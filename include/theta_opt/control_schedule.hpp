#pragma once
#include "system.hpp"

namespace theta_opt {

// A finite vector of optimization coefficients can parameterize a control
// trajectory. The callbacks expose the schedule Jacobian and its directional
// derivative so compose_dynamics preserves exact tangent/adjoint/Hessian
// derivatives through the time-varying control map.
template<std::size_t C,std::size_t P> struct ControlParameterization {
    std::function<Vec<C>(double,const Vec<P>&)> value;
    std::function<Mat<C,P>(double,const Vec<P>&)> jacobian;
    std::function<Mat<C,P>(double,const Vec<P>&,const Vec<P>&)> directional_jacobian;
};

template<std::size_t Channels,std::size_t Knots>
class PiecewiseLinearControl {
    static_assert(Channels>0&&Knots>=2);
    std::array<double,Knots> time_{};
    std::pair<std::size_t,double> interval(double t) const {
        if(!std::isfinite(t))throw std::invalid_argument("Nonfinite control time");
        if(t<=time_.front())return {0,0.};
        if(t>=time_.back())return {Knots-2,1.};
        auto upper=std::upper_bound(time_.begin(),time_.end(),t);
        const auto i=static_cast<std::size_t>(upper-time_.begin()-1);
        return {i,(t-time_[i])/(time_[i+1]-time_[i])};
    }
public:
    static constexpr std::size_t parameter_count=Channels*Knots;
    explicit PiecewiseLinearControl(std::array<double,Knots> time):time_(time) {
        for(std::size_t i=0;i<Knots;++i)
            if(!std::isfinite(time_[i])||(i&&!(time_[i]>time_[i-1])))
                throw std::invalid_argument("Control knot times must be finite and increasing");
    }
    Vec<Channels> value(double t,const Vec<parameter_count>& coefficient) const {
        const auto [i,s]=interval(t);Vec<Channels> result{};
        for(std::size_t c=0;c<Channels;++c)
            result[c]=(1.-s)*coefficient[c*Knots+i]+s*coefficient[c*Knots+i+1];
        return result;
    }
    Mat<Channels,parameter_count> jacobian(double t,const Vec<parameter_count>&) const {
        const auto [i,s]=interval(t);Mat<Channels,parameter_count> result{};
        for(std::size_t c=0;c<Channels;++c){result[c][c*Knots+i]=1.-s;result[c][c*Knots+i+1]=s;}
        return result;
    }
    Mat<Channels,parameter_count> directional_jacobian(double,const Vec<parameter_count>&,
                                                         const Vec<parameter_count>&) const {
        return {};
    }
    ControlParameterization<Channels,parameter_count> parameterization() const {
        auto copy=*this;ControlParameterization<Channels,parameter_count> map;
        map.value=[copy](double t,const Vec<parameter_count>& p){return copy.value(t,p);};
        map.jacobian=[copy](double t,const Vec<parameter_count>& p){return copy.jacobian(t,p);};
        map.directional_jacobian=[copy](double t,const Vec<parameter_count>& p,const Vec<parameter_count>& d){return copy.directional_jacobian(t,p,d);};
        return map;
    }
};

// A positive, sum-to-one schedule for compositions. Species 0 is the softmax
// reference; the other species use log-ratio coefficients at each knot. Thus
// every interpolated composition lies strictly inside the simplex.
template<std::size_t Species,std::size_t Knots>
class SoftmaxCompositionControl {
    static_assert(Species>=2&&Knots>=2);
    static constexpr std::size_t P=(Species-1)*Knots;
    std::array<double,Knots> time_{};
    std::pair<std::size_t,double> interval(double t) const {
        if(!std::isfinite(t))throw std::invalid_argument("Nonfinite composition time");
        if(t<=time_.front())return {0,0.};
        if(t>=time_.back())return {Knots-2,1.};
        auto upper=std::upper_bound(time_.begin(),time_.end(),t);
        const auto i=static_cast<std::size_t>(upper-time_.begin()-1);
        return {i,(t-time_[i])/(time_[i+1]-time_[i])};
    }
    Vec<Species> composition(double t,const Vec<P>& z) const {
        const auto [k,s]=interval(t);Vec<Species> logits{};double maximum=0.;
        for(std::size_t species=1;species<Species;++species){
            const auto offset=(species-1)*Knots;
            logits[species]=(1.-s)*z[offset+k]+s*z[offset+k+1];
            maximum=std::max(maximum,logits[species]);
        }
        double sum=0.;for(auto& q:logits){q=std::exp(q-maximum);sum+=q;}
        for(auto& q:logits)q/=sum;
        return logits;
    }
public:
    static constexpr std::size_t parameter_count=P;
    explicit SoftmaxCompositionControl(std::array<double,Knots> time):time_(time) {
        for(std::size_t i=0;i<Knots;++i)
            if(!std::isfinite(time_[i])||(i&&!(time_[i]>time_[i-1])))
                throw std::invalid_argument("Composition knot times must be finite and increasing");
    }
    Vec<Species> value(double t,const Vec<P>& z) const {return composition(t,z);}
    Mat<Species,P> jacobian(double t,const Vec<P>& z) const {
        const auto [k,s]=interval(t);const auto y=composition(t,z);Mat<Species,P> result{};
        for(std::size_t species=1;species<Species;++species){
            const auto offset=(species-1)*Knots;
            for(std::size_t knot_side=0;knot_side<2;++knot_side){
                const std::size_t knot=k+knot_side;const double weight=knot_side?s:1.-s;
                for(std::size_t i=0;i<Species;++i)
                    result[i][offset+knot]=weight*y[i]*((i==species?1.:0.)-y[species]);
            }
        }
        return result;
    }
    Mat<Species,P> directional_jacobian(double t,const Vec<P>& z,const Vec<P>& dz) const {
        const auto [k,s]=interval(t);const auto y=composition(t,z);Vec<Species> dlogit{};
        for(std::size_t species=1;species<Species;++species){
            const auto offset=(species-1)*Knots;
            dlogit[species]=(1.-s)*dz[offset+k]+s*dz[offset+k+1];
        }
        double mean=0.;for(std::size_t i=0;i<Species;++i)mean+=y[i]*dlogit[i];
        Vec<Species> dy{};for(std::size_t i=0;i<Species;++i)dy[i]=y[i]*(dlogit[i]-mean);
        Mat<Species,P> result{};
        for(std::size_t species=1;species<Species;++species){
            const auto offset=(species-1)*Knots;
            for(std::size_t side=0;side<2;++side){
                const std::size_t knot=k+side;const double weight=side?s:1.-s;
                for(std::size_t i=0;i<Species;++i){
                    const double q=(i==species?1.:0.)-y[species];
                    result[i][offset+knot]=weight*(dy[i]*q-y[i]*dy[species]);
                }
            }
        }
        return result;
    }
    ControlParameterization<Species,P> parameterization() const {
        auto copy=*this;ControlParameterization<Species,P> map;
        map.value=[copy](double t,const Vec<P>& z){return copy.value(t,z);};
        map.jacobian=[copy](double t,const Vec<P>& z){return copy.jacobian(t,z);};
        map.directional_jacobian=[copy](double t,const Vec<P>& z,const Vec<P>& d){return copy.directional_jacobian(t,z,d);};
        return map;
    }
};

// Chain a model whose independent control is u(t) through a finite-dimensional
// schedule z -> u(t;z). State and control Jacobians, including directional
// second derivatives, are composed analytically.
template<std::size_t N,std::size_t C,std::size_t P>
Dynamics<N,P> compose_dynamics(Dynamics<N,C> model,ControlParameterization<C,P> schedule) {
    if(!model.value||!model.derivatives||!model.directional_derivatives||
       !schedule.value||!schedule.jacobian||!schedule.directional_jacobian)
        throw std::invalid_argument("Incomplete scheduled dynamics callbacks");
    Dynamics<N,P> result;
    result.value=[model,schedule](double t,const Vec<N>& x,const Vec<P>& z){return model.value(t,x,schedule.value(t,z));};
    result.derivatives=[model,schedule](double t,const Vec<N>& x,const Vec<P>& z){
        const auto u=schedule.value(t,z);const auto d=model.derivatives(t,x,u);const auto J=schedule.jacobian(t,z);
        Derivatives<N,P> out;out.x=d.x;
        for(std::size_t i=0;i<N;++i)for(std::size_t j=0;j<P;++j)
            for(std::size_t c=0;c<C;++c)out.control[i][j]+=d.control[i][c]*J[c][j];
        return out;
    };
    result.directional_derivatives=[model,schedule](double t,const Vec<N>& x,const Vec<P>& z,
                                                     const Vec<N>& dx,const Vec<P>& dz){
        const auto u=schedule.value(t,z);const auto J=schedule.jacobian(t,z);const auto du=mv(J,dz);
        const auto directional=model.directional_derivatives(t,x,u,dx,du);
        const auto dJ=schedule.directional_jacobian(t,z,dz);const auto first=model.derivatives(t,x,u);
        Derivatives<N,P> out;out.x=directional.x;
        for(std::size_t i=0;i<N;++i)for(std::size_t j=0;j<P;++j)
            for(std::size_t c=0;c<C;++c)
                out.control[i][j]+=directional.control[i][c]*J[c][j]+first.control[i][c]*dJ[c][j];
        return out;
    };
    result.admissible=[model,schedule](double t,const Vec<N>& x,const Vec<P>& z){
        return model.admissible(t,x,schedule.value(t,z));
    };
    return result;
}

// Chain an observation, running/terminal term, or path inequality that depends
// on the instantaneous control through the same schedule map.
template<std::size_t N,std::size_t C,std::size_t P>
ScalarFunction<N,P> compose_scalar_function(ScalarFunction<N,C> function,
                                             ControlParameterization<C,P> schedule) {
    if(!function.evaluate||!function.directional_gradient||!schedule.value||
       !schedule.jacobian||!schedule.directional_jacobian)
        throw std::invalid_argument("Incomplete scheduled scalar callbacks");
    ScalarFunction<N,P> result;
    result.evaluate=[function,schedule](double t,const Vec<N>& x,const Vec<P>& z){
        const auto value=function.evaluate(t,x,schedule.value(t,z));
        const auto J=schedule.jacobian(t,z);Vec<P> control_gradient{};
        for(std::size_t j=0;j<P;++j)for(std::size_t c=0;c<C;++c)
            control_gradient[j]+=value.control[c]*J[c][j];
        return ScalarValue<N,P>{value.value,value.x,control_gradient};
    };
    result.directional_gradient=[function,schedule](double t,const Vec<N>& x,const Vec<P>& z,
                                                     const Vec<N>& dx,const Vec<P>& dz){
        const auto u=schedule.value(t,z);const auto J=schedule.jacobian(t,z);
        const auto du=mv(J,dz);const auto first=function.evaluate(t,x,u);
        const auto directional=function.directional_gradient(t,x,u,dx,du);
        const auto dJ=schedule.directional_jacobian(t,z,dz);ScalarDirection<N,P> out;
        out.x=directional.x;
        for(std::size_t j=0;j<P;++j)for(std::size_t c=0;c<C;++c)
            out.control[j]+=directional.control[c]*J[c][j]+first.control[c]*dJ[c][j];
        return out;
    };
    return result;
}
}
