#pragma once
#include "linalg.hpp"
#include <functional>
#include <utility>
#include <vector>

namespace theta_opt {
struct Grid {
    std::vector<double> time,theta;
    Grid(std::vector<double> t,std::vector<double> th):time(std::move(t)),theta(std::move(th)){
        if(time.size()<2||theta.size()+1!=time.size())throw std::invalid_argument("Grid sizes");
        for(std::size_t n=0;n<time.size();++n){if(!std::isfinite(time[n])||(n&&(!(time[n]>time[n-1])||!std::isfinite(time[n]-time[n-1]))))throw std::invalid_argument("Times and positive intervals must be finite");}
        for(auto x:theta)if(!std::isfinite(x)||x<0||x>1)throw std::invalid_argument("Theta must lie in [0,1]");
    }
    static Grid uniform(double t0,double tf,std::size_t steps,double th){
        if(!steps||!std::isfinite(t0)||!std::isfinite(tf)||!(tf>t0))throw std::invalid_argument("Uniform grid");
        std::vector<double> t(steps+1);for(std::size_t n=0;n<=steps;++n)t[n]=t0+(tf-t0)*static_cast<double>(n)/static_cast<double>(steps);
        return Grid(std::move(t),std::vector<double>(steps,th));
    }
    std::size_t steps() const{return theta.size();}
    std::vector<double> weights() const{
        std::vector<double> w(time.size(),0);
        for(std::size_t n=0;n<steps();++n){double h=time[n+1]-time[n];w[n]+=h*(1-theta[n]);w[n+1]+=h*theta[n];}return w;
    }
};
struct NewtonOptions {double absolute_tolerance=1e-13,relative_tolerance=1e-13;std::size_t max_iterations=30,max_backtracks=30;};
template<std::size_t N,std::size_t P> struct Derivatives {Mat<N,N> x{};Mat<N,P> control{};};
template<std::size_t N,std::size_t P> struct ScalarValue {double value{};Vec<N> x{};Vec<P> control{};};
template<std::size_t N,std::size_t P> struct ScalarDirection {Vec<N> x{};Vec<P> control{};};

// All directions differentiate at fixed time and fixed prescribed grid.
template<std::size_t N,std::size_t P> struct Dynamics {
    std::function<Vec<N>(double,const Vec<N>&,const Vec<P>&)> value;
    std::function<Derivatives<N,P>(double,const Vec<N>&,const Vec<P>&)> derivatives;
    std::function<Derivatives<N,P>(double,const Vec<N>&,const Vec<P>&,const Vec<N>&,const Vec<P>&)> directional_derivatives;
    std::function<bool(double,const Vec<N>&,const Vec<P>&)> admissible=[](double,const Vec<N>&,const Vec<P>&){return true;};
};
template<std::size_t N,std::size_t P> struct ScalarFunction {
    std::function<ScalarValue<N,P>(double,const Vec<N>&,const Vec<P>&)> evaluate=[](double,const Vec<N>&,const Vec<P>&){return ScalarValue<N,P>{};};
    // Optional for first order, required explicitly for Hessian products.
    std::function<ScalarDirection<N,P>(double,const Vec<N>&,const Vec<P>&,const Vec<N>&,const Vec<P>&)> directional_gradient;
    static ScalarFunction zero(){ScalarFunction f;f.directional_gradient=[](double,const Vec<N>&,const Vec<P>&,const Vec<N>&,const Vec<P>&){return ScalarDirection<N,P>{};};return f;}
};
template<std::size_t N,std::size_t P> struct InitialCondition {
    std::function<Vec<N>(const Vec<P>&)> value;
    std::function<Mat<N,P>(const Vec<P>&)> jacobian;
    std::function<Mat<N,P>(const Vec<P>&,const Vec<P>&)> directional_jacobian;
};
template<std::size_t P> struct Penalty {
    std::function<double(const Vec<P>&)> value=[](const Vec<P>&){return 0.;};
    std::function<Vec<P>(const Vec<P>&)> gradient=[](const Vec<P>&){return Vec<P>{};};
    std::function<Vec<P>(const Vec<P>&,const Vec<P>&)> hessian_vector;
    static Penalty zero(){Penalty p;p.hessian_vector=[](const Vec<P>&,const Vec<P>&){return Vec<P>{};};return p;}
};
template<std::size_t N,std::size_t P> struct Problem {
    Dynamics<N,P> dynamics;InitialCondition<N,P> initial;
    ScalarFunction<N,P> running=ScalarFunction<N,P>::zero(),terminal=ScalarFunction<N,P>::zero();
    Penalty<P> penalty=Penalty<P>::zero();
    // Isolated measurements are nodal costs, not running-cost densities.
    std::vector<std::pair<std::size_t,ScalarFunction<N,P>>> nodal;
};
template<std::size_t N> struct Trajectory {std::vector<Vec<N>> state;double max_residual{};std::size_t newton_updates{};};
template<std::size_t N,std::size_t P> struct Evaluation {
    Trajectory<N> trajectory;double cost{};Vec<P> gradient{};
    std::vector<Vec<N>> nodal_adjoint,step_multiplier;
};

template<std::size_t N,std::size_t P> class System {
    static_assert(N>0&&P>0);
    Problem<N,P> problem_;Grid grid_;NewtonOptions options_;
    auto deriv(double t,const Vec<N>& x,const Vec<P>& p) const{return problem_.dynamics.derivatives(t,x,p);}
    std::pair<Mat<N,N>,Mat<N,N>> blocks(std::size_t n,const std::vector<Vec<N>>& x,const Vec<P>& p) const{
        double h=grid_.time[n+1]-grid_.time[n],th=grid_.theta[n];
        return {add(identity<N>(),scaled(-h*th,deriv(grid_.time[n+1],x[n+1],p).x)),
                add(identity<N>(),scaled(h*(1-th),deriv(grid_.time[n],x[n],p).x))};
    }
    Mat<N,P> forcing(std::size_t n,const std::vector<Vec<N>>& x,const Vec<P>& p) const{
        double h=grid_.time[n+1]-grid_.time[n],th=grid_.theta[n];
        return add(scaled(h*(1-th),deriv(grid_.time[n],x[n],p).control),scaled(h*th,deriv(grid_.time[n+1],x[n+1],p).control));
    }
public:
    System(Problem<N,P> problem,Grid grid,NewtonOptions options={}):problem_(std::move(problem)),grid_(std::move(grid)),options_(options){
        if(!problem_.dynamics.value||!problem_.dynamics.derivatives||!problem_.dynamics.admissible||!problem_.initial.value||!problem_.initial.jacobian||
           !problem_.running.evaluate||!problem_.terminal.evaluate||!problem_.penalty.value||!problem_.penalty.gradient)throw std::invalid_argument("Missing first-order callback");
        if(!std::isfinite(options_.absolute_tolerance)||!std::isfinite(options_.relative_tolerance)||options_.absolute_tolerance<=0||options_.relative_tolerance<0||!options_.max_iterations||!options_.max_backtracks)throw std::invalid_argument("Newton options");
        for(const auto& [node,cost]:problem_.nodal)
            if(node>=grid_.time.size()||!cost.evaluate)throw std::invalid_argument("Invalid nodal cost");
    }
    const Grid& grid() const{return grid_;}
    Trajectory<N> forward(const Vec<P>& p) const{
        (void)norm(p);Trajectory<N> out;out.state.resize(grid_.time.size());out.state[0]=problem_.initial.value(p);(void)norm(out.state[0]);
        if(!problem_.dynamics.admissible(grid_.time[0],out.state[0],p))throw std::domain_error("Inadmissible initial state");
        for(std::size_t n=0;n<grid_.steps();++n){
            double h=grid_.time[n+1]-grid_.time[n],th=grid_.theta[n],tn=grid_.time[n+1];
            const auto old=out.state[n],fold=problem_.dynamics.value(grid_.time[n],old,p);auto next=old;
            auto residual=[&](const Vec<N>& y){return y-old-h*((1-th)*fold+th*problem_.dynamics.value(tn,y,p));};
            bool converged=false;
            for(std::size_t k=0;k<=options_.max_iterations;++k){
                auto r=residual(next);double error=norm(r),tol=options_.absolute_tolerance+options_.relative_tolerance*std::max(norm(old),norm(next));
                if(error<=tol){
                    if(!problem_.dynamics.admissible(tn,next,p))throw std::domain_error("Inadmissible converged state");
                    converged=true;out.max_residual=std::max(out.max_residual,error);break;
                }
                if(k==options_.max_iterations)break;
                auto bn=add(identity<N>(),scaled(-h*th,deriv(tn,next,p).x));auto change=LU<N>(bn).solve(-1.*r);
                double alpha=1.;bool accepted=false;
                for(std::size_t j=0;j<options_.max_backtracks;++j){auto trial=next+alpha*change;
                    if(problem_.dynamics.admissible(tn,trial,p)&&norm(residual(trial))<=(1-1e-4*alpha)*error){next=trial;accepted=true;break;}alpha*=.5;
                }
                if(!accepted)throw std::runtime_error("Forward Newton line search failed");
                ++out.newton_updates;
            }
            if(!converged)throw std::runtime_error("Forward Newton failed to converge");
            out.state[n+1]=next;
        }return out;
    }
    std::vector<Vec<N>> tangent(const Vec<P>& p,const Trajectory<N>& trajectory,const Vec<P>& direction) const{
        const auto& x=trajectory.state;if(x.size()!=grid_.time.size())throw std::invalid_argument("Trajectory size");
        std::vector<Vec<N>> s(x.size());s[0]=mv(problem_.initial.jacobian(p),direction);
        for(std::size_t n=0;n<grid_.steps();++n){auto [bn,cn]=blocks(n,x,p);s[n+1]=LU<N>(bn).solve(mv(cn,s[n])+mv(forcing(n,x,p),direction));}return s;
    }
    double cost(const Vec<P>& p,const Trajectory<N>& trajectory) const{
        const auto& x=trajectory.state;
        if(x.size()!=grid_.time.size())throw std::invalid_argument("Trajectory size");
        auto w=grid_.weights();double value=problem_.penalty.value(p);
        for(std::size_t n=0;n<x.size();++n)value+=w[n]*problem_.running.evaluate(grid_.time[n],x[n],p).value;
        value+=problem_.terminal.evaluate(grid_.time.back(),x.back(),p).value;
        for(const auto& [node,cost]:problem_.nodal)value+=cost.evaluate(grid_.time[node],x[node],p).value;
        if(!std::isfinite(value))throw std::runtime_error("Nonfinite cost");
        return value;
    }
    Evaluation<N,P> evaluate(const Vec<P>& p) const{return adjoint(p,forward(p));}
    Evaluation<N,P> adjoint(const Vec<P>& p,Trajectory<N> trajectory) const{
        if(trajectory.state.size()!=grid_.time.size())throw std::invalid_argument("Trajectory size");
        Evaluation<N,P> e;e.trajectory=std::move(trajectory);const auto& x=e.trajectory.state;auto w=grid_.weights();
        e.nodal_adjoint.resize(x.size());e.step_multiplier.resize(grid_.steps());e.cost=problem_.penalty.value(p);e.gradient=problem_.penalty.gradient(p);
        for(std::size_t n=0;n<x.size();++n){auto c=problem_.running.evaluate(grid_.time[n],x[n],p);e.cost+=w[n]*c.value;e.nodal_adjoint[n]=w[n]*c.x;e.gradient=e.gradient+w[n]*c.control;}
        auto phi=problem_.terminal.evaluate(grid_.time.back(),x.back(),p);e.cost+=phi.value;e.nodal_adjoint.back()=e.nodal_adjoint.back()+phi.x;e.gradient=e.gradient+phi.control;
        for(const auto& [node,cost]:problem_.nodal){auto v=cost.evaluate(grid_.time[node],x[node],p);e.cost+=v.value;e.nodal_adjoint[node]=e.nodal_adjoint[node]+v.x;e.gradient=e.gradient+v.control;}
        for(std::size_t n=grid_.steps();n-->0;){auto [bn,cn]=blocks(n,x,p);auto z=LU<N>(transpose(bn)).solve(e.nodal_adjoint[n+1]);e.step_multiplier[n]=z;
            e.nodal_adjoint[n]=e.nodal_adjoint[n]+mv(transpose(cn),z);e.gradient=e.gradient+mv(transpose(forcing(n,x,p)),z);
        }
        e.gradient=e.gradient+mv(transpose(problem_.initial.jacobian(p)),e.nodal_adjoint[0]);
        if(!std::isfinite(e.cost))throw std::runtime_error("Nonfinite cost");
        (void)norm(e.gradient);return e;
    }
    // Pull back arbitrary nodal state sources through the SAME discrete model.
    // Useful for goal gradients and J^T R^-1 J information actions.
    Vec<P> state_pullback(const Vec<P>& p,const Trajectory<N>& trajectory,
                          std::vector<Vec<N>> sources) const{
        const auto& x=trajectory.state;
        if(x.size()!=grid_.time.size()||sources.size()!=x.size())throw std::invalid_argument("Pullback size");
        Vec<P> gradient{};
        for(std::size_t n=grid_.steps();n-->0;){auto [bn,cn]=blocks(n,x,p);
            auto z=LU<N>(transpose(bn)).solve(sources[n+1]);
            sources[n]=sources[n]+mv(transpose(cn),z);
            gradient=gradient+mv(transpose(forcing(n,x,p)),z);
        }
        return gradient+mv(transpose(problem_.initial.jacobian(p)),sources[0]);
    }
    Vec<P> goal_gradient(const Vec<P>& p,const Trajectory<N>& trajectory,
                         std::size_t node,const Vec<N>& state_gradient) const{
        if(node>=grid_.time.size())throw std::invalid_argument("Goal node");
        std::vector<Vec<N>> sources(grid_.time.size());sources[node]=state_gradient;
        return state_pullback(p,trajectory,std::move(sources));
    }
    // Converts the paper's nodal adjoint to the remaining-cost convention.
    // Only the fixed-grid running-cost endpoint term is removed.
    std::vector<Vec<N>> cost_to_go_adjoint(const Vec<P>& p,const Evaluation<N,P>& e) const{
        if(e.nodal_adjoint.size()!=grid_.time.size()||e.trajectory.state.size()!=grid_.time.size())throw std::invalid_argument("Evaluation size");
        auto lam=e.nodal_adjoint;
        for(std::size_t n=1;n<lam.size();++n){double previous=(grid_.time[n]-grid_.time[n-1])*grid_.theta[n-1];lam[n]=lam[n]-previous*problem_.running.evaluate(grid_.time[n],e.trajectory.state[n],p).x;}return lam;
    }
    Vec<P> hessian_vector(const Vec<P>& p,const Vec<P>& direction) const{
        if(!problem_.dynamics.directional_derivatives||!problem_.initial.directional_jacobian||!problem_.running.directional_gradient||!problem_.terminal.directional_gradient||!problem_.penalty.hessian_vector)throw std::invalid_argument("Hessian-vector product requires all second-order callbacks");
        auto e=evaluate(p);const auto& x=e.trajectory.state;auto s=tangent(p,e.trajectory,direction);auto w=grid_.weights();
        std::vector<Vec<N>> dl(x.size());auto hv=problem_.penalty.hessian_vector(p,direction);
        for(std::size_t n=0;n<x.size();++n){auto dc=problem_.running.directional_gradient(grid_.time[n],x[n],p,s[n],direction);dl[n]=w[n]*dc.x;hv=hv+w[n]*dc.control;}
        auto dp=problem_.terminal.directional_gradient(grid_.time.back(),x.back(),p,s.back(),direction);dl.back()=dl.back()+dp.x;hv=hv+dp.control;
        for(const auto& [node,cost]:problem_.nodal){
            if(!cost.directional_gradient)throw std::invalid_argument("Nodal cost requires second-order callback");
            auto dc=cost.directional_gradient(grid_.time[node],x[node],p,s[node],direction);
            dl[node]=dl[node]+dc.x;hv=hv+dc.control;
        }
        for(std::size_t n=grid_.steps();n-->0;){double h=grid_.time[n+1]-grid_.time[n],th=grid_.theta[n];auto [bn,cn]=blocks(n,x,p);
            auto old=problem_.dynamics.directional_derivatives(grid_.time[n],x[n],p,s[n],direction),next=problem_.dynamics.directional_derivatives(grid_.time[n+1],x[n+1],p,s[n+1],direction);
            auto db=scaled(-h*th,next.x),dc=scaled(h*(1-th),old.x);
            auto dz=LU<N>(transpose(bn)).solve(dl[n+1]-mv(transpose(db),e.step_multiplier[n]));
            dl[n]=dl[n]+mv(transpose(dc),e.step_multiplier[n])+mv(transpose(cn),dz);
            auto dk=add(scaled(h*(1-th),old.control),scaled(h*th,next.control));
            hv=hv+mv(transpose(dk),e.step_multiplier[n])+mv(transpose(forcing(n,x,p)),dz);
        }
        hv=hv+mv(transpose(problem_.initial.directional_jacobian(p,direction)),e.nodal_adjoint[0])+mv(transpose(problem_.initial.jacobian(p)),dl[0]);(void)norm(hv);return hv;
    }
};
}
