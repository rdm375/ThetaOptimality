#pragma once
#include "linalg.hpp"
#include <optional>
#include <vector>

namespace theta_opt {
// Each input is a persistent qualification time under the same criterion and
// horizon. All goals qualify from the latest individual time, if every goal
// qualifies. A missing time means the joint requirement is not established.
template<std::size_t Goals>
std::optional<double> all_goals_qualification_time(const std::array<std::optional<double>,Goals>& times){
    static_assert(Goals>0);
    std::optional<double> joint;
    for(auto time:times){
        if(!time)return std::nullopt;
        if(!std::isfinite(*time))throw std::invalid_argument("Nonfinite goal qualification time");
        joint=joint?std::max(*joint,*time):*time;
    }
    return joint;
}
// A simultaneous noise box for m independent N(0,sd^2) errors. Its nonlinear
// model-consistency preimage is a confidence set without a linearization.
inline double independent_gaussian_box_radius(double sd,std::size_t measurements,double joint_probability){
    if(!std::isfinite(sd)||sd<=0||!measurements||!std::isfinite(joint_probability)||joint_probability<=0||joint_probability>=1)
        throw std::invalid_argument("Gaussian noise-box arguments");
    double target=.5*(1+std::exp(std::log(joint_probability)/static_cast<double>(measurements)));
    if(target>=1)throw std::invalid_argument("Gaussian box probability rounds to one");
    auto cdf=[](double z){return .5*std::erfc(-z/std::sqrt(2.));};
    double lo=0.,hi=1.;while(cdf(hi)<target)hi*=2.;
    for(int k=0;k<100;++k){double mid=.5*(lo+hi);if(cdf(mid)<target)lo=mid;else hi=mid;}
    return sd*.5*(lo+hi);
}
// Local linear uncertainty indicator, not a nonlinear coverage guarantee.
template<std::size_t P>
double local_goal_variance(const Vec<P>& gradient,const Mat<P,P>& covariance){
    double variance=dot(gradient,mv(covariance,gradient));
    if(!std::isfinite(variance)||variance< -1e-14)throw std::domain_error("Invalid goal variance");
    return std::max(0.,variance);
}
template<std::size_t P>
double box_support(const Vec<P>& gradient,const Vec<P>& radii){
    double support=0.;
    for(std::size_t i=0;i<P;++i){if(!std::isfinite(radii[i])||radii[i]<0)throw std::domain_error("Invalid box radius");support+=std::abs(gradient[i])*radii[i];}
    return support;
}
// Earliest node for which every subsequent sampled ratio is <= 1.
// This is a finite-horizon NODE criterion; interpolation is not a certificate.
inline std::optional<double> persistent_node_qualification(
    const std::vector<double>& times,const std::vector<double>& ratios){
    if(times.empty()||times.size()!=ratios.size())throw std::invalid_argument("Qualification sizes");
    std::optional<double> answer;
    for(std::size_t n=0;n<times.size();++n){
        if(!std::isfinite(times[n])||(n&&times[n]<=times[n-1])||!std::isfinite(ratios[n])||ratios[n]<0)throw std::invalid_argument("Invalid qualification data");
    }
    for(std::size_t n=times.size();n-->0;){if(ratios[n]>1)break;answer=times[n];}
    return answer;
}
// Continuous error comparison E' <= mu*E + D with constant bounds per step.
// mu may be positive: asymptotic stability does not imply norm contraction.
inline double propagate_defect_bound(double error,double defect,double step,double logarithmic_norm_bound){
    if(!std::isfinite(error)||!std::isfinite(defect)||!std::isfinite(step)||!std::isfinite(logarithmic_norm_bound)||error<0||defect<0||step<=0)throw std::invalid_argument("Defect comparison arguments");
    double z=logarithmic_norm_bound*step;
    double integral=logarithmic_norm_bound==0?step:std::expm1(z)/logarithmic_norm_bound;
    double result=std::exp(z)*error+integral*defect;if(!std::isfinite(result))throw std::overflow_error("Defect comparison overflow");return result;
}
inline double quadratic_absolute_maximum(double c0,double c1,double c2){
    if(!std::isfinite(c0)||!std::isfinite(c1)||!std::isfinite(c2))throw std::invalid_argument("Finite polynomial coefficients required");
    double result=std::max(std::abs(c0),std::abs(c0+c1+c2));if(c2!=0){double s=-c1/(2*c2);if(s>0&&s<1)result=std::max(result,std::abs(c0+c1*s+c2*s*s));}return result;
}
template<std::size_t N>
Vec<N> linear_quadratic_reconstruction(const Mat<N,N>& A,const Vec<N>& x,const Vec<N>& next,double step,double fraction){
    if(!std::isfinite(step)||step<=0||!std::isfinite(fraction)||fraction<0||fraction>1)throw std::invalid_argument("Linear reconstruction arguments");
    return x+step*fraction*mv(A,x)+fraction*fraction*(next-x-step*mv(A,x));
}
// Residual bound for a quadratic reconstruction of a constant linear theta step.
// The caller supplies the chosen physical norm; no contraction is assumed here.
template<std::size_t N,class Norm>
double linear_theta_quadratic_defect_bound(const Mat<N,N>& A,const Vec<N>& x,const Vec<N>& next,double step,double theta,Norm physical_norm){
    if(!std::isfinite(step)||step<=0||!std::isfinite(theta)||theta<0||theta>1)throw std::invalid_argument("Linear theta defect arguments");
    auto norm_value=[&](const Vec<N>& v){(void)norm(v);double result=physical_norm(v);if(!std::isfinite(result)||result<0)throw std::invalid_argument("Invalid physical norm");return result;};
    auto delta=next-x,Ax=mv(A,x),Anext=mv(A,next),Adelta=mv(A,delta),r=delta-step*((1-theta)*Ax+theta*Anext);
    double roundoff=64*std::numeric_limits<double>::epsilon()*(norm_value(x)+norm_value(next)+step*(norm_value(Ax)+norm_value(Anext)));
    return std::abs(2*theta-1)*norm_value(Adelta)+step*theta*.25*norm_value(mv(A,Adelta))+.25*norm_value(mv(A,r))+2*(norm_value(r)+roundoff)/step;
}
template<std::size_t Goals> struct IntervalQualification {
    std::array<std::optional<double>,Goals> first_below{},persistent{};std::optional<double> joint;
};
// Bounds must cover each WHOLE interval [edges[n],edges[n+1]], not sampled nodes.
template<std::size_t Goals>
IntervalQualification<Goals> qualify_intervals(const std::vector<double>& edges,const std::vector<std::array<double,Goals>>& bounds,const std::array<double,Goals>& tolerances){
    if(edges.size()<2||bounds.size()+1!=edges.size())throw std::invalid_argument("Interval qualification sizes");
    for(std::size_t n=0;n<edges.size();++n)if(!std::isfinite(edges[n])||(n&&edges[n]<=edges[n-1]))throw std::invalid_argument("Increasing finite interval edges required");
    IntervalQualification<Goals> result;std::vector<double> starts(edges.begin(),edges.end()-1);
    for(std::size_t j=0;j<Goals;++j){if(!std::isfinite(tolerances[j])||tolerances[j]<=0)throw std::invalid_argument("Positive finite goal tolerances required");std::vector<double> ratios;for(std::size_t n=0;n<bounds.size();++n){if(!std::isfinite(bounds[n][j])||bounds[n][j]<0)throw std::invalid_argument("Finite nonnegative interval bounds required");ratios.push_back(bounds[n][j]/tolerances[j]);if(ratios.back()<=1&&!result.first_below[j])result.first_below[j]=edges[n];}result.persistent[j]=persistent_node_qualification(starts,ratios);}result.joint=all_goals_qualification_time(result.persistent);return result;
}

}
