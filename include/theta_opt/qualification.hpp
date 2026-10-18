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
}
