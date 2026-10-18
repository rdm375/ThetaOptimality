#pragma once
#include "system.hpp"
#include <type_traits>
#include <string>

namespace theta_opt {
struct OptimizeOptions {double gradient_tolerance=1e-11;std::size_t max_iterations=40,max_backtracks=40;double cost_roundoff_slack=0.;};
template<std::size_t P> struct Iteration {Vec<P> control{};double cost{},gradient_norm{},step_length{};};
template<std::size_t P> struct Optimum {Vec<P> control{};std::vector<Iteration<P>> trace;};
// A small dense Newton driver for examples. Exact H columns come from H*v.
template<std::size_t N,std::size_t P>
Optimum<P> minimize(const System<N,P>& system,Vec<P> p,OptimizeOptions options={},
                  std::type_identity_t<std::function<bool(const Vec<P>&)>> admissible=[](const Vec<P>&){return true;}){
    if(!std::isfinite(options.gradient_tolerance)||options.gradient_tolerance<=0||!options.max_iterations||!options.max_backtracks||!std::isfinite(options.cost_roundoff_slack)||options.cost_roundoff_slack<0)throw std::invalid_argument("Optimizer options");
    if(!admissible(p))throw std::domain_error("Inadmissible initial control");
    Optimum<P> out;
    for(std::size_t k=0;k<=options.max_iterations;++k){auto e=system.evaluate(p);out.trace.push_back({p,e.cost,norm(e.gradient),0.});
        if(norm(e.gradient)<=options.gradient_tolerance){out.control=p;return out;}
        if(k==options.max_iterations)break;
        Mat<P,P> H{};
        for(std::size_t j=0;j<P;++j){Vec<P> unit{};unit[j]=1.;auto col=system.hessian_vector(p,unit);for(std::size_t i=0;i<P;++i)H[i][j]=col[i];}
        Vec<P> step=-1.*e.gradient;
        try {auto trial=LU<P>(H).solve(-1.*e.gradient);if(dot(e.gradient,trial)<0)step=trial;}catch(const std::runtime_error&){/* Gradient fallback. */}
        bool accepted=false;double alpha=1.;
        for(std::size_t j=0;j<options.max_backtracks;++j){auto trial=p+alpha*step;
            if(admissible(trial)){
                // Near stationarity the predicted decrease can be below one ulp.
                double slack=options.cost_roundoff_slack+8*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(e.cost));
                try {if(system.evaluate(trial).cost<=e.cost+1e-4*alpha*dot(e.gradient,step)+slack){p=trial;accepted=true;break;}}
                catch(const std::domain_error&){/* Try a shorter admissible step. */}
            }alpha*=.5;
        }
        if(!accepted)throw std::runtime_error("Optimizer line search failed at iteration "+std::to_string(k)+" gradient "+std::to_string(norm(e.gradient))+" cost "+std::to_string(e.cost));
        out.trace.back().step_length=alpha;
    }throw std::runtime_error("Optimizer failed to converge");
}
}
