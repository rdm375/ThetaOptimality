// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include "system.hpp"
#include <algorithm>
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

template<std::size_t P> struct BoxOptimum {
    Vec<P> control{},gradient{},lower_multiplier{},upper_multiplier{};
    std::array<bool,P> active_lower{},active_upper{};
    std::vector<Iteration<P>> trace;
    double kkt_residual{};
};

// Dense active-set Newton method for simple bound constraints l <= p <= u.
// Iteration.gradient_norm records the projected-gradient norm; the returned
// multipliers use stationarity g - lambda_lower + lambda_upper = 0.
template<std::size_t N,std::size_t P>
BoxOptimum<P> minimize_box_constrained(const System<N,P>& system,Vec<P> p,
        const Vec<P>& lower,const Vec<P>& upper,OptimizeOptions options={}) {
    if(!std::isfinite(options.gradient_tolerance)||options.gradient_tolerance<=0||
       !options.max_iterations||!options.max_backtracks||
       !std::isfinite(options.cost_roundoff_slack)||options.cost_roundoff_slack<0)
        throw std::invalid_argument("Box optimizer options");
    for(std::size_t j=0;j<P;++j) {
        if(!std::isfinite(lower[j])||!std::isfinite(upper[j])||!(lower[j]<upper[j])||
           !std::isfinite(p[j])||p[j]<lower[j]||p[j]>upper[j])
            throw std::invalid_argument("Invalid box or initial control");
    }
    auto project=[&](const Vec<P>& x){Vec<P> y{};for(std::size_t j=0;j<P;++j)y[j]=std::clamp(x[j],lower[j],upper[j]);return y;};
    auto projected_gradient=[&](const Vec<P>& x,const Vec<P>& g){return x-project(x-g);};
    BoxOptimum<P> out;
    for(std::size_t k=0;k<=options.max_iterations;++k) {
        auto evaluation=system.evaluate(p);auto pg=projected_gradient(p,evaluation.gradient);double pg_norm=norm(pg);
        out.trace.push_back({p,evaluation.cost,pg_norm,0.});
        if(pg_norm<=options.gradient_tolerance) {
            out.control=p;out.gradient=evaluation.gradient;
            double stationarity=0.,complementarity=0.,primal=0.;
            for(std::size_t j=0;j<P;++j) {
                const double tol=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(lower[j]),std::abs(upper[j])});
                out.active_lower[j]=p[j]-lower[j]<=tol;out.active_upper[j]=upper[j]-p[j]<=tol;
                if(out.active_lower[j])out.lower_multiplier[j]=std::max(0.,evaluation.gradient[j]);
                if(out.active_upper[j])out.upper_multiplier[j]=std::max(0.,-evaluation.gradient[j]);
                stationarity=std::max(stationarity,std::abs(evaluation.gradient[j]-out.lower_multiplier[j]+out.upper_multiplier[j]));
                complementarity=std::max(complementarity,std::max(std::abs(out.lower_multiplier[j]*(p[j]-lower[j])),std::abs(out.upper_multiplier[j]*(upper[j]-p[j]))));
                primal=std::max(primal,std::max(lower[j]-p[j],p[j]-upper[j]));
            }
            out.kkt_residual=std::max({stationarity,complementarity,primal,pg_norm});return out;
        }
        if(k==options.max_iterations)break;
        Mat<P,P> H{};
        for(std::size_t j=0;j<P;++j){Vec<P> unit{};unit[j]=1.;auto col=system.hessian_vector(p,unit);for(std::size_t i=0;i<P;++i)H[i][j]=col[i];}
        Vec<P> step{};Mat<P,P> reduced{};Vec<P> rhs{};std::array<bool,P> fixed{};
        for(std::size_t j=0;j<P;++j){const double tol=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(lower[j]),std::abs(upper[j])});fixed[j]=(p[j]-lower[j]<=tol&&evaluation.gradient[j]>0.)||(upper[j]-p[j]<=tol&&evaluation.gradient[j]<0.);}
        for(std::size_t i=0;i<P;++i)for(std::size_t j=0;j<P;++j)reduced[i][j]=fixed[i]||fixed[j]?0.:H[i][j];
        for(std::size_t j=0;j<P;++j){if(fixed[j])reduced[j][j]=1.;else rhs[j]=-evaluation.gradient[j];}
        try {step=LU<P>(reduced).solve(rhs);if(!(dot(evaluation.gradient,step)<0.))step=-1.*pg;}catch(const std::runtime_error&){step=-1.*pg;}
        bool accepted=false;double alpha=1.;
        for(std::size_t j=0;j<options.max_backtracks;++j){auto trial=project(p+alpha*step);auto displacement=trial-p;
            if(norm(displacement)>0&&dot(evaluation.gradient,displacement)<0.) {
                const double slack=options.cost_roundoff_slack+8*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(evaluation.cost));
                try {if(system.evaluate(trial).cost<=evaluation.cost+1e-4*dot(evaluation.gradient,displacement)+slack){p=trial;accepted=true;break;}}
                catch(const std::domain_error&){/* Backtrack to an admissible state. */}
            }alpha*=.5;
        }
        if(!accepted)throw std::runtime_error("Box optimizer line search failed at iteration "+std::to_string(k)+" projected gradient "+std::to_string(pg_norm));
        out.trace.back().step_length=alpha;
    }
    throw std::runtime_error("Box optimizer failed to converge");
}

struct PathConstraintOptions {
    double feasibility_tolerance=1e-9;
    double complementarity_tolerance=1e-8;
    double initial_penalty=10.;
    double penalty_growth=10.;
    std::size_t max_outer_iterations=20;
};
template<std::size_t P> struct PathConstrainedOptimum {
    Vec<P> control{},gradient{},lower_multiplier{},upper_multiplier{};
    std::array<bool,P> active_lower{},active_upper{};
    std::vector<double> path_multiplier,path_value;
    std::vector<Iteration<P>> trace;
    double objective{},stationarity_residual{},feasibility_residual{},complementarity_residual{},kkt_residual{};
    std::size_t outer_iterations{};
};

// Augmented-Lagrangian solve for nodal path inequalities g_i(x_{n_i},p)<=0,
// combined with the existing simple control bounds. Each path term is inserted
// as a nodal cost, so state dependence is differentiated by the discrete
// adjoint and its Hessian by the incremental adjoint. The constraints are only
// enforced at the selected transcription nodes.
template<std::size_t N,std::size_t P>
PathConstrainedOptimum<P> minimize_path_constrained(const System<N,P>& system,Vec<P> p,
        const Vec<P>& lower,const Vec<P>& upper,OptimizeOptions optimize_options={},
        PathConstraintOptions constraint_options={}) {
    const std::size_t count=system.path_inequality_count();
    if(!count)throw std::invalid_argument("Path-constrained optimizer needs at least one path inequality");
    if(!std::isfinite(constraint_options.feasibility_tolerance)||constraint_options.feasibility_tolerance<=0.||
       !std::isfinite(constraint_options.complementarity_tolerance)||constraint_options.complementarity_tolerance<=0.||
       !std::isfinite(constraint_options.initial_penalty)||constraint_options.initial_penalty<=0.||
       !std::isfinite(constraint_options.penalty_growth)||constraint_options.penalty_growth<=1.||
       !constraint_options.max_outer_iterations)
        throw std::invalid_argument("Path-constraint optimizer options");
    std::vector<double> multiplier(count,0.);double penalty=constraint_options.initial_penalty;
    double previous_feasibility=std::numeric_limits<double>::infinity();
    PathConstrainedOptimum<P> result;
    for(std::size_t outer=0;outer<constraint_options.max_outer_iterations;++outer){
        auto augmented=system.with_path_lagrangian(multiplier,penalty);
        auto inner=minimize_box_constrained(augmented,p,lower,upper,optimize_options);
        p=inner.control;auto trajectory=system.forward(p);
        auto value=system.path_inequality_values(p,trajectory);
        double feasibility=0.;
        std::vector<double> updated(count);
        for(std::size_t i=0;i<count;++i){
            feasibility=std::max(feasibility,std::max(0.,value[i]));
            updated[i]=std::max(0.,multiplier[i]+penalty*value[i]);
        }
        auto lagrangian=system.with_path_lagrangian(updated,0.);
        auto evaluation=lagrangian.evaluate(p);
        Vec<P> projected{};
        for(std::size_t j=0;j<P;++j)
            projected[j]=p[j]-std::clamp(p[j]-evaluation.gradient[j],lower[j],upper[j]);
        double complementarity=0.;
        for(std::size_t i=0;i<count;++i)complementarity=std::max(complementarity,std::abs(updated[i]*value[i]));
        result.control=p;result.gradient=evaluation.gradient;result.path_multiplier=updated;
        result.path_value=value;result.objective=system.cost(p,evaluation.trajectory);
        result.stationarity_residual=norm(projected);result.feasibility_residual=feasibility;
        result.complementarity_residual=complementarity;
        result.kkt_residual=std::max({result.stationarity_residual,feasibility,complementarity});
        result.trace.insert(result.trace.end(),inner.trace.begin(),inner.trace.end());
        result.outer_iterations=outer+1;
        for(std::size_t j=0;j<P;++j){
            const double tol=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(lower[j]),std::abs(upper[j])});
            result.active_lower[j]=p[j]-lower[j]<=tol;result.active_upper[j]=upper[j]-p[j]<=tol;
            if(result.active_lower[j])result.lower_multiplier[j]=std::max(0.,evaluation.gradient[j]);
            if(result.active_upper[j])result.upper_multiplier[j]=std::max(0.,-evaluation.gradient[j]);
            complementarity=std::max(complementarity,std::max(
                std::abs(result.lower_multiplier[j]*(p[j]-lower[j])),
                std::abs(result.upper_multiplier[j]*(upper[j]-p[j]))));
        }
        result.complementarity_residual=complementarity;
        result.kkt_residual=std::max({result.stationarity_residual,result.feasibility_residual,complementarity});
        if(result.stationarity_residual<=optimize_options.gradient_tolerance&&
           feasibility<=constraint_options.feasibility_tolerance&&
           complementarity<=constraint_options.complementarity_tolerance)return result;
        multiplier=std::move(updated);
        if(feasibility>.5*previous_feasibility)penalty*=constraint_options.penalty_growth;
        previous_feasibility=feasibility;
        if(!std::isfinite(penalty)||penalty>1e14)break;
    }
    throw std::runtime_error("Path-constrained optimizer failed to satisfy discrete KKT tolerances; residual="+
                             std::to_string(result.kkt_residual));
}
}
