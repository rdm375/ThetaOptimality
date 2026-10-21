// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include "decay_model.hpp"
#include <theta_opt/qualification.hpp>
#include <optional>
#include <string>

namespace decay_qualification {
using namespace theta_opt;
inline constexpr std::array<double,9> elapsed{0.,.25,.5,1.,2.,4.,6.,8.,10.};
// Retained synthetic records from the Python benchmark, embedded to avoid
// implementation-dependent C++ random-distribution sequences.
inline constexpr std::array<double,9> gaussian{
3.023319070661288688,2.475663321902065839,1.997103570151420593,
1.480409391827697707,.7845566751614213574,.2861883603767760853,
.06383756281227698226,.03825846087639751120,.009863706548190513679};
inline constexpr std::array<double,9> bounded{
3.026954999610935726,2.463975848002926572,2.090975787271590836,
1.463207653679871045,.7633181737440541959,.2943523657622753031,
.1109271660351506816,.03968597501298244679,.002966809496667933282};
inline constexpr double truth=3.,sigma=.03,delta=.04;
inline constexpr Vec<2> tolerance{.02,.025};
inline constexpr NewtonOptions solver_options{2e-14,0.,40,40};
struct Exposure {double elapsed=2.,value=decay::exact(2.,truth).state[1]+.002,bound=.005;};
struct Interval {double lower{},upper{};};

inline Grid mesh(std::size_t subdivisions,double theta){
    std::vector<double> times;
    for(std::size_t j=0;j+1<elapsed.size();++j)
        for(std::size_t k=0;k<subdivisions;++k)
            times.push_back(decay::t0+elapsed[j]+(elapsed[j+1]-elapsed[j])*static_cast<double>(k)/static_cast<double>(subdivisions));
    times.push_back(decay::tf);return Grid(times,std::vector<double>(times.size()-1,theta));
}
inline std::size_t node_at(const Grid& grid,double offset){
    for(std::size_t n=0;n<grid.time.size();++n)
        if(std::abs(grid.time[n]-decay::t0-offset)<1e-11)return n;
    throw std::invalid_argument("Observation must coincide with a node");
}
inline ScalarFunction<2,1> measurement(std::size_t component,double value,double sd){
    auto cost=ScalarFunction<2,1>::zero();
    cost.evaluate=[=](double,const Vec<2>& x,const Vec<1>&){
        double residual=x[component]-value;ScalarValue<2,1> v;
        v.value=.5*residual*residual/(sd*sd);v.x[component]=residual/(sd*sd);return v;
    };
    cost.directional_gradient=[=](double,const Vec<2>&,const Vec<1>&,const Vec<2>& s,const Vec<1>&){
        ScalarDirection<2,1> v;v.x[component]=s[component]/(sd*sd);return v;
    };return cost;
}
inline Problem<2,1> assimilation(const Grid& grid,const std::array<double,9>& data,
                               std::optional<Exposure> extra=std::nullopt){
    auto problem=decay::problem();
    problem.running=ScalarFunction<2,1>::zero();problem.terminal=ScalarFunction<2,1>::zero();problem.penalty=Penalty<1>::zero();
    for(std::size_t j=0;j<elapsed.size();++j)
        problem.nodal.emplace_back(node_at(grid,elapsed[j]),measurement(0,data[j],sigma));
    if(extra)problem.nodal.emplace_back(node_at(grid,extra->elapsed),measurement(1,extra->value,extra->bound));
    return problem;
}
inline System<2,1> system(const Grid& grid,const std::array<double,9>& data,
                        std::optional<Exposure> extra=std::nullopt){
    return System<2,1>(assimilation(grid,data,extra),grid,solver_options);
}
inline Optimum<1> fit(const System<2,1>& system){
    return minimize(system,Vec<1>{1.},OptimizeOptions{1e-8,40,40},
                    [](const Vec<1>& xi){return xi[0]>0;});
}
inline std::array<double,3> analytical_objective(double xi,const std::array<double,9>& data,
                                               std::optional<Exposure> extra=std::nullopt){
    std::array<double,3> out{};
    for(std::size_t j=0;j<elapsed.size();++j){auto ex=decay::exact(elapsed[j],xi);double r=ex.state[0]-data[j];
        out[0]+=.5*r*r/(sigma*sigma);out[1]+=r*ex.tangent[0]/(sigma*sigma);
        out[2]+=(ex.tangent[0]*ex.tangent[0]+r*ex.second[0])/(sigma*sigma);
    }
    if(extra){auto ex=decay::exact(extra->elapsed,xi);double r=ex.state[1]-extra->value,variance=extra->bound*extra->bound;
        out[0]+=.5*r*r/variance;out[1]+=r*ex.tangent[1]/variance;out[2]+=(ex.tangent[1]*ex.tangent[1]+r*ex.second[1])/variance;}
    return out;
}
inline double analytical_fit(const std::array<double,9>& data,std::optional<Exposure> extra=std::nullopt){
    return decay::root([&](double xi){return analytical_objective(xi,data,extra)[1];});
}

// Intersect monotone observation preimages. Forward calls below use the C++
// library; the same routine also accepts the independent analytical oracle.
template<class Forward>
std::optional<Interval> consistent_interval(Forward forward,const Grid& grid,
    const std::array<double,9>& data,const std::array<double,9>& bounds,
    std::optional<Exposure> extra=std::nullopt){
    Interval interval{std::max(1e-10,data[0]-bounds[0]),data[0]+bounds[0]};
    auto intersect=[&](std::size_t node,std::size_t component,double value,double radius){
        auto f=[&](double xi){return forward(xi)[node][component];};
        double a=f(interval.lower),b=f(interval.upper),lo=value-radius,hi=value+radius;
        if(lo>b+1e-12||hi<a-1e-12)return false;
        if(lo>a)interval.lower=decay::root([&](double xi){return f(xi)-lo;},interval.lower,interval.upper);
        if(hi<b)interval.upper=decay::root([&](double xi){return f(xi)-hi;},interval.lower,interval.upper);
        return interval.lower<=interval.upper;
    };
    for(std::size_t j=0;j<elapsed.size();++j)if(!intersect(node_at(grid,elapsed[j]),0,data[j],bounds[j]))return std::nullopt;
    if(extra&&!intersect(node_at(grid,extra->elapsed),1,extra->value,extra->bound))return std::nullopt;
    return interval;
}
inline std::vector<Vec<2>> exact_path(const Grid& grid,double xi){
    std::vector<Vec<2>> x;for(double t:grid.time)x.push_back(decay::exact(t-decay::t0,xi).state);return x;
}
inline std::array<double,9> constant_bounds(double bound=delta){std::array<double,9> v;v.fill(bound);return v;}
inline Interval exact_interval(std::optional<Exposure> extra=std::nullopt){
    const auto grid=mesh(1,.5);return consistent_interval([&](double xi){return exact_path(grid,xi);},grid,bounded,constant_bounds(),extra).value();
}
inline std::vector<Vec<2>> envelope(const std::vector<Vec<2>>& nominal,
                                       const std::vector<Vec<2>>& lower,const std::vector<Vec<2>>& upper){
    std::vector<Vec<2>> out(nominal.size());
    for(std::size_t n=0;n<out.size();++n)for(std::size_t j=0;j<2;++j)
        out[n][j]=std::max(std::abs(nominal[n][j]-lower[n][j]),std::abs(nominal[n][j]-upper[n][j]));
    return out;
}
inline std::array<std::optional<double>,2> qualify(const Grid& grid,const std::vector<Vec<2>>& error){
    std::array<std::optional<double>,2> out;
    for(std::size_t j=0;j<2;++j){std::vector<double> ratio;for(auto v:error)ratio.push_back(v[j]/tolerance[j]);
        out[j]=persistent_node_qualification(grid.time,ratio);if(out[j])*out[j]-=decay::t0;}
    return out;
}
inline double output(std::optional<double> t){return t.value_or(std::numeric_limits<double>::quiet_NaN());}

// Decay-specific analytic defect bound, computed from library trajectories.
// It is uniform on the initial interval; it is NOT an exact-solution oracle.
// Floating-point evaluation is not directed-rounding interval certification.
inline std::vector<Vec<2>> physical_allowance(const System<2,1>& system,Interval prior){
    const auto& grid=system.grid();auto high=system.forward({prior.upper});auto low=system.forward({prior.lower});
    std::vector<Vec<2>> allowance(grid.time.size());double E=0.,integral=0.,quadrature=0.;
    // Uniform Newton residual allowance on this positive, decreasing branch.
    const double residual_limit=solver_options.absolute_tolerance;
    for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n],th=grid.theta[n];
        double x0=high.state[n][0],x1=high.state[n+1][0];
        if(1-h*(1-th)*(decay::a+2*decay::b*x0)<0||x1>x0||low.state[n+1][0]<=0)
            throw std::runtime_error("Refine grid: positive monotone branch required for bound");
        double g0=decay::a*x0+decay::b*x0*x0,g1=decay::a*x1+decay::b*x1*x1;
        double D=std::max(th,1-th)*std::abs(g0-g1)+residual_limit/h;
        double decay_factor=std::exp(-decay::a*h),phi=-std::expm1(-decay::a*h)/decay::a;
        integral+=E*phi+D/decay::a*(h-phi);E=decay_factor*E+D*phi;
        quadrature+=h*std::abs(th-.5)*std::abs(x1-x0)+residual_limit;
        allowance[n+1]={E,integral+quadrature};
    }return allowance;
}
inline std::vector<Vec<2>> curvature_allowance(const System<2,1>& system,Interval hull){
    auto lo=system.forward({hull.lower}),hi=system.forward({hull.upper});const auto& grid=system.grid();
    std::vector<Vec<2>> bounds(grid.time.size());double S=1.,Z=0.,Wq=0.;
    for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n],th=grid.theta[n];
        double Cmin=1-h*(1-th)*(decay::a+2*decay::b*hi.state[n][0]);
        double Cmax=1-h*(1-th)*(decay::a+2*decay::b*lo.state[n][0]);
        double Bmin=1+h*th*(decay::a+2*decay::b*lo.state[n+1][0]);
        if(Cmin<0)throw std::runtime_error("Nonmonotone branch in curvature bound");
        double T=Cmax/Bmin,Sn=T*S,Zn=T*Z+2*h*decay::b/Bmin*((1-th)*S*S+th*Sn*Sn);
        Wq+=h*((1-th)*Z+th*Zn);bounds[n+1]={Zn,Wq};S=Sn;Z=Zn;
    }return bounds;
}
}
