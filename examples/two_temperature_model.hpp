// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include <theta_opt/theta_opt.hpp>
#include <algorithm>
#include <array>
#include <optional>

namespace thermal {
using namespace theta_opt;
inline constexpr double t0=2.,tf=12.,sigma=.1;
inline constexpr Vec<2> equilibrium{30.,20.},truth{35.,15.};
inline constexpr Mat<2,2> dynamics_matrix{{{-2.,1.},{.25,-.5}}};
inline constexpr std::array<Vec<2>,2> goals{{{1.,0.},{60.,240.}}};
inline constexpr Vec<2> tolerance{.2,60.};
inline constexpr std::array<double,8> elapsed{0.,.5,1.,2.,4.,6.,8.,10.};
inline constexpr std::array<double,8> observations{
35.07773023553762926,30.53933840230146401,28.97603352688368972,
28.95529246360522535,29.36570771254467971,29.77228485594519469,
29.75098075498287997,29.94019717082843357};
inline constexpr NewtonOptions solver_options{2e-13,0.,10,40};
inline constexpr double slow=(-5.+3.605551275463989293)/4.;
inline constexpr double fast=(-5.-3.605551275463989293)/4.;
inline Mat<2,2> slow_projector(){return scaled(1./(slow-fast),add(dynamics_matrix,scaled(-fast,identity<2>())));}
inline Mat<2,2> propagator(double tau){
    auto S=slow_projector(),F=add(identity<2>(),scaled(-1.,S));
    return add(scaled(std::exp(slow*tau),S),scaled(std::exp(fast*tau),F));
}
inline Vec<2> exact(double tau,const Vec<2>& initial){return equilibrium+mv(propagator(tau),initial-equilibrium);}
inline Grid mesh(std::size_t subdivisions,double theta){
    if(!subdivisions)throw std::invalid_argument("Thermal subdivisions");
    std::vector<double> times;
    for(std::size_t j=0;j+1<elapsed.size();++j)for(std::size_t k=0;k<subdivisions;++k)
        times.push_back(t0+elapsed[j]+(elapsed[j+1]-elapsed[j])*static_cast<double>(k)/static_cast<double>(subdivisions));
    times.push_back(tf);return Grid(times,std::vector<double>(times.size()-1,theta));
}
inline std::size_t node_at(const Grid& grid,double offset){
    for(std::size_t n=0;n<grid.time.size();++n)if(std::abs(grid.time[n]-t0-offset)<1e-11)return n;
    throw std::invalid_argument("Thermal observation must be a node");
}
inline double wall_observation(){return exact(1.,truth)[1]+.03;}
inline ScalarFunction<2,2> measurement(std::size_t component,double value){
    auto cost=ScalarFunction<2,2>::zero();
    cost.evaluate=[=](double,const Vec<2>& x,const Vec<2>&){ScalarValue<2,2> v;double r=x[component]-value;v.value=.5*r*r/(sigma*sigma);v.x[component]=r/(sigma*sigma);return v;};
    cost.directional_gradient=[=](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>&){ScalarDirection<2,2> v;v.x[component]=s[component]/(sigma*sigma);return v;};return cost;
}
inline Problem<2,2> problem(const Grid& grid,bool wall_sensor=false){
    Problem<2,2> p;
    p.initial.value=[](const Vec<2>& xi){return xi;};
    p.initial.jacobian=[](const Vec<2>&){return identity<2>();};
    p.initial.directional_jacobian=[](const Vec<2>&,const Vec<2>&){return Mat<2,2>{};};
    // M^-1(K T + forcing), with M=diag(1,4). This is algebraically
    // equivalent to the original constant-mass energy balance theta step.
    p.dynamics.value=[](double,const Vec<2>& x,const Vec<2>&){return mv(dynamics_matrix,x-equilibrium);};
    p.dynamics.derivatives=[](double,const Vec<2>&,const Vec<2>&){Derivatives<2,2> d;d.x=dynamics_matrix;return d;};
    p.dynamics.directional_derivatives=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>&,const Vec<2>&){return Derivatives<2,2>{};};
    for(std::size_t j=0;j<elapsed.size();++j)p.nodal.emplace_back(node_at(grid,elapsed[j]),measurement(0,observations[j]));
    if(wall_sensor)p.nodal.emplace_back(node_at(grid,1.),measurement(1,wall_observation()));
    return p;
}
inline System<2,2> system(const Grid& grid,bool wall=false){return System<2,2>(problem(grid,wall),grid,solver_options);}
inline Optimum<2> fit(const System<2,2>& sys){return minimize(sys,equilibrium,OptimizeOptions{1e-8,10,40});}
inline Mat<2,2> inverse(const Mat<2,2>& a){Mat<2,2> result{};LU<2> lu(a);for(std::size_t j=0;j<2;++j){Vec<2> unit{};unit[j]=1.;auto col=lu.solve(unit);for(std::size_t i=0;i<2;++i)result[i][j]=col[i];}return result;}
inline Mat<2,2> outer(const Vec<2>& v){return {{{v[0]*v[0],v[0]*v[1]},{v[1]*v[0],v[1]*v[1]}}};}
struct Reference {Vec<2> fit{};Mat<2,2> covariance{};};
inline Reference analytical_fit(bool wall=false){
    Mat<2,2> information{};Vec<2> rhs{};
    for(std::size_t j=0;j<elapsed.size();++j){auto h=propagator(elapsed[j])[0];information=add(information,scaled(1./(sigma*sigma),outer(h)));rhs=rhs+((observations[j]-equilibrium[0])/(sigma*sigma))*h;}
    if(wall){auto h=propagator(1.)[1];information=add(information,scaled(1./(sigma*sigma),outer(h)));rhs=rhs+((wall_observation()-equilibrium[1])/(sigma*sigma))*h;}
    return {equilibrium+LU<2>(information).solve(rhs),inverse(information)};
}
template<class Function> double root(Function f,double lo,double hi){
    double a=f(lo),b=f(hi);if(a==0)return lo;if(b==0)return hi;if(a*b>0)throw std::runtime_error("Thermal root not bracketed");
    for(int k=0;k<100;++k){double mid=.5*(lo+hi),v=f(mid);if(mid==lo||mid==hi||v==0)return mid;if(a*v<=0)hi=mid;else{lo=mid;a=v;}}return .5*(lo+hi);
}
inline double ellipse_bound(double tau,const Mat<2,2>& covariance,const Vec<2>& goal){
    auto g=mv(transpose(propagator(tau)),goal);return std::sqrt(-2*std::log(.05))*std::sqrt(local_goal_variance(g,covariance));
}
inline std::optional<double> analytical_qualification(const Mat<2,2>& covariance,std::size_t goal){
    auto S=slow_projector(),F=add(identity<2>(),scaled(-1.,S));
    auto s=mv(transpose(S),goals[goal]),f=mv(transpose(F),goals[goal]);
    double p=dot(s,mv(covariance,s)),q=2*dot(s,mv(covariance,f)),r=dot(f,mv(covariance,f));
    double aa=2*fast*r,bb=(slow+fast)*q,cc=2*slow*p;
    std::vector<double> critical{0.,10.};
    auto add_root=[&](double z){if(z>0){double t=std::log(z)/(fast-slow);if(t>0&&t<10)critical.push_back(t);}};
    if(std::abs(aa)<1e-25){if(std::abs(bb)>1e-25)add_root(-cc/bb);}
    else{double disc=bb*bb-4*aa*cc;if(disc>=0){add_root((-bb-std::sqrt(disc))/(2*aa));add_root((-bb+std::sqrt(disc))/(2*aa));}}
    std::sort(critical.begin(),critical.end());auto e=[&](double t){return ellipse_bound(t,covariance,goals[goal])-tolerance[goal];};
    if(e(10)>0)return std::nullopt;
    double last=0.;
    for(std::size_t k=0;k+1<critical.size();++k)if(e(critical[k])*e(critical[k+1])<0)last=root(e,critical[k],critical[k+1]);
    return last;
}

// Convex confidence set: physical initial domain intersected with affine
// observation strips. Polygon clipping is performed entirely in C++20.
using Polygon=std::vector<Vec<2>>;
inline Polygon prior_polygon(){return {{0.,0.},{60.,0.},{60.,60.},{0.,60.}};}
inline Polygon clip(Polygon polygon,const Vec<2>& normal,double bound){
    Polygon out;if(polygon.empty())return out;
    Vec<2> previous=polygon.back();double fp=dot(normal,previous)-bound;
    for(auto current:polygon){double fc=dot(normal,current)-bound;bool inside=fc<=0,was_inside=fp<=0;
        if(inside!=was_inside){double alpha=fp/(fp-fc);out.push_back(previous+alpha*(current-previous));}
        if(inside)out.push_back(current);
        previous=current;fp=fc;
    }return out;
}
inline bool contains(const Polygon& p,const Vec<2>& x,double tol=1e-8){
    if(p.size()<3)return false;
    for(std::size_t i=0;i<p.size();++i){auto a=p[i],b=p[(i+1)%p.size()];if((b[0]-a[0])*(x[1]-a[1])-(b[1]-a[1])*(x[0]-a[0])< -tol)return false;}return true;
}
struct AffineTrajectory {std::vector<Vec<2>> base;std::array<std::vector<Vec<2>>,2> tangent;};
inline AffineTrajectory affine(const System<2,2>& sys){
    auto tr=sys.forward(equilibrium);return {tr.state,{sys.tangent(equilibrium,tr,{1.,0.}),sys.tangent(equilibrium,tr,{0.,1.})}};
}
inline Vec<2> affine_goal_gradient(const AffineTrajectory& path,std::size_t node,const Vec<2>& goal){return {dot(goal,path.tangent[0][node]),dot(goal,path.tangent[1][node])};}
inline std::vector<double> physical_allowance(const System<2,2>& sys){
    const auto& grid=sys.grid();std::vector<double> D(grid.steps(),0.),E(grid.time.size(),0.);
    for(auto corner:prior_polygon()){auto tr=sys.forward(corner);
        for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n],th=grid.theta[n];
            double defect=std::max(th,1-th)*norm(mv(dynamics_matrix,tr.state[n+1]-tr.state[n]))+solver_options.absolute_tolerance/h;
            D[n]=std::max(D[n],defect);}}
    // mu_infinity(A)=-0.25, so all model errors contract at least at this rate.
    for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n];E[n+1]=std::exp(-.25*h)*E[n]+(-std::expm1(-.25*h)/.25)*D[n];}return E;
}
inline Polygon confidence_polygon(const Grid& grid,const AffineTrajectory& path,double radius,
                                  const std::vector<double>& allowance,bool wall=false){
    auto polygon=prior_polygon();
    auto observation=[&](std::size_t node,std::size_t component,double value){
        Vec<2> h{path.tangent[0][node][component],path.tangent[1][node][component]};
        double offset=path.base[node][component]-dot(h,equilibrium),bound=radius+allowance[node];
        polygon=clip(std::move(polygon),h,value+bound-offset);
        polygon=clip(std::move(polygon),-1.*h,-value+bound+offset);
    };
    for(std::size_t j=0;j<elapsed.size();++j)observation(node_at(grid,elapsed[j]),0,observations[j]);
    if(wall)observation(node_at(grid,1.),1,wall_observation());
    if(polygon.empty())throw std::runtime_error("Empty thermal confidence polygon");
    return polygon;
}
inline Polygon analytical_polygon(const Grid& grid,double radius,bool wall=false){
    AffineTrajectory path;
    for(double time:grid.time){auto E=propagator(time-t0);path.base.push_back(equilibrium);path.tangent[0].push_back({E[0][0],E[1][0]});path.tangent[1].push_back({E[0][1],E[1][1]});}
    return confidence_polygon(grid,path,radius,std::vector<double>(grid.time.size(),0.),wall);
}
inline Mat<2,2> analytical_running_hessian(){
    auto S=slow_projector(),F=add(identity<2>(),scaled(-1.,S));
    auto integral=[](double rate){return std::expm1(10.*rate)/rate;};
    Mat<2,2> q=add(scaled(integral(2*slow),outer(S[0])),scaled(integral(2*fast),outer(F[0])));
    for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j)q[i][j]+=integral(slow+fast)*(S[0][i]*F[0][j]+F[0][i]*S[0][j]);
    return q;
}
inline Problem<2,2> running_problem(const Grid& grid){
    auto p=problem(grid);p.nodal.clear();
    p.running.evaluate=[](double,const Vec<2>& x,const Vec<2>&){double r=x[0]-equilibrium[0];return ScalarValue<2,2>{.5*r*r,{r,0.},{0.,0.}};};
    p.running.directional_gradient=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>&){return ScalarDirection<2,2>{{s[0],0.},{0.,0.}};};return p;
}
}
