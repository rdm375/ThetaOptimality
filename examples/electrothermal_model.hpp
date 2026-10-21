// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include <theta_opt/theta_opt.hpp>
#include <complex>
namespace electrothermal {
using namespace theta_opt;
inline constexpr double t0=2.,tf=10.,horizon=tf-t0,R=1.,L=1.,C=1.,heat_capacity=2.,cooling=.4,decay=cooling/heat_capacity,beta=R/heat_capacity,gamma=R/(2*L),ambient=20.,sigma=.01;
inline const double omega=std::sqrt(1/(L*C)-gamma*gamma);
inline constexpr Vec<2> truth{1.2,.7},probe{.9,.5};
inline constexpr Mat<2,2> A{{{0.,-1/C},{1/L,-R/L}}};
inline constexpr Box2 prior{{-1.6,-1.4},{1.6,1.4}};
inline Mat<2,2> outer(Vec<2> a,Vec<2> b){return {{{a[0]*b[0],a[0]*b[1]},{a[1]*b[0],a[1]*b[1]}}};}
inline Mat<2,2> propagator(double tau){double e=std::exp(-gamma*tau),s=std::sin(omega*tau)/omega,c=std::cos(omega*tau);return {{{e*(c+gamma*s),-e*s/C},{e*s/L,e*(c-gamma*s)}}};}
inline Mat<2,2> temperature_matrix(double tau){double z=decay-2*gamma;double f0=std::expm1(z*tau)/z;std::complex<double> w(z,2*omega);auto f=(std::exp(w*tau)-1.)/w;Vec<2> a{0.,1.},b{1/(L*omega),-gamma/omega};return scaled(beta*std::exp(-decay*tau),add(add(scaled(.5*(f0+f.real()),outer(a,a)),scaled(.5*(f0-f.real()),outer(b,b))),scaled(.5*f.imag(),add(outer(a,b),outer(b,a)))));}
inline Vec<3> exact(double tau,Vec<2> xi){auto e=mv(propagator(tau),xi);return {e[0],e[1],dot(xi,mv(temperature_matrix(tau),xi))};}
inline Mat<3,2> exact_jacobian(double tau,Vec<2> xi){auto F=propagator(tau);return {F[0],F[1],2.*mv(temperature_matrix(tau),xi)};}
inline Mat<3,3> state_flow_jacobian(double tau,Vec<3> x){auto F=propagator(tau);auto g=2.*mv(temperature_matrix(tau),Vec<2>{x[0],x[1]});return {{{F[0][0],F[0][1],0.},{F[1][0],F[1][1],0.},{g[0],g[1],std::exp(-decay*tau)}}};}
struct Observation {double tau;std::size_t component;double value,scale;};
// design: 0 temperature-only; 1 informative signed voltage; 2 weak signed voltage.
inline std::vector<Observation> observations(int design=0,bool noisy=true){std::vector<Observation> r;std::array<double,6> times{.5,1.,2.,3.,5.,7.},noise{.1,-.15,.2,-.1,.12,-.08};for(std::size_t j=0;j<times.size();++j)r.push_back({times[j],2,exact(times[j],truth)[2]+(noisy?sigma*noise[j]:0.),sigma});if(design)r.push_back({0.,0,truth[0]+(noisy?.002:0.),design==1?.02:3.});return r;}
inline std::size_t node_at(const Grid& grid,double tau){for(std::size_t n=0;n<grid.time.size();++n)if(std::abs(grid.time[n]-t0-tau)<1e-12)return n;throw std::invalid_argument("Electrothermal observation must be a node");}
inline System<3,2> system(std::size_t steps,double theta,int design=0,bool running=false){Problem<3,2> p;auto grid=Grid::uniform(t0,tf,steps,theta);p.initial.value=[](auto xi){return Vec<3>{xi[0],xi[1],0.};};p.initial.jacobian=[](auto){return Mat<3,2>{{{1.,0.},{0.,1.},{0.,0.}}};};p.initial.directional_jacobian=[](auto,auto){return Mat<3,2>{};};p.dynamics.value=[](double,auto x,auto){return Vec<3>{-x[1]/C,(x[0]-R*x[1])/L,beta*x[1]*x[1]-decay*x[2]};};p.dynamics.derivatives=[](double,auto x,auto){Derivatives<3,2> d;d.x={{{0.,-1/C,0.},{1/L,-R/L,0.},{0.,2*beta*x[1],-decay}}};return d;};p.dynamics.directional_derivatives=[](double,auto,auto,auto s,auto){Derivatives<3,2> d;d.x[2][1]=2*beta*s[1];return d;};
 if(running){p.running.evaluate=[](double,auto x,auto){return ScalarValue<3,2>{R*x[1]*x[1],{0.,2*R*x[1],0.},{}};};p.running.directional_gradient=[](double,auto,auto,auto s,auto){return ScalarDirection<3,2>{{0.,2*R*s[1],0.},{}};};}
 else {for(auto o:observations(design)){auto cost=ScalarFunction<3,2>::zero();cost.evaluate=[o](double,auto x,auto){ScalarValue<3,2> v;double r=x[o.component]-o.value;v.value=.5*r*r/(o.scale*o.scale);v.x[o.component]=r/(o.scale*o.scale);return v;};cost.directional_gradient=[o](double,auto,auto,auto s,auto){ScalarDirection<3,2> d;d.x[o.component]=s[o.component]/(o.scale*o.scale);return d;};p.nodal.emplace_back(node_at(grid,o.tau),cost);}}
 return System<3,2>(p,grid);}
struct ReferenceCost {double cost{};Vec<2> gradient{};Mat<2,2> hessian{},gauss_newton{};};
inline ReferenceCost reference_cost(Vec<2> xi,int design=0){ReferenceCost r;for(auto o:observations(design)){auto state=exact(o.tau,xi);auto g=exact_jacobian(o.tau,xi)[o.component];double residual=state[o.component]-o.value,weight=1/(o.scale*o.scale);r.cost+=.5*weight*residual*residual;r.gradient=r.gradient+(weight*residual)*g;auto GN=scaled(weight,outer(g,g));r.gauss_newton=add(r.gauss_newton,GN);r.hessian=add(r.hessian,GN);if(o.component==2)r.hessian=add(r.hessian,scaled(2*weight*residual,temperature_matrix(o.tau)));}return r;}
// Independent Newton solve of closed analytical objective, never a forward theta solve.
inline Vec<2> analytical_fit(Vec<2> xi,int design=0){for(int k=0;k<50;++k){auto r=reference_cost(xi,design);if(norm(r.gradient)<1e-8)return xi;auto step=LU<2>(r.hessian).solve(-1.*r.gradient);if(dot(step,r.gradient)>=0)step=-1.*r.gradient;double alpha=1.;bool accepted=false;for(int j=0;j<50;++j){auto candidate=xi+alpha*step;if(reference_cost(candidate,design).cost<=r.cost+1e-4*alpha*dot(r.gradient,step)+1e-13){xi=candidate;accepted=true;break;}alpha*=.5;}if(!accepted)throw std::runtime_error("Analytical reference optimization failed");}throw std::runtime_error("Analytical optimum not stationary");}
inline Optimum<2> fit(const System<3,2>& sys,Vec<2> seed){return minimize(sys,seed,OptimizeOptions{1e-7,60,50,1e-12});}
inline Vec<3> nodal_costate(double tau,Vec<2> xi,int design=0){Vec<3> l{};auto current=exact(tau,xi);for(auto o:observations(design))if(o.tau+1e-12>=tau){Vec<3> source{};source[o.component]=(exact(o.tau,xi)[o.component]-o.value)/(o.scale*o.scale);l=l+mv(transpose(state_flow_jacobian(o.tau-tau,current)),source);}return l;}
inline double energy(Vec<3> x){return .5*(C*x[0]*x[0]+L*x[1]*x[1]);}
inline Mat<2,2> energy_matrix(double tau){auto F=propagator(tau);return add(scaled(.5*C,outer(F[0],F[0])),scaled(.5*L,outer(F[1],F[1])));}
inline Mat<2,2> dissipated_matrix(double duration){return add(energy_matrix(0.),scaled(-1.,energy_matrix(duration)));}
inline Vec<3> running_costate(double tau,Vec<2> xi){auto x=exact(tau,xi);auto g=2.*mv(dissipated_matrix(horizon-tau),Vec<2>{x[0],x[1]});return {g[0],g[1],0.};}
}
