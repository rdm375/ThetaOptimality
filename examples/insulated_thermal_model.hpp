// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include <theta_opt/theta_opt.hpp>
namespace insulated_thermal {
using namespace theta_opt;
inline constexpr double t0=2.,tf=10.,horizon=tf-t0,capacity1=2.,capacity2=3.,conductance=.6;
inline constexpr double w1=capacity1/(capacity1+capacity2),w2=1-w1,rate=conductance*(1/capacity1+1/capacity2),sigma=.1;
inline constexpr Mat<2,2> A{{{-conductance/capacity1,conductance/capacity1},{conductance/capacity2,-conductance/capacity2}}};
inline constexpr Vec<2> truth{32.,20.},probe{30.,19.},seed{24.,24.},difference{1.,-1.},mean{w1,w2},common{1.,1.};
inline Mat<2,2> outer(const Vec<2>& x){return {{{x[0]*x[0],x[0]*x[1]},{x[0]*x[1],x[1]*x[1]}}};}
inline Mat<2,2> propagator(double tau){double e=std::exp(-rate*tau);return {{{w1+w2*e,w2*(1-e)},{w1*(1-e),w2+w1*e}}};}
inline Vec<2> exact(double tau,const Vec<2>& x){double m=dot(mean,x),d=dot(difference,x)*std::exp(-rate*tau);return {m+w2*d,m-w1*d};}
inline double energy(const Vec<2>& x){return capacity1*x[0]+capacity2*x[1];} // kJ relative to 0 degrees C
struct Observation {double tau;Vec<2> goal;double value,scale;};
inline std::vector<Observation> observations(bool added=false,bool noisy=true){std::vector<Observation> r;std::array<double,4> times{0.,1.,2.,4.},noise{-.2,.3,-.4,.1};for(std::size_t j=0;j<times.size();++j)r.push_back({times[j],difference,dot(difference,exact(times[j],truth))+(noisy?sigma*noise[j]:0.),sigma});if(added)r.push_back({1.,{1.,0.},exact(1.,truth)[0]+(noisy?.01:0.),.05});return r;}
inline std::size_t node_at(const Grid& grid,double tau){for(std::size_t n=0;n<grid.time.size();++n)if(std::abs(grid.time[n]-t0-tau)<1e-12)return n;throw std::invalid_argument("Thermal observations must be grid nodes");}
inline Vec<2> observation_row(const Observation& o){return mv(transpose(propagator(o.tau)),o.goal);}
inline Mat<2,2> information(bool added){Mat<2,2> H{};for(auto o:observations(added))H=add(H,scaled(1/(o.scale*o.scale),outer(observation_row(o))));return H;}
// Independent modal normal equations: the absolute datum is matched exactly
// since the difference data leave the mean unconstrained.
inline Vec<2> analytical_fit(bool added=false,Vec<2> start=seed,bool noisy=true){double numerator=0.,denominator=0.;auto data=observations(added,noisy);for(std::size_t j=0;j<4;++j){double e=std::exp(-rate*data[j].tau);numerator+=e*data[j].value/(sigma*sigma);denominator+=e*e/(sigma*sigma);}double d=numerator/denominator;if(!added){double s=.5*(start[0]+start[1]);return {s+.5*d,s-.5*d};}auto o=data.back();double m=o.value-w2*d*std::exp(-rate*o.tau);return {m+w2*d,m-w1*d};}
inline ScalarFunction<2,2> measurement(Observation o){auto c=ScalarFunction<2,2>::zero();c.evaluate=[o](double,const Vec<2>& x,const Vec<2>&){double r=dot(o.goal,x)-o.value;return ScalarValue<2,2>{.5*r*r/(o.scale*o.scale),(r/(o.scale*o.scale))*o.goal,{}};};c.directional_gradient=[o](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>&){return ScalarDirection<2,2>{(dot(o.goal,s)/(o.scale*o.scale))*o.goal,{}};};return c;}
inline System<2,2> system(std::size_t steps,double theta,bool added=false,bool running=false,bool noisy=true){Problem<2,2> p;auto grid=Grid::uniform(t0,tf,steps,theta);p.initial.value=[](auto x){return x;};p.initial.jacobian=[](auto){return identity<2>();};p.initial.directional_jacobian=[](auto,auto){return Mat<2,2>{};};p.dynamics.value=[](double,auto x,auto){return mv(A,x);};p.dynamics.derivatives=[](double,auto,auto){Derivatives<2,2> r;r.x=A;return r;};p.dynamics.directional_derivatives=[](double,auto,auto,auto,auto){return Derivatives<2,2>{};};
 if(running){p.running.evaluate=[](double,auto x,auto){double d=dot(difference,x);return ScalarValue<2,2>{.5*d*d,d*difference,{}};};p.running.directional_gradient=[](double,auto,auto,auto s,auto){return ScalarDirection<2,2>{dot(difference,s)*difference,{}};};}else for(auto o:observations(added,noisy))p.nodal.emplace_back(node_at(grid,o.tau),measurement(o));return System<2,2>(p,grid);}
inline QuadraticOptimum2 fit(const System<2,2>& sys,Vec<2> start=seed){return minimize_quadratic_psd_2d(sys,start);}
inline double exact_cost(Vec<2> x,bool added){double J=0.;for(auto o:observations(added)){double r=dot(o.goal,exact(o.tau,x))-o.value;J+=.5*r*r/(o.scale*o.scale);}return J;}
inline Vec<2> exact_gradient(Vec<2> x,bool added){Vec<2> g{};for(auto o:observations(added))g=g+((dot(o.goal,exact(o.tau,x))-o.value)/(o.scale*o.scale))*observation_row(o);return g;}
inline Vec<2> nodal_costate(double tau,Vec<2> x,bool added){Vec<2> l{};for(auto o:observations(added))if(o.tau+1e-12>=tau)l=l+((dot(o.goal,exact(o.tau,x))-o.value)/(o.scale*o.scale))*mv(transpose(propagator(o.tau-tau)),o.goal);return l;}
inline Mat<2,2> gramian(double duration){double integral=-std::expm1(-2*rate*duration)/(2*rate);return scaled(integral,outer(difference));}
inline Vec<2> running_costate(double tau,Vec<2> x){return mv(gramian(horizon-tau),exact(tau,x));}
}
