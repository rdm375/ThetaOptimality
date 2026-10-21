// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include <theta_opt/theta_opt.hpp>

namespace two_tank {
using namespace theta_opt;
inline constexpr double t0=2.,tf=12.,horizon=tf-t0,flow=1.,fast=4.,slow=.2;
inline constexpr double volume1=flow/fast,volume2=flow/slow,sigma=.003;
inline constexpr double mixing=fast/(fast-slow);
inline constexpr Vec<2> truth{.15,1.},scales{1.,1.}; // mg/L for both states
inline constexpr Mat<2,2> A{{{-fast,fast},{0.,-slow}}};
inline constexpr std::array<double,4> elapsed{0.,6.,8.,9.},noise{-.2,.5,-.4,.3};
inline Mat<2,2> propagator(double tau){double a=std::exp(-fast*tau),d=std::exp(-slow*tau);return {{{a,mixing*(d-a)},{0.,d}}};}
inline Vec<2> exact(double tau,const Vec<2>& xi){return mv(propagator(tau),xi);}
inline double integral(double rate,double duration){return -std::expm1(-rate*duration)/rate;}
inline Mat<2,2> gramian(double duration){double a=integral(2*fast,duration),d=integral(2*slow,duration),ad=integral(fast+slow,duration);return {{{a,mixing*(ad-a)},{mixing*(ad-a),mixing*mixing*(a+d-2*ad)}}};}
inline Mat<2,2> outer(const Vec<2>& v){return {{{v[0]*v[0],v[0]*v[1]},{v[1]*v[0],v[1]*v[1]}}};}
inline Mat<2,2> inverse(const Mat<2,2>& H){Mat<2,2> result{};for(int j=0;j<2;++j){Vec<2> e{};e[j]=1.;auto x=LU<2>(H).solve(e);for(int i=0;i<2;++i)result[i][j]=x[i];}return result;}
struct Observation {double tau;std::size_t component;double value,scale;};
inline std::vector<Observation> observations(bool added=false,bool noisy=true){std::vector<Observation> result;for(std::size_t j=0;j<elapsed.size();++j)result.push_back({elapsed[j],0,exact(elapsed[j],truth)[0]+(noisy?sigma*noise[j]:0.),sigma});if(added)result.push_back({.5,1,exact(.5,truth)[1]+(noisy?.0002:0.),.001});return result;}
inline double observation(std::size_t j,bool noisy=true){return observations(false,noisy)[j].value;}
inline double inventory(const Vec<2>& x){return volume1*x[0]+volume2*x[1];}
inline Mat<2,2> information(bool added=false){Mat<2,2> H{};for(auto o:observations(added)){auto h=propagator(o.tau)[o.component];H=add(H,scaled(1/(o.scale*o.scale),outer(h)));}return H;}
inline Vec<2> analytical_fit(bool added=false,bool noisy=true){Vec<2> b{};for(auto o:observations(added,noisy))b=b+(o.value/(o.scale*o.scale))*propagator(o.tau)[o.component];return LU<2>(information(added)).solve(b);}
inline std::size_t node_at(const Grid& grid,double tau){for(std::size_t n=0;n<grid.time.size();++n)if(std::abs(grid.time[n]-t0-tau)<1e-12)return n;throw std::invalid_argument("Two-tank observation must be a grid node");}
inline ScalarFunction<2,2> measurement(const Observation& o){auto cost=ScalarFunction<2,2>::zero();cost.evaluate=[o](double,const Vec<2>& x,const Vec<2>&){double r=x[o.component]-o.value;ScalarValue<2,2> v;v.value=.5*r*r/(o.scale*o.scale);v.x[o.component]=r/(o.scale*o.scale);return v;};cost.directional_gradient=[o](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>&){ScalarDirection<2,2> d;d.x[o.component]=s[o.component]/(o.scale*o.scale);return d;};return cost;}
inline System<2,2> system(std::size_t steps,double theta,bool running=false,bool added=false,bool noisy=true){
    auto grid=Grid::uniform(t0,tf,steps,theta);Problem<2,2> p;
    p.initial.value=[](const Vec<2>& xi){return xi;};p.initial.jacobian=[](const Vec<2>&){return identity<2>();};p.initial.directional_jacobian=[](const Vec<2>&,const Vec<2>&){return Mat<2,2>{};};
    p.dynamics.value=[](double,const Vec<2>& x,const Vec<2>&){return mv(A,x);};p.dynamics.derivatives=[](double,const Vec<2>&,const Vec<2>&){Derivatives<2,2> d;d.x=A;return d;};p.dynamics.directional_derivatives=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>&,const Vec<2>&){return Derivatives<2,2>{};};
    p.dynamics.admissible=[](double,const Vec<2>& x,const Vec<2>&){return x[0]>=-1e-12&&x[1]>=-1e-12;};
    if(running){p.running.evaluate=[](double,const Vec<2>& x,const Vec<2>&){return ScalarValue<2,2>{.5*x[0]*x[0],{x[0],0.},{0.,0.}};};p.running.directional_gradient=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>&){return ScalarDirection<2,2>{{s[0],0.},{0.,0.}};};}
    else for(auto o:observations(added,noisy))p.nodal.emplace_back(node_at(grid,o.tau),measurement(o));
    return System<2,2>(p,grid);
}
inline Optimum<2> fit(const System<2,2>& sys){return minimize(sys,Vec<2>{.2,.8},OptimizeOptions{1e-8,6,40});}
inline double exact_cost(const Vec<2>& xi){double J=0.;for(auto o:observations()){double r=exact(o.tau,xi)[o.component]-o.value;J+=.5*r*r/(o.scale*o.scale);}return J;}
inline Vec<2> exact_gradient(const Vec<2>& xi){Vec<2> g{};for(auto o:observations())g=g+((exact(o.tau,xi)[o.component]-o.value)/(o.scale*o.scale))*propagator(o.tau)[o.component];return g;}
inline Vec<2> nodal_costate(double tau,const Vec<2>& xi){Vec<2> lambda{};for(auto o:observations())if(o.tau+1e-12>=tau)lambda=lambda+((exact(o.tau,xi)[o.component]-o.value)/(o.scale*o.scale))*propagator(o.tau-tau)[o.component];return lambda;}
inline Vec<2> running_costate(double tau,const Vec<2>& xi){return mv(gramian(horizon-tau),exact(tau,xi));}
inline double analytical_gain(double tau){auto P=propagator(tau);double tr=P[0][0]*P[0][0]+P[0][1]*P[0][1]+P[1][1]*P[1][1],det=std::exp(-2*(fast+slow)*tau);return std::sqrt(.5*(tr+std::sqrt(std::max(0.,tr*tr-4*det))));}
inline double gain_derivative(double tau){double ea=std::exp(-fast*tau),ed=std::exp(-slow*tau),off=mixing*(ed-ea),offdot=mixing*(-slow*ed+fast*ea);double tr=ea*ea+ed*ed+off*off,dt=-2*fast*ea*ea-2*slow*ed*ed+2*off*offdot,det=std::exp(-2*(fast+slow)*tau),dd=-2*(fast+slow)*det;return .5*(dt+(tr*dt-2*dd)/std::sqrt(tr*tr-4*det));}
inline double peak_time(){double lo=1e-5,hi=horizon;for(int j=0;j<100;++j){double mid=.5*(lo+hi);if(gain_derivative(mid)>0)lo=mid;else hi=mid;}return .5*(lo+hi);}
}
