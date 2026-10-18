#pragma once
#include <theta_opt/optimize.hpp>
#include <array>
#include <cmath>

namespace decay {
using namespace theta_opt;
inline constexpr double a=.5,b=.1,kappa=.1,gamma=1.,beta=.4,prior=3.2,t0=2.,tf=12.;
struct Exact {Vec<2> state{},tangent{},second{};double cost{},cost_tangent{},cost_second{};};
inline Exact exact(double elapsed,double xi){
    double r=std::exp(-a*elapsed),k=b/a*(-std::expm1(-a*elapsed)),d=1+k*xi;
    Exact e;e.state={xi*r/d,std::log1p(k*xi)/b};e.tangent={r/(d*d),k/(b*d)};e.second={-2*k*r/(d*d*d),-k*k/(b*d*d)};
    e.cost=kappa*(xi-e.state[0]-a*e.state[1])/(2*b);
    e.cost_tangent=kappa*(1-e.tangent[0]-a*e.tangent[1])/(2*b);
    e.cost_second=kappa*(-e.second[0]-a*e.second[1])/(2*b);return e;
}
inline double target(){return exact(tf-t0,3.).state[1];}
inline std::array<double,3> reference(double xi){auto e=exact(tf-t0,xi);double r=e.state[1]-target();return {e.cost+.5*gamma*r*r+.5*beta*(xi-prior)*(xi-prior),e.cost_tangent+gamma*r*e.tangent[1]+beta*(xi-prior),e.cost_second+gamma*(e.tangent[1]*e.tangent[1]+r*e.second[1])+beta};}
inline Vec<2> costate(double t,double xi){auto now=exact(t-t0,xi),future=exact(tf-t,now.state[0]);double lq=gamma*(exact(tf-t0,xi).state[1]-target());return {future.cost_tangent+lq*future.tangent[1],lq};}
inline Problem<2,1> problem(){
    Problem<2,1> p;
    p.initial.value=[](const Vec<1>& xi){return Vec<2>{xi[0],0.};};
    p.initial.jacobian=[](const Vec<1>&){return Mat<2,1>{{{1.},{0.}}};};
    p.initial.directional_jacobian=[](const Vec<1>&,const Vec<1>&){return Mat<2,1>{};};
    p.dynamics.value=[](double,const Vec<2>& x,const Vec<1>&){return Vec<2>{-a*x[0]-b*x[0]*x[0],x[0]};};
    p.dynamics.derivatives=[](double,const Vec<2>& x,const Vec<1>&){Derivatives<2,1> d;d.x={{{-a-2*b*x[0],0.},{1.,0.}}};return d;};
    p.dynamics.directional_derivatives=[](double,const Vec<2>&,const Vec<1>&,const Vec<2>& v,const Vec<1>&){Derivatives<2,1> d;d.x[0][0]=-2*b*v[0];return d;};
    p.dynamics.admissible=[](double,const Vec<2>& x,const Vec<1>&){return x[0]>0;};
    p.running.evaluate=[](double,const Vec<2>& x,const Vec<1>&){return ScalarValue<2,1>{.5*kappa*x[0]*x[0],{kappa*x[0],0.},{0.}};};
    p.running.directional_gradient=[](double,const Vec<2>&,const Vec<1>&,const Vec<2>& v,const Vec<1>&){return ScalarDirection<2,1>{{kappa*v[0],0.},{0.}};};
    p.terminal.evaluate=[](double,const Vec<2>& x,const Vec<1>&){double r=x[1]-target();return ScalarValue<2,1>{.5*gamma*r*r,{0.,gamma*r},{0.}};};
    p.terminal.directional_gradient=[](double,const Vec<2>&,const Vec<1>&,const Vec<2>& v,const Vec<1>&){return ScalarDirection<2,1>{{0.,gamma*v[1]},{0.}};};
    p.penalty.value=[](const Vec<1>& xi){return .5*beta*(xi[0]-prior)*(xi[0]-prior);};
    p.penalty.gradient=[](const Vec<1>& xi){return Vec<1>{beta*(xi[0]-prior)};};
    p.penalty.hessian_vector=[](const Vec<1>&,const Vec<1>& d){return beta*d;};return p;
}
template<class Function> double root(Function g,double lo=.1,double hi=6.){
    double gl=g(lo),gh=g(hi);if(gl==0)return lo;if(gh==0)return hi;if(gl*gh>0)throw std::runtime_error("Root not bracketed");
    for(int k=0;k<100;++k){double mid=.5*(lo+hi);if(mid==lo||mid==hi)return mid;double gm=g(mid);if(gm==0)return mid;if(gl*gm<=0)hi=mid;else{lo=mid;gl=gm;}}return .5*(lo+hi);
}
}
