// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include "insulated_thermal_model.hpp"
namespace insulated_thermal::qualification {
inline constexpr std::array<Vec<2>,4> goals{{{1.,0.},{0.,1.},difference,mean}};
inline constexpr std::array<double,4> tolerances{.3,.3,.05,.2};
inline constexpr std::array<const char*,4> names{"temperature1","temperature2","difference","weighted_mean"};
inline constexpr std::array<const char*,4> units{"degC","degC","K_difference","degC"};
inline ConvexPolygon2 prior(){return rectangle_polygon({10.,10.},{40.,40.});}
inline double dual(Vec<2> c){return std::abs(c[0])+std::abs(c[1]);}
inline std::vector<double> allowance(const System<2,2>& sys){auto& g=sys.grid();std::vector<double> D(g.steps()),E(g.time.size());for(auto v:prior()){auto tr=sys.forward(v);for(std::size_t n=0;n<g.steps();++n)D[n]=std::max(D[n],linear_theta_quadratic_defect_bound(A,tr.state[n],tr.state[n+1],g.time[n+1]-g.time[n],g.theta[n],[](auto x){return norm(x);}));}for(std::size_t n=0;n<g.steps();++n)E[n+1]=propagate_defect_bound(E[n],D[n],g.time[n+1]-g.time[n],0.);return E;}
inline ConvexPolygon2 feasible(const System<2,2>& sys,const std::vector<Mat<2,2>>& J,const std::vector<double>& E,bool added,bool analytical){auto p=prior();for(auto o:observations(added)){auto n=node_at(sys.grid(),o.tau);auto row=analytical?observation_row(o):mv(transpose(J[n]),o.goal);double width=o.scale+(analytical?0.:dual(o.goal)*E[n]);p=clip_halfspace(std::move(p),row,o.value+width);p=clip_halfspace(std::move(p),-1.*row,-o.value+width);}if(p.size()<3)throw std::runtime_error("Empty thermal feasible polygon");return p;}
inline double exact_bound(double tau,const ConvexPolygon2& p,Vec<2> center,Vec<2> c){return polygon_goal_bound(p,mv(transpose(propagator(tau)),c),center);}
// For each vertex the goal error is a constant plus one exponential. Its
// absolute value is bounded on an interval by its endpoints.
inline std::optional<double> exact_persistent(const ConvexPolygon2& p,Vec<2> center,Vec<2> c,double tolerance){if(exact_bound(horizon,p,center,c)>tolerance)return {};if(exact_bound(0.,p,center,c)<=tolerance)return 0.;double lo=0.,hi=horizon;for(int j=0;j<90;++j){double mid=.5*(lo+hi);if(exact_bound(mid,p,center,c)>tolerance)lo=mid;else hi=mid;}return hi;}
struct Result {Vec<2> center{};ConvexPolygon2 polygon,reference;std::vector<double> error;std::vector<std::array<double,4>> whole;IntervalQualification<4> decision;std::array<std::optional<double>,4> analytical_time;};
template<class Check> Result experiment(const System<2,2>& sys,bool added,Vec<2> start,Check&& check){Result r;r.center=fit(sys,start).optimum.control;auto tr=sys.forward(r.center);auto J=control_to_state_jacobians(sys,r.center,tr);auto& grid=sys.grid();r.error=allowance(sys);r.polygon=feasible(sys,J,r.error,added,false);r.reference=feasible(sys,J,r.error,added,true);
 check("exact feasible set contains truth",polygon_contains(r.reference,truth)?0.:1.,0.);check("numerical feasible set contains truth",polygon_contains(r.polygon,truth)?0.:1.,0.);for(auto v:r.reference)check("outer polygon contains analytical vertex",polygon_contains(r.polygon,v)?0.:1.,0.);
 double goal_pairing=0.;for(std::size_t n:std::array<std::size_t,3>{0,grid.steps()/2,grid.steps()})for(auto goal:goals)goal_pairing=std::max(goal_pairing,norm(sys.goal_gradient(r.center,tr,n,goal)-mv(transpose(J[n]),goal)));check("goal adjoints agree with tangent maps",goal_pairing,2e-11);
 r.whole.resize(grid.steps());double enclosure=0.,error_violation=0.;for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n],tau=grid.time[n]-t0;for(std::size_t k=0;k<4;++k){double spread=0.;for(auto v:r.polygon){auto x=mv(J[n],v-r.center),next=mv(J[n+1],v-r.center);spread=std::max(spread,quadratic_absolute_maximum(dot(goals[k],x),h*dot(goals[k],mv(A,x)),dot(goals[k],next-x-h*mv(A,x))));}r.whole[n][k]=spread+dual(goals[k])*std::max(r.error[n],r.error[n+1]);}
 for(double s:{0.,.25,.5,.75,1.}){auto center=linear_quadratic_reconstruction(A,tr.state[n],tr.state[n+1],h,s);for(auto v:r.reference){auto delta=exact(tau+h*s,v)-center;for(std::size_t k=0;k<4;++k)enclosure=std::max(enclosure,(std::abs(dot(goals[k],delta))-r.whole[n][k])/tolerances[k]);}}
 }
 for(auto corner:prior()){auto path=sys.forward(corner);for(std::size_t n=0;n<grid.steps();++n)for(double s:{0.,.5,1.}){double h=grid.time[n+1]-grid.time[n];error_violation=std::max(error_violation,norm(exact(grid.time[n]-t0+h*s,corner)-linear_quadratic_reconstruction(A,path.state[n],path.state[n+1],h,s))-std::max(r.error[n],r.error[n+1]));}}
 check("quadratic defect allowance versus exact corners",error_violation,2e-11);check("whole interval envelopes enclose analytical feasible probes",enclosure,2e-10);
 std::vector<double> times;for(double t:grid.time)times.push_back(t-t0);r.decision=qualify_intervals(times,r.whole,tolerances);auto reference_center=analytical_fit(added,start);for(std::size_t k=0;k<4;++k){r.analytical_time[k]=exact_persistent(r.reference,reference_center,goals[k],tolerances[k]);if(r.decision.persistent[k]){double maxratio=0.;for(std::size_t n=0;n<grid.steps();++n)if(times[n]>=*r.decision.persistent[k]-1e-12)maxratio=std::max(maxratio,r.whole[n][k]/tolerances[k]);check("qualified remaining whole intervals",maxratio,1.);}}
 return r;}
}
