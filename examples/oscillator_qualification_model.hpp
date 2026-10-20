#pragma once
#include "oscillator_model.hpp"
#include <optional>

namespace oscillator::qualification {
using Polygon=ConvexPolygon2;
inline constexpr std::array<Vec<2>,3> goals{{{1.,0.},{0.,1.},{stiffness,damping}}};
inline constexpr std::array<double,3> tolerances{.001,.002,.04};
inline constexpr std::array<const char*,3> names{"displacement","velocity","support_force"};
inline constexpr std::array<const char*,3> units{"m","m/s","N"};
inline constexpr Vec<2> prior_radii{.08,.30};
inline Polygon prior(){return rectangle_polygon(-1.*prior_radii,prior_radii);}
inline Vec<2> transform(const Vec<2>& x){return {frequency*x[0],x[1]+decay*x[0]};}
inline double energy_norm(const Vec<2>& x){auto z=transform(x);return std::hypot(z[0],z[1]);}
inline double goal_dual_norm(const Vec<2>& c){return std::hypot((c[0]-decay*c[1])/frequency,c[1]);}
inline Polygon clip(Polygon p,const Vec<2>& normal,double bound){return clip_halfspace(std::move(p),normal,bound);}
inline bool contains(const Polygon& p,const Vec<2>& x,double tolerance=1e-13){return polygon_contains(p,x,tolerance);}
struct Affine {std::array<std::vector<Vec<2>>,2> tangent;};
inline Affine affine(const System<2,2>& sys){auto tr=sys.forward({0.,0.});return {{sys.tangent({0.,0.},tr,{1.,0.}),sys.tangent({0.,0.},tr,{0.,1.})}};}
inline Vec<2> apply(const Affine& path,std::size_t node,const Vec<2>& xi){return xi[0]*path.tangent[0][node]+xi[1]*path.tangent[1][node];}
inline Vec<2> goal_gradient(const Affine& path,std::size_t node,const Vec<2>& c){return {dot(c,path.tangent[0][node]),dot(c,path.tangent[1][node])};}
// Quadratic reconstruction matches both state endpoints. At theta=.5 it also
// matches the model derivative at both ends when the algebraic residual is zero.
inline Vec<2> reconstruct(const Vec<2>& x,const Vec<2>& next,double h,double s){return linear_quadratic_reconstruction(A,x,next,h,s);}
inline double quadratic_peak(double c0,double c1,double c2){return quadratic_absolute_maximum(c0,c1,c2);}
inline std::vector<double> allowance(const System<2,2>& sys){
    const auto& grid=sys.grid();std::vector<double> D(grid.steps()),E(grid.time.size());
    for(auto corner:prior()){auto tr=sys.forward(corner);
        for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n],theta=grid.theta[n];auto x=tr.state[n],next=tr.state[n+1];
            double defect=linear_theta_quadratic_defect_bound(A,x,next,h,theta,energy_norm);
            D[n]=std::max(D[n],defect);
        }
    }
    // In z=(omega*q,v+a*q), z'=(-a I+omega J)z and ||Phi||=exp(-a*t).
    for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n];E[n+1]=propagate_defect_bound(E[n],D[n],h,-decay);}return E;
}
inline Polygon observation_set(const Grid& grid,const Affine& path,double radius,const std::vector<double>& E,bool analytical=false){
    auto p=prior();for(std::size_t j=0;j<elapsed.size();++j){auto n=node_at(grid,elapsed[j]);auto g=analytical?propagator(elapsed[j])[0]:goal_gradient(path,n,{1.,0.});double b=radius+(analytical?0.:goal_dual_norm({1.,0.})*E[n]);p=clip(std::move(p),g,observation(j)+b);p=clip(std::move(p),-1.*g,-observation(j)+b);}if(p.size()<3)throw std::runtime_error("Empty oscillator uncertainty polygon");return p;
}
inline double analytical_spread(double tau,const Polygon& set,const Vec<2>& center,const Vec<2>& goal){auto g=mv(transpose(propagator(tau)),goal);double maximum=0.;for(auto v:set)maximum=std::max(maximum,std::abs(dot(g,v-center)));return maximum;}
inline std::optional<double> analytical_persistent_time(const Polygon& set,const Vec<2>& center,const Vec<2>& goal,double tolerance){
    double latest=0.;for(auto vertex:set){auto delta=vertex-center;double b=dot(goal,delta),c=dot(goal,(1/frequency)*(mv(A,delta)+decay*delta));
        auto error=[&](double tau){return std::exp(-decay*tau)*std::abs(b*std::cos(frequency*tau)+c*std::sin(frequency*tau))-tolerance;};
        if(error(horizon)>0)return std::nullopt;
        std::vector<double> critical{0.,horizon};double U=-decay*b+frequency*c,V=-decay*c-frequency*b;
        if(U!=0||V!=0){double phase=std::atan2(-U,V);for(int k=-2;k<20;++k){double tau=(phase+k*std::numbers::pi)/frequency;if(tau>0&&tau<horizon)critical.push_back(tau);}}
        // Include zero crossings, so each absolute-value segment is monotone.
        if(b!=0||c!=0){double phase=std::atan2(-b,c);for(int k=-2;k<20;++k){double tau=(phase+k*std::numbers::pi)/frequency;if(tau>0&&tau<horizon)critical.push_back(tau);}}
        std::sort(critical.begin(),critical.end());
        for(std::size_t j=0;j+1<critical.size();++j)if(error(critical[j])>0&&error(critical[j+1])<=0){double lo=critical[j],hi=critical[j+1];for(int k=0;k<70;++k){double mid=.5*(lo+hi);if(error(mid)>0)lo=mid;else hi=mid;}latest=std::max(latest,hi);}
    }return latest;
}
struct Result {
    Vec<2> center{};Polygon polygon,reference_polygon;std::vector<double> numerical;
    std::vector<std::array<double,3>> nodes,intervals,local95;
    std::array<std::optional<double>,3> node_time,continuous_time,analytical_time,first_interval_time;
    std::optional<double> joint;double radius{};
};
template<class Check> Result experiment(const System<2,2>& sys,double radius,Check&& check){
    Result result;result.radius=radius;auto opt=fit(sys);result.center=opt.control;auto tr=sys.forward(result.center);auto path=affine(sys);const auto& grid=sys.grid();result.numerical=allowance(sys);
    result.polygon=observation_set(grid,path,radius,result.numerical);result.reference_polygon=observation_set(grid,path,radius,result.numerical,true);
    check("analytical set contains synthetic truth",contains(result.reference_polygon,truth)?0.:1.,0.);
    check("numerical set contains synthetic truth",contains(result.polygon,truth)?0.:1.,0.);
    for(auto v:result.reference_polygon)check("corrected numerical set contains analytical vertex",contains(result.polygon,v)?0.:1.,0.);
    double bound_check=0.;for(auto corner:prior()){auto trajectory=sys.forward(corner);for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n];for(double s:{0.,.25,.5,.75,1.}){auto y=reconstruct(trajectory.state[n],trajectory.state[n+1],h,s);double E=std::max(result.numerical[n],result.numerical[n+1]);bound_check=std::max(bound_check,energy_norm(y-exact(grid.time[n]-t0+s*h,corner))-E);}}}check("uniform reconstruction allowance versus exact solution",bound_check,2e-13);
    Mat<2,2> I{};for(std::size_t j=0;j<elapsed.size();++j){auto node=node_at(grid,elapsed[j]);auto g=sys.goal_gradient(result.center,tr,node,{1.,0.});I=add(I,scaled(1/(sigma*sigma),outer(g)));}auto covariance=inverse(I);
    result.nodes.resize(grid.time.size());result.local95.resize(grid.time.size());result.intervals.resize(grid.steps());
    double adjoint_error=0.,prior_support_error=0.,hessian_error=0.;
    for(int j=0;j<2;++j){Vec<2> unit{};unit[j]=1.;hessian_error=std::max(hessian_error,norm(sys.hessian_vector(result.center,unit)-mv(I,unit)));}check("goal pullback information equals incremental adjoint",hessian_error/std::max(norm(I[0]),norm(I[1])),2e-13);
    for(std::size_t n=0;n<grid.time.size();++n)for(std::size_t k=0;k<3;++k){auto c=goals[k],g=goal_gradient(path,n,c);
        if(n==0||n==grid.steps()/2||n==grid.steps()){auto adj=sys.goal_gradient(result.center,tr,n,c);adjoint_error=std::max(adjoint_error,norm(adj-g));}
        double spread=0.;for(auto v:result.polygon){spread=std::max(spread,std::abs(dot(g,v-result.center)));prior_support_error=std::max(prior_support_error,std::abs(dot(g,v))-box_support(g,prior_radii));}
        result.nodes[n][k]=spread+goal_dual_norm(c)*result.numerical[n];
        result.local95[n][k]=std::sqrt(-2*std::log(.05))*std::sqrt(local_goal_variance(g,covariance));
    }
    check("goal adjoints equal analytical discrete tangents",adjoint_error,2e-11);check("prior box support encloses polygon projections",prior_support_error,2e-12);
    for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n];for(std::size_t k=0;k<3;++k){double spread=0.;for(auto v:result.polygon){auto x=apply(path,n,v-result.center),next=apply(path,n+1,v-result.center);double c0=dot(goals[k],x),c1=h*dot(goals[k],mv(A,x)),c2=dot(goals[k],next-x-h*mv(A,x));spread=std::max(spread,quadratic_peak(c0,c1,c2));}result.intervals[n][k]=spread+goal_dual_norm(goals[k])*std::max(result.numerical[n],result.numerical[n+1]);}}
    std::vector<double> times;for(double time:grid.time)times.push_back(time-t0);
    auto qualified=qualify_intervals(times,result.intervals,tolerances);result.first_interval_time=qualified.first_below;result.continuous_time=qualified.persistent;result.joint=qualified.joint;
    for(std::size_t k=0;k<3;++k){std::vector<double> node_ratio;for(auto v:result.nodes)node_ratio.push_back(v[k]/tolerances[k]);result.node_time[k]=persistent_node_qualification(times,node_ratio);result.analytical_time[k]=analytical_persistent_time(result.reference_polygon,analytical_fit(),goals[k],tolerances[k]);
        if(result.continuous_time[k]){double last=0.;for(std::size_t n=0;n<grid.steps();++n)if(times[n]>=*result.continuous_time[k]-1e-12)last=std::max(last,result.intervals[n][k]/tolerances[k]);check("qualified suffix covers every whole interval",last,1.);}
    }
    double physical_violation=0.;
    for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n];for(double s:{0.,.25,.5,.75,1.}){auto center=reconstruct(tr.state[n],tr.state[n+1],h,s);double tau=grid.time[n]-t0+h*s;for(auto vertex:result.reference_polygon){auto error=exact(tau,vertex)-center;for(std::size_t k=0;k<3;++k)physical_violation=std::max(physical_violation,(std::abs(dot(goals[k],error))-result.intervals[n][k])/tolerances[k]);}}}
    check("whole-interval envelopes enclose analytical physical-set probes",physical_violation,2e-12);
    return result;
}
}
