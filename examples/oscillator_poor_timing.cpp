// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#include "oscillator_model.hpp"
#include "oscillator_qualification_model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

// A second oscillator experiment: all displacement observations are taken at
// exact half periods. The fit uses the actual theta-discrete objective; the
// engineering feasible set is formed from the exact physical observation map
// and enlarged by a certified theta reconstruction allowance.
namespace {
using namespace oscillator;
namespace oq=oscillator::qualification;
using Schedule=std::array<double,4>;
inline constexpr Schedule poor_times{
    std::numbers::pi/2.,std::numbers::pi,3*std::numbers::pi/2.,2*std::numbers::pi};
inline constexpr Schedule proper_times{.3,.9,1.8,3.2};
inline constexpr std::array<double,4> noise{.4,-.6,.2,.8};
inline constexpr Vec<2> control_scale{.02,.12};
inline constexpr Vec<2> prior_u{4.,2.5};
inline constexpr double radius=sigma;

double obs(std::size_t j,const Schedule& schedule){return exact(schedule[j],truth)[0]+sigma*noise[j];}
std::size_t node(const Grid& g,double tau){for(std::size_t n=0;n<g.time.size();++n)if(std::abs(g.time[n]-t0-tau)<2e-12)return n;throw std::logic_error("half-period time missing from grid");}
Grid mesh(std::size_t n,double theta,const Schedule& schedule){
    std::vector<double> t; t.reserve(n+schedule.size()+2);
    for(std::size_t j=0;j<=n;++j)t.push_back(t0+horizon*static_cast<double>(j)/static_cast<double>(n));
    for(double tau:schedule)t.push_back(t0+tau);
    std::sort(t.begin(),t.end());t.erase(std::unique(t.begin(),t.end(),[](double a,double b){return std::abs(a-b)<2e-14;}),t.end());
    const auto intervals=t.size()-1;return Grid(std::move(t),std::vector<double>(intervals,theta));
}
System<2,2> make_system(std::size_t n,double theta,const Schedule& schedule){
    auto g=mesh(n,theta,schedule);Problem<2,2> p;
    p.initial.value=[](const Vec<2>& u){return Vec<2>{control_scale[0]*u[0],control_scale[1]*u[1]};};
    p.initial.jacobian=[](const Vec<2>&){return Mat<2,2>{{{control_scale[0],0.},{0.,control_scale[1]}}};};
    p.initial.directional_jacobian=[](const Vec<2>&,const Vec<2>&){return Mat<2,2>{};};
    p.dynamics.value=[](double,const Vec<2>& x,const Vec<2>&){return mv(A,x);};
    p.dynamics.derivatives=[](double,const Vec<2>&,const Vec<2>&){Derivatives<2,2> d;d.x=A;return d;};
    p.dynamics.directional_derivatives=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>&,const Vec<2>&){return Derivatives<2,2>{};};
    for(std::size_t j=0;j<schedule.size();++j){auto c=ScalarFunction<2,2>::zero();const double y=obs(j,schedule);
        c.evaluate=[y](double,const Vec<2>& x,const Vec<2>&){ScalarValue<2,2> v;const double r=x[0]-y;v.value=.5*r*r/(sigma*sigma);v.x[0]=r/(sigma*sigma);return v;};
        c.directional_gradient=[](double,const Vec<2>&,const Vec<2>&,const Vec<2>& s,const Vec<2>&){ScalarDirection<2,2> v;v.x[0]=s[0]/(sigma*sigma);return v;};
        p.nodal.emplace_back(node(g,schedule[j]),c);
    }
    return System<2,2>(std::move(p),std::move(g));
}
Vec<2> physical(const Vec<2>& u){return {control_scale[0]*u[0],control_scale[1]*u[1]};}
Mat<2,2> information(const System<2,2>& s,const Vec<2>& u){auto e=s.evaluate(u);Mat<2,2> H{};for(int j=0;j<2;++j){Vec<2> d{};d[j]=1.;auto col=s.hessian_vector(u,d);for(int i=0;i<2;++i)H[i][j]=col[i];}return H;}
std::vector<double> allowance(const System<2,2>& sys){
    const auto& g=sys.grid();std::vector<double> D(g.steps()),E(g.time.size());
    for(double uq:{-prior_u[0],prior_u[0]})for(double uv:{-prior_u[1],prior_u[1]}){auto tr=sys.forward({uq,uv});for(std::size_t n=0;n<g.steps();++n){double h=g.time[n+1]-g.time[n];D[n]=std::max(D[n],linear_theta_quadratic_defect_bound(A,tr.state[n],tr.state[n+1],h,g.theta[n],oq::energy_norm));}}
    for(std::size_t n=0;n<g.steps();++n)E[n+1]=propagate_defect_bound(E[n],D[n],g.time[n+1]-g.time[n],-decay);
    return E;
}
struct Qualification {ConvexPolygon2 exact_set,numerical_set;Vec<2> center{};std::array<std::optional<double>,3> times;std::optional<double> joint;};
Qualification qualify(const System<2,2>& sys,double error_radius,const Schedule& schedule){
    Qualification q;const auto& g=sys.grid();auto path=oq::affine(sys);auto E=allowance(sys);q.exact_set=q.numerical_set=rectangle_polygon(-1.*prior_u,prior_u);
    for(std::size_t j=0;j<schedule.size();++j){auto row=propagator(schedule[j])[0];Vec<2> a{row[0]*control_scale[0],row[1]*control_scale[1]};q.exact_set=clip_halfspace(std::move(q.exact_set),a,obs(j,schedule)+error_radius);q.exact_set=clip_halfspace(std::move(q.exact_set),-1.*a,-obs(j,schedule)+error_radius);
        auto n=node(g,schedule[j]);auto b=oq::goal_gradient(path,n,{1.,0.});double r=error_radius+oq::goal_dual_norm({1.,0.})*E[n];q.numerical_set=clip_halfspace(std::move(q.numerical_set),b,obs(j,schedule)+r);q.numerical_set=clip_halfspace(std::move(q.numerical_set),-1.*b,-obs(j,schedule)+r);}
    if(q.exact_set.size()<3||q.numerical_set.size()<3)throw std::runtime_error("Empty half-period feasible set");
    // Pick an explicit, feasible representative for goal errors. The discrete
    // least-squares fit can sit outside the exact bounded-data set because its
    // phase error acts like false information; do not center qualification on it.
    for(auto u:q.exact_set)q.center=q.center+(1./static_cast<double>(q.exact_set.size()))*u;
    auto tr=sys.forward(q.center);
    std::vector<std::array<double,3>> intervals(g.steps());std::vector<double> edges;edges.reserve(g.time.size());for(double t:g.time)edges.push_back(t-t0);
    for(std::size_t n=0;n<g.steps();++n){double h=g.time[n+1]-g.time[n];for(std::size_t k=0;k<3;++k){double bound=0.;for(auto u:q.numerical_set){auto d=u-q.center;auto x=oq::apply(path,n,d),next=oq::apply(path,n+1,d);double c0=dot(oq::goals[k],x),c1=h*dot(oq::goals[k],mv(A,x)),c2=dot(oq::goals[k],next-x-h*mv(A,x));bound=std::max(bound,quadratic_absolute_maximum(c0,c1,c2));}intervals[n][k]=bound+oq::goal_dual_norm(oq::goals[k])*std::max(E[n],E[n+1]);}}
    auto result=qualify_intervals(edges,intervals,oq::tolerances);q.times=result.persistent;q.joint=result.joint;return q;
}
}

int main(int argc,char** argv){try{
    const bool proper=argc>2&&std::string(argv[2])=="proper";
    const Schedule& schedule=proper?proper_times:poor_times;
    const std::string case_name=proper?"proper":"poor";
    const std::filesystem::path out=argc>1?argv[1]:(proper?"results/oscillator_proper_timing":"results/oscillator_poor_timing");std::filesystem::create_directories(out);
    std::ofstream summary(out/"summary.csv"),observations(out/"observations.csv"),states(out/"reconstruction.csv"),qualification(out/"qualification.csv"),sets(out/"feasible_sets.csv"),orders(out/"convergence.csv");
    summary<<std::setprecision(17)<<"quantity,value\n"<<"sampling_design,"<<case_name<<'\n';observations<<std::setprecision(17)<<"elapsed_s,observed_m,exact_truth_m,measurement_bound_m\n";for(std::size_t j=0;j<schedule.size();++j)observations<<schedule[j]<<','<<obs(j,schedule)<<','<<exact(schedule[j],truth)[0]<<','<<radius<<'\n';
    qualification<<"goal,units,tolerance,persistent_time_s,qualified,joint_time_s,joint_qualified\n";sets<<std::setprecision(17)<<"set,vertex,q0_m,v0_m_s\n";orders<<std::setprecision(17)<<"steps,h_s,q0_fit_m,v0_fit_m_s,scaled_control_error,fixed_state_error,scaled_information_small_eigenvalue\n";
    std::size_t checks=0;
    auto add_check=[&](const std::string& name,double error,double tol){++checks;if(!std::isfinite(error)||error>tol)throw std::runtime_error(name+" failed: "+std::to_string(error));};
    double exact_sensitivity=0.;Mat<2,2> analytical_H{};Vec<2> analytical_rhs{};
    for(std::size_t j=0;j<schedule.size();++j){double tau=schedule[j];auto phi=propagator(tau);exact_sensitivity=std::max(exact_sensitivity,std::abs(phi[0][1]));Vec<2> r{phi[0][0]*control_scale[0],phi[0][1]*control_scale[1]};analytical_H=add(analytical_H,scaled(1/(sigma*sigma),outer(r)));analytical_rhs=analytical_rhs+((obs(j,schedule)/(sigma*sigma))*r);}
    auto analytic_rank=psd_spectrum_2d(analytical_H,1e-11);
    add_check("analytical observation rank",std::abs(static_cast<double>(analytic_rank.rank)-(proper?2.:1.)),0.);
    Vec<2> exact_u_fit{};if(proper)exact_u_fit=LU<2>(analytical_H).solve(analytical_rhs);else {double qnumerator=0.,qdenominator=0.;for(std::size_t j=0;j<schedule.size();++j){double a=propagator(schedule[j])[0][0];qnumerator+=a*obs(j,schedule);qdenominator+=a*a;}exact_u_fit[0]=(qnumerator/qdenominator)/control_scale[0];}
    const auto exact_fit=physical(exact_u_fit);
    std::cout<<(proper?"Properly":"Poorly")<<" timed oscillator: 4 displacement samples; analytical scaled observation rank "<<analytic_rank.rank<<" / 2.\n";
    std::cout<<"Maximum exact initial-velocity sensitivity: "<<std::setprecision(10)<<exact_sensitivity<<" s.\n";
    std::cout<<"Analytical least-squares representative: (q0,v0)=("<<exact_fit[0]<<", "<<exact_fit[1]<<") [m, m/s].\n";
    double previous_state_error=0.,previous_control_error=0.;std::vector<std::array<double,5>> fit_rows;
    for(std::size_t n:{64u,128u,256u,512u}){auto sys=make_system(n,.5,schedule);auto spectrum=psd_spectrum_2d(information(sys,{0.,0.}),1e-14);auto opt=minimize(sys,Vec<2>{0.,0.},OptimizeOptions{1e-8,8,40});auto xi=physical(opt.control);double h=horizon/n;auto fixed=sys.forward({truth[0]/control_scale[0],truth[1]/control_scale[1]});double fixed_error=0.;for(std::size_t k=0;k<fixed.state.size();++k){double tau=sys.grid().time[k]-t0;auto e=fixed.state[k]-exact(tau,truth);fixed_error=std::max(fixed_error,std::max(std::abs(e[0])/control_scale[0],std::abs(e[1])/control_scale[1]));}
        // The full discrete objective may weakly lift the exact physical null
        // direction through phase error. Report that tiny eigenvalue, but use
        // the continuous exact feasible set—not that artifact—for qualification.
        double control_error=norm(Vec<2>{(xi[0]-exact_fit[0])/control_scale[0],(xi[1]-exact_fit[1])/control_scale[1]});
        orders<<n<<','<<h<<','<<xi[0]<<','<<xi[1]<<','<<control_error<<','<<fixed_error<<','<<spectrum.eigenvalues[1]<<'\n';
        fit_rows.push_back({static_cast<double>(n),xi[0],xi[1],spectrum.eigenvalues[1],fixed_error});
        if(previous_state_error>0){double order=std::log(previous_state_error/fixed_error)/std::log(2.);add_check("midpoint forward state convergence order",std::abs(order-2.),.08);summary<<"fixed_state_order_N_"<<(n/2)<<"_to_"<<n<<','<<order<<'\n';}
        if(proper&&previous_control_error>0){double order=std::log(previous_control_error/control_error)/std::log(2.);summary<<"proper_fit_order_N_"<<(n/2)<<"_to_"<<n<<','<<order<<'\n';if(n>=128)add_check("properly sampled recovered-control convergence",std::abs(order-2.),.25);}
        previous_state_error=fixed_error;previous_control_error=control_error;
        add_check("discrete optimizer stationarity",norm(sys.evaluate(opt.control).gradient),1e-8);
        if(n==128){auto ev=sys.evaluate({.1,-.2});auto d=Vec<2>{.03,-.04};auto tangent=sys.tangent({.1,-.2},ev.trajectory,d);double pair=0.;for(std::size_t j=0;j<schedule.size();++j){auto k=node(sys.grid(),schedule[j]);pair+=(ev.trajectory.state[k][0]-obs(j,schedule))*tangent[k][0]/(sigma*sigma);}add_check("theta tangent-adjoint pairing",std::abs(pair-dot(ev.gradient,d)),2e-8);
            auto eps=1e-5;auto ep=sys.evaluate(Vec<2>{.1+eps*d[0],-.2+eps*d[1]}),em=sys.evaluate(Vec<2>{.1-eps*d[0],-.2-eps*d[1]});add_check("discrete objective adjoint gradient",std::abs((ep.cost-em.cost)/(2*eps)-dot(ev.gradient,d)),2e-5);
            auto hp=sys.hessian_vector({.1,-.2},d);add_check("incremental-adjoint Hessian",norm(hp-(1/(2*eps))*(ep.gradient-em.gradient)),2e-5);
        }
        if(n==512){auto result=qualify(sys,radius,schedule);auto truth_u=Vec<2>{truth[0]/control_scale[0],truth[1]/control_scale[1]};add_check("analytical feasible set contains truth",polygon_contains(result.exact_set,truth_u)?0.:1.,0.);double containment=0.;for(auto vertex:result.exact_set)if(!polygon_contains(result.numerical_set,vertex))containment=1.;add_check("numerical outer set contains analytical feasible set",containment,0.);double velocity_radius=0.;for(auto vertex:result.exact_set)velocity_radius=std::max(velocity_radius,std::abs(vertex[1]-result.center[1]));if(proper)add_check("proper schedule narrows velocity uncertainty",velocity_radius<.5*prior_u[1]?0.:1.,0.);else add_check("poor schedule retains full prior velocity range",std::abs(velocity_radius-prior_u[1]),2e-13);summary<<"discrete_rank_indicator,"<<spectrum.rank<<'\n'<<"discrete_small_information_eigenvalue,"<<spectrum.eigenvalues[1]<<'\n'<<"exact_rank,"<<analytic_rank.rank<<'\n'<<"analytical_velocity_sensitivity_max,"<<exact_sensitivity<<'\n'<<"analytical_fit_q0_m,"<<exact_fit[0]<<'\n'<<"analytical_fit_v0_m_s,"<<exact_fit[1]<<'\n'<<"discrete_fit_q0_m,"<<xi[0]<<'\n'<<"discrete_fit_v0_m_s,"<<xi[1]<<'\n'<<"qualification_representative_q0_m,"<<physical(result.center)[0]<<'\n'<<"qualification_representative_v0_m_s,"<<physical(result.center)[1]<<'\n'<<"retained_velocity_half_width_m_s,"<<velocity_radius*control_scale[1]<<'\n';
            auto trajectory=sys.forward(opt.control);auto qualification_trajectory=sys.forward(result.center);auto path=oq::affine(sys);auto E=allowance(sys);states<<std::setprecision(17)<<"elapsed_s,true_q_m,true_v_m_s,discrete_q_m,discrete_v_m_s,qualification_center_q_m,qualification_center_v_m_s,displacement_bound_m,velocity_bound_m_s,numerical_allowance\n";
            for(std::size_t k=0;k<sys.grid().time.size();++k){double tau=sys.grid().time[k]-t0;auto xt=exact(tau,truth),xd=trajectory.state[k],xc=qualification_trajectory.state[k];double bq=0.,bv=0.;auto gq=oq::goal_gradient(path,k,{1.,0.}),gv=oq::goal_gradient(path,k,{0.,1.});for(auto u:result.numerical_set){bq=std::max(bq,std::abs(dot(gq,u-result.center)));bv=std::max(bv,std::abs(dot(gv,u-result.center)));}states<<tau<<','<<xt[0]<<','<<xt[1]<<','<<xd[0]<<','<<xd[1]<<','<<xc[0]<<','<<xc[1]<<','<<bq+oq::goal_dual_norm({1.,0.})*E[k]<<','<<bv+oq::goal_dual_norm({0.,1.})*E[k]<<','<<E[k]<<'\n';}
            auto write_set=[&](const char* name,const ConvexPolygon2& p){for(std::size_t j=0;j<p.size();++j){auto v=physical(p[j]);sets<<name<<','<<j<<','<<v[0]<<','<<v[1]<<'\n';}};write_set("analytical",result.exact_set);write_set("numerical_outer",result.numerical_set);
            for(std::size_t k=0;k<3;++k){auto t=result.times[k];if(proper)add_check(std::string(oq::names[k])+" qualifies with proper timing",t&&std::abs(*t)<1e-14?0.:1.,0.);else add_check(std::string(oq::names[k])+" remains unqualified",t?1.:0.,0.);qualification<<oq::names[k]<<','<<oq::units[k]<<','<<oq::tolerances[k]<<',';if(t)qualification<<*t;qualification<<','<<static_cast<bool>(t)<<',';if(result.joint)qualification<<*result.joint;qualification<<','<<static_cast<bool>(result.joint)<<'\n';std::cout<<"  "<<oq::names[k]<<": "<<(t?"QUALIFIED":"NOT QUALIFIED")<< (t?" from "+std::to_string(*t)+" s":" over 8 s")<<'\n';}
            if(proper)add_check("proper schedule establishes joint qualification at zero",result.joint?*result.joint:1.,0.);else add_check("joint qualification remains unestablished",result.joint?1.:0.,0.);
            std::cout<<"  joint decision: "<<(result.joint?"QUALIFIED":"NOT QUALIFIED")<<'\n';
        }
    }
    for(const auto& r:fit_rows)std::cout<<"  N="<<static_cast<std::size_t>(r[0])<<": discrete fit (q0,v0)=("<<r[1]<<", "<<r[2]<<"), small scaled eigenvalue="<<r[3]<<", fixed-state error="<<r[4]<<'\n';
    summary<<"checks_passed,"<<checks<<'\n';std::cout<<"Checks passed: "<<checks<<". All qualification outputs are C++20 calculations.\n"<<"Convergence table: "<<(out/"convergence.csv")<<"\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
