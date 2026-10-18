#include "decay_qualification_model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace dq=decay_qualification;
using namespace theta_opt;
namespace {
struct Checks {
    std::ofstream file;std::size_t count=0;
    explicit Checks(const std::filesystem::path& path):file(path){file<<std::setprecision(17)<<"name,error,tolerance,passed\n";}
    void check(const std::string& name,double error,double tol){
        bool passed=std::isfinite(error)&&error<=tol;file<<name<<','<<error<<','<<tol<<','<<passed<<'\n';++count;
        if(!passed)throw std::runtime_error("Qualification check failed: "+name+" error="+std::to_string(error));
    }
};
struct BoundedResult {
    double fit{};dq::Interval outer{};std::optional<dq::Interval> raw;
    std::vector<Vec<2>> nominal,allowance,raw_error,physical_error;
    std::vector<Vec<2>> continuous_error;
    std::array<std::optional<double>,2> raw_time,physical_time,continuous_time;
};
BoundedResult bounded_experiment(const System<2,1>& system,std::optional<dq::Exposure> extra={}){
    BoundedResult result;const auto& grid=system.grid();
    auto optimum=dq::fit(system);result.fit=optimum.control[0];result.nominal=system.forward(optimum.control).state;
    const dq::Interval prior{dq::bounded[0]-dq::delta,dq::bounded[0]+dq::delta};
    result.allowance=dq::physical_allowance(system,prior);
    auto forward=[&](double xi){return system.forward({xi}).state;};
    result.raw=dq::consistent_interval(forward,grid,dq::bounded,dq::constant_bounds(),extra);
    auto error_bounds=dq::constant_bounds();
    for(std::size_t j=0;j<dq::elapsed.size();++j)error_bounds[j]+=result.allowance[dq::node_at(grid,dq::elapsed[j])][0];
    auto inflated_extra=extra;
    if(inflated_extra)inflated_extra->bound+=result.allowance[dq::node_at(grid,inflated_extra->elapsed)][1];
    auto outer=dq::consistent_interval(forward,grid,dq::bounded,error_bounds,inflated_extra);
    if(!outer)throw std::runtime_error("Error-aware data set unexpectedly empty");
    result.outer=*outer;
    result.physical_error=dq::envelope(result.nominal,forward(outer->lower),forward(outer->upper));
    for(std::size_t n=0;n<grid.time.size();++n)result.physical_error[n]=result.physical_error[n]+result.allowance[n];
    result.physical_time=dq::qualify(grid,result.physical_error);
    // Enclose the full interval, relative to the piecewise-linear nominal
    // reconstruction. Node-only qualification does not cover between nodes.
    auto discrete_spread=dq::envelope(result.nominal,forward(outer->lower),forward(outer->upper));
    auto high=forward(prior.upper);result.continuous_error.resize(grid.steps());
    for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n],th=grid.theta[n];
        double x0=high[n][0],x1=high[n+1][0],g0=decay::a*x0+decay::b*x0*x0,g1=decay::a*x1+decay::b*x1*x1;
        double D=std::max(th,1-th)*std::abs(g0-g1)+dq::solver_options.absolute_tolerance/h;
        double phi=-std::expm1(-decay::a*h)/decay::a;
        double integral=result.allowance[n][0]*phi+D/decay::a*(h-phi);
        double between_q=result.allowance[n][1]+integral+h*std::max(th,1-th)*std::abs(x1-x0)+dq::solver_options.absolute_tolerance;
        result.continuous_error[n]={
            std::max(discrete_spread[n][0],discrete_spread[n+1][0])+std::max(result.allowance[n][0],result.allowance[n+1][0]),
            std::max(discrete_spread[n][1],discrete_spread[n+1][1])+between_q};
    }
    std::vector<double> starts(grid.time.begin(),grid.time.end()-1);
    for(std::size_t j=0;j<2;++j){std::vector<double> ratios;for(auto e:result.continuous_error)ratios.push_back(e[j]/dq::tolerance[j]);
        result.continuous_time[j]=persistent_node_qualification(starts,ratios);if(result.continuous_time[j])*result.continuous_time[j]-=decay::t0;}
    if(result.raw){result.raw_error=dq::envelope(result.nominal,forward(result.raw->lower),forward(result.raw->upper));result.raw_time=dq::qualify(grid,result.raw_error);}
    return result;
}
double subset_error(dq::Interval outer,dq::Interval inner){return std::max({0.,outer.lower-inner.lower,inner.upper-outer.upper});}
double information(const System<2,1>& system,double fit,Checks& checks,const std::string& label){
    auto trajectory=system.forward({fit});auto tangent=system.tangent({fit},trajectory,{1.});
    std::vector<Vec<2>> sources(trajectory.state.size());double I=0.;
    for(double offset:dq::elapsed){auto n=dq::node_at(system.grid(),offset);double s=tangent[n][0];I+=s*s/(dq::sigma*dq::sigma);sources[n][0]+=s/(dq::sigma*dq::sigma);}
    double adjoint_I=system.state_pullback({fit},trajectory,sources)[0];
    checks.check(label+" information tangent versus adjoint",std::abs(adjoint_I-I),1e-8);return adjoint_I;
}
double qualification_exact(double center,dq::Interval interval,std::size_t component){
    auto err=[&](double elapsed){auto nominal=decay::exact(elapsed,center).state,lo=decay::exact(elapsed,interval.lower).state,hi=decay::exact(elapsed,interval.upper).state;
        return std::max(std::abs(lo[component]-nominal[component]),std::abs(hi[component]-nominal[component]));};
    if(component==1)return err(10.)<=dq::tolerance[1]?0.:std::numeric_limits<double>::quiet_NaN();
    if(err(0.)<=dq::tolerance[0])return 0.;
    if(err(10.)>dq::tolerance[0])return std::numeric_limits<double>::quiet_NaN();
    return decay::root([&](double t){return err(t)-dq::tolerance[0];},0.,10.);
}
void validate_bounded(const System<2,1>& system,const BoundedResult& r,dq::Interval exact,
                      Checks& checks,const std::string& label){
    checks.check(label+" corrected set contains exact set",subset_error(r.outer,exact),2e-9);
    double excess=0.,bound_excess=0.;const auto& grid=system.grid();
    for(int k=0;k<5;++k){double initial=dq::bounded[0]-dq::delta+2*dq::delta*k/4.;auto path=system.forward({initial});
        for(std::size_t n=0;n<grid.time.size();++n){auto reference=decay::exact(grid.time[n]-decay::t0,initial).state;
            for(std::size_t j=0;j<2;++j)bound_excess=std::max(bound_excess,std::abs(path.state[n][j]-reference[j])-r.allowance[n][j]);}}
    for(std::size_t n=0;n<grid.time.size();++n){auto reference=decay::exact(grid.time[n]-decay::t0,dq::truth).state;
        for(std::size_t j=0;j<2;++j)excess=std::max(excess,std::abs(r.nominal[n][j]-reference[j])-r.physical_error[n][j]);}
    checks.check(label+" physical bound against analytical paths",bound_excess,2e-10);
    checks.check(label+" full goal allowance covers truth",excess,2e-9);
    double continuous_excess=0.;
    for(std::size_t n=0;n<grid.steps();++n)for(double fraction:{.25,.5,.75}){
        double t=(1-fraction)*grid.time[n]+fraction*grid.time[n+1];auto reference=decay::exact(t-decay::t0,dq::truth).state;
        auto reconstruction=(1-fraction)*r.nominal[n]+fraction*r.nominal[n+1];
        for(std::size_t j=0;j<2;++j)continuous_excess=std::max(continuous_excess,std::abs(reconstruction[j]-reference[j])-r.continuous_error[n][j]);
    }
    checks.check(label+" continuous interval envelope against analytical truth",continuous_excess,2e-9);
}
}

int main(int argc,char** argv){try{
    namespace fs=std::filesystem;fs::path output=argc>1?argv[1]:"results/qualification";fs::create_directories(output);
    Checks checks(output/"checks.csv");std::ofstream summary(output/"summary.csv");summary<<std::setprecision(17)<<"quantity,value\n";
    auto report=[&](const std::string& key,double value){summary<<key<<','<<value<<'\n';std::cout<<std::setprecision(12)<<key<<": "<<value<<'\n';};
    const double exact_gaussian_fit=dq::analytical_fit(dq::gaussian),exact_bounded_fit=dq::analytical_fit(dq::bounded);
    const auto exact_set=dq::exact_interval();const dq::Exposure extra;
    const auto added_exact_set=dq::exact_interval(extra);const double added_exact_fit=dq::analytical_fit(dq::bounded,extra);
    report("analytical_gaussian_fit",exact_gaussian_fit);report("analytical_bounded_fit",exact_bounded_fit);
    report("exact_interval_lower",exact_set.lower);report("exact_interval_upper",exact_set.upper);
    report("analytical_qualification_x_min",qualification_exact(exact_bounded_fit,exact_set,0));
    report("analytical_qualification_q_min",qualification_exact(exact_bounded_fit,exact_set,1));
    checks.check("exact bounded interval contains truth",std::max({0.,exact_set.lower-dq::truth,dq::truth-exact_set.upper}),1e-12);
    // Baseline primary grid matches the earlier nonlinear qualification experiment.
    auto primary_grid=dq::mesh(16,.65);auto gaussian_system=dq::system(primary_grid,dq::gaussian);
    auto gaussian_opt=dq::fit(gaussian_system);double gaussian_fit=gaussian_opt.control[0];
    auto gaussian_evaluation=gaussian_system.evaluate({gaussian_fit});
    auto gaussian_tangent=gaussian_system.tangent({gaussian_fit},gaussian_evaluation.trajectory,{1.});
    double I=information(gaussian_system,gaussian_fit,checks,"primary");double P=1./I;
    report("discrete_gaussian_fit",gaussian_fit);report("local_initial_variance",P);
    report("gaussian_gradient_norm",norm(gaussian_evaluation.gradient));
    checks.check("Gaussian adjoint optimizer stationarity",norm(gaussian_evaluation.gradient),1e-8);
    double gaussian_radius=independent_gaussian_box_radius(dq::sigma,dq::elapsed.size(),.95);
    auto gaussian_exact_set=dq::consistent_interval([&](double xi){return dq::exact_path(primary_grid,xi);},primary_grid,dq::gaussian,dq::constant_bounds(gaussian_radius)).value();
    auto gaussian_numerical_allowance=dq::physical_allowance(gaussian_system,{dq::gaussian[0]-gaussian_radius,dq::gaussian[0]+gaussian_radius});
    auto gaussian_bounds=dq::constant_bounds(gaussian_radius);
    for(std::size_t j=0;j<dq::elapsed.size();++j)gaussian_bounds[j]+=gaussian_numerical_allowance[dq::node_at(primary_grid,dq::elapsed[j])][0];
    auto gaussian_outer=dq::consistent_interval([&](double xi){return gaussian_system.forward({xi}).state;},primary_grid,dq::gaussian,gaussian_bounds).value();
    auto finite95_error=dq::envelope(gaussian_evaluation.trajectory.state,gaussian_system.forward({gaussian_outer.lower}).state,gaussian_system.forward({gaussian_outer.upper}).state);
    for(std::size_t n=0;n<finite95_error.size();++n)finite95_error[n]=finite95_error[n]+gaussian_numerical_allowance[n];
    checks.check("finite-noise confidence outer set contains analytical set",subset_error(gaussian_outer,gaussian_exact_set),2e-9);
    report("joint95_noise_radius",gaussian_radius);report("joint95_exact_lower",gaussian_exact_set.lower);report("joint95_exact_upper",gaussian_exact_set.upper);
    report("joint95_outer_lower",gaussian_outer.lower);report("joint95_outer_upper",gaussian_outer.upper);
    auto bounded_system=dq::system(primary_grid,dq::bounded);auto primary=bounded_experiment(bounded_system);
    validate_bounded(bounded_system,primary,exact_set,checks,"primary");
    report("discrete_bounded_fit",primary.fit);report("raw_interval_lower",primary.raw->lower);report("raw_interval_upper",primary.raw->upper);
    report("outer_interval_lower",primary.outer.lower);report("outer_interval_upper",primary.outer.upper);
    report("raw_node_qualification_x_min",dq::output(primary.raw_time[0]));
    report("physical_node_qualification_x_min",dq::output(primary.physical_time[0]));report("physical_node_qualification_q_min",dq::output(primary.physical_time[1]));
    report("physical_continuous_qualification_x_min",dq::output(primary.continuous_time[0]));report("physical_continuous_qualification_q_min",dq::output(primary.continuous_time[1]));
    checks.check("concentration-only data do not qualify accumulated exposure",primary.continuous_time[1]?1.:0.,0.);
    checks.check("raw discrete interval excludes truth on this grid",primary.raw->upper<dq::truth?0.:1.,0.);

    std::ofstream curves(output/"curves.csv");curves<<std::setprecision(17)<<"elapsed,x_nominal,q_nominal,truth_x,truth_q,local95_x,local95_q,exact_bound_x,exact_bound_q,raw_bound_x,raw_bound_q,physical_bound_x,physical_bound_q,numerical_allowance_x,numerical_allowance_q,linear_bound_x,linear_bound_q,taylor_bound_x,taylor_bound_q,joint95_physical_bound_x,joint95_physical_bound_q\n";
    auto exact_nominal=dq::exact_path(primary_grid,exact_bounded_fit),exact_lower=dq::exact_path(primary_grid,exact_set.lower),exact_upper=dq::exact_path(primary_grid,exact_set.upper);
    auto exact_error=dq::envelope(exact_nominal,exact_lower,exact_upper);
    std::vector<Vec<2>> local95(primary_grid.time.size());double goal_difference=0.,curvature_excess=0.;
    auto bounded_trajectory=bounded_system.forward({primary.fit});auto bounded_tangent=bounded_system.tangent({primary.fit},bounded_trajectory,{1.});
    double radius=std::max(std::abs(primary.raw->lower-primary.fit),std::abs(primary.raw->upper-primary.fit));
    auto curvature=dq::curvature_allowance(bounded_system,{std::min(primary.raw->lower,primary.fit),std::max(primary.raw->upper,primary.fit)});
    for(std::size_t n=0;n<primary_grid.time.size();++n){
        Vec<2> confidence{},linear{},taylor{};
        for(std::size_t j=0;j<2;++j){Vec<2> component{};component[j]=1.;
            double goal=gaussian_system.goal_gradient({gaussian_fit},gaussian_evaluation.trajectory,n,component)[0];
            goal_difference=std::max(goal_difference,std::abs(goal-gaussian_tangent[n][j]));
            confidence[j]=1.959963984540054*std::sqrt(local_goal_variance(Vec<1>{goal},Mat<1,1>{{{P}}}));
            linear[j]=box_support(Vec<1>{bounded_tangent[n][j]},Vec<1>{radius});
            taylor[j]=linear[j]+.5*curvature[n][j]*radius*radius;
            curvature_excess=std::max(curvature_excess,primary.raw_error[n][j]-taylor[j]);
        }
        local95[n]=confidence;auto truth=decay::exact(primary_grid.time[n]-decay::t0,dq::truth).state;
        curves<<primary_grid.time[n]-decay::t0;
        for(const auto& pair:{gaussian_evaluation.trajectory.state[n],truth,confidence,exact_error[n],primary.raw_error[n],primary.physical_error[n],primary.allowance[n],linear,taylor,finite95_error[n]})curves<<','<<pair[0]<<','<<pair[1];
        curves<<'\n';
    }
    checks.check("all goal adjoints versus tangent",goal_difference,1e-11);
    checks.check("nonlinear interval envelope covered by linear plus curvature",curvature_excess,2e-10);
    auto local_time=dq::qualify(primary_grid,local95);
    auto joint95_time=dq::qualify(primary_grid,finite95_error);
    report("local_gaussian_node_qualification_x_min",dq::output(local_time[0]));report("local_gaussian_node_qualification_q_min",dq::output(local_time[1]));
    // A real additional exposure observation: reoptimize and reconstruct its
    // nonlinear consistency set. Refine enough to make its error budget useful.
    auto added_grid=dq::mesh(512,.5);auto added_system=dq::system(added_grid,dq::bounded,extra);
    auto added=bounded_experiment(added_system,extra);validate_bounded(added_system,added,added_exact_set,checks,"added exposure");
    report("added_analytical_fit",added_exact_fit);report("added_discrete_fit",added.fit);
    report("added_exact_lower",added_exact_set.lower);report("added_exact_upper",added_exact_set.upper);
    report("added_outer_lower",added.outer.lower);report("added_outer_upper",added.outer.upper);
    report("added_analytical_qualification_x_min",qualification_exact(added_exact_fit,added_exact_set,0));report("added_analytical_qualification_q_min",qualification_exact(added_exact_fit,added_exact_set,1));
    report("added_physical_node_qualification_x_min",dq::output(added.physical_time[0]));report("added_physical_node_qualification_q_min",dq::output(added.physical_time[1]));
    report("added_physical_continuous_qualification_x_min",dq::output(added.continuous_time[0]));report("added_physical_continuous_qualification_q_min",dq::output(added.continuous_time[1]));
    checks.check("added exposure retains synthetic truth",std::max({0.,added_exact_set.lower-dq::truth,dq::truth-added_exact_set.upper}),1e-12);
    checks.check("added exposure concentration qualified over whole horizon",added.continuous_time[0]&&*added.continuous_time[0]==0.?0.:1.,0.);
    checks.check("added exposure accumulated goal qualified over whole horizon",added.continuous_time[1]&&*added.continuous_time[1]==0.?0.:1.,0.);
    checks.check("added exposure adjoint optimizer stationarity",norm(added_system.evaluate({added.fit}).gradient),1e-8);
    std::ofstream added_curves(output/"added_exposure.csv");added_curves<<std::setprecision(17)<<"elapsed,physical_bound_x,physical_bound_q,exact_bound_x,exact_bound_q\n";
    auto added_exact_error=dq::envelope(dq::exact_path(added_grid,added_exact_fit),dq::exact_path(added_grid,added_exact_set.lower),dq::exact_path(added_grid,added_exact_set.upper));
    for(std::size_t n=0;n<added_grid.time.size();++n)added_curves<<added_grid.time[n]-decay::t0<<','<<added.physical_error[n][0]<<','<<added.physical_error[n][1]<<','<<added_exact_error[n][0]<<','<<added_exact_error[n][1]<<'\n';
    std::ofstream continuous_curves(output/"continuous_envelopes.csv");continuous_curves<<std::setprecision(17)<<"case,start,end,physical_bound_x,physical_bound_q\n";
    for(bool augmented:{false,true}){const auto& grid=augmented?added_grid:primary_grid;const auto& result=augmented?added:primary;
        for(std::size_t n=0;n<grid.steps();++n)continuous_curves<<(augmented?"added":"baseline")<<','<<grid.time[n]-decay::t0<<','<<grid.time[n+1]-decay::t0<<','<<result.continuous_error[n][0]<<','<<result.continuous_error[n][1]<<'\n';}

    // All derivative checks use nonstationary controls. Added exposure sources
    // must be differentiated in both the q and concentration adjoint channels.
    for(bool augmented:{false,true}){auto test_system=augmented?dq::system(primary_grid,dq::bounded,extra):gaussian_system;double xi=2.4,eps=1e-4;
        auto e=test_system.evaluate({xi}),ep=test_system.evaluate({xi+eps}),em=test_system.evaluate({xi-eps});double H=test_system.hessian_vector({xi},{1.})[0];
        checks.check(augmented?"added observation adjoint gradient FD":"nodal observation adjoint gradient FD",std::abs((ep.cost-em.cost)/(2*eps)-e.gradient[0])/std::max(1.,std::abs(e.gradient[0])),2e-8);
        checks.check(augmented?"added observation incremental Hessian FD":"nodal observation incremental Hessian FD",std::abs((ep.gradient[0]-em.gradient[0])/(2*eps)-H)/std::max(1.,std::abs(H)),2e-8);
        checks.check("standalone observation quadrature cost",std::abs(test_system.cost({xi},e.trajectory)-e.cost),1e-10);
    }

    std::ofstream refinement(output/"refinement.csv");refinement<<std::setprecision(17)<<"theta,subdivisions,h,state_error,cost_error,tangent_error,goal_error,gradient_error,hessian_error,control_error,raw_lower,raw_upper,outer_lower,outer_upper,raw_qualification_x,physical_qualification_x,physical_qualification_q\n";
    for(double theta:{.5,.65}){
        std::array<double,7> previous{},last{};
        for(std::size_t k:{1U,2U,4U,8U,16U,32U,64U,128U}){
            auto grid=dq::mesh(k,theta);auto sys=dq::system(grid,dq::gaussian);auto opt=dq::fit(sys);double xi=2.4;
            auto e=sys.evaluate({xi});auto s=sys.tangent({xi},e.trajectory,{1.});auto reference=dq::analytical_objective(xi,dq::gaussian);
            double state_error=0.,tangent_error=0.;for(std::size_t n=0;n<grid.time.size();++n){auto ex=decay::exact(grid.time[n]-decay::t0,xi);state_error=std::max(state_error,norm(e.trajectory.state[n]-ex.state));tangent_error=std::max(tangent_error,norm(s[n]-ex.tangent));}
            auto goal=sys.goal_gradient({xi},e.trajectory,grid.steps(),{0.,1.});double goal_error=std::abs(goal[0]-decay::exact(10.,xi).tangent[1]);
            previous=last;last={state_error,std::abs(e.cost-reference[0]),tangent_error,goal_error,std::abs(e.gradient[0]-reference[1]),std::abs(sys.hessian_vector({xi},{1.})[0]-reference[2]),std::abs(opt.control[0]-exact_gaussian_fit)};
            auto bs=dq::system(grid,dq::bounded);auto r=bounded_experiment(bs);validate_bounded(bs,r,exact_set,checks,"refined "+std::to_string(theta)+"/"+std::to_string(k));
            checks.check("refined observation stationarity",norm(sys.evaluate(opt.control).gradient),1e-8);
            double h=0.;for(std::size_t n=0;n<grid.steps();++n)h=std::max(h,grid.time[n+1]-grid.time[n]);
            refinement<<theta<<','<<k<<','<<h;for(double err:last)refinement<<','<<err;
            double nan=std::numeric_limits<double>::quiet_NaN();refinement<<','<<(r.raw?r.raw->lower:nan)<<','<<(r.raw?r.raw->upper:nan)<<','<<r.outer.lower<<','<<r.outer.upper<<','<<dq::output(r.raw_time[0])<<','<<dq::output(r.physical_time[0])<<','<<dq::output(r.physical_time[1])<<'\n';
        }
        for(std::size_t j=0;j<last.size();++j){double p=std::log2(previous[j]/last[j]);report("order_"+std::to_string(theta)+"_component_"+std::to_string(j),p);checks.check("observation problem temporal order",std::abs(p-(theta==.5?2.:1.)),.13);}
    }
    std::ofstream observations(output/"observations.csv");observations<<std::setprecision(17)<<"elapsed,truth_x,gaussian,bounded\n";
    for(std::size_t j=0;j<dq::elapsed.size();++j)observations<<dq::elapsed[j]<<','<<decay::exact(dq::elapsed[j],dq::truth).state[0]<<','<<dq::gaussian[j]<<','<<dq::bounded[j]<<'\n';
    report("qualification_checks_passed",static_cast<double>(checks.count));
    // One authoritative C++ report shared by the console and plot renderer.
    // Absence of a time is a qualification outcome, not a numerical failure.
    struct QualificationRow {
        std::string name,criterion;
        std::array<std::optional<double>,2> time;
    };
    auto optional_time=[](double t)->std::optional<double>{return std::isfinite(t)?std::optional<double>{t}:std::nullopt;};
    const std::vector<QualificationRow> rows{
        {"baseline_analytical","continuous",{optional_time(qualification_exact(exact_bounded_fit,exact_set,0)),optional_time(qualification_exact(exact_bounded_fit,exact_set,1))}},
        {"baseline_raw_discrete","nodes_uncorrected",primary.raw_time},
        {"baseline_physical_nodes","nodes",primary.physical_time},
        {"baseline_physical_continuous","continuous",primary.continuous_time},
        {"local_gaussian","nodes_local95",local_time},
        {"joint95_gaussian","nodes_joint95",joint95_time},
        {"added_analytical","continuous",{optional_time(qualification_exact(added_exact_fit,added_exact_set,0)),optional_time(qualification_exact(added_exact_fit,added_exact_set,1))}},
        {"added_physical_nodes","nodes",added.physical_time},
        {"added_physical_continuous","continuous",added.continuous_time}
    };
    std::ofstream times(output/"qualification_times.csv");
    if(!times)throw std::runtime_error("Cannot open qualification-time report");
    times<<std::setprecision(17)<<"case,criterion,goal,status,elapsed_min,absolute_time_min\n";
    std::cout<<"\nPersistent state qualification over 10 elapsed minutes\n"
             <<"Times are elapsed from t0; unqualified means no persistent qualification.\n"
             <<std::left<<std::setw(31)<<"Case"<<std::setw(21)<<"Criterion"
             <<std::setw(17)<<"x time (min)"<<"q time (min)\n";
    for(const auto& row:rows){
        std::cout<<std::left<<std::setw(31)<<row.name<<std::setw(21)<<row.criterion;
        for(std::size_t j=0;j<2;++j){
            times<<row.name<<','<<row.criterion<<','<<(j==0?"x":"q")<<','
                 <<(row.time[j]?"qualified":"unqualified")<<',';
            if(row.time[j])times<<*row.time[j]<<','<<decay::t0+*row.time[j];else times<<',';
            times<<'\n';
            if(row.time[j])std::cout<<std::setw(17)<<std::fixed<<std::setprecision(6)<<*row.time[j];
            else std::cout<<std::setw(17)<<"unqualified";
        }
        std::cout<<'\n';
    }
    std::cout<<"Raw discrete times ignore numerical model error. Local95 is a linear indicator.\n";
    // Physical determinations use only the continuous, numerical-error-aware
    // envelopes. Analytical references and local/node-only diagnostics above
    // do not vote on this verdict.
    std::ofstream decisions(output/"qualification_decisions.csv");
    if(!decisions)throw std::runtime_error("Cannot open qualification decisions");
    decisions<<std::setprecision(17)<<"case,role,concentration_status,concentration_elapsed_min,exposure_status,exposure_elapsed_min,all_goals_status,all_goals_elapsed_min,blocking_goals\n";
    std::cout<<"\nPHYSICAL QUALIFICATION DETERMINATION\n"
             <<"Requested goals: concentration error <= 0.02; exposure error <= 0.025 min.\n"
             <<"Both cases use bounded-error data and include numerical model error.\n"
             <<"Criterion: persistent from the reported time through elapsed 10 min, including between nodes.\n"
             <<"Times below are elapsed from t0.\n";
    for(bool augmented:{false,true}){
        const auto& time=augmented?added.continuous_time:primary.continuous_time;
        auto joint=all_goals_qualification_time(time);
        const std::string name=augmented?"added_exposure":"baseline";
        std::string blockers;
        if(!time[0])blockers="concentration";
        if(!time[1])blockers+=(blockers.empty()?"":"+")+std::string("exposure");
        decisions<<name<<",physical_continuous,";
        for(auto t:time){decisions<<(t?"qualified":"not_qualified")<<',';if(t)decisions<<*t;decisions<<',';}
        decisions<<(joint?"qualified":"not_qualified")<<',';if(joint)decisions<<*joint;
        decisions<<','<<blockers<<'\n';
        std::cout<<'\n'<<(augmented?"ADDED EXPOSURE MEASUREMENT (refined grid)":"BASELINE: CONCENTRATION MEASUREMENTS ONLY")<<'\n';
        for(std::size_t j=0;j<2;++j){
            std::cout<<"  "<<(j==0?"Concentration":"Exposure")<<": ";
            if(time[j])std::cout<<"QUALIFIED from "<<std::fixed<<std::setprecision(6)<<*time[j]<<" min through the end of the horizon.\n";
            else std::cout<<"NOT QUALIFIED over the requested horizon.\n";
        }
        std::cout<<"  ALL REQUESTED GOALS: ";
        if(joint)std::cout<<"QUALIFIED from "<<*joint<<" min through the end of the horizon.\n";
        else std::cout<<"NOT QUALIFIED. Blocking goal: "<<blockers<<".\n";
    }
    std::cout<<"\nQualification test checks: PASS ("<<checks.count<<" checks).\n"
             <<"Passing tests verify the expected decisions; they do not imply every case qualifies.\n"
             <<"NOT QUALIFIED means the stated data and error budgets do not establish the tolerance.\n"
             <<"Added-measurement qualification is retrospective, using all observations.\n";
    std::cout<<"C++20 state qualification complete. Results: "<<output<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
