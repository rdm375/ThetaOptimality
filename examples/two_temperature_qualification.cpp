// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#include "two_temperature_model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace theta_opt;
namespace {
struct Checks {
    std::ofstream csv;std::size_t count=0;
    explicit Checks(const std::filesystem::path& p):csv(p){if(!csv)throw std::runtime_error("Cannot write thermal checks");csv<<std::setprecision(17)<<"name,error,tolerance,passed\n";}
    void check(const std::string& name,double error,double tol){bool pass=std::isfinite(error)&&error<=tol;csv<<name<<','<<error<<','<<tol<<','<<pass<<'\n';++count;if(!pass)throw std::runtime_error("Thermal check failed: "+name+" error="+std::to_string(error));}
};
double matrix_error(const Mat<2,2>& a,const Mat<2,2>& b){return std::max(norm(a[0]-b[0]),norm(a[1]-b[1]));}
struct Case {
    Vec<2> fit{};Mat<2,2> covariance{};thermal::Polygon polygon;
    std::vector<Vec<2>> state,ellipse,physical,continuous;
    std::vector<double> numerical;
    std::array<std::optional<double>,2> ellipse_time,physical_time;
    std::vector<Iteration<2>> optimizer_trace;
};
std::array<std::optional<double>,2> qualify(const std::vector<double>& time,const std::vector<Vec<2>>& error){
    std::array<std::optional<double>,2> result;
    for(std::size_t j=0;j<2;++j){std::vector<double> ratios;for(auto e:error)ratios.push_back(e[j]/thermal::tolerance[j]);result[j]=persistent_node_qualification(time,ratios);if(result[j])*result[j]-=thermal::t0;}return result;
}
Case experiment(const System<2,2>& sys,bool wall,Checks& checks,const std::string& label,bool check_all_adjoint=false){
    Case result;auto opt=thermal::fit(sys);result.fit=opt.control;result.optimizer_trace=opt.trace;auto evaluation=sys.evaluate(opt.control);result.state=evaluation.trajectory.state;const auto& grid=sys.grid();
    checks.check(label+" optimized adjoint stationarity",norm(evaluation.gradient),1e-8);
    checks.check(label+" one quadratic Newton update",std::abs(static_cast<double>(opt.trace.size())-2.),0.);
    Mat<2,2> info{};std::array<std::vector<Vec<2>>,2> s;
    for(std::size_t j=0;j<2;++j){Vec<2> unit{};unit[j]=1.;s[j]=sys.tangent(opt.control,evaluation.trajectory,unit);
        std::vector<Vec<2>> sources(grid.time.size());
        for(double t:thermal::elapsed){auto n=thermal::node_at(grid,t);sources[n][0]+=s[j][n][0]/(thermal::sigma*thermal::sigma);}
        if(wall){auto n=thermal::node_at(grid,1.);sources[n][1]+=s[j][n][1]/(thermal::sigma*thermal::sigma);}
        auto col=sys.state_pullback(opt.control,evaluation.trajectory,sources),hcol=sys.hessian_vector(opt.control,unit);
        checks.check(label+" information action equals Hessian",norm(col-hcol),1e-10);
        for(std::size_t i=0;i<2;++i)info[i][j]=col[i];
    }
    checks.check(label+" information symmetry",std::abs(info[0][1]-info[1][0]),1e-10);
    result.covariance=thermal::inverse(info);
    result.numerical=thermal::physical_allowance(sys);
    auto affine=thermal::affine(sys);
    double radius=independent_gaussian_box_radius(thermal::sigma,wall?9:8,.95);
    result.polygon=thermal::confidence_polygon(grid,affine,radius,result.numerical,wall);
    auto exact_set=thermal::analytical_polygon(grid,radius,wall);
    checks.check(label+" exact confidence polygon contains truth",thermal::contains(exact_set,thermal::truth)?0.:1.,0.);
    for(auto vertex:exact_set)checks.check(label+" numerical polygon contains analytical vertex",thermal::contains(result.polygon,vertex)?0.:1.,0.);
    result.ellipse.resize(grid.time.size());result.physical.resize(grid.time.size());
    std::vector<Vec<2>> spread(grid.time.size());double goal_error=0.;
    for(std::size_t n=0;n<grid.time.size();++n)for(std::size_t j=0;j<2;++j){
        Vec<2> g=thermal::affine_goal_gradient(affine,n,thermal::goals[j]);
        if(check_all_adjoint||n==0||n==grid.steps()/2||n==grid.steps()){
            auto adj=sys.goal_gradient(opt.control,evaluation.trajectory,n,thermal::goals[j]);goal_error=std::max(goal_error,norm(adj-g));g=adj;
        }
        result.ellipse[n][j]=std::sqrt(-2*std::log(.05))*std::sqrt(local_goal_variance(g,result.covariance));
        for(auto vertex:result.polygon)spread[n][j]=std::max(spread[n][j],std::abs(dot(g,vertex-result.fit)));
        double gain=std::abs(thermal::goals[j][0])+std::abs(thermal::goals[j][1]);
        result.physical[n][j]=spread[n][j]+gain*result.numerical[n];
    }
    checks.check(label+" goal adjoints versus tangents",goal_error,2e-10);
    result.ellipse_time=qualify(grid.time,result.ellipse);
    result.continuous.resize(grid.steps());
    for(std::size_t n=0;n<grid.steps();++n)for(std::size_t j=0;j<2;++j){double gain=std::abs(thermal::goals[j][0])+std::abs(thermal::goals[j][1]);result.continuous[n][j]=std::max(spread[n][j],spread[n+1][j])+gain*std::max(result.numerical[n],result.numerical[n+1]);}
    result.physical_time=qualify(std::vector<double>(grid.time.begin(),grid.time.end()-1),result.continuous);
    double envelope_excess=0.,defect_excess=0.,balance=0.;
    for(auto initial:thermal::prior_polygon()){auto path=sys.forward(initial);
        for(std::size_t n=0;n<grid.time.size();++n)defect_excess=std::max(defect_excess,norm(path.state[n]-thermal::exact(grid.time[n]-thermal::t0,initial))-result.numerical[n]);}
    for(std::size_t n=0;n<grid.steps();++n){double h=grid.time[n+1]-grid.time[n],th=grid.theta[n];
        double lhs=dot(thermal::goals[1],result.state[n+1]-result.state[n]);
        double old=result.state[n][0]+result.state[n][1]-50.,next=result.state[n+1][0]+result.state[n+1][1]-50.;
        balance=std::max(balance,std::abs(lhs+60*h*((1-th)*old+th*next)));
        for(double fraction:{.25,.5,.75}){double tau=(1-fraction)*grid.time[n]+fraction*grid.time[n+1]-thermal::t0;auto exact=thermal::exact(tau,thermal::truth),nominal=(1-fraction)*result.state[n]+fraction*result.state[n+1];
            for(std::size_t j=0;j<2;++j)envelope_excess=std::max(envelope_excess,std::abs(dot(thermal::goals[j],exact-nominal))-result.continuous[n][j]);}
    }
    checks.check(label+" conservative theta thermal-energy balance",balance,2e-9);
    checks.check(label+" computable defect bound versus analytical corners",defect_excess,2e-10);
    checks.check(label+" whole-step physical envelope covers truth",envelope_excess,2e-8);
    return result;
}
double value(std::optional<double> t){return t.value_or(std::numeric_limits<double>::quiet_NaN());}
}

int main(int argc,char** argv){try{
    namespace fs=std::filesystem;fs::path output=argc>1?argv[1]:"results/two_temperature";fs::create_directories(output);Checks checks(output/"checks.csv");
    auto reference=thermal::analytical_fit();std::array<std::optional<double>,2> exact_time;
    for(std::size_t j=0;j<2;++j)exact_time[j]=thermal::analytical_qualification(reference.covariance,j);
    auto coarse_grid=thermal::mesh(16,.65),fine_grid=thermal::mesh(256,.5);
    auto coarse_sys=thermal::system(coarse_grid),fine_sys=thermal::system(fine_grid),added_sys=thermal::system(fine_grid,true);
    auto coarse=experiment(coarse_sys,false,checks,"primary",true),fine=experiment(fine_sys,false,checks,"refined"),added=experiment(added_sys,true,checks,"wall sensor");
    checks.check("refined thermal case qualifies both goals",all_goals_qualification_time(fine.physical_time)?0.:1.,0.);
    checks.check("wall observation qualifies both goals",all_goals_qualification_time(added.physical_time)?0.:1.,0.);
    checks.check("wall observation advances joint qualification",*all_goals_qualification_time(added.physical_time)<*all_goals_qualification_time(fine.physical_time)?0.:1.,0.);
    // Independent conservative constant-mass step, matching the original
    // Python model rather than calling the library's normalized ODE update.
    const Mat<2,2> mass{{{1.,0.},{0.,4.}}},K{{{-2.,1.},{1.,-2.}}};Vec<2> independent=coarse.fit;double mass_error=0.;
    for(std::size_t n=0;n<coarse_grid.steps();++n){double h=coarse_grid.time[n+1]-coarse_grid.time[n],theta=coarse_grid.theta[n];
        auto B=add(mass,scaled(-h*theta,K)),C=add(mass,scaled(h*(1-theta),K));
        independent=LU<2>(B).solve(mv(C,independent)+h*Vec<2>{40.,10.});mass_error=std::max(mass_error,norm(independent-coarse.state[n+1]));}
    checks.check("library theta equals original conservative mass step",mass_error,2e-11);
    // Independent covariance update for the actually assimilated wall sensor.
    auto baseline_path=fine_sys.forward(fine.fit);auto wall_g=fine_sys.goal_gradient(fine.fit,baseline_path,thermal::node_at(fine_grid,1.),{0.,1.});
    auto v=mv(fine.covariance,wall_g);double denominator=thermal::sigma*thermal::sigma+dot(wall_g,v);auto update=add(fine.covariance,scaled(-1./denominator,thermal::outer(v)));
    checks.check("wall covariance update versus new information solve",matrix_error(update,added.covariance),2e-12);
    std::ofstream summary(output/"summary.csv");summary<<std::setprecision(17)<<"quantity,value\n";
    auto report=[&](const std::string& key,double v){summary<<key<<','<<v<<'\n';std::cout<<std::setprecision(12)<<key<<": "<<v<<'\n';};
    report("analytical_fit_fluid_C",reference.fit[0]);report("analytical_fit_wall_C",reference.fit[1]);
    report("discrete_fit_fluid_C",coarse.fit[0]);report("discrete_fit_wall_C",coarse.fit[1]);
    report("analytical_ellipse_fluid_qualification_min",value(exact_time[0]));report("analytical_ellipse_energy_qualification_min",value(exact_time[1]));
    report("primary_discrete_ellipse_fluid_node_min",value(coarse.ellipse_time[0]));report("primary_discrete_ellipse_energy_node_min",value(coarse.ellipse_time[1]));
    for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j)report("primary_covariance_"+std::to_string(i)+std::to_string(j),coarse.covariance[i][j]);
    checks.check("reference initial fit fluid matches retained benchmark",std::abs(reference.fit[0]-35.0831075),1e-7);
    checks.check("reference initial fit wall matches retained benchmark",std::abs(reference.fit[1]-14.70559968),1e-7);
    // Check the full observation-cost adjoint at a nonoptimal two-control point.
    Vec<2> control{34.,17.},direction{.3,-.8};double eps=1e-4;
    auto e=coarse_sys.evaluate(control),ep=coarse_sys.evaluate(control+eps*direction),em=coarse_sys.evaluate(control-eps*direction);
    checks.check("two-control adjoint objective finite difference",std::abs((ep.cost-em.cost)/(2*eps)-dot(e.gradient,direction)),2e-6);
    checks.check("two-control incremental adjoint finite difference",norm((.5/eps)*(ep.gradient-em.gradient)-coarse_sys.hessian_vector(control,direction)),2e-6);
    checks.check("nodal thermal cost matches library cost API",std::abs(e.cost-coarse_sys.cost(control,e.trajectory)),1e-10);

    std::ofstream curves(output/"curves.csv");curves<<std::setprecision(17)<<"elapsed,fluid_C,wall_C,truth_fluid_C,truth_wall_C,energy_deviation_kJ,exact_ellipse_fluid_K,exact_ellipse_energy_kJ,discrete_ellipse_fluid_K,discrete_ellipse_energy_kJ\n";
    for(std::size_t n=0;n<coarse_grid.time.size();++n){double tau=coarse_grid.time[n]-thermal::t0;auto truth=thermal::exact(tau,thermal::truth);
        curves<<tau<<','<<coarse.state[n][0]<<','<<coarse.state[n][1]<<','<<truth[0]<<','<<truth[1]<<','<<dot(thermal::goals[1],coarse.state[n]-thermal::equilibrium)<<','<<thermal::ellipse_bound(tau,reference.covariance,thermal::goals[0])<<','<<thermal::ellipse_bound(tau,reference.covariance,thermal::goals[1])<<','<<coarse.ellipse[n][0]<<','<<coarse.ellipse[n][1]<<'\n';}
    std::ofstream envelopes(output/"continuous_envelopes.csv");envelopes<<std::setprecision(17)<<"case,start,end,fluid_bound_K,energy_bound_kJ\n";
    std::ofstream decisions(output/"qualification_decisions.csv");decisions<<std::setprecision(17)<<"case,criterion,fluid_status,fluid_elapsed_min,energy_status,energy_elapsed_min,all_goals_status,all_goals_elapsed_min,blocking_goals\n";
    std::ofstream trace_csv(output/"optimizer_trace.csv");trace_csv<<std::setprecision(17)<<"case,iteration,initial_fluid_C,initial_wall_C,cost,gradient_norm,alpha\n";
    std::ofstream polygons(output/"confidence_polygons.csv");polygons<<std::setprecision(17)<<"case,fluid_initial_C,wall_initial_C\n";
    std::cout<<"\nTWO-TEMPERATURE PHYSICAL QUALIFICATION DETERMINATION\n"
             <<"Goals: fluid error <= 0.2 K; stored-energy error <= 60 kJ.\n"
             <<"95% simultaneous independent-Gaussian noise box; known initial domain [0,60] C squared.\n"
             <<"Numerical error included; entire intervals covered; elapsed horizon [0,10] min.\n";
    for(std::size_t k=0;k<3;++k){const auto& result=k==0?coarse:(k==1?fine:added);const auto& grid=k==0?coarse_grid:fine_grid;const std::string name=k==0?"primary":(k==1?"refined":"added_wall_sensor");
        for(std::size_t j=0;j<result.optimizer_trace.size();++j){const auto& row=result.optimizer_trace[j];trace_csv<<name<<','<<j<<','<<row.control[0]<<','<<row.control[1]<<','<<row.cost<<','<<row.gradient_norm<<','<<row.step_length<<'\n';}
        for(auto vertex:result.polygon)polygons<<name<<','<<vertex[0]<<','<<vertex[1]<<'\n';
        for(std::size_t n=0;n<grid.steps();++n)envelopes<<name<<','<<grid.time[n]-thermal::t0<<','<<grid.time[n+1]-thermal::t0<<','<<result.continuous[n][0]<<','<<result.continuous[n][1]<<'\n';
        auto joint=all_goals_qualification_time(result.physical_time);std::string blockers;
        std::cout<<'\n'<<name<<'\n';decisions<<name<<",continuous_joint95_box,";
        for(std::size_t j=0;j<2;++j){auto t=result.physical_time[j];decisions<<(t?"qualified":"not_qualified")<<',';if(t)decisions<<*t;decisions<<',';
            std::cout<<"  "<<(j==0?"Fluid temperature":"Stored energy")<<": ";
            if(t)std::cout<<"QUALIFIED from "<<*t<<" min through the horizon.\n";else{std::cout<<"NOT QUALIFIED.\n";blockers+=(blockers.empty()?"":"+")+std::string(j==0?"fluid":"energy");}}
        decisions<<(joint?"qualified":"not_qualified")<<',';if(joint)decisions<<*joint;decisions<<','<<blockers<<'\n';
        std::cout<<"  ALL REQUESTED GOALS: ";if(joint)std::cout<<"QUALIFIED from "<<*joint<<" min.\n";else std::cout<<"NOT QUALIFIED; blocking goal: "<<blockers<<".\n";
    }
    // Repeat complete two-control inverse solves on 16 nested meshes.
    std::ofstream convergence(output/"convergence.csv");convergence<<std::setprecision(17)<<"theta,subdivisions,h,state_error,control_error,tangent_error,goal_error,covariance_error,running_cost_error,running_gradient_error\n";
    const auto running_Q=thermal::analytical_running_hessian();auto running_gradient=mv(running_Q,thermal::truth-thermal::equilibrium);double running_value=.5*dot(thermal::truth-thermal::equilibrium,running_gradient);
    for(double theta:{.5,.65}){Vec<7> previous{},last{};
        for(std::size_t subdivisions:{1U,2U,4U,8U,16U,32U,64U,128U}){auto grid=thermal::mesh(subdivisions,theta);auto sys=thermal::system(grid);auto fit=thermal::fit(sys);auto trajectory=sys.forward(thermal::truth);
            auto s0=sys.tangent(thermal::truth,trajectory,{1.,0.}),s1=sys.tangent(thermal::truth,trajectory,{0.,1.});double state_error=0.,tangent_error=0.;
            for(std::size_t n=0;n<grid.time.size();++n){auto E=thermal::propagator(grid.time[n]-thermal::t0);state_error=std::max(state_error,norm(trajectory.state[n]-thermal::exact(grid.time[n]-thermal::t0,thermal::truth)));tangent_error=std::max({tangent_error,norm(s0[n]-Vec<2>{E[0][0],E[1][0]}),norm(s1[n]-Vec<2>{E[0][1],E[1][1]})});}
            Mat<2,2> information{};for(std::size_t j=0;j<2;++j){Vec<2> unit{};unit[j]=1.;auto col=sys.hessian_vector(fit.control,unit);for(std::size_t i=0;i<2;++i)information[i][j]=col[i];}
            auto goal=sys.goal_gradient(thermal::truth,trajectory,grid.steps(),thermal::goals[1]);auto exact_goal=mv(transpose(thermal::propagator(10.)),thermal::goals[1]);
            System<2,2> running_system(thermal::running_problem(grid),grid,thermal::solver_options);auto running=running_system.evaluate(thermal::truth);
            if(subdivisions==16){auto weights=grid.weights();double sum=0.;for(std::size_t n=0;n<weights.size();++n){double r=running.trajectory.state[n][0]-thermal::equilibrium[0];sum+=.5*weights[n]*r*r;}
                checks.check("thermal theta running quadrature",std::abs(sum-running.cost),1e-12);
                auto hp=running_system.hessian_vector(thermal::truth,direction);double re=1e-4;auto pp=running_system.evaluate(thermal::truth+re*direction),pm=running_system.evaluate(thermal::truth-re*direction);
                checks.check("thermal running adjoint gradient FD",std::abs((pp.cost-pm.cost)/(2*re)-dot(running.gradient,direction)),2e-7);
                checks.check("thermal running Hessian FD",norm((.5/re)*(pp.gradient-pm.gradient)-hp),2e-7);}
            previous=last;last={state_error,norm(fit.control-reference.fit),tangent_error,norm(goal-exact_goal),matrix_error(thermal::inverse(information),reference.covariance),std::abs(running.cost-running_value),norm(running.gradient-running_gradient)};
            double h=0.;for(std::size_t n=0;n<grid.steps();++n)h=std::max(h,grid.time[n+1]-grid.time[n]);
            convergence<<theta<<','<<subdivisions<<','<<h;for(double err:last)convergence<<','<<err;convergence<<'\n';
            checks.check("refined thermal inverse stationarity",norm(sys.evaluate(fit.control).gradient),1e-8);
        }
        for(std::size_t j=0;j<7;++j){double p=std::log2(previous[j]/last[j]);report("order_"+std::to_string(theta)+"_component_"+std::to_string(j),p);checks.check("thermal temporal order",std::abs(p-(theta==.5?2.:1.)),.15);}
    }
    std::ofstream observation_csv(output/"observations.csv");observation_csv<<std::setprecision(17)<<"elapsed,truth_fluid_C,observation_fluid_C\n";
    for(std::size_t j=0;j<thermal::elapsed.size();++j)observation_csv<<thermal::elapsed[j]<<','<<thermal::exact(thermal::elapsed[j],thermal::truth)[0]<<','<<thermal::observations[j]<<'\n';
    report("thermal_checks_passed",static_cast<double>(checks.count));
    std::cout<<"\nThermal numerical tests: PASS. Test success and engineering qualification are separate.\n"
             <<"Ellipse times above are reference/model diagnostics; final verdicts use whole-step physical confidence envelopes.\n"
             <<"All inference and qualification decisions ran in C++20. Results: "<<output<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
