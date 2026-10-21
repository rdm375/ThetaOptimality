// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Richard Myers
#include "peng_robinson_model.hpp"
#include "peng_robinson_bounds.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace gas_inverse;
int main(int argc,char** argv){try{
    std::filesystem::path out=argc>1?argv[1]:"results/peng_robinson";std::filesystem::create_directories(out);
    std::size_t count=0;std::ofstream checks(out/"checks.csv");checks<<"check,passed\n";
    auto check=[&](const std::string& name,bool ok){checks<<name<<','<<ok<<'\n';if(!ok)throw std::runtime_error("Failed: "+name);++count;};
    auto data=observations();auto sys=system(data);auto optimum=fit(sys);auto x=optimum.control;auto ev=sys.evaluate(x);
    std::ofstream summary(out/"summary.csv");summary<<std::setprecision(17)<<"quantity,value\n";
    auto report=[&](std::string name,double v){summary<<name<<','<<v<<'\n';std::cout<<std::setprecision(12)<<name<<": "<<v<<'\n';};
    report("fit_ethane",x[0]);report("fit_nitrogen",x[1]);report("fit_methane",1-x[0]-x[1]);report("cost",ev.cost);report("gradient_norm",norm(ev.gradient));
    for(auto start:std::array<Vec<2>,3>{{{.08,.04},{.12,.06},{.09,.08}}})check("multistart_same_stationary_point",norm(fit(sys,start).control-x)<1e-8);
    auto exact=fit(system(observations(false)));check("noiseless_truth_recovery",norm(exact.control-truth)<1e-9);
    Vec<2> v{.3,-.2};double eps=2e-5;auto p=x+eps*v,m=x-eps*v;
    for(const auto& o:data){auto r=density(o.state,x),rp=density(o.state,p),rm=density(o.state,m);
        check("unique_root_residual",std::abs(r.residual)<1e-12&&r.Z>0&&std::abs(r.root_slope)>1e-5);
        check("implicit_density_direction",std::abs((rp.value-rm.value)/(2*eps)-dot(r.gradient,v))<1e-6);
        check("implicit_density_hessian",norm((1/(2*eps))*(rp.gradient-rm.gradient)-mv(r.hessian,v))<1e-5);
    }
    auto trajectory=sys.forward(x);auto tangent=sys.tangent(x,trajectory,v);
    for(auto s:tangent)check("constant_composition_tangent",norm(s-v)<1e-14);
    auto ep=sys.evaluate(p),em=sys.evaluate(m);auto hv=sys.hessian_vector(x,v);
    check("adjoint_vs_cost_difference",std::abs((ep.cost-em.cost)/(2*eps)-dot(ev.gradient,v))<1e-3);
    check("incremental_adjoint_vs_gradient_difference",norm((1/(2*eps))*(ep.gradient-em.gradient)-hv)<.02);
    // Differentiate away from stationarity so cancellation cannot mask a wrong gradient.
    Vec<2> probe{.08,.07};double de=1e-6;auto probe_eval=sys.evaluate(probe);
    double fd=(sys.evaluate(probe+de*v).cost-sys.evaluate(probe-de*v).cost)/(2*de);
    check("off_solution_adjoint_direction",std::abs(fd-dot(probe_eval.gradient,v))<2e-5);
    double he=1e-5;auto dh=(1/(2*he))*(sys.evaluate(probe+he*v).gradient-sys.evaluate(probe-he*v).gradient);
    auto exact_hv=sys.hessian_vector(probe,v);
    std::ofstream derivative(out/"derivative_check.csv");derivative<<std::setprecision(17)<<"perturbation,hessian_direction_error\n";
    for(double step:{1e-3,1e-4,1e-5,1e-6,1e-7}){auto diff=(1/(2*step))*(sys.evaluate(probe+step*v).gradient-sys.evaluate(probe-step*v).gradient);derivative<<step<<','<<norm(diff-exact_hv)<<'\n';}
    check("off_solution_incremental_adjoint",norm(dh-exact_hv)<2e-4);
    auto probe_traj=sys.forward(probe);auto ts=sys.tangent(probe,probe_traj,v);
    double tangent_cost=0.;for(std::size_t j=0;j<data.size();++j){auto r=density(data[j].state,probe_traj.state[j]);tangent_cost+=(r.value-data[j].value)*dot(r.gradient,ts[j])/(data[j].sigma*data[j].sigma);}
    check("tangent_adjoint_duality",std::abs(tangent_cost-dot(probe_eval.gradient,v))<1e-9);
    for(const auto& o:data){auto r=thermogpu::evaluate_mixture(mixture(x),o.state);double RT=thermogpu::gas_constant*o.state.temperature_K,volume=r.molar_mass_kg_per_mol/r.density_kg_per_m3;
        double den=volume*(volume+r.b_m)+r.b_m*(volume-r.b_m);
        double pressure=RT/(volume-r.b_m)-r.a_m/den;
        check("density_pressure_consistency",std::abs(pressure-o.state.pressure_Pa)/o.state.pressure_Pa<1e-12);
        double dpdv=-RT/((volume-r.b_m)*(volume-r.b_m))+r.a_m*(2*volume+2*r.b_m)/(den*den);
        check("fitted_gas_mechanically_stable",dpdv<0);
    }
    auto bad_options=OptimizeOptions{};bad_options.cost_roundoff_slack=-1.;bool rejected=false;
    try{minimize(sys,x,bad_options,admissible);}catch(const std::invalid_argument&){rejected=true;}
    check("negative_oracle_slack_rejected",rejected);
    Vec<2> direct{};for(auto o:data){auto r=density(o.state,x);direct=direct+((r.value-o.value)/(o.sigma*o.sigma))*r.gradient;}
    check("adjoint_vs_independent_chain_rule",norm(ev.gradient-direct)<1e-9);
    check("theta_static_invariance",norm(fit(system(data,.65)).control-x)<1e-12);
    auto I=information(sys,data,x);double det=I[0][0]*I[1][1]-I[0][1]*I[1][0];check("four_observation_information_full_rank",det>0);
    double tr=I[0][0]+I[1][1],gap=std::hypot(I[0][0]-I[1][1],2*I[0][1]);double large=.5*(tr+gap),small=det/large;report("information_condition_number",large/small);
    Mat<2,2> covariance{};for(int j=0;j<2;++j){Vec<2> unit{};unit[j]=1.;auto col=LU<2>(I).solve(unit);for(int i=0;i<2;++i)covariance[i][j]=col[i];}
    report("local_ethane_standard_deviation",std::sqrt(covariance[0][0]));report("local_nitrogen_standard_deviation",std::sqrt(covariance[1][1]));
    auto single=data;single.resize(1);auto g=density(data[0].state,x).gradient;Vec<2> null{g[1],-g[0]};check("one_density_has_unobservable_direction",std::abs(dot(g,null))<1e-10);
    std::ofstream obs(out/"observations.csv");obs<<std::setprecision(17)<<"index,temperature_K,pressure_Pa,observed_density_kg_m3,sigma_kg_m3,bound_kg_m3,fitted_density_kg_m3,true_density_kg_m3,Z\n";
    for(std::size_t j=0;j<data.size();++j){auto o=data[j];auto r=density(o.state,x);obs<<j<<','<<o.state.temperature_K<<','<<o.state.pressure_Pa<<','<<o.value<<','<<o.sigma<<','<<o.bound<<','<<r.value<<','<<density(o.state,truth).value<<','<<r.Z<<'\n';}
    std::ofstream trace(out/"optimizer_trace.csv");trace<<std::setprecision(17)<<"iteration,cost,gradient_norm,ethane,nitrogen,step_length\n";
    for(std::size_t k=0;k<optimum.trace.size();++k){auto it=optimum.trace[k];trace<<k<<','<<it.cost<<','<<it.gradient_norm<<','<<it.control[0]<<','<<it.control[1]<<','<<it.step_length<<'\n';}
    // Exact reduced Hessian need not equal Gauss-Newton information with noise.
    check("hessian_symmetry",std::abs(sys.hessian_vector(x,{1.,0.})[1]-sys.hessian_vector(x,{0.,1.})[0])<1e-7);
    std::ofstream qualification(out/"qualification.csv"),cells(out/"outer_cells.csv");qualification<<std::setprecision(17)<<"case,goal,lower,upper,estimate,tolerance,worst_error,qualified\n";cells<<std::setprecision(17)<<"case,ethane_lower,ethane_upper,nitrogen_lower,nitrogen_upper\n";
    std::ofstream decisions(out/"qualification_decisions.csv"),experiment_obs(out/"experiment_observations.csv");
    decisions<<"case,joint_qualified,blocking_goals\n";
    experiment_obs<<std::setprecision(17)<<"case,index,temperature_K,pressure_Pa,density_kg_m3,bound_kg_m3\n";
    std::array<std::string,3> names{"methane","ethane","nitrogen"};
    double lo=0.,hi=.3;for(int k=0;k<60;++k){double mid=.5*(lo+hi);if(density(single[0].state,{mid,.08}).value<single[0].value)lo=mid;else hi=mid;}
    Vec<2> single_representative{.5*(lo+hi),.08};
    check("single_representative_fits",std::abs(density(single[0].state,single_representative).value-single[0].value)<1e-10);
    for(int scenario=0;scenario<3;++scenario){bool all=scenario>0;std::string name=scenario==0?"one_density":scenario==1?"four_densities":"four_precise_densities";auto selected=scenario==0?single:scenario==1?data:observations(true,.0002);auto outer=outer_set(selected,scenario==2?.000125:.0005);check(name+"_outer_nonempty",!outer.cells.empty());check(name+"_truth_retained",possible({truth,truth},selected));
        auto range=ranges(outer);
        for(int a=0;a<5;++a)for(int b=0;b<5;++b){Vec<2> candidate{.09+.005*a,.04+.005*b};bool feasible=true;for(auto o:selected)feasible=feasible&&std::abs(density(o.state,candidate).value-o.value)<=o.bound;if(feasible)check(name+"_feasible_probe_retained",possible({candidate,candidate},selected));}
        bool contains_truth=false;for(auto c:outer.cells){if(c.lo[0]<=truth[0]&&c.hi[0]>=truth[0]&&c.lo[1]<=truth[1]&&c.hi[1]>=truth[1])contains_truth=true;cells<<name<<','<<c.lo[0]<<','<<c.hi[0]<<','<<c.lo[1]<<','<<c.hi[1]<<'\n';}check(name+"_truth_in_union",contains_truth);
        for(std::size_t j=0;j<selected.size();++j){auto o=selected[j];experiment_obs<<name<<','<<j<<','<<o.state.temperature_K<<','<<o.state.pressure_Pa<<','<<o.value<<','<<o.bound<<'\n';}
        auto estimate=all?fractions(scenario==1?x:fit(system(selected)).control):fractions(single_representative);
        bool joint=true;std::string blockers;for(int j=0;j<3;++j){double error=up(std::max(std::abs(range[j][0]-estimate[j]),std::abs(range[j][1]-estimate[j])));bool qualified=error<=.01;joint=joint&&qualified;if(!qualified){if(!blockers.empty())blockers+="+";blockers+=names[j];}std::cout<<name<<" / "<<names[j]<<": "<<(qualified?"QUALIFIED":"NOT QUALIFIED")<<"; outer ["<<range[j][0]<<", "<<range[j][1]<<"], worst error "<<error<<'\n';qualification<<name<<','<<names[j]<<','<<range[j][0]<<','<<range[j][1]<<','<<estimate[j]<<",0.01,"<<error<<','<<qualified<<'\n';}
        decisions<<name<<','<<joint<<','<<(blockers.empty()?"none":blockers)<<'\n';
        std::cout<<name<<": "<<(joint?"QUALIFIED":"NOT QUALIFIED")<<" for all three mole fractions (absolute tolerance 0.01); retained "<<outer.cells.size()<<" cells, visited "<<outer.visited<<'\n';
        check(name+"_expected_verdict",joint==(scenario==2));
        report(name+"_cells",static_cast<double>(outer.cells.size()));
    }
    // Objective surface for plotting only; all inverse computations remain C++.
    std::ofstream surface(out/"objective_surface.csv");surface<<std::setprecision(17)<<"ethane,nitrogen,cost\n";
    for(int a=0;a<=100;++a)for(int b=0;b<=100;++b){Vec<2> q{.02+.0018*a,.001+.00149*b};surface<<q[0]<<','<<q[1]<<','<<sys.evaluate(q).cost<<'\n';}
    report("checks_passed",static_cast<double>(count));std::cout<<"Numerical validation passed; qualification verdicts are separate.\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
