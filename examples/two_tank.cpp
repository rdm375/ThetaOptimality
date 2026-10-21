// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#include "two_tank_model.hpp"
#include "two_tank_qualification_model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace two_tank;
namespace {
inline constexpr Vec<2> scale{1.,1.},probe{.25,.85};
double normalized(const Vec<2>& e){return std::max(std::abs(e[0])/scale[0],std::abs(e[1])/scale[1]);}
double matrix_error(const Mat<2,2>& x,const Mat<2,2>& y){return std::max(norm(x[0]-y[0]),norm(x[1]-y[1]));}
struct Checks {
    std::ofstream csv;std::size_t count=0;
    explicit Checks(const std::filesystem::path& path):csv(path){csv<<"name,error,tolerance,passed\n"<<std::setprecision(17);}
    void operator()(const std::string& name,double error,double tol){++count;bool pass=std::isfinite(error)&&error<=tol;csv<<name<<','<<error<<','<<tol<<','<<pass<<'\n';if(!pass)throw std::runtime_error(name+" failed: "+std::to_string(error));}
};
}
int main(int argc,char** argv){try{
    std::filesystem::path out=argc>1?argv[1]:"results/two_tank";std::filesystem::create_directories(out);Checks check(out/"checks.csv");
    auto reference=analytical_fit();auto H=information(),Q=gramian(horizon);auto covariance=inverse(H);
    std::ofstream summary(out/"summary.csv");summary<<std::setprecision(17)<<"quantity,value\n";
    auto report=[&](const std::string& name,double value){summary<<name<<','<<value<<'\n';std::cout<<std::setprecision(12)<<name<<": "<<value<<'\n';};
    report("analytical_initial_c1_mg_L",reference[0]);report("analytical_initial_c2_mg_L",reference[1]);
    report("truth_c1_mg_L",truth[0]);report("truth_c2_mg_L",truth[1]);
    report("local_initial_c1_sd_mg_L",std::sqrt(covariance[0][0]));report("local_initial_c2_sd_mg_L",std::sqrt(covariance[1][1]));
    check("analytical noiseless recovery",normalized(analytical_fit(false,false)-truth),2e-14);
    auto p03=propagator(.3),p10=propagator(1.);Mat<2,2> product{};for(int j=0;j<2;++j){Vec<2> col{p03[0][j],p03[1][j]};auto result=mv(p10,col);for(int i=0;i<2;++i)product[i][j]=result[i];}check("analytical semigroup",matrix_error(propagator(1.3),product),1e-15);
    std::ofstream observations_file(out/"observations.csv");observations_file<<std::setprecision(17)<<"elapsed_min,observed_c1_mg_L,true_c1_mg_L,sigma_mg_L\n";
    for(std::size_t j=0;j<elapsed.size();++j)observations_file<<elapsed[j]<<','<<observation(j)<<','<<exact(elapsed[j],truth)[0]<<','<<sigma<<'\n';
    constexpr std::size_t components=13;
    const std::array<std::string,components> names{"state_c1","state_c2","tangent","nodal_adjoint","control_c1","control_c2","running_cost","running_gradient","running_costate","running_hessian","nodal_cost","nodal_gradient","nodal_hessian"};
    std::ofstream controls(out/"control_convergence.csv");controls<<std::setprecision(17)<<"theta,steps,h_min,discrete_c10_mg_L,discrete_c20_mg_L,analytical_c10_mg_L,analytical_c20_mg_L\n";
    std::ofstream convergence(out/"convergence.csv"),orders(out/"orders.csv");convergence<<std::setprecision(17)<<"theta,steps,h_min";for(auto name:names)convergence<<','<<name;convergence<<'\n';orders<<std::setprecision(17)<<"theta,component,order\n";
    for(double theta:{.5,.65}){
        std::array<double,components> previous{},errors{};
        for(std::size_t steps:{40u,80u,160u,320u,640u,1280u}){
            auto sys=system(steps,theta);auto ev=sys.evaluate(probe);auto optimum=fit(sys);auto fitted=sys.evaluate(optimum.control);
            controls<<theta<<','<<steps<<','<<horizon/steps<<','<<optimum.control[0]<<','<<optimum.control[1]<<','<<reference[0]<<','<<reference[1]<<'\n';
            if(steps==1280){report("finest_c10_theta_"+std::to_string(theta),optimum.control[0]);report("finest_c20_theta_"+std::to_string(theta),optimum.control[1]);}
            auto run=system(steps,theta,true);auto re=run.evaluate(probe);auto costate=run.cost_to_go_adjoint(probe,re);
            auto tangent=sys.tangent(probe,ev.trajectory,scale);
            std::array<Vec<2>,2> hcols,runcols;Mat<2,2> Hd{},Qd{};
            for(int j=0;j<2;++j){Vec<2> e{};e[j]=1.;hcols[j]=sys.hessian_vector(probe,e);runcols[j]=run.hessian_vector(probe,e);for(int i=0;i<2;++i){Hd[i][j]=hcols[j][i];Qd[i][j]=runcols[j][i];}}
            check("optimized discrete stationarity",norm(fitted.gradient),1e-8);
            check("one Newton correction for quadratic objective",std::abs(static_cast<double>(optimum.trace.size())-2),0.);
            check("forward residual",ev.trajectory.max_residual,5e-13);
            check("nodal Hessian symmetry",std::abs(Hd[0][1]-Hd[1][0])/std::max(norm(Hd[0]),norm(Hd[1])),5e-14);
            double sq=0.,sv=0.,st=0.,sa=0.,sr=0.,adjscale=0.,runscale=0.;
            double balance=0.;
            for(std::size_t n=0;n<=steps;++n){double tau=sys.grid().time[n]-t0;auto x=exact(tau,probe),s=exact(tau,scale),a=nodal_costate(tau,probe),ra=running_costate(tau,probe);
                sq=std::max(sq,std::abs(ev.trajectory.state[n][0]-x[0])/scale[0]);sv=std::max(sv,std::abs(ev.trajectory.state[n][1]-x[1])/scale[1]);st=std::max(st,normalized(tangent[n]-s));
                sa=std::max(sa,norm(Vec<2>{scale[0]*(ev.nodal_adjoint[n][0]-a[0]),scale[1]*(ev.nodal_adjoint[n][1]-a[1])}));
                sr=std::max(sr,norm(Vec<2>{scale[0]*(costate[n][0]-ra[0]),scale[1]*(costate[n][1]-ra[1])}));
                adjscale=std::max(adjscale,norm(Vec<2>{scale[0]*a[0],scale[1]*a[1]}));runscale=std::max(runscale,norm(Vec<2>{scale[0]*ra[0],scale[1]*ra[1]}));
                if(n){auto old=ev.trajectory.state[n-1],next=ev.trajectory.state[n],mid=(1-theta)*old+theta*next;double h=sys.grid().time[n]-sys.grid().time[n-1];double identity_error=inventory(next)-inventory(old)+h*flow*mid[0];balance=std::max(balance,std::abs(identity_error));}
            }
            // Closed-form 2x2 theta step, independent of the forward Newton solve.
            Vec<2> independent=probe;double step_error=0.,h=horizon/steps;
            double b00=1+h*theta*fast,b01=-h*theta*fast,b10=0.,b11=1+h*theta*slow,det=b00*b11-b01*b10;
            for(std::size_t n=1;n<=steps;++n){auto rhs=independent+h*(1-theta)*mv(A,independent);independent={(b11*rhs[0]-b01*rhs[1])/det,(-b10*rhs[0]+b00*rhs[1])/det};step_error=std::max(step_error,normalized(independent-ev.trajectory.state[n]));}
            check("Newton versus explicit rational theta map",step_error,2e-12);
            check("theta tracer mass balance",balance,2e-14);
            check("terminal running costate zero",norm(costate.back()),1e-15);
            previous=errors;
            double Jr=.5*dot(probe,mv(Q,probe));auto gr=mv(Q,probe),g=exact_gradient(probe);
            errors={sq,sv,st,sa/adjscale,std::abs(optimum.control[0]-reference[0])/scale[0],std::abs(optimum.control[1]-reference[1])/scale[1],std::abs(re.cost-Jr)/Jr,norm(re.gradient-gr)/norm(gr),sr/runscale,matrix_error(Qd,Q)/std::max(norm(Q[0]),norm(Q[1])),std::abs(ev.cost-exact_cost(probe))/exact_cost(probe),norm(ev.gradient-g)/norm(g),matrix_error(Hd,H)/std::max(norm(H[0]),norm(H[1]))};
            convergence<<theta<<','<<steps<<','<<horizon/steps;for(auto error:errors)convergence<<','<<error;convergence<<'\n';
            if(theta==.5&&steps==160){
                double epsilon=1e-4;Vec<2> direction{.03,-.04};auto ep=sys.evaluate(probe+epsilon*direction),em=sys.evaluate(probe-epsilon*direction);
                check("adjoint cost directional difference",std::abs((ep.cost-em.cost)/(2*epsilon)-dot(ev.gradient,direction)),2e-6);
                check("incremental adjoint directional difference",norm((1/(2*epsilon))*(ep.gradient-em.gradient)-sys.hessian_vector(probe,direction)),5e-6);
                double paired=0.;auto ts=sys.tangent(probe,ev.trajectory,direction);for(std::size_t j=0;j<elapsed.size();++j){auto n=node_at(sys.grid(),elapsed[j]);paired+=(ev.trajectory.state[n][0]-observation(j))*ts[n][0]/(sigma*sigma);}check("tangent-adjoint duality",std::abs(paired-dot(ev.gradient,direction)),2e-9);
            }
            if(theta==.5&&steps==160){
                std::ofstream curves(out/"trajectories.csv");curves<<std::setprecision(17)<<"elapsed_min,true_c1_mg_L,true_c2_mg_L,analytical_fit_c1_mg_L,analytical_fit_c2_mg_L,discrete_c1_mg_L,discrete_c2_mg_L,true_inventory_mg,discrete_inventory_mg,discrete_probe_running_costate_c1,exact_probe_running_costate_c1,discrete_probe_running_costate_c2,exact_probe_running_costate_c2\n";
                for(std::size_t n=0;n<=steps;++n){double tau=sys.grid().time[n]-t0;auto xt=exact(tau,truth),xa=exact(tau,reference),xd=fitted.trajectory.state[n],ra=running_costate(tau,probe);curves<<tau<<','<<xt[0]<<','<<xt[1]<<','<<xa[0]<<','<<xa[1]<<','<<xd[0]<<','<<xd[1]<<','<<inventory(xt)<<','<<inventory(xd)<<','<<costate[n][0]<<','<<ra[0]<<','<<costate[n][1]<<','<<ra[1]<<'\n';}
                std::ofstream trace(out/"optimizer_trace.csv");trace<<std::setprecision(17)<<"iteration,cost,gradient_norm,c10_mg_L,c20_mg_L\n";for(std::size_t j=0;j<optimum.trace.size();++j){auto it=optimum.trace[j];trace<<j<<','<<it.cost<<','<<it.gradient_norm<<','<<it.control[0]<<','<<it.control[1]<<'\n';}
                report("primary_discrete_c1_mg_L",optimum.control[0]);report("primary_discrete_c2_mg_L",optimum.control[1]);
            }
        }
        for(std::size_t j=0;j<components;++j){
            if(j==4){check("direct initial observation control remains at roundoff scale",errors[j],2e-10);orders<<theta<<','<<names[j]<<",nan\n";continue;}
            double order=std::log2(previous[j]/errors[j]);orders<<theta<<','<<names[j]<<','<<order<<'\n';check("temporal order "+names[j],std::abs(order-(theta==.5?2.:1.)),theta==.5?.07:.22);std::cout<<"order "<<theta<<" / "<<names[j]<<": "<<order<<'\n';}

    }
    // Gain benchmark: analytical scalar formula vs library singular gain and
    // a tangent/adjoint power iteration through the actual discrete model.
    double peak=peak_time(),maximum_gain=analytical_gain(peak);report("analytical_peak_time_min",peak);report("analytical_peak_gain",maximum_gain);
    report("euclidean_logarithmic_norm",.5*(-fast-slow+std::hypot(fast-slow,fast)));
    check("analytical gain stationary at peak",std::abs(gain_derivative(peak)),2e-13);
    check("transient gain exceeds unity",maximum_gain>1.2?0.:1.,0.);
    auto gain_sys=system(1280,.5);auto gain_tr=gain_sys.forward(truth);auto maps=control_to_state_jacobians(gain_sys,truth,gain_tr);
    std::ofstream gain_csv(out/"amplification.csv");gain_csv<<std::setprecision(17)<<"elapsed_min,analytical_gain,discrete_gain,direction_c1,direction_c2\n";
    double gain_error=0.;std::size_t peak_node=0;double discrete_peak=0.;
    for(std::size_t n=0;n<maps.size();++n){double tau=gain_sys.grid().time[n]-t0;auto g=worst_amplification_2d(maps[n]);gain_csv<<tau<<','<<analytical_gain(tau)<<','<<g.gain<<','<<g.scaled_input[0]<<','<<g.scaled_input[1]<<'\n';gain_error=std::max(gain_error,std::abs(worst_amplification_2d(propagator(tau)).gain-analytical_gain(tau)));if(g.gain>discrete_peak){discrete_peak=g.gain;peak_node=n;}}
    check("closed analytical gain versus library SVD",gain_error,5e-14);
    report("discrete_peak_time_min",gain_sys.grid().time[peak_node]-t0);report("discrete_peak_gain",discrete_peak);
    check("discrete peak time versus analytical",std::abs(gain_sys.grid().time[peak_node]-t0-peak),horizon/1280);
    check("discrete peak gain versus analytical",std::abs(discrete_peak-maximum_gain),3e-5);
    Vec<2> direction{1.,.3};direction=(1/std::hypot(direction[0],direction[1]))*direction;
    for(int iteration=0;iteration<50;++iteration){auto s=gain_sys.tangent(truth,gain_tr,direction);auto pullback=gain_sys.goal_gradient(truth,gain_tr,peak_node,s[peak_node]);direction=(1/std::hypot(pullback[0],pullback[1]))*pullback;}
    auto optimal=worst_amplification_2d(maps[peak_node]);check("tangent-adjoint worst direction",norm(direction-optimal.scaled_input),2e-11);
    auto sensitivity=gain_sys.tangent(truth,gain_tr,direction);check("worst direction attains discrete gain",std::abs(std::hypot(sensitivity[peak_node][0],sensitivity[peak_node][1])-optimal.gain),2e-12);
    report("worst_direction_c1",direction[0]);report("worst_direction_c2",direction[1]);
    auto gain_pullback=gain_sys.goal_gradient(truth,gain_tr,peak_node,sensitivity[peak_node]);check("gain tangent-adjoint pairing",std::abs(dot(direction,gain_pullback)-optimal.gain*optimal.gain),2e-12);
    std::ofstream worst_path(out/"worst_direction.csv");worst_path<<std::setprecision(17)<<"elapsed_min,perturbation_c1,perturbation_c2,euclidean_gain,weighted_mass_error_bound\n";
    double mass_increase=0.;for(std::size_t n=0;n<sensitivity.size();++n){auto x=sensitivity[n];if(n)mass_increase=std::max(mass_increase,inventory(x)-inventory(sensitivity[n-1]));worst_path<<gain_sys.grid().time[n]-t0<<','<<x[0]<<','<<x[1]<<','<<std::hypot(x[0],x[1])<<','<<volume1*std::abs(x[0])+volume2*std::abs(x[1])<<'\n';}
    check("positive worst perturbation loses total mass while concentration norm grows",mass_increase,1e-13);
    namespace tq=two_tank::qualification;
    std::ofstream estimates(out/"estimates.csv"),extra_data(out/"experiment_observations.csv");
    estimates<<std::setprecision(17)<<"case,discrete_c10,discrete_c20,analytical_c10,analytical_c20\n";extra_data<<std::setprecision(17)<<"case,elapsed_min,component,value_mg_L,scale_mg_L\n";
    std::ofstream decisions(out/"qualification_decisions.csv"),envelopes(out/"qualification_envelopes.csv"),polygons(out/"qualification_polygons.csv");
    decisions<<std::setprecision(17)<<"case,goal,units,tolerance,first_below_min,persistent_min,analytical_min,qualified,joint_qualified,joint_min,blocking_goals\n";
    envelopes<<std::setprecision(17)<<"case,elapsed_min,interval_end_min,goal,tolerance,node_bound,whole_bound,analytical_bound,local95_indicator,numerical_allowance\n";
    polygons<<std::setprecision(17)<<"case,set,vertex,c10,c20\n";
    struct Specification {const char* name;std::size_t steps;bool added,gaussian;};
    std::cout<<"\nTWO-TANK ENGINEERING QUALIFICATION\nKnown prior: c10 in [0,0.8], c20 in [0,1.6] mg/L.\nTolerances: downstream 0.005 mg/L, upstream 0.012 mg/L, inventory 0.04 mg.\nWhole intervals through elapsed 10 min; all-data retrospective decisions.\n";
    for(auto spec:std::array<Specification,4>{{{"bounded_primary",80,false,false},{"bounded_refined",1280,false,false},{"added_upstream",1280,true,false},{"joint95_refined",1280,false,true}}}){
        auto qs=system(spec.steps,.5,false,spec.added);auto result=tq::experiment(qs,spec.added,spec.gaussian,[&](const std::string& name,double error,double tol){check(std::string(spec.name)+": "+name,error,tol);});
        bool expected_joint=std::string(spec.name)!="bounded_primary";check(std::string(spec.name)+" expected engineering verdict",static_cast<bool>(result.decision.joint)==expected_joint?0.:1.,0.);
        if(std::string(spec.name)=="bounded_refined"){bool rebound=false;for(std::size_t n=0;n<spec.steps;++n)if(result.whole[n][0]>tq::tolerance[0])rebound=true;check("initially acceptable downstream goal subsequently amplifies beyond tolerance",result.decision.first_below[0]==0.&&result.decision.persistent[0].value_or(0)>0.&&rebound?0.:1.,0.);}
        if(spec.added)check("actual added measurement restores whole-horizon qualification",result.decision.joint.value_or(-1),0.);
        auto ref=analytical_fit(spec.added);estimates<<spec.name<<','<<result.center[0]<<','<<result.center[1]<<','<<ref[0]<<','<<ref[1]<<'\n';for(auto o:observations(spec.added))extra_data<<spec.name<<','<<o.tau<<','<<o.component<<','<<o.value<<','<<o.scale<<'\n';
        auto write_time=[](std::ostream& o,const std::optional<double>& t){if(t)o<<*t;};std::string blockers;
        for(std::size_t k=0;k<3;++k)if(!result.decision.persistent[k]){if(!blockers.empty())blockers+='+';blockers+=tq::names[k];}
        if(blockers.empty())blockers="none";
        std::cout<<spec.name<<": ";if(result.decision.joint)std::cout<<"ALL GOALS QUALIFIED from "<<*result.decision.joint<<" min.\n";else std::cout<<"NOT QUALIFIED; blockers: "<<blockers<<'\n';
        for(std::size_t k=0;k<3;++k){decisions<<spec.name<<','<<tq::names[k]<<','<<tq::units[k]<<','<<tq::tolerance[k]<<',';write_time(decisions,result.decision.first_below[k]);decisions<<',';write_time(decisions,result.decision.persistent[k]);decisions<<',';write_time(decisions,result.analytical_time[k]);decisions<<','<<static_cast<bool>(result.decision.persistent[k])<<','<<static_cast<bool>(result.decision.joint)<<',';write_time(decisions,result.decision.joint);decisions<<','<<blockers<<'\n';}
        for(std::size_t n=0;n<qs.grid().time.size();++n)for(std::size_t k=0;k<3;++k){double tau=qs.grid().time[n]-t0;envelopes<<spec.name<<','<<tau<<','<<(n<spec.steps?qs.grid().time[n+1]-t0:tau)<<','<<tq::names[k]<<','<<tq::tolerance[k]<<','<<result.node[n][k]<<',';if(n<spec.steps)envelopes<<result.whole[n][k];envelopes<<','<<tq::exact_spread(tau,result.reference,analytical_fit(spec.added),tq::goals[k])<<','<<result.local95[n][k]<<','<<result.error[n]<<'\n';}
        for(bool exact_set:{false,true}){const auto& p=exact_set?result.reference:result.polygon;for(std::size_t j=0;j<p.size();++j)polygons<<spec.name<<','<<(exact_set?"analytical":"numerical_outer")<<','<<j<<','<<p[j][0]<<','<<p[j][1]<<'\n';}
    }
    report("checks_passed",static_cast<double>(check.count));std::cout<<"Two-tank validation PASS; engineering decisions are reported separately.\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
