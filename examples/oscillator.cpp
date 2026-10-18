#include "oscillator_model.hpp"
#include "oscillator_qualification_model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace oscillator;
namespace {
inline constexpr Vec<2> scale{.02,.12},probe{.013,.085};
double normalized(const Vec<2>& e){return std::max(std::abs(e[0])/scale[0],std::abs(e[1])/scale[1]);}
double matrix_error(const Mat<2,2>& x,const Mat<2,2>& y){return std::max(norm(x[0]-y[0]),norm(x[1]-y[1]));}
struct Checks {
    std::ofstream csv;std::size_t count=0;
    explicit Checks(const std::filesystem::path& path):csv(path){csv<<"name,error,tolerance,passed\n"<<std::setprecision(17);}
    void operator()(const std::string& name,double error,double tol){++count;bool pass=std::isfinite(error)&&error<=tol;csv<<name<<','<<error<<','<<tol<<','<<pass<<'\n';if(!pass)throw std::runtime_error(name+" failed: "+std::to_string(error));}
};
}
int main(int argc,char** argv){try{
    std::filesystem::path out=argc>1?argv[1]:"results/oscillator";std::filesystem::create_directories(out);Checks check(out/"checks.csv");
    auto reference=analytical_fit();auto H=information(),Q=gramian(horizon);auto covariance=inverse(H);
    std::ofstream summary(out/"summary.csv");summary<<std::setprecision(17)<<"quantity,value\n";
    auto report=[&](const std::string& name,double value){summary<<name<<','<<value<<'\n';std::cout<<std::setprecision(12)<<name<<": "<<value<<'\n';};
    report("analytical_initial_displacement_m",reference[0]);report("analytical_initial_velocity_m_s",reference[1]);
    report("truth_displacement_m",truth[0]);report("truth_velocity_m_s",truth[1]);
    report("local_initial_displacement_sd_m",std::sqrt(covariance[0][0]));report("local_initial_velocity_sd_m_s",std::sqrt(covariance[1][1]));
    check("analytical noiseless recovery",normalized(analytical_fit(false)-truth),2e-14);
    auto p03=propagator(.3),p10=propagator(1.);Mat<2,2> product{};for(int j=0;j<2;++j){Vec<2> col{p03[0][j],p03[1][j]};auto result=mv(p10,col);for(int i=0;i<2;++i)product[i][j]=result[i];}check("analytical semigroup",matrix_error(propagator(1.3),product),1e-15);
    // Deliberately poor design: every observation at an integer half-period.
    double blind=0.;Mat<2,2> poor{};std::ofstream design(out/"measurement_design.csv");design<<std::setprecision(17)<<"design,elapsed_s,q0_sensitivity,v0_sensitivity,true_displacement_m\n";
    for(int j=1;j<=4;++j){double tau=j*std::numbers::pi/frequency;auto row=propagator(tau)[0];blind=std::max(blind,std::abs(row[1]));poor=add(poor,outer(row));design<<"half_period,"<<tau<<','<<row[0]<<','<<row[1]<<','<<exact(tau,truth)[0]<<'\n';}
    for(double tau:elapsed){auto row=propagator(tau)[0];design<<"well_spaced,"<<tau<<','<<row[0]<<','<<row[1]<<','<<exact(tau,truth)[0]<<'\n';}
    check("half-period design blind to initial velocity",blind,1e-15);
    check("well-spaced information positive definite",eigen_min(H)>0?0.:1.,0.);
    report("half_period_velocity_sensitivity_max",blind);
    Mat<2,2> scaledH{};for(int i=0;i<2;++i)for(int j=0;j<2;++j)scaledH[i][j]=scale[i]*H[i][j]*scale[j];
    double small=eigen_min(scaledH),large=scaledH[0][0]+scaledH[1][1]-small;report("scaled_information_condition",large/small);
    std::ofstream observations_file(out/"observations.csv");observations_file<<std::setprecision(17)<<"elapsed_s,observed_displacement_m,true_displacement_m,sigma_m\n";
    for(std::size_t j=0;j<elapsed.size();++j)observations_file<<elapsed[j]<<','<<observation(j)<<','<<exact(elapsed[j],truth)[0]<<','<<sigma<<'\n';
    constexpr std::size_t components=13;
    const std::array<std::string,components> names{"state_q","state_v","tangent","nodal_adjoint","control_q","control_v","running_cost","running_gradient","running_costate","running_hessian","nodal_cost","nodal_gradient","nodal_hessian"};
    std::ofstream controls(out/"control_convergence.csv");controls<<std::setprecision(17)<<"theta,steps,h_s,discrete_q0_m,discrete_v0_m_s,analytical_q0_m,analytical_v0_m_s\n";
    std::ofstream convergence(out/"convergence.csv"),orders(out/"orders.csv");convergence<<std::setprecision(17)<<"theta,steps,h_s";for(auto name:names)convergence<<','<<name;convergence<<'\n';orders<<std::setprecision(17)<<"theta,component,order\n";
    for(double theta:{.5,.65}){
        std::array<double,components> previous{},errors{};
        for(std::size_t steps:{32u,64u,128u,256u,512u,1024u}){
            auto sys=system(steps,theta);auto ev=sys.evaluate(probe);auto optimum=fit(sys);auto fitted=sys.evaluate(optimum.control);
            controls<<theta<<','<<steps<<','<<horizon/steps<<','<<optimum.control[0]<<','<<optimum.control[1]<<','<<reference[0]<<','<<reference[1]<<'\n';
            if(steps==1024){report("finest_q0_theta_"+std::to_string(theta),optimum.control[0]);report("finest_v0_theta_"+std::to_string(theta),optimum.control[1]);}
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
                if(n){auto old=ev.trajectory.state[n-1],next=ev.trajectory.state[n],delta=next-old,mid=(1-theta)*old+theta*next;double h=sys.grid().time[n]-sys.grid().time[n-1];double identity_error=energy(next)-energy(old)+h*damping*mid[1]*mid[1]+(theta-.5)*(stiffness*delta[0]*delta[0]+mass*delta[1]*delta[1]);balance=std::max(balance,std::abs(identity_error));}
            }
            // Closed-form 2x2 theta step, independent of the forward Newton solve.
            Vec<2> independent=probe;double step_error=0.,h=horizon/steps;
            double b00=1.,b01=-h*theta,b10=h*theta*stiffness/mass,b11=1+h*theta*damping/mass,det=b00*b11-b01*b10;
            for(std::size_t n=1;n<=steps;++n){auto rhs=independent+h*(1-theta)*mv(A,independent);independent={(b11*rhs[0]-b01*rhs[1])/det,(-b10*rhs[0]+b00*rhs[1])/det};step_error=std::max(step_error,normalized(independent-ev.trajectory.state[n]));}
            check("Newton versus explicit rational theta map",step_error,2e-12);
            check("theta mechanical energy balance",balance,2e-14);
            check("terminal running costate zero",norm(costate.back()),1e-15);
            previous=errors;
            double Jr=.5*dot(probe,mv(Q,probe));auto gr=mv(Q,probe),g=exact_gradient(probe);
            errors={sq,sv,st,sa/adjscale,std::abs(optimum.control[0]-reference[0])/scale[0],std::abs(optimum.control[1]-reference[1])/scale[1],std::abs(re.cost-Jr)/Jr,norm(re.gradient-gr)/norm(gr),sr/runscale,matrix_error(Qd,Q)/std::max(norm(Q[0]),norm(Q[1])),std::abs(ev.cost-exact_cost(probe))/exact_cost(probe),norm(ev.gradient-g)/norm(g),matrix_error(Hd,H)/std::max(norm(H[0]),norm(H[1]))};
            convergence<<theta<<','<<steps<<','<<horizon/steps;for(auto error:errors)convergence<<','<<error;convergence<<'\n';
            if(theta==.5&&steps==128){
                double epsilon=1e-4;Vec<2> direction{.03,-.04};auto ep=sys.evaluate(probe+epsilon*direction),em=sys.evaluate(probe-epsilon*direction);
                check("adjoint cost directional difference",std::abs((ep.cost-em.cost)/(2*epsilon)-dot(ev.gradient,direction)),2e-6);
                check("incremental adjoint directional difference",norm((1/(2*epsilon))*(ep.gradient-em.gradient)-sys.hessian_vector(probe,direction)),5e-6);
                double paired=0.;auto ts=sys.tangent(probe,ev.trajectory,direction);for(std::size_t j=0;j<elapsed.size();++j){auto n=node_at(sys.grid(),elapsed[j]);paired+=(ev.trajectory.state[n][0]-observation(j))*ts[n][0]/(sigma*sigma);}check("tangent-adjoint duality",std::abs(paired-dot(ev.gradient,direction)),2e-9);
            }
            if(theta==.5&&steps==128){
                std::ofstream curves(out/"trajectories.csv");curves<<std::setprecision(17)<<"elapsed_s,true_q_m,true_v_m_s,analytical_fit_q_m,analytical_fit_v_m_s,discrete_q_m,discrete_v_m_s,true_energy_J,discrete_energy_J,discrete_probe_running_costate_q,exact_probe_running_costate_q,discrete_probe_running_costate_v,exact_probe_running_costate_v\n";
                for(std::size_t n=0;n<=steps;++n){double tau=sys.grid().time[n]-t0;auto xt=exact(tau,truth),xa=exact(tau,reference),xd=fitted.trajectory.state[n],ra=running_costate(tau,probe);curves<<tau<<','<<xt[0]<<','<<xt[1]<<','<<xa[0]<<','<<xa[1]<<','<<xd[0]<<','<<xd[1]<<','<<energy(xt)<<','<<energy(xd)<<','<<costate[n][0]<<','<<ra[0]<<','<<costate[n][1]<<','<<ra[1]<<'\n';}
                std::ofstream trace(out/"optimizer_trace.csv");trace<<std::setprecision(17)<<"iteration,cost,gradient_norm,q0_m,v0_m_s\n";for(std::size_t j=0;j<optimum.trace.size();++j){auto it=optimum.trace[j];trace<<j<<','<<it.cost<<','<<it.gradient_norm<<','<<it.control[0]<<','<<it.control[1]<<'\n';}
                report("primary_discrete_displacement_m",optimum.control[0]);report("primary_discrete_velocity_m_s",optimum.control[1]);
            }
        }
        for(std::size_t j=0;j<components;++j){double order=std::log2(previous[j]/errors[j]);orders<<theta<<','<<names[j]<<','<<order<<'\n';check("temporal order "+names[j],std::abs(order-(theta==.5?2.:1.)),theta==.5?.06:.20);std::cout<<"order "<<theta<<" / "<<names[j]<<": "<<order<<'\n';}
    }
    std::ofstream blind_curves(out/"hidden_mode.csv");blind_curves<<std::setprecision(17)<<"elapsed_s,q_vminus_m,q_vtrue_m,q_vplus_m\n";
    for(int n=0;n<=800;++n){double tau=horizon*n/800.;blind_curves<<tau<<','<<exact(tau,{truth[0],-.03})[0]<<','<<exact(tau,truth)[0]<<','<<exact(tau,{truth[0],.25})[0]<<'\n';}
    namespace oq=oscillator::qualification;
    std::ofstream decision(out/"qualification_decisions.csv"),envelopes(out/"qualification_envelopes.csv"),polygons(out/"qualification_polygons.csv");
    decision<<std::setprecision(17)<<"case,goal,units,tolerance,radius,node_time_s,first_interval_time_s,continuous_time_s,analytical_time_s,qualified,joint_qualified,joint_time_s,blocking_goals\n";
    envelopes<<std::setprecision(17)<<"case,elapsed_s,interval_end_s,goal,tolerance,node_bound,whole_interval_bound,analytical_bound,local95_indicator,numerical_state_allowance\n";
    polygons<<std::setprecision(17)<<"case,set,vertex,q0_m,v0_m_s\n";
    struct Specification {const char* name;std::size_t steps;double theta;bool gaussian;};
    std::cout<<"\nOSCILLATOR ENGINEERING STATE-QUALIFICATION DETERMINATION\n"
             <<"Goals: displacement <= 0.001 m; velocity <= 0.002 m/s; support-force error <= 0.04 N.\n"
             <<"Known initial domain: |q0| <= 0.08 m, |v0| <= 0.30 m/s.\n"
             <<"Whole intervals through elapsed 8 s; numerical reconstruction error included.\n";
    double primary_allowance=0.,refined_allowance=0.;
    for(auto spec:std::array<Specification,7>{{{"bounded_primary",128,.5,false},{"bounded_refined",1024,.5,false},{"bounded_first_order",128,.65,false},{"bounded_first_order_refined",1024,.65,false},{"bounded_first_order_4096",4096,.65,false},{"bounded_first_order_8192",8192,.65,false},{"joint95_refined",1024,.5,true}}}){
        double radius=spec.gaussian?independent_gaussian_box_radius(sigma,elapsed.size(),.95):sigma;
        auto qs=system(spec.steps,spec.theta);
        auto result=oq::experiment(qs,radius,[&](const std::string& name,double error,double tolerance){check(std::string(spec.name)+": "+name,error,tolerance);});
        bool expected_joint=std::string(spec.name)=="bounded_refined"||std::string(spec.name)=="joint95_refined"||std::string(spec.name)=="bounded_first_order_8192";
        check(std::string(spec.name)+" expected joint engineering decision",static_cast<bool>(result.joint)==expected_joint?0.:1.,0.);
        if(std::string(spec.name)=="bounded_refined")check("refined bounded-data qualification begins at elapsed zero",result.joint.value_or(-1.),0.);
        if(spec.gaussian){bool rebound=false;for(std::size_t n=0;n<spec.steps;++n){double tau=qs.grid().time[n]-t0;if(tau>=result.first_interval_time[0].value_or(horizon)&&tau<result.continuous_time[0].value_or(0.)&&result.intervals[n][0]>oq::tolerances[0])rebound=true;}check("oscillatory first crossing is not persistent qualification",rebound?0.:1.,0.);}
        std::string blockers;for(std::size_t k=0;k<3;++k)if(!result.continuous_time[k]){if(!blockers.empty())blockers+='+';blockers+=oq::names[k];}if(blockers.empty())blockers="none";
        auto write_time=[](std::ostream& stream,const std::optional<double>& value){if(value)stream<<*value;};
        std::cout<<"\n"<<spec.name<<(spec.gaussian?" (95% simultaneous Gaussian noise box)":" (deterministic measurement-error bounds)")<<"\n";
        for(std::size_t k=0;k<3;++k){decision<<spec.name<<','<<oq::names[k]<<','<<oq::units[k]<<','<<oq::tolerances[k]<<','<<radius<<',';write_time(decision,result.node_time[k]);decision<<',';write_time(decision,result.first_interval_time[k]);decision<<',';write_time(decision,result.continuous_time[k]);decision<<',';write_time(decision,result.analytical_time[k]);decision<<','<<static_cast<bool>(result.continuous_time[k])<<','<<static_cast<bool>(result.joint)<<',';write_time(decision,result.joint);decision<<','<<blockers<<'\n';
            std::cout<<"  "<<oq::names[k]<<": ";if(result.continuous_time[k])std::cout<<"QUALIFIED from "<<*result.continuous_time[k]<<" s through the horizon.\n";else std::cout<<"NOT QUALIFIED over the requested horizon.\n";
            if(result.analytical_time[k]){double maximum=0.;for(int j=0;j<=2000;++j){double tau=*result.analytical_time[k]+(horizon-*result.analytical_time[k])*j/2000.;maximum=std::max(maximum,oq::analytical_spread(tau,result.reference_polygon,analytical_fit(),oq::goals[k])/oq::tolerances[k]);}check(std::string(spec.name)+" analytical persistent-time suffix",maximum,1.+2e-12);}
        }
        std::cout<<"  ALL REQUESTED GOALS: ";if(result.joint)std::cout<<"QUALIFIED from "<<*result.joint<<" s.\n";else std::cout<<"NOT QUALIFIED; blocking goals: "<<blockers<<".\n";
        for(std::size_t n=0;n<qs.grid().time.size();++n)for(std::size_t k=0;k<3;++k){double tau=qs.grid().time[n]-t0;envelopes<<spec.name<<','<<tau<<','<<(n<spec.steps?qs.grid().time[n+1]-t0:tau)<<','<<oq::names[k]<<','<<oq::tolerances[k]<<','<<result.nodes[n][k]<<',';if(n<spec.steps)envelopes<<result.intervals[n][k];envelopes<<','<<oq::analytical_spread(tau,result.reference_polygon,analytical_fit(),oq::goals[k])<<','<<result.local95[n][k]<<','<<result.numerical[n]<<'\n';}
        for(bool analytical:{false,true}){const auto& p=analytical?result.reference_polygon:result.polygon;for(std::size_t j=0;j<p.size();++j)polygons<<spec.name<<','<<(analytical?"analytical":"numerical_outer")<<','<<j<<','<<p[j][0]<<','<<p[j][1]<<'\n';}
        if(std::string(spec.name)=="bounded_primary")primary_allowance=*std::max_element(result.numerical.begin(),result.numerical.end());
        if(std::string(spec.name)=="bounded_refined")refined_allowance=*std::max_element(result.numerical.begin(),result.numerical.end());
    }
    report("qualification_numerical_allowance_refinement_order",std::log(primary_allowance/refined_allowance)/std::log(8.));
    check("midpoint reconstruction allowance improves quadratically",std::abs(std::log(primary_allowance/refined_allowance)/std::log(8.)-2.),.08);
    report("checks_passed",static_cast<double>(check.count));std::cout<<"Oscillator numerical and qualification checks PASS. Final decisions above are separate from test success.\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
