#include "oscillator_model.hpp"
#include "oscillator_qualification_model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

// A theta-consistent optimality-system convergence test on a nested,
// observation-aligned, nonuniform mesh family.
namespace {
using namespace oscillator;
namespace oq=oscillator::qualification;
inline constexpr Vec<2> probe{.013,.085};
inline constexpr Vec<2> tangent_direction{.007,-.031};
inline constexpr std::array<double,8> anchors{
    t0,t0+elapsed[0],t0+elapsed[1],t0+elapsed[2],
    t0+elapsed[3],t0+elapsed[4],t0+elapsed[5],tf};

Grid graded_grid(std::size_t subdivisions,double theta){
    if(!subdivisions)throw std::invalid_argument("Positive subdivision count required");
    std::vector<double> time;time.reserve((anchors.size()-1)*subdivisions+1);
    time.push_back(anchors.front());
    for(std::size_t segment=0;segment+1<anchors.size();++segment){
        const double left=anchors[segment],length=anchors[segment+1]-left;
        for(std::size_t k=1;k<=subdivisions;++k){
            const double s=static_cast<double>(k)/static_cast<double>(subdivisions);
            time.push_back(k==subdivisions?anchors[segment+1]:left+length*s*s);
        }
    }
    std::vector<double> weights(time.size()-1,theta);
    return Grid(std::move(time),std::move(weights));
}

double state_error(const std::vector<Vec<2>>& numerical,const Grid& grid,
                   const Vec<2>& initial){
    double error=0.;
    for(std::size_t n=0;n<grid.time.size();++n){
        const double tau=grid.time[n]-t0;
        const auto difference=numerical[n]-exact(tau,initial);
        error=std::max(error,std::max(std::abs(difference[0])/.02,
                                     std::abs(difference[1])/.12));
    }
    return error;
}

double matrix_error(const Mat<2,2>& a,const Mat<2,2>& b){
    return std::max(norm(a[0]-b[0]),norm(a[1]-b[1]));
}

struct Errors {
    double state{},tangent{},adjoint{},cost{},gradient{},hessian{},fit{};
    double running_cost{},running_gradient{},running_costate{},running_hessian{};
};

void verify(const std::string& name,double error,double tolerance,
            std::ofstream& csv,std::size_t& checks){
    const bool passed=std::isfinite(error)&&error<=tolerance;
    csv<<name<<','<<std::setprecision(17)<<error<<','<<tolerance<<','<<passed<<'\n';
    ++checks;
    if(!passed)throw std::runtime_error(name+" failed: "+std::to_string(error));
}

Errors evaluate_family(const System<2,2>& system,const Vec<2>& exact_fit){
    const auto& grid=system.grid();
    const auto evaluation=system.evaluate(probe);
    const auto tangent=system.tangent(probe,evaluation.trajectory,tangent_direction);
    const auto reference_gradient=exact_gradient(probe);
    const auto reference_hessian=information();
    Errors e;
    e.state=state_error(evaluation.trajectory.state,grid,probe);
    double tangent_scale=0.,adjoint_scale=0.;
    for(std::size_t n=0;n<grid.time.size();++n){
        const double tau=grid.time[n]-t0;
        const auto tangent_reference=mv(propagator(tau),tangent_direction);
        e.tangent=std::max(e.tangent,std::max(
            std::abs(tangent[n][0]-tangent_reference[0])/.02,
            std::abs(tangent[n][1]-tangent_reference[1])/.12));
        const auto adjoint_reference=nodal_costate(tau,probe);
        const auto adjoint_delta=evaluation.nodal_adjoint[n]-adjoint_reference;
        e.adjoint=std::max(e.adjoint,norm(Vec<2>{.02*adjoint_delta[0],
                                                  .12*adjoint_delta[1]}));
        adjoint_scale=std::max(adjoint_scale,norm(Vec<2>{.02*adjoint_reference[0],
                                                         .12*adjoint_reference[1]}));
        tangent_scale=std::max(tangent_scale,norm(Vec<2>{tangent_reference[0]/.02,
                                                         tangent_reference[1]/.12}));
    }
    e.tangent/=std::max(tangent_scale,1e-30);
    e.adjoint/=std::max(adjoint_scale,1e-30);
    const double reference_cost=exact_cost(probe);
    e.cost=std::abs(evaluation.cost-reference_cost)/std::max(1.,std::abs(reference_cost));
    e.gradient=norm(evaluation.gradient-reference_gradient)/std::max(1.,norm(reference_gradient));
    Mat<2,2> discrete_hessian{};
    for(std::size_t j=0;j<2;++j){
        Vec<2> basis{};basis[j]=1.;const auto column=system.hessian_vector(probe,basis);
        for(std::size_t i=0;i<2;++i)discrete_hessian[i][j]=column[i];
    }
    e.hessian=matrix_error(discrete_hessian,reference_hessian)/
        std::max(norm(reference_hessian[0]),norm(reference_hessian[1]));
    const auto optimum=fit(system);
    e.fit=norm(optimum.control-exact_fit)/std::max(1.,norm(exact_fit));

    auto running_system=System<2,2>(problem(grid,true),grid);
    const auto running_eval=running_system.evaluate(probe);
    const auto running_costates=running_system.cost_to_go_adjoint(probe,running_eval);
    const auto gramian=oscillator::gramian(horizon);
    const double exact_running_cost=.5*dot(probe,mv(gramian,probe));
    const auto exact_running_gradient=mv(gramian,probe);
    e.running_cost=std::abs(running_eval.cost-exact_running_cost)/std::max(1.,std::abs(exact_running_cost));
    e.running_gradient=norm(running_eval.gradient-exact_running_gradient)/
        std::max(1.,norm(exact_running_gradient));
    double running_costate_scale=0.;
    for(std::size_t n=0;n<grid.time.size();++n){
        const double tau=grid.time[n]-t0;
        const auto reference=running_costate(tau,probe);
        const auto difference=running_costates[n]-reference;
        e.running_costate=std::max(e.running_costate,
            norm(Vec<2>{.02*difference[0],.12*difference[1]}));
        running_costate_scale=std::max(running_costate_scale,
            norm(Vec<2>{.02*reference[0],.12*reference[1]}));
    }
    e.running_costate/=std::max(running_costate_scale,1e-30);
    Mat<2,2> running_hessian{};
    for(std::size_t j=0;j<2;++j){
        Vec<2> basis{};basis[j]=1.;const auto column=running_system.hessian_vector(probe,basis);
        for(std::size_t i=0;i<2;++i)running_hessian[i][j]=column[i];
    }
    e.running_hessian=matrix_error(running_hessian,gramian)/
        std::max(norm(gramian[0]),norm(gramian[1]));
    return e;
}

double observed_order(double coarse,double fine,double h_coarse,double h_fine){
    return std::log(coarse/fine)/std::log(h_coarse/h_fine);
}
}

int main(int argc,char** argv){try{
    const std::filesystem::path output=argc>1?argv[1]:"results/nonuniform_oscillator";
    std::filesystem::create_directories(output);
    std::ofstream convergence(output/"convergence.csv"),orders(output/"orders.csv"),
        checks_csv(output/"checks.csv"),mesh_csv(output/"mesh.csv"),
        summary(output/"summary.csv"),qualification_csv(output/"qualification.csv"),
        uncertainty_csv(output/"state_uncertainty.csv"),
        envelopes_csv(output/"qualification_envelopes.csv");
    convergence<<std::setprecision(17)<<"theta,subdivisions,steps,h_min_s,h_max_s,mesh_ratio"
        ",state,tangent,nodal_adjoint,point_cost,point_gradient,point_hessian,"
        "recovered_initial_state,running_cost,running_gradient,running_costate,running_hessian\n";
    orders<<"theta,quantity,observed_order,expected_order,passed\n";
    checks_csv<<"check,error,tolerance,passed\n";
    mesh_csv<<std::setprecision(17)<<"theta,subdivisions,node,elapsed_s,step_s\n";
    summary<<std::setprecision(17)<<"quantity,value\n";
    qualification_csv<<std::setprecision(17)
        <<"theta,goal,units,tolerance,qualified,persistent_time_s,joint_qualified,joint_time_s\n";
    uncertainty_csv<<std::setprecision(17)
        <<"theta,elapsed_s,true_displacement_m,center_displacement_m,displacement_bound_m,"
        <<"true_velocity_m_s,center_velocity_m_s,velocity_bound_m_s,"
        <<"true_support_force_N,center_support_force_N,support_force_bound_N\n";
    envelopes_csv<<std::setprecision(17)
        <<"theta,goal,units,tolerance,elapsed_s,interval_end_s,node_bound,whole_interval_bound\n";
    std::size_t checks=0;
    const Vec<2> analytic_fit=analytical_fit();
    const std::array<std::string,11> names{"state","tangent","nodal_adjoint",
        "point_cost","point_gradient","point_hessian","recovered_initial_state",
        "running_cost","running_gradient","running_costate","running_hessian"};
    const std::array<std::size_t,4> subdivisions{8,16,32,64};
    for(double theta:{.5,.65}){
        std::vector<Errors> errors;
        std::vector<double> hmaxes;
        std::vector<double> previous_mesh;
        for(auto m:subdivisions){
            const auto grid=graded_grid(m,theta);
            double hmin=std::numeric_limits<double>::infinity(),hmax=0.;
            for(std::size_t n=0;n<grid.steps();++n){
                const double h=grid.time[n+1]-grid.time[n];hmin=std::min(hmin,h);hmax=std::max(hmax,h);
                mesh_csv<<theta<<','<<m<<','<<n<<','<<grid.time[n]-t0<<','<<h<<'\n';
            }
            mesh_csv<<theta<<','<<m<<','<<grid.steps()<<','<<grid.time.back()-t0<<",\n";
            verify("mesh is genuinely nonuniform theta="+std::to_string(theta)+" m="+std::to_string(m),
                   hmax/hmin>10.?0.:1.,0.,checks_csv,checks);
            for(double old_node:previous_mesh){
                const auto it=std::lower_bound(grid.time.begin(),grid.time.end(),old_node);
                const bool retained=it!=grid.time.end()&&std::abs(*it-old_node)<2e-14;
                verify("nested mesh retains coarse node theta="+std::to_string(theta)+" m="+std::to_string(m),
                       retained?0.:1.,0.,checks_csv,checks);
            }
            previous_mesh=grid.time;
            for(double observation_time:elapsed){
                const auto node=node_at(grid,observation_time);
                verify("fixed observation node theta="+std::to_string(theta)+" m="+std::to_string(m),
                    std::abs(grid.time[node]-(t0+observation_time)),2e-14,checks_csv,checks);
            }
            const auto system=System<2,2>(problem(grid),grid);
            auto e=evaluate_family(system,analytic_fit);
            verify("forward residual theta="+std::to_string(theta)+" m="+std::to_string(m),
                   system.forward(probe).max_residual,2e-12,checks_csv,checks);
            const auto fitted=fit(system);
            verify("discrete fit stationarity theta="+std::to_string(theta)+" m="+std::to_string(m),
                   norm(system.evaluate(fitted.control).gradient),2e-8,checks_csv,checks);
            const std::array<double,11> values{e.state,e.tangent,e.adjoint,e.cost,e.gradient,
                e.hessian,e.fit,e.running_cost,e.running_gradient,e.running_costate,e.running_hessian};
            convergence<<theta<<','<<m<<','<<grid.steps()<<','<<hmin<<','<<hmax<<','<<hmax/hmin;
            for(double value:values)convergence<<','<<value;
            convergence<<'\n';errors.push_back(e);hmaxes.push_back(hmax);
            if(m==64){
                const auto qualified=oq::experiment(system,sigma,
                    [&](const std::string& name,double error,double tolerance){
                        verify("qualification: "+name+" theta="+std::to_string(theta),
                               error,tolerance,checks_csv,checks);
                    });
                const auto reconstructed=system.forward(qualified.center);
                for(std::size_t k=0;k<oq::goals.size();++k){
                    const auto time=qualified.continuous_time[k];
                    qualification_csv<<theta<<','<<oq::names[k]<<','<<oq::units[k]<<','
                        <<oq::tolerances[k]<<','<<static_cast<bool>(time)<<',';
                    if(time)qualification_csv<<*time;
                    qualification_csv<<','<<static_cast<bool>(qualified.joint)<<',';
                    if(qualified.joint)qualification_csv<<*qualified.joint;
                    qualification_csv<<'\n';
                    std::cout<<"  qualification theta="<<theta<<" / "<<oq::names[k]<<": ";
                    if(time)std::cout<<"QUALIFIED from "<<*time<<" s\n";
                    else std::cout<<"NOT QUALIFIED over the horizon\n";
                }
                for(std::size_t n=0;n<grid.time.size();++n){
                    const double tau=grid.time[n]-t0;
                    const auto truth_state=exact(tau,truth);
                    const auto center_state=reconstructed.state[n];
                    const double true_force=stiffness*truth_state[0]+damping*truth_state[1];
                    const double center_force=stiffness*center_state[0]+damping*center_state[1];
                    uncertainty_csv<<theta<<','<<tau<<','<<truth_state[0]<<','<<center_state[0]<<','
                        <<qualified.nodes[n][0]<<','<<truth_state[1]<<','<<center_state[1]<<','
                        <<qualified.nodes[n][1]<<','<<true_force<<','<<center_force<<','
                        <<qualified.nodes[n][2]<<'\n';
                    for(std::size_t k=0;k<oq::goals.size();++k){
                        envelopes_csv<<theta<<','<<oq::names[k]<<','<<oq::units[k]<<','
                            <<oq::tolerances[k]<<','<<tau<<',';
                        if(n<grid.steps())envelopes_csv<<grid.time[n+1]-t0;
                        else envelopes_csv<<tau;
                        envelopes_csv<<','<<qualified.nodes[n][k]<<',';
                        if(n<grid.steps())envelopes_csv<<qualified.intervals[n][k];
                        envelopes_csv<<'\n';
                    }
                }
            }
            if(m==16){
                const Vec<2> direction{.03,-.04};const double epsilon=1e-5;
                const auto ev=system.evaluate(probe);
                const auto plus=system.evaluate(probe+epsilon*direction);
                const auto minus=system.evaluate(probe-epsilon*direction);
                double tangent_pairing=0.;const auto tangent=system.tangent(probe,ev.trajectory,direction);
                for(std::size_t j=0;j<elapsed.size();++j){
                    const auto n=node_at(grid,elapsed[j]);
                    tangent_pairing+=(ev.trajectory.state[n][0]-observation(j))
                        *tangent[n][0]/(sigma*sigma);
                }
                verify("nonuniform tangent-adjoint pairing theta="+std::to_string(theta),
                    std::abs(tangent_pairing-dot(ev.gradient,direction)),2e-9,checks_csv,checks);
                verify("nonuniform discrete objective gradient theta="+std::to_string(theta),
                    std::abs((plus.cost-minus.cost)/(2*epsilon)-dot(ev.gradient,direction)),
                    3e-6,checks_csv,checks);
                verify("nonuniform incremental-adjoint Hessian theta="+std::to_string(theta),
                    norm(system.hessian_vector(probe,direction)-
                        (1/(2*epsilon))*(plus.gradient-minus.gradient)),8e-6,checks_csv,checks);
                const auto run=System<2,2>(problem(grid,true),grid);
                const auto re=run.evaluate(probe);
                verify("nonuniform running costate terminal zero theta="+std::to_string(theta),
                    norm(run.cost_to_go_adjoint(probe,re).back()),1e-14,checks_csv,checks);
            }
        }
        const double expected=theta==.5?2.:1.;
        for(std::size_t j=0;j<names.size();++j){
            const auto member=[j](const Errors& e){
                const std::array<double,11> a{e.state,e.tangent,e.adjoint,e.cost,e.gradient,
                    e.hessian,e.fit,e.running_cost,e.running_gradient,e.running_costate,e.running_hessian};
                return a[j];
            };
            const double order=observed_order(member(errors[2]),member(errors[3]),hmaxes[2],hmaxes[3]);
            const bool passed=std::isfinite(order)&&std::abs(order-expected)<(expected==2.?0.35:0.3);
            orders<<theta<<','<<names[j]<<','<<order<<','<<expected<<','<<passed<<'\n';
            summary<<"order_theta_"<<theta<<'_'<<names[j]<<','<<order<<'\n';
            verify("nonuniform observed order theta="+std::to_string(theta)+" "+names[j],
                std::abs(order-expected),expected==2.?0.35:0.3,checks_csv,checks);
            std::cout<<"theta="<<theta<<" nonuniform order / "<<names[j]<<": "<<order<<" (expected "<<expected<<")\n";
        }
    }
    summary<<"checks_passed,"<<checks<<'\n';
    std::cout<<"Nonuniform oscillator optimality-system test PASS: "<<checks
             <<" checks. Fixed observation times; nested quadratic mesh grading.\n"
             <<"Results: \""<<output.string()<<"\"\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
