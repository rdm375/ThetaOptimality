#include <theta_opt/theta_opt.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Lotka--Volterra inverse problems, solved in three stages: initial populations,
// kinetic rates, and all six quantities jointly. The independent RK4 reference
// is also checked against the model's conserved quantity.
namespace {
using namespace theta_opt;
inline constexpr double t0=0.,tf=10.,reference_step=2.e-4;
inline constexpr std::array<double,2> truth_initial{1.6,.6};
inline constexpr std::array<double,4> truth_rates{1.0,.8,.6,.7};
inline constexpr std::array<double,2> assay_sd{.025,.02};
inline constexpr std::size_t sample_count=20;

using State=std::array<double,2>;
using Rates=std::array<double,4>;
State rhs(const State& x,const Rates& p) {
    const auto [a,b,c,d]=p;
    return {a*x[0]-b*x[0]*x[1],c*x[0]*x[1]-d*x[1]};
}
State add_scaled(const State& x,double h,const State& v) {
    return {x[0]+h*v[0],x[1]+h*v[1]};
}
State rk4_step(const State& x,const Rates& p,double h) {
    const auto k1=rhs(x,p);
    const auto k2=rhs(add_scaled(x,.5*h,k1),p);
    const auto k3=rhs(add_scaled(x,.5*h,k2),p);
    const auto k4=rhs(add_scaled(x,h,k3),p);
    return {x[0]+h*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6.,
            x[1]+h*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6.};
}
double invariant(const State& x,const Rates& p) {
    const auto [a,b,c,d]=p;
    return c*x[0]-d*std::log(x[0])+b*x[1]-a*std::log(x[1]);
}
State interpolate(const std::vector<State>& values,double time) {
    const double coordinate=(time-t0)/reference_step;
    const auto lower=static_cast<std::size_t>(std::floor(coordinate));
    if(lower>=values.size()-1)return values.back();
    const double fraction=coordinate-static_cast<double>(lower);
    return {(1.-fraction)*values[lower][0]+fraction*values[lower+1][0],
            (1.-fraction)*values[lower][1]+fraction*values[lower+1][1]};
}
struct Reference {std::vector<State> dense,observed;};
Reference make_reference() {
    const std::size_t count=static_cast<std::size_t>((tf-t0)/reference_step);
    if(std::abs(reference_step*static_cast<double>(count)-(tf-t0))>1e-12)
        throw std::logic_error("Reference step must divide horizon");
    Reference result;result.dense.resize(count+1);result.dense[0]={truth_initial[0],truth_initial[1]};
    for(std::size_t k=0;k<count;++k)result.dense[k+1]=rk4_step(result.dense[k],truth_rates,reference_step);
    const std::size_t stride=static_cast<std::size_t>(.5/reference_step);
    if(std::abs(reference_step*static_cast<double>(stride)-.5)>1e-12)
        throw std::logic_error("Reference step must divide assay spacing");
    for(std::size_t j=1;j<=sample_count;++j)result.observed.push_back(result.dense[j*stride]);
    return result;
}
State assay(std::size_t j,const Reference& reference) {
    const double phase=static_cast<double>(j+1);
    return {reference.observed.at(j)[0]+assay_sd[0]*(.35*std::sin(1.7*phase)-.15*std::cos(.6*phase)),
            reference.observed.at(j)[1]+assay_sd[1]*(.30*std::cos(1.3*phase)+.12*std::sin(.8*phase))};
}
std::size_t sample_node(const Grid& grid,std::size_t j) {
    const double t=.5*static_cast<double>(j+1);
    for(std::size_t n=0;n<grid.time.size();++n)
        if(std::abs(grid.time[n]-t)<1e-11)return n;
    throw std::logic_error("Assay time is not a theta-grid node");
}

template<std::size_t P>
Problem<2,P> make_problem(const Grid& grid,int mode,const Reference& reference) {
    // mode 0: initial populations unknown; mode 1: rates unknown;
    // mode 2: all initial populations and rates unknown.
    const bool variable_initial=mode!=1;
    const bool variable_rates=mode!=0;
    const std::size_t rate_offset=mode==1?0:2;
    Problem<2,P> problem;
    problem.initial.value=[variable_initial](const Vec<P>& z){
        return variable_initial?Vec<2>{z[0],z[1]}:Vec<2>{truth_initial[0],truth_initial[1]};
    };
    problem.initial.jacobian=[variable_initial](const Vec<P>&){
        Mat<2,P> J{};
        if(variable_initial){J[0][0]=1.;J[1][1]=1.;}
        return J;
    };
    problem.initial.directional_jacobian=[](const Vec<P>&,const Vec<P>&){return Mat<2,P>{};};

    const auto rates=[variable_rates,rate_offset](const Vec<P>& z){
        Rates p=truth_rates;
        if(variable_rates)for(std::size_t j=0;j<4;++j)p[j]=z[rate_offset+j];
        return p;
    };
    problem.dynamics.value=[rates](double,const Vec<2>& x,const Vec<P>& z){
        const auto p=rates(z);const auto f=rhs({x[0],x[1]},p);return Vec<2>{f[0],f[1]};
    };
    problem.dynamics.derivatives=[variable_rates,rate_offset,rates](double,const Vec<2>& x,const Vec<P>& z){
        const auto p=rates(z);const auto [a,b,c,d]=p;Derivatives<2,P> D;
        D.x={{{a-b*x[1],-b*x[0]},{c*x[1],c*x[0]-d}}};
        if(variable_rates){
            D.control[0][rate_offset]=x[0];D.control[0][rate_offset+1]=-x[0]*x[1];
            D.control[1][rate_offset+2]=x[0]*x[1];D.control[1][rate_offset+3]=-x[1];
        }
        return D;
    };
    problem.dynamics.directional_derivatives=[variable_rates,rate_offset,rates](double,const Vec<2>& x,const Vec<P>& z,
                                                                                 const Vec<2>& dx,const Vec<P>& dz){
        const auto p=rates(z);const auto [a,b,c,d]=p;double da=0.,db=0.,dc=0.,dd=0.;
        if(variable_rates){da=dz[rate_offset];db=dz[rate_offset+1];dc=dz[rate_offset+2];dd=dz[rate_offset+3];}
        Derivatives<2,P> D;
        D.x={{{da-db*x[1]-b*dx[1],-db*x[0]-b*dx[0]},
              {dc*x[1]+c*dx[1],dc*x[0]+c*dx[0]-dd}}};
        if(variable_rates){
            D.control[0][rate_offset]=dx[0];
            D.control[0][rate_offset+1]=-dx[0]*x[1]-x[0]*dx[1];
            D.control[1][rate_offset+2]=dx[0]*x[1]+x[0]*dx[1];
            D.control[1][rate_offset+3]=-dx[1];
        }
        (void)a;(void)c;
        return D;
    };
    problem.dynamics.admissible=[](double,const Vec<2>& x,const Vec<P>&){
        return std::isfinite(x[0])&&std::isfinite(x[1])&&x[0]>0.&&x[1]>0.;
    };

    for(std::size_t j=0;j<sample_count;++j)for(std::size_t species=0;species<2;++species) {
        const double y=assay(j,reference)[species],variance=assay_sd[species]*assay_sd[species];
        auto observation=ScalarFunction<2,P>::zero();
        observation.evaluate=[species,y,variance](double,const Vec<2>& x,const Vec<P>&){
            ScalarValue<2,P> value;const double residual=x[species]-y;
            value.value=.5*residual*residual/variance;value.x[species]=residual/variance;return value;
        };
        observation.directional_gradient=[species,variance](double,const Vec<2>&,const Vec<P>&,
                                                              const Vec<2>& dx,const Vec<P>&){
            ScalarDirection<2,P> value;value.x[species]=dx[species]/variance;return value;
        };
        problem.nodal.emplace_back(sample_node(grid,j),std::move(observation));
    }
    return problem;
}

template<std::size_t P>
System<2,P> make_system(std::size_t steps,double theta,int mode,const Reference& reference) {
    auto grid=Grid::uniform(t0,tf,steps,theta);auto problem=make_problem<P>(grid,mode,reference);
    return System<2,P>(std::move(problem),std::move(grid));
}

template<std::size_t P>
struct Fit {Vec<P> control;Evaluation<2,P> evaluation;std::size_t iterations{};};

template<std::size_t P>
Fit<P> fit_mode(std::size_t steps,double theta,int mode,const Reference& reference) {
    auto system=make_system<P>(steps,theta,mode,reference);Vec<P> initial{},lower{},upper{};
    if constexpr(P==2){initial={1.35,.78};lower={.2,.1};upper={3.,2.5};}
    else if constexpr(P==4){initial={.9,.72,.52,.82};lower={.2,.15,.1,.15};upper={2.,1.6,1.5,1.5};}
    else {initial={1.35,.78,.9,.72,.52,.82};lower={.2,.1,.2,.15,.1,.15};upper={3.,2.5,2.,1.6,1.5,1.5};}
    OptimizeOptions options;options.gradient_tolerance=2e-7;options.max_iterations=80;
    options.cost_roundoff_slack=1e-10;
    auto optimum=minimize_box_constrained(system,initial,lower,upper,options);
    return {optimum.control,system.evaluate(optimum.control),optimum.trace.size()};
}

template<std::size_t P>
std::vector<State> continuous_fit_trajectory(const Vec<P>& z) {
    State initial{truth_initial[0],truth_initial[1]};Rates p=truth_rates;
    if constexpr(P==2){initial={z[0],z[1]};}
    else if constexpr(P==4){for(std::size_t j=0;j<4;++j)p[j]=z[j];}
    else {initial={z[0],z[1]};for(std::size_t j=0;j<4;++j)p[j]=z[j+2];}
    const std::size_t count=static_cast<std::size_t>((tf-t0)/reference_step);
    std::vector<State> result(count+1);result[0]=initial;
    for(std::size_t k=0;k<count;++k)result[k+1]=rk4_step(result[k],p,reference_step);
    return result;
}
double trajectory_error(const std::vector<Vec<2>>& discrete,const std::vector<State>& continuous,
                        std::size_t steps) {
    double error=0.;
    for(std::size_t n=0;n<discrete.size();++n){
        const double coordinate=static_cast<double>(n)*static_cast<double>(continuous.size()-1)/static_cast<double>(steps);
        const auto lower=static_cast<std::size_t>(std::floor(coordinate));
        const double fraction=coordinate-static_cast<double>(lower);
        const auto& lo=continuous[lower];const auto& hi=continuous[std::min(lower+1,continuous.size()-1)];
        const double dx=discrete[n][0]-((1.-fraction)*lo[0]+fraction*hi[0]);
        const double dy=discrete[n][1]-((1.-fraction)*lo[1]+fraction*hi[1]);
        error=std::max(error,std::hypot(dx,dy));
    }
    return error;
}

template<std::size_t P>
double gradient_fd_error(const System<2,P>& system,const Vec<P>& p) {
    auto exact=system.evaluate(p);double error=0.;
    for(std::size_t j=0;j<P;++j){const double h=2.e-6*std::max(1.,std::abs(p[j]));auto plus=p,minus=p;plus[j]+=h;minus[j]-=h;
        const double fd=(system.cost(plus,system.forward(plus))-system.cost(minus,system.forward(minus)))/(2*h);
        error=std::max(error,std::abs(fd-exact.gradient[j])/(1.+std::abs(fd)));
    }
    return error;
}
template<std::size_t P>
double tangent_fd_error(const System<2,P>& system,const Vec<P>& p) {
    Vec<P> direction{};for(std::size_t j=0;j<P;++j)direction[j]=.2/static_cast<double>(j+1);
    const auto trajectory=system.forward(p);const auto tangent=system.tangent(p,trajectory,direction);
    const double epsilon=2.e-6;const auto plus_p=p+epsilon*direction,minus_p=p-epsilon*direction;
    const auto plus=system.forward(plus_p),minus=system.forward(minus_p);double error=0.;
    for(std::size_t n=0;n<tangent.size();++n)for(std::size_t i=0;i<2;++i){
        const double finite_difference=(plus.state[n][i]-minus.state[n][i])/(2.*epsilon);
        error=std::max(error,std::abs(finite_difference-tangent[n][i])/(1.+std::abs(finite_difference)));
    }
    return error;
}
template<std::size_t P>
PsdSpectrum<P> information_spectrum(const System<2,P>& system,const Vec<P>& p,
                                    const Trajectory<2>& trajectory) {
    std::array<std::vector<Vec<2>>,P> sensitivities;
    for(std::size_t j=0;j<P;++j){Vec<P> direction{};direction[j]=1.;sensitivities[j]=system.tangent(p,trajectory,direction);}
    Mat<P,P> information{};
    for(std::size_t sample=0;sample<sample_count;++sample)for(std::size_t species=0;species<2;++species){
        const auto node=sample_node(system.grid(),sample);const double variance=assay_sd[species]*assay_sd[species];
        for(std::size_t i=0;i<P;++i)for(std::size_t j=0;j<P;++j)
            information[i][j]+=sensitivities[i][node][species]*sensitivities[j][node][species]/variance;
    }
    return psd_spectrum(information);
}
template<std::size_t P>
double hessian_fd_error(const System<2,P>& system,const Vec<P>& p) {
    double error=0.;
    for(std::size_t j=0;j<P;++j){Vec<P> e{};e[j]=1.;const auto hv=system.hessian_vector(p,e);const double h=3.e-5*std::max(1.,std::abs(p[j]));auto plus=p,minus=p;plus[j]+=h;minus[j]-=h;
        const auto gp=system.evaluate(plus).gradient,gm=system.evaluate(minus).gradient;
        for(std::size_t i=0;i<P;++i){const double fd=(gp[i]-gm[i])/(2*h);error=std::max(error,std::abs(fd-hv[i])/(1.+std::abs(fd)));}
    }
    return error;
}

struct CheckLog {
    std::ofstream csv;std::size_t count{};
    explicit CheckLog(const std::filesystem::path& path):csv(path) {
        if(!csv)throw std::runtime_error("Cannot open predator-prey checks CSV");
        csv<<"check,error,tolerance,passed\n";
    }
    void require(const std::string& name,double error,double tolerance) {
        const bool passed=std::isfinite(error)&&error<=tolerance;
        csv<<name<<','<<std::setprecision(17)<<error<<','<<tolerance<<','<<(passed?1:0)<<'\n';
        ++count;if(!passed)throw std::runtime_error("Predator-prey check failed: "+name+" error="+std::to_string(error));
    }
};
double order(double coarse,double fine) {
    if(!(coarse>0.)||!(fine>0.))return std::numeric_limits<double>::quiet_NaN();
    return std::log(coarse/fine)/std::log(2.);
}

void write_outputs(const std::filesystem::path& dir,const Reference& reference,CheckLog& checks) {
    std::filesystem::create_directories(dir);
    std::ofstream observations(dir/"observations.csv");observations<<"time,prey_observed,predator_observed,prey_sd,predator_sd\n";
    for(std::size_t j=0;j<sample_count;++j){auto y=assay(j,reference);observations<<.5*(j+1)<<','<<y[0]<<','<<y[1]<<','<<assay_sd[0]<<','<<assay_sd[1]<<'\n';}
    std::ofstream summary(dir/"fits.csv");summary<<"case,theta,steps,prey_initial,predator_initial,prey_growth,predation,predator_growth,predator_mortality,cost,gradient_norm,iterations\n";
    std::ofstream convergence(dir/"convergence.csv");convergence<<"case,theta,steps,h,state_error,fit_error,tangent_error,gradient_error,hessian_error\n";
    std::ofstream identifiability(dir/"identifiability.csv");identifiability<<"case,rank,controls,smallest_retained_eigenvalue,largest_eigenvalue,condition_number\n";
    const std::array<std::size_t,4> meshes{40,80,160,320};
    std::array<std::vector<Vec<2>>,3> representative_trajectories;
    for(double theta:{.5,.65})for(int mode=0;mode<3;++mode){
        if(mode==0) {
            std::vector<Fit<2>> fits;for(auto n:meshes)fits.push_back(fit_mode<2>(n,theta,mode,reference));
            const auto fine=fits.back();std::array<double,4> state_errors{};for(std::size_t k=0;k<meshes.size();++k){const double h=(tf-t0)/meshes[k];
                auto system=make_system<2>(meshes[k],theta,mode,reference);auto trajectory=system.forward(fits[k].control);double se=0.;
                const auto continuous=continuous_fit_trajectory(fits[k].control);
                se=trajectory_error(trajectory.state,continuous,meshes[k]);state_errors[k]=se;
                if(theta==.5&&k+1==meshes.size())representative_trajectories[0]=trajectory.state;
                const double fe=std::hypot(fits[k].control[0]-fine.control[0],fits[k].control[1]-fine.control[1]);
                convergence<<"initial_state,"<<theta<<','<<meshes[k]<<','<<h<<','<<se<<','<<fe<<','<<tangent_fd_error(system,fits[k].control)<<','<<gradient_fd_error(system,fits[k].control)<<','<<hessian_fd_error(system,fits[k].control)<<'\n';
                auto p=fits[k].control;summary<<"initial_state,"<<theta<<','<<meshes[k]<<','<<p[0]<<','<<p[1]<<','<<truth_rates[0]<<','<<truth_rates[1]<<','<<truth_rates[2]<<','<<truth_rates[3]<<','<<fits[k].evaluation.cost<<','<<norm(fits[k].evaluation.gradient)<<','<<fits[k].iterations<<'\n';
            }
            auto system=make_system<2>(160,theta,mode,reference);
            checks.require("initial_state_tangent_fd_theta_"+std::to_string(theta),tangent_fd_error(system,fits[2].control),1e-6);
            checks.require("initial_state_gradient_fd_theta_"+std::to_string(theta),gradient_fd_error(system,fits[2].control),1e-4);
            checks.require("initial_state_hessian_fd_theta_"+std::to_string(theta),hessian_fd_error(system,fits[2].control),1e-4);
            const double state_order=order(state_errors[2],state_errors[3]);
            const double d1=norm(fits[1].control-fits[2].control),d2=norm(fits[2].control-fits[3].control);
            const double fit_order=order(d1,d2);
            checks.require("initial_state_state_order_theta_"+std::to_string(theta),std::abs(state_order-(theta==.5?2.:1.)),.25);
            checks.require("initial_state_fit_order_theta_"+std::to_string(theta),std::abs(fit_order-(theta==.5?2.:1.)),.4);
            std::cout<<"  initial-state theta="<<theta<<" observed orders: trajectory="<<state_order<<", fitted initial populations="<<fit_order<<"\n";
            if(theta==.5) {
                checks.require("initial_state_recovery",norm(fits.back().control-Vec<2>{truth_initial[0],truth_initial[1]}),.003);
                const auto spec=information_spectrum(make_system<2>(meshes.back(),theta,mode,reference),fits.back().control,fits.back().evaluation.trajectory);
                identifiability<<"initial_state,"<<spec.rank<<",2,"<<spec.eigenvalues[spec.rank-1]<<','<<spec.eigenvalues[0]<<','<<spec.eigenvalues[0]/spec.eigenvalues[spec.rank-1]<<'\n';
                checks.require("initial_state_information_full_rank",static_cast<double>(2-spec.rank),0.);
                std::cout<<"    N=320 estimate (prey0,predator0)=("<<fits.back().control[0]<<", "<<fits.back().control[1]<<")\n";
            }
        } else if(mode==1) {
            std::vector<Fit<4>> fits;for(auto n:meshes)fits.push_back(fit_mode<4>(n,theta,mode,reference));
            const auto fine=fits.back();std::array<double,4> state_errors{};for(std::size_t k=0;k<meshes.size();++k){auto system=make_system<4>(meshes[k],theta,mode,reference);const double fe=norm(fits[k].control-fine.control);
                const double se=trajectory_error(system.forward(fits[k].control).state,continuous_fit_trajectory(fits[k].control),meshes[k]);state_errors[k]=se;
                if(theta==.5&&k+1==meshes.size())representative_trajectories[1]=system.forward(fits[k].control).state;
                convergence<<"parameters,"<<theta<<','<<meshes[k]<<','<<(tf-t0)/meshes[k]<<','<<se<<','<<fe<<','<<tangent_fd_error(system,fits[k].control)<<','<<gradient_fd_error(system,fits[k].control)<<','<<hessian_fd_error(system,fits[k].control)<<'\n';
                const auto& p=fits[k].control;summary<<"parameters,"<<theta<<','<<meshes[k]<<','<<truth_initial[0]<<','<<truth_initial[1]<<','<<p[0]<<','<<p[1]<<','<<p[2]<<','<<p[3]<<','<<fits[k].evaluation.cost<<','<<norm(fits[k].evaluation.gradient)<<','<<fits[k].iterations<<'\n';}
            auto system=make_system<4>(160,theta,mode,reference);
            checks.require("parameter_tangent_fd_theta_"+std::to_string(theta),tangent_fd_error(system,fits[2].control),1e-6);
            checks.require("parameter_gradient_fd_theta_"+std::to_string(theta),gradient_fd_error(system,fits[2].control),1e-4);
            checks.require("parameter_hessian_fd_theta_"+std::to_string(theta),hessian_fd_error(system,fits[2].control),1e-4);
            const double state_order=order(state_errors[2],state_errors[3]);
            const double fit_order=order(norm(fits[1].control-fits[2].control),norm(fits[2].control-fits[3].control));
            checks.require("parameter_state_order_theta_"+std::to_string(theta),std::abs(state_order-(theta==.5?2.:1.)),.25);
            checks.require("parameter_fit_order_theta_"+std::to_string(theta),std::abs(fit_order-(theta==.5?2.:1.)),.4);
            std::cout<<"  parameter theta="<<theta<<" observed orders: trajectory="<<state_order<<", fitted rates="<<fit_order<<"\n";
            if(theta==.5) {
                checks.require("parameter_recovery",norm(fits.back().control-Vec<4>{truth_rates[0],truth_rates[1],truth_rates[2],truth_rates[3]}),.003);
                const auto spec=information_spectrum(make_system<4>(meshes.back(),theta,mode,reference),fits.back().control,fits.back().evaluation.trajectory);
                identifiability<<"parameters,"<<spec.rank<<",4,"<<spec.eigenvalues[spec.rank-1]<<','<<spec.eigenvalues[0]<<','<<spec.eigenvalues[0]/spec.eigenvalues[spec.rank-1]<<'\n';
                checks.require("parameter_information_full_rank",static_cast<double>(4-spec.rank),0.);
                std::cout<<"    N=320 estimate (a,b,c,d)=("<<fits.back().control[0]<<", "<<fits.back().control[1]<<", "<<fits.back().control[2]<<", "<<fits.back().control[3]<<")\n";
            }
        } else {
            std::vector<Fit<6>> fits;for(auto n:meshes)fits.push_back(fit_mode<6>(n,theta,mode,reference));
            const auto fine=fits.back();std::array<double,4> state_errors{};for(std::size_t k=0;k<meshes.size();++k){auto system=make_system<6>(meshes[k],theta,mode,reference);const double fe=norm(fits[k].control-fine.control);
                const double se=trajectory_error(system.forward(fits[k].control).state,continuous_fit_trajectory(fits[k].control),meshes[k]);state_errors[k]=se;
                if(theta==.5&&k+1==meshes.size())representative_trajectories[2]=system.forward(fits[k].control).state;
                convergence<<"joint,"<<theta<<','<<meshes[k]<<','<<(tf-t0)/meshes[k]<<','<<se<<','<<fe<<','<<tangent_fd_error(system,fits[k].control)<<','<<gradient_fd_error(system,fits[k].control)<<','<<hessian_fd_error(system,fits[k].control)<<'\n';
                const auto& p=fits[k].control;summary<<"joint,"<<theta<<','<<meshes[k]<<','<<p[0]<<','<<p[1]<<','<<p[2]<<','<<p[3]<<','<<p[4]<<','<<p[5]<<','<<fits[k].evaluation.cost<<','<<norm(fits[k].evaluation.gradient)<<','<<fits[k].iterations<<'\n';}
            auto system=make_system<6>(160,theta,mode,reference);
            checks.require("joint_tangent_fd_theta_"+std::to_string(theta),tangent_fd_error(system,fits[2].control),1e-6);
            checks.require("joint_gradient_fd_theta_"+std::to_string(theta),gradient_fd_error(system,fits[2].control),1e-4);
            checks.require("joint_hessian_fd_theta_"+std::to_string(theta),hessian_fd_error(system,fits[2].control),1e-4);
            const double state_order=order(state_errors[2],state_errors[3]);
            const double fit_order=order(norm(fits[1].control-fits[2].control),norm(fits[2].control-fits[3].control));
            checks.require("joint_state_order_theta_"+std::to_string(theta),std::abs(state_order-(theta==.5?2.:1.)),.25);
            checks.require("joint_fit_order_theta_"+std::to_string(theta),std::abs(fit_order-(theta==.5?2.:1.)),.4);
            std::cout<<"  joint theta="<<theta<<" observed orders: trajectory="<<state_order<<", joint estimate="<<fit_order<<"\n";
            if(theta==.5) {
                checks.require("joint_initial_recovery",norm(Vec<2>{fits.back().control[0]-truth_initial[0],fits.back().control[1]-truth_initial[1]}),.003);
                checks.require("joint_parameter_recovery",norm(Vec<4>{fits.back().control[2]-truth_rates[0],fits.back().control[3]-truth_rates[1],fits.back().control[4]-truth_rates[2],fits.back().control[5]-truth_rates[3]}),.003);
                const auto spec=information_spectrum(make_system<6>(meshes.back(),theta,mode,reference),fits.back().control,fits.back().evaluation.trajectory);
                identifiability<<"joint,"<<spec.rank<<",6,"<<spec.eigenvalues[spec.rank-1]<<','<<spec.eigenvalues[0]<<','<<spec.eigenvalues[0]/spec.eigenvalues[spec.rank-1]<<'\n';
                checks.require("joint_information_full_rank",static_cast<double>(6-spec.rank),0.);
                std::cout<<"    local information rank="<<spec.rank<<"/6, condition="<<spec.eigenvalues[0]/spec.eigenvalues[spec.rank-1]<<"\n";
                std::cout<<"    N=320 estimate (x0,y0,a,b,c,d)=("<<fits.back().control[0]<<", "<<fits.back().control[1]<<", "<<fits.back().control[2]<<", "<<fits.back().control[3]<<", "<<fits.back().control[4]<<", "<<fits.back().control[5]<<")\n";
            }
        }
    }
    std::ofstream trajectories(dir/"reconstruction.csv");
    trajectories<<"time,prey_reference,predator_reference,prey_initial_only,predator_initial_only,prey_parameters_only,predator_parameters_only,prey_joint,predator_joint\n";
    for(std::size_t n=0;n<=meshes.back();++n){const double time=tf*static_cast<double>(n)/static_cast<double>(meshes.back());const auto truth=interpolate(reference.dense,time);
        trajectories<<time<<','<<truth[0]<<','<<truth[1];
        for(const auto& trajectory:representative_trajectories){if(trajectory.size()!=meshes.back()+1)throw std::logic_error("Missing reconstruction trajectory");trajectories<<','<<trajectory[n][0]<<','<<trajectory[n][1];}
        trajectories<<'\n';}
}

} // namespace

int main(int argc,char** argv) {
    try {
        const std::filesystem::path output=argc>1?argv[1]:"results/predator_prey";
        const auto reference=make_reference();
        std::filesystem::create_directories(output);
        std::cout<<"Predator-prey inverse problem (Lotka-Volterra), synthetic bounded-noise assays.\n";
        double drift=0.;const double initial_invariant=invariant(reference.dense.front(),truth_rates);
        for(const auto& x:reference.dense)drift=std::max(drift,std::abs(invariant(x,truth_rates)-initial_invariant));
        std::cout<<"Reference invariant max drift: "<<std::setprecision(10)<<drift<<"\n";
        CheckLog checks(output/"checks.csv");checks.require("rk4_reference_invariant",drift,2e-12);
        write_outputs(output,reference,checks);
        std::cout<<"Checks passed: "<<checks.count<<"\n";
        std::cout<<"Results: \""<<output.string()<<"\"\n";
        return 0;
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
