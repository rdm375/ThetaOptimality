// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#include <theta_opt/theta_opt.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

// Batch-culture biomass reconstruction with known logistic-growth parameters.
// This isolates state-dependent first and second derivatives while retaining
// a closed-form continuous state, tangent, objective gradient, and Hessian.
namespace {
using namespace theta_opt;
inline constexpr double t0=0.,tf=4.,growth_rate=.8,carrying_capacity=8.;
inline constexpr double truth_xi=1.6,measurement_sd=.08;
inline constexpr double assay_bound=.04,goal_tolerance=.02;
inline constexpr double control_lower=.02,control_upper=.60;
inline constexpr std::array<double,5> sample_times{.25,.75,1.5,2.5,3.75};
inline constexpr std::array<double,5> standardized_noise{-.30,.25,-.20,.10,.15};

double exact_state(double tau,double xi) {
    const double e=std::exp(growth_rate*tau);
    return carrying_capacity*xi*e/(carrying_capacity+xi*std::expm1(growth_rate*tau));
}
double exact_sensitivity(double tau,double xi) {
    const double e=std::exp(growth_rate*tau);
    const double d=carrying_capacity+xi*std::expm1(growth_rate*tau);
    return carrying_capacity*carrying_capacity*e/(d*d);
}
double exact_second_sensitivity(double tau,double xi) {
    const double e=std::exp(growth_rate*tau);
    const double em1=std::expm1(growth_rate*tau);
    const double d=carrying_capacity+xi*em1;
    return -2.*carrying_capacity*carrying_capacity*e*em1/(d*d*d);
}
double observation(std::size_t j) {
    return exact_state(sample_times.at(j),truth_xi)+measurement_sd*standardized_noise.at(j);
}
std::size_t observation_node(const Grid& grid,double tau) {
    for(std::size_t n=0;n<grid.time.size();++n)
        if(std::abs(grid.time[n]-t0-tau)<1e-12)return n;
    throw std::logic_error("Observation time is not a grid node");
}
double inverse_state(double tau,double state) {
    return carrying_capacity*state/
        (state+(carrying_capacity-state)*std::exp(growth_rate*tau));
}

Problem<1,1> make_problem(const Grid& grid) {
    Problem<1,1> problem;
    problem.initial.value=[](const Vec<1>& u){return Vec<1>{carrying_capacity*u[0]};};
    problem.initial.jacobian=[](const Vec<1>&){return Mat<1,1>{{{carrying_capacity}}};};
    problem.initial.directional_jacobian=[](const Vec<1>&,const Vec<1>&){return Mat<1,1>{};};
    problem.dynamics.value=[](double,const Vec<1>& x,const Vec<1>&){
        return Vec<1>{growth_rate*x[0]*(1.-x[0]/carrying_capacity)};
    };
    problem.dynamics.derivatives=[](double,const Vec<1>& x,const Vec<1>&){
        Derivatives<1,1> d;d.x[0][0]=growth_rate*(1.-2.*x[0]/carrying_capacity);return d;
    };
    problem.dynamics.directional_derivatives=[](double,const Vec<1>&,const Vec<1>&,
                                                  const Vec<1>& dx,const Vec<1>&){
        Derivatives<1,1> d;d.x[0][0]=-2.*growth_rate*dx[0]/carrying_capacity;return d;
    };
    problem.dynamics.admissible=[](double,const Vec<1>& x,const Vec<1>&){
        return x[0]>0.&&x[0]<carrying_capacity;
    };
    for(std::size_t j=0;j<sample_times.size();++j) {
        const double y=observation(j);
        auto term=ScalarFunction<1,1>::zero();
        term.evaluate=[y](double,const Vec<1>& x,const Vec<1>&){
            ScalarValue<1,1> value;const double residual=(x[0]-y)/measurement_sd;
            value.value=.5*residual*residual;
            value.x[0]=(x[0]-y)/(measurement_sd*measurement_sd);return value;
        };
        term.directional_gradient=[](double,const Vec<1>&,const Vec<1>&,
                                      const Vec<1>& dx,const Vec<1>&){
            ScalarDirection<1,1> value;value.x[0]=dx[0]/(measurement_sd*measurement_sd);return value;
        };
        problem.nodal.emplace_back(observation_node(grid,sample_times[j]),std::move(term));
    }
    return problem;
}
System<1,1> make_system(std::size_t steps,double theta) {
    auto grid=Grid::uniform(t0,tf,steps,theta);
    auto problem=make_problem(grid);
    return System<1,1>(std::move(problem),std::move(grid));
}
double exact_cost(double u) {
    const double xi=carrying_capacity*u;double cost=0.;
    for(std::size_t j=0;j<sample_times.size();++j) {
        const double residual=exact_state(sample_times[j],xi)-observation(j);
        cost+=.5*residual*residual/(measurement_sd*measurement_sd);
    }
    return cost;
}
double exact_gradient(double u) {
    const double xi=carrying_capacity*u;double gradient=0.;
    for(std::size_t j=0;j<sample_times.size();++j) {
        const double residual=exact_state(sample_times[j],xi)-observation(j);
        const double sensitivity=carrying_capacity*exact_sensitivity(sample_times[j],xi);
        gradient+=residual*sensitivity/(measurement_sd*measurement_sd);
    }
    return gradient;
}
double exact_hessian(double u) {
    const double xi=carrying_capacity*u;double hessian=0.;
    for(std::size_t j=0;j<sample_times.size();++j) {
        const double residual=exact_state(sample_times[j],xi)-observation(j);
        const double first=carrying_capacity*exact_sensitivity(sample_times[j],xi);
        const double second=carrying_capacity*carrying_capacity*
                            exact_second_sensitivity(sample_times[j],xi);
        hessian+=(first*first+residual*second)/(measurement_sd*measurement_sd);
    }
    return hessian;
}
double exact_fit() {
    double lo=control_lower,hi=control_upper;
    if(!(exact_gradient(lo)<0.&&exact_gradient(hi)>0.))
        throw std::runtime_error("Analytical optimum is not bracketed");
    for(int k=0;k<100;++k) {
        const double mid=.5*(lo+hi);
        if(exact_gradient(mid)>0.)hi=mid;else lo=mid;
    }
    return .5*(lo+hi);
}
double exact_costate(double tau,double xi) {
    double lambda=0.;
    for(std::size_t j=0;j<sample_times.size();++j)if(sample_times[j]>tau+1e-12) {
        const double residual=(exact_state(sample_times[j],xi)-observation(j))/
                              (measurement_sd*measurement_sd);
        lambda+=residual*exact_sensitivity(sample_times[j],xi)/exact_sensitivity(tau,xi);
    }
    return lambda;
}
struct CheckLog {
    std::ofstream csv;std::size_t count{};
    explicit CheckLog(const std::filesystem::path& path):csv(path) {
        if(!csv)throw std::runtime_error("Cannot open checks CSV");
        csv<<"check,error,tolerance,passed\n";
    }
    void require(const std::string& name,double error,double tolerance) {
        const bool passed=std::isfinite(error)&&error<=tolerance;
        csv<<name<<','<<std::setprecision(17)<<error<<','<<tolerance<<','<<(passed?1:0)<<'\n';
        ++count;if(!passed)throw std::runtime_error("Check failed: "+name+" error="+std::to_string(error));
    }
};
struct MeshRow {double theta;std::size_t steps;double h;std::array<double,7> error;double fit;};
std::optional<double> observed_order(double coarse,double fine) {
    if(!(coarse>0.)||!(fine>0.)||!std::isfinite(coarse)||!std::isfinite(fine))return std::nullopt;
    return std::log(coarse/fine)/std::log(2.);
}
std::pair<double,double> feasible_initial_interval() {
    double lower=.4,upper=4.5;
    for(std::size_t j=0;j<sample_times.size();++j) {
        const double y=observation(j),lo_state=std::max(0.,y-assay_bound);
        const double hi_state=std::min(carrying_capacity,y+assay_bound);
        lower=std::max(lower,inverse_state(sample_times[j],lo_state));
        upper=std::min(upper,inverse_state(sample_times[j],hi_state));
    }
    if(lower>upper)throw std::runtime_error("No initial state fits the bounded observations");
    return {lower,upper};
}
double pair_difference(double tau,double a,double b) {
    return std::abs(exact_state(tau,a)-exact_state(tau,b));
}
double max_pair_difference(double a_time,double b_time,double a,double b) {
    if(std::abs(a-b)<1e-15)return 0.;
    const double sum_a=exact_state(a_time,a)+exact_state(a_time,b);
    const double sum_b=exact_state(b_time,a)+exact_state(b_time,b);
    if(sum_a>=carrying_capacity)return pair_difference(a_time,a,b);
    if(sum_b<=carrying_capacity)return pair_difference(b_time,a,b);
    double lo=a_time,hi=b_time;
    for(int k=0;k<80;++k) {
        const double mid=.5*(lo+hi);
        if(exact_state(mid,a)+exact_state(mid,b)<carrying_capacity)lo=mid;else hi=mid;
    }
    return pair_difference(.5*(lo+hi),a,b);
}
std::vector<double> theta_error_bound(const System<1,1>& system,const Vec<1>& control,
                                      const Trajectory<1>& trajectory) {
    const auto& grid=system.grid();std::vector<double> error(grid.time.size(),0.);
    for(std::size_t n=0;n<grid.steps();++n) {
        const double h=grid.time[n+1]-grid.time[n];
        const double exact_step=exact_state(h,trajectory.state[n][0]);
        const double defect=std::abs(exact_step-trajectory.state[n+1][0]);
        // The exact logistic flow map has derivative <= exp(r h) on x >= 0.
        error[n+1]=std::exp(growth_rate*h)*error[n]+defect;
    }
    (void)control;return error;
}
struct Qualification {double fit_xi{},lower{},upper{},numerical_allowance{};std::optional<double> time;};
Qualification qualify_case(std::size_t steps,std::ofstream& curves,std::ofstream& intervals) {
    auto system=make_system(steps,.5);
    const auto optimum=minimize(system,{.2},OptimizeOptions{1e-11,40,50},
        [](const Vec<1>& u){return u[0]>control_lower&&u[0]<control_upper;});
    const auto trajectory=system.forward(optimum.control);const double fit=carrying_capacity*optimum.control[0];
    const auto [lo,hi]=feasible_initial_interval();
    const auto errors=theta_error_bound(system,optimum.control,trajectory);
    const double model_allowance=*std::max_element(errors.begin(),errors.end());
    const auto& grid=system.grid();std::vector<double> edges;std::vector<std::array<double,1>> bounds;
    for(std::size_t n=0;n<grid.steps();++n) {
        const double a=grid.time[n]-t0,b=grid.time[n+1]-t0;
        const double spread=std::max(max_pair_difference(a,b,lo,fit),max_pair_difference(a,b,hi,fit));
        bounds.push_back({spread+model_allowance});edges.push_back(grid.time[n]-t0);
        intervals<<steps<<','<<a<<','<<b<<','<<spread<<','<<model_allowance<<','<<spread+model_allowance<<'\n';
    }
    edges.push_back(tf-t0);
    const auto result=qualify_intervals<1>(edges,bounds,{goal_tolerance});
    for(std::size_t n=0;n<grid.time.size();++n) {
        const double tau=grid.time[n]-t0;
        curves<<steps<<','<<tau<<','<<exact_state(tau,truth_xi)<<','<<trajectory.state[n][0]<<','
              <<exact_state(tau,fit)<<','<<exact_state(tau,lo)-model_allowance<<','
              <<exact_state(tau,hi)+model_allowance<<','<<model_allowance<<'\n';
    }
    return {fit,lo,hi,model_allowance,result.persistent[0]};
}
}

int main(int argc,char** argv) {
    try {
        using namespace theta_opt;namespace fs=std::filesystem;
        const fs::path output=argc>1?argv[1]:"results/logistic_growth";fs::create_directories(output);
        CheckLog checks(output/"checks.csv");
        std::ofstream convergence(output/"convergence.csv"),curves(output/"qualification_curves.csv"),
                      intervals(output/"qualification_intervals.csv"),summary(output/"summary.csv"),
                      optimality(output/"optimality_curves.csv"),observations(output/"observations.csv");
        if(!convergence||!curves||!intervals||!summary||!optimality||!observations)
            throw std::runtime_error("Cannot open output CSV");
        convergence<<"theta,steps,h,state,tangent,adjoint,cost,gradient,hessian,recovered_initial\n";
        curves<<"steps,elapsed_h,truth_g_L,theta_fit_g_L,analytic_fit_g_L,feasible_lower_g_L,feasible_upper_g_L,theta_error_allowance_g_L\n";
        intervals<<"steps,interval_start_h,interval_end_h,feasible_spread_g_L,theta_allowance_g_L,total_bound_g_L\n";
        summary<<"steps,discrete_fit_g_L,analytical_fit_g_L,feasible_lower_g_L,feasible_upper_g_L,theta_allowance_g_L,qualification_time_h\n";
        optimality<<"control_u,exact_cost,exact_gradient,exact_hessian,cost_N32,gradient_N32,hessian_N32,cost_N128,gradient_N128,hessian_N128\n";
        observations<<"elapsed_h,measured_biomass_g_L,assay_bound_g_L\n";
        for(std::size_t j=0;j<sample_times.size();++j)
            observations<<sample_times[j]<<','<<observation(j)<<','<<assay_bound<<'\n';

        const double fit_reference=exact_fit(),probe=.2,direction=.37;
        const double eps=2e-6;
        const double fd_gradient=(exact_cost(probe+eps)-exact_cost(probe-eps))/(2.*eps);
        const double fd_hessian=(exact_gradient(probe+eps)-exact_gradient(probe-eps))/(2.*eps);
        checks.require("analytical objective gradient finite difference",std::abs(fd_gradient-exact_gradient(probe)),2e-6);
        checks.require("analytical objective Hessian finite difference",std::abs(fd_hessian-exact_hessian(probe)),2e-5);
        checks.require("bounded observations contain the truth",std::max({0.,feasible_initial_interval().first-truth_xi,
            truth_xi-feasible_initial_interval().second}),1e-12);

        std::vector<MeshRow> rows;
        for(double theta:{.5,.65})for(std::size_t steps:{32u,64u,128u,256u,512u,1024u}) {
            auto system=make_system(steps,theta);const Vec<1> p{probe};
            const auto trajectory=system.forward(p);const auto tangent=system.tangent(p,trajectory,{direction});
            const auto evaluation=system.adjoint(p,trajectory);
            const auto costate=system.cost_to_go_adjoint(p,evaluation);
            const auto hessian=system.hessian_vector(p,{1.});
            Optimum<1> optimum;
            try {
                optimum=minimize(system,{.2},OptimizeOptions{1e-10,40,50,1e-13},
                    [](const Vec<1>& u){return u[0]>control_lower&&u[0]<control_upper;});
            } catch(const std::exception& error) {
                throw std::runtime_error("logistic fit theta="+std::to_string(theta)+
                    " N="+std::to_string(steps)+": "+error.what());
            }
            double state_error=0.,tangent_error=0.,adjoint_error=0.;
            double pairing=0.;
            for(std::size_t j=0;j<sample_times.size();++j) {
                const auto n=observation_node(system.grid(),sample_times[j]);
                const double residual=(trajectory.state[n][0]-observation(j))/(measurement_sd*measurement_sd);
                pairing+=residual*tangent[n][0];
            }
            for(std::size_t n=0;n<=steps;++n) {
                const double tau=system.grid().time[n]-t0;
                state_error=std::max(state_error,std::abs(trajectory.state[n][0]-exact_state(tau,carrying_capacity*probe)));
                tangent_error=std::max(tangent_error,std::abs(tangent[n][0]-direction*carrying_capacity*exact_sensitivity(tau,carrying_capacity*probe)));
                bool is_observation=false;for(double sample:sample_times)if(std::abs(sample-tau)<1e-12)is_observation=true;
                if(!is_observation)adjoint_error=std::max(adjoint_error,
                    std::abs(costate[n][0]-exact_costate(tau,carrying_capacity*probe)));
            }
            const double objective_error=std::abs(evaluation.cost-exact_cost(probe));
            const double gradient_error=std::abs(evaluation.gradient[0]-exact_gradient(probe));
            const double hessian_error=std::abs(hessian[0]-exact_hessian(probe));
            const double fit_error=std::abs(carrying_capacity*optimum.control[0]-carrying_capacity*fit_reference);
            const std::array<double,7> errors{state_error,tangent_error,adjoint_error,objective_error,
                                               gradient_error,hessian_error,fit_error};
            const double h=(tf-t0)/static_cast<double>(steps);
            convergence<<std::setprecision(17)<<theta<<','<<steps<<','<<h;
            for(double error:errors)convergence<<','<<error;
            convergence<<'\n';rows.push_back({theta,steps,h,errors,carrying_capacity*optimum.control[0]});
            if(steps==64) {
                checks.require("tangent-adjoint directional pairing theta="+std::to_string(theta),
                    std::abs(pairing-evaluation.gradient[0]*direction),2e-9);
                checks.require("discrete fit stationarity theta="+std::to_string(theta),
                    std::abs(system.evaluate(optimum.control).gradient[0]),2e-8);
                const double fd_discrete_gradient=(system.evaluate({probe+eps}).cost-
                    system.evaluate({probe-eps}).cost)/(2.*eps);
                const double fd_discrete_hessian=(system.evaluate({probe+eps}).gradient[0]-
                    system.evaluate({probe-eps}).gradient[0])/(2.*eps);
                checks.require("discrete objective gradient finite difference theta="+
                    std::to_string(theta),std::abs(fd_discrete_gradient-evaluation.gradient[0]),2e-5);
                checks.require("incremental-adjoint Hessian finite difference theta="+
                    std::to_string(theta),std::abs(fd_discrete_hessian-hessian[0]),2e-4);
            }
        }
        for(double theta:{.5,.65}) {
            std::array<const MeshRow*,6> group{};std::size_t count=0;
            for(const auto& row:rows)if(row.theta==theta)group.at(count++)=&row;
            if(count!=group.size())throw std::logic_error("Unexpected refinement count");
            constexpr std::array<const char*,7> quantity{"state","tangent","adjoint","cost","gradient","Hessian","fitted initial state"};
            std::cout<<"theta="<<theta<<" observed finest-pair orders:";
            for(std::size_t q=0;q<7;++q) {
                const auto order=observed_order(group[4]->error[q],group[5]->error[q]);
                const double expected=theta==.5?2.:1.;
                std::cout<<' '<<quantity[q]<<'='<<(order?*order:std::numeric_limits<double>::quiet_NaN());
                checks.require("theta="+std::to_string(theta)+" order for "+quantity[q],
                    order?std::abs(*order-expected):std::numeric_limits<double>::infinity(),
                    theta==.5?.30:.35);
            }
            std::cout<<"\n";
        }
        const auto system32=make_system(32,.5),system128=make_system(128,.5);
        for(std::size_t j=0;j<=100;++j) {
            const double u=control_lower+(control_upper-control_lower)*static_cast<double>(j)/100.;
            const auto e32=system32.evaluate({u}),e128=system128.evaluate({u});
            optimality<<std::setprecision(17)<<u<<','<<exact_cost(u)<<','<<exact_gradient(u)<<','
                      <<exact_hessian(u)<<','<<e32.cost<<','<<e32.gradient[0]<<','
                      <<system32.hessian_vector({u},{1.})[0]<<','<<e128.cost<<','
                      <<e128.gradient[0]<<','<<system128.hessian_vector({u},{1.})[0]<<'\n';
        }
        for(std::size_t steps:{32u,128u}) {
            auto q=qualify_case(steps,curves,intervals);
            summary<<steps<<','<<q.fit_xi<<','<<carrying_capacity*fit_reference<<','<<q.lower<<','<<q.upper<<','
                   <<q.numerical_allowance<<',';
            if(q.time)summary<<*q.time;else summary<<"unqualified";
            summary<<'\n';
            std::cout<<"N="<<steps<<": fitted initial biomass="<<q.fit_xi<<" g/L; exact bounded-data interval=["
                     <<q.lower<<", "<<q.upper<<"] g/L; theta allowance="<<q.numerical_allowance
                     <<" g/L; qualification="<<(q.time?std::to_string(*q.time)+" h":"NOT QUALIFIED")<<'\n';
            checks.require("bounded interval contains the fitted analytical minimizer",
                std::max({0.,q.lower-carrying_capacity*fit_reference,
                           carrying_capacity*fit_reference-q.upper}),1e-10);
        }
        std::cout<<"Analytical fit: "<<carrying_capacity*fit_reference<<" g/L; checks passed: "
                 <<checks.count<<". Results: "<<output<<'\n';
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
