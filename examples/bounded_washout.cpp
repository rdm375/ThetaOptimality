// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#include <theta_opt/theta_opt.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>

// A pharmacokinetic washout with a physically bounded initial concentration.
// Blank-corrected assay noise can favor a negative concentration near the
// detection limit, while a dose prepared at the stock limit can favor a value
// above that limit. The optimizer must return the constrained KKT solution.
namespace {
using namespace theta_opt;
inline constexpr double t0=2.,tf=8.,rate=.45,stock_limit=10.,baseline=.2;
inline constexpr double assay_sd=.02,assay_bound=.05,goal_tolerance=.02,theta=.5;
inline constexpr Mat<1,1> dynamics{{{-rate}}};
inline constexpr std::array<double,6> sample_time{0.,.5,1.5,3.,5.,6.};
struct Case {const char* name;double truth;std::array<double,6> noise;};
inline constexpr std::array<Case,3> cases{{
    {"near_detection_lower",.02,{-2.4,-2.0,-1.8,-1.5,-1.2,-1.0}},
    {"interior",4.0,{-.4,.1,-.2,.3,-.1,.2}},
    {"stock_limit_upper",10.,{1.2,1.0,.9,.8,.6,.5}}
}};
struct Interval {double lower{},upper{};};

double exact(double tau,double c0){return c0*std::exp(-rate*tau);}
double datum(const Case& c,std::size_t j){return baseline+exact(sample_time[j],c.truth)+assay_sd*c.noise[j];}
std::size_t node_at(const Grid& grid,double tau){for(std::size_t n=0;n<grid.time.size();++n)if(std::abs(grid.time[n]-t0-tau)<1e-12)return n;throw std::logic_error("assay time absent from grid");}
Grid make_grid(std::size_t n){return Grid::uniform(t0,tf,n,theta);}
System<1,1> make_system(std::size_t n,const Case& c){auto grid=make_grid(n);Problem<1,1> p;
    p.initial.value=[](const Vec<1>& u){return Vec<1>{stock_limit*u[0]};};
    p.initial.jacobian=[](const Vec<1>&){return Mat<1,1>{{{stock_limit}}};};
    p.initial.directional_jacobian=[](const Vec<1>&,const Vec<1>&){return Mat<1,1>{};};
    p.dynamics.value=[](double,const Vec<1>& x,const Vec<1>&){return Vec<1>{-rate*x[0]};};
    p.dynamics.derivatives=[](double,const Vec<1>&,const Vec<1>&){Derivatives<1,1> d;d.x=dynamics;return d;};
    p.dynamics.directional_derivatives=[](double,const Vec<1>&,const Vec<1>&,const Vec<1>&,const Vec<1>&){return Derivatives<1,1>{};};
    p.dynamics.admissible=[](double,const Vec<1>& x,const Vec<1>&){return x[0]>=-1e-12;};
    for(std::size_t j=0;j<sample_time.size();++j){const double y=datum(c,j);auto term=ScalarFunction<1,1>::zero();
        term.evaluate=[y](double,const Vec<1>& x,const Vec<1>&){ScalarValue<1,1> v;const double r=baseline+x[0]-y;v.value=.5*r*r/(assay_sd*assay_sd);v.x[0]=r/(assay_sd*assay_sd);return v;};
        term.directional_gradient=[](double,const Vec<1>&,const Vec<1>&,const Vec<1>& dx,const Vec<1>&){ScalarDirection<1,1> v;v.x[0]=dx[0]/(assay_sd*assay_sd);return v;};
        p.nodal.emplace_back(node_at(grid,sample_time[j]),term);
    }
    return System<1,1>(std::move(p),std::move(grid));
}
double exact_unconstrained_fit(const Case& c){double numerator=0.,denominator=0.;for(std::size_t j=0;j<sample_time.size();++j){double s=std::exp(-rate*sample_time[j]),y=datum(c,j)-baseline;numerator+=s*y;denominator+=s*s;}return numerator/denominator;}
double discrete_reference_control(const Case& c,std::size_t n){const double dt=(tf-t0)/static_cast<double>(n),r=(1-rate*(1-theta)*dt)/(1+rate*theta*dt);double numerator=0.,denominator=0.;for(std::size_t j=0;j<sample_time.size();++j){const auto steps=static_cast<std::size_t>(std::llround(sample_time[j]/dt));const double sensitivity=stock_limit*std::pow(r,static_cast<double>(steps));numerator+=sensitivity*(datum(c,j)-baseline);denominator+=sensitivity*sensitivity;}return std::clamp(numerator/denominator,0.,1.);}
double exact_cost(const Case& c,double c0){double value=0.;for(std::size_t j=0;j<sample_time.size();++j){double r=baseline+exact(sample_time[j],c0)-datum(c,j);value+=.5*r*r/(assay_sd*assay_sd);}return value;}
double exact_gradient(const Case& c,double u){double value=0.;for(std::size_t j=0;j<sample_time.size();++j){double s=stock_limit*std::exp(-rate*sample_time[j]);double r=baseline+u*s-datum(c,j);value+=r*s/(assay_sd*assay_sd);}return value;}
double exact_hessian(){double value=0.;for(double tau:sample_time){double s=stock_limit*std::exp(-rate*tau);value+=s*s/(assay_sd*assay_sd);}return value;}
std::vector<double> numerical_allowance(const System<1,1>& sys){const auto& g=sys.grid();auto trajectory=sys.forward({1.});std::vector<double> defect(g.steps()),error(g.time.size());
    for(std::size_t n=0;n<g.steps();++n){const double h=g.time[n+1]-g.time[n];defect[n]=linear_theta_quadratic_defect_bound(dynamics,trajectory.state[n],trajectory.state[n+1],h,g.theta[n],[](const Vec<1>& x){return std::abs(x[0]);});error[n+1]=propagate_defect_bound(error[n],defect[n],h,-rate);}return error;}
Interval feasible_interval(const Case& c,const Grid& g,const std::vector<double>& error){Interval result{0.,stock_limit};for(std::size_t j=0;j<sample_time.size();++j){const double tau=sample_time[j],s=std::exp(-rate*tau),y=datum(c,j)-baseline,r=assay_bound+error[node_at(g,tau)];result.lower=std::max(result.lower,(y-r)/s);result.upper=std::min(result.upper,(y+r)/s);}if(result.lower>result.upper)throw std::runtime_error("empty observation-consistent concentration interval");return result;}
struct Checks {std::ofstream out;std::size_t count{};explicit Checks(const std::filesystem::path& p):out(p){out<<"check,error,tolerance,passed\n";}void require(const std::string& name,double error,double tol){bool ok=std::isfinite(error)&&error<=tol;out<<name<<','<<std::setprecision(17)<<error<<','<<tol<<','<<(ok?1:0)<<'\n';++count;if(!ok)throw std::runtime_error("Check failed: "+name+" error="+std::to_string(error));}};
struct Run {BoxOptimum<1> optimum;Interval feasible;std::vector<double> allowance;double qualified_time{-1.};};
Run run_case(const Case& c,std::size_t n,Checks& checks,bool validate){auto sys=make_system(n,c);auto opt=minimize_box_constrained(sys,{.5},{0.},{1.},OptimizeOptions{2e-10,40,50});auto eval=sys.evaluate(opt.control);const double exact_unconstrained=exact_unconstrained_fit(c);
    const bool at_lower=opt.active_lower[0],at_upper=opt.active_upper[0];if(validate){const bool expected_lower=exact_unconstrained<0.,expected_upper=exact_unconstrained>stock_limit;
        checks.require(std::string(c.name)+" expected active set",(at_lower==expected_lower&&at_upper==expected_upper)?0.:1.,0.);
        checks.require(std::string(c.name)+" KKT residual",opt.kkt_residual,2e-9);
        double state_violation=0.;for(const auto& x:eval.trajectory.state)state_violation=std::max(state_violation,-x[0]);checks.require(std::string(c.name)+" nonnegative discrete trajectory",std::max(0.,state_violation),1e-12);
        checks.require(std::string(c.name)+" discrete fit matches analytical discrete constrained minimizer",std::abs(opt.control[0]-discrete_reference_control(c,n)),2e-10);
        if(expected_lower)checks.require(std::string(c.name)+" lower multiplier sign",opt.lower_multiplier[0]>0.&&std::abs(opt.gradient[0]-opt.lower_multiplier[0])<2e-9?0.:1.,0.);
        else if(expected_upper)checks.require(std::string(c.name)+" upper multiplier sign",opt.upper_multiplier[0]>0.&&std::abs(opt.gradient[0]+opt.upper_multiplier[0])<2e-9?0.:1.,0.);
        else checks.require(std::string(c.name)+" interior stationarity",std::abs(opt.gradient[0]),2e-8);
        const double u=.37,eps=1e-4;const double exact_grad_fd=(exact_cost(c,stock_limit*(u+eps))-exact_cost(c,stock_limit*(u-eps)))/(2*eps);const double exact_hess_fd=(exact_gradient(c,u+eps)-exact_gradient(c,u-eps))/(2*eps);checks.require(std::string(c.name)+" analytical exact-cost gradient",std::abs(exact_grad_fd-exact_gradient(c,u)),2e-6);checks.require(std::string(c.name)+" analytical exact-cost Hessian",std::abs(exact_hess_fd-exact_hessian()),2e-6);
    }
    auto allowance=numerical_allowance(sys);auto feasible=feasible_interval(c,sys.grid(),allowance);if(validate){checks.require(std::string(c.name)+" bounded-data set contains truth",std::max({0.,feasible.lower-c.truth,c.truth-feasible.upper}),1e-12);}
    std::vector<double> edges;std::vector<std::array<double,1>> bounds;for(std::size_t i=0;i<sys.grid().steps();++i){double a=sys.grid().time[i]-t0,b=sys.grid().time[i+1]-t0,h=b-a;edges.push_back(a);double spread=.5*(feasible.upper-feasible.lower)*std::exp(-rate*a);bounds.push_back({spread+std::max(allowance[i],allowance[i+1])});(void)h;}edges.push_back(tf-t0);auto q=qualify_intervals<1>(edges,bounds,{goal_tolerance});Run result{std::move(opt),feasible,std::move(allowance),q.persistent[0].value_or(-1.)};return result;
}
}

int main(int argc,char** argv){try{const std::filesystem::path out=argc>1?argv[1]:"results/bounded_washout";std::filesystem::create_directories(out);Checks checks(out/"checks.csv");
    std::ofstream summary(out/"summary.csv"),refinement(out/"refinement.csv"),curves(out/"trajectory_qualification.csv");summary<<std::setprecision(17)<<"case,unconstrained_exact_fit_mg_L,constrained_exact_fit_mg_L,discrete_fit_mg_L,raw_gradient,lower_multiplier,upper_multiplier,kkt_residual,active_bound,feasible_lower_mg_L,feasible_upper_mg_L,qualification_time_h\n";refinement<<std::setprecision(17)<<"case,steps,dt_h,fit_mg_L,analytic_fit_mg_L,control_error,maximum_state_error,maximum_tangent_error,objective\n";curves<<std::setprecision(17)<<"case,elapsed_h,truth_mg_L,discrete_fit_mg_L,exact_fit_mg_L,feasible_center_mg_L,uncertainty_bound_mg_L,tolerance_mg_L\n";
    for(const auto& c:cases){const double unconstrained=exact_unconstrained_fit(c),reference=std::clamp(unconstrained,0.,stock_limit);double previous_error=0.,previous_state_error=0.,previous_tangent_error=0.;Run finest;
        for(std::size_t n:{24u,48u,96u,192u}){auto result=run_case(c,n,checks,n==96);auto sys=make_system(n,c);double fit=stock_limit*result.optimum.control[0],control_error=std::abs(fit-reference);auto trajectory=sys.forward(result.optimum.control);auto fixed=sys.forward({reference/stock_limit});double state_error=0.,tangent_error=0.;auto tangent=sys.tangent({reference/stock_limit},fixed,{1.});
            for(std::size_t j=0;j<sys.grid().time.size();++j){double tau=sys.grid().time[j]-t0;state_error=std::max(state_error,std::abs(fixed.state[j][0]-exact(tau,reference)));tangent_error=std::max(tangent_error,std::abs(tangent[j][0]-stock_limit*std::exp(-rate*tau)));}
            refinement<<c.name<<','<<n<<','<<(tf-t0)/static_cast<double>(n)<<','<<fit<<','<<reference<<','<<control_error<<','<<state_error<<','<<tangent_error<<','<<sys.evaluate(result.optimum.control).cost<<'\n';
            if(n==96){if(previous_state_error>0.&&state_error>0.){double order=std::log(previous_state_error/state_error)/std::log(2.);checks.require(std::string(c.name)+" midpoint state convergence order",std::abs(order-2.),.12);}else checks.require(std::string(c.name)+" zero-state equilibrium remains exact",state_error,1e-16);double tangent_order=std::log(previous_tangent_error/tangent_error)/std::log(2.);checks.require(std::string(c.name)+" midpoint tangent convergence order",std::abs(tangent_order-2.),.12);if(c.truth==4.){double control_order=std::log(previous_error/control_error)/std::log(2.);checks.require(std::string(c.name)+" interior fitted-state convergence order",std::abs(control_order-2.),.2);} }
            if(n==192){if(previous_state_error>0.&&state_error>0.){double order=std::log(previous_state_error/state_error)/std::log(2.);checks.require(std::string(c.name)+" midpoint state order on finest pair",std::abs(order-2.),.15);}else checks.require(std::string(c.name)+" zero-state equilibrium remains exact on finest grid",state_error,1e-16);double tangent_order=std::log(previous_tangent_error/tangent_error)/std::log(2.);checks.require(std::string(c.name)+" midpoint tangent order on finest pair",std::abs(tangent_order-2.),.15);}previous_error=control_error;previous_state_error=state_error;previous_tangent_error=tangent_error;if(n==192)finest=std::move(result);
        }
        const char* bound=finest.optimum.active_lower[0]?"lower":(finest.optimum.active_upper[0]?"upper":"interior");summary<<c.name<<','<<unconstrained<<','<<reference<<','<<stock_limit*finest.optimum.control[0]<<','<<finest.optimum.gradient[0]<<','<<finest.optimum.lower_multiplier[0]<<','<<finest.optimum.upper_multiplier[0]<<','<<finest.optimum.kkt_residual<<','<<bound<<','<<finest.feasible.lower<<','<<finest.feasible.upper<<','<<finest.qualified_time<<'\n';
        auto sys=make_system(192,c);auto trajectory=sys.forward(finest.optimum.control);const double center=.5*(finest.feasible.lower+finest.feasible.upper);auto center_path=sys.forward({center/stock_limit});for(std::size_t j=0;j<sys.grid().time.size();++j){double tau=sys.grid().time[j]-t0;double bound_width=.5*(finest.feasible.upper-finest.feasible.lower)*std::exp(-rate*tau)+finest.allowance[j];curves<<c.name<<','<<tau<<','<<exact(tau,c.truth)<<','<<trajectory.state[j][0]<<','<<exact(tau,reference)<<','<<center_path.state[j][0]<<','<<bound_width<<','<<goal_tolerance<<'\n';}
        if(c.truth==4.){Vec<1> d{.37},u{.41};auto ev=sys.evaluate(u);auto tv=sys.tangent(u,ev.trajectory,d);double pairing=dot(ev.gradient,d),sum=0.;for(std::size_t j=0;j<sample_time.size();++j){auto m=node_at(sys.grid(),sample_time[j]);const double r=ev.trajectory.state[m][0]+baseline-datum(c,j);sum+=r*tv[m][0]/(assay_sd*assay_sd);}checks.require("tangent-adjoint scalar pairing",std::abs(pairing-sum),2e-9);double eps=1e-5;auto ep=sys.evaluate({u[0]+eps*d[0]}),em=sys.evaluate({u[0]-eps*d[0]});checks.require("discrete adjoint objective gradient",std::abs((ep.cost-em.cost)/(2*eps)-pairing),2e-6);checks.require("incremental-adjoint Hessian",std::abs(sys.hessian_vector(u,d)[0]-(ep.gradient[0]-em.gradient[0])/(2*eps)),2e-6);checks.require("theta discrete objective versus exact discrete reference",std::abs(ev.cost-([&]{double jv=0.;for(std::size_t j=0;j<sample_time.size();++j){double h=(1-rate*(1-theta)*(tf-t0)/192)/(1+rate*theta*(tf-t0)/192);double s=stock_limit*std::pow(h,sample_time[j]/((tf-t0)/192));double r=s*u[0]+baseline-datum(c,j);jv+=.5*r*r/(assay_sd*assay_sd);}return jv;})()),2e-8);}
    }
    std::cout<<"BOUNDED INITIAL-STATE WASHOUT CASE\n"<<"Known one-compartment elimination rate: "<<rate<<" h^-1; initial concentration is constrained to [0,"<<stock_limit<<"] mg/L.\n"<<"Cases exercise interior, lower-bound, and upper-bound optima. KKT multipliers and numerical error are included.\n";
    for(const auto& c:cases){auto r=run_case(c,192,checks,false);std::cout<<"  "<<c.name<<": c0="<<stock_limit*r.optimum.control[0]<<" mg/L, active="<<(r.optimum.active_lower[0]?"lower":r.optimum.active_upper[0]?"upper":"none")<<", gradient="<<r.optimum.gradient[0]<<", multipliers (lower,upper)=("<<r.optimum.lower_multiplier[0]<<','<<r.optimum.upper_multiplier[0]<<"), KKT="<<r.optimum.kkt_residual<<", qualification="<<(r.qualified_time<0?"NOT QUALIFIED":std::to_string(r.qualified_time)+" h")<<'\n';}
    std::cout<<"Checks passed: "<<checks.count<<". Results: \""<<out.string()<<"\"\n";return 0;
}catch(const std::exception& e){std::cerr<<"bounded_washout: "<<e.what()<<'\n';return 1;}}
