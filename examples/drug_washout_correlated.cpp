#include <theta_opt/theta_opt.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace theta_opt;

constexpr double k_elimination = 0.45; // 1/hour
constexpr double final_time = 8.0; // hours
constexpr double true_initial = 10.0; // mg/L
constexpr double prior_lower = 0.0;
constexpr double prior_upper = 20.0;
constexpr double sensor_sd = 0.020; // mg/L, marginal standard deviation
constexpr double correlation = 0.70; // AR(1), one-sample lag
constexpr double chi_square_6_95 = 12.591587243743977;
constexpr double concentration_tolerance = 0.020; // mg/L
constexpr double auc_tolerance = 0.060; // mg hour/L
constexpr std::array<double, 6> sample_times{0.5, 1.0, 2.0, 4.0, 6.0, 8.0};
constexpr std::array<double, 6> innovations{0.7, -0.4, 0.6, -0.2, -0.5, 0.3};
constexpr std::array<std::size_t, 5> convergence_steps{32, 64, 128, 256, 512};
constexpr std::array<std::size_t, 3> qualification_steps{32, 128, 512};

struct Checks {
    std::ofstream file;
    std::size_t count{};
    explicit Checks(const std::filesystem::path& path) : file(path) {
        file << std::setprecision(17) << "check,error,tolerance,passed\n";
        if (!file) throw std::runtime_error("Cannot open checks.csv");
    }
    void check(const std::string& name, double error, double tolerance) {
        const bool passed = std::isfinite(error) && error <= tolerance;
        file << name << ',' << error << ',' << tolerance << ',' << passed << '\n';
        ++count;
        if (!passed) throw std::runtime_error(name + " failed; error=" + std::to_string(error));
    }
};

struct Measurement { std::array<double, sample_times.size()> value{}; };
struct ReferenceFit {
    double estimate{};
    double lower{};
    double upper{};
    double information{};
    double minimum_quadratic{};
    double true_quadratic{};
};
struct ErrorComponents { double state{}, tangent{}, cost{}, gradient{}, hessian{}, fit{}, auc{}; };
struct Budget {
    std::vector<double> measurement_c, numerical_c, total_c;
    std::vector<double> measurement_auc, numerical_auc, total_auc;
    double measurement_c_time{}, numerical_c_time{}, total_c_time{};
    double measurement_auc_time{}, numerical_auc_time{}, total_auc_time{};
};

double concentration(double time, double initial) {
    return initial * std::exp(-k_elimination * time);
}
double remaining_auc(double time, double initial) {
    return initial * (std::exp(-k_elimination * time) -
                      std::exp(-k_elimination * final_time)) / k_elimination;
}

Measurement make_observations() {
    Measurement data;
    double previous_noise = sensor_sd * innovations[0];
    data.value[0] = concentration(sample_times[0], true_initial) + previous_noise;
    const double innovation_scale = sensor_sd * std::sqrt(1.0 - correlation * correlation);
    for (std::size_t j = 1; j < sample_times.size(); ++j) {
        previous_noise = correlation * previous_noise + innovation_scale * innovations[j];
        data.value[j] = concentration(sample_times[j], true_initial) + previous_noise;
    }
    return data;
}

std::vector<std::vector<double>> precision_matrix(bool use_correlation) {
    const std::size_t m = sample_times.size();
    std::vector<std::vector<double>> W(m, std::vector<double>(m, 0.0));
    if (!use_correlation) {
        for (std::size_t i = 0; i < m; ++i) W[i][i] = 1.0 / (sensor_sd * sensor_sd);
        return W;
    }
    const double scale = 1.0 / (sensor_sd * sensor_sd * (1.0 - correlation * correlation));
    W[0][0] = scale;
    W[m - 1][m - 1] = scale;
    for (std::size_t i = 1; i + 1 < m; ++i) W[i][i] = scale * (1.0 + correlation * correlation);
    for (std::size_t i = 0; i + 1 < m; ++i) {
        W[i][i + 1] = -scale * correlation;
        W[i + 1][i] = -scale * correlation;
    }
    return W;
}

double quadratic(const std::array<double, 6>& residual,
                 const std::vector<std::vector<double>>& W) {
    double value = 0.0;
    for (std::size_t i = 0; i < residual.size(); ++i)
        for (std::size_t j = 0; j < residual.size(); ++j)
            value += residual[i] * W[i][j] * residual[j];
    return value;
}

ReferenceFit analytical_fit(const Measurement& data,
                            const std::vector<std::vector<double>>& W) {
    std::array<double, 6> sensitivity{}, y{};
    for (std::size_t j = 0; j < sample_times.size(); ++j) {
        sensitivity[j] = std::exp(-k_elimination * sample_times[j]);
        y[j] = data.value[j];
    }
    double information = 0.0, rhs = 0.0;
    for (std::size_t i = 0; i < y.size(); ++i)
        for (std::size_t j = 0; j < y.size(); ++j) {
            information += sensitivity[i] * W[i][j] * sensitivity[j];
            rhs += sensitivity[i] * W[i][j] * y[j];
        }
    const double estimate = rhs / information;
    std::array<double, 6> residual{};
    for (std::size_t i = 0; i < y.size(); ++i) {
        residual[i] = sensitivity[i] * estimate - y[i];
    }
    const double minimum = quadratic(residual, W);
    const double discriminant = chi_square_6_95 - minimum;
    if (discriminant < 0.0) throw std::runtime_error("95% joint data ellipsoid excludes its fitted point");
    const double radius = std::sqrt(discriminant / information);
    ReferenceFit result{estimate, std::max(prior_lower, estimate - radius),
                        std::min(prior_upper, estimate + radius), information,
                        minimum, 0.0};
    std::array<double, 6> true_residual{};
    for (std::size_t j = 0; j < y.size(); ++j)
        true_residual[j] = sensitivity[j] * true_initial - y[j];
    result.true_quadratic = quadratic(true_residual, W);
    return result;
}

std::size_t sample_node(const Grid& grid, double time) {
    for (std::size_t n = 0; n < grid.time.size(); ++n)
        if (std::abs(grid.time[n] - time) < 1e-12) return n;
    throw std::runtime_error("Observation time is missing from grid");
}

Problem<1, 1> washout_problem(const Grid& grid, const Measurement& data,
                              const std::vector<std::vector<double>>& W) {
    Problem<1, 1> problem;
    problem.initial.value = [](const Vec<1>& c0) { return c0; };
    problem.initial.jacobian = [](const Vec<1>&) { return identity<1>(); };
    problem.initial.directional_jacobian = [](const Vec<1>&, const Vec<1>&) { return Mat<1, 1>{}; };
    problem.dynamics.value = [](double, const Vec<1>& x, const Vec<1>&) {
        return Vec<1>{-k_elimination * x[0]};
    };
    problem.dynamics.derivatives = [](double, const Vec<1>&, const Vec<1>&) {
        Derivatives<1, 1> d; d.x[0][0] = -k_elimination; return d;
    };
    problem.dynamics.directional_derivatives = [](double, const Vec<1>&, const Vec<1>&,
                                                   const Vec<1>&, const Vec<1>&) {
        return Derivatives<1, 1>{};
    };
    CorrelatedObservationTerm<1, 1> observations;
    observations.precision = W;
    for (std::size_t j = 0; j < sample_times.size(); ++j) {
        observations.node.push_back(sample_node(grid, sample_times[j]));
        const double observed = data.value[j];
        auto residual = ScalarFunction<1, 1>::zero();
        residual.evaluate = [observed](double, const Vec<1>& x, const Vec<1>&) {
            ScalarValue<1, 1> v; v.value = x[0] - observed; v.x[0] = 1.0; return v;
        };
        residual.directional_gradient = [](double, const Vec<1>&, const Vec<1>&,
                                           const Vec<1>&, const Vec<1>&) {
            return ScalarDirection<1, 1>{};
        };
        observations.residual.push_back(std::move(residual));
    }
    problem.correlated_observations.push_back(std::move(observations));
    return problem;
}

Problem<1, 1> exposure_problem() {
    Problem<1, 1> problem;
    problem.initial.value = [](const Vec<1>& c0) { return c0; };
    problem.initial.jacobian = [](const Vec<1>&) { return identity<1>(); };
    problem.initial.directional_jacobian = [](const Vec<1>&, const Vec<1>&) { return Mat<1, 1>{}; };
    problem.dynamics.value = [](double, const Vec<1>& x, const Vec<1>&) {
        return Vec<1>{-k_elimination * x[0]};
    };
    problem.dynamics.derivatives = [](double, const Vec<1>&, const Vec<1>&) {
        Derivatives<1, 1> d; d.x[0][0] = -k_elimination; return d;
    };
    problem.dynamics.directional_derivatives = [](double, const Vec<1>&, const Vec<1>&,
                                                   const Vec<1>&, const Vec<1>&) {
        return Derivatives<1, 1>{};
    };
    problem.running.evaluate = [](double, const Vec<1>& x, const Vec<1>&) {
        ScalarValue<1, 1> v; v.value = x[0]; v.x[0] = 1.0; return v;
    };
    problem.running.directional_gradient = [](double, const Vec<1>&, const Vec<1>&,
                                               const Vec<1>&, const Vec<1>&) {
        return ScalarDirection<1, 1>{};
    };
    return problem;
}

double discrete_closed_fit(const System<1, 1>& system, const Grid& grid,
                           const Measurement& data,
                           const std::vector<std::vector<double>>& W) {
    const auto unit = system.forward({1.0});
    std::array<double, 6> s{}, y{};
    for (std::size_t j = 0; j < sample_times.size(); ++j) {
        s[j] = unit.state[sample_node(grid, sample_times[j])][0];
        y[j] = data.value[j];
    }
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t i = 0; i < y.size(); ++i)
        for (std::size_t j = 0; j < y.size(); ++j) {
            numerator += s[i] * W[i][j] * y[j];
            denominator += s[i] * W[i][j] * s[j];
        }
    return numerator / denominator;
}

double exp_line_error(double left, double h, double line_left, double line_right,
                      double amplitude, double offset = 0.0) {
    const double slope = (line_right - line_left) / h;
    auto error = [&](double t) {
        const double alpha = (t - left) / h;
        const double line = (1.0 - alpha) * line_left + alpha * line_right;
        return std::abs(line - amplitude * std::exp(-k_elimination * t) - offset);
    };
    double maximum = std::max(error(left), error(left + h));
    if (amplitude > 0.0 && slope < 0.0) {
        const double ratio = -slope / (k_elimination * amplitude);
        if (ratio > 0.0) {
            const double critical = -std::log(ratio) / k_elimination;
            if (critical > left && critical < left + h) maximum = std::max(maximum, error(critical));
        }
    }
    return maximum;
}

double persistent(const std::vector<double>& starts,const std::vector<double>& bounds,double tolerance) {
    for (std::size_t i=0;i<bounds.size();++i) {
        bool suffix=true;
        for (std::size_t j=i;j<bounds.size();++j) suffix=suffix&&bounds[j]<=tolerance;
        if (suffix) return starts[i];
    }
    return std::numeric_limits<double>::quiet_NaN();
}

Budget make_budget(const Grid& grid,const ReferenceFit& fit,double discrete_fit,
                   const Trajectory<1>& unit_path) {
    Budget b;
    const double radius=std::max(std::abs(fit.lower-fit.estimate),std::abs(fit.upper-fit.estimate));
    std::vector<double> auc_unit(grid.time.size(),0.0);
    for (std::size_t n=grid.steps();n-->0;) {
        const double h=grid.time[n+1]-grid.time[n],th=grid.theta[n];
        auc_unit[n]=auc_unit[n+1]+h*((1-th)*unit_path.state[n][0]+th*unit_path.state[n+1][0]);
    }
    std::vector<double> starts;
    for(std::size_t n=0;n<grid.steps();++n){
        const double left=grid.time[n],h=grid.time[n+1]-left;
        const double mc=radius*std::exp(-k_elimination*left);
        const double nc=exp_line_error(left,h,discrete_fit*unit_path.state[n][0],
                                       discrete_fit*unit_path.state[n+1][0],fit.estimate);
        const double ma=radius*(std::exp(-k_elimination*left)-std::exp(-k_elimination*final_time))/k_elimination;
        const double na=exp_line_error(left,h,discrete_fit*auc_unit[n],discrete_fit*auc_unit[n+1],
                                       fit.estimate/k_elimination,
                                       -fit.estimate*std::exp(-k_elimination*final_time)/k_elimination);
        b.measurement_c.push_back(mc);b.numerical_c.push_back(nc);b.total_c.push_back(mc+nc);
        b.measurement_auc.push_back(ma);b.numerical_auc.push_back(na);b.total_auc.push_back(ma+na);
        starts.push_back(left);
    }
    b.measurement_c_time=persistent(starts,b.measurement_c,concentration_tolerance);
    b.numerical_c_time=persistent(starts,b.numerical_c,concentration_tolerance);
    b.total_c_time=persistent(starts,b.total_c,concentration_tolerance);
    b.measurement_auc_time=persistent(starts,b.measurement_auc,auc_tolerance);
    b.numerical_auc_time=persistent(starts,b.numerical_auc,auc_tolerance);
    b.total_auc_time=persistent(starts,b.total_auc,auc_tolerance);
    return b;
}

double analytic_cost(double c0,const Measurement& data,const std::vector<std::vector<double>>& W) {
    std::array<double,6> r{};
    for(std::size_t j=0;j<r.size();++j)r[j]=concentration(sample_times[j],c0)-data.value[j];
    return .5*quadratic(r,W);
}
double analytic_gradient(double c0,const Measurement& data,const std::vector<std::vector<double>>& W) {
    std::array<double,6> r{},s{};
    for(std::size_t j=0;j<r.size();++j){s[j]=std::exp(-k_elimination*sample_times[j]);r[j]=s[j]*c0-data.value[j];}
    double g=0.;for(std::size_t i=0;i<r.size();++i)for(std::size_t j=0;j<r.size();++j)g+=s[i]*W[i][j]*r[j];
    return g;
}
double analytic_hessian(const std::vector<std::vector<double>>& W) {
    std::array<double,6>s{};for(std::size_t j=0;j<s.size();++j)s[j]=std::exp(-k_elimination*sample_times[j]);
    double h=0.;for(std::size_t i=0;i<s.size();++i)for(std::size_t j=0;j<s.size();++j)h+=s[i]*W[i][j]*s[j];
    return h;
}

ErrorComponents convergence_errors(std::size_t steps,double theta,const Measurement& data,
                                  const std::vector<std::vector<double>>& W) {
    const Grid grid=Grid::uniform(0.,final_time,steps,theta);
    const System<1,1> system(washout_problem(grid,data,W),grid);
    const Vec<1> truth{true_initial};
    const auto trajectory=system.forward(truth);
    const auto tangent=system.tangent(truth,trajectory,{1.});
    const Vec<1> diagnostic{true_initial+1.0};
    const auto diagnostic_path=system.forward(diagnostic);
    const auto evaluation=system.adjoint(diagnostic,diagnostic_path);
    const auto hessian=system.hessian_vector(diagnostic,{1.});
    const auto fit=minimize(system,Vec<1>{8.},OptimizeOptions{1e-9,40,40},
                            [](const Vec<1>& c0){return c0[0]>0.;});
    const System<1,1> exposure(exposure_problem(),grid);
    const auto auc_path=exposure.forward(truth);
    ErrorComponents e;
    for(std::size_t n=0;n<=steps;++n){
        e.state=std::max(e.state,std::abs(trajectory.state[n][0]-concentration(grid.time[n],true_initial)));
        e.tangent=std::max(e.tangent,std::abs(tangent[n][0]-std::exp(-k_elimination*grid.time[n])));
    }
    const auto exact_fit=analytical_fit(data,W);
    e.cost=std::abs(evaluation.cost-analytic_cost(diagnostic[0],data,W));
    e.gradient=std::abs(evaluation.gradient[0]-analytic_gradient(diagnostic[0],data,W));
    e.hessian=std::abs(hessian[0]-analytic_hessian(W));
    e.fit=std::abs(fit.control[0]-exact_fit.estimate);
    e.auc=std::abs(exposure.cost(truth,auc_path)-remaining_auc(0.,true_initial));
    return e;
}

void run_convergence(const std::filesystem::path& output,const Measurement& data,
                     const std::vector<std::vector<double>>& correlated_W,
                     const std::vector<std::vector<double>>& diagonal_W,Checks& checks) {
    std::ofstream csv(output/"convergence.csv"),orders(output/"orders.csv");
    csv<<std::setprecision(17)<<"model,theta,steps,h,state_error,tangent_error,cost_error,gradient_error,hessian_error,fit_error,auc_error\n";
    orders<<std::setprecision(17)<<"model,theta,quantity,observed_order,expected_order,passed\n";
    const std::array<const char*,7> names{"state","tangent","cost","gradient","hessian","fit","auc"};
    for(const auto& model:{std::pair{"correlated",&correlated_W},std::pair{"diagonal_assumption",&diagonal_W}}){
        for(double theta:{.5,.65}){
            std::array<ErrorComponents,convergence_steps.size()> errors{};
            for(std::size_t i=0;i<convergence_steps.size();++i){
                const auto N=convergence_steps[i];errors[i]=convergence_errors(N,theta,data,*model.second);
                const double h=final_time/static_cast<double>(N);const auto&e=errors[i];
                csv<<model.first<<','<<theta<<','<<N<<','<<h<<','<<e.state<<','<<e.tangent<<','<<e.cost<<','<<e.gradient<<','<<e.hessian<<','<<e.fit<<','<<e.auc<<'\n';
            }
            const double expected=theta==.5?2.:1.;
            for(std::size_t q=0;q<names.size();++q){
                auto component=[q](const ErrorComponents&e){return std::array<double,7>{e.state,e.tangent,e.cost,e.gradient,e.hessian,e.fit,e.auc}[q];};
                const double coarse=component(errors[errors.size()-2]),fine=component(errors.back());
                if(!(coarse>0.&&fine>0.))throw std::runtime_error("Convergence quantity reached roundoff");
                const double observed=std::log(coarse/fine)/std::log(2.);
                const double tolerance=theta==.5?.30:.35;
                const double err=std::abs(observed-expected);
                const bool passed=std::isfinite(observed)&&err<=tolerance;
                orders<<model.first<<','<<theta<<','<<names[q]<<','<<observed<<','<<expected<<','<<passed<<'\n';
                checks.check(std::string(model.first)+" theta="+std::to_string(theta)+" "+names[q]+" order",err,tolerance);
                std::cout<<"  "<<model.first<<" theta="<<theta<<" / "<<names[q]<<" order="<<observed<<" (expected "<<expected<<")\n";
            }
        }
    }
}

void run_case(const char* model,const std::vector<std::vector<double>>& W,
              const Measurement& data,const ReferenceFit& reference,std::size_t steps,
              Checks& checks,std::ofstream& comparisons,std::ofstream& decisions,
              std::ofstream& budgets,std::ofstream& paths) {
    const Grid grid=Grid::uniform(0.,final_time,steps,.5);
    const System<1,1> system(washout_problem(grid,data,W),grid);
    const auto optimum=minimize(system,Vec<1>{8.},OptimizeOptions{1e-9,40,40},
                                [](const Vec<1>& c0){return c0[0]>0.;});
    const auto evaluation=system.evaluate(optimum.control);
    const double closed=discrete_closed_fit(system,grid,data,W);
    const auto unit=system.forward({1.});
    checks.check(std::string(model)+" discrete fit equals closed form",std::abs(optimum.control[0]-closed),2e-8);
    checks.check(std::string(model)+" discrete adjoint stationarity",std::abs(evaluation.gradient[0]),2e-8);
    checks.check(std::string(model)+" exact state lies in 95% joint ellipsoid",
                 std::max(0.,reference.true_quadratic-chi_square_6_95),1e-12);
    const auto budget=make_budget(grid,reference,optimum.control[0],unit);
    const System<1,1> exposure(exposure_problem(),grid);
    const auto exp_eval=exposure.evaluate(optimum.control);
    const double auc_gradient=(1.-std::exp(-k_elimination*final_time))/k_elimination;
    checks.check(std::string(model)+" theta-quadrature exposure gradient",
                 std::abs(exp_eval.gradient[0]-auc_gradient),steps==512?2e-5:.02);

    comparisons<<model<<','<<steps<<','<<reference.estimate<<','<<optimum.control[0]<<','
               <<reference.lower<<','<<reference.upper<<','<<(reference.upper-reference.lower)<<','
               <<reference.information<<','<<reference.minimum_quadratic<<','
               <<reference.true_quadratic<<','<<std::max(0.,reference.true_quadratic-chi_square_6_95)<<'\n';
    decisions<<model<<','<<steps<<",concentration,"<<concentration_tolerance<<','
             <<budget.measurement_c_time<<','<<budget.numerical_c_time<<','<<budget.total_c_time<<'\n';
    decisions<<model<<','<<steps<<",remaining_auc,"<<auc_tolerance<<','
             <<budget.measurement_auc_time<<','<<budget.numerical_auc_time<<','<<budget.total_auc_time<<'\n';
    std::cout<<"\n"<<model<<" covariance model, N="<<steps<<"\n"
             <<"  continuous weighted fit="<<reference.estimate<<" mg/L; discrete fit="<<optimum.control[0]<<" mg/L\n"
             <<(std::string(model)=="correlated" ? "  nominal 95% joint-ellipsoid interval=[" :
                                                    "  interval under diagonal covariance assumption=[")
             <<reference.lower<<", "<<reference.upper<<"] mg/L\n"
             <<"  concentration qualification: measurement="<<budget.measurement_c_time
             <<" h; numerical="<<budget.numerical_c_time<<" h; combined="<<budget.total_c_time<<" h\n"
             <<"  remaining-AUC qualification: measurement="<<budget.measurement_auc_time
             <<" h; numerical="<<budget.numerical_auc_time<<" h; combined="<<budget.total_auc_time<<" h\n";
    std::vector<double> auc_unit(grid.time.size(),0.);
    for(std::size_t n=grid.steps();n-->0;){double h=grid.time[n+1]-grid.time[n];auc_unit[n]=auc_unit[n+1]+h*.5*(unit.state[n][0]+unit.state[n+1][0]);}
    const auto discrete_path=system.forward(optimum.control);
    const double radius=std::max(std::abs(reference.lower-reference.estimate),std::abs(reference.upper-reference.estimate));
    for(std::size_t n=0;n<grid.steps();++n){
        for(const auto& goal:{std::pair{"concentration",concentration_tolerance},std::pair{"remaining_auc",auc_tolerance}}){
            const bool is_c=std::string(goal.first)=="concentration";
            const double m=is_c?budget.measurement_c[n]:budget.measurement_auc[n];
            const double nu=is_c?budget.numerical_c[n]:budget.numerical_auc[n];
            const double total=is_c?budget.total_c[n]:budget.total_auc[n];
            budgets<<model<<','<<steps<<','<<grid.time[n]<<','<<grid.time[n+1]<<','<<goal.first<<','<<goal.second<<','<<m<<','<<nu<<','<<total<<'\n';
        }
    }
    for(std::size_t n=0;n<=steps;++n){
        const double t=grid.time[n];
        const double analytic_c=concentration(t,reference.estimate);
        const double discrete_c=discrete_path.state[n][0];
        const double analytic_a=remaining_auc(t,reference.estimate);
        const double discrete_a=optimum.control[0]*auc_unit[n];
        const double mc=radius*std::exp(-k_elimination*t);
        const double ma=radius*(std::exp(-k_elimination*t)-std::exp(-k_elimination*final_time))/k_elimination;
        paths<<model<<','<<steps<<','<<t<<','<<concentration(t,true_initial)<<','<<analytic_c<<','<<discrete_c<<','
             <<mc<<','<<std::abs(analytic_c-discrete_c)<<','<<mc+std::abs(analytic_c-discrete_c)<<','
             <<remaining_auc(t,true_initial)<<','<<analytic_a<<','<<discrete_a<<','<<ma<<','
             <<std::abs(analytic_a-discrete_a)<<','<<ma+std::abs(analytic_a-discrete_a)<<'\n';
    }
}

} // namespace

int main(int argc,char**argv){
    try{
        namespace fs=std::filesystem;
        const fs::path output=argc>1?argv[1]:"results/drug_washout_correlated";
        fs::create_directories(output);
        Checks checks(output/"checks.csv");
        const auto data=make_observations();
        const auto correlated_W=precision_matrix(true),diagonal_W=precision_matrix(false);
        const auto correlated_fit=analytical_fit(data,correlated_W);
        const auto independent_fit=analytical_fit(data,diagonal_W);
        const double innovation_norm=std::inner_product(innovations.begin(),innovations.end(),innovations.begin(),0.0);
        checks.check("fixed standardized innovation vector is inside 95% chi-square-6 ellipsoid",
                     std::max(0.,innovation_norm-chi_square_6_95),1e-12);

        std::ofstream comparisons(output/"covariance_comparison.csv");
        std::ofstream decisions(output/"decisions.csv");
        std::ofstream budgets(output/"error_budgets.csv");
        std::ofstream paths(output/"reconstruction.csv");
        comparisons<<std::setprecision(17)<<"model,steps,analytical_fit_mg_L,discrete_fit_mg_L,confidence_lower_mg_L,confidence_upper_mg_L,confidence_width_mg_L,information,minimum_quadratic,true_quadratic,true_quadratic_excess\n";
        decisions<<std::setprecision(17)<<"model,steps,goal,tolerance,measurement_only_time_h,numerical_only_time_h,combined_time_h\n";
        budgets<<std::setprecision(17)<<"model,steps,elapsed_h,interval_end_h,goal,tolerance,measurement_component,numerical_component,total_bound\n";
        paths<<std::setprecision(17)<<"model,steps,elapsed_h,true_concentration_mg_L,analytical_fit_mg_L,discrete_fit_mg_L,measurement_concentration_bound,numerical_concentration_bound,total_concentration_bound,true_remaining_auc_mg_h_L,analytical_fit_remaining_auc_mg_h_L,discrete_fit_remaining_auc_mg_h_L,measurement_auc_bound,numerical_auc_bound,total_auc_bound\n";

        std::cout<<"CORRELATED ASSAY ERRORS IN DRUG WASHOUT\n"
                 <<"One-compartment model; known k="<<k_elimination<<" 1/h; six fixed sample times.\n"
                 <<"Marginal assay SD="<<sensor_sd<<" mg/L; AR(1) correlation rho="<<correlation<<".\n"
                 <<"Same fixed data are fitted with the correct covariance and with its diagonal approximation.\n"
                 <<"95% data-consistency set: residual^T R^-1 residual <= chi-square(6,0.95)="<<chi_square_6_95<<".\n"
                 <<"The diagonal-covariance comparison is misspecified for these AR(1)-generated data; its interval is not a calibrated 95% claim.\n"
                 <<"The fixed standardized innovation norm squared is "<<innovation_norm<<".\n";
        std::cout<<"Analytical fits: correlated="<<correlated_fit.estimate<<" mg/L; diagonal assumption="
                 <<independent_fit.estimate<<" mg/L.\n"
                 <<"95% interval widths: correlated="<<correlated_fit.upper-correlated_fit.lower
                 <<" mg/L; diagonal assumption="<<independent_fit.upper-independent_fit.lower<<" mg/L.\n";

        run_convergence(output,data,correlated_W,diagonal_W,checks);
        for(std::size_t N:qualification_steps){
            run_case("correlated",correlated_W,data,correlated_fit,N,checks,comparisons,decisions,budgets,paths);
            run_case("diagonal_assumption",diagonal_W,data,independent_fit,N,checks,comparisons,decisions,budgets,paths);
        }
        checks.check("correct correlated model has lower information than diagonal approximation",
                     std::max(0.,correlated_fit.information-independent_fit.information),1e-12);
        checks.check("correlated confidence interval is wider than diagonal approximation",
                     std::max(0.,(independent_fit.upper-independent_fit.lower)-
                                  (correlated_fit.upper-correlated_fit.lower)),1e-12);
        std::cout<<"  covariance-model fit shift="
                 <<correlated_fit.estimate-independent_fit.estimate
                 <<" mg/L; confidence-width ratio="
                 <<(correlated_fit.upper-correlated_fit.lower)/
                   (independent_fit.upper-independent_fit.lower)<<".\n";
        std::cout<<"Checks passed: "<<checks.count<<". Correlated-observation experiment PASS.\n"
                 <<"Results: \""<<output.string()<<"\"\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
