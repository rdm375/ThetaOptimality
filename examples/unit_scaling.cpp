#include <theta_opt/theta_opt.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace {
using namespace theta_opt;

constexpr double k_hour = 0.45;
constexpr double c0_hour_mg_l = 10.0;
constexpr double horizon_hour = 8.0;
constexpr double concentration_scale = 1000.0; // microgram/L per mg/L
constexpr double time_scale = 60.0;             // minutes per hour
constexpr double sensor_sd_mg_l = 0.05;
constexpr double concentration_goal_mg_l = 0.012;
constexpr double auc_goal_mg_h_l = 0.040;
constexpr std::array<double, 6> sample_hour{0.5, 1.0, 2.0, 4.0, 6.0, 8.0};
constexpr std::array<double, 6> offset_mg_l{0.0010, -0.0015, 0.0007,
                                            -0.0008, 0.0012, -0.0004};
constexpr std::array<std::size_t, 5> meshes{32, 64, 128, 256, 512};

struct Units {
    const char* name;
    const char* concentration_label;
    const char* time_label;
    double concentration_factor;
    double time_factor;
    double rate;
    double horizon;
    double true_initial;
    double observation_sd;
    double concentration_tolerance;
    double auc_tolerance;
};

const Units mg_hour{"mg_L_hour", "mg/L", "h", 1.0, 1.0, k_hour,
                    horizon_hour, c0_hour_mg_l, sensor_sd_mg_l,
                    concentration_goal_mg_l, auc_goal_mg_h_l};
const Units ug_min{"ug_L_min", "microgram/L", "min", concentration_scale,
                   time_scale, k_hour / time_scale, horizon_hour * time_scale,
                   c0_hour_mg_l * concentration_scale,
                   sensor_sd_mg_l * concentration_scale,
                   concentration_goal_mg_l * concentration_scale,
                   auc_goal_mg_h_l * concentration_scale * time_scale};

struct Checks {
    std::ofstream out;
    std::size_t count{};
    explicit Checks(const std::filesystem::path& path) : out(path) {
        out << std::setprecision(17) << "check,error,tolerance,passed\n";
        if (!out) throw std::runtime_error("Cannot open checks.csv");
    }
    void close(double error, double tolerance, const std::string& name) {
        const bool passed = std::isfinite(error) && error <= tolerance;
        out << name << ',' << error << ',' << tolerance << ',' << passed << '\n';
        ++count;
        if (!passed) throw std::runtime_error(name + " failed; error=" +
                                               std::to_string(error));
    }
    void equal(double a, double b, double tolerance, const std::string& name) {
        close(std::abs(a - b) / std::max({1.0, std::abs(a), std::abs(b)}),
              tolerance, name);
    }
};

double exact_state(double time, const Units& unit, double initial) {
    return initial * std::exp(-unit.rate * time);
}

double exact_auc(double time, const Units& unit, double initial) {
    return initial * (std::exp(-unit.rate * time) -
                      std::exp(-unit.rate * unit.horizon)) / unit.rate;
}

std::array<double, sample_hour.size()> observations(const Units& unit) {
    std::array<double, sample_hour.size()> data{};
    for (std::size_t j = 0; j < data.size(); ++j) {
        const double time = sample_hour[j] * unit.time_factor;
        data[j] = exact_state(time, unit, unit.true_initial) +
                  offset_mg_l[j] * unit.concentration_factor;
    }
    return data;
}

std::size_t node_for(const Grid& grid, double time) {
    for (std::size_t n = 0; n < grid.time.size(); ++n)
        if (std::abs(grid.time[n] - time) < 1e-10) return n;
    throw std::runtime_error("Sample time is missing from grid");
}

Problem<1, 1> observation_problem(const Grid& grid, const Units& unit,
                                  const std::array<double, 6>& data) {
    Problem<1, 1> problem;
    problem.initial.value = [](const Vec<1>& xi) { return xi; };
    problem.initial.jacobian = [](const Vec<1>&) { return identity<1>(); };
    problem.initial.directional_jacobian = [](const Vec<1>&, const Vec<1>&) {
        return Mat<1, 1>{};
    };
    problem.dynamics.value = [&unit](double, const Vec<1>& x, const Vec<1>&) {
        return Vec<1>{-unit.rate * x[0]};
    };
    problem.dynamics.derivatives = [&unit](double, const Vec<1>&, const Vec<1>&) {
        Derivatives<1, 1> d;
        d.x[0][0] = -unit.rate;
        return d;
    };
    problem.dynamics.directional_derivatives = [](double, const Vec<1>&,
        const Vec<1>&, const Vec<1>&, const Vec<1>&) { return Derivatives<1, 1>{}; };
    const double variance = unit.observation_sd * unit.observation_sd;
    for (std::size_t j = 0; j < data.size(); ++j) {
        const double observed = data[j];
        auto term = ScalarFunction<1, 1>::zero();
        term.evaluate = [observed, variance](double, const Vec<1>& x, const Vec<1>&) {
            const double residual = x[0] - observed;
            ScalarValue<1, 1> value;
            value.value = 0.5 * residual * residual / variance;
            value.x[0] = residual / variance;
            return value;
        };
        term.directional_gradient = [variance](double, const Vec<1>&,
            const Vec<1>&, const Vec<1>& direction, const Vec<1>&) {
            ScalarDirection<1, 1> value;
            value.x[0] = direction[0] / variance;
            return value;
        };
        problem.nodal.emplace_back(node_for(grid, sample_hour[j] * unit.time_factor),
                                   std::move(term));
    }
    return problem;
}

Problem<1, 1> auc_problem(const Units& unit) {
    Problem<1, 1> problem;
    problem.initial.value = [](const Vec<1>& xi) { return xi; };
    problem.initial.jacobian = [](const Vec<1>&) { return identity<1>(); };
    problem.initial.directional_jacobian = [](const Vec<1>&, const Vec<1>&) {
        return Mat<1, 1>{};
    };
    problem.dynamics.value = [&unit](double, const Vec<1>& x, const Vec<1>&) {
        return Vec<1>{-unit.rate * x[0]};
    };
    problem.dynamics.derivatives = [&unit](double, const Vec<1>&, const Vec<1>&) {
        Derivatives<1, 1> d;
        d.x[0][0] = -unit.rate;
        return d;
    };
    problem.dynamics.directional_derivatives = [](double, const Vec<1>&,
        const Vec<1>&, const Vec<1>&, const Vec<1>&) { return Derivatives<1, 1>{}; };
    problem.running.evaluate = [](double, const Vec<1>& x, const Vec<1>&) {
        ScalarValue<1, 1> value;
        value.value = x[0];
        value.x[0] = 1.0;
        return value;
    };
    problem.running.directional_gradient = [](double, const Vec<1>&,
        const Vec<1>&, const Vec<1>&, const Vec<1>&) { return ScalarDirection<1, 1>{}; };
    return problem;
}

double continuous_fit(const Units& unit, const std::array<double, 6>& data) {
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t j = 0; j < data.size(); ++j) {
        const double time = sample_hour[j] * unit.time_factor;
        const double sensitivity = std::exp(-unit.rate * time);
        numerator += sensitivity * data[j];
        denominator += sensitivity * sensitivity;
    }
    return numerator / denominator;
}

std::pair<double, double> exact_data_interval(const Units& unit,
                                               const std::array<double, 6>& data) {
    double lower = 0.0;
    double upper = 20.0 * unit.concentration_factor;
    for (std::size_t j = 0; j < data.size(); ++j) {
        const double time = sample_hour[j] * unit.time_factor;
        const double sensitivity = std::exp(-unit.rate * time);
        lower = std::max(lower, (data[j] - unit.observation_sd) / sensitivity);
        upper = std::min(upper, (data[j] + unit.observation_sd) / sensitivity);
    }
    if (!(lower <= upper)) throw std::runtime_error("Empty bounded-data set");
    return {lower, upper};
}

double line_exponential_error(double left_time, double h, double left,
                              double right, double amplitude, double rate,
                              double offset) {
    const double slope = (right - left) / h;
    auto error = [&](double time) {
        const double alpha = (time - left_time) / h;
        const double line = (1.0 - alpha) * left + alpha * right;
        return line - (amplitude * std::exp(-rate * time) + offset);
    };
    double maximum = std::max(std::abs(error(left_time)), std::abs(error(left_time + h)));
    if (amplitude > 0.0 && slope < 0.0) {
        const double ratio = -slope / (rate * amplitude);
        if (ratio > 0.0) {
            const double critical = -std::log(ratio) / rate;
            if (critical > left_time && critical < left_time + h)
                maximum = std::max(maximum, std::abs(error(critical)));
        }
    }
    return maximum;
}

double persistent_time(const std::vector<double>& bound, double tolerance,
                       const Grid& grid) {
    for (std::size_t n = 0; n < bound.size(); ++n) {
        bool suffix = true;
        for (std::size_t j = n; j < bound.size(); ++j) suffix = suffix && bound[j] <= tolerance;
        if (suffix) return grid.time[n];
    }
    return std::numeric_limits<double>::infinity();
}

struct Metrics {
    double state{}, tangent{}, cost{}, gradient{}, hessian{}, fit{}, auc{};
};

Metrics convergence_metrics(const Units& unit, const std::array<double, 6>& data,
                           std::size_t steps, double theta) {
    const auto grid = Grid::uniform(0.0, unit.horizon, steps, theta);
    const System<1, 1> fit_system(observation_problem(grid, unit, data), grid);
    const Vec<1> truth{unit.true_initial};
    const auto state = fit_system.forward(truth);
    const auto tangent = fit_system.tangent(truth, state, {1.0});
    Metrics error;
    for (std::size_t n = 0; n <= steps; ++n) {
        error.state = std::max(error.state,
            std::abs(state.state[n][0] - exact_state(grid.time[n], unit, unit.true_initial)) /
            unit.concentration_factor);
        error.tangent = std::max(error.tangent,
            std::abs(tangent[n][0] - std::exp(-unit.rate * grid.time[n])));
    }
    const Vec<1> diagnostic{11.0 * unit.concentration_factor};
    const auto diagnostic_trajectory = fit_system.forward(diagnostic);
    const auto evaluation = fit_system.adjoint(diagnostic, diagnostic_trajectory);
    const auto hessian = fit_system.hessian_vector(diagnostic, {1.0});
    double exact_cost = 0.0, exact_gradient = 0.0, exact_hessian = 0.0;
    for (std::size_t j = 0; j < data.size(); ++j) {
        const double time = sample_hour[j] * unit.time_factor;
        const double sensitivity = std::exp(-unit.rate * time);
        const double residual = exact_state(time, unit, diagnostic[0]) - data[j];
        const double variance = unit.observation_sd * unit.observation_sd;
        exact_cost += 0.5 * residual * residual / variance;
        exact_gradient += residual * sensitivity / variance;
        exact_hessian += sensitivity * sensitivity / variance;
    }
    error.cost = std::abs(evaluation.cost - exact_cost);
    error.gradient = std::abs(evaluation.gradient[0] - exact_gradient) *
                     unit.concentration_factor;
    error.hessian = std::abs(hessian[0] - exact_hessian) *
                    unit.concentration_factor * unit.concentration_factor;
    const auto optimum = minimize(fit_system, Vec<1>{8.0 * unit.concentration_factor},
        OptimizeOptions{1e-10, 40, 40}, [](const Vec<1>& p) { return p[0] > 0.0; });
    error.fit = std::abs(optimum.control[0] - continuous_fit(unit, data)) /
                unit.concentration_factor;

    const System<1, 1> exposure(auc_problem(unit), grid);
    const auto exposure_evaluation = exposure.evaluate(truth);
    error.auc = std::abs(exposure_evaluation.cost - exact_auc(0.0, unit, unit.true_initial)) /
                (unit.concentration_factor * unit.time_factor);
    return error;
}

std::vector<double> remaining_discrete_auc(const Grid& grid, const Trajectory<1>& state) {
    std::vector<double> value(grid.time.size(), 0.0);
    for (std::size_t n = grid.steps(); n-- > 0;) {
        const double h = grid.time[n + 1] - grid.time[n];
        const double theta = grid.theta[n];
        value[n] = value[n + 1] + h * ((1.0 - theta) * state.state[n][0] +
                                       theta * state.state[n + 1][0]);
    }
    return value;
}

struct Qualification { double concentration{}, auc{}; };

Qualification qualification(const Units& unit, const Grid& grid,
                            const std::array<double, 6>& data,
                            double discrete_fit, const Trajectory<1>& unit_state) {
    const auto feasible = exact_data_interval(unit, data);
    const double radius = std::max(std::abs(feasible.first - continuous_fit(unit, data)),
                                   std::abs(feasible.second - continuous_fit(unit, data)));
    const auto auc_nodes_unit = remaining_discrete_auc(grid, unit_state);
    std::vector<double> concentration_bound, auc_bound;
    for (std::size_t n = 0; n < grid.steps(); ++n) {
        const double time = grid.time[n], h = grid.time[n + 1] - time;
        const double measurement_c = radius * std::exp(-unit.rate * time);
        const double exact_fit = continuous_fit(unit, data);
        const double numerical_c = line_exponential_error(time, h,
            discrete_fit * unit_state.state[n][0],
            discrete_fit * unit_state.state[n + 1][0], exact_fit, unit.rate, 0.0);
        const double measurement_auc = radius * (std::exp(-unit.rate * time) -
            std::exp(-unit.rate * unit.horizon)) / unit.rate;
        const double auc_offset = -exact_fit *
            std::exp(-unit.rate * unit.horizon) / unit.rate;
        const double numerical_auc = line_exponential_error(time, h,
            discrete_fit * auc_nodes_unit[n], discrete_fit * auc_nodes_unit[n + 1],
            exact_fit / unit.rate, unit.rate, auc_offset);
        concentration_bound.push_back(measurement_c + numerical_c);
        auc_bound.push_back(measurement_auc + numerical_auc);
    }
    return {persistent_time(concentration_bound, unit.concentration_tolerance, grid),
            persistent_time(auc_bound, unit.auc_tolerance, grid)};
}

void write_csv(const std::filesystem::path& directory, Checks& checks) {
    const auto data_a = observations(mg_hour), data_b = observations(ug_min);
    checks.equal(ug_min.true_initial / concentration_scale, mg_hour.true_initial,
                 1e-15, "initial-state unit conversion");
    checks.equal(ug_min.rate * time_scale, mg_hour.rate, 1e-15,
                 "rate unit conversion");
    checks.equal(ug_min.observation_sd / concentration_scale, mg_hour.observation_sd,
                 1e-15, "measurement standard-deviation conversion");
    checks.equal(ug_min.concentration_tolerance / concentration_scale,
                 mg_hour.concentration_tolerance, 1e-15,
                 "concentration goal conversion");
    checks.equal(ug_min.auc_tolerance / (concentration_scale * time_scale),
                 mg_hour.auc_tolerance, 1e-15, "AUC goal conversion");
    checks.equal(ug_min.horizon / time_scale, mg_hour.horizon, 1e-15,
                 "horizon time conversion");
    for (std::size_t j = 0; j < data_a.size(); ++j) {
        checks.equal(data_a[j], data_b[j] / concentration_scale, 2e-15,
                     "observation conversion j=" + std::to_string(j));
        checks.equal(sample_hour[j] * mg_hour.time_factor,
                     sample_hour[j] * ug_min.time_factor / time_scale,
                     1e-15, "sample-time conversion j=" + std::to_string(j));
    }
    std::ofstream comparison(directory / "unit_comparison.csv");
    comparison << std::setprecision(17)
        << "steps,quantity,mg_L_hour,ug_L_min,ug_L_min_converted_back,relative_difference\n";
    auto row = [&](std::size_t steps, const std::string& quantity, double base, double alternate,
                   double converted) {
        comparison << steps << ',' << quantity << ',' << base << ',' << alternate << ',' << converted << ','
                   << std::abs(base - converted) / std::max(1.0, std::abs(base)) << '\n';
    };

    std::ofstream trajectories(directory / "trajectory_comparison.csv");
    trajectories << std::setprecision(17)
        << "time_h,mg_L,ug_L_converted_to_mg_L,absolute_difference_mg_L\n";
    std::ofstream qualifications(directory / "qualification_comparison.csv");
    qualifications << std::setprecision(17)
        << "steps,goal,mg_L_hour_time_h,ug_L_min_time_min,converted_time_h,difference_h\n";
    for (const std::size_t steps : {std::size_t{32}, std::size_t{128}, std::size_t{512}}) {
        const auto grid_a = Grid::uniform(0.0, mg_hour.horizon, steps, 0.5);
        const auto grid_b = Grid::uniform(0.0, ug_min.horizon, steps, 0.5);
        const System<1, 1> system_a(observation_problem(grid_a, mg_hour, data_a), grid_a);
        const System<1, 1> system_b(observation_problem(grid_b, ug_min, data_b), grid_b);
        const Vec<1> initial_a{mg_hour.true_initial}, initial_b{ug_min.true_initial};
        const auto state_a = system_a.forward(initial_a), state_b = system_b.forward(initial_b);
        const auto tangent_a = system_a.tangent(initial_a, state_a, {1.0});
        const auto tangent_b = system_b.tangent(initial_b, state_b, {1.0});
        const auto opt_a = minimize(system_a, Vec<1>{8.0}, OptimizeOptions{1e-10, 40, 40},
                                    [](const Vec<1>& x) { return x[0] > 0.0; });
        const auto opt_b = minimize(system_b, Vec<1>{8000.0}, OptimizeOptions{1e-7, 40, 40},
                                    [](const Vec<1>& x) { return x[0] > 0.0; });
        const auto eval_a = system_a.evaluate({11.0});
        const auto eval_b = system_b.evaluate({11000.0});
        const auto hess_a = system_a.hessian_vector({11.0}, {1.0});
        const auto hess_b = system_b.hessian_vector({11000.0}, {1.0});
        const auto auc_a = System<1, 1>(auc_problem(mg_hour), grid_a);
        const auto auc_b = System<1, 1>(auc_problem(ug_min), grid_b);
        const auto auc_eval_a = auc_a.evaluate(initial_a), auc_eval_b = auc_b.evaluate(initial_b);
        const auto response_a = system_a.forward({1.0});
        const auto response_b = system_b.forward({1.0});
        const auto q_a = qualification(mg_hour, grid_a, data_a, opt_a.control[0], response_a);
        const auto q_b = qualification(ug_min, grid_b, data_b, opt_b.control[0], response_b);

        row(steps, "fit_initial_concentration", opt_a.control[0], opt_b.control[0],
            opt_b.control[0] / concentration_scale);
        row(steps, "objective_at_diagnostic_state", eval_a.cost, eval_b.cost, eval_b.cost);
        row(steps, "scaled_initial_gradient", eval_a.gradient[0],
            eval_b.gradient[0] * concentration_scale,
            eval_b.gradient[0] * concentration_scale);
        row(steps, "scaled_initial_hessian", hess_a[0], hess_b[0] * concentration_scale * concentration_scale,
            hess_b[0] * concentration_scale * concentration_scale);
        row(steps, "total_auc", auc_eval_a.cost, auc_eval_b.cost,
            auc_eval_b.cost / (concentration_scale * time_scale));
        row(steps, "auc_initial_sensitivity", auc_eval_a.gradient[0], auc_eval_b.gradient[0] / time_scale,
            auc_eval_b.gradient[0] / time_scale);

        for (std::size_t n = 0; n <= steps; ++n) {
            checks.equal(state_a.state[n][0], state_b.state[n][0] / concentration_scale,
                         2e-12, "trajectory conversion N=" + std::to_string(steps));
            checks.equal(tangent_a[n][0], tangent_b[n][0], 2e-12,
                         "tangent conversion N=" + std::to_string(steps));
            checks.equal(eval_a.nodal_adjoint[n][0],
                         eval_b.nodal_adjoint[n][0] * concentration_scale, 2e-11,
                         "nodal adjoint covector scaling N=" + std::to_string(steps));
            if (steps == 512) {
                const double mg = state_a.state[n][0];
                const double converted = state_b.state[n][0] / concentration_scale;
                trajectories << grid_a.time[n] << ',' << mg << ',' << converted << ','
                             << std::abs(mg - converted) << '\n';
            }
        }
        checks.equal(eval_a.cost, eval_b.cost, 2e-12,
                     "dimensionless objective invariance N=" + std::to_string(steps));
        checks.equal(eval_a.gradient[0], eval_b.gradient[0] * concentration_scale, 2e-10,
                     "initial gradient covector scaling N=" + std::to_string(steps));
        checks.equal(hess_a[0], hess_b[0] * concentration_scale * concentration_scale,
                     2e-10, "initial Hessian scaling N=" + std::to_string(steps));
        checks.equal(opt_a.control[0], opt_b.control[0] / concentration_scale, 2e-10,
                     "optimized initial concentration conversion N=" + std::to_string(steps));
        checks.equal(auc_eval_a.cost, auc_eval_b.cost /
                     (concentration_scale * time_scale), 2e-12,
                     "AUC unit conversion N=" + std::to_string(steps));
        checks.equal(auc_eval_a.gradient[0], auc_eval_b.gradient[0] / time_scale,
                     2e-12, "AUC sensitivity unit conversion N=" + std::to_string(steps));
        const double qbc_h = q_b.concentration / time_scale;
        const double qba_h = q_b.auc / time_scale;
        checks.equal(q_a.concentration, qbc_h, 2e-12,
                     "concentration qualification time conversion N=" + std::to_string(steps));
        checks.equal(q_a.auc, qba_h, 2e-12,
                     "AUC qualification time conversion N=" + std::to_string(steps));
        for (const auto& [goal, ta, tb] : {
            std::tuple<std::string, double, double>{"concentration", q_a.concentration, q_b.concentration},
            std::tuple<std::string, double, double>{"remaining_auc", q_a.auc, q_b.auc}}) {
            qualifications << steps << ',' << goal << ',' << ta << ',' << tb << ','
                           << tb / time_scale << ',' << std::abs(ta - tb / time_scale) << '\n';
        }
    }

    std::ofstream convergence(directory / "convergence.csv");
    convergence << std::setprecision(17)
        << "units,theta,steps,state,tangent,cost,gradient,hessian,fit,auc\n";
    std::ofstream orders(directory / "orders.csv");
    orders << std::setprecision(17) << "units,theta,quantity,observed_order,expected_order\n";
    const std::array<const Units*, 2> units{&mg_hour, &ug_min};
    const std::array<const char*, 7> names{"state", "tangent", "cost", "gradient",
                                           "hessian", "fit", "auc"};
    for (const auto* unit : units) {
        for (double theta : {0.5, 0.65}) {
            std::array<Metrics, meshes.size()> metrics{};
            for (std::size_t j = 0; j < meshes.size(); ++j) {
                metrics[j] = convergence_metrics(*unit, observations(*unit), meshes[j], theta);
                const auto& e = metrics[j];
                convergence << unit->name << ',' << theta << ',' << meshes[j] << ','
                    << e.state << ',' << e.tangent << ',' << e.cost << ',' << e.gradient << ','
                    << e.hessian << ',' << e.fit << ',' << e.auc << '\n';
            }
            const std::array<double, 7> first{metrics[3].state, metrics[3].tangent,
                metrics[3].cost, metrics[3].gradient, metrics[3].hessian,
                metrics[3].fit, metrics[3].auc};
            const std::array<double, 7> last{metrics[4].state, metrics[4].tangent,
                metrics[4].cost, metrics[4].gradient, metrics[4].hessian,
                metrics[4].fit, metrics[4].auc};
            for (std::size_t j = 0; j < names.size(); ++j) {
                const double observed = std::log(first[j] / last[j]) / std::log(2.0);
                const double expected = theta == 0.5 ? 2.0 : 1.0;
                orders << unit->name << ',' << theta << ',' << names[j] << ','
                       << observed << ',' << expected << '\n';
                checks.close(std::abs(observed - expected), 0.10,
                    std::string(unit->name) + " convergence order " + names[j] +
                    " theta=" + std::to_string(theta));
            }
        }
    }
    if (!comparison || !trajectories || !qualifications || !convergence || !orders)
        throw std::runtime_error("Failed writing unit-scaling CSV output");
}

} // namespace

int main(int argc, char** argv) {
    try {
        const std::filesystem::path directory = argc > 1 ? argv[1] : "results/unit_scaling";
        std::filesystem::create_directories(directory);
        Checks checks(directory / "checks.csv");
        write_csv(directory, checks);
        std::cout << "UNIT AND SCALING INVARIANCE: DRUG WASHOUT\n"
                  << "Same physical data and theta grid in mg/L-hours and microgram/L-minutes.\n"
                  << "c conversion factor=1000; time conversion factor=60.\n"
                  << "Dimensionless objective is invariant; gradients, Hessians, and AUC covectors transform by their units.\n"
                  << "At N=512, concentration and remaining-AUC qualification times agree after converting minutes to hours.\n"
                  << "Convergence checked for theta=0.5 (order 2) and theta=0.65 (order 1).\n"
                  << "Checks passed: " << checks.count << ". Unit-scaling experiment PASS.\n"
                  << "Results: \"" << directory.string() << "\"\n";
    } catch (const std::exception& error) {
        std::cerr << "unit_scaling: " << error.what() << '\n';
        return 1;
    }
}
