#include <theta_opt/theta_opt.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace {
using namespace theta_opt;

constexpr double k_elimination = 0.45; // 1/hour
constexpr double t0 = 0.0;
constexpr double tf = 8.0; // hours after an IV bolus
constexpr double true_initial = 10.0; // mg/L
constexpr double prior_lower = 0.0;
constexpr double prior_upper = 20.0;
constexpr double concentration_tolerance = 0.012; // mg/L
constexpr double remaining_auc_tolerance = 0.040; // mg hour/L
constexpr std::array<double, 6> sample_times{0.5, 1.0, 2.0, 4.0, 6.0, 8.0};
constexpr std::array<double, 6> fixed_sensor_offsets{
    0.0010, -0.0015, 0.0007, -0.0008, 0.0012, -0.0004}; // mg/L
constexpr std::array<std::size_t, 5> convergence_steps{32, 64, 128, 256, 512};
constexpr std::array<std::size_t, 3> qualification_steps{32, 128, 512};
constexpr std::array<double, 2> sensor_bounds{0.005, 0.050}; // mg/L

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
        if (!passed)
            throw std::runtime_error(name + " failed: error=" + std::to_string(error) +
                                     ", tolerance=" + std::to_string(tolerance));
    }
};

struct Measurement {
    std::array<double, sample_times.size()> value{};
};

struct Budget {
    std::vector<double> measurement_concentration;
    std::vector<double> numerical_concentration;
    std::vector<double> total_concentration;
    std::vector<double> measurement_remaining_auc;
    std::vector<double> numerical_remaining_auc;
    std::vector<double> total_remaining_auc;
    std::optional<double> measurement_concentration_time;
    std::optional<double> numerical_concentration_time;
    std::optional<double> total_concentration_time;
    std::optional<double> measurement_auc_time;
    std::optional<double> numerical_auc_time;
    std::optional<double> total_auc_time;
};

struct Scenario {
    std::size_t steps{};
    double error_bound{};
    double theta{};
    const char* name{};
};

struct ConvergenceErrors {
    double state{};
    double tangent{};
    double cost{};
    double gradient{};
    double hessian{};
    double recovered_initial{};
    double auc{};
};

double exact_concentration(double time, double initial) {
    return initial * std::exp(-k_elimination * (time - t0));
}

double exact_remaining_auc(double time, double initial) {
    return initial * (std::exp(-k_elimination * (time - t0)) -
                      std::exp(-k_elimination * (tf - t0))) / k_elimination;
}

Measurement observations() {
    Measurement data;
    for (std::size_t j = 0; j < sample_times.size(); ++j)
        data.value[j] = exact_concentration(sample_times[j], true_initial) + fixed_sensor_offsets[j];
    return data;
}

std::size_t sample_node(const Grid& grid, double sample_time) {
    for (std::size_t n = 0; n < grid.time.size(); ++n)
        if (std::abs(grid.time[n] - sample_time) < 1e-12) return n;
    throw std::runtime_error("A sample time is not a node of the prescribed grid");
}

Problem<1, 1> observation_problem(const Grid& grid, const Measurement& data,
                                  double error_bound) {
    Problem<1, 1> problem;
    problem.initial.value = [](const Vec<1>& xi) { return xi; };
    problem.initial.jacobian = [](const Vec<1>&) { return identity<1>(); };
    problem.initial.directional_jacobian = [](const Vec<1>&, const Vec<1>&) {
        return Mat<1, 1>{};
    };
    problem.dynamics.value = [](double, const Vec<1>& x, const Vec<1>&) {
        return Vec<1>{-k_elimination * x[0]};
    };
    problem.dynamics.derivatives = [](double, const Vec<1>&, const Vec<1>&) {
        Derivatives<1, 1> derivative;
        derivative.x[0][0] = -k_elimination;
        return derivative;
    };
    problem.dynamics.directional_derivatives = [](double, const Vec<1>&, const Vec<1>&,
                                                   const Vec<1>&, const Vec<1>&) {
        return Derivatives<1, 1>{};
    };
    const double variance = error_bound * error_bound;
    for (std::size_t j = 0; j < sample_times.size(); ++j) {
        const double value = data.value[j];
        auto cost = ScalarFunction<1, 1>::zero();
        cost.evaluate = [value, variance](double, const Vec<1>& x, const Vec<1>&) {
            const double residual = x[0] - value;
            ScalarValue<1, 1> result;
            result.value = 0.5 * residual * residual / variance;
            result.x[0] = residual / variance;
            return result;
        };
        cost.directional_gradient = [variance](double, const Vec<1>&, const Vec<1>&,
                                                const Vec<1>& direction, const Vec<1>&) {
            ScalarDirection<1, 1> result;
            result.x[0] = direction[0] / variance;
            return result;
        };
        problem.nodal.emplace_back(sample_node(grid, sample_times[j]), std::move(cost));
    }
    return problem;
}

Problem<1, 1> auc_problem() {
    Problem<1, 1> problem;
    problem.initial.value = [](const Vec<1>& xi) { return xi; };
    problem.initial.jacobian = [](const Vec<1>&) { return identity<1>(); };
    problem.initial.directional_jacobian = [](const Vec<1>&, const Vec<1>&) {
        return Mat<1, 1>{};
    };
    problem.dynamics.value = [](double, const Vec<1>& x, const Vec<1>&) {
        return Vec<1>{-k_elimination * x[0]};
    };
    problem.dynamics.derivatives = [](double, const Vec<1>&, const Vec<1>&) {
        Derivatives<1, 1> derivative;
        derivative.x[0][0] = -k_elimination;
        return derivative;
    };
    problem.dynamics.directional_derivatives = [](double, const Vec<1>&, const Vec<1>&,
                                                   const Vec<1>&, const Vec<1>&) {
        return Derivatives<1, 1>{};
    };
    problem.running.evaluate = [](double, const Vec<1>& x, const Vec<1>&) {
        ScalarValue<1, 1> result;
        result.value = x[0];
        result.x[0] = 1.0;
        return result;
    };
    problem.running.directional_gradient = [](double, const Vec<1>&, const Vec<1>&,
                                               const Vec<1>&, const Vec<1>&) {
        return ScalarDirection<1, 1>{};
    };
    return problem;
}

double analytic_fit(const Measurement& data) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t j = 0; j < sample_times.size(); ++j) {
        const double sensitivity = std::exp(-k_elimination * sample_times[j]);
        numerator += sensitivity * data.value[j];
        denominator += sensitivity * sensitivity;
    }
    return numerator / denominator;
}

double analytic_cost(double initial, const Measurement& data, double sd) {
    double result = 0.0;
    for (std::size_t j = 0; j < sample_times.size(); ++j) {
        const double residual = exact_concentration(sample_times[j], initial) - data.value[j];
        result += 0.5 * residual * residual / (sd * sd);
    }
    return result;
}

double analytic_gradient(double initial, const Measurement& data, double sd) {
    double result = 0.0;
    for (std::size_t j = 0; j < sample_times.size(); ++j) {
        const double sensitivity = std::exp(-k_elimination * sample_times[j]);
        const double residual = initial * sensitivity - data.value[j];
        result += residual * sensitivity / (sd * sd);
    }
    return result;
}

double analytic_hessian(double sd) {
    double result = 0.0;
    for (double time : sample_times) {
        const double sensitivity = std::exp(-k_elimination * time);
        result += sensitivity * sensitivity / (sd * sd);
    }
    return result;
}

double linear_exponential_maximum(double t_left, double step, double line_left,
                                  double line_right, double amplitude,
                                  double offset = 0.0) {
    const double slope = (line_right - line_left) / step;
    auto difference = [&](double time) {
        const double fraction = (time - t_left) / step;
        const double line = (1.0 - fraction) * line_left + fraction * line_right;
        return line - (amplitude * std::exp(-k_elimination * time) + offset);
    };
    double maximum = std::max(std::abs(difference(t_left)),
                              std::abs(difference(t_left + step)));
    if (amplitude > 0.0 && slope < 0.0) {
        const double ratio = -slope / (k_elimination * amplitude);
        if (ratio > 0.0) {
            const double critical = -std::log(ratio) / k_elimination;
            if (critical > t_left && critical < t_left + step)
                maximum = std::max(maximum, std::abs(difference(critical)));
        }
    }
    return maximum;
}

double persistent_time(const std::vector<double>& interval_starts,
                      const std::vector<double>& interval_bounds, double tolerance) {
    for (std::size_t i = 0; i < interval_bounds.size(); ++i) {
        bool suffix_qualified = true;
        for (std::size_t j = i; j < interval_bounds.size(); ++j)
            suffix_qualified = suffix_qualified && interval_bounds[j] <= tolerance;
        if (suffix_qualified) return interval_starts[i];
    }
    return std::numeric_limits<double>::quiet_NaN();
}

std::pair<double, double> measurement_interval(const Measurement& data, double error_bound) {
    double lower = prior_lower;
    double upper = prior_upper;
    for (std::size_t j = 0; j < sample_times.size(); ++j) {
        const double sensitivity = std::exp(-k_elimination * sample_times[j]);
        lower = std::max(lower, (data.value[j] - error_bound) / sensitivity);
        upper = std::min(upper, (data.value[j] + error_bound) / sensitivity);
    }
    if (!(lower <= upper)) throw std::runtime_error("Analytical measurement set is empty");
    return {lower, upper};
}

std::vector<double> auc_remaining_unit(const Grid& grid,
                                       const Trajectory<1>& unit_trajectory) {
    std::vector<double> remaining(grid.time.size(), 0.0);
    for (std::size_t n = grid.steps(); n-- > 0;) {
        const double h = grid.time[n + 1] - grid.time[n];
        const double theta = grid.theta[n];
        remaining[n] = remaining[n + 1] + h *
            ((1.0 - theta) * unit_trajectory.state[n][0] +
             theta * unit_trajectory.state[n + 1][0]);
    }
    return remaining;
}

Budget compute_budget(const Grid& grid, const std::pair<double, double>& measured_set,
                      double continuous_fit, double discrete_fit,
                      const Trajectory<1>& unit_trajectory) {
    Budget result;
    const double measurement_radius = std::max(std::abs(measured_set.first - continuous_fit),
                                               std::abs(measured_set.second - continuous_fit));
    const auto auc_unit = auc_remaining_unit(grid, unit_trajectory);
    const auto fit_trajectory = [&]() {
        std::vector<double> value(grid.time.size());
        for (std::size_t n = 0; n < value.size(); ++n)
            value[n] = discrete_fit * unit_trajectory.state[n][0];
        return value;
    }();
    const double auc_offset = -continuous_fit * std::exp(-k_elimination * tf) / k_elimination;
    std::vector<double> starts;
    for (std::size_t n = 0; n < grid.steps(); ++n) starts.push_back(grid.time[n]);

    for (std::size_t n = 0; n < grid.steps(); ++n) {
        const double left = grid.time[n];
        const double h = grid.time[n + 1] - left;
        const double measurement_c = measurement_radius * std::exp(-k_elimination * left);
        const double numerical_c = linear_exponential_maximum(
            left, h, fit_trajectory[n], fit_trajectory[n + 1], continuous_fit);
        const double measurement_auc = measurement_radius *
            (std::exp(-k_elimination * left) - std::exp(-k_elimination * tf)) / k_elimination;
        const double numerical_auc = linear_exponential_maximum(
            left, h, discrete_fit * auc_unit[n], discrete_fit * auc_unit[n + 1],
            continuous_fit / k_elimination, auc_offset);
        result.measurement_concentration.push_back(measurement_c);
        result.numerical_concentration.push_back(numerical_c);
        result.total_concentration.push_back(measurement_c + numerical_c);
        result.measurement_remaining_auc.push_back(measurement_auc);
        result.numerical_remaining_auc.push_back(numerical_auc);
        result.total_remaining_auc.push_back(measurement_auc + numerical_auc);
    }

    result.measurement_concentration_time = persistent_time(
        starts, result.measurement_concentration, concentration_tolerance);
    result.numerical_concentration_time = persistent_time(
        starts, result.numerical_concentration, concentration_tolerance);
    result.total_concentration_time = persistent_time(
        starts, result.total_concentration, concentration_tolerance);
    result.measurement_auc_time = persistent_time(
        starts, result.measurement_remaining_auc, remaining_auc_tolerance);
    result.numerical_auc_time = persistent_time(
        starts, result.numerical_remaining_auc, remaining_auc_tolerance);
    result.total_auc_time = persistent_time(
        starts, result.total_remaining_auc, remaining_auc_tolerance);
    return result;
}

ConvergenceErrors convergence_errors(std::size_t steps, double theta,
                                     const Measurement& data) {
    const Grid grid = Grid::uniform(t0, tf, steps, theta);
    const double fit_sd = sensor_bounds[0];
    const System<1, 1> assimilation(observation_problem(grid, data, fit_sd), grid);
    const Vec<1> truth{true_initial};
    const auto trajectory = assimilation.forward(truth);
    const auto tangent = assimilation.tangent(truth, trajectory, {1.0});
    // Measure cost and derivative convergence away from the near-zero-residual
    // truth. At the truth, the deliberately small fixed sensor offsets can
    // cancel the leading objective error and make an asymptotic order estimate
    // misleading. This fixed diagnostic state keeps the same observations and
    // objective while avoiding that accidental cancellation.
    const Vec<1> diagnostic{true_initial + 1.0};
    const auto diagnostic_trajectory = assimilation.forward(diagnostic);
    const auto evaluation = assimilation.adjoint(diagnostic, diagnostic_trajectory);
    const auto hessian = assimilation.hessian_vector(diagnostic, {1.0});
    const auto optimum = minimize(assimilation, Vec<1>{8.0}, OptimizeOptions{1e-9, 40, 40},
                                  [](const Vec<1>& xi) { return xi[0] > 0.0; });
    const System<1, 1> auc(auc_problem(), grid);
    const auto auc_trajectory = auc.forward(truth);
    const double discrete_auc = auc.cost(truth, auc_trajectory);

    ConvergenceErrors errors;
    for (std::size_t n = 0; n <= steps; ++n) {
        const double time = grid.time[n];
        errors.state = std::max(errors.state,
            std::abs(trajectory.state[n][0] - exact_concentration(time, true_initial)));
        errors.tangent = std::max(errors.tangent,
            std::abs(tangent[n][0] - std::exp(-k_elimination * time)));
    }
    errors.cost = std::abs(evaluation.cost - analytic_cost(diagnostic[0], data, fit_sd));
    errors.gradient = std::abs(evaluation.gradient[0] - analytic_gradient(diagnostic[0], data, fit_sd));
    errors.hessian = std::abs(hessian[0] - analytic_hessian(fit_sd));
    errors.recovered_initial = std::abs(optimum.control[0] - analytic_fit(data));
    errors.auc = std::abs(discrete_auc - exact_remaining_auc(t0, true_initial));
    return errors;
}

void write_convergence(const std::filesystem::path& output, const Measurement& data,
                       Checks& checks) {
    std::ofstream curves(output / "convergence.csv");
    std::ofstream orders(output / "orders.csv");
    curves << std::setprecision(17)
        << "theta,steps,h,state_error,tangent_error,cost_error,gradient_error,"
        << "hessian_error,recovered_initial_error,auc_error\n";
    orders << std::setprecision(17)
        << "theta,quantity,observed_order,expected_order,passed\n";
    const std::array<const char*, 7> names{
        "state", "tangent", "cost", "gradient", "hessian", "recovered_initial", "auc"};
    for (double theta : {0.5, 0.65}) {
        std::array<ConvergenceErrors, convergence_steps.size()> errors{};
        for (std::size_t i = 0; i < convergence_steps.size(); ++i) {
            const std::size_t steps = convergence_steps[i];
            errors[i] = convergence_errors(steps, theta, data);
            const double h = (tf - t0) / static_cast<double>(steps);
            const auto& e = errors[i];
            curves << theta << ',' << steps << ',' << h << ',' << e.state << ',' << e.tangent
                   << ',' << e.cost << ',' << e.gradient << ',' << e.hessian << ','
                   << e.recovered_initial << ',' << e.auc << '\n';
        }
        const double expected = theta == 0.5 ? 2.0 : 1.0;
        for (std::size_t q = 0; q < names.size(); ++q) {
            auto get = [q](const ConvergenceErrors& e) {
                const std::array<double, 7> v{e.state, e.tangent, e.cost, e.gradient,
                    e.hessian, e.recovered_initial, e.auc};
                return v[q];
            };
            const double coarse = get(errors[errors.size() - 2]);
            const double fine = get(errors.back());
            if (!(coarse > 0.0 && fine > 0.0))
                throw std::runtime_error("Convergence quantity reached roundoff: " +
                                         std::string(names[q]));
            const double observed = std::log(coarse / fine) / std::log(2.0);
            const double tolerance = theta == 0.5 ? 0.28 : 0.32;
            const double rate_error = std::abs(observed - expected);
            const bool passed = std::isfinite(observed) && rate_error <= tolerance;
            orders << theta << ',' << names[q] << ',' << observed << ',' << expected << ','
                   << passed << '\n';
            checks.check(std::string("theta=") + std::to_string(theta) + " " + names[q] +
                         " observed convergence order", rate_error, tolerance);
            std::cout << "  theta=" << theta << " / " << names[q] << " order="
                      << observed << " (expected " << expected << ")\n";
        }
    }
}

void run_scenario(const Scenario& scenario, const Measurement& data, Checks& checks,
                  std::ofstream& decisions, std::ofstream& budgets,
                  std::ofstream& reconstructions) {
    const Grid grid = Grid::uniform(t0, tf, scenario.steps, scenario.theta);
    for (double time : sample_times) (void)sample_node(grid, time);
    const System<1, 1> system(observation_problem(grid, data, scenario.error_bound), grid);
    const auto optimum = minimize(system, Vec<1>{8.0}, OptimizeOptions{1e-9, 40, 40},
                                  [](const Vec<1>& xi) { return xi[0] > 0.0; });
    const auto evaluation = system.evaluate(optimum.control);
    const auto unit_trajectory = system.forward({1.0});
    const double continuous_fit = analytic_fit(data);
    const double discrete_fit_closed = [&]() {
        double numerator = 0.0;
        double denominator = 0.0;
        for (std::size_t j = 0; j < sample_times.size(); ++j) {
            const double sensitivity = unit_trajectory.state[sample_node(grid, sample_times[j])][0];
            numerator += sensitivity * data.value[j];
            denominator += sensitivity * sensitivity;
        }
        return numerator / denominator;
    }();
    checks.check(std::string(scenario.name) + " discrete optimum matches closed-form discrete LS",
                 std::abs(optimum.control[0] - discrete_fit_closed), 2e-8);
    checks.check(std::string(scenario.name) + " discrete adjoint gradient stationarity",
                 std::abs(evaluation.gradient[0]), 2e-8);
    checks.check(std::string(scenario.name) + " analytic center lies in measurement set",
                 std::max({0.0, measurement_interval(data, scenario.error_bound).first - continuous_fit,
                           continuous_fit - measurement_interval(data, scenario.error_bound).second}),
                 1e-12);

    const auto measured_set = measurement_interval(data, scenario.error_bound);
    checks.check(std::string(scenario.name) + " true initial state lies in measurement set",
                 std::max({0.0, measured_set.first - true_initial,
                           true_initial - measured_set.second}), 1e-12);
    const auto budget = compute_budget(grid, measured_set, continuous_fit,
                                       optimum.control[0], unit_trajectory);
    const auto auc_unit = auc_remaining_unit(grid, unit_trajectory);
    const System<1, 1> auc_system(auc_problem(), grid);
    const auto fitted_path = system.forward(optimum.control);
    const auto auc_evaluation = auc_system.evaluate(optimum.control);
    const double exact_auc_gradient =
        (1.0 - std::exp(-k_elimination * (tf - t0))) / k_elimination;
    checks.check(std::string(scenario.name) + " theta AUC adjoint equals quadrature sensitivity",
                 std::abs(auc_evaluation.gradient[0] -
                          auc_system.cost({1.0}, auc_system.forward({1.0}))), 2e-11);
    checks.check(std::string(scenario.name) + " AUC gradient converges to analytical sensitivity",
                 std::abs(auc_evaluation.gradient[0] - exact_auc_gradient),
                 scenario.steps == 512 ? 2e-5 : 0.02);
    checks.check(std::string(scenario.name) + " numerical concentration allowance encloses truth",
                 [&]() {
                     double violation = 0.0;
                     for (std::size_t n = 0; n < grid.steps(); ++n) {
                         const double left = grid.time[n];
                         const double h = grid.time[n + 1] - left;
                         const double error = linear_exponential_maximum(
                             left, h, fitted_path.state[n][0], fitted_path.state[n + 1][0],
                             true_initial);
                         const double measured = std::max(
                             std::abs(measured_set.first - continuous_fit),
                             std::abs(measured_set.second - continuous_fit)) *
                             std::exp(-k_elimination * left);
                         const double numerical = budget.numerical_concentration[n];
                         violation = std::max(violation, error - measured - numerical - 1e-11);
                     }
                     return std::max(0.0, violation);
                 }(), 2e-10);
    checks.check(std::string(scenario.name) + " measurement and numerical AUC allowances enclose truth",
                 [&]() {
                     double violation = 0.0;
                     for (std::size_t n = 0; n < grid.steps(); ++n) {
                         const double left = grid.time[n];
                         const double h = grid.time[n + 1] - left;
                         const double error = linear_exponential_maximum(
                             left, h, discrete_fit_closed * auc_unit[n],
                             discrete_fit_closed * auc_unit[n + 1],
                             true_initial / k_elimination,
                             -true_initial * std::exp(-k_elimination * tf) / k_elimination);
                         const double measured = std::max(
                             std::abs(measured_set.first - continuous_fit),
                             std::abs(measured_set.second - continuous_fit)) *
                             (std::exp(-k_elimination * left) -
                              std::exp(-k_elimination * tf)) / k_elimination;
                         const double numerical = budget.numerical_remaining_auc[n];
                         violation = std::max(violation, error - measured - numerical - 1e-11);
                     }
                     return std::max(0.0, violation);
                 }(), 2e-10);

    const double joint_time = std::max(*budget.total_concentration_time,
                                       *budget.total_auc_time);
    const bool concentration_is_measurement_dominated =
        budget.measurement_concentration.front() >= budget.numerical_concentration.front();
    const bool auc_is_measurement_dominated =
        budget.measurement_remaining_auc.front() >= budget.numerical_remaining_auc.front();
    decisions << scenario.name << ',' << scenario.theta << ',' << scenario.steps << ','
              << scenario.error_bound << ',' << measured_set.first << ',' << measured_set.second
              << ',' << optimum.control[0] << ',' << continuous_fit << ','
              << concentration_tolerance << ',' << *budget.measurement_concentration_time << ','
              << *budget.numerical_concentration_time << ',' << *budget.total_concentration_time
              << ',' << (concentration_is_measurement_dominated ? "measurement" : "numerical")
              << ',' << remaining_auc_tolerance << ',' << *budget.measurement_auc_time << ','
              << *budget.numerical_auc_time << ',' << *budget.total_auc_time << ','
              << (auc_is_measurement_dominated ? "measurement" : "numerical")
              << ',' << joint_time << '\n';

    std::cout << "\n" << scenario.name << ": theta=" << scenario.theta << ", N="
              << scenario.steps << ", sensor bound=+/-" << scenario.error_bound << " mg/L\n"
              << "  fitted c0: discrete=" << optimum.control[0]
              << ", analytical=" << continuous_fit << " mg/L\n"
              << "  concentration qualification: measurement-only from "
              << *budget.measurement_concentration_time << " h; numerical-only from "
              << *budget.numerical_concentration_time << " h; combined from "
              << *budget.total_concentration_time << " h\n"
              << "  remaining-AUC qualification: measurement-only from "
              << *budget.measurement_auc_time << " h; numerical-only from "
              << *budget.numerical_auc_time << " h; combined from "
              << *budget.total_auc_time << " h\n"
              << "  joint qualification from " << joint_time << " h\n";

    std::vector<double> discrete_remaining_auc(grid.time.size());
    for (std::size_t n = 0; n < grid.time.size(); ++n)
        discrete_remaining_auc[n] = optimum.control[0] * auc_unit[n];
    for (std::size_t n = 0; n < grid.steps(); ++n) {
        budgets << scenario.name << ',' << scenario.theta << ',' << scenario.steps << ','
                << scenario.error_bound << ',' << grid.time[n] << ',' << grid.time[n + 1]
                << ",concentration," << concentration_tolerance << ','
                << budget.measurement_concentration[n] << ','
                << budget.numerical_concentration[n] << ','
                << budget.total_concentration[n] << '\n';
        budgets << scenario.name << ',' << scenario.theta << ',' << scenario.steps << ','
                << scenario.error_bound << ',' << grid.time[n] << ',' << grid.time[n + 1]
                << ",remaining_auc," << remaining_auc_tolerance << ','
                << budget.measurement_remaining_auc[n] << ','
                << budget.numerical_remaining_auc[n] << ','
                << budget.total_remaining_auc[n] << '\n';
    }
    const auto observed_fit = system.forward(optimum.control);
    const double radius = std::max(std::abs(measured_set.first - continuous_fit),
                                   std::abs(measured_set.second - continuous_fit));
    for (std::size_t n = 0; n <= grid.steps(); ++n) {
        const double at = grid.time[n];
        double meas_c = radius * std::exp(-k_elimination * at);
        double num_c = std::abs(continuous_fit * std::exp(-k_elimination * at) -
                                observed_fit.state[n][0]);
        double rem_auc = continuous_fit *
            (std::exp(-k_elimination * at) - std::exp(-k_elimination * tf)) / k_elimination;
        const double meas_auc = radius *
            (std::exp(-k_elimination * at) - std::exp(-k_elimination * tf)) / k_elimination;
        const double num_auc = std::abs(rem_auc - discrete_remaining_auc[n]);
        reconstructions << scenario.name << ',' << scenario.theta << ',' << scenario.steps << ','
                        << scenario.error_bound << ',' << at << ','
                        << exact_concentration(at, true_initial) << ','
                        << exact_concentration(at, continuous_fit) << ','
                        << observed_fit.state[n][0] << ',' << meas_c << ',' << num_c << ','
                        << meas_c + num_c << ',' << exact_remaining_auc(at, true_initial) << ','
                        << rem_auc << ',' << discrete_remaining_auc[n] << ',' << meas_auc << ','
                        << num_auc << ',' << meas_auc + num_auc << '\n';
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        namespace fs = std::filesystem;
        const fs::path output = argc > 1 ? argv[1] : "results/drug_washout";
        fs::create_directories(output);
        Checks checks(output / "checks.csv");
        const Measurement data = observations();

        std::ofstream decisions(output / "decisions.csv");
        std::ofstream budgets(output / "error_budgets.csv");
        std::ofstream reconstructions(output / "reconstruction.csv");
        decisions << std::setprecision(17)
            << "case,theta,steps,sensor_error_bound_mg_L,measurement_set_lower_mg_L,"
            << "measurement_set_upper_mg_L,discrete_fit_mg_L,analytical_fit_mg_L,"
            << "concentration_tolerance_mg_L,measurement_only_time_h,numerical_only_time_h,"
            << "combined_time_h,dominant_initial_budget,remaining_auc_tolerance_mg_h_L,"
            << "measurement_only_auc_time_h,numerical_only_auc_time_h,combined_auc_time_h,"
            << "dominant_initial_auc_budget,joint_time_h\n";
        budgets << std::setprecision(17)
            << "case,theta,steps,sensor_error_bound_mg_L,elapsed_h,interval_end_h,goal,"
            << "tolerance,measurement_component,numerical_component,total_bound\n";
        reconstructions << std::setprecision(17)
            << "case,theta,steps,sensor_error_bound_mg_L,elapsed_h,true_concentration_mg_L,"
            << "analytical_fit_mg_L,discrete_fit_mg_L,measurement_bound_mg_L,"
            << "numerical_bound_mg_L,total_bound_mg_L,true_remaining_auc_mg_h_L,"
            << "analytical_fit_remaining_auc_mg_h_L,discrete_fit_remaining_auc_mg_h_L,"
            << "measurement_auc_bound_mg_h_L,numerical_auc_bound_mg_h_L,"
            << "total_auc_bound_mg_h_L\n";

        std::cout << "DRUG WASHOUT: MEASUREMENT ACCURACY VERSUS NUMERICAL ACCURACY\n"
                  << "One-compartment IV-bolus washout; known elimination rate "
                  << k_elimination << " 1/h; unknown initial concentration.\n"
                  << "The exact same six observations are used for every grid and error bound.\n"
                  << "Theta=0.5 qualification cases; theta=0.65 is included as a first-order convergence control.\n";

        std::cout << "\nTemporal convergence at fixed observations and initial state:\n";
        write_convergence(output, data, checks);

        for (std::size_t i = 0; i < qualification_steps.size(); ++i) {
            for (double sensor_bound : sensor_bounds) {
                const std::string name = std::string("N") +
                    std::to_string(qualification_steps[i]) +
                    (sensor_bound < 0.01 ? "_precise" : "_coarse_sensors");
                run_scenario({qualification_steps[i], sensor_bound, 0.5, name.c_str()},
                             data, checks, decisions, budgets, reconstructions);
            }
        }

        std::cout << "\nChecks passed: " << checks.count
                  << ". Drug-washout measurement/model-error experiment PASS.\n"
                  << "Results: \"" << output.string() << "\"\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
