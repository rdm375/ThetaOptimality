#include "stiff_reaction_model.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace theta_opt;
namespace {
using namespace stiff_reaction;

struct Checks {
    std::ofstream out;
    std::size_t count{};
    explicit Checks(const std::filesystem::path& path) : out(path) {
        if (!out) throw std::runtime_error("Cannot write stiff-reaction checks");
        out << std::setprecision(17) << "name,error,tolerance,passed\n";
    }
    void require(const std::string& name, double error, double tolerance) {
        const bool passed = std::isfinite(error) && error <= tolerance;
        out << name << ',' << error << ',' << tolerance << ',' << (passed ? 1 : 0) << '\n';
        ++count;
        if (!passed) throw std::runtime_error("Stiff-reaction validation failed: " + name +
                                             " error=" + std::to_string(error));
    }
};

std::vector<double> as_vector(const auto& times) { return {times.begin(), times.end()}; }

double observed_product(double time, std::size_t index) {
    return species(time, truth, 2) + sensor_sd * noise_unit[index];
}

Vec<2> exact_gradient(const std::vector<double>& times, const Vec<2>& xi) {
    Vec<2> gradient{};
    for (std::size_t j = 0; j < times.size(); ++j) {
        const auto p = propagator(times[j]);
        const Vec<2> sensitivity{-(p[0][0] + p[1][0]), -(p[0][1] + p[1][1])};
        const double residual = species(times[j], xi, 2) - observed_product(times[j], j);
        gradient = gradient + (residual / (sensor_sd * sensor_sd)) * sensitivity;
    }
    return gradient;
}

double exact_cost(const std::vector<double>& times, const Vec<2>& xi) {
    double value = 0.0;
    for (std::size_t j = 0; j < times.size(); ++j) {
        const double residual = species(times[j], xi, 2) - observed_product(times[j], j);
        value += 0.5 * residual * residual / (sensor_sd * sensor_sd);
    }
    return value;
}

double matrix_error(const Mat<2, 2>& lhs, const Mat<2, 2>& rhs) {
    return std::max(norm(lhs[0] - rhs[0]), norm(lhs[1] - rhs[1]));
}

double information_condition(const PsdSpectrum2& spectrum) {
    return spectrum.rank == 2 ? spectrum.eigenvalues[0] / spectrum.eigenvalues[1]
                              : std::numeric_limits<double>::infinity();
}

Vec<2> feasible_reference(const ConvexPolygon2& polygon) {
    Vec<2> center{};
    for (const auto& point : polygon) center = center + point;
    return (1.0 / static_cast<double>(polygon.size())) * center;
}

ConvexPolygon2 discrete_confidence_set(const System<2, 2>& sys,
                                       const std::vector<double>& times,
                                       const std::string& label,
                                       Checks& checks) {
    const auto& grid = sys.grid();
    const auto zero = sys.forward({0.0, 0.0});
    const auto e0 = sys.forward({1.0, 0.0});
    const auto e1 = sys.forward({0.0, 1.0});
    const auto prior = prior_polygon();
    auto feasible = prior;
    for (std::size_t j = 0; j < times.size(); ++j) {
        const auto n = node_at(grid, times[j]);
        auto product = [](const Vec<2>& ab) { return total - ab[0] - ab[1]; };
        const double offset = product(zero.state[n]);
        const Vec<2> sensitivity{product(e0.state[n]) - offset,
                                 product(e1.state[n]) - offset};
        const double observed = observed_product(times[j], j);
        double model_allowance = 0.0;
        for (const auto& vertex : prior) {
            const double numerical = product(sys.forward(vertex).state[n]);
            model_allowance = std::max(model_allowance,
                std::abs(numerical - species(times[j], vertex, 2)));
        }
        const double radius = bounded_noise + model_allowance;
        feasible = clip_halfspace(std::move(feasible), sensitivity, observed + radius - offset);
        feasible = clip_halfspace(std::move(feasible), -1.0 * sensitivity,
                                  -observed + radius + offset);
    }
    checks.require(label + " confidence set retains truth",
                   polygon_contains(feasible, truth, 1e-9) ? 0.0 : 1.0, 0.0);
    if (feasible.size() < 3) throw std::runtime_error(label + " confidence set is degenerate");
    return feasible;
}

struct Qualification { IntervalQualification<3> intervals; std::vector<std::array<double, 3>> bounds; };

Qualification qualify(const ConvexPolygon2& feasible, const Vec<2>& center,
                      const std::vector<double>& sample_times) {
    Qualification result;
    std::vector<double> edges;
    constexpr std::size_t per_measurement_interval = 32;
    for (std::size_t j = 0; j + 1 < sample_times.size(); ++j)
        for (std::size_t k = 0; k < per_measurement_interval; ++k)
            edges.push_back(sample_times[j] + (sample_times[j + 1] - sample_times[j]) *
                            static_cast<double>(k) / per_measurement_interval);
    edges.push_back(sample_times.back());
    result.bounds.resize(edges.size() - 1);
    for (std::size_t n = 0; n + 1 < edges.size(); ++n)
        for (std::size_t goal = 0; goal < 3; ++goal)
            for (const auto& vertex : feasible)
                result.bounds[n][goal] = std::max(result.bounds[n][goal],
                    max_abs_goal_difference(edges[n], edges[n + 1], goal, vertex, center));
    result.intervals = qualify_intervals<3>(edges, result.bounds, goal_tolerance);
    return result;
}

void run_stiffness_diagnostic(Checks& checks, std::ofstream& diagnostics,
                              std::ofstream& diagnostic_path) {
    const double h = 0.01;
    const double end = 0.1;
    const std::array<double, 3> theta_values{0.5, 0.65, 1.0};
    const std::size_t steps = static_cast<std::size_t>(std::round(end / h));
    diagnostics << "theta,step_min,kf_h,fast_mode_factor,exact_fast_mode_factor,minimum_species_mol_L,negative_species_nodes,maximum_mass_defect_mol_L\n";
    diagnostic_path << "theta,elapsed_min,A_mol_L,B_mol_L,C_mol_L,exact_A_mol_L,exact_B_mol_L,exact_C_mol_L\n";
    for (double theta : theta_values) {
        const double z = -k_fast * h;
        const double factor = (1.0 + (1.0 - theta) * z) / (1.0 - theta * z);
        const auto grid = Grid::uniform(0.0, end, steps, theta);
        auto p = Problem<2, 2>{};
        p.initial.value = [](const Vec<2>& xi) { return xi; };
        p.initial.jacobian = [](const Vec<2>&) { return identity<2>(); };
        p.initial.directional_jacobian = [](const Vec<2>&, const Vec<2>&) { return Mat<2, 2>{}; };
        p.dynamics.value = [](double, const Vec<2>& x, const Vec<2>&) { return mv(A, x); };
        p.dynamics.derivatives = [](double, const Vec<2>&, const Vec<2>&) {
            Derivatives<2, 2> d; d.x = A; return d;
        };
        p.dynamics.directional_derivatives = [](double, const Vec<2>&, const Vec<2>&,
                                                 const Vec<2>&, const Vec<2>&) {
            return Derivatives<2, 2>{};
        };
        const System<2, 2> system(std::move(p), grid, newton_options);
        const auto path = system.forward(truth);
        double minimum = total;
        std::size_t negatives = 0;
        double mass_defect = 0.0;
        for (std::size_t n = 0; n < path.state.size(); ++n) {
            const double tau = grid.time[n];
            const std::array<double, 3> discrete{path.state[n][0], path.state[n][1],
                                                total - path.state[n][0] - path.state[n][1]};
            const auto exact_ab = exact(tau, truth);
            const std::array<double, 3> exact_all{exact_ab[0], exact_ab[1],
                                                 total - exact_ab[0] - exact_ab[1]};
            for (double v : discrete) {
                minimum = std::min(minimum, v);
                if (v < 0.0) ++negatives;
            }
            mass_defect = std::max(mass_defect,
                std::abs(discrete[0] + discrete[1] + discrete[2] - total));
            diagnostic_path << theta << ',' << tau;
            for (double v : discrete) diagnostic_path << ',' << v;
            for (double v : exact_all) diagnostic_path << ',' << v;
            diagnostic_path << '\n';
        }
        const double exact_factor = std::exp(-k_fast * h);
        diagnostics << theta << ',' << h << ',' << k_fast * h << ',' << factor << ','
                    << exact_factor << ',' << minimum << ',' << negatives << ',' << mass_defect << '\n';
        checks.require("theta stiffness diagnostic mass balance theta=" + std::to_string(theta),
                       mass_defect, 2e-15);
        if (theta < 1.0)
            checks.require("unresolved stiff step exposes negative species theta=" + std::to_string(theta),
                           minimum < 0.0 ? 0.0 : 1.0, 0.0);
        if (theta == 1.0)
            checks.require("backward Euler coarse diagnostic remains nonnegative",
                           minimum >= -1e-14 ? 0.0 : 1.0, 0.0);
    }
}

void run_schedule(const std::string& name, const std::vector<double>& times,
                  double theta, std::size_t subdivisions, Checks& checks,
                  std::ofstream& summary, std::ofstream& observations,
                  std::ofstream& trajectories, std::ofstream& qualification) {
    const auto grid = mesh(times, subdivisions, theta);
    const auto sys = system(grid, times);
    const auto optimum = minimize(sys, Vec<2>{0.5, 0.3}, OptimizeOptions{1e-3, 20, 40, 1e-6});
    const auto eval = sys.evaluate(optimum.control);
    const auto exact_optimum = analytical_fit(times);
    const auto information_matrix = information(times);
    const auto spectrum = psd_spectrum_2d(information_matrix);
    const auto set = discrete_confidence_set(sys, times, name, checks);
    const auto center = feasible_reference(set);
    checks.require(name + " feasible-polygon reference lies in confidence set",
                   polygon_contains(set, center, 1e-10) ? 0.0 : 1.0, 0.0);
    const auto decisions = qualify(set, center, times);
    const std::array<Vec<2>, 3> initial_goals{{{1.0, 0.0}, {0.0, 1.0}, {-1.0, -1.0}}};
    std::array<double, 3> initial_bounds{};
    for (std::size_t j = 0; j < initial_goals.size(); ++j)
        initial_bounds[j] = polygon_goal_bound(set, initial_goals[j], center);

    const Vec<2> direction{0.3, -0.7};
    constexpr double fd_step = 1e-2; // linear dynamics and quadratic cost make a wider difference stable
    const auto tangent = sys.tangent(optimum.control, eval.trajectory, direction);
    const auto plus_path = sys.forward(optimum.control + fd_step * direction);
    const auto minus_path = sys.forward(optimum.control - fd_step * direction);
    double tangent_fd_error = 0.0;
    for (std::size_t n = 0; n < grid.time.size(); ++n) {
        const auto finite_difference = (1.0 / (2.0 * fd_step)) *
                                       (plus_path.state[n] - minus_path.state[n]);
        tangent_fd_error = std::max(tangent_fd_error, norm(tangent[n] - finite_difference));
    }
    checks.require(name + " tangent against state finite differences", tangent_fd_error, 2e-8);
    std::vector<Vec<2>> observation_sources(grid.time.size());
    double forward_observation_pairing = 0.0;
    for (std::size_t j = 0; j < times.size(); ++j) {
        const std::size_t node = node_at(grid, times[j]);
        const double weight = 0.2 + 0.1 * static_cast<double>(j);
        observation_sources[node] = observation_sources[node] + weight * Vec<2>{-1.0, -1.0};
        forward_observation_pairing += dot(weight * Vec<2>{-1.0, -1.0}, tangent[node]);
    }
    const auto observation_pullback = sys.state_pullback(optimum.control, eval.trajectory,
                                                          std::move(observation_sources));
    checks.require(name + " tangent-adjoint observation pairing",
                   std::abs(forward_observation_pairing - dot(observation_pullback, direction)), 2e-8);
    const Vec<2> gradient_probe{0.70, 0.15};
    const auto probe_eval = sys.evaluate(gradient_probe);
    const auto probe_plus_path = sys.forward(gradient_probe + fd_step * direction);
    const auto probe_minus_path = sys.forward(gradient_probe - fd_step * direction);
    const double directional_fd = (sys.cost(gradient_probe + fd_step * direction, probe_plus_path) -
                                   sys.cost(gradient_probe - fd_step * direction, probe_minus_path)) /
                                  (2.0 * fd_step);
    checks.require(name + " discrete adjoint gradient directional difference",
                   std::abs(directional_fd - dot(probe_eval.gradient, direction)) /
                       std::max(1.0, std::abs(directional_fd)), 2e-6);
    const auto hv = sys.hessian_vector(optimum.control, direction);
    const auto plus_eval = sys.evaluate(optimum.control + fd_step * direction);
    const auto minus_eval = sys.evaluate(optimum.control - fd_step * direction);
    const auto hv_fd = (1.0 / (2.0 * fd_step)) * (plus_eval.gradient - minus_eval.gradient);
    checks.require(name + " incremental-adjoint Hessian finite difference",
                   norm(hv - hv_fd) / std::max(1.0, norm(hv)), 2e-6);
    std::array<Vec<2>, 2> hessian_columns{};
    for (std::size_t j = 0; j < 2; ++j) {
        Vec<2> basis{}; basis[j] = 1.0;
        hessian_columns[j] = sys.hessian_vector(optimum.control, basis);
    }
    const Mat<2, 2> discrete_hessian{{{hessian_columns[0][0], hessian_columns[1][0]},
                                      {hessian_columns[0][1], hessian_columns[1][1]}}};
    checks.require(name + " discrete Hessian symmetry",
                   std::abs(discrete_hessian[0][1] - discrete_hessian[1][0]) /
                       std::max(1.0, norm(discrete_hessian[0])), 1e-12);
    for (std::size_t j = 0; j < 2; ++j) {
        Vec<2> basis{}; basis[j] = 1.0;
        checks.require(name + " Hessian versus analytical information column " + std::to_string(j),
            norm(hessian_columns[j] - mv(information_matrix, basis)) /
                std::max(1.0, norm(mv(information_matrix, basis))), 0.08);
    }
    checks.require(name + " initial A is not qualified at t0",
                   initial_bounds[0] > goal_tolerance[0] ? 0.0 : 1.0, 0.0);
    checks.require(name + " initial B is not qualified at t0",
                   initial_bounds[1] > goal_tolerance[1] ? 0.0 : 1.0, 0.0);
    checks.require(name + " A and B eventually qualify after fast transient",
                   decisions.intervals.persistent[0] && decisions.intervals.persistent[1] &&
                   *decisions.intervals.persistent[0] > 0.0 &&
                   *decisions.intervals.persistent[1] > 0.0 ? 0.0 : 1.0, 0.0);

    checks.require(name + " discrete gradient stationarity", norm(eval.gradient), 2e-4);
    checks.require(name + " discrete estimate versus analytical estimate",
                   norm(optimum.control - exact_optimum), 0.08);
    checks.require(name + " exact data are inside declared prior",
                   truth[0] >= prior_lower[0] && truth[0] <= prior_upper[0] &&
                   truth[1] >= prior_lower[1] && truth[1] <= prior_upper[1] &&
                   truth[0] + truth[1] <= 0.9 ? 0.0 : 1.0, 0.0);

    const double t0 = grid.time.front(), tf = grid.time.back();
    summary << name << ',' << theta << ',' << subdivisions << ',' << times.size() << ','
            << spectrum.rank << ',' << information_condition(spectrum) << ',' << eval.cost << ','
            << optimum.control[0] << ',' << optimum.control[1] << ','
            << total - optimum.control[0] - optimum.control[1] << ','
            << exact_optimum[0] << ',' << exact_optimum[1] << ',' << set.size() << ','
            << initial_bounds[0] << ',' << initial_bounds[1] << ',' << initial_bounds[2] << ','
            << eval.gradient[0] << ',' << eval.gradient[1] << '\n';
    for (std::size_t j = 0; j < times.size(); ++j)
        observations << name << ',' << times[j] << ',' << observed_product(times[j], j) << ','
                     << species(times[j], truth, 2) << ',' << sensor_sd << '\n';
    const auto fitted_path = sys.forward(optimum.control);
    for (std::size_t n = 0; n < grid.time.size(); ++n) {
        const double tau = grid.time[n] - t0;
        const auto exact_ab = exact(tau, optimum.control);
        trajectories << name << ',' << grid.time[n] << ','
                     << fitted_path.state[n][0] << ',' << fitted_path.state[n][1] << ','
                     << total - fitted_path.state[n][0] - fitted_path.state[n][1] << ','
                     << species(tau, optimum.control, 0) << ','
                     << species(tau, optimum.control, 1) << ','
                     << species(tau, optimum.control, 2) << ','
                     << fitted_path.state[n][0] - exact_ab[0] << ','
                     << fitted_path.state[n][1] - exact_ab[1] << '\n';
    }
    for (std::size_t j = 0; j < 3; ++j) {
        const auto time = decisions.intervals.persistent[j];
        qualification << name << ',' << goal_names[j] << ',' << goal_tolerance[j] << ','
                      << (time ? *time : std::numeric_limits<double>::quiet_NaN()) << ','
                      << (time ? "qualified" : "unqualified") << ',' << set.size() << '\n';
        std::cout << "  " << name << " / species " << goal_names[j] << " (tolerance "
                  << goal_tolerance[j] << " mol/L): ";
        if (time)
            std::cout << "QUALIFIED persistently from " << std::fixed << std::setprecision(8)
                      << *time << " min through " << (tf - t0) << " min.\n";
        else
            std::cout << "NOT QUALIFIED over the requested horizon.\n";
    }
}

void convergence_study(Checks& checks, std::ofstream& convergence,
                       std::ofstream& orders) {
    const auto times = as_vector(early_times);
    const std::array<std::size_t, 11> refinements{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    const std::array<const char*, 8> names{"state", "tangent", "goal_adjoint", "initial_fit",
        "objective", "gradient", "hessian", "fitted_final_state"};
    for (double theta : {0.5, 0.65}) {
        std::array<double, 8> previous{}, current{};
        double previous_h = 0.0;
        bool first = true;
        for (std::size_t sub : refinements) {
            const auto grid = mesh(times, sub, theta);
            const auto sys = system(grid, times);
            const auto fit = minimize(sys, Vec<2>{0.5, 0.3}, OptimizeOptions{1e-3, 20, 40, 1e-6});
            const auto eval = sys.evaluate(fit.control);
            const auto exact_fit = analytical_fit(times);
            const auto exact_path = sys.forward(truth);
            const auto s0 = sys.tangent(truth, exact_path, {1.0, 0.0});
            const auto s1 = sys.tangent(truth, exact_path, {0.0, 1.0});
            double state_error = 0.0, tangent_error = 0.0, adjoint_error = 0.0;
            for (std::size_t n = 0; n < grid.time.size(); ++n) {
                const auto exact_ab = exact(grid.time[n], truth);
                state_error = std::max(state_error, norm(exact_path.state[n] - exact_ab));
                const auto exact_p = propagator(grid.time[n]);
                tangent_error = std::max({tangent_error,
                    norm(s0[n] - Vec<2>{exact_p[0][0], exact_p[1][0]}),
                    norm(s1[n] - Vec<2>{exact_p[0][1], exact_p[1][1]})});
                if (n == 0 || n + 1 == grid.time.size() || n == grid.time.size() / 2) {
                    const auto adj = sys.goal_gradient(truth, exact_path, n, {1.0, 0.0});
                    const auto exact_adjoint = mv(transpose(exact_p), Vec<2>{1.0, 0.0});
                    adjoint_error = std::max(adjoint_error, norm(adj - exact_adjoint));
                }
            }
            const auto hessian = information(times);
            const auto discrete_hessian = Mat<2, 2>{{
                {sys.hessian_vector(fit.control, {1.0, 0.0})[0], sys.hessian_vector(fit.control, {0.0, 1.0})[0]},
                {sys.hessian_vector(fit.control, {1.0, 0.0})[1], sys.hessian_vector(fit.control, {0.0, 1.0})[1]}}};
            const Vec<2> objective_probe{0.70, 0.15};
            const double objective_error = std::abs(sys.cost(objective_probe,
                                                    sys.forward(objective_probe)) -
                                                    exact_cost(times, objective_probe));
            const double gradient_error = norm(eval.gradient - exact_gradient(times, fit.control));
            const double final_error = norm(sys.forward(fit.control).state.back() - exact(tf - t0, exact_fit));
            current = {state_error, tangent_error, adjoint_error, norm(fit.control - exact_fit),
                       objective_error, gradient_error, matrix_error(discrete_hessian, hessian), final_error};
            double h = 0.0;
            for (std::size_t n = 0; n < grid.steps(); ++n)
                h = std::max(h, grid.time[n + 1] - grid.time[n]);
            convergence << theta << ',' << sub << ',' << h;
            for (double error : current) convergence << ',' << error;
            convergence << '\n';
            checks.require("convergence optimizer stationarity theta=" + std::to_string(theta) +
                           " sub=" + std::to_string(sub), norm(eval.gradient), 1e-3);
            if (!first) {
                const double ratio = previous_h / h;
                std::array<double, 8> observed{};
                for (std::size_t j = 0; j < observed.size(); ++j) {
                    observed[j] = std::log(previous[j] / current[j]) / std::log(ratio);
                    orders << theta << ',' << names[j] << ',' << observed[j] << ','
                           << (theta == 0.5 ? 2.0 : 1.0) << '\n';
                }
                if (sub == refinements.back()) {
                    const double expected = theta == 0.5 ? 2.0 : 1.0;
                    for (std::size_t j = 0; j < observed.size(); ++j)
                        checks.require(std::string("theta=") + std::to_string(theta) +
                                       " observed order " + names[j],
                                       std::abs(observed[j] - expected), 0.30);
                    std::cout << "last-pair orders theta=" << theta << ':';
                    for (double order : observed) std::cout << ' ' << order;
                    std::cout << '\n';
                }
            }
            previous = current;
            previous_h = h;
            first = false;
        }
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        namespace fs = std::filesystem;
        const fs::path output = argc > 1 ? argv[1] : "results/stiff_reaction";
        fs::create_directories(output);
        Checks checks(output / "checks.csv");
        std::ofstream summary(output / "summary.csv"), observations(output / "observations.csv"),
            trajectories(output / "trajectories.csv"), diagnostics(output / "stiffness_diagnostic.csv"),
            diagnostic_path(output / "stiffness_trajectories.csv"), convergence(output / "convergence.csv"),
            orders(output / "orders.csv"), qualification(output / "qualification.csv");
        for (auto* stream : {&summary, &observations, &trajectories, &diagnostics,
                             &diagnostic_path, &convergence, &orders, &qualification}) {
            if (!*stream) throw std::runtime_error("Cannot write stiff-reaction outputs");
            *stream << std::setprecision(17);
        }
        summary << "schedule,theta,subdivisions,observation_count,information_rank,information_condition,cost,fit_A,fit_B,fit_C,exact_fit_A,exact_fit_B,exact_confidence_vertices,initial_A_bound,initial_B_bound,initial_C_bound,gradient_A,gradient_B\n";
        observations << "schedule,elapsed_min,observed_C_mol_L,exact_C_truth_mol_L,sensor_sd_mol_L\n";
        trajectories << "schedule,elapsed_min,theta_A,theta_B,theta_C,exact_A,exact_B,exact_C,A_error,B_error\n";
        convergence << "theta,subdivisions,max_step_min,state_error,tangent_error,goal_adjoint_error,initial_fit_error,objective_error,gradient_error,hessian_error,fitted_final_state_error\n";
        orders << "theta,quantity,observed_order,expected_order\n";
        qualification << "schedule,goal,tolerance_mol_L,persistent_elapsed_min,decision,confidence_polygon_vertices\n";

        std::cout << "Stiff reaction: A -> B -> C; k_fast/k_slow=" << k_fast / k_slow
                  << "; C is inferred from the conserved total concentration.\n";
        std::cout << "Observed convergence orders (last refinement pair):\n";
        run_stiffness_diagnostic(checks, diagnostics, diagnostic_path);
        const auto early = as_vector(early_times);
        run_schedule("early_sampled", early, 0.5, 16, checks, summary, observations,
                     trajectories, qualification);
        run_schedule("late_sampled", as_vector(late_times), 0.5, 16, checks, summary,
                     observations, trajectories, qualification);
        const auto early_information = psd_spectrum_2d(information(early));
        const auto late_information = psd_spectrum_2d(information(as_vector(late_times)));
        std::cout << "information condition: early=" << information_condition(early_information)
                  << ", late=" << information_condition(late_information)
                  << "; minimum eigenvalue early=" << early_information.eigenvalues[1]
                  << ", late=" << late_information.eigenvalues[1] << '\n';
        checks.require("early-sampled and late-sampled information matrices are full rank",
                       early_information.rank == 2 && late_information.rank == 2 ? 0.0 : 1.0, 0.0);
        checks.require("fast-transient samples increase the weakest information eigenvalue",
                       early_information.eigenvalues[1] > late_information.eigenvalues[1] ? 0.0 : 1.0, 0.0);
        convergence_study(checks, convergence, orders);
        std::ofstream report(output / "validation.txt");
        report << "stiff_reaction_checks_passed=" << checks.count << '\n'
               << "reaction=A->B->C\nk_fast_min^-1=" << k_fast
               << "\nk_slow_min^-1=" << k_slow
               << "\nall forward, fit, derivative, and qualification calculations ran in C++20\n";
        std::cout << "Stiff reaction state-finding complete: " << checks.count
                  << " checks passed. Results: " << output << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
