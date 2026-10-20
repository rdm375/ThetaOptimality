#include "weak_coupling_model.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace theta_opt;
namespace {
struct Checks {
    std::ofstream out;
    std::size_t count{};
    explicit Checks(const std::filesystem::path& path) : out(path) {
        if (!out) throw std::runtime_error("Cannot write weak-coupling checks");
        out << std::setprecision(17) << "name,error,tolerance,passed\n";
    }
    void require(const std::string& name, double error, double tolerance) {
        const bool passed = std::isfinite(error) && error <= tolerance;
        out << name << ',' << error << ',' << tolerance << ',' << (passed ? 1 : 0) << '\n';
        ++count;
        if (!passed) throw std::runtime_error("Weak-coupling validation failed: " + name +
                                             " error=" + std::to_string(error));
    }
};

double matrix_error(const Mat<2, 2>& a, const Mat<2, 2>& b) {
    return std::max(norm(a[0] - b[0]), norm(a[1] - b[1]));
}

Mat<2, 2> information_matrix(const System<2, 2>& sys, const Vec<2>& xi) {
    Mat<2, 2> result{};
    for (std::size_t j = 0; j < 2; ++j) {
        Vec<2> e{};
        e[j] = 1.0;
        const auto col = sys.hessian_vector(xi, e);
        for (std::size_t i = 0; i < 2; ++i) result[i][j] = col[i];
    }
    return result;
}

Vec<2> analytical_gradient(double contact, const Vec<2>& xi) {
    const auto eq = weak_thermal::equilibrium(contact);
    const auto y = weak_thermal::measurements(contact);
    Vec<2> result{};
    for (std::size_t j = 0; j < weak_thermal::elapsed.size(); ++j) {
        const auto e = weak_thermal::propagator(contact, weak_thermal::elapsed[j]);
        const Vec<2> h{e[1][0], e[1][1]};
        const double residual = eq[1] + dot(h, xi - eq) - y[j];
        result = result + (residual / (weak_thermal::sensor_sd * weak_thermal::sensor_sd)) * h;
    }
    return result;
}

double analytical_cost(double contact, const Vec<2>& xi) {
    const auto eq = weak_thermal::equilibrium(contact);
    const auto y = weak_thermal::measurements(contact);
    double result = 0.0;
    for (std::size_t j = 0; j < weak_thermal::elapsed.size(); ++j) {
        const auto e = weak_thermal::propagator(contact, weak_thermal::elapsed[j]);
        const Vec<2> h{e[1][0], e[1][1]};
        const double residual = eq[1] + dot(h, xi - eq) - y[j];
        result += 0.5 * residual * residual /
                  (weak_thermal::sensor_sd * weak_thermal::sensor_sd);
    }
    return result;
}

ConvexPolygon2 analytical_polygon(double contact, double radius) {
    auto polygon = rectangle_polygon(weak_thermal::prior_lower, weak_thermal::prior_upper);
    const auto eq = weak_thermal::equilibrium(contact);
    const auto y = weak_thermal::measurements(contact);
    for (std::size_t j = 0; j < weak_thermal::elapsed.size(); ++j) {
        const auto e = weak_thermal::propagator(contact, weak_thermal::elapsed[j]);
        const Vec<2> h{e[1][0], e[1][1]};
        const double offset = eq[1] - dot(h, eq);
        polygon = clip_halfspace(std::move(polygon), h, y[j] + radius - offset);
        polygon = clip_halfspace(std::move(polygon), -1.0 * h, -y[j] + radius + offset);
    }
    return polygon;
}

struct PhysicalError {
    std::vector<double> weighted;
    std::vector<double> component[2];
};

PhysicalError numerical_error(const System<2, 2>& sys, double contact) {
    const auto& grid = sys.grid();
    PhysicalError result;
    result.weighted.assign(grid.time.size(), 0.0);
    for (auto& component : result.component) component.assign(grid.time.size(), 0.0);
    std::vector<double> defect(grid.steps(), 0.0);
    const std::array<Vec<2>, 4> corners{{weak_thermal::prior_lower,
        {weak_thermal::prior_upper[0], weak_thermal::prior_lower[1]},
        weak_thermal::prior_upper,
        {weak_thermal::prior_lower[0], weak_thermal::prior_upper[1]}}};
    const auto a = weak_thermal::dynamics_matrix(contact);
    for (const auto& corner : corners) {
        const auto path = sys.forward(corner);
        for (std::size_t n = 0; n < grid.steps(); ++n) {
            const double h = grid.time[n + 1] - grid.time[n];
            const double theta = grid.theta[n];
            const auto delta = path.state[n + 1] - path.state[n];
            const double interpolation_defect = std::max(theta, 1.0 - theta) *
                                                weak_thermal::energy_norm(mv(a, delta));
            const auto f0 = mv(a, path.state[n] - weak_thermal::equilibrium(contact));
            const auto f1 = mv(a, path.state[n + 1] - weak_thermal::equilibrium(contact));
            const auto residual = delta - h * ((1.0 - theta) * f0 + theta * f1);
            const double solve_defect = weak_thermal::energy_norm(residual) / h;
            defect[n] = std::max(defect[n], interpolation_defect + solve_defect);
        }
    }
    const double mu = weak_thermal::slow_rate(contact);
    for (std::size_t n = 0; n < grid.steps(); ++n) {
        const double h = grid.time[n + 1] - grid.time[n];
        result.weighted[n + 1] = propagate_defect_bound(result.weighted[n], defect[n], h, mu);
    }
    for (std::size_t n = 0; n < grid.time.size(); ++n) {
        result.component[0][n] = result.weighted[n] / std::sqrt(weak_thermal::winding_capacity);
        result.component[1][n] = result.weighted[n] / std::sqrt(weak_thermal::housing_capacity);
    }
    return result;
}

ConvexPolygon2 numerical_confidence_set(const System<2, 2>& sys, double contact,
                                        const PhysicalError& error, double radius) {
    const auto& grid = sys.grid();
    const Vec<2> reference{115.0, 55.0};
    const auto base = sys.forward(reference);
    const auto s0 = sys.tangent(reference, base, {1.0, 0.0});
    const auto s1 = sys.tangent(reference, base, {0.0, 1.0});
    auto polygon = rectangle_polygon(weak_thermal::prior_lower, weak_thermal::prior_upper);
    const auto y = weak_thermal::measurements(contact);
    for (std::size_t j = 0; j < weak_thermal::elapsed.size(); ++j) {
        const std::size_t n = weak_thermal::node_at(grid, weak_thermal::elapsed[j]);
        const Vec<2> h{s0[n][1], s1[n][1]};
        const double offset = base.state[n][1] - dot(h, reference);
        const double bound = radius + error.component[1][n];
        polygon = clip_halfspace(std::move(polygon), h, y[j] + bound - offset);
        polygon = clip_halfspace(std::move(polygon), -1.0 * h, -y[j] + bound + offset);
    }
    if (polygon.size() < 3) throw std::runtime_error("Winding/housing confidence set is empty or degenerate");
    return polygon;
}

struct Decision {
    std::array<std::optional<double>, 4> persistent{};
    std::array<std::optional<double>, 4> first{};
    std::optional<double> joint;
    std::vector<std::array<double, 4>> envelopes;
};

Decision qualify(const System<2, 2>& sys, const Vec<2>& fit, const ConvexPolygon2& set,
                 const PhysicalError& error) {
    const auto& grid = sys.grid();
    const auto path = sys.forward(fit);
    Decision result;
    std::vector<std::array<double, 4>> nodal(grid.time.size());
    for (std::size_t n = 0; n < grid.time.size(); ++n) {
        for (std::size_t j = 0; j < weak_thermal::goals.size(); ++j) {
            const auto goal = weak_thermal::goals[j];
            const auto gradient = sys.goal_gradient(fit, path, n, goal);
            double data_bound = 0.0;
            for (const auto& vertex : set)
                data_bound = std::max(data_bound, std::abs(dot(gradient, vertex - fit)));
            const double numerical = weak_thermal::goal_error_factor(goal) * error.weighted[n];
            nodal[n][j] = data_bound + numerical;
        }
    }
    result.envelopes.resize(grid.steps());
    for (std::size_t n = 0; n < grid.steps(); ++n) {
        for (std::size_t j = 0; j < weak_thermal::goals.size(); ++j)
            result.envelopes[n][j] = std::max(nodal[n][j], nodal[n + 1][j]);
    }
    const auto intervals = qualify_intervals<4>(grid.time, result.envelopes, weak_thermal::goal_tolerance);
    result.persistent = intervals.persistent;
    result.first = intervals.first_below;
    result.joint = intervals.joint;
    for (auto& t : result.persistent) if (t) *t -= weak_thermal::t0;
    for (auto& t : result.first) if (t) *t -= weak_thermal::t0;
    if (result.joint) *result.joint -= weak_thermal::t0;
    return result;
}

std::string label(double contact, const std::string& suffix) {
    return "K=" + std::to_string(contact) + " " + suffix;
}

std::string case_name(double contact) {
    for (std::size_t j = 0; j < weak_thermal::conductances.size(); ++j)
        if (contact == weak_thermal::conductances[j]) return weak_thermal::case_names[j];
    throw std::invalid_argument("Unknown conductance case");
}

double val(const std::optional<double>& x) {
    return x.value_or(std::numeric_limits<double>::quiet_NaN());
}

void run_case(double contact, Checks& checks,
              std::ofstream& summary, std::ofstream& decisions, std::ofstream& envelopes,
              std::ofstream& trajectories, std::ofstream& observations, std::ofstream& adjoints) {
    const auto grid = weak_thermal::mesh(256, 0.5);
    const auto sys = weak_thermal::system(grid, contact);
    const auto reference = weak_thermal::analytical_fit(contact);
    const auto optimum = weak_thermal::fit(sys);
    const auto eval = sys.evaluate(optimum.control);
    const auto spectrum = psd_spectrum_2d(reference.information);
    const auto scaled_information = Mat<2, 2>{{
        {reference.information[0][0] * 45.0 * 45.0, reference.information[0][1] * 45.0 * 20.0},
        {reference.information[1][0] * 45.0 * 20.0, reference.information[1][1] * 20.0 * 20.0}}};
    const auto scaled_spectrum = psd_spectrum_2d(scaled_information);
    const std::string prefix = weak_thermal::case_names[contact == weak_thermal::conductances[0] ? 0 :
                                                        contact == weak_thermal::conductances[1] ? 1 : 2];

    checks.require(label(contact, "forward state at t0"), norm(eval.trajectory.state.front() - optimum.control), 1e-13);
    checks.require(label(contact, "discrete optimizer stationarity"), norm(eval.gradient), 2e-8);
    checks.require(label(contact, "discrete fit versus analytical optimum"), norm(optimum.control - reference.initial), 1.5);
    checks.require(label(contact, "Hessian symmetry"), std::abs(reference.information[0][1] - reference.information[1][0]), 2e-9);
    checks.require(label(contact, "information matrix is full rank"), spectrum.rank == 2 ? 0.0 : 1.0, 0.0);
    const Vec<2> direction{0.7, -0.4};
    const double epsilon = 1e-4;
    const auto fd_plus = sys.evaluate(optimum.control + epsilon * direction);
    const auto fd_minus = sys.evaluate(optimum.control - epsilon * direction);
    checks.require(label(contact, "discrete adjoint objective gradient finite difference"),
                   std::abs((fd_plus.cost - fd_minus.cost) / (2.0 * epsilon) - dot(eval.gradient, direction)), 2e-7);

    const auto tangent = sys.tangent(optimum.control, eval.trajectory, direction);
    const auto plus = sys.forward(optimum.control + epsilon * direction);
    const auto minus = sys.forward(optimum.control - epsilon * direction);
    double tangent_fd = 0.0;
    for (std::size_t n = 0; n < grid.time.size(); ++n)
        tangent_fd = std::max(tangent_fd, norm((1.0 / (2.0 * epsilon)) *
                                               (plus.state[n] - minus.state[n]) - tangent[n]));
    checks.require(label(contact, "tangent finite difference"), tangent_fd, 3e-7);

    const auto arbitrary = sys.forward(optimum.control);
    std::vector<Vec<2>> source(grid.time.size());
    for (std::size_t n : {std::size_t{0}, grid.steps() / 3, 2 * grid.steps() / 3, grid.steps()})
        source[n] = source[n] + Vec<2>{0.3 + 0.001 * static_cast<double>(n), -0.2};
    const auto pullback = sys.state_pullback(optimum.control, arbitrary, source);
    double pairing = 0.0;
    for (std::size_t n = 0; n < grid.time.size(); ++n) pairing += dot(source[n], tangent[n]);
    checks.require(label(contact, "tangent-adjoint pairing"), std::abs(pairing - dot(pullback, direction)), 2e-9);

    const auto hessian_direction = sys.hessian_vector(optimum.control, direction);
    const auto eval_plus = sys.evaluate(optimum.control + epsilon * direction);
    const auto eval_minus = sys.evaluate(optimum.control - epsilon * direction);
    checks.require(label(contact, "incremental-adjoint Hessian finite difference"),
                   norm((1.0 / (2.0 * epsilon)) * (eval_plus.gradient - eval_minus.gradient) - hessian_direction), 3e-7);
    const auto discrete_information = information_matrix(sys, optimum.control);
    checks.require(label(contact, "discrete Hessian symmetry"),
                   std::abs(discrete_information[0][1] - discrete_information[1][0]), 2e-9);
    double pullback_error = 0.0;
    for (std::size_t j = 0; j < 2; ++j) {
        Vec<2> unit{};
        unit[j] = 1.0;
        const auto sensitivity = sys.tangent(optimum.control, eval.trajectory, unit);
        std::vector<Vec<2>> sources(grid.time.size());
        for (double tau : weak_thermal::elapsed) {
            const auto n = weak_thermal::node_at(grid, tau);
            sources[n][1] += sensitivity[n][1] / (weak_thermal::sensor_sd * weak_thermal::sensor_sd);
        }
        const auto column = sys.state_pullback(optimum.control, eval.trajectory, std::move(sources));
        const auto hcolumn = sys.hessian_vector(optimum.control, unit);
        pullback_error = std::max(pullback_error, norm(column - hcolumn));
    }
    checks.require(label(contact, "observation tangent-adjoint information action"), pullback_error, 2e-9);

    const auto errors = numerical_error(sys, contact);
    const double radius = independent_gaussian_box_radius(weak_thermal::sensor_sd,
                                                          weak_thermal::elapsed.size(), 0.95);
    const auto feasible = numerical_confidence_set(sys, contact, errors, radius);
    const auto exact_feasible = analytical_polygon(contact, radius);
    checks.require(label(contact, "bounded set contains synthetic truth"),
                   polygon_contains(feasible, weak_thermal::truth) ? 0.0 : 1.0, 0.0);
    double cover_error = 0.0;
    for (const auto& vertex : exact_feasible)
        if (!polygon_contains(feasible, vertex, 2e-7)) cover_error = 1.0;
    checks.require(label(contact, "numerical set contains analytical set vertices"), cover_error, 0.0);
    const auto result = qualify(sys, optimum.control, feasible, errors);

    const auto a = weak_thermal::dynamics_matrix(contact);
    const auto eq = weak_thermal::equilibrium(contact);
    const auto exact_final = weak_thermal::exact(contact, 600.0, weak_thermal::truth);
    const auto fit_final = weak_thermal::exact(contact, 600.0, reference.initial);
    (void)a;
    summary << contact << ',' << prefix << ',' << weak_thermal::slow_rate(contact) << ','
            << weak_thermal::fast_rate(contact) << ','
            << std::log(2.0) / std::abs(weak_thermal::fast_rate(contact)) << ','
            << std::log(2.0) / std::abs(weak_thermal::slow_rate(contact)) << ','
            << spectrum.eigenvalues[0] << ',' << spectrum.eigenvalues[1] << ','
            << spectrum.eigenvalues[0] / spectrum.eigenvalues[1] << ','
            << scaled_spectrum.eigenvalues[0] / scaled_spectrum.eigenvalues[1] << ','
            << reference.initial[0] << ',' << reference.initial[1] << ','
            << optimum.control[0] << ',' << optimum.control[1] << ',' << norm(eval.gradient) << ','
            << exact_final[0] << ',' << exact_final[1] << ',' << fit_final[0] << ',' << fit_final[1] << ','
            << eq[0] << ',' << eq[1] << ',' << feasible.size() << ',' << result.joint.value_or(std::numeric_limits<double>::quiet_NaN()) << '\n';

    for (std::size_t j = 0; j < weak_thermal::goal_names.size(); ++j) {
        decisions << contact << ',' << prefix << ',' << weak_thermal::goal_names[j] << ','
                  << weak_thermal::goal_tolerance[j] << ',' << val(result.first[j]) << ','
                  << val(result.persistent[j]) << ',' << (result.persistent[j] ? 1 : 0) << '\n';
    }
    for (std::size_t n = 0; n < grid.steps(); ++n)
        for (std::size_t j = 0; j < weak_thermal::goal_names.size(); ++j)
            envelopes << contact << ',' << prefix << ',' << grid.time[n] - weak_thermal::t0 << ','
                      << weak_thermal::goal_names[j] << ',' << result.envelopes[n][j] << ','
                      << weak_thermal::goal_tolerance[j] << '\n';

    const auto costate = sys.cost_to_go_adjoint(optimum.control, eval);
    const auto truth_path = sys.forward(weak_thermal::truth);
    const auto exact_fit_path = sys.forward(reference.initial);
    for (std::size_t n = 0; n < grid.time.size(); ++n) {
        const double tau = grid.time[n] - weak_thermal::t0;
        const auto exact_truth = weak_thermal::exact(contact, tau, weak_thermal::truth);
        const auto exact_fit = weak_thermal::exact(contact, tau, reference.initial);
        const auto exact_discrete_fit = weak_thermal::exact(contact, tau, optimum.control);
        const double exact_measurement = std::find(weak_thermal::elapsed.begin(), weak_thermal::elapsed.end(), tau) != weak_thermal::elapsed.end()
            ? exact_truth[1] : std::numeric_limits<double>::quiet_NaN();
        trajectories << contact << ',' << prefix << ',' << tau << ','
                     << exact_truth[0] << ',' << exact_truth[1] << ','
                     << truth_path.state[n][0] << ',' << truth_path.state[n][1] << ','
                     << exact_fit[0] << ',' << exact_fit[1] << ','
                     << exact_fit_path.state[n][0] << ',' << exact_fit_path.state[n][1] << ','
                     << eval.trajectory.state[n][0] << ',' << eval.trajectory.state[n][1] << ','
                     << exact_discrete_fit[0] << ',' << exact_discrete_fit[1] << ','
                     << errors.component[0][n] << ',' << errors.component[1][n] << ',' << exact_measurement << '\n';
        adjoints << contact << ',' << prefix << ',' << tau << ',' << costate[n][0] << ',' << costate[n][1] << '\n';
    }
    const auto y = weak_thermal::measurements(contact);
    for (std::size_t j = 0; j < weak_thermal::elapsed.size(); ++j)
        observations << contact << ',' << prefix << ',' << weak_thermal::elapsed[j] << ','
                     << weak_thermal::exact(contact, weak_thermal::elapsed[j], weak_thermal::truth)[1]
                     << ',' << y[j] << ',' << radius << '\n';

    // Independent physical balance: the exact continuous model satisfies the
    // winding/housing first-law equations at every sampled state.
    double physical_balance = 0.0;
    for (std::size_t n = 0; n < grid.time.size(); n += std::max<std::size_t>(1, grid.steps() / 50)) {
        const auto x = weak_thermal::exact(contact, grid.time[n] - weak_thermal::t0, weak_thermal::truth);
        const auto f = mv(a, x - eq);
        const double winding_balance = weak_thermal::winding_capacity * f[0] -
            (weak_thermal::motor_loss - contact * (x[0] - x[1]));
        const double housing_balance = weak_thermal::housing_capacity * f[1] -
            (contact * (x[0] - x[1]) - weak_thermal::ambient_conductance * (x[1] - weak_thermal::ambient));
        physical_balance = std::max({physical_balance, std::abs(winding_balance), std::abs(housing_balance)});
    }
    checks.require(label(contact, "independent thermal balance"), physical_balance, 2e-11);
}

struct RunningReference { double cost{}; Vec<2> gradient{}; Mat<2, 2> hessian{}; };

RunningReference running_reference(double contact, const Vec<2>& xi) {
    const auto a = weak_thermal::dynamics_matrix(contact);
    const double trace = a[0][0] + a[1][1];
    const double determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    const double delta = std::sqrt(trace * trace - 4.0 * determinant);
    const double slow = 0.5 * (trace + delta), fast = 0.5 * (trace - delta);
    const auto ps = scaled(1.0 / (slow - fast), add(a, scaled(-fast, identity<2>())));
    const auto pf = add(identity<2>(), scaled(-1.0, ps));
    const Vec<2> rs{ps[1][0], ps[1][1]}, rf{pf[1][0], pf[1][1]};
    auto integral = [](double rate) { return std::expm1(rate * (weak_thermal::tf - weak_thermal::t0)) / rate; };
    auto q = add(scaled(integral(2.0 * slow), weak_thermal::outer(rs)),
                 scaled(integral(2.0 * fast), weak_thermal::outer(rf)));
    for (std::size_t i = 0; i < 2; ++i) for (std::size_t j = 0; j < 2; ++j)
        q[i][j] += integral(slow + fast) * (rs[i] * rf[j] + rf[i] * rs[j]);
    const auto d = xi - weak_thermal::equilibrium(contact);
    return {0.5 * dot(d, mv(q, d)), mv(q, d), q};
}

void run_quadrature_test(double contact, Checks& checks, std::ofstream& csv) {
    const auto grid = weak_thermal::mesh(256, 0.5);
    auto p = weak_thermal::problem(grid, contact);
    p.nodal.clear();
    p.running.evaluate = [contact](double, const Vec<2>& x, const Vec<2>&) {
        ScalarValue<2, 2> out;
        const double r = x[1] - weak_thermal::equilibrium(contact)[1];
        out.value = 0.5 * r * r;
        out.x[1] = r;
        return out;
    };
    p.running.directional_gradient = [](double, const Vec<2>&, const Vec<2>&,
                                         const Vec<2>& s, const Vec<2>&) {
        ScalarDirection<2, 2> out;
        out.x[1] = s[1];
        return out;
    };
    System<2, 2> sys(std::move(p), grid, weak_thermal::solver_options);
    const Vec<2> xi{125.0, 62.0};
    const auto evaluated = sys.evaluate(xi);
    const auto exact = running_reference(contact, xi);
    // An independent fine composite trapezoid integrates the exact analytical
    // trajectory, providing a check separate from theta endpoint quadrature.
    constexpr std::size_t count = 2000000;
    double reference = 0.0;
    const double h = (weak_thermal::tf - weak_thermal::t0) / static_cast<double>(count);
    for (std::size_t j = 0; j <= count; ++j) {
        const double t = h * static_cast<double>(j);
        const double r = weak_thermal::exact(contact, t, xi)[1] - weak_thermal::equilibrium(contact)[1];
        reference += (j == 0 || j == count ? 0.5 : 1.0) * 0.5 * r * r * h;
    }
    checks.require(label(contact, "theta running quadrature versus analytic integral"),
                   std::abs(evaluated.cost - exact.cost), 0.1);
    checks.require(label(contact, "running adjoint gradient versus analytic gradient"),
                   norm(evaluated.gradient - exact.gradient), 0.1);
    checks.require(label(contact, "theta running quadrature versus fine independent integral"),
                   std::abs(exact.cost - reference), 5e-6);
    const Vec<2> direction{0.4, -0.7};
    const double eps = 1e-4;
    const auto plus = sys.evaluate(xi + eps * direction), minus = sys.evaluate(xi - eps * direction);
    checks.require(label(contact, "running adjoint gradient finite difference"),
                   std::abs((plus.cost - minus.cost) / (2.0 * eps) - dot(evaluated.gradient, direction)), 2e-6);
    const auto hv = sys.hessian_vector(xi, direction);
    checks.require(label(contact, "running Hessian-vector versus analytic Hessian"),
                   norm(hv - mv(exact.hessian, direction)), 0.1);
    checks.require(label(contact, "running incremental-adjoint finite difference"),
                   norm((1.0 / (2.0 * eps)) * (plus.gradient - minus.gradient) - hv), 2e-6);
    csv << contact << ',' << evaluated.cost << ',' << exact.cost << ',' << reference << ','
        << norm(evaluated.gradient) << ',' << norm(hv) << '\n';
}

void convergence_study(double contact, double theta, Checks& checks, std::ofstream& csv, std::ofstream& order_csv) {
    const auto reference = weak_thermal::analytical_fit(contact);
    const std::array<std::size_t, 11> refinements{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    std::array<double, 8> previous{};
    std::array<double, 8> current{};
    std::array<double, 8> orders{};
    double previous_h = 0.0;
    bool first = true;
    for (const auto sub : refinements) {
        const auto grid = weak_thermal::mesh(sub, theta);
        const auto sys = weak_thermal::system(grid, contact);
        const auto fit = weak_thermal::fit(sys);
        const auto eval = sys.evaluate(fit.control);
        double state_error = 0.0, tangent_error = 0.0, adjoint_error = 0.0;
        const auto state = sys.forward(weak_thermal::truth);
        const auto s0 = sys.tangent(weak_thermal::truth, state, {1.0, 0.0});
        const auto s1 = sys.tangent(weak_thermal::truth, state, {0.0, 1.0});
        for (std::size_t n = 0; n < grid.time.size(); ++n) {
            const double tau = grid.time[n] - weak_thermal::t0;
            const auto exact_state = weak_thermal::exact(contact, tau, weak_thermal::truth);
            const auto e = weak_thermal::propagator(contact, tau);
            state_error = std::max(state_error, norm(state.state[n] - exact_state));
            tangent_error = std::max({tangent_error, norm(s0[n] - Vec<2>{e[0][0], e[1][0]}),
                                      norm(s1[n] - Vec<2>{e[0][1], e[1][1]})});
            if (n == 0 || n == grid.steps() / 3 || n == 2 * grid.steps() / 3 || n == grid.steps()) {
                const auto goal_adj = sys.goal_gradient(weak_thermal::truth, state, n, {1.0, 0.0});
                const auto exact_adj = mv(transpose(e), Vec<2>{1.0, 0.0});
                adjoint_error = std::max(adjoint_error, norm(goal_adj - exact_adj));
            }
        }
        const auto exact_h = reference.information;
        const auto discrete_h = information_matrix(sys, fit.control);
        const double hessian_error = matrix_error(discrete_h, exact_h);
        const double cost_error = std::abs(analytical_cost(contact, weak_thermal::truth) -
                                           sys.cost(weak_thermal::truth, sys.forward(weak_thermal::truth)));
        const double gradient_error = norm(eval.gradient - analytical_gradient(contact, fit.control));
        const double controls_error = norm(fit.control - reference.initial);
        const auto fitted_path = sys.forward(fit.control);
        const double final_state_error = norm(fitted_path.state.back() -
                                              weak_thermal::exact(contact, 600.0, reference.initial));
        current = {state_error, tangent_error, adjoint_error, controls_error, cost_error,
                   gradient_error, hessian_error, final_state_error};
        double max_step = 0.0;
        for (std::size_t n = 0; n < grid.steps(); ++n)
            max_step = std::max(max_step, grid.time[n + 1] - grid.time[n]);
        csv << contact << ',' << case_name(contact) << ',' << theta << ',' << sub << ',' << max_step;
        for (double e : current) csv << ',' << e;
        csv << '\n';
        checks.require(label(contact, "convergence optimizer stationarity"), norm(eval.gradient), 2e-8);
        if (!first) {
            const double ratio_h = previous_h / max_step;
            for (std::size_t j = 0; j < orders.size(); ++j)
                orders[j] = (previous[j] > 0.0 && current[j] > 0.0) ?
                    std::log(previous[j] / current[j]) / std::log(ratio_h) :
                    std::numeric_limits<double>::quiet_NaN();
        }
        previous = current;
        previous_h = max_step;
        first = false;
    }
    const double expected = theta == 0.5 ? 2.0 : 1.0;
    const std::array<const char*, 8> names{{"state", "tangent", "goal_adjoint", "initial_fit",
        "objective", "gradient", "hessian", "fitted_final_state"}};
    for (std::size_t j = 0; j < orders.size(); ++j) {
        order_csv << contact << ',' << theta << ',' << names[j] << ',' << orders[j] << ',' << expected << '\n';
        checks.require(label(contact, std::string("theta=") + std::to_string(theta) + " observed order " + names[j]),
                       std::abs(orders[j] - expected), 0.22);
    }
    std::cout << std::fixed << std::setprecision(4)
              << contact << ',' << theta;
    for (double order : orders) std::cout << ',' << order;
    std::cout << '\n';
    (void)previous;
    (void)previous_h;
}
}

int main(int argc, char** argv) {
    try {
        namespace fs = std::filesystem;
        const fs::path output = argc > 1 ? argv[1] : "results/weak_coupling";
        fs::create_directories(output);
        Checks checks(output / "checks.csv");
        std::ofstream summary(output / "summary.csv"), decisions(output / "qualification_decisions.csv"),
                      envelopes(output / "qualification_envelopes.csv"), trajectories(output / "trajectories.csv"),
                      observations(output / "observations.csv"), adjoints(output / "costates.csv"),
                      convergence(output / "convergence.csv"), orders(output / "orders.csv"),
                      quadrature(output / "quadrature.csv");
        for (auto* stream : {&summary, &decisions, &envelopes, &trajectories, &observations,
                             &adjoints, &convergence, &orders, &quadrature}) {
            if (!*stream) throw std::runtime_error("Cannot write weak-coupling output CSV");
            *stream << std::setprecision(17);
        }
        summary << "contact_kW_per_K,case,slow_rate_per_min,fast_rate_per_min,fast_half_life_min,slow_half_life_min,information_eigenvalue_max,information_eigenvalue_min,information_condition,scaled_information_condition,analytical_winding_fit_C,analytical_housing_fit_C,discrete_winding_fit_C,discrete_housing_fit_C,gradient_inf_norm,truth_final_winding_C,truth_final_housing_C,analytical_fit_final_winding_C,analytical_fit_final_housing_C,equilibrium_winding_C,equilibrium_housing_C,feasible_polygon_vertices,joint_qualification_elapsed_min\n";
        decisions << "contact_kW_per_K,case,goal,tolerance,first_crossing_elapsed_min,persistent_elapsed_min,qualified\n";
        envelopes << "contact_kW_per_K,case,elapsed_min,goal,bound,tolerance\n";
        trajectories << "contact_kW_per_K,case,elapsed_min,truth_exact_winding_C,truth_exact_housing_C,truth_theta_winding_C,truth_theta_housing_C,analytical_fit_winding_C,analytical_fit_housing_C,discrete_trajectory_at_analytical_fit_winding_C,discrete_trajectory_at_analytical_fit_housing_C,optimized_theta_winding_C,optimized_theta_housing_C,optimized_exact_winding_C,optimized_exact_housing_C,winding_numerical_allowance_K,housing_numerical_allowance_K,observation_reference_C\n";
        observations << "contact_kW_per_K,case,elapsed_min,truth_housing_C,observed_housing_C,joint_noise_box_halfwidth_K\n";
        adjoints << "contact_kW_per_K,case,elapsed_min,lambda_winding,lambda_housing\n";
        convergence << "contact_kW_per_K,case,theta,subdivisions,max_step_min,state_error_inf_K,tangent_error_inf,goal_adjoint_error_inf,initial_fit_error_K,objective_error_at_truth,gradient_error,hessian_error,fitted_final_state_error_K\n";
        orders << "contact_kW_per_K,theta,quantity,observed_order,expected_order\n";
        quadrature << "contact_kW_per_K,theta_cost,analytic_integral,fine_exact_integral,gradient_inf_norm,hessian_vector_inf_norm\n";

        std::cout << "Observed convergence orders from the last mesh-refinement pair\n"
                  << "Expected order: theta=0.5 -> 2; theta=0.65 -> 1\n"
                  << "contact_K,theta,state,tangent,goal_adjoint,initial_fit,objective,gradient,hessian,fitted_final_state\n";
        for (double contact : weak_thermal::conductances) {
            run_case(contact, checks, summary, decisions, envelopes, trajectories, observations, adjoints);
            run_quadrature_test(contact, checks, quadrature);
            for (double theta : {0.5, 0.65}) convergence_study(contact, theta, checks, convergence, orders);
        }
        std::ofstream report(output / "validation.txt");
        report << "weak_coupling_checks_passed=" << checks.count << '\n'
               << "cases=3\nconvergence_theta_values=0.5,0.65\n"
               << "all forward, inverse, derivative, balance, and qualification decisions executed in C++20\n";
        std::cout << "Weak-coupling winding/housing qualification complete: " << checks.count << " checks passed.\n"
                  << "Numerical decisions and reference comparisons ran in C++20. Results: " << output << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
