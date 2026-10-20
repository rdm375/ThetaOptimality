#pragma once

#include <theta_opt/theta_opt.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace stiff_reaction {
using namespace theta_opt;

// A -> B -> C, with total concentration conserved. Time is in minutes;
// concentrations are mol/L. C is eliminated using the known total, leaving
// the two independent dynamic states A and B.
inline constexpr double k_fast = 1000.0; // min^-1
inline constexpr double k_slow = 1.0;    // min^-1
inline constexpr double total = 1.0;     // mol/L
inline constexpr double t0 = 0.0;
inline constexpr double tf = 2.0;
inline constexpr double sensor_sd = 5e-6; // high-precision lab concentration assay
inline constexpr double bounded_noise = 2.0 * sensor_sd;
inline constexpr Vec<2> truth{0.55, 0.25}; // A(t0), B(t0); C(t0)=0.20
inline constexpr Vec<2> prior_lower{0.05, 0.05};
inline constexpr Vec<2> prior_upper{0.80, 0.80};
inline constexpr std::array<double, 9> early_times{
    0.0, 0.0005, 0.001, 0.002, 0.01, 0.1, 0.5, 1.0, 2.0};
inline constexpr std::array<double, 7> late_times{
    0.0, 0.02, 0.1, 0.25, 0.5, 1.0, 2.0};
inline constexpr std::array<double, 9> noise_unit{
    0.20, -0.35, 0.10, 0.40, -0.25, 0.15, -0.30, 0.05, 0.25};
inline constexpr std::array<double, 3> goal_tolerance{0.005, 0.005, 0.005};
inline constexpr std::array<const char*, 3> goal_names{"A", "B", "C"};
inline constexpr Mat<2, 2> A{{{-k_fast, 0.0}, {k_fast, -k_slow}}};
inline constexpr NewtonOptions newton_options{2e-13, 0.0, 8, 30};

inline Mat<2, 2> propagator(double tau) {
    const double ef = std::exp(-k_fast * tau);
    const double es = std::exp(-k_slow * tau);
    const double coupling = k_fast / (k_fast - k_slow);
    return {{{ef, 0.0}, {coupling * (es - ef), es}}};
}

inline Vec<2> exact(double tau, const Vec<2>& xi) {
    return mv(propagator(tau), xi);
}

inline double species(double tau, const Vec<2>& xi, std::size_t which) {
    const auto ab = exact(tau, xi);
    if (which == 0) return ab[0];
    if (which == 1) return ab[1];
    if (which == 2) return total - ab[0] - ab[1];
    throw std::invalid_argument("Unknown species index");
}

inline std::size_t node_at(const Grid& grid, double tau) {
    for (std::size_t n = 0; n < grid.time.size(); ++n)
        if (std::abs(grid.time[n] - tau) < 2e-12) return n;
    throw std::invalid_argument("Reaction observation is not a mesh node");
}

inline Grid mesh(const std::vector<double>& observations, std::size_t subdivisions,
                 double theta) {
    if (!subdivisions || observations.size() < 2)
        throw std::invalid_argument("Invalid reaction mesh");
    std::vector<double> time;
    for (std::size_t j = 0; j + 1 < observations.size(); ++j)
        for (std::size_t k = 0; k < subdivisions; ++k)
            time.push_back(observations[j] + (observations[j + 1] - observations[j]) *
                           static_cast<double>(k) / static_cast<double>(subdivisions));
    time.push_back(observations.back());
    return Grid(time, std::vector<double>(time.size() - 1, theta));
}

inline ScalarFunction<2, 2> product_observation(double observed) {
    auto cost = ScalarFunction<2, 2>::zero();
    cost.evaluate = [=](double, const Vec<2>& x, const Vec<2>&) {
        ScalarValue<2, 2> result;
        const double residual = total - x[0] - x[1] - observed;
        result.value = 0.5 * residual * residual / (sensor_sd * sensor_sd);
        result.x = {-residual / (sensor_sd * sensor_sd),
                    -residual / (sensor_sd * sensor_sd)};
        return result;
    };
    cost.directional_gradient = [](double, const Vec<2>&, const Vec<2>&,
                                   const Vec<2>& direction, const Vec<2>&) {
        ScalarDirection<2, 2> result;
        const double value = (direction[0] + direction[1]) / (sensor_sd * sensor_sd);
        result.x = {value, value};
        return result;
    };
    return cost;
}

inline Problem<2, 2> problem(const Grid& grid, const std::vector<double>& observations) {
    Problem<2, 2> p;
    p.initial.value = [](const Vec<2>& xi) { return xi; };
    p.initial.jacobian = [](const Vec<2>&) { return identity<2>(); };
    p.initial.directional_jacobian = [](const Vec<2>&, const Vec<2>&) { return Mat<2, 2>{}; };
    p.dynamics.value = [](double, const Vec<2>& x, const Vec<2>&) { return mv(A, x); };
    p.dynamics.derivatives = [](double, const Vec<2>&, const Vec<2>&) {
        Derivatives<2, 2> d;
        d.x = A;
        return d;
    };
    p.dynamics.directional_derivatives = [](double, const Vec<2>&, const Vec<2>&,
                                             const Vec<2>&, const Vec<2>&) {
        return Derivatives<2, 2>{};
    };
    for (std::size_t j = 0; j < observations.size(); ++j) {
        const double observed = species(observations[j], truth, 2) + sensor_sd * noise_unit[j];
        p.nodal.emplace_back(node_at(grid, observations[j]), product_observation(observed));
    }
    return p;
}

inline System<2, 2> system(const Grid& grid, const std::vector<double>& observations) {
    return System<2, 2>(problem(grid, observations), grid, newton_options);
}

inline Mat<2, 2> information(const std::vector<double>& observations) {
    Mat<2, 2> hessian{};
    for (double time : observations) {
        const auto row = mv(transpose(propagator(time)), Vec<2>{-1.0, -1.0});
        hessian = add(hessian, scaled(1.0 / (sensor_sd * sensor_sd),
                                      Mat<2, 2>{{{row[0] * row[0], row[0] * row[1]},
                                                 {row[1] * row[0], row[1] * row[1]}}}));
    }
    return hessian;
}

inline Vec<2> analytical_fit(const std::vector<double>& observations) {
    Mat<2, 2> hessian{};
    Vec<2> rhs{};
    for (std::size_t j = 0; j < observations.size(); ++j) {
        const auto p = propagator(observations[j]);
        const Vec<2> row{p[0][0] + p[1][0], p[0][1] + p[1][1]};
        const double target = total - (species(observations[j], truth, 2) +
                              sensor_sd * noise_unit[j]);
        hessian = add(hessian, scaled(1.0 / (sensor_sd * sensor_sd),
                                       Mat<2, 2>{{{row[0] * row[0], row[0] * row[1]},
                                                  {row[1] * row[0], row[1] * row[1]}}}));
        rhs = rhs + (target / (sensor_sd * sensor_sd)) * row;
    }
    return LU<2>(hessian).solve(rhs);
}

inline std::vector<Vec<2>> prior_polygon() {
    auto polygon = rectangle_polygon(prior_lower, prior_upper);
    return clip_halfspace(std::move(polygon), Vec<2>{1.0, 1.0}, 0.9);
}

inline double goal_coeff_fast(std::size_t species_index, const Vec<2>& delta) {
    if (species_index == 0) return delta[0];
    if (species_index == 1) return -k_fast / (k_fast - k_slow) * delta[0];
    if (species_index == 2) return -delta[0] + k_fast / (k_fast - k_slow) * delta[0];
    throw std::invalid_argument("Unknown species index");
}
inline double goal_coeff_slow(std::size_t species_index, const Vec<2>& delta) {
    if (species_index == 0) return 0.0;
    if (species_index == 1)
        return delta[1] + k_fast / (k_fast - k_slow) * delta[0];
    if (species_index == 2) return -delta[1] - k_fast / (k_fast - k_slow) * delta[0];
    throw std::invalid_argument("Unknown species index");
}

inline double max_abs_goal_difference(double begin, double end, std::size_t which,
                                      const Vec<2>& vertex, const Vec<2>& center) {
    const Vec<2> delta = vertex - center;
    const double cf = goal_coeff_fast(which, delta);
    const double cs = goal_coeff_slow(which, delta);
    auto value = [&](double t) { return cf * std::exp(-k_fast * t) + cs * std::exp(-k_slow * t); };
    double maximum = std::max(std::abs(value(begin)), std::abs(value(end)));
    if (cf != 0.0) {
        const double ratio = -k_slow * cs / (k_fast * cf);
        if (ratio > 0.0) {
            const double critical = -std::log(ratio) / (k_fast - k_slow);
            if (critical > begin && critical < end)
                maximum = std::max(maximum, std::abs(value(critical)));
        }
    }
    return maximum;
}
} // namespace stiff_reaction
