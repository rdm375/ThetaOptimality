// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once

#include <theta_opt/theta_opt.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace weak_thermal {
using namespace theta_opt;

inline constexpr double t0 = 3.0;
inline constexpr double tf = 603.0;
inline constexpr double ambient = 25.0;       // deg C
inline constexpr double winding_capacity = 18.0; // kJ/K
inline constexpr double housing_capacity = 120.0; // kJ/K
inline constexpr double ambient_conductance = 2.0; // kW/K
inline constexpr double motor_loss = 2.0;        // kW, known and constant
inline constexpr double sensor_sd = 0.20;        // K
inline constexpr Vec<2> truth{130.0, 55.0};       // winding, housing; deg C at t0
inline constexpr Vec<2> prior_lower{70.0, 35.0};
inline constexpr Vec<2> prior_upper{160.0, 75.0};
inline constexpr std::array<double, 7> elapsed{0.0, 20.0, 60.0, 120.0, 240.0, 420.0, 600.0};
inline constexpr std::array<double, 7> noise_in_sd{0.10, -0.15, 0.05, 0.20, -0.10, 0.15, -0.05};
inline constexpr std::array<double, 4> goal_tolerance{3.0, 1.0, 4.0, 1000.0};
inline constexpr std::array<Vec<2>, 4> goals{{{1.0, 0.0}, {0.0, 1.0}, {1.0, -1.0},
                                               {winding_capacity, housing_capacity}}};
inline constexpr std::array<const char*, 4> goal_names{{"winding_temperature", "housing_temperature",
                                                        "winding_housing_difference", "stored_heat"}};
inline constexpr NewtonOptions solver_options{2e-13, 0.0, 8, 30};

inline constexpr std::array<double, 3> conductances{0.03, 0.30, 3.0}; // kW/K
inline constexpr std::array<const char*, 3> case_names{{"weak", "moderate", "strong"}};

inline Mat<2, 2> dynamics_matrix(double contact) {
    if (!std::isfinite(contact) || contact <= 0.0) throw std::invalid_argument("Positive contact conductance required");
    return {{{-contact / winding_capacity, contact / winding_capacity},
             {contact / housing_capacity, -(contact + ambient_conductance) / housing_capacity}}};
}

inline Vec<2> equilibrium(double contact) {
    const double housing = ambient + motor_loss / ambient_conductance;
    return {housing + motor_loss / contact, housing};
}

inline Mat<2, 2> propagator(double contact, double tau) {
    const auto a = dynamics_matrix(contact);
    const double trace = a[0][0] + a[1][1];
    const double determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    const double discriminant = std::sqrt(trace * trace - 4.0 * determinant);
    const double slow = 0.5 * (trace + discriminant);
    const double fast = 0.5 * (trace - discriminant);
    Mat<2, 2> pslow = scaled(1.0 / (slow - fast), add(a, scaled(-fast, identity<2>())));
    Mat<2, 2> pfast = add(identity<2>(), scaled(-1.0, pslow));
    return add(scaled(std::exp(slow * tau), pslow), scaled(std::exp(fast * tau), pfast));
}

inline Vec<2> exact(double contact, double tau, const Vec<2>& initial) {
    const auto eq = equilibrium(contact);
    return eq + mv(propagator(contact, tau), initial - eq);
}

inline double slow_rate(double contact) {
    const auto a = dynamics_matrix(contact);
    const double trace = a[0][0] + a[1][1];
    const double determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    return 0.5 * (trace + std::sqrt(trace * trace - 4.0 * determinant));
}

inline double fast_rate(double contact) {
    const auto a = dynamics_matrix(contact);
    const double trace = a[0][0] + a[1][1];
    const double determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    return 0.5 * (trace - std::sqrt(trace * trace - 4.0 * determinant));
}

inline std::vector<double> measurements(double contact) {
    std::vector<double> values;
    values.reserve(elapsed.size());
    for (std::size_t j = 0; j < elapsed.size(); ++j)
        values.push_back(exact(contact, elapsed[j], truth)[1] + sensor_sd * noise_in_sd[j]);
    return values;
}

inline Grid mesh(std::size_t subdivisions, double theta) {
    if (subdivisions == 0) throw std::invalid_argument("Positive thermal mesh subdivisions required");
    std::vector<double> time;
    for (std::size_t j = 0; j + 1 < elapsed.size(); ++j) {
        for (std::size_t k = 0; k < subdivisions; ++k)
            time.push_back(t0 + elapsed[j] + (elapsed[j + 1] - elapsed[j]) *
                           static_cast<double>(k) / static_cast<double>(subdivisions));
    }
    time.push_back(tf);
    return Grid(time, std::vector<double>(time.size() - 1, theta));
}

inline std::size_t node_at(const Grid& grid, double tau) {
    for (std::size_t n = 0; n < grid.time.size(); ++n)
        if (std::abs(grid.time[n] - t0 - tau) < 1e-10) return n;
    throw std::invalid_argument("Thermal observation is not a mesh node");
}

inline ScalarFunction<2, 2> housing_observation(double value) {
    auto cost = ScalarFunction<2, 2>::zero();
    cost.evaluate = [=](double, const Vec<2>& x, const Vec<2>&) {
        ScalarValue<2, 2> result;
        const double residual = x[1] - value;
        result.value = 0.5 * residual * residual / (sensor_sd * sensor_sd);
        result.x[1] = residual / (sensor_sd * sensor_sd);
        return result;
    };
    cost.directional_gradient = [](double, const Vec<2>&, const Vec<2>&,
                                   const Vec<2>& direction, const Vec<2>&) {
        ScalarDirection<2, 2> result;
        result.x[1] = direction[1] / (sensor_sd * sensor_sd);
        return result;
    };
    return cost;
}

inline Problem<2, 2> problem(const Grid& grid, double contact) {
    Problem<2, 2> p;
    const auto a = dynamics_matrix(contact);
    const auto eq = equilibrium(contact);
    p.initial.value = [](const Vec<2>& xi) { return xi; };
    p.initial.jacobian = [](const Vec<2>&) { return identity<2>(); };
    p.initial.directional_jacobian = [](const Vec<2>&, const Vec<2>&) { return Mat<2, 2>{}; };
    p.dynamics.value = [=](double, const Vec<2>& x, const Vec<2>&) { return mv(a, x - eq); };
    p.dynamics.derivatives = [=](double, const Vec<2>&, const Vec<2>&) {
        Derivatives<2, 2> d;
        d.x = a;
        return d;
    };
    p.dynamics.directional_derivatives = [](double, const Vec<2>&, const Vec<2>&,
                                            const Vec<2>&, const Vec<2>&) { return Derivatives<2, 2>{}; };
    const auto y = measurements(contact);
    for (std::size_t j = 0; j < elapsed.size(); ++j)
        p.nodal.emplace_back(node_at(grid, elapsed[j]), housing_observation(y[j]));
    return p;
}

inline System<2, 2> system(const Grid& grid, double contact) {
    return System<2, 2>(problem(grid, contact), grid, solver_options);
}

inline Optimum<2> fit(const System<2, 2>& sys, const Vec<2>& seed = {115.0, 55.0}) {
    return minimize(sys, seed, OptimizeOptions{1e-9, 8, 30});
}

inline Mat<2, 2> inverse(const Mat<2, 2>& a) {
    Mat<2, 2> result{};
    LU<2> lu(a);
    for (std::size_t j = 0; j < 2; ++j) {
        Vec<2> unit{};
        unit[j] = 1.0;
        const auto column = lu.solve(unit);
        for (std::size_t i = 0; i < 2; ++i) result[i][j] = column[i];
    }
    return result;
}

inline Mat<2, 2> outer(const Vec<2>& x) {
    return {{{x[0] * x[0], x[0] * x[1]}, {x[1] * x[0], x[1] * x[1]}}};
}

struct AnalyticalFit { Vec<2> initial{}; Mat<2, 2> information{}, covariance{}; };
inline AnalyticalFit analytical_fit(double contact) {
    AnalyticalFit result;
    const auto eq = equilibrium(contact);
    const auto y = measurements(contact);
    for (std::size_t j = 0; j < elapsed.size(); ++j) {
        const auto e = propagator(contact, elapsed[j]);
        const Vec<2> sensitivity{e[1][0], e[1][1]};
        result.information = add(result.information, scaled(1.0 / (sensor_sd * sensor_sd), outer(sensitivity)));
        result.initial = result.initial + ((y[j] - eq[1] + dot(sensitivity, eq)) /
                                           (sensor_sd * sensor_sd)) * sensitivity;
    }
    result.initial = mv(inverse(result.information), result.initial);
    result.covariance = inverse(result.information);
    return result;
}

inline Vec<2> physical_direction(const Vec<2>& v) {
    return {std::sqrt(winding_capacity) * v[0], std::sqrt(housing_capacity) * v[1]};
}
inline double energy_norm(const Vec<2>& v) {
    const auto w = physical_direction(v);
    return std::hypot(w[0], w[1]);
}
inline double goal_error_factor(const Vec<2>& goal) {
    return std::sqrt(goal[0] * goal[0] / winding_capacity + goal[1] * goal[1] / housing_capacity);
}
}
