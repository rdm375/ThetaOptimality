// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include "amplification.hpp"
#include "identifiability.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace theta_opt {

template<std::size_t N>
Mat<N,N> state_transition_product(const std::vector<Mat<N,N>>& transition) {
    Mat<N,N> product=identity<N>();
    for(const auto& step:transition){
        for(const auto& row:step)(void)norm(row);
        Mat<N,N> next{};
        for(std::size_t i=0;i<N;++i)
            for(std::size_t j=0;j<N;++j)
                for(std::size_t k=0;k<N;++k)
                    next[i][j]+=step[i][k]*product[k][j];
        product=next;
    }
    return product;
}

template<std::size_t N>
struct StateAmplification {
    double gain{};
    Vec<N> input_direction{};
    Vec<N> state_response{};
    std::array<double,N> singular_values{};
};

// Largest finite-horizon gain of D_out^-1 Phi D_in, computed from the
// eigenvalues of its symmetric Gram matrix. Scales define physical units.
template<std::size_t N>
StateAmplification<N> worst_state_amplification(
    const Mat<N,N>& propagator,const Vec<N>& input_scales,
    const Vec<N>& output_scales,double relative_cutoff=1e-12) {
    Mat<N,N> scaled_map{};
    for(std::size_t i=0;i<N;++i){
        if(!std::isfinite(input_scales[i])||input_scales[i]<=0.||
           !std::isfinite(output_scales[i])||output_scales[i]<=0.)
            throw std::invalid_argument("Positive finite state scales required");
        for(std::size_t j=0;j<N;++j){
            scaled_map[i][j]=propagator[i][j]*input_scales[j]/output_scales[i];
            if(!std::isfinite(scaled_map[i][j]))
                throw std::invalid_argument("Nonfinite scaled state propagator");
        }
    }
    Mat<N,N> gram{};
    for(std::size_t i=0;i<N;++i)
        for(std::size_t j=0;j<N;++j)
            for(std::size_t k=0;k<N;++k)
                gram[i][j]+=scaled_map[k][i]*scaled_map[k][j];
    const auto spectrum=psd_spectrum(gram,relative_cutoff);
    StateAmplification<N> result;
    for(std::size_t i=0;i<N;++i)
        result.singular_values[i]=std::sqrt(std::max(0.,spectrum.eigenvalues[i]));
    result.input_direction=spectrum.directions[0];
    if(spectrum.eigenvalues[0]>0.)
        result.gain=result.singular_values[0];
    else {
        result.input_direction={};
        result.input_direction[0]=1.;
    }
    // Choose a deterministic sign for reproducible diagnostics.
    for(double component:result.input_direction){
        if(component!=0.){
            if(component<0.)result.input_direction=-1.*result.input_direction;
            break;
        }
    }
    Vec<N> physical_input{};
    for(std::size_t i=0;i<N;++i)
        physical_input[i]=input_scales[i]*result.input_direction[i];
    result.state_response=mv(propagator,physical_input);
    (void)norm(result.state_response);
    return result;
}

template<std::size_t N>
StateAmplification<N> finite_horizon_amplification(
    const std::vector<Mat<N,N>>& transition,const Vec<N>& input_scales,
    const Vec<N>& output_scales,double relative_cutoff=1e-12) {
    return worst_state_amplification(state_transition_product(transition),
                                     input_scales,output_scales,relative_cutoff);
}

template<std::size_t N,std::size_t P>
struct LinearizedObservation {
    std::size_t node{};
    Vec<N> state_gradient{};
    Vec<P> direct_control_gradient{};
    double standard_deviation{1.};
};

// Each output row is d y_i / d xi, assembled with the system's discrete
// adjoint and including any direct observation dependence on the controls.
template<std::size_t N,std::size_t P>
std::vector<Vec<P>> observation_sensitivities(
    const System<N,P>& system,const Vec<P>& xi,const Trajectory<N>& trajectory,
    const std::vector<LinearizedObservation<N,P>>& observations) {
    std::vector<Vec<P>> rows;
    rows.reserve(observations.size());
    for(const auto& observation:observations){
        if(observation.node>=trajectory.state.size())
            throw std::invalid_argument("Observation node out of range");
        (void)norm(observation.state_gradient);
        (void)norm(observation.direct_control_gradient);
        if(!std::isfinite(observation.standard_deviation)||
           observation.standard_deviation<=0.)
            throw std::invalid_argument("Positive observation standard deviation required");
        rows.push_back(system.goal_gradient(xi,trajectory,observation.node,
                                            observation.state_gradient)+
                       observation.direct_control_gradient);
    }
    return rows;
}

template<std::size_t P>
Mat<P,P> fisher_information(const std::vector<Vec<P>>& sensitivity_rows,
                            const std::vector<double>& standard_deviations) {
    if(sensitivity_rows.size()!=standard_deviations.size())
        throw std::invalid_argument("Sensitivity/noise dimensions differ");
    Mat<P,P> information{};
    for(std::size_t k=0;k<sensitivity_rows.size();++k){
        (void)norm(sensitivity_rows[k]);
        const double sigma=standard_deviations[k];
        if(!std::isfinite(sigma)||sigma<=0.)
            throw std::invalid_argument("Positive observation standard deviation required");
        const double weight=1./(sigma*sigma);
        if(!std::isfinite(weight))throw std::overflow_error("Observation precision overflow");
        for(std::size_t i=0;i<P;++i)
            for(std::size_t j=0;j<P;++j)
                information[i][j]+=weight*sensitivity_rows[k][i]*sensitivity_rows[k][j];
    }
    for(const auto& row:information)(void)norm(row);
    return information;
}

// Combine data information and a quadratic-prior precision in shared control
// coordinates. The result may still be rank deficient if the prior is weak.
template<std::size_t P>
Mat<P,P> add_information_precision(const Mat<P,P>& data_information,
                                   const Mat<P,P>& prior_precision) {
    Mat<P,P> result{};
    for(std::size_t i=0;i<P;++i)
        for(std::size_t j=0;j<P;++j)
            result[i][j]=data_information[i][j]+prior_precision[i][j];
    (void)psd_spectrum(data_information);
    (void)psd_spectrum(prior_precision);
    (void)psd_spectrum(result);
    return result;
}

// A covariance is returned only when the supplied precision is numerically
// full rank. Unidentified directions are not assigned finite uncertainty.
template<std::size_t P>
Mat<P,P> covariance_from_precision(const Mat<P,P>& precision,
                                   double relative_cutoff=1e-12) {
    const auto spectrum=psd_spectrum(precision,relative_cutoff);
    if(spectrum.rank!=P)
        throw std::domain_error("Precision is rank deficient; covariance is unbounded in null directions");
    Mat<P,P> covariance{};
    for(std::size_t k=0;k<P;++k)
        for(std::size_t i=0;i<P;++i)
            for(std::size_t j=0;j<P;++j)
                covariance[i][j]+=spectrum.directions[k][i]*spectrum.directions[k][j]/
                                  spectrum.eigenvalues[k];
    for(const auto& row:covariance)(void)norm(row);
    return covariance;
}

template<std::size_t P>
Vec<P> standard_errors(const Mat<P,P>& covariance) {
    Vec<P> result{};
    for(std::size_t i=0;i<P;++i){
        for(std::size_t j=0;j<P;++j)
            if(!std::isfinite(covariance[i][j])||
               std::abs(covariance[i][j]-covariance[j][i])>
                    1e-10*std::max({1.,std::abs(covariance[i][j]),std::abs(covariance[j][i])}))
                throw std::invalid_argument("Finite symmetric covariance required");
        if(!(covariance[i][i]>0.))
            throw std::invalid_argument("Positive covariance diagonal required");
        result[i]=std::sqrt(covariance[i][i]);
    }
    return result;
}

template<std::size_t P>
Mat<P,P> correlation_matrix(const Mat<P,P>& covariance) {
    const auto sd=standard_errors(covariance);
    Mat<P,P> result{};
    for(std::size_t i=0;i<P;++i)
        for(std::size_t j=0;j<P;++j){
            result[i][j]=covariance[i][j]/(sd[i]*sd[j]);
            if(!std::isfinite(result[i][j]))throw std::overflow_error("Correlation overflow");
        }
    return result;
}

template<std::size_t P>
double sensitivity_condition_number(const Mat<P,P>& information,
                                    double relative_cutoff=1e-12) {
    const auto spectrum=psd_spectrum(information,relative_cutoff);
    if(spectrum.rank<P)return std::numeric_limits<double>::infinity();
    return std::sqrt(spectrum.eigenvalues[0]/spectrum.eigenvalues[P-1]);
}

} // namespace theta_opt
