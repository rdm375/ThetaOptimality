// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include "identifiability.hpp"
#include "system.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace theta_opt {

// Convex quadratic penalty 1/2 (p-center)^T precision (p-center).
template<std::size_t P>
struct QuadraticRegularizer {
    Mat<P,P> precision{};
    Vec<P> center{};

    QuadraticRegularizer(const Mat<P,P>& matrix,const Vec<P>& reference,
                         double relative_cutoff=1e-12)
        : precision(matrix),center(reference) {
        (void)norm(center);
        (void)psd_spectrum(precision,relative_cutoff);
    }
    double value(const Vec<P>& p) const {
        const auto d=p-center;
        const double result=.5*dot(d,mv(precision,d));
        if(!std::isfinite(result))throw std::overflow_error("Regularizer value overflow");
        return result;
    }
    Vec<P> gradient(const Vec<P>& p) const { return mv(precision,p-center); }
    Vec<P> hessian_vector(const Vec<P>&,const Vec<P>& v) const {
        return mv(precision,v);
    }
};

template<std::size_t P>
QuadraticRegularizer<P> prior_precision_regularizer(
    const Mat<P,P>& prior_precision,const Vec<P>& prior_mean,
    double relative_cutoff=1e-12) {
    return QuadraticRegularizer<P>(prior_precision,prior_mean,relative_cutoff);
}

template<std::size_t P>
QuadraticRegularizer<P> independent_prior_regularizer(
    const Vec<P>& prior_mean,const Vec<P>& standard_deviations) {
    Mat<P,P> precision{};
    for(std::size_t i=0;i<P;++i){
        if(!std::isfinite(standard_deviations[i])||standard_deviations[i]<=0.)
            throw std::invalid_argument("Positive finite prior standard deviations required");
        precision[i][i]=1./(standard_deviations[i]*standard_deviations[i]);
        if(!std::isfinite(precision[i][i]))throw std::overflow_error("Prior precision overflow");
    }
    return QuadraticRegularizer<P>(precision,prior_mean);
}

template<std::size_t P>
Penalty<P> as_penalty(const QuadraticRegularizer<P>& regularizer) {
    Penalty<P> result;
    result.value=[regularizer](const Vec<P>& p){return regularizer.value(p);};
    result.gradient=[regularizer](const Vec<P>& p){return regularizer.gradient(p);};
    result.hessian_vector=[regularizer](const Vec<P>& p,const Vec<P>& v){
        return regularizer.hessian_vector(p,v);
    };
    return result;
}

template<std::size_t P>
Penalty<P> add_penalty(Penalty<P> existing,const Penalty<P>& additional) {
    if(!existing.value||!existing.gradient||!additional.value||!additional.gradient)
        throw std::invalid_argument("Penalty callbacks are required");
    const auto old_value=std::move(existing.value);
    const auto old_gradient=std::move(existing.gradient);
    const auto old_hessian=std::move(existing.hessian_vector);
    existing.value=[old_value,additional](const Vec<P>& p){return old_value(p)+additional.value(p);};
    existing.gradient=[old_gradient,additional](const Vec<P>& p){return old_gradient(p)+additional.gradient(p);};
    if(old_hessian&&additional.hessian_vector)
        existing.hessian_vector=[old_hessian,additional](const Vec<P>& p,const Vec<P>& v){
            return old_hessian(p,v)+additional.hessian_vector(p,v);
        };
    else existing.hessian_vector={};
    return existing;
}

// A projector penalty for an orthonormal basis of weakly observed directions.
template<std::size_t P,std::size_t K>
QuadraticRegularizer<P> subspace_regularizer(
    const std::array<Vec<P>,K>& basis,double weight,const Vec<P>& center,
    double orthogonality_tolerance=1e-10) {
    if(!std::isfinite(weight)||weight<0.||!std::isfinite(orthogonality_tolerance)||
       orthogonality_tolerance<=0.)throw std::invalid_argument("Invalid subspace penalty settings");
    Mat<P,P> precision{};
    for(std::size_t k=0;k<K;++k){
        (void)norm(basis[k]);
        for(std::size_t j=0;j<K;++j){
            const double expected=k==j?1.:0.;
            if(std::abs(dot(basis[k],basis[j])-expected)>orthogonality_tolerance)
                throw std::invalid_argument("Subspace basis must be orthonormal");
        }
        for(std::size_t i=0;i<P;++i)
            for(std::size_t j=0;j<P;++j)
                precision[i][j]+=weight*basis[k][i]*basis[k][j];
    }
    return QuadraticRegularizer<P>(precision,center);
}

// Spacing-aware H1/H2 smoothness precision for nodal values. H1 uses
// integral |u'|^2; H2 uses integral |u''|^2 and leaves affine profiles unpenalized.
template<std::size_t P>
QuadraticRegularizer<P> spatial_difference_regularizer(
    unsigned order,double spacing,double weight,const Vec<P>& reference) {
    if((order!=1&&order!=2)||!std::isfinite(spacing)||spacing<=0.||
       !std::isfinite(weight)||weight<0.)
        throw std::invalid_argument("H1/H2 require order 1 or 2, positive spacing, nonnegative weight");
    Mat<P,P> precision{};
    if(order==1){
        for(std::size_t k=0;k+1<P;++k){
            const double c=weight/spacing;
            precision[k][k]+=c; precision[k+1][k+1]+=c;
            precision[k][k+1]-=c; precision[k+1][k]-=c;
        }
    }else if constexpr(P>=3){
        const double c=weight/(spacing*spacing*spacing);
        for(std::size_t k=0;k+2<P;++k){
            const std::array<double,3> d{1.,-2.,1.};
            for(std::size_t i=0;i<3;++i)
                for(std::size_t j=0;j<3;++j)
                    precision[k+i][k+j]+=c*d[i]*d[j];
        }
    }else if(order==2){
        throw std::invalid_argument("H2 regularization requires at least three nodes");
    }
    return QuadraticRegularizer<P>(precision,reference);
}

template<std::size_t N,std::size_t Steps>
struct WeakConstraintResult {
    double value{};
    std::array<Vec<N>,Steps+1> gradient{};
    std::array<Vec<N>,Steps> defect{};
};

// Penalize defects in deviations from a reference trajectory:
// r_k=(x_{k+1}-xbar_{k+1})-M_k(x_k-xbar_k), J=.5 sum r_k^T Q_k^-1 r_k.
template<std::size_t N,std::size_t Steps>
WeakConstraintResult<N,Steps> weak_constraint_trajectory_penalty(
    const std::array<Vec<N>,Steps+1>& state,
    const std::array<Vec<N>,Steps+1>& reference,
    const std::array<Mat<N,N>,Steps>& transition,
    const std::array<Mat<N,N>,Steps>& process_precision) {
    WeakConstraintResult<N,Steps> result;
    for(const auto& x:state)(void)norm(x);
    for(const auto& x:reference)(void)norm(x);
    for(std::size_t k=0;k<Steps;++k){
        const auto dx=state[k]-reference[k];
        const auto next=state[k+1]-reference[k+1];
        result.defect[k]=next-mv(transition[k],dx);
        const auto weighted=mv(process_precision[k],result.defect[k]);
        (void)psd_spectrum(process_precision[k]);
        result.value+=.5*dot(result.defect[k],weighted);
        result.gradient[k]=result.gradient[k]-mv(transpose(transition[k]),weighted);
        result.gradient[k+1]=result.gradient[k+1]+weighted;
    }
    if(!std::isfinite(result.value))throw std::overflow_error("Weak-constraint cost overflow");
    for(const auto& g:result.gradient)(void)norm(g);
    return result;
}

} // namespace theta_opt
