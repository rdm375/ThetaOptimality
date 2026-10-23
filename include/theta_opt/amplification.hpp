// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Richard Myers

#pragma once
#include "system.hpp"

namespace theta_opt {
// Propagate a state perturbation through discrete theta-step TLM maps.
template<std::size_t N>
std::vector<Vec<N>> propagate_state_tangent(
    const Vec<N>& initial,const std::vector<Mat<N,N>>& transition) {
    (void)norm(initial);
    std::vector<Vec<N>> result(transition.size()+1);
    result[0]=initial;
    for(std::size_t k=0;k<transition.size();++k)
        result[k+1]=mv(transition[k],result[k]);
    return result;
}

// Pull a terminal state gradient backward through the same maps. This is the
// transpose of propagate_state_tangent under the Euclidean inner product.
template<std::size_t N>
std::vector<Vec<N>> propagate_state_adjoint(
    const Vec<N>& terminal,const std::vector<Mat<N,N>>& transition) {
    (void)norm(terminal);
    std::vector<Vec<N>> result(transition.size()+1);
    result.back()=terminal;
    for(std::size_t k=transition.size();k-->0;)
        result[k]=mv(transpose(transition[k]),result[k+1]);
    return result;
}

// Jacobians of the implemented state trajectory with respect to its controls.
// Includes the initial-condition map and direct dynamics/control dependence.
template<std::size_t N,std::size_t P>
std::vector<Mat<N,P>> control_to_state_jacobians(const System<N,P>& system,const Vec<P>& xi,const Trajectory<N>& trajectory){
    std::vector<Mat<N,P>> result(trajectory.state.size());
    for(std::size_t j=0;j<P;++j){Vec<P> direction{};direction[j]=1.;auto tangent=system.tangent(xi,trajectory,direction);for(std::size_t n=0;n<result.size();++n)for(std::size_t i=0;i<N;++i)result[n][i][j]=tangent[n][i];}return result;
}
struct Amplification2D {double gain{};Vec<2> scaled_input{},input_direction{},state_response{};};
// Largest singular gain of D_out^-1 J D_in; scales carry physical units.
// A deliberately small, closed-form 2D utility, not a general SVD backend.
inline Amplification2D worst_amplification_2d(const Mat<2,2>& J,const Vec<2>& input_scales={1.,1.},const Vec<2>& output_scales={1.,1.}){
    Mat<2,2> K{};double maximum=0.;
    for(int i=0;i<2;++i){if(!std::isfinite(input_scales[i])||!std::isfinite(output_scales[i])||input_scales[i]<=0||output_scales[i]<=0)throw std::invalid_argument("Positive finite amplification scales required");}
    for(int i=0;i<2;++i)for(int j=0;j<2;++j){K[i][j]=J[i][j]*input_scales[j]/output_scales[i];if(!std::isfinite(K[i][j]))throw std::invalid_argument("Nonfinite scaled amplification matrix");maximum=std::max(maximum,std::abs(K[i][j]));}
    Amplification2D result;result.scaled_input={1.,0.};
    if(maximum>0){K=scaled(1/maximum,K);double a=K[0][0]*K[0][0]+K[1][0]*K[1][0],b=K[0][0]*K[0][1]+K[1][0]*K[1][1],d=K[0][1]*K[0][1]+K[1][1]*K[1][1];double lambda=.5*(a+d+std::hypot(a-d,2*b));Vec<2> v{b,lambda-a},other{lambda-d,b};if(std::hypot(other[0],other[1])>std::hypot(v[0],v[1]))v=other;double length=std::hypot(v[0],v[1]);if(length>0)result.scaled_input=(1/length)*v;result.gain=maximum*std::sqrt(lambda);}
    if(result.scaled_input[0]<0||(result.scaled_input[0]==0&&result.scaled_input[1]<0))result.scaled_input=-1.*result.scaled_input;
    for(int j=0;j<2;++j){result.input_direction[j]=input_scales[j]*result.scaled_input[j];}
    result.state_response=mv(J,result.input_direction);
    if(!std::isfinite(result.gain))throw std::overflow_error("Amplification gain overflow");
    return result;
}
}
