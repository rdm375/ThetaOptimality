# Transient flow state and parameter analysis

The inverse problem is defined by a finite theta-method trajectory. This
document adds diagnostics derived from that exact discrete model: state
perturbation propagation, finite-horizon amplification, observation
sensitivities, and local parameter uncertainty. These diagnostics complement
the existing adjoint gradient, identifiability, and qualification tools.

## Discrete state transition maps

Write one implicit theta step as
\[
 R_k(x_{k+1},x_k,p)=x_{k+1}-x_k-h_k\{(1-\theta_k)f(t_k,x_k,p)
 +\theta_kf(t_{k+1},x_{k+1},p)\}=0.
\]
Linearization at the computed trajectory gives
\[
 B_k\,\delta x_{k+1}=C_k\,\delta x_k,\qquad
 M_k=B_k^{-1}C_k,
\]
where \(B_k=R_{x_{k+1}}\), \(C_k=-R_{x_k}\). The new
`System::state_transition_operators(trajectory,p)` computes each \(M_k\)
using the same step Jacobians and linear solver as the TLM. It is local to that
trajectory and fixed parameter vector. It is not a nonlinear stability proof.

`propagate_state_tangent` applies \(\delta x_{k+1}=M_k\delta x_k\), and
`propagate_state_adjoint` applies the reverse transpose. Their inner-product
duality is validated. Over several steps the ordered propagator is
\(\Phi_{n,0}=M_{n-1}\cdots M_0\), built by `state_transition_product`.
Multiplication order matters because the maps generally do not commute.

## Finite-horizon amplification

The largest transient gain is the largest singular value of the propagator
after scaling input and output states into comparable dimensionless units:
\[
 G=D_{\rm out}^{-1}\Phi_{n,0}D_{\rm in},\qquad
 g_{\max}=\sigma_{\max}(G).
\]
`finite_horizon_amplification` returns that gain, all singular values, the
dominant scaled input direction, and its physical state response. Positive
scales are mandatory: mixing pressure, flow, temperature, or composition
without scaling makes the Euclidean norm physically arbitrary.

Large gain means some small initial-state perturbation can become large over
the selected horizon. It does not imply an unstable equilibrium, nor does a
small gain imply parameter identifiability. For non-normal or time-varying
systems, eigenvalue magnitudes or each step's spectral radius alone can miss
finite-time growth. Inspect gains for several horizons and use the scales,
mesh, trajectory, and norm in the report. These dense fixed-size matrices are
intended for small/medium teaching and diagnostics problems.

The theta value is part of the discrete model. A value such as 0.5 gives the
trapezoidal rule and second-order accuracy for sufficiently smooth solutions
with consistent nonlinear solves, but it is not universally energy
conserving or unconditionally stable for every nonlinear pipe formulation.
Check temporal refinement and physical invariants for the actual equations.

## Observation sensitivity and observability

For a scalar observation \(y_i=h_i(x_{k_i},p)\), its row in the local
parameter-sensitivity matrix is
\[
 S_i= h_{i,x}\,\frac{\partial x_{k_i}}{\partial p}+h_{i,p}.
\]
`observation_sensitivities` computes the state-mediated part by the discrete
adjoint `System::goal_gradient` and adds the direct-control derivative. Each
observation supplies its node, state gradient, direct control gradient, and
standard deviation. The resulting rows are derivatives of the implemented
discrete observations, not of a separate continuous approximation.

Scale parameter columns and observation rows before interpreting singular
values. A small singular value indicates a weakly observed parameter
combination relative to the chosen scales and local linearization. Its right
singular vector names that combination; it does not say each named parameter
is individually unidentifiable. Numerical rank depends on the cutoff, scales,
sampling schedule, and trajectory. Recompute it after changing any of these.

## Fisher information, covariance, and correlation

For independent, zero-mean observation errors with known standard deviations
\(\sigma_i\), the local Fisher/Gauss–Newton information is
\[
 F=S^T\Sigma^{-1}S=\sum_i S_i^TS_i/\sigma_i^2.
\]
`fisher_information` assembles this matrix. If it is full rank and the local
linear-Gaussian approximation is defensible, `covariance_from_precision(F)`
returns \(F^{-1}\), and helpers return marginal standard errors and the
parameter correlation matrix. A rank-deficient matrix raises an error rather
than reporting finite covariance in an unobserved direction. For correlated
measurement errors, whiten the sensitivity rows with the full covariance (or
extend the assembly accordingly); passing only marginal standard deviations
assumes independence.

For a quadratic Gaussian prior with precision \(K_0\),
`add_information_precision(F,K0)` forms the local posterior precision
\(F_{\rm post}=F+K_0\). Its inverse is a prior-conditioned covariance, not a
data-only observability result; report both separately.

The condition number describes anisotropy after scaling, not absolute
accuracy. High off-diagonal correlation indicates that changes in two
parameters compensate in the local fit. Large standard errors indicate weak
information under the assumed noise model. These are local approximations:
they omit model discrepancy, nonlinear posterior shape, bounds, multi-modal
solutions, and numerical discretization error unless those are explicitly
included. A regularizer can make a posterior precision invertible; report
data-only information separately so prior information is not mistaken for
observability.

## How to use the diagnostics

1. Verify the forward solve and refine time and space until outputs and
   sensitivities are numerically stable at the intended scale.
2. Define physical state and parameter scales; construct state maps and inspect
   finite-horizon gains for operating schedules of interest.
3. Define observation operators, sampling times, and a defensible error
   covariance; compute sensitivity rows and scaled singular directions.
4. Assemble Fisher information and report rank, cutoff, condition number,
   covariance, standard errors, and correlations with the assumptions.
5. If weak directions remain, add informative sensors/excitation or an
   explicit prior/smoothness assumption. Then distinguish data information
   from regularization and validate on changed schedules or independent data.

The numerical validation suite checks the transition-product/TLM agreement,
tangent-adjoint duality, known amplification, observation pullback, Fisher
covariance, rank-deficient rejection, regularizer derivatives, and weak
constraint gradients. See `tests/validation.cpp` and the regularizer overview
in [matrix-based regularizers](matrix_based_regularizers.md).
