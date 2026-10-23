# Matrix-based regularizers

This note describes additive regularization utilities that operate on the
finite-dimensional control vector or on a candidate state trajectory. They
use the same controls and discrete trajectory as the existing theta-method
forward solve, tangent linear model (TLM), and discrete adjoint. Regularization
adds prior or smoothness information; it does not create information in the
observations.

## Quadratic prior and its derivatives

For a control vector \(p\), reference value \(p_0\), and symmetric positive
semidefinite precision matrix \(K\), the penalty is

\[
  J_{\rm prior}(p)=\tfrac12(p-p_0)^T K(p-p_0),\qquad
  \nabla J_{\rm prior}=K(p-p_0),\qquad
  H_{\rm prior}v=Kv.
\]

If the prior covariance \(C_0\) is positive definite, take \(K=C_0^{-1}\).
`prior_precision_regularizer` accepts \(K\) directly;
`independent_prior_regularizer` constructs a diagonal precision from standard
deviations. `as_penalty` converts the result to the existing `Penalty<P>`
callback type, and `add_penalty` combines it with the problem's current
penalty. The exact Hessian-vector callback is retained only when both summed
penalties provide one.

Covariance and precision must be expressed in the same parameter coordinates
as the optimizer. For example, a covariance in physical diameter units cannot
be applied directly to log-diameter controls. Transform the covariance (or
linearize the coordinate transform locally) before constructing the prior.
A positive semidefinite precision may intentionally leave directions
unpenalized; it must not be described as a proper probability density unless
it is positive definite on the modeled parameter space.

## Sensitivity-nullspace penalty

Let \(S\) be a scaled observation-sensitivity matrix and let columns of
\(V_0\) be an orthonormal basis for its numerically weak right-singular
subspace. `subspace_regularizer` forms

\[
 K_{\rm null}=\lambda V_0V_0^T,
 \qquad J_{\rm null}=\tfrac12\lambda
 \|V_0^T(p-p_0)\|_2^2.
\]

This selects a solution near a stated center in directions the current data
resolve poorly. It does not make those directions observed, and the returned
uncertainty should still report the data-only nullspace separately from the
regularized posterior. The basis supplied to the helper must be orthonormal;
derive it from a sensitivity SVD after choosing meaningful parameter and
observation scales.

## Spatial H1 and H2 penalties

For equally spaced nodal values with interval \(h\),
`spatial_difference_regularizer< P >(1,...)` discretizes the H1 seminorm

\[
  J_{H1}=\tfrac12\lambda\sum_i (p_{i+1}-p_i)^2/h,
\]

and the order-2 option discretizes the H2 seminorm

\[
  J_{H2}=\tfrac12\lambda\sum_i
  (p_{i+2}-2p_{i+1}+p_i)^2/h^3.
\]

H1 does not penalize a constant offset. H2 does not penalize an affine trend.
These nullspaces are deliberate; add an appropriate prior or boundary
information if those modes must be controlled. These helpers assume uniform
spacing and nodal values, not arbitrary finite-element bases or nonuniform
grids. The reference vector is the center of the penalty, so a nonzero
reference means smoothness is imposed on deviations from that profile.

## Weak-constraint trajectory penalty

The strong-constraint `System` continues to enforce every discrete model
equation exactly. For state estimation with model discrepancy, the standalone
`weak_constraint_trajectory_penalty` instead evaluates a candidate trajectory
against a reference trajectory and local TLM maps:

\[
 e_k=x_k-\bar x_k,\quad r_k=e_{k+1}-M_ke_k,\quad
 J_{\rm wc}=\tfrac12\sum_k r_k^TQ_k^{-1}r_k.
\]

Its state gradient is assembled as
\[
 \nabla_{x_k}J_{\rm wc}\;{+}{=}-M_k^TQ_k^{-1}r_k,\qquad
 \nabla_{x_{k+1}}J_{\rm wc}\;{+}{=}Q_k^{-1}r_k.
\]
The helper accepts process precision matrices (the API name is
`process_precision`, corresponding to \(Q_k^{-1}\)). They must be symmetric
positive semidefinite. The reference trajectory is essential: without it,
nonzero nominal dynamics would be incorrectly penalized as model error.
This utility returns cost, defects, and per-node state gradients; callers add
these to their observation and prior terms in a weak-constraint estimator.
It does not change the existing strong-constraint integrator or claim a
complete weak-constraint optimization interface.

## Interpretation and checks

Regularization strength has units implied by the parameter coordinates and
the cost scaling. Choose it from a stated prior uncertainty, cross-validation,
or a documented discrepancy model rather than tuning it only to make a fit
look smooth. Report data-only sensitivity rank and uncertainty beside the
regularized solution. Tests check the analytic quadratic derivatives, the
nullspaces of H1/H2, and the weak-constraint state gradient against centered
finite differences. See `tests/validation.cpp`.
