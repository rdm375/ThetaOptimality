# Time-varying controls and discrete path-constraint KKT conditions

This library extension lets a finite set of optimization coefficients represent
a time-varying control and lets state/control inequalities participate in the
same theta-discrete optimization system. It is a small, dense implementation
for the repository's verification examples.

## Discrete problem

The state still satisfies the chosen theta scheme on a fixed grid. The reduced
problem has coefficient vector $\boldsymbol z\in\mathbb R^P$, control history
$\boldsymbol u(t;\boldsymbol z)$, and selected nodal inequalities

$$
g_i(t_{n_i},\boldsymbol x_{n_i}(\boldsymbol z),\boldsymbol z)\le 0,
\qquad i=1,\ldots,m.
$$

For box bounds $\boldsymbol z_L\le\boldsymbol z\le\boldsymbol z_U$, the
discrete KKT conditions are

$$
\nabla_{\boldsymbol z}J(\boldsymbol z)
+\sum_i\mu_i\nabla_{\boldsymbol z}g_i(\boldsymbol z)
-\boldsymbol\lambda_L+\boldsymbol\lambda_U=0,
$$

$$
g_i\le0,\quad \mu_i\ge0,\quad \mu_i g_i=0,
\qquad \boldsymbol\lambda_L,\boldsymbol\lambda_U\ge0,
$$

with the corresponding box feasibility and complementarity conditions. The
gradients of state-dependent $g_i$ are reduced through the discrete adjoint.
Their second derivatives enter the augmented objective through the
incremental-adjoint callback.

## Library interface

Add each inequality to `Problem<N,P>::path_inequalities` as
`PathInequality<N,P>{node, function}`. Its scalar function follows the same
callback contract as a nodal objective term and must provide
`directional_gradient` when second-order optimization is used. A path callback
may depend on the state and on the coefficient vector. To express a constraint
on an instantaneous scheduled control, first use `compose_scalar_function`.

Call:

```cpp
auto result = minimize_path_constrained(
    system, initial_coefficients, coefficient_lower, coefficient_upper,
    optimize_options, path_options);
```

The solver uses an augmented-Lagrangian iteration. It returns the fitted
coefficients, path values and multipliers, box multipliers, and residuals for
stationarity, primal feasibility, and complementarity. The test suite includes
an analytically solvable dynamic example: minimize
$\frac12(x(t_f)-1)^2$ subject to $x'=z$, $x(t_0)=0$, and
$x(t_0+0.5)\le0.4$. The exact optimum is $z=0.8$ and the path multiplier is
$0.4$. The check validates the state-path gradient through the adjoint and the
path Hessian through the incremental adjoint.

The constraints are imposed at selected theta-grid nodes. They do not certify
the inequality between nodes. A continuous-time claim needs mesh refinement,
additional collocation points, or a separately justified interval bound.

## Time-varying controls

`PiecewiseLinearControl<C,K>` uses the values of $C$ channels at $K$ fixed
knots as optimization coefficients. `SoftmaxCompositionControl<S,K>` uses
$S-1$ log-ratios at each knot and returns positive fractions summing to one at
every time. It therefore represents a changing gas composition without adding
a separate equality constraint for the mole-fraction sum.

`compose_dynamics` chains the model's derivatives with respect to the
instantaneous control through the schedule. `compose_scalar_function` performs
the same chain rule for running/terminal terms, observations, and path
inequalities. Both propagate directional second derivatives, preserving the
library's tangent, discrete adjoint, and incremental-adjoint calculations for
the chosen finite-dimensional control basis.

For a dynamic gas-mixture example, the schedule can supply mole fractions to a
mixing/residence model and to a thermodynamic property observation such as
density or compressibility. The C++ callback must provide derivatives of that
property with respect to composition (and state, when applicable). The current
extension supplies the schedule and theta-optimality chain rule; it does not
yet implement a transient gas-mixture experiment or automatically differentiate
the ThermoGPU equation-of-state routines.

## Scope and verification

The tests check the softmax simplex, schedule Jacobian and directional
Jacobian, composed dynamics tangent, discrete-adjoint gradient, and
incremental-adjoint Hessian by finite differences. They also check the dynamic
path-constrained optimum and KKT residuals against its analytic solution.

The current solver is dense and intended for small problems. It supports
nodal inequalities and simple coefficient box bounds; general equality
constraints, sparse large-scale SQP/interior-point methods, adaptive control
knots, and between-node path certification remain future work. The finite
control basis is part of the model: a recovered knot schedule is not evidence
that an arbitrary continuous control history is uniquely determined.
