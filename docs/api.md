# API and callback contracts

Include `<theta_opt/theta_opt.hpp>` and link no compiled library. Public symbols live in `theta_opt`. `N` is state dimension and `P` is control dimension. Both must be positive. `Vec<N>` and `Mat<N,P>` are `std::array` aliases; import `theta_opt` operators explicitly or use `using namespace theta_opt` when writing vector arithmetic outside the namespace, because the underlying arrays belong to `std` for argument-dependent lookup.

## Problem definition

`Problem<N,P>` owns these callback groups by value:

- `dynamics.value(t,x,xi)` returns the ODE right-hand side.
- `dynamics.derivatives(t,x,xi)` returns `Derivatives<N,P>`. Its `.x` is the state Jacobian and `.control` is the direct control Jacobian, both evaluated with time fixed.
- `dynamics.directional_derivatives(t,x,xi,s,v)` returns the directional derivatives of those two Jacobians along state direction `s` and control direction `v`. It is needed only by second-order operations.
- `dynamics.admissible(t,x,xi)` returns whether a trial state is acceptable. The default accepts every state; applications may enforce positivity. It is a domain predicate, not an inequality-constrained optimality condition.
- `initial.value(xi)` and `initial.jacobian(xi)` return the initial state and its `N` by `P` derivative. `initial.directional_jacobian(xi,v)` differentiates that matrix along `v`.
- `running.evaluate(t,x,xi)` and `terminal.evaluate(t,x,xi)` return `ScalarValue<N,P>{value, x_gradient, control_gradient}`. Control gradients hold the state fixed. The terminal function is evaluated only at the final node.
- `running.directional_gradient(t,x,xi,s,v)` and the terminal equivalent return `ScalarDirection<N,P>{delta_x_gradient, delta_control_gradient}`.
- `penalty.value(xi)`, `penalty.gradient(xi)`, and `penalty.hessian_vector(xi,v)` depend directly on controls.
- `nodal` is a vector of `(node, ScalarFunction<N,P>)` pairs for isolated measurements or other discrete costs. They use the same scalar-function callbacks, and are added without time weights. Multiple costs may share a node, including the initial and final node. First- and second-order sources accumulate in the appropriate sweeps.
- `path_inequalities` is a vector of `PathInequality<N,P>{node,function}` entries, each imposing `function(t_node,x_node,xi) <= 0` at that mesh node. Add one entry for every constrained node and inequality component. These are discretized path constraints; node feasibility alone is not a certificate between nodes.

The defaults for `Problem::running`, `terminal`, and `penalty` are mathematically zero and include zero second derivatives. A newly constructed generic `ScalarFunction` or `Penalty` does not silently assume second derivatives when nonzero first-order callbacks have been assigned. Prefer defining the complete callback group before constructing a `System`. Missing first-order callbacks are rejected by its constructor; missing second-order callbacks are rejected by `hessian_vector`.

Callbacks must be pure with respect to a solve: the same arguments must give the same output. The code does not finite-difference or automatically differentiate user callbacks. Analytical callbacks, or future AD adapters, must return consistent derivatives themselves. No derivative of the time grid or theta values is taken.

## Grid and numerical options

`Grid(time,theta)` takes at least two finite strictly increasing times, and one theta value per interval in `[0,1]`. `Grid::uniform(t0,tf,steps,theta)` constructs equally spaced nodes. `Grid::weights()` returns the collected running-cost quadrature weights.

`NewtonOptions` contains absolute and relative residual tolerances, maximum Newton updates per step, and maximum backtracks per Newton update. The default convergence test uses the infinity norm of the residual and the state scale. A converged residual and an admissible state are both required. Scaling is the application's responsibility: variables with very different physical magnitudes should be nondimensionalized or otherwise scaled before use.

The built-in linear solver is dense LU with partial pivoting. It rejects nonfinite matrices and pivots small relative to the matrix scale. This is appropriate for the small teaching examples and is not a condition-number certificate. There is no sparse factorization, mass-matrix specialization, or matrix-free state solver.

## Numerical operations

`System<N,P>(problem,grid,options)` owns copies of the problem and grid. It exposes its grid through a const reference.

| Method | Result |
|---|---|
| `forward(xi)` | `Trajectory<N>` with states, largest converged step residual, and total Newton updates |
| `cost(xi,trajectory)` | Running theta quadrature plus terminal cost and penalty |
| `tangent(xi,trajectory,v)` | One state perturbation per node, including the initial-map derivative |
| `adjoint(xi,trajectory)` | `Evaluation<N,P>` containing that trajectory, cost, reduced gradient, raw nodal adjoints, and step multipliers |
| `evaluate(xi)` | Convenience call: `adjoint(xi,forward(xi))` |
| `cost_to_go_adjoint(xi,evaluation)` | Nodal adjoint with the preceding running-cost endpoint contribution removed |
| `hessian_vector(xi,v)` | Exact incremental-adjoint reduced Hessian product; recomputes its own trajectory and adjoint |
| `state_pullback(xi,trajectory,sources)` | Pulls arbitrary nodal state gradients back to controls through the discrete model, including its direct control dependence |
| `goal_gradient(xi,trajectory,node,state_gradient)` | A single goal pullback with a source at the requested node |
| `state_transition_operators(trajectory,xi)` | Local state-to-state maps \(M_k=B_k^{-1}C_k\) for every discrete theta step |

Pass a trajectory generated by the same system and controls. Size mismatches are rejected, but trajectories do not carry an ownership token and the API does not prove that arbitrary supplied states solve the intended problem. `adjoint` takes its trajectory by value so an rvalue can be moved into the result; passing an lvalue makes a copy. Derivative operations store a full time history. This first version prioritizes a clear numerical contract over checkpointing and factorization reuse.

## Optimization

`minimize(system,guess,options,admissible)` is a small dense reduced-space Newton driver. It assembles Hessian columns using exact Hessian-vector products. It requires second-order callbacks. If the Newton direction cannot be solved or is not a descent direction, it tries negative gradient. An Armijo line search safeguards step lengths and records accepted alpha values. The cost comparison includes a small floating-point allowance near stationarity.

`OptimizeOptions` specifies gradient tolerance, maximum iterations, and maximum backtracks. `Optimum<P>` returns the final controls and an iteration trace containing controls, cost, Euclidean gradient norm, and accepted step length. The final converged record has step length zero because no further step is taken. Convergence is based on the reduced gradient norm; exhausting iterations or failing a line search throws.

`minimize` targets interior unconstrained stationarity. Its control-admissibility predicate keeps trial points in the model domain; it does not make that method a constrained KKT solver, and a boundary optimum may fail the gradient test. Objective and callback exceptions are not silently treated as successful convergence. Domain errors in trial controls can trigger backtracking; other numerical failures terminate the solve.

`minimize_box_constrained(system,guess,lower,upper,options)` handles finite
componentwise bounds `lower[j] <= control[j] <= upper[j]`. It uses a dense
active-set Newton solve over free components, projection during line search, and
a projected-gradient fallback. The initial guess must be finite and feasible;
each lower bound must be strictly less than its upper bound. `BoxOptimum<P>`
returns the final control, raw reduced gradient, active lower/upper flags,
nonnegative multipliers, iteration trace, and KKT residual. The multiplier
convention is

$$
\nabla J-\lambda_L+\lambda_U=0,
\qquad \lambda_L,\lambda_U\ge0.
$$

At a lower-bound optimum, the raw gradient can be positive and is balanced by
the lower multiplier; at an upper-bound optimum, the raw gradient can be
negative and is balanced by the upper multiplier. The trace's `gradient_norm`
is the projected-gradient norm for this solver. See
[the bounded washout case](bounded_washout.md).

`minimize_path_constrained(system,guess,lower,upper,optimizer_options,
path_options)` handles a finite set of the system's nodal path inequalities and
control box bounds. It uses an augmented-Lagrangian outer iteration. Each path
term is inserted as a nodal cost, so the existing discrete adjoint computes its
reduced gradient and the incremental adjoint computes its exact Hessian action;
the path callback must provide `directional_gradient`. `PathConstrainedOptimum<P>`
returns nonnegative path multipliers, path values, box-bound multipliers, and
primal-feasibility, stationarity, complementarity, and combined KKT residuals.
This is a dense small-problem solver. A node-constrained result does not enforce
the inequality between nodes; refine the mesh and/or provide an interval
certificate for a continuous-time guarantee.

## Time-varying controls

The optimizer controls remain a fixed-size vector `Vec<P>`, but its entries can
represent coefficients of a prescribed temporal basis rather than constant
physical controls. `PiecewiseLinearControl<C,K>` interpolates `C` channels at
`K` strictly increasing times. `SoftmaxCompositionControl<S,K>` interpolates
`S-1` log-ratios at each knot and maps them through a stable softmax, so every
composition is positive and sums to one. It uses `P=(S-1)*K` optimized
coefficients.

`compose_dynamics(model,schedule.parameterization())` chains the instantaneous
control Jacobian and directional second derivatives through the schedule.
`compose_scalar_function` does the same for a running/terminal term, observation,
or path function that depends on the instantaneous control. Model callbacks
differentiate with respect to the instantaneous control; the schedule applies
the chain rule to the finite coefficient vector. Knot times and the time grid
are fixed during differentiation. A scheduled gas composition can therefore
use the same theta forward, tangent, discrete adjoint, and incremental-adjoint
machinery. Thermodynamic property callbacks still need accurate derivatives
with respect to state and instantaneous composition.

## Failure reporting and scope

Invalid grids, options, callback requirements, or trajectory sizes throw `std::invalid_argument`. Inadmissible initial/converged states and controls throw `std::domain_error`. Singular/poorly scaled linear systems, nonfinite numerical results, or failed iterations throw `std::runtime_error`. User callbacks may propagate their own exceptions. There is no global mutable solver state; callback thread safety remains the caller's responsibility.

Use `running` for a continuous cost density and `nodal` for isolated measurements. `cost_to_go_adjoint` removes only the preceding running-cost endpoint term; nodal sources remain, as appropriate for costs with observation jumps. General coupled equality constraints, general mass matrices, and DAEs remain extensions rather than implicit capabilities.

`qualification.hpp` supplies `local_goal_variance(gradient,covariance)`, `box_support(gradient,radii)`, `independent_gaussian_box_radius(sd,count,joint_probability)`, and `persistent_node_qualification(times,ratios)`. The last method returns an optional earliest sampled time after which all remaining ratios are at most one. It makes no between-node claim by itself. The decay example separately constructs whole-step physical envelopes and applies the suffix criterion to their interval start times, which yields a finite-horizon continuous-time statement under its stated model assumptions.

`all_goals_qualification_time(times)` combines individual persistent qualification times under a common horizon and criterion. It returns the latest time if every requested goal qualifies, including a valid zero time; it returns no time if any goal is unqualified. It does not combine node-only indicators with continuous physical envelopes or otherwise strengthen the input criterion.

`transient_analysis.hpp` provides state tangent/adjoint propagation, ordered
finite-horizon propagators, scaled worst-gain analysis, observation sensitivity
rows, Fisher information, and full-rank local covariance summaries. See
[`transient_flow_state_parameter_analysis.md`](transient_flow_state_parameter_analysis.md)
for assumptions and interpretation. `regularization.hpp` provides quadratic
prior, sensitivity-subspace, spatial H1/H2, and weak-constraint trajectory
penalty utilities; see
[`matrix_based_regularizers.md`](matrix_based_regularizers.md). These helpers
are additive to the existing strong-constraint system and examples.


## Amplification and bounded two-state qualification (0.7.0)

Include `<theta_opt/theta_opt.hpp>` or the individual headers.

- `control_to_state_jacobians(system, xi, trajectory)` assembles the state/control
  Jacobian at every node using the actual tangent solver. It includes nonlinear
  initial-condition maps and direct control dependence of the dynamics.
- `worst_amplification_2d(J, input_scales, output_scales)` computes the largest
  singular value of `Dout^-1 J Din` and returns the normalized scaled input,
  physical input direction, and physical response. Scales must be positive and
  finite. This is a two-dimensional utility, not a general SVD implementation.
- `ConvexPolygon2`, `rectangle_polygon`, `clip_halfspace`, `polygon_contains`,
  and `polygon_goal_bound` support a convex, counterclockwise two-dimensional
  initial-state domain. Containment requires a nondegenerate polygon; support
  takes the largest absolute linear goal deviation from a supplied center.
- `linear_quadratic_reconstruction(A,x,next,h,s)` reconstructs a constant-linear
  step at fraction `s` in [0,1], matching both endpoints and the left derivative.
- `linear_theta_quadratic_defect_bound(A,x,next,h,theta,norm)` bounds that
  reconstruction's residual in a caller-supplied physical norm, including the
  computed step residual and a floating-point allowance. Its contract is a
  constant linear system; it does not certify nonlinear dynamics.
- `propagate_defect_bound(E,D,h,mu)` applies the scalar comparison
  `exp(mu*h)*E + expm1(mu*h)/mu*D`, using `E+h*D` at zero. The caller must
  establish the logarithmic-norm upper bound; positive `mu` is supported.
- `quadratic_absolute_maximum(c0,c1,c2)` checks endpoints and an interior
  stationary point on [0,1].
- `qualify_intervals(edges,bounds,tolerances)` returns first-below, persistent,
  and joint times. Each supplied bound must cover its **entire interval**.
  A missing persistent time means qualification was not established through the
  final edge. The function cannot turn node samples into interval certificates.

These utilities use ordinary double precision. Conservative analytical bounds
and validation probes are not a formal outward-rounded interval proof. The
oscillator and two-tank examples share these geometry and qualification tools.


## Rank-deficient quadratic inverse problems (0.8.0)

`identifiability.hpp` supplies general fixed-size PSD spectral diagnostics and
a specialized two-control quadratic optimizer:

- `psd_spectrum<P>(H, relative_cutoff=1e-12)` now also supplies a dense
  fixed-dimension Jacobi eigensolver for symmetric positive-semidefinite
  information matrices of any compile-time control dimension `P`. It returns
  descending eigenvalues, orthonormal directions, numerical rank, and the
  absolute cutoff. Symmetry and material indefiniteness are checked. The
  threshold is relative to the largest eigenvalue; users remain responsible
  for scaling heterogeneous control coordinates. This reports numerical local
  information and is not a global-identifiability proof.
- `PsdSpectrum<P>::solve_compatible(rhs, relative_tolerance=1e-10)` applies the
  inverse on retained eigen-directions and rejects unresolved null-space
  components. As with the two-dimensional API, it does not interpret null
  directions as zero uncertainty.

- `psd_spectrum_2d(H, relative_cutoff=1e-12)` returns descending eigenvalues,
  orthonormal directions, numerical rank, and the absolute cutoff. It rejects
  nonfinite, nonsymmetric, and materially indefinite matrices. Coordinate
  scaling is the caller's responsibility; rank is relative to the largest
  information eigenvalue. A zero matrix has rank zero.
- `PsdSpectrum2::solve_compatible(rhs, relative_tolerance=1e-10)` applies the
  inverse on retained directions. It rejects right-hand sides with unresolved
  null components instead of silently discarding them. This operation does not
  turn the pseudoinverse into a covariance with known zero null uncertainty.
- `minimize_quadratic_psd_2d(system, seed, gradient_tolerance=1e-8,
  relative_cutoff=1e-12)` returns `QuadraticOptimum2`, containing an ordinary
  `Optimum<2>` plus spectral diagnostics. Its contract is an unconstrained convex
  quadratic reduced objective. It uses the system adjoint gradient and exact
  incremental-adjoint Hessian columns, takes a compatible pseudoinverse Newton
  correction, and preserves the seed's Euclidean null component. It verifies
  stationarity, nonincreasing cost within roundoff, and endpoint Hessian
  agreement. Those checks do not establish global quadratic structure for an
  arbitrary supplied function. Admissibility is governed by the forward model;
  this routine does not solve boundary KKT conditions.

No null-space penalty is added. A representative selected by the seed is not
additional measurement information. Use bounded feasible sets or an explicitly
stated prior to qualify goals affected by missing directions. See
[the insulated thermal experiment](insulated_thermal.md).


## Disconnected bounded feasible sets (0.9.0)

`bounded_sets.hpp` adds ordinary double-precision conservative range tools for
small two-control problems. These are not outward-rounded interval arithmetic.

- `Box2 {lower,upper}` describes a closed finite box; `validate_box` and
  `contains` check its ordering and membership.
- `Quadratic2 {constant,linear,quadratic}` represents `c+bᵀx+xᵀQx`.
  `evaluate` evaluates it. `quadratic_range(q,box)` uses a centered expansion
  and absolute coefficient bounds, with roundoff padding. Convexity is not
  required; ranges can overestimate the actual image.
- `polynomial_quadratic_range(coefficients,box)` bounds the full polynomial
  `sum s^j*q_j(x)` for `s` in [0,1], not just its values at time nodes.
- `subdivide_feasible_boxes(prior,target_width,possibly_feasible,budget)` returns
  `BoxCover2`, including boxes, visitation/rejection counts, and a
  `resolution_complete` flag. The callback must return true for any box that
  may contain an admissible point. Rejection soundness is the caller's
  responsibility. Unresolved boxes at a budget or floating-point split limit
  are retained and flagged incomplete. Empty output means every box was
  rejected under that callback's range argument, not that Newton failed.
- `cover_contains` tests membership in the retained union.
  `cover_quadratic_range` and `cover_polynomial_quadratic_range` bound goals
  over **all** retained boxes, without choosing an optimizer branch. Empty
  covers have no goal enclosure and cause these union range functions to throw.
- `absolute_bound` converts an ordered range to a worst absolute magnitude.

Model-specific numerical allowances and observation construction remain the
caller's responsibility. The electrothermal example supplies electrical and
thermal defect propagation, widens observation bounds, then uses the shared
union/time-polynomial bounds for persistent qualification. A selected-region
calculation is retained only as a diagnostic of the error caused by dropping a
feasible polarity. Full-rank local information does not establish a single
feasible region. See [electrothermal.md](electrothermal.md).
