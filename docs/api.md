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

Pass a trajectory generated by the same system and controls. Size mismatches are rejected, but trajectories do not carry an ownership token and the API does not prove that arbitrary supplied states solve the intended problem. `adjoint` takes its trajectory by value so an rvalue can be moved into the result; passing an lvalue makes a copy. Derivative operations store a full time history. This first version prioritizes a clear numerical contract over checkpointing and factorization reuse.

## Optimization

`minimize(system,guess,options,admissible)` is a small dense reduced-space Newton driver. It assembles Hessian columns using exact Hessian-vector products. It requires second-order callbacks. If the Newton direction cannot be solved or is not a descent direction, it tries negative gradient. An Armijo line search safeguards step lengths and records accepted alpha values. The cost comparison includes a small floating-point allowance near stationarity.

`OptimizeOptions` specifies gradient tolerance, maximum iterations, and maximum backtracks. `Optimum<P>` returns the final controls and an iteration trace containing controls, cost, gradient infinity norm, and accepted step length. The final converged record has step length zero because no further step is taken. Convergence is based on the reduced gradient norm; exhausting iterations or failing a line search throws.

The optimizer targets interior unconstrained stationarity. A control-admissibility predicate keeps trial points in the model domain; it does not make the method a constrained KKT solver, and a boundary optimum may fail the gradient test. Objective and callback exceptions are not silently treated as successful convergence. Domain errors in trial controls can trigger backtracking; other numerical failures terminate the solve.

## Failure reporting and scope

Invalid grids, options, callback requirements, or trajectory sizes throw `std::invalid_argument`. Inadmissible initial/converged states and controls throw `std::domain_error`. Singular/poorly scaled linear systems, nonfinite numerical results, or failed iterations throw `std::runtime_error`. User callbacks may propagate their own exceptions. There is no global mutable solver state; callback thread safety remains the caller's responsibility.

Use `running` for a continuous cost density and `nodal` for isolated measurements. `cost_to_go_adjoint` removes only the preceding running-cost endpoint term; nodal sources remain, as appropriate for costs with observation jumps. Differentiable mesh adaptation, inequality constraints, general mass matrices, and DAEs remain extensions rather than implicit capabilities.

`qualification.hpp` supplies `local_goal_variance(gradient,covariance)`, `box_support(gradient,radii)`, `independent_gaussian_box_radius(sd,count,joint_probability)`, and `persistent_node_qualification(times,ratios)`. The last method returns an optional earliest sampled time after which all remaining ratios are at most one. It makes no between-node claim by itself. The decay example separately constructs whole-step physical envelopes and applies the suffix criterion to their interval start times, which yields a finite-horizon continuous-time statement under its stated model assumptions.

`all_goals_qualification_time(times)` combines individual persistent qualification times under a common horizon and criterion. It returns the latest time if every requested goal qualifies, including a valid zero time; it returns no time if any goal is unqualified. It does not combine node-only indicators with continuous physical envelopes or otherwise strengthen the input criterion.
