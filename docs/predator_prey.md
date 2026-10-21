# Predator–prey inverse problem: initial state, rates, and joint recovery

This example applies the same theta-consistent discrete optimality system to a
two-state nonlinear model in three progressively harder inverse problems. It
asks which quantities the observations determine when (1) the kinetic rates
are known, (2) the initial populations are known, and (3) both are unknown.
All three fits use the C++20 library's theta forward solve, nodal observation
cost, tangent linear model, discrete adjoint, exact incremental-adjoint Hessian
products, and box-constrained optimizer.

## Model and analytical check

The dimensionless Lotka–Volterra equations are

$$
\dot{x}=a x-bxy,\qquad \dot{y}=cxy-dy,
$$

where $x$ and $y$ are normalized prey and predator populations. The synthetic
truth is

$$
(x(t_0),y(t_0))=(1.6,0.6),\qquad (a,b,c,d)=(1.0,0.8,0.6,0.7),
\quad t\in[0,10].
$$

For positive populations and positive rates, the model remains in the positive
quadrant. It also has the conserved quantity

$$
H(x,y)=cx-d\log x+by-a\log y.
$$

The reference trajectory is generated independently with fixed-step RK4 at
$h=2\times10^{-4}$. Its maximum drift in $H$ is printed and tested. This is a
strong analytical invariant check, although the classical Lotka–Volterra model
does not have an elementary explicit formula for $x(t)$ and $y(t)$.

## Observation experiment

Both populations are observed every 0.5 time units, giving 20 paired assays.
Small fixed perturbations are added to the reference states. The inverse cost
is the dimensionless least-squares objective

$$
J=\frac12\sum_{j=1}^{20}\left[
 \frac{(x(t_j)-x_j^{\rm obs})^2}{\sigma_x^2}+
 \frac{(y(t_j)-y_j^{\rm obs})^2}{\sigma_y^2}\right],
\qquad \sigma_x=0.025,\quad \sigma_y=0.020.
$$

These are nodal observation costs; they are not multiplied by running-cost
quadrature weights. Synthetic perturbations are fixed across mesh refinements,
so changes in the fitted controls measure the effect of the discrete model,
not a regenerated data sample.

## Three inverse problems

| Case | Unknown control vector | Quantities treated as known |
|---|---|---|
| Initial-state estimate | $(x(t_0),y(t_0))$ | $(a,b,c,d)$ |
| Rate estimation | $(a,b,c,d)$ | $(x(t_0),y(t_0))$ |
| Joint estimate | $(x(t_0),y(t_0),a,b,c,d)$ | Observation schedule and model form |

Each unknown is constrained by simple positive componentwise bounds, handled by
the library's box-constrained solver. The positivity of the Lotka–Volterra
trajectory follows from positive initial populations and rates. The dynamics'
admissibility callback rejects a nonpositive Newton state; it is a domain
safeguard, not a path-constraint KKT solver. General path-constraint KKT
machinery is therefore not needed to run these three estimates. A future case
with an active minimum-biomass, capacity, or harvest constraint would be an
appropriate reason to add it.

## Derivative and convergence checks

For every inverse problem, the example compares the tangent trajectory with
finite differences of the forward solve, the discrete adjoint gradient with
central finite differences of the reduced objective, and the
incremental-adjoint Hessian with finite differences of the adjoint gradient.
The tests exercise state-state, state-rate, and rate-rate second derivatives
in a coupled two-population system.

The mesh sequence is $N=40,80,160,320$, with all assay times on mesh nodes. The
example reports trajectory errors against an RK4 integration using the same
fitted controls, and changes in the fitted controls relative to the finest
mesh. It runs both midpoint ($\theta=0.5$) and a first-order theta value
($\theta=0.65$). The expected observed orders are two and one, respectively;
the C++ program checks those rates for all three inverse problems.

At $N=320$, midpoint recovers the synthetic truth to within the fixed assay
noise and discretization error. The six-control joint fit also has full rank in
the local Gauss–Newton information matrix; its eigenvalue condition number is
reported as a scaling-dependent local diagnostic, not as a global uniqueness
claim. The exact estimates and all derivative, invariant, rank, and convergence
checks are written to `results/predator_prey/`.

## Build and reproduce

```bash
make predator-prey
make predator-prey-plots
make test
```

The C++ executable performs all simulation, fitting, differentiation, and
validation. The Python script only plots its CSV outputs:

- `state_reconstruction.png`: both populations and observations for all three
  fits;
- `phase_portrait.png`: coupled prey–predator trajectories and paired assays;
- `convergence_and_derivatives.png`: theta refinement and derivative checks.

## What this adds to the test program

The logistic-growth example isolates one nonlinear state with elementary exact
time-dependent references. This example moves to a coupled two-state system
with an invariant and one shared set of simultaneous population observations.
It demonstrates that initial-state-only, parameter-only, and joint estimation
can all use the same discrete optimality library. The optimizer and derivative
system already supported six controls. The case did require a new reusable
library feature: `psd_spectrum<P>` now diagnoses local information rank and
conditioning beyond the previous two-control limit. General path constraints
remain a separate extension until an experiment has a meaningful active path
bound.
