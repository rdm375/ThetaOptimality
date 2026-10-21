# Nonlinear dynamics: logistic biomass growth

This experiment extends the state-finding workflow from linear dynamics to a
nonlinear state equation whose exact solution is known. Its purpose is to test
whether the library differentiates the **discrete theta solution and discrete
objective correctly when the dynamics have a state-dependent Jacobian**. We
compare the state, tangent, adjoint, gradient, Hessian, and recovered initial
biomass with analytical references. Midpoint (`theta = 0.5`) should converge at
second order; `theta = 0.65` supplies the expected first-order comparison.

The model is a batch-culture biomass concentration. Its growth rate and carrying
capacity are known; only the initial biomass is inferred. The observations are
synthetic, fixed in physical time, and generated from the exact solution with a
small fixed assay perturbation. The theta solver—not the analytical solution—
produces the states used by the optimizer.

Run the case and make its figures with:

```bash
make logistic-growth-plots
```

## 1. Continuous model and exact solution

For $t\in[t_0,t_f]=[0,4]$ h, biomass concentration $x(t)$ in g/L obeys

$$
\dot x=f(x)=r x\left(1-\frac{x}{K}\right),
\qquad x(t_0)=\xi,
$$

with known specific growth rate $r=0.8\;\mathrm{h}^{-1}$ and carrying
capacity $K=8\;\mathrm{g/L}$. For elapsed time $\tau=t-t_0$,

$$
x(t_0+\tau;\xi)
=\frac{K\xi e^{r\tau}}{K+\xi(e^{r\tau}-1)}.
$$

This solution is nonlinear in the initial state. Its first two initial-state
sensitivities are

$$
s(\tau;\xi)=\frac{\partial x}{\partial\xi}
=\frac{K^2e^{r\tau}}{[K+\xi(e^{r\tau}-1)]^2},
$$

$$
z(\tau;\xi)=\frac{\partial^2x}{\partial\xi^2}
=-\frac{2K^2e^{r\tau}(e^{r\tau}-1)}
        {[K+\xi(e^{r\tau}-1)]^3}.
$$

The state Jacobian and its directional derivative are

$$
f_x(x)=r\left(1-\frac{2x}{K}\right),
\qquad Df_x(x)[v]=-\frac{2r}{K}v.
$$

The second formula is the nonlinear-dynamics term used by the incremental
adjoint. A constant-Jacobian model would not test that term.

## 2. Observations and inverse problem

The unknown is normalized initial biomass $u_0=\xi/K$, constrained to the
interior interval $0.02<u_0<0.60$. Five biomass assays are taken at elapsed
times $0.25$, $0.75$, $1.5$, $2.5$, and $3.75$ h. Each value is generated as

$$
y_j=x(t_0+\tau_j;\xi_{\rm true})+\sigma\eta_j,
\qquad \xi_{\rm true}=1.6\;\mathrm{g/L},\quad
\sigma=0.08\;\mathrm{g/L},
$$

where the fixed standardized perturbations are
$(-0.30,0.25,-0.20,0.10,0.15)$. The same data and scales are used on every
mesh. No assay is at $t_0$, so the initial state is inferred through nonlinear
growth rather than copied directly from an initial observation.

The analytical reference objective is the isolated-observation least-squares
functional

$$
J(\xi)=\frac12\sum_{j=1}^{5}
  \left(\frac{x(t_0+\tau_j;\xi)-y_j}{\sigma}\right)^2.
$$

The library represents the same quantity as nodal costs. Nodal costs are not
multiplied by time-step weights: they represent individual assays, not a
continuous measurement density. The continuous analytical optimum is found by
bracketing the scalar root $J'(\xi)=0$; the discrete solve uses the library's
forward theta method, adjoint gradient, incremental-adjoint Hessian, Newton
optimizer, and line search.

Writing $r_j=x(t_0+\tau_j;\xi)-y_j$,

$$
J'(\xi)=\sum_j \frac{r_j s_j}{\sigma^2},\qquad
J''(\xi)=\sum_j\frac{s_j^2+r_jz_j}{\sigma^2}.
$$

Because the library control is $u_0=\xi/K$, the control derivatives used in
the checks are $dJ/du_0=KJ'(\xi)$ and $d^2J/du_0^2=K^2J''(\xi)$.

For times strictly between observation times, the continuous nodal-observation
costate is

$$
\lambda(t)=\sum_{j:t<t_j}
 \frac{x(t_j;\xi)-y_j}{\sigma^2}
 \frac{s(t_j-t_0;\xi)}{s(t-t_0;\xi)}.
$$

The strict inequality avoids comparing different one-sided values at an
observation jump. With no terminal observation, $\lambda(t_f)=0$.

## 3. Discrete optimality checks

The code evaluates a fixed probe state on six nested uniform meshes,
$N=32,64,128,256,512,1024$, for each of $\theta=0.5$ and $\theta=0.65$.
Observation times lie exactly on every grid. It compares:

1. The theta trajectory with the exact logistic state.
2. The discrete tangent with $K s(t;\xi)$ times the chosen control direction.
3. The discrete adjoint away from its point-source jumps with the analytical
   costate.
4. The discrete cost, gradient, and Hessian-vector product with their
   analytical objective references.
5. The recovered discrete minimizer with the continuous analytical minimizer.

The tangent/adjoint pairing is independently checked. Central differences also
check the analytical objective gradient/Hessian and the discrete gradient and
incremental-adjoint Hessian. The output reports observed orders from the finest
mesh pair for all seven quantities. The expected pattern is order two for
midpoint and order one for $\theta=0.65$.

## 4. Bounded-data qualification and numerical allowance

Each assay is also interpreted with a deterministic error bound of
$0.04\;\mathrm{g/L}$ and a prior interval $\xi\in[0.4,4.5]\;\mathrm{g/L}$.
The logistic formula can be inverted analytically:

$$
\xi=\frac{Kx}{x+(K-x)e^{r\tau}}.
$$

Inverting every assay band and intersecting the resulting intervals gives one
nonlinear data-consistent initial-state interval. The exact logistic flow maps
its endpoints to the full-time feasible state envelope because the solution is
monotone in $\xi$.

For the theta fit, the code computes a forward-error allowance from the
one-step discrepancy between the discrete step and the exact logistic flow
started at the discrete node. The exact flow map has derivative bounded by
$e^{rh}$ for nonnegative states, so these local discrepancies are propagated
through the mesh. The largest nodal allowance is added conservatively to each
whole-interval feasible-state spread. The exact difference between two
logistic trajectories has at most one interior maximum; the code locates it
from $x(t;\xi_a)+x(t;\xi_b)=K$ and passes whole-step bounds to the shared
persistent qualification routine.

At $N=32$ the combined state bound does not meet the $0.02\;\mathrm{g/L}$
tolerance over the requested horizon. At $N=128$, it qualifies persistently
from 2.75 h. This is a controlled illustration of temporal error affecting an
engineering conclusion even when the fitted initial biomass is already close
to the analytical optimum.

## 5. Retained outputs and scope

The C++ executable writes convergence, objective-derivative, bounded-data,
qualification-interval, observation, and check CSV files under `results/logistic_growth/`.
Python reads those files and generates:

- `state_qualification.png`: truth, discrete reconstruction, feasible envelope,
  and qualification time for coarse and refined meshes.
- `theta_convergence.png`: state, tangent, adjoint, cost, gradient, Hessian,
  and recovered-initial-state errors.
- `objective_derivatives.png`: analytical and theta-discrete objective,
  gradient, and Hessian curves.

This case tests nonlinear **dynamics** with a single unknown initial state and
known parameters. It does not estimate growth rate or carrying capacity, and it
does not impose path constraints. A predator-prey model belongs after this
single-state nonlinear validation: it adds coupled nonlinear states, phase
interaction, and potentially weak or ambiguous observation designs. We should
keep its parameters known initially, compare against a high-accuracy or
analytical reference where available, and only then consider estimating more
than the initial state.
