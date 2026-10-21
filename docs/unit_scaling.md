# Unit and scaling invariance in drug-washout state reconstruction

## Purpose

This experiment checks that the discrete optimality system describes the same
physical inference problem when its variables are expressed in different
units. We use the analytically solvable one-compartment washout problem and run
the same fixed observations in two representations:

| Representation | State and initial control | Time | Elimination rate |
|---|---|---|---|
| A | mg/L | h | $0.45\;\mathrm{h}^{-1}$ |
| B | microgram/L | min | $0.0075\;\mathrm{min}^{-1}$ |

The conversions are $1\;\mathrm{mg/L}=1000\;\mathrm{microgram/L}$ and
$1\;\mathrm{h}=60\;\mathrm{min}$. This is a numerical consistency test, not a
new data set: observations, their error bounds, prior, and engineering goals
are all transformed from one physical experiment. The test exercises the
forward theta solution, nodal observation objective, tangent, discrete adjoint,
Hessian-vector product, optimization, running-cost AUC, uncertainty bounds, and
persistent qualification times.

## Physical model and its transformed form

The model in representation A is

$$
\frac{dc}{dt}=-k c,\qquad c(t_0)=c_0,
\quad k=0.45\;\mathrm{h}^{-1},\quad t_f=8\;\mathrm{h}.
$$

Let $a=1000$ be the concentration scale and $b=60$ the time scale, and define
$\widetilde c=a c$ and $\tau=b t$. Then

$$
\frac{d\widetilde c}{d\tau}=-\widetilde k\widetilde c,
\qquad \widetilde k=\frac{k}{b}=0.0075\;\mathrm{min}^{-1},
\qquad \widetilde c(\tau_0)=a c_0.
$$

The analytic solution is $c(t)=c_0e^{-kt}$ in either representation. The
discrete theta update depends on the product $h k$; because
$\widetilde h\widetilde k=(b h)(k/b)=hk$, both unit systems have the same
dimensionless update and the same number of steps.

The observations are fixed at $0.5,1,2,4,6,8$ hours, with the same deterministic
offsets as the original washout experiment. In representation B, sample times
are multiplied by 60, concentrations and standard deviations by 1000, and the
same prior is expressed as $[0,20000]$ microgram/L. The concentration goal
tolerance is multiplied by 1000. The remaining-AUC tolerance and values are
multiplied by $1000\cdot60$, since

$$
\widetilde A(\tau)=\int_{\tau}^{\tau_f}\widetilde c(s)\,ds
=ab\,A(t).
$$

## How the objective and derivatives must transform

For the same six observations $y_j$, standard deviation $\sigma$, and model
predictions $c_j$, the objective is

$$
J(c_0)=\frac12\sum_j\left(\frac{c_j-y_j}{\sigma}\right)^2.
$$

Under the unit conversion, residuals and standard deviations are both
multiplied by $a$, so every standardized residual and the entire objective are
unchanged. The parameter changes from $c_0$ to $\widetilde c_0=a c_0$.
Consequently, gradients and Hessians are covectors and transform as

$$
\frac{\partial J}{\partial\widetilde c_0}
=\frac1a\frac{\partial J}{\partial c_0},
\qquad
\frac{\partial^2J}{\partial\widetilde c_0^2}
=\frac1{a^2}\frac{\partial^2J}{\partial c_0^2}.
$$

The nodal adjoint has units of objective per state, so it also scales by
$1/a$. The tangent $\partial c(t)/\partial c_0$ is dimensionless and should
remain unchanged. For the AUC, its value scales by $ab$, while its derivative
with respect to the initial concentration scales by $b$:

$$
\widetilde A=ab A,\qquad
\frac{\partial\widetilde A}{\partial\widetilde c_0}
=b\frac{\partial A}{\partial c_0}.
$$

These transformations are checked directly against the library's discrete
forward, adjoint, and incremental-adjoint results. Comparing raw gradient or
Hessian numbers without applying their unit transformations would be an
incorrect invariance test.

## Validation and results

For $N=32,128,512$ midpoint steps, the C++ program compares converted physical
trajectories and fitted initial concentrations, checks objective invariance,
checks gradient/Hessian/costate scaling, and compares AUC and its initial-state
sensitivity after conversion. Qualification uses the same bounded observations
and error decomposition in both representations; times from minutes are
converted back to hours before comparing decisions.

At $N=512$, the fitted initial concentration is $10.0000695$ mg/L in
representation A and $10000.0695$ microgram/L in representation B. The
dimensionless objective, scaled gradient and Hessian, AUC, and AUC sensitivity
agree after conversion to within roundoff (the largest reported relative
difference is below $3\times10^{-15}$). Both unit systems give concentration
qualification at $3.71875$ h and remaining-AUC qualification at $2.609375$ h.
The converted minute-based times are $223.125$ min and $156.5625$ min,
respectively.

The same convergence study is performed in both representations, with errors
normalized back to mg/L-hours. At $\theta=0.5$, the observed orders for state,
tangent, cost, gradient, Hessian, fit, and AUC are approximately 2. At
$\theta=0.65$, all are approximately 1. The unit change therefore preserves
the discretization errors and their convergence rates, not just the final
trajectory.

The executable performs more than 2,000 checks. Many compare every time node
across three meshes; the rest verify the analytic unit conversions and both
theta-method convergence regimes. These are algebraic and numerical
consistency checks for this controlled scalar model. They do not guarantee
that arbitrary poorly scaled problems will be well conditioned for a finite
precision optimizer. In practice, optimization variables and residuals should
still be nondimensionalized or scaled to sensible magnitudes.

## Reproduce and inspect

```bash
make unit-scaling-plots
```

The C++ executable writes the converted values, errors, observed orders, and
qualification times to `results/unit_scaling/`. Python reads those CSV files
and only makes plots; it performs no modeling, fitting, differentiation, or
qualification.

Figures:

- `trajectory_unit_conversion.png`: physical state curves after unit conversion
  and their roundoff-scale difference;
- `unit_transformed_quantities.png`: invariance of the fitted state, objective,
  scaled gradient and Hessian, AUC, and AUC sensitivity;
- `unit_scaling_convergence.png`: observed order for both unit systems and both
  theta values;
- `unit_scaling_qualification.png`: identical physical qualification times
  over mesh refinement.
