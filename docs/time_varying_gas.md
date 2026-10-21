# Time-varying gas composition in a mixing vessel

This example asks whether thermodynamic measurements can recover a changing
gas composition when the inlet composition varies with time. It extends the
static Peng–Robinson composition fit into a transient state-finding problem.
The synthetic experiment is small enough to have an analytical reference
solution while exercising the C++20 discrete optimality system and the
ThermoGPU Peng–Robinson implementation.

## Physical model

Consider an ideally mixed vessel with constant total inventory and turnover
rate $k=0.7\,\mathrm{min}^{-1}$. Its outlet composition equals its well-mixed
state. The inlet contains a known nitrogen fraction of 0.05; ethane is the
unknown varying feed fraction, and methane makes up the remainder. Let
$x_E(t)$ and $x_N(t)$ be the vessel ethane and nitrogen mole fractions, and
$u_E(t)$ the ethane mole fraction in the feed. Then

$$
\dot{x}_E=k(u_E-x_E),\qquad
\dot{x}_N=k(0.05-x_N),\qquad
x_M=1-x_E-x_N.
$$

The unknown initial state and three piecewise-linear feed knots at $0,5,10$
min form the five fitted coefficients:

$$
\boldsymbol{z}=(x_E(0),x_N(0),u_E(0),u_E(5),u_E(10)).
$$

The synthetic truth is $(0.08,0.04,0.06,0.18,0.08)$. Twenty-one grab
measurements, spaced by 0.5 min over $[0,10]$ min, provide mixture density.
For measurement diversity, sample conditions cycle through four known
pressure/temperature pairs: (1 MPa, 300 K), (5 MPa, 300 K), (10 MPa, 350 K),
and (20 MPa, 400 K). Density is evaluated using the existing ThermoGPU
Peng–Robinson implementation. Each synthetic observation has a stated
relative scale of 0.1%; the values are noiseless so discretization and
inversion errors can be inspected without random noise.

## Exact reference and discrete estimation

On an interval where the feed is linear, $u_E(t+s)=u_0+rs$, ethane has the
closed-form update

$$
x_E(t+h)=u_0+r\left(h-\frac{1}{k}\right)
+\left(x_E(t)-u_0+\frac{r}{k}\right)e^{-kh}.
$$

Nitrogen has the exact update

$$
x_N(t+h)=0.05+(x_N(t)-0.05)e^{-kh}.
$$

These formulas generate the reference trajectory and synthetic density
observations. The inverse fit uses the library's discrete theta method, with
observation costs attached at their matching time nodes. The same discrete
system supplies the tangent linear model and adjoint gradient; the gradient
is checked by finite differences. Refinement with $\theta=0.5$ should show
second-order state convergence, while $\theta=0.65$ should show first-order
convergence. Fitted coefficient errors show how discretization error enters
the recovered initial composition and feed schedule. At the finest grid, the
local information matrix has rank 5 of 5 and condition number about 390 in
the reported coordinates. This is a local, linearized diagnostic; it does
not establish a global uniqueness or a bounded-data qualification result.

The C++ program performs the forward solve, thermodynamic evaluation,
optimization, derivatives, and convergence calculations. Python only renders
the CSV files produced by C++.

## Interpretation and limits

This is a dynamic composition-reconstruction test, not a full pipeline or
compressible vessel simulator. Constant total inventory and turnover rate are
assumed known. Only the inlet ethane fraction changes; nitrogen in the feed is
fixed. Density depends nonlinearly on composition through Peng–Robinson, and
is the sole measured property in this experiment. The three-knot schedule is
a modeling assumption that makes the control finite dimensional; recovery
does not imply that one thermodynamic property can identify an arbitrary
feed history.

The 0.1% scale weights the least-squares objective. The example does not
compute conservative state-qualification bounds, account for uncertain
vessel parameters, or model measurement-noise realizations. Optimization
success and convergence checks are numerical validation, not an engineering
qualification verdict.

## Run

Run: make time-varying-gas-plots

This writes CSV results and creates:

- composition_reconstruction.png: true and reconstructed vessel mole fractions.
- feed_reconstruction.png: true and fitted inlet ethane fraction for both theta values.
- convergence_and_parameter_error.png: state convergence and joint coefficient error.
- density_residuals.png: fitted density residuals scaled by the observation scale.
- information_spectrum.png: local observation information across the five fitted coefficients.
