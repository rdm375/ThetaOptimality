# Drug washout: measurement accuracy versus numerical accuracy

## Purpose of the experiment

This example tests a distinction that is easy to lose in a state reconstruction:
how much uncertainty comes from the sensors, and how much comes from solving the
model on a finite time grid? We hold six plasma-concentration observations fixed,
refine the theta-method grid, and then widen the allowed measurement errors
without changing those observations. We ask when two engineering quantities are
established to their tolerances:

1. plasma concentration, in mg/L;
2. remaining exposure (AUC from the current time to the end of the horizon), in
   mg h/L.

The example is deliberately small. It has a closed-form state, sensitivities,
least-squares estimate, and AUC, so every component of the discrete optimality
system can be compared with an independent analytical reference. It is a
verification and teaching case, not a pharmacokinetic dosing model: elimination
is first order, the elimination rate is known, and there is no absorption,
distribution between compartments, or parameter uncertainty.

## Physical model and exact reference

After an idealized intravenous bolus at elapsed time $t_0=0$ h, the one-compartment
concentration satisfies

$$
\dot c(t)=-k c(t),\qquad c(0)=c_0,\qquad
k=0.45\;\mathrm{h}^{-1},\quad 0\le t\le t_f=8\;\mathrm{h}.
$$

Here $c_0$ is the unknown initial concentration. The analytic state and its
initial-state sensitivity are

$$
c(t;c_0)=c_0e^{-kt},\qquad
\frac{\partial c(t;c_0)}{\partial c_0}=e^{-kt}.
$$

The engineering exposure goal is the AUC still remaining after time $t$:

$$
A(t;c_0)=\int_t^{t_f} c(s;c_0)\,ds
=\frac{c_0}{k}\left(e^{-kt}-e^{-kt_f}\right).
$$

Its initial-state sensitivity is

$$
\frac{\partial A(t;c_0)}{\partial c_0}
=\frac{1}{k}\left(e^{-kt}-e^{-kt_f}\right).
$$

These formulas are the references for the C++ forward state, tangent, objective,
gradient, Hessian, theta-quadrature exposure, and recovered initial state.

## Fixed observations and measurement model

The six observation times are $0.5,1,2,4,6,8$ h. A fixed data vector is
generated from the analytical state at the true value $c_0=10$ mg/L with small,
deterministic sensor offsets:

| Time (h) | Offset (mg/L) |
|---:|---:|
| 0.5 | +0.0010 |
| 1 | -0.0015 |
| 2 | +0.0007 |
| 4 | -0.0008 |
| 6 | +0.0012 |
| 8 | -0.0004 |

Exactly this same vector is used for every grid and for both measurement-error
scenarios. The two deterministic per-observation error bounds are
$\varepsilon=0.005$ mg/L and $\varepsilon=0.050$ mg/L. Both contain every
offset above. The prior is $0\le c_0\le20$ mg/L.

The reconstruction objective is the isolated-observation least-squares cost

$$
J_h(c_0)=\frac12\sum_{j=1}^{6}
\left(\frac{c_{h}(t_j;c_0)-y_j}{\varepsilon}\right)^2,
$$

where $c_h$ is produced by the library's theta method. Each observation is a
**nodal cost**: it is not multiplied by a time-step or theta-quadrature weight.
Changing $\varepsilon$ changes the objective scale but, because all six bounds
are the same and the objective is unconstrained at its interior solution, does
not change the fitted value for this fixed dataset. It does change the
observation-consistent set and therefore the engineering uncertainty.

For the continuous analytical observation model, write $s_j=e^{-kt_j}$. The
least-squares minimizer and exact bounded-data interval are

$$
\widehat c_0=\frac{\sum_j s_jy_j}{\sum_j s_j^2},\qquad
c_{0,\min}=\max\!\left(0,\max_j\frac{y_j-\varepsilon}{s_j}\right),\qquad
c_{0,\max}=\min\!\left(20,\min_j\frac{y_j+\varepsilon}{s_j}\right).
$$

The C++ program computes these exact analytical references and independently
checks the theta-discrete least-squares minimizer against its closed-form
discrete counterpart.

## Discrete optimality system exercised

For each fixed grid, the state is advanced by the library's theta method. The
observation objective is assembled from nodal costs at the six sample nodes.
The library then supplies the discrete tangent, reverse discrete adjoint,
reduced gradient, and incremental-adjoint Hessian action. A separate running
cost system integrates $c(t)$ with the same theta weights to compute discrete
total exposure; summing its remaining interval contributions gives the discrete
remaining-AUC values used by this goal.

The validation checks compare the forward state and tangent to their exact
formulas; the objective, adjoint gradient, Hessian action, and recovered initial
state to analytical references; and the discrete adjoint AUC sensitivity to the
theta-quadrature result. The convergence study uses the same six observations at
every resolution. It measures cost and its derivatives at a fixed diagnostic
initial value $c_0=11$ mg/L. This avoids an accidental near-zero-residual
cancellation that occurred when measuring objective error only at the true
initial state. The state and tangent errors are still measured from the true
initial state.

Midpoint ($\theta=0.5$) should converge at second order. The non-midpoint
comparison $\theta=0.65$ should converge at first order. Current output reports
the observed final-pair orders for state, tangent, cost, gradient, Hessian,
recovered initial concentration, and theta-quadrature AUC.

## Separating the accuracy budgets

At each interval the C++ calculation retains two contributions:

- **Measurement contribution:** the exact analytical state/AUC spread induced
  by the feasible interval for $c_0$ from the bounded sensor data.
- **Numerical contribution:** the discrepancy between the analytical fit and
  the theta-discrete fit, including the piecewise-linear reconstruction between
  time nodes used for the whole-interval decision.

The plotted and reported total is the sum of these two bounds. This is a
controlled analytic-reference decomposition for this benchmark. It is not a
general-purpose certified discretization estimator: exact knowledge of this
model provides the numerical reference. The interval between nodes is assessed
using linear interpolation of the nodal theta solution, and the maximum
discrepancy against the exponential reference is checked over each interval,
including any interior extremum.

Qualification is persistent: a goal qualifies at time $t$ only if its summed
bound remains at or below tolerance on every subsequent interval through 8 h.
The concentration tolerance is 0.012 mg/L; the remaining-AUC tolerance is
0.040 mg h/L. The joint decision requires both goals to qualify, so its start
time is the later of their individual persistent times.

The experiment contrasts three effects:

1. On a coarse grid with precise sensors, numerical error can dominate even
   though the data uncertainty is small.
2. Refining the grid reduces the numerical contribution while leaving the data
   and measurement contribution fixed.
3. With coarse sensor bounds, measurement uncertainty limits qualification
   after the numerical contribution has become small.

These cases demonstrate why increasing the numerical resolution cannot repair
poor measurements, and why improving sensors does not repair an under-resolved
numerical model.

## Reproduce and inspect

```bash
make drug-washout
make drug-washout-plots
```

`make test` includes the executable. The C++ program writes its numerical and
qualification results to `results/drug_washout/`. Python only renders those
already-computed CSV outputs; it performs no fitting, error analysis, or
qualification decisions.

The figures are:

- `reconstruction_and_qualification.png`: analytical truth, analytical fit,
  theta-discrete fit, total uncertainty allowance, and qualification times for
  concentration and remaining AUC;
- `error_budget_decomposition.png`: measurement, numerical, and summed error
  contributions in units of each engineering tolerance;
- `convergence.png`: the seven observed convergence quantities for both theta
  values;
- `qualification_times.png`: persistent qualification times across grid and
  sensor-accuracy choices.
