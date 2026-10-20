# Feature-based roadmap for the discrete optimality and state-qualification tests

The purpose of this series is to test **what can go wrong in engineering state
finding**, using small problems whose behavior we can understand independently.
Each example should exercise a specific feature as well as the common C++20
numerical workflow. The oscillator's hidden velocity is the model for this
approach: an unmeasured state is recoverable through dynamics, but particular
measurement times can remove the information needed to recover it.

This is a test roadmap, not a claim that all proposed features are implemented.
Current coverage refers to the examples shipped with this repository. New examples
should estimate initial states under known inputs and known model parameters.
Identification of an arbitrary time-dependent forcing or inlet-composition
history belongs to the deferred time-dependent-control work.

## 1. Features to test

| Feature | Small engineering example | Question and required demonstration | Current coverage |
|---|---|---|---|
| **Hidden state** | Damped mechanical oscillator: observe displacement, infer initial velocity | Can the observations identify a state component that is never directly measured? Recover the full initial state and compare its sensitivities, errors, and engineering qualification with analytical references. | Implemented in [oscillator.md](oscillator.md). |
| **Poor measurement timing / sampling blindness** | Compare four oscillator displacement samples at half-periods and at well-spaced times | Show how half-period samples leave initial velocity unobservable, while proper sampling restores rank, stable recovery, and engineering qualification. Do not interpret numerical phase error as physical information. | Implemented as matched theta-discrete inverse problems with exact rank-one/rank-two references, adjoint checks, phase-leakage refinement study, bounded feasible-state sets, and persistent qualification in [oscillator.md](oscillator.md). |
| **Multiple decay timescales** | Motor winding and housing with variable thermal contact | Separate fast and slow modes. Show that the measured housing can look settled while uncertainty in the hidden winding state still prevents that goal from qualifying. Compare goal-specific qualification times with modal predictions. | Implemented in [weak_coupling.md](weak_coupling.md): three known contact conductances, housing-only sensing, exact modal reference, full discrete optimality validation, and bounded persistent qualification. |
| **Transient error amplification** | Two coupled storage states with one-way transfer, or a damped coupled mechanical system | Despite asymptotic stability, an initial perturbation can grow temporarily in a specified physical norm or engineering goal. Predict the largest gain and its initial direction, then verify the tangent, adjoint, uncertainty envelope, and delayed persistent qualification. | Implemented in [two_tank.md](two_tank.md), including worst-direction tangent/adjoint checks, analytical convergence, and delayed persistent qualification. |
| **Persistent unobservable mode** | Two insulated thermal masses; measure only their temperature difference | Recover or qualify the difference while a common temperature offset remains invisible and does not decay. With a prior domain broader than the absolute-temperature tolerance, correctly leave absolute-temperature goals unqualified. | Implemented in [insulated_thermal.md](insulated_thermal.md), including singular optimization, persistent feasible directions, and actual added absolute data. |
| **Nonlinear observation ambiguity** | Infer RLC voltage/current from resistor temperature driven by Joule heating $Ri^2$ | Exhibit distinct initial states with identical observations; local derivative rank must not be mistaken for global uniqueness. Preserve both feasible branches in qualification, then test how an added signed observation changes the result. | Implemented in [electrothermal.md](electrothermal.md): nonlinear dynamics, two full-rank local optima, global bounded box covers, invariant/signed goals, and informative versus insufficient signed measurements. |
| **Stiff dynamics** | Consecutive first-order reaction $A\to B\to C$, with $k_f/k_s=1000$ | Distinguish stability, damping, positivity, and accuracy. Show how an unresolved fast reaction mode affects concentrations, reconstruction, and qualification; compare midpoint, an intermediate theta, and backward Euler. Verify expected convergence after resolving both timescales. | Implemented in [stiff_reaction.md](stiff_reaction.md): exact solution, two product-observation schedules, discrete optimality checks, positivity diagnostic, bounded persistent qualification, and measured rates for eight workflow quantities. |

For every proposed feature, choose physical units and scaling before assessing
conditioning, perturbation size, amplification, or engineering tolerances.
Changing units must not manufacture an apparent physical difficulty.

## 2. What would count as a successful feature test?

### Hidden state

- Use observations of only part of the state.
- Identify the initial-state directions to which the observation record responds.
- Recover the unmeasured component without supplying its true initial value to
  the optimizer or adding an undisclosed observation.
- Compare the numerical optimum with the optimum of the analytical observation
  model for the **same retained noisy data**.
- Qualify goals involving the hidden state, including measurement and numerical
  uncertainty. A visually accurate reconstruction is insufficient by itself.

The oscillator already demonstrates this with displacement, velocity, and signed
support-force goals, including qualified and intentionally unqualified grids.

### Poor measurement timing

The number of observations alone does not determine their information content.
For the oscillator, the initial-velocity displacement sensitivity contains
$\sin(\omega(t-t_0))$, which vanishes at integer half-periods.

A complete follow-on test should:

- Verify the exact missing sensitivity direction and singular information matrix.
- Demonstrate that distinct initial states match the same poorly timed data.
- Avoid applying an ordinary inverse to a singular information matrix.
- Retain the unobservable direction in a bounded uncertainty set.
- Show which engineering goals remain unqualified and how moving or adding a
  measurement restores information.

Regularization may select one initial state, but must be reported as prior
information rather than evidence supplied by the measurements.

### Multiple decay timescales

Implemented in [the motor winding/housing example](weak_coupling.md). Vary a
known physical contact conductance to create a controlled separation between
fast and slow modes. Keep model parameters known during each state-finding
experiment.

Compare analytical modal decay, discrete modal behavior, and goal uncertainty.
The weak-contact case has a 423-minute slow-mode half-life. Housing-only data
leave the winding-temperature bound above tolerance through the 600-minute
horizon even though the housing and total stored-heat goals qualify at elapsed
zero. At moderate contact, the winding goal qualifies after 103.59 minutes. At
strong contact, the housing bound first passes and later fails before it
qualifies persistently. The study reports goal-specific decisions, modal rates,
the scaled information condition, exact references, and full observed rates.
It exercises the existing library without requiring a new generic API; it
also demonstrates that a full-rank Hessian can still correspond to weak
practical information.

### Transient amplification

Implemented in [the two-tank example](two_tank.md), following the oscillator. The concentration norm amplifies even though tracer inventory decreases. The refined bounded-data case initially meets the downstream tolerance, subsequently exceeds it, and qualifies persistently only from 2.515625 elapsed minutes. Shared library utilities now support tangent-map gains, polygon uncertainty sets, defect propagation, and whole-interval decisions.

For further amplification examples, specify the physical
state scaling and norm before computing a gain such as

$$
G(t)=\max_{\delta\boldsymbol\xi\ne0}
\frac{\|\delta\boldsymbol x(t)\|_W}
     {\|\delta\boldsymbol\xi\|_W},
\qquad t\in[t_0,t_f].
$$

For a linear example, use an exact matrix exponential as the independent
reference. Verify the maximum gain, its time, and the corresponding initial
perturbation direction. The tangent and adjoint should reproduce those results.

Separate growth of a particular downstream quantity from growth of a chosen
whole-state norm. Redistribution can increase one engineering quantity without
increasing total stored energy. Negative eigenvalue real parts alone do not rule
out transient amplification in a prescribed norm.

The qualification test must retain a possible temporary rise in uncertainty.
It must not infer persistent qualification from the first threshold crossing or
assume that asymptotic decay implies monotonically decreasing error.

### Persistent unobservable mode

Implemented in [the insulated thermal example](insulated_thermal.md). It uses the shared rank-aware quadratic optimizer and whole-interval qualification. Difference data alone qualify the difference on the refined midpoint grid while absolute temperatures remain unqualified. An actual added absolute datum restores rank and enables joint qualification.

For two insulated thermal masses with internal heat exchange, adding the same
constant to both initial temperatures leaves their temperature difference
unchanged. That common offset persists because there is no heat exchange with an
external reference temperature.

Use a finite, explicit initial-state domain. Demonstrate that the data constrain
the difference but do not narrow the common-offset direction. With appropriate
goal tolerances, the correct result is a qualified difference goal and unqualified
absolute-temperature goals. Passing this test means reporting the missing
information correctly, not forcing every goal to qualify.

### Nonlinear observation ambiguity

Implemented in [the electrothermal discharge example](electrothermal.md).
Joule heating supplies the physical nonlinear observation mechanism: opposite
initial capacitor voltages and inductor currents produce identical resistor
temperatures. Each branch has full-rank local curvature, yet both must remain
in qualification. Shared box-subdivision and union range routines search the
stated prior rather than treating multistart optimization as a completeness proof.

Temperature and stored-energy goals qualify with both polarities retained;
signed electrical goals remain unqualified on the stated horizon. An actual
accurate signed datum removes the opposite region and enables joint
qualification. A weak signed datum changes local cost preference while leaving
both regions admissible. The original squared-displacement oscillator remains
a possible alternative illustration, not an implemented separate example.

### Stiff dynamics

Implemented in [the stiff consecutive-reaction example](stiff_reaction.md).
The exact $A\to B\to C$ mechanism uses $k_f/k_s=1000$ and observes only product
$C$. It tests the forward model, discrete optimality derivatives, reconstruction,
and bounded qualification. On a step ten times larger than the fast timescale,
midpoint is stable but leaves an alternating fast mode and produces negative
species; backward Euler damps the same mode much more strongly. Refinement
recovers the measured second-order rate for $\theta=0.5$ and first-order rates
for $\theta=0.65$ throughout the state-finding workflow.

The case also makes clear that numerical damping does not create measurement
information. The prior and product data are retained in a bounded feasible set,
and the species goals qualify only after the hidden fast-composition uncertainty
has decayed below tolerance.

## 3. Shared test contract for every dynamic example

Each new feature example should include the following, with exceptions stated
explicitly rather than silently omitted.

1. **Standalone problem statement.** State what is unknown, what is measured,
   why the feature matters, all units, known inputs, model assumptions, and
   $t\in[t_0,t_f]$ with $\boldsymbol x(t_0)=\boldsymbol\xi$.
2. **Independent analytical reference.** Prefer exact states and sensitivities.
   Derive objective, adjoint, gradient, Hessian, and reference-optimum formulas
   wherever tractable. Clearly identify any reference that instead uses separate
   high-accuracy numerical integration.
3. **Actual theta forward solve.** Use the library's numerical solver. Analytical
   propagation may generate observations and validation references, but must not
   replace the forward solve used by optimization.
4. **Matching discrete objective.** Use theta-consistent quadrature for running
   costs. Isolated observations use nodal costs without timestep weights. Retain
   the same physical observations and error scales across grids.
5. **Tangent and discrete adjoint.** Verify directional state sensitivities,
   tangent–adjoint pairing, and the gradient of the implemented discrete cost.
   Explain any costate convention needed for a continuous-reference comparison.
   With no terminal source, use $\boldsymbol\lambda(t_f)=0$.
6. **Incremental adjoint and optimization.** Check Hessian-vector products and
   stationarity. Distinguish convergence of the optimizer from temporal
   convergence and from identifiability. Test nearby starts or multiple branches
   where they are relevant.
7. **Physical balance or invariant.** Include an independent check of conservation,
   dissipation, or another model-specific identity when available.
8. **Temporal refinement.** Compare state, tangent, adjoint, objective derivatives,
   and reconstructed controls. Test second order for theta = 0.5 and the expected
   order for comparison methods once the smooth problem is resolved. Report
   observed rates, not just an expected reference slope on a plot.
9. **Engineering qualification.** Specify goals, tolerances, uncertainty
   assumptions, and prior domains. Include numerical model error, check the
   between-node behavior, and report individual and joint persistent times or
   explicit unqualified decisions. Keep local covariance indicators distinct
   from the final uncertainty-set argument.
10. **C++20 decisions and reproducible outputs.** All forward, inverse, reference,
    error-budget, and qualification logic runs in C++20. Python renders retained
    outputs. Supply machine-readable results, labeled plots, and a tutorial that
    explains why each step is performed.

A separate exactly integrable running-cost verification functional is acceptable
when the physical inverse problem contains only point measurements, as in the
oscillator. It must not be presented as additional sensor information used to
recover a hidden state.

Qualification is finite-horizon and persistent: a goal must satisfy its tolerance
at every remaining time through $t_f$. When all observations are used to qualify
an earlier trajectory, label the decision **retrospective**, not online.

## 4. Cross-cutting features to extend across the series

These complement the seven physical/observability features above. Their dedicated
coverage should be tracked independently of whether an API already accepts the
relevant input.

| Feature | Required test |
|---|---|
| Variable timesteps and nonuniform schedules | Keep observation times fixed; verify the forward/cost/tangent/adjoint consistency on nonuniform meshes and define the refinement family explicitly. Implemented in [nonuniform_mesh.md](nonuniform_mesh.md), with a nested graded oscillator mesh and observed orders for eleven optimality-system quantities. |
| Goal dependence | On the same reconstructed state and uncertainty set, show goals that qualify at different times and cases where one blocked goal prevents joint qualification. Implemented as a dedicated two-case oscillator cross-cutting test in [goal_dependence.md](goal_dependence.md). |
| Measurement quality versus numerical accuracy | Refine the mesh with data fixed, then separately change measurement bounds. Demonstrate which budget term limits qualification. |
| Observation weighting and correlated errors | Specify measurement covariance independently of goal-adjoint sensitivity. A dedicated correlated-noise example needs an appropriate joint uncertainty construction; the independent-Gaussian noise-box formula must not be reused unchanged. |
| Units and scaling | Verify consistent results under unit conversion and use physically stated scales for conditioning and gain comparisons. |
| State constraints and boundary optima | Test admissibility and positivity where relevant. Boundary KKT optimization is an explicit future solver extension, not something the current interior Newton driver already implements. |
| Weak information and near singularity | Vary measurement design continuously toward a blind direction; distinguish a small numerical gradient from an accurately identified initial state. |
| Nonlinear dynamics | Retain exact or independently verified references; test state-dependent Jacobians and the second-order terms required by the incremental adjoint. Existing nonlinear decay supplies a foundation. |
| Added measurement experiments | Add an actual observation, reoptimize, and recompute the uncertainty set and qualification. Do not merely shrink covariance by assumption. |

## 5. Suggested sequence and completion checklist

The existing examples form the baseline: nonlinear decay, two-temperature state
finding, static Peng–Robinson composition inference, and the hidden-state
oscillator. The static EOS example does not count as a transient theta test.

A proposed continuation is:

1. **Completed:** transient amplification with an exact linear reference; see [two_tank.md](two_tank.md).
2. **Completed:** persistent unobservable thermal mode with qualified and unqualified goals; see [insulated_thermal.md](insulated_thermal.md).
3. **Completed:** nonlinear observation ambiguity and branch-resolving measurements in [electrothermal.md](electrothermal.md).
4. **Completed:** a practical weak-coupling timescale sweep in the motor winding/housing example; see [weak_coupling.md](weak_coupling.md).
5. **Completed:** a stiff, analytically solvable chemical reaction with bounded
   species qualification; see [stiff_reaction.md](stiff_reaction.md).
6. **Completed:** matched poorly and properly timed oscillator reconstruction and qualification experiments; see [oscillator.md](oscillator.md#15-poorly-timed-sampling-is-a-full-inverse-problem-test).

The ordering can change as the series develops. For each new example, complete
and record:

- [ ] The feature and the failure it should expose.
- [ ] The engineering model, physical scaling, and independent references.
- [ ] The observation design and identifiability argument.
- [ ] Forward, objective, tangent, adjoint, incremental-adjoint, and optimization checks.
- [ ] Temporal-convergence results at fixed controls and for recovered initial states.
- [ ] Engineering goals, measurement assumptions, numerical allowances, and persistent decisions.
- [ ] At least one intentionally difficult or unqualified case with an explained result.
- [ ] C++20 outputs, plots, build/test integration, and standalone documentation.

A test succeeds when it reaches the correct conclusion—including a justified
unqualified or nonidentifiable result. The purpose is to determine which
engineering predictions the available information supports, and when it supports
them, rather than simply to produce an initial state that fits the observations.
