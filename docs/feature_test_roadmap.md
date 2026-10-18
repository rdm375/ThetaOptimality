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
| **Poor measurement timing / sampling blindness** | Observe the oscillator only at integer half-periods | Show a nonzero initial-state direction that leaves every observation unchanged. Contrast this with a well-spaced schedule. Do not interpret numerical phase error as genuine physical observability. | Analytical timing diagnostic implemented in the oscillator; a full ill-posed reconstruction and qualification case remains to be added. |
| **Multiple decay timescales** | Two thermal masses with weak coupling | Separate fast and slow modes. Show that a state can look settled while uncertainty in a slow mode still prevents a selected goal from qualifying. Compare goal-specific qualification times with modal predictions. | Two-temperature model provides an analytical two-mode foundation; a dedicated weak-coupling / scale-separation sweep remains to be added. |
| **Transient error amplification** | Two coupled storage states with one-way transfer, or a damped coupled mechanical system | Despite asymptotic stability, an initial perturbation can grow temporarily in a specified physical norm or engineering goal. Predict the largest gain and its initial direction, then verify the tangent, adjoint, uncertainty envelope, and delayed persistent qualification. | Proposed dedicated example. |
| **Persistent unobservable mode** | Two insulated thermal masses; measure only their temperature difference | Recover or qualify the difference while a common temperature offset remains invisible and does not decay. With a prior domain broader than the absolute-temperature tolerance, correctly leave absolute-temperature goals unqualified. | Proposed dedicated example. |
| **Nonlinear observation ambiguity** | Observe squared oscillator displacement, $y=q^2$, without a sign-sensitive measurement | Exhibit distinct initial states with identical observations; local derivative rank must not be mistaken for global uniqueness. Preserve both feasible branches in qualification, then test how an added signed observation changes the result. | Proposed dedicated example. Static PR composition inference provides related nonlinear inverse-problem experience, but does not establish this sign-ambiguity test. |
| **Stiff dynamics** | Coupled fast and slow thermal relaxation | Distinguish stability from accuracy. Show how unresolved fast modes affect state reconstruction and qualification; test midpoint's weak damping of very stiff modes as well as theta values with stronger damping. Recover expected convergence once the relevant scales are resolved. | Proposed dedicated example. Existing theta comparisons are not a dedicated stiffness test. |

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

Vary a known physical coupling to create a controlled separation between fast and
slow modes. Keep model parameters known during each state-finding experiment.

Compare analytical modal decay, discrete modal behavior, and goal uncertainty.
A temperature difference may qualify much earlier than a goal dominated by the
slow mean-temperature mode. Include measurement schedules that poorly observe
the slow direction, and an added measurement that targets it.

### Transient amplification

This is a proposed next example, following the oscillator. Specify the physical
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

For a homogeneous oscillator and observations $q(t)^2$, the initial states
$\boldsymbol\xi$ and $-\boldsymbol\xi$ generate identical observation records.
Their local sensitivity matrices can each have full rank away from degeneracies;
that does not remove the global ambiguity.

The test should retain both branches over the stated initial-state domain.
Compare signed displacement/velocity goals with sign-invariant goals. A small
local covariance around one optimum must not be used as a global uncertainty
certificate. Add a sign-sensitive measurement and determine whether it actually
separates the branches at the chosen measurement time and error level.

### Stiff dynamics

Construct an exact linear reference with well-separated decay rates. Include
initial conditions that excite the fast mode, rather than choosing an initial
state that hides it from the experiment.

For theta = 0.5, distinguish second-order convergence at resolved timesteps from
poor damping of an unresolved fast mode. Compare against theta = 1 and an
intermediate value such as 0.65 where useful. Stability, physical admissibility,
observation fit, and engineering qualification are separate checks.

Temporal-order claims should use a refinement range that reaches the asymptotic
regime. A stable coarse calculation need not be accurate enough to reconstruct
or qualify the state.

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
| Variable timesteps and nonuniform schedules | Keep observation times fixed; verify the forward/cost/tangent/adjoint consistency on nonuniform meshes and define the refinement family explicitly. |
| Goal dependence | On the same reconstructed state and uncertainty set, show goals that qualify at different times and cases where one blocked goal prevents joint qualification. |
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

1. Transient amplification with an exact linear reference.
2. A persistent unobservable thermal mode with qualified and unqualified goals.
3. Nonlinear observation ambiguity and branch-resolving measurements.
4. A weak-coupling timescale sweep, followed by a deliberately stiff case.
5. A complete poorly timed oscillator reconstruction/qualification experiment.

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
