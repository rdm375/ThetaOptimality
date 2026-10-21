# ThetaOptimality

**A small C++20 library for building and testing a numerically consistent
discrete optimality system for time-dependent inverse problems.**

The project develops one coupled discrete construction: a theta-method state
trajectory, a matching theta-quadrature objective, its tangent linear model,
the exact discrete adjoint, reduced gradients, and incremental-adjoint
Hessian-vector products. It also includes state-estimation and engineering
qualification experiments that compare discrete results with analytical
solutions, bounded observation sets, and numerical error estimates.

The purpose is educational and methodological: make the pieces inspectable on
small problems, test their derivatives and convergence, and expose cases where
good numerical optimization does not imply that the physical state is
identifiable or sufficiently known. This is a research and teaching codebase,
not a general-purpose sparse PDE or production estimation package.

## Why this project

For a time-discrete model, the objective and its derivatives should belong to
the same discrete problem as the forward state. Differentiating a different
continuous or discretized problem can produce gradients that do not match the
numerical objective actually optimized. ThetaOptimality makes that connection
explicit and tests it with directional finite differences, tangent-adjoint
pairings, Hessian checks, analytical references, and mesh refinement.

The examples also separate three questions that are easy to conflate:

1. **Did the numerical optimization converge?**
2. **Do the observations determine the state or parameters?**
3. **Does the available information establish the requested engineering
   tolerance, including measurement and numerical uncertainty?**

An unqualified result means the stated assumptions and error budgets do not
establish the tolerance. It does not prove that the actual state error exceeds
it. Local information diagnostics are not global uniqueness certificates.

## The discrete workflow

For a model \(\dot{x}=f(t,x,\xi)\) with unknown finite-dimensional vector
\(\xi\), the library supplies these parts of one discrete optimality system:

| Component | Role in the experiment |
|---|---|
| Theta-method forward solve | Defines the numerical trajectory used in the inverse problem. |
| Matching theta-quadrature objective | Evaluates running and terminal costs on that trajectory. |
| Tangent linear model | Propagates perturbations in initial conditions or model controls. |
| Discrete adjoint | Computes objective sensitivities through the discrete equations. |
| Reduced gradient and Hessian-vector products | Drive optimization and test second-order derivatives. |
| State qualification | Compares whole-interval state uncertainty and numerical error with engineering tolerances. |

The complete equations and conventions are in
[`docs/formulation.md`](docs/formulation.md); callback contracts and usage are
in [`docs/api.md`](docs/api.md). The library is header-only in
[`include/theta_opt/`](include/theta_opt/).

## Experiments and documentation

Each experiment has a standalone write-up describing its model, assumptions,
numerical formulation, checks, findings, and limitations. The table links to
every experiment document so readers can go directly from the overview to the
derivation.

| Experiment | What it tests | Documentation |
|---|---|---|
| Nonlinear decay and accumulated exposure | Full forward, cost, tangent, adjoint, gradient, Hessian, and analytical convergence workflow. | [Decay example](docs/decay.md) · [State qualification](docs/state_qualification.md) |
| Electrical state from resistor temperature | Nonlinear observation ambiguity, local curvature, and disconnected feasible state sets. | [Electrothermal experiment](docs/electrothermal.md) |
| Persistent unobservable thermal mode | Difference observations identify a decaying mode but not the common temperature offset. | [Insulated thermal experiment](docs/insulated_thermal.md) |
| Poorly timed oscillator samples | Half-period sampling loses the initial-velocity direction; moving sample times restores it. | [Oscillator experiment](docs/oscillator.md) · [Weak-information sweep](docs/weak_information.md) |
| Nonuniform oscillator mesh | Tests the full discrete system and convergence on nested, nonuniform time grids. | [Nonuniform-mesh experiment](docs/nonuniform_mesh.md) |
| Goal-dependent qualification | Shows why different physical goals can qualify at different times and one can block a joint decision. | [Goal-dependence experiment](docs/goal_dependence.md) |
| Transient amplification in two tanks | A stable coupled system can temporarily amplify concentration error and delay qualification. | [Two-tank experiment](docs/two_tank.md) |
| Hidden-state oscillator | Reconstructs an unmeasured initial velocity from displacement observations. | [Oscillator experiment](docs/oscillator.md) |
| Static Peng–Robinson composition | Estimates gas-mixture fractions from thermodynamic density observations and compares identifiability. | [Peng–Robinson experiment](docs/peng_robinson.md) |
| Two-temperature state finding | Infers two coupled temperatures from one thermometer and tests engineering qualification. | [Two-temperature experiment](docs/two_temperature.md) |
| Time-varying heater schedule | Estimates a finite-dimensional control schedule from thermal observations. | [Thermal schedule experiment](docs/time_varying_thermal.md) · [Control and path-constraint theory](docs/time_varying_controls_path_constraints.md) |
| Weakly coupled winding and housing | Tests separated thermal time scales, hidden modes, and goal-specific qualification. | [Weak-coupling experiment](docs/weak_coupling.md) |
| Stiff sequential reaction | Uses a closed-form chemical system to examine stiffness, positivity, sampling, and convergence. | [Stiff-reaction experiment](docs/stiff_reaction.md) |
| State constraints and boundary optima | Checks active bounds, KKT multipliers, stationarity, and qualification at constrained optima. | [Bounded washout experiment](docs/bounded_washout.md) |
| Nonlinear logistic growth | Estimates an initial state for nonlinear dynamics and checks discrete derivative convergence. | [Logistic-growth experiment](docs/logistic_growth.md) |
| Joint predator–prey recovery | Separately estimates initial conditions, model parameters, and both together. | [Predator–prey experiment](docs/predator_prey.md) |
| Measurement versus numerical accuracy | Separates sensor-error and discretization contributions in a drug-washout model. | [Drug-washout experiment](docs/drug_washout.md) |
| Correlated observation errors | Tests the effect of correlated assay errors on fit, information, and qualification. | [Correlated-washout experiment](docs/drug_washout_correlated.md) |
| Unit and scaling invariance | Checks that changing concentration and time units transforms derivatives and decisions consistently. | [Unit-scaling experiment](docs/unit_scaling.md) |
| Time-varying gas composition | Reconstructs a changing feed composition and vessel state from Peng–Robinson density observations. | [Dynamic gas experiment](docs/time_varying_gas.md) |
| Roadmap and cross-cutting tests | Tracks the theory features tested across multiple examples and future cases. | [Feature-test roadmap](docs/feature_test_roadmap.md) |

These cases deliberately include both qualified and **NOT QUALIFIED** results.
An unqualified result means that the stated data and error assumptions do not
establish the requested tolerance; it does not prove that the actual state
error exceeds that tolerance. Local information or a small optimization
gradient is not, by itself, evidence of global uniqueness or engineering
sufficiency.

## Start here

The nonlinear decay example is the smallest complete optimality-system
verification. Its continuous solution, objective, tangent, adjoint, gradient,
and Hessian have analytical references. Numerical optimization still uses the
library's discrete forward solve and objective; analytical formulas are used
as independent comparisons, not as substitutes for the numerical calculation.

```bash
make example       # Run the basic decay optimality-system example
make validation    # Check library identities and derivatives
make test          # Build and run the C++ validation and experiment suite
```

To run the full set of examples and regenerate their figures:

```bash
make qualification-plots
```

Requirements: a C++20 compiler, GNU Make, and the C++ standard library. Python
3 with Matplotlib is only needed to render figures; simulations, fits,
derivatives, convergence estimates, and qualification decisions are computed
in C++20. The optional Peng–Robinson examples use the pinned scalar CPU
ThermoGPU source included in `third_party/ThermoGPU`; they do not require CUDA.
See its [`PROVENANCE.md`](third_party/ThermoGPU/PROVENANCE.md) for provenance
and licensing details.

Install the plotting dependency, if needed, with:

```bash
python3 -m pip install matplotlib
```

The Peng–Robinson examples use the pinned scalar CPU ThermoGPU implementation
in [`third_party/ThermoGPU/`](third_party/ThermoGPU/). No CUDA is required.
Its provenance and separate licensing are documented in
[`third_party/ThermoGPU/PROVENANCE.md`](third_party/ThermoGPU/PROVENANCE.md).

CMake is also supported:

```bash
cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build-cmake --parallel
ctest --test-dir build-cmake --output-on-failure
```

## Example gallery

Each example documents its model, assumptions, theory exercised, checks, and
limitations in a standalone page under [`docs/`](docs/). Figures below are
committed outputs from the corresponding C++ experiment; their Python scripts
only render the saved results.

| Case | What it tests | Experiment notes and representative figure |
|---|---|---|
| **Nonlinear decay and exposure** | Complete forward, quadrature, tangent, adjoint, gradient, Hessian, and second-order midpoint convergence workflow. | [State qualification](docs/state_qualification.md) · [qualification figure](results/qualification/qualification_verdict.png) |
| **Hidden-state oscillator** | Recover an unmeasured initial velocity; compare analytical and discrete fits, derivatives, and convergence. | [Oscillator](docs/oscillator.md) · [reconstruction](results/oscillator/reconstruction.png) |
| **Poor versus proper sampling** | Show exact rank loss at half-period sample times and recovery when the same number of samples are placed informatively. | [Sampling experiment](docs/oscillator.md#15-poorly-timed-sampling-is-a-full-inverse-problem-test) · [comparison](results/oscillator_poor_timing/sampling_comparison.png) |
| **State qualification** | Retain bounded data-consistent uncertainty and include a numerical-model allowance in whole-interval decisions. | [Qualification theory](docs/state_qualification.md) · [decision plot](results/qualification/qualification_verdict.png) |
| **Time-varying thermal control** | Estimate a finite-dimensional heater schedule using a discrete optimality system. | [Thermal schedule](docs/time_varying_thermal.md) · [reconstruction](results/time_varying_thermal/reconstruction.png) |
| **Dynamic gas composition** | Infer initial vessel composition and a changing feed from Peng–Robinson density observations. | [Gas schedule](docs/time_varying_gas.md) · [reconstruction](results/time_varying_gas/composition_reconstruction.png) |
| **Stiff reaction kinetics** | Test a closed-form sequential reaction, stiff time scales, positivity, observation timing, and qualification. | [Stiff reaction](docs/stiff_reaction.md) · [reconstruction](results/stiff_reaction/stiffness_and_reconstruction.png) |
| **Nonlinear population dynamics** | Estimate initial state, model parameters, and both jointly for Lotka–Volterra dynamics. | [Predator–prey](docs/predator_prey.md) · [phase portrait](results/predator_prey/phase_portrait.png) |
| **Ambiguity and missing modes** | Preserve multiple feasible branches or persistent unobservable directions rather than reporting a unique-looking fit. | [Electrothermal](docs/electrothermal.md) · [Insulated thermal](docs/insulated_thermal.md) |
| **Cross-cutting cases** | Nonuniform meshes, transient amplification, correlated observations, units/scaling, path constraints, and goal-dependent qualification. | [Feature roadmap](docs/feature_test_roadmap.md) |

The complete set of derivations and experiments is indexed in
[`docs/`](docs/feature_test_roadmap.md). To regenerate figures after running
the relevant C++ cases, use `make qualification-plots` (requires Matplotlib).

## Library use

For a model with \(N\) states and \(P\) unknown controls, define a
`Problem<N,P>` and a time `Grid`:

```cpp
#include <theta_opt/theta_opt.hpp>
#include "decay_model.hpp"

using namespace theta_opt;

System<2, 1> system(decay::problem(),
    Grid::uniform(decay::t0, decay::tf, 160, 0.5));
Vec<1> xi{2.4};

auto trajectory = system.forward(xi);
auto tangent = system.tangent(xi, trajectory, Vec<1>{1.0});
auto adjoint = system.adjoint(xi, trajectory);
auto gradient = adjoint.gradient;
auto hessian_vector = system.hessian_vector(xi, Vec<1>{1.0});
auto evaluation = system.evaluate(xi); // forward + cost + adjoint + gradient
```

Second-order callbacks are explicit. A Hessian request without the required
directional derivatives throws; there is no silent finite-difference Hessian.
See [`docs/api.md`](docs/api.md) for the callback contract and
[`docs/formulation.md`](docs/formulation.md) for the full discrete equations.

## Scope and interpretation

- The core uses dense, compile-time-sized vectors and matrices for small ODE
  systems; it is not a sparse PDE or large-scale optimization library.
- Time-varying controls are finite-dimensional schedules or bases, not
  unrestricted functions of time.
- Qualification conclusions depend on the stated model, prior, observation
  assumptions, tolerances, and numerical-error estimates. They are not formal
  directed-rounding certificates or safety certifications.
- The examples are synthetic or idealized method demonstrations, not
  validated clinical, industrial, or safety-critical tools.
- See the [feature-test roadmap](docs/feature_test_roadmap.md) for current
  coverage, limitations, and proposed cases.

## License

ThetaOptimality project code is licensed under
[`AGPL-3.0-only`](LICENSE). ThermoGPU source and the connected composition
examples carry their separately identified GPL-3.0-or-later terms; see
[`third_party/ThermoGPU/PROVENANCE.md`](third_party/ThermoGPU/PROVENANCE.md).
The AGPL version is available under the terms of `LICENSE`. Separate
commercial licensing, if offered by the rights holders, requires a separate
written agreement. See [`COPYRIGHT`](COPYRIGHT) for the licensing notice.
