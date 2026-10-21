CXX ?= c++

CPPFLAGS ?= -Iinclude
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror
HEADERS := $(wildcard include/theta_opt/*.hpp) $(wildcard examples/*_model.hpp)

.PHONY: all test validation example plots qualification qualification-plots decay-plots thermal thermal-plots clean logistic-growth logistic-growth-plots
all: build/theta_decay build/theta_validation build/decay_qualification build/two_temperature_qualification build/weak_coupling build/peng_robinson_inverse build/oscillator build/oscillator_poor_timing build/oscillator_information_sweep build/bounded_washout build/nonuniform_oscillator build/goal_dependence build/two_tank build/insulated_thermal build/electrothermal build/stiff_reaction build/drug_washout build/drug_washout_correlated build/unit_scaling build/logistic_growth

build:
	mkdir -p build

build/theta_decay: examples/decay.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/theta_validation: tests/validation.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/decay_qualification: examples/decay_qualification.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/two_temperature_qualification: examples/two_temperature_qualification.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/weak_coupling: examples/weak_coupling.cpp examples/weak_coupling_model.hpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

test: validation qualification thermal weak-coupling gas oscillator poor-oscillator proper-oscillator information-sweep bounded-washout nonuniform-oscillator goal-dependence two-tank insulated-thermal electrothermal stiff-reaction drug-washout correlated-washout unit-scaling logistic-growth predator-prey

validation: build/theta_validation
	./build/theta_validation

example: build/theta_decay
	./build/theta_decay results

plots: example
	python3 tools/plot_results.py

qualification: build/decay_qualification
	./build/decay_qualification results/qualification

qualification-plots: decay-plots thermal-plots weak-coupling-plots gas-plots oscillator-plots information-sweep-plots bounded-washout-plots nonuniform-oscillator-plots goal-dependence-plots two-tank-plots insulated-thermal-plots electrothermal-plots stiff-reaction-plots drug-washout-plots correlated-washout-plots unit-scaling-plots logistic-growth-plots predator-prey-plots

decay-plots: qualification
	python3 tools/plot_qualification.py

thermal: build/two_temperature_qualification
	./build/two_temperature_qualification results/two_temperature

thermal-plots: thermal
	python3 tools/plot_two_temperature.py

.PHONY: weak-coupling weak-coupling-plots
weak-coupling: build/weak_coupling
	./build/weak_coupling results/weak_coupling

weak-coupling-plots: weak-coupling
	python3 tools/plot_weak_coupling.py

build/stiff_reaction: examples/stiff_reaction.cpp examples/stiff_reaction_model.hpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/drug_washout: examples/drug_washout.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: drug-washout drug-washout-plots
drug-washout: build/drug_washout
	./build/drug_washout results/drug_washout

drug-washout-plots: drug-washout
	python3 tools/plot_drug_washout.py

build/drug_washout_correlated: examples/drug_washout_correlated.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: correlated-washout correlated-washout-plots
correlated-washout: build/drug_washout_correlated
	./build/drug_washout_correlated results/drug_washout_correlated

correlated-washout-plots: correlated-washout
	python3 tools/plot_drug_washout_correlated.py

build/unit_scaling: examples/unit_scaling.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: unit-scaling unit-scaling-plots
unit-scaling: build/unit_scaling
	./build/unit_scaling results/unit_scaling

unit-scaling-plots: unit-scaling
	python3 tools/plot_unit_scaling.py

.PHONY: stiff-reaction stiff-reaction-plots
stiff-reaction: build/stiff_reaction
	./build/stiff_reaction results/stiff_reaction

stiff-reaction-plots: stiff-reaction
	python3 tools/plot_stiff_reaction.py

clean:
	rm -rf build

# Pinned unmodified ThermoGPU scalar CPU backend (no CUDA required).
THERMOGPU_ROOT ?= third_party/ThermoGPU
build/peng_robinson_inverse: examples/peng_robinson_inverse.cpp examples/peng_robinson_model.hpp examples/peng_robinson_bounds.hpp $(HEADERS) $(THERMOGPU_ROOT)/src/cpu/pr_scalar.cpp $(wildcard $(THERMOGPU_ROOT)/include/thermogpu/*.hpp) | build
	$(CXX) $(CPPFLAGS) -I$(THERMOGPU_ROOT)/include $(CXXFLAGS) $< $(THERMOGPU_ROOT)/src/cpu/pr_scalar.cpp -o $@

.PHONY: gas gas-plots
gas: build/peng_robinson_inverse
	./build/peng_robinson_inverse results/peng_robinson

gas-plots: gas
	python3 tools/plot_peng_robinson.py

build/oscillator: examples/oscillator.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/oscillator_poor_timing: examples/oscillator_poor_timing.cpp examples/oscillator_model.hpp examples/oscillator_qualification_model.hpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/oscillator_information_sweep: examples/oscillator_information_sweep.cpp examples/oscillator_model.hpp examples/oscillator_qualification_model.hpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/bounded_washout: examples/bounded_washout.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/logistic_growth: examples/logistic_growth.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

logistic-growth: build/logistic_growth
	./build/logistic_growth results/logistic_growth

logistic-growth-plots: logistic-growth
	python3 tools/plot_logistic_growth.py

build/nonuniform_oscillator: examples/nonuniform_oscillator.cpp examples/oscillator_model.hpp examples/oscillator_qualification_model.hpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/predator_prey: examples/predator_prey.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

predator-prey: build/predator_prey
	./build/predator_prey results/predator_prey

predator-prey-plots: predator-prey
	python3 tools/plot_predator_prey.py

build/goal_dependence: examples/goal_dependence.cpp examples/oscillator_model.hpp examples/oscillator_qualification_model.hpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: oscillator oscillator-plots
oscillator: build/oscillator
	./build/oscillator results/oscillator

oscillator-plots: oscillator poor-oscillator proper-oscillator
	python3 tools/plot_oscillator.py

.PHONY: poor-oscillator proper-oscillator poor-oscillator-plots sampling-comparison
poor-oscillator: build/oscillator_poor_timing
	./build/oscillator_poor_timing results/oscillator_poor_timing

proper-oscillator: build/oscillator_poor_timing
	./build/oscillator_poor_timing results/oscillator_proper_timing proper

sampling-comparison: poor-oscillator proper-oscillator
	python3 tools/plot_oscillator.py

poor-oscillator-plots: sampling-comparison

.PHONY: information-sweep information-sweep-plots
information-sweep: build/oscillator_information_sweep
	./build/oscillator_information_sweep results/oscillator_information_sweep

information-sweep-plots: information-sweep
	python3 tools/plot_oscillator_information_sweep.py

.PHONY: bounded-washout bounded-washout-plots
bounded-washout: build/bounded_washout
	./build/bounded_washout results/bounded_washout

bounded-washout-plots: bounded-washout
	python3 tools/plot_bounded_washout.py

.PHONY: nonuniform-oscillator nonuniform-oscillator-plots
nonuniform-oscillator: build/nonuniform_oscillator
	./build/nonuniform_oscillator results/nonuniform_oscillator

nonuniform-oscillator-plots: nonuniform-oscillator
	python3 tools/plot_nonuniform_oscillator.py

.PHONY: goal-dependence goal-dependence-plots
goal-dependence: build/goal_dependence
	./build/goal_dependence results/goal_dependence

goal-dependence-plots: goal-dependence
	python3 tools/plot_goal_dependence.py

build/two_tank: examples/two_tank.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: two-tank two-tank-plots
two-tank: build/two_tank
	./build/two_tank results/two_tank

two-tank-plots: two-tank
	python3 tools/plot_two_tank.py

build/insulated_thermal: examples/insulated_thermal.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: insulated-thermal insulated-thermal-plots
insulated-thermal: build/insulated_thermal
	./build/insulated_thermal results/insulated_thermal

insulated-thermal-plots: insulated-thermal
	python3 tools/plot_insulated_thermal.py

build/electrothermal: examples/electrothermal.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

.PHONY: electrothermal electrothermal-plots
electrothermal: build/electrothermal
	./build/electrothermal results/electrothermal

electrothermal-plots: electrothermal
	python3 tools/plot_electrothermal.py
