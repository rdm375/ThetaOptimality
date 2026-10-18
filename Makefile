CXX ?= c++
CPPFLAGS ?= -Iinclude
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror
HEADERS := $(wildcard include/theta_opt/*.hpp) $(wildcard examples/*_model.hpp)

.PHONY: all test validation example plots qualification qualification-plots decay-plots thermal thermal-plots clean
all: build/theta_decay build/theta_validation build/decay_qualification build/two_temperature_qualification build/peng_robinson_inverse build/oscillator

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

test: validation qualification thermal gas oscillator

validation: build/theta_validation
	./build/theta_validation

example: build/theta_decay
	./build/theta_decay results

plots: example
	python3 tools/plot_results.py

qualification: build/decay_qualification
	./build/decay_qualification results/qualification

qualification-plots: decay-plots thermal-plots gas-plots oscillator-plots

decay-plots: qualification
	python3 tools/plot_qualification.py

thermal: build/two_temperature_qualification
	./build/two_temperature_qualification results/two_temperature

thermal-plots: thermal
	python3 tools/plot_two_temperature.py

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

.PHONY: oscillator oscillator-plots
oscillator: build/oscillator
	./build/oscillator results/oscillator

oscillator-plots: oscillator
	python3 tools/plot_oscillator.py
