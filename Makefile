CXX ?= c++
CPPFLAGS ?= -Iinclude
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror
HEADERS := $(wildcard include/theta_opt/*.hpp) examples/decay_model.hpp

.PHONY: all test example plots clean
all: build/theta_decay build/theta_validation

build:
	mkdir -p build

build/theta_decay: examples/decay.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/theta_validation: tests/validation.cpp $(HEADERS) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

test: build/theta_validation
	./build/theta_validation

example: build/theta_decay
	./build/theta_decay results

plots: example
	python3 tools/plot_results.py

clean:
	rm -rf build
