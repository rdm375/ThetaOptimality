# ThermoGPU scalar backend

Source: https://github.com/rdm375/ThermoGPU

Pinned revision: `5a1a4174231a389b1724571eb646d48bfddb98e7`

The six `include/thermogpu/*.hpp` files, `src/cpu/pr_scalar.cpp`, the two data CSVs,
and `LICENSE` are unmodified copies from that revision. Only the scalar CPU
backend is included. There is no CUDA requirement. The full upstream repository
contains the other backends, tests, documentation, and benchmarks.

Copyright (C) 2026 Richard Myers. SPDX-License-Identifier: GPL-3.0-or-later.
The composition example and its adapters also use GPL-3.0-or-later. The optional
linked example is distinct from the header-only theta library; users distributing
that executable must observe ThermoGPU's included license.

The example's three component constants are copied verbatim from the included
components.csv. All nine binary interaction entries are explicitly set to zero:
this is an experiment assumption, not an empirical interaction dataset. The
upstream CSV supplies diagonal entries only. The upstream backend uses
Omega_a = 0.4572355289213822 and Omega_b = 0.0777960739038884.
