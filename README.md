# HAMMERHEAD

[![CI](https://github.com/jmurphy6895/HAMMERHEAD.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jmurphy6895/HAMMERHEAD.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/jmurphy6895/HAMMERHEAD.jl/branch/main/graph/badge.svg?token=47G4OLV6PD)](https://codecov.io/gh/jmurphy6895/HAMMERHEAD.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

## Description

This project (**H**igh-**A**ccuracy **M**odelling **E**nvironment **R**esourcing **H**ighly-**E**fficient Solvers **A**nd **D**ifferentiability) implements high-fidelity Equations of Motion using Zonal Harmonics, Solar Radiation Pressure, Drag, and Third-Body Perturbations based on packages developed in the Satellite Toolbox.

These forces are then called from a number of high performance integrator methodologies. Some of these are ported and implemented from the THALASSA library [CITE]. A full list can be found below:

1. Cowell
2. EDromo
3. Kustaanheimo-Stiefel
4. Stiefel-Scheifel
5. Unified State Model
6. Gauss Variational Equations

## Installation

```julia
julia> using Pkg
julia> Pkg.add("HAMMERHEAD")
```

## Documentation

For more information, see the [documentation][docs-stable-url].

[comment]: <>  UPDATE WITH OUR DOCS

[docs-dev-url]: https://jmurphy6895.github.io/HAMMERHEAD.jl/stable/
[docs-stable-url]: https://jmurphy6895.github.io/HAMMERHEAD.jl/stable/
