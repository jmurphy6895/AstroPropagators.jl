# AstroPropagators

[![CI](https://github.com/jmurphy6895/AstroPropagators.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jmurphy6895/AstroPropagators.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/jmurphy6895/AstroPropagators.jl/branch/main/graph/badge.svg?token=47G4OLV6PD)](https://codecov.io/gh/jmurphy6895/AstroPropagators.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

## Description

This project implements several propagation methods for the **HAMMERHEAD.jl** ecosystem using equations of motion for different orbital element sets. The force models can either be user-supplied or built up with the **AstroForceModels.jl** package. This package was inspired by the THALASSA library [CITE]. A full list of implemented propagators can be found below:

- [ ] Cowell
- [ ] EDromo
- [ ] Kustaanheimo-Stiefel
- [ ] Stiefel-Scheifel
- [ ] Unified State Model
- [ ] Gauss Variational Equations

## Installation

```julia
julia> using Pkg
julia> Pkg.add("AstroPropagators")
```

## Documentation

For more information, see the [documentation][docs-stable-url].

[comment]: <>  UPDATE WITH OUR DOCS

[docs-dev-url]: https://jmurphy6895.github.io/AstroPropagators.jl/stable/
[docs-stable-url]: https://jmurphy6895.github.io/AstroPropagators.jl/stable/
