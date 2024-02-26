# HAMMERHEAD

[![Build Status](https://github.com/jmurphy6895/HAMMERHEAD.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jmurphy6895/HAMMERHEAD.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## Description

This project (**H**igh-**A**ccuracy **M**odelling **E**nvironment **R**esourcing **H**ighly-**E**fficient Solvers **A**nd **D**ifferentiability) implements high-fidelity Equations of Motion using Zonal Harmonics, Solar Radiation Pressure, Drag, and Third-Body Perturbations based on packages developed in the Satellite Toolbox.

These forces are then called from a number of high performance integrator methodologies. Some of these are ported and implemented from the THALASSA library [CITE]. A full list can be found below:

1. Cowell
2. EDromo
3. Kustaanheimo-Stiefel
4. Stiefel-Scheifel
5. Unified State Model
6. Gauss Variational Equations