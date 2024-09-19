using AstroPropagators
using Aqua
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqCore
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SpaceIndices
using Test

@testset "AstroPropagators.jl" begin
    include("propagators/test_cowell.jl")
end

@testset "Aqua.jl" begin
    Aqua.test_all(AstroForceModels; ambiguities=(recursive = false))
end
