using AstroPropagators
using ComponentArrays
using OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqCore
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SpaceIndices
using Test

@testset "AstroPropagators.jl" begin
    include("propagators/test_cowell.jl")
end
