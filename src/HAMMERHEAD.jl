module HAMMERHEAD

using ValSplit
using StaticArrays
using LinearAlgebra, AngleBetweenVectors
using SatelliteToolboxBase, SatelliteToolboxGravityModels, SatelliteToolboxTransformations
using SatelliteToolboxAtmosphereModels, SpaceIndices

############################################################################################
#                                         Includes
############################################################################################
include("accelerations/non_potential_based_accel.jl")
include("accelerations/potential_based_accel.jl")

include("propagators/Cowell.jl")
include("propagators/EDROMO.jl")
include("propagators/Kustaanheimo-Stiefel.jl")
include("propagators/Stiefel-Scheifel.jl")
include("propagators/USM.jl")

include("api.jl")

end
