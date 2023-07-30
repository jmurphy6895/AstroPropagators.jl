module HAMMERHEAD

using ValSplit
using StaticArrays
using LinearAlgebra, AngleBetweenVectors
using SatelliteToolboxBase, SatelliteToolboxGravityModels, SatelliteToolboxTransformations
using SatelliteToolboxAtmosphereModels, SpaceIndices

############################################################################################
#                                         Includes
############################################################################################

include("force_models/third_body/third_body_models.jl")

include("force_models/solar_radiation_pressure/shadow_models.jl")
include("force_models/solar_radiation_pressure/srp_model.jl")

include("force_models/drag/drag_models.jl")
include("force_models/drag/density_calculator.jl")

include("accelerations/non_potential_based_accel.jl")
include("accelerations/potential_based_accel.jl")

include("propagators/Cowell.jl")
include("propagators/EDROMO.jl")
include("propagators/Kustaanheimo-Stiefel.jl")
include("propagators/Stiefel-Scheifel.jl")
include("propagators/USM.jl")

include("api.jl")

end
