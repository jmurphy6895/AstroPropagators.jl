# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Solar Radiation Pressure
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] https://ai-solutions.com/_freeflyeruniversityguide/solar_radiation_pressure.htm
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export srp_accel

"""
    srp_accel(u::AbstractArray, sun_pos::AbstractArray, R_Sun::Number, R_Earth::Number, Ψ::Number, RC::Number, t::Number, [SRPModel::Symbol]; ShadowModel::Symbol)

Compute the Acceleration from Solar Radiaiton Pressure

Radiation from the Sun reflects off the satellite's surface and transfers momentum perturbing the satellite's trajectory. This
force can be computed using the a Cannonball model with teh following equation

                𝐚 = F * RC * Ψ * (AU/(R_sc_Sun))^2 * R̂_sc_Sun


!!! note
    Currently only Cannonball SRP is supported, to use a higher fidelity drag either use a state varying function or compute
    the ballsitic coeffient further upstream

# Arguments

- `u::AbstractArray`: The current state of the spacecraft in the central body's inertial frame.
- `sun_pos::AbstractArray`: The current position of the Sun.
- `R_Sun::Number`: The radius of the Sun.
- `R_Earth::Number`: The radius of the Earth.
- `Ψ::Number`: Solar Constant at 1 Astronomical Unit.
- `RC::Number`: The solar ballistic coeffient of the satellite -- (Area/mass) * Reflectivity Coefficient [kg/m^2].
- `t::Number`: The current time of the Simulation

# Optional Arguments

- `SRPModel::Symbol`: SRP Model to use. Current Options -- :Cannonball, :None
- `ShadowModel::Symbol`: SRP Earth Shadow Model to use. Current Options -- :Conical, :Conical_Simplified, :Cylinderical

# Returns

- `SVector{3}{Number}`: Inertial acceleration from the 3rd body
"""

function srp_accel(
    u::AbstractArray,
    sun_pos::AbstractArray,
    R_Sun::Number,
    R_Earth::Number,
    Ψ::Number,
    RC::Number,
    t::Number,
    SRPModel::Val{:Cannonball};
    ShadowModel::Symbol=:Conical)

    sat_pos = @view(u[1:3])

    # Compute the lighting factor
    F = shadow_model(sat_pos, sun_pos, R_Sun, R_Earth, t, ShadowModel)

    # Compute the Vector Between the Satellite and Sun
    R_spacecraft_Sun = sat_pos - sun_pos

    #Compute the SRP Force
    return SVector{3}((F * RC * Ψ * (ASTRONOMICAL_UNIT/norm(R_spacecraft_Sun))^2 * R_spacecraft_Sun/norm(R_spacecraft_Sun)) ./ 1E3)

end

@inline function srp_accel(
    u::AbstractArray,
    sun_pos::AbstractArray,
    Ψ::Number,
    RC::Number,
    t::Number,
    DragModel::Val{:None};
    ShadowModel::Symbol=:Conical)

    return @SVector zeros(3)

end


@inline function srp_accel(
    u::AbstractArray,
    sun_pos::AbstractArray,
    Ψ::Number,
    RC::Number,
    t::Number;
    ShadowModel::Symbol=:Conical) 
    
    return srp_accel(u, sun_pos, Ψ, RC, t, :Cannonball; ShadowModel=ShadowModel)

end

@valsplit function srp_accel(
    u::AbstractArray,
    sun_pos::AbstractArray,
    Ψ::Number,
    RC::Number,
    t::Number,
    Val(SRPModel::Symbol);
    ShadowModel::Symbol=:Conical)

    error("Model Not Defined for $SRPModel")

end