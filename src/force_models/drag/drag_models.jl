# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from Drag
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#TODO: REFERENCE
#   [1] 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export drag_accel

"""
    drag_accel(u::AbstractArray, rho::Number, BC::Number, ω_vec::AbstractArray, t::Number, [DragModel]) -> SVector{3}{Number}

Compute the Acceleration Atmospheric Drag

The atmosphere is treated as a solid revolving with the Earth and the apparent velocity of the satellite is computed
using the transport theorem

                𝐯_app = 𝐯 - 𝛚 x 𝐫

The acceleration from drag is then computed with a cannonball model as
                
                𝐚 = 1/2 * ρ * BC * |𝐯_app|₂^2 * v̂


!!! note
    Currently only Cannonball Drag is supported, to use a higher fidelity drag either use a state varying function or compute
    the ballsitic coeffient further upstream

# Arguments

- `u::AbstractArray`: The current state of the spacecraft in the central body's inertial frame.
- `rho::Number`: Atmospheric Density at (t, u) [kg/m^3].
- `BC::Number`: The ballistic coeffient of the satellite -- (Area/mass) * Drag Coefficient [kg/m^2].
- `ω_vec::AbstractArray`: The angular velocity vector of Earth. Typically appozimated as [0.0; 0.0; ω_Earth]
- `t::Number`: Current Time of the Simulation.

# Optional Arguments

- `DragModel::Symbol`: Drag Model to use. Current Options -- :Cannonball, :None

# Returns

- `SVector{3}{Number}`: Inertial acceleration from Drag
"""
@inline function drag_accel(
    u::AbstractArray,
    rho::Number,
    BC::Number,
    ω_vec::AbstractArray,
    t::Number,
    DragModel::Val{:Cannonball})

    # Compute Apparent Velocity w.r.t the Atmosphere using the Transport Theorem
    apparent_vel = @view(u[4:6]) - cross(ω_vec, @view(u[1:3]))

    # Scaled by 1E3 to convert to km/s
    return SVector{3}((-.5 * BC * rho * norm(apparent_vel) * apparent_vel) .* 1E3)

end

@inline function drag_accel(
    u::AbstractArray,
    rho::Number,
    BC::Number,
    ω_vec::AbstractArray,
    t::Number,
    DragModel::Val{:None})

    return @SVector zeros(3)

end

@inline function drag_accel(
    u::AbstractArray,
    rho::Number,
    BC::Number,
    ω_vec::AbstractArray,
    t::Number) 
    
    return drag_accel(u, rho, BC, ω_vec, t, :Cannonball)

end

@valsplit function drag_accel(
    u::AbstractArray,
    rho::Number,
    BC::Number,
    ω_vec::AbstractArray,
    t::Number,
    Val(DragModel::Symbol))

    error("Model Not Defined for $DragModel")

end