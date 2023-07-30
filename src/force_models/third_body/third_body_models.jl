# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Acceleration from a 3rd Body Represented as a Point Mass
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/propagation_setup/acceleration_models/third_body_acceleration.html
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export third_body_accel

"""
    third_body_accel(u::AbstractArray, μ_body::Number, body_pos::AbstractArray, h::Number) -> SVector{3}{Number}

Compute the Acceleration from a 3rd Body Represented as a Point Mass

Since the central body is also being acted upon by the third body, the acceleration of body 𝐁 acting on 
spacecraft 𝐀 in the orbiting body's 𝐂 is part of the force not acting on the central body

                a = ∇UB(rA) - ∇UB(rC)

# Arguments

- `u::AbstractArray`: The current state of the spacecraft in the central body's inertial frame.
- `μ_body`: Gravitation Parameter of the 3rd body.
- `body_pos::AbstractArray`: The current position of the 3rd body in the central body's inertial frame.

# Returns

- `SVector{3}{Number}`: Inertial acceleration from the 3rd body
"""
@inline function third_body_accel(
    u::AbstractArray,
    μ_body::Number,
    body_pos::AbstractArray) 

    # Compute Position Vectors for the Spacecraft w.r.t the Central and 3rd Body Respectively
    sat_pos = @view(u[1:3])
    r_spacecraft_to_body = body_pos - sat_pos

    # Calculate and Return the Acceleration from the Difference in Potential
    return SVector{3}(μ_body * (r_spacecraft_to_body/(norm(r_spacecraft_to_body)^3)) -
        μ_body * (body_pos/(norm(body_pos)^3)))

end