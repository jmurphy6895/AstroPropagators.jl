# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Different Shadow Models used Mainly in SRP Calculation
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#TODO: REFERENCE
#   [1]
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export shadow_model

"""
    shadow_model(sat_pos::AbstractArray, sun_pos::AbstractArray, R_Sun::Number, R_Earth::Number, t::Number, ShadowModel::Symbol)
Computes the Lighting Factor of the Sun occur from the Umbra and Prenumbra of Earth's Shadow

# Arguments

- `sat_pos::AbstractArray`: The current satellite position.
- `sun_pos::AbstractArray`: The current Sun position.
- `R_Sun::Number`: The radius of the Sun.
- `R_Earth::Number`: The radius of the Earth.
- `t::Number`: The current time of the simulation

# Optional Arguments 

- `ShadowModel::Symbol`: The Earth shadow model to use. Current Options -- :Cylindrical, :Conical, :Conical_Simplified, :Cylinder

# Returns

- `SVector{3}{Number}`: Inertial acceleration from the 3rd body
"""

#TODO: Double Check Math
@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    R_Sun::Number,
    R_Earth::Number,
    t::Number,
    ShadowModel::Val{:Cylinder})

    sat_pos = @view(sat_pos[1:3])

    # Compute dot product between sun and satellite positions
    dp_sun_sat = dot(normalize(sun_pos), sat_pos)

    return (dp_sun_sat > 0.0 || norm(sat_pos - dp_sun_sat * normalize(sun_pos)) > R_Earth) ? 1.0 : 0.0

end

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    R_Sun::Number,
    R_Earth::Number,
    t::Number,
    ShadowModel::Val{:Conical_Simplified})

    R_spacecraft_Sun = sat_pos - sun_pos

    con_a = asin(R_Sun/norm(R_spacecraft_Sun))
    con_b = asin(R_Earth/norm(sat_pos))
    con_c = AngleBetweenVectors.angle(normalize(@view(R_spacecraft_Sun[1:3])), normalize(@view(sat_pos[1:3])))

    return if con_c ≥ (con_b + con_a)
        1.0
    elseif con_c < (con_b - con_a)
        0.0
    else
        .5 + (con_c - con_b) / (2.0 * con_a)
    end

end

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    R_Sun::Number,
    R_Earth::Number,
    t::Number,
    ShadowModel::Val{:Conical})

    R_spacecraft_Sun = sat_pos - sun_pos

    con_a = asin(R_Sun/norm(R_spacecraft_Sun))
    con_b = asin(R_Earth/norm(sat_pos))
    con_c = AngleBetweenVectors.angle(normalize(@view(R_spacecraft_Sun[1:3])), normalize(@view(sat_pos[1:3])))

    return if con_c ≥ (con_b + con_a)
        1.0
    elseif con_c < (con_b - con_a)
        0.0
    elseif con_c < (con_a - con_b)
        1.0 - (con_b/con_a)^2
    else
        x = (con_c^2 + con_a^2 - con_b^2)/(2.0 * con_c)
        y = √(con_a^2 - x^2)
        area = con_a^2*acos(x/con_a) + con_b^2*acos((con_c - x)/con_b) - con_c*y
        1.0 - area/(π * con_a^2)
    end

end

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    R_Sun::Number,
    R_Earth::Number,
    t::Number,
    ShadowModel::Val{:None})

    return 1.0

end

@inline function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    R_Sun::Number,
    R_Earth::Number,
    t::Number)

    return shadow_model(sat_pos, sun_pos, R_Sun, R_Earth, t, :Conical)

end

@valsplit function shadow_model(
    sat_pos::AbstractArray,
    sun_pos::AbstractArray,
    R_Sun::Number,
    R_Earth::Number,
    t::Number,
    Val(ShadowModel::Symbol))

    error("Shadow Model is not Defined $ShadowModel")
    
end

#TODO: ADD 3RD BODY OCCULSION