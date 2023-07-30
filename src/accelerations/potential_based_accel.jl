function potential_accel(
    u::AbstractArray,
    p::ComponentArray,
    grav_model::AbstractGravityModel{T},
    eop_data::EOPData_IAU1980,
    t::Number;
    max_degree::Int=-1,
    max_order::Int=-1) where T

    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), p.JD + t/86400.0, eop_data)

    ecef_pos = R_J20002ITRF * @view(u[1:3])

    #TODO: FIGURE OUT HOW TO SUPPORT TIME VARYING COEFFS
    return R_J20002ITRF' * gravitational_acceleration(
        grav_model, ecef_pos;
        max_degree=max_degree, max_order=max_order) ./ 1E3 

end

function keplerian_accel(u, μ_Earth)

    r = @view(u[1:3])

    return SVector{3}(-μ_Earth/(norm(r)^3) * r)

end

function potential_perturb(
    u::AbstractArray,
    p::ComponentArray,
    grav_model::AbstractGravityModel{T},
    eop_data::EOPData_IAU1980,
    t::Number;
    max_degree::Int=-1,
    max_order::Int=-1) where T

    return potential_accel(
        u, p, grav_model, eop_data, t; 
        max_order=max_order, 
        max_degree=max_degree) - keplerian_accel(u, gravity_constant(grav_model))

end