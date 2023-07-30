function EOM!(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    t::Number,
    grav_coeffs::AbstractGravityModel,
    eop_data::EOPData_IAU1980;
    max_order::Int=-1,
    max_degree::Int=-1,
    atmosphere_type::Symbol=:,
    drag_model::Symbol=:Cannonball,
    srp_model::Symbol=:Cannonball,
    shadow_model::Symbol=:Conical,
    lunar_3rd_body::Bool=true,
    solar_3rd_body::Bool=true)

    du .= [@view(u[4:6]);
        potential_accel(u, p, grav_coeffs, eop_data, t; 
            max_order=max_order, 
            max_degree=max_degree) +
        non_potential_accel(u, p, t, eop_data; 
            atmosphere_type=atmosphere_type, 
            srp_model=srp_model, 
            drag_model=drag_model, 
            shadow_model=shadow_model, 
            lunar_3rd_body=lunar_3rd_body, 
            solar_3rd_body=solar_3rd_body)]

    return nothing

end