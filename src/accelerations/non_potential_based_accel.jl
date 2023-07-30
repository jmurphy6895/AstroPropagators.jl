function non_potential_accel(
    u::AbstractArray,
    p::ComponentVector,
    t::Number,
    eop_data::EOPData_IAU1980;
    drag_model::Symbol=:Cannonball,
    srp_model::Symbol=:Cannonball,
    shadow_model::Symbol=:Conical,
    lunar_3rd_body::Bool=true,
    solar_3rd_body::Bool=true)

    JD = p.JD + t/86400.0

    #* Third Body Accelerations
    R_MOD2J2000 = r_eci_to_eci(MOD(), J2000(), JD, eop_data)

    moon_pos = R_MOD2J2000 * moon_position_i(JD) ./ 1E3
    moon_3rd_body_accel = lunar_3rd_body * third_body_accel(u, p.μ_Moon, moon_pos)

    sun_pos = R_MOD2J2000 * sun_position_i(JD) ./ 1E3
    sun_3rd_body_accel = solar_3rd_body * third_body_accel(u, p.μ_Sun, sun_pos)

    #* Drag Acceleration
    rho = density_calculator()
    drag_accel = drag_accel(u, rho, p.BC, [0.0; 0.0; ω_Earth], t, drag_model)

    #* SRP Acceleration
    srp_accel = srp_accel(u, sun_pos, p.R_Sun, p.R_Earth, p.Ψ, p.RC, t, srp_model; ShadowModel=shadow_model)

    return SVector{3}(moon_3rd_body_accel + sun_3rd_body_accel + drag_accel + srp_accel)

end