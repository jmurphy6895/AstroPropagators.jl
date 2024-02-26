export USM7_ODE

function USM7_ODE!(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    t::Number,
    grav_coeffs::AbstractGravityModel,
    eop_data::EOPData_IAU1980;
    max_order::Int=-1,
    max_degree::Int=-1,
    atmosphere_type::Symbol=:JR1971,
    drag_model::Symbol=:Cannonball,
    srp_model::Symbol=:Cannonball,
    shadow_model::Symbol=:Conical,
    lunar_3rd_body::Bool=true,
    solar_3rd_body::Bool=true)

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u
    μ = gravity_constant(grav_coeffs)

    sinλ = (2*ϵO3*η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
        
    l = (ϵO1*ϵO3 - ϵO2*η0) / (ϵO3^2 + η0^2)

    ve2 = C - Rf1*sinλ + Rf2*cosλ

    ω3 = (C*ve2^2)/μ

    ρ = C / ve2

    u_cart = USM72cart(u, μ)

    fe = RTN_frame(u_cart) * 
        (potential_perturb(
            u_cart, p, grav_coeffs, eop_data, t;
            max_degree=max_degree, max_order=max_order) +
        non_potential_accel(
            u_cart, p, t, eop_data; 
            atmosphere_type=atmosphere_type, 
            srp_model=srp_model, 
            drag_model=drag_model, 
            shadow_model=shadow_model, 
            lunar_3rd_body=lunar_3rd_body, 
            solar_3rd_body=solar_3rd_body))

    ω1 = fe[3] / ve2    

    du .= [-ρ*fe[2];
        fe[1]*cosλ - fe[2]*(1. + ρ)*sinλ - fe[3]*l*(Rf2/ve2);
        fe[1]*sinλ + fe[2]*(1. + ρ)*cosλ + fe[3]*l*(Rf1/ve2);
        .5*(ω3*ϵO2 + ω1*η0);
        .5*(-ω3*ϵO1 + ω1*ϵO3);
        .5*(-ω1*ϵO2 + ω3*η0);
        .5*(-ω1*ϵO1 - ω3*ϵO3)]

    return nothing
    
end