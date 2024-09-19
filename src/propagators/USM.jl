export USM7_ODE, USM7_ODE!

function USM7_ODE(
    u::AbstractArray,
    p::ComponentVector,
    t::Number,
    models::NTuple{N, AstroForceModels.AbstractAstroForceModel}) where N

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u
    μ = 

    sinλ = (2*ϵO3*η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    l = (ϵO1*ϵO3 - ϵO2*η0) / (ϵO3^2 + η0^2)
    ve2 = C - Rf1*sinλ + Rf2*cosλ
    ω3 = (C*ve2^2)/μ
    ρ = C / ve2

    u_cart = USM72cart(u, μ)

    fe = RTN_frame(u_cart) * build_dynamics_model(u, p, t, models)

    ω1 = fe[3] / ve2    

    return SVector{7}([-ρ*fe[2];
        fe[1]*cosλ - fe[2]*(1. + ρ)*sinλ - fe[3]*l*(Rf2/ve2);
        fe[1]*sinλ + fe[2]*(1. + ρ)*cosλ + fe[3]*l*(Rf1/ve2);
        .5*(ω3*ϵO2 + ω1*η0);
        .5*(-ω3*ϵO1 + ω1*ϵO3);
        .5*(-ω1*ϵO2 + ω3*η0);
        .5*(-ω1*ϵO1 - ω3*ϵO3)])
    
end

function USM7_ODE!(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    t::Number,
    models::NTuple{N, AstroForceModels.AbstractAstroForceModel}) where N

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u
    μ = 

    sinλ = (2*ϵO3*η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    l = (ϵO1*ϵO3 - ϵO2*η0) / (ϵO3^2 + η0^2)
    ve2 = C - Rf1*sinλ + Rf2*cosλ
    ω3 = (C*ve2^2)/μ
    ρ = C / ve2

    u_cart = USM72cart(u, μ)

    fe = RTN_frame(u_cart) * build_dynamics_model(u, p, t, models)

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