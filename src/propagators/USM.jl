export USM7_ODE

function USM7_ODE(
    u::AbstractVector, 
    p::ComponentVector, 
    t::Number)

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u
    μ = p[1]

    sinλ = (2*ϵO3*η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
        
    l = (ϵO1*ϵO3 - ϵO2*η0) / (ϵO3^2 + η0^2)

    ve2 = C - Rf1*sinλ + Rf2*cosλ

    ω3 = (C*ve2^2)/μ

    ρ = C / ve2

    fe = 

    ω1 = fe[3] / ve2    

    dC = -ρ*fe[2]
    dRf1 = fe[1]*cosλ - fe[2]*(1. + ρ)*sinλ - fe[3]*l*(Rf2/ve2)
    dRf2 = fe[1]*sinλ + fe[2]*(1. + ρ)*cosλ + fe[3]*l*(Rf1/ve2)
    dϵO1 = .5*(ω3*ϵO2 + ω1*η0)
    dϵO2 = .5*(-ω3*ϵO1 + ω1*ϵO3)
    dϵO3 = .5*(-ω1*ϵO2 + ω3*η0)
    dη0 = .5*(-ω1*ϵO1 - ω3*ϵO3)

    return [dC; dRf1; dRf2; dϵO1; dϵO2; dϵO3; dη0]
    
end