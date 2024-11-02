export USM7_EOM, USM7_EOM!

function USM7_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u

    μ = p.μ

    sinλ = (2 * ϵO3 * η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    l = (ϵO1 * ϵO3 - ϵO2 * η0) / (ϵO3^2 + η0^2)
    ve2 = C - Rf1 * sinλ + Rf2 * cosλ
    ω3 = (C * ve2^2) / μ
    ρ = C / ve2

    u_cart = Cartesian(USM7(u), p.μ)

    fe =
        RTN_frame(u_cart) * (
            build_dynamics_model(u_cart, p, t, models) -
            acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=p.μ))
        )

    ω1 = fe[3] / ve2

    return SVector{7}(
        [
            -ρ * fe[2]
            fe[1] * cosλ - fe[2] * (1.0 + ρ) * sinλ - fe[3] * l * (Rf2 / ve2)
            fe[1] * sinλ + fe[2] * (1.0 + ρ) * cosλ + fe[3] * l * (Rf1 / ve2)
            0.5 * (ω3 * ϵO2 + ω1 * η0)
            0.5 * (-ω3 * ϵO1 + ω1 * ϵO3)
            0.5 * (-ω1 * ϵO2 + ω3 * η0)
            0.5 * (-ω1 * ϵO1 - ω3 * ϵO3)
        ],
    )
end

function USM7_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    du .= USM7_EOM(u, p, t, models)

    return nothing
end

export USM6_EOM, USM6_EOM!
"""
USM6 ODE System

Arguments:
-'u::AbstractArray{AbstractFloat, 1}': USM6 State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
-'p::AbstractArray{AbstractFloat, 1}': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the orbital frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': USM6 ODE Change in State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
"""
function USM6_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    σ = SVector{3}(u[4], u[5], u[6])
    σ_norm = √(sum(abs2.(σ)))

    #TODO: SHADOW SET SHOULD PROBABLY BE EVENT
    if σ_norm > 1.0
        u[4:6] .= -u[4:6] / σ_norm
    end

    C, Rf1, Rf2, σ1, σ2, σ3 = u
    μ = p.μ

    σ = SVector{3}(u[4], u[5], u[6])
    σ_norm = √(sum(abs2.(σ)))

    _, _, _, ϵO1, ϵO2, ϵO3, η0 = USM7(USM6(u), μ)
    u_cart = Cartesian(USM6(u), μ)

    sinλ = (4.0 * σ3 * (1.0 - σ_norm^2)) / (4.0 * σ3^2 + (1 - σ_norm^2)^2)
    cosλ = ((1.0 - σ_norm^2)^2 - 4.0 * σ3^2) / (4.0 * σ3^2 + (1 - σ_norm^2)^2)

    l = (ϵO1 * ϵO3 - ϵO2 * η0) / (ϵO3^2 + η0^2)

    ve2 = C - Rf1 * sinλ + Rf2 * cosλ

    ω3 = (C * ve2^2) / μ

    ρ = C / ve2

    fe =
        RTN_frame(u_cart) * (
            build_dynamics_model(u_cart, p, t, models) -
            acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=μ))
        )

    ω1 = fe[3] / ve2

    dC = -ρ * fe[2]
    dRf1 = fe[1] * cosλ - fe[2] * (1.0 + ρ) * sinλ - fe[3] * l * (Rf2 / ve2)
    dRf2 = fe[1] * sinλ + fe[2] * (1.0 + ρ) * cosλ + fe[3] * l * (Rf1 / ve2)
    dσ1 = 0.25 * ((1.0 - σ_norm^2 + 2.0 * σ1^2) * ω1 + 2.0 * (σ1 * σ3 + σ2) * ω3)
    dσ2 = 0.25 * (2.0 * (σ2 * σ1 + σ3) * ω1 + 2.0 * (σ2 * σ3 - σ1) * ω3)
    dσ3 = 0.25 * (2.0 * (σ3 * σ1 - σ2) * ω1 + (1.0 - σ_norm^2 + 2.0 * σ3^2) * ω3)

    return SVector{6}(dC, dRf1, dRf2, dσ1, dσ2, dσ3)
end

"""
USM6 ODE System

Arguments:
-'u::AbstractArray{AbstractFloat, 1}': USM6 State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
-'p::AbstractArray{AbstractFloat, 1}': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the orbital frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': USM6 ODE Change in State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
"""
function USM6_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    du .= USM6_EOM(u, p, t, models)

    return nothing
end

export USMEM_EOM, USMEM_EOM!
"""
USMEM ODE System

Arguments:
-'u::AbstractArray{AbstractFloat, 1}': USM6 State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
-'p::AbstractArray{AbstractFloat, 1}': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the orbital frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': USM6 ODE Change in State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
"""
function USMEM_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    Φ_tol::Float64=1E-8,
) where {N}
    C, Rf1, Rf2, a1, a2, a3 = u

    a = SVector{3}(a1, a2, a3)
    Φ = √(sum(abs2.(a)))
    a_cross = skew_sym(a)

    μ = p.μ

    (_, _, _, ϵO1, ϵO2, ϵO3, η0) = USM7(USMEM(u), μ)
    u_cart = Cartesian(USMEM(u), μ)

    sinλ = (2 * ϵO3 * η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    l = (ϵO1 * ϵO3 - ϵO2 * η0) / (ϵO3^2 + η0^2)
    ve2 = C - Rf1 * sinλ + Rf2 * cosλ
    ω3 = (C * ve2^2) / μ
    ρ = C / ve2

    fe =
        RTN_frame(u_cart) * (
            build_dynamics_model(u_cart, p, t, models) -
            acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=p.μ))
        )

    ω1 = fe[3] / ve2

    ω = SVector{3}(ω1, 0.0, ω3)

    dC = -ρ * fe[2]
    dRf1 = fe[1] * cosλ - fe[2] * (1.0 + ρ) * sinλ - fe[3] * l * (Rf2 / ve2)
    dRf2 = fe[1] * sinλ + fe[2] * (1.0 + ρ) * cosλ + fe[3] * l * (Rf1 / ve2)
    if Φ > Φ_tol
        da =
            (
                I(3) +
                a_cross / 2.0 +
                1.0 / (Φ^2.0) * (1.0 - Φ / 2.0 * cot(Φ / 2.0)) * a_cross * a_cross
            ) * ω
    else
        da =
            0.5 * (
                (12.0 - Φ^2.0) / 6.0 * ω - cross(ω, a) -
                dot(ω, a) * ((60.0 + Φ^2.0) / 360.0) * a
            )
    end

    return SVector{6}(dC, dRf1, dRf2, da[1], da[2], da[3])
end

"""
USM6 ODE System

Arguments:
-'u::AbstractArray{AbstractFloat, 1}': USM6 State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
-'p::AbstractArray{AbstractFloat, 1}': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the orbital frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': USM6 ODE Change in State Vector [C; Rf1; Rf2; σ1; σ2; σ3]
"""
function USMEM_EOM!(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    Φ_tol::Float64=1E-8,
) where {N}
    du .= USMEM_EOM(u, p, t, models; Φ_tol=Φ_tol)

    return nothing
end
