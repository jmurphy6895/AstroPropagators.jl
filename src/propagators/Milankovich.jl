export Milankovich_EOM, Milankovich_EOM!
"""
Milankovich ODE System

Arguments:
-'u::AbstractArray': USM State Vector [H; e; L]
-'p::AbstractArray': Parameter Vector
-'t::AbstractFloat': Time
-'accel::Function': Function computing perturbing acceleration in the inertial frame 

Returns:
-'du::AbstractArray{AbstractFloat, 1}': USM ODE Change in State Vector [dH; de; dL]
"""
function Milankovich_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    Hx, Hy, Hz, _, _, _, _ = u
    μ = p.μ

    u_cart = AstroCoords.Mil2cart(u, μ)

    r = SVector{3}(u_cart[1], u_cart[2], u_cart[3])
    v = SVector{3}(u_cart[4], u_cart[5], u_cart[6])

    ad =
        build_dynamics_model(u_cart, p, t, models) -
        acceleration(u_cart, p, t, KeplerianGravityAstroModel(μ))

    H = SVector{3}(Hx, Hy, Hz)
    ẑ = SVector{3}(0.0, 0.0, 1.0)

    dH = skew_sym(r) * ad
    de = 1 / μ * (skew_sym(v) * skew_sym(r) - skew_sym(H)) * ad
    dL =
        (dot(ẑ, r) / (norm(H) * (norm(H) + dot(ẑ, H)))) * dot(H, ad) +
        norm(H) / (norm(r)^2)

    return SVector{7}(dH[1], dH[2], dH[3], de[1], de[2], de[3], dL)
end

function Milankovich_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    du .= Milankovich_EOM(u, p, t, models)

    return nothing
end
