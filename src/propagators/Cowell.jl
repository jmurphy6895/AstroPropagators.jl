export Cowell_EOM, Cowell_EOM!

function Cowell_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    accel = build_dynamics_model(u, p, t, models)
    return SVector{6}(u[4], u[5], u[6], accel[1], accel[2], accel[3])
end

function Cowell_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    du .= Cowell_EOM(u, p, t, models)

    return nothing
end
