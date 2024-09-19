export Cowell_EOM, Cowell_EOM!

function Cowell_EOM!(
    u::AbstractArray,
    p::ComponentVector,
    t::Number,
    models::NTuple{N, AstroForceModels.AbstractAstroForceModel}) where N

    return SVector{6}([@view(u[4:6]); build_dynamics_model(u, p, t, models)])

    return nothing

end

function Cowell_EOM!(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    t::Number,
    models::NTuple{N, AstroForceModels.AbstractAstroForceModel}) where N

    du[1:3] .= u[4:6]
    du[4:6] .= build_dynamics_model(u, p, t, models)

    return nothing

end