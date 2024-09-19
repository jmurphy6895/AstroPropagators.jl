abstract type AbstractPropType end
struct Cowell <: AbstractPropType end

function propagate(
    u0::AbstractArray,
    p::ComponentArray,
    models::NTuple{N, AstroForceModels.AbstractAstroForceModel},
    tspan::Tuple{Number, Number};
    tsteps=nothing,
    ODE_solver=VCABM(),
    abstol=1E-13,
    reltol=1E-13,
    output_file::Union{String, Nothing}=nothing) where N

    EOM!(du, u, p, t) = 
        if prop_type == Cowell()
            Cowell_EOM!(du, u, p, t, models)
        end

    #TODO: Add in rest of propagators

    ODE_prob = ODEProblem(EOM!, u0, tspan, p)

    sol = solve(ODE_prob, ODE_solver, reltol=reltol, abstol=abstol)

    if output_file !== nothing
        times = (tsteps === nothing) ? sol.t : tsteps
        states = Array(sol(times))

        #TODO: print to file
        #TODO: Support Ephemeris Writing
    end

    return sol
end