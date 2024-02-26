function propagate(
    u0::AbstractArray,
    p::ComponentArray,
    tspan::Tuple{Number, Number};
    tsteps=nothing,
    ODE_solver=VCABM(),
    abstol=1E-13,
    reltol=1E-13,
    grav_model::AbstractGravityModel=load(IcgemFile, fetch_icgem_file(:EGM2008)),
    eop_data::EOPData_IAU1980=get_iers_eop(),
    prop_type::Symbol=:DROMO,
    max_order::Int=72,
    max_degree::Int=72,
    atmosphere_type::Symbol=:JR1971,
    drag_model::Symbol=:Cannonball,
    srp_model::Symbol=:Cannonball,
    shadow_model::Symbol=:Conical,
    lunar_3rd_body::Bool=true,
    solar_3rd_body::Bool=true,
    output_file::Union{String, Nothing}=nothing)

    if (atmosphere_type != :None || atmosphere_type != :ExpAtmo)
        SpaceIndices.init()
    end

    EOM!(du, u, p, t) = 
        if prop_type == :Cowell
            Cowell_EOM!(du, u, p, t, grav_model, eop_data;
                max_degree=max_degree,
                max_order=max_order,
                atmosphere_type=atmosphere_type,
                drag_model=drag_model,
                srp_model=srp_model,
                shadow_model=shadow_model,
                lunar_3rd_body=lunar_3rd_body,
                solar_3rd_body=solar_3rd_body)
        elseif prop_type == :EDROMO
            EDROMO_EOM!(du, u, p, t, grav_model, eop_data;
                max_degree=max_degree,
                max_order=max_order,
                atmosphere_type=atmosphere_type,
                drag_model=drag_model,
                srp_model=srp_model,
                shadow_model=shadow_model,
                lunar_3rd_body=lunar_3rd_body,
                solar_3rd_body=solar_3rd_body)
        elseif prop_type == :KS
            KS_EOM!(du, u, p, t, grav_model, eop_data;
                max_degree=max_degree,
                max_order=max_order,
                atmosphere_type=atmosphere_type,
                drag_model=drag_model,
                srp_model=srp_model,
                shadow_model=shadow_model,
                lunar_3rd_body=lunar_3rd_body,
                solar_3rd_body=solar_3rd_body)
        elseif prop_type == :STI_SCHE
            STI_SCHE_EOM!(du, u, p, t, grav_model, eop_data;
                max_degree=max_degree,
                max_order=max_order,
                atmosphere_type=atmosphere_type,
                drag_model=drag_model,
                srp_model=srp_model,
                shadow_model=shadow_model,
                lunar_3rd_body=lunar_3rd_body,
                solar_3rd_body=solar_3rd_body)
        elseif prop_type == :USM
            USM7_EOM!(du, u, p, t, grav_model, eop_data;
                max_degree=max_degree,
                max_order=max_order,
                atmosphere_type=atmosphere_type,
                drag_model=drag_model,
                srp_model=srp_model,
                shadow_model=shadow_model,
                lunar_3rd_body=lunar_3rd_body,
                solar_3rd_body=solar_3rd_body)
        end

        #TODO: CONVERT u0

        ODE_prob = ODEProblem(EOM!, u0, tspan, p)

        sol = solve(ODE_prob, ODE_solver, reltol=reltol, abstol=abstol)

        if output_file !== nothing
            times = (tsteps === nothing) ? sol.t : tsteps
            states = Array(sol(times))

            #TODO: print to file
        end

        return sol

end