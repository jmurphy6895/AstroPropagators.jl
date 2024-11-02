@testset "Milankovich Propagator Keplerian" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    grav_model = KeplerianGravityAstroModel()
    p = ComponentVector(; JD=JD, μ=grav_model.μ)

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    u0_Mil = Array(AstroCoords.cart2Mil(u0, p.μ))

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = Milankovich_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_Mil, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)
    states = mapreduce(permutedims, vcat, sol.u)
    NRG = [
        (norm(state[4:6])^2.0) / 2.0 - AstroForceModels.μ_EARTH / norm(state[1:3]) for
        state in AstroCoords.Mil2cart.(sol.u, p.μ)
    ]

    @test NRG[1] ≈ NRG[end]

    # Comparison Against Cowell
    expected_end = [
        29447.829228927185,
        21027.31807398145,
        -1675.1455650449682,
        0.1548633780399451,
        2.381456403719347,
        0.14019776429122285,
    ]
    @test AstroCoords.Mil2cart(sol.u[end], p.μ) ≈ expected_end
end

@testset "Milankovich Propagator High-Fidelity" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    p = ComponentVector(;
        JD=JD, μ=GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
    )

    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    satellite_srp_model = CannonballFixedSRP(0.2)
    srp_model = SRPAstroModel(satellite_srp_model, sun_third_body, eop_data, Conical())

    satellite_drag_model = CannonballFixedDrag(0.2)
    drag_model = DragAstroModel(satellite_drag_model, JB2008(), eop_data)

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s
    u0_Mil = Array(Milankovich(Cartesian(u0), p.μ))

    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = Milankovich_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_Mil, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

    # Comparison Against Cowell
    expected_end = [
        29212.218059568793,
        22213.569774646894,
        -1540.5554989791624,
        0.021304477578929348,
        2.305181821713249,
        0.15097194461767285,
    ]
    @test Cartesian(Milankovich(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end

@testset "Milankovich Propagator High-Fidelity 2" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    p = ComponentVector(;
        JD=JD, μ=GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
    )

    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    satellite_srp_model = CannonballFixedSRP(0.5)
    srp_model = SRPAstroModel(satellite_srp_model, sun_third_body, eop_data, Conical())

    satellite_drag_model = CannonballFixedDrag(0.2)
    drag_model = DragAstroModel(satellite_drag_model, JB2008(), eop_data)

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        8.956857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    u0_Mil = Array(Milankovich(Cartesian(u0), p.μ))

    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

    tspan = (0.0, 3 * 86400.0)

    EOM!(du, u, p, t) = Milankovich_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_Mil, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

    # Comparison Against Cowell
    expected_end = [
        -6786.287820442294
        -1796.8353557785467
        578.0942989518618
        4.233306215217419
        -8.074542523286715
        -1.0229403931458645
    ]
    @test Cartesian(Milankovich(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end
