@testset "USM7 Propagator Keplerian" begin
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

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = USM7_EOM!(du, u, p, t, model_list)

    u0_USM7 = Array(USM7(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_USM7, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

    NRG = orbitalNRG.(USM7.(sol.u), p.μ)

    @test NRG[1] ≈ NRG[end]

    # Comparison Against Cowell
    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test Cartesian(USM7(sol.u[end]), p.μ) ≈ expected_end
end

@testset "USM7 Propagator High-Fidelity" begin
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

    u0_USM7 = Array(USM7(Cartesian(u0), p.μ))

    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = USM7_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_USM7, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

    # Comparison Against Cowell
    expected_end = [
        29212.62374256986
        22212.76353663187
        -1540.6511103523974
        0.021438031925788435
        2.3052533789977283
        0.15096291364452943
    ]
    @test Cartesian(USM7(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end

@testset "USM7 Propagator High-Fidelity 2" begin
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
        8.956857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    u0_USM7 = Array(USM7(Cartesian(u0), p.μ))

    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

    tspan = (0.0, 3.0 * 86400.0)

    EOM!(du, u, p, t) = USM7_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_USM7, tspan, p)
    sol = solve(prob, VCABM(); adaptive=true, abstol=1e-13, reltol=1e-13)

    # Regression Test
    # TODO: Cowell and USM diverge after a while, need better tests
    expected_end = [
        5568.572766761747
        -6677.611439306551
        -1061.0170203608875
        7.488631356803176
        2.1801419660996624
        -0.6222432860620839
    ]
    @test Cartesian(USM7(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end

@testset "USM6 Propagator Keplerian" begin
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

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = USM6_EOM!(du, u, p, t, model_list)

    u0_USM6 = Array(USM6(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_USM6, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

    NRG = orbitalNRG.(USM6.(sol.u), p.μ)

    @test NRG[1] ≈ NRG[end]

    # Comparison Against Cowell
    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test Cartesian(USM6(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end

@testset "USM6 Propagator High-Fidelity" begin
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

    u0_USM6 = Array(USM6(Cartesian(u0), p.μ))

    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = USM6_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_USM6, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

    # Comparison Against Cowell
    expected_end = [
        29212.62374256986
        22212.76353663187
        -1540.6511103523974
        0.021438031925788435
        2.3052533789977283
        0.15096291364452943
    ]
    @test Cartesian(USM6(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end

@testset "USM6 Propagator High-Fidelity 2" begin
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
        8.956857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    u0_USM6 = Array(USM6(Cartesian(u0), p.μ))

    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

    tspan = (0.0, 3 * 86400.0)

    EOM!(du, u, p, t) = USM6_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_USM6, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

    # Regression Test
    # TODO: Cowell and USM diverge after a while, need better tests
    # Uses USM7 result
    expected_end = [
        5568.572766761747
        -6677.611439306551
        -1061.0170203608875
        7.488631356803176
        2.1801419660996624
        -0.6222432860620839
    ]
    @test Cartesian(USM6(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-3
end

@testset "USMEM Propagator Keplerian" begin
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

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = USMEM_EOM!(du, u, p, t, model_list)

    u0_USMEM = Array(USMEM(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_USMEM, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

    NRG = orbitalNRG.(USMEM.(sol.u), p.μ)

    @test NRG[1] ≈ NRG[end]

    #Comparison Against Cowell
    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test Cartesian(USMEM(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-6
end

@testset "USMEM Propagator High-Fidelity" begin
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

    u0_USMEM = Array(USMEM(Cartesian(u0), p.μ))

    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = USMEM_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_USMEM, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

    # Regression Test
    expected_end = [
        29212.62374256986
        22212.76353663187
        -1540.6511103523974
        0.021438031925788435
        2.3052533789977283
        0.15096291364452943
    ]
    @test Cartesian(USMEM(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end

@testset "USMEM Propagator High-Fidelity 2" begin
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
        8.956857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    u0_USMEM = Array(USMEM(Cartesian(u0), p.μ))

    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

    tspan = (0.0, 3 * 86400.0)

    EOM!(du, u, p, t) = USMEM_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_USMEM, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-15, reltol=1e-15)

    # Regression Test
    # TODO: Cowell and USM diverge after a while, need better tests
    # Uses USM7 result
    expected_end = [
        5568.572766761747
        -6677.611439306551
        -1061.0170203608875
        7.488631356803176
        2.1801419660996624
        -0.6222432860620839
    ]
    @test Cartesian(USMEM(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end
