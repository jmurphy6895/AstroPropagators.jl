@inline function compute_density(
    JD::T,
    u::AbstractArray,
    eop_data::EOPData_IAU1980,
    AtmosphereType::Val{:JB2008})

    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), p.JD + t/86400.0, eop_data)
    ecef_pos = R_J20002ITRF * @view(u[1:3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return (geodetic_pos[3] < 1000E3) * jb2008(JD, geodetic_pos...).total_density

end

@inline function compute_density(
    JD::T,
    u::AbstractArray,
    eop_data::EOPData_IAU1980,
    AtmosphereType::Val{:JR1971})

    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), p.JD + t/86400.0, eop_data)
    ecef_pos = R_J20002ITRF * @view(u[1:3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return (geodetic_pos[3] < 2500E3) * jr1971(JD, geodetic_pos...)

end

@inline function compute_density(
    JD::T,
    u::AbstractArray,
    eop_data::EOPData_IAU1980,
    AtmosphereType::Val{:MSIS2000})

    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), p.JD + t/86400.0, eop_data)
    ecef_pos = R_J20002ITRF * @view(u[1:3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return (geodetic_pos[3] < 1000E3) * nrlmsis00(JD, geodetic_pos...).total_density

end

@inline function compute_density(
    JD::T,
    u::AbstractArray,
    eop_data::EOPData_IAU1980,
    AtmosphereType::Val{:ExpAtmo})

    R_J20002ITRF = r_eci_to_ecef(DCM, J2000(), ITRF(), p.JD + t/86400.0, eop_data)
    ecef_pos = R_J20002ITRF * @view(u[1:3])
    geodetic_pos = ecef_to_geodetic(ecef_pos .* 1E3)

    return exponential(geodetic_pos[3])
    

end

@inline function compute_density(
    JD::T,
    u::AbstractArray,
    eop_data::EOPData_IAU1980,
    AtmosphereType::Val{:None})

    return 0.0

end

@valsplit function compute_density(
    JD::T,
    u::AbstractArray,
    eop_data::EOPData_IAU1980,
    Val(AtmosphereType::Symbol))

    error("Atmosphere Type Not Defined $AtmosphereType")

end