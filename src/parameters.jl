# Generic physical constants that adapt to the precision type
function get_constants(::Type{Tf}) where {Tf <: AbstractFloat}
    return (
        cp = Tf(1005),         # Specific heat capacity of air (J/K/kg)
        eps_fsm = Tf(0.622),   # Ratio of molecular weights of water and dry air
        e0 = Tf(610.78),       # Saturation vapour pressure at Tm (Pa)
        grav = Tf(9.81),       # Acceleration due to gravity (m/s^2)
        hcap_ice = Tf(2100),   # Specific heat capacity of ice (J/K/kg)
        hcap_wat = Tf(4180),   # Specific heat capacity of water (J/K/kg)
        hcon_air = Tf(0.025),  # Thermal conductivity of air (W/m/K)
        hcon_clay = Tf(1.16),  # Thermal conductivity of clay (W/m/K)
        hcon_ice = Tf(2.24),   # Thermal conducivity of ice (W/m/K)
        hcon_sand = Tf(1.57),  # Thermal conductivity of sand (W/m/K)
        hcon_wat = Tf(0.56),   # Thermal conductivity of water (W/m/K)
        I0 = Tf(1367),         # Solar constant (W/m^2)
        Lf = Tf(0.334e6),      # Latent heat of fusion (J/kg)
        Lv = Tf(2.501e6),      # Latent heat of vapourisation (J/kg)
        Ls = Tf(0.334e6 + 2.501e6), # Latent heat of sublimation (J/kg)
        pi = Tf(3.14159),      # pi
        Rair = Tf(287),        # Gas constant for air (J/K/kg)
        Rwat = Tf(462),        # Gas constant for water vapour (J/K/kg)
        rho_ice = Tf(917),     # Density of ice (kg/m^3)
        rho_wat = Tf(1000),    # Density of water (kg/m^3)
        sb = Tf(5.67e-8),      # Stefan-Boltzmann constant (W/m^2/K^4)
        em_snow = Tf(0.99),    # Emissivity snow for Stefan-Boltzmann
        em_soil = Tf(0.9),    # Emissivity soil for Stefan-Boltzmann
        Tm = Tf(273.15),       # Melting point (K)
        vkman = Tf(0.4),        # Von Karman constant
    )
end

# Convenience macro to unpack constants in functions
macro unpack_constants(Tf)
    return esc(
        quote
            constants = get_constants($Tf)
            cp = constants.cp
            eps_fsm = constants.eps_fsm
            e0 = constants.e0
            grav = constants.grav
            hcap_ice = constants.hcap_ice
            hcap_wat = constants.hcap_wat
            hcon_air = constants.hcon_air
            hcon_clay = constants.hcon_clay
            hcon_ice = constants.hcon_ice
            hcon_sand = constants.hcon_sand
            hcon_wat = constants.hcon_wat
            I0 = constants.I0
            Lf = constants.Lf
            Lv = constants.Lv
            Ls = constants.Ls
            pi = constants.pi
            Rair = constants.Rair
            Rwat = constants.Rwat
            rho_ice = constants.rho_ice
            rho_wat = constants.rho_wat
            sb = constants.sb
            em_snow = constants.em_snow
            em_soil = constants.em_soil
            Tm = constants.Tm
            vkman = constants.vkman
        end
    )
end
