# Snow albedo parameterizations

@kwdef struct DiagnosticAlbedo{Tf} <: AbstractAlbedo{Tf}
    amin::Tf = 0.6                       # Minimum albedo for melting snow (-)
    amax::Tf = 0.86                      # Maximum albedo for fresh snow (-)
    Talb::Tf = -2                        # Albedo decay temperature threshold (C)
end

@kwdef struct DecayAlbedo{Tf} <: AbstractAlbedo{Tf}
    amin::Tf = 0.6                       # Minimum albedo for melting snow (-)
    tcld::Tf = 3600 * 1000               # Cold snow albedo decay time scale (s)
    tmlt::Tf = 3600 * 100                # Melting snow albedo decay time scale (s)
    adfs::Tf = 3                         # Albedo adjustment, shortwave (-)
    adfl::Tf = 2                         # Albedo adjustment, longwave (-)
    Sfmin::Tf = 10                       # Min 24h snowfall to refresh albedo (kg/m^2)
end

@kwdef struct PrognosticAlbedo{Tf} <: AbstractAlbedo{Tf}
    ALRADT::Bool = true                  # Aspect-dependent decay tuning
    adm::Tf = 100                        # Melting snow albedo decay time (h)
    amin::Tf = 0.6                       # Minimum albedo for melting snow (-)
    Sfmin::Tf = 10                       # Min 24h snowfall to refresh albedo (kg/m^2)
end

DiagnosticAlbedo{Tf}(Nx, Ny; kwargs...) where {Tf} = DiagnosticAlbedo{Tf}(; kwargs...)
DecayAlbedo{Tf}(Nx, Ny; kwargs...) where {Tf} = DecayAlbedo{Tf}(; kwargs...)
PrognosticAlbedo{Tf}(Nx, Ny; kwargs...) where {Tf} = PrognosticAlbedo{Tf}(; kwargs...)

"""
    snow_albedo!(scheme, i, j, state, landuse, meteo, params, summer_decay)

Update snow albedo `albs[i, j]` for cell `(i, j)` implemented for every 
`AbstractAlbedo`.
"""
function snow_albedo! end

@inline function snow_albedo!(c::DiagnosticAlbedo{Tf}, i, j, state, landuse, meteo, params, summer_decay) where {Tf}
    @unpack_constants(Tf)
    (; albs, Tsrf) = state
    afs_loc = c.amax
    a = c.amin + (afs_loc - c.amin) * (Tsrf[i, j] - Tm) / c.Talb
    a = max(a, min(afs_loc, c.amin))
    a = min(a, max(afs_loc, c.amin))
    albs[i, j] = a
    return nothing
end

@inline function snow_albedo!(c::DecayAlbedo{Tf}, i, j, state, landuse, meteo, params, summer_decay) where {Tf}
    @unpack_constants(Tf)
    (; albs, Tsrf) = state
    (; fveg, trcn, fsky, afs) = landuse
    (; Sdir, Sdif, Sf, Tv) = meteo
    (; dt) = params
    afs_loc = afs[i, j]

    tau = c.tcld
    if (Tsrf[i, j] >= Tm)
        tau = c.tmlt
    end
    
    if summer_decay
        tau = Tf(70.0) * Tf(3600.0)
    end

    if fveg[i, j] > Tf(0) && Sdir[i, j] > eps(Tf)
        tau = tau / ((Tf(1) - trcn[i, j] * fsky[i, j]) * (Tf(1) + c.adfl * Tv[i, j]) + c.adfs * Tv[i, j])
    elseif fveg[i, j] > Tf(0) && Sdif[i, j] > eps(Tf)
        tau = tau / ((Tf(1) - trcn[i, j] * fsky[i, j]) + c.adfs * trcn[i, j] * fsky[i, j])
    elseif (fveg[i, j] > Tf(0) && (Sdir[i, j] + Sdif[i, j] <= eps(Tf)))
        tau = tau / (Tf(2.0) - trcn[i, j] * fsky[i, j])
    end

    rt = Tf(1) / tau + Sf[i, j] / c.Sfmin
    alim = (c.amin / tau + Sf[i, j] * afs_loc / c.Sfmin) / rt
    a = alim + (albs[i, j] - alim) * exp(-rt * dt)
    if (a < min(afs_loc, c.amin))
        a = min(afs_loc, c.amin)
    end
    if (a > max(afs_loc, c.amin))
        a = max(afs_loc, c.amin)
    end
    albs[i, j] = a
    return nothing
end

@inline function snow_albedo!(c::PrognosticAlbedo{Tf}, i, j, state, landuse, meteo, params, summer_decay) where {Tf}
    @unpack_constants(Tf)
    (; albs, Tsrf, Sice, Sliq) = state
    (; afs, adc) = landuse
    (; Sdir, Sdird, Sf, Sf24h) = meteo
    (; dt) = params
    adc_loc = adc[i, j]
    adm_loc = c.adm
    afs_loc = afs[i, j]

    SWEtmp = zero(Tf)
    for si in 1:size(Sice, 1)
        SWEtmp += Sice[si, i, j] + Sliq[si, i, j]
    end

    # Aspect-dependent albedo tuning
    if c.ALRADT
        if ((Sdir[i, j] > eps(Tf)) && (Sdird[i, j] < Sdir[i, j]))
            adm_loc = adm_loc * (Sdird[i, j]) / (Sdir[i, j])
            adc_loc = adc_loc * (Sdird[i, j]) / (Sdir[i, j])
            if (adm_loc < eps(Tf))
                adm_loc = eps(Tf)
            end
            if (adc_loc < eps(Tf))
                adc_loc = eps(Tf)
            end
        end
    end

    # Temperature dependent albedo update
    a = albs[i, j]
    if (Tsrf[i, j] >= Tm)
        a = (a - c.amin) * exp(-(dt / Tf(3600)) / adm_loc) + c.amin
    else
        a = a - (dt / Tf(3600)) / adc_loc
    end

    # Reduce albedo for thin and patchy snow cover
    if (SWEtmp < Tf(75.0))
        afs_loc *= Tf(0.8)
    end

    # Reset to fresh snow albedo
    if ((Sf[i, j] * dt) > Tf(0.0) && Sf24h[i, j] > c.Sfmin)
        a = afs_loc
    else
        a = a + (afs_loc - a) * Sf[i, j] * dt / c.Sfmin
    end

    ## End Adjustments
    if (a > afs_loc)
        a = afs_loc
    end
    if (a < c.amin)
        a = c.amin
    end
    albs[i, j] = a
    return nothing
end

# Canopy radiative transfer

"""
    solar_radiation!(canopy, i, j, state, diag, landuse, meteo)

Surface albedo and shortwave transmission for cell `(i, j)`: fills `diag.alb`,
`diag.asrf_out`, `diag.SWsrf`, `diag.SWveg` and `diag.SWsci`. Expects `state.albs`
to already hold the bare-ground albedo where the snow has gone.
"""
function solar_radiation! end

@inline function solar_radiation!(c::NoCanopy{Tf}, i, j, state, diag, landuse, meteo) where {Tf}
    (; albs) = state
    (; alb, asrf_out, SWveg, SWsrf, SWsci) = diag
    (; Sdif, Sdir) = meteo

    asrf = albs[i, j]
    alb[i, j] = asrf
    asrf_out[i, j] = asrf
    SWveg[i, j] = Tf(0)
    SWsrf[i, j] = (Tf(1) - asrf) * (Sdir[i, j] + Sdif[i, j])
    SWsci[i, j] = Sdif[i, j] + Sdir[i, j]
    return nothing
end

@inline function solar_radiation!(c::OneLayerCanopy{Tf}, i, j, state, diag, landuse, meteo) where {Tf}
    (; albs, fsnow, Sveg) = state
    (; alb, asrf_out, SWveg, SWsrf, SWsci) = diag
    (; fveg, fsky, fsky_terr, scap, trcn) = landuse
    (; Sdif, Sdir, Tv) = meteo

    asrf = albs[i, j]
    if (fsnow[i, j] > eps(Tf))
        asrf *= Tf(1) - fveg[i, j] * canopy_fsar(c)
    end

    fcans = Tf(0.0)
    if (scap[i, j] > eps(Tf))
        fcans = Sveg[i, j] / scap[i, j]
    end
    aveg = (Tf(1) - fcans) * canopy_avg0(c) + fcans * canopy_avgs(c)
    acan = fveg[i, j] * aveg

    asrf_out[i, j] = fveg[i, j] * aveg + (Tf(1) - fveg[i, j]) * asrf

    Sdif_aux = fsky[i, j] / fsky_terr[i, j] * Sdif[i, j]
    tdif = trcn[i, j]
    tdir = Tv[i, j]
    alb[i, j] = acan + (Tf(1) - acan) * asrf * tdif^Tf(2)
    if (Sdif_aux + Sdir[i, j] > eps(Tf))
        alb[i, j] = (acan * (Sdif_aux + tdir * Sdir[i, j]) + asrf * tdif * (tdif * Sdif_aux + tdir * Sdir[i, j])) / (Sdif_aux + Sdir[i, j])
    end
    SWsrf[i, j] = (Tf(1) - asrf) * (tdif * Sdif_aux + tdir * Sdir[i, j])
    SWveg[i, j] = ((Tf(1) - tdif) * (Tf(1) - aveg) + tdif * asrf * (Tf(1) - tdif)) * Sdif_aux + (tdir * fveg[i, j] * (Tf(1) - aveg) + tdir * asrf * (Tf(1) - tdif)) * Sdir[i, j]   # local SWR absorption by vegetation correlates with local tdir
    SWsci[i, j] = tdif * Sdif_aux + tdir * Sdir[i, j]
    return nothing
end

"""
    thermal_radiation!(canopy, i, j, diag, landuse, meteo)

Effective incoming longwave `diag.LWeff` for cell `(i, j)`. Without canopy the terrain
emission is computed here, while in forested cells the process is accounted for in `ebalfor!`.
"""
function thermal_radiation! end

@inline function thermal_radiation!(c::NoCanopy{Tf}, i, j, diag, landuse, meteo) where {Tf}
    @unpack_constants(Tf)
    (; LWeff) = diag
    (; fsky_terr) = landuse
    (; LW, Ta) = meteo

    LWeff[i, j] = fsky_terr[i, j] * LW[i, j] + (Tf(1) - fsky_terr[i, j]) * sb * Ta[i, j]^Tf(4)
    return nothing
end

@inline function thermal_radiation!(c::OneLayerCanopy{Tf}, i, j, diag, landuse, meteo) where {Tf}
    (; LWeff) = diag
    (; LW) = meteo

    LWeff[i, j] = LW[i, j]
    return nothing
end

"""
    radiation!(fsm, meteo, t)

Snow albedo calculations, surface and canopy net shortwave radiation, terrain correction of
longwave radiation for open terrain.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions
- `t`: Current simulation time
"""
function radiation!(fsm::FSM{Tf, Ti}, meteo::MET{Tf, Ti}, t) where {Tf <: Real, Ti <: Integer}

    (; CANOPY, ALBEDO) = fsm.physics

    summer_decay = Dates.value(Month(t)) > 4 && Dates.value(Month(t)) < 10

    backend = get_backend(fsm.state.albs)
    kernel! = radiation_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.landuse, fsm.params, meteo,
        CANOPY, ALBEDO, summer_decay;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function radiation_kernel!(
        state, diag, landuse, params::Parameters{Tf}, meteo,
        CANOPY::AbstractCanopy{Tf},
        ALBEDO::AbstractAlbedo{Tf}, summer_decay::Bool,
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; tthresh) = params
    (; tilefrac, alb0) = landuse
    (; albs, fsnow) = state

    if (tilefrac[i, j] >= tthresh)

        # Snow albedo
        snow_albedo!(ALBEDO, i, j, state, landuse, meteo, params, summer_decay)

        # Bare ground shows through once the snow has gone
        if (fsnow[i, j] <= eps(Tf))
            albs[i, j] = alb0[i, j]
        end

        # Surface albedo, shortwave transmission and thermal emission from surroundings
        solar_radiation!(CANOPY, i, j, state, diag, landuse, meteo)
        thermal_radiation!(CANOPY, i, j, diag, landuse, meteo)

    end
end
