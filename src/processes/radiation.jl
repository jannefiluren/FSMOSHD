# Radiation process: snow-albedo application, canopy shortwave transmission, terrain longwave

"""
    radiation!(fsm, meteo, t)

Snow albedo calculations, surface and canopy net shortwave radiation, terrain correction of
longwave radiation for open terrain.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions
- `t`: Current simulation time
"""
function radiation!(fsm::FSM{Tf}, meteo::MET{Tf}, t) where {Tf <: Real}

    (; canopy, snow_albedo) = fsm.physics

    summer_decay = Dates.value(Month(t)) > 4 && Dates.value(Month(t)) < 10

    backend = get_backend(fsm.state.albs)
    kernel! = radiation_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.params, meteo,
        canopy, snow_albedo, summer_decay;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function radiation_kernel!(
        state, diag, surface, params::Parameters{Tf}, meteo,
        canopy::AbstractCanopy{Tf},
        snow_albedo::AbstractAlbedo{Tf}, summer_decay::Bool,
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; tthresh) = params
    (; tilefrac, alb0) = surface
    (; albs, fsnow) = state

    if (tilefrac[i, j] >= tthresh)

        # Snow albedo
        snow_albedo!(snow_albedo, i, j, state, surface, meteo, params, summer_decay)

        # Bare ground shows through once the snow has gone
        if (fsnow[i, j] <= eps(Tf))
            albs[i, j] = alb0[i, j]
        end

        # Surface albedo, shortwave transmission and thermal emission from surroundings
        solar_radiation!(canopy, i, j, state, diag, surface, meteo)
        thermal_radiation!(canopy, i, j, diag, surface, meteo)

    end
end

# Canopy radiative transfer

"""
    solar_radiation!(canopy, i, j, state, diag, surface, meteo)

Surface albedo and shortwave transmission for cell `(i, j)`: fills `diag.alb`,
`diag.asrf_out`, `diag.SWsrf`, `diag.SWveg` and `diag.SWsci`. Expects `state.albs`
to already hold the bare-ground albedo where the snow has gone.
"""
function solar_radiation! end

@inline function solar_radiation!(c::NoCanopy{Tf}, i, j, state, diag, surface, meteo) where {Tf}
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

@inline function solar_radiation!(c::OneLayerCanopy{Tf}, i, j, state, diag, surface, meteo) where {Tf}
    (; albs, fsnow, Sveg) = state
    (; alb, asrf_out, SWveg, SWsrf, SWsci) = diag
    (; fveg, fsky, fsky_terr, scap, trcn) = surface
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
    thermal_radiation!(canopy, i, j, diag, surface, meteo)

Effective incoming longwave `diag.LWeff` for cell `(i, j)`. Without canopy the terrain
emission is computed here, while in forested cells `energy_balance!` accounts for it.
"""
function thermal_radiation! end

@inline function thermal_radiation!(c::NoCanopy{Tf}, i, j, diag, surface, meteo) where {Tf}
    @unpack_constants(Tf)
    (; LWeff) = diag
    (; fsky_terr) = surface
    (; LW, Ta) = meteo

    LWeff[i, j] = fsky_terr[i, j] * LW[i, j] + (Tf(1) - fsky_terr[i, j]) * sb * Ta[i, j]^Tf(4)
    return nothing
end

@inline function thermal_radiation!(c::OneLayerCanopy{Tf}, i, j, diag, surface, meteo) where {Tf}
    (; LWeff) = diag
    (; LW) = meteo

    LWeff[i, j] = LW[i, j]
    return nothing
end
