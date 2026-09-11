# Snow hydraulics parameterizations.

struct FreeDrainingHydrology{Tf} <: AbstractHydrology{Tf} end

@kwdef struct BucketHydrology{Tf} <: AbstractHydrology{Tf}
    Wirr::Tf = 0.03             # Irreducible liquid water content of snow (-)
end

struct DensityBucketHydrology{Tf} <: AbstractHydrology{Tf} end

FreeDrainingHydrology{Tf}(grid::Grid; kwargs...) where {Tf} = FreeDrainingHydrology{Tf}()
BucketHydrology{Tf}(grid::Grid; kwargs...) where {Tf} = BucketHydrology{Tf}(; kwargs...)
DensityBucketHydrology{Tf}(grid::Grid; kwargs...) where {Tf} = DensityBucketHydrology{Tf}()

"""
    snow_hydrology!(scheme, i, j, state, diag, params)

Route liquid water through the snow column at cell `(i, j)`: update `Sliq`,
`Sice`, `Tsnow`, `histowet` and the runoff/meltflux diagnostics in place, for
every layer.
"""
function snow_hydrology! end

# Free-draining snow
@inline function snow_hydrology!(::FreeDrainingHydrology{Tf}, i, j, state, diag, params) where {Tf}
    (; Sliq, Nsnow) = state
    (; Roff_snow, meltflux_out) = diag
    meltflux_out[i, j] = Tf(0)
    for k in 1:Nsnow[i, j]
        Roff_snow[i, j] = Roff_snow[i, j] + Sliq[k, i, j]
        meltflux_out[i, j] = meltflux_out[i, j] + Sliq[k, i, j]
        Sliq[k, i, j] = Tf(0)
    end
    return nothing
end

# Bucket storage
@inline function snow_hydrology!(c::BucketHydrology{Tf}, i, j, state, diag, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, Tsnow, histowet, fsnow, Nsnow) = state
    (; Roff_snow, meltflux_out) = diag
    (; Wirr) = c
    for k in 1:Nsnow[i, j]
        phi = Tf(0.0)
        if (Ds[k, i, j] > eps(Tf))
            phi = Tf(1) - Sice[k, i, j] / (rho_ice * Ds[k, i, j] * fsnow[i, j])
        end
        SliqMax = fsnow[i, j] * rho_wat * Ds[k, i, j] * phi * Wirr
        Sliq[k, i, j] = Sliq[k, i, j] + Roff_snow[i, j]
        Roff_snow[i, j] = Tf(0)
        if (Sliq[k, i, j] > SliqMax)   # Liquid capacity exceeded and drain to next layer
            Roff_snow[i, j] = Sliq[k, i, j] - SliqMax
            Sliq[k, i, j] = SliqMax
            histowet[k, i, j] = Tf(1.0)
        end
        # Rescale areal heat capacity of snow after mass updates
        csnow = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
        coldcont = csnow * (Tm - Tsnow[k, i, j])
        if (coldcont > Tf(0))          # Liquid water can freeze
            dSice = min(Sliq[k, i, j], fsnow[i, j] * coldcont / Lf)
            Sliq[k, i, j] = Sliq[k, i, j] - dSice
            Sice[k, i, j] = Sice[k, i, j] + dSice
            meltflux_out[i, j] = meltflux_out[i, j] - dSice
            Tsnow[k, i, j] = Tsnow[k, i, j] + Lf * dSice / csnow / fsnow[i, j]
        end
    end

    if (meltflux_out[i, j] < Tf(0))
        meltflux_out[i, j] = Tf(0)
    end
    return nothing
end

# Density-dependent bucket storage
@inline function snow_hydrology!(::DensityBucketHydrology{Tf}, i, j, state, diag, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, Tsnow, histowet, fsnow, Nsnow) = state
    (; Roff_snow, meltflux_out) = diag
    for k in 1:Nsnow[i, j]
        SliqCap = Tf(0.0)
        if (Ds[k, i, j] > eps(Tf))
            rhos = Sice[k, i, j] / Ds[k, i, j] / fsnow[i, j]
            SliqCap = Tf(0.03) + Tf(0.07) * (Tf(1) - rhos / Tf(200))
            SliqCap = max(SliqCap, Tf(0.03))
        end
        SliqMax = SliqCap * Sice[k, i, j]
        Sliq[k, i, j] = Sliq[k, i, j] + Roff_snow[i, j]
        Roff_snow[i, j] = Tf(0)
        if (Sliq[k, i, j] > SliqMax)   # Liquid capacity exceeded and drain to next layer
            Roff_snow[i, j] = Sliq[k, i, j] - SliqMax
            Sliq[k, i, j] = SliqMax
            histowet[k, i, j] = Tf(1.0)
        end
        # Rescale areal heat capacity of snow after mass updates
        csnow = (Sice[k, i, j] * hcap_ice + Sliq[k, i, j] * hcap_wat) / fsnow[i, j]
        coldcont = csnow * (Tm - Tsnow[k, i, j])
        if (coldcont > eps(Tf))        # Liquid water can freeze
            dSice = min(Sliq[k, i, j], fsnow[i, j] * coldcont / Lf)
            Sliq[k, i, j] = Sliq[k, i, j] - dSice
            Sice[k, i, j] = Sice[k, i, j] + dSice
            # Account for refreezing for melt and snow temperature
            meltflux_out[i, j] = meltflux_out[i, j] - dSice
            Tsnow[k, i, j] = Tsnow[k, i, j] + Lf * dSice / csnow / fsnow[i, j]
        end
    end

    if (meltflux_out[i, j] < Tf(0))
        meltflux_out[i, j] = Tf(0)
    end
    return nothing
end
