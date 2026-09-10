# Snow thermal conductivity parameterizations.

@kwdef struct FixedConductivity{Tf} <: AbstractConductivity{Tf}
    kfix::Tf = 0.24        # Fixed thermal conductivity of snow (W/m/K)
end

@kwdef struct DensityConductivity{Tf} <: AbstractConductivity{Tf}
    bthr::Tf = 2           # Snow thermal conductivity exponent (-)
end

FixedConductivity{Tf}(grid::Grid; kwargs...) where {Tf} = FixedConductivity{Tf}(; kwargs...)
DensityConductivity{Tf}(grid::Grid; kwargs...) where {Tf} = DensityConductivity{Tf}(; kwargs...)

"""
    snow_conductivity!(scheme, i, j, state, diag, params)

Fill the snow thermal conductivity `ksnow[1:Nsnow, i, j]` for cell `(i, j)`,
implemented for every `AbstractConductivity`.
"""
function snow_conductivity! end

@inline function snow_conductivity!(c::FixedConductivity, i, j, state, diag, params)
    (; Nsnow) = state
    (; ksnow) = diag
    for k in 1:Nsnow[i, j]
        ksnow[k, i, j] = c.kfix
    end
    return nothing
end

@inline function snow_conductivity!(c::DensityConductivity{Tf}, i, j, state, diag, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow) = state
    (; ksnow) = diag
    (; rhof) = params
    for k in 1:Nsnow[i, j]
        rhos = rhof
        if ((Ds[k, i, j] > eps(Tf)) && fsnow[i, j] > eps(Tf))
            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
        end
        ksnow[k, i, j] = hcon_ice * (rhos / rho_ice)^c.bthr
    end
    return nothing
end
