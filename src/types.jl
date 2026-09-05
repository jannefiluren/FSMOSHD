@with_kw struct Grid{Ti, VF <: AbstractVector}
    Dzsnow::VF = [0.1, 0.2, 0.4]                     # Maximum snow layer thicknesses (m)
    Dzsoil::VF = [0.1, 0.2, 0.4, 0.8]                # Maximum soil layer thicknesses (m)
    Nsmax::Ti = length(Dzsnow)                       # Number of snow layers
    Nsoil::Ti = length(Dzsoil)                       # Number of soil layers
    Nx::Ti = 1                                       # First array dimension (rows)
    Ny::Ti = 1                                       # Second array dimension (columns)
end

"""
    check_layer_thicknesses(grid)

Verify that the first snow layer can grow at least as thick as the top soil layer.
The surface layer in `thermal!` blends snow and soil over a depth of `Dzsoil[1]`, so a
thinner first snow layer leaves `Ts1` and `ks1` contaminated by soil however deep the
snowpack gets.
"""
function check_layer_thicknesses(grid::Grid)
    grid.Dzsnow[1] >= grid.Dzsoil[1] || throw(
        ArgumentError(
            "Dzsnow[1] = $(grid.Dzsnow[1]) must be at least Dzsoil[1] = $(grid.Dzsoil[1]), " *
                "otherwise the surface layer in thermal! stays blended with the soil under deep snow"
        )
    )
    return nothing
end

@with_kw struct Parameters{Tf, Ti}
    dt::Tf = 3600                                    # Time step (s)
    zT::Tf = 10                                      # Temperature measurement height (m)
    zU::Tf = 10                                      # Wind speed measurement height (m)
    zRH::Tf = 10                                     # Relative humidity measurement height (m)
    tthresh::Tf = 0.1                                # Tile threshold
    Nitr::Ti = 4                                     # Iterations for surface energy balance
    cvai::Tf = 4.4                                   # Canopy snow capacity per unit vegetation area index (kg/m^2)
    Gcn1::Tf = 0.5                                   # Leaf angle distribution parameter (-)
    Gcn2::Tf = 0                                     # Leaf angle distribution parameter (-)
    gsnf::Tf = 0                                     # Snow-free vegetation moisture conductance (m/s)
    kdif::Tf = 0.5                                   # Diffuse radiation extinction coefficient (-)
    kveg::Tf = 1                                     # Canopy cover coefficient (-)
    tcnc::Tf = 3600 * 240                            # Canopy unloading time scale for cold snow (s)
    tcnm::Tf = 3600 * 48                             # Canopy unloading time scale for melting snow (s)
    pmultf_for::Tf = 0.5                             # Multiplier for snowfall in forest (-)
    a_eta::Tf = 0.1                                  # Temperature factor for Crocus B92 compaction (K^-1)
    b_eta::Tf = 0.023                                # First density factor for Crocus B92 compaction (m^3/kg)
    c_eta::Tf = 250                                  # Second density factor for Crocus B92 compaction (kg/m^3)
    eta0::Tf = 3.7e7                                 # Reference snow viscosity (Pa s)
    eta1::Tf = 7.62237e6                             # Reference snow viscosity for Crocus B92 compaction (Pa s)
    Tsnow_min::Tf = -Inf                             # Floor on snow layer temperature (K); -Inf disables
    hfsn::Tf = 0.1                                   # Snowcover fraction depth scale (m)
    rho0::Tf = 300                                   # Fixed snow density (kg/m^3)
    rhob::Tf = 6                                     # Temperature factor in fresh snow density (kg/m^3/K)
    rhoc::Tf = 26                                    # Wind factor in fresh snow density (kg s^0.5/m^3.5)
    rhof::Tf = 109                                   # Fresh snow density (kg/m^3)
    rhos_min::Tf = 50                                # Minimum snow density (kg/m^3)
    rhos_max::Tf = 750                               # Maximum snow density (kg/m^3)
    rcld::Tf = 300                                   # Maximum density for cold snow (kg/m^3)
    rgr0::Tf = 5.0e-5                                # Fresh snow grain radius (m)
    rmlt::Tf = 500                                   # Maximum density for melting snow (kg/m^3)
    Salb::Tf = 10                                    # Albedo decay constant (kg/m^2)
    snda::Tf = 2.8e-6                                # Thermal metamorphism parameter (1/s)
    trho::Tf = 3600 * 200                            # Snow compaction time scale (s)
    Wirr::Tf = 0.03                                  # Irreducible liquid water content of snow (-)
    Ds_min::Tf = 0.01                                # Minimum possible snow layer thickness (m)
    Ds_surflay::Tf = 0.5                             # Maximum thickness of surface fine snow layering (m)
    gsat::Tf = 0.01                                  # Surface conductance for saturated soil (m/s)
    zsub::Tf = 2                                     # Sub-canopy reference height (m)
    fsat::Tf = 0.5                                   # Initial soil moisture as fraction of saturation
    Tprof::Tf = 285                                  # Initial soil layer temperatures (K)
end

@with_kw struct Landuse{GT, MF <: AbstractMatrix{<:AbstractFloat}, MF64 <: AbstractMatrix{Float64}}
    grid::GT
    z0_snow::MF = 0.002 * ones(grid.Nx, grid.Ny)    # Roughness length of snow (m)
    alb0::MF = 0.2 * ones(grid.Nx, grid.Ny)         # Snow-free ground albedo (-)
    afs::MF = 0.86 * ones(grid.Nx, grid.Ny)         # Maximum albedo for fresh snow (-)  [Decay/Prognostic albedo]
    adc::MF = 1000 * ones(grid.Nx, grid.Ny)         # Cold snow albedo decay time (h)  [Prognostic albedo]
    z0sf::MF = 0.2 * ones(grid.Nx, grid.Ny)         # Snow-free roughness length (m)
    fcly::MF = 0.3 * ones(grid.Nx, grid.Ny)         # Soil clay fraction (-)
    fsnd::MF = 0.6 * ones(grid.Nx, grid.Ny)         # Soil sand fraction (-)
    VAI::MF = zeros(grid.Nx, grid.Ny)               # Vegetation area index (-)
    vfhp::MF = fill(NaN, grid.Nx, grid.Ny)          # Hemispherical sky-view fraction incl. canopy (-)
    canh::MF = fill(NaN, grid.Nx, grid.Ny)          # Canopy heat capacity (J/K/m^2)
    fsky::MF = ones(grid.Nx, grid.Ny)               # Sky view fraction (-)
    fveg::MF = zeros(grid.Nx, grid.Ny)              # Canopy cover fraction (-)  = 1-exp(-kveg*VAI), VAI=0
    fves::MF = zeros(grid.Nx, grid.Ny)              # Stand-scale canopy cover fraction (-)
    hcan::MF = zeros(grid.Nx, grid.Ny)              # Canopy height (m)
    lai::MF = fill(NaN, grid.Nx, grid.Ny)           # Leaf area index (-)
    pmultf::MF = fill(NaN, grid.Nx, grid.Ny)        # Precipitation multiplier reverting open-area correction (-)
    scap::MF = fill(NaN, grid.Nx, grid.Ny)          # Canopy snow capacity (kg/m^2)
    trcn::MF = ones(grid.Nx, grid.Ny)               # Canopy transmissivity (-)  = exp(-kdif*VAI), VAI=0
    slopemu::MF = fill(NaN, grid.Nx, grid.Ny)       # Slope parameter (-)
    xi::MF = fill(NaN, grid.Nx, grid.Ny)            # Terrain correlation length (m)
    Ld::MF = fill(NaN, grid.Nx, grid.Ny)            # Grid cell size (m)
    fsky_terr::MF = fill(NaN, grid.Nx, grid.Ny)     # Sky view fraction terrain (-)
    dem::MF = fill(NaN, grid.Nx, grid.Ny)           # Grid elevation (m)
    tilefrac::MF = ones(grid.Nx, grid.Ny)           # Tile fraction (-)
    prec_multi::MF64 = fill(NaN, grid.Nx, grid.Ny)  # Precipitation multiplier (-)  TODO float64 legacy
    b::MF = zeros(grid.Nx, grid.Ny)                 # Clapp-Hornberger exponent (-)
    hcap_soil::MF = zeros(grid.Nx, grid.Ny)         # Volumetric heat capacity of dry soil (J/K/m^3)
    hcon_soil::MF = zeros(grid.Nx, grid.Ny)         # Thermal conductivity of dry soil (W/m/K)
    sathh::MF = zeros(grid.Nx, grid.Ny)             # Saturated soil water pressure (m)
    Vsat::MF = zeros(grid.Nx, grid.Ny)              # Volumetric soil moisture at saturation (-)
    Vcrit::MF = zeros(grid.Nx, grid.Ny)             # Volumetric soil moisture at critical point (-)
end

@with_kw struct State{GT, MF <: AbstractMatrix{<:AbstractFloat}, MI <: AbstractMatrix{<:Integer}, AF <: AbstractArray{<:AbstractFloat, 3}}
    grid::GT
    albs::MF = 0.85 * ones(grid.Nx, grid.Ny)                 # Snow albedo (-)
    Ds::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)                 # Snow layer thicknesses (m)
    Nsnow::MI = zeros(Int, grid.Nx, grid.Ny)                      # Number of snow layers
    Qcan::MF = zeros(grid.Nx, grid.Ny)                           # Canopy air space humidity (kg/kg)
    Sice::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)               # Ice content of snow layers (kg/m^2)
    Sliq::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)               # Liquid content of snow layers (kg/m^2)
    Sveg::MF = zeros(grid.Nx, grid.Ny)                           # Snow mass on vegetation (kg/m^2)
    Tcan::MF = 285 * ones(grid.Nx, grid.Ny)                  # Canopy air space temperature (K)
    theta::AF = zeros(grid.Nsoil, grid.Nx, grid.Ny)              # Volumetric moisture content of soil layers (-)
    Tsnow::AF = 273.15 * ones(grid.Nsmax, grid.Nx, grid.Ny)  # Snow layer temperatures (K)
    Tsoil::AF = 285 * ones(grid.Nsoil, grid.Nx, grid.Ny)     # Soil layer temperatures (K)
    Tsrf::MF = 285 * ones(grid.Nx, grid.Ny)                  # Surface skin temperature (K)
    fsnow::MF = zeros(grid.Nx, grid.Ny)                          # Snow cover fraction (-)
    Tveg::MF = 285 * ones(grid.Nx, grid.Ny)                  # Vegetation temperature (K)
    snowdepthmin::MF = zeros(grid.Nx, grid.Ny)                   # Min snow depth at time of swemin (m)
    snowdepthmax::MF = zeros(grid.Nx, grid.Ny)                   # Max snow depth at time of swemax (m)
    snowdepthhist::AF = zeros(14, grid.Nx, grid.Ny)              # Snow depth over last 14 days (m)
    swemin::MF = zeros(grid.Nx, grid.Ny)                         # Minimum SWE during the season (kg/m^2)
    swemax::MF = zeros(grid.Nx, grid.Ny)                         # Maximum SWE during the season (kg/m^2)
    swehist::AF = zeros(14, grid.Nx, grid.Ny)                    # SWE over last 14 days (kg/m^2)
    histowet::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)           # Historical past wetting of a layer (-)
end

@with_kw struct Diagnostics{GT, MF <: AbstractMatrix{<:AbstractFloat}, AF <: AbstractArray{<:AbstractFloat, 3}}
    grid::GT
    # drive
    es::MF = zeros(grid.Nx, grid.Ny)                 # Saturation vapour pressure (Pa)
    Qa::MF = zeros(grid.Nx, grid.Ny)                 # Specific humidity (kg/kg)
    Uaeff::MF = zeros(grid.Nx, grid.Ny)              # Wind speed with lower bound applied (m/s)
    Sfeff::MF = zeros(grid.Nx, grid.Ny)              # Snowfall reaching the surface (kg/m^2/s)
    # radiation
    alb::MF = zeros(grid.Nx, grid.Ny)                # Albedo (-)
    asrf_out::MF = zeros(grid.Nx, grid.Ny)           # Surface albedo (-)
    SWveg::MF = zeros(grid.Nx, grid.Ny)              # Net shortwave absorbed by vegetation (W/m^2)
    SWsrf::MF = zeros(grid.Nx, grid.Ny)              # Net shortwave absorbed by the surface (W/m^2)
    SWsci::MF = zeros(grid.Nx, grid.Ny)              # Subcanopy incoming shortwave (W/m^2)
    LWeff::MF = zeros(grid.Nx, grid.Ny)              # Incoming longwave used in the energy balance (W/m^2)
    # thermal
    ksnow::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)  # Thermal conductivity of snow (W/m/K)
    csoil::AF = zeros(grid.Nsoil, grid.Nx, grid.Ny)  # Areal heat capacity of soil (J/K/m^2)
    ksoil::AF = zeros(grid.Nsoil, grid.Nx, grid.Ny)  # Thermal conductivity of soil (W/m/K)
    gs1::MF = zeros(grid.Nx, grid.Ny)                # Surface moisture conductance (m/s)
    Ds1::MF = zeros(grid.Nx, grid.Ny)                # Surface layer thickness (m)
    Ts1::MF = zeros(grid.Nx, grid.Ny)                # Surface layer temperature (K)
    ks1::MF = zeros(grid.Nx, grid.Ny)                # Surface thermal conductivity (W/m/K)
    Tveg0::MF = zeros(grid.Nx, grid.Ny)              # Vegetation temperature at start of timestep (K)
    # surface_exchange_coefficients
    KH::MF = zeros(grid.Nx, grid.Ny)                 # Eddy diffusivity for heat to the atmosphere (m/s)
    KHa::MF = zeros(grid.Nx, grid.Ny)                # Eddy diffusivity from the canopy air space (m/s)
    KHg::MF = zeros(grid.Nx, grid.Ny)                # Eddy diffusivity for heat from the ground (m/s)
    KHv::MF = zeros(grid.Nx, grid.Ny)                # Eddy diffusivity for heat from vegetation (m/s)
    KWg::MF = zeros(grid.Nx, grid.Ny)                # Eddy diffusivity for water from the ground (m/s)
    KWv::MF = zeros(grid.Nx, grid.Ny)                # Eddy diffusivity for water from vegetation (m/s)
    Usc::MF = zeros(grid.Nx, grid.Ny)                # Wind speed in canopy layer (m/s)
    # surface_energy_balance
    Esrf::MF = zeros(grid.Nx, grid.Ny)               # Moisture flux from the surface (kg/m^2/s)
    Eveg::MF = zeros(grid.Nx, grid.Ny)               # Moisture flux from vegetation (kg/m^2/s)
    G::MF = zeros(grid.Nx, grid.Ny)                  # Heat flux into the surface (W/m^2)
    H::MF = zeros(grid.Nx, grid.Ny)                  # Sensible heat flux to the atmosphere (W/m^2)
    Hsrf::MF = zeros(grid.Nx, grid.Ny)               # Sensible heat flux from the surface (W/m^2)
    LE::MF = zeros(grid.Nx, grid.Ny)                 # Latent heat flux to the atmosphere (W/m^2)
    LEsrf::MF = zeros(grid.Nx, grid.Ny)              # Latent heat flux from the surface (W/m^2)
    LWsci::MF = zeros(grid.Nx, grid.Ny)              # Subcanopy incoming longwave (W/m^2)
    LWveg::MF = zeros(grid.Nx, grid.Ny)              # Net longwave absorbed by vegetation (W/m^2)
    Melt::MF = zeros(grid.Nx, grid.Ny)               # Surface melt rate (kg/m^2/s)
    Rnet::MF = zeros(grid.Nx, grid.Ny)               # Net radiation (W/m^2)
    Rsrf::MF = zeros(grid.Nx, grid.Ny)               # Net radiation at surface (W/m^2)
    # canopy
    intcpt::MF = zeros(grid.Nx, grid.Ny)             # Canopy interception (kg/m^2)
    Sbveg::MF = zeros(grid.Nx, grid.Ny)              # Sublimation from vegetation (kg/m^2)
    unload::MF = zeros(grid.Nx, grid.Ny)             # Snow mass unloaded from canopy (kg/m^2)
    # snow
    Gsoil::MF = zeros(grid.Nx, grid.Ny)              # Heat flux into soil (W/m^2)
    Roff::MF = zeros(grid.Nx, grid.Ny)               # Total runoff (kg/m^2)
    meltflux_out::MF = zeros(grid.Nx, grid.Ny)       # Runoff from snowmelt at base of snow (kg/m^2)
    Sbsrf::MF = zeros(grid.Nx, grid.Ny)              # Sublimation from the snow surface (kg/m^2)
    Roff_bare::MF = zeros(grid.Nx, grid.Ny)          # Bare soil runoff (kg/m^2)
    Roff_snow::MF = zeros(grid.Nx, grid.Ny)          # Runoff at base of snow (kg/m^2)
    snowdepth0::MF = zeros(grid.Nx, grid.Ny)         # Snow depth at start of timestep (m)
    Sice0::MF = zeros(grid.Nx, grid.Ny)              # Ice content at start of timestep (kg/m^2)
    # snow_layering
    Ds0::MF = zeros(grid.Nx, grid.Ny)                # Snow layer thickness at start of timestep (m)
end

# Adapt (@adapt_structure) and on_architecture rebuild these structs by calling the
# type positionally with the converted fields. Tf is the declared type of no field -
# it appears only in the bounds on MF/MI/AF - so Julia generates no such constructor.
# Spell them out, taking each parameter from a representative field.
mutable struct FSM{Tf, Ti, G, P, L, S, D, PH}
    grid::G
    params::P
    landuse::L
    state::S
    diag::D
    physics::PH
end

# Positional constructor used by on_architecture
function FSM(grid::Grid, params::Parameters{Tf, Ti}, landuse::Landuse, state::State,
        diag::Diagnostics, physics) where {Tf, Ti}
    return FSM{Tf, Ti, typeof(grid), typeof(params), typeof(landuse), typeof(state),
        typeof(diag), typeof(physics)}(grid, params, landuse, state, diag, physics)
end

# No getproperty/setproperty! forwarding: fields are reached explicitly through
# their sub-struct (fsm.state.Ds, fsm.diag.KH, fsm.params.dt, ...). Each
# sub-struct is its own namespace, so a field name may repeat across them without
# ambiguity.

# FSM{Tf, Ti}(; Nx, Ny, schemes...) builds the model on plain CPU Arrays
function (::Type{FSM{Tf, Ti}})(;
        Nx = 1, Ny = 1,
        ALBEDO = PrognosticAlbedo{Tf}(Nx, Ny),
        CANOPY = NoCanopy{Tf}(),
        SUBSTR = SoilSubstrate{Tf}(),
        CONDCT = DensityConductivity{Tf}(),
        reference_height = AboveGround{Tf}(),
        surface_layer = OpenSurfaceLayer{Tf}(; stability = LouisStabilityCorrection{Tf}()),
        FSNRHO = ElevationFreshSnowDensity{Tf}(),
        COMPACT = CrocusCompaction{Tf}(),
        HYDROL = DensityBucketHydrology{Tf}(),
        LAYERING = OriginalLayering{Tf}(),
        SNFRAC = PointSnowFraction{Tf}()) where {Tf, Ti}

    grid    = Grid{Ti, Vector{Tf}}(; Nx = Nx, Ny = Ny)
    check_layer_thicknesses(grid)
    GT      = typeof(grid)
    params  = Parameters{Tf, Ti}()
    landuse = Landuse{GT, Matrix{Tf}, Matrix{Float64}}(; grid = grid)
    state   = State{GT, Matrix{Tf}, Matrix{Ti}, Array{Tf, 3}}(; grid = grid)
    diag    = Diagnostics{GT, Matrix{Tf}, Array{Tf, 3}}(; grid = grid)
    physics = (ALBEDO = ALBEDO, CANOPY = CANOPY, SUBSTR = SUBSTR, CONDCT = CONDCT,
        reference_height = reference_height, surface_layer = surface_layer,
        FSNRHO = FSNRHO, COMPACT = COMPACT, HYDROL = HYDROL, LAYERING = LAYERING,
        SNFRAC = SNFRAC)

    all(s -> s isa AbstractParameterization{Tf}, values(physics)) ||
        throw(ArgumentError("physics scheme precision does not match model Tf = $Tf"))

    return FSM(grid, params, landuse, state, diag, physics)
end

@with_kw mutable struct MET{
        Tf, Ti,
        MF <: AbstractMatrix{Tf}, MF64 <: AbstractMatrix{Float64},
        AF64_3 <: AbstractArray{Float64, 3},
    }

    # Domain size

    Nx::Ti = 1                                     # Size of first array dimension (rows)
    Ny::Ti = 1                                     # Size of second array dimension (columns)

    # Meteorological variables

    Sdir::MF = fill(NaN, Nx, Ny)                   # Direct shortwave radiation per inclined surface area (W/m^2)
    Sdif::MF = fill(NaN, Nx, Ny)                   # Diffuse shortwave radiation (W/m^2)
    Sdird::MF = fill(NaN, Nx, Ny)                  # Direct shortwave radiation per horizontal surface area (W/m^2)
    LW::MF = fill(NaN, Nx, Ny)                     # Incoming longwave radiation (W/m^2)
    Sf::MF = fill(NaN, Nx, Ny)                     # Snowfall rate (kg/m^2/s)
    Rf::MF = fill(NaN, Nx, Ny)                     # Rainfall rate (kg/m^2/s)
    Sf24h::MF = fill(NaN, Nx, Ny)                  # Total snowfall over 24h (kg/m^2)
    Ta::MF = fill(NaN, Nx, Ny)                     # Air temperature (K)
    RH::MF = fill(NaN, Nx, Ny)                     # Relative humidity (%)
    Ua::MF = fill(NaN, Nx, Ny)                     # Wind speed (m/s)
    Ps::MF = fill(NaN, Nx, Ny)                     # Surface air pressure (Pa)
    Tv::MF = fill(NaN, Nx, Ny)                     # Time-varying transmissivity for direct shortwave radiation (-)

    # Snowfall tracking variables

    Sf24h_f64::MF64 = zeros(Nx, Ny)                # Total snowfall over 24h (kg/m^2)  TODO intermediate variable using Float64 to match matlab/fortran code - remove later
    Sf_history_f64::AF64_3 = zeros(Nx, Ny, 24)     # History of snowfall over the last 24h (kg/m^2)  TODO using Float64 to match matlab/fortran code - change precision later

end

function (::Type{MET{Tf, Ti}})(; kwargs...) where {Tf, Ti}
    return MET{Tf, Ti, Matrix{Tf}, Matrix{Float64}, Array{Float64, 3}}(; kwargs...)
end

function MET(Nx::Integer, fields...)
    nt = NamedTuple{fieldnames(MET)}((Nx, fields...))
    return MET{eltype(nt.Sdir), typeof(nt.Nx), typeof(nt.Sdir), typeof(nt.Sf24h_f64),
        typeof(nt.Sf_history_f64)}(nt...)
end

# Let the array-holding structs cross into a kernel: Adapt rewrites each array
# field to the device array type at launch (a no-op on the CPU). Parameters and
# the physics schemes are isbits and need no adaptor.
@adapt_structure Grid
@adapt_structure Landuse
@adapt_structure State
@adapt_structure Diagnostics
@adapt_structure MET
