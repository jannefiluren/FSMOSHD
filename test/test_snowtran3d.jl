using FlexibleSnowModelOSHD
using Test
using Libdl

# Cross-validation of the Fortran (ccall) and pure-Julia SnowTran3D implementations: both run
# on identical synthetic inputs and their outputs must match within Float32 transport
# tolerances. The Fortran comparison is skipped (with a warning) when the shared library has
# not been built; the Julia routine still runs so the CPU path and the physical sanity checks
# are exercised on every machine.

const HAVE_LIBSNOWTRAN3D =
    isfile(joinpath(pkgdir(FlexibleSnowModelOSHD), "deps", "libsnowtran3d." * Libdl.dlext))

HAVE_LIBSNOWTRAN3D ||
    @warn "libsnowtran3d not built (run Pkg.build with gfortran) - skipping the Fortran vs Julia comparison"

"""
Set up a synthetic snow transport test case: spatially varying wind (some
pixels below the transport threshold of 4 m/s, others well above), cold dry
soft snow on most pixels, some wet/previously-wetted pixels (hard snow), and
some bare pixels. Layer densities lie within [rhos_min, rhos_max].

With `forest = true`, a forest fraction pattern is added: mostly open pixels,
some partly forested pixels that reduce the saltation flux, and dense forest
pixels (> 0.9) that block it. `tiled` sets the tiled trans run flag, which
changes how the forest weighting is applied.

Returns the model `fsm`, the meteo `met`, and the transport workspace `w`.
"""
function setup_snowtran3d_case(Udir0; forest = false, tiled = false)

    Nx, Ny = 12, 10

    fsm = FSM{Float32, Int32}(Nx = Nx, Ny = Ny)
    met = MET{Float32, Int32}(Nx = Nx, Ny = Ny)
    w = SnowTransport{Float32, Int32}(Nx = Int32(Nx), Ny = Int32(Ny))
    w.tiled_trans_run = tiled

    state = fsm.state

    for j in 1:Ny, i in 1:Nx

        if forest
            ff = (3 * i + 5 * j) % 10
            w.forestfrac[i, j] = ff < 6 ? 0.0f0 : (ff < 9 ? 0.3f0 : 0.95f0)
        end

        # Wind field: speed varies over the domain (2-16 m/s; suspension needs
        # friction velocities above ~0.7 m/s), direction wiggles around the
        # base direction Udir0
        met.Ua[i, j] = 2.0f0 + 14.0f0 * Float32(i + j) / Float32(Nx + Ny)
        met.Udir[i, j] = Udir0 + 5.0f0 * ((i + 2 * j) % 5 - 2)

        # Cold and fairly dry air
        met.Ta[i, j] = 260.0f0 + i % 5
        met.RH[i, j] = 60.0f0 + 3.0f0 * (j % 10)

        # Terrain / surface properties
        fsm.surface.Ld[i, j] = 100.0f0                    # Grid cell size (m)
        w.vegsnowd_xy[i, j] = 0.1f0
        fsm.surface.z0_snow[i, j] = 0.001f0 + 0.001f0 * (j % 5)

        # Snowpack with 0-2 layers (bare pixels included)
        nsnow = (i + 2 * j) % 3
        state.Nsnow[i, j] = nsnow
        state.fsnow[i, j] = nsnow == 0 ? 0.0f0 : (isodd(i + j) ? 0.8f0 : 1.0f0)

        # Some pixels have a wet / previously wetted (hard) snowpack
        wet = (i + j) % 4 == 0

        for k in 1:nsnow
            rho = 150.0f0 + 50.0f0 * k + 10.0f0 * (i % 3)  # Layer density (kg/m^3)
            state.Ds[k, i, j] = 0.3f0 + 0.05f0 * k + 0.01f0 * (j % 4)
            if wet
                state.Sice[k, i, j] = 0.95f0 * rho * state.Ds[k, i, j] * state.fsnow[i, j]
                state.Sliq[k, i, j] = 0.05f0 * rho * state.Ds[k, i, j] * state.fsnow[i, j]
                state.histowet[k, i, j] = 1.0f0
            else
                state.Sice[k, i, j] = rho * state.Ds[k, i, j] * state.fsnow[i, j]
                state.Sliq[k, i, j] = 0.0f0
                state.histowet[k, i, j] = 0.0f0
            end
            state.Tsnow[k, i, j] = 262.0f0 + k
        end

    end

    # The transport routines read the wind from fsm.diag.Uaeff, which drive!
    # computes from met.Ua; the tests call the routines directly, so mirror drive! here
    fsm.diag.Uaeff .= max.(met.Ua, 0.1f0)

    return fsm, met, w

end

total_swe_tran(fsm, Sice0) = sum(fsm.state.Sice) + sum(fsm.state.Sliq) + sum(Sice0)

# Base wind directions: westerly, easterly, southerly, northerly, diagonal,
# plus a swirling field with sign changes in both wind components; each with
# and without forest, and the forest cases also as a tiled run
@testset "Fortran vs Julia comparison (Udir = $Udir0, forest = $forest, tiled = $tiled)" for
    Udir0 in (270.0f0, 90.0f0, 180.0f0, 0.0f0, 225.0f0, :swirl),
        (forest, tiled) in ((false, false), (true, false), (true, true))

    fsm_f, met, w_f = setup_snowtran3d_case(Udir0 === :swirl ? 270.0f0 : Udir0; forest, tiled)
    Nx, Ny = fsm_f.grid.Nx, fsm_f.grid.Ny
    if Udir0 === :swirl
        for j in 1:Ny, i in 1:Nx
            met.Udir[i, j] = Float32((i * 37 + j * 61) % 360)
        end
    end
    fsm_j, w_j = deepcopy(fsm_f), deepcopy(w_f)

    swe_before = total_swe_tran(fsm_f, zeros(Float32, Nx, Ny))

    # Two consecutive calls: the second starts from the wind-compacted and
    # eroded snowpack left by the first
    for call in 1:2

        # The caller zeroes these before the call
        snowdepth0_f = zeros(Float32, Nx, Ny)
        Sice0_f = zeros(Float32, Nx, Ny)
        dSWE_salt_f = zeros(Float32, Nx, Ny)
        dSWE_susp_f = zeros(Float32, Nx, Ny)
        dSWE_subl_f = zeros(Float32, Nx, Ny)

        snowdepth0_j = zeros(Float32, Nx, Ny)
        Sice0_j = zeros(Float32, Nx, Ny)
        dSWE_salt_j = zeros(Float32, Nx, Ny)
        dSWE_susp_j = zeros(Float32, Nx, Ny)
        dSWE_subl_j = zeros(Float32, Nx, Ny)

        snowtran3d_julia!(fsm_j, met, w_j, snowdepth0_j, Sice0_j, dSWE_salt_j, dSWE_susp_j, dSWE_subl_j)

        if HAVE_LIBSNOWTRAN3D
            snowtran3d!(fsm_f, met, w_f, snowdepth0_f, Sice0_f, dSWE_salt_f, dSWE_susp_f, dSWE_subl_f)

            @testset "State comparison after call $call" begin
                atol = 1.0f-5
                @test isapprox(snowdepth0_j, snowdepth0_f; rtol = 1.0f-4, atol = atol)
                @test isapprox(Sice0_j, Sice0_f; rtol = 1.0f-4, atol = atol)
                @test isapprox(dSWE_salt_j, dSWE_salt_f; rtol = 1.0f-4, atol = atol)
                @test isapprox(dSWE_susp_j, dSWE_susp_f; rtol = 1.0f-4, atol = atol)
                @test isapprox(dSWE_subl_j, dSWE_subl_f; rtol = 1.0f-4, atol = atol)
                @test isapprox(w_j.dSWE_tot_salt, w_f.dSWE_tot_salt; rtol = 1.0f-4, atol = atol)
                @test isapprox(w_j.dSWE_tot_susp, w_f.dSWE_tot_susp; rtol = 1.0f-4, atol = atol)
                @test isapprox(w_j.dSWE_tot_subl, w_f.dSWE_tot_subl; rtol = 1.0f-4, atol = atol)
                @test isapprox(fsm_j.state.Ds, fsm_f.state.Ds; rtol = 1.0f-4, atol = atol)
                @test isapprox(fsm_j.state.Sice, fsm_f.state.Sice; rtol = 1.0f-4, atol = atol)
                @test isapprox(fsm_j.state.Sliq, fsm_f.state.Sliq; rtol = 1.0f-4, atol = atol)
                @test isapprox(fsm_j.state.Tsnow, fsm_f.state.Tsnow; rtol = 1.0f-4, atol = atol)
                @test fsm_j.state.histowet == fsm_f.state.histowet
                @test fsm_j.state.Nsnow == fsm_f.state.Nsnow
            end
        end

        if call == 1
            @testset "Sanity checks" begin

                # Snow was actually transported (the comparison is not trivial)
                @test any(!=(0), dSWE_salt_j)
                @test any(!=(0), dSWE_susp_j)
                @test any(>(0), Sice0_j)

                # Sublimation is a net loss
                @test sum(dSWE_subl_j) <= 0

                # Mass balance: the change in total SWE equals the net SWE
                # change from saltation, suspension and sublimation (transport
                # over the downwind domain edge leaves the domain). In a tiled
                # run the transferred SWE is intentionally re-weighted for the
                # later tile combination, so mass is not conserved within the run.
                if !tiled
                    swe_after = total_swe_tran(fsm_j, Sice0_j)
                    dswe_net = sum(dSWE_salt_j) + sum(dSWE_susp_j) + sum(dSWE_subl_j)
                    @test isapprox(swe_after - swe_before, dswe_net; rtol = 1.0f-3, atol = 0.1f0)
                end

            end
        end

    end

end

@testset "No transport below wind_min" begin

    # Wind everywhere below wind_min = 4 m/s: no transport, and (since the
    # upstream code dropped the surface_snow wind compaction) the snowpack
    # state stays untouched
    fsm_f, met, w_f = setup_snowtran3d_case(270.0f0)
    met.Ua .= 3.0f0
    fsm_f.diag.Uaeff .= 3.0f0
    fsm_j, w_j = deepcopy(fsm_f), deepcopy(w_f)
    Ds_before = copy(fsm_f.state.Ds)

    Nx, Ny = fsm_f.grid.Nx, fsm_f.grid.Ny
    args_j = [zeros(Float32, Nx, Ny) for _ in 1:5]

    snowtran3d_julia!(fsm_j, met, w_j, args_j...)

    # The snowpack is left unchanged and nothing was transported
    @test fsm_j.state.Ds == Ds_before
    @test all(a == zeros(Float32, Nx, Ny) for a in args_j)

    if HAVE_LIBSNOWTRAN3D
        args_f = [zeros(Float32, Nx, Ny) for _ in 1:5]
        snowtran3d!(fsm_f, met, w_f, args_f...)
        @test isapprox(fsm_j.state.Ds, fsm_f.state.Ds; rtol = 1.0f-4, atol = 1.0f-5)
    end

end

@testset "Erosion down to the holding-depth floor (swe_from_hs loss path)" begin

    # Shallow, dry, fully-soft snow under a strong wind: erosion drives cells down to the
    # vegetation snow-holding floor, exercising the getnewdepth / sublimation loss branches that
    # convert an eroded depth back to SWE via swe_from_hs. Regression guard - this path is not
    # reached by the moderate-wind cases above (it surfaced only on a full operational domain).
    fsm_f, met, w_f = setup_snowtran3d_case(270.0f0)
    Nx, Ny = fsm_f.grid.Nx, fsm_f.grid.Ny

    met.Ua .*= 1.8f0                              # up to ~29 m/s, with the existing gradient
    fsm_f.diag.Uaeff .= max.(met.Ua, 0.1f0)

    st = fsm_f.state
    st.Ds .= 0.0f0; st.Sice .= 0.0f0; st.Sliq .= 0.0f0; st.histowet .= 0.0f0; st.Tsnow .= 263.0f0
    for j in 1:Ny, i in 1:Nx
        st.Nsnow[i, j] = 1
        st.fsnow[i, j] = 1.0f0
        st.Ds[1, i, j] = 0.15f0
        st.Sice[1, i, j] = 120.0f0 * 0.15f0      # dry, low-density (soft) snow
    end

    fsm_j, w_j = deepcopy(fsm_f), deepcopy(w_f)
    aj = [zeros(Float32, Nx, Ny) for _ in 1:5]
    snowtran3d_julia!(fsm_j, met, w_j, aj...)

    # Some cells eroded below the initial 0.15 m (guards against the swe_from_hs MethodError)
    @test any(<(0.15f0 - 1.0f-6), fsm_j.state.Ds[1, :, :])

    if HAVE_LIBSNOWTRAN3D
        af = [zeros(Float32, Nx, Ny) for _ in 1:5]
        snowtran3d!(fsm_f, met, w_f, af...)
        @test isapprox(fsm_j.state.Ds, fsm_f.state.Ds; rtol = 1.0f-4, atol = 1.0f-5)
        @test isapprox(fsm_j.state.Sice, fsm_f.state.Sice; rtol = 1.0f-4, atol = 1.0f-5)
        for (a, b) in zip(aj, af)
            @test isapprox(a, b; rtol = 1.0f-4, atol = 1.0f-5)
        end
    end

end
