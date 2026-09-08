using FlexibleSnowModelOSHD
using Test
using Libdl

# Cross-validation of the Fortran (ccall) and pure-Julia SnowSlide implementations: both run
# on identical synthetic inputs and their outputs must match. The Fortran comparison is
# skipped (with a warning) when the shared library has not been built (Pkg.build needs
# gfortran); the Julia routine still runs so the CPU path and the physical sanity checks are
# exercised on every machine.

const HAVE_LIBSNOWSLIDE =
    isfile(joinpath(pkgdir(FlexibleSnowModelOSHD), "deps", "libsnowslide." * Libdl.dlext))

HAVE_LIBSNOWSLIDE ||
    @warn "libsnowslide not built (run Pkg.build with gfortran) - skipping the Fortran vs Julia comparison"

"""
Set up a synthetic snow slide test case on a steep inclined plane with a
cross-slope tilt and a flat run-out zone in the lowest rows. Snow depths are
well above the snow holding depth so that plenty of snow slides, edge pixels
are steep to exercise the edge-dumping branch, deposits trigger the dynamic
reduction of the snow holding depth, and bare pixels in the steep zone receive
deposits that partly re-slide without an underlying snowpack.

With `forest = true`, a forest fraction pattern is added: mostly open pixels,
some partly forested pixels that reduce the slide transfer, pixels above the
50% no-slide threshold, and dense forest. `tiled` sets the tiled trans run
flag, which changes how the forest weighting is applied.

Returns the model `fsm` and its transport workspace `w`.
"""
function setup_snowslide_case(; forest = false, tiled = false)

    Nx, Ny = 12, 10

    fsm = FSM{Float32, Int32}(Nx = Nx, Ny = Ny)
    w = SnowTransport{Float32, Int32}(Nx = Int32(Nx), Ny = Int32(Ny))
    w.tiled_trans_run = tiled

    state = fsm.state

    for j in 1:Ny, i in 1:Nx

        if forest
            ff = (2 * i + 3 * j) % 8
            w.forestfrac[i, j] =
                ff < 5 ? 0.0f0 : (ff == 5 ? 0.3f0 : (ff == 6 ? 0.6f0 : 0.95f0))
        end

        # Terrain descending from north (high i) to south (low i), with a small
        # cross-slope tilt so that diagonal neighbours also receive snow
        fsm.surface.dem[i, j] = 1000.0f0 + 50.0f0 * i + 2.0f0 * j

        # Steep in the upper part, flat run-out zone in the lowest rows
        w.slope[i, j] = i >= 4 ? 45.0f0 : 10.0f0

        # Small snow holding depth so that a lot of the snow slides
        w.Shd[i, j] = 0.05f0 + 0.01f0 * ((i + j) % 5)

        # Snowpack with 0-2 layers (bare pixels included) and layer densities
        # within [rhos_min, rhos_max]
        nsnow = (i + 2 * j) % 3
        state.Nsnow[i, j] = nsnow
        state.fsnow[i, j] = nsnow == 0 ? 0.0f0 : (isodd(i + j) ? 0.8f0 : 1.0f0)
        for k in 1:nsnow
            rho = 150.0f0 + 50.0f0 * k + 10.0f0 * (i % 3)  # Layer density (kg/m^3)
            state.Ds[k, i, j] = 0.2f0 + 0.05f0 * k + 0.01f0 * (j % 4)
            state.Sice[k, i, j] = 0.95f0 * rho * state.Ds[k, i, j] * state.fsnow[i, j]
            state.Sliq[k, i, j] = 0.05f0 * rho * state.Ds[k, i, j] * state.fsnow[i, j]
            state.Tsnow[k, i, j] = 265.0f0 + k
            state.histowet[k, i, j] = Float32((i + j + k) % 2)
        end

    end

    # Processing order from highest to lowest pixel
    FlexibleSnowModelOSHD.sort_dem_indices!(w.index_sorted_dem, fsm.surface.dem)

    return fsm, w

end

total_swe(fsm, Sice0) = sum(fsm.state.Sice) + sum(fsm.state.Sliq) + sum(Sice0)

@testset "Fortran vs Julia comparison (forest = $forest, tiled = $tiled)" for (forest, tiled) in
    ((false, false), (true, false), (true, true))

    # Identical initial states for the two implementations
    fsm_f, w_f = setup_snowslide_case(; forest, tiled)
    fsm_j, w_j = deepcopy(fsm_f), deepcopy(w_f)

    Nx, Ny = fsm_f.grid.Nx, fsm_f.grid.Ny

    snowdepth0_f = zeros(Float32, Nx, Ny)
    Sice0_f = zeros(Float32, Nx, Ny)
    dSWE_slide_f = zeros(Float32, Nx, Ny)

    snowdepth0_j = zeros(Float32, Nx, Ny)
    Sice0_j = zeros(Float32, Nx, Ny)
    dSWE_slide_j = zeros(Float32, Nx, Ny)

    swe_before = total_swe(fsm_j, Sice0_j)
    swe_injected = 0.0f0

    # Two calls: the first slides snow from the original snowpack, the second
    # starts from a fresh avalanche deposit injected on entry (as if deposited
    # by a slide from outside the local neighbourhood) and re-slides it
    for call in 1:2

        if call == 2
            for j in 4:6, i in 7:9
                for (snowdepth0, Sice0) in ((snowdepth0_f, Sice0_f), (snowdepth0_j, Sice0_j))
                    snowdepth0[i, j] += 1.5f0
                    Sice0[i, j] += 1.5f0 * w_f.rho_deposit
                end
                swe_injected += 1.5f0 * w_f.rho_deposit
            end
        end

        snowslide_julia!(fsm_j, w_j, snowdepth0_j, Sice0_j, dSWE_slide_j)

        if HAVE_LIBSNOWSLIDE
            snowslide!(fsm_f, w_f, snowdepth0_f, Sice0_f, dSWE_slide_f)

            @testset "State comparison after call $call" begin
                @test snowdepth0_j ≈ snowdepth0_f rtol = 1.0f-5
                @test Sice0_j ≈ Sice0_f rtol = 1.0f-5
                @test dSWE_slide_j ≈ dSWE_slide_f rtol = 1.0f-5
                @test w_j.dSWE_tot_slide ≈ w_f.dSWE_tot_slide rtol = 1.0f-5
                @test fsm_j.state.Ds ≈ fsm_f.state.Ds rtol = 1.0f-5
                @test fsm_j.state.Sice ≈ fsm_f.state.Sice rtol = 1.0f-5
                @test fsm_j.state.Sliq ≈ fsm_f.state.Sliq rtol = 1.0f-5
                @test fsm_j.state.Tsnow ≈ fsm_f.state.Tsnow rtol = 1.0f-5
                @test fsm_j.state.histowet == fsm_f.state.histowet
                @test fsm_j.state.Nsnow == fsm_f.state.Nsnow
            end
        end

    end

    @testset "Sanity checks" begin

        # Snow was actually redistributed (the comparison is not trivial)
        @test any(!=(0), dSWE_slide_j)
        @test any(>(0), Sice0_j)

        # Mass balance: the change in total SWE (apart from the injected
        # deposit) equals the net dSWE_slide, which is negative because edge
        # pixels dump snow out of the domain. In a tiled run the transferred
        # SWE is intentionally re-weighted for the later tile combination, so
        # mass is not conserved within the run.
        if !tiled
            swe_after = total_swe(fsm_j, Sice0_j)
            @test swe_after - swe_before - swe_injected ≈ sum(w_j.dSWE_tot_slide) rtol = 1.0f-4
            @test sum(w_j.dSWE_tot_slide) < 0
        end

        # Snow moves downhill: deposits only on pixels with a higher neighbour
        dem = fsm_j.surface.dem
        for j in 1:Ny, i in 1:Nx
            if Sice0_j[i, j] > 0
                @test any(
                    dem[i + di, j + dj] > dem[i, j]
                        for di in -1:1, dj in -1:1
                        if checkbounds(Bool, dem, i + di, j + dj)
                )
            end
        end

    end

end
