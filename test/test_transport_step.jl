using FlexibleSnowModelOSHD
using Test
using Dates

# transport! is the step!-level snow-transport step. These tests check three things on a small
# multi-cell open-tile domain: (1) transport! reproduces exactly the documented
# operator-then-relayer sequence, (2) step! runs transport! when a workspace is passed and is a
# no-op otherwise, (3) transport's relayer passes do not roll the 14-day history (update_hist is
# false there). The Julia operators are used throughout so the tests need no Fortran build; the
# operators themselves are cross-validated against Fortran in test_snowslide/test_snowtran3d.

const Nx, Ny = 8, 6
const t_test = DateTime(2024, 1, 15, 12)      # not a history-roll hour

"""
Build a small open-tile model with a seeded snowpack, strong wind and a sloped DEM, plus its
transport workspace. `wind`/`slide` select the processes. Returns `(fsm, met, w)`.
"""
function make_transport_step_case(; wind = false, slide = false)

    elevation = Float64[1000 + 60 * i + 3 * j for i in 1:Nx, j in 1:Ny]
    slope = Float64[i >= 3 ? 50.0 : 8.0 for i in 1:Nx, j in 1:Ny]

    lus = Dict(
        "elevation" => Dict("data" => elevation),
        "slope" => Dict("data" => slope),
        "forest" => Dict("data" => zeros(Nx, Ny)),
        "skyvf" => Dict("data" => ones(Nx, Ny)),
        "slopemu" => Dict("data" => ones(Nx, Ny)),
        "xi" => Dict("data" => ones(Nx, Ny)),
        "Ld" => Dict("data" => fill(100.0, Nx, Ny)),
        "prec_multi" => Dict("data" => ones(Nx, Ny)),
    )

    settings = Dict("tile" => "open", "config" => Dict("SNFRAC" => 0))
    fsm = setup(Float32, lus, Nx, Ny, settings)
    w = setup_transport(fsm, lus; wind = wind, slide = slide, use_fortran = false)

    # Seed a snowpack (0-2 layers per cell, densities within [rhos_min, rhos_max])
    state = fsm.state
    for j in 1:Ny, i in 1:Nx
        nsnow = (i + 2 * j) % 3
        state.Nsnow[i, j] = nsnow
        state.fsnow[i, j] = nsnow == 0 ? 0.0f0 : (isodd(i + j) ? 0.8f0 : 1.0f0)
        for k in 1:nsnow
            rho = 150.0f0 + 50.0f0 * k + 10.0f0 * (i % 3)
            state.Ds[k, i, j] = 0.6f0 + 0.3f0 * k
            state.Sice[k, i, j] = rho * state.Ds[k, i, j] * state.fsnow[i, j]
            state.Sliq[k, i, j] = 0.0f0
            state.Tsnow[k, i, j] = 263.0f0 + k
            state.histowet[k, i, j] = 0.0f0
        end
    end

    # Full, benign forcing so step! runs cleanly; strong cold wind so transport fires
    met = MET{Float32}(Nx = Nx, Ny = Ny)
    met.Sdir .= 0.0f0
    met.Sdif .= 0.0f0
    met.Sdird .= 0.0f0
    met.LW .= 250.0f0
    met.Rf .= 0.0f0
    met.Sf .= 0.0f0
    met.Sf24h .= 0.0f0
    met.Ta .= 263.0f0
    met.RH .= 80.0f0
    met.Ua .= 10.0f0
    met.Ps .= 80000.0f0
    met.Tv .= 1.0f0
    met.Udir .= 270.0f0
    fsm.diag.Uaeff .= max.(met.Ua, 0.1f0)

    return fsm, met, w
end

state_fields(fsm) = (
    fsm.state.Ds, fsm.state.Sice, fsm.state.Sliq, fsm.state.Tsnow,
    fsm.state.Nsnow, fsm.state.fsnow, fsm.state.histowet,
)

@testset "transport! reproduces the operator + relayer sequence" begin

    @testset "wind" begin
        fsm, met, w = make_transport_step_case(; wind = true)

        fsm_a, w_a = deepcopy(fsm), deepcopy(w)
        transport!(fsm_a, met, w_a, t_test)

        # Manual: the exact documented sequence
        fsm_b, w_b = deepcopy(fsm), deepcopy(w)
        d = fsm_b.diag
        fill!(d.snowdepth0, 0.0f0); fill!(d.Sice0, 0.0f0)
        fill!(w_b.dSWE_salt, 0.0f0); fill!(w_b.dSWE_susp, 0.0f0); fill!(w_b.dSWE_subl, 0.0f0)
        snowtran3d_julia!(fsm_b, met, w_b, d.snowdepth0, d.Sice0, w_b.dSWE_salt, w_b.dSWE_susp, w_b.dSWE_subl)
        relayer!(fsm_b, met, t_test; update_hist = false)

        for (a, b) in zip(state_fields(fsm_a), state_fields(fsm_b))
            @test a == b
        end
        @test w_a.dSWE_tot_salt == w_b.dSWE_tot_salt
        @test w_a.dSWE_tot_susp == w_b.dSWE_tot_susp
        @test w_a.dSWE_tot_subl == w_b.dSWE_tot_subl

        # The comparison is non-trivial: transport actually moved snow
        @test any(!=(0), w_a.dSWE_tot_salt)
    end

    @testset "slide" begin
        fsm, met, w = make_transport_step_case(; slide = true)

        fsm_a, w_a = deepcopy(fsm), deepcopy(w)
        transport!(fsm_a, met, w_a, t_test)

        fsm_b, w_b = deepcopy(fsm), deepcopy(w)
        d = fsm_b.diag
        fill!(d.snowdepth0, 0.0f0); fill!(d.Sice0, 0.0f0)
        fill!(w_b.dSWE_slide, 0.0f0)
        snowslide_julia!(fsm_b, w_b, d.snowdepth0, d.Sice0, w_b.dSWE_slide)
        relayer!(fsm_b, met, t_test; update_hist = false)

        for (a, b) in zip(state_fields(fsm_a), state_fields(fsm_b))
            @test a == b
        end
        @test w_a.dSWE_tot_slide == w_b.dSWE_tot_slide
        @test any(!=(0), w_a.dSWE_tot_slide)
    end

end

@testset "step! runs transport! only when a workspace is passed" begin

    fsm, met, w = make_transport_step_case(; wind = true, slide = true)

    fsm_none = deepcopy(fsm)
    step!(fsm_none, met, t_test)                              # no transport

    fsm_tr = deepcopy(fsm)
    step!(fsm_tr, met, t_test; transport = deepcopy(w))       # with transport

    # Transport redistributes snow, so the two snowpacks must differ
    @test fsm_tr.state.Ds != fsm_none.state.Ds

    # And a no-transport step! is exactly a step! with transport = nothing (default), i.e. the
    # kwarg default changes nothing
    fsm_default = deepcopy(fsm)
    step!(fsm_default, met, t_test; transport = nothing)
    @test fsm_default.state.Ds == fsm_none.state.Ds
    @test fsm_default.state.Sice == fsm_none.state.Sice

end

@testset "transport does not roll the 14-day history" begin

    # At a history-roll hour, snow!'s new-snow pass rolls the 14-day buffers once; transport!'s
    # relayer passes must not roll them again (update_hist = false)
    t_roll = DateTime(2024, 1, 15, 5)
    fsm, met, w = make_transport_step_case(; wind = true, slide = true)

    swehist_before = copy(fsm.state.swehist)
    snowdepthhist_before = copy(fsm.state.snowdepthhist)

    transport!(fsm, met, w, t_roll)

    @test fsm.state.swehist == swehist_before
    @test fsm.state.snowdepthhist == snowdepthhist_before

end
