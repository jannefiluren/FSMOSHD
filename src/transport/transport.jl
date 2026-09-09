# The snow-transport step of the model, called from step! after snow!. Transport is a
# neighbour-coupled, grid-global, CPU-only operator (SnowTran3D wind transport + SnowSlide),
# kept separate from the per-cell snow! kernel. Each enabled process writes its redistributed
# ("incoming") mass into the diag deposit arrays snowdepth0/Sice0 and erodes the pack directly;
# a relayer! pass then accumulates the deposit into the snowpack - the same three-step idiom the
# pre-refactor snow! used, but as its own step driven by a run-level SnowTransport workspace.

# Fall back to the Julia port (with a warning) when use_fortran is set but the shared library
# was not built, so a run never hard-fails on a missing gfortran toolchain.
function transport_use_fortran(w::SnowTransport, lib::String)
    w.use_fortran || return false
    if !isfile(joinpath(@__DIR__, "..", "..", "deps", lib * "." * Libdl.dlext))
        @warn "use_fortran = true but $lib is not built (run Pkg.build with gfortran); " *
            "falling back to the Julia implementation" maxlog = 1
        return false
    end
    return true
end

"""
    transport!(fsm, met, w::SnowTransport, t)

Redistribute snow horizontally for one time step: wind transport (SnowTran3D) and/or snow
slides (SnowSlide), as selected by `w.wind` / `w.slide`, using the Fortran or Julia
implementation per `w.use_fortran`. Each process zeroes the deposit arrays
`fsm.diag.snowdepth0` / `Sice0`, runs the operator (which erodes the pack and writes incoming
mass into the deposit), then relayers the deposit with [`relayer!`](@ref).

Transport is **CPU-only**; `fsm` must hold host `Array`s. Called from `step!` when a workspace
is passed; a normal `step!` (no `transport`) does not run it.
"""
function transport!(fsm::FSM{Tf}, met::MET{Tf}, w::SnowTransport{Tf}, t) where {Tf}

    fsm.state.Ds isa Array ||
        throw(ArgumentError("snow transport is CPU-only; `fsm` must hold host Arrays (got $(typeof(fsm.state.Ds)))"))

    (w.wind || w.slide) || return nothing

    (; snowdepth0, Sice0) = fsm.diag

    if w.wind
        fill!(snowdepth0, zero(Tf))
        fill!(Sice0, zero(Tf))
        fill!(w.dSWE_salt, zero(Tf))
        fill!(w.dSWE_susp, zero(Tf))
        fill!(w.dSWE_subl, zero(Tf))
        if transport_use_fortran(w, "libsnowtran3d")
            snowtran3d!(fsm, met, w, snowdepth0, Sice0, w.dSWE_salt, w.dSWE_susp, w.dSWE_subl)
        else
            snowtran3d_julia!(fsm, met, w, snowdepth0, Sice0, w.dSWE_salt, w.dSWE_susp, w.dSWE_subl)
        end
        relayer!(fsm, met, t; update_hist = false)
    end

    if w.slide
        fill!(snowdepth0, zero(Tf))
        fill!(Sice0, zero(Tf))
        fill!(w.dSWE_slide, zero(Tf))
        if transport_use_fortran(w, "libsnowslide")
            snowslide!(fsm, w, snowdepth0, Sice0, w.dSWE_slide)
        else
            snowslide_julia!(fsm, w, snowdepth0, Sice0, w.dSWE_slide)
        end
        relayer!(fsm, met, t; update_hist = false)
    end

    return nothing
end
