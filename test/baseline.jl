#
# Output baseline manifest — generate one, or check the current tree against a frozen one.
#
#   julia --project=. test/baseline.jl generate out.txt
#   julia --project=. test/baseline.jl check            # against test/baseline_6c4dda5.txt
#   julia --project=. test/baseline.jl check other.txt
#   julia --project=. test/baseline.jl matrix out.txt   # per-flag matrix (slow, ~30 min)
#
# The manifest is 3 tile configurations x 19 output variables. `hash` is the authority;
# min/max/mean are recorded so that a differing hash can be triaged by eye.
#
# Why this exists alongside test/test_regression.jl: that test compares against
# test/simulation_results.jls, which is *gitignored* and regenerated whenever it is absent.
# A refactor can therefore be validated against a reference produced by the refactor itself.
# The manifest here is committed, so the baseline cannot drift.
#
# NOTE: the helper functions (load_domain_data, interpolate_meteo, run_simulations) live at the
# top of test_regression.jl, above its testsets. They are read out by slicing the file rather
# than `include`ing it, which would run the full regression suite. Moving them into their own
# file would be cleaner and is worth doing at some point.
#

using FlexibleSnowModelOSHD
using NCDatasets
using Dates
using Printf

const PROJDIR = pkgdir(FlexibleSnowModelOSHD)
const DEFAULT_BASELINE = joinpath(PROJDIR, "test", "baseline_6c4dda5.txt")

let src = read(joinpath(PROJDIR, "test", "test_regression.jl"), String)
    marker = findfirst("# Test data paths", src)
    marker === nothing && error("could not locate the helper/testset boundary in test_regression.jl")
    include_string(Main, src[1:first(marker) - 1], "test_regression_helpers.jl")
end

# The three tile configurations, matching test/test_regression.jl
const TILE_SETTINGS = [
    ("open", Dict("tile" => "open", "config" => Dict("SNFRAC" => 0))),
    (
        "forest", Dict(
            "tile" => "forest",
            "config" => Dict("CANMOD" => 1, "EXCHNG" => 2, "SNFRAC" => 4, "ZOFFST" => 1),
            "params" => Dict("hfsn" => 0.3, "z0_snow" => 0.01),
        ),
    ),
    ("glacier", Dict("tile" => "glacier", "config" => Dict("SNFRAC" => 0))),
]

# Per-flag baseline matrix. Every configuration flag that the type-domain refactor will convert,
# at every value it accepts, run over all three tiles. A converted parameterization must reproduce
# the rows for its old integer value bit-for-bit.
#
# Generated from a *pinned* commit (6c4dda5) rather than from the tree preceding each conversion,
# so the target cannot drift stage by stage.
#
# SNTRAN / SNSLID are absent: enabling them needs a 'slope' field that data/domain_data.nc does
# not carry, so transport cannot be baselined with the current fixture at all.
#
# NOTE the `Vector{Any}` values: written as a plain literal, Julia promotes the element type across
# all the pairs, so `[false, true]` silently becomes `[0, 1]`. That matters because `setfield!` on a
# `Bool` field *rejects* an integer ("non-boolean (Int32) used in boolean context") rather than
# converting it, so HN_ON must be fed genuine Bools.
const FLAG_MATRIX = Pair{String, Vector{Any}}[
    "ALBEDO" => [0, 1, 2],
    "CONDCT" => [0, 1],
    "DENSTY" => [1, 2, 3],   # constant-density (0) option removed in the AbstractCompaction refactor
    "HYDROL" => [0, 1, 2],
    "SNFRAC" => [0, 1, 2, 3, 4],
    "EXCHNG" => [0, 1, 2],
    "ZOFFST" => [0, 1],
    "FSNRHO" => [0, 1, 2],
    "SNOLAY" => [0, 1],
    "HN_ON" => [false, true],
]

"""
    generate_matrix(path)

Write a per-flag baseline matrix to `path`: each flag in [`FLAG_MATRIX`] at each of its values,
over all three tiles. A configuration that raises is recorded as an ERROR row rather than
aborting the run, since an unsupported combination is itself worth knowing about.
"""
function generate_matrix(path, only_flags = String[])

    selected = isempty(only_flags) ? FLAG_MATRIX :
        [pr for pr in FLAG_MATRIX if first(pr) in only_flags]
    isempty(selected) && error("no flags matched $only_flags")

    commit = try
        strip(read(`git -C $PROJDIR rev-parse HEAD`, String))
    catch
        "unknown"
    end

    total = sum(length(vals) for (_, vals) in selected) * length(TILE_SETTINGS)
    done = 0

    open(path, "w") do io
        println(io, "# FSM per-flag baseline matrix")
        println(io, "# source commit : ", commit)
        println(io, "# hash is the authority; min/max/mean are for triage. See baseline.jl.")
        println(io, "#")
        @printf(io, "# %-16s %-8s %-10s %-18s %8s %14s %14s %14s\n",
            "config", "tile", "variable", "hash", "nonfinite", "min", "max", "mean")

        for (flag, values) in selected, value in values
            label = string(flag, "=", value)
            for (name, base) in TILE_SETTINGS
                done += 1
                settings = deepcopy(base)
                # `Dict("SNFRAC" => 0)` is a Dict{String,Int64}; assigning a Bool into it would
                # convert to 0/1, and setfield! then rejects an integer for a Bool field.
                settings["config"] = Dict{String, Any}(settings["config"])
                settings["config"][flag] = value
                print(stderr, "[", done, "/", total, "] ", label, " ", name, "\n")
                try
                    results = run_simulations(settings, Float32, Int32)
                    for var in sort(collect(keys(results)))
                        var == "timestamps" && continue
                        data = Float64.(results[var])
                        @printf(io, "%-18s %-8s %-10s %-18s %8d %14.6g %14.6g %14.6g\n",
                            label, name, var, string(hash(results[var]), base = 16),
                            count(!isfinite, data),
                            minimum(data), maximum(data), sum(data) / length(data))
                    end
                catch err
                    msg = first(split(replace(sprint(showerror, err), r"\s+" => " "), " Stacktrace"))
                    @printf(io, "%-18s %-8s %-10s ERROR: %s\n", label, name, "-", first(msg, 140))
                end
                flush(io)
            end
        end
    end

    return path

end

"""
    generate(path)

Run all three tile configurations and write an output manifest to `path`.
"""
function generate(path)

    commit = try
        strip(read(`git -C $PROJDIR rev-parse HEAD`, String))
    catch
        "unknown"
    end

    landuse = load_domain_data()
    Nx, Ny = size(landuse["elevation"]["data"])

    open(path, "w") do io
        println(io, "# FSM output baseline manifest")
        println(io, "# source commit : ", commit)
        println(io, "# domain        : ", Nx, " x ", Ny, "  (data/domain_data.nc)")
        println(io, "# precision     : Float32 / Int32")
        println(io, "# tile configs  : as in test/test_regression.jl")
        println(io, "#")
        println(io, "# hash is Base.hash of the full output array and is the authority.")
        println(io, "# min/max/mean are for triage when a hash differs; they are not the check.")
        println(io, "#")
        @printf(io, "# %-8s %-10s %-18s %8s %14s %14s %14s\n",
            "tile", "variable", "hash", "nonfinite", "min", "max", "mean")

        for (name, settings) in TILE_SETTINGS
            results = run_simulations(settings, Float32, Int32)
            for var in sort(collect(keys(results)))
                var == "timestamps" && continue
                data = Float64.(results[var])
                @printf(io, "%-10s %-10s %-18s %8d %14.6g %14.6g %14.6g\n",
                    name, var, string(hash(results[var]), base = 16),
                    count(!isfinite, data),
                    minimum(data), maximum(data), sum(data) / length(data))
            end
        end
    end

    return path

end

rows(path) = [l for l in eachline(path) if !startswith(l, "#") && !isempty(strip(l))]

"""
    check(baseline_path)

Generate a manifest from the current tree and compare its data rows against `baseline_path`.
Returns `true` when every row matches. Comment lines (including the source commit) are ignored.
"""
function check(baseline_path = DEFAULT_BASELINE)

    isfile(baseline_path) || error("baseline not found: $baseline_path")

    current_path = tempname() * ".txt"
    generate(current_path)

    expected, actual = rows(baseline_path), rows(current_path)

    println("baseline : ", baseline_path)
    println("rows     : ", length(expected), " expected, ", length(actual), " produced")

    if length(expected) != length(actual)
        println("MISMATCH : row counts differ")
        return false
    end

    bad = [(e, a) for (e, a) in zip(expected, actual) if e != a]

    if isempty(bad)
        println("RESULT   : IDENTICAL — all ", length(expected), " rows match")
        return true
    end

    println("RESULT   : ", length(bad), " of ", length(expected), " rows differ\n")
    for (e, a) in bad
        println("  expected  ", e)
        println("  actual    ", a)
    end
    return false

end

if abspath(PROGRAM_FILE) == @__FILE__
    action = isempty(ARGS) ? "check" : ARGS[1]
    if action == "generate"
        length(ARGS) >= 2 || error("usage: baseline.jl generate <path>")
        generate(ARGS[2])
        println("wrote ", ARGS[2])
    elseif action == "matrix"
        length(ARGS) >= 2 || error("usage: baseline.jl matrix <path> [FLAG ...]")
        generate_matrix(ARGS[2], ARGS[3:end])
        println("wrote ", ARGS[2])
    elseif action == "check"
        ok = length(ARGS) >= 2 ? check(ARGS[2]) : check()
        exit(ok ? 0 : 1)
    else
        error("unknown action '$action'; expected 'generate' or 'check'")
    end
end
