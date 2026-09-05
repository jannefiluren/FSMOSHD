---
paths:
  - src/**/*.jl
---

# Kernel point-function rules

FSM physics runs in KernelAbstractions kernels launched over the 2-D `(i, j)`
grid. Each process is a thin launcher + a `@kernel` that dispatches to small
`@inline` "point functions". These rules keep those functions uniform, testable,
and GPU-safe.

## Kernel structure

- **Launcher** `foo!(fsm, ...)` (exported): pick the backend, launch, synchronize.
  No physics. Pass the sub-structs the kernel needs — `fsm.state`, `fsm.diag`,
  `fsm.landuse`, `fsm.grid`, `fsm.params` (and `meteo`) — plus isbits scheme
  objects and any `Val`/computed scalars. Do not unpack arrays in the launcher.
- **Kernel** `foo_kernel!` (`@kernel`): `i, j = @index(Global, NTuple)`,
  `@unpack_constants(Tf)`, destructure the fields it uses from the passed
  sub-structs (`(; Tsrf, Ds) = state`), then dispatch to point functions.
- **Point functions** (`@inline`): the physics. See the rules below.

## Point-function rules (array-update functions)

For dispatched physics that reads/writes the model arrays:

- **A — no ad-hoc bundles.** Take the model sub-structs (`state`, `diag`,
  `landuse`, `meteo`, `params`) + `i, j`, and index inside (`state.Tsrf[i, j]`).
  Do **not** build re-bundled NamedTuples in the kernel (e.g.
  `alb_states = (; Tsrf, ...)`). The scheme object is passed for dispatch.
  Physical constants come from `@unpack_constants(Tf)` **inside** the function
  (precedent: `qsat`). A lone computed scalar that lives in no struct (e.g.
  `summer_decay`) is passed explicitly.
- **B — may mutate in place.** Write outputs into `state`/`diag` arrays and
  `return nothing`. FSM is operator-splitting, not a tendency machine, so in-place
  is the idiom. `landuse` is read-only after `setup`, and **`meteo` (MET) is never
  written** — it is meteorological forcing (input only). Never assign to a MET
  field from a kernel or point function; the `MET Immutability` tests enforce this.
- **C — take `i, j`, never `k`.** A point function owns the whole column at a
  cell: loop the layers internally and update all of them. Kernels stay 2-D
  because `Nsnow` is a runtime per-cell height, so the layer loop cannot be a
  kernel dimension.

Canonical examples: `snow_conductivity!` (thermal.jl), `snow_albedo!`
(radiation.jl), `snowcoverfraction_point!` (snowcoverfraction.jl).

### What the rules do NOT cover

Rules A–C apply to **array-update** point functions. Two other *kinds* of helper
are exempt — classify a new function by its shape, not by name:

- **Scalar-transformation functions** — take scalars (or bundles of *derived
  intermediates*) and return values or bundles; they do no `[i, j]` indexing and
  touch no model array. Chaining them, where each stage returns its own bundle for
  the next, is an encouraged pattern (it keeps each stage's intermediates separate
  and independently testable). *Example:* the surface exchange chain —
  `reference_heights` → `surface_layer_state` → `stability_factor` /
  `canopy_richardson` → `eddy_diffusivities`, where `surface_layer_state` returns
  an `S` bundle and `eddy_diffusivities` a `K` bundle. Scheme accessors
  (`canopy_fsar`) and `fresh_snow_density` are simpler members of this kind.
- **Index-agnostic numerical utilities** — solvers, reductions, and math on plain
  buffers/scalars, independent of the grid. *Examples:* `tridiag!`, `ludcmp!`,
  `qsat`, `column_sum`, `first_argmin`/`first_argmax`.

## GPU / correctness (all kernel code)

- **`Tf`-purity:** no literal `0.0`/`1.0`/`2` in kernels or point functions — use
  `Tf(0)`, `Tf(0.1)`, `zero(Tf)`. A Float64 literal silently promotes the whole
  expression to Float64.
- **Physics schemes are isbits** (scalars only): they cross into a kernel by
  value. Per-cell array maps (e.g. `afs`, `adc`) live on `Landuse`, not on the
  scheme.
- **Passing a sub-struct into a kernel** requires `@adapt_structure` on that
  struct (see types.jl); Adapt rewrites its array fields to the device array type
  at launch (a no-op on CPU). `Parameters` and the schemes are isbits and need no
  adaptor.
- **`@kernel inbounds = true`** for kernels with `MVector` scratch — never a raw
  `@inbounds` block inside a `@kernel` body (it miscompiles with the KA CPU
  aliasscope on Julia ≥ 1.11). Keep helper functions `@inline` so their `MVector`
  scratch stays off the heap.

## Bit-identity

Every refactor must stay **bit-identical** to the physics baseline. Verify with
`Pkg.test`, `test/baseline.jl check` (IDENTICAL vs commit `6c4dda5`), and — for
schemes the baseline matrix cannot cover (ALBEDO/CONDCT are schemes, not integer
flags) — a git A/B hash of a synthetic run across scheme combinations.
