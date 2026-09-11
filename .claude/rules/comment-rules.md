---
paths:
  - src/**/*.jl
  - test/**/*.jl
  - script/**/*.jl
---

# Comment rules

Be minimal with comments. The code should speak for itself through clear names
and obvious structure; heavy commentary is a sign the code itself needs work,
not more prose around it. These rules are always in effect, in every file and
every scope.

- **Do not comment an individual obvious line.** If a reader who knows Julia and
  the domain can follow a statement from its names and shape alone, a comment on
  it is noise — leave it out.
- **Do put a single-line banner above a block of related lines to name what it
  computes**, e.g. `# Roughness lengths, friction velocity and canopy wind
  profile` above the block that computes them. A banner states a block's overall
  purpose in domain terms and groups a long function or file into labeled stages;
  it never restates the mechanics of one line.
- **Comment genuinely non-obvious logic**, where a one-line hint saves the reader
  real time: a non-obvious invariant, a subtle index trick, a workaround for a
  specific upstream bug, a numerical-stability detail, a sign or unit convention
  that contradicts intuition, or *why* a bit-identity/GPU constraint forces an
  unusual form.
- **Keep such comments surgical:** one line, placed exactly at the confusing
  step. Never restate what the next line obviously does.
- **No multi-line block comments (`#= =#`) or multi-paragraph prose inside
  function bodies.** The single-line banner above is the only structural comment.
- **Do not narrate the task or history:** no "added for X", "fixes Y", "was
  CANMOD==0 before", and no references to callers. A `TODO` is tolerated during
  active development but should be removed before merging.
- Docstrings on public functions are separate and encouraged — this rule governs
  in-body and inline comments, not docstrings.

Examples:

- ❌ `# increment counter` above `i += 1` — restates one obvious line
- ❌ `# loop over all cells` above `for i in 1:Nx` — restates control flow
- ✅ `# Roughness lengths, friction velocity and canopy wind profile` above the
  block that computes them — names the block's purpose
- ✅ `# 0.05 m snow-depth threshold stabilises tuning-point runs (vs fsnow)`
- ✅ `# keep Usc here: dropping it would break bit-identity vs the baseline`
