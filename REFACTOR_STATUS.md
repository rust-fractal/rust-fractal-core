# Refactor status

This branch started from the one-commit PR refactor (`7a539bf`) and was
audited against `master` before implementation. The refactor keeps the
numerical core, native scheduling, I/O, colouring, root finding, and fuzz
manifests separate; the old renderer/UI API is not retained.

## Preserved from `master`

- Mandelbrot and Burning Ship, powers 2 through 16.
- Smooth iteration, distance estimation, stripes, validity, supersampling,
  palettes, raw channels, image export, root finding, and legacy TOML loading.
- Extended-range arithmetic, perturbation, series approximation, glitch
  detection, native tile scheduling, progressive events, cancellation, and
  deterministic raw output.
- High-precision reference arithmetic through Malachite, now integrated into
  the normal `RenderSession` path rather than being benchmark-only.

## Implemented in this pass

- `ComplexCoordinate` and `Zoom` preserve validated decimal components/scales,
  including values such as `1E8000`; shallow values still use an `f64` hint.
- `RenderSession` selects precision from `PrecisionPolicy`, parses the exact
  reference once, stores compact extended points plus sparse Malachite
  checkpoints, and rejects incomplete deep references instead of falling back
  to invalid absolute `f64` iteration.
- Deep viewport transforms operate on center-relative extended offsets. Tests
  verify adjacent e8000 pixels remain distinct and deep production frames are
  nontrivial.
- Quadratic and cubic perturbation have specialized extended kernels; generic
  power-N and Burning Ship paths remain explicit fallbacks. Deep glitches
  rebase into extended direct recurrence and report correction metadata.
- Series order is honored, coefficient states are checkpointed (interval 32
  for long references), and deep evaluation can start from a validated sparse
  state. Unsupported formula/order combinations return typed errors.
- Runtime delivery is coordinator-driven and bounded. Completion-order streams
  provide immediate delivery; the existing synchronous API keeps deterministic
  tile order with a small reorder buffer. A full-frame path fills disjoint row
  views without tile allocation/merge copies.
- Cancellation is checked between pixels and periodically inside iteration
  loops. The simple token uses relaxed atomic ordering.
- Benchmarks cover the production deep render path, reference generation,
  scalar recurrence alternatives, allocation counts, runtime scaling, and
  full-frame versus progressive execution. CI now targets Ubuntu and Windows.

## Removed intentionally

- The old `FractalRenderer`/UI/progress/data-export API and its `rug`-based
  public implementation.
- The attempted intra-iteration Rayon reference prototype. On this machine it
  was slower than scalar arithmetic through 65,536-bit operands because clone,
  scheduling, and synchronization costs dominated. Production reference
  arithmetic therefore remains single-threaded and the dependency chain is
  respected.

## Previous losses corrected

The initial PR could preserve deep strings only in `fractal-io`, ignored the
configured series order, stored dense per-iteration coefficient vectors,
used an invalid direct-f64 glitch fallback, recomputed viewport invariants in
the pixel path, materialized all tiles before delivery, and checked
cancellation only between pixels. Those paths are now covered by the exact
coordinate, series, glitch, runtime, transform, and cancellation tests.

## Remaining limitations

- The production reference keeps compact points and a shallow `f64` view for
  compatibility with series preparation; it does not retain every Malachite
  value. Sparse checkpoints are currently stored for future rebase work but
  deep correction uses the compact reference directly, so independent
  concurrent rebase references are not needed by the current kernel.
- Series acceleration is intentionally limited to quadratic Mandelbrot, where
  it has the measured benefit. Other formulas use direct specialized kernels
  unless approximation is disabled.
- The old `fractal-arbitrary` crate remains as a small Malachite/root-finding
  oracle and regression helper. The production render and location benchmarks
  no longer depend on its duplicate orbit wrapper; merging it would be a
  separate API cleanup rather than a numerical requirement.
- No SIMD path or alternate Rug/MPFR backend is retained: representative
  scalar/runtime measurements did not establish a worthwhile portable SIMD or
  backend change, and no unsupported dependency was added.
- Allocation counts are measured by the standalone reference probe; peak RSS
  and allocator-specific live-byte statistics are not claimed.
