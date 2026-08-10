# Arbitrary-Precision Backend Report

## Scope

This report records the current backend decision for the deterministic core
refactor. Malachite is the arbitrary-precision reference oracle for every
checked reference location.

## Current candidates

| Backend | Current state | Decision |
| --- | --- | --- |
| `malachite` 0.10.x | `fractal-numerics/malachite` implementation | Use as the only arbitrary-precision arithmetic |
| `fractal-core::generate_reference_exact` | Production exact-coordinate recurrence and compact orbit | Use for normal rendering |
| `fractal-arbitrary::MalachiteReference` | Independent Malachite/root-finding oracle | Keep for differential tests and legacy experiments |

No project-owned arbitrary/fixed-point backend is maintained. This is an
explicit scope decision: Malachite provides the arbitrary-precision arithmetic,
while this repository owns the fractal recurrence, precision policy, and
backend-facing tests.

## Evidence collected

- Malachite reference points are compared against an `f64` quadratic orbit for
  32 iterations at 256 bits; the current test requires a maximum component
  error below `1e-12`.
- The production wrapper reserves compact reference storage once and performs
  arbitrary arithmetic through Malachite `Float` values; sparse checkpoints
  replace a Malachite clone for every iteration.
- The wrapper is covered by unit tests and Criterion benchmarks for complex
  square-add and 100-iteration reference generation.
- Every reference-location benchmark uses the raw decimal coordinates and the
  production exact-reference wrapper. No reference-location calculation falls
  back to `f64`.

The finite render benchmark measures the real locations used by
`benches/core_render.rs`. The benchmark keeps each case bounded at its declared
pixel dimensions and iteration cap, while retaining the location's center,
zoom, and rotation. The `arbitrary_reference` benchmark calls
`fractal_core::generate_reference_exact`, parsing each location's decimal
coordinates once at the production precision derived from its zoom, so
`flake`, `wfs`, `e8000`, and the other deep locations use the same exact
reference pipeline as rendering.

| Case | Median range reported by Criterion |
| --- | --- |
| extended multiply | 6.123–6.459 ns |
| Malachite complex square-add | 649.81–660.58 ns |
| Malachite 100-iteration reference | 41.189–41.776 µs |
| `locations/e14.toml` (64×64, 512 iterations) | 250.35–251.16 ms |
| `locations/burning_ship_test.toml` (128×128, 256 iterations) | 2.972–3.003 ms |
| `locations/e50.toml` (64×64, 512 iterations) | 194.02–194.34 ms |
| `locations/flake.toml` (64×64, 256 iterations) | 2.1739–2.1773 ms |
| `locations/cubic_test.toml` (power 3, 64×64, 256 iterations) | 27.513–27.818 ms |
| arbitrary `e14` (302 bits, 512 iterations) | 522.66–529.44 µs |
| arbitrary `e50` (422 bits, 512 iterations) | 554.78–557.08 µs |
| arbitrary `flake` (777 bits, 512 iterations) | 758.46–762.95 µs |
| arbitrary `wfs` (5,305 bits, 512 iterations) | 7.249–7.352 ms |
| arbitrary `e8000` (26,832 bits, 256 iterations) | 29.375–29.427 ms |
| arbitrary `e1086` (3,863 bits, 256 iterations) | 1.786–1.865 ms |
| arbitrary `e1121` (3,979 bits, 256 iterations) | 1.789–1.795 ms |
| arbitrary `e1000` (3,578 bits, 256 iterations) | 2.002–2.003 ms |

These are historical Windows local baselines, not cross-backend comparisons or
a claim of final production speedup. Final production timings and allocation
measurements are in `PERFORMANCE.md`. The runs are repeatable with:

```text
cargo bench --bench numerics --features malachite --locked -- --quick
cargo bench --bench core_render --locked -- --quick
cargo bench --bench arbitrary_reference --features malachite -- --quick
```

## Malachite baseline

The `malachite` feature resolves to Malachite 0.10.x. It uses
Malachite's arbitrary-precision `Float` type and does not require a native GMP,
MPFR, MPC, or `m4` installation, which keeps the reference benchmark portable
to Windows. Malachite documents `Float` as implemented and tested but not yet
feature-complete, so this backend is retained as a correctness/reference
oracle rather than silently treating it as the production render backend.

The Malachite dependency remains feature-gated in `fractal-numerics`, while
`fractal-core` enables it for the production exact-reference path. Both
wrappers parse raw decimal coordinates at the selected precision; the
benchmark does not use an `f64` escape check.

## Adoption decision

Malachite is the selected arbitrary-precision implementation for the current
scope. A project-owned custom backend is explicitly rejected: it would
duplicate a general arbitrary-precision library and expand maintenance without
evidence that it improves the real reference workloads. Allocation measurements
and end-to-end exact-render integration are complete; peak RSS and alternate
backend comparisons remain outside the supported artifact.
