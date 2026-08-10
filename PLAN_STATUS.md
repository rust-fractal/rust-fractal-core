# Refactor plan status

This is the implementation status for the refactor plan in `PLAN.md` after the
scope decisions recorded at its top.

Completed in the current Windows-native scope:

- deterministic `fractal-core` request, frame, tile, reference, perturbation,
  series, root, colouring, supersampling, and raw-frame APIs;
- immutable `RenderSession` state shared across tiles;
- pixel-boundary cancellation and ordered runtime events;
- JSON/TOML settings, legacy location parsing that preserves raw decimal
  strings, palette conversion, PNG/JPEG/EXR export boundaries, and fuzz targets;
- Malachite-only arbitrary reference orbits and root finding;
- reference-location benchmarks using `e14`, `e50`, `flake`, `wfs`, `e8000`,
  `e1086`, `e1121`, and `e1000`, with raw decimal input and no `f64` fallback;
- Windows Criterion baselines recorded in `ARBITRARY_PRECISION_REPORT.md`;
- no GMP, MPFR, MPC, or project-owned custom/fixed arbitrary backend in
  the active source or dependency graph.

Explicitly deferred or rejected by project decision:

- WASM targets, browser harnesses, and WASM-specific benchmarking;
- a project-owned fixed/custom arbitrary backend;
- replacing the production `f64`/perturbation render path with Malachite. The
  Malachite path is the reference-location oracle and is not silently used as
  the production renderer.

Windows acceptance commands:

```text
cargo fmt --all -- --check
cargo check --workspace --all-targets --locked
cargo test --workspace --all-targets --all-features --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo bench --bench numerics --features malachite --locked -- --quick
cargo bench --bench core_render --locked -- --quick
cargo bench --bench arbitrary_reference --features malachite --locked -- --quick
```
