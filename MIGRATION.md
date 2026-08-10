# Migration from the historical renderer

The old `FractalRenderer` combined configuration parsing, reference creation,
iteration, threading, colouring, and file output. The refactored workspace
separates those responsibilities.

| Historical responsibility | New boundary |
| --- | --- |
| `config::Config` and split TOML files | `fractal_io::request_from_legacy_toml` |
| viewport and render flags | `fractal_core::RenderRequest` and `RenderSettings` |
| raw iteration/data buffers | `fractal_core::RawFrame` |
| palette and RGBA mapping | `fractal_core::ColourPipeline` |
| Rayon, cancellation, and progress | `fractal_runtime::RenderExecutor` |
| PNG/EXR/filesystem output | `fractal_io` |
| high-precision production reference baseline | `fractal_core::generate_reference_exact` |
| independent arbitrary-precision oracle/root helpers | `fractal_arbitrary::MalachiteReference` |

Minimal synchronous rendering now looks like:

```rust
use fractal_core::{
    ComplexCoordinate, RenderRequest, RenderSettings, Viewport, Zoom, compute_frame,
};

let request = RenderRequest {
    generation: 1,
    viewport: Viewport {
        center: ComplexCoordinate::new(-0.75, 0.0),
        zoom: Zoom::unit(),
        rotation_radians: 0.0,
        width: 800,
        height: 600,
    },
    settings: RenderSettings::default(),
};
let raw_frame = compute_frame(&request)?;
```

For native progressive rendering, pass the request to
`fractal_runtime::RenderExecutor` and consume `RenderEvent`s. Consumers should
discard events whose generation is no longer current; the runtime provides
`GenerationGate` for this policy.

The historical renderer source has been removed from the build after its
responsibilities were split across the workspace. The new API is intentionally
breaking; no compatibility feature is provided for the old renderer.
