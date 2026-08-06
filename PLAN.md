# rust-fractal-core Refactor and Arbitrary-Precision Evaluation Plan

**Repository:** `rust-fractal/rust-fractal-core`  
**Plan date:** 2026-08-06  
**Primary goal:** turn the existing renderer into a clear, fast, well-tested, UI-independent computational library that can support both a native Bevy application and a later WebAssembly target.

### Execution decisions

The implementation scope was narrowed by later project decisions:

- Malachite is the sole arbitrary-precision implementation. The project-owned
  custom/fixed-point backend experiment is explicitly rejected and must not be
  introduced as a second arbitrary backend.
- Reference locations are always calculated through Malachite from their raw
  decimal strings; there is no `f64` reference-location fallback.
- WebAssembly targeting is deferred. No active WASM target, browser harness, or
  WASM-specific dependency is required until that decision is revisited.

---

## 1. Objectives

This work has two linked objectives:

1. **Refactor the core renderer**
   - improve clarity and maintainability;
   - reduce production lines of code and repeated logic;
   - remove dead code, obsolete comments, and unnecessary abstraction;
   - make ownership, state transitions, and numerical assumptions explicit;
   - improve native speed and memory efficiency;
   - separate algorithms from configuration, threading, UI, file I/O, and image export;
   - make the computational core buildable for `wasm32-unknown-unknown`.

2. **Evaluate and, where justified, build arbitrary-precision implementations**
   - retain `rug`/MPFR/MPC as the initial correctness and native-performance baseline;
   - benchmark realistic fractal kernels, not generic big-number demonstrations;
   - compare pure-Rust candidates for native and WASM use;
   - prototype a project-owned, fractal-specialized arbitrary/fixed-precision backend;
   - investigate portable, architecture-specific, SIMD, and carefully thresholded multithreaded kernels;
   - select the backend or backend combination using measured accuracy, runtime, allocation, build, portability, and maintenance results.

The work must preserve the mathematical behaviour of perturbation iteration, series approximation, glitch correction, root/period finding, colouring inputs, and deep-zoom coordinate handling.

---

## 2. Current-state observations

The current implementation has strong algorithms, but several concerns are mixed together.

### 2.1 Core/application coupling

`FractalRenderer` currently owns or directly coordinates:

- parsed `config::Config`;
- viewport and render settings;
- arbitrary-precision reference data;
- series approximation;
- period finding;
- render-index generation;
- progress atomics;
- shared `DataExport`;
- Rayon execution;
- channels and native threads;
- palette generation;
- output/export behaviour.

This makes the renderer difficult to test in isolation and unsuitable as the direct dependency of a browser frontend.

### 2.2 Numerical backend coupling

The public aliases are currently direct aliases:

```rust
pub type ComplexArbitrary = rug::Complex;
pub type FloatArbitrary = rug::Float;
```

`Reference` and root finding directly construct and mutate `rug` values. This prevents a WASM-friendly alternative without editing the algorithms themselves.

### 2.3 Shared mutable state

The core currently uses `Arc<Mutex<DataExport>>`, atomics, channels, and threads inside the rendering implementation. This creates:

- lock contention and opaque ownership;
- difficult cancellation behaviour;
- coupling to the native threading model;
- awkward progressive updates;
- poor separation between deterministic computation and scheduling.

### 2.4 Large mixed-purpose modules

Examples include:

- `renderer.rs`: configuration, orchestration, algorithms, threading, progress, palette setup, and output;
- `data_export.rs`: render data storage, colouring, lighting, PNG/JPEG writing, EXR writing, and buffer lifecycle;
- `series_approximation.rs`: coefficient generation, probe construction, validation, interpolation, progress, and diagnostics.

These modules should be split by responsibility, but not fragmented into excessive one-function files.

### 2.5 Hot-path inefficiencies visible in the current code

Items to verify with benchmarks include:

- repeated cloning of arbitrary-precision values;
- repeated allocation and cloning of series coefficient vectors;
- per-pixel temporary `Vec` allocation for stripe colouring;
- complete buffer replacement in `clear_buffers`;
- cloning whole output channels before EXR writing;
- frequent `SeqCst` atomics where weaker ordering or non-atomic local counters may suffice;
- repeated locking of shared output buffers;
- runtime branches nested inside heavily iterated kernels;
- `Vec<f64>` arguments where slices or fixed small arrays are sufficient;
- linear searches for extended-iteration boundaries;
- parsing zero values from strings rather than constructing them directly;
- use of native C `frexp`/`ldexp` FFI when a safe portable implementation is preferable;
- mixed compile-time and runtime fractal dispatch.

No optimization is accepted solely because it “looks faster”; every hot-path change must be measured.

### 2.6 Correctness and maintainability concerns

The current code contains:

- `unwrap`/`expect` in library-facing paths;
- `unimplemented!()` branches in extended perturbation modes;
- commented-out alternatives and debugging output;
- implicit iteration-index conventions;
- implicit precision selection;
- no obvious comprehensive test suite;
- I/O and configuration errors that cannot be handled cleanly by callers.

---

## 3. Non-goals

This phase will **not**:

- build the Bevy GUI;
- move the perturbation algorithm to the GPU;
- redesign the visual appearance of palettes;
- change the mathematical method merely to reduce code;
- promise identical floating-point bits across all arbitrary-precision backends;
- force a single arbitrary-precision implementation before benchmarks exist;
- add broad generic numeric traits to every hot loop;
- build a general-purpose replacement for GMP, MPFR, or MPC;
- parallelize inherently sequential reference iterations merely to claim multithreading;
- retain old APIs solely for compatibility.

Breaking changes are acceptable. Dead APIs should be removed rather than deprecated.

---

## 4. Target repository structure

Convert the repository into a small Cargo workspace. Keep the number of crates limited and responsibility-driven.

```text
rust-fractal-core/
├── Cargo.toml
├── crates/
│   ├── fractal-core/
│   │   ├── domain types
│   │   ├── extended-range arithmetic
│   │   ├── perturbation
│   │   ├── series approximation
│   │   ├── reference generation
│   │   ├── glitch correction
│   │   ├── root/period finding
│   │   └── deterministic tile/frame computation
│   │
│   ├── fractal-numerics/
│   │   ├── arbitrary-precision backend interface
│   │   ├── rug backend
│   │   ├── pure-Rust candidate backends
│   │   └── backend conformance tests
│   │
│   ├── fractal-runtime/
│   │   ├── native executor
│   │   ├── Rayon scheduling
│   │   ├── cancellation
│   │   └── progressive render events
│   │
│   └── fractal-io/
│       ├── settings serialization
│       ├── PNG/JPEG export
│       ├── EXR export
│       └── native filesystem integration
│
├── benches/
│   ├── arithmetic/
│   ├── kernels/
│   ├── rendering/
│   └── backends/
│
├── tests/
│   ├── fixtures/
│   ├── regression/
│   └── differential/
│
└── fuzz/
```

A smaller alternative is acceptable if the implementation shows that `fractal-numerics` can remain a module of `fractal-core`. `fractal-io` and the native executor must remain separable from the WASM-capable core.

---

## 5. Target public API

The core API should be based on explicit domain types rather than `config::Config`.

```rust
pub struct Viewport {
    pub center: ComplexCoordinate,
    pub zoom: Zoom,
    pub rotation_radians: f64,
    pub width: u32,
    pub height: u32,
}

pub struct RenderSettings {
    pub fractal: FractalFormula,
    pub iteration_limit: u32,
    pub approximation: ApproximationSettings,
    pub glitch: GlitchSettings,
    pub sampling: SamplingSettings,
    pub output: OutputChannels,
}

pub struct RenderRequest {
    pub generation: u64,
    pub viewport: Viewport,
    pub settings: RenderSettings,
}

pub enum RenderEvent {
    Started(RenderStarted),
    Progress(RenderProgress),
    TileReady(RenderTile),
    Finished(RenderResult),
    Cancelled { generation: u64 },
    Failed(RenderError),
}
```

### Required API properties

- `RenderRequest` is immutable after submission.
- A render generation uniquely identifies all results.
- A late result can be rejected without mutating the renderer.
- The core can synchronously compute a reference, series state, tile, or complete frame.
- Scheduling and message delivery are supplied by `fractal-runtime`.
- Rendering can be deterministic regardless of worker count and tile order.
- Raw numerical channels are separate from colour mapping.
- Recolouring does not rerun fractal iteration.
- Export is outside the numerical core.
- Invalid settings return typed errors.
- Public types support serialization where appropriate, but core algorithms do not depend on a particular file format.

---

## 6. Numerical architecture

### 6.1 Avoid dynamic dispatch in hot loops

Do not place `dyn ArbitraryFloat` or `dyn ArbitraryComplex` inside reference iteration. Instead use one of:

- feature-selected concrete backend types;
- generic backend implementations monomorphized at compile time;
- backend-specific adapters that expose a common higher-level operation set.

The abstraction should sit around coarse kernels, not around every add and multiply.

### 6.2 Required arbitrary-precision operations

Derive the backend requirements from actual code. At minimum:

#### Float operations

- construction at a requested bit precision;
- zero/one/integer construction without parsing strings;
- clone and in-place assignment;
- add/subtract/multiply/divide;
- square;
- `2^n` construction or exponent adjustment;
- precision query and precision change;
- finite/NaN/infinite checks;
- comparison;
- conversion to `f64`;
- conversion to mantissa plus binary exponent;
- string parse/format for persisted coordinates.

#### Complex operations

Implement a thin project-owned complex type over backend floats where a candidate has no native complex type:

```rust
pub struct BigComplex<F> {
    pub re: F,
    pub im: F,
}
```

Required operations are intentionally narrow:

- add/subtract;
- square in place;
- multiply;
- divide for nucleus finding;
- norm squared;
- precision management;
- conversion to fixed and extended representations.

Do not require a large generic complex-number ecosystem unless the renderer actually uses it.

### 6.3 Precision policy

Move precision selection into an explicit policy:

```rust
pub struct PrecisionPolicy {
    pub minimum_bits: u32,
    pub guard_bits: u32,
    pub max_bits: Option<u32>,
}
```

The policy should compute the required precision from viewport radius/zoom and algorithmic needs. Tests must verify monotonicity and sufficient guard precision.

### 6.4 Extended-range types

Keep the useful mantissa/exponent representation, but refactor it as a small independent numerical module.

Requirements:

- safe, portable mantissa/exponent decomposition;
- no C FFI;
- invariants documented;
- normalization/reduction centralised;
- `#[must_use]` where ignored results would be suspicious;
- no silently invalid zero exponent states;
- property tests for arithmetic and ordering;
- clear distinction between value range extension and significand precision.


### 6.5 Project-owned specialized precision backend

The arbitrary-precision study must include a custom backend, tentatively named `fractal-arbitrary`. It should be designed around the actual fractal kernels rather than around compatibility with a general-purpose floating-point API.

The first prototype should target **reference-orbit generation for bounded escape-time fractals**, especially:

```text
z_{n+1} = z_n² + c
```

For this workload, the most promising representation is a fixed-scale binary value rather than a fully general normalized arbitrary float.

#### Proposed reference representation

```rust
pub struct FixedBig {
    sign: Sign,
    limbs: Box<[u64]>,
}

pub struct FixedComplex {
    re: FixedBig,
    im: FixedBig,
}

pub struct FixedContext {
    fractional_bits: usize,
    guard_bits: usize,
    limb_count: usize,
    rounding: RoundingMode,
}
```

The exact representation should be selected through prototypes, but it should have these properties:

- all values in one reference calculation use one known binary scale;
- the precision and guard-bit count are fixed for the duration of a calculation;
- values are stored in contiguous little-endian limbs;
- storage is allocated once and reused;
- no per-iteration heap allocation;
- no decimal parsing in the iteration loop;
- no per-operation precision negotiation;
- no per-operation exponent normalization;
- explicit overflow/escape detection;
- deterministic rounding;
- portable fallback implementation;
- optional architecture-specific optimized kernels.

Because the reference orbit remains bounded until escape, fixed-point arithmetic may avoid much of the machinery required by MPFR-style general floats. This hypothesis must be benchmarked and verified rather than assumed.

#### Narrow operation set

The first custom backend should implement only what the reference kernel requires:

- fixed-scale signed add and subtract;
- absolute comparison;
- left/right shifts with defined rounding;
- specialized square;
- fixed-scale multiply followed by scale restoration;
- fused complex square and add;
- conversion from a high-precision parsed coordinate;
- conversion to `f64` and the project extended-range type;
- exact or conservatively bounded escape tests;
- reusable scratch storage.

A possible fused complex-square kernel is:

```text
real = (x + y)(x - y) + c_re
imag = 2xy + c_im
```

This uses two large multiplications rather than independently expanding both components. Benchmark it against the conventional formulation because operand growth, signs, carry propagation, and truncation may change the practical result.

Do not add arbitrary division, transcendental functions, generic decimal formatting, or a broad trait hierarchy to the first prototype. Root/nucleus finding may continue to use `rug` or the selected general backend until a dedicated custom root kernel is separately justified.

#### Multiplication algorithm ladder

Implement algorithms incrementally and retain only measured winners:

1. portable basecase/schoolbook multiplication;
2. specialized basecase squaring;
3. Comba-style accumulation where it is faster and safe;
4. Karatsuba above a measured limb threshold;
5. Toom variants only if the project precision distribution justifies them;
6. FFT/NTT multiplication only for very large precisions where it wins end to end;
7. truncated or middle-product algorithms only after a written error proof and differential tests.

GMP uses different multiplication and squaring algorithms at size-dependent thresholds. The custom backend should likewise use benchmark-derived thresholds rather than a single algorithm at every size.

#### CPU specialization

Evaluate isolated implementations for:

- portable `u64` limbs with `u128` accumulation;
- x86-64 BMI2 multiplication;
- x86-64 ADX carry chains;
- combined BMI2+ADX paths;
- AArch64 multiply-high/carry operations;
- WASM32 portable scalar;
- WASM SIMD128 where it benefits independent limb work or batched kernels.

Requirements:

- runtime feature detection on native targets;
- a correct portable fallback;
- all `unsafe` code isolated in a small architecture module;
- every optimized kernel tested against the portable implementation;
- no architecture-specific path accepted without end-to-end reference-orbit improvement.

SIMD is expected to help more with multiple independent values or batched work than with a single carry-dependent limb chain. Benchmark both arrangements.

#### Memory and scratch design

Create caller-owned reusable workspaces:

```rust
pub struct FixedScratch {
    product: Vec<u64>,
    temporary_a: Vec<u64>,
    temporary_b: Vec<u64>,
    carry: Vec<u64>,
}
```

The final design should:

- reserve all required capacity before iteration;
- reuse output and temporary buffers;
- support swap-based iteration;
- avoid zeroing memory that will be fully overwritten;
- keep hot limbs contiguous and aligned where useful;
- avoid reference-counted limb storage in hot paths;
- report workspace memory as part of benchmarks.

#### Rounding and error budget

The custom backend must have a documented numerical contract.

At minimum:

- define the fixed binary scale;
- define guard bits;
- define rounding after multiplication;
- track or bound truncation error;
- prove that escape and glitch decisions remain valid under the chosen bound;
- compare every reference checkpoint with a higher-precision `rug` oracle;
- increase precision or fall back to the general backend when the error budget is exhausted.

A fast backend that silently changes orbit classification is unacceptable.

### 6.6 Multithreading policy for arbitrary precision

Use parallelism at the coarsest profitable level first.

#### Level 1 — independent high-precision jobs

Prioritize parallel execution of:

- glitch-correction reference orbits;
- independent candidate/reference calculations;
- independent root searches;
- separate frames;
- independent benchmark cases.

This usually offers more parallel work and less synchronization than splitting one arithmetic operation.

#### Level 2 — parallel large multiplication

Prototype operation-level parallelism only above measured thresholds:

- parallel Karatsuba branches;
- parallel Toom evaluation/interpolation stages;
- parallel FFT/NTT transforms;
- independent real/imaginary products in a complex kernel when the fused two-product form is not used.

The central reference recurrence itself is sequential: iteration `n + 1` depends on iteration `n`. It must not be artificially partitioned across threads.

#### Thread-budget ownership

Prevent oversubscription:

```rust
pub struct ArithmeticThreadBudget {
    pub max_workers: usize,
    pub minimum_parallel_limbs: usize,
}
```

- the outer render executor owns the global thread budget;
- a custom arithmetic operation requests workers rather than creating its own unconstrained pool;
- when tiles or correction references already saturate the CPU, large-number operations remain single-threaded;
- nested Rayon pools are avoided;
- thread thresholds are benchmarked separately for each multiplication algorithm and architecture.

#### Web policy

For the first browser release:

- run the custom backend single-threaded inside a Web Worker;
- parallelize independent jobs across workers if useful;
- keep shared-memory operation-level threading optional;
- require explicit benchmarks before enabling WASM atomics/shared-memory code;
- document the hosting headers required by browser shared memory.

Rust exposes WASM atomics and SIMD, but multithreaded `wasm32-unknown-unknown` still requires deliberate toolchain and browser setup. The custom backend must not make shared-memory WASM a baseline requirement.

---

## 7. Work phases

## Phase 0 — Freeze behaviour and establish baselines

Before structural changes:

- record current commit and dependency lockfile;
- build the current release binary;
- capture production LOC by module;
- capture binary size and clean build time;
- record representative render times and peak memory;
- save canonical raw numerical outputs, not only PNG images;
- record reference trajectories and root-finding results;
- document known unsupported branches and current failures.

Create a fixture manifest containing:

- viewport centre;
- zoom;
- rotation;
- image size;
- formula and power;
- iteration count;
- approximation settings;
- expected raw-channel hashes;
- expected selected pixel values;
- expected reference and period statistics.

Suggested fixture classes:

1. ordinary Mandelbrot exterior;
2. ordinary boundary region;
3. interior region;
4. moderate perturbation zoom;
5. very deep zoom;
6. known glitch-heavy location;
7. series approximation enabled and disabled;
8. distance estimate;
9. stripe and distance-stripe colouring inputs;
10. Burning Ship;
11. power 3;
12. generic power greater than 3;
13. root finding at short period;
14. root finding at long period.

**Gate:** no refactor starts until the baseline suite can reproduce the current accepted outputs.

---

## Phase 1 — Modernize the package without changing algorithms

- move to the current stable Rust edition used by the selected Bevy version;
- remove nightly-only requirements unless measured and justified;
- update dependencies in controlled groups;
- add `rust-toolchain.toml`;
- add `rustfmt.toml`;
- enable strict Clippy in CI;
- define MSRV policy;
- replace obsolete configuration APIs;
- remove generated/build metadata that is not required at runtime;
- remove dead code and commented-out implementations;
- convert print debugging to structured tracing in outer crates only;
- add a typed error enum.

Suggested CI matrix:

- Windows x86-64;
- Linux x86-64;
- macOS where practical;
- `wasm32-unknown-unknown` check for WASM-capable crates;
- stable Rust;
- optional nightly Miri/fuzz jobs.

**Gate:** all baseline fixtures remain unchanged.

---

## Phase 2 — Introduce typed domain configuration

Replace `config::Config` inside the core with typed structures.

Tasks:

- define `Viewport`, `RenderSettings`, and nested setting types;
- validate width, height, zoom, precision, iteration limits, powers, and sampling;
- move TOML/config parsing to `fractal-io`;
- provide explicit defaults through `Default`;
- provide conversion code from the current `start.toml`;
- remove string-based mode dispatch;
- replace integer magic values for data/fractal/export modes with enums;
- document units and angle conventions;
- make iteration indexing explicit and consistently zero- or one-based.

**Gate:** `fractal-core` has no dependency on `config`.

---

## Phase 3 — Refactor output data and colouring

Split the current `DataExport` responsibilities into:

```text
RawFrame
├── iterations
├── smooth
├── distance
├── stripes
└── validity/glitch metadata

ColourPipeline
├── palette
├── palette mapping
└── lighting

ImageFrame
└── RGBA or RGB pixels

Exporters
├── PNG/JPEG
└── EXR
```

Tasks:

- allocate only channels requested by `OutputChannels`;
- replace per-pixel temporary vectors with fixed arrays or direct scalar operations;
- reuse buffers rather than replacing entire vectors;
- support tile updates into a frame;
- make recolouring a separate deterministic operation;
- remove file writing from raw frame structures;
- avoid cloning full channels during export where exporter APIs permit borrowing;
- make GUI output use RGBA8 unless a measured reason favours RGB8;
- return export errors rather than panic.

**Gate:** recolouring fixtures match the current accepted outputs and no per-pixel heap allocation remains.

---

## Phase 4 — Extract deterministic numerical stages

Split `FractalRenderer` into explicit stages:

1. viewport geometry;
2. arbitrary-precision reference generation;
3. series coefficient generation;
4. approximation validity/probe calculation;
5. initial pixel/tile setup;
6. perturbation iteration;
7. glitch detection;
8. glitch correction;
9. raw-channel finalization;
10. colour mapping.

Each stage should:

- accept explicit inputs;
- return explicit outputs;
- avoid hidden global/shared state;
- be testable synchronously;
- expose timings through returned statistics rather than printing;
- use cancellation at bounded checkpoints;
- avoid locking.

Suggested orchestrator:

```rust
pub struct RenderSession<B: ArbitraryBackend> {
    request: RenderRequest,
    reference: Reference<B>,
    approximation: ApproximationState,
    frame: RawFrame,
    statistics: RenderStatistics,
}
```

`RenderSession` should not own threads or UI handles.

**Gate:** a single-threaded complete render works using only deterministic core APIs.

---

## Phase 5 — Simplify and optimize numerical kernels

Perform measured changes in small commits.

### Reference generation

- construct zeros directly;
- reuse temporaries;
- use in-place operations where backend APIs support them;
- reduce clones in `z = z² + c`;
- benchmark native complex arithmetic against project-owned complex pairs;
- preallocate reference arrays from iteration limits;
- make high-precision checkpoint storage explicit;
- verify the checkpoint indexing used for glitch correction.

### Series approximation

- replace `Vec<Vec<T>>` if a flat contiguous layout improves access;
- eliminate full coefficient-vector cloning each iteration by swapping reusable buffers;
- preallocate coefficient storage;
- isolate coefficient generation from probe validation;
- avoid repeated `contains` scans for selected probes;
- replace diagnostic print methods with optional debug helpers;
- make storage interval assumptions explicit.

### Perturbation

- reduce nested mode branches by selecting a specialized kernel before the pixel loop;
- retain compile-time specialization only where benchmarks prove it helps;
- use slices rather than `&Vec<T>`;
- avoid repeating extended-iteration searches per pixel;
- benchmark binary search, cursor propagation, and precomputed lookup tables;
- replace global atomics inside inner loops with local counts reduced at chunk completion;
- lower atomic ordering where correctness permits;
- validate bounds-check elimination through profiling rather than unsafe indexing;
- preserve dedicated fast paths for power 2 and power 3;
- decide whether generic powers warrant a separate slower kernel.

### Buffers and memory

- reserve exact or estimated capacities;
- reuse render workspaces between renders where ownership remains clear;
- avoid false sharing between workers;
- group pixel fields by access pattern only after cache profiling;
- compare array-of-structs and structure-of-arrays for the perturbation kernel;
- report peak resident memory and allocations.

**Performance gate:** no accepted optimization may regress the relevant canonical workload by more than 3%. End-to-end performance should be at least equal to baseline, with a target improvement of 15–30% where the current architecture causes avoidable allocation/locking overhead.

---


## Phase 5A — Prototype and evaluate `fractal-arbitrary`

This phase is an explicit research gate, not an automatic commitment to maintain a custom bignum library.

### Step 1 — kernel characterization

Instrument real workloads and record:

- precision in bits;
- limb count;
- reference iteration count;
- number of arbitrary squares, multiplies, additions, divisions, conversions, and clones;
- proportion of total render time in arbitrary arithmetic;
- value ranges encountered;
- frequency of root finding versus reference generation;
- number and size of correction references;
- thread availability during each high-precision stage.

This determines whether the custom work targets the actual bottleneck.

### Step 2 — portable fixed-point prototype

Implement a correctness-first backend with:

- `u64` limbs;
- `u128` intermediates;
- reusable scratch;
- fixed binary scale;
- fixed guard precision;
- basecase add/subtract/multiply/square;
- fused complex square-add;
- conversion from and to the oracle backend.

Use no architecture-specific code initially.

### Step 3 — differential validation

For every iteration in canonical reference orbits:

- compare real and imaginary components with `rug` at extra precision;
- record absolute and relative error;
- verify rounding/error bounds;
- compare escape iteration;
- compare generated fixed and extended reference buffers;
- verify downstream perturbation and raw-frame results.

Fuzz precision, coordinates, signs, carries, and limb boundaries.

### Step 4 — optimize one layer at a time

In separate benchmarked commits:

- eliminate remaining allocations;
- specialize squaring;
- test fused complex formulas;
- add truncated products if proven safe;
- add Karatsuba;
- add native CPU kernels;
- add batch/SIMD paths;
- tune algorithm thresholds;
- add prefetch/alignment only if profiling supports it.

Retain a portable reference implementation permanently.

### Step 5 — multithreading experiment

Benchmark:

- multiple independent custom references;
- parallel real/imaginary products;
- parallel Karatsuba branches;
- parallel Toom/FFT only at sufficiently high precision;
- interaction with outer render threads;
- 1, 2, 4, 8, and available-core worker limits.

Record speedup, efficiency, latency, synchronization, and peak memory. Disable operation-level multithreading below measured thresholds.

### Step 6 — root-finding decision

Only after the reference backend proves valuable, decide whether to implement:

- arbitrary complex division;
- Newton update kernels;
- adaptive precision;
- root residual evaluation.

Keep root finding on `rug` or the selected general backend if its runtime contribution is small or a custom implementation would materially expand scope.

### Adoption gate

Adopt `fractal-arbitrary` as a production backend only if all of these hold:

- it passes the complete correctness and regression suite;
- its error contract is documented;
- it is materially faster on important real renders, with a target of at least 1.5× for the high-precision reference stage or a comparably strong WASM advantage;
- it improves end-to-end workloads enough to justify maintenance;
- its production code remains smaller and narrower than a general arbitrary-float library;
- architecture-specific code has a portable fallback;
- multithreading does not regress common precisions;
- one maintainer can understand and modify the limb kernels using the documentation and tests.

If the adoption gate is not met, retain the prototype and benchmark report as research evidence, and use the best maintained external backend.

---

## Phase 6 — Introduce scheduling outside the core

Implement `fractal-runtime` after deterministic rendering works.

### Native executor

- use Rayon or Bevy-compatible compute pools through an executor adapter;
- use tile or chunk tasks with owned output;
- send `RenderEvent`s through a bounded channel;
- throttle progress updates;
- use generation-based cancellation;
- do not share the entire renderer behind a mutex;
- cap worker count through configuration;
- verify deterministic output across 1, 2, and all available worker threads.

### Cancellation

A cancellation token should be cheap to check:

```rust
pub struct CancellationToken {
    generation: u64,
    cancelled: Arc<AtomicBool>,
}
```

Use acquire/release or relaxed ordering only after documenting the memory-safety and visibility requirements. Avoid `SeqCst` by default.

### Progressive rendering

Emit rectangular tiles or scanline blocks:

```rust
pub struct RenderTile {
    pub generation: u64,
    pub rect: PixelRect,
    pub pixels: Box<[u8]>,
    pub raw: Option<RawTile>,
}
```

The UI can discard stale generations. The core never needs a UI callback.

**Gate:** cancellation, repeated zooming, and stale-result rejection have integration tests.

---

## Phase 7 — Separate native I/O

Move into `fractal-io`:

- TOML/RON settings;
- project/session serialization;
- palette file loading;
- PNG/JPEG export;
- EXR export;
- filesystem path handling;
- native command-line support.

Feature-gate large optional exporters:

```toml
[features]
default = ["png"]
png = ["dep:image"]
exr = ["dep:exr"]
```

`fractal-core` must not depend on `image`, `exr`, filesystem APIs, or native configuration crates.

---

## Phase 8 — Testing programme

## 8.1 Unit tests

### Extended arithmetic

Test:

- normalization invariant;
- zero handling;
- sign handling;
- very large and very small exponents;
- addition with widely separated exponents;
- subtraction/cancellation;
- multiplication/division;
- powers;
- norm and norm squared;
- conversion to/from `f64`;
- ordering;
- formatting round trips.

### Geometry and viewport

Test:

- pixel-to-complex mapping;
- centre pixel;
- corners;
- aspect ratio;
- rotation;
- zoom interpolation;
- deep-exponent coordinates.

### Colouring

Test:

- palette wrapping;
- cyclic and non-cyclic palettes;
- smooth iteration;
- stripe interpolation;
- distance lighting;
- black interior;
- invalid/glitch display;
- channel allocation.

## 8.2 Property tests

Use `proptest` for:

- `reduce()` preserving represented value;
- addition and multiplication agreeing with a high-precision oracle within tolerance;
- parse/format round trips;
- zoom conversion monotonicity;
- pixel mapping invertibility within numeric limits;
- deterministic output under varied chunk sizes;
- no panic for valid settings.

Persist failing seeds.

## 8.3 Differential numerical tests

Use `rug` at additional guard precision as the oracle.

Test:

- each arbitrary backend’s primitive operations;
- complex square/add kernel;
- reference iteration;
- conversions to extended precision;
- root-finding residuals;
- period/nucleus calculations.

For perturbation:

1. calculate a high-precision direct orbit for selected pixels;
2. calculate the perturbation result;
3. compare escape iteration, smooth value, and final orbit difference;
4. repeat over boundary, deep-zoom, and glitch-heavy locations.

For series approximation:

1. evaluate direct perturbation without approximation;
2. evaluate with approximation;
3. compare skip decision and reconstructed delta;
4. test around the validity boundary.

## 8.4 Golden regression tests

Store compact raw fixtures rather than only large images.

For each canonical render:

- selected pixel records;
- histogram of escape iterations;
- min/max channel values;
- reference count;
- approximation skip range;
- glitch count;
- deterministic raw-channel hash when bitwise stability is expected;
- tolerance-based comparisons otherwise.

PNG screenshot fixtures are secondary diagnostics, not the only correctness test.

## 8.5 Concurrency tests

Test:

- 1-thread versus N-thread equality;
- chunk sizes from 1 to large blocks;
- cancellation during every render stage;
- rapid sequence of generations;
- no stale tile accepted by the runtime;
- progress monotonicity;
- no deadlocks;
- bounded channel backpressure.

Use Loom only for small synchronization primitives where it adds value; do not attempt to model the full renderer.

## 8.6 Fuzzing

Add `cargo-fuzz` targets for:

- coordinate parsing;
- zoom parsing;
- project/settings deserialization;
- palette parsing;
- extended-number operations;
- arbitrary backend conversion;
- EXR metadata parsing if custom parsing remains;
- random but bounded render settings at tiny image sizes.

Every discovered crash becomes a normal regression test.

## 8.7 Static and dynamic checks

CI or scheduled jobs should include:

- `cargo fmt --check`;
- `cargo clippy --workspace --all-targets --all-features`;
- `cargo test --workspace --all-features`;
- `cargo check --target wasm32-unknown-unknown`;
- Miri for numerical utility modules where supported;
- sanitizer/fuzz smoke tests;
- dependency audit;
- minimum-version or MSRV check if maintained.

---

## 9. Benchmark infrastructure

Use multiple tools because they answer different questions.

### Criterion

Use Criterion for:

- wall-clock latency;
- throughput;
- scaling by precision, image size, and iteration count;
- comparison plots;
- local regression analysis.

### Iai-Callgrind or its maintained successor

Use instruction/cache-oriented measurements on Linux CI for stable structural comparisons. Do not interpret instruction count as a substitute for wall-clock performance.

### Allocation and memory profiling

Use an allocation profiler such as DHAT plus OS-level peak-memory measurement for:

- reference generation;
- series coefficient generation;
- complete frames;
- export;
- backend comparisons.

### Sampling profiler

Use Samply, `perf`, or an equivalent platform profiler on Criterion profile runs. Preserve scripts that reproduce important profiles.

### Benchmark hygiene

- run release builds;
- record compiler, target CPU, dependency versions, backend versions, OS, and thread count;
- disable unrelated applications where possible;
- use fixed CPU power settings on the primary benchmark machine;
- separate single-thread and multi-thread results;
- retain raw Criterion data;
- do not compare results from materially different hardware as direct regressions;
- benchmark both default and native-CPU builds;
- report median, confidence interval, allocations, and peak memory;
- benchmark warm and cold-start paths separately when parsing/initialization matters.

---

## 10. Core performance benchmark suite

## 10.1 Microbenchmarks

### Extended-range arithmetic

- normalize/reduce;
- add/subtract with similar exponents;
- add/subtract with exponent gaps;
- multiply/divide;
- norm/norm squared;
- conversion to/from fixed precision;
- formatting.

### Arbitrary-precision primitives

At 64, 128, 256, 512, 1,024, 2,048, and 4,096 bits:

- construct zero/one/from `f64`;
- parse representative decimal coordinates;
- format representative coordinates;
- clone;
- precision change;
- add/subtract;
- multiply;
- square;
- divide;
- `2^n`;
- conversion to `f64`;
- conversion to mantissa/exponent;
- finite checks.

### Complex kernels

At the same precisions:

- `z = z² + c`;
- complex multiply;
- complex divide;
- norm squared;
- 1,000 chained reference iterations;
- 100,000 chained reference iterations where practical.

Compare:

- native `rug::Complex`;
- two `rug::Float`s in project-owned `BigComplex`;
- each pure-Rust float in project-owned `BigComplex`.

This determines whether MPC convenience costs or saves time for the actual workload.

## 10.2 Algorithm benchmarks

- reference generation;
- reference checkpoint storage;
- series coefficient generation by order and iteration count;
- probe validation;
- perturbation for a fixed tile;
- glitch detection;
- glitch correction;
- raw-channel finalization;
- recolouring;
- root/nucleus finding at several periods;
- period estimation.

## 10.3 End-to-end benchmarks

Use canonical render classes:

| Class | Purpose |
|---|---|
| Shallow interactive | UI responsiveness and overhead |
| Boundary 1080p | ordinary CPU throughput |
| Moderate deep zoom | reference + perturbation balance |
| Very deep zoom | arbitrary-precision cost |
| Glitch-heavy | correction scheduling |
| High iteration | long reference and perturbation loops |
| Distance/lighting | derivative and channel cost |
| Power 3 / generic power | non-default formula cost |
| Root finding | arbitrary complex division and Newton loops |

Measure:

- first visible tile latency;
- total render time;
- CPU utilization;
- peak memory;
- bytes allocated;
- number of reference recalculations;
- number of completed/cancelled tasks;
- scaling from 1 to N threads.

---

## 11. Arbitrary-precision candidate study

## 11.1 Candidates

### A. `rug` 1.30.x — baseline

Properties:

- MPFR/MPC/GMP-backed;
- correct rounding;
- native arbitrary complex type;
- mature API;
- native C dependencies;
- awkward browser/WASM portability.

Evaluate both `rug::Complex` and a project-owned complex pair of `rug::Float`.

### B. `dashu-float` 0.5.x

Properties to test:

- pure Rust;
- arbitrary precision and rounding;
- `no_std` support;
- no native complex type required for this project;
- browser build suitability.

Do not assume transcendental-function performance matters; prioritize the exact operations used by the renderer.

### C. `astro-float` 0.9.x

Properties to test:

- pure Rust arbitrary precision;
- explicit precision and rounding;
- configurable exponent range;
- correctly rounded operations according to its documentation;
- custom complex wrapper required.

Evaluate setup/context overhead as well as arithmetic.

### D. `malachite-float` 0.9.x

Properties to test:

- pure Rust;
- algorithms partially derived from MPFR;
- extensive built-in benchmark/test infrastructure;
- float support currently documented as experimental and incomplete;
- custom complex wrapper required.

It is a benchmark candidate, not an assumed production choice.


### E. Project-owned `fractal-arbitrary` backend

Develop and benchmark a narrowly specialized backend in two variants:

1. **Portable fixed-point reference backend**
   - pure Rust;
   - fixed scale and fixed precision per reference;
   - `u64` limbs;
   - reusable scratch;
   - specialized square and fused complex square-add;
   - intended for native and WASM.

2. **Optimized native backend**
   - architecture-specific multiplication/carry kernels;
   - size-dependent algorithm dispatch;
   - optional thresholded operation-level parallelism;
   - identical external numerical contract.

The first version is not required to support general arbitrary floats or complex division. It may coexist with `rug` for parsing and root finding.

Benchmark the custom backend against:

- `rug::Complex`;
- two `rug::Float`s;
- the best external pure-Rust backend;
- the portable custom implementation;
- every architecture-specific implementation;
- single-threaded and multithreaded modes.

Reject it if the end-to-end improvement does not justify its maintenance cost.

### F. Optional tiered precision experiment

After the arbitrary backends are working, test whether moderate zooms benefit from a tiered policy:

- `f64` or current extended-range arithmetic;
- double-double or quad-double;
- custom fixed-point reference arithmetic;
- general arbitrary precision only beyond a measured threshold.

This is optional and must not complicate the core unless it materially improves real renders.

---

## 11.2 Benchmark data sets

Use fixed decimal coordinate strings and preserve them in the repository.

Include:

- values near zero;
- ordinary Mandelbrot coordinates;
- coordinates with heavy cancellation;
- deep coordinates requiring 128–4,096 bits;
- stress cases at 8,192, 16,384, 65,536, and higher bit counts where practical, specifically to identify multiplication and threading thresholds;
- extremely large negative zoom exponents;
- nucleus guesses for several periods;
- values encountered in actual reference trajectories, not random numbers alone.

Use deterministic random seeds for supplementary generated inputs.

---

## 11.3 Correctness protocol

`rug`/MPFR/MPC at target precision plus at least 64 guard bits is the primary oracle.

For each candidate:

1. parse the same decimal input;
2. run the same primitive operation sequence;
3. compare to the high-precision oracle after rounding;
4. run reference iteration;
5. compare intermediate checkpoints;
6. run nucleus finding;
7. compare final residual;
8. run selected full renders;
9. compare raw numerical channels.

Record:

- maximum relative error;
- maximum absolute error;
- ULP-style error where representation permits;
- first iteration where trajectories diverge beyond tolerance;
- incorrect escape decisions;
- incorrect period/root convergence;
- parse/format loss;
- NaN/infinity behaviour;
- fixed-point guard-bit consumption;
- custom-backend overflow and fallback events;
- equality between portable and architecture-specific custom kernels;
- equality between single-threaded and multithreaded custom kernels.

For `fractal-arbitrary`, compare selected orbits after **every iteration**, not only at the final result. Validate carries, rounding, truncation, sign handling, scale restoration, fused complex squaring, and all algorithm-dispatch boundaries.

A backend is rejected if it is fast but changes accepted escape decisions or fails deep-zoom/root cases without a clear, correctable precision-policy explanation.

---

## 11.4 Native performance protocol

For every backend and precision:

- reuse allocated values where the API permits;
- separately benchmark allocation-heavy naïve code and optimized in-place code;
- measure direct complex and float-pair complex implementations;
- report time per reference iteration;
- report allocations per iteration;
- report memory per stored checkpoint;
- benchmark single thread;
- benchmark realistic complete reference generation;
- benchmark independent-reference parallelism;
- benchmark operation-level parallelism at increasing precision;
- record the crossover precision where threading becomes beneficial;
- measure interaction with the outer render worker pool;
- measure compile time and release binary size.

Do not compare one backend’s optimized in-place path against another backend’s naïve operator-overload path.


### Custom-backend threading matrix

For the custom backend, report a matrix over:

- precision/limb count;
- multiplication algorithm;
- number of arithmetic workers;
- number of outer render workers;
- batch size;
- reference iteration count;
- wall-clock time;
- CPU efficiency;
- peak memory;
- synchronization time.

At minimum compare:

1. one reference, one arithmetic thread;
2. multiple independent references, one thread each;
3. one reference with parallel large multiplication;
4. multiple references with nested parallel multiplication;
5. a globally budgeted hybrid.

The expected default is independent-job parallelism. Nested operation-level parallelism becomes the default only where measured results show a clear benefit without oversubscription.


---

## 11.5 WASM protocol

For every pure-Rust candidate:

- `cargo check --target wasm32-unknown-unknown`;
- compile an executable benchmark harness;
- run arithmetic and reference kernels in a browser;
- measure Chrome/Edge and Firefox where practical;
- record WASM size before and after `wasm-opt`;
- record initialization time;
- record JS/WASM boundary overhead separately;
- test single-threaded baseline first;
- test worker-based execution later;
- test parse/format behaviour;
- run the same correctness fixtures.

`rug` should also receive a feasibility investigation, but no substantial porting effort should be spent unless an existing maintainable MPFR/MPC WASM path is identified.

---

## 11.6 Decision matrix

Score each backend from 1–5, with evidence linked from benchmark output.

| Criterion | Weight |
|---|---:|
| Numerical correctness and convergence | 25% |
| Native reference/root performance | 25% |
| WASM build and browser performance | 20% |
| Allocation and memory behaviour | 10% |
| API fit and ability to reuse temporaries | 8% |
| Maintenance/activity/documentation | 7% |
| Compile time and binary/WASM size | 5% |

For the custom backend, add a separate maintenance gate covering unsafe-code size, algorithm count, documentation quality, cross-platform coverage, and required specialist knowledge. Performance alone cannot override a failed maintenance gate.

### Possible outcomes

1. **One backend everywhere**  
   Adopt one pure-Rust backend if it is accurate and close enough to `rug` natively.

2. **Native `rug`, web pure Rust**  
   Keep `rug` for desktop and use the best pure-Rust backend for WASM behind a coarse backend abstraction.

3. **Native `rug`, web worker/server fallback**  
   If pure-Rust browser arithmetic is inadequate, retain native rendering and provide a restricted web precision range or remote render path. This is a last resort, not the initial design.

4. **Custom specialized reference backend plus general root backend**  
   Use `fractal-arbitrary` for reference-orbit generation and retain `rug` or a pure-Rust general backend for parsing, division, and root finding.

5. **Custom backend everywhere it is justified**  
   Expand the custom backend only if later operations independently meet correctness, performance, and maintenance gates.

6. **Tiered precision**  
   Use fixed precision for moderate cases and arbitrary precision only when required.

The benchmark report must recommend one outcome and list the rejected alternatives with evidence.

---

## 12. LOC reduction and clarity targets

Measure production LOC before and after with generated files, tests, and benchmarks excluded.

Targets:

- remove all dead/commented-out implementation blocks;
- remove obsolete compatibility code;
- remove duplicate mode/config conversion logic;
- eliminate the monolithic `FractalRenderer` and mixed `DataExport`;
- reduce production LOC by **at least 20%**, unless the final report demonstrates that additional explicit error handling or portability code provides greater value;
- keep test and benchmark LOC unconstrained;
- no source file should remain large merely because unrelated responsibilities were left together;
- avoid excessive file splitting: cohesive numerical kernels may remain substantial;
- document public APIs and non-obvious numerical invariants;
- prefer straightforward loops over clever abstractions;
- use macros only where they remove verified repetition without obscuring control flow.

---

## 13. Performance and quality acceptance criteria

The core refactor is complete when:

- [ ] all baseline fixtures pass;
- [ ] root/nucleus residual tests pass;
- [ ] perturbation agrees with direct high-precision orbits within documented tolerances;
- [ ] series approximation agrees with the non-approximated path;
- [ ] deterministic results are obtained across worker counts and chunk sizes;
- [ ] no per-pixel heap allocation remains in hot loops;
- [ ] the deterministic core contains no `Arc<Mutex<FractalRenderer>>` or UI callbacks;
- [ ] `fractal-core` has no dependency on `config`, `image`, `exr`, native filesystem APIs, or native threads;
- [ ] the WASM-capable core has no C FFI;
- [ ] invalid public inputs return typed errors rather than panic;
- [ ] `cargo check --target wasm32-unknown-unknown` passes for the core and selected web numerical backend;
- [ ] end-to-end native rendering is not slower than the baseline outside the 3% noise/regression gate;
- [ ] important workloads have Criterion baselines;
- [ ] allocation and peak-memory results are recorded;
- [ ] the arbitrary-precision comparison report is committed;
- [ ] the custom `fractal-arbitrary` prototype has been implemented or explicitly rejected after kernel characterization;
- [ ] portable, architecture-specific, and multithreaded custom variants have been compared where applicable;
- [ ] operation-level multithreading has measured thresholds and cannot oversubscribe the outer executor;
- [ ] the selected backend strategy is justified by measured results;
- [ ] public documentation includes a minimal render example;
- [ ] the Bevy frontend can consume the core through request/event/tile APIs without accessing internal renderer state.

---

## 14. Recommended pull-request sequence

1. **Baseline fixtures and benchmark harness**
2. **Toolchain/dependency modernization**
3. **Typed domain settings and errors**
4. **Portable extended arithmetic**
5. **Raw frame / colouring / export split**
6. **Reference generation extraction**
7. **Series approximation extraction**
8. **Perturbation and glitch pipeline extraction**
9. **Deterministic render session**
10. **Native executor and progressive tile events**
11. **Allocation and hot-loop optimization**
12. **Arbitrary backend adapters**
13. **Custom fixed-point reference prototype**
14. **Custom kernel correctness and fuzzing**
15. **Custom native optimization and multiplication thresholds**
16. **Custom multithreading experiment**
17. **Backend benchmark implementation**
18. **WASM benchmark harness**
19. **Backend decision and cleanup**
20. **Final documentation and Bevy-ready release**

Each PR should be independently reviewable and keep the regression suite green.

---

## 15. Deliverables

- refactored Cargo workspace;
- deterministic, UI-independent render API;
- native progressive executor;
- typed configuration and errors;
- separated raw frame, colour pipeline, and exporters;
- comprehensive unit/property/differential/regression tests;
- fuzz targets;
- Criterion and instruction/cache benchmark suites;
- allocation and memory profiles;
- arbitrary-precision backend adapters;
- portable `fractal-arbitrary` reference prototype;
- optional architecture-specific custom kernels;
- arithmetic thread-budget and threshold results;
- native and browser benchmark results;
- `ARBITRARY_PRECISION_REPORT.md`;
- migration notes from the old public API;
- a tagged core release suitable for the Bevy frontend.

---

## 16. Reference material snapshot

Versions and APIs must be rechecked at implementation time.

- Rug documentation: <https://docs.rs/rug/latest/rug/>
- Dashu Float documentation: <https://docs.rs/dashu-float/latest/dashu_float/>
- Astro Float documentation: <https://docs.rs/astro-float/latest/astro_float/>
- Malachite Float documentation: <https://docs.rs/malachite-float/latest/malachite_float/>
- Criterion.rs guide: <https://bheisler.github.io/criterion.rs/book/>
- Iai-Callgrind documentation: <https://docs.rs/iai-callgrind/latest/iai_callgrind/>
- Proptest book: <https://proptest-rs.github.io/proptest/proptest/>
- Rust Fuzz Book: <https://rust-fuzz.github.io/book/>
- GMP multiplication algorithms: <https://gmplib.org/manual/Multiplication-Algorithms.html>
- FLINT threading model: <https://flintlib.org/doc/threading.html>
- Rust WASM32 atomics and SIMD: <https://doc.rust-lang.org/core/arch/wasm32/index.html>
- `ibig` pure-Rust big integers: <https://docs.rs/ibig/latest/ibig/>
- `malachite-nz` limb-based integers: <https://docs.rs/malachite-nz/latest/malachite_nz/>
