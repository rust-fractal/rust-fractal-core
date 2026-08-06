# Changelog

Release notes for `rust-fractal-core` are listed below in descending version order.

## Unreleased

- Replaced the legacy native feature and root renderer entry point with a small workspace.
- Added typed render domain settings, deterministic tile/frame computation, raw numerical channels, and regression fixtures.
- Added separate native runtime, I/O, and numerical crates.
- Updated the workspace to Rust 1.97.1 and current direct dependencies.
- Replaced the native C-backed arbitrary-precision experiment with the Windows-compatible Malachite reference implementation.
- Added arbitrary-precision reference benchmarks for finite and deep locations, including `flake`, `wfs`, and `e8000`.

## 0.14.0 - 2020-12-01

- Added progress counters and support for stopping long-running renders.
- Improved iteration-limit handling and distance-estimation colouring.
- Restored the default settings and removed unused code and noisy output.

## 0.13.0 - 2020-11-19

- Added GUI-oriented rendering options and improved GUI interaction support.
- Improved data export, zoom storage, palette offsets, and palette division.
- Fixed low-resolution series-approximation overskipping and location updates.

## 0.12.0 - 2020-10-22

- Added analytic distance-estimation and derivative support.
- Added jitter for moiré-pattern reduction.
- Extended EXR export with distance-estimation and direction channels.
- Fixed rotation handling for additional locations and EXR keyframe export.

## 0.11.0 - 2020-10-17

- Fixed tiled series approximation and made it the default in `high.toml`.
- Improved conservative iteration handling for high-quality renders.

## 0.10.0 - 2020-10-05

- Stabilised probe-based series approximation for rotated images.
- Added configurable data-storage intervals and refreshed the default configs.
- Fixed probe regressions and the test runner script.

## 0.9.1 - 2020-10-02

- Fixed regressions involving high series-approximation term counts and low iteration counts.

## 0.9.0 - 2020-10-02

- Added experimental tiled series approximation and probe improvements.
- Added configurable glitch tolerances and experimental options.
- Fixed data-interval handling and rare regressions in glitch-reference creation.

## 0.8.0 - 2020-09-29

- Added configurable colour palettes and integrated palette colouring.
- Added experimental series-approximation options and new example locations.
- Fixed extended-float range handling and refactored data-export initialisation.

## 0.7.2 - 2020-09-18

- Increased series-approximation terms at low zoom levels to avoid rare errors.
- Added default quality files and fixed rotation settings in location files.

## 0.7.1 - 2020-09-15

- Optimised probe perturbation for series approximation.
- Stopped tracking generated segment data.

## 0.7.0 - 2020-09-11

- Added probe sampling and automatic reduction of series-approximation terms.
- Added frame offsets, centre-removal controls, and improved EXR segment data.
- Added export-type options and automatic maximum-keyframe checks.

## 0.6.0 - 2020-09-05

- Added KFR/KFB output support and improved video assembly support.
- Added central-part removal and automatic iteration experiments.
- Improved memory use and optimised the core perturbation loop.

## 0.5.0 - 2020-08-27

- Added rotation parameters and more location examples.
- Reorganised dependencies and expanded export options.

## 0.4.0 - 2020-08-26

- Added EXR iteration-data export and reusable references for image sequences.
- Improved colouring, image output, zoom-out configuration, and multithreading.
- Added more location and palette examples.

## 0.3.0 - 2020-08-23

- Added automatic series-approximation term calculation.
- Improved series-approximation packing, optimisation, and glitch correction.

## 0.2.0 - 2020-08-21

- Reworked the perturbation algorithm for deeper zoom levels.
- Improved series approximation reliability and performance.
- Simplified dependencies and improved GNU/Linux compilation support.

## 0.1.0 - 2020-05-05

- Initial public release of the Mandelbrot renderer.
- Added perturbation rendering, series approximation, multithreading, smooth and histogram colouring, palettes, supersampling, and configurable zoom locations.
