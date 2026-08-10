# Performance report

These measurements were taken on the local Windows runner and are intended to
be reproducible, not portable performance guarantees.

## Environment

| Item | Value |
| --- | --- |
| CPU | AMD Ryzen 7 5800X3D, 8 physical cores / 16 logical processors |
| OS/target | Windows 11 Core, `x86_64-pc-windows-msvc` |
| Rust | `rustc 1.97.1`, LLVM 22.1.6; Cargo 1.97.1 |
| Criterion | 0.8.2 |
| Malachite | 0.10.0, `std`, `floats`, `naturals_and_integers`, `rationals` (the last feature is required by this registry release’s `Float` re-export) |
| Rayon | 1.12.0 |
| Other relevant locked crates | serde 1.0.229, image 0.25.10, exr 1.74.2, toml 1.1.4 |

Relevant commands:

```text
cargo bench --bench core_render --all-features --locked -- --sample-size 10 --warm-up-time 0.4 --measurement-time 0.4
cargo bench --bench arbitrary_reference --all-features --locked -- --sample-size 10 --warm-up-time 0.5 --measurement-time 0.5
cargo bench --bench numerics --all-features --locked -- --sample-size 10 --warm-up-time 0.5 --measurement-time 0.5
cargo bench --bench runtime_scaling --all-features --locked -- --sample-size 10 --warm-up-time 0.5 --measurement-time 0.5
cargo bench --bench reference_allocations --all-features --locked
```

Criterion results below are medians from ten samples. The benchmark binaries
were built with the repository release profile (`debug = true`); no
architecture-specific `target-cpu` setting was used.

## End-to-end render comparison

The PR baseline is a clean run of commit `7a539bf` in a temporary worktree.
`master` has a different renderer/API and no equivalent Criterion target, so a
same-input number is not meaningful and is shown as unavailable. Final cases
use the actual exact-coordinate production pipeline; the baseline did not have
the added deep cases.

| Case | `master` | PR #45 before | final production | Notes |
| --- | ---: | ---: | ---: | --- |
| e14, 64×64 | n/a | 367.75 ms | 34.01 ms | comparable shallow/deep render |
| e50, 64×64 | n/a | 301.72 ms | 39.37 ms | final uses exact center-relative offsets |
| flake, 64×64 | n/a | 3.18 ms | 32.74 ms | baseline was an invalid collapsed-`f64` deep render; no speedup claim |
| Burning Ship, 128×128 | n/a | 4.07 ms | 3.77 ms | comparable native scalar path |
| cubic, 64×64 | n/a | 43.73 ms | 8.41 ms | specialized power-3 perturbation |
| e1000, 16×16 | n/a | n/a | 2.92 ms | new real deep production case |
| e1086, 16×16 | n/a | n/a | 3.15 ms | new real deep production case |
| e1121, 16×16 | n/a | n/a | 3.01 ms | new real deep production case |
| e8000, 8×8 | n/a | n/a | 13.22 ms | new real deep production case |

The Criterion “change” column is not used for these claims because its
persistent comparison data was changed by the refactor. The PR baseline was
measured independently with the same ten-sample settings.

## Production reference generation

`arbitrary_reference` now calls `fractal_core::generate_reference_exact` with
the raw decimal locations. Precision is selected by the production
`PrecisionPolicy`; it is not an `f64` approximation.

| Location | Precision bits | Iterations | Final median |
| --- | ---: | ---: | ---: |
| e14 | 143 | 512 | 367 µs |
| e50 | 263 | 512 | 657 µs |
| flake | 619 | 512 | 1.06 ms |
| wfs | 5,149 | 512 | 8.94 ms |
| e1000 | 3,418 | 256 | 1.68 ms |
| e1086 | 3,704 | 256 | 3.88 ms |
| e1121 | 3,822 | 256 | 2.69 ms |
| e8000 | 26,673 | 256 | 46.44 ms |

The initial PR’s separate wrapper benchmark reported approximately 526 µs,
556 µs, 760 µs, 7.30 ms, 2.00 ms, 1.83 ms, 1.79 ms, and 29.4 ms for those
locations respectively, but used a different precision policy and retained
dense Malachite points. Those are historical reference numbers, not a direct
speedup comparison.

## Series variants

Measured on e50, 16×16, 128 iterations. The disabled case performs direct
extended perturbation; enabled cases use the configured order.

| Series mode | Total render | Series generation |
| --- | ---: | ---: |
| disabled | 20.48 ms | — |
| order 4 | 737 µs | 15.8 µs |
| order 8 | 1.74 ms | 65.9 µs |
| order 16 | 2.25 ms | 131 µs |
| order 32 | 4.29 ms | 735 µs |
| order 64 | 11.25 ms | 2.75 ms |

For references longer than 256 iterations, coefficient states are stored at
interval 32; shorter references retain interval 1. The current benchmark
records total work and generation cost, but does not retain a per-pixel average
skip counter in the raw output API. The large disabled/enabled difference is
the end-to-end evidence used for the default series path.

## Runtime scaling

Measured on a 192×128, 256-iteration render. Progressive delivery uses the
bounded coordinator; full-frame uses one preallocated channel set and disjoint
row views.

| Workers | Progressive | Speedup | Full frame | Speedup |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 4.09 ms | 1.00× | 3.70 ms | 1.00× |
| 2 | 2.39 ms | 1.71× | 2.25 ms | 1.65× |
| 4 | 1.59 ms | 2.57× | 1.64 ms | 2.26× |
| 8 (physical) | 1.26 ms | 3.24× | 1.59 ms | 2.34× |
| 16 (logical) | 1.20 ms | 3.41× | 1.59 ms | 2.34× |

The progressive path scales better at 8–16 workers for this small workload;
the full-frame path wins at one to four workers by avoiding tile copies. The
runtime does not create a second pool per render.

## Allocation probe

The standalone `reference_allocations` binary wraps the system allocator and
reports allocation/reallocation events and total bytes requested during one
reference generation. `allocated_bytes` is cumulative allocation traffic, not
peak live memory.

| Location | Precision | Alloc events | Allocated bytes | Wall time |
| --- | ---: | ---: | ---: | ---: |
| e14 | 143 | 3,618 | 134 KiB | 523 µs |
| e50 | 263 | 7,708 | 450 KiB | 809 µs |
| flake | 619 | 7,708 | 840 KiB | 1.20 ms |
| wfs | 5,149 | 7,718 | 10.6 MiB | 5.93 ms |
| e1000 | 3,418 | 3,866 | 3.82 MiB | 1.92 ms |
| e1086 | 3,704 | 6,420 | 5.14 MiB | 2.55 ms |
| e1121 | 3,822 | 3,870 | 4.04 MiB | 2.17 ms |
| e8000 | 26,673 | 7,766 | 23.6 MiB | 35.2 ms |

Production pixels do not allocate arbitrary-precision values. The remaining
reference traffic is dominated by Malachite temporaries and the compact/f64
orbit vectors; sparse checkpoints avoid retaining a Malachite clone for every
iteration. Peak RSS was not measured.

## Numerical kernel experiments

The scalar extended multiply benchmark was approximately 7.95 ns. Current
Malachite square-add versus difference-of-squares medians were:

| Precision | Square-add | Difference-of-squares |
| ---: | ---: | ---: |
| 256 bits | 1.23 µs | 1.02 µs |
| 512 bits | 1.56 µs | 1.24 µs |
| 1,024 bits | 2.16 µs | 1.86 µs |
| 2,048 bits | 1.43 µs | 1.25 µs |
| 4,096 bits | 1.77 µs | 1.71 µs |
| 8,192 bits | 2.49 µs | 2.17 µs |

An intra-iteration Rayon prototype was tested at 256 through 65,536 bits and
1, 2, 4, 8, and 16 workers using a pool created outside the measured loop.
It did not beat scalar arithmetic: at 16,384 bits, for example, the scalar
product was about 2.52 µs while the two-worker prototype was about 3.00 µs;
at 65,536 bits it was about 10.75 µs versus 14.17 µs. No production threshold
was therefore selected—the effective threshold is “disabled”—and the
prototype/dependency was removed.

## Removed or unsuccessful experiments

- Dense arbitrary point storage was replaced by compact points and sparse
  checkpoints.
- Direct absolute-`f64` deep fallback was removed in favor of typed reference
  errors or extended rebase correction.
- Intra-iteration Rayon reference arithmetic was removed after the measured
  overhead result above.
- No SIMD, Rug/MPFR backend, or project-owned fixed-point backend was kept;
  representative end-to-end evidence did not justify their portability and
  maintenance costs.
- Release LTO/codegen-unit/`target-cpu=native` tuning was not imposed on the
  portable artifact. The benchmark profile remains portable and keeps debug
  symbols for profiling.
