# Rust support policy

The workspace is developed and tested with Rust 1.97.1, pinned in
`rust-toolchain.toml`. The manifests declare Rust 1.97 as the minimum supported
version for the current release line. Dependency updates must continue to pass
the pinned toolchain before they are merged.

The current CI target is native Windows. WebAssembly checks are intentionally
deferred until the frontend/backend decision is made.
