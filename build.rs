extern crate ispc;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    let mut cfg = ispc::Config::new();

    if cfg!(windows) {
        cfg.debug(false);
    }

    let ispc_files = vec!["src/math/perturbation_double.ispc"];

    for s in &ispc_files[..] {
        cfg.file(*s);
    }

    cfg.target_isas(vec![
        ispc::opt::TargetISA::SSE2i32x4,
        ispc::opt::TargetISA::SSE4i32x4,
        ispc::opt::TargetISA::AVX1i32x8,
        ispc::opt::TargetISA::AVX2i32x8,
        ispc::opt::TargetISA::AVX512KNLi32x16,
        ispc::opt::TargetISA::AVX512SKXi32x16
    ]);

    cfg.compile("mandelbrot");
}