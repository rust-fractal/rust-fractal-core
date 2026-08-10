use criterion::{Criterion, criterion_group, criterion_main};
use fractal_core::{
    FractalFormula, OutputChannels, PrecisionPolicy, RenderRequest, compute_frame,
    generate_reference_exact, generate_series,
};
use fractal_io::request_from_legacy_toml;

const LEGACY_OPTIONS: &str = include_str!("../default.toml");

struct LocationCase {
    name: &'static str,
    source: &'static str,
    formula: FractalFormula,
    width: u32,
    height: u32,
    iteration_cap: u32,
}

const LOCATION_CASES: &[LocationCase] = &[
    LocationCase {
        name: "e14_mandelbrot",
        source: include_str!("../locations/e14.toml"),
        formula: FractalFormula::Mandelbrot { power: 2 },
        width: 64,
        height: 64,
        iteration_cap: 512,
    },
    LocationCase {
        name: "burning_ship_test",
        source: include_str!("../locations/burning_ship_test.toml"),
        formula: FractalFormula::BurningShip { power: 2 },
        width: 128,
        height: 128,
        iteration_cap: 256,
    },
    LocationCase {
        name: "e50_mandelbrot",
        source: include_str!("../locations/e50.toml"),
        formula: FractalFormula::Mandelbrot { power: 2 },
        width: 64,
        height: 64,
        iteration_cap: 512,
    },
    LocationCase {
        name: "flake",
        source: include_str!("../locations/flake.toml"),
        formula: FractalFormula::Mandelbrot { power: 2 },
        width: 64,
        height: 64,
        iteration_cap: 256,
    },
    LocationCase {
        name: "cubic_test",
        source: include_str!("../locations/cubic_test.toml"),
        formula: FractalFormula::Mandelbrot { power: 3 },
        width: 64,
        height: 64,
        iteration_cap: 256,
    },
    LocationCase {
        name: "e1000_production_deep",
        source: include_str!("../locations/e1000.toml"),
        formula: FractalFormula::Mandelbrot { power: 2 },
        width: 16,
        height: 16,
        iteration_cap: 128,
    },
    LocationCase {
        name: "e1086_production_deep",
        source: include_str!("../locations/e1086.toml"),
        formula: FractalFormula::Mandelbrot { power: 2 },
        width: 16,
        height: 16,
        iteration_cap: 128,
    },
    LocationCase {
        name: "e1121_production_deep",
        source: include_str!("../locations/e1121.toml"),
        formula: FractalFormula::Mandelbrot { power: 2 },
        width: 16,
        height: 16,
        iteration_cap: 128,
    },
    LocationCase {
        name: "e8000_production_deep",
        source: include_str!("../locations/e8000.toml"),
        formula: FractalFormula::Mandelbrot { power: 2 },
        width: 8,
        height: 8,
        iteration_cap: 64,
    },
];

fn request(case: &LocationCase) -> RenderRequest {
    let mut request = request_from_legacy_toml(LEGACY_OPTIONS, case.source, 1)
        .unwrap_or_else(|error| panic!("invalid benchmark location {}: {error}", case.name));
    request.viewport.width = case.width;
    request.viewport.height = case.height;
    request.settings.fractal = case.formula;
    if !matches!(case.formula, FractalFormula::Mandelbrot { power: 2 }) {
        request.settings.approximation.enabled = false;
        request.settings.approximation.order = 0;
    }
    request.settings.iteration_limit = request.settings.iteration_limit.min(case.iteration_cap);
    request.settings.output =
        OutputChannels { smooth_iteration: true, distance: true, ..request.settings.output };
    request.validate().expect("benchmark request should be valid");
    request
}

fn benchmark_render(c: &mut Criterion) {
    let mut group = c.benchmark_group("location_render");
    for case in LOCATION_CASES {
        let request = request(case);
        group.bench_function(case.name, |b| b.iter(|| compute_frame(&request).unwrap()));
    }
    group.finish();
}

fn benchmark_series(c: &mut Criterion) {
    let mut base = request(&LOCATION_CASES[2]);
    base.viewport.width = 16;
    base.viewport.height = 16;
    base.settings.iteration_limit = 128;

    let mut render_group = c.benchmark_group("series_render_variants");
    for (name, order) in [
        ("disabled", 0),
        ("order4", 4),
        ("order8", 8),
        ("order16", 16),
        ("order32", 32),
        ("order64", 64),
    ] {
        let mut request = base.clone();
        request.settings.approximation.enabled = order != 0;
        request.settings.approximation.order = order;
        render_group.bench_function(name, |b| {
            b.iter(|| std::hint::black_box(compute_frame(&request).unwrap()))
        });
    }
    render_group.finish();

    let precision = PrecisionPolicy::default().required_bits(base.viewport.zoom.clone()).unwrap();
    let reference = generate_reference_exact(
        &base.viewport.center,
        base.settings.fractal,
        precision,
        base.settings.iteration_limit,
    )
    .unwrap();
    let mut generation_group = c.benchmark_group("series_generation_variants");
    for order in [4, 8, 16, 32, 64] {
        generation_group.bench_function(format!("order{order}"), |b| {
            b.iter(|| std::hint::black_box(generate_series(&reference, order).unwrap()))
        });
    }
    generation_group.finish();
}

criterion_group!(benches, benchmark_render, benchmark_series);
criterion_main!(benches);
