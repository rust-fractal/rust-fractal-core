use criterion::{Criterion, criterion_group, criterion_main};
use fractal_core::{FractalFormula, OutputChannels, RenderRequest, compute_frame};
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
];

fn request(case: &LocationCase) -> RenderRequest {
    let mut request = request_from_legacy_toml(LEGACY_OPTIONS, case.source, 1)
        .unwrap_or_else(|error| panic!("invalid benchmark location {}: {error}", case.name));
    request.viewport.width = case.width;
    request.viewport.height = case.height;
    request.settings.fractal = case.formula;
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

criterion_group!(benches, benchmark_render);
criterion_main!(benches);
