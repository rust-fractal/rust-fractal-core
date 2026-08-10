use fractal_core::{
    ApproximationSettings, ComplexCoordinate, FractalFormula, RenderRequest, RenderSettings,
    Viewport, Zoom, compute_frame,
};
use serde::Deserialize;

#[derive(Deserialize)]
struct FixtureManifest {
    version: u32,
    hash_algorithm: String,
    fixture: Vec<FixtureSpec>,
}

#[derive(Deserialize)]
struct FixtureSpec {
    name: String,
    center: [f64; 2],
    zoom: f64,
    rotation_radians: f64,
    width: u32,
    height: u32,
    formula: String,
    power: u32,
    iteration_limit: u32,
    approximation: bool,
    approximation_order: Option<u32>,
    expected_iterations: Option<Vec<u32>>,
    expected_raw_hash: u64,
}

#[test]
fn canonical_manifest_matches_current_raw_outputs() {
    let manifest: FixtureManifest = toml::from_str(include_str!("fixtures/manifest.toml")).unwrap();
    assert_eq!(manifest.version, 1);
    assert!(manifest.hash_algorithm.contains("fnv1a64"));
    assert_eq!(manifest.fixture.len(), 5);

    for fixture in manifest.fixture {
        let fractal = match fixture.formula.as_str() {
            "mandelbrot" => FractalFormula::Mandelbrot { power: fixture.power },
            "burning_ship" => FractalFormula::BurningShip { power: fixture.power },
            other => panic!("unknown fixture formula: {other}"),
        };
        let request = RenderRequest {
            generation: 1,
            viewport: Viewport {
                center: ComplexCoordinate::new(fixture.center[0], fixture.center[1]),
                zoom: Zoom::new(fixture.zoom).unwrap(),
                rotation_radians: fixture.rotation_radians,
                width: fixture.width,
                height: fixture.height,
            },
            settings: RenderSettings {
                fractal,
                iteration_limit: fixture.iteration_limit,
                approximation: ApproximationSettings {
                    enabled: fixture.approximation,
                    order: fixture.approximation_order.unwrap_or(0),
                },
                ..RenderSettings::default()
            },
        };
        let result = compute_frame(&request).unwrap();
        assert_eq!(result.raw_hash().unwrap(), fixture.expected_raw_hash, "{}", fixture.name);
        if let Some(expected) = fixture.expected_iterations {
            assert_eq!(result.iterations, expected, "{}", fixture.name);
        }
    }
}
