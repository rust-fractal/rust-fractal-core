use fractal_core::{
    ComplexCoordinate, FractalFormula, OutputChannels, RenderRequest, RenderSettings, Viewport,
    Zoom, compute_frame,
};

#[test]
fn four_by_four_mandelbrot_fixture_is_stable() {
    let request = RenderRequest {
        generation: 1,
        viewport: Viewport {
            center: ComplexCoordinate::new(0.0, 0.0),
            zoom: Zoom::unit(),
            rotation_radians: 0.0,
            width: 4,
            height: 4,
        },
        settings: RenderSettings {
            fractal: FractalFormula::Mandelbrot { power: 2 },
            iteration_limit: 32,
            ..RenderSettings::default()
        },
    };

    let expected = include_str!("fixtures/mandelbrot_4x4.iterations")
        .split_whitespace()
        .map(|value| value.parse::<u32>().unwrap())
        .collect::<Vec<_>>();
    let actual = compute_frame(&request).unwrap();
    assert_eq!(actual.iterations, expected);
    assert_eq!(actual.raw_hash().unwrap(), 783574553016584100);
}

#[test]
fn optional_channels_are_aligned_with_iteration_data() {
    let request = RenderRequest {
        generation: 2,
        viewport: Viewport {
            center: ComplexCoordinate::new(-0.75, 0.0),
            zoom: Zoom::unit(),
            rotation_radians: 0.0,
            width: 8,
            height: 8,
        },
        settings: RenderSettings {
            output: OutputChannels {
                smooth_iteration: true,
                distance: true,
                ..OutputChannels::default()
            },
            ..RenderSettings::default()
        },
    };

    let result = compute_frame(&request).unwrap();
    assert_eq!(result.iterations.len(), 64);
    assert_eq!(result.smooth_iteration.as_ref().unwrap().len(), 64);
    assert_eq!(result.distance.as_ref().unwrap().len(), 64);
    assert!(result.smooth_iteration.unwrap().into_iter().all(f64::is_finite));
    assert!(result.distance.unwrap().into_iter().all(f64::is_finite));
}

#[test]
fn approximation_path_matches_direct_path_for_a_moderate_zoom() {
    let direct = RenderRequest {
        generation: 3,
        viewport: Viewport {
            center: ComplexCoordinate::new(-0.75, 0.0),
            zoom: Zoom::new(8.0).unwrap(),
            rotation_radians: 0.0,
            width: 16,
            height: 12,
        },
        settings: RenderSettings { iteration_limit: 128, ..RenderSettings::default() },
    };
    let mut approximated = direct.clone();
    approximated.settings.approximation.enabled = true;
    approximated.settings.approximation.order = 8;
    let direct = compute_frame(&direct).unwrap();
    let approximated = compute_frame(&approximated).unwrap();
    assert_eq!(approximated.iterations, direct.iterations);
}
