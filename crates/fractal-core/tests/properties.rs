use fractal_core::{
    ComplexCoordinate, FractalFormula, RenderRequest, RenderSettings, Viewport, Zoom,
    compute_frame, pixel_coordinate,
};
use proptest::prelude::*;

proptest! {
    #[test]
    fn bounded_valid_requests_do_not_panic(
        width in 1_u32..33,
        height in 1_u32..33,
        iteration_limit in 1_u32..128,
        center_real in -2.0_f64..2.0,
        center_imag in -2.0_f64..2.0,
        zoom in 0.25_f64..1000.0,
    ) {
        let request = RenderRequest {
            generation: 1,
            viewport: Viewport {
                center: ComplexCoordinate::new(center_real, center_imag),
                zoom: Zoom::new(zoom).unwrap(),
                rotation_radians: 0.25,
                width,
                height,
            },
            settings: RenderSettings {
                fractal: FractalFormula::Mandelbrot { power: 2 },
                iteration_limit,
                ..RenderSettings::default()
            },
        };
        let frame = compute_frame(&request).unwrap();
        prop_assert_eq!(frame.iterations.len(), (width * height) as usize);
        prop_assert!(frame.validate().is_ok());
    }

    #[test]
    fn pixel_mapping_stays_finite_for_bounded_viewports(
        width in 1_u32..4096,
        height in 1_u32..4096,
        x in 0_u32..4095,
        y in 0_u32..4095,
        zoom in 0.25_f64..1000.0,
    ) {
        prop_assume!(x < width && y < height);
        let request = RenderRequest {
            generation: 1,
            viewport: Viewport {
                center: ComplexCoordinate::new(-0.75, 0.1),
                zoom: Zoom::new(zoom).unwrap(),
                rotation_radians: 0.75,
                width,
                height,
            },
            settings: RenderSettings::default(),
        };
        let point = pixel_coordinate(&request, x, y);
        prop_assert!(point.real.is_finite() && point.imag.is_finite());
    }
}
