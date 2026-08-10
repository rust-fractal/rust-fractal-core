use criterion::{Criterion, criterion_group, criterion_main};
use fractal_core::{
    ComplexCoordinate, OutputChannels, RenderRequest, RenderSettings, Viewport, Zoom,
};
use fractal_runtime::{CancellationToken, RenderExecutor};
use std::hint::black_box;

fn request() -> RenderRequest {
    RenderRequest {
        generation: 1,
        viewport: Viewport {
            center: ComplexCoordinate::new(-0.75, 0.0),
            zoom: Zoom::unit(),
            rotation_radians: 0.0,
            width: 192,
            height: 128,
        },
        settings: RenderSettings {
            iteration_limit: 256,
            output: OutputChannels {
                smooth_iteration: true,
                distance: true,
                ..OutputChannels::default()
            },
            ..RenderSettings::default()
        },
    }
}

fn benchmark_runtime_scaling(c: &mut Criterion) {
    let request = request();
    let mut group = c.benchmark_group("runtime_scaling");
    let logical = std::thread::available_parallelism().map_or(1, |count| count.get());
    let mut workers = vec![1, 2, 4, 8, logical];
    workers.sort_unstable();
    workers.dedup();
    for worker_count in workers {
        let executor = RenderExecutor::with_workers(32, 32, worker_count).unwrap();
        group.bench_function(format!("progressive_{worker_count}_workers"), |b| {
            b.iter(|| {
                let token = CancellationToken::for_generation(1);
                black_box(executor.render(&request, &token))
            })
        });
        group.bench_function(format!("full_frame_{worker_count}_workers"), |b| {
            b.iter(|| {
                let token = CancellationToken::for_generation(1);
                black_box(executor.compute_frame_full(&request, &token).unwrap())
            })
        });
    }
    group.finish();
}

criterion_group!(benches, benchmark_runtime_scaling);
criterion_main!(benches);
