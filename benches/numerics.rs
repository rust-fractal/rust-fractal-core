use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use fractal_arbitrary::MalachiteReference;
use fractal_numerics::{Extended, malachite::MalachiteComplex};

fn benchmark_numerics(c: &mut Criterion) {
    let value = Extended::from_f64(1.2345e200).unwrap();
    c.bench_function("extended_multiply", |b| b.iter(|| value.checked_mul(value).unwrap()));

    let constant = MalachiteComplex::from_f64(-0.75, 0.0, 256).unwrap();
    c.bench_function("malachite_complex_square_add", |b| {
        b.iter_batched(
            || MalachiteComplex::from_f64(0.25, -0.5, 256).unwrap(),
            |mut complex| {
                complex.square_add(&constant).unwrap();
                std::hint::black_box(complex.to_f64());
            },
            BatchSize::SmallInput,
        )
    });

    c.bench_function("malachite_reference_100_iterations", |b| {
        b.iter(|| MalachiteReference::generate("-0.75", "0", 256, 100).unwrap())
    });
}

criterion_group!(benches, benchmark_numerics);
criterion_main!(benches);
