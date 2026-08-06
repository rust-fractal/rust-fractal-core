use criterion::{Criterion, criterion_group, criterion_main};
use fractal_arbitrary::MalachiteReference;
use std::hint::black_box;

struct LocationCase {
    name: &'static str,
    source: &'static str,
    iteration_cap: u32,
}

const LOCATION_CASES: &[LocationCase] = &[
    LocationCase { name: "e14", source: include_str!("../locations/e14.toml"), iteration_cap: 512 },
    LocationCase { name: "e50", source: include_str!("../locations/e50.toml"), iteration_cap: 512 },
    LocationCase {
        name: "flake",
        source: include_str!("../locations/flake.toml"),
        iteration_cap: 512,
    },
    LocationCase { name: "wfs", source: include_str!("../locations/wfs.toml"), iteration_cap: 512 },
    LocationCase {
        name: "e8000",
        source: include_str!("../locations/e8000.toml"),
        iteration_cap: 256,
    },
    LocationCase {
        name: "e1086",
        source: include_str!("../locations/e1086.toml"),
        iteration_cap: 256,
    },
    LocationCase {
        name: "e1121",
        source: include_str!("../locations/e1121.toml"),
        iteration_cap: 256,
    },
    LocationCase {
        name: "e1000",
        source: include_str!("../locations/e1000.toml"),
        iteration_cap: 256,
    },
];

fn field<'a>(source: &'a str, name: &str) -> &'a str {
    source
        .lines()
        .find_map(|line| {
            let (key, value) = line.split_once('=')?;
            (key.trim() == name).then(|| value.trim().trim_matches('"'))
        })
        .unwrap_or_else(|| panic!("location is missing {name}"))
}

fn arbitrary_reference(source: &str, iteration_cap: u32) -> (u32, u32) {
    let zoom = field(source, "zoom");
    let reference = MalachiteReference::generate_for_zoom(
        field(source, "real"),
        field(source, "imag"),
        zoom,
        iteration_cap,
    )
    .unwrap();
    (reference.statistics.precision_bits, reference.statistics.escaped_at.unwrap_or(0))
}

fn benchmark_arbitrary_references(c: &mut Criterion) {
    let mut group = c.benchmark_group("arbitrary_reference_locations");
    for case in LOCATION_CASES {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(arbitrary_reference(case.source, case.iteration_cap)))
        });
    }
    group.finish();
}

criterion_group!(benches, benchmark_arbitrary_references);
criterion_main!(benches);
