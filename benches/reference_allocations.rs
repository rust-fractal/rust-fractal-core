use fractal_core::{
    ComplexCoordinate, FractalFormula, PrecisionPolicy, Zoom, generate_reference_exact,
};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

struct CountingAllocator;

static ALLOCATIONS: AtomicU64 = AtomicU64::new(0);
static ALLOCATED_BYTES: AtomicU64 = AtomicU64::new(0);

#[global_allocator]
static GLOBAL: CountingAllocator = CountingAllocator;

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
            ALLOCATED_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        unsafe { System.dealloc(pointer, layout) };
    }

    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let resized = unsafe { System.realloc(pointer, layout, new_size) };
        if !resized.is_null() {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
            ALLOCATED_BYTES.fetch_add(new_size as u64, Ordering::Relaxed);
        }
        resized
    }
}

const LOCATIONS: &[(&str, &str, u32)] = &[
    ("e14", include_str!("../locations/e14.toml"), 512),
    ("e50", include_str!("../locations/e50.toml"), 512),
    ("flake", include_str!("../locations/flake.toml"), 512),
    ("wfs", include_str!("../locations/wfs.toml"), 512),
    ("e8000", include_str!("../locations/e8000.toml"), 256),
    ("e1086", include_str!("../locations/e1086.toml"), 256),
    ("e1121", include_str!("../locations/e1121.toml"), 256),
    ("e1000", include_str!("../locations/e1000.toml"), 256),
];

fn field<'a>(source: &'a str, name: &str) -> &'a str {
    source
        .lines()
        .find_map(|line| {
            let (key, value) = line.split_once('=')?;
            (key.trim() == name).then(|| value.trim().trim_matches('"'))
        })
        .unwrap_or("")
}

fn main() {
    for &(name, source, iteration_limit) in LOCATIONS {
        let zoom = Zoom::from_decimal(field(source, "zoom")).unwrap();
        let center =
            ComplexCoordinate::from_decimal(field(source, "real"), field(source, "imag")).unwrap();
        let precision_bits = PrecisionPolicy::default().required_bits(zoom).unwrap();
        ALLOCATIONS.store(0, Ordering::Relaxed);
        ALLOCATED_BYTES.store(0, Ordering::Relaxed);
        let start = Instant::now();
        let reference = generate_reference_exact(
            &center,
            FractalFormula::Mandelbrot { power: 2 },
            precision_bits,
            iteration_limit,
        )
        .unwrap();
        let elapsed = start.elapsed();
        println!(
            "location={name} precision_bits={precision_bits} generated={} elapsed_us={} allocations={} allocated_bytes={}",
            reference.statistics.generated_iterations,
            elapsed.as_micros(),
            ALLOCATIONS.load(Ordering::Relaxed),
            ALLOCATED_BYTES.load(Ordering::Relaxed),
        );
    }
}
