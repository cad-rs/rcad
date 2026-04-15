//! IGES I/O benchmark stub.

use criterion::{criterion_group, criterion_main, Criterion};

fn stub_benchmark(_c: &mut Criterion) {
    // Placeholder benchmark
}

criterion_group!(benches, stub_benchmark);
criterion_main!(benches);
