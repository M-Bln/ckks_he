use bnum::types::I1024;
use ckks::algebra::big_int::BigInt;
use ckks::algebra::complex::{Complex, C64};
use ckks::encoding::Encoder;
use ckks::keys::client_key::to_plaintext;
use ckks::random_distributions::generate_random_vector;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::Rng;

const PLAINTEXT_BOUND: f64 = 4.0;
const Q_EXPONENT: u32 = 30;
const LEVEL_MAX: u32 = 5;

fn bench_encode(c: &mut Criterion) {
    let dimension_exponents = vec![15];
    for &dimension_exponent in &dimension_exponents {
        let q = 1 << Q_EXPONENT;
        let modulus = I1024::from(q).fast_exp(LEVEL_MAX);
        let encoder = Encoder::new(dimension_exponent, modulus, q as f64);
        let real_plaintext = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
        let plaintext = to_plaintext(&real_plaintext);
        c.bench_function(
            &format!("encode dimension_exponent={}", dimension_exponent),
            |b| b.iter(|| encoder.encode(black_box(&plaintext))),
        );
    }
}

fn bench_decode(c: &mut Criterion) {
    let dimension_exponents = vec![15];
    for &dimension_exponent in &dimension_exponents {
        let q = 1 << Q_EXPONENT;
        let modulus = I1024::from(q).fast_exp(LEVEL_MAX);
        let encoder = Encoder::new(dimension_exponent, modulus, q as f64);
        let real_plaintext = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
        let plaintext = to_plaintext(&real_plaintext);
        let ciphertext = encoder.encode(&plaintext);

        c.bench_function(
            &format!("decode dimension_exponent={}", dimension_exponent),
            |b| b.iter(|| encoder.decode(black_box(&ciphertext.clone()))),
        );
    }
}

criterion_group!(benches, bench_encode, bench_decode);
criterion_main!(benches);
