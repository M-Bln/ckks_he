use bnum::types::I1024;
use ckks_he::algebra::big_int::BigInt;
use ckks_he::algebra::complex::C64;
use ckks_he::algebra::polynomial::Polynomial;
use ckks_he::encoding::Encoder;
use ckks_he::keys::client_key::{calculate_relative_error, to_plaintext};
use ckks_he::keys::key_generator::generate_pair_keys_default;
use ckks_he::random_distributions::generate_random_vector;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

const PLAINTEXT_BOUND: f64 = 4.0;
const Q_EXPONENT: u32 = 30;
const LEVEL_MAX: u32 = 5;
const EXPONENT: u32 = 2;

fn bench_encode(c: &mut Criterion) {
    let dimension_exponents = vec![EXPONENT];
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
    let dimension_exponents = vec![EXPONENT];
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

fn bench_encrypt_decrypt(c: &mut Criterion) {
    let dimension_exponents = vec![EXPONENT]; // Try smaller exponents first
    for &dimension_exponent in &dimension_exponents {
        // Generate keys
        let (mut client_key, _server_key) =
            generate_pair_keys_default::<I1024>(dimension_exponent, LEVEL_MAX);

        // Generate plaintext
        let real_plaintext = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
        let plaintext = to_plaintext(&real_plaintext);

        // Encrypt the plaintext
        c.bench_function(&format!("encrypt_dim_{}", dimension_exponent), |b| {
            b.iter(|| {
                client_key
                    .encrypt(black_box(&plaintext), PLAINTEXT_BOUND)
                    .unwrap()
            })
        });

        let ciphertext = client_key.encrypt(&plaintext, PLAINTEXT_BOUND).unwrap();

        // Decrypt the ciphertext
        c.bench_function(&format!("decrypt_dim_{}", dimension_exponent), |b| {
            b.iter(|| {
                let _decrypted = client_key.decrypt(black_box(&ciphertext));
            })
        });

        let decrypted = client_key.decrypt(&ciphertext);
        let relative_error = calculate_relative_error(&plaintext, &decrypted);
        println!(
            "Relative error after decryption (dim {}): {}",
            dimension_exponent, relative_error
        );
    }
}

fn bench_encrypt_add_decrypt(c: &mut Criterion) {
    let dimension_exponents = vec![EXPONENT];
    for &dimension_exponent in &dimension_exponents {
        // Generate keys
        let (mut client_key, server_key) =
            generate_pair_keys_default::<I1024>(dimension_exponent, LEVEL_MAX);

        // Generate ciphertexts
        let real_plaintext1 = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
        let plaintext1 = to_plaintext(&real_plaintext1);
        let ciphertext1 = client_key.encrypt(&plaintext1, PLAINTEXT_BOUND).unwrap();

        let real_plaintext2 = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
        let plaintext2 = to_plaintext(&real_plaintext2);
        let ciphertext2 = client_key.encrypt(&plaintext2, PLAINTEXT_BOUND).unwrap();

        c.bench_function(
            &format!("homomorphic add dimension: 2^{}", dimension_exponent),
            |b| {
                b.iter(|| {
                    server_key.add(black_box(&ciphertext1), black_box(&ciphertext2));
                })
            },
        );

        let expected: Vec<C64> = plaintext1
            .iter()
            .zip(plaintext2.iter())
            .map(|(c1, c2)| *c1 + *c2)
            .collect();
        let result = &server_key.add(&ciphertext1, &ciphertext2);
        let decrypted = client_key.decrypt(&result);
        let expected_error = client_key.rescaled_error(&result);
        let relative_error = calculate_relative_error(&expected, &decrypted);
        println!("Expected error: {}", expected_error);
        println!(
            "Relative error after decryption (dim {}): {}",
            dimension_exponent, relative_error
        );
    }
}

fn bench_encrypt_mul_decrypt(c: &mut Criterion) {
    let dimension_exponents = vec![EXPONENT];
    for &dimension_exponent in &dimension_exponents {
        // Generate keys
        let (mut client_key, server_key) =
            generate_pair_keys_default::<I1024>(dimension_exponent, LEVEL_MAX);

        // Generate ciphertexts
        let real_plaintext1 = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
        let plaintext1 = to_plaintext(&real_plaintext1);
        let ciphertext1 = client_key
            .encrypt(&plaintext1, 100.0 * PLAINTEXT_BOUND)
            .unwrap();

        let real_plaintext2 = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
        let plaintext2 = to_plaintext(&real_plaintext2);
        let ciphertext2 = client_key
            .encrypt(&plaintext2, 100.0 * PLAINTEXT_BOUND)
            .unwrap();

        c.bench_function(
            &format!(
                "homomorphic multiplication, dimension: 2^{}",
                dimension_exponent
            ),
            |b| {
                b.iter(|| {
                    server_key
                        .mul(black_box(&ciphertext1), black_box(&ciphertext2))
                        .unwrap();
                })
            },
        );

        let expected: Vec<C64> = plaintext1
            .iter()
            .zip(plaintext2.iter())
            .map(|(c1, c2)| *c1 * c2)
            .collect();
        let result = &server_key.mul(&ciphertext1, &ciphertext2).unwrap();
        let decrypted = client_key.decrypt(&result);
        let expected_error = client_key.rescaled_error(&result);
        let relative_error = calculate_relative_error(&expected, &decrypted);
        println!("Expected error: {}", expected_error);
        println!(
            "Relative error after decryption (dim {}): {}",
            dimension_exponent, relative_error
        );
    }
}

fn bench_apply_polynomial(c: &mut Criterion) {
    let dimension_exponents = vec![EXPONENT];
    for &dimension_exponent in &dimension_exponents {
        // Generate keys
        let (mut client_key, server_key) =
            generate_pair_keys_default::<I1024>(dimension_exponent, LEVEL_MAX);

        // Generate ciphertexts
        let real_plaintext = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
        let plaintext = to_plaintext(&real_plaintext);
        let ciphertext = client_key
            .encrypt(&plaintext, 100.0 * PLAINTEXT_BOUND)
            .unwrap();

        let polynomial = Polynomial::<I1024>::new(vec![
            I1024::from(2),
            I1024::from(1),
            I1024::from(2),
            I1024::from(1),
        ]);
        let complex_polynomial = polynomial.to_c64();
        let expected: Vec<C64> = plaintext
            .iter()
            .map(|c| complex_polynomial.eval(*c))
            .collect();
        let result = server_key
            .apply_polynomial(black_box(&polynomial), black_box(&ciphertext))
            .unwrap();
        let decrypted = client_key.decrypt(&result);
        let expected_error = client_key.rescaled_error(&result);
        let relative_error = calculate_relative_error(&expected, &decrypted);
        c.bench_function(
            &format!("Polynomial evaluation, dimension: 2^{}", dimension_exponent),
            |b| {
                b.iter(|| {
                    server_key
                        .apply_polynomial(black_box(&polynomial), black_box(&ciphertext))
                        .unwrap();
                })
            },
        );

        println!("Expected error: {}", expected_error);
        println!(
            "Relative error after decryption (dim {}): {}",
            dimension_exponent, relative_error
        );
    }
}

// criterion_group!(benches, bench_encode, bench_decode, bench_encrypt_decrypt, bench_encrypt_add_decrypt, bench_encrypt_mul_decrypt, bench_apply_polynomial);

// Create a Criterion configuration with a reduced sample size
fn configure_criterion() -> Criterion {
    Criterion::default().sample_size(10) // Adjust the sample size as needed
}

criterion_group! {
    name = benches;
    config = configure_criterion();
    targets = bench_encode, bench_decode, bench_encrypt_decrypt, bench_encrypt_add_decrypt, bench_encrypt_mul_decrypt, bench_apply_polynomial
}

criterion_main!(benches);
