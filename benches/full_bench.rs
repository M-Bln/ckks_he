use bnum::types::I1024;
use ckks::algebra::big_int::BigInt;
use ckks::algebra::complex::{Complex, C64};
use ckks::encoding::Encoder;
use ckks::keys::client_key::to_plaintext;
use ckks::keys::key_generator::{generate_pair_keys_default};
use ckks::random_distributions::generate_random_vector;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::Rng;

const PLAINTEXT_BOUND: f64 = 4.0;
const Q_EXPONENT: u32 = 30;
const LEVEL_MAX: u32 = 5;

fn bench_encode(c: &mut Criterion) {
    let dimension_exponents = vec![8];
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
    let dimension_exponents = vec![8];
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
    let dimension_exponents = vec![8]; // Try smaller exponents first
    for &dimension_exponent in &dimension_exponents {
        // Generate keys
        let (mut client_key, _server_key) = generate_pair_keys_default::<I1024>(dimension_exponent, LEVEL_MAX);

        // Generate plaintext
        let real_plaintext = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
	let plaintext=to_plaintext(&real_plaintext);

        // Encrypt the plaintext
        c.bench_function(&format!("encrypt_dim_{}", dimension_exponent), |b| {
            b.iter(|| {
                client_key.encrypt(black_box(&plaintext), PLAINTEXT_BOUND).unwrap()
            })
        });

        let ciphertext = client_key.encrypt(&plaintext, PLAINTEXT_BOUND).unwrap();

        // Decrypt the ciphertext
        c.bench_function(&format!("decrypt_dim_{}", dimension_exponent), |b| {
            b.iter(|| {
                let decrypted = client_key.decrypt(black_box(&ciphertext));
                // calculate_relative_error(&plaintext, &decrypted);
            })
        });

        let decrypted = client_key.decrypt(&ciphertext);
        let relative_error = calculate_relative_error(&plaintext, &decrypted);
        println!("Relative error after decryption (dim {}): {}", dimension_exponent, relative_error);
    }
}

fn bench_encrypt_add_decrypt(c: &mut Criterion) {
    let dimension_exponents = vec![8];
    for &dimension_exponent in &dimension_exponents {
        // Generate keys
        let (mut client_key, server_key) = generate_pair_keys_default::<I1024>(dimension_exponent, LEVEL_MAX);

        // Generate ciphertexts
        let real_plaintext1 = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
	let plaintext1=to_plaintext(&real_plaintext1);
        let ciphertext1 = client_key.encrypt(&plaintext1, PLAINTEXT_BOUND).unwrap();
	
        let real_plaintext2 = generate_random_vector(
            1 << (dimension_exponent - 1),
	    -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
	let plaintext2=to_plaintext(&real_plaintext2);
	let ciphertext2 = client_key.encrypt(&plaintext2, PLAINTEXT_BOUND).unwrap();
	

        c.bench_function(&format!("homomorphic add dimension: 2^{}", dimension_exponent), |b| {
            b.iter(|| {
                server_key.add(black_box(&ciphertext1), black_box(&ciphertext2));
            })
        });
        


        // // Decrypt the ciphertext
        // c.bench_function(&format!("decrypt_dim_{}", dimension_exponent), |b| {
        //     b.iter(|| {
        //         let decrypted = client_key.decrypt(black_box(&ciphertext));
        //         // calculate_relative_error(&plaintext, &decrypted);
        //     })
        // });
	let expected: Vec::<C64> = plaintext1.iter().zip(plaintext2.iter()).map(|(c1,c2)| *c1+*c2).collect();
	let result = &server_key.add(&ciphertext1,&ciphertext2);
        let decrypted = client_key.decrypt(&result);
	let expected_error = client_key.rescaled_error(&result);
        let relative_error = calculate_relative_error(&expected, &decrypted);
	println!("Expected error: {}", expected_error);
        println!("Relative error after decryption (dim {}): {}", dimension_exponent, relative_error);
    }
}


fn bench_encrypt_mul_decrypt(c: &mut Criterion) {
    let dimension_exponents = vec![8];
    for &dimension_exponent in &dimension_exponents {
        // Generate keys
        let (mut client_key, server_key) = generate_pair_keys_default::<I1024>(dimension_exponent, LEVEL_MAX);

        // Generate ciphertexts
        let real_plaintext1 = generate_random_vector(
            1 << (dimension_exponent - 1),
            -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
	let plaintext1=to_plaintext(&real_plaintext1);
        let ciphertext1 = client_key.encrypt(&plaintext1, 100.0*PLAINTEXT_BOUND).unwrap();
	
        let real_plaintext2 = generate_random_vector(
            1 << (dimension_exponent - 1),
	    -PLAINTEXT_BOUND,
            PLAINTEXT_BOUND,
        );
	let plaintext2=to_plaintext(&real_plaintext2);
	let ciphertext2 = client_key.encrypt(&plaintext2, 100.0*PLAINTEXT_BOUND).unwrap();
	

        c.bench_function(&format!("homomorphic multiplication, dimension: 2^{}", dimension_exponent), |b| {
            b.iter(|| {
                server_key.mul(black_box(&ciphertext1), black_box(&ciphertext2)).unwrap();
            })
        });
        


        // // Decrypt the ciphertext
        // c.bench_function(&format!("decrypt_dim_{}", dimension_exponent), |b| {
        //     b.iter(|| {
        //         let decrypted = client_key.decrypt(black_box(&ciphertext));
        //         // calculate_relative_error(&plaintext, &decrypted);
        //     })
        // });
	let expected: Vec::<C64> = plaintext1.iter().zip(plaintext2.iter()).map(|(c1,c2)| *c1*c2).collect();
	let result = &server_key.mul(&ciphertext1,&ciphertext2).unwrap();
        let decrypted = client_key.decrypt(&result);
	let expected_error = client_key.rescaled_error(&result);
        let relative_error = calculate_relative_error(&expected, &decrypted);
	println!("Expected error: {}", expected_error);
        println!("Relative error after decryption (dim {}): {}", dimension_exponent, relative_error);
    }
}


fn calculate_relative_error(original: &[C64], decrypted: &[C64]) -> f64 {
    original.iter()
        .zip(decrypted.iter())
        .map(|(o, d)| {
            let error = (*o-*d).magnitude();
	    let relative_error = error / (o.magnitude());
	    println!("expected: {}", o);
	    println!("decrypted: {}", d);
	    println!("error: {}", error);
	    println!("relative error: {}", relative_error);
	    relative_error
        })
        .fold(0.0, |max_error, current_error| max_error.max(current_error))
}

criterion_group!(benches, bench_encode, bench_decode, bench_encrypt_decrypt, bench_encrypt_add_decrypt, bench_encrypt_mul_decrypt);
criterion_main!(benches);
