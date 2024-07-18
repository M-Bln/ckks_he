mod algebra;
mod ciphertext;
mod encoding;
mod keys;
mod random_distributions;

use crate::algebra::polynomial::Polynomial;
use crate::keys::client_key::{calculate_relative_error, to_plaintext, Plaintext};
use crate::keys::key_generator::generate_pair_keys_toy;
use bnum::types::I1024;

fn main() {
    // Generate keys with toy parameters
    let (mut client_key, server_key) = generate_pair_keys_toy();

    // Define some toy parameters
    let plaintext_bound = 10.0;

    // Chose a frist plaintext
    // Adjust length of the vector to client_key.dimension_plaintext() (default is 8)
    let real_plaintext1 = vec![0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let plaintext1 = to_plaintext(&real_plaintext1);

    // Chose a second plaintext
    let real_plaintext2 = vec![7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.1];
    let plaintext2 = to_plaintext(&real_plaintext2);

    // Encrypt the plaintexts
    let ciphertext1 = client_key.encrypt(&plaintext1, plaintext_bound).unwrap();
    let ciphertext2 = client_key.encrypt(&plaintext2, plaintext_bound).unwrap();

    // Perform homomorphic multiplication
    let result_mul = server_key.mul(&ciphertext1, &ciphertext2).unwrap();

    // Decrypt the result
    let decrypted_mul = client_key.decrypt(&result_mul);

    // Calculate error
    let expected_mul: Plaintext = plaintext1
        .iter()
        .zip(plaintext2.iter())
        .map(|(c1, c2)| *c1 * *c2)
        .collect();
    let error_mul = calculate_relative_error(&expected_mul, &decrypted_mul);
    println!(" Multiplication:");
    for (expected, obtained) in expected_mul.iter().zip(decrypted_mul.iter()) {
        println!("expected:{}, obtained:{}", expected, obtained);
    }
    println!("Multiplication relative error: {}", error_mul);

    // Polynomial for evaluation (2 + X + X^2 + X^3)
    let polynomial = Polynomial::<I1024>::new(vec![
        I1024::from(2),
        I1024::from(1),
        I1024::from(2),
        I1024::from(1),
    ]);

    let complex_polynomial = polynomial.to_c64();

    // Apply the polynomial to the ciphertext
    let result_poly = server_key
        .apply_polynomial(&polynomial, &ciphertext1)
        .unwrap();

    // Decrypt the result
    let decrypted_poly = client_key.decrypt(&result_poly);

    // Calculate error
    let expected_poly: Plaintext = plaintext1
        .iter()
        .map(|c| complex_polynomial.eval(*c))
        .collect();
    let error_poly = calculate_relative_error(&expected_poly, &decrypted_poly);
    println!("\n Polynomial evalutation:");
    for (expected, obtained) in expected_poly.iter().zip(decrypted_poly.iter()) {
        println!("expected:{}, obtained:{}", expected, obtained);
    }
    println!("Polynomial evaluation relative error: {}", error_poly);
}
