# CKKS Homomorphic Encryption

A rust implementation of the homomorphic encryption scheme from Cheon–Kim–Kim–Song [Homomorphic encryption for arithmetic of approximate numbers](https://eprint.iacr.org/2016/421.pdf).

## Quick start
Build with `cargo build`. Test with `cargo test`.
You can code directly in [/src/main.rs](https://github.com/M-Bln/ckks_he/blob/master/src/main.rs) and modify the provided example.
For quick tests, generate a pair of key with
```
let (mut client_key, server_key) = generate_pair_keys_toy();
```
Create a real plaintext (default dimension is 8), and specify an upper bound, large than any coefficient of the plaintext
```
let real_plaintext1 = vec![0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
let plaintext1_bound = 10.0;
let plaintext1 = to_plaintext(&real_plaintext1);
```
Encrypt with
```
let ciphertext1 = client_key.encrypt(&plaintext1, plaintext1_bound).unwrap();
```
Perform homomorphic operations such as multiplication
```
let result_mul = server_key.mul(&ciphertext1, &ciphertext1).unwrap();
```
Decrypt the result
```
let decrypted_mul = client_key.decrypt(&result_mul);
```
Run with `cargo run`.

## Advanced features
### Customed parameters
The dimension of the space of ciphertext is $2^{dimension_exponent}$, and the dimension of the space of plaintext $2^{dimension_exponent -1}$.
The initial encryption of a plaintext is a polynomial with coefficients modulo $q^level_max$. The higher the dimension the higher the security
but the slowest the computation. Specify the dimension of the space of ciphertext and the maximum level of ciphertext:
```
let level_max = 5;
let dimension_exponent = 13; // For these parameters the security reaches 86 bits
let (mul client_key, server_key) = generate_pair_keys_default(dimension_exponent, level_max);
```
