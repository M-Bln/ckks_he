# CKKS Homomorphic Encryption

A rust implementation of the homomorphic encryption scheme from Cheon–Kim–Kim–Song [Homomorphic encryption for arithmetic of approximate numbers](https://eprint.iacr.org/2016/421.pdf).

## Quick start
Build with `cargo build`. Test with `cargo test`.
You can code directly in [/src/main.rs](https://github.com/M-Bln/ckks_he/blob/master/src/main.rs) and modify the provided example.
For quick tests, generate a pair of key with
```rust
let (mut client_key, server_key) = generate_pair_keys_toy();
```
Create a real plaintext (default dimension is 8), and specify an upper bound, large than any coefficient of the plaintext
```rust
let real_plaintext1 = vec![0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
let plaintext1_bound = 10.0;
let plaintext1 = to_plaintext(&real_plaintext1);
```
Encrypt with
```rust
let ciphertext1 = client_key.encrypt(&plaintext1, plaintext1_bound).unwrap();
```
Perform homomorphic operations such as multiplication
```rust
let result_mul = server_key.mul(&ciphertext1, &ciphertext1).unwrap();
```
Decrypt the result
```rust
let decrypted_mul = client_key.decrypt(&result_mul);
```
Run with `cargo run`.

## Advanced features
### Customed parameters
For a dimension exponent $h$ the dimension of the space of ciphertext is $2^h$ and the dimension of the space of plaintext is $2^{h-1}$.
For a level max $L$, the initial encryption of a plaintext is a polynomial with coefficients in  $\mathbb{Z}/q^level_max\mathbb{Z}$.
The higher the dimension the higher the security but the slowest the computation. Specify the dimension of the space of ciphertext and the maximum level of ciphertext with
```
let level_max = 5;
let dimension_exponent = 13; // For these parameters the security reaches 86 bits
let (mul client_key, server_key) = generate_pair_keys_default(dimension_exponent, level_max);
```
### Polynomial evaluation
The CKKS scheme allows for the homomorphic evaluation of polynomials on encrypted data.
```rust
// Polynomial for evaluation (2 + X + X^2 + X^3)
let polynomial = Polynomial::<I1024>::new(vec![
I1024::from(2),
I1024::from(1),
I1024::from(2),
I1024::from(1),
]);

// Apply the polynomial to the ciphertext
let result_poly = server_key.apply_polynomial(&polynomial, &ciphertext1).unwrap();
```
## Benchmark
Run the benchmark with
```
cargo bench
```
Benchmark various parameters by modifying the `const` in `/benches/full_bench.rs`
```rust
const PLAINTEXT_BOUND: f64 = 4.0;
const Q_EXPONENT: u32 = 30;
const LEVEL_MAX: u32 = 5;
const EXPONENT: u32 = 2;
```
To bench various exponent at once modify the vector
```rust
let dimension_exponents = vec![EXPONENT];	
```
in each bench. Notice than with `dimension_exponent = 13` the bench taked roughly half an hour.

## Build and browse documentation
Build the documentation with `cargo doc --no-deps`. Browse it by opening the resulting `/target/doc/ckks_he/index.html` in a web browser.