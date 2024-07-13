use crate::algebra::complex::{C64, Complex};

pub fn identity_matrix(n: usize) -> Vec<Vec<C64>> {
    let mut matrix = vec![vec![C64::new(0.0, 0.0); n]; n];
    for i in 0..n {
        matrix[i][i] = C64::new(1.0, 0.0);
    }
    matrix
}


pub fn multiply_matrices(a: &[Vec<C64>], b: &[Vec<C64>]) -> Vec<Vec<C64>> {
    let n = a.len();
    let mut result = vec![vec![C64::new(0.0, 0.0); n]; n];
    
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                result[i][j] = result[i][j] + a[i][k] * b[k][j];
            }
        }
    }

    result
}

