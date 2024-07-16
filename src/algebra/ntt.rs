use crate::algebra::big_int::BigInt;
use crate::algebra::arithmetic::RingMod;

fn ntt<T: BigInt>(a: &mut [RingMod<T>], root: RingMod<T>) {
    let n = a.len();
    assert!(n.is_power_of_two());

    for i in (1..n.trailing_zeros()).rev() {
        let m = 1 << i;
        let root_m = root.pow((n / m) as u64);
        for j in (0..n).step_by(m << 1) {
            let mut w = RingMod::new(T::from(1), root.modulus);
            for k in 0..m {
                let t = w.mul(a[j + k + m]);
                a[j + k + m] = a[j + k].sub(t);
                a[j + k] = a[j + k].add(t);
                w = w.mul(root_m);
            }
        }
    }

    let mut j = 0;
    for i in 1..n - 1 {
        let mut k = n >> 1;
        while j & k != 0 {
            j ^= k;
            k >>= 1;
        }
        j ^= k;
        if i < j {
            a.swap(i, j);
        }
    }
}

fn intt<T: BigInt>(a: &mut [RingMod<T>], root: RingMod<T>) {
    let n = a.len();
    assert!(n.is_power_of_two());
    let inv_n = RingMod::new(T::from(n as u64), root.modulus).inv();

    for i in (1..n.trailing_zeros()).rev() {
        let m = 1 << i;
        let root_m = root.pow((n / m) as u64).inv();
        for j in (0..n).step_by(m << 1) {
            let mut w = RingMod::new(T::from(1), root.modulus);
            for k in 0..m {
                let t = w.mul(a[j + k + m]);
                a[j + k + m] = a[j + k].sub(t);
                a[j + k] = a[j + k].add(t);
                w = w.mul(root_m);
            }
        }
    }

    let mut j = 0;
    for i in 1..n - 1 {
        let mut k = n >> 1;
        while j & k != 0 {
            j ^= k;
            k >>= 1;
        }
        j ^= k;
        if i < j {
            a.swap(i, j);
        }
    }

    for x in a.iter_mut() {
        *x = x.mul(inv_n);
    }
}
