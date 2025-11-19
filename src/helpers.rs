

pub fn poly_size_n() -> impl Fn(usize, usize) -> Vec<usize> {
    |n, degree| {
        let mut buf = Vec::with_capacity(degree);
        for i in 0..degree {
            buf.push(n.pow(i as u32));
        }
        buf
    }
}