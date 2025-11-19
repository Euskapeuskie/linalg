use crate::types::Matrix;


impl<T, const M: usize, const N: usize> Matrix<T, M, N> {

    /// Compress the given matrix using FFT, keeps only the most dominant frequencies
    pub fn fft_compress(&self, compression_rate: f64) -> () {
        todo!()
    }

    /// Compress the given vector using haar wavelets
    pub fn hwt_compress(&self, compression_rate: f64) -> () {

    }

    /// Compress the given matrix using SVD
    pub fn svd_compress(&self, compression_rate: f64) -> () {
        todo!()
    }
}