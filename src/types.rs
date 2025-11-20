use crate::complex_impl::ToComplex;

use num_traits;
use num_complex::Complex;



// general 2d-matrix
#[derive(Debug, Clone)]
pub struct Matrix<T, const M: usize, const N: usize> {
    // M rows, N columns
    pub data: [[T; N]; M],
}


// column vector
pub type Vector<T, const M: usize> = Matrix<T, M, 1>;


// General rectangular matrices
impl<T, const M: usize, const N: usize> Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
{
    /// Returns the number of rows (M) of the matrix
    /// 
    /// # Example
    /// ```rust
    /// use linalg::types::Matrix;
    /// 
    /// let c: Matrix<usize, 10, 3> = Matrix::ones();
    /// let n_rows = c.n_rows();
    /// assert_eq!(n_rows, 10);
    /// ```
    pub fn n_rows(&self) -> usize {
        self.data.len()
    }

    /// Returns the number of columns (N) of the matrix
    /// 
    /// # Example
    /// ```rust
    /// use linalg::types::Matrix;
    /// 
    /// let c: Matrix<usize, 10, 3> = Matrix::ones();
    /// let n_rows = c.n_cols();
    /// assert_eq!(n_rows, 3);
    /// ```
    pub fn n_cols(&self) -> usize {
        let maybe_n_cols = self.data.first();
        if let Some(n_cols) = maybe_n_cols {
            return n_cols.len();
        }
        0
    }

    /// create a matrix of size MxN with just 0s
    /// 
    /// # Example
    /// ```rust
    /// use linalg::types::Matrix;
    /// 
    /// let c: Matrix<usize, 10, 3> = Matrix::zeros();
    /// assert_eq!(c[2], [0, 0, 0]);
    /// ```
    pub fn zeros() -> Self {
        let data = [[T::zero(); N]; M];
        Matrix { data: data }
    }

    /// create a matrix of size MxN with just 1s
    /// 
    /// # Example
    /// ```rust
    /// use linalg::types::Matrix;
    /// 
    /// let c: Matrix<usize, 10, 3> = Matrix::ones();
    /// assert_eq!(c[2], [1, 1, 1]);
    /// ```
    pub fn ones() -> Self {
        let data = [[T::one(); N]; M];
        Matrix { data: data }
    }


    /// Returns the rank of the matrix
    /// 
    /// # Example
    /// ```rust
    /// use linalg::types::Matrix;
    /// 
    /// let a: Matrix<f64, 5, 5> = Matrix::identity();
    /// assert_eq!(a.rank(), 5);
    /// ```
    pub fn rank(&self) -> usize
    where 
        T: num_traits::Float,
    {
        let rref = self.row_echelon();
        let mut rank = 0;

        for i in 0..M {
            if rref[i][i].abs() > T::epsilon() {
                rank += 1;
            }
            else {
                break;
            }
        }
        rank
    }

    /// create a matrix from a function that generates each row
    /// The resulting matrix has to be of size f().len() x n
    /// 
    /// # Example
    /// ```rust
    /// use linalg::types::Matrix;
    /// 
    /// let c: Matrix<usize, 10, 3> = Matrix::from_function(|x: usize| [1, x, x.pow(2)]);
    /// assert_eq!(c[2], [1, 2, 4]);
    /// ```
    pub fn from_function<F>(f: F) -> Self
    where 
        F: Fn(usize) -> [T; N],
    {
        let mut data = [[T::zero(); N]; M];
        for i in 0..M {
            data[i] = f(i);
        }
        Matrix { data: data }
    }


    /// create a matrix from two functions - one that generates the rows, one that generates the columns
    /// 
    /// 
    pub fn from_function_xy<F, G>(f: F, g: G) -> Self
    where 
        F: Fn(usize) -> [T; N],
        G: Fn(usize) -> [T; M],
    {
        let data = [[T::zero(); N]; M];
        let mut ans = Matrix { data: data };
        for i in 0..M {
            ans[i] = f(i);
        }

        let mut ans_t = ans.transpose();
        for j in 0..N {
            let col = g(j);
            for k in 0..M {
                ans_t[j][k] = ans_t[j][k] * col[k];
            }
        }
        let ans = ans_t.transpose();
        ans
    }


    /// As slice
    pub fn as_slice(&self) -> &[T] {
        &self[..M]
    }


    /// transposing is just swapping rows and columns, so
    /// AT_ji = A_ij
    pub fn transpose(&self) -> Matrix<T, N, M> {
        let mut ans = Matrix::<T, N, M>::zeros();
        // row
        for i in 0..M {
            // col
            for j in 0..N {
                ans[j][i] = self[i][j];
            }
        }
        ans
    }


    /// Tranpose of a complex matrix (Hermetian)
    pub fn hermetian<U>(&self) -> Matrix<Complex<U>, N, M>
    where 
        T: ToComplex<U>,
        U: Copy + num_traits::Num + num_traits::Float,
    {
        let mut ans = [[U::zero().to_complex(); M]; N];
        for i in 0..M {
            for j in 0..N {
                let x = self[i][j].to_complex();
                ans[i][j] = x.conj();
            }
        }
        let ans = Matrix::from(ans);
        ans

    }


    /// Bring a matrix to row echelon form
    pub fn row_echelon(&self) -> Self 
    where
        T: num_traits::Float

    {
        let mut u = self.clone();
        let max_iter = M.min(N);

        // 1) Permutation - put the biggest element by absolute value in the upper most position as pivot (numeric stability)
        // 2) Elimination - row by row
        for i in 0..max_iter {
            // find the row index where that column has the maximum value
            let maybe_row_max_index = u.data.iter().enumerate().skip(i) // only look at rows from i downwards
                .map(|(row_i, row)| (row_i, row[i])) // look at the i_th element in each row
                .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).expect("Unable to order rows"))
                .map(|(index, _)| index);

            // swap rows if necessary
            if let Some(row_max_index) = maybe_row_max_index {
                // if the largest element of column i isn't in row i, then swap the rows
                u.data.swap(row_max_index, i);
            }

            // eliminate rows
            let (top, bottom) = u.data.split_at_mut(i+1);
            let pivot_row = &top[i];
            let pivot_factor = pivot_row[i];

            // Check if pivot is a 0 -> if yes go to the next row
            if pivot_factor.is_zero() {
                continue
            }

            for j in 0..bottom.len() {
                let row_j = &mut bottom[j];
                let factor = row_j[i] / pivot_factor;
                for k in i..N {
                    row_j[k] = row_j[k] - factor*pivot_row[k]
                }
            }
        }
        u
    }


    /// Compute the QR decomposition of the matrix that satisfies A = QR using Gram-Schmidt
    /// - Q: Orthonomormal basis of the column space of A. Q is of size m x n
    /// - R: Upper triangular matrix. R is of size n x n
    /// Returns (q, r)
    /// 
    /// Internal workings:
    /// 1) choose the first non-zero column as my initial basis
    /// 2) calculate for each following column...
    pub fn qr_decomposition(&self) -> (Self, Matrix<T, N, N>)
    where 
        T: num_traits::Float
    {

        // transpose columns into rows just more convenient to work with
        let mut q = Self::zeros().transpose();

        // Build the Q-Matrix row by row
        // q[i] is supposed to be the first column of the later q, that's why transposing after the for loop
        for i in 0..self.n_cols() {
            q[i] = {
                // b is the vector we want to orthonormalize to q
                let mut b = Vector::from(self.transpose()[i].map(|x| [x]));

                // subtract the components of b in the other orthonormal directions to get only the new component of b
                for j in 0..i {
                    let a = Vector::from(q[j].map(|x| [x]));
                    b  = &b - &(&a * (((&a.transpose() * &b)[0][0]) / ((&a.transpose() * &a)[0][0])));
                }
                // normalize b to unit length
                b = &b * (T::one()/b.magnitude());
                b.to_array()
            };
        }


        // Build the R-Matrix:
        // A = QR with Q^T = Q^-1 -> R = Q^T*A
        let r = &q * self;
        
        // transpose Q so orthonormal basis are now column vectors
        (q.transpose(), r)
    }


    pub fn svd(&self) -> () {
        todo!()
    }

    pub fn pinv(&self) -> Self {
        todo!()
    }
}



// Square matrices
impl<T, const N: usize> Matrix<T, N, N>
where
    T: Copy + num_traits::Num
{
    /// Create an Identity matrix of size NxN
    pub fn identity() -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = T::one();
        }
        Self { data: data }
    }


    /// Check if the matrix is symmetric
    pub fn is_sym(&self) -> bool {
        if self == self.transpose() {
            true
        }
        else {
            false
        }
    }

    /// Performs the LU factorization on the matrix A to satisfy the equation PA = LU and returns a tuple
    /// 
    /// (P, L, U)
    ///  
    /// - P: Permutation matrix
    /// - L: Lower triangular matrix
    /// - U: Upper traingular matrix
    /// This only works for floating point type matrices.
    pub fn lu_decomposition(&self) -> Result<(Self, Self, Self), String>
    where T: num_traits::Float,
    {
        let mut p: Matrix<T, N, N> = Self::identity();
        let mut l = Self::zeros();
        let mut u = self.clone();
        

        // 1) Permutation - put the biggest element by absolute value in the upper most position as pivot (numeric stability)
        // 2) Elimination - row by row
        for i in 0..N {
            // find the row index where that column has the maximum value
            let maybe_row_max_index = u.data.iter().enumerate().skip(i) // only look at rows from i downwards
                .map(|(row_i, row)| (row_i, row[i])) // look at the i_th element in each row
                .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).expect("Unable to order rows"))
                .map(|(index, _)| index);
            
            // swap rows if necessary
            if let Some(row_max_index) = maybe_row_max_index {
                // if the largest element of column i isn't in row i, then swap the rows
                u.data.swap(row_max_index, i);
                l.data.swap(row_max_index, i);
                // calculate new permutation matrix
                let mut p_i: Matrix<T, N, N> = Self::identity();
                p_i.data.swap(row_max_index, i);
                p = &p_i * &p;
            }

            // eliminate rows
            let (top, bottom) = u.data.split_at_mut(i+1);
            let pivot_row = &top[i];
            let pivot_factor = pivot_row[i];
            // Check if division by 0
            if pivot_factor.is_zero() {
                return Err(String::from("Singular matrix cannot be factorized"));
            }

            for j in 0..bottom.len() {
                let row_j = &mut bottom[j];
                let factor = row_j[i] / pivot_factor;
                // update L matrix: stores information
                // "how did I eliminate element [i][j]"
                l.data[i+1+j][i] = factor;
                // for every column in row j
                for k in i..N {
                    row_j[k] = row_j[k] - factor*pivot_row[k]
                }
            }

        }
        // add the ones for the main diagonal of L
        l = &l + &Self::identity();
        Ok((p, l, u))
    }


    /// Returns the inverse of the matrix or an error if the matrix is singular
    /// uses LU-factorization for calculating the inverse
    /// A * A.inv() = I | with PA = LU
    /// LU * A.inv() = P*I | with Y = U * A.inv()
    /// LY = P -> solve for Y
    /// Resubstitution: U*A.inv() = Y -> solve for A.inv()
    pub fn inv(&self) -> Result<Self, String>
    where
        T: num_traits::Float,
    {
        let (p, l, u) = self.lu_decomposition()?;

        // 1) Solve for LY = P
        // L is a lower triangular matrix with only 1s on the main diagonal!
        let mut y = Self::zeros();
        // i: row we're acting on of Y
        for i in 0..N {
            // j: column of Y
            for j in 0..N {
                let mut sum = T::zero();
                for k in 0..i {
                    sum = sum + l[i][k] * y[k][j];
                }
                y[i][j] = (p[i][j] - sum) / l[i][i] // l[i][i] is always 1
            }
        }

        // Solve for U*A.inv() = Y
        // U is a upper triangular matrix, so start with last row
        let mut a_inv = Self::zeros();
        // start with last row because U = [0, 0, ...., x]
        for i in (0..N).rev() {
            for j in 0..N {
                let mut ans = T::zero();
                for k in i+1..N {
                    ans = ans + u[i][k] * a_inv[k][j];
                }
                a_inv[i][j] = (y[i][j] - ans) / u[i][i];
            }
        }
        Ok(a_inv)
    }


    /// Calculates the determinant of a square matrix
    /// 1) Bring matrix in upper triangular form (keep track of row exchanges) !!!! CANNOT USE ROW_ECHOLONG METHOD BECAUSE THIS METHOD SETS PIVOTS TO 1
    /// 2) Determinant = (-1)^n_row_changes * product of pivots
    pub fn det(&self) -> T 
    where
        T: num_traits::Float
    {
        todo!();
    }

    pub fn eig(&self) -> () {
        todo!()
    }

    pub fn eigvals(&self) -> () {
        todo!()
    }

    /// Generate a N x N fourier matrix
    pub fn fourier_matrix() -> Matrix<Complex<T>, N, N>
    where
        T: num_traits::Float
    {
        let mut f = [[T::zero().to_complex(); N]; N];
        for i in 0..N {
            for j in 0..N {
                // use f64 for accuracy, convert to T (which might be f32) later
                let theta = -2.0 * std::f64::consts::PI * (i as f64) * (j as f64) / (N as f64);
                f[i][j] = Complex::new(T::zero(), T::from(theta).unwrap()).exp();
            }
        }
        let f = Matrix::from(f);
        f
    }
}



// Vectors
impl<T, const M: usize> Vector<T, M>
where
    T: Copy + num_traits::Num
{

    /// Returns the magnitude or length of a vector:
    /// sqrt(x_transpose * x)
    pub fn magnitude(&self) -> T
    where
        T: num_traits::Float,
    {
        let x = &(self.transpose()) * self;
        x[0][0].sqrt()
    }


    /// Vector into Array
    pub fn to_array(self) -> [T; M] {
        self.transpose()[0]
    }


    /// DFT transform -> projects time into frequency domain
    /// runs in O(n²)
    pub fn dft<U>(&self) -> Matrix<Complex<U>, M, 1>
    where
        T: ToComplex<U>,
        U: Copy + num_traits::Num + num_traits::Float,
    {

        let f = Matrix::<U, M, M>::fourier_matrix();

        // Convert self to a complex vector if it isn't already complex
        let x = Vector::from(self.data.map(|row| [row[0].to_complex()]));
        let x = f * x;
        x
    }


    /// Inverse DFT Transform -> projects frequency domain in time domain
    pub fn idft<U>(&self) -> Matrix<Complex<U>, M, 1>
    where
        T: ToComplex<U>,
        U: Copy + num_traits::Num + num_traits::Float,
    {
        let f = Matrix::<U, M, M>::fourier_matrix().hermetian();
        // Convert self to a complex vector if it isn't already complex
        let x = Vector::from(self.data.map(|row| [row[0].to_complex()]));
        let mut x = f * x;

        for i in 0..M {
            x[i][0] = x[i][0] * U::one()/U::from(M).unwrap();
        }
        x
    }


    /// FFT (not in place), runs in O(n*log(n)), recursive
    /// Uses FFT algorithm if input vector is a power of 2
    /// Else falls back to DFT algorithm (stack overflow guaranteed if input size big)
    pub fn fft(&self) -> Matrix<Complex<T>, M, 1>
    where
        T: num_traits::Float,
    {
        if !(self.n_rows().is_power_of_two()) {
            return self.dft();
        }
        let mut data = self.transpose()[0].map(|x| Complex::new(x, T::zero()));

        _fft(&mut data);

        /// Helper function for doing the fft algorithm recursively
        fn _fft<T>(x: &mut [Complex<T>]) -> ()
        where
            T: num_traits::Float
        {
            let n = x.len();

            // Base case
            if n < 2 {
                return;
            }

            // 1) Split x into x_even and x_odd components
            let mut x_even: Vec<_> = x.iter().step_by(2).cloned().collect();
            let mut x_odd: Vec<_> = x.iter().skip(1).step_by(2).cloned().collect();

            // 2) Recursive call for breaking down the xs further
            // once we hit length of x_even == 1 the function returns empty
            // -> we build the fourier matrix for size 2 first
            _fft(&mut x_even);
            _fft(&mut x_odd);

            // do the fourier transform on each half (x_even and x_odd) of the components
            // cannot use matrix unfortunately because non constant initialization
            for k in 0..n/2 {
                // twiddle factor in row k: d_k = exp(i*2*pi*k/n) --> Twiddle factor matrix D is diagonal
                let d_k = Complex::new(0.0, -2.0 * std::f64::consts::PI * (k as f64) / (n as f64)).exp();
                let d_k = Complex::new(T::from(d_k.re).unwrap(), T::from(d_k.im).unwrap());
                let t = d_k * x_odd[k];
                x[k] = x_even[k] + t;
                x[k + n/2] = x_even[k] - t;
            }
        }

        Matrix::from(data.map(|x| [x]))
    }
}




// TESTS
#[cfg(test)]
mod test {
    use crate::types::Matrix;
    use crate::types::Vector;

    #[test]
    fn create_matrix() {
        // zeros
        let a: Matrix<usize, 3, 2> = Matrix::zeros();
        assert_eq!(a[0], [0, 0]);

        // ones
        let b: Matrix<usize, 3, 2> = Matrix::ones();
        assert_eq!(b[0], [1, 1]);

        // identity
        let c: Matrix<usize, 3, 3> = Matrix::identity();
        assert_eq!(c[0], [1, 0, 0]);

        // from array
        let arr = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        let _m = Matrix::from(arr);

        // from function
        let c: Matrix<usize, 10, 3> = Matrix::from_function(|x| [1, x, x.pow(2)]);
        assert_eq!(c[2], [1, 2, 4]);

        // from function xy
        let d: Matrix<usize, 4, 4> = Matrix::from_function_xy(|x| [1, x, x.pow(2), x.pow(3)], |y| [1, y, y.pow(2), y.pow(3)]);
        println!("{d}");
    }

    #[test]
    fn print_matrix() {
        let a: Matrix<f64, 3, 2> = Matrix::ones();
        println!("{a}");
    }

    #[test]
    fn add_sub_matrices() {
        let a: Matrix<usize, 5, 2> = Matrix::ones();
        let b: Matrix<usize, 5, 2> = Matrix::zeros();
        let _c = &a + &b;
        let _d = &a - &b;
    }

    #[test]
    fn mul_matrices() {
        let a: Matrix<f64, 3, 3> = Matrix::ones();
        let b: Matrix<f64, 3, 3> = Matrix::ones();
        let d: Matrix<f64, 3, 3> = Matrix::ones();
        let _c = &(&a * &b) * &d;
        let _e = &a + &b;

        let a: Matrix<f64, 3, 7> = Matrix::ones();
        let b: Matrix<f64, 7, 3> = Matrix::ones();
        let _c = &a * &b;
    }

    #[test]
    fn mul_matrix_scalar() {
        // Matrix
        let a: Matrix<f64, 3, 2> = Matrix::ones();
        let s: f64 = 64.3;
        let _b = &a * s;

        // Vector
        let v = Vector::<f64, 5>::ones();
        let _w = &v * s;
    }

    #[test]
    fn transpose() {
        let a: Matrix<f64, 5, 2> = Matrix::ones();
        let _b = a.transpose();
    }

    #[test]
    fn lu_factors() {
        let arr = [[1., 0., 7.4], [0., 1., 0.], [1., 0., 0.]];
        let a = Matrix::from(arr);
        let (_p, _l, _u) = a.lu_decomposition().unwrap();
    }

    #[test]
    fn index() {
        let mut a = Matrix::<f64, 5, 5>::ones();
        let _x = a[3][4]; // index borrow
        let y = &mut a[3][4]; // index mut
        *y = 4.0;
        println!("{a}"); // index borrow again
    }

    #[test]
    fn inverse() {
        let a = [
            [1., 1., 1.],
            [2., 3., 3.],
            [3., 4., 5.],
        ];
        let a = Matrix::from(a);
        let inv = a.inv().unwrap();
        println!("Inverse:");
        println!("{inv}");
    }

    #[test]
    fn row_echelon() {
        let a = Matrix::<f64, 4, 4>::identity();
        let a = &a+&a;
        let t = a.row_echelon();
        println!("{t}");
    }

    #[test]
    fn qr_decomposition() {
        let a = [
            [1., 2., 4.],
            [0., 0., 5.],
            [0., 3., 6.],
            [0., 1., 2.],
            [1., 7., -1.]
        ];
        let a = Matrix::from(a);
        let (q, r) = a.qr_decomposition();
        println!("{q}");
        println!("{r}");
    }


    #[test]
    fn vector() {
        let a = [1, 2, 3, 4];
        let a = Vector::from(a.map(|x| [x]));
        let a = &a[..4];
        println!("{a:?}");
    }

    #[test]
    fn dft() {
        let a = (0..5).map(|x| [x as f64]).collect::<Vec<_>>();
        let arr: [[f64; 1]; 5] = a.try_into().expect("unable");
        println!("{arr:?}");
        let a = Vector::from(arr);

        let f = a.dft();
        println!("{f}");
        let g = f.idft::<f64>();
        println!("{g:?}");
    }

    #[test]
    fn fft() {
        let a = (0..23).map(|x| [x as f64]).collect::<Vec<_>>();
        let arr: [[f64; 1]; 23] = a.try_into().expect("unable");
        let a = Box::new(Vector::from(arr));
        let f = a.fft();
        println!("{f}");
    }

    #[test]
    fn symmetry() {
        let a = [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]
        ];
        let b = [
            [1, 2, 3],
            [2, 1, 4],
            [3, 4, 1],
        ];
        let a = Matrix::from(a);
        let b = Matrix::from(b);
        let t_1 = a.is_sym();
        let t_2 = b.is_sym();
        assert_eq!(t_1, false);
        assert_eq!(t_2, true);
    }
}