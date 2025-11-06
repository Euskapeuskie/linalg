use crate::types::{Matrix, Vector};



// Linear system of equations Ax = b
pub struct LinearSystem<T, const M: usize, const N: usize> {
    A: Matrix<T, M, N>,
    b: Vector<T, M>,
    x: Option<Vector<T, N>>,
    total_error: Option<T>,
}

impl<T, const M: usize, const N: usize> LinearSystem<T, M, N>
where 
    T: Copy + num_traits::Num,
{
    pub fn new(A: Matrix<T, M, N>, b: Vector<T, M>) -> Self {
        Self {
            A: A,
            b: b,
            x: None,
            total_error: None,
        }
    }


    /// Linear least squares approximation to fit the linear system to the given bs
    /// Returns the resulting best approximation x_hat
    /// 
    /// Also updates internally:
    /// self.x = x_hat
    /// self.total_error = (b - A*x_hat).magnitude().powi(2)
    /// 
    /// Uses QR-decomposition as a numerically stable solution
    /// A = QR
    /// A * x_hat = b
    /// -> QR * x_hat = b | *Q^T (with Q^T = Q^-1 for orthonormal basis)
    /// -> R * x_hat = Q^T * b
    /// R is upper triangular
    pub fn least_sq(&mut self) -> Result<Vector<T, N>, String> 
    where
        T: num_traits::Float,
    {
        let mut x_hat = Vector::<T, N>::zeros();
        let (q, r) = self.A.qr_decomposition();

        let mut res = &q.transpose() * &self.b;
        for i in (0..N).rev() {
            for j in i..N {
                res[i][0] = res[i][0] - r[i][j]*x_hat[j][0]
            }
            x_hat[i][0] = res[i][0] / r[i][i];
        }

        self.x = Some(x_hat.clone());

        // total error = ||(measurement_vec - estimate_vec)||²
        self.total_error = Some((&self.b - &(&self.A * &x_hat)).magnitude().powi(2));

        Ok(x_hat)
    }
}



impl<T, const N: usize> LinearSystem<T, N, N>
where 
    T: Copy + num_traits::Num,
{
    /// Solves the linear system
    /// Returns a Result that contains either the Vector that solves the system or an Error if the system cannot be solved.
    /// Assumes matrix A is invertible
    /// Updates the variable x of the LinearSystem
    /// 
    /// 
    /// This is done using PLU decomposition (Ax = b -> LUx = Pb)
    /// 1) Ly = Pb
    /// 2) Ux = y
    pub fn solve(&mut self) -> Result<Vector<T, N>, String>
    where
        T: num_traits::Float
    {
        // plu-factorization
        let (p, l, u) = self.A.lu_decomposition()?;

        // Calculate y from Ly = Pb
        let pb = &p * &self.b;
        let mut y = Vector::<T, N>::zeros();
        for i in 0..N {
            // 1. Eintrag kann man direkt ablesen da lower triangular mit 1 auf Hauptdiagonalen
            if i == 0 {
                y[i][0] = pb[i][0]
            }
            // Ansonsten: y[i] = pb[i] - summe
            else {
                let mut ans = T::zero();
                for j in 0..i {
                    ans = ans + l[i][j] * y[j][0];
                }
                y[i][0] = pb[i][0] - ans;
            }
        }

        // Calculate x from Ux = y
        let mut x = Vector::<T, N>::zeros();
        
        for i in (0..N).rev() {
            // Letzten Eintrag kann man direkt ablesen da upper triangular (ABER keine 1 auf Hauptdiagonalen!!!)
            if i == N-1 {
                x[i][0] = y[i][0] / u[i][i];
            }
            else {
                let mut ans = T::zero();
                for j in i+1..N {
                    ans = ans + u[i][j] * x[j][0];
                }
                x[i][0] = (y[i][0] - ans) / u[i][i];
            }
        }
        self.x = Some(x.clone());
        self.total_error = Some(T::zero());
        Ok(x)
    }
}



mod test {
    use crate::{systems::LinearSystem, types::Matrix, types::Vector};

    #[test]
    fn solve_full_rank() {
        let _a = [
            [1., 2., 0., 13.],
            [0., -1., 0., 0.],
            [0., 0., -1., 4.],
            [0., 0., 4., 1.],
        ];
        let a = [
            [4., 3., 2.],
            [3., 2., 1.],
            [2., 1., 3.],
        ];
        let a = [
            [1., 0., 0.],
            [0., 1., 0.],
            [1., 1., 0.],
        ];
        let a = [
            [2., -1., 0., -1.],
            [-1., 3., -2., 0.],
            [0., -2., 4., -2.],
            [-1., 0., -2., 3.],
        ];
        let b = [[1.], [0.], [-1.], [0.]];
        // let a = Matrix::<f64, 3, 3>::ones();
        let a = Matrix::from(a);
        let b = Vector::from(b);
        let mut sys = LinearSystem::new(a,b);
        sys.solve();
        if let Some(x) = sys.x {
            println!("Solution to system:\n{x}");
        }
    }

    #[test]
    fn least_sq() {
        let a = [
            [1., 1.],
            [1., 2.],
            [1., 3.],
        ];
        let a = Matrix::from(a);
        let b = [[1., 2., 3.]];
        let b = Matrix::from(b).transpose();
        let mut sys = LinearSystem::new(a, b);
        let x = sys.least_sq().unwrap();
        println!("{x}");
        let error = sys.total_error.unwrap();
        println!("{error}");
    }
}