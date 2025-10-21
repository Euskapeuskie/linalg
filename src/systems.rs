use crate::types::{Matrix, Vector};



// Linear system of equations Ax = b
pub struct LinearSystem<T, const M: usize, const N: usize> {
    A: Matrix<T, M, N>,
    x: Option<Vector<T, N>>,
    b: Vector<T, N>
}

impl<T, const M: usize, const N: usize> LinearSystem<T, M, N>
where 
    T: Copy + num_traits::Num,
{
    pub fn new(A: Matrix<T, M, N>, b: Vector<T, N>) -> Self {
        Self {
            A: A,
            x: None,
            b: b
        }
    }

    pub fn least_sq(&mut self) -> Vector<T, N> {
        todo!()
    }
}



impl<T, const N: usize> LinearSystem<T, N, N>
where 
    T: Copy + num_traits::Num,
{
    /// Solve the linear system with PLU decomposition (Ax = b -> LUx = Pb)
    /// 1) Ly = Pb
    /// 2) Ux = y
    pub fn solve(&mut self) -> Result<(), String>
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
        self.x = Some(x);
        Ok(())
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
        // let a = Matrix::<f64, 3, 3>::ones();
        let a = Matrix::from(a);
        let b = Vector::<f64, 3>::ones();
        let mut sys = LinearSystem::new(a,b);
        sys.solve();
        if let Some(x) = sys.x {
            println!("Solution to system:\n{x}");
        }
    }
}