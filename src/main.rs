extern crate blas;

use blas::c::*;

fn main() {
    let n = 5;
    let k = 10;
    let alpha = 1.0;
    let beta = 0.0;

    let Sigma = vec![1.0; (n * n) as usize];
    let H = vec![1.0; (k * n) as usize];
    let mu = vec![1.0; n as usize];
    let data = vec![1.0; k as usize];
    let var_11 = vec![1.0; n as usize];
    let var_5 = vec![1.0; (n * n) as usize];

    let mut var_20 = vec![1.0; (n * n) as usize];
    let mut INFO = -1;
    let mut muvar_2 = vec![1.0; (n * n) as usize];
    let mut Sigmavar_2 = vec![1.0; (n * n) as usize];
    let mut var_17 = vec![1.0; (n * k) as usize];
    //TODO is this R right?
    let mut R = vec![1.0; (k * k) as usize];
    let mut var_19 = vec![1.0; (n * n) as usize];
    let mut Hvar_2 = vec![1.0; (k * n) as usize];

    dcopy(n * n, &var_5, 1, &mut var_20, 1);

    dsymm(Layout::ColumnMajor, Side::Left, Part::Upper,
          n, k, alpha,
          &Sigma, n, &H, k, beta, &mut var_17, n);

    dgemm(Layout::ColumnMajor, Transpose::None, Transpose::None,
          k, k, n, alpha,
          &H, k, &var_17, n, beta, &mut R, k);

    dcopy(n * n, &var_5, 1, &mut var_19, 1);

    dcopy(n, &mu, 1, &mut muvar_2, 1);

    dcopy(n * n, &Sigma, 1, &mut Sigmavar_2, 1);

    dcopy(k * n, &H, 1, &mut Hvar_2, 1);

    //TODO:
    /*
    call dgemm('N', 'N', k, 1, n, 1, H, k, mu, n, -1, data, k)
    call dposv('U', k, n, R, k, Hvar_2, k, INFO)
    call dposv('U', k, 1, R, k, data, k, INFO)
    call dgemm('N', 'N', n, n, k, 1, H, k, Hvar_2, k, 0, var_19, n)
    call dgemm('N', 'N', n, 1, k, 1, H, k, data, k, 0, var_11, n)
    call dsymm('L', 'U', n, n, 1, var_19, n, Sigma, n, 0, var_20, n)
    call dsymm('L', 'U', n, 1, 1, Sigma, n, var_11, n, 1, muvar_2, n)
    call dsymm('L', 'U', n, n, -1, Sigmavar_2, n, var_20, n, 1, Sigmavar_2, n)
*/
}