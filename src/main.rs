extern crate blas;

use blas::c::*;

fn main() {

    let (m, n, k) = (2, 4, 3);
    let a = vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0];
    let b = vec![1.0, 5.0, 9.0, 2.0, 6.0, 10.0, 3.0, 7.0, 11.0, 4.0, 8.0, 12.0];
    let mut c = vec![2.0, 7.0, 6.0, 2.0, 0.0, 7.0, 4.0, 2.0];

    dgemm(Layout::ColumnMajor, Transpose::None, Transpose::None,
          m, n, k, 1.0, &a, m, &b, k, 1.0, &mut c, m);

    assert_eq!(&c, &vec![40.0, 90.0, 50.0, 100.0, 50.0, 120.0, 60.0, 130.0]);

    let n = 10;
    let k = 1;
    

}


fn kalman_f90() {
    /*
    subroutine f(mu, Sigma, H, INFO, R, Sigmavar_2, data, muvar_2, k, n)
    implicit none

    integer, intent(in) :: k
    integer, intent(in) :: n
    real*8, intent(in) :: Sigma(n, n)        !  Sigma
    real*8, intent(in) :: H(k, n)            !  H
    real*8, intent(in) :: mu(n)              !  mu
    real*8, intent(in) :: R(k, k)            !  R, H*Sigma*H' + R
    real*8, intent(in) :: data(k)            !  (H*Sigma*H' + R)^-1*((-1)*data + H*mu), data, (-1)*   data + H*mu
    integer, intent(out) :: INFO             !  INFO
    real*8, intent(out) :: muvar_2(n)        !  mu, Sigma*H'*(H*Sigma*H' + R)^-1*((-1)*data + H*  mu) + mu
    real*8, intent(out) :: Sigmavar_2(n, n)  !  Sigma, (-1)*Sigma*H'*(H*Sigma*H' + R)^-1*H* Sigma + Sigma
    real*8 :: var_17(n, k)                   !  Sigma*H', 0
    real*8 :: Hvar_2(k, n)                   !  (H*Sigma*H' + R)^-1*H, H
    real*8 :: var_11(n)                      !  0, H'*(H*Sigma*H' + R)^-1*((-1)*data + H*mu)
    real*8 :: var_19(n, n)                   !  0, H'*(H*Sigma*H' + R)^-1*H
    real*8 :: var_5(n, n)                    !  0
    real*8 :: var_20(n, n)                   !  H'*(H*Sigma*H' + R)^-1*H*Sigma, 0

    call dcopy(n**2, var_5, 1, var_20, 1)
    call dsymm('L', 'U', n, k, 1, Sigma, n, H, k, 0, var_17, n)
    call dgemm('N', 'N', k, k, n, 1, H, k, var_17, n, 1, R, k)
    call dcopy(n**2, var_5, 1, var_19, 1)
    call dcopy(n, mu, 1, muvar_2, 1)
    call dcopy(n**2, Sigma, 1, Sigmavar_2, 1)
    call dcopy(k*n, H, 1, Hvar_2, 1)
    call dgemm('N', 'N', k, 1, n, 1, H, k, mu, n, -1, data, k)
    call dposv('U', k, n, R, k, Hvar_2, k, INFO)
    call dposv('U', k, 1, R, k, data, k, INFO)
    call dgemm('N', 'N', n, n, k, 1, H, k, Hvar_2, k, 0, var_19, n)
    call dgemm('N', 'N', n, 1, k, 1, H, k, data, k, 0, var_11, n)
    call dsymm('L', 'U', n, n, 1, var_19, n, Sigma, n, 0, var_20, n)
    call dsymm('L', 'U', n, 1, 1, Sigma, n, var_11, n, 1, muvar_2, n)
    call dsymm('L', 'U', n, n, -1, Sigmavar_2, n, var_20, n, 1, Sigmavar_2, n)

    RETURN
    END
    */

}