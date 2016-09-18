extern crate blas;

use blas::c::*;

fn main() {

    let n = 5;
    let k = 10;
    let alpha = 1.0;
    let beta = 0.0;

    let Sigma =             vec![1.0; (n * n) as usize];
    let H =                 vec![1.0; (k * n) as usize];
    let mu =                vec![1.0; n as usize];
    let R =                 vec![1.0; (k * k) as usize];
    let data =              vec![1.0; k as usize];
    let Hvar_2 =            vec![1.0; (k * n) as usize];
    let var_11 =            vec![1.0; n as usize];
    let var_19 =            vec![1.0; (n * n) as usize];
    let var_5 =             vec![1.0; (n * n) as usize];

    let mut var_20 =        vec![1.0; (n * n) as usize];
    let mut INFO =          -1;
    let mut muvar_2 =       vec![1.0; (n * n) as usize];
    let mut Sigmavar_2 =    vec![1.0; (n * n) as usize];
    let mut var_17 =        vec![1.0; (n * k) as usize];

    dcopy(
        n*n,
        &var_5,
        1,
        &mut var_20,
        1);

    dsymm(
        Layout::ColumnMajor,
        Side::Left,
        Part::Upper,
        n,
        k,
        alpha,
        &Sigma,
        n,
        &H,
        k,
        beta,
        &mut var_17,
        n);

    
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