! Created by Benjamin Alheit on 2021/06/14.


!module Increments
!
!contains
!
!    subroutine IncrementAndPrintReal(data)
!        character(len=1) :: data(:)
!        real             :: r
!        r = transfer(data, r)
!        r = r + 1.0
!        print *,r
!        data = transfer(r, data)
!    end subroutine
!
!    subroutine IncrementAndPrintInteger(data)
!        character(len=1) :: data(:)
!        integer          :: i
!        i = transfer(data, i)
!        i = i + 1
!        print *,i
!        data = transfer(i, data)
!    end subroutine
!
!    subroutine IncrementTenTimes(incrFunc, data)
!        character(len=1) :: data(:)
!        integer :: i
!        interface
!            subroutine incrFunc(data)
!                character(len=1) :: data(:)
!            end subroutine
!        end interface
!        do i = 1, 10
!            call incrFunc(data)
!        enddo
!    end subroutine
!
!end module



!subroutine call_back_test()
!    print *, "In call back test"
!end
!
!
!subroutine call_back(sub)
!    interface
!        subroutine sub()
!        end subroutine sub
!    end interface
!    call sub()
!end

!TODO Compare lin viscoelastic model with Abaqus

module la
contains
    function inv(A) result(Ainv)
        implicit none
        doubleprecision, intent(in) :: A(:, :)
        doubleprecision :: Ainv(size(A, 1), size(A, 2))
        doubleprecision :: work(size(A, 1))            ! work array for LAPACK
        integer :: n, info, ipiv(size(A, 1))     ! pivot indices

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A, 1)
        ! SGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        !    call SGETRF(n,n,Ainv,n,ipiv,info)
        call DGETRF(n, n, Ainv, n, ipiv, info)
        !    if (info.ne.0) stop 'Matrix is numerically singular!'
        ! SGETRI computes the inverse of a matrix using the LU factorization
        ! computed by SGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)
        if (info.ne.0) stop 'Matrix inversion failed!'
    end function inv


    function invreal(A) result(Ainv)
        implicit none
        real, intent(in) :: A(:, :)
        real :: Ainv(size(A, 1), size(A, 2))
        real :: work(size(A, 1))            ! work array for LAPACK
        integer :: n, info, ipiv(size(A, 1))     ! pivot indices

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A, 1)
        ! SGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call SGETRF(n, n, Ainv, n, ipiv, info)
        if (info.ne.0) stop 'Matrix is numerically singular!'
        ! SGETRI computes the inverse of a matrix using the LU factorization
        ! computed by SGETRF.
        call SGETRI(n, Ainv, n, ipiv, work, n, info)
        if (info.ne.0) stop 'Matrix inversion failed!'
    end function invreal

end module la


module solid_mechanics


contains



subroutine lin_elastic(sig_s_val, sig_s_params, eps_s, n_s)
    integer, intent(in) :: n_s
    doubleprecision, intent(in) :: eps_s(6), sig_s_params(n_s)
    doubleprecision :: mu, lambda
    doubleprecision, intent(out) :: sig_s_val(6)

    lambda = sig_s_params(1)
    mu = sig_s_params(2)

    sig_s_val = 0d0
    sig_s_val(1:3) = sum(eps_s(1:3)) * lambda
    sig_s_val = sig_s_val + 2.d0 * mu * eps_s

end subroutine lin_elastic

subroutine iso_lin_elastic(sig_s_val, sig_s_params, eps_s, n_s)
    integer, intent(in) :: n_s
    doubleprecision, intent(in) :: eps_s(6), sig_s_params(n_s)
    doubleprecision :: mu
    doubleprecision, intent(out) :: sig_s_val(6)

    mu = sig_s_params(1)

    sig_s_val = 2.d0 * mu * eps_s

end subroutine iso_lin_elastic

subroutine iso_lin_elastic_sig_d_depss(d_sig_s_depss_val, sig_s_params, eps_s, n_s)
    integer, intent(in) :: n_s
    doubleprecision, intent(in) :: eps_s(6), sig_s_params(n_s)
    doubleprecision :: mu
    integer :: m
    doubleprecision, intent(out) :: d_sig_s_depss_val(6, 6)

    mu = sig_s_params(1)
    d_sig_s_depss_val = 0d0
    forall (m = 1:6) d_sig_s_depss_val(m, m) = 2 * mu

end subroutine iso_lin_elastic_sig_d_depss

subroutine lin_visc(sig_d_val, sig_d_params, eps_d_dot_bar, n_d)
    integer, intent(in) :: n_d
    doubleprecision, intent(in) :: eps_d_dot_bar(6), sig_d_params(n_d)
    doubleprecision :: eta
    doubleprecision, intent(out) :: sig_d_val(6)

    eta = sig_d_params(1)
    sig_d_val = eta * eps_d_dot_bar

end subroutine lin_visc

subroutine lin_viscd_sig_d_depsd_dot(lin_viscd_sig_d_depsd_dot_val, sig_d_params, eps_d_dot_bar, n_d)
    integer, intent(in) :: n_d
    doubleprecision, intent(in) :: eps_d_dot_bar(6), sig_d_params(n_d)
    doubleprecision :: eta
    integer :: m
    doubleprecision, intent(out) :: lin_viscd_sig_d_depsd_dot_val(6, 6)

    eta = sig_d_params(1)
    lin_viscd_sig_d_depsd_dot_val = 0d0
    forall (m = 1:6) lin_viscd_sig_d_depsd_dot_val(m, m) = eta

end subroutine lin_viscd_sig_d_depsd_dot

subroutine increment_eps_bar(sig_s, dsig_s_depse, sig_s_params, n_s, &
        sig_d, dsig_d_depsd_dot, sig_d_params, n_d, &
        eps_n1, eps_d_n, eps_d_n1, sig_dn1, dt)
    use la

    integer, intent(in) :: n_s, n_d
    doubleprecision, intent(in) :: eps_n1(6), eps_d_n(6), sig_s_params(n_s), sig_d_params(n_d), dt
!    doubleprecision, intent(in) :: eps_d_n(6)

    doubleprecision :: res(6), res_norm, eps_d_dot(6)

    doubleprecision :: d_sig_d_depsd_dot_val(6, 6), d_sig_s_depss_val(6, 6), drdeps(6, 6)

    doubleprecision :: sig_s_val(6)

    doubleprecision, intent(out) :: eps_d_n1(6), sig_dn1(6)

    external :: sig_s, dsig_s_depse, sig_d, dsig_d_depsd_dot

!    interface
!        subroutine sig_s(sig_s_val, sig_s_params, eps_e, n_s)
!            integer, intent(in) :: n_s
!            doubleprecision, intent(in) :: eps_e(6), sig_s_params(n_s)
!            doubleprecision, intent(out) :: sig_s_val(6)
!        end subroutine sig_s
!        subroutine dsig_s_depse(dsig_s_depse_val, sig_s_params, eps_e, n_s)
!            integer, intent(in) :: n_s
!            doubleprecision, intent(in) :: eps_e(6), sig_s_params(n_s)
!            doubleprecision, intent(out) :: dsig_s_depse_val(6, 6)
!        end subroutine dsig_s_depse
!        subroutine sig_d(sig_d_val, sig_d_params, eps_d_dot, n_d)
!            integer, intent(in) :: n_d
!            doubleprecision, intent(in) :: eps_d_dot(6), sig_d_params(n_d)
!            doubleprecision, intent(out) :: sig_d_val(6)
!        end subroutine sig_d
!        subroutine dsig_d_depsd_dot(dsig_d_depse_val, sig_d_params, eps_d_dot, n_d)
!            integer, intent(in) :: n_d
!            doubleprecision, intent(in) :: eps_d_dot(6), sig_d_params(n_d)
!            doubleprecision, intent(out) :: dsig_d_depse_val(6, 6)
!        end subroutine dsig_d_depsd_dot
!!        function inv(A) result(Ainv)
!!            implicit none
!!            doubleprecision, intent(in) :: A(:, :)
!!            doubleprecision :: Ainv(size(A, 1), size(A, 2))
!!            doubleprecision :: work(size(A, 1))            ! work array for LAPACK
!!            integer :: n, info, ipiv(size(A, 1))     ! pivot indices
!!        end function inv
!    end interface

    !    res_norm = 1
    eps_d_n1 = eps_n1

    eps_d_dot = (eps_d_n1 - eps_d_n) / dt
    call sig_s(sig_s_val, sig_s_params, eps_n1 - eps_d_n1, n_s)
    call sig_d(sig_dn1, sig_d_params, eps_d_dot, n_d)

    !    sig_dn1 = eta * eps_d_dot_bar
    !    sig_s_val = mu * (eps_n1 - eps_d_n1)

    res = sig_dn1 - sig_s_val
    res_norm = dot_product(res, res)

    do while (res_norm > 1e-6)

        !        call lin_viscd_sig_d_depsd_dot(d_sig_d_depsd_dot_val, [eta], eps_d_dot, 1)
        !        call iso_lin_elastic_sig_d_depss(d_sig_s_depss_val, [mu], eps_d_dot_bar, 1)
        call dsig_d_depsd_dot(d_sig_d_depsd_dot_val, sig_d_params, eps_d_dot, n_d)
        call dsig_s_depse(d_sig_s_depss_val, sig_s_params, eps_n1 - eps_d_n1, n_s)

        drdeps = d_sig_d_depsd_dot_val / dt + d_sig_s_depss_val
        eps_d_n1 = eps_d_n1 - matmul(inv(drdeps), res)

        eps_d_dot = (eps_d_n1 - eps_d_n) / dt
        call sig_d(sig_dn1, sig_d_params, eps_d_dot, n_d)
        call sig_s(sig_s_val, sig_s_params, eps_n1 - eps_d_n1, n_s)
        !        sig_dn1 = eta * eps_d_dot_bar
        !        sig_s_val = mu * (eps_n1 - eps_d_n1)

        res = sig_dn1 - sig_s_val
        res_norm = dot_product(res, res)

    end do

end

subroutine increment_eps_bar_lin(eta, mu, eps_n1, eps_d_n, eps_d_n1, sig_dn1, dt)

    use la
    !    integer, intent(in) :: n_s, n_d
    !    doubleprecision, intent(in) :: eps_n1(6), sig_s_params(n_s), sig_d_params(n_d), dt
    doubleprecision, intent(in) :: eta, mu, dt, eps_n1(6)
    doubleprecision, intent(in) :: eps_d_n(6)

    doubleprecision :: res(6), res_norm, eps_d_dot_bar(6)

    doubleprecision :: d_sig_d_depsd_dot_val(6, 6), d_sig_s_depss_val(6, 6), drdeps(6, 6), drdeps_inv(6, 6)

    doubleprecision :: sig_s_iso(6)

    doubleprecision, intent(out) :: eps_d_n1(6), sig_dn1(6)

!    interface
!        function inv(A) result(Ainv)
!            implicit none
!            doubleprecision, intent(in) :: A(:, :)
!            doubleprecision :: Ainv(size(A, 1), size(A, 2))
!            doubleprecision :: work(size(A, 1))            ! work array for LAPACK
!            integer :: n, info, ipiv(size(A, 1))     ! pivot indices
!        end function inv
!    end interface

    !    res_norm = 1
    eps_d_n1 = eps_n1

    eps_d_dot_bar = (eps_d_n1 - eps_d_n) / dt
    sig_dn1 = eta * eps_d_dot_bar
    sig_s_iso = 2*mu * (eps_n1 - eps_d_n1)
    res = sig_dn1 - sig_s_iso
    res_norm = dot_product(res, res)

    do while (res_norm > 1e-6)

        call lin_viscd_sig_d_depsd_dot(d_sig_d_depsd_dot_val, [eta], eps_d_dot_bar, 1)
        call iso_lin_elastic_sig_d_depss(d_sig_s_depss_val, [mu], eps_d_dot_bar, 1)
        drdeps = d_sig_d_depsd_dot_val / dt + d_sig_s_depss_val
        drdeps_inv = inv(drdeps)
        eps_d_n1 = eps_d_n1 - matmul(drdeps_inv, res)

        eps_d_dot_bar = (eps_d_n1 - eps_d_n) / dt
        sig_dn1 = eta * eps_d_dot_bar
        sig_s_iso = 2*mu * (eps_n1 - eps_d_n1)
        res = sig_dn1 - sig_s_iso
        res_norm = dot_product(res, res)

    end do

end

end module

subroutine eye_voigt(I)
    doubleprecision, intent(inout) :: I(6)
    I(1:3) = 1d0
    I(4:6) = 0d0
end subroutine eye_voigt

function trace_voigt(voigt_vec) result (out)
    doubleprecision, intent(in) :: voigt_vec(6)
    integer :: i
    doubleprecision :: out

    out = 0d0
    forall (i = 1:3) out = out + voigt_vec(i)

end function trace_voigt

function norm_voigt(voigt_vec) result (out)
    doubleprecision, intent(in) :: voigt_vec(6)
    integer :: i
    doubleprecision :: out

    out = 0d0
    forall (i = 1:3) out = out + voigt_vec(i) * voigt_vec(i)
    forall (i = 4:6) out = out + 2 * voigt_vec(i) * voigt_vec(i)

    out = out ** 0.5d0

end function norm_voigt


subroutine print_matrix_real(mat, rows, collumns)
    integer :: rows, collumns
    real :: mat(rows, collumns)
    integer :: i, j

    do i = 1, rows
        print*, mat(i, :)
    enddo
end

subroutine print_matrix(mat, rows, collumns)
    integer :: rows, collumns
    double precision :: mat(rows, collumns)
    integer :: i, j

    do i = 1, rows
        print*, mat(i, :)
    enddo
end


program main
    !    use test_callback
    use solid_mechanics

    doubleprecision :: eps_d_n1(6), sig_dn1(6), a
!
!
!
!    !    character(len=1), allocatable :: data(:)
!    !    integer                       :: lengthData
!    !    real                          :: r = 5.0
!    !    integer                       :: i = 10
!
!    !    call increment_eps_bar_lin(sig_s, dsig_s_depse, sig_s_params, n_s, &
!    !    sig_d, dsig_d_depsd_dot, sig_d_params, n_d, &
!    !    eps_n1, eps_d_n, eps_d_n1, sig_dn1, dt)
!
        call increment_eps_bar_lin(1.d0, 2.d0, &
                [2.d0, 2.d0, 2.d0, 0.d0, 0.d0, 0.d0], &
                [1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0], &
                eps_d_n1, &
                sig_dn1, &
                0.1d0)

        a=0
!
        call increment_eps_bar(iso_lin_elastic, iso_lin_elastic_sig_d_depss, [2d0], 1, &
                lin_visc, lin_viscd_sig_d_depsd_dot, [1d0], 1, &
                [2.d0, 2.d0, 2.d0, 0.d0, 0.d0, 0.d0], &
                [1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0], &
                eps_d_n1, sig_dn1, 0.1d0)

        a = 0
!!    call call_back(call_back_test)
!
end program
