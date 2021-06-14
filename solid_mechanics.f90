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

!program main
!    use Increments
!    character(len=1), allocatable :: data(:)
!    integer                       :: lengthData
!    real                          :: r = 5.0
!    integer                       :: i = 10
!
!    lengthData = size(transfer(r, data))
!    allocate(data(lengthData))
!    data = transfer(r, data)
!    call IncrementTenTimes(IncrementAndPrintReal, data)
!    deallocate(data)
!
!    lengthData = size(transfer(i, data))
!    allocate(data(lengthData))
!    data = transfer(i, data)
!    call IncrementTenTimes(IncrementAndPrintInteger, data)
!
!end program

subroutine lin_elastic(sig_s_val, sig_s_params, eps_s, n_s)
    integer, intent(in) :: n_s
    doubleprecision, intent(in) :: eps_s(6), sig_s_params(n_s)
    doubleprecision :: mu, lambda
    doubleprecision, intent(out) :: sig_s_val(6)

    lambda = sig_s_params(1)
    mu = sig_s_params(2)

    sig_s_val = 0d0
    sig_s_val(1:3) = sum(eps_s(1:3)) * lambda
    sig_s_val =  sig_s_val + 2.d0 * mu * eps_s

end subroutine lin_elastic

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
    forall (m=1:6) lin_viscd_sig_d_depsd_dot_val(m, m) = eta

end subroutine lin_viscd_sig_d_depsd_dot

subroutine increment_eps_bar(sig_s, dsig_s_depse, sig_s_params, n_s, &
                             sig_d, dsig_d_depsd_dot, sig_d_params, n_d, &
                             eps_n1, eps_d_n, eps_d_n1, sig_dn1)

    integer, intent(in) :: n_s, n_d
    doubleprecision, intent(in) :: eps_n1(6), sig_s_params(n_s), sig_d_params(n_d)
    doubleprecision, intent(inout) :: eps_d_n(6)

    doubleprecision, intent(out) :: eps_d_n1(6), sig_dn1(6)

    interface
        subroutine sig_s(sig_s_val, sig_s_params, eps_e, n_s)
            integer, intent(in) :: n_s
            doubleprecision, intent(in) :: eps_e(6), sig_s_params(n_s)
            doubleprecision, intent(out) :: sig_s_val(6)
        end subroutine sig_s
        subroutine dsig_s_depse(dsig_s_depse_val, sig_s_params, eps_e, n_s)
            integer, intent(in) :: n_s
            doubleprecision, intent(in) :: eps_e(6), sig_s_params(n_s)
            doubleprecision, intent(out) :: dsig_s_depse_val(6, 6)
        end subroutine dsig_s_depse
        subroutine sig_d(sig_d_val, sig_d_params, eps_d_dot, n_d)
            integer, intent(in) :: n_d
            doubleprecision, intent(in) :: eps_d_dot(6), sig_d_params(n_d)
            doubleprecision, intent(out) :: sig_d_val(6)
        end subroutine sig_d
        subroutine dsig_d_depse(dsig_d_depse_val, sig_d_params, eps_d_dot, n_d)
            integer, intent(in) :: n_d
            doubleprecision, intent(in) :: eps_d_dot(6), sig_d_params(n_d)
            doubleprecision, intent(out) :: dsig_d_depse_val(6, 6)
        end subroutine dsig_d_depse
    end interface

end subroutine increment_eps_bar

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
    forall (i=1:3) out = out + voigt_vec(i)

end function trace_voigt