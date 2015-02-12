module mydefs

    implicit none

    integer, parameter :: rk = selected_real_kind(15,300)
    integer, parameter :: ik = selected_int_kind(5)
    integer, parameter :: iw1=15, iw2=16, iw3=17, iw4=18
    real(rk), parameter :: pi = 3.141592653589793238462643383279502884197_rk    
    real(rk), parameter :: autoaa = 0.529177_rk

contains

    real(rk) function rad(ang)

    real(rk), intent(inout) :: ang
    real(rk) :: den

    den = 180.000000_rk
    rad = ang * pi / den

    end function rad

    real(rk) function deg(ang)

    real(rk), intent(inout) :: ang
    real(rk) :: den

    den = 180.000000_rk
    deg = ang * den / pi

    end function deg


end module mydefs
