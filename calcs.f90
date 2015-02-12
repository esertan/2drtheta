module calcs
    use mydefs

    implicit none

contains

    subroutine make_stretch(r1, r2, r01, r02, xmin, step, ssize)

    real(rk), intent(in) :: r01, r02, xmin, ssize
    real(rk), allocatable, intent(inout) :: r1(:,:), r2(:,:)
    real(rk) :: k, l
    integer(ik) :: step, i, j

    do i=1,step
        k=i*ssize
        do j=1,step
            l=j*ssize

            r1(j,i) = r01+xmin+l
            r2(j,i) = r02+xmin+k
            write(*,*) r1(j,i), r2(j,i), l, k

        end do
    end do
    
    end subroutine make_stretch

    subroutine make_bend(a, a0, amin, astep, asize)

    real(rk), intent(in) :: a0, amin, asize
    real(rk), allocatable, intent(inout) :: a(:)
    real(rk) :: k
    integer(ik) :: i, astep

    do i=1,astep
        k=i*asize

        a(i) = a0+amin+k
        write(*,*) a(i), deg(a(i)) 

    end do
    
    end subroutine make_bend
        

end module calcs
