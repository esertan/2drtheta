module writes
    use mydefs

    implicit none

contains

    subroutine write_rgrid(r1, r2, step)

    !Writes bond distances to file in atomic units, each step has its own grid file.
    
    real(rk), allocatable, intent(in) :: r1(:,:), r2(:,:)
    integer(ik), intent(in) :: step
    integer(ik) :: i,j
    character(len=4) :: qfile

    !open(unit=iw1,form='formatted', status='unknown', file='Q')
    do i=1,step

        select case(i)
            case(1:9)
                write(qfile, '(a,i1)') 'Q00', i
            case(10:99)
                write(qfile, '(a,i2)') 'Q0', i
            case(100:999)
                write(qfile, '(a,i3)') 'Q', i
        end select
        open(unit=iw1,form='formatted', status='unknown', file=qfile)

        do j=1,step
            write(iw1, '(f6.3,3x,f6.3)') r2(j,i), r1(j,i)
        end do
        close(iw1)
    end do

    !close(iw1)
    
    end subroutine write_rgrid

    subroutine write_stretch(r1,r2,a0, maxatoms, step, aname)

    !Writes bond displaced coordinates to file. Output in cartesian format.

    integer(ik) :: maxatoms, step
    real(rk), allocatable, intent(in) :: r1(:,:), r2(:,:)
    real(rk) :: a0
    character(len=3), allocatable, intent(in) :: aname(:)
    character(len=11) :: cfile
    character(len=7) :: dfile
    integer(ik) :: imol, jmol, fc=0

    do jmol=1,step
        do imol=1,step
            fc=fc+1
!            write(*,*) "imol=", imol, "fc=", fc
            select case(fc)
                case(1:9)
                    write(cfile, '(a,i1)') '2DRTHETA00', fc
                case(10:99)
                    write(cfile, '(a,i2)') '2DRTHETA0', fc
                case(100:999)
                    write(cfile, '(a,i3)') '2DRTHETA', fc
            end select

!            write(*,*) 'Writing file:', cfile

            open(unit=iw2,form='formatted',status='unknown',file=cfile)

            !Write coordinates to file in cartesian format

            write(iw2, '(a2,2x,f9.6,2x,f9.6,2x,f9.6)') aname(1), 0.d0, r1(imol,jmol)*sin(0.5_rk*a0)*autoaa, &
                    -r1(imol,jmol)*cos(0.5_rk*a0)*autoaa

            write(iw2,'(a2,2x,f9.6,2x,f9.6,2x,f9.6)') aname(2), 0d0, 0d0, 0d0

            write(iw2,'(a2,2x,f9.6,2x,f9.6,2x,f9.6)') aname(3), 0.d0, -r2(imol,jmol)*sin(0.5_rk*a0)*autoaa,&
                    -r2(imol,jmol)*cos(0.5_rk*a0)*autoaa

            write(iw2,'(a1,3x,f9.6,2x,f9.6,2x,f9.6)') "X", 0.d0, 0.d0, 0.00005d0

            close(iw2)

            !Write coordinates to trajectory file (for append to work remove these files before running the code)

            select case(imol)
                case(1:9)
                    write(dfile, '(a,i1)') 'TRAJ00', imol
                case(10:99)
                    write(dfile, '(a,i2)') 'TRAJ0', imol
                case(100:999)
                    write(dfile, '(a,i3)') 'TRAJ', imol

            end select

            open(unit=iw3,form='formatted',status='unknown',file=dfile, position='append')

            write(iw3,'(i1)') maxatoms
            write(iw3,'(i1)')

            write(iw3,'(a1,2x,f9.6,2x,f9.6,2x,f9.6)') aname(1), 0.d0, r1(imol,jmol)*sin(0.5_rk*a0)*autoaa, &
                    -r1(imol,jmol)*cos(0.5_rk*a0)*autoaa
            write(iw3,'(a1,2x,f9.6,2x,f9.6,2x,f9.6)') aname(2), 0.d0, 0.d0, 0.d0
            write(iw3,'(a1,2x,f9.6,2x,f9.6,2x,f9.6)') aname(3), 0.d0, -r2(imol,jmol)*sin(0.5_rk*a0)*autoaa,&
                    -r2(imol,jmol)*cos(0.5_rk*a0)*autoaa

            close(iw3)
        end do
    end do

    end subroutine write_stretch

end module writes
