module inits
   
     use mydefs

    implicit none

contains

    integer function read_key(checkkey)

    character(len=*), intent(in) :: checkkey
    character(len=5) :: key    

    read(*,*) key
    write(*,*) key
    
    if(key==checkkey) then
        read_key=1
    else
        read_key=0
    end if

        
    end function read_key


    subroutine allocate_arrays(xyz0, r1, r2, a, maxatoms, step, astep, aname)

    ! Subroutine allocates the arrays for atom positions and dispacements.

    integer :: istat
    integer(ik) :: maxatoms, step, astep
    real(rk), allocatable :: xyz0(:,:)
    real(rk), allocatable :: r1(:,:)
    real(rk), allocatable :: r2(:,:)
    real(rk), allocatable :: a(:)
    character(len=3), allocatable :: aname(:)

    allocate(xyz0(1:maxatoms,1:3), STAT=istat)
    allocate(r1(1:step,1:step), STAT=istat)
    allocate(r2(1:step,1:step), STAT=istat)
    allocate(a(1:astep), STAT=istat)
    allocate(aname(1:maxatoms), STAT=istat)

    end subroutine allocate_arrays

    subroutine deallocate_arrays(xyz0, r1, r2, a, aname)

    ! Subroutine deallocates arrays that were allocated in the allocate_arrays subroutine.

    integer :: istat
    real(rk), allocatable :: xyz0(:,:)
    real(rk), allocatable :: r1(:,:)
    real(rk), allocatable :: r2(:,:)
    real(rk), allocatable :: a(:)
    character(len=3), allocatable :: aname(:)

    deallocate(xyz0, STAT=istat)
    deallocate(r1, STAT=istat)
    deallocate(r2, STAT=istat)
    deallocate(a, STAT=istat)      
    deallocate(aname, STAT=istat)
    ! TODO: Check istat !!

    end subroutine deallocate_arrays    
!
    subroutine get_molinfo(maxatoms)

    ! Subroutine reads in nr atoms from stdin.

    integer(ik), intent(inout) :: maxatoms

    read(*,*) maxatoms

    write(*,*) maxatoms

    end subroutine get_molinfo

    subroutine get_atypes(aname, maxatoms)

    character(len=3), allocatable :: aname(:)
    integer(ik), intent(in) :: maxatoms
    integer(ik) :: i

    read(*,*) (aname(i), i=1,maxatoms)
    write(*,*) (aname(i), i=1,maxatoms)

    end subroutine get_atypes
!
!    subroutine init_coords(xyz0, maxatoms, r01, r02, a0)
    
    ! Subroutine reads in initial geometry from stdin in cartesian coordinates (Angstrom). Converts coordinates to (r1, r2, theta) format.
!    integer(ik), intent(inout) :: maxatoms
!    real(rk), allocatable, intent(inout) :: xyz0(:,:)
!    integer :: i,k
!
!    write(*,*) 'Reading input geometry (a.u)'

!    do i=1, maxatoms
!        read(*,*) (xyz0(i,k), k=1,3)
!        xyz0(i,:) = xyz0(i,:)/autoaa
!!        write(*,*) (xyz0(i,k), k=1,3)
!    end do
!    
!    r01 = sqrt(sum((xyz0(1,:)-xyz0(2,:))**2))
!    r02 = sqrt(sum((xyz0(3,:)-xyz0(2,:))**2))
!    a0 = acos(sum((xyz0(3,:)-xyz0(2,:))*(xyz0(1,:)-xyz0(2,:)))/r01/r02)
!
!
!    end subroutine init_coords

    subroutine init_fixed(r01, r02, a0)

    !Subroutine reads in initial distances (a.u) and angle. Works only for N=3 (no dihedral specified) 
   
    real(rk), intent(inout) :: r01, r02, a0 

    read(*,*) r01, r02, a0

    a0=rad(a0)

    write(*,*) r01, r02, a0

    end subroutine init_fixed

    subroutine init_grid(xmin, step, ssize, amin, astep, asize)

    ! Reads grid parameters from stdin

    integer(ik), intent(inout) :: step, astep
    real(rk), intent(inout) :: xmin, ssize, amin, asize

    read(*,*) xmin, step, ssize
    write(*,*) xmin, step, ssize
    read(*,*) amin, astep, asize
    !Convert astep to radians
    asize=rad(asize)
    
    write(*,*) amin, astep, asize


    end subroutine init_grid

end module inits
