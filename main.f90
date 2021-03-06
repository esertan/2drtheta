program main

    use mydefs
    use inits
    use calcs
    use writes

    implicit none

    integer(ik) :: maxatoms, step, astep
    real(rk) :: xmin, amin, ssize, asize
    real(rk) :: a0, r01, r02
    real(rk), allocatable :: xyz0(:,:)
    real(rk), allocatable :: r1(:,:)
    real(rk), allocatable :: r2(:,:)
    real(rk), allocatable :: a(:)

    character(len=3), allocatable :: aname(:)
    character(len=5) :: grid='GRIDI', mol='MOLIN', atype='ATOMT', geo='GEOME', fix='FIXED'
    integer :: iflag

    

    iflag = read_key(mol)

    if(iflag==1) then
        
        call get_molinfo(maxatoms)

    end if

    iflag = read_key(grid)

    if(iflag==1) then

        call init_grid(xmin, step, ssize, amin, astep, asize)

    end if


    call allocate_arrays(xyz0, r1, r2, a, maxatoms, step, astep, aname)

    iflag = read_key(atype)

    if(iflag==1) then
        
        call get_atypes(aname, maxatoms)

    end if


    iflag = read_key(geo)

    if(iflag==1) then

       call init_coords(xyz0, maxatoms, r01, r02, a0)

       call make_stretch(r1, r2, r01, r02, xmin, step, ssize)

       call write_rgrid(r1,r2, step)

       call write_stretch(r1, r2, a0, maxatoms, step, aname)

       call make_bend(a, amin, astep, asize)

       call write_agrid(a, astep)

       call write_bend(a, r01, r02, maxatoms, astep, aname)
    end if

    call deallocate_arrays(xyz0, r1, r2, a, aname)


    
end program main
