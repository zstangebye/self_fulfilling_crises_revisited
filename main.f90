!  Code for:
!  'Rollover Crises Revisited: The Art of the Desperate Deal'
!  Authors:
!  Satyajit Chatterjee, Harold Cole, Mark Aguiar, and Zach Stangebye
!  September 2018
!
!  In this code, zm corresponds to exp(g) in the paper, and
!  shock s corresponds to shock z in the paper.

program main

    use globvar
    use funcall

    implicit none

    integer, parameter::NN=4
    real(doub),dimension(NN)::funcinp
    real(doub):: res

    ! If running in parallel, initialize relevant MPI variables: Else, set them to serial values:

    if(mpi_on) then
        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,id,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)
    else
        id = 0
        nproc = 1
    end if

    open(200,file="n_procs.txt")

    write(200,*) nproc

    close(200)

    ! Initiate random matrices used in simulations:

    call random_number(shockzw) ! Output/wealth shock
    call random_number(shockh)  ! Shock determines re-inclusion into financial markets
    call random_number(shocks)  ! Sovereign BC shock to help with computation
    call random_number(shockp)  ! Sunspot shock

    i1 = 0
    do iz=1,Z
        do iw  = 1,W
            i1 = i1+1
            zwloop_z(i1) = iz
            zwloop_w(i1) = iw
            iizw(iz,iw)  = i1
        end do
    end do

    !smthresh: the end points for the intervals of s

    dels        = 2*m2/S*sds
    smthresh(1) = smean-m2*sds

    do i = 2,S+1
        smthresh(i) = smthresh(i-1)+dels
    end do

    do i = 1,S
        sm(i) = (smthresh(i)+smthresh(i+1))*0.5
    end do

    do i = 1,S
        smdiff(i) = smthresh(i+1)-smthresh(i)
    end do

    !pdfs: the probability of drawing s (completely independent of z and w)

    do j = 1,S
        call cdfnor(xx,x1,non,dble(smthresh(j+1)),dble(0),dble(sds),sta,boun)
        call cdfnor(xx,x2,non,dble(smthresh(j)),dble(0),dble(sds),sta,boun)
        pdfs(j) = x1-x2
    end do

    divids = spread(sum(pdfs),DIM=1,NCOPIES=S)
    pdfs   = pdfs/divids

    ! The following loop vectorizes the state space into one long vector (assets, shocks, wealth)

    i1 = 0
    do ia = 1,A
        do iz = 1,Z
            do iw = 1,W
                i1 = i1+1
                izloop(i1) = iz
                ialoop(i1) = ia
                iwloop(i1) = iw
                iiv(ia,iz,iw) = i1
            end do
        end do
    end do


    nii  = int(real(nzaw-1,8)/real(nproc,8))+1
    itop = id*nii+1
    iend = min((id+1)*nii,nzaw)


!    write (temp_name,'(i2)') id
!
!    open(50,file=trim(temp_name))
!
!    write(50,*) itop
!    write(50,*) iend
!
!    close(50)


    ! Allocate objects that may be too big to initiate on compilation


    allocate(val1(nii), val1agg(nii*nproc), val2(nii),val2agg(nii*nproc), &
        val_int1(nii*S), val_int1agg(nii*S*nproc), val_int2(nii*S), val_int2agg(nii*S*nproc), val_int3(nii*S), val_int3agg(nii*S*nproc))
    allocate(val3(nii*A*S), val3agg(nii*A*S*nproc), val4(nii*A*S), val4agg(nii*A*S*nproc), val5(nii*A*S), val5agg(nii*A*S*nproc) )

    allocate(q(Z,W,A),q2(Z,W,A))

    aminvec = aminset

    kkk = 1

    funcinp(1) = beta_true
    funcinp(2) = defp0
    funcinp(3) = w_base

    ticker = 0

    if(GMM_on) then ! If we are matching moments, solve the model many times
        n  = 4 !dimension of problem
        mm = 6

        iprint = 1 ! specify display

        factr = 1.0e7 ! specify tolerance
        pgtol = 1.0e-5
        grad_tol = 1.0e-3

        nbd(1) = 2 ! Specify bounds on variables
        nbd(2) = 2
        nbd(3) = 2
        nbd(4) = 2

        l(1) = 0.8
        l(2) = .06
        l(3) = 2.0
        l(4) = 0.05

        u(1) = 0.85
        u(2) = 0.08
        u(3) = 3.5
        u(4) = 0.95


        x(1:4) = funcinp(1:4)


        task = 'START' ! initialize algorithm

        111  continue

        !    This is the call to the L-BFGS-B code.

        call setulb(n,mm,x,l,u,nbd,f,gg,factr,pgtol,wa,iwa,task,iprint,&
        csave,lsave,isave,dsave)

        if (task(1:2) .eq. 'FG') then
            !        the minimization routine has returned to request the
            !        function f and gradient g values at the current x

            !        Compute function value f
            ticker = ticker + 1
            f = funcmin(x)

            !        Compute gradient g numerically

            do min_ind = 1,n
                ticker  = ticker + 1
                xtemp   = x
                xtemp(min_ind) = xtemp(min_ind) + grad_tol
                gg(min_ind)    = (funcmin(xtemp)-f)/grad_tol
            end do

            temp_str = 'current_x_'
            temp_name = trim(temp_str) // tagg
            open(200,file=trim(temp_name))

            write(200,'(f12.9,2X,f12.9,2X,f12.9,2x,f12.9)') x(1:4)

            close(200)

            !          go back to the minimization routine.
            goto 111
        endif
        !
        if (task(1:5) .eq. 'NEW_X')  goto 111
        !        the minimization routine has returned with a new iterate,
        !         and we have opted to continue the iteration.

        !           ---------- the end of the loop -------------

        !     If task is neither FG nor NEW_X we terminate execution.

    else ! otherwise, solve it once
        res = funcmin(funcinp)
    end if

    ! write the solutions (only relevant if matching moments)

    temp_str = 'solution_'
    temp_name = trim(temp_str) // tagg
    open(20,file=trim(temp_name))

    write(20,'(f12.9,2X,f12.9,2X,f12.9,2x,f12.9)') x(1:4)

    close(20)

    if(mpi_on) then
        call mpi_finalize(ierr)
    end if

    100 format (10000(1x, f12.8))
    110 format (1000(1x,i4))
    120 format (1000(1x, f7.3))


end program main
