
module funcall

    use globvar

    implicit none

    contains

    real(doub) function funcmin (xx0)
        implicit none


        real(doub), intent(in), dimension(4) :: xx0

        beta_true=xx0(1)
        defp0 =xx0(2)
        w_base = xx0(3)


        ! Discretize the state-space with Tauchen's method

        lrsdz = sdz               ! stdev of invariant distribution of log prod
        delz  = 2*m/(Z-1)*lrsdz   ! length of intervals between z in discretizations (equal distance)
        zm(1) = muz-m*lrsdz

        do  i = 2,Z
            zm(i) = zm(i-1)+delz
        end do

        INFIN(1)  = 2
        INFIN(2)  = 2
        correl(1) = sigma_zw
        correl(2) = 1.0

        if(W > 1) then

            lrsdw = sdw             !stdev of invariant distribution of
            num   = (W-1)
            delw  = 2*m/num*lrsdw   !length of intervals between w in discretizations (equal distance)
            wm(1) = muw-m*lrsdw

            do i = 2,W
                wm(i) = wm(i-1)+delw
            end do

        else

            lrsdw = 10
            delw  = 10
            wm(1) = muw

        end if

        do i1 = 1,NZW                                       ! Joint state of (z,w) today
            iz   = zwloop_z(i1)
            iw   = zwloop_w(i1)

            do i2 = 1,NZW                                     ! Joint state of (z',w') tomorrow
                iz2   = zwloop_z(i2)
                iw2   = zwloop_w(i2)

                lb(1) = ( (zm(iz2) - delz/2.0) - rhoz*zm(iz) - (1-rhoz)*muz )/(sdz*sqrt(1-rhoz**2)) ! MVNDST gives results for standard normal; must normalize here
                ub(1) = ( (zm(iz2) + delz/2.0) - rhoz*zm(iz) - (1-rhoz)*muz )/(sdz*sqrt(1-rhoz**2)) ! Divide by stdev(error) = stdev(uncond)*sqrt(1-rhoz**2) < stdev(uncond)
                lb(2) = ( (wm(iw2) - delw/2.0) -  rhow*wm(iw) - (1-rhow)*muw )/max(sdw*sqrt(1-rhow**2),1e-10)
                ub(2) = ( (wm(iw2) + delw/2.0) -  rhow*wm(iw) - (1-rhow)*muw )/max(sdw*sqrt(1-rhow**2),1e-10)

                call MVNDST(DN,lb,ub,INFIN,correl,MXPNTS,ABSEPS,RELEPS,ERR,int_val,INFORM)
                pdfzw(i1,i2) = int_val ! Assume independence in sunspot process
            end do
        end do


        dividzw = spread(sum(pdfzw,DIM=2),DIM=2,NCOPIES=NZW)
        pdfzw   = pdfzw/dividzw

        logzm = zm
        zm = exp(zm)

        if (beta_true*zm(1)**(1-gamma) .GE. 1.0) then
            write(*,*) 'Beta too large! Likely not stationary'
        end if

        ! Used in taking of expectations
        do i = 1,W
            beta_aug((i-1)*Z+1:i*Z) = beta_true*zm**(1-gamma)
        end do

        logwm = wm
        wm    = w_base + exp(wm) - 1

        !cdfzw: cdf of z and w

        cdfzw(:,1) = pdfzw(:,1)
        do i = 2,NZW
            cdfzw(:,i) = cdfzw(:,i-1)+pdfzw(:,i)
        end do
        cdfzw(:,NZW) = one


        amin = aminset

        ! output shock z under default
        ! Comes from non-normalized: C = Y(1-defp0-defp1*exp(g))
        ! Divide by Y to derive:     c = 1-defp0-defp1*exp(g)
        ! i.e. Higher g implies relatively lower consumption post-default
        do iz=1,Z
            defz(iz)=zbase - defp0 - defp1*zm(iz) ! Asymmetric default costs
        end do


        ! flow utility under default
        do iz=1,Z
            do iw = 1,W

                izw = iizw(iz,iw)

                do is=1,S
                    def_fix(is,izw)= sov_utility(sm(is)+defz(iz)) ! Not a function of w, but makes it easier to take expectations later
                end do

            end do

        end do

        ! asset grids


        !        dela1=(amax-amax_grid_break)/real(course_grid_A-1,doub)
        !        am(A-course_grid_A+1)=amax_grid_break
        !        am(A)=amax
        !        do i=A-course_grid_A+2,A-1
        !            am(i)=am(i-1)+dela1
        !        end do
        !
        !        dela1=(amax_grid_break-amin)/real(A-course_grid_A+1,doub)
        !        am(1)=amin
        !        do i=2,A-course_grid_A
        !            am(i)=am(i-1)+dela1
        !        end do
        !
        !        ambound=A
        !
        !        if (id .EQ. 0) then
        !
        !            ! Print debt grid
        !            temp_str = 'debt_grid'
        !            temp_name = trim(temp_str) // tagg
        !            open(20,file=trim(temp_name))
        !
        !            do i1=1,A
        !                write(20,formatter) am(i1)
        !            end do
        !
        !            close(20)
        !
        !        end if


        dela1=(amax-amin)/real(A-1,doub)
        am(1)=amin
        am(A)=amax
        do i=2,A-1
            am(i)=am(i-1)+dela1
        end do

        ! HERE, WE INTERPRET DEBT AS SCALED BY SIZE OF GOVERNMENT (only do this
        ! if zbase = 1.0)
          am = am
          aminset = aminset

        ambound=A

        ! epsilon grids

        epsilon_delta =(epsilon_upper_bar-epsilon_lower_bar)/real(EPS_CHECK-1,doub)
        epsilon_grid(1)=epsilon_lower_bar
        epsilon_grid(EPS_CHECK) = epsilon_upper_bar
        do i=2,EPS_CHECK-1
            epsilon_grid(i)=epsilon_grid(i-1)+epsilon_delta
        end do

        ! Here, allow for V-sunspot to be increasing in b
        ! spt(:,2) = am*spt_prob2/amin

        do iz=1,Z
            do ia=1,A
                zamvec(iz,ia)=zbase+(lambda+coup)*am(ia)/zm(iz) ! Output less debt obligations (consumption absent borrowing and s-shock)
            end do
        end do


        temp_str = 'debt_grid'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))
        do ia=1,A

            write(20,formatter) zamvec(:,ia)

        end do
        close(20)

        amnx=A
        do ia=1,A
            num=(one-lambda)*am(ia)
            i=A-1
            do ia2=A-1,1,-1
                if (num>=am(ia2)) then
                    i=ia2
                    exit
                end if
            end do
            amnx(ia)=i
            ampr(ia)=(am(i+1)-num)/(am(i+1)-am(i))
        end do

        ! start with risk-free price
        qbase=(lambda+coup)/(rbase+lambda) ! Always equals 1 when rbase=coup

        q=qbase
        !q=0.0


        errtot_q=0
        errtot_v=0

        vn = -5.0
        vds= -5.0
        vns= -5.0


        do t=1,time

            ! Begin iteration by computing value of default and continuation value of repayment
            ! given last iteration's computed repayment valuations

            ! First compute the value of default
            do izw=1,NZW
                Vd(izw)=dot_product(pdfs,Vds(:,izw)) ! State-by-state value of default
            end do

            Vd_fut=matmul(pdfzw,beta_aug*Vd) ! Expected value of default (continuation)

            do ia=1,A
                do iz=1,Z
                    do iw = 1,W

                        i1 = iiv(ia,iz,iw)
                        izw = iizw(iz,iw)

                        zvec(izw) = Vn(i1) ! State-by-state continuation repayment value (all states tomorrow)

                    end do
                end do
                do iz=1,Z
                    do iw = 1,W

                        izw = iizw(iz,iw)
                        Vn_fut(ia,iz,iw) = dot_product(pdfzw(izw,:),beta_aug*zvec) ! Expected continuation value

                    end do
                end do
            end do

            ! Full value of default, including contemporaneous consumption
            do iz = 1,Z
                do iw = 1,W

                    izw = iizw(iz,iw)

                    num=ent*Vn_fut(A,iz,iw)+(one-ent)*Vd_fut(izw) ! Full expected continuation value of default

                    do is=1,S
                        Vds(is,izw)=def_fix(is,izw)+num ! Full value of default (for a shock s; every subsequent period)
                    end do

                end do
            end do


            ! Now, compute optimal choice of contemporaneous borrowing

            do i1=itop,iend ! Begins parallel execution (itop & iend differ by core)

                shortS=thrmin ! Vector of thresholds between optimal asset choices
                defthr=thrmin ! Vector of default thresholds

                iz=izloop(i1) ! Output today
                ia=ialoop(i1) ! Debt level today
                iw=iwloop(i1) ! Creditor wealth today

                iins=i1-itop+1

                do ia2=1,A
                    paym(ia2) = -q(iz,iw,ia2)*(am(ia2)-(1.0-lambda)*am(ia)/zm(iz)) ! This is the EG price
                end do

                izw = iizw(iz,iw)


                i_sa_counter = 1
                ! optimization routine (also pins down price schedules across beliefs)
                do isind=1,S

                    stop_mix_fix1 = .FALSE.
                    stop_mix_fix2 = .FALSE.
                    do ia2=1,A

                        if( paym(ia2) > 0) then
                            if(  ( sov_utility(sm(isind) + zamvec(iz,ia) + paym(ia2)) + Vn_fut(ia2,iz,iw) ) .GE. ( sov_utility(sm(isind) + zbase + (1.0-lend_share)*paym(ia2) ) + Vd_fut(izw) + epsilon_upper_bar ) ) then
                                ! always repay current epsilon-shock does not matter in this case
                                ! Checked EG price here
                                util_vecOPT(ia2) = sov_utility( sm(isind) + zamvec(iz,ia) + paym(ia2) ) + Vn_fut(ia2,iz,iw)
                                epsOPT(isind,ia2) = epsilon_upper_bar
                                qOPT(isind,ia2) = q(iz,iw,ia2)


                                if( (sov_utility( sm(isind) + zamvec(iz,ia) ) + Vn_fut(ia2,iz,iw)) < ( sov_utility( sm(isind) + zbase ) + Vd_fut(izw) + epsilon_lower_bar ) ) then
                                    ! Now check if Cole-Kehoe can happen regardless of epsilon-shock (CASE 1)
                                    ! Checked 0 price here
                                    util_vecPES(ia2) = sov_utility( sm(isind) + zbase ) + Vd_fut(izw) + (epsilon_upper_bar + epsilon_lower_bar)/2.0
                                    epsPES(isind,ia2) = epsilon_lower_bar
                                    qPES(isind,ia2) = 0.0

                                    ! In this case, we also know that concerned equilibrium exists (under uniform, it will be unique)
                                    call interval_bisection_eps(epsilon_lower_bar,epsilon_upper_bar)

                                    epsCON(isind,ia2) = eps_star
                                    util_vecCON(ia2) = cdf_epsilon(eps_star)*(sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2))+Vn_fut(ia2,iz,iw)) + &
                                    (1.0-cdf_epsilon(eps_star))*( sov_utility(sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) ) + (eps_star+epsilon_upper_bar)/2.0*(epsilon_upper_bar-eps_star)/(epsilon_upper_bar-epsilon_lower_bar)
                                    qCON(isind,ia2) = FPfun(eps_star,ia,iz,ia2)*q(iz,iw,ia2)

                                    con_prob_vec(ia2) = cdf_epsilon(eps_star)
                                    con_repay_val_vec(ia2) = sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2))+Vn_fut(ia2,iz,iw)
                                    con_def_val_vec(ia2) = sov_utility(sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw)
                                    con_eps_shock_cont(ia2) = (eps_star+epsilon_upper_bar)/2.0*(epsilon_upper_bar-eps_star)/(epsilon_upper_bar-epsilon_lower_bar)


                                else ! In this case, by concavity we know two worst equilibria are same as optimistic (CASE 2)
                                    ! Case in which pess=conc and are all zero is a knife-edge case

                                    util_vecPES(ia2) = util_vecOPT(ia2)
                                    epsPES(isind,ia2) = epsOPT(isind,ia2)
                                    qPES(isind,ia2) = qOPT(isind,ia2)

                                    util_vecCON(ia2) = util_vecOPT(ia2)
                                    epsCON(isind,ia2) = epsOPT(isind,ia2)
                                    qCON(isind,ia2) = qOPT(isind,ia2)

                                    con_prob_vec(ia2) = 1.0
                                    con_repay_val_vec(ia2) = sov_utility( sm(isind) + zamvec(iz,ia) + paym(ia2) ) + Vn_fut(ia2,iz,iw)
                                    con_def_val_vec(ia2) = sov_utility(sm(isind) + zbase + (1.0-lend_share)*paym(ia2) ) + Vd_fut(izw)
                                    con_eps_shock_cont(ia2) = .5*(epsilon_lower_bar+epsilon_upper_bar)


                                end if

                            elseif(  (sov_utility(sm(isind) + zamvec(iz,ia) + paym(ia2)) + Vn_fut(ia2,iz,iw)) < ( sov_utility(sm(isind) + zbase + (1.0-lend_share)*paym(ia2)) + Vd_fut(izw) + epsilon_lower_bar ) ) then
                                ! always default (price=zero under all beliefs) (CASE 3)
                                ! Checked EG price here
                                util_vecOPT(ia2) = sov_utility(sm(isind) + zbase ) + Vd_fut(izw) + .5*(epsilon_lower_bar+epsilon_upper_bar)
                                epsOPT(isind,ia2) = epsilon_lower_bar
                                qOPT(isind,ia2) = 0.0

                                util_vecPES(ia2) = util_vecOPT(ia2)
                                epsPES(isind,ia2) = epsOPT(isind,ia2)
                                qPES(isind,ia2) = qOPT(isind,ia2)

                                util_vecCON(ia2) = util_vecOPT(ia2)
                                epsCON(isind,ia2) = epsOPT(isind,ia2)
                                qCON(isind,ia2) = qOPT(isind,ia2)

                                con_prob_vec(ia2) = 0.0
                                con_repay_val_vec(ia2) = sov_utility(sm(isind) + zamvec(iz,ia) ) + Vn_fut(ia2,iz,iw)
                                con_def_val_vec(ia2) = sov_utility(sm(isind) + zbase ) + Vd_fut(izw)
                                con_eps_shock_cont(ia2) = .5*(epsilon_lower_bar+epsilon_upper_bar)


                            ! Here, we do not repay for every epsilon-shock regardless of beliefs, but we also do not default for every epsilon-shock regardless of beliefs
                            elseif( ( sov_utility( sm(isind) + zamvec(iz,ia) ) + Vn_fut(ia2,iz,iw) ) < ( sov_utility( sm(isind) + zbase ) + Vd_fut(izw) + epsilon_lower_bar ) ) then

                                ! In this case, we have the pure-strategy rollover for sure; other two equilibria may or may not be rollover (CASES 4.1 and 4.2)
                                ! Case in which optimistic and concerned are same is degenerate
                                util_vecPES(ia2) = sov_utility(sm(isind) + zbase ) + Vd_fut(izw) + (epsilon_upper_bar + epsilon_lower_bar)/2.0 ! mean of epsilon-shock
                                epsPES(isind,ia2) = epsilon_lower_bar
                                qPES(isind,ia2) = 0.0

                                ! Assume they all result in default; check and see if they ever don't
                                util_vecCON(ia2) = sov_utility(sm(isind) + zbase ) + Vd_fut(izw) + (epsilon_upper_bar + epsilon_lower_bar)/2.0 ! mean of epsilon-shock
                                epsCON(isind,ia2) = epsilon_lower_bar
                                qCON(isind,ia2) = 0.0

                                util_vecOPT(ia2) = sov_utility(sm(isind) + zbase ) + Vd_fut(izw) + (epsilon_upper_bar + epsilon_lower_bar)/2.0 ! mean of epsilon-shock
                                epsOPT(isind,ia2) = epsilon_lower_bar
                                qOPT(isind,ia2) = 0.0

                                con_prob_vec(ia2) = 0.0
                                con_repay_val_vec(ia2) = sov_utility(sm(isind) + zamvec(iz,ia) )+Vn_fut(ia2,iz,iw)
                                con_def_val_vec(ia2) = sov_utility(sm(isind) + zbase ) + Vd_fut(izw)
                                con_eps_shock_cont(ia2) = .5*(epsilon_lower_bar+epsilon_upper_bar)


                                ! Check over all epsilons to see when/if they switch sides
                                eps_check_not_done = .TRUE.
                                ieps = 1


                                current_eps = epsilon_grid(ieps)
                                do while( (ieps <= EPS_CHECK) .AND. eps_check_not_done)

                                    eps_flag = sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(current_eps,ia,iz,ia2)*paym(ia2)) + Vn_fut(ia2,iz,iw) - &
                                    ( sov_utility( sm(isind) + zbase + (1.0-lend_share)*FPfun(current_eps,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) + current_eps)
                                    if( eps_flag > 0.0 ) then


                                        if(ieps > 1) then

                                            call interval_bisection_eps(epsilon_grid(ieps-1),current_eps)

                                        else

                                            call interval_bisection_eps(epsilon_grid(ieps)/100.0,current_eps) ! don't want to start at lowest bound (since this is an equilibrium for sure)

                                        end if

                                        epsCON(isind,ia2) = eps_star
                                        util_vecCON(ia2) = cdf_epsilon(eps_star)*(sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2)) + Vn_fut(ia2,iz,iw)) + &
                                        (1.0-cdf_epsilon(eps_star))*( sov_utility(sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) ) + (eps_star+epsilon_upper_bar)/2.0*(epsilon_upper_bar-eps_star)/(epsilon_upper_bar-epsilon_lower_bar)
                                        qCON(isind,ia2) = FPfun(eps_star,ia,iz,ia2)*q(iz,iw,ia2)

                                        con_prob_vec(ia2) = cdf_epsilon(eps_star)
                                        con_repay_val_vec(ia2) = sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2))+Vn_fut(ia2,iz,iw)
                                        con_def_val_vec(ia2) = sov_utility(sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw)
                                        con_eps_shock_cont(ia2) = (eps_star+epsilon_upper_bar)/2.0*(epsilon_upper_bar-eps_star)/(epsilon_upper_bar-epsilon_lower_bar)


                                        call interval_bisection_eps(current_eps,epsilon_upper_bar)

                                        epsOPT(isind,ia2) = eps_star
                                        util_vecOPT(ia2) = cdf_epsilon(eps_star)*(sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2)) + Vn_fut(ia2,iz,iw)) + &
                                        (1.0-cdf_epsilon(eps_star))*( sov_utility(sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) ) + (eps_star+epsilon_upper_bar)/2.0*(epsilon_upper_bar-eps_star)/(epsilon_upper_bar-epsilon_lower_bar)
                                        qOPT(isind,ia2) = FPfun(eps_star,ia,iz,ia2)*q(iz,iw,ia2)

                                        eps_check_not_done = .FALSE.

                                    end if

                                    ieps = ieps + 1

                                end do

                            else
                                ! If we cannot have a pure-strategy rollover, then by concavity there is a unique interior solution for all three (CASE 5)
                                ! Checked zero price here

                                call interval_bisection_eps(epsilon_lower_bar,epsilon_upper_bar)

                                epsOPT(isind,ia2) = eps_star
                                util_vecOPT(ia2) = cdf_epsilon(eps_star)*(sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2))+Vn_fut(ia2,iz,iw)) + &
                                (1.0-cdf_epsilon(eps_star))*( sov_utility(sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) ) + (eps_star+epsilon_upper_bar)/2.0*(epsilon_upper_bar-eps_star)/(epsilon_upper_bar-epsilon_lower_bar)
                                qOPT(isind,ia2) = FPfun(eps_star,ia,iz,ia2)*q(iz,iw,ia2)

                                util_vecPES(ia2) = util_vecOPT(ia2)
                                epsPES(isind,ia2) = epsOPT(isind,ia2)
                                qPES(isind,ia2) = qOPT(isind,ia2)

                                util_vecCON(ia2) = util_vecOPT(ia2)
                                epsCON(isind,ia2) = epsOPT(isind,ia2)
                                qCON(isind,ia2) = qOPT(isind,ia2)

                                con_prob_vec(ia2) = cdf_epsilon(eps_star)
                                con_repay_val_vec(ia2) = sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2))+Vn_fut(ia2,iz,iw)
                                con_def_val_vec(ia2) = sov_utility(sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw)
                                con_eps_shock_cont(ia2) = (eps_star+epsilon_upper_bar)/2.0*(epsilon_upper_bar-eps_star)/(epsilon_upper_bar-epsilon_lower_bar)

                            end if

                        else
                            ! In this case, we have a buyback, which we know implies a unique equilibrium. Just need to see what
                            ! threshold eps_star is and then we know everything

                            if( (sov_utility( sm(isind) + zamvec(iz,ia) + paym(ia2) ) + Vn_fut(ia2,iz,iw) ) >= &
                            (sov_utility( sm(isind) + zbase + (1.0-lend_share)*paym(ia2) ) + Vd_fut(izw) + epsilon_upper_bar) ) then

                                ! Here, even the best default shock induces repayment
                                eps_star = epsilon_upper_bar

                            elseif( (sov_utility(sm(isind) + zamvec(iz,ia)) + Vn_fut(ia2,iz,iw) ) < &
                            (sov_utility( sm(isind) + zbase) + Vd_fut(izw) + epsilon_lower_bar) ) then

                                ! Here, even the worst default shock induces default
                                eps_star = epsilon_lower_bar

                            else
                                ! Here, we have the interior-equilibrium case (know it must be unique by monotonicity)

                                call interval_bisection_eps(epsilon_lower_bar, epsilon_upper_bar)

                            end if

                            epsPES(isind,ia2) = eps_star
                            util_vecPES(ia2) = cdf_epsilon(eps_star)*(sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2))+Vn_fut(ia2,iz,iw)) + &
                            (1.0-cdf_epsilon(eps_star))*( sov_utility(sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) ) + (eps_star+epsilon_upper_bar)/2.0*(epsilon_upper_bar-eps_star)/(epsilon_upper_bar-epsilon_lower_bar)
                            qPES(isind,ia2) = FPfun(eps_star,ia,iz,ia2)*q(iz,iw,ia2)

                            epsCON(isind,ia2) = epsPES(isind,ia2)
                            util_vecCON(ia2) = util_vecPES(ia2)
                            qCON(isind,ia2) = qPES(isind,ia2)

                            epsOPT(isind,ia2) = epsPES(isind,ia2)
                            util_vecOPT(ia2) = util_vecPES(ia2)
                            qOPT(isind,ia2) = qPES(isind,ia2)

                        end if


                        ! Here, set `corner' prices to EG under assumption that epsilon-shock small: Fixes grid issues

                        if( mixing_fix .AND. (qCON(isind,ia2) < qOPT(isind,ia2)) .AND. (qCON(isind,ia2) > qPES(isind,ia2)) ) then

                            if( ia2 > 1 ) then ! First check overborrowing case (if not borrowing maximum amount already)

                                if( (qCON(isind,ia2-1) .LE. 1.0e-5) .OR. (qCON(isind,ia2-1) .GE. q(iz,iw,ia2-1) - 1.0e-5) ) then

                                    if(.NOT. stop_mix_fix1) then ! Make sure to only do this once for all possible debt choices

                                        ! Close enough; set to EG price
                                        epsCON(isind,ia2) = epsilon_upper_bar
                                        qCON(isind,ia2) = q(iz,iw,ia2)
                                        util_vecCON(ia2) = util_vecOPT(ia2)
                                        stop_mix_fix1 = .TRUE.

                                    end if

                                elseif( ia2 < A ) then  ! If not overborrowing, check consolidation case (if not borrowing zero already)

                                    if( (qCON(isind,ia2+1) .LE. 1.0e-5) .OR. (qCON(isind,ia2+1) .GE. q(iz,iw,ia2+1) - 1.0e-5) ) then

                                        if(.NOT. stop_mix_fix2) then ! Make sure to only do this once for all possible debt choices

                                            ! Close enough; set to EG price
                                            epsCON(isind,ia2) = epsilon_upper_bar
                                            qCON(isind,ia2) = q(iz,iw,ia2)
                                            util_vecCON(ia2) = util_vecOPT(ia2)
                                            stop_mix_fix2 = .TRUE.

                                        end if

                                    end if

                                end if

                            end if

                        end if
                        ! End mixing-fix routine; comment out if not wanted (or set mixing_fix to .FALSE.)

                        val3( (iins-1)*S*A+i_sa_counter) = qOPT(isind,ia2)
                        val4( (iins-1)*S*A+i_sa_counter) = qPES(isind,ia2)
                        val5( (iins-1)*S*A+i_sa_counter) = qCON(isind,ia2)

                        i_sa_counter = i_sa_counter + 1

                    end do

                    dummy = maxloc(util_vecOPT)

                    apolOPT(isind) = dummy(1)
                    ucons(isind) = util_vecOPT(dummy(1))

                    ! Determine how much was actually transferred from lenders to sovereign
                    if(am(apolOPT(isind)) < (1.0-lambda)*am(ia)/zm(iz)) then
                    auc_rev = cdf_epsilon(epsOPT(isind,apolOPT(isind)) )*paym(ia2)
                    else
                        auc_rev = 0.0
                    end if

                    q_fut(isind) = cdf_epsilon(epsOPT(isind,apolOPT(isind)) )*(lambda+coup+(one-lambda)*(qOPT(isind,apolOPT(isind)))) + &
                    (1.0 - cdf_epsilon(epsOPT(isind,apolOPT(isind)) ) )*( lend_share*(lambda+coup)*auc_rev/(-(lambda + coup)*am(ia)/zm(iz)-am(apolOPT(isind))) + (1.0-lambda)*qOPT(isind,apolOPT(isind)) ) ! Account for pro-rata recovery

                    dummy = maxloc(util_vecPES)

                    apolPES(isind) = dummy(1)
                    ucons2(isind) = util_vecPES(dummy(1))

                    ! In cases where default is certain, restrict issuance (even though indifferent)
                    if( cdf_epsilon(epsPES(isind,apolPES(isind)) ) .EQ. 0.0 ) then
                        apolPES(isind) = A
                        ucons2(isind) = util_vecPES(A)
                    end if

                    ! Determine how much was actually transferred from lenders to sovereign
                    if(am(apolPES(isind)) < (1.0-lambda)*am(ia)/zm(iz)) then
                    auc_rev = cdf_epsilon(epsPES(isind,apolPES(isind)) )*paym(ia2)
                    else
                        auc_rev = 0.0
                    end if

                    q_fut2(isind) = cdf_epsilon(epsPES(isind,apolPES(isind)) )*(lambda+coup+(one-lambda)*(qPES(isind,apolPES(isind)))) + &
                    (1.0 - cdf_epsilon(epsPES(isind,apolPES(isind)) ) )*( lend_share*(lambda+coup)*auc_rev/(-(lambda + coup)*am(ia)/zm(iz)-am(apolPES(isind))) + (1.0-lambda)*qPES(isind,apolPES(isind)) )

                    dummy = maxloc(util_vecCON)

                    apolCON(isind) = dummy(1)
                    ucons3(isind) = util_vecCON(dummy(1))

                    ! Check if super-close to indifference is throwing off maximization
                    if( show_obj .AND. (ia .EQ. 70) .AND. (iz .EQ. 8) .AND. (isind .EQ. 6) ) then

                        temp_str = 'con_objective_'
                        temp_name = trim(temp_str) // tagg
                        open(20,file=trim(temp_name))

                        write(20,formatter) util_vecCON
                        write(20,formatter) con_prob_vec
                        write(20,formatter) con_repay_val_vec
                        write(20,formatter) con_def_val_vec
                        write(20,formatter) con_eps_shock_cont

                        close(20)

                    end if

                    ! Determine how much was actually transferred from lenders to sovereign
                    if(am(apolCON(isind)) < (1.0-lambda)*am(ia)/zm(iz)) then
                    auc_rev = cdf_epsilon(epsCON(isind,apolCON(isind)) )*paym(ia2)
                    else
                        auc_rev = 0.0
                    end if
                    !
                    q_fut3(isind) = cdf_epsilon(epsCON(isind,apolCON(isind)) )*(lambda+coup+(one-lambda)*(qCON(isind,apolCON(isind)))) + &
                    (1.0 - cdf_epsilon(epsCON(isind,apolCON(isind)) ) )*( lend_share*(lambda+coup)*auc_rev/(-(lambda + coup)*am(ia)/zm(iz)-am(apolCON(isind))) + (1.0-lambda)*qCON(isind,apolCON(isind)) )


                end do

                qpr = opt_bel*dot_product(pdfs,q_fut) + pes_bel*dot_product(pdfs,q_fut2) + con_bel*dot_product(pdfs,q_fut3)

                uexpt = opt_bel*dot_product(pdfs,ucons) + pes_bel*dot_product(pdfs,ucons2) + con_bel*dot_product(pdfs,ucons3)


                ! Recall iins=i1-itop+1

                ! Write value functions and policy thresholds into MPI-readable vectors
                val1(iins) = uexpt
                val2(iins) = qpr

                val_int1( (iins-1)*S+1:iins*S ) = apolOPT
                val_int2( (iins-1)*S+1:iins*S ) = apolPES
                val_int3( (iins-1)*S+1:iins*S ) = apolCON

            end do

            ! Make sure all processes stop to wait here
            if(mpi_on) then

                call MPI_BARRIER(mpi_comm_world, ierr)

            end if

            ! Aggregate price and value functions from cluster operations

            ! only price and value function will be called during solution of the model
            if(mpi_on) then
                call mpi_allgather(val1(1),nii,mpi_double_precision,&
                val1agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
                Vn2=val1agg(1:nzaw)

            else
                Vn2 = val1
            end if


            if(mpi_on) then
                call mpi_allgather(val2(1),nii,mpi_double_precision,&
                val2agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
                qsum=val2agg(1:nzaw)

            else
                qsum = val2
            end if

            q2 = 0

            ! Take an expectation to get the updated EG price
            do ia=1,A ! Asset choice for tomorrow
                do iz=1,Z ! Today's output
                    do iw = 1,W ! Today's wealth

                        izw = iizw(iz,iw)
                        qdenom = 0.0
                        q_num = 0.0

                        do iz2 = 1,Z ! Tomorrow's potential output
                            do iw2 = 1,W ! Tomorrow's potential wealth

                                i1=iiv(ia,iz2,iw2)
                                izw2 = iizw(iz2,iw2)

                                q_num = q_num + pdfzw(izw,izw2)*qsum(i1)

                            end do
                        end do

                        q2(iz,iw,ia) = q_num/(1.0+rbase)

                    end do
                end do
            end do


            ! Each of these has the same number of elements; Vn is a long vector and q is a 3D array
            errq=maxval(abs(q-q2))
            errv=maxval(abs(Vn-Vn2))

            if(errv < converge_tol) then
                qsmooth = qsmooth2
            end if

            q=qsmooth*q+(one-qsmooth)*q2
            Vn=vsmooth*Vn+(one-vsmooth)*Vn2

            errtot_q(t)=errq
            errtot_v(t)=errv

            ! Write the state to check as it runs
            temp_str = 'tracker_'
            temp_name = trim(temp_str) // tagg
            open(20,file=trim(temp_name)) !20 is the file number (can be anything)
            write(20,'(i4)') t
            write(20,'(3f8.5)') errtot_q(t), errtot_v(t)
            close(20)

            if (t > tmin .AND. errq < converge_tol .AND. errv < converge_tol) exit   ! no-deals convergence check

        end do

        ! only EG price and value function will be called during solution of the model
        ! Now, we pull out policy and other price schedules
        if(mpi_on) then

            call mpi_allgather(val_int1(1),nii*S,mpi_integer,&
            val_int1agg(1),nii*S,mpi_integer,mpi_comm_world,ierr)
            apolVEC=val_int1agg(1:NZAWS)

        else
            apolVEC = val_int1
        end if


        do iz=1,Z
            do iw=1,W
                do ia=1,A

                    i1 = iiv(ia,iz,iw)

                    apolOPTfull(iz,iw,ia,:) = apolVEC( (i1-1)*S+1:i1*S )

                end do
            end do
        end do


        if(mpi_on) then

            call mpi_allgather(val_int2(1),nii*S,mpi_integer,&
            val_int2agg(1),nii*S,mpi_integer,mpi_comm_world,ierr)
            apolVEC=val_int2agg(1:NZAWS)

        else
            apolVEC = val_int2
        end if

        do iz=1,Z
            do iw=1,W
                do ia=1,A

                    i1 = iiv(ia,iz,iw)

                    apolPESfull(iz,iw,ia,:) = apolVEC( (i1-1)*S+1:i1*S )

                end do
            end do
        end do

        if(mpi_on) then

            call mpi_allgather(val_int3(1),nii*S,mpi_integer,&
            val_int3agg(1),nii*S,mpi_integer,mpi_comm_world,ierr)
            apolVEC=val_int3agg(1:NZAWS)

        else
            apolVEC = val_int3
        end if

        do iz=1,Z
            do iw=1,W
                do ia=1,A

                    i1 = iiv(ia,iz,iw)

                    apolCONfull(iz,iw,ia,:) = apolVEC( (i1-1)*S+1:i1*S )

                end do
            end do
        end do



        ! Now for the pricing schedules
        if(mpi_on) then

            call mpi_allgather(val3(1),nii*S*A,mpi_double_precision,&
            val3agg(1),nii*S*A,mpi_double_precision,mpi_comm_world,ierr)
            qlongVEC=val3agg(1:NZAWSA)

        else
            qlongVEC = val3
        end if

        do iz=1,Z
            do iw=1,W
                do ia=1,A

                    i1 = iiv(ia,iz,iw)

                    i2 = 1

                    do isind=1,S
                        do ia2 = 1,A
                            qOPTfull(iz,iw,ia,isind,ia2) = qlongVEC( (i1-1)*S*A+i2 )
                            i2 = i2+1
                        end do

                    end do

                end do
            end do
        end do


        if(mpi_on) then

            call mpi_allgather(val4(1),nii*S*A,mpi_double_precision,&
            val4agg(1),nii*S*A,mpi_double_precision,mpi_comm_world,ierr)
            qlongVEC=val4agg(1:NZAWSA)

        else
            qlongVEC = val4
        end if

        do iz=1,Z
            do iw=1,W
                do ia=1,A

                    i1 = iiv(ia,iz,iw)

                    i2 = 1

                    do isind=1,S
                        do ia2 = 1,A
                            qPESfull(iz,iw,ia,isind,ia2) = qlongVEC( (i1-1)*S*A+i2 )
                            i2 = i2+1

                        end do

                    end do

                end do
            end do
        end do


        if(mpi_on) then

            call mpi_allgather(val5(1),nii*S*A,mpi_double_precision,&
            val5agg(1),nii*S*A,mpi_double_precision,mpi_comm_world,ierr)
            qlongVEC=val5agg(1:NZAWSA)

        else
            qlongVEC = val5
        end if

        do iz=1,Z
            do iw=1,W
                do ia=1,A

                    i1 = iiv(ia,iz,iw)

                    i2 = 1

                    do isind=1,S
                        do ia2 = 1,A
                            qCONfull(iz,iw,ia,isind,ia2) = qlongVEC( (i1-1)*S*A+i2 )
                            i2 = i2+1
                        end do

                    end do

                end do
            end do
        end do

        if( id .EQ. 0) then

            call write_functions

        end if

        if( id .EQ. 0) then
            ! Here, we write out all the constituent elements that determine the concerned price for a particular state
            i1 = 1
            do ia2=1,A

                do iz=1,Z

                    zvec(iz) = Vn(i1) ! State-by-state continuation repayment value (all states tomorrow)
                    i1 = i1 + 1

                end do

                do iz=1,Z

                    Vn_fut(ia2,iz,1) = dot_product(pdfzw(iz,:),beta_aug*zvec) ! Expected continuation value

                end do

            end do


            ia = A - 120 + 1
            iz = int(Z/2)+1
            isind = int(S/2)+1

            do ia2=1,A

                const_vec1(ia2) = sov_utility(sm(isind) + zamvec(iz,ia) - qCONfull(iz,1,ia,isind,ia2)*(am(ia2)-(1.0-lambda)*am(ia)/zm(iz)) )
                const_vec2(ia2) = sov_utility(sm(isind) + zbase - (1.0-lend_share)*qCONfull(iz,1,ia,isind,ia2)*(am(ia2)-(1.0-lambda)*am(ia)/zm(iz)) )
                const_vec3(ia2) = Vn_fut(ia2,iz,1)
                const_vec4(ia2) = (epsilon_upper_bar-epsilon_lower_bar)*qCONfull(iz,1,ia,isind,ia2)/q(iz,1,ia2) + epsilon_lower_bar

            end do

            temp_str = 'constituent_parts_'
            temp_name = trim(temp_str) // tagg
            open(20,file=trim(temp_name))

            write(20,formatter) qCONfull(iz,1,ia,isind,:)
            write(20,formatter) const_vec1
            write(20,formatter) const_vec2
            write(20,formatter) const_vec3
            write(20,formatter) const_vec4

            close(20)

        end if

        maxT = t
        if(maxT < time) then

            conv_yes(kkk) = .TRUE.

        else

            conv_yes(kkk) = .FALSE.

        end if


        ! debt service
        do iz=1,Z
            do j=1,A
                amser(iz,j)=-am(j)*(lambda+coup)/zm(iz) ! Sum together the fraction of debt matured and the fraction of debt with only coupons
            end do
        end do

        ! debt buyback
        do iz = 1,Z
            do i=1,A
                do j=1,A
                    amdiff(iz,i,j)=am(j)-(one-lambda)*am(i)/zm(iz) ! New debt less debt still yet to mature
                end do
            end do
        end do

        ! Simulate the sample economy

        call model_simul

        ! This subroutine will print out key model moments int he funcreg vector as follows

        !funcreg(kkk,1)=average debtFV/GDP
        !funcreg(kkk,2)=average spread (annual)
        !funcreg(kkk,3)= crisis percentile debtPmt/TaxRev
        !funcreg(kkk,4)=annual default frequency
        !funcreg(kkk,5) = spread volatility
        !funcreg(kkk,6) = C/Y volatility
        !funcreg(kkk,7) = (trade balance volatility)/(output volatility)
        !funcreg(kkk,8) = correlation(output,consumption)
        !funcreg(kkk,9) = correlation(output,trade balance)
        !funcreg(kkk,10) = correlation(output,spreads)

        !!!!!!!!!!!Output!!!!!!!!!!!


        temp_str = 'accs_results_v2_'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))

        write(20,'(A,i3,2x,A,i3,2x,A,i3,2x,A,i3)') 'Z=', Z, 'A=', A, 'W=',W,'ticker= ', ticker
        write(20,'(A,2X,f7.4)') 'g_mean= ', muz
        write(20,'(A,2X,f6.3)') 'g_vol = ', sdz
        write(20,'(A,2X,f6.3)') 'g_pers = ', rhoz
        write(20,'(A,2X,f6.3)') 'sds = ', sds
        write(20,'(A,2X,f7.4)') 're-entry probability= ', ent
        write(20,'(A,2X,f12.9)') 'def_lev= ', defp0
        write(20,'(A,2X,f12.9)') 'def_curv= ', defp1
        write(20,'(A,2X,f12.9)') 'disc_fact= ', beta_true
        write(20,'(A,2X,2f12.9)') 'epsilon_bnds= ', epsilon_lower_bar, epsilon_upper_bar
        write(20,'(A,2X,f12.9)') 'default_lender_loss_share= ', lend_share
        write(20,'(A,2X,f12.9)') 'optimistic_belief_prob= ', opt_bel
        write(20,'(A,2X,f12.9)') 'pessimistic_belief_prob= ', pes_bel
        write(20,'(A,2X,f12.9)') 'concerned_belief_prob= ', con_bel


        write(20,'(2X)')

        write(20,'(A)') 'Maturity structure'
        write(20,'(1(A,f8.4))') 'lambda=', lambda
        write(20,'(1(A,f8.4))') 'coupon=', coup

        write(20,'(2X)')

        write(20,'(A)') 'Convergence Tools'
        write(20,'(1(A,f8.4))') 'qsmooth=', qsmooth
        write(20,'(1(A,f8.4))') 'vsmooth=', vsmooth
        write(20,'(1(A,f8.4))') 'FV/QuarterlyTax UpperBnd=', -aminset

        write(20,'(2X)')

        Write(20,'(A)') 'Results'
        Write(20,'(2X)')
        i = 1
        if(conv_yes(i)) then
            write(20,'(A,2X,i6.4,X,A)') 'Successful Convergence:', maxT, 'iterations'

            !            funcmin=(abs(funcreg(i,2) - spreaddata) + abs(funcreg(i,3) - debtdata) + &
            !            abs(funcreg(i,4) - def_freq_data) + abs( funcreg(i,5) - sprvoldata ) )/GMM_coarseness !Weight volatility more heavily
            funcmin = ( abs(funcreg(i,3) - debtdata) + abs( funcreg(i,5) - sprvoldata ) )/GMM_coarseness !Weight volatility more heavily

            !		funcmin=(abs(funcreg(i,2) - 0.034) + abs(funcreg(i,3) - 0.656) + &
            !            abs(funcreg(i,4) - 0.02) + abs( funcreg(i,11) - 0.26 ) + &
            !		abs( funcreg(i,5) - 0.025 ))/GMM_coarseness

        if(isnan(funcmin)) funcmin = 1000
        else
            write(20,'(A)') 'No Convergence'
            ! Return a large value if it does not converge
            funcmin=1000
        end if
        write(20,'(A,2X,f8.5)') 'Average Spread:', funcreg(i,2)
        write(20,'(A,2X,f8.5)') 'Average FV/QuarterlyGDP:', funcreg(i,1)*tax_rev_share
        write(20,'(A,2X,f8.5)') 'FV/AnnualTaxRev 97%ile::', funcreg(i,3)/(4*(lambda + coup))
        write(20,'(A,2X,f8.5)') 'DebtPmt/QuarterlyTaxRev 97%ile:', funcreg(i,3)
        write(20,'(A,2X,f8.5)') 'Default Frequency:', funcreg(i,4)
        write(20,'(A,2X,f8.5)') 'Std Dev of Spread:', funcreg(i,5)
        write(20,'(A,2X,f8.5)') 'stdev(c)/stdev(y):', funcreg(i,6)
        write(20,'(A,2X,f8.5)') 'stdev(NX/y)/stdev(y):', funcreg(i,7)
        write(20,'(A,2X,f8.5)') 'corr(y,c):', funcreg(i,8)
        write(20,'(A,2X,f8.5)') 'corr(NX/y,y):', funcreg(i,9)
        write(20,'(A,2X,f8.5)') 'corr(r-r_f,y):', funcreg(i,10)
        write(20,'(2X)')
        write(20,'(A,2X,f13.10)') 'Absolute Distance from Targets:', funcmin*GMM_coarseness

        close(20)


    end function funcmin



    ! This subroutine simply prints policy and value functions

    subroutine write_functions
        implicit none

        ! Print value function
        temp_str = 'CEC'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))

        do ia=1,A
            do iz=1,Z
                do iw = 1,W

                    i1 = iiv(ia,iz,iw)

                    write(20,formatter) ( (1-gamma)*(1.0-beta_true)*Vn(i1) )**(1.0/(1.0 - gamma))

                end do
            end do
        end do

        close(20)

        ! Print value function
        temp_str = 'CEC_default'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))

        do izw=1,NZW

            write(20,formatter) ( (1-gamma)*(1.0-beta_true)*Vd(izw) )**(1.0/(1.0 - gamma))

        end do

        close(20)


        ! Print benchmark and artificial pricing schedules
        temp_str = 'price_function_EG_'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))

        do i1=1,NZW

            iz = zwloop_z(i1)
            iw = zwloop_w(i1)
            write(20,formatter) q(iz,iw,:)
        end do

        close(20)

        temp_str = 'price_function_OPT_'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))
        do ia=1,A

            write(20,formatter) qOPTfull(int(Z/2)+1,1,ia,int(S/2)+1,:)

        end do
        close(20)

        temp_str = 'price_function_PES_'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))
        do ia=1,A

            write(20,formatter) qPESfull(int(Z/2)+1,1,ia,int(S/2)+1,:)

        end do
        close(20)

        temp_str = 'price_function_CON_'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))
        do ia=1,A

            write(20,formatter) qCONfull(int(Z/2)+1,1,ia,int(S/2)+1,:)

        end do
        close(20)

        temp_str = 'policy_function_OPT_'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))
        do iz=1,Z
            do iw=1,W
                do isind=1,S
                    do ia=1,A

                        if( (apolOPTfull(iz,iw,ia,isind) .GE. 1) .AND. (apolOPTfull(iz,iw,ia,isind) .LE. A) ) then

                            apolTEMP(ia) = am( apolOPTfull(iz,iw,ia,isind) )

                            if(apolTEMP(ia) < aminset + 1e-3) then
                                apolTEMP(ia) = 0.0
                            end if

                        else

                            apolTEMP(ia) = 0.0

                        end if

                    end do

                    write(20,formatter) apolTEMP

                end do
            end do
        end do
        close(20)

        temp_str = 'policy_function_PES_'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))
        do iz=1,Z
            do iw=1,W
                do isind=1,S
                    do ia=1,A

                        if( (apolPESfull(iz,iw,ia,isind) .GE. 1) .AND. (apolPESfull(iz,iw,ia,isind) .LE. A) ) then

                            apolTEMP(ia) = am(apolPESfull(iz,iw,ia,isind) )

                            if(apolTEMP(ia) < aminset + 1e-3) then
                                apolTEMP(ia) = 0.0
                            end if

                        else

                            apolTEMP(ia) = 0.0

                        end if

                    end do

                    write(20,formatter) apolTEMP

                end do
            end do
        end do
        close(20)

        temp_str = 'policy_function_CON_'
        temp_name = trim(temp_str) // tagg
        open(20,file=trim(temp_name))
        do iz=1,Z
            do iw=1,W
                do isind=1,S
                    do ia=1,A

                        if( (apolCONfull(iz,iw,ia,isind) .GE. 1) .AND. (apolCONfull(iz,iw,ia,isind) .LE. A) ) then

                            apolTEMP(ia) = am(apolCONfull(iz,iw,ia,isind) )

                            if(apolTEMP(ia) < aminset + 1e-3) then
                                apolTEMP(ia) = 0.0
                            end if

                        else

                            apolTEMP(ia) = 0.0

                        end if

                    end do

                    write(20,formatter) apolTEMP

                end do
            end do
        end do
        close(20)

    end subroutine write_functions

    ! This subroutine uses the global policy and pricing functions to simulate the model
    ! It alters a global vector of moments (funcreg) that will be printed and can be used to match
    ! moments. It assumes the s-shock is part of the endowment process

    subroutine model_simul
        implicit none

        DOUBLE PRECISION dinvnr
        EXTERNAL dinvnr

        !! Simulate the economy and generate sample moments
        shock_path(TB) = NZW/2 ! governs which belief gets selected
        zpath(TB)=zwloop_z(shock_path(TB)) ! initiate zpath.
        wpath(TB)=zwloop_w(shock_path(TB)) ! initiate wpath.
        crisis_path(TB) = 1 ! 1 is optimistic, 2 is concerned, 3 is pessimistic
        apath(TB) = A
        history(TB) = 1
        first_pd_def(TB) = 0

        def_RO = 0
        def_total = 0

        temp_str = 'sample_sim_filter_'
        temp_name = trim(temp_str) // tagg
        open(30,file=trim(temp_name))
        do ss=1,sim

            ! shock paths
            call random_number(shockone)
            do i = 1,3

                if( (shockone .GE. 0.0) .AND. (shockone < opt_bel) ) then

                    crisis_path(1) = 1

                elseif( shockone < (opt_bel + pes_bel) ) then

                    crisis_path(1) = 2

                else

                    crisis_path(1) = 3

                end if

            end do

            shock_path(1) = shock_path(TB)
            zpath(1) = zpath(TB)
            wpath(1) = wpath(TB)
            crisis_path(1) = crisis_path(TB)
            do t = 2,TB
                temp = cdfzw(shock_path(t-1),:) !!!!!!!!!!!!!!!!!
                do izw = 1,NZW
                    if (shockzw(t,ss) .le. temp(izw)) then

                        shock_path(t) = izw

                        zpath(t)=zwloop_z(izw)
                        wpath(t)=zwloop_w(izw)
                        exit

                    end if
                end do

                call random_number(shockone)

                if( (shockone .GE. 0.0) .AND. (shockone < opt_bel) ) then

                    ! optimistic beliefs
                    crisis_path(t) = 1

                elseif( shockone < opt_bel + pes_bel ) then

                    ! pessimistic beliefs
                    crisis_path(t) = 2

                else

                    ! concerned beliefs
                    crisis_path(t) = 3

                end if


            end do

            ! Create spath (iid) using an inverse normal cdf to generate random normal deviates
            do i = 1,TB
                num=shocks(i,ss)
                spath(i)=dinvnr(shocks(i,ss),one-shocks(i,ss))*sds
                if ( spath(i)>smthresh(S+1) .OR. spath(i)<smthresh(1) ) then ! If they fall outside the boundary, pick new ones until they don't
                    defpr=.true.
                    do while (defpr)
                        call random_number(shockone)
                        spath(i)=dinvnr(shockone,one-shockone)*sds
                        if ( spath(i)<=smthresh(S+1) .AND. spath(i)>=smthresh(1) )then
                            defpr=.false.
                        end if
                    end do
                end if
            end do

            apath(1)=apath(TB) !initiating assets (no debt initially)
            history(1)=history(TB) ! history: 1 if included in financial markets, 0 if excluded from financial markets.
            first_pd_def(1) = first_pd_def(TB)

            do t=1,TB-1

                ! First, locate nearest s-shock and assume that for policy functions
                dummy = minloc( (sm - spath(t))**2.0 )
                spath_ind(t) = dummy(1)


                if (history(t)==1 .OR. (history(t)==0 .AND. shockh(t,ss)<ent .AND. first_pd_def(t)==0)) then    ! if included in fin markets or redeemed into fin markets.

                    ! Beliefs matter for price here since in fin markets
                    call random_number(shockone) ! This is the epsilon shock (transformed through cdf into [0,1])

                    defpr = .FALSE. ! Assume no default; check if that's correct (consolidates code a bit)

                    if(crisis_path(t) .EQ. 1) then

                        ! HERE, WE HAVE OPTIMISTIC BELIEFS

                        if(apolOPTfull(zpath(t),wpath(t),apath(t),spath_ind(t)) .EQ. A) then

                            !  Here, creditors offer zero price even for zero issuance (necessitates default today)
                            rep_prob = 0.0
                            apath(t+1) = A
                            qpath(t) = qbase

                        else

                            ! Ratio of prices gives repayment probability
                            rep_prob = qOPTfull(zpath(t),wpath(t),apath(t),spath_ind(t),apolOPTfull(zpath(t),wpath(t),apath(t),spath_ind(t)))/q(zpath(t),wpath(t),apolOPTfull(zpath(t),wpath(t),apath(t),spath_ind(t)))
                            apath(t+1) = apolOPTfull(zpath(t),wpath(t),apath(t),spath_ind(t))
                            qpath(t)= max(qOPTfull(zpath(t),wpath(t),apath(t),spath_ind(t),apath(t+1)),eps)

                        end if


                        !                        ! Eliminate interior defaults (if we're rigging the game)
                        !                        if( (rep_prob > 1.0e-6) .AND. (rep_prob < 1.0- 1.0e-6) .AND. load_the_die ) then
                        !
                        !                            rep_prob = 1.0
                        !
                        !                        end if

                        !                        apath(t+1) = apolOPTfull(zpath(t),wpath(t),apath(t),spath_ind(t))
                        !                        qpath(t)= max(qOPTfull(zpath(t),wpath(t),apath(t),spath_ind(t),apath(t+1)),eps)

                        if( shockone .LE. rep_prob ) then
                            ! In this case, the sovereign has repaid its debts

                            if(t > 1) then

                                grsim(t) = zm(zpath(t))
                                hat_ysim(t)=zm(zpath(t))*hat_ysim(t-1)*(zbase)

                                ysim(t)=hat_ysim(t)*(zbase + spath(t))
                                grsim_real(t) = log(ysim(t)/ysim(t-1))

                            else

                                grsim(t) = zm(zpath(t))
                                hat_ysim(t) = zbase

                                ysim(t)= hat_ysim(t)*(zbase + spath(t))
                                grsim_real(t) = log( zm(zpath(t))*(zbase + spath(t)) )

                            end if
                            csim(t) = hat_ysim(t)*( spath(t) + zamvec(zpath(t),apath(t)) - qpath(t)*(am(apath(t+1))-(1.0-lambda)*am(apath(t))/zm(zpath(t))) )

                            defpath(t)=0
                            history(t+1)=1
                            history(t)=1
                            nx(t) = ysim(t) - csim(t)

                            if (amdiff(zpath(t),apath(t), apath(t+1))>0) then  ! if there is buyback of debt

                                debts(t)=amser(zpath(t),apath(t))+qpath(t)*amdiff(zpath(t),apath(t), apath(t+1)) ! debt-service

                            else

                                debts(t)=amser(zpath(t),apath(t))

                            end if

                            first_pd_def(t+1) = 0
                            if (t==1) then

                                inc(t)=1
                                inc2(t)=1

                            else

                                inc(t)=inc(t-1)+1
                                inc2(t)=inc2(t-1)+1

                            end if

                        else

                            ! In this case, the sovereign has defaulted on its debts (first period of default)
                            defpr = .TRUE.

                        end if

                    elseif(crisis_path(t) .EQ. 2) then

                        ! HERE, WE HAVE PESSIMISTIC BELIEFS

                        if(apolPESfull(zpath(t),wpath(t),apath(t),spath_ind(t)) .EQ. A) then

                            !  Here, creditors offer zero price even for zero issuance (necessitates default today)
                            rep_prob = 0.0
                            apath(t+1) = A
                            qpath(t) = qbase

                        else

                            ! Ratio of prices gives repayment probability
                            rep_prob = qPESfull(zpath(t),wpath(t),apath(t),spath_ind(t),apolPESfull(zpath(t),wpath(t),apath(t),spath_ind(t)))/q(zpath(t),wpath(t),apolPESfull(zpath(t),wpath(t),apath(t),spath_ind(t)))
                            apath(t+1) = apolPESfull(zpath(t),wpath(t),apath(t),spath_ind(t))
                            qpath(t)= max(qPESfull(zpath(t),wpath(t),apath(t),spath_ind(t),apath(t+1)),eps)

                        end if

                        !                        ! Eliminate interior defaults (if we're rigging the game)
                        !                        if( (rep_prob > 1.0e-6) .AND. (rep_prob < 1.0- 1.0e-6) .AND. load_the_die) then
                        !
                        !                            rep_prob = 1.0
                        !
                        !                        end if

                        !                        apath(t+1) = apolPESfull(zpath(t),wpath(t),apath(t),spath_ind(t))
                        !                        qpath(t)= max(qPESfull(zpath(t),wpath(t),apath(t),spath_ind(t),apath(t+1)),eps)

                        if( shockone .LE. rep_prob ) then
                            ! In this case, the sovereign has repaid its debts

                            if(t > 1) then

                                grsim(t) = zm(zpath(t))
                                hat_ysim(t)=zm(zpath(t))*hat_ysim(t-1)*(zbase)

                                ysim(t)=hat_ysim(t)*(zbase + spath(t))
                                grsim_real(t) = log(ysim(t)/ysim(t-1))

                            else

                                grsim(t) = zm(zpath(t))
                                hat_ysim(t) = zbase

                                ysim(t)= hat_ysim(t)*(zbase + spath(t))
                                grsim_real(t) = log( zm(zpath(t))*(zbase + spath(t)) )

                            end if
                            csim(t) = hat_ysim(t)*( spath(t) + zamvec(zpath(t),apath(t)) - qpath(t)*(am(apath(t+1))-(1.0-lambda)*am(apath(t))/zm(zpath(t))) )

                            defpath(t)=0
                            history(t+1)=1
                            history(t)=1
                            nx(t) = ysim(t) - csim(t)

                            if (amdiff(zpath(t),apath(t), apath(t+1))>0) then  ! if there is buyback of debt

                                debts(t)=amser(zpath(t),apath(t))+qpath(t)*amdiff(zpath(t),apath(t), apath(t+1)) ! debt-service

                            else

                                debts(t)=amser(zpath(t),apath(t))

                            end if

                            first_pd_def(t+1) = 0
                            if (t==1) then

                                inc(t)=1
                                inc2(t)=1

                            else

                                inc(t)=inc(t-1)+1
                                inc2(t)=inc2(t-1)+1

                            end if

                        else

                            ! In this case, the sovereign has defaulted on its debts (first period of default)
                            defpr = .TRUE.

                        end if

                    else

                        ! HERE, WE HAVE CONCERNED BELIEFS

                        if(apolCONfull(zpath(t),wpath(t),apath(t),spath_ind(t)) .EQ. A) then

                            !  Here, creditors offer zero price even for zero issuance (necessitates default today)
                            rep_prob = 0.0
                            apath(t+1) = A
                            qpath(t) = qbase

                        else

                            ! Ratio of prices gives repayment probability
                            rep_prob = qCONfull(zpath(t),wpath(t),apath(t),spath_ind(t),apolCONfull(zpath(t),wpath(t),apath(t),spath_ind(t)))/q(zpath(t),wpath(t),apolCONfull(zpath(t),wpath(t),apath(t),spath_ind(t)))
                            apath(t+1) = apolCONfull(zpath(t),wpath(t),apath(t),spath_ind(t))
                            qpath(t)= max(qCONfull(zpath(t),wpath(t),apath(t),spath_ind(t),apath(t+1)),eps)

                        end if

                        ! Eliminate only borderline interior defaults (if we're rigging the game)
                        if( (rep_prob > 1.0e-6) .AND. (rep_prob < 1.0- 1.0e-6) .AND. load_the_die .AND. ( ( qCONfull(zpath(t),wpath(t),apath(t),spath_ind(t),apath(t+1)-1) < 1.0e-6  ) .OR. ( qCONfull(zpath(t),wpath(t),apath(t),spath_ind(t),apath(t+1)-1) > q(zpath(t),wpath(t),apath(t+1)-1)-1.0e-6  ) ) ) then

                            rep_prob = 1.0

                        end if

                        !                        ! Eliminate interior defaults (if we're rigging the game)
                        !                        if( (rep_prob > 1.0e-6) .AND. (rep_prob < 1.0- 1.0e-6) .AND. load_the_die) then
                        !
                        !                            rep_prob = 1.0
                        !
                        !                        end if

                        !                        apath(t+1) = apolCONfull(zpath(t),wpath(t),apath(t),spath_ind(t))
                        !                        qpath(t)= max(qCONfull(zpath(t),wpath(t),apath(t),spath_ind(t),apath(t+1)),eps)
                        !
                        if( shockone .LE. rep_prob ) then
                            ! In this case, the sovereign has repaid its debts

                            if(t > 1) then

                                grsim(t) = zm(zpath(t))
                                hat_ysim(t)=zm(zpath(t))*hat_ysim(t-1)*(zbase)

                                ysim(t)=hat_ysim(t)*(zbase + spath(t))
                                grsim_real(t) = log(ysim(t)/ysim(t-1))

                            else

                                grsim(t) = zm(zpath(t))
                                hat_ysim(t) = zbase

                                ysim(t)= hat_ysim(t)*(zbase + spath(t))
                                grsim_real(t) = log( zm(zpath(t))*(zbase + spath(t)) )

                            end if
                            csim(t) = hat_ysim(t)*( spath(t) + zamvec(zpath(t),apath(t)) - qpath(t)*(am(apath(t+1))-(1.0-lambda)*am(apath(t))/zm(zpath(t))) )

                            defpath(t)=0
                            history(t+1)=1
                            history(t)=1
                            nx(t) = ysim(t) - csim(t)

                            if (amdiff(zpath(t),apath(t), apath(t+1))>0) then  ! if there is buyback of debt

                                debts(t)=amser(zpath(t),apath(t))+qpath(t)*amdiff(zpath(t),apath(t), apath(t+1)) ! debt-service

                            else

                                debts(t)=amser(zpath(t),apath(t))

                            end if

                            first_pd_def(t+1) = 0
                            if (t==1) then

                                inc(t)=1
                                inc2(t)=1

                            else

                                inc(t)=inc(t-1)+1
                                inc2(t)=inc2(t-1)+1

                            end if

                        else

                            ! In this case, the sovereign has defaulted on its debts (first period of default)
                            defpr = .TRUE.

                        end if

                    end if

                    ! In this case, sovereign has defaulted. Beliefs don't matter now after price and issuance known and set (qpath(t), apath(t+1))
                    if(defpr) then

                        def_total = def_total+1

                        !!!!! inc and inc2 keep track of periods of exclusion for purposes of computing moments: Exclude these when
                        !!!!! and 20 quarters after re-entry when computing key moments. inc and inc2 differ by only 1, generally.
                        !!!!! inc is used to take average of spreads (excludes period of default) and inc2 used to compute default
                        !!!!! frequency (which includes period of default)
                        first_pd_def(t+1) = 1

                        if (t==1) then

                            inc2(t)=1

                        else

                            inc2(t)=inc2(t-1)+1

                        end if

                        ! Don't pay default costs in first period:
                        if(t > 1) then

                            grsim(t) = zm(zpath(t))
                            hat_ysim(t)=zm(zpath(t))*hat_ysim(t-1)*(zbase)

                            ysim(t)=hat_ysim(t)*(spath(t) + zbase) ! average s-shock in first pd of default
                            grsim_real(t) = log(ysim(t)/ysim(t-1))

                        else

                            grsim(t) = zm(zpath(t))
                            hat_ysim(t)= zbase

                            ysim(t)= hat_ysim(t)*( spath(t) + zbase )
                            grsim_real(t) = log( zm(zpath(t))*(spath(t) + zbase ) )

                        end if
                        csim(t)= hat_ysim(t)*( spath(t) + zbase - qpath(t)*(am(apath(t+1))-(1.0-lambda)*am(apath(t))/zm(zpath(t))) ) ! Get to run away with auction revenue
                        apath(t+1)= A  ! reset to zero debt level

                        if(qpath(t) < .01) then
                            qpath(t) = qbase ! reset prices if default is known prior to realization of epsilon to keep moments finite
                        end if


                        defpath(t)=1
                        history(t+1)=0
                        history(t)=1
                        debts(t)=0 ! Did not service debt here
                        nx(t)=ysim(t) - csim(t)

                    end if


                else ! Case of exclusion from credit markets
                    first_pd_def(t+1) = 0
                    inc2(t)=0
                    inc(t)=0
                    if(t > 1) then

                        grsim(t) = zm(zpath(t))
                        hat_ysim(t)=zm(zpath(t))*hat_ysim(t-1)*(zbase)

                        ysim(t)=hat_ysim(t)*(defz(zpath(t))+ spath(t))
                        grsim_real(t) = log(ysim(t)/ysim(t-1))

                    else

                        grsim(t) = zm(zpath(t))
                        hat_ysim(t)= zbase

                        ysim(t)= hat_ysim(t)*(defz(zpath(t))+ spath(t))
                        grsim_real(t) = log( zm(zpath(t))*(defz(zpath(t))+ spath(t)) )

                    end if

                    csim(t)=ysim(t)
                    defpath(t)=0
                    history(t+1)=0
                    apath(t+1)=A
                    qpath(t)=qbase
                    debts(t)=0
                    nx(t)=0
                end if
            end do

            ! Compute annual yields from quarterly prices
            rpath=(lambda+coup-lambda*qpath(1:TB-1))/qpath(1:TB-1)  !interest rate

            spreadpath=rpath-rbase


            spreadpatha=(1+rpath)**4-(1+rbase)**4 ! annual spreads


            assets=am(apath)/(zm(zpath)*(1+spath(t))) ! Debt-to-GDP = b/y = am/(zm*(1+sm))
            cred_wealth=wm(wpath) ! Foreign wealth

            dummy(:)=minloc(apath(1:TB))  ! finding minimum asset choice to see if it is binding
            if (ss==1) then

                aminloc=apath(dummy(1))

            elseif (aminloc>apath(dummy(1))) then

                aminloc=apath(dummy(1))

            end if

            num=0
            num2=0
            do i=1,TB-1
                if (inc(i)>maxinc) then

                    num=num+1
                    samp1(i)=.TRUE.

                else

                    samp1(i)=.FALSE.

                end if

                if (inc2(i)>maxinc) then

                    num2=num2+1

                end if
            end do

            spreadt=spreadpatha(TB-nos:TB-1)

            if((num > 0.0) .AND. (num2 > 0.0)) then

                meanspread(ss)=sum(spreadpatha, MASK=samp1)/num
                meandebt(ss) = sum(-assets(2:TB), MASK=samp1)/num !sum(am(apath(2:TB)), MASK=samp1)/num ! new debt/gdp = b'/y = am

                debt_burden = PACK( -assets(2:TB)*(lambda+coup) ,MASK=samp1 )
                temp_size = INT(size(debt_burden)*doub_percentile)

                meandebtserv(ss)= quantile( temp_size, debt_burden )

                num3=sum(defpath(1:TB-1),MASK=inc2(1:TB-1)>maxinc)/num2 ! quarterly default probability
                defaultpc(ss)=one-(one-num3)**4.0  ! yearly default probability
                call moment(spreadt, stdall1(3,ss))

            else
                meanspread(ss) = 1.0/eps
                meandebt(ss) = 0.0
                meandebtserv(ss) = 0.0
                defaultpc(ss) = 1.0
                stdall1(3,ss) = 1.0/eps

            end if


            !turning quarterly debt service data to yearly
            t=0
            j=0
            do i=1,TB-1
                j=j+1
                ds1(j)=debts(i)
                ds2(j)=ysim(i)
                if (samp1(i)) then
                dsint(j)=1
                else
                    dsint(j)=0
                end if
                if (j==4) then
                    j=0
                    if (sum(dsint)>0) then
                        t=t+1
                        debtssim(t)=dot_product(dsint,ds1)/dot_product(dsint,ds2)
                    end if
                end if
            end do


            logy=log(ysim(TB-nos:TB-1))
            nxy = nx(TB-nos:TB-1)/ysim(TB-nos:TB-1)
            logc = log(csim(TB-nos:TB-1))
            cgrowth(1) = muz
            cgrowth(2:TB-1) = logc(TB-nos+1:TB-1)-logc(TB-nos:TB-2)

            call HPFILT(logy,ydevs,vscrap,nos,real(1600,doub),1)
            call HPFILT(logc,cdevs,vscrap,nos,real(1600,doub),1)
            call HPFILT(nxy,nxy_devs,vscrap,nos,real(1600,doub),1)

            call moment(ydevs, stdall1(1,ss))
            call moment(cdevs, stdall1(2,ss))
            call moment(nxy_devs, stdall1(4,ss))

            ! Write relevant entries into the simulation text file (only one processor needs to do this)
            if( id .EQ. 0) then

                do int1 = 1,TB-1
                    write(30,'(f8.5,f8.5,2X,f8.5,f8.5,f8.5,i6.3,2X,f8.5,f8.5,2X,i8.5,i8.5,f8.5)') &
                    ydevs(int1), -assets(int1), nxy_devs(int1), cgrowth(int1), spreadpatha(int1), crisis_path(int1),  &
                    grsim(int1), spath(int1), zpath(int1), apath(int1), grsim_real(int1)
                end do

            end if

            spread_devs = spreadpatha(TB-nos:TB-1)

            corrall1(1,ss)=corr(ydevs, cdevs)
            corrall1(2,ss)=corr(ydevs, spread_devs)
            corrall1(3,ss)=corr(ydevs, nxy_devs)
            corrall1(4,ss)=corr(cdevs, spread_devs)
            corrall1(5,ss)=corr(cdevs, nxy_devs)
            corrall1(6,ss)=corr(spread_devs, nxy_devs)

        end do
        close(30)

        corr1=sum(corrall1,DIM=2)/real(sim)
        stdn1=sum(stdall1,DIM=2)/real(sim)

        meanmeandef=sum(defaultpc)/real(sim)
        meanmeanspr=sum(meanspread)/real(sim)
        meanmeandebt= sum(meandebt)/real(sim)
        meanmeansv=sum(meandebtserv)/real(sim)
        meandebtsvol=sum(debtsvol)/real(sim)


        funcreg(kkk,1)= meanmeandebt !welfc
        funcreg(kkk,2)=meanmeanspr
        funcreg(kkk,3)=meanmeansv
        funcreg(kkk,4)=meanmeandef
        funcreg(kkk,5) = stdn1(3)
        funcreg(kkk,6) = stdn1(2)/stdn1(1)
        funcreg(kkk,7) = stdn1(4)/stdn1(1)
        funcreg(kkk,8) = corr1(1)
        funcreg(kkk,9) = corr1(3)
        funcreg(kkk,10) = corr1(2)

    end subroutine model_simul


    real(doub) function FPfun(curr_eps,ia_curr, iz_curr, ia_fut)
        implicit none

        integer, intent(in) :: ia_curr, ia_fut, iz_curr
        real(doub), intent(in)::curr_eps

        FPfun = cdf_epsilon(curr_eps)/(1.0-(1.0-cdf_epsilon(curr_eps))*lend_share*max(-am(ia_fut)+(1.0-lambda)*am(ia_curr)/zm(iz_curr),0.0)/( -(lambda+coup)*am(ia_curr)/zm(iz_curr)-am(ia_fut) )  )

        return
    end function FPfun


    subroutine interval_bisection_eps(epsLin,epsHin)
        implicit none

        real(doub), intent(in) :: epsLin, epsHin
        real(doub):: epsL, epsH, HL, HH, conv_tol_bis

        conv_tol_bis = 1.0e-12

        epsL = epsLin
        epsH = epsHin


        HL = sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(epsL,ia,iz,ia2)*paym(ia2))+Vn_fut(ia2,iz,iw) - &
        ( sov_utility( sm(isind) + zbase + (1.0-lend_share)*FPfun(epsL,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) + epsL)

        HH = sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(epsH,ia,iz,ia2)*paym(ia2))+Vn_fut(ia2,iz,iw) - &
        ( sov_utility( sm(isind) + zbase + (1.0-lend_share)*FPfun(epsH,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) + epsH)



        if( HH*HL > 0.0) then

            open(50,file='big_problem.txt') !20 is the file number (can be anything)
            write(50,*) iz, iw, ia, isind
            write(50,'(i4)') ia2
            write(50,*) paym(ia2), HL, HH
            write(50,*) FPfun(epsL,ia,iz,ia2), FPfun(epsH,ia,iz,ia2), epsL, epsH
            close(50)

        else

            if( HL > 0.0 ) then ! In this case, the H-function is decreasing; otherwise it's increasing


                do while(epsH-epsL > conv_tol_bis)

                    eps_star = (epsL+epsH)/2.0

                    eps_gap = sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2)) + Vn_fut(ia2,iz,iw) - &
                    ( sov_utility( sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) + eps_star)

                    if(eps_gap .GE. 0) then

                        epsL = eps_star

                    else

                        epsH = eps_star

                    end if

                end do



            else

                do while(epsH-epsL > conv_tol_bis)

                    eps_star = (epsL+epsH)/2.0

                    eps_gap = sov_utility(sm(isind) + zamvec(iz,ia) + FPfun(eps_star,ia,iz,ia2)*paym(ia2)) + Vn_fut(ia2,iz,iw) - &
                    ( sov_utility( sm(isind) + zbase + (1.0-lend_share)*FPfun(eps_star,ia,iz,ia2)*paym(ia2) ) + Vd_fut(izw) + eps_star)

                    if(eps_gap .GE. 0) then

                        epsH = eps_star

                    else

                        epsL = eps_star

                    end if

                end do

            end if

        end if


        eps_star = (epsH+epsL)/2.0

    end subroutine interval_bisection_eps


    ! Subroutines moment and corr simply compute volatilities and correlation coefficients
    ! for vector inputs

    subroutine moment(xx, std)
        implicit none
        integer::mm
        real(doub),dimension(:), intent(in)::xx
        real(doub), intent(out):: std
        real(doub):: mean

        mm=sum(history(TB-nos:TB-1),MASK=samp1(TB-nos:TB-1))
        mean=sum(xx,MASK=samp1(TB-nos:TB-1))/mm
        std=sum((xx-mean)**2,MASK=samp1(TB-nos:TB-1))/mm
        std=sqrt(std)

        return
    end subroutine moment

    real(doub) function corr(xx1,xx2)
        implicit none
        integer::mm
        real(doub)::mean1, mean2, std1, std2
        real(doub), dimension(:), intent(in)::xx1, xx2

        mm=sum(history(TB-nos:TB-1),MASK=samp1(TB-nos:TB-1))
        mean1=sum(xx1,MASK=samp1(TB-nos:TB-1))/mm
        mean2=sum(xx2,MASK=samp1(TB-nos:TB-1))/mm
        std1=sqrt(sum((xx1-mean1)**2,MASK=samp1(TB-nos:TB-1))/mm)
        std2=sqrt(sum((xx2-mean2)**2,MASK=samp1(TB-nos:TB-1))/mm)
        corr=sum((xx1-mean1)*(xx2-mean2),MASK=samp1(TB-nos:TB-1))/mm
        corr=corr/(std1*std2)

        return
    end function corr

    ! Prescott's original HP-filter routine: For use in simulations
    ! ----------------------------------------------------------------------
    !  SR: hpfilt
    !  Kalman smoothing routine for HP filter written by E Prescott.
    !   y=data series, d=deviations from trend, t=trend, n=no. obs,
    !   s=smoothing parameter (eg, 1600 for std HP).
    !   Array v is scratch area and must have dimension at least 3n.
    !   If IOPT=1 and n and s are the same as for the previous call,
    !   the numbers in v are not recomputed.  This reduces execution
    !   time by about 30 percent.  Note that if this option is exercised,
    !   v cannot be used for other purposes between calls.
    !   This version does NOT release the trend in order to save memory.
    ! ----------------------------------------------------------------------
    SUBROUTINE HPFILT(Y,D,V,N,S,IOPT)
        INTEGER*4 N
        REAL(doub) Y(N),T(N),V(N,3),D(N),S
        REAL(doub) M1,M2,V11,V12,V22,X,Z,B11,B12,B22,DET,E1,E2,SS
        INTEGER*4 IOPT,NN,I,I1,IB
        DATA SS,NN/0.D0,0/
        !
        !     compute sequences of covariance matrix for f[x(t),x(t-1) | y(<t)]
        !
        IF(IOPT.NE.1.OR.NN.NE.N.OR.S.NE.SS)  THEN
            SS=S
            NN=N
            V11=1.D0
            V22=1.D0
            V12=0.D0
            DO 5 I=3,N
            X=V11
            Z=V12
            V11=1.D0/S + 4.D0*(X-Z) + V22
            V12=2.D0*X - Z
            V22=X
            DET=V11*V22-V12*V12
            V(I,1)=V22/DET
            V(I,3)=V11/DET
            V(I,2)=-V12/DET
            X=V11+1.D0
            Z=V11
            V11=V11-V11*V11/X
            V22=V22-V12*V12/X
            V12=V12-Z*V12/X
            5    CONTINUE
        ENDIF
        !
        !     this is the forward pass
        !
        M1=Y(2)
        M2=Y(1)
        DO 10 I=3,N
        X=M1
        M1=2.0*M1-M2
        M2=X
        T(I-1)= V(I,1)*M1+V(I,2)*M2
        D(I-1)= V(I,2)*M1+V(I,3)*M2
        DET=V(I,1)*V(I,3)-V(I,2)*V(I,2)
        V11=V(I,3)/DET
        V12=-V(I,2)/DET
        Z=(Y(I)-M1)/(V11+1.D0)
        M1=M1+V11*Z
        M2=M2+V12*Z
        10     CONTINUE
        T(N)=M1
        T(N-1)=M2
        !
        !       this is the backward pass
        !
        M1=Y(N-1)
        M2=Y(N)
        DO 15 I=N-2,1,-1
        I1=I+1
        IB=N-I+1
        X=M1
        M1=2.D0*M1 - M2
        M2=X
        !
        !           combine info for y(.lt.i) with info for y(.ge.i)
        !
        IF(I.GT.2)                 THEN
            E1=V(IB,3)*M2 + V(IB,2)*M1 + T(I)
            E2=V(IB,2)*M2 + V(IB,1)*M1 + D(I)
            B11=V(IB,3)+V(I1,1)
            B12=V(IB,2)+V(I1,2)
            B22=V(IB,1)+V(I1,3)
            DET=B11*B22-B12*B12
            T(I)=(-B12*E1+B11*E2)/DET
        ENDIF
        !
        !           end of combining
        !
        DET=V(IB,1)*V(IB,3)-V(IB,2)*V(IB,2)
        V11=V(IB,3)/DET
        V12=-V(IB,2)/DET
        Z=(Y(I)-M1)/(V11+1.D0)
        M1=M1+V11*Z
        M2=M2+V12*Z
        15     CONTINUE
        T(1)=M1
        T(2)=M2
        DO I=1,N
              D(I)=Y(I)-T(I)
        END DO
        RETURN
    END SUBROUTINE HPFILT
    !***********************************************************************

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(doub) function cdf_epsilon(sample_epsilon)
        implicit none
        real(doub)::sample_epsilon

        if(sample_epsilon > epsilon_upper_bar) then

            cdf_epsilon = 1.0

        elseif(sample_epsilon < epsilon_lower_bar) then

            cdf_epsilon = 0.0

        else

            cdf_epsilon = (sample_epsilon-epsilon_lower_bar)/(epsilon_upper_bar-epsilon_lower_bar)

        end if

        return
    end function cdf_epsilon

    real(doub) function sov_utility(c)
        implicit none
        real(doub):: c

        if(c > 0) then

            sov_utility = c**(1.0-gamma)/(1.0-gamma)

        else

            sov_utility = -1e10

        end if

    end function sov_utility

        ! Efficient quantile finder for simulations. Hoare's algorithm, e.g., median = quantile( size(a)/2, a)

        recursive function quantile( k, a) result( value)
          integer,             intent (in)  :: k          ! position in array
          real(8), dimension (:), intent (in)  :: a
          real(8)                              :: value      ! output value of quantile
          integer                           :: j
          real(8)                              :: ak

          ak = a( k)
          j = count( a < ak)  ! how many a(:) < ak

          if( j >= k) then
              value = quantile( k, pack( a, a < ak))
          else
      	    j = count( a > ak) + k - size( a)
      	    if( j > 0) then
      		    value = quantile( j, pack( a, a > ak))
      	    else
      		    value = ak
            end if
          end if

          end function quantile

end module funcall
