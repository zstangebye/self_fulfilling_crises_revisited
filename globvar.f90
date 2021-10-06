
module globvar
implicit none
! IF RUNNING IN PARALLEL, uncomment " include 'mpif.h' ", set mpi_on = .TRUE.,
! and comment the definition of the four MPI variables (which are also defined in mpif.h)

 include 'mpif.h'
!use mpi

logical, parameter :: mpi_on = .TRUE.
!integer :: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, MPI_LOGICAL, MPI_INTEGER ! These variables are not used if mpi is not on but have to be defined, nevertheless

! Throughout the code, 'zm' typically refers to exp(g), where g is the economy's growth rate

integer, parameter::doub=8, sing=4
! Model parameters (relevant)
real(doub), parameter:: sdz = .0138298 ! (MX) ! .0082939 ! (IT) !         GROWTH RATE VOLATILITY (unconditional)
real(doub), parameter:: sds = .0113592 ! (MX) !  .0096982 !(IT) ! .07 !     s-SHOCK VOLATILITY
real(doub), parameter:: rhoz = 0.3308571 ! (MX) !  .579109 ! (IT) !     GROWTH RATE PERSISTENCE
real(doub), parameter:: muz= .0033107 ! (MX) ! .0022693 ! (IT) !     AVERAGE GROWTH RATE
real(doub), parameter:: lambdaset = 1.0 ! .0638 ! (MXLT) ! .0286 !(ITLT) !  1.0 !    1.0 !0.125 ! MATURITY OF DEBT (FRACTION THAT MATURES)
real(doub), parameter:: rbase=0.01 ! RISK-FREE RATE
real(doub), parameter:: ent=0.125 ! RE-ENTRY PROBABILITY SHOCK
real(doub), parameter:: lend_share = 0.0 ! 0.0 ! Share of auction revenue returned to lenders pro-rata in default
real(doub), parameter:: tax_rev_share = 0.176 ! (MX) ! 0.29 ! (IT) !    Tax revenue share of GDP (set 1.0 if treating output as GDP; .106 if output is tax rev)
real(doub), parameter:: doub_percentile = .97 ! (MX) ! .93 ! (IT) ! Used to match fraction of crisis periods in calibration

real(doub), parameter :: epsilon_upper_bar = .0001 ! 0.0001 ! intratemporal default shock (non-degenerate)
real(doub), parameter :: epsilon_lower_bar = .0000 ! 0.0  ! uniform distribution necessary with zero mean


real(doub) :: aminset =  -.75 ! (MXST) ! -8.0 ! (MXLTEG) ! -15.0 ! (IT)  ! -1.0 ! (ITST) !  -5.0 ! (MXSTEG) ! -8.0 ! (MX) !  Bound on asset grid


! COUPON PAYMENTS (FRACTION OF DEBT)
real(doub), parameter :: coupset = .03

logical, parameter :: GMM_on = .FALSE. ! SWITCH THAT DETERMINES WHETHER RUNNING ONCE OR MATCHING MOMENTS

! If computing once, these parameters are fixed;
! During GMM, these parameters are the initial guess
real(doub) :: beta_true = 0.8 ! (MX) ! 0.95 ! (IT) !     SOVEREIGN DISCOUNT FACTOR
real(doub) :: defp0 =  0.176 ! (MXST) ! 0.467 !(MXLT) !  .665 ! (ITLT) !  .143 ! (ITST) ! 0.455 ! (MXLTEG) !     PROPORTIONAL OUTPUT COST OF DEFAULT
real(doub), parameter :: pes_bel = 0.0, con_bel = 0.0079 ! 0.0025 !  .0044 ! 0.003 ! (ITLT) ! 0.0044 !(MXLT) ! !! Belief probabilities (must sum to one)
real(doub), parameter :: opt_bel = 1.0-pes_bel-con_bel

real(doub) :: sdw = 0.0001 ! 2.640 !.653 ! LENDER WEALTH VOLATILITY (unconditional)


real(doub) :: defp1 = 0.0 ! CURVATURE PARAMETER (IF ZERO, NO CURVATURE)
character(LEN=100) :: tagg = 'revision2_mex_ST.txt', temp_str, temp_name

! output shock z under default
! Comes from non-normalized: C = Y*(1-defp0) +defp1*Y(-1)
! Divide by Y to derive:     c = 1-defp0+defp1/exp(g)
! i.e. Higher g implies relatively lower consumption post-default

! Set integer parameters used as grid sizes
integer, parameter::S=11, Z=15, A=500, W = 1, EPS_CHECK = 100 ! EPS_CHECK = 50  ! GRID SIZES 11 25 350 1
character(LEN=10), parameter :: formatter = '(500f9.6)' ! Change this when changing size of 'A'
character(LEN=9), parameter :: formatter2 = '(11f12.7)'  ! Change this when changing size of 'S'
character(LEN=9), parameter :: formatter3 = '(500i6.3)'  ! Change this when changing size of 'A'

! COMPUTATIONAL TOOLS (SMOOTHING PARAMETERS)
real(doub) :: qsmooth = 0.5, qsmooth2 = 0.5 
real(doub) :: vsmooth = 0.5
real(doub) :: converge_tol =  0.00001 ! (ST case) !


! Parameters used for robustness/sensitivity (IRRELEVANT IN BENCHMARK MODEL)
real(doub), parameter:: rhow = 0.91 ! LENDER WEALTLH PERSISTENCE 
real(doub), parameter:: sigma_zw = 0.0 ! CORRELATION (g_z,w)
real(doub) :: u_ben_def = 0.0
real(doub) :: smean = 0.0, sdef
logical, parameter :: mixing_fix = .FALSE., show_obj = .TRUE., load_the_die = .FALSE.
logical :: stop_mix_fix1, stop_mix_fix2
real(doub) :: w_base =  3.75  ! !!!!!!!!!! WEALTH OF LENDERS (MULTIPLE OF ENDOWMENT)
! No-deals case is zero
real(doub), parameter:: deal_frac_debt = 0.0 !0.9*(1.0-lambdaset) !(2.0-lambdaset)/2.0 ! Fraction of debt issued/bought back during desperate deal
real(doub) :: ppismooth = 1.0  ! no deals case
!real(doub) :: ppismooth =  0.5 ! benchmark case

logical, parameter :: no_deals = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TECHNICAL DEFINITIONS
real(doub) :: eps_star, eps_gap 
real(doub)::coup = coupset ! technical re-defining of coupon (equal to risk-free rate)
real(doub)::lambda = lambdaset ! technical re-defining of maturity
real(doub), parameter:: muw = 0.0 ! Keep at zero (normalization)
integer, parameter::gamma=2 ! RISK-AVERSION COEFFICIENT (CANNOT BE CHANGED)
integer :: def_total, def_RO, num_eq
real(doub) :: RO_fraction, eps_check_num

real(doub), parameter :: GMM_coarseness = 2.0 !100
integer :: ticker
real(doub), parameter:: mxspt=0.0

! Used in computation of multivariate normal CDF (MVNDST)
real(doub), parameter :: ABSEPS = 0.00005, RELEPS = 0.0
integer, parameter::  DN = 2
integer :: INFIN(DN), maxT
integer :: MXPNTS = 2*1000, INFORM = 0
real(doub) :: ERR, int_val, correl(2)

integer, parameter :: NZW = Z*W
integer, parameter:: NZA=Z*A
integer, parameter:: NZAW = Z*A*W
integer, parameter:: NZAWS = Z*A*W*S
integer, parameter:: NZAWSA = Z*A*W*S*A
integer, parameter:: NZAA=Z*A*A
integer, parameter :: NZZA = Z*Z*A
integer, parameter :: NZZWWA = Z*Z*W*W*A
integer, parameter:: DA=A

real(doub), dimension(NZW) :: beta_aug
real(doub), parameter :: zbase = 1.0 ! tax_rev_share ! 1.0
integer, parameter :: tmin = 10 ! Typically 250

! Used to discretize the grids
real(doub), parameter :: amax = 0.0, amax_grid_break = -.3
integer, parameter :: course_grid_A = 10
real(doub):: amin, epsilon_delta, eps_flag, eps_flag2
real(doub):: dela1, dela2
real(doub), dimension(2) :: lb, ub
real(doub), dimension(A):: am
integer, dimension(A) :: ambound
real(doub), dimension(S) :: sm, smdiff
real(doub), dimension(S+1) :: smthresh ! Used in expectation
real(doub), parameter :: thrmin=-10, thrmax=10 ! Used in simulation

! Used with Tauchen's method to turn continuous Markov processes into discrete ones
! Some used in computation of normal cdf (CDFNOR)
real(doub), parameter::m=3.0, m2=3.0
real(doub):: lrsdz, delz, dels, mu_temp, lrsdw, delw, mu_temp_1pd_RA
real(doub), dimension(NZW,NZW)::pdfzw, cdfzw, dividzw
real(doub),dimension(Z)::zm, logzm
real(doub), dimension(W)::wm, logwm
real(doub), dimension(S)::pdfs, divids
integer,parameter::xx=1
real(doub)::x1,x2, non, boun, x01, x02, auc_rev
integer:: sta

! Used in the value function iteration
integer, parameter:: time= 800 !800 !10000
real(doub)::qbase, current_eps, current_eps2, rep_prob
real(doub), dimension(A) :: const_vec1, const_vec2, const_vec3, const_vec4, const_vec5
integer:: srcmax
logical :: eps_check_not_done
real(doub), dimension(EPS_CHECK):: epsilon_grid
real(doub), dimension(NZAW)::Vn, Vn2 ! Value function; s is not a state variable here; we take an expectation over it in each state to get the value
real(doub), dimension(NZAWS) :: apolVEC
real(doub), dimension(NZAWSA) :: qlongVEC
real(doub), dimension(:), allocatable::qsum, mu_lend, qsum2, mu_lend2, rep_bin_full, qsum_1pd_RA, mu_lend_1pd_RA, debt_burden
real(doub),dimension(A)::shorts
real(doub):: uexpt, qpr, q_num, qdenom, q_short, q_num_1pd_RA, qdenom_1pd_RA
real(doub), dimension(S,NZW):: def_fix
real(doub), dimension(Z):: def_fixn
real(doub), dimension(S,NZW):: Vds
real(doub), dimension(NZW):: vd_fut, Vd
real(doub), dimension(Z,A):: zamvec
real(doub), dimension(A,Z,W) :: Vn_fut
real(doub), dimension(:,:,:), allocatable :: q, q2, q_RN, q_1pd, q_RN_LT
real(doub), dimension(S,A) :: qOPT, qPES, qCON, epsOPT, epsPES, epsCON
integer, dimension(S) :: apolOPT, apolPES, apolCON
integer, dimension(Z,W,A,S) :: apolOPTfull, apolPESfull, apolCONfull
real, dimension(Z,W,A,S,A) :: qOPTfull, qPESfull, qCONfull
logical, dimension(Z,W,A) :: crisis_zone
real(doub), dimension(Z,W,A) :: crisis_thresh
real(doub),dimension(A)::fixnum, util_vecOPT, util_vecPES, util_vecCON, apolTEMP, con_prob_vec, con_repay_val_vec, con_def_val_vec, con_eps_shock_cont
real(doub),dimension(nzw)::zvec, zvec2
real(doub),dimension(S)::ucons, mu_lend_cons, mu_lend_cons2, mu_lend_cons3, mu_lend_cons_1pd_RA, mu_lend_cons_1pd_RA2
real(doub):: latval
real(doub), dimension(NZW) :: mu_exp, rep_bin, qpr_1pd_RA, mu_exp_1pd_RA
real(doub), dimension(S)::q_fut, q_fut_1pd_RA, q_fut_1pd_RA2
logical, dimension(A)::rela
real(doub), dimension(Z,W):: Vdsnone, V_ro ! Value of defaulting and not being in the credit market (for sure) next period; no asset holdings and does not depend on w
real(doub), dimension(A):: defthr
real(doub), dimension(Z):: defz

real(doub):: errq, errv, errPpi
real(doub), dimension(time)::errtot_v, errtot_q, errtot_ppi
real(doub),dimension(:),allocatable:: val1, val2, val3, val4, val5, val1agg, val2agg, val3agg, val4agg, val5agg
integer,dimension(:),allocatable:: val_int1, val_int2, val_int3, val_int1agg, val_int2agg, val_int3agg

real(doub), dimension(A)::paym, cons1
integer::intnc
real(doub):: defthnc, defthc

! Temporary policy functions
integer:: combint
real(doub), dimension(DA)::sthcomb, thcomb
integer, dimension(DA)::chass, asscomb
logical, dimension(DA)::chdef
real(doub), dimension(DA)::thcomb2 ! The '2' attached to the policy function represents the confidence crisis sunspot equilibrium
integer, dimension(DA)::asscomb2
logical, dimension(DA)::chdef2
integer:: combint2
real(doub),dimension(S)::ucons2, ucons3, q_fut2, q_fut3, rep_binary, rep_binary2
real(doub), dimension(A,Z,W,S)::Vns, Ppi, Ppi_updated
real(doub)::welfc

! Policy functions
integer, dimension(:,:,:,:), allocatable :: assfin, assfin2
real(doub), dimension(:,:,:,:), allocatable :: thfin, thfin2
real(doub), dimension(A) :: q_eq
real(doub), dimension(A) :: tempA
integer, dimension(Z,4) :: IC
integer, dimension(A,Z,W)::endint, endint2
real(doub),dimension(A,Z,W,A)::fixnumall
logical :: crisis_light, temp_log
integer :: max_s1, max_s2, max_s, temp_size

! Used in the simulation
integer, parameter:: sim=500, TB = 500 !sim=1500, TB=1000
real(doub)::shockone
real(doub),dimension(Z,A)::amser
real(doub),dimension(Z,A,A)::amdiff
real(doub),dimension(TB, sim)::shockzw, shockh, shocks
integer, dimension(TB)::zpath, apath, wpath, shock_path, history, crisis_path, spath_ind
real(doub), dimension(TB):: sptpathprob
real(doub), dimension(TB):: spath, defpath, csim, ysim, grsim, grsim_real, hat_ysim, nx, debts, qpath,qpath1,qpathRN, qpathRNLT, assets,cred_wealth
real(doub), dimension(TB-1):: sprpatha, rpath, spreadpath, spreadpatha, rpath1, spreadpath1a, spreadpath1,yield_curve, risk_spread, risk_spreadLT, &
    rpathRN, spreadpathaRN, rpathSTRN, spreadpathSTRN, spreadpath_predict, mixpath
logical, dimension(TB-1)::samp1

real(doub),dimension(sim)::defaultpc, debtsvol
real(doub),dimension(sim)::meanspread, meandebt, meandebtserv, Rsquare

real(doub):: meanmeandef, meandebtsvol
real(doub):: meanmeanspr, meanmeandebt, meanmeansv
integer, parameter:: TBd=TB/4
real(doub),dimension(TBd)::debtssim
real(doub),dimension(4)::ds1, ds2, dsint

integer,parameter:: nos=TB-1
real(doub), dimension(nos) :: logy, logc, nxy, ydevs, cdevs, nxy_devs, cgrowth
real(doub), dimension(nos,3) :: vscrap
real(doub), dimension(nos)::spreadt, yfilt,cfilt, spread_devs
real(doub),dimension(4, sim):: stdall1
real(doub),dimension(10,sim):: corrall1
real(doub), dimension(10):: corr1
real(doub), dimension(4)::stdn1

integer, dimension(TB)::inc, inc2, first_pd_def
integer, parameter::maxinc=20

integer::aminloc

! Used in welfare analysis (ergodic distributions)
real(doub), dimension(nzw)::invpdfzw, invpdfzw2
real(doub), dimension(Z) :: invpdfz
real(doub), dimension(W) :: invpdfw
real(doub)::Ez, Ew

! Data moments used as targets
real(doub),parameter:: spreaddata=0.034, debtdata=0.70, sprvoldata=0.025, def_freq_data = .017

real(doub), parameter::one=1.0
integer, parameter::two=2
real(doub),parameter::eps=1e-15

! Used by MPI and in parallelization
integer,parameter:: idmaster=0
integer:: ierr,id,nproc
integer:: itop, iend, nii
integer::i0, i1, i2, i0pr, i1_u

integer, dimension(NZAW):: izloop, ialoop, iwloop
integer::iins, iins2, iins3
integer, dimension(A,Z,W)::iiv
integer, dimension(Z,W) :: iizw
integer, dimension(NZW) :: zwloop_z, zwloop_w

! Used in sunspot analysis
integer, parameter::LL=6
real(doub), dimension(LL)::sptvec, aminvec
logical, dimension(LL) :: conv_yes
real(doub), dimension(LL,11)::funcreg

real(doub), dimension(A) :: spt
real(doub), dimension(A)::ampr
integer, dimension(A)::amnx
real(doub),dimension(TB, sim)::shockp

! Temporary numbers for loops, indices, and computation
integer:: t, i, j, j2, ij, k, int1, int2, int3, int4, int5, int_temp, ieps, i_sa_counter
integer::is, iz, is2, iz2, ia, ia2, ia3, ia4, iw, iw2, izw, izw2, iaa_m_temp, iss_temp, iz_u, ia_u, iw_u, is_u, a_alt
integer::kkk, ss, isind
logical:: defpr, spread_spike
real(doub)::num, num1, num2, num3, num4, num5, Kons1, yy, CRRA_temp, bbeta_0, bbeta_1, qtemp, qtemp_u, xxtemp, stemp, return_u_def, return_u_rep, temp_cons, temp_ppi, num_temp, lender_ret_def
real(doub), dimension(NZW) :: temp
integer, dimension(1:1)::dummy

! Used in the minimization routine for GMM
integer, parameter :: nmax = 1024, mmax = 17
character*60 ::    task, csave
logical, dimension(4) :: lsave
integer  :: n, mm, iprint,&
    nbd(nmax), iwa(3*nmax), isave(44), min_ind
real(doub) ::f, factr, pgtol, &
    x(nmax), xtemp(nmax), l(nmax), u(nmax), gg(nmax), dsave(29), &
    wa(2*mmax*nmax + 5*nmax + 11*mmax*mmax + 8*mmax), grad_tol


end module globvar
