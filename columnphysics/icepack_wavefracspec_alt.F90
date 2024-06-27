


     module icepack_wavefracspec_alt

      use icepack_kinds
      use icepack_parameters, only: p01, p5, c0, c1, c2, c3, c4, c10
      use icepack_parameters, only: bignum, puny, gravit, pi, rhow, rhoi
      use icepack_tracers, only: nt_fsd
      use icepack_warnings, only: warnstr, icepack_warnings_add,  icepack_warnings_aborted
      use icepack_fsd

      implicit none
      private
      public :: icepack_step_wavefracture_alt

      real (kind=dbl_kind), parameter  :: &
         young_mod  = 10e9, &          ! Youngs Modulus for ice (Pa)
         straincrit = 3.e-5_dbl_kind, & ! critical strain
         dx = c1 ! domain spacing

!=======================================================================

      contains


!=======================================================================
!autodocument_start icepack_step_wavefracture
!
!  Given fracture histogram computed from local wave spectrum, evolve
!  the floe size distribution
!
!  authors: 2018 Lettie Roach, NIWA/VUW
!
      subroutine icepack_step_wavefracture_alt(wave_spec_type,   &
                  dt,            ncat,            nfsd,      &
                  nfreq,                                     &
                  aice,          vice,            aicen,     &
                  floe_rad_l,    floe_rad_c,                 &
                  wave_spectrum, wavefreq,        dwavefreq, &
                  trcrn,         d_afsd_wave)


      character (len=char_len), intent(in) :: &
         wave_spec_type   ! type of wave spectrum forcing

      integer (kind=int_kind), intent(in) :: &
         nfreq,        & ! number of wave frequency categories
         ncat,         & ! number of thickness categories
         nfsd            ! number of floe size categories

      real (kind=dbl_kind), intent(in) :: &
         dt,           & ! time step
         aice,         & ! ice area fraction
         vice            ! ice volume per unit area

      real (kind=dbl_kind), dimension(ncat), intent(in) :: &
         aicen           ! ice area fraction (categories)

      real(kind=dbl_kind), dimension(:), intent(in) ::  &
         floe_rad_l,   & ! fsd size lower bound in m (radius)
         floe_rad_c      ! fsd size bin centre in m (radius)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         wavefreq,     & ! wave frequencies (s^-1)
         dwavefreq       ! wave frequency bin widths (s^-1)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         wave_spectrum   ! ocean surface wave spectrum as a function of frequency
                         ! power spectral density of surface elevation, E(f) (units m^2 s)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcrn           ! tracer array

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         d_afsd_wave     ! change in fsd due to waves

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         d_afsdn_wave    ! change in fsd due to waves, per category

!autodocument_end
      ! local variables
      integer (kind=int_kind) :: &
         n, k,  &
         nsubt ! number of subcycles

      real (kind=dbl_kind), dimension (nfsd, nfsd) :: &
         frac, fracture_hist

      real (kind=dbl_kind) :: &
         hbar         , & ! mean ice thickness
         elapsed_t    , & ! elapsed subcycling time
         subdt        , & ! subcycling time step
         cons_error       ! area conservation error

      real (kind=dbl_kind), dimension (nfsd) :: &
         omega, loss, gain, &
         afsd_init    , & ! tracer array
         afsd_tmp     , & ! tracer array
         d_afsd_tmp       ! change
   character(len=*),parameter :: &
         subname='(icepack_step_wavefracture_alt)'


      !------------------------------------

      ! initialize
      d_afsd_wave    (:)   = c0
      d_afsdn_wave   (:,:) = c0
      fracture_hist  (:,:)   = c0

!      hbar = c1
!      wave_spectrum = (/0.000000000000000E+000,  4.048108530696481E-005, &
!  5.282969796098769E-004,  1.064894371666014E-003,  1.249741762876511E-003, &
!  1.229491783306003E-003,  1.080903923138976E-003,  8.635559934191406E-004, &
!  9.837691904976964E-004,  1.176746212877333E-003,  2.027775160968304E-003, &
!  4.147783387452364E-003,  8.442047052085400E-003,  3.563777357339859E-002, &
!  5.805501341819763E-002,  2.729533798992634E-002,  7.663844618946314E-003, &
!  1.658817403949797E-003,  7.883401121944189E-004,  4.551284946501255E-004, &
!  4.689317429438233E-004,  9.280506637878716E-004,  5.240151658654213E-004, &
!  5.421090172603726E-004,  5.024557467550039E-004/)

!      print *, 'wave_spec ',wave_spectrum
!      print *, 'hbar = ',hbar

      ! if all ice is not in first floe size category
      if (.NOT. ALL(trcrn(nt_fsd,:).ge.c1-puny)) then
      if ((aice > p01).and.(MAXVAL(wave_spectrum(:)) > puny)) then

          hbar = vice/aice ! note- average thickness
          DO k = 2, nfsd
              if (.NOT. ALL(trcrn(nt_fsd+k-1,:).ge.c1-puny)) then
                  call solve_yt_for_strain(nfsd, nfreq, & 
                               floe_rad_l, floe_rad_c, &
                               wavefreq, dwavefreq, &
                               c2*floe_rad_c(k), &
                               hbar, wave_spectrum, & 
                               fracture_hist(k,:))
 
              end if          
          END DO
          if (MAXVAL(fracture_hist) > puny) then
            ! protect against small numerical errors
            call icepack_cleanup_fsd (ncat, nfsd, trcrn(nt_fsd:nt_fsd+nfsd-1,:) )
            if (icepack_warnings_aborted(subname)) return

             DO n = 1, ncat

              afsd_init(:) = trcrn(nt_fsd:nt_fsd+nfsd-1,n)

              ! if there is ice, and a FSD, and not all ice is the smallest floe size
              if ((aicen(n) > puny) .and. (SUM(afsd_init(:)) > puny) &
                                    .and.     (afsd_init(1) < c1)) then

                  afsd_tmp =  afsd_init
                  loss(:) = c0
                  gain(:) = c0
                  omega(:) = c0
                  DO k = 1, nfsd
                      omega(k) = afsd_tmp(k)*SUM(fracture_hist(k,1:k))
                      loss(k) = omega(k)
                  END DO

                  DO k = 1, nfsd
                      gain(k) = SUM(omega(:)*fracture_hist(:,k))
                  END DO
                  afsd_tmp = afsd_tmp + gain -loss

                  ! update trcrn
                  trcrn(nt_fsd:nt_fsd+nfsd-1,n) = afsd_tmp/SUM(afsd_tmp)
                  call icepack_cleanup_fsd (ncat, nfsd, trcrn(nt_fsd:nt_fsd+nfsd-1,:) )
                  if (icepack_warnings_aborted(subname)) return

                  ! for diagnostics
                  d_afsdn_wave(:,n) = afsd_tmp(:) - afsd_init(:)
                  d_afsd_wave (:)   = d_afsd_wave(:) + aicen(n)*d_afsdn_wave(:,n)
 
              end if
             END DO

      end if
      end if
      end if

      end subroutine icepack_step_wavefracture_alt


!===========================================================================
!
      subroutine four_by_four_matrix_solver(a,b,c)

      real (kind=dbl_kind), dimension(4,4), intent(in) :: &
         a ! four by four matrix

      real (kind=dbl_kind), dimension(4), intent(in) :: &
         b ! four-element column vector

      real (kind=dbl_kind), dimension(4), intent(out) :: &
         c ! four-element column vector

      !--- local
      
      real (kind=dbl_kind) :: &
          deno, deno1, deno2, detaa

      real (kind=dbl_kind), dimension (2,2) :: &
          aa, aainv

      real (kind=dbl_kind), dimension (2) :: &
          bb


      deno  = a(2,2) - a(1,2)*a(2,1)/a(1,1) 
      deno1 = a(1,1) * deno
      deno2 = a(1,1) * deno1 
    
      aa(1,1) = -a(3,1)*a(1,2)*a(2,1)*a(1,3)/deno2 + a(3,1)*a(1,2)*a(2,3)/deno1 &
      - a(3,1)*a(1,3)/a(1,1) + a(3,2)*a(2,1)*a(1,3)/deno1 &
      - a(3,2)*a(2,3)/deno + a(3,3)

      aa(1,2) = -a(3,1)*a(1,2)*a(2,1)*a(1,4)/deno2 + a(3,1)*a(1,2)*a(2,4)/deno1 &
      - a(3,1)*a(1,4)/a(1,1) + a(3,2)*a(2,1)*a(1,4)/deno1 &
      - a(3,2)*a(2,4)/deno + a(3,4)

      aa(2,1) = -a(4,1)*a(1,2)*a(2,1)*a(1,3)/deno2 + a(4,1)*a(1,2)*a(2,3)/deno1 &
      - a(4,1)*a(1,3)/a(1,1) + a(4,2)*a(2,1)*a(1,3)/deno1 &
      - a(4,2)*a(2,3)/deno + a(4,3)

      aa(2,2) = -a(4,1)*a(1,2)*a(2,1)*a(1,4)/deno2 + a(4,1)*a(1,2)*a(2,4)/deno1 &
      - a(4,1)*a(1,4)/a(1,1) + a(4,2)*a(2,1)*a(1,4)/deno1 &
      - a(4,2)*a(2,4)/deno + a(4,4)

      bb = (/b(3) - a(3,1)*b(1)/a(1,1) + a(3,1)*a(1,2)*b(2)/deno1 &
     - a(3,1)*a(1,2)*a(2,1)*b(1)/deno2 - a(3,2)*b(2)/deno &
     + a(3,2)*a(2,1)*b(1)/deno1, &
     b(4) - a(4,1)*b(1)/a(1,1) + a(4,1)*a(1,2)*b(2)/deno1 &
     - a(4,1)*a(1,2)*a(2,1)*b(1)/deno2 - a(4,2)*b(2)/deno &
     + a(4,2)*a(2,1)*b(1)/deno1 /)


      detaa = aa(1,1)*aa(2,2)-aa(1,2)*aa(2,1)
      aainv(1,1) = aa(2,2) 
      aainv(1,2) = -aa(1,2) 
      aainv(2,1) = -aa(2,1) 
      aainv(2,2) = aa(1,1)
      aainv = aainv / detaa

      c(3:4) = MATMUL(aainv,bb)
      
      c(1) = (b(1) - a(1,2)/deno * (b(2) - a(2,1)/a(1,1) * &
                (b(1) - a(1,3)*c(3) - a(1,4)*c(4)) &
              - a(2,3)*c(3) - a(2,4)*c(4)) - a(1,3)*c(3) &
              - a(1,4)*c(4)) / a(1,1) 

      c(2) = (b(2) - a(2,1)/a(1,1) * (b(1) - a(1,3)*c(3) - a(1,4)*c(4)) &
              - a(2,3)*c(3) - a(2,4)*c(4)) / deno 

      end subroutine four_by_four_matrix_solver

!===========================================================================
!
      subroutine alt_get_fraclengths(nfsd, floe_rad_c, floe_rad_l, &
                                     x, strain, frac_local)
 
      integer (kind=int_kind), intent(in) :: &
          nfsd

      real (kind=dbl_kind), dimension(:), intent(in) :: &
          x, strain, floe_rad_c, floe_rad_l

      real(kind=dbl_kind), dimension(:), intent(inout) :: &
          frac_local ! binned histogram of fractures

      ! local
      logical (kind=log_kind), dimension (:), allocatable :: &
          exceed_crit_pos, &
          exceed_crit_neg

      integer (kind=int_kind), dimension(:), allocatable :: &
          extremelocs !minlocs, maxlocs

      real (kind = dbl_kind), dimension(:), allocatable :: &
          fraclengths

      integer (kind=int_kind) :: &
        nx, n_exceed, j_beg, j_end, j, jj, k, nfrac
 
      print *, 'BEGIN get fraclengths'

      nx = SIZE(strain)
      allocate(exceed_crit_pos (nx))
      allocate(exceed_crit_neg (nx))
      exceed_crit_pos(:) = .false.
      exceed_crit_neg(:) = .false.

      WHERE (strain.gt.straincrit) exceed_crit_pos = .true.
      WHERE (strain.lt.-straincrit) exceed_crit_neg = .true.
      n_exceed = COUNT(exceed_crit_pos) + COUNT(exceed_crit_neg)
      allocate(extremelocs(n_exceed))

      j_beg = 0
      j_end = 0
      j = 1
      k = 0
      DO WHILE (j<nx)
          if (exceed_crit_neg(j)) then
              j_beg = j
              j_end = j

              DO jj = 1, nx-j
                  if (exceed_crit_neg(j+jj)) then
                      j_end = j+jj
                  else
                      EXIT
                  end if
              END DO
              k = k + 1
              extremelocs(k) = MINLOC(strain(j_beg:j_end),DIM=1)+j_beg-1
              j = j_end + 1 ! skip to end of segment
          else if (exceed_crit_pos(j)) then
              j_beg = j
              j_end = j

              DO jj = 1, nx-j
                  if (exceed_crit_pos(j+jj)) then
                      j_end = j+jj
                  else
                      EXIT
                  end if
              END DO
              k = k + 1
              extremelocs(k) = MAXLOC(strain(j_beg:j_end),DIM=1)+j_beg-1
              j = j_end + 1 ! skip to end of segment

          else
              j = j + 1 ! move to next point
          end if

      END DO

      nfrac = COUNT(extremelocs>0)
      if (nfrac.eq.0) stop 'need to deal with 0 fracture case'
      allocate(fraclengths(nfrac+1))

      fraclengths(1) = X(extremelocs(1)) - X(1) 
      do k = 2, nfrac
          fraclengths(k) = X(extremelocs(k)) - X(extremelocs(k-1))
      end do
      fraclengths(nfrac+1) = X(nx) - X(extremelocs(nfrac))
      print *, 'fraclengths',fraclengths

      frac_local(:) = c0

      ! convert from diameter to radii
      fraclengths(:) = fraclengths(:)/c2

      if (.not. ALL(fraclengths.lt.floe_rad_l(1))) then
          ! bin into FS cats
          do j = 1, size(fraclengths)
              if (fraclengths(j).gt.floe_rad_l(1)) then
                   do k = 1, nfsd-1
                       if ((fraclengths(j) >= floe_rad_l(k)) .and. &
                        (fraclengths(j) < floe_rad_l(k+1))) then
                            frac_local(k) = frac_local(k) + 1
                       end if
                   end do
                if (fraclengths(j)>floe_rad_l(nfsd)) frac_local(nfsd) = frac_local(nfsd) + 1
              end if
          end do

          do k = 1, nfsd
               frac_local(k) = floe_rad_c(k)*frac_local(k)
          end do

          ! normalize
          if (SUM(frac_local) /= c0) frac_local(:) = frac_local(:) / SUM(frac_local(:))

      end if
      print *, 'frac_local', frac_local


      end subroutine alt_get_fraclengths

!===========================================================================
!
      subroutine solve_yt_for_strain(nfsd, nfreq, &
                                     floe_rad_l, floe_rad_c, &
                                     wavefreq, dwavefreq, &
                                     L, hbar, spec_efreq, &
                                     frac_local)
      
      integer(kind=int_kind), intent(in) :: &
           nfreq, & ! number of wave frequencies
           nfsd     ! number of floe size categories

      real (kind = dbl_kind), intent (in) :: &
           L, & ! floe diameter (m)
           hbar ! floe thickness (m)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         wavefreq,   & ! wave frequencies (s^-1)
         dwavefreq,  & ! wave frequency bin widths (s^-1)
         spec_efreq, & ! wave spectrum (m^2 s)
         floe_rad_c, & ! floe radius center (m)
         floe_rad_l    ! floe radius lower bin edge (m)


      real (kind=dbl_kind), dimension(:), intent(inout) :: &
           frac_local

      real (kind=dbl_kind), dimension(:),allocatable :: &!, intent(out) :: &
           strain, strain_yP

      ! local variables
      integer (kind = int_kind) :: &
           Lint, nx, & ! length, number of points in domain
           j ! index

      real (kind = dbl_kind), dimension (:), allocatable :: &
           x, xp,  & ! spatial domain 
           yH, yP, &  ! homogenous, particular SSH solution
           yppH, yppP, & ! second deriv SSH for each solution
           ypp           ! second deriv SSH for total solution

      real (kind = dbl_kind) :: &
           I, & ! moment of inertia [m^3]
           Lambda, & ! characteristic length scale [m]
           gamm, &   ! non-dimensional number
           m, b, ap, bp

       real (kind = dbl_kind), dimension(nfreq) :: &
           lamdai, &  ! wavelengths [m]
           spec_coeff, & ! spectral coefficients
           langi, &  ! rescaled wavelength
           AAmi, &      ! rescaled spectral coefficients
           PHIi         ! phase for SSH

       real (kind=dbl_kind), dimension(4,4) :: &
           aa ! 4x4 matrix

       real (kind=dbl_kind), dimension(4) :: &
           bb, cc ! column vectors

       real (kind=dbl_kind), dimension(:,:), allocatable :: &
           arg

       Lint = NINT(L)
       nx = NINT(Lint/dx+dx)  
       allocate(x(nx))
       allocate(xp(nx))
       allocate(yP(nx))
       allocate(yH(nx))
       allocate(ypp(nx))
       allocate(strain(nx))
       allocate(strain_yP(nx))


       ! dispersion relation
       lamdai (:) = gravit/(c2*pi*wavefreq (:)**2)

       ! spectral coefficients
       spec_coeff = sqrt(c2*spec_efreq*dwavefreq)

       
       DO j=1,nx
           x(j) = -Lint/c2 + (j-1)*dx
       END DO

       ! this should be the same each run
       ! and for restarts
       CALL RANDOM_NUMBER(PHIi)
       PHIi = c2*pi*PHIi

!       PHIi = c2*pi*(/3.920868194323862E-007,  2.548044275764261E-002, &
!  0.352516161261067,       0.666914481524251,       0.963055531894656, &     
!  0.838288203465982,       0.335355043646496,       0.915327203368213, &     
!  0.795863676652503,       0.832693143644796,       0.345042693116063, &     
!  0.871183932316783,       8.991835668825542E-002,  0.888283839684037, &     
!  0.700978902440147,       0.734552583860683,       0.300175817923128, &     
!  4.971772349719251E-002,  0.908189377373128,       9.765859753870422E-002, &
!  4.031338096905369E-002,  8.502479466940610E-002,  0.558820973383161, &     
!  0.926451747654190,       7.564077406631106E-002/)

       I=hbar**3/12
       Lambda = (young_mod*I/(rhow*gravit))**(0.25_dbl_kind)

       gamm = L/(c2*SQRT(c2)*Lambda)
       langi = lamdai/(c2*pi)
       AAmi = spec_coeff*langi**4/(Lambda**4 + langi**4)
       xp = x/(SQRT(c2)*Lambda)

       ! floating line
       m = 6./(L**2)*SUM(spec_coeff*lamdai/pi*sin(PHIi)*(-cos(pi*L/lamdai)+lamdai/(pi*L)*sin(pi*L/lamdai)))
       b = (c1/L)*SUM(spec_coeff*lamdai/pi*sin(pi*L/lamdai)*cos(PHIi))- rhoi/rhow*hbar
       bp = -b - rhoi/rhow*hbar
       ap = -m

       aa(1,1) = EXP(-gamm)*SIN(gamm)
       aa(1,2) = EXP(-gamm)*COS(gamm)
       aa(1,3) = -EXP(gamm)*SIN(gamm)
       aa(1,4) = -EXP(gamm)*COS(gamm)

       aa(2,1) = -EXP(gamm)*SIN(gamm)
       aa(2,2) = EXP(gamm)*COS(gamm)
       aa(2,3) = EXP(-gamm)*SIN(gamm)
       aa(2,4) = -EXP(-gamm)*COS(gamm)

       aa(3,1) = SIN(gamm)*COSH(gamm) + COS(gamm)*SINH(gamm)
       aa(3,2) = -COS(gamm)*SINH(gamm) + SIN(gamm)*COSH(gamm)
       aa(3,3) = SIN(gamm)*COSH(gamm) + COS(gamm)*SINH(gamm)
       aa(3,4) = -SIN(gamm)*COSH(gamm) + COS(gamm)*SINH(gamm)

       aa(4,1) = gamm*SIN(gamm)*SINH(gamm) + gamm*COS(gamm)*COSH(gamm) - SIN(gamm)*COSH(gamm)
       aa(4,2) = gamm*SIN(gamm)*SINH(gamm) - gamm*COS(gamm)*COSH(gamm) + COS(gamm)*SINH(gamm)
       aa(4,3) = -gamm*COS(gamm)*COSH(gamm) - gamm*SIN(gamm)*SINH(gamm) + SIN(gamm)*COSH(gamm)
       aa(4,4) = -gamm*COS(gamm)*COSH(gamm) + gamm*SIN(gamm)*SINH(gamm) + COS(gamm)*SINH(gamm)

       bb(1) = Lambda**2*SUM(AAmi/(langi**2)*COS(L/(c2*langi)+PHIi))
       bb(2) = Lambda**2*SUM(AAmi/(langi**2)*COS(L/(c2*langi)-PHIi))
       
       bb(3) = - SQRT(c2)/(c2*Lambda)*(SUM(AAmi*langi* ( SIN(L/(c2*langi)-PHIi) +&
               SIN(L/(c2*langi)+PHIi) )) + bp*L)

       bb(4) = - c1/(c2*Lambda**2)*(SUM(AAmi*langi* &
               ( L/c2*(SIN(L/(c2*langi)-PHIi) - SIN(L/(c2*langi)+PHIi)) + &
                 langi*(COS(L/(c2*langi)-PHIi) - COS(L/(c2*langi)+PHIi)))) + ap*L**3/12)


       call four_by_four_matrix_solver(aa,bb,cc)

       yH = EXP(xp)*(cc(1)*COS(xp)+cc(2)*SIN(xp)) + EXP(-xp)*(cc(3)*COS(xp)+cc(4)*SIN(xp))
       yppH = (EXP(xp)*(-cc(1)*SIN(xp)+cc(2)*COS(xp)) - EXP(-xp)*(-cc(3)*SIN(xp)+cc(4)*COS(xp)))/Lambda**2

       allocate(arg(nfreq,nx))
       DO j=1,nx
           arg(:,j) = x(j)/langi(:) - PHIi(:)
       END DO

       yP = MATMUL(AAmi,COS(arg))+(ap*x+bp)
       yppP = -MATMUL(AAmi/langi**2,COS(arg))

       ypp = yppP + yppH

       strain = hbar*ypp/c2
       strain_yP = hbar*yppP/c2

       ! only consider particular solution for floes>300m
       if (L.gt.300_dbl_kind) then
              print *, 'max strain yp=',MAXVAL(ABS(strain_yP))

              strain = strain_yP
       end if

       print *, 'max strain=',MAXVAL(ABS(strain))

       if (MAXVAL(ABS(strain)).gt.straincrit) then
           print *, 'condition true'
           call alt_get_fraclengths(nfsd, floe_rad_c, floe_rad_l, &
                                    x,strain, frac_local)
       end if


       end subroutine solve_yt_for_strain

!=======================================================================


!=======================================================================

      end module icepack_wavefracspec_alt

!=======================================================================


