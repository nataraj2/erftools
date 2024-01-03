!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
    program mbm
    implicit none

!-----------------------------------------------------------------------
!
!  mbm - a fortran90 program to generate the initial conditions for the
!        "moist benchmark" case of Bryan and Fritsch (2002)
!
!  Version 1.00                           Last modified:  18 October 2008
!
!  Author:  George H. Bryan
!           Mesoscale and Microscale Meteorology Division
!           National Center for Atmospheric Research
!           Boulder, Colorado, USA
!           gbryan@ucar.edu
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!
!  Reference:  Bryan and Fritsch (2002, MWR, p. 2917)
!
!-----------------------------------------------------------------------
!
!  Special instructions:  This code is designed to generate initial conditions
!  for the Bryan Cloud Model, CM1 (http://www.mmm.ucar.edu/people/bryan/cm1/).
!  Thus, the default thermodynamic equation, formulation of saturation vapor
!  pressure, physical constants, and grid staggering are all appropriate for
!  that code.
!
!  This code should be modified for input to other numerical models.  The
!  following is a list of the sections of code that may need to be changed:
!
!  1) physical constants:  (see "physical constants" section below)
!
!  2) thermodynamic equation:  (see section labeled  "Code your model's 
!     thermodynamic equation here")
!
!  3) function for saturation vapor pressure:  (see section labeled 
!     "Code your model's equation for saturation vapor pressure here)
!
!  4) grid staggering:  a staggered grid is assumed here.  The first
!     scalar level is located one-half of a grid increment above the 
!     ground.  If you need to change this, you must do so in two places:
!     a) in the "base-state profile" section of code, and b) in the
!     two-dimensional section of code. 
!
!  Output is in GrADS format.  This should be easy to change, if another
!  format is desired.
!
!  The default "user settings" below are used to generate the control moist 
!  benchmark case in Bryan and Fritsch (2002):  that is, theta_e = 320 K
!  and r_t = 0.020
!
!-----------------------------------------------------------------------
!  user settings:

    integer, parameter :: ni = 200    ! number of grid points in x direction
    integer, parameter :: nk = 100    ! number of grid points in z direction

    real, parameter :: dx = 100.0     ! grid spacing in x direction (m)
    real, parameter :: dz = 100.0     ! grid spacing in z direction (m)

    real, parameter :: th_sfc = 289.8486  ! potential temperature at surface (K)
    real, parameter :: p_sfc  = 100000.0  ! pressure at surface (Pa)
    real, parameter :: qt_mb  = 0.020     ! total water mixing ratio (kg/kg)

!-----------------------------------------------------------------------
!  physical constants:

    real, parameter :: g      = 9.81
    real, parameter :: to     = 273.15
    real, parameter :: rd     = 287.04
    real, parameter :: rv     = 461.5
    real, parameter :: cp     = 1005.7
    real, parameter :: cpv    = 1870.0
    real, parameter :: p00    = 1.0e5
    real, parameter :: rp00   = 1.0/p00
    real, parameter :: eps    = rd/rv
    real, parameter :: reps   = rv/rd
    real, parameter :: xlv    = 2501000.0
    real, parameter :: cpl    = 4190.0
    real, parameter :: lv1    = xlv+(cpl-cpv)*to
    real, parameter :: lv2    = cpl-cpv

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    integer :: i,k,n,orec
    real :: pi_sfc,t_sfc,qv_sfc,ql_sfc,thv_sfc,thex,t1,th1,qv1,ql1,   &
            p1,pi1,thv1,delz,thlast,th2,thv2,pi2,p2,t2,qv2,ql2,tbar,   &
            qvbar,qlbar,lhv,rm,cpm,diff,xh,zh,beta,thpert
    real, dimension(nk) :: t0,th0,thv0,prs0,pi0,qv0,ql0
    real, dimension(ni,nk) :: tha,prs,qv,ql
    real :: getqvs
    logical :: converged

!-----------------------------------------------------------------------
!  misc constants

    real, parameter :: converge  = 0.0001
    real, parameter :: converge2 = 0.0002
    real, parameter :: pi        = 3.14159265358979323

!-----------------------------------------------------------------------
!  get base state profile:

    print *,'  getting base state profile ...'

    ! surface values:
    pi_sfc  = (p_sfc/p00)**(rd/cp)
    t_sfc   = th_sfc*pi_sfc
    qv_sfc  = getqvs(eps,p_sfc,t_sfc)
    ql_sfc  = qt_mb-qv_sfc
    thv_sfc = th_sfc*(1.0+reps*qv_sfc)/(1.0+qv_sfc+ql_sfc)

    t1   = t_sfc
    th1  = th_sfc
    qv1  = qv_sfc
    ql1  = ql_sfc
    p1   = p_sfc
    pi1  = pi_sfc
    thv1 = th1*(1.0+reps*qv1)/(1.0+qv1+ql1)

    DO k=1,nk

      ! staggered grid:  first grid levels is 1/2 dz above ground
      if(k.eq.1)then
        delz = 0.5*dz
      else
        delz = dz
      endif

      n = 0
      converged = .false.

      thlast = th1
      th2    = th1
      thv2   = thv1

      do while( .not. converged )
        n = n + 1

        pi2  = pi1 - delz*g/(cp*0.5*(thv1+thv2))
        p2   = p00*(pi2**(cp/rd))
        t2   = thlast*pi2
        qv2  = getqvs(eps,p2,t2)
        ql2  = qt_mb - qv2
        thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2)

        tbar  = 0.5*(t1+t2)
        qvbar = 0.5*(qv1+qv2)
        qlbar = 0.5*(ql1+ql2)

        lhv = lv1-lv2*tbar
        rm  = rd+rv*qvbar
        cpm = cp+cpv*qvbar+cpl*qlbar

!-----------------------------------------------------------------------
!  Code your model's thermodynamic equation here:
!-----------------------------------------------------------------------
!  an exact governing equation (e.g., for CM1)
        th2 = th1*exp( -lhv*(qv2-qv1)/(cpm*tbar)     &
                       +(rm/cpm-rd/cp)*alog(p2/p1) )
!-----------------------------------------------------------------------
!  a commonly used approximate governing equation (e.g., for MM5, WRF, etc)
!!!        th2 = th1 - (qv2-qv1)*lhv/(cp*pi2)
!-----------------------------------------------------------------------

        diff = th2-thlast
!!!        print *,n,thlast,th2,diff

        if( abs(diff).gt.0.0001 )then
          thlast = thlast + 0.30*(th2-thlast)
        else
          converged = .true.
        endif

        if(n.gt.100)then
          print *,'  Error:  the solution is not converging'
          stop 1111
        endif

      enddo

      t2   = th2*pi2
      qv2  = getqvs(eps,p2,t2)
      ql2  = qt_mb - qv2
      thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2)

      if(ql2.lt.0.0)then
        print *,'  Warning:  ql < 0'
        stop 1113
      endif

      t0(k) = t2
      th0(k) = th2
      thv0(k) = thv2
      prs0(k) = p2
      pi0(k) = pi2
      qv0(k) = qv2
      ql0(k) = ql2

      t1 = t2
      th1 = th2
      thv1 = thv2
      p1 = p2
      pi1 = pi2
      qv1 = qv2
      ql1 = ql2

    ENDDO

!-----------------------------------------------------------------------
!  get the two-dimensional fields:

    print *,'  getting two-dimensional fields ...'

    do k=1,nk
    do i=1,ni
      prs(i,k) = prs0(k)

      ! staggered grid:
      xh = -0.5*ni*dx + (i-1)*dx + 0.5*dx
      zh = (k-1)*dz + 0.5*dz

      beta=sqrt( ((xh-   0.0)/2000.0)**2    &
                +((zh-2000.0)/2000.0)**2)
      if(beta.lt.1.0)then
        thpert=2.0*(cos(0.5*pi*beta)**2)
        pi2 = (prs(i,k)*rp00)**(rd/cp)
        qv(i,k) = qv0(k)
        ql(i,k) = ql0(k)
        thlast = th0(k)
        converged = .false.
        n = 0
        do while( .not. converged )
          n = n + 1
          tha(i,k)=( (thpert/300.0)+(1.0+qt_mb)/(1.0+qv(i,k)) )  &
                     *thv0(k)*(1.0+qv(i,k))/(1.0+reps*qv(i,k))
          t2 = tha(i,k)*pi2
          qv(i,k) = getqvs(eps,prs(i,k),t2)
          ql(i,k) = qt_mb - qv(i,k)
          if(ql(i,k).lt.0.0)then
            print *,'  Warning:  ql < 0'
            stop 1112
          endif
!!!          print *,n,thlast,tha(i,k),abs(tha(i,k)-thlast)
          if( abs(tha(i,k)-thlast).lt.converge2 )then
            converged = .true.
          else
            thlast = tha(i,k)
          endif
          if(n.gt.100)then
            print *,'  Error:  the solution is not converging'
            stop 1115
          endif
        enddo
      else
        tha(i,k) = th0(k)
         qv(i,k) = qv0(k)
         ql(i,k) = ql0(k)
      endif
    enddo
    enddo

!-----------------------------------------------------------------------
!  write GrADS data file

    print *,'  writing output file ...'

    open(unit=20,file='mbm.dat',form='unformatted',access='direct',recl=4)
    orec = 1

    do n=1,4
    do k=1,nk
    do i=1,ni
      if(n.eq.1) write(20,rec=orec) tha(i,k)
      if(n.eq.2) write(20,rec=orec) prs(i,k)
      if(n.eq.3) write(20,rec=orec)  qv(i,k)
      if(n.eq.4) write(20,rec=orec)  ql(i,k)
      orec=orec+1
    enddo
    enddo
    enddo

!-----------------------------------------------------------------------
!  write GrADS descriptor file

    open(unit=21,file='mbm.ctl')
    write(21,101)
    write(21,102)
    write(21,103)
    write(21,104) ni,0.001*0.5*dx,0.001*dx
    write(21,105)  1,0.001*0.0   ,0.001*dx
    write(21,106) nk,0.001*0.5*dx,0.001*dz
    write(21,107)
    write(21,108) 4
    write(21,109) 'th      ',nk,'potential temperature (K)          '
    write(21,109) 'prs     ',nk,'pressure (Pa)                      '
    write(21,109) 'qv      ',nk,'water vapor mixing ratio (kg/kg)   '
    write(21,109) 'ql      ',nk,'liquid water mixing ratio (kg/kg)  '
    write(21,110)

101 format('dset ^mbm.dat')
102 format('title moist benchmark initial conditions')
103 format('undef -99999999.')
104 format('xdef ',i5,' linear ',f14.6,1x,f14.6)
105 format('ydef ',i5,' linear ',f14.6,1x,f14.6)
106 format('zdef ',i5,' linear ',f14.6,1x,f14.6)
107 format('tdef          1 linear 00Z03JUL2000     1MN')
108 format('vars ',i4)
109 format(a8,1x,i5,' 99 ',a35)
110 format('endvars')

!-----------------------------------------------------------------------

   
    print *,'  ... code completed successfully'
    stop 99999
    end program mbm

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  Code your model's equation for saturation vapor pressure here:
!  This formulation is from Bolton (1980, MWR, p. 1047).

    real function getqvs(eps,p,t)
    implicit none

    real :: eps,p,t,es

    es = 611.2*exp(17.67*(t-273.15)/(t-29.65))
    getqvs = eps*es/(p-es)

    return
    end function getqvs

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
