    program main
    implicit none

    ! param
    integer,parameter   :: r_size=4
    integer,parameter   :: ystr=<yearStart>
    integer,parameter   :: yend=<yearEnd>
    integer,parameter   :: ijmax=<ijmax>
    integer,parameter   :: nz=<nz>
    integer,parameter   :: yymax=yend-ystr+1
    real(r_size),parameter   :: vmiss=-999.d0
    ! input
    real(r_size)  :: prcp  (ijmax,12,yymax)   !! precipitation (mm/d)
    real(r_size)  :: prcp01(ijmax,12,yymax)   !! precipitation isotope ratio (ratio * mm/d, SMOW)
    real(r_size)  :: qv    (ijmax,nz,12,yymax)!! humidity (mm/d)
    real(r_size)  :: qv01  (ijmax,nz,12,yymax)!! humidity isotope ratio (ratio * mm/d, SMOW)
    real(r_size)  :: tmp2m (ijmax,12,yymax)   !! 2m temperature (K)
    real(r_size)  :: rh2m  (ijmax,12,yymax)   !! relative humidity (ND)
    real(r_size)  :: ps    (ijmax,12,yymax)   !! surface pressure (hPa)
    real(r_size)  :: npp   (ijmax,12)         !! net primary production
    ! output
    real(r_size)  :: celldw(ijmax,yymax)
    ! internal work
    real(r_size) :: esat, ei, ea, et, es, wi, wa, ws, bp
    real(r_size) :: ef, rsw, ratm, lho, ldho, celld, dv, sw
    real(r_size) :: t, tmp, atmp, rh
    real(r_size) :: wgt, pr
    integer :: mm, yy, it, jt, ij, iz
    integer :: mm2, yy2
    character*128 :: inFilePrcp, inFilePrcpIso, inFileQv, inFileQvIso
    character*128 :: inFileT2, inFileRH2, inFilePs, inFileNPP
    character*128 :: outFile
    real(r_size) :: gtw !! total conductance (stomatal plus boundary-layer) to water vapor (mol m-2 s-1)
    real(r_size) :: gsw !! stomatal conductance to water (mol m-2 s-1)
    real(r_size) :: gbc !! boundary layer conductance to water (mol m-2 s-1)
    real(r_size) :: dt  !! temperature difference b/w air and leaf
    real(r_size) :: ak  !! kinetic frac associated with diffusion in air
    real(r_size) :: akb !! kinetic frac associated with diffusion thru boundary layer
    real(r_size) :: ro  !! standard oxygen isotope ratio (VSMOW)
    real(r_size) :: ff  !! isotope elevation during photosynthesis
    real(r_size) :: fre !! fraction of carbon-bounded oxygen undergoes exchange with medium water
    real(r_size) :: ts  !! steam temperature at standard pressure
    integer :: tcnv     !! e-folding time of source water (month)
    data    ts  / 373.16    /
    data    gsw / 0.0001    /   !<-- isn't it too small?? --> little sensitivity
    data    gbc / 0.0004    /   !<-- isn't it too small?? --> little sensitivity
    data    dt  / 1.d0      /
    !data    ak  / 1.032     /   ! kei
    !data    akb / 1.021     /   ! kei
    data    ak  / 1.0285    /   ! roden et al., 2000
    data    akb / 1.0189    /   ! roden et al., 2000
    data    ro  / 0.0020052 /
    data    ff  / 27.d0     /
    data    fre / 0.42d0    /
    data    tcnv/ 3         /   

    ! get args
    ! check # files
    if( iargc() /= 9 )then
      write(*,*) ' @@@ # of input files are not correct'
      write(*,*) ' ==> Aborted.'
      stop
    endif
    ! get arguments
    call getarg(1,inFilePrcp)
    call getarg(2,inFilePrcpIso)
    call getarg(3,inFileQv)
    call getarg(4,inFileQvIso)
    call getarg(5,inFileT2)
    call getarg(6,inFileRH2)
    call getarg(7,inFilePs)
    call getarg(8,inFileNPP)
    call getarg(9,outFile)

    !! read input files
    open(11,file=inFilePrcp,access='direct',recl=4*ijmax)
    open(12,file=inFilePrcpIso,access='direct',recl=4*ijmax)
    open(13,file=inFileQv,access='direct',recl=4*ijmax*nz)
    open(14,file=inFileQvIso,access='direct',recl=4*ijmax*nz)
    open(15,file=inFileT2,access='direct',recl=4*ijmax)
    open(16,file=inFileRH2,access='direct',recl=4*ijmax)
    open(17,file=inFilePs,access='direct',recl=4*ijmax)
    open(18,file=inFileNPP,access='direct',recl=4*ijmax)
    do yy=1,yymax
      do mm=1,12
        read(13,rec=(yy-1)*12+mm) ((qv  (ij,iz,mm,yy),ij=1,ijmax),iz=1,nz)
        read(14,rec=(yy-1)*12+mm) ((qv01(ij,iz,mm,yy),ij=1,ijmax),iz=1,nz)
      enddo
    enddo
    do yy=1,yymax
      do mm=1,12
        read(11,rec=(yy-1)*12+mm) (prcp  (ij,mm,yy),ij=1,ijmax)
        read(12,rec=(yy-1)*12+mm) (prcp01(ij,mm,yy),ij=1,ijmax)
        read(15,rec=(yy-1)*12+mm) (tmp2m (ij,mm,yy),ij=1,ijmax)
        read(16,rec=(yy-1)*12+mm) (rh2m  (ij,mm,yy),ij=1,ijmax)
        read(17,rec=(yy-1)*12+mm) (ps    (ij,mm,yy),ij=1,ijmax)
      enddo
    enddo
    do mm=1,12
      read(18,rec=mm) (npp(ij,mm),ij=1,ijmax)
    enddo
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)

    !! Calculate leaf conductance of water vapor (gtw) in mol*m-2*s-1
    gtw = 1.d0 / (1.d0/gsw + 1.d0/gbc)

    !! Calculate monthly mean cellulose isotope ratio
    celldw = 0.d0
    do ij = 1, ijmax
      do yy = 1, yymax
        wgt = 0.d0
        do mm = 1, 12
          it = (yy-1) * 12 + mm

          !!! <--- future
          !!! define first month of growing season
          !gsbn = 6
          !!! define last month of growing season
          !gsen = 9
          !!! definition of "water year", which spans from the month after the end of the growing season to the end of the next growing season
          !wybn = gsen - 11
          !wyen = gsbn
          !!! ---> future

          !! define air temperature (K)
          tmp = tmp2m(ij,mm,yy)
          !! define temperature at leaf surface (K)
          atmp = tmp + dt
          !! define growing season relative humidity (ND)
          rh = max(min(rh2m(ij,mm,yy),1.d0),0.d0)
          !! define growing season barometric pressure (Pa)
          bp = ps(ij,mm,yy) * 100
          !! define atmospheric water d18O value
          !! find bottom layer
          do iz = 1, nz
            if ( qv(ij,iz,mm,yy) /= vmiss ) exit
          enddo
          dv = qv01(ij,iz,mm,yy) / qv(ij,iz,mm,yy)
          !! define weighted average source water isotopic value
          sw = 0.d0
          pr = 0.d0
          do jt = it, it-tcnv, -1
            if (jt <= 0 ) cycle
            mm2 = mod(jt-1,12)+1
            yy2 = int((jt-1)/12)+1
            sw = sw + prcp01(ij,mm2,yy2)
            pr = pr + prcp(ij,mm2,yy2)
          enddo
          if ( pr /= 0.d0 ) then
            sw = sw / pr
          else
            cycle
          endif

          !! Calculate Leaf vapor pressure (ei) in kPa (Richards, 1971)
          !t = ts / tk
          t = 1.d0 - ts / atmp
          ei=(101325*EXP((((-0.1299*t-0.6445)*t-1.976)*t+13.3185)*t))/1000
          !! Calculate saturation vapor pressure (esat) in kPa
          t = 1.d0 - ts / tmp
          esat=(101325*EXP((((-0.1299*t-0.6445)*t-1.976)*t+13.3185)*t))/1000
          !! calculate ambient vapor pressure (ea) in kPa
          ea=rh*esat
          !! calculate leaf water vapor fraction (wi) as a molar fraction
          wi=ei/bp
          !! calculate ambient water vapor (wa) as a molar fraction
          wa=ea/bp
          !! calculate leaf transpiration (et) in mol*m-2*s-1
          et=(wi-wa)*gtw
          !! calculate leaf surface water vapor (ws) as a molar fraction
          ws=((gsw*wi)-et*(1-wi/2))/(gsw-et/2)
          !! calculate vapor pressure at leaf surface (es) in kPa
          es=ws*bp
          !! calculate temperature-dependent equilibrium fractionation factor between vapor and liquid water at leaf surface (ef)
          ef=EXP((1.137*1000/atmp/atmp)-(0.4156/atmp)-0.0020667)
          !! calculate isotopic composition of source water as an isotopic ratio, R-value (rsw)
          rsw=ro*sw
          !! calculate isotopic composition of atmospheric water vapor as an isotopic ratio, R-value, (ratm)
          ratm=ro*dv
          !! calculate isotopic composition of leaf water as an isotopic ratio, R-value (lho)
          lho=ef*((ak*rsw*((ei-es)/ei))+(akb*rsw*(es-ea)/ei)+(ratm*ea/ei))
          !! calculate delta value for leaf water (ldho)
          ldho=((lho/ro)-1)*1000
          !! calculate delta value for annual integrated cellulose (celld)
          celld=fre*(sw+ff)+(1-fre)*(ldho+ff)
          !! weight by NPP
          celldw(ij,yy) = celldw(ij,yy) + celld * npp(ij,mm)
          wgt = wgt + npp(ij,mm)
        enddo
        if ( wgt /= 0.d0 ) then
          celldw(ij,yy) = celldw(ij,yy) / wgt
        else
          celldw(ij,yy) = vmiss
        endif
      enddo
    enddo

    !! check nan
    do ij=1,ijmax
      do yy=1,yymax
        if(celldw(ij,yy) /= celldw(ij,yy)) celldw(ij,yy) = vmiss
      enddo
    enddo

    !! output file
    open(21,file=outFile,access='direct',recl=4*ijmax)
    do yy=1,yymax
      write(21,rec=yy)(celldw(ij,yy),ij=1,ijmax)
    enddo
    close(21)

    end

