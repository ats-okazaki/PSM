    ! 
    ! Description ;;
    ! This program is based on the following paper :
    ! Liu, G., Kojima, K., Yoshimura, K., Okai, T., Suzuki, A., Oki, T., Siringan, F. P., Yoneda, M., and Kawahata, H. (2013), A model-based test of accuracy of seawater oxygen isotope ratio record derived from a coral dual proxy method at southeastern Luzon Island, the Philippines, J. Geophys. Res. Biogeosci., 118, 853- 859, doi:10.1002/jgrg.20074.
    ! Liu, G., Kojima, K., Yoshimura, K., and Oka, A. (2014), Proxy interpretation of coral-recorded seawater δ18O using 1-D model forced by isotope-incorporated GCM in tropical oceanic regions, J. Geophys. Res. Atmos., 119, 12,021- 12,033, doi:10.1002/2014JD021583.
    !
    ! 
    ! Inputs ;;
    ! - monthly mean 
    !   -- precipitation
    !   -- precipitation isotope
    !   -- evaporation
    !   -- evaporation isotope
    !   -- sea surface temperature
    ! - climatology
    !   -- surface sea water isotope
    !   -- deeper sea water isotope
    ! 
    ! Outputs ;;
    ! - annual mean isotope ratio in coral calcite
    !
    ! Notes ;;
    ! - resolution of input files must be the same
    !
    program main
    implicit none

    ! param
    integer,parameter   :: r_size=4
    integer,parameter   :: ystr=<yearStart>
    integer,parameter   :: yend=<yearEnd>
    integer,parameter   :: ijmax=<ijmax>
    integer,parameter   :: yymax=yend-ystr+1
    real(r_size),parameter   :: vmiss=-999.d0
    ! input
    real(r_size)  :: prcp  (ijmax,12,yymax)   !! precipitation (mm/d)
    real(r_size)  :: prcp01(ijmax,12,yymax)   !! precipitation isotope ratio (ratio * mm/d, SMOW)
    real(r_size)  :: evap  (ijmax,12,yymax)   !! evaporation (W/m2)
    real(r_size)  :: evap01(ijmax,12,yymax)   !! evaporation isotope ratio (ratio * W/m2, SMOW)
    real(r_size)  :: sst   (ijmax,12,yymax)   !! sea surface temperature (deg-C)
    real(r_size)  :: dow   (ijmax)            !! surface ocean water isotope (permil, SMOW)
    real(r_size)  :: dowd  (ijmax)            !! deep ocean water isotope (permil, SMOW)
    ! output
    real(r_size)  :: d_coral(ijmax,12,yymax)  !! monthly isotope ration in coral (permil, VPDB)
    real(r_size)  :: d_coral2d(ijmax,yymax)   !! annual mean isotope ration in speleothem (permil, VPDB)
    ! internal work
    real(r_size)  :: pterm,eterm,dterm
    real(r_size),parameter :: q=20000         !! depth of upper layer [mm]
    real(r_size),parameter :: d=0.4*q         !! mixing from the lower layer [mm/mon]
    real(r_size),parameter :: a=-0.22         !! dT/d(d_coral)
    real(r_size)  :: d_ow,d_owd
    integer :: t,yy,mm,ij,dend
    character*128 :: inFilePrcp, inFilePrcpIso, inFileEvap, inFileEvapIso
    character*128 :: inFileSST, inFileSSWIso, inFileSDWIso
    character*128 :: outFile

    !! get args
    ! check # files
    if( iargc() /= 8 )then
      write(*,*) ' @@@ # of input files are not correct'
      write(*,*) ' ==> Aborted.'
      stop
    endif
    ! get arguments
    call getarg(1,inFilePrcp)
    call getarg(2,inFilePrcpIso)
    call getarg(3,inFileEvap)
    call getarg(4,inFileEvapIso)
    call getarg(5,inFileSST)
    call getarg(6,inFileSSWIso)
    call getarg(7,inFileSDWIso)
    call getarg(8,outFile)

    !! read input files
    open(11,file=inFilePrcp,access='direct',recl=4*ijmax)
    open(12,file=inFilePrcpIso,access='direct',recl=4*ijmax)
    open(13,file=inFileEvap,access='direct',recl=4*ijmax)
    open(14,file=inFileEvapIso,access='direct',recl=4*ijmax)
    open(15,file=inFileSST,access='direct',recl=4*ijmax)
    open(16,file=inFileSSWIso,access='direct',recl=4*ijmax)
    open(17,file=inFileSDWIso,access='direct',recl=4*ijmax)
    do yy=1,yymax
      do mm=1,12
        read(11,rec=(yy-1)*12+mm) (prcp(ij,mm,yy),ij=1,ijmax)
        read(12,rec=(yy-1)*12+mm) (prcp01(ij,mm,yy),ij=1,ijmax)
        read(13,rec=(yy-1)*12+mm) (evap(ij,mm,yy),ij=1,ijmax)
        read(14,rec=(yy-1)*12+mm) (evap01(ij,mm,yy),ij=1,ijmax)
        read(15,rec=(yy-1)*12+mm) (sst(ij,mm,yy),ij=1,ijmax)
      enddo
    enddo
    read(16,rec=1) (dow(ij),ij=1,ijmax)
    read(17,rec=1) (dowd(ij),ij=1,ijmax)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    !! change unit
    do yy=1,yymax
      do mm=1,12
        do ij=1,ijmax
          prcp(ij,mm,yy)=prcp(ij,mm,yy)*24*60*60    !! kg/m2/s --> mm/d
          prcp01(ij,mm,yy)=prcp01(ij,mm,yy)*24*60*60    !! kg/m2/s --> mm/d
          evap(ij,mm,yy)=evap(ij,mm,yy)/2500000*24*60*60    !! W/m2 --> mm/d
          evap01(ij,mm,yy)=evap01(ij,mm,yy)/2500000*24*60*60    !! W/m2 --> mm/d
        enddo
      enddo
    enddo

    !do ij=1,ijmax
    !  dow(ij)  = dow(ij) / 1000.+1  !! permil --> ratio
    !  dowd(ij) = dowd(ij)/ 1000.+1  !! permil --> ratio
    !enddo

    do ij=1,ijmax
      d_ow=dow(ij)
      d_owd=dowd(ij)
      !d_coral(ij,:)=0.d0
      if ( d_ow == vmiss .or. sst(ij,1,1) == vmiss ) then
        d_coral2d(ij,:) = vmiss
        cycle
      endif
      do yy=1,yymax
        do mm=1,12
          call getday(dend,mm,yy)
          pterm=((prcp01(ij,mm,yy)-prcp(ij,mm,yy))*1000.d0-d_ow*prcp(ij,mm,yy))/q*dend
          eterm=((evap01(ij,mm,yy)-evap(ij,mm,yy))*1000.d0-d_ow*evap(ij,mm,yy))/q*dend
          dterm=(d_owd-d_ow)*d/q
          d_ow=d_ow+pterm-eterm+dterm
          d_coral(ij,mm,yy) = d_ow + a * sst(ij,mm,yy)
          !d_coral(ij,yy)=d_coral(ij,yy)+d_ow+a*(sst(ij,yy,mm)-273.15)
        enddo
        !if (dow(ij).ne.vmiss) then
        !  d_coral(ij,yy)=d_coral(ij,yy)/12
        !else
        !  d_coral(ij,yy)=vmiss
        !endif
        d_coral2d(ij,yy) = 0.d0
        do mm=1,12
          d_coral2d(ij,yy) = d_coral2d(ij,yy) + d_coral(ij,mm,yy)
        enddo
        d_coral2d(ij,yy) = d_coral2d(ij,yy) / 12
      enddo
    enddo
    do ij=1,ijmax
      if (dow(ij).eq.vmiss.or.dowd(ij).eq.vmiss) then
        do yy=1,yymax
          d_coral2d(ij,yy)=vmiss
        enddo
      endif
    enddo

    !! check nan
    do ij=1,ijmax
      do yy=1,yymax
        if(d_coral2d(ij,yy) /= d_coral2d(ij,yy)) d_coral2d(ij,yy) = vmiss
      enddo
    enddo


    !! output file
    open(21,file=outFile,access='direct',recl=4*ijmax)
    do yy=1,yymax
      write(21,rec=yy)(d_coral2d(ij,yy),ij=1,ijmax)
    enddo
    close(21)

    end
        
!**********************************************************************
    subroutine getday(dd,mm,yy)
    implicit none
    !input
    integer :: mm,yy
    !output
    integer :: dd

    if(mm.eq.1.or.mm.eq.3.or.mm.eq.5.or.mm.eq.7.or.mm.eq.8.or.mm.eq.10.or.mm.eq.12)then
      dd=31
    elseif(mm.eq.4.or.mm.eq.6.or.mm.eq.9.or.mm.eq.11)then
      dd=30
    else
      if(mod(yy,4).eq.0)then
        if(mod(yy,100).eq.0)then
          if(mod(yy,400).eq.0)then
            dd=29
          else
            dd=28
          endif
        else
          dd=29
        endif
      else
        dd=28
      endif
    endif

    return
    end subroutine
