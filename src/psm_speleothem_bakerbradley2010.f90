    ! 
    ! Description ;;
    ! This program is based on the following paper :
    ! - Baker, A., Bradley, C. (2010), Modern stalagmite δ18O: Instrumental calibration and forward modelling, GPC, 71, 3-4, 201-206. doi: https://doi.org/10.1016/j.gloplacha.2009.05.002
    !
    ! 
    ! Inputs ;;
    ! - monthly mean 
    !   -- precipitation
    !   -- precipitation isotope
    !   -- evaporation
    !   -- evaporation isotope
    !   -- surface temperature
    ! 
    ! Outputs ;;
    ! - annual mean isotope ratio in speleothem
    !
    ! Notes ;;
    ! - resolution of input files must be the same
    !
    program main
    implicit none

    ! param
    integer,parameter   :: r_size = 4
    integer,parameter   :: ystr=<yearStart>
    integer,parameter   :: yend=<yearEnd>
    integer,parameter   :: ijmax=<ijmax>
    integer,parameter   :: yymax=yend-ystr+1
    real(r_size),parameter   :: vmiss=-999.d0
    ! input
    real(r_size)  :: prcp(ijmax,12,yymax)
    real(r_size)  :: prcp01(ijmax,12,yymax)
    real(r_size)  :: evap(ijmax,12,yymax)
    real(r_size)  :: evap01(ijmax,12,yymax)
    real(r_size)  :: t2(ijmax,12,yymax)
    ! output
    real(r_size)  :: d_spl(ijmax,yymax) !! annual mean of isotope ration in speleothem (VPDB)
    ! internal work
    real(r_size)  :: avpr,avt2
    real(r_size)  :: inf,r_inf    !! infiltration
    real(r_size)  :: storage,r_str    !! storage
    real(r_size)  :: drip,r_drip  !! drip water
    real(r_size)  :: alpha        !! equilibrium fractionation factor
    real(r_size)  :: wgt
    real(r_size)  :: strmax
    real(r_size)  :: k1,k2,init
    real(r_size)  :: r_ini
    integer :: yy,mm,ij,dend
    character*10    :: arg1,arg2,arg3
    character*128 :: inFilePrcp, inFilePrcpIso, inFileEvap, inFileEvapIso
    character*128 :: inFileT2
    character*128 :: outFile
    data    k1  / 6.0 /     ! maximum storage (storage * annual prcp)
    data    k2  / 0.1 /     ! outlet amout (outlet * storage)
    data    init/ 0.66/     ! portion of storage initially full
    data    r_ini/ 0.985 /  ! initil iso in soil

    !! get args
    ! check # files
    if( iargc() < 6 )then
      write(*,*) ' @@@ # of input files are not correct'
      write(*,*) ' ==> Aborted.'
      stop
    endif
    ! get arguments
    call getarg(1,inFilePrcp)
    call getarg(2,inFilePrcpIso)
    call getarg(3,inFileEvap)
    call getarg(4,inFileEvapIso)
    call getarg(5,inFileT2)
    call getarg(6,outFile)
    if ( iargc() == 9 ) then
      call getarg(7,arg1)
      call getarg(8,arg2)
      call getarg(9,arg3)
      read(arg1,'(f10.5)') k1
      read(arg2,'(f10.5)') k2
      read(arg3,'(f10.5)') init
    endif

    open(11,file=inFilePrcp,access='direct',recl=4*ijmax)
    open(12,file=inFilePrcpIso,access='direct',recl=4*ijmax)
    open(13,file=inFileEvap,access='direct',recl=4*ijmax)
    open(14,file=inFileEvapIso,access='direct',recl=4*ijmax)
    open(15,file=inFileT2,access='direct',recl=4*ijmax)
    do yy=1,yymax
      do mm=1,12
        read(11,rec=(yy-1)*12+mm) (prcp(ij,mm,yy),ij=1,ijmax)
        read(12,rec=(yy-1)*12+mm) (prcp01(ij,mm,yy),ij=1,ijmax)
        read(13,rec=(yy-1)*12+mm) (evap(ij,mm,yy),ij=1,ijmax)
        read(14,rec=(yy-1)*12+mm) (evap01(ij,mm,yy),ij=1,ijmax)
        read(15,rec=(yy-1)*12+mm) (t2(ij,mm,yy),ij=1,ijmax)
      enddo
    enddo
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    !! unit
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

    do ij=1,ijmax
      !! set param & initialization
      r_str=0.
      avpr=0.
      avt2=0.
      do yy=1,yymax
        do mm=1,12
          avpr=avpr+prcp(ij,mm,yy)
          avt2=avt2+t2(ij,mm,yy)
          r_str=r_str+prcp01(ij,mm,yy)
        enddo
      enddo
      avpr=avpr/yymax*365
      avt2=avt2/(12*yymax)
      r_str=(r_str/yymax*365)/avpr
      if(avpr.le.0.d0)then
        r_str=r_ini
        write(*,*) 'R_str is set to R_ini'
      endif
      strmax=k1*avpr
      storage=init*strmax

      do yy=1,yymax
        wgt=0.d0
        do mm=1,12
          call getday(dend,mm,yy)
          !infiltration
          if (prcp(ij,mm,yy).gt.evap(ij,mm,yy).and.storage.lt.strmax) then
            inf=(prcp(ij,mm,yy)-evap(ij,mm,yy))*dend    !! mm/d --> mm/mon
            r_inf=(prcp01(ij,mm,yy)-evap01(ij,mm,yy))/inf*dend

            storage=storage+inf
            r_str=(storage*r_str+inf*r_inf)/(storage+inf)

          ! drip
            drip=k2*storage
            r_drip=r_str
            storage=min(max(storage-drip,0.d0),strmax)

          ! isotopic equilibrium
          ! Kim and O'Neil (1997)
            alpha=exp(18.03/avt2-0.03242)
          ! Friedman and O'Neil (1977)
          !  alpha=exp(2780/avt2/avt2-0.00289) 
            d_spl(ij,yy)=d_spl(ij,yy)+drip*(r_drip/1.03092*alpha-1.)*1000.
            wgt=wgt+drip
          endif
        enddo
        if (wgt.gt.0.d0) then
          d_spl(ij,yy)=d_spl(ij,yy)/wgt
        else
          d_spl(ij,yy)=vmiss
        endif
      enddo
    enddo

    open(21,file=outFile,access='direct',recl=4*ijmax)
    do yy=1,yymax
      write(21,rec=yy)(d_spl(ij,yy),ij=1,ijmax)
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
