    ! 
    ! Description ;;
    ! This program is based on the following paper :
    ! Dee et al., 2015, https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015MS000447
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
    ! - annual mean isotope ratio in speleothem in VPDB
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
    real*4  :: prcp(ijmax,12,yymax)
    real*4  :: prcp01(ijmax,12,yymax)
    real*4  :: evap(ijmax,12,yymax)
    real*4  :: evap01(ijmax,12,yymax)
    real*4  :: t2(ijmax,12,yymax)
    ! output
    real*4  :: d_spl(ijmax,yymax) !! annual mean of isotope ration in speleothem (VPDB)
    ! internal work
    real*4  :: avt(yymax)   !! annual mean temperature
    real*4  :: inf,r_inf(yymax) !! infiltration
    real*4  :: drip,r_drip  !! drip water
    real*4  :: alpha        !! equilibrium fractionation factor
    real*4  :: wgt
    real*4  :: tau          !! mean transit time [years]
    integer :: t,yy,mm,ij,dend
    character*10    :: arg1
    character*128 :: inFilePrcp, inFilePrcpIso, inFileEvap, inFileEvapIso
    character*128 :: inFileT2
    character*128 :: outFile
    data    tau  / 0.5 /    !! mean transit time 

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
    if ( iargc() == 7 ) then
      call getarg(7,arg1)
      read(arg1,'(f10.5)') tau
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
      avt=0.d0
      !! set param & initialization
      do yy=1,yymax
        do mm=1,12
          avt(yy)=avt(yy)+t2(ij,mm,yy)
        enddo
      enddo
      avt=avt/12

      ! annual average of infilitration
      do yy=1,yymax
        inf=0.d0
        r_inf(yy)=0.d0
        do mm=1,12
          call getday(dend,mm,yy)
          if( prcp(ij,mm,yy)-evap(ij,mm,yy) > 0.d0 )then
            inf = inf + (prcp(ij,mm,yy)-evap(ij,mm,yy)) * dend    !! mm/d --> mm/mon
            r_inf(yy) = r_inf(yy) + (prcp01(ij,mm,yy)-evap01(ij,mm,yy)) * dend
          endif
        enddo
        if (inf.gt.0.d0) then
          r_inf(yy)=r_inf(yy)/inf
        else
          r_inf(yy)=vmiss
        endif
      enddo


      do yy=1,yymax
        ! convolution
        r_drip=0.d0
        wgt=0.d0
        do t=1,yy
          if (r_inf(t).ne.vmiss) then
            r_drip=r_drip+r_inf(t)*exp(-1*(yy-t+0.5d0)/tau)/tau
            wgt=wgt+exp(-1*(yy-t+0.5d0)/tau)/tau
          endif
        enddo
        if (wgt.gt.0.d0) then
          ! isotopic fractionation (Wackerbarth et al., 2010)
          r_drip=r_drip/wgt
          alpha=exp(2780/avt(yy)/avt(yy)-0.00289) 
          d_spl(ij,yy)=(r_drip*alpha/1.03086-1.d0) * 1000.d0
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
