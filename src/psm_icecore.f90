    program main
    implicit none

    !! param
    integer,parameter :: r_size = 4
    integer,parameter :: ystr = <yearStart>
    integer,parameter :: yend = <yearEnd>
    integer,parameter :: ijmax = <ijmax>
    integer,parameter :: yymax = yend - ystr + 1
    real(r_size),parameter :: vmiss = -999.d0
    !! input
    real(r_size) :: prcp(ijmax,12,yymax)
    real(r_size) :: prcp01(ijmax,12,yymax)
    !! output
    real(r_size) :: ice(ijmax,yymax)
    !! internal work
    real(r_size) :: wgt
    integer :: ij, mm, yy
    character*128 :: inFilePrcp, inFilePrcpIso
    character*128 :: outFile

    !! 
    !! input ::
    !! 
    !! output 
    !! 
    !!  ice :: annual mean precipitation (permil)
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
    !! 
 

    !! get args
    if ( iargc() /= 3 ) then
      write(*,*) ' @@@ # of input files are not correct'
      write(*,*) ' ==> Aborted.'
      stop
    endif
    call getarg(1,inFilePrcp)
    call getarg(2,inFilePrcpIso)
    call getarg(3,outFile)

    !! read input data
    open(11,file=inFilePrcp,access='direct',recl=4*ijmax)
    open(12,file=inFilePrcpIso,access='direct',recl=4*ijmax)
    do yy=1,yymax
      do mm=1,12
        read(11,rec=(yy-1)*12+mm) (prcp(ij,mm,yy),ij=1,ijmax)
        read(12,rec=(yy-1)*12+mm) (prcp01(ij,mm,yy),ij=1,ijmax)
      enddo
    enddo
    close(11)
    close(12)

    !! weight by precipitation amount
    ice = 0.d0
    do ij=1,ijmax
      do yy=1,yymax
        wgt = 0.d0
        do mm=1,12
          ice(ij,yy) = ice(ij,yy) + prcp01(ij,mm,yy)
          wgt = wgt + prcp(ij,mm,yy)
        enddo
        if ( wgt > 0.d0 ) then
          ice(ij,yy) = ( ice(ij,yy) / wgt - 1.d0 ) * 1000.d0
        else
          ice(ij,yy) = vmiss
        endif
      enddo
    enddo


    !! archive model

    !! check nan
    do ij=1,ijmax
      do yy=1,yymax
        if(ice(ij,yy) /= ice(ij,yy)) ice(ij,yy) = vmiss
      enddo
    enddo

    !! output file
    open(21,file=outFile,access='direct',recl=4*ijmax)
    do yy=1,yymax
      write(21,rec=yy)(ice(ij,yy),ij=1,ijmax)
    enddo
    close(21)

    end

