module akmod
  logical :: ak_initialized =.false.

  real(kind=8),target               :: depak(136),vpak(136),vsak(136),dak(136)
  real(kind=8),dimension(:),pointer :: vak

contains
 
  subroutine ak_init
    
    if (ak_initialized) return

! read in the ak135 velocity model
    open(1,file='ak135.dat')
    do i=1,136
       read(1,*) depak(i),vpak(i),vsak(i),dak(i)
    end do
    print *,'ak135 models read in'
    ak_initialized=.true.

  end subroutine ak_init

end module akmod



! the functions below return the velocity in each region

function vel1(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel1,dep,lat,long,deppth,deppth2,deppth3,deppth4
  real(kind=8)  :: v1p,v2p,v3p,v4p,v5p,v1s,v2s,v3s,v4s,v5s,delt,po
  integer        :: vtype
  REAL :: str,dip,dpd,la,lo

  call ak_init

!-----------------------------!
!    SELECTION OF VP OR VS    !
!-----------------------------!

  select case (vtype)
     case(1)
        vak=>vpak
     case(2)
        vak=>vsak
     case default
        stop 'vtype can only be 1 or 2'
  end select

!-----------------------!
!    SYNTHETIC CASES    !
!-----------------------!

  str=000*3.1415/180
  dip1=10*3.1415/180 !15
  dip2=20*3.1415/180 !15
  dip3=30*3.1415/180 !15
  dpd=111.32
  la=lat-1   !38
  lo=long-1  !23
  delt=10

  v1p = 8.0 !5.8
  v1s = 5.0 !3.5
  deppth1=45
  v2p = 8.8 !8.5
  v2s = 5.5 !4.8
  deppth2=070 + SIN(str)*la*dpd*TAN(dip1) + COS(str)*lo*dpd*TAN(dip1)
  v3p = 8.0 !5.8
  v3s = 5.0 !3.5
  deppth3=105 + SIN(str)*la*dpd*TAN(dip2) + COS(str)*lo*dpd*TAN(dip2)
  v4p = 8.8 !9.5
  v4s = 5.5 !5.2
  deppth4=150 + SIN(str)*la*dpd*TAN(dip3) + COS(str)*lo*dpd*TAN(dip3)
  v5p = 8.0 !8.5
  v5s = 5.0 !4.8

  !deppth1=100 + SIN(str)*la*dpd*TAN(dip) + COS(str)*lo*dpd*TAN(dip)
  !deppth2=158+(lat-5)*111.320*TAN(30*3.14159265/180)
  !deppth3=169.5+(lat-5)*111.320*TAN(30*3.14159265/180)
  !deppth4=227.2+(lat-5)*111.320*TAN(30*3.14159265/180)


  select case (vtype)
     case(1)

        vel1 = v1p

        if ( ABS(dep-deppth1) .LT. delt ) then
           vel1 = (deppth1+delt-dep)/(2*delt) * vel1 + &
                   (dep-deppth1+delt)/(2*delt) * v2p
        else if ( dep .GT. deppth1 ) then
           vel1 = v2p
        endif

        if ( ABS(dep-deppth2) .LT. delt ) then
           vel1 = (deppth2+delt-dep)/(2*delt) * vel1 + &
                   (dep-deppth2+delt)/(2*delt) * v3p
        else if ( dep .GT. deppth2 ) then
           vel1 = v3p
        endif

        if ( ABS(dep-deppth3) .LT. delt ) then
           vel1 = (deppth3+delt-dep)/(2*delt) * vel1 + &
                   (dep-deppth3+delt)/(2*delt) * v4p
        else if ( dep .GT. deppth3 ) then
           vel1 = v4p
        endif

        if ( ABS(dep-deppth4) .LT. delt ) then
           vel1 = (deppth4+delt-dep)/(2*delt) * vel1 + &
                   (dep-deppth4+delt)/(2*delt) * v5p
        else if ( dep .GT. deppth4 ) then
           vel1 = v5p
        endif

     case(2)

        vel1 = v1s

        if ( ABS(dep-deppth1) .LT. delt ) then
           vel1 = (deppth1+delt-dep)/(2*delt) * vel1 + &
                   (dep-deppth1+delt)/(2*delt) * v2s
        else if ( dep .GT. deppth1 ) then
           vel1 = v2s
        endif

        if ( ABS(dep-deppth2) .LT. delt ) then
           vel1 = (deppth2+delt-dep)/(2*delt) * vel1 + &
                   (dep-deppth2+delt)/(2*delt) * v3s
        else if ( dep .GT. deppth2 ) then
           vel1 = v3s
        endif
         
        if ( ABS(dep-deppth3) .LT. delt ) then
           vel1 = (deppth3+delt-dep)/(2*delt) * vel1 + &
                   (dep-deppth3+delt)/(2*delt) * v4s
        else if ( dep .GT. deppth3 ) then
           vel1 = v4s
        endif
         
        if ( ABS(dep-deppth4) .LT. delt ) then
           vel1 = (deppth4+delt-dep)/(2*delt) * vel1 + &
                   (dep-deppth4+delt)/(2*delt) * v5s
        else if ( dep .GT. deppth4 ) then
           vel1 = v5s
        endif

  end select

  return

end function vel1


function vel2(dep,lat,long,vtype)
  use akmod
  real(kind=8)   :: vel2,dep,lat,long
  real(kind=8)   :: DEPTH,THICK,VTOP,VBOT
  integer        :: vtype,ii

  call ak_init
  select case (vtype)
     case(1)
        vak=>vpak
        VTOP = 8.1
     case(2)
        vak=>vsak
        VTOP = 4.5
     case default
        stop 'vtype can only be 1 or 2'
  end select

! Linear interpolation between bottom ak135 and top regional

  DEPTH = 200!500
  THICK = 10
  VBOT  = VTOP

  do ii=1,56!26
     if (DEPTH >= depak(ii+1)) then
       VBOT = vak(ii+1)!+(vak(ii+1)-vak(ii))*(DEPTH-depak(ii))/(depak(ii+1)-depak(ii))
     endif
  enddo

  if (dep >= DEPTH) then
     vel2=VBOT
     return
  endif
  if (dep <= (DEPTH-THICK)) then
     vel2=VTOP
     return
  endif

  vel2 = VTOP+((VBOT-VTOP)*(dep-DEPTH+THICK)/THICK)

  return


! Old definition

!  if (dep >= depak(14)) then
!     vel2=vak(14)
!     return
!  endif
!  if (dep <= depak(1)) then
!     vel2=vak(1)
!     return
!  endif
!  do ii=1,135
!     if (dep <= depak(ii+1)) then
!        vel2 = vak(ii)+(vak(ii+1)-vak(ii))*(dep-depak(ii))/(depak(ii+1)-depak(ii))
!        return
!     endif
!  end do
!  stop 'vel2 error'

end function vel2




!-------------------------------------
! below this functions are not used





function vel3(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel3,dep,lat,long
  integer        :: vtype,ii

  call ak_init
  select case (vtype)
     case(1)
        vel3=3.0d0
     case(2)
        vel3=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

end function vel3

function vel4(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel4,dep,lat,long
  integer        :: vtype,ii

  call ak_init
  select case (vtype)
     case(1)
        vel4=3.0d0
     case(2)
        vel4=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

end function vel4

function vel5(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel5,dep,lat,long
  integer        :: vtype,ii

  call ak_init
  select case (vtype)
     case(1)
        vel5=3.0d0
     case(2)
        vel5=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

end function vel5

function vel6(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel6,dep,lat,long
  integer        :: vtype,ii

  call ak_init
  select case (vtype)
     case(1)
        vel6=3.0d0
     case(2)
        vel6=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

end function vel6

function vel7(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel7,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel7=3.0d0
     case(2)
        vel7=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel7

function vel8(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel8,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel8=3.0d0
     case(2)
        vel8=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel8

function vel9(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel9,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel9=3.0d0
     case(2)
        vel9=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel9

function vel10(dep,lat,long,vtype)
  use akmod
  real(kind=8)  :: vel10,dep,lat,long
  integer        :: vtype

  select case (vtype)
     case(1)
        vel10=3.0d0
     case(2)
        vel10=1.8d0
     case default
        stop 'vtype can only be 1 or 2'
  end select

  return

end function vel10

!function vel1(dep,lat,long,vtype)
!  use akmod
!  real(kind=8)  :: vel1,dep,lat,long
!  integer        :: vtype,ii
!
!  call ak_init
!  select case (vtype)
!     case(1)
!        vak=>vpak
!     case(2)
!        vak=>vsak
!     case default
!        stop 'vtype can only be 1 or 2'
!  end select
!
!
!! we have to define the velocity outside the region to be constant
!! and equal to that on the boundary, to ensure the proper velocity jump
!! at the interface. Just interpolating could destroy the discontinuity
!! if an interface node does not lie exactly on the discontinuity because
!! of rounding errors
!
!  if (dep >= depak(13)) then
!     vel1=vak(13)   ! at the bottom of region 3 and below the velocity is constant
!!!     print *, 'BELOW BOTTOM'
!     return
!  endif
!  if (dep <= depak(1)) then
!    vel1=vak(1)     ! at the top of region 3 and above the velocity is constant
!     return
!  endif
!
!  do ii=1,135  ! inbetween do linear interpolation 
!     if (dep <= depak(ii+1)) then
!        vel1 = vak(ii)+(vak(ii+1)-vak(ii))*(dep-depak(ii))/(depak(ii+1)-depak(ii))
!        return
!     endif
!  end do
!  stop 'vel1 error'
!end function vel1







!--------------------------!
!    1D REFERENCE MODEL    !
!--------------------------!

!  do ii=1,135
!     if (dep <= depak(ii+1)) then
!        vel1 = vak(ii)+(vak(ii+1)-vak(ii))*(dep-depak(ii))/(depak(ii+1)-depak(ii))
!        return
!     endif
!  end do
!  stop 'vel1 error'

!  select case (vtype)
!     case(1)
!        vel1=5.8d0!6.2d0!5.8d0
!        if (dep .GT. deppth1) then
!            vel1=8.5d0!6.8d0!8.5d0
!        endif
!        if (dep .GT. deppth2) then
!            vel1=5.8d0!7.6d0!9.0d0
!        endif
!        if (dep .GT. deppth3) then
!            vel1=9.5d0
!        endif
!        if (dep .GT. deppth4) then
!           vel1=8.5d0!9.0d0
!        endif
!     case(2)
!        vel1=3.5d0!3.6d0!3.5d0
!        if (dep .GT. deppth1) then
!            vel1=4.8d0!3.8d0!4.8d0
!        endif
!        if (dep .GT. deppth2) then
!            vel1=3.5d0!4.2d0!4.9d0
!        endif
!        if (dep .GT. deppth3) then
!            vel1=5.2d0
!        endif
!        if (dep .GT. deppth4) then
!            vel1=4.8d0!4.9d0
!        endif
!  end select


!function vel2(dep,lat,long,vtype,j,k)
!  use akmod
!  real(kind=8)  :: vel1,dep,lat,long
!  integer        :: vtype
!
!  select case (vtype)
!     case(1)
!        vel2=7.0d0
!     case(2)
!        vel2=6.0d0
!     case default
!        stop 'vtype can only be 1 or 2'
!  end select
!
!  return
!
!end function vel2
