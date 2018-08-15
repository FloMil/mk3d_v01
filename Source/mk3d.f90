      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!                                                         !!
      !!                                                         !!
      !!      MULTI-MODE 3D KIRCHHOFF MIGRATION OF RF DATA       !!
      !!                                                         !!
      !!                                                         !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!                            !                            !!
      !!      Florian Millet        !         UCB Lyon 1         !!
      !!                            !         Ui Bergen          !!
      !!                            !                       2018 !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM mk3d
      IMPLICIT NONE


                   !*******************************!
                   !*                             *!
!==================!*    Define the variables     *!===================!
                   !*                             *!
                   !*******************************!

!       Parameters
      INCLUDE 'cst_param.finc'
!       Variables
      INTEGER :: nsrc, nrec                 ! Number of src & rec
      INTEGER :: lentt,lenwf                ! Length of tt & wf files
      INTEGER :: srat, Pdel, mute           ! wf sampling rate & delay
      REAL    :: sfrq                       ! wf sampling frequency
      REAL    :: tmig                       ! ts+tp2-tp1 for grid
      INTEGER :: twf                        ! ts+tp2-tp1 for RF
      REAL    :: tppp,tpps,tpss             ! IDEM multiples for grid
      INTEGER :: twf2,twf3,twf4             ! IDEM multiples for RF
      REAL    :: maxT                       ! Duration of RF
      INTEGER :: maxLat,maxLon,maxDep       ! Migration grid maxima
      REAL    :: LatStep,LonStep,DepStep    ! Lat, Lon & Dep resolution
      REAL    :: rlat,rlon                  ! for bilinear interpolation
      REAL    :: val1,val2                  ! Bilinear interpolation
      REAL    :: dist,sprd                  ! Geometrical spreading
      REAL    :: WF,wfR,wfT,wfZ,NT          ! Tot. & 3C RF and nthroot
      COMPLEX :: PH,phR,phT,phZ,phC         ! Tot. & 3C Instant. phase
      REAL    :: nval                       ! 'N' of Nth root
      REAL    :: teta_dps                   ! arr-dep & dep-rad angles
      REAL    :: teta_spp,teta_sps          ! arr-dep (surface) angles
      REAL    :: teta_mpp,teta_mps          ! arr-dep (multiples) angles
      REAL    :: teta_msp,teta_mss          ! arr-dep (multiples) angles
      REAL    :: PuS_b                      ! P-to-S / Vs scat pat
      REAL    :: PuP_r,PuS_r                ! P-to-(p,s) / rho scat pat
      REAL    :: PuP_s,PuS_s                ! P-to-(p,s) / beta scat pat
      REAL    :: PuP_2,PuS_2,SuS_1,SuS_2    ! (p,s)to(P,S) / Vs scat pat
      REAL    :: spr                        ! Vs/Vp ratio
      REAL    :: svc,shc                    ! Ps(S) Sv and Sh components
      REAL    :: weight                     ! Total weight of each wf
      INTEGER :: pow                        ! Focus power (cos^pow)
      REAL    :: PhRa                       ! +- Phase
      REAL(kind=8)      :: renorm           ! Vector renormalization
      REAL,DIMENSION(3) :: AmpCor           ! Amplitude correction (mig)
      REAL,DIMENSION(5) :: focus            ! Focus energy under station
      REAL,DIMENSION(5) :: phaw             ! Phases' weights
      INTEGER,DIMENSION(8) :: rn            ! QC renorm
!       Arrays
      INTEGER,DIMENSION(:), ALLOCATABLE :: grdinfo    ! Sizes box/grid
      REAL,DIMENSION(:,:),  ALLOCATABLE :: tp1        ! First P arrival
      REAL,DIMENSION(:,:),  ALLOCATABLE :: RECinfo    ! Receivers info
      REAL,DIMENSION(:,:),  ALLOCATABLE :: SRCinfo    ! Sources info
      INTEGER,DIMENSION(:,:), ALLOCATABLE :: Corr     ! Src-Rec-Wf corr
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: wfrms    ! Waveforms
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: Emig     ! Migrated energy
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: Fmig     ! Migratd nrj PWS
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: Nmig     ! Migrated Nth rt
      COMPLEX,DIMENSION(:,:,:,:),ALLOCATABLE :: ansig ! Analytical Sig.
      COMPLEX,DIMENSION(:,:,:,:),ALLOCATABLE :: Pmig  ! Migrated phase
      REAL,DIMENSION(36)                  :: Emax     ! Max energy/ph
      INTEGER,DIMENSION(:,:),ALLOCATABLE  :: CN       ! Closest neigh'
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: Ptt      ! P  travel times
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: Stt      ! xS travel times
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: Att      ! Pp travel times
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: Btt      ! Ps travel times
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: Ctt      ! xP travel times
      REAL,DIMENSION(:,:),ALLOCATABLE   :: bsip       ! baz,slw,inc,pol
      REAL(kind=8),DIMENSION(3) :: arvl,dptr,rDir,tDir,vDir! ar&dp vctrs
      REAL(kind=8),DIMENSION(3) :: fsri,fsro          ! free srf vctrs
      REAL(kind=8),DIMENSION(3) :: mulp,muls,finp     ! mltpls   vctrs
      REAL,DIMENSION(3) :: gama_ps,beta_ps            ! Scat Pat vectors
      REAL,DIMENSION(3) :: gama_pps,beta_pps          ! Scat Pat vectors
      REAL,DIMENSION(3) :: gama_psx,beta_psx          ! Scat Pat vectors
      REAL,DIMENSION(3) :: gama_pss,beta_pss          ! Scat Pat vectors
      REAL,DIMENSION(3) :: delta_psx,delta_pss        ! Scat Pat vectors
      REAL,DIMENSION(:,:),ALLOCATABLE :: Projection   ! (P,S) > (R,T,Z)
!       File names
      CHARACTER(LEN=30) :: infile           ! Input file name
      CHARACTER(LEN=30) :: outfile          ! Output file name
      CHARACTER(LEN=30) :: Tname            ! Traveltime files
      CHARACTER(LEN=60) :: Tfile            ! Traveltime files FFORT
      CHARACTER(LEN=30),DIMENSION(4) :: WFname   ! Waveform file
      CHARACTER(LEN=60),DIMENSION(4) :: WFfile   ! Waveform file FFORT
      CHARACTER(LEN=30) :: Infoname         ! Information file
      CHARACTER(LEN=60) :: Infofile         ! Information file FFORT
      CHARACTER(LEN=30) :: COname           ! src-rec correlation file
      CHARACTER(LEN=60) :: COfile           ! src-rec correlation FFORT
!       Counters, flags, dummies and timing
      INTEGER :: i, j, k, l                 ! Counters
      INTEGER :: ii,jj,kk,ll                ! Counters
      INTEGER :: lP,lS                      ! Correspondance counters
      INTEGER :: dum                        ! Integer dummy
      REAL    :: dumm                       ! Real dummy
      INTEGER :: recount                    ! Records counter (emax)
      REAL    :: ranu                       ! Random Number
      INTEGER :: t0,t1,t2,t3,t4,rate        ! Timing routine
      INTEGER :: nbhrs,nbmns                ! Timings
      INTEGER :: splmnt                     ! File number indicative
      splmnt  = 0000


                    !******************************!
                    !*                            *!
!===================!*    Read the input files    *!===================!
                    !*                            *!
                    !******************************!

!       Write input file name
      infile = '../Input/mk3d.in'      ! Name of input file
      ALLOCATE(grdinfo(15))
      ALLOCATE(CN(4,2))

!-----------------------------------------------------------!
!     Read the input parameters for selection of traces     !
!-----------------------------------------------------------!
      OPEN(101,file=infile,form="formatted")
      READ(101,*) !#
      READ(101,*) !#
      READ(101,*) !#
      READ(101,*) ! Waveforms
      READ(101,*) WFname(1)
      READ(101,*) WFname(2)
      READ(101,*) WFname(3)
      READ(101,*) WFname(4)
      READ(101,*) ! Src-Rec
      READ(101,*) Infoname
      READ(101,*) ! Corresp
      READ(101,*) COname
      READ(101,*) !#
      READ(101,*) !#
      READ(101,*) !#
      READ(101,*) ! Weights of phases
      READ(101,*) phaw(1),phaw(2),phaw(3),phaw(4),phaw(5)
      READ(101,*) !#
      READ(101,*) !#
      READ(101,*) !#
      READ(101,*) ! Size of boxes
      READ(101,*) grdinfo(1)
      READ(101,*) grdinfo(2)
      READ(101,*) grdinfo(3) 
      READ(101,*) ! Size of boxes
      READ(101,*) grdinfo(13)
      READ(101,*) grdinfo(14)
      READ(101,*) grdinfo(15) 
      READ(101,*) ! Size of propagation grid
      READ(101,*) grdinfo(4) 
      READ(101,*) grdinfo(5) 
      READ(101,*) grdinfo(6) 
      READ(101,*) ! Size of migration grid
      READ(101,*) grdinfo(7) 
      READ(101,*) grdinfo(8) 
      READ(101,*) grdinfo(9) 
      READ(101,*) ! Origin of migration grid
      READ(101,*) grdinfo(10) 
      READ(101,*) grdinfo(11) 
      READ(101,*) grdinfo(12) 
      CLOSE(101)
      lentt   = grdinfo(4) * grdinfo(5) * grdinfo(6)
      maxDep  = grdinfo(7) + grdinfo(10)
      maxLat  = grdinfo(8) + grdinfo(11)
      maxLon  = grdinfo(9) + grdinfo(12)
      DepStep = REAL(grdinfo(1)) / REAL((grdinfo(4)-1))
      LatStep = REAL(grdinfo(2)) / REAL((grdinfo(5)-1))
      LonStep = REAL(grdinfo(3)) / REAL((grdinfo(6)-1))


              !******************************************!
              !*                                        *!
!=============!*    Read the preliminary information    *!=============!
              !*                                        *!
              !******************************************!

!------------------------------------------!
!     Prepare all the files and arrays     !
!------------------------------------------!
!       Write input files names
      Infofile  = '../Data/' // Infoname    ! Src-Rec info
      WFfile(1) = '../Data/' // WFname(1)   ! Waveform
      WFfile(2) = '../Data/' // WFname(2)   ! Waveform
      WFfile(3) = '../Data/' // WFname(3)   ! Waveform
      WFfile(4) = '../Data/' // WFname(4)   ! Waveform

!       Open the files
      OPEN(113,file=Infofile, form="formatted")
      OPEN(114,file=WFfile(1),form="formatted")

!------------------------------!
!     Read the information     !
!------------------------------!
!       Number of Sources and Receivers
      READ(113,*) nsrc
      DO i=1,nsrc
        READ(113,*)
      ENDDO
      READ(113,*) nrec
      CLOSE(113)
      OPEN(113,file=Infofile,form="formatted")

!       Waveform information
      READ(114,*)
      READ(114,*) dumm, lenwf, sfrq, dumm, dumm
      srat=INT(1/sfrq)
      Pdel=INT(dumm)
      CLOSE(114)

!       Allocate memory
      ALLOCATE(RECinfo(nrec,4),SRCinfo(nsrc,3))
      ALLOCATE(Ptt(nsrc,grdinfo(6),grdinfo(5),  grdinfo(4)))
      ALLOCATE(Stt(nrec,grdinfo(6),grdinfo(5),  grdinfo(4)))
      ALLOCATE(Att(nsrc,grdinfo(6),grdinfo(5),  grdinfo(4)))
      ALLOCATE(Btt(nsrc,grdinfo(6),grdinfo(5),  grdinfo(4)))
      ALLOCATE(Ctt(nrec,grdinfo(6),grdinfo(5),  grdinfo(4)))
      ALLOCATE(Emig(7,grdinfo(9)+1,grdinfo(8)+1,grdinfo(7)+1))
      ALLOCATE(Pmig(7,grdinfo(9)+1,grdinfo(8)+1,grdinfo(7)+1))
      ALLOCATE(Fmig(7,grdinfo(9)+1,grdinfo(8)+1,grdinfo(7)+1))
      ALLOCATE(Nmig(7,grdinfo(9)+1,grdinfo(8)+1,grdinfo(7)+1))
      ALLOCATE(wfrms(nsrc*nrec,lenwf,3,4))
      ALLOCATE(ansig(nsrc*nrec,lenwf,3,4))
      ALLOCATE(Corr(nsrc*nrec,8))
      ALLOCATE(tp1(nrec,nsrc))
      ALLOCATE(bsip(nrec*nsrc,4))
      ALLOCATE(Projection(5,3))

!        Talk to people :)
      WRITE(*,12) ''
      WRITE(*,12) '+-------------------------------+'
      WRITE(*,12) '|                               |'
      WRITE(*,11) '|   Number of sources   : ',nsrc,'   |'
      WRITE(*,11) '|   Number of receivers : ',nrec,'   |'
      WRITE(*,12) '|                               |'
      WRITE(*,12) '+-------------------------------+'
      WRITE(*,12) ''
11    FORMAT(A46,I3,A4)
12    FORMAT(A53)


                 !***********************************!
                 !*                                 *!
!================!*    Read the wfields & wforms    *!=================!
                 !*                                 *!
                 !***********************************!

!----------------------------------------!
!     Put all the data in the memory     !
!----------------------------------------!
      CALL system_clock(t0,rate)
      WRITE(*,*) 
      WRITE(*,*) '  PUT DATA IN MEMORY'

!       Read Source and Receiver information
      READ(113,*)
      DO i=1,nsrc
        READ(113,*) SRCinfo(i,1),SRCinfo(i,2),SRCinfo(i,3)
      ENDDO
      READ(113,*)
      DO i=1,nrec
        READ(113,*) RECinfo(i,1),RECinfo(i,2),RECinfo(i,3),RECinfo(i,4)
      ENDDO
      CLOSE(113)
      Ptt = -100
      Stt = -100

!       Work with the wavefields
      WRITE(*,*)
      WRITE(*,*) '    Sources (P-waves)'
      DO l=1,nsrc                                   ! Nb src
        !IF (l .NE. 03) CYCLE                       ! Only one source
        ll = l + splmnt
        WRITE(Tfile,13) 'fm3d/data/srcs/s', ll, ".dat"  ! Wvflds
        OPEN(111,file=Tfile,form="formatted")
!         Direct P wave
        DO i=1,5
          READ(111,*)                               ! Header
        ENDDO
        DO i=1,grdinfo(6)                           ! Longitude
          DO j=1,grdinfo(5)                         ! Latitude
            DO k=1,grdinfo(4)                       ! Depth
              READ(111,*) Ptt(l,i,j,k)
            ENDDO
          ENDDO
        ENDDO
!         Pp wave
        READ(111,*)                                 ! Header
        DO i=1,grdinfo(6)                           ! Longitude
          DO j=1,grdinfo(5)                         ! Latitude
            DO k=1,grdinfo(4)                       ! Depth
              READ(111,*) Att(l,i,j,k)
            ENDDO
          ENDDO
        ENDDO
!         Ps wave
        READ(111,*)                                 ! Header
        DO i=1,grdinfo(6)                           ! Longitude
          DO j=1,grdinfo(5)                         ! Latitude
            DO k=1,grdinfo(4)                       ! Depth
              READ(111,*) Btt(l,i,j,k)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(111)
        IF (MOD(l,10).EQ.1) WRITE(*,"(A5)",advance='no') ''
        WRITE(*,"(i6)",advance='no') l
        IF (MOD(l,10).EQ.0) WRITE(*,*)
      ENDDO

      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) '    Receivers (S-waves)'
      DO l=1,nrec                                   ! Nb rec
        !IF (l .NE. 5) CYCLE                        ! Only one receiver
        ll = l + splmnt
        WRITE(Tfile,13) 'fm3d/data/recs/r', ll, ".dat"  ! Wvflds
        OPEN(111,file=Tfile,form="formatted")
!         xS wave
        DO i=1,5
          READ(111,*)                               ! Header
        ENDDO
        DO i=1,grdinfo(6)                           ! Longitude
          DO j=1,grdinfo(5)                         ! Latitude
            DO k=1,grdinfo(4)                       ! Depth
              READ(111,*) Stt(l,i,j,k)
            ENDDO
          ENDDO
        ENDDO
!         xP wave
        READ(111,*)                                 ! Header
        DO i=1,grdinfo(6)                           ! Longitude
          DO j=1,grdinfo(5)                         ! Latitude
            DO k=1,grdinfo(4)                       ! Depth
              READ(111,*) Ctt(l,i,j,k)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(111)
        IF (MOD(l,10).EQ.1) WRITE(*,"(A5)",advance='no') ''
        WRITE(*,"(i6)",advance='no') l
        IF (MOD(l,10).EQ.0) WRITE(*,*)
      ENDDO
13    FORMAT(A16,I4.4,A4)

!       Correspondance Src-Rec-Wf
      WRITE(*,*) 
      WRITE(*,*) 
      WRITE(*,*) '    Quality Control...'
      WRITE(*,*) 
      COfile = '../Data/' // COname     ! Src-Rec info
      OPEN(121,file=COfile,   form="formatted")
      DO i=1,nsrc*nrec
        READ(121,*) Corr(i,1), Corr(i,2), Corr(i,3), Corr(i,4), &
                    Corr(i,5), Corr(i,6), Corr(i,7), Corr(i,8)
      ENDDO
      CLOSE(121)

!       Renorm given QC values for final stacking
      rn(4) = SUM(Corr(:,4))                ! PS
      rn(5) = SUM(Corr(:,5))                ! PpP
      rn(6) = SUM(Corr(:,6))                ! PpS
      rn(7) = SUM(Corr(:,7))                ! PsSv
      rn(8) = SUM(Corr(:,8))                ! PsSh
      WRITE(*,15) 'PS  :',rn(4)
      WRITE(*,15) 'PpP :',rn(5)
      WRITE(*,15) 'PpS :',rn(6)
      WRITE(*,15) 'PsSv:',rn(7)
      WRITE(*,15) 'PsSh:',rn(8)
15    FORMAT(A12,I5)

!         Work with the waveforms
      WRITE(*,*) 
      WRITE(*,*) 
      WRITE(*,*) '    Waveforms...'
      WRITE(*,*) 
      DO j=1,4
        IF (j.EQ.1) dum = 114
        IF (j.EQ.2) dum = 115
        IF (j.EQ.3) dum = 116
        IF (j.EQ.4) dum = 117
        OPEN(dum,file=WFfile(j),form="formatted")
        READ(dum,*)                         ! File header 
        READ(dum,*)
        DO l=1,nsrc*nrec
          READ(dum,*)                       ! Trace header 
          READ(dum,*)
          READ(dum,*)
          IF (Corr(l,3).NE.0) THEN
            DO i=1,lenwf
              READ(dum,*) wfrms(l,i,1,j),wfrms(l,i,2,j),wfrms(l,i,3,j),&
                          ansig(l,i,1,j),ansig(l,i,2,j),ansig(l,i,3,j)
              IF (ISNAN(wfrms(l,i,1,j))) wfrms(l,i,1,j) = 0
              IF (ISNAN(wfrms(l,i,2,j))) wfrms(l,i,2,j) = 0
              IF (ISNAN(wfrms(l,i,3,j))) wfrms(l,i,3,j) = 0
              IF (ISNAN(REAL(ansig(l,i,1,j)))) ansig(l,i,1,j) = 0
              IF (ISNAN(REAL(ansig(l,i,2,j)))) ansig(l,i,2,j) = 0
              IF (ISNAN(REAL(ansig(l,i,3,j)))) ansig(l,i,3,j) = 0
            ENDDO
          ENDIF
        ENDDO
        CLOSE(dum)
      ENDDO

      CALL system_clock(t1,rate)
      WRITE(*,*) 
      WRITE(*,*) 
      WRITE(*,*) 'Time to load data :      ',(t1-t0)/rate,' seconds'

!       Talk to people :)
      WRITE(*,*) 
      WRITE(*,*) 
      WRITE(*,*) '    All timefields read'

               !***************************************!
               !*                                     *!
!==============!*    Start working with the traces    *!===============!
               !*                                     *!
               !***************************************!

!------------------------------------!
!     Find tp1 (first P arrival)     !
!------------------------------------!
      DO i=1,nsrc                           ! Number of sources
        DO j=1,nrec                         ! Number of receivers

!           Find closest neighbours
          rlat = (RECinfo(j,2)-grdinfo(14))/LatStep+1
          rlon = (RECinfo(j,3)-grdinfo(15))/LonStep+1
          lP = Corr((i-1)*nrec+j,1)
          lS = Corr((i-1)*nrec+j,2)
          CN(1,:) = (/ FLOOR(rlat),  FLOOR(rlon) /)
          CN(2,:) = (/ FLOOR(rlat),  FLOOR(rlon)+1 /)
          CN(3,:) = (/ FLOOR(rlat)+1,FLOOR(rlon)+1 /)
          CN(4,:) = (/ FLOOR(rlat)+1,FLOOR(rlon) /)

!           Bilinear interpolation
          val1 = (CN(2,2)-rlon)*Ptt(i,CN(1,2),CN(1,1),maxDep) &
               + (rlon-CN(1,2))*Ptt(i,CN(2,2),CN(2,1),maxDep)
          val2 = (CN(3,2)-rlon)*Ptt(i,CN(4,2),CN(4,1),maxDep) &
               + (rlon-CN(4,2))*Ptt(i,CN(3,2),CN(3,1),maxDep)
          tp1(lS,lP) = (CN(4,1)-rlat) * val1 &
                     + (rlat-CN(1,1)) * val2
        ENDDO
      ENDDO

!-----------------------------------------!
!     Initialize the values and stuff     !
!-----------------------------------------!
!       Initialize Emig
      Emig = ZERO
      Emax = ZERO
      Pmig = ZERO
      Fmig = ZERO
      Nmig = ZERO

!       Find max time of RF
      maxT=(lenwf/srat)-Pdel

!       Grab BAZ and SLWN information
      OPEN(102,file='../Run/raysum/geo',form='formatted')
      READ(102,*)
      DO i=1,nsrc*nrec
        READ(102,*) bsip(i,1), bsip(i,2)
        bsip(i,1) = bsip(i,1) * pi / 180
        bsip(i,3) = ASIN(bsip(i,2)*5800)! * 180 / pi
      ENDDO
      CLOSE(102)
!       Define Vp/Vs ratio
      spr = 0.54 

!----------------------------!
!     Migrate the traces     !
!----------------------------!
      WRITE(*,*) 
      WRITE(*,*) 
      WRITE(*,*) '  MIGRATE THE TRACES'
      recount = 0
      DO l=1,nsrc*nrec
        lP = Corr(l,1)
        lS = Corr(l,2)
        !IF  (lP .NE. 01) CYCLE
        !IF  (lS .NE. 01) CYCLE
        !IF ((lP .NE. 01) .AND. ( lP .NE. 01 )) CYCLE

!     Write the progression of the migration
        IF (MOD(l,10).EQ.1) WRITE(*,"(A5)",advance='no') ''
        IF (Corr(l,3).NE.0) WRITE(*,"(i6)",advance='no') l
        IF (Corr(l,3).EQ.0) WRITE(*,"(A6)",advance='no') 'X'
        IF (MOD(l,10).EQ.0) WRITE(*,*)
        IF (Corr(l,3).EQ.0) CYCLE
        recount = recount + 1
        CALL random_number(PhRa)
        PhRa = SIGN(1.,(PhRa-.5))

!     Start migration loops
        rlat = (RECinfo(lS,2)-grdinfo(14))/LatStep+1
        rlon = (RECinfo(lS,3)-grdinfo(15))/LonStep+1
        DO i=grdinfo(12),maxLon             ! Longitude
          ii = i - grdinfo(12) + 1
          DO j=grdinfo(11),maxLat           ! Latitude
            jj = j - grdinfo(11) + 1
            DO k=grdinfo(10),maxDep         ! Depth
              kk = k - grdinfo(10) + 1

!     Skip the edges of the box (gradient)
              IF ((ii.EQ.1).OR.(i.EQ.maxLon).OR. &
                  (jj.EQ.1).OR.(j.EQ.maxLat).OR. &
                  (kk.EQ.1).OR.(k.EQ.maxDep)) THEN
                CYCLE
              ENDIF

!     Find migration time and associate waveform
              tmig = Stt(lS,i,j,k) + Ptt(lP,i,j,k) - tp1(lS,lP)
              twf  = INT((tmig + Pdel) * srat)     ! PS conversion
              tppp = Ctt(lS,i,j,k) + Att(lP,i,j,k) - tp1(lS,lP)
              twf2 = INT((tppp + Pdel) * srat)     ! PpP multiple
              tpps = Stt(lS,i,j,k) + Att(lP,i,j,k) - tp1(lS,lP)
              twf3 = INT((tpps + Pdel) * srat)     ! PpS multiple
              tpss = Stt(lS,i,j,k) + Btt(lP,i,j,k) - tp1(lS,lP)
              twf4 = INT((tpss + Pdel) * srat)     ! PsS multiple

!     Correct for geometrical spreading
              dist = SQRT(((j-rlat) * LatStep * dkm)**2 &
                           + ((i-rlon) * LonStep * dkm)**2 &
                             + ((maxDep-k) * DepStep)**2)
              sprd = 1 / (MAX(dist,30.)**(1.0))
              !sprd = 1 / MAX(tmig,10.)                     

!     Find arrival-to-departure angle
!       Direct P wave
              arvl(1) = (Ptt(lP,i,j+1,k) - Ptt(lP,i,j-1,k) ) / 2 !+-
              arvl(2) = (Ptt(lP,i+1,j,k) - Ptt(lP,i-1,j,k) ) / 2 !-+
              arvl(3) = (Ptt(lP,i,j,k-1) - Ptt(lP,i,j,k+1) ) / 2 !-+
              renorm = arvl(1)**2 + arvl(2)**2 + arvl(3)**2
              arvl = arvl / DSQRT((renorm))
!       Downgoing p wave
              mulp(1) = (Att(lP,i,j+1,k) - Att(lP,i,j-1,k) ) / 2
              mulp(2) = (Att(lP,i+1,j,k) - Att(lP,i-1,j,k) ) / 2
              mulp(3) = (Att(lP,i,j,k-1) - Att(lP,i,j,k+1) ) / 2
              renorm = mulp(1)**2 + mulp(2)**2 + mulp(3)**2
              mulp = mulp / DSQRT((renorm))
!       Downgoing s wave
              muls(1) = (Btt(lP,i,j+1,k) - Btt(lP,i,j-1,k) ) / 2
              muls(2) = (Btt(lP,i+1,j,k) - Btt(lP,i-1,j,k) ) / 2
              muls(3) = (Btt(lP,i,j,k-1) - Btt(lP,i,j,k+1) ) / 2
              renorm = muls(1)**2 + muls(2)**2 + muls(3)**2
              muls = muls / DSQRT((renorm))
!       Final S wave
              dptr(1) = (Stt(lS,i,j+1,k) - Stt(lS,i,j-1,k) ) / 2 !-+
              dptr(2) = (Stt(lS,i+1,j,k) - Stt(lS,i-1,j,k) ) / 2 !-+
              dptr(3) = (Stt(lS,i,j,k-1) - Stt(lS,i,j,k+1) ) / 2 !+-
              renorm = dptr(1)**2 + dptr(2)**2 + dptr(3)**2
              dptr = dptr / DSQRT((renorm))
!       Final P wave
              finp(1) = (Ctt(lS,i,j+1,k) - Ctt(lS,i,j-1,k) ) / 2
              finp(2) = (Ctt(lS,i+1,j,k) - Ctt(lS,i-1,j,k) ) / 2
              finp(3) = (Ctt(lS,i,j,k-1) - Ctt(lS,i,j,k+1) ) / 2
              renorm = finp(1)**2 + finp(2)**2 + finp(3)**2
              finp = finp / DSQRT((renorm))
!       Uniform background model for free surface reflection
              fsri(1) = -COS(bsip(l,1)) * SIN(bsip(l,3))
              fsri(2) = -SIN(bsip(l,1)) * SIN(bsip(l,3))
              fsri(3) = -COS(bsip(l,3))
              fsro(1) =  COS(bsip(l,1)) * SIN(bsip(l,3))
              fsro(2) =  SIN(bsip(l,1)) * SIN(bsip(l,3))
              fsro(3) = -COS(bsip(l,3))
!       Compute the angles
              teta_dps = ACOS(REAL(DOT_PRODUCT(arvl,dptr)))
              teta_spp = ACOS(REAL(DOT_PRODUCT(fsri,fsro)))
              teta_sps = ACOS(REAL(DOT_PRODUCT(fsri,fsro)))
              teta_mps = ACOS(REAL(DOT_PRODUCT(mulp,dptr)))
              teta_mpp = ACOS(REAL(DOT_PRODUCT(mulp,finp)))
              teta_mss = ACOS(REAL(DOT_PRODUCT(muls,dptr)))

!     Apply 3-D scattering patterns (Beylkin & Burridge 1990)
              ! gama_ps = cross_product( arvl , dptr )       !-+
              gama_ps(1) = arvl(2)*dptr(3) - arvl(3)*dptr(2) !-+
              gama_ps(2) = arvl(3)*dptr(1) - arvl(1)*dptr(3) !-+
              gama_ps(3) = arvl(1)*dptr(2) - arvl(2)*dptr(1) !-+
              renorm = gama_ps(1)**2 + gama_ps(2)**2 + gama_ps(3)**2
              gama_ps = gama_ps / DSQRT((renorm))
              ! beta_ps = cross_product( dptr , gama_ps )           !++
              beta_ps(1) = dptr(2)*gama_ps(3) - dptr(3)*gama_ps(2) !+-
              beta_ps(2) = dptr(3)*gama_ps(1) - dptr(1)*gama_ps(3) !+-
              beta_ps(3) = dptr(1)*gama_ps(2) - dptr(2)*gama_ps(1) !+-
              renorm = beta_ps(1)**2 + beta_ps(2)**2 + beta_ps(3)**2
              beta_ps = beta_ps / DSQRT((renorm))

              ! gama_pps
              gama_pps(1) = mulp(2)*dptr(3) - mulp(3)*dptr(2)
              gama_pps(2) = mulp(3)*dptr(1) - mulp(1)*dptr(3)
              gama_pps(3) = mulp(1)*dptr(2) - mulp(2)*dptr(1)
              renorm = gama_pps(1)**2 + gama_pps(2)**2 + gama_pps(3)**2
              gama_pps = gama_pps / DSQRT((renorm))
              ! beta_pps
              beta_pps(1) = dptr(2)*gama_pps(3) - dptr(3)*gama_pps(2)
              beta_pps(2) = dptr(3)*gama_pps(1) - dptr(1)*gama_pps(3)
              beta_pps(3) = dptr(1)*gama_pps(2) - dptr(2)*gama_pps(1)
              renorm = beta_pps(1)**2 + beta_pps(2)**2 + beta_pps(3)**2
              beta_pps = beta_pps / DSQRT((renorm))

              ! gama_psx
              gama_psx(1) = fsri(2)*fsro(3) - fsri(3)*fsro(2)
              gama_psx(2) = fsri(3)*fsro(1) - fsri(1)*fsro(3)
              gama_psx(3) = fsri(1)*fsro(2) - fsri(2)*fsro(1)
              renorm = gama_psx(1)**2 + gama_psx(2)**2 + gama_psx(3)**2
              gama_psx = gama_psx / DSQRT((renorm))
              ! beta_psx
              beta_psx(1) = fsro(2)*gama_psx(3) - fsro(3)*gama_psx(2)
              beta_psx(2) = fsro(3)*gama_psx(1) - fsro(1)*gama_psx(3)
              beta_psx(3) = fsro(1)*gama_psx(2) - fsro(2)*gama_psx(1)
              renorm = beta_psx(1)**2 + beta_psx(2)**2 + beta_psx(3)**2
              beta_psx = beta_psx / DSQRT((renorm))
              ! delta_psx
              delta_psx(1) = -muls(2)*gama_psx(3) + muls(3)*gama_psx(2)
              delta_psx(2) = -muls(3)*gama_psx(1) + muls(1)*gama_psx(3)
              delta_psx(3) = -muls(1)*gama_psx(2) + muls(2)*gama_psx(1)
              renorm = delta_psx(1)**2+delta_psx(2)**2+delta_psx(3)**2
              delta_psx = delta_psx / DSQRT((renorm))

              ! gama_pss
              gama_pss(1) = muls(2)*dptr(3) - muls(3)*dptr(2)
              gama_pss(2) = muls(3)*dptr(1) - muls(1)*dptr(3)
              gama_pss(3) = muls(1)*dptr(2) - muls(2)*dptr(1)
              renorm = gama_pss(1)**2 + gama_pss(2)**2 + gama_pss(3)**2
              gama_pss = gama_pss / DSQRT((renorm))
              ! beta_pss
              beta_pss(1) = dptr(2)*gama_pss(3) - dptr(3)*gama_pss(2)
              beta_pss(2) = dptr(3)*gama_pss(1) - dptr(1)*gama_pss(3)
              beta_pss(3) = dptr(1)*gama_pss(2) - dptr(2)*gama_pss(1)
              renorm = beta_pss(1)**2 + beta_pss(2)**2 + beta_pss(3)**2
              beta_pss = beta_pss / DSQRT((renorm))
              ! delta_pss
              delta_pss(1) = muls(2)*gama_pss(3) - muls(3)*gama_pss(2)
              delta_pss(2) = muls(3)*gama_pss(1) - muls(1)*gama_pss(3)
              delta_pss(3) = muls(1)*gama_pss(2) - muls(2)*gama_pss(1)
              renorm = delta_pss(1)**2+delta_pss(2)**2+delta_pss(3)**2
              delta_pss = delta_pss / DSQRT((renorm))

              ! Sv and Sh component of Ps(S)
              svc = DOT_PRODUCT(delta_psx,delta_pss)
              shc = DOT_PRODUCT(delta_psx,gama_pss)

              ! P-to-S forward scattering pattern with vs/beta
              PuS_b = SIN(2*teta_dps)
              PuS_b = -PuS_b ! Decrease upwards

              ! P-to-P at free surface with vs/beta
              PuP_s = 1                   ! vp/alpha
              PuP_s = COS(2*teta_spp) - 1
              PuP_s = -PuP_s ! Decrease upwards

              ! P-to-S at free surface with vs/beta
              PuS_s = SIN(2*teta_sps)
              PuS_s = -PuS_s ! Decrease upwards

              ! P-to-P back scattering with vs/beta
              PuP_2 = 1                   ! vp/alpha
              PuP_2 = COS(2*teta_mpp) - 1

              ! P-to-S back scattering with vs/beta
              PuS_2 = SIN(2*teta_mps)

              ! S-to-S back scattering with vs/beta (Sv-Sv)
              SuS_1 = COS(2*teta_mss)
              ! S-to-S back scattering with vs/beta (Sh-Sh)
              SuS_2 = COS(teta_mss)

!     Apply surface geometrical weighting (R, T and Z rfs)
              rDir(1) = -COS(bsip(l,1))
              rDir(2) = -SIN(bsip(l,1))
              rDir(3) = 0

              tDir(1) =  SIN(bsip(l,1))
              tDir(2) = -COS(bsip(l,1))
              tDir(3) = 0

              vDir(1) = 0
              vDir(2) = 0
              vDir(3) = -1

              Projection(1,1) = DOT_PRODUCT(rDir,beta_ps)
              Projection(1,2) = DOT_PRODUCT(tDir,beta_ps)
              Projection(1,3) = DOT_PRODUCT(vDir,beta_ps)
              Projection(2,1) = DOT_PRODUCT(rDir,-finp)
              Projection(2,2) = DOT_PRODUCT(tDir,-finp)
              Projection(2,3) = DOT_PRODUCT(vDir,-finp)
              Projection(3,1) = DOT_PRODUCT(rDir,beta_pps)
              Projection(3,2) = DOT_PRODUCT(tDir,beta_pps)
              Projection(3,3) = DOT_PRODUCT(vDir,beta_pps)
              Projection(4,1) = DOT_PRODUCT(rDir,beta_pss)
              Projection(4,2) = DOT_PRODUCT(tDir,beta_pss)
              Projection(4,3) = DOT_PRODUCT(vDir,beta_pss)
              Projection(5,1) = DOT_PRODUCT(rDir,gama_pss)
              Projection(5,2) = DOT_PRODUCT(tDir,gama_pss)
              Projection(5,3) = DOT_PRODUCT(vDir,gama_pss)

!     Combine all the amplitude corrections
              AmpCor = 1                    ! Initialize
              WF     = 0                    ! Initialize
              pow    = 4                    ! Focus power (see next)
              nval   = 2                    ! Nth root stack
              mute   = 3                    ! Hide initial pulse (sec)

!     Focus energy under station (Similar to Cheng)
              focus(1) = ABS(DOT_PRODUCT(-vDir,dptr)**pow)
              focus(2) = ABS(DOT_PRODUCT(-vDir,finp)**pow)
              focus(3) = ABS(DOT_PRODUCT(-vDir,dptr)**pow)
              focus(4) = ABS(DOT_PRODUCT(-vDir,dptr)**pow)
              focus(5) = ABS(DOT_PRODUCT(-vDir,dptr)**pow)

!     Re-initialize for tests
              sprd     = 1
              !focus    = 1

!     Migrate + Cover initial pulse
              IF ((tmig.LT.maxT-2).AND.(twf.GT.(mute+Pdel)*srat) &
                                  .AND.(Corr(l,4).NE.0))  THEN     ! PS
                AmpCor(:) = Projection(1,:) * PuS_b * sprd
                WF  = DOT_PRODUCT(wfrms(l,twf,:,1),AmpCor(:))* focus(1)
                phC = DOT_PRODUCT(ansig(l,twf,:,1),AmpCor(:))
                PH  = EXP(PhRa*yi*ATAN2(AIMAG(phC),REAL(phC)))
                  PH  = PH * focus(1)
                NT  = ABS(WF)**(1/nval) * SIGN(1.,WF)        * focus(1)
                weight = rn(4) / phaw(1) / RECinfo(lS,4)
                Emig(1,ii,jj,kk) = Emig(1,ii,jj,kk) + WF / weight 
                Pmig(1,ii,jj,kk) = Pmig(1,ii,jj,kk) + PH / weight 
                Nmig(1,ii,jj,kk) = Nmig(1,ii,jj,kk) + NT / weight 
              ELSEIF (Corr(l,4).NE.0) THEN
                CALL random_number(ranu)
                PH  = EXP(yi*1000*ranu)
                Pmig(1,ii,jj,kk) = Pmig(1,ii,jj,kk) + PH / weight
              ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              IF ((tppp.LT.maxT-2).AND.(twf2.GT.(mute+Pdel)*srat) &
                                  .AND.(Corr(l,5).NE.0))  THEN     ! PpP
                AmpCor(:) = Projection(2,:) * PuP_s * PuP_2 * sprd
                WF  = DOT_PRODUCT(wfrms(l,twf2,:,2),AmpCor(:))* focus(2)
                phC = DOT_PRODUCT(ansig(l,twf2,:,2),AmpCor(:))
                PH  = EXP(PhRa*yi*ATAN2(AIMAG(phC),REAL(phC)))
                  PH  = PH * focus(2)
                NT  = ABS(WF)**(1/nval) * SIGN(1.,WF)         * focus(2)
                weight = rn(5) / phaw(2) / RECinfo(lS,4)
                Emig(2,ii,jj,kk) = Emig(2,ii,jj,kk) + WF / weight
                Pmig(2,ii,jj,kk) = Pmig(2,ii,jj,kk) + PH / weight
                Nmig(2,ii,jj,kk) = Nmig(2,ii,jj,kk) + NT / weight
              ELSEIF (Corr(l,5).NE.0) THEN
                CALL random_number(ranu)
                PH  = EXP(yi*1000*ranu)
                Pmig(2,ii,jj,kk) = Pmig(2,ii,jj,kk) + PH / weight
              ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              IF ((tpps.LT.maxT-2).AND.(twf3.GT.(mute+Pdel)*srat) &
                                  .AND.(Corr(l,6).NE.0))  THEN     ! PpS
                AmpCor(:) = Projection(3,:) * PuP_s * PuS_2 * sprd
                WF  = DOT_PRODUCT(wfrms(l,twf3,:,3),AmpCor(:))* focus(3)
                phC = DOT_PRODUCT(ansig(l,twf3,:,3),AmpCor(:))
                PH  = EXP(PhRa*yi*ATAN2(AIMAG(phC),REAL(phC)))
                  PH  = PH * focus(3)
                NT  = ABS(WF)**(1/nval) * SIGN(1.,WF)         * focus(3)
                weight = rn(6) / phaw(3) / RECinfo(lS,4)
                Emig(3,ii,jj,kk) = Emig(3,ii,jj,kk) + WF / weight
                Pmig(3,ii,jj,kk) = Pmig(3,ii,jj,kk) + PH / weight
                Nmig(3,ii,jj,kk) = Nmig(3,ii,jj,kk) + NT / weight
              ELSEIF (Corr(l,6).NE.0) THEN
                CALL random_number(ranu)
                PH  = EXP(yi*1000*ranu)
                Pmig(3,ii,jj,kk) = Pmig(3,ii,jj,kk) + PH / weight
              ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              IF ((tpss.LT.maxT-2).AND.(twf4.GT.(mute+Pdel)*srat) &
                                  .AND.(Corr(l,7).NE.0))  THEN     ! PsSv
                AmpCor(:) = Projection(4,:) * PuS_s * SuS_1 * svc * sprd
                WF  = DOT_PRODUCT(wfrms(l,twf4,:,4),AmpCor(:))* focus(4)
                phC = DOT_PRODUCT(ansig(l,twf4,:,4),AmpCor(:))
                PH  = EXP(PhRa*yi*ATAN2(AIMAG(phC),REAL(phC)))
                  PH  = PH * focus(4)
                NT  = ABS(WF)**(1/nval) * SIGN(1.,WF)         * focus(4)
                weight = rn(7) / phaw(4) / RECinfo(lS,4)
                Emig(4,ii,jj,kk) = Emig(4,ii,jj,kk) + WF / weight
                Pmig(4,ii,jj,kk) = Pmig(4,ii,jj,kk) + PH / weight
                Nmig(4,ii,jj,kk) = Nmig(4,ii,jj,kk) + NT / weight
              ELSEIF (Corr(l,7).NE.0) THEN
                CALL random_number(ranu)
                PH  = EXP(yi*1000*ranu)
                Pmig(4,ii,jj,kk) = Pmig(4,ii,jj,kk) + PH / weight
              ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              IF ((tpss.LT.maxT-2).AND.(twf4.GT.(mute+Pdel)*srat) &
                                  .AND.(Corr(l,8).NE.0))  THEN     ! PsSh
                AmpCor(:) = Projection(5,:) * PuS_s * SuS_2 * shc * sprd
                WF  = DOT_PRODUCT(wfrms(l,twf4,:,4),AmpCor(:))* focus(5)
                phC = DOT_PRODUCT(ansig(l,twf4,:,4),AmpCor(:))
                PH  = EXP(PhRa*yi*ATAN2(AIMAG(phC),REAL(phC)))
                  PH  = PH * focus(5)
                NT  = ABS(WF)**(1/nval) * SIGN(1.,WF)         * focus(5)
                weight = rn(8) / phaw(5) / RECinfo(lS,4)
                Emig(5,ii,jj,kk) = Emig(5,ii,jj,kk) + WF / weight
                Pmig(5,ii,jj,kk) = Pmig(5,ii,jj,kk) + PH / weight
                Nmig(5,ii,jj,kk) = Nmig(5,ii,jj,kk) + NT / weight
              ELSEIF (Corr(l,8).NE.0) THEN
                CALL random_number(ranu)
                PH  = EXP(yi*1000*ranu)
                Pmig(5,ii,jj,kk) = Pmig(5,ii,jj,kk) + PH / weight
              ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ENDDO
          ENDDO
        ENDDO
      ENDDO
      WRITE(*,*) 
      Emig(6,:,:,:) = 0 
      Pmig(6,:,:,:) = 0 
      Nmig(6,:,:,:) = 0 

!     Sum the fields for the joint migration (PWS)
      Emig(6,:,:,:) = Emig(1,:,:,:) + Emig(2,:,:,:) + &
                        Emig(3,:,:,:) + Emig(4,:,:,:) + Emig(5,:,:,:)
      Pmig(6,:,:,:) = Pmig(1,:,:,:) + Pmig(2,:,:,:) + &
                        Pmig(3,:,:,:) + Pmig(4,:,:,:) + Pmig(5,:,:,:)
      Nmig(6,:,:,:) = Nmig(1,:,:,:) + Nmig(2,:,:,:) + &
                        Nmig(3,:,:,:) + Nmig(4,:,:,:) + Nmig(5,:,:,:)

!     Take the norm of the complex phase vector (PWS)
      Pmig(1,:,:,:) = ABS(Pmig(1,:,:,:))
      Pmig(2,:,:,:) = ABS(Pmig(2,:,:,:))
      Pmig(3,:,:,:) = ABS(Pmig(3,:,:,:))
      Pmig(4,:,:,:) = ABS(Pmig(4,:,:,:))
      Pmig(5,:,:,:) = ABS(Pmig(5,:,:,:))
      Pmig(6,:,:,:) = ABS(Pmig(6,:,:,:))

!     Compute the phase weighted stack F=E*P^v
      Fmig(1,:,:,:) = Emig(1,:,:,:) * (REAL(Pmig(1,:,:,:))**1)
      Fmig(2,:,:,:) = Emig(2,:,:,:) * (REAL(Pmig(2,:,:,:))**1)
      Fmig(3,:,:,:) = Emig(3,:,:,:) * (REAL(Pmig(3,:,:,:))**1)
      Fmig(4,:,:,:) = Emig(4,:,:,:) * (REAL(Pmig(4,:,:,:))**1)
      Fmig(5,:,:,:) = Emig(5,:,:,:) * (REAL(Pmig(5,:,:,:))**1)
      Fmig(6,:,:,:) = Emig(6,:,:,:) * (REAL(Pmig(6,:,:,:))**1)

!     Take the Nth power for the Nth root stacking
      DO i=grdinfo(12),maxLon             ! Longitude
        ii = i - grdinfo(12) + 1
        DO j=grdinfo(11),maxLat           ! Latitude
          jj = j - grdinfo(11) + 1
          DO k=grdinfo(10),maxDep         ! Depth
            kk = k - grdinfo(10) + 1
            IF (Nmig(1,ii,jj,kk).LT.0) THEN
              Nmig(1,ii,jj,kk) = -ABS(Nmig(1,ii,jj,kk))**(nval)
            ELSE
              Nmig(1,ii,jj,kk) =  ABS(Nmig(1,ii,jj,kk))**(nval)
            ENDIF
            IF (Nmig(2,ii,jj,kk).LT.0) THEN
              Nmig(2,ii,jj,kk) = -ABS(Nmig(2,ii,jj,kk))**(nval)
            ELSE
              Nmig(2,ii,jj,kk) =  ABS(Nmig(2,ii,jj,kk))**(nval)
            ENDIF
            IF (Nmig(3,ii,jj,kk).LT.0) THEN
              Nmig(3,ii,jj,kk) = -ABS(Nmig(3,ii,jj,kk))**(nval)
            ELSE
              Nmig(3,ii,jj,kk) =  ABS(Nmig(3,ii,jj,kk))**(nval)
            ENDIF
            IF (Nmig(4,ii,jj,kk).LT.0) THEN
              Nmig(4,ii,jj,kk) = -ABS(Nmig(4,ii,jj,kk))**(nval)
            ELSE
              Nmig(4,ii,jj,kk) =  ABS(Nmig(4,ii,jj,kk))**(nval)
            ENDIF
            IF (Nmig(5,ii,jj,kk).LT.0) THEN
              Nmig(5,ii,jj,kk) = -ABS(Nmig(5,ii,jj,kk))**(nval)
            ELSE
              Nmig(5,ii,jj,kk) =  ABS(Nmig(5,ii,jj,kk))**(nval)
            ENDIF
            IF (Nmig(6,ii,jj,kk).LT.0) THEN
              Nmig(6,ii,jj,kk) = -ABS(Nmig(6,ii,jj,kk))**(nval)
            ELSE
              Nmig(6,ii,jj,kk) =  ABS(Nmig(6,ii,jj,kk))**(nval)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!-----------------------------------------------------!
!     Renormalize the migration for visualization     ! 
!-----------------------------------------------------!
      DO i=grdinfo(12),maxLon             ! Longitude
        ii = i - grdinfo(12) + 1
        DO j=grdinfo(11),maxLat           ! Latitude
          jj = j - grdinfo(11) + 1
          DO k=grdinfo(10),maxDep         ! Depth
            kk = k - grdinfo(10) + 1
              Emax(1) =  MAX( ABS(Emig(1,ii,jj,kk)) , Emax(1)  )
              Emax(2) =  MAX( ABS(Emig(2,ii,jj,kk)) , Emax(2)  )
              Emax(3) =  MAX( ABS(Emig(3,ii,jj,kk)) , Emax(3)  )
              Emax(4) =  MAX( ABS(Emig(4,ii,jj,kk)) , Emax(4)  )
              Emax(5) =  MAX( ABS(Emig(5,ii,jj,kk)) , Emax(5)  )
              Emax(6) =  MAX( ABS(Emig(6,ii,jj,kk)) , Emax(6)  )

              Emax(11) = MAX( ABS(Pmig(1,ii,jj,kk)) , Emax(11) )
              Emax(12) = MAX( ABS(Pmig(2,ii,jj,kk)) , Emax(12) )
              Emax(13) = MAX( ABS(Pmig(3,ii,jj,kk)) , Emax(13) )
              Emax(14) = MAX( ABS(Pmig(4,ii,jj,kk)) , Emax(14) )
              Emax(15) = MAX( ABS(Pmig(5,ii,jj,kk)) , Emax(15) )
              Emax(16) = MAX( ABS(Pmig(6,ii,jj,kk)) , Emax(16) )

              Emax(21) = MAX( ABS(Fmig(1,ii,jj,kk)) , Emax(21) )
              Emax(22) = MAX( ABS(Fmig(2,ii,jj,kk)) , Emax(22) )
              Emax(23) = MAX( ABS(Fmig(3,ii,jj,kk)) , Emax(23) )
              Emax(24) = MAX( ABS(Fmig(4,ii,jj,kk)) , Emax(24) )
              Emax(25) = MAX( ABS(Fmig(5,ii,jj,kk)) , Emax(25) )
              Emax(26) = MAX( ABS(Fmig(6,ii,jj,kk)) , Emax(26) )

              Emax(31) = MAX( ABS(Nmig(1,ii,jj,kk)) , Emax(31) )
              Emax(32) = MAX( ABS(Nmig(2,ii,jj,kk)) , Emax(32) )
              Emax(33) = MAX( ABS(Nmig(3,ii,jj,kk)) , Emax(33) )
              Emax(34) = MAX( ABS(Nmig(4,ii,jj,kk)) , Emax(34) )
              Emax(35) = MAX( ABS(Nmig(5,ii,jj,kk)) , Emax(35) )
              Emax(36) = MAX( ABS(Nmig(6,ii,jj,kk)) , Emax(36) )
          ENDDO
        ENDDO
      ENDDO

      Emig(1,:,:,:) = Emig(1,:,:,:) / Emax(1)
      Emig(2,:,:,:) = Emig(2,:,:,:) / Emax(2)
      Emig(3,:,:,:) = Emig(3,:,:,:) / Emax(3)
      Emig(4,:,:,:) = Emig(4,:,:,:) / Emax(4)
      Emig(5,:,:,:) = Emig(5,:,:,:) / Emax(5)
      Emig(6,:,:,:) = Emig(6,:,:,:) / Emax(6)

      Pmig(1,:,:,:) = Pmig(1,:,:,:) / Emax(11)
      Pmig(2,:,:,:) = Pmig(2,:,:,:) / Emax(12)
      Pmig(3,:,:,:) = Pmig(3,:,:,:) / Emax(13)
      Pmig(4,:,:,:) = Pmig(4,:,:,:) / Emax(14)
      Pmig(5,:,:,:) = Pmig(5,:,:,:) / Emax(15)
      Pmig(6,:,:,:) = Pmig(6,:,:,:) / Emax(16)

      Fmig(1,:,:,:) = Fmig(1,:,:,:) / Emax(21)
      Fmig(2,:,:,:) = Fmig(2,:,:,:) / Emax(22)
      Fmig(3,:,:,:) = Fmig(3,:,:,:) / Emax(23)
      Fmig(4,:,:,:) = Fmig(4,:,:,:) / Emax(24)
      Fmig(5,:,:,:) = Fmig(5,:,:,:) / Emax(25)
      Fmig(6,:,:,:) = Fmig(6,:,:,:) / Emax(26)

      Nmig(1,:,:,:) = Nmig(1,:,:,:) / Emax(31)
      Nmig(2,:,:,:) = Nmig(2,:,:,:) / Emax(32)
      Nmig(3,:,:,:) = Nmig(3,:,:,:) / Emax(33)
      Nmig(4,:,:,:) = Nmig(4,:,:,:) / Emax(34)
      Nmig(5,:,:,:) = Nmig(5,:,:,:) / Emax(35)
      Nmig(6,:,:,:) = Nmig(6,:,:,:) / Emax(36)

!       Talk to people :)
      CALL system_clock(t2,rate)
      nbhrs = ((t2-t1)/rate/3600)
      nbmns = ((t2-t1)/rate/60) - (nbhrs*60)
      WRITE(*,*) 
      WRITE(*,*) 'Time to migrate :      ',(t2-t1)/rate,' seconds'
      WRITE(*,*) '                      (',nbhrs,' hours',nbmns,' mins)'

      WRITE(*,*) 
      WRITE(*,*) '    All traces migrated' 

                      !***************************!
                      !*                         *!
!=====================!*    Write the outputs    *!====================!
                      !*                         *!
                      !***************************!

!       Write output file name
      OPEN(250,file='../Output/Migrated/E_mig_ps',  form='formatted')
      OPEN(251,file='../Output/Migrated/E_mig_ppp', form='formatted')
      OPEN(252,file='../Output/Migrated/E_mig_pps', form='formatted')
      OPEN(253,file='../Output/Migrated/E_mig_pssv',form='formatted')
      OPEN(254,file='../Output/Migrated/E_mig_pssh',form='formatted')
      OPEN(255,file='../Output/Migrated/E_mig_all', form='formatted')
      !OPEN(255,file='../Output/Migrated/E_mig_lin',form='formatted')

      OPEN(260,file='../Output/Migrated/P_mig_ps',  form='formatted')
      OPEN(261,file='../Output/Migrated/P_mig_ppp', form='formatted')
      OPEN(262,file='../Output/Migrated/P_mig_pps', form='formatted')
      OPEN(263,file='../Output/Migrated/P_mig_pssv',form='formatted')
      OPEN(264,file='../Output/Migrated/P_mig_pssh',form='formatted')
      OPEN(265,file='../Output/Migrated/P_mig_all', form='formatted')
      !OPEN(265,file='../Output/Migrated/P_mig_lin',form='formatted')

      OPEN(270,file='../Output/Migrated/F_mig_ps',  form='formatted')
      OPEN(271,file='../Output/Migrated/F_mig_ppp', form='formatted')
      OPEN(272,file='../Output/Migrated/F_mig_pps', form='formatted')
      OPEN(273,file='../Output/Migrated/F_mig_pssv',form='formatted')
      OPEN(274,file='../Output/Migrated/F_mig_pssh',form='formatted')
      OPEN(275,file='../Output/Migrated/F_mig_all', form='formatted')
      !OPEN(275,file='../Output/Migrated/F_mig_lin',form='formatted')

      OPEN(280,file='../Output/Migrated/N_mig_ps',  form='formatted')
      OPEN(281,file='../Output/Migrated/N_mig_ppp', form='formatted')
      OPEN(282,file='../Output/Migrated/N_mig_pps', form='formatted')
      OPEN(283,file='../Output/Migrated/N_mig_pssv',form='formatted')
      OPEN(284,file='../Output/Migrated/N_mig_pssh',form='formatted')
      OPEN(285,file='../Output/Migrated/N_mig_all', form='formatted')
      !OPEN(285,file='../Output/Migrated/F_mig_lin',form='formatted')

!       Normalize max amplitude by number of migrated traces
      Emax = Emax / recount
      WRITE(250,*) Emax(1)
      WRITE(251,*) Emax(2)
      WRITE(252,*) Emax(3)
      WRITE(253,*) Emax(4)
      WRITE(254,*) Emax(5)
      WRITE(255,*) Emax(6)
      WRITE(260,*) Emax(11)
      WRITE(261,*) Emax(12)
      WRITE(262,*) Emax(13)
      WRITE(263,*) Emax(14)
      WRITE(264,*) Emax(15)
      WRITE(265,*) Emax(16)
      WRITE(270,*) Emax(21)
      WRITE(271,*) Emax(22)
      WRITE(272,*) Emax(23)
      WRITE(273,*) Emax(24)
      WRITE(274,*) Emax(25)
      WRITE(275,*) Emax(26)
      WRITE(280,*) Emax(31)
      WRITE(281,*) Emax(32)
      WRITE(282,*) Emax(33)
      WRITE(283,*) Emax(34)
      WRITE(284,*) Emax(35)
      WRITE(285,*) Emax(36)

!       Maximum amplitude on a selection on images
      WRITE(*,*) 
      WRITE(*,*) '    Emax LIN',   Emax(6)
      WRITE(*,*) '      - PS    ', Emax(1)
      WRITE(*,*) '      - PpP   ', Emax(2)
      WRITE(*,*) '      - PpS   ', Emax(3)
      WRITE(*,*) '      - PsSv  ', Emax(4)
      WRITE(*,*) '      - PsSh  ', Emax(5)
      WRITE(*,*) '    Emax PWS',   Emax(26)
      WRITE(*,*) '    Emax NTH',   Emax(36)

!       Write to the output
      DO i=grdinfo(12),maxLon               ! Longitude
        ii = i - grdinfo(12) + 1
        DO j=grdinfo(11),maxLat             ! Latitude
          jj = j - grdinfo(11) + 1
          DO k=grdinfo(10),maxDep           ! Depth
            kk = k - grdinfo(10) + 1

            WRITE(250,*) i, j, k, Emig(1,ii,jj,kk)
            WRITE(251,*) i, j, k, Emig(2,ii,jj,kk)
            WRITE(252,*) i, j, k, Emig(3,ii,jj,kk)
            WRITE(253,*) i, j, k, Emig(4,ii,jj,kk)
            WRITE(254,*) i, j, k, Emig(5,ii,jj,kk)
            WRITE(255,*) i, j, k, Emig(6,ii,jj,kk)

            !WRITE(260,*) i, j, k, REAL(Pmig(1,ii,jj,kk))
            !WRITE(261,*) i, j, k, REAL(Pmig(2,ii,jj,kk))
            !WRITE(262,*) i, j, k, REAL(Pmig(3,ii,jj,kk))
            !WRITE(263,*) i, j, k, REAL(Pmig(4,ii,jj,kk))
            !WRITE(264,*) i, j, k, REAL(Pmig(5,ii,jj,kk))
            WRITE(265,*) i, j, k, REAL(Pmig(6,ii,jj,kk))

            !WRITE(270,*) i, j, k, Fmig(1,ii,jj,kk)
            !WRITE(271,*) i, j, k, Fmig(2,ii,jj,kk)
            !WRITE(272,*) i, j, k, Fmig(3,ii,jj,kk)
            !WRITE(273,*) i, j, k, Fmig(4,ii,jj,kk)
            !WRITE(274,*) i, j, k, Fmig(5,ii,jj,kk)
            WRITE(275,*) i, j, k, Fmig(6,ii,jj,kk)

            !WRITE(280,*) i, j, k, Nmig(1,ii,jj,kk)
            !WRITE(281,*) i, j, k, Nmig(2,ii,jj,kk)
            !WRITE(282,*) i, j, k, Nmig(3,ii,jj,kk)
            !WRITE(283,*) i, j, k, Nmig(4,ii,jj,kk)
            !WRITE(284,*) i, j, k, Nmig(5,ii,jj,kk)
            WRITE(285,*) i, j, k, Nmig(6,ii,jj,kk)

          ENDDO
        ENDDO
      ENDDO
      CALL system_clock(t3,rate)
      WRITE(*,*) 
      WRITE(*,*) 'Time to write   :      ',(t3-t2)/rate,' seconds'

!       Close the relevant files
      CLOSE(250)
      CLOSE(251)
      CLOSE(252)
      CLOSE(253)
      CLOSE(254)
      CLOSE(255)
      CLOSE(260)
      CLOSE(261)
      CLOSE(262)
      CLOSE(263)
      CLOSE(264)
      CLOSE(265)
      CLOSE(270)
      CLOSE(271)
      CLOSE(272)
      CLOSE(273)
      CLOSE(274)
      CLOSE(275)
      CLOSE(280)
      CLOSE(281)
      CLOSE(282)
      CLOSE(283)
      CLOSE(284)
      CLOSE(285)

!       Deallocate memory
      DEALLOCATE(Ptt,Stt,Att,Btt,Ctt)
      DEALLOCATE(RECinfo,SRCinfo)
      DEALLOCATE(Emig,Pmig,Fmig,Nmig)
      DEALLOCATE(wfrms,Corr,tp1,bsip)
      DEALLOCATE(ansig)
      DEALLOCATE(grdinfo,CN)

!       Talk to people :)
      WRITE(*,*) 
      WRITE(*,*) '    Output written'
      WRITE(*,*) 


                           !*****************!
                           !*               *!
!==========================!*    The end    *!=========================!
                           !*               *!
                           !*****************!
      END
