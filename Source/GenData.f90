      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!                                                         !!
      !!                                                         !!
      !!        GENERATE THE DATA FOR THE MK3D MIGRATION         !!
      !!                                                         !!
      !!                                                         !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!                            !                            !!
      !!      Florian Millet        !         UCB Lyon 1         !!
      !!                            !         Ui Bergen          !!
      !!                            !                       2018 !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM GenData
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
      INTEGER :: kase                       ! Real vs synth geometry
      INTEGER :: kase2                      ! "Slice.f90" 2D projection
      REAL    :: P1,P2                      ! for lat
      REAL    :: L1,L2,L3                   ! for lon
      REAL    :: A1,A2,A3                   ! for dist
      REAL    :: B0,B1,B2,B3,B4             ! for bear/baz
      REAL    :: DIST                       ! Distance btw src & rec
      REAL    :: BEAR                       ! Bearing
      INTEGER :: oric,ranc,orid,rand        ! Sources Geometry
      INTEGER :: depp,dvar,mode             ! Sources Geometry
      REAL    :: alat,alon                  ! Array center
      REAL    :: asiz,angs,ratio            ! Array dimension & geom
      INTEGER :: dim1,dim2                  ! Size of rect. array
      REAL    :: nrow                       ! Nb of rows in array
      REAL    :: sspa                       ! Angular source spacing
      REAL    :: ranu                       ! Random number
      REAL    :: TOangle      !not used!    ! Take-off angle (ttimes)
      INTEGER :: prec                       ! Precision of projection
      REAL    :: latp0,lonp0,latp1,lonp1    ! Projection coordinates
      REAL    :: d1,d2,d3                   ! Rec - proj distance
      REAL,DIMENSION(2) :: curpt            ! Current proj point
      REAL,DIMENSION(2) :: LLC              ! Lower Left Corner 
      INTEGER :: splmnt                     ! Value for numerotation
      REAL    :: sdr                        ! Station Density Radius
      INTEGER :: nwfs                       ! Number of waveforms
      INTEGER :: nspl                       ! Number of sampling point
      REAL    :: dt, shft                   ! Sampl rate & align t shft 
!       Arrays
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: corr      ! Src-Rec table
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: sel       ! Selection table
      REAL,DIMENSION(:,:),   ALLOCATABLE :: RECinfo   ! Receivers info
      REAL,DIMENSION(:,:),   ALLOCATABLE :: REC_2     ! Proj Rec info
      INTEGER,DIMENSION(:),  ALLOCATABLE :: nsir      ! Nb Sta In Rad
      REAL,DIMENSION(:,:),   ALLOCATABLE :: sd        ! Station Density
      REAL,DIMENSION(:,:),   ALLOCATABLE :: SRCinfo   ! Sources info
      REAL,DIMENSION(:),     ALLOCATABLE :: BAZ       ! Back-azimuth
      REAL,DIMENSION(:),     ALLOCATABLE :: SLWN      ! Slowness
      REAL,DIMENSION(:),     ALLOCATABLE :: NSshift   ! N-S shift (RS)
      REAL,DIMENSION(:),     ALLOCATABLE :: EWshift   ! E-W shift (RS)
      REAL,DIMENSION(:),     ALLOCATABLE :: diff      ! r-p distance
!           ( analytical signal )
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: wf         ! Waveforms
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: Hilb       ! Hilb. Trans.
      COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: AS       ! Anal. Sig.
      COMPLEX,DIMENSION(:) ,ALLOCATABLE :: ipR,ipT,ipZ ! Fourier PHA
      COMPLEX,DIMENSION(:) ,ALLOCATABLE :: sgn         ! Sign function
      REAL,DIMENSION(:)    ,ALLOCATABLE :: r,t,z       ! RTZ in time
      COMPLEX,DIMENSION(:) ,ALLOCATABLE :: Fr,Ft,Fz    ! RTZ in Four
      REAL,DIMENSION(:)    ,ALLOCATABLE :: fft1,fft2,fft3
!       File names
      CHARACTER(LEN=30) :: infile           ! Input file name
      CHARACTER(LEN=120):: TMPcmd           ! Temp command (TauP)
      CHARACTER(LEN=30) :: RSname           ! Raysum input file
      CHARACTER(LEN=30) :: RSfile           ! Raysum input file
      CHARACTER(LEN=30) :: COname           ! SRC-REC corresp file
      CHARACTER(LEN=30) :: COfile           ! SRC-REC corresp file
      CHARACTER(LEN=30) :: SRname           ! SrcRecLoc geometry file
      CHARACTER(LEN=30) :: SRfile           ! SrcRecLoc geometry file
      CHARACTER(LEN=30) :: WFname           ! Waveforms file
      CHARACTER(LEN=30) :: WFfile           ! Waveforms file
      CHARACTER(LEN=30) :: FMfile           ! FM3D input file
      CHARACTER(LEN=30) :: srcnam,recnam,cornam,wfsnam,selnam ! Field
      CHARACTER(LEN=40) :: Rsrc  ,Rrec  ,Rcor  ,Rwfs,  Rsel   ! Field
!       Counters, flags & dummies
      INTEGER :: i, j, k, l, ii             ! Counters
      INTEGER :: records                    ! Number of records present
      CHARACTER(LEN=30) :: dummy            ! Char dummy
      INTEGER :: dum                        ! Integer dummy
      REAL    :: dumm                       ! Real dummy


                    !******************************!
                    !*                            *!
!===================!*    Read the input files    *!===================!
                    !*                            *!
                    !******************************!

!------------------------------------------!
!     Prepare all the files and arrays     !
!------------------------------------------!
!       Write input files names
      infile = '../Input/GenData.in'       ! Name of input file
!       Open the files
      OPEN(114,file=infile,form="formatted")

!--------------------------------!
!     Prepare the parameters     !
!--------------------------------!
!       Read the information
      READ(114,*) !#
      READ(114,*) !# Case
      READ(114,*) !#
      READ(114,*) kase
      READ(114,*) !#
      READ(114,*) !# Numbers
      READ(114,*) !#
      READ(114,*) nsrc
      READ(114,*) oric
      READ(114,*) ranc
      READ(114,*) orid
      READ(114,*) rand
      READ(114,*) depp
      READ(114,*) dvar
      READ(114,*) nrec
      READ(114,*) mode, angs, dim1
      READ(114,*) asiz
      READ(114,*) alat, alon
      READ(114,*) !#
      READ(114,*) !# Files
      READ(114,*) !#
      READ(114,*) RSname
      READ(114,*) COname
      READ(114,*) SRname
      READ(114,*) WFname
      READ(114,*) !#
      READ(114,*) !# Input files
      READ(114,*) !#
      READ(114,*) srcnam
      READ(114,*) recnam
      READ(114,*) cornam
      READ(114,*) wfsnam
      READ(114,*) selnam
      READ(114,*) !#
      READ(114,*) !# Projection
      READ(114,*) !#
      READ(114,*) kase2
      READ(114,*) prec
      READ(114,*) splmnt
      CLOSE(114)

!        If field data, find nsrc & nrec
      IF ((kase .EQ. 2).OR.(kase .EQ. 3)) THEN
        Rsrc = '../Data/' // srcnam        ! Sources
        Rrec = '../Data/' // recnam        ! Receivers
        OPEN(301,file=Rsrc, form="formatted")
        OPEN(302,file=Rrec, form="formatted")
        READ(301,*) nsrc
        READ(302,*) nrec
        CLOSE(301)
        CLOSE(302)
      ENDIF

!        Talk to people :)
      WRITE(*,*)  ''
      WRITE(*,*)  '   +-------------------------------+'
      WRITE(*,*)  '   |                               |'
      WRITE(*,11) '|   Number of sources   : ',nsrc,'   |'
      WRITE(*,11) '|   Number of receivers : ',nrec,'   |'
      WRITE(*,*)  '   |                               |'
      WRITE(*,*)  '   +-------------------------------+'
      WRITE(*,*)  ''
11    FORMAT(A30,I3,A4)

!        Talk to people :)
      IF (kase.EQ.1) WRITE(*,*)  ' Case ~ Type: Synthetics'
      IF (kase.EQ.2) WRITE(*,*)  ' Case ~ Type: Real data'
      IF (kase.EQ.3) WRITE(*,*)  ' Case ~ Type: Resolution test'
      IF (kase2.EQ.1) WRITE(*,*) '  Projection: No'
      IF (kase2.EQ.2) WRITE(*,*) '  Projection: Yes'
      WRITE(*,*)


                  !***********************************!
                  !*                                 *!
!=================!*    Create synthetic geometry    *!================!
                  !*                                 *!
                  !***********************************!

!       Check case scenario
      IF (kase .EQ. 1) THEN

!       Allocate memory
        ALLOCATE(SRCinfo(nsrc,3))
        ALLOCATE(RECinfo(nrec,4))
        ALLOCATE(corr(   nrec*nsrc,8))
        ALLOCATE(BAZ(    nrec*nsrc),SLWN(   nrec*nsrc))
        ALLOCATE(NSshift(nrec*nsrc),EWshift(nrec*nsrc))

!-----------------------------------!
!     Create the array geometry     !
!-----------------------------------!
!         Check if number of receiver is correct
        IF (mode .EQ. 1) THEN
          nrow = SQRT(FLOAT(nrec))     ! Square
          IF (MOD(nrow,1.).NE.0) THEN
            WRITE(*,*) 'INPUT ERROR'
            WRITE(*,*) 
            WRITE(*,*) '  >>  Incorrect number of receivers  <<  '
            WRITE(*,*) 
            STOP
          ENDIF
        ELSE IF (mode .EQ. 3) THEN     ! Rectangle
          IF (MOD((FLOAT(nrec)/FLOAT(dim1)),1.).NE.0) THEN
            WRITE(*,*) 'INPUT ERROR'
            WRITE(*,*) 
            WRITE(*,*) '  >>  Incorrect number of receivers  <<  '
            WRITE(*,*) 
            STOP
          ENDIF
        ENDIF

!         -- square --
        IF (mode .EQ. 1) THEN
          LLC(1) = alat - asiz / 2               ! Lower Left Corner
          LLC(2) = alon - asiz / 2               ! Lower Left Corner
          IF(nrow.NE.1) THEN
            DO i=1,INT(nrow)
              DO j=1,INT(nrow)
                k = (i-1) * nrow + j
                RECinfo(k,1) = LLC(1) + (i-1) * asiz / (nrow-1)
                RECinfo(k,2) = LLC(2) + (j-1) * asiz / (nrow-1)
                RECinfo(k,3) = ZERO              ! Depth
                RECinfo(k,4) = ONE               ! Weight
              ENDDO
            ENDDO
          ELSE
            RECinfo(1,1) = alat                  ! Latitude
            RECinfo(1,2) = alon                  ! Longitude
            RECinfo(1,3) = ZERO                  ! Depth
            RECinfo(1,4) = ONE                   ! Weight
          ENDIF

!         -- linear --
        ELSE IF (mode .EQ. 2) THEN
          LLC(1) = alat - (asiz / 2 * COS(angs*3.14/180))
          LLC(2) = alon - (asiz / 2 * SIN(angs*3.14/180))
          DO i=1,nrec
            RECinfo(i,1) = LLC(1)+(i-1)*asiz/(nrec-1)*COS(angs*3.14/180)
            RECinfo(i,2) = LLC(2)+(i-1)*asiz/(nrec-1)*SIN(angs*3.14/180)
            RECinfo(i,3) = ZERO                  ! Depth
            RECinfo(i,4) = ONE                   ! Weight
          ENDDO

!         -- rectangle --
        ELSE IF (mode .EQ. 3) THEN
          ratio  = nrec/FLOAT(dim1)/FLOAT(dim1)  ! Length vs. width
          LLC(1) = alat &
                    - asiz/2 * COS(angs*3.14/180) * ratio &
                    - asiz/2 * SIN(angs*3.14/180)
          LLC(2) = alon &
                    - asiz/2 * SIN(angs*3.14/180) * ratio &
                    + asiz/2 * COS(angs*3.14/180)
          dim2 = INT(nrec)/dim1
          DO i=1,dim1
            DO j=1,dim2!INT(nrec)/dim1
              k = (i-1) * nrec/dim1 + j
              RECinfo(k,1) = LLC(1) &
                  + (j-1)*asiz/(dim2-1)*COS(angs*3.14/180) * ratio &
                  + (i-1)*asiz/(dim1-1)*SIN(angs*3.14/180)
              RECinfo(k,2) = LLC(2) &
                  + (j-1)*asiz/(dim2-1)*SIN(angs*3.14/180) * ratio &
                  - (i-1)*asiz/(dim1-1)*COS(angs*3.14/180)
              RECinfo(k,3)= ZERO                 ! Depth
              RECinfo(k,4)= ONE                  ! Weight
            ENDDO
          ENDDO
        ENDIF

!          Talk to people :)
        WRITE(*,*) ' Array geometry created'
        WRITE(*,*)

!------------------------------------!
!     Create the source geometry     !
!------------------------------------!
        sspa = ranc / nsrc
        DO i=1,nsrc

!           Bearing and distance center-source
          CALL random_number(ranu)
          BEAR = ((i-1) * sspa + oric) * dtr
          DIST = (orid + ranu * rand) * dtr

!           Latitude
          P1 = SIN(alat * dtr) * COS(DIST)
          P2 = COS(alat * dtr) * SIN(DIST) * COS(BEAR)
          SRCinfo(i,1) = ASIN(P1+P2) / dtr

!           Longitude
          L1 = SIN(BEAR) * SIN(DIST) * COS(alat * dtr)
          L2 = COS(DIST)
          L3 = SIN(alat * dtr) * SIN(SRCinfo(i,1) * dtr)
          SRCinfo(i,2) = alon + ATAN2(L1,L2-L3) / dtr

!           Depth
          CALL random_number(ranu)
          !SRCinfo(i,3) = depp + ((2*ranu)-1) * dvar
          SRCinfo(i,3) = depp

        ENDDO

!          Talk to people :)
        WRITE(*,*) ' Source geometry created'
        WRITE(*,*)

!---------------------------------!
!     Compute the 3D geometry     !
!---------------------------------!
!         Start writing file for ttimes
        OPEN(201,file='input_for_ttimes',form='formatted')
        WRITE(201,'(A1)') 'Y'
        WRITE(201,'(A8)') 'ak135_tt'
        WRITE(201,'(A1)') 'P'
        WRITE(201,'(A1)') ' '
        WRITE(201,'(A1)') 'N'

        DO i=1,nsrc
          DO j=1,nrec
            k = (i-1) * nrec + j

!         Distance
            A1 = SIN(SRCinfo(i,1) * dtr) * SIN(RECinfo(j,1) * dtr)
            A2 = COS(SRCinfo(i,1) * dtr) * COS(RECinfo(j,1) * dtr)
            A3 = COS(SRCinfo(i,2) * dtr - RECinfo(j,2) * dtr)
            DIST = ACOS(A1+A2*A3)/dtr

!         Back-azimuth
            B0 = COS(SRCinfo(i,2) * dtr - RECinfo(j,2) * dtr)
            B1 = SIN(SRCinfo(i,2) * dtr - RECinfo(j,2) * dtr)
            B2 = COS(SRCinfo(i,1) * dtr)
            B3 = COS(RECinfo(j,1) * dtr) * SIN(SRCinfo(i,1) * dtr)
            B4 = SIN(RECinfo(j,1) * dtr) * COS(SRCinfo(i,1) * dtr)
            BAZ(k) = ATAN2(B1*B2,B3-B4*B0) / dtr
            IF (BAZ(k).LT.0) BAZ(k) = BAZ(k) + 360

!         Incidence angle with ttimes (write here, read later)
            WRITE(201,'(F6.3)') SRCinfo(i,3)
            WRITE(201,'(F6.3)') DIST
            WRITE(201,'(A2)') '-1'

!         NS and EW shifts (m)
            NSshift(k) = (RECinfo(j,1) - alat) * dkm * 1000
            EWshift(k) = (RECinfo(j,2) - alon) * dkm * 1000

!         Correspondance table
            corr(k,1) = i                ! Source
            corr(k,2) = j                ! Receiver
            corr(k,3) = 1                ! Data altogether
            corr(k,4) = 1                ! PS
            corr(k,5) = 1                ! PpP
            corr(k,6) = 1                ! PpS
            corr(k,7) = 1                ! PsSv
            corr(k,8) = 1                ! PsSh

            IF (k.EQ.1) WRITE(*,*) ' Number of records :',nrec*nsrc
            IF (k.EQ.1) WRITE(*,*) 

          ENDDO
        ENDDO
        CLOSE(201)

!         Slownesses with ttimes
        WRITE(*,*) '  > ignore next 4 lines <'
        WRITE(*,*) '---------------------------'
        CALL SYSTEM('./ttimes < input_for_ttimes > output_from_ttimes')
        WRITE(*,*) '---------------------------'
        WRITE(*,*) ! WARNING: Needs to be in the correct range to work 
        OPEN(902,file='ttimes.lst',form='formatted')
        DO i=1,nsrc
          DO j=1,nrec
            k = (i-1) * nrec + j
            READ(902,*) !
            READ(902,*) !
            READ(902,*) ! Source Depth
            READ(902,*) !
            READ(902,*) ! delta [...]
            READ(902,*) dumm, dum, dummy, dumm, TOangle, SLWN(k)
            READ(902,*) ! PKiKP
            SLWN(k) = SLWN(k) / (dkm * 1000)     ! Slowness
          ENDDO
        ENDDO
        CLOSE(902)

!         Talk to people :)
        WRITE(*,*) ' Source-Receiver couples geometries computed'
        WRITE(*,*) 


                  !**********************************!
                  !*                                *!
!=================!*    Read real world geometry    *!=================!
                  !*                                *!
                  !**********************************!

!       Check case scenario
      ELSEIF ((kase .EQ. 2).OR.(kase .EQ. 3)) THEN

!---------------------------!
!     Read the geometry     !
!---------------------------!
!         Open the files
        Rsrc = '../Data/' // srcnam        ! Sources
        Rrec = '../Data/' // recnam        ! Receivers
        Rcor = '../Data/' // cornam        ! Correspondances
        Rwfs = '../Data/' // wfsnam        ! Waveforms
        Rsel = '../Data/' // selnam        ! Quality Control
        OPEN(301,file=Rsrc, form="formatted")
        OPEN(302,file=Rrec, form="formatted")
        OPEN(303,file=Rcor, form="formatted")
        OPEN(305,file=Rsel, form="formatted")

!         Read number of sources and receivers
        READ(301,*) nsrc
        READ(302,*) nrec

!         Allocate memory
        ALLOCATE(SRCinfo(nsrc,3))
        ALLOCATE(RECinfo(nrec,4))
        ALLOCATE(REC_2(  nrec,4))
        ALLOCATE(nsir(nrec),sd(nrec,nrec))
        ALLOCATE(corr(nrec*nsrc,8))
        ALLOCATE(sel(nsrc,5))
        ALLOCATE(BAZ(    nrec*nsrc),SLWN(   nrec*nsrc))
        ALLOCATE(NSshift(nrec*nsrc),EWshift(nrec*nsrc))
        ALLOCATE(diff(nrec))

!         Read geometry information
        corr = 0 
        DO i=1,nsrc
          READ(301,*) SRCinfo(i,3), SRCinfo(i,1), SRCinfo(i,2)
        ENDDO
        DO i=1,nrec
          READ(302,*) RECinfo(i,3), RECinfo(i,1), RECinfo(i,2)
        ENDDO
        DO i=1,nsrc
          READ(305,*) dum,sel(i,1),sel(i,2),sel(i,3),sel(i,4),sel(i,5)
        ENDDO
        DO i=1,nsrc
          DO j=1,nrec
            k = (i-1)*nrec + j
            READ(303,*) corr(k,1), corr(k,2), corr(k,3)
            IF (corr(k,3).NE.0) THEN
              corr(k,4) = 1!sel(i,1)
              corr(k,5) = 1!sel(i,2)
              corr(k,6) = 1!sel(i,3)
              corr(k,7) = 1!sel(i,4)
              corr(k,8) = 1!sel(i,5)
            ENDIF
          ENDDO
        ENDDO

!------------------------------------------------!
!     Project if needed (params in Slice.in)     !
!------------------------------------------------!
!         Open the files
        OPEN(401,file="../Input/Slice.in",form="formatted")
        DO i=1,8
          READ(401,*)
        ENDDO
        READ(401,*) latp0
        READ(401,*) lonp0
        READ(401,*) !
        READ(401,*) !
        READ(401,*) latp1
        READ(401,*) lonp1
        CLOSE(401)

!         Project
        diff = 1000
        REC_2 = RECinfo
        DO i=1,prec
          curpt(1) = (latp1-latp0)*(i-1)/prec + latp0
          curpt(2) = (lonp1-lonp0)*(i-1)/prec + lonp0
          DO j=1,nrec
            d1 = (curpt(1)-RECinfo(j,1))**2
            d2 = (curpt(2)-RECinfo(j,2))**2
            d3 = SQRT(d1+d2)
            IF (d3 .LT. diff(j)) THEN
              diff(j)    = d3               ! Distance proj-original
              REC_2(j,1) = curpt(1)         ! New latitude
              REC_2(j,2) = curpt(2)         ! New longitude
            ENDIF
          ENDDO
        ENDDO

!         Write output file
        OPEN(402,file="../Data/proj_rec",form="formatted")
        WRITE(402,*) nrec
        DO i=1,nrec
          WRITE(402,*) RECinfo(i,3), REC_2(i,1), REC_2(i,2), diff(i)
        ENDDO
        CLOSE(402)

!         Attribute REC_2 to RECinfo if needed
        IF (kase2 .EQ. 2) THEN
          RECinfo = REC_2
        ENDIF

!---------------------------------!
!     Compute the 3D geometry     !
!---------------------------------!
!         Start writing ttimes file
        OPEN(201,file='input_for_ttimes',form='formatted')
        WRITE(201,'(A1)') 'Y'
        WRITE(201,'(A8)') 'ak135_tt'
        WRITE(201,'(A1)') 'P'
        WRITE(201,'(A1)') ' '
        WRITE(201,'(A1)') 'N'
        DO i=1,nsrc
          DO j=1,nrec
            k = (i-1) * nrec + j

!         Distance
            A1 = SIN(SRCinfo(i,1) * dtr) * SIN(RECinfo(j,1) * dtr)
            A2 = COS(SRCinfo(i,1) * dtr) * COS(RECinfo(j,1) * dtr)
            A3 = COS(SRCinfo(i,2) * dtr - RECinfo(j,2) * dtr)
            DIST = ACOS(A1+A2*A3)/dtr

!         Back-azimuth
            B0 = COS(SRCinfo(i,2) * dtr - RECinfo(j,2) * dtr)
            B1 = SIN(SRCinfo(i,2) * dtr - RECinfo(j,2) * dtr)
            B2 = COS(SRCinfo(i,1) * dtr)
            B3 = COS(RECinfo(j,1) * dtr) * SIN(SRCinfo(i,1) * dtr)
            B4 = SIN(RECinfo(j,1) * dtr) * COS(SRCinfo(i,1) * dtr)
            BAZ(k) = ATAN2(B1*B2,B3-B4*B0) / dtr
            IF (BAZ(k).LT.0) BAZ(k) = BAZ(k) + 360

!         NS and EW shifts (m)
            NSshift(k) = (RECinfo(j,1) - alat) * dkm * 1000
            EWshift(k) = (RECinfo(j,2) - alon) * dkm * 1000

!         correct for >100km depth
            IF ((kase.EQ.3).AND.(SRCinfo(i,3).GT.99)) THEN
              corr(k,3) = 0
            ENDIF

!         ttimes input
            WRITE(201,'(F6.3)') REAL(MIN(99,INT(SRCinfo(i,3))))
            WRITE(201,'(F6.3)') DIST
            WRITE(201,'(A2)') '-1'
          ENDDO
        ENDDO
        CLOSE(201)

!         Run ttimes for slownesses
        WRITE(*,*) '  > ignore next 4 lines <'
        WRITE(*,*) '---------------------------'
        CALL SYSTEM('./ttimes < input_for_ttimes > out_from_ttimes')
        WRITE(*,*) '---------------------------'
        WRITE(*,*) 
        OPEN(902,file='ttimes.lst',form='formatted')
        DO i=1,nsrc
          DO j=1,nrec
            k = (i-1) * nrec + j
            READ(902,*) !
            READ(902,*) !
            READ(902,*) ! Source Depth
            READ(902,*) !
            READ(902,*) ! delta [...]
            READ(902,*) dumm, dum, dummy, dumm, TOangle, SLWN(k)
            READ(902,*) ! PKiKP
            SLWN(k) = SLWN(k) / (dkm * 1000)
          ENDDO
        ENDDO
        CLOSE(902)

!--------------------------------------------------!
!     Compute station weight based on density      !
!--------------------------------------------------!
!         Compute distance from every receiver to all others
        DO i=1,nrec
          DO j=1,nrec
            sd(i,j) = SQRT(((RECinfo(i,1)-RECinfo(j,1))*dkm)**2 + &
                            ((RECinfo(i,2)-RECinfo(j,2))*dkm)**2 )
          ENDDO
        ENDDO

!         Compute average inter-station distance
        sdr = 0
        DO i=1,nrec
          dumm = 999999
          DO j=1,nrec
            IF (sd(i,j).EQ.0.0) CYCLE       ! Stations at same location
            dumm = MIN(dumm,sd(i,j))
          ENDDO
          sdr = sdr + dumm / nrec
        ENDDO

!         Find number of closest stations to every station
        nsir = ZERO
        DO i=1,nrec
          DO j=1,nrec
            IF (sd(i,j).LT.sdr) THEN
              nsir(i) = nsir(i) + 1
            ENDIF
          ENDDO
        ENDDO

!         Compute weight
        DO i=1,nrec
          RECinfo(i,4) = ( 1 / FLOAT(nsir(i)) ) ** 1
        ENDDO

!         Close files and talk to people
        CLOSE(301)
        CLOSE(302)
        CLOSE(303)
        CLOSE(305)
        WRITE(*,*) ' Source-Receiver geometry read'
        WRITE(*,*)

      ENDIF


                    !******************************!
                    !*                            *!
!===================!*    Write and Run Raysum    *!===================!
                    !*                            *!
                    !******************************!

!       Prepare the file
      RSfile  = '../Run/raysum/' // RSname  ! Raysum input file
      OPEN(203,file=RSfile, form="formatted")

!       Write to the file
      WRITE(203,*) '# BAZ (deg), SLWN (s/m), NS shift (m), EW shift (m)'
      DO i=1,nsrc*nrec
        WRITE(203,*) BAZ(i), SLWN(i), NSshift(i), EWshift(i)
      ENDDO

!       Close the relevant files
      CLOSE(203)

!       Talk to people :)
      WRITE(*,*) ' Inputs for Raysum prepared'
      WRITE(*,*) 

!       Do it
      IF ((kase .EQ. 1).OR.(kase .EQ. 3)) THEN
        WRITE(TMPcmd,*) './raysum/seis-spread raysum/mod raysum/geo', &
                         ' raysum/ph raysum/arr raysum/tr > raysum/out'
        CALL SYSTEM(TMPcmd)

!       Talk to people :)
        WRITE(*,*) ' Raysum run'
        WRITE(*,*) 
      ENDIF


               !****************************************!
               !*                                      *!
!==============!*    Compute the Analytical Signals    *!==============!
               !*                                      *!
               !****************************************!

!-----------------------------!
!     Prepare everything      !
!-----------------------------!
!       Read the input parameters
      IF ((kase .EQ. 1).OR.(kase .EQ. 3)) THEN
        OPEN(101,file='raysum/tr',form="formatted")
      ELSEIF (kase .EQ. 2) THEN
        OPEN(101,file=Rwfs,form="formatted")
      ENDIF
      READ(101,*) !#
      READ(101,*) nwfs, nspl, dt, dumm, shft

!       Allocate memory
      ALLOCATE(wf(nwfs,nspl,3))
      ALLOCATE(AS(nwfs,nspl,3))
      ALLOCATE(Hilb(nwfs,nspl,3))
      ALLOCATE(ipR(nspl),ipT(nspl),ipZ(nspl))
      ALLOCATE(sgn(nspl))
      ALLOCATE( r(nspl), t(nspl), z(nspl))
      ALLOCATE(Fr(nspl),Ft(nspl),Fz(nspl))
      ALLOCATE(fft1(2*nspl),fft2(2*nspl),fft3(2*nspl))
      wf      = 0
      records = 0

!       Read the waveforms
      DO i=1,nwfs
        READ(101,*) !--------!
        READ(101,*) ! Header !
        READ(101,*) !--------!
        IF (corr(i,3).NE.0) THEN
          records = records + 1
          DO j=1,nspl
            READ(101,*) wf(i,j,1), wf(i,j,2), wf(i,j,3)
          ENDDO
        ELSE
          IF (kase .EQ. 3) THEN        ! Remove unavailable data
            DO j=1,nspl                ! for resolution test
              READ(101,*)
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      CLOSE(101)

!-----------------------------!
!     Do the calculations     !
!-----------------------------!
!       Start the loop
      DO i=1,nwfs
        r = 0
        t = 0
        z = 0
        DO j=1,nspl
          CALL random_number(ranu)
          r(j) = wf(i,j,1)
          CALL random_number(ranu)
          t(j) = wf(i,j,2)
          CALL random_number(ranu)
          z(j) = wf(i,j,3)
        ENDDO

!         Go to Fourier domain
        Fr = 0
        Ft = 0
        Fz = 0
        fft1 = 0
        fft2 = 0
        fft3 = 0
        DO j=1,nspl
          fft1(2*j-1) = r(j)
          fft1(2*j)   = 0
          fft2(2*j-1) = t(j)
          fft2(2*j)   = 0
          fft3(2*j-1) = z(j)
          fft3(2*j)   = 0
        ENDDO
        CALL four1(fft1,nspl,1)
        CALL four1(fft2,nspl,1)
        CALL four1(fft3,nspl,1)
        DO j=1,nspl
          Fr(j) = CMPLX(fft1(2*j-1),fft1(2*j))
          Ft(j) = CMPLX(fft2(2*j-1),fft2(2*j))
          Fz(j) = CMPLX(fft3(2*j-1),fft3(2*j))
        ENDDO

!         Compute Hilbert Transform of RF (PWS)
        sgn = 0 
        ipR = 0
        ipT = 0
        ipZ = 0
        DO j=1,nspl
          IF (j .LT. nspl/2) sgn(j) = +1
          IF (j .GT. nspl/2) sgn(j) = -1
        ENDDO
        DO j=1,nspl                         ! Apply to original wfs
          ipR(j) = yi * Fr(j) * sgn(j)      ! > r [wf(1)]
          ipT(j) = yi * Ft(j) * sgn(j)      ! > t [wf(2)]
          ipZ(j) = yi * Fz(j) * sgn(j)      ! > z [wf(3)]
        ENDDO
        fft1 = 0
        fft2 = 0
        fft3 = 0
        DO j=1,nspl
          fft1(2*j-1) = REAL( ipR(j))
          fft1(2*j)   = AIMAG(ipR(j))
          fft2(2*j-1) = REAL( ipT(j))
          fft2(2*j)   = AIMAG(ipT(j))
          fft3(2*j-1) = REAL( ipZ(j))
          fft3(2*j)   = AIMAG(ipZ(j))
        ENDDO
        CALL four1(fft1,nspl,-1)
        CALL four1(fft2,nspl,-1)
        CALL four1(fft3,nspl,-1)
        DO j=1,nspl
          Hilb(i,j,1) = fft1(2*j-1) / nspl
          Hilb(i,j,2) = fft2(2*j-1) / nspl
          Hilb(i,j,3) = fft3(2*j-1) / nspl
        ENDDO

!         Compute analytical signal
        DO j=1,nspl
          AS(i,j,1) = wf(i,j,1) + yi*Hilb(i,j,1)
          AS(i,j,2) = wf(i,j,2) + yi*Hilb(i,j,2)
          AS(i,j,3) = wf(i,j,3) + yi*Hilb(i,j,3)
        ENDDO

      ENDDO

!       Talk to people :)
      WRITE(*,*) ' Analytical signals computed'
      WRITE(*,*)


                      !***************************!
                      !*                         *!
!=====================!*    Write the outputs    *!====================!
                      !*                         *!
                      !***************************!

!-------------------------------!
!     Prepare all the files     !
!-------------------------------!
      WFfile  = '../Data/' // WFname        ! RF output file
      COfile  = '../Data/' // COname        ! Src-Rec corresp file
      SRfile  = '../Data/' // SRname        ! SrcRecLoc input file
      OPEN(202,file=WFfile,form="formatted")
      OPEN(204,file=COfile, form="formatted")
      OPEN(205,file=SRfile, form="formatted")

!------------------------------!
!     Write the FM3D files     !
!------------------------------!
!       Sources
      WRITE(205,*) nsrc
      DO i=1,nsrc
        ii = i + splmnt
        WRITE(205,*)    SRCinfo(i,3), SRCinfo(i,1), SRCinfo(i,2)

!         Inputs for the sources in FM3D
        WRITE(FMfile,13) '../Output/FMins/s', i, '.in'
        OPEN( 201,file=FMfile, form="formatted")
        WRITE(201,*) '1'                    ! Teleseismic
        WRITE(201,*) 'P'                    ! P wave coda
        WRITE(201,*) SRCinfo(i,3), SRCinfo(i,1), SRCinfo(i,2)
        WRITE(201,*) '3'                    ! 3 time fields P Pp Ps
        WRITE(201,*) '2'                    !  2 legs
        WRITE(201,*) '3 2  2 1'             !  Bottom > Top
        WRITE(201,*) '1    1'               !  P wave
        WRITE(201,*) '3'                    !   3 legs
        WRITE(201,*) '3 2  2 1  1 2'        !   Bottom > Top > Back
        WRITE(201,*) '1    1    1'          !   Pp
        WRITE(201,*) '3'                    !    3 legs 
        WRITE(201,*) '3 2  2 1  1 2'        !    Bottom > Top > Back
        WRITE(201,*) '1    1    2'          !    Ps
        WRITE(201,14) 's', ii               ! Output file name
        CLOSE(201)
      ENDDO

!       Receivers
      WRITE(205,*) nrec
      DO i=1,nrec
        ii = i + splmnt
        WRITE(205,*) RECinfo(i,3),RECinfo(i,1),RECinfo(i,2),RECinfo(i,4)

!         Inputs for the receivers in FM3D
        WRITE(FMfile,13) '../Output/FMins/r', i, '.in'
        OPEN( 201,file=FMfile, form="formatted")
        WRITE(201,*) '0'                    ! Local source
        WRITE(201,*) RECinfo(i,3), RECinfo(i,1), RECinfo(i,2)
        WRITE(201,*) '2'                    ! 2 times fields S P
        WRITE(201,*) '2'                    !  2 legs
        WRITE(201,*) '0 2  2 1'             !  Top > Bottom
        WRITE(201,*) '2    2'               !  S wave
        WRITE(201,*) '2'                    !   2 legs
        WRITE(201,*) '0 2  2 1'             !   Top > Bottom
        WRITE(201,*) '1    1'               !   P wave
        WRITE(201,14) 'r', ii               ! Output file name
        CLOSE(201)
      ENDDO
13    FORMAT(A17,I4.4,A3)
14    FORMAT(A1,I4.4)

!---------------------------------------!
!     Write the correspondance file     !
!---------------------------------------!
      DO i=1,nsrc*nrec
        WRITE(204,*) corr(i,1), corr(i,2), corr(i,3), corr(i,4), &
                     corr(i,5), corr(i,6), corr(i,7), corr(i,8) 
      ENDDO

!----------------------------------!
!     Write the waveforms file     !
!----------------------------------!
      WRITE(202,*) '#traces #samples dt(s) align shift(s)'
      WRITE(202,*) nwfs, nspl, dt, ONE, shft
      DO i=1,nwfs
        WRITE(202,*) '#--------'
        WRITE(202,*) '# trace',i
        WRITE(202,*) '#--------'
        IF (corr(i,3).NE.0) THEN
          DO j=1,nspl
            WRITE(202,*) wf(i,j,1), wf(i,j,2), wf(i,j,3), &
                         AS(i,j,1), AS(i,j,2), AS(i,j,3)!, &
          ENDDO
        ENDIF
      ENDDO

!       Close the relevant files
      CLOSE(201)
      CLOSE(202)
      CLOSE(204)
      CLOSE(205)

!       Deallocate memory
      DEALLOCATE(wf)
      DEALLOCATE( r, t, z)
      DEALLOCATE(Fr,Ft,Fz)
      DEALLOCATE(fft1,fft2,fft3)
      DEALLOCATE(AS)

!       Talk to people :)
      WRITE(*,*) ' Records present :',records,'/',nwfs
      WRITE(*,*) 
      WRITE(*,*) ' Outputs written'
      WRITE(*,*) 


                           !*****************!
                           !*               *!
!==========================!*    The end    *!=========================!
                           !*               *!
                           !*****************!
      END



!========================!
!                        !
!     Function four1     !
!                        !
!========================!
      SUBROUTINE four1(data,nn,isign)
!---------------------------------------!
!     Define variables for the code     !
!---------------------------------------!
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
!------------------!
!     The code     !
!------------------!
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.
