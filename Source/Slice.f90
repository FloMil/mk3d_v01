      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!                                                         !!
      !!                                                         !!
      !!      CODE FOR SLICING 3D DATA FOR 2D VISUALISATION      !!
      !!                                                         !!
      !!                                                         !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!                            !                            !!
      !!      Florian Millet        !         UCB Lyon 1         !!
      !!                            !         Ui Bergen          !!
      !!                            !                       2018 !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM Slice
      IMPLICIT NONE


                   !*******************************!
                   !*                             *!
!==================!*    Define the variables     *!===================!
                   !*                             *!
                   !*******************************!

!       Parameters
      INCLUDE 'cst_param.finc'
!       Variables
      INTEGER :: nlat,nlon,ndep                  ! Nb of pts
      REAL    :: dlat,dlon,ddep                  ! Spacing
      REAL    :: slat,slon,sdep                  ! Start of grid
      INTEGER :: n1,n2                           ! Nb of pts
      INTEGER :: nbr                             ! Nb of fields
      REAL    :: lat0,lon0,dep0                  ! Origin of slice
      REAL    :: lat1,lon1,dep1                  ! End of slice
      CHARACTER(LEN=30) :: ftr,ftp               ! Files to read/print
      REAL    :: latr,lonr,depr                  ! Pt pos in real world
      REAL    :: latg,long,depg                  ! Pt pos on grid
!       Arrays
      REAL,ALLOCATABLE,DIMENSION(:)       :: hdrs
      REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: Emig
      REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: Epro
      REAL,ALLOCATABLE,DIMENSION(:)       :: Emax
      REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: Cpro
      CHARACTER(LEN=30),ALLOCATABLE,DIMENSION(:) :: filnams,fnams
      INTEGER,DIMENSION(8,3) :: CN               ! Closest Neighbours
      REAL,DIMENSION(7)      :: val              ! Trilinear interpol
      REAL,DIMENSION(5)      :: Ctrb
!       Counters
      INTEGER :: i,j,k,l,m,n
      INTEGER :: dum
      REAL    :: dumm


                   !********************************!
                   !*                              *!
!==================!*    Load the data and info    *!==================!
                   !*                              *!
                   !********************************!

!       Read input file
      OPEN(101,file='../Input/Slice.in',form='formatted')
      READ(101,*) !#
      READ(101,*) nlat, dlat, slat
      READ(101,*) nlon, dlon, slon
      READ(101,*) ndep, ddep, sdep
      READ(101,*) !#
      READ(101,*) n1
      READ(101,*) n2
      READ(101,*) !#
      READ(101,*) lat0
      READ(101,*) lon0
      READ(101,*) dep0
      READ(101,*) !#
      READ(101,*) lat1
      READ(101,*) lon1
      READ(101,*) dep1
      READ(101,*) !#
      READ(101,"(A20)") ftr
      READ(101,"(A20)") ftp
      CLOSE(101)

!       Read migrated files names
      OPEN(102,file=ftr,form='formatted')
      READ(102,*) nbr
      ALLOCATE(filnams(nbr),fnams(nbr),hdrs(nbr))
      DO i=1,nbr
        READ(102,"(A30)") filnams(i)
      ENDDO

!       Allocate the data
      ALLOCATE(Emig(nbr,nlon,nlat,ndep),Epro(nbr,n1,n2))
      ALLOCATE(Cpro(n1,n2,1))
      ALLOCATE(Emax(nbr))
      Emax = 0

!       Load the data from the migration
      DO i=1,nbr
        j = 200 + i
        OPEN(j,file=TRIM(filnams(i)),form='formatted')
        READ(j,*) hdrs(i)
      ENDDO
      WRITE(*,*) ''
      WRITE(*,*) '  Loading the data...'
      WRITE(*,*) ''
      WRITE(*,*) '          0',&
                 '                      50',&
                 '                      100'
      WRITE(*,"(A11)",advance='no') ''

      DO i=1,nlon
        DO j=1,nlat
          DO k=1,ndep
            DO l=1,nbr
              m = 200 + l
              READ(m,*) dum,dum,dum,Emig(l,i,j,k)
            ENDDO
          ENDDO
        ENDDO
        IF(MOD(REAL(i),REAL(nlon-1)/50).LT.1) THEN
          WRITE(*,"(A1)",advance='no') '#'
        ENDIF
      ENDDO
      WRITE(*,*) ''
      WRITE(*,*) ''

      Epro=0

      DO i=1,nbr
        j = 200 + i
        CLOSE(j)
      ENDDO


                      !**************************!
                      !*                        *!
!=====================!*    Project the data    *!==================!
                      !*                        *!
                      !**************************!

!       Project data on the plane
      WRITE(*,*) ''
      WRITE(*,*) '  Projecting the data...'
      WRITE(*,*) ''
      WRITE(*,*) '          0',&
                 '                      50',&
                 '                      100'
      WRITE(*,"(A11)",advance='no') ''
      DO i=1,n1
        DO j=1,n2

!       Find position of the point on the 2D plan
          latr = lat0 + ((i-1)/REAL(n1))*(lat1-lat0)
          lonr = lon0 + ((i-1)/REAL(n1))*(lon1-lon0)
          depr = dep0 + ((j-1)/REAL(n2))*(dep1-dep0)

!       Transform positions to grid numbers
          latg = (latr-slat)/(dlat) + 1
          IF (latg.LT.0.5) CYCLE   
          IF (latg.GT.nlat) CYCLE
          long = (lonr-slon)/(dlon) + 1
          IF (long.LT.0.5) CYCLE
          IF (long.GT.nlon) CYCLE
          depg = (sdep-depr)/(ddep) + 1
          IF (depg.LT.0.5) CYCLE
          IF (depg.GT.ndep) CYCLE

!       Find 8 neighbouring points
          CN(1,:) = (/FLOOR(long),  FLOOR(latg),  FLOOR(depg)  /) ! FFF
          CN(2,:) = (/FLOOR(long),  FLOOR(latg)+1,FLOOR(depg)  /) ! FCF
          CN(3,:) = (/FLOOR(long)+1,FLOOR(latg)+1,FLOOR(depg)  /) ! CCF
          CN(4,:) = (/FLOOR(long)+1,FLOOR(latg),  FLOOR(depg)  /) ! CFF
          CN(5,:) = (/FLOOR(long),  FLOOR(latg),  FLOOR(depg)+1/) ! FFC
          CN(6,:) = (/FLOOR(long),  FLOOR(latg)+1,FLOOR(depg)+1/) ! FCC
          CN(7,:) = (/FLOOR(long)+1,FLOOR(latg)+1,FLOOR(depg)+1/) ! CCC
          CN(8,:) = (/FLOOR(long)+1,FLOOR(latg),  FLOOR(depg)+1/) ! CFC

          IF (ANY(CN .EQ. 0)       .OR. ANY(CN(:,1).EQ.nlon).OR. &
              ANY(CN(:,2).EQ.nlat) .OR. ANY(CN(:,3).EQ.ndep)) CYCLE

!       Trilinear interpolation
          DO k=1,nbr
            val(1) = (CN(5,3)-depg)*Emig(k,CN(1,1),CN(1,2),CN(1,3)) &
                     + (depg-CN(1,3))*Emig(k,CN(5,1),CN(5,2),CN(5,3))
            val(2) = (CN(6,3)-depg)*Emig(k,CN(2,1),CN(2,2),CN(2,3)) &
                     + (depg-CN(2,3))*Emig(k,CN(6,1),CN(6,2),CN(6,3))
            val(3) = (CN(7,3)-depg)*Emig(k,CN(3,1),CN(3,2),CN(3,3)) &
                     + (depg-CN(3,3))*Emig(k,CN(7,1),CN(7,2),CN(7,3))
            val(4) = (CN(8,3)-depg)*Emig(k,CN(4,1),CN(4,2),CN(4,3)) &
                     + (depg-CN(4,3))*Emig(k,CN(8,1),CN(8,2),CN(8,3))
            val(5) = (CN(4,1)-long)*val(1) + (long-CN(1,1))*val(4)
            val(6) = (CN(3,1)-long)*val(2) + (long-CN(2,1))*val(3)
            val(7) = (CN(2,2)-latg)*val(5) + (latg-CN(1,2))*val(6)
            Epro(k,i,j) = val(7)
            Emax(k)     = MAX(ABS(Epro(k,i,j)),Emax(k))
          ENDDO

        ENDDO
        IF(MOD(REAL(i),REAL(n1)/50).LT.1) THEN
          WRITE(*,"(A1)",advance='no') '#'
        ENDIF
      ENDDO

      DO k=1,nbr
        Epro(k,:,:) = Epro(k,:,:)/Emax(k)
      ENDDO

      WRITE(*,*) ''
      WRITE(*,*) ''


                !**************************************!
                !*                                    *!
!===============!*    Compare the amplitude images    *!===============!
                !*                                    *!
                !**************************************!

      DO i=1,n1
        DO j=1,n2
          IF (ABS(Epro(1,i,j)).LT.0.1) Ctrb(1) = 0
          IF (ABS(Epro(1,i,j)).GE.0.1) Ctrb(1) = SIGN(1.,Epro(1,i,j))
          IF (ABS(Epro(2,i,j)).LT.0.1) Ctrb(2) = 0
          IF (ABS(Epro(2,i,j)).GE.0.1) Ctrb(2) = SIGN(1.,Epro(2,i,j))
          IF (ABS(Epro(3,i,j)).LT.0.1) Ctrb(3) = 0
          IF (ABS(Epro(3,i,j)).GE.0.1) Ctrb(3) = SIGN(1.,Epro(3,i,j))
          IF (ABS(Epro(4,i,j)).LT.0.1) Ctrb(4) = 0
          IF (ABS(Epro(4,i,j)).GE.0.1) Ctrb(4) = SIGN(1.,Epro(4,i,j))
          IF (ABS(Epro(5,i,j)).LT.0.1) Ctrb(5) = 0
          IF (ABS(Epro(5,i,j)).GE.0.1) Ctrb(5) = SIGN(1.,Epro(5,i,j))
          Cpro(i,j,1) = ABS( Ctrb(1)+Ctrb(2)+Ctrb(3)+Ctrb(4)+Ctrb(5) )
        ENDDO
      ENDDO


                      !***************************!
                      !*                         *!
!=====================!*    Write the outputs    *!====================!
                      !*                         *!
                      !***************************!

!       Do stuff
      OPEN(103,file=ftp,form='formatted')
      READ(103,*)
      DO i=1,nbr
        READ(103,"(A30)") fnams(i)
      ENDDO

      DO i=1,nbr
        j = 300 + i
        OPEN(j,file=fnams(i),form='formatted')
        WRITE(j,*) hdrs(i)
      ENDDO

      OPEN(400,file='../Output/Sliced/C_pro',form='formatted')
      WRITE(400,*) 5

!       Do more stuff
      WRITE(*,*) ''
      WRITE(*,*) '  Writing the outputs...'
      WRITE(*,*) ''
      WRITE(*,*) '          0',&
                 '                      50',&
                 '                      100'
      WRITE(*,"(A11)",advance='no') ''
      DO i=1,n1
        DO j=1,n2
          WRITE(400,*) i, j, Cpro(i,j,1)
          DO k=1,nbr
            l = 300 + k
            WRITE(l,*) i, j, Epro(k,i,j)
          ENDDO
        ENDDO
        IF(MOD(REAL(i),REAL(n1)/50).LT.1) THEN
          WRITE(*,"(A1)",advance='no') '#'
        ENDIF
      ENDDO
      WRITE(*,*) ''
      WRITE(*,*) ''

      DO i=1,nbr
        j = 300 + i
        CLOSE(j)
      ENDDO

      !CLOSE(400)

      DEALLOCATE(Emig)
      DEALLOCATE(Epro)
      DEALLOCATE(Cpro)


                           !*****************!
                           !*               *!
!==========================!*    The end    *!=========================!
                           !*               *!
                           !*****************!
      END
