PROGRAM MAIN
 IMPLICIT NONE
 INTEGER      , PARAMETER  :: r8=8
 INTEGER      , PARAMETER  :: r4=4
 INTEGER                   :: N
 REAL                      :: x0,k
 CHARACTER(LEN=100)        :: FileName='trajetorias'
 LOGICAL                   :: CtrlWriteDataFile=.TRUE.
 INTEGER                   :: UnitData=10
 INTEGER                   :: UnitCtl=11
 INTEGER                   :: status
 REAL, ALLOCATABLE         :: a1(:)
 REAL, ALLOCATABLE         :: a2(:)
 REAL, ALLOCATABLE         :: a3(:)

 !case  0 < k <1
 N =100 
 k =3.7
 status = SchemeWriteCtl(N+1)
 ALLOCATE(a1(0:N))
 ALLOCATE(a2(0:N))
 ALLOCATE(a3(0:N))
 x0=0.6
 a1     = recorrencia(k,x0,N)
 x0=0.60001
 a2     = recorrencia(k,x0,N)
 x0=0.59999
 a3     = recorrencia(k,x0,N)
 status = SchemeWriteData(a1,a2,a3,N)
 IF(ALLOCATED(a1))DEALLOCATE(a1)
 IF(ALLOCATED(a2))DEALLOCATE(a2)
 IF(ALLOCATED(a3))DEALLOCATE(a3)



 STOP
CONTAINS 
 FUNCTION recorrencia(k,x0,N)  RESULT (x_vec)
  IMPLICIT NONE
  !
  !  x(n) = k*x(n-1)*(1-x(n-1))
  !
  REAL   , INTENT(IN   ) :: k
  REAL   , INTENT(IN   ) :: x0
  INTEGER, INTENT(IN   ) :: N
  REAL                   :: x_vec(0:N)
  INTEGER                :: i
  x_vec(0)=x0
  DO i=1,N
    x_vec(i) = k*x_vec(i-1)*(1-x_vec(i-1))
  END DO
  
 END FUNCTION recorrencia


 FUNCTION SchemeWriteData(vars1,vars2,vars3,ndata)  RESULT (ok)
    IMPLICIT NONE
    INTEGER       , INTENT (IN   ) :: ndata
    REAL (KIND=r4), INTENT (INOUT) :: vars1(0:ndata)
    REAL (KIND=r4), INTENT (INOUT) :: vars2(0:ndata)
    REAL (KIND=r4), INTENT (INOUT) :: vars3(0:ndata)
    INTEGER        :: ok
    INTEGER        :: i
    INTEGER        :: irec
    INTEGER        :: lrec
    REAL (KIND=r4) :: Yout(0:ndata)
    IF(CtrlWriteDataFile)INQUIRE (IOLENGTH=lrec) Yout(0)
    IF(CtrlWriteDataFile)OPEN(UnitData,FILE=TRIM(FileName)//'.bin',&
    FORM='UNFORMATTED', ACCESS='DIRECT', STATUS='UNKNOWN', &
    ACTION='WRITE',RECL=lrec)
    ok=1
    irec=0
    CtrlWriteDataFile=.FALSE.
    DO i=0,ndata
       Yout(i)=REAL(vars1(i),KIND=r4)
       irec=irec+1
       WRITE(UnitData,rec=irec)Yout(i)
       Yout(i)=REAL(vars2(i),KIND=r4)
       irec=irec+1
       WRITE(UnitData,rec=irec)Yout(i)
       Yout(i)=REAL(vars3(i),KIND=r4)
       irec=irec+1
       WRITE(UnitData,rec=irec)Yout(i)
    END DO
    ok=0
 END FUNCTION SchemeWriteData

 FUNCTION SchemeWriteCtl(nrec)  RESULT (ok)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nrec
    INTEGER              :: ok,i
    INTEGER              :: Idim=1
    ok=1
   OPEN(UnitCtl,FILE=TRIM(FileName)//'.ctl',FORM='FORMATTED', &
   ACCESS='SEQUENTIAL',STATUS='UNKNOWN',ACTION='WRITE')
    WRITE (UnitCtl,'(A6,A           )')'dset ^',TRIM(FileName)//'.bin'
    WRITE (UnitCtl,'(A                 )')'title  EDO'
    WRITE (UnitCtl,'(A                 )')'undef  -9999.9'
    WRITE (UnitCtl,'(A6,I8,A12   )')'xdef  ',Idim,' levels 0   '
    WRITE (UnitCtl,'(A6,I8,A12   )')'ydef  ',Idim,' levels 0   '
    WRITE (UnitCtl,'(A                  )')'ydef  1 linear  -1.27 1'
    WRITE (UnitCtl,'(A6,I6,A25   )')'tdef  ',nrec,' linear  00z01jan0001 1hr'
    WRITE (UnitCtl,'(A20             )')'zdef  1 levels 1000 '
    WRITE (UnitCtl,'(A           )')'vars 3'
    WRITE (UnitCtl,'(A           )')'x1 0 99 trajetoria da 1 '
    WRITE (UnitCtl,'(A           )')'x2 0 99 trajetoria da 2 '
    WRITE (UnitCtl,'(A           )')'x3 0 99 trajetoria da 3 '
    WRITE (UnitCtl,'(A           )')'endvars'
    CLOSE (UnitCtl,STATUS='KEEP')
    ok=0
 END FUNCTION SchemeWriteCtl
END PROGRAM MAIN 
