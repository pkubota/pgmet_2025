!  $Author: pkubota 						$
!  $Date: 2008/09/23 17:51:54 					$
!  $Revision: 1.9 						$
!  $Revisions are currently made by the class's students.	$
!  $Update Date: 01/11/2023 10:19 AM				$
!  
!  Implementações: 
!  	1) Colocar INIT e FINALIZE - OK
! 	7) Amortecimento das condições de contorno - OK
!	8) Colocar todos os campos em "WRITE FIELDS" e salvar um bin - OK

MODULE Class_Module_Fields
 USE Constants, Only: r8, r4, i4, pi, Deg2Rad, r_earth, omega, nfprt,unit_namelist

  IMPLICIT NONE
  PRIVATE       

  ! Selecting Unit

  INTEGER,PUBLIC, PARAMETER :: unitsurp=50
  INTEGER,PUBLIC, PARAMETER :: unituvel=51
  INTEGER,PUBLIC, PARAMETER :: unitvvel=52
  INTEGER,PUBLIC, PARAMETER :: unitomeg=56
  INTEGER,PUBLIC, PARAMETER :: unittemp=53
  INTEGER,PUBLIC, PARAMETER :: unitumes=54
  INTEGER,PUBLIC, PARAMETER :: unitzgeo=55
  INTEGER,PUBLIC, PARAMETER :: unitoutp=60
  
  INTEGER,PUBLIC, PARAMETER :: unitumodel=61  
  INTEGER,PUBLIC, PARAMETER :: unittmodel=62
  INTEGER,PUBLIC, PARAMETER :: unitmmodel=63
  INTEGER,PUBLIC, PARAMETER :: unitqmodel=64    

  INTEGER                    :: irec_local
  INTEGER,PUBLIC             :: nLon     =  93!161
  INTEGER,PUBLIC             :: nLat     =  61!161
  REAL(KIND=r4)              :: InitLon  = 302.0 !0  - 360 
  REAL(KIND=r4)              :: InitLat  = -35.0 !-90   90 
  REAL(KIND=r8)              :: DeltaLon =  0.25
  REAL(KIND=r8)              :: DeltaLat =  0.25
  INTEGER,PUBLIC             :: nLev     =37
  REAL(KIND=8), PUBLIC       :: vis      =  2.5e3_r8! 1.5e-10 !1.5e-5        !  viscosity
  REAL(KIND=8), PUBLIC       :: vis_q    =  3.0e4_r8! 1.5e-10 !1.5e-5        !  viscosity
  REAL(KIND=8), PUBLIC       :: taul     =  1800_r8

  INTEGER                    :: lrec2D
  INTEGER                    :: lrec3D

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: Plevs(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: CoordLat(:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: CoordLon(:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: FcorPar (:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: DeltaLamda(:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: DeltaTheta(:,:) 

  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: u_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: v_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: w_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: t_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: q_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: z_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: p_ref(:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var2P_A(:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var2P_B(:,:) 
  
  REAL(KIND=r4), ALLOCATABLE :: var3U_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3U_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3V_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3V_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3W_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3W_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3T_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3T_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3Q_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3Q_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3Z_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3Z_B(:,:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: U_N(:,:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: U_C(:,:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: V_N(:,:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: V_C(:,:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: T_N(:,:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: T_C(:,:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: Q_N(:,:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: Q_C(:,:,:) 
  
!!!!!!!!!!!!!!!!!SUBROUTINES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_LWRH(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_SWRH(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_LHCV(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_MSCV(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_LGLH(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_LGMS(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_SCVH(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_SCVM(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_PBLT(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_PBLQ(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_GDTZ(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: r_GDTM(:)
   
  PUBLIC :: Init_Class_Module_Fields, ReadFields, WriteFields, Finalize_Class_Module_Fields
  
CONTAINS

 SUBROUTINE Init_Class_Module_Fields() !INIT
  IMPLICIT NONE
  INTEGER :: i,j,k
  REAL(KIND=r8), ALLOCATABLE :: r_q1(:,:)
  REAL(KIND=r8), ALLOCATABLE :: r_q2(:,:)
  INTEGER :: ForcLW
  INTEGER :: ForcSW
  INTEGER :: ForcDC
  INTEGER :: ForcLG
  INTEGER :: ForcSC
  INTEGER :: ForcPB
  NAMELIST /ForcingControl/ForcLW,ForcSW,ForcDC,ForcLG,ForcSC,ForcPB
  
  NAMELIST /MeshesControl/ nLon,nLat,InitLon,InitLat,DeltaLon,DeltaLat,nLev,&
                           vis,vis_q,taul 
  rewind(unit_namelist)
  read(unit_namelist, nml=MeshesControl)
  rewind(unit_namelist)
  read(unit_namelist, nml=ForcingControl)
  
  ALLOCATE(Plevs(nLev));Plevs=0.0
  ALLOCATE(CoordLat(nLon,nLat));CoordLat=0.0
  ALLOCATE(CoordLon(nLon,nLat));CoordLon=0.0
  ALLOCATE(FcorPar(nLon,nLat));FcorPar=0.0
  ALLOCATE(DeltaLamda(nLon,nLat));DeltaLamda=0.0
  ALLOCATE(DeltaTheta(nLon,nLat));DeltaTheta=0.0
  
  ALLOCATE(var2P_A(nLon,nLat));var2P_A=0.0
  ALLOCATE(var2P_B(nLon,nLat));var2P_B=0.0

  ALLOCATE(var3U_A(nLon,nLat,nLev));var3U_A=0.0
  ALLOCATE(var3U_B(nLon,nLat,nLev));var3U_B=0.0

  ALLOCATE(var3V_A(nLon,nLat,nLev));var3V_A=0.0
  ALLOCATE(var3V_B(nLon,nLat,nLev));var3V_B=0.0

  ALLOCATE(var3W_A(nLon,nLat,nLev));var3W_A=0.0
  ALLOCATE(var3W_B(nLon,nLat,nLev));var3W_B=0.0

  ALLOCATE(var3T_A(nLon,nLat,nLev));var3T_A=0.0
  ALLOCATE(var3T_B(nLon,nLat,nLev));var3T_B=0.0

  ALLOCATE(var3Q_A(nLon,nLat,nLev));var3Q_A=0.0
  ALLOCATE(var3Q_B(nLon,nLat,nLev));var3Q_B=0.0

  ALLOCATE(var3Z_A(nLon,nLat,nLev));var3Z_A=0.0
  ALLOCATE(var3Z_B(nLon,nLat,nLev));var3Z_B=0.0

  ALLOCATE(U_N(nLon,nLat,nLev));U_N=0.0
  ALLOCATE(U_C(nLon,nLat,nLev));U_C=0.0

  ALLOCATE(V_N(nLon,nLat,nLev));V_N=0.0
  ALLOCATE(V_C(nLon,nLat,nLev));V_C=0.0

  ALLOCATE(T_N(nLon,nLat,nLev));T_N=0.0
  ALLOCATE(T_C(nLon,nLat,nLev));T_C=0.0

  ALLOCATE(Q_N(nLon,nLat,nLev));Q_N=0.0
  ALLOCATE(Q_C(nLon,nLat,nLev));Q_C=0.0
  
  ALLOCATE(p_ref(nLon,nLat));p_ref=0.0  
  ALLOCATE(u_ref(nLon,nLat,nLev));u_ref=0.0
  ALLOCATE(v_ref(nLon,nLat,nLev));v_ref=0.0
  ALLOCATE(w_ref(nLon,nLat,nLev));w_ref=0.0
  ALLOCATE(t_ref(nLon,nLat,nLev));t_ref=0.0
  ALLOCATE(q_ref(nLon,nLat,nLev));q_ref=0.0
  ALLOCATE(z_ref(nLon,nLat,nLev));z_ref=0.0

  ALLOCATE(r_LWRH(1:nLev));r_LWRH=0.0_r8
  ALLOCATE(r_SWRH(1:nLev));r_SWRH=0.0_r8
  ALLOCATE(r_LHCV(1:nLev));r_LHCV=0.0_r8
  ALLOCATE(r_MSCV(1:nLev));r_MSCV=0.0_r8
  ALLOCATE(r_LGLH(1:nLev));r_LGLH=0.0_r8
  ALLOCATE(r_LGMS(1:nLev));r_LGMS=0.0_r8
  ALLOCATE(r_SCVH(1:nLev));r_SCVH=0.0_r8
  ALLOCATE(r_SCVM(1:nLev));r_SCVM=0.0_r8
  ALLOCATE(r_PBLT(1:nLev));r_PBLT=0.0_r8
  ALLOCATE(r_PBLQ(1:nLev));r_PBLQ=0.0_r8
  ALLOCATE(r_GDTZ(1:nLev));r_GDTZ=0.0_r8
  ALLOCATE(r_GDTM(1:nLev));r_GDTM=0.0_r8
 
  INQUIRE(IOLENGTH=lrec2D)var2P_A
  INQUIRE(IOLENGTH=lrec3D)var3U_A

  OPEN(unit=unitzgeo,FILE='../datain/GeoPotential.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unittemp,FILE='../datain/Temperature.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitumes,FILE='../datain/SpecificHumidy.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unituvel,FILE='../datain/ZonalWind.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitvvel,FILE='../datain/MeridionalWind.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitomeg,FILE='../datain/Omega.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitsurp,FILE='../datain/SurfacePressure.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec2D,ACTION='READ',STATUS='OLD') 
       
!SAIDAS DO MODELO SAO SALVAS EM:
  OPEN(unit=unitumodel,FILE='../dataout/ModelU.bin',&
      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='WRITE',STATUS='UNKNOWN')
      
  OPEN(unit=unittmodel,FILE='../dataout/ModelT.bin',&
      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='WRITE',STATUS='UNKNOWN') 
      
  OPEN(unit=unitmmodel,FILE='../dataout/ModelM.bin',&
      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='WRITE',STATUS='UNKNOWN') 
      
!  OPEN(unit=unitqmodel,FILE='ModelQ.bin',&
!      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec2D,ACTION='WRITE',STATUS='UNKNOWN')   
  OPEN(unit=unitqmodel,FILE='../dataout/ModelQ.bin',&
      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='WRITE',STATUS='UNKNOWN')   

!!!!!!!!!!!!!!!!!!!!DEFINICAO DA MALHA (QUADRADA)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO j=1,nLat
     CoordLon(1,j) = InitLon*Deg2Rad
     DO i=2,nLon
        CoordLon(i,j) = CoordLon(i-1,j) + (DeltaLon*Deg2Rad)
     END DO
  END DO
  
  DO i=1,nLon
     CoordLat(i,1) = InitLat*Deg2Rad
     DO j=2,nLat
        CoordLat(i,j) = CoordLat(i,j-1) + (DeltaLat*Deg2Rad)
     END DO
  END DO

  DO j=1,nLat
     DO i=1,nLon
        FcorPar(i,j) = 2.0_r8*omega*sin(CoordLat(i,j))
     END DO
  END DO

  DO j=1,nLat
     DO i=2,nLon-1
        DeltaLamda(i,j) = (CoordLon(i+1,j)-CoordLon(i,j))
     END DO
        DeltaLamda(1,j)    = (DeltaLamda(2,j))
        DeltaLamda(nLon,j) = (DeltaLamda(nLon-1,j))
  END DO

  DO i=1,nLon
     DO j=2,nLat-1
        DeltaTheta(i,j) = (CoordLat(i,j+1)-CoordLat(i,j))
     END DO
        DeltaTheta(i,1) = DeltaTheta(i,2)
        DeltaTheta(i,nLat) = DeltaTheta(i,nLat-1)
  END DO
  Plevs(1:nLev)=(/1000.0_r8,975.0_r8,950.0_r8,925.0_r8,900.0_r8,875.0_r8,850.0_r8,825.0_r8,800.0_r8,775.0_r8,&
                   750.0_r8,700.0_r8,650.0_r8,600.0_r8,550.0_r8,500.0_r8,450.0_r8,400.0_r8,350.0_r8,300.0_r8,&
                   250.0_r8,225.0_r8,200.0_r8,175.0_r8,150.0_r8,125.0_r8,100.0_r8, 70.0_r8, 50.0_r8, 30.0_r8,&
                    20.0_r8, 10.0_r8,  7.0_r8,  5.0_r8,  3.0_r8,  2.0_r8,  1.0_r8/)             

  Plevs=Plevs*100.0_r8
  
  ALLOCATE(r_q1(1:nLev,1:6));r_q1=0.0_r8
  ALLOCATE(r_q2(1:nLev,1:4));r_q2=0.0_r8
  r_q1 = RESHAPE( (/&
! r_LWRH(k) , r_SWRH(k) , r_LHCV(k) , r_LGLH(k) , r_SCVH(k) , r_PBLT(k)
  0.1518E+00, 0.2005E-01, 0.7377E-01, 0.3375E+00, 0.0000E+00, 0.4168E+00,&
  0.3944E+00, 0.3373E-01, 0.9684E-01, 0.4621E+00, 0.0000E+00, 0.1288E-01,&
  0.4087E+00, 0.3984E-01, 0.7604E-01, 0.3658E+00, 0.0000E+00, 0.1097E+00,&
  0.3974E+00, 0.5922E-01, 0.5273E-01, 0.3819E+00, 0.0000E+00, 0.1088E+00,&
  0.3659E+00, 0.1017E+00, 0.1915E-01, 0.4840E+00, 0.0000E+00, 0.2925E-01,&
  0.2421E+00, 0.1111E+00, 0.8334E-01, 0.4612E+00, 0.0000E+00, 0.1023E+00,&
  0.2802E+00, 0.1005E+00, 0.1017E+00, 0.4156E+00, 0.0000E+00, 0.1020E+00,&
  0.3296E+00, 0.1104E+00, 0.1368E+00, 0.3708E+00, 0.0000E+00, 0.5244E-01,&
  0.3555E+00, 0.1334E+00, 0.1992E+00, 0.2895E+00, 0.0000E+00, 0.2246E-01,&
  0.3341E+00, 0.1353E+00, 0.2531E+00, 0.1852E+00, 0.0000E+00, 0.9229E-01,&
  0.2819E+00, 0.1469E+00, 0.3302E+00, 0.1146E+00, 0.0000E+00, 0.1262E+00,&
  0.2384E+00, 0.1583E+00, 0.3807E+00, 0.8140E-01, 0.0000E+00, 0.1412E+00,&
  0.3104E+00, 0.1814E+00, 0.3087E+00, 0.5525E-01, 0.0000E+00, 0.1442E+00,&
  0.3835E+00, 0.1830E+00, 0.2917E+00, 0.4227E-01, 0.0000E+00, 0.9952E-01,&
  0.4805E+00, 0.1319E+00, 0.2904E+00, 0.9293E-01, 0.0000E+00, 0.4252E-02,&
  0.4605E+00, 0.1158E+00, 0.2624E+00, 0.9486E-01, 0.0000E+00, 0.6647E-01,&
  0.4074E+00, 0.1508E+00, 0.2592E+00, 0.5118E-01, 0.0000E+00, 0.1314E+00,&
  0.4332E+00, 0.1918E+00, 0.2488E+00, 0.8094E-03, 0.0000E+00, 0.1254E+00,&
  0.4372E+00, 0.2039E+00, 0.2236E+00, 0.6251E-01, 0.0000E+00, 0.7277E-01,&
  0.4443E+00, 0.2136E+00, 0.2011E+00, 0.1061E+00, 0.0000E+00, 0.3488E-01,&
  0.4183E+00, 0.2412E+00, 0.1779E+00, 0.1462E+00, 0.0000E+00, 0.1641E-01,&
  0.4319E+00, 0.2616E+00, 0.1544E+00, 0.1307E+00, 0.0000E+00, 0.2141E-01,&
  0.4282E+00, 0.3202E+00, 0.1233E+00, 0.1000E+00, 0.0000E+00, 0.2822E-01,&
  0.5115E+00, 0.3717E+00, 0.5096E-01, 0.2992E-01, 0.0000E+00, 0.3592E-01,&
  0.5525E+00, 0.3956E+00, 0.1254E-01, 0.3796E-02, 0.0000E+00, 0.3555E-01,&
  0.4865E+00, 0.4155E+00, 0.4267E-02, 0.4971E-01, 0.0000E+00, 0.4408E-01,&
  0.3521E+00, 0.4688E+00, 0.2058E-02, 0.1043E+00, 0.0000E+00, 0.7271E-01,&
  0.2907E+00, 0.5950E+00, 0.6246E-03, 0.8891E-02, 0.0000E+00, 0.1047E+00,&
  0.2899E+00, 0.5899E+00, 0.1449E-03, 0.5083E-01, 0.0000E+00, 0.6922E-01,&
  0.3086E+00, 0.6526E+00, 0.3578E-04, 0.3716E-01, 0.0000E+00, 0.1602E-02,&
  0.2737E+00, 0.6176E+00, 0.8161E-05, 0.1343E-01, 0.0000E+00, 0.9529E-01,&
  0.2157E+00, 0.6167E+00, 0.3623E-05, 0.7204E-02, 0.0000E+00, 0.1604E+00,&
  0.1152E+00, 0.7213E+00, 0.3406E-05, 0.6771E-02, 0.0000E+00, 0.1568E+00,&
  0.7200E-01, 0.8329E+00, 0.1955E-05, 0.3887E-02, 0.0000E+00, 0.9124E-01,&
  0.2202E+00, 0.7549E+00, 0.5127E-06, 0.1019E-02, 0.0000E+00, 0.2393E-01,&
  0.3547E+00, 0.6398E+00, 0.1113E-06, 0.2213E-03, 0.0000E+00, 0.5194E-02,&
  0.3974E+00, 0.6004E+00, 0.4616E-07, 0.9177E-04, 0.0000E+00, 0.2154E-02/),(/nLev,6/), order=(/2,1/))

  DO k=1,nLev
    r_LWRH(k) = r_q1(k,1)
    r_SWRH(k) = r_q1(k,2)
    r_LHCV(k) = r_q1(k,3)
    r_LGLH(k) = r_q1(k,4)
    r_SCVH(k) = r_q1(k,5)
    r_PBLT(k) = r_q1(k,6)    
  END DO 
  IF(ForcLW == 0)r_LWRH=0.0_r8
  IF(ForcSW == 0)r_SWRH=0.0_r8
  IF(ForcDC == 0)r_LHCV=0.0_r8
  IF(ForcLG == 0)r_LGLH=0.0_r8
  IF(ForcSC == 0)r_SCVH=0.0_r8
  IF(ForcPB == 0)r_PBLT=0.0_r8

  !WRITE(1,'(A)')'MSCV   37 99 CONVECTIVE MOISTURE SOURCE  	    (1/Sec	     )'
  !WRITE(1,'(A)')'LGMS   37 99 LARGE SCALE MOISTURE SOURCE 	    (1/Sec	     )'
  !WRITE(1,'(A)')'SCVM   37 99 SHALLOW CONV. MOISTURE SOURCE	    (1/Sec	     )'
  !WRITE(1,'(A)')'PBLQ   37 99 VERTICAL DIFF. MOISTURE SOURCE	    (1/Sec	     )'
  r_q2 = RESHAPE( (/&
!   r_SCVM(k) , r_MSCV(k) , r_LGMS(k) , r_PBLQ(k)
    0.0000E+00, 0.5001E-01, 0.3444E+00, 0.6056E+00,&
    0.0000E+00, 0.4076E-01, 0.3377E+00, 0.6216E+00,&
    0.0000E+00, 0.1077E+00, 0.2979E+00, 0.5944E+00,&
    0.0000E+00, 0.2034E+00, 0.2599E+00, 0.5367E+00,&
    0.0000E+00, 0.3024E+00, 0.2334E+00, 0.4642E+00,&
    0.0000E+00, 0.3415E+00, 0.2815E+00, 0.3770E+00,&
    0.0000E+00, 0.3187E+00, 0.3970E+00, 0.2843E+00,&
    0.0000E+00, 0.3796E+00, 0.5227E+00, 0.9766E-01,&
    0.0000E+00, 0.4690E+00, 0.3905E+00, 0.1404E+00,&
    0.0000E+00, 0.5336E+00, 0.2364E+00, 0.2300E+00,&
    0.0000E+00, 0.5836E+00, 0.1273E+00, 0.2892E+00,&
    0.0000E+00, 0.6843E+00, 0.8767E-01, 0.2281E+00,&
    0.0000E+00, 0.6170E+00, 0.9126E-01, 0.2917E+00,&
    0.0000E+00, 0.5496E+00, 0.1241E+00, 0.3263E+00,&
    0.0000E+00, 0.3128E+00, 0.3334E+00, 0.3537E+00,&
    0.0000E+00, 0.2090E+00, 0.2666E+00, 0.5245E+00,&
    0.0000E+00, 0.3812E+00, 0.1228E+00, 0.4960E+00,&
    0.0000E+00, 0.5039E+00, 0.2374E-02, 0.4938E+00,&
    0.0000E+00, 0.4381E+00, 0.2458E+00, 0.3160E+00,&
    0.0000E+00, 0.3347E+00, 0.4963E+00, 0.1691E+00,&
    0.0000E+00, 0.3566E+00, 0.5819E+00, 0.6147E-01,&
    0.0000E+00, 0.3370E+00, 0.6500E+00, 0.1303E-01,&
    0.0000E+00, 0.2058E+00, 0.6570E+00, 0.1373E+00,&
    0.0000E+00, 0.3440E+00, 0.5462E+00, 0.1098E+00,&
    0.0000E+00, 0.6231E+00, 0.2409E+00, 0.1360E+00,&
    0.0000E+00, 0.2406E-01, 0.4686E+00, 0.5074E+00,&
    0.0000E+00, 0.5444E-02, 0.4591E+00, 0.5355E+00,&
    0.0000E+00, 0.2350E-01, 0.5590E+00, 0.4175E+00,&
    0.0000E+00, 0.6233E-03, 0.3628E+00, 0.6365E+00,&
    0.0000E+00, 0.2288E-03, 0.3945E+00, 0.6053E+00,&
    0.0000E+00, 0.1483E-03, 0.4050E+00, 0.5948E+00,&
    0.0000E+00, 0.1252E-03, 0.4133E+00, 0.5866E+00,&
    0.0000E+00, 0.1248E-03, 0.4119E+00, 0.5880E+00,&
    0.0000E+00, 0.1248E-03, 0.4119E+00, 0.5880E+00,&
    0.0000E+00, 0.1248E-03, 0.4119E+00, 0.5880E+00,&
    0.0000E+00, 0.1248E-03, 0.4119E+00, 0.5880E+00,&
    0.0000E+00, 0.1248E-03, 0.4119E+00, 0.5880E+00/),(/nLev,4/), order=(/2,1/))
  DO k=1,nLev
    r_SCVM(k) = r_q2(k,1) 
    r_MSCV(k) = r_q2(k,2)
    r_LGMS(k) = r_q2(k,3)
    r_PBLQ(k) = r_q2(k,4)
  END DO 
  IF(ForcDC == 0)r_MSCV=0.0
  IF(ForcLG == 0)r_LGMS=0.0
  IF(ForcSC == 0)r_SCVM=0.0
  IF(ForcPB == 0)r_PBLQ=0.0
 END SUBROUTINE Init_Class_Module_Fields

 SUBROUTINE ReadFields(irec,TimeIncrSeg) 
  IMPLICIT NONE 
  INTEGER      , INTENT(IN   ) :: irec
  REAL(KIND=r8), INTENT(IN   ) :: TimeIncrSeg
  REAL(KIND=r8) :: w1
  REAL(KIND=r8) :: w2
  INTEGER       :: i,j,k

  w1=1.0_r8-mod(TimeIncrSeg,3600.0_r8)/3600.0_r8
  w2=mod(TimeIncrSeg,3600.0_r8)/3600.0_r8

  IF(irec > irec_local)THEN
    READ(unituvel,rec=irec)var3U_A
    READ(unitvvel,rec=irec)var3V_A
    READ(unitomeg,rec=irec)var3W_A
    READ(unittemp,rec=irec)var3T_A
    READ(unitumes,rec=irec)var3Q_A
    READ(unitzgeo,rec=irec)var3Z_A
    READ(unitsurp,rec=irec)var2P_A

    READ(unituvel,rec=irec+1)var3U_B
    READ(unitvvel,rec=irec+1)var3V_B
    READ(unitomeg,rec=irec+1)var3W_B
    READ(unittemp,rec=irec+1)var3T_B
    READ(unitumes,rec=irec+1)var3Q_B
    READ(unitzgeo,rec=irec+1)var3Z_B
    READ(unitsurp,rec=irec+1)var2P_B
    irec_local=irec
  END IF

  p_ref= var2P_A *w1 + w2*var2P_B
  u_ref= var3U_A *w1 + w2*var3U_B
  v_ref= var3V_A *w1 + w2*var3V_B
  w_ref= var3W_A *w1 + w2*var3W_B
  t_ref= var3T_A *w1 + w2*var3T_B
  q_ref= var3Q_A *w1 + w2*var3Q_B
  z_ref= var3Z_A *w1 + w2*var3Z_B

  DO k=1,nLev

     z_ref(1   ,1:nLat,k)=0.25_r8*(z_ref(1   ,1:nLat,k) + z_ref(2     ,1:nLat,k) + z_ref(3     ,1:nLat,k) + z_ref(2     ,1:nLat,k) )
     z_ref(nLon,1:nLat,k)=0.25_r8*(z_ref(nLon,1:nLat,k) + z_ref(nLon-1,1:nLat,k) + z_ref(nLon-2,1:nLat,k) + z_ref(nLon-1,1:nLat,k) )

     z_ref(1:nLon,1   ,k)=0.25_r8*(z_ref(1:nLon,1   ,k) + z_ref(1:nLon,2     ,k) + z_ref(1:nLon,3     ,k) + z_ref(1:nLon,2     ,k) )
     z_ref(1:nLon,nLat,k)=0.25_r8*(z_ref(1:nLon,nLat,k) + z_ref(1:nLon,nLat-1,k) + z_ref(1:nLon,nLat-2,k) + z_ref(1:nLon,nLat-1,k) )

     w_ref(1   ,1:nLat,k)=0.25_r8*(w_ref(1   ,1:nLat,k) + w_ref(2     ,1:nLat,k) + w_ref(3     ,1:nLat,k) + w_ref(2     ,1:nLat,k) )
     w_ref(nLon,1:nLat,k)=0.25_r8*(w_ref(nLon,1:nLat,k) + w_ref(nLon-1,1:nLat,k) + w_ref(nLon-2,1:nLat,k) + w_ref(nLon-1,1:nLat,k) )

     w_ref(1:nLon,1   ,k)=0.25_r8*(w_ref(1:nLon,1   ,k) + w_ref(1:nLon,2     ,k) + w_ref(1:nLon,3     ,k) + w_ref(1:nLon,2     ,k) )
     w_ref(1:nLon,nLat,k)=0.25_r8*(w_ref(1:nLon,nLat,k) + w_ref(1:nLon,nLat-1,k) + w_ref(1:nLon,nLat-2,k) + w_ref(1:nLon,nLat-1,k) )
  END DO
  DO k=1,nLev
     DO j=2,nLat-1
        DO i=2,nLon-1
           z_ref(i   ,j,k)=(1.0_r8/5.0_r8)*(z_ref(i+1,j,k) + z_ref(i,j,k) + z_ref(i-1,j,k)+ z_ref(i,j+1,k)+ z_ref(i,j-1,k))
           w_ref(i   ,j,k)=(1.0_r8/5.0_r8)*(w_ref(i+1,j,k) + w_ref(i,j,k) + w_ref(i-1,j,k)+ w_ref(i,j+1,k)+ w_ref(i,j-1,k))
        END DO
     END DO
  END DO

  DO k=1,nLev-1
     u_ref(:,:,k)= 0.5_r8*(u_ref(:,:,k+1)+u_ref(:,:,k))
     v_ref(:,:,k)= 0.5_r8*(v_ref(:,:,k+1)+v_ref(:,:,k))
     w_ref(:,:,k)= 0.5_r8*(w_ref(:,:,k+1)+w_ref(:,:,k))
     t_ref(:,:,k)= 0.5_r8*(t_ref(:,:,k+1)+t_ref(:,:,k))
     q_ref(:,:,k)= 0.5_r8*(q_ref(:,:,k+1)+q_ref(:,:,k))
     z_ref(:,:,k)= 0.5_r8*(z_ref(:,:,k+1)+z_ref(:,:,k))
  END DO
 END SUBROUTINE ReadFields
 
 ! NAO ENTENDI A LOGICA MAS VOU SEGUI-LA
 
 SUBROUTINE WriteFields(irec)
  IMPLICIT NONE 
  INTEGER      , INTENT(INOUT) :: irec
  
  irec=irec+1
  

  WRITE(unittmodel,rec=irec)REAL(T_C,kind=r4) !Temperature
  WRITE(unitqmodel,rec=irec)REAL(Q_C,kind=r4) !SpecificHumidy
  WRITE(unitumodel,rec=irec)REAL(U_C,kind=r4) !ZonalWind
  WRITE(unitmmodel,rec=irec)REAL(V_C,kind=r4) !MeridionalWind

    
 END SUBROUTINE WriteFields

 SUBROUTINE Finalize_Class_Module_Fields() !FINALIZE
  IMPLICIT NONE 

  DEALLOCATE(Plevs)
  DEALLOCATE(CoordLat)
  DEALLOCATE(CoordLon)
  DEALLOCATE(FcorPar)
  DEALLOCATE(DeltaLamda)
  DEALLOCATE(DeltaTheta)
  DEALLOCATE(var2P_A)
  DEALLOCATE(var2P_B)
  DEALLOCATE(var3U_A)
  DEALLOCATE(var3U_B)
  DEALLOCATE(var3V_A)
  DEALLOCATE(var3V_B)
  DEALLOCATE(var3W_A)
  DEALLOCATE(var3W_B)
  DEALLOCATE(var3T_A)
  DEALLOCATE(var3T_B)
  DEALLOCATE(var3Q_A)
  DEALLOCATE(var3Q_B)
  DEALLOCATE(var3Z_A)
  DEALLOCATE(var3Z_B)
  DEALLOCATE(U_N)   !tempo futuro
  DEALLOCATE(U_C)
  DEALLOCATE(V_N)
  DEALLOCATE(V_C)
  DEALLOCATE(T_N)
  DEALLOCATE(T_C)
  DEALLOCATE(Q_N)
  DEALLOCATE(Q_C)
  DEALLOCATE(p_ref)  
  DEALLOCATE(u_ref)
  DEALLOCATE(v_ref)
  DEALLOCATE(w_ref)
  DEALLOCATE(t_ref)
  DEALLOCATE(q_ref)
  DEALLOCATE(z_ref)

  CLOSE(unit=unitzgeo, STATUS='KEEP')
  CLOSE(unit=unittemp, STATUS='KEEP')
  CLOSE(unit=unitumes, STATUS='KEEP')
  CLOSE(unit=unituvel, STATUS='KEEP')
  CLOSE(unit=unitvvel, STATUS='KEEP')
  CLOSE(unit=unitomeg, STATUS='KEEP')
  CLOSE(unit=unitsurp, STATUS='KEEP')
  CLOSE(unit=unitoutp, STATUS='KEEP')
  
  CLOSE(unit=unitumodel, STATUS='KEEP')
  CLOSE(unit=unittmodel, STATUS='KEEP')
  CLOSE(unit=unitmmodel, STATUS='KEEP')
  CLOSE(unit=unitqmodel, STATUS='KEEP')
 


 END SUBROUTINE Finalize_Class_Module_Fields

END MODULE Class_Module_Fields
