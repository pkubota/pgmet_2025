!  $Author: pkubota 						$
!  $Date: 2008/09/23 17:51:54 					$
!  $Revision: 1.9 						$
!  $Revisions are currently made by the class's students.	$
!  $Update Date: 01/11/2023 10:19 AM				$
!  
!  ImplementaÃ§Ãµes: 
!  	1) Colocar INIT e FINALIZE - OK 
! 	2) Colocar as equaÃ§Ãµes dos campos no "Class_Module_Dynamics" - OK
! 	7) Amortecimento das condiÃ§Ãµes de contorno - OK (funcionando parcialmente)

MODULE Class_Module_Dynamics
 USE Constants, Only: r8,r4,i4,pi,Deg2Rad,Rd,Cp,kappa,r_earth,omega,CTv,nfprt, Eps,g_cte
  
 USE Class_Module_LapaceOperator, Only:sor_solve,&
                               multigrid_solver
 
 
 USE Class_Module_Fields, Only: u_ref, v_ref,w_ref, t_ref,q_ref,z_ref,p_ref,&
                                U_N,U_C,  V_N,V_C,  W_N,W_C, T_N,T_C,  Q_N,Q_C,  Z_N,Z_C,&
                                Plevs,CoordLat,CoordLon,FcorPar ,DeltaLamda,DeltaTheta,&
                                termA_QG,termB_QG,termC_QG,termD_QG,termQ_QG,termS_QG,total_forcing,&
                                thetas ,ug_surf,vg_surf,ug_iLND,vg_iLND,kMinLev
  IMPLICIT NONE
  PRIVATE       

  INTEGER ::  Idim
  INTEGER ::  Jdim
  INTEGER ::  Kdim
  REAL(KIND=r8) :: DeltaT
  REAL(KIND=r8), PUBLIC   ,parameter        :: vis    =  2.0e2_r8! 1.5e-10 !1.5e-5        !  viscosity
  REAL(KIND=r8), PUBLIC   ,parameter        :: taul   = 3600_r8
  ! Selecting Unit
  PUBLIC :: Init_Class_Module_Dynamics, Finalize_Class_Module_Dynamics, RunDynamics
  
CONTAINS

 SUBROUTINE Init_Class_Module_Dynamics (nLat,nLon,nLev,dt_step) !INIT
  IMPLICIT NONE
  INTEGER , INTENT(IN   ) :: nLat
  INTEGER , INTENT(IN   ) :: nLon
  INTEGER , INTENT(IN   ) :: nLev
  REAL(KIND=r8), INTENT(IN   )  :: dt_step
  Idim=nLon
  Jdim=nLat
  Kdim=nLev
  DeltaT=dt_step

 END SUBROUTINE Init_Class_Module_Dynamics

 SUBROUTINE RunDynamics (it) 
  IMPLICIT NONE
  INTEGER     , INTENT(IN   ) :: it
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: ku1, uo
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: ku2, ku3, ku4
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kv1, vo
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kv2, kv3, kv4
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kt1, to
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kt2, kt3, kt4
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kq1, qo
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kq2, kq3, kq4
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kz1, zo
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kz2, kz3, kz4
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kw1, wo
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kw2, kw3, kw4

  REAL(KIND=r8) :: dt6
  REAL(KIND=r8) :: dt2
  REAL(KIND=r8) :: dt
  INTEGER      :: i,j,k
  dt=DeltaT;  dt6 = DeltaT/6.0; dt2=DeltaT*0.5

  IF(it==1)THEN
     DO k=1,Kdim
        DO j=1,Jdim
           DO i=1,  Idim
              uo(i,j,k) =  u_ref(i,j,k)
              vo(i,j,k) =  v_ref(i,j,k)
              wo(i,j,k) =  w_ref(i,j,k)
              to(i,j,k) =  t_ref(i,j,k)
              qo(i,j,k) =  q_ref(i,j,k)
              zo(i,j,k) =  z_ref(i,j,k)

              U_C(i,j,k)=  u_ref(i,j,k)
              V_C(i,j,k)=  v_ref(i,j,k)
              W_C(i,j,k)=  w_ref(i,j,k)
              T_C(i,j,k)=  t_ref(i,j,k)
              Q_C(i,j,k)=  q_ref(i,j,k)
              Z_C(i,j,k)=  z_ref(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     DO k=1,Kdim
        DO j=1,Jdim
           DO i=1,  Idim
              uo(i,j,k)=  U_C(i,j,k)
              vo(i,j,k)=  V_C(i,j,k) 
              wo(i,j,k)=  W_C(i,j,k) 
              to(i,j,k)=  T_C(i,j,k)
              qo(i,j,k)=  Q_C(i,j,k)
              zo(i,j,k)=  Z_C(i,j,k)
           END DO
        END DO
     END DO
  END IF
 !       PRINT*,'h eta',MAXVAL(h),MINVAL(h),'u eta',MAXVAL(u),MINVAL(u),'v eta',MAXVAL(v),MINVAL(v)
 
 call Pettersen_Sutcliffe_Model(uo        , vo            , to         , qo        , zo        , wo        ,t_ref)
   
 ! final step and time marching / new values for RK4
 DO k=1,Kdim-1
    DO j=2,Jdim-1
       DO i=2,  Idim-1
          U_N(i,j,k) =  u_ref(i,j,k)!uo(i,j,k) + (ku1(i,j,k) + 2.0*ku2(i,j,k) + 2.0*ku3(i,j,k) + ku4(i,j,k))*dt6
          V_N(i,j,k) =  v_ref(i,j,k)!vo(i,j,k) + (kv1(i,j,k) + 2.0*kv2(i,j,k) + 2.0*kv3(i,j,k) + kv4(i,j,k))*dt6
          W_N(i,j,k) =  w_ref(i,j,k)!vo(i,j,k) + (kv1(i,j,k) + 2.0*kv2(i,j,k) + 2.0*kv3(i,j,k) + kv4(i,j,k))*dt6
          T_N(i,j,k) =  t_ref(i,j,k)!to(i,j,k) + (kt1(i,j,k) + 2.0*kt2(i,j,k) + 2.0*kt3(i,j,k) + kt4(i,j,k))*dt6
          Q_N(i,j,k) =  q_ref(i,j,k)!qo(i,j,k) + (kq1(i,j,k) + 2.0*kq2(i,j,k) + 2.0*kq3(i,j,k) + kq4(i,j,k))*dt6
          Z_N(i,j,k) =  z_ref(i,j,k)!zo(i,j,k) + (kz1(i,j,k) + 2.0*kz2(i,j,k) + 2.0*kz3(i,j,k) + kz4(i,j,k))*dt6
       END DO
    END DO
 END DO

 ! updating the data
 DO k=1,Kdim-1
    DO j=2,Jdim-1
       DO i=2,  Idim-1
          U_C(i,j,k) = U_N(i,j,k) 
          V_C(i,j,k) = V_N(i,j,k)
          W_C(i,j,k) = W_N(i,j,k)
          T_C(i,j,k) = T_N(i,j,k)
          Q_C(i,j,k) = Q_N(i,j,k)
          Z_C(i,j,k) = Z_N(i,j,k)
       END DO
    END DO
 END DO
 PRINT*,it,MAXVAL(ug_surf(2:Idim,2:Jdim)), MINVAL(ug_surf(2:Idim,2:Jdim)),& !TESTAR
           MAXVAL(vg_surf(2:Idim,2:Jdim)), MINVAL(vg_surf(2:Idim,2:Jdim))
 END SUBROUTINE RunDynamics
 
 
 SUBROUTINE  Pettersen_Sutcliffe_Model(u_in, v_in, t_in,q_in,z_in,w_in,t_new)
  IMPLICIT NONE
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:kdim), intent(in) ::  u_in,  v_in,  t_in , q_in, z_in,w_in,t_new
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:kdim) ::  u,  v,  t , q, z, w,zeta_out, div_out
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:kdim) ::  ug_out, vg_out
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:kdim) ::  LapOperOmega_qg
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:kdim) ::  SigmaLapOperOmega_qg
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:kdim) ::  LapOperHeat_qg
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,1:kdim) ::  termA,Q_heating
  REAL(KIND=r8), dimension(1:Idim, 1:Jdim,3,1:kdim) ::  GradOperZeta
  REAL(KIND=r8) :: udux, udvx, udqx, vdTx   
  REAL(KIND=r8) :: vduy, vdvy, vdqy, udTy
  REAL(KIND=r8) :: wduz, wdvz, wdqz, wdTz
  REAL(KIND=r8) :: fcov,minabs,pLND
  REAL(KIND=r8) :: hh  
  REAL(KIND=r8) :: dPdx, dPdy
  REAL(KIND=r8) :: lnPs_xf
  REAL(KIND=r8) :: lnPs_xb
  REAL(KIND=r8) :: Ps_xf
  REAL(KIND=r8) :: Ps_xc
  REAL(KIND=r8) :: Tv, kTv, kTv0,f0
  REAL(KIND=r8) :: vis2dudx, vis2dvdx, vis2dqdx
  REAL(KIND=r8) :: vis2dudy, vis2dvdy, vis2dqdy 
  REAL(KIND=r8) :: factor4,factor2
  REAL(KIND=r8) :: term1,uadvc,vadvc,wadvc
  REAL(KIND=r8) :: TermUNewton, TermVNewton, TermTNewton, TermQNewton
  REAL(KIND=r8) :: factor_dlambda
  REAL(KIND=r8) :: factor_dtheta 
  real(kind=r8) :: dx(Idim)
  real(kind=r8) :: dy(Jdim)
  real(kind=r8) :: div_profile(1:kdim)
  real(kind=r8) :: sigma_ref
  logical :: have_theta
  INTEGER :: i,j,k,iLND
  INTEGER :: xb,xc,xf
  INTEGER :: yb,yc,yf

  real(kind=r8) :: AdvcOperThichGeo(1:Idim, 1:Jdim)
  real(kind=r8) :: SigmaLapOperOmegaMed(1:Idim, 1:Jdim)
  real(kind=r8) :: LapOperHeatMed(1:Idim, 1:Jdim)
  ! QG-omega solver arrays
  real(kind=r8) :: omega_qg(1:Idim, 1:Jdim,1:kdim)
  real(kind=r8) :: rhs_qg  (1:Idim, 1:Jdim,1:kdim)
  REAL(KIND=r8) :: w_med(1:Idim, 1:Jdim)
  ! solver control params (tune these!)
      real(r8), parameter :: omegapar = 1.7_r8   ! SOR relaxation
  real(kind=r8) :: nu            ! horizontal diffusion coefficient (m^2 / s) effective scale (tunable)
  integer :: maxiter
  real(kind=r8) :: tol_sor
  integer :: smooth_freq, smooth_passes
  ug_out=0.0_r8; vg_out=0.0_r8;div_out=0.0_r8;zeta_out=0.0_r8;div_profile=0.0_r8;GradOperZeta=0.0_r8
  termA=0.0_r8; omega_qg=0.0_r8;termB_QG=0.0_r8;termC_QG=0.0_r8;termD_QG=0.0_r8;ug_surf=0.0_r8;vg_surf=0.0_r8
  ug_iLND=0.0_r8;vg_iLND=0.0_r8
  u=u_in;v=v_in;t=t_in;q=q_in;z=z_in;w=w_in
  
  !O QUE ISSO FAZ?
  
  u(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(u_ref(1:Idim   ,1:Jdim,Kdim) + u(1:Idim,1:Jdim,Kdim-1))

  v(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(v_ref(1:Idim   ,1:Jdim,Kdim) + v(1:Idim,1:Jdim,Kdim-1))

  t(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(t_ref(1:Idim   ,1:Jdim,Kdim))!+ t(1:Idim,1:Jdim,Kdim-1))

  q(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(q_ref(1:Idim   ,1:Jdim,Kdim))!+ q(1:Idim,1:Jdim,Kdim-1))
 
  z(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(z_ref(1:Idim   ,1:Jdim,Kdim))!+ q(1:Idim,1:Jdim,Kdim-1))

  w(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(w_ref(1:Idim   ,1:Jdim,Kdim))!+ q(1:Idim,1:Jdim,Kdim-1))

 !IMPLEMENTAR AQUI AMORTECIMENTO COM 4 PONTOS

  DO k=1,Kdim

     u(1   ,1:Jdim,k)=0.25_r8*(u_ref(1   ,1:Jdim,k) + u(2     ,1:Jdim,k) + u(3     ,1:Jdim,k) + u_ref(2     ,1:Jdim,k) )
     u(Idim,1:Jdim,k)=0.25_r8*(u_ref(Idim,1:Jdim,k) + u(Idim-1,1:Jdim,k) + u(Idim-2,1:Jdim,k) + u_ref(Idim-1,1:Jdim,k) )

     u(1:Idim,1   ,k)=0.25_r8*(u_ref(1:Idim,1   ,k) + u(1:Idim,2     ,k) + u(1:Idim,3     ,k) + u_ref(1:Idim,2     ,k) )
     u(1:Idim,Jdim,k)=0.25_r8*(u_ref(1:Idim,Jdim,k) + u(1:Idim,Jdim-1,k) + u(1:Idim,Jdim-2,k) + u_ref(1:Idim,Jdim-1,k) )

     v(1   ,1:Jdim,k)=0.25_r8*(v_ref(1   ,1:Jdim,k) + v(2     ,1:Jdim,k) + v(3     ,1:Jdim,k) + v_ref(2     ,1:Jdim,k) )
     v(Idim,1:Jdim,k)=0.25_r8*(v_ref(Idim,1:Jdim,k) + v(Idim-1,1:Jdim,k) + v(Idim-2,1:Jdim,k) + v_ref(Idim-1,1:Jdim,k) )

     v(1:Idim,1   ,k)=0.25_r8*(v_ref(1:Idim,1   ,k) + v(1:Idim,2     ,k) + v(1:Idim,3     ,k) + v_ref(1:Idim,2     ,k) )
     v(1:Idim,Jdim,k)=0.25_r8*(v_ref(1:Idim,Jdim,k) + v(1:Idim,Jdim-1,k) + v(1:Idim,Jdim-2,k) + v_ref(1:Idim,Jdim-1,k) )

     t(1   ,1:Jdim,k)=0.25_r8*(t_ref(1   ,1:Jdim,k) + t(2     ,1:Jdim,k) + t(3     ,1:Jdim,k) + t_ref(2     ,1:Jdim,k))
     t(Idim,1:Jdim,k)=0.25_r8*(t_ref(Idim,1:Jdim,k) + t(Idim-1,1:Jdim,k) + t(Idim-2,1:Jdim,k) + t_ref(Idim-1,1:Jdim,k))

     t(1:Idim,1   ,k)=0.25_r8*(t_ref(1:Idim,1   ,k) + t(1:Idim,2     ,k) + t(1:Idim,3     ,k) + t_ref(1:Idim,2     ,k))
     t(1:Idim,Jdim,k)=0.25_r8*(t_ref(1:Idim,Jdim,k) + t(1:Idim,Jdim-1,k) + t(1:Idim,Jdim-2,k) + t_ref(1:Idim,Jdim-1,k))

     q(1   ,1:Jdim,k)=0.25_r8*(q_ref(1   ,1:Jdim,k) + q(2     ,1:Jdim,k) + q(3     ,1:Jdim,k) + q_ref(2     ,1:Jdim,k))
     q(Idim,1:Jdim,k)=0.25_r8*(q_ref(Idim,1:Jdim,k) + q(Idim-1,1:Jdim,k) + q(Idim-2,1:Jdim,k) + q_ref(Idim-1,1:Jdim,k))

     q(1:Idim,1   ,k)=0.25_r8*(q_ref(1:Idim,1   ,k) + q(1:Idim,2     ,k) + q(1:Idim,3     ,k) + q_ref(1:Idim,2     ,k))
     q(1:Idim,Jdim,k)=0.25_r8*(q_ref(1:Idim,Jdim,k) + q(1:Idim,Jdim-1,k) + q(1:Idim,Jdim-2,k) + q_ref(1:Idim,Jdim-1,k))

     z(1   ,1:Jdim,k)=0.25_r8*(z_ref(1   ,1:Jdim,k) + z(2     ,1:Jdim,k) + z(3     ,1:Jdim,k) + z_ref(2     ,1:Jdim,k))
     z(Idim,1:Jdim,k)=0.25_r8*(z_ref(Idim,1:Jdim,k) + z(Idim-1,1:Jdim,k) + z(Idim-2,1:Jdim,k) + z_ref(Idim-1,1:Jdim,k))

     z(1:Idim,1   ,k)=0.25_r8*(z_ref(1:Idim,1   ,k) + z(1:Idim,2     ,k) + z(1:Idim,3     ,k) + z_ref(1:Idim,2     ,k))
     z(1:Idim,Jdim,k)=0.25_r8*(z_ref(1:Idim,Jdim,k) + z(1:Idim,Jdim-1,k) + z(1:Idim,Jdim-2,k) + z_ref(1:Idim,Jdim-1,k))

     w(1   ,1:Jdim,k)=0.25_r8*(w_ref(1   ,1:Jdim,k) + w(2     ,1:Jdim,k) + w(3     ,1:Jdim,k) + w_ref(2     ,1:Jdim,k))
     w(Idim,1:Jdim,k)=0.25_r8*(w_ref(Idim,1:Jdim,k) + w(Idim-1,1:Jdim,k) + w(Idim-2,1:Jdim,k) + w_ref(Idim-1,1:Jdim,k))

     w(1:Idim,1   ,k)=0.25_r8*(w_ref(1:Idim,1   ,k) + w(1:Idim,2     ,k) + w(1:Idim,3     ,k) + w_ref(1:Idim,2     ,k))
     w(1:Idim,Jdim,k)=0.25_r8*(w_ref(1:Idim,Jdim,k) + w(1:Idim,Jdim-1,k) + w(1:Idim,Jdim-2,k) + w_ref(1:Idim,Jdim-1,k))

  END DO

  DO j=1,Jdim
     DO i=1,Idim
        factor_dlambda=(r_earth)*(cos(CoordLat(i,j))**2)
        factor_dtheta =(r_earth)*(cos(CoordLat(i,j))**1)
        dx(i) = (factor_dlambda*DeltaLamda(i,j))
        dy(j) = (factor_dtheta *DeltaTheta(i,j))
     END DO
  END DO 
        
  !PARTE ESPACIAL
  DO k=1,Kdim-1 
      ! compute geostrophic wind, vorticity, divergence
      CALL geostrophic_from_phi(z(1:Idim, 1:Jdim,k),dx,dy,ug_out(1:Idim, 1:Jdim,k), vg_out(1:Idim, 1:Jdim,k))
      CALL compute_vort_div(ug_out(1:Idim, 1:Jdim,k), vg_out(1:Idim, 1:Jdim,k), dx, dy, zeta_out(1:Idim, 1:Jdim,k), div_out(1:Idim, 1:Jdim,k))
  END DO
  ! find LND: search for level where divergence crosses zero between consecutive levels
  iLND = -1
  DO k=1,Kdim-1 
      ! check sign change somewhere in lat-lon (we consider global search per column later)
      ! Here we compute mean over domain to find global LND index (approx)
      div_profile(k) = sum(abs(div_out(1:Idim, 1:Jdim,k)))/real(Jdim*Idim,r8)
  END DO
  ! find level where abs(div_profile) minimum
  minabs = 1.0e30_r8
  DO k = kMinLev, Kdim-1
      IF (abs(div_profile(k)) < minabs .and.Plevs(k) >= 50000.0) then
        minabs = abs(div_profile(k))
        iLND = k
      ENDIF
  END DO
  !div_profile(1:Kdim)=abs(Plevs(1:Kdim)-50000.0_r8)
  !iLND=minloc(div_profile,dim=1)
  IF (iLND < 1) iLND = 1
  pLND = Plevs(iLND)/100.0_r8
  
  ! choose a representative f0 (Coriolis at domain mean latitude)
  f0 = sum(FcorPar(1,1:Jdim))/Jdim
  if (abs(f0) < 1e-6_r8) f0 = 1.0e-4_r8

  !PARTE ESPACIAL
  DO k=kMinLev,Kdim-1 
      ! compute gradient , vorticity, zeta
      CALL gradient(zeta_out(1:Idim, 1:Jdim,k),&
                    dx                        ,&
                    dy                        ,&
                    GradOperZeta(1:Idim, 1:Jdim,1:3,k))

      CALL compute_termA_from_zeta(GradOperZeta(1:Idim, 1:Jdim,1:3,k), &
                                   ug_out(1:Idim, 1:Jdim,k), &
                                   vg_out(1:Idim, 1:Jdim,k), &
                                   f0                      , &
                                   termA(1:Idim, 1:Jdim,k))

  END DO
  ! Choose a representative sigma for denominator f0^2 / sigma
  do yc=1,Jdim
     do xc=1,Idim
        thetas(xc,yc) = 0.5_r8*(t(xc,yc,kMinLev)*((Plevs(1)/Plevs(kMinLev))**(287.0_r8/1004.6_r8)) + &
                                t(xc,yc,kMinLev+1)*((Plevs(1)/Plevs(kMinLev+1))**(287.0_r8/1004.6_r8))  )
     end do
  end do

  ! if theta available compute approximate sigma; else use a typical midtropospheric value:
  have_theta=.true.  
  if (have_theta) then
      ! crude estimate: sigma_ref ~ (p/theta) dtheta/dp averaged over mid-layer
      sigma_ref = estimate_sigma_from_theta(t, Plevs, Kdim)
  else
      sigma_ref = 2.5e-6_r8   ! typical value (tunable)
  end if

  if (sigma_ref <= 0.0_r8 ) sigma_ref = 2.5e-6_r8  
  DO yc=1,Jdim
     DO xc=1,Idim
        ug_surf (xc,yc) =  vertical_mean_p(Plevs(kMinLev:kMinLev+1), ug_out(xc, yc,kMinLev:kMinLev+1))
        vg_surf (xc,yc) =  vertical_mean_p(Plevs(kMinLev:kMinLev+1), vg_out(xc, yc,kMinLev:kMinLev+1))
     END DO
  END DO

  do yc=1,Jdim
      do xc=1,Idim

        termA_QG(xc,yc) = - (1.0_r8/sigma_ref*(Plevs(iLND)-Plevs(kMinLev)))*termA(xc,yc,iLND)
        termQ_QG(xc,yc) = 0.0_r8 
        termS_QG(xc,yc) = 0.0_r8 


        ug_iLND (xc,yc) = ug_out(xc,yc,iLND)
        vg_iLND (xc,yc) = vg_out(xc,yc,iLND)
        !
        ! total_forcing is the sum of A + Q + S at this stage (S=0 until omega computed)
        total_forcing(xc,yc) = termA_QG(xc,yc) + termQ_QG(xc,yc) + termS_QG(xc,yc)
        ! Build RHS for QG-omega solver:
        ! The QG-omega eqn has form: Lap_h(omega) + (f0^2/sigma) d2omega/dp2 = F  (RHS forcing)
        ! Here we use F = - total_forcing * f0^2  (this is schematic ? replace by full physical form)
        ! For demonstration use a safe scaling so magnitudes are reasonable:
        DO k=kMinLev,Kdim-1 
           rhs_qg(xc,yc,k) = - total_forcing(xc,yc) * (f0*f0)
        END DO
     END DO
  END DO
  nu            = 1.0e4_r8          ! ~1e4 m^2/s as example (tune). Acts like viscosity smoothing small scales.
  maxiter       = 8000
  tol_sor       = 1.0e-8_r8
  !
  !     compute inverse Laplacian to get quase geostrophyc omega
  !
  DO k=kMinLev,Kdim-1 
     call multigrid_solver(rhs_qg(1:Idim,1:Jdim,k), omega_qg(1:Idim,1:Jdim,k), Idim,Jdim, dx, dy, tol_sor)
  END DO 

  call Advection(z(1:Idim, 1:Jdim,iLND)-z(1:Idim, 1:Jdim,kMinLev),&
                 ug_surf (1:Idim, 1:Jdim)                  ,&
                 vg_surf (1:Idim, 1:Jdim)                  ,&
                 dx                                        ,&
                 dy                                        ,&
                 AdvcOperThichGeo                           )

  CALL laplacian(AdvcOperThichGeo(1:Idim, 1:Jdim), dx, dy, termB_QG(1:Idim, 1:Jdim))

  DO j=1,Jdim
     DO i=1,Idim
        termA_QG(i, j) = + termA(i, j,iLND)
        termB_QG(i, j) = - termB_QG(i, j)

        w_med(i, j) = vertical_mean_p(Plevs(kMinLev:iLND), omega_qg(i,j,kMinLev:iLND))
     END DO
  END DO


  DO k=kMinLev,Kdim-1 
      !CALL laplacian(w(1:Idim, 1:Jdim,k), dx, dy, LapOperOmega_qg(1:Idim, 1:Jdim,k))
      !CALL laplacian(omega_qg(1:Idim, 1:Jdim,k), dx, dy, LapOperOmega_qg(1:Idim, 1:Jdim,k))
      CALL laplacian(w_med(1:Idim, 1:Jdim), dx, dy, LapOperOmega_qg(1:Idim, 1:Jdim,k))
     SigmaLapOperOmega_qg(1:Idim, 1:Jdim,k) = sigma_ref*LapOperOmega_qg(1:Idim, 1:Jdim,k)
  END DO 

  DO j=1,Jdim
     DO i=1,Idim
        SigmaLapOperOmegaMed(i, j) = vertical_mean_p(Plevs(kMinLev:iLND), SigmaLapOperOmega_qg(i,j,kMinLev:iLND))


        termC_QG(i, j) = -SigmaLapOperOmegaMed(i, j)*(Plevs(iLND)-Plevs(kMinLev))

     END DO
  END DO
  

  DO k=kMinLev,Kdim
      DO j=1,Jdim
         DO i=1,Idim
            Q_heating(i,j,k) =(t_new(i,j,k) - t(i,j,k))/DeltaT
         END DO
      END DO
  END DO
  DO j=1,Jdim
     DO i=1,Idim
        total_forcing(i, j) = vertical_mean_p(Plevs(kMinLev:iLND), Q_heating(i,j,kMinLev:iLND))
     END DO
  END DO

  DO k=kMinLev,Kdim-1 
     CALL laplacian(Q_heating(1:Idim, 1:Jdim,k), dx, dy, LapOperHeat_qg(1:Idim, 1:Jdim,k))
     LapOperHeat_qg(1:Idim, 1:Jdim,k) = (Rd/Cp)*LapOperHeat_qg(1:Idim, 1:Jdim,k)
  END DO 

  DO j=1,Jdim
     DO i=1,Idim
        LapOperHeatMed(i, j) = vertical_mean_p(Plevs(kMinLev:iLND), LapOperHeat_qg(i,j,kMinLev:iLND))

        termD_QG(i, j) = -LapOperHeatMed(i, j)*(log(Plevs(iLND))-log(Plevs(kMinLev)))

     END DO
  END DO


 END SUBROUTINE  Pettersen_Sutcliffe_Model



  !==============================================================
  !  Calcula média ponderada simples por pressão
  !==============================================================
  pure function vertical_mean_p(p, x) result(xmean)
    implicit none
    real(kind=r8), intent(in) :: p(:), x(:)
    real(kind=r8)              :: xmean
    real(kind=r8)              :: dp(size(p,dim=1))
    real(kind=r8)              ::  num, den
    integer :: k,kdim
    xmean=0.0_r8
    kdim=size(p,dim=1)
    do k = 1, kdim-1
       dp(k) = abs(p(k+1) - p(k))
    end do
    dp(kdim) = dp(kdim-1)

    num = 0.0_r8
    den = 0.0_r8

    do k = 1, kdim
       num = num + x(k) * dp(k)
       den = den + dp(k)
    end do

    if (den /= 0.0) then
       xmean = num / den
    else
       xmean = 0.0
    end if   
  end function vertical_mean_p

!-------------------------------------------------------
  subroutine apply_3x3_smoother2d(field, passes)
    implicit none
    real(kind=r8), intent(inout) :: field(:,:)
    integer, intent(in) :: passes
    real(kind=r8) :: temp(size(field,1),size(field,2))
    real(kind=r8) :: wcenter, wneigh
    integer  :: pass, k, i, j,yb,yc,yf,xb,xc,xf, Jdim, Idim

    Jdim = size(field,2); Idim = size(field,1)
    wcenter = 0.5_r8
    wneigh  = 0.125_r8  ! 4 neighbors * 0.125 + center 0.5 = 1.0
    do pass = 1, passes
       temp = field
          do j = 2, Jdim-1
             CALL index(j,Jdim,yb,yc,yf)
             DO i=2,Idim-1
                CALL index(i,Idim,xb,xc,xf)
                field(xc,yc) = wcenter*temp(xc,yc) + wneigh*( temp(xf,yc) + temp(xb,yc) + temp(xc,yf) + temp(xc,yb) )
             end do
          end do
    end do
    ! Condições simples nas bordas (Neumann: derivada normal = 0)
    do xc = 1, Idim
       field(xc,   1)  = field(xc,     2)
       field(xc,Jdim)  = field(xc,Jdim-1)
    end do
    do yc = 1, Jdim
       field(1   ,yc)  = field(2     ,yc)
       field(Idim,yc)  = field(Idim-1,yc)
    end do
  end subroutine apply_3x3_smoother2d
  !-------------------------------------------------------
!-------------------------------------------------------
  subroutine apply_3x3_smoother(field, passes)
    implicit none
    real(kind=r8), intent(inout) :: field(:,:,:)
    integer, intent(in) :: passes
    real(kind=r8) :: temp(size(field,1),size(field,2),size(field,3))
    real(kind=r8) :: wcenter, wneigh
    integer  :: pass, k, i, j,yb,yc,yf,xb,xc,xf, kdim, Jdim, Idim

    kdim = size(field,3); Jdim = size(field,2); Idim = size(field,1)
    wcenter = 0.5_r8
    wneigh  = 0.125_r8  ! 4 neighbors * 0.125 + center 0.5 = 1.0
    do pass = 1, passes
       temp = field
       do k = 1, kdim
          do j = 2, Jdim-1
             CALL index(j,Jdim,yb,yc,yf)
             DO i=2,Idim-1
                CALL index(i,Idim,xb,xc,xf)
                field(xc,yc,k) = wcenter*temp(xc,yc,k) + wneigh*( temp(xf,yc,k) + temp(xb,yc,k) + temp(xc,yf,k) + temp(xc,yb,k) )
             end do
          end do
       end do
    end do
    do k = 1, kdim
    ! Condições simples nas bordas (Neumann: derivada normal = 0)
       do xc = 1, Idim
         field(xc,   1,k)  = field(xc,     2,k)
         field(xc,Jdim,k)  = field(xc,Jdim-1,k)
       end do
       do yc = 1, Jdim
         field(1   ,yc,k)  = field(2     ,yc,k)
         field(Idim,yc,k)  = field(Idim-1,yc,k)
       end do
    end do
  end subroutine apply_3x3_smoother
  !-------------------------------------------------------
  ! Estimate sigma from theta (crude layer-mean)
  pure function estimate_sigma_from_theta(temp, levels, Kdim) result(sigmean)
    implicit none
    integer      , intent(in) :: Kdim
    real(kind=r8), intent(in) :: temp(:,:,:), levels(:)
    real(kind=r8) :: sigmean
    real(kind=r8) :: accum
    real(kind=r8) :: theta(size(temp,dim=1),size(temp,dim=2),size(temp,dim=3))
    real(kind=r8) :: locsum, dp_level
    integer :: k,i,j,xc,yc,Jdim,Idim, count
    
    accum = 0.0_r8; count = 0
    Idim=size(temp,dim=1);Jdim=size(temp,dim=2)
    do k=1,Kdim
       do yc=1,Jdim
          do xc=1,Idim
             theta(xc,yc,k) = temp(xc,yc,k)*((levels(1)/levels(k))**(287.0_r8/1004.6_r8))
          end do
       end do
    end do


    do k=2,Kdim-1
       ! approximate p/theta * dtheta/dp at level k, average over domain
       locsum   = 0.0_r8
       dp_level = (levels(k+1)-levels(k-1))
       do yc=1,Jdim
          do xc=1,Idim
             locsum = locsum + (levels(k)/ max(theta(xc,yc,k),1e-6_r8)) * &
                               ( (theta(xc,yc,k+1)-theta(xc,yc,k-1)) / dp_level )
          end do
       end do
       accum = accum + locsum / real(Idim*Jdim,r8)
       count = count + 1
    end do
    if (count > 0) then
       sigmean = accum / real(count,r8)
    else 
       sigmean = 2.5e-6_r8
    end if
  end function estimate_sigma_from_theta
  
  
!-------------------------------------------------------
  subroutine compute_termA_from_zeta(GradOperZeta, ug, vg, f0, termA)
    implicit none
    real(kind=r8), intent(in) :: GradOperZeta(:,:,:), ug(:,:), vg(:,:)
    real(kind=r8), intent(in) :: f0
    real(kind=r8), intent(out) :: termA(size(GradOperZeta,1),size(GradOperZeta,2))
    integer :: Idim,Jdim, i, j, yb,yc,yf,xb,xc,xf,versor_i,versor_j
    Idim = size(GradOperZeta,1); Jdim = size(GradOperZeta,2)
    versor_i=1;versor_j=2
    do yc=1,Jdim
       do xc=1,Idim
          termA(xc,yc) = -  f0 * ( -ug(xc,yc)*GradOperZeta(xc,yc,versor_i) - vg(xc,yc)*GradOperZeta(xc,yc,versor_j) )
       end do
    end do
  end subroutine compute_termA_from_zeta
  
 !-------------------------------------------------------
subroutine laplacian(phi, dx, dy, lap)
  implicit none
  real(kind=r8), intent(in)  :: dx(:), dy(:)
  real(kind=r8), intent(in)  :: phi(:, :)
  real(kind=r8), intent(out) :: lap(:, :)
  integer :: Idim,Jdim, i, j, yb,yc,yf,xb,xc,xf
  Idim = size(phi,1); Jdim = size(phi,2)

  ! Zera saida
  lap = 0.0_r8

  ! Diferencas centrais no interior
  DO j=2,Jdim-1
     CALL index(j,Jdim,yb,yc,yf)
      DO i=2,Idim-1
         CALL index(i,Idim,xb,xc,xf)
         lap(xc,yc) = (phi(xf, yc) - 2.0_r8*phi(xc,yc) + phi(xb,yc)) / (dx(xc)*dx(xc)) + &
                      (phi(xc ,yf) - 2.0_r8*phi(xc,yc) + phi(xc,yb)) / (dy(yc)*dy(yc))
     end do
  end do

  ! Condições simples nas bordas (Neumann: derivada normal = 0)
  do xc = 1, Idim
      lap(xc,   1)  = lap(xc,     2)
      lap(xc,Jdim)  = lap(xc,Jdim-1)
  end do
  do yc = 1, Jdim
      lap(1   ,yc)  = lap(2     ,yc)
      lap(Idim,yc)  = lap(Idim-1,yc)
  end do
  call apply_3x3_smoother2d(lap, passes=1)

end subroutine laplacian
!-------------------------------------------------------
 
subroutine gradient(phi, dx, dy, GradOper)
  implicit none
  real(kind=r8), intent(in ) :: dx(:), dy(:)
  real(kind=r8), intent(in ) :: phi(:, :)
  real(kind=r8), intent(out) :: GradOper(size(phi,1), size(phi,2),3)
  real(kind=r8) :: dphidx(size(phi,1), size(phi,2))
  real(kind=r8) :: dphidy(size(phi,1), size(phi,2))
  integer :: Idim,Jdim, i, j, yb,yc,yf,xb,xc,xf
  Idim = size(phi,1); Jdim = size(phi,2)

  

  ! Zera bordas (condicoes simples)

  dphidx = 0.0_r8
  dphidy = 0.0_r8

  ! Diferencas centrais no interior

  DO j=2,Jdim-1
     CALL index(j,Jdim,yb,yc,yf)
      DO i=2,Idim-1
         CALL index(i,Idim,xb,xc,xf)
         dphidx(xc,yc) = (phi(xf,yc) - phi(xb,yc)) / (2.0*dx(xc))
         dphidy(xc,yc) = (phi(xc,yf) - phi(xc,yb)) / (2.0*dy(yc))
      END DO
  END DO

  ! Condicoes nas bordas (diferenças para frente/para tras)
  do yc = 1, Jdim
      dphidx(1   ,yc)  = (phi(2   ,yc) - phi(1     ,yc)) / dx(1   )
      dphidx(Idim,yc)  = (phi(Idim,yc) - phi(Idim-1,yc)) / dx(Idim)
  end do

  do xc = 1, Idim
      dphidy(xc,   1)  = (phi(xc,   2) - phi(xc,     1)) / dy(1   )
      dphidy(xc,Jdim)  = (phi(xc,Jdim) - phi(xc,Jdim-1)) / dy(Jdim)
  end do
  
  call apply_3x3_smoother2d(dphidy, passes=1)

  call apply_3x3_smoother2d(dphidy, passes=1)


  do yc = 1, Jdim
     do xc = 1, Idim
        GradOper(xc,yc,1) = dphidx(xc,yc)
        GradOper(xc,yc,2) = dphidy(xc,yc)
        GradOper(xc,yc,3) = 0.0_r8
     end do
  end do

end subroutine gradient

!-------------------------------------------------------
 
subroutine Advection(phi,ug2d, vg2d,  dx, dy, AdvcOper)
  implicit none
  real(kind=r8), intent(in ) :: dx(:), dy(:)
  real(kind=r8), intent(in ) :: phi(:, :),ug2d(:, :), vg2d(:, :)
  real(kind=r8), intent(out) :: AdvcOper(size(phi,1), size(phi,2))
  real(kind=r8) :: VectorGradOper(size(phi,1), size(phi,2),3)
  real(kind=r8) :: Vector2AuxOper(size(phi,1), size(phi,2),3)
  real(kind=r8) :: div   (size(phi,1), size(phi,2))
  real(kind=r8) :: Scaler(size(phi,1), size(phi,2))
  
  integer :: Idim,Jdim, i, j, yb,yc,yf,xb,xc,xf
  Idim = size(phi,1); Jdim = size(phi,2)
  Scaler=phi
  call gradient(Scaler, dx, dy, VectorGradOper)
  DO  yc= 1, Jdim
    do xc = 1, Idim
       Vector2AuxOper(xc,yc,1)=ug2d(xc,yc)
       Vector2AuxOper(xc,yc,2)=vg2d(xc,yc)
    end do
  END DO  

  CALL divergence(Vector2AuxOper, VectorGradOper, dx, dy, div)
 
  ! Zera bordas (condicoes simples)

  AdvcOper = 0.0_r8

  ! Diferencas centrais no interior

  DO j=2,Jdim-1
     CALL index(j,Jdim,yb,yc,yf)
      DO i=2,Idim-1
         CALL index(i,Idim,xb,xc,xf)
         AdvcOper(xc,yc) = div(xc,yc)
      END DO
  END DO

  ! Condicoes nas bordas (diferenças para frente/para tras)

  do yc = 1, Jdim
      AdvcOper(1   ,yc)  = AdvcOper(2     ,yc)
      AdvcOper(Idim,yc)  = AdvcOper(Idim-1,yc)
  end do

  do xc = 1, Idim
      AdvcOper(xc,   1)  = AdvcOper(xc,     2)
      AdvcOper(xc,Jdim)  = AdvcOper(xc,Jdim-1)
  end do

end subroutine Advection
 
 
 subroutine divergence(vector1, vector2, dx, dy, div)
  !==========================================================
  ! Calcula divergencia 2D: naBla·V = du/dx + dv/dy
  !----------------------------------------------------------
  ! Entradas:
  !   u(nx,ny), v(nx,ny) : componentes do vetor
  !   dx, dy             : espacamentos em metros
  ! Saída:
  !   div(nx,ny)         : divergencia 
  !==========================================================
  implicit none
  real(kind=r8)   , intent(in) :: dx(:), dy(:)
  real(kind=r8)   , intent(in) :: vector1(:,:,:), vector2(:,:,:)
  real(kind=r8)   , intent(out) :: div(size(dx,dim=1), size(dy,dim=1))
  integer :: Idim,Jdim, i, j, yb,yc,yf,xb,xc,xf

  Idim=size(dx,dim=1);Jdim=size(dy,dim=1)

  div = 0.0_r8

  ! Diferenças centrais no interior

  DO j=2,Jdim-1
     CALL index(j,Jdim,yb,yc,yf)
      DO i=2,Idim-1
         CALL index(i,Idim,xb,xc,xf)
         div(xc,yc)      =  vector1(xc,yc,1) * vector2(xc,yc,1) + &
                            vector1(xc,yc,2) * vector2(xc,yc,2)

      END DO
  END DO

  ! Condicoes nas bordas (diferenças para frente/para tras)
  do yc = 1, Jdim
      div(1   ,yc)  = div(2     ,yc)
      div(Idim,yc)  = div(Idim-1,yc)
  end do

  do xc = 1, Idim
      div(xc,   1)  = div(xc,     2)
      div(xc,Jdim)  = div(xc,Jdim-1)
  end do
end subroutine divergence

 !-------------------------------------------------------
  subroutine compute_vort_div(ug2d, vg2d, dx, dy, zeta_out, div_out)
    implicit none
    real(kind=r8), intent(in ) :: ug2d(:,:), vg2d(:,:), dx(:), dy(:)
    real(kind=r8), intent(out) :: zeta_out(size(ug2d,1),size(ug2d,2))
    real(kind=r8), intent(out) :: div_out (size(ug2d,1),size(ug2d,2))
    integer :: i,j,yb,yc,yf,xb,xc,xf
    Idim = size(ug2d,1); Jdim = size(ug2d,2)
      DO j=2,Jdim-1
         CALL index(j,Jdim,yb,yc,yf)
         DO i=2,Idim-1
            CALL index(i,Idim,xb,xc,xf)
            zeta_out(xc,yc) = ((vg2d(xf,yc)-vg2d(xb,yc))/(2.0_r8*dx(xc))) - ((ug2d(xc,yf)-ug2d(xc,yb))/(2.0_r8*dy(yc))) 
            div_out (xc,yc) = ((ug2d(xf,yc)-ug2d(xb,yc))/(2.0_r8*dx(xc))) + ((vg2d(xc,yf)-vg2d(xc,yb))/(2.0_r8*dy(yc)))
         END DO
      END DO
      zeta_out(1   ,1:Jdim)=0.8_r8*( zeta_out(2     ,1:Jdim) )
      zeta_out(Idim,1:Jdim)=0.8_r8*( zeta_out(Idim-1,1:Jdim) )

      zeta_out(1:Idim,1   )=0.8_r8*( zeta_out(1:Idim,2     ) )
      zeta_out(1:Idim,Jdim)=0.8_r8*( zeta_out(1:Idim,Jdim-1) )

      div_out(1   ,1:Jdim)=0.8_r8*( div_out(2     ,1:Jdim) )
      div_out(Idim,1:Jdim)=0.8_r8*( div_out(Idim-1,1:Jdim) )

      div_out(1:Idim,1   )=0.8_r8*( div_out(1:Idim,2     ) )
      div_out(1:Idim,Jdim)=0.8_r8*( div_out(1:Idim,Jdim-1) )

  end subroutine compute_vort_div
 !-------------------------------------------------------
  subroutine geostrophic_from_phi(phi2d, dx, dy, ug_out, vg_out)
    implicit none
    real(kind=r8), intent(in ) :: phi2d(:,:), dx(:), dy(:)
    real(kind=r8), intent(out) :: ug_out(size(dx,dim=1),size(dy,dim=1))
    real(kind=r8), intent(out) :: vg_out(size(dx,dim=1),size(dy,dim=1))
    integer :: i,j,Idim,Jdim,xb,xc,xf,yb,yc,yf
    Idim=size(dx,dim=1)
    Jdim=size(dy,dim=1)
    !DO k=1,Kdim-1 
       DO j=2,Jdim-1
          CALL index(j,Jdim,yb,yc,yf)
          DO i=2,Idim-1
             CALL index(i,Idim,xb,xc,xf)
             !
             !          1       D FHI
             ! ug = - ----- * ---------
             !          f0       D y
             !
             ug_out(xc,yc) = - (1.0_r8 / FcorPar(xc,yc)) * ((phi2d(xc,yf)-phi2d(xc,yb))/(2.0*dy(yc))) !d_dy(phi2d, j, i, dy)
             !
             !          1       D FHI
             ! vg = + ----- * ---------
             !          f0       D x
             !
             vg_out(xc,yc) =   (1.0_r8 / FcorPar(xc,yc)) * ((phi2d(xf,yc)-phi2d(xb,yc))/(2.0*dx(xc))) ! d_dx(phi2d, j, i, dx)
          END DO
       END DO
      ug_out(1   ,1:Jdim)=( ug_out(2	 ,1:Jdim) )
      ug_out(Idim,1:Jdim)=( ug_out(Idim-1,1:Jdim) )

      ug_out(1:Idim,1   )=( ug_out(1:Idim,2	) )
      ug_out(1:Idim,Jdim)=( ug_out(1:Idim,Jdim-1) )

      vg_out(1   ,1:Jdim)=( vg_out(2	 ,1:Jdim) )
      vg_out(Idim,1:Jdim)=( vg_out(Idim-1,1:Jdim) )

      vg_out(1:Idim,1   )=( vg_out(1:Idim,2	) )
      vg_out(1:Idim,Jdim)=( vg_out(1:Idim,Jdim-1) )
 
    !END DO

  end subroutine geostrophic_from_phi
  
   SUBROUTINE index(i,Idim,xb,xc,xf)
      IMPLICIT NONE
      INTEGER, INTENT(IN   ) :: i
      INTEGER, INTENT(IN   ) :: Idim
      INTEGER, INTENT(OUT  ) :: xb,xc,xf
      IF(i==1) THEN
        xb=Idim
        xc=i
        xf=i+1
      ELSE IF(i==Idim)THEN
        xb=Idim-1
        xc=Idim
        xf=1
      ELSE
        xb=i-1
        xc=i
        xf=i+1
      END IF
   END SUBROUTINE index
   
   SUBROUTINE Finalize_Class_Module_Dynamics() !FINALIZE
   	IMPLICIT NONE
   	
   END SUBROUTINE Finalize_Class_Module_Dynamics
   
END MODULE Class_Module_Dynamics
