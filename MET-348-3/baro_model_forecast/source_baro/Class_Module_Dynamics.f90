!  $Author: pkubota 						$
!  $Date: 2008/09/23 17:51:54 					$
!  $Revision: 1.9 						$
!  $Revisions are currently made by the class's students.	$
!  $Update Date: 01/11/2023 10:19 AM				$
!  
!  Implementações: 
!  	1) Colocar INIT e FINALIZE - OK 
! 	2) Colocar as equações dos campos no "Class_Module_Dynamics" - OK
! 	7) Amortecimento das condições de contorno - OK (funcionando parcialmente)

MODULE Class_Module_Dynamics
 USE Constants, Only: r8,r4,i4,pi,Deg2Rad,Rd,Cp,kappa,r_earth,omega,CTv,nfprt, Eps
 
 USE Class_Module_LapaceOperator, Only:sor_solve,&
                               multigrid_solver
 
 USE Class_Module_Fields, Only: u_ref, v_ref,w_ref, t_ref,q_ref,z_ref,p_ref,s_ref,k_start,k_finish,&
                                U_N,U_C,  V_N,V_C,  T_N,T_C,  Q_N,Q_C,Z_C,Z_N,S_C,S_N,&
                                Plevs,CoordLat,CoordLon,FcorPar ,BETAPar,DeltaLamda,DeltaTheta
  IMPLICIT NONE
  PRIVATE       

  INTEGER ::  Idim
  INTEGER ::  Jdim
  INTEGER ::  Kdim
  REAL(KIND=r8) :: DeltaT
  REAL(KIND=8), PUBLIC   ,parameter        :: vis    =  2.0e2_r8! 1.5e-10 !1.5e-5        !  viscosity
  REAL(KIND=8), PUBLIC   ,parameter        :: taul   = 3600_r8
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
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: ku1, uo
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: ku2, ku3, ku4
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kv1, vo
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kv2, kv3, kv4
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kt1, to
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kt2, kt3, kt4
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kq1, qo
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kq2, kq3, kq4
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kz1, zo
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kz2, kz3, kz4
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: ks1, so
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: ks2, ks3, ks4
  real(kind=r8) :: dx(Idim)
  real(kind=r8) :: dy(Jdim)
  real(kind=r8) :: factor_dlambda
  real(kind=r8) :: factor_dtheta 
  REAL(KIND=8) :: dt6
  REAL(KIND=8) :: dt2
  REAL(KIND=8) :: dt
  INTEGER      :: i,j,k
  dt=DeltaT;  dt6 = DeltaT/6.0; dt2=DeltaT*0.5

  IF(it==1)THEN
     DO k=1,Kdim
        DO j=1,Jdim
           DO i=1,  Idim
              uo(i,j,k) =  u_ref(i,j,k)
              vo(i,j,k) =  v_ref(i,j,k)
              to(i,j,k) =  t_ref(i,j,k)
              qo(i,j,k) =  q_ref(i,j,k)
              zo(i,j,k) =  z_ref(i,j,k)
              so(i,j,k) =  s_ref(i,j,k)

              U_C(i,j,k)=  u_ref(i,j,k)
              V_C(i,j,k)=  v_ref(i,j,k)
              T_C(i,j,k)=  t_ref(i,j,k)
              Q_C(i,j,k)=  q_ref(i,j,k)
              Z_C(i,j,k)=  z_ref(i,j,k)
              S_C(i,j,k)=  s_ref(i,j,k)

           END DO
        END DO
     END DO
  ELSE
     DO k=1,Kdim
        DO j=1,Jdim
           DO i=1,  Idim
              uo(i,j,k)=  U_C(i,j,k)
              vo(i,j,k)=  V_C(i,j,k) 
              to(i,j,k)=  T_C(i,j,k)
              qo(i,j,k)=  Q_C(i,j,k)
              zo(i,j,k)=  Z_C(i,j,k)
              so(i,j,k)=  S_C(i,j,k)
           END DO
        END DO
     END DO
  END IF
 !       PRINT*,'h eta',MAXVAL(h),MINVAL(h),'u eta',MAXVAL(u),MINVAL(u),'v eta',MAXVAL(v),MINVAL(v)
 
 call Solve_Forward_Beta_plane_Grid_A(ku1, kv1, kt1,kq1,ks1, uo        , vo		, to	     , qo	 , so	     )
 call Solve_Forward_Beta_plane_Grid_A(ku2, kv2, kt2,kq2,ks2, uo+ku1*dt2, vo+kv1*dt2	, to+kt1*dt2 , qo+kq1*dt2, so+ks1*dt2)
 call Solve_Forward_Beta_plane_Grid_A(ku3, kv3, kt3,kq3,ks3, uo+ku2*dt2, vo+kv2*dt2	, to+kt2*dt2 , qo+kq2*dt2, so+ks2*dt2)
 call Solve_Forward_Beta_plane_Grid_A(ku4, kv4, kt4,kq4,ks4, uo+ku3*dt , vo+kv3*dt	, to+kt3*dt  , qo+kq3*dt , so+ks3*dt )
   
 ! final step and time marching / new values for RK4
 DO k=k_start,k_finish
    DO j=2,Jdim-1
       DO i=2,  Idim-1
          S_N(i,j,k) = so(i,j,k) + (ks1(i,j,k) + 2.0*ks2(i,j,k) + 2.0*ks3(i,j,k) + ks4(i,j,k))*dt6
       END DO
    END DO
 END DO

 ! updating the data
 DO k=k_start,k_finish
    DO j=2,Jdim-1
       DO i=2,  Idim-1
          S_C(i,j,k) = S_N(i,j,k)
       END DO
    END DO
 END DO

 ! updating the data
 DO k=k_start,k_finish
    DO j=3,Jdim-3
       DO i=3,  Idim-3
          factor_dlambda=(r_earth)*(cos(CoordLat(i,j))**2)
          factor_dtheta =(r_earth)*(cos(CoordLat(i,j))**1)
          dx(i) = (factor_dlambda*DeltaLamda(i,j))
          dy(j) = (factor_dtheta *DeltaTheta(i,j))
          U_N(i,j,k) = -(S_N(i  ,j+1,k) - S_N(i  ,j-1,k))/(2*dy(j))
          V_N(i,j,k) =  (S_N(i+1,j  ,k) - S_N(i-1,j  ,k))/(2*dx(i))
          U_C(i,j,k) = U_N(i,j,k)
          V_C(i,j,k) = V_N(i,j,k)
       END DO
    END DO
 END DO
 PRINT*,it,MAXVAL(S_C  (2:Idim,2:Jdim,k_start:k_finish)), MINVAL(S_C  (2:Idim,2:Jdim,k_start:k_finish)),& !TESTAR
           MAXVAL(w_ref(2:Idim,2:Jdim,k_start:k_finish)), MINVAL(w_ref(2:Idim,2:Jdim,k_start:k_finish))
 END SUBROUTINE RunDynamics
 
 
 SUBROUTINE  Solve_Forward_Beta_plane_Grid_A(TermEqMomU, TermEqMomV,TermEqConT, TermEqConQ,TermEqConS, u_in, v_in, t_in,q_in,s_in)
  IMPLICIT NONE
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:kdim), intent(in) ::  u_in,  v_in,  t_in , q_in, s_in
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:kdim), intent(out) :: TermEqMomU,TermEqMomV, TermEqConT, TermEqConQ,TermEqConS
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:kdim) ::  u,  v,  t , q, z, s, vtx,dsdt
  REAL(KIND=8) :: udux, udvx, udqx, vdTx   
  REAL(KIND=8) :: vduy, vdvy, vdqy, udTy
  REAL(KIND=8) :: wduz, wdvz, wdqz, wdTz
  REAL(KIND=8) :: fcov
  REAL(KIND=8) :: hh  
  REAL(KIND=8) :: dPdx, dPdy
  REAL(KIND=8) :: lnPs_xf
  REAL(KIND=8) :: lnPs_xb
  REAL(KIND=8) :: Ps_xf
  REAL(KIND=8) :: Ps_xc
  REAL(KIND=8) :: Tv, kTv, kTv0
  REAL(KIND=8) :: vis2dudx, vis2dvdx, vis2dqdx
  REAL(KIND=8) :: vis2dudy, vis2dvdy, vis2dqdy 
  REAL(KIND=8) :: factor4,factor2
  REAL(KIND=8) :: term1,uadvc,vadvc,wadvc
  REAL(KIND=8) :: TermUNewton, TermVNewton, TermTNewton, TermQNewton
  INTEGER :: i,j,k
  INTEGER :: xb,xc,xf
  INTEGER :: yb,yc,yf
  TermEqMomU=0.0_r8;TermEqMomV=0.0_r8; TermEqConT=0.0_r8; TermEqConQ=0.0_r8;dsdt=0.0_r8
  u=u_in;v=v_in;t=t_in;q=q_in;s=s_in
  
  !O QUE ISSO FAZ?
  
  u(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(u_ref(1:Idim   ,1:Jdim,Kdim) + u(1:Idim,1:Jdim,Kdim-1))

  v(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(v_ref(1:Idim   ,1:Jdim,Kdim) + v(1:Idim,1:Jdim,Kdim-1))

  t(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(t_ref(1:Idim   ,1:Jdim,Kdim))!+ t(1:Idim,1:Jdim,Kdim-1))

  q(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(q_ref(1:Idim   ,1:Jdim,Kdim))!+ q(1:Idim,1:Jdim,Kdim-1))

  z(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(z_ref(1:Idim   ,1:Jdim,Kdim))!+ q(1:Idim,1:Jdim,Kdim-1))

  s(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(s_ref(1:Idim   ,1:Jdim,Kdim))!+ q(1:Idim,1:Jdim,Kdim-1))

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

     s(1   ,1:Jdim,k)=0.25_r8*(s_ref(1   ,1:Jdim,k) + s(2     ,1:Jdim,k) + s(3     ,1:Jdim,k) + s_ref(2     ,1:Jdim,k))
     s(Idim,1:Jdim,k)=0.25_r8*(s_ref(Idim,1:Jdim,k) + s(Idim-1,1:Jdim,k) + s(Idim-2,1:Jdim,k) + s_ref(Idim-1,1:Jdim,k))

     s(1:Idim,1   ,k)=0.25_r8*(s_ref(1:Idim,1   ,k) + s(1:Idim,2     ,k) + s(1:Idim,3     ,k) + s_ref(1:Idim,2     ,k))
     s(1:Idim,Jdim,k)=0.25_r8*(s_ref(1:Idim,Jdim,k) + s(1:Idim,Jdim-1,k) + s(1:Idim,Jdim-2,k) + s_ref(1:Idim,Jdim-1,k))

  END DO
  !PARTE ESPACIAL
  !
  !     a) compute vorticity from streamfunction
  !
  !print*,"mkvort"
  call mkvort(s,vtx,Idim,Jdim,Kdim) 
  !
  !     b) compute streamfunction tendencies
  !
  !print*,"tendency"

  call tendency(s,vtx,dsdt,Idim,Jdim,Kdim)

  !PARTE ESPACIAL
  DO k=k_start,k_finish
      DO j=2,Jdim-1
         CALL index(j,Jdim,yb,yc,yf)
         DO i=2,Idim-1
            CALL index(i,Idim,xb,xc,xf)
            factor4=(r_earth**2)*(cos(CoordLat(xc,yc))**4)
            factor2=(r_earth**2)*(cos(CoordLat(xc,yc))**2)

            TermQNewton = (s(xc,yc,k) - s_ref(xc,yc,k))/taul !Amortecimento para as forçantes
            !      --       --            !                 
            ! d(s)            |d(d(s))    |
            ! -----  - Neta * |--------   | = 0
            ! dt              |dxdx       |
            !                  -         -
            vis2dqdx= -vis*(1/factor4)*((s(xf,yc,k) - 2.0*s(xc,yc,k) + s(xb,yc,k))/(DeltaLamda(xc,yc)*DeltaLamda(xc,yc)))
            !
            !                  -       -
            ! d(s)             |d(d(s)) |
            ! -----   - Neta * |--------| = 0
            ! dt               |dydy    |
            !                  -       -

            vis2dqdy= -vis*(1.0/factor2)*((s(xc,yf,k) - 2.0*s(xc,yc,k) + s(xc,yb,k))/(DeltaTheta(xc,yc)*DeltaTheta(xc,yc)))

            TermEqMomU(xc,yc,k) = dsdt(xc,yc,k) - (vis2dqdx + vis2dqdy + TermQNewton)
            TermEqMomV(xc,yc,k) = dsdt(xc,yc,k) - (vis2dqdx + vis2dqdy + TermQNewton)
            TermEqConT(xc,yc,k) = dsdt(xc,yc,k) - (vis2dqdx + vis2dqdy + TermQNewton)
            TermEqConQ(xc,yc,k) = dsdt(xc,yc,k) - (vis2dqdx + vis2dqdy + TermQNewton)
            TermEqConS(xc,yc,k) = dsdt(xc,yc,k) - (vis2dqdx + vis2dqdy + TermQNewton)
         END DO
      END DO
  END DO


 END SUBROUTINE  Solve_Forward_Beta_plane_Grid_A

!
!----------------------------------------------------------------------
!
 subroutine tendency(s,vtx,dsdt,Idim,Jdim,Kdim)
      !real, INTENT(inout) :: dsdt(0:NX+1,0:NY+1) ! streamfct tendnecy (ds/dt) [m**2/s**2] 
      integer, intent(in) :: Idim,Jdim,Kdim
      real(kind=r8)    :: s  (Idim,Jdim,kdim)     ! input: stream function
      real(kind=r8)    :: vtx(Idim,Jdim,kdim)    ! output: vorticity
      real(kind=r8)    :: dsdt (Idim,Jdim,kdim)  ! ! outpout: inverse Laplacian of input
      integer :: i,j,k ,yb,yc,yf,xb,xc,xf
!
!     subroutine tendency computes the streamfunction tendency 
!
      real(kind=r8) :: zdt(Idim,Jdim,kdim)  ! vorticity tendency
      real(kind=r8) :: zj (Idim,Jdim,kdim)  ! jacobian (advection of rel vort.)
      real(kind=r8) :: zv (Idim,Jdim,kdim)  ! y-velocity v (ds/dx)

      real(kind=r8) :: st (Idim,Jdim,kdim)  ! y-velocity v (ds/dx)
      real(kind=r8) :: vo (Idim,Jdim,kdim)  ! y-velocity v (ds/dx)
      real(kind=r8) :: dif(Idim,Jdim,kdim) 
      real(kind=r8) :: dx(Idim)
      real(kind=r8) :: dy(Jdim)
      real(kind=r8) :: factor_dlambda
      real(kind=r8) :: factor_dtheta 
      real(kind=r8), parameter :: tol = 1.0e-8_r8
      integer, parameter :: maxiter_sor = 20000
      real(r8), parameter :: omega = 1.7_r8   ! SOR relaxation

      zdt=0.0_r8;      zj =0.0_r8;      zv =0.0_r8;      st =0.0_r8;      vo =0.0_r8

      DO k=k_start,k_finish
         DO j=1,Jdim
            DO i=1,Idim
               factor_dlambda=(r_earth)*(cos(CoordLat(i,j))**2)
               factor_dtheta =(r_earth)*(cos(CoordLat(i,j))**1)
               dx(i) = (factor_dlambda*DeltaLamda(i,j))
               dy(j) = (factor_dtheta *DeltaTheta(i,j))
               st(i,j,k)=s(i,j,k)
               vo(i,j,k)=vtx(i,j,k)
               dsdt=vtx(i,j,k)
            END DO
         END DO 
      END DO 
      !
      !     compute jacobian J(s,vo) (adv. of rel. vort.)
      !
      call jacobi(st,vo,zj,Idim,Jdim,kdim)
      !
      !     compute  (ds/dx)
      !
      call mkdfdx(st,zv,Idim,Jdim,kdim)
      !
      !     compute ni * (d2s/d2x)
      !
      !call mkd2fd2x(vo,dif,Idim,Jdim,kdim)
      !
      !     add advection of rel. and planetray vorticity to get vort. tendency
      !
      !print*,"advection"

      DO k=k_start,k_finish
         DO j=2,Jdim-1
            CALL index(j,Jdim,yb,yc,yf)
            DO i=2,Idim-1
               CALL index(i,Idim,xb,xc,xf)
               zdt(xc,yc,k)= - zj(xc,yc,k) - BETAPar(xc,yc)*zv(xc,yc,k) 
            END DO
         END DO 
      END DO 
      !
      !     compute inverse Laplacian to get streamfunction tendency
      !
      DO k=k_start,k_finish
!      call sor_solve(zdt, dsdt, Idim,Jdim, dx, dy, tol, maxiter_sor, omega)   

         call multigrid_solver(zdt(1:Idim,1:Jdim,k), dsdt(1:Idim,1:Jdim,k), Idim,Jdim, dx, dy, tol)
      END DO 

      return
 end subroutine tendency

 subroutine sor(pdf,pf,Idim,Jdim,kdim)
  implicit none
!
!     subroutine sor compute the inverse Laplacian from a given field
!     by using the Successive OverRelaxation method (SOR)
!
      real,parameter :: wsor = 1.50 ! overrelaxaton factor
      real,parameter :: eps = 1.E-4 ! minimum reduction of error to stop iter.
!
      integer :: Idim                ! x dimension,Idim,Jdim,kdim
      integer :: Jdim                ! y dimension
      integer :: kdim                !  z dimension
      real(kind=r8)    :: pdf(Idim,Jdim,kdim)   ! input: field 
      real(kind=r8)    :: pf (Idim,Jdim,kdim)    ! outpout: inverse Laplacian of input
!
      integer :: niter             ! maximum number of iterations
      real(kind=r8) :: dx                  ! x grid point distance
      real(kind=r8) :: dy                  ! y grid point distance
      real(kind=r8) :: zn(Idim,Jdim)    ! work array: new result
      real(kind=r8) :: zo(Idim,Jdim)    ! work array: old result
      real(kind=r8) :: zfac                 ! help
      real(kind=r8) :: zzres                ! local error
      real(kind=r8) :: zresf                ! initial error
      real(kind=r8) :: zres                 ! actual error
      real(kind=r8) :: factor_dlambda
      real(kind=r8) :: factor_dtheta 
      real(kind=r8) :: factor4       
      real(kind=r8) :: factor2       
      integer :: i,j,jiter         ! loop indizes
      integer :: k,yb ,yc,yf,xb,xc,xf               ! loop indizes

!
!     set maximum number of iteration 
!     (note: sor should need about sqrt(kx*ky)/3*log10(1/eps) iterations to 
!      decrease the initial error by factor eps... but it takes more(?))
!
      niter=(Idim*Jdim*kdim)*NINT(log10(1./eps)+0.5)/3 
!
!
!     copy first guess to work (here first guess is 0.)
!
      zo(:,:)=0. ! pf(:,:)
      zn(:,:)=0. ! pf(:,:)
!
!     compute initial error
!
!
!     df/dx=(f(i+1)-f(i-1))/2dx
!
      DO k=k_start,k_finish
         !
         !     copy first guess to work (here first guess is 0.)
         !
         zo(:,:)=0. ! pf(:,:)
         zn(:,:)=0. ! pf(:,:)
         zresf=0.
         DO j=2,Jdim-1
            CALL index(j,Jdim,yb,yc,yf)
            DO i=2,Idim-1
               CALL index(i,Idim,xb,xc,xf)
               factor_dlambda=(r_earth)*(cos(CoordLat(xc,yc))**2)
               factor_dtheta =(r_earth)*(cos(CoordLat(xc,yc))**1)
               factor4       =(r_earth**2)*(cos(CoordLat(xc,yc))**4)
               factor2       =(r_earth**2)*(cos(CoordLat(xc,yc))**2)
               dx = (factor_dlambda*DeltaLamda(xc,yc))
               dy = (factor_dtheta *DeltaTheta(xc,yc))
               zfac=2./(dx*dx)+2./(dy*dy)
               !
               zzres = - zfac*zn(xc,yc) + (zn(xb,yc) + zn(xf,yc))/(dx*dx) &
     &                 +                  (zn(xc,yb) + zn(xc,yf))/(dy*dy) - pdf(xc,yc,k)
             
               zresf = zresf + abs(zzres)/real(Idim*Jdim)

           enddo
         enddo
!        do the iteration
!
         do jiter=1,niter
            zo(:,:)=zn(:,:)
            zres=0.
            DO j=2,Jdim-1
                CALL index(j,Jdim,yb,yc,yf)
                DO i=2,Idim-1
                   CALL index(i,Idim,xb,xc,xf)
                   factor_dlambda=(r_earth)*(cos(CoordLat(xc,yc))**2)
                   factor_dtheta =(r_earth)*(cos(CoordLat(xc,yc))**1)
                   factor4       =(r_earth**2)*(cos(CoordLat(xc,yc))**4)
                   factor2       =(r_earth**2)*(cos(CoordLat(xc,yc))**2)
                   dx = (factor_dlambda*DeltaLamda(xc,yc))
                   dy = (factor_dtheta *DeltaTheta(xc,yc))
                   zfac=2./(dx*dx)+2./(dy*dy)

                   zzres=-zfac*zo(xc,yc)+(zn(xb,yc)+zo(xf,yc))/(dx*dx)            &
                                        +(zn(xc,yb)+zo(xc,yf))/(dy*dy)-pdf(xc,yc,k)

                   zn(xc,yc)=zo(xc,yc)+wsor*zzres/zfac
                   zres=zres + abs(zzres)/real(Idim*Jdim)
                enddo
            enddo
            !      ps(0,:)    = ps(kx,:)
            !      ps(kx+1,:) = ps(1 ,:)
            !      ps(:,0)    = 0.0
            !      ps(:,ky+1) = 0.0
            zn(1,     :) = 0.5*(zn(2    ,:) + zn(3      ,:))
            zn(Idim,:) = 0.5*(zn(Idim-1 ,:) + zn(Idim-2 ,:))
            zn(:,     1) = 0.5*(zn(: ,2   ) + zn(: ,3)     )
            zn(:,Jdim) = 0.5*(zn(: ,Jdim-1) + zn(: ,Jdim-2))

            if(zres <= eps*zresf) exit
         enddo
!
         if(jiter >= niter) then 
            print*,'MAXIMUM iterations needed',niter
            print*,'error, initial error, eps*ie= ',zres,zresf,eps*zresf
         endif
!
!     copy result to output
!
         DO j=2,Jdim-1
             CALL index(j,Jdim,yb,yc,yf)
             DO i=2,Idim-1
                CALL index(i,Idim,xb,xc,xf)
                 pf(xc,yc,k)=zn(xc,yc)
             enddo
         enddo

!
      enddo

!      zresf=0.
!      do j=1,ky
!       do i=1,kx
!        zzres=-zfac*zn(i,j)+(zn(i-1,j)+zn(i+1,j))/(pdx*pdx)            &
!     &       +(zn(i,j-1)+zn(i,j+1))/(pdy*pdy)-pdf(i,j)
!        zresf=zresf+abs(zzres)/real(kx*ky)
!       enddo
!      enddo
!
!     do the iteration
!
!     do jiter=1,niter
!      zo(:,:)=zn(:,:)
!      zres=0.
!      do j=1,ky
!       do i=1,kx
!        zzres=-zfac*zo(i,j)+(zn(i-1,j)+zo(i+1,j))/(pdx*pdx)            &
!    &        +(zn(i,j-1)+zo(i,j+1))/(pdy*pdy)-pdf(i,j)
!        zn(i,j)=zo(i,j)+wsor*zzres/zfac
!        zres=zres+abs(zzres)/real(kx*ky)
!       enddo
!      enddo
!      call boundary(zn,kx,ky)
!      if(zres <= eps*zresf) exit
!     enddo
!
!     if(jiter >= niter) then 
!      print*,'MAXIMUM iterations needed',niter
!      print*,'error, initial error, eps*ie= ',zres,zresf,eps*zresf
!     endif
!
!     copy result to output
!
!     pf(:,:)=zn(:,:)
!
      return
      end subroutine sor

!-----------------------------------------------------------------------
!
      subroutine mkd2fd2x(q,pdfdx,Idim,Jdim,kdim)
      implicit none
!
!     subroutine mkdfdx computes x derivation from a field 
!     using central differences
!
      integer :: Idim                       ! x-dimension
      integer :: Jdim                       ! y-dimension
      integer :: kdim                       ! z-dimension
      real(kind=r8)    :: q   (Idim,Jdim,kdim)      ! input: field
      real(kind=r8)    :: pdfdx(Idim,Jdim,kdim)  ! output: dfield/dx
      real(kind=r8)    :: dfdx  ! output: dfield/dx
      real(kind=r8)    :: factor4
      real(kind=r8)    :: factor2
      real(kind=r8)    :: factor_dlambda
      real(kind=r8)    :: factor_dtheta
      real(kind=r8)    :: dx
      real(kind=r8)    :: dy
      real(kind=r8)    :: zx  ,vis2dqdx      ,vis2dqdy             ! 2.*dx
      integer :: i,j,k,yb ,yc,yf,xb,xc,xf               ! loop indizes
!
!     df/dx=(f(i+1)-f(i-1))/2dx
!
      DO k=k_start,k_finish
         DO j=2,Jdim-1
            CALL index(j,Jdim,yb,yc,yf)
            DO i=2,Idim-1
               CALL index(i,Idim,xb,xc,xf)
               factor_dlambda=(r_earth)*(cos(CoordLat(xc,yc))**2)
               factor_dtheta =(r_earth)*(cos(CoordLat(xc,yc))**1)
               factor4=(r_earth**2)*(cos(CoordLat(xc,yc))**4)
               factor2=(r_earth**2)*(cos(CoordLat(xc,yc))**2)
               dx = (factor_dlambda*DeltaLamda(xc,yc))
               dy = (factor_dtheta *DeltaTheta(xc,yc))
               zx=2.0*dx
               !      --       --            !                 
               ! d(q)            |d(d(q))    |
               ! -----  - Neta * |--------   | = 0
               ! dt              |dxdx       |
               !                  -         -
               vis2dqdx= -vis*(1/factor4)*((q(xf,yc,k) - 2.0*q(xc,yc,k) + q(xb,yc,k))/(DeltaLamda(xc,yc)*DeltaLamda(xc,yc)))
               !
               !                  -       -
               ! d(q)             |d(d(q)) |
               ! -----   - Neta * |--------| = 0
               ! dt               |dydy    |
               !                  -       -

               vis2dqdy= -vis*(1.0/factor2)*((q(xc,yf,k) - 2.0*q(xc,yc,k) + q(xc,yb,k))/(DeltaTheta(xc,yc)*DeltaTheta(xc,yc)))



               !
               pdfdx(i,j,k)=vis2dqdx + vis2dqdy
           enddo
         enddo
      enddo
!
!
      return
      end subroutine mkd2fd2x
!

!
!-----------------------------------------------------------------------
!
      subroutine mkdfdx(pf,pdfdx,Idim,Jdim,kdim)
      implicit none
!
!     subroutine mkdfdx computes x derivation from a field 
!     using central differences
!
      integer :: Idim                       ! x-dimension
      integer :: Jdim                       ! y-dimension
      integer :: kdim                       ! z-dimension
      real(kind=r8)    :: pf   (Idim,Jdim,kdim)      ! input: field
      real(kind=r8)    :: pdfdx(Idim,Jdim,kdim)  ! output: dfield/dx
      real(kind=r8)    :: dfdx  ! output: dfield/dx
      real(kind=r8)    :: factor4
      real(kind=r8)    :: factor2
      real(kind=r8)    :: factor_dlambda
      real(kind=r8)    :: factor_dtheta
      real(kind=r8)    :: dx
      real(kind=r8)    :: dy
      real(kind=r8)    :: zx                     ! 2.*dx
      integer :: i,j,k,yb ,yc,yf,xb,xc,xf               ! loop indizes
!
!     df/dx=(f(i+1)-f(i-1))/2dx
!
      DO k=k_start,k_finish
         DO j=2,Jdim-1
            CALL index(j,Jdim,yb,yc,yf)
            DO i=2,Idim-1
               CALL index(i,Idim,xb,xc,xf)
               factor_dlambda=(r_earth)*(cos(CoordLat(xc,yc))**2)
               factor_dtheta =(r_earth)*(cos(CoordLat(xc,yc))**1)
               factor4=(r_earth**2)*(cos(CoordLat(xc,yc))**4)
               factor2=(r_earth**2)*(cos(CoordLat(xc,yc))**2)
               dx = (factor_dlambda*DeltaLamda(xc,yc))
               dy = (factor_dtheta *DeltaTheta(xc,yc))
               zx=2.0*dx
               !
               pdfdx(i,j,k)=(pf(xf,yc,k)-pf(xb,yc,k))/zx
           enddo
         enddo
      enddo
!
!
      return
      end subroutine mkdfdx
!
!
!----------------------------------------------------------------------
!
      subroutine mkvort(s,vtx,Idim,Jdim,kdim)
      implicit none
!
!     subroutine mkvort computes the vorticity from a given streamfunction
!
      integer :: Idim                 ! x-dimension
      integer :: Jdim                 ! y-dimension
      integer :: kdim                 ! z-dimension
      real(kind=r8)    :: s (Idim,Jdim,kdim)     ! input: stream function
      real(kind=r8)    :: vtx(Idim,Jdim,kdim)    ! output: vorticity
!      integer :: i,j                ! loop indizes
!
!     compute Laplacian
!
      call laplace(s,vtx,Idim,Jdim,kdim)
!
!     make boundary conditions
!
!      call boundary(pvo,kx,ky)
!
      return
      end subroutine mkvort

!
!-----------------------------------------------------------------------
!
      subroutine laplace(s,vtx,Idim,Jdim,kdim)
      implicit none
!
!     subroutine laplace computes the laplacian from a field 
!
      integer :: Idim                 ! x-dimension
      integer :: Jdim                 ! y-dimension
      integer :: kdim                 ! z-dimension
      real(kind=r8) :: d2qd2x,d2qd2y,factor4,factor2
      real(kind=r8) :: s (Idim,Jdim,kdim)     ! input: field
      real(kind=r8) :: vtx(Idim,Jdim,kdim)    ! output: Laplacian
      real(kind=r8) :: q(Idim,Jdim,kdim)    ! output: Laplacian
      integer :: i,j,k,yb ,yc,yf,xb,xc,xf               ! loop indizes

      DO k=k_start,k_finish
         DO j=1,Jdim
            DO i=1,Idim
               q(i,j,k)=s(i,j,k)
               vtx(i,j,k)=0.0
            END DO
         END DO 
      END DO 

!      DO k=1,kdim
!         q(0     ,0:Jdim+1,k)=0.5*(s(1   ,1:Jdim,k) + s_ref(1     ,1:Jdim,k))
!         q(Idim+1,0:Jdim+1,k)=0.5*(s(Idim,1:Jdim,k) + s_ref(Idim  ,1:Jdim,k))

!         q(0:Idim+1,0     ,k)=0.5*(s(1:Idim,1   ,k) + s_ref(1:Idim,1     ,k))
!         q(0:Idim+1,Jdim+1,k)=0.5*(s(1:Idim,Jdim,k) + s_ref(1:Idim,Jdim  ,k))

!      END DO 
!
      DO k=k_start,k_finish
         DO j=2,Jdim-1
            CALL index(j,Jdim,yb,yc,yf)
            DO i=2,Idim-1
               CALL index(i,Idim,xb,xc,xf)
               factor4=(r_earth**2)*(cos(CoordLat(xc,yc))**4)
               factor2=(r_earth**2)*(cos(CoordLat(xc,yc))**2)

               !                    --        --
               ! dd(q)             |d(d(q))    |
               ! -----       =     |--------   | 
               ! dxdx              |dxdx       |
               !                   -          -
               d2qd2x= (1/factor4)*((q(xf,yc,k) - 2.0*q(xc,yc,k) + q(xb,yc,k))/(DeltaLamda(xc,yc)*DeltaLamda(xc,yc)))
               !
               !                    --        --
               ! dd(q)             |d(d(q))    |
               ! -----       =     |--------   | 
               ! dydy              |dydy       |
               !                   -          -

               d2qd2y= (1.0/factor2)*((q(xc,yf,k) - 2.0*q(xc,yc,k) + q(xc,yb,k))/(DeltaTheta(xc,yc)*DeltaTheta(xc,yc)))

               vtx(i,j,k)=d2qd2x  + d2qd2y
           enddo
         enddo
      enddo
!
      return
      end subroutine laplace

      subroutine jacobi(st,vo,pj,Idim,Jdim,kdim)
      implicit none
!
!     subroutine jacobi computes the jacobi operator according to 
!     Arakawa (1966)
!
!     jacobi(1,2)=(j1+j2+j3)/3. after Arakawa 1966
!,
      integer :: Idim              ! x-dimension
      integer :: Jdim              ! y-dimension
      integer :: kdim              ! y-dimension
      real(kind=r8) :: st(Idim,Jdim,kdim)  ! input: field 1  real :: p1(Idim,Jdim,kdim)  ! input: field 1
      real(kind=r8) :: vo(Idim,Jdim,kdim)  ! input: field 2  real :: p2(Idim,Jdim,kdim)  ! input: field 2
      real(kind=r8) :: pj(Idim,Jdim,kdim)  ! output: jacobi(1,2)
      real(kind=r8) :: pdx                ! x grid distance
      real(kind=r8) :: pdy                ! y grid distance
      real(kind=r8) :: zx                 ! 2.*dx
      real(kind=r8) :: zy                 ! 2.*dy
      real(kind=r8) :: DxXDy
      real(kind=r8) :: zj1(Idim,Jdim)         ! 1st jacobi
      real(kind=r8) :: zj2(Idim,Jdim)         ! 2nd jacobi
      real(kind=r8) :: zj3(Idim,Jdim)         ! 3rd jacobi 
      real(kind=r8) :: var_a(Idim,Jdim,kdim)  ! input: field 1
      real(kind=r8) :: var_b(Idim,Jdim,kdim)  ! input: field 2
      real(kind=r8) :: factor_dlambda
      real(kind=r8) :: factor_dtheta 
      integer :: i,j,k,yb ,yc,yf,xb,xc,xf               ! loop indizes
zj1=0.0_r8;var_a=0.0_r8
zj2=0.0_r8;var_b=0.0_r8
zj3=0.0_r8
      var_a=st
      var_b=vo
!
!    _____________ ________________________
!   |            |            |            |
!   | a(i-1,j+1) | a(i  ,j+1) | a(i+1,j+1) |
!   | b(i-1,j+1) | b(i  ,j+1) | b(i+1,j+1) |
!   |____________|____________|____________|
!   |            |            |            |
!   | a(i-1,j  ) | a(i  ,j  ) | a(i+1,j  ) |
!   | b(i-1,j  ) | b(i  ,j  ) | b(i+1,j  ) |
!   |____________|_____ ______|_____ ______|
!   |            |            |            |
!   | a(i-1,j-1) | a(i  ,j-1) | a(i+1,j-1) |
!   | b(i-1,j-1) | b(i  ,j-1) | b(i+1,j-1) |
!   |____________|____________|____________|
!
      DO k=k_start,k_finish
         DO j=2,Jdim-1
            CALL index(j,Jdim,yb,yc,yf)
            DO i=2,Idim-1
               CALL index(i,Idim,xb,xc,xf)
               factor_dlambda=(r_earth)*(cos(CoordLat(xc,yc))**2)
               factor_dtheta =(r_earth)*(cos(CoordLat(xc,yc))**1)
               zx = (factor_dlambda*DeltaLamda(xc,yc))
               zy = (factor_dtheta *DeltaTheta(xc,yc))
               DxXDy = 4.0*(zx * zy)
               !
               !                     --                                                      
               !             1       |                                                       
               ! J1     = ---------  | ( a1(i+1,j) -a3(i-1,j) ) *( b2(i,j+1) -b4(i,j-1) )  - 
               !           4*Dx*Dy   |                                                       
               !                     -                                                       
               !
               !                                                                          --
               !                                                                            |
               !                       ( a2(i,j+1) -a4(i,j-1) ) *( b1(i+1,j) - b3(i-1,j) )  | 
               !                                                                            |
               !                                                                          --
               !
               zj1(xc,yc) = 1.0/(DxXDy) * ((var_a(xf,yc,k) - var_a(xb,yc,k)) * (var_b(xc,yf,k)-var_b(xc,yb,k))  &
                                        -  (var_a(xc,yf,k) - var_a(xc,yb,k)) * (var_b(xf,yc,k)-var_b(xb,yc,k)))
               !
               !                     --                                                                                             
               !             1       |                                                                                               
               ! J2     = ---------  | ( a1(i+1,j)) *( b5(i+1,j+1) - b8(i+1,j-1) )  -   ( a3(i-1,j)) *( b6(i-1,j+1) - b7(i-1,j-1) )  -
               !           4*Dx*Dy   |                                                                                               
               !                     -                                                                                               
               !
               !                                                                                                                     --
               !                                                                                                                      |
               !                       ( a2(i,j+1)) *( b5(i+1,j+1) - b6(i-1,j+1) )  -   ( a4(i,j-1) ) *( b8(i+1,j-1) - b7(i-1,j-1) )  |
               !                                                                                                                      |
               !                                                                                                                     -
               !
               zj2(xc,yc)= 1.0/(DxXDy) * ((var_b(xc,yf,k)*(var_a(xf,yf,k)-var_a(xb,yf,k))  -  var_b(xc,yb,k)*(var_a(xf,yb,k)-var_a(xb,yb,k)))   &
                                       -  (var_b(xf,yc,k)*(var_a(xf,yf,k)-var_a(xf,yb,k))  -  var_b(xb,yc,k)*(var_a(xb,yf,k)-var_a(xb,yb,k))) )

               !
               !                     --                                                                                               
               !             1       |                                                                                                
               ! J3     = ---------  | ( b2(i,j+1)) *( a5(i+1,j+1) - a6(i-1,j+1) )  -   ( b4(i,j-1) ) *( a8(i+1,j-1) - a7(i-1,j-1) ) - 
               !           4*Dx*Dy   |                                                                                                
               !                     -                                                                                                
               !
               !                                                                                                                      --
               !                                                                                                                       |
               !                       ( b1(i+1,j) ) *( a5(i+1,j+1) - a8(i+1,j-1) )  -   ( b3(i-1,j) ) *( a6(i-1,j+1) - a7(i-1,j-1) )  | 
               !                                                                                                                       |
               !                                                                                                                     --
               !
               zj3(xc,yc)= 1.0/(DxXDy) * ((var_a(xf,yc,k) * (var_b(xf,yf,k)-var_b(xf,yb,k))  - var_a(xb,yc,k) * (var_b(xb,yf,k)-var_b(xb,yb,k)))   &
                                       -  (var_a(xc,yf,k) * (var_b(xf,yf,k)-var_b(xb,yf,k))  - var_a(xc,yb,k) * (var_b(xf,yb,k)-var_b(xb,yb,k))) )


               pj(xc,yc,k)=(zj1(xc,yc)+zj2(xc,yc)+zj3(xc,yc))/3.0
           enddo
         enddo
      enddo

      return
      end subroutine jacobi
      
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

  SUBROUTINE index2(i,Idim,xb,xc,xf)
      IMPLICIT NONE
      INTEGER, INTENT(IN   ) :: i
      INTEGER, INTENT(IN   ) :: Idim
      INTEGER, INTENT(OUT  ) :: xb,xc,xf
      IF(i==1) THEN
        xb=i-1
        xc=i
        xf=i+1
      ELSE IF(i==Idim)THEN
        xb=Idim-1
        xc=Idim
        xf=Idim+1
      ELSE
        xb=i-1
        xc=i
        xf=i+1
      END IF
   END SUBROUTINE index2   
   SUBROUTINE Finalize_Class_Module_Dynamics() !FINALIZE
   IMPLICIT NONE
   
   END SUBROUTINE Finalize_Class_Module_Dynamics
   
END MODULE Class_Module_Dynamics
