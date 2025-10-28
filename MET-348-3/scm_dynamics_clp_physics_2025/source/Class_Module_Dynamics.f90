!  $Author: pkubota                                                 $
!  $Date: 2008/09/23 17:51:54                                         $
!  $Revision: 1.9                                                 $
!  $Revisions are currently made by the class's students.        $
!  $Update Date: 01/11/2023 10:19 AM                                $
!  
!  Implementações: 
!          1) Colocar INIT e FINALIZE - OK 
!         2) Colocar as equações dos campos no "Class_Module_Dynamics" - OK
!         7) Amortecimento das condições de contorno - OK (funcionando parcialmente)

MODULE Class_Module_Dynamics
 USE Constants, Only: r8,r4,i4,pi,Deg2Rad,Rd,Cp,kappa,r_earth,omega,CTv,nfprt, Eps
  
 
 USE Class_Module_Fields, Only: u_ref, v_ref,w_ref, t_ref,q_ref,z_ref,p_ref,&
                                vis,vis_q,taul,&
                                U_N,U_C,  V_N,V_C,  T_N,T_C,  Q_N,Q_C,&
                                r_LWRH,r_SWRH,r_LHCV,r_MSCV,r_LGLH,r_LGMS,&
                                r_SCVH,r_SCVM,r_PBLT,r_PBLQ,r_GDTZ,r_GDTM,&
                                Plevs,CoordLat,CoordLon,FcorPar ,DeltaLamda,DeltaTheta
  IMPLICIT NONE
  PRIVATE       

  INTEGER ::  Idim
  INTEGER ::  Jdim
  INTEGER ::  Kdim
  REAL(KIND=r8) :: DeltaT
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
              U_C(i,j,k)=  u_ref(i,j,k)
              V_C(i,j,k)=  v_ref(i,j,k)
              T_C(i,j,k)=  t_ref(i,j,k)
              Q_C(i,j,k)=  q_ref(i,j,k)
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
           END DO
        END DO
     END DO
  END IF
 !       PRINT*,'h eta',MAXVAL(h),MINVAL(h),'u eta',MAXVAL(u),MINVAL(u),'v eta',MAXVAL(v),MINVAL(v)
 
 call Solve_Forward_Beta_plane_Grid_A(ku1, kv1, kt1, kq1, uo        , vo             , to         , qo        )
 call Solve_Forward_Beta_plane_Grid_A(ku2, kv2, kt2, kq2, uo+ku1*dt2, vo+kv1*dt2     , to+kt1*dt2 , qo+kq1*dt2)
 call Solve_Forward_Beta_plane_Grid_A(ku3, kv3, kt3, kq3, uo+ku2*dt2, vo+kv2*dt2     , to+kt2*dt2 , qo+kq2*dt2)
 call Solve_Forward_Beta_plane_Grid_A(ku4, kv4, kt4, kq4, uo+ku3*dt , vo+kv3*dt      , to+kt3*dt  , qo+kq3*dt )
   
 ! final step and time marching / new values for RK4
 DO k=1,Kdim-1
    DO j=2,Jdim-1
       DO i=2,  Idim-1
          U_N(i,j,k) = uo(i,j,k) + (ku1(i,j,k) + 2.0*ku2(i,j,k) + 2.0*ku3(i,j,k) + ku4(i,j,k))*dt6
          V_N(i,j,k) = vo(i,j,k) + (kv1(i,j,k) + 2.0*kv2(i,j,k) + 2.0*kv3(i,j,k) + kv4(i,j,k))*dt6
          T_N(i,j,k) = to(i,j,k) + (kt1(i,j,k) + 2.0*kt2(i,j,k) + 2.0*kt3(i,j,k) + kt4(i,j,k))*dt6
          Q_N(i,j,k) = qo(i,j,k) + (kq1(i,j,k) + 2.0*kq2(i,j,k) + 2.0*kq3(i,j,k) + kq4(i,j,k))*dt6
       END DO
    END DO
 END DO

 ! updating the data
 U_C=u_ref
 V_C=v_ref
 T_C=t_ref
 Q_C=q_ref
 DO k=1,Kdim-1
    DO j=2,Jdim-1
       DO i=2,  Idim-1
          U_C(i,j,k) = U_N(i,j,k) 
          V_C(i,j,k) = V_N(i,j,k)
          T_C(i,j,k) = T_N(i,j,k)
          Q_C(i,j,k) = MAX(Q_N(i,j,k),1e-12)
       END DO
    END DO
 END DO
 PRINT*,it,MAXVAL(sqrt(U_C(2:Idim,2:Jdim,2:Kdim-1)**2 +V_C(2:Idim,2:Jdim,2:Kdim-1)**2)),&
           MAXVAL(t_C(2:Idim,2:Jdim,2:Kdim-1)),MINVAL(T_C(2:Idim,2:Jdim,2:Kdim-1)),& !TESTAR
           MAXVAL(Q_C(2:Idim,2:Jdim,2:Kdim-1)),MINVAL(Q_C(2:Idim,2:Jdim,2:Kdim-1))
! DO k=1,Kdim-1
!
! PRINT*,it,MAXVAL(T_C(2:Idim,2:Jdim,k)),MINVAL(T_C(2:Idim,2:Jdim,k)),&
!           MAXVAL(Q_C(2:Idim,2:Jdim,k)),MINVAL(Q_C(2:Idim,2:Jdim,k))
! END DO

 END SUBROUTINE RunDynamics
 
 
 SUBROUTINE  Solve_Forward_Beta_plane_Grid_A(TermEqMomU, TermEqMomV,TermEqConT, TermEqConQ, u_in, v_in, t_in,q_in)
  IMPLICIT NONE
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:kdim), intent(in) ::  u_in,  v_in,  t_in , q_in
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:kdim), intent(out) :: TermEqMomU,TermEqMomV, TermEqConT, TermEqConQ
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:kdim) ::  u,  v,  t , q
  REAL(KIND=8) :: udux, udvx, udqx, vdTx   
  REAL(KIND=8) :: vduy, vdvy, vdqy, udTy
  REAL(KIND=8) :: wduz, wdvz, wdqz, wdTz
  REAL(KIND=8) :: fcov
  REAL(KIND=8) :: hh  
  REAL(KIND=8) :: dPdx, dPdy ,dGeop1, dGeop2
  REAL(KIND=8) :: lnPs_xf
  REAL(KIND=8) :: lnPs_xb
  REAL(KIND=8) :: Ps_xf, q_zf
  REAL(KIND=8) :: Ps_xc, q_zb
  REAL(KIND=8) :: Tv, Tv_xf,Tv_xb, kTv, kTv0
  REAL(KIND=8) :: vis2dudx, vis2dvdx, vis2dqdx, vis2dTdx
  REAL(KIND=8) :: vis2dudy, vis2dvdy, vis2dqdy, vis2dTdy
  REAL(KIND=8) :: factor4,factor2
  REAL(KIND=8) :: term1,uadvc,vadvc,wadvc
  REAL(KIND=8) :: TermUNewton, TermVNewton, TermTNewton, TermQNewton, TermTNewtonAux, TermQNewtonAux
  INTEGER :: i,j,k
  INTEGER :: xb,xc,xf
  INTEGER :: yb,yc,yf
  TermEqMomU=0.0_r8;TermEqMomV=0.0_r8; TermEqConT=0.0_r8; TermEqConQ=0.0_r8
   TermTNewtonAux=0.0_r8; TermQNewtonAux=0.0_r8
  u=u_in;v=v_in;t=t_in;q=q_in
  
  !O QUE ISSO FAZ?
  
  u(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(u_ref(1:Idim   ,1:Jdim,Kdim) + u(1:Idim,1:Jdim,Kdim-1))

  v(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(v_ref(1:Idim   ,1:Jdim,Kdim) + v(1:Idim,1:Jdim,Kdim-1))

  t(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(t_ref(1:Idim   ,1:Jdim,Kdim)+ t(1:Idim,1:Jdim,Kdim-1))

  q(1:Idim   ,1:Jdim,Kdim)=0.5_r8*(q_ref(1:Idim   ,1:Jdim,Kdim)+ q(1:Idim,1:Jdim,Kdim-1))

  u(1:Idim   ,1:Jdim,1)=0.5_r8*(u_ref(1:Idim   ,1:Jdim,1) + u(1:Idim,1:Jdim,2))

  v(1:Idim   ,1:Jdim,1)=0.5_r8*(v_ref(1:Idim   ,1:Jdim,1) + v(1:Idim,1:Jdim,2))

  t(1:Idim   ,1:Jdim,1)=0.5_r8*(t_ref(1:Idim   ,1:Jdim,1) + t(1:Idim,1:Jdim,2))

  q(1:Idim   ,1:Jdim,1)=0.5_r8*(q_ref(1:Idim   ,1:Jdim,1) + q(1:Idim,1:Jdim,2))

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



     DO j=1,Jdim
        if(u(1   ,j,k) < 0.0)THEN
           q(1   ,j,k)=q_ref(1   ,j,k)                  - u(1   ,j,k)*&
           ((1.0_r8/(r_earth*(cos(CoordLat(1   ,j))**2))/(2*DeltaLamda(1,j)))&
                     *( q(3   ,j,k) - q_ref(1   ,j,k))) - v(1   ,j,k)*&
                     ((1.0_r8/(r_earth*(cos(CoordLat(1,j))))/(2*DeltaTheta(1,j)))&
                     *( q_ref(1   ,j,k) - q(3   ,j,k) ))
        else 
           q(1   ,j,k)=q_ref(1   ,j,k)                  - u(1   ,j,k)*&
           ((1.0_r8/(r_earth*(cos(CoordLat(1   ,j))**2))/2*DeltaLamda(1,j))&
                     *(q_ref(1   ,j,k) - q(3   ,j,k)))  - v(1   ,j,k)*&
                     ((1.0_r8/(r_earth*(cos(CoordLat(1  ,j))))/(2*DeltaTheta(1,j)))&
                     *( q(3   ,j,k)-q_ref(1   ,j,k)))
        end if
        if(u(Idim,j,k) <=0)THEN
           q(Idim,j,k)=q_ref(Idim,j,k)                   - u(Idim,j,k)*&
           ((1.0_r8/(r_earth*(cos(CoordLat(Idim,j))**2))/(2*DeltaLamda(Idim,j)))&
                     *(q_ref(Idim,j,k) - q(Idim-2,j,k))) - v(idim,j,k)*&
                     ((1.0_r8/(r_earth*(cos(CoordLat(idim,j))))/(2*DeltaTheta(idim,j)))&
                     *( q(idim-2,j,k)-q_ref(idim,j,k) ))
        else
           q(Idim,j,k)=q_ref(Idim,j,k)                  - u(Idim,j,k)*&
           ((1.0_r8/(r_earth*(cos(CoordLat(Idim,j))**2))/(2*DeltaLamda(Idim,j)))&
                     *( q(Idim-2,j,k)-q_ref(Idim,j,k))) - v(idim,j,k)*&
                     ((1.0_r8/(r_earth*(cos(CoordLat(idim,j))))/2*DeltaTheta(idim,j))&
                     *( q_ref(idim,j,k)-q(idim-2,j,k)))
        endif 
     END DO
     DO i=1,Idim
        if(v(i   ,1,k) > 0.0)THEN
           q(i ,1   ,k)=q_ref(i, 1   ,k)               - u(i   ,1,k)*&
           ((1.0_r8/(r_earth*(cos(CoordLat(i,1))**2))/(2*DeltaLamda(i,1)))&
                     *( q_ref(i,1,k) - q(i,3,k)))      - v(i   ,1,k)* &
                     ((1.0_r8/(r_earth*(cos(CoordLat(i,1))))/2*DeltaTheta(i,1))&
                     *( q_ref(i,1,k) - q(i,3,k)))
        else
           q(i ,1   ,k)=q_ref(i, 1   ,k)               - u(i   ,1,k)* &
           ((1.0_r8/(r_earth*(cos(CoordLat(i,1))**2))/(2*DeltaLamda(i,1)))&
                     *(q(i,3,k) - q_ref(i,1,k))) - v(i   ,1,k)*&
                      ((1.0_r8/(r_earth*(cos(CoordLat(i,1))))/(2*DeltaTheta(i,1)))&
                     *(q(i,3,k) - q_ref(i,1,k)))
        end if  
        if(v(i   ,jdim,k) > 0.0)THEN
           q(i ,jdim   ,k)=q_ref(i, jdim   ,k)           - u(i   ,jdim,k)*&
            ((1.0_r8/(r_earth*(cos(CoordLat(i,jdim))**2))/(2*DeltaLamda(i,jdim)))&
                     *(q_ref(i,jdim,k)-q(i,jdim-2,k)))   - v(i   ,jdim,k)*&
                      ((1.0_r8/(r_earth*(cos(CoordLat(i,jdim))))/(2*DeltaTheta(i,jdim)))&
                     *(q(i,jdim-2,k)-q_ref(i,jdim,k)))
        else
           q(i ,jdim   ,k)=q_ref(i, jdim   ,k)         - u(i,jdim,k)*&
           ((1.0_r8/(r_earth*(cos(CoordLat(i,jdim))**2))/2*DeltaLamda(i,jdim))&
                     *( q(i,jdim-2,k)-q_ref(i,jdim,k)))- v(i,jdim,k)*&
                     ((1.0_r8/(r_earth*(cos(CoordLat(i,jdim))))/(2*DeltaTheta(i,jdim)))&
                     *( q_ref(i,jdim,k)-q(i,jdim-2,k)))
        end if  
     END DO


     q(1   ,1:Jdim,k)=0.25_r8*(q_ref(1   ,1:Jdim,k) + q(2     ,1:Jdim,k) + q(3     ,1:Jdim,k) + q_ref(2     ,1:Jdim,k))
     q(Idim,1:Jdim,k)=0.25_r8*(q_ref(Idim,1:Jdim,k) + q(Idim-1,1:Jdim,k) + q(Idim-2,1:Jdim,k) + q_ref(Idim-1,1:Jdim,k))

     q(1:Idim,1   ,k)=0.25_r8*(q_ref(1:Idim,1   ,k) + q(1:Idim,2     ,k) + q(1:Idim,3     ,k) + q_ref(1:Idim,2     ,k))
     q(1:Idim,Jdim,k)=0.25_r8*(q_ref(1:Idim,Jdim,k) + q(1:Idim,Jdim-1,k) + q(1:Idim,Jdim-2,k) + q_ref(1:Idim,Jdim-1,k))


  END DO
  
  !PARTE ESPACIAL
  DO k=1,Kdim-1 
      DO j=2,Jdim-1
         CALL index(j,Jdim,yb,yc,yf)
         DO i=2,Idim-1
            CALL index(i,Idim,xb,xc,xf)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*            
            !
            !TERMO DA ZONAL
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*            

            TermUNewton = (u(xc,yc,k) - u_ref(xc,yc,k))/taul !Amortecimento para as forçantes
            !
            !                       --               --   
            !                      |                   |  
            !            1         |          du       |  
            !udux= ----------------|   U * ----------  |  
            !       a*cos^2(theta) |         d lambda  |  
            !                      |                   |  
            !                       --               --   
            uadvc = (1.0_r8/6.0_r8) * (u(xf,yc,k) + u(xc,yc,k) + u(xc,yf,k) + u(xc,yb,k) + u(xc,yc,k) + u(xb,yc,k))
   
            udux = (1.0_r8/(r_earth*(cos(CoordLat(xc,yc))**2))) * &
                     (uadvc *((u(xf,yc,k) - u(xb,yc,k))/(2_r8*DeltaLamda(xc,yc)))) 

            !
            !                       --                --
            !                       |                   |
            !             1         |          du       | 
            !vduy=  ----------------|   V * ----------  | 
            !        a*cos(theta)   |         d theta   | 
            !                       |                   |
            !                       --                 --
            vadvc = (1.0_r8/6.0_r8)*(v(xf,yc,k)+v(xc,yc,k)+v(xc,yf,k)+v(xc,yb,k)+v(xc,yc,k)+v(xb,yc,k))

            vduy  = (1.0_r8/(r_earth*cos(CoordLat(xc,yc))))  * &
                    (vadvc*((u(xc,yf,k) - u(xc,yb,k))/(2.0_r8*DeltaTheta(xc,yc))))

            !      --                --
            !      |                   |
            !      |          du       | 
            !wduz= |   w * ----------  | 
            !      |         d P       | 
            !      |                   |
            !      --                 --
            wadvc = 0.05_r8*(0.25_r8*(w_ref (xf,yc,k  ) + w_ref (xb,yc,k  ) +w_ref (xc,yf,k  ) + w_ref (xc,yb,k  )) + &
                             0.25_r8*(w_ref (xf,yc,k+1) + w_ref (xb,yc,k+1) +w_ref (xc,yf,k+1) + w_ref (xc,yb,k+1)))
            Ps_xf=Plevs(k+1)!(p_ref(xc,yc)*(Plevs(k+1)/100000.0_r8))
            Ps_xc=Plevs(k+0)!(p_ref(xc,yc)*(Plevs(k+0)/100000.0_r8))
            wduz = wadvc * ((u(xc,yc,k+1) - u(xc,yc,k))/(Ps_xf-Ps_xc))

            !      --       --
            !      |         |
            !      |         | 
            !fcov= | f  * v  | 
            !      |         | 
            !      |         |
            !      --      --
    
            fcov  = - FcorPar(xc,yc)*v(xc,yc,k)
            !
            !           --                            --
            !          |                                |
            !        1 |       dZGeo           dln(P)   | 
            !dPdx=  ---|    ---------- + Rd*Tv--------  | 
            !        a |     dLambda           dLambda  | 
            !          |                                |
            !          --                             --
            Tv_xf=t(xf,yc,k)*(1.0_r8 + CTv*q(xf,yc,k))
            Tv   =t(xc,yc,k)*(1.0_r8 + CTv*q(xc,yc,k))
            Tv_xb=t(xb,yc,k)*(1.0_r8 + CTv*q(xb,yc,k))
            lnPs_xf=log(Plevs(k+0))! (p_ref(xf,yc)*(Plevs(k+0)/100000.0_r8)) )  
            lnPs_xb=log(Plevs(k+0))! (p_ref(xb,yc)*(Plevs(k+0)/100000.0_r8)) )
            dGeop2 = z_ref(xf,yc,k)*(1.0_r8 + 0.0001*Tv_xf)
            dGeop1 = z_ref(xb,yc,k)*(1.0_r8 + 0.0001*Tv_xb)
            dPdx =(1.0_r8/(r_earth)) * (((dGeop2 - dGeop1)/(2.0*DeltaLamda(xc,yc))) + &
                                          Rd*Tv*((lnPs_xf-lnPs_xb)/(DeltaLamda(xc,yc))))
          
            !Testar com e sem a difusão
            !                  -         -
            ! d(u)            |d(d(u))    |
            ! -----  - Neta * |--------   | = 0
            ! dt              |dxdx       |
            !                  -         -
            factor4=(r_earth**2)*(cos(CoordLat(xc,yc))**4)
            vis2dudx= -vis*(1/factor4)*((u(xf,yc,k) - 2.0*u(xc,yc,k) + u(xb,yc,k)) &
                            /(DeltaLamda(xc,yc)*DeltaLamda(xc,yc)))
            !
            !                  -       -
            ! d(u)             |d(d(u)) |
            ! -----   - Neta * |--------| = 0
            ! dt               |dydy    |
            !                  -       -
            factor2=(r_earth**2)*(cos(CoordLat(xc,yc))**2)
            vis2dudy= -vis*(1.0/factor2)*((u(xc,yf,k) - 2.0*u(xc,yc,k) + u(xc,yb,k)) & 
                            /(DeltaTheta(xc,yc)*DeltaTheta(xc,yc)))

            TermEqMomU(xc,yc,k) = -( udux + vduy + wduz + fcov + dPdx + vis2dudx + vis2dudy + TermUNewton)
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*            
            !
            !TERMO DA MERIDIONAL
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*            
            
            TermVNewton = (v(xc,yc,k) - v_ref(xc,yc,k))/taul !Amortecimento para as forçantes           
            !                       --               --   
            !                      |                   |  
            !            1         |          dv       |  
            !udvx= ----------------|   U * ----------  |  
            !       a*cos^2(theta) |         d lambda  |  
            !                      |                   |  
            !                       --               --      
            udvx = (1.0_r8/(r_earth*(cos(CoordLat(xc,yc))**2))) * &
                     (uadvc *((v(xf,yc,k) - v(xb,yc,k))/(2_r8*DeltaLamda(xc,yc)))) 
            !                       --                --
            !                       |                   |
            !             1         |          dv       | 
            !vdvy=  ----------------|   V * ----------  | 
            !        a*cos(theta)   |         d theta   | 
            !                       |                   |
            !                       --                 --
            vadvc = (1.0_r8/6.0_r8)*(v(xf,yc,k) + v(xc,yc,k) + v(xc,yf,k) + &
                            v(xc,yb,k) + v(xc,yc,k) + v(xb,yc,k))

            vdvy  = (1.0_r8/(r_earth*cos(CoordLat(xc,yc))))  * &
                    (vadvc*((v(xc,yf,k) - v(xc,yb,k))/(2.0_r8*DeltaTheta(xc,yc))))

            !      --                --
            !      |                   |
            !      |          dv       | 
            !wdvz= |   w * ----------  | 
            !      |         d P       | 
            !      |                   |
            !      --                 --
            wadvc = 0.05_r8*(0.25_r8*(w_ref (xf,yc,k  ) + w_ref (xb,yc,k  ) +w_ref (xc,yf,k  ) + w_ref (xc,yb,k  )) + &
                             0.25_r8*(w_ref (xf,yc,k+1) + w_ref (xb,yc,k+1) +w_ref (xc,yf,k+1) + w_ref (xc,yb,k+1)))
            Ps_xf=Plevs(k+1)!(p_ref(xc,yc)*(Plevs(k+1)/100000.0_r8))
            Ps_xc=Plevs(k+0)!(p_ref(xc,yc)*(Plevs(k+0)/100000.0_r8))

            wdvz = wadvc * ((v(xc,yc,k+1) - v(xc,yc,k))/(Ps_xf-Ps_xc))

            !      --       --
            !      |         |
            !      |         | 
            !fcov= | f  * u  | 
            !      |         | 
            !      |         |
            !      --      --
    
            fcov  =  FcorPar(xc,yc)*u(xc,yc,k)
            !
            !           --                            --
            !                   |                                |
            !        cos(theta) |       dZGeo           dln(P)   | 
            !dPdy=   ---------  |    ---------- + Rd*Tv--------  | 
            !            a      |     d theta           d theta  | 
            !                   |                                |
            !          --                             --
            Tv_xf=t(xc,yf,k)*(1.0_r8 + CTv*q(xc,yf,k))
            Tv   =t(xc,yc,k)*(1.0_r8 + CTv*q(xc,yc,k))
            Tv_xb=t(xc,yb,k)*(1.0_r8 + CTv*q(xc,yb,k))
            lnPs_xf=log(Plevs(k+0))! (p_ref(xc,yf)*(Plevs(k+0)/100000.0_r8)) )  
            lnPs_xb=log(Plevs(k+0))! (p_ref(xc,yb)*(Plevs(k+0)/100000.0_r8)) )
            dGeop2 = z_ref(xc,yf,k)*(1.0_r8 + 0.0001*Tv_xf)
            dGeop1 = z_ref(xc,yb,k)*(1.0_r8 + 0.0001*Tv_xb)
            dPdy =(cos(CoordLat(xc,yc))/(r_earth)) * (((dGeop2 - dGeop1)/(2.0*DeltaTheta(xc,yc))) + &
                                          Rd*Tv*((lnPs_xf-lnPs_xb)/(2.0_r8*DeltaTheta(xc,yc)))) 
            !                  -         -
            ! HH=    sin(theta)      |      U^2 + V^2   |
            !        -----------     |      ---------   | = 0
            !      a*cos^2(theta)    |          1       |
            !                  -         -
        
            hh=(sin(CoordLat(xc,yc))/(r_earth*cos(CoordLat(xc,yc))**2)) * ((u(xc,yc,k)**2) + (v(xc,yc,k)**2)) 
            
            !Testar com e sem a difusão
            
            !                  -         -
            ! d(v)            |d(d(v))    |
            ! -----  - Neta * |--------   | = 0
            ! dt              |dxdx       |
            !                  -         -         
            vis2dvdx= -vis*(1/factor4)*((v(xf,yc,k) - 2.0*v(xc,yc,k) + v(xb,yc,k))/(DeltaLamda(xc,yc)*DeltaLamda(xc,yc)))
            !
            !                  -       -
            ! d(v)             |d(d(v)) |
            ! -----   - Neta * |--------| = 0
            ! dt               |dydy    |
            !                  -       -
            vis2dvdy= -vis*(1.0/factor2)*((v(xc,yf,k) - 2.0*v(xc,yc,k) + v(xc,yb,k))/(DeltaTheta(xc,yc)*DeltaTheta(xc,yc)))

            TermEqMomV(xc,yc,k) = -( udvx + vdvy + wdvz + fcov + dPdy + hh + vis2dvdx + vis2dvdy + TermVNewton)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*            
            !
            !TERMO DA TEMPERATURA:
            !
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*            

            TermTNewtonAux = (t(xc,yc,k) - t_ref(xc,yc,k))/taul !Amortecimento para as forçantes
            
            !TermTNewton = 0.4*TermTNewtonAux 

            TermTNewton =  (r_LWRH(k))*TermTNewtonAux + (r_SWRH(k))*TermTNewtonAux + &
                           (r_LHCV(k))*TermTNewtonAux + (r_LGLH(k))*TermTNewtonAux + &
                           (r_SCVH(k))*TermTNewtonAux + (r_PBLT(k))*TermTNewtonAux !Amortecimento para as forçantes
            !                       --               --   
            !                      |                   |  
            !            1         |          dT       |  
            !udT= ----------------|   U * ----------  |  
            !       a*cos^2(theta) |         d lambda  |  
            !                      |                   |  
            !                       --               --   
   
            udTy = (1.0_r8/(r_earth*(cos(CoordLat(xc,yc))**2))) * &
                           (uadvc *((t(xf,yc,k) - t(xb,yc,k))/(2_r8*DeltaLamda(xc,yc)))) 
  
            !
            !                       --                --
            !                       |                   |
            !             1         |          dT       | 
            !vdTx=  ----------------|   V * ----------  | 
            !        a*cos(theta)   |         d theta   | 
            !                       |                   |
            !                       --                 --
  
            vdTx  = (1.0_r8/(r_earth*cos(CoordLat(xc,yc))))  * &
                    (vadvc*((t(xc,yf,k) - t(xc,yb,k))/(2.0_r8*DeltaTheta(xc,yc))))
  
            !      --                --
            !      |                   |
            !      |          dT       | 
            !wdTz= |   w * ----------  | 
            !      |         d P       | 
            !      |                   |
            !      --                 --
  
            wdTz = wadvc * ((t(xc,yc,k+1) - t(xc,yc,k))/(Ps_xf-Ps_xc))
            
            !                                   --                
            !                                     |  
            !            kappa*T_v*omega          |            
            !kTv= - -------------------------     |    
            !           (1 + (delta-1)q)p         |          
            !                                     |                    
            !--                                  --  
            kTv0 = (1 + (Eps - 1)*q(xc,yc,k) )*Plevs(k)
            kTv  = - Tv*kappa*w_ref(i,j,k) / kTv0

            !      --       --            !                 
            ! d(T)            |d(d(T))    |
            ! -----  - Neta * |--------   | = 0
            ! dt              |dxdx       |
            !                  -         -
            vis2dtdx= -vis*(1/factor4)*((t(xf,yc,k) - 2.0*t(xc,yc,k) + t(xb,yc,k))/(DeltaLamda(xc,yc)*DeltaLamda(xc,yc)))
            !
            !                  -       -
            ! d(T)             |d(d(T)) |
            ! -----   - Neta * |--------| = 0
            ! dt               |dydy    |
            !                  -       -

            vis2dtdy= -vis*(1.0/factor2)*((t(xc,yf,k) - 2.0*t(xc,yc,k) + t(xc,yb,k))/(DeltaTheta(xc,yc)*DeltaTheta(xc,yc)))

            TermEqConT(xc,yc,k) = -( vdTx + udTy + wdTz + kTv + vis2dtdx + vis2dtdy + TermTNewton)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*            
            !
            !TERMO DA UMIDADE
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*            
                                             
            TermQNewtonAux = (q(xc,yc,k) - q_ref(xc,yc,k))/taul !Amortecimento para as forcantes
            
            TermQNewtonAux = TermQNewtonAux
            TermQNewton =  (r_SCVM(k))*TermQNewtonAux +  (r_MSCV(k))*TermQNewtonAux + &
                           (r_LGMS(k))*TermQNewtonAux +  (r_PBLQ(k))*TermQNewtonAux

            !
            !                       --               --   
            !                      |                   |  
            !            1         |          dq       |  
            !udqx= ----------------|   U * ----------  |  
            !       a*cos^2(theta) |         d lambda  |  
            !                      |                   |  
            !                       --               --   
            udqx = (1.0_r8/(r_earth*(cos(CoordLat(xc,yc))**2))) * &
                     (uadvc *((q(xf,yc,k) - q(xb,yc,k))/(2_r8*DeltaLamda(xc,yc)))) 

            !
            !                       --                --
            !                       |                   |
            !             1         |          dq       | 
            !vdqy=  ----------------|   V * ----------  | 
            !        a*cos(theta)   |         d theta   | 
            !                       |                   |
            !                       --                 --
            vdqy  = (1.0_r8/(r_earth*cos(CoordLat(xc,yc))))  * &
                    (vadvc*((q(xc,yf,k) - q(xc,yb,k))/(2.0_r8*DeltaTheta(xc,yc))))

            !      --                --
            !      |                   |
            !      |          dq       | 
            !wdvz= |   w * ----------  | 
            !      |         d P       | 
            !      |                   |
            !      --                 --
            q_zb=0.25_r8*(q (xf,yc,k  ) + q (xb,yc,k  ) + q (xc,yf,k  ) + q (xc,yb,k  ))
            q_zf=0.25_r8*(q (xf,yc,k+1) + q (xb,yc,k+1) + q (xc,yf,k+1) + q (xc,yb,k+1))
            wdqz = 0.25*wadvc * (( q_zf - q_zb)/(Ps_xf-Ps_xc))
            
            !      --       --            !                 
            ! d(q)            |d(d(q))    |
            ! -----  - Neta * |--------   | = 0
            ! dt              |dxdx       |
            !                  -         -
            vis2dqdx= -vis_q*(1/factor4)*((q(xf,yc,k) - 2.0*q(xc,yc,k) + q(xb,yc,k))/(DeltaLamda(xc,yc)*DeltaLamda(xc,yc)))
            !
            !                  -       -
            ! d(q)             |d(d(q)) |
            ! -----   - Neta * |--------| = 0
            ! dt               |dydy    |
            !                  -       -

            vis2dqdy= -vis_q*(1.0/factor2)*((q(xc,yf,k) - 2.0*q(xc,yc,k) + q(xc,yb,k))/(DeltaTheta(xc,yc)*DeltaTheta(xc,yc)))

            TermEqConQ(xc,yc,k) = -(udqx + vdqy + wdqz + vis2dqdx + vis2dqdy + TermQNewton)
                        
         END DO
      END DO
  END DO


 END SUBROUTINE  Solve_Forward_Beta_plane_Grid_A

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
