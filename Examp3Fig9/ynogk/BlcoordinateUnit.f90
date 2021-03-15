



!*******************************************************************************
      module BLcoordinate
!*******************************************************************************
!*     PURPOSE: This module aims on computing 4 Boyer-Lindquist coordinates (r,\theta,\phi,t)
!*              and affine parameter \sigam.    
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012 
!***********************************************************************
      use constants
      use rootsfinding
      use ellfunction      
      implicit none

      contains
!******************************************************************************************** 
      SUBROUTINE YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        radi,mu,phi,time,sigma,sign_pr,sign_pth) 
!********************************************************************************************
!*     PURPOSE:  Computes four Boyer-Lindquist coordinates (r,\mu,\phi,t) and affine parameter 
!*               \sigma as functions of parameter p, i.e. functions r(p), \mu(p), \phi(p), t(p)
!*               and \sigma(p). Cf. discussions in Yang & Wang (2012).    
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f1234----------array of p_1, p_2, p_3, p_4, which are the components of four-
!*                              momentum of a photon measured under the LNRF frame. This array 
!*                              can be computed by subroutine lambdaq(...), see below.     
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or initialposition of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  radi-----------value of function r(p). 
!*               mu-------------value of function \mu(p). 
!*               phi------------value of function \phi(p).
!*               time-----------value of function t(p).
!*               sigma----------value of function \sigma(p).
!*               sign_pr--------sign of r component of 4-momentum of the photon.
!*               sign_pth-------sign of \theta component of 4-momentum of the photon.
!*               tm1,tm2--------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively. 
!*               tr1,tr2--------number of times of photon meets turning points r_tp1 and r_tp2
!*                              respectively.            
!*     ROUTINES CALLED: INTRPART, INTTPART.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ********************************************************************
        USE constants 
        IMPLICIT NONE 
        DOUBLE PRECISION :: f1234(4),lambda,q,sinobs,muobs,a_spin,robs,scal,&
                phi_r,time_r,aff_r,phi_t,time_t,&
                rp,mup,p_int,mu_cos,r_coord,radi,mu,time,phi,sigma,p,&
                mu_tp,mu_tp2,Rab,sign_pr,sign_pth
        CHARACTER :: varble
        LOGICAL :: rotate,err,mobseqmtp
        INTEGER :: tm1,tm2,tr1,tr2,reals,del,t1,t2 

        !************************************************************************************
! call integrat_r_part to evaluate t_r,\phi_r,\sigma_r, and function r(p) (here is r_coord).
        call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,&
                        scal,phi_r,time_r,aff_r,r_coord,tr1,tr2,sign_pr) 
! call integrat_theta_part to evaluate t_\mu,\phi_\mu,\sigma_\mu, and function \mu(p) (here is mu_cos).
        call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,&
                        scal,phi_t,time_t,mu_cos,tm1,tm2,sign_pth) 
!write(*,*)'ynogk===', tm1, tm2
        radi=r_coord
        mu=mu_cos
!time coordinate value, equation (74) of Yang & Wang (2012).
        time=time_r+time_t
!affine parameter value, equation (74) of Yang & Wang (2012).
! write(*,*)'ynogk=',aff_r,time_t,time_r!,time_t
        sigma=aff_r+time_t
!phi coordinate value.
        rotate=.false.
        err=.false.
!write(*,*)'phi2=',p,phi_r,phi_t,tm1,tm2,time_r,time_t,lambda,f1234(3)
        IF(DABS(muobs).NE.ONE)THEN
! equation (74) of Yang & Wang (2012).
            phi=-(phi_r+phi_t)
            IF(f1234(3).EQ.zero)THEN 
                phi=phi+(tm1+tm2)*PI
            ENDIF
            phi=DMOD(phi,twopi)
            IF(phi.LT.zero)THEN
                phi=phi+twopi
            ENDIF
        ELSE 
! equation (74) of Yang & Wang (2012).
            phi=-(phi_t+phi_r+(tm1+tm2)*PI)
            Rab=dsqrt(f1234(3)**two+f1234(2)**two)
            IF(phi.NE.zero)THEN
                rotate=.TRUE.
            ENDIF
            IF(Rab.NE.zero)THEN
! a muobs was multiplied to control the rotate direction
                if((f1234(3).ge.zero).and.(f1234(2).gt.zero))then
                    phi=muobs*phi+dasin(f1234(2)/Rab)  
                endif
                if((f1234(3).lt.zero).and.(f1234(2).ge.zero))then
                    phi=muobs*phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).le.zero).and.(f1234(2).lt.zero))then
                    phi=muobs*phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).gt.zero).and.(f1234(2).le.zero))then
                    phi=muobs*phi+twopi+dasin(f1234(2)/Rab)
                endif
            ELSE
                phi=zero
            ENDIF
            IF(rotate)THEN
                phi=dMod(phi,twopi)
                IF(phi.LT.zero)THEN
                   phi=phi+twopi
                ENDIF
            ENDIF
        ENDIF        
        RETURN
      END SUBROUTINE YNOGK

!============================================================================================
      Function mucos(p,f12343,f12342,lambda,q,sinobs,muobs,a_spin,scal, sign_pth)
!============================================================================================
!*     PURPOSE:  Computes function \mu(p) defined by equation (32) in Yang & Wang (2012). That is
!*               \mu(p)=b0/(4*\wp(p+PI0;g_2,g_3)-b1)+\mu_tp1. \wp(p+PI0;g_2,g_3) is the Weierstrass'
!*               elliptic function.  
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f12342---------p_2, \theta component of four momentum of a photon measured under a LNRF.
!*               f12343---------p_3, \phi component of four momentum of a photon measured under a LNRF..
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  mucos----------\mu coordinate of photon corresponding to a given p. 
!*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision mucos,f12343,f12342,p,sinobs,muobs,a_spin,q,lambda,bc,cc,dc,ec,mu_tp,&
                        zero,b0,b1,b2,b3,g2,g3,tinf,fzero,a4,b4,AA,BB,two,delta,scal,four,&
                        mu_tp2,one,three,integ4(4),rff_p, pp0, u
      Double precision :: half_period_wp, period_wp
      Double precision, intent(out) :: sign_pth
      Double precision f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1
      complex*16 dd(3)
      integer ::  reals,p4,index_p4(4),del,cases,count_num=1
      logical :: mobseqmtp
      save  f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1,&
                mu_tp, mu_tp2,b0,b1,b2,b3,g2,g3,dd,fzero,count_num,AA,BB,del, &
                half_period_wp, period_wp, pp0
      parameter (zero=0.D0,two=2.0D0,four=4.D0,one=1.D0,three=3.D0)

10    continue
      If(count_num.eq.1)then
          f12343_1=f12343
          f12342_1=f12342
          lambda_1=lambda
          q_1=q
          muobs_1=muobs
          sinobs_1=sinobs
          a_spin_1=a_spin 
          scal_1=scal
   !*****************************************************************************************************
          If(f12343.eq.zero.and.f12342.eq.zero.and.dabs(muobs).eq.one)then
              sign_pth = zero
              mucos=muobs     !this is because that mu==1 for ever,this because that Theta_mu=-a^2(1-mu^2)
              return          !so,mu must =+1 or -1 for ever. q=-a^2, X=lambda/sin(theta)=0
          endif               !so Theta_mu=q+a^2mu^2-X^2mu^4=-a^2(1-mu^2) 
! spin is zero.
          If(a_spin.eq.zero)then
              if(q.gt.zero)then        
                  AA=dsqrt((lambda**two+q)/q)
                  BB=dsqrt(q)
                  If(f12342.lt.zero)then
                      u = dasin(muobs*AA)+p*BB*AA  
                      !mucos=dsin(dasin(muobs*AA)+p*BB*AA)/AA
                      mucos = dsin( u ) / AA
                      u = dmod(u + halfpi, twopi)
                      If ( zero <= u .AND. u <= pi ) then
                          sign_pth = -one
                      else
                          sign_pth = one
                      endif
                  else
                      If(f12342.eq.zero)then
                          u = p*AA*BB
                          mucos=dcos(u)*muobs
                          u = dmod(u, twopi)
                          If ( zero <= u .AND. u <= pi ) then
                              sign_pth = dsign(one, muobs)
                          else
                              sign_pth = -dsign(one, muobs)
                          endif  
                      else
                          u = dasin(muobs*AA)-p*BB*AA                     
                          !mucos=dsin(dasin(muobs*AA)-p*AA*BB)/AA        
                          mucos = dsin( u ) / AA 
                          u = dmod(u - halfpi, twopi)
                          If ( -pi <= u .AND. u <= zero ) then
                              sign_pth = one
                          else
                              sign_pth = -one
                          endif       
                      endif        
                  endif
              else
                  sign_pth = zero
                  mucos=muobs
              endif
          else
! Equatorial plane motion.
              If(muobs.eq.zero.and.q.eq.zero)then
                  sign_pth = zero
                  mucos=zero
                  return                        
              endif
              call mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,&
                                       mu_tp,mu_tp2,reals,mobseqmtp)
              a4=zero
              b4=one
              p4=0
! equations (26)-(29) in Yang & Wang (2012).
              b0=-four*a_spin**2*mu_tp**3+two*mu_tp*(a_spin**2-lambda**2-q)
              b1=-two*a_spin**2*mu_tp**2+one/three*(a_spin**2-lambda**2-q)
              b2=-four/three*a_spin**2*mu_tp
              b3=-a_spin**2
! equation (31) in Yang & Wang (2012).
              g2=three/four*(b1**2-b0*b2)
              g3=one/16.D0*(three*b0*b1*b2-two*b1**3-b0**2*b3)

              call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
              index_p4(1)=0        
              cases=1
! equation (33) in Yang & Wang (2012).
              If(muobs.ne.mu_tp)then        
                  tinf=b0/( - four*dabs(muobs-mu_tp))+b1/four
                  call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,&
                                          index_p4,rff_p,integ4,cases)
                  fzero=integ4(1)   
                  tinf=b0/(four*(mu_tp2-mu_tp))+b1/four
                  call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,&
                                          index_p4,rff_p,integ4,cases)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                  half_period_wp = integ4(1)
                  !pp0 = halfperiodwp(g2, g3, dd, del) 
                  period_wp = two * half_period_wp
         !write(unit = *, fmt = *)'************************************************************'  
         !write(unit = *, fmt = 2001) half_period_wp, pp0, half_period_wp - pp0
         !write(unit = *, fmt = *)'************************************************************' 
         !2001 Format(3F50.40)
         !3001 Format(2F50.40) 
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              else
                  fzero=zero                
              endif
              If(f12342.lt.zero)then
                  fzero = - fzero         
              endif
! equation (32) in Yang & Wang (2012).
              u = p + fzero
              mucos=mu_tp+b0/(four*weierstrassP(u, g2,g3,dd,del)-b1)
              if ( u <= zero ) then
                  sign_pth = -one
              else
                  u = dmod(u, period_wp)
                  if ( u <= half_period_wp ) then
                      sign_pth = one
                  else
                      sign_pth = - one
                  endif
              endif
                  if(scal == -two)then
                      write(unit = *, fmt = *)'**uuuuuu1== ', p + fzero, u, half_period_wp
                      write(unit = *, fmt = *)'**uuuuuu2== ', b0, b1, mu_tp, mucos, sign_pth
                  endif
                  if(mucos .lt. mu_tp2)then
                      write(unit = *, fmt = *)'**mucos ddd== === ', mucos, mu_tp2, mu_tp2 - mucos
                  endif
                  if (scal == -one)then  
      write(unit = *, fmt = *)'************************************************************'  
      write(unit = *, fmt = *)'**mucos dddd Savedd computing ===== === ', &
            b0, b1, p, fzero  
      write(unit = *, fmt = *)'**mucos dddd Savedd computing ===== === ', &
            b0/(four*(muobs-mu_tp))+b1/four, mu_tp, mu_tp2, muobs, f12343,f12342
      write(unit = *, fmt = *)'**mucos dddd Savedd computing ===== === ',&
           halfperiodwp(g2,g3,dd,del), pp0, weierstrassP((p+fzero),g2,g3,dd,del), &
             weierstrassP(halfperiodwp(g2,g3,dd,del),g2,g3,dd,del)
      write(unit = *, fmt = *)'**mu_tp2 ===== === ', &  
            p+fzero, halfperiodwp(g2,g3,dd,del), p+fzero - halfperiodwp(g2,g3,dd,del)
      write(unit = *, fmt = *)'**mu_tp2 ===== === ', &  
           mucos, mu_tp + b0 / (four*weierstrassP(halfperiodwp(g2,g3,dd,del),g2,g3,dd,del)-b1), mu_tp2 
      write(unit = *, fmt = *)'************************************************************'
                  tinf=b0/(four*dabs(muobs-mu_tp))+b1/four
                  call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,&
                                          index_p4,rff_p,integ4,cases)
      write(*,*)'55==', tinf, infinity,dd,del,a4,b4,&
                                          index_p4,rff_p,integ4(1),cases
                  endif
              !write(*,*)'mu=',weierstrassP(p+fzero,g2,g3,dd,del),b0,mu_tp,b1!tinf,infinity,g2,g3,a4,b4,p4,fzero
              ! If muobs eq 0,q eq 0,and mu_tp eq 0,so b0 eq 0,
              ! so mucos eq mu_tp eq 0.
              count_num=count_num+1
          endif
 !**************************************************************************       
      else
          If(f12343.eq.f12343_1.and.f12342.eq.f12342_1.and.lambda.eq.lambda_1.and.q.eq.q_1.and.&
              sinobs.eq.sinobs_1.and.muobs.eq.muobs_1.and.a_spin.eq.a_spin_1&
              .and.scal.eq.scal_1)then
              !******************************************************************        
              If(f12343.eq.zero.and.f12342.eq.zero.and.dabs(muobs).eq.one)then
                  sign_pth = zero
                  mucos=muobs                 !this is because that mu==1 for ever,
                                              !this is because that Theta_mu=-a^2(1-mu^2)
                  return                 !so,mu must =+1 or -1 for ever. q=-a^2, X=lambda/sin(theta)=0
              endif
              If(a_spin.eq.zero)then
                  if(q.gt.zero)then        
                      If(f12342.lt.zero)then
                          u = dasin(muobs*AA) + p*BB*AA
                          !mucos = dsin(dasin(muobs*AA)+p*BB*AA)/AA
                          mucos = dsin( u ) / AA
                          u = dmod(u + halfpi, twopi)
                          If ( zero <= u .AND. u <= pi ) then
                              sign_pth = -one
                          else
                              sign_pth = one
                          endif
                      else
                          If(f12342.eq.zero)then
                              u = p*AA*BB
                              !mucos=dcos(p*AA*BB)*muobs
                              mucos=dcos(u)*muobs
                              u = dmod(u, twopi)
                              If ( zero <= u .AND. u <= pi ) then
                                  sign_pth = dsign(one, muobs)
                              else
                                  sign_pth = -dsign(one, muobs)
                              endif  
                          else
                              u = dasin(muobs*AA) - p*BB*AA
                              !mucos=dsin(dasin(muobs*AA)-p*AA*BB)/AA        
                              mucos = dsin( u ) / AA 
                              u = dmod(u - halfpi, twopi)
                              If ( -pi <= u .AND. u <= zero ) then
                                  sign_pth = one
                              else
                                  sign_pth = -one
                              endif  
                          endif        
                      endif
                  else
                      mucos=muobs
                  endif
              else
                  If(muobs.eq.zero.and.q.eq.zero)then
                      mucos=zero
                      sign_pth = zero
                      return                        
                  endif    
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                  !write(unit = *, fmt = *)'************************************************************'  
                  !write(unit = *, fmt = 2001) half_period_wp, pp0, half_period_wp -pp0
                  !write(unit = *, fmt = *)'************************************************************' 
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! equation (32) in Yang & Wang (2012).    
                  u = p + fzero
                  mucos = mu_tp + b0 / (four*weierstrassP(u, g2,g3,dd,del)-b1) 
                  if ( u <= zero ) then
                      sign_pth = -one
                  else
                      u = dmod(u, period_wp)
                      if ( u <= half_period_wp ) then
                          sign_pth = one
                      else
                          sign_pth = - one
                      endif
                  endif
                  if(scal == -two)then
                      write(unit = *, fmt = *)'**uuuuuu1== ', p + fzero, u, half_period_wp
                      write(unit = *, fmt = *)'**uuuuuu2== ', b0, b1, mu_tp, mucos, sign_pth
                  endif
                  if(mucos .lt. mu_tp2)then
                      write(unit = *, fmt = *)'**mucos ddd== === ', mucos, mu_tp2, mu_tp2- mucos
                  endif
                  if (scal == -one)then 
      write(unit = *, fmt = *)'************************************************************'  
      write(unit = *, fmt = *)'**mucos dddd Savedd computing ===== === ', &
            b0, b1, p, fzero
      write(unit = *, fmt = *)'**mucos dddd Savedd computing ===== === ', &
            b0/(four*(muobs-mu_tp))+b1/four, mu_tp, mu_tp2, muobs, f12343,f12342
      write(unit = *, fmt = *)'**mucos dddd Savedd computing ===== === ',&
           halfperiodwp(g2,g3,dd,del), pp0, weierstrassP((p+fzero),g2,g3,dd,del), &
             weierstrassP(halfperiodwp(g2,g3,dd,del),g2,g3,dd,del)
      write(unit = *, fmt = *)'**mu_tp2 ===== === ', &  
            p+fzero, halfperiodwp(g2,g3,dd,del), p+fzero - halfperiodwp(g2,g3,dd,del)
      write(unit = *, fmt = *)'**mu_tp2 ===== === ', &  
           mucos, mu_tp + b0 / (four*weierstrassP(halfperiodwp(g2,g3,dd,del),g2,g3,dd,del)-b1), mu_tp2 
      write(unit = *, fmt = *)'************************************************************'
                  tinf=b0/(four*dabs(muobs-mu_tp))+b1/four
                  call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,&
                                          index_p4,rff_p,integ4,cases)
      write(*,*)'55==', tinf, infinity,dd,del,a4,b4,&
                                          index_p4,rff_p,integ4(1),cases
      write(unit = *, fmt = *)'************************************************************'
                  endif
                  !write(*,*)'mu=',mucos,weierstrassP(p+fzero,g2,g3,dd,del),mu_tp,b0,b1        
              endif 
              !******************************************************************** 
          else
              count_num=1        
              goto 10        
          endif        
      endif
      return
      end Function mucos

!********************************************************************************************
      Function radius(p,f1234r,lambda,q,a_spin,robs,scal,sign_pr)
!============================================================================================
!*     PURPOSE:  Computes function r(p) defined by equation (41) and (49) in Yang & Wang (2012). That is
!*               r(p)=b0/(4*\wp(p+PIr;g_2,g_3)-b1)+r_tp1. \wp(p+PIr;g_2,g_3) is the Weierstrass'
!*               elliptic function; Or r=r_+, r=r_-. 
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f1234r---------p_1, r components of four momentum of a photon measured under a LNRF. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.  
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  radius---------radial coordinate of photon corresponding to a given p.
!*               sign_pr--------the sign of r component of 4-momentum of the photon. 
!*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision, intent(in) :: p, f1234r, lambda, q, a_spin, robs, scal
      Double precision, intent(out) :: sign_pr
      Double precision :: up, periodwp, half_periodwp, p_ini_r_tp2
      Double precision radius,rhorizon ,zero,integ4(4),&
             bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tinf,tinf1,PI0,PI1, cr,dr,integ04(4),&
             u,v,w,L1,L2,thorizon,m2,pinf,sn,cn,dn,a4,b4,one,two,four,PI2,ttp,sqt3,&
             integ14(4),three,six,nine,r_tp1,r_tp2,tp2,tp,t_inf,PI0_total,&
             PI0_inf_obs,PI0_obs_hori,PI01,PI0_total_2,rff_p
      Double precision f1234r_1,lambda_1,q_1,p_1,a_spin_1,robs_1,scal_1
      parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0,six=6.D0,nine=9.D0)
      complex*16 bb(1:4),dd(3)
      integer ::  reals,i,p4,cases_int,del,index_p4(4),cases,count_num=1
      logical :: robs_eq_rtp,indrhorizon
      save :: f1234r_1,lambda_1,q_1,a_spin_1,robs_1,scal_1,r_tp1,r_tp2,reals,&
                robs_eq_rtp,indrhorizon,cases,bb,rhorizon,b0,b1,b2,b3,g2,g3,dd,del,cc,tinf,tp2,&
                thorizon,tinf1,PI0,u,w,v,L1,L2,m2,t_inf,pinf,a4,b4,PI0_total,PI0_inf_obs,PI0_obs_hori,&
                PI0_total_2, half_periodwp, periodwp, p_ini_r_tp2

 20   continue
      If(count_num.eq.1)then
          f1234r_1=f1234r 
          lambda_1=lambda
          q_1=q 
          a_spin_1=a_spin
          robs_1=robs
          scal_1=scal          
    !*********************************************************************************************          
          rhorizon=one+dsqrt(one-a_spin**2)
          a4=zero
          b4=one
          cc=a_spin**2-lambda**2-q
          robs_eq_rtp=.false.
          indrhorizon=.false.
          call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
                             reals,robs_eq_rtp,indrhorizon,cases,bb)
          !write(*,*)'hererererer = = = ', p
          If(reals.ne.0)then
! equations (35)-(38) in Yang & Wang (2012).
              b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
              b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
              b2=four/three*r_tp1
              b3=one
              g2=three/four*(b1**2-b0*b2)
              g3=one/16.D0*(3*b0*b1*b2-2*b1**3-b0**2*b3)
! equation (39) in Yang & Wang (2012).
              If(robs-r_tp1.ne.zero)then          
                  tinf=b0/four/(robs-r_tp1)+b1/four
              else
                  tinf=infinity
              endif 
              If(rhorizon-r_tp1.ne.zero)then
                  thorizon=b1/four+b0/four/(rhorizon-r_tp1)
              else
                  thorizon=infinity           
              endif
              tp2=b0/four/(r_tp2-r_tp1)+b1/four   
              tinf1=b1/four

              call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)          
              index_p4(1)=0          
              cases_int=1
! equation (42) in Yang & Wang (2012).
              call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,&
                                 index_p4,rff_p,integ04,cases_int)          
              PI0=integ04(1)
              select case(cases)
              case(1)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then  
                        call weierstrass_int_j3(tinf1,infinity,dd,del,a4,b4,&
                                            index_p4,rff_p,integ14,cases_int)        
                        PI0_total=PI0+integ14(1)
                        If(p.lt.PI0_total)then   
! equation (41) in Yang & Wang (2012).                          
                            radius=r_tp1+b0/(four*weierstrassP&
                                        (p-PI0,g2,g3,dd,del)-b1)
                            If ( p < PI0 ) then
                                sign_pr = -one
                            else
                                sign_pr = one
                            endif
                        else
                            radius=infinity  !Goto infinity, far away.
                            sign_pr = one
                        endif        
                    else
                        call weierstrass_int_J3(tinf1,tinf,dd,del,a4,b4,&
                                        index_p4,rff_p,integ04,cases_int)
                        PI0_inf_obs=integ04(1)
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif
                        sign_pr = one
                    endif
                else
                    If(f1234r.lt.zero)then
                        call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,&
                                           index_p4,rff_p,integ04,cases_int)        
                        PI0_obs_hori=integ04(1) 
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif
                        sign_pr = -one
                    else
                        call weierstrass_int_J3(tinf1,tinf,dd,del,a4,b4,&
                                         index_p4,rff_p,integ04,cases_int)
                        PI0_inf_obs=integ04(1)        
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                             radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif
                        sign_pr = one
                    endif                 
                endif
              case(2)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then
                        PI01=-PI0
                    else
                        PI01=PI0        
                    endif
! equation (41) in Yang & Wang (2012). 
                    radius=r_tp1+b0/(four*weierstrassP(p+PI01,g2,g3,dd,del)-b1)
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! this part of the code aims on obtaining the sign of the r component of  !|
               ! 4-momentum of the photon.                                               !|
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    call weierstrass_int_J3(tp2,infinity,dd,del,a4,b4,&                  !|
                                           index_p4,rff_p,integ14,cases_int)             !|
                    half_periodwp = integ14(1)                                           !|
                    periodwp = two * half_periodwp                                       !|
                    If (f1234r > zero) then                                              !|
                        up = dmod(p + PI01, periodwp)                                    !|
                        If ( up <= half_periodwp ) then                                  !|
                            sign_pr = one                                                !|
                        else                                                             !|
                            sign_pr = - one                                              !|
                        endif                                                            !|
                    Else If (f1234r < zero) then                                         !|
                        up = dmod(p + PI01 + half_periodwp, periodwp)                    !|
                        If (up < half_periodwp) then                                     !|
                            sign_pr = -one                                               !|
                        Else                                                             !|
                            sign_pr = one                                                !|
                        endif                                                            !|
                    Else                                                                 !|
                        If (robs == r_tp1) then                                          !|
                            If ( dmod(p, periodwp) <= half_periodwp ) then               !|
                                sign_pr = one                                            !|
                            else                                                         !|
                                sign_pr = -one                                           !|
                            endif                                                        !|
                        else                                                             !| 
                            If ( dmod(p + half_periodwp, periodwp) &                     !|
                                                   <= half_periodwp ) then               !|
                                sign_pr = one                                            !|
                            else                                                         !|
                                sign_pr = -one                                           !|
                            endif                                                        !| 
                        endif                                                            !|
                    Endif                                                                !|
               !============================================================================
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                else        
                    If(f1234r.le.zero)then
                        call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,&
                                           index_p4,rff_p,integ14,cases_int)        
                        PI0_obs_hori = integ14(1) 
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif
                        sign_pr = -one                     
                    else
                        call weierstrass_int_J3(tp2,thorizon,dd,del,a4,b4,index_p4,&
                                                             rff_p,integ14,cases_int)                        
                        call weierstrass_int_J3(tp2,tinf,dd,del,a4,b4,&
                                       index_p4,rff_p,integ4,cases_int)
                        p_ini_r_tp2 = integ4(1)
                        PI0_total_2=integ14(1)+integ4(1)
                        If(p.lt.PI0_total_2)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        If (p <= p_ini_r_tp2) then              !|
                            sign_pr = one                       !|
                        else                                    !|
                            sign_pr = -one                      !|
                        endif                                   !|
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
                    endif                              
                endif
              end select                    
              If(a_spin.eq.zero)then
                If(cc.eq.zero)then
                    If(f1234r.lt.zero)then
                        If(p.lt.one/rhorizon-one/robs)then
                            radius=robs/(robs*p+one)
                        else
                            radius=rhorizon                  
                        endif
                        sign_pr = -one
                    else
                        If(p.lt.one/robs)then
                            radius=robs/(one-robs*p)
                        else
                            radius=infinity          
                        endif  
                        sign_pr = one                      
                    endif        
                endif
                If(cc.eq.-27.D0)then
                    sqt3=dsqrt(three)        
                    If(f1234r.lt.zero)then
                        cr=-three*dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)&
                                           /(three-robs))*dexp(three*sqt3*p)-sqt3
                        dr=-dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                          (three-robs))*dexp(three*sqt3*p)+two/sqt3
                        If(p.ne.zero)then        
                            radius=(three+cr*dr+dsqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=robs!infinity
                        endif
                        sign_pr = - one
                    else        
                        cr=-three*dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                      (three-robs))*dexp(-three*sqt3*p)-sqt3
                        dr=-dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                      (three-robs))*dexp(-three*sqt3*p)+two/sqt3
                        PI0=dLog(dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                      (robs-three)))/three/sqt3-dLog(one+two/sqt3)/three/sqt3
                        If(p.lt.PI0)then        
                            radius=(three+cr*dr+dsqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=infinity
                        endif  
                        sign_pr = one                      
                    endif                
                endif
              endif        
          else
            u=real(bb(4))
            w=dabs(aimag(bb(4)))
            v=dabs(aimag(bb(2)))
            If(u.ne.zero)then
! equation (45) in Yang & Wang (2012). 
                L1=(four*u**2+w**2+v**2+dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                L2=(four*u**2+w**2+v**2-dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012). 
                thorizon=dsqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/&
                                  (L1-one))/dsqrt((rhorizon-u)**2+w**2)
! equation (48) in Yang & Wang (2012). 
                m2=(L1-L2)/L1
                tinf=dsqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/&
                          (L1-one))/dsqrt((robs-u)**two+w**two)
                t_inf=dsqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012). 
                pinf=EllipticF(tinf,m2)/w/dsqrt(L1)
                call sncndn(p*w*dsqrt(L1)+dsign(one,f1234r)*&
                            pinf*w*dsqrt(L1),one-m2,sn,cn,dn)
                If(f1234r.lt.zero)then
                    PI0=pinf-EllipticF(thorizon,m2)/(w*dsqrt(L1))        
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
                        radius=u+(-two*u+w*(L1-L2)*sn*dabs(cn))/((L1-L2)*sn**two-(L1-one))
                    else
                        radius=rhorizon
                    endif   
                    sign_pr = - one                 
                else
                    PI0=EllipticF(t_inf,m2)/(w*dsqrt(L1))-pinf
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
                        radius=u+(-two*u-w*(L1-L2)*sn*dabs(cn))/((L1-L2)*sn**two-(L1-one))        
                    else
                        radius=infinity
                    endif
                    sign_pr = one
                endif
            else
                If(f1234r.lt.zero)then
                    if(p.lt.(datan(robs/w)-datan(rhorizon/w))/w)then
                        radius=w*dtan(datan(robs/w)-p*w)        
                    else
                        radius=rhorizon
                    endif
                    sign_pr = - one
                else
                    if(p.lt.(PI/two-datan(robs/w))/w)then
                        radius=w*dtan(datan(robs/w)+p*w)        
                    else
                        radius=infinity
                    endif 
                    sign_pr = one               
                endif
            endif                        
          endif
          count_num=count_num+1
      else
        If(f1234r.eq.f1234r_1.and.lambda.eq.lambda_1.and.q.eq.q_1.and.&
        a_spin.eq.a_spin_1.and.robs.eq.robs_1.and.scal.eq.scal_1)then
     !***************************************************************************************************
        !write(*,*)'I am here, how are you dears. p =', p, reals, cases
        If(reals.ne.0)then        
            index_p4(1)=0        
            cases_int=1
            select case(cases)
            case(1)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then  
                        If(p.lt.PI0_total)then 
! equation (41) in Yang & Wang (2012).                                
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0, &
                                                      g2,g3,dd,del)-b1)
                            If ( p < PI0 ) then
                                sign_pr = -one
                            else
                                sign_pr = one
                            endif                                  
                        else
                            radius=infinity  !Goto infinity, far away.
                            sign_pr = one
                        endif                
                    else
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,&
                                                     g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif      
                        sign_pr = one                  
                    endif
                else
                    If(f1234r.lt.zero)then 
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,&
                                                     g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif    
                        sign_pr = -one                            
                    else
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,&
                                                     g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif  
                        sign_pr = one      
                    endif                 
                endif
            case(2)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then
                        PI01=-PI0
                    else
                        PI01=PI0
                    endif
! equation (41) in Yang & Wang (2012). 
                    radius=r_tp1+b0/(four*weierstrassP(&
                                p+PI01,g2,g3,dd,del)-b1) 
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! this part of the code aims on obtaining the sign of the r component of  !|
               ! 4-momentum of the photon.                                               !|
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                    If (f1234r > zero) then                                              !|
                        up = dmod(p + PI01, periodwp)                                    !|
                        If ( up <= half_periodwp ) then                                  !|
                            sign_pr = one                                                !|
                        else                                                             !|
                            sign_pr = - one                                              !|
                        endif                                                            !|
                    Else If (f1234r < zero) then                                         !|
                        up = dmod(p + PI01 + half_periodwp, periodwp)                    !|
                        If (up < half_periodwp) then                                     !|
                            sign_pr = -one                                               !|
                        Else                                                             !|
                            sign_pr = one                                                !|
                        endif                                                            !|
                    Else                                                                 !|
                        If (robs == r_tp1) then                                          !|
                            If ( dmod(p, periodwp) <= half_periodwp ) then               !|
                                sign_pr = one                                            !|
                            else                                                         !|
                                sign_pr = -one                                           !|
                            endif                                                        !|
                        else                                                             !| 
                            If ( dmod(p + half_periodwp, periodwp) &                     !|
                                                   <= half_periodwp ) then               !|
                                sign_pr = one                                            !|
                            else                                                         !|
                                sign_pr = -one                                           !|
                            endif                                                        !| 
                        endif                                                            !|
                    Endif                                                                !|
               !============================================================================
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                   
                else        
                    If(f1234r.le.zero)then         
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(&
                                       p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif   
                        sign_pr = -one                        
                    else
                        If(p.lt.PI0_total_2)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(&
                                         p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif 
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        If (p <= p_ini_r_tp2) then              !|
                            sign_pr = one                       !|
                        else                                    !|
                            sign_pr = -one                      !|
                        endif                                   !|
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
                    endif                              
                endif
            end select                    
            If(a_spin.eq.zero)then
                If(cc.eq.zero)then
                    If(f1234r.lt.zero)then
                        If(p.lt.one/rhorizon-one/robs)then
                            radius=robs/(robs*p+one)
                        else
                            radius=rhorizon                  
                        endif
                        sign_pr = -one
                    else
                        If(p.lt.one/robs)then
                            radius=robs/(one-robs*p)
                        else
                            radius=infinity          
                        endif 
                        sign_pr = one                       
                    endif        
                endif
                If(cc.eq.-27.D0)then
                    sqt3=dsqrt(three)        
                    If(f1234r.lt.zero)then
                        cr=-three*dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                              (three-robs))*dexp(three*sqt3*p)-sqt3
                        dr=-dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                              (three-robs))*dexp(three*sqt3*p)+two/sqt3
                        If(p.ne.zero)then        
                            radius=(three+cr*dr+dsqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=robs!infinity
                        endif
                        sign_pr = - one
                    else        
                        cr=-three*dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                            (three-robs))*dexp(-three*sqt3*p)-sqt3
                        dr=-dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                            (three-robs))*dexp(-three*sqt3*p)+two/sqt3
                        PI0=dLog(dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                              (robs-three)))/three/sqt3-dLog(one+two/sqt3)/three/sqt3
                        If(p.lt.PI0)then        
                            radius=(three+cr*dr+dsqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=infinity
                        endif 
                        sign_pr = one                       
                    endif                
                endif
            endif        
        else
            If(u.ne.zero)then
                call sncndn(p*w*dsqrt(L1)+dsign(one,f1234r)*pinf*w*dsqrt(L1),one-m2,sn,cn,dn)
                If(f1234r.lt.zero)then        
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
                        radius=u+(-two*u+w*(L1-L2)*sn*dabs(cn))/((L1-L2)*sn**two-(L1-one))
                    else
                        radius=rhorizon
                    endif
                    sign_pr = - one                   
                else
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
                        radius=u+(-two*u-w*(L1-L2)*sn*dabs(cn))/((L1-L2)*sn**two-(L1-one))        
                    else
                        radius=infinity
                    endif
                    sign_pr = one 
                endif
            else
                If(f1234r.lt.zero)then
                    if(p.lt.(datan(robs/w)-datan(rhorizon/w))/w)then
                        radius=w*dtan(datan(robs/w)-p*w)        
                    else
                        radius=rhorizon
                    endif
                    sign_pr = - one 
                else
                    if(p.lt.(PI/two-datan(robs/w))/w)then
                        radius=w*dtan(datan(robs/w)+p*w)        
                    else
                        radius=infinity
                    endif 
                    sign_pr = one                
                endif
            endif                        
        endif
      !***************************************************************************************************
        else
            count_num=1        
            goto  20
        endif
      endif        
      return                 
      End function radius

!********************************************************************************************
      Function phi(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)
!********************************************************************************************
!*     PURPOSE:  Computes function \phi(p). 
!*     INPUTS:   p--------------independent variable, which must be nonnegative.  
!*               f1234----------array of p_1, p_2, p_3, p_4, which are the components of four-
!*                              momentum of a photon measured under the LNRF frame. This array 
!*                              can be computed by subroutine lambdaq(...), see below.   
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer. 
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phi------------value of function \phi(p). 
!*     ROUTINES CALLED: INTRPART, INTTPART.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ******************************************
        Implicit none
        Double precision p,sinobs,muobs,a_spin,phi,phi_r,phi_c,twopi,mu,f1234(4),timet,&
                         two,arc,Rab,robs,scal,zero,one,lambda,q,p_axis,mu_tp,mu_tp2,phi_c1,&
                         time_r,aff_r,mu_cos,r_coord,sign_pr,sign_pth
        parameter (zero=0.D0,one=1.D0,two=2.D0)
        logical :: rotate,clines,err,mobseqmtp
        integer  tm1,tm2,reals,tr1,tr2

        twopi=two*PI
        rotate=.false.
        err=.false.
        IF(a_spin.EQ.ZERO.and.(f1234(2).NE.zero.or.muobs.NE.zero))THEN
!When spin a=0, and the motion of photon is not confined in the equatorial plane, then \phi_r = 0.
            phi_r=zero
        ELSE
            call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,&
                                scal,phi_r,time_r,aff_r,r_coord,tr1,tr2,sign_pr)
        ENDIF         
! Call integrat_theta_part to evaluate \phi_\mu.
        call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,&
                     a_spin,scal,phi_c,timet,mu_cos,tm1,tm2,sign_pth)
        !write(*,*)'ss=', p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,scal
        !write(*,*)'phi2=',phi_r,phi_c,tm1,tm2
        If(dabs(muobs).ne.one)then
! When the observer is on the equatorial plane, and p_\theta (f1234(2)) = 0, then the photon is
! confined in the equatorial plane.
            If(muobs.eq.zero.and.f1234(2).eq.zero)then
                phi_c=zero
            endif
!Equation (74) of Yang & Wang (2012).
            phi=-(phi_r+phi_c)        
            If(f1234(3).eq.zero)then
                phi=phi+(tm1+tm2)*PI
            endif        
            phi=dMod(phi,twopi)
            If(phi.lt.zero)then
                phi=phi+twopi
            Endif
        else       
!Equation (74) of Yang & Wang (2012). 
            phi=-(phi_c+phi_r+(tm1+tm2)*PI)
            Rab=dsqrt(f1234(3)**two+f1234(2)**two)
            If(phi.ne.zero)then
                rotate=.true.
            endif
            If(Rab.ne.zero)then
! a muobs was multiplied to control the rotate direction
                if((f1234(3).ge.zero).and.(f1234(2).gt.zero))then
                    phi=phi+dasin(f1234(2)/Rab)  
                endif
                if((f1234(3).lt.zero).and.(f1234(2).ge.zero))then
                    phi=phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).le.zero).and.(f1234(2).lt.zero))then
                    phi=phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).gt.zero).and.(f1234(2).le.zero))then
                    phi=phi+twopi+dasin(f1234(2)/Rab)
                endif
            else
                phi=zero
            endif
            If(rotate)then
                phi=dMod(phi,twopi)
                If(phi.lt.zero)then
                    phi=phi+twopi
                Endif
            Endif        
        endif
        return
      End Function phi

!********************************************************************************************
     SUBROUTINE GEOKERR(p_int,rp,mup,varble,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        tr1,tr2,tm1,tm2,radi,mu,time,phi,sigma,sign_pr,sign_pth) 
!********************************************************************************************
!*     PURPOSE:  Computes four Boyer-Lindquist coordinates (r,\mu,\phi,t) and affine parameter 
!*               \sigma as functions of parameter p, i.e. functions r(p), \mu(p), \phi(p), t(p)
!*               and \sigma(p). Cf. discussions in Yang & Wang (2012).    
!*     INPUTS:   p_int----------this parameter will be taken as independent variable, if 
!*                              varble='p', which must be nonnegative.
!*               rp-------------this parameter will be taken as independent variable, if 
!*                              varble='r'.
!*               mup------------this parameter will be taken as independent variable, if 
!*                              varble='mu'.
!*               varble---------Tell the routine which parameter to be as independent variable,
!*                              r, mu or p.
!*               f12342---------array of f_1, f_2, f_3, f_4, which was defined by equation (102)-(105) 
!*                              in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or initialposition of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               tm1,tm2--------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively. If varble='mu', these two parameter must be provided. 
!*               tr1,tr2--------number of times of photon meets turning points r_tp1 and r_tp2
!*                              respectively. If varble='r', these two parameter must be provided.        
!*     OUTPUTS:  radi-----------value of function r(p). 
!*               mu-------------value of function \mu(p). 
!*               phi------------value of function \phi(p).
!*               time-----------value of function t(p).
!*               sigma----------value of function \sigma(p).          
!*               tm1,tm2--------number of times of the photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively. 
!*               tr1,tr2--------number of times of the photon meets turning points r_tp1 and r_tp2
!*                              respectively.
!*     ROUTINES CALLED: INTRPART, INTTPART.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
        IMPLICIT NONE 
        DOUBLE PRECISION f1234(4),lambda,q,sinobs,muobs,a_spin,robs,scal,radi,mu,time,phi,sigma,&
                zero,one,two,three,four,phi_r,time_r,aff_r,phi_t,time_t,p,mu_tp,mu_tp2,Rab,&
                rp,mup,p_int,mu_cos,r_coord,sign_pth, sign_pr
        CHARACTER varble
        PARAMETER(zero=0.D0, one=1.D0, two=2.D0, three=3.D0, four=4.D0)
        LOGICAL rotate,err,mobseqmtp
        INTEGER tm1,tm2,tr1,tr2,reals,del
 
        SELECT CASE(varble)
        CASE('r')
            radi=rp
            p=r2p(f1234(1),rp,lambda,q,a_spin,robs,scal,tr1,tr2)
            mu=mucos(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal,sign_pth) 
        CASE('mu')
            mu=mup
            p=mu2p(f1234(3),f1234(2),lambda,q,mup,sinobs,muobs,a_spin,tm1,tm2,scal)
            radi=radius(p,f1234(1),lambda,q,a_spin,robs,scal,sign_pr)
        CASE('p')
            p=p_int 
        END SELECT 
        !************************************************************************************

! Call integrate_r_part to evaluate t_r,\phi_r,\sigma_r, and function r(p)=r_coord.
        call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,&
                              scal,phi_r,time_r,aff_r,r_coord,tr1,tr2,sign_pr)
! Call integrate_theta_part to evaluate t_\mu,\phi_\mu,\sigma_\mu, and function \mu(p)=mu_cos.
        call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal,&
                              phi_t,time_t,mu_cos,tm1,tm2,sign_pth) 
        radi=r_coord
        mu=mu_cos
!time coordinate value **************************************************************
        time=time_r+time_t
!affine parameter value *************************************************************
        sigma=aff_r+time_t  
!phi coordinate value ***************************************************************
        rotate=.false.
        err=.false.
        IF(ABS(muobs).NE.ONE)THEN
! Equation (74) of Yang & Wang (2012).
            phi=-(phi_r+phi_t)
            IF(f1234(3).EQ.zero)THEN 
                phi=phi+(tm1+tm2)*PI
            ENDIF
            phi=DMOD(phi,twopi)
            IF(phi.LT.zero)THEN
                phi=phi+twopi
            ENDIF
        ELSE 
! Equation (74) of Yang & Wang (2012).
            phi=-(phi_t+phi_r+(tm1+tm2)*PI)
            Rab=dsqrt(f1234(3)**two+f1234(2)**two)
            IF(phi.NE.zero)THEN
                rotate=.TRUE.
            ENDIF
            IF(Rab.NE.zero)THEN
! a muobs was multiplied to control the rotate direction
                if((f1234(3).ge.zero).and.(f1234(2).gt.zero))then
                    phi=muobs*phi+dasin(f1234(2)/Rab)  
                endif
                if((f1234(3).lt.zero).and.(f1234(2).ge.zero))then
                    phi=muobs*phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).le.zero).and.(f1234(2).lt.zero))then
                    phi=muobs*phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).gt.zero).and.(f1234(2).le.zero))then
                    phi=muobs*phi+twopi+dasin(f1234(2)/Rab)
                endif
            ELSE
                phi=zero
            ENDIF
            IF(rotate)THEN
                phi=dMod(phi,twopi)
                IF(phi.LT.zero)THEN
                    phi=phi+twopi
                ENDIF
            ENDIF
        ENDIF      
        RETURN
      END SUBROUTINE GEOKERR

!**************************************************
      Function rms(a_spin)                       
!**************************************************
!*     PURPOSE: Computes inner most stable circular orbit r_{ms}. 
!*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
!*     OUTPUTS:  radius of inner most stable circular orbit: r_{ms}
!*     ROUTINES CALLED: root4
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ******************************************************* 
      implicit none
      Double precision rms,a_spin,b,c,d,e
      complex*16 rt(1:4)
      integer  reals,i
      If(a_spin.eq.0.D0)then
            rms=6.D0
            return
      endif
      b=0.D0
      c=-6.D0
      d=8.D0*a_spin
      e=-3.D0*a_spin**2
        ! Bardeen et al. (1972) 
      call root4(b,c,d,e,rt(1),rt(2),rt(3),rt(4),reals)
      Do i=4,1,-1
         If(aimag(rt(i)).eq.0.D0)then
            rms=real(rt(i))**2
             return                     
         endif             
      enddo
      end function rms

!*********************************************************** 
      Function rph(a_spin)
!***********************************************************
!*     PURPOSE: Computes photon orbit of circluar orbits: r_{ph}. 
!*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
!*     OUTPUTS:  radius of photon orbit: r_{ph}
!*     ROUTINES CALLED: NONE
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision rph,a_spin
        ! Bardeen et al. (1972) 
      rph=2.D0*(1.D0+dcos(2.D0/3.D0*dacos(-a_spin)))
      End function  rph      

!************************************************************* 
      Function rmb(a_spin)
!*************************************************************
!*     PURPOSE: Computes marginally bound orbit of circluar orbits: r_{mb}. 
!*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
!*     OUTPUTS:  radius of marginally bound orbit: r_{mb}
!*     ROUTINES CALLED: NONE
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision rmb,a_spin
        ! Bardeen et al. (1972)  
      rmb=2.D0-a_spin+2.D0*dsqrt(1.D0-a_spin)
      End function  rmb      

!********************************************************************************************
      subroutine mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp1,mu_tp2,reals,mobseqmtp)
!********************************************************************************************
!*     PURPOSE: Returns the coordinates of turning points \mu_tp1 and \mu_tp2 of poloidal motion, judges
!*                whether the initial poloidal angle \theta_{obs} is one of turning points, if 
!*                it is true then mobseqmtp=.TRUE..  
!*     INPUTS:   f12342---------p_2, the \theta component of four momentum of the photon measured 
!*                              under the LNRF, see equation (84) in Yang & Wang (2012).
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*     OUTPUTS:  mu_tp1, mu_tp2----the turning points, between which the poloidal motion of 
!*                                 the photon was confined, and mu_tp2 <= mu_tp1. 
!*               reals------number of real roots of equation \Theta_\mu(\mu)=0.
!*               mobseqmtp---If mobseqmtp=.TRUE., then muobs equals to be one of the turning points. 
!*     ROUTINES CALLED: NONE
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision f1234(4),sinobs,muobs,a_spin,lambda,q,zero,one,two,four,&
                  mu_tp1,mu_tp2,delta,mutemp,Bprime,f12342,f12343 
      integer  reals
      logical :: mobseqmtp
      parameter (zero=0.D0,two=2.0D0,four=4.D0,one=1.D0)
 
      mobseqmtp=.false.
        If(a_spin .eq. zero)then 
          If(f12342.ne.zero)then
              mu_tp1=dsqrt(q/(lambda**two+q))
               mu_tp2=-mu_tp1      
          else
               mu_tp1=dabs(muobs)
               mu_tp2=-mu_tp1
               mobseqmtp=.true.
          endif
          reals=2    
        ELSE
          If(lambda.ne.zero)then
            delta=(a_spin**2-lambda**2-q)**2+four*a_spin**2*q
            mu_tp1=dsqrt( dabs((dsqrt(delta)-(lambda**2+q-a_spin**2))/two) )/dabs(a_spin) 
            !write(*,*)mu_tp1, muobs, f12342
            If( dsqrt(delta)+(lambda**two+q-a_spin**two).le.zero )then
                mu_tp2=dsqrt(-(dsqrt(delta)+(lambda**two+q-a_spin**two))/two)/dabs(a_spin)
                If(f12342.eq.zero)then
                  If(dabs(muobs-mu_tp1).le.1.D-4)then
                      mu_tp1=dabs(muobs)
                  else       
                      mu_tp2=dabs(muobs)
                  endif
                  mobseqmtp=.true.
                endif      
                reals=4
            else
                If(f12342.ne.zero)then      
                  mu_tp2=-mu_tp1
                else
                  mu_tp1=dabs(muobs)
                  mu_tp2=-mu_tp1
                  mobseqmtp=.true.
                endif      
                reals=2
            endif      
          else 
            If(dabs(muobs).ne.one)then
                If(q.le.zero)then
                    If(f12342.ne.zero)then
                        mu_tp2=dsqrt(-q)/dabs(a_spin)
                    else
                        mu_tp2=dabs(muobs)!a=B=zero.
                        mobseqmtp=.true.
                    endif
                    reals=4
                else
                    mu_tp2=-one
                    reals=2
                endif 
                mu_tp1=one
            else
                mu_tp1=one
                If(q.le.zero.and.f12342*f12342+f12343*f12343.ne.zero)then
                    mu_tp2=dsqrt(-q)/dabs(a_spin)
                    reals=4
                else
                    mu_tp2=-one
                    reals=2
                endif
            endif
          endif
        ENDIF
      If(dabs(muobs).eq.one)mobseqmtp=.true.
      If(muobs.lt.zero.and.reals.eq.4)then
            mutemp=mu_tp1
            mu_tp1=-mu_tp2
            mu_tp2=-mutemp
      endif
      return
      end subroutine mutp      

!============================================================================================
      Subroutine radiustp(f12341,a_spin,robs,lambda,q,r_tp1,&
                        r_tp2,reals,robs_eq_rtp,indrhorizon,cases,bb)
!********************************************************************************************
!*     PURPOSE: Returns the coordinates of turning points r_tp1 and r_tp2 of radial motion, judges
!*                whether the initial radius robs is one of turning points, if 
!*                it is true then robs_eq_rtp=.TRUE.. And if r_tp1 less or equal r_horizon,
!*                then indrhorizon=.TRUE. Where r_horizon is the radius of the event horizon.  
!*     INPUTS:   f12341---------p_r, the r component of four momentum of the photon measured 
!*                              under the LNRF, see equation (83) in Yang & Wang (2012).
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               robs-----------radial coordinate of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*     OUTPUTS:  r_tp1, r_tp2----the turning points, between which the radial motion of 
!*                                 the photon was confined, and r_tp2 >= r_tp1.
!*               bb(1:4)----roots of equation R(r)=0.                
!*               reals------number of real roots of equation R(r)=0.
!*               robs_eq_rtp---If robs_eq_rtp=.TRUE., then robs equal to be one of turning points. 
!*               cases-------If r_tp2=infinity, then cases=1, else cases=2.
!*               indrhorizon----if r_tp1 less or equals r_horizon, indrhorizon=.TRUE.. 
!*     ROUTINES CALLED: root4
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
       Double precision f12341,a_spin,robs,lambda,q,r_tp1,r_tp2,&
                  zero,one,two,four,b1,c1,d1,e1,r1(2),rhorizon
      integer  reals,i,j,cases
      logical :: robs_eq_rtp,indrhorizon
      complex*16 bb(4)
      parameter (zero=0.D0,two=2.0D0,four=4.D0,one=1.D0)       

      rhorizon=one+dsqrt(one-a_spin**two)
      robs_eq_rtp=.false.
      indrhorizon=.false.
      b1=zero
      c1=a_spin**two-lambda**two-q
      d1=two*(q+(lambda-a_spin)**two)
      e1=-q*a_spin**two
      call root4(b1,c1,d1,e1,bb(1),bb(2),bb(3),bb(4),reals) 
      !write(*,*)'sss=', b1,c1,d1,e1,bb(1),bb(2),bb(3),bb(4),reals
          SELECT CASE(reals) 
          CASE(4)   
              IF(f12341.eq.zero)THEN
                  IF(dabs(robs-real(bb(4))) .LE. 1.D-4)THEN
                      r_tp1=robs
                      r_tp2=infinity
                      cases=1
                  ENDIF 
                  IF(dabs(robs-real(bb(2))) .LE. 1.D-4)THEN
                      r_tp1=robs  
                      r_tp2=real(bb(3))
                      cases=2
                  ENDIF
                  IF(dabs(robs-real(bb(3))) .LE. 1.D-4)THEN
                      r_tp1=real(bb(2))       
                      r_tp2=robs
                      cases=2
                  ENDIF
                  IF(dabs(robs-real(bb(1))) .LE. 1.D-4)THEN
                      r_tp1=-infinity
                      r_tp2=robs 
                      cases=3
                      write(*,*)'radiustp(): wrong! 4 roots, cases = 3'
                      stop  
                  ENDIF  
                  robs_eq_rtp = .TRUE. 
              ELSE 
                  If( robs .gt. real(bb(4)) )then      
                      r_tp1=real(bb(4))      
                      r_tp2=infinity
                      cases=1 
                  else if ( robs > real(bb(3)) ) then  
                      r_tp1=robs!real(bb(4))      
                      r_tp2=infinity
                      cases=1  
                      robs_eq_rtp = .TRUE. 
                      write(*,*)'p_r_LNRF === ', f12341, robs, bb
                      !stop
                  else
                      If( (robs.ge.real(bb(2)) .and. robs.le.real(bb(3))) .or. &
                      dabs( robs - real(bb(3)) ) <=1.D-10 )then      
                          r_tp1=real(bb(2))      
                          r_tp2=real(bb(3))
                          if( r_tp2 < robs )r_tp2 = robs
                          cases=2
                      else
                          IF( real(bb(1)) .GT. rhorizon .AND.  robs .LE. real(bb(1))  )THEN
                              write(*,*)'radiustp(): wrong! 4 roots,cases = 3'
                              stop 
                              r_tp2=real(bb(1))  
                              r_tp1=-infinity 
                          ELSE   
                              write(*,*)'radiustp(): wrong! 4 roots',robs,bb
                              stop
                          ENDIF
                      endif
                  endif
              ENDIF          
          CASE(2)
              j=1
              Do  i=1,4
                  If (aimag(bb(i)).eq.zero) then
                      r1(j)=real(bb(i))
                      j=j+1      
                  endif
              Enddo
              IF(f12341.eq.zero)THEN 
                  IF(dabs(robs-r1(2)) .LE. 1.D-4)THEN
                      r_tp1=robs 
                      r_tp2=infinity
                      cases = 1
                  ENDIF
                  IF(dabs(robs-r1(1)) .LE. 1.D-4)THEN
                      r_tp1=-infinity 
                      r_tp2=robs 
                      write(*,*)'radiustp(): wrong! 2 roots, cases = 3'
                      stop
                  ENDIF
                  robs_eq_rtp=.TRUE.
              ELSE IF(dabs(robs-r1(2)) .LE. 1.D-10)THEN
                  r_tp1=robs 
                  r_tp2=infinity
                  cases = 1
              ELSE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  If( robs.ge.r1(2) )then      
                      r_tp1=r1(2)      
                      r_tp2=infinity
                      cases=1
                  else  
                      If( r1(1).ge.rhorizon .and. robs.le.r1(1) )then
                          write(*,*)'radiustp(): wrong! 2 roots, cases = 3'
                          stop      
                      endif
                  endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ENDIF
          CASE(0)
              r_tp1=zero
              r_tp2=infinity
              cases=1       
          END SELECT  
 
      IF(rhorizon.ge.r_tp1 .and. rhorizon.le.r_tp2)then
          indrhorizon=.true.
      Endif
      End Subroutine radiustp

!********************************************************************************************  
      Function mu2p(f12343,f12342,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,scal)
!********************************************************************************************
!*     PURPOSE:  Computes the value of parameter p from \mu coordinate. In other words, to compute 
!*               the \mu part of integral of equation (24), using formula (54) in Yang & Wang (2012).
!*               (54) is: p=-sign(p_\theta)*p_0+2*t1*p_2+2*t2*p_2. where p_\theta is initial \theta
!*               component of 4 momentum of photon.
!*     INPUTS:   f12342---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f12343---------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               t1,t2----------Number of photon meets the turning points \mu_tp1 and \mu_tp2
!*                              respectively.
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.     
!*               mu-------------\mu coordinate of photon.     
!*     OUTPUTS:  value of \mu part of integral of (24). 
!*     ROUTINES CALLED: mu2p_schwartz, mutp, root3, weierstrass_int_J3
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision mu2p,f12342,f12343,mu,sinobs,muobs,a_spin,lambda,q,mu_tp,tposition,tp2,four,&
                   bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tinf,p1,p2,pp,a4,b4,delta,two,mu_tp2,&
                   scal,zero,mutemp,one,integ4(4),three,rff_p
      parameter (zero=0.D0,two=2.D0,four=4.D0,one=1.D0,three=3.D0)
      integer  t1,t2,reals,p4,index_p4(4),del,cases
      complex*16 bb(1:4),dd(3)
      logical :: mobseqmtp  

      If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
          mu2p=zero!-one
          return            
      endif
      If(a_spin.eq.zero)then
            call mu2p_schwartz(f12343,f12342,lambda,q,mu,sinobs,muobs,t1,t2,mu2p,scal)
            return      
      endif
      
      a4=zero
      b4=one
      p4=0
      mobseqmtp=.false.
      call mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp,mu_tp2,reals,mobseqmtp)
! equatorial plane motion.
      If(mu_tp.eq.zero)then
            mu2p=zero
            return
      endif
! equations (26)-(29) in Yang & Wang (2012).
      b0=-four*a_spin**2*mu_tp**3+two*mu_tp*(a_spin**2-lambda**2-q)
      b1=-two*a_spin**2*mu_tp**2+one/three*(a_spin**2-lambda**2-q)
      b2=-four/three*a_spin**2*mu_tp
      b3=-a_spin**2
      g2=three/four*(b1**2-b0*b2)
      g3=one/16.D0*(three*b0*b1*b2-two*b1**3-b0**2*b3)
! equation (30) in Yang & Wang (2012).
      If(dabs(mu-mu_tp).ne.zero)then
           tposition=b0/(four*(mu-mu_tp))+b1/four
      else
           tposition=infinity      
      endif
      If(muobs.ne.mu_tp)then      
           tinf=b0/four/(muobs-mu_tp)+b1/four
      else
           tinf=infinity
      endif      

      call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
      index_p4(1)=0
      cases=1 

        If(mu.gt.mu_tp.or.mu.lt.mu_tp2)then
              mu2p=-one
              return
        endif      
! equation (30) in Yang & Wang (2012).
              tp2=b0/four/(mu_tp2-mu_tp)+b1/four
            If(t1.eq.0)then
               p1=zero
            else
! equation (53) in Yang & Wang (2012).
               call weierstrass_int_J3(tposition,infinity,dd,del,a4,b4,index_p4,rff_p,integ4,cases)
               p1=integ4(1)
            endif
            If(t2.eq.0)then
               p2=zero
            else
               call weierstrass_int_J3(tp2,tposition,dd,del,a4,b4,index_p4,rff_p,integ4,cases)      
               p2=integ4(1)
            endif
               call weierstrass_int_J3(tinf,tposition,dd,del,a4,b4,index_p4,rff_p,integ4,cases)            
               pp=integ4(1) 

! equation (54) in Yang & Wang (2012).
      If(mobseqmtp)then
          If(muobs.eq.mu_tp)then  
              mu2p=-pp+two*(t1*p1+t2*p2)            
          else
              mu2p=pp+two*(t1*p1+t2*p2)            
          endif       
      else
          If(f12342.lt.zero)then
            mu2p=pp+two*(t1*p1+t2*p2)
          endif
          If(f12342.gt.zero)then            
            mu2p=-pp+two*(t1*p1+t2*p2)
          endif      
      endif
      return
      end Function mu2p

!============================================================================================
      subroutine mu2p_schwartz(f12343,f12342,lambda,q,mu,sinobs,muobs,t1,t2,mu2p,scal)
!********************************************************************************************
!*     PURPOSE:  Computes the value of parameter p from \mu coordinate. In other words, to compute 
!*               the \mu part of integral of equation (24), using formula (54) in Yang & Wang (2012).
!*               (54) is: p=-sign(p_\theta)*p_0+2*t1*p_2+2*t2*p_2. where p_\theta is initial \theta
!*               component of 4 momentum of photon.
!*               And black hole spin is zero. 
!*     INPUTS:   f12342---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f12343---------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               t1,t2----------Number of photon meets the turning points \mu_tp1 and \mu_tp2
!*                              respectively.
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.     
!*               mu-------------\mu coordinate of photon.     
!*     OUTPUTS:  mu2p-----------value of \mu part of integral of (24). 
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision f12343,f12342,mu,sinobs,muobs,mu2p,pp,p1,p2,AA,BB,two,DD,&
                   lambda,q,a_spin,scal,zero,one,mu_tp,mu_tp2
      integer :: t1,t2      
      parameter(two=2.D0,zero=0.D0,one=1.D0)
      logical :: mobseqmtp
      If(f12343.eq.zero.and.f12342.eq.zero)then !this Theta=q(1-mu^2),so if B=0,then q=0.
          mu2p=-two            !so Theta_mu=0 for ever.But we do not need to 
          return  !consider it,for q=0,so the next part means that it will return
      endif !zero value.
      mobseqmtp=.false.
      If(q.gt.zero)then      
          BB=dsqrt(q)
          If(f12342.ne.zero)then
              mu_tp=dsqrt(q/(lambda**two+q))
              mu_tp2=-mu_tp      
          else
              mu_tp=muobs
              mu_tp2=-mu_tp
              mobseqmtp=.true.
          endif 
          If(dabs(muobs).eq.one)mobseqmtp=.true.            
          pp=(dasin(mu/mu_tp)-dasin(muobs/mu_tp))*mu_tp/BB      
          If(t1.eq.0)then
              p1=zero
          else      
              p1=(PI/two-dasin(mu/mu_tp))*mu_tp/BB      
          endif
          If(t2.eq.0)then
              p2=zero
          else      
              p2=(dasin(mu/mu_tp)+PI/two)*mu_tp/BB
            endif
          If(mobseqmtp)then
              If(muobs.eq.mu_tp)then  
                  mu2p=-pp+two*(t1*p1+t2*p2)            
              else
                  mu2p=pp+two*(t1*p1+t2*p2)            
              endif             
          else
              mu2p=dsign(one,-f12342)*pp+two*(t1*p1+t2*p2)       
          endif
      else      
          mu2p=zero
      endif            
      return
      end subroutine mu2p_schwartz

!********************************************************************************************
      Function r2p(f1234r,rend,lambda,q,a_spin,robs,scal,t1,t2)
!============================================================================================
!*     PURPOSE:  Computes the value of parameter p from radial coordinate. In other words, to compute 
!*               the r part of integral of equation (24), using formula (58) in Yang & Wang (2012).
!*               (58) is: p=-sign(p_r)*p_0+2*t1*p_2+2*t2*p_2. where p_r is initial radial
!*               component of 4 momentum of photon. 
!*     INPUTS:   f1234r---------f_1, which was defined by equation (106) in Yang & Wang (2012). 
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               robs-----------radial coordinate of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
!*               scal-----------a dimentionless parameter to control the size of the images. 
!*                              Which is usually be set to 1.D0.   
!*               t1,t2----------Number of photon meets the turning points r_tp1 and r_tp2
!*                              respectively in radial motion.
!*     OUTPUTS:  r2p------------value of r part of integral (24) in Yang & Wang (2012).
!*     ROUTINES CALLED: radiustp, root3, weierstrass_int_j3, EllipticF.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision r2p,p,a_spin,rhorizon,q,lambda,scal,zero,integ4(4),&
                    bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tinf,tinf1,PI0,PI1,robs,cr,dr,integ04(4),&
                    u,v,w,L1,L2,thorizon,m2,pinf,sn,cn,dn,a4,b4,one,two,four,PI2,ttp,sqrt3,&
                    integ14(4),three,six,nine,r_tp1,r_tp2,f1234r,tp2,tp,t_inf,PI0_total,&
                    PI0_inf_obs,PI0_obs_hori,PI01,PI0_total_2,pp,p1,p2,rend,rff_p
      Double precision f1234r_1,lambda_1,q_1,p_1,a_spin_1,robs_1,scal_1
      parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0,six=6.D0,nine=9.D0)
      complex*16 bb(1:4),dd(3)
      integer  reals,i,p4,cases_int,del,index_p4(4),cases,t1,t2
      logical :: robs_eq_rtp,indrhorizon
      
      rhorizon=one+dsqrt(one-a_spin**2)
      a4=zero
      b4=one
      cc=a_spin**2-lambda**2-q
      robs_eq_rtp=.false.
      indrhorizon=.false.
      call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,&
                    r_tp2,reals,robs_eq_rtp,indrhorizon,cases,bb)
      If(reals.ne.0)then
          If(rend.lt.r_tp1.or.rend.gt.r_tp2)then
              r2p=-one
              return
          endif
! equations (35)-(38) in Yang & Wang (2012).
          b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
          b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
          b2=four/three*r_tp1
          b3=one
          g2=three/four*(b1**2-b0*b2)
          g3=one/16.D0*(3*b0*b1*b2-2*b1**3-b0**2*b3)
! equation (39) in Yang & Wang (2012).
          If(robs-r_tp1.ne.zero)then      
              tinf=b0/four/(robs-r_tp1)+b1/four
          else
              tinf=infinity
          endif 
          If(rhorizon-r_tp1.ne.zero)then
              thorizon=b1/four+b0/four/(rhorizon-r_tp1)
          else
              thorizon=infinity       
          endif
          If(rend-r_tp1.ne.zero)then
              tp=b1/four+b0/four/(rend-r_tp1)
          else
              tp=infinity       
          endif
          tp2=b0/four/(r_tp2-r_tp1)+b1/four   
          tinf1=b1/four

          call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)      
          index_p4(1)=0      
          cases_int=1
! equation (42) in Yang & Wang (2012).
          call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,&
                             index_p4,rff_p,integ04,cases_int)      
          PI0=integ04(1)
          select case(cases)
          case(1)
              If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then  
                        call weierstrass_int_j3(tinf,tp,dd,del,a4,b4,&
                                     index_p4,rff_p,integ14,cases_int)      
                        pp=integ14(1)
                        If(t1.ne.zero)then     
! equation (57) in Yang & Wang (2012).                   
                            call weierstrass_int_j3(tp,infinity,dd,del,&
                                 a4,b4,index_p4,rff_p,integ14,cases_int)
                            p1=integ14(1)              
                        else
                            p1=zero  !Goto infinity, far away.
                        endif 
! equation (58) in Yang & Wang (2012).
                        r2p=pp+two*p1*t1
                    else
                        call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,&
                                     index_p4,rff_p,integ04,cases_int)
                        pp=integ04(1)
                        r2p=-pp        
                    endif
                else
                    If(f1234r.lt.zero)then
                        If(rend.le.rhorizon)then
                            tp=thorizon
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,&
                                          index_p4,rff_p,integ04,cases_int)        
                            r2p=integ04(1)
                        else
                            call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,&
                                          index_p4,rff_p,integ04,cases_int)        
                            r2p=integ04(1)        
                        endif
                    else
                        If(rend.lt.infinity)then
                            call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,&
                                         index_p4,rff_p,integ04,cases_int)        
                            r2p=-pp
                        else
                            call weierstrass_int_J3(tinf,tinf1,dd,del,a4,b4,&
                                            index_p4,rff_p,integ04,cases_int)        
                            r2p=-pp        
                        endif
                    endif                 
                endif
            case(2)
                If(.not.indrhorizon)then   
! equation (57) in Yang & Wang (2012).                 
                        call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,&
                                      index_p4,rff_p,integ4,cases_int)
                        pp=integ4(1)
                        If(t1.eq.zero)then
                            p1=zero
                        else
                            call weierstrass_int_J3(tp,infinity,dd,del,a4,b4,&
                                              index_p4,rff_p,integ4,cases_int)
                            p1=integ4(1)        
                        endif        
                        If(t2.eq.zero)then
                            p2=zero
                        else
                            call weierstrass_int_J3(tp2,tp,dd,del,a4,b4,&
                                         index_p4,rff_p,integ4,cases_int)
                            p2=integ4(1)        
                        endif
                        If(f1234r.ne.zero)then
                            r2p=dsign(one,-f1234r)*pp+two*(t1*p1+t2*p2)
                        else
! equation (58) in Yang & Wang (2012).
                            If(robs.eq.r_tp1)then
                                r2p=-pp+two*(t1*p1+t2*p2)
                            else
                                r2p=pp+two*(t1*p1+t2*p2)
                            endif        
                        endif                    
                else        
                    If(f1234r.le.zero)then
                        If(rend.le.rhorizon)then
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,&
                                          b4,index_p4,rff_p,integ4,cases_int)
                            pp=integ4(1)                            
                        else
                            call weierstrass_int_J3(tinf,tp,dd,del,a4,&
                                    b4,index_p4,rff_p,integ4,cases_int)
                            pp=integ4(1)
                        endif
                    else
                        call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,&
                                      index_p4,rff_p,integ4,cases_int)
                        pp=integ4(1)
                        If(t2.eq.zero)then
                            p2=zero
                        else
                            call weierstrass_int_J3(tp2,tp,dd,del,a4,b4,&
                                         index_p4,rff_p,integ4,cases_int)
                            p2=integ4(1)
                        endif
! equation (58) in Yang & Wang (2012).
                        r2p=-pp+two*t2*p2
                    endif                              
                endif
            end select                    
            If(a_spin.eq.zero)then
                If(cc.eq.zero)then
                    If(f1234r.lt.zero)then
                        If(rend.le.rhorizon)then
                            r2p=one/rhorizon-one/robs
                        else
                            r2p=one/rend-one/robs         
                        endif
                    else
                        If(rend.lt.infinity)then
                            r2p=one/robs-one/rend
                        else
                            r2p=one/robs
                        endif 
                    endif 
                endif
                If(cc.eq.-27.D0)then        
                    sqrt3=dsqrt(three)        
                    If(f1234r.lt.zero)then
                        If(rend.gt.rhorizon)then
                            r2p=-dlog(dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/(sqrt3))/&
                                                         (robs-three)))/(three*sqrt3)+&
                                 dlog(dabs((dsqrt(rend*(rend+6.D0))+(three+two*rend)/&
                                                         (sqrt3))/(rend-three)))/(three*sqrt3)
                        else
                            r2p=-dlog(dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/&
                                            (sqrt3))/(robs-three)))/(three*sqrt3)+&
                                 dlog(dabs((dsqrt(rhorizon*(rhorizon+6.D0))+(three+two*rhorizon)&
                                            /(sqrt3))/(rhorizon-three)))/(three*sqrt3)
                        endif
                    else        
                        If(rend.lt.infinity)then
                            r2p=-dlog(dabs((dsqrt(rend*(rend+6.D0))+(three+two*rend)/(sqrt3))/&
                                (rend-three)))/(three*sqrt3)+dlog(dabs((dsqrt(robs*(robs+6.D0))+&
                                (three+two*robs)/(sqrt3))/(robs-three)))/(three*sqrt3)
                        else
                            r2p=-dlog(one+two/sqrt3)/three/sqrt3+&
                                 dlog(dabs((dsqrt(rend*(rend+6.D0))+(three+two*rend)/&
                                    (sqrt3))/(rend-three)))/(three*sqrt3)
                        endif                                                
                    endif                
                endif
            endif        
        else
! equation (44) in Yang & Wang (2012).
            u=real(bb(4))
            w=dabs(aimag(bb(4)))
            v=dabs(aimag(bb(2)))
            If(u.ne.zero)then
! equation (45) in Yang & Wang (2012).
                L1=(four*u**2+w**2+v**2+dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                L2=(four*u**2+w**2+v**2-dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012).
                thorizon=dsqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/&
                               (L1-one))/dsqrt((rhorizon-u)**2+w**2)
                tp=dsqrt((L1-one)/(L1-L2))*(rend-u*(L1+one)/&
                         (L1-one))/dsqrt((rend-u)**2+w**2)
! equation (48) in Yang & Wang (2012).
                m2=(L1-L2)/L1
                tinf=dsqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/&
                           (L1-one))/dsqrt((robs-u)**two+w**two)
                t_inf=dsqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012).
                pinf=EllipticF(tinf,m2)/w/dsqrt(L1) 
                If(f1234r.lt.zero)then
                    If(rend.le.rhorizon)then
                        r2p=pinf-EllipticF(thorizon,m2)/(w*dsqrt(L1))
                    else
                        r2p=pinf-EllipticF(tp,m2)/(w*dsqrt(L1))
                    endif                            
                else
                    If(rend.lt.infinity)then
                        r2p=EllipticF(tp,m2)/(w*dsqrt(L1))-pinf
                    else
                        r2p=EllipticF(t_inf,m2)/(w*dsqrt(L1))-pinf
                    endif
                endif
            else
                If(f1234r.lt.zero)then
                    If(rend.le.rhorizon)then
                        r2p=(datan(robs/w)-datan(rhorizon/w))/w
                    else
                        r2p=(datan(robs/w)-datan(rend/w))/w
                    endif        
                else
                    if(rend.lt.infinity)then
                        r2p=(datan(rend/w)-datan(robs/w))/w
                    else
                        r2p=(PI/two-datan(robs/w))/w
                    endif                
                endif
            endif                        
      endif        
      return                 
      End function r2p 

!********************************************************************************************
      SUBROUTINE INTTPART(p,f12343,f12342,lambda,q,sinobs,muobs,&
                          a_spin,scal,phyt,timet,mucos,t1,t2,sign_pth)    
!********************************************************************************************
!*     PURPOSE:  Computes \mu part of integrals in coordinates \phi, t and affine parameter \sigma,
!*               expressed by equation (71) and (72) in Yang & Wang (2012).    
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f12342---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f12343---------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phyt-----------value of integral \phi_\theta expressed by equation (72) in
!*                              Yang & Wang (2012).  
!*               timet----------value of integral \t_\theta expressed by equation (71) in
!*                              Yang & Wang (2012). And \sigma_\theta=time_\theta.
!*               mucos----------value of function \mu(p).
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.
!*               sign_pth-------sign of \theta component of 4-momentum of the photon.         
!*     ROUTINES CALLED: mutp, root3, weierstrass_int_J3 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
      USE constants
      IMPLICIT NONE
      Double precision :: phyt,timet,f12343,f12342,p,sinobs,muobs,a_spin,lambda,q,&
             mu_tp1,tposition,tp2,mu,tmu,p1J2,&
             bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tinf,p1,p2,pp,Wmup,Wmum,tplus,tminus,p1J1,&
             p1I0,a4,b4 ,delta,mu_tp2,scal,mutemp ,integ4(4),integ(4),rff_p,&
             integ14(4),pp2,f1234(4),PI0,integ04(4),fzero,mu2p,PI01,h,p1_t,p2_t,pp_t,p1_phi,&
             p2_phi,pp_phi,radius,mtp1,mtp2,mucos,sqt3,difference,p_mt1_mt2,&
             PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_p,PI2_p,p1_temp,p2_temp,sign_pth,&
             half_period_wp, period_wp, u
      Double precision :: f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1 
      integer :: t1,t2,i,j,reals,cases,p4,index_p4(4),del,cases_int,&
                 count_num=1,N_temp, tt1,tt2
      complex*16 bb(1:4),dd(3)
      logical :: err,mobseqmtp
      save  f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,&
             scal_1,a4,b4,mu_tp1,mu_tp2,reals,mobseqmtp,b0,b1,b2,b3,&
             g2,g3,dd,del,PI0,Wmup,Wmum,tplus,tminus,tp2,tinf,h,p_mt1_mt2,&
             PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_p,PI2_p,PI01, &
             half_period_wp, period_wp

30      continue
        IF(count_num.eq.1)then        
            f12343_1=f12343
            f12342_1=f12342
            lambda_1=lambda
            q_1=q
            muobs_1=muobs
            sinobs_1=sinobs
            a_spin_1=a_spin 
            scal_1=scal
            t1=0
            t2=0
            !************************************************************************
            If(f12343.eq.zero.and.f12342.eq.zero.and.dabs(muobs).eq.one)then
                sign_pth = zero
                mucos = dsign(one,muobs)
                timet=zero             !this is because that mu==1 for ever
                phyt=zero              !this is because that mu==1 for ever,this 
                count_num=count_num+1
                return        !because that Theta_mu=-a^2(1-mu^2), so,mu must =+1 or -1 for ever.
            endif
            If(muobs.eq.zero.and.(dabs(lambda).lt.dabs(a_spin)).and.q.eq.zero)then
                sign_pth = zero
                timet=zero
                phyt=zero
                mucos=zero 
                count_num=count_num+1
                return                        
            endif
            mobseqmtp=.false.
            call mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,&
                                    mu_tp1,mu_tp2,reals,mobseqmtp)        
            If(mu_tp1.eq.zero)then
            !photons are confined in the equatorial plane, so the integrations about \theta are valished.
                sign_pth = zero
                timet=zero
                phyt=zero
                mucos=zero
                count_num=count_num+1
                return
            endif
            !**************************************************************************
            If(a_spin.eq.zero)then
                timet=zero
                CALL phyt_schwatz(p,f12343,f12342,lambda,q,&
                         sinobs,muobs,scal,phyt,mucos,t1,t2,sign_pth)
                count_num=count_num+1
                return
            endif
            a4=zero
            b4=one 
! equations (26)-(29) in Yang & Wang (2012). 
            b0=-four*a_spin**2*mu_tp1**3+two*mu_tp1*(a_spin**2-lambda**2-q)
            b1=-two*a_spin**2*mu_tp1**2+one/three*(a_spin**2-lambda**2-q)
            b2=-four/three*a_spin**2*mu_tp1
            b3=-a_spin**2
            g2=three/four*(b1**2-b0*b2)
            g3=one/16.D0*(three*b0*b1*b2-two*b1**3-b0**2*b3)
            call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
! equation (30) in Yang & Wang (2012). 
            If(muobs.ne.mu_tp1)then        
                tinf=b0/four/(muobs-mu_tp1)+b1/four
            else
                tinf=infinity
            endif
            If(mu_tp1-one.ne.zero)then
! equation (64) in Yang & Wang (2012). 
                Wmum=b0/(eight*(-one-mu_tp1)**2)
                Wmup=b0/(eight*(one-mu_tp1)**2) 
                tminus=b0/four/(-one-mu_tp1)+b1/four
                tplus=b0/four/(one-mu_tp1)+b1/four
            endif
            index_p4(1)=0
            cases_int=1
            call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,&
                               index_p4,rff_p,integ04,cases_int)
            PI0=integ04(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            tp2=b0/four/(mu_tp2-mu_tp1)+b1/four  
            call weierstrass_int_J3(tp2,infinity,dd,del,a4,&
                            b4,index_p4,rff_p,integ14,cases_int)
            half_period_wp = integ14(1) 
            period_wp = two * half_period_wp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! equation (34) in Yang & Wang (2012).    
            If(f12342.lt.zero)then
                PI01=-PI0        
            else
                PI01=PI0
            endif
            tmu=weierstrassP(p+PI01,g2,g3,dd,del)
! equation (32) in Yang & Wang (2012). 
            mucos = mu_tp1+b0/(four*tmu-b1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Now we get the value of parameter sign_pth
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            u = p+PI01
            if ( u <= zero ) then
                sign_pth = -one
            else
                u = dmod(u, period_wp)
                if ( u <= half_period_wp ) then
                    sign_pth = one
                else
                    sign_pth = - one
                endif
            endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            h=-b1/four        
            !to get number of turning points of t1 and t2.
            !111111111*********************************************************** 
                call weierstrass_int_J3(tinf,tmu,dd,del,a4,b4,&
                                index_p4,rff_p,integ4,cases_int) 
            !write(*,*)tp2,tinf,tmu,mu_tp1,mu_tp2
! equation (51) in Yang & Wang (2012).        
                p_mt1_mt2=integ14(1)
                PI1_p = PI0
                PI2_p=p_mt1_mt2-PI0
                pp=integ4(1)
                p1=PI0-pp
                p2=p_mt1_mt2-p1
                PI1_phi=zero
                PI2_phi=zero
                PI1_time=zero
                PI2_time=zero 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Now we determine parameters: t1 and t2 which are the times that the photon 
! meets the two turnning points, mu_tp1 and mu_tp2 respectively.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                !Do j=0,10
                !    Do i=j,j+1 
                !        If(mobseqmtp)then
                !            If(muobs.eq.mu_tp1)then
                !                t1=j
                !                t2=i
! equation (52) in Yang & Wang (2012). 
                !                mu2p=-pp+two*(t1*p1+t2*p2)
                !            else
                !                t1=i
                !                t2=j  
                !                mu2p=pp+two*(t1*p1+t2*p2)
                !            endif
                !        else        
                !            If(f12342.lt.zero)then        
                !                t1=i
                !                t2=j
! equation (52) in Yang & Wang (2012). 
                !                mu2p=pp+two*(t1*p1+t2*p2)
                !            endif
                !            If(f12342.gt.zero)then        
                !                t1=j
                !                t2=i      
! equation (52) in Yang & Wang (2012).    
                !                mu2p=-pp+two*(t1*p1+t2*p2)                              
                !            endif
                !        endif 
                !        !write(*,*)p,mu2p,abs(p-mu2p),pp,p1,p2! 
                !        If(dabs(p-mu2p).lt.1.D-4)goto 400
                !    enddo
                !enddo
                !11111111***************************************************** 
            !400 continue
            !tt1 = t1
            !tt2 = t2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! To determine the number of time_0 N_t1, N_t2 that the particle meets the
! two turn points mu_tp1, mu_tp2 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      If(mobseqmtp)then
                          p1_temp = zero
                          p2_temp = p_mt1_mt2
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(muobs.eq.mu_tp1)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_mt1_mt2
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          Else If(muobs.eq.mu_tp2)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_mt1_mt2
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          Endif 
                      Else
                          p1_temp = zero
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(f12342.gt.zero)then
                              p2_temp = PI2_p
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_mt1_mt2
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          ENDIF
                          If(f12342.lt.zero)then
                              p2_temp = PI1_p
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_mt1_mt2
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          ENDIF
                      Endif   
!========================================================================
          !write(unit = *, fmt = *)'************************************************************'
          !write(*,*)'sssssssssssssssssstttttt1= ',t1,tt1
          !write(*,*)'sssssssssssssssssstttttt2= ',t2,tt2
          !If (t1/=tt1 .or. t2/=tt2) then 
          !!    write(*,*)'sssssstttttt= ', p1_temp, p ,p2_temp, mu2p,pp
          !    write(*,*)'sssssstttttt= ', mobseqmtp,f12342, N_temp, p_mt1_mt2
          !    write(*,*)'sssssstttttt= ', p, -pp+two*(t1*p1+t2*p2),pp+two*(t1*p1+t2*p2)
          !    stop
          !endif
          !write(unit = *, fmt = *)'************************************************************'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            index_p4(1)=-1
            index_p4(2)=-2
            index_p4(3)=0
            index_p4(4)=-4
           !*****pp part***************************************
            If(lambda.ne.zero)then        
                cases_int=2
                call weierstrass_int_J3(tinf,tmu,dd,del,-tplus,&
                          b4,index_p4,dabs(pp),integ4,cases_int)
                call weierstrass_int_J3(tinf,tmu,dd,del,-tminus,&
                          b4,index_p4,dabs(pp),integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                pp_phi=lambda*(pp/(one-mu_tp1*mu_tp1)+&
                        integ4(2)*Wmup-integ14(2)*Wmum)
            else 
                pp_phi=zero                 
            endif
            cases_int=4 
            call weierstrass_int_J3(tinf,tmu,dd,del,h,b4,&
                        index_p4,dabs(pp),integ,cases_int)
            !write(*,*)'ynogk=',integ,cases_int!tinf,tmu,dd,del,h,b4,index_p4,abs(pp) 
! equation (71) in Yang & Wang (2012). 
            pp_t=a_spin**two*(pp*mu_tp1**two+integ(2)*mu_tp1*&
                               b0/two+integ(4)*b0**two/sixteen)
            !write(*,*)'ynogk=', mu_tp1,mu_tp2,b0,integ(2),integ(4)
           !*****p1 part***************************************
            If(t1.eq.0)then        
                p1_phi=zero
                p1_t=zero
            else  
                If(lambda.ne.zero)then  
                    IF(PI1_phi .EQ. zero)THEN
                        cases_int=2        
                        call weierstrass_int_J3(tinf,infinity,dd,del,&
                              -tplus,b4,index_p4,PI0,integ4,cases_int)
                        call weierstrass_int_J3(tinf,infinity,dd,del,&
                            -tminus,b4,index_p4,PI0,integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                        PI1_phi=lambda*(PI0/(one-mu_tp1**two)+&
                                integ4(2)*Wmup-integ14(2)*Wmum)
                    ENDIF 
! equation (51) in Yang & Wang (2012). 
                    p1_phi=PI1_phi-pp_phi 
                else 
                    p1_phi=zero             
                endif 
                IF(PI1_time .EQ. zero)THEN  
                    cases_int=4  
                    call weierstrass_int_J3(tinf,infinity,dd,del,&
                                h,b4,index_p4,PI0,integ,cases_int)
! equation (62) in Yang & Wang (2012). 
                    PI1_time=a_spin**two*(PI0*mu_tp1**two+integ(2)*&
                             mu_tp1*b0/two+integ(4)*b0**two/sixteen) 
                ENDIF
! equation (51) in Yang & Wang (2012). 
                p1_t=PI1_time-pp_t 
            endif 
          !*****p2 part***************************************
            If(t2.eq.0)then
                p2_phi=zero
                p2_t=zero
            else
                IF(lambda.ne.zero)then  
                    IF(PI2_phi .EQ. zero)THEN  
                        cases_int=2        
                        call weierstrass_int_J3(tp2,tinf,dd,del,&
                            -tplus,b4,index_p4,PI2_p,integ4,cases_int)
                        call weierstrass_int_J3(tp2,tinf,dd,del,&
                            -tminus,b4,index_p4,PI2_p,integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                        PI2_phi=lambda*(PI2_p/(one-mu_tp1*mu_tp1)+&
                                    integ4(2)*Wmup-integ14(2)*Wmum) 
                    ENDIF
! equation (51) in Yang & Wang (2012). 
                    p2_phi=PI2_phi+pp_phi
                ELSE
                    p2_phi=zero                
                ENDIF 

                IF(PI2_time .EQ. zero)THEN  
                    cases_int=4  
                    call weierstrass_int_J3(tp2,tinf,dd,del,h,b4,&
                                   index_p4,PI2_p,integ,cases_int)
! equation (71) in Yang & Wang (2012).  
                    PI2_time=a_spin**two*(PI2_p*mu_tp1**two+integ(2)*&
                               mu_tp1*b0/two+integ(4)*b0**two/sixteen) 
                ENDIF   
! equation (51) in Yang & Wang (2012).      
                p2_t=PI2_time+pp_t   
                !write(*,*)'ynogk=',tp2,tinf,h,dd
            endif   
        !**************************************************************
! equation (52) in Yang & Wang (2012). 
            !write(*,*)'ynogk=',pp_t,p1_t,p2_t,t1,t2
            If(mobseqmtp)then
                If(muobs.eq.mu_tp1)then  
                    phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    timet=-pp_t+two*(t1*p1_t+t2*p2_t)                
                else
                    phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)        
                    timet=pp_t+two*(t1*p1_t+t2*p2_t)                
                endif 
            else
                If(f12342.lt.zero)then
                    phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)        
                    timet=pp_t+two*(t1*p1_t+t2*p2_t)
                endif
                If(f12342.gt.zero)then                
                    phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    timet=-pp_t+two*(t1*p1_t+t2*p2_t)        
                endif
            endif
            count_num=count_num+1
        ELSE 
            If(f12343_1.eq.f12343.and.f12342_1.eq.f12342.and.lambda_1&
                       .eq.lambda.and.q_1.eq.q.and.sinobs_1.eq.sinobs&
            .and.muobs_1.eq.muobs.and.a_spin_1.eq.a_spin.and.scal_1.eq.scal)then  
        !***************************************************************************
                    t1=0
                    t2=0        
                    If(f12343.eq.zero.and.f12342.eq.zero.and.dabs(muobs).eq.one)then
                        sign_pth = zero
                        mucos = dsign(one,muobs)
                        timet=zero      !this is because that mu==1 for ever
                        phyt=zero       !this is because that mu==1 for ever,this is &
                                        !because that Theta_mu=-a^2(1-mu^2)
                        return          !so,mu must =+1 or -1 for ever.
                    endif
                    If(muobs.eq.zero.and.(dabs(lambda).lt.dabs(a_spin)).and.q.eq.zero)then
                        timet=zero
                        phyt=zero
                        mucos=zero
                        sign_pth = zero
                        return                        
                    endif          
                    If(mu_tp1.eq.zero)then
                        !photons are confined in the equatorial plane, &
                        !so the integrations about \theta are valished.
                        timet=zero
                        phyt=zero
                        mucos=zero
                        sign_pth = zero
                        return
                    endif
                    If(a_spin.eq.zero)then
                        timet=zero
                        CALL phyt_schwatz(p,f12343,f12342,lambda,q,&
                                 sinobs,muobs,scal,phyt,mucos,t1,t2,sign_pth)
                        return
                    endif
         
                    tmu=weierstrassP(p+PI01,g2,g3,dd,del)
                    mucos=mu_tp1+b0/(four*tmu-b1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Now we get the value of parameter: sign_pth
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    if ( u <= zero ) then
                        sign_pth = -one
                    else
                        u = dmod(u, period_wp)
                        if ( u <= half_period_wp ) then
                            sign_pth = one
                        else
                            sign_pth = - one
                        endif
                    endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    !get numbers of turn points of t1 and t2.
                    !111111111************************************************************  
                        index_p4(1)=0
                        cases_int=1    
                        call weierstrass_int_J3(tinf,tmu,dd,del,a4,b4,&
                                       index_p4,rff_p,integ4,cases_int)                
                        pp=integ4(1)
                        p1=PI0-pp
                        p2=p_mt1_mt2-p1
                        !p1=zero
                        !p2=zero
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        !Do j=0,10 
                        !    Do i=j,j+1
                        !        If(mobseqmtp)then
                        !            If(muobs.eq.mu_tp1)then
                        !                t1=j
                        !                t2=i
! equation (54) in Yang & Wang (2012). 
                        !                mu2p=-pp+two*(t1*p1+t2*p2) 
                        !            else
                        !                t1=i
                        !                t2=j 
                        !                mu2p=pp+two*(t1*p1+t2*p2) 
                        !            endif
                        !        else        
                        !            If(f12342.lt.zero)then        
                        !                t1=i
                        !                t2=j
! equation (54) in Yang & Wang (2012). 
                        !                mu2p=pp+two*(t1*p1+t2*p2)
                        !            endif
                        !            If(f12342.gt.zero)then        
                        !                t1=j
                        !                t2=i     
! equation (54) in Yang & Wang (2012).      
                        !                mu2p=-pp+two*(t1*p1+t2*p2)                             
                        !            endif
                        !        endif  
                        !        !write(*,*)p,mu2p,t1,t2
                        !        If(dabs(p-mu2p).lt.1.D-4)goto 410
                        !    enddo
                        !enddo
                        !11111111********************************************* 
                    !410 continue
                    !tt1 = t1
                    !tt2 = t2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! To determine the number of time_0 N_t1, N_t2 that the particle meets the
! two turn points mu_tp1, mu_tp2 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      If(mobseqmtp)then
                          p1_temp = zero
                          p2_temp = p_mt1_mt2
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(muobs.eq.mu_tp1)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_mt1_mt2
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          Else If(muobs.eq.mu_tp2)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_mt1_mt2
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          Endif 
                      Else
                          p1_temp = zero
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(f12342.gt.zero)then
                              p2_temp = PI2_p
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_mt1_mt2
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          ENDIF
                          If(f12342.lt.zero)then
                              p2_temp = PI1_p
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_mt1_mt2
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          ENDIF
                      Endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    index_p4(1)=-1
                    index_p4(2)=-2
                    index_p4(3)=0
                    index_p4(4)=-4 
                    !*****pp parts************************************
                    If(lambda.ne.zero)then
                        cases_int=2        
                        call weierstrass_int_J3(tinf,tmu,dd,del,&
                             -tplus,b4,index_p4,dabs(pp),integ4,cases_int)
                        call weierstrass_int_J3(tinf,tmu,dd,del,&
                             -tminus,b4,index_p4,dabs(pp),integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                        pp_phi=lambda*(pp/(one-mu_tp1**two)+&
                              integ4(2)*Wmup-integ14(2)*Wmum) 
                    else 
                        pp_phi=zero         
                    endif
                    cases_int=4
                    call weierstrass_int_J3(tinf,tmu,dd,del,h,&
                           b4,index_p4,dabs(pp),integ,cases_int)
! equation (71) in Yang & Wang (2012). 
                    pp_t=a_spin**two*(pp*mu_tp1**two+integ(2)*&
                        mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                    !*****p1 parts************************************
                    If(t1.eq.0)then        
                        p1_phi=zero
                        p1_t=zero
                    else  
                        If(lambda.ne.zero)then  
                            IF(PI1_phi .EQ. zero)THEN
                                cases_int=2        
                                call weierstrass_int_J3(tinf,infinity,dd,del,-tplus,&
                                                    b4,index_p4,PI0,integ4,cases_int)
                                call weierstrass_int_J3(tinf,infinity,dd,del,-tminus,&
                                                    b4,index_p4,PI0,integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                                PI1_phi=lambda*(PI0/(one-mu_tp1**two)+&
                                        integ4(2)*Wmup-integ14(2)*Wmum)
                            ENDIF
! equation (51) in Yang & Wang (2012). 
                            p1_phi=PI1_phi-pp_phi 
                        else 
                            p1_phi=zero            
                        endif  
                        IF(PI1_time .EQ. zero)THEN
                            cases_int=4  
                            call weierstrass_int_J3(tinf,infinity,dd,&
                                del,h,b4,index_p4,PI0,integ,cases_int)
! equation (71) in Yang & Wang (2012). 
                            PI1_time=a_spin**two*(PI0*mu_tp1**two+integ(2)&
                                   *mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                        ENDIF
! equation (51) in Yang & Wang (2012). 
                        p1_t=PI1_time-pp_t
                    endif 
                    !*****p2 parts************************************
                    If(t2.eq.0)then
                        p2_phi=zero
                        p2_t=zero
                    else        
                        If(lambda.ne.zero)then
                            IF(PI2_phi .EQ. zero)THEN 
                               cases_int=2
                               call weierstrass_int_J3(tp2,tinf,dd,del,-tplus,&
                                            b4,index_p4,PI2_p,integ4,cases_int)
                               call weierstrass_int_J3(tp2,tinf,dd,del,-tminus,&
                                            b4,index_p4,PI2_p,integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                               PI2_phi=lambda*(PI2_p/(one-mu_tp1**two)+&
                                         integ4(2)*Wmup-integ14(2)*Wmum)
                            ENDIF
! equation (51) in Yang & Wang (2012). 
                            p2_phi=PI2_phi+pp_phi 
                        else 
                            p2_phi=zero                     
                        endif
                        IF(PI2_time .EQ. zero)THEN 
                            cases_int=4
                            call weierstrass_int_J3(tp2,tinf,dd,del,h,b4,&
                                           index_p4,PI2_p,integ,cases_int)
! equation (71) in Yang & Wang (2012). 
                            PI2_time=a_spin**two*(PI2_p*mu_tp1**two+&
                                     integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                        ENDIF
! equation (51) in Yang & Wang (2012). 
                        p2_t=PI2_time+pp_t 
                    endif 
                !**************************************************************
! equation (52) in Yang & Wang (2012). 
                    If(mobseqmtp)then
                        If(muobs.eq.mu_tp1)then  
                            phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timet=-pp_t+two*(t1*p1_t+t2*p2_t)                
                        else
                            phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)        
                            timet=pp_t+two*(t1*p1_t+t2*p2_t)                
                        endif 
                    else
                        If(f12342.lt.zero)then
                            phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)        
                            timet=pp_t+two*(t1*p1_t+t2*p2_t)
                        endif
                        If(f12342.gt.zero)then                
                            phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timet=-pp_t+two*(t1*p1_t+t2*p2_t)
                        endif
                    endif                
            else
                count_num=1
                goto 30          
            endif                               
        ENDIF
        !write(*,*)'ff1=',phyt,timet,pp_phi,p1_phi,p2_phi,t1,t2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !write(unit = *, fmt = *)'************************************************************'
          !write(*,*)'sssssssssssssssssstttttt1= ',t1,tt1
          !write(*,*)'sssssssssssssssssstttttt2= ',t2,tt2
          !If (t1/=tt1 .or. t2/=tt2) then
          !    write(*,*)'sssssstttttt= ', p1_temp, p ,p2_temp, mu2p,pp
          !    write(*,*)'sssssstttttt= ', mobseqmtp,f12342, N_temp, p_mt1_mt2
          !    write(*,*)'sssssstttttt= ', p, -pp+two*(t1*p1+t2*p2),pp+two*(t1*p1+t2*p2)
          !    write(*,*)'mmmmmmmm======',-pp_t+two*(t1*p1_t+t2*p2_t),pp_t+two*(t1*p1_t+t2*p2_t)
          !    write(*,*)'mmmmmmmm======',-pp_t+two*(tt1*p1_t+tt2*p2_t),pp_t+two*(tt1*p1_t+tt2*p2_t)
          !    write(*,*)'sssss',-pp_phi+two*(t1*p1_phi+t2*p2_phi),pp_phi+two*(t1*p1_phi+t2*p2_phi)
          !    write(*,*)'sssss',-pp_phi+two*(tt1*p1_phi+tt2*p2_phi),pp_phi+two*(tt1*p1_phi+tt2*p2_phi)
          !    stop
          !endif
          !write(unit = *, fmt = *)'************************************************************'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RETURN
      END SUBROUTINE INTTPART

!********************************************************************************************
      SUBROUTINE phyt_schwatz(p,f3,f2,lambda,q,sinobs,muobs,&
                         scal,phyc_schwatz,mucos,t1,t2,sign_pth)
!******************************************************************************************** 
!*     PURPOSE:  Computes \mu part of integrals in coordinates \phi, expressed by equation (72) 
!*               in Yang & Wang (2012) with zero spin of black hole.    
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f2-------------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f3-------------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phyc_schwatz-----------value of integral \phi_\theta expressed by equation (71) in
!*                              Yang & Wang (2012).   
!*               mucos----------value of function \mu(p) with zero spin.
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.            
!*     ROUTINES CALLED: schwatz_int
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
        USE constants
        implicit none
        Double precision phyc_schwatz,f3,f2,p,sinobs,muobs,phyr_schwatz,pp,p1,p2,mu,&
                    lambda,AA,BB,AAi,AAim,scal,q,mu1,mu2,mu_tp1,mu_tp2,PI2,sign_pth,&
                    mu2p,mucos,Pt,PI1,PI1_phi,PI2_phi,f3_1,f2_1,lambda_1,q_1,u,&
                    sinobs_1,muobs_1,scal_1,pp_phi,p1_phi,p2_phi,p1_temp,p2_temp
        integer :: t1,t2,cases,i,j,count_num=1,N_temp,tt1,tt2
        logical :: err,mobseqmtp
        save :: PI1,PI2,PI1_phi,PI2_phi,Pt,f3_1,f2_1,lambda_1,q_1,pp_phi,&
                p1_phi,p2_phi,sinobs_1,muobs_1,scal_1,mobseqmtp,AA,BB,&
                mu_tp1,mu_tp2 

60 continue 
      IF(count_num .EQ. 1)THEN
          f3_1=f3
          f2_1=f2
          lambda_1=lambda
          q_1=q
          muobs_1=muobs
          sinobs_1=sinobs 
          scal_1=scal
          t1=0
          t2=0         
          mobseqmtp=.false.
          If(q.gt.zero)then        
              AA=dsqrt((lambda**two+q)/q)
              BB=dsqrt(q)
        !*****************************************************
              If(f2.lt.zero)then
                  u = dasin(muobs*AA)+p*BB*AA
                  !mu=dsin(dasin(muobs*AA)+p*BB*AA)/AA
                  mu=dsin( u )/AA
                  u = dmod(u + halfpi, twopi)
                  If ( zero <= u .AND. u <= pi ) then
                      sign_pth = -one
                  else
                      sign_pth = one
                  endif
              else
                  If(f2.eq.zero)then
                      u = p*AA*BB
                      !mu=dcos(p*AA*BB)*muobs
                      mu=dcos(u)*muobs
                      u = dmod(u, twopi)
                      If ( zero <= u .AND. u <= pi ) then
                          sign_pth = dsign(one, muobs)
                      else
                          sign_pth = -dsign(one, muobs)
                      endif  
                  else 
                      u = dasin(muobs*AA)-p*BB*AA                     
                      !mu=dsin(dasin(muobs*AA)-p*AA*BB)/AA        
                      mu = dsin( u ) / AA 
                      u = dmod(u - halfpi, twopi)
                      If ( -pi <= u .AND. u <= zero ) then
                          sign_pth = one
                      else
                          sign_pth = -one
                      endif       
                  endif        
              endif
              mucos = mu  
        !****************************************************
              If(f2.ne.zero)then
                  mu_tp1=dsqrt(q/(lambda**two+q))
                  mu_tp2=-mu_tp1        
              else
                  mu_tp1=dabs(muobs)
                  mu_tp2=-mu_tp1
                  mobseqmtp=.true.
              endif
              If(dabs(muobs).eq.one)mobseqmtp=.true.        

              If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, 
              !so the integrations about !\theta are valished.
                  sign_pth=zero
                  phyc_schwatz=zero
                  return
              endif

              !***************************************************
              PI1=(PI/two-dasin(muobs/mu_tp1))*mu_tp1/BB        
              Pt=PI*mu_tp1/BB
              PI2=Pt - PI1      
              pp=(dasin(mu/mu_tp1)-dasin(muobs/mu_tp1))*mu_tp1/BB        
              p1=PI1-pp
              p2=Pt-p1 
              PI1_phi=zero
              PI2_phi=zero
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              Do j=0,100 
!                  Do i=j,j+1
!                      If(mobseqmtp)then
!                          If(muobs.eq.mu_tp1)then
!                              t1=j
!                              t2=i
!                              mu2p=-pp+two*(t1*p1+t2*p2)
!                          else
!                              t1=i
!                              t2=j
!                              mu2p=pp+two*(t1*p1+t2*p2)
!                          endif
!                      else        
!                          If(f2.lt.zero)then        
!                              t1=i
!                              t2=j
!                              mu2p=pp+two*(t1*p1+t2*p2)
!                          endif
!                          If(f2.gt.zero)then        
!                              t1=j
!                              t2=i    
!                              mu2p=-pp+two*(t1*p1+t2*p2)                            
!                          endif
!                      endif  
!                      If(dabs(p-mu2p).lt.1.D-4)goto 300
!                  enddo
!              enddo
!              !*************************************************************** 
!              300 continue
!              tt1 = t1
!              tt2 = t2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! To determine the number of time_0 N_t1, N_t2 that the particle meets the
! two turn points mu_tp1, mu_tp2 respectively.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      If(mobseqmtp)then
                          p1_temp = zero
                          p2_temp = pt
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(muobs.eq.mu_tp1)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + pt
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          Else If(muobs.eq.mu_tp2)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + pt
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          Endif 
                      Else
                          p1_temp = zero
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(f2.gt.zero)then
                              p2_temp = PI2
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + pt
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          ENDIF
                          If(f2.lt.zero)then
                              p2_temp = PI1
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + pt
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          ENDIF
                      Endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              If(lambda.eq.zero)then 
                  phyc_schwatz = zero
                  return
              endif
              pp_phi=lambda*schwatz_int(muobs,mu,AA)/BB
              If(t1.eq.0)then
                  p1_phi=zero
              else
                  IF(PI1_phi .eq. zero)THEN
                      PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB
                  ENDIF
                  p1_phi=PI1_phi-pp_phi       
              endif
              If(t2.eq.0)then
                  p2_phi=zero
              else
                  IF(PI2_phi .EQ. zero)THEN
                      PI2_phi=lambda*schwatz_int(mu_tp2,muobs,AA)/BB
                  ENDIF
                  p2_phi=PI2_phi+pp_phi 
              endif
              If(mobseqmtp)then
                  If(muobs.eq.mu_tp1)then  
                      phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)                
                  else
                      phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)                
                  endif         
              else
                 If(f2.lt.zero)then
                      phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                 endif         
                 If(f2.gt.zero)then
                      phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                 endif          
              endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !write(unit = *, fmt = *)'************************************************************'
          !write(*,*)'sssssssssssssssssstttttt1= ',t1,tt1
          !write(*,*)'sssssssssssssssssstttttt2= ',t2,tt2
          !If (t1/=tt1 .or. t2/=tt2) then
          !    write(*,*)'sssssstttttt1= ', p1_temp, p ,p2_temp, mu2p,pp
          !    write(*,*)'sssssstttttt2= ', mobseqmtp,f2, N_temp, pt
          !    write(*,*)'sssssstttttt3= ', p, -pp+two*(t1*p1+t2*p2),pp+two*(t1*p1+t2*p2) 
          !    write(*,*)'sssss',-pp_phi+two*(t1*p1_phi+t2*p2_phi),pp_phi+two*(t1*p1_phi+t2*p2_phi)
          !    write(*,*)'sssss',-pp_phi+two*(tt1*p1_phi+tt2*p2_phi),pp_phi+two*(tt1*p1_phi+tt2*p2_phi)
          !    stop
          !endif
          !write(unit = *, fmt = *)'************************************************************'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else
              !write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is an affending',&
              !                'value, the program should be',&  
              !                'stoped! and q = ',q
              !stop
              mucos=muobs 
              t1 = 0
              t2 = 0
              phyc_schwatz = zero
              sign_pth = zero
          endif        
      ELSE
          IF(f3_1.eq.f3.and.f2_1.eq.f2.and.lambda_1.eq.lambda.and.q_1.eq.q.and.sinobs_1.eq.sinobs&
          .and.muobs_1.eq.muobs.and.scal_1.eq.scal)THEN
              If(q.gt.zero)then         
        !*****************************************************
                  If(f2.lt.zero)then 
                      u = dasin(muobs*AA) + p*BB*AA
                      !mu=dsin(dasin(muobs*AA)+p*BB*AA)/AA
                      mu = dsin( u ) / AA
                      u = dmod(u + halfpi, twopi)
                      If ( zero <= u .AND. u <= pi ) then
                          sign_pth = -one
                      else
                          sign_pth = one
                      endif        
                  else
                      If(f2.eq.zero)then 
                          u = p*AA*BB
                          !mu=dcos(p*AA*BB)*muobs
                          mu=dcos(u)*muobs
                          u = dmod(u, twopi)
                          If ( zero <= u .AND. u <= pi ) then
                              sign_pth = dsign(one, muobs)
                          else
                              sign_pth = -dsign(one, muobs)
                          endif 
                      else 
                          u = dasin(muobs*AA) - p*BB*AA
                          !mu=dsin(dasin(muobs*AA)-p*AA*BB)/AA        
                          mu = dsin( u ) / AA 
                          u = dmod(u - halfpi, twopi)
                          If ( -pi <= u .AND. u <= zero ) then
                              sign_pth = one
                          else
                              sign_pth = -one
                          endif        
                      endif        
                  endif
                  mucos = mu  
        !****************************************************  
                  If(mu_tp1.eq.zero)then
                  !photons are confined in the equatorial plane, 
                  !so the integrations about !\theta are valished.
                      t1=0
                      t2=0
                      sign_pth=zero
                      phyc_schwatz=zero
                      return
                  endif

                  !***************************************************  
                  pp=(dasin(mu/mu_tp1)-dasin(muobs/mu_tp1))*mu_tp1/BB        
                  p1=PI1-pp
                  p2=Pt-p1  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  Do j=0,100 
!                      Do i=j,j+1
!                      If(mobseqmtp)then
!                          If(muobs.eq.mu_tp1)then
!                              t1=j
!                              t2=i
!                              mu2p=-pp+two*(t1*p1+t2*p2) 
!                          else
!                              t1=i
!                              t2=j
!                              mu2p=pp+two*(t1*p1+t2*p2) 
!                          endif
!                      else        
!                          If(f2.lt.zero)then        
!                              t1=i
!                              t2=j
!                              mu2p=pp+two*(t1*p1+t2*p2)
!                          endif
!                          If(f2.gt.zero)then        
!                              t1=j
!                              t2=i      
!                              mu2p=-pp+two*(t1*p1+t2*p2)                            
!                          endif
!                      endif   
!                      If(dabs(p-mu2p).lt.1.D-4)goto 310
!                      enddo
!                  enddo
!                  !************************************************************ 
!                  310 continue
!                  tt1 = t1
!                  tt2 = t2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! To determine the number of time_0 N_t1, N_t2 that the particle meets the
! two turn points mu_tp1, mu_tp2 respectively.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      If(mobseqmtp)then
                          p1_temp = zero
                          p2_temp = pt
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(muobs.eq.mu_tp1)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + pt
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          Else If(muobs.eq.mu_tp2)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + pt
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          Endif 
                      Else
                          p1_temp = zero
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(f2.gt.zero)then
                              p2_temp = PI2
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + pt
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          ENDIF
                          If(f2.lt.zero)then
                              p2_temp = PI1
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + pt
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          ENDIF
                      Endif  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  If(lambda.eq.zero)then 
                      phyc_schwatz = zero
                      return
                  endif
                  pp_phi=lambda*schwatz_int(muobs,mu,AA)/BB
                  If(t1.eq.0)then
                      p1_phi=zero
                  else
                      IF(PI1_phi .eq. zero)THEN
                          PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB
                      ENDIF
                      p1_phi=PI1_phi-pp_phi       
                  endif
                  If(t2.eq.0)then
                      p2_phi=zero
                  else
                      IF(PI2_phi .EQ. zero)THEN
                          PI2_phi=lambda*schwatz_int(mu_tp2,muobs,AA)/BB
                      ENDIF
                      p2_phi=PI2_phi+pp_phi 
                  endif
                  If(mobseqmtp)then
                      If(muobs.eq.mu_tp1)then  
                          phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                      else
                          phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                      endif         
                  else
                      If(f2.lt.zero)then
                          phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                      endif         
                      If(f2.gt.zero)then
                          phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                      endif          
                  endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !write(unit = *, fmt = *)'************************************************************'
          !write(*,*)'sssssssssssssssssstttttt1= ',t1,tt1
          !write(*,*)'sssssssssssssssssstttttt2= ',t2,tt2
          !If (t1/=tt1 .or. t2/=tt2) then
          !    write(*,*)'sssssstttttt1= ', p1_temp, p ,p2_temp, mu2p,pp
          ! !   write(*,*)'sssssstttttt2= ', mobseqmtp,f2, N_temp, pt
          !    write(*,*)'sssssstttttt3= ', p, -pp+two*(t1*p1+t2*p2),pp+two*(t1*p1+t2*p2) 
          !    write(*,*)'sssss',-pp_phi+two*(t1*p1_phi+t2*p2_phi),pp_phi+two*(t1*p1_phi+t2*p2_phi)
          !    write(*,*)'sssss',-pp_phi+two*(tt1*p1_phi+tt2*p2_phi),pp_phi+two*(tt1*p1_phi+tt2*p2_phi)
          !    stop
          !endif
          !write(unit = *, fmt = *)'************************************************************'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              else
                  !write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is a affending',&
                  !                'value, the program should be',&  
                  !                'stoped! and q = ',q
                  !stop
                  mucos=muobs 
                  t1 = 0
                  t2 = 0
                  phyc_schwatz = zero
                  sign_pth=zero
              endif                  
          ELSE
              count_num=1
              goto 60
          ENDIF
      ENDIF                                
      return
      End SUBROUTINE phyt_schwatz 
!************************************************************************* 
      Function schwatz_int(y,x,AA) 
!************************************************************************* 
!*     PURPOSE:  Computes \int^x_y dt/(1-t^2)/sqrt(1-AA^2*t^2) and AA .gt. 1  
!*     INPUTS:   components of above integration.      
!*     OUTPUTS:  valve of integral.             
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
      USE constants
      implicit none
      Double precision y,x,yt,xt,AA,schwatz_int,ppx,ppy,A2,tp  
      logical :: inverse

      xt=x
      yt=y        
      If(yt.eq.xt)then
          schwatz_int=0.D0
          return        
      endif
      If(dabs(AA).ne.one)then
          A2=AA*AA
          ppx=datan(dsqrt(A2-one)*xt/dsqrt(dabs(one-A2*xt*xt)))
          ppy=datan(dsqrt(A2-one)*yt/dsqrt(dabs(one-A2*yt*yt)))
          schwatz_int=(ppx-ppy)/dsqrt(A2-one) 
      ELse
          If(dabs(xt).eq.one)then
              schwatz_int=infinity
          Else
              If(dabs(yt).eq.one)then
                  schwatz_int=-infinity
              Else
                  ppx=xt/dsqrt(dabs(one-xt*xt))
                  ppy=yt/dsqrt(dabs(one-yt*yt))
                  schwatz_int=ppx-ppy 
              endif        
          Endif           
      Endif
      return
      End Function schwatz_int   

!********************************************************************************************
      SUBROUTINE INTRPART(p,f1234r,f1234t,lambda,q,sinobs,muobs,a_spin,&
                            robs,scal,phyr,timer,affr,r_coord,t1,t2,sign_pr)
!******************************************************************************************** 
!*     PURPOSE:  Computes r part of integrals in coordinates \phi, t and affine parameter \sigma,
!*               expressed by equations (62), (63), (65), (67), (68) and (69) in Yang & Wang (2012).    
!*     INPUTS:   p--------------independent variable, which must be nonnegative. 
!*               f1234r---------p_r, which is the r component of four momentum of a photon 
!*                              measured under the LNRF, see equation (83) in Yang & Wang (2012).
!*               f1234t---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phyr-----------value of integral \phi_r expressed by equation (65) or (69) in
!*                              Yang & Wang (2012).  
!*               affr-----------value of integral \sigma_r expressed by equation (62) or (67) in
!*                              Yang & Wang (2012).  
!*               timer----------value of integral t_r expressed by equation (63) or (68) in
!*                              Yang & Wang (2012).  
!*               r_coord--------value of function r(p).
!*               t1,t2----------number of times of photon meets turning points r_tp1 and r_tp2
!*                              respectively.
!*               sign_pr--------the sign of r component of 4-momentum of the photon.
!*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
        USE constants
        IMPLICIT NONE
        DOUBLE PRECISION phyr,radius,re,a,B,p,sinobs,muobs,a_spin,rhorizon,q,lambda,integ,integ4(4),&
                          bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tobs,tp,pp,p1,p2,PI0,p1I0,p1J1,p1J2,E_add,E_m,&
                          u,v,w,L1,L2,thorizon,m2,pinf,sn,cn,dn,r_add,r_m,B_add,B_m,D_add,D_m,&
                         y,x,f1,g1,h1,f2,h2,a5,b5,a4,b4,integ0,integ1,integ2,robs,ttp,sign_pr,&
                         PI1,PI2,scal,tinf,integ04(4),integ14(4),integ5(5),integ15(5),pp2,up,&
                         r_tp1,r_tp2,t_inf,tp2,f1234r,f1234t,p_temp,PI0_obs_inf,PI0_total,PI0_obs_hori,&
                         PI0_obs_tp2,PI01,timer,affr,r_coord,cr,dr,rff_p,p_t1_t2,half_periodwp,&
                         Ap,Am,h,wp,wm,wbarp,wbarm,hm,hp,pp_time,pp_phi,pp_aff,p1_phi,p1_time,p1_aff,&
                         p2_phi,p2_time,p2_aff,time_temp,sqt3,p_tp1_tp2,PI2_p,PI1_p,periodwp,&
                         PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_aff,PI2_aff,p1_temp,p2_temp
        DOUBLE PRECISION f1234r_1,f1234t_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,scal_1
        COMPLEX*16 bb(1:4),dd(3)
        INTEGER :: reals,i,j,t1,t2,p5,p4,index_p4(4),index_p5(5),del,cases_int,&
                   cases,count_num=1,N_temp,tt1,tt2
        LOGICAL :: robs_eq_rtp,indrhorizon
        SAVE :: f1234r_1,f1234t_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,&
                robs_1,scal_1,rhorizon,r_add,r_m,a4,b4,B_add,half_periodwp,periodwp,&
                B_m,robs_eq_rtp,indrhorizon,r_tp1,r_tp2,reals,cases,bb,b0,b1,b2,b3,g2,g3,tobs,thorizon,&
                tp2,tinf,dd,E_add,E_m,D_add,D_m,PI0,PI0_obs_inf,PI0_total,PI0_obs_hori,PI0_obs_tp2,del,&
                u,v,w,L1,L2,m2,t_inf,pinf,f1,g1,h1,f2,h2,b5,Ap,Am,h,wp,wm,wbarp,wbarm,hm,hp,a5,cc,&
                PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_aff,PI2_aff,PI2_p,PI1_p,p_tp1_tp2,sqt3          

  40 continue        
        If(count_num.eq.1)then
                f1234r_1=f1234r        
                f1234t_1=f1234t
                lambda_1=lambda
                q_1=q
                muobs_1=muobs
                sinobs_1=sinobs
                a_spin_1=a_spin
                robs_1=robs
                scal_1=scal
            !************************************************************************************
                rhorizon=one+dsqrt(one-a_spin**two)
! equation (64) in Yang & Wang (2012).
                r_add=rhorizon
                r_m=one-dsqrt(one-a_spin**two)      
! equation (64) in Yang & Wang (2012).  
                B_add=(two*r_add-a_spin*lambda)/(r_add-r_m)        
                B_m=(two*r_m-a_spin*lambda)/(r_add-r_m)      
! equation (64) in Yang & Wang (2012).
                Ap=(r_add*(four-a_spin*lambda)-two*a_spin**two)/dsqrt(one-a_spin**two)
                Am=(r_m*(four-a_spin*lambda)-two*a_spin**two)/dsqrt(one-a_spin**two)
                b4=one
                a4=zero
                cc=a_spin**2-lambda**2-q
                robs_eq_rtp=.false.
                indrhorizon=.false.
                call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
                                  reals,robs_eq_rtp,indrhorizon,cases,bb)
! equation (55) in Yang & Wang (2012).
                PI1_phi=zero
                PI2_phi=zero
                PI1_time=zero
                PI2_time=zero
                PI1_aff=zero
                PI2_aff=zero 
!** R(r)=0 has real roots and turning points exists in radial r.
                If(reals.ne.0)then  
! equations (35)-(38) in Yang & Wang (2012).
                    b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+&
                                              two*(q+(lambda-a_spin)**2)
                    b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
                    b2=four/three*r_tp1
                    b3=one
                    g2=three/four*(b1**2-b0*b2)
                    g3=one/16.D0*(three*b0*b1*b2-two*b1**three-b0**two*b3)
! equation (39) in Yang & Wang (2012).
                    If(robs-r_tp1.ne.zero)then        
                        tobs=b0/four/(robs-r_tp1)+b1/four
                    else
                        tobs=infinity
                    endif 
                    If(rhorizon-r_tp1.ne.zero)then
                        thorizon=b1/four+b0/four/(rhorizon-r_tp1)
                    else
                        thorizon=infinity         
                    endif
                    tp2=b0/four/(r_tp2-r_tp1)+b1/four   
                    tinf=b1/four
                    h=-b1/four        
! equation (64), (66) and (70) in Yang & Wang (2012).
                    call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)        
                        E_add=b0/(four*(r_add-r_tp1))+b1/four
                          E_m=b0/(four*(r_m-r_tp1))+b1/four
                        D_add=b0/(four*(r_tp1-r_add)**2)
                          D_m=b0/(four*(r_tp1-r_m)**2)

                               wp=one/(r_tp1-r_add)
                               wm=one/(r_tp1-r_m)
                            wbarp=b0/four/(r_tp1-r_add)**two
                            wbarm=b0/four/(r_tp1-r_m)**two
                               hp=b0/four/(r_add-r_tp1)+b1/four
                               hm=b0/four/(r_m-r_tp1)+b1/four

                    index_p4(1)=0
                    cases_int=1
                    call weierstrass_int_J3(tobs,infinity,dd,del,a4,&
                                 b4,index_p4,rff_p,integ04,cases_int) 
! equation (42) in Yang & Wang (2012).
                    PI0=integ04(1)   
                    select case(cases)
                    CASE(1)
                        If(f1234r .ge. zero)then !**photon will goto infinity.
                            index_p4(1)=0
                            cases_int=1
                            call weierstrass_int_J3(tinf,tobs,dd,del,a4,&
                                     b4,index_p4,rff_p,integ04,cases_int)
                            PI0_obs_inf=integ04(1)
                            If(p.lt.PI0_obs_inf)then    
! equation (41) in Yang & Wang (2012).    
                                tp=weierstrassP(p+PI0,g2,g3,dd,del) 
                                r_coord = r_tp1+b0/(four*tp-b1)
                                pp=-p                                  
                            else
                                tp=tinf! !Goto infinity, far away. 
                                r_coord = infinity
                                pp=-PI0_obs_inf 
                            endif
                            t1=0
                            t2=0   
                            sign_pr = one       
                        ELSE 
                            If(.not.indrhorizon)then
                                index_p4(1)=0
                                cases_int=1 
                                call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,&
                                                   index_p4,rff_p,integ14,cases_int)
                                PI0_total=PI0+integ14(1)
                                t2=0
                                If(p.le.PI0)then
                                    t1=0
                                    pp=p  
! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    sign_pr = - one
                                else
                                    t1=1
                                    PI1_p=PI0        
                                    If(p.lt.PI0_total)then   
! equation (41) in Yang & Wang (2012).     
                                        tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=two*PI0-p
                                        p1=dabs(p-PI0)
                                    else        
                                        tp=tinf !Goto infinity, far away.
                                        r_coord = infinity
                                        pp=-PI0_total+two*PI0 
                                        p1=pI0_total-PI0
                                    endif 
                                    sign_pr = one       
                                endif        
                            ELSE
!f1234r<0, photon will fall into black hole unless something encountered. 
                                index_p4(1)=0                
                                cases_int=1
                                call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,&
                                                   index_p4,rff_p,integ04,cases_int)
                                PI0_obs_hori=integ04(1)
                                If(p.lt.PI0_obs_hori)then     
! equation (41) in Yang & Wang (2012).   
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)        
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p                                  
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0  
                                sign_pr = - one
                            ENDIF
                        ENDIF         
                    CASE(2)
                        If(.not.indrhorizon)then
                            If(f1234r.lt.zero)then
                                PI01=-PI0
                            else
                                PI01=PI0        
                            endif
! equation (41) in Yang & Wang (2012).
                            tp=weierstrassP(p+PI01,g2,g3,dd,del)
                            r_coord = r_tp1+b0/(four*tp-b1)     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
                            index_p4(1)=0
                            cases_int=1 
                            call weierstrass_int_J3(tp2,infinity,dd,del,a4,&
                                      b4,index_p4,rff_p,integ14,cases_int)      
                            call weierstrass_int_J3(tobs,tp,dd,del,a4,b4,&
                                          index_p4,rff_p,integ4,cases_int)
                            half_periodwp = integ14(1)
                            periodwp = two * half_periodwp
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! this part of the code aims on obtaining the sign of the r component of  !|
               ! 4-momentum of the photon.                                               !|
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    If (f1234r > zero) then                                              !|
                        up = dmod(p + PI01, periodwp)                                    !|
                        If ( up <= half_periodwp ) then                                  !|
                            sign_pr = one                                                !|
                        else                                                             !|
                            sign_pr = - one                                              !|
                        endif                                                            !|
                    Else If (f1234r < zero) then                                         !|
                        up = dmod(p + PI01 + half_periodwp, periodwp)                    !|
                        If (up < half_periodwp) then                                     !|
                            sign_pr = -one                                               !|
                        Else                                                             !|
                            sign_pr = one                                                !|
                        endif                                                            !|
                    Else                                                                 !|
                        If (robs == r_tp1) then                                          !|
                            If ( dmod(p, periodwp) <= half_periodwp ) then               !|
                                sign_pr = one                                            !|
                            else                                                         !|
                                sign_pr = -one                                           !|
                            endif                                                        !|
                        else                                                             !| 
                            If ( dmod(p + half_periodwp, periodwp) &                     !|
                                                   <= half_periodwp ) then               !|
                                sign_pr = one                                            !|
                            else                                                         !|
                                sign_pr = -one                                           !|
                            endif                                                        !| 
                        endif                                                            !|
                    Endif                                                                !| 
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                            pp=integ4(1)
! equation (57) in Yang & Wang (2012).
                            p_tp1_tp2=integ14(1) 
                            PI2_p=p_tp1_tp2-PI0 
                            PI1_p=PI0 
                            p1=PI0-pp
                            p2=p_tp1_tp2-p1 
                            !p1=zero
                            !p2=zero
                        !*************************************************************************************
! equation (58) in Yang & Wang (2012).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! To determine t1 and t2, which are the times of the photon meets the two
! turnning points: r_tp1 and r_tp2 respectively.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            Do j=0,100
                                Do i=j,j+1
                                    If(robs_eq_rtp)then        
                                        If(robs.eq.r_tp1)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        else
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    else
                                        If(f1234r.gt.zero)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        endif
                                        If(f1234r.lt.zero)then
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    endif  
                                    If(dabs(p-p_temp).lt.1.D-4)goto 200
                                Enddo
                            Enddo
!*************************************************************************************
                        200     continue
                        tt1 = t1
                        tt2 = t2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This code section is added here at 2017--11--16 am 11:23
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      If(robs_eq_rtp)then
                          p1_temp = zero
                          p2_temp = p_tp1_tp2
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(robs.eq.r_tp1)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_tp1_tp2
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          Else If(robs.eq.r_tp2)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_tp1_tp2
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          Endif 
                      Else
                          p1_temp = zero
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(f1234r.gt.zero)then
                              p2_temp = PI2_p
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_tp1_tp2
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          ENDIF
                          If(f1234r.lt.zero)then
                              p2_temp = PI1_p
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_tp1_tp2
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          ENDIF
                      Endif  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        else  !photon has probability to fall into black hole.
                            If(f1234r.le.zero)then
                                index_p4(1)=0
                                cases_int=1
                                call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,&
                                                   index_p4,rff_p,integ04,cases_int)
                                PI0_obs_hori=integ04(1)
                                If(p.lt.PI0_obs_hori)then   
! equation (41) in Yang & Wang (2012).     
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)        
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p                                  
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0
                                sign_pr = -one 
                            ELSE  !p_r>0, photon will meet the r_tp2 &
                                  !turning point and turn around then goto vevnt horizon.     
                                index_p4(1)=0
                                cases_int=1        
                                call weierstrass_int_J3(tp2,tobs,dd,del,a4,b4,&
                                              index_p4,rff_p,integ04,cases_int)  
                                call weierstrass_int_j3(tp2,thorizon,dd,del,a4,b4,&
                                              index_p4,rff_p,integ14,cases_int)
                                PI0_obs_tp2=integ04(1)        
                                PI2_p=PI0_obs_tp2
                                PI0_total=integ14(1)+PI0_obs_tp2
                                If(p.le.PI0_obs_tp2)then
                                    t1=0
                                    t2=0
                                    pp=-p 
! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    sign_pr = one
                                else
                                    t1=0
                                    t2=1
                                    If(p.lt.PI0_total)then  
! equation (41) in Yang & Wang (2012).      
                                        tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=p-two*PI0_obs_tp2
                                        p2=p-PI0_obs_tp2
                                    else        
                                        tp=thorizon !Fall into black hole. 
                                        r_coord = rhorizon
                                        pp=PI0_total-two*PI0_obs_tp2
                                        p2=PI0_total-PI0_obs_tp2
                                    endif   
                                    sign_pr = -one     
                                endif        
                            ENDIF
                        ENDIF                             
                    END SELECT  
              !****************************************************************** 
                    index_p4(1)=-1
                    index_p4(2)=-2
                    index_p4(3)=0
                    index_p4(4)=-4
                !pp part ***************************************************        
                    cases_int=4
                    call weierstrass_int_J3(tobs,tp,dd,del,h,b4,&
                              index_p4,dabs(pp),integ4,cases_int)
! equation (62) in Yang & Wang (2012).
                    pp_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pp*r_tp1**two  
! equation (63) in Yang & Wang (2012).
                    pp_time=integ4(2)*b0/two+pp*(two*r_tp1+four+Ap*wp)+pp_aff
                    time_temp=pp*(-Am*wm)           
        
                    cases_int=2        
                    call weierstrass_int_J3(tobs,tp,dd,del,-E_add,&
                             b4,index_p4,dabs(pp),integ4,cases_int)
! equation (63) in Yang & Wang (2012).
                    pp_time=pp_time-Ap*wbarp*integ4(2)         
                    IF(a_spin.NE.zero)THEN
                        call weierstrass_int_J3(tobs,tp,dd,del,-E_m,&
                              b4,index_p4,dabs(pp),integ14,cases_int)
! equation (63) in Yang & Wang (2012).
                        pp_time=pp_time+Am*wbarm*integ14(2)+time_temp     
! equation (65) in Yang & Wang (2012).               
                        pp_phi=pp*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                 -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                    ELSE
                        pp_phi=zero
                    ENDIF                
                    IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).
                        pp_phi=pp_phi+pp*lambda
                    ENDIF        
                    !p1 part *******************************************************
                    IF(t1 .EQ. 0)THEN
                        p1_phi=ZERO
                        p1_time=ZERO
                        p1_aff=ZERO
                    ELSE
                        IF(PI1_aff .EQ. zero .AND. PI1_time .EQ. zero)THEN
                            cases_int=4
                            call weierstrass_int_J3(tobs,infinity,dd,del&
                                     ,h,b4,index_p4,PI0,integ4,cases_int)
! equation (62) in Yang & Wang (2012).
                            PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*&
                                           b0*r_tp1/two+PI0*r_tp1**two    
! equation (63) in Yang & Wang (2012).     
                            PI1_time=integ4(2)*b0/two+PI0*(two*r_tp1+four+Ap*wp)+PI1_aff        
                            time_temp=PI0*(-Am*wm)
        
                            cases_int=2        
                            call weierstrass_int_J3(tobs,infinity,dd,del,-E_add,&
                                                b4,index_p4,PI0,integ4,cases_int)
! equation (63) in Yang & Wang (2012).
                            PI1_time=PI1_time-Ap*wbarp*integ4(2) 
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tobs,infinity,dd,del,-E_m,&
                                                 b4,index_p4,PI0,integ14,cases_int)
! equation (63) in Yang & Wang (2012).       
                                PI1_time=PI1_time+Am*wbarm*integ14(2)+time_temp   
! equation (65) in Yang & Wang (2012).
                                PI1_phi=PI0*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                           -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI1_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).       
                                PI1_phi=PI1_phi+PI0*lambda
                            ENDIF
                        ENDIF 
! equation (55) in Yang & Wang (2012).       
                        p1_aff=PI1_aff-pp_aff
                        p1_time=PI1_time-pp_time
                        P1_phi=PI1_phi-pp_phi 
                    ENDIF
                !p2 part *******************************************************
                    IF(t2.EQ.ZERO)THEN
                        p2_phi=ZERO
                        p2_time=ZERO
                        p2_aff=ZERO
                    ELSE
                        IF(PI2_aff .EQ. zero .AND. PI2_time .EQ. zero)THEN
                            cases_int=4
                            call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,&
                                          index_p4,PI2_p,integ4,cases_int)
! equation (62) in Yang & Wang (2012).       
                            PI2_aff=integ4(4)*b0**two/sixteen+integ4(2)*&
                                            b0*r_tp1/two+PI2_p*r_tp1**two  
! equation (63) in Yang & Wang (2012).       
                            PI2_time=integ4(2)*b0/two+PI2_p*(two*r_tp1+four+Ap*wp)+PI2_aff        
                            time_temp=PI2_p*(-Am*wm)

                            cases_int=2        
                            call weierstrass_int_J3(tp2,tobs,dd,del,-E_add,&
                                         b4,index_p4,PI2_p,integ4,cases_int)
! equation (63) in Yang & Wang (2012).
                            PI2_time=PI2_time-Ap*wbarp*integ4(2)!+Am*wbarm*integ14(2) 
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tp2,tobs,dd,del,-E_m,&
                                          b4,index_p4,PI2_p,integ14,cases_int)
! equation (63) in Yang & Wang (2012).       
                                PI2_time=PI2_time+Am*wbarm*integ14(2)+time_temp   
! equation (65) in Yang & Wang (2012).
                                PI2_phi=PI2_p*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                        -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI2_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).       
                                PI2_phi=PI2_phi+PI2_p*lambda
                            ENDIF
                        ENDIF
! equation (55) in Yang & Wang (2012).       
                        p2_aff=PI2_aff+pp_aff
                        p2_time=PI2_time+pp_time
                        p2_phi=PI2_phi+pp_phi
                    ENDIF
                    !phi, aff,time part *******************************************************
! equation (56) in Yang & Wang (2012).       
                    If(f1234r.ne.zero)then
                        phyr=dsign(one,-f1234r)*pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timer=dsign(one,-f1234r)*pp_time+two*(t1*p1_time+t2*p2_time)
                        affr=dsign(one,-f1234r)*pp_aff+two*(t1*p1_aff+t2*p2_aff)
                    else
                        If(robs.eq.r_tp1)then
                            phyr=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=-pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=-pp_aff+two*(t1*p1_aff+t2*p2_aff)        
                        else
                            phyr=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=pp_aff+two*(t1*p1_aff+t2*p2_aff)
                        endif        
                    endif 
!************************************************************************************************ 
                    If(a_spin.eq.zero)then
                        If(cc.eq.zero)then
                            If(f1234r.lt.zero)then
                                If(p.lt.one/rhorizon-one/robs)then
                                    radius=robs/(robs*p+one)
                                else
                                    radius=rhorizon                  
                                endif
                                sign_pr = -one
                            else
                                If(p.lt.one/robs)then
                                    radius=robs/(one-robs*p)
                                else
                                    radius=infinity          
                                endif  
                                sign_pr = one
                            endif        
                        endif
                        If(cc.eq.-27.D0)then
                                sqt3=dsqrt(three)        
                            If(f1234r.lt.zero)then
                                cr=-three*dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/&
                                                sqt3)/(three-robs))*dexp(three*sqt3*p)-sqt3
                                dr=-dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                                (three-robs))*dexp(three*sqt3*p)+two/sqt3
                                If(p.ne.zero)then        
                                    radius=(three+cr*dr+dsqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=robs!infinity
                                endif
                                sign_pr = - one
                            else        
                                cr=-three*dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                                   (three-robs))*dexp(-three*sqt3*p)-sqt3
                                dr=-dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                          (three-robs))*dexp(-three*sqt3*p)+two/sqt3
                                PI0=dLog(dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                         (robs-three)))/three/sqt3-dLog(one+two/sqt3)/three/sqt3
                                If(p.lt.PI0)then        
                                    radius=(three+cr*dr+dsqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=infinity
                                endif 
                                sign_pr = one                       
                            endif                
                        endif
                    endif                                   
                ELSE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
! equation (44) in Yang & Wang (2012). !equation R(r)=0 has no real roots. we use the Legendre elliptic 
                    u=real(bb(4))      !integrations and functions to compute the calculations.
                    w=dabs(aimag(bb(4)))
                    v=dabs(aimag(bb(2)))
                    If(u.ne.zero)then
! equation (45) in Yang & Wang (2012).       
                        L1=(four*u**2+w**2+v**2+dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                        L2=(four*u**2+w**2+v**2-dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012).       
                        thorizon=dsqrt((L1-one)/(L1-L2))*(rhorizon-u*&
                                 (L1+one)/(L1-one))/dsqrt((rhorizon-u)**2+w**2)
! equation (48) in Yang & Wang (2012).   
                        m2=(L1-L2)/L1
                        tinf=dsqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/&
                             dsqrt((robs-u)**two+w**two)
                        t_inf=dsqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012).       
                        pinf=EllipticF(tinf,m2)/(w*dsqrt(L1))
                        call sncndn(p*w*dsqrt(L1)+dsign(one,f1234r)*&
                                 pinf*w*dsqrt(L1),one-m2,sn,cn,dn)
                        f1=u**2+w**2
                        g1=-two*u
                        h1=one
                        f2=u**2+v**2
                        g2=-g1
                        h2=one
                        a5=zero
                        b5=one
                        index_p5(1)=-1
                        index_p5(2)=-2
                        index_p5(3)=2
                        index_p5(4)=-4
                        index_p5(5)=4
                        IF(f1234r.lt.zero)THEN
                            PI0=pinf-EllipticF(thorizon,m2)/(w*dsqrt(L1))                                
                            if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012).       
                                y=u+(-two*u+w*(L1-L2)*sn*dabs(cn))/((L1-L2)*sn**2-(L1-one))
                                r_coord = y 
                                pp=p  
                            else
                                y=rhorizon
                                r_coord = y
                                pp=PI0
                            endif         
                            x=robs 
                            sign_pr = - one
                        ELSE
                            PI0=EllipticF(t_inf,m2)/(w*dsqrt(L1))-pinf
                            if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012).       
                                x=u+(-two*u-w*(L1-L2)*sn*dabs(cn))/((L1-L2)*sn**2-(L1-one))
                                r_coord = x
                                pp=p
                            else
                                x=infinity 
                                r_coord = x
                                pp=PI0
                            endif
                            y=robs   
                            sign_pr = one      
                        ENDIF 
 !affine parameter part integration **********************************************
                        cases_int=5    
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,&
                             h2,a5,b5,index_p5,dabs(pp),integ5,cases_int)
! equation (67) in Yang & Wang (2012).   
                        affr=integ5(5)
! equation (68) in Yang & Wang (2012).   
                        timer=two*integ5(3)+four*pp+affr
                        cases_int=2    
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,&
                             h2,-r_add,b5,index_p5,dabs(pp),integ5,cases_int)
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,&
                             h2,-r_m,b5,index_p5,dabs(pp),integ15,cases_int)
                        !phy part**************************************************************************
! equation (68) in Yang & Wang (2012).   
                        timer=timer+Ap*integ5(2)-Am*integ15(2)
! equation (69) in Yang & Wang (2012).   
                        phyr=a_spin*(B_add*integ5(2)-B_m*integ15(2))        
                        IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).   
                            phyr=phyr+pp*lambda
                        ENDIF 
                    ELSE
                        If(f1234r.lt.zero)then
                            PI0=(datan(robs/w)-datan(rhorizon/w))/w        
                            If(p.lt.PI0)then
                                radius=w*dtan(datan(robs/w)-p*w)
                                r_coord = radius 
                            else
                                radius=rhorizon
                                r_coord = radius
                            endif
                            !timer part ****************************************    
                            y=radius
                            x=robs
                            sign_pr = - one
                        ELSE
                            PI0=(PI/two-datan(robs/w))/w        
                            If(p.lt.PI0)then
                                radius=w*dtan(datan(robs/w)+p*w)
                                r_coord = radius
                            else
                                radius=infinity
                                r_coord = radius
                            endif
                            !timer part ***************************************
                            y=robs
                            x=radius
                            sign_pr = one
                        ENDIF
                        pp_time=(x-y)+datan(x/w)*(-w+four/w-r_add*Ap/w/&
                                  (w**two+r_add**two)+r_m*Am/w/(w**two+r_m**two))-&
                                  datan(y/w)*(-w+four/w-r_add*Ap/w/(w**two+&
                                  r_add**two)+r_m*Am/w/(w**two+r_m**two))
                        pp_time=pp_time+dlog(x**two+w**two)*(one-Ap/two/&
                                (w**two+r_add**two)+Am/two/(w**two+r_m**two))-&
                                (dlog(y**two+w**two)*(one-Ap/two/(w**two+r_add**two)+&
                                Am/two/(w**two+r_m**two)))        
                        timer=pp_time+Ap*dlog(dabs(x-r_add))/(w**two+r_add**two)-&
                              Am*dlog(dabs(x-r_m))/(w**two+r_m**two)-&
                              (Ap*dlog(dabs(y-r_add))/(w**two+r_add**two)-&
                              Am*dlog(dabs(y-r_m))/(w**two+r_m**two))        
                        !affine parameter part **************************************
                        affr=(x-y)-w*datan(x/w)+w*datan(y/w)          
                        !phy part ***************************************************
                        IF(a_spin .NE. zero)THEN
                            phyr=(-B_add*r_add/w/(r_add**two+w**two)+B_m*r_m/w/(r_m**two+w**two))&
                                               *(datan(x/w)-datan(y/w))+&
                                 dlog(dabs(x-r_add)/dsqrt(x**two+w**two))*B_add/(r_add*two+w**two)-&
                                 dlog(dabs(y-r_add)/dsqrt(y**two+w**two))*B_add/(r_add*two+w**two)-&
                                 dlog(dabs(x-r_m)/dsqrt(x**two+w**two))*B_m/(r_m*two+w**two)+&
                                 dlog(dabs(y-r_m)/dsqrt(y**two+w**two))*B_m/(r_m*two+w**two)
                            phyr=phyr*a_spin 
                        ELSE
                            phyr=zero
                        ENDIF        
                        If(muobs.eq.zero.and.f1234t.eq.zero)then         
                            phyr=phyr+lambda*(datan(x/w)-datan(y/w))/w
                        ENDIF   
                    ENDIF                        
                ENDIF
                count_num=count_num+1
        !*****************************************************************************************
        else 
        !*****************************************************************************************
            If(f1234r.eq.f1234r_1.and.f1234t.eq.f1234t_1.and.&
               lambda.eq.lambda_1.and.q.eq.q_1.and.&
               sinobs.eq.sinobs_1.and.muobs.eq.muobs_1.and.a_spin.eq.&
               a_spin_1.and.robs.eq.robs_1.and.scal.eq.scal_1)then
                !***********************************************************************
                If(reals.ne.0)then  !** R(r)=0 has real roots and turning points exists in radial r.          
                    select case(cases)
                    CASE(1)
                        If(f1234r .ge. zero)then !**photon will goto infinity. 
                            If(p.lt.PI0_obs_inf)then  
! equation (41) in Yang & Wang (2012).             
                                tp=weierstrassP(p+PI0,g2,g3,dd,del)        
                                r_coord = r_tp1+b0/(four*tp-b1)
                                pp=-p                                  
                            else
                                tp=tinf! !Goto infinity, far away. 
                                r_coord = infinity
                                pp=-PI0_obs_inf 
                            endif
                            t1=0
                            t2=0
                            sign_pr = one
                        ELSE 
                            If(.not.indrhorizon)then 
                                t2=0
                                If(p.le.PI0)then
                                    t1=0
                                    pp=p  
! equation (41) in Yang & Wang (2012).       
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    sign_pr = -one
                                else
                                    t1=1        
                                    If(p.lt.PI0_total)then   
! equation (41) in Yang & Wang (2012).            
                                        tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=two*PI0-p
                                        p1=dabs(p-PI0)
                                    else        
                                        tp=tinf !Goto infinity, far away.
                                        r_coord = infinity
                                        pp=-PI0_total+two*PI0 
                                        p1=pI0_total-PI0
                                    endif
                                    sign_pr = one
                                endif        
                            ELSE     
!f1234r<0, photon will fall into black hole unless something encountered.
                                If(p.lt.PI0_obs_hori)then  
! equation (41) in Yang & Wang (2012).             
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)        
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p                                  
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0
                                sign_pr = -one
                            ENDIF
                        ENDIF         
                    CASE(2)
                        If(.not.indrhorizon)then
                            If(f1234r.lt.zero)then
                                PI01=-PI0
                            else
                                PI01=PI0        
                            endif
! equation (41) in Yang & Wang (2012).       
                            tp=weierstrassP(p+PI01,g2,g3,dd,del)
                            r_coord = r_tp1+b0/(four*tp-b1)
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! this part of the code aims on obtaining the sign of the r component of  !|
               ! 4-momentum of the photon.                                               !|
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    If (f1234r > zero) then                                              !|
                        up = dmod(p + PI01, periodwp)                                    !|
                        If ( up <= half_periodwp ) then                                  !|
                            sign_pr = one                                                !|
                        else                                                             !|
                            sign_pr = - one                                              !|
                        endif                                                            !|
                    Else If (f1234r < zero) then                                         !|
                        up = dmod(p + PI01 + half_periodwp, periodwp)                    !|
                        If (up < half_periodwp) then                                     !|
                            sign_pr = -one                                               !|
                        Else                                                             !|
                            sign_pr = one                                                !|
                        endif                                                            !|
                    Else                                                                 !|
                        If (robs == r_tp1) then                                          !|
                            If ( dmod(p, periodwp) <= half_periodwp ) then               !|
                                sign_pr = one                                            !|
                            else                                                         !|
                                sign_pr = -one                                           !|
                            endif                                                        !|
                        else                                                             !| 
                            If ( dmod(p + half_periodwp, periodwp) &                     !|
                                                   <= half_periodwp ) then               !|
                                sign_pr = one                                            !|
                            else                                                         !|
                                sign_pr = -one                                           !|
                            endif                                                        !| 
                        endif                                                            !|
                    Endif                                                                !|
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                            index_p4(1)=0
                            cases_int=1        
                            call weierstrass_int_J3(tobs,tp,dd,del,a4,&
                                    b4,index_p4,rff_p,integ4,cases_int) 
! equation (57) in Yang & Wang (2012).       
                            pp=integ4(1) 
                            p1=PI0-pp
                            p2=p_tp1_tp2-p1  
                        !******************************************************
! equation (58) in Yang & Wang (2012).       
                            Do j=0,100
                                Do i=j,j+1
                                    If(robs_eq_rtp)then        
                                        If(robs.eq.r_tp1)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        else
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    else
                                        If(f1234r.gt.zero)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        endif
                                        If(f1234r.lt.zero)then
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    endif  
                                    If(dabs(p-p_temp).lt.1.D-4)goto 210
                                Enddo
                            Enddo
                        !**********************************************************
                            210  continue
                            tt1 = t1
                            tt2 = t2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! To determine t1 and t2, which are the times of the photon meets the two
! turnning points: r_tp1 and r_tp2 respectively.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      If(robs_eq_rtp)then
                          p1_temp = zero
                          p2_temp = p_tp1_tp2
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(robs.eq.r_tp1)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_tp1_tp2
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          Else If(robs.eq.r_tp2)then
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_tp1_tp2
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO
                          Endif 
                      Else
                          p1_temp = zero
                          t1 = 0
                          t2 = 0
                          N_temp = 0
                          If(f1234r.gt.zero)then
                              p2_temp = PI2_p
                              Do While(.TRUE.)
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_tp1_tp2
                                      t1 = t2
                                      t2 = N_temp-t1
                                  ENDIF
                              ENDDO
                          ENDIF
                          If(f1234r.lt.zero)then 
                              p2_temp = PI1_p
                              Do While(.TRUE.)   
                                  IF(p1_temp.le.p .and. p.le.p2_temp)then 
                                      EXIT
                                  ELSE
                                      N_temp = N_temp + 1
                                      p1_temp = p2_temp
                                      p2_temp = p2_temp + p_tp1_tp2
                                      t2 = t1
                                      t1 = N_temp-t2
                                  ENDIF
                              ENDDO 
                          ENDIF
                      Endif   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                  
                        else  !photon has probability to fall into black hole.
                            If(f1234r.le.zero)then 
                                If(p.lt.PI0_obs_hori)then    
! equation (41) in Yang & Wang (2012).        
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)        
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p                                  
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0
                                sign_pr = -one
                            ELSE  !p_r>0, photon will meet the r_tp2 turning &
                                  !point and turn around then goto vevnt horizon.  
                                If(p.le.PI0_obs_tp2)then
                                    t1=0
                                    t2=0
                                    pp=-p 
! equation (41) in Yang & Wang (2012).    
                                    tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    sign_pr = one
                                else
                                    t1=0
                                    t2=1
                                    If(p.lt.PI0_total)then    
! equation (41) in Yang & Wang (2012).        
                                        tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=p-two*PI0_obs_tp2
                                        p2=p-PI0_obs_tp2
                                    else        
                                        tp=thorizon !Fall into black hole. 
                                        r_coord = rhorizon
                                        pp=PI0_total-two*PI0_obs_tp2
                                        p2=PI0_total-PI0_obs_tp2
                                    endif
                                    sign_pr = -one
                                endif        
                            ENDIF
                        ENDIF                             
                    END SELECT 
              !****************************************************************** 
                    index_p4(1)=-1
                    index_p4(2)=-2
                    index_p4(3)=0
                    index_p4(4)=-4
                !pp part ***************************************************        
                    cases_int=4
                    call weierstrass_int_J3(tobs,tp,dd,del,h,&
                           b4,index_p4,dabs(pp),integ4,cases_int)
! equation (62) in Yang & Wang (2012).    
                    pp_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pp*r_tp1**two  
! equation (63) in Yang & Wang (2012).           
                    pp_time=integ4(2)*b0/two+pp*(two*r_tp1+four+Ap*wp)+pp_aff
                    time_temp=pp*(-Am*wm)         
        
                    cases_int=2        
                    call weierstrass_int_J3(tobs,tp,dd,del,-E_add,b4,&
                                   index_p4,dabs(pp),integ4,cases_int)
! equation (63) in Yang & Wang (2012).    
                    pp_time=pp_time-Ap*wbarp*integ4(2)         
                    IF(a_spin.NE.zero)THEN
                        call weierstrass_int_J3(tobs,tp,dd,del,-E_m,&
                               b4,index_p4,dabs(pp),integ14,cases_int)
! equation (63) in Yang & Wang (2012).    
                        pp_time=pp_time+Am*wbarm*integ14(2)+time_temp     
! equation (65) in Yang & Wang (2012).                   
                        pp_phi=pp*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                    -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                    ELSE
                        pp_phi=zero
                    ENDIF                
                    IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).    
                        pp_phi=pp_phi+pp*lambda
                    ENDIF        
                 !p1 part *******************************************************
                    IF(t1 .EQ. 0)THEN
                        p1_phi=ZERO
                        p1_time=ZERO
                        p1_aff=ZERO
                    ELSE
                        IF(PI1_aff .EQ. zero .AND. PI1_time .EQ. zero)THEN
                            cases_int=4
                            call weierstrass_int_J3(tobs,infinity,dd,del,&
                                       h,b4,index_p4,PI0,integ4,cases_int)
! equation (62) in Yang & Wang (2012).    
                            PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*&
                                               b0*r_tp1/two+p1*r_tp1**two   
! equation (63) in Yang & Wang (2012).          
                            PI1_time=integ4(2)*b0/two+PI0*(two*r_tp1+four+Ap*wp)+PI1_aff        
                            time_temp=PI0*(-Am*wm)
        
                            cases_int=2        
                            call weierstrass_int_J3(tobs,infinity,dd,del,-E_add,&
                                                b4,index_p4,PI0,integ4,cases_int)
! equation (63) in Yang & Wang (2012).    
                            PI1_time=PI1_time-Ap*wbarp*integ4(2) 
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tobs,infinity,dd,del,-E_m,&
                                                 b4,index_p4,PI0,integ14,cases_int)
! equation (65) in Yang & Wang (2012).    
                                PI1_time=PI1_time+Am*wbarm*integ14(2)+time_temp   
                                PI1_phi=PI0*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                           -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI1_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).    
                                PI1_phi=PI1_phi+PI0*lambda
                            ENDIF
                        ENDIF 
! equation (55) in Yang & Wang (2012).    
                        p1_aff=PI1_aff-pp_aff
                        p1_time=PI1_time-pp_time
                        P1_phi=PI1_phi-pp_phi 
                    ENDIF
                !p2 part *******************************************************
                    IF(t2.EQ.ZERO)THEN
                        p2_phi=ZERO
                        p2_time=ZERO
                        p2_aff=ZERO
                    ELSE
                        IF(PI2_aff .EQ. zero .AND. PI2_time .EQ. zero)THEN
                            cases_int=4
                            call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,&
                                          index_p4,PI2_p,integ4,cases_int)
! equation (62) in Yang & Wang (2012).    
                            PI2_aff=integ4(4)*b0**two/sixteen+integ4(2)*&
                                            b0*r_tp1/two+PI2_p*r_tp1**two   
! equation (63) in Yang & Wang (2012).          
                            PI2_time=integ4(2)*b0/two+PI2_p*(two*r_tp1+four+Ap*wp)+PI2_aff        
                            time_temp=PI2_p*(-Am*wm)

                            cases_int=2        
                            call weierstrass_int_J3(tp2,tobs,dd,del,-E_add,&
                                         b4,index_p4,PI2_p,integ4,cases_int)
! equation (63) in Yang & Wang (2012).    
                            PI2_time=PI2_time-Ap*wbarp*integ4(2)!+Am*wbarm*integ14(2) 
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tp2,tobs,dd,del,-E_m,&
                                          b4,index_p4,PI2_p,integ14,cases_int)
! equation (63) in Yang & Wang (2012).    
                                PI2_time=PI2_time+Am*wbarm*integ14(2)+time_temp   
! equation (65) in Yang & Wang (2012).    
                                PI2_phi=PI2_p*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                           -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI2_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).    
                                PI2_phi=PI2_phi+PI2_p*lambda
                            ENDIF
                        ENDIF
! equation (55) in Yang & Wang (2012).    
                        p2_aff=PI2_aff+pp_aff
                        p2_time=PI2_time+pp_time
                        p2_phi=PI2_phi+pp_phi
                    ENDIF                    
                    !phi, aff,time part *******************************************************
! equation (56) in Yang & Wang (2012).    
                    If(f1234r.ne.zero)then
                        phyr=dsign(one,-f1234r)*pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timer=dsign(one,-f1234r)*pp_time+two*(t1*p1_time+t2*p2_time)
                        affr=dsign(one,-f1234r)*pp_aff+two*(t1*p1_aff+t2*p2_aff)
                    else
                        If(robs.eq.r_tp1)then
                            phyr=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=-pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=-pp_aff+two*(t1*p1_aff+t2*p2_aff)        
                        else
                            phyr=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=pp_aff+two*(t1*p1_aff+t2*p2_aff)
                        endif        
                    endif 
!************************************************************************************************ 
                    If(a_spin.eq.zero)then
                        If(cc.eq.zero)then
                            If(f1234r.lt.zero)then
                                If(p.lt.one/rhorizon-one/robs)then
                                    radius=robs/(robs*p+one)
                                else
                                    radius=rhorizon                  
                                endif
                                sign_pr = -one
                            else
                                If(p.lt.one/robs)then
                                    radius=robs/(one-robs*p)
                                else
                                    radius=infinity          
                                endif 
                                sign_pr = one
                            endif        
                        endif
                        If(cc.eq.-27.D0)then
                                sqt3=dsqrt(three)        
                            If(f1234r.lt.zero)then
                                cr=-three*dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/&
                                              sqt3)/(three-robs))*dexp(three*sqt3*p)-sqt3
                                dr=-dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                              (three-robs))*dexp(three*sqt3*p)+two/sqt3
                                If(p.ne.zero)then        
                                    radius=(three+cr*dr+dsqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=robs!infinity
                                endif
                                sign_pr = -one
                            else        
                                cr=-three*dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/&
                                              sqt3)/(three-robs))*dexp(-three*sqt3*p)-sqt3
                                dr=-dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/&
                                                (three-robs))*dexp(-three*sqt3*p)+two/sqt3
                                PI0=dLog(dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)&
                                          /(robs-three)))/three/sqt3-dLog(one+two/sqt3)/three/sqt3
                                If(p.lt.PI0)then        
                                    radius=(three+cr*dr+dsqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=infinity
                                endif
                                sign_pr = one
                            endif                
                        endif
                    endif                                   
                ELSE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! equation R(r)=0 has no real roots. we use the Legendre elliptic 
! integrations and functions to compute the calculations. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    If(u.ne.zero)then 
                        call sncndn(p*w*dsqrt(L1)+dsign(one,f1234r)*&
                                    pinf*w*dsqrt(L1),one-m2,sn,cn,dn) 
                        index_p5(1)=-1
                        index_p5(2)=-2
                        index_p5(3)=2
                        index_p5(4)=-4
                        index_p5(5)=4
                        IF(f1234r.lt.zero)THEN                                 
                            if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012).    
                                y=u+(-two*u+w*(L1-L2)*sn*dabs(cn))/&
                                             ((L1-L2)*sn**2-(L1-one))
                                r_coord = y 
                                pp=p  
                            else
                                y=rhorizon
                                r_coord = y
                                pp=PI0
                            endif         
                            x=robs
                            sign_pr = -one
                        ELSE 
                            if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012).    
                                x=u+(-two*u-w*(L1-L2)*sn*dabs(cn))/&
                                             ((L1-L2)*sn**2-(L1-one))
                                r_coord = x
                                pp=p
                            else
                                x=infinity 
                                r_coord = x
                                pp=PI0
                            endif
                            y=robs
                            sign_pr = one
                        ENDIF 
!affine parameter part integration **********************************************
                        cases_int=5
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,a5,&
                                        b5,index_p5,dabs(pp),integ5,cases_int)
! equation (67) in Yang & Wang (2012).    
                        affr=integ5(5)
! equation (68) in Yang & Wang (2012).    
                        timer=two*integ5(3)+four*pp+affr
                        cases_int=2
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,&
                                -r_add,b5,index_p5,dabs(pp),integ5,cases_int)
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,&
                                -r_m,b5,index_p5,dabs(pp),integ15,cases_int)
!phy part**************************************************************************
! equation (68) in Yang & Wang (2012).    
                        timer=timer+Ap*integ5(2)-Am*integ15(2)
! equation (69) in Yang & Wang (2012).    
                        phyr=a_spin*(B_add*integ5(2)-B_m*integ15(2))        
                        IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).    
                            phyr=phyr+pp*lambda        
                        ENDIF 
                    ELSE
                        If(f1234r.lt.zero)then         
                            If(p.lt.PI0)then
                                radius=w*dtan(datan(robs/w)-p*w)
                                r_coord = radius 
                            else
                                radius=rhorizon
                                r_coord = radius
                            endif
!timer part ******************************************************************        
                            y=radius
                            x=robs
                            sign_pr = -one
                        ELSE 
                            If(p.lt.PI0)then
                                radius=w*dtan(datan(robs/w)+p*w)
                                r_coord = radius
                            else
                                radius=infinity
                                r_coord = radius
                            endif
!timer part ************************************************************************
                            y=robs
                            x=radius
                            sign_pr = one
                        ENDIF
                        pp_time=(x-y)+datan(x/w)*(-w+four/w-r_add*Ap/w/(w**two+&
                                  r_add**two)+r_m*Am/w/(w**two+r_m**two))-&
                                  datan(y/w)*(-w+four/w-r_add*Ap/w/(w**two+&
                                  r_add**two)+r_m*Am/w/(w**two+r_m**two))
                        pp_time=pp_time+dlog(x**two+w**two)*(one-Ap/two/(w**two+&
                                  r_add**two)+Am/two/(w**two+r_m**two))-&
                                  (dlog(y**two+w**two)*(one-Ap/two/(w**two+&
                                  r_add**two)+Am/two/(w**two+r_m**two)))        
                        timer=pp_time+Ap*dlog(dabs(x-r_add))/(w**two+&
                                   r_add**two)-Am*dlog(dabs(x-r_m))/(w**two+r_m**two)-&
                                   (Ap*dlog(dabs(y-r_add))/(w**two+r_add**two)-&
                                   Am*dlog(dabs(y-r_m))/(w**two+r_m**two))        
!affine parameter part *****************************************************************
                        affr=(x-y)-w*datan(x/w)+w*datan(y/w)          
!phy part ******************************************************************************
                        IF(a_spin .NE. zero)THEN
                            phyr=(-B_add*r_add/w/(r_add**two+w**two)+B_m*r_m/w/(r_m**two+w**two))&
                                               *(datan(x/w)-datan(y/w))+&
                                 dlog(dabs(x-r_add)/dsqrt(x**two+w**two))*B_add/(r_add*two+w**two)-&
                                 dlog(dabs(y-r_add)/dsqrt(y**two+w**two))*B_add/(r_add*two+w**two)-&
                                 dlog(dabs(x-r_m)/dsqrt(x**two+w**two))*B_m/(r_m*two+w**two)+&
                                 dlog(dabs(y-r_m)/dsqrt(y**two+w**two))*B_m/(r_m*two+w**two)
                            phyr=phyr*a_spin 
                        ELSE
                            phyr=zero
                        ENDIF        
                        If(muobs.eq.zero.and.f1234t.eq.zero)then         
                            phyr=phyr+lambda*(datan(x/w)-datan(y/w))/w
                        ENDIF   
                    ENDIF                        
                ENDIF
            ELSE
                count_num=1
                goto 40
            endif
        endif                
        RETURN
     END SUBROUTINE INTRPART 

!********************************************************************************************
      Function  Pemdisk(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,mu,rout,rin) 
!********************************************************************************************
!*     PURPOSE:  Solves equation \mu(p)=mu, i.e. to search the value p_{em} of
!*               parameter p, corresponding to the intersection point of geodesic with with 
!*               disk, i.e. a surface has a constants inclination angle with respect to 
!*               equatorial plane of black hole.  
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are defined by equation 
!*                              (82)-(85) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon.
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               mu-------------mu=\cos(\pi/2-\theta_disk), where \theta_disk is the inclination 
!*                              angle of disk surface with respect to the equatorial plane of 
!*                              black hole.    
!*               rin, rout------inner and outer radius of disk.        
!*     OUTPUTS:  pemdisk--------value of root of equation \mu(p)= mu for p.  
!*                              pemdisk=-1.D0, if the photon goto infinity.
!*                              pemdisk=-2.D0, if the photon fall into event horizon.       
!*     REMARKS:                 This routine just search the intersection points of geodesic with 
!*                              up surface of disk. Following routine Pemdisk_all will searches 
!*                              intersection points of geodesic with up and down surface of disk.      
!*     ROUTINES CALLED: mutp, mu2p.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision Pemdisk,f1234(4),f3,f2,p,sinobs,muobs,a_spin,lambda,q,mu_tp,tposition,tp2,two,&
                         bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tinf,p1,p2,pp,mu,scal,zero,robs,&
                         mutemp,mu_tp2,delta,four,one,pm,rout,rin,re, sign_pr
        parameter(zero=0.D0,two=2.D0,four=4.D0,one=1.D0)
        integer  t1,t2,reals,i,j
        complex*16 bb(1:4)
        logical :: err,mobseqmtp 

        f3 = f1234(3)
        f2 = f1234(2) 
        mobseqmtp=.false. 
        call mutp(f2,f3,sinobs,muobs,a_spin,lambda,q,mu_tp,mu_tp2,reals,mobseqmtp) 
 
        If(reals.eq.2)then
            If(mobseqmtp)then        
                t1=0
                t2=0
            else
                If(muobs.gt.zero)then
                    If(f2.lt.zero)then        
                        t1=1
                        t2=0
                    endif
                    If(f2.gt.zero)then        
                        t1=0
                        t2=0                                
                    endif
                else
                    If(muobs.eq.zero)then
                        If(f2.lt.zero)then        
                            t1=1
                            t2=0
                        endif
                        If(f2.gt.zero)then        
                            t1=0
                            t2=1                                
                        endif                        
                    else
                        If(f2.lt.zero)then        
                            t1=0
                            t2=0
                        endif
                        If(f2.gt.zero)then        
                            t1=0
                            t2=1                                
                        endif
                    endif        
                endif        
            endif
            !write(*,*)'mu=',mu,B,t1,t2
            pm=mu2p(f3,f2,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,scal)  
            re = radius(pm,f1234(1),lambda,q,a_spin,robs,scal,sign_pr)  
            IF(re .le. rout .and. re .ge. rin)THEN
                Pemdisk = pm
                return 
            ENDIF
            IF(re .gt. rout)THEN
                Pemdisk = -one 
                RETURN
            ENDIF 
            IF(re .lt. rin)THEN
                Pemdisk = -two
                return   
            ENDIF
        else
            Pemdisk=-one                 
        endif 
        return
      End Function Pemdisk 
!********************************************************************* 
      Function  Pemdisk_all(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,mu,rout,rin) 
!********************************************************************* 
!*     PURPOSE:  Solves equation \mu(p)=mu, where \mu(p)=\mu(p), i.e. to search the value p_{em} of
!*               parameter p, corresponding to the intersection point of geodesic with 
!*               disk, i.e. a surface has a constants inclination angle with respect to 
!*               equatorial plane of black hole.
!* 
!*     INPUTS:   f1234(1:4)-----array of f_1, f_2, f_3, f_0, which are defined by equation 
!*                              (106)-(109) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               mu-------------mu=\cos(\pi/2-\theta_disk), where \theta_disk is the inclination 
!*                              angle of disk surface with respect to the equatorial plane of 
!*                              black hole.  
!*               rin, rout------inner and outer radius of disk.        
!*     OUTPUTS:  pemdisk_all--------value of root of equation \mu(p)= mu for p.  
!*                              pemdisk=-1.D0, if the photon goto infinity.
!*                              pemdisk=-2.D0, if the photon fall into event horizon.        
!*     REMARKS:                 This routine will searches intersection points of 
!*                              geodesic with double surfaces of disk.      
!*     ROUTINES CALLED: mutp, mu2p, radius.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision Pemdisk_all,f3,f2,p,sinobs,muobs,a_spin,lambda,q,mu_tp1,tposition,tp2,two,&
                         bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tinf,p1,p2,pp,mu,scal,robs,zero,&
                         mutemp,mu_tp2,delta,four,one,pm,f1234(4),rout,rin,re,rhorizon,sign_pr
        parameter(zero=0.D0,two=2.D0,four=4.D0,one=1.D0)
        integer  t1,t2,reals,i,j
        complex*16 bb(1:4)
        logical :: err,mobseqmtp 

        f3 = f1234(3)
        f2 = f1234(2) 
        mobseqmtp=.false.
        rhorizon = one + dsqrt(one-a_spin*a_spin) 
        call mutp(f2,f3,sinobs,muobs,a_spin,lambda,q,mu_tp1,mu_tp2,reals,mobseqmtp)         

        If(reals.eq.2)then
            Do j = 0,10
                Do i=j,j+1    
                    If(mobseqmtp)then        
                        If(muobs.eq.mu_tp1)then
                            t1=j
                            t2=i
                        else
                            t1=i
                            t2=j
                        endif
                    else
                        If(muobs.gt.zero)then
                            If(f2.lt.zero)then        
                                t1=i
                                t2=j
                            endif
                            If(f2.gt.zero)then        
                                t1=j
                                t2=i                                
                            endif
                        else
                            If(muobs.eq.zero)then
                                If(f2.lt.zero)then        
                                    t1=i
                                    t2=j
                                endif
                                If(f2.gt.zero)then        
                                    t1=j
                                    t2=i                                
                                endif                        
                            else
                                If(f2.lt.zero)then        
                                    t1=i
                                    t2=j
                                endif
                                If(f2.gt.zero)then        
                                    t1=j
                                    t2=i                                
                                endif
                            endif        
                        endif        
                    endif
                    pm=mu2p(f3,f2,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,scal) 
                    IF(pm .le. zero)cycle
                    re = radius(pm,f1234(1),lambda,q,a_spin,robs,scal,sign_pr)  
!write(*,*)'mu=',re,pm,t1,t2,mucos(pm,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal),f2
                    IF(re .le. rout .and. re .ge. rin)THEN
                        Pemdisk_all = pm
                        return
                    ELSE
                        IF(re .ge. infinity)THEN
                            Pemdisk_all = -one 
                            RETURN
                        ELSE
                            IF(re .le. rhorizon)THEN
                                Pemdisk_all = -two
                                return  
                            ENDIF
                        ENDIF 
                    ENDIF 
                ENDDO
            ENDDO
        else
            pm=-two                        
        endif
        Pemdisk_all=pm
        return
      End Function Pemdisk_all
!*****************************************************************************************************
      subroutine metricg(robs,sinobs,muobs,a_spin,somiga,expnu,exppsi,expmu1,expmu2)
!*****************************************************************************************************
!*     PURPOSE:  Computes Kerr metric, exp^\nu, exp^\psi, exp^mu1, exp^\mu2, and omiga at position:
!*               r_obs, \theta_{obs}.     
!*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).         
!*     OUTPUTS:  somiga,expnu,exppsi,expmu1,expmu2------------Kerr metrics under Boyer-Lindquist coordinates.
!*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision robs,theta,a_spin,Delta,bigA,two,sinobs,muobs,one,sigma
        Double precision ,optional :: somiga,expnu,exppsi,expmu1,expmu2        

        two=2.D0        
        one=1.D0
! equations (1) and (2) in Yang & Wang (2012).
        Delta=robs**two-two*robs+a_spin**two
        sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*Delta
        somiga=two*a_spin*robs/bigA
        expnu=dsqrt(sigma*Delta/bigA)
        exppsi=sinobs*dsqrt(bigA/sigma)
        expmu1=dsqrt(sigma/Delta)
        expmu2=dsqrt(sigma)
        return        
      End subroutine metricg        

!********************************************************************************************
      Subroutine lambdaq_old(alpha,beta,robs,sinobs,muobs,a_spin,scal,velocity,f1234,lambda,q)
!********************************************************************************************
!*     PURPOSE:  Computes constants of motion from impact parameters alpha and beta by using 
!*               formulae (110) and (112) in Yang & Wang (2012).    
!*     INPUTS:   alpha,beta-----Impact parameters.
!*               robs-----------radial coordinate of observer or the initial position of photon.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               velocity(1:3)--Array of physical velocity of observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  f1234(1:4)-----array of f_1, f_2, f_3, f_4, which was defined by equation 
!*                              (106)-(109) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision f1234(4),robs,sinobs,muobs,a_spin,lambda,q,As,Bs,A1,a2,b1,b2,c1,c2,d1,d2,&
                         Delta,zero,one,two,four,gff,Sigma,Ac,Bc,Cc,at,Bt,lambdap,lambdam,qp,qm,&
                         RaB2,scal,Vr,Vt,Vp,gama,gama_tp,gama_p,F1,F2,NN1,NN2,DD1,DD2,Aab,expnu2,&
                         eppsi2,epnu2,epmu12,epmu22,bigA,G1,G2,KF,three,alpha,beta,&
                         velocity(3),lambda2,somiga,q1,q2,q3
        parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)

        if(dabs(beta).lt.1.D-7)beta=zero
        if(dabs(alpha).lt.1.D-7)alpha=zero
! equations (113), (114) in Yang & Wang (2012).
        at=alpha/scal/robs
        Bt=beta/scal/robs        
        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
! equation (92) in Yang & Wang (2012).
        Aab=-one/dsqrt(one+at**two+Bt**two)
        gama=one/dsqrt(one-(Vr**two+Vt**two+Vp**two))
        gama_tp=one/dsqrt(one-(Vt**two+Vp**two))
        gama_p=one/dsqrt(one-Vp**two)
! equation (106)-(109) in Yang & Wang (2012).
        f1234(1)=(-gama*Vr+gama/gama_tp*Aab)
        f1234(2)=-((-gama*Vt+gama*gama_tp*Vr*Vt*Aab)*robs*scal+gama_tp/gama_p*beta*Aab)
        f1234(3)=-((-gama*Vp+gama*gama_tp*Vr*Vp*Aab)*robs*scal+gama_tp*&
                                  gama_p*Vt*Vp*beta*Aab+gama_p*alpha*Aab)
        f1234(4)=((gama-gama*gama_tp*Vr*Aab)*robs*scal-gama_tp*gama_p*Vt*&
                                  beta*Aab-gama_p*Vp*alpha*Aab)/robs/scal
        If(dabs(f1234(1)).lt.1.D-7)f1234(1)=zero
        If(dabs(f1234(2)).lt.1.D-7)f1234(2)=zero
        If(dabs(f1234(3)).lt.1.D-7)f1234(3)=zero
        If(dabs(f1234(4)).lt.1.D-7)f1234(4)=zero 
! equations (1), (2) in Yang & Wang (2012).
        Delta=robs**two-two*robs+a_spin**two
        Sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*Delta
        somiga=two*a_spin*robs/bigA
        expnu2=Sigma*Delta/bigA
        eppsi2=sinobs**two*bigA/Sigma
        epmu12=Sigma/Delta
        epmu22=Sigma 
! equations (110), (112) in Yang & Wang (2012).
        A1 = f1234(3)/(dsqrt(Delta)*Sigma/bigA*f1234(4)*robs*scal+f1234(3)*somiga*sinobs) 
        lambda = A1*sinobs 
        q=(A1*A1-a_spin*a_spin)*muobs*muobs+(f1234(2)/robs/scal/f1234(4)*&
                    (one-lambda*somiga))**two*bigA/Delta 
        !write(*,*)'kk=',f1234(3),dsqrt(Delta)*Sigma/bigA*f1234(4)*robs*scal,f1234(3)*somiga*sinobs
       return
       End subroutine lambdaq_old 

!********************************************************************************************
      Subroutine lambdaq(alpha,beta,robs,sinobs,muobs,a_spin,scal,velocity,f1234,lambda,q)
!********************************************************************************************
!*     PURPOSE:  Computes constants of motion from impact parameters alpha and beta by using 
!*               formulae (86) and (87) in Yang & Wang (2012).    
!*     INPUTS:   alpha,beta-----Impact parameters.
!*               robs-----------radial coordinate of observer or the initial position of photon.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               velocity(1:3)--Array of physical velocities of the observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are the components of 
!*                              four momentum of a photon measured under the LNRF frame, and 
!*                              defined by equations (82)-(85) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision f1234(4),robs,sinobs,muobs,a_spin,lambda,q,As,Bs,A1,a2,b1,b2,c1,c2,d1,d2,&
                         Delta,zero,one,two,four,gff,Sigma,Ac,Bc,Cc,at,Bt,lambdap,lambdam,qp,qm,&
                         RaB2,scal,Vr,Vt,Vp,gama,gama_tp,gama_p,F1,F2,NN1,NN2,DD1,DD2,Aab,expnu2,&
                         eppsi2,epnu2,epmu12,epmu22,bigA,G1,G2,KF,three,alpha,beta,&
                         velocity(3),lambda2,somiga,q1,q2,q3,prt,ptt,ppt
        parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)

        if(dabs(beta).lt.1.D-7)beta=zero
        if(dabs(alpha).lt.1.D-7)alpha=zero
! equations (94), (95) in Yang & Wang (2012).
        at=alpha/scal/robs
        Bt=beta/scal/robs        
        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
! equation (90) in Yang & Wang (2012). 
        gama=one/dsqrt(one-(Vr**two+Vt**two+Vp**two)) 

! equations (97), (98), (99) in Yang & Wang (2012). 
        prt=-one/dsqrt(one+at**two+Bt**two)
        ptt=Bt*prt
        ppt=at*prt
! equations (89), (90) and (91) in Yang & Wang (2012).
        f1234(1)=( gama*Vr-prt*(one+gama*gama*Vr*Vr/(one+gama))-&
                 ptt*gama*gama*Vr*Vt/(one+gama)-ppt*gama*gama*Vr*Vp/(one+gama) )*robs*scal 
        f1234(2)=( gama*Vt-prt*gama*gama*Vt*Vr/(one+gama)-ptt*(one+&
                 gama*gama*Vt*Vt/(one+gama))-ppt*gama*gama*Vt*Vp/(one+gama) )*robs*scal 
        f1234(3)=( gama*Vp-prt*gama*gama*Vp*Vr/(one+gama)-ptt*gama*gama*Vp*Vt/(one+gama)-&
                 ppt*(one+gama*gama*Vp*Vp/(one+gama)) )*robs*scal
        f1234(4)=gama*(one-prt*Vr-ptt*Vt-ppt*Vp) 

! Keep r component p_r of four momentum to be negative, so the photon will go
! to the central black hole.
        f1234(1)=-f1234(1) 
        If(dabs(f1234(1)).lt.1.D-6)f1234(1)=zero
        If(dabs(f1234(2)).lt.1.D-6)f1234(2)=zero
        If(dabs(f1234(3)).lt.1.D-6)f1234(3)=zero
        If(dabs(f1234(4)).lt.1.D-6)f1234(4)=zero 
! equations (1), (2) in Yang & Wang (2012).
        Delta=robs**two-two*robs+a_spin**two
        Sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*Delta
        somiga=two*a_spin*robs/bigA
        expnu2=Sigma*Delta/bigA
        eppsi2=sinobs**two*bigA/Sigma
        epmu12=Sigma/Delta
        epmu22=Sigma 
! equations (86) and (87) in Yang & Wang (2012).
        A1 = f1234(3)/(dsqrt(Delta)*Sigma/bigA*f1234(4)*robs*scal+f1234(3)*somiga*sinobs) 
        lambda = A1*sinobs 
        q=(A1*A1-a_spin*a_spin)*muobs*muobs+(f1234(2)/f1234(4)/robs/scal*&
                    (one-lambda*somiga))**two*bigA/Delta  
        !write(*,*)'kk=',f1234(3)!,dsqrt(Delta)*Sigma/bigA*f1234(4)*robs*scal,f1234(3)*somiga*sinobs
       return
       End subroutine lambdaq

!********************************************************************************************
      Subroutine initialdirection(pr,ptheta,pphi,sinobs,&
                              muobs,a_spin,robs,velocity,lambda,q,f1234)
!********************************************************************************************
!*     PURPOSE:  Computes constants of motion from components of initial 4 momentum 
!*               of photon measured by emitter in its local rest frame, by using 
!*               formulae (86) and (87) in Yang & Wang (2012).    
!*     INPUTS:   pr,ptheta,pphi-----components of initial 4 momentum of photon measured by 
!*               emitter in its local rest frame.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or the initial position of photon.
!*               velocity(1:3)--Array of physical velocities of the observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are the components of 
!*                              four momentum of a photon measured under the LNRF frame, and 
!*                              defined by equations (82)-(85) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision lambda,q,sinobs,muobs,a_spin,robs,zero,one,two,three,four,&
                        pt,pp,cost,sint,cosp,sinp,Ac,Bc,Cc,moment(4),velocity(3),Vr,Vt,Vp,&
                        a1,b1,c1,d1,a2,b2,c2,d2,gama,gama_tp,gama_p,f1234(4),KF,F1,F2,&
                        Delta,Er,Sigma,bigA,eppsi2,epmu12,epmu22,NN1,NN2,DD1,p_E,somiga,&
                        lambdap,lambdam,pr,ptheta,pphi,theta,phi,pr1,ptheta1,pphi1,&
                        lambda1,q1
        parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)
        optional pr,ptheta,pphi 
        integer cases

        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
! equation (92) in Yang & Wang (2012).
        gama=one/dsqrt(one-(Vr**two+Vt**two+Vp**two))
        gama_tp=one/dsqrt(one-(Vt**two+Vp**two))
        gama_p=one/dsqrt(one-Vp**two)
! equation (106)-(109) in Yang & Wang (2012).
        f1234(1)=-(-gama*Vr+gama/gama_tp*pr)
        f1234(2)=-(-gama*Vt+gama*gama_tp*Vr*Vt*pr+gama_tp/gama_p*ptheta)
        f1234(3)=-(-gama*Vp+gama*gama_tp*Vr*Vp*pr+gama_tp*gama_p*Vt*Vp*ptheta+gama_p*pphi)
        f1234(4)=(gama-gama*gama_tp*Vr*pr-gama_tp*gama_p*Vt*ptheta-gama_p*Vp*pphi)
        If(dabs(f1234(1)).lt.1.D-7)f1234(1)=zero
        If(dabs(f1234(2)).lt.1.D-7)f1234(2)=zero
        If(dabs(f1234(3)).lt.1.D-7)f1234(3)=zero
        If(dabs(f1234(4)).lt.1.D-7)f1234(4)=zero
! equations (1), (2) in Yang & Wang (2012).
        Delta=robs**two-two*robs+a_spin**two
        Sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*Delta
        somiga=two*a_spin*robs/bigA
        eppsi2=sinobs**two*bigA/Sigma
        epmu12=Sigma/Delta
        epmu22=Sigma
! equations (86) and (87) in Yang & Wang (2012).
        A1 = f1234(3)/(dsqrt(Delta)*Sigma/bigA*f1234(4)+f1234(3)*somiga*sinobs) 
        lambda = A1*sinobs 
        q=(A1*A1-a_spin*a_spin)*muobs*muobs+(f1234(2)/f1234(4)*&
                   (one-lambda*somiga))**two*bigA/Delta   
      return
      End subroutine initialdirection

!********************************************************************************************        
      Subroutine center_of_image(robs,sinobs,muobs,a_spin,scal,velocity,alphac,betac)
!********************************************************************************************
!*     PURPOSE:  Solves equations f_3(alphac,betac)=0, f_2(alphac,betac)=0, of (100)  
!*               and (101) in Yang & Wang (2012). alphac, betac are the coordinates of 
!*               center point of images on the screen of observer.    
!*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.  
!*               velocity(1:3)--Array of physical velocity of observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  alphac,betac---coordinates of center point of images on the screen of observer.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  9 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision robs,sinobs,muobs,a_spin,scal,velocity(3),alphac,betac,zero,one,two,four,&
                        Vr,Vt,Vp,gama,a1,b1,c1,alphap,alpham,betap,betam
        parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0)

        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
! equation (90) in Yang & Wang (2012).
        gama=one/dsqrt(one-(Vr*Vr+Vt*Vt+Vp*Vp))       

        If(Vt.ne.zero)then
            If(Vp.ne.zero)then   
                a1=(one+gama*gama*(Vt*Vt+Vp*Vp)/(one+gama))**two-gama*gama*(Vp*Vp+Vt*Vt)  
                b1=two*gama*gama*Vt*Vr*(one+gama+gama*gama*(Vt*Vt+Vp*Vp))/(one+gama)**two
                c1=(gama*gama*Vt*Vr/(one+gama))**two-gama*gama*Vt*Vt
                betap=(-b1+dsqrt(b1**two-four*a1*c1))/two/a1
                betam=(-b1-dsqrt(b1**two-four*a1*c1))/two/a1                         
                If(betap*Vp.lt.zero)then
                    betac=betap
                Else
                    betac=betam                
                Endif
                alphac=Vp/Vt*betac
            Else
                alphac=zero
                a1=(one+gama*gama*Vt*Vt/(one+gama))**two-gama*gama*Vt*Vt
                b1=two*gama*gama*Vt*Vr*(one+gama+gama*gama*Vt*Vt)/(one+gama)**two
                c1=(gama*gama*Vt*Vr/(one+gama))**two-gama*gama*Vt*Vt    
                betap=(-b1+dsqrt(b1**two-four*a1*c1))/two/a1
                betam=(-b1-dsqrt(b1**two-four*a1*c1))/two/a1        
                If(betap*Vt.lt.zero)then
                    betac=betap
                Else
                    betac=betam                
                Endif
            Endif  
        else
            betac=zero        
            If(Vp.ne.zero)then
                a1=(one+gama*gama*Vp*Vp/(one+gama))**two-gama*gama*Vp*Vp
                b1=two*gama*gama*Vp*Vr*(one+gama+gama*gama*Vp*Vp)/(one+gama)**two
                c1=(gama*gama*Vp*Vr/(one+gama))**two-gama*gama*Vp*Vp   
                alphap=(-b1+dsqrt(b1**two-four*a1*c1))/two/a1        
                alpham=(-b1-dsqrt(b1**two-four*a1*c1))/two/a1        
                If(alphap*Vp.lt.zero)then
                    alphac=alphap        
                Else
                    alphac=alpham        
                Endif
            Else
                alphac=zero
            Endif
        endif
        alphac=alphac*robs*scal
        betac=betac*robs*scal
        return
      End Subroutine center_of_image

!********************************************************************************************        
      Subroutine center_of_image_old(robs,sinobs,muobs,a_spin,scal,velocity,alphac,betac)
!********************************************************************************************
!*     PURPOSE:  Solves equations f_3(alphac,betac)=0, f_2(alphac,betac)=0, of (119)  
!*               and (120) in Yang & Wang (2012). alphac, betac are the coordinates of 
!*               center point of images on the screen of observer.    
!*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.  
!*               velocity(1:3)--Array of physical velocity of observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  alphac,betac---coordinates of center point of images on the screen of observer.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  9 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision robs,sinobs,muobs,a_spin,scal,velocity(3),alphac,betac,zero,one,two,four,&
                        Vr,Vt,Vp,gama,gama_tp,gama_p,a1,b1,c1,alphap,alpham,betap,betam
        parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0)

        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
! equation (92) in Yang & Wang (2012).
        gama=one/dsqrt(one-(Vr**two+Vt**two+Vp**two))
        gama_tp=one/dsqrt(one-(Vt**two+Vp**two))
        gama_p=one/dsqrt(one-Vp**two)        

        If(Vt.ne.zero)then
          If(Vp.ne.zero)then         
            a1=(gama_tp/gama_p)**two-((Vp/gama_tp)**two+Vt**two)*gama**two
            b1=two*gama*gama_tp*Vr*Vp/gama_p
            c1=-Vp**two
            alphap=(-b1+dsqrt(b1**two-four*a1*c1))/two/a1
            alpham=(-b1-dsqrt(b1**two-four*a1*c1))/two/a1                         
            If(alphap*Vp.lt.zero)then
                alphac=alphap
            Else
                alphac=alpham                
            Endif
            betac=alphac*Vt*gama_tp/Vp
          Else
            alphac=zero
            a1=one-(gama*Vt/gama_tp)**two
            b1=two*gama*Vr*Vt
            c1=(gama*Vt)**two*(Vr**two-one)        
            betap=(-b1+dsqrt(b1**two-four*a1*c1))/two/a1
            betam=(-b1-dsqrt(b1**two-four*a1*c1))/two/a1        
            If(betap*Vt.lt.zero)then
                betac=betap
            Else
                betac=betam                
            Endif
          Endif  
        else
            betac=zero        
            If(Vp.ne.zero)then
                a1=gama_p**two-(gama*Vp)**two
                b1=two*gama*gama_tp*gama_p*Vr*Vp
                c1=-(gama_tp*Vp)**two
                alphap=(-b1+dsqrt(b1**two-four*a1*c1))/two/a1        
                alpham=(-b1-dsqrt(b1**two-four*a1*c1))/two/a1        
                If(alphap*Vp.lt.zero)then
                    alphac=alphap        
                Else
                    alphac=alpham        
                Endif
            Else
                alphac=zero
            Endif
        endif
        alphac=alphac*robs*scal
        betac=betac*robs*scal
        return
      End Subroutine center_of_image_old

!********************************************************************************************
      FUNCTION p_total(f1234r,lambda,q,sinobs,muobs,a_spin,robs,scal)
!******************************************************************************************** 
!*     PURPOSE:  Computes the integral value of \int^r dr (R)^{-1/2}, from the starting position to
!*               the termination----either the infinity or the event horizon.   
!*     INPUTS:   f1234r---------p_r, the r component of four momentum of the photon measured 
!*                              under the LNRF, see equation (83) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  p_total--------which is the value of integrals \int^r dr (R)^{-1/2}, along a 
!*                              whole geodesic, that is from the starting position to either go to
!*                              infinity or fall in to black hole.          
!*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ******************************************   
        IMPLICIT NONE
        DOUBLE PRECISION phyr,re,a,B,sinobs,muobs,a_spin,rhorizon,q,lambda,integ,integ4(4),&
                          bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tobs,tp,pp,p1,p2,PI0,p1I0,p1J1,p1J2,E_add,E_m,&
                          u,v,w,L1,L2,thorizon,m2,pinf,r_add,r_m,B_add,B_m,D_add,D_m,&
                         y,x,f1,g1,h1,f2,h2,a5,b5,a4,b4,integ0,integ1,integ2,robs,ttp,&
                         PI1,PI2,scal,tinf,integ04(4),integ14(4),integ5(5),integ15(5),pp2,&
                         r_tp1,r_tp2,t_inf,tp2,f1234r,p_temp,PI0_obs_inf,PI0_total,PI0_obs_hori,&
                         PI0_obs_tp2,PI01,rff_p,p_t1_t2,p_total,p_tp1_tp2,PI2_p,PI1_p,sqt3         
        !PARAMETER(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0)
        COMPLEX*16 bb(1:4),dd(3)
        INTEGER ::  reals,i,j,t1,t2,p5,p4,index_p4(4),index_p5(5),del,cases_int,cases,count_num=1
        LOGICAL :: robs_eq_rtp,indrhorizon                

                rhorizon=one+dsqrt(one-a_spin**two)
! equation (64) in Yang & Wang (2012).
                r_add=rhorizon
                r_m=one-dsqrt(one-a_spin**two)      
! equation (64) in Yang & Wang (2012).   
                b4=one
                a4=zero
                cc=a_spin**2-lambda**2-q
                robs_eq_rtp=.false.
                indrhorizon=.false. 
                call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
                                  reals,robs_eq_rtp,indrhorizon,cases,bb) 
!** R(r)=0 has real roots and turning points exists in radial r.
                If(reals.ne.0)then  
! equations (35)-(38) in Yang & Wang (2012).
                    b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*&
                                    r_tp1+two*(q+(lambda-a_spin)**2)
                    b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
                    b2=four/three*r_tp1
                    b3=one
                    g2=three/four*(b1**2-b0*b2)
                    g3=one/16.D0*(three*b0*b1*b2-two*b1**three-b0**two*b3)
! equation (39) in Yang & Wang (2012).
                    If(robs-r_tp1.ne.zero)then        
                        tobs=b0/four/(robs-r_tp1)+b1/four
                    else
                        tobs=infinity
                    endif 
                    If(rhorizon-r_tp1.ne.zero)then
                        thorizon=b1/four+b0/four/(rhorizon-r_tp1)
                    else
                        thorizon=infinity         
                    endif
                    tp2=b0/four/(r_tp2-r_tp1)+b1/four   
                    tinf=b1/four       
! equation (64), (66) and (70) in Yang & Wang (2012).
                    call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)       

                    index_p4(1)=0
                    cases_int=1
                    call weierstrass_int_J3(tobs,infinity,dd,del,a4,&
                                 b4,index_p4,rff_p,integ04,cases_int) 
! equation (42) in Yang & Wang (2012).
                    PI0=integ04(1)   
                    select case(cases)
                    CASE(1)
                        If(f1234r .ge. zero)then !**photon will goto infinity.
                            index_p4(1)=0
                            cases_int=1
                            call weierstrass_int_J3(tinf,tobs,dd,del,a4,b4,&
                                           index_p4,rff_p,integ04,cases_int)
                            p_total = integ04(1)      
                        ELSE 
                            If(.not.indrhorizon)then
                                index_p4(1)=0
                                cases_int=1 
                                call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,&
                                                   index_p4,rff_p,integ14,cases_int)
                                p_total = PI0+integ14(1)      
                            ELSE     !f1234r<0, photon will fall into 
                                     !black hole unless something encountered.        
                                index_p4(1)=0                
                                cases_int=1
                                call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,&
                                                   index_p4,rff_p,integ04,cases_int)
                                p_total = integ04(1)  
                            ENDIF
                        ENDIF         
                    CASE(2)
                        If(.not.indrhorizon)then
                            write(*,*)'we come here!'  
                            If(f1234r.lt.zero)then
                                PI01=-PI0
                            else
                                PI01=PI0        
                            endif
! equation (41) in Yang & Wang (2012).    
                            index_p4(1)=0
                            cases_int=1         
                            call weierstrass_int_J3(tp2,infinity,dd,del,a4,b4,&
                                              index_p4,rff_p,integ14,cases_int)
                            pp=zero
! equation (57) in Yang & Wang (2012).
                            p_tp1_tp2=integ14(1) 
                            PI2_p=p_tp1_tp2-PI0 
                            PI1_p=PI0 
                            p1=PI0-pp
                            p2=p_tp1_tp2-p1   
! equation (58) in Yang & Wang (2012). 
                            t1 = 2
                            t2 = 2
                            If(robs_eq_rtp)then         
                                p_total=dabs(pp)+two*(t1*p1+t2*p2)
                            else
                                If(f1234r.gt.zero)then 
                                    p_total=-pp+two*(t1*p1+t2*p2)
                                endif
                                If(f1234r.lt.zero)then 
                                    p_total=pp+two*(t1*p1+t2*p2)
                                endif
                            endif   
                        !*************************************************************************************
                        200     continue                    
                        else  !photon has probability to fall into black hole.
                            If(f1234r.le.zero)then
                                index_p4(1)=0
                                cases_int=1
                                call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,&
                                                   index_p4,rff_p,integ04,cases_int)
                                p_total = integ04(1) 
                            ELSE  !p_r>0, photon will meet the r_tp2 turning &
                                  !point and turn around then goto vevnt horizon.     
                                index_p4(1)=0
                                cases_int=1        
                                call weierstrass_int_J3(tp2,tobs,dd,del,a4,b4,&
                                              index_p4,rff_p,integ04,cases_int)  
                                call weierstrass_int_j3(tp2,thorizon,dd,del,a4,b4,&
                                              index_p4,rff_p,integ14,cases_int) 
                                p_total = integ14(1)+integ04(1)        
                            ENDIF
                        ENDIF                             
                    END SELECT       
!************************************************************************************************ 
                    If(a_spin.eq.zero)then
                        If(cc.eq.zero)then
                            If(f1234r.lt.zero)then
                                p_total = one/rhorizon-one/robs
                            else
                                p_total = one/robs                       
                            endif        
                        endif
                        If(cc.eq.-27.D0)then
                            sqt3=dsqrt(three)        
                            If(f1234r.lt.zero)then
                                p_total=dLog(dabs((dsqrt(robs*(robs+6.D0))+&
                                (three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                                -dLog(one+two/sqt3)/three/sqt3 
                            else        
                                p_total=dLog(dabs((dsqrt(robs*(robs+6.D0))+&
                                (three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                                -dLog(one+two/sqt3)/three/sqt3                     
                            endif                
                        endif
                    endif                                   
                ELSE            
! equation (44) in Yang & Wang (2012).   !equation R(r)=0 has no real roots. we use the Legendre elliptic 
                    u=real(bb(4))        !integrations and functions to compute the calculations.
                    w=dabs(aimag(bb(4)))
                    v=dabs(aimag(bb(2)))
                    If(u.ne.zero)then
! equation (45) in Yang & Wang (2012).       
                        L1=(four*u**2+w**2+v**2+dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                        L2=(four*u**2+w**2+v**2-dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012).       
                        thorizon=dsqrt((L1-one)/(L1-L2))*(rhorizon-&
                                 u*(L1+one)/(L1-one))/dsqrt((rhorizon-u)**2+w**2)
! equation (48) in Yang & Wang (2012).   
                        m2=(L1-L2)/L1
                        tinf=dsqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/&
                             (L1-one))/dsqrt((robs-u)**two+w**two)
                        t_inf=dsqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012).       
                        pinf=EllipticF(tinf,m2)/(w*dsqrt(L1))
                        IF(f1234r.lt.zero)THEN
                            p_total = pinf-EllipticF(thorizon,m2)/(w*dsqrt(L1))             
                        ELSE
                            p_total = EllipticF(t_inf,m2)/(w*dsqrt(L1))-pinf        
                        ENDIF  
                    ELSE
                        If(f1234r.lt.zero)then
                            p_total = (datan(robs/w)-datan(rhorizon/w))/w         
                        ELSE
                            p_total = (PI/two-datan(robs/w))/w         
                        ENDIF   
                    ENDIF                        
                ENDIF             
        RETURN
     END FUNCTION p_total

!********************************************************************************************
      SUBROUTINE r2p_new( f1234r, r_outer, r_inner, lambda, q, a_spin, robs, &
                           p_out1, p_out2, p_in1, p_in2, InterSection_Cases, r_tp1, r_tp2 )
!============================================================================================
!*     PURPOSE:  Computes the value of parameter p from radial coordinate. In other words, to compute 
!*               the r part of integral of equation (24), using formula (58) in Yang & Wang (2012).
!*               (58) is: p=-sign(p_r)*p_0+2*t1*p_2+2*t2*p_2. where p_r is initial radial
!*               component of 4 momentum of photon. 
!*     INPUTS:   f1234r---------f_1, which was defined by equation (106) in Yang & Wang (2012). 
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               robs-----------radial coordinate of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               t1,t2----------Number of photon meets the turning points r_tp1 and r_tp2
!*                              respectively in radial motion.
!*     OUTPUTS:  r2p------------value of r part of integral (24) in Yang & Wang (2012).
!*     ROUTINES CALLED: radiustp, root3, weierstrass_int_j3, EllipticF.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision, intent(in) :: f1234r, lambda, q, a_spin, robs, r_outer, r_inner
      Double precision, intent(inout) :: p_out1, p_out2, p_in1, p_in2
      integer, intent(out) :: InterSection_Cases
      Double precision r2p_new1,p,rhorizon,zero,integ4(4),&
                    bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tinf,tinf1,PI0,PI1,cr,dr,integ04(4),&
                    u,v,w,L1,L2,thorizon,m2,pinf,sn,cn,dn,a4,b4,one,two,four,PI2,ttp,sqrt3,&
                    integ14(4),three,six,nine,r_tp1,r_tp2,tp2,tp,t_inf,PI0_total,&
                    PI0_inf_obs,PI0_obs_hori,PI01,PI0_total_2,pp,p1,p2,rff_p, &
                    t_in, t_out
      parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0,six=6.D0,nine=9.D0)
      complex*16 bb(1:4),dd(3)
      integer  reals,i,p4,cases_int,del,index_p4(4),cases,t1,t2
      logical :: robs_eq_rtp,indrhorizon
      
      rhorizon=one+dsqrt(one-a_spin**2)
      a4=zero
      b4=one
      cc=a_spin**2-lambda**2-q
      robs_eq_rtp=.false.
      indrhorizon=.false.
      call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,&
                    r_tp2,reals,robs_eq_rtp,indrhorizon,cases,bb)
      !write(*,*)'sss=',reals, cases, f1234r, indrhorizon, bb, robs
      If(reals.ne.0)then
          !If(r_outer.lt.r_tp1.or.r_outer.gt.r_tp2)then
              !r2p_new1 = -one
              !return
          !endif
! equations (35)-(38) in Yang & Wang (2012).
          b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
          b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
          b2=four/three*r_tp1
          b3=one
          g2=three/four*(b1**2-b0*b2)
          g3=one/16.D0*(3*b0*b1*b2-2*b1**3-b0**2*b3)
! equation (39) in Yang & Wang (2012).
          If(robs-r_tp1.ne.zero)then      
              tinf=b0/four/(robs-r_tp1)+b1/four
          else
              tinf=infinity
          endif 
          If(rhorizon-r_tp1.ne.zero)then
              thorizon=b1/four+b0/four/(rhorizon-r_tp1)
          else
              thorizon=infinity       
          endif
          If(r_outer-r_tp1.ne.zero)then
              t_out = b1/four+b0/four/(r_outer-r_tp1)
          else
              t_out = infinity       
          endif
          If(r_inner-r_tp1.ne.zero)then
              t_in = b1/four+b0/four/(r_inner-r_tp1)
          else
              t_in = infinity       
          endif
          tp2=b0/four/(r_tp2-r_tp1)+b1/four   
          tinf1=b1/four

          call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)      
          index_p4(1)=0      
          cases_int=1
! equation (42) in Yang & Wang (2012).
          call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,&
                             index_p4,rff_p,integ04,cases_int)      
          PI0=integ04(1) 
          select case(cases)
          case(1)
              If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then  
                        if( r_tp1 >= r_inner )then
                            call weierstrass_int_j3(tinf,t_out,dd,del,a4,b4,&
                                     index_p4,rff_p,integ14,cases_int)      
                            pp=integ14(1)   
                            call weierstrass_int_j3(t_out,infinity,dd,del,&
                                 a4,b4,index_p4,rff_p,integ14,cases_int)
                            p1=integ14(1)
                            p_out1 = pp + two * p1
                            InterSection_Cases = 1
                        else if( r_tp1 < r_inner )then    
                            call weierstrass_int_j3(tinf,t_in,dd,del,&
                                 a4,b4,index_p4,rff_p,integ14,cases_int)
                            !pp=integ14(1)
                            p_in1 = integ14(1)
                            call weierstrass_int_j3(t_in,infinity,dd,del,&
                                 a4,b4,index_p4,rff_p,integ14,cases_int)
                            p1=integ14(1)
                            p_in2 = p_in1 + two * p1
                            call weierstrass_int_j3(t_out,t_in,dd,del,&
                                 a4,b4,index_p4,rff_p,integ14,cases_int)
                            p_out1 = p_in2 + integ14(1)
                            InterSection_Cases = 3
                        endif
                    else
                        call weierstrass_int_J3(tinf,t_out,dd,del,a4,b4,&
                                     index_p4,rff_p,integ04,cases_int)
                        pp=integ04(1)
                        p_out1 = - pp
                        InterSection_Cases = 1
                    endif
                else
                    If(f1234r.lt.zero)then
                        If(r_inner .le. rhorizon)then 
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,&
                                          index_p4,rff_p,integ04,cases_int)        
                            p_in1 = integ04(1)
                        else
                            call weierstrass_int_J3(tinf,t_in,dd,del,a4,b4,&
                                          index_p4,rff_p,integ04,cases_int)        
                            p_in1 = integ04(1)
                        endif
                        !write(*,*)'r2p_new=', p_in1, tinf, t_in, thorizon
                        InterSection_Cases = - 1
                    else
                        If(r_outer .lt. infinity)then
                            call weierstrass_int_J3(tinf,t_out,dd,del,a4,b4,&
                                         index_p4,rff_p,integ04,cases_int) 
                            p_out1 = - integ04(1)
                        else
                            call weierstrass_int_J3(tinf,tinf1,dd,del,a4,b4,&
                                            index_p4,rff_p,integ04,cases_int) 
                            p_out1 = - integ04(1)
                        endif
                        !write(*,*)'r2p_new =', tinf, tp, tinf1, integ04(1)
                        InterSection_Cases = 1
                    endif                 
                endif
            case(2)
                If(.not.indrhorizon)then    
                        write(*,*)'r2p_new stops!! case 2 occured, it is wrong!'
                        stop         
                        call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,&
                                      index_p4,rff_p,integ4,cases_int)
                        pp=integ4(1)
                        If(t1.eq.zero)then
                            p1=zero
                        else
                            call weierstrass_int_J3(tp,infinity,dd,del,a4,b4,&
                                              index_p4,rff_p,integ4,cases_int)
                            p1=integ4(1)
                        endif        
                        If(t2.eq.zero)then
                            p2=zero
                        else
                            call weierstrass_int_J3(tp2,tp,dd,del,a4,b4,&
                                         index_p4,rff_p,integ4,cases_int)
                            p2=integ4(1)        
                        endif
                        If(f1234r.ne.zero)then
                            r2p_new1=dsign(one,-f1234r)*pp+two*(t1*p1+t2*p2)
                        else
! equation (58) in Yang & Wang (2012).
                            If(robs.eq.r_tp1)then
                                r2p_new1=-pp+two*(t1*p1+t2*p2)
                            else
                                r2p_new1=pp+two*(t1*p1+t2*p2)
                            endif        
                        endif                    
                else        
                    If(f1234r.le.zero)then
                        If(r_inner.le.rhorizon)then
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,&
                                          b4,index_p4,rff_p,integ4,cases_int)
                            p_in1 = integ4(1)                            
                        else
                            call weierstrass_int_J3(tinf,t_in,dd,del,a4,&
                                    b4,index_p4,rff_p,integ4,cases_int)
                            p_in1 = integ4(1)
                        endif
                        InterSection_Cases = - 1
                    else
                        call weierstrass_int_J3(tinf,tp2,dd,del,a4,b4,&
                                      index_p4,rff_p,integ4,cases_int)
                        pp = integ4(1) 
                        call weierstrass_int_J3(tp2,t_in,dd,del,a4,b4,&
                                         index_p4,rff_p,integ4,cases_int)
                        p2 = integ4(1)  
                        p_in1 = -pp + p2
                        InterSection_Cases = - 1
                        !write(*, *)'heress ==', pp, p2, p_in1
                        !stop
                    endif                              
                endif
            end select                    
            If(a_spin.eq.zero)then
                If(cc.eq.zero)then
                    write(*,*)'r2p_new stops!! cc==0 occured, it is wrong!'
                    stop 
                    If(f1234r.lt.zero)then
                        If(r_inner.le.rhorizon)then
                            r2p_new1=one/rhorizon-one/robs
                        else
                            r2p_new1=one/r_inner-one/robs         
                        endif
                    else
                        If(r_outer.lt.infinity)then
                            r2p_new1=one/robs-one/r_outer
                        else
                            r2p_new1=one/robs
                        endif 
                    endif 
                endif
                If(cc.eq.-27.D0)then  
                    write(*,*)'r2p_new stops!! cc==-27.D0 occured, it is wrong!'
                    stop       
                    sqrt3=dsqrt(three)        
                    If(f1234r.lt.zero)then
                        If(r_inner.gt.rhorizon)then
                            r2p_new1=-dlog(dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/(sqrt3))/&
                                                         (robs-three)))/(three*sqrt3)+&
                                 dlog(dabs((dsqrt(r_inner*(r_inner+6.D0))+(three+two*r_inner)/&
                                                         (sqrt3))/(r_inner-three)))/(three*sqrt3)
                        else
                            r2p_new1=-dlog(dabs((dsqrt(robs*(robs+6.D0))+(three+two*robs)/&
                                            (sqrt3))/(robs-three)))/(three*sqrt3)+&
                                 dlog(dabs((dsqrt(rhorizon*(rhorizon+6.D0))+(three+two*rhorizon)&
                                            /(sqrt3))/(rhorizon-three)))/(three*sqrt3)
                        endif
                    else        
                        If(r_outer.lt.infinity)then
                            r2p_new1=-dlog(dabs((dsqrt(r_outer*(r_outer+6.D0))+(three+two*r_outer)/(sqrt3))/&
                                (r_outer-three)))/(three*sqrt3)+dlog(dabs((dsqrt(robs*(robs+6.D0))+&
                                (three+two*robs)/(sqrt3))/(robs-three)))/(three*sqrt3)
                        else
                            r2p_new1=-dlog(one+two/sqrt3)/three/sqrt3+&
                                 dlog(dabs((dsqrt(r_outer*(r_outer+6.D0))+(three+two*r_outer)/&
                                    (sqrt3))/(r_outer-three)))/(three*sqrt3)
                        endif                                                
                    endif                
                endif
            endif        
        else
! equation (44) in Yang & Wang (2012).
            !write(*,*)'r2p_new stops!! R(r)==0 has no real roots occured, it is wrong!'
            !stop 
            u=real(bb(4))
            w=dabs(aimag(bb(4)))
            v=dabs(aimag(bb(2)))
            !write(*,*)'r2p_new stops!! R(r)==0 has no real roots occured, it is wrong!', u, w, v
            If(u.ne.zero)then
! equation (45) in Yang & Wang (2012).
                L1=(four*u**2+w**2+v**2+dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                L2=(four*u**2+w**2+v**2-dsqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012).
                thorizon=dsqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/&
                               (L1-one))/dsqrt((rhorizon-u)**2+w**2)
                t_in=dsqrt((L1-one)/(L1-L2))*(r_inner-u*(L1+one)/&
                         (L1-one))/dsqrt((r_inner-u)**2+w**2)
                t_out=dsqrt((L1-one)/(L1-L2))*(r_outer-u*(L1+one)/&
                         (L1-one))/dsqrt((r_outer-u)**2+w**2)
! equation (48) in Yang & Wang (2012).
                m2=(L1-L2)/L1
                tinf=dsqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/&
                           (L1-one))/dsqrt((robs-u)**two+w**two)
                t_inf=dsqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012).
                pinf=EllipticF(tinf,m2)/w/dsqrt(L1) 
                If(f1234r.lt.zero)then
                    If(r_inner.le.rhorizon)then
                        p_in1=pinf-EllipticF(thorizon,m2)/(w*dsqrt(L1))
                    else
                        p_in1=pinf-EllipticF(t_in,m2)/(w*dsqrt(L1))
                    endif   
                    InterSection_Cases = - 1                         
                else
                    If(r_outer.lt.infinity)then
                        p_out1=EllipticF(t_out,m2)/(w*dsqrt(L1))-pinf
                    else
                        p_out1=EllipticF(t_inf,m2)/(w*dsqrt(L1))-pinf
                    endif
                    InterSection_Cases = 1
                endif
            else
                If(f1234r.lt.zero)then
                    If(r_inner.le.rhorizon)then
                        p_in1=(datan(robs/w)-datan(rhorizon/w))/w
                    else
                        p_in1=(datan(robs/w)-datan(r_inner/w))/w
                    endif 
                    InterSection_Cases = - 1       
                else
                    if(r_outer.lt.infinity)then
                        p_out1=(datan(r_outer/w)-datan(robs/w))/w
                    else
                        p_out1=(PI/two-datan(robs/w))/w
                    endif  
                    InterSection_Cases = 1              
                endif
            endif                        
      endif        
      return                 
      End SUBROUTINE r2p_new

!********************************************************************************************        
      end module BLcoordinate 




