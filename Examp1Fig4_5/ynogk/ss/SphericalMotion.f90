!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      module sphericalmotion
!********************************************************************
!* This module aim on the calculation of the geodesic trajectory of a photon
!* doing spherical motion.
!********************************************************************
      USE blcoordinate
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      contains 
!********************************************************************************************
      SUBROUTINE SPHERICALMOTION_BL(p,kp,kt,lambda,q,sinobs,muobs,a_spin,&
                               robs,thetamax,phyt,timet,sigmat,mucos,t1,t2)    
!******************************************************************************************** 
!*     PURPOSE:   This routine computs the four B_L coordiants r,\mu,\phi,t and affine 
!*                parameter \sigma of the spherical motion. 
!*     INPUTS:   p--------------the independent variable.
!*               kt-------------p_\theta, the \theta component of four momentum of the photon.
!*               kp-------------p_\phi,  the \phi component of four momentum of the photon.
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{ini}), muobs=cos(\theta_{ini}), where 
!*                              \theta_{ini} is the initial theta angle of the photon.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               thetamax-------the maximum or minimum value of the \theta coordinate of the geodesic.
!*                              which also the \theta turning point of the spherical motion.       
!*     OUTPUTS:  phyt-----------the \phi coordinat of the photon.
!*               timet----------the t coordinats of the photon.
!*               sigmat---------the affine parameter \sigma.
!*               mucos----------the \theta coordinates of the photon, and mucos=cos(\theta).  
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.      
!*     ROUTINES CALLED: metricg, circ_mb, weierstrass_int_J3, weierstrassP
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
      IMPLICIT NONE
      Double precision phyt,timet,kp,kt,p,sinobs,muobs,a_spin,lambda,q,mu_tp1,tposition,tp2,mu,tmu,p1J2,&
             bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tobs,p1,p2,pp,c_add,c_m,a_add,a_m,p1J1,come,&
             p1I0,a4,b4,delta,mu_tp2,robs,mutemp ,integ5(5),integ(5),rff_p,tp1,&
             integ15(5),pp2,f1234(4),PI0,integ05(5),fzero,p_mu,PI01,h,p1_t,p2_t,pp_t,p1_phi,&
             p2_phi,pp_phi,radius,mtp1,mtp2,mucos,sqt3,difference,p_mt1_mt2,PI1_sig,PI2_sig,&
             PI1_phi,PI2_phi,PI1_time,PI2_time,PI2_p,bigT,p1_sig,p2_sig,pp_sig,sigmat 
      Double precision kp_1,kt_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,&
             c_phi,c_time,thetamax,thetamax_1,sinmax,mumax,Omega,somiga,&
             expnu,exppsi,expmu1,expmu2,c_tau 
      integer ::  t1,t2,i,j,reals,cases,p4,index_p5(5),del,cases_int,count_num=1
      complex*16 bb(1:4),dd(3)
      logical :: err,mobseqmtp
      save  kp_1,kt_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,a4,b4,mu_tp1,mu_tp2,reals,&
            mobseqmtp,b0,b1,b2,b3,g2,g3,dd,del,PI0,c_m,c_add,a_m,a_add,tp2,tobs,h,p_mt1_mt2,&
            PI1_phi,PI2_phi,PI1_time,PI2_time,PI2_p,come,tp1,bigT,Delta,c_phi,c_time,&
            PI01,sinmax,mumax,thetamax_1,Omega,c_tau

30    continue
      IF(count_num.eq.1)then	 
          kp_1=kp
          kt_1=kt
          lambda_1=lambda
          q_1=q
          muobs_1=muobs
          sinobs_1=sinobs
          a_spin_1=a_spin
          robs_1=robs
          thetamax_1 = thetamax
          t1=0
          t2=0  
          IF(thetamax.ne.90.D0 .and. thetamax.ne.180.D0)THEN
              sinmax = dsin(thetamax*dtor)
              mumax = dcos(thetamax*dtor)
          Else
              IF(thetamax.eq.90.D0)THEN
                  sinmax = one
                  mumax = zero
              ENDIF
              IF(thetamax.eq.180.D0)THEN
                  sinmax = zero
                  mumax = one
              ENDIF
          ENDIF
          mu_tp1 = abs(mumax)
          mu_tp2 = -mu_tp2
     !***********************************************************************    
          If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, so the integrations about \theta are valished.
              Omega = one/(robs**(three/two)+a_spin)
              mucos=zero
              phyt=p
              timet=phyt/Omega
              call metricg(robs,sinobs,muobs,a_spin,somiga,expnu,exppsi,expmu1,expmu2)
              c_tau = dsqrt(-expnu*expnu+somiga*somiga*exppsi*exppsi-two*exppsi*somiga*Omega+&
                            exppsi*exppsi*Omega*Omega)
              sigmat=c_tau*timet
              count_num=count_num+1 
              return
          endif
     !**************************************************************************
          If(a_spin.eq.zero)then
              timet=zero
              CALL SPINZERO(p,kp,kt,lambda,q,sinobs,muobs,a_spin,&
                           robs,thetamax,phyt,timet,sigmat,mucos,t1,t2) 
              count_num=count_num+1
              return
          endif
          a4=zero
          b4=one
!      equations (26)-(29) of Yang and Wang (2013).
          come = -one  
          b0=four*a_spin**2*mu_tp1**3*come-two*mu_tp1*(a_spin**2*come+lambda**2+q)
          b1=two*a_spin**2*mu_tp1**2*come-(a_spin**2*come+lambda**2+q)/three
          b2=four/three*a_spin**2*mu_tp1*come
          b3=a_spin**2*come
!      equations (31) of Yang and Wang (2013).
          g2=three/four*(b1**2-b0*b2)
          g3=(three*b0*b1*b2-two*b1**3-b0**2*b3)/16.D0  
          call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)

!      equations (30) of Yang and Wang (2013).
          If(muobs.ne.mu_tp1)then	
              tobs=b0/four/(muobs-mu_tp1)+b1/four
          else
              tobs=infinity
          endif
          tp1=infinity
          tp2=b0/four/(mu_tp2-mu_tp1)+b1/four
!      equations (72)-(73) of Yang and Wang (2013).
          If(mu_tp1-one.ne.zero)then
               c_m=b0/(four*(-one-mu_tp1)**2)
             c_add=b0/(four*(one-mu_tp1)**2) 
               a_m=b0/four/(-one-mu_tp1)+b1/four
             a_add=b0/four/(one-mu_tp1)+b1/four
          endif
          index_p5(1)=0
          cases_int=1
!      equations (53) of Yang and Wang (2013).
          call weierstrass_int_J3(tobs,tp1,dd,del,a4,b4,index_p5,rff_p,integ05,cases_int)
          PI0=integ05(1)	
          If(kt.lt.zero)then
              PI01=-PI0	
          else
              PI01=PI0
          endif
          tmu=weierstrassP(p+PI01,g2,g3,dd,del)
!      equations (32) of Yang and Wang (2013).
          mucos = mu_tp1+b0/(four*tmu-b1)
          h=-b1/four	
          !to get number of turn points of t1 and t2.
          !111111111*****************************************************************************************
          !mu=mu_tp+b0/(four*tmu-b1)
          call weierstrass_int_J3(tp2,tp1,dd,del,a4,b4,index_p5,rff_p,integ15,cases_int)
          call weierstrass_int_J3(tobs,tmu,dd,del,a4,b4,index_p5,rff_p,integ5,cases_int)
!      equations (53) of Yang and Wang (2013).	
          p_mt1_mt2=integ15(1)
          PI2_p=p_mt1_mt2-PI0
          pp=integ5(1)
          p1=PI0-pp
!      equations (53) of Yang and Wang (2013).
          p2=p_mt1_mt2-p1
          PI1_phi=zero
          PI2_phi=zero
          PI1_sig=zero
          PI2_sig=zero
          PI1_time=zero
          PI2_time=zero 
!      equations (54) of Yang and Wang (2013).
          Do j=0,10**6
              Do i=j,j+1
                  If(mobseqmtp)then
                      If(muobs.eq.mu_tp1)then
                          t1=j
                          t2=i
                          p_mu=-pp+two*(t1*p1+t2*p2)
                      else
                          t1=i
                          t2=j
                          p_mu=pp+two*(t1*p1+t2*p2)
                      endif
                  else	
                      If(kt.lt.zero)then	
                          t1=i
                          t2=j
                          p_mu=pp+two*(t1*p1+t2*p2)
                      endif
                      If(kt.gt.zero)then	
                          t1=j
                          t2=i	
                          p_mu=-pp+two*(t1*p1+t2*p2)			
                      endif
                  endif    
                  If(dabs(p-p_mu).lt.1.D-3)goto 400
              enddo
          enddo
          !11111111*************************************************************** 
400       continue
!      equations (71)-(73) of Yang and Wang (2013).
          Delta=robs*robs+a_spin*a_spin-two*robs
          c_phi = a_spin*(robs*(two)-(lambda*a_spin))/Delta
          c_time = ( (two)*robs**three+(two*a_spin*(a_spin-lambda))*robs )/Delta
          index_p5(1)=-1
          index_p5(2)=-2
          index_p5(3)=0
          index_p5(4)=-4
          index_p5(5)=0
          !*****pp part***************************************
          If(lambda.ne.zero)then	
              cases_int=2
              call weierstrass_int_J3(tobs,tmu,dd,del,-a_add,b4,index_p5,abs(pp),integ5,cases_int)
              call weierstrass_int_J3(tobs,tmu,dd,del,-a_m,b4,index_p5,abs(pp),integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
              pp_phi=(pp/(one-mu_tp1**two)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*pp 
          else 
              pp_phi=c_phi*pp 		
          endif
          cases_int=4
          call weierstrass_int_J3(tobs,tmu,dd,del,h,b4,index_p5,abs(pp),integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
          pp_sig=(pp*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                  sixteen)*a_spin*a_spin+robs*robs*pp
!      equations (22) of Yang and Wang (2013).
          pp_t=pp_sig+c_time*pp  
          !*****p1 part***************************************
          If(t1.eq.0)then	
              p1_phi=zero
              p1_sig=zero 
              p1_t=zero
          else  
              If(lambda.ne.zero)then  
                  IF(PI1_phi .EQ. zero)THEN
                      cases_int=2	
                      call weierstrass_int_J3(tobs,infinity,dd,del,-a_add,b4,index_p5,PI0,integ5,cases_int)
                      call weierstrass_int_J3(tobs,infinity,dd,del,-a_m,b4,index_p5,PI0,integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
                      PI1_phi=(PI0/(one-mu_tp1**two)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*PI0
                  ENDIF 
                  p1_phi=PI1_phi-pp_phi 
              else 
                  IF(PI1_phi.eq.zero)PI1_phi=c_phi*PI0
                  p1_phi=PI1_phi-pp_phi	    
              endif 
              IF(PI1_time .EQ. zero .or. PI1_sig.eq.zero)THEN  
                  cases_int=4  
                  call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p5,PI0,integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
                  PI1_sig=(PI0*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                           sixteen)*a_spin*a_spin+robs*robs*PI0
!      equations (22) of Yang and Wang (2013).
                  PI1_time=PI1_sig+c_time*PI0  
              ENDIF
              p1_sig=PI1_sig-pp_sig
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
                      call weierstrass_int_J3(tp2,tobs,dd,del,-a_add,b4,index_p5,PI2_p,integ5,cases_int)
                      call weierstrass_int_J3(tp2,tobs,dd,del,-a_m,b4,index_p5,PI2_p,integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
                      PI2_phi=(PI2_p/(one-mu_tp1*mu_tp1)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*PI2_p
                  ENDIF
                  p2_phi=PI2_phi+pp_phi
              ELSE
                  IF(PI2_phi.eq.zero)PI2_phi=c_phi*PI2_p  
                  p2_phi=PI2_phi+pp_phi             
              ENDIF 
                
              IF(PI2_time .EQ. zero)THEN  
                  cases_int=4  
                  call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p5,PI2_p,integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
                  PI2_sig=(PI2_p*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                                  sixteen)*a_spin*a_spin+robs*robs*PI2_p
!      equations (22) of Yang and Wang (2013).
                  PI2_time=PI2_sig+c_time*PI2_p  
              ENDIF	
              p2_sig=PI2_sig+pp_sig
              p2_t=PI2_time+pp_t   
          endif   
          !write(*,*)pp_phi,p1_phi,p2_phi,t1,t2
	!**************************************************************
          If(mobseqmtp)then 
              If(muobs .eq. mu_tp1)then
!      equations (52) of Yang and Wang (2013).
                  phyt= -pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet= -pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat= -pp_sig+two*(t1*p1_sig+t2*p2_sig) 
              else
!      equations (52) of Yang and Wang (2013).
                  phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig)  
              endif
          else
              If(kt.lt.zero)then
!      equations (52) of Yang and Wang (2013).
                  phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)	
                  timet=pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig)
              endif
              If(kt.gt.zero)then
!      equations (52) of Yang and Wang (2013).		
                  phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=-pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=-pp_sig+two*(t1*p1_sig+t2*p2_sig)	
              endif
          endif
          IF(mu_tp1.eq.one)phyt = phyt+(t1+t2)*PI
          !phyt = mod(phyt,twopi)
          !If(phyt .lt. zero)phyt=phyt+twopi
          count_num=count_num+1
      ELSE  
          If(kp_1.eq.kp.and.kt_1.eq.kt.and.lambda_1.eq.lambda.and.q_1.eq.q.and.&
          sinobs_1.eq.sinobs.and.muobs_1.eq.muobs.and.a_spin_1.eq.a_spin.and.robs_1.eq.robs&
          .and.thetamax_1.eq.thetamax)then   
	!***************************************************************************
          t1=0
          t2=0  
        !***********************************************************************   	
          If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, so the integrations about \theta are valished. 
              mucos=zero
              phyt=p
              timet=phyt/Omega 
              sigmat=c_tau*timet 
              return
          endif
          !**************************************************************************
          If(a_spin.eq.zero)then 
              CALL SPINZERO(p,kp,kt,lambda,q,sinobs,muobs,a_spin,&
                           robs,thetamax,phyt,timet,sigmat,mucos,t1,t2) 
              return
          endif   
          tmu=weierstrassP(p+PI01,g2,g3,dd,del)
          mucos = mu_tp1+b0/(four*tmu-b1) 	
          !to get number of turn points of t1 and t2.
          !111111111*****************************************************************************************
          !mu=mu_tp+b0/(four*tmu-b1) 
          index_p5(1)=0
          cases_int=1
          call weierstrass_int_J3(tobs,tmu,dd,del,a4,b4,index_p5,rff_p,integ5,cases_int)  
          pp=integ5(1)
          p1=PI0-pp
          p2=p_mt1_mt2-p1 
          Do j=0,10**6
              Do i=j,j+1
                  If(mobseqmtp)then
                      If(muobs.eq.mu_tp1)then
                          t1=j
                          t2=i
                          p_mu=-pp+two*(t1*p1+t2*p2)
                      else
                          t1=i
                          t2=j
                          p_mu=pp+two*(t1*p1+t2*p2)
                      endif
                  else	
                      If(kt.lt.zero)then	
                          t1=i
                          t2=j
                          p_mu=pp+two*(t1*p1+t2*p2)
                      endif
                      If(kt.gt.zero)then	
                          t1=j
                          t2=i	
                          p_mu=-pp+two*(t1*p1+t2*p2)			
                      endif
                  endif  
                  !write(*,*)p,p_mu,abs(p-p_mu),t1,t2!(a,B,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,robs),t1,t2
                  If(abs(p-p_mu).lt.1.D-3)goto 410
              enddo
          enddo
          !11111111*************************************************************** 
410       continue
          index_p5(1)=-1
          index_p5(2)=-2
          index_p5(3)=0
          index_p5(4)=-4
          index_p5(5)=0
          !*****pp part***************************************
          If(lambda.ne.zero)then	
              cases_int=2
              call weierstrass_int_J3(tobs,tmu,dd,del,-a_add,b4,index_p5,abs(pp),integ5,cases_int)
              call weierstrass_int_J3(tobs,tmu,dd,del,-a_m,b4,index_p5,abs(pp),integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
              pp_phi=(pp/(one-mu_tp1**two)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*pp 
          else 
              pp_phi=c_phi*pp		
          endif
          cases_int=4
          call weierstrass_int_J3(tobs,tmu,dd,del,h,b4,index_p5,abs(pp),integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
          pp_sig=(pp*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                    sixteen)*a_spin*a_spin+robs*robs*pp
!      equations (22) of Yang and Wang (2013).
          pp_t=pp_sig+c_time*pp
          !*****p1 part***************************************
          If(t1.eq.0)then	
              p1_phi=zero
              p1_t=zero
          else  
              If(lambda.ne.zero)then  
                  IF(PI1_phi .EQ. zero)THEN
                      cases_int=2	
                      call weierstrass_int_J3(tobs,infinity,dd,del,-a_add,b4,index_p5,PI0,integ5,cases_int)
                      call weierstrass_int_J3(tobs,infinity,dd,del,-a_m,b4,index_p5,PI0,integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
                      PI1_phi=(PI0/(one-mu_tp1**two)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*PI0
                  ENDIF 
                  p1_phi=PI1_phi-pp_phi 
              else 
                  IF(PI1_phi.eq.zero)PI1_phi=c_phi*PI0	
                  p1_phi=PI1_phi-pp_phi     
              endif 
              IF(PI1_time .EQ. zero)THEN  
                  cases_int=4  
                  call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p5,PI0,integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
                  PI1_sig=(PI0*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                                  sixteen)*a_spin*a_spin+robs*robs*PI0
!      equations (22) of Yang and Wang (2013).
                  PI1_time=PI1_sig+c_time*PI0  
              ENDIF
              p1_sig=PI1_sig-pp_sig
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
                      call weierstrass_int_J3(tp2,tobs,dd,del,-a_add,b4,index_p5,PI2_p,integ5,cases_int)
                      call weierstrass_int_J3(tp2,tobs,dd,del,-a_m,b4,index_p5,PI2_p,integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
                      PI2_phi=(PI2_p/(one-mu_tp1*mu_tp1)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*PI2_p 
                  ENDIF
                  p2_phi=PI2_phi+pp_phi
              ELSE
                  IF(PI2_phi.eq.zero)PI2_phi=c_phi*PI2_p 
                  p2_phi=PI2_phi+pp_phi           
              ENDIF 
                
              IF(PI2_time .EQ. zero)THEN  
                  cases_int=4  
                  call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p5,PI2_p,integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
                  PI2_sig=(PI2_p*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                                  sixteen)*a_spin*a_spin+robs*robs*PI2_p
!      equations (22) of Yang and Wang (2013).
                  PI2_time=PI2_sig+c_time*PI2_p   
              ENDIF	
              p2_sig=PI2_sig+pp_sig 
              p2_t=PI2_time+pp_t   
          endif   
          !write(*,*)'kkk=',pp_phi,p1_phi,p2_phi,t1,t2
	!**************************************************************
          If(mobseqmtp)then 
              If(muobs .eq. mu_tp1)then
!      equations (52) of Yang and Wang (2013).
                  phyt= -pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet= -pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat= -pp_sig+two*(t1*p1_sig+t2*p2_sig) 
              else
!      equations (52) of Yang and Wang (2013).
                  phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig)  
              endif
          else
              If(kt.lt.zero)then
!      equations (52) of Yang and Wang (2013).
                  phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)	
                  timet=pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig)
              endif
              If(kt.gt.zero)then
!      equations (52) of Yang and Wang (2013).		
                  phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=-pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=-pp_sig+two*(t1*p1_sig+t2*p2_sig)	
              endif
          endif 
          IF(mu_tp1.eq.one)phyt = phyt+(t1+t2)*PI 
       !***************************************************** 
          else
              count_num=1
              goto 30  	
          endif				
      ENDIF 
      RETURN
      END SUBROUTINE SPHERICALMOTION_BL
!********************************************************************************************
      SUBROUTINE SPINZERO(p,kp,kt,lambda,q,sinobs,muobs,a_spin,robs,&
                         thetamax,phyt,timet,sigmat,mucos,t1,t2) 
!******************************************************************************************** 
!*     PURPOSE:   This routine computs the four B_L coordiants r,\mu,\phi,t and affine 
!*                parameter \sigma of the spherical motion when black hole spin a is zero. 
!*     INPUTS:   p--------------the independent variable.
!*               kt-------------p_\theta, the \theta component of four momentum of the photon.
!*               kp-------------p_\phi,  the \phi component of four momentum of the photon.
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               thetamax-------the maximum or minimum value of the \theta coordinate of the geodesic.
!*                              which also the \theta turning point of the spherical motion.       
!*     OUTPUTS:  phyt-----------the \phi coordinat of the photon.
!*               timet----------the t coordinats of the photon.
!*               sigmat---------the affine parameter \sigma.
!*               mucos----------the \theta coordinates of the photon, and mucos=cos(\theta).  
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.      
!*     ROUTINES CALLED: metricg, circ_mb, weierstrass_int_J3, weierstrassP
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
      use constants 
      IMPLICIT NONE
      DOUBLE PRECISION kp,kt,lambda,q,p,sinobs,muobs,a_spin,robs,phyt,timet,mucos,&
             kp_1,kt_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,sigmat,AA,BB,&
             mu_tp1,mu_tp2,PI1,Ptotal,pp,p1,p2,PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_sigma,&
             PI2_sigma,pp_sigma,p1_sigma,p2_sigma,pp_time,p1_time,p2_time,pp_phi,p1_phi,&
             p2_phi,p_mu,Delta,c_time,c_phi,PI2,mu,thetamax,thetamax_1
      integer :: t1,t2,count_num=1,i,j
      save :: kp_1,kt_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,PI1_phi,PI2_phi,&
              PI1_time,PI2_time,PI1_sigma,PI2_sigma,Delta,c_time,c_phi,PI1,PI2,mobseqmtp,&
              AA,BB,mu_tp1,mu_tp2,Ptotal,thetamax_1
      logical :: mobseqmtp 

30    continue
      IF(count_num .eq. 1)THEN
          kp_1 = kp
          kt_1 = kt
          lambda_1 = lambda
          q_1 =q
          sinobs_1 = sinobs
          muobs_1 = muobs
          a_spin_1 = a_spin
          robs_1 = robs
          thetamax_1 = thetamax
          t1=0
          t2=0
          mobseqmtp = .false.
 
          IF(q.gt.zero)THEN
              AA = dsqrt((q+lambda*lambda)/q)
              BB = dsqrt(q)
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              If(kt.gt.zero)then
                  mu=sin(asin(muobs*AA)-p*BB*AA)/AA	
              else
                  If(kt.eq.zero)then
                      mu=cos(p*AA*BB)*muobs
                  else			      
                      mu=sin(asin(muobs*AA)+p*AA*BB)/AA	
                  endif	
              endif
              mucos = mu  
          !****************************************************
              If(kt.ne.zero)then
                  mu_tp1=sqrt(q/(lambda**two+q))
                  mu_tp2=-mu_tp1	
              else
                  mu_tp1=abs(muobs)
                  mu_tp2=-mu_tp1
                  mobseqmtp=.true.
              endif
              If(abs(muobs).eq.one)mobseqmtp=.true.
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
              If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, 
              !so the integrations about !\theta are valished.
                  timet = zero
                  sigmat = zero
                  phyt = zero
                  return
              endif  
          !***************************************************
              PI1=(PI/two-asin(muobs/mu_tp1))*mu_tp1/BB	
              Ptotal=PI*mu_tp1/BB
              PI2=Ptotal-PI1	
              pp=(asin(mu/mu_tp1)-asin(muobs/mu_tp1))*mu_tp1/BB	
              p1=PI1-pp
              p2=Ptotal-p1 
              PI1_phi=zero
              PI2_phi=zero
              PI1_time=zero
              PI2_time=zero
              PI1_sigma=zero
              PI2_sigma=zero
              Do j=0,10**6
                  Do i=j,j+1
                      If(mobseqmtp)then
                          If(muobs.eq.mu_tp1)then
                              t1=j
                              t2=i
                              p_mu=-pp+two*(t1*p1+t2*p2)
                          else
                              t1=i
                              t2=j
                              p_mu=pp+two*(t1*p1+t2*p2)
                          endif 
                      else	
                          If(kt.lt.zero)then	
                              t1=i
                              t2=j
                              p_mu=pp+two*(t1*p1+t2*p2)
                          endif
                          If(kt.gt.zero)then	
                              t1=j
                              t2=i	
                              p_mu=-pp+two*(t1*p1+t2*p2)			
                          endif
                      endif   
                      If(abs(p-p_mu).lt.1.D-6)goto 300
                  enddo
              enddo
300           continue
              !p0 part!
              Delta = robs*robs+a_spin*a_spin-two*robs
              c_phi = a_spin*(robs*(two)-(lambda*a_spin))/Delta
              c_time = ( (two)*robs**three+(two*a_spin*(a_spin-lambda))*robs )/Delta

              pp_sigma = a_spin*a_spin*mveone_int(muobs,mu,one/AA)/BB+robs*robs*pp
              pp_time = pp_sigma+c_time*pp
              pp_phi = lambda*schwatz_int(muobs,mu,AA)/BB+c_phi*pp
              !******p1 part***********************************
              IF(t1 .eq. 0)THEN
                  p1_sigma=zero
                  p1_time=zero
                  p1_phi=zero
              ELSE
                  IF(PI1_time .eq. zero .or. PI1_sigma .eq. zero)THEN
                      PI1_sigma = a_spin*a_spin*mveone_int(muobs,mu_tp1,one/AA)/BB+robs*robs*PI1
                      PI1_time = PI1_sigma+c_time*PI1 
                  ENDIF
                  IF(PI1_phi .eq. zero)THEN
                      PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB+c_phi*PI1
                  ENDIF
                  p1_sigma = PI1_sigma-pp_sigma
                  p1_time = PI1_time-pp_time
                  p1_phi = PI1_phi-pp_phi
              ENDIF  
              !******p2 part***********************************
              IF(t2 .eq. 0)THEN
                  p2_sigma=zero
                  p2_time=zero
                  p2_phi=zero
              ELSE
                  IF(PI2_time .eq. zero .or. PI2_sigma .eq. zero)THEN
                      PI2_sigma = a_spin*a_spin*mveone_int(mu_tp2,muobs,one/AA)/BB+robs*robs*PI2
                      PI2_time = PI2_sigma+c_time*PI2 
                  ENDIF
                  IF(PI2_phi .eq. zero)THEN
                      PI2_phi = lambda*schwatz_int(mu_tp2,muobs,AA)/BB+c_phi*PI2
                  ENDIF
                  p2_sigma = PI2_sigma+pp_sigma
                  p2_time = PI2_time+pp_time
                  p2_phi = PI2_phi+pp_phi
              ENDIF
              !**********************************************  
              If(mobseqmtp)then
                  If(muobs.eq.mu_tp1)then  
!      equations (52) of Yang and Wang (2013).	
                      sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = -pp_time+two*(t1*p1_time+t2*p2_time)	
                      phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi)	
                  else
!      equations (52) of Yang and Wang (2013).	
                      sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = pp_time+two*(t1*p1_time+t2*p2_time)	
                      phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi)		
                  endif  
              else
                 If(kt.lt.zero)then
!      equations (52) of Yang and Wang (2013).	
                      sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = pp_time+two*(t1*p1_time+t2*p2_time)	
                      phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi)	
                 endif 	
                 If(kt.gt.zero)then
!      equations (52) of Yang and Wang (2013).	
                      sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = -pp_time+two*(t1*p1_time+t2*p2_time)	
                      phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi)	
                 endif	  
              endif
              If(thetamax.eq.zero.or.thetamax.eq.180.D0)phyt = phyt+(t1+t2)*PI
          ELSE
              !write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is a affending',&
              !		'value, the program should be',&  
              !		'stoped! and q = ',q
              !stop
              mucos=muobs 
              t1 = 0
              t2 = 0
              phyt = zero 
              timet = zero
              sigmat = zero
          ENDIF
      ELSE 
          IF(kp_1 .eq. kp.and.kt_1 .eq. kt.and.lambda_1 .eq. lambda.and.q_1 .eq.q.and.&
          sinobs_1 .eq. sinobs.and.muobs_1 .eq. muobs.and.a_spin_1 &
          .eq. a_spin.and.robs_1 .eq. robs .and.thetamax_1.eq.thetamax)THEN
      !*******************************************************************************
            IF(q.gt.zero)THEN 
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              If(kt.gt.zero)then
                  mu=sin(asin(muobs*AA)-p*BB*AA)/AA	
              else
                  If(kt.eq.zero)then
                      mu=cos(p*AA*BB)*muobs
                  else			      
                      mu=sin(asin(muobs*AA)+p*AA*BB)/AA	
                  endif	
              endif
              mucos = mu   
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
              If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, 
              !so the integrations about !\theta are valished.
                  timet = zero
                  sigmat = zero
                  phyt = zero
                  return
              endif  
          !*************************************************** 
              pp=(asin(mu/mu_tp1)-asin(muobs/mu_tp1))*mu_tp1/BB	
              p1=PI1-pp
              p2=Ptotal-p1  
              Do j=0,10**6
                  Do i=j,j+1
                      If(mobseqmtp)then
                          If(muobs.eq.mu_tp1)then
                              t1=j
                              t2=i
                              p_mu=-pp+two*(t1*p1+t2*p2)
                          else
                              t1=i
                              t2=j
                              p_mu=pp+two*(t1*p1+t2*p2)
                          endif 
                      else	
                          If(kt.lt.zero)then	
                              t1=i
                              t2=j
                              p_mu=pp+two*(t1*p1+t2*p2)
                          endif
                          If(kt.gt.zero)then	
                              t1=j
                              t2=i	
                              p_mu=-pp+two*(t1*p1+t2*p2)			
                          endif
                      endif   
                      If(abs(p-p_mu).lt.1.D-6)goto 310
                  enddo
              enddo
310           continue
              !p0 part!   
              pp_sigma = a_spin*a_spin*mveone_int(muobs,mu,one/AA)/BB+robs*robs*pp
              pp_time = pp_sigma+c_time*pp
              pp_phi = lambda*schwatz_int(muobs,mu,AA)/BB+c_phi*pp
              !******p1 part***********************************
              IF(t1 .eq. 0)THEN
                  p1_sigma=zero
                  p1_time=zero
                  p1_phi=zero
              ELSE
                  IF(PI1_time .eq. zero .or. PI1_sigma .eq. zero)THEN
                      PI1_sigma = a_spin*a_spin*mveone_int(muobs,mu_tp1,one/AA)/BB+robs*robs*PI1
                      PI1_time = PI1_sigma+c_time*PI1 
                  ENDIF
                  IF(PI1_phi .eq. zero)THEN
                      PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB+c_phi*PI1
                  ENDIF
                  p1_sigma = PI1_sigma-pp_sigma
                  p1_time = PI1_time-pp_time
                  p1_phi = PI1_phi-pp_phi
              ENDIF  
              !******p2 part***********************************
              IF(t2 .eq. 0)THEN
                  p2_sigma=zero
                  p2_time=zero
                  p2_phi=zero
              ELSE
                  IF(PI2_time .eq. zero .or. PI2_sigma .eq. zero)THEN
                      PI2_sigma = a_spin*a_spin*mveone_int(mu_tp2,muobs,one/AA)/BB+robs*robs*PI2
                      PI2_time = PI2_sigma+c_time*PI2 
                  ENDIF
                  IF(PI2_phi .eq. zero)THEN
                      PI2_phi = lambda*schwatz_int(mu_tp2,muobs,AA)/BB+c_phi*PI2
                  ENDIF
                  p2_sigma = PI2_sigma+pp_sigma
                  p2_time = PI2_time+pp_time
                  p2_phi = PI2_phi+pp_phi
              ENDIF
              !**********************************************  
              If(mobseqmtp)then
                  If(muobs.eq.mu_tp1)then  
!      equations (52) of Yang and Wang (2013).	
                      sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = -pp_time+two*(t1*p1_time+t2*p2_time)	
                      phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi)	
                  else
!      equations (52) of Yang and Wang (2013).	
                      sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = pp_time+two*(t1*p1_time+t2*p2_time)	
                      phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi)		
                  endif  
              else
                 If(kt.lt.zero)then
!      equations (52) of Yang and Wang (2013).	
                      sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = pp_time+two*(t1*p1_time+t2*p2_time)	
                      phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi)	
                 endif 	
                 If(kt.gt.zero)then
!      equations (52) of Yang and Wang (2013).	
                      sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = -pp_time+two*(t1*p1_time+t2*p2_time)	
                      phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi)	
                 endif	  
              endif
              If(thetamax.eq.zero.or.thetamax.eq.180.D0)phyt = phyt+(t1+t2)*PI
            ELSE 
              mucos=muobs 
              t1 = 0
              t2 = 0
              phyt = zero 
              timet = zero
              sigmat = zero
            ENDIF
          ELSE
              count_num = 1
              goto 30
          ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SPINZERO

!********************************************************************************************
      FUNCTION mveone_int(x,y,z0)
!******************************************************************************************** 
!*     PURPOSE:   To compute the integral \int^y_x z^2/sqrt(1-(z/z0)^2)dz, and z0<1. 
!*     INPUTS:    x,y-----------the integral limits.    
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ******************************************  
      USE constants
      implicit none
      double precision x,y,z0,mveone_int,xt,yt,Iy,Ix

      xt = x
      yt = y
      If(xt .eq. yt)THEN
          mveone_int = zero
          return
      ENDIF 
      Iy = z0**three*half*(dasin(yt/z0)-y/z0/z0*dsqrt(z0*z0-y*y))
      Ix = z0**three*half*(dasin(xt/z0)-x/z0/z0*dsqrt(z0*z0-x*x))
      mveone_int = Iy-Ix  
      return
      END FUNCTION mveone_int

!*****************************************************************************************************
      Subroutine lambdaq_sphericalm(r_sm,a_spin,lambda,q,theta_min)
!******************************************************************************************** 
!*     PURPOSE:   To compute the constants of motion for the spherical motion of a photon. 
!*     INPUTS:    r_sm--------------the radius of the spherical motion.
!*                theta_max---------the maximum or minimum value the \theta coordinate of the 
!*                                  spherical motion, which also the turning point of the motion
!*                                  in theta coordinate.  
!*                signs-------------which determine the direction of the motion of a photon with respect 
!*                                  to the black hole. signs>0 is for co-rotating/prograde orbiting,
!*                                  signs< is for counter-rotating/retrograde orbiting.
!*                a_spin------------the black hole spin.
!*     OUTPUTS:   lamda,q-----------constants of motion.    
!*                theta_min---------the minimum value of the theta coordinate the orbiting.                         
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ******************************************  
      implicit none
      Double precision :: r_sm,a_spin,zetaw,epsilonw,qw,qw1,zero,one,two,three,four,&
             sinthemax,costhemax,cotmax,dtor,AL,AE,DD,Sigma,Delta,sqrts,signs,lambda,q,mve,ep,&
             a1,b1,c1,mu_starp,mu_starm,r_max,r_min,c_temp,mutp2,sintp2,theta_min
      parameter(dtor=asin(1.D0)/90.D0,zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)
	
      r_max = radiusofsphericalmotion(a_spin,zero)
      c_temp = 90.D0
      r_min = radiusofsphericalmotion(a_spin,c_temp)
      If(a_spin.lt.zero)then
          c_temp = r_max
          r_max = r_min
          r_min = c_temp
      endif
      write(*,*)'r_min=',r_min,'   r_max=',r_max
      
      If (a_spin.ne.zero) then   
          If (r_sm.lt.r_min .or. r_sm.gt.r_max .or. r_sm-r_min.le.-1.D-6 .or. r_sm-r_max.ge.1.D-6) then
              write(*,*)'lambdaq_sphericalm(): For spin=',a_spin,' the allowed range for radial coordinate '&
              ,'of the spherical motion of a photon is between ',r_min,' and ',r_max,'. The radius you input is '&
              ,'out of the range and the code shall stop.'
              stop
          endif
      else
          If(r_sm.ne.three)then
              write(*,*)'For a=0, the radius of the spherical motion of a photon is 3 for all inclination anges of '&
                       ,'the orbit with respect to the equatorial plane.'
              r_sm = three
          endif
      endif
!*******************************************************************
      If (a_spin.ne.zero) then
          a1 = a_spin**four*(r_sm-one)**two
          b1 = two*a_spin**two*r_sm*(two*a_spin*a_spin-three*r_sm+r_sm**three)
          c1 = r_sm**three*(-four*a_spin*a_spin+r_sm*(r_sm-three)**two)
    
          mutp2 = (-b1+dsqrt(b1*b1-four*a1*c1))/two/a1 
          sintp2 = one-mutp2
          mu_starm = (-b1-dsqrt(b1*b1-four*a1*c1))/two/a1  
          theta_min = dacos(dsqrt(mutp2))/dtor
          write(*,*)'Theta_min is:',dacos(dsqrt(mutp2))/dtor 
      else
          write(*,*)'For a=0, we need you to input the minimum of the theta coordinate:'
          read(unit=5,fmt=*)theta_min      
          If(theta_min.ne.90.D0)then
              sintp2=sin(theta_min*dtor)**two
              mutp2=cos(theta_min*dtor)**two
          else
              sintp2=one
              mutp2=zero
          endif
      endif
!*******************************************************************
      cotmax=costhemax/sinthemax
      Sigma=r_sm**two+(a_spin)**two*mutp2
      Delta=r_sm**two-two*r_sm+a_spin**two

      AL=Delta*(three*r_sm**four+(a_spin*r_sm)**two*&
                    (mutp2+one)-a_spin**four*mutp2)
      AE=Delta*(r_sm**two-(a_spin)**two*mutp2)

      lambda = dsign(one,a_spin)*dsqrt(sintp2*AL/AE)
      q = mutp2*(AL/AE-a_spin*a_spin) 
      return	 	
      End Subroutine lambdaq_sphericalm

!********************************************************************************
      function Denomenator_r(r,theta,a_spin)
!********************************************************************************
      use constants
      implicit none
      double precision r,theta,a_spin,Denomenator_r,sinthe,costhe

      If(theta.ne.90.D0)then
          sinthe=dsin(theta*dtor)
          costhe=dcos(theta*dtor)
      else
          sinthe=one
          costhe=zero
      endif
      Denomenator_r = ( -(  -three*r*r+(r+one)*a_spin*a_spin*costhe**two+two*a_spin*sinthe*& 
               dsqrt( r*(r*r-a_spin*a_spin*costhe**two) )   ) )**(one/three) 
      return
!********************************************************************************
      end function

!********************************************************************************
      function radiusofsphericalmotion(a,c)
!********************************************************************************
      implicit none
      double precision a,c,radiusofsphericalmotion,r1,r2,Dr

      r1 = 30.D0
      r2 = Denomenator_r(r1,c,a)
      Dr = r2-r1
      r1 = r2
      do while(dabs(Dr).ge.1.D-10)
          r2 = Denomenator_r(r1,c,a)
          Dr = r2-r1
          r1 = r2
      enddo
      radiusofsphericalmotion = r1
      return
      end function
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      end module sphericalmotion

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

 

