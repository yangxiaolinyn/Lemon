

!********************************************************************************************
      Function radius(p,f1234r,lambda,q,a_spin,robs,scal)
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
!*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision radius,p,a_spin,rhorizon,q,lambda,scal,zero,integ4(4),&
             bc,cc,dc,ec,b0,b1,b2,b3,g2,g3,tinf,tinf1,PI0,PI1,robs,cr,dr,integ04(4),&
             u,v,w,L1,L2,thorizon,m2,pinf,sn,cn,dn,a4,b4,one,two,four,PI2,ttp,sqt3,&
             integ14(4),three,six,nine,r_tp1,r_tp2,f1234r,tp2,tp,t_inf,PI0_total,&
             PI0_inf_obs,PI0_obs_hori,PI01,PI0_total_2,rff_p
      Double precision f1234r_1,lambda_1,q_1,p_1,a_spin_1,robs_1,scal_1
      parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0,six=6.D0,nine=9.D0)
      complex*16 bb(1:4),dd(3)
      integer ::  reals,i,p4,cases_int,del,index_p4(4),cases,count_num=1
      logical :: robs_eq_rtp,indrhorizon
      save  f1234r_1,lambda_1,q_1,a_spin_1,robs_1,scal_1,r_tp1,r_tp2,reals,&
                robs_eq_rtp,indrhorizon,cases,bb,rhorizon,b0,b1,b2,b3,g2,g3,dd,del,cc,tinf,tp2,&
                thorizon,tinf1,PI0,u,w,v,L1,L2,m2,t_inf,pinf,a4,b4,PI0_total,PI0_inf_obs,PI0_obs_hori,&
                PI0_total_2

 20   continue
      If(count_num.eq.1)then
          f1234r_1=f1234r 
          lambda_1=lambda
          q_1=q 
          a_spin_1=a_spin
          robs_1=robs
          scal_1=scal          
    !*********************************************************************************************          
          rhorizon=one+sqrt(one-a_spin**2)
          a4=zero
          b4=one
          cc=a_spin**2-lambda**2-q
          robs_eq_rtp=.false.
          indrhorizon=.false.
          call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
                             reals,robs_eq_rtp,indrhorizon,cases,bb)
          write(*,*)'hererererer = = = ', p
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
              call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)          
              PI0=integ04(1)
              select case(cases)
              case(1)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then  
                        call weierstrass_int_j3(tinf1,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)        
                        PI0_total=PI0+integ14(1)
                        If(p.lt.PI0_total)then   
! equation (41) in Yang & Wang (2012).                             
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                                  
                        else
                            radius=infinity  !Goto infinity, far away.
                        endif        
                    else
                        call weierstrass_int_J3(tinf1,tinf,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        PI0_inf_obs=integ04(1)
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif                
                    endif
                else
                    If(f1234r.lt.zero)then
                        call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)        
                        PI0_obs_hori=integ04(1) 
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif
                    else
                        call weierstrass_int_J3(tinf1,tinf,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        PI0_inf_obs=integ04(1)        
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                             radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif        
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
                else        
                    If(f1234r.le.zero)then
                        call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)        
                        PI0_obs_hori = integ14(1) 
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif                        
                    else
                        call weierstrass_int_J3(tp2,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)                        
                        call weierstrass_int_J3(tp2,tinf,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                        PI0_total_2=integ14(1)+integ4(1)
                        If(p.lt.PI0_total_2)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif        
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
                    else
                        If(p.lt.one/robs)then
                            radius=robs/(one-robs*p)
                        else
                            radius=infinity          
                        endif                        
                    endif        
                endif
                If(cc.eq.-27.D0)then
                        sqt3=sqrt(three)        
                    If(f1234r.lt.zero)then
                        cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                        dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                        If(p.ne.zero)then        
                            radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=robs!infinity
                        endif
                    else        
                        cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                        dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                        PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                        -Log(one+two/sqt3)/three/sqt3
                        If(p.lt.PI0)then        
                            radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=infinity
                        endif                        
                    endif                
                endif
              endif        
          else
            u=real(bb(4))
            w=abs(aimag(bb(4)))
            v=abs(aimag(bb(2)))
            If(u.ne.zero)then
! equation (45) in Yang & Wang (2012). 
                L1=(four*u**2+w**2+v**2+sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                L2=(four*u**2+w**2+v**2-sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012). 
                thorizon=sqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/(L1-one))/sqrt((rhorizon-u)**2+w**2)
! equation (48) in Yang & Wang (2012). 
                m2=(L1-L2)/L1
                tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/sqrt((robs-u)**two+w**two)
                t_inf=sqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012). 
                pinf=EllipticF(tinf,m2)/w/sqrt(L1)
                call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn)
                If(f1234r.lt.zero)then
                    PI0=pinf-EllipticF(thorizon,m2)/(w*sqrt(L1))        
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
                        radius=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))
                    else
                        radius=rhorizon
                    endif                    
                else
                    PI0=EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
                        radius=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))        
                    else
                        radius=infinity
                    endif
                endif
            else
                If(f1234r.lt.zero)then
                    if(p.lt.(atan(robs/w)-atan(rhorizon/w))/w)then
                        radius=w*tan(atan(robs/w)-p*w)        
                    else
                        radius=rhorizon
                    endif
                else
                    if(p.lt.(PI/two-atan(robs/w))/w)then
                        radius=w*tan(atan(robs/w)+p*w)        
                    else
                        radius=infinity
                    endif                
                endif
            endif                        
          endif
          count_num=count_num+1
      else
        If(f1234r.eq.f1234r_1.and.lambda.eq.lambda_1.and.q.eq.q_1.and.&
        a_spin.eq.a_spin_1.and.robs.eq.robs_1.and.scal.eq.scal_1)then
     !***************************************************************************************************
        write(*,*)'I am here, how are you dears. p =', p, reals, cases
        If(reals.ne.0)then        
            index_p4(1)=0        
            cases_int=1
            select case(cases)
            case(1)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then  
                        If(p.lt.PI0_total)then 
! equation (41) in Yang & Wang (2012).                                
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                                  
                        else
                            radius=infinity  !Goto infinity, far away.
                        endif                
                    else
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif                        
                    endif
                else
                    If(f1234r.lt.zero)then 
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif                                
                    else
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif        
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
                else        
                    If(f1234r.le.zero)then         
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif                        
                    else
                        If(p.lt.PI0_total_2)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif        
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
                    else
                        If(p.lt.one/robs)then
                            radius=robs/(one-robs*p)
                        else
                            radius=infinity          
                        endif                        
                    endif        
                endif
                If(cc.eq.-27.D0)then
                        sqt3=sqrt(three)        
                    If(f1234r.lt.zero)then
                        cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                        dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                        If(p.ne.zero)then        
                            radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=robs!infinity
                        endif
                    else        
                        cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                        dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                        PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                        -Log(one+two/sqt3)/three/sqt3
                        If(p.lt.PI0)then        
                            radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=infinity
                        endif                        
                    endif                
                endif
            endif        
        else
            If(u.ne.zero)then
                call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn)
                If(f1234r.lt.zero)then        
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
                        radius=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))
                    else
                        radius=rhorizon
                    endif                    
                else
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
                        radius=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))        
                    else
                        radius=infinity
                    endif
                endif
            else
                If(f1234r.lt.zero)then
                    if(p.lt.(atan(robs/w)-atan(rhorizon/w))/w)then
                        radius=w*tan(atan(robs/w)-p*w)        
                    else
                        radius=rhorizon
                    endif
                else
                    if(p.lt.(PI/two-atan(robs/w))/w)then
                        radius=w*tan(atan(robs/w)+p*w)        
                    else
                        radius=infinity
                    endif                
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




