!********************************************************************************************
      module pemfinding
!*******************************************************************************
!*     PURPOSE: This module aims on solving more general equation f(p)=0. For 
!*              detail definition of f(p), cf. Yang & Wang (2012).     
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  15 Jan 2012. 
!*******************************************************************************        
      use constants
      use parameters
      use rootsfinding
      use ellfunction
      USE BLcoordinate                
      implicit none

      contains
!******************************************************************************
      function pemfind(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&                              
                rin,rout,muup,mudown,phy1,phy2,caserange,Fp,paras,orir,oricosth,bisection)   
!******************************************************************************
!*     PURPOSE:  Searches for minimum root pem of equations f(p)=0.  s   
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_\phi, p_0. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               rin,rout-------Inner and outer radius of emission region or emission objects.
!*               muup,mudown----Boundary \mu coordinates of emission region or objects.
!*                              Or the maximum and minimum of \mu coordinates of emission region or objects.
!*                              These parameter can be specified void value by setting
!*                              values of parameter caserange. 
!*               phy1,phy2------Boundary \phi coordinates of emission region or objects.
!*                              Or the maximum and minimum of \phi coordinates of emission region or objects.
!*                              These parameter can be specified void value by setting
!*                              values of parameter caserange. 
!*               caserange------Tell routine whether muup, mudown and phy1, phy2 are provided.
!*                              caserange = 1, muup, mudown and phy1, phy2 are provided.
!*                              caserange = 2, muup, mudown are provided, but phy1, phy2 not.
!*                              caserange = 3, muup, mudown are not provided, phy1, phy2 are provided.
!*                              caserange = 4, muup, mudown and phy1, phy2 not are provided. 
!*                              Provided means corresponding parameters have been set specific value.
!*                              Not provided means corresponding parameters have not been set specific value,
!*                              but these parameter should also be given as dummy parameter.
!*               Fp-------------name of function f(p). This routine to compute f(p) should be prvided by
!*                              user, and the dummy variable of Fp should have following form:
!*                              Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras).
!*               paras(1:10)----Array of parameters to descirbe function f(p).  
!*               bisection------Logical variable, if TURE, then use bisection method to search the 
!*                              root of equation f(p)=0, else use Newton-Raphson method. 
!*               NN-------------In the subroutine pemfind there is an important parameter: NN, which is the number
!*                              of sections of the interval (p1 , p2 ) or (p3 , p4 ) has been divided when 
!*                              searching the roots. One needs to set a proper value for NN, since if 
!*                              NN is too small, the roots exit on the interval (p1 , p2 ) or (p3 , p4 ) 
!*                              maybe omitted, if NN is too big, the time consumed by the code will be
!*                              large.
!*     OUTPUTS:  pemfind--------value of root of equation f(p)=0 for p.  
!*                              pemdisk=-1.D0, if the photon goto infinity.
!*                              pemdisk=-2.D0, if the photon fall into event horizon.         
!*     REMARKS:  This routine will search root between interval (p1, p2). We will chose NN points on 
!*               this interval, and check one by one to wheter f(p) changing its sign, if so, the root
!*               must be on interval (p_{i-1}, p_{i}). Then we use Bisection or Newton-Raphson method 
!*               to find the roots. One should set NN propriately to guarantee no root missing and
!*               also quickly find the root, thus saving the running time of the code.
!*     ROUTINES CALLED: radiustp, r2p, rootfind, Sectionp.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        use BLcoordinate
        implicit none
        Double precision pemfind,a,B,lambda,q,sinobs,muobs,a_spin,&
                        robs,scal,rhorizon,NN,r_tp1,r_tp2,mu_tp,mu_tp2,&
                        deltax,p_rout,p_rout2,&
                        p_rin,phya1,phya2,p1,p2,p,f_p,paras(10),f1234(4),&
                        r1,r2
        Double precision :: rin,rout,muup,mudown,phy1,phy2, orir(N+1), oricosth(N+1)
        optional rin,rout,muup,mudown,phy1,phy2
        complex*16  bb(1:4)        
        integer  i,j,k,reals,cases,cases_of_tp,NNf,caserange,tr1,tr2,tm1,tm2,t1,t2
        Double precision,external :: Fp!,Bisection,rootfind,Sectionp
        parameter(deltax=5.D-5)
        logical :: err,mobseqmtp,clines,robs_eq_rtp,indrhorizon,bisection
        character varble 

        rhorizon=one+sqrt(one-a_spin**two)        
        call radiustp(f1234(1),a_spin,robs,lambda,q,r_tp1,r_tp2,&
                reals,robs_eq_rtp,indrhorizon,cases_of_tp,bb)
!NN is important paramerter in the searching of roods by Bisection or Newton-Raphson algorithm.
!Which divides the internal [p1,p2] or [p3,p4] into NN parts. If NN is too small then one may
!miss the roots, while if NN is too big, the speed of code will very slow. Thus one should to
!chose an apropriate one. 
        IF(r_tp1.ge. 10.D0)THEN
            NN=100.D0
        ELSE
            NN=100.D0
        ENDIF 
! To classify the geodesic according its relationship with a shell.
! case 1, 2, 3, 4 correspond to A, B, C, D discussed in Yang & Wang (2012).
        If(present(rout).and.present(rin))then
            If(rin.gt.rhorizon)then
                If(r_tp1.ge.rout)then 
                    cases=1 
                    pemfind=-one 
                else
                    If(r_tp1.gt.rin)then
                        cases=2
                    else
                        If(r_tp1.gt.rhorizon)then
                            cases=3
                        else
                            cases=4        
                        endif
                    endif        
                endif
            else
                cases=5 
            endif
        else
            cases=5
        endif
        !write(*,*)'cases=',cases,r_tp1,rout,rin,f1234(2)
        select case(cases)
        case(1)               
        case(2)
            tr1=0
            tr2=0
            p1=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)
            tr1=1        
            p2=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)   
            paras(1) = mucos(p1,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)
!write(*,*)'p1,p2=',p1,p2 
            If(caserange.eq.4)then        
                pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                        robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)        
            else        
                pemfind=Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        tr1,tr2,p1,p2,muup,mudown,phy1,phy2,caserange,NN,Fp,paras,orir,oricosth,bisection)
            endif 
        case(3)
            tr1=0
            tr2=0        
            p1=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)        
            p2=r2p(f1234(1),rin,lambda,q,a_spin,robs,scal,tr1,tr2) 
            paras(1) = mucos(p1,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)     
            If(caserange.eq.4)then        
                pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)        
            else        
                pemfind=Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        tr1,tr2,p1,p2,muup,mudown,phy1,phy2,caserange,NN,Fp,paras,orir,oricosth,bisection)
                !write(*,*)'ff=',pemfind
            endif
            If(pemfind.eq.-one)then
                tr1=1
                p1=r2p(f1234(1),rin,lambda,q,a_spin,robs,scal,tr1,tr2)        
                p2=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)  
                p2 = p2!+1.D-2	 
                paras(1) = mucos(p1,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)  
                If(caserange.eq.4)then        
                    pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)        
                else        
                    pemfind=Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,tr1,tr2,p1,&
                               p2,muup,mudown,phy1,phy2,caserange,NN,Fp,paras,orir,oricosth,bisection)
                    !write(*,*)'ff=',p1,p2,pemfind
                endif                
            endif                        
        case(4)        
            tr1=0
            tr2=0
            p1=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)        
            p2=r2p(f1234(1),rin,lambda,q,a_spin,robs,scal,tr1,tr2) 
            paras(1) = mucos(p1,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal) 
            If(caserange.eq.4)then        
                pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)        
            else        
                pemfind=Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,tr1,tr2,p1,p2,&
                          muup,mudown,phy1,phy2,caserange,NN,Fp,paras,orir,oricosth,bisection)
            endif
!the photon will fall into the black hole.
            If(pemfind.eq.-one)then 
                pemfind=-two
            endif        
        case(5)
            tr1=0        
            tr2=0           
            p1=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)        
            If(r_tp1.gt.rhorizon)then
                tr1=1
                p2=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)
            else
                p2=r2p(f1234(1),rhorizon,lambda,q,a_spin,robs,scal,tr1,tr2)
            endif
            paras(1) = mucos(p1,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)
            pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                           robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)
        end select
        return        
      end function pemfind
!*****************************************************************************************************
      Function Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                t1,t2,p1,p2,muup,mudown,phy1,phy2,caserange,NN,Fp,paras,orir,oricosth,bisection)  
!*****************************************************************************************************
!*     PURPOSE:  To judge whether a geodesic intersects with the surface of object or emission region
!*               which are described by function f(p). The space range of surface of object or
!*               emission region is (rin, rout) in radius, (mudown,muup) in poloidal and 
!*               (phy1, phy2) in azimuthal. If geodesic has no intersections on those intervals
!*               then it will no intersects with the surface and emission region, a special       
!*               value -1.D0 will be returned.
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t which are the four momentum of the photon. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               t1,t2----------Number of photon meets turning points r_tp1 and r_tp2 respectively.
!*               p1,p2----------roots may hidden inside interval (p1,p2).
!*               muup,mudown,phy1,phy2,
!*               caserange,Fp,paras,bisection--------see instructions in pemfind.  
!*
!*     OUTPUTS:  Sectionp-------value of root of equation f(p)=0 for p.  
!*                              If Sectionp=-1.D0, No roots were found no interval (p1, p2).     
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        use BLcoordinate
        implicit none
        Double precision Sectionp,a,B,lambda,q,sinobs,muobs,a_spin,robs,&
                        scal,rhorizon,NN,mu1,mu2,paras(10),&
                        deltax,p_rout,p_rout2, orir(N+1), oricosth(N+1),&
                        p_rin,phya1,phya2,p1,p2,p,f_p,f1234(4)
        Double precision ,optional :: muup,mudown,phy1,phy2
        integer  i,j,k,t11,t21,t12,t22,tt,reals,cases,NNf,caserange,t1,t2
        Double precision,external :: Fp!,Bisection,rootfind
        parameter(deltax=5.D-5)
        logical :: err,mobseqmtp,clines,bisection

            p1=p1-deltax
            p2=p2+deltax
            If(caserange.eq.1 .or.caserange.eq.3)then        
                phya1=phi(p1,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)
                phya2=phi(p2,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)
            endif                
            If(caserange.eq.1 .or.caserange.eq.2)then        
                mu1=mucos(p1,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)
                mu2=mucos(p2,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal) 
                If((mu1-muup)*(mu1-mudown).lt.zero.or.(mu2-muup)*(mu2-mudown).lt.zero.or.&
                (muup-mu1)*(muup-mu2).lt.zero.or.(mudown-mu1)*(mudown-mu2).lt.zero)then
                    If(caserange.eq.1 .or.caserange.eq.3)then 
                        If((phya1-phy1)*(phya1-phy2).lt.zero.or.(phya2-phy1)*(phya2-phy2).lt.zero.or.&
                        (phy1-phya1)*(phy1-phya2).lt.zero.or.(phy2-phya1)*(phy2-phya2).lt.zero)then
! the geodesic intersecte the zone defined by r1,r2,mu1,mu2,phy1,phy2,so it has the 
!possibility to hit the surface of the object.        
                            Sectionp=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                        robs,scal,t1,t2,p1,p2,NN,Fp,paras,orir,oricosth,bisection) !(1,1)
!write(*,*)'ff=',p1,p2,Sectionp
                        else 
!the (phy1,phy2) and (phya1,phya2) has no public point,so we don't consider it.
                            Sectionp=-one  !(1,1)                
                        endif
                    else
                        Sectionp=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,t1,t2,p1,p2,NN,Fp,paras,orir,oricosth,bisection) !(1,0)
                    endif
                else 
!the internal of (mu1,mu2) and (muup,mudown) does not overfold each other,so the geodesic will not hit the
!surface of the object at the internal (p_rout,p_rout2),which also means the geodesic will never hit the object.
                    Sectionp=-one          ! so nothing needed to be done.                        
                endif
            else
                If(caserange.eq.1 .or.caserange.eq.3)then
                    If((phya1-phy1)*(phya1-phy2).lt.zero.or.(phya2-phy1)*(phya2-phy2).lt.zero.or.&
                    (phy1-phya1)*(phy1-phya2).lt.zero.or.(phy2-phya1)*(phy2-phya2).lt.zero)then
                        Sectionp=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,p1,p2,NN,Fp,paras,orir,oricosth,bisection) !(0,1)
                    else 
!the (phy1,phy2) and (phya1,phya2) has no public points,so we don't consider it.
                        Sectionp=-one  !(0,1)                
                    endif
                else              
                    Sectionp=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,t1,t2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)        ! (0,0)        
                endif        
            endif
        return
      End Function Sectionp


!********************************************************************************* 
      Function rootfind(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                           t1,t2,p1,p2,NN,Fp,paras,orir,oricosth,bisection)
!********************************************************************************* 
!*     PURPOSE:  To search roots on interval (p1, p2). If no roots were found 
!*               a special value -1.D0 will return. 
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t which are the four momentum of the photon. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               t1,t2----------Number of photon meets turning points r_tp1 and r_tp2 respectively.
!*               p1,p2----------roots may hidden inside interval (p1,p2). 
!*               Fp,paras,bisection--------see instructions in pemfind.  
!*
!*     OUTPUTS:  rootfind-------value of root of equation f(p)=0 for p.  
!*                              If rootfind=-1.D0, No roots were found no interval (p1, p2).     
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        use BLcoordinate
        implicit none
        Double precision :: rootfind,a,B,lambda,q,sinobs,muobs,a_spin,&
                        robs,scal,p1,p2,NN,sp1,sp2,f_p_old,&
                        f_p,p,paras(10),f1234(4),deltap,dp, orir(N+1), oricosth(N+1)
        Double precision ,external :: Fp
        parameter (dp=1.D-5)
        integer NNf,k,t1,t2
        logical :: bisection
 
        !p1 = p1 + dp
        deltap=(p2-p1)/NN
        NNf=floor(NN)
        p=p1
        f_p=Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)
        !write(*,*)'NNf=',f_p
        !write(unit=6,fmt=*)p1,p2,f_p
        If(f_p.eq.0.D0)then
                rootfind=p
                return
        endif 
        If(f_p.lt.zero)then
                k=0
                Do while(.true.)			
                        If(f_p.eq.zero)then
                            rootfind=p
                            return
                        endif 
                        If(f_p.gt.zero .or.k.gt.NNf)exit 
                        k=k+1
                        p=p1+deltap*k
                        f_p=Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)
                Enddo                                        
        else
                k=0        
                Do while(.true.)
                        If(f_p.eq.zero)then
                            rootfind=p
                            return
                        endif
                        If(f_p.lt.zero.or.k.gt.NNf)exit 
                        k=k+1
                        p=deltap*k+p1
                        f_p=Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)        
                Enddo
        endif
 !write(unit=6,fmt=*)'f_p=',f_p,k,deltap,p-deltap,p
                If(k.le.NNf)then
                    sp1=p-deltap
                    sp2=p      
! Using bisection or Newton Raphson method to find roots on interval (sp1, sp2).   
                    IF(bisection)THEN
                        rootfind=Bisectionp(f1234,lambda,q,sinobs,muobs,&
                                a_spin,robs,scal,t1,t2,sp1,sp2,Fp,paras,orir,oricosth)
                    else
                        rootfind=NewRapson(f1234,lambda,q,sinobs,muobs,&
                                a_spin,robs,scal,t1,t2,sp1,sp2,Fp,paras,orir,oricosth) 
                    endif
                else
!In (p1,p2) no roots were found!
                        rootfind=-1.D0        
                endif
        return
      End Function rootfind
!*****************************************************************************************************
      Function Bisectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,p1,p2,Fp,paras,orir,oricosth)
!*****************************************************************************************************
!*     PURPOSE:  To search roots on interval (p1, p2) by bisenction method for 
!*               equation Fp(p)=0.
!*     INPUTS:   Parameters descriptions ---------see instructions in pemfind.  
!*
!*     OUTPUTS:  Bisectionp-------value of root of equation f(p)=0 for p.        
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        use BLcoordinate
        implicit none
        Double precision Bisectionp,a,B,lambda,q,p1,p2,sinobs,muobs,&
                        a_spin,robs,scal,pc,f1, orir(N+1), oricosth(N+1),&
                        f2,fc,rout,rin,counter,paras(10),f1234(4)
        Double precision,external :: Fp
        integer t1,t2
        
        pc=(p1+p2)/two
        f1=Fp(p1,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)
        f2=Fp(p2,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)
        fc=Fp(pc,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)        
!write(unit=6,fmt=*)'fc=',f1,f2,fc
        counter=0
        Do while(abs(p2-p1).gt.1.D-4)
             If(f1*fc.gt.zero)then
                p1=pc
                f1=fc
                pc=(p1+p2)/two 
                fc=Fp(pc,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)
             Endif        
             If(f2*fc.gt.zero)then
                p2=pc
                f2=fc
                pc=(p1+p2)/two
                fc=Fp(pc,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)
             Endif
             IF(fc.eq.zero)Then 
                 exit
             Endif
             If(counter.gt.200)exit
!write(unit=6,fmt=*)'fc=',f1,fc,f2,p1,pc,p2
             counter=counter+1                
        Enddo
        Bisectionp=pc        
        return
      End Function Bisectionp 
!*****************************************************************************************************
      Function NewRapson(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,p1,p2,Fp,paras,orir,oricosth)
!*****************************************************************************************************
!*     PURPOSE:  To search roots on interval (p1, p2) by Newton Raphson method for 
!*               equation Fp(p)=0.
!*     INPUTS:   Parameters descriptions ---------see instructions in pemfind.  
!*
!*     OUTPUTS:  NewRapson-------value of root of equation f(p)=0 for p.        
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*****************************************************************
        use BLcoordinate
        implicit none
        Double precision NewRapson,rootfind,lambda,q,sinobs,muobs,a_spin,robs,&
               scal,p1,p2,deltap,NN,sp1,sp2,f_p,f_pj,p,paras(10),f1234(4),dfp,&
               dp,pacc,EPS,temp,h,pj,ptn, orir(N+1), oricosth(N+1)
        Double precision ,external :: Fp
        integer NNf,k,kMax,t1,t2
        parameter (kmax=30, deltap=1.D-6, pacc=1.D-2, EPS=1.D-5)

        ptn=(p2+p1)/two
        Do k=1,kmax
            temp=ptn
            h=EPS*abs(temp)
            if(h.eq.zero)h=EPS        
            pj=temp+h
            h=pj-temp
            f_p=Fp(temp,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)
            f_pj=Fp(pj,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras,orir,oricosth)
            dfp=(f_pj-f_p)/h        
            dp=f_p/dfp
            ptn=ptn-dp
!If((ptn-p1)*(p2-ptn).lt.zero)then
!write(unit=6,fmt=*)'ptn jumps out of brakets [p1, p2]!'
            If(abs(dp).lt.pacc)then
                NewRapson=ptn
                return
            endif  
        Enddo
      End Function NewRapson
!******************************************************* 
      end module pemfinding 


