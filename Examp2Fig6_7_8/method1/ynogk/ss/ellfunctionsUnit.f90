
     module ellfunction
!********************************************************************
!*    PURPOSE:  This module includes supporting functions and subroutines to compute 
!*              Weierstrass' and Jacobi's elliptical integrals and functions by Carlson's 
!*              integral method. Those codes mainly come from Press (2007) and geokerr.f of
!*              Dexter & Agol (2009).   
!*    AUTHOR:     Yang, Xiao-lin & Wang, Jian-cheng (2012)
!*    DATE WRITTEN:  1 Jan 2012 
!********************************************************************
      use constants
      use rootsfinding
      implicit none

      contains
!***************************************************************************
      Double precision function weierstrassP(z,g2,g3,r1,del)
!***************************************************************************
!*    PURPOSE:   to compute Weierstrass' elliptical function \wp(z;g_2,g_3) and all of 
!*               this function involved are real numbers.   
!*    INPUTS:    z----------the independent variable value.
!*               g_2, g_3---two parameters.
!*               r1(1:3)----an array which is the roots of equation W(t)=4t^3-g_2t-g_3=0.
!*               del--------number of real roots among r1(1:3). 
!*    RETURN:    weierstrassP----the value of function \wp(z;g_2,g_3).  
!*    ROUTINES CALLED:  sncndn
!*    AUTHOR:     Yang, Xiao-lin & Wang, Jian-cheng (2012)
!*    DATE WRITTEN:  1 Jan 2012 
!********************************************************************
      implicit none
      Double precision, intent(in) :: z, g2, g3
      Double precision :: z1,e1,e2,e3,k2,u,sn,alp,bet,sig2,lamb,cn,dn,realp,&
             imagep,rfx,rfy,rfz,EK,two,four,zero, weierstrassP1, &
             weierstrassP2, weierstrassP3, U1, U2
      parameter(two=2.D0,four=4.D0,zero=0.D0)                
      complex*16, intent(in) :: r1(1:3)        
      integer, intent(in) :: del
        
      If(z.eq.zero)then
          weierstrassP=infinity
      else                 
          z1=dabs(z)     
          if(del.eq.3)then
              e1=real(r1(1))
              e2=real(r1(2))
              e3=real(r1(3))
              k2=(e2-e3)/(e1-e3)
              u=z1*dsqrt(e1-e3)
              call sncndn(u,1.D0-k2,sn,cn,dn)
              weierstrassP=e3+(e1-e3)/sn**2
          else
              !alp=-real(r1(1))/two
              alp = real(r1(2))
              bet = dabs(aimag(r1(2)))
              !sig=(9.D0*alp**2+bet**2)**(one/four)
              sig2 = dsqrt( 9.D0*alp**2+bet**2 )
              lamb = (sig2-three*alp)/bet
              !k2=0.5D0+1.5D0*alp/sig**2
              k2 = one / ( one + lamb**2 )
              u = two * dsqrt(sig2) * z1
              call sncndn(u,1.D0-k2,sn,cn,dn)
              if(cn.gt.zero)then                
                  weierstrassP = two * sig2 * ( one + cn ) / sn**2 - two * alp - sig2
                  !write(*,*)'ssssssss = ', weierstrassP, sn
              else
                  weierstrassP = two * sig2 / ( one - cn ) - two * alp - sig2 
                  !If(dabs(sn).gt.1.D-7)then
                  !    weierstrassP=two*sig2*(1.D0+cn)/sn**2-two*alp-sig2
                  !else
                  !    weierstrassP=-two*alp + sig2 * sn**2 / four
                  !endif  
                  !U1 = tempfuncU1( sn, sig2 ) 
                  !U2 = tempfuncU2( U1, sn, sig2 )
                  !weierstrassP1 = tempfuncU1( sn, sig2 ) - two * alp - sig2
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               !   U1 = tempfuncNewton( sn, sig2 )
               !   !write(*,*)'U1 = ', U1, sig2, sn
               !   If (dabs(U1) < sig2) then
               !       !write(*,*)'U11111 = ', U1, sig2
               !       U1 = sig2
               !       stop
               !   Endif
               !   weierstrassP2 = U1 - two * alp - sig2
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  !weierstrassP2 = U2 - two * alp - sig2
                  !write(*,*)'tttt nnnn1 = == ',weierstrassP , sn
                  !write(*,*)'tttt nnnn2 = == ',weierstrassP2, cn
                  !write(*,*)'tttt nnnn3 = == ',weierstrassP2 - weierstrassP
                  !write(*,*)'tttt nnnnd = == ',weierstrassP-weierstrassP1!, weierstrassP- weierstrassP1
                  !If(dabs(weierstrassP - weierstrassP) > 1.D-15)&
                  !    write(*,*)'tttt nnnn3 = == ',dabs(weierstrassPtt - weierstrassP), sn
              endif     
          endif
      endif        
      return
      end function weierstrassP

!*************************************************************************************************
      Function tempfuncU1( t, sig2 )
!*************************************************************************************************
      Double precision, intent(in) :: t, sig2
      Double precision :: tempfuncU1, y1, u1, u2, y2
      integer(kind=8) :: i

      u1 = sig2
      y1 = t * u1
      i = 0
      Do
          u2 = y1**2/sig2/four + sig2
          !write(*,*)'sss = ',u1,u2,dabs(u1 -u2)
          If ( dabs(u2-u1) <= 1.D-26 ) exit
          y1 = t * u2
          u1 = u2
          i = i + 1
      Enddo
      write(*,*)'iissskjjjj = ',i
      tempfuncU1 = u2
      return
      End function tempfuncU1

!*************************************************************************************************
      Function tempfuncU2( u1, t, sig2 )
!*************************************************************************************************
      Double precision, intent(in) :: u1, t, sig2
      Double precision :: tempfuncU2, y1, v1, v2, y2
      integer(kind=8) :: i

      v1 = u1 + 1.D0
      y1 = two * dsqrt( sig2 * ( v1 - sig2 ) )
      i = 0
      Do
          v2 = y1 / t
          !write(*,*)'sss = ',v1,v2,sig2
          If ( dabs(v2-v1) <= 1.D-10 ) exit
          y1 = two * dsqrt( sig2 * ( v2 - sig2 ) )
          v1 = v2
          i = i + 1
          if (i > 100) exit
      Enddo
      tempfuncU2 = v2
      return
      End function tempfuncU2

!*************************************************************************************************
      Function tempfuncNewton( t, sig2 )
!*************************************************************************************************
      Double precision, intent(in) :: t, sig2
      Double precision :: tempfuncNewton, u1, u2, X, d1, d2
      integer(kind=8) :: i

      u1 = zero
      u2 = u1 - ( t**2 * u1**2 - four * sig2 * u1 + four * sig2**2 ) / &
               ( two * t**2 * u1 - four * sig2 )
      !d1 = dabs(u1-u2)
      !d2 = zero
      i = 0
      Do
          !write(*,*)'sssffffssdf === ', d1,d2
          If ( dabs(u1-u2) <= 1.D-14 ) exit
          !d2 = dabs(u1-u2)
          u1 = u2
          u2 = u1 - ( t**2 * u1**2 - four * sig2 * u1 + four * sig2**2 ) / &
                   ( two * t**2 * u1 - four * sig2 )
          !d1 = dabs(u1-u2)
          i = i + 1
          if (i > 50) exit
      Enddo
      !write(*,*)'iiiiiiii = ',i
      tempfuncNewton = u2
      return
      End function tempfuncNewton

!*************************************************************************************************
      Function halfperiodwp(g2,g3,r1,del)
!************************************************************************************************* 
!*    PURPOSE:   to compute the semi period of Weierstrass' elliptical function \wp(z;g_2,g_3) and all of 
!*               this function involved are real numbers.   
!*    INPUTS:    g_2, g_3---two parameters. 
!*               r1(1:3)----an array which is the roots of equation W(t)=4t^3-g_2t-g_3=0.
!*               del--------number of real roots among r1(1:3). 
!*    RETURN:    halfperiodwp----the semi period of function \wp(z;g_2,g_3).  
!*    ROUTINES CALLED:  rf
!*    AUTHOR:    Yang, Xiao-lin & Wang, Jian-cheng (2012)
!*    DATE WRITTEN:  1 Jan 2012 
!********************************************************************************************
      implicit none
      Double precision halfperiodwp,g2,g3,g1,e1,e2,e3,zero,one,two,&
                           three,four,EK,alp,bet,sig2,lamb,k2
      parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0) 
      complex*16 r1(3)      
      integer  i,del
 
      if(del.eq.3)then
         e1=real(r1(1))
         e2=real(r1(2))
         e3=real(r1(3))
         k2=(e2-e3)/(e1-e3)
         EK=rf(zero,one-k2,one)
         halfperiodwp=EK/dsqrt(e1-e3)
      else
         alp=real(r1(2))
         bet=dabs(aimag(r1(2)))
         sig2=dsqrt(9.D0*alp**two+bet**two)
         k2=one/two+three/two*alp/sig2
         EK=rf(zero,one-k2,one)
         halfperiodwp=EK/dsqrt(sig2)
      endif
      return
      End Function halfperiodwp
!*************************************************************************************************
      subroutine sncndn(uu,emmc,sn,cn,dn)
!************************************************************************************************* 
!*    PURPOSE:   Returns the Jacobian elliptic functions sn(u|k^2), cn(u|k^2), 
!*            and dn(u|k^2). Here uu=u, while emmc=1-k^2. 
!*    RETURN:    sn, cn, dn----Jacobian elliptic functions sn(u|k^2), cn(u|k^2), dn(u|k^2). 
!*    AUTHOR:    Press et al. (2007) 
!********************************************************************
      implicit none      
      Double precision uu,emmc,sn,CA
      Double precision ,optional :: cn,dn
      !parameter (CA=3.D0-8)
      parameter (CA=3.D0-8)
      integer i,ii,l
      Double precision a,b,c,d,emc,u,em(13),en(13)
      logical bo

      emc=emmc
      u=uu
      if(emc .ne. 0.D0)then
          bo=(emc .lt. 0.D0)        
          if(bo)then
              d=1.D0-emc   !t'=t*k, u'=k*u, k'^2=1./k^2,  
              emc= - emc / d
              d=dsqrt(d)
              u=d*u
          end if 
          a=1.D0
          dn=1.D0
          l1: do i=1,13 
              l=i
              em(i)=a
              emc=dsqrt(emc)
              en(i)=emc
              c=0.5D0*(a+emc)
              if(dabs(a-emc).le.CA*a) exit
              emc=a*emc
              a=c
          end do l1
          u=c*u
          sn=dsin(u)  
          cn=dcos(u)                
          if(sn.eq.0.D0)then
              if(bo)then
                  a=dn
                  dn=cn
                  cn=a
                  sn=sn/d
              end if
              return        
          endif 
          a=cn/sn
          c=a*c
          l2: do ii=l,1,-1
              b=em(ii)
              a=c*a
              c=dn*c
              dn=(en(ii)+a)/(b+a)
              a=c/b
          enddo l2
          a=1.D0/dsqrt(c**2+1.D0)
          if(sn.lt.0.D0)then
              sn=-a
          else
              sn=a
          endif
              cn=c*sn
          return           
      else
          cn=1.D0/dcosh(u)
          dn=cn                   
          sn=dtanh(u)
          return
      end if
      end subroutine sncndn
!*************************************************************************************************
      Double precision FUNCTION rf(x,y,z) 
!*************************************************************************************************
!*     PURPOSE: Compute Carlson fundamental integral RF
!*              R_F=1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2)
!*     ARGUMENTS: Symmetric arguments x,y,z
!*     ROUTINES CALLED:  None.
!*     ALGORITHM: Due to B.C. Carlson.
!*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
!*     REMARKS:  
!*     AUTHOR:  Press et al (2007).
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!***********************************************************************
      implicit none
      Double precision x,y,z,ERRTOL,TINY1,BIG,THIRD,C1,C2,C3,C4,delta,zero
      parameter (ERRTOL=0.0025D0,TINY1=1.5D-38,BIG=3.D37,THIRD=1.D0/3.D0,C1=1.D0/24.D0,&
                          C2=0.1D0,C3=3.D0/44.D0,C4=1.D0/14.D0,zero=0.D0)
      Double precision alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      !if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY1.or.max(x,y,z).gt.BIG)then
      !      rf=0.D0
      !      return
      !endif
      IF(x.lt.zero)x=zero
      IF(y.lt.zero)y=zero
      IF(z.lt.zero)z=zero
      xt=x
      yt=y
      zt=z
      sqrtx=dsqrt(xt)
      sqrty=dsqrt(yt)
      sqrtz=dsqrt(zt)
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      xt=0.25D0*(xt+alamb)
      yt=0.25D0*(yt+alamb)
      zt=0.25D0*(zt+alamb)
      ave=THIRD*(xt+yt+zt)
      delx=(ave-xt)/ave
      dely=(ave-yt)/ave
      delz=(ave-zt)/ave      
      delta=max(dabs(delx),dabs(dely),dabs(delz))
      Do while(delta.gt.ERRTOL)
            sqrtx=dsqrt(xt)
            sqrty=dsqrt(yt)
            sqrtz=dsqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            ave=THIRD*(xt+yt+zt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
            delta=max(dabs(delx),dabs(dely),dabs(delz))
      enddo
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.D0+(C1*e2-C2-C3*e3)*e2+C4*e3)/dsqrt(ave)
      return      
      end Function rf
!************************************************************************ 
      Double precision FUNCTION rj(x,y,z,p) 
!************************************************************************ 
!*     PURPOSE: Compute Carlson fundamental integral RJ
!*     RJ(x,y,z,p) = 3/2 \int_0^\infty dt
!*                      (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1)
!*     ARGUMENTS: x,y,z,p
!*     ROUTINES CALLED:  RF, RC.
!*     ALGORITHM: Due to B.C. Carlson.
!*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
!*     REMARKS:  
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!*********************************************************************** 
       implicit none
      Double precision x,y,z,p,ERRTOL,TINY1,BIG,C1,C2,C3,C4,C5,C6,C7,C8,delta,zero
      Double precision a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,fac,pt,&
      rcx,rho,sqrtx,sqrty,sqrtz,sum1,tau,xt,yt,zt
      Parameter(ERRTOL=0.0015D0,TINY1=2.5D-13,BIG=9.D11,C1=3.D0/14.D0,&
                     C2=1.D0/3.D0,C3=3.D0/22.D0,C4=3.D0/26.D0,C5=0.75D0*C3,&
                           C6=1.5D0*C4,C7=0.5D0*C2,C8=C3+C3,zero=0.D0)
      
      !if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY1.or.max(x,y,z,abs(p)).gt.BIG)then
      !   rj=0.D0            
      !   write(*,*)'ffsdfa'
      !   return
      !endif
      IF(x.lt.zero)x=zero
      IF(y.lt.zero)y=zero
      IF(z.lt.zero)z=zero
      sum1=0.D0
      fac=1.D0
      If(p.gt.0.D0)then
            xt=x
            yt=y
            zt=z
            pt=p
      else
            xt=min(x,y,z)
            zt=max(x,y,z)
            yt=x+y+z-xt-zt
            a=1.D0/(yt-p)
            b=a*(zt-yt)*(yt-xt)
            pt=yt+b
            rho=xt*zt/yt
            tau=p*pt/yt
            rcx=rc(rho,tau)
      endif
            sqrtx=dsqrt(xt)
            sqrty=dsqrt(yt)
            sqrtz=dsqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
            beta=pt*(pt+alamb)**2
            sum1=sum1+fac*rc(alpha,beta)
            fac=0.25D0*fac
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            pt=0.25D0*(pt+alamb)
            ave=0.2D0*(xt+yt+zt+pt+pt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
            delp=(ave-pt)/ave
            delta=max(dabs(delx),dabs(dely),dabs(delz),dabs(delp))
      Do while(delta.gt.ERRTOL)
            sqrtx=dsqrt(xt)
            sqrty=dsqrt(yt)
            sqrtz=dsqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
            beta=pt*(pt+alamb)**2
            sum1=sum1+fac*rc(alpha,beta)
            fac=0.25D0*fac
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            pt=0.25D0*(pt+alamb)
            ave=0.2D0*(xt+yt+zt+pt+pt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
            delp=(ave-pt)/ave
            delta=max(dabs(delx),dabs(dely),dabs(delz),dabs(delp))
      enddo  
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3.D0*ec
      ee=eb+2.D0*delp*(ea-ec)
      rj=3.D0*sum1+fac*(1.D0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+&
                     delp*(-C8+delp*C4))+delp*ea*(C2-delp*C3)-&
      C2*delp*ec)/(ave*dsqrt(ave))      
      If(p.le.0.D0)rj=a*(b*rj+3.D0*(rcx-rf(xt,yt,zt)))
      return
      end Function rj
!************************************************************************ 
      FUNCTION rc(x,y)
!************************************************************************
!*     PURPOSE: Compute Carlson degenerate integral RC
!*              R_C(x,y)=1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1)
!*     ARGUMENTS: x,y
!*     ROUTINES CALLED:  None.
!*     ALGORITHM: Due to B.C. Carlson.
!*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
!*     REMARKS:  
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!***********************************************************************
      implicit none
      Double precision rc,x,y,ERRTOL,TINY1,sqrtNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2,C3,C4
      PARAMETER(ERRTOL=0.0012D0,TINY1=1.69D-38,sqrtNY=1.3D-19,BIG=3.D37,&
      TNBG=TINY1*BIG,COMP1=2.236D0/sqrtNY,COMP2=TNBG*TNBG/25.D0,&
      THIRD=1.D0/3.D0,C1=0.3D0,C2=1.D0/7.D0,C3=0.375D0,C4=9.D0/22.D0)
      Double precision alamb,ave,s,w,xt,yt 

      !if((x.lt.0.D0).or.(y.eq.0.D0).or.((x+abs(y)).lt.TINY1).or.((x+abs(y)).gt.BIG).or.&
      !                  ((y.lt.-COMP1).and.(x.gt.0.D0).and.(x.lt.COMP2)))then
      !    rc=0.D0            
      !    return
      !endif
      if(y.gt.0.D0)then
            xt=x
            yt=y
            w=1.D0
      else
            xt=x-y
            yt=-y
            w=dsqrt(x)/dsqrt(xt)
      endif
            alamb=2.D0*dsqrt(xt)*dsqrt(yt)+yt
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            ave=THIRD*(xt+yt+yt)
            s=(yt-ave)/ave
      Do While(dabs(s).gt.ERRTOL)
            alamb=2.D0*dsqrt(xt)*dsqrt(yt)+yt
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            ave=THIRD*(xt+yt+yt)
            s=(yt-ave)/ave
      ENDdo
      rc=w*(1.D0+s*s*(C1+s*(C2+s*(C3+s*C4))))/dsqrt(ave)
      return
      END FUNCTION rc
!**********************************************************************
      FUNCTION rd(x,y,z)
!**********************************************************************
!*     PURPOSE: Compute Carlson degenerate integral RD
!*              R_D(x,y,z)=3/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2)
!*     ARGUMENTS: x,y,z
!*     ROUTINES CALLED:  None.
!*     ALGORITHM: Due to B.C. Carlson.
!*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
!*     REMARKS:  
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!***********************************************************************
      implicit none
      Double precision rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,zero
        PARAMETER (ERRTOL=0.0015D0,TINY=1.D-25,BIG=4.5D21,C1=3.D0/14.D0,C2=1.D0/6.D0,&
                         C3=9.D0/22.D0,C4=3.D0/26.D0,C5=0.25D0*C3,C6=1.5D0*C4,zero=0.D0) 
      Double precision alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt

      !if(min(x,y).lt.0.D0.or.min(x+y,z).lt.TINY.or. max(x,y,z).gt.BIG)then
      !      rd=0.D0      
      !      return
      !endif
      IF(x.lt.zero)x=zero
      IF(y.lt.zero)y=zero      
            xt=x
            yt=y
            zt=z
            sum=0.D0
            fac=1.D0
      
            sqrtx=dsqrt(xt)
            sqrty=dsqrt(yt)
            sqrtz=dsqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            sum=sum+fac/(sqrtz*(zt+alamb))
            fac=0.25D0*fac
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            ave=0.2D0*(xt+yt+3.D0*zt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
      DO While(max(dabs(delx),dabs(dely),dabs(delz)).gt.ERRTOL)
            sqrtx=dsqrt(xt)
            sqrty=dsqrt(yt)
            sqrtz=dsqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            sum=sum+fac/(sqrtz*(zt+alamb))
            fac=0.25D0*fac
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            ave=0.2D0*(xt+yt+3.D0*zt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
      End DO
            ea=delx*dely
            eb=delz*delz
            ec=ea-eb
            ed=ea-6.D0*eb
            ee=ed+ec+ec
            rd=3.D0*sum+fac*(1.D0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+&
                                 delz*(-C3*ec+delz*C4*ea)))/(ave*dsqrt(ave))
            return
      END function rd
!******************************************************************* 
      Function EllipticF(t,k2)
!******************************************************************* 
!*     PURPOSE: calculate Legendre's first kind elliptic integral: 
!*              F(t,k2)=\int_0^t dt/sqrt{(1-t^2)*(1-k2*t^2)}.  
!*     ARGUMENTS: t, k2
!*     ROUTINES CALLED:  RF  
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!***********************************************************************
      implicit none
      Double precision t,k2,EllipticF,x1,y1,z1  
      
      x1=1.D0-t*t
      y1=1.D0-k2*t*t
      z1=1.D0
!Press et al. 2007 (6.12.19)
      EllipticF=t*rf(x1,y1,z1)
      return
      end function EllipticF       
!********************************************************************* 
      Function EllipticE(t,k2)
!********************************************************************* 
!*     PURPOSE: calculate Legendre's second kind elliptic integrals: 
!*              E(t,k2)=\int_0^t sqrt{1-k2*t^2}/sqrt{(1-t^2)}dt.
!*     ARGUMENTS: t, k2
!*     ROUTINES CALLED:  RF, RD 
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!********************************************************************** 
      implicit none
      Double precision t,k2,EllipticE,x1,y1,z1 
      
      x1=1.D0-t*T
      y1=1.D0-k2*t*t
      z1=1.D0
!Press et al. 2007 (6.12.20)
      EllipticE=t*rf(x1,y1,z1)-1.D0/3.D0*k2*t**3*rd(x1,y1,z1)
      return
      end function EllipticE
!*********************************************************************** 
      Function EllipticPI(t,n,k2)
!*********************************************************************** 
!*     PURPOSE: calculate Legendre's third kind elliptic integrals: 
!*              PI(t,n,k2)=\int_0^t /(1+nt^2)/sqrt{(1-k2*t^2)(1-t^2)}dt. 
!*     ARGUMENTS: t, k2
!*     ROUTINES CALLED:  RF, RJ
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!********************************************************************** 
      implicit none
      Double precision t,k2,EllipticPI,x1,y1,z1,w1,n 
      
      x1=1.D0-t*t
      y1=1.D0-k2*t*t
      z1=1.D0
      w1=1.D0+n*t*t
!Press et al. 2007 (6.12.20)
      EllipticPI=t*rf(x1,y1,z1)-1.D0/3.D0*n*t*t*t*rj(x1,y1,z1,w1)
      return
      end function EllipticPI
!*************************************************************************************************
      subroutine weierstrass_int_J3(y,x,bb,del,a4,b4,p4,rff_p,integ,cases)
!*************************************************************************************************
!*     PURPOSE: Computes integrals: J_k(h)=\int^x_y (b4*t+a4)^(k/2)*(4*t^3-g_2*t-g_3)^(-1/2)dt.
!*              Where integer index k can be 0, -2, -4 and 2. (75) and (76) of Yang & Wang (2012).   
!*     INPUTS:  x,y -- limits of integral.
!*              bb(1:3) -- Roots of equation 4*t^3-g_2*t-g_3=0 solved by routine root3.
!*              del -- Number of real roots in bb(1:3).
!*              p4,rff_p,integ -- p4(1:4) is an array which specifies the value of index k of J_k(h). 
!*                 If p4(1)=0, then J_0 was computed and sent to integ(1).
!*                 If p4(1)!=0, then J_0 was replaced by parameter p, and rff_p=p.  
!*                 If p4(2)=-2, then J_{-2} was computed and sent to integ(2).                      
!*                 If p4(3)=2, then J_{2} was computed and sent to integ(3).
!*                 If p4(4)=-4, then J_{-4} was computed and sent to integ(4).
!*              cases -- If cases=1, then only J_0 was computed.
!*                       If cases=2, then only J_0 and J_{-2} are computed.
!*                       If cases=3, then only J_0, J_{-2} and J_{2} are computed.    
!*                       If cases=4, then J_0, J_{-2}, J_{2} and J_{-4} are computed.            
!*     OUTPUTS: integ -- is an array saved the results of J_k(h). 
!*     ROUTINES CALLED:  ellcubicreals, ellcubiccomplexs     
!*     ACCURACY:   Machine.
!*     REMARKS: Based on Yang & Wang (2012).
!*     AUTHOR:     Yang & Wang (2012).
!*     DATE WRITTEN:  4 Jan 2012 
!***********************************************************************
      implicit none
      Double precision y,x,yt,xt,h4,g2,g3,a1,b1,a2,b2,a3,b3,a4,b4,rff_p,integ(4),b,three,two,one,&
                  tempt,f,g,h,a44,b44,sign_h
      parameter  (three=3.0D0,two=2.0D0,one=1.D0)
      integer  del,p4(4),i,cases
      complex*16 bb(1:3)
      logical :: inverse,neg

      xt=x
      yt=y
      a44=a4
      b44=b4
      inverse=.false.
      neg=.false.
      If(dabs(xt-yt).le.1.D-9)then
         integ=0.D0
         return      
      endif
      if(yt.gt.xt)then
         tempt=xt
         xt=yt
         yt=tempt 
         inverse=.true.                  
      endif
      b=0.D0
      sign_h=dsign(one,b44*xt+a44)
!write(*,*)'weierstrass_int_J3=',del
      If(del.eq.3)then
! equation (75) of Yang & Wang (2012).
            a44=sign_h*a44
            b44=sign_h*b44
            a1=-real(bb(1))
            a2=-real(bb(2))
            a3=-real(bb(3))
            b1=1.D0
            b2=1.D0
            b3=1.D0
            call ellcubicreals(p4,a1,b1,a2,b2,a3,b3,a44,b44,yt,xt,rff_p*two,integ,cases)
            if(inverse)then
                  integ(1)=-integ(1)
                  integ(2)=-integ(2)
                  integ(3)=-integ(3)
                  integ(4)=-integ(4)
            endif
            Do i=1,4
                integ(i)=integ(i)/two      
                integ(i)=integ(i)*(sign_h)**(-p4(i)/2)
            Enddo
!write(*,*)'weierstrass_int_J3',integ
      else
! equation (76) of Yang & Wang (2012).
            a44=sign_h*a44
            b44=sign_h*b44
            a1=-real(bb(1))
            b1=one
            f=real(bb(2))**2+aimag(bb(2))**2
            g=-two*real(bb(2))
            h=one
               call ellcubiccomplexs(p4,a1,b1,a44,b44,f,g,h,yt,xt,rff_p*two,integ,cases)
            if(inverse)then
                  integ(1)=-integ(1)
                  integ(2)=-integ(2)
                  integ(3)=-integ(3)
                  integ(4)=-integ(4)
            endif
            Do i=1,4
                integ(i)=integ(i)/two      
                integ(i)=integ(i)*(sign_h)**(-p4(i)/2)
            Enddo
!write(*,*)'weierstrass_int_J3',integ
      endif 
      return
      end subroutine weierstrass_int_J3
!**********************************************************************************************
      subroutine carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,a5,b5,p5,rff_p,integ,cases)
!**********************************************************************************************    
!*     PURPOSE: Computes integrals: J_k(h)=\int^x_y (b5*r+a5)^(k/2)*[(h1*r^2+g1*r+f1)(h2*r^2+g2*r+f2)]^(-1/2)dr.
!*              Where integer index k can be 0, -2, -4, 2 and 4. (77) of Yang & Wang (2012).   
!*     INPUTS:  x,y -- limits of integral.  
!*              p5,rff_p,integ -- p5(1:5) is an array which specifies the value of index k of J_k(h). 
!*                 If p5(1)=0, then J_0 was computed and sent to integ(1).
!*                 If p5(1)!=0, then J_0 was replaced by parameter p, and rff_p=p.  
!*                 If p5(2)=-2, then J_{-2} was computed and sent to integ(2).                      
!*                 If p5(3)=2, then J_{2} was computed and sent to integ(3).
!*                 If p5(4)=-4, then J_{-4} was computed and sent to integ(4).
!*                 If p5(4)=4, then J_{4} was computed and sent to integ(5).
!*              cases -- If cases=1, then only J_0 will be computed.
!*                       If cases=2, then only J_0 and J_{-2} will be computed.
!*                       If cases=3, then only J_0, J_{-2} and J_{2} will be computed.
!*                       If cases=4, then J_0, J_{-2}, J_{2} and J_{-4} will be computed. 
!*                       If cases=5, then J_0, J_{-2}, J_{2} and J_{4} will be computed.     
!*     OUTPUTS: integ -- is an array saved the results of J_k(h). 
!*     ROUTINES CALLED:  elldoublecomplexs
!*     ACCURACY:   Machine.
!*     REMARKS: Based on Yang & Wang (2012).
!*     AUTHOR:     Yang & Wang (2012).
!*     DATE WRITTEN:  4 Jan 2012 
!***********************************************************************
      implicit none
      Double precision y,x,xt,yt,f1,g1,h1,f2,g2,h2,a5,b5,rff,integ(1:5),b,tempt,f,g,h,a55,&
                  b55,sign_h,one,zero,rff_p
      parameter(one=1.D0,zero=0.D0)
      integer  reals,p5(1:5),cases,i
      logical :: inverse,neg
      xt=x
      yt=y
      a55=a5
      b55=b5
      inverse=.false.
      neg=.false.
      If(dabs(xt-yt).le.1.D-9)then
         integ=zero
         return      
      endif
      if(yt.gt.xt)then
         tempt=xt
         xt=yt
         yt=tempt 
         inverse=.true.                  
      endif
      sign_h=dsign(one,b55*xt+a55)      
      a55=sign_h*a55
      b55=sign_h*b55
! equation (77) of Yang & Wang (2012).
         call elldoublecomplexs(p5,f1,g1,h1,f2,g2,h2,a55,b55,yt,xt,rff_p,integ,cases)
      if(inverse)then
            integ(1)=-integ(1)
            integ(2)=-integ(2)
            integ(3)=-integ(3)
            integ(4)=-integ(4)
            integ(5)=-integ(5)
      endif      
      Do i=1,5 
            integ(i)=integ(i)*(sign_h)**(-p5(i)/2)
      Enddo
      return
      end subroutine carlson_doublecomplex5
!*******************************************************************************
      subroutine ellcubicreals(index_p4,a1,b1,a2,b2,a3,b3,a4,b4,y,x,rff_p,integ,cases)
!*******************************************************************************
!*     PURPOSE: Computes J_k(h)=\int_y^x dt (b4*t+a4)^(k/2)[(b1*t+a1)*(b2*t+a2)*(b3*t+a3)]^{-1/2}. 
!*              It is the case of equation W(t)=4*t^3-g_2*t-g_3=0 has three real roots. 
!*              Equation (61) of Yang & Wang (2012).   
!*     INPUTS:  Arguments for above integral. If index_p4(1)=0, then J_0 will be computed, 
!*              else J_0 will be replaced by parameter p, and rff_p=p. 
!*     OUTPUTS:  Value of integral J_k(h).
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.
!*     AUTHOR:     Dexter & Agol (2009)
!*     MODIFIED:   Yang & Wang (2012)
!*     DATE WRITTEN:  4 Mar 2009
!*     REVISIONS: 
!***********************************************************************
      implicit NONE
      Double precision zero,one,half,two,three,ellcubic,d12,d13,d14,d24,d34,X1,X2,X3,X4,&
                       Y1,Y2,Y3,Y4,U1c,U32,U22,W22,U12,Q22,P22,I1c,I3c,r12,r13,r24i,r34i,&
                       I2c,J2c,K2c,a1,b1,a2,b2,a3,b3,a4,b4,y,x,rff_p,r14,r24,r34,rff,integ(4)
      integer  index_p4(4),cases
      PARAMETER ( ZERO=0.D0, ONE=1.D0, TWO=2.D0, HALF=0.5d0, THREE=3.D0 )
      ellcubic=0.d0
      !c (2.1) Carlson (1989)
      d12=a1*b2-a2*b1
      d13=a1*b3-a3*b1
      d14=a1*b4-a4*b1
      d24=a2*b4-a4*b2
      d34=a3*b4-a4*b3
      r14=a1/b1-a4/b4
      r24=a2/b2-a4/b4
      r34=a3/b3-a4/b4            
      !c (2.2) Carlson (1989)
      X1=dsqrt(dabs(a1+b1*x))
      X2=dsqrt(dabs(a2+b2*x))
      X3=dsqrt(dabs(a3+b3*x))
      X4=dsqrt(dabs(a4+b4*x))
      Y1=dsqrt(dabs(a1+b1*y))
      Y2=dsqrt(dabs(a2+b2*y))
      Y3=dsqrt(dabs(a3+b3*y))
      Y4=dsqrt(dabs(a4+b4*y))
      !c! (2.3) Carlson (1989)
      If(x.lt.infinity)then      
          U1c=(X1*Y2*Y3+Y1*X2*X3)/(x-y)
          U12=U1c**2
          U22=((X2*Y1*Y3+Y2*X1*X3)/(x-y))**2
          U32=((X3*Y1*Y2+Y3*X1*X2)/(x-y))**2
      else
          U1c=dsqrt(dabs(b2*b3))*Y1
          U12=U1c**2
          U22=b1*b3*Y2**2
          U32=b2*b1*Y3**2      
      endif             
      !c (2.4) Carlson (1989)
      W22=U12
      W22=U12-b4*d12*d13/d14
      ! (2.5) Carlson (1989)
      If(x.lt.infinity)then            
          Q22=(X4*Y4/X1/Y1)**2*W22 
      else
          Q22=b4/b1*(Y4/Y1)**2*W22       
      endif 
      P22=Q22+b4*d24*d34/d14
       
      !c Now, compute the three integrals we need [-1,-1,-1],[-1,-1,-1,-2], and 
      !c  [-1,-1,-1,-4]:we need to calculate the [-1,-1,-1,2] integral,we add it in this part.
      if(index_p4(1).eq.0) then
        !c (2.21) Carlson (1989)
        rff=rf(U32,U22,U12)
        integ(1)=two*rff
        if(cases.eq.1)return
      else
        rff=rff_p/two
      endif
          !c (2.12) Carlson (1989)
          I1c=rff_p!two*rff
      If(index_p4(3).eq.2) then
            !c (2.13) Carlson (1989)
            I2c=two/three*d12*d13*rd(U32,U22,U12)+two*X1*Y1/U1c
            !  (2.39) Carlson (1989)
            integ(3)=(b4*I2c-d14*I1c)/b1      
            if(cases.eq.3)return
      endif
      If(X1*Y1.ne.zero) then
      !c (2.14) Carlson (1989)
          I3c=two*rc(P22,Q22)-two*b1*d12*d13/three/d14*rj(U32,U22,U12,W22)
      Else
       ! One can read the paragraph between (2.19) and (2.20) of Carlson (1989).
        I3c=-two*b1*d12*d13/three/d14*rj(U32,U22,U12,W22)
      Endif
        if(index_p4(2).eq.-2) then
      !c (2.49) Carlson (1989)
          integ(2)=(b4*I3c-b1*I1c)/d14
        if(cases.eq.2)return
        endif

        If(index_p4(4).eq.-4)then
            !c (2.1)  Carlson (1989)
               r12=a1/b1-a2/b2
                r13=a1/b1-a3/b3
                r24i=b2*b4/(a2*b4-a4*b2)
                r34i=b3*b4/(a3*b4-a4*b3)      
             If(x.lt.infinity)then
            !c (2.17) Carlson (1989)
               J2c=two/three*d12*d13*rd(U32,U22,U12)+two*d13*X2*Y2/X3/Y3/U1c
            !c (2.59) & (2.6) Carlson (1989)
               K2c=b3*J2c-two*d34*(X1*X2/X3/X4**2-Y1*Y2/Y3/Y4**2)
           else
               J2c=two/three*d12*d13*rd(U32,U22,U12)+two*d13*Y2/b3/Y3/Y1
               K2c=b3*J2c+two*d34*Y1*Y2/Y3/Y4**2
           endif      
               !c (2.62) Carlson (1989)
               integ(4)=-I3c*half/d14*(one/r14+one/r24+one/r34)+&
               half*b4/d14/d24/d34*K2c+(b1/d14)**2*(one-half*r12*r13*r24i*r34i)*I1c
                !!write(*,*)I3c,d14,d24,d34,r12,r13,r14,r24,r34,K2c,r24i,r34i,I1c,integ(4)
        endif  
      return
      end  subroutine ellcubicreals
!*******************************************************************************
      subroutine ellcubiccomplexs(index_p4,a1,b1,a4,b4,f,g,h,y,x,rff_p,integ,cases)
!*******************************************************************************
!*     PURPOSE: Computes J_k(h)=\int_y^x dt (b4*t+a4)^(k/2)[(b1*t+a1)*(h*t^2+g*t+f)]^{-1/2}. 
!*              It is the case of equation W(t)=4*t^3-g_2*t-g_3=0 has one real root. 
!*              Equation (62) of Yang & Wang (2012).  
!*     INPUTS:  Arguments for above integral.
!*     OUTPUTS:  Value of integral. If index_p4(1)=0, then J_0 will be computed, 
!*               else J_0 will be replaced by parameter p, and rff_p=p. 
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.
!*     AUTHOR:     Dexter & Agol (2009)
!*     MODIFIED:   Yang & Wang (2012)
!*     DATE WRITTEN:  4 Mar 2009
!*     REVISIONS: 
!*********************************************************************** 
      Double precision a1,b1,a4,b4,f,g,h,y,x,X1,X4,Y1,Y4,d14,d34,beta1,beta4,a11,c44,integ(4),&
             a142,xi,eta,M2,Lp2,Lm2,I1c,U,U2,Wp2,W2,Q2,P2,rho,I3c,I2c,r24xr34,r12xr13,&
             N2c,K2c,ellcubic,zero,one,two,four,three,half,six,rff,rdd,r14,Ay1,rff_p,By1
      integer   index_p4(4),cases
      PARAMETER ( ONE=1.D0, TWO=2.D0, HALF=0.5d0, THREE=3.D0, FOUR=4.d0, SIX=6.d0, ZERO=0.D0 )
      ellcubic=0.d0
      X1=dsqrt(dabs(a1+b1*x))
      X4=dsqrt(dabs(a4+b4*x))      
      Y1=dsqrt(dabs(a1+b1*y))      
      Y4=dsqrt(dabs(a4+b4*y))      
      r14=a1/b1-a4/b4      
      d14=a1*b4-a4*b1
      !c (2.2) Carlson (1991)
      beta1=g*b1-two*h*a1
      beta4=g*b4-two*h*a4
      !c (2.3) Carlson (1991)
      a11=dsqrt(two*f*b1*b1-two*g*a1*b1+two*h*a1*a1)
      c44=dsqrt(two*f*b4*b4-two*g*a4*b4+two*h*a4*a4)
      a142=two*f*b1*b4-g*(a1*b4+a4*b1)+two*h*a1*a4
      !c (2.4) Carlson (1991)
      xi=dsqrt(f+g*x+h*x*x)
      eta=dsqrt(f+g*y+h*y*y)
      !write(*,*)'ellcubiccomplexs= ',f,g,x,f+g*x+h*x*x,eta
      !c (3.1) Carlson (1991):
      if(x.lt.infinity)then
         M2=((X1+Y1)*dsqrt((xi+eta)**two-h*(x-y)**two)/(x-y))**two
      else 
         M2=b1*(two*dsqrt(h)*eta+g+two*h*y)
      endif            
      !c (3.2) Carlson (1991):
      Lp2=M2-beta1+dsqrt(two*h)*a11
      Lm2=M2-beta1-dsqrt(two*h)*a11

      if(index_p4(1).eq.0) then
      !c (1.2)   Carlson (1991)
        rff=rf(M2,Lm2,Lp2)
        integ(1)=four*rff
        if(cases.eq.1)return
      else
        rff=rff_p/four  
      endif
       !c (3.8)  1991
         I1c=rff_p!four*rff
       !c (3.3) 1991
          if(x.lt.infinity)then      
              U=(X1*eta+Y1*xi)/(x-y)
          else
              U=dsqrt(h)*Y1
              !If(Y1.eq.zero)U=1.D-9
          endif 
          
        !c (3.5) 1991
          rho=dsqrt(two*h)*a11-beta1
          rdd=rd(M2,Lm2,Lp2)
        If(index_p4(3).eq.2)  then
           !  (3.9) Carlson (1991)
           I2c=a11*dsqrt(two/h)/three*(four*rho*rdd-six*rff+three/U)+two*X1*Y1/U
           !  (2.39) Carlson (1989)
           integ(3)=(b4*I2c-d14*I1c)/b1      
           if(cases.eq.3)return            
        endif

         U2=U*U
         Wp2=M2-b1*(a142+a11*c44)/d14
         W2=U2-a11**two*b4/two/d14
       !c (3.4) 1991
          if(x.lt.infinity)then        
                Q2=(X4*Y4/X1/Y1)**two*W2
          else
                Q2=(b4/b1)*(Y4/Y1)**two*W2
         endif       
        P2=Q2+c44**two*b4/two/d14
      !c!! (3.9) 1991
      If(X1*Y1.ne.0.D0) then
        I3c=(two*a11/three/c44)*((-four*b1/d14)*(a142+a11*c44)*rj(M2,Lm2,Lp2,Wp2)&
                        -six*rff+three*rc(U2,W2))+two*rc(P2,Q2)
      Else
        I3c=(two*a11/three/c44)*((-four*b1/d14)*(a142+a11*c44)*rj(M2,Lm2,Lp2,Wp2)&
                        -six*rff+three*rc(U2,W2))
      Endif

        if(index_p4(2).eq.-2) then
        !c (2.49) Carlson (1989)  
          integ(2)=(b4*I3c-b1*I1c)/d14
        if(cases.eq.2)return
        endif

        If(index_p4(4).eq.-4)  then
           ! (2.19) Carlson (1991)
             r24Xr34=half*c44**two/h/b4**two
             r12Xr13=half*a11**two/h/b1**two
           !c (3.11) Carlson (1991)                  
             N2c=two/three*dsqrt(two*h)/a11*(four*rho*rdd-six*rff)!+three/U)!+two/X1/Y1/U
           If(Y1.eq.zero)then
             Ay1=-two*b1*xi/X1+dsqrt(two*h)*a11/U
             If(x.ge.infinity)Ay1=zero!two*d14*eta/Y4**two/U+dsqrt(two*h)*a11/U 
           else
             Ay1=(two*d14*eta/Y4**two+a11**two/X1/U)/Y1+dsqrt(two*h)*a11/U
           endif      
           !c (2.5) & (3.12) Carlson (1991)                  
             K2c=half*a11**two*N2c-two*d14*(xi/X1/X4**two)+Ay1 
           !c (2.62) Carlson (1989)
              integ(4)=-I3c*half/d14*(one/r14+two*b4*beta4/c44**two)+&
              half/d14/(h*b4*r24Xr34)*K2c+(b1/d14)**two*(one-half*r12Xr13/r24Xr34)*I1c  
              !write(*,*)'herer = ',a11,N2c,d14,xi,X1,X4,Ay1 
        endif 
      return
      end  subroutine ellcubiccomplexs
!********************************************************************************************
      subroutine  elldoublecomplexs(index_p5,f1,g1,h1,f2,g2,h2,a5,b5,y,x,rff_p,integ,cases)
!*******************************************************************************
!*     PURPOSE: Computes J_k(h)=\int_y^x dt (f_1+g_1t+h_1t^2)^{p_1/2} 
!*                       (f_2+g_2t+h_2t^2)^{p_2/2} (a_5+b_5t)^{p_5/2}. 
!*     INPUTS:  Arguments for above integral.
!*     OUTPUTS:  Value of integral. If index_p5(1)=0, then J_0 will be computed, 
!*               else J_0 will be replaced by parameter p, and rff_p=p. 
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Dexter & Agol (2009)
!*     MODIFIED:   Yang & Wang (2012)
!*     DATE WRITTEN:  4 Mar 2009
!*     REVISIONS: [-1,-1,-1,-1,0]=integ(1),[-1,-1,-1,-1,-2]=integ(2),
!*                [-1,-1,-1,-1,-4]=integ(4),[-1,-1,-1,-1,2]=integ(3),
!*                [-1,-1,-1,-1,4]=integ(5)
!***********************************************************************
      implicit NONE 
      Double precision one,half,two,three,four,six,f1,g1,h1,f2,g2,h2,a5,b5,y,x,xi1,xi2,eta1,eta2,integ(5),&
                       theta1,theta2,zeta1,zeta2,M,M2,delta122,delta112,delta222,delta,deltap,lp2,lm2,&
                       deltam,rff,ellquartic,U,U2,alpha15,beta15,alpha25,beta25,lambda,omega2,psi,xi5,eta5,&
                       gamma1,gamma2,Am111m1,A1111m4,XX,S,mu,T,V2,b2,a2,H,A1111m2,xi1p,B,G,Sigma,Lambda0,&
                       Omega02,psi0,X0,mu0,b02,a02,H0,S2,T2,eta1p,T02,V02,psi2,T0,A1111,rff_p
      integer  index_p5(5),cases
      one=1.D0
      half=0.5D0
      two=2.D0
      three=3.D0
      four=4.D0
      six=6.D0
      !c (2.1) Carlson (1992)
      If(x.lt.infinity)then
          xi1=dsqrt(f1+g1*x+h1*x**two)
          xi2=dsqrt(f2+g2*x+h2*x**two)
      else
          xi1=x*dsqrt(h1)
          xi2=x*dsqrt(h2)      
      endif            
      eta1=dsqrt(f1+g1*y+h1*y*y)
      eta2=dsqrt(f2+g2*y+h2*y*y)
      !c (2.4) Carlson (1992)
      If(x.lt.infinity)then
          theta1=two*f1+g1*(x+y)+two*h1*x*y
          theta2=two*f2+g2*(x+y)+two*h2*x*y
      else
          theta1=(g1+two*h1*y)*x
          theta2=(g2+two*h2*y)*x
      endif      
      !c (2.5) Carlson (1992)
      zeta1=dsqrt(two*xi1*eta1+theta1)
      zeta2=dsqrt(two*xi2*eta2+theta2)
      !c (2.6) Carlson (1992)
      If(x.lt.infinity)then
          M=zeta1*zeta2/(x-y)
          M2=M*M
      else
         M2=(two*dsqrt(h1)*eta1+g1+two*h1*y)*(two*dsqrt(h2)*eta2+g2+two*h2*y)
      endif

      !c (2.7) Carlson (1992)
      delta122=two*f1*h2+two*f2*h1-g1*g2
      delta112=four*f1*h1-g1*g1
      delta222=four*f2*h2-g2*g2
      Delta=dsqrt(delta122*delta122-delta112*delta222)
      !c (2.8) Carlson (1992)
      Deltap=delta122+Delta
      Deltam=delta122-Delta
      Lp2=M2+Deltap
      Lm2=M2+Deltam
 
      if(index_p5(1).eq.0) then
        rff=rf(M2,Lm2,Lp2)
        !c (2.36) Carlson (1992)
        integ(1)=four*rff
        if(cases.eq.1)return
      else
        rff=rff_p/four 
      endif
      !c (2.6) Carlson (1992)
          If(x.lt.infinity)then
             U=(xi1*eta2+eta1*xi2)/(x-y)
             U2=U*U
         else
             U=dsqrt(h1)*eta2+dsqrt(h2)*eta1
             U2=U*U      
         endif
        
      !c (2.11) Carlson (1992)
        alpha15=two*f1*b5-g1*a5 
        alpha25=two*f2*b5-g2*a5
        beta15=g1*b5-two*h1*a5 
        beta25=g2*b5-two*h2*a5
      !c (2.12) Carlson (1992)
        gamma1=half*(alpha15*b5-beta15*a5)
        gamma2=half*(alpha25*b5-beta25*a5)
      !c (2.13) Carlson (1992)
        Lambda=delta112*gamma2/gamma1
        Omega2=M2+Lambda
        psi=half*(alpha15*beta25-alpha25*beta15)
        psi2=psi*psi
      !c (2.15) Carlson (1992)
        xi5=a5+b5*x
        eta5=a5+b5*y
      !c (2.16) Carlson (1992)
          If(x.lt.infinity)then
              Am111m1=one/xi1*xi2-one/eta1*eta2
              A1111m4=xi1*xi2/xi5**two-eta1*eta2/eta5**two
              A1111m2=xi1*xi2/xi5-eta1*eta2/eta5      
              XX = (xi5*eta5*theta1*half*Am111m1-xi1*xi2*eta5**two+&
                                  eta1*eta2*xi5**two)/(x-y)**two
              mu=gamma1*xi5*eta5/xi1/eta1
          else
              Am111m1=dsqrt(h2/h1)-one/eta1*eta2
              A1111m4=dsqrt(h1*h2)/b5**two-eta1*eta2/eta5**two      
              A1111m2=dsqrt(h1*h2)*x/b5-eta1*eta2/eta5
              XX=b5*eta5*h1*y*Am111m1-eta5**two*dsqrt(h1*h2)+eta1*eta2*b5**two
              mu=gamma1*b5*eta5/dsqrt(h1)/eta1
         endif

      !c (2.17) Carlson (1992)
    
      !c (2.18) Carlson (1992)
        S=half*(M2+delta122)-U2
        S2=S*S
      !c (2.19) Carlson (1992)
        T=mu*S+two*gamma1*gamma2
        T2=T*T
        V2=mu**two*(S2+Lambda*U2)
      !c (2.20) Carlson (1992)
        b2=Omega2**two*(S2/U2+Lambda)
        a2=b2+Lambda**two*psi2/gamma1/gamma2
      !c (2.22) Carlson (1992)
        H=delta112*psi*(rj(M2,Lm2,Lp2,Omega2)/three+&
             half*rc(a2,b2))/gamma1**two-XX*rc(T2,V2)
      If(index_p5(3).eq.2 .or. index_p5(5).eq.4)then
            !(2.23)--(2.29) Carlson (1992)
            psi=g1*h2-g2*h1
            Lambda=delta112*h2/h1
            Omega2=M2+Lambda
            Am111m1=one/xi1*xi2-one/eta1*eta2
            A1111=xi1*xi2-eta1*eta2
            XX=(theta1*half*Am111m1-A1111)/(x-y)**two
            b2=Omega2**two*(S2/U2+Lambda)
            a2=b2+Lambda**two*psi**two/h1/h2
            mu=h1/xi1/eta1
            T=mu*S+two*h1*h2
            T2=T*T
            V2=mu**two*(S2+Lambda*U2)
            H0=delta112*psi*(rj(M2,Lm2,Lp2,Omega2)/three+&
                          half*rc(a2,b2))/h1**two-XX*rc(T2,V2)
          If(index_p5(3).eq.2)then
            !(2.42) Carlson (1992)
            integ(3)=two*b5*H0-two*beta15*rff/h1
            if(cases.eq.3)return
          endif
          If(index_p5(5).eq.4)then      
              !c (2.2) Carlson (1992)
                If(x.lt.infinity)then
                 xi1p=half*(g1+two*h1*x)/xi1
                else
                 xi1p=dsqrt(h1)
                endif
                eta1p=half*(g1+two*h1*y)/eta1
              !c (2.3) Carlson (1992)
                B=xi1p*xi2-eta1p*eta2
              !c (2.9) Carlson (1992)
                If(x.lt.infinity)then
                   G=two/three*Delta*Deltap*rd(M2,Lm2,Lp2)+half*Delta/U+&
                       (delta122*theta1-delta112*theta2)/four/xi1/eta1/U  
                else
                   G=two/three*Delta*Deltap*rd(M2,Lm2,Lp2)+half*Delta/U+&
                  (delta122*(g1+two*h1*y)-delta112*(g2+two*h2*y))/four/dsqrt(h1)/eta1/U  
                endif
              !c (2.10) Carlson (1992)  
                  Sigma=G-Deltap*rff+B
              !(2.44) Carlson (1992)
                  integ(5)=-b5*(beta15/h1+beta25/h2)*H0+b5**two*&
                            Sigma/h1/h2+beta15**two*rff/h1**two
               if(cases.eq.5)return
          endif             
      endif
            if (index_p5(2).eq.-2) then
              !c (2.39) Carlson (1992)
                integ(2)=-two*(b5*H+beta15*rff/gamma1)
            if(cases.eq.2)return
            endif
          If(index_p5(4).eq.-4)then
              !c (2.2) Carlson (1992)
                If(x.lt.infinity)then
                 xi1p=half*(g1+two*h1*x)/xi1
                else
                 xi1p=dsqrt(h1)
                endif
                eta1p=half*(g1+two*h1*y)/eta1
              !c (2.3) Carlson (1992)
                B=xi1p*xi2-eta1p*eta2
              !c (2.9) Carlson (1992)
                If(x.lt.infinity)then
                G=two/three*Delta*Deltap*rd(M2,Lm2,Lp2)+half*Delta/U+&
                         (delta122*theta1-delta112*theta2)/four/xi1/eta1/U  
                else
                 G=two/three*Delta*Deltap*rd(M2,Lm2,Lp2)+half*Delta/U+&
                  (delta122*(g1+two*h1*y)-delta112*(g2+two*h2*y))/four/dsqrt(h1)/eta1/U  
                endif
              !c (2.10) Carlson (1992)  
                Sigma=G-Deltap*rff+B
              !c (2.41) Carlson (1992)
                integ(4)=b5*(beta15/gamma1+beta25/gamma2)*H+beta15**two*rff/gamma1**two+&
                                   b5**two*(Sigma-b5*A1111m2)/gamma1/gamma2
          endif
      integ=ellquartic
      return
      end  subroutine  elldoublecomplexs
!*************************************************************************** 
      end module ellfunction


