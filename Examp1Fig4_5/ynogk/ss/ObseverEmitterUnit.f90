!*******************************************************************************
      module obsemitter
!*******************************************************************************
!*     PURPOSE: This module aims on determining geodesic connecting observer
!*              and emitter, i.e. to solve equations (121)-(123) in  
!*              Yang & Wang (2012) by Newton Raphson method.     
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  15 Jan 2012. 
!*******************************************************************************  
      use constants
      use rootsfinding
      use ellfunction
      use pemfinding
      use BLcoordinate         
      implicit none
!*****************************************************************************************************
      contains
!*****************************************************************************************************
      SUBROUTINE alphabetap(a0,B0,alen,Blen,sinobs,muobs,a_spin,robs,scal,&
                              obs_V,rmuphy_em,abp,func1,func2,func3)
      IMPLICIT NONE
!*******************************************************************************
!*    PURPOSE:   This subroutine aim to determine motion constants of a geodesic which connecting the 
!*               emmiter and the observer. The Boyer-Lindquist coordinates of them (the emmiter and 
!*               observer) have given.
!************************************************************************
!*     INPUTS:    (a0,B0)--------is the coordinates of one corner of the rectangle on the 
!*                               screen plane of the observer.
!*                alen,Blen------size of the rectangle. 
!*                sinobs,muobs---sinobs=sin(theta_obs) and muobs=cos(theta_obs), theta_obs is the 
!*                               theta coordinate of the observer.
!*                a_spin---------spin of black hole.                       
!*                robs-----------The radiual coordinate of the observer.
!*                scal-----------a parameter to control the size of the image.   
!*                               The value of which usually was set to be 1.D0.
!*                obs_V(1:3)-----array of the velocity of the observer respect to the ZAMO( or LNRF).    
!*                rmuphy_em------array of Boyer-Lindquist coordinates of the emmiter. 
!*                func1,func2,func3-----Name of functions: r(p), \mu(p), \phi(p).
!*     OUTPUTS:   abp(1:3)-------array of roots of equations (121)-(123) in Yang & Wang (2012),  
!*                               \alpha_em, \beta_em, p_em.
!*     ROUTINES CALLED: lambdaq, mnewt2.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012       
!*******************************************************************************  
      INTEGER n,ntrial,i,j,k,i1,j1,ntriend
      INTEGER tr1,tr2,reals,del,ax,Bx
      PARAMETER(n=3,ntrial=1000)
      DOUBLE PRECISION a0,B0,a_spin,cobs,robs,scal,x(n),xem(n),tolx,tolf,sinobs,muobs,a1,B1,&
                   Delta,sigma,bigA,somiga,expmu,zero,one,two,three,four,rmss,&
                   fvec(n),rem,muem,phyem,p,mu_tp,mu_tp2,r_tp,lambda,q,pin,&         
                   g,re,ut,bomiga,deltaa0,deltaB0,abp(3),rmuphy_em(3),pem1,pem2,&
                   mu1,mu2,phy1,phy2,f1234(4),obs_V(3),alphac,Betac,alpha1,Beta1,&
                   alen,Blen,adel,Bdel        
      PARAMETER(two=2.D0,three=3.D0,four=4.D0,one=1.D0,zero=0.D0)
      DOUBLE PRECISION,EXTERNAL :: func1,func2,func3
      LOGICAL :: mobseqmtp,err

      p=one
      tr1=0
      tr2=0
      rmss=rms(a_spin) 
      tolx=1.D-4
      tolf=1.D-4
      xem(1)=rmuphy_em(1)
      xem(2)=rmuphy_em(2)
      rmuphy_em(3)=DMOD(rmuphy_em(3),twopi)
      If(rmuphy_em(3).lt.zero)rmuphy_em(3)=rmuphy_em(3)+twopi
      xem(3)=rmuphy_em(3)
      adel=0.1D0
      Bdel=0.1D0         
        B1=-sqrt(xem(1)**two+a_spin**two)*cos(xem(3))*zero
        a1=-sqrt(xem(1)**two+a_spin**two)*sin(xem(3))*zero
        Do i=0,floor(Blen/Bdel)
            Beta1=dabs(B0)-i*Bdel
            Do j=0,floor(alen/adel)
            alpha1=-dabs(a0)+j*adel
                x(1)=alpha1
                x(2)=Beta1
                tr1=0
                tr2=0
                call lambdaq(x(1),x(2),robs,sinobs,muobs,a_spin,scal,obs_V,f1234,lambda,q)
                x(3)=r2p(f1234(1),xem(1),lambda,q,a_spin,robs,scal,tr1,tr2)
                If(x(3).lt.zero)cycle        
                call mnewt2(ntrial,x,n,tolx,tolf,sinobs,muobs,a_spin,robs,scal,&
                                                 obs_V,xem,err,func1,func2,func3) 
                !write(*,*)'i,j',i,j,alpha1,Beta1,lambda,q,x,err
                If(.not.err)then
                !write(*,*)'ggg',x,err
                        abp(1)=x(1)
                        abp(2)=x(2)
                        abp(3)=x(3)
                        goto 100
                endif           
            Enddo   
        Enddo   
100     continue
        return
      END subroutine alphabetap
!********************************************************************************* 
      SUBROUTINE mnewt2(ntrial,x,n,tolx,tolf,sinobs,muobs,a_spin,robs,&
                                 scal,obs_V,xend,err,func1,func2,func3) 
!************************************************************************************** 
!*    PURPOSE:    Given an initial guess x for a root in n dimensions,take ntrial 
!*                Newton-Raphson steps to improve the root.Stop if the root converges 
!*                in either summed absolute variable imcrements tolx
!*                or summed absolute function values tolf.
!*         
!*     INPUTS:    
!*     ROUTINES CALLED:  lubksb,ludcmp,usrfun.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Press et al. (2007)  
!*     DATE WRITTEN:   
!*******************************************************************************  
        IMPLICIT NONE
        INTEGER n,ntrial,NP
        INTEGER tr1,tr2
        DOUBLE PRECISION tolx,tolf,x(n),zero,xend(n)
        PARAMETER(NP=3,zero=0.D0)     !Up to NP variables.
        INTEGER i,k,indx(NP)
        DOUBLE PRECISION d,errf,errx,fjac(NP,NP),fvec(NP),p(NP),lambda,q,&
                         a,B,p_null,sinobs,muobs,a_spin,robs,scal,f1234(4),obs_V(3)
        DOUBLE PRECISION, EXTERNAL :: func1,func2,func3
        LOGICAL err,err1

        do k=1,ntrial 
            call usrfun2(x,n,NP,fvec,fjac,sinobs,muobs,a_spin,robs,scal,&
                                          obs_V,xend,err,func1,func2,func3) 
            If(err)then
                return
            endif
!User subroutine supplies function valuse at x in fvec and
!Jacobian matrix in fjac. Check function convergence 
            errf=zero                           
            do i=1,n                            
                errf=errf+abs(fvec(i))
            enddo 
            if(errf.le.tolf)return
            do i=1,n
                p(i)=-fvec(i)
            enddo        
!Solve linear equations using LU decomposition.
            call ludcmp(fjac,n,NP,indx,d,err)  
                If(err)then
!write(unit=6,fmt=*)'ludcmp(): err=',err,x(3)
                    return        
                endif
            call lubksb(fjac,n,NP,indx,p) 
!Check root convergence.
            errx=zero                                 
            do i=1,n
                errx=errx+abs(p(i))
                x(i)=x(i)+p(i)
            enddo
            If(x(3).lt.zero)then
                err=.true.
                return
            endif        
!write(*,*)'ntrials=',x
            if(errx.le.tolx)return
        enddo
        return
      END SUBROUTINE mnewt2
!********************************************************************* 
      SUBROUTINE usrfun2(x,n,NP,fvec,fjac,sinobs,muobs,a_spin,robs,&
                                   scal,obs_V,xend,err,func1,func2,func3)
!*********************************************************************
!*    PURPOSE:   This subroutine aim to compute array: fvec and fjac.
!************************************************************************
!*     INPUTS:     
!*     OUTPUTS:  fvec and fjac 
!*     ROUTINES CALLED: lambdaq, mnewt2.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012       
!*******************************************************************************  
      USE BLcoordinate 
      IMPLICIT NONE                
      INTEGER  n,NP,i,j
      DOUBLE PRECISION x(n),fjac(NP,NP),fvec(NP),sinobs,muobs,a_spin,robs,scal,&
                   ra,mua,phya,deltax,ra1,mua1,phya1,f(NP),h,EPS,&
                   temp,xend(n),rhorizon,lambda,q,f1234(4),obs_V(3),&
                   temp_aB(2),a,B
      PARAMETER(deltax=1.D-6,EPS=1.D-5)
      DOUBLE PRECISION,EXTERNAL :: func1,func2,func3
      LOGICAL err
        
      err=.false.
      rhorizon=one+sqrt(one-a_spin**two)
      call lambdaq(x(1),x(2),robs,sinobs,muobs,a_spin,scal,obs_V,f1234,lambda,q)
      fvec(1)=func1(x(3),f1234(1),lambda,q,a_spin,robs,scal)-xend(1)
      fvec(2)=func2(x(3),f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)-xend(2)
      fvec(3)=func3(x(3),f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)-xend(3)

      If(fvec(1)+xend(1).ge.robs.or.fvec(1)+xend(1).le.rhorizon)then
          err=.true.
          return
      Endif
                
      Do j=1,n
          temp=x(j)
          h=EPS*ABS(temp)
          IF(h.eq.zero)h=EPS        
          x(j)=temp+h
          call lambdaq(x(1),x(2),robs,sinobs,muobs,a_spin,scal,obs_V,f1234,lambda,q)
          h=x(j)-temp
          f(1)=func1(x(3),f1234(1),lambda,q,a_spin,robs,scal)-xend(1)
          f(2)=func2(x(3),f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)-xend(2)
          f(3)=func3(x(3),f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)-xend(3)
          x(j)=temp
          Do i=1,n
              fjac(i,j)=(f(i)-fvec(i))/h
          Enddo
      Enddo        

      END SUBROUTINE usrfun2
!********************************************************************* 
      SUBROUTINE ludcmp(a,n,np,indx,d,err) 
!***********************************************************************************************
!*    PURPOSE:     PGiven a matrix a(1:n,1:n), with physical dimension np by np, this routine 
!*                 replaces it by the LU decomposition of a rowwise permutation of itself. a 
!*                 and n are input. a is output, arranged as in equation (2.3.14) above; 
!*                 indx(1:n) is an output vector that records the row permutation effected 
!*                 by the partial pivoting; d is output as +/-1 depending on whether
!*                 the number of row interchanges was even or odd, respectively. This 
!*                 routine is used in combination with lubksb to solve linear equations 
!*                 or invert a matrix. 
!*
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Press et al. (2007)  
!*     DATE WRITTEN:  
!***********************************************************************************************
      IMPLICIT NONE
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY 
      PARAMETER (NMAX=500,TINY=1.0D-20)     !Largest expected n, and a small number.
      INTEGER i,imax,j,k
      LOGICAL  err
      DOUBLE PRECISION aamax,dum,sums,vv(NMAX)   !vv stores the implicit scaling of each row.
      d=one                            !No row interchanges yet.
      err=.false. 
      do i=1,n                         !Loop over rows to get the implicit scaling information.
          aamax=zero
          do j=1,n
              if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
              !write(*,*)'aa==',aamax
          enddo
          if (aamax.eq.zero)then
              !write(unit=6,fmt=*)'pause singular matrix in ludcmp'
              !write(*,*)'ss==',a(i,1),a(i,2),a(i,3)
              err=.true.
              return
              !stop
          endif         !pause 'singular matrix in ludcmp'   !No nonzero largest element.
          vv(i)=one/aamax                                   !Save the scaling.
      enddo
      do  j=1,n                       !This is the loop over columns of Crout's method.
          do  i=1,j-1                 !This is equation (2.3.12) except for i = j.
              sums=a(i,j)
              do  k=1,i-1
                  sums=sums-a(i,k)*a(k,j)
              enddo 
              a(i,j)=sums
          enddo 
          aamax=zero      !Initialize for the search for largest pivot element.
          do  i=j,n            !This is i = j of equation (2.3.12) and i = j+1: ::N
              sums=a(i,j)                 !of equation (2.3.13).
              do  k=1,j-1
                  sums=sums-a(i,k)*a(k,j)
              enddo 
              a(i,j)=sums
              dum=vv(i)*abs(sums)                 !Figure of merit for the pivot.
              if (dum.ge.aamax) then         !Is it better than the best so far?
                  imax=i
                  aamax=dum
              endif
          enddo 
          if (j.ne.imax)then                 !Do we need to interchange rows?
              do  k=1,n                         !Yes, do so...
                  dum=a(imax,k)
                  a(imax,k)=a(j,k)
                  a(j,k)=dum
              enddo 
              d=-d                         !...and change the parity of d.
              vv(imax)=vv(j)                 !Also interchange the scale factor.
          endif
          indx(j)=imax
          if(a(j,j).eq.zero)a(j,j)=TINY
!*  If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!*  For some applications on singular matrices, it is desirable to substitute TINY
!*  for zero.
          if(j.ne.n)then
!Now, finally, divide by the pivot element.                                
              dum=one/a(j,j)
              do  i=j+1,n
                  a(i,j)=a(i,j)*dum
              enddo 
          endif
!Go back for the next column in the reduction.
      enddo   
      return
      END  SUBROUTINE ludcmp
!**************************************************************************************
      SUBROUTINE lubksb(a,n,np,indx,b) 
!**************************************************************************************
!*    PURPOSE:   Solves the set of n linear equations A \dot X = B. Here a is input, 
!*               not as the matrix A but rather as its LU decomposition, determined 
!*               by the routine ludcmp. indx is input as the permutation vector returned 
!*               by ludcmp. b(1:n) is input as the right-hand side vector B, and returns 
!*               with the solution vector X. a, n, np, and indx are not modified by this 
!*               routine and can be left in place for successive calls with diffierent 
!*               right-hand sides b. This routine takes into account the possibility that 
!*               b will begin with many zero elements, so it is efficient for use in 
!*               matrix inversion.
!*
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Press et al. (2007)  
!*     DATE WRITTEN:  
!**************************************************************************************
      IMPLICIT NONE
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n),zero,one,sums
      INTEGER  i,ii,j,ll
      ii=0      !When ii is set to a positive value, it will become the index
                !of the first nonvanishing element of b. We now do
                !the forward substitution, equation (2.3.6). The only new
                !wrinkle is to unscramble the permutation as we go.
      do  i=1,n
          ll=indx(i)
          sums=b(ll)
          b(ll)=b(i)
          if (ii.ne.0)then
              do  j=ii,i-1
                  sums=sums-a(i,j)*b(j)
              enddo 
          else 
!A nonzero element was encountered, so from now on we will
!have to do the sums in the loop above.
              if (sums.ne.zero) then
                  ii=i                
              endif                
          endif
          b(i)=sums
      enddo 
!Now we do the backsubstitution, equation (2.3.7).
      do  i=n,1,-1                        
          sums=b(i)
          do  j=i+1,n
          sums=sums-a(i,j)*b(j)
          enddo 
!Store a component of the solution vector X.
          b(i)=sums/a(i,i)         
      enddo 
!All done!
      return                                 
      END  SUBROUTINE lubksb
!************************************************************ 
      ENd module obsemitter



