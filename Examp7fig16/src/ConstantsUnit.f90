
!************************************************************
      module constants
!********************************************************************
!*    This module defines many constants often uesd in our code.
!*    One can use these constants through a command "use constants" in their
!*    own subroutines or functions. 
!********************************************************************
        Implicit None
        integer :: rand_inst = 0
        integer, parameter :: krand = KIND(1.d0)
        integer, parameter :: mcp= KIND(1.d0)
        integer :: Rand_Feedback = 1
        double precision, parameter :: zero = 0.D0, &
              one = 1.D0, &
              two = 2.D0, &
              three = 3.D0, &
              four = 4.D0, &
              five = 5.D0, &
              six = 6.D0, &
              seven = 7.D0, &
              eight = 8.D0, &
              nine = 9.D0, &
              ten = 10.D0, &
              sixteen = 16.D0, &
              sqrt2 = dsqrt(2.D0), &
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              pi = 3.1415926535897932384828D0,&
              twopi = 6.28318530717958647692D0, &
              halfpi = 1.57079632679489661923132169163D0, & 
              !sqrtpi = 1.772453850905516027298167483341145182797549456122387128213...
              sqrtpi = 1.77245385090551602729D0, &
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              half = 0.5D0, &
              half2 = 0.25D0, &
              sigma_T = 6.65D-25, &
              barn = 1.D-24, &
              mec2 = 0.511D0, &
              Boltzman_Costant_k = 1.380662D-16, &
              electron_mass = 9.109534D-28, &
              electron_Charge = 4.803242D-10, &
              ln_e = 2.718281828459045235360D0, &
              infinity=1.D100, & 
              Boltzman_k_ergK = 1.380662D-16, &
              Boltzman_k_evK = 1.380662D-16 / 1.6021892D-12, &
              dtors=3.1415926535897932384828D0/180.D0, &
              dtor=3.1415926535897932384828D0/180.D0, &
              mh=1.6726231D-24, &
              hbar = 1.0545887D-27, & ! h / 2 pi
              planck_h=6.626178D-27, & !erg * s
              h_ev = 4.13566743D-15, & !eV·s
              pho_v=2.99792458D10, & 
              erg_of_one_ev = 1.6021892D-12, &
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Gnewton = 6.6720D-8, &
              Msun = 1.9891D33, &
              Cv = 2.99792458D10, &
              rg_SUN = 6.6720D-8 * 1.9891D33 / 2.99792458D10**2
!********************************************************************************************
      end module constants 
!********************************************************************************************
 
