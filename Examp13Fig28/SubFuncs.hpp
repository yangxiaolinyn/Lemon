#include "RandUtils.hpp"
 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
    double cosm_phi_sampling( const int& m, const double& Phi1, const double& rnum )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    { 
        double temp_phi1, m2, xi; 
 
        //cout << " sss3 = " << m << endl; 
        if (m == 0) 
            return Phi1;
        else
        {   
            xi = rnum ;  
            if (  xi * two <= one + cos( m*Phi1 ) ) 
            {
                temp_phi1 = Phi1; 
            }
            else
            {
                m2 = floor( Phi1 / ( pi/m ) );
                temp_phi1 = (two * m2 + one) * pi / m - Phi1; 
            }  
            return temp_phi1;
        } 
    } 
