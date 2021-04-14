#include "SubFuncs.hpp"
#include <mpi.h>
#include <fstream>


int main( )
{ 
    //~~~~~~~~~~~~~ MPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int myid, np;
    int namelen;
    char processor_name[ MPI_MAX_PROCESSOR_NAME ];
    int argc;
    char** argv;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &np );
    MPI_Comm_rank (MPI_COMM_WORLD, &myid );
    MPI_Get_processor_name(processor_name, &namelen );

    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << " My id is :" << myid << " np = " << np << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

    long long mydutyphot, Total_Smple_Num;
    Total_Smple_Num = 5.0 * pow(10, 7);
    if( myid == np-1 ) 
        mydutyphot = Total_Smple_Num / np + ( Total_Smple_Num % np );
    else
        mydutyphot = Total_Smple_Num / np;

    //~~~~~~~~~~~~~ MPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    RNumGenerate rnums;
    rnums.initRandom();

    //~~~ Set initial parameters and conditions ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    double Normal = {zero};
    long const m1 = 20;
    long const m2 = 10;
    long const nn = 500;
    double am[m1 + 1] = {zero};
    double bm[m2 + 1] = {zero};
    double prob[m1 + 1] = {zero};
    double Pr[m1 + 1] = {zero}; 
    double xi[nn+2] = {zero};
 
    double delta_x = twopi / nn;
    double F_sin_mphi( long n, double& x_var, double bm[] );
    double F_cos_mphi( long n, double& x_var, double bm[] );

    if( myid == np-1 )
    { 
        for( int i = 1; i <= m1; i++) 
        { 
            am[i] = rnums.RANMAR(); 
            Normal += am[i]; 
        };

        for( int i = 1; i <= m2; i++) 
        { 
            bm[i] = rnums.RANMAR();
            Normal += bm[i]; 
        };

        Normal *= 1.10;
        am[0] = one; // / twopi;
        for( int i = 1; i <= m1; i++) 
        { 
            am[i] /= Normal; 
        };

        for( int i = 1; i <= m2; i++) 
        { 
            bm[i] /= Normal; 
        };

        prob[0] = one;
        for( int i = 1; i <= m1; i++) 
        { 
            prob[0] -= am[i];
            prob[i] = am[i];
        }

        double pro = zero;
        pro = prob[0]; 
        Pr[0] = pro;
        for( int i = 1; i <= m1; i++) 
        {
            pro += prob[i];
            Pr[i] = pro; 
        } 
        //MPI_Barrier( MPI_COMM_WORLD ); 
        //cout << " ass = " << am[m1] << endl;  
    } 
   
    MPI_Bcast( am, m1 + 1, MPI_DOUBLE_PRECISION, np-1, MPI_COMM_WORLD );
    MPI_Bcast( bm, m2 + 1, MPI_DOUBLE_PRECISION, np-1, MPI_COMM_WORLD ); 
    MPI_Bcast( prob, m1 + 1, MPI_DOUBLE_PRECISION, np-1, MPI_COMM_WORLD ); 
    MPI_Bcast( Pr, m1 + 1, MPI_DOUBLE_PRECISION, np-1, MPI_COMM_WORLD );
 

    double Phi1; 
    long long i_count = 0;
    do
    { 
        Phi1 = rnums.RANMAR() * twopi;
          
        double p_1 = zero;
        double temp_phi;
        p_1 = rnums.RANMAR();
        for( int i = 0; i <= m1; i++)
        { 
            if( p_1 <= Pr[i] ) 
            {  
                temp_phi = cosm_phi_sampling( i, Phi1, rnums.RANMAR() ); 
                //cout << " P2 = "<< Pr[i] << endl; 
                break;
            } 
        }
 
        double F_phi = zero, H_phi = zero, x; 
        for( int k = 0; k <= m1; k++ )
            F_phi += am[k] * cos(k * temp_phi);
        for( int k = 1; k <= m2; k++ )
            H_phi += bm[k] * sin(k * temp_phi); 
 

        if( (rnums.RANMAR() - one) <= (H_phi / F_phi) )
            x = temp_phi;
        else
            x = twopi - temp_phi; 
 
        int i;
        i = floor(x / delta_x) + 1; 

        //cout << " ss iii = " << i << "   dex = " << delta_x << endl;
        //cout << " x = " << x << endl;
        if ( i > nn )
            continue; 

        xi[i] += 1.0; 

        i_count += 1;
        //cout << " My id is = " << myid << " i_count = " << i_count <<  endl;
        //cout << " i is = " << i << " i_count = " << xi[i] <<  endl;
        if( i_count >= mydutyphot ) 
            break;

        if( myid == np-1 && i_count % 500000 == 0 ) 
            cout << " J is = " << i_count << " i_count = " << mydutyphot << endl;

    }while(true);
 
    // collect data from all other MPI processors.
     if ( myid != np-1 )
     {
         if ( np-2 >= 0 )  
         {
             long send_num, send_tag;
             send_num  = nn + 2;
             send_tag = 1;
             MPI_Send( xi, send_num, MPI_DOUBLE_PRECISION, np-1, \
                        send_tag, MPI_COMM_WORLD );
             cout << "MPI Processor: " << myid << " has send xi array to Processor " << np-1 << endl;
          } 
     }
     else
     {
         if (np-2 >= 0)
         { 
             for( int RECV_SOURCE = 0; RECV_SOURCE < np - 1; RECV_SOURCE++ ) 
             {
                 long recv_num, recv_tag;
                 double xi_RECV[nn + 2]; 
                 recv_num = nn + 2;
                 recv_tag = 1;
                 MPI_Recv( xi_RECV, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, \
                           recv_tag, MPI_COMM_WORLD, &status );
                 for( int i = 0; i < nn + 2; i++ )
                      xi[i] += xi_RECV[i];
             }
         }

         //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ofstream myFile;
         myFile.open("./fig/dataFF_xx.dat", ios_base::out); 
         if (myFile.is_open())
         { 
             cout << "File open successful" << endl;  
             cout << "Finished writing to file, will close now" << endl; 
             double Px, x_var;
             for( int i = 0; i < nn + 2; i++ )
             {
                 x_var = i * delta_x;
                 Px = F_cos_mphi( m1, x_var, am ) + F_sin_mphi( m2, x_var, bm );
                 myFile << xi[i] << "      " << x_var << "      " << Px << endl;
             }
             myFile.close();
         }
         //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         /*ofstream myFile2;
         myFile2.open("./dataXX.txt", ios_base::out); 
         if (myFile2.is_open())
         {  
             for( int i = 0; i < nn + 2; i++ )
             {
                 myFile2 << (i * delta_x) << endl;
             }
             myFile2.close();
         }*/
         //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     }

     
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MPI_Finalize(); 
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return 0;
}

 
double F_cos_mphi( long m, double& phi, double am[] )
{
    double temp_v; 

    temp_v = zero;
    for(long i = 0; i <= m; i++)
        temp_v += am[i] * cos(i * phi);
    return temp_v;
} 
 
double F_sin_mphi( long m, double& phi, double am[] )
{
    double temp_v; 

    temp_v = zero;
    for(long i = 0; i <= m; i++)
        temp_v += am[i] * sin(i * phi);
    return temp_v;
} 








