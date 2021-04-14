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
    Total_Smple_Num = 2.5 * pow(10, 8);
    if( myid == np-1 ) 
        mydutyphot = Total_Smple_Num / np + ( Total_Smple_Num % np );
    else
        mydutyphot = Total_Smple_Num / np;

    //~~~~~~~~~~~~~ MPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    RNumGenerate rnums;
    rnums.initRandom();

    //~~~ Set initial parameters and conditions ~~~~~~~~~~~~~~~~ 
    long const nn = 500; 

    double xi[nn+2] = {zero};
 
    double x_max = 3.0;
    double Sigma_max = x_max*x_max;
    double N_normal =  one - exp(-Sigma_max * x_max); 

    double delta_x = x_max / nn; 
    double func_x( double& x );
  
    double x, x1; 
    long long i_counts = 0;
    long long num_count = 0;
    bool continu;
    do
    { 
        x = zero;
        i_counts = 0;
        continu = false;
       
        do
        {
            x1 = - log(one - rnums.RANMAR() * N_normal) / Sigma_max; 
            i_counts += 1;
            x = x + x1;
            if (x > x_max)
            {
                continu = true;
                break;
            } else if ( rnums.RANMAR() <= ( func_x( x ) / Sigma_max ) )
            {
                if( i_counts == 1 )
                    break;
                else if( rnums.RANMAR() <= N_normal )
                    break;
                else
                {
                    continu = true;
                    break;
                } 
            } else 
            { 
                if( i_counts == 1 )
                {}
                else
                {
                    if( rnums.RANMAR() > N_normal )
                    {
                        continu = true;
                        break;
                    }
                }
            }
        }while(true);
        
        if(continu)
            continue; 

        int i;
        i = floor(x / delta_x) + 1; 
 
        if ( i > nn )
            continue; 

        xi[i] += 1.0; 

        num_count += 1; 
        if( num_count >= mydutyphot ) 
            break;

        if( myid == np-1 && num_count % 500000 == 0 ) 
            cout << " J is = " << num_count << " i_count = " << mydutyphot << endl;

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
         myFile.open("./fig/dataFF_xx.txt", ios_base::out); 
         if (myFile.is_open())
         { 
             cout << "File open successful" << endl;  
             cout << "Finished writing to file, will close now" << endl; 
             double Px, x_var;
             for( int i = 0; i < nn + 2; i++ )
             {
                 x_var = i * delta_x; 
                 Px = pow(x_var, 2) * exp( - pow(x_var, 3) / 3.0 );
                 myFile << xi[i] << "      " << x_var << "      " << Px << endl;
             }
             myFile.close();
         } 
     }

     
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MPI_Finalize(); 
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return 0;
}
 
 
double func_x( double& x )
{ 
    return x * x;
} 
 







