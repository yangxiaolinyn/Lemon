#include <iostream> 
#include <ctime>
#include <ratio>
#include <chrono>
#include <assert.h>
#include <cmath>
#include <math.h>
#include <iomanip>
using namespace std;
using std::chrono::system_clock;

double const zero = 0.0L;
double const one = 1.0L;
double const two = 2.0L;
double const three = 3.0L;
double const four = 4.0L;
double const five = 5.0L;
double const six = 6.0L;
double const seven = 7.0L;
double const eight = 8.0L;
double const nine = 9.0L;
double const ten = 10.0L;

double const pi = 3.14159265358979323846264380L;
double const twopi = 3.14159265358979323846264380L * 2.0L;
 
struct RNumGenerate
{
private: 
    long rand_inst = 0; 
    long Rand_Feedback = 1;
    double U[98], C, CD, CM, S, T;
    long I97, J97;


public: 
    void initRandom( const long& i = -1, const long& i2 = 9373 )
    {  
        long seed_in, kl, ij;
        string fred;
        //double klr;
 
        seed_in = i;
        //cout << i2 << endl;

        if ( seed_in != -1 )
        {
            kl = i2; 
            if (i2 > 30081)
                cout << "initRandom:second seed too large" << endl ;
            ij = i;
        }
        else
        { 
            //system_clock(count=ij);
            auto begin = std::chrono::steady_clock::now();
            system_clock::time_point today = system_clock::now();
            std::time_t tt;
            tt = system_clock::to_time_t ( today );
            //cout << "today is: " << ctime(&tt) << endl ;
            //cout << "today is: " << tt << endl ; 
            ij = (tt + rand_inst*100) % 31328;
 
            auto end = std::chrono::steady_clock::now();
            auto diff = (end - begin).count();//end-begin得到一个duration类型
            //cout<<  "diff is: " <<  diff << endl;

            //date_and_time( time = fred );
            //read (fred,'(e10.3)') klr;
            kl = ( int(diff*1000) ) % 30081;
            //cout << "ij is: " << ij << endl ;
            //cout << "kl is: " << kl << endl ;
        };
       
        RMARIN( ij, kl );
    };

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    void RMARIN( const long& IJ, const long& KL )
    { 
        if( (IJ < 0) || (IJ > 31328) || (KL < 0) || (KL > 30081) )
        { 
            cout << "The first random number seed must have a value between 0 and 31328." << endl;
            cout << "The second seed must have a value between 0 and 30081. " << endl;
            cout << "T " << IJ << "  " <<  KL << "  " << endl;
            exit(100);
        };

        long I, J, K, L, M;
        I = ( IJ/177) % 177 + 2;
        J = IJ % 177 + 2;
        K = (KL/169) % 178 + 1;
        L = KL % 169;

        //double S, T;
        for(int II = 1; II < 98; II++)
        {
            S = 0.0;
            T = 0.5;
            for(int JJ = 1; JJ < 25; JJ++)
            {  
                M = ( ( (I*J) % 179 )*K ) % 179;
                I = J;
                J = K;
                K = M;
                L = (53*L+1) % 169;
                if ( ( (L*M) % 64 ) > 31 ) 
                    S += T;
                T = 0.5 * T;
            }; 
            U[II] = S; 
        };
        C = 362436.0 / 16777216.0;
        CD = 7654321.0 / 16777216.0;
        CM = 16777213.0 /16777216.0;
        I97 = 97;
        J97 = 33;
        //U[0] = zero; 
        //for( int i=1; i < 98; i++)
        //   cout << U[i] << endl;
    };

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    double RANMAR()
    {
    // This is the random number generator proposed by George Marsaglia in
    // Florida State University Report: FSU-SCRI-87-50
    // It was slightly modified by F. James to produce an array of pseudorandom
    // numbers.
        //double U[97], C, CD, CM;
        //long I97, J97;
        double UNI;
 
        //     INTEGER IVEC
        UNI = U[I97] - U[J97];
        if( UNI < 0.0 )
            UNI += 1.0;
        U[I97] = UNI;
        I97 = I97 - 1;
        if(I97 == 0)
            I97 = 97;
        J97 = J97 - 1;
        if(J97 == 0)
            J97 = 97;
        C = C - CD;
        if( C < 0.0 )
            C = C + CM;
        UNI = UNI - C;
        if( UNI < 0.0 )
            UNI += 1.0; // bug?
        if( UNI == 0.0 || UNI == 1.0 )
            UNI = 0.2459428724394560; 
        return UNI;
    };
};



