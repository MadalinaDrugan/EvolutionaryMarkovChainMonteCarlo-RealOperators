#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <iostream.h>

#include "PARAMETERS.h"

#ifndef RANDOM_FILE_
#define RANDOM_FILE_

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1()
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2()
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3()
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53() 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/*set seed on time clock*/
void setSeed()
{ int i;
  time_t t1;

  (void) time(&t1);
  
   srand((unsigned long) t1); /* use time in seconds to set seed */
   
   init_genrand(t1);
}

//box muller
//return a value normal distribution from this mean with this standard deviation
double gaussian()
{
    double number;
    static int iset = 0;
    static float gset;
    double v1, v2, rsq = 0, fac;

        if(iset ==0){
      //here it is a truncation of the generated individuals to fit in the 3 * searcnig space
	  //do{ 
	    do{
	      v1= 2.0*genrand_real1()-1.0; 
	      v2= 2.0*genrand_real1()-1.0;
	      rsq=v1*v1+v2*v2;
	      //cout << "rsq(" << rsq << "," << v1 << "," << v2 << ") \n";
	    } while (rsq > 1.0 || rsq < 0.000000000001);
	    
	    fac = sqrt( log(rsq) *(-2.0 / rsq));
	    gset = fac * v1;
	    iset = 1;
	    number = v2 * fac;
	    //cout << "{" << number << "," << gset << ",{" << rsq << "," << v1 << "," << v2 << "," << fac <<"}}\t";
	    // }while(number >= 1.0 || number <= -1.0 || gset >= 1.0 || gset <= - 1.0); 
	} else { 
	iset = 0;
	number = gset;
	//cout << "{"<< gset <<"}\n";
	}
    //cout << "accepted {" << number << "," << gset << ",{" << rsq << "," << v1 << "," << v2 << "," << fac <<"}}\n";
    //make it a cycle on 0..1
    //if((number * std) + mean < 0 || (number * std) + mean >= 1)
    //return 1 - ((number * std) + mean);
    //else 
    return number;
}


//#ifdef MIXTURE_REALS
//scale - convert from 0..1 interval to 0..scale interval 
//convert to binary
void convertBinary(double* reals, unsigned long* integers, long size){
  
    //cout << "Convert reals from ";
    // for(int k =0; k < size; k++)
    // cout << reals[k] << " \t";
    // cout << "\n";

  //initialization
  for(short i = 0; i < size*BLOCKsize; i++)
      integers[i] = 0;

  for(unsigned long k = 0; k < size; k++){
      unsigned long tempReal = (unsigned long)((reals[k]/scale) * pow(2.0,(double)BLOCKsize));
      //cout << "for k = " << k << " -> " << reals[k]/scale << " conv " << tempReal << "\n";
      unsigned long t = 0, rest;
      while(tempReal != 0 && t < BLOCKsize){
	  rest = tempReal % 2;
	  integers[(k+1)*BLOCKsize-t-1] = rest;
	  //  cout << " rest " << rest << " at pos " << (k+1)*BLOCKsize - t-1 
	  //   << " ramas "<< (tempReal - rest)/2 << "\n";
	  t++;
	  tempReal = (tempReal-rest) / 2;
      }
  }

  //cout << "\n to binary ";
  //for(int j = 0; j < size*BLOCKsize; j++)
  //    cout << integers[j];
  //    cout << "\n";
 }

//convert to reals
 void convertReal(double* reals, unsigned long* integers, long size){
     /////////test convert reals
     // cout << "from binary ";
     //for(int j = 0; j < size*BLOCKsize; j++)
     //	 cout << integers[j];
     //   cout << "\n";

  for(unsigned long k = 0; k < size; k++){
      double tempReal = 0;
      for(int t = 0; t < BLOCKsize; t++){
	  //if(integers[k*BLOCKsize + t] != 0)
	  tempReal += integers[k*BLOCKsize + t] * pow(2.0,(double)BLOCKsize-t-1);
	  //tempReal += integers[(k+1)*BLOCKsize - t-1] * pow(1.0/2.0,BLOCKsize);
	  //   cout << " tempreal " << tempReal << " at pos " << (k+1)*BLOCKsize - t-1 << " from "<< t 
	  // << " and interger" << integers[(k+1)*BLOCKsize - t-1] << "\n";
      }
      reals[k] =(double)tempReal/pow(2.0,(double)BLOCKsize);
  }
  
  //cout << "Convert reals to ";
  //for(int k =0; k < size; k++)
  //    cout << scale*reals[k] << " \t";
  //cout << "\n";
  //scale the reals from [0,1) to [0,2)
  for(int k =0; k < size; k++)
      reals[k] *= scale;
 }

//#endif //mixture reals

#endif

