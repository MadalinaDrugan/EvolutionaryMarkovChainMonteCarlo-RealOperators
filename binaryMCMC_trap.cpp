#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
                                     // newmatap.h will get include.h
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <math.h>

#include "MatrixUtil\newmatap.h"                // need matrix applications
#include "MatrixUtil\newmatio.h"                // need matrix output routines

#include "PARAMETERS.h"
#include "binaryMCMC_trap.h"
//#include "binaryPopulation_trap.h"
//#include "random.cpp"

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

//extern Random;

extern long unsigned int genrand_int32();
extern double genrand_real2();
extern double genrand_real1();
//extern void  setSeed();
extern void init_genrand(unsigned long);
extern double gaussian();
extern void convertReal(double*,unsigned long*, long);
extern void convertBinary(double*,unsigned long*, long);
//extern static double** positionMixture;
//extern static double*** covarianceMixture;

#ifdef BINOMIAL
const int binomial_table[2 * BLOCKsize - 1] = {1,4,6,4,1};
#ifdef BERNOULLI
const float bernoulli_table[2*BLOCKsize - 1] = {0.4,0.03,0.01,0.03,0.3};
#endif
#endif

#ifndef MIXTURE_REALS
binaryMCMC::binaryMCMC(unsigned long mySize, int* myGenome, double myTemperature){
    //genrand_real1 = new UniformDistribution(0,1);

    index = 0;
	size = mySize;

	sampleSize = sampleSize_Const;

#ifdef LIST
	next = NULL;
#else
	nextValues = new double[sampleSize_Const + 4];
	nextGENOME = new int*[sampleSize_Const + 4];
	for(unsigned long i = 0; i < sampleSize_Const + 4; i++)
	  nextGENOME[i] = new int[size];
	nextTemperatures = new double[sampleSize_Const + 4];
#endif //LIST

	runs = 0;
	genome = new int[size];
	for(unsigned long i=0;i<size;i++)
		genome[i] = myGenome[i]; 
	temperature = myTemperature;	
	value = fitness();
	nrOfGenerations = 0;
	recombination = RECOMBINATION_PROB;
#ifdef MULTIPLE_RUNS
#ifndef MUTATION_VARIES
#ifdef IsRecombination
	PROB_nPOINT = 1;
	UNIF_RECOMB_PARAM = 0.5;
#endif //IsRecombination
#else
	//	mutation = (1.0 / (double) size) * mutation_multi;

#endif //#ifdef MUTATION_VARIES
#else
	//mutation = (1.0 / (double) size) * mutation_multi;
#endif //MULTIPLE_RUNS

	accepted = 0;

#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	proposalDistributionVector = new double*[size];
	for(unsigned long i=0; i < size; i++){
	  proposalDistributionVector[i] = new double[2];
	  for(unsigned long j = 0; j < 2; j++)
	    proposalDistributionVector[i][j] = 1.0;
	}
#endif //Metropolis_Algorithm
}

binaryMCMC::binaryMCMC(unsigned long mySize, int* myGenome, 
		       double myTemperature, double myValue){
    //genrand_real1 = new UniformDistribution(0,1);
	index = 0;
	size = mySize;

	sampleSize = sampleSize_Const;

#ifdef LIST
	next = NULL;
#else
	nextValues = new double[sampleSize_Const + 4];
	nextGENOME = new int*[sampleSize_Const + 4];
	for(unsigned long i = 0; i < sampleSize_Const + 4; i++)
	  nextGENOME[i] = new int[size];
	nextTemperatures = new double[sampleSize_Const + 4];
#endif //LIST

	runs = 0;
	genome = new int[size];
	for(unsigned long i=0;i<size;i++)
		genome[i] = myGenome[i]; 
	temperature = myTemperature;	
	value = myValue;
	nrOfGenerations = 0;
	recombination = RECOMBINATION_PROB;
#ifdef MULTIPLE_RUNS
#ifndef MUTATION_VARIES
#ifdef IsRecombination
	PROB_nPOINT = 1;
	UNIF_RECOMB_PARAM = 0.5;
#endif //IsRecombination
#else
	//	mutation = (1.0 / (double) size) * mutation_multi;

#endif //#ifdef MUTATION_VARIES
#else
	//	mutation = (1.0 / (double) size) * mutation_multi;
#endif //MULTIPLE_RUNS

	accepted = 0;
#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	proposalDistributionVector = new double*[size];
	for(unsigned long i=0; i < size; i++){
	  proposalDistributionVector[i] = new double[2];
	  for(unsigned long j = 0; j < 2; j++)
	    proposalDistributionVector[i][j] = 1.0;
	}
#endif //Metropolis_algorithm
}


binaryMCMC::binaryMCMC(unsigned long mySize, double myTemperature, int myRuns){
    //genrand_real1 = new UniformDistribution(0,1);
	index = 0;
	size = mySize;

	sampleSize = sampleSize_Const;

#ifdef LIST
	next = NULL;
#else
	nextValues = new double[sampleSize_Const + 4];
	nextGENOME = new int*[sampleSize_Const + 4];
	for(unsigned long i = 0; i < sampleSize_Const + 4; i++)
	  nextGENOME[i] = new int[size];
	nextTemperatures = new double[sampleSize_Const + 4];
#endif //LIST

	runs = myRuns;
	genome = new int[size];

#ifdef NORMAL_RUNS 
	for(unsigned long i=0;i<size;i++){
	    genome[i] = genrand_int32() % 2;}
#else
#ifdef BINOMIAL
    for(int j = 0; j < (size/(double)BLOCKsize); j++){
	for(unsigned long i = 0; i<BLOCKsize/2; i++){
	    genome[i + j*BLOCKsize] = 0;
	    //   cout << genome[i + j*BLOCKsize];
	}
	for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++){
	    genome[i + j*BLOCKsize] = 1;
	    //cout << genome[i + j*BLOCKsize];
	}
    }
//	  if(i % BLOCKsize >= BLOCKsize/2)
//	      genome[i] = 1;
//	  else genome[i] = 0;
#else
#ifdef TWO_ATTRACTORS_FUNCTION
	double count1 = 0;
	for(unsigned long i=0;i<size;i++){
	  if(i - count1 >= size/2) {
	      genome[i] = 1;
	      count1++;
	  } else if(count1 <= size/2){
	      //srand(time(0));
	      genome[i] = (unsigned long)(genrand_int32() % 2);
	      //genome[i] = genrand_int32() % 2;
	      if(genome[i] == 1)
		  count1++;
	  }
	  else genome[i] = 0;  
	}
#else
	for(unsigned long i=0;i<size;i++){
	    //if(runs == 0) 
	    genome[i] = 0;
	    /*else if(runs == 1) genome[i] = 1;
	  else if(runs == 2) { 
	    if(i < size/2) genome[i] = 0;
	    else genome[i] = 1; 
	  } else if(runs == 3) {
	    if(i < size/3) genome[i] = 0;
	    else genome[i] = 1;
	  } else if(runs == 4) {
	    if(i < size/4) genome[i] = 0;
	    else genome[i] = 1;
	  } else if(runs == 5) {
	    if(i < BLOCKsize) genome[i] = 0;
	    else genome[i] = 1;
	  }
	  else genome[i] = genrand_int32() % 2;*/
	}
#endif
#endif
#endif
//	  cout << genome[i] << "||";
//	cout << "\n";	

	temperature = myTemperature;	
	value = fitness();
       	nrOfGenerations = 0;
	recombination = RECOMBINATION_PROB;
#ifdef MULTIPLE_RUNS
#ifndef MUTATION_VARIES
#ifdef IsRecombination
	PROB_nPOINT = 1;
	UNIF_RECOMB_PARAM = 0.5;
#endif //IsRecombination
#else
	//mutation = (1.0 / (double) size) * mutation_multi;

#endif //#ifdef MUTATION_VARIES
#else
	//mutation = (1.0 / (double) size) * mutation_multi;
#endif //MULTIPLE_RUNS

	accepted = 0;
#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	proposalDistributionVector = new double*[size];
	for(unsigned long i=0; i < size; i++){
	  proposalDistributionVector[i] = new double[2];
	  for(long j = 0; j < 2; j++)
	    proposalDistributionVector[i][j] = 1.0;
	}
#endif //Metropolis_Algorithm
}


#ifdef SINGLE_vs_MULTIPLE
binaryMCMC::binaryMCMC(unsigned long mySize, double myTemperature, 
		       int nr_chain, int sizeMCMCstring){
    //genrand_real1 = new UniformDistribution(0,1);
	index = 0;
	size = mySize;

#ifdef LIST
	next = NULL;
#else

	sampleSize = sampleSize_Const*((int)(log(sizeMCMCstring)/log(2.0)) +1); 

	nextValues = new double[sampleSize + 4];
	nextGENOME = new int*[sampleSize + 4];
	for(unsigned long i = 0; i < sampleSize + 4; i++)
	  nextGENOME[i] = new int[size];
	nextTemperatures = new double[sampleSize + 4];
#endif //LIST

	runs = 0;
	genome = new int[size];
	double count1 = 0;

	for(unsigned long i=0;i<size;i++){
#ifdef NORMAL_RUNS
	  genome[i] = genrand_int32() % 2;
#else
#ifdef BINOMIAL
	  if(i % BLOCKsize >= BLOCKsize/2)
	      genome[i] = 1;
	  else genome[i] = 0;
#else
#ifdef TWO_ATTRACTORS_FUNCTION
	  if(i - count1 >= size/2) {
	      genome[i] = 1;
	      count1++;
	  } else if(count1 <= size/2){
	      //srand(time(0));
	      genome[i] = (unsigned long)(genrand_int32() % 2);
	      //genome[i] = genrand_int32() % 2;
	      if(genome[i] == 1)
		  count1++;
	  }
	  else genome[i] = 0;  
#else
	  genome[i] = 0;
#endif
#endif
#endif
	  cout << genome[i];
	}
	cout << "\n";	

	temperature = myTemperature;	
	value = fitness();
       	nrOfGenerations = 0;
	recombination = RECOMBINATION_PROB;
#ifdef MULTIPLE_RUNS
#ifndef MUTATION_VARIES
#ifdef IsRecombination
	PROB_nPOINT = 1;
	UNIF_RECOMB_PARAM = 0.5;
#endif //IsRecombination
#else
	//mutation = (1.0 / (double) size) * mutation_multi;

#endif //#ifdef MUTATION_VARIES
#else
	mutation = (1.0 / (double) size) * mutation_multi;
#endif //MULTIPLE_RUNS

	accepted = 1;
#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	proposalDistributionVector = new double*[size];
	for(unsigned long i=0; i < size; i++){
	  proposalDistributionVector[i] = new double[2];
	  for(long j = 0; j < 2; j++)
	    proposalDistributionVector[i][j] = 1.0;
	}
#endif //Metropolis_Algorithm

}
#endif

#ifdef QBP
binaryMCMC::binaryMCMC(unsigned long mySize, double myTemperature, double** qbp){
    //genrand_real1 = new UniformDistribution(0,1);
	index = 0;
	size = mySize;

	sampleSize = sampleSize_Const;

#ifdef LIST
	next = NULL;
#else
	nextValues = new double[sampleSize_Const + 4];
	nextGENOME = new int*[sampleSize_Const + 4];
	for(unsigned long i = 0; i < sampleSize_Const + 4; i++)
	  nextGENOME[i] = new int[size];
	nextTemperatures = new double[sampleSize_Const + 4];
#endif //LIST

	runs = 0;
	genome = new int[size];
	double count1 = 0;

#ifdef NORMAL_RUNS
	for(unsigned long i=0;i<size;i++){
	  genome[i] = genrand_int32() % 2;
	}
#else
#ifdef BINOMIAL
     for(int j = 0; j < (size/(double)BLOCKsize); j++){
	 for(unsigned long i = 0; i<BLOCKsize/2; i++){
	     //if(chain %2 == 0)
		 genome[i + j*BLOCKsize] = 0;
		 //else 
		 //genome[i + j*BLOCKsize] = 1;
	     cout << genome[i+ j*BLOCKsize] << "||";
	 }
	 for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++){
	     //if(chain %2 == 0)
		 genome[i + j*BLOCKsize] = 1;
		 //else 
		 //genome[i + j*BLOCKsize] = 0;
	     cout << genome[i+ j*BLOCKsize] << "||";
	 }
     }
     cout << " = "<< fitness() << "\n";	
#else
#ifdef TWO_ATTRACTORS_FUNCTION
	for(unsigned long i=0;i<size;i++){
	  if(i - count1 >= size/2) {
	      genome[i] = 1;
	      count1++;
	  } else if(count1 <= size/2){
	      //srand(time(0));
	      //genome[i] = (unsigned long)(rand() % 2);
	      genome[i] = genrand_int32() % 2;
	      if(genome[i] == 1)
		  count1++;
	  }
	  else genome[i] = 0;  
	  cout << genome[i] << "||";
	}
	cout <<"\n";	
#else
	for(unsigned long i=0;i<size;i++){
	  genome[i] = 0;
	}
#endif
#endif
#endif
	fitness();	

	temperature = myTemperature;	
	value = fitness();
       	nrOfGenerations = 0;
	mutation = (1.0 / (double) size) * mutation_multi;
	recombination = RECOMBINATION_PROB;
#ifdef MULTIPLE_RUNS
#ifdef IsRecombination
	PROB_nPOINT = 1;
	UNIF_RECOMB_PARAM = 0.5;
#endif //IsRecombination
#endif //MULTIPLE_RUNS
	accepted = 0;
#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	proposalDistributionVector = new double*[size];
	for(unsigned long i=0; i < size; i++){
	  proposalDistributionVector[i] = new double[2];
	  for(long j = 0; j < 2; j++)
	    proposalDistributionVector[i][j] = 1.0;
	}
#endif //Metropolis_Algorithm

	if(qbp != NULL){
	    matrixQBP = new double*[size];
	    for(int i = 0; i < size; i++){
		matrixQBP[i] = new double[size];
		for(int j = 0; j < size; j++)
		    matrixQBP[i][j] = qbp[i][j];
	    }
	}else{
	    //genereaza matricea
	    cout << "this is not correct: matrix QBP uninitialised in creator stringMCMC\n";
	    exit(1);
	    for(int i = 0; i < size; i++)
		for(int j = 0; j < size; j++)
		    matrixQBP[i][j] =  UniformDistribution(0,1);
	    
	}
}
#endif //QBP

binaryMCMC::binaryMCMC(unsigned long mySize, double myTemperature, unsigned long chain){
    //genrand_real1 = new UniformDistribution(0,1);
	index = 0;
	size = mySize;

	sampleSize = sampleSize_Const;

#ifdef LIST
	next = NULL;
#else
	nextValues = new double[sampleSize_Const + 4];
	nextGENOME = new int*[sampleSize_Const + 4];
	for(unsigned long i = 0; i < sampleSize_Const + 4; i++)
	  nextGENOME[i] = new int[size];
	nextTemperatures = new double[sampleSize_Const + 4];
#endif //LIST

	runs = 0;
	genome = new int[size];
	double count1 = 0;

#ifdef NORMAL_RUNS
	for(unsigned long i=0;i<size;i++){
	  genome[i] = genrand_int32() % 2;
	}
#else
#ifdef BINOMIAL
     for(int j = 0; j < (size/(double)BLOCKsize); j++){
	 for(unsigned long i = 0; i<BLOCKsize/2; i++){
	     //if(chain %2 == 0)
		 genome[i + j*BLOCKsize] = 0;
		 //else 
		 //genome[i + j*BLOCKsize] = 1;
	     cout << genome[i+ j*BLOCKsize] << "||";
	 }
	 for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++){
	     //if(chain %2 == 0)
		 genome[i + j*BLOCKsize] = 1;
		 //else 
		 //genome[i + j*BLOCKsize] = 0;
	     cout << genome[i+ j*BLOCKsize] << "||";
	 }
     }
     cout << " = "<< fitness() << "\n";	
#else
#ifdef TWO_ATTRACTORS_FUNCTION
	for(unsigned long i=0;i<size;i++){
	  if(i - count1 >= size/2) {
	      genome[i] = 1;
	      count1++;
	  } else if(count1 <= size/2){
	      //srand(time(0));
	      //genome[i] = (unsigned long)(rand() % 2);
	      genome[i] = genrand_int32() % 2;
	      if(genome[i] == 1)
		  count1++;
	  }
	  else genome[i] = 0;  
	  cout << genome[i] << "||";
	}
	cout <<"\n";	
#else
	for(unsigned long i=0;i<size;i++){
	  genome[i] = 0;
	}
#endif
#endif
#endif
	fitness();	

	temperature = myTemperature;	
	value = fitness();
       	nrOfGenerations = 0;
	recombination = RECOMBINATION_PROB;
#ifdef MULTIPLE_RUNS
#ifndef MUTATION_VARIES
#ifdef IsRecombination
	PROB_nPOINT = 1;
	UNIF_RECOMB_PARAM = 0.5;
#endif //IsRecombination
#else
	//mutation = (1.0 / (double) size) * mutation_multi;

#endif //#ifdef MUTATION_VARIES
#else
	//mutation = (1.0 / (double) size) * mutation_multi;
#endif //MULTIPLE_RUNS

	accepted = 0;
#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	proposalDistributionVector = new double*[size];
	for(unsigned long i=0; i < size; i++){
	  proposalDistributionVector[i] = new double[2];
	  for(long j = 0; j < 2; j++)
	    proposalDistributionVector[i][j] = 1.0;
	}
#endif //Metropolis_Algorithm
}

binaryMCMC::binaryMCMC(unsigned long mySize, double myTemperature){
    //genrand_real1 = new UniformDistribution(0,1);
	index = 0;
	size = mySize;

	sampleSize = sampleSize_Const;

#ifdef LIST
	next = NULL;
#else
	nextValues = new double[sampleSize_Const + 4];
	nextGENOME = new int*[sampleSize_Const + 4];
	for(unsigned long i = 0; i < sampleSize_Const + 4; i++)
	  nextGENOME[i] = new int[size];
	nextTemperatures = new double[sampleSize_Const + 4];
#endif //LIST

	runs = 0;
	genome = new int[size];

#ifdef NORMAL_RUNS
	for(unsigned long i=0;i<size;i++){
	  genome[i] = genrand_int32() % 2;
	}
#else
#ifdef BINOMIAL 
     for(int j = 0; j < (size/(double)BLOCKsize); j++){
	 for(unsigned long i = 0; i<BLOCKsize/2; i++)
	     genome[i + j*BLOCKsize] = 0;
	 for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++)
	     genome[i + j*BLOCKsize] = 1;
     }
#else
#ifdef TWO_ATTRACTORS_FUNCTION
	double count1 = 0;
	for(unsigned long i=0;i<size;i++){
	    if(i - count1 >= size/2) {
		genome[i] = 1;
		count1++;
	    } else if(count1 <= size/2){
		//srand(time(0));
		genome[i] = (unsigned long)(genrand_int32() % 2);
		if(genome[i] == 1)
		    count1++;
	    }
	    else genome[i] = 0;  
	    
	    cout << genome[i] << "||";
	}
	cout << "\n";	
#else
	for(unsigned long i=0;i<size;i++){
	    genome[i] = 0;
	    
	    cout << genome[i] << "||";
	}
	cout << "\n";	
#endif //TWO_Attractors
#endif //binomial
#endif //NORMAL_RUNS

	temperature = myTemperature;	
	value = fitness();
       	nrOfGenerations = 0;
	recombination = RECOMBINATION_PROB;
#ifdef MULTIPLE_RUNS
#ifndef MUTATION_VARIES
#ifdef IsRecombination
	PROB_nPOINT = 1;
	UNIF_RECOMB_PARAM = 0.5;
#endif //IsRecombination
#else
	//	mutation = (1.0 / (double) size) * mutation_multi;

#endif //#ifdef MUTATION_VARIES
#else
	//mutation = (1.0 / (double) size) * mutation_multi;
#endif //MULTIPLE_RUNS

	accepted = 0;
#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	proposalDistributionVector = new double*[size];
	for(unsigned long i=0; i < size; i++){
	  proposalDistributionVector[i] = new double[2];
	  for(long j = 0; j < 2; j++)
	    proposalDistributionVector[i][j] = 1.0;
	}
#endif //Metropolis_Algorithm

}

#else //MIXTURE_REALS
binaryMCMC::binaryMCMC(unsigned long mySize, double myTemperature){
    //genrand_real1 = new UniformDistribution(0,1);
    //UniformDistribution genrand_real1(0,1);
    index = 0;
    size = mySize;
    
    sampleSize = sampleSize_Const;

#ifdef MIXTURE_BIVARIATE
    read_parametri();
#endif //MIXTURE_BIVARIATE
    
    //gaussian = new CDataDistrib();
#ifdef LIST
    next = NULL;
#else
    nextValues = new double[sampleSize_Const + 4];
    nextGENOME = new double*[sampleSize_Const + 4];
    for(unsigned long i = 0; i < sampleSize_Const + 4; i++)
	nextGENOME[i] = new double[size];
    nextTemperatures = new double[sampleSize_Const + 4];
#endif //LIST
    
    runs = 0;
    genome = new double[size];
    //genomeRecomb = new int[size*((int) pow(2,BLOCKsize))];
 
    /*#ifndef BIVARIATE   
//#ifdef NORMAL_RUNS
    unsigned long integers[size*BLOCKsize];
    for(unsigned long i=0;i<size*BLOCKsize;i++){
	integers[i] = genrand_int32() % 2;
	//cout << "\t" << integers[i];
	}
    //cout << "\n";
    convertReal(genome,integers,size);
	//cout << "\n";
//#endif //NORMAL_RUNS
#else*/
    bool true_bit = true;
    while(true_bit){
      cout << "proposed genome \t";
     for(unsigned long i=0;i<size;i++){
#ifndef INIT_HALF
#ifndef INIT_CENTER
#ifndef INIT_SIDE
	genome[i] = genrand_real1() * scale;
	genome[i] = genome[i] - scale/2.0;
#else
	genome[i] = genrand_real1() * scale/init_coef;
	if(i == 0)
	  genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
	else genome[i] = genome[i] - (init_coef - 1.0) *scale/2.0/init_coef;
	//if(i == 0)
	//  genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
	//else genome[i] = genome[i] - (init_coef - 1.0) *scale/2.0/init_coef;
#endif //init side
#else
	genome[i] = genrand_real1() * scale/init_coef;
	genome[i] = genome[i] - scale/2.0/init_coef;
#endif //init_center
#else
      if(i == 0)
	genome[i] = genrand_real1() * scale - scale/2.0;
      else genome[i] = genrand_real1() * scale/2.0;
      //	if(i == 0)
      //  genome[i] = genrand_real1() * scale/2.0;
      //else genome[i] = genrand_real1() * scale;
#endif 
	//if(genome[i] < scale/2.0)
	//  genome[i] = - genome[i];
	//else genome[i] -= scale/2;
	//cout << genome[i] << "\t";
     }
     //cout << "\n";

     if(fitness(genome,size) > MIN_FITNESS){
     for(unsigned long i=0;i<size;i++)
       cout << genome[i] << "\t";
     cout << fitness(genome,size) << "\n";
       true_bit = false;
     } //else cout << " rejected fitness " << fitness(genome,size) << "\n";
    }
    //#endif //bivariate

	temperature = myTemperature;	
	value = fitness();
       	nrOfGenerations = 0;
	recombination = RECOMBINATION_PROB;
#ifdef MULTIPLE_RUNS
#ifndef MUTATION_VARIES
#ifdef IsRecombination
	PROB_nPOINT = 1;
	UNIF_RECOMB_PARAM = 0.5;
#endif //IsRecombination
#else
	mutation = (1.0 / (double) size) * mutation_multi;

#endif //#ifdef MUTATION_VARIES
#else
	mutation = (1.0 / (double) size) * mutation_multi;
#endif //MULTIPLE_RUNS

	accepted = 0;
#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	proposalDistributionVector = new double*[size];
	for(unsigned long i=0; i < size; i++){
	  proposalDistributionVector[i] = new double[2];
	  for(long j = 0; j < 2; j++)
	    proposalDistributionVector[i][j] = 1.0;
	}
#endif //Metropolis_Algorithm
}

binaryMCMC::binaryMCMC(unsigned long mySize, double* myGenome, 
		       double myTemperature, double myValue){
    //UniformDistribution genrand_real1(0,1);
    //genrand_real1 = new UniformDistribution(0,1);
	index = 0;
	size = mySize;

	sampleSize = sampleSize_Const;

	//gaussian = new CDataDistrib();

#ifdef MIXTURE_BIVARIATE
	read_parametri();
#endif //MIXTURE_BIVARIATE

#ifdef LIST
	next = NULL;
#else
	nextValues = new double[sampleSize_Const + 4];
	nextGENOME = new double*[sampleSize_Const + 4];
	for(unsigned long i = 0; i < sampleSize_Const + 4; i++)
	  nextGENOME[i] = new double[size];
	nextTemperatures = new double[sampleSize_Const + 4];
#endif //LIST

	runs = 0;
	genome = new double[size];
	for(unsigned long i=0;i<size;i++)
		genome[i] = myGenome[i]; 
	temperature = myTemperature;	
	value = myValue;
	nrOfGenerations = 0;
	recombination = RECOMBINATION_PROB;
#ifdef MULTIPLE_RUNS
#ifndef MUTATION_VARIES
#ifdef IsRecombination
	PROB_nPOINT = 1;
	UNIF_RECOMB_PARAM = 0.5;
#endif //IsRecombination
#else
	mutation = (1.0 / (double) size) * mutation_multi;

#endif //#ifdef MUTATION_VARIES
#else
	mutation = (1.0 / (double) size) * mutation_multi;
#endif //MULTIPLE_RUNS

	accepted = 0;
#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	proposalDistributionVector = new double*[size];
	for(unsigned long i=0; i < size; i++){
	  proposalDistributionVector[i] = new double[2];
	  for(unsigned long j = 0; j < 2; j++)
	    proposalDistributionVector[i][j] = 1.0;
	}
#endif //Metropolis_algorithm
}

#endif //REALS

binaryMCMC::~binaryMCMC()
{
#ifdef LIST
  if(next != NULL) delete next;
#else
  delete[] nextValues;
  delete[] nextGENOME;
  delete[] nextTemperatures;
#endif

  delete[] genome;
#ifndef METROPOLIS_ALGORITHM
  delete[] proposalDistributionVector;
#endif //Metropolis_algorithm

#ifdef QBP
  if(matrixQBP != NULL) delete[] matrixQBP;
#endif

#ifdef MIXTURE_BIVARIATE
	delete_parametri();
#endif //MIXTURE_BIVARIATE
}


#ifdef TRAP_FUNCTION
//fitness for multiple trap functions
//BLOCKsize - size of a trap function
//size - size of a genome
double binaryMCMC::fitness(){
	double tempValue = 0;
	value = 0;

#ifdef DEBUG_
  cout << "fitness is for(";
  for(unsigned long i = 0; i<size; i++)
    if(i % BLOCKsize == BLOCKsize - 1) cout <<genome[i]<<"-"; 
      else cout <<genome[i]<<",";
  cout<<") == "; 
#endif

  for(unsigned long j = 0; j<(size/BLOCKsize); j++)
    {
      tempValue = 0;
      for(unsigned long i = 0; i<BLOCKsize; i++)
	tempValue += genome[i + j * BLOCKsize];
      if(BLOCKsize == 3){
	  if(tempValue == 3) tempValue = a;
	  else if(tempValue == 2) tempValue = 0;
	  else if (tempValue == 1) tempValue = c;
	  else tempValue = b;
      } else if(BLOCKsize == 4){
	  if(tempValue == 4) tempValue = a;
	  else tempValue = b - tempValue; 
      } else {
	  if(tempValue <= z)  tempValue = (b/z)*(z-tempValue);
	  else tempValue = (a/(BLOCKsize-z))*(tempValue-z);
      }
#ifdef LINEAR_TRAP_FUNCTION
      value += tempValue;	
#ifdef DEBUG_
    cout << "("<<j+1<<","<<value<<")\t";
#endif
#else
      value += pow(j+1,ORDER_BLOCK)*tempValue;
#ifdef DEBUG_
    cout << "("<<j+1<<","<<value<<")\t";
#endif
#endif
    }	
#ifdef DEBUG_
  cout << "\n";
#endif	
  if(value == 0) value = AwayFromZero;
  return value;
}

double binaryMCMC::fitness(int* tempGenome, unsigned long mySize){
	double tempValue = 0; 
	double myValue = 0;

#ifdef DEBUG_
  cout << "fitness is for(";
  for(unsigned long i = 0; i<size; i++)
    if(i % BLOCKsize == BLOCKsize - 1) cout <<tempGenome[i]<<"-"; 
      else cout <<tempGenome[i]<<",";
  cout<<") == "; 
#endif

	for(unsigned long j = 0; j<(mySize / BLOCKsize); j++) 
	{
		tempValue = 0;
		for(unsigned long i = 0; i<BLOCKsize; i++) 
			tempValue += tempGenome[i + j * BLOCKsize];
		if(BLOCKsize == 3){
		    if(tempValue == 3) tempValue = a;
		    else if(tempValue == 2) tempValue = 0;
		    else if (tempValue == 1) tempValue = c;
		    else tempValue = b;
		} else if(BLOCKsize == 4){
		    if(tempValue == 4) tempValue = a;
		    else tempValue = b - tempValue;
		} else {
		    if(tempValue <= z)  tempValue = (b/z)*(z-tempValue);
		    else tempValue = (a/(BLOCKsize-z))*(tempValue-z);
		}		 
#ifdef LINEAR_TRAP_FUNCTION
		myValue += tempValue;	
#ifdef DEBUG_
    cout << "("<<j+1<<","<<myValue<<")\t";
#endif
#else
		myValue += pow(j+1,ORDER_BLOCK)*tempValue;
#endif
	}	
#ifdef DEBUG_
  cout << "\n";
#endif		
  if(myValue == 0) return AwayFromZero;
  return myValue;
}

int binaryMCMC::stop(){
#ifdef  LINEAR_TRAP_FUNCTION
	if(value == a* (size / BLOCKsize)) return 1;
	if(value == b* (size / BLOCKsize)) return 2;	
#else
	unsigned long sizeM = size / BLOCKsize;
	if(value == a * (sizeM * (sizeM + 1 ) / 2)) return 1;
	if(value == b * (sizeM * (sizeM + 1 ) / 2)) return 2;
#endif
	return 0;
}

void binaryMCMC::getMax(int* maxGenome){
    for(int i = 0; i < size; i++)
	maxGenome[i] = 1;
}

void binaryMCMC::getMax0(int* maxGenome){
    for(int i = 0; i < size; i++)
	maxGenome[i] = 0;
}

#ifdef RECOMB_MUT
void binaryMCMC::setMax0(){
    for(int i = 0; i < size; i++)
	genome[i] = 0;
}

void binaryMCMC::setMax1(){
    for(int i = 0; i < size; i++)
	genome[i] = 1;
}
#endif //RECOMB_MUT

#else 
#ifdef ONEMAX_FUNCTION
double binaryMCMC::fitness(){
	double tempValue = 0;
	value = 0;
	for(unsigned long i = 0; i<size; i++)
		tempValue += genome[i];
	value = tempValue;	
	return value;
}

double binaryMCMC::fitness(int* tempGenome, unsigned long mySize){
	double tempValue = 0; 
	for(unsigned long i = 0; i<mySize; i++) 
		tempValue += tempGenome[i];
	return tempValue;
}

int binaryMCMC::stop(){
	if(value == size) return 1;
	return 0;
}
#else
#ifdef MULTI_PICK_TRAP_FUNCTION
//fitness for 1 trap function with 2 picks
double binaryMCMC::fitness(int* tempGenome, unsigned long mySize){
	double tempValue = 0; 
	double myValue = 0;
	int tempBLOCKsize = BLOCKsize;

	if(BLOCKsize%2 == 1) tempBLOCKsize = BLOCKsize - 1;

	for(unsigned long j = 0; j<(mySize / BLOCKsize); j++){
	  tempValue = 0;

	  for(unsigned long i = 0; i<BLOCKsize; i++) 
		tempValue += tempGenome[i + j*BLOCKsize];

	  if(tempValue <= z)  tempValue = (b/z)*(z-tempValue);
	  else 
	      if(tempValue < tempBLOCKsize/2) 
		  tempValue = (a/(tempBLOCKsize/2 - z))*(tempValue-z);		 
	      else 
		  if(BLOCKsize - tempValue > z) 
		      tempValue = (a/(tempBLOCKsize/2 - z))*(tempBLOCKsize - tempValue-z);
		  else tempValue = (b/z)*(z-tempBLOCKsize+tempValue);
	  
	  myValue += tempValue;
	}

	return myValue;
}

double binaryMCMC::fitness(){
  double tempValue = 0; 
#ifdef DEBUG_
  cout << "fitness is for(";
  for(unsigned long i = 0; i<size; i++)
    if(i % BLOCKsize == BLOCKsize - 1) cout <<genome[i]<<"-"; 
      else cout <<genome[i]<<",";
  cout<<") == "; 
#endif
  value = 0;

  int tempBLOCKsize = BLOCKsize;
  if(BLOCKsize%2 == 1) tempBLOCKsize = BLOCKsize- 1;

  for(unsigned long j = 0; j<(size / BLOCKsize); j++){
    tempValue = 0;
    for(unsigned long i = 0; i<BLOCKsize; i++)
      tempValue += genome[i + j*BLOCKsize];
 
    if(tempValue <= z)  tempValue = (b/z)*(z-tempValue);
    else if(tempValue < tempBLOCKsize/2) tempValue = (a/(tempBLOCKsize/2-z))*(tempValue-z);		 
    else if(BLOCKsize - tempValue > z) tempValue = (a/(tempBLOCKsize/2 - z))*(tempBLOCKsize - tempValue-z);
    else tempValue = (b/z)*(z-tempBLOCKsize+tempValue);
    value += tempValue;
#ifdef DEBUG_
    cout << "("<<j<<","<<value<<")\t";
#endif
  }	
#ifdef DEBUG_
  cout << "\n";
#endif
  return value;
}

int binaryMCMC::stop(){
	if(value == a * size / BLOCKsize) return 1;
	if(value == b * size / BLOCKsize) return 2;
	return 0;
}

void binaryMCMC::getMax(int* maxGenome){
    for(int i = 0; i < size; i++)
	maxGenome[i] = 1;
}

#else
#ifdef TWO_ATTRACTORS_FUNCTION
//fitness 2 attractors and 1 platou between them
double binaryMCMC::fitness(int* tempGenome, unsigned long mySize){
  double tempValue = 0; 
  double myValue = 0;
  int zero = 0;
  int one = 0;
  for(unsigned long i = 0; i<mySize; i++){
    if(tempGenome[i] == 0) zero++;
    else 
      if(tempGenome[i] == 1) one++;
      else {
	  cout << "Genome eronat tempGenome[" <<i<<"]="
	       <<tempGenome[i]<<"in binaryMCMC_trap.cpp:467\n";
	  exit(1);
      }
  }

  if(zero >= wrongBits && one >= wrongBits) return AwayFromZero;
  
  int prim = 0;
  if(wrongBits > one) {
      /*for(unsigned long i = 0; i<mySize; i++)
	if(tempGenome[i] == 0 && (prim < one && prim != 0)) return 0;
	else if(tempGenome[i] == 1) prim++;*/
    tempValue = wrongBits - one;
  }
  if(wrongBits > zero) {
      /*for(unsigned long i = 0; i<mySize; i++)
	if(tempGenome[i] == 1 && (prim < zero && prim != 0)) return 0;
	else if(tempGenome[i] == 0) prim++; */ 
      tempValue = wrongBits - zero;
  } 
#ifdef DEBUG_
  if(tempValue > 0){
    cout << "Greater then 0";
    for(unsigned long i = 0; i<size; i++)
      cout <<tempGenome[i];    
    cout << tempValue << "\n";
  }
#endif
  return tempValue;
}

double binaryMCMC::fitness(){
#ifdef DEBUG_
  cout << "fitness is for(";
  for(unsigned long i = 0; i<size; i++)
       cout <<genome[i]<<",";
  cout<<") == "; 
#endif
  int zero = 0;
  int one = 0;
  for(unsigned long i = 0; i<size; i++){
    if(genome[i] == 0) zero++;
    else 
      if(genome[i] == 1) one++;
      else {
	  cout << "Genome eronat binaryMCMC_trap.cpp:648\n";
	  exit(1);
      }
  }

  if(zero >= wrongBits && one >= wrongBits) {
    value = AwayFromZero;
    return value;
  }

  int prim = 0;
  if(wrongBits > one) {
    /*for(unsigned long i = 0; i<size; i++)
      if(genome[i] == 0 && (prim < one && prim != 0)) return 0;
	else if(genome[i] == 1) prim++;*/
    value = wrongBits - one;
  }
  if(wrongBits > zero) {
    /*for(unsigned long i = 0; i<size; i++)
      if(genome[i] == 1 && (prim < zero && prim != 0)) return 0;
      else if(genome[i] == 0) prim++;*/ 
    value = wrongBits - zero;
  }   	
#ifdef DEBUG_
  if(value > 0){
    cout << "Greater then 0";
    for(unsigned long i = 0; i<size; i++)
      cout <<genome[i];    
    cout << value << "\n";
  }
#endif
  return value;
}

int binaryMCMC::stop(){
  int zero = 0;
  int one = 0;
  for(unsigned long i = 0; i<size; i++){
    if(genome[i] == 0) zero++;
    else 
      if(genome[i] == 1) one++;
      else {
	  cout << "Genome eronat \n";
	  exit(1);
      }
  }
  if(zero == size) return 1;
  else if(one == size) return 2;
  else return 0;
}

int binaryMCMC::stop(int* tempGenome, unsigned long mySize){
  int zero = 0;
  int one = 0;
  for(unsigned long i = 0; i<mySize; i++){
    if(tempGenome[i] == 0) zero++;
    else 
      if(tempGenome[i] == 1) one++;
      else {
	  cout << "Genome eronat \n";
	  exit(1);
      }
  }
  if(zero == size) return 1;
  else if(one == size) return 2;
  else return 0;
}

void binaryMCMC::getMax(int* maxGenome){
    for(int i = 0; i < size; i++)
	maxGenome[i] = 1;
}

void binaryMCMC::getMax0(int* maxGenome){
    for(int i = 0; i < size; i++)
	maxGenome[i] = 0;
}

#else
#ifdef BINOMIAL
#ifdef BERNOULLI
double binaryMCMC::fitness(int* tempGenome, unsigned long mySize){
  int tempValue = 0; 
  double myValue = 0;

  for(unsigned long j = 0; j<(mySize / BLOCKsize); j++){
    tempValue = 0;

    for(unsigned long i = 0; i<BLOCKsize; i++) 
      tempValue += tempGenome[i + j*BLOCKsize];
    
    myValue += binomial_table[tempValue] * pow(p_binomial,tempValue) * pow(1 - p_binomial,tempValue);
  }

  return myValue;
}

double binaryMCMC::fitness(){
  int tempValue = 0; 
  value = 0;

  for(unsigned long j = 0; j<(size / BLOCKsize); j++){
    tempValue = 0;

    for(unsigned long i = 0; i<BLOCKsize; i++)
      tempValue += genome[i + j*BLOCKsize];
 
    value += binomial_table[tempValue] * pow(p_binomial,tempValue) * pow(1 - p_binomial,tempValue);
  }	
  return value;
}

int binaryMCMC::stop(){
	if(value == size / BLOCKsize) return 1;
	return 0;
}

void binaryMCMC::getMax(int* maxGenome){
    for(int i = 0; i < size; i++)
	maxGenome[i] = 1;
}

void binaryMCMC::getMax0(int* maxGenome){
    for(int i = 0; i < size; i++)
	maxGenome[i] = 0;
}

#else
//#ifndef BERNOULLI
// are doua pickuri - maxim la 11..10..00, unde 11..1 sint jumatate din 
//lungimea totala a stringului
// al doilea pick 00..011..1 unde 11..1 sint din nou jumatate din string
double binaryMCMC::fitness(int* tempGenome, unsigned long mySize){
  int tempValue = 0; 
  double myValue = 0;

  for(unsigned long j = 0; j<(mySize / BLOCKsize); j++){
    tempValue = 0;

    for(unsigned long i = 0; i<BLOCKsize/2; i++) 
      if(tempGenome[i + j*BLOCKsize] == 0) tempValue ++;

    for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++) 
      if(tempGenome[i + j*BLOCKsize] == 1) tempValue ++;
    
    //myValue += bernoulli_table[tempValue];
    if(tempValue < BLOCKsize/2) tempValue = 2*(BLOCKsize/2 - tempValue);
    else if(tempValue == BLOCKsize/2) tempValue = 0;
    else if(tempValue <= BLOCKsize) tempValue = 2*(tempValue - BLOCKsize/2)-1;
    else {
	cout << "wrong fitness value" << tempValue << " see binaryMCMC_trap.cpp:613\n";
	exit(1);
    }

#ifdef LINEAR_TRAP_FUNCTION
		myValue += tempValue;	
#ifdef DEBUG_
    cout << "("<<j+1<<","<<myValue<<")\t";
#endif
#else
		myValue += pow(j+1,ORDER_BLOCK)*tempValue;
#endif

  }
  
  if(myValue == 0) myValue = AwayFromZero;
  return myValue;
}

double binaryMCMC::fitness(){
  int tempValue = 0; 
  value = 0;

  for(unsigned long j = 0; j<(size / BLOCKsize); j++){
    tempValue = 0;

    for(unsigned long i = 0; i<BLOCKsize/2; i++)
      if(genome[i + j*BLOCKsize] == 0) tempValue++;
 
     for(unsigned long i = BLOCKsize/2; i < BLOCKsize; i++)
      if(genome[i + j*BLOCKsize] == 1) tempValue++;

     //value += bernoulli_table[tempValue];
     if(tempValue < BLOCKsize/2) tempValue = 2*(BLOCKsize/2 - tempValue);
     else if(tempValue == BLOCKsize/2) tempValue = 0;
     else if(tempValue <= BLOCKsize) tempValue = 2*(tempValue - BLOCKsize/2) - 1;
     else {
	 cout << "wrong fitness value" << tempValue  << " with block size "
	      << BLOCKsize << " see binaryMCMC_trap.cpp:639\n";
	// cout << "wrong fitness value, see binaryMCMC_trap.cpp:613\n";
	 exit(1);
     }
#ifdef LINEAR_TRAP_FUNCTION
		value += tempValue;	
#ifdef DEBUG_
    cout << "("<<j+1<<","<<myValue<<")\t";
#endif
#else
		value += pow(j+1,ORDER_BLOCK)*tempValue;
#endif

  }	

  if(value == 0) value = AwayFromZero;
  return value;
}

int binaryMCMC::stop(){
//	if(value == size / BLOCKsize * bernoulli_table[0]) return 1;
	if(value == (size / BLOCKsize) * BLOCKsize) return 1;
	if(value == (size / BLOCKsize) * (BLOCKsize-1)) return 2;
	return 0;
}

void binaryMCMC::getMax(int* maxGenome){
    for(int j = 0; j < (size/(double)BLOCKsize); j++){
	for(unsigned long i = 0; i<BLOCKsize/2; i++)
	    maxGenome[i + j*BLOCKsize] = 0;
	for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++)
	    maxGenome[i + j*BLOCKsize] = 1;
    }
}

void binaryMCMC::getMax0(int* maxGenome){
    for(int j = 0; j < (size/(double)BLOCKsize); j++){
	for(unsigned long i = 0; i<BLOCKsize/2; i++)
	    maxGenome[i + j*BLOCKsize] = 1;
	for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++)
	    maxGenome[i + j*BLOCKsize] = 0;
    }
}

#ifdef RECOMB_MUT
void binaryMCMC::setMax0(){
    for(int j = 0; j < (size/(double)BLOCKsize); j++){
	for(unsigned long i = 0; i<BLOCKsize/2; i++){
	    genome[i + j*BLOCKsize] = 0;
	    cout << genome[i + j*BLOCKsize];
	}
	for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++){
	    genome[i + j*BLOCKsize] = 1;	
	    cout << genome[i + j*BLOCKsize];
	}
    }
    //cout << " = " << 
    fitness();
// << "\n";
}

void binaryMCMC::setMax1(){
    for(int j = 0; j < (size/(double)BLOCKsize); j++){
	for(unsigned long i = 0; i<BLOCKsize/2; i++){
	    genome[i + j*BLOCKsize] = 1;
	    cout << genome[i + j*BLOCKsize];
	}
	for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++){
	    genome[i + j*BLOCKsize] = 0;
	    cout << genome[i + j*BLOCKsize];
	}
    }
    //cout << " = " << 
    fitness();
    // << "\n";
}
#endif //RECOMB_MUT
#endif //ifndef BERNOULLI
#else
#ifdef QBP
double binaryMCMC::fitness(int* tempGenome, unsigned long mySize){
  int tempValue = 0; 

  for(unsigned long j = 0; j<mySize; j++)
      for(unsigned long i = 0; i<mySize; i++){ 
	  tempValue += matrixQBP[i][j] * tempGenome[i] * tempGenome[j];
      }
  
  if(tempValue == 0) tempValue = AwayFromZero;
  return tempValue;
}

double binaryMCMC::fitness(){
  int tempValue = 0; 
  value = 0;

  for(unsigned long j = 0; j<size; j++)
      for(unsigned long i = 0; i<size; i++){ 
	  tempValue += matrixQBP[i][j] * genome[i] * genome[j];
      }

  value = tempValue;
  if(value == 0) value = AwayFromZero;
  return value;
}

int binaryMCMC::stop(){
	if(value == size) return 1;
	if(value == size) return 2;
	return 0;
}

#else 
#ifdef MIXTURE_REALS
//if the reals values are in-between [0..2)
double binaryMCMC::fitness(){
#ifdef UNIFORM_DISTRIBUTION
  value = 1;
  return value;
#endif 

  double tempValue = 0; 
  //  value = 0;

////////
///// was the reals values in between [0..1)  
//  for(unsigned long i = 0; i<size; i++){ 
//      tempValue += 0.2*n(genome[i],0.1,0.2) + 0.4*n(genome[i],0.4,0.3) + 
//	  0.1*n(genome[i],0.6,0.1) + 0.3*n(genome[i],0.9,0.2); /// maxBLOCK;
///  }
///////
#ifndef BIVARIATE
#ifdef ROSENBROCK
  for(unsigned long i = 0; i<size - 1; i++){ 
    tempValue +=  (genome[i] * genome[i] - genome[i+1]) * (genome[i] * genome[i] - genome[i+1]) + 0.1* (genome[i] - 1) * (genome[i] - 1);
  }
#endif //rosenbrobk
#endif //BIVARIATE

#ifdef BIVARIATE
     tempValue = fitness_gauss(genome); /// maxBLOCK;
#endif 

  value = tempValue;
  //if(value < AwayFromZero) value = AwayFromZero;
  //cout << "big fitness" << tempValue << " for " << genome[0] << "\t" << genome[1] << "\n";
  return value;
}

double binaryMCMC::fitness(double* tempGenome, unsigned long mySize){
  double tempValue = 0; 

#ifdef UNIFORM_DISTRIBUTION
  return 1;
#endif 
///////
/// when reals valeas are in-between [0,1)  
//  for(unsigned long i = 0; i<mySize; i++){ 
//      tempValue += 0.2*n(tempGenome[i],0.1,0.2) + 0.4*n(tempGenome[i],0.4,0.3) + 
//	  0.1*n(tempGenome[i],0.6,0.1) + 0.3*n(tempGenome[i],0.9,0.2); /// maxBLOCK;
///  }
///////
///real value [0,2)
#ifndef BIVARIATE
#ifdef ROSENBROCK
  for(unsigned long i = 0; i<size - 1; i++){ 
    tempValue +=  (tempGenome[i] * tempGenome[i] - tempGenome[i+1]) * (tempGenome[i] * tempGenome[i] - tempGenome[i+1]) + 0.1*(tempGenome[i] - 1) * (tempGenome[i] - 1);
  }
#endif //rosenbrock
#endif //BIVARIATE


#ifdef BIVARIATE
      tempValue = fitness_gauss(tempGenome); /// maxBLOCK;
#endif //bivariate  

      //if(tempValue < AwayFromZero) tempValue = AwayFromZero; 
  //cout << "big fitness" << tempValue << " for " << tempGenome[0] << "\t" << tempGenome[1] << "\n";
  return tempValue;
}

double binaryMCMC::fitness(unsigned long* tempGenome, unsigned long mySize){
  double tempValue = 0; 
#ifdef UNIFORM_DISTRIBUTION
  return 1;
#endif 

  double recomb[mySize/BLOCKsize];
  convertReal(recomb,tempGenome,mySize/BLOCKsize);
  //for(unsigned long k = 0; k < mySize/BLOCKsize; k++){
  //    double tempReal = 0;
  //    for(int t = 0; t < BLOCKsize; t++)
//	  tempReal += tempGenome[k*BLOCKsize + t] * pow(2.0,(double)BLOCKsize - t - 1);
  //     recomb[k] = scale * tempReal / pow(2.0,(double)BLOCKsize);
  //}

#ifndef BIVARIATE
#ifdef ROSENBROCK
  for(unsigned long i = 0; i < size; i++)
    recomb[i] -= scale/2;

  for(unsigned long i = 0; i<size - 1; i++){ 
    tempValue +=  (recomb[i] * recomb[i] - recomb[i+1]) * (recomb[i] * recomb[i] - recomb[i+1]) + 0.1*(recomb[i] - 1) * (recomb[i] - 1);
  }
#endif //rosenbrock
#endif //BIVARIATE

#ifdef BIVARIATE
      tempValue =  fitness_gauss(recomb); // / maxBLOCK;
#endif // bivariate  

      //if(tempValue < AwayFromZero) tempValue = AwayFromZero;
  //cout << "big fitness" << tempValue << " for " << recomb[0] << "\t" << recomb[1] << "\n";
  
  return tempValue;
}

void binaryMCMC::getMax(double* maxGenome){
    for(unsigned long j = 0; j < size; j++){
	maxGenome[j] = scale* 0.5;
    }
}

void binaryMCMC::getMax0(double* maxGenome){
    for(unsigned long j = 0; j < size; j++){
	maxGenome[j] = 0;
    }
}

void binaryMCMC::getMax1(double* maxGenome){
    for(unsigned long j = 0; j < size; j++){
	maxGenome[j] = scale;
    }
}

#ifdef RECOMB_MUT
void binaryMCMC::setMax0(){
    //unsigned long tempGenome[size*BLOCKsize];
    //for(unsigned long j = 0; j < size*BLOCKsize; j++){
//	tempGenome[j] = 1;
    //  }
    //convertReal(genome,tempGenome,size);
    //for(unsigned long j = 0; j < size; j++){
//	genome[j] = 3.0*scale/4.0;
    //  }
    genome[0] = 0.6 * scale;
    genome[1] = 0.3 * scale;
}

void binaryMCMC::setMax1(){
    //for(unsigned long j = 0; j < size; j++){
//	genome[j] = scale/4.0;
    //  }
    genome[0] = 0.4 * scale;
    genome[1] = 0.7 * scale;
}

void binaryMCMC::setMax2(){
    //for(unsigned long j = 0; j < size; j++){
//	genome[j] = scale/4.0;
    //  }
    genome[0] = 0.4 * scale;
    genome[1] = 0.2 * scale;
}
#endif //recomb mut

double binaryMCMC::getMaxValue(){
  double maximVal = scale * 0.5;
    return  2.0*fitness_mixture(maximVal);
}

double binaryMCMC::n(double x, double mean, double sigma){
    return (1.0/(sigma*sqrt(2*pi))) *pow(e,-(x - mean)*(x-mean)/(2.0* sigma * sigma));
}

double binaryMCMC::fitness_mixture(double temp){
 return 2*n(temp,0.1,0.02) + 2.2*n(temp,0.3,0.03) + 2.5*n(temp,0.6,0.035) + 2.1*n(temp,0.9,0.02);
}

double binaryMCMC::fitness_bivariate(double* x,double s1, double s2, double m1, double m2, double r1){
  double z = (x[0] - m1)*(x[0] - m1)/(s1 * s1) - 2 * r1 * (x[0] - m1)*(x[1] - m2)/(s1 * s2) + (x[1] - m2)*(x[1] - m2)/(s2 * s2);
 
  double res = (1.0/(s1 * s2 * 2 * pi * sqrt(1 - r1 * r1))) * pow(e,- z /(2.0* (1 - r1 * r1)) );
  return res;
}

//sum over all mixtures
double binaryMCMC::fitness_multivariate(double* x){
  //cout << " one fitness evaluation \n";
  double value = 0;

  /* if(size == 2){
     for(int i = 0; i < mixture; i++){
  	 double rho = covarianceMixture[i].element(0,1) / sqrt(covarianceMixture[i].element(0,0) * covarianceMixture[i].element(1,1));
  	 double z = pow(x[0] - positionMixture[i].element(0,0),2.0)/covarianceMixture[i].element(0,0) + pow(x[1] - positionMixture[i].element(1,0),2.0)/covarianceMixture[i].element(1,1) - 2 * rho * (x[0] - positionMixture[i].element(0,0)) * (x[1] - positionMixture[i].element(1,0))/ sqrt(covarianceMixture[i].element(0,0) * covarianceMixture[i].element(1,1));
  	 value += exp(-z/(2 * (1 - rho * rho))) / (2 * pi * sqrt(covarianceMixture[i].element(0,0) * covarianceMixture[i].element(1,1) * ( 1 - rho * rho)));	 
     }
     } else {*/
  //cout << mixture << "\n";
      Matrix x_current(size,1);
      Matrix multiM(1,1); 
      
      for(int j = 0; j < size; j++)
	  x_current.element(j,0) = x[j];
      //cout << "current position " << x_current.t() << "\n";
      
      for(int i = 0; i < mixture; i++){	  
	  //cout << " position mixture " << positionMixture[i].t() << " difference with current " << (x_current - positionMixture[i]).t() 
	  //	 << "half matrix " << (x_current - positionMixture[i]).t() * covarianceMixtureInverse[i] << "\n";
	  
	  // double far = 0;
	  //for(int j = 0; j < size; j++)
	  //  if(!(x_current.element(j,0) - positionMixture[i].element(j,0) < 1 && x_current.element(j,0) - positionMixture[i].element(j,0) > -1))
	  //	far = 1;
	  // else {
	  ////cout << "not to far " << x_current.element(j,0) - positionMixture[i].element(j,0) << "\n"; 
	  //}
	  	  //if(far != 1){
	  multiM = (x_current - positionMixture[i]).t() * covarianceMixtureInverse[i] * (x_current - positionMixture[i]);
	  //cout << " for mixture " << i << " multiple with " << multiM << " determinant " << abs(determinantMixture[i]); 
	  value += pow(2.0 * pi,-((double)size)/2.0) * pow(abs(determinantMixture[i]),-1.0/2.0) * exp(-0.5 * multiM.element(0,0));
	  
	  //cout << " determinant^-1/2" << pow(abs(determinantMixture[i]),-1.0/2.0) << " size " << -((double)size)/2.0 
	  //<< " pi " << pow(2.0 * pi,-((double) size)/2.0) << " exp " 
	  //   << exp(-0.5 * multiM.element(0,0)) << "value " << value << "\n";
	  // } else {
	  //cout << " !!!!!!!!!!!!!!! too far away " << (x_current - positionMixture[i]).t();
	  //}
      }
      //}
  //cout << "fitness value " << value << "\n";
   value *= multiplication;
#ifdef AwayFromZero
   if(value < AwayFromZero) value = AwayFromZero;
#endif
  return value;
}

#ifdef BIVARIATE
double binaryMCMC::fitness_gauss(double* x){
#ifdef UNIFORM_DISTRIBUTION
  return 1;
#else

#ifndef MIXTURE_BIVARIATE
#ifndef TWO_PICKS
  //bivariate
  double s1 = sigma_1;
  double s2 = sigma_2;
  double m1 = (scale)/2.0;
  double m2 = (scale)/2.0;

  return fitness_bivariate(x,s1,s2,m1,m2,rho);

#else //two picks

  double m = scale/2.0;

  double s11 = sigma_11;
  double s12 = sigma_12;

  double s21 = sigma_21;
  double s22 = sigma_22;

  return coef1 * fitness_bivariate(x,s11,s12,m-0.25 * scale,m-0.25*scale,rho_1) + coef2 * fitness_bivariate(x,s21,s22,m + 0.25*scale,m + 0.25* scale,rho_2);

#endif //TWO_PICKS
#else 
#ifdef NINE_BIVARIATE
  //9 bivariates

  double m = scale/2.0;

  double s1 = sigma_1;
  double s2 = sigma_2;

  double s3 = sigma_3;
  double s4 = sigma_4;

  double s5 = sigma_5;
  double s6 = sigma_6;

  double* xNew = new double[size];
  xNew[0] = x[0] * 1;
  xNew[1] = x[1] * 1;

  return coef1 * fitness_bivariate(xNew,s1,s2,m-0.375*scale,m-0.375*scale,rho1) + coef1 * fitness_bivariate(xNew,s1,s2,m,m,rho1) + 
      coef1* fitness_bivariate(xNew,s1,s2,m+0.375*scale,m+0.375*scale,rho1)+ 
      coef3 * fitness_bivariate(xNew,s1,s2,m-0.375*scale,m,rho3) + coef3 * fitness_bivariate(xNew,s1,s2,m+0.375*scale,m,rho3)+
      coef3 * fitness_bivariate(xNew,s1,s2,m,m-0.375*scale,rho5) + coef3 * fitness_bivariate(xNew,s1,s2,m,m+0.375*scale,rho5)
      + coef5 * fitness_bivariate(xNew,s1,s2,m-0.375*scale,m+0.375*scale,rho1) + coef5 * fitness_bivariate(xNew,s1,s2,m+0.375*scale,m-0.375*scale,rho1);
#else 
  //cout << "I should enter here \n";
  return fitness_multivariate(x);
#endif //9_BIVARIATE
#endif //mixture bivariate
#endif // 
}
#endif //bivariate

unsigned long binaryMCMC::stop(){
    for(unsigned long i = 0; i < size; i++)
	if(genome[i] != scale/2) return 0;
    return 1;
}

#endif //MIXTURE_REALS   
#endif //QBP
#endif //BINOMIAL
#endif
#endif
#endif
#endif

#ifndef RANDOM_WALK
#ifndef MIXTURE_REALS
#ifdef METROPOLIS_ALGORITHM 
int binaryMCMC::acceptance(int* tempGenome, double tempValue){
#else
int binaryMCMC::acceptance(int* tempGenome, double tempValue, double tempProposalDistribution)
{
#endif //metropolis algorithm
#else //mixture reals
unsigned long binaryMCMC::acceptance(double* tempGenome, double tempValue){
#endif //mixture reals
    //cout << " tempvalue " << tempValue << " value " << value << "\n";
    //cout << "2 individuals tempIndiv = ";
    //for(int i = 0; i < size; i++)
//	cout << tempGenome[i] << "\t";
    //cout << "\n actualIndiv = ";
    //for(int i = 0; i < size; i++)
//	cout << genome[i] << "\t";
    //  cout << "\n";

  double prop = 0;

#ifdef Boltzmann
#ifdef Normal_Boltzmann
	if(tempValue == 0 && value == 0)
	    prop = 1;
	else {
	    if(tempValue == 0) tempValue = AwayFromZero;
	    if(value == 0) value = AwayFromZero;
	    //prop = exp((tempValue - value)/temperature);
	    prop = exp((double)(tempValue - value)/(double)DEFAULT_TEMPERATURE);
	}
#else 
	if(tempValue == 0 && value == 0)
	  prop = 1;
	else {
	  if(tempValue == 0) tempValue = AwayFromZero;
	  if(value == 0) value = AwayFromZero;
	  prop = exp((log(tempValue) - log(value))/(double)temperature);
	}
#endif //normal boltzman
#else //boltzmann
#ifndef MIXTURE_REALS
#ifdef FLOOR
	tempValue = floor(tempValue);
	value = floor(value);
#endif
#endif //mixture reals
	//if(tempValue == 0 && value == 0)
	//  prop = 1;
	//else{
#ifdef AwayFromZero
	  if(tempValue == 0) tempValue = AwayFromZero;
	  if(value == 0) value = AwayFromZero;	  
#endif
	  //atentie temporar
	  //prop = pow(tempValue/value,1.0/temperature);
	  //prop = tempValue/value;
	  /*  if(temperature != 0)
	      prop = pow(tempValue/value,1.0/(double)temperature);
	      else {*/
	      //cout << "temperature is 0\n"; 
	  prop = (double)tempValue/(double)value;
	  //cout << " probab acc " << prop << "\n";
	  //prop = ((double)tempValue/(double)value);
	      //}
	  //}
#endif 

#ifndef METROPOLIS_ALGORITHM
#ifndef MIXTURE_REALS
	ProposalDistribution();
	prop *= tempProposalDistribution/proposalDistribution;
#endif
#endif //METROPOLUS_HASTINGS_ALGORITHM

	//cout << "accetance with temp =" <<tempProposalDistribution << ", prop =" << proposalDistribution << ", frac =" << tempProposalDistribution/proposalDistribution << "\n";
#ifdef DEBUG_
	cout << "probability of acceptance " << prop << "\n";
#endif
	if(prop > 1) {
	  accepted = 1;
	  //cout << " accepted \n";
	  return 1;
	} else if(prop == 1){
	  return 1;
	}
	//UniformDistribution genrand_real1(0,1);
	double temp = genrand_real1();
	if(temp > prop) {
	    // cout << " rejected \n";
#ifdef DEBUG_		
	    cout << "Regection because random number = " << temp << "\n";
	      //     cout << ": generation = " << nrOfGenerations << "\n";
	      //cout << "proposed genome: |";
	      //for(int i=0; i<size; i++)
	      //cout << tempGenome[i] << " | ";	
	      //cout << " == "<<tempValue<<"\n";	
	      //for(int i=0; i<size; i++)
	      //cout << genome[i] << " | ";
	      //cout << " == "<<value<<"\n \n";
#endif
	  accepted = 0;
	  return 0;
	}

	//cout << " accepted \n";
	accepted = 1;
	return 1;
}

#endif //RANDOM_WALK

#ifndef MIXTURE_REALS
int* binaryMCMC::MutationGenome(int* newGenome){
    UniformDistribution genrand_real1(0,1);

  //  cout << "size:" << size << " mutation:" << mutation <<"\n";
  for(unsigned long i=0; i<size; i++)
    newGenome[i] = genome[i];
  
  //cout << "Next generation\n";
  /*  if(mutation == -1) {		
      //flip un bit din genome
    unsigned long first = genrand_int32() % size;
    newGenome[first] = (newGenome[first] + 1) % 2;
  } else {
  */    //mutatia standard - cu probabilitatea 1/l
    for(unsigned long i=0; i<size; i++){

//#ifdef INDEPENDENT_SAMPLER
	//algorithmul lui Laskey
//      if(2*mutation*newProposalDistribution[i][(newGenome[i] + 1)%2] > genrand_real2()){
//	newGenome[i] = (newGenome[i] + 1) % 2;
	//cout << "Mutation " << mutation << "am mutat " << i << "\n";
//      }
//#else
//#ifdef METROPOLIS_ALGORITHM 
      if(mutation > genrand_real1()){
	newGenome[i] = (newGenome[i] + 1) % 2;
	//cout << "Mutation " << mutation << "am mutat " << i << "\n";
	//  }
//#endif 
//#endif

      }
    }	

  return newGenome;		
}
#else
double binaryMCMC::bounds(double number){
  //cout << "enter bounds ";
    double temp = number;
#ifndef BIVARIATE
#ifdef CYCLE_SPACE
    if(temp < 0) temp = scale + temp;
    if(temp > scale) temp = temp -scale;
#else 
    cout << "torus";
    //# ifdef TORUS SPACE
    if(temp < -scale/2) {
      //           cout << "temp " << temp << "ceil(" << ceil((temp-scale/2)/scale) << ") ";
      temp = -scale/2.0 -  ((temp+scale/2)/scale - ceil((double)(temp-scale/2)/scale));
      //cout << "temp " << temp << "\n";
      // temp = - temp;
    }
    if(temp > scale/2) {
      //   cout << "temp " << temp << "ceil(" << ceil((temp-scale/2)/scale) << ") ";
      temp = scale/2.0 -  ((temp-scale/2)/scale - ceil((double)(temp-scale/2)/scale));
      //cout << "temp " << temp << "\n";
 //temp = 2* scale - temp;	    
    }
#endif //CYCLE SPACE
#else //bivariate
#ifdef CYCLE_SPACE
    if(temp < 0) temp = scale + temp;
    if(temp > scale) temp = temp -scale;
#else 
# ifdef TORUS_SPACE
    if(temp < -scale/2){
      ///      cout << "temp " << temp  << " temp+scale/2 " << temp+scale/2.0  << ", "<< (temp+scale/2.0)/scale<< " ceil (" << ceil((double)(temp+scale/2.0)/scale) << ") ";
      temp = -scale/2.0 -  ((temp+scale/2)/scale - (long)(temp+scale/2.0)/scale);
      //cout << "temp " << temp << "\n";
      //temp = scale + temp;
    }
    if(temp > scale/2) {
      //cout << "temp " << temp << " temp-scale/2 " << temp-scale/2.0  << ", "<< (temp-scale/2.0)/scale << "ceil(" << ceil((double)(temp-scale/2.0)/scale) 
      // << ") = " << ((temp-scale/2)/scale - (long)(temp-scale/2)/scale);
      temp = scale/2.0 -  ((temp-scale/2)/scale - (long)(temp-scale/2)/scale);
      //cout << " temp " << temp << "\n";
      //temp = temp - scale;
    }
#else //indef closed_space
    cout << "not inplemented here \n";
    exit(1);
#endif //torus space
#endif //CYCLE SPACE    
#endif //bivariate
    return temp;
}

double* binaryMCMC::MutationGenome(double* newGenome){
  //   cout << "size:" << size << " mutation:" << mutation_sigma << " mutation " << mutation << "\n";
  //cout << "mutation all genemoe ";
  //UniformDistribution genrand_real1(0,1);
  for(unsigned long i=0; i<size; i++){
    newGenome[i] = genome[i];
    //cout << genome[i] << "\t";
  }
  double tempValue; /* = fitness(genome,size);
  if(tempValue < MIN_FITNESS){
    cout << " Wrong individ transl to mutation simple" << genome[0] << "\t" << genome[1] << "\t" << tempValue << "\n";
    exit(1);
    } else {*/
    /*    cout << "mutate 1 ";
     for(unsigned long i=0;i<size;i++)
       cout << genome[i] << "\t";
       cout << tempValue << " into ";*/
  //}
    //cout << " mutation in mutationGenome binaryMCMC_trap.cpp:1762 ";
    //cout << " \n mutation scale" << scale << "\t";
    for(unsigned long i=0; i<size; i++){
      //if(mutation > genrand_real1()){
	    //NormalDistribution genrand_norm(newGenome[i],mutation_sigma);
	    //  cout << " at i=" << i << " old " << genome[i] << " new ";
	    double normal = gaussian();
#ifdef RANDOM_MUTATION
	   newGenome[i] = genrand_real1() * scale/2.0 * normal + newGenome[i]; 
#else
#ifdef OPTIMAL_MUTATION
	   newGenome[i] = normal * sqrt(covarianceMixture[0].element(i,i)) + newGenome[i];
#else
	   newGenome[i] = mutation_sigma * normal + newGenome[i];
#endif //optimal mutation
#endif //random mutation
	    //if genone is generated ustide the reange - it is a cycle
	    //cout << "(" << genome[i] << "->" << newGenome[i] << "," << normal << ",";
	    //  cout <<  newGenome[i] << ")\t";
	   /*#ifndef CLOSED_EQ_SPACE 
	    newGenome[i] = bounds(newGenome[i]);
	    #endif */
	    //}
	//else cout << "("<<  newGenome[i] << ")\t";
    }	
    //cout << "\n";
    tempValue = fitness(newGenome,size);
    //#ifdef CLOSED_EQ_SPACE 
    if(tempValue < MIN_FITNESS){
      /* cout << "failed ";
      for(unsigned long i=0;i<size;i++)
	cout << newGenome[i] << "\t";
      cout << tempValue << "\t return ";
      for(unsigned long i=0;i<size;i++)
	cout << genome[i] << "\t";
	cout << "\n";*/
      for(unsigned long i=0;i<size;i++)
	newGenome[i] = genome[i];
      accepted = 0;
    } else accepted = 1;
    //#endif 
    //ut << " succes "; 
    /*for(unsigned long i=0;i<size;i++)
      cout << newGenome[i] << "\t";
      cout << tempValue << " accept " << accepted << "\n";*/
    return newGenome;		
}

double* binaryMCMC::MutationOtherGenome(double* newGenome){
  //  cout << "size:" << size << " mutation:" << mutation_sigma << " mutation " << mutation << "\n";
    //UniformDistribution genrand_real1(0,1);
  //cout << " mutation old geneme ";
  double oldGenome[size];
  for(unsigned long i=0; i<size; i++)
    oldGenome[i] = newGenome[i];

  double tempValue; /* = fitness(oldGenome,size);
  if(tempValue < MIN_FITNESS){
    cout << " Wrong individ transl to mutation " << oldGenome[0] << " \t" << oldGenome[1] << "\t" << tempValue << "\n";
    exit(1);
  } else {
    cout << "mutate ";
     for(unsigned long i=0;i<size;i++)
       cout << oldGenome[i] << "\t";
     cout << tempValue << " into ";
     }*/
  //cout << newGenome[i] << "\t";
  //   }
  //  cout << " mutation in mutationOtherGenome binaryMCMC_trap.cpp:1798 ";*/
    //cout << " \n mutation scale" << scale << "\t";
    for(unsigned long i=0; i<size; i++){
      //if(mutation > genrand_real1()){
	    //NormalDistribution genrand_norm(newGenome[i],mutation_sigma);
	    //  cout << " at i=" << i << " old " << genome[i] << " new ";
	    double normal = gaussian();
#ifdef RANDOM_MUTATION
	   newGenome[i] = genrand_real1() * scale/2.0 * normal + newGenome[i]; 
#else
#ifdef OPTIMAL_MUTATION
	   newGenome[i] = normal * sqrt(covarianceMixture[0].element(i,i)) + newGenome[i];
#else
	    newGenome[i] = mutation_sigma * normal + newGenome[i];
#endif //optimal mutation
#endif //random mutation
	    //if genone is generated ustide the reange - it is a cycle
	    //cout << "(" << genome[i] << "->" << newGenome[i] << "," << normal << ",";
	    /*#ifndef CLOSED_EQ_SPACE 
	    newGenome[i] = bounds(newGenome[i]);
	    #endif*/
	    //cout <<  newGenome[i] << "\t";
	}
	//else cout << "("<<  newGenome[i] << ")\t";
    //}	
    //cout << "\n";
    tempValue = fitness(newGenome,size);
    //#ifdef CLOSED_EQ_SPACE 
    if(tempValue < MIN_FITNESS){
      /*cout << " failed ";
      for(unsigned long i=0;i<size;i++)
	cout << newGenome[i] << "\t";
	cout << tempValue << "\n";*/
      for(int i = 0; i < size; i++)
	newGenome[i] = oldGenome[i];
      accepted = 0;
    } else accepted = 1;
    
    //#endif 
    /*cout << " mutate ";
    for(unsigned long i=0;i<size;i++)
      cout << newGenome[i] << "\t";
      cout << tempValue << "\n";*/
    return newGenome;		
}

#endif //MIXTURE_REALS

#ifdef POPULATION_MCMC 
int* binaryMCMC::MutationGenome(int* newGenome, double** newProposalDistribution, unsigned long myRandJ){
  for(unsigned long i=0; i<size; i++)
    newGenome[i] = genome[i];

  randJ = myRandJ;
    if(2*mutation*newProposalDistribution[randJ][(newGenome[randJ] + 1)%2] > genrand_real1()){
	newGenome[randJ] = (newGenome[randJ] + 1) % 2;
	//cout << "Mutation " << mutation << "am mutat " << i << "\n";
      }
  return newGenome;		
}
#endif

#ifdef IsRecombination
#ifdef  BINARY_RECOMB
// this operation is implemented only for sizeCoupling = 2
unsigned long** binaryMCMC::RecombinationGenome(unsigned long** tableGenomes, unsigned long Recom_Type){
//    UniformDistribution genrand_real1(0,1);
#ifdef MIXTURE_REALS
    size = size * BLOCKsize;
    //    double* tempGenomeRecomb = new double[size/BLOCKsize];
    //   double* newGenome = new double[size/BLOCKsize];
#endif //mixture reals

  for(unsigned int i = 0; i < sizeCoupling; i++)
    for(unsigned int j = 0; j <size; j++)
      tableGenomes[i+sizeCoupling][j] = tableGenomes[i][j];

  //cout << "Next generation recombination " << Recom_Type<<"recombination number "<< recombination<<"\n";
  if(recombination == 0 || Recom_Type == 0) {
    //cout << "No recombination \n";
    return tableGenomes;
  } else if(recombination >= genrand_real1())
    {
      //cout << "Bigger than a random number\n";
      if(Recom_Type == 1) //one point cross over
	{
	  unsigned long randPos = genrand_int32() % size;
	  for(unsigned long i=0; i<size; i++)
	    if(randPos >= i){
	      tableGenomes[0+sizeCoupling][i] = tableGenomes[1][i];
	      tableGenomes[1+sizeCoupling][i] = tableGenomes[0][i];
	    }
	} else
	  if(Recom_Type == 2)
	    {
	      unsigned long randPos1 = genrand_int32() % size;
	      unsigned long randPos2 = genrand_int32() % size;
	      for(unsigned long i=0; i<size; i++)
		if((randPos1 > i && randPos2 < i) || (randPos1 < i && randPos2 > i)){
		  tableGenomes[0+sizeCoupling][i] = tableGenomes[1][i];
		  tableGenomes[1+sizeCoupling][i] = tableGenomes[0][i];
		}
	    }
	  else
	    if(Recom_Type == 3)
	      {
		//Uniform recombination
		//#ifdef DEBUG_
		//cout << "Uniform recombination switches with the probability"<< UNIF_RECOMB_PARAM <<"\n";
		//#endif
		for(unsigned long i=0; i<size; i++)
		  if(UNIF_RECOMB_PARAM > genrand_real1()){
		    tableGenomes[0+sizeCoupling][i] = tableGenomes[1][i]; 
		    tableGenomes[1+sizeCoupling][i] = tableGenomes[0][i]; 
		    //#ifdef DEBUG_
		    //cout << i << "\t";
		    //#endif
		  }
		//cout << "\n";
	      }
	    else
#ifndef METROPOLIS_ALGORITHM
		if(Recom_Type == 4){
		    //choose at random a parent to be replaced 
		    //cout << "Recombinare probabilistica \n";
		    //int parentReplaced = 0;
#ifndef TWO_CHILDREN
		  /*if(0.5 < genrand_real1()) {
			parentReplaced = 0;
		    
			for(unsigned long i=0; i<size; i++){
			    if(tableGenomes[0][i] != tableGenomes[1][i]){
				//if(newProposalDistribution[i][0] > genrand_real1())
				if(0.5 >= genrand_real1())
				    tableGenomes[0 + sizeCoupling][i] = 0;
				else tableGenomes[0 +sizeCoupling][i] = 1;
			    } else {
				if(mutation >= genrand_real1()) //mutation
				tableGenomes[0+sizeCoupling][i] = notAllele(tableGenomes[0][i]);
			    }
			    tableGenomes[1 + sizeCoupling][i] = tableGenomes[1][i];
			    //    cout << "from " << i << "  " << tableGenomes[0][i] << tableGenomes[1][i] << 
			    //" to " << tableGenomes[0+sizeCoupling][i] << tableGenomes[1+sizeCoupling][i] << 
			    //" " << parentReplaced << " 0 copied \n";
			}
			}else{*/
		  parentReplaced = 1;
		  
		  for(unsigned long i=0; i<size; i++){
		    if(tableGenomes[0][i] != tableGenomes[1][i]){
		      if(0.5 >= genrand_real1())
			tableGenomes[1 + sizeCoupling][i] = 0;
		      else tableGenomes[1 +sizeCoupling][i] = 1;
		    } else {
		      //if(mutation >= genrand_real1()) 
			tableGenomes[1+sizeCoupling][i] =tableGenomes[1][i];
		    }
		    tableGenomes[0 + sizeCoupling][i] = tableGenomes[0][i];
		    //cout << "from " << i << "  " << tableGenomes[0][i] << tableGenomes[1][i] << 
		    //" to " << tableGenomes[0+sizeCoupling][i] << tableGenomes[1+sizeCoupling][i] <<
		    //" " << parentReplaced <<" 1 copied \n";
		  }
#ifdef MIXTURE_REALS
		  double* tempGenome[size/BLOCKsize];
		  double* newGenome[size/BLOCKsize];
		  convertReal(tempGenome,tableGenomes[0],size/BLOCKsize);
		  MutationOtherGenome(tempGenome,newGenome);
		  convertBinary(newGenome,tableGenomes[0+sizeCoupling],size/BLOCKsize);			
#endif // mixture reals
			// MutationOtherGenome(tableGenomes[0], tableGenomes[0 + sizeCoupling]);
	    
		  //}
#else
		    for(unsigned long i=0; i<size; i++){
			if(tableGenomes[0][i] != tableGenomes[1][i]){
			    //if(newProposalDistribution[i][0] > genrand_real1())
			    if(0.5 >= genrand_real1())
				tableGenomes[0 + sizeCoupling][i] = 0;
			    else tableGenomes[0 +sizeCoupling][i] = 1;
			} else {
			    if(mutation >= genrand_real1()) 
				tableGenomes[0+sizeCoupling][i] = notAllele(tableGenomes[0][i]);
			}
		    }		    
		    for(unsigned long i=0; i<size; i++){
			if(tableGenomes[0][i] != tableGenomes[1][i]){
			    if(0.5 >= genrand_real1())
				tableGenomes[1 + sizeCoupling][i] = 0;
			    else tableGenomes[1 +sizeCoupling][i] = 1;
			} else {
			    if(mutation >= genrand_real1()) 
				tableGenomes[1+sizeCoupling][i] = notAllele(tableGenomes[1][i]);
			}
		    }
		    
#endif // TWO_CHILDREN 
		} else if(Recom_Type == 6){
			// protejeaza zonele comune
		  /*if(0.5 < genrand_real1()) {
			parentReplaced = 0;

			for(int i = 0; i < size; i++)
			if(tableGenomes[1][i] == tableGenomes[0][i]){
			    //do nothing
			    if(1.0/(double)size > genrand_real1()){
				tableGenomes[0+sizeCoupling][i] = (tableGenomes[0][i] + 1) % 2;
			    }
			    //if(mutation/(double)3.0 > genrand_real1()){
			     // tableGenomes[1+sizeCoupling][i] = (tableGenomes[0][i] + 1) % 2;
			     // }
			} else {
			    //mutate first parent with mutation rate and don not do anything else to the second parent
			    if(mutation > genrand_real1()){
				tableGenomes[0+sizeCoupling][i] = (tableGenomes[0][i] + 1) % 2;
			    }
			    //if(mutation > genrand_real1()){
			    // tableGenomes[1+sizeCoupling][i] = (tableGenomes[1][i] + 1) % 2;
			     // }
			}
    } else {*/
			parentReplaced = 1;

			for(int i = 0; i < size; i++)
			    if(tableGenomes[1][i] == tableGenomes[0][i]){
				//do nothing
			      //	if(1.0/(double)size > genrand_real1()){
			      //	tableGenomes[1+sizeCoupling][i] = (tableGenomes[1][i] + 1) % 2;
			      //}
			    /*if(mutation/(double)3.0 > genrand_real1()){
			      tableGenomes[1+sizeCoupling][i] = (tableGenomes[0][i] + 1) % 2;
			      }*/
			} else {
			    //mutate first parent with mutation rate and don not do anything else to the second parent
			    if(mutation > genrand_real1()){
				tableGenomes[1+sizeCoupling][i] = (tableGenomes[1][i] + 1) % 2;
			    }
			    /*if(mutation > genrand_real1()){
			      tableGenomes[1+sizeCoupling][i] = (tableGenomes[1][i] + 1) % 2;
			      }*/
			}
#ifdef MIXTURE_REALS
			double* tempGenome = new double[size/BLOCKsize];
			double* newGenome = new double[size/BLOCKsize];
			convertReal(tempGenome,tableGenomes[0],size/BLOCKsize);
			MutationOtherGenome(tempGenome,newGenome);
			convertBinary(newGenome,tableGenomes[0+sizeCoupling],size/BLOCKsize);			
#endif // mixture reals
			//}
		}
		else
#endif //METROPOLIS_ALGORITHM	  
		    if(Recom_Type == 12)	
		    {
			//	cout << "recombinare cu parametri " << PROB_nPOINT << "\n";
			unsigned long randPos[PROB_nPOINT+1];
			int i = 0;
			while(i < PROB_nPOINT){
			     randPos[i] = genrand_int32() % size;
			     int gasit = 0;
			     for(int j = 0; j < i; j++)
				 if(randPos[i] == randPos[j])
				     gasit = 1;
			     if(gasit == 0)
				 i++;
			}
			randPos[PROB_nPOINT] = size; // punct de schimbare
			
			//sorteaza noul vector - bbuble sort?
			for(i = 0; i < PROB_nPOINT; i++)
			    for(int j = 0; j < PROB_nPOINT; j++)
				if((randPos[i] > randPos[j] && i < j) || (randPos[i] < randPos[j] && i > j)){
				    unsigned long temp = randPos[i];
				    randPos[i] = randPos[j];
				    randPos[j] = temp;
				}
			//cout << "swap prob ";
			/*for(i = 0; i < PROB_nPOINT+1; i++){
			    cout << "(" << i << "," << randPos[i] << ")\t";
			}
			cout << "\n swap at position ";*/

		        // pozitia in vector
			unsigned long j = 0; i = 0;
			while(i < size && j <= PROB_nPOINT){
			    if((j == 0 && randPos[j] > i)||(j %2 == 0 && randPos[j] > i)){
				tableGenomes[0+sizeCoupling][i] = tableGenomes[1][i]; 
				tableGenomes[1+sizeCoupling][i] = tableGenomes[0][i]; 
				//cout << i << "\t";
			    }
			    i++;
			    if(randPos[j] <= i){
				j++;
				//cout << "increase j " << j << "\n";
			    }
			}
			//cout << " I " << i << " j " << j << "\n";
		    }
		    else if(Recom_Type == 6){
			// protejeaza zonele comune
			//for(unsigned long i = 0; i < size; i++)
			//if(tableGenomes[1][i] == tableGenomes[0][i]){
			    //do nothing
			    /*if(mutation/(double)3.0 > genrand_real1()){
				tableGenomes[0+sizeCoupling][i] = (tableGenomes[0][i] + 1) % 2;
				}*/
			    /*if(mutation/(double)3.0 > genrand_real1()){
				tableGenomes[1+sizeCoupling][i] = (tableGenomes[0][i] + 1) % 2;
				}*/
		      //} else {
			    //mutate first parent with mutation rate and don not do anything else to the second parent
		      //    if(mutation > genrand_real1()){
		      //	tableGenomes[0+sizeCoupling][i] = (tableGenomes[0][i] + 1) % 2;
		      //	    }
			    /*if(mutation > genrand_real1()){
				tableGenomes[1+sizeCoupling][i] = (tableGenomes[1][i] + 1) % 2;
				}*/
		      //	}
		      //parentReplaced = 1;

			for(int i = 0; i < size; i++)
			    if(tableGenomes[1][i] == tableGenomes[0][i]){
				//do nothing
			    } else {
			      //mutate first parent with mutation rate and don not do anything else to the second parent
			      if(mutation > genrand_real1()){
				tableGenomes[1+sizeCoupling][i] = (tableGenomes[1][i] + 1) % 2;
			      }
			      /*if(mutation > genrand_real1()){
				tableGenomes[1+sizeCoupling][i] = (tableGenomes[1][i] + 1) % 2;
				}*/
			    }
			
			//#ifdef MIXTURE_REALS
			//convertReal(tempGenomeRecomb,tableGenomes[0],size/BLOCKsize);
			//MutationOtherGenome(tempGenomeRecomb);
			//convertBinary(tempGenomeRecomb,tableGenomes[0+sizeCoupling],size/BLOCKsize);
			//#endif // mixture reals
			//}

		    }
		    else 
			cout << "Tip inexistent recomb " << Recom_Type << "\n";
    }
	/*cout << "Initial \n";
	for(unsigned long i = 0; i < sizeCoupling; i++){
	  for(unsigned int j = 0; j <size; j++)
	    cout << tableGenomes[i][j] << "|";
	  cout << "\n";
	}
	cout << "Proposed \n";
	for(unsigned long i = 0; i < sizeCoupling; i++){
	  for(unsigned int j = 0; j <size; j++)
	    cout << tableGenomes[i+sizeCoupling][j] << "|";
	  cout << "\n";
	  }*/
	
#ifdef MIXTURE_REALS
  size/= BLOCKsize;
  //if(Recom_Type == 6){  
  //delete[] tempGenomeRecomb;
  //}
#endif
	return tableGenomes;		
}

#ifndef METROPOLIS_ALGORITHM
#ifdef POPULATION_MCMC
unsigned long** binaryMCMC::RecombinationGenome(unsigned long** tableGenomes, unsigned long Recom_Type, double** newProposalDistribution){
  for(unsigned unsigned long i = 0; i < sizeCoupling; i++)
    for(unsigned unsigned long j = 0; j <size; j++)
      tableGenomes[i+sizeCoupling][j] = tableGenomes[i][j];
  
  //cout << "Next generation recombination " << Recom_Type<<"recombination number "<< recombination<<"\n";
  if(recombination == 0 || Recom_Type == 0) {		
    //cout << "No recombination \n";
    return tableGenomes;
  } else if(recombination >= genrand_real1())
    {
      //cout << "Bigger than a random number\n";
	//recombination si mutation in acelasi operator cu generarea unui singur copil
	if(Recom_Type == 4){
	    //choose at random a parent to be replaced 

	    unsigned long parentReplaced = 0;
	    if(0.5 < genrand_real1()) parentReplaced = 1;

	    for(unsigned long i = 0; i < size; i++)
		if(tableGenomes[0][i] != tableGenomes[1][i]){
		    //if(newProposalDistribution[i][0] > genrand_real1())
		    if(0.5 > genrand_real1()) 
			tableGenomes[(parentReplaced + 1)%2 + sizeCoupling][i] = 0;
		    else tableGenomes[(parentReplaced + 1)%2 +sizeCoupling][i] = 1;
		    
		    tableGenomes[parentReplaced+sizeCoupling][i] = tableGenomes[parentReplaced][i];
		    
		} else {
		    if(mutation > genrand_real1()) 
			tableGenomes[0+sizeCoupling][i] = (tableGenomes[0+sizeCoupling][i] +1)%2;
		    if(mutation > genrand_real1()) 
			tableGenomes[1+sizeCoupling][i] = (tableGenomes[1+sizeCoupling][i] +1)%2;
		}
	    //cout << "Am terminat selection \n";
	} else 
	    if(Recom_Type == 5){
		for(unsigned long i = 0; i < size; i++)
		    if(tableGenomes[0][i] != tableGenomes[1][i]){
			if(newProposalDistribution[i][0] > genrand_real1())
			    //if(0.5 > genrand_real1()) 
			    tableGenomes[0+sizeCoupling][i] = 0;
			else tableGenomes[0+sizeCoupling][i] = 1;
			
			if(newProposalDistribution[i][0] > genrand_real1()) 
			    //if(0.5 > genrand_real1()) 
			    tableGenomes[1+sizeCoupling][i] = 0;
			else tableGenomes[1+sizeCoupling][i] = 1;
			
			/*double probZero = 1;
			  if(tableGenomes[0][i] == 0) probZero++;
			  if(tableGenomes[1][i] == 0) probZero++;
			  probZero /=(double)(2+2);
			  if(probZero >= genrand_real1()) tableGenomes[0 + sizeCoupling][i] = 0;
			  else tableGenomes[0 + sizeCoupling][i] = 1;
			  if(probZero >= genrand_real1()) tableGenomes[1 + sizeCoupling][i] = 0;
			  else tableGenomes[1 + sizeCoupling][i] = 1;*/
		    } else {
			if(mutation/2.0 > genrand_real1())
			    tableGenomes[0+sizeCoupling][i] = (tableGenomes[0+sizeCoupling][i] +1)%2;
			if(mutation/2.0 > genrand_real1())
			    tableGenomes[1+sizeCoupling][i] = (tableGenomes[1+sizeCoupling][i] +1)%2;
		    }
	    }
	    else {
		cout << "Error undefinit recombination number =" << Recom_Type << "\n";
	    }
    }
  /*cout << "Initial \n";
    for(unsigned long i = 0; i < sizeCoupling; i++){
    for(unsigned int j = 0; j <size; j++)
    cout << tableGenomes[i][j] << "|";
    cout << "\n";
    }
    cout << "Proposed \n";
    for(unsigned long i = 0; i < sizeCoupling; i++){
    for(unsigned int j = 0; j <size; j++)
    cout << tableGenomes[i+sizeCoupling][j] << "|";
    cout << "\n";
    }*/

  return tableGenomes;
}
#endif //POPULATION_MCMC
#endif //METROPOLIS_ALGORITHM

#else //real recombination
double** binaryMCMC::RecombinationGenome(double** tableGenomes, unsigned long Recom_Type){
  //cout << "proposed \n";
  for(unsigned int i = 0; i < sizeCoupling; i++){
    //   cout << "(";
    for(unsigned int j = 0; j <size; j++){
      tableGenomes[i+sizeCoupling][j] = tableGenomes[i][j];
      //cout << tableGenomes[i][j] << "\t";
    }
    //cout << ")\t";
    /*if(fitness(tableGenomes[i],size) < MIN_FITNESS){
      cout << " wrong individual transmeded to recombination " << fitness(tableGenomes[i],size) << " \n";
      exit(1);
      }*/
  }
  //cout << "\n";
  /*cout << "Next generation recombination " << Recom_Type <<" with sizeCoupling " << sizeCoupling << "\n";
  cout << " given individuals \n";
  for(unsigned int i = 0; i < sizeCoupling; i++){
    for(unsigned int j = 0; j <size; j++)
      cout << " " << tableGenomes[i][j];
    cout << "\n";
  }
  cout << "\n";*/

  if(Recom_Type == PARENT_CENTRIC) {		
    // aply parent centric always on the first position
    //since the parents are randomly picked, this does not have any effect

    //1.0 sample on a direction according with the position of the center of all individuals
    //1.1compute the center
    double g[size];
    double euclidian[sizeCoupling];
    //#ifdef PCTX
#ifdef wPCTX
    double fit[sizeCoupling];
    for(int i = 0; i < sizeCoupling; i++)
      fit[i] = fitness(tableGenomes[i+sizeCoupling],size);
#endif //wpctx 
    //   cout << " compute the center of the indiv";
    for(int i = 0; i < size; i++){     
      g[i] = 0;
#ifndef wPCTX
      for(unsigned long j = 0; j < sizeCoupling; j++){
	g[i]+= tableGenomes[j+sizeCoupling][i];    
      }
      g[i] /= (double) sizeCoupling;
#else
      //the center of the mass contains the fitness function
      double jos = 0;
      for(unsigned long j = 0; j < sizeCoupling; j++){
	g[i]+= tableGenomes[j + sizeCoupling][i] * fit[j];
	jos += fit[j];
      }
      //if(g[i] != 0 && jos != 0) 
      g[i] /= jos;
	  //else if(g[i] == 0 && jos == 0) g[i] = 0;
	  //else {
	  //cout << " Erroare in computation in binaryMCMC_trap.cpp::2345 \n";
	  //exit(1);
	  //}
#endif 
      //cout << " " << g[i];
    }
    //cout << "\n";
    //1.2 compute the distances to the center of the parent for recombination
    /*double d[size];
    //cout << " compute the distance to the center of the first individual ";
    for(int i = 0; i < size; i++){     
      d[i] = tableGenomes[0 + sizeCoupling][i] - g[i];
      //cout << " " << d[i];
    }
    //cout << "\n";
    */
    //1.3 compute the eucidian distance from the center to the parents
    //cout << " compute the euclidian distance to the center of the two individuals";
    for(int j = 0; j < sizeCoupling; j++){
      euclidian[j] = 0;
      for(int i =0; i < size; i++)
	euclidian[j] += pow(tableGenomes[j + sizeCoupling][i] - g[i],2.0);
      euclidian[j] = sqrt(euclidian[j]);
      //cout << " j= " << j << "\t" << euclidian[j];
    }
    //#endif //PCTX

    //cout << "\n";
    //1.4 only if it is respectful
    /*int same = 0;
    for(int i = 0; i < size; i++)
      for(int j = 0; j < size; j++)
	if(i != j && euclidian[i] != euclidian[j])
	  same ==1;
    if(same == 0){ 
      cout << " coupling = " << sizeCoupling << "all point in the same place \n";
      return tableGenomes;
      }*/

    //2.0 copy in an intermediate buffer
    double** tempGenome = new double*[sizeCoupling];
    for(int j = 0; j < sizeCoupling; j++){
      tempGenome[j] = new double[size];
      for(int i = 0; i < size; i++)
	tempGenome[j][i]= tableGenomes[j + sizeCoupling][i];
    }

    //2.1 operatiile de translation, rotatie si scalare

#ifdef TRANSL_ROT

    //0. assymmetrical operator
#ifdef wPCTX
    //if(mix_mixture == 1){
#ifdef TRANSLATION
      translation(tempGenome,g,euclidian);
      //mix_mixture = 1;
      //#else //rotation
      //rotation(tempGenome,g,euclidian);
      //mix_mixture = 0;
#endif //translation    
      /*    } else {
#ifdef ROTATION
      rotation(tempGenome,g,euclidian);
      mix_mixture = 0;
#else //rotation
      translation(tempGenome,g,euclidian);
      mix_mixture = 1;
#endif //translation    
}*/
#else //wPCTX
    //#ifdef wPCTX
    //#ifdef TRANSLATION
    // translation(tempGenome,g,euclidian);
    //#endif //translation    
    //#endif //wPCTX

     double randValue = genrand_real1(); 

#ifdef MIXTURE_TRANSL_SCALING
      if(randValue < MIXTURE_TRANSL_SCALING){
#ifdef TRANSLATION
	translation(tempGenome,g,euclidian);
#endif //translation
      } else if(randValue < MIXTURE_TRANSL_SCALING + 0.5)  {
#ifdef ROTATION
	rotation(tempGenome,g,euclidian);
#endif  // rotation
     } else {
#ifdef SCALING
	scaling_dim(tempGenome,g,euclidian);
#endif //SCALING
}  
#endif

    //1. tested algorithm; without mutation
      //cout << " before mixture \n";
#ifdef MIXTURE_TRANSL_ROT_SCALE
      //cout << "transl rot scale mixture\n";
      if(1.0/3.0 > randValue){
	//cout << "translation \n";
#ifdef TRANSLATION
	translation(tempGenome,g,euclidian);
#endif //translation
	
      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "trans(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}
          
      } else if(2.0 / 3.0 > randValue){
	//cout << "rotation \n";
#ifdef ROTATION
	rotation(tempGenome,g,euclidian);
#endif  // rotation
	
	// for(int j = 0; j < sizeCoupling; j++)
	//for(int i = 0; i < size; i++){
	//if(abs(tempGenome[j][i]) > scale/2) cout << "rotation(" << i << "," <<  tempGenome[j][i] << ")\n";
	//}

      } else {
	//cout << "scaling \n";
#ifdef SCALING
	scaling_dim(tempGenome,g,euclidian);
#endif //SCALING
      }
#endif // mixture without mutation

    //1.1. tested algorithm; without mutation
#ifdef CYCLE
      //cout << "cycle 2 posibilities \n";
      if(0.5 > randValue){
#ifdef TRANSLATION
      translation(tempGenome,g,euclidian);
#endif //translation

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "trans(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}
    
#ifdef SCALING
      scaling_dim(tempGenome,g,euclidian);
#endif //SCALING

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "scale(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef ROTATION
      rotation(tempGenome,g,euclidian);
#endif  // rotation
      
      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "rotation(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}
      }else{
#ifdef ROTATION
      rotation(tempGenome,g,euclidian);
#endif  // rotation
      
      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "rotation(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}
    
#ifdef SCALING
      scaling_dim(tempGenome,g,euclidian);
#endif //SCALING

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "scale(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef TRANSLATION
      translation(tempGenome,g,euclidian);
#endif //translation

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "trans(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}
      } 
#endif 

#ifdef CYCLE_TRANSL_ROT_SCALE
      //cout << " cycle transl rot scale 6 posibilities \n";
    if(1.0/6.0 > randValue){

#ifdef TRANSLATION
      //cout << "here \n";
      translation(tempGenome,g,euclidian);
#endif //translation

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "trans(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}
    
#ifdef SCALING
      scaling_dim(tempGenome,g,euclidian);
#endif //SCALING

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "scale(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef ROTATION
      rotation(tempGenome,g,euclidian);
#endif  // rotation
      
      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "rotation(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}
      
      } else if(2.0 / 6.0 > randValue){
#ifdef ROTATION
	rotation(tempGenome,g,euclidian);
#endif  // rotation

	// for(int j = 0; j < sizeCoupling; j++)
	//for(int i = 0; i < size; i++){
	//if(abs(tempGenome[j][i]) > scale/2) cout << "rotation(" << i << "," <<  tempGenome[j][i] << ")\n";
	//}

#ifdef SCALING
	scaling_dim(tempGenome,g,euclidian);
#endif //SCALING
	
 	//for(int j = 0; j < sizeCoupling; j++)
	//for(int i = 0; i < size; i++){
	//if(abs(tempGenome[j][i]) > scale/2) cout << "scale(" << i << "," <<  tempGenome[j][i] << ")\n";
	//}

#ifdef TRANSLATION
	translation(tempGenome,g,euclidian);
#endif //translation
	
	//for(int j = 0; j < sizeCoupling; j++)
	//for(int i = 0; i < size; i++){
	//if(abs(tempGenome[j][i]) > scale/2) cout << "transl(" << i << "," <<  tempGenome[j][i] << ")\n";
	//}

      } else if(3.0 / 6.0 > randValue){
    
#ifdef SCALING
      scaling_dim(tempGenome,g,euclidian);
#endif //SCALING

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "scale(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef ROTATION
      rotation(tempGenome,g,euclidian);
#endif  // rotation

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "trans(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef TRANSLATION
      translation(tempGenome,g,euclidian);
#endif //translation

      } else if(4.0 / 6.0 > randValue){
    
#ifdef TRANSLATION
      translation(tempGenome,g,euclidian);
#endif //translation

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "scale(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef ROTATION
      rotation(tempGenome,g,euclidian);
#endif  // rotation

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "trans(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef SCALING
      scaling_dim(tempGenome,g,euclidian);
#endif //SCALING

      } else if(5.0 / 6.0 > randValue){
    
#ifdef ROTATION
      rotation(tempGenome,g,euclidian);
#endif  // rotation

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "scale(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef TRANSLATION
      translation(tempGenome,g,euclidian);
#endif //translation

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "trans(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef SCALING
      scaling_dim(tempGenome,g,euclidian);
#endif //SCALING

      } else {
    
#ifdef SCALING
      scaling_dim(tempGenome,g,euclidian);
#endif //SCALING

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "trans(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef TRANSLATION
      translation(tempGenome,g,euclidian);
#endif //translation

      //for(int j = 0; j < sizeCoupling; j++)
      //for(int i = 0; i < size; i++){
      //if(abs(tempGenome[j][i]) > scale/2) cout << "scale(" << i << "," <<  tempGenome[j][i] << ")\n";
      //}

#ifdef ROTATION
      rotation(tempGenome,g,euclidian);
#endif  // rotation

      }
#endif //cycle transl rotation scaling

    //2. tested algorithm: without mutation
#ifdef MIXTURE_TRANSL_DIFF_ROT
    //cout << " before diff ";
    if(MIXTURE_TRANSL_DIFF_ROT < genrand_real1()){
#ifdef ROTATION
      rotation(tempGenome,g,euclidian);
#endif  // rotation
    } else {
#ifdef TRANSLATION
      translation(tempGenome,g,euclidian);
#endif //translation
    }
#endif //mixture_transl_diff_rot
#endif //wPCTX
    //3. snooker recomb, symmetrical, one parent generated

    //4. sum
#ifdef SUM_TRANSL_ROT_SCALE
#ifdef TRANSLATION
      translation(tempGenome,g,euclidian);
#endif //translation

    double** tempGenomeR = new double*[sizeCoupling];
    for(int j = 0; j < sizeCoupling; j++){
      tempGenomeR[j] = new double[size];
      for(int i = 0; i < size; i++)
	tempGenomeR[j][i]= tableGenomes[j + sizeCoupling][i];
    }

#ifdef ROTATION
    rotation(tempGenomeR,g,euclidian);
#endif //rotation

    for(int j = 0; j < sizeCoupling; j++){
      for(int i = 0; i < size; i++)
	tempGenome[j][i] += tempGenomeR[j][i];
    }

    for(int j = 0; j < sizeCoupling; j++){
      for(int i = 0; i < size; i++)
	tempGenomeR[j][i]= tableGenomes[j + sizeCoupling][i];
    }

#ifdef SCALING
    scaling_dim(tempGenomeR,g,euclidian);
#endif //scaling

    for(int j = 0; j < sizeCoupling; j++){
      for(int i = 0; i < size; i++)
	tempGenome[j][i] += tempGenomeR[j][i];
    }

    delete[] tempGenomeR;
#endif //sum transl rot scale

#endif //rot scale

    //actualizate the candidate individuals
#ifdef PCTX
    for(int j = 0; j < sizeCoupling; j++){
#else
    for(int j = 0; j < 1; j++){
#endif
      for(int i = 0; i < size; i++)
	tableGenomes[j + sizeCoupling][i] = tempGenome[j][i];
    }
    
    delete[] tempGenome;
  } 
  //cout << " scale " << scale << "\n";
#ifdef DEBUG
      cout << "resulting individuals \t";
  for(int j = 0; j < sizeCoupling; j++){
    cout << j << " =( ";
    for(int i = 0; i < size; i++)
      cout << "(" << tableGenomes[j][i] << "," << tableGenomes[j+sizeCoupling][i] <<") ";
    cout << ")\t";
  }
  cout << "\n";
#endif 
  //fix the values in the bounds
  /*#ifndef CLOSED_EQ_SPACE 
  bool stop = false;
  for(int j = 0; j < sizeCoupling; j++)
    for(int i = 0; i < size; i++){
      //if(stop == false)
      //cout << "bounds (" << j << "," << i << "," << tableGenomes[j + sizeCoupling][i] << ")\n";
      tableGenomes[j + sizeCoupling][i] = bounds(tableGenomes[j + sizeCoupling][i]);
      if(abs(tableGenomes[j + sizeCoupling][i]) > scale/2) {
	cout << "w(" << i << "," <<  tableGenomes[j + sizeCoupling][i] << ")\t";
	stop = true;
      }
    }
  if(stop == true) exit(0);
  #else*/
#ifdef SIMPLEX
    for(int j =0; j < sizeCoupling; j++){
      if(fitness(tableGenomes[j + sizeCoupling],size) < MIN_FITNESS){
	for(int t = 0; t < sizeCoupling; t++)
	  for(int t2 = 0; t2 < size; t2++)
	    tableGenomes[t + sizeCoupling][t2] = tableGenomes[t][t2];      
	accepted = 0;
	//cout << " not accept ANY \n";
	return tableGenomes;
      } 
    }
    //cout << "accept \n";
    accepted = 1;
#else //simplex
#ifdef PCTX
#ifdef SINGLE_INDIV_ROT
    for(int j =0; j < 1; j++){
       if(fitness(tableGenomes[j + sizeCoupling],size) < MIN_FITNESS){
	  //for(int t = 0; t < sizeCoupling; t++)
	  for(int t2 = 0; t2 < size; t2++)
	      tableGenomes[0 + sizeCoupling][t2] = tableGenomes[0][t2];      
	  accepted = 0;
	  //cout << " not accept ANY \n";
	  return tableGenomes;
      } 
    }
    accepted = 1;
#else //single indiv rot
#ifndef wPCTX
  for(int j =0; j < sizeCoupling; j++){
#else
  for(int j =0; j < 1; j++){
#endif
    if(fitness(tableGenomes[j + sizeCoupling],size) < MIN_FITNESS){
      for(int t = 0; t < sizeCoupling; t++)
	for(int t2 = 0; t2 < size; t2++)
	  tableGenomes[t + sizeCoupling][t2] = tableGenomes[t][t2];      
      accepted = 0;
      //cout << " not accept ANY \n";
      return tableGenomes;
    } 
  }
  accepted = 1;
  // cout << "accept all \n";
#endif //single indiv rot
#else //ifndef PCTX
#ifdef HALF_PCTX
  if(mix_rotation != 0){
    accepted = 0;
    for(int j =0; j < sizeCoupling; j++){
      if(fitness(tableGenomes[j + sizeCoupling],size) < MIN_FITNESS){
	for(int t2 = 0; t2 < size; t2++)
	    tableGenomes[j + sizeCoupling][t2] = tableGenomes[j][t2];      
	//accepted+= 0;
	//cout << " not accept ANY \n";
      } else accepted++;
    }
    return tableGenomes;
  } else {
  accepted = 0;
  double tempValueFit = fitness(tableGenomes[0 + sizeCoupling],size);
    if(tempValueFit < MIN_FITNESS){
      for(int t2 = 0; t2 < size; t2++)
	tableGenomes[0 + sizeCoupling][t2] = tableGenomes[0][t2];      
      accepted++;
      //cout << " NOT accept " << tableGenomes[0 + sizeCoupling][0] << " , " << tableGenomes[0 + sizeCoupling][1] << " , " << tempValueFit << " \n";
      //return tableGenomes;
    } //else cout << " ACCEPT " << tableGenomes[0 + sizeCoupling][0] << " , " << tableGenomes[0 + sizeCoupling][1] << " , " << tempValueFit << " \n";
    //}
    
    for(int j =1; j < sizeCoupling; j++){
      if(fitness(tableGenomes[j + sizeCoupling],size) < MIN_FITNESS){
	for(int t = 1; t < sizeCoupling; t++)
	  for(int t2 = 0; t2 < size; t2++)
	    tableGenomes[t + sizeCoupling][t2] = tableGenomes[t][t2];      
	//accepted+= 0;
	//cout << " not accept ANY \n";
	return tableGenomes;
      }
    }
 
    accepted += sizeCoupling -1;
  //cout << " accept all \n";
  }
#else
  //for(int j =0; j < 1; j++){
  double tempValueFit = fitness(tableGenomes[0 + sizeCoupling],size);
    if(tempValueFit < MIN_FITNESS){
      for(int t = 0; t < 1; t++)
	for(int t2 = 0; t2 < size; t2++)
	  tableGenomes[t + sizeCoupling][t2] = tableGenomes[t][t2];      
      accepted = 0;
      //cout << " NOT accept " << tableGenomes[0 + sizeCoupling][0] << " , " << tableGenomes[0 + sizeCoupling][1] << " , " << tempValueFit << " \n";
      return tableGenomes;
    } //else cout << " ACCEPT " << tableGenomes[0 + sizeCoupling][0] << " , " << tableGenomes[0 + sizeCoupling][1] << " , " << tempValueFit << " \n";
    //}
  accepted = 1;
  //cout << " accept all \n";
#endif //half_PCTX
#endif //pctx
#endif //simplex

  //#endif 
#ifdef DEBUG  
  cout << "resulting individuals after bounds";
  for(int j = 0; j < sizeCoupling; j++){
    cout << j << " =( ";
    for(int i = 0; i < size; i++)
      cout << tableGenomes[j+sizeCoupling][i] <<",";
    cout << ")\t";
  }
  cout << "\n\n\n";
#endif 

  return tableGenomes;
}

#ifdef TWO_PARENTS_SIMPLE

#else //two parents
#ifdef TRANSLATION
void binaryMCMC::translation(double** tempGenome, double* g, double* euclidian){
#ifdef SINGLE_INDIV_ROT
  mix_rotation = 0; 
  //translationDiff(tempGenome, g, euclidian);
  //return;
#endif 

  //cout << "transl\t";
#ifdef SIMPLEX
  //double fac;
  double randValue = genrand_real1(); 

  //if(randValue < 0.5){
    fac = -1.0;
    translationSimplex(tempGenome, g, euclidian);
    //} else{
    //positionChanged = 0;
    //translationDiff(tempGenome, g, euclidian);
    //mix_rotation = -2;
    //}
    /*else if(randValue < 0.75)
    fac = 2.0;
  else 
  fac = 0.5;*/
    
#else
#ifdef SYMMETRICAL_SNOOKER
    translationSnooker(tempGenome, g);
    cout << "should not enter here \n";
#else
#ifdef PCTX
/*   translationDiff(tempGenome, g, euclidian);

   //cycles on individuals between PCTX and Diff
   translationPCTX((tempGenome + 1), (g+1), (euclidian + 1), sizeCoupling - 1);
*/
    translationPCTX(tempGenome, g, euclidian, sizeCoupling);
#ifdef SUM_TRANSL 
   for(int i = 0; i < size; i++)
     tempGenome[0][i] = tempGenome[0][i] + gaussian() * 0.001;
#endif 
   //cout << "shound not enter here \n";
#else
   //   cout << "should enter here\n";
   // difference recombination on three individuals
   //first individual is translated
   translationDiff(tempGenome, g, euclidian);

#ifdef HALF_PCTX
   //cycles on individuals between PCTX and Diff
   translationPCTX((tempGenome + 1), (g+1), (euclidian + 1), sizeCoupling - 1);
#endif //half PCTX
#endif //PCTX
#endif //symmetrical snooker
#endif //simplex

}

#ifdef SIMPLEX
#ifdef METROPOLIS_HASTINGS_ALGORITHM
 void binaryMCMC::translationSimplex(double** tempGenome, double* g, double* euclidian){

   double tempSimplex[sizeCoupling][size];
   unsigned long position[sizeCoupling];
   for(unsigned long i = 0; i < sizeCoupling; i++){
     position[i] = i;
     for(unsigned long j = 0; j < size; j++)
       tempSimplex[i][j] = tempGenome[i][j];
   }

   //sort the individuals after their fitness
   double fitnessValues[sizeCoupling];
   for(unsigned long i = 0; i < sizeCoupling; i++)
     fitnessValues[i] = fitness(tempGenome[i],size);

   
   for(unsigned long i = 0; i < sizeCoupling; i++)
     for(unsigned long j = 0; j < sizeCoupling; j++)
       if(fitnessValues[i] < fitnessValues[j]){
	 double temp = fitnessValues[i];
	 fitnessValues[i] = fitnessValues[j];
	 fitnessValues[j] = temp; 
	 
	 temp = position[i];
	 position[i] = position[j];
	 position[j] = temp; 
	 
	 for(unsigned long k = 0; k < size; k++){
	   temp = tempSimplex[i][k];
	   tempSimplex[i][k] = tempSimplex[j][k];
	   tempSimplex[j][k] = temp; 
	 }
       }
   
   //REFLECTION CONTRACTION
   double  randValue = genrand_real1(); 
   
   if(randValue < 0.25){
     double transl1 = fac; //* gaussian();
     double fac1 = (1.0 - transl1)/size;
     double fac2 = fac1 - transl1;
     //transl1 = gaussian();
     randValue = genrand_real1(); 
     
     //
     for( unsigned long i =0; i < sizeCoupling; i++){
       if(randValue <= proportions[i]){
	 double sumDimension[size];
	 for(unsigned long j = 0; j < size; j++)
	   sumDimension[j] = 0;
	 
	 for(unsigned long k = 0; k < sizeCoupling; k++)
	   for(unsigned long j = 0; j < size; j++)
	     //if(k != i)
	     sumDimension[j] += tempSimplex[k][j];
	 
	 //make the transformation
	 for(unsigned long j = 0; j < size; j++)
	   tempSimplex[i][j] = sumDimension[j] * fac1 - tempSimplex[i][j] * fac2;
	 
	 for(unsigned long k = 0; k < sizeCoupling; k++)
	   for(unsigned long j = 0; j < size; j++)
	     //if(k != i)
	     sumDimension[j] += tempSimplex[k][j];
	 
	 randValue = genrand_real1(); 
	 
	 if(randValue < 0.5){
	   for(unsigned long j = 0; j < size; j++)
	     tempGenome[position[i]][j] = tempSimplex[i][j]; //* transl1;
	   
	   positionChanged = position[i];
	   rank = i;
	   
	   /*cout << " order { ";
	     for(unsigned long j = 0; j < sizeCoupling; j++)
	     cout << " ( " << position[j] << " , " << fitnessValues[j] << " )\t";
	     cout << " } =  positionChanged " << positionChanged << " rank = " << rank << " \n"; */
	   
	   //only one point is changed
	   return; 
	 }
	 else if(randValue < 0.75){
	   fac = 2;
	 } else {
	   fac = 0.5;
	 }
	 fac1 = (1.0 - transl1)/size;
	 fac2 = fac1 - transl1;
	 //make the transformation
	 for(unsigned long j = 0; j < size; j++)
	   tempSimplex[i][j] = sumDimension[j] * fac1 - tempSimplex[i][j] * fac2;
	 
	 for(unsigned long j = 0; j < size; j++)
	   tempGenome[position[i]][j] = tempSimplex[i][j]; //* transl1;
	 
	 positionChanged = position[i];
	 rank = i;
	 
	 /*cout << " order { ";
	   for(unsigned long j = 0; j < sizeCoupling; j++)
	   cout << " ( " << position[j] << " , " << fitnessValues[j] << " )\t";
	   cout << " } =  positionChanged " << positionChanged << " rank = " << rank << " \n"; */
	 
	 //only one point is changed
	 return; 
       }
     }
   } else if(randValue < 0.75){
     for( unsigned long i =0; i < sizeCoupling; i++){
       if(randValue <= proportions[i]){

	 for(unsigned long k1 = 0; k1 < sizeCoupling; k1++)
	   for(unsigned long k2 = k1 + 1; k2 < sizeCoupling; k2++)
	     if(i != k1 && i != k2){
	       double transl = gaussian() * sigma_transl;
	       
	       for(unsigned long j = 0; j < size; j++)
		 tempSimplex[i][j] += transl * (tempGenome[k1][j] - tempGenome[k2][j]);;
	     }
	 
	 for(unsigned long j = 0; j < size; j++)
	   tempGenome[position[i]][j] = tempSimplex[i][j]; //* transl1;
	 
	 positionChanged = position[i];
	 rank = i;

	 return;
       }
     }     
   } else {
     //rotation
     for( unsigned long j =0; j < sizeCoupling; j++){
       if(randValue <= proportions[j]){

	 //
	 Matrix intermediate(size,1);
	 Matrix y(size,1);
	 Matrix rotate(size - 1,1);
//circumcernter
      double mean_00 = (tempGenome[0][0] + tempGenome[1][0])/2;
      double mean_01 = (tempGenome[0][1] + tempGenome[1][1])/2;
      double slope_0 = -(tempGenome[0][1] - tempGenome[1][1])/(tempGenome[0][0] + tempGenome[1][0]);
      double mean_10 = (tempGenome[0][0] + tempGenome[2][0])/2;
      double mean_11 = (tempGenome[0][1] + tempGenome[2][1])/2;
      double slope_1 = -(tempGenome[0][1] - tempGenome[2][1])/(tempGenome[0][0] + tempGenome[2][0]);

      g[0] = (mean_00/slope_0 - mean_10/slope_1 - mean_01 + mean_11)/(1/slope_0 - 1/slope_1); 
      g[1] = (mean_01*slope_0 - mean_11*slope_1 - mean_00 + mean_10)/(slope_0 - slope_1); 
      for(int i = 2; i < size; i++){
	  double slope = -(tempGenome[0][i] - tempGenome[i+1][i])/(tempGenome[0][0] + tempGenome[i+1][0]);
	  double mean_2 = (tempGenome[0][i] + tempGenome[i+1][i])/2;
	  double mean_1 = (tempGenome[0][0] + tempGenome[i+1][0])/2;
	  g[i]=mean_2 - slope * (g[0] - mean_1);
      }
	 
	 for(int i = 0; i < size; i++)
	   y.element(i,0) = g[i];
  
	 double transl;
	 //what are parameters for rotation ?!?
	 //cout << "rotation angles \n";
	 for(int i =0; i < size - 1;i++){
	   //  if(size == 2)
	   //transl = gaussian() * pi * sigma_rotate;
	   //else 
	   transl = gaussian() * pi/2.0 * sigma_rotate;
	   //transl = gaussian() * pi * sigma_rotate;
	   rotate.element(i,0) = transl;
	 }
  
	 for(int i = 0; i < size; i++)
	   intermediate.element(i,0) = tempSimplex[j][i];
    
	 intermediate = CartesianToPolar(intermediate, y);
	 
	 intermediate = rotationPolar(intermediate,rotate);
    
	 intermediate = PolarToCartesian(intermediate,y);

	 for(int i = 0; i < size; i++){
#ifndef SUM_ROT
	   tempSimplex[j][i] = intermediate.element(i,0);
#else
	   tempSimpex[j][i] = intermediate.element(i,0) + gaussian() * 0.1;
#endif //sum rotation
	 }
       
	 for(unsigned long i = 0; i < size; i++)
	   tempGenome[position[j]][i] = tempSimplex[j][i]; //* transl1;
	 
	 positionChanged = position[j];
	 rank = j;

	 return;
       }
     }
   }
 } 
#else //metropolis hastihn
 void binaryMCMC::translationSimplex(double** tempGenome, double* g, double* euclidian){

   double tempSimplex[sizeCoupling][size];
   unsigned long position[sizeCoupling];
   for(unsigned long i = 0; i < sizeCoupling; i++){
     position[i] = i;
     for(unsigned long j = 0; j < size; j++)
       tempSimplex[i][j] = tempGenome[i][j];
   }
   
   //REFLECTION CONTRACTION
   double  randValue = genrand_real1(); 
   
   if(randValue < 0.05){
     double transl1 = fac; //* gaussian();
     double fac1 = (1.0 - transl1)/size;
     double fac2 = fac1 - transl1;
     //transl1 = gaussian();
     randValue = genrand_real1(); 
     
     //
     for( unsigned long i =0; i < sizeCoupling; i++){
       if(randValue <= proportions[i]){
	 double sumDimension[size];
	 for(unsigned long j = 0; j < size; j++)
	   sumDimension[j] = 0;
	 
	 for(unsigned long k = 0; k < sizeCoupling; k++)
	   for(unsigned long j = 0; j < size; j++)
	     //if(k != i)
	     sumDimension[j] += tempSimplex[k][j];
	 
	 //make the transformation
	 for(unsigned long j = 0; j < size; j++)
	   tempSimplex[i][j] = sumDimension[j] * fac1 - tempSimplex[i][j] * fac2;
	 
	 for(unsigned long k = 0; k < sizeCoupling; k++)
	   for(unsigned long j = 0; j < size; j++)
	     //if(k != i)
	     sumDimension[j] += tempSimplex[k][j];
	 
	 randValue = genrand_real1(); 
	 
	 if(randValue < 0.5){
	   for(unsigned long j = 0; j < size; j++)
	     tempGenome[position[i]][j] = tempSimplex[i][j]; //* transl1;
	   
	   positionChanged = position[i];
	   rank = i;
	   
	   /*cout << " order { ";
	     for(unsigned long j = 0; j < sizeCoupling; j++)
	     cout << " ( " << position[j] << " , " << fitnessValues[j] << " )\t";
	     cout << " } =  positionChanged " << positionChanged << " rank = " << rank << " \n"; */
	   
	   //only one point is changed
	   return; 
	 }
	 else if(randValue < 0.25){
	   fac = 2;
	 } else {
	   fac = 0.5;
	 }
	 fac1 = (1.0 - transl1)/size;
	 fac2 = fac1 - transl1;
	 //make the transformation
	 for(unsigned long j = 0; j < size; j++)
	   tempSimplex[i][j] = sumDimension[j] * fac1 - tempSimplex[i][j] * fac2;
	 
	 for(unsigned long j = 0; j < size; j++)
	   tempGenome[position[i]][j] = tempSimplex[i][j]; //* transl1;
	 
	 positionChanged = position[i];
	 rank = i;
	 
	 /*cout << " order { ";
	   for(unsigned long j = 0; j < sizeCoupling; j++)
	   cout << " ( " << position[j] << " , " << fitnessValues[j] << " )\t";
	   cout << " } =  positionChanged " << positionChanged << " rank = " << rank << " \n"; */
	 
	 //only one point is changed
	 return; 
       }
     }
   } else if(randValue < 0.5){
     for( unsigned long i =0; i < sizeCoupling; i++){
       if(randValue <= proportions[i]){

	 for(unsigned long k1 = 0; k1 < sizeCoupling; k1++)
	   for(unsigned long k2 = k1 + 1; k2 < sizeCoupling; k2++)
	     if(i != k1 && i != k2){
	       double transl = gaussian() * sigma_transl;
	       
	       for(unsigned long j = 0; j < size; j++)
		 tempSimplex[i][j] += transl * (tempGenome[k1][j] - tempGenome[k2][j]);;
	     }
	 
	 for(unsigned long j = 0; j < size; j++)
	   tempGenome[position[i]][j] = tempSimplex[i][j]; //* transl1;
	 
	 positionChanged = position[i];
	 rank = i;

	 return;
       }
     }     
   } else {
     //rotation
     for( unsigned long j =0; j < sizeCoupling; j++){
       if(randValue <= proportions[j]){

	 //
	 Matrix intermediate(size,1);
	 Matrix y(size,1);
	 Matrix rotate(size - 1,1);
	 
//circumcenter here 
      double mean_00 = (tempGenome[0][0] + tempGenome[1][0])/2;
      double mean_01 = (tempGenome[0][1] + tempGenome[1][1])/2;
      double slope_0 = -(tempGenome[0][1] - tempGenome[1][1])/(tempGenome[0][0] + tempGenome[1][0]);
      double mean_10 = (tempGenome[0][0] + tempGenome[2][0])/2;
      double mean_11 = (tempGenome[0][1] + tempGenome[2][1])/2;
      double slope_1 = -(tempGenome[0][1] - tempGenome[2][1])/(tempGenome[0][0] + tempGenome[2][0]);

      g[0] = (mean_00/slope_0 - mean_10/slope_1 - mean_01 + mean_11)/(1/slope_0 - 1/slope_1); 
      g[1] = (mean_01*slope_0 - mean_11*slope_1 - mean_00 + mean_10)/(slope_0 - slope_1); 
      for(int i = 2; i < size; i++){
	  double slope = -(tempGenome[0][i] - tempGenome[i+1][i])/(tempGenome[0][0] + tempGenome[i+1][0]);
	  double mean_2 = (tempGenome[0][i] + tempGenome[i+1][i])/2;
	  double mean_1 = (tempGenome[0][0] + tempGenome[i+1][0])/2;
	  g[i]=mean_2 - slope * (g[0] - mean_1);
      }
	 for(int i = 0; i < size; i++)
	   y.element(i,0) = g[i];
  
	 double transl;
	 //what are parameters for rotation ?!?
	 //cout << "rotation angles \n";
	 for(int i =0; i < size - 1;i++){
	   //  if(size == 2)
	   //transl = gaussian() * pi * sigma_rotate;
	   //else 
	   transl = gaussian() * pi/2.0 * sigma_rotate;
	   //transl = gaussian() * pi * sigma_rotate;
	   rotate.element(i,0) = transl;
	 }
  
	 for(int i = 0; i < size; i++)
	   intermediate.element(i,0) = tempSimplex[j][i];
    
	 intermediate = CartesianToPolar(intermediate, y);
	 
	 intermediate = rotationPolar(intermediate,rotate);
    
	 intermediate = PolarToCartesian(intermediate,y);

	 for(int i = 0; i < size; i++){
#ifndef SUM_ROT
	   tempSimplex[j][i] = intermediate.element(i,0);
#else
	   tempSimpex[j][i] = intermediate.element(i,0) + gaussian() * 0.1;
#endif //sum rotation
	 }
       
	 for(unsigned long i = 0; i < size; i++)
	   tempGenome[position[j]][i] = tempSimplex[j][i]; //* transl1;
	 
	 positionChanged = position[j];
	 rank = j;

	 return;
       }
     }
   }
 } 

#endif
#endif 

 //only if there are three parents --> translate according with the last two parents that are then translated with PCTX
 void binaryMCMC::translationDiff(double** tempGenome, double* g, double* euclidian){
   double x_transl[size];

   //cout << "diff recomb \t"; 
   // the translation distance
   //double transl = 0;
   //for(int i = 0; i < size; i++)
   //  transl += pow(tempGenome[1][i] - tempGenome[2][i],2.0);
   //transl = sqrt(transl);

   //the translation variable
   //#ifndef TRANSLATION_COMPLEX
#ifndef TerBraak
   double transl = gaussian() * sigma_transl;
   // the translation
#ifdef HALF_PCTX
   //if(mix_rotation == 0){
     for(int i = 0; i < size; i++)
       tempGenome[0][i] = tempGenome[0][i] + transl * (tempGenome[1][i] - tempGenome[2][i]);
     //} 
     /*else if(mix_mutation == -1){
     double m[size-1];
     for(int i = 0; i < siz1; i++)
       m[i] = (tempGenome[1][i+1] - tempGenome[2][i+1])/(tempGenome[1][0] - tempGenome[2][0]);
     
   } else {
   }*/
#else
     for(int j = 2; j < sizeCoupling; j++){
       transl = gaussian() * sigma_transl;
       for(int i = 0; i < size; i++)
	 tempGenome[0][i] += transl * (tempGenome[1][i] - tempGenome[j][i]);
     }
     // transl = gaussian() * sigma_transl;
     //for(int i = 0; i < size; i++)
     //  tempGenome[0][i] = tempGenome[0][i] + transl * (tempGenome[1][i] - tempGenome[3][i]);

#endif //HALF_PCTX
#else
   //double transl = constTerBraak;

   for(int i = 0; i < size; i++)
     tempGenome[0][i] += constTerBraak * (tempGenome[1][i] - tempGenome[2][i]) + genrand_real1() * 0.0002 - 0.0001;
#endif 
   /*#else
   double transl = gaussian() * sigma_transl;
   double x[size], y[size];

   // the translation
   for(int i = 0; i < size; i++)
     y[i] = transl * (tempGenome[1][i] - tempGenome[2][i]);

   double m[size-1];
   for(int i = 0; i < size-1; i++)
     m[i] = (tempGenome[1][i+1] - tempGenome[2][i+1])/(tempGenome[1][0] - tempGenome[2][0]);

   x[0] = m[0]/(pow(m[0],2.0)+1) *(tempGenome[0][1] + 1.0/m[0] * tempGenome[0][0] - tempGenome[1][1] + m[0] * tempGenome[1][0]);

   for(int i = 1; i < size; i++)
     x[i] = tempGenome[1][i] + m[i-1] * (x[0] - tempGenome[1][0]);

   transl = gaussian() * sigma_transl_oriz;

   for(int i = 0; i < size; i++)
     tempGenome[0][i] = tempGenome[0][i] + y[i] + transl * (x[i] - tempGenome[0][i]);

     #endif*/ //perpendicular
 }

 void binaryMCMC::translationSnooker(double** tempGenome, double* g){
   double m[size-1];

   for(int j = 0; j < sizeCoupling; j++){
     for(int iteratorCycle = 1; iteratorCycle < size; iteratorCycle++)
       m[iteratorCycle-1] = (tempGenome[j][iteratorCycle] - g[iteratorCycle])/(tempGenome[j][0] - g[0]);

     if(genrand_real1() > 0.5)
       tempGenome[j][0] = tempGenome[j][0] + sqrt(m[0] * m[0] / (m[0] * m[0] + 1.0));
     else
       tempGenome[j][0] = tempGenome[j][0] - sqrt(m[0] * m[0] / (m[0] * m[0] + 1.0));

     for(int i = 1; i < size; i++)
       if(genrand_real1() > 0.5)
	 tempGenome[j][i] = tempGenome[j][i] + sqrt(1.0 / (m[i-1] * m[i-1] + 1.0));
       else
	 tempGenome[j][i] = tempGenome[j][i] - sqrt(1.0 / (m[i-1] * m[i-1] + 1.0));
   }
 }

 void binaryMCMC::translationPCTX(double** tempGenome, double* g, double* euclidian, double sizeCouplingLocal){
  // if there are different number of parents recompute the center and euclidian
  //cout << "translation ";
   //#ifdef DEBUG
  if(sizeCoupling != sizeCouplingLocal){
#ifndef wPCTX
    for(int i = 0; i < size; i++){
      g[i] = 0;
      for(unsigned long j = 0; j < sizeCouplingLocal; j++){
	g[i]+= tempGenome[j][i];
      }
      g[i] /= (double) sizeCouplingLocal;
      //cout << " " << g[i];
    }
#endif //wPCTX

    for(int j = 0; j < sizeCouplingLocal; j++){
      euclidian[j] = 0;
      for(int i =0; i < size; i++)
	euclidian[j] += pow(tempGenome[j][i] - g[i],2.0);
      euclidian[j] = sqrt(euclidian[j]);
      //cout << " j= " << j << "\t" << euclidian[j];
    }
  } else {
    //only for DEBUG
    // for(int i = 0; i < size; i++){
    // cout << " " << g[i];
    //}
    //cout << "\n";
  }
  //#endif
    double x_transl[size];

    //0.0 mutation with the same distance in all directions
    double transl;

    /*    for(int i = 0; i < size; i++){
      transl = gaussian() * mutation_sigma;
      for(int j = 0; j < sizeCouplingLocal; j++)
	tempGenome[j][i] += transl;
      g[i]+=transl;
      }*/

    //1.0 translation towars the center of the mass of the triangle
    //slopes computed for the first individual; the other have the same slopes
#ifdef EUCLIDIAN
    transl = gaussian() * sigma_transl * euclidian[0]; // + gaussian() * mutation_sigma; // * euclidian[0];
#else 
    transl = gaussian() * sigma_transl;
#endif 

//#ifdef DIFFERENT_DIMENSIONS
//    cout << "Not implemented different dimension recombination yet \n";
    //   exit(1);
//#else
//#ifndef TRANSLATION_AXES

    for(int i = 0; i < size; i++){
      x_transl[i] = transl * (tempGenome[0][i] - g[i]);
      //cout << "c(" << i << "," << x_transl[i] << "), ";
    }

    /*   for(int i = 0; i < size; i++){
#ifndef wPCTX
      for(int j = 0; j < sizeCouplingLocal; j++){
#else
	for(int j = 0; j < 1; j++){
#endif
	  tempGenome[j][i] = tempGenome[j][i] + x_transl[i];
	}

	//move the center of the point as well
	//for(int i = 0; i < size; i++)

#ifndef wPCTX
	g[i] = g[i] + x_transl[i];
#else
	g[i] = g[i] + x_transl[i]/size;
#endif //wPCTX
      }
    */
      //cout << "transl ("<< tempGenome[0][0] << "," << tempGenome[0][1] << ") , (" << g[0] << "," << g[1] <<") , "
      //   << euclidian[0] << "\n";

#ifdef TRANSLATION_COMPLEX
      double m[size-1];
      //double y[size];
      //if(abs(tempGenome[0][0] - g[0]) >= AwayFromZero){
      for(int iteratorCycle = 1; iteratorCycle < size; iteratorCycle++)
	m[iteratorCycle-1] = (tempGenome[0][iteratorCycle] - g[iteratorCycle])/(tempGenome[0][0] - g[0]);

      //DEBUG purposes
      //2.0 translation perpendicular on the direction towards the center of the mass, d_i,
      //cout << " | perpendiculat to the center \n";
      //reset the translation to 0
      //for(int i = 0; i < size; i++)
      //x_transl[i] = 0;

#ifdef TWO_PARENTS
      if(size == 2){
#ifdef EUCLIDIAN
	transl = gaussian() * sigma_transl_oriz * euclidian[0]; // * euclidian[0];
#else
	transl = gaussian() * sigma_transl_oriz; // * euclidian[0];
#endif 
	//double sumTemp = 0;
	//for(int k = 0; k < size-1; k++)
	//sumTemp +=  1.0/pow(m[k],2.0);
	x_transl[0] += transl * (tempGenome[0][1] - g[1]);//sqrt(abs(transl)/(1 + sumTemp));
	//if(genrand_real1() < 0.5)
	// x_transl[0] = - x_transl[0];
	//for(int k = 1; k < size; k++)
	x_transl[1] += - transl * (tempGenome[0][0] - g[0]); //x_transl[0]/m[k-1];
      }
#else
      if(size == 2){
#ifdef SIMPLE_PERP_TRANSL
      //if(size == 2){
      double intermP = m[0]/(pow(m[0],2.0)+1) *(tempGenome[1][1] + 1.0/m[0] * tempGenome[1][0] - tempGenome[0][1] + m[0] * tempGenome[0][0]);
      //intermP[1] = recomb_table[1][1] - 1.0/m * (intermP[0] - recomb_table[1][0]);
#ifdef EUCLIDIAN 
      double d1 = sqrt(1 + 1.0/pow(m[0],2)) * sqrt(pow(tempGenome[1][0] - intermP,2.0));
      
#ifdef wPCTX
      intermP = m[0]/(pow(m[0],2.0)+1) *(tempGenome[2][1] + 1.0/m[0] * tempGenome[2][0] - tempGenome[0][1] + m[0] * tempGenome[0][0]);
      double d2 = sqrt(1 + 1.0/pow(m[0],2)) * sqrt(pow(tempGenome[2][0] - intermP,2.0));
      transl = gaussian() * sigma_transl_oriz * (d1 +  d2)/2.0; // * euclidian[0];
#else
      transl = gaussian() * sigma_transl_oriz * d1; // * euclidian[0];
#endif 
      x_transl[0] += transl * (tempGenome[0][1] - g[1]); 
// * euclidian[0];//sqrt(abs(transl)/(1 + sumTemp));
      x_transl[1] += - transl * (tempGenome[0][0] - g[0]); 
// * euclidian[0]; //x_transl[0]/m[k-1];
      //x_transl[2] = - transl * (tempGenome[0][0] - g[0]) * (euclidian[1] +  euclidian[2])/2.0; // * euclidian[0]; //x_transl[0]/m[k-1];
      //}
#else
      transl = gaussian() * sigma_transl_oriz;
#ifndef wPCTX
      x_transl[0] += transl * (intermP - tempGenome[1][0]);
      x_transl[1] += -transl * (intermP - tempGenome[1][0])/m[0];
#else
      double intermP1 = m[0]/(pow(m[0],2.0)+1) *(tempGenome[2][1] + 1.0/m[0] * tempGenome[2][0] - tempGenome[0][1] + m[0] * tempGenome[0][0]);
 
      x_transl[0] += transl * (tempGenome[1][0]-intermP + tempGenome[2][0] -intermP1);
      x_transl[1] += -transl * (tempGenome[1][0]-intermP + tempGenome[2][0] -intermP1)/m[0];
#endif 
#endif 
#else
      for(int i = 0; i < size; i++)
	x_transl[i] += transl * (tempGenome[1][i] - tempGenome[2][i]);

#endif  //simple perp transl
      } else{
	//compute for each dimension the perpendicular
	double intermP[size];
	
	for(int inter = 1; inter < sizeCouplingLocal; inter++){  	  
	  transl = gaussian() * sigma_transl_oriz;
	  
  	  intermP[0] = m[0]/(pow(m[0],2.0)+1) *(tempGenome[inter][1] + 1.0/m[0] * tempGenome[inter][0] - tempGenome[0][1] + m[0] * tempGenome[0][0]);

	  for(int inter_x = 1; inter_x < size; inter_x++)
	    intermP[inter_x] = tempGenome[inter][inter_x] - 1 / m[inter_x-1] * (intermP[0] - tempGenome[inter][0]);

	  for(int inter_x = 0; inter_x < size; inter_x++)
	    x_transl[inter_x] += transl * (intermP[inter_x] - tempGenome[inter][inter_x]);
	} 
      }
/*#else
//2.0 translatarea perpendiculara
//incepe de la al doilea individ
for(int j = 1; j < sizeCouplingLocal; j++){
//
	transl = gaussian() * sigma_transl_oriz; // * euclidian[0];

	//2.0.1 compute the perpendicular distance from a parent x_j to the distance d_i
	if(abs(m[0]) < AwayFromZero){
	  x_transl[0] += transl * (g[0] - tempGenome[j][0]);
	  for(int i = 1; i < size; i++)
	    x_transl[i] += 0;
	} else if(m[0] - (tempGenome[j][1] - tempGenome[0][1])/(tempGenome[j][0] - tempGenome[0][0]) < AwayFromZero){
	  //are only two parents and/or the center of mass is on a line between two parents
	  //a free dimension
	  for(int k = 0; k < size; k++){
	    if(abs((tempGenome[j][0] - g[0]) / (tempGenome[j][k] - g[k])) > AwayFromZero && abs(tempGenome[j][k] - g[k]) > AwayFromZero &&
	       abs(tempGenome[j][0] - g[0]) > AwayFromZero && abs((tempGenome[j][k] - g[k]) / (tempGenome[j][0] - g[0])) > AwayFromZero)
	      x_transl[k] += - transl * (tempGenome[j][0] - g[0]) / (tempGenome[j][k] - g[k]);
	    //else 0
	    //cout << "a free dimension " << (tempGenome[j][0] - g[0]) << "\t "<< (tempGenome[j][k] - g[k]) << "\t"
	    //	 << abs((tempGenome[j][0] - g[0]) / (tempGenome[j][k] - g[k])) << "\n";
	  }
	} else {
	  x_transl[0] += transl * (- m[0] * m[0]/(1 + m[0] * m[0]) * (tempGenome[j][0] - g[0]) + m[0]/(1 + m[0] * m[0]) * (tempGenome[j][1] - g[1]));
	  for(int i = 1; i < size; i++){
	    //cout << "m[" << i-1 << "]="<<m[i-1]<< "\n";
	    x_transl[i] += transl * ( m[i-1] / (1 + m[i-1] * m[i-1]) * (tempGenome[j][0] - g[0]) - (tempGenome[j][i] - g[i])/(1 + m[i-1] * m[i-1]) );
	  }
	}
      }
      #endif*/ 
#endif
#endif //translation complex
      //cout << " transl vector (";
      //for(int i = 0; i < size; i++)
      //cout << x_transl[i] << ",";
      //cout << ")\n";

      //translates the candidates
    for(int i = 0; i < size; i++)
#ifndef wPCTX
      for(int j = 0; j < sizeCouplingLocal; j++){
#else
	for(int j = 0; j < 1; j++){
#endif
	  tempGenome[j][i] += x_transl[i];
	}

	//translate the center
	  for(int i = 0; i < size; i++)
	    g[i] = g[i] + x_transl[i];

	  /*	  cout << "actual parameters after second translation";
	  for(int j = 0; j < sizeCouplingLocal; j++){
	    cout << j << " =( ";
	    for(int i = 0; i < size; i++)
	      cout << tempGenome[j][i] <<",";
	    cout << ")\t";
	    }*/
/*#else //TRANSLATION_AXES

	//translation parallel with axes
    double x_transl[size];
    for(int i = 0; i < size; i++){
      x_transl[i] = gaussian() * sigma_transl * (tempGenome[0][i] - g[i]);
    }
    
    //cout << " translation with " << x_transl[0] << "\t" << x_transl[1] << "\n"; 
 
    for(int i = 0; i < size; i++)
#ifndef wPCTX
      for(int j = 0; j < sizeCoupling; j++)
#else
	for(int j =0; j < 1; j++)
#endif //PCTX
	tempGenome[j][i] = tempGenome[j][i] + x_transl[i];
      

    //move the center of the point as well
    //double g1[size];
    for(int i = 0; i < size; i++)
      g[i] = g[i] + x_transl[i];

    //cout << " new coordinates after translation (" << tempGenome[0][0] << "," << tempGenome[0][1] << ")\t(" 
    //	 << tempGenome[1][0] << "," << tempGenome[1][1] << ")\t center ("
    //	 << g[0] << "," << g[1] << ")\n";

#endif // TRANSLATION_AXES 
#endif //different dimensions*/
}
#endif //translation
 
    //#ifdef ROTATION
void binaryMCMC::rotation(double** tempGenome, double* g, double* euclidian){
#ifdef SINGLE_INDIV_ROT
  mix_rotation = 1; 
#endif 
  //cout << "transl\t";
   //   cout << "should enter here\n";
   // difference recombination on three individuals
   //first individual is translated
  /*#ifdef HALF_PCTX
   translationDiff(tempGenome, g, euclidian);

   //cycles on individuals between PCTX and Diff
   rotation_((tempGenome + 1), (g+1), (euclidian + 1), sizeCoupling - 1);
   #else*/
   rotation(tempGenome, g, euclidian, sizeCoupling);
   //#endif //half PCTX
   //#endif //PCTX
   //#endif //symmetrical snooker
}

//rotate the individuals; ATENTION!!! only for 2 dimensions
// we assume that the distances remains the same
// the center remains the same
//initilize matrices for the roation 
Matrix binaryMCMC::CartesianToPolar(Matrix& rotate, Matrix around){
  //compute the euclidian distance between the individual and the center
  Matrix temp(size-1,1);
  Matrix tempValueMatrix(1,1);
  tempValueMatrix = (rotate - around).t() * (rotate - around);
  //cout << " Cartezian =(" << rotate.t() << ") arround =(" << around.t() << ") (rotate - arround) = (" << (rotate - around).t() << ") square = (" 
  //   << tempValueMatrix << ")\n";  

  tempValueMatrix.element(0,0) = sqrt(tempValueMatrix.element(0,0));

  if(size != 2){
  if(rotate.element(0,0) - around.element(0,0) < 0)
    tempValueMatrix.element(0,0) = - tempValueMatrix.element(0,0);

  //cout << " compute temp \n ";
  for(int i = 0; i < size-1; i++){
    double tempCos = 1;
    for(int j = 0; j < i; j++)
      tempCos *= cos(temp.element(j,0));
    temp.element(i,0) = asin((rotate.element(size - i - 1,0) - around.element(size - i - 1,0))/(tempValueMatrix.element(0,0) * tempCos));
    //cout << "("  << size << ",index = ("<< i  << "," << size - i - 1 << ")," << temp.element(i,0) << ",(" 
    // << rotate.element(size - i - 1,0) - around.element(size - i - 1,0) << " < " << tempValueMatrix.element(0,0) << "), cos = " << tempCos << ", sin =" 
    // << (rotate.element(size - i - 1,0) - around.element(size - i - 1,0))/(tempValueMatrix.element(0,0) * tempCos)<< ")\t";
  }
  //cout << "\n";
  } else {
    temp.element(0,0) = acos((rotate.element(0,0)- around.element(0,0))/tempValueMatrix.element(0,0));
    if(rotate.element(1,0) - around.element(1,0) < 0) temp.element(0,0) = - temp.element(0,0);
  }
  //cout << " transform from cartesian to polar \n";
  rotate.element(0,0) = tempValueMatrix.element(0,0);
  //cout << "( 0, " << rotate.element(0,0) << ")\t";
  for(int i = 0; i < size-1; i++){
    rotate.element(i+1,0) = temp.element(i,0);
    //if(rotate.element(i+1,0) < 0) rotate.element(i+1,0) = 2 * pi + rotate.element(i+1,0);
    //cout << " ( " << i+1 << "," << rotate.element(i+1,0) << "),sin " << sin(rotate.element(i+1,0)) << ",cos " << cos(rotate.element(i+1,0)) << ")\t";
  }
  //cout << "\n";

  return rotate;
}

/*Matrix binaryMCMC::CartesianToPolar(Matrix rotate, Matrix around){
  //compute the euclidian distance between the individual and the center
  Matrix temp(size-1,1);
  Matrix tempValueMatrix(1,1);
  tempValueMatrix = (rotate - around).t() * (rotate - around);
  //cout << " rotate " << rotate << " arround " << around << " rotate - arround " << rotate - around << " square " << tempValueMatrix << "\n";  

  tempValueMatrix.element(0,0) = sqrt(tempValueMatrix.element(0,0));

  if(rotate.element(0,0) - around.element(0,0) < 0)
    tempValueMatrix.element(0,0) = - tempValueMatrix.element(0,0);

  //cout << " compute temp \n ";
  for(int i = 0; i < size-1; i++){
    double tempCos = 1;
    for(int j = 0; j < i; j++)
      tempCos *= cos(temp.element(j,0));
    temp.element(i,0) = asin((rotate.element(size - i - 1,0) - around.element(size - i - 1,0))/(tempValueMatrix.element(0,0) * tempCos));
    //cout << "("  << size << ",index = ("<< i  << "," << size - i - 1 << ")," << temp.element(i,0) << ",(" 
    // << rotate.element(size - i - 1,0) - around.element(size - i - 1,0) << " < " << tempValueMatrix.element(0,0) << ")," << tempCos << "," 
    // << (rotate.element(size - i - 1,0) - around.element(size - i - 1,0))/(tempValueMatrix.element(0,0) * tempCos)<< ")\t";
  }
  //cout << "\n";

  // cout << " transform from cartesian to polar \n";
  rotate.element(0,0) = tempValueMatrix.element(0,0);
  //cout << "( 0, " << rotate.element(0,0) << ")\t";
  for(int i = 0; i < size-1; i++){
    rotate.element(i+1,0) = temp.element(i,0);
    //cout << " ( " << i << "," << rotate.element(i+1,0) << ")\t";
  }
  //cout << "\n";

  return rotate;
  }*/
    
Matrix binaryMCMC::rotationPolar(Matrix& rotate, Matrix angle){
  //cout << " old rotation angles \n";
  //for(int i = 0; i < size; i++)
  //cout << " ( " << i << "," << rotate.element(i,0) << ")\t";
  //cout << "\n";

  for(int i = 1; i < size; i++)
    rotate.element(i,0) += angle.element(i-1,0); 

  //cout << " new rotate angle \n";
  //for(int i = 0; i < size; i++)
  //if(i != 0)
  //  cout << " ( " << i << "," << rotate.element(i,0) << "," << angle.element(i-1,0)<< ")\t";
  //else
  //  cout << " ( " << i << "," << rotate.element(i,0) << ")\t";
  //cout << "\n";

  return rotate;
}

Matrix binaryMCMC::PolarToCartesian(Matrix& rotate, Matrix around){
  //compute the euclidian distance between the individual and the center
  Matrix temp(size,1);

  if(size != 2){
    double tempCos = 1;
    for(int i = 0; i < size - 1; i++)
      tempCos *= cos(rotate.element(i + 1,0));
    temp.element(0,0) = around.element(0,0) + rotate.element(0,0) * tempCos;
    
    // cout << " ( 0, " << tempCos << " , " << temp.element(0,0) << ") \t";
    
    for(int i = 1; i < size; i++){
      tempCos = 1;
      for(int j = 1; j < size - i; j++)
	tempCos *= cos(rotate.element(j,0));
      temp.element(i,0) = around.element(i,0) + rotate.element(0,0) * tempCos * sin(rotate.element(size - i,0));
    }
  } else {
    temp.element(0,0) = rotate.element(0,0) * cos(rotate.element(1,0));
    temp.element(1,0) = rotate.element(0,0) * sin(rotate.element(1,0));
  }
  //cout << "polar to cartezian \t";
  //for(int i = 0; i < size; i++)
  // cout << "(" << i << " , " << temp.element(i,0) << ") \t";
  //cout << "\n";
  rotate = temp + around;
  
  return rotate;
}

void binaryMCMC::rotation(double** tempGenome, double* g, double* euclidian, unsigned long sizeCoup){
#ifdef SINGLE_INDIV_ROT
  mix_rotation = 1; 
#endif 
  if(sizeCoupling != sizeCoup){
#ifndef wPCTX
      //for(int i = 0; i < size; i++){
      // g[i] = 0;
      //  for(unsigned long j = 0; j < sizeCoup; j++){
      //	g[i]+= tempGenome[j][i];
      //}
      //g[i] /= (double) sizeCoup;
      //cout << " " << g[i];
//circumcenter
      double mean_00 = (tempGenome[0][0] + tempGenome[1][0])/2;
      double mean_01 = (tempGenome[0][1] + tempGenome[1][1])/2;
      double slope_0 = -(tempGenome[0][1] - tempGenome[1][1])/(tempGenome[0][0] + tempGenome[1][0]);
      double mean_10 = (tempGenome[0][0] + tempGenome[2][0])/2;
      double mean_11 = (tempGenome[0][1] + tempGenome[2][1])/2;
      double slope_1 = -(tempGenome[0][1] - tempGenome[2][1])/(tempGenome[0][0] + tempGenome[2][0]);

      g[0] = (mean_00/slope_0 - mean_10/slope_1 - mean_01 + mean_11)/(1/slope_0 - 1/slope_1); 
      g[1] = (mean_01*slope_0 - mean_11*slope_1 - mean_00 + mean_10)/(slope_0 - slope_1); 
      for(int i = 2; i < size; i++){
	  double slope = -(tempGenome[0][i] - tempGenome[i+1][i])/(tempGenome[0][0] + tempGenome[i+1][0]);
	  double mean_2 = (tempGenome[0][i] + tempGenome[i+1][i])/2;
	  double mean_1 = (tempGenome[0][0] + tempGenome[i+1][0])/2;
	  g[i]=mean_2 - slope * (g[0] - mean_1);
      }
#else
    cout << " EROR in binaryMCMC_trap.cpp::3621 \n";
    exit(1);
#endif //wPCTX

    for(int j = 0; j < sizeCoup; j++){
      euclidian[j] = 0;
      for(int i =0; i < size; i++)
	euclidian[j] += pow(tempGenome[j][i] - g[i],2.0);
      euclidian[j] = sqrt(euclidian[j]);
      //cout << " j= " << j << "\t" << euclidian[j];
    }
  } else {
     double mean_00 = (tempGenome[0][0] + tempGenome[1][0])/2;
      double mean_01 = (tempGenome[0][1] + tempGenome[1][1])/2;
      double slope_0 = -(tempGenome[0][1] - tempGenome[1][1])/(tempGenome[0][0] + tempGenome[1][0]);
      double mean_10 = (tempGenome[0][0] + tempGenome[2][0])/2;
      double mean_11 = (tempGenome[0][1] + tempGenome[2][1])/2;
      double slope_1 = -(tempGenome[0][1] - tempGenome[2][1])/(tempGenome[0][0] + tempGenome[2][0]);

      g[0] = (mean_00/slope_0 - mean_10/slope_1 - mean_01 + mean_11)/(1/slope_0 - 1/slope_1); 
      g[1] = (mean_01*slope_0 - mean_11*slope_1 - mean_00 + mean_10)/(slope_0 - slope_1); 
      for(int i = 2; i < size; i++){
	  double slope = -(tempGenome[0][i] - tempGenome[i+1][i])/(tempGenome[0][0] + tempGenome[i+1][0]);
	  double mean_2 = (tempGenome[0][i] + tempGenome[i+1][i])/2;
	  double mean_1 = (tempGenome[0][0] + tempGenome[i+1][0])/2;
	  g[i]=mean_2 - slope * (g[0] - mean_1);
      }
    //only for DEBUG
    //for(int i = 0; i < size; i++){
      //cout << " " << g[i];
      //}
      //cout << "\n";

  }

  Matrix intermediate(size,1);
  Matrix y(size,1);
  Matrix rotate(size - 1,1);

  for(int i = 0; i < size; i++)
    y.element(i,0) = g[i];
  
  double transl;
  //what are parameters for rotation ?!?
  //cout << "rotation angles \n";
  for(int i =0; i < size - 1;i++){
    //  if(size == 2)
    //transl = gaussian() * pi * sigma_rotate;
    //else 
    transl = gaussian() * pi/2.0 * sigma_rotate;
    //transl = gaussian() * pi * sigma_rotate;
#ifndef ROTATION_PROP 
#ifdef EUCLIDIAN
    rotate.element(i,0) = transl * euclidian[i];
#else
    rotate.element(i,0) = transl;
#endif 
#else
#ifdef EUCLIDIAN
    if(abs((tempGenome[0][0] - g[0])/(tempGenome[0][i+1] - g[i+1])) <= 1)
      rotate.element(i,0) = transl * (1 - sqrt(pow((tempGenome[0][0] - g[0])/(tempGenome[0][i+1] - g[i+1]),2.0)))* euclidian[i];
    else 
      rotate.element(i,0) = transl * (1 - sqrt(pow((tempGenome[0][i+1] - g[i+1])/(tempGenome[0][0] - g[0]),2.0)))* euclidian[i];
#else
    if(abs((tempGenome[0][0] - g[0])/(tempGenome[0][i+1] - g[i+1])) <= 1)
      rotate.element(i,0) = transl * (1 - sqrt(pow((tempGenome[0][0] - g[0])/(tempGenome[0][i+1] - g[i+1]),2.0)));
    else 
      rotate.element(i,0) = transl * (1 - sqrt(pow((tempGenome[0][i+1] - g[i+1])/(tempGenome[0][0] - g[0]),2.0)));
#endif 
    /*if(abs((tempGenome[0][0] - g[0])/(tempGenome[0][i+1] - g[i+1])) <= 1)
      rotate.element(i,0) = transl * sigma_rotate * sqrt(pow((tempGenome[0][0] - g[0])/(tempGenome[0][i+1] - g[i+1]),2.0)) * euclidian[i];
    else 
    rotate.element(i,0) = transl * sigma_rotate * sqrt(pow((tempGenome[0][i+1] - g[i+1])/(tempGenome[0][0] - g[0]),2.0)) * euclidian[i];*/
#endif //rotation prop
    //cout << "rotation with " << rotate.element(i,0) << "\n"; 
  }
  
  //if(sizeCoupling != 3){
  //#ifndef SINGLE_INDIV_ROT
#ifdef SINGLE_INDIV_ROT
  for(int j = 0; j < 1; j++){
#else
  for(int j = 0; j < sizeCoup; j++){
#endif
    //#else
    //for(int j = 0; j < 1; j++){
    //#endif 
    for(int i = 0; i < size; i++)
      intermediate.element(i,0) = tempGenome[j][i];
    
    //if(size != 2){
    /*cout << " intermediate in cartesian \n";
    for(int i = 0; i < size; i++)
    cout << "(" << i << "," << intermediate.element(i,0) << ")\t";
    cout << "\n";*/

    intermediate = CartesianToPolar(intermediate, y);

    intermediate = rotationPolar(intermediate,rotate);

    /*cout << " intermediate in polar \n";
    for(int i = 0; i < size; i++)
    cout << "(" << i << "," << intermediate.element(i,0) << ")\t";
    cout << "\n";*/
    //for each dimension rotate separatly
    //if(sizeCoupling == 3 && j < 0)
    
    intermediate = PolarToCartesian(intermediate,y);

    /* cout << " intermediate in cartezian again  \n";
    for(int i = 0; i < size; i++)
    cout << "(" << i << "," << intermediate.element(i,0) << ")\t";
    cout << "\n";*/
    /*} else {
      //if(sizeCoupling != 3){
#ifndef TWO_PARENTS
      double transl = gaussian() * pi * sigma_rotate;
#endif 
      intermediate.element(0,0) = g[0] + cos(transl*euclidian[0])*(tempGenome[j][0] - g[0])+sin(transl*euclidian[0]) * (tempGenome[j][1] - g[1]);
	
      intermediate.element(1,0) = g[1] - sin(transl*euclidian[0])*(tempGenome[j][0] - g[0]) + cos(transl*euclidian[0]) * (tempGenome[j][1] - g[1]);
	//} else {
	
	//}
	}*/
    for(int i = 0; i < size; i++){
#ifndef SUM_ROT
      tempGenome[j][i] = intermediate.element(i,0);
#else
      tempGenome[j][i] = intermediate.element(i,0) + gaussian() * 0.001;
#endif //sum rotation
    }
  }
}
  //#endif // rotation

  //#ifdef SCALING
void binaryMCMC::scaling_dim(double** tempGenome, double* g, double* euclidian){
#ifdef SINGLE_INDIV_ROT
  mix_rotation = -1; 
#endif 
  //cout << "transl\t";
   //   cout << "should enter here\n";
   // difference recombination on three individuals
   //first individual is translated
  /*#ifdef HALF_PCTX
   translationDiff(tempGenome, g, euclidian);

   //cycles on individuals between PCTX and Diff
   scaling_dim_((tempGenome + 1), (g+1), (euclidian + 1), sizeCoupling - 1);
   #else*/
   scaling_dim(tempGenome, g, euclidian, sizeCoupling);
   //#endif //half PCTX

}

void binaryMCMC::scaling_dim(double** tempGenome, double* g, double* euclidian, unsigned long sizeCoup){
    //scaling
    //cout << "resulting individuals before scale(" << tempGenome[0][0] << "," << tempGenome[0][1] << ") ; (" 
  // <<  tempGenome[1][0] << "," << tempGenome[1][1] << "); (" << tempGenome[2][0] << "," << tempGenome[2][1] << ")\n";
#ifdef SINGLE_INDIV_ROT
  mix_rotation = -1; 
#endif 

  if(sizeCoupling != sizeCoup){
#ifndef wPCTX
    for(int i = 0; i < size; i++){
      g[i] = 0;
      for(unsigned long j = 0; j < sizeCoup; j++){
	g[i]+= tempGenome[j][i];
      }
      g[i] /= (double) sizeCoup;
      //cout << " " << g[i];
    }
#endif //wPCTX

    for(int j = 0; j < sizeCoup; j++){
      euclidian[j] = 0;
      for(int i =0; i < size; i++)
	euclidian[j] += pow(tempGenome[j][i] - g[i],2.0);
      euclidian[j] = sqrt(euclidian[j]);
      //cout << " j= " << j << "\t" << euclidian[j];
    }
  } else {
    //only for DEBUG
    // for(int i = 0; i < size; i++){
    // cout << " " << g[i];
    //}
    //cout << "\n";
  }

  //#ifdef SCALING_NO_SYMETRICAL
  //double scale_ = gaussian() * sigma_scale * euclidian[0];
    //cout << "scaling factor " << scale_ * (tempGenome[j][i] - g[i]) / euclidian[i]
  //#else
  //#ifdef SCALING_MUTATION_SAME
    double scale_ = sigma_scale * gaussian();
    //#endif //scaling_mutation_same
    //#endif //scaling

    double m[size-1];
    //double y[size];
    //if(abs(tempGenome[0][0] - g[0]) >= AwayFromZero){
    for(int iteratorCycle = 1; iteratorCycle < size; iteratorCycle++)
      m[iteratorCycle-1] = (tempGenome[0][iteratorCycle] - g[iteratorCycle])/(tempGenome[0][0] - g[0]);
    
#ifdef SINGLE_INDIV_ROT
    for(int j = 0; j < 1; j++){
#else
    for(int j = 0; j < sizeCoup; j++){
#endif
      //double tang = (tableGenomes[j + sizeCoupling][1] - g[1]) / (tableGenomes[j + sizeCoupling][0] - g[0]);
      tempGenome[j][0] = tempGenome[j][0] + scale_ * (m[0] + 1);
      for(int i = 1; i < size; i++)
	tempGenome[j][i] = tempGenome[j][i] + m[i-1] * (tempGenome[j][0] - tempGenome[j][i]); // / euclidian[j];
    }

    //cout << " scale " << scale << " with scale_ " << scale_ << "\n";
    //cout << "resulting individuals after scale(" << tempGenome[0][0] << "," << tempGenome[0][1] << ") ; (" 
    //	 <<  tempGenome[1][0] << "," << tempGenome[1][1] << "); (" << tempGenome[2][0] << "," << tempGenome[2][1] << ")\n";
}
    //#endif //SCALING
#endif //two parents

#endif // real recombination
#endif //IsRecombination


/*#ifndef METROPOLIS_ALGORITHM
double binaryMCMC::ProposalDistribution(int* tempGenome, unsigned long size, double** tempProposalDistributionVector){
  double tempProposalDistribution = 0;
#ifdef INDEPENDENT_SAMPLER
  //cout << "Calculate independent sampler \n";
    tempProposalDistribution = 2*mutation*tempProposalDistributionVector[randJ][tempGenome[randJ]];
  //cout << "\n Final results" << tempProposalDistribution <<", "<<exp(tempProposalDistribution) << "\n";
#else
  tempProposalDistribution = mutation;
#endif
  return tempProposalDistribution;
}

double binaryMCMC::ProposalDistribution(double** tempProposalDistributionVector){
  double tempProposalDistribution = 0;
#ifdef INDEPENDENT_SAMPLER
    tempProposalDistribution = 2*mutation*tempProposalDistributionVector[randJ][genome[randJ]];
    //cout << "\n Final results" << tempProposalDistribution <<", "<<exp(tempProposalDistribution) << "\n";
#else
  tempProposalDistribution = mutation;
#endif
  return tempProposalDistribution;
}

double binaryMCMC::ProposalDistribution(){
  double tempProposalDistribution = 0;
#ifdef INDEPENDENT_SAMPLER
    tempProposalDistribution = 2*mutation*proposalDistributionVector[randJ][genome[randJ]];
  //cout << "\n Final results" << tempProposalDistribution <<", "<<exp(tempProposalDistribution) << "\n";
#else
  tempProposalDistribution = mutation;
#endif
  proposalDistribution = tempProposalDistribution;
  return proposalDistribution;
}

void binaryMCMC::AddaptationMethods(double** tempProposalDistributionVector){
	for(unsigned long i = 0; i < size; i++)
	  for(unsigned long j = 0; j < 2; j++)
	    tempProposalDistributionVector[i][j] = proposalDistributionVector[i][j];
}
#endif // Metropolis_algorithm
*/
unsigned long binaryMCMC::nextGeneration(){
#ifndef MIXTURE_REALS
	unsigned long* newGenome = new unsigned long[size];
#else
	double* newGenome = new double[size];
#endif
	for(unsigned long i = 0; i < size; i++)
	    newGenome[i] = -1;
#ifndef METROPOLIS_ALGORITHM
	double** newProposalDistributionVector = new double*[size];
	for(unsigned long i = 0; i < size; i++)
	  newProposalDistributionVector[i]= new double[2];
#endif

#ifdef DEBUG_
	cout << "Next generation with mutation \n";
#endif
		
#ifdef METROPOLIS_ALGORITHM
	MutationGenome(newGenome);
#else
#ifdef POPULATION_MCMC 
	AddaptationMethods(newProposalDistributionVector);
	MutationGenome(newGenome,newProposalDistributionVector);
#else
	MutationGenome(newGenome);
#endif
#endif //METROPOLIS_ALGORITHM
	
	double newValue = fitness(newGenome,size);
#ifndef METROPOLIS_ALGORITHM
#ifndef MIXTURE_REALS
	double newProposalDistribution = ProposalDistribution(newGenome,size,newProposalDistributionVector);
#endif
#endif

#ifndef RANDOM_WALK
#ifdef METROPOLIS_ALGORITHM
	if(acceptance(newGenome,newValue) == 1) {
#else
#ifndef MIXTURE_REALS
	if(acceptance(newGenome,newValue,newProposalDistribution) == 1) {
#else
	if(acceptance(newGenome,newValue) == 1) {
#endif
#endif
#ifdef DEBUG_
	    cout << "Accept proposed individual with acceptance " << accepted << "\n";
#endif
#endif //RANDOM_WALK
	  //nrOfGenerations++;
	  for(unsigned long i=0; i<size; i++)
	    genome[i] = newGenome[i];
	  value = newValue;
#ifndef RANDOM_WALK 
	} else {
	  
#ifdef DEBUG_
	    cout << "reject proposed individual with acceptance " << accepted << "\n";
	    //for(unsigned long i = 0; i<size; i++)
	    //cout << newGenome[i] <<"|";
	    //cout << " value: "<<newValue <<"\n";	 
#endif
	  //nrOfGenerations++;
	}	
#endif //RANDOM_WALK	
	free(newGenome);
#ifndef METROPOLIS_ALGORITHM
	delete[] newProposalDistributionVector;
#endif
	return 1;
}

#ifndef REALS 
unsigned long binaryMCMC::nextGeneration(unsigned long* nextGenome, double nextValue, unsigned long nextSize){
#else
unsigned long binaryMCMC::nextGeneration(double* nextGenome, double nextValue, unsigned long nextSize){
#endif
	if(size != nextSize) {
		delete genome;
#ifndef REALS 
		genome = new unsigned long[nextSize];
#else
		genome = new double[nextSize];
#endif
		size = nextSize;
	}	
	for(unsigned long i = 0; i<size; i++)
		genome[i] = nextGenome[i];
	value = nextValue;
	/*if(fitness() != nextValue) {
	  cout << "\n!!!!!!!!Error in calculation of fitness" << fitness() <<" = " << nextValue << " generation " << nrOfGenerations << "binaryMCMC_trap.cpp:2758 \n";
	  for(unsigned long i = 0; i < size; i++)
	    cout << genome[i] << ",";
	  cout << "\n";
	  exit(1);
	  value = fitness();
	  }*/
	/*#ifndef METROPOLIS_ALGORITHM
	proposalDistribution = mutation;
	#endif*/
	//nrOfGenerations++;	
	return 1;
}

unsigned long binaryMCMC::nextGenerationBlank(){
  //nrOfGenerations++;	
	return 1;
}

void binaryMCMC::program(char* outFile){
	for(unsigned long i= 0; i<1; i++) {
		nextGeneration();
		nrOfGenerations++;	
	}	
	nextGeneration();	
	nrOfGenerations++;	
	
	for(unsigned long i=0; i<sampleSize-1; i++){	
#ifdef DEBUG_
		printState();
#endif
		addElement();
		for(unsigned long j=0; j<1; j++){
			nextGeneration();
			nrOfGenerations++;	
		}		
	}

	ofstream myFile(outFile);		
	if(!myFile.is_open())
	    cout << "Error opening file "<< outFile << " binaryMCMC.trap:1362\n";
	else print(myFile);
	myFile.close();
}

#ifdef LIST
unsigned long binaryMCMC::addElement(){
  if(genome == NULL) {
    cout << "Eroare obiect \n"; return -1;
    }
  if(index == 0) {
    unsigned long myGenome[size];
    for(unsigned long i = 0; i<size; i++)
      myGenome[i] = genome[i];
    
    double myValue = value;
    double myTemperature = temperature;
    unsigned long myGenerations = nrOfGenerations;

    if(next == NULL){
      //cout << "add element to this level " << index << "\n";

      binaryMCMC* temp = new binaryMCMC(size,myGenome,myTemperature,myValue);
      temp->index = index + 1;
      temp->nrOfGenerations = myGenerations;
#ifdef DEBUG_
      for(unsigned long i = 0; i < size; i++)
	cout << "||" <<temp->genome[i];
      cout << "value "<< temp->value << "\n";	
#endif
      temp->next = NULL;
      next = temp;
      return 1;
    }
    else {
      //cout << "add element to the next level \n";
      return next->addElement(size,myGenome,myTemperature,myGenerations);
    }
  }
  else {
    cout << "Error in adding an element at index="<<index<<"\n";
    return -1;
  }
  return 1;
}

unsigned long binaryMCMC::addElement( unsigned long mySize, unsigned long* myGenome, double myTemperature, unsigned long mynrOfGenerations){
 if(genome == NULL) {
    cout << "Eroare obiect \n"; return -1;
    }
  
  if(next == NULL) 
    {
      //cout << "add element at the level "<< index << "\n";

      binaryMCMC* temp = new binaryMCMC(mySize,myGenome,myTemperature);
      temp->index = index + 1;
      temp->nrOfGenerations = mynrOfGenerations;
#ifdef DEBUG_
      for(unsigned long i = 0; i < mySize; i++)
	cout << "||" <<temp->genome[i];
      cout << "value "<< temp->value << "\n";	
#endif
      temp->next = NULL;
      next = temp;		
      return 1;
    } else {
      //cout << "add element at the next index \n"; 
      return next->addElement(mySize,myGenome,myTemperature,mynrOfGenerations);
    }
  return 1;
}

double binaryMCMC::getElementValue(unsigned long i){
  if(genome == NULL || i <= 0)
    {
      cout << "Eroare obiect " << index << " because genome" << genome << " or i"  << i << "\n";
      return -1;
    }
  if(i == index) return value;
  else 
    if(i != index && next != NULL)
      return next -> getElementValue(i);
    else 
      //if(next == NULL && i != index)
      {
	cout << "Obiect unfound in getElementValue for i="<<i<<" index="<<index<<"\n";
	return 0;
      }
  return -1;
}

double binaryMCMC::getElementValue(){
  if(genome == NULL)
    {
      cout << "Eroare obiect " << index << " because genome" << genome << "\n";
      return -1;
    }
  if(next != NULL && index > 0){
    //cout << "  index =" << index;     
    return next -> getElementValue() + value;
  }
  else if(next == NULL && index > 0){
    //cout << "  index =" << index; 
    return value;
  } else if(next != NULL){
    return next -> getElementValue();
  } else {
    cout << "Error no stored element in this\n";
    return -1;
  }
}

double binaryMCMC::getElementValue(unsigned long n, unsigned long burn){
  if(genome == NULL)
    {
      cout << "Eroare obiect " << index << " because genome" << genome << "\n";
      return -1;
    }
  double tempValue = 0;

  if(next != NULL && index >= burn && index < n){
    //cout << "  index =" << index<< " value " << value << " n=" << n << " burn=" << burn << "\n";     
    return next -> getElementValue(n,burn) + value;
  }else if(index == n){
    //cout << "  index =" << index << " value " << value << " n=" << n << " burn=" << burn << "\n"; 
    return value;
  } else if(next != NULL && index < burn){
    //cout << "   notI ="<< index << " n=" << n << " burn=" << burn <<" \n"; 
    return next -> getElementValue(n,burn);
  } else {
    cout << "Error no stored element in this\n";
    return -1;
  }
  return tempValue;
}

double binaryMCMC::getElementValue(unsigned long n, unsigned long burn, double* teta){
  if(genome == NULL)
    {
      cout << "Eroare obiect " << index << " because genome" << genome << "\n";
      return -1;
    }

  if(next != NULL && index >= burn && index < n){
    teta[index - 1] = value;
    //cout << "  index =" << index<< " value " << value << " n=" << n << " burn=" << burn << "\n";     
    return next -> getElementValue(n,burn, teta) + value;
  }else if(index == n){
    teta[index - 1] = value;
    //cout << "  index =" << index << " value " << value << " n=" << n << " burn=" << burn << "\n"; 
    return value;
  } else if(next != NULL && index < burn){
    //cout << "   notI ="<< index << " n=" << n << " burn=" << burn <<" \n"; 
    return next -> getElementValue(n,burn,teta);
  } else {
    cout << "Error no stored element in this\n";
    return -1;
  }
}

double binaryMCMC::getElementValue(double* teta){
  if(genome == NULL)
    {
      cout << "Eroare obiect " << index << " because genome" << genome << "\n";
      return -1;
    }

  if(next != NULL && index != 0){
    teta[index - 1] = value;
    //cout << "  index =" << index<< " value " << value << " n=" << n << " burn=" << burn << "\n";     
    return next -> getElementValue(teta) + value;
  } else if(next != NULL && index == 0){
    //cout << "   notI ="<< index << " n=" << n << " burn=" << burn <<" \n"; 
    return next -> getElementValue(teta);
  } else if(next == NULL && index != 0){
    teta[index - 1] = value;
    return value;
  } else {
    cout << "Error no stored element in this\n";
    exit(1);
    return -1;
  }
}

double binaryMCMC::getVarianceValue(unsigned long n, unsigned long burn, double mean){
  if(genome == NULL)
    {
      cout << "Eroare obiect " << index << " because genome" << genome << "\n";
      return -1;
    }

  if(next != NULL && index >= burn && index < n){
    //cout << "  index =" << index << " tempValue " << tempValue << "\n";     
    return (value - mean)*(value - mean)+ next -> getVarianceValue(n,burn,mean);
  }else if(index == n){
    //cout << "  index =" << index << " tempValue " << tempValue << "\n"; 
    return (value - mean)*(value - mean);
  } else if(next != NULL && index < burn){
    //cout << "   notI ="<< index << "\n"; 
    return next -> getVarianceValue(n,burn,mean);
  } else {
    cout << "Error no stored element in this\n";
    return -1;
  }
}

double binaryMCMC::getVarianceValue(unsigned long n, unsigned long burn, double mean, unsigned long s){
  if(genome == NULL)
    {
      cout << "Eroare obiect " << index << " because genome" << genome << "\n";
      return -1;
    }

  if(next != NULL && index >= burn && index < n){
    //cout << "  index =" << index << " tempValue " << tempValue << "\n"; 
    if(s % 2 == 0) return pow(value - mean,s) + next -> getVarianceValue(n,burn,mean,s);
    else if(value < mean)
      return -pow(value - mean,s) + next -> getVarianceValue(n,burn,mean,s);
    else return pow(value - mean,s) + next -> getVarianceValue(n,burn,mean,s);
  }else if(index == n){
    //cout << "  index =" << index << " tempValue " << tempValue << "\n"; 
    if(s % 2 == 0) return pow(value - mean, s);
    else if(value < mean) return -pow(value - mean,s);
    else return pow(value - mean,s);
  } else if(next != NULL && index < burn){
    //cout << "   notI ="<< index << "\n"; 
    return next -> getVarianceValue(n,burn,mean,s);
  } else {
    cout << "Error no stored element in this\n";
    return -1;
  }
}

binaryMCMC* binaryMCMC::getElement(unsigned long i){
  if(genome == NULL || i <= 0)
    {
      cout << "Eroare obiect \n";
      return NULL;
    }
  if(i == index) {
    binaryMCMC* temp = new binaryMCMC(size,genome,temperature,value);
    temp -> nrOfGenerations = nrOfGenerations;
    return temp;
  }
  if(i != index && next != NULL)
    return next -> getElement(i);
  if(next == NULL && i != index)
    {
      cout << "Obiect unfound in getElement for i ="<<i<< " index="<<index<<"\n";
      return NULL;
    }
  return NULL;
}

void binaryMCMC::reset(unsigned long myRuns){
  index = 0;
  runs = myRuns;
  if(next != NULL) {
    delete next;
    next = NULL;	
  }
  nrOfGenerations = 0;

  setSeed();

#ifdef MIXTURE_REALS
#ifdef BIVARIATE
#ifndef MIXTURE_BIVARIATE
  unsigned long integers[size*BLOCKsize];
  for(unsigned long i=0;i<size*BLOCKsize;i++)
	  integers[i] = genrand_int32() % 2;
  convertReal(genome,integers,size);
#else //mixture bivariate
  bool true_bit = true;
  while(true_bit){
    for(unsigned long i=0;i<size;i++){
#ifndef INIT_HALF
#ifndef INIT_CENTER
#ifndef INIT_SIDE
      genome[i] = genrand_real1() * scale; 
      genome[i] += -scale/2.;
#else
	genome[i] = genrand_real1() * scale/init_coef;
	if(i == 0)
	  genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
	else genome[i] = genome[i] - (init_coef - 1.0) *scale/2.0/init_coef;
	//
	//genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
#endif //init side
#else
      genome[i] = genrand_real1() * scale/init_coef;
      genome[i] += -scale/2.0/init_coef;
#endif //init_center
#else
      if(i == 0)
	genome[i] = genrand_real1() * scale - scale/2.0;
      else genome[i] = genrand_real1() * scale/2.0;
#endif 
    }
    if(fitness(genome,size) > MIN_FITNESS)
      true_bit = false;
  }
#endif //mixture bivariate
#else //bivariate
  for(unsigned long i=0;i<size;i++){
      genome[i] = genrand_real1() * scale; 
      if(genome[i] < scale/2)
	  genome[i] = -genome[i];
      else genome[i] -= scale/2;
  }
#endif //bivariate
#else //mixture reals 
#ifdef NORMAL_RUNS
  for(unsigned long i=0;i<size;i++){
	  genome[i] = genrand_int32() % 2;
  }
#else
#ifdef BINOMIAL
     for(unsigned long j = 0; j < (size/(double)BLOCKsize); j++){
	 for(unsigned long i = 0; i<BLOCKsize/2; i++)
	     genome[i + j*BLOCKsize] = 0;
	 for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++)
	     genome[i + j*BLOCKsize] = 1;
     } 
#else
#ifdef TWO_ATTRACTORS_FUNCTION
  double count1 = 0;
  for(unsigned long i=0;i<size;i++){
	  if(i - count1 >= size/2) {
	      genome[i] = 1;
	      count1++;
	  } else if(count1 <= size/2){
	      genome[i] = rand() % 2;
	      if(genome[i] == 1)
		  count1++;
	  }
	  else genome[i] = 0;  
  }
#else
  for(unsigned long i=0;i<size;i++){
	  if(runs == 0) genome[i] = 0;
	  else if(runs == 1) genome[i] = 1;
	  else if(runs == 2) { 
	    if(i < size/2) genome[i] = 0;
	    else genome[i] = 1; 
	  } else if(runs == 3) {
	    if(i < size/3) genome[i] = 0;
	    else genome[i] = 1;
	  } else if(runs == 4) {
	    if(i < size/4) genome[i] = 0;
	    else genome[i] = 1;
	  } else if(runs == 5) {
	    if(i < BLOCKsize) genome[i] = 0;
	    else genome[i] = 1;
	  }
	  else genome[i] = genrand_int32() % 2;
  }
#endif
#endif 
#endif
#endif //mixture reals
    //cout << genome[i] << "||";
  //cout << "\n";	
  value = fitness();
}

void binaryMCMC::reset(unsigned long myRuns, unsigned long chain){
  index = 0;
  runs = myRuns;
  if(next != NULL) {
    delete next;
    next = NULL;	
  }
  nrOfGenerations = 0;
  
  setSeed();
#ifdef MIXTURE_REALS
#ifdef BIVARIATE
#ifndef MIXTURE_BIVARIATE
  unsigned long integers[size*BLOCKsize];
  for(unsigned long i=0;i<size*BLOCKsize;i++)
	  integers[i] = genrand_int32() % 2;
  convertReal(genome,integers,size);
  // for(unsigned long i=0;i<size;i++)
//	  genome[i] = genrand_real1();
#else //mixture bivariate
  bool true_bit = true;
  while(true_bit){
    for(unsigned long i=0;i<size;i++){
#ifndef INIT_HALF
#ifndef INIT_CENTER
#ifndef INIT_SIDE
      genome[i] = genrand_real1() * scale; 
      genome[i] += -scale/2;
#else
      genome[i] = genrand_real1() * scale/init_coef;
      if(i == 0)
	genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
      else genome[i] = genome[i] - (init_coef - 1.0) *scale/2.0/init_coef;
      //genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
#endif //init side
#else
      genome[i] = genrand_real1() * scale/init_coef;
      genome[i] += -scale/2.0/init_coef;
#endif //init_center
#else
      if(i == 0)
	genome[i] = genrand_real1() * scale - scale/2.0;
      else genome[i] = genrand_real1() * scale/2.0;
      //      genome[i] = genrand_real1() * scale/2.0;
#endif 
    }
    if(fitness(genome,size) > MIN_FITNESS)
      true_bit = false;
  }
#endif //mixture bivariate
#else //bivariate
  for(unsigned long i=0;i<size;i++){
      genome[i] = genrand_real1() * scale; 
      if(genome[i] < scale/2)
	  genome[i] = -genome[i];
      else genome[i] -= scale/2;
  }
#endif //bivariate
#else //mixture reals
#ifdef NORMAL_RUNS
  for(unsigned long i=0;i<size;i++){
	  genome[i] = genrand_int32() % 2;}
#else
#ifdef BINOMIAL
     for(unsigned long j = 0; j < (size/(double)BLOCKsize); j++){
	 for(unsigned long i = 0; i<BLOCKsize/2; i++)
	     genome[i + j*BLOCKsize] = 0;
	 for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++)
	     genome[i + j*BLOCKsize] = 1;
     } 
//	  if(chain == 0){
//	      if(i % BLOCKsize >= BLOCKsize/2)
//		  genome[i] = 1;
//	      else genome[i] = 0;
//	  } else
//	      genome[i] = genrand_int32() % 2;   
#else
#ifdef TWO_ATTRACTORS_FUNCTION
  double count1 = 0;
  for(unsigned long i=0;i<size;i++){
	  if(i - count1 >= size/2) {
	      genome[i] = 1;
	      count1++;
	  } else if(count1 <= size/2){
	      genome[i] = rand() % 2;
	      //genome[i] = genrand_int32() % 2;
	      if(genome[i] == 1)
		  count1++;
	  }
	  else genome[i] = 0;  
  }
#else
  for(unsigned long i=0;i<size;i++){
      if(runs == 0) genome[i] = 0;
	  else if(runs == 1) genome[i] = 1;
	  else if(runs == 2) { 
	    if(i < size/2) genome[i] = 0;
	    else genome[i] = 1; 
	  } else if(runs == 3) {
	    if(i < size/3) genome[i] = 0;
	    else genome[i] = 1;
	  } else if(runs == 4) {
	    if(i < size/4) genome[i] = 0;
	    else genome[i] = 1;
	  } else if(runs == 5) {
	    if(i < BLOCKsize) genome[i] = 0;
	    else genome[i] = 1;
	  }
	  else genome[i] = genrand_int32() % 2;
  }
#endif
#endif
#endif 
#endif //mixture reals
    //cout << genome[i] << "||";
  //cout << "\n";	
  value = fitness();
}
#else // no LIST

long binaryMCMC::addElement(){
  if(genome == NULL || nrOfGenerations > sampleSize + 3) {
    cout << "Eroare obiect,  nrOfGenerations==" << nrOfGenerations << " sample size==" << sampleSize << "\n"; 
    return -1;
    }

 nextTemperatures[nrOfGenerations] = temperature;
 for(unsigned long i = 0; i < size; i++)
   nextGENOME[nrOfGenerations][i] = genome[i];
 nextValues[nrOfGenerations] = value;

 return 1;
}

#ifdef MIXTURE_REALS 
 long binaryMCMC::addElement( unsigned long mySize, double* myGenome, double myTemperature, unsigned long mynrOfGenerations){
#else
long binaryMCMC::addElement( unsigned long mySize, unsigned long* myGenome, double myTemperature, unsigned long mynrOfGenerations){
#endif
 if(genome == NULL || sampleSize + 3 < nrOfGenerations) {
    cout << "Eroare obiect ,  nrOfGenerations==" << nrOfGenerations << " sample size==" << sampleSize << "\n"; 
    return -1;
 }
  
 nextTemperatures[nrOfGenerations] = myTemperature;
 for(unsigned long i = 0; i < mySize; i++)
   nextGENOME[nrOfGenerations][i] = myGenome[i];
 nextValues[nrOfGenerations] = fitness(myGenome,mySize);
 return 1;
}

double binaryMCMC::getElementValue(unsigned long i){
  if(genome == NULL || i <= 0 || i > nrOfGenerations + 1 ){
      cout << "Eroare obiect " << index << " because genome" << genome << " or i"  << i << "\n";
      return -1;
    }

  return nextValues[i-1];
}

double binaryMCMC::getElementValue(){
  if(genome == NULL)
    {
      cout << "Eroare obiect " << index << " because genome" << genome << "\n";
      return -1;
    }

  double temp = 0;
  for(unsigned long i = 0; i < nrOfGenerations; i++){
    temp += nextValues[i];
  }

  return temp;
}

double binaryMCMC::getElementValue(unsigned long n, unsigned long burn){
  if(genome == NULL || n > nrOfGenerations + 1 || n < burn)
    {
      cout << "Eroare obiect " << index << " because genome" << genome << "\n";
      return -1;
    }

  double temp = 0;
  for(unsigned long i = burn; i <= n; i++){
    temp += nextValues[i-1];
  }

  return temp;
}

double binaryMCMC::getElementValue(unsigned long n, unsigned long burn, double* teta){
  if(genome == NULL || n > nrOfGenerations + 1 || n < burn){
      cout << "Eroare obiect " << index << " because genome" << genome << "\n";
      return -1;
    }

  double temp = 0;
  for(unsigned long i = burn; i <= n; i++){
     teta[i - 1] = nextValues[i-1];
    temp += teta[i-1];
  }

  return temp;
}

double binaryMCMC::getElementValue(double* teta){
  if(genome == NULL || nrOfGenerations > sampleSize + 2){
      cout << "Eroare obiect nrOfGenerations =" << nrOfGenerations << " because genome" << genome << "\n";
      return -1;
    }

  double temp = 0;
  for(unsigned long i = 0; i <= sampleSize; i++){
    teta[i] = nextValues[i];
    temp += teta[i];
  }

  return temp;
}

double binaryMCMC::getVarianceValue(unsigned long n, unsigned long burn, double mean){
  if(genome == NULL || n > nrOfGenerations + 1 || n < burn)
    {
      cout << "Eroare obiect n=" << n << " because genome" << genome << "\n";
      return -1;
    }

  double temp = 0;
  for(unsigned long i = burn; i <= n; i++)
      temp += pow(nextValues[i-1] - mean,2);

  return temp;
}

double binaryMCMC::getVarianceValue(unsigned long n, unsigned long burn, double mean, unsigned long s){
  if(genome == NULL || n > nrOfGenerations + 1 || n < burn)
    {
      cout << "Eroare obiect n=" << n << " because genome" << genome << "\n";
      return -1;
    }

  double temp = 0;
  for(unsigned long i = burn; i <= n; i++)
    if(s % 2 == 0) temp += pow((double)nextValues[i-1] - mean,(double) s);
    else if(nextValues[i-1] < mean)
      temp += -pow((double)nextValues[i-1] - mean,(double)s);
    else temp += pow((double)nextValues[i-1] - mean,(double)s);

  return temp;
}

binaryMCMC* binaryMCMC::getElement(unsigned long i){
  if(genome == NULL || i <= 0 || i > nrOfGenerations + 1)
    {
      cout << "Eroare obiect i=" << i << "\n";
      return NULL;
    }

    binaryMCMC* temp = new binaryMCMC(size,nextGENOME[i-1],nextTemperatures[i-1],nextValues[i-1]);
    temp -> nrOfGenerations = i-1;
    return temp;
}

void binaryMCMC::reset(unsigned long myRuns){
//    UniformDistribution genrand_real1(0,1);
  index = 0;
  runs = myRuns;
  nrOfGenerations = 0;
  
  //init_genrand(time(0));
  //setSeed();
#ifdef MIXTURE_REALS
#ifdef BIVARIATE
#ifndef MIXTURE_BIVARIATE
  unsigned long integers[size*BLOCKsize];
  for(unsigned long i=0;i<size*BLOCKsize;i++)
	  integers[i] = genrand_int32() % 2;
  convertReal(genome,integers,size);
//  for(unsigned long i=0;i<size;i++){
//	  genome[i] = genrand_real1();
//  }
#else //mixture bivariate
  bool true_bit = true;
  while(true_bit){
    for(unsigned long i=0;i<size;i++){
#ifndef INIT_HALF
#ifndef INIT_CENTER
#ifndef INIT_SIDE
      genome[i] = genrand_real1() * scale; 
      genome[i] += -scale/2;
#else
	genome[i] = genrand_real1() * scale/init_coef;
	if(i == 0)
	  genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
	else genome[i] = genome[i] - (init_coef - 1.0) *scale/2.0/init_coef;
	//	genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
#endif //init side
#else
	genome[i] = genrand_real1() * scale/init_coef;
	genome[i] = genome[i] - scale/2.0/init_coef;
#endif //init_center
#else
      if(i == 0)
	genome[i] = genrand_real1() * scale - scale/2.0;
      else genome[i] = genrand_real1() * scale/2.0;
      //	genome[i] = genrand_real1() * scale/2.0;
#endif 
    }
    if(fitness(genome,size) > MIN_FITNESS)
      true_bit = false;
  }
#endif //mixture bivariate

#else
  bool true_bit = true;
  while(true_bit){
    for(unsigned long i=0;i<size;i++){
#ifndef INIT_HALF
#ifndef INIT_CENTER
#ifndef INIT_SIDE
      genome[i] = genrand_real1() * scale; 
      genome[i] += -scale/2;
#else
      genome[i] = genrand_real1() * scale/init_coef;
      if(i == 0)
	genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
      else genome[i] = genome[i] - (init_coef - 1.0) *scale/2.0/init_coef;
	//      genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
#endif //init side
#else
	genome[i] = genrand_real1() * scale/init_coef;
	genome[i] = genome[i] - scale/2.0/init_coef;
#endif //init_center
#else
      if(i == 0)
	genome[i] = genrand_real1() * scale - scale/2.0;
      else genome[i] = genrand_real1() * scale/2.0;
      //	genome[i] = genrand_real1() * scale/2.0;
#endif 
    }
    if(fitness(genome,size) > MIN_FITNESS)
      true_bit = false;
  }
#endif //bivariate
#else //mixture reals
#ifdef NORMAL_RUNS
  for(unsigned long i=0;i<size;i++){
	  genome[i] = genrand_int32() % 2;
  }
#else
#ifdef BINOMIAL
  if(runs == 1){
      for(unsigned long j = 0; j < (size/(double)BLOCKsize); j++){
	  for(unsigned long i = 0; i<BLOCKsize/2; i++){
	      genome[i + j*BLOCKsize] = 1;
	      //   cout << genome[i + j*BLOCKsize];
	  }
	  for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++){
	      genome[i + j*BLOCKsize] = 0;
	      //cout << genome[i + j*BLOCKsize];
	  }
      }
  } else {
      for(unsigned long i=0;i<size;i++){
	  genome[i] = genrand_int32() % 2;
      }
  }
#else
#ifdef TWO_ATTRACTORS_FUNCTION
  double count1 = 0;
  for(unsigned long i=0;i<size;i++){
	  if(i - count1 >= size/2) {
	      genome[i] = 1;
	      count1++;
	  } else if(count1 <= size/2){
	      genome[i] = rand() % 2;
	      //genome[i] = genrand_int32() % 2;
	      if(genome[i] == 1)
		  count1++;
	  }
	  else genome[i] = 0;  
  }
#else
  for(unsigned long i=0;i<size;i++){
	  if(runs == 0) genome[i] = 0;
	  else if(runs == 1) genome[i] = 1;
	  else 
	      /*if(runs == 2) { 
	    if(i < size/2) genome[i] = 0;
	    else genome[i] = 1; 
	  } else if(runs == 3) {
	    if(i < size/3) genome[i] = 0;
	    else genome[i] = 1;
	  }
	  else*/ genome[i] = genrand_int32() % 2;
  }
#endif
#endif
#endif 
#endif //mixture reals
//	  cout << genome[i];
//   cout << "\n";	
  value = fitness();
}

void binaryMCMC::reset(unsigned long myRuns, unsigned long chain){
    //  UniformDistribution genrand_real1(0,1);
  index = 0;
  runs = myRuns;
  nrOfGenerations = 0;

  //init_genrand(time(0));
  //setSeed();
  double count1 = 0;
#ifdef MIXTURE_REALS
#ifdef BIVARIATE
#ifndef MIXTURE_BIVARIATE
  unsigned long integers[size*BLOCKsize];
  for(unsigned long i=0;i<size*BLOCKsize;i++)
	  integers[i] = genrand_int32() % 2;
  convertReal(genome,integers,size);
//  for(unsigned long i=0;i<size;i++){
//	  genome[i] = genrand_real1();
//  }
#else //mixture bivariate
  bool true_bit = true;
  while(true_bit){
    for(unsigned long i=0;i<size;i++){
#ifndef INIT_HALF
#ifndef INIT_CENTER
#ifndef INIT_SIDE
      genome[i] = genrand_real1() * scale; 
      genome[i] += -scale/2;
#else
	genome[i] = genrand_real1() * scale/init_coef;
	if(i == 0)
	  genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
	else genome[i] = genome[i] - (init_coef - 1.0) *scale/2.0/init_coef;
	//	genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
#endif //init side
#else
	genome[i] = genrand_real1() * scale/init_coef;
	genome[i] = genome[i] - scale/2.0/init_coef;
#endif //init_center
#else
      if(i == 0)
	genome[i] = genrand_real1() * scale - scale/2.0;
      else genome[i] = genrand_real1() * scale/2.0;
      //	genome[i] = genrand_real1() * scale/2.0;
#endif 
    }
    if(fitness(genome,size) > MIN_FITNESS)
      true_bit = false;
  }
#endif //mixture bivariate
#else //bivariate
  bool true_bit = true;
  while(true_bit){
    for(unsigned long i=0;i<size;i++){
#ifndef INIT_HALF
#ifndef INIT_CENTER
#ifndef INIT_SIDE
      genome[i] = genrand_real1() * scale; 
      genome[i] += -scale/2;
#else
      genome[i] = genrand_real1() * scale/init_coef;
      if(i == 0)
	genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
      else genome[i] = genome[i] - (init_coef - 1.0) *scale/2.0/init_coef;
      //genome[i] = genome[i] - (init_coef - 1.0) * scale/2.0/init_coef;
#endif //init side
#else
      genome[i] = genrand_real1() * scale/init_coef;
      genome[i] = genome[i] - scale/2.0/init_coef;
#endif //init_center
#else
      if(i == 0)
	genome[i] = genrand_real1() * scale - scale/2.0;
      else genome[i] = genrand_real1() * scale/2.0;
      //genome[i] = genrand_real1() * scale/2.0;
#endif 
    }
    if(fitness(genome,size) > MIN_FITNESS)
      true_bit = false;
  }
#endif //bivariate
#else
#ifdef NORMAL_RUNS
  for(unsigned long i=0;i<size;i++){
	  genome[i] = genrand_int32() % 2;}
#else
#ifdef BINOMIAL
  if(runs == 1){
    for(unsigned long j = 0; j < (size/(double)BLOCKsize); j++){
	for(unsigned long i = 0; i<BLOCKsize/2; i++){
	    genome[i + j*BLOCKsize] = 1;
	    //    cout << genome[i + j*BLOCKsize];
	}
	for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++){
	    genome[i + j*BLOCKsize] = 0;
	    //cout << genome[i + j*BLOCKsize];
	}
    }
  } else {
      for(unsigned long i=0;i<size;i++)
	  genome[i] = genrand_int32() % 2;  
  }
//	  if(chain == 0){
//	      if(i % BLOCKsize >= BLOCKsize/2)
//		  genome[i] = 1;
//	      else genome[i] = 0;
//	  } else
//	      genome[i] = genrand_int32() % 2;   
#else
#ifdef TWO_ATTRACTORS_FUNCTION
  for(unsigned long i=0;i<size;i++){
	  if(i - count1 >= size/2) {
	      genome[i] = 1;
	      count1++;
	  } else if(count1 <= size/2){
	      //genome[i] = drand48() * 2;
	      genome[i] = genrand_int32() % 2;
	      if(genome[i] == 1)
		  count1++;
	  }
	  else genome[i] = 0;  
  }
#else
  for(unsigned long i=0;i<size;i++){
	  if(runs == 0) genome[i] = 0;
	  else if(runs == 1) genome[i] = 1;
	  else /*if(runs == 2) { 
	    if(i < size/2) genome[i] = 0;
	    else genome[i] = 1; 
	  } else if(runs == 3) {
	    if(i < size/3) genome[i] = 0;
	    else genome[i] = 1;
	  }
	  else */ genome[i] = genrand_int32() % 2;
  }
#endif
#endif
#endif
#endif //mixture reals 
//	  cout << genome[i];
//  cout << "\n";	
  value = fitness();
}
#endif //LIST


#ifdef TRAP_LIKE_FUNCTION
void binaryMCMC::print(){
        unsigned long correct = 0;
	unsigned long contor = 0;
	unsigned long wrong = 0;
	unsigned long contorW = 0;

	for(unsigned long i = 0; i<size; i++) 
	  if(i % BLOCKsize == BLOCKsize - 1) { 
	    cout << genome[i] << "-";
	    if(genome[i] == 1) contor++;
	    else contorW++;
	    if (contor == BLOCKsize) correct ++;
	    if (contorW == BLOCKsize) wrong++;
	    contor = 0;
	    contorW = 0;
	  } 
	  else {
	    cout << genome[i] << "|";
	    if(genome[i] == 1) contor++;
	    contorW++;
	  }
	cout << "\tF=" <<value;
	cout << "\tT=" <<temperature;
	cout << "\tCB=" << correct;
	cout << "\tWB=" << wrong <<"\n";
	//if(next != NULL) next -> print();
}

void binaryMCMC::print(ofstream& myFile){
        unsigned long correct = 0;
	unsigned long wrong = 0;
	unsigned long contor = 0;
	unsigned long contorW = 0;
	for(unsigned long i = 0; i<size; i++) 
	  if(i % BLOCKsize == BLOCKsize - 1) { 
	    myFile << genome[i] << "-";
	    if(genome[i] == 1) contor++;
	    else contorW ++;
	    if (contor == BLOCKsize) correct ++;
	    if (contorW == BLOCKsize) wrong ++;
	    contor = 0;
	    contorW = 0;
	  } 
	  else {
	    myFile << genome[i] << "|";
	    if(genome[i] == 1) contor++;
	    else contorW++;
	  }
	myFile << "\tF=" <<value;
	myFile << "\tT=" <<temperature;
	myFile << "\tCB=" << correct; 
	myFile << "\tWB=" <<wrong<<"\n";	
#ifdef SAMPLE
#ifndef MULTIPLE_RUN_SAMPLE
	if(index+1 > sampleSize)
	  if(next != NULL) next -> printLEVEL(index+1,myFile);
#endif
#endif
}

void binaryMCMC::printLEVEL(unsigned long level, ofstream& myFile){
  unsigned long correct = 0;
  unsigned long wrong = 0;
  unsigned long contor = 0;
  unsigned long contorW = 0;
#ifdef LIST
  if(level != index && next != NULL) {next -> printLEVEL(level,myFile);return;}
  if(level != index && next == NULL) {
    cout << "ERROR in finding the print level " << level << "\n";
    return;
  }
#endif
  for(unsigned long i = 0; i<size; i++) 
    if(i % BLOCKsize == BLOCKsize - 1) { 
      myFile << genome[i] << "-";
      if(genome[i] == 1) contor++;
      else contorW ++;
      if (contor == BLOCKsize) correct ++;
      if (contorW == BLOCKsize) wrong ++;
      contor = 0;
      contorW = 0;
    } 
    else {
      myFile << genome[i] << "|";
      if(genome[i] == 1) contor++;
      else contorW++;
    }
  myFile << "\tF=" <<value;
  myFile << "\tT=" <<temperature;
  myFile << "\tCB=" << correct; 
  myFile << "\tWB=" <<wrong<<"\n";	
}
 
#else 
void binaryMCMC::printLEVEL(unsigned long level, ofstream& myFile){
  unsigned long correct = 0;
  unsigned long wrong = 0;
#ifdef LIST
  if(level != index && next != NULL) {next -> printLEVEL(level,myFile);return;}
  if(level != index && next == NULL) {
    cout << "ERROR in finding the print level " << level << "\n";
    return;
  }
#endif
  for(unsigned long i = 0; i<size; i++) { 
    myFile << genome[i];
  } 
  myFile << "\tF=" <<value;
  myFile << "\tT=" <<temperature;
  myFile << "\tCB=" << correct; 
  myFile << "\tWB=" <<wrong<<"\n";	
}

void binaryMCMC::print(){
  for(unsigned long i = 0; i<size; i++) 
    cout << genome[i] << " | "; 
  cout << "fitness== " <<value;
  cout << "temperature" <<temperature<<"\n";	
#ifdef LIST
  if(next != NULL) next -> print();
#endif
}

void binaryMCMC::print(ofstream& myFile){
  for(unsigned long i = 0; i<size; i++) 
    myFile << genome[i] << " | "; 
  myFile <<" ==" <<value << " || T = " <<temperature;
  myFile <<"|| generations == " <<nrOfGenerations<< '\n'; 
#ifdef LIST
  if(next!=NULL) next -> print(myFile);
#endif
}

#endif

void binaryMCMC::printState(){
	for(unsigned long i = 0; i<size; i++) 
		cout << genome[i] << " | "; 
	cout << "fitnessBLOCKsize== " <<value <<"|| generations == "<<nrOfGenerations <<"\n";	
}



/*void binaryMCMC::setMutation(double myMutation)
{
	mutation = myMutation * mutation_multi;
}*/

void binaryMCMC::setRecombination(double myRecombination)
{
  recombination = myRecombination;
}

void binaryMCMC::setRandomGenome(){
	for(unsigned long i= 0; i<size; i++)
		genome[i] = genrand_int32() % 2;
}

void binaryMCMC::setSolutionGenome1(){
	for(unsigned long i= 0; i<size; i++)
		genome[i] = 1;
}

void binaryMCMC::setSolutionGenome0(){
	for(unsigned long i= 0; i<size; i++)
		genome[i] = 0;
}

double binaryMCMC::getValue(){
	return value;
}

double binaryMCMC::getValue(unsigned long level){
    unsigned long i = 0;
    binaryMCMC* top = this;
#ifdef LIST
    while(top != NULL && i < level){
	i++;
	top = top->next;
    }
    if(i!= level) return -1;
    else retun top->value;
#else
	return nextValues[level];
#endif //LIST
}

double binaryMCMC::getTEMPERATURE(){
	return temperature;
}

void binaryMCMC::setTEMPERATURE(double mytemp){
	temperature = mytemp;
}

#ifndef MIXTURE_REALS
unsigned long* binaryMCMC::getGENOME(unsigned long* myGenome)
{
	for(unsigned long i = 0; i<size; i++)
		myGenome[i] = genome[i];
	return myGenome;
}

unsigned long* binaryMCMC::getGENOME(unsigned long i, unsigned long* myGenome)
{
  if(genome == NULL || i <= 0 || i > nrOfGenerations + 1)
    {
      cout << "Eroare obiect i=" << i << "\n";
      return NULL;
    }
#ifndef LIST
  for(unsigned long j = 0; j<size; j++)
      myGenome[j] = nextGENOME[i][j];
#else
#endif //LIST  
  return myGenome;
}

void binaryMCMC::setGENOME(unsigned long* myGenome)
{
	for(unsigned long i = 0; i<size; i++)
		genome[i] = myGenome[i];
	fitness();
}
#else
double* binaryMCMC::getGENOME(double* myGenome)
{
	for(unsigned long i = 0; i<size; i++)
		myGenome[i] = genome[i];
	return myGenome;
}

double* binaryMCMC::getGENOME(unsigned long i, double* myGenome)
{
  if(genome == NULL || i <= 0 || i > nrOfGenerations + 1)
    {
      cout << "Eroare obiect i=" << i << "\n";
      return NULL;
    }
#ifndef LIST
  for(unsigned long j = 0; j<size; j++)
      myGenome[j] = nextGENOME[i][j];
#else
#endif //LIST  
  return myGenome;
}

void binaryMCMC::setGENOME(double* myGenome)
{
	for(unsigned long i = 0; i<size; i++)
		genome[i] = myGenome[i];
	fitness();
}
#endif //MIXTURE_REALS

unsigned long binaryMCMC::getACCEPTED(){
    if(accepted != 1 && accepted != 0) {
	cout << "wrong acceptance=" <<accepted<<" in binaryMCMC.cpp:1653\n";
	exit(1);
    }
	return accepted;
}

void binaryMCMC::setACCEPTED(unsigned long myacc){
	accepted = myacc;
}

#ifndef METROPOLIS_ALGORITHM
void binaryMCMC::setProposalDistribution(double myProposal){
  proposalDistribution = myProposal;
}

double binaryMCMC::getProposalDistribution(){
  return proposalDistribution;
}

void binaryMCMC::setProposalDistributionVector(double** tempDistribution){
  /*if(proposalDistributionVector == NULL){
    proposalDistributionVector = new double*[size];
    for(long i=0; i < size; i++){
      proposalDistributionVector[i] = new double[2];
      for(long j = 0; j < 2; j++)
	proposalDistributionVector[i][j] = 1.0;
    }
  }

  for(long i=0; i < size; i++){
    for(long j = 0; j < 2; j++)
      cout << proposalDistributionVector[i][j] << " " << tempDistribution[i][j] << " " << i << " " << j<<"\n";
      } */  
  for(unsigned long i = 0; i < size; i++){
    for(unsigned long j = 0; j < 2; j++){
      proposalDistributionVector[i][j] = tempDistribution[i][j];
    }}
}

double** binaryMCMC::getProposalDistributionVector(double** tempDistribution){
  for(unsigned long i = 0; i < size; i++)
    for(unsigned long j = 0; j < 2; j++)
      tempDistribution[i][j] = proposalDistributionVector[i][j];
  return tempDistribution;
}
#endif

unsigned long binaryMCMC::getNrRuns(){
	return runs;
}

void binaryMCMC::setNrRuns(unsigned long myRuns){
       runs = myRuns;
}

void binaryMCMC::generation_Update(){
  nrOfGenerations++;
}

void binaryMCMC::runs_Update(){
  runs++;
}

unsigned long binaryMCMC::equal(binaryMCMC* string2){
#ifndef MIXTURE_REALS
  unsigned long* genome2 = new unsigned long[size];
#else
  double* genome2 = new double[size];
#endif //MIXTURE_REALS
  
  genome2 = string2->getGENOME(genome2);
  for(unsigned long i = 0; i < size; i++)
    if(genome[i] != genome2[i]){
      delete genome2;
      return 0;
    }
  delete genome2;
  
  genome2 = string2->getGENOME(genome2);
  for(unsigned long i = 0; i < size; i++)
    if(genome[i] != genome2[i]){
      delete genome2;
      return 0;
    }
  delete genome2;
  return 1;
}

#ifndef MIXTURE_REALS
unsigned long binaryMCMC::getAllele(unsigned long i){
#else
double binaryMCMC::getAllele(unsigned long i){
#endif //MIXTURE_REALS
    if(i >= 0 && i < size)
	return genome[i];
    cout << "Error, no such allele " << i << " in binaryMCMCtrap.cpp:2241 \n";
    return -1;
}

unsigned long binaryMCMC::notAllele(unsigned long i){
    if(i ==0) return 1;
    else return 0;
}

#ifdef MIXTURE_REALS
//read file "parametri_mixture.file" 
void binaryMCMC::read_parametri(){
//parametri de la mixture
  char buffer[500];

  ifstream myFile("parametri_mixture.file");		
  if(!myFile.is_open()){
    cout << "Error opening file \n";
    exit(1);
  }
 
//type of lanscape
  myFile >> buffer;
  typeLandscape = atoi(buffer);
  cout << buffer << " should be" << typeLandscape << "\n";
  myFile.getline(buffer,500);
  //cout << buffer << "\n";
  cout << " discard " << buffer << "\n";

//number of components
//primel cimp este explicative
  myFile.getline(buffer,500);
  cout << " discard " << buffer << "\n";
  myFile.getline(buffer,500);
  mixture = atoi(buffer);
  cout << buffer << " should be " << mixture << "\n";

//number of dimensions
  myFile.getline(buffer,500);
  cout << " discard " << buffer << "\n";
  myFile.getline(buffer,500);
  cout << buffer << "should be" << atoi(buffer) << "\n";
  //mixture = atoi(buffer);
  
//the bounds of the landscape
  myFile.getline(buffer,500);
  cout << " discard " << buffer << "\n";
  myFile.getline(buffer,500);
  cout << " discard " << buffer << "should be " <<  atoi(buffer);
  //initialization of the position and covariance matrices

  //cout << "mixture " << mixture << "\n";
  positionMixture = new Matrix[mixture];
  determinantMixture = new double[mixture]; 
  covarianceMixture = new Matrix[mixture];
  covarianceMixtureInverse = new Matrix[mixture];
  for(int i = 0; i < mixture; i++){
    positionMixture[i].ReSize(size,1);
    for(int j = 0; j < size; j++)
      positionMixture[i].element(j,0) = 0;
    covarianceMixture[i].ReSize(size,size);
    covarianceMixtureInverse[i].ReSize(size,size);
  }

//if we have correation in a matrix
  if(typeLandscape == 1 || typeLandscape == 2) {
    //position is always in the middle of the space 0
    //cout << " enter correctly on typeLandscape " << typeLandscape << " position " << positionMixture[0] << "\n";
    for(int i = 0; i < size; i++){
      //cout << "position " << i << "\n";
      positionMixture[0].element(i,0) = 0;
    }
    
    //read the title of the covariance matrix
    myFile.getline(buffer,500);
    //cout << "discard last before matrix" << buffer << "\n";

      //read the non-correlated covariance matrix
      for(int i = 0; i < size; i++){
	  for(int j = 0; j < size; j++){
	      myFile >> buffer;
	      if(typeLandscape == 2){
		covarianceMixture[0].element(i,j) = atof(buffer);
		//	cout << "(" << i << "," << j << "," << covarianceMixture[0].element(i,j) << ")\n";
	      } //else 
		//cout << "(" << i << "," << j << "," << buffer << ")\n";
	  }
	  //cout << "\n";
      }

      //read the angles title and angles
      myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";
      myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";

      //read the title of the corelated covariance matrix
      myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";
      myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";

      for(int i = 0; i < size; i++)
	  for(int j = 0; j < size; j++){
	      myFile >> buffer;
	      if(typeLandscape == 1){
		  covarianceMixture[0].element(i,j) = atof(buffer);
		  //	  cout << "(" << i << "," << j << "," << covarianceMixture[0].element(i,j) << ")\n";
	      } //else 
		//cout << "(" << i << "," << j << "," << buffer << ")\n";
	  }

      //write the inverse and the determinant here
      covarianceMixtureInverse[0] = covarianceMixture[0].i();

      determinantMixture[0] = covarianceMixture[0].Determinant();
      //cout << " type lanscape " << typeLandscape << " is " << covarianceMixture[0] << "\n";
  } else if(typeLandscape == 3){
     //position is always in the middle of the space 0
    //cout << " enter correctly on typeLandscape " << typeLandscape << " position " << positionMixture[0] << "\n";
    for(int i = 0; i < size; i++){
      //cout << "position " << i << "\n";
      positionMixture[0].element(i,0) = 0;
    }
    
    //read the title of the covariance matrix
    myFile.getline(buffer,500);
    //cout << "discard last before correlation matrix" << buffer << "\n";

      //read the correlation matrix
      for(int i = 0; i < size; i++){
	  for(int j = 0; j < size; j++){
	      myFile >> buffer;
	      //  cout << "(" << i << "," << j << "," << buffer << ")\t";
	      //if(typeLandscape == 2){
	      //	covarianceMixture[0].element(i,j) = atof(buffer);
		//	cout << "(" << i << "," << j << "," << covarianceMixture[0].element(i,j) << ")\n";
	      //    } //else 
	  	//cout << "(" << i << "," << j << "," << buffer << ")\n";
	  }
	  //cout << "\n";
      }

      //read the angles title and angles
      //myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";
      //myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";

      //read the title of the corelated covariance matrix
      myFile.getline(buffer,500);
      //cout << " discard last buffer before covariance matrix" << buffer << "\n";
      myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";

      for(int i = 0; i < size; i++){
	  for(int j = 0; j < size; j++){
	      myFile >> buffer;
	      //if(typeLandscape == 1){
	      covarianceMixture[0].element(i,j) = atof(buffer);
	      //cout << "(" << i << "," << j << "," << covarianceMixture[0].element(i,j) << ")\t";
	      //} //else 
		//cout << "(" << i << "," << j << "," << buffer << ")\n";
	  }
	  //cout << "\n";
      }
      //write the inverse and the determinant here
      covarianceMixtureInverse[0] = covarianceMixture[0].i();

      determinantMixture[0] = covarianceMixture[0].Determinant();
      //cout << " type lanscape " << typeLandscape << " is " << covarianceMixture[0] << "\n"; 
  } else {
    // extract which size
    cout << mixture << "\n";
    for(int i = 0; i < mixture; i++){
      myFile >> buffer;
      myFile.getline(buffer,100); 
      //cout << " numerical position[" << i << "]= "  << buffer << "\t";
      for(int k = 0; k < size; k++){
	myFile >> buffer;
	positionMixture[i].element(k,0) = atof(buffer);
	//cout << buffer << " , ";
      }
      //cout << "\n";
      
      //read the title of the covariance matrix
      myFile.getline(buffer,500);
      //cout << "discard last before correlation matrix" << buffer << "\n";
      myFile >> buffer;
      //read the correlation matrix
      for(int i1 = 0; i1 < size; i1++){
	for(int j = 0; j < size; j++){
	  myFile >> buffer;
	  //  cout << "(" << i1 << "," << j << "," << buffer << ")\t";
	  //if(typeLandscape == 2){
	  //	covarianceMixture[0].element(i,j) = atof(buffer);
	  //	cout << "(" << i << "," << j << "," << covarianceMixture[0].element(i,j) << ")\n";
	  //    } //else 
	  //cout << "(" << i << "," << j << "," << buffer << ")\n";
	}
	//cout << "\n";
      }
      
      //cout << "\n";
      //read the angles title and angles
      //myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";
      //myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";

      //read the title of the corelated covariance matrix
      myFile.getline(buffer,500);
      //cout << " discard last buffer before covariance matrix" << buffer << "\n";
      myFile.getline(buffer,500);
      //cout << " discard " << buffer << "\n";

      for(int i1 = 0; i1 < size; i1++){
	  for(int j = 0; j < size; j++){
	      myFile >> buffer;
	      //if(typeLandscape == 1){
	      covarianceMixture[i].element(i1,j) = atof(buffer);
	      //cout << "(" << i1 << "," << j << "," << covarianceMixture[i].element(i1,j) << ")\t";
	      //} //else 
		//cout << "(" << i << "," << j << "," << buffer << ")\n";
	  }
	  //cout << "\n";
      }
      //cout << "\n";
      //write the inverse and the determinant here
      covarianceMixtureInverse[i] = covarianceMixture[i].i();

      determinantMixture[i] = covarianceMixture[i].Determinant();
      //cout << " type lanscape " << typeLandscape << " is " << covarianceMixture[0] << "\n"; 
    }
  }

  //---
  // cout << "covariance \n" << covarianceMixture[0] << "position \n" << positionMixture[0] << " covariance Inverse \n" << covarianceMixtureInverse[0] << "\n"; 
  //cout << " determinant " << determinantMixture[0] << "\n";
  myFile.close();
}

void binaryMCMC::delete_parametri(){
  delete[] positionMixture;
  delete[] covarianceMixture;
  delete[] covarianceMixtureInverse;
  delete[] determinantMixture;
}
#endif //mixture reals
