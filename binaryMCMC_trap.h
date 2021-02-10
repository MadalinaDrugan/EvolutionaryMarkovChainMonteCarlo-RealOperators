#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <math.h>

#include "MatrixUtil\newmatap.h"                // need matrix applications
#include "MatrixUtil\newmatio.h"                // need matrix output routines

#include "PARAMETERS.h"
//#include "gaussian.cpp"

class binaryMCMC{
    //list of samples
    // static UniformDistribution* genrand_real1;

 public:
#ifdef LIST
    binaryMCMC* next;
#else

#ifndef MIXTURE_REALS
  unsigned long** nextGENOME;
#else

  double** nextGENOME;
#ifdef wPCTX
  double mix_mixture;
#endif //wPCTX
#endif //MIXTURE_REALS

  double* nextTemperatures;
  double* nextValues;
#endif

  unsigned long index;
  unsigned long runs;
  
  //GA part
#ifndef MIXTURE_REALS
  unsigned long* genome;
#else
  double* genome;
  //int* genomeREALS;
#endif

  unsigned long nrOfGenerations; // the number of generation in the population  
  unsigned long size; //dimensiunea genotipului
  
  //MCMC and GA
  double value; //fitness of genome
  double mutation;
//#ifdef REALS
//  double mutation_sigma;
//#endif
  double recombination;

#ifdef MULTIPLE_RUNS
#ifndef MUTATION_VARIES
#ifdef IsRecombination
	unsigned long PROB_nPOINT;
	double UNIF_RECOMB_PARAM;
#endif //IsRecombination
#endif
#endif //MULTIPLE_RUNS
 
  //MCMC part
  unsigned long accepted; //indicator if it is accepted or not in the next population; defalt = 0;
  double temperature; //tempreatura MCMC

#ifndef METROPOLIS_ALGORITHM
  double proposalDistribution;// proposal distribution for MCMC
  double** proposalDistributionVector;
#ifndef POPULATION_MCMC
  unsigned long parentReplaced;
#endif
#endif //Metropolis_algorithm

  /*#ifdef BINOMIAL
 static const int binomial_table[2 * BLOCKsize - 1] = {1,4,6,4,1};
 #endif*/

#ifndef MIXTURE_REALS
  binaryMCMC(unsigned long,double);
  binaryMCMC(unsigned long,double,unsigned long);
  binaryMCMC(unsigned long,double,unsigned long);
  binaryMCMC(unsigned long,unsigned long*,double);
  binaryMCMC(unsigned long,unsigned long*,double,double);
#ifdef SINGLE_vs_MULTIPLE
  binaryMCMC(unsigned long,double,unsigned long,unsigned long);
#endif
#else
  binaryMCMC(unsigned long,double*,double,double);
  binaryMCMC(unsigned long,double);

  double n(double,double,double);
  double fitness_mixture(double);
  double fitness_gauss(double*);
  double fitness_bivariate(double*,double, double, double, double, double);
#endif

  ~binaryMCMC();
 
  unsigned long sampleSize;
 
  double fitness();
  unsigned long stop();

#ifndef MIXTURE_REALS
  double fitness(unsigned long*,unsigned long);
  unsigned long stop(unsigned long*,unsigned long);
#else
  double fitness(double*,unsigned long);
  double fitness(unsigned long*,unsigned long);
  unsigned long stop(double*,unsigned long);
#endif

#ifndef MIXTURE_REALS 
  unsigned long* MutationGenome(unsigned long*);
#else
  double* MutationGenome(double*);
  double* MutationOtherGenome(double*);
#endif

#ifdef POPULATION_MCMC
  unsigned long* MutationGenome(unsigned long*,double**,unsigned long);
  unsigned long randJ;
#endif

#ifndef METROPOLIS_ALGORITHM
  double ProposalDistribution(unsigned long*,unsigned long,double**);
  double ProposalDistribution(double**);  
  double ProposalDistribution();

  void AddaptationMethods(double**);
#endif //Metropolis_algorithm

#ifdef BINARY_RECOMB
  unsigned long** RecombinationGenome(unsigned long**,unsigned long);
#ifndef METROPOLIS_ALGORITHM
  unsigned long** RecombinationGenome(unsigned long**,unsigned long,double**);
#endif //Metropolis_algorithm
#else
  double** RecombinationGenome(double**,unsigned long);
#endif // 

#ifdef MIXTURE_REALS
  unsigned long acceptance(double*,double);
#else
#ifdef METROPOLIS_ALGORITHM
  unsigned long acceptance(unsigned long*,double);
#else
  unsigned long acceptance(unsigned long*,double,double);
#endif
#endif //MIXTURE_REALS

  unsigned long nextGeneration();		
  unsigned long nextGenerationBlank();
#ifndef MIXTURE_REALS		
  unsigned long nextGeneration(unsigned long[],double,unsigned long);
#else
  unsigned long nextGeneration(double[],double,unsigned long);
#endif
  void generation_Update();
  void runs_Update();
  void program(char*);
	
  void reset(unsigned long);
  void reset(unsigned long,unsigned long);
  void print();
  void printState();
  void print(ofstream&);
  void printLEVEL(unsigned long,ofstream&);
	
  long addElement();

#ifndef MIXTURE_REALS
  long addElement(unsigned long, unsigned long[], double, unsigned long);
#else
  long addElement(unsigned long, double[], double, unsigned long);
#endif

  double getElementValue(unsigned long);
  double getElementValue(unsigned long,unsigned long);
  double getElementValue(unsigned long,unsigned long, double*);
  double getElementValue(double*);
  double getVarianceValue(unsigned long,unsigned long, double);
  double getVarianceValue(unsigned long,unsigned long, double, unsigned long);
  double getElementValue();
  binaryMCMC* getElement(unsigned long);

  void setRandomGenome();
  void setSolutionGenome1();
  void setSolutionGenome0();
  void setMutation(double);
  void setRecombination(double);
  double getValue();
  double getValue(unsigned long);
  double getTEMPERATURE();
  void setTEMPERATURE(double);

#ifndef MIXTURE_REALS
  unsigned long* getGENOME(unsigned long*);
  unsigned long* getGENOME(unsigned long,unsigned long*);
  void setGENOME(unsigned long*);
#else
  double* getGENOME(double*);
  double* getGENOME(unsigned long,double*);
  void setGENOME(double*);
#endif //MIXTURE_REALS

  void setACCEPTED(unsigned long);
  unsigned long getACCEPTED();
#ifndef METROPOLIS_ALGORITHM
  double getProposalDistribution();
  void setProposalDistribution(double);
  void setProposalDistributionVector(double**);
  double** getProposalDistributionVector(double**); 
#endif //Metropolis_algorithm

  void setNrRuns(unsigned long);
  unsigned long getNrRuns();

  unsigned long equal(binaryMCMC*);

#ifndef MIXTURE_REALS
  void getMax(unsigned long*);
  void getMax0(unsigned long*);
#else
  void getMax0(double*);
  void getMax(double*);
  double getMaxValue();
  void getMax1(double*);
  /*#ifdef TWO_PARENTS

#ifdef  TRANSLATION
  void translation(double**,double*,double*);
  void translationPCTX(double**,double*,double*,double);
  void translationDiff(double**,double*,double*);
  void translationSnooker(double**, double*);
#endif 

#ifdef ROTATION
  Matrix CartesianToPolar(Matrix&, Matrix);
  Matrix PolarToCartesian(Matrix&, Matrix);
  Matrix rotationPolar(Matrix, Matrix);
  void rotation(double**,double*,double*);
#ifdef SINGLE_INDIV_ROT
  int mix_rotation; 
#endif 
  //void rotationMatrix(Matrix*, double);
#endif 

#ifdef SCALING
  void scaling_dim(double**,double*,double*);
#endif
#else*/ //two parents
#ifdef  TRANSLATION
  void translation(double**,double*,double*);
  void translationDiff(double**,double*,double*);
  void translationPCTX(double**,double*,double*,double);
  void translationSnooker(double**, double*);
#ifdef SIMPLEX
  void translationSimplex(double**, double*, double*);
  double fac; 
  double* proportions;
  unsigned long positionChanged;
  unsigned long rank;
#endif //simplex
#endif 

  //#ifdef ROTATION
  void rotation(double**,double*,double*);
  void rotation(double**,double*,double*,unsigned long);
  Matrix CartesianToPolar(Matrix&, Matrix);
  Matrix PolarToCartesian(Matrix&, Matrix);
  Matrix rotationPolar(Matrix&, Matrix);
  //void rotationMatrix(Matrix*, double);
#ifdef SINGLE_INDIV_ROT
  int mix_rotation; 
#endif 
  //#endif 

  //#ifdef SCALING
  void scaling_dim(double**,double*,double*);
  void scaling_dim(double**,double*,double*,unsigned long);
  //#endif
  //#endif //two parents
  
  double bounds(double);

#endif //mixture reals

#ifdef RECOMB_MUT
  void setMax0();
  void setMax1();
  void setMax2();
#endif

#ifdef QBP
  double** matrixQBP;
  void readQBP(unsigned long,double**);
  binaryMCMC(unsigned long,double,double**);
#endif

#ifndef MIXTURE_REALS
  unsigned long getAllele(unsigned long);
#else
  double getAllele(unsigned long);

  double fitness_multivariate(double*);

  Matrix* positionMixture;
  Matrix* covarianceMixture;
  double* determinantMixture;
  Matrix* covarianceMixtureInverse;

  void read_parametri();
  void delete_parametri();

  int typeLandscape;
  int mixture;
#endif

  unsigned long notAllele(unsigned long);

};





