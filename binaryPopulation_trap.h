//fileName binaryPopulation.cpp
//Version 0.4
//Date 13.11.2002
//Description - population of EMCMC with different Metropolis algorithms
//	- one individual - one MCMC - one binary chromozon
//	- each individual has a differenct temperatureaccording with its fitness
//	- best and worst fitness of individuals are stored
// news
//      - implementation for elitist coupling wit regeneration

#include <iostream.h>

#include "Utils\SIMUL.h"
#include "HashtablesUtil\Hashtable.h"
#include "HashtablesUtil\Hashtable_histogram.h"

//#include "random/randoma.h"
#include "binaryMCMC_trap.h"

class binaryPopulation{//: public population{
//generation of uniform numbers
    // U01_Distribution genrand_real1;
//generation of normal numbers for the component functions
/*#ifdef MIXTURE_REALS
    static NormalDistribution* n1;
    static NormalDistribution* n2;
    static NormalDistribution* n3;
    static NormalDistribution* n4;
//generation of normal functions (0,1) for the alleator number generator
    static NormalDistribution* nn; 
#endif //mixture reals
*/

#ifdef MIXTURE
  int mix_mutation;
#endif

public:
  binaryMCMC** MCMCstring;
  binaryMCMC** proposedMCMCstring;
   
  char* fileData;
  unsigned long generation;
  unsigned long sizeMCMCstring;
  unsigned long sizeGENOME;
  
  int burn_in;
  int sampleSize;  

  //new
#ifndef SINGLE_vs_MULTIPLE
  SIMULATIONS *simul;
#else
  SIMULATIONS **simul;
#endif

  double averageGeneration();
  double geometricGeneration();
  double averageFitness;
  double bestFitnessP();

//#ifdef MIXTURE_REALS
//  CDataDistrib* gaussian;
//#endif

#ifdef QBP
  double** matrixQBP;
  fstream myQBP;
  void generateQBP(unsigned long,double**);
  void readQBP(unsigned long,double**);
#endif

#ifdef SINGLE_vs_MULTIPLE
#ifdef TWO_ATTRACTORS_FUNCTION
  double* meanFitness;
  double* diffFitness;
#else
  double meanFitness;
  double diffFitness;
#endif
#else
  double meanFitness;
  double diffFitness;
#endif
  double mean;
  double total;

#ifndef SINGLE_vs_MULTIPLE
#ifdef TWO_ATTRACTORS_FUNCTION
  double** bestFitnessVector;
#else
  double* bestFitnessVector;
#endif 
  double* secondBestFitnessVector;
#else //no single vs. 2 chains
#ifdef TWO_ATTRACTORS_FUNCTION
  double*** bestFitnessVector;
#else
  double** bestFitnessVector;
#endif
#endif

  double bestFitness;
  double worstFitness;

  double bestTemperature;
  double worstTemperature;

  // best of worst of MCMC are updated
  binaryMCMC* bestMCMC;
  binaryMCMC* worstMCMC;

#ifndef METROPOLIS_ALGORITHM
#ifdef POPULATION_MCMC
  double ** PropVector;
  unsigned long randI;
  unsigned long randJ;
#else
  double* proposalDistribution;
  double* proposedProposalDistribution;
#endif
#endif


#ifdef ELITIST_ACCEPTANCE_WITH_REGENERATION
  //temporal data
  double** tempResults;

#ifdef COMPLICAT_TEMPERING
  Hashtable* myHashtables2;
#endif
#endif //ELITIST_ACCEPTANCE_WITH_REGENERATION

  //private functions
  void assignTemperature();
#ifndef MIXTURE_REALS 
  Hashtable* myHashtables;
#else 
#ifdef HISTOGRAM
  Hashtable* myHashtables;  
#else //histogram
  Hashtable_histogram** myHashtables;
#endif //histogram  
#endif //mixture reals

  unsigned long *burn_inT;
  unsigned long interSamplesize;

#ifdef GELMAN_RUBIN

  double **teta;
  double **s_2;

#ifdef NORMAL_DISRIBUTION
  void varianceCalculation();
  double varianceCalculation(unsigned long);
  double meanCalculation(unsigned long);
  void printReductionGNUFILE();
#endif
  void varianceHigerMomentsUnivariateCalculation(unsigned long, unsigned long);
  double varianceCalculation(unsigned long,unsigned long);
  double meanCalculation(unsigned long,unsigned long);
  void printReductionGNUFILE(unsigned long);

  void printReductionMultivariateGNUFILE();
  double ***tetaM;
  void varianceMultivariateCalculation();
  double varianceMultivariateCalculation(unsigned long, double**);
  double meanMultivariateCalculation(unsigned long, double**);
#endif //GELMAN_RUBIN

#ifdef FILE_BOOKEEPING
    fstream diversityFile;
#else
    ofstream diversityFile;
    void printDiversity();
    void printDiversityPopulation();
#endif

    unsigned long runs;

  double f(double);

#ifdef ACCEPTANCE_RATIO

  void printAcceptance(ofstream&);
  double acceptanceAverage();

#ifdef SINGLE_vs_MULTIPLE
#ifdef TWO_ATTRACTORS_FUNCTION
  double* meanAcceptance;
#else
  double meanAcceptance;
#endif
  double** acceptance_ratio;
#else
  double* acceptance_ratio;
  double meanAcceptance;
#endif
#endif

#ifdef HISTOGRAM
#ifdef KULLBACK_INFORMATION
#ifdef TIME_PERFORMANCE
  double** timeKLdistance;
#endif
#ifdef EXPAND_KULLBACK_INFORMATION
    void expandTrueUnivariateDistribution();
    void expandTrueMultivariateDistribution();
#endif
    double **Kullback;
    double mixingTimeUnivariateCalculation();
#ifndef HISTOGRAM
    double *trueDistribution;
#else
    double *trueDistribution;
    double *diversityTrueDistribution;
#endif
    double printMixingGNUFILE();
#ifdef MULTIPLE_RUNS
    void resetKullback();
#endif
#endif
#endif //histogram

#ifdef CUSUM
#ifndef SINGLE_vs_MULTIPLE
#ifndef multiple_restarts
    double* Cusum;
    double* Hairness;
#else
    double **Cusum;
    double **Hairness;
#endif //nr_runs
#else
    double ***Cusum;
    double ***Hairness;
#endif
    double CusumCalculation();
    double printCusumGNUFILE();
#ifdef MULTIPLE_RUNS
    void resetCusum();
#endif
#endif

//performance measure that referes to diversity
#ifdef HISTOGRAM
    unsigned long diversityMeasure;
#endif //histogram

//total number of found solutions 
    double procentRightSolutions;
    double weightRightSolutions;
#ifndef SINGLE_vs_MULTIPLE
    double* globalRightSolutions;
#else
    double** globalRightSolutions;
#endif

//diversity in each run
#ifndef SINGLE_vs_MULTIPLE
#ifndef MIXTURE_REALS
    double** globalDiversity;
    double** globalDiversityRightSolutions;
#else
#ifdef HISTOGRAM
    double*** globalDiversity; //number of individuals in a bin
    double*** globalDiversityRightSolutions;//weight of the individuals in a bin
#endif //histogram
#endif
#else
    double*** globalDiversity;
    double*** globalDiversityRightSolutions;
#endif

//diversity in a run - do I need them is I have global diversity
#ifndef HISTOGRAM
    double* diversity;
    double* diversityTrueDistribution;
    double* diversityRightSolutions;
#else
    //do not use them; use only the global variable
    //double** diversity;
    //double** diversityTrueDistribution;
    //double** diversityRightSolutions;
#endif

    double* transform_table;
    double printDiversityGNUFILE();
    double printDiversityRightSolutions();

    double scoreDiversity();

    void resetGlobalVectors();

#ifdef MIXTURE_REALS
#ifndef HISTOGRAM
    void diversityBookeeping(unsigned long);
    void diversityFileBookeeping(unsigned long);
#endif
#else
    void diversityBookeeping(unsigned long);
    void diversityFileBookeeping(unsigned long);
#endif //mixture_reals
    //cite valori posibile ale functiei avem

#ifdef TRAP_LIKE_FUNCTION
//intentioned fro multiblock functions
    double** globalDiversityDetail;
    double** diversityDetail;
    double printDiversityDetailGNUFILE();
    double scoreDiversityDetail();
#endif

    void resetRuns(unsigned long);

#ifdef HISTOGRAM
//#ifndef EXPAND_KULLBACK_INFORMATION
     void readHistogram();
//#endif
#endif //histogram

#ifdef TWO_ATTRACTORS_FUNCTION
     double* histogram;
#ifdef MULTIPLE_RUNS
     void resetHistogram();
#endif
#ifndef SINGLE_vs_MULTIPLE
    double** ActualHistogram;
    double* globalHistogram;
#else
    double*** ActualHistogram;
    double** globalHistogram;
#endif
#else
#ifdef BINOMIAL
/*    double* histogram;
    void readHistogram();
    double** ActualHistogram;
    double* globalHistogram;
*/
#endif
#endif

  double factorial(unsigned long);
  double combinatorial(unsigned long,unsigned long);

#ifdef ELITIST_ACCEPTANCE
  void proportionalSelection(double[],unsigned long*);
#else
  void proportionalSelection(unsigned long*);
#endif

  void proposedGeneration(unsigned long**, unsigned long);
  void init();
  void mutation();
  void recombination(unsigned long**,unsigned long,const unsigned long);
  void mutation(unsigned long**,unsigned long);

  void randomCOUPLING(unsigned long**);
  void neighboursCOUPLING(unsigned long**);
  double neighboursCOUPLING_NO_OVERLAP(unsigned long**);

#ifndef METROPOLIS_ALGORITHM 
#ifdef POPULATION
  void proposalDistributionCalculation();
  void CalculateNewProposalDistribution(unsigned long);  
  void CalculateNewProposalDistributionRecombination(unsigned long, unsigned long);  
#else
  void proposalDistributionCalculation(unsigned long**, unsigned long);
#endif
#endif

  void acceptanceProposed(unsigned long**, unsigned long);
  unsigned long acceptancePopulation(unsigned long*, unsigned long);
  unsigned long acceptanceCoupledPopulation(unsigned long*);
  unsigned long postProcessing();

public:
  binaryPopulation(unsigned long, unsigned long, char*);
  //binaryPopulation(unsigned long, unsigned long, double, char*);
  ~binaryPopulation();

  void evolution();
#ifdef MULTIPLE_RUNS  
  void evolution(ofstream&);
#endif

  void reset();
  void noCoupling();
//print output for further analyses
#ifdef PLOTS_MCMC
  void printAllChains(ofstream&);
  void printHashtable(ofstream&);
  void printOneChain(ofstream&,unsigned long);
  void printMultiValuedChain(ofstream&);
  void printAllChains_I(ofstream&);
  void printOneChain_I(ofstream&,unsigned long);
  void printMultiValuedChain_I(ofstream&);
#else
#ifdef DEBUG_
  void print();	
  void print(char*);	
  void printLEVEL(ofstream&);
#endif
#endif

#ifdef MULTIPLE_RUNS
  void printMultiRuns(ofstream&);
#endif

#ifdef SINGLE_vs_MULTIPLE
  double updateHashtable();
#endif

#ifndef MIXTURE_REALS
#ifdef DISTANCE
  double HammingDistance(unsigned long*);
  double HammingDistance(unsigned long i);
#endif
#else
#ifdef DISTANCE
  double HammingDistance(double*);
  double HammingDistance(unsigned long i);
  double averageDistanceGoodIndiv();
  double HammingDistance(double*,double*);
#endif
#endif

//studierea effectului recombinarii de la un capat la altul
#ifdef RECOMB_MUT
  double recombRW();
  double mutationRW();
#endif

#ifdef NESTED_MH
  double proposalLocus2_2(unsigned long, unsigned long, unsigned long, unsigned long);
  double proposalDistributionCalculation(double***, int);
#endif //NESTED_MH
//#ifdef MIXTURE_REALS
//  void convertBinary(double*,int*,int);
//  void convertReal(double*,int*,int);
//#endif

#ifdef MIXTURE_REALS
#ifdef BIVARIATE
  unsigned long nrHighIndivGlobal;

#ifdef MIXTURE_BIVARIATE
  //static int dimension;
  //static float boundMax;
  //static float boundMin;
  //static int mixture;
#endif //mixture_bivariate

#endif //BIVARIATE
#endif //MIXTURE_REALS

};

