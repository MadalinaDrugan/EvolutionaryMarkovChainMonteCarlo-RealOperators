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
#include <string>

#include "binaryPopulation_trap.h"
#include "PARAMETERS.h"
using namespace std;

//#include "convergence.cpp"

extern SIMULATIONS *insert_simulation(SIMULATIONS*, long int, double);
extern void show_simulations(SIMULATIONS*);
extern void write_to_file(char*,SIMULATIONS*, SIMULATIONS*);
extern void write_to_file_all(char*,SIMULATIONS*, SIMULATIONS*);
extern void setSeed();
extern long unsigned int genrand_int32();
extern double genrand_real2();
extern double genrand_real1();
extern void convertReal(double*,unsigned long*, long);
extern void convertBinary(double*,unsigned long*, long);

binaryPopulation::binaryPopulation(unsigned long size, unsigned long genSize, char* myFileName){
//random generator 
    //genrand  UniformDistribution genrand_real1(0,1.0);
/*#ifdef MIXTURE_REALS 
    n1 = new NormalDistribution(0.1,0.2);
    n2 = new NormalDistribution(0.4,0.3);
    n3 = new NormalDistribution(0.6,0.1);
    n4 = new NormalDistribution(0.9,0.3);
    nn = new NormalDistribution(0.0,1.0);
#endif //MIXTURE_REALS
*/
//  mix_mutation = 1;
#ifdef BIVARIATE
  nrHighIndivGlobal = 0;
#endif //BIVARIATE

	sizeMCMCstring = size;

#ifndef SINGLE_vs_MULTIPLE
        simul = NULL;

	sampleSize = sampleSize_Const;

	globalRightSolutions = new double[nr_runs];
	for(unsigned long i = 0; i < nr_runs; i++)
	  globalRightSolutions[i] = 0;
#else
	simul = new SIMULATIONS *[(long)(log(sizeMCMCstring)/log(2.0)) + 1];
	for(unsigned long i = 0; i < (long)(log(sizeMCMCstring)/log(2.0)) + 1; i++)
	    simul[i] = NULL;

	sampleSize = sampleSize_Const*sizeMCMCstring;

	globalRightSolutions = new double*[(int)(log(sizeMCMCstring)/log(2.0)) + 1];
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) +1; i++){
	    globalRightSolutions[i] = new double[nr_runs];
	    for(unsigned long j = 0; j < nr_runs; j++)
		globalRightSolutions[i][j] = 0;
	}
#endif

	sizeGENOME = genSize;
	burn_in = 0;
	interSamplesize = 1;
	setSeed();
	averageFitness = 0;
	mean = 0;

#ifdef QBP
	myQBP.open(FileQBP);
	if(!myQBP.is_open()){
	  cout << "Error in opening the qbp File" << FileQBP << "\n";
	}
#ifdef generatedQBP
	matrixQBP = new double*[sizeGENOME];
	for(unsigned long i = 0; i < sizeGENOME; i++){
	    matrixQBP[i] = new double[sizeGENOME];
	    for(unsigned long j = 0; j < sizeGENOME; j++)
		matrixQBP[i][j] = 0;
	}

	generateQBP(sizeGENOME,matrixQBP);
#else
	myQBP >> sizeGENOME;
	readQBP(sizeGENOME,matrixQBP);
#endif //generatedQBP
#endif //QBP

#ifndef SINGLE_vs_MULTIPLE
#ifdef TWO_ATTRACTORS_FUNCTION
	bestFitnessVector = new double*[nr_runs];
	for(unsigned long i =0; i < nr_runs; i++){
	    //two maximum from two picks
	    bestFitnessVector[i] = new double[2];
	    bestFitnessVector[i][0] = 0;
	    bestFitnessVector[i][1] = 0;
	}
#else
	bestFitnessVector = new double[nr_runs];
	for(unsigned long i =0; i < nr_runs; i++)
	    bestFitnessVector[i] = 0;

	secondBestFitnessVector = new double[nr_runs];
	for(unsigned long t =0; t < nr_runs; t++)
	    secondBestFitnessVector[t] = 0;

#endif
#else
#ifdef TWO_ATTRACTORS_FUNCTION
	bestFitnessVector = new double**[(long)(log(sizeMCMCstring)/log(2.0)) + 1];
	for(unsigned long j =0; j < (long)(log(sizeMCMCstring)/log(2.0)) + 1; j++){
	    bestFitnessVector[j] = new double*[nr_runs];
	    for(unsigned long i =0; i < nr_runs; i++){
		//two maximum from two picks
		bestFitnessVector[j][i] = new double[2];
		bestFitnessVector[j][i][0] = 0;
		bestFitnessVector[j][i][1] = 0;
	    }
	}
#else
	bestFitnessVector = new double*[(long)(log(sizeMCMCstring)/log(2.0)) + 1];
	for(unsigned long j =0; j < (long)(log(sizeMCMCstring)/log(2.0)) + 1; j++){
	    bestFitnessVector[j] = new double[nr_runs];
	    for(unsigned long i =0; i < nr_runs; i++)
		bestFitnessVector[j][i] = 0;
	}

#endif
#endif
	cout << "Nr chains:" << size << " sizeGENOME:" << sizeGENOME 
	     << " output:" << myFileName << "\n";

#ifndef FILE_BOOKEEPING
#ifndef MIXTURE_REALS 
	myHashtables = new Hashtable(sizeGENOME);
	myHashtables -> listCurrent = new Hashtable*[sizeGENOME];
	myHashtables -> secondListCurrent = new Hashtable*[sizeGENOME];
	myHashtables -> indexCurrent = new unsigned long[sizeGENOME];
	myHashtables -> secondIndexCurrent = new unsigned long[sizeGENOME];
	unsigned long* temp = new unsigned long[sizeGENOME];
#else
	double* temp = new double[sizeGENOME];
#ifdef HISTOGRAM 
	myHashtables = new Hashtable(sizeGENOME*BLOCKsize);
	myHashtables -> listCurrent = new Hashtable*[sizeGENOME*BLOCKsize];
	myHashtables -> secondListCurrent = new Hashtable*[sizeGENOME*BLOCKsize];
	myHashtables -> indexCurrent = new unsigned long[sizeGENOME*BLOCKsize];
	myHashtables -> secondIndexCurrent = new unsigned long[sizeGENOME*BLOCKsize];
#else // histogram
    //aditional hashtables
	//cout << " am fost in constructor \n";
	if(nr_runs == 1){
	  myHashtables = new Hashtable_histogram*[(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
	  for(unsigned long i = 0; i < (unsigned long)(sampleSize_Const/SEE_SAMPLE)+1; i++){
	    myHashtables[i] = new Hashtable_histogram(sizeGENOME,NrBins);
#ifdef FIRST_ITERATOR
	    myHashtables[i] -> listCurrent = new Hashtable_histogram*[sizeGENOME];
	    myHashtables[i] -> indexCurrent = new unsigned long[sizeGENOME];
#endif //first iterator
#ifdef SECOND_ITERATOR
	    myHashtables[i] -> listSecondCurrent = new Hashtable_histogram*[sizeGENOME];
	    myHashtables[i] -> indexSecondCurrent = new unsigned long[sizeGENOME];
#endif // second iterator
	  }
	} else {
	  myHashtables = new Hashtable_histogram*[nr_runs];
	  for(unsigned long i = 0; i < nr_runs; i++){
	    myHashtables[i] = new Hashtable_histogram(sizeGENOME,NrBins);
#ifdef FIRST_ITERATOR
	    myHashtables[i] -> listCurrent = new Hashtable_histogram*[sizeGENOME];
	    myHashtables[i] -> indexCurrent = new unsigned long[sizeGENOME];
#endif // first iterator
#ifdef SECOND_ITERATOR
	    myHashtables[i] -> listSecondCurrent = new Hashtable_histogram*[sizeGENOME];
	    myHashtables[i] -> indexSecondCurrent = new unsigned long[sizeGENOME];
#endif // second iterator
	  }      
	}

#endif //histogram
#endif //mixture reals
#endif //file bookeeping

	MCMCstring = new binaryMCMC*[sizeMCMCstring];
	for(unsigned long i =0; i<sizeMCMCstring; i++){
	  //cout << "generaza " << i << " with size:" << sizeGENOME<<"\n";
#ifndef SINGLE_vs_MULTIPLE
#ifdef BINOMIAL
	    MCMCstring[i] = new binaryMCMC(sizeGENOME,DEFAULT_TEMPERATURE,i);
#else
#ifdef QBP
	    MCMCstring[i] = new binaryMCMC(sizeGENOME,DEFAULT_TEMPERATURE,matrixQBP);
#else  //not BINOMIAL or QbP
	    MCMCstring[i] = new binaryMCMC(sizeGENOME,DEFAULT_TEMPERATURE);
#ifndef FILE_BOOKEEPING
	    MCMCstring[i]->getGENOME(temp);
	    //myHashtables->store_individual(temp,sizeGENOME);
#endif //FILE_BOOKEEPOING
	     averageFitness += MCMCstring[i]->getValue(); 
#endif //QBP
#endif //BINOMIAL
#else
	    MCMCstring[i] = new binaryMCMC(sizeGENOME,DEFAULT_TEMPERATURE,i,sizeMCMCstring);
#endif //SINGLE_vs_MULTIPLE
	}	

	averageFitness /= sizeMCMCstring;

#ifndef FILE_BOOKEEPING
	delete[] temp;
#endif

	generation = 0;
	if(myFileName != NULL) {
		fileData = new char[strlen(myFileName)+1];
		strcpy(fileData,myFileName);
	} 
	
	proposedMCMCstring = new binaryMCMC*[sizeMCMCstring];
	for(unsigned long i =0; i<sizeMCMCstring; i++){
	  //cout << "generaza proposed independent not count" << i <<"\n";
		proposedMCMCstring[i] = new binaryMCMC(sizeGENOME,1);
	}

	runs = 0;
	procentRightSolutions = 0;
	weightRightSolutions = 0;

#ifndef METROPOLIS_ALGORITHM
#ifdef POPULATION_MCMC
	PropVector = new double*[sizeGENOME];
	for(unsigned long i = 0; i < sizeGENOME; i++){
	    PropVector[i] = new double[2];
	    PropVector[i][0] = 1;
	    PropVector[i][1] = 1;
	}
#else
	proposalDistribution = new double[sizeMCMCstring+2];
	proposedProposalDistribution = new double[sizeMCMCstring+2];
	for(unsigned long k = 0; k < sizeMCMCstring; k++){
	  proposalDistribution[k] = 1;
	  proposedProposalDistribution[k] = 1;
	}
#endif // POPULATION_MCMC
#endif //Metropolis_algorithm
	
#ifdef  ELITIST_ACCEPTANCE_WITH_REGENERATION	
	tempResults = new double*[sizeCoupling];
	for(unsigned long i =0; i<sizeCoupling; i++)
	    tempResults[i] = new double[2*sizeCoupling];
#endif //ELITIST_ACCEPTANCE_WITH_REGENERATION

////////////////////////
//data to measure the diversity of the generated individuals
//////////////////////

#ifdef FILE_BOOKEEPING
	diversityFile.open(diversityTempFILE,ios::trunc|ios::out|ios::in);
	if(!diversityFile.is_open()){
	  cout << "Error in opening the diveristy File" <<diversityTempFILE<< "\n";
	}
#else 
	diversityFile.open(diversityTempFILE);
	if(!diversityFile.is_open()){
	  cout << "Error in opening the diveristy File" <<diversityTempFILE<< "\n";
	}	
#endif // FILE_BOOKEEPING

	//despre forma functiilor
	//numarul valorilor solutiilor posibile
#ifdef TRAP_LIKE_FUNCTION

#ifdef TRAP_FUNCTION
//ATENTIE!!! a si b trebuie sa aiba valori intregi altfel trebuie modificata formula
#ifdef LINEAR_TRAP_FUNCTION
	diversityMeasure = a*sizeGENOME/BLOCKsize + 1;
	transform_table = new double[diversityMeasure];
      	//transform_table[0] = 0.1;
	for(unsigned long t = 0; t < diversityMeasure; t++)
	    transform_table[t] = t;

	//diversityMeasure = sizeGENOME;
#else
	if(ORDER_BLOCK != 1){
	    diversityMeasure = 0;
	    for(unsigned long i = 0; i < (unsigned long)(sizeGENOME/BLOCKsize); i++)
		diversityMeasure += (unsigned long) pow(i+1,ORDER_BLOCK);
	    diversityMeasure = a * diversityMeasure + 1;
	}
	else
	    diversityMeasure = a*(sizeGENOME/BLOCKsize)*(sizeGENOME/BLOCKsize + 1)/2 + 1;
	transform_table = new double[diversityMeasure];
      	//transform_table[0] = 0.1;
	for(unsigned long t = 0; t < diversityMeasure; t++)
	    transform_table[t] = t;
#endif // LINEAR_TRAP_FUNCTION
#endif //TRAP_FUNCTION

	//#else

#ifdef MULTI_PICK_TRAP_FUNCTION
	diversityMeasure = b*sizeGENOME/BLOCKsize + 1;
	//cout << "diversityMeasure" <<diversityMeasure <<"\n";
#endif //MULTI_PICK_TRAP_FUNCTION

	diversityDetail = new double*[sizeGENOME/BLOCKsize];
	globalDiversityDetail = new double*[sizeGENOME/BLOCKsize];
	for(unsigned long i = 0; i < sizeGENOME/BLOCKsize; i++){
	  diversityDetail[i] = new double[BLOCKsize+1];
	  globalDiversityDetail[i] = new double[BLOCKsize+1];
	  for(unsigned long j = 0; j < BLOCKsize+1; j++){
	    diversityDetail[i][j] = 0;
	    globalDiversityDetail[i][j] = 0;
	  }
	}

#endif //TRAP_LIKE_FUNCTION


#ifdef TWO_ATTRACTORS_FUNCTION
	diversityMeasure = wrongBits+1;
	cout << "\n Diversity" << diversityMeasure << "\n";
	transform_table = new double[diversityMeasure];
      	//transform_table[0] = 0.1;
	for(unsigned long t = 0; t < diversityMeasure; t++)
	    transform_table[t] = t;

	histogram = new double[sizeGENOME+1];
	//readHistogram();

#ifndef SINGLE_vs_MULTIPLE
	ActualHistogram = new double*[nr_runs];
	for(unsigned long i = 0; i < nr_runs; i++){
	    ActualHistogram[i] = new double[sizeGENOME+1];
	    for(unsigned long j = 0; j < sizeGENOME+1; j++)
		ActualHistogram[i][j] = 0;
	}
	
	globalHistogram = new double[sizeGENOME+1];
	for(unsigned long j = 0; j < sizeGENOME+1; j++)
	    globalHistogram[j] = 0;
#else
	meanFitness = new double[(unsigned long)(log(sizeMCMCstring)/log(2.0))+1];
	for(unsigned long k = 0; k < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; k++)
	    meanFitness[k] = 0;

	diffFitness = new double[(unsigned long)(log(sizeMCMCstring)/log(2.0))+1];
	for(unsigned long k = 0; k < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; k++)
	    diffFitness[k] = 0;

	ActualHistogram = new double**[(unsigned long)(log(sizeMCMCstring)/log(2.0))+1];
	for(unsigned long k = 0; k < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; k++){
	    ActualHistogram[k] = new double*[nr_runs];
	    for(unsigned long i = 0; i < nr_runs; i++){
		ActualHistogram[k][i] = new double[sizeGENOME+1];
		for(unsigned long j = 0; j < sizeGENOME+1; j++)
		    ActualHistogram[k][i][j] = 0;
	    }
	}
	
	globalHistogram = new double*[(unsigned long)(log(sizeMCMCstring)/log(2.0))+1];
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++){
	    globalHistogram[i] = new double[sizeGENOME+1];
	    for(unsigned long j = 0; j < sizeGENOME+1; j++)
		globalHistogram[i][j] = 0;
	}
#endif
#endif // TWO_ATTRACTORS_FUNCTION

#ifdef ONEMAX_FUNCTION
	diversityMeasure = sizeGENOME + 1;
#endif //ONEMAX_FUNCTION

#ifdef BINOMIAL
#ifdef BERNOULLI
	diversityMeasure = (2 * BLOCKsize - 1) * sizeGENOME / BLOCKsize + 1;
#else 
#ifdef LINEAR_TRAP_FUNCTION
	diversityMeasure = sizeGENOME+1;
#else
	if(ORDER_BLOCK != 1){
	    diversityMeasure = 0;
	    for(unsigned long i = 0; i < (unsigned long)(sizeGENOME/BLOCKsize); i++)
		diversityMeasure += (unsigned long) pow(i+1,ORDER_BLOCK);
	    diversityMeasure = BLOCKsize * diversityMeasure + 1;
	}
	else 
	    diversityMeasure = BLOCKsize * (sizeGENOME/BLOCKsize) * (sizeGENOME/BLOCKsize + 1) / 2 + 1;
#endif
	transform_table = new double[diversityMeasure];
      	//transform_table[0] = 0.1;
	for(unsigned long t = 0; t < diversityMeasure; t++)
	    transform_table[t] = t;


	/*histogram = new double[sizeGENOME+1];

	ActualHistogram = new double*[nr_runs];
	for(unsigned long i = 0; i < nr_runs; i++){
	    ActualHistogram[i] = new double[sizeGENOME+1];
	    for(unsigned long j = 0; j < sizeGENOME+1; j++)
		ActualHistogram[i][j] = 0;
	}
	
	globalHistogram = new double[sizeGENOME+1];
	for(unsigned long j = 0; j < sizeGENOME+1; j++)
	    globalHistogram[j] = 0;
	*/
#endif
#endif //BINOMIAL

#ifdef MIXTURE_REALS
//each individual is diverse in the limit of the binary individuals
/////???????????????????? should be diversity + 1
#ifndef HISTOGRAM
	//diversityMeasure = sizeGENOME*BLOCKsize;
#else
	diversityMeasure = (unsigned long) pow((double)hist_points,(double)sizeGENOME);
#endif //HISTOGRAM
#endif //MIXTURE_REALS

#ifdef HISTOGRAM
	//if(burn_intT != NULL) delete[] burn_inT;
	//proprietati pentru diveristyMeasure
	if(diversityMeasure == 0) {
	  cout << "Error in diversityMeasure =0 \n";
	  exit(1);
	}
#endif //histogram

#ifndef HISTOGRAM
	/*
	diversity = new double[diversityMeasure];
	diversityRightSolutions = new double[diversityMeasure];
	diversityTrueDistribution = new double[diversityMeasure];
	for(unsigned long i = 0; i <  diversityMeasure; i++){
	    diversityRightSolutions[i] = 0;
	    diversity[i] = 0;
	    diversityTrueDistribution[i] = 0;
	}
	*/
#else
	
//	diversityRightSolutions = new double*[diversityMeasure];
	cout << " diversity measure " << diversityMeasure << "\n";
	diversityTrueDistribution = new double[diversityMeasure];
//	diversity = new double*[diversityMeasure]; //number of individuals in a bin
	for(unsigned long i = 0; i < diversityMeasure; i++){
//		diversityRightSolutions[j][i] = 0;
		diversityTrueDistribution[i] = 0;
//		diversity[j][i] = 0;
	}
	
#endif //histogram
	      
#ifndef SINGLE_vs_MULTIPLE
#ifndef MIXTURE_REALS
	globalDiversityRightSolutions = new double*[nr_runs];
	globalDiversity = new double*[nr_runs];
	for(unsigned long i = 0; i <  nr_runs; i++){
	    globalDiversityRightSolutions[i] = new double[diversityMeasure];
	    globalDiversity[i] = new double[diversityMeasure];
	    for(unsigned long j = 0; j < diversityMeasure; j++){
		globalDiversityRightSolutions[i][j] = 0;
		globalDiversity[i][j] = 0;
	    }
	}

#else
#ifdef HISTOGRAM
	globalDiversityRightSolutions = new double**[(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
	globalDiversity = new double**[(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
	for(unsigned long t = 0; t < (unsigned long)(sampleSize_Const/SEE_SAMPLE)+1; t++){
	  globalDiversityRightSolutions[t] = new double*[nr_runs];
	  globalDiversity[t] = new double*[nr_runs];
	  for(unsigned long i = 0; i < nr_runs; i++){
	    globalDiversityRightSolutions[t][i] = new double[diversityMeasure];
	    globalDiversity[t][i] = new double[diversityMeasure];
	    for(unsigned long j = 0; j < diversityMeasure; j++){
		globalDiversityRightSolutions[t][i][j] = 0;
		globalDiversity[t][i][j] = 0;
	    }
	  }
	}
#endif //histogram
#endif
#else
	globalDiversityRightSolutions = new double**[(unsigned long)(log(sizeMCMCstring)/log(2.0))+1];
	globalDiversity = new double**[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1];
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++){
	    globalDiversityRightSolutions[i] = new double*[nr_runs];
	    globalDiversity[i] = new double*[nr_runs];
	    for(unsigned long j = 0; j < nr_runs; j++){
		globalDiversityRightSolutions[i][j] = new double[diversityMeasure];
		globalDiversity[i][j] = new double[diversityMeasure];
		for(unsigned long k = 0; k < diversityMeasure; k++){
		    globalDiversityRightSolutions[i][j][k] = 0;
		    globalDiversity[i][j][k] = 0;
		}
	    }
	}
#endif


//////////////////////
//date pentru masurarea convergentei
//////////////////////
#ifdef ACCEPTANCE_RATIO
#ifdef SINGLE_vs_MULTIPLE
#ifdef TWO_ATTRACTORS_FUNCTION
	meanAcceptance = new double[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1];
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++)
	    meanAcceptance[i] = 0;
#endif
	acceptance_ratio = new double*[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1];
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++){
	    acceptance_ratio[i] = new double[nr_runs];
	    for(unsigned long j = 0; j < nr_runs; j++)
		acceptance_ratio[i][j] = 0;
	}
#else
#ifndef multiple_restarts
	acceptance_ratio = new double[(long)((double)sampleSize_Const/(double)SEE_SAMPLE) + 1];
	for(unsigned long i = 0; i < (long)((double)sampleSize_Const/(double)SEE_SAMPLE) + 1; i++)
	  acceptance_ratio[i] = 0;	  
#else
	//if(nr_runs != 1){
	acceptance_ratio = new double[nr_runs];
	for(unsigned long i = 0; i < nr_runs; i++)
	  acceptance_ratio[i] = 0;
#endif
#endif

#endif

#ifdef CUSUM
	burn_inT = new unsigned long[nr_runs];
	for(unsigned long i = 0; i < nr_runs; i++){
	  burn_inT[i] = 0;
	}

#ifndef SINGLE_vs_MULTIPLE
#ifndef multiple_restarts
	  Cusum = new double[sampleSize_Const+4];
	  Hairness = new double[sampleSize_Const+4];
	  for(unsigned long j = 0 ; j < sampleSize_Const+4; j++){
	      Cusum[j] = 0; 
	      Hairness[j] = 0;
	  }
#else  
	  Cusum = new double*[nr_runs];
	  Hairness = new double*[nr_runs];
	  for(unsigned long i = 0; i < nr_runs; i++){
	    Cusum[i] = new double[sampleSize_Const+4];
	    Hairness[i] = new double[sampleSize_Const+4];
	    for(unsigned long j = 0 ; j < sampleSize_Const+4; j++){
	      Cusum[i][j] = 0; 
	      Hairness[i][j] = 0;
	    } 
	  }
#endif //nr_runs
#else
	Cusum = new double**[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1];
	Hairness = new double**[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1];
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1; i++){
	    Cusum[i] = new double*[nr_runs];
	    Hairness[i] = new double*[nr_runs];
	    for(unsigned long j = 0; j < nr_runs; j++){
		Cusum[i][j] = new double[sampleSize+4];
		Hairness[i][j] = new double[sampleSize+4];
		for(unsigned long k = 0 ; k < sampleSize+4; k++){
		    Cusum[i][j][k] = 0; 
		    Hairness[i][j][k] = 0;
		} 
	    }
	}
#endif
#endif //CUSUM

#ifdef HISTOGRAM
#ifdef KULLBACK_INFORMATION
	Kullback = new double*[nr_runs+1];
	for(unsigned long i = 0; i < nr_runs+1; i++){
	  Kullback[i] = new double[(unsigned long)(sampleSize_Const/SEE_SAMPLE) + 1]; 
	}

#ifndef HISTOGRAM
	trueDistribution = new double[diversityMeasure];
	for(unsigned long i = 0; i < diversityMeasure; i++)
	  trueDistribution[i] = 0; 
#else
	trueDistribution = new double[(unsigned long) pow((double)diversityMeasure,(double)sizeGENOME)];
	for(unsigned long i = 0; i < (unsigned long) pow((double)diversityMeasure,(double)sizeGENOME); i++){
	      trueDistribution[i] = 0; 
	} 
#endif

#ifdef TIME_PERFORMANCE
	timeKLdistance = new double*[nr_runs+1];
	for(unsigned long i = 0; i < nr_runs+1; i++)
	  timeKLdistance[i] = new double[(unsigned long)(sampleSize_Const/SEE_SAMPLE) + 1]; 
#endif
#endif //KULLBACK_INFORMATION
#endif //histogram

#ifdef GELMAN_RUBIN
	teta = new double*[nr_runs];
	s_2 = new double*[nr_runs];
	for(unsigned long i = 0; i< nr_runs; i++){
	  teta[i] = new double[sampleSize + 4];
	  for(unsigned long j = 0; j < sampleSize + 1; j++)
	    teta[i][j] = 0;
#ifdef HALF_GELMAN
	  s_2[i] = new double[(sampleSize + 4)/2];
	  for(unsigned long j = 0; j < (sampleSize+1)/2; j++)
	    s_2[i][j] = 0;
#else
	  s_2[i] = new double[sampleSize+4];
	  for(unsigned long j = 0; j < sampleSize; j++)
	    s_2[i][j] = 0;
#endif
	}

	tetaM = new double**[nr_runs];
	for(unsigned long i = 0; i < nr_runs; i++){
	  tetaM[i] = new double*[sizeMCMCstring];
	  for(unsigned long j = 0; j < sizeMCMCstring; j++)
	    tetaM[i][j] = new double[sampleSize + 1];
	}

#endif //GELMAN_RUBIN
       	//cout << "sizeGenome 0:" << MCMCstring[0]->size << "\n";
}


binaryPopulation::~binaryPopulation()
{
    //  delete genrand_real1;
/*#ifdef MIXTURE_REALS
    delete n1;
    delete n2;
    delete n3;
    delete n4;
    delete nn;
#endif
*/
#ifdef FILE_BOOKEEPING
  diversityFile.close();  
#endif

#ifndef MIXTURE_REALS
  if(myHashtables!= NULL) myHashtables->free_table(myHashtables);
  //delete myHashtables;
#else
#ifdef HISTOGRAM
  if(myHashtables!= NULL) myHashtables->free_table(myHashtables);
  //delete myHashtables;
#else
  //for(int i = 0; i < )
    //if(myHashtables!= NULL) myHashtables->free_table(myHashtables);
  delete[] myHashtables;
#endif //histogram
#endif //mixture reals

#ifndef HISTOGRAM
  if(diversity != NULL)   delete[] diversity;
  if(diversityRightSolutions != NULL)  delete[] diversityRightSolutions;
  if(diversityTrueDistribution != NULL)  delete[] diversityTrueDistribution;
#endif

#ifdef MIXTURE_REALS
#ifdef HISTOGRAM
  if(globalDiversityRightSolutions != NULL)  delete[] globalDiversityRightSolutions;
  if(globalRightSolutions != NULL) delete[] globalRightSolutions;
  if(globalDiversity!= NULL) delete[] globalDiversity;
#endif 
#else
  if(globalDiversityRightSolutions != NULL)  delete[] globalDiversityRightSolutions;
  if(globalRightSolutions != NULL) delete[] globalRightSolutions;
  if(globalDiversity!= NULL) delete[] globalDiversity;
#endif //mixture reals

#ifdef TRAP_LIKE_FUNCTION
  if(diversityDetail != NULL) delete[] diversityDetail;
  if(globalDiversityDetail != NULL) delete[] globalDiversityDetail;
#endif

///////////////////////
///masurarea convergentei
//////////////////////
#ifdef ACCEPTANCE_RATIO
   if(acceptance_ratio != NULL)  delete[] acceptance_ratio; 
#ifdef SINGLE_vs_MULTIPLE
#ifdef TWO_ATTRACTORS_FUNCTION
   if(meanAcceptance != NULL) delete[] meanAcceptance;
#endif
#endif
#endif

#ifdef HISTOGRAM
#ifdef KULLBACK_INFORMATION
      if(Kullback != NULL) delete Kullback;
      if(trueDistribution != NULL) delete[] trueDistribution;
#ifdef TIME_PERFORMANCE
      if(timeKLdistance != NULL) delete[] timeKLdistance;
#endif
#endif
#endif //histogram

      if(bestFitnessVector != NULL) delete[] bestFitnessVector;
      if(secondBestFitnessVector != NULL) delete[] secondBestFitnessVector;

#ifdef GELMAN_RUBIN
  delete[] teta;
  delete[] s_2;
  delete[] tetaM;
#endif

#ifdef CUSUM
  if(Cusum != NULL) delete[] Cusum;
  if(Hairness != NULL) delete[] Hairness;
  //if(burn_intT != NULL) delete[] burn_inT;
#endif

  if(fileData != NULL) delete[] fileData;
  delete[] MCMCstring;
  delete[] proposedMCMCstring;

#ifndef METROPOLIS_ALGORITHM
#ifdef POPULATION_MCMC
  if(PropVector != NULL) delete[] PropVector;
#else
  if(proposalDistribution != NULL) delete[] proposalDistribution; 
  if(proposedProposalDistribution != NULL) delete[] proposedProposalDistribution; 
#endif
#endif //METROPOLIS_ALGORITHM

#ifdef  ELITIST_ACCEPTANCE_WITH_REGENERATION
  delete[] tempResults; 
#endif

  delete[] transform_table;

#ifdef TWO_ATTRACTORS_FUNCTION
    delete[] histogram;
    delete[] ActualHistogram;
    delete[] globalHistogram;
#ifdef SINGLE_vs_MULTIPLE
    delete[] meanFitness;
    delete[] diffFitness;
#endif
#endif

#ifdef QBP
  if(matrixQBP != NULL) delete[] matrixQBP;
#endif

}

// -------- Useful functions ----------

void binaryPopulation::reset(){
//reset the seed
//    unsigned long seed 1812433253;
//    Ran002 gen(seed);

    generation = 0;

    for(unsigned long i =0; i<sizeMCMCstring; i++){
#ifdef BINOMIAL
	MCMCstring[i]->reset(runs,i);
#else
	MCMCstring[i] -> reset(runs);
#endif
    }	

#ifdef  ELITIST_ACCEPTANCE_WITH_REGENERATION
    for(unsigned long i =0; i<sizeCoupling; i++)
	for(unsigned long j = 0; j < 2*sizeCoupling; j++)
	    tempResults[i][j] = 0;	
#endif
#ifdef HISTOGRAM
    if(diversityMeasure != 0){
#ifndef HISTOGRAM
	if( diversity != NULL && diversityRightSolutions !=NULL)
	    for(unsigned long i = 0; i < diversityMeasure; i++){
		diversity[i] = 0;
		diversityRightSolutions[i] = 0;
	    }
#endif
    }
#endif //histogram

#ifndef MIXTURE_REALS
  if(myHashtables != NULL){
    myHashtables->reset(sizeGENOME);
    //delete myHashtables;
    //myHashtables = new Hashtable(sizeGENOME);
  }
#else 
#ifdef HISTOGRAM
  if(myHashtables != NULL){
    myHashtables->reset(sizeGENOME*BLOCKsize);
#else
    if(nr_runs == 1){
      for(unsigned long i = 0; i < (unsigned long)(sampleSize_Const/SEE_SAMPLE)+1; i++){
	myHashtables[i] -> reset();
      }
    } else {
      /*for(unsigned long i = 0; i < nr_runs; i++){
	myHashtables[i] -> reset();
	} */     
    }
#endif //histogram
#endif //mixture reals

#ifdef TRAP_LIKE_FUNCTION
  for(unsigned long i = 0; i < sizeGENOME/BLOCKsize; i++)
    for(unsigned long j = 0; j < BLOCKsize+1; j++)
      diversityDetail[i][j] = 0;
#endif

/*#ifdef KULLBACK_INFORMATION
  for(unsigned long i = 0; i < sampleSize; i++)
    Kullback[i] = 0; 
    #endif  */

  procentRightSolutions = 0;
  weightRightSolutions = 0;

#ifdef FILE_BOOKEEPING
  diversityFile.flush();
  diversityFile.close();
  diversityFile.open(diversityTempFILE,ios::trunc|ios::out|ios::in);
  if(!diversityFile.is_open()){
    cout << "Error in opening after reset the diversityFile" << diversityTempFILE << "\n";
  }
#endif

#ifndef METROPOLIS_ALGORITHM
#ifdef POPULATION_MCMC
  for(unsigned long i = 0; i < sizeGENOME; i++){
      PropVector[i][0] = 1;
      PropVector[i][1] = 1;
  }
#else
  for(unsigned long j = 0; j < sizeMCMCstring; j++){
	  proposalDistribution[j] = 1; 
	  proposedProposalDistribution[j] = 1; 
      }
#endif
#endif

  averageFitness = 0;
}

#ifdef MULTIPLE_RUNS
void binaryPopulation::resetRuns(unsigned long size){
    runs = 0;
    sizeMCMCstring = size;
    reset();
    //reseteaza tablourile globale
#ifdef CUSUM
    resetCusum();
#endif
#ifdef KULLBACK_INFORMATION
    resetKullback();
#endif
    resetGlobalVectors();
#ifdef TWO_ATTRACTORS_FUNCTION
    resetHistogram();
#endif
}

#ifndef MIXTURE_REALS
void binaryPopulation::resetGlobalVectors(){
#ifndef SINGLE_vs_MULTIPLE
        simul = NULL;

	sampleSize = sampleSize_Const;

	for(unsigned long i = 0; i < nr_runs; i++)
	  globalRightSolutions[i] = 0;
#else
	for(unsigned long i = 0; i < (long)(log(sizeMCMCstring)/log(2.0)) + 1; i++)
	    simul[i] = NULL;

	sampleSize = sampleSize_Const*sizeMCMCstring;

	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) +1; i++)
	    for(unsigned long j = 0; j < nr_runs; j++)
		globalRightSolutions[i][j] = 0;
	
#endif

	burn_in = 0;
	interSamplesize = 1;
	setSeed();
	averageFitness = 0;
	mean = 0;
	
#ifndef SINGLE_vs_MULTIPLE
#ifdef TWO_ATTRACTORS_FUNCTION
	for(unsigned long i =0; i < nr_runs; i++){
	    //two maximum from two picks
	    bestFitnessVector[i][0] = 0;
	    bestFitnessVector[i][1] = 0;
	}
#else
	for(unsigned long i =0; i < nr_runs; i++)
	    bestFitnessVector[i] = 0;
#endif
#else
#ifdef TWO_ATTRACTORS_FUNCTION
	for(unsigned long j =0; j < (long)(log(sizeMCMCstring)/log(2.0)) + 1; j++)
	    for(unsigned long i =0; i < nr_runs; i++){
		//two maximum from two picks
		bestFitnessVector[j][i][0] = 0;
		bestFitnessVector[j][i][1] = 0;
	    }
#else
	for(unsigned long j =0; j < (long)(log(sizeMCMCstring)/log(2.0)) + 1; j++)
	    for(unsigned long i =0; i < nr_runs; i++)
		bestFitnessVector[j][i] = 0;
#endif
#endif
	procentRightSolutions = 0;
	weightRightSolutions = 0;

	for(unsigned long i = 0; i < nr_runs; i++){
	  burn_inT[i] = 0;
	}

#ifdef TRAP_LIKE_FUNCTION
	for(unsigned long i = 0; i < sizeGENOME/BLOCKsize; i++){
	  for(unsigned long j = 0; j < BLOCKsize+1; j++)
	    globalDiversityDetail[i][j] = 0;
	}
#endif //TRAP_FUNCTION

#ifdef BINOMIAL
	for(unsigned long t =0; t < nr_runs; t++)
	    secondBestFitnessVector[t] = 0;
#endif

	  for(unsigned long i = 0; i <  diversityMeasure; i++)
	      diversityRightSolutions[i] = 0;
	      
#ifndef SINGLE_vs_MULTIPLE
	  for(unsigned long i = 0; i <  nr_runs; i++)
	      for(unsigned long j = 0; j < nr_runs; j++)
		  globalDiversityRightSolutions[i][j] = 0;
#else
	  for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++)
	      for(unsigned long j = 0; j < nr_runs; j++)
		  for(unsigned long k = 0; k < diversityMeasure; k++)
		      globalDiversityRightSolutions[i][j][k] = 0;
#endif

#ifndef SINGLE_vs_MULTIPLE
	for(unsigned long i = 0; i <  nr_runs; i++)
	    for(unsigned long j = 0; j < diversityMeasure; j++)
		globalDiversity[i][j] = 0;
#else
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1; i++)
	    for(unsigned long j = 0; j < nr_runs; j++)
		  for(unsigned long k = 0; k < diversityMeasure; k++)
		      globalDiversity[i][j][k] = 0;
#endif
	for(unsigned long i = 0; i < diversityMeasure; i++)
	    diversityTrueDistribution[i] = 0;

#ifdef ACCEPTANCE_RATIO
#ifdef SINGLE_vs_MULTIPLE
#ifdef TWO_ATTRACTORS_FUNCTION
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++)
	    meanAcceptance[i] = 0;
#endif
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++){
	    for(unsigned long j = 0; j < nr_runs; j++)
		acceptance_ratio[i][j] = 0;
	}
#else
	for(unsigned long i = 0; i < nr_runs; i++){
	    acceptance_ratio[i] = 0;
	}
#endif

#endif

#ifndef METROPOLIS_ALGORITHM
#ifdef POPULATION_MCMC
#else
	for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    proposalDistribution[i] = 1;
	    proposedProposalDistribution[i] = 1;
	}
#endif
#endif // METROPOLIS_ALGORITHM
	averageFitness = 0;
}
#endif //mixture reals

#ifdef GELMAN
void binaryPopulation::resetGelman(){

	for(unsigned long i = 0; i< nr_runs; i++)
	  for(unsigned long j = 0; j < sampleSize + 1; j++)
	    teta[i][j] = 0;

#ifdef HALF_GELMAN
	 
	  for(unsigned long j = 0; j < (sampleSize+1)/2; j++)
	    s_2[i][j] = 0;
#else
	  
	  for(unsigned long j = 0; j < sampleSize; j++)
	    s_2[i][j] = 0;
#endif
}
#endif //GELMAN_RUBIN

#ifdef KULLBACK_INFORMATION
void binaryPopulation::resetKullback(){
	for(unsigned long j = 0; j < (unsigned long)(sampleSize/SEE_SAMPLE)+1; j++)
	    for(unsigned long i = 0; i < nr_runs; i++)
	    Kullback[i][j] = 0; 
	for(unsigned long i = 0; i < diversityMeasure; i++)
	  trueDistribution[i] = 0;
} 
#endif //KULLBACK_INFORMATION

 
#ifdef CUSUM
void binaryPopulation::resetCusum(){
#ifndef SINGLE_vs_MULTIPLE
	for(unsigned long i = 0; i < nr_runs; i++){
	  for(unsigned long j = 0 ; j < sampleSize; j++){
	    Cusum[i][j] = 0; 
	    Hairness[i][j] = 0;
	  } 
	}
#else
	for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1; i++){
	    for(unsigned long j = 0; j < nr_runs; j++){
		for(unsigned long k = 0 ; k < sampleSize+2; k++){
		    Cusum[i][j][k] = 0; 
		    Hairness[i][j][k] = 0;
		} 
	    }
	}
#endif
}
#endif //CUSUM

#ifdef TWO_ATTRACTORS_FUNCTION
void binaryPopulation::resetHistogram(){
#ifndef SINGLE_vs_MULTIPLE
    for(unsigned long i = 0; i < nr_runs; i++){
	for(unsigned long j = 0; j < sizeGENOME+1; j++)
	    ActualHistogram[i][j] = 0;
    }

    for(unsigned long j = 0; j < sizeGENOME+1; j++)
	globalHistogram[j] = 0;
#else
    for(unsigned long k = 0; k < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; k++){
	    meanFitness[k] = 0;
	    diffFitness[k] = 0;
    }
    for(unsigned long k = 0; k < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; k++)
	for(unsigned long i = 0; i < nr_runs; i++)
	    for(unsigned long j = 0; j < sizeGENOME+1; j++)
		ActualHistogram[k][i][j] = 0;

    for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++)
	for(unsigned long j = 0; j < sizeGENOME+1; j++)
	    globalHistogram[i][j] = 0;
    
#endif
}
#endif // TWO_ATTRACTORS_FUNCTION

#endif //MULTIPLE RUNS

//////////------

//factorial - (N)!
double binaryPopulation::factorial(unsigned long nr){
     if(nr < 1) return 1;
	return nr * factorial(nr-1);
}

//combinatorial C_{n}^{i}
double binaryPopulation::combinatorial(unsigned long size, unsigned long i)
{
	return factorial(size) / (factorial(i) * factorial(size - i));
}

#ifdef PARALLEL_SIMULATED_ANNEALING
//the information about temperature is stored in bestTemperature parameter
void binaryPopulation::assignTemperature()
{
  double bestTemperature = bestTEMPERATURE;
  double worstTemperature = worstTEMPERATURE;

  //linear schedual
#ifdef LINEAR_COOLING
  double sliceT = (bestTemperature - worstTemperature) / (double) (burn_in + sampleSize*interSamplesize);
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    MCMCstring[i] ->setTEMPERATURE(bestTemperature - generation*sliceT);
  //cout << "Temperature at "<<generation <<" is " << (bestTemperature - generation*sliceT) << "\n";
#else
#ifdef LOGARITHMIC_COOLING
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    MCMCstring[i] ->setTEMPERATURE(bestTemperature / log(generation + 2));
#else
#ifdef GEOMETRIC_COOLING
  double base = exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)(burn_in + sampleSize*interSamplesize)));
  //  cout << "base " << base;
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    if(generation == 0) MCMCstring[i] ->setTEMPERATURE(bestTemperature);
      else MCMCstring[i] ->setTEMPERATURE(bestTemperature * pow(base,generation));
#endif
#endif
#endif
}
#else 
#ifdef PARALLEL_TEMPERING
void binaryPopulation::assignTemperature(){

#ifdef LINEAR_COOLING  
  double sliceT = (bestTemperature - worstTemperature) / (double)sizeMCMCstring;
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(worstTemperature + sliceT*(i+0.5));
#else
#ifdef LOGARITHMIC_COOLING
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(bestTemperature / log(i+2));
#else
#ifdef GEOMETRIC_COOLING
  double bestTemperature = bestTEMPERATURE;
  double worstTemperature = worstTEMPERATURE;
  double base = exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)sizeMCMCstring));

  for(unsigned long i = 0; i < sizeMCMCstring; i++){
#ifdef MULTIPLE_RUNS
	MCMCstring[i] -> setTEMPERATURE(bestTemperature);
#else
    MCMCstring[i] -> setTEMPERATURE(bestTemperature * pow(base,i));
#endif
  }
#else
#ifdef NO_COOLING
  for(unsigned long i = 0; i < sizeMCMCstring; i++){
#ifdef MULTIPLE_RUNS
	MCMCstring[i] -> setTEMPERATURE(bestTemperature);
#else
    MCMCstring[i] -> setTEMPERATURE(DEFAULT_TEMPERATURE);
    //MCMCstring[i] -> setTEMPERATURE(bestTemperature * pow(base,i));
#endif
  }
#endif //PARALLEL_TEMPERING
#endif
#endif
#endif
}
#else
#ifdef NO_COUPLING
void binaryPopulation::assignTemperature(){
  double bestTemperature = bestTEMPERATURE;
  double worstTemperature = worstTEMPERATURE;

#ifdef LINEAR_COOLING  
  double sliceT = (bestTemperature - worstTemperature) / (double)sizeMCMCstring;
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(worstTemperature + sliceT*(i+0.5));
#else
#ifdef LOGARITHMIC_COOLING
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(bestTemperature / log(i+2));
#else
#ifdef GEOMETRIC_COOLING
  double base = exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)sizeMCMCstring));
  //cout << "base =" <<base << "\n";
  for(unsigned long i = 0; i < sizeMCMCstring; i++){
    MCMCstring[i] -> setTEMPERATURE(bestTemperature * pow(base,i));
    //cout << "i =(" <<i << "," << bestTemperature * pow(base,i) <<")\n";
  }
#else
#ifdef NO_COOLING
  for(unsigned long i = 0; i < sizeMCMCstring; i++){
#ifdef MULTIPLE_RUNS
	MCMCstring[i] -> setTEMPERATURE(bestTemperature);
#else
	MCMCstring[i] -> setTEMPERATURE(DEFAULT_TEMPERATURE);
#endif
  }
#endif
#endif
#endif
#endif
}
#else
#ifdef RANDOM_COUPLING
void binaryPopulation::assignTemperature(){
  double bestTemperature = bestTEMPERATURE;
  double worstTemperature = worstTEMPERATURE;

#ifdef LINEAR_COOLING  
  double sliceT = (bestTemperature - worstTemperature) / (double)sizeMCMCstring;
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(worstTemperature + sliceT*(i+0.5));
#else
#ifdef LOGARITHMIC_COOLING
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(bestTemperature / log(i+2));
#else
#ifdef GEOMETRIC_COOLING
  double base = 
      exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)sizeMCMCstring));
  //cout << "base =" <<base << "\n";
  for(unsigned long i = 0; i < sizeMCMCstring; i++){
    MCMCstring[i] -> setTEMPERATURE(bestTemperature * pow(base,i));
    //cout << "i =(" <<i << "," << bestTemperature * pow(base,i) <<")\n";
  }
#else
#ifdef NO_COOLING
    for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(DEFAULT_TEMPERATURE);
#endif
#endif
#endif
#endif
}
#else 
#ifdef POPULATION_MCMC
void binaryPopulation::assignTemperature(){
  double bestTemperature = bestTEMPERATURE;
  double worstTemperature = worstTEMPERATURE;

#ifdef LINEAR_COOLING  
  double sliceT = (bestTemperature - worstTemperature) / (double)sizeMCMCstring;
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(worstTemperature + sliceT*(i+0.5));
#else
#ifdef LOGARITHMIC_COOLING
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(bestTemperature / log(i+2));
#else
#ifdef GEOMETRIC_COOLING
  double base = exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)sizeMCMCstring));
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(bestTemperature * pow(base,i));
#else
#ifdef NO_COOLING
    for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(DEFAULT_TEMPERATURE);
#endif
#endif
#endif
#endif
}
#else 
#ifdef  WITHOUT_MEM
//procedure to assign to each state a temperature
void binaryPopulation::assignTemperature()
{
  unsigned long * tempTemp = new unsigned long[sizeMCMCstring];
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    tempTemp[i] = i;
	
  //sort the vector after the values of the fitness
  //bulesort
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    for(unsigned long j =0; j<sizeMCMCstring; j++) 	
      if(i!=j) {
	double temp_i = MCMCstring[tempTemp[i]] -> getValue();
	double temp_j = MCMCstring[tempTemp[j]] -> getValue();
	if(temp_i < temp_j && i < j)
	  {
	    unsigned long temp = tempTemp[i];
	    tempTemp[i] = tempTemp[j];
	    tempTemp[j] = temp;
	  }
	
	if(temp_j < temp_i && j < i)
	  {
	    unsigned long temp = tempTemp[i];
	    tempTemp[i] = tempTemp[j];
	    tempTemp[j] = temp;
	  }
      }
  
  //for(unsigned long i =0; i<sizeMCMCstring; i++)
  //  cout << "tempTemp[" <<i<<"]="<<tempTemp[i]<< " for " << MCMCstring[tempTemp[i]] -> getValue() << "\n";

  //the first and the last element are memorized if neccesary
  bestFitness = MCMCstring[tempTemp[0]] -> getValue();
  worstFitness = MCMCstring[tempTemp[sizeMCMCstring - 1]] -> getValue();	
  
#ifdef COMPLICAT_TEMPERING
  myHashtables2 = new Hashtable;
  
  unsigned long* temp = new unsigned long[sizeGENOME];
  unsigned long diferit = 1, ckeckI = -1;

  double** diferitValues = new double*[sizeMCMCstring];
  for(long i = 0; i<sizeMCMCstring; i++)
    diferitValues[i] = new double[4];

  for(long i =0; i<sizeMCMCstring; i++)
  { diferitValues[i][0] = 0; diferitValues[i][1] = 0; diferitValues[i][2] = 0; diferitValues[i][3] = 0; }

  diferitValues[diferit-1][0] = MCMCstring[tempTemp[0]]->getValue();
  diferitValues[diferit-1][1]++;
  diferitValues[diferit-1][2] = 1;
  diferitValues[diferit-1][3] = 0;
  MCMCstring[tempTemp[0]]->getGENOME(temp);
  myHashtables2 -> store_individual(temp,sizeGENOME);

  for(unsigned long i = 1; i < sizeMCMCstring; i++)
    if(MCMCstring[tempTemp[i-1]]->getValue()!= MCMCstring[tempTemp[i]]->getValue()) {
      diferitValues[diferit][0] = MCMCstring[tempTemp[i]]->getValue();
      diferitValues[diferit][1] = 1;
      diferitValues[diferit][2] = 1;
      diferitValues[diferit][3] = i;
      MCMCstring[tempTemp[i]]->getGENOME(temp);
      myHashtables2 -> store_individual(temp,sizeGENOME);
      //cout << diferit << " temp "<< temp << " hash " << myHashtables[diferit] -> check_individual(temp,sizeGENOME)<< "\n";
      diferit++;
    } else {
      MCMCstring[tempTemp[i]]->getGENOME(temp);
      if(myHashtables2 -> check_individual(temp,sizeGENOME) == 0){
	diferitValues[diferit-1][1]++;
	myHashtables2 -> store_individual(temp,sizeGENOME);
        //cout << "a "<< diferit-1 << " hash " << myHashtables[diferit-1] -> check_individual(temp,sizeGENOME)<< "\n";
      }
      else {
	myHashtables2 -> store_individual(temp,sizeGENOME);
	//cout << "b" << diferit-1 << " hash " << myHashtables[diferit-1] -> check_individual(temp,sizeGENOME)<< "\n";
      }
      double valueH = myHashtables2 -> check_individual(temp,sizeGENOME);
      if(valueH > diferitValues[diferit-1][2]) diferitValues[diferit-1][2] = valueH;
    }
  
  //for(unsigned long i = 1; i < diferit; i++)
  //cout << "i " << diferitValues[i][0] << " " << diferitValues[i][1] << " " << diferitValues[i][2] << " " <<diferitValues[i][3] << "\n";

#else
  unsigned long diferit = 1;
  for(unsigned long i = 1; i<sizeMCMCstring; i++)
    if(MCMCstring[tempTemp[i-1]]->getValue()!= MCMCstring[tempTemp[i]]->getValue()) diferit++;
#endif

  //temperatures are seted apriori
  bestTemperature = bestTEMPERATURE;
  worstTemperature = worstTEMPERATURE;

#ifdef LINEAR_COOLING
  double sliceF = (double)(bestFitness - worstFitness) / (double)diferit;
  
  if(bestFitness == worstFitness) sliceF = 1 / (double)sizeMCMCstring;

#ifdef PARALLEL_SIMULATED_ANNEALING_TEMPERING
  double baseTemp = exp(log((double)bestworstTEMPERATURE/(double)bestTEMPERATURE)*(1.0/(double)(burn_in+sampleSize*interSamplesize)));
  bestTemperature = bestTEMPERATURE * pow(baseTemp,runs);
#endif

  double sliceT = (bestTemperature - worstTemperature) / (double)diferit;

  if(bestFitness == worstFitness) sliceT =  (bestTemperature - worstTemperature)/ (double)sizeMCMCstring;
  
  //cout << "best fitness"<<bestFitness <<" worst fitness"<<worstFitness<< " diferit value "<<diferit;
  //cout << " first" << tempTemp[0] << " last" << tempTemp[sizeMCMCstring - 1];
  //cout << " slice T "<< sliceT << " slice F " << sliceF <<"\n";
  
  //assign the temperature for each state
  // this is a linear assignments (BT -WT)/sizePopulation for (BF - WF)/sizePopulation
  for(unsigned long i =0; i<sizeMCMCstring; i++){
    double tempX = MCMCstring[i] -> getValue();
    //find out in which interval is the given MCMC
    
    unsigned long j = 0;
    if(bestFitness != worstFitness){
      while(j < diferit){
	if(tempX < worstFitness + (j+1) * sliceF  &&  tempX >= worstFitness + j * sliceF) 
	  {
	    MCMCstring[i] ->setTEMPERATURE(worstTemperature + (diferit - j - 0.5)*sliceT);
	    //    cout << "assign " << i <<  "to "<< worstTemperature + (diferit - j-0.5)*sliceT << "\n";
	    j = sizeMCMCstring+1;
	  }	
	j++;	
      }
      if(j == diferit) {
	MCMCstring[i] -> setTEMPERATURE(worstTemperature + (0.5)*sliceT);
	//cout << "assign " << i <<  "to "<< worstTemperature + (0.5)*sliceT << "\n";
      } 
    }
    else {
      MCMCstring[i] ->setTEMPERATURE(worstTemperature + (i + 0.5)*sliceT);
      //cout << "equal assign " << i <<  "to "<< worstTemperature + (i+0.5)*sliceT << "\n";
    }  
  }	
#else
#ifdef GEOMETRIC_COOLING

#ifdef COMPLICAT_TEMPERING
  double diferitValueTotal = 0;
  for(unsigned long i = 0; i < diferit; i++)
    diferitValueTotal += diferitValues[i][2];
  //diferitValueTotal --;

  double sliceF;
  sliceF = (double)(bestFitness - worstFitness) / (double)diferitValueTotal;
  if(bestFitness == worstFitness) sliceF = 1 / (double)sizeMCMCstring;
#else 
  double sliceF;
  sliceF = (double)(bestFitness - worstFitness) / (double)diferit;
  if(bestFitness == worstFitness) sliceF = 1 / (double)sizeMCMCstring;
#endif

#ifdef PARALLEL_SIMULATED_ANNEALING_TEMPERING
  double baseTemp = exp(log((double)bestworstTEMPERATURE/(double)bestTEMPERATURE)*(1.0/(double)(burn_in+sampleSize*interSamplesize)));
  bestTemperature = bestTEMPERATURE * pow(baseTemp,generation);
#endif

#ifdef COMPLICAT_TEMPERING  
  double base = exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)diferitValueTotal));
#else
  double base = exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)diferit));
#endif

  if(bestFitness == worstFitness) base = exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)sizeMCMCstring));

  /*if(generation % 10 == 9){
    cout << "best fitness"<<bestFitness <<" worst fitness"<<worstFitness<< " diferit value "<<diferit;
    cout << " first MCMC= "<< tempTemp[0] << " last MCMC=" << tempTemp[sizeMCMCstring - 1];
    cout << " base "<< base << " slice F " << sliceF ;
    cout << " bestTemperature " << bestTemperature << " worstTemperature "<<worstTemperature
	 << " diferitTotal " <<diferitValueTotal << "\n";
    //cout << " generation " <<generation << " baseTemp "<< baseTemp;
    //cout << " bestworstTEMPERATURE"<<bestworstTEMPERATURE<<" bestTEMPERATURE"<< bestTEMPERATURE<<"\n"; 
    }*/
  
  for(unsigned long i =0; i<sizeMCMCstring; i++){
    double tempX = MCMCstring[i] -> getValue();
    
#ifdef COMPLICAT_TEMPERING   
    MCMCstring[i]->getGENOME(temp);
    unsigned long positionDif = -1, /*position = 0,*/ positionValDif = 0;
    
    for(unsigned long j = 0; j < diferit; j++)
      if(tempX == diferitValues[j][0]) {
	positionDif = j;
	break;
      } 
    
    for(unsigned long j = 0; j < positionDif; j++)
      positionValDif += diferitValues[j][2];
#else
    unsigned long j = 0;
#endif

    if(bestFitness != worstFitness){
#ifdef COMPLICAT_TEMPERING      
      double valueRamas = myHashtables2->check_individual(temp,sizeGENOME);

#ifndef FOARTE_COMPLICAT      
      if(valueRamas < 1)
	MCMCstring[i] ->setTEMPERATURE(bestTemperature*pow(base,diferitValueTotal - positionValDif+1));
      else {
	//double worstTemperature_i = bestTemperature*pow(base,diferitValueTotal - diferitValues[positionDif][3]);
	//double base_i = exp(log((double)worstTemperature_i/(double)bestTemperature)*(1.0/(double)diferitValues[positionDif][2]));
	double tempTempT = bestTemperature*pow(base,diferitValueTotal - valueRamas + 1 - positionValDif);
	MCMCstring[i] ->setTEMPERATURE(bestTemperature*pow(base,diferitValueTotal - valueRamas - positionValDif));
	//MCMCstring[i] ->setTEMPERATURE(bestTemperature*pow(base,positionValDif + valueRamas - 2));
      }
#else
      if(valueRamas == 1)
	MCMCstring[i] ->setTEMPERATURE(bestTemperature*pow(base,sizeMCMCstring - diferitValues[positionDif][3]));
      else {
	double worstTemperature_i = bestTemperature*pow(base,sizeMCMCstring - diferitValues[positionDif][3]);
	double base_i = exp(log((double)worstTemperature_i/(double)bestTemperature)*(1.0/(double)diferitValues[positionDif][2]));
	MCMCstring[i] ->setTEMPERATURE(bestTemperature*pow(base,sizeMCMCstring - valueRamas + 1 - diferitValues[positionDif][3]));
      }
#endif
      if(valueRamas > 0)
       myHashtables2->delete_individual(temp,sizeGENOME);

      /*if(generation % 10 == 9){
	valueRamas = myHashtables[positionDif]->check_individual(temp,sizeGENOME);
	cout << "after " << valueRamas;
	cout << "\n";
	}*/
#else
      while(j < diferit){
	if(tempX < worstFitness + (j+1) * sliceF  &&  tempX >= worstFitness + j * sliceF) 
	  //if(tempX >= bestFitness * pow(sliceF,diferit - j+1)   &&  tempX < bestFitness * pow(sliceF,diferit - j)) 
	  {
	    MCMCstring[i] ->setTEMPERATURE(bestTemperature*pow(base,j));
	    j = sizeMCMCstring+1;
	  } 
	//else 
	//cout << "tempX" << tempX << "best " <<bestFitness * pow(sliceF,diferit - j+1) << "worst " <<bestFitness * pow(sliceF,diferit - j) << "\n";
	j++;	
      }
      if(j == diferit) {
	MCMCstring[i] -> setTEMPERATURE(bestTemperature * pow(base,diferit));
	//cout << "assign " << i <<  "to "<< bestTemperature * pow(base,diferit) << " j= "<< j << "for temperature "<< tempX<<"\n";
      } 
#endif      
    }
    else {
      MCMCstring[i] ->setTEMPERATURE(bestTemperature * pow(base,i*sliceF));
      //cout << "equal assign " << i <<  "to "<< bestTemperature * pow(base,i*sliceF) << "for temperature "<< tempX<<"\n";
    }  
  }
  /*if(generation % 10 == 9){
    char t;
    cin >> t;
    }*/
#else
#ifdef NO_COOLING
  for(unsigned long i =0; i<sizeMCMCstring; i++){
    MCMCstring[i] ->setTEMPERATURE(DEFAULT_TEMPERATURE);
  }
#endif
#endif
#endif
#ifdef COMPLICAT_TEMPERING
  myHashtables2 -> free_table(myHashtables2); 
  //delete myHashtables2;
  delete[] diferitValues;
  delete[] temp;
#endif
}
#else
//parallel tempering
#ifdef NORMAL_TEMPERING
//procedure to assign to each state a temperature
void binaryPopulation::assignTemperature(){
  double bestTemperature = bestTEMPERATURE;
  double worstTemperature = worstTEMPERATURE;
  
#ifdef LINEAR_COOLING  
  double sliceT = (bestTemperature - worstTemperature) / (double)sizeMCMCstring;
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(worstTemperature + sliceT*(i+0.5));
#else
#ifdef LOGARITHMIC_COOLING
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(bestTemperature / log(i+2));
#else
#ifdef GEOMETRIC_COOLING
  double base = exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)sizeMCMCstring));
  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    MCMCstring[i] -> setTEMPERATURE(bestTemperature * pow(base,i));
#else
#ifdef NO_COOLING
    for(unsigned long i = 0; i < sizeMCMCstring; i++)
      MCMCstring[i] -> setTEMPERATURE(DEFAULT_TEMPERATURE);
#endif
#endif
#endif
#endif
}
#else
#ifdef PARALLEL_SIMULATED_ANNEALING_TEMPERING
//the information about temperature is stored in bestTemperature parameter
void binaryPopulation::assignTemperature()
{
  double bestTemperature = bestTEMPERATURE;
  double worstTemperature = worstTEMPERATURE;

  //linear schedual
#ifdef LINEAR_COOLING
  double sliceT = (bestTemperature - worstTemperature) / (double) (burn_in + sampleSize*interSamplesize);
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    MCMCstring[i] ->setTEMPERATURE(bestTemperature - generation*sliceT);
  //cout << "Temperature at "<<generation <<" is " << (bestTemperature - generation*sliceT) << "\n";
#else
#ifdef LOGARITHMIC_COOLING
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    MCMCstring[i] ->setTEMPERATURE(bestTemperature / log(generation + 2));
#else
#ifdef GEOMETRIC_COOLING
  double base = exp(log((double)worstTemperature/(double)bestTemperature)*(1.0/(double)(burn_in + sampleSize*interSamplesize)));
  //  cout << "base " << base;
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    if(generation == 0) MCMCstring[i] ->setTEMPERATURE(bestTemperature);
      else MCMCstring[i] ->setTEMPERATURE(bestTemperature * pow(base,generation));
#endif
#endif
#endif
}
#else
//procedure to assign to each state a temperature
void binaryPopulation::assignTemperature(){
  unsigned long * tempTemp = new unsigned long[sizeMCMCstring];
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    tempTemp[i] = i;

  //sort the vector after the values of the fitness
  //bulesort
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    for(unsigned long j =0; j<sizeMCMCstring; j++) 	
      if(i!=j) {
	double temp_i = MCMCstring[tempTemp[i]] -> getValue();
	double temp_j = MCMCstring[tempTemp[j]] -> getValue();
	if(temp_i < temp_j && i < j)
	  {
	    unsigned long temp = tempTemp[i];
	    tempTemp[i] = tempTemp[j];
	    tempTemp[j] = temp;
	  }
	
	if(temp_j < temp_i && j < i)
	  {
	    unsigned long temp = tempTemp[i];
	    tempTemp[i] = tempTemp[j];
	    tempTemp[j] = temp;
	  }
      }
  
  //the first and the last element are memorized if neccesary
  if(generation == 0) {
    bestFitness = MCMCstring[tempTemp[0]] -> getValue();
    worstFitness = MCMCstring[tempTemp[sizeMCMCstring - 1]] -> getValue();	
  } else {
    if(bestFitness < MCMCstring[tempTemp[0]] -> getValue())
      bestFitness = MCMCstring[tempTemp[0]] -> getValue();
    if(worstFitness > MCMCstring[tempTemp[sizeMCMCstring - 1]] -> getValue())
      worstFitness = MCMCstring[tempTemp[sizeMCMCstring - 1]] -> getValue();	
  }
  unsigned long diferit = 1;
  for(unsigned long i = 1; i<sizeMCMCstring; i++)
    if(MCMCstring[tempTemp[i-1]]->getValue()!= MCMCstring[tempTemp[i]]->getValue()) diferit++;

  double sliceF = (double)(bestFitness - worstFitness) / (double)diferit;
  
  if(bestFitness == worstFitness) sliceF = 1 / (double)sizeMCMCstring;
  
  //temperatures are seted apriori
  bestTemperature = bestTEMPERATURE;
  worstTemperature = worstTEMPERATURE;
  double sliceT = (bestTemperature - worstTemperature) / (double)diferit;

  if(bestFitness == worstFitness) sliceT =  (bestTemperature - worstTemperature)/ (double)sizeMCMCstring;
    //cout << "slice T "<< sliceT << " slice F " << sliceF <<"\n";
  
  //assign the temperature for each state
  // this is a linear assignments (BT -WT)/sizePopulation for (BF - WF)/sizePopulation
  for(unsigned long i =0; i<sizeMCMCstring; i++)
    {
      double tempX = MCMCstring[i] -> getValue();
      //find out in which interval is the given MCMC
      
      unsigned long j = 0;
      
      if(bestFitness != worstFitness){
	while(j < diferit){
	  if(tempX < worstFitness + (j+1) * sliceF  &&  tempX >= worstFitness + j * sliceF) 
	    {
	      MCMCstring[i] ->setTEMPERATURE(worstTemperature + (diferit - j - 0.5)*sliceT);
	      //cout << "assign " << i <<  "to "<< worstTemperature + (diferit - j-0.5)*sliceT << "\n";
	      j = sizeMCMCstring+1;
	    }	
	  j++;	
	}
	if(j == diferit) {
	  MCMCstring[i] -> setTEMPERATURE(worstTemperature + (0.5)*sliceT);
	} 
      }
      else MCMCstring[i] ->setTEMPERATURE(worstTemperature + (i + 0.5)*sliceT);  
      
    }	

}
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

void binaryPopulation::randomCOUPLING(unsigned long** table){
  unsigned long tempSize = sizeMCMCstring / sizeCoupling;
  unsigned long tempTableR[sizeMCMCstring];
  unsigned long index = 0;
  unsigned long mytrue;

  //setSeed();

#ifdef DEBUG_
  cout << "random coupling is begin with" << tempSize <<" "<<sizeCoupling <<"\n";
#endif

  for(unsigned long i = 0; i<tempSize ; i++)
      for(unsigned long j = 0; j < sizeCoupling ; j++)
	{
	  mytrue = 0;
	  while(mytrue != 1){
	    //srand(time(0));
	    table[i][j] = genrand_int32() % sizeMCMCstring;
	    mytrue = 1;
	    if(index != 0)
	      for(unsigned long k =0; k < index; k++)
		if(table[i][j] == tempTableR[k]) {
		  mytrue = 0;
		  //cout << "reject " <<table[i][j] << "cu " <<i <<","<<j<<"\n";
		}
	  }
	  tempTableR[index++] = table[i][j];
	}
  
#ifdef DEBUG_
  cout << "random coupling is finished \n";
  for(unsigned long i= 0; i<tempSize; i++){
    for(unsigned long j = 0; j<sizeCoupling; j++)
      cout << table[i][j] << " ";
    cout << "\n";
  }
#endif
}

void binaryPopulation::neighboursCOUPLING(unsigned long** table){
  unsigned long tempSize = sizeMCMCstring / sizeCoupling;
  unsigned long tempTableR[sizeMCMCstring];
  unsigned long index = 0;
  unsigned long mytrue;

  for(unsigned long i = 0; i<tempSize ; i++){
	  mytrue = 0;
	  while(mytrue != 1){
	    table[i][0] = genrand_int32() % sizeMCMCstring;
	    mytrue = 1;
	    if(index != 0)
	    for(unsigned long k =0; k < index; k++)
	      if(table[i][0] == tempTableR[k]) {
		    mytrue = 0;
		    //cout << "reject " <<table[i][j] << "cu " <<i <<","<<j<<"\n";
		  }
	  }
	  tempTableR[index++] = table[i][0];

	  if(table[i][0] == 0) table[i][1] = 1;
	   else if(table[i][0] == sizeMCMCstring - 1) table[i][1] = sizeMCMCstring - 2;
	     else if(genrand_int32() % 2 == 0) table[i][1] = table[i][0] + 1;
	       else table[i][1] = table[i][0] - 1;
    
  }
}

/*ATENTIE-Dimensiunea tabelei nu este fixa
returneaza dimensiunea tabelei*/
double binaryPopulation::neighboursCOUPLING_NO_OVERLAP(unsigned long** table){
  unsigned long tempTableR[sizeMCMCstring];
  unsigned long index = 0;
  unsigned long mytrue;
  unsigned long posible = 1;
  unsigned long dimension = 0;
  unsigned long vecinis, vecinid;
 
  if(sizeMCMCstring < sizeCoupling) return dimension;

  while(posible){
    mytrue = 0;
    while(mytrue != 1){
      table[dimension+1][0] = genrand_int32() % sizeMCMCstring;
      mytrue = 1;
      vecinis = 0; vecinid = 0;
      if(table[dimension + 1][0] == 0) vecinis++;
      if(table[dimension + 1][0] == sizeMCMCstring - 1) vecinid++;
      if(index != 0)
	for(unsigned long k =0; k < index; k++){
	      if(table[dimension+1][0] == tempTableR[k]) {
		mytrue = 0;
		//cout << "reject " <<table[i][j] << "cu " <<i <<","<<j<<"\n";
	      }
	      if(table[dimension + 1][0] - 1 == tempTableR[k])
		vecinis++;
	      if(table[dimension + 1][0] + 1 == tempTableR[k])
		vecinid++;
	}
      if(vecinis > 0 && vecinid > 0) mytrue = 0;
    }

    tempTableR[index++] = table[dimension + 1][0];

    if(table[dimension + 1][0] == 0) table[dimension + 1][1] = 1;
    else if(table[dimension + 1][0] == sizeMCMCstring - 1) table[dimension + 1][1] = sizeMCMCstring - 2;
    else if(vecinis == 1) table[dimension + 1][1] = table[dimension + 1][0] + 1;
    else if(vecinid == 1) table[dimension + 1][1] = table[dimension + 1][0] - 1;
    else if(genrand_int32() % 2 == 0) table[dimension + 1][1] = table[dimension + 1][0] + 1;
    else table[dimension + 1][1] = table[dimension + 1][0] - 1;

    tempTableR[index++] = table[dimension + 1][1];
    
    if(dimension + 1 == sizeMCMCstring/2 - 1) posible = 0;

    //exista cel putin un individ faca cuplaj si fara vecini
    unsigned long i = 0; 
    posible = 0;
    while(i < sizeMCMCstring && posible == 0){
      unsigned k = 0; 
      
      posible = 1;
      while(k < index && posible != 0){
	if(tempTableR[k] == i) posible = 0;
	k++;
      }
      
      if(posible == 1){
	k = 0; vecinis = 0; vecinid = 0;
	if(i == 0) vecinis = 1;
	if(i == sizeMCMCstring - 1) vecinid = 1;
	while(k < index && !(vecinis == 1 && vecinid == 1)){
	  if(tempTableR[k] == i - 1) vecinis = 1;
	  if(tempTableR[k] == i + 1) vecinid = 1;
	  k++;
	}
	if(vecinis ==1 && vecinid == 1)
	  posible = 0;
      }
      i++;
    }
    dimension++;
  }
  return dimension + 1;
}

void binaryPopulation::proposedGeneration(unsigned long** table, unsigned long dimension){
#ifdef BINARY_RECOM
#ifdef METROPOLIS_ALGORITHM
  init();

#ifndef MIXTURE_REALS
#ifdef RECOMBINATION_TEST
#ifdef MIXTURE
  if(MIXTURE <= genrand_real1()){
      recombination(table, dimension, RECOMBINATION_TEST_PARAM);
      mutation(table,dimension);
      //recombination(table, dimension, UNIFORM_RECOMB);
      //cout << "recombinatie ";
  }
  else{
      mutation();
      //cout << "mutatie ";
  }
#else  //cout << "recomb test" << RECOMBINATION_TEST_PARAM << "\n";
  recombination(table, dimension, RECOMBINATION_TEST_PARAM);
  mutation(table,dimension);
#endif //MIXTURE
#else // recombination test
#ifdef CYCLE
  if(generation % 2 == 1){
      recombination(table, dimension, ONE_POINT_RECOMB);
      //mutation();
      //recombination(table, dimension, UNIFORM_RECOMB);
      //cout << "recombinatie ";
  }
  else{
      mutation();
      //cout << "mutatie ";
  }
#else //cycle
#ifdef MIXTURE
  if(MIXTURE <= genrand_real1()){
      recombination(table, dimension, ONE_POINT_RECOMB);
      //mutation();
      //recombination(table, dimension, UNIFORM_RECOMB);
      //cout << "recombinatie ";
  }
  else{
      mutation();
      //cout << "mutatie ";
  }
#else //mixture
#ifdef IsMutation
  mutation();
#endif
#ifdef IsRecombination
#ifndef POPULATION_MCMC
  //recombination(table, dimension, ONE_POINT_RECOMB);
  //recombination(table, dimension, N_POINT_RECOMB);
  //recombination(table, dimension, TWO_POINTS_RECOMB);
  recombination(table, dimension, UNIFORM_RECOMB);
  //recombination(table, dimension, UNIFORM_RECOMB_PROBABIL);
  //recombination(table, dimension, NO_RECOMB);
  //cout << "Recombination" << UNIFORM_RECOMB <<"\n";
#else
  //recombination(table, dimension, NO_RECOMB);
  //recombination(table, UNIFORM_RECOMB_PROBABIL);
#endif
#endif //IsRecombination
#endif //Mixture
#endif //Cycle
#endif //Recombination test
#else /////////////////////mixture reals and metropolis
 
 //cout << "verification that it is entering in the right place \n";

/////////%%%%%%%%%%%%%%%
#ifndef RECOMBINATION_TEST
#ifndef MIXTURE
#ifdef IsMutation
  //cout << "Selected mutation \n";
  mutation();
#endif

#ifdef IsRecombination
  //cout << "Selected recombination \n";
  recombination(table, dimension, UNIFORM_RECOMB);
#endif

#else // mixture

  if(MIXTURE <= genrand_real1()){
#ifdef IsRecombination
    recombination(table, dimension, UNIFORM_RECOMB);
#else
    mutation();
#endif
    mix_mutation = 0;
  }
  else{
    mutation();
    mix_mutation = 1;
  }
  
#endif //mixture

#else //recomb test param
  recombination(table, dimension, RECOMBINATION_TEST_PARAM);
  mutation(table,dimension);
#endif  //recomb test
////////%%%%%%

#endif //mixture reals with symmetry

#else //metropolis
#ifdef RECOMBINATION_TEST
  init();
  //cout << "recomb test" << RECOMBINATION_TEST_PARAM << "\n";
  recombination(table, dimension, RECOMBINATION_TEST_PARAM);
  proposalDistributionCalculation(table,dimension);
#else
#ifdef POPULATION_MCMC
  proposalDistributionCalculation();
  init();
  recombination(table, dimension, UNIFORM_RECOMB_PROBABIL);
#else
  // cout << "Propose generation \n";
  init();
  recombination(table, dimension, UNIFORM_RECOMB_PROBABIL);
  proposalDistributionCalculation(table,dimension);
#endif //POPULATION_MCMC
#endif //RECOMBINATION_TEST
#endif //METROPOLIS_ALGORITHM

#else //not binary recombination

  init();

#ifdef MIXTURE
  //cout << "\mixture \n";
  if(MIXTURE >= genrand_real1()){
    //cout << "recombination \t";
#ifdef IsRecombination
    //cout << " recombniae \n";
    recombination(table, dimension, PARENT_CENTRIC);
#ifdef METROPOLIS_HASTINGS_ALGORITHM
    proposalDistributionCalculation(table,dimension);
#endif //metropolis
    mix_mutation = 0;
#else
    //cout << " mutate after recombine \n"; 
    mutation();
    mix_mutation = 1;
#endif
  }
  else{
#ifdef IsMutation
    //cout << "mutation \n";
    mutation();
    mix_mutation = 1;
#else
    recombination(table, dimension, PARENT_CENTRIC);
#ifdef METROPOLIS_HASTINGS_ALGORITHM
    proposalDistributionCalculation(table,dimension);
#endif //metropolis
    mix_mutation = 0;
#endif
  }
#else //mixture
#ifdef CYCLE
  if(generation % 2 == 0){
#ifdef IsMutation
    mutation();
    mix_mutation = 1;
#else
#ifdef IsRecombination
      recombination(table, dimension, PARENT_CENTRIC);
#else
      cout << "not implemented \n"
	exit(1);
#endif 
#endif//mutation
      //mutation();
      //recombination(table, dimension, UNIFORM_RECOMB);
      //cout << "recombinatie ";
  }
  else{
#ifdef IsRecombination
      recombination(table, dimension, PARENT_CENTRIC);
      mix_mutation = 0;
#else
#ifdef IsMutation
      mutation();
#else
      cout << "not implemented \n";
      exit(1);
#endif
#endif //is recombination
      //cout << "mutatie ";
  }
#else //cycle
#ifdef SYMMETRIC_PROP

#ifdef IsRecombination
    recombination(table, dimension, PARENT_CENTRIC);
#endif

#ifdef IsMutation
  //cout << "Selected mutation \n";
  mutation();
#endif 

#else
    if(0.5 <= genrand_real1()){
#ifdef IsMutation
      //cout << "Selected mutation \n";
      mutation();
#endif 
      
#ifdef IsRecombination
      recombination(table, dimension, PARENT_CENTRIC);
#endif
    } else {
#ifdef IsRecombination
      recombination(table, dimension, PARENT_CENTRIC);
#endif
      
#ifdef IsMutation
      //cout << "Selected mutation \n";
      mutation();
#endif 
    }
#endif //symmetric proposal distribution
#endif //cycle
#endif // mixture

#endif // binary recomb
}

#ifndef METROPOLIS_ALGORITHM
#ifdef RECOMBINATION_TEST 
void binaryPopulation::proposalDistributionCalculation(unsigned long** table, int dimension){

  unsigned long tempSize = dimension;
  //int temp[sizeGENOME];

  int** recomb_table = new int*[sizeCoupling];
  int** prop_recomb_table = new int*[sizeCoupling];
  for(unsigned long i = 0; i<sizeCoupling; i++){
    recomb_table[i] = new int[sizeGENOME];
    prop_recomb_table[i] = new int[sizeGENOME];
  }

  for(unsigned long i = 0; i < sizeMCMCstring; i++){
      proposalDistribution[i] = 1;
      proposedProposalDistribution[i] = 1;
  }

  for(unsigned long i = 0; i<tempSize ; i++){
      // cout << "pair (";
      for(unsigned long j = 0; j<sizeCoupling; j++){
  	MCMCstring[table[i][j]] -> getGENOME(recomb_table[j]);
  	proposedMCMCstring[table[i][j]] -> getGENOME(prop_recomb_table[j]);
	//cout << table[i][j] << ",";
	//for(int k = 0; k < sizeGENOME; k++)
	//    cout << recomb_table[j][k];
	//cout << ",";
	//for(int k = 0; k < sizeGENOME; k++)
	//    cout << prop_recomb_table[j][k];
	//    cout << ",";
      }
     
#ifndef TWO_CHILDREN 
      //cout << ") \n Proposed distribution =(";
      if(MCMCstring[table[i][0]]->parentReplaced == 1){
	  proposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  proposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  //  cout << "parent repl. 0 " << proposalDistribution[table[i][0]] << ",(";
	  for(int k = 0; k < sizeGENOME; k++){
	      if(recomb_table[0][k] == recomb_table[1][k]){
		  if(prop_recomb_table[1][k] == recomb_table[0][k])
		      proposalDistribution[table[i][1]] *= (double)(sizeGENOME-1)/(double)sizeGENOME;
		  else proposalDistribution[table[i][1]] *= 1.0/(double)sizeGENOME;
	      } else {
		  if(prop_recomb_table[1][k] == recomb_table[1][k])
		      proposalDistribution[table[i][1]] *= 1 - mutation_multi/(double)sizeGENOME; 
		  else proposalDistribution[table[i][1]] *= mutation_multi/(double)sizeGENOME;
	      }
	      //  cout << proposalDistribution[table[i][1]] <<" ";
	  }
	  // cout << "=" << proposalDistribution[table[i][1]] << ")";
      } else {
	  proposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  proposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  for(int k = 0; k < sizeGENOME; k++)
	      if(recomb_table[0][k] == recomb_table[1][k]){
		  if(prop_recomb_table[0][k] == recomb_table[0][k])
		      proposalDistribution[table[i][0]] *=  (double)(sizeGENOME-1)/(double)sizeGENOME;
		  else proposalDistribution[table[i][0]] *= 1.0/(double)sizeGENOME;
	      } else {
		  if(prop_recomb_table[0][k] == recomb_table[0][k])
		      proposalDistribution[table[i][0]] *= 1.0 - mutation_multi/(double)sizeGENOME;
		  else   proposalDistribution[table[i][0]] *= mutation_multi/(double)sizeGENOME;
	      }
	  //cout << proposalDistribution[table[i][0]] << "," << proposalDistribution[table[i][1]];
      }

      // cout << ") proposed proposal distribution (";
      if(MCMCstring[table[i][0]]->parentReplaced == 1){
	  proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  for(int k = 0; k < sizeGENOME; k++)
	      if(prop_recomb_table[0][k] == prop_recomb_table[1][k]){
		  if(recomb_table[1][k] == prop_recomb_table[1][k])
		      proposedProposalDistribution[table[i][1]] *=  (double)(sizeGENOME-1)/(double)sizeGENOME;
		      else proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeGENOME;
		  } else {
		      if(recomb_table[1][k] == prop_recomb_table[1][k])
			  proposedProposalDistribution[table[i][1]] *= 1.0 - mutation_multi/(double)sizeGENOME;
		      else proposedProposalDistribution[table[i][1]] *= mutation_multi/(double)sizeGENOME;
		  }
	  //  cout << proposedProposalDistribution[table[i][0]] << "," << proposedProposalDistribution[table[i][1]];
      } else { 
	  proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  for(int k = 0; k < sizeGENOME; k++)
	      if(prop_recomb_table[0][k] == prop_recomb_table[1][k]){
		  if(recomb_table[0][k] == prop_recomb_table[0][k])
		      proposedProposalDistribution[table[i][0]] *= (double)(sizeGENOME-1)/(double)sizeGENOME;
		      else proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeGENOME;
		  } else {
		      if(recomb_table[0][k] == prop_recomb_table[0][k])
			  proposedProposalDistribution[table[i][0]] *= 1.0 - mutation_multi/(double)sizeGENOME;
		      else 
			  proposedProposalDistribution[table[i][0]] *= mutation_multi/(double)sizeGENOME;
		  }
	  //cout << proposedProposalDistribution[table[i][0]] << "," << proposedProposalDistribution[table[i][1]];
      }
      //cout << ")\n";
#else
      //cout << ") \n Proposed distribution =(";
      for(int k = 0; k < sizeGENOME; k++){
	  if(recomb_table[0][k] == recomb_table[1][k]){
	      if(prop_recomb_table[1][k] == recomb_table[1][k])
		  proposalDistribution[table[i][1]] *= (double)(sizeGENOME-1)/(double)sizeGENOME;
	      else proposalDistribution[table[i][1]] *= 1.0/(double)sizeGENOME;
	  } else {
	      if(prop_recomb_table[1][k] == recomb_table[1][k])
		  proposalDistribution[table[i][1]] *= 1.0 - mutation_multi/(double)sizeGENOME;
	      else 
		  proposalDistribution[table[i][1]] *= mutation_multi/(double)sizeGENOME;
	  }
	  //  cout << proposalDistribution[table[i][1]] <<" ";
      }
      // cout << "=" << proposalDistribution[table[i][1]] << ")";
      for(int k = 0; k < sizeGENOME; k++)
	  if(recomb_table[0][k] == recomb_table[1][k]){
	      if(prop_recomb_table[0][k] == recomb_table[0][k])
		  proposalDistribution[table[i][0]] *=  (double)(sizeGENOME-1)/(double)sizeGENOME;
	      else proposalDistribution[table[i][0]] *= 1.0/(double)sizeGENOME;
	  } else {
	      if(prop_recomb_table[0][k] == recomb_table[0][k]) 
		  proposalDistribution[table[i][0]] *= 1.0 - mutation_multi/(double)sizeGENOME;
	      else 
		  proposalDistribution[table[i][0]] *= mutation_multi/(double)sizeGENOME;
	  }
      //cout << proposalDistribution[table[i][0]] << "," << proposalDistribution[table[i][1]];
  

      // cout << ") proposed proposal distribution (";
      for(int k = 0; k < sizeGENOME; k++)
	  if(prop_recomb_table[0][k] == prop_recomb_table[1][k]){
	      if(recomb_table[1][k] == prop_recomb_table[1][k])
		  proposedProposalDistribution[table[i][1]] *=  (double)(sizeGENOME-1)/(double)sizeGENOME;
	      else proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeGENOME;
	  } else {
	      if(recomb_table[1][k] == prop_recomb_table[1][k])
		  proposedProposalDistribution[table[i][1]] *= 1.0 - mutation_multi/(double)sizeGENOME;
	      else 
		  proposedProposalDistribution[table[i][1]] *= mutation_multi/(double)sizeGENOME;
	  }
      //  cout << proposedProposalDistribution[table[i][0]] << "," << proposedProposalDistribution[table[i][1]];
      for(int k = 0; k < sizeGENOME; k++)
	  if(prop_recomb_table[0][k] == prop_recomb_table[1][k]){
	      if(recomb_table[0][k] == prop_recomb_table[0][k])
		      proposedProposalDistribution[table[i][0]] *= (double)(sizeGENOME-1)/(double)sizeGENOME;
	      else proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeGENOME;
	  } else {
	      if(recomb_table[0][k] == prop_recomb_table[0][k])
		  proposedProposalDistribution[table[i][0]] *= 1.0 - mutation_multi/(double)sizeGENOME;
	      else 
		  proposedProposalDistribution[table[i][0]] *= mutation_multi/(double)sizeGENOME;
	  }
      //cout << proposedProposalDistribution[table[i][0]] << "," << proposedProposalDistribution[table[i][1]];
      
      //cout << ")\n";
#endif // TWO_CHILDREN
  }
  delete[] recomb_table;
  delete[] prop_recomb_table;
}
#else
#ifdef POPULATION_MCMC
void binaryPopulation::proposalDistributionCalculation(){
#ifdef INDEPENDENT_SAMPLER
  for(unsigned long i = 0; i < sizeGENOME; i++)
    for(unsigned long j = 0; j < 2; j++)
      proposalDistribution[i][j] = 1;

  int temp[sizeGENOME];
  for(unsigned long j = 0; j < sizeMCMCstring; j++){
      MCMCstring[j]->getGENOME(temp);
      for(unsigned long i = 0; i < sizeGENOME; i++)
	  proposalDistribution[i][temp[i]] ++;
  }
  
  for(unsigned long i = 0; i < sizeGENOME; i++){
      for(unsigned long j = 0; j < 2; j++)
#ifdef EXPLOITATION_INDEPENDENT_SAMPLER
	  proposalDistribution[i][j] /= (sizeMCMCstring + 2);
#else
#ifdef METROPOLIS_HASTINGS_ALGORITHM
#ifdef POPULATION_MCMC
// pentru operatorul same loci mutation and recombination
      proposalDistribution[i][j] = (double)(sizeMCMCstring+2-proposalDistribution[i][j])/(sizeMCMCstring+2);
#else
      //trebuie sa stiu care este partenerul
      
#endif
#else
//#ifdef METROPOLIS_ALGORITHM
    proposalDistribution[i][j] = 1;
#endif
#endif
    }
#else
    for(unsigned long i = 0; i < sizeGENOME; i++)
      for(unsigned long j = 0; j < 2; j++)
	  proposalDistribution[i][j] = 1;    
#endif //independent sampler
}

void binaryPopulation::CalculateNewProposalDistribution(unsigned long kint){
#ifdef INDEPENDENT_SAMPLER
  for(unsigned long i = 0; i < sizeGENOME; i++)
    for(unsigned long j = 0; j < 2; j++)
      PropVector[i][j] = 1;
  
  int temp[sizeGENOME];
    for(unsigned long j = 0; j < sizeMCMCstring ; j++)
      if (kint != j){
	MCMCstring[j]->getGENOME(temp);
	for(unsigned long i = 0; i < sizeGENOME; i++)
	  PropVector[i][temp[i]] ++;
      } else {
	proposedMCMCstring[j]->getGENOME(temp);
	for(unsigned long i = 0; i < sizeGENOME; i++)
	  PropVector[i][temp[i]] ++;
      }
    for(unsigned long i = 0; i < sizeGENOME; i++){
      for(unsigned long j = 0; j < 2; j++)
#ifdef EXPLOITATION_INDEPENDENT_SAMPLER
	PropVector[i][j] /= (sizeMCMCstring + 2);
#else
#ifdef METROPOLIS_HASTINGS_ALGORITHM
        PropVector[i][j] = (double)(sizeMCMCstring+2-PropVector[i][j])/(sizeMCMCstring+2);
#else
	PropVector[i][j] = 1;
#endif
#endif    
    }
#else 
    //   for(unsigned long i = 0; i < sizeGENOME; i++)
    //  for(unsigned long j = 0; j < 2; j++)
    //	  PropVector[i][j] = 1;
#endif
}

void binaryPopulation::CalculateNewProposalDistributionRecombination(unsigned long kint, unsigned long rint){
#ifdef INDEPENDENT_SAMPLER
  for(unsigned long i = 0; i < sizeGENOME; i++)
    for(unsigned long j = 0; j < 2; j++)
      PropVector[i][j] = 1;
  
  int temp[sizeGENOME];
    for(unsigned long j = 0; j < sizeMCMCstring ; j++)
      if (kint != j && rint != j){
	MCMCstring[j]->getGENOME(temp);
	for(unsigned long i = 0; i < sizeGENOME; i++)
	  PropVector[i][temp[i]] ++;
      } else {
	proposedMCMCstring[j]->getGENOME(temp);
	for(unsigned long i = 0; i < sizeGENOME; i++)
	  PropVector[i][temp[i]] ++;
      }
    for(unsigned long i = 0; i < sizeGENOME; i++){
      for(unsigned long j = 0; j < 2; j++)
#ifdef EXPLOITATION_INDEPENDENT_SAMPLER
	PropVector[i][j] /= (sizeMCMCstring + 2);
#else
#ifdef METROPOLIS_HASTINGS_ALGORITHM
        PropVector[i][j] = (double)(sizeMCMCstring+2-PropVector[i][j])/(sizeMCMCstring+2);
#else
	PropVector[i][j] = 1;
#endif
#endif   
    } 
#else // it is not independent sampler
    //for(unsigned long i = 0; i < sizeGENOME; i++){
    //  for(unsigned long j = 0; j < 2; j++)
//	  PropVector[i][j] = 1;
#endif
}
#else //POPULATION_MCMC

#ifdef BINARY_RECOMB
void binaryPopulation::proposalDistributionCalculation(unsigned long** table, int dimension){

  unsigned long tempSize = dimension;
  //int temp[sizeGENOME];

  int** recomb_table = new int*[sizeCoupling];
  int** prop_recomb_table = new int*[sizeCoupling];
  for(unsigned long i = 0; i<sizeCoupling; i++){
    recomb_table[i] = new int[sizeGENOME];
    prop_recomb_table[i] = new int[sizeGENOME];
  }

  for(unsigned long i = 0; i < sizeMCMCstring; i++){
      proposalDistribution[i] = 1;
      proposedProposalDistribution[i] = 1;
  }

  for(unsigned long i = 0; i<tempSize ; i++){
      // cout << "pair (";
      for(unsigned long j = 0; j<sizeCoupling; j++){
  	MCMCstring[table[i][j]] -> getGENOME(recomb_table[j]);
  	proposedMCMCstring[table[i][j]] -> getGENOME(prop_recomb_table[j]);
	//cout << table[i][j] << ",";
	//for(int k = 0; k < sizeGENOME; k++)
	//    cout << recomb_table[j][k];
	//cout << ",";
	//for(int k = 0; k < sizeGENOME; k++)
	//    cout << prop_recomb_table[j][k];
	//    cout << ",";
      }
     
#ifndef TWO_CHILDREN 
      //cout << ") \n Proposed distribution =(";
      if(MCMCstring[table[i][0]]->parentReplaced == 1){
	  proposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  proposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  //  cout << "parent repl. 0 " << proposalDistribution[table[i][0]] << ",(";
	  for(int k = 0; k < sizeGENOME; k++){
	      if(recomb_table[0][k] == recomb_table[1][k]){
		  if(prop_recomb_table[1][k] == recomb_table[0][k])
		      proposalDistribution[table[i][1]] *= (double)(sizeGENOME-1)/(double)sizeGENOME;
		  else proposalDistribution[table[i][1]] *= 1.0/(double)sizeGENOME;
	      } else proposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	      //  cout << proposalDistribution[table[i][1]] <<" ";
	  }
	  // cout << "=" << proposalDistribution[table[i][1]] << ")";
      } else {
	  proposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  proposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  for(int k = 0; k < sizeGENOME; k++)
	      if(recomb_table[0][k] == recomb_table[1][k]){
		  if(prop_recomb_table[0][k] == recomb_table[0][k])
		      proposalDistribution[table[i][0]] *=  (double)(sizeGENOME-1)/(double)sizeGENOME;
		  else proposalDistribution[table[i][0]] *= 1.0/(double)sizeGENOME;
	      } else proposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  //cout << proposalDistribution[table[i][0]] << "," << proposalDistribution[table[i][1]];
      }

      // cout << ") proposed proposal distribution (";
      if(MCMCstring[table[i][0]]->parentReplaced == 1){
	  proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  for(int k = 0; k < sizeGENOME; k++)
	      if(prop_recomb_table[0][k] == prop_recomb_table[1][k]){
		  if(recomb_table[1][k] == prop_recomb_table[0][k])
		      proposedProposalDistribution[table[i][1]] *=  (double)(sizeGENOME-1)/(double)sizeGENOME;
		      else proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeGENOME;
		  } else proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  //  cout << proposedProposalDistribution[table[i][0]] << "," << proposedProposalDistribution[table[i][1]];
      } else { 
	  proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  for(int k = 0; k < sizeGENOME; k++)
	      if(prop_recomb_table[0][k] == prop_recomb_table[1][k]){
		  if(recomb_table[0][k] == prop_recomb_table[0][k])
		      proposedProposalDistribution[table[i][0]] *= (double)(sizeGENOME-1)/(double)sizeGENOME;
		      else proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeGENOME;
		  } else proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
	  //cout << proposedProposalDistribution[table[i][0]] << "," << proposedProposalDistribution[table[i][1]];
      }
      //cout << ")\n";
#else
      //cout << ") \n Proposed distribution =(";
      for(int k = 0; k < sizeGENOME; k++){
	  if(recomb_table[0][k] == recomb_table[1][k]){
	      if(prop_recomb_table[1][k] == recomb_table[0][k])
		  proposalDistribution[table[i][1]] *= (double)(sizeGENOME-1)/(double)sizeGENOME;
	      else proposalDistribution[table[i][1]] *= 1.0/(double)sizeGENOME;
	  } else proposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
	  //  cout << proposalDistribution[table[i][1]] <<" ";
      }
      // cout << "=" << proposalDistribution[table[i][1]] << ")";
      for(int k = 0; k < sizeGENOME; k++)
	  if(recomb_table[0][k] == recomb_table[1][k]){
	      if(prop_recomb_table[0][k] == recomb_table[0][k])
		  proposalDistribution[table[i][0]] *=  (double)(sizeGENOME-1)/(double)sizeGENOME;
	      else proposalDistribution[table[i][0]] *= 1.0/(double)sizeGENOME;
	  } else proposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
      //cout << proposalDistribution[table[i][0]] << "," << proposalDistribution[table[i][1]];
  

      // cout << ") proposed proposal distribution (";
      for(int k = 0; k < sizeGENOME; k++)
	  if(prop_recomb_table[0][k] == prop_recomb_table[1][k]){
	      if(recomb_table[1][k] == prop_recomb_table[0][k])
		  proposedProposalDistribution[table[i][1]] *=  (double)(sizeGENOME-1)/(double)sizeGENOME;
	      else proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeGENOME;
	  } else proposedProposalDistribution[table[i][1]] *= 1.0/(double)sizeCoupling;
      //  cout << proposedProposalDistribution[table[i][0]] << "," << proposedProposalDistribution[table[i][1]];
      for(int k = 0; k < sizeGENOME; k++)
	  if(prop_recomb_table[0][k] == prop_recomb_table[1][k]){
	      if(recomb_table[0][k] == prop_recomb_table[0][k])
		      proposedProposalDistribution[table[i][0]] *= (double)(sizeGENOME-1)/(double)sizeGENOME;
	      else proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeGENOME;
	  } else proposedProposalDistribution[table[i][0]] *= 1.0/(double)sizeCoupling;
      //cout << proposedProposalDistribution[table[i][0]] << "," << proposedProposalDistribution[table[i][1]];
      
      //cout << ")\n";
#endif // TWO_CHILDREN
  }
  delete[] recomb_table;
  delete[] prop_recomb_table;
}
#else //not binary_recomb
#ifdef wPCTX
void binaryPopulation::proposalDistributionCalculation(unsigned long** table, unsigned long dimension){

 unsigned long tempSize = dimension;
  //int temp[sizeGENOME];

  double** recomb_table = new double*[sizeCoupling];
  double** intermediate = new double*[sizeGENOME];
  double** prop_recomb_table = new double*[sizeCoupling];
  for(unsigned long i = 0; i<sizeCoupling; i++){
    recomb_table[i] = new double[sizeGENOME];
    prop_recomb_table[i] = new double[sizeGENOME];
  }
  for(unsigned long i = 0; i<sizeGENOME; i++){
    intermediate[i] = new double[sizeGENOME];
  }

  double* center = new double[sizeGENOME];
  double* fit = new double[sizeCoupling];
  double* prop_fit = new double[sizeCoupling];
  
  for(unsigned long i = 0; i < sizeMCMCstring; i++){
      proposalDistribution[i] = 1;
      proposedProposalDistribution[i] = 1;
  }

  //cout << "propocal dsitrib comp \n";

  for(unsigned long i = 0; i<tempSize ; i++){
      // cout << "pair (";
      for(unsigned long j = 0; j<sizeCoupling; j++){
	  MCMCstring[table[i][j]] -> getGENOME(recomb_table[j]);
	  proposedMCMCstring[table[i][j]] -> getGENOME(prop_recomb_table[j]);

	  fit[j] = MCMCstring[table[i][j]]->getValue();
	  prop_fit[j] = proposedMCMCstring[table[i][j]]->fitness();
	  //cout << table[i][j] << ",";
	  //for(int k = 0; k < sizeGENOME; k++)
	  //    cout << recomb_table[j][k];
	  //cout << ",";
	  //for(int k = 0; k < sizeGENOME; k++)
	  //    cout << prop_recomb_table[j][k];
	  //    cout << ",";
      }

//translation
      //if(MCMCstring[table[i][0]] -> mix_mixture == 0){
	  //compute the translation probability based on the fitness of the proposed and of the originals
	
	//1.0 comute the centers of the originals
	for(unsigned long j = 0; j < sizeGENOME; j++){
	  center[j] = 0;
	  double jos = 0;
	  for(unsigned long k = 0; k < sizeCoupling; k++){
	    center[j] += recomb_table[k][j] * fit[k];
	    jos += fit[k];
	  }
	  center[j] /= jos;
	}
	  
	if(sizeGENOME == 2){
	//for(unsigned long j = 0; j < sizeCoupling; j++){
	//for(unsigned long k = 1; k < sizeGENOME; k++ ){
	//learn the first intermediate which is in the same direction as the center-point
	//  cout << "start\n";
	double m = (recomb_table[0][1] - center[1])/(recomb_table[0][0] - center[0]);
#ifdef EUCLIDIAN
	double euclid = sqrt(pow(recomb_table[0][1] - center[1],2.0)+ pow(recomb_table[0][0] - center[0],2.0));
#endif 

#ifndef TRANSLATION_COMPLEX

#ifndef EUCLIDIAN
	proposalDistribution[table[i][0]] *= exp(-pow(prop_recomb_table[0][0] - recomb_table[0][0],2) / (2.0 * pow((recomb_table[0][0] - center[0])*sigma_transl,2.0)))/(sqrt(pow(recomb_table[0][0] - center[0],2.0))*sigma_transl*sqrt(2 * pi)) ;
	proposalDistribution[table[i][0]] *= exp(-pow(prop_recomb_table[0][1] - recomb_table[0][1],2) / (2.0 * pow((recomb_table[0][1] - center[1])*sigma_transl,2.0)))/(sqrt(pow(recomb_table[0][1] - center[1],2.0)) * sigma_transl * sqrt(2 * pi));
#else
	//cout << "proposal(" << proposalDistribution[table[i][0]] << ",";
	proposalDistribution[table[i][0]] *= exp(-pow(prop_recomb_table[0][0] - recomb_table[0][0],2) / (2.0 * pow((recomb_table[0][0] - center[0])*sigma_transl*euclid,2.0)))/(sqrt(pow(recomb_table[0][0] - center[0],2.0))*sigma_transl*euclid*sqrt(2 * pi));
	//cout << proposalDistribution[table[i][0]] << ",";
	proposalDistribution[table[i][0]] *= exp(-pow(prop_recomb_table[0][1] - recomb_table[0][1],2) / (2.0 * pow((recomb_table[0][1] - center[1])*sigma_transl*euclid,2.0)))/(sqrt(pow(recomb_table[0][1] - center[1],2.0)) * sigma_transl * euclid * sqrt(2 * pi));
	//cout << proposalDistribution[table[i][0]] << ") \t";
#endif
	
#else
	  // in the dirrection x - g
	//	  if(k == 1) {
	intermediate[0][0] = m / (pow(m,2) + 1) * (prop_recomb_table[0][1] + 
		   1.0 /m * prop_recomb_table[0][0] + m * recomb_table[0][0] - recomb_table[0][1]);

	intermediate[0][1] = recomb_table[0][1] + m * (intermediate[0][0] - recomb_table[0][0]); 

	intermediate[1][0] = m / (pow(m,2) + 1) *(recomb_table[0][1] + 
		   1.0 / m * recomb_table[0][0] - prop_recomb_table[0][1] + m * prop_recomb_table[0][0]);

	intermediate[1][1] = prop_recomb_table[0][1] + m * (intermediate[1][0] - prop_recomb_table[0][0]); 

#ifndef EUCLIDIAN
	proposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(recomb_table[0][0] - center[0],2.0))*sigma_transl*sqrt(2 * pi)) *
	  exp(-pow(intermediate[0][0] - recomb_table[0][0],2) / (2.0 * pow((recomb_table[0][0] - center[0])*sigma_transl,2.0)));

	proposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(recomb_table[0][1] - center[1],2.0)) * sigma_transl * sqrt(2 * pi)) *
	    exp(-pow(intermediate[0][1] - recomb_table[0][1],2) / (2.0 * pow((recomb_table[0][1] - center[1])*sigma_transl,2.0)));
#else
	proposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(recomb_table[0][0] - center[0],2.0))*sigma_transl*euclid*sqrt(2 * pi)) *
	  exp(-pow(intermediate[0][0] - recomb_table[0][0],2) / (2.0 * pow((recomb_table[0][0] - center[0])*
									   sigma_transl*euclid,2.0)));
	proposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(recomb_table[0][1] - center[1],2.0)) * sigma_transl * euclid * sqrt(2 * pi)) *
	    exp(-pow(intermediate[0][1] - recomb_table[0][1],2) / (2.0 * pow((recomb_table[0][1] - center[1])*sigma_transl*euclid,2.0)));
#endif
	//} 
	  
	//cout << "m=" << m << "euclid " << euclid << "(" << intermediate[0][0] << "," << intermediate[0][1] << ")(" << recomb_table[0][0]
	//cout << "first prop (" << i  << "," << proposalDistribution[i] << ", " << 1.0/(sqrt(pow(recomb_table[0][1] - center[1],2.0)) * sigma_transl * euclid * sqrt(2 * pi)) *
	//exp(-pow(intermediate[0][1] - recomb_table[0][1],2) / (2.0 * pow((recomb_table[0][1] - center[1])*sigma_transl*euclid,2.0))) << "," << euclid << ","<< recomb_table[0][1] - center[1] << "\n";
	// in the direction perp x - g in the dimension 
	//for(int j = 1; j < sizeGENOME; j++){
	//if(k == 1){
#ifdef TWO_PARENTS
	//	if(sizeCoupling == 2){
#ifndef EUCLIDIAN
	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(recomb_table[0][1] - center[1],2.0)) * sigma_transl_oriz * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][0] - recomb_table[0][0],2) / (2.0 * pow((recomb_table[0][1] - center[1])*sigma_transl_oriz,2.0)));
	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(recomb_table[0][0] - center[0],2.0)) * sigma_transl_oriz * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][1] - recomb_table[0][1],2) / (2.0 * pow((recomb_table[0][0] - center[0])*sigma_transl_oriz,2.0)));
#else
	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(recomb_table[0][1] - center[1],2.0)) * sigma_transl_oriz * euclid * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][0] - recomb_table[0][0],2) / (2.0 * pow((recomb_table[0][1] - center[1])*sigma_transl_oriz*euclid,2.0)));
	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(recomb_table[0][0] - center[0],2.0)) * sigma_transl_oriz * euclid * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][1] - recomb_table[0][1],2) / (2.0 * pow((recomb_table[0][0] - center[0])*sigma_transl_oriz*euclid,2.0)));
#endif
#else
	  //} else if(sizeCoupling == 3){
	  double intermP;
	  intermP = m/(pow(m,2.0)+1) *(recomb_table[1][1] + 1.0/m * recomb_table[1][0] - recomb_table[0][1] + m * recomb_table[0][0]);
#ifdef EUCLIDIAN
	  //intermP[1] = recomb_table[1][1] - 1.0/m * (intermP[0] - recomb_table[1][0]);
	  double d1 = sqrt(1 + 1.0/pow(m,2)) * sqrt(pow(recomb_table[1][0] - intermP,2.0));

	  intermP = m/(pow(m,2.0)+1) *(recomb_table[2][1] + 1.0/m * recomb_table[2][0] - recomb_table[0][1] + m * recomb_table[0][0]);
	  double d2 = sqrt(1 + 1.0/pow(m,2)) * sqrt(pow(recomb_table[2][0] - intermP,2.0));

	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(recomb_table[0][1] - center[1],2.0)) * sigma_transl_oriz * (d1+d2)/2.0 * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][0] - recomb_table[0][0],2) / (2.0 * pow((recomb_table[0][1] - center[1])*sigma_transl_oriz* (d1+d2)/2.0,2.0)));

	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(recomb_table[0][0] - center[0],2.0)) * sigma_transl_oriz * (d1+d2)/2.0 * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][1] - recomb_table[0][1],2) / (2.0 * pow((recomb_table[0][0] - center[0])*sigma_transl_oriz* (d1+d2)/2.0,2.0)));
#else
	  double intermP1 = m/(pow(m,2.0)+1) *(recomb_table[2][1] + 1.0/m * recomb_table[2][0] - recomb_table[0][1] + m * recomb_table[0][0]);

	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(0.5 *((recomb_table[1][0] - intermP) + (recomb_table[2][0] - intermP1)),2.0)) * 
						     sigma_transl_oriz * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][0] - recomb_table[0][0],2) / (2.0 * pow(0.5*((recomb_table[1][0] - intermP) + (recomb_table[2][0] - intermP1))*sigma_transl_oriz,2.0)));
	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(-0.5 *((recomb_table[1][0] - intermP) + (recomb_table[2][0] - intermP1))/m,2.0))
						     * sigma_transl_oriz * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][1] - recomb_table[0][1],2)/(2.0 * pow(-0.5 *((recomb_table[1][0] - intermP) + (recomb_table[2][0] - intermP1))/m*sigma_transl_oriz,2.0)));
	  
#endif //euclidian
#endif 
	  //} else {
	  //cout << " size sizecoupling is not implemented for wPCTX\n";
	  //exit(1);
	  //}
	//}
	  
#endif //translation complex

	//cout << " (" << i << " , (" << recomb_table[0][0] << " , " << recomb_table[0][1]<< ") , (" << 
	//   prop_recomb_table[0][0] << " , " << prop_recomb_table[0][1] << "), (" << center[0] << " , "
	//    << center[1] << "), (" << euclid << ") , "<< proposalDistribution[table[i][0]] << "," 
	//  << proposedProposalDistribution[table[i][0]] << ")\t";

     
	//1.0 comute the centers of the originals
	for(unsigned long j = 0; j < sizeGENOME; j++){
	  center[j] = 0;
	  double jos = 0;
	  for(unsigned long k = 0; k < sizeCoupling; k++){
	    center[j] += prop_recomb_table[k][j] * prop_fit[k];
	    jos += prop_fit[k];
	  }
	  center[j] /= jos;
	}
	  
	m = (prop_recomb_table[0][1] - center[1])/(prop_recomb_table[0][0] - center[0]);
	//for(unsigned long j = 0; j < sizeCoupling; j++){
#ifdef EUCLIDIAN
	euclid = sqrt(pow(prop_recomb_table[0][1] - center[1],2.0)+ pow(prop_recomb_table[0][0] - center[0],2.0));
#endif
 
#ifndef TRANSLATION_COMPLEX
#ifndef EUCLIDIAN
	proposedProposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(prop_recomb_table[0][0] - center[0],2.0))*sigma_transl*sqrt(2 * pi)) *
	  exp(-pow(recomb_table[0][0] - prop_recomb_table[0][0],2) / (2.0 * pow((prop_recomb_table[0][0] - center[0])*sigma_transl,2.0)));
	proposedProposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(prop_recomb_table[0][1] - center[1],2.0)) * sigma_transl * sqrt(2 * pi)) *
	    exp(-pow(recomb_table[0][1] - prop_recomb_table[0][1],2) / (2.0 * pow((prop_recomb_table[0][1] - center[1])*sigma_transl,2.0)));
#else
	proposedProposalDistribution[table[i][0]] *= exp(-pow(recomb_table[0][0] - prop_recomb_table[0][0],2) / (2.0 * pow((prop_recomb_table[0][0] - center[0])*sigma_transl*euclid,2.0)))/(sqrt(pow(prop_recomb_table[0][0] - center[0],2.0))*sigma_transl*euclid*sqrt(2 * pi));
	proposedProposalDistribution[table[i][0]] *= exp(-pow(recomb_table[0][1] - prop_recomb_table[0][1],2) / (2.0 * pow((prop_recomb_table[0][1] - center[1])*sigma_transl*euclid,2.0)))/(sqrt(pow(prop_recomb_table[0][1] - center[1],2.0)) * sigma_transl * euclid * sqrt(2 * pi));
#endif

#else
	//m = (prop_recomb_table[0][1] - center[1])/(prop_recomb_table[0][0] - center[0]);
	  // in the dirrection x - g
	//	  if(k == 1) {
	intermediate[0][0] = m / (pow(m,2) + 1) * (recomb_table[0][1] + 
		       1.0 /m * recomb_table[0][0] + m * prop_recomb_table[0][0] - prop_recomb_table[0][1]);

	intermediate[0][1] = prop_recomb_table[0][1] + m * (intermediate[0][0] - prop_recomb_table[0][0]); 

	intermediate[1][0] = m / (pow(m,2) + 1) *(prop_recomb_table[0][1] + 
		       1.0 / m * prop_recomb_table[0][0]+ m * recomb_table[0][0] - recomb_table[0][1] );

	intermediate[1][1] = recomb_table[0][1] + m * (intermediate[1][0] - recomb_table[0][0]); 

#ifndef EUCLIDIAN
	proposedProposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(prop_recomb_table[0][0] - center[0],2.0))*sigma_transl*sqrt(2 * pi)) *
	  exp(-pow(intermediate[0][0] - prop_recomb_table[0][0],2) / (2.0 * pow((prop_recomb_table[0][0] - center[0])*sigma_transl,2.0)));
#else
	proposedProposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(prop_recomb_table[0][0] - center[0],2.0))*sigma_transl*euclid*sqrt(2 * pi)) *
	  exp(-pow(intermediate[0][0] - prop_recomb_table[0][0],2) / (2.0 * pow((prop_recomb_table[0][0] - center[0])*sigma_transl*euclid,2.0)));
#endif
	//} 
	  
#ifndef EUCLIDIAN
	proposedProposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(prop_recomb_table[0][1] - center[1],2.0)) * sigma_transl * sqrt(2 * pi)) *
	    exp(-pow(intermediate[0][1] - prop_recomb_table[0][1],2) / (2.0 * pow((prop_recomb_table[0][1] - center[1])*sigma_transl,2.0)));
#else
	proposedProposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(prop_recomb_table[0][1] - center[1],2.0)) * sigma_transl * euclid * sqrt(2 * pi)) *
	    exp(-pow(intermediate[0][1] - prop_recomb_table[0][1],2) / (2.0 * pow((prop_recomb_table[0][1] - center[1])*sigma_transl*euclid,2.0)));
#endif

	// in the direction perp x - g in the dimension 
	//for(int j = 1; j < sizeGENOME; j++){
	//if(k == 1){
#ifdef TWO_PARENTS 
	//	if(sizeCoupling == 2){
#ifndef EUCLIDIAN
	  proposedProposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(prop_recomb_table[0][1] - center[1],2.0)) * sigma_transl_oriz * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][0] - prop_recomb_table[0][0],2) / (2.0 * pow((prop_recomb_table[0][1] - center[1])*sigma_transl_oriz,2.0)));
	  proposedProposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(prop_recomb_table[0][0] - center[0],2.0)) * sigma_transl_oriz * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][1] - prop_recomb_table[0][1],2) / (2.0 * pow((prop_recomb_table[0][0] - center[0])*sigma_transl_oriz,2.0)));
#else
	  proposedProposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(prop_recomb_table[0][1] - center[1],2.0)) * sigma_transl_oriz * euclid * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][0] - prop_recomb_table[0][0],2) / (2.0 * pow((prop_recomb_table[0][1] - center[1])*sigma_transl_oriz*euclid,2.0)));
	  proposedProposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(prop_recomb_table[0][0] - center[0],2.0)) * sigma_transl_oriz * euclid * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][1] - prop_recomb_table[0][1],2) / (2.0 * pow((prop_recomb_table[0][0] - center[0])*sigma_transl_oriz*euclid,2.0)));
#endif
#else
	  //} else if(sizeCoupling == 3){
	  //double intermP;
	  intermP = m/(pow(m,2.0)+1) *(prop_recomb_table[1][1] + 1.0/m * prop_recomb_table[1][0] - prop_recomb_table[0][1] + m * prop_recomb_table[0][0]);
#ifdef EUCLIDIAN
	  //intermP[1] = recomb_table[1][1] - 1.0/m * (intermP[0] - recomb_table[1][0]);
	  d1 = sqrt(1 + 1.0/pow(m,2)) * sqrt(pow(prop_recomb_table[1][0] - intermP,2.0));

	  intermP = m/(pow(m,2.0)+1) *(prop_recomb_table[2][1] + 1.0/m * prop_recomb_table[2][0] - prop_recomb_table[0][1] + m * prop_recomb_table[0][0]);
	  d2 = sqrt(1 + 1.0/pow(m,2)) * sqrt(pow(prop_recomb_table[2][0] - intermP,2.0));

	  proposedProposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(prop_recomb_table[0][1] - center[1],2.0)) * sigma_transl_oriz * (d1+d2)/2.0 * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][0] - prop_recomb_table[0][0],2) / (2.0 * pow((prop_recomb_table[0][1] - center[1])*sigma_transl_oriz* (d1+d2)/2.0,2.0)));

	  proposedProposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(prop_recomb_table[0][0] - center[0],2.0)) * sigma_transl_oriz * (d1+d2)/2.0 * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][1] - prop_recomb_table[0][1],2) / (2.0 * pow((prop_recomb_table[0][0] - center[0])*sigma_transl_oriz* (d1+d2)/2.0,2.0)));
#else
	  intermP1 = m/(pow(m,2.0)+1) *(prop_recomb_table[2][1] + 1.0/m * prop_recomb_table[2][0] - prop_recomb_table[0][1] + m * prop_recomb_table[0][0]);

	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(0.5 * ((prop_recomb_table[1][0] - intermP) + (prop_recomb_table[2][0] - intermP1)),2.0)) 
						     * sigma_transl_oriz * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][0] - prop_recomb_table[0][0],2) / (2.0 * pow(0.5 * ((recomb_table[1][0] - intermP) + (recomb_table[2][0] - intermP1))
										  *sigma_transl_oriz,2.0)));
	  proposalDistribution[table[i][0]] *=  1.0/(sqrt(pow(-0.5 * ((prop_recomb_table[1][0] - intermP) + (prop_recomb_table[2][0] - intermP1))/m,2.0))
						     * sigma_transl_oriz * sqrt(2 * pi)) *
	    exp(-pow(intermediate[1][1] - prop_recomb_table[0][1],2) / (2.0 * pow(-0.5 *((prop_recomb_table[1][0] - intermP) + (prop_recomb_table[2][0] - intermP1))/m
										  *sigma_transl_oriz,2.0)));
#endif
#endif  
	  //} else {
	  //cout << " size sizecoupling is not implemented for wPCTX\n";
	  //exit(1);
	  //}
	//}
	  
#endif //TRANSLATION_COMPLEX
	
//	cout << " (" << i << " , (" << recomb_table[0][0] << " , " << recomb_table[0][1]<< ") , (" << 
	//   prop_recomb_table[0][0] << " , " << prop_recomb_table[0][1] << "), (" << center[0] << " , "
	//   << center[1] << "), (" << euclid << ") , "<< proposalDistribution[table[i][0]] << "," 
	//<< proposedProposalDistribution[table[i][0]] << ")\t";
	}
	else{
	//1.0 comute the centers of the originals
	for(unsigned long j = 0; j < sizeGENOME; j++){
	  center[j] = 0;
	  double jos = 0;
	  for(unsigned long k = 0; k < sizeCoupling; k++){
	    center[j] += recomb_table[k][j] * fit[k];
	    jos += fit[k];
	  }
	  center[j] /= jos;
	}

	  for(int j = 0; j < sizeGENOME; j++)
	    proposalDistribution[table[i][0]] *= exp(-pow(prop_recomb_table[i][j] - recomb_table[i][j],2) / (2.0 * pow((recomb_table[i][j] - center[j])*sigma_transl,2.0)))/(sqrt(pow(recomb_table[i][j] - center[j],2.0))*sigma_transl*sqrt(2 * pi)) ;

	//1.0 comute the centers of the originals
	for(unsigned long j = 0; j < sizeGENOME; j++){
	  center[j] = 0;
	  double jos = 0;
	  for(unsigned long k = 0; k < sizeCoupling; k++){
	    center[j] += prop_recomb_table[k][j] * prop_fit[k];
	    jos += prop_fit[k];
	  }
	  center[j] /= jos;
	}

	for(int j = 0; j < sizeGENOME; j++)
	proposedProposalDistribution[table[i][0]] *= 1.0/(sqrt(pow(prop_recomb_table[i][j] - center[j],2.0))*sigma_transl*sqrt(2 * pi)) *
	  exp(-pow(recomb_table[i][j] - prop_recomb_table[i][j],2) / (2.0 * pow((prop_recomb_table[i][j] - center[j])*sigma_transl,2.0)));

	}
	  //}
	  
	/*   } else { 
//rotation
	  // for(int k = 1; k < sizeGENOME; k++){
	  //  double m = acos(1 - (pow(recomb_table[0][k] - prop_recomb_table[0][k],2.0) + pow(recomb_table[0][0] - 
	  //     prop_recomb_table[0][0],2.0))/(2*pow(recomb_table[0][k] - center[k],2.0) + 
	  //				    2*pow(recomb_table[0][0] - center[0],2.0)));
	  proposalDistribution[i] *= 1;
	  proposedProposalDistribution[i] *= 1;
	      //}
	      }*/
  }
  //cout << "\n";
  delete[] recomb_table;
  delete[] intermediate;
  delete[] center;
  delete[] prop_recomb_table;  
  delete[] fit;
  delete[] prop_fit;
}
#else 
#ifdef SIMPLEX
void binaryPopulation::proposalDistributionCalculation(unsigned long** table, unsigned long dimension){

 unsigned long tempSize = dimension;
  //int temp[sizeGENOME];

  for(unsigned long i = 0; i < sizeMCMCstring; i++){
      proposalDistribution[i] = 1;
      proposedProposalDistribution[i] = 1;
  }

  double rest = 1 - MCMCstring[0]->proportions[sizeCoupling - 2];

  
  /* cout << " Correction \n";

  cout << "probabilities { ";
  for(unsigned long i = 0; i < sizeCoupling; i++)
    cout << " ( " << i << " , " <<  MCMCstring[0]->proportions[i] << " ) ";
  cout << " } \n ";
  */

  for(unsigned long i = 0; i < sizeMCMCstring/sizeCoupling; i++){
    if(MCMCstring[table[i][0]] -> rank != sizeCoupling -1)
      proposalDistribution[table[i][0]] *= pow(0.5,(double)MCMCstring[table[i][0]] -> rank)/sizeGENOME;
    else proposalDistribution[table[i][0]] *= rest;
    
    //cout << "i " << i << " rank " << MCMCstring[table[i][0]] -> rank << " proposal " << proposalDistribution[table[i][0]] << " \n";
 }

 for(unsigned long i = 0; i < sizeMCMCstring/sizeCoupling; i++){
   int small = 0;
   double fitnessChanged = proposedMCMCstring[table[i][MCMCstring[table[i][0]] -> positionChanged]] -> fitness();
   for(unsigned long j = 0; j < sizeCoupling; j++)
     if(j != MCMCstring[table[i][0]] -> positionChanged && proposedMCMCstring[table[i][j]] -> fitness() > fitnessChanged)
       small++;
   if(small != sizeCoupling - 1)
     proposedProposalDistribution[table[i][0]] *= pow(1.0/2.0,(double)small)/sizeGENOME;
   else proposedProposalDistribution[table[i][0]] *= rest;
   
   // cout << "i " << i << " new fitness " << fitnessChanged << " inverse rank " << small << " inverse proposal " << proposedProposalDistribution[table[i][0]] << " \n";
 }
 
}
#endif //simplex 
#endif //wPCTX
#endif //BINARY_RECOMB
#endif //POPULATION_MCMC
#endif //Recombination test
#endif //METROPOLIS_ALGORITHM

/////////////////////
//////////// proposal distribution for NESTED MH algorithm
////////////////////

//reprezinte s_r(x3,x2 \mid x1,x4)/s_r(x1,x2 \mid x3,x4)
#ifdef NESTED_MH
double binaryPopulation::proposalLocus2_2(unsigned long x1,unsigned long x2, unsigned long x3, unsigned long x4){
  //aceasta recombinare este numai pentru x4 from x2 and x1
  if(x1 == x2 && x2 == x3 && x3 == x4)
    return 1.0;
  //1. if x2 and x4 there are the same then no bit were changed 
  if(x2 == x4 && x3 == x1)
    return 1.0;
  //2. if x2 and x4 difera in acelasi fel de parteneri -> atunci propoertia este 1
  if(x3 == x1) 
    if(x2 == x1)
      return (mutation_multi*BLOCKsize);
    else return 1.0/(mutation_multi*BLOCKsize); 
  //3. 
  return 1.0/(mutation_multi*BLOCKsize);

}

double binaryPopulation::proposalDistributionCalculation(double*** table, int dimension){

  unsigned long tempSize = dimension;
  //int temp[sizeGENOME];

  unsigned long*** recomb_table = new unsigned long**[sizeCoupling];
  for(unsigned long i = 0; i < sizeCoupling; i++){
    recomb_table[i] = new unsigned long*[sizeCoupling];
    for(unsigned long j = 0; j < sizeCoupling; j++)
      recomb_table[i][j] = new unsigned long[sizeGENOME*BLOCKsize];
  }
  double proposal = 1;
  
  //#ifndef TWO_CHILDREN 
  
  //1.convert all individuals from reals to integers
  //cout << ") \n Proposed distribution =(";
  convertBinary(table[0][0],recomb_table[0][0],sizeGENOME);
  convertBinary(table[0][1],recomb_table[0][1],sizeGENOME);
  convertBinary(table[1][0],recomb_table[1][0],sizeGENOME);
  convertBinary(table[1][1],recomb_table[1][1],sizeGENOME);
  
  //compute the proposal probabilites for the recombination that is second child
  //we assume that for S_{r,mk}(table[0][0], table[0][1] -> table[1][1])
  for(int k = 0; k < sizeGENOME; k++)
    proposal *= proposalLocus2_2(recomb_table[0][0][k],recomb_table[0][1][k],recomb_table[1][0][k],recomb_table[1][1][k]);

  //cout << ")\n";
  //#endif // TWO_CHILDREN
  
  delete[] recomb_table;
  return proposal;
}
#endif // nested MH

/////////////////
///// initialize the proposed individuals with the current individuals
////////////////
void binaryPopulation::init()
{
  unsigned long tempSize = sizeMCMCstring / sizeCoupling;
#ifndef MIXTURE_REALS
  unsigned long temp[sizeGENOME];  
#else
  double temp[sizeGENOME];  
#endif
 
#ifdef DEBUG_
  cout << "Copy parents to children and initialize them \n";
#endif

  for(unsigned long i = 0; i<sizeMCMCstring ; i++){
	MCMCstring[i] -> getGENOME(temp);
	proposedMCMCstring[i] -> setGENOME(temp);
	proposedMCMCstring[i] -> fitness();
	//proposedMCMCstring[i] -> setTEMPERATURE(MCMCstring[i]->getTEMPERATURE());
	//proposedMCMCstring[i] -> setProposalDistributionVector(proposalDistribution);
	//proposedMCMCstring[i] -> ProposalDistribution();
  }
}

void binaryPopulation::mutation()
{
    // cout << "start mutation \n";
  unsigned long tempSize = sizeMCMCstring / sizeCoupling;
#ifndef REALS
  unsigned long temp[sizeGENOME];  
#else
  double temp[sizeGENOME];
#endif

#ifdef DEBUG_ 
#ifndef METROPOLIS_ALGORITHM 
  cout << "proposalDistribution \n"; 
  for(unsigned long k = 0; k < sizeGENOME; k++)
    cout << "k=" <<k << "=("<<proposalDistribution[k][0] << "," << proposalDistribution[k][1] << ")\t";
  cout <<"\n";
#endif //METROPOLIS_ALGORITHM
  cout << "Mutate children and initialize them \n";
#endif

#ifndef POPULATION_MCMC
  for(unsigned long i = 0; i< sizeMCMCstring ; i++){
    ///???? has to be here?
    //proposedMCMCstring[i] -> getGENOME(temp);
      //cout << "before mutation";
      //for(int k =0; k < sizeGENOME; k++)
      //  cout << temp[k] << "\t";
      //cout << "\n";
      proposedMCMCstring[i] -> MutationGenome(temp);
      //cout << "after mutation";
      //for(int k =0; k < sizeGENOME; k++)
      //	  cout << temp[k] << "\t";
      //    cout << "\n";
      proposedMCMCstring[i] -> setGENOME(temp);
      proposedMCMCstring[i] -> fitness();
      //proposedMCMCstring[i] -> setTEMPERATURE(MCMCstring[i]->getTEMPERATURE());
#ifndef METROPOLIS_ALGORITHM
#ifdef POPULATION_MCMC
      CalculateNewProposalDistribution(i);
#endif
#ifdef INDEPENDENT_SAMPLER
      proposedMCMCstring[i] -> setProposalDistributionVector(PropVector);
      //cout << "i,j" << i<< "," << j << " " <<proposedMCMCstring[table[i][j]] -> ProposalDistribution() << "\n";
#endif //INDEPENDENT_SAMPLER
#endif //n METROPOLIS_ALGORITHM
  }
  
#else //POPULATION_MCMC
  randI = genrand_int32() % sizeMCMCstring;  
  randJ = genrand_int32() % sizeGENOME;
  proposedMCMCstring[randI] -> getGENOME(temp);
  MCMCstring[randI] -> MutationGenome(temp,proposalDistribution,randJ);
  
  proposedMCMCstring[randI] -> setGENOME(temp);
  proposedMCMCstring[randI] -> fitness();
  proposedMCMCstring[randI] -> setTEMPERATURE(MCMCstring[randI]->getTEMPERATURE());
  CalculateNewProposalDistribution(randI);
#ifdef INDEPENDENT_SAMPLER
  proposedMCMCstring[randI] -> setProposalDistributionVector(PropVector);
  //cout << "i,j" << i<< "," << j << " " <<proposedMCMCstring[randI] -> ProposalDistribution() << "\n";
#endif //INDEPENDENT_SAMPLER
#endif // not POPULATION_MCMC
  
#ifdef DEBUG_
  cout << "mutation is finished \n";
#endif
  
}

#ifdef IsRecombination
void binaryPopulation::recombination(unsigned long** table,unsigned long dimension, const unsigned long rECOMB)
{
    
//1. Initilization of the data
  unsigned long tempSize = dimension;
  //int temp[sizeGENOME];
#ifndef MIXTURE_REALS 
  unsigned long** recomb_table = new unsigned long*[2*sizeCoupling];
  for(unsigned long i = 0; i<2*sizeCoupling; i++)
    recomb_table[i] = new unsigned long[sizeGENOME];
#else
  //cout << "recombination is made in the real space \n";
#ifdef BINARY_RECOMB
  double** recomb_reals = new double*[2*sizeCoupling];
  unsigned long** recomb_table = new unsigned long*[2*sizeCoupling];
  for(unsigned long i = 0; i<2*sizeCoupling; i++){
    recomb_reals[i] = new double[sizeGENOME];
    recomb_table[i] = new unsigned long[sizeGENOME*BLOCKsize];
    for(unsigned long j =0; j < sizeGENOME*BLOCKsize; j++)
	recomb_table[i][j] = 0;
  }
#else 
  //cout << " real recombination initialization data\n";
  double** recomb_table = new double*[2*sizeCoupling];
  for(unsigned long i = 0; i<2*sizeCoupling; i++)
    recomb_table[i] = new double[sizeGENOME];
#endif // binary recomb
 
#endif

#ifdef DEBUG_
  cout << "begin recombination \n"; 
#endif

  for(unsigned long i = 0; i<tempSize ; i++){

  //2. prepare the recombination table
#ifdef DEBUG_
      cout << "with intermidare genomes for i=" <<i<<"\n";      
      for(unsigned long j = 0; j<sizeCoupling; j++)
	proposedMCMCstring[table[i][j]] -> print();
#endif 

#ifndef MIXTURE_REALS
      for(unsigned long j = 0; j<sizeCoupling; j++)
	proposedMCMCstring[table[i][j]] -> getGENOME(recomb_table[j]);
#else
      //cout << "genome to recombine, pair " << i << "\n";
      for(unsigned long j = 0; j<sizeCoupling; j++){
#ifdef BINARY_RECOM
	  proposedMCMCstring[table[i][j]] -> getGENOME(recomb_reals[j]);
          //convert everything to integers
	  //cout << " real genome "<< j << "\t"; 
	  //for(int k = 0; k < sizeGENOME; k++)
	  //    cout << recomb_reals[j][k] << "\t";
	  //cout << "\n";
	  convertBinary(recomb_reals[j],recomb_table[j],sizeGENOME);
	  //cout << " genome "<< j << "\t"; 
	  //for(int k = 0; k < sizeGENOME*BLOCKsize; k++)
	  //    cout << recomb_table[j][k];
	  ///cout << "\n";
#else
	  //cout << " real recombination initialization proposed\n";	  
	  proposedMCMCstring[table[i][j]] -> getGENOME(recomb_table[j]);  
#endif // binary_recomb
      }
#endif

   //3. Recombination 
//#ifdef IsRecombination
#ifdef METROPOLIS_ALGORITHM
     recomb_table = MCMCstring[table[i][0]] -> RecombinationGenome(recomb_table,rECOMB);
     //cout << " accepted imediat (" << i << "," << table[i][0] << "," << MCMCstring[table[i][0]] -> accepted << ")\n";
     /*cout << "old genomes ";
     for(int j = 0; j < sizeCoupling; j++){
       cout << "("<< j << ","; 
       for(int k = 0; k < sizeGENOME; k++)
	 cout << "," << recomb_table[j][k];
       cout << ") \t";
     }
     cout << "\n";
     cout << " new genomes ";
     for(int j = 0; j < sizeCoupling; j++){
       cout << "("<< j << ","; 
       for(int k = 0; k < sizeGENOME; k++)
	 cout << "," << recomb_table[j + sizeCoupling][k];
       cout << ") \t";
     }
     cout << "\n";*/
#else
#ifndef POPULATION_MCMC
     MCMCstring[table[i][0]] -> RecombinationGenome(recomb_table,rECOMB);
     //cout << "parent replace" << MCMCstring[table[i][0]] ->parentReplaced << "\n";
#else
     MCMCstring[table[i][0]] -> RecombinationGenome(recomb_table,rECOMB,proposalDistribution);
#endif 
#endif //Metropolis
//#endif //IsRecombination

     //4. compute the proposal distribution
#ifndef METROPOLIS_ALGORITHM 
#ifdef POPULATION_MCMC
     CalculateNewProposalDistributionRecombination(table[i][0],table[i][1]);
#endif
#endif //METROPOLIS_ALGORITHM

     //5. set the new individuals as the proposal individuals
#ifndef MIXTURE_REALS
     for(unsigned long j = 0; j<sizeCoupling; j++){
	 proposedMCMCstring[table[i][j]] -> setGENOME(recomb_table[j+sizeCoupling]);
	 proposedMCMCstring[table[i][j]] -> fitness();
	 //	 proposedMCMCstring[table[i][j]] -> setTEMPERATURE(MCMCstring[table[i][j]]->getTEMPERATURE());
#ifndef METROPOLIS_ALGORITHM 
#ifdef POPULATION_MCMC
#ifdef INDEPENDENT_SAMPLER
	 proposedMCMCstring[table[i][j]] -> setProposalDistributionVector(PropVector);
#endif
#endif //METROPOLIS_ALGORITHM 
#endif //POPULATION_MCMC
       }
#else //mixture reals
      //cout << "after recombination, pair " << i << "\n";
      for(unsigned long j = 0; j<sizeCoupling; j++){
	//cout << " genome "<< j << "\t"; 
	//    for(int k = 0; k < sizeGENOME*BLOCKsize; k++)
	//cout << recomb_table[j+sizeCoupling][k];
	       //   cout << "\n";
#ifdef BINARY_RECOMB	     
	     convertReal(recomb_reals[j+sizeCoupling],recomb_table[j+sizeCoupling],sizeGENOME);
	     //cout << " real genome "<< j << "\t"; 
	     //for(int k = 0; k < sizeGENOME; k++)
	     //  cout << recomb_reals[j+sizeCoupling][k] << "\t";
	     //cout << "\n";

	     proposedMCMCstring[table[i][j]] -> setGENOME(recomb_reals[j+sizeCoupling]);	 
#else
	     /*cout << " real genome "<< j << "\t"; 
	     for(int k = 0; k < sizeGENOME; k++)
	       cout << recomb_table[j+sizeCoupling][k] << "\t";
	       cout << "\n";*/
	     proposedMCMCstring[table[i][j]] -> setGENOME(recomb_table[j+sizeCoupling]);
#endif //binary recomb	     
	     proposedMCMCstring[table[i][j]] -> fitness();
	     //	 proposedMCMCstring[table[i][j]] -> setTEMPERATURE(MCMCstring[table[i][j]]->getTEMPERATURE());
#ifndef METROPOLIS_ALGORITHM 
#ifdef POPULATION_MCMC
#ifdef INDEPENDENT_SAMPLER
	 proposedMCMCstring[table[i][j]] -> setProposalDistributionVector(PropVector);
#endif
#endif //METROPOLIS_ALGORITHM 
#endif //POPULATION_MCMC
       }
#endif //mixture reals
	

#ifdef DEBUG_
     cout << "with final genomes \n";      
     for(unsigned long j = 0; j<sizeCoupling; j++)
       proposedMCMCstring[table[i][j]] -> print();
#endif 

    }

#ifdef DEBUG_
  cout << "recombination is finished \n";
#endif

  delete[] recomb_table;
}

//#ifndef MIXTURE_REALS
void binaryPopulation::mutation(unsigned long** table,unsigned long dimension)
{
  unsigned long tempSize = dimension;
#ifndef MIXTURE_REALS
  unsigned long temp[sizeGENOME];
#else
  double temp[sizeGENOME];
  //  unsigned long tempInt[sizeGENOME*BLOCKsize];
#endif  

  //#ifdef DEBUG_
  /*if(generation % SEE_SAMPLE == 0){
    cout << "begin mutation " << MCMCstring[0]->mutation << "\n"; 
    }*/
 //#endif
  
  for(unsigned long i = 0; i<tempSize ; i++){
      
      proposedMCMCstring[table[i][0]] -> getGENOME(temp);

#ifdef MIXTURE_REALS
      proposedMCMCstring[table[i][0]] -> MutationGenome(temp);
#else // mixture reals
      proposedMCMCstring[table[i][0]] -> MutationGenome(temp);
#endif // mixture reals
      proposedMCMCstring[table[i][0]] -> setGENOME(temp);
      proposedMCMCstring[table[i][0]] -> fitness();
      
      
      /* if(generation % SEE_SAMPLE == 0)
    for(unsigned long j = 0; j<sizeCoupling; j++){
      for(int k = 0; k < sizeGENOME; k++)
	cout << proposedMCMCstring[table[i][j]] -> genome[k];
      //cout << "\n";
      }*/
      //#endif 
      
  }
  
  //#ifdef DEBUG_
  //if(generation % SEE_SAMPLE == 0)  cout << "mutation is finished \n";
  //#endif
  
}
//#endif //MIXTURE_REALS
#endif

unsigned long binaryPopulation::postProcessing(){
    unsigned long stop = 0;
    
    for(unsigned long i = 0; i<sizeMCMCstring; i++) {	
#ifdef RECOMB_MUT
	if(i == 0){
#endif
#ifdef SAMPLE
//    if(generation >= burn_in && (generation - burn_in)%interSamplesize == 0){
	//cout << "verifica" << i <<"   " << generation << " diversity" << diversityMeasure << "\n";
	MCMCstring[i] -> addElement();
	MCMCstring[i]->generation_Update(); //?
//    }
	diversityFileBookeeping(i);
#ifndef MIXTURE_REALS
	diversityBookeeping(i);
#endif 
	
#else	//SAMPLE
	if((MCMCstring[i]) -> stop() == 1) {
	    stop++;
	    if(stop %1000==999) {
		cout << "Maximum in here \n"; 
	    }
	}
#endif //SAMPLE
#ifdef RECOMB_MUT
	}
#endif
    }
    
#ifndef SAMPLE
    if(stop != 0) print();
    
#ifdef RUN_SEQUENCES
    if(generation % 1000 == 999){
	print();
	char t;
	cin >> t;
	if(t == 's') stop =sizeMCMCstring;
    }
#endif //RUN_SEQUENCES
#endif //SAMPLE
    return 0;
}

#ifndef MIXTURE_REALS
void binaryPopulation::diversityBookeeping(unsigned long i){
    //diversity Bookeeping
    double tempValue;
#ifdef THRESHOLD
#ifndef THRESHOLD_DISTRIBUTION
  tempValue = f(MCMCstring[i]->getValue());
  if(tempValue <= 0) return;
#else
  tempValue = MCMCstring[i]->getValue();
#endif
#else 
  tempValue = MCMCstring[i]->getValue();
#endif
#ifdef FLOOR
  tempValue =  floor(tempValue);
#endif
    
    //diversity bookeeping
    unsigned long count = 0;
    for(unsigned long t=0;t< diversityMeasure;t++) 
	if(tempValue == transform_table[t]){
	    diversity[t]++;
	    count = 1;
	}
    
    if(count == 0) {
	cout << "Eroare in diversity in fisier binaryPopulation_trap.cpp:1494\n";
	cout << "error in chain " << i << "which reports a value " << tempValue << 
	  " at generation " << runs <<"\n"; 
      exit(1);
  }
  
#ifdef TRAP_LIKE_FUNCTION
    //detailed diversity bookeeping
    unsigned long temp[sizeGENOME];
    MCMCstring[i]->getGENOME(temp);

#ifdef DEBUG_
    for(unsigned long j= 0; j < sizeGENOME; j++)
	cout << temp[j];
    cout <<"\n";
#endif

    for(unsigned long j= 0; j < sizeGENOME/BLOCKsize; j++){
      unsigned long nr_temp = 0; 
      for(unsigned long k = 0; k < BLOCKsize; k++)
	nr_temp += temp[j*BLOCKsize+k];
      //cout << "("<<j<<","<<nr_temp<<")\t";
      diversityDetail[j][nr_temp]++;
    }
#endif
}
#endif //mixture reals

void binaryPopulation::diversityFileBookeeping(unsigned long i){
    double tempValue;
#ifndef MIXTURE_REALS
#ifdef THRESHOLD
#ifndef THRESHOLD_DISTRIBUTION
  tempValue = f(MCMCstring[i]->getValue());
  if(tempValue <= 0) return;
#else
  tempValue = MCMCstring[i]->getValue();
#endif
#else 
  tempValue = MCMCstring[i]->getValue();
#endif
#ifdef FLOOR
  tempValue =  floor(tempValue);
#endif
#else //mixture reals
///thresold
#ifndef THRESHOLD
  tempValue = MCMCstring[i]->getValue();
#else
  tempValue = f(MCMCstring[i]->fitness());
  //tempValue = f(MCMCstring[i]->getValue());
  if(tempValue <= 0) return;  
#endif
#endif // MIXTURE_REALS

#ifdef MULTI_PICK_TRAP_FUNCTION
   if(MCMCstring[i] -> stop() == 0) return; 
#endif

  unsigned long count = 0;
#ifndef MIXTURE_REALS
  unsigned long trans_value = 0;
  for(unsigned long t=0;t< diversityMeasure;t++) 
      if(tempValue == transform_table[t]){
	  trans_value = t;
	  count = 1;
      }
    
  if(count == 0) {
      cout << "Eroare in diversity in fisier binaryPopulation_trap.cpp:1564\n";
      cout << "error in chain " << i << "which reports a value " << 
	  tempValue << " at generation " << runs << " when diversity measure is" 
	   << diversityMeasure << "\n"; 
      exit(0);
  }
#else
  double trans_value = 0;
  //the same value, not transformed
  trans_value = tempValue;
  count = 1;
#endif //MIXTURE_REALS

#ifndef MIXTURE_REALS
  unsigned long temp[sizeGENOME];
  MCMCstring[i]->getGENOME(temp);
#else
  double recomb_reals[sizeGENOME];

  MCMCstring[i] -> getGENOME(recomb_reals);
  //cout << "\n from real ";
  //for(int j = 0; j < sizeGENOME; j++)
  //    cout << recomb_reals[j] << "\t";
  //cout << "\n";
  unsigned long temp[sizeGENOME*BLOCKsize];
  convertBinary(recomb_reals,temp,sizeGENOME);

  //cout << "\n to binary ";
  //for(int j = 0; j < sizeGENOME*BLOCKsize; j++)
  //    cout << (unsigned long)temp[j];
  //cout << "\n";
  //cout << "\n from real ";
  //for(int j = 0; j < sizeGENOME; j++)
  //    cout << recomb_reals[j] << "\t";
  //cout << "\n";

#endif 

#ifdef FILE_BOOKEEPING
  // diversity bookeeping
  if (diversityMeasure == 0) return;
  diversityFile.seekg(0,ios::beg);
  unsigned long lenghtFI = diversityFile.tellg();
  diversityFile.seekg(0,ios::end);
  unsigned long lenghtFE = diversityFile.tellg();
  
  char buffer[100];
  char* tempGenome = new char[sizeGENOME+2];
  for(unsigned long j = 0; j < sizeGENOME; j++)
      if(temp[j] == 0) tempGenome[j] = '0';
      else if(temp[j] == 1) tempGenome[j] = '1';
  tempGenome[sizeGENOME] = '\0';
  
  unsigned long count = 0;
  
  diversityFile.seekg(0,ios::beg);
  //cout << "To be compared "<< tempGenome << "\n";
  while((lenghtFE-lenghtFI) > sizeGENOME && count == 0){
      diversityFile.getline(buffer,100);  
      //cout << "is compared with "<< buffer << "\n";
      if(!strcmp(buffer,tempGenome)) count ++;
      diversityFile.seekg(0,ios::cur);
      lenghtFI = diversityFile.tellg();
  }
  
  if(count == 0) {
      diversityFile.flush();
      diversityFile.close();
      diversityFile.open(diversityTempFILE,ios::in|ios::out);
      diversityFile.seekp(0,ios::end);
      //cout << "Position " << diversityFile.tellp();
      tempGenome[sizeGENOME+1] = '\n';
      if(!diversityFile.write(tempGenome,sizeGENOME+2))
	  cout << "error of writing\n";
      //cout << "tempGenome is writed "<< tempGenome <<"\n";
      //diversityFile.flush();
      
      procentRightSolutions++;
      
      diversityRightSolutions[trans_value] ++;
      
      weightRightSolutions += tempValue;
  }     
  myHashtables -> store_individual(temp,sizeGENOME);
  
  delete[] tempGenome;
#else

#ifndef MIXTURE_REALS
  unsigned long size = sizeGENOME;
#else
  unsigned long size = sizeGENOME * BLOCKsize;
#endif

  //the individual introduced
  //cout << "we have for hashtable ";
  //for(int i = 0; i < sizeGENOME*BLOCKsize; i++)
  //    cout << temp[i];
  //cout << "\t";

#ifdef HISTOGRAM
  if(myHashtables -> check_individual((unsigned long*)temp,size) == 0){
      procentRightSolutions++;

      //  cout << "found\n";

#ifndef MIXTURE_REALS
      diversityRightSolutions[trans_value]++;
#endif
      
      weightRightSolutions += tempValue;

#ifndef MIXTURE_REALS
#ifdef TWO_ATTRACTORS_FUNCTION
      unsigned long onesM = 0;
      for(unsigned long ones = 0; ones < size; ones++){
	  if(temp[ones] != 0 && temp[ones] != 1){
	      cout << "Error in genome binaryPopulation_trap.cpp:1674\n";
	      exit(1);
	  }
	  else onesM += temp[ones];
      }
#ifndef SINGLE_vs_MULTIPLE
      ActualHistogram[runs-1][onesM]++;
#else
      //atentie aceasta conditie este buna numai pentru 2 chains
      /*if(i < sizeMCMCstring/2.0){
       for(unsigned long chains = 0; chains < (int)(log(sizeMCMCstring)/log(2.0))+1;chains++)
	  ActualHistogram[chains][runs-1][onesM]++;
      }
      else ActualHistogram[(int)(log(sizeMCMCstring)/log(2.0))][runs-1][onesM]++; 
      */
      ActualHistogram[(unsigned long)(log(sizeMCMCstring)/log(2))][runs-1][onesM]++;
#endif //ifndef SINgle vs multiple
#endif //ifdef two attractorfunctions
#endif //mixture reals

  } //else cout << "not found \n";
#endif //histogram

  //verificare ce pune in Hashtable
  //cout << "insert hashtable \t";
  //for(unsigned long inter =0; inter < size; inter++)
  //    cout << temp[inter];
  //cout << " = " << tempValue << "\n";
#ifdef HISTOGRAM
  myHashtables -> store_individual(temp,size);
  //verificare cit este reals
  //  if(generation % SEE_SAMPLE == 0){
    //numara in care dintre regiuni sa pun numarul
    // in care din generatii sa pun 
    unsigned long gen = (unsigned long) (sampleSize_Const/SEE_SAMPLE);
    int nrGeneration = (unsigned long) ((sampleSize_Const - generation)/SEE_SAMPLE);
    //cout << "introduce ot generation " << generation << " nr table " << nrGeneration << "\n";
    //place in the histogram
    unsigned long indexIK = 0;
#ifndef BIVARIATE
    for(int ik = 0; ik < sizeGENOME; ik++){
      indexIK += ((unsigned long) recomb_reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
    }
#else //bivariate
    for(int ik = 0; ik < sizeGENOME; ik++){
      indexIK += ((unsigned long) (recomb_reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
    }
#endif //bivariate

    while(nrGeneration >= 0){
      //cout << " (" << nrGeneration << "," << gen - (unsigned long)nrGeneration << ") \t";

      globalDiversity[gen - (unsigned long)nrGeneration][runs-1][indexIK]++;
      globalDiversityRightSolutions[gen - (unsigned long)nrGeneration][runs-1][indexIK] += 1.0/tempValue;

      nrGeneration--;
    }
    //cout << "\n";
  //}

#endif //histogram
  
#endif //ifdef file bookipingf
}

#ifndef TRACE_RUN
#ifndef SINGLE_vs_MULTIPLE
#ifndef MIXTURE_REALS 
double binaryPopulation::printDiversityGNUFILE(){
/*  ifstream mySample(sampleFile);
    if(!mySample.is_open()){
    cout <<"Error in opening the distribution file\n";
    exit(1);
    }
    //unsigned long dataFromSampleFile;
    char temp[100];
*/
    ofstream myData(gnuplotFILE,ios::out | ios:: app);
    if(!myData.is_open()){
	cout <<"Error in opening the diversity file\n";
	exit(1);
    }
    
    double performance[diversityMeasure];
    for(unsigned long i = 0; i < diversityMeasure; i++)
	performance[i] = 0;
    
    double performance_runs[diversityMeasure];
    for(unsigned long i = 0; i < diversityMeasure; i++)
	performance_runs[i] = 0;
    
    /*double sampleData[diversityMeasure];
      for(int i = 0; i < diversityMeasure; i++){
      mySample.getline(temp, 100);
      char* r = temp; while(r[0] != '\t') r++;
      r++;
      sampleData[i] = (unsigned long) atof(r);
      //cout << "vede " << temp << "\t " << sampleData[i] << "\n";
      }*/
    
    for(unsigned long i = 0; i < nr_runs; i++){
	//myData << "Nr runs = " << i << "\n";
	for(unsigned long j= 0; j < diversityMeasure; j++){
	    if (diversityTrueDistribution[j] != 0 & globalDiversityRightSolutions[i][j] != 0) {
		performance[j] += 
		    log((double)globalDiversityRightSolutions[i][j]/(double)diversityTrueDistribution[j]);
		performance_runs[j]++;
/*		myData << j <<"\t"<<
		    (double)globalDiversity[i][j]<<"\t" << 
		    globalDiversityRightSolutions[i][j] << "\t" <<
		    log((double)globalDiversityRightSolutions[i][j]/(double)diversityTrueDistribution[j])
		    <<"\n";*/
	    }
	    else {/*
		if (diversityTrueDistribution[j] == 0) 
		    myData << j <<"\t"<<(double)globalDiversity[i][j]<<"\t" << 0 
			   << "\t" << "\n";
		else if(globalDiversityRightSolutions[i][j] == 0) 
		    myData << j <<"\t"<<(double)globalDiversity[i][j]<<"\t" << 0 
		    << "\t" << -pow(10,30) << "\n";*/
	    }
	} 
    }
    
#ifdef MULTIPLE_RUNS
#ifdef MUTATION_VARIES
    myData << "#Mutation - Diversity mean " << sizeGENOME << "\t" << sampleSize << "\t" << MCMCstring[0]->mutation <<"\n";
#else
    //myData << "#\n Recombination - Diversity mean " << sizeGENOME << "\t" << sampleSize << "\t" 
//	   << MCMCstring[0]->PROB_nPOINT <<"\n";
    myData << "#\n Recombination - Diversity mean " << sizeGENOME << "\t" << sampleSize << "\t" 
	   << MCMCstring[0]->UNIF_RECOMB_PARAM <<"\n";
#endif
#endif
//    for(int i = 0; i < diversityMeasure; i++)
//	myData << i << " " << performance[i]/(double)performance_runs[i] << "\n";
    
//    myData << "Diversity mean histrogram \n";
    double st_dev_perf[diversityMeasure];
    for(unsigned long i = 0; i < diversityMeasure; i++){
	performance[i] = 0;
	for(unsigned long j = 0; j < nr_runs; j++){
	    performance[i] += (double)globalDiversityRightSolutions[j][i]/(double)diversityTrueDistribution[i];
	}
	performance[i] /= (double)nr_runs;
	st_dev_perf[i] = 0;
	for(unsigned long j = 0; j < nr_runs; j++){
	    st_dev_perf[i] += pow(performance[i] - (double)globalDiversityRightSolutions[j][i]/(double)diversityTrueDistribution[i],2.0);
      }
      myData << i << "\t" << performance[i] << "\t" << 
	  sqrt(st_dev_perf[i]/(double)nr_runs) << "\n";
    }
    
/*#ifdef ACCEPTANCE_RATIO
    myData << "\n acceptance ratio \n";
    double accepted = 0;
    for(unsigned long i = 0; i < nr_runs; i++){
	myData << "run = " << i  <<" acceptance= " <<
	    acceptance_ratio[i]/(sampleSize) << "\n";
	accepted += acceptance_ratio[i];
    }
    accepted /= (sampleSize*nr_runs);
    meanAcceptance = accepted;
    myData << "overall acceptance" << accepted << "\n";
    #endif */
    
#ifdef TWO_ATTRACTORS_FUNCTION
    double firstBestMean =0;
    double secondBestMean = 0;
    double minBestMean = 0;
    double meanBestMean = 0;
    //myData << " mean best " << "\n";
    for(unsigned long i = 0; i < nr_runs; i++){
	myData << i << "\t"<< bestFitnessVector[i][0] << "\t"<< bestFitnessVector[i][1] << "\n";
	if(bestFitnessVector[i][0] != 0)
	    firstBestMean += bestFitnessVector[i][0];
	else firstBestMean += sizeMCMCstring * sampleSize;
	if(bestFitnessVector[i][1] != 0)
	    secondBestMean+= bestFitnessVector[i][1];
	else secondBestMean+= sizeMCMCstring * sampleSize;
	meanBestMean += abs(bestFitnessVector[i][0] - bestFitnessVector[i][1]);
	if(bestFitnessVector[i][0]>bestFitnessVector[i][1] && bestFitnessVector[i][1] != 0)
	    minBestMean += bestFitnessVector[i][1];
	else if(bestFitnessVector[i][0] != 0)
	    minBestMean += bestFitnessVector[i][0];
	else minBestMean += sizeMCMCstring * sampleSize;
    }
    myData << " means  " << firstBestMean/nr_runs << "\t" << secondBestMean/nr_runs << 
	"\t" << minBestMean/nr_runs << "\t" << meanBestMean/nr_runs << "\n";
    meanFitness = minBestMean/nr_runs;
    diffFitness = meanBestMean/nr_runs;
#else 
    double firstbestFitness = 0;
    double firstExists = 0, meanFirstExists = 0; 
    //myData << " mean first best " << "\n";
    for(unsigned long i = 0; i < nr_runs; i++)
	if(bestFitnessVector[i] != 0){
	    //    myData << i << "\t"<< bestFitnessVector[i] << "\n";
	    firstbestFitness += bestFitnessVector[i];
	    firstExists ++;
	    meanFirstExists += bestFitnessVector[i];
	} else {
	    //myData << i << "\t"<< sizeMCMCstring * sampleSize << "\n";
	    firstbestFitness += sizeMCMCstring * sampleSize;
	}
    
    meanFitness = firstbestFitness/nr_runs;
    meanFirstExists /= firstExists;

    double std_best = 0;
    double std_exists = 0;
    for(unsigned long i = 0; i < nr_runs; i++)
	if(bestFitnessVector[i] != 0){
	    std_best += pow(bestFitnessVector[i] - meanFitness,2.0);
	    std_exists += pow(bestFitnessVector[i] - meanFirstExists,2.0);
	} else {
	    std_best += pow(sizeMCMCstring * sampleSize - meanFitness,2.0);
	}

    std_best = pow(std_best/(double)nr_runs,0.5);
    std_exists = pow(std_exists/(double)nr_runs,0.5);

    double meansecond = 0;
    //myData << " mean second best " << "\n";
    for(unsigned long i = 0; i < nr_runs; i++)
	if(secondBestFitnessVector[i] != 0){
	    //    myData << i << "\t"<< secondBestFitnessVector[i] << "\n";
	    meansecond += secondBestFitnessVector[i];
	} else {
	    //myData << i << "\t"<< sizeMCMCstring * sampleSize << "\n";
	    meansecond += sizeMCMCstring * sampleSize;
	}
    
    meansecond /= nr_runs;

    double std_second = 0;
    for(unsigned long i = 0; i < nr_runs; i++)
	if(secondBestFitnessVector[i] != 0){
	    std_second += pow(secondBestFitnessVector[i] - meansecond,2.0);
	} else {
	    std_second += pow(sizeMCMCstring * sampleSize - meansecond,2.0);
	}

    std_second += pow(std_second/(double)nr_runs,0.5);

    //myData << "#mean second = " << meansecond << "\n";
    myData << "#" << sizeMCMCstring << "\t" <<  sampleSize << "\t" << meanFitness << "\t" << std_best << "\t" 
	   << meanFirstExists << "\t" << std_exists << "\t" << meansecond << "\t" << std_second << "\n";

#endif
    
    myData.flush();
    myData.close();
    return 1;
}

#else // mixture reals
#ifdef HISTOGRAM
#ifndef MIXTURE_BIVARIATE
double binaryPopulation::printDiversityGNUFILE(){

  cout <<"print diversity \n";
  
  double** performance = new double*[diversityMeasure];
  for(unsigned long j = 0; j < diversityMeasure; j++){
    performance[j] = new double[diversityMeasure];
  }
  
  
  for(unsigned long gen = 0; gen <= (unsigned long)(sampleSize_Const/SEE_SAMPLE); gen++){
    
    char* fileN = new char[50];
    char* fileNI = new char[50];
    char* fileNII = new char[50];
    char* fileNIII = new char[50];
    char* fileNIV = new char[50];
#ifdef TWO_PICKS
    char* fileNV = new char[50];
#endif

    fileN = strcpy (fileN, "gnuFILE_mut\0");
    fileNI = strcpy (fileNI, "gnuFILE_mut_I\0");
    fileNII = strcpy (fileNII, "gnuFILE_mut_II\0");
    fileNIII = strcpy (fileNIII, "gnuFILE_mut_III\0");
    fileNIV = strcpy (fileNIV, "gnuFILE_mut_IV\0");
#ifdef TWO_PICKS
    fileNV = strcpy (fileNV, "gnuFILE_mut_V\0");
#endif 

    if(gen == 0){
      fileN = strcat (fileN, "_0.txt\0");
      fileNI = strcat(fileNI,"_0.txt\0");
      fileNII = strcat(fileNII,"_0.txt\0");
      fileNIII = strcat(fileNIII,"_0.txt\0");
      fileNIV = strcat(fileNIV,"_0.txt\0");
#ifdef TWO_PICKS 
      fileNV = strcat(fileNV,"_0.txt\0");
#endif
    }
    else if(gen == 1){ 
      fileN = strcat (fileN, "_500.txt\0");
      fileNI = strcat(fileNI,"_500.txt\0");
      fileNII = strcat(fileNII,"_500.txt\0");
      fileNIII = strcat(fileNIII,"_500.txt\0");
      fileNIV = strcat(fileNIV,"_500.txt\0");
#ifdef TWO_PICKS 
      fileNV = strcat(fileNV,"_500.txt\0");
#endif
    }
    else if(gen == 2){
      fileN = strcat (fileN, "_1000.txt\0");
      fileNI = strcat(fileNI,"_1000.txt\0");
      fileNII = strcat(fileNII,"_1000.txt\0");
      fileNIII = strcat(fileNIII,"_1000.txt\0");
      fileNIV = strcat(fileNIV,"_1000.txt\0");
#ifdef TWO_PICKS 
      fileNV = strcat(fileNV,"_1000.txt\0");
#endif
    }  
    else if(gen == 3){	
      fileN = strcat (fileN, "_1500.txt\0");
      fileNI = strcat(fileNI,"_1500.txt\0");
      fileNII = strcat(fileNII,"_1500.txt\0");
      fileNIII = strcat(fileNIII,"_1500.txt\0");
      fileNIV = strcat(fileNIV,"_1500.txt\0");
#ifdef TWO_PICKS 
      fileNV = strcat(fileNV,"_1500.txt\0");
#endif
    }
    else if(gen == 4){
      fileN = strcat (fileN, "_2000.txt\0");
      fileNI = strcat(fileNI,"_2000.txt\0");
      fileNII = strcat(fileNII,"_2000.txt\0");
      fileNIII = strcat(fileNIII,"_2000.txt\0");
      fileNIV = strcat(fileNIV,"_2000.txt\0");
#ifdef TWO_PICKS 
      fileNV = strcat(fileNV,"_2000.txt\0");
#endif
    }
    else {
      fileN = strcat (fileN, "_2500.txt\0");
      fileNI = strcat(fileNI,"_2500.txt\0");
      fileNII = strcat(fileNII,"_2500.txt\0");
      fileNIII = strcat(fileNIII,"_2500.txt\0");
      fileNIV = strcat(fileNIV,"_2500.txt\0");
#ifdef TWO_PICKS 
      fileNV = strcat(fileNV,"_2500.txt\0");
#endif
    }

    cout << "file name" << fileN << " because " << gen << "\n";

    ofstream myData(fileN);
    ofstream myDataI(fileNI);
    ofstream myDataII(fileNII);
    ofstream myDataIII(fileNIII);
    ofstream myDataIV(fileNIV);
#ifdef TWO_PICKS 
    ofstream myDataV(fileNV);
#endif

    delete[] fileN;
    delete[] fileNI;
    delete[] fileNII;
    delete[] fileNIII;
    delete[] fileNIV;
#ifdef TWO_PICKS 
    delete[] fileNV;
#endif 

    if(!myData.is_open()){
	cout <<"Error in opening the diversity file\n";
	exit(1);
    }

#ifdef TWO_PICKS
    //cout << "chiar a intrat aici?\n";
    char* best = new char[50];
    best = strcpy (best, "best_19\0");
    ofstream file_best_37(best,ios::out | ios:: app);

    best = strcpy (best, "best_137\0");
    ofstream file_best_112(best,ios::out | ios:: app);


    delete[] best;
#else
    
    ofstream myBest("best_fitness",ios::out | ios:: app);
    if(!myBest.is_open()){
      cout <<"Error in opening the diversity file\n";
      exit(1);
    }
#endif

    //#ifdef MULTIPLE_RUNS
#ifndef IsRecombination
    myData << "#Mutation - Diversity mean " << sizeGENOME << "\t" << sampleSize << "\t" << MCMCstring[0]->mutation <<"\n";
#else
    //myData << "#\n Recombination - Diversity mean " << sizeGENOME << "\t" << sampleSize << "\t" 
//	   << MCMCstring[0]->PROB_nPOINT <<"\n";
    myData << "#\n Recombination - Diversity mean " << sizeGENOME << "\t" << sampleSize << "\t" 
#ifdef BINARY_RECOMB
	   << UNIF_RECOMB_PARAM << "\t"<< " and mutation rate " << MCMCstring[0]->mutation <<"\n";
#else
	   << PARENT_CENTRIC << "\t"<< " and mutation rate " << MCMCstring[0]->mutation <<"\n";
#endif //binary recomb
#endif
    //#endif //multiple runs
    
    //if (diversityTrueDistribution[j][k] != 0 & globalDiversityRightSolutions[i][j][k] != 0) {
    /////bin hight at a point in sapce
    for(unsigned long j = 0; j < diversityMeasure; j++)
      for(unsigned long k = 0; k < diversityMeasure; k++)
	performance[j][k] = 0;
    
    double tempDiv = 0;
    
    //measure the mean of the histogram
	for(unsigned long j= 0; j < diversityMeasure; j++)
	  for(unsigned long k= 0; k < diversityMeasure; k++){
	    for(unsigned long i = 0; i < nr_runs; i++)
	      if(globalDiversity[gen][i][j][k] != 0 && globalDiversityRightSolutions[gen][i][j][k] != 0)
		performance[j][k] += (double)globalDiversity[gen][i][j][k]/(double)globalDiversityRightSolutions[gen][i][j][k];	  
	    performance[j][k] /= (double)nr_runs;
	  }
    
    double st_dev_perf[diversityMeasure][diversityMeasure];
    for(unsigned long k = 0; k < diversityMeasure; k++){
	for(unsigned long j = 0; j < diversityMeasure; j++){
	  myData << k << "\t" << j << "\t";
	    st_dev_perf[k][j] = 0;
	    for(unsigned long i = 0; i < nr_runs; i++){
		/////bin hight at a point in sapce
	      if(globalDiversity[gen][i][k][j] != 0 && globalDiversityRightSolutions[gen][i][k][j] != 0){
	        tempDiv = globalDiversity[gen][i][k][j]/globalDiversityRightSolutions[gen][i][k][j];
		st_dev_perf[k][j] += pow(performance[k][j] - tempDiv,2.0);
	      }
	    }
	    myData << performance[k][j] << "\t" << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";

#ifndef TWO_PICKS
	    if(k == diversityMeasure/2)
	      myDataI << k << "\t" << j << "\t" << performance[k][j] << "\t" 
		      << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";
	    if(j == diversityMeasure/2)
	      myDataII << k << "\t" << j << "\t" << performance[k][j] << "\t" 
		      << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";
	    if(j == k)
	      myDataIII << k << "\t" << j << "\t" << performance[k][j] << "\t" 
		      << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";
	    if(k == diversityMeasure - 1 - j)
	      myDataIV << k << "\t" << j << "\t" << performance[k][j] << "\t" 
		      << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";
#else
	    if(j == k)
	      myDataI << k << "\t" << j << "\t" << performance[k][j] << "\t" 
		      << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";

	    if(k == 37) //diversityMeasure/4)
	      myDataII << k << "\t" << j << "\t" << performance[k][j] << "\t" 
		      << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";

	    if(k == 112) //3* diversityMeasure/4)
	      myDataIII << k << "\t" << j << "\t" << performance[k][j] << "\t" 
		      << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";

	    if(j == 37) //diversityMeasure/4)
	      myDataIV << k << "\t" << j << "\t" << performance[k][j] << "\t" 
		      << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";

	    if(j == 112) //3* diversityMeasure/4)
	      myDataV << k << "\t" << j << "\t" << performance[k][j] << "\t" 
		      << sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";

#endif //two picks

#ifndef TWO_PICKS
	    if(k == diversityMeasure/2 && j == diversityMeasure/2)
	      myBest << gen*SEE_SAMPLE << "\t" << k << "\t" << j << "\t" << performance[k][j] << "\t" << 
		sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j] << "\n";
#else
	    if(k == 37 && j == 37)
	      file_best_37 << gen*SEE_SAMPLE << "\t" << k << "\t" << j << "\t" << performance[k][j] << "\t" << 
		sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j];

	    if(k == 112 && j == 112)
	      file_best_112 << gen*SEE_SAMPLE << "\t" << k << "\t" << j << "\t" << performance[k][j] << "\t" << 
		sqrt(st_dev_perf[k][j]/(double)nr_runs) << "\t" << trueDistribution[k][j];

#endif
	}
	myData << "\n";
	//if(j == diversityMeasure/2) myDataII << "\n";
	//if(k == diversityMeasure/2) myDataI << "\n";
	//if(j == k) myDataIII << "\n";
	//if(k == diversityMeasure - 1 - j) myDataIV << "\n";
    }
    myData << "\n";
    myDataI << "\n";
    myDataII << "\n";
    myDataIII << "\n";
    myDataIV << "\n";
#ifdef TWO_PICKS
    myDataV << "\n";
#endif 

#ifndef TWO_PICKS
    myBest << "\n";
#else
    file_best_37 << "\n";
    file_best_112 << "\n";
#endif //

#ifdef KNOWN_MAXIM    
#ifdef TWO_ATTRACTORS_FUNCTION
    double firstBestMean =0;
    double secondBestMean = 0;
    double minBestMean = 0;
    double meanBestMean = 0;
    //myData << " mean best " << "\n";
    for(unsigned long i = 0; i < nr_runs; i++){
	myData << i << "\t"<< bestFitnessVector[i][0] << "\t"<< bestFitnessVector[i][1] << "\n";
	if(bestFitnessVector[i][0] != 0)
	    firstBestMean += bestFitnessVector[i][0];
	else firstBestMean += sizeMCMCstring * sampleSize;
	if(bestFitnessVector[i][1] != 0)
	    secondBestMean+= bestFitnessVector[i][1];
	else secondBestMean+= sizeMCMCstring * sampleSize;
	meanBestMean += abs(bestFitnessVector[i][0] - bestFitnessVector[i][1]);
	if(bestFitnessVector[i][0]>bestFitnessVector[i][1] && bestFitnessVector[i][1] != 0)
	    minBestMean += bestFitnessVector[i][1];
	else if(bestFitnessVector[i][0] != 0)
	    minBestMean += bestFitnessVector[i][0];
	else minBestMean += sizeMCMCstring * sampleSize;
    }
    myData << " means  " << firstBestMean/nr_runs << "\t" << secondBestMean/nr_runs << 
	"\t" << minBestMean/nr_runs << "\t" << meanBestMean/nr_runs << "\n";
    meanFitness = minBestMean/nr_runs;
    diffFitness = meanBestMean/nr_runs;
#else 
#ifndef BIVARIATE 
    double firstbestFitness = 0;
    double firstExists = 0, meanFirstExists = 0; 
    //myData << " mean first best " << "\n";
    for(unsigned long i = 0; i < nr_runs; i++)
	if(bestFitnessVector[i] != 0){
	    //    myData << i << "\t"<< bestFitnessVector[i] << "\n";
	    firstbestFitness += bestFitnessVector[i];
	    firstExists ++;
	    meanFirstExists += bestFitnessVector[i];
	} else {
	    //myData << i << "\t"<< sizeMCMCstring * sampleSize << "\n";
	    firstbestFitness += sizeMCMCstring * sampleSize;
	}
    
    meanFitness = firstbestFitness/nr_runs;
    meanFirstExists /= firstExists;

    double std_best = 0;
    double std_exists = 0;
    for(unsigned long i = 0; i < nr_runs; i++)
	if(bestFitnessVector[i] != 0){
	    std_best += pow(bestFitnessVector[i] - meanFitness,2.0);
	    std_exists += pow(bestFitnessVector[i] - meanFirstExists,2.0);
	} else {
	    std_best += pow(sizeMCMCstring * sampleSize - meanFitness,2.0);
	}

    std_best = pow(std_best/(double)nr_runs,0.5);
    std_exists = pow(std_exists/(double)nr_runs,0.5);

    myData << "#" << sizeMCMCstring << "\t" <<  sampleSize << "\t" << meanFitness << "\t" << std_best << "\t" 
	   << meanFirstExists << "\t" << std_exists << "\n";
#else
#endif //BIVARIATE
#endif
#else //known maxim
#endif // known maxim    

#ifndef TWO_PICKS
    myBest.flush();
    myBest.close();
#else
    file_best_37.flush();
    file_best_37.close();

    file_best_112.flush();
    file_best_112.close();
#endif 

    myData.flush();
    myData.close();

    myDataI.flush();
    myDataI.close();

    myDataII.flush();
    myDataII.close();

    myDataIII.flush();
    myDataIII.close();

    myDataIV.flush();
    myDataIV.close();

#ifdef TWO_PICKS
    myDataV.flush();
    myDataV.close();
#endif //two picks
  }

  delete[] performance;
  return 1;
}
#else //MIXTURE_BIVARIATE
double binaryPopulation::printDiversityGNUFILE(){

  cout <<"print diversity \n";
  
  double* performance = new double[diversityMeasure];
  
  
  for(unsigned long gen = 0; gen <= (unsigned long)(sampleSize_Const/SEE_SAMPLE); gen++){
    
    char* fileN = new char[50];

    fileN = strcpy (fileN, "gnuFILE_mut\0");

    if(gen == 0)
      fileN = strcat (fileN, "_0.txt\0");
    else if(gen == 1) 
      fileN = strcat (fileN, "_500.txt\0");
    else if(gen == 2)
      fileN = strcat (fileN, "_1000.txt\0");
    else if(gen == 3)
      fileN = strcat (fileN, "_1500.txt\0");
    else if(gen == 4)
      fileN = strcat (fileN, "_2000.txt\0");
    else {
      fileN = strcat (fileN, "_2500.txt\0");
    }

    cout << "file name" << fileN << " because " << gen << "\n";

    ofstream myData(fileN);
    delete[] fileN;

    if(!myData.is_open()){
	cout <<"Error in opening the diversity file\n";
	exit(1);
    }

    //#ifdef MULTIPLE_RUNS
#ifndef IsRecombination
    myData << "#Mutation - Diversity mean " << sizeGENOME << "\t" << sampleSize << "\t" << MCMCstring[0]->mutation <<"\n";
#else
    //myData << "#\n Recombination - Diversity mean " << sizeGENOME << "\t" << sampleSize << "\t" 
//	   << MCMCstring[0]->PROB_nPOINT <<"\n";
    myData << "#\n Recombination - Diversity mean " << sizeGENOME << "\t" << sampleSize << "\t" 
#ifdef BINARY_RECOMB
	   << UNIF_RECOMB_PARAM << "\t"<< " and mutation rate " << MCMCstring[0]->mutation <<"\n";
#else
	   << PARENT_CENTRIC << "\t"<< " and mutation rate " << MCMCstring[0]->mutation <<"\n";
#endif //binary recomb
#endif
    //#endif //multiple runs
    
    //if (diversityTrueDistribution[j][k] != 0 & globalDiversityRightSolutions[i][j][k] != 0) {
    /////bin hight at a point in sapce
    for(unsigned long j = 0; j < diversityMeasure; j++)
	performance[j] = 0;
    
    double tempDiv = 0;
    
    //measure the mean of the histogram
    for(unsigned long j= 0; j < diversityMeasure; j++){
      double performance = 0;
      for(unsigned long i = 0; i < nr_runs; i++)
	if(globalDiversity[gen][i][j] != 0 && globalDiversityRightSolutions[gen][i][j] != 0)
	  performance += (double)globalDiversity[gen][i][j]/(double)globalDiversityRightSolutions[gen][i][j];	  
      performance /= (double)nr_runs;
      
      double st_dev_perf = 0;
      for(unsigned long i = 0; i < nr_runs; i++){
	/////bin hight at a point in sapce
	if(globalDiversity[gen][i][j] != 0 && globalDiversityRightSolutions[gen][i][j] != 0){
	  tempDiv = globalDiversity[gen][i][j]/globalDiversityRightSolutions[gen][i][j];
		st_dev_perf += pow(performance - tempDiv,2.0);
	}
      }

      //split the data "j" in dimensions
      int k = sizeGENOME - 1;
      double rest[sizeGENOME];
      double cit;
      while(k >= -1){
	cit = ((unsigned long)j)/100;
	rest[k] = j - cit;
	j = cit;
	k--;
      }
      
      if(rest[sizeGENOME - 1] == 0) myData << "\n";
      
      for(k =0; k < sizeGENOME; k++)
	myData << rest[k] << "\t";

      myData << performance << "\t" << sqrt(st_dev_perf/(double)nr_runs) << "\t" << trueDistribution[j] << "\n";
    }
    myData << "\n";

    myData.flush();
    myData.close();
  }

  delete[] performance;
  return 1;
}
#endif //MIXTURE_BIVARIATE
#endif //histogram
#endif // mixture reals

double binaryPopulation::printDiversityRightSolutions(){
    ofstream myData(gnuplotDRSFILE, ios::out|ios::app);
    if(!myData.is_open()){
	cout <<"Error in opening the diversity file\n";
	return -1;
    }
    for(unsigned long j= 0; j < nr_runs; j++){
	cout << "RightSolution "<< j << " is " << (double)globalRightSolutions[j] << "\n";
	myData << j <<"\t"<<((double)globalRightSolutions[j]/(double)(sampleSize_Const*sizeMCMCstring))<<"\n";
    } 
    myData.flush();
    myData.close();
    return 1;
}
#endif //SINGLE_vs_MULTIPLE
#endif //TRACE_RUN 

#ifndef MIXTURE_REALS
double binaryPopulation::scoreDiversity(){
    double temp = 0;
#ifndef SINGLE_vs_MULTIPLE
    for(unsigned long j= 0; j < diversityMeasure; j++){ 
      if(diversity[j] != 0) {
	temp++;
      } else cout << "diversity is 0 for j = " << j << "\n";
      globalDiversity[runs-1][j] = diversity[j];
      globalDiversityRightSolutions[runs-1][j] = diversityRightSolutions[j];
    }
    temp /= (double)diversityMeasure;
#else 
    unsigned long i = (unsigned long)(log(sizeMCMCstring)/log(2.0)); 
    for(unsigned long j= 0; j < diversityMeasure; j++){ 
      if(diversity[j] != 0) {
	temp++;
      }
      globalDiversity[i][runs-1][j] = diversity[j];
      globalDiversityRightSolutions[i][runs-1][j] = diversityRightSolutions[j];
    }
    temp /= (double)diversityMeasure;
#endif //SINGLE_vs_MULTIPLE
    return temp;
}
#endif //mixture reals

#ifdef TRAP_LIKE_FUNCTION
double binaryPopulation::scoreDiversityDetail(){
    double temp = 0;
    for(unsigned long j= 0; j < sizeGENOME/BLOCKsize; j++) 
      for(unsigned long k = 0; k < BLOCKsize+1; k++){
	if(diversityDetail[j][k] != 0) {
	  temp++;
	}
	globalDiversityDetail[j][k] += diversityDetail[j][k] / (sampleSize_Const * sizeMCMCstring);
      }
    temp /= (double) (sizeGENOME/BLOCKsize)*(BLOCKsize+1);
    return temp;
}

double binaryPopulation::printDiversityDetailGNUFILE(){
    ofstream myData(gnuplotDDFILE, ios::out|ios::app);
    if(!myData.is_open()){
	cout <<"Error in opening the diversity file\n";
	return -1;
    }
    myData << "\n sizeGENOME " << sizeGENOME << "\t BLOCKsize " << BLOCKsize << "\n";
    for(unsigned long j= 0; j < sizeGENOME/BLOCKsize; j++) {
	myData << j;
	for(unsigned long k = 0; k < BLOCKsize+1; k++)
	    myData <<"\t"<<((double)globalDiversityDetail[j][k]/(double)(nr_runs));
	myData <<"\n";
    }
  myData.flush();
  myData.close();
  return 1;
}
#endif //TRAP_LIKE_FUNCTION

double binaryPopulation::f(double temp)
{
#ifndef MIXTURE_REALS
#ifdef THRESHOLD
#ifndef THRESHOLD_DISTRIBUTION
#ifndef HALF_MAX_FITNESS
#ifdef TRAP_LIKE_FUNCTION
  if(temp > (diversityMeasure -1)* b / a) return temp;
#else
#ifdef BINOMIAL 
  if(temp > (diversityMeasure -1)* (BLOCKsize-1) / BLOCKsize) return temp;
  return 0;
#endif //BINOMIAL
#endif //trap like function
#else
  if(temp >= (diversityMeasure -1) / 2) return temp;
#endif //HALF_MAX_FITNESS
#else
  return temp;
#endif //THRESHOLD_DISTRIBUTION
#else 
  return temp;
#endif //Thresold
#else //mixture reals
  if(temp > 1.0) return temp;
#endif //mixture reals
  return -1;
}

double binaryPopulation::averageGeneration(){
    double temp = 0;
    //cout << "Average to ";
    for(unsigned long i = 0; i<sizeMCMCstring; i++){
	temp += MCMCstring[i] ->fitness();
	//cout << MCMCstring[i] -> value << "\t";
    }
    temp /= sizeMCMCstring;
    //cout << " is " << temp << "\n";
    //averageFitness += temp;
    return temp;
}

double binaryPopulation::geometricGeneration(){
    double temp = 1;
    //cout << "Geometric average ";
    for(unsigned long i = 0; i<sizeMCMCstring; i++){
	temp *= MCMCstring[i] ->getValue();
	//cout << MCMCstring[i] -> value << "\t";
    }
    //temp = sqrt(temp);
    //cout << " is " << temp << "\n";
    temp = pow(temp,1.0/(double)sizeMCMCstring);
    //averageFitness = temp;
    return temp;
}

double binaryPopulation::bestFitnessP(){
    double count = 0;
    double count2 = 0;

    for(unsigned long i = 0; i<sizeMCMCstring; i++){
	if((MCMCstring[i]) -> stop() == 2) 
	    count2 = 2;  
	if((MCMCstring[i]) -> stop() == 1) 
	    count = 1;
    }
    return count2+count;
}

//----------- evolution of mixed strategies --------
/*
  void binaryPopulation::evolution(unsigned long myGenerations){
  int stopBLOCKS = 0; 
  setSeed();
  while(generation != myGenerations && stopBLOCKS!= 1)
    {	
      //distribution program
      unsigned int tempValue = generation % (sizeGENOME / 2 + 1);
      //unsigned int tempValue = generation % (sizeGENOME + 1);
      if(tempValue < sizeGENOME / 2) {
	noCoupling();
      } else //temperingCoupling(2);
	elitistRandomCoupling(2);
      //randomCoupling(2);
      
      for(unsigned long i = 0; i<sizeMCMCstring; i++) {
	if((MCMCstring[i]) -> stopBLOCKS() == 1) 
	  stopBLOCKS = 1;
      }	
      //print();
    }	
}
*/

#ifdef ACCEPTANCE_RATIO
//include si average
void binaryPopulation::printAcceptance(ofstream& myAcceptance){
    myAcceptance << sizeGENOME << "\t" << sizeMCMCstring << "\t" << 
	generation << "\t" << runs - 1 << "\t" << acceptanceAverage() 
		 << "\t" << geometricGeneration() << "\t" << 
	averageGeneration() << "\t" << acceptance_ratio[runs-1] << "\t";
    myAcceptance << averageFitness/(generation+1) << "\t";
#ifdef DISTANCE
//	    Kullback[generation] += temp;
    //if(generation % (SEE_SAMPLE) == 0){
#ifndef MIXTURE_REALS
	    unsigned long temp[sizeGENOME];
#else
	    double temp[sizeGENOME];
#endif
	    MCMCstring[0]->getMax(temp);
	    myAcceptance << HammingDistance(temp);	    
//simul = insert_simulation(simul,generation,distH);
	    //show_simulations(simul);
	    //write_to_file_all(simul_all_proc_right,simul,NULL);
	    //}
#endif
    myAcceptance << "\n";
    if (generation >= burn_in + sampleSize*interSamplesize){
	myAcceptance << "\n";
	if(runs < nr_runs)
	    myAcceptance << "\n";
    }
}


double binaryPopulation::acceptanceAverage(){
    double temp = 0;
    for(unsigned long i = 0; i<sizeMCMCstring; i++)
	temp += MCMCstring[i] ->accepted;
    temp /= sizeMCMCstring;
    return temp;
}
#endif

#ifndef FILE_BOOKEEPING
void binaryPopulation::printDiversity(){
#ifndef MIXTURE_REALS
    unsigned long* temp = new unsigned long[sizeGENOME];
#else
    double* temp = new double[sizeGENOME];
#endif
    for(unsigned long i = 0; i < sizeMCMCstring; i++){
      diversityFile << MCMCstring[i]->getValue() << "\t";
      MCMCstring[i]->getGENOME(temp);
      for(unsigned long j = 0; j < sizeGENOME; j++)
	  diversityFile << temp[j];
      diversityFile << "\n";
      if (generation >= burn_in + sampleSize*interSamplesize){
	  diversityFile << "\n";
	  if(runs < nr_runs)
	      diversityFile << "\n";
      }
    }
    delete temp;
} 

void binaryPopulation::printDiversityPopulation(){
    double tempValue = 0;
    for(unsigned long j = 0; j < sizeMCMCstring; j++){
	  tempValue += MCMCstring[j] -> getValue();
    }
    tempValue /= sizeMCMCstring;
    diversityFile << tempValue << "\n";
    if (generation >= burn_in + sampleSize*interSamplesize){
	diversityFile << "\n";
	if(runs < nr_runs)
	    diversityFile << "\n";
    }
} 

#endif //filw bookiping

//----------REC_MUT RANDOM_WALK-------
//numai pentru recombinare cu 2 chainuri
#ifdef RECOMB_MUT
void binaryPopulation::evolution(){
  unsigned long** table = new unsigned long*[1];
  for(unsigned long i = 0; i< 1; i++)
    table[i] = new unsigned long[sizeCoupling];

  cout << "RW recomb vs. mut with sizeCoupling "<< sizeCoupling<< "\n";

#ifdef IsMutation 
  ofstream myFile("rw_file_mut.data");
#else
  ofstream myFile("rw_file_recomb.data");
#endif		
  if(!myFile.is_open()){
    cout << "Error opening file rw \n";
    return;
  }

#ifndef MIXTURE_REALS 
  double tempDistance[diversityMeasure];
  double myDistance[diversityMeasure];
  double count[diversityMeasure];
  for(unsigned long i = 0; i < diversityMeasure; i++){
      myDistance[i] = 0;
      tempDistance[i] = 0;
      count[i] = 0;
  }
#else
  double tempDistance[diversityMeasure][diversityMeasure];
  double myDistance[diversityMeasure][diversityMeasure];
  double count[diversityMeasure][diversityMeasure];
  for(unsigned long i = 0; i < diversityMeasure; i++)
      for(unsigned long j = 0; j < diversityMeasure; j++)
      {
	  myDistance[i][j] = 0;
	  tempDistance[i][j] = 0;
	  count[i][j] = 0;
      }

  // the point with which it is compared
  unsigned long temp1[sizeGENOME*BLOCKsize];
#endif

#ifndef MIXTURE_REALS	
	unsigned long temp[sizeGENOME];
#else
	double temp[sizeGENOME];

#ifdef BOTH_GOOD
	double temp_recomb[sizeGENOME];

#ifndef TWO_PARENTS
	//if(sizeCoupling == 3){ 
	  double temp_recomb2[sizeGENOME];
	  //}
#endif
#endif
#endif //mixture reals

	  //cout << "Diversity measure " << diversityMeasure << "\t nr_runs" << nr_runs << "\n";
  while(runs < nr_runs){

    runs++;
    
    MCMCstring[0]->setMax0();
    MCMCstring[1]->setMax1();
#ifndef TWO_PARENTS
    //if(sizeCoupling == 3){
      MCMCstring[2]->setMax2();
      table[0][2] = 2;
      // }
#endif
    table[0][0] = 0;
    table[0][1] = 1;

    MCMCstring[0]->getGENOME(temp);
    double tempValue = MCMCstring[0]->getValue();

#ifdef BOTH_GOOD
    MCMCstring[1]->getGENOME(temp_recomb);
    double tempValue_recomb = MCMCstring[1]->getValue();

#ifndef TWO_PARENTS
    //if(sizeCoupling == 3){
      MCMCstring[2]->getGENOME(temp_recomb2);
      double tempValue_recomb2 = MCMCstring[2]->getValue();      
      //}
#endif
#endif // both good    

#ifndef MIXTURE_REALS	
    unsigned long div = (unsigned long)HammingDistance(temp);
    tempDistance[div] += tempValue;
    //cout << " gen Val" << tempValue << "\n";
    count[div]++;
    //cout << "div " << div << " generation " << generation << "\n";
#else
    unsigned long dim1 = (unsigned long) (temp[0] * hist_points/scale);
    unsigned long dim2 = (unsigned long) (temp[1] * hist_points/scale);
    tempDistance[dim1][dim2] += tempValue;
    count[dim1][dim2]++;

    //cout << "at " << generation << " dist("<< temp[0]*hist_points/scale<< "," <<
    //  temp[1]*hist_points/scale << ") val=" 
    //	 << tempValue << " with " << tempDistance[dim1][dim2]<< " no " << count[dim1][dim2] << "\n";

#ifdef BOTH_GOOD
    dim1 = (unsigned long) (temp_recomb[0] * hist_points/scale);
    dim2 = (unsigned long) (temp_recomb[1] * hist_points/scale);
    tempDistance[dim1][dim2] += tempValue_recomb;
    count[dim1][dim2]++;

    //cout << "at " << generation << " dist("<< temp_recomb[0]*hist_points/scale<< "," <<
    //  temp_recomb[1]*hist_points/scale << ") val=" 
    //	 << tempValue_recomb << " with " << tempDistance[dim1][dim2]<< " no " << count[dim1][dim2] << "\n";

#ifndef TWO_PARENTS
    if(sizeCoupling == 3){
      dim1 = (unsigned long) (temp_recomb2[0] * hist_points/scale);
      dim2 = (unsigned long) (temp_recomb2[1] * hist_points/scale);
      tempDistance[dim1][dim2] += tempValue_recomb2;
      count[dim1][dim2]++;
    }
#endif 
    //cout << "at " << generation << " dist("<< temp_recomb2[0]*hist_points/scale<< "," <<
    //  temp_recomb2[1]*hist_points/scale << ") val=" 
    //	 << tempValue_recomb2 << " with " << tempDistance[dim1][dim2]<< " no " << count[dim1][dim2] << "\n";
#endif // both good

#endif
    
    
    //cout << "to binary first genome";
    //convertBinary(MCMCstring[0]->getGENOME(temp),temp1,sizeGENOME);
    //for(unsigned long i = 0; i < sizeGENOME*BLOCKsize; i++)
    //cout << temp1[i];
    //cout << "\n";
    //
    //cout << "to binary second genome";
    //convertBinary(MCMCstring[1]->getGENOME(temp),temp1,sizeGENOME);
    //for(unsigned long i = 0; i < sizeGENOME*BLOCKsize; i++)
    //	cout << temp1[i];
    //cout << "\n";
    
    while(generation <= burn_in + sampleSize*interSamplesize){	
	//randomCOUPLING(table);
	proposedGeneration(table,1);
	//acceptanceProposed(table,dimension);

	proposedMCMCstring[0]->getGENOME(temp);
	//MCMCstring[1]->getGENOME(temp1);
	cout << "gen" << generation << " genome " << temp[0] << "\t" << temp[1] << "\n";
	//cout << "to binary ";
	//convertBinary(temp,temp1,sizeGENOME);
	//for(unsigned long i = 0; i < sizeGENOME*BLOCKsize; i++)
	//    cout << temp1[i];
	//cout << "\n";
	//cout << "gen" << generation << " genome " << temp1[0] << "\t" << temp1[1] << "\n";
 
	tempValue = proposedMCMCstring[0]->getValue();

#ifdef BOTH_GOOD
	proposedMCMCstring[1]->getGENOME(temp_recomb);
	//MCMCstring[1]->getGENOME(temp1);
	cout << "gen" << generation << " genome " << temp_recomb[0] << "\t" << temp_recomb[1] << "\n";
	//cout << "to binary ";
	//convertBinary(temp,temp1,sizeGENOME);
	//for(unsigned long i = 0; i < sizeGENOME*BLOCKsize; i++)
	//    cout << temp1[i];
	//cout << "\n";
	//cout << "gen" << generation << " genome " << temp1[0] << "\t" << temp1[1] << "\n";
 
	tempValue_recomb = proposedMCMCstring[1]->getValue();

#ifndef TWO_PARENTS
	if(sizeCoupling == 3){
	  proposedMCMCstring[2]->getGENOME(temp_recomb2);
	  //MCMCstring[1]->getGENOME(temp1);
	  cout << "gen" << generation << " genome " << temp_recomb2[0] << "\t" << temp_recomb2[1] << "\n";
	  //cout << "to binary ";
	  //convertBinary(temp,temp1,sizeGENOME);
	  //for(unsigned long i = 0; i < sizeGENOME*BLOCKsize; i++)
	  //    cout << temp1[i];
	  //cout << "\n";
	  //cout << "gen" << generation << " genome " << temp1[0] << "\t" << temp1[1] << "\n";
	  
	  tempValue_recomb2 = proposedMCMCstring[2]->getValue();
	  }
#endif
#endif // both good

#ifndef MIXTURE_REALS	
	unsigned long div = (unsigned long)HammingDistance(temp);
	tempDistance[div] += tempValue;
	//cout << " gen Val" << tempValue << "\n";
	count[div]++;
	//cout << "div " << div << " generation " << generation << "\n";
#else
	unsigned long dim1 = (unsigned long) (temp[0] * hist_points/scale);
	unsigned long dim2 = (unsigned long) (temp[1] * hist_points/scale);
	tempDistance[dim1][dim2] += tempValue;
	count[dim1][dim2]++;

	//cout << "at " << generation << " dist("<< temp[0]*hist_points/scale<< "," <<
	//    temp[1]*hist_points/scale << ") val=" 
	//     << tempValue << " with " << tempDistance[dim1][dim2]<< " no " << count[dim1][dim2] << "\n";

#ifdef BOTH_GOOD
	dim1 = (unsigned long) (temp_recomb[0] * hist_points/scale);
	dim2 = (unsigned long) (temp_recomb[1] * hist_points/scale);
	tempDistance[dim1][dim2] += tempValue_recomb;
	count[dim1][dim2]++;

	//cout << "at " << generation << " dist("<< temp_recomb[0]*hist_points/scale<< "," <<
	//    temp_recomb[1]*hist_points/scale << ") val=" 
	//     << tempValue_recomb << " with " << tempDistance[dim1][dim2]<< " no " << count[dim1][dim2] << "\n";

#ifndef TWO_PARENTS
	if(sizeCoupling == 3){
	  dim1 = (unsigned long) (temp_recomb2[0] * hist_points/scale);
	  dim2 = (unsigned long) (temp_recomb2[1] * hist_points/scale);
	  tempDistance[dim1][dim2] += tempValue_recomb2;
	  count[dim1][dim2]++;
	  
	  //cout << "at " << generation << " dist("<< temp_recomb2[0]*hist_points/scale<< "," <<
	  // temp_recomb2[1]*hist_points/scale << ") val=" 
	  //    << tempValue_recomb2 << " with " << tempDistance[dim1][dim2]<< " no " << count[dim1][dim2] << "\n";
	}
#endif
#endif // both good

#endif

#ifndef RANDOM_WALK
	if(MCMCstring[0] -> acceptance(temp,tempValue) == 1){
#endif
	  MCMCstring[0] -> nextGeneration(temp,tempValue,sizeGENOME);
	  MCMCstring[0] -> setACCEPTED(1);
#ifndef RANDOM_WALK
	} else{
	  MCMCstring[0] -> nextGenerationBlank();
	  MCMCstring[0] -> setACCEPTED(0);
	}
#endif
	 
	postProcessing();
	generation++;
    }
    
#ifndef MIXTURE_REALS	
    for(unsigned long i = 0; i < diversityMeasure; i++){
	if(count[i] != 0) {
	    myDistance[i] += tempDistance[i]/count[i];
	    cout << " Distance " << tempDistance[i]/count[i] << " i " << i << 
	    	" runs " << runs << "\n"; 
	}
	count[i] = 0;
	tempDistance[i] = 0;
    }
#else
    //cout << "run" << runs << "\n";
    for(unsigned long i = 0; i < diversityMeasure; i++)
	for(unsigned long j = 0; j < diversityMeasure; j++)
	{
	    if(count[i][j] != 0) {
#ifndef UNIFORM_DISTRIBUTION
		myDistance[i][j] += tempDistance[i][j]/count[i][j];
#else
		myDistance[i][j] += count[i][j] ;
#endif //UNIFORM_DISTRIBUTION
		//cout << " Distance " << tempDistance[i][j] << "\t" << count[i][j] << " i " << i << 
		//  ", j " << j << " runs " << runs << " myDistance" << myDistance[i][j] << "\n"; 
	    }
	    count[i][j] = 0;
	    tempDistance[i][j] = 0;
	}
#endif //mixture reals

    if(runs != nr_runs) reset();
    
  }
  
#ifndef MIXTURE_REALS	
  for(unsigned long i = 0; i < diversityMeasure; i++){
      myFile << i << "\t" << myDistance[i]/runs << "\n";
  }
#else
  for(unsigned long i = 0; i < diversityMeasure; i++){
    int one = 0;
      for(unsigned long j = 0; j < diversityMeasure; j++)
	if(myDistance[i][j] != 0){
	  myFile << i << "\t" << j << "\t" << myDistance[i][j]/runs << "\n";
	  one = 1;
      }
     if(one == 1)  myFile << "\n";
  }
#endif // mixture reals
  
  myFile.flush();
  myFile.close();
}
#endif
// -------- No coupling - individual evolution ---------

#ifdef NO_COUPLING
void binaryPopulation::noCoupling(){
    for(unsigned long i = 0; i<sizeMCMCstring ; i++){
      MCMCstring[i] -> nextGeneration();		
    }
}

#ifndef SINGLE_vs_MULTIPLE
#ifdef MULTIPLE_RUNS
void binaryPopulation::evolution(ofstream& myMultiFile){
#else		
void binaryPopulation::evolution(){
#endif
  cout << "NO COUPLING \n";

#ifndef MIXTURE_REALS
  unsigned long tempGENOME[sizeGENOME];
#else
  double tempGENOME[sizeGENOME];
#endif //mixture reals

#ifndef MULTIPLE_RUNS
#ifndef TRACE_RUN
#ifdef HISTOGRAM
  ofstream myFileFG("err_mut_fit_secondGen");		
  if(!myFileFG.is_open()){
      cout << "Error opening file for data err_mut_dist_firstGen \n";
    return;
  }

  ofstream myFileFit("err_mut_fit_firstGen");		
  if(!myFileFit.is_open()){
      cout << "Error opening file for data err_mut_fit_firstGen \n";
    return;
  }

 ofstream myFile(fileData);		
  if(!myFile.is_open()){
      cout << "Error opening file write in " << fileData << " \n";
    return;
  }

  ofstream myTempFile("temporar.txt");		
  if(!myTempFile.is_open()){
      cout << "Error opening file for data " << fileData << "\n";
    return;
  }

  ofstream myFileHash("err_mut_hash_firstGen");		
  if(!myFileHash.is_open()){
      cout << "Error opening file for data err_mut_hash_firstGen \n";
      return;
  }
#else //histogram
  
  int records;
  double mse, cor, s1, s2, jump1, jumpGen;
  double tempPrev[sizeMCMCstring][sizeGENOME];
  double tempValuePrev[sizeMCMCstring];
  if(nr_runs == 1)
    records = (int)((double)sampleSize_Const/ (double)SEE_SAMPLE) + 1;
  else records = nr_runs;

  double averageFitnessIndivEach = 0;
  double sqrAverageFitnessIndivEach = 0;
  double averageVariablesIndivEach[sizeGENOME];
  double sqrAverageVariablesIndivEach[sizeGENOME];
  for(int j = 0; j < sizeGENOME; j++){
    averageVariablesIndivEach[j] = 0;
    sqrAverageVariablesIndivEach[j] = 0;
  }

  double averageFitnessIndivRuns[records];
  double averageFitnessRuns[records];
  double averageVariablesIndivRuns[records][sizeGENOME];

  double corVariablesIndivRuns[records][sizeGENOME * (sizeGENOME - 1) / 2];
  double sqrCorVariablesIndivRuns[records][sizeGENOME * (sizeGENOME - 1) / 2];

  double sqrAverageFitnessIndivRuns[records];
  double sqrAverageFitnessRuns[records];
  double sqrAverageVariablesIndivRuns[records][sizeGENOME];
  for(int j = 0; j < records; j++){
    averageFitnessIndivRuns[j] = 0;
    averageFitnessRuns[j] = 0;
    for(int k = 0; k < sizeGENOME; k++)
      averageVariablesIndivRuns[j][k] = 0;
    sqrAverageFitnessIndivRuns[j] = 0;
    sqrAverageFitnessRuns[j] = 0;
    for(int k = 0; k < sizeGENOME; k++)
      sqrAverageVariablesIndivRuns[j][k] = 0;
    for(int k = 0; k < sizeGENOME * (sizeGENOME - 1) / 2; k++){
      corVariablesIndivRuns[j][k] = 0;
      sqrCorVariablesIndivRuns[j][k] = 0;
    }
  }
  
  double averageFitnessIndiv[records][sizeMCMCstring];
  double averageVariablesIndiv[records][sizeMCMCstring][sizeGENOME];
  double corVariablesIndiv[records][sizeMCMCstring][sizeGENOME * (sizeGENOME - 1) / 2];
  double sqrCorVariablesIndiv[records][sizeMCMCstring][sizeGENOME * (sizeGENOME - 1) / 2];
  double averageJumpIndiv[records][sizeMCMCstring];
  double sqrAverageJumpIndiv[records][sizeMCMCstring];
  double averageJumpIndivRuns[records];
  double sqrAverageJumpIndivRuns[records];
  double averageExpJumpIndivRuns[records];
  double sqrAverageExpJumpIndivRuns[records];
  double averageJump[records];
  double sqrAverageJump[records];
  for(int j = 0; j < records; j++){
    for(int k = 0; k < sizeMCMCstring; k++){
      averageFitnessIndiv[j][k] = 0;
      for(int i = 0; i < sizeGENOME; i++)
	averageVariablesIndiv[j][k][i] = 0;
      for(int i = 0; i < sizeGENOME * (sizeGENOME - 1) / 2; i++){
	corVariablesIndiv[j][k][i] = 0;
	sqrCorVariablesIndiv[j][k][i] = 0;
      }
      averageJumpIndiv[j][k] = 0;
      sqrAverageJumpIndiv[j][k] = 0;
    }
    averageJumpIndivRuns[j] = 0;
    sqrAverageJumpIndivRuns[j] = 0;
    averageExpJumpIndivRuns[j] = 0;
    sqrAverageExpJumpIndivRuns[j] = 0;
    averageJump[j] = 0;
    sqrAverageJump[j] = 0;
  }

  double sqrAverageFitnessIndiv[records][sizeMCMCstring];
  double sqrAverageVariablesIndiv[records][sizeMCMCstring][sizeGENOME];
  for(int j = 0; j < records; j++)
    for(int k = 0; k < sizeMCMCstring; k++){
      sqrAverageFitnessIndiv[j][k] = 0;
      for(int i = 0; i < sizeGENOME; i++)
	sqrAverageVariablesIndiv[j][k][i] = 0;
    }
  

  double averageFitnessIndivAll[sizeMCMCstring];
  double averageVariablesIndivAll[sizeMCMCstring][sizeGENOME];
  for(int j = 0; j < sizeMCMCstring; j++){
    averageFitnessIndivAll[j] = 0;
    for(int i = 0; i < sizeGENOME; i++)
      averageVariablesIndivAll[j][i] = 0;
  }

  double sqrAverageFitnessIndivAll[sizeMCMCstring];
  double sqrAverageVariablesIndivAll[sizeMCMCstring][sizeGENOME];
  for(int j = 0; j < sizeMCMCstring; j++){
    sqrAverageFitnessIndivAll[j] = 0;
    for(int i = 0; i < sizeGENOME; i++)
      sqrAverageVariablesIndivAll[j][i] = 0;
  }
  
  ofstream* myFileAllSplit = new ofstream[records];
  char nameFileSplit[30] = "outMutAllChainSplit";		
  for(int i = 0; i < records; i++){
    nameFileSplit[19] = 48 + (int)(i / 10);
    nameFileSplit[20] = 48 + i - 10 * (int)(i / 10);
    nameFileSplit[21] = '\0'; 
    myFileAllSplit[i].open(nameFileSplit);
    if(!myFileAllSplit[i].is_open()){
      cout << "Error opening file for data output_mut_all_chains_splits \n";
      return;
    } else {
      //write the first row
      /*for(int k = 0; k < sizeMCMCstring; k++){
	for(int j = 0; j < sizeGENOME;j++)
	  myFileAllSplit[i] << "(C" << k << ",V" << j <<")\t(C" << k << ",MV" << j <<")\t";
	myFileAllSplit[i] << "(C"<< k << ",F)\t(C"<< k << ",MF)\t";
      }
      myFileAllSplit[i] << "(A) \t (MA)\n";
      */
    }
  }

  /*nameFileSplit[19] = 'A'; nameFileSplit[20] = 'l';  nameFileSplit[21] = 'l';  nameFileSplit[22] = '\0'; 
  ofstream tempFileAllSplit(nameFileSplit);		
  if(!tempFileAllSplit.is_open()){
    cout << "Error opening file for data err_mut_dist_firstGen \n";
    return;
    }*/
  
  ofstream myFileAllNoSplit("outMutAllChainNoSplit");		
  if(!myFileAllNoSplit.is_open()){
      cout << "Error opening file for data err_mut_fit_firstGen \n";
    return;
  } else {
      //write the first row
    /*  for(int i = 0; i < records; i++){
      for(int k = 0; k < sizeMCMCstring; k++){
	for(int j = 0; j < sizeGENOME;j++)
	  myFileAllNoSplit << "(R" << i  << ",C" << k << ",V" << j << ")\t(R" << i << ",C" << k << ",MV" << j << ")\t";
	myFileAllNoSplit << "(R" << i << ",C"<< k << ",F)\t(R" << i << ",C"<< k << ",MF)\t";
      }
      myFileAllNoSplit << "(R" << i << ",A)\t(R" << i << ",MA)\t";
    }
    myFileAllNoSplit << "\n"; 
    */
  }

  ofstream* myFile1ChainSplit = new ofstream[records];
  nameFileSplit[6] = 'O'; nameFileSplit[7] = 'n'; nameFileSplit[8] = 'e'; 		
  for(int i = 0; i < records; i++){
    nameFileSplit[19] = 48 + (int)(i / 10);
    nameFileSplit[20] = 48 + i - 10 * (int)(i / 10);
    nameFileSplit[21] = '\0'; 
    myFile1ChainSplit[i].open(nameFileSplit);
    if(!myFile1ChainSplit[i].is_open()){
      cout << "Error opening file for data output_mut_all_chains_splits \n";
      return;
    } else {
      /*for(int j = 0; j < sizeGENOME; j++)
	myFile1ChainSplit[i] << "(V" << j << ")\t(MV" << j << ")\t";
	myFile1ChainSplit[i] << "(F)\t(MF)\n";*/
    }
  }

  ofstream myFile1ChainNoSplit("outMutOneChainNoSplit");		
  if(!myFile1ChainNoSplit.is_open()){
      cout << "Error opening file for data err_mut_fit_firstGen \n";
    return;
  } else {
    /*for(int j = 0; j < sizeGENOME; j++)
      myFile1ChainNoSplit << "(V" << j << ")\t(MV" << j << ")\t";
    myFile1ChainNoSplit << "(F)\t(MF)\n";
    */
  }


  /*nameFileSplit[19] = 'A'; nameFileSplit[20] = 'l';  nameFileSplit[21] = 'l';  nameFileSplit[22] = '\0'; 
  fstream tempFile1ChainSplit(nameFileSplit);		
  if(!tempFile1ChainSplit.is_open()){
      cout << "Error opening file for data err_mut_dist_firstGen \n";
    return;
    }*/
 
  ofstream myFileHash("outMutHashAllChainNoSplit");		
  if(!myFileHash.is_open()){
    cout << "Error opening file for data err_mut_hash_firstGen \n";
    return;
  }

#endif //histogram
#else // trace run

    ofstream myFile("err_mut_0.1_500");
    ofstream myFileI("err_mut_0.1_1000");
    ofstream myFileII("err_mut_0.1_1500");
    ofstream myFileIII("err_mut_0.1_2000");
    ofstream myFileIV("err_mut_0.1_2500");

#endif //trace runs
#endif //MULTIPLE RUNS

#ifndef TRACE_RUN
#ifdef ACCEPTANCE_RATIO
//de calculat average
  cout << "sample/see " << (unsigned long)(sampleSize_Const/SEE_SAMPLE) << "\n";
#ifndef multiple_restarts
  //  if(nr_runs == 1){
      double myAverage[(unsigned long)(sampleSize_Const/SEE_SAMPLE) + 1];
      double myAcceptanceT[(unsigned long)(sampleSize_Const/SEE_SAMPLE) + 1];
      
      double nrHighIndiv[(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
      unsigned long distHighIndiv[(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];

      double maxIndiv[(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
      for(unsigned long j = 0; j < (unsigned long)(sampleSize_Const/SEE_SAMPLE)+1; j++){
	  maxIndiv[j] = - 2500;
	  distHighIndiv[j] = 0;
      }
      //} else {
#else 
  double myAverage[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
  double myAcceptanceT[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];

  double nrHighIndiv[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
  unsigned long distHighIndiv[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];

  double maxIndiv[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
  for(int i=0; i < nr_runs; i++)
    for(unsigned long j = 0; j < (unsigned long)(sampleSize_Const/SEE_SAMPLE)+1; j++){
	  maxIndiv[i][j] = - 2500;
	  distHighIndiv[i][j] = 0;
    }
  //}
#endif 
  double maximI = -2500;
#endif //acceptance ration

#ifndef MIXTURE_REALS
#ifdef DISTANCE
    double myDistance[nr_runs][sizeGENOME+1];
    double count[sizeGENOME+1], tempDistance[sizeGENOME+1];
    for(unsigned long i = 0; i < sizeGENOME+1; i++){
	for(unsigned long j = 0; j < nr_runs; j++){
	    myDistance[j][i] = 0;
	}
	count[i] = 0;
	tempDistance[i] = 0;
    }
#endif
#else //mixture reals
#ifndef HISTOGRAM
    unsigned long *discreteGenome = new unsigned long[sizeGENOME];
#endif //histogram
#endif //MIXTURE_REALS

#ifdef TRAP_LIKE_FUNCTION
    //store 1111 blocks
    double detailBLOCKS[nr_runs][sizeGENOME/BLOCKsize][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
#endif

#ifdef HISTOGRAM
#ifdef EXPAND_KULLBACK_INFORMATION
    expandTrueUnivariateDistribution();
#else
    readHistogram();
#endif
#endif //histogram
#endif //NOT trance run

/*#ifdef MULTIPLE_RUNS 
   double averageKL[nr_runs];
   #endif*/

  while(runs < nr_runs){

    runs++;
#ifdef multiple_restarts
  for(int j = 0; j < sizeGENOME; j++){
    averageVariablesIndivEach[j] = 0;
    sqrAverageVariablesIndivEach[j] = 0;
  }

  for(int j = 0; j < records; j++){
    averageFitnessIndivRuns[j] = 0;
    averageFitnessRuns[j] = 0;
    for(int k = 0; k < sizeGENOME; k++)
      averageVariablesIndivRuns[j][k] = 0;
    sqrAverageFitnessIndivRuns[j] = 0;
    sqrAverageFitnessRuns[j] = 0;
    for(int k = 0; k < sizeGENOME; k++)
      sqrAverageVariablesIndivRuns[j][k] = 0;
  }
  
  for(int j = 0; j < records; j++)
    for(int k = 0; k < sizeMCMCstring; k++){
      averageFitnessIndiv[j][k] = 0;
      for(int i = 0; i < sizeGENOME; i++)
	averageVariablesIndiv[j][k][i] = 0;
    }

  for(int j = 0; j < records; j++)
    for(int k = 0; k < sizeMCMCstring; k++){
      sqrAverageFitnessIndiv[j][k] = 0;
      for(int i = 0; i < sizeGENOME; i++)
	sqrAverageVariablesIndiv[j][k][i] = 0;
    }
  

  for(int j = 0; j < sizeMCMCstring; j++){
    averageFitnessIndivAll[j] = 0;
    for(int i = 0; i < sizeGENOME; i++)
      averageVariablesIndivAll[j][i] = 0;
  }

  for(int j = 0; j < sizeMCMCstring; j++){
    sqrAverageFitnessIndivAll[j] = 0;
    for(int i = 0; i < sizeGENOME; i++)
      sqrAverageVariablesIndivAll[j][i] = 0;
  }
#endif 

    //setSeed();

#ifndef TRACE_RUN
    double firstBestFitness = 0;
    double secondBestFitness = 0;

    diversityFile << "\n \n runs " << runs << "\n";
#endif //trace run
    double genIndiv = 0;
    
    while(generation <= burn_in + sampleSize*interSamplesize)
      {	
#ifndef multiple_restarts
	if(generation % SEE_SAMPLE = 0){
	  for(int j = 0; j < sizeGENOME; j++){
	    averageVariablesIndivEach[j] = 0;
	    sqrAverageVariablesIndivEach[j] = 0;
	  }
	  
	  for(int j = 0; j < records; j++){
	    averageFitnessIndivRuns[j] = 0;
	    averageFitnessRuns[j] = 0;
	    for(int k = 0; k < sizeGENOME; k++)
	      averageVariablesIndivRuns[j][k] = 0;
	    sqrAverageFitnessIndivRuns[j] = 0;
	    sqrAverageFitnessRuns[j] = 0;
	    for(int k = 0; k < sizeGENOME; k++)
	      sqrAverageVariablesIndivRuns[j][k] = 0;
	  }
	  
	  for(int j = 0; j < records; j++)
	    for(int k = 0; k < sizeMCMCstring; k++){
	      averageFitnessIndiv[j][k] = 0;
	      for(int i = 0; i < sizeGENOME; i++)
		averageVariablesIndiv[j][k][i] = 0;
	    }
	  
	  for(int j = 0; j < records; j++)
	    for(int k = 0; k < sizeMCMCstring; k++){
	      sqrAverageFitnessIndiv[j][k] = 0;
	      for(int i = 0; i < sizeGENOME; i++)
		sqrAverageVariablesIndiv[j][k][i] = 0;
	    }
	   
	  for(int j = 0; j < sizeMCMCstring; j++){
	    averageFitnessIndivAll[j] = 0;
	    for(int i = 0; i < sizeGENOME; i++)
	      averageVariablesIndivAll[j][i] = 0;
	  }
	  
	  for(int j = 0; j < sizeMCMCstring; j++){
	    sqrAverageFitnessIndivAll[j] = 0;
	    for(int i = 0; i < sizeGENOME; i++)
	      sqrAverageVariablesIndivAll[j][i] = 0;
	  }
	}
#endif 

	  if(generation > 0) noCoupling();

#ifdef TRACE_RUN
	  if(generation <= SEE_SAMPLE ){
	    //myFile << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      myFile << generation << "\t";
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFile << tempGENOME[nr_lenght] << "\t";
		//cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFile << "\n";
	      //cout << "\n";
	    }
	    //myFile << "\n";
	  } 
	  if(generation <= 2 * SEE_SAMPLE ){
	    //myFileI << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      myFileI << generation << "\t";
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFileI << tempGENOME[nr_lenght] << "\t";
		//cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFileI << "\n";
	      //cout << "\n";
	    }
	    //myFileI << "\n";
	  }
	  if(generation <= 3 * SEE_SAMPLE ){
	    //myFileII << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      myFileII << generation << "\t";
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFileII << tempGENOME[nr_lenght] << "\t";
		//cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFileII << "\n";
	      //cout << "\n";
	    }
	    //myFileII << "\n";
	  }
	  if(generation <= 4 * SEE_SAMPLE ){
	    //myFileIII << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      myFileIII << generation << "\t";
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFileIII << tempGENOME[nr_lenght] << "\t";
		//cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFileIII << "\n";
	      //cout << "\n";
	    }
	    //myFileIII << "\n";
	  }
	  if(generation <= 5 * SEE_SAMPLE ){
	    //myFileIV << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      myFileIV << generation << "\t";
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFileIV << tempGENOME[nr_lenght] << "\t";
		//cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFileIV << "\n";
	      //cout << "\n";
	    }
	    //myFileIV << "\n";
	  }


#else //trace run

	  postProcessing();

	  averageFitness += averageGeneration();

	  //for MIXTURE_REALS write the samples in a file just like for TRACE_RUN
	  //overwrite the file myFileFit
#ifdef MIXTURE_REALS
	  //write them in one long file
#ifdef BURN_IN_GEN
	  if(generation > throw_gen){
#endif 
	    genIndiv++;

	  jumpGen = 0;
	  for(unsigned long i = 0; i < sizeMCMCstring; i++){

	    MCMCstring[i] -> getGENOME(tempGENOME);
	    double tempValue = MCMCstring[i] -> getValue();

	    if(tempValue != MCMCstring[i] -> fitness()){
	      cout << "true fitness not equal with the read one " << tempValue << " != " << MCMCstring[i] -> fitness() << " in gen " << generation << "\n";
	      exit(1);
	    }
	    
	    if(tempValue < MIN_FITNESS){
	      cout << " smaller then thresold " << tempValue << " != " << tempGENOME[0] << " , " << tempGENOME[1] << " , " << generation << "\n";
	      
	      for(int i1 = 0; i1 < sizeMCMCstring; i1++)
		tempGENOME[i1] = tempPrev[i][i1];
	      tempValue = tempValuePrev[i];
	      MCMCstring[i] -> accepted = 0;
	      //exit(1);
	    }

	    //write the values of the chains
#ifndef HISTOGRAM
	    //discretize before introduce it in one of the hashtable
	    //position in histogram
	    //cout << " dicrete";
	    for(unsigned long i1 =0; i1 < sizeGENOME; i1++ ){
	      discreteGenome[i1] = (unsigned long)((tempGENOME[i1] + scale/2.0)/((double)scale/(double)NrBins));

	      if(discreteGenome[i1] >= NrBins) {
		cout << "WRONG(" << i1 << "," << discreteGenome[i1] << "," << tempGENOME[i1] << ")"; 
		exit(1);
	      }
	    }
	    if(nr_runs ==1){
	      //split in more consecutive and independent regions;
	      myHashtables[(int)((double)generation/ (double)SEE_SAMPLE)]->store_individual(discreteGenome,tempValue);
	    } else {
	      //more runs and for each run a hashtable
	      myHashtables[runs-1]->store_individual(discreteGenome,tempValue);
	    }
#endif //not histogram

	    for(unsigned long j = 0; j < sizeGENOME; j++){
	      //no splited file
	      averageVariablesIndivAll[i][j] += tempGENOME[j];
	      sqrAverageVariablesIndivAll[i][j] += pow(tempGENOME[j],2.0);
	      mse = sqrAverageVariablesIndivAll[i][j]/(generation + 1) - pow(averageVariablesIndivAll[i][j]/(generation + 1),2.0);
	      myFileAllNoSplit << tempGENOME[j] << "\t" << averageVariablesIndivAll[i][j]/(generation + 1) << "\t"<< mse << "\t";

	      averageVariablesIndivEach[j] += tempGENOME[j];
	      sqrAverageVariablesIndivEach[j] += pow(tempGENOME[j],2.0);	      
	      mse = sqrAverageVariablesIndivEach[j]/(generation*sizeMCMCstring + i  + 1) - pow(averageVariablesIndivEach[j]/(generation*sizeMCMCstring + i + 1),2.0);
	      myFile1ChainNoSplit << tempGENOME[j] << "\t" << averageVariablesIndivEach[j]/(generation*sizeMCMCstring + i + 1) << "\t" << mse << "\t";

	      //splited files
	      if(nr_runs == 1){
		averageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j] += tempGENOME[j];
		sqrAverageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j]/ (1 + (int)((double)generation / (double)SEE_SAMPLE)) - 
		  pow(averageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j]/ (1 + (int)((double)generation / (double)SEE_SAMPLE)),2.0);
		myFileAllSplit[(int)((double)generation / (double)SEE_SAMPLE)] << tempGENOME[j] << "\t" << 
		  averageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j]/ (1 + (int)((double)generation / (double)SEE_SAMPLE))<<"\t"<<mse<<"\t";

		averageVariablesIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)][j] += tempGENOME[j];
		sqrAverageVariablesIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)][j]/(1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeMCMCstring + i)
	     - pow(averageVariablesIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)][j]/(1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeMCMCstring + i),2.0);
		myFile1ChainSplit[(int)((double)generation /(double) SEE_SAMPLE)] << tempGENOME[j] << "\t" << 
		  averageVariablesIndivRuns[(int)((double)generation/(double)SEE_SAMPLE)][j]/(1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeMCMCstring + i) 
										  << "\t" << mse <<"\t";
	      } else {
		averageVariablesIndiv[runs-1][i][j] += tempGENOME[j];
		sqrAverageVariablesIndiv[runs-1][i][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndiv[runs-1][i][j]/(generation -throw_gen+ 1) - pow(averageVariablesIndiv[runs-1][i][j]/(generation-throw_gen + 1),2.0);
		if(mse < -0.0001){
		  cout << " mse of average variable runs " << runs << " genome " << i << " generation " << generation << " averageVariable "
		       << averageVariablesIndiv[runs-1][i][j] <<" sqr " << sqrAverageVariablesIndiv[runs-1][i][j] << " mse " << mse << " in negative \n";
		  //exit(1);
		}
		if(mse<0) mse = -mse;
		mse = sqrt(mse);
		myFileAllSplit[runs-1] << tempGENOME[j] << "\t" << averageVariablesIndiv[runs-1][i][j]/(generation-throw_gen + 1) << "\t" << mse << "\t";

		averageVariablesIndivRuns[runs-1][j] += tempGENOME[j];
		sqrAverageVariablesIndivRuns[runs-1][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndivRuns[runs-1][j]/((generation-throw_gen) * sizeMCMCstring + i + 1) 
		  - pow(averageVariablesIndivRuns[runs-1][j]/((generation-throw_gen) * sizeMCMCstring + i + 1),2.0);
		if(mse < -0.0001){
		  cout << " mse of average variable runs all indiv" << runs << " genome " << i << " generation " << generation << " averageVariable "
		       << averageVariablesIndivRuns[runs-1][j] << " sqr " << sqrAverageVariablesIndivRuns[runs-1][j] << " mse " << mse << " in negative \n";
		  //exit(1);
		}
		if(mse<0) mse = -mse;
		mse = sqrt(mse);
		myFile1ChainSplit[runs-1] << tempGENOME[j] << "\t" << averageVariablesIndivRuns[runs-1][j]/ ((generation-throw_gen)*sizeMCMCstring+i+1)<< "\t" << mse << "\t";
	      }
	      //	      cout << tempGENOME[j] << "\t";
	    }

	    
	    //no splited files
	    averageFitnessIndivAll[i] += tempValue;
	    sqrAverageFitnessIndivAll[i] += pow(tempValue,2.0);
	    mse = sqrAverageFitnessIndivAll[i]/(generation-throw_gen + 1) - pow(averageFitnessIndivAll[i]/(generation-throw_gen + 1),2.0);
	    myFileAllNoSplit << tempValue << "\t" << averageFitnessIndivAll[i]/(generation-throw_gen + 1) << "\t" << mse << "\t";

	    averageFitnessIndivEach += tempValue;
	    sqrAverageFitnessIndivEach += pow(tempValue,2.0);
	    mse = sqrAverageFitnessIndivEach/((generation-throw_gen)*sizeMCMCstring + i + 1) - pow(averageFitnessIndivEach/((generation-throw_gen)*sizeMCMCstring+i+1),2.0);
	    myFile1ChainNoSplit << tempValue << "\t" << averageFitnessIndivEach/((generation-throw_gen)*sizeMCMCstring + i + 1) << "\t" << mse << "\n";
	    //cout << tempValue << "\n";

	    //splited files
	    if(nr_runs == 1){
	      averageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i] += tempValue;
	      sqrAverageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i] += pow(tempValue,2.0);
	      mse = sqrAverageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i]/(1 + (int)((double)generation / (double)SEE_SAMPLE)) - 
		pow(averageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i]/ (1 + (int)((double)generation / (double)SEE_SAMPLE)),2.0);
	      myFileAllSplit[(int)((double)generation / (double)SEE_SAMPLE)] << tempValue << "\t" << 
		averageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i]/ (1 + (int)((double)generation / (double)SEE_SAMPLE)) << "\t" << mse << "\t";

	      averageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)] += tempValue;
	      sqrAverageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)] += pow(tempValue,2.0);
	      mse = sqrAverageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)]/(1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeMCMCstring + i) - 
		pow(averageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)]/ (1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeMCMCstring + i),2.0); 
	      myFile1ChainSplit[(int)((double)generation /(double)SEE_SAMPLE)] << tempValue << "\t" << 
		averageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)]/ (1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeMCMCstring + i) << "\t"
									       << mse << "\n";
	    } else {
	      for(int i1 = 0; i1 < sizeGENOME; i1++)
		for(int i2 = i1+1; i2 < sizeGENOME; i2++){
		  int tempIndex = i1 * (sizeGENOME -1) - i1 * (i1-1)/2 + i2 - i1 - 1;
		  corVariablesIndiv[runs-1][i][tempIndex] += tempGENOME[i1]* tempGENOME[i2];
		  sqrCorVariablesIndiv[runs-1][i][tempIndex]+=pow(tempGENOME[i1]*tempGENOME[i2]-tempGENOME[i1]*averageVariablesIndiv[runs-1][i][i2]/(generation-throw_gen+1) - 
		    tempGENOME[i2] * averageVariablesIndiv[runs-1][i][i1]/(generation -throw_gen+ 1) + 
							     averageVariablesIndiv[runs-1][i][i1]*averageVariablesIndiv[runs-1][i][i2]/pow(generation -throw_gen+ 1,2.0),2.0);
		  mse = sqrCorVariablesIndiv[runs-1][i][tempIndex]/(generation-throw_gen +1) - pow(corVariablesIndiv[runs-1][i][tempIndex]/(generation-throw_gen+1) 
					    - averageVariablesIndiv[runs-1][i][i1]*averageVariablesIndiv[runs-1][i][i2]/pow(generation-throw_gen + 1,2.0),2.0);
		  cor = corVariablesIndiv[runs-1][i][tempIndex]/(generation-throw_gen+1) - 
		    averageVariablesIndiv[runs-1][i][i1]*averageVariablesIndiv[runs-1][i][i2]/pow(generation-throw_gen + 1,2.0);
		  s1 = sqrAverageVariablesIndiv[runs-1][i][i1]/(generation-throw_gen + 1) - pow(averageVariablesIndiv[runs-1][i][i1]/(generation-throw_gen + 1),2.0);
		  s2 = sqrAverageVariablesIndiv[runs-1][i][i2]/(generation-throw_gen + 1) - pow(averageVariablesIndiv[runs-1][i][i2]/(generation-throw_gen + 1),2.0);
		  if(mse < -0.0001){
		    cout << " mse of average correlation runs " << runs << " genome " << i << " generation " << generation << " averageCor "
			 << cor << " sqr " << sqrCorVariablesIndiv[runs-1][i][tempIndex] << " mse " << mse << " in negative \n";
		    //exit(1);
		  }
		  if(mse<0) mse = -mse;
		  mse = sqrt(mse);
		  if(s1 > 0 && s2 > 0)
		    myFileAllSplit[runs-1] << cor << "\t" << mse << "\t" << cor/sqrt(s1 * s2) << "\t";
		  else 
		    if(s1 < -0.0001 || s2 < -0.0001){
		      cout << " s1 " << s1 << " s2 " << s2 << "\n";
		      exit(1);
		    }
		    myFileAllSplit[runs-1] << cor << "\t" << mse << "\t" << "0" << "\t";		    
		}

	      jump1 = 0;
	      if(generation != 0){
		for(int i1 = 0; i1 < sizeGENOME; i1++)
		  jump1 += pow(tempGENOME[i1] - tempPrev[i][i1],2.0);
		jump1 = jump1;
		
		//averge jumping for each chains
		averageJumpIndiv[runs-1][i] += sqrt(jump1/sizeGENOME/scale);
		jumpGen += sqrt(jump1/sizeGENOME/scale);
		sqrAverageJumpIndiv[runs-1][i] += jump1/sizeGENOME/scale;
		mse = sqrAverageJumpIndiv[runs-1][i]/(generation-throw_gen + 1) - pow(averageJumpIndiv[runs-1][i]/(generation-throw_gen + 1),2.0); 
		if(mse < -0.0001){
		  cout << " mse of average jumps runs " << runs << " genome " << i << " generation " << generation << " averageJump "
		       << averageJumpIndiv[runs-1][i] << " sqr " << sqrAverageJumpIndiv[runs-1][i] << " mse " << mse << " in negative \n";
		  //exit(1);
		}
		// mse = sqrt(abs(mse));
		if(mse<0) mse = -mse;
		mse = sqrt(mse);
	      }
	      if(generation != 0)
		myFileAllSplit[runs-1] << averageJumpIndiv[runs-1][i]/(generation-throw_gen + 1) << "\t" << mse << "\t";
	      else 
		myFileAllSplit[runs-1] << "0\t0\t";		
	      
	      averageFitnessIndiv[runs-1][i] += tempValue;
	      sqrAverageFitnessIndiv[runs-1][i] += pow(tempValue,2.0);
	      mse = sqrAverageFitnessIndiv[runs-1][i]/(generation-throw_gen + 1) - pow(averageFitnessIndiv[runs-1][i]/(generation-throw_gen + 1),2.0);
	      myFileAllSplit[runs-1] << tempValue << "\t" << averageFitnessIndiv[runs-1][i]/(generation-throw_gen + 1) << "\t" << mse << "\t"; 

	      for(int i1 = 0; i1 < sizeGENOME; i1++)
		for(int i2 = i1+1; i2 < sizeGENOME; i2++){
		  int tempIndex = i1 * (sizeGENOME -1) - i1 * (i1-1)/2 + i2 - i1 - 1;
		  corVariablesIndivRuns[runs-1][tempIndex] += tempGENOME[i1]* tempGENOME[i2];
		  double ave1 = averageVariablesIndivRuns[runs-1][i1]/((generation- throw_gen) * sizeMCMCstring + i+ 1);
		  double ave2 = averageVariablesIndivRuns[runs-1][i2]/((generation- throw_gen) * sizeMCMCstring + i+ 1);
		  sqrCorVariablesIndivRuns[runs-1][tempIndex] += pow((tempGENOME[i1] - ave1) * (tempGENOME[i2] - ave2),2.0);
		  cor = corVariablesIndivRuns[runs-1][tempIndex]/((generation - throw_gen)* sizeMCMCstring + i + 1) - ave1 * ave2;

		  mse = sqrCorVariablesIndivRuns[runs-1][tempIndex]/((generation - throw_gen)* sizeMCMCstring + i +1) - pow(cor,2.0);

		  s1 = sqrAverageVariablesIndivRuns[runs-1][i1]/((generation- throw_gen)*sizeMCMCstring+i+1) - pow(ave1,2.0);
		  s2 = sqrAverageVariablesIndivRuns[runs-1][i2]/((generation- throw_gen)*sizeMCMCstring+i+1) - pow(ave2,2.0);

		  if(mse < -0.0001){
		    cout << " mse of average jumps runs " << runs << " genome " << i << " generation " << generation << " averageJump "
			 << cor << " sqr " << sqrCorVariablesIndivRuns[runs-1][tempIndex]/((generation - throw_gen)* sizeMCMCstring + i +1) << " mse " << mse 
			 << " in negative \n";
		    //exit(1);
		  }
		  // mse = sqrt(abs(mse));
		  if(mse<0) mse = -mse;
		  mse = sqrt(mse);

		  if(s1 > 0 && s2 > 0)
		    myFile1ChainSplit[runs-1] << cor << "\t" << mse << "\t" << cor/sqrt(s1 * s2) << "\t";
		  else 
		    myFile1ChainSplit[runs-1] << cor << "\t" << mse << "\t" << 0 << "\t";
		}

	      // if(generation != 0){
		averageJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		sqrAverageJumpIndivRuns[runs-1] += pow(jump1/sizeGENOME/scale,2.0);
		mse = sqrAverageJumpIndivRuns[runs-1]/((generation-throw_gen)*sizeMCMCstring+i+1)
		  -pow(averageJumpIndivRuns[runs-1]/((generation-throw_gen)*sizeMCMCstring+i+1),2.0); 
		//}

	      if(generation != 0)
		myFile1ChainSplit[runs-1] << averageJumpIndivRuns[runs-1]/((generation- throw_gen)*sizeMCMCstring+i+1) << "\t" << mse << "\t";
	      else 
		myFile1ChainSplit[runs-1] << "0\t0\t";		

	      if(generation != 0){
		if(tempValuePrev[i] < tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		}
		else{ 
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation-throw_gen)*sizeMCMCstring+i+1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation-throw_gen)*sizeMCMCstring+i+1),2.0); 
	      }

	      if(generation != 0)
		myFile1ChainSplit[runs-1] << averageExpJumpIndivRuns[runs-1]/((generation- throw_gen)*sizeMCMCstring+i+1) << "\t" << mse << "\t";
	      else 
		myFile1ChainSplit[runs-1] << "0\t0\t";		

	      for(int i1 = 0; i1 < sizeGENOME; i1++)
		tempPrev[i][i1] = tempGENOME[i1];
	      tempValuePrev[i] = tempValue;

	      averageFitnessIndivRuns[runs-1] += tempValue;
	      sqrAverageFitnessIndivRuns[runs-1] += pow(tempValue,2.0);
	      mse = sqrAverageFitnessIndivRuns[runs-1]/ ((generation- throw_gen)*sizeMCMCstring + i + 1) 
		- pow(averageFitnessIndivRuns[runs-1]/((generation- throw_gen) * sizeMCMCstring + i + 1),2.0);
	      myFile1ChainSplit[runs-1] << tempValue << "\t" << averageFitnessIndivRuns[runs-1]/((generation- throw_gen)*sizeMCMCstring + i + 1) << "\t" << mse << "\n"; 
	    }
	      //myString.insert(i*(2*sizeGENOME + 2) + 2*(sizeGENOME - 1) + 1,);

	  }

	  //myFileFit << "\n";
	  myFileAllNoSplit << averageGeneration()  << "\t" << averageFitness/(generation+1)<< "\n";

	  if(nr_runs == 1){
	    averageFitnessRuns[(int)((double)generation / (double)SEE_SAMPLE)] += averageGeneration();
	    myFileAllSplit[(int)((double)generation / (double)SEE_SAMPLE)] << averageGeneration() << "\t"
		    << averageFitnessRuns[(int)((double)generation / (double)SEE_SAMPLE)]/(1 + (int)((double)generation / (double)SEE_SAMPLE)) << "\n";
	  } else {
	    averageFitnessRuns[runs-1] += averageGeneration();
	    myFileAllSplit[runs-1] << averageGeneration()  << "\t" << averageFitnessRuns[runs]/(generation + 1) << "\t"; 

	  //compute the jumps for all chains
	    averageJump[runs-1] += jumpGen/sizeMCMCstring;
	    sqrAverageJump[runs-1] += pow(jumpGen/sizeMCMCstring,2.0);
	    mse = sqrAverageJump[runs-1]/(generation -throw_gen+ 1) - pow(averageJump[runs-1]/(generation-throw_gen + 1),2.0);
	    myFileAllSplit[runs-1] << averageJump[runs-1]/(generation -throw_gen+ 1) << "\t" << mse << "\n";   
	  }


	  //each generation is separated by a witte line 	
	  //how to introduce individuals in hashtable when it is not a diversity measure ?!?
#ifdef BURN_IN_GEN
	  } else if(generation == throw_gen) {
	    for(unsigned long i =0 ; i < sizeMCMCstring; i++){
	      MCMCstring[i] -> getGENOME(tempGENOME);
	      double tempValue = MCMCstring[i] -> fitness();
	      /*   //cout << "transform dicrete" << i << "(" << tempGENOME[0] << "," << tempGENOME[1] << ")\n";
	      for(unsigned long jack =0; jack < sizeGENOME; jack++ ){
		discreteGenome[jack] = (unsigned long)((tempGENOME[jack] + scale/2.0)/((double)scale/(double)NrBins));
		//cout << "(" << i << " , " << tempGENOME[i] << " , " << discreteGenome[i] << ")\t"; 
		if(discreteGenome[i] >= NrBins) cout << "wrong individual ? (" << i << "," << discreteGenome[i] << "," << tempGENOME[i] << ")"; 
	      }
	      //cout << "\n";
	      if(nr_runs ==1){
		//split in more consecutive and independent regions;
		myHashtables[(int)((double)generation/ (double)SEE_SAMPLE)]->store_individual(discreteGenome,tempValue);
	      } else {
		//more runs and for each run a hashtable
		myHashtables[runs-1]->store_individual(discreteGenome,tempValue);
		//if(generation % sampleSize_Const == sampleSize_Const - 1)
		//	myFileFit << "\n";	      
		}*/
	      for(unsigned long i =0; i < sizeMCMCstring; i++){
		for(int i1 = 0; i1 < sizeGENOME; i1++)
		  tempPrev[i][i1] = tempGENOME[i1];
		tempValuePrev[i] = tempValue;
	      }
	    }
	  }
		  
#endif
#else // is not mixture reals

#ifndef MULTIPLE_RUNS
	if(runs == 1){
	  for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    MCMCstring[i] -> getGENOME(tempGENOME);
	    myFileFit << MCMCstring[i] -> getValue() << "\t";
	  }
	  myFileFit << "\n";
	} else if(runs == 2){
	  for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    MCMCstring[i] -> getGENOME(tempGENOME);
	    myFileFG << MCMCstring[i] -> getValue() << "\t";
	  }
	  myFileFG << "\n";
	}
	myTempFile << MCMCstring[0]->getValue() << "\n" << averageGeneration() << "\n";
#endif	

#endif //mixture reals	  

#ifndef FILE_BOOKEEPING
	if(generation % SEE_SAMPLE == 0){
	    printDiversityPopulation();
	}
#endif

#ifndef MIXTURE_REALS
#ifdef DISTANCE
	//int tempGenome[sizeGENOME];
	//MCMCstring[0]->getMax(tempGenome);
	double averageDistance = 0;
	for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    unsigned long ate = (unsigned long)HammingDistance(i);
	    myDistance[runs-1][ate] += MCMCstring[i]->fitness();
	    count[ate] ++; 
	    averageDistance += ate;
#ifndef MULTIPLE_RUNS
	if(i == 0)
	    myTempFile << ate << "\n";
#endif	
	}
	averageDistance/=sizeMCMCstring;
#ifndef MULTIPLE_RUNS
	//if(runs == 1)
	myTempFile << averageDistance << "\n";
#endif //multiple runs
#endif //DISTANCE
#else //mixture reals

#ifdef HISTOGRAM
	double averageDistance = 0;
	for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    double ate = (double)HammingDistance(i);
	    averageDistance += ate;
#ifndef MULTIPLE_RUNS
	if(i == 0)
	    myTempFile << ate << "\n";
#endif	
	}
	averageDistance/=sizeMCMCstring;
#ifndef MULTIPLE_RUNS
	//if(runs == 1)
	myTempFile << averageDistance << "\n";
#endif
#endif //histogram
#endif //mixture reals

#ifdef HISTOGRAM
#ifdef KULLBACK_INFORMATION
	      //cout << "KL diff" << mixingTimeUnivariateCalculation()<< "\n";
	      if(generation % SEE_SAMPLE == 0){	
		  double temp = mixingTimeUnivariateCalculation();
		  //Kullback[generation] += temp;
		  simul = insert_simulation(simul,generation,temp);
/*#ifdef MULTIPLE_RUNS
		  if(generation == sampleSize)
		      averageKL[runs-1] = temp;
#endif*/
		  show_simulations(simul);
		  write_to_file_all(simul_all_proc_right,simul,NULL);
		  //mixingTimeUnivariateCalculation();
	      }      
#endif //KULLBACK_INFORMATION
#endif //histogram

#ifdef ACCEPTANCE_RATIO
	for(int i = 0; i < sizeMCMCstring; i++)
	    if(maximI < MCMCstring[i]->getValue()){
		maximI = MCMCstring[i]->getValue();
	    }
	//if(nr_runs == 1){
#ifndef multiple_restarts
	  acceptance_ratio[(int)((double)generation / (double)SEE_SAMPLE)] += acceptanceAverage();
	  if(generation % SEE_SAMPLE == 0){
	    //printAcceptance(myAcceptance);
	    myAverage[(unsigned long)(generation/SEE_SAMPLE)] = averageGeneration();	
	    myAcceptanceT[(unsigned long)(generation/SEE_SAMPLE)] = acceptanceAverage();	
	    nrHighIndiv[(unsigned long)(generation/SEE_SAMPLE)] = myHashtables[runs-1]->nrElem;
	    //	    distHighIndiv[runs-1][(int)(generation/SEE_SAMPLE)] = nrHighIndivGlobal;
	    //distHighIndiv[runs-1][(int)(generation/SEE_SAMPLE)] = averageDistanceGoodIndiv();
	    maxIndiv[(int)(generation/SEE_SAMPLE)] = maximI;
	  }
	  //} else {
#else
	  acceptance_ratio[runs-1] += acceptanceAverage();
	  averageFitness += averageGeneration();
	  if(generation % SEE_SAMPLE == 0){
	    //printAcceptance(myAcceptance);
	    myAverage[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = averageGeneration();	
	    myAcceptanceT[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = acceptanceAverage();	
	    
	    //if(nr_runs == 1)
	    nrHighIndiv[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = myHashtables[runs-1]->nrElem;	
	    //distHighIndiv[runs-1][(int)(generation/SEE_SAMPLE)] += nrHighIndivGlobal;
	      //averageDistanceGoodIndiv();

	    maxIndiv[runs-1][(int)(generation/SEE_SAMPLE)] = maximI;
	  } 
	  //}
#endif //nr_runs
#endif

#ifdef TRAP_LIKE_FUNCTION
 	if(generation % SEE_SAMPLE == 0){
	  
	  for(unsigned long i = 0 ; i < sizeGENOME/BLOCKsize; i++)
#ifdef BINOMIAL
	    detailBLOCKS[runs-1][i][(unsigned long)(sampleSize_Const/SEE_SAMPLE)] = diversityDetail[i][BLOCKsize];
#else
	    detailBLOCKS[runs-1][i][(unsigned long)(sampleSize_Const/SEE_SAMPLE)] = diversityDetail[i][a];
#endif
	}
#endif //TRAP_LIKE_FUNCTION

//print the outout for further analyses
#ifdef PLOTS_MCMC
	//printAllChains_I(myFile);
#endif
	      //functii pentru hitting time
#ifdef KNOWN_MAXIM
#ifndef TWO_ATTRACTORS_FUNCTION
	for(unsigned long i=0; i < sizeMCMCstring; i++){
	    if(MCMCstring[i]->stop() == 1 && firstBestFitness == 0){
		bestFitness += generation;
		firstBestFitness = 1;
		bestFitnessVector[runs-1] = generation; 
	    } else if(MCMCstring[i]->stop() == 2 && secondBestFitness == 0){
		//bestFitness += generation;
		secondBestFitness = 1;
		secondBestFitnessVector[runs-1] = generation; 
	    } 
	}
#else
//#ifdef HISTOGRAM
	for(unsigned long i=0; i < sizeMCMCstring; i++){
	    if(MCMCstring[i]->stop() == 1 && firstBestFitness == 0){
		bestFitness += generation;
		firstBestFitness = 1;
		bestFitnessVector[runs-1][0] = generation; 
	    } else if(MCMCstring[i]->stop() == 2 && secondBestFitness == 0){
		//bestFitness += generation;
		secondBestFitness = 1;
		bestFitnessVector[runs-1][1] = generation; 
	    }	
	}

#endif
#endif //known matrix

#endif //NOT trace run
	generation++;	
      }

#ifndef TRACE_RUN
#ifdef HISTOGRAM 
    cout << "Begin convergence diagnostic \n";
#ifdef PLOTS_MCMC
/*    if(runs == 1){
    myFile << "run=" << nr_runs<<"\n";
    printAllChains(myFile);
    }*/
    // printHashtable(myFile);
#else //plot mcmc
#ifdef DEBUG_ 
    printLEVEL(myFile);
#endif //debug
#endif //plot mcmc

    myFileHash << "max " << MCMCstring[0]->getMaxValue() << "\n";
    printHashtable(myFileHash);

#ifndef MIXTURE_REALS
    scoreDiversity();
#endif 

#ifdef TRAP_LIKE_FUNCTION
    scoreDiversityDetail();
#endif
    
    globalRightSolutions[runs-1] = procentRightSolutions;
        
#ifdef GELMAN_RUBIN
    varianceHigerMomentsUnivariateCalculation(GELMAN_Moment,sizeMCMCstring - 1);
    varianceMultivariateCalculation();
#endif
    cout << "End convergence diagnostic \n";

#ifndef MIXTURE_REALS
#ifdef DISTANCE
	for(unsigned long i = 0; i < sizeGENOME+1; i++){
	    if(count[i] != 0){
		myDistance[runs-1][i] /= count[i];
		count[i] = 0;
	    }
	}
#endif
#endif //MIXTURE_REALS
#endif //histogram

#ifdef CUSUM
    CusumCalculation();
#endif

#endif //trace run

    cout << " indiv mut " << genIndiv << "\n";
	
	if(runs != nr_runs) reset();
  }
  
#ifndef TRACE_RUN

  /// write in the file the hashtables -> in the one that supose to write the values of the second run
  //cout << "arive here 2 \n";
#ifdef MIXTURE_REALS
#ifndef HISTOGRAM
  //write the general statistical data
  //cout << "arive here 3 \n";
  
  double populationVariance, populationMean, centerMass;
  cout << " runs " << nr_runs << " sampleSize " << sampleSize_Const << " sample " << SEE_SAMPLE << " \n";
    if(nr_runs ==1){
    //split in more consecutive and independent regions;
    myHashtables[0]->printStd(myFileHash,myHashtables,(int)((double)sampleSize_Const/ (double)SEE_SAMPLE) + 1,centerMass,populationVariance);
  } else {
    //more runs and for each run a hashtable
    myHashtables[0]->printStd(myFileHash,myHashtables,nr_runs,centerMass,populationVariance);
    }  
  //write the data in a file for histogram
  ofstream myHistogram(fileHistogram,ios::app);		
  if(!myHistogram.is_open()){
      cout << "Error opening histogram file \n";
      return;
  }
  
  myHistogram << "# no coupling \n";
  myHistogram << sizeGENOME << "\t" << sampleSize_Const << "\t" << nr_runs << "\t" << (int)((double)sampleSize_Const/ (double)SEE_SAMPLE) + 1<< "\t"
	      << centerMass << "\t" << populationMean << "\t" << populationVariance << "\n";

  myHistogram.flush();
  myHistogram.close();
#endif //histogram
#endif //mixutre reals

#ifdef HISTOGRAM  
#ifdef TRAP_LIKE_FUNCTION
  printDiversityDetailGNUFILE();
#endif
  printDiversityRightSolutions();

#ifdef GELMAN_RUBIN
  printReductionGNUFILE(GELMAN_Moment);
  printReductionMultivariateGNUFILE();
#endif

#ifdef TWO_ATTRACTORS_FUNCTION
  ofstream myHistogram(fileHistogram);		
  if(!myHistogram.is_open()){
      cout << "Error opening histogram file \n";
      return;
  }

  for(unsigned long k = 0; k < sizeGENOME+1; k++){
      for(unsigned long j = 0; j < nr_runs; j++)
	  globalHistogram[k] += ActualHistogram[k][j];
      globalHistogram[k] /= nr_runs;
      myHistogram << globalHistogram[k];
  }
  
  myHistogram.flush();
  myHistogram.close();
#endif
  
#ifdef DEBUG_
  myFile.close();
#endif

#ifndef FILE_BOOKEEPING
  diversityFile.flush();
  diversityFile.close();
#endif

#endif //histogram

/*#ifdef KULLBACK_INFORMATION
  printMixingGNUFILE();
  #endif*/

#ifdef CUSUM
  double hair = printCusumGNUFILE();
#endif


#ifdef HISTOGRAM
  printDiversityGNUFILE();
#endif //histogram

#ifdef ACCEPTANCE_RATIO
  ofstream myAcceptance(fileAcceptance,ios::out|ios::app);
  if(!myAcceptance.is_open()){
      cout << "Error opening file acceptance with no coupling \n";
      return;
  } 
  double cumA = 0, cumV = 0, cumN = 0;
  double meanAv, std_meanAv;

#ifdef multiple_restarts 
  for(unsigned long i = 0; i < (unsigned long)(sampleSize_Const/SEE_SAMPLE)+1; i++){
      double meanA = 0, std_devA = 0;
      double meanV = 0, std_devV = 0;
      double meanN = 0, std_devN = 0;
      double cusum_mean = 0, hairness_mean = 0;
      double std_dev_c = 0, std_dev_h = 0;
      double meanD = 0, std_devD = 0;
      double meanM = 0, std_devM = 0;

      for(unsigned long j = 0; j < nr_runs; j++){
	  meanA += myAcceptanceT[j][i];
	  meanV += myAverage[j][i];
	  meanN += nrHighIndiv[j][i];
	  meanD += distHighIndiv[j][i];
	  if(i !=sampleSize_Const/SEE_SAMPLE){
	    cusum_mean += Cusum[j][i*SEE_SAMPLE];
	    hairness_mean += Hairness[j][i*SEE_SAMPLE];
	  } else {
	    cusum_mean += Cusum[j][i*SEE_SAMPLE-3];
	    hairness_mean += Hairness[j][i*SEE_SAMPLE-3];
	  }
	  meanM += maxIndiv[j][i];
      }      

      meanA /= nr_runs;
      meanV /= nr_runs;
      meanN /= nr_runs;
      cusum_mean /= nr_runs;
      hairness_mean /= nr_runs;
      meanD /= nr_runs;
      meanM /= nr_runs;

      for(unsigned long j = 0; j < nr_runs; j++){
	  std_devA += pow(myAcceptanceT[j][i] - meanA,2);
	  std_devV += pow(myAverage[j][i] - meanV,2);
	  std_devN += pow(nrHighIndiv[j][i] - meanN,2);
	  std_devD += pow(distHighIndiv[j][i] - meanD,2);
	  if(i != sampleSize_Const/SEE_SAMPLE){
	    std_dev_c += pow(cusum_mean - Cusum[j][i*SEE_SAMPLE],2);
	    std_dev_h += pow(hairness_mean - Hairness[j][i*SEE_SAMPLE],2);
	  } else {
	    std_dev_c += pow(cusum_mean - Cusum[j][i*SEE_SAMPLE-3],2);
	    std_dev_h += pow(hairness_mean - Hairness[j][i*SEE_SAMPLE-3],2);
	  }
	  std_devM += pow(maxIndiv[j][i] - meanM,2);
      }

      std_devA = sqrt(std_devA/nr_runs);
      std_devV = sqrt(std_devV/nr_runs);
      std_devN = sqrt(std_devN/nr_runs);
      std_devD = sqrt(std_devD/nr_runs);
      std_dev_c = sqrt(std_dev_c/nr_runs);
      std_dev_h = sqrt(std_dev_h/nr_runs);
      std_devM = sqrt(std_devM/nr_runs);

      cumA += meanA;
      cumV += meanV;
      cumN += meanN;

      myAcceptance << i*SEE_SAMPLE << "\t" << meanA << "\t" << std_devA << "\t" << cumA << "\t"<< meanV << "\t" << std_devV << "\t" << cumV  		    
		   << "\t" << meanN << "\t" << std_devN << "\t" << cumN << "\t" << meanM << "\t" << std_devM << "\t" << meanD << "\t" << std_devD 
		   << "\t" << cusum_mean << "\t" << std_dev_c  << "\t" << hairness_mean << "\t" << std_dev_h 
		   << "\n";

      //meanAv += meanA;
//	std_meanAv += pow(meanA,2.0)/nr_runs;

      }
      for(unsigned long j = 0; j < nr_runs; j++){
	  meanAv += acceptance_ratio[j]/generation;
	  std_meanAv += pow(acceptance_ratio[j]/generation,2.0);
      }
      
      meanAv /= nr_runs;
      std_meanAv = std_meanAv/nr_runs - pow(meanAv,2.0);
      myAcceptance << meanAv << "\t" << std_meanAv << "\n";
#else 
	double meanA = 0, std_devA = 0;
	double meanV = 0, std_devV = 0;
	double meanN = 0, std_devN = 0;
	double cusum_mean = 0, hairness_mean = 0;
	double std_dev_c = 0, std_dev_h = 0;
	double meanD = 0, std_devD = 0;
	double meanM = 0, std_devM = 0;

       for(unsigned long i = 0; i < (sampleSize_Const/SEE_SAMPLE) + 1; i++){
	  meanA += myAcceptanceT[i];
	  meanV += myAverage[i];
	  meanN += nrHighIndiv[i];
	  meanD += distHighIndiv[i];
	  if(i !=sampleSize_Const/SEE_SAMPLE){
	    cusum_mean += Cusum[i*SEE_SAMPLE];
	    hairness_mean += Hairness[i*SEE_SAMPLE];
	  } else {
	    cusum_mean += Cusum[i*SEE_SAMPLE-3];
	    hairness_mean += Hairness[i*SEE_SAMPLE-3];
	  }
 	  meanM += maxIndiv[i];
       }

	cusum_mean /= (sampleSize_Const/SEE_SAMPLE) + 1;
	hairness_mean /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanA /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanV /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanN /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanD /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanM /= (sampleSize_Const/SEE_SAMPLE) + 1;
	
	for(unsigned long j = 0; j < (sampleSize_Const/SEE_SAMPLE) + 1; j++){
	  std_devA += pow(myAcceptanceT[j] - meanA,2);
	  std_devV += pow(myAverage[j] - meanV,2);
	  std_devN += pow(nrHighIndiv[j] - meanN,2);
	  std_devD += pow(distHighIndiv[j] - meanD,2);
	  if(j != sampleSize_Const/SEE_SAMPLE){
	    std_dev_c += pow(cusum_mean - Cusum[j*SEE_SAMPLE],2);
	    std_dev_h += pow(hairness_mean - Hairness[j*SEE_SAMPLE],2);
	  } else {
	    std_dev_c += pow(cusum_mean - Cusum[j*SEE_SAMPLE-3],2);
	    std_dev_h += pow(hairness_mean - Hairness[j*SEE_SAMPLE-3],2);
	  }
	  std_devM += pow(maxIndiv[j] - meanM,2);      
	}      
	
	std_devA = sqrt(std_devA/(sampleSize_Const/SEE_SAMPLE + 1));
	std_devV = sqrt(std_devV/(sampleSize_Const/SEE_SAMPLE + 1));
	std_devN = sqrt(std_devN/(sampleSize_Const/SEE_SAMPLE + 1));
	std_devD = sqrt(std_devD/(sampleSize_Const/SEE_SAMPLE + 1));
	std_dev_c = sqrt(std_dev_c/(sampleSize_Const/SEE_SAMPLE + 1));
	std_dev_h = sqrt(std_dev_h/(sampleSize_Const/SEE_SAMPLE + 1));
	std_devM = sqrt(std_devM/(sampleSize_Const/SEE_SAMPLE + 1));
	
	cumA += meanA;
	cumV += meanV;
	cumN += meanN;
	
	myAcceptance << sampleSize_Const*SEE_SAMPLE << "\t" << meanA << "\t" << std_devA << "\t" << meanA * ((sampleSize_Const/SEE_SAMPLE) + 1) << "\t" 
		     << meanV << "\t" << std_devV << "\t" << meanV * ((sampleSize_Const/SEE_SAMPLE) + 1)   		    
		     << "\t" << meanN << "\t" << std_devN << "\t" << meanN*((sampleSize_Const/SEE_SAMPLE) + 1) << "\t" 
		     << meanM << "\t" << std_devM << "\t" << meanD << "\t" << std_devD 
		     << "\t" << cusum_mean << "\t" << std_dev_c  << "\t" << hairness_mean << "\t" << std_dev_h 
		     << "\n";
#endif //nr_runs  

  myAcceptance.flush();
  myAcceptance.close();
#endif

#ifndef MIXTURE_REALS
#ifdef DISTANCE
    ofstream myFileDistance("distance_no_coupling.txt");
  cout << "average on this distance \n";
  double std_dev[sizeGENOME+1];
  for(unsigned long i = 0; i < sizeGENOME+1; i++){
      tempDistance[i] = 0;
      std_dev[i] = 0;

      for(unsigned long j = 0; j < nr_runs; j++)
	  tempDistance[i] += myDistance[j][i];

      tempDistance[i] /= nr_runs;

      for(unsigned long j = 0; j < nr_runs; j++)
	  std_dev[i] += pow(myDistance[j][i] - tempDistance[i],2);

      std_dev[i] = sqrt(std_dev[i]/nr_runs);

      myFileDistance << i << "\t" << tempDistance[i] << "\t" <<  std_dev[i] << "\n";
  }
  
  myFileDistance.flush();
  myFileDistance.close();
#endif
#else //mixture reals
#ifdef HISTOGRAM
  //delete[] discreteGenome;
#endif //histogram
#endif //MIXTURE_REALS

#ifndef MULTIPLE_RUNS
#ifdef HISTOGRAM
  myTempFile.flush();
  myTempFile.close();

  myFileFG.flush();
  myFileFG.close();

  myFileFit.flush();
  myFileFit.close();

  myFileHash.flush();
  myFileHash.close();


  ofstream myFileD("err_mut_dist");		
  if(!myFileD.is_open()){
      cout << "Error opening file for data " << fileData << "\n";
    return;
  }

  ifstream myTempFile1("temporar.txt");		
  if(!myTempFile1.is_open()){
      cout << "Error opening file for data " << fileData << "\n";
    return;
  }

  double tempV[sampleSize+2][nr_runs+1];
  double tempD[sampleSize+2][nr_runs+1];
  double tempV1[sampleSize+2][nr_runs+1];
  double tempD1[sampleSize+2][nr_runs+1];

  for(unsigned long i = 0; i < nr_runs; i++)
      for(unsigned long j = 0; j < sampleSize+1;j++){
	  myTempFile1 >> tempV[j][i];
	  myTempFile1 >> tempV1[j][i];
	  myTempFile1 >> tempD[j][i];
	  myTempFile1 >> tempD1[j][i];
      }

  for(unsigned long j = 0; j < sampleSize+1;j++){
      for(unsigned long i = 0; i < nr_runs; i++){
	  myFile << tempV[j][i] << "\t";
	  myFile << tempV1[j][i] << "\t";
	  myFileD << tempD[j][i] << "\t";
	  myFileD << tempD1[j][i] << "\t";
      }
      myFile << "\n";
      myFileD << "\n";
  }
  myTempFile1.close();

//#ifdef DEBUG_
  myFile.flush();
  myFile.close();
//#endif
  myFileD.flush();
  myFileD.close();

#else //histogram
  for(int i = 0; i < records; i++){
    myFileAllSplit[i].flush();
    myFileAllSplit[i].close();

    myFile1ChainSplit[i].flush();
    myFile1ChainSplit[i].close();
  }

  myFileAllNoSplit.flush();
  myFileAllNoSplit.close();

  myFile1ChainNoSplit.flush();
  myFile1ChainNoSplit.close();

  delete[] myFileAllSplit;
  delete[] myFile1ChainSplit;

  myFileHash.flush();
  myFileHash.close();
#endif //histogram
#endif //MULTIPLE_RUNS

#else

  myFile.flush();
  myFile.close();

  myFileI.flush();
  myFileI.close();

  myFileII.flush();
  myFileII.close();

  myFileIII.flush();
  myFileIII.close();

  myFileIV.flush();
  myFileIV.close();

#endif //trace run
}
#else
////////////////////
//SINGLE_vs_MULTIPLE
////////////////////
void binaryPopulation::evolution(){
    cout << "NO COUPLING \n";
    
    long beginMCMCstring = sizeMCMCstring;
    
    ofstream myFile(fileData);		
    if(!myFile.is_open()){
	cout << "Error opening file for data " << myFile << "\n";
	return;
    }
    
    ifstream mySample(sampleFile);
    if(!mySample.is_open()){
	cout <<"Error in opening the distribution file\n";
	exit(1);
    }
    
    ofstream myData(gnuplotFILE);
    if(!myData.is_open()){
	cout <<"Error in opening the diversity file\n";
	exit(1);
    }

    ofstream myHitting("hitting_time.txt",ios::out|ios::app);
    if(!myHitting.is_open()){
	cout << "Error in opening the hitting time file \n";
	exit(1);
    }

#ifdef CUSUM
      ofstream myCUSUM(CusumFILE);
      if(!myCUSUM.is_open()){
      cout <<"Error in opening the mixing file "<< CusumFILE <<"\n";
      exit(1);
      }
#endif    

    ofstream myMixing(MixingFILE, ios::out | ios:: app);
    if(!myMixing.is_open()){
	cout <<"Error in opening the mixing file "<< MixingFILE <<"\n";
	exit(1);
    }

#ifdef ACCEPTANCE_RATIO
    ofstream myAcceptance(fileAcceptance, ios::out | ios::app);
    if(!myAcceptance.is_open()){
	cout << "Error in opening the acceptance file \n";
	exit(1);
    }
#endif

    //unsigned long dataFromSampleFile;
    double performance[(unsigned long)(log(sizeMCMCstring)/log(2.0))+1][diversityMeasure];
    for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++)
	for(unsigned long j = 0; j < diversityMeasure; j++)
	    performance[i][j] = 0;
    
    double performance_runs[(unsigned long)(log(sizeMCMCstring)/log(2.0))+1][diversityMeasure];
    for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++)
	for(unsigned long j = 0; j < diversityMeasure; j++)
	    performance_runs[i][j] = 0;
    
#ifdef EXPAND_KULLBACK_INFORMATION
    expandTrueUnivariateDistribution();
#else
    readHistogram();
#endif
    
    for(unsigned long i =0; i< sizeMCMCstring; i++)
	MCMCstring[i]->temperature = bestTemperature;

    while(runs < nr_runs){
	
	runs++;
	setSeed();
	
	//assignTemperature();     
	sizeMCMCstring = beginMCMCstring;
	sampleSize = sampleSize_Const;
	
	generation = 0;
	bestFitness = 0;
	double firstBestFitness[(unsigned long)(log(sizeMCMCstring)/log(2))+1];
	for(unsigned long i =0; i<=(unsigned long)(log(sizeMCMCstring)/log(2)); i++)
	    firstBestFitness[i] = 0;
	
#ifdef TWO_ATTRACTORS_FUNCTION
	double secondBestFitness[(unsigned long)(log(sizeMCMCstring)/log(2))+1];
	for(unsigned long i =0; i<=(unsigned long)(log(sizeMCMCstring)/log(2)); i++)
	    secondBestFitness[i] = 0;
#endif
	
	for(unsigned long i =0; i<sizeMCMCstring; i++){
	    MCMCstring[i] -> nrOfGenerations = 0;
	    MCMCstring[i] -> reset(runs);
	}
	
	while((unsigned long)(log(sizeMCMCstring)/log(2.0)) >= 0)
	{	
	    
	    cout << "Runs:" << runs << " nrOfchains " << sizeMCMCstring << 
		" generation " << generation << " sample size:" << sampleSize  << " \n"; 
	    
	    while(generation <= sampleSize){
		
		//if(generation )
		//cout << "generation:"<< generation << " diversity:" << diversity[0] << "\n";
		
		noCoupling();
		postProcessing();
		//cout << "generation punct:"<< generation << " diversity:" << diversity[0] << "\n";
		
#ifdef ACCEPTANCE_RATIO
		//aceasta este gata numai pentru sizeMCMCstring = 2
		for(unsigned long i = 0; i < sizeMCMCstring; i++){
		    if(sizeMCMCstring == 2 && i == 1)
			acceptance_ratio[1][runs-1] += MCMCstring[i]->accepted;
		    else if(sizeMCMCstring == 2){
			acceptance_ratio[0][runs-1] += MCMCstring[i]->accepted;
			acceptance_ratio[1][runs-1] += MCMCstring[i]->accepted;
		    } else acceptance_ratio[0][runs-1] += MCMCstring[i]->accepted;
		} 
#endif
		
//mean best fitness
      for(unsigned long i = 0; i < sizeMCMCstring; i++){
	  if(MCMCstring[i]->stop() == 1 && firstBestFitness[(unsigned long)(log(sizeMCMCstring)/log(2.0))] != 1){
	      bestFitness += generation;
	      firstBestFitness[(unsigned long)(log(sizeMCMCstring)/log(2.0))] = 1;
#ifndef TWO_ATTRACTORS_FUNCTION
	      bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1] = generation; 
#else
	      bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][0] = generation; 
	      cout << "first 0 best " << i << " size chain size " << sizeMCMCstring <<
		  " in generation " << generation << " in run" << (runs-1) << "\n";
#endif
	  }
#ifdef TWO_ATTRACTORS_FUNCTION
	  if(MCMCstring[i]->stop() == 2 && secondBestFitness[(unsigned long)(log(sizeMCMCstring)/log(2.0))] != 1){
	      bestFitness += generation;
	      secondBestFitness[(unsigned long)(log(sizeMCMCstring)/log(2.0))] = 1;
	      bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][1] = generation;
	      cout << "second 1 best " << i << " size chain size " << sizeMCMCstring <<
		  " in generation " << generation << "in run " << (runs-1) << "\n";
	  }    
#endif
      }
	    
      
#ifdef KULLBACK_INFORMATION
		//cout << "KL diff" << mixingTimeUnivariateCalculation()<< "\n";
		if(generation % 10000 == 0){
		    unsigned long tempSize = (unsigned long)(log(sizeMCMCstring)/log(2.0));
		    simul[tempSize] = 
		    insert_simulation(simul[tempSize],generation, 
				      mixingTimeUnivariateCalculation());
		show_simulations(simul[tempSize]);
		
		if(tempSize == 0){
		    //char* tempSimul = new char[strlen(simul_all_proc_right) + 3];
		    //tempSimul = strcat(,(char)tempSize);
		    write_to_file_all("simul_all_proc_right_no_coupling_0.txt",
				      simul[tempSize],NULL);
		}
		else if(tempSize == 1){
		    write_to_file_all("simul_all_proc_right_no_coupling_1.txt"
				      ,simul[tempSize],NULL);
		    //write_to_file_all("simul_all_proc_right_no_coupling_0.txt",
		    //		 simul[tempSize],NULL);
		}
		else 
		    write_to_file_all(simul_all_proc_right,simul[tempSize],NULL);
		//mixingTimeUnivariateCalculation();
	    }
  
#endif //KULLBACK_INFORMATION
	    
//print the outout for further analyses
#ifdef PLOTS_MCMC
	    //printAllChains_I(myFile);
#endif
	    generation++;	
	} //while(generation <= sampleSize)	
	
	cout << "Begin convergence diagnostic for nr chains:" << sizeMCMCstring <<
	    "which run: " << sampleSize <<"\n";
	
	//scoreDiversity();
	
#ifdef KULLBACK_INFORMATION
	//myMixing << "Nr runs:" << runs << " nr chains:" << sizeMCMCstring << "\n";
	//for(unsigned long i = 0 ; i < sampleSize; i++)
	//    myMixing << Kullback[i] << "\n";
#endif //KULLBACK_INFORMATION
	
	//scrie in fisierul gnufile
	  myData << "Diversity mean" << mean << " nr chains:" << sizeMCMCstring << " \n";
	  myData << "rightSolutions: "<< procentRightSolutions/(sizeMCMCstring*sampleSize)
		 << "\n";	  
	  myData << "Acceptance:" << 
	      acceptance_ratio[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1]/
	      (sampleSize*sizeMCMCstring) <<"\n";
	  
	  for(unsigned long j= 0; j < diversityMeasure; j++){
	      if (diversityTrueDistribution[j] != 0 & diversityRightSolutions[j] != 0) {
		  myData << j <<"\t"<< (double)diversity[j]<<"\t"<< 
		      diversityRightSolutions[j] << "\t" <<
		      log((double)diversityRightSolutions[j]/
			  (double)diversityTrueDistribution[j]) <<"\n";
	      }
	      else {
		  if (diversityTrueDistribution[j] == 0) 
		      myData << j <<"\t"<<(double)diversity[j]<<"\t" << 0 
			     << "\t" << 0 << "\n";
		  else if(diversityRightSolutions[j] == 0) 
		      myData << j <<"\t"<<(double)diversity[j]<<"\t" << 0 
			     << "\t" << -pow(10,30) << "\n";
	      }
	  }
	  
	  //scrie in fisier best fitness din aceasta rulare
	  myData << "\n best Fitness" << bestFitness/nr_runs << "\n";
	  //scrie in fisierul err
	  myFile << "nr chains" << sizeMCMCstring << "\n";
	  
	  printHashtable(myFile);
	  
       	  scoreDiversity();
	  
	  //globalRightSolutions[runs-1] = procentRightSolutions;
	  
#ifdef CUSUM
	  CusumCalculation();

//atentie la prelucrare - scriere pe orinzital- o linie - rulare; 
/*#ifdef BURN_IN
	  unsigned long maxim = burn_inT[0];
	  //for(unsigned long j = 1; j < nr_runs; j++)
	  if(burn_inT[j-1] > maxim) maxim = burn_inT[runs-1];
	  myCUSUM << maxim << "\n";
	  
	  for(unsigned long i = maxim; i < sampleSize; i++){
	      //for(unsigned long j = 0; j < nr_runs; j++)
	      myCUSUM << Cusum[runs-1][i]/(double)(i+1) << "\t" << Hairness[runs-1][i] <<"\t";
	      myCUSUM << "\n";
	  }
#else
	  for(unsigned long i = 0; i < sampleSize; i++){
	      //double tempCusum =  Cusum[runs-1][i]/(double)(i+1);
	      myCUSUM << Cusum[runs-1][i]/(double)(i+1) << "\t";
	  }
	  myCUSUM << "\n";
	  for(unsigned long i = 0; i < sampleSize; i++)
	      myCUSUM << Hairness[runs-1][i] <<"\t";
	  myCUSUM << "\n";
	 
#endif //BURN_IN
*/	  
#endif //CUSUM
	  
#ifdef GELMAN_RUBIN
	  varianceHigerMomentsUnivariateCalculation(GELMAN_Moment,sizeMCMCstring - 1);
	  varianceMultivariateCalculation();
#endif
	  sizeMCMCstring = sizeMCMCstring / 2;
	  sampleSize *= 2;
	  //transforma Hashtable in Hashtable in care sint MCMC ramase
	  updateHashtable();

	  if(sizeMCMCstring != 0){
	      if(bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][0]!=0)
		  firstBestFitness[(unsigned long)(log(sizeMCMCstring)/log(2.0))] = 1;
	      if(bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][1]!=0)
		  secondBestFitness[(unsigned long)(log(sizeMCMCstring)/log(2.0))] = 1;	      
	  }
	  //    myFile << "After cutting half"<< " nr chains" << sizeMCMCstring << "\n";
	  //    
	  //    printHashtable(myFile);
	  //}
      }     //while(generation <= beginSampleSize*beginMCMCstring)
      
      cout << "End convergence diagnostic \n";
      if(runs != nr_runs) {
	  sizeMCMCstring = beginMCMCstring;
	  sampleSize = sampleSize_Const;
	  reset();
      }
      else cout << "Runnnig finished:" << runs << "\t" << MCMCstring[0]->temperature << " \n";
      //myHashtables -> reset();
    }
    
//  printDiversityGNUFILE();
    
    sizeMCMCstring = beginMCMCstring;
    sampleSize = sampleSize_Const;
    double minimumPerformance[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1][diversityMeasure];
    double maximumPerformance[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1][diversityMeasure];
    double stdDev[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1][diversityMeasure];
    for(unsigned long i= 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1; i++)
      for(unsigned long k= 0; k < diversityMeasure; k++){
	  stdDev[i][k] = 0;
	  minimumPerformance[i][k] = 1;
	  maximumPerformance[i][k] = -1000000;
      }

  myData << "\n Statistical data over all runs \n";
  for(unsigned long i= 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1; i++)
      for(unsigned long j= 0; j < nr_runs; j++)
	  for(unsigned long k= 0; k < diversityMeasure; k++)
	      if (diversityTrueDistribution[k] != 0 & 
		  globalDiversityRightSolutions[i][j][k] != 0) {
		  performance[i][k] += log((double)globalDiversityRightSolutions[i][j][k]/ 
					   (double)diversityTrueDistribution[k]);
		  performance_runs[i][k]++;
		  if(log((double)globalDiversityRightSolutions[i][j][k]
			 /(double)diversityTrueDistribution[k]) > maximumPerformance[i][k])
		      maximumPerformance[i][k] = 
			  log((double)globalDiversityRightSolutions[i][j][k]/ 
			      (double)diversityTrueDistribution[k]);

		  if(log((double)globalDiversityRightSolutions[i][j][k]
			 /(double)diversityTrueDistribution[k]) < minimumPerformance[i][k])
		      minimumPerformance[i][k] = 
			  log((double)globalDiversityRightSolutions[i][j][k]/ 
			      (double)diversityTrueDistribution[k]);

		  //  stdDev[i][k] += pow(log((double)globalDiversityRightSolutions[i][j][k]
		  // /(double)sampleData[k]),2);
		  //myData << j <<"\t"<<(double)globalDiversity[i][j]<<"\t" << 
		  //      globalDiversityRightSolutions[runs-1][j] << "\t" <<
		  //      log((double)globalDiversityRightSolutions[runs-1][j]/
		  //	  (double)sampleData[j]) <<"\n";
	      }
	      else {
		  if (diversityTrueDistribution[k] == 0) {}
		  //myData << j <<"\t"<<(double)globalDiversity[runs-1][j]<<"\t" << 0 
		  //	 << "\t" << 0 << "\n";
		  else if(globalDiversityRightSolutions[i][j][k] == 0)
		      minimumPerformance[i][k] = -1000000;
		  //   myData << j <<"\t"<<(double)globalDiversity[runs-1][j]<<"\t" << 0 
		  //	     << "\t" << -pow(10,30) << "\n";
	      }
  
  
  for(unsigned long i= 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1; i++)
      for(unsigned long j= 0; j < nr_runs; j++)
	  for(unsigned long k= 0; k < diversityMeasure; k++)
	      if (diversityTrueDistribution[k] != 0 & 
		  globalDiversityRightSolutions[i][j][k] != 0) {
		  stdDev[i][k] += pow(performance[i][k]/(double)performance_runs[i][k] - 
				      log((double)globalDiversityRightSolutions[i][j][k]
					  /(double)diversityTrueDistribution[k]),2);
	      }

  for(unsigned long i= 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) +1; i++){
      myData << "nr chains:" << pow(2,i) << " \n";
      for(unsigned long k = 0; k < diversityMeasure; k++){
	  //myData << "diversity: " << k << " with performance " <<
	  myData << k << "\t" << performance[i][k]/(double)performance_runs[i][k] << 
		"\t" << sqrt(stdDev[i][k]/(double)performance_runs[i][k]) << "\t" 
		<< maximumPerformance[i][k] << "\t" << minimumPerformance[i][k] << "\n";
	  myMixing << pow(2,i) << "\t" << DEFAULT_TEMPERATURE << "\t" 
		   << k << "\t" << performance[i][k]/(double)performance_runs[i][k] << 
		"\t" << sqrt(stdDev[i][k]/(double)performance_runs[i][k]) << "\t" 
		<< maximumPerformance[i][k] << "\t" << minimumPerformance[i][k] << "\n";
      }
      myMixing << "\n";
  }

  //cout << "Temperature " << MCMCstring[0]->temperature << "\n";
#ifdef ACCEPTANCE_RATIO
  myData << "\n acceptance ratio \n";
  for(unsigned long i=  0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++){
      double accepted = 0;
      myData << "nr chains:" << pow(2,i) << " \n";
      for(unsigned long k = 0; k < nr_runs; k++){
	  myData << "run = " << k  <<" acceptance= " <<
	      acceptance_ratio[i][k]/(sampleSize*sizeMCMCstring) << "\n";
	  accepted += acceptance_ratio[i][k];
      }
      accepted /= (sampleSize*nr_runs*sizeMCMCstring);
#ifdef TWO_ATTRACTORS_FUNCTION
      meanAcceptance[i] = accepted;
#else
      meanAcceptance = accepted;
#endif
      myData << "overall acceptance" << accepted << "\n";

      double st_dev = 0;
      for(unsigned long k = 0; k < nr_runs; k++)
	  st_dev += pow(acceptance_ratio[i][k]/(sampleSize*sizeMCMCstring - accepted),2);

      myAcceptance << sizeGENOME << "\t" << pow(2,i) << "\t" << (sampleSize * sizeMCMCstring) 
		   << "\t" << DEFAULT_TEMPERATURE << "\t" << accepted << "\t" <<
		   sqrt(st_dev/nr_runs) << "\n";
  }
  myAcceptance << "\n";

#endif

//scrie in fisier hitting time best fitness
  myData << "bestFitness \n";
#ifndef TWO_ATTRACTORS_FUNCTION
  for(unsigned long i= 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++){
      double tempBest  = 0;
      for(unsigned long k = 0; k < nr_runs; k++)
	  tempBest += bestFitnessVector[i][k];
      myData << "nr chains:" << pow(2,i) << " \t";
      myData << tempBest/nr_runs << "\n"; 
  }
#else
  for(unsigned long i= 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++){
      myData << "nr chains:" << pow(2,i) << " \n";
      double firstBest = 0;
      double difBest = 0;
      double tempBest  = 0;
      double tempBest2 = 0;
      double minim = sampleSize * pow(2,i);
      double maxim = -10000;
      double st_dev = 0;
      double firstBestVector[nr_runs];
      for(unsigned long k = 0; k < nr_runs; k++){
	  double best = sampleSize*pow(2,i);
	  myData << bestFitnessVector[i][k][0] << "\t" << bestFitnessVector[i][k][1] << "\n";
	  if(bestFitnessVector[i][k][0] != 0){
	      tempBest += bestFitnessVector[i][k][0];
	      if(minim > bestFitnessVector[i][k][0])
		  minim = bestFitnessVector[i][k][0];
	      if(maxim < bestFitnessVector[i][k][0])
		  maxim = bestFitnessVector[i][k][0];
	      best = bestFitnessVector[i][k][0]; 
	  }
	  else tempBest += sampleSize * pow(2,i);
	  if(bestFitnessVector[i][k][1] != 0){
	      tempBest2 += bestFitnessVector[i][k][1];
	      if(minim > bestFitnessVector[i][k][1])
		  minim = bestFitnessVector[i][k][1];
	      if(maxim < bestFitnessVector[i][k][1])
		  maxim = bestFitnessVector[i][k][1];
	      if(best > bestFitnessVector[i][k][1])
		  best = bestFitnessVector[i][k][1];
	  }
	  else tempBest2 += sampleSize * pow(2,i);

	  if(bestFitnessVector[i][k][0] >= bestFitnessVector[i][k][1] && bestFitnessVector[i][k][1] != 0){
	      firstBest += bestFitnessVector[i][k][1];
	      difBest += bestFitnessVector[i][k][0] - bestFitnessVector[i][k][1];
	  } else if(bestFitnessVector[i][k][0] <= bestFitnessVector[i][k][1] && bestFitnessVector[i][k][0] != 0){
	      firstBest += bestFitnessVector[i][k][0];
	      difBest += bestFitnessVector[i][k][1] - bestFitnessVector[i][k][0];
	  } else if(bestFitnessVector[i][k][0] > bestFitnessVector[i][k][1]){
	      firstBest += bestFitnessVector[i][k][0];
	      difBest += sampleSize * sizeMCMCstring - bestFitnessVector[i][k][0];
	 }else if(bestFitnessVector[i][k][0] < bestFitnessVector[i][k][1]){
	      firstBest += bestFitnessVector[i][k][1];
	      difBest += sampleSize*sizeMCMCstring - bestFitnessVector[i][k][1];
	  } else { // do not fine any picks
	      firstBest += sampleSize*sizeMCMCstring;
	      difBest += sampleSize*sizeMCMCstring; 
	  }
	  firstBestVector[k] = best;
      }
      myData << "mean  " << tempBest/nr_runs << ",\t" << tempBest2/nr_runs<< 
	  ",\t" << firstBest/nr_runs << ",\t" << difBest/nr_runs << "\n";

      for(unsigned long k =0; k < nr_runs; k++)
	  st_dev += pow(firstBestVector[k] - firstBest/nr_runs,2); 

      myHitting << sizeGENOME << "\t" << pow(2,i) << "\t" << (sampleSize*sizeMCMCstring) <<
	  "\t" << DEFAULT_TEMPERATURE << "\t" << firstBest/nr_runs << "\t" << 
	  sqrt(st_dev/nr_runs) << "\t" << minim << "\t" << maxim << "\n";

      meanFitness[i] = firstBest/nr_runs;
      diffFitness[i] = difBest/nr_runs;
      //else
//	  meanseconfFitness = 
  }
  myHitting << "\n";
#endif
/*  for(unsigned long j= 0; j < nr_runs; j++){
    cout << "RightSolution "<< j << " is " << (double)globalRightSolutions[j] << "\n";
    myData << j <<"\t"<<((double)globalRightSolutions[j]/(double)(sampleSize*sizeMCMCstring))<<"\n";
    } */
  
/*#ifdef TRAP_LIKE_FUNCTION
  printDiversityDetailGNUFILE();
  #endif
  printDiversityRightSolutions();
  
  #ifdef KULLBACK_INFORMATION
  printMixingGNUFILE();
  #endif*/

  sizeMCMCstring = beginMCMCstring;
  sampleSize = sampleSize_Const;
#ifdef CUSUM
  printCusumGNUFILE();
#endif

//print histogram

  double minimumHistogram[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1][sizeGENOME+1];
  double maximumHistogram[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1][sizeGENOME+1];
  double stdDevHistogram[(unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1][sizeGENOME+1];
  for(unsigned long i= 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1; i++)
      for(unsigned long k= 0; k < sizeGENOME+1; k++){
	  stdDevHistogram[i][k] = 0;
	  minimumHistogram[i][k] = 1;
	  maximumHistogram[i][k] = -1000000;
      }

  //char histogramName1[100], histogramName2[100];
/*  ofstream myHistogram0("histogram_0.file");		
  if(!myHistogram0.is_open()){
      cout << "Error opening histogram file \n";
      return;
  }

  ofstream myHistogram1("histogram_1.file");		
  if(!myHistogram.is_open()){
      cout << "Error opening histogram file \n";
      return;
      }*/

  ofstream myHistogram(fileHistogram, ios::out | ios:: app);		
  if(!myHistogram.is_open()){
      cout << "Error opening histogram file \n";
      return;
  }
  
  myData << "Histogram is = \n";
  
  for(unsigned long i =0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++){
      myData << "nr chains \t" << pow(2,i) << "\n";
      for(unsigned long k = 0; k < sizeGENOME+1; k++){
	  for(unsigned long j = 0; j < nr_runs; j++){
	      globalHistogram[i][k] += ActualHistogram[i][j][k]/histogram[k];
	      if(maximumHistogram[i][k] > ActualHistogram[i][j][k]/histogram[k])
		  maximumHistogram[i][k] = ActualHistogram[i][j][k]/histogram[k];
	      if(minimumHistogram[i][k] < ActualHistogram[i][j][k]/histogram[k])
		  minimumHistogram[i][k] = ActualHistogram[i][j][k]/histogram[k];
	      //if(i==0)
	      //cout <<  pow(2,i) << "\t" << j << "\t" << k << "\t" 
	      //   << ActualHistogram[i][j][k] << "\t" << histogram[k] << "\n";
	      //cout << "nr chains="<< pow(2,i) << "\t diversity=" <<k << "\t runs="
	      //	   << j << " actualHistogram=" << ActualHistogram[i][j][k] << 
	      //	  "real histogram=" << histogram[k] << "\n";
	  }
	  globalHistogram[i][k] /= nr_runs;
	  myData << k << "\t" << globalHistogram[i][k] << "\n";
	  cout << "global histogram[" << i << "][" << k << "]=" 
	       << globalHistogram[i][k] << "\n";
      }
  }
  
  cout << "temperature " << MCMCstring[0]->temperature << "\n";

  for(unsigned long i =0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0))+1; i++){
      for(unsigned long k = 0; k < sizeGENOME+1; k++){
	  for(unsigned long j = 0; j < nr_runs; j++)
	      stdDevHistogram[i][k] += pow(globalHistogram[i][k] - 
					   ActualHistogram[i][j][k]/histogram[k],2);
	  stdDevHistogram[i][k] = sqrt(stdDevHistogram[i][k]/nr_runs);
	  myHistogram << sizeGENOME << "\t" << DEFAULT_TEMPERATURE << "\t" << 
	      (sampleSize * sizeMCMCstring) << "\t"<< k << "\t" << pow(2,i) << "\t" 
		      << globalHistogram[i][k]
		      << "\t" << stdDevHistogram[i][k] << "\t" << maximumHistogram[i][k]
		      << "\t" << minimumHistogram[i][k] << "\n";
      }
      myHistogram << "\n";
  }
  
  myHistogram.flush();
  myHistogram.close();
  
/*#ifdef GELMAN_RUBIN
  printReductionGNUFILE(GELMAN_Moment);
  printReductionMultivariateGNUFILE();
  #endif //GELMAN_RUBIN
*/  
//#ifdef DEBUG_
  myFile.flush();
  myFile.close();
//#endif //DEBUG_
  
  myData.flush();
  myData.close();

  myHitting.flush();
  myHitting.close();

#ifdef CUSUM
  myCUSUM.flush();
  myCUSUM.close();
#endif
  
  myMixing.flush();
  myMixing.close();

  myAcceptance.flush();
  myAcceptance.close();

#ifndef FILE_BOOKEEPING
  diversityFile.flush();
  diversityFile.close();
#endif
}

#endif //SINGLE_vs_MULTIPLE

#else 

#ifdef RANDOM_COUPLING
// ------------Random coupling -----------

void binaryPopulation::evolution(){
  unsigned long tempSize = sizeMCMCstring / sizeCoupling;
  unsigned long** table = new unsigned long*[tempSize + 1];
  for(unsigned long i = 0; i< tempSize + 1; i++)
    table[i] = new unsigned long[sizeCoupling];

  cout << "Random coupling \n";

#ifdef DEBUG_  
  ofstream myFile(fileData);		
  if(!myFile.is_open())
    cout << "Error opening file";
#endif  
 
#ifdef EXPAND_KULLBACK_INFORMATION
    expandTrueUnivariateDistribution();
#else
    readHistogram();
#endif

  while(runs < nr_runs){
    runs++;
    setSeed();

    while(generation <=  burn_in + sampleSize*interSamplesize)
      {	
	assignTemperature();
	unsigned long dimension = sizeMCMCstring/sizeCoupling;
#ifndef NEAR_CHANGING  
	randomCOUPLING(table);
#else
	dimension = neighborsCOUPLING_NO_OVERLAPING(table);
#endif
	proposedGeneration(table,dimension);
	acceptanceProposed(table,dimension);
	postProcessing();

	simul = insert_simulation(simul,generation,procentRightSolutions);
	show_simulations(simul);
	write_to_file_all(simul_all_proc_right,simul,NULL);

#ifdef KULLBACK_INFORMATION
	Kullback[generation] += mixingTimeUnivariateCalculation();      
#endif //KULLBACK_INFORMATION

	generation++;	
      }	

#ifdef DEBUG_
    printLEVEL(myFile);
#else
    scoreDiversity();
#endif

    globalRightSolutions[runs-1] = procentRightSolutions;
#ifdef CUSUM
	CusumCalculation();
#endif

#ifdef GELMAN_RUBIN
  varianceHigerMomentsUnivariateCalculation(GELMAN_Moment,sizeMCMCstring - 1);
  varianceMultivariateCalculation();
#endif

    reset();
  }

#ifdef GELMAN_RUBIN
  printReductionGNUFILE(GELMAN_Moment);
  printReductionMultivariateGNUFILE();
#endif

  printDiversityGNUFILE();

#ifdef TRAP_LIKE_FUNCTION
  printDiversityDetailGNUFILE();
#endif
  printDiversityRightSolutions();

#ifdef CUSUM
	printCusumGNUFILE();
#endif

#ifdef DEBUG_
  myFile.close();
#endif

#ifndef FILE_BOOKEEPING
  diversityFile.flush();
  diversityFile.close();
#endif

  delete[] table;
}


void binaryPopulation::acceptanceProposed(unsigned long** table, unsigned long dimension)
{
  unsigned long tempSize = dimension;
  unsigned long temp[sizeGENOME];
  double tempValue;

#ifdef COUPLED_ACCEPTANCE
  for(unsigned long i = 0; i<tempSize ; i++){
    for(unsigned long j = 0; j<sizeCoupling; j++){
      proposedMCMCstring[table[i][j]]->setACCEPTED(0); 
      MCMCstring[table[i][j]]->setACCEPTED(0);
    }
    if(acceptancePopulation(table[i]) == 1)
      for(unsigned long j = 0; j<sizeCoupling; j++){
	proposedMCMCstring[table[i][j]]->getGENOME(temp);
	double tempValue = proposedMCMCstring[table[i][j]]->getValue();
	MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	MCMCstring[table[i][j]]->setACCEPTED(1);
      }
    //cout << "am acceptat ceva \n";	
    else {
      for(unsigned long j = 0; j<sizeCoupling; j++)
	MCMCstring[table[i][j]] -> nextGenerationBlank();
      //cout << "nu am acceptat ceva \n";	
    }      
  }
#else
  for(unsigned long i =0; i <sizeMCMCstring; i++){
    proposedMCMCstring[i] -> getGENOME(temp);
    tempValue =  proposedMCMCstring[i] -> fitness();
    double tempPropDistrValue = proposedMCMCstring[i] -> ProposalDistribution(temp,sizeGENOME,proposalDistribution);
    proposedMCMCstring[i] -> setProposalDistribution(tempPropDistrValue);

    if(MCMCstring[i] -> acceptance(temp,tempValue,tempPropDistrValue) == 1)
	MCMCstring[i] -> nextGeneration(temp,tempValue,sizeGENOME);
      else 
	MCMCstring[i] -> nextGenerationBlank();
  }
  //MCMCstring[i] -> nextGeneration();
#endif // coupled acceptance
}

unsigned long binaryPopulation::acceptancePopulation(unsigned long* table){
  double prop = 1;
  for(unsigned long i = 0; i<sizeCoupling; i++) 
    prop = prop * 
      exp((MCMCstring[table[i]]->getValue() - 
	   proposedMCMCstring[table[i]]-> getValue())/MCMCstring[table[i]]->getTEMPERATURE());
  
  if(prop > 1) return 1;
  double temp = genrand_int32();
  if(temp > prop) {
    return 0;
  }
  return 1;
}

#else 
#ifdef PARALLEL_TEMPERING
//-------------Parallel Tempering ------------

#ifdef MULTIPLE_RUNS
void binaryPopulation::evolution(ofstream& myMultiFile){
#else		
void binaryPopulation::evolution(){
#endif
  unsigned long tempSize = sizeMCMCstring / sizeCoupling;
    unsigned long** table = new unsigned long*[tempSize + 1];
    for(unsigned long i = 0; i< tempSize + 1; i++)
	table[i] = new unsigned long[sizeCoupling];
    //int table[tempSize + 1][sizeCoupling];

#ifndef MIXTURE_REALS
  unsigned long tempGENOME[sizeGENOME];
#else
  double tempGENOME[sizeGENOME];
#endif
  cout << "PARALLEL_TEMPERING\n";

#ifndef TRACE_RUN 
#ifndef MULTIPLE_RUNS
#ifdef HISTOGRAM
  ofstream myFileFG("err_recomb_fit_secondGen");		
  if(!myFileFG.is_open()){
      cout << "Error opening file for data err_mut_dist_firstGen \n";
    return;
  }

  ofstream myFileFit("err_recomb_fit_firstGen");		
  if(!myFileFit.is_open()){
      cout << "Error opening file for data err_mut_fit_firstGen \n";
    return;
  }

  ofstream myFileHash("err_recomb_hash_firstGen");		
  if(!myFileHash.is_open()){
    cout << "Error opening file for data err_mut_hash_firstGen \n";
    return;
  }

//#ifdef DEBUG_ 
  ofstream myFile(fileData);		
  if(!myFile.is_open()){
      cout << "Error opening file for data " << fileData << "\n";
    return;
  }
//#endif  

  ofstream myTempFile("temporar.txt");		
  if(!myTempFile.is_open()){
      cout << "Error opening file for data " << fileData << "\n";
    return;
  }

  ofstream myFileHash("err_recomb_hash_firstGen");		
  if(!myFileHash.is_open()){
    cout << "Error opening file for data err_mut_hash_firstGen \n";
    return;
  }
#else //histogram

  int records;
  double mse, cor, s1, s2, jump1, jumpGen;
  double tempPrev[sizeMCMCstring][sizeGENOME];
  double tempValuePrev[sizeMCMCstring];
  if(nr_runs == 1)
    records = (int)((double)sampleSize_Const/ (double)SEE_SAMPLE) + 1;
  else records = nr_runs;
  
  double averageFitnessIndivEach = 0;
  double sqrAverageFitnessIndivEach = 0;
  double averageVariablesIndivEach[sizeGENOME];
  double sqrAverageVariablesIndivEach[sizeGENOME];
  for(int j = 0; j < sizeGENOME; j++){
    averageVariablesIndivEach[j] = 0;
    sqrAverageVariablesIndivEach[j] = 0;
  }

  double averageFitnessIndivRuns[records];
  double averageFitnessRuns[records];
  double averageVariablesIndivRuns[records][sizeGENOME];
  double corVariablesIndivRuns[records][sizeGENOME * (sizeGENOME - 1) / 2];
  double sqrCorVariablesIndivRuns[records][sizeGENOME * (sizeGENOME - 1) / 2];
  
  double sqrAverageFitnessIndivRuns[records];
  double sqrAverageFitnessRuns[records];

  double sqrAverageVariablesIndivRuns[records][sizeGENOME];
  for(int j = 0; j < records; j++){
    averageFitnessIndivRuns[j] = 0;
    averageFitnessRuns[j] = 0;
    for(int k = 0; k < sizeGENOME; k++)
      averageVariablesIndivRuns[j][k] = 0;
    sqrAverageFitnessIndivRuns[j] = 0;
    sqrAverageFitnessRuns[j] = 0;
    for(int k = 0; k < sizeGENOME; k++)
      sqrAverageVariablesIndivRuns[j][k] = 0;
    for(int k = 0; k < sizeGENOME * (sizeGENOME - 1) / 2; k++){
      corVariablesIndivRuns[j][k] = 0;
      sqrCorVariablesIndivRuns[j][k] = 0;
    }
  }
  
  double averageFitnessIndiv[records][sizeMCMCstring];
  double averageVariablesIndiv[records][sizeMCMCstring][sizeGENOME];
  double corVariablesIndiv[records][sizeMCMCstring][sizeGENOME * (sizeGENOME - 1) / 2];
  double sqrCorVariablesIndiv[records][sizeMCMCstring][sizeGENOME * (sizeGENOME - 1) / 2];
  double averageJumpIndiv[records][sizeMCMCstring], sqrAverageJumpIndiv[records][sizeMCMCstring];
  double averageJumpIndivRuns[records], sqrAverageJumpIndivRuns[records];
  double averageJump[records], sqrAverageJump[records];
  //double averageAcceptanceIndivRuns[records], sqrAverageAcceptanceIndivRuns[records];
   double averageExpJumpIndivRuns[records];
  double sqrAverageExpJumpIndivRuns[records];
   
  for(int j = 0; j < records; j++){
    for(int k = 0; k < sizeMCMCstring; k++){
      averageFitnessIndiv[j][k] = 0;
      for(int i = 0; i < sizeGENOME; i++)
	averageVariablesIndiv[j][k][i] = 0;
      for(int i = 0; i < sizeGENOME * (sizeGENOME - 1) / 2; i++){
	corVariablesIndiv[j][k][i] = 0;
	sqrCorVariablesIndiv[j][k][i] = 0;
      }
      averageJumpIndiv[j][k] = 0;
      sqrAverageJumpIndiv[j][k] = 0;
    }
    averageJumpIndivRuns[j] = 0;
    sqrAverageJumpIndivRuns[j] = 0;
    averageExpJumpIndivRuns[j] = 0;
    sqrAverageExpJumpIndivRuns[j] = 0;
    averageJump[j] = 0;
    sqrAverageJump[j] = 0;
  }
  double sqrAverageFitnessIndiv[records][sizeMCMCstring];
  double sqrAverageVariablesIndiv[records][sizeMCMCstring][sizeGENOME];
  for(int j = 0; j < records; j++)
    for(int k = 0; k < sizeMCMCstring; k++){
      sqrAverageFitnessIndiv[j][k] = 0;
      for(int i = 0; i < sizeGENOME; i++)
	sqrAverageVariablesIndiv[j][k][i] = 0;
    }

  double averageFitnessIndivAll[sizeMCMCstring];
  double averageVariablesIndivAll[sizeMCMCstring][sizeGENOME];
  for(int j = 0; j < sizeMCMCstring; j++){
    averageFitnessIndivAll[j] = 0;
    for(int i = 0; i < sizeGENOME; i++)
      averageVariablesIndivAll[j][i] = 0;
  }

  double sqrAverageFitnessIndivAll[sizeMCMCstring];
  double sqrAverageVariablesIndivAll[sizeMCMCstring][sizeGENOME];
  for(int j = 0; j < sizeMCMCstring; j++){
    sqrAverageFitnessIndivAll[j] = 0;
    for(int i = 0; i < sizeGENOME; i++)
      sqrAverageVariablesIndivAll[j][i] = 0;
  }

  ofstream* myFileAllSplit = new ofstream[records];
  char nameFileSplit[30] = "outRecombAllChainSplit";		
  for(int i = 0; i < records; i++){
    nameFileSplit[22] = 48 + (int)(i / 10);
    nameFileSplit[23] = 48 + i - 10 * (int)(i / 10);
    nameFileSplit[24] = '\0'; 
    myFileAllSplit[i].open(nameFileSplit);
    if(!myFileAllSplit[i].is_open()){
      cout << "Error opening file for data output_mut_all_chains_splits \n";
      return;
    } 
  }
 
  ofstream myFileAllNoSplit("outRecombAllChainNoSplit");		
  if(!myFileAllNoSplit.is_open()){
      cout << "Error opening file for data err_mut_fit_firstGen \n";
    return;
  }

  ofstream* myFile1ChainSplit = new ofstream[records];
  nameFileSplit[9] = 'O'; nameFileSplit[10] = 'n'; nameFileSplit[11] = 'e'; 		
  for(int i = 0; i < records; i++){
    nameFileSplit[22] = 48 + (int)(i / 10);
    nameFileSplit[23] = 48 + i - 10 * (int)(i / 10);
    nameFileSplit[24] = '\0'; 
    myFile1ChainSplit[i].open(nameFileSplit);
    if(!myFile1ChainSplit[i].is_open()){
      cout << "Error opening file for data output_mut_all_chains_splits \n";
      return;
    } 
  }

  ofstream myFile1ChainNoSplit("outRecombOneChainNoSplit");		
  if(!myFile1ChainNoSplit.is_open()){
      cout << "Error opening file for data err_mut_fit_firstGen \n";
    return;
  }

  ofstream myFileHash("outRecombHashAllChainNoSplit");		
  if(!myFileHash.is_open()){
    cout << "Error opening file for data err_mut_hash_firstGen \n";
    return;
  }
#endif //histogram

#endif //MULTIPLE_RUNS
#else //trace run

    ofstream myFile("err_recomb_2p_0.1_500");
    ofstream myFileI("err_recomb_2p_0.1_1000");
    ofstream myFileII("err_recomb_2p_0.1_1500");
    ofstream myFileIII("err_recomb_2p_0.1_2000");
    ofstream myFileIV("err_recomb_2p_0.1_2500");

#endif //TRACE_RUN

#ifdef ACCEPTANCE_RATIO
//de calculat average
#ifndef multiple_restarts
      double myAverage[(unsigned long)(sampleSize_Const/SEE_SAMPLE) + 1];
      double myAcceptanceT[(unsigned long)(sampleSize_Const/SEE_SAMPLE) + 1];
      
      double nrHighIndiv[(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
      unsigned long distHighIndiv[(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];

      double maxIndiv[(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
      for(unsigned long j = 0; j < (unsigned long)(sampleSize_Const/SEE_SAMPLE)+1; j++){
	  maxIndiv[j] = - 2500;
	  distHighIndiv[j] = 0;
      }
#else 
      double myAverage[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE) + 1];
      double myAcceptanceT[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE) + 1];
      
      double nrHighIndiv[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
      unsigned long distHighIndiv[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];

      double maxIndiv[nr_runs][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
      for(int i=0; i < nr_runs; i++)
	for(unsigned long j = 0; j < (unsigned long)(sampleSize_Const/SEE_SAMPLE)+1; j++){
	  maxIndiv[i][j] = - 2500;
	  distHighIndiv[i][j] = 0;
	}
#endif //nr_runs
    double maximI = -2500;
#endif

#ifndef TRACE_RUN 

#ifndef MIXTURE_REALS
#ifdef DISTANCE
    double myDistance[nr_runs][sizeGENOME+1];
    double count[sizeGENOME+1], tempDistance[sizeGENOME+1];
    for(unsigned long i = 0; i < sizeGENOME+1; i++){
	for(unsigned long j = 0; j < nr_runs; j++){
	    myDistance[j][i] = 0;
	}
	count[i] = 0;
	tempDistance[i] = 0;
    }
    ofstream myFileDistance("distance_par_temp.txt");
#endif
#else //mixture reals
#ifndef HISTOGRAM
    unsigned long *discreteGenome = new unsigned long[sizeGENOME];
#endif //histogram
#endif //mixture reals

#ifdef TRAP_LIKE_FUNCTION
    //store 1111 blocks
    double detailBLOCKS[nr_runs][sizeGENOME/BLOCKsize][(unsigned long)(sampleSize_Const/SEE_SAMPLE)+1];
#endif

#ifdef HISTOGRAM
#ifdef EXPAND_KULLBACK_INFORMATION
  //cout << "Expand distribution \n";  
  expandTrueUnivariateDistribution();
#else
    readHistogram();
#endif
#endif //histogram

#endif //TRACE_RUN

#ifdef SIMPLEX
    for(unsigned long i = 0; i < sizeMCMCstring; i++){
      MCMCstring[i] -> proportions = new double[sizeCoupling];
      double sumProp = 0;
      for(unsigned long j = 0; j < sizeCoupling - 1; j++){
	sumProp += pow(1.0/2.0,(double)j)/sizeGENOME; 
	MCMCstring[i] -> proportions[j] = sumProp;
      }
      MCMCstring[i] -> proportions[sizeCoupling - 1] = 1;
    }
#endif //simplex 
 
  cout << "Start runs\n";
  while(runs < nr_runs){
    runs++;

    generation=0;

#ifndef  TRACE_RUN
    //assignTemperature();
    //setSeed();
 
    double firstBestFitness = 0;
    double secondBestFitness = 0;

    diversityFile << "\n \n runs " << runs << "\n";
#endif //TRACE_RUN

    unsigned long genIndiv = 0;
    unsigned long mutPar = 0, recombPar = 0;
    while(generation <=  sampleSize)
    {	
	// GA process
	//cout << generation << "\n";
	if(generation != 0){
	    unsigned long dimension1 = sizeMCMCstring/sizeCoupling;

#ifndef NEAR_CHANGING  
	    randomCOUPLING(table);
#else
	    dimension1 = neighboursCOUPLING_NO_OVERLAP(table);
#endif

	    /* cout << " genemeome before proposal  first ";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++)
		cout << tempGENOME[nr_lenght] << "\t";
	      cout << " \t\t second \t";
	    }
	    cout << "\n"; */

	    proposedGeneration(table,dimension1);

	    /*cout << " genemeome before acceptance first ";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++)
		cout << tempGENOME[nr_lenght] << "\t";
	      cout << " \t\tsecond \t";
	    }
	    cout << "\n";*/

	    acceptanceProposed(table,dimension1);
	}

#ifdef TRACE_RUN
	
	  if(generation <= SEE_SAMPLE ){
	      //cout << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      myFile << generation << "\t";
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFile << tempGENOME[nr_lenght] << "\t";
		cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFile << "\n";
	      //cout << "\n";
	    }
	    //myFile << "\n";
	  } 
	  if(generation <= 2 *  SEE_SAMPLE ){
	    //myFileI << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      myFileI << generation << "\t";
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFileI << tempGENOME[nr_lenght] << "\t";
		//cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFileI << "\n";
	      //cout << "\n";
	    }
	    //myFileI << "\n";
	  }
	  if(generation <= 3 * SEE_SAMPLE ){
	    //myFileII << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      myFileII << generation << "\t";
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFileII << tempGENOME[nr_lenght] << "\t";
		//cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFileII << "\n";
	      //cout << "\n";
	    }
	    //myFileII << "\n";
	  }
	  if(generation <= 4 * SEE_SAMPLE ){
	    //myFileIII << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      myFileIII << generation << "\t";
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFileIII << tempGENOME[nr_lenght] << "\t";
		//cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFileIII << "\n";
	      //cout << "\n";
	    }
	    //myFileIII << "\n";
	  }
	  if(generation <= 5 * SEE_SAMPLE ){
	    //myFileIV << generation << "\t";
	    for(unsigned long nr_indiv = 0; nr_indiv < sizeMCMCstring; nr_indiv++){ 
	      myFileIV << generation << "\t";
	      MCMCstring[nr_indiv] ->getGENOME(tempGENOME);
	      //cout << "generation" << generation << "\n";
	      for(unsigned long nr_lenght = 0; nr_lenght < sizeGENOME; nr_lenght++){
		myFileIV << tempGENOME[nr_lenght] << "\t";
		//cout << tempGENOME[nr_lenght] << "\t";
	      }
	      myFileIV << "\n";
		//cout << "\n";
	    }
	    //myFileIV << "\n";
	  }
	  
#else //trace run

	postProcessing();

	averageFitness += averageGeneration();
	
#ifdef MIXTURE_REALS
	
	//write them in one long file
	//cout << "runs=" << runs << " generation=" << generation << "\t";
#ifdef BURN_IN_GEN
	if(generation > throw_gen){
#endif
	  //cout << "generation " << generation << "\n";
	  jumpGen = 0;
	  for(unsigned long ter = 0; ter < sizeMCMCstring/sizeCoupling; ter++)
	    for(unsigned long ter_j = 0; ter_j < sizeCoupling; ter_j++){
	      //find the coupling and its position
	      double sizeString = sizeMCMCstring;
	      unsigned long i = table[ter][ter_j];
	      MCMCstring[i] -> getGENOME(tempGENOME);
	      double tempValue = MCMCstring[i] -> getValue(); //fitness();
	      
#ifdef IsRecombination
#ifdef SIMPLEX
#ifdef IsMutation
	      if(mix_mutation != 1 && MCMCstring[table[ter][0]] -> mix_rotation == 0 && ter_j != MCMCstring[table[ter][0]] -> positionChanged){
		//cout << "skip " << generation << " indiv (" << tempGENOME[0] << "," << tempGENOME[1] << "," << tempValue << ")  prev (" 
		  //   << tempPrev[i][0] << " , " << tempPrev[i][1] << " , " << tempValuePrev[i] << " ) \n";
		continue; 
	      }
	      else if(mix_mutation != 1 && MCMCstring[table[ter][0]] -> mix_rotation == 0){
		sizeString = sizeMCMCstring/sizeCoupling;
		//i = table[ter][0];
		genIndiv++;
		recombPar++;
	      } else if(mix_mutation != 1 && ter_j != 0){
		//cout << "skip " << generation << " indiv (" << tempGENOME[0] << "," << tempGENOME[1] << "," << tempValue << ")  prev (" 
		//   << tempPrev[i][0] << " , " << tempPrev[i][1] << " , " << tempValuePrev[i] << " ) \n";
		continue; 
	      }
	      else if(mix_mutation != 1){
		sizeString = sizeMCMCstring/sizeCoupling;
		//i = table[ter][0];
		genIndiv++;
		recombPar++;
	      } else 
#ifdef OneMut 
		if(mix_mutation == 1 && ter_j != 0)
		  continue; else if(mix_mutation == 1){
		    mutPar ++; 
		    genIndiv++;
		  }
#else
	      if(mix_mutation == 1){
		mutPar ++; 
		genIndiv++;
	      }

#endif //one mut
#else
	      if(MCMCstring[table[ter][0]] -> mix_rotation == 0 && ter_j != MCMCstring[table[ter][0]] -> positionChanged){
		//cout << "skip " << generation << " indiv (" << tempGENOME[0] << "," << tempGENOME[1] << "," << tempValue << ")  prev (" 
		// << tempPrev[i][0] << " , " << tempPrev[i][1] << " , " << tempValuePrev[i] << " ) \n";
		continue;
	      } else if(MCMCstring[table[ter][0]] -> mix_rotation != 0 && ter_j != 0){
		continue;
	      }
	      sizeString = sizeMCMCstring/sizeCoupling;
	      //i = table[ter][0];
#endif //IsMutation
#else //simplex 
#ifdef SINGLE_INDIV_ROT
#ifdef IsMutation
	      if(mix_mutation != 1 && ter_j != 0){
		//cout << "skip " << generation << " indiv (" << tempGENOME[0] << "," << tempGENOME[1] << "," << tempValue << ")  prev (" 
		//   << tempPrev[i][0] << " , " << tempPrev[i][1] << " , " << tempValuePrev[i] << " ) \n";
		continue; 
	      }
	      else if(mix_mutation != 1){
		sizeString = sizeMCMCstring/sizeCoupling;
		//i = table[ter][0];
		recombPar ++; 
		genIndiv++;
	      } else {
		mutPar ++;
		genIndiv++;
	      }
#else
	      if(ter_j != 0){
		//cout << "skip " << generation << " indiv (" << tempGENOME[0] << "," << tempGENOME[1] << "," << tempValue << ")  prev (" 
		// << tempPrev[i][0] << " , " << tempPrev[i][1] << " , " << tempValuePrev[i] << " ) \n";
		continue;
	      }
	      sizeString = sizeMCMCstring/sizeCoupling;
	      //i = table[ter][0];
#endif 
#endif //single indiv rot
#endif //simplex

#ifdef wPCTX
	      if(mix_mutation != 1 && ter_j != 0)
		continue; 
	      else if(mix_mutation != 1){
		sizeString = sizeMCMCstring/sizeCoupling;
		//i = table[ter][0];
	      }
#endif //wPCTX

#ifndef PCTX
#ifndef SIMPLEX
#ifndef SINGLE_INDIV_ROT
#ifndef HALF_PCTX 
#ifdef IsMutation
	      if(mix_mutation != 1 && ter_j != 0){
		//cout << "skip " << generation << " indiv (" << tempGENOME[0] << "," << tempGENOME[1] << "," << tempValue << ")  prev (" 
		//   << tempPrev[i][0] << " , " << tempPrev[i][1] << " , " << tempValuePrev[i] << " ) \n";
		continue; 
	      }
	      else if(mix_mutation != 1){
		sizeString = sizeMCMCstring/sizeCoupling;
		//i = table[ter][0];
		genIndiv++;
	      } else genIndiv+=sizeCoupling;
#else
	      if(ter_j != 0){
		//cout << "skip " << generation << " indiv (" << tempGENOME[0] << "," << tempGENOME[1] << "," << tempValue << ")  prev (" 
		// << tempPrev[i][0] << " , " << tempPrev[i][1] << " , " << tempValuePrev[i] << " ) \n";
		continue;
	      }
	      sizeString = sizeMCMCstring/sizeCoupling;
	      //i = table[ter][0];
#endif 
#endif //HALF_PCTX
#endif // SINGLE_INDIV_ROT
#endif //simplex
#endif //not PCTX
#else // only mut
	      genIndiv++;
	      mutPar ++;
#endif //recombination

#ifdef DEBUG_
	      if(tempValue != MCMCstring[i] -> fitness()){
		cout << "true fitness not equal with the read one " << tempValue << " != " << MCMCstring[i] -> fitness() << " in gen " << generation << "\n";
		exit(1);
	      }

	      if(tempValue < MIN_FITNESS){
		cout << " smaller then thresold " << tempValue << " != " << tempGENOME[0] << " , " << tempGENOME[1] << " , " << generation << "\n";
		
		for(int i1 = 0; i1 < sizeMCMCstring; i1++)
		  tempGENOME[i1] = tempPrev[i][i1];
		tempValue = tempValuePrev[i];
		if(tempValue < MIN_FITNESS){
		  cout << " AFTER smaller then thresold " << tempValue << " != " << tempGENOME[0] << " , " << tempGENOME[1] << " , " << generation << "\n";
		  exit(1);
		}
		
		for(int i1 = 0; i1 < sizeCoupling; i1++)
		  MCMCstring[table[ter][i1]] -> accepted = 0;
		//exit(1);
	      } //else {
		//cout << "modify " << generation << " indiv (" << tempGENOME[0] << "," << tempGENOME[1] << "," << tempValue << ")  prev (" 
		//<< tempPrev[i][0] << " , " << tempPrev[i][1] << " , " << tempValuePrev[i] << " ) \n";
	      //}
	      
	      /*for(int i1 = 1; i1 < sizeCoupling; i1++){
		if(tempValuePrev[table[ter][i1]] < MIN_FITNESS){
		  cout << "temp value wrong in wrong indiviudal"<<tempValuePrev[table[ter][i1]]<<" != "<<tempPrev[i1][0]<< " , " << tempPrev[i1][1]<<" , "<<generation << "\n";
		  exit(1);
		}
		MCMCstring[table[ter][i1]] -> getGENOME(tempGENOME);
		if(tempValuePrev[table[ter][i1]] != MCMCstring[table[ter][i1]] ->value || tempPrev[table[ter][i1]][0] !=  tempGENOME[0] ||  tempPrev[table[ter][i1]][1] !=  tempGENOME[1]){
		  cout << " AFTER temp wrong " <<  table[ter][i1] << " , " << tempValuePrev[table[ter][i1]]<<" != "
		       <<tempPrev[table[ter][i1]][0]<< " , " << tempPrev[table[ter][i1]][1]<<" , "<<generation << "\n";
		  cout << "all values \n";
		  for(int i2 = 0; i2 < sizeCoupling; i2++){
		    MCMCstring[table[ter][i2]] -> getGENOME(tempGENOME);
		    cout << "(" << table[ter][i2] << ",(" << tempGENOME[0] << " , " << tempGENOME[1] << ") != " 
			 << tempPrev[table[ter][i2]][0] << " , " << tempPrev[table[ter][i2]][1] << ") \n ";   
		  }
		  exit(1);
		}
	      }

	      MCMCstring[i] -> getGENOME(tempGENOME);
	      double tempValue = MCMCstring[i] -> getValue(); //fitness();
	      */

#endif //debug_
	  //myString.insert(i*(2*sizeGENOME + 2) + 2*(sizeGENOME - 1) + 1,);

#ifndef HISTOGRAM
	    //discretize before introduce it in one of the hashtable
	    //position in histogram
	  //cout << "transform dicrete";
	  for(unsigned long i1 =0; i1 < sizeGENOME; i1++ ){
	    // if(tempGENOME[i] + scale/2.0 < 0 || tempGENOME[i] + scale/2.0 > scale){
	    // cout << "wrong individual ? (" << i << "," << discreteGenome[i] << "," << tempGENOME[i] << ") = " << MCMCstring[0]->fitness(tempGENOME,sizeGENOME) << "\n"; 
	    // exit(0);
	    //}
	    discreteGenome[i1] = (unsigned long)((tempGENOME[i1] + scale/2.0)/((double)scale/(double)NrBins));
	    //cout << "(" << i << " , " << tempGENOME[i] << " , " << discreteGenome[i] << ")\t"; 
	    if(discreteGenome[i1] >= NrBins) {
	      cout << "wrong individual ? (" << i1 << ","<< generation << "," << discreteGenome[i1] << "," << tempGENOME[i1] << "," << MCMCstring[i]->getValue()<< ")"; 
	      exit(1);
	    }
	  }
	  //cout << "\n";
	  if(nr_runs ==1){
	    //split in more consecutive and independent regions;
	    myHashtables[(int)((double)generation/ (double)SEE_SAMPLE)]->store_individual(discreteGenome,tempValue);
	  } else {
	    //more runs and for each run a hashtable
	    myHashtables[runs-1]->store_individual(discreteGenome,tempValue);
	    //if(generation % sampleSize_Const == sampleSize_Const - 1)
	    //	myFileFit << "\n";	      
	  }
#endif //not histogram

#ifdef WRITE_OUT_FILE
	      //cout << " gene " << generation << " i " << i << " (" << tempGENOME[0] << "," << tempGENOME[1] << "," << tempGENOME[2] << ") " << tempValue << "\n";  
	      for(unsigned long j = 0; j < sizeGENOME; j++){

	      if(nr_runs == 1){
	      //no splited file
		averageVariablesIndivAll[i][j] += tempGENOME[j];
		sqrAverageVariablesIndivAll[i][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndivAll[i][j]/(generation + 1) - pow(averageVariablesIndivAll[i][j]/(generation + 1),2.0);
		myFileAllNoSplit << tempGENOME[j] << "\t" << averageVariablesIndivAll[i][j]/(generation + 1) << "\t"<< mse << "\t";
		
		averageVariablesIndivEach[j] += tempGENOME[j];
		sqrAverageVariablesIndivEach[j] += pow(tempGENOME[j],2.0);	      
		mse = sqrAverageVariablesIndivEach[j]/(generation*sizeString + i  + 1) - pow(averageVariablesIndivEach[j]/(generation*sizeString + i + 1),2.0);
		myFile1ChainNoSplit << tempGENOME[j] << "\t" << averageVariablesIndivEach[j]/(generation*sizeString + i + 1) << "\t" << mse << "\t";
		
	      //splited files
		averageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j] += tempGENOME[j];
		sqrAverageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j]/ (1 + (int)((double)generation / (double)SEE_SAMPLE)) 
		  - pow(averageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j]/ (1 + (int)((double)generation / (double)SEE_SAMPLE)),2.0);
		myFileAllSplit[(int)((double)generation / (double)SEE_SAMPLE)] << tempGENOME[j] << "\t" << 
		  averageVariablesIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i][j]/ (1 + (int)((double)generation / (double)SEE_SAMPLE)) << "\t" << mse << "\t";

		averageVariablesIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)][j] += tempGENOME[j];
		sqrAverageVariablesIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)][j]/(1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeString + i)
	     - pow(averageVariablesIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)][j]/(1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeString + i),2.0);
		myFile1ChainSplit[(int)((double)generation /(double) SEE_SAMPLE)] << tempGENOME[j] << "\t" << 
		  averageVariablesIndivRuns[(int)((double)generation/(double)SEE_SAMPLE)][j]/(1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeString + i) 
										  << "\t" << mse <<"\t";
	      } else {
	      //splited files
		averageVariablesIndiv[runs-1][i][j] += tempGENOME[j];
		sqrAverageVariablesIndiv[runs-1][i][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndiv[runs-1][i][j]/(generation  - throw_gen+ 1) - pow(averageVariablesIndiv[runs-1][i][j]/(generation  - throw_gen+ 1),2.0);
		if(mse < -0.0001){
		  cout << " mse of average variable runs " << runs << " genome " << i << " generation " << generation << " averageVariable "
		       << averageVariablesIndiv[runs-1][i][j] <<" sqr " << sqrAverageVariablesIndiv[runs-1][i][j] << " mse " << mse << " in negative \n";
		  //exit(1);
		}
		if(mse<0) mse = -mse;
		mse = sqrt(mse);

		myFileAllSplit[runs-1] << tempGENOME[j] << "\t" << averageVariablesIndiv[runs-1][i][j]/(generation  - throw_gen+ 1) << "\t" << mse << "\t";

		averageVariablesIndivRuns[runs-1][j] += tempGENOME[j];
		sqrAverageVariablesIndivRuns[runs-1][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndivRuns[runs-1][j]/((generation -1 - throw_gen)*sizeString+i+1)-pow(averageVariablesIndivRuns[runs-1][j]
												    /((generation -1 - throw_gen)*sizeString + i+1),2.0);
		//if(mse < -0.0001){
		//  cout << " mse of average variable runs all indiv" << runs << " genome " << i << " generation " << generation << " averageVariable "
		//       << averageVariablesIndivRuns[runs-1][j] << " sqr " << sqrAverageVariablesIndivRuns[runs-1][j] << " mse " << mse << " in negative \n";
		//  //exit(1);
		//}
		if(mse<0) mse = -mse;
		mse = sqrt(mse);

		myFile1ChainSplit[runs-1] << tempGENOME[j] << "\t" << averageVariablesIndivRuns[runs-1][j]/ ((generation  - throw_gen)* sizeString + i) 
					  << "\t" << mse << "\t";

		//no splited file
		averageVariablesIndivAll[i][j] += tempGENOME[j];
		sqrAverageVariablesIndivAll[i][j] += pow(tempGENOME[j],2.0);
		mse = sqrAverageVariablesIndivAll[i][j]/(generation  - throw_gen+ 1) - pow(averageVariablesIndivAll[i][j]/(generation - throw_gen + 1),2.0);
		myFileAllNoSplit << tempGENOME[j] << "\t" << averageVariablesIndivAll[i][j]/(generation - throw_gen + 1) << "\t"<< mse << "\t";
		
		averageVariablesIndivEach[j] += tempGENOME[j];
		sqrAverageVariablesIndivEach[j] += pow(tempGENOME[j],2.0);	      
		mse = sqrAverageVariablesIndivEach[j]/((generation - throw_gen)*sizeString + i  + 1) - 
		  pow(averageVariablesIndivEach[j]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0);
		myFile1ChainNoSplit << tempGENOME[j] << "\t" << averageVariablesIndivEach[j]/((generation - 1 - throw_gen)*sizeString + i + 1) << "\t" << mse << "\t";
		
	      }
	      //	      cout << tempGENOME[j] << "\t";
	    }
	    
	    //no splited files
	    averageFitnessIndivAll[i] += tempValue;
	    sqrAverageFitnessIndivAll[i] += pow(tempValue,2.0);
	    mse = sqrAverageFitnessIndivAll[i]/(generation  - throw_gen+ 1) - pow(averageFitnessIndivAll[i]/(generation  - throw_gen+ 1),2.0);
	    myFileAllNoSplit << tempValue << "\t" << averageFitnessIndivAll[i]/(generation  - throw_gen+ 1) << "\t" << mse << "\t";

	    averageFitnessIndivEach += tempValue;
	    sqrAverageFitnessIndivEach += pow(tempValue,2.0);
	    mse = sqrAverageFitnessIndivEach/((generation - 1 - throw_gen)*sizeString + i + 1) - 
	      pow(averageFitnessIndivEach/((generation - 1 - throw_gen)*sizeString + i + 1),2.0);
	    myFile1ChainNoSplit << tempValue << "\t" << averageFitnessIndivEach/((generation - 1 - throw_gen)*sizeString + i + 1) << "\t" << mse << "\n";
	    //cout << tempValue << "\n";

	    //splited files
	    if(nr_runs == 1){
	      averageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i] += tempValue;
	      sqrAverageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i] += pow(tempValue,2.0);
	      mse = sqrAverageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i]/(1 + (int)((double)generation / (double)SEE_SAMPLE)) - 
		pow(averageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i]/ (1 + (int)((double)generation / (double)SEE_SAMPLE)),2.0);
	      myFileAllSplit[(int)((double)generation / (double)SEE_SAMPLE)] << tempValue << "\t" << 
		averageFitnessIndiv[(int)((double)generation / (double)SEE_SAMPLE)][i]/ (1 + (int)((double)generation / (double)SEE_SAMPLE)) << "\t" << mse << "\t";

	      averageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)] += tempValue;
	      sqrAverageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)] += pow(tempValue,2.0);
	      mse = sqrAverageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)]/(1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeString + i) - 
		pow(averageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)]/ (1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeString + i),2.0); 
	      myFile1ChainSplit[(int)((double)generation /(double)SEE_SAMPLE)] << tempValue << "\t" << 
		averageFitnessIndivRuns[(int)((double)generation / (double)SEE_SAMPLE)]/ (1 + (int)((double)generation / (double)SEE_SAMPLE) * sizeString + i) << "\t"
									       << mse << "\n";
	    } else {
	      //file all split
	      //corelation computation
	      for(int i1 = 0; i1 < sizeGENOME; i1++)
		for(int i2 = i1+1; i2 < sizeGENOME; i2++){
		  int tempIndex = i1 * (sizeGENOME -1) - i1 * (i1-1)/2 + i2 - i1 - 1;
		  corVariablesIndiv[runs-1][i][tempIndex] += tempGENOME[i1]* tempGENOME[i2];
		  double ave1 = averageVariablesIndiv[runs-1][i][i1]/(generation - throw_gen+ 1);
		  double ave2 = averageVariablesIndiv[runs-1][i][i2]/(generation  - throw_gen+ 1);
		  sqrCorVariablesIndiv[runs-1][i][tempIndex] += pow((tempGENOME[i1] - ave1)* (tempGENOME[i2] - ave2),2.0);
		  cor = corVariablesIndiv[runs-1][i][tempIndex]/(generation - throw_gen+1) - ave1 * ave2;
		  mse = sqrCorVariablesIndiv[runs-1][i][tempIndex]/(generation - throw_gen+1) - pow(cor,2.0);

		  s1 = sqrAverageVariablesIndiv[runs-1][i][i1]/(generation  - throw_gen+ 1) - pow(averageVariablesIndiv[runs-1][i][i1]/(generation  - throw_gen+ 1),2.0);
		  s2 = sqrAverageVariablesIndiv[runs-1][i][i2]/(generation - throw_gen+ 1) - pow(averageVariablesIndiv[runs-1][i][i2]/(generation - throw_gen+ 1),2.0);

		  //if(mse < -0.0001){
		  //cout << " mse of average correlation runs " << runs << " genome " << i << " generation " << generation << " averageCor "
		  // << cor << " sqr " << sqrCorVariablesIndiv[runs-1][i][tempIndex] << " mse " << mse << " in negative \n";
		    //exit(1);
		  //}
		  if(mse<0) mse = -mse;
		  mse = sqrt(mse);
		  if(s1 > 0 && s2 > 0)
		    myFileAllSplit[runs-1] << cor << "\t" << mse << "\t" << cor/sqrt(s1 * s2) << "\t";
		  else {
		    if(s1 < -0.0001 || s2 < -0.0001){
		      cout << " s1 " << s1 << " s2 " << s2 << "\n";
		      exit(1);
		    }
		    myFileAllSplit[runs-1] << cor << "\t" << mse << "\t" << 0 << "\t";
		  }
		}
	      
	      jump1 = 0;
	      if(generation != 0){
		//cout << "gen = " << i << " acc " << MCMCstring[i]->accepted << " ";
		for(int i1 = 0; i1 < sizeGENOME; i1++){
		  jump1 += pow(tempGENOME[i1] - tempPrev[i][i1],2.0);
		  //cout << "("<<tempGENOME[i1] << "," << tempPrev[i][i1] << "," << pow(tempGENOME[i1] - tempPrev[i][i1],2.0) << ","<< jump1 << ")";
		}
		//cout << "\n";
		jump1 = jump1;

		//averge jumping for each chains
		averageJumpIndiv[runs-1][i] += sqrt(jump1/sizeGENOME/scale);
		jumpGen += sqrt(jump1/sizeGENOME/scale);
		sqrAverageJumpIndiv[runs-1][i] += jump1/sizeGENOME/scale;
		mse = sqrAverageJumpIndiv[runs-1][i]/(generation- throw_gen + 1) - pow(averageJumpIndiv[runs-1][i]/(generation - throw_gen+ 1),2.0); 
		//if(mse < -0.0001){
		//cout << " mse of average jumps runs " << runs << " genome " << i << " generation " << generation << " averageJump "
		//     << averageJumpIndiv[runs-1][i] << " sqr " << sqrAverageJumpIndiv[runs-1][i] << " mse " << mse << " in negative \n";
		  //exit(1);
		//}
		// mse = sqrt(abs(mse));
		if(mse<0) mse = -mse;
		mse = sqrt(mse);
	      }
	      if(generation != 0)
		myFileAllSplit[runs-1] << averageJumpIndiv[runs-1][i]/(generation - throw_gen+ 1) << "\t" << mse << "\t";
	      else 
		myFileAllSplit[runs-1] << "0\t0\t";		

	      averageFitnessIndiv[runs-1][i] += tempValue;
	      sqrAverageFitnessIndiv[runs-1][i] += pow(tempValue,2.0);
	      mse = sqrAverageFitnessIndiv[runs-1][i]/(generation - throw_gen+ 1) - pow(averageFitnessIndiv[runs-1][i]/(generation- throw_gen  +1),2.0);
	      myFileAllSplit[runs-1] << tempValue << "\t" << averageFitnessIndiv[runs-1][i]/(generation - throw_gen +1) << "\t" << mse << "\t"; 

	      //file one chain split
	      for(int i1 = 0; i1 < sizeGENOME; i1++)
		for(int i2 = i1+1; i2 < sizeGENOME; i2++){
		  int tempIndex = i1 * (sizeGENOME -1) - i1 * (i1-1)/2 + i2 - i1 - 1;
		  corVariablesIndivRuns[runs-1][tempIndex] += tempGENOME[i1]* tempGENOME[i2];
		  double ave1 = averageVariablesIndivRuns[runs-1][i1]/((generation - 1 - throw_gen)*sizeString + i + 1);
		  double ave2 = averageVariablesIndivRuns[runs-1][i2]/((generation - 1 - throw_gen)*sizeString + i + 1);
		  sqrCorVariablesIndivRuns[runs-1][tempIndex] += pow((tempGENOME[i1] - ave1) * (tempGENOME[i2] - ave2),2.0);

		  mse = sqrCorVariablesIndivRuns[runs-1][tempIndex]/((generation - 1 - throw_gen)*sizeString + i + 1) - 
		    pow(corVariablesIndivRuns[runs-1][tempIndex]/((generation - 1 - throw_gen)*sizeString + i + 1) - ave1 * ave2,2.0);

		  cor = corVariablesIndivRuns[runs-1][tempIndex]/((generation - 1 - throw_gen)*sizeString + i + 1) - ave1 * ave2;
		  s1 = sqrAverageVariablesIndivRuns[runs-1][i1]/((generation - 1 - throw_gen)*sizeString + i + 1) - pow(ave1,2.0);
		  s2 = sqrAverageVariablesIndivRuns[runs-1][i2]/((generation - 1 - throw_gen)*sizeString + i + 1) - pow(ave2,2.0);

		  //if(mse < -0.0001){
		  //cout << " mse of average jumps runs " << runs << " genome " << i << " geenration " << generation << " averageJump "
		  // << corVariablesIndivRuns[runs-1][tempIndex] << " sqr " << sqrCorVariablesIndivRuns[runs-1][tempIndex] << " mse " << mse << " in negative \n";
		    //exit(1);
		  //}
		  // mse = sqrt(abs(mse));
		  if(mse<0) mse = -mse;
		  mse = sqrt(mse);

		  if(s1 > 0 && s2 > 0)
		    myFile1ChainSplit[runs-1] << cor << "\t" << mse << "\t" << cor/sqrt(s1 * s2) << "\t";
		  else 
		    myFile1ChainSplit[runs-1] << cor << "\t" << mse << "\t" << 0 << "\t";
		}

	      if(generation != 0){
		averageJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		sqrAverageJumpIndivRuns[runs-1] += pow(jump1/sizeGENOME/scale,2.0);
		mse = sqrAverageJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
	      }
	      if(generation != 0)
		myFile1ChainSplit[runs-1] << averageJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1) << "\t" << mse << "\t";
	      else 
		myFile1ChainSplit[runs-1] << "0\t0\t";		

#ifdef MIXTURE
	      if(mix_mutation == 1){
#ifdef IsMutation
	      if(generation != 0){
		if(tempValuePrev[i] < tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		}
		else{ 
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
	      }
#else 
#ifdef IsRecombination
		  int theOthers[sizeCoupling - 1];
		  for(int coup = 0; coup < sizeString/sizeCoupling; coup++)
		    for(int dim = 0; dim < sizeCoupling; dim++)
		      if(table[coup][dim] == i){
			for(int inter = 1; inter < sizeCoupling; inter++)
			  theOthers[inter-1] = table[coup][(dim+inter)%sizeCoupling];
			break;
		      }
		  double theOtherValue = 1;
		  for(int inter = 0; inter< sizeCoupling-1; inter++)
		    theOtherValue *= MCMCstring[theOthers[inter]]->fitness();
		  
		  double theOtherPrevValue = 1;
		  
		  for(int inter = 0; inter< sizeCoupling-1; inter++)
		    theOtherPrevValue *= tempValuePrev[inter];
		  
		  if(tempValuePrev[i] * theOtherPrevValue < tempValue * theOtherValue){
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME); 
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME;
		  }
		  else{ 
#ifndef JUMP_SINGLE
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue / tempValuePrev[i];
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue/tempValuePrev[i], 2.0);
#else
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue * theOtherValue /(tempValuePrev[i] * theOtherPrevValue);
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue* theOtherValue/(tempValuePrev[i]* theOtherPrevValue), 2.0);
#endif 
		  }
		  mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		    -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
		  
#endif //is recombination
#endif //is recombination		
	      } else{
#ifdef IsRecombination 
#ifndef PCTX
		if(tempValuePrev[i] < tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		}
		else { 
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
#else
#ifdef SINGLE_INDIV_ROT
		if(tempValuePrev[i] < tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		}
		else{ 
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
#else
#ifdef wPCTX
#ifndef SUM_DISTRIB
		if(tempValuePrev[i] * (double)proposalDistribution[table[ter][0]] < (double)proposedProposalDistribution[table[ter][0]] * tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		} else{
		  averageExpJumpIndivRuns[runs-1] += 
		    sqrt(jump1/sizeGENOME/scale)*tempValue/tempValuePrev[i]*(double)proposedProposalDistribution[table[ter][0]]/(double)proposalDistribution[table[ter][0]];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * 
		    pow(tempValue/tempValuePrev[i] * (double)proposedProposalDistribution[table[ter][0]]/(double)proposalDistribution[table[ter][0]], 2.0);
		}

		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
#else
		double theOtherValue = 1;
		double theOtherPrevValue = 1;

		for(int inter = 0; inter< sizeCoupling; inter++)
		  theOtherValue += MCMCstring[table[ter][inter]]->fitness();		
		
		for(int inter = 0; inter< sizeCoupling; inter++)
		  theOtherPrevValue += tempValuePrev[table[ter][inter]];

		if(theOtherPrevValue * (double)proposalDistribution[table[ter][0]] < (double)proposedProposalDistribution[table[ter][0]] * theOtherValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME;
		}
		  else{ 
		    averageExpJumpIndivRuns[runs-1] += 
		      sqrt(jump1/sizeGENOME)*theOtherValue/theOtherPrevValue*(double)proposedProposalDistribution[table[ter][0]]/(double)proposalDistribution[table[ter][0]];
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * 
		      pow(theOtherValue/theOtherPrevValue * (double)proposedProposalDistribution[table[ter][0]]/(double)proposalDistribution[table[ter][0]], 2.0);
		  }
		  mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		    -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 		  
#endif //sum prop 
#else
		double theOtherValue = 1;
		double theOtherPrevValue = 1;

#ifndef SUM_DISTRIB

		for(int inter = 0; inter< sizeCoupling; inter++)
		  theOtherValue *= MCMCstring[table[ter][inter]]->fitness();		
		
		for(int inter = 0; inter< sizeCoupling; inter++)
		  theOtherPrevValue *= tempValuePrev[table[ter][inter]];
#else
		for(int inter = 0; inter< sizeCoupling; inter++)
		  theOtherValue += MCMCstring[table[ter][inter]]->fitness();		
		
		for(int inter = 0; inter< sizeCoupling; inter++)
		  theOtherPrevValue += tempValuePrev[table[ter][inter]];

#endif 		  
		if(theOtherPrevValue < theOtherValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME;
		}
		else{ 
#ifndef JUMP_SINGLE
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue / tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue/tempValuePrev[i], 2.0);
#else
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * theOtherValue /theOtherPrevValue;
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(theOtherValue/theOtherPrevValue, 2.0);
#endif
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 		  
#endif  //rot
#endif //wPCTX
#endif //ndef PCTX
#else 
#ifdef IsMutation
		if(tempValuePrev[i] < tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		}
		else{ 
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
#endif 
#endif 
	    }
#else //mixture
#ifdef CYCLE
	      if(generation % 2 == 0) {
#ifdef IsMutation
	      if(generation != 0){
		if(tempValuePrev[i] < tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		}
		else{ 
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
	      }
#else 
#ifdef IsRecombination
	      if(generation != 0){
		  int theOthers[sizeCoupling - 1];
		  for(int coup = 0; coup < sizeString/sizeCoupling; coup++)
		    for(int dim = 0; dim < sizeCoupling; dim++)
		      if(table[coup][dim] == i){
			for(int inter = 1; inter < sizeCoupling; inter++)
			  theOthers[inter-1] = table[coup][(dim+inter)%sizeCoupling];
			break;
		      }
		  double theOtherValue = 1;
		  for(int inter = 0; inter< sizeCoupling-1; inter++)
		    theOtherValue *= MCMCstring[theOthers[inter]]->fitness();
		  
		  double theOtherPrevValue = 1;
		  
		  for(int inter = 0; inter< sizeCoupling-1; inter++)
		    theOtherPrevValue *= tempValuePrev[inter];
		  
		  if(tempValuePrev[i] * theOtherPrevValue < tempValue * theOtherValue){
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME); 
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME;
		  }
		  else{ 
#ifndef JUMP_SINGLE
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue / tempValuePrev[i];
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue/tempValuePrev[i], 2.0);
#else
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue * theOtherValue /(tempValuePrev[i] * theOtherPrevValue);
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue* theOtherValue/(tempValuePrev[i]* theOtherPrevValue), 2.0);
#endif 
		  }
		  mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		    -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
	      }
#endif //is recombination
#endif //is recombination		
	      } else 
#ifdef IsRecombination 
		if(generation != 0){
		  int theOthers[sizeCoupling - 1];
		  for(int coup = 0; coup < sizeString/sizeCoupling; coup++)
		    for(int dim = 0; dim < sizeCoupling; dim++)
		      if(table[coup][dim] == i){
			for(int inter = 1; inter < sizeCoupling; inter++)
			  theOthers[inter-1] = table[coup][(dim+inter)%sizeCoupling];
			break;
		      }
		  double theOtherValue = 1;
		  for(int inter = 0; inter< sizeCoupling-1; inter++)
		    theOtherValue *= MCMCstring[theOthers[inter]]->fitness();
		  
		  double theOtherPrevValue = 1;
		  
		  for(int inter = 0; inter< sizeCoupling-1; inter++)
		    theOtherPrevValue *= tempValuePrev[inter];
		  
		  if(tempValuePrev[i] * theOtherPrevValue < tempValue * theOtherValue){
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME); 
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME;
		  }
		  else{ 
#ifndef JUMP_SINGLE
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue / tempValuePrev[i];
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue/tempValuePrev[i], 2.0);
#else
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue * theOtherValue /(tempValuePrev[i] * theOtherPrevValue);
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue* theOtherValue/(tempValuePrev[i]* theOtherPrevValue), 2.0);
#endif 
		  }
		  mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		    -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
		}
#else 
#ifdef IsMutation
	      if(generation != 0){
		if(tempValuePrev[i] < tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		}
		else{ 
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
	      }
#endif 
#endif 
#else 
	      //cycle in proposal dsitribution 
	      //cout << "should enter here\n";
#ifdef IsRecombination
#ifdef PCTX 
#ifndef SINGLE_INDIV_ROT
		if(generation != 0){
		  int theOthers[sizeCoupling - 1];
		  for(int coup = 0; coup < sizeString/sizeCoupling; coup++)
		    for(int dim = 0; dim < sizeCoupling; dim++)
		      if(table[coup][dim] == i){
			for(int inter = 1; inter < sizeCoupling; inter++)
			  theOthers[inter-1] = table[coup][(dim+inter)%sizeCoupling];
			break;
		      }
		  double theOtherValue = 1;
		  for(int inter = 0; inter< sizeCoupling-1; inter++)
		    theOtherValue *= MCMCstring[theOthers[inter]]->getValue();
		  
		  double theOtherPrevValue = 1;
		  
		  for(int inter = 0; inter< sizeCoupling-1; inter++)
		    theOtherPrevValue *= tempValuePrev[inter];
		  
		  if(tempValuePrev[i] * theOtherPrevValue < tempValue * theOtherValue){
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME); 
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME;
		  }
		  else{ 
#ifndef JUMP_SINGLE
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue / tempValuePrev[i];
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue/tempValuePrev[i], 2.0);
#else
		    averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue * theOtherValue /(tempValuePrev[i] * theOtherPrevValue);
		    sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue* theOtherValue/(tempValuePrev[i]* theOtherPrevValue), 2.0);
#endif
		  }
		  mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		    -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
		}
#else //#ifdef SINGLE_INDIV_ROT
		if(MCMCstring[i] -> mix_rotation == 1){
		  if(generation != 0){
		    if(tempValuePrev[i] < tempValue){
		      averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		      sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		    }
		    else{ 
		      averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		      sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		    }
		    mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		      -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
		  }
		} else {
		  if(generation != 0){
		    int theOthers[sizeCoupling - 1];
		    for(int coup = 0; coup < sizeString/sizeCoupling; coup++)
		      for(int dim = 0; dim < sizeCoupling; dim++)
			if(table[coup][dim] == i){
			  for(int inter = 1; inter < sizeCoupling; inter++)
			    theOthers[inter-1] = table[coup][(dim+inter)%sizeCoupling];
			  break;
			}
		    double theOtherValue = 1;
		    for(int inter = 0; inter< sizeCoupling-1; inter++)
		      theOtherValue *= MCMCstring[theOthers[inter]]->getValue();
		    
		    double theOtherPrevValue = 1;
		    
		    for(int inter = 0; inter< sizeCoupling-1; inter++)
		      theOtherPrevValue *= tempValuePrev[inter];
		    
		    if(tempValuePrev[i] * theOtherPrevValue < tempValue * theOtherValue){
		      averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME); 
		      sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME;
		    }
		    else{ 
#ifndef JUMP_SINGLE
		      averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue / tempValuePrev[i];
		      sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue/tempValuePrev[i], 2.0);
#else
		      averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME) * tempValue * theOtherValue /(tempValuePrev[i] * theOtherPrevValue);
		      sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME * pow(tempValue* theOtherValue/(tempValuePrev[i]* theOtherPrevValue), 2.0);
#endif
		    }
		    mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		      -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
		  }
		}
#endif //#ifdef SINGLE_INDIV_ROT

#else //not PCTX
	      if(generation != 0){
		if(tempValuePrev[i] < tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		}
		else{ 
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
	      }
#endif //PCTX
#else 
#ifdef IsMutation
	      if(generation != 0){
		if(tempValuePrev[i] < tempValue){
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale); 
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale;
		}
		else{ 
		  averageExpJumpIndivRuns[runs-1] += sqrt(jump1/sizeGENOME/scale) * tempValue /tempValuePrev[i];
		  sqrAverageExpJumpIndivRuns[runs-1] += jump1/sizeGENOME/scale * pow(tempValue /tempValuePrev[i], 2.0);
		}
		mse = sqrAverageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1)
		  -pow(averageExpJumpIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0); 
	      }
#endif //mutation
#endif //recombniation
#endif //cycle
#endif //mixture
  
	      if(generation != 0)
		myFile1ChainSplit[runs-1] << averageExpJumpIndivRuns[runs-1]/((generation- throw_gen)*sizeString+i+1) << "\t" << mse << "\t";
	      else 
		myFile1ChainSplit[runs-1] << "0\t0\t";		

	      averageFitnessIndivRuns[runs-1] += tempValue;
	      sqrAverageFitnessIndivRuns[runs-1] += pow(tempValue,2.0);
	      mse = sqrAverageFitnessIndivRuns[runs-1]/ ((generation - 1 - throw_gen)*sizeString + i + 1) 
		- pow(averageFitnessIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1),2.0);
	      myFile1ChainSplit[runs-1] << tempValue << "\t" << averageFitnessIndivRuns[runs-1]/((generation - 1 - throw_gen)*sizeString + i + 1) << "\t" << mse << "\n"; 
	    }
#endif //#ifdef WRITE_OUT_FILE
	    }
	
	  for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    MCMCstring[i] ->getGENOME(tempGENOME);
	    for(int i1 = 0; i1 < sizeGENOME; i1++)
	      tempPrev[i][i1] = tempGENOME[i1];
	    tempValuePrev[i] = MCMCstring[i]->getValue();
	  }

#ifdef WRITE_OUT_FILE

	myFileAllNoSplit << averageGeneration()  << "\t" << averageFitness/(generation+1)<< "\n";
	
	if(nr_runs == 1){
	  averageFitnessRuns[(int)((double)generation / (double)SEE_SAMPLE)] += averageGeneration();
	  myFileAllSplit[(int)((double)generation / (double)SEE_SAMPLE)] << averageGeneration() << "\t"
	       << averageFitnessRuns[(int)((double)generation / (double)SEE_SAMPLE)]/(1 + (int)((double)generation / (double)SEE_SAMPLE)) << "\n";
	} else {
	  averageFitnessRuns[runs-1] += averageGeneration();
	  myFileAllSplit[runs-1] << averageGeneration()  << "\t" << averageFitnessRuns[runs-1]/(generation - throw_gen + 1)<< "\n";
	  //with the means and the std until now make a 

	  //compute the jumps for all chains
	    averageJump[runs-1] += jumpGen/sizeMCMCstring;
	    sqrAverageJump[runs-1] += pow(jumpGen/sizeMCMCstring,2.0);
	    mse = sqrAverageJump[runs-1]/(generation - throw_gen+ 1) - pow(averageJump[runs-1]/(generation - throw_gen+ 1),2.0);
	    myFileAllSplit[runs-1] << averageJump[runs-1]/(generation - throw_gen+ 1) << "\t" << mse << "\n";   
	}
#endif //#ifdef WRITE_OUT_FILE    	
	//how to introduce individuals in hashtable when it is not a diversity measure ?!?
#ifdef BURN_IN_GEN
	}else if(generation == throw_gen) {
	  /*for(unsigned long i =0 ; i < sizeMCMCstring; i++){
	      MCMCstring[i] -> getGENOME(tempGENOME);
	      double tempValue = MCMCstring[i] -> fitness();
	      //cout << "transform dicrete" << i << "(" << tempGENOME[0] << "," << tempGENOME[1] << ")\n";
	      for(unsigned long jack =0; jack < sizeGENOME; jack++ ){
		discreteGenome[jack] = (unsigned long)((tempGENOME[jack] + scale/2.0)/((double)scale/(double)NrBins));
		//cout << "(" << i << " , " << tempGENOME[i] << " , " << discreteGenome[i] << ")\t"; 
		if(discreteGenome[jack] >= NrBins) cout << "wrong individual ? (" << jack << "," << discreteGenome[jack] << "," << tempGENOME[jack] << ")"; 
	      }
	      //cout << "\n";
	      if(nr_runs ==1){
		//split in more consecutive and independent regions;
		myHashtables[(int)((double)generation/ (double)SEE_SAMPLE)]->store_individual(discreteGenome,tempValue);
	      } else {
		//more runs and for each run a hashtable
		myHashtables[runs-1]->store_individual(discreteGenome,tempValue);
		//if(generation % sampleSize_Const == sampleSize_Const - 1)
		//	myFileFit << "\n";	      
	      }
	      }*/
	  for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    MCMCstring[i] ->getGENOME(tempGENOME);
	    for(int i1 = 0; i1 < sizeGENOME; i1++)
	      tempPrev[i][i1] = tempGENOME[i1];
	    tempValuePrev[i] = MCMCstring[i]->getValue();
	  }
	}
	
#endif
#else // is not mixture reals

#ifndef MULTIPLE_RUNS
	//scrie fitness
	if(runs == 1){
	  for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    MCMCstring[i] -> getGENOME(tempGENOME);
	    myFileFit << MCMCstring[i] -> getValue() << "\t";
	  }
	  myFileFit << "\n";
	} else if(runs == 2){
	  for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    MCMCstring[i] -> getGENOME(tempGENOME);
	    myFileFG << MCMCstring[i] -> getValue() << "\t";
	  }
	  myFileFG << "\n";
	  }
	myTempFile << MCMCstring[0]->getValue() << "\n" << averageGeneration() << "\n";
	
	//scrie distance dintre indivizi buni descoperiti pina acum
	//descrie hashtable
	
#endif	
#endif //mixture reals

#ifndef FILE_BOOKEEPING
	if(generation % SEE_SAMPLE == 0)
	    printDiversityPopulation();
#endif

#ifdef HISTOGRAM
#ifdef KULLBACK_INFORMATION
//	    Kullback[generation] += temp;
	if(generation % SEE_SAMPLE == 0){
	    double temp = mixingTimeUnivariateCalculation();
	    //Kullback[(int)(generation/SEE_SAMPLE)] = temp;
	    simul = insert_simulation(simul,generation,temp);
	    show_simulations(simul);
	    write_to_file_all(simul_all_proc_right,simul,NULL);
	}
#endif //KULLBACK_INFORMATION
#endif //histogram

#ifdef ACCEPTANCE_RATIO
	for(unsigned long i = 0; i < sizeMCMCstring; i++)
	    if(maximI < MCMCstring[i]->getValue()){
		maximI = MCMCstring[i]->getValue();
	    }
#ifndef multiple_restarts
	  acceptance_ratio[(int)((double)generation / (double)SEE_SAMPLE)] += acceptanceAverage();
	  //cout << " aceptance ratio " << acceptance_ratio[(int)((double)generation / (double)SEE_SAMPLE)] << "\n";
	  //averageFitness += averageGeneration();
	  if(generation % SEE_SAMPLE == 0){
	    //printAcceptance(myAcceptance);
	    myAverage[(unsigned long)(generation/SEE_SAMPLE)] = averageGeneration();	
	    myAcceptanceT[(unsigned long)(generation/SEE_SAMPLE)] = acceptanceAverage();	
	    nrHighIndiv[(unsigned long)(generation/SEE_SAMPLE)] = myHashtables[runs-1]->nrElem;
	    //	    distHighIndiv[runs-1][(int)(generation/SEE_SAMPLE)] = nrHighIndivGlobal;
	    //distHighIndiv[runs-1][(int)(generation/SEE_SAMPLE)] = averageDistanceGoodIndiv();
	    maxIndiv[(int)(generation/SEE_SAMPLE)] = maximI;
	  }
#else
	  acceptance_ratio[runs-1] += acceptanceAverage();
	  averageFitness += averageGeneration();
	  if(generation % SEE_SAMPLE == 0){
	    //printAcceptance(myAcceptance);
	    myAverage[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = averageGeneration();	
	    myAcceptanceT[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = acceptanceAverage();	
	    nrHighIndiv[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = myHashtables[runs-1]->nrElem;
	    //	    distHighIndiv[runs-1][(int)(generation/SEE_SAMPLE)] = nrHighIndivGlobal;
	    //distHighIndiv[runs-1][(int)(generation/SEE_SAMPLE)] = averageDistanceGoodIndiv();
	    maxIndiv[runs-1][(int)(generation/SEE_SAMPLE)] = maximI;
	  } 
#endif //nr_runs
#endif
	
#ifdef HISTOGRAM
#ifdef DISTANCE
	//int tempGenome[sizeGENOME];
	//MCMCstring[0]->getMax(tempGenome);
	double averageDistance = 0;
	for(unsigned long i = 0; i < sizeMCMCstring; i++){
#ifndef MIXTURE_REALS
	    unsigned long ate = (unsigned long)HammingDistance(i);
	    myDistance[runs-1][ate] += MCMCstring[i]->fitness();
	    count[ate]++;
#else
	    double ate = (double)HammingDistance(i);
#endif // MIXTURE_REALS
	    averageDistance += ate;
#ifndef MULTIPLE_RUNS
	if(i == 0)
	    myTempFile << ate << "\n";
#endif	
	}
	averageDistance/=sizeMCMCstring;
#ifndef MULTIPLE_RUNS
	//if(runs == 1)
	myTempFile << averageDistance << "\n";
#endif
#endif //DISTANCE
#endif //histogram

#ifdef TRAP_LIKE_FUNCTION
 	if(generation % SEE_SAMPLE == 0){
	  
           for(unsigned long i = 0 ; i < sizeGENOME/BLOCKsize; i++)
#ifdef BINOMIAL

	    detailBLOCKS[runs-1][i][(unsigned long)(sampleSize_Const/SEE_SAMPLE)] = diversityDetail[i][BLOCKsize];
#else
	    detailBLOCKS[runs-1][i][(unsigned long)(sampleSize_Const/SEE_SAMPLE)] = diversityDetail[i][a];
#endif
	}
#endif //TRAP_LIKE_FUNCTION

#ifdef KNOWN_MAXIM
#ifndef TWO_ATTRACTORS_FUNCTION
	for(unsigned long i=0; i < sizeMCMCstring; i++){
	    if(MCMCstring[i]->stop() == 1 && firstBestFitness == 0){
		bestFitness += generation;
		firstBestFitness = 1;
		bestFitnessVector[runs-1] = generation; 
	    } else if(MCMCstring[i]->stop() == 2 && secondBestFitness == 0){
		//bestFitness += generation;
		secondBestFitness = 1;
		secondBestFitnessVector[runs-1] = generation; 
	    } 
	}
#else
	for(unsigned long i=0; i < sizeMCMCstring; i++){
	    if(MCMCstring[i]->stop() == 1 && firstBestFitness == 0){
		bestFitness += generation;
		firstBestFitness = 1;
		bestFitnessVector[runs-1][0] = generation; 
	    } else if(MCMCstring[i]->stop() == 2 && secondBestFitness == 0){
		//bestFitness += generation;
		secondBestFitness = 1;
		bestFitnessVector[runs-1][1] = generation; 
	    }	
	}
#endif //BINOMIAL
#endif //KNOWN_MAXIM

#endif //TRACE_RUN	
	generation++;	
    }	

#ifndef TRACE_RUN
#ifdef HISTOGRAM
   
#ifdef GELMAN_RUBIN
    varianceHigerMomentsUnivariateCalculation(GELMAN_Moment,sizeMCMCstring - 1);
    varianceMultivariateCalculation();
#endif
  
#ifdef DEBUG_
  //printLEVEL(myFile);
#else 
    /*if(runs == 1){
      myFile << "\n nr chains" << sizeMCMCstring << " runs " << runs-1 <<  
	  " lenght MCMC " << generation << "\n";
      
//  printHashtable(myFile);
      
      printAllChains(myFile);
      }*/
    myFileHash << "max " << MCMCstring[0]->getMaxValue() << "\n";
    printHashtable(myFileHash);
    
//    scoreDiversity();
#endif

#ifdef TRAP_LIKE_FUNCTION
    scoreDiversityDetail();
#endif
  
  globalRightSolutions[runs-1] = procentRightSolutions;

#ifndef MIXTURE_REALS
#ifdef DISTANCE
	for(unsigned long i = 0; i < sizeGENOME+1; i++){
	    if(count[i] != 0){
		myDistance[runs-1][i] /= count[i];
		count[i] = 0;
	    } else{
		if(myDistance[runs-1][i] != 0)
		    cout<<"Greseala in distanta "<<myDistance[runs-1][i]<<"i =" << i << " \n";
	    }
	}
#endif
#endif //mixture reals

//	firstBestFitness = 0;
//	secondBestFitness = 0;

#endif //histogram

#ifdef CUSUM
  CusumCalculation();
#endif

#endif //NOT trace run
  cout << " generated individuals " << genIndiv << "\t" << " mut " << mutPar << "\t" << " recomb " << recombPar << "\n"; 
  genIndiv = 0;
  mutPar = 0;
  recombPar = 0;
  if(runs != nr_runs) reset();
  }
  
#ifndef TRACE_RUN

  /// write in the file the hashtables -> in the one that supose to write the values of the second run
  //cout << "arive here 2 \n";
#ifdef MIXTURE_REALS
#ifndef HISTOGRAM
  //write the general statistical data
  //cout << " 3 \n";
  
  double populationVariance, populationMean, centerMass;
  cout << " runs " << nr_runs << " sampleSize " << sampleSize_Const << " sample " << SEE_SAMPLE << " \n";
    if(nr_runs ==1){
    //split in more consecutive and independent regions;
    myHashtables[0]->printStd(myFileHash,myHashtables,(int)((double)sampleSize_Const/ (double)SEE_SAMPLE) + 1,centerMass,populationVariance);
  } else {
    //more runs and for each run a hashtable
    myHashtables[0]->printStd(myFileHash,myHashtables,nr_runs,centerMass,populationVariance);
    }  
    //cout << "Is here \n";
    //double temp = mixingTimeUnivariateCalculation();
  //write the data in a file for histogram
  ofstream myHistogram(fileHistogram,ios::app);		
  if(!myHistogram.is_open()){
      cout << "Error opening histogram file \n";
      return;
  }
  
  myHistogram << "# par tempering \n";
  myHistogram << sizeGENOME << "\t" << sampleSize_Const << "\t" << nr_runs << "\t" << (int)((double)sampleSize_Const/ (double)SEE_SAMPLE) + 1 << "\t"
	      << centerMass << "\t" << populationVariance << "\n";

  myHistogram.flush();
  myHistogram.close();
#endif //histogram
#endif //mixutre reals

#ifdef HISTOGRAM
#ifdef GELMAN_RUBIN
  printReductionGNUFILE(GELMAN_Moment);
  printReductionMultivariateGNUFILE();
#endif
  
  printDiversityGNUFILE();

#ifdef TRAP_LIKE_FUNCTION
  printDiversityDetailGNUFILE();
#endif

  printDiversityRightSolutions();
  
//#ifdef KULLBACK_INFORMATION
//  printMixingGNUFILE();
//#endif
    
#ifndef FILE_BOOKEEPING
  diversityFile.flush();
  diversityFile.close();
#endif

#endif //histogram  
 
#ifdef CUSUM
    double hair = printCusumGNUFILE();
#endif

#ifdef ACCEPTANCE_RATIO
    ofstream myAcceptance(fileAcceptance, ios::out | ios::app);
    if(!myAcceptance.is_open()){
	cout << "Error in opening the acceptance file \n";
	exit(1);
    }

    double cumA = 0, cumV = 0, cumN = 0;;
    double meanAv = 0, std_meanAv = 0;

#ifdef multiple_restarts
      for(unsigned long i = 1; i < (sampleSize_Const/SEE_SAMPLE) + 1; i++){
	double meanA = 0, std_devA = 0;
	double meanV = 0, std_devV = 0;
	double meanN = 0, std_devN = 0;
	double cusum_mean = 0, hairness_mean = 0;
	double std_dev_c = 0, std_dev_h = 0;
	double meanD = 0, std_devD = 0;
	double meanM = 0, std_devM = 0;
	
	for(unsigned long j = 0; j < nr_runs; j++){
	  meanA += myAcceptanceT[j][i];
	  meanV += myAverage[j][i];
	  meanN += nrHighIndiv[j][i];
	  meanD += distHighIndiv[j][i];
	  if(i !=sampleSize_Const/SEE_SAMPLE){
	    cusum_mean += Cusum[j][i*SEE_SAMPLE];
	    hairness_mean += Hairness[j][i*SEE_SAMPLE];
	  } else {
	    cusum_mean += Cusum[j][i*SEE_SAMPLE-3];
	    hairness_mean += Hairness[j][i*SEE_SAMPLE-3];
	  }
 	  meanM += maxIndiv[j][i];
	}

	cusum_mean /= nr_runs;
	hairness_mean /= nr_runs;
	meanA /= nr_runs;
	meanV /= nr_runs;
	meanN /= nr_runs;
	meanD /= nr_runs;
	meanM /= nr_runs;

	for(unsigned long j = 0; j < nr_runs; j++){
	  std_devA += pow(myAcceptanceT[j][i] - meanA,2);
	  std_devV += pow(myAverage[j][i] - meanV,2);
	  std_devN += pow(nrHighIndiv[j][i] - meanN,2);
	  std_devD += pow(distHighIndiv[j][i] - meanD,2);
	  if(i != sampleSize_Const/SEE_SAMPLE){
	    std_dev_c += pow(cusum_mean - Cusum[j][i*SEE_SAMPLE],2);
	    std_dev_h += pow(hairness_mean - Hairness[j][i*SEE_SAMPLE],2);
	  } else {
	    std_dev_c += pow(cusum_mean - Cusum[j][i*SEE_SAMPLE-3],2);
	    std_dev_h += pow(hairness_mean - Hairness[j][i*SEE_SAMPLE-3],2);
	  }
	  std_devM += pow(maxIndiv[j][i] - meanM,2);      
	}      
	
	std_devA = sqrt(std_devA/nr_runs);
	std_devV = sqrt(std_devV/nr_runs);
	std_devN = sqrt(std_devN/nr_runs);
	std_devD = sqrt(std_devD/nr_runs);
	std_dev_c = sqrt(std_dev_c/nr_runs);
	std_dev_h = sqrt(std_dev_h/nr_runs);
	std_devM = sqrt(std_devM/nr_runs);
	
	cumA += meanA;
	cumV += meanV;
	cumN += meanN;
	
	myAcceptance << i*SEE_SAMPLE << "\t" << meanA << "\t" << std_devA << "\t" << cumA/(i+1.) << "\t"<< meanV << "\t" << std_devV << "\t" << cumV/(i+1.)  		    
		     << "\t" << meanN << "\t" << std_devN << "\t" << cumN/(i+1.) << "\t" << meanM << "\t" << std_devM << "\t" << meanD << "\t" << std_devD 
		     << "\t" << cusum_mean << "\t" << std_dev_c  << "\t" << hairness_mean << "\t" << std_dev_h 
		     << "\n";

	//meanAv += meanA;
	//std_meanAv += pow(meanA,2.0)/nr_runs;

      }

      for(unsigned long j = 0; j < nr_runs; j++){
	  meanAv += acceptance_ratio[j]/generation;
	  std_meanAv += pow(acceptance_ratio[j]/generation,2.0);
      }
      
      meanAv /= nr_runs;
      std_meanAv = std_meanAv/nr_runs - pow(meanAv,2.0);
      myAcceptance << meanAv << "\t" << std_meanAv << "\n";
#else 
	double meanA = 0, std_devA = 0;
	double meanV = 0, std_devV = 0;
	double meanN = 0, std_devN = 0;
	double cusum_mean = 0, hairness_mean = 0;
	double std_dev_c = 0, std_dev_h = 0;
	double meanD = 0, std_devD = 0;
	double meanM = 0, std_devM = 0;

       for(unsigned long i = 1; i < (sampleSize_Const/SEE_SAMPLE) + 1; i++){
	  meanA += myAcceptanceT[i];
	  meanV += myAverage[i];
	  meanN += nrHighIndiv[i];
	  meanD += distHighIndiv[i];
	  if(i !=sampleSize_Const/SEE_SAMPLE){
	    cusum_mean += Cusum[i*SEE_SAMPLE];
	    hairness_mean += Hairness[i*SEE_SAMPLE];
	  } else {
	    cusum_mean += Cusum[i*SEE_SAMPLE-3];
	    hairness_mean += Hairness[i*SEE_SAMPLE-3];
	  }
 	  meanM += maxIndiv[i];
       }

	cusum_mean /= (sampleSize_Const/SEE_SAMPLE) + 1;
	hairness_mean /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanA /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanV /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanN /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanD /= (sampleSize_Const/SEE_SAMPLE) + 1;
	meanM /= (sampleSize_Const/SEE_SAMPLE) + 1;
	
	for(unsigned long j = 0; j < (sampleSize_Const/SEE_SAMPLE) + 1; j++){
	  std_devA += pow(myAcceptanceT[j] - meanA,2);
	  std_devV += pow(myAverage[j] - meanV,2);
	  std_devN += pow(nrHighIndiv[j] - meanN,2);
	  std_devD += pow(distHighIndiv[j] - meanD,2);
	  if(j != sampleSize_Const/SEE_SAMPLE){
	    std_dev_c += pow(cusum_mean - Cusum[j*SEE_SAMPLE],2);
	    std_dev_h += pow(hairness_mean - Hairness[j*SEE_SAMPLE],2);
	  } else {
	    std_dev_c += pow(cusum_mean - Cusum[j*SEE_SAMPLE-3],2);
	    std_dev_h += pow(hairness_mean - Hairness[j*SEE_SAMPLE-3],2);
	  }
	  std_devM += pow(maxIndiv[j] - meanM,2);      
	}      
	
	std_devA = sqrt(std_devA/(sampleSize_Const/SEE_SAMPLE + 1));
	std_devV = sqrt(std_devV/(sampleSize_Const/SEE_SAMPLE + 1));
	std_devN = sqrt(std_devN/(sampleSize_Const/SEE_SAMPLE + 1));
	std_devD = sqrt(std_devD/(sampleSize_Const/SEE_SAMPLE + 1));
	std_dev_c = sqrt(std_dev_c/(sampleSize_Const/SEE_SAMPLE + 1));
	std_dev_h = sqrt(std_dev_h/(sampleSize_Const/SEE_SAMPLE + 1));
	std_devM = sqrt(std_devM/(sampleSize_Const/SEE_SAMPLE + 1));
	
	cumA += meanA;
	cumV += meanV;
	cumN += meanN;
	
	myAcceptance << sampleSize_Const*SEE_SAMPLE << "\t" << meanA << "\t" << std_devA << "\t" << meanA * ((sampleSize_Const/SEE_SAMPLE) + 1) << "\t" 
		     << meanV << "\t" << std_devV << "\t" << meanV * ((sampleSize_Const/SEE_SAMPLE) + 1)   		    
		     << "\t" << meanN << "\t" << std_devN << "\t" << meanN*((sampleSize_Const/SEE_SAMPLE) + 1) << "\t" 
		     << meanM << "\t" << std_devM << "\t" << meanD << "\t" << std_devD 
		     << "\t" << cusum_mean << "\t" << std_dev_c  << "\t" << hairness_mean << "\t" << std_dev_h 
		     << "\n";
#endif //nr_runs

  myAcceptance.flush();
  myAcceptance.close();
#endif //aceptance ratio

#ifndef MIXTURE_REALS
#ifdef DISTANCE
  cout << "average on this distance \n";
  double std_dev[sizeGENOME+1];
  for(unsigned long i = 0; i < sizeGENOME+1; i++){
      tempDistance[i] = 0;
      std_dev[i] = 0;

      for(unsigned long j = 0; j < nr_runs; j++)
	  tempDistance[i] += myDistance[j][i];

      tempDistance[i] /= nr_runs;

      for(unsigned long j = 0; j < nr_runs; j++)
	  std_dev[i] += pow(myDistance[j][i] - tempDistance[i],2);

      std_dev[i] = sqrt(std_dev[i]/nr_runs);

      myFileDistance << i << "\t" << tempDistance[i] << "\t" <<  std_dev[i] << "\n";
  }
  
  myFileDistance.flush();
  myFileDistance.close();
#endif
#endif //mixture reals

#ifndef MULTIPLE_RUNS
#ifdef HISTOGRAM
  myTempFile.flush();
  myTempFile.close();

  myFileFG.flush();
  myFileFG.close();

  myFileFit.flush();
  myFileFit.close();

  myFileHash.flush();
  myFileHash.close();

  ofstream myFileD("err_recomb_dist");		
  if(!myFileD.is_open()){
      cout << "Error opening file for data " << fileData << "\n";
    return;
  }

  ifstream myTempFile1("temporar.txt");		
  if(!myTempFile1.is_open()){
      cout << "Error opening file for data " << fileData << "\n";
    return;
  }

  double tempV[sampleSize+2][nr_runs+1];
  double tempD[sampleSize+2][nr_runs+1];
  double tempV1[sampleSize+2][nr_runs+1];
  double tempD1[sampleSize+2][nr_runs+1];

  for(unsigned long i = 0; i < nr_runs; i++)
      for(unsigned long j = 0; j < sampleSize+1;j++){
	  myTempFile1 >> tempV[j][i];
	  myTempFile1 >> tempV1[j][i];
	  myTempFile1 >> tempD[j][i];
	  myTempFile1 >> tempD1[j][i];
      }

  for(unsigned long j = 0; j < sampleSize+1;j++){
      for(unsigned long i = 0; i < nr_runs; i++){
	  myFile << tempV[j][i] << "\t";
	  myFile << tempV1[j][i] << "\t";
	  myFileD << tempD[j][i] << "\t";
	  myFileD << tempD1[j][i] << "\t";
      }
      myFile << "\n";
      myFileD << "\n";
  }
  myTempFile1.close();

//#ifdef DEBUG_
  myFile.flush();
  myFile.close();
//#endif
  myFileD.flush();
  myFileD.close();
#else //histogram
  for(int i = 0; i < records; i++){
    myFileAllSplit[i].flush();
    myFileAllSplit[i].close();

    myFile1ChainSplit[i].flush();
    myFile1ChainSplit[i].close();
  }
  myFileAllNoSplit.flush();
  myFileAllNoSplit.close();

  myFile1ChainNoSplit.flush();
  myFile1ChainNoSplit.close();

  delete[] myFileAllSplit;
  delete[] myFile1ChainSplit;

  myFileHash.flush();
  myFileHash.close();
#endif //histogram
#endif //MULTIPLE_RUNS

#else //not trace runs

  myFile.flush();
  myFile.close();

  myFileI.flush();
  myFileI.close();

  myFileII.flush();
  myFileII.close();

  myFileIII.flush();
  myFileIII.close();

  myFileIV.flush();
  myFileIV.close();

#endif // TRACE_RUN
  delete[] table; 
}

#ifndef RANDOM_WALK
unsigned long binaryPopulation::acceptancePopulation(unsigned long* table, unsigned long sizeCoup){
    double prop = 1;
//    UniformDistibution genrand_real1(0,1);

#ifdef DEBUG_
    double* temp1 = new double[sizeGENOME];
    for(unsigned long i = 0; i<sizeCoup; i++) {
	MCMCstring[table[i]] -> getGENOME(temp1);
	cout << "parent " << i << "\t";
	for(unsigned long j = 0; j < sizeGENOME; j++)
	    cout << temp1[j] ;
	cout << "\t" << MCMCstring[table[i]] -> getValue() << "\t table" << table[i] << "\n";
	proposedMCMCstring[table[i]] -> getGENOME(temp1);
	cout << "child " << i << "\t";
	for(unsigned long j = 0; j < sizeGENOME; j++)
	    cout << temp1[j];
	cout << "\t" << proposedMCMCstring[table[i]] -> getValue() << "\t table" 
	     << table[i]<< "\n";
    }  
    delete[] temp1;
#endif //DEBUG_
    
    //#ifndef PCTX
    //for(unsigned long i = 0; i<1; i++) {
    //#else
#ifndef SUM_DISTRIB
    for(unsigned long i = 0; i<sizeCoup; i++) {
      //#endif 

#ifdef Boltzmann
#ifdef Normal_Boltzmann
	double child = proposedMCMCstring[table[i]]->getValue();
	double parent = MCMCstring[table[i]]->getValue();
	//prop *= exp((child - parent)/MCMCstring[table[i]]->getTEMPERATURE());
	prop *= exp((double)(child - parent)/(double)DEFAULT_TEMPERATURE);
#else 
	prop = prop * 
	    exp((log(MCMCstring[table[(i+1)%sizeCoup]]->getValue()) - log(MCMCstring[table[i]]->getValue()))
		/MCMCstring[table[i]]->getTEMPERATURE());
#endif //NORMAL_BOLTZMAN
#else //boltzman
	
#ifndef MIXTURE_REALS
	double child = floor(proposedMCMCstring[table[i]]->getValue());
	double parent = floor(MCMCstring[table[i]]->getValue());
	if(floor(child) == 0)
	    child = AwayFromZero;
	if(floor(parent) == 0)
	    parent = AwayFromZero;
#else
	double child = proposedMCMCstring[table[i]]->getValue();
	double parent = MCMCstring[table[i]]->getValue();
#endif //FLOOR

	//prop *= pow(child / parent, 1/DEFAULT_TEMPERATURE);
	prop *= child / parent;
	
#ifdef DEBUG_
	cout << "prop (" << i << "," << child << "," << parent << ") = \t" << prop << "\n";
#endif
	//prop *= pow(MCMCstring[table[(i+1)%sizeCoupling]]->getValue() / 
	//	MCMCstring[table[i]]->getValue(),1.0/MCMCstring[table[i]]->getTEMPERATURE());
#endif //non - Boltzman

    }
	
#else //sum distrib
    prop = 0;
    for(unsigned long i = 0; i<sizeCoup; i++)
      prop+= MCMCstring[table[i]]->getValue();

    double prop_child = 0;
    for(unsigned long i = 0; i<sizeCoup; i++)
      prop_child += proposedMCMCstring[table[i]]->getValue();

    prop = prop_child/prop;

    
#endif //sum distrib

#ifndef METROPOLIS_ALGORITHM
#ifndef POPULATION_MCMC
	//cout << " proposed distrib (" << proposedProposalDistribution[table[0]] << "," << proposedProposalDistribution[table[1]] << "," << proposedProposalDistribution[table[2]]<<
	//") , proposal (" <<proposalDistribution[table[0]] << "," <<proposalDistribution[table[1]] << ","<<proposalDistribution[table[2]]<< "\n"; 
    if(MCMCstring[table[0]] -> mix_rotation == 0){
      prop *= (double)proposedProposalDistribution[table[0]]/(double)proposalDistribution[table[0]];
    }
#endif // METROPOLIS_ALGORITHM 
#endif  

    if(prop >= 1) return 1;
    double temp = genrand_real1();
    if(temp > prop) {
      //cout << "Rejected at generation " << generation << " in run " << runs-1 << "\n"; 
	return 0;
    }
    return 1;
}

#endif //RANDOM_WALK

void binaryPopulation::acceptanceProposed(unsigned long** table, unsigned long dimension){
    unsigned long tempSize = dimension;
#ifndef MIXTURE_REALS
    unsigned long temp[sizeGENOME];
#else
#ifndef NESTED_MH
    double temp[sizeGENOME];
#else
    double*** temp = new double**[sizeCoupling];
    for(unsigned long i = 0; i < sizeCoupling; i++){
      temp[i] = new double*[sizeCoupling];
      for(unsigned long j = 0; j < sizeCoupling; j++)
	temp[i][j] = new double[sizeGENOME];
    }
#endif //nested mh
#endif
    
//#ifdef METROPOLIS_HASTINGS_ALGORITHM
    //  cout << "proposed distrbution";
    // proposalDistributionCalculation(table, dimension);
//#endif

    for(unsigned long i = 0; i<tempSize ; i++){
      //for each pair
      //cout << "pair = " << i << "\n";
 
      //      for(unsigned long j = 0; j<sizeCoupling; j++){
      //proposedMCMCstring[table[i][j]]->setACCEPTED(0); 
      //	MCMCstring[table[i][j]]->setACCEPTED(0);
      //}
      
#ifndef NESTED_MH
#ifdef CYCLE
      if(generation % 2 == 0){
#else
#ifdef MIXTURE
/*#ifdef SINGLE_INDIV_ROT
      if(COUPLING_PROBABILITY == 0 || mix_mutation == 1 || MCMCstring[table[i][0]]->mix_rotation != 0){
      #else*/
      if(COUPLING_PROBABILITY == 0 || mix_mutation == 1){
//#endif
	//cout << "mixture = " << mix_mutation << " and coupling" << COUPLING_PROBABILITY << " binaryPopulation_trap.cpp:6389 \n"; 
#else
      if(COUPLING_PROBABILITY == 0){
	//cout << "not coupled \t";
#endif //MIXTURE
#endif //cycle

/*#ifdef SINGLE_INDIV_ROT
	  int acept_number = 0;
	  if(MCMCstring[table[i][0]] -> accepted != 0)
	  #endif*/
//#define OneMut
#ifdef OneMut
	for(unsigned long j = 0; j < 1; j++){
#else
	for(unsigned long j = 0; j < sizeCoupling; j++){
#endif
	  proposedMCMCstring[table[i][j]]->getGENOME(temp);
	  double tempValue = proposedMCMCstring[table[i][j]]->getValue();
	  //cout << " i" << i << " j " << j << " (" << temp[0] << " , "<< temp[1] << " , " <<temp[2] << ") accept " << MCMCstring[table[i][j]] -> accepted << " , " << proposedMCMCstring[table[i][j]] -> accepted << "\n";
	  /*#ifndef METROPOLIS_ALGORITHM
#ifdef METROPOLIS_HASTINGS_ALGORITHM
	  double tempProposalDistribution = MCMCstring[table[i][j]]->getProposalDistribution();
#else
	  double tempProposalDistribution = proposedMCMCstring[table[i][j]]->getProposalDistribution();     
#endif
	
	  
#ifndef MIXTURE_REALS
	  if(MCMCstring[table[i][j]] -> acceptance(temp,tempValue,tempProposalDistribution) == 1){
#else
	  if(MCMCstring[table[i][j]] -> acceptance(temp,tempValue) == 1){
#endif
	    //if(MCMCstring[table[i][0]] -> accepted != 0)
	    // MCMCstring[table[i][j]]->acceted = 1;
	    //else MCMCstring[table[i][j]]->accepted = 0;
	    //cout << "accepted old (" << MCMCstring[table[i][j]] -> genome[0] << ", " 
	    //<< MCMCstring[table[i][j]] -> genome[1] << "," << MCMCstring[table[i][j]]->value <<")\t";
	    MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	    cout << "accepted new (" << MCMCstring[table[i][j]] -> genome[0] << ", " << 
	      MCMCstring[table[i][j]] -> genome[1] << "," <<  MCMCstring[table[i][j]] -> genome[2] << "," << MCMCstring[table[i][j]]->value <<")\n";
	    acept_number++;
	  }
	  else {
	    MCMCstring[table[i][j]] -> nextGenerationBlank();
	    cout << "rejected old (" << MCMCstring[table[i][j]] -> genome[0] << ", " 
		 << MCMCstring[table[i][j]] -> genome[1] << "," <<  MCMCstring[table[i][j]] -> genome[2] << "," << MCMCstring[table[i][j]]->value <<")\t";
	    MCMCstring[table[i][j]]->accepted = 0; 
	  }
	  #else*/
#ifndef RANDOM_WALK
//#ifndef SINGLE_INDIV_ROT
	  if(proposedMCMCstring[table[i][j]] -> accepted != 0){
//#endif
	    if(MCMCstring[table[i][j]] -> acceptance(temp,tempValue) == 1){
#endif //RANDOM_WALK
	      //cout << "accepted single \n";
	      MCMCstring[table[i][j]]->accepted = 1;
	      //cout<<"accepted old ("<<MCMCstring[table[i][j]] -> genome[0] << ", " << MCMCstring[table[i][j]] -> genome[1] << "," << MCMCstring[table[i][j]]->value <<")\t";
	      MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	      //cout<<"accepted new ("<<MCMCstring[table[i][j]] -> genome[0] << ", " << MCMCstring[table[i][j]] -> genome[1] << "," << MCMCstring[table[i][j]]->value <<")\n";
	    }
#ifndef RANDOM_WALK
	    else {
	      //cout << "not accepted single\n";
	    MCMCstring[table[i][j]] -> nextGenerationBlank();
	    MCMCstring[table[i][j]]->accepted = 0; 
	    }
#endif //RANDOM_WALK
//#ifndef SINGLE_INDIV_ROT
	  }else {
	    //cout << "mut not accepted \n";
	    MCMCstring[table[i][j]] -> nextGenerationBlank();
	    MCMCstring[table[i][j]]->accepted = 0;
	    //MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	  }
	  //#endif //METROPOLIS
//#endif
/*#ifdef SINGLE_INDIV_ROT
	if(acept_number > MCMCstring[table[i][0]] -> accepted){
	  MCMCstring[table[i][0]] -> accepted = acept_number;
	  for(int j = 1; j < sizeCoupling; j++)
	    MCMCstring[table[i][j]] -> accepted = 0;
	} else{
	  for(int j = 1; j < sizeCoupling; j++)
	    MCMCstring[table[i][j]] -> accepted = 0;
	}
	#endif*/
	  }
      } else {
	//#ifdef COUPLED_ACCEPTANCE
	//cout << " i" << i << " (" << temp[0] << " , "<< temp[1] << " , " <<temp[2] << ") accept " << MCMCstring[table[i][0]] -> accepted << " , " << proposedMCMCstring[table[i][0]] -> accepted << "\n";
#ifdef HALF_PCTX
#ifndef RANDOM_WALK
	if(MCMCstring[table[i][0]] -> acceptance(temp,tempValue) == 1){
#endif //RANDOM_WALK
	  //cout << "accepted single \n";
	  MCMCstring[table[i][0]]->accepted = 1;
	  //cout<<"accepted old ("<<MCMCstring[table[i][j]] -> genome[0] << ", " << MCMCstring[table[i][j]] -> genome[1] << "," << MCMCstring[table[i][j]]->value <<")\t";
	  MCMCstring[table[i][0]] -> nextGeneration(temp,tempValue,sizeGENOME);
	  //cout<<"accepted new ("<<MCMCstring[table[i][j]] -> genome[0] << ", " << MCMCstring[table[i][j]] -> genome[1] << "," << MCMCstring[table[i][j]]->value <<")\n";
	}
#ifndef RANDOM_WALK
	else {
	  //cout << "not accepted single\n";
	    MCMCstring[table[i][0]] -> nextGenerationBlank();
	    MCMCstring[table[i][0]]->accepted = 0; 
	}
#endif //RANDOM_WALK
#ifndef RANDOM_WALK
	if(acceptancePopulation((table[i]+1),sizeCoupling -1) == 1){
	    // cout << "accepted coupled\n";
#endif
	    for(unsigned long j = 1; j<sizeCoupling; j++){
	      proposedMCMCstring[table[i][j]]->getGENOME(temp);
	      double tempValue = proposedMCMCstring[table[i][j]]->fitness();
	      
	      //cout << " accept (" << i << "," << table[i][0] << "," << MCMCstring[table[i][0]] -> accepted 
	      //   << " , " << MCMCstring[table[i][1]] -> accepted << " , " << MCMCstring[table[i][2]] -> accepted 
	      //   << " , " << proposedMCMCstring[table[i][0]] -> accepted << " , " 
	      ///   << proposedMCMCstring[table[i][1]] -> accepted << " , " 
	      //  << proposedMCMCstring[table[i][2]] -> accepted<< ")\n";
	      MCMCstring[table[i][j]]->accepted = 1;
	      MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	      //cout << "set YES " << temp[0] << "," << temp[1] << "," << tempValue << "\n";
	    }
	  } else {
	      //cout << "set NO (" << i << "," << temp[0] << "," << temp[1] << ")\n"; 
	    for(unsigned long j = 0; j<sizeCoupling; j++){
	      MCMCstring[table[i][j]]->accepted = 0;
	      MCMCstring[table[i][j]] -> nextGenerationBlank();
	    }
	  }
	    //}
	    //MCMCstring[table[i][j]]->setACCEPTED(1);
	    //}
      }
#else //half_PCTX
#ifndef RANDOM_WALK
#ifdef IsRecombination
	if(MCMCstring[table[i][0]] -> accepted != 0){
	  //cout << "candidate for acceptance \n";
#endif 
	  if(acceptancePopulation(table[i],sizeCoupling) == 1){
	    // cout << "accepted coupled\n";
#endif
	    for(unsigned long j = 0; j<sizeCoupling; j++){
	      proposedMCMCstring[table[i][j]]->getGENOME(temp);
	      double tempValue = proposedMCMCstring[table[i][j]]->fitness();
	      
	      //cout << " accept (" << i << "," << table[i][0] << "," << MCMCstring[table[i][0]] -> accepted 
	      //   << " , " << MCMCstring[table[i][1]] -> accepted << " , " << MCMCstring[table[i][2]] -> accepted 
	      //   << " , " << proposedMCMCstring[table[i][0]] -> accepted << " , " 
	      ///   << proposedMCMCstring[table[i][1]] -> accepted << " , " 
	      //  << proposedMCMCstring[table[i][2]] -> accepted<< ")\n";
	      MCMCstring[table[i][j]]->accepted = 1;
	      MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	      //cout << "set YES " << temp[0] << "," << temp[1] << "," << tempValue << "\n";
	    }
	  } else {
	      //cout << "set NO (" << i << "," << temp[0] << "," << temp[1] << ")\n"; 
	    for(unsigned long j = 0; j<sizeCoupling; j++){
	      MCMCstring[table[i][j]]->accepted = 0;
	      MCMCstring[table[i][j]] -> nextGenerationBlank();
	    }
	  }
	    //}
	    //MCMCstring[table[i][j]]->setACCEPTED(1);
	    //}
#ifndef RANDOM_WALK
#ifdef IsRecombination
	}else{
	  //cout << "not accepted coupled" << i << " \n";
	      for(unsigned long j = 0; j<sizeCoupling; j++){
		MCMCstring[table[i][j]] -> nextGenerationBlank();
		MCMCstring[table[i][j]]->accepted = 0;
	      }
	}
#endif
#endif //RANDOM_WALK
	    //#endif //COUPLIED ACCEPTANCE
      }
#endif //half PCTX
#else //ifndef nested MH
	//.1. first each individual compets agains one parent
	double tempValue[sizeCoupling];
	unsigned long tempAccepted[sizeCoupling];

	//to remember: 0 parents; 1 children in the temp data-structure
	//1.0 first read de parents in temp
	unsigned long diff = 0;
	for(unsigned long j = 0; j<sizeCoupling; j++){
	  proposedMCMCstring[table[i][j]]->getGENOME(temp[1][j]);
	  tempValue[j] = proposedMCMCstring[table[i][j]]->getValue();
	  if(MCMCstring[table[i][j]] -> acceptance(temp[1][j],tempValue[j]) == 1){
	    tempAccepted[j] = 1;
	    diff++;
	  } else 
	    tempAccepted[j] = 0;
	}
	
	//2. daca amindoi individu aceptati sau rejectati acceptai si rejecteaza si de data asta
	if(diff == 0 || diff == sizeCoupling){
	  for(unsigned long j = 0; j<sizeCoupling; j++){
	    if(diff == sizeCoupling){
	      MCMCstring[table[i][j]] -> nextGeneration(temp[1][j],tempValue[j],sizeGENOME);
	      MCMCstring[table[i][j]]->setACCEPTED(1);
	    } else {
	      MCMCstring[table[i][j]] -> nextGenerationBlank();
	      MCMCstring[table[i][j]]->setACCEPTED(0);
	    }
	  }	  
	} else {
	  //.3. daca un individ aceptat rejectat atunci trebuie sa ii acepti si cu proposal dsitribution
	  for(unsigned long j = 0; j < sizeCoupling; j++)
	    MCMCstring[table[i][j]] -> getGENOME(temp[0][j]);
	  double prop = proposalDistributionCalculation(temp,dimension);
	  //accept
	  if(prop >= 1 || genrand_real1() > prop){
	    for(unsigned long j = 0; j < sizeCoupling; j++){
	      if(tempAccepted[j] == 1){
		MCMCstring[table[i][j]] -> nextGeneration(temp[1][j],tempValue[j],sizeGENOME);
		MCMCstring[table[i][j]]->setACCEPTED(1);
	      } else {
		MCMCstring[table[i][j]] -> nextGenerationBlank();
		MCMCstring[table[i][j]]->setACCEPTED(0);
	      }
	    }
	  } else {
	    //reject everithing
	    for(unsigned long j = 0; j < sizeCoupling; j++){
	      MCMCstring[table[i][j]] -> nextGenerationBlank();
	      MCMCstring[table[i][j]]->setACCEPTED(0);
	    }
	  }
	}
#endif //nested MH

	
#ifdef SWAP_CHAIN
	if(SWAP > genrand_int32()){
	    if(acceptancePopulation(table[i]) == 1){
		unsigned long tempTable[sizeCoupling][sizeGENOME];
		
		for(unsigned long j = 0; j<sizeCoupling; j++)
		    MCMCstring[table[i][(j+1)%sizeCoupling]]->getGENOME(tempTable[j]);
		
		for(unsigned long j = 0; j<sizeCoupling; j++)
		    MCMCstring[table[i][j]] -> setGENOME(tempTable[j]);
		
		//cout << "I accepted sometihng for " << i <<"\n";
		continue;
	    }
	}
	//   cout << "I reject for "<< i <<"\n";
#endif
	}

}


#else
#ifdef ELITIST_ACCEPTANCE
// --------- Elitist Evolution ---------------
//
//myGenerations - number of generations to run
//
void binaryPopulation::evolution(){
    unsigned long tempSize = sizeMCMCstring / sizeCoupling;
    unsigned long** table = new unsigned long*[tempSize + 1];
    for(unsigned long i = 0; i< tempSize + 1; i++)
	table[i] = new unsigned long[sizeCoupling];
    
    cout << "ELITIST_ACCEPTANCE \n";
    
#ifdef DEBUG_
    ofstream myFile(fileData);		
    if(!myFile.is_open())
	cout << "Error opening file";
#endif
    
#ifdef EXPAND_KULLBACK_INFORMATION
    expandTrueUnivariateDistribution();
#else
    readHistogram();
#endif
    
    while(runs < nr_runs){
	runs++;
	setSeed();
	
	while(generation <=  burn_in + sampleSize*interSamplesize)
	{		
	    assignTemperature();
	    
	    unsigned long dimension = sizeMCMCstring/sizeCoupling;
#ifndef NEAR_CHANGING  
	    randomCOUPLING(table);
#else
	    dimension = neighborsCOUPLING_NO_OVERLAPING(table);
#endif
	    proposedGeneration(table,dimension);
	    acceptanceProposed(table,dimension);
	    postProcessing();
	    
	    simul = insert_simulation(simul,generation,procentRightSolutions);
	    //show_simulations(simul);
	    write_to_file_all(simul_all_proc_right,simul,NULL);
	    
#ifdef KULLBACK_INFORMATION
	    mixingTimeUnivariateCalculation();      
#endif //KULLBACK_INFORMATION
	    
	    generation++;
	}	
	
#ifdef DEBUG_
	printLEVEL(myFile);
#else 
	scoreDiversity();
#endif
	
	globalRightSolutions[runs-1] = procentRightSolutions;
	
#ifdef CUSUM
	CusumCalculation();
#endif
	
#ifdef GELMAN_RUBIN
	varianceHigerMomentsUnivariateCalculation(GELMAN_Moment,sizeMCMCstring - 1);
	varianceMultivariateCalculation();
#endif
	
	reset();
    }
    
#ifdef GELMAN_RUBIN
    printReductionGNUFILE(GELMAN_Moment);
    printReductionMultivariateGNUFILE();
#endif
    
    printDiversityGNUFILE();
    printDiversityDetailGNUFILE();  
    printDiversityRightSolutions();
    
#ifdef KULLBACK_INFORMATION
    printMixingGNUFILE();
#endif
    
#ifdef CUSUM
    printCusumGNUFILE();
#endif
    
#ifdef DEBUG_
    myFile.close();
#endif
    
#ifndef FILE_BOOKEEPING
    diversityFile.flush();
    diversityFile.close();
#endif
    
    delete[] table;
}


//
//tempTable - temporal table with the values acceptance for (x',y) and (x',y) and (x',y') and (x,y)
//prop - store the result of acceptance operation
//tempStringValue - table with the fitness functiob for the involved chains
//tempTEMPvalue - table with the temperature for the involved chains
//size - size of coupling
unsigned long binaryPopulation::acceptancePopulation(unsigned long* table){
    double prop = 1;
    double tempTable[sizeCoupling];
    for(unsigned long i = 0; i<sizeCoupling; i++) {
	proposedMCMCstring[table[i]]->setACCEPTED(0);
	tempTable[i] = exp( (proposedMCMCstring[table[i]]->getValue() - MCMCstring[table[i]]->getValue()) 
			    / MCMCstring[table[i]] -> getTEMPERATURE());
	if(tempTable[i] > 1){
	  prop = prop * tempTable[i]; 
	  proposedMCMCstring[table[i]]->setACCEPTED(1);
	}
#ifdef DEBUG_
	cout << " tempTable[" <<i<<"] = "<<tempTable[i];
	cout << "value i = " << (MCMCstring[table[i]]) ->getValue()<<" tempValue i "<<tempStringValue[i];		
	cout << "temperaturei " << MCMCstring[table[i]] -> getTEMPERATURE() << " tempTemperature " << tempTEMPvalue[i] << "\t";
	cout << "prop " << prop << "\n"; 
#endif
}
	if(prop > 1) return 1;
	
	//else pick up some combination given some algorithm
	proportionalSelection(tempTable,table);
	return 1;
}
	
	
void binaryPopulation::proportionalSelection(double* tempResults, unsigned long* table)
{
  double tempTable[(unsigned long)pow(2,sizeCoupling)-1];
  for(unsigned long i = 0; i < sizeCoupling; i++) 
    tempTable[i] = 1;
  
  unsigned long index = 0;
  double norm = 0;

  //compute tempTable
  for(unsigned long i = 0; i < sizeCoupling; i++) {
    //compute C_{i}^{sizeCoupling}
    if(i == 0) {
      for(unsigned long k = 0; k<sizeCoupling;k++)
	tempTable[index] *= tempResults[k];
      norm += tempResults[index];
    }
    else 
      if(i == 1) {
	for(unsigned long j =0; j < (unsigned long)combinatorial(sizeCoupling,i); j++){
	  tempResults[index+j] = tempResults[j];
	  norm += tempTable[index + j];
	}
      }
      else {
	//recursively resolve the problem
	//for(unsigned long j =0; j < (unsigned long)combinatorial(sizeCoupling,i); j++){
	  //tempTable[i+j] = tempResults[j];
	//}
      }
    index++;
  }

  //find out the maximum and multiply by it
  double max = -1;
  index = 0;
  for(unsigned long i = 0; i < (unsigned long)pow(2,sizeCoupling)-1; i++)
    if(max < tempTable[i]) {
      max = tempTable[i];
      index = i;
    }

  //normalize the results from tempTable
  for(unsigned long i = 0; i < (unsigned long)pow(2,sizeCoupling)-1; i++)
    tempTable[i] *= max / norm; 

  //print 
#ifdef DEBUG_
  cout << " \n value Table for rejection \n";
  for(unsigned long i = 0; i < (unsigned long)pow(2,sizeCoupling)-1; i++)
    cout << tempTable[i] << " || "; 
    cout << " norm = "<< norm << " maxim = " << max;
#endif

  double temp = genrand_real1();
#ifdef DEBUG_  
  cout << " \n    With ramdom number = "<<temp <<"combinari (1,"<<sizeCoupling<<") = "<<combinatorial(sizeCoupling,1) <<"\n";
#endif

  double leftNorm = 0;
  for(unsigned long i = 0; i < (unsigned long)pow(2,sizeCoupling)-1; i++) {
		if(temp > leftNorm + tempTable[i]) {
		  leftNorm += tempTable[i]; 
		  continue;
		}
		else {
		  //fill in the results vector
		  if(i == 0 ) {
		    for(unsigned long k = 0; k < sizeCoupling ; k++)
		      proposedMCMCstring[table[k]]->setACCEPTED(1);
		    return;
		  }
		  if(i - 1 < combinatorial(sizeCoupling,1)) {
		    proposedMCMCstring[table[i - 1]] ->setACCEPTED(1);
		    return;
		  } else {
		    //recursively find out the place, otherwise
		  }
		}
	}
  return;
}

//
//sizeCoupling - number of chains to be coupled randomly
//tempSize - number of couples formed in population
//table - the table with the individuals that couples
//tempString - contained for the proposed individuals
//tempAcceptedProposal - table which contains the accepted proposal
//
void binaryPopulation::acceptanceProposed(unsigned long** table, unsigned long dimension){
	unsigned long tempSize = dimension;
	unsigned long temp[sizeGENOME];
	for(unsigned long i = 0; i<tempSize ; i++)
	{
	  for(unsigned long j = 0; j<sizeCoupling; j++){
	    proposedMCMCstring[table[i][j]]->setACCEPTED(0); 
	    MCMCstring[table[i][j]]->setACCEPTED(0);
	  }
	  unsigned long index = acceptancePopulation(table[i]);
		
#ifdef DEBUG_
	  for(unsigned long j = 0; j<sizeCoupling; j++)
	    cout <<proposedMCMCstring[table[i][j]]->getACCEPTED() << "||";
	  cout << " table pf acceptance \n";
#endif
	  for(unsigned long j = 0; j<sizeCoupling; j++)
	    if(proposedMCMCstring[table[i][j]]->getACCEPTED() == 1){
	      proposedMCMCstring[table[i][j]]->getGENOME(temp);
	      double tempValue = proposedMCMCstring[table[i][j]]->getValue();
	      MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	    }
	    else MCMCstring[table[i][j]] -> nextGenerationBlank();
	  
	}			
}
#else

#ifdef ELITIST_ACCEPTANCE_WITH_REGENERATION
// ------ Elitisit coupling with regeneration -------------
//
//sizeCoupling - number of chains to be coupled randomly
//tempSize - number of couples formed in population
//table - the table with the individuals that couples
//tempString - contained for the proposed individuals
//tempAcceptedProposal - table which contains the accepted proposal
//

void binaryPopulation::evolution(){
  unsigned long tempSize = sizeMCMCstring / sizeCoupling;
  unsigned long** table = new unsigned long*[tempSize + 1];
  for(unsigned long i = 0; i< tempSize + 1; i++)
    table[i] = new unsigned long[sizeCoupling];
  
  cout << "Elitist acceptance with regeneration \n";

#ifdef DEBUG_
  ofstream myFile(fileData);		
  if(!myFile.is_open())
    cout << "Error opening file";
#endif

#ifdef EXPAND_KULLBACK_INFORMATION
    expandTrueUnivariateDistribution();
#else
    readHistogram();
#endif

  while(runs < nr_runs){
    runs++;
    setSeed();
    
    while(generation <=  burn_in + sampleSize*interSamplesize)
      {	
	assignTemperature();
	unsigned long dimension = sizeMCMCstring/sizeCoupling;
#ifndef NEAR_CHANGING  
	randomCOUPLING(table);
#else
	dimension = neighborsCOUPLING_NO_OVERLAPING(table);
#endif
	proposedGeneration(table,dimension);
	acceptanceProposed(table,dimension);
	postProcessing();

	simul = insert_simulation(simul,generation,procentRightSolutions);
	//show_simulations(simul);
	write_to_file_all(simul_all_proc_right,simul,NULL);
	
#ifdef KULLBACK_INFORMATION
	mixingTimeUnivariateCalculation();      
#endif //KULLBACK_INFORMATION

 	//cout << "generation" << generation++ << "\n";
        generation++;	
      }

#ifdef DEBUG_
    printLEVEL(myFile);
#else 
    scoreDiversity();
#endif

    globalRightSolutions[runs-1] = procentRightSolutions;
#ifdef CUSUM
	CusumCalculation();
#endif

#ifdef GELMAN_RUBIN
  varianceHigerMomentsUnivariateCalculation(GELMAN_Moment,sizeMCMCstring - 1);
  varianceMultivariateCalculation();
#endif

    reset();
  }

#ifdef GELMAN_RUBIN
  printReductionGNUFILE(GELMAN_Moment);
  printReductionMultivariateGNUFILE();
#endif

  printDiversityGNUFILE();
#ifdef TRAP_LIKE_FUNCTION
  printDiversityDetailGNUFILE();
#endif
  printDiversityRightSolutions();

#ifdef KULLBACK_INFORMATION
	printMixingGNUFILE();
#endif

#ifdef CUSUM
	printCusumGNUFILE();
#endif

#ifdef DEBUG_
  myFile.close();
#endif

#ifndef FILE_BOOKEEPING
  diversityFile.flush();
  diversityFile.close();
#endif

   delete[] table;
}


void binaryPopulation::acceptanceProposed(unsigned long** table, unsigned long dimension){
	unsigned long tempSize = dimension;
	unsigned long temp[sizeGENOME]; 

	for(unsigned long i = 0; i<tempSize ; i++)
	{
	  //cout <<" the" << i<<"-th couple is= (";  
	  for(unsigned long j = 0; j<sizeCoupling; j++){
	    proposedMCMCstring[table[i][j]]->setACCEPTED(0); 
	    MCMCstring[table[i][j]]->setACCEPTED(0);
	    //cout << table[i][j] << ",";
	  }
	  if(COUPLING_PROBABILITY > genrand_real1()){
	    //cout << ")\n";
	    unsigned long index = acceptancePopulation(table[i]);
	    
	    for(unsigned long j = 0; j<sizeCoupling; j++){
	      if(MCMCstring[table[i][j]]->getACCEPTED() == 1 && proposedMCMCstring[table[i][j]]->getACCEPTED() == 0){
#ifdef DEBUG_
		cout << "j = " << j << "is further! \n";
#endif
		MCMCstring[table[i][j]]->nextGenerationBlank();
	      }
	      else if(MCMCstring[table[i][j]]->getACCEPTED() == 0 && proposedMCMCstring[table[i][j]]->getACCEPTED() == 1){
#ifdef DEBUG_
		cout << "j = " << j << "is replaced :( \n";
#endif
		proposedMCMCstring[table[i][j]]->getGENOME(temp);
		double tempValue = proposedMCMCstring[table[i][j]]->getValue();
		MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
		MCMCstring[table[i][j]]-> setACCEPTED(1);
	      } 
	      else if(MCMCstring[table[i][j]]->getACCEPTED() == 1 && proposedMCMCstring[table[i][j]]->getACCEPTED() == 1){
		unsigned long k = 0;
		while(k < sizeCoupling && MCMCstring[table[i][k]]->getACCEPTED() == 1) k++;
		if(k == sizeCoupling) cout << "error in acceptanceProposed \n";
		else {
		  MCMCstring[table[i][j]]->nextGenerationBlank();
		  proposedMCMCstring[table[i][j]]->getGENOME(temp);
		  double tempValue = proposedMCMCstring[table[i][j]]->getValue();
		  MCMCstring[table[i][k]] -> nextGeneration(temp,tempValue,sizeGENOME);
		}
#ifdef DEBUG_
		cout << "j = " << j << "is replaced and further :) helped by " << k << "\n";
#endif	      
	      } 
	      else {
#ifdef DEBUG_
		cout << "j = " << j << "is dead :((( \n";
#endif	      	    
	      }
	    }
	  }
	  else {
	    for(unsigned long j = 0; j < sizeCoupling; j++){
	      proposedMCMCstring[table[i][j]]->getGENOME(temp);
	      double tempValue = proposedMCMCstring[table[i][j]]->getValue();
#ifdef METROPOLIS_ALGORITHM
	      double tempProposalDistribution = MCMCstring[table[i][j]]->getProposalDistribution();
#else
	      double tempProposalDistribution = proposedMCMCstring[table[i][j]]->getProposalDistribution();     
#endif
	      if(MCMCstring[table[i][j]] -> acceptance(temp,tempValue,tempProposalDistribution) == 1)
		MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	      else MCMCstring[table[i][j]] -> nextGenerationBlank();
	    }  
	  }
	}
#ifdef DEBUG_
	cout << "I finished the acceptance \n";
#endif
}


// !!!!this function do not consider different temperatures
//
//tempTable - temporal table with the values acceptance for (x',y) and (x',y) and (x',y') and (x,y)
//prop - store the result of acceptance operation
//tempStringValue - table with the fitness functiob for the involved chains
//tempTEMPvalue - table with the temperature for the involved chains
//size - size of coupling
//trece - table for the wining strategies

unsigned long binaryPopulation::acceptancePopulation(unsigned long* table){
	double prop = 1;
	double max = -1;

#ifndef GA_SELECTION
#ifdef IRREDUCTIBIL
	//cazul in care nu exista nici o modificare in generatia urmatoare
	if(MCMCstring[table[0]]->equal(proposedMCMCstring[table[0]]) && 
	   MCMCstring[table[1]]->equal(proposedMCMCstring[table[1]])){
	  MCMCstring[table[0]]->setACCEPTED(1);
	  MCMCstring[table[1]]->setACCEPTED(1);
	  return 1;
	}
	
	if(MCMCstring[table[0]]->equal(proposedMCMCstring[table[1]]) && 
	   MCMCstring[table[1]]->equal(proposedMCMCstring[table[0]])){
	  MCMCstring[table[0]]->setACCEPTED(1);
	  MCMCstring[table[1]]->setACCEPTED(1);
	  return 1;
	  }
#endif
#endif

	//calculeaza raporturile P(y)/P(x)
	for(unsigned long i = 0; i<sizeCoupling; i++) 
	  for(unsigned long j = 0; j<sizeCoupling*2; j++) {
	    if(i == j) tempResults[i][j] = 1;
	    else if(j < sizeCoupling) {
#ifdef Boltzmann
#ifdef Normal_Boltzmann
	      tempResults[i][j] = exp((double)(MCMCstring[table[j]] ->getValue() -  MCMCstring[table[i]] ->getValue()) / (double)MCMCstring[table[i]] -> getTEMPERATURE());		
#else
	      if(MCMCstring[table[j]] ->getValue() == 0 && MCMCstring[table[i]] ->getValue() == 0)
		tempResults[i][j] = 0;			
	      else{
		  if(MCMCstring[table[j]] ->getValue() == 0) MCMCstring[table[j]] ->getValue() == 0.01;
		  if(MCMCstring[table[i]] ->getValue() == 0) MCMCstring[table[i]] ->getValue() == 0.01;		
		  tempResults[i][j] = exp((log(MCMCstring[table[j]] ->getValue()) - log(MCMCstring[table[i]] ->getValue())) / (double)MCMCstring[table[i]] -> getTEMPERATURE());
	      }		
#endif
#else
	      if(MCMCstring[table[j]] ->getValue() == 0 && MCMCstring[table[i]] ->getValue() == 0)
		 tempResults[i][j] = 1;			
	      else{
		    if(MCMCstring[table[j]] ->getValue() == 0) MCMCstring[table[j]] ->getValue() == 0.01;
		    if(MCMCstring[table[i]] ->getValue() == 0) MCMCstring[table[i]] ->getValue() == 0.01;
		    tempResults[i][j] = (double)MCMCstring[table[j]] ->getValue() / (double)MCMCstring[table[i]] ->getValue();
	      }
#endif
	} else {
#ifdef Boltzmann 
#ifdef Normal_Boltzmann
	  tempResults[i][j] = exp((double)(proposedMCMCstring[table[j - sizeCoupling]]->getValue() -  MCMCstring[table[i]] ->getValue()) / (double) MCMCstring[table[i]] -> getTEMPERATURE());
#else
	  if(proposedMCMCstring[table[j - sizeCoupling]] ->getValue() == 0 && proposedMCMCstring[table[i]] ->getValue() == 0)
		tempResults[i][j] = 0;			
	  else{
	    if(proposedMCMCstring[table[j - sizeCoupling]] ->getValue() == 0) 
	      proposedMCMCstring[table[j - sizeCoupling]] ->getValue() == 0.01;
	    if(proposedMCMCstring[table[i]] ->getValue() == 0) 
	      proposedMCMCstring[table[i]] ->getValue() == 0.01;		 
	    tempResults[i][j] = exp((log(proposedMCMCstring[table[j - sizeCoupling]]->getValue()) - log(MCMCstring[table[i]] ->getValue())) / (double) MCMCstring[table[i]] -> getTEMPERATURE());
	  }
#endif
#else
	  if(proposedMCMCstring[table[j - sizeCoupling]] ->getValue() == 0 && proposedMCMCstring[table[i]] ->getValue() == 0)
	    tempResults[i][j] = 0;			
	  else{
	    if(proposedMCMCstring[table[j - sizeCoupling]] ->getValue() == 0) 
	      proposedMCMCstring[table[j - sizeCoupling]] ->getValue() == 0.01;
	    if(proposedMCMCstring[table[i]] ->getValue() == 0) 
	      proposedMCMCstring[table[i]] ->getValue() == 0.01;   
	    tempResults[i][j] = (double)proposedMCMCstring[table[j - sizeCoupling]]->getValue() / (double)MCMCstring[table[i]] ->getValue();
	  }
#endif
	}
	  }

#ifdef DETAILED_BALANCE	  
	  prop = (tempResults[0][2]*tempResults[1][3]);
	  if(prop < 1 && genrand_real1() > prop){ 
	    MCMCstring[table[0]]->setACCEPTED(1);
	    MCMCstring[table[1]]->setACCEPTED(1);
	    return 1;
	  }
#endif //DETAILED_BALANCE
	
	//maximum calculation
	double maxIndex[6][3];
	maxIndex[0][0] = tempResults[0][2]*tempResults[1][3]; maxIndex[0][1] = 2; maxIndex[0][2] = 3;
	maxIndex[1][0] = tempResults[0][2]                  ; maxIndex[1][1] = 2; maxIndex[1][2] = 1;
	maxIndex[2][0] = tempResults[1][3]                  ; maxIndex[2][1] = 0; maxIndex[2][2] = 3;
	maxIndex[3][0] = tempResults[0][2]*tempResults[1][0]; maxIndex[3][1] = 2; maxIndex[3][2] = 0;
	maxIndex[4][0] = tempResults[0][1]*tempResults[1][3]; maxIndex[4][1] = 1; maxIndex[4][2] = 3;
	maxIndex[5][0] = 1;                                   maxIndex[5][1] = 0; maxIndex[5][2] = 1;
       
	//ordering procedure
	for(unsigned long i = 0; i<6; i++) 
	  for(unsigned long j = 0; j<6; j++)
	    if(i != j)
	    if(maxIndex[i][0] < maxIndex[j][0] && i < j){
	      double temp = maxIndex[i][0];
	      maxIndex[i][0] = maxIndex[j][0];
	      maxIndex[j][0] = temp;

	      temp = maxIndex[i][1];
	      maxIndex[i][1] = maxIndex[j][1];
	      maxIndex[j][1] = temp;

	      temp = maxIndex[i][2];
	      maxIndex[i][2] = maxIndex[j][2];
	      maxIndex[j][2] = temp;
	    } 
	    else if(maxIndex[i][0] > maxIndex[j][0] && i > j){
	      double temp = maxIndex[i][0];
	      maxIndex[i][0] = maxIndex[j][0];
	      maxIndex[j][0] = temp;
	      
	      temp = maxIndex[i][1];
	      maxIndex[i][1] = maxIndex[j][1];
	      maxIndex[j][1] = temp;
	      
	      temp = maxIndex[i][2];
	      maxIndex[i][2] = maxIndex[j][2];
	      maxIndex[j][2] = temp;
	    }	

#ifdef DEBUG_
	cout << "tempResults is ";
	for(unsigned long i = 0; i<sizeCoupling; i++) {
	  cout << "\n for i=" << i << " | "; 
	  for(unsigned long j = 0; j<6; j++) 
	    cout << tempResults[i][j] << " | ";
	}
	cout << "\n";

	cout << "results for maximization are\n";
	for(unsigned long i = 0; i<6; i++) {
		cout << "|(" <<maxIndex[i][0] <<","<<maxIndex[i][1]<<","<<maxIndex[i][2]<<")";
	}
	cout << "\n";
	
	cout << "generation of trece: ";
#endif

#ifdef GA_SELECTION
	//cout << "trece afisare \n";
	for(unsigned long i = 1; i <= sizeCoupling; i++){
	  //cout <<  maxIndex[0][i] << "\t";
	  if(maxIndex[0][i] >= sizeCoupling)
	    proposedMCMCstring[table[(unsigned long)maxIndex[0][i] - sizeCoupling]]->setACCEPTED(1);
	  else MCMCstring[table[(unsigned long)maxIndex[0][i]]]->setACCEPTED(1);
	}
	//cout << "\n";
	return 1;
#endif

	double contor = 0;
	for(unsigned long i = 0; i < sizeCoupling; i++){
	  if(maxIndex[0][i + 1] < sizeCoupling) contor++;
#ifdef DEBUG_
	  cout << "i = " << maxIndex[0][i] << "\t"; 
#endif
	}

#ifdef DEBUG_
	cout<< " with contor = " << contor << "\n";
#endif	

	//if there are some new better states;
	if (contor < sizeCoupling) {
	    
	  for(unsigned long i = 0; i < sizeCoupling; i++){
	    if(maxIndex[0][i + 1] >= sizeCoupling){
	      proposedMCMCstring[table[(unsigned long)maxIndex[0][i+1] - sizeCoupling]]->setACCEPTED(1);
	    }
	    else{ 
		MCMCstring[table[(unsigned long)maxIndex[0][i + 1]]]->setACCEPTED(1);
	    }
	  }
#ifdef DEBUG_
	  cout << " some was accepted \n ";
	  for(unsigned long i = 0; i < sizeCoupling; i++){
	    if(proposedMCMCstring[(unsigned long)maxIndex[0][i + 1]]->getACCEPTED() == 1)
	      cout << " prop. i ==" <<i;
	    else cout << " not prop. i ==" <<i;
	    if(MCMCstring[(unsigned long)maxIndex[0][i + 1]]->getACCEPTED() == 1)
	      cout << " i ==" <<i;
	    else cout <<" not i ==" <<i;
	  }
	  cout << "\n";
#endif
	  return 1;
	}

#ifdef DEBUG_
	cout << " nobady was accepted so generates random ";
#endif

	//return -1;
	//else pick up some random combination
#ifdef MAX_SELECTION
	prop = maxIndex[1][0];
	if(genrand_real1() <= prop){ 
	  //cout << "accept the children \n";
	  for(unsigned long i = 0; i < sizeCoupling; i++){
	    if(maxIndex[1][i + 1] >= sizeCoupling)
	      proposedMCMCstring[table[(unsigned long)maxIndex[1][i + 1] - sizeCoupling]]->setACCEPTED(1);
	    else MCMCstring[table[(unsigned long)maxIndex[1][i + 1]]]->setACCEPTED(1);
	  }
	} else {
	  //cout << "accept the parents\n";
	  for(unsigned long i = 0; i < sizeCoupling; i++)
	    MCMCstring[table[(unsigned long)maxIndex[0][i + 1]]]->setACCEPTED(1);
	}
	return 1;
#else
#ifdef SIMPLE_SELECTION
	prop = 1;
	for(unsigned long i = 0; i < sizeCoupling; i++){
	  //prop *= tempResults[i][(i+1)%sizeCoupling];
	  prop *= tempResults[i][i+sizeCoupling];
	}
	//cout << "prop =" <<prop <<"\n";
	//if(genrand_real2() <= 5*prop){ 
	if(genrand_real1() <= prop){ 
	  //cout << "accept the childrens\n";
	  for(unsigned long i = 0; i < sizeCoupling; i++)
	    proposedMCMCstring[table[i]]->setACCEPTED(1);
	}
	else{
	  //cout << "accept the parents\n";
	  for(unsigned long i = 0; i < sizeCoupling; i++)
	    MCMCstring[table[i]]->setACCEPTED(1);
	}
	return 1;
#else
	proportionalSelection(table);
#endif
#endif
	return 1;
}

#ifdef WITH_PARENTS
//it is working only for m = 2 
void binaryPopulation::proportionalSelection(unsigned long* table)
{
  double tempTable[3];
  for(unsigned long i = 0; i < 3; i++) 
    tempTable[i] = -1;
  
  unsigned long index = 0;
  
  tempTable[0] = tempResults[0][2]*tempResults[1][3]; 
  tempTable[1] = tempResults[0][2]; 
  tempTable[2] = tempResults[1][3];

  double norm = 0;
  //compute tempTempResults
  for(unsigned long i = 0; i < 3; i++) 
	norm += tempTable[i];
 
  //find out the maximum and multiply by it
  double max = -1;
  index = 0;
  for(unsigned long i = 0; i < 3; i++)
    if(max < tempTable[i]) {
      max = tempTable[i];
      index = i;
    }

  //normalize the results from tempTable
  for(unsigned long i = 0; i < 3; i++)
    tempTable[i] *= max / norm; 

  //print 

#ifdef DEBUG_
  cout << " \n value Table for rejection \n";
  for(unsigned long i = 0; i < 3; i++)
    cout << tempTable[i] << " || "; 
    cout << " norm = "<< norm << " maxim = " << max;
#endif

  double temp = genrand_real1();

#ifdef DEBUG_
  cout << " \n    With ramdom number = "<<temp <<"\n";
#endif

  double leftNorm = 0;
  for(unsigned long i = 0; i < 3; i++) {
    if(temp > leftNorm + tempTable[i]) {
      leftNorm += tempTable[i]; 
      continue;
    }
    else {
      //fill in the results vector
      unsigned long g = 0, t = 0;
      if(i == 0)      {g = 2; t = 3;}
      else if(i == 1) {g = 2; t = 1;}
      else if(i == 2) {g = 3; t = 0;}
      if(t >= sizeCoupling)
	proposedMCMCstring[table[t - sizeCoupling]]->setACCEPTED(1);
      else MCMCstring[table[t]] -> setACCEPTED(1);
      if(g >= sizeCoupling)
	proposedMCMCstring[table[g- sizeCoupling]]->setACCEPTED(1);
      else MCMCstring[table[g]] -> setACCEPTED(1);
      return;
    }
  }
  
  return;
}

#else 
//it is working only for m = 2 
void binaryPopulation::proportionalSelection(unsigned long* table)
{
  double tempTable[5];
  for(unsigned long i = 0; i < 5; i++) 
    tempTable[i] = -1;
  
  unsigned long index = 0;

  tempTable[0] = tempResults[0][2]*tempResults[1][3]; 
  tempTable[1] = tempResults[0][2]; 
  tempTable[2] = tempResults[1][3];
  tempTable[3] = tempResults[0][2]*tempResults[1][0]; 
  tempTable[4] = tempResults[0][1]*tempResults[1][3];

  double norm = 0;
  //compute tempTempResults
  for(unsigned long i = 0; i < 5; i++) 
	norm += tempTable[i];
 
  //find out the maximum and multiply by it
  double max = -1;
  index = 0;
  for(unsigned long i = 0; i < 5; i++)
    if(max < tempTable[i]) {
      max = tempTable[i];
      index = i;
    }

#ifdef DEBUG_
  cout << " \n value Table for rejection \n";
  for(unsigned long i = 0; i < 5; i++)
    cout << tempTable[i] << "||"; 
    cout << "\n";
#endif
  
  //normalize the results from tempTable
  for(unsigned long i = 0; i < 5; i++)
    tempTable[i] *= max / norm; 

#ifdef DEBUG_
  //print 
  cout << " \n value Table for rejection after normalization\n";
  for(unsigned long i = 0; i < 5; i++)
    cout << tempTable[i] << "||"; 
  cout << " norm = "<< norm << " maxim = " << max;
#endif

  double temp = genrand_real1();

  //#ifdef DEBUG_
  // cout << " \n    With ramdom number = "<<temp << " max "<< max <<"\n";
  //#endif

  if(temp > max) {
    return;
  }

  double leftNorm = 0;
  for(unsigned long i = 0; i < 5; i++) {
    if(temp >= leftNorm + tempTable[i]) {
      leftNorm += tempTable[i]; 
      continue;
    }
    else {
      //fill in the results vector
      unsigned long g = 0, t = 0;
      if(i == 0)      {g = 2; t = 3;}
      else if(i == 1) {g = 2; t = 1;}
      else if(i == 2) {g = 3; t = 0;}
      else if(i == 3) {g = 2; t = 0;}
      else if(i == 4) {g = 1; t = 2;}
      if(t >= sizeCoupling)
	proposedMCMCstring[table[t - sizeCoupling]]->setACCEPTED(1);
      else MCMCstring[table[t]] -> setACCEPTED(1);
      if(g >= sizeCoupling)
	proposedMCMCstring[table[g-sizeCoupling]]->setACCEPTED(1);
      else MCMCstring[table[g]] -> setACCEPTED(1);
      return;
    }
  }
  
  cout << "erroare \n";
  return;
}
#endif

#else

#ifdef PARALLEL_SIMULATED_ANNEALING
//------------------Parallel Simulated annealing----------------

void binaryPopulation::evolution()
{

  unsigned long tempSize = sizeMCMCstring / sizeCoupling;
  unsigned long** table = new unsigned long*[tempSize + 1];
  for(unsigned long i = 0; i< tempSize + 1; i++)
    table[i] = new unsigned long[sizeCoupling];
  
  cout << "Parallel simulated annealing ";
  //cin >> bestTemperature;

#ifdef DEBUG_
	ofstream myFile(fileData);		
	if(!myFile.is_open())
	  cout << "Error opening file";
#endif	

#ifdef EXPAND_KULLBACK_INFORMATION
    expandTrueUnivariateDistribution();
#else
    readHistogram();
#endif

	while(runs < nr_runs){
	  runs++;
	  
	  setSeed();
	  while(generation <=  burn_in + sampleSize*interSamplesize){	
	      assignTemperature();	      
	      unsigned long dimension = sizeMCMCstring/sizeCoupling;
#ifndef NEAR_CHANGING  
	      randomCOUPLING(table);
#else
	      dimension = neighborsCOUPLING_NO_OVERLAPING(table);
#endif
	      proposedGeneration(table,dimension);
	      acceptanceProposed(table,dimension);
	      postProcessing();

	      simul = insert_simulation(simul,generation,procentRightSolutions);
	      //show_simulations(simul);
	      write_to_file_all(simul_all_proc_right,simul,NULL);

#ifdef KULLBACK_INFORMATION
	mixingTimeUnivariateCalculation();      
#endif //KULLBACK_INFORMATION

	      generation++;	
	    }	

#ifdef DEBUG_
    printLEVEL(myFile,tempData);
#else 
    scoreDiversity();
#endif

    globalRightSolutions[runs-1] = procentRightSolutions;
#ifdef CUSUM
	CusumCalculation();
#endif

#ifdef GELMAN_RUBIN
  varianceHigerMomentsUnivariateCalculation(GELMAN_Moment,sizeMCMCstring - 1);
  varianceMultivariateCalculation();
#endif

    reset();
  }

#ifdef GELMAN_RUBIN
	printReductionGNUFILE(GELMAN_Moment);
	printReductionMultivariateGNUFILE();
#endif
       
	printDiversityGNUFILE();
#ifdef TRAP_LIKE_FUNCTION
	printDiversityDetailGNUFILE();
#endif
	printDiversityRightSolutions();

#ifdef KULLBACK_INFORMATION
	printMixingGNUFILE();
#endif

#ifdef CUSUM
	printCusumGNUFILE();
#endif

#ifdef DEBUG_
  myFile.close();
#endif

#ifndef FILE_BOOKEEPING
  diversityFile.flush();
  diversityFile.close();
#endif

  delete[] table;	
}

#ifdef COUPLED_ACCEPTANCE
unsigned long binaryPopulation::acceptancePopulation(unsigned long* table){
  double prop = 1;
  for(unsigned long i = 0; i<sizeCoupling; i++) {
    prop = prop * 
      exp((proposedMCMCstring[table[i]]->getValue() - MCMCstring[table[i]]->getValue())
	  /MCMCstring[table[i]]->getTEMPERATURE());
  }	
  if(prop > 1) return 1;
  double temp = genrand_real1();
  if(temp > prop) {
    return 0;
  }
  return 1;
}

void binaryPopulation::acceptanceProposed(unsigned long** table, unsigned long dimension)
{
  unsigned long tempSize = dimension;
  unsigned long temp[sizeGENOME];
  for(unsigned long i = 0; i<tempSize ; i++){
    for(unsigned long j = 0; j<sizeCoupling; j++){
      proposedMCMCstring[table[i][j]]->setACCEPTED(0); 
      MCMCstring[table[i][j]]->setACCEPTED(0);
    }
    
    if(COUPLING_PROBABILITY > genrand_real1()){
      if(acceptancePopulation(table[i]) == 1){
	for(unsigned long j = 0; j<sizeCoupling; j++){
	  proposedMCMCstring[table[i][j]]->getGENOME(temp);
	  double tempValue = proposedMCMCstring[table[i][j]]->getValue();
	  MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	}
	//cout << "I accepted sometihng for " << i <<"\n";
	continue;
      } else {
	//   cout << "I reject for "<< i <<"\n";
	for(unsigned long j = 0; j < sizeCoupling; j++){
	  MCMCstring[table[i][j]] -> nextGenerationBlank();
	}
      }
    }
    else {
      for(unsigned long j = 0; j < sizeCoupling; j++){
        proposedMCMCstring[table[i][j]]->getGENOME(temp);
	double tempValue = proposedMCMCstring[table[i][j]]->getValue();
#ifdef METROPOLIS_ALGORITHM
	double tempProposalDistribution = MCMCstring[table[i][j]]->getProposalDistribution();
#else
	double tempProposalDistribution = proposedMCMCstring[table[i][j]]->getProposalDistribution();     
#endif
	if(MCMCstring[table[i][j]] -> acceptance(temp,tempValue,tempProposalDistribution) == 1)
	  MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
	else MCMCstring[table[i][j]] -> nextGenerationBlank();
      }  
    }
  }
}
#else
//single Boltzman competitions
void binaryPopulation::acceptanceProposed(unsigned long** table, unsigned long dimension)
{
  unsigned long tempSize = dimension;
  unsigned long temp[sizeGENOME];
  
  for(unsigned long i = 0; i<tempSize ; i++)
    for(unsigned long j = 0; j<sizeCoupling; j++){
      proposedMCMCstring[table[i][j]]->setACCEPTED(0); 
      MCMCstring[table[i][j]]->setACCEPTED(0);
      
      proposedMCMCstring[table[i][j]]->getGENOME(temp);
      double tempValue = proposedMCMCstring[table[i][j]]->getValue();           
#ifdef METROPOLIS_ALGORITHM
      double tempProposalDistribution = MCMCstring[table[i][j]]->getProposalDistribution();
#else
      double tempProposalDistribution = proposedMCMCstring[table[i][j]]->getProposalDistribution();     
#endif
      if(MCMCstring[table[i][j]]->acceptance(temp,tempValue,tempProposalDistribution) == 1)
	MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
      else 
	MCMCstring[table[i][j]] -> nextGenerationBlank();
    }
}
#endif
#else

//--------------------POPULATION of MCMC------------------------------
#ifdef POPULATION_MCMC
void binaryPopulation::evolution()
{

  unsigned long tempSize = sizeMCMCstring / sizeCoupling;
  unsigned long** table = new unsigned long*[tempSize + 1];
  for(unsigned long i = 0; i< tempSize + 1; i++)
    table[i] = new unsigned long[sizeCoupling];

#ifdef DEBUG_
	ofstream myFile(fileData);		
	if(!myFile.is_open())
	  cout << "Error opening file";
#endif
	
#ifdef EXPAND_KULLBACK_INFORMATION
    expandTrueUnivariateDistribution();
#else
    readHistogram();
#endif

	while(runs < nr_runs){
	  runs++;
	  setSeed();
	  assignTemperature();
	  while(generation <=  (burn_in + sampleSize*interSamplesize)*sizeMCMCstring)
	    {		      
#ifdef DEBUG_
	      cout << "generation " << generation << "\n";
	      print();
#endif	      
	unsigned long dimension = sizeMCMCstring/sizeCoupling;
#ifndef NEAR_CHANGING  
	randomCOUPLING(table);
#else
	dimension = neighborsCOUPLING_NO_OVERLAPING(table);
#endif
	proposedGeneration(table,dimension);
	acceptanceProposed(table,dimension);
	postProcessing();
	
	simul = insert_simulation(simul,generation,procentRightSolutions);
	//show_simulations(simul);
	write_to_file_all(simul_all_proc_right,simul,NULL);
	
#ifdef KULLBACK_INFORMATION
	mixingTimeUnivariateCalculation();      
#endif //KULLBACK_INFORMATION

	generation++;	

	    }	
#ifdef DEBUG_
	  printLEVEL(myFile);
#else 
    scoreDiversity();
#endif

	  globalRightSolutions[runs-1] = procentRightSolutions;

#ifdef CUSUM
	    CusumCalculation();
#endif

#ifdef GELMAN_RUBIN
  varianceHigerMomentsUnivariateCalculation(GELMAN_Moment,sizeMCMCstring  -1);
  varianceMultivariateCalculation();
#endif

	  reset();
	}
	
	printDiversityGNUFILE();

#ifdef TRAP_LIKE_FUNCTION
	printDiversityDetailGNUFILE();
#endif
	printDiversityRightSolutions();
	
#ifdef KULLBACK_INFORMATION
	printMixingGNUFILE();
#endif

#ifdef CUSUM
	printCusumGNUFILE();
#endif

#ifdef GELMAN_RUBIN
	printReductionGNUFILE(GELMAN_Moment);
	printReductionMultivariateGNUFILE();
#endif

#ifdef DEBUG_
  myFile.close();
#endif

#ifndef FILE_BOOKEEPING
  diversityFile.flush();
  diversityFile.close();
#endif

  delete[] table;	
}

void binaryPopulation::acceptanceProposed(unsigned long** table, unsigned long dimension)
{
  unsigned long tempSize = dimension;
  unsigned long temp[sizeGENOME];

  for(unsigned long i = 0; i<tempSize ; i++)
    for(unsigned long j = 0; j<sizeCoupling; j++){
      proposedMCMCstring[table[i][j]]->setACCEPTED(0); 
      MCMCstring[table[i][j]]->setACCEPTED(0);
      
      proposedMCMCstring[table[i][j]]->getGENOME(temp);
      double tempValue = proposedMCMCstring[table[i][j]]->getValue();           
#ifdef METROPOLIS_ALGORITHM
      double tempProposalDistribution = MCMCstring[table[i][j]]->getProposalDistribution();
#else
      double tempProposalDistribution = proposedMCMCstring[table[i][j]]->getProposalDistribution();     
      //double tempProposalDistribution = MCMCstring[table[i][j]]->getProposalDistribution();     
#endif
      if(MCMCstring[table[i][j]]->acceptance(temp,tempValue,tempProposalDistribution) == 1){
	  MCMCstring[table[i][j]] -> nextGeneration(temp,tempValue,sizeGENOME);
#ifdef DEBUG_
	  cout << "I accepted " << table[i][j] <<"\n";
#endif
      }
      else 
	  MCMCstring[table[i][j]] -> nextGenerationBlank();
    }
}
#endif
#endif
#endif
#endif
#endif
#endif
#endif

/**********/
//Extra functions
#ifdef MULTIPLE_RUNS
void binaryPopulation::printMultiRuns(ofstream& myFile){
#ifdef SINGLE_vs_MULTIPLE
#ifdef TWO_ATTRACTORS_FUNCTION
    myFile << sizeGENOME << "\t" << MCMCstring[0]->temperature << "\t" <<
	sampleSize * sizeMCMCstring << "\t";
    for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2))+1 ;i++){
	myFile << pow(2,i) << "\t" << meanFitness[i] << "\t" << diffFitness[i] 
	       << "\t" << meanAcceptance[i] << "\t"; 
    }
    myFile << "\n";
#else
    myFile << sizeGENOME << "\t" << MCMCstring[0]->temperature << "\t" <<
	sampleSize * sizeMCMCstring  << "\t" << 
	meanFitness << "\t" << diffFitness << "\t" << meanAcceptance << "\n";
#endif
#else
    myFile << sizeGENOME << "\t" << MCMCstring[0]->temperature << "\t" <<
#ifdef IsRecombination
	UNIFORM_RECOMB << "\t" <<
#else
	0 << "\t" << 
#endif 
	sampleSize * sizeMCMCstring  << "\t" << 
	meanFitness << "\t" << diffFitness << "\t" << meanAcceptance << "\n";
#endif
}
#endif

#ifdef PLOTS_MCMC
void binaryPopulation::printAllChains(ofstream& myFile){
    if(MCMCstring != NULL){
	for(unsigned long level = 1; level <= sampleSize; level++){
	    for(unsigned long i =0; i< sizeMCMCstring; i++)
		myFile << MCMCstring[i]->getValue(level) << "\t";
	    myFile << "\n";
	}
    }
}

void binaryPopulation::printHashtable(ofstream& myFile){
   long nrIndiv = 0;
//    cout << "verification hashtable \n";
#ifndef MIXTURE_REALS
	unsigned long tempBuffer[sizeGENOME+1];
#else
	unsigned long tempBuffer[sizeGENOME*BLOCKsize];
	double reals[sizeGENOME];
#endif

    if(MCMCstring != NULL){
	myFile << "mutation: " << MCMCstring[0]->mutation << "\n";
	myFile << "size population: " << sizeMCMCstring << "\n";
	myFile <<"Nr runs = " << runs << "\n Nr. of discovered indiv " 
	       << myHashtables[runs] ->nrElem << "\n";
	
#ifndef MIXTURE_REALS
	nrIndiv = myHashtables -> iteratorInit(tempBuffer,sizeGENOME); 
	double tempValue = MCMCstring[0] -> fitness(tempBuffer,sizeGENOME);

	for(unsigned long i = 0; i < sizeGENOME; i++) 
	    myFile << tempBuffer[i];
	myFile << "\t" << tempValue<<"\t" << nrIndiv << "\n";
#else
#ifdef HISTOGRAM
	nrIndiv = myHashtables -> iteratorInit((unsigned long*)tempBuffer,sizeGENOME*BLOCKsize); 
	//for(int i = 0; i < sizeGENOME*BLOCKsize; i++) {
	//    cout << tempBuffer[i] << "\t";
	//}
	//cout << "\n";
	convertReal(reals,tempBuffer,sizeGENOME);
	double tempValue = MCMCstring[0] -> fitness(reals,sizeGENOME);

	for(unsigned long i = 0; i < sizeGENOME; i++) {
	    myFile << reals[i] << "\t(";
	    for(unsigned long j = i*BLOCKsize; j < (i+1)*BLOCKsize; j++)
		myFile << tempBuffer[j];
	    myFile << ")";
	    //    cout << reals[i] << "\t";
	}
	myFile << "\t = " << tempValue<<"\t" << nrIndiv << "\n";
	//cout << "\t = " << tempValue<<"\t" << nrIndiv << "\n";
#endif //histogram
#endif
		
#ifndef MIXTURE_REALS
	while((nrIndiv = myHashtables ->iterator(tempBuffer,sizeGENOME)) != -1){
	    //  for(int i = 0; i < sizeGENOME; i++) 
	    //cout << tempBuffer[i];
	    //cout << "\n";
	    double tempValue = MCMCstring[0] -> fitness(tempBuffer,sizeGENOME);
	    //if(tempValue >= 3){
	    for(unsigned long i = 0; i < sizeGENOME; i++){ 
		myFile << tempBuffer[i];
	    }
	    myFile <<  "\t" << tempValue << "\t" << nrIndiv << "\n";
	    //}
	}
	myHashtables->iteratorFree(sizeGENOME);
	myFile << "\n";
#else
#ifdef HISTOGRAM
	while((nrIndiv = myHashtables ->iterator((unsigned long*)tempBuffer,sizeGENOME*BLOCKsize)) != -1){
	    //for(int i = 0; i < sizeGENOME*BLOCKsize; i++) {
	    //cout << tempBuffer[i] << "\t";
	    //}
	    //cout << "\n";
	    convertReal(reals,tempBuffer,sizeGENOME);
	    double tempValue = MCMCstring[0] -> fitness(reals,sizeGENOME);
	    for(unsigned long i = 0; i < sizeGENOME; i++){ 
		myFile << reals[i] << "\t(";
		for(unsigned long j = i*BLOCKsize; j < (i+1)*BLOCKsize; j++)
		    myFile << tempBuffer[j];
		myFile << ")";
		//cout << reals[i] << "\t";
	    }
	    myFile <<  "\t = " << tempValue << "\t" << nrIndiv << "\n";
	    //cout <<  "\t = " << tempValue << "\t" << nrIndiv << "\n";
	}
	myHashtables->iteratorFree(sizeGENOME*BLOCKsize);
	myFile << "\n";
#endif //histogram
#endif //MIXTURE_REALS
    }
}

void binaryPopulation::printOneChain(ofstream& myFile, unsigned long i){
  if(MCMCstring != NULL){
    for(unsigned long level = 1; level <= sampleSize; level++)
	    myFile << MCMCstring[i]->getValue(level) << "\n";
  }
}

void binaryPopulation::printMultiValuedChain(ofstream& myFile){
  if(MCMCstring != NULL){
      for(unsigned long level = 1; level <= sampleSize; level++){
	  for(unsigned long i =0; i< sizeMCMCstring; i++){
	    myFile << MCMCstring[i]->getValue(level) << "\t";
	  }
	  myFile << "\n";
      }
  }
}

void binaryPopulation::printAllChains_I(ofstream& myFile){
  if(MCMCstring != NULL){
	for(unsigned long i =0; i< sizeMCMCstring; i++)
	    myFile << MCMCstring[i]->getValue() << "\n";
  }
}

void binaryPopulation::printOneChain_I(ofstream& myFile, unsigned long i){
    if(MCMCstring != NULL){
	myFile << MCMCstring[i]->getValue() << "\n";
    }
}

void binaryPopulation::printMultiValuedChain_I(ofstream& myFile){
  if(MCMCstring != NULL){
      for(unsigned long i =0; i< sizeMCMCstring; i++){
	    myFile << MCMCstring[i]->getValue() << "\t";
      }
      myFile << "\n";
  }
}

#else
#ifdef DEBUG_
void binaryPopulation::printLEVEL(ofstream& myFile){
  double tempValue = 0;
  if(MCMCstring != NULL){
    myFile << "nr. of runs= " << runs << "\n";

    for(unsigned long level = 1; level <= sampleSize; level++){
      myFile << "level=" << level << "\n";
      for(unsigned long i=0;i<sizeMCMCstring;i++) {
	myFile << i << "-th MCMC chain: ";
	//cout << "MCMC chain i="<<i<<"\t";
	//MCMCstring[i] -> printLEVEL(level,myFile,tempData);	
  
	double t = MCMCstring[i]->getElementValue(level);
	if (t == -1)
	  cout << "Error in finding the level " << level << " of the element " << i << "\n";
	else tempValue += t;

	MCMCstring[i] -> printLEVEL(level,myFile);	
      }
    }

    myFile << generation << " generations  \n";
    myFile << bestFitness << "best fitness " <<worstFitness <<" worstFitness \n";
    myFile << "procent hitted values " << procentRightSolutions << "values of it" << weightRightSolutions << "\n";
    myFile << "Diversity score=" << scoreDiversity(); 
#ifdef TRAP_FUNCTION
    myFile << " Detailed diversity score=" << scoreDiversityDetail();
#endif
    myFile << "\nDiversivity tables \n";
    for(unsigned long i = 0; i<diversityMeasure; i++)
      myFile << "("<<i<<","<<diversity[i]<<")\t";
    myFile <<"\n";
#ifdef TRAP_FUNCTION
    for(unsigned long i = 0; i<sizeGENOME/BLOCKsize; i++){
      myFile << "Block " << i;
      for(unsigned long j = 0; j<BLOCKsize+1; j++)
	myFile << "("<<i<<","<<j<<","<<diversityDetail[i][j]<<")\t";
      myFile << "\n";
    }
#endif
  } 
  else cout << "error MCMC empty \n";	
  //cout << "I finished the print for 1 chain \n";
}

void binaryPopulation::print(){
  cout << " Population \n";
  if(MCMCstring != NULL)
    for(unsigned long i=0;i<sizeMCMCstring;i++) {
      cout << i << "-th MCMC chain: \n";
      MCMCstring[i] -> print();
    } else cout << "error MCMC empty \n";
  cout << " NEW Population \n";
  if(proposedMCMCstring != NULL)
    for(unsigned long i=0;i<sizeMCMCstring;i++) {
      cout << i << "-th MCMC chain: \n";
      proposedMCMCstring[i] -> print();
    } else cout << "error proposed MCMC empty \n";

  //cout << bestTemperature << "best temperature " <<worstTemperature <<" worstTemperature \n";
  cout << bestFitness << "best fitness " <<worstFitness <<" worstFitness \n";
  cout << generation << "generations \n";	
}

void binaryPopulation::print(char* outFile){
  ofstream myFile(outFile);		
  if(!myFile.is_open()){
    cout << "Error opening file";
    return;
  }
  if(MCMCstring != NULL) 
    {
      for(unsigned long i=0;i<sizeMCMCstring;i++) {
	myFile << i << "-th MCMC chain: \n";
	MCMCstring[i] -> print(myFile);	
      }
      myFile << generation << " generations  \n";
      myFile << bestFitness << "best fitness " <<worstFitness <<" worstFitness \n";
    } else cout << "error MCMC empty \n";		
  myFile.close();
}
#endif 
#endif

#ifdef SINGLE_vs_MULTIPLE
double binaryPopulation::updateHashtable(){

    for(unsigned long i =0; i < diversityMeasure; i++){
	diversity[i] = 0;
	diversityRightSolutions[i] = 0;
    }

    procentRightSolutions = 0;
    weightRightSolutions = 0;

#ifdef TRAP_FUNCTION
    for(unsigned long j1= 0; j1 < sizeGENOME/BLOCKsize; j1++)
	for(unsigned long k = 0; k < BLOCKsize; k++)
	    diversityDetail[j1][k] = 0;
#endif    

//#ifdef TRAP_FUNCTION
    //detailed diversity bookeeping
    unsigned long temp[sizeGENOME];
//#endif
    
    cout <<" update the hashtable size population:" << sizeMCMCstring 
	 << " sample size:" << sampleSize << " \n";
    myHashtables -> reset();
    
    //test the existing histogram
    double checkHistogram[sizeGENOME+1];
    for(unsigned long i= 0; i < sizeGENOME+1; i++)
	checkHistogram[i]= 0;

    double firstBestFitness = 0;
    double secondBestFitness = 0;

    for(unsigned long i= 0; i < sizeMCMCstring; i++)
	for(unsigned long j=1; j < sampleSize/2; j++){
		MCMCstring[i]->getGENOME(j,temp);
		
		double tempValue = MCMCstring[0] -> fitness(temp,sizeGENOME);

		if(myHashtables -> check_individual(temp,sizeGENOME) == 0){
		    //cout << "sample j:" << j << " chain:" << i << "\t";
		    //for(int i = 0; i < sizeGENOME; i++) 
		    //cout << temp[i];
		    procentRightSolutions++;
		    
		    for(unsigned long t =0; t < diversityMeasure; t++)
			if(tempValue == transform_table[t])
			    diversityRightSolutions[t]++;
		    
		    weightRightSolutions += tempValue;
		    
		    unsigned long oneC = 0;
		    for(unsigned long k = 0; k < sizeGENOME; k++)
			oneC += temp[k];
		    ActualHistogram[(unsigned long)(log(sizeMCMCstring)/log(2))][runs-1][oneC]++;
		    double stop = MCMCstring[0] -> stop(temp,sizeGENOME);

		    if(stop == 1 && firstBestFitness == 0){ 
#ifndef TWO_ATTRACTORS_FUNCTION		    
			firstBestFitness = 1;
			bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1] = j;
			//cout << " best 1 in best" << 
			//    bestFitnessVector[(int)(log(sizeMCMCstring)/log(2.0))][runs-1] 
			//     << " with chains" << (int)(log(sizeMCMCstring)/log(2.0)) << 
			//    " at runs" << (runs -1);
#else		  
		    firstBestFitness = 1;
		    bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][0] = j;
		    cout << " best 1 in best" << 
		   	bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][0] 
		    	 << " with chains" << (unsigned long)(log(sizeMCMCstring)/log(2.0)) << 
		    	" at runs" << (runs -1) << "in chain 0\n";
		    }
		    if(stop == 2 && secondBestFitness == 0){
		      secondBestFitness = 1;
		      bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][1]=j;
		      	cout << " best 2 in best" << 
		        bestFitnessVector[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][1] 
		           << " with chains" << (unsigned long)(log(sizeMCMCstring)/log(2.0)) << 
			    " at runs" << (runs -1) << "in chain 0\n";
#endif
		    }
		    //cout << "\n";
		}
		
		myHashtables -> store_individual(temp,sizeGENOME);

		//myHashtables->store_individual(tempValue,sizeGENOME);
		for(unsigned long t =0; t < diversityMeasure; t++)
		    if(tempValue == transform_table[t]){
			diversity[t]++;
		    }
#ifdef TRAP_FUNCTION
		for(unsigned long j1= 0; j1 < sizeGENOME/BLOCKsize; j1++){
		    unsigned long nr_temp = 0; 
		    for(unsigned long k = 0; k < BLOCKsize; k++)
			nr_temp += temp[j*BLOCKsize+k];
		    //cout << "("<<j<<","<<nr_temp<<")\t";
		    diversityDetail[j1][nr_temp]++;
		}
#endif
		
	}
//verifica daca histogramele sint aceasi
/*    if(sizeMCMCstring == 1)
    for(int k = 0; k < sizeGENOME+1; k++)
	if(checkHistogram[k] != ActualHistogram[0][runs-1][k]){
	    cout << "Histograma diferita binaryPopulation_trap.cpp:4456 " << k <<
		" with chack " << checkHistogram[k]<< " and real" << 
		ActualHistogram[0][runs-1][k]<< "\n";
	    exit(1);
	}
*/
    cout << "Terminat update \n";
    return 1;
}

#endif
