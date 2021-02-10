//test program for binaryPopulation
//#include <stdio.h>
#include<iostream.h>

#include "binaryPopulation_trap.h"
//#include "sample.h"

int main ()
{
   unsigned long nrVariables,nrVariablesMax,incrementVariables;	
    unsigned long size;
    char buffer[100];
    unsigned long sampleSizeA;
    srand(time(0));
    //  srand48(time(0));

//#ifdef MULTIPLE_RUNS  
//    cout << "To number of variables: ";
//    cin >> nrVariablesMax;

//    cout << "from how many to how many: ";
//    cin >> incrementVariables;

    double lowTemperature;
    double highTemperature;
    double incrementTemperature;

    //   cout << " Lowest temperature:";
    //cin >> lowTemperature;
    
    //cout << "Highest temperature";
    //cin >> highTemperature;

    //cout << "Increment temperature";
    //cin >> incrementTemperature;
//#else 
    cout << "Enter number of variables: ";
    cin >> nrVariables;

    cout << "Enter number of chains: ";
    cin >> size;
    
    cout << "Enter the output file name: ";
    cin >> buffer;
//#endif

    //sample* s = new sample(nrVariables);
    //s -> run();
    //delete s;

    //#ifndef MULTIPLE_RUNS
     binaryPopulation* p1 = new binaryPopulation(size,nrVariables,buffer);  
    p1 -> evolution();
#ifdef PLOTS_MCMC
    //p1->print(buffer);
#else
#ifdef DEBUG_
    p1->print(buffer);
#endif
#endif
    delete p1;    
    
    /*#else //ifdef MULTIPLE_RUNS
    //fisier in care se sciu datele de la rulari
        ofstream myMultiFile("test_2chain",ios::out | ios:: app);		
	  if(!myMultiFile.is_open()){
	  cout << "Error opening file \n";
	  exit(1);
	  }
	  
	  size = 2;
	  #ifdef TEMPERATURE_CHANGE
	  unsigned long sizeGENOMEVector[SIZEGENOME_CT] = {15};
	  double temperatureVector[TEMPERATURE_CT] = {1.0};
	  //   for(int r = 0; r <= 1; r++){
	  //	if(r == 0){
	  //#define IsRecombination
	  //	} else{
	  //#undef IsRecombination
	  //	}
	  //double temperatureVector[1] = {8};
	  for(int i = 0; i < SIZEGENOME_CT; i++){

	  for(int j = 0; j< TEMPERATURE_CT; j++){
	  if(j == 0){
	  #define EXPAND_KULLBACK_INFORMATION
	  }else if(j==1){
	  #undef EXPAND_KULLBACK_INFORMATION
	  //if(j != 0) p1->resetRuns(sizeGENOMEVector[i]);
	  }
	  binaryPopulation* p1 = new binaryPopulation(size,sizeGENOMEVector[i],"err");
	  p1->bestTemperature = temperatureVector[j];
	  p1->worstTemperature = temperatureVector[j];
	  for(int k = 0; k < size; k++)
	  p1->MCMCstring[k]->temperature = temperatureVector[j];

	  p1->evolution();
	  p1->printMultiRuns(myMultiFile);
	  myMultiFile.flush();
	  cout << "it is finished the runs for the size of the genome "<<
	  p1->sizeGENOME <<
	  " and temperature" << p1->MCMCstring[0]->temperature << "\n";
	  //cin >> r;
	  delete p1;
	  }

	  }
	  //}
	  #else // other
*/
    //fisier in care se sciu datele de la rulari
/*    ofstream myMultiFile("bin_recomb_mut",ios::out | ios:: app);
    if(!myMultiFile.is_open()){
	cout << "Error opening file \n";
	exit(1);
    }

	binaryPopulation* p1 = new binaryPopulation(size,nrVariables,buffer);
#ifdef MUTATION_VARIES
	myMultiFile << "\n #mutation " << size << "\t" << sampleSize_Const << "\n" ;
#else
	myMultiFile << "\n #Recomb unif " << nrVariables << "\t" << size
		    << "\t" << sampleSize_Const << "\n" ;
#endif
	for(int i = 8; i <= 10; i++){
	    //   for(int i = 6; i <= (int)((double)nrVariables/2.0); i++){
	//for(int i = 1; i <= ; i++){
//		if(i == 1){
//#define EXPAND_KULLBACK_INFORMATION
//		}else if(i==2){
//#undef EXPAND_KULLBACK_INFORMATION
		    //if(j != 0) p1->resetRuns(sizeGENOMEVector[i]);
//		}
#ifdef MUTATION_VARIES
	  //#undef mutation_multi
#define mutation_multi i
		myMultiFile << mutation_multi << "\t";

		for(int k = 0; k < p1->sizeMCMCstring; k++){
		    p1->MCMCstring[k]->mutation = (1.0 / (double) nrVariables) * i;
		    p1->proposedMCMCstring[k]->mutation =  (1.0 / (double) nrVariables) * i;
		    //p1->MCMCstring[k]->UNIF_RECOMB_PARAM = 0.5;
		}
#else

//#undef PROB_nPOINT
//#define PROB_nPOINT i

		myMultiFile << (double)i/(double)nrVariables << "\t";
		for(int k = 0; k < p1->sizeMCMCstring; k++)
		    p1->MCMCstring[k]->UNIF_RECOMB_PARAM = (double)i/(double)nrVariables;
		    //p1->MCMCstring[k]->PROB_nPOINT = i;

#endif 
       		p1->evolution(myMultiFile);
       		//p1->evolution();
		p1->resetRuns(size);
		
		//myMultiFile << "\n";
		myMultiFile.flush();
	}
	delete p1;
	myMultiFile.flush();
	myMultiFile.close();
	
	//#endif // MULTI_RUNS
//#endif
*/
	return 0;
}

