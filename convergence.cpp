#include "binaryPopulation_trap.h"
#include "PARAMETERS.h"

#include "Utils\long_doubleList.h"

extern double genrand_real1();
extern void convertReal(double*,unsigned long*, long);
extern void convertBinary(double*,unsigned long*, long);

//metoda lui Gelman si Rubin cu univariate distribution
#ifdef GELMAN_RUBIN

#ifdef NORMAL_DISTRIBUTION
void binaryPopulation::varianceCalculation(){
  
  double tempVal = 0;
  
#ifdef HALF_GELMAN
  for(unsigned long j = 1; j <= (sampleSize + 1)/constBurn; j++){
    teta[runs-1][j -1] = 0;
    for(unsigned long i = 0; i < sizeMCMCstring; i++){
      tempVal = MCMCstring[i]->getElementValue(j*constBurn,j*constBurn/2 + 1);
      if(tempVal != -1)
	teta[runs-1][j -1] += tempVal;
      else {
	cout << "It is an error in varaince Calculation \n";
	return;
      }
    }
    //cout << teta[runs-1][j - localBurn - 1] << " for runs=" << runs 
    //	   << " and for j=" << j << "and local burn " << localBurn << " after teta";
    teta[runs-1][j - 1] /= j*constBurn*sizeMCMCstring/2;
    //cout << teta[runs-1][j - 1] << " for runs=" << runs << " and for j=" << j << "\n";
  }
#else
    for(unsigned long j = 2; j <= sampleSize + 1; j++){
      teta[runs-1][j - 2] =0;
      for(unsigned long i = 0; i < sizeMCMCstring; i++){
	tempVal = MCMCstring[i]->getElementValue(j,1);
	if(tempVal != -1)
	  teta[runs-1][j - 2] += tempVal;
	else {
	  cout << "It is an error in varaince Calculation \n";
	  return;
	}
      }
      //cout << teta[runs-1][j - localBurn - 1] << " for runs=" << runs 
      //	   << " and for j=" << j << "and local burn " << localBurn << " after teta";
      teta[runs-1][j - 2] /= j*sizeMCMCstring;
      //cout << teta[runs-1][j - localBurn - 1] << "\n";
    }
#endif
     

#ifdef HALF_GELMAN
    for(unsigned long j = 1; j <= (sampleSize + 1)/constBurn; j++){
      s_2[runs-1][j - 1] = 0;
      for(unsigned long i = 0; i < sizeMCMCstring; i++){
	tempVal = MCMCstring[i]->getVarianceValue(j*constBurn, j*constBurn/2 + 1, teta[runs-1][j - 1]);
	if(tempVal != -1)
	    s_2[runs-1][j - 1] += tempVal;
	  else {
	    cout << "It is an error in the variance calculation\n";
	    return;
	  }
	}
	s_2[runs-1][j - 1] /= j*sizeMCMCstring*constBurn/2 - 1;
	//cout << s_2[runs-1][j - 1] << " for runs=" << runs << " and for j=" << j << "\n";
      }
#else
      for(unsigned long j = 2; j <= sampleSize + 1; j++){
	s_2[runs-1][j-2]= 0;
	for(unsigned long i = 0; i < sizeMCMCstring; i++){
	  tempVal = MCMCstring[i]->getVarianceValue(j,1, teta[runs-1][j - 2]);
	  if(tempVal != -1)
	    s_2[runs-1][j - 2] += tempVal;
	  else {
	    cout << "It is an error in the variance calculation\n";
	    return;
	  }
	}
	s_2[runs-1][j - 2] /= j*sizeMCMCstring - 1;
	//cout << s_2[runs-1][j - localBurn - 1] << " for runs=" << runs 
	//   << " and for j=" << j << "and local burn " << localBurn << "\n";
      }
#endif
}

double binaryPopulation::varianceCalculation(unsigned long mySamplesize){
  double teta_i= 0;
  double B = 0;

  for(unsigned long i = 0; i < runs; i++){
      teta_i += teta[i][mySamplesize];
    //cout << teta[i][mySamplesize] << " for i=" << i << " and for mySamplesize=" << mySamplesize << "\n";
 }
  teta_i /= (double)runs;
  
  
  //cout << "teta = "<<teta << " teta=";
  for(unsigned long i = 0; i < runs; i++){
    //cout << "("<<i<<","<<teta[i]<<")\t";
    B += (teta[i][mySamplesize]-teta_i)*(teta[i][mySamplesize]-teta_i);
  }
  B /= (double)runs -1;
  //cout << "B=" << B << "\t";
   
  return B;
}

double binaryPopulation::meanCalculation(unsigned long mySamplesize){
  double W = 0;

  for(unsigned long i = 0; i < runs; i++){
    W += s_2[i][mySamplesize];
    //cout << s_2[i][mySamplesize] << " for i=" << i << " and for mySamplesize=" << mySamplesize << "\n";
  }  
  W /= (double) runs;
  //cout << "W=" << W << "\n";

  return W;
}

void binaryPopulation::printReductionGNUFILE(){
   ofstream myReduction(reductionData);
   if(!myReduction.is_open()){
     cout << "Error opening file reduction in convergence.cpp:122 \n";
     exit(1);
   }
  //cout << "Calculeaza " << mySamplesize  << "\n";

  double temp = 0;

#ifdef HALF_GELMAN
   cout << "Take only half \n";
   for(unsigned long i = 1; i <= (sampleSize + 1)/constBurn; i ++ ){
     double variance = varianceCalculation(i -1);
     double mean = meanCalculation(i -1);

     temp = (double) (i*constBurn*sizeMCMCstring/2 -1)/(double) (i*constBurn*sizeMCMCstring/2);
     //cout << "variance= "<< variance << " mean=" << mean <<" \t" ;  
     if(variance >= 0 && mean > 0)
       temp += (runs + 1) * variance / (runs * mean);
     else {
       cout << "Error in reduction \n";
     }
     //temp /= mean;
     temp = pow(temp*(sizeMCMCstring*sizeGENOME + 3)/(sizeMCMCstring*sizeGENOME + 1),0.5);
     //cout << "reduction " << temp << "\n";
     
     myReduction << i*constBurn <<"\t"<<temp<<"\n";
   }
#else
   cout << "Take all \n";
   for(unsigned long i = 1; i <= sampleSize; i++){
     double variance = varianceCalculation(i - 1);
     double mean = meanCalculation(i - 1);
     
     temp = (double) (i*sizeMCMCstring -1)/(double) (i*sizeMCMCstring);
     //cout << "variance= "<< variance << " mean=" << mean <<" \t" ;  
     if(variance >= 0 && mean > 0)
       temp += (runs + 1) * variance / (runs * mean);
     else {
       cout << "Error in reduction \n";
     }
     //temp /= mean;
     temp = pow(temp*(sizeMCMCstring*sizeGENOME + 3)/(sizeMCMCstring*sizeGENOME + 1),0.5);
     //cout << "reduction " << temp << "\n";
     myReduction<<i<<"\t"<<temp<<"\n";
   }
#endif
   myReduction.close();
}
#endif //NORMAL_DISTRIBUTION

void binaryPopulation::varianceHigerMomentsUnivariateCalculation(unsigned long s, unsigned long chain){
  
    MCMCstring[chain]->getElementValue(teta[runs - 1]);

#ifdef HALF_GELMAN
  double tempVal[(sampleSize + 1)/constBurn];
  double sum = 0;

  for(unsigned long j = 1; j <= (sampleSize + 1)/constBurn; j++){
    for(unsigned long i = j*constBurn/2 + 1; i <= j*constBurn; i++)
      sum += teta[runs-1][i-1];
    tempVal[j - 1] = 2*sum/(j*constBurn);
  }

  double tempValT = 0;
  for(unsigned long j = 1; j <= (sampleSize + 1)/constBurn; j++){
    tempValT = MCMCstring[chain]->getVarianceValue(j*constBurn, j*constBurn/2 + 1, tempVal[j - 1],s);
    if(tempValT != -1)
      s_2[runs-1][j - 1] = tempValT;
    else {
      cout << "It is an error in the variance calculation\n";
      return;
    }
    s_2[runs-1][j - 1] /= j*constBurn/2 - 1;
    //cout << s_2[runs-1][j - 1] << " for runs=" << runs << " and for j=" << j << "\n";
  }
#else
  double tempVal[sampleSize  +1];
  double sum = teta[runs-1][0];

  for(unsigned long j = 2; j <= sampleSize + 1; j++){
      sum += teta[runs-1][j-1];
    tempVal[j - 1] = sum/(j-1);
  }

  double tempValT = 0;
    for(unsigned long j = 2; j <= sampleSize + 1; j++){
      tempValT = MCMCstring[chain]->getVarianceValue(j,1, tempVal[j - 2],s);
      if(tempValT != -1)
	s_2[runs-1][j - 2] = tempValT;
      else {
	cout << "It is an error in the variance calculation\n";
	return;
      }
      s_2[runs-1][j - 2] /= j - 2;
    }
#endif
}

double binaryPopulation::varianceCalculation(unsigned long mySamplesize, unsigned long s){
  double teta_i= 0;
  double B = 0;
  
#ifdef HALF_GELMAN
  for(unsigned long i = 0; i < runs; i++){
    double tempTeta = 0;

    for(unsigned long j = (mySamplesize+1)*constBurn/2 + 1; j <= (mySamplesize + 1)*constBurn; j++)
      tempTeta += teta[i][j - 1];

    teta_i += tempTeta * 2/((mySamplesize + 1)*constBurn);

  }
  teta_i /= (double)runs;

  //cout << "teta_i = "<<teta_i << " \n";
  for(unsigned long i = 0; i < runs; i++)
    for(unsigned long j = (mySamplesize+1)*constBurn/2 + 1; j <= (mySamplesize + 1)*constBurn; j++)
      if(teta[i][j-1] < teta_i) 
	B += -pow(teta[i][j-1]-teta_i, s);
      else B += pow(teta[i][j-1] - teta_i, s);
  
  B /= (double)(mySamplesize + 1)* constBurn * nr_runs/2 - 1;
#else
  for(unsigned long i = 0; i < runs; i++){
    double tempTeta = 0;
    for(unsigned long j = 0; j < mySamplesize; j++)
      tempTeta += teta[i][j];
    teta_i += tempTeta / (double)mySamplesize;
  }
  teta_i /= (double)runs;

  //cout << "teta_i = "<<teta_i << " \n";
  for(unsigned long i = 0; i < runs; i++)
    for(unsigned long j = 0; j < mySamplesize; j++)
      if(teta[i][j] < teta_i) 
	B += -pow(teta[i][j]-teta_i, s);
      else B += pow(teta[i][j] - teta_i, s);

  B /= (double)mySamplesize * runs - 1;
#endif

  //cout << "B=" << B << "\t";
  
  return B;
}

double binaryPopulation::meanCalculation(unsigned long mySamplesize, unsigned long s){
  double W = 0;

  for(unsigned long i = 0; i < runs; i++){
    W += s_2[i][mySamplesize];
    //cout << s_2[i][mySamplesize] << " for i=" << i << " and for mySamplesize=" << mySamplesize << "\n";
  }  
  W /= (double) runs;
  //cout << "W=" << W << "\n";
  
  return W;
}
  
void binaryPopulation::printReductionGNUFILE(unsigned long s){
   ofstream myReduction(reductionData);
   if(!myReduction.is_open()){
     cout << "Error opening file reduction in convergence.cpp:283 \n";
     exit(1);
   }

  double temp = 0;

#ifdef HALF_GELMAN
   cout << "Take only half \n";
   for(unsigned long i = 1; i <= (sampleSize + 1)/constBurn; i = i+SEE_SAMPLE ){
     double variance = varianceCalculation(i - 1, s);
     double mean = meanCalculation(i -1, s);

     if(variance >= 0 && mean > 0)
       temp = variance / mean;
     else {
       cout << "Error in reduction \n";
     }
     temp = pow(temp,1.0/(double)s);
     //cout << "reduction " << temp << "\n";
     
     myReduction << i*constBurn <<"\t"<<temp<<"\n";
   }
#else
   cout << "Take all \n";
   for(unsigned long i = 1; i <= sampleSize; i++){
     double variance = varianceCalculation(i - 1, s);
     double mean = meanCalculation(i - 1, s);
     
     if(variance >= 0 && mean > 0)
       temp = variance / mean;
     else {
       cout << "Error in reduction \n";
     }
     temp = pow(temp,1.0/(double)s);
     //cout << "reduction " << temp << "\n";
     myReduction<<i<<"\t"<<temp<<"\n";
   }
#endif
   myReduction.close();
}

/*---Multivariate Calculation*/

void binaryPopulation::varianceMultivariateCalculation(){
    for(unsigned long i = 0; i < sizeMCMCstring; i++)
      MCMCstring[i]->getElementValue(tetaM[runs - 1][i]);
}

double binaryPopulation::varianceMultivariateCalculation(unsigned long mySamplesize, double** B){
  double meanM[sizeMCMCstring];
  double meanL[sizeMCMCstring][nr_runs];

  for(unsigned long i = 0; i < sizeMCMCstring; i++){
    meanM[i] = 0;
    for(unsigned long r = 0; r < nr_runs; r++ ){
#ifdef HALF_GELMAN
    double sum = 0;
    for(unsigned long j = (mySamplesize + 1)*constBurn/2 + 1; j <= (mySamplesize+1)*constBurn; j++)
      sum += tetaM[r][i][j-1];

    meanM[i] += sum/((double) (mySamplesize + 1)*constBurn/2);
    meanL[i][r] = sum/((double)(mySamplesize + 1)*constBurn/2);
#else
    double sum = 0;
    for(unsigned long j = 0; j <= mySamplesize; j++)
      sum += tetaM[r][i][j];

    meanM[i] += sum/(double) mySamplesize;
    meanL[i][r] = sum/(double)mySamplesize;
#endif    
    }

    meanM[i] /= (double)(runs -1);
  }

  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    for(unsigned long j = i; j < sizeMCMCstring; j++){
      B[i][j] = 0;

      for(unsigned long r = 0; r < nr_runs; r++)
	B[i][j] += (meanL[i][r] - meanM[i])*(meanL[j][r] - meanM[j]);
      
      B[i][j] /= (double)(nr_runs - 1);
      if (i != j) B[j][i] = B[i][j];
    }

  //delete[] meanL;
//  delete[] meanM;
  return 0;
}

double binaryPopulation::meanMultivariateCalculation(unsigned long mySamplesize, double** W){
  double meanL[sizeMCMCstring][nr_runs];

  for(unsigned long i = 0; i < sizeMCMCstring; i++){
    for(unsigned long r = 0; r < nr_runs; r++ ){
#ifdef HALF_GELMAN
    double sum = 0;
    for(unsigned long j = (mySamplesize + 1)*constBurn/2 + 1; j <= (mySamplesize+1)*constBurn; j++)
      sum += tetaM[r][i][j - 1];

    meanL[i][r] = sum/((double)(mySamplesize + 1)*constBurn/2);
    //cout << meanL[i][r] << "," << i << "," << r << "\n";
#else
    double sum = 0;
    for(unsigned long j = 0; j <= mySamplesize; j++)
      sum += tetaM[r][i][j];
    meanL[i][r] = sum/(double)mySamplesize;
#endif    
    }
  }

  for(unsigned long i = 0; i < sizeMCMCstring; i++)
    for(unsigned long j = i; j < sizeMCMCstring; j++){
      W[i][j] = 0;

#ifdef HALF_GELMAN
      for(unsigned long r = 0; r < nr_runs; r++)
	for(unsigned long k = (mySamplesize + 1)*constBurn/2 + 1; k <= (mySamplesize+1)*constBurn; k++ )
	  W[i][j] += (tetaM[r][i][k-1] - meanL[i][r])*(tetaM[r][j][k-1] - meanL[j][r]); 
	
      //cout << W[i][j] << "," << i << "," << j << "," << (mySamplesize+1)*constBurn/2 - 1 << "\n";

      W[i][j] /= (double)(nr_runs*((mySamplesize+1)*constBurn/2 - 1));

      //cout << W[i][j] << "," << i << "," << j << "\n";
#else
      for(unsigned long r = 0; r < nr_runs; r++)
	for(unsigned long k = 1; k <= mySamplesize; k++)
	  W[i][j] += (tetaM[r][i][k] - meanL[i][r])*(tetaM[r][j][k] - meanL[j][r]); 

      W[i][j] /= (double)(nr_runs*(mySamplesize - 1));
#endif


      if(i != j) W[j][i] = W[i][j];
    }

  //delete[] meanL;
  return 0;
}

void binaryPopulation::printReductionMultivariateGNUFILE(){
  ofstream myReduction(reductionMultivariateData);
  if(!myReduction.is_open()){
    cout << "Error opening file reduction in convergence.cpp:431 \n";
    exit(1);
  }

  double** B;
  double** W;
  B = new double*[sizeMCMCstring];
  W = new double*[sizeMCMCstring];
  for(unsigned long j = 0; j < sizeMCMCstring; j++){
    B[j] = new double[sizeMCMCstring]; 
    W[j] = new double[sizeMCMCstring];
  }
//  double B[sizeMCMCstring][sizeMCMCstring];
//  double W[sizeMCMCstring][sizeMCMCstring];
  
#ifdef HALF_GELMAN
  cout << "Take only half \n";
  for(unsigned long i = 1; i <= (sampleSize + 1)/constBurn; i ++ ){
    varianceMultivariateCalculation(i - 1, B);
    meanMultivariateCalculation(i -1, W);
    
     //scrierea W in fisier
    myReduction << "W" << i <<"=[";
    for(unsigned long ci = 0; ci < sizeMCMCstring; ci++){
      for(unsigned long g = 0; g < sizeMCMCstring - 1; g++)
	myReduction << W[ci][g] << ",";
      if(ci != sizeMCMCstring - 1){
	myReduction << W[ci][sizeMCMCstring-1] << ";";
	//cout << W[ci][sizeMCMCstring-1] << "," << ci << "," << sizeMCMCstring -1 << "\n";
      }
      else myReduction << W[ci][sizeMCMCstring-1] << "];";
    }
    
    //scrierea B in fisier
    myReduction << "\n" << "B" << i <<"=[";
    for(unsigned long ci = 0; ci < sizeMCMCstring; ci++){
      for(unsigned long g = 0; g < sizeMCMCstring - 1; g++)
	myReduction << B[ci][g] <<",";
      if(ci != sizeMCMCstring - 1) 
	myReduction << B[ci][sizeMCMCstring-1] << ";";
      else myReduction << B[ci][sizeMCMCstring-1] << "];";
    }
    
    //scrierea operatiilor in fisier
    myReduction <<"Lmax=sort(abs(eig(inverse(W"<<i <<")*B"<<i<<")));\n";
    myReduction <<"l1= Lmax/max(Lmax);\n";
    myReduction <<"l= l1(" << sizeMCMCstring - 2 << ");\n";
    myReduction << "Rp(" << i << ")=" << (i*constBurn/2 - 1)/(double)(i*constBurn/2) 
		<< "+" << (nr_runs + 1)/(double)nr_runs << "*l;";   
  }
  myReduction << "\nsave \"gnuREDMultiFILE" << specific << " Rp \n";
 
#else
  cout << "Take all \n";
   for(unsigned long i = 1; i <= sampleSize+1; i=i+SEE_SAMPLE){
     varianceMultivariateCalculation(i - 1, B);
     meanMultivariateCalculation(i - 1, W);

     //scrierea W in fisier
     myReduction << "W=[";
     for(unsigned long c1 = 0; c1 < sizeMCMCstring; c1++){
       for(unsigned long g = 0; g < sizeMCMCstring - 1; g++)
	 myReduction << W[c1][g] <<",";
       if(c1 < sizeMCMCstring - 1) 
	 myReduction << W[c1][sizeMCMCstring-1] << ";";
       else myReduction << W[c1][sizeMCMCstring-1] << "];\n";
     }

     //scrierea B in fisier
     myReduction << "B=[";
     for(unsigned long c1 = 0; c1 < sizeMCMCstring; c1++){
       for(unsigned long g = 0; g < sizeMCMCstring - 1; g++)
	 myReduction << B[c1][g] <<",";
       if(c1 < sizeMCMCstring - 1) 
	 myReduction << B[c1][sizeMCMCstring-1] << ";";
       else myReduction << B[c1][sizeMCMCstring-1] << "];\n";
     }
     
     //scrierea operatiilor in fisier
    myReduction <<"Lmax=sort(abs(eig(inverse(W)*B)));\n";
    myReduction <<"l1= Lmax/max(Lmax);\n";
    myReduction <<"l= l1(" << sizeMCMCstring - 1 << ");\n";
    myReduction << "Rp(" << i << ")=" << (i - 1)/(double)i 
		<< "+" << (nr_runs + 1)/(double)nr_runs << "*l;\n";   
  }
  myReduction << "\nsave \"gnuREDMultiFILE" << specific << "\" Rp \n";

#endif
  delete[] B;
   delete[] W;
   myReduction.close();
}

#endif //Gelman and Rubin convergence

#ifdef HISTOGRAM
#ifdef KULLBACK_INFORMATION
#ifdef EXPAND_KULLBACK_INFORMATION
//calculeaza distributia ideala
void binaryPopulation::expandTrueUnivariateDistribution(){

//1. initialize the data
  ofstream distributionFile;
  distributionFile.open(distributionFILE);
  if(!distributionFile.is_open()){
    cout << "Error in opening for writting the distribution file " 
	 << distributionFILE << "in convergence.cpp:529\n"; 
    exit(1);
  }

#ifndef MIXTURE_REALS  
#ifdef RECOMB_MUT
  double distance[diversityMeasure];
  double count[diversityMeasure];
  for(unsigned long i = 0; i < diversityMeasure; i++){
    distance[i] = 0;
    count[i] = 0;
  }
#endif

  const unsigned long number = (unsigned long) pow(2,sizeGENOME); // dimensiunea spatiului de cautare
  //double diversityTrueDistribution[diversityMeasure];
  unsigned long tempExpGenome[sizeGENOME];
  for(unsigned long i = 0; i < sizeGENOME; i++)
    tempExpGenome[i] = 0;
  for(unsigned long i = 0; i < diversityMeasure; i++){
    diversityTrueDistribution[i] = 0;
  }

#else
  unsigned long number = (unsigned long) pow(2.0,(double)sizeGENOME*BLOCKsize); // precizia spatiului de cautare
  //double diversityTrueDistribution[diversityMeasure];
  unsigned long tempExpGenome[sizeGENOME*BLOCKsize];
  double reals[sizeGENOME];
#ifdef HISTOGRAM
  // for(unsigned long i = 0; i < diversityMeasure; i++)
  //    for(unsigned long j = 0; j < diversityMeasure; j++){
//	  diversityTrueDistribution[i][j] = 0;
  //     }
#endif // histogram
  //cout << "sizetemp" << sizeGENOME*BLOCKsize;
  for(unsigned long i = 0; i < sizeGENOME*BLOCKsize; i++){
    tempExpGenome[i] = 0;
    //cout << " (i " << i << "," << tempExpGenome[i] << ")\t";
  }
  //cout << "\n";
#endif

  double meanFitness = 0;

  //2. expand the distribution
  unsigned long long impartitor;
  unsigned long rest;
  long double total = 0;

  for(unsigned long long i = 0; i < number; i++){
    impartitor = i;
    unsigned long long j = 0;
    while(impartitor != 0){
      rest = impartitor % 2;
      tempExpGenome[j++] = rest;
      impartitor = (impartitor-rest) / 2;
    }

#ifndef MIXTURE_REALS
#ifndef FLOOR 
    unsigned long value = (unsigned long) MCMCstring[0]->fitness(tempExpGenome,sizeGENOME);
#else
    unsigned long value = floor(MCMCstring[0]->fitness(tempExpGenome,sizeGENOME));
#endif
    if(value > AwayFromZero)
	total += value;
    else total += AwayFromZero;

    diversityTrueDistribution[value] ++;

#ifdef TWO_ATTRACTORS_FUNCTION
    unsigned long ones = 0;
    for(unsigned long k =0; k < sizeGENOME; k++)
	ones += tempExpGenome[k];
    histogram[ones]++;
    //cout << "histohram" << histogram[ones] << "  with" << ones << "\n";
#else
#ifdef BINOMIAL
#ifdef RECOMB_MUT
    unsigned long div = (unsigned long) HammingDistance(tempExpGenome);
    distance[div] += value;
    count[div]++;
#endif
/*    if(value % 2 == 0)
	histogram[BLOCKsize/2-value/2]++;
    else
	histogram[BLOCKsize/2+(value+1)/2]++;
*/
#endif
#endif // two attractor function
#else //mixtire reals
    
    //////////////////test convert reals
    //cout << "sizetemp" << sizeGENOME*BLOCKsize;
    //for(unsigned long i = 0; i < sizeGENOME*BLOCKsize; i++){
      //  tempExpGenome[i] = 0;
    //cout << " (i " << i << "," << tempExpGenome[i] << ")\t";
    //}
    //cout << "\n";
    //cout << "realtemp" << sizeGENOME;
    //for(unsigned long i = 0; i < sizeGENOME; i++){
      //  tempExpGenome[i] = 0;
    //cout << " (i " << i << "," << reals[i] << ")\t";
    //}
    //cout << "\n";

    convertReal(reals,tempExpGenome,sizeGENOME);

    ////////////////test convertReals
    // cout << "sizetemp" << sizeGENOME*BLOCKsize;
    //for(unsigned long i = 0; i < sizeGENOME*BLOCKsize; i++){
      //  tempExpGenome[i] = 0;
    //cout << " (i " << i << "," << tempExpGenome[i] << ")\t";
    /// }
    //cout << "\n";
    //cout << "realtemp" << sizeGENOME;
    //for(unsigned long i = 0; i < sizeGENOME; i++){
    //  //  tempExpGenome[i] = 0;
    //cout << " (i " << i << "," << reals[i] << ")\t";
    //}
    //cout << "\n";

    double value = MCMCstring[0]->fitness(reals,sizeGENOME);
    if(value > AwayFromZero)
	total += value;
    else total += AwayFromZero;

    unsigned long indexIK = 0;
#ifndef BIVARIATE
    for(int ik = 0; ik < sizeGENOME; ik++){
      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
    }
#else //bivariate
    for(int ik = 0; ik < sizeGENOME; ik++){
      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
    }
#endif //bivariate

    diversityTrueDistribution[indexIK]++;
    trueDistribution[indexIK] += value;
#endif //mixture reals

   //cout << "Total value for i=" << i << " is " << total << " with genome ";
    //for(j = 0; j < sizeGENOME; j++)
    //  cout << tempExpGenome[j];
    //cout << "\n"; 
  }   

 //3. write the histogram for the discrete function in a separate file 
  //3'. for the reals function, write the histogram for ploting
#ifndef MIXTURE_REALS 
#ifdef TWO_ATTRACTORS_FUNCTION
  ofstream myHistogram(histogramFILE);		
  if(!myHistogram.is_open())
      cout << "Error opening file " << histogramFILE << " convergence.cpp:602 \n";

  //cout << "Diversity true distribution \n";
  for(unsigned long i =0; i< sizeGENOME + 1; i++){
    myHistogram << i <<"\t" << histogram[i]<<"\n";
    // cout << i << "\t" <<histogram[i]<<"\n";
  }

  myHistogram.flush();
  myHistogram.close();
#else
#ifdef BINOMIAL
#ifdef RECOMB_MUT
  ofstream myDistance(histogramFILE);		
  if(!myDistance.is_open())
      cout << "Error opening file " << histogramFILE << "convergence.cpp:616\n";

  for(unsigned long i =0; i< diversityMeasure; i++){
    myDistance << i <<"\t" << distance[i]/count[i]<<"\n";
    // cout << i << "\t" <<histogram[i]<<"\n";
  }

  myDistance.flush();
  myDistance.close();
#endif //recomb mut
/*  ofstream myHistogram(histogramFILE);		
  if(!myHistogram.is_open())
    cout << "Error opening file \n";

  //cout << "Diversity true distribution \n";
  for(unsigned long i =0; i< sizeGENOME + 1; i++){
    myHistogram << i <<"\t" << histogram[i]<<"\n";
    // cout << i << "\t" <<histogram[i]<<"\n";
  }
  myHistogram.flush();
  myHistogram.close();
*/
#endif //binomial
#endif //two attractor
#else
#ifdef HISTOGRAM
  ofstream myDistance(histogramFILE);		
  if(!myDistance.is_open())
      cout << "Error opening file " << histogramFILE << "convergence.cpp:616\n";

  for(unsigned long i =0; i < diversityMeasure; i++){
    myDistance << i <<"\t";
    myDistance << trueDistribution[i]/diversityTrueDistribution[i] << "n";
    // cout << i << "\t" <<histogram[i]<<"\n";
  }

  myDistance.flush();
  myDistance.close();
 
#endif
#endif //mixture reals

////////////////////////
////5. compute the mean of the normalized distribution; how much bigger, how mach smaller then the mean
///////////////////////
  double mean = total/number;
  double mean_unnorm = 0;
  double lessMean = 0, greaterMean = 0;
  double lessNoMean = 0, greaterNoMean = 0;
//calcularea mean
  for(unsigned long long i = 0; i < number; i++){
    impartitor = i;
    unsigned long long j = 0;
    while(impartitor != 0){
      rest = impartitor % 2;
      tempExpGenome[j++] = rest;
      impartitor = (impartitor-rest) / 2;
    }

#ifndef MIXTURE_REALS
#ifndef FLOOR 
    unsigned long value = (unsigned long) MCMCstring[0]->fitness(tempExpGenome,sizeGENOME);
#else
    unsigned long value = floor(MCMCstring[0]->fitness(tempExpGenome,sizeGENOME));
#endif
#else	
    convertReal(reals,tempExpGenome,sizeGENOME);	
    double value = MCMCstring[0]->fitness(reals,sizeGENOME);
#endif //MIXTURE_REALS
    if(value != 0){
	//double fitnessValue = MCMCstring[0]->fitness(tempExpGenome,sizeGENOME);
	if(value < (double)total/number){
	    mean_unnorm = mean_unnorm - 1;
	    lessNoMean++;
	    lessMean += value;
	}
	else if (value > (double)total/number){ 
		 mean_unnorm=mean_unnorm + 1;
		 greaterNoMean++;
		 greaterMean += value;
	}
    }
    else {
	mean_unnorm = mean_unnorm - 1;
	lessNoMean++;
	lessMean += AwayFromZero;
    }
    //cout << "Total value for i=" << i << " is " << total << " with genome ";
    //for(j = 0; j < sizeGENOME; j++)
    //  cout << tempExpGenome[j];
    //cout << "\n"; 

  }   

  //distributionFile << "Total value: " << total << "\t" 
  //<< "Mean: " << mean << "\t" << "total invid: " << number << "\n\n";
  
  distributionFile << total << "\n"; // << total << "\n";
  //distributionFile << total/number << "\n";
  distributionFile << mean << "\n";
 //cout << "Total distribution " << total << "mean" << "\n";

//6. print the histogram: for discrete spaces -> based on the number of individuals for a fitness value
// for reals spaces -> the discretized histogram 
#ifndef MIXTURE_REALS
  for(unsigned long i = 0; i < diversityMeasure; i++){
    if (diversityTrueDistribution[i] != 0)
	//distributionFile << ((double)diversityTrueDistribution[i]/number) << "\n";
	if(i != 0)
	    distributionFile << i/total << "\n";
		//((double)i/total) << "\n";
	else  distributionFile << AwayFromZero/total << "\n";
    else distributionFile << 0.0 << "\n";
  }
#else
#ifdef HISTOGRAM
  for(unsigned long i = 0; i < diversityMeasure; i++){
    trueDistribution[i] = trueDistribution[i]/diversityTrueDistribution[i];
    //diversityTrueDistribution[i][j] = trueDistribution[i][j];
    distributionFile << trueDistribution[i] << "\n";
  }
#endif
#endif //MIXTURE_REALS

//5'. print what we found at 5; these data are not read in a file
  distributionFile << "\n" << total/(double)number << "\n" << mean_unnorm/(double)number << "\n";
  distributionFile << "\n" << lessMean/(double)lessNoMean << "\t" << lessMean << "\n";
  distributionFile << "\n" << greaterMean/(double)greaterNoMean << "\t" << greaterMean << "\n";
  distributionFile << "\n" << greaterMean/(double)greaterNoMean - lessMean/(double)lessNoMean << 
    "\t" << (greaterMean - lessMean)/(greaterNoMean + lessNoMean) <<  "\t" << 
    (greaterMean - lessMean)/(total) << "\t" << (greaterNoMean - lessNoMean)/(number) << "\n";

  distributionFile.close();

//7. print number of individuals in each bin of the histogram
#ifndef MIXTURE_REALS
  ofstream mySampleFile(sampleFile);		
  if(!mySampleFile.is_open())
      cout << "Error opening file "<< sampleFile << "convergence.cpp:688 \n";

  //cout << "Diversity true distribution \n";
  for(unsigned long i =0; i< diversityMeasure; i++){
    mySampleFile <<diversityTrueDistribution[i]<<"\n";
    //cout << i << "\t" <<diversityTrueDistribution[i]<<"\n";
  }

  mySampleFile.close();
#else
#ifdef HISTOGRAM
  ofstream mySampleFile(sampleFile);		
  if(!mySampleFile.is_open())
      cout << "Error opening file "<< sampleFile << "convergence.cpp:790 \n";

  //cout << "Diversity true distribution \n";
  for(unsigned long i =0; i< diversityMeasure; i++)
    mySampleFile << diversityTrueDistribution[i]<<"\n";
  //cout << i << "\t" <<diversityTrueDistribution[i]<<"\n";
  

  mySampleFile.close();
#endif //HISTOGRAM
#endif //MIXTURE_REALS
}
/*
void binaryPopulation::expandTrueMultivariateDistribution(){
   const unsigned long number = (unsigned long) pow(2.0,(double)sizeGENOME*sizeMCMCstring); 
// dimensiunea spatiului de cautar
   
  doubleList* diversityTrueDistribution = new doubleList();

  unsigned long tempPopulation[sizeGENOME*sizeMCMCstring];
  for(unsigned long i = 0; i < sizeGENOME*sizeMCMCstring; i++)
    tempPopulation[i] = 0;

  unsigned long tempGENOME[sizeGENOME];

  unsigned long impartitor;
  unsigned long rest;
  long double total = 0;

  for(unsigned long i = 0; i < number; i++){
    impartitor = i;
    unsigned long j = 0;
    while(impartitor != 0){
      rest = impartitor % 2;
      tempPopulation[j++] = rest;
      impartitor = (impartitor-rest) / 2;
    }

    unsigned long value = 1;

    for(unsigned long k = 0; k < sizeMCMCstring; k++){
      
      for(unsigned long t = 0; t < sizeGENOME; t++)
	tempGENOME[t] = tempPopulation[k*sizeGENOME+t];
      
      value *= (unsigned long) MCMCstring[k]->fitness(tempGENOME,sizeGENOME);
    }

    diversityTrueDistribution -> addElement((double)value);

    total += value;
    mean = total /(double) number;
  }   

  //cout << "Total distribution " << total << "\n";
  ofstream distributionFile;
  distributionFile.open(distributionFILE);
  if(!distributionFile.is_open()){
    cout << "Error in opening for writting the distribution file " << distributionFILE 
	 << " in convergence.cpp:630\n";
    exit(1);
  }

  for(unsigned long i = 0; i < diversityTrueDistribution->getSize(); i++){
    distributionFile << diversityTrueDistribution->getObject(i) << "\t" 
		     << ((double)diversityTrueDistribution->getObject(i)/total) << "\t"
		     << (double)diversityTrueDistribution->getNumbers(i) << "\n";
  } 

  distributionFile.close();
  }*/
#endif //EXPAND_KULLBACK_INFORMATION

//univariate
//compute the difference between the true distribution and the actual distribution	
double binaryPopulation::mixingTimeUnivariateCalculation(){
    //   cout << "compute difference between individuals mixingTimeUnivariateCalculation \n";
 double tempValue = 0;
  double difDistr = 0;
#ifndef MIXTURE_REALS
  unsigned long tempBuffer[sizeGENOME+1];
#else
  unsigned long tempBuffer[sizeGENOME*BLOCKsize];
  double reals[sizeGENOME];
#ifdef HISTOGRAM
  double actualDistrib[diversityMeasure];
  double actualFitness[diversityMeasure];
  for(unsigned long i = 0; i < diversityMeasure; i++){
      actualDistrib[i] = 0;
      actualFitness[i] = 0;
  }

#ifdef BIVARIATE
  nrHighIndivGlobal = 0;
  double actualHighDistrib[diversityMeasure];
  double actualHighFitness[diversityMeasure];
  for(unsigned long i = 0; i < diversityMeasure; i++){
      actualHighDistrib[i] = 0;
      actualHighFitness[i] = 0;
  }
#endif //bivariate

#endif //histogram  
#endif //MIXTURE_REALS

//1. if no hashtable exists 
#ifdef FILE_BOOKEEPING
  tempValue = 0;
  char buffer[sizeGENOME+1];
  if (diversityMeasure == 0) return -1;
  diversityFile.seekg(0,ios::beg);
  unsigned long lenghtFI = diversityFile.tellg();
  diversityFile.seekg(0,ios::end);
  unsigned long lenghtFE = diversityFile.tellg();

  // calcularea constantei de normalizare pentru statiu expoarat
  //Kullback information
  diversityFile.seekg(0,ios::beg);
  while((lenghtFE-lenghtFI) > sizeGENOME){
    diversityFile.getline(buffer,100);  
   for(unsigned long i = 0; i<sizeGENOME; i++)
      if(buffer[i] == '0') tempBuffer[i] = 0;
      else tempBuffer[i] = 1;

    tempValue += MCMCstring[0]->fitness(tempBuffer,sizeGENOME);

    //calcularea valori approximative din rulari
    double valueT = (double)myHashtables->check_individual(tempBuffer,sizeGENOME)/(double)(sizeMCMCstring*((generation - burn_in)/interSamplesize+ 1));
    //verification += valueT;

    double valueA = (double)trueDistribution[(unsigned long) MCMCstring[0]->fitness(tempBuffer,sizeGENOME)];

    //calcularea KL
    if(valueT  > AwayFromZero && valueA > AwayFromZero)
    Kullback[unsigned long ((generation - burn_in)/interSamplesize)] += - valueT * log((double)valueT/(double)valueA);
    
    //for(int i = 0; i < sizeGENOME; i++)
    // cout << tempBuffer[i]; 
    //cout<< "= "<< buffer << " frecventa=(" << valueT <<","<<myHashtables->check_individual(tempBuffer,sizeGENOME) 
    //<< "," << sizeMCMCstring*((generation-burn_in)/interSamplesize + 1)<< ")\t distributie=(" << valueA << "," 
    //<< MCMCstring[0]->fitness(tempBuffer,sizeGENOME) << ")\t KL[" << int ((generation -burn_in)/interSamplesize) <<"]=" 
    //<< Kullback[int ((generation -burn_in)/interSamplesize)] << "from " 
    //<< (- valueT * log((double)valueT/(double)valueA)) <<"\n";

    diversityFile.seekg(0,ios::cur);
    lenghtFI = diversityFile.tellg();
  }
  cout << "total= " << total << "total current " << tempValue << "\n";

  //calcularea diferentei dintre distributtii
  difDistr = 0;
  diversityFile.seekg(0,ios::beg);
  while((lenghtFE-lenghtFI) > sizeGENOME){
    diversityFile.getline(buffer,100);  
    unsigned long tempBuffer[sizeGENOME];
    for(unsigned long i = 0; i<sizeGENOME; i++)
      if(buffer[i] == '0') tempBuffer[i] = 0;
      else tempBuffer[i] = 1;

    double tempFitnessValue = MCMCstring[0]->fitness(tempBuffer,sizeGENOME);
    difDistr += (double) tempFitnessValue * ((double)1/(double)total - (double)1/(double)tempValue) * (log((double)tempFitnessValue/(double)total) - log((double)tempFitnessValue/(double)tempFitnessValue));
    
    diversityFile.seekg(0,ios::cur);
    lenghtFI = diversityFile.tellg();
  }

#else 
//2. fara scriere in fisier; cu hashtable

  //2.0 initializarea datelor
   tempValue = (generation+1) * sizeMCMCstring;

#ifdef THRESHOLD_DISTRIBUTION
#ifndef MIXTURE_REALS
   long double totalNr = pow(2,sizeGENOME); // - tempValue;
#else
   long double totalNr = pow(2,sizeGENOME*BLOCKsize); // - tempValue;
#endif
   unsigned long mari = 0;
   unsigned long limitaMari = 0;
#endif

  difDistr = 0;

  cout << "calcularea KL distribution \n";

  //2.1 citirea din hashtable a distributiei care a ramas
  //2.1.1 initilaizarea interatorului

#ifndef MIXTURE_REALS
  double tempNrValue = myHashtables -> iteratorInit(tempBuffer,sizeGENOME);
#else
  //cout << "first elelemnt in Hashtable \n";
  double tempNrValue = myHashtables -> iteratorInit(tempBuffer,sizeGENOME*BLOCKsize);
  //cout << " value" << tempNrValue << " for genome ";
  //for(unsigned long inter = 0; inter < sizeGENOME*BLOCKsize; inter++)
  //    cout << tempBuffer[inter];
  //cout << "\n";
#endif //MIXTURE_REALS

/// 2.1.1.1 if exist something in the hashtable
  if(tempNrValue == -1){
    cout << " it is nothing in the hashtable \n";
  }
  else{ 
      //primul element din Hashtable
      if(tempValue == 0) {
	difDistr = 0;
	cout << "nothing left there \n";
      }
      else{
	//cout << "first element in hashtable \n";
#ifndef MIXTURE_REALS
#ifndef FLOOR
	  unsigned long tempFitness = MCMCstring[0]->fitness(tempBuffer,sizeGENOME);
#else
	  unsigned long tempFitness = floor(MCMCstring[0]->fitness(tempBuffer,sizeGENOME));
#endif

#ifdef THRESHOLD_DISTRIBUTION
#ifndef HALF_MAX_FITNESS
#ifdef BINOMIAL
	  if(tempFitness > (diversityMeasure - 1) * (BLOCKsize-1)/BLOCKsize){
#else
#ifdef  TRAP_LIKE_FUNCTION
	      if(tempFitness > (diversityMeasure - 1) * b/a){
#endif
#endif //BINOMIAL
#else
		  if(tempFitness >= (diversityMeasure - 1) /2){
#endif //HALF_MAX_FITNESS

	  tempDiversity[tempFitness] = tempDiversity[tempFitness] - 1;
	  //if(tempFitness == 0) tempFitness = AwayFromZero;
 
	  //for(int i = 0; i < sizeGENOME; i++)
	  //  cout << tempBuffer[i];
	  //cout << " = " << tempNrValue/(double)tempValue << "\t" << tempNrValue << "\t" << (double)tempValue 
	  // << "\t" << trueDistribution[tempFitness] << "\n";
	  	  
#ifndef KL_DIFFERENCE	  
	  if(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue) > 0){
	      difDistr = trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue);
	  }
	  else {
	      difDistr = -(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue));
	      mari++;
	  }
#else
	  difDistr = (trueDistribution[tempFitness] - 
	   (double)tempNrValue/(double)(tempValue))*(log(trueDistribution[tempFitness]) - 
	   log((double)tempNrValue/(double)(tempValue)));
	  if(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue) < 0){
	      mari++;
	  }
#endif //KL_DIFFERENCE
	  limitaMari++;
	  }
	  else difDistr = 0;
#else // THRESHOLD_DISTRIBUTION
		  tempDiversity[tempFitness] = tempDiversity[tempFitness] - 1;
	  //if(tempFitness == 0) tempFitness = AwayFromZero;
 
	  //for(int i = 0; i < sizeGENOME; i++)
	  //  cout << tempBuffer[i];
	  //cout << " = " << tempNrValue/(double)tempValue << "\t" << tempNrValue << "\t" << (double)tempValue 
	  // << "\t" << trueDistribution[tempFitness] << "\n";
	  
#ifndef KL_DIFFERENCE	  
	  if(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue) > 0)
	      difDistr = trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue);
	      else difDistr = -(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue));
#else
	  difDistr = (trueDistribution[tempFitness] - 
	   (double)tempNrValue/(double)(tempValue))*(log(trueDistribution[tempFitness]) - 
	   log((double)tempNrValue/(double)(tempValue)));
#endif
#endif // THRESHOLD_DISTRIBUTION
#else //mixture reals
	  //cout << " arrived on the right way ..................\n";
	  convertReal(reals,tempBuffer,sizeGENOME);
	  double tempFitness = MCMCstring[0]->fitness(reals,sizeGENOME);
	  //cout << "find in hashtable \t";
	  //for(unsigned long inter = 0; inter < sizeGENOME*BLOCKsize; inter++)
	  //  cout << tempBuffer[inter];
	  //cout << "= " << tempFitness << "\n";

#ifndef HISTOGRAM
	  double trueD = (double)tempFitness/(double)total;
	  double actualD = (double)tempNrValue/(double)tempValue;
	  
#ifndef KL_DIFFERENCE
	  if(trueD > actualD)
	    difDistr = trueD - actualD;
	  else difDistr = - trueD + actualD;
#else
	  difDistr = (trueD - actualD)* (log(trueD) - log(actualD));
#endif //kl difference

#else // histogram
	  //discreteze first
	  //cout << "arived on histogram way \n";
	  unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate

	  actualDistrib[indexIK] += tempNrValue;
	  actualFitness[indexIK] += tempNrValue/(double)tempFitness;

#ifdef UNIFORM_DISTRIBUTION
	  //nrHighIndivGlobal = nrHighIndivGlobal + 1;
	    //actualHighDistrib[(unsigned long)(reals[0]*hist_points/scale)][(unsigned long)(reals[1]*hist_points/scale)]+= tempNrValue;
	    //actualHighFitness[(unsigned long)(reals[0]*hist_points/scale)][(unsigned long)(reals[1]*hist_points/scale)]+= tempNrValue/(double)tempFitness;
#else //histogram
#ifdef BIVARIATE	  
#ifdef MIXTURE_BIVARIATE
#ifdef NINE_BIVARIATE
	  if((reals[0] < sigma_1 + scale/2. && reals[0] > -sigma_1 + scale/2. && reals[1] < sigma_2  + scale/2. && reals[1] > -sigma_2 + scale/2.) ||
	     (reals[0] < sigma_1 + scale/2. + 0.375*scale && reals[0] > -sigma_1 + scale/2. + 0.375*scale && reals[1] < sigma_2  + scale/2. + 0.375*scale && reals[1] > -sigma_2 + scale/2. + 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. + 0.375*scale && reals[0] > -sigma_1 + scale/2. + 0.375*scale && reals[1] < sigma_2  + scale/2. - 0.375*scale && reals[1] > -sigma_2 + scale/2. - 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. - 0.375*scale && reals[0] > -sigma_1 + scale/2. - 0.375*scale && reals[1] < sigma_2  + scale/2. + 0.375*scale && reals[1] > -sigma_2 + scale/2. + 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. - 0.375*scale && reals[0] > -sigma_1 + scale/2. - 0.375*scale && reals[1] < sigma_2  + scale/2. - 0.375*scale && reals[1] > -sigma_2 + scale/2. - 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. + 0.375*scale && reals[0] > -sigma_1 + scale/2. + 0.375*scale && reals[1] < sigma_2  + scale/2. && reals[1] > -sigma_2 + scale/2.) ||
	     (reals[0] < sigma_1 + scale/2. && reals[0] > -sigma_1 + scale/2. && reals[1] < sigma_2  + scale/2. + 0.375*scale && reals[1] > -sigma_2 + scale/2. + 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. - 0.375*scale && reals[0] > -sigma_1 + scale/2. - 0.375*scale && reals[1] < sigma_2  + scale/2. && reals[1] > -sigma_2 + scale/2.) ||
	     (reals[0] < sigma_1 + scale/2. && reals[0] > -sigma_1 + scale/2. && reals[1] < sigma_2  + scale/2. - 0.375*scale && reals[1] > -sigma_2 + scale/2. - 0.375*scale)){
	    nrHighIndivGlobal = nrHighIndivGlobal + 1;

	    unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate
	    
	    actualHighDistrib[indexIK]+= tempNrValue;
	    actualHighFitness[indexIK]+= tempNrValue/(double)tempFitness;
	    //   cout << " fitness " << tempFitness << " nr Indiv " << tempNrValue << " for (" << reals[0] << "," << reals[1] << ") \n"; 
	  }
#else
	  //????????????????????????????????????????????????????????????????????????????
	  //write the variance in each point
	  //??????????????????? not inplemented
	  
	  double good = 1;
	  for(int iDimension = 0; iDimension < sizeGENOME; iDimension++)
	  for(int iMixture = 0; iMixture < mixture; iMixture++)
	    if(!(reals[iDimension] > MCMCstring[0]->positionMixture[iMixture].element(iDimension,0) - 
		 sqrt(MCMCstring[0]->covarianceMixture[iMixture].element(iDimension,iDimension)) 
		 && reals[iDimension] < MCMCstring[0]->positionMixture[iMixture].element(iDimension,0) + 
		 sqrt(MCMCstring[0]->covarianceMixture[iMixture].element(iDimension,iDimension))))
	      good = 0;
	  
	  if(good == 1){
	    nrHighIndivGlobal = nrHighIndivGlobal + 1;

	    unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate

	    actualHighDistrib[indexIK] += tempNrValue;
	    actualHighFitness[indexIK] += tempNrValue/(double)tempFitness;
	  }
#endif //9_bivariate
#else // mixture bivariate
#ifndef TWO_PICKS
	  if(reals[0] < sigma_1 * scale + scale/2. && reals[0] > -sigma_1 * scale + scale/2. && 
	     reals[1] < sigma_2 * scale + scale/2. && reals[1] > -sigma_2 * scale + scale/2.){
	    nrHighIndivGlobal = nrHighIndivGlobal + 1;

	    unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate

	    actualHighDistrib[indexIK]+= tempNrValue;
	    actualHighFitness[indexIK]+= tempNrValue/(double)tempFitness;
	    //cout << " fitness " << tempFitness << " nr Indiv " << tempNrValue << 
	    //  " for (" << reals[0] << "," << reals[1] << ") = " << MCMCstring[0]->fitness(reals,sizeGENOME) << "\n"; 
	  }
#else
	  if((reals[0] < sigma_11 + scale/2. - 1.5 && reals[0] > -sigma_11 + scale/2. - 1.5 && 
	      reals[1] < sigma_12  + scale/2. - 1.5 && reals[1] > -sigma_12 + scale/2.- 1.5) ||
	     (reals[0] < sigma_21 + scale/2. + 1.5 && reals[0] > -sigma_21 + scale/2. + 1.5 && 
	      reals[1] < sigma_22  + scale/2.+ 1.5 && reals[1] > -sigma_22 + scale/2. + 1.5)){
	    nrHighIndivGlobal = nrHighIndivGlobal + 1;

	    unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate

	    actualHighDistrib[indexIK] += tempNrValue;
	    actualHighFitness[indexIK] += tempNrValue/(double)tempFitness;
	  }
#endif //two picks
#endif //mixture bivariate
#endif //bivariate
#endif  //uniform distribution

#endif //histogram


#endif //MIXTURE_REALS

///iteration of the rest of the list
#ifndef MIXTURE_REALS
	  while((tempNrValue = myHashtables ->iterator(tempBuffer,sizeGENOME)) != -1){
#ifndef FLOOR
	      unsigned long tempFitness = MCMCstring[0]->fitness(tempBuffer,sizeGENOME);
#else
	      unsigned long tempFitness =  floor(MCMCstring[0]->fitness(tempBuffer,sizeGENOME));
#endif
#else
	      tempFitness =  MCMCstring[0]->fitness(tempBuffer,sizeGENOME*BLOCKsize);
#endif //MIXTURE_REALS

#ifndef MIXTURE_REALS

#ifdef THRESHOLD_DISTRIBUTION
#ifndef HALF_MAX_FITNESS
#ifdef BINOMIAL
	      if(tempFitness > (diversityMeasure - 1) * (BLOCKsize-1)/BLOCKsize){
#else
#ifdef  TRAP_LIKE_FUNCTION
		  if(tempFitness > (diversityMeasure - 1) * b/a){
#endif //trap like function
#endif //binomial
#else //half max fitness
		      if(tempFitness > (diversityMeasure - 1)/ 2){
#endif //halfmax fitness
	      //for(int i = 0; i < sizeGENOME; i++)
	      //cout << tempBuffer[i];
	      //cout << " = " << tempNrValue/(double)tempValue << "\t" << tempNrValue << "\t" << (double)tempValue 
	      //	   << "\t" << trueDistribution[tempFitness] << "\n";
	      tempDiversity[tempFitness] = tempDiversity[tempFitness] - 1;
	      
#ifndef KL_DIFFERENCE
	      if(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue) > 0){
		  difDistr += trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue);
	      } else {
		  difDistr += -(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue));
		  mari++;
	      }
#else //kl difference
	      difDistr += (trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue + 1)) * 
		(log(trueDistribution[tempFitness]) -  log((double)tempNrValue/(double)(tempValue + 1)));
	      if(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue) < 0){
		  mari++;
	      }
#endif //kl difference
	      limitaMari++;
	  }
#else //thresold distribution
	      //for(int i = 0; i < sizeGENOME; i++)
	      //cout << tempBuffer[i];
	      //cout << " = " << tempNrValue/(double)tempValue << "\t" << tempNrValue << "\t" << (double)tempValue 
	      //	   << "\t" << trueDistribution[tempFitness] << "\n";
	      tempDiversity[tempFitness] = tempDiversity[tempFitness] - 1;
	      
#ifndef KL_DIFFERENCE
	      if(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue) > 0)
		  difDistr += trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue);
	      else 
	      difDistr += -(trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue));
#else 
	      difDistr += (trueDistribution[tempFitness] - (double)tempNrValue/(double)(tempValue)) * 
		(log(trueDistribution[tempFitness]) -  log((double)tempNrValue/(double)(tempValue)));
#endif
#endif //THRESHOLD_DISTRIBUTION
	  }
	  myHashtables->iteratorFree(sizeGENOME);
#else //MIXTURE_REALS
	  //cout << "arrived before loop in histohgram\n";
	  while((tempNrValue = myHashtables ->iterator(tempBuffer,sizeGENOME*BLOCKsize)) != -1){
	    //cout << "find in hashtable \t";
	    //for(unsigned long inter = 0; inter < sizeGENOME*BLOCKsize; inter++)
	    //cout << tempBuffer[inter];
	    //cout << "= " << tempFitness << "\n";
	  convertReal(reals,tempBuffer,sizeGENOME);
	  tempFitness = MCMCstring[0]->fitness(reals,sizeGENOME);
#ifndef HISTOGRAM
	  double trueD = (double)tempFitness/(double)total;
	  double actualD = (double)tempNrValue/(double)tempValue;

#ifndef KL_DIFFERENCE
	  if(trueD > actualD)
	    difDistr += trueD - actualD;
	  else difDistr += - trueD + actualD;
#else
	  difDistr += (trueD - actualD)* (log(trueD) - log(actualD));
#endif
#else // histogram
	  //discreteze first
	      //for(int i = 0; i < sizeGENOME; i++)
	      //cout << tempBuffer[i];
	      //cout << " = " << tempNrValue/(double)tempValue << "\t" << tempNrValue << "\t" << (double)tempValue 
	      //	   << "\t" << trueDistribution[tempFitness] << "\n";

	  //convertReal(reals,tempBuffer,sizeGENOME);
	  unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate

	  actualDistrib[indexIK] += tempNrValue;
	  actualFitness[indexIK] += tempNrValue/(double)tempFitness;

#ifdef UNIFORM_DISTRIBUTION
	  //  nrHighIndivGlobal = nrHighIndivGlobal + 1;
	  //  actualHighDistrib[(unsigned long)(reals[0]*hist_points/scale)][(unsigned long)(reals[1]*hist_points/scale)]+= tempNrValue;
	  //  actualHighFitness[(unsigned long)(reals[0]*hist_points/scale)][(unsigned long)(reals[1]*hist_points/scale)]+= tempNrValue/(double)tempFitness;
#else	  
#ifdef BIVARIATE
#ifdef MIXTURE_BIVARIATE
#ifdef NINE_BIVARIATE
	  if((reals[0] < sigma_1 + scale/2. && reals[0] > -sigma_1 + scale/2. && reals[1] < sigma_2  + scale/2. && reals[1] > -sigma_2 + scale/2.) ||
	     (reals[0] < sigma_1 + scale/2. + 0.375*scale && reals[0] > -sigma_1 + scale/2. + 0.375*scale && reals[1] < sigma_2  + scale/2. + 0.375*scale && reals[1] > -sigma_2 + scale/2. + 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. + 0.375*scale && reals[0] > -sigma_1 + scale/2. + 0.375*scale && reals[1] < sigma_2  + scale/2. - 0.375*scale && reals[1] > -sigma_2 + scale/2. - 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. - 0.375*scale && reals[0] > -sigma_1 + scale/2. - 0.375*scale && reals[1] < sigma_2  + scale/2. + 0.375*scale && reals[1] > -sigma_2 + scale/2. + 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. - 0.375*scale && reals[0] > -sigma_1 + scale/2. - 0.375*scale && reals[1] < sigma_2  + scale/2. - 0.375*scale && reals[1] > -sigma_2 + scale/2. - 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. + 0.375*scale && reals[0] > -sigma_1 + scale/2. + 0.375*scale && reals[1] < sigma_2  + scale/2. && reals[1] > -sigma_2 + scale/2.) ||
	     (reals[0] < sigma_1 + scale/2. && reals[0] > -sigma_1 + scale/2. && reals[1] < sigma_2  + scale/2. + 0.375*scale && reals[1] > -sigma_2 + scale/2. + 0.375*scale) ||
	     (reals[0] < sigma_1 + scale/2. - 0.375*scale && reals[0] > -sigma_1 + scale/2. - 0.375*scale && reals[1] < sigma_2  + scale/2. && reals[1] > -sigma_2 + scale/2.) ||
	     (reals[0] < sigma_1 + scale/2. && reals[0] > -sigma_1 + scale/2. && reals[1] < sigma_2  + scale/2. - 0.375*scale && reals[1] > -sigma_2 + scale/2. - 0.375*scale)){
	    nrHighIndivGlobal = nrHighIndivGlobal + 1;

	    unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate

	    actualHighDistrib[indexIK]+= tempNrValue;
	    actualHighFitness[indexIK]+= tempNrValue/(double)tempFitness;
	    //cout << " fitness " << tempFitness << " nr Indiv " << tempNrValue << " for (" << reals[0] << "," << reals[1] << ") \n"; 
	  }
#else
	  //????????????????????????????????????????????????????????
	  //not impleemnted
	  //??????????????????????????????????????????????????????????????

	  double good = 1;
	  for(int iDimension = 0; iDimension < sizeGENOME; iDimension++)
	  for(int iMixture = 0; iMixture < mixture; iMixture++)
	    if(!(reals[iDimension] > MCMCstring[0]->positionMixture[iMixture].element(iDimension,0) - 
		 sqrt(MCMCstring[0]->covarianceMixture[iMixture].element(iDimension,iDimension)) 
		 && reals[iDimension] < MCMCstring[0]->positionMixture[iMixture].element(iDimension,0) + 
		 sqrt(MCMCstring[0]->covarianceMixture[iMixture].element(iDimension,iDimension))))
	      good = 0;
	  
	  if(good == 1){
	    nrHighIndivGlobal = nrHighIndivGlobal + 1;

	    unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate

	    actualHighDistrib[indexIK] += tempNrValue;
	    actualHighFitness[indexIK] += tempNrValue/(double)tempFitness;
	  }
#endif //9_bivariate
#else // mixture bivariate
#ifndef TWO_PICKS
	  if(reals[0] < sigma_1 * scale + scale/2. && reals[0] > -sigma_1 * scale + scale/2. && 
	     reals[1] < sigma_2 * scale + scale/2. && reals[1] > -sigma_2 * scale + scale/2.){
	    nrHighIndivGlobal = nrHighIndivGlobal + 1;
	    unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate

	    actualHighDistrib[indexIK]+= tempNrValue;
	    actualHighFitness[indexIK]+= tempNrValue/(double)tempFitness;
	    //cout << " fitness " << tempFitness << " nr Indiv " << tempNrValue << 
	    //  " for (" << reals[0] << "," << reals[1] << ") = " << MCMCstring[0]->fitness(reals,sizeGENOME) << "\n"; 
	  }
#else
	  if((reals[0] < sigma_11 + scale/2. - 1.5 && reals[0] > -sigma_11 + scale/2. - 1.5 && 
	      reals[1] < sigma_12  + scale/2.- 1.5 && reals[1] > -sigma_12 + scale/2.- 1.5) ||
	     (reals[0] < sigma_21 + scale/2. + 1.5 && reals[0] > -sigma_21 + scale/2. + 1.5 && 
	      reals[1] < sigma_22  + scale/2.+ 1.5 && reals[1] > -sigma_22 + scale/2. + 1.5)){
	    nrHighIndivGlobal = nrHighIndivGlobal + 1;

	    unsigned long indexIK = 0;
#ifndef BIVARIATE
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) reals[sizeGENOME - 1 - ik]*hist_points) * pow(hist_points,ik); 
	    }
#else //bivariate
	    for(int ik = 0; ik < sizeGENOME; ik++){
	      indexIK += ((unsigned long) (reals[sizeGENOME - 1 - ik]+scale/2)*hist_points) * pow((double)hist_points,(double)ik); 
	    }
#endif //bivariate

	    actualHighDistrib[indexIK]+= tempNrValue;
	    actualHighFitness[indexIK]+= tempNrValue/(double)tempFitness;
	  }
#endif //two picks
#endif //mixture bivariate
#endif //bivariate
#endif //uniform distribution

#endif  //histogram
	  }
	  myHashtables->iteratorFree(sizeGENOME*BLOCKsize);
#endif //MIXTURE_REALS
      }
  }

///if sometihng is not discovered add a penalty
#ifndef MIXTURE_REALS
#ifdef PENALITY
#ifndef THRESHOLD
  for(unsigned long k = 0; k < diversityMeasure; k++)
#else
#ifndef HALF_MAX_FITNESS
#ifdef BINOMIAL
  for(unsigned long k = (diversityMeasure - 1) * (BLOCKsize-1)/BLOCKsize+1; k < diversityMeasure; k++)
#else
#ifdef TRAP_LIKE_FUNCTION
  for(unsigned long k = (diversityMeasure - 1) * b/a+1; k < diversityMeasure; k++)
#endif
#endif
#else 
  for(unsigned long k = (diversityMeasure - 1) /2 + 1; k < diversityMeasure; k++)
#endif //HALF_MAX_FITNESS
#endif
      if(tempDiversity[k] > 0)
#ifndef KL_DIFFERENCE
	      difDistr += tempDiversity[k] * trueDistribution[k];
#else
	      difDistr += tempDiversity[k] * (trueDistribution[k] - 1.0/totalNr)*
		  (log(trueDistribution[k]) - log(1.0/totalNr));
#endif
      else if(tempDiversity[k] < 0)
	  cout << "Erorare in diversity " << tempDiversity[k] << " at " << k << "\n";

  cout << mari << " from " << tempValue << "\n";
  Kullback[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = mari/tempValue;
  timeKLdistance[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = limitaMari;

#else //PENALITY
  double notFound = 0, notFoundNr = 0;

#ifndef THRESHOLD
  for(unsigned long k = 0; k < diversityMeasure; k++)
#else
#ifndef HALF_MAX_FITNESS
#ifdef BINOMIAL
  for(unsigned long k = (diversityMeasure - 1) * (BLOCKsize-1)/BLOCKsize+1; k < diversityMeasure; k++)
#else
#ifdef TRAP_LIKE_FUNCTION
  for(unsigned long k = (diversityMeasure - 1) * b/a+1; k < diversityMeasure; k++)
#endif
#endif
#else 
  for(unsigned long k = (diversityMeasure - 1) /2 + 1; k < diversityMeasure; k++)
#endif //HALF_MAX_FITNESS
#endif
      if(tempDiversity[k] > 0){
	  notFoundNr += tempDiversity[k];
#ifdef KL_DIFFERENCE
	  notFound += tempDiversity[k] * trueDistribution[k] * log(trueDistribution[k]);
#endif
      }
      else if(tempDiversity[k] < 0)
	  cout << "Erorare in diversity " << tempDiversity[k] << " at " << k << "\n";

  cout << mari << " from " << tempValue << "\n";
  Kullback[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = notFound;
  timeKLdistance[runs-1][(unsigned long)(generation/SEE_SAMPLE)] = (double)limitaMari/(double)(notFoundNr+limitaMari);

#endif //PENALITY
#endif //ifndef MIXTURE_REALS

      /*cout << "enum indivizi \n f(";
  for(unsigned long t = 0; t < sizeGENOME; t++)
      cout << tempBuffer[t];
      cout << ") = " << tempNrValue << " difDistr =" << difDistr <<" \n";*/

  //cout << "verifica KL procedure  difDistr="<< difDistr << "fitness " <<tempFitnessValue<< "\n";

#endif //FILE_BOOKIPEENG 

#ifdef HISTOGRAM
  difDistr = 0;
  for(unsigned long i = 0; i < diversityMeasure; i++){
      double dif = 0;
      if(actualDistrib[i] != 0 && actualFitness[i] != 0) 
	dif = trueDistribution[i] - actualDistrib[i]/actualFitness[i];
      else dif = trueDistribution[i];
      if(dif < 0) dif = -dif;
      difDistr += dif;
    }

  //print the information about the best point from the distribution 
#ifdef BIVARIATE
  difDistr = 0;
  cout << " cum se calculeaza distribution \n";
  double diff = 0, true_d = 0;
  int count_true = 0;
  for(unsigned long i = 0; i < diversityMeasure; i++){
      double dif = 0;
      if(actualHighDistrib[i] != 0 && actualHighFitness[i] != 0){ 
	//	cout << i << "\t" << j << "\t" << actualHighDistrib[i][j] << "\t" << actualHighFitness[i][j] << "\t" << 
	//  actualHighDistrib[i][j]/actualHighFitness[i][j] << "\t" << trueDistribution[i][j] << "\n"; 
	difDistr += actualHighDistrib[i]/actualHighFitness[i];
	diff += actualHighDistrib[i];
	true_d += trueDistribution[i];
	count_true++;
      }
    }
  if(nrHighIndivGlobal != 0){
    cout << " mean most " << diff << "   " << difDistr << " true mean " << true_d << " " << count_true << "\n";
    difDistr = (difDistr/count_true)/(double)(true_d/count_true);
  }
  else difDistr = 0;

#endif //BIVARIATE

#endif  
  return difDistr;
}

double binaryPopulation::printMixingGNUFILE(){
  ofstream myMixing(MixingFILE);
  if(!myMixing.is_open()){
    cout <<"Error in opening the mixing file "<< MixingFILE <<"\n";
    return -1;
  }

  for(unsigned long i = 0 ; i < sampleSize / SEE_SAMPLE; i++)
    myMixing << (Kullback[i][0]/runs) << "\n";

  myMixing.flush();
  myMixing.close();
  return 0;

}
           
#ifndef TWO_ATTRACTORS_FUNCTION
void binaryPopulation::readHistogram(){
//1.read the true distribution from the file
  ifstream distributionFile;
  distributionFile.open(distributionFILE,ios::in);
  if(!distributionFile.is_open()){
      cout << "Error in opening the distribution File" <<distributionFILE<< "\n";
  }

  ifstream mySampleFile(sampleFile);		
  if(!mySampleFile.is_open())
      cout << "Error opening file "<< sampleFile << "convergence.cpp:793\n";


  distributionFile >> total;
  distributionFile >> mean;

  unsigned long k = 0;

#ifndef MIXTURE_REALS  
  //cout << "Write true distribution = ";
  while(!distributionFile.eof() && k < diversityMeasure){
	    if(k >= diversityMeasure) {
		cout << "Error in diversity measure convergence.cpp:724\n";
		exit(1);
	    } else{
		distributionFile >> trueDistribution[k];
		//	cout << trueDistribution[k] << "\t";
		k++;
	    }
  }
#else //mixture reals
#ifdef HISTOGRAM
  while(!distributionFile.eof() && k < diversityMeasure){
	    if(k >= diversityMeasure) {
		cout << "Error in diversity measure convergence.cpp:891\n";
		exit(1);
	    } else{
	      distributionFile >> trueDistribution[k];
	      //trueDistribution[k][j] = diversityTrueDistribution[k][j];
	      //	cout << trueDistribution[k] << "\t";
	      k++;
	    }
  }
#endif //histogram
#endif //MIXTURE_REALS
  //cout << "\n";
  distributionFile.close();

//2. read a sampled file
#ifndef MIXTURE_REALS
  k = 0;
  //cout << "Diversity true distribution \n";
  while(!mySampleFile.eof() && k < diversityMeasure){
	    if(k >= diversityMeasure) {
		cout << "Error in diversity measure convergence.cpp:724\n";
		exit(1);
	    } else{
		mySampleFile >> diversityTrueDistribution[k];
		k++;
		//cout << i << "\t" <<diversityTrueDistribution[i]<<"\t" << tempDiversity[k] << "\n";
	    }
  }
#else //mixture reals
#ifdef HISTOGRAM
  while(!mySampleFile.eof() && k < diversityMeasure){
	    if(k >= diversityMeasure) {
		cout << "Error in diversity measure convergence.cpp:724\n";
		exit(1);
	    } else{
	      mySampleFile >> diversityTrueDistribution[k];
	      k++;
	      //cout << i << "\t" <<diversityTrueDistribution[i]<<"\t" << tempDiversity[k] << "\n";
	    }
  }
#endif //histogram
#endif //MIXTURE_REALS
  mySampleFile.close();
 }
#endif //two attractors distribution
#endif //KULLBACK_INFORMATION
#endif //histogram

#ifdef CUSUM
//metoda grafica de masurare a convergentei CUSUM, 
//univariata
double binaryPopulation::CusumCalculation(){
  double medie = 0;
  double trecut = 0, acum = 0, viitor = 0;
  double grafic = 0;
  unsigned long burn = 0, set_burn = 0;
//#ifdef BURN_IN
  double BurnT[constBurn];
  for(unsigned long i = 0; i < constBurn; i++)
    BurnT[i] = 0;
//#endif

  if(sampleSize < 3) return 0;

#ifndef MEAN_UNKNOWN
  cout << "We have the mean " << mean <<"\n";
#endif

#ifdef BURN_IN 
  //Calculeaza burn-in period
  //calcularea secventelor de Hairness
  for(unsigned long k = 1; k <= constBurn; k++){
    //calcularea burn-in
    unsigned long n = ((sampleSize+ 1)/constBurn) * k;
    burn = n / 2;
    grafic = 0; medie = 0; trecut = 0; acum = 0; viitor = 0;

    //calcularea mediei
    //cout << "k=" << k << "\t n=" << n << "\t burn=" << burn << "\n";

    for(unsigned long j = 0; j < sizeMCMCstring; j++) {
      //cout << "\n Pentru ("<< j << ","<<n << "," << burn << ")\t";     
      medie += MCMCstring[j] -> getElementValue(n,burn + 1);
    }
    medie /= (n - burn_in) * sizeMCMCstring;
    //cout << "medie " << medie << " over "<< n* sizeMCMCstring/2<< "\n";
    
    //primul trecut si acum
    for(unsigned long i = 0; i < sizeMCMCstring; i++){
      trecut += MCMCstring[i] -> getElementValue(burn) - medie;
      acum += MCMCstring[i] -> getElementValue(burn + 1) - medie;
    }

    //calcularea grafic
    double acelasi = -1;

    for(unsigned long i = burn + 2; i <= n; i++){
      for(unsigned long j = 0; j < sizeMCMCstring; j++)
	viitor += MCMCstring[j] -> getElementValue(i) - medie;
    
#ifdef STICKY
      if ((trecut > acum && acum < viitor) || 
	  (trecut < acum && acum > viitor) || 
	  (trecut = acum && acelasi != -1 && ((acelasi > acum && acum < viitor) || (acelasi <acum && acum > viitor)))) 
	{grafic += 1;acelasi = -1;}
      else if(trecut == acum && viitor == acum) {grafic += 0.5;}
      
      if(viitor == acum) acelasi = trecut;
      else acelasi = -1;
#else
      if ((trecut > acum && acum < viitor) || (trecut < acum && acum > viitor))
	grafic += 1;
#endif

      trecut = acum;
      acum = viitor;
      viitor = 0;
    }
    BurnT[k-1] = (double)grafic/(double)burn;

    //cout << "BurnT[" << k-1 << "]=" << BurnT[k-1] << "\n";
  }

  //calcularea burn-in prin stabilizare
  set_burn = constBurn;

#ifdef BURN_IN_STABLE
#ifndef BURN_IN_NORMAL

  for(unsigned long i = 0; i < constBurn - 2; i++)
    for(unsigned long j = i+1; j < constBurn-1; j++)
        if(BurnT[i] - BurnT[j] > epsilon || BurnT[i] - BurnT[j] < -epsilon) {
	  set_burn = i + 1;
	  //cout << "BurnT[" << i << "]=" << BurnT[i] << "- BurnT[" << j << "]=" << BurnT[j] << " is " 
	  //   << abs(BurnT[i] - BurnT[j]) << " which is >" << epsilon << "\n"; 
	}
#else
 for(unsigned long i = 0; i < constBurn-1; i++)
   if(BurnT[i] < 0.5 + epsilon && BurnT[i] > 0.5 - epsilon) continue;
   else      
     set_burn = i+1;
#endif //burn_in normal
#else //burn_in_stable
  //se presupune ca distributia este simetrica fata de medie
  for(unsigned long i = 0; i < constBurn-1; i++){
    double simetricalC1 = 0.5 + CI*sqrt(1.0/(4*( (sampleSize + 1)/constBurn *(i+1)*sizeMCMCstring - (sampleSize + 1)/constBurn *(i+1)*sizeMCMCstring/2 -2) ));
    double simetricalC2 = 0.5 - CI*sqrt(1.0/(4*( (sampleSize + 1)/constBurn *(i+1)*sizeMCMCstring - (sampleSize + 1)/constBurn *(i+1)*sizeMCMCstring /2 -2) ));
    
    //cout << "BurnT[" << i << "]=" << BurnT[i] << " with up limit=" << simetricalC1 << " and down limit=" << simetricalC2 << "\n";
    
    if(BurnT[i] < simetricalC1 && BurnT[i] > simetricalC2) continue;
    else { 
      set_burn = i+1;
      //cout << "the burn in is set further " << i << "\n";
    }
  }
#endif //BURN_IN_STABLE

  if(set_burn < constBurn - 1){
    burn_in = (unsigned long) ((sampleSize + 1)/constBurn)*set_burn/2;
    cout << "Sa gasit burn_in" << burn_in << "\n";
  } else {
    cout << "Nu este burn in in aceastta secventa alegere default burn_in\n";
    burn_in = (sampleSize+1)/2;
  }

  burn_inT[runs-1] = burn_in;
#else //no burn in
  burn_in = 0;
  //for(unsigned long i = 0; i < constBurn; i++)
  //    BurnT[i] = sampleSize/2;
  //burn_in = sampleSize/2;
#endif //BURN_IN

  //Calculeaza partial sums and hairness

    grafic = 0;  trecut = 0; acum = 0; viitor = 0;
#ifdef MEAN_UNKNOWN
//#ifndef KULLBACK_INFORMATION
   //calcularea mediei
    medie = 0;
    for(unsigned long j = 0; j < sizeMCMCstring; j++) {
       double sumChain = MCMCstring[j] -> getElementValue(sampleSize,burn_in + 1);
       cout << "chain " << j << "has mean " << sumChain/sampleSize << "\n";
       medie += sumChain;
    }
    medie /= (sampleSize - burn_in)* sizeMCMCstring;
    cout << "mean unknown " << medie << " bucati "<< sampleSize* sizeMCMCstring << "\n";
#else 
    medie = mean;
    cout << "mean known" << medie << " bucati "
	 << (sampleSize-burn_in)* sizeMCMCstring << "\n";
#endif //MEAN_UNKNOWN    

#ifndef SINGLE_vs_MULTIPLE
    //primul trecut si acum
    for(unsigned long i = 0; i < sizeMCMCstring; i++){
      trecut += MCMCstring[i] -> getElementValue(burn_in + 1) - medie;
      acum += MCMCstring[i] -> getElementValue(burn_in + 2) - medie;
    }

    double convergence = 0;

    double acelasi = -1;

    for(unsigned long i = burn_in+2; i < sampleSize; i++){
      for(unsigned long j = 0; j < sizeMCMCstring; j++)
	viitor += MCMCstring[j] -> getElementValue(i+1) - medie;
      
#ifdef STICKY
      if ((trecut > acum && acum < viitor) || 
	  (trecut < acum && acum > viitor) || 
	  (trecut = acum && acelasi != -1 && ((acelasi > acum && acum < viitor) || (acelasi <acum && acum > viitor)))) 
	{grafic += 1;acelasi = -1;}
      else if(trecut == acum && viitor == acum) {grafic += 0.5;}
      
      if(viitor == acum) acelasi = trecut;
      else acelasi = -1;
#else
      if ((trecut > acum && acum < viitor) || (trecut < acum && acum > viitor))
	grafic += 1;
#endif
      
#ifndef multiple_restarts
      Hairness[i-2] = grafic/(double)(i- (burn_in + 1));
      //cout << "Cusum[" << runs-1 << "][" << i-2 << "]=" <<Cusum[runs-1][i-2] << "\n";
      Cusum[i-2] = acum;
#else //nr_runs 1
      Hairness[runs-1][i-2] = grafic/(double)(i- (burn_in + 1));
      //cout << "Cusum[" << runs-1 << "][" << i-2 << "]=" <<Cusum[runs-1][i-2] << "\n";
      Cusum[runs-1][i-2] += acum;
#endif //nr_runs
      
      //cout << "(" << (i - 2) << ", " << grafic << "," << Cusum[runs-1][i - 2] << ","<< Hairness[runs-1][i - 2] <<")\n"; 
#ifdef CONVERGENCE_CUSUM
#ifdef NORMAL_FUNCTION
#ifndef BURN_IN
      double simetricalC1 = mean + CI*sqrt(mean*(1-mean))/sqrt(sampleSize*(i-1));
      double simetricalC2 = mean - CI*sqrt(mean*(1-mean))/sqrt(sampleSize*(i-1));
   
#ifndef multiple_restarts
      if((Hairness[i-2] < simetricalC1 && Hairness[i-2] > simetricalC2)) 
#else
      if((Hairness[runs-1][i-2] < simetricalC1 && Hairness[runs-1][i-2] > simetricalC2)) 
#endif 
	{
	  // cout << "Hainness[" << (runs-1) << "][" << (i-2) << "]=" << Hairness[runs-1][i-2]
	  //   << " where \mu+e = " << simetricalC1 << " \mu-e=" << simetricalC2 << "\n";
	  continue;
      }
      else {
	convergence = i;
	//cout << "Convergence is not present at step " << i -1 << "\n";
      }
#endif
#else
      double meanT = BurnT[set_burn];
      double simetricalC1 = meanT + CI*sqrt(meanT*(1.0-meanT)/( (sampleSize + 1)/constBurn *(i+1) - (sampleSize + 1)/constBurn *(i+1)*0.5 -1));
      double simetricalC2 = meanT - CI*sqrt(meanT*(1.0-meanT)/((sampleSize + 1)/constBurn *(i+1) - (sampleSize + 1)/constBurn *(i+1)*0.5 -1));
   
#ifndef multiple_restarts
      if((Hairness[i-2] < simetricalC1 && Hairness[i-2] > simetricalC2)) 
#else
      if((Hairness[runs-1][i-2] < simetricalC1 && Hairness[runs-1][i-2] > simetricalC2)) 
#endif //nr_runs
	  continue;
     else {
	convergence = i - 1;
	//cout << "Convergence is not present at step " << i -1 << "\n";
      }
      
#endif //NORMAL_FUNCTION
#endif //CONVERGENCE CUSUM

      trecut = acum;
      acum = viitor;
      viitor = 0;
    
#else //single_vs_multiple
    
    //primul trecut si acum
    for(unsigned long j = burn_in+2; j < sampleSize; j++){
	for(unsigned long i = 0; i < sizeMCMCstring; i++){
	    trecut += MCMCstring[i] -> getElementValue(burn_in + 1) - medie;
	    acum += MCMCstring[i] -> getElementValue(burn_in + 2) - medie;
	    
	    double convergence = 0;
	    
	    double acelasi = -1;
	    
	    viitor += MCMCstring[i] -> getElementValue(j+1) - medie;
	    
#ifdef STICKY
	    if ((trecut > acum && acum < viitor) || 
		(trecut < acum && acum > viitor) || 
		(trecut = acum && acelasi != -1 && ((acelasi > acum && acum < viitor) || (acelasi <acum && acum > viitor)))) 
	    {grafic += 1;acelasi = -1;}
	    else if(trecut == acum && viitor == acum) {grafic += 0.5;}
	    
	    if(viitor == acum) acelasi = trecut;
	    else acelasi = -1;
#else
	    if ((trecut > acum && acum < viitor) || (trecut < acum && acum > viitor))
		grafic += 1;
#endif
	    
	    Hairness[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][(j-2)*sizeMCMCstring+i] = 
		grafic/(double)(j - (burn_in + 1));
	    //cout << "Cusum[" << runs-1 << "][" << i-2 << "]=" <<
	    //Cusum[runs-1][i-2] << "\n";
	    Cusum[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][(j-2)*sizeMCMCstring+i] 
		+= acum;
	    
	    //cout << "(" << (i - 2) << ", " << grafic << "," << Cusum[runs-1][i - 2] << ","<< Hairness[runs-1][i - 2] <<")\n"; 
#ifdef NORMAL_FUNCTION
#ifndef BURN_IN
	    double simetricalC1 = mean + CI*sqrt(mean*(1-mean))/sqrt(sampleSize*(i-1));
	    double simetricalC2 = mean - CI*sqrt(mean*(1-mean))/sqrt(sampleSize*(i-1));
	    
	    if((Hairness[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][(j-2)*sizeMCMCstring+i]< simetricalC1&& 
		Hairness[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][(j-2)*sizeMCMCstring+i] > simetricalC2)){
		// cout << "Hainness[" << (runs-1) << "][" << (i-2) << "]=" 
		//<< Hairness[runs-1][i-2] << " where \mu+e = " << simetricalC1 
		//<< " \mu-e=" << simetricalC2 << "\n";
		continue;
	    }
	    else {
		convergence = i;
		//cout << "Convergence is not present at step " << i -1 << "\n";
	    }
#endif
#else
	    double meanT = BurnT[set_burn];
	    double simetricalC1 = meanT + CI*sqrt(meanT*(1.0-meanT)/
		       ( (sampleSize + 1)/constBurn *((j+1)*sizeMCMCstring+i) - 
		       (sampleSize + 1)/constBurn *((j+1)*sizeMCMCstring+i)*0.5 -1));
	    double simetricalC2 = meanT - CI*sqrt(meanT*(1.0-meanT)/
		       ((sampleSize + 1)/constBurn *((j+1)*sizeMCMCstring+i) - 
		       (sampleSize + 1)/constBurn *((j+1)*sizeMCMCstring+i)*0.5 -1));
	    
	    if((Hairness[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][(j-2)*sizeMCMCstring+i] 
		< simetricalC1 && 
		Hairness[(unsigned long)(log(sizeMCMCstring)/log(2.0))][runs-1][(j-2)*sizeMCMCstring+i] 
		> simetricalC2)) 
		continue;
	    else 
		convergence = i - 1;
	    //cout << "Convergence is not present at step " << i -1 << "\n";
	    
#endif
	    
	    trecut = acum;
	    acum = viitor;
	    viitor = 0;
#endif
    
	    if(convergence < sampleSize-1){
		//cout << "Convergenta sa produs la pasul" << convergence << " with burn-in " 
		//    << burn_in << "\n";
	    } else {
		//cout << "Convergenta not achived \n";
	    }
	}
#ifdef SINGLE_vs_MULTIPLE
    }
    }
#endif
    return 0;
}

double binaryPopulation::printCusumGNUFILE(){
    ofstream myCUSUM(CusumFILE);
    if(!myCUSUM.is_open()){
	cout <<"Error in opening the mixing file "<< CusumFILE <<"\n";
	return -1;
    }
    
#ifdef BURN_IN
    unsigned long maxim = burn_inT[0];
    for(unsigned long j = 1; j < nr_runs; j++)
//#ifndef SINGLE_vs_MULTIPLE

	if(burn_inT[j] > maxim) maxim = burn_inT[j];
    myCUSUM << maxim << "\n";
    
  for(unsigned long i = maxim; i < sampleSize-1; i++){
      for(unsigned long j = 0; j < nr_runs; j++)
	  myCUSUM << Cusum[j][i]/(double)(i+1) << "\t" << Hairness[j][i] <<"\t";
      myCUSUM << "\n";
  }
#else
#ifndef SINGLE_vs_MULTIPLE
//printeaza o suma a cusum si st dev de la burn_in
  double cusum_mean, hairness_mean;
  double std_dev_c, std_dev_h;
  //cout << " sampleSize " << sampleSize << " sampleSize - 2" << sampleSize - 2 << " burn_in " << burn_in << "\n"; 
  if(sampleSize - 2 >= burn_in){
    for(int i = burn_in; i < sampleSize-2; i++){
#ifndef multiple_restarts
      myCUSUM << i << "\t" << Cusum[i]/(double)(i+1) << "\t" << Hairness[i] << "\n";
#else
      cusum_mean = 0; hairness_mean = 0;
      std_dev_c = 0; std_dev_h = 0;
      for(unsigned long j = 0; j < nr_runs; j++){
        //myCUSUM << "\n";
	cusum_mean += Cusum[j][i];
	hairness_mean += Hairness[j][i];
      }
      cusum_mean /= nr_runs;
      hairness_mean /= nr_runs;
      for(unsigned long j = 0; j < nr_runs; j++){
	//myCUSUM << Cusum[j][i]/(double)(i+1) << "\t" << Hairness[j][i] <<"\t";
	//myCUSUM << "\n";
	std_dev_c += pow(cusum_mean - Cusum[j][i],2);
	std_dev_h += pow(hairness_mean - Hairness[j][i],2);
      }      
      myCUSUM << i << "\t" << Cusum[0][i]/(double)(i+1) <<
	"\t" << cusum_mean << "\t" << pow(std_dev_c/nr_runs,0.5) 
	      << "\t" << hairness_mean << "\t" << pow(std_dev_h/nr_runs,0.5);
      //myCUSUM << "\t";
      //for(unsigned long j = 0; j < nr_runs; j++){
      //  myCUSUM << Cusum[j][i] << "\t" << Hairness[j][i] <<"\t";
      // }
      myCUSUM << "\n";
    
#endif //nr_runs
    }
  }

#else
  //int tempSize = (int)(log(sizeMCMCstring)/log(2.0));
  for(unsigned long j = 0; j < (unsigned long)sampleSize*sizeMCMCstring; j++){
      for(unsigned long k = 0; k < nr_runs; k++)
	  for(unsigned long i = 0; i < (unsigned long)(log(sizeMCMCstring)/log(2.0)) + 1; i++)
	      myCUSUM << Cusum[i][k][j]/(double)(j+1)<<"\t"<<Hairness[i][k][j]<<"\t";
      myCUSUM << "\n";
      // if(j > (int)sampleSize*pow(tempSize-i,2))
      //	  myCUSUM << "-" << "\t" << "-" << "\t";
      //else
  }
  
#endif // ifndef single vs multiple
#endif //BURN_IN
  
  myCUSUM.flush();
  myCUSUM.close();
  
  return std_dev_h;
}

#endif //CUSUM
#ifndef EXPAND_KULLBACK_INFORMATION
#ifdef TWO_ATTRACTORS_FUNCTION
#ifndef MIXTURE_REALS
void binaryPopulation::readHistogram(){
  ifstream histogramFile;
  histogramFile.open(histogramFILE,ios::in);
  if(!histogramFile.is_open()){
      cout << "Error in opening the distribution File" <<distributionFILE<< "\n";
  }

  unsigned long k = 0;
  unsigned long fInt = 0;
//cout << "Write true distribution = ";
  while(!histogramFile.eof() && k <= sizeGENOME){
	    if(k > sizeGENOME) {
		cout << "Error in diversity measure convergence.cpp:670\n";
		//exit(1);
	    } else{
		histogramFile >> fInt;
		histogramFile >> histogram[k];
		//	cout << trueDistribution[k] << "\t";
		k++;
	    }
  }
}
#else 
#ifdef BINOMIAL
/*void binaryPopulation::readHistogram(){
  ifstream histogramFile;
  histogramFile.open(histogramFILE,ios::in);
  if(!histogramFile.is_open()){
      cout << "Error in opening the distribution File" <<distributionFILE<< "\n";
  }

  unsigned long k = 0;
  int fInt = 0;
//cout << "Write true distribution = ";
  while(!histogramFile.eof() && k <= sizeGENOME){
	    if(k > sizeGENOME) {
		cout << "Error in diversity measure convergence.cpp:670\n";
		//exit(1);
	    } else{
		histogramFile >> fInt;
		histogramFile >> histogram[k];
		//	cout << trueDistribution[k] << "\t";
		k++;
	    }
  }
  }*/
#endif
#endif
#endif
#endif

#ifdef DISTANCE
#ifndef MIXTURE_REALS
#ifndef RECOMB_MUT
 double binaryPopulation::HammingDistance(unsigned long i){
     double temp = 0;
     unsigned long tempExpGenome[sizeGENOME];
     MCMCstring[i]->getGENOME(tempExpGenome);
//#ifdef BINOMIAL
     unsigned long opGENOME[sizeGENOME];
     MCMCstring[0]->getMax0(opGENOME);
     for(unsigned long i = 0; i < sizeGENOME; i++){
	 if(opGENOME[i] != tempExpGenome[i]){
	    temp++;
	    //cout << opGenome[i] << "||";
	 }
    }
//#endif
     //temp = temp/(double)sizeMCMCstring;
     //cout << " = " << temp << "\n";
    return temp;
}

double binaryPopulation::HammingDistance(unsigned long* opGENOME){
     double temp = 0;
    for(unsigned long j = 0; j < sizeMCMCstring; j++){
	unsigned long tempExpGenome[sizeGENOME];
	MCMCstring[j]->getGENOME(tempExpGenome);
	for(unsigned long i = 0; i < sizeGENOME; i++)
	    if(opGENOME[i] != tempExpGenome[i])
		temp++;
    }
    temp /= sizeMCMCstring;
    return temp;
} 
#else
 double binaryPopulation::HammingDistance(unsigned long* opGENOME){
     double temp = 0;
     for(unsigned long j = 0; j < (sizeGENOME/(double)BLOCKsize); j++){
#ifdef BINOMIAL 
	 for(unsigned long i = 0; i<BLOCKsize/2; i++)
	     if(opGENOME[i + j*BLOCKsize] != 0){
		 temp++;
//		cout << opGENOME[i + j*BLOCKsize];
	     }
	 for(unsigned long i = BLOCKsize/2; i<BLOCKsize; i++)
	     if(opGENOME[i + j*BLOCKsize] != 1){
		 temp++;
		 //	cout << opGENOME[i + j*BLOCKsize];
	     }
#else
	 for(unsigned long i = 0; i<BLOCKsize; i++)
	     if(opGENOME[i + j*BLOCKsize] != 0){
		 temp++;
	     }
#endif
     }
     // cout << " dist " << temp;
     return temp;
 } 
 
double binaryPopulation::HammingDistance(unsigned long i){
     double temp = 0;
     unsigned long tempExpGenome[sizeGENOME];
     MCMCstring[i]->getGENOME(tempExpGenome);
//#ifdef BINOMIAL
     unsigned long opGENOME[sizeGENOME];
     MCMCstring[0]->getMax0(opGENOME);
     for(unsigned long i = 0; i < sizeGENOME; i++){
	 if(opGENOME[i] != tempExpGenome[i])
	     temp++;
     }
     //temp /= sizeMCMCstring;
//#endif
     return temp;
}
#endif //recomb mut
#else //mixture reals
/////////////////////////////
////// it is actualy the euclidian distance
////////////////////////////
double binaryPopulation::HammingDistance(unsigned long i){
     double temp = 0;
     double tempExpGenome[sizeGENOME];
     MCMCstring[i]->getGENOME(tempExpGenome);
     double opGENOME[sizeGENOME];
     MCMCstring[0]->getMax0(opGENOME);
     for(unsigned long i = 0; i < sizeGENOME; i++)
	 temp += pow(opGENOME[i] - tempExpGenome[i],2);
     temp = sqrt(temp); 
//#endif
     //temp = temp/(double)sizeMCMCstring;
     //cout << " = " << temp << "\n";
    return temp;
}

double binaryPopulation::HammingDistance(double* opGENOME){
    double tempTotal = 0;
    double tempExpGenome[sizeGENOME];
    
    for(unsigned long j = 0; j < sizeMCMCstring; j++){
	double temp =0;
	MCMCstring[j]->getGENOME(tempExpGenome);
	for(unsigned long i = 0; i < sizeGENOME; i++)
	    temp += pow(opGENOME[i] - tempExpGenome[i],2);
	tempTotal += sqrt(temp); 
    }
    tempTotal /= sizeMCMCstring;
    return tempTotal;
} 
 
double binaryPopulation::HammingDistance(double* opGENOME, double* secondGENOME){
     double temp = 0;

     for(unsigned long i = 0; i < sizeGENOME; i++){
	 temp += pow(opGENOME[i]- secondGENOME[i],2);
     }
     return sqrt(temp);
 } 

#ifdef HISTOGRAM
double binaryPopulation::averageDistanceGoodIndiv(){
  double tempValue =  myHashtables -> nrElem;
  double total = tempValue;
  double difDistr = 0;
  double realBuffer[sizeGENOME+1], realsecondBuffer[sizeGENOME+1];
  unsigned long tempBuffer[sizeGENOME*BLOCKsize+1], secondBuffer[sizeGENOME*BLOCKsize+1];

  // cout << "calcularea distantei \n";
  double tempNrValue = myHashtables -> iteratorInit((unsigned long*)tempBuffer,sizeGENOME*BLOCKsize);

  if(tempNrValue != -1){
	//calculeaza distanta cu ceilalti indivizi din hashtable
    double tempNrSecondValue =  myHashtables -> iteratorSecondInit((unsigned long*)tempBuffer,sizeGENOME*BLOCKsize);
    //unsigned long ate = (unsigned long) HammingDistance(tempBuffer,secondBuffer);
    double localDist = 0; //ate;
    //localDistr += ate;
    total--;
    /*for(unsigned long i = 0; i < sizeGENOME; i++)
      cout << tempBuffer[i];
      cout << " = " << tempNrValue/(double)tempValue << "\t" << trueDistribution[ate] << "\n";*/
    
    while((tempNrSecondValue = myHashtables ->iteratorSecond((unsigned long*)secondBuffer,sizeGENOME*BLOCKsize)) != -1){
	convertReal(realBuffer,tempBuffer,sizeGENOME);
	convertReal(realsecondBuffer,secondBuffer,sizeGENOME);
	double ate = (double) HammingDistance(realBuffer,realsecondBuffer);
	localDist += ate;
    }    
    if(total != 0)
      difDistr += localDist/total;

    while((tempNrValue = myHashtables ->iterator((unsigned long*)tempBuffer,sizeGENOME*BLOCKsize)) != -1){

      double tempNrSecondValue =  myHashtables -> iteratorSecondInit((unsigned long*)tempBuffer,sizeGENOME*BLOCKsize);
      //unsigned long ate = (unsigned long) HammingDistance(tempBuffer,secondBuffer);
      double localDist = 0; //ate;
      //localDistr += ate;
      total--;
      //	      unsigned long ate = (unsigned long) HammingDistance(tempBuffer);
      //	      double localDist = ate;
      //	      localDistr += ate;
      
      /*for(unsigned long i = 0; i < sizeGENOME; i++)
	cout << tempBuffer[i];
	cout << " = " << tempNrValue/(double)tempValue << "\t" << trueDistribution[ate] << "\n";*/
      
      while((tempNrSecondValue = myHashtables ->iteratorSecond((unsigned long*)secondBuffer,sizeGENOME*BLOCKsize)) != -1){
	convertReal(realBuffer,tempBuffer,sizeGENOME);
	convertReal(realsecondBuffer,secondBuffer,sizeGENOME);
	double ate = (double) HammingDistance(realBuffer,realsecondBuffer);
	localDist += ate;
      }    

      if(total != 0)
	difDistr += localDist/total;
      
    }
    myHashtables->iteratorFree(sizeGENOME);
  }
  
  return difDistr;
}
#endif //histogram 
#endif //MIXTURE_REALS

#endif //DISTANCE
 
#ifdef QBP 
void binaryPopulation::generateQBP(unsigned long sizeGENOME, double** matrixQBP){
    myQBP << sizeGENOME << "\n";
    for(unsigned long i = 0; i < sizeGENOME; i++)
	for(unsigned long j = 0; j < sizeGENOME; j++){
	    matrixQBP[j][i] = genrand_real1();
	    myQBP << matrixQBP[j][i] << "\n";
	}
}
void binaryPopulation::readQBP(unsigned long sizeGENOME, double** matrixQBP){
    //myQBP << sizeGENOME << "\n";
    for(unsigned long i = 0; i < sizeGENOME; i++)
	for(unsigned long j = 0; j < sizeGENOME; j++){
	    myQBP >> matrixQBP[j][i];
	}
}
#endif //QBP
 
