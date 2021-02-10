//write the histogram of the

#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
                                     // newmatap.h will get include.h
#include<iostream.h>
#include <fstream.h>
#include <stdlib.h>
//#include <math.h>

#include "MatrixUtil\newmatap.h"                // need matrix applications
#include "MatrixUtil\newmatio.h"                // need matrix output routines

#include"random.cpp"
//#include "newmat.h"                // need matrix applications

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

//maximum fitness
double maxFitness(Matrix* covarianceMatrix, int dimensionMy, int mixtureMy){
  double value = 0; 
  Matrix multiM(1,1);

  //compute the values
  for(int i = 0; i < mixtureMy; i++){
    //cout << " mixiture " << i << "\n";

    double tempValue = 0;

    tempValue = pow(2.0 * pi,-((double)dimensionMy)/2.0) * pow(abs(covarianceMatrix[i].Determinant()),-1.0/2.0);

    value += tempValue;

    //cout << " exp " << exp(-0.5 * multiM.element(0,0)) << "value " << tempValue << "\n\n";
  }

  return value;
}

//compute the fitness of the covariance matrix

double fitness(Matrix x_current, Matrix* covarianceMatrix, Matrix* covarianceMatrixInverse, Matrix* positionMatrix, int dimensionMy, int mixtureMy){
  double value = 0; 
  Matrix multiM(1,1);

  //compute the values
  for(int i = 0; i < mixtureMy; i++){
    //cout << " mixiture " << i << "\n";

    double tempValue = 0;

    multiM = (x_current - positionMatrix[i]).t() * covarianceMatrixInverse[i] * (x_current - positionMatrix[i]);

    //cout << " position mixture " << positionMatrix[i].t() << " difference with current " << (x_current - positionMatrix[i]).t() 
    //	 << "\n half matrix " << (x_current - positionMatrix[i]).t() * covarianceMatrixInverse[i]
    // << " multiple with " << multiM << " determinant " << abs(covarianceMatrix[i].Determinant()) << "\n"; 

    tempValue = pow(2.0 * pi,-((double)dimensionMy)/2.0) * pow(abs(covarianceMatrix[i].Determinant()),-1.0/2.0) * exp(-0.5 * multiM.element(0,0));
    
    value += tempValue;

    //cout << " exp " << exp(-0.5 * multiM.element(0,0)) << "value " << tempValue << "\n\n";
  }
  value *= multiplication;
#ifdef AwayFromZero
  if(value < AwayFromZero)  value = AwayFromZero;
#endif
  return value;
}
 
//global variables
int NrBins_h=120; // nr of individuals in a bin
int NrDivBin=10; // nr of splits in a bin


double histogram(Matrix* covarianceMatrix, Matrix* positionMatrix, int dimension, int lenght, int mixture){
  Matrix x_current(dimension,1);
  cout << "multiply the covariance with scale - lenght \n ";
  for(int i = 0; i < mixture; i++){
    //for(int j = 0; j < dimension; j++)
    //  for(int k = 0; k < dimension; k++)
    //covarianceMatrix[i].element(j,k) *= lenght/2;
    cout << " covariance\n" << covarianceMatrix[i]; 
  }
  Matrix* covarianceInverse = new Matrix[mixture];
  for(int i = 0; i < mixture; i++){
    covarianceInverse[i].ReSize(dimension, dimension);
  }
  for(int i = 0; i < mixture; i++){
    covarianceInverse[i] = covarianceMatrix[i].i();
  }

  double kl_distance = 0;

  cout << " trect de asta \n";
  double total = 0;
  int point[dimension];
  for(int i = 0; i < dimension; i++)
    point[i] = 0;

  //statistical tables
  ofstream myFile("distribution.file");		
  if(!myFile.is_open()){
    cout << "Error opening file \n";
    exit(1);
  }
  ofstream histFile("histogram.file");		
  if(!histFile.is_open()){
    cout << "Error opening file \n";
    exit(1);
  }
  double maxDiff= 0;
  double cov12 = 0;
  double mean1 = 0, mean2 = 0;
  double sqrTotal = 0, sqrBin = 0;
  unsigned long totalHigh = 0;
  //for all possible points that you can think of
  //for each of this points compute other points in the bin
  for(unsigned long long i = 0; i < pow((double)NrBins_h,(double)dimension); i++){
    
    int fPoint[dimension];
    for(int k = 0; k < dimension; k++)
      fPoint[k] = 0;
    
    //the bin value
    double binHight = 0;
    double sqrBin = 0;
    double maxBin = 0, minBin = 32000; //sper ca maximum nu este atit !!!!!!!!!!!!!!!!!!!!

    for(unsigned long long iBin = 0; iBin < pow((double)NrDivBin,(double)dimension); iBin++){

      //compute the fitness of the function 
      for(int index = 0; index < dimension; index++)
	x_current.element(index,0) = -(double)lenght/2.0 + point[index] * (double)lenght/(double)NrBins_h + fPoint[index] * (double) lenght / (double)(NrBins_h * NrDivBin);
      
      /* cout << "bin " << i << " pointBin (";
      for(int indexTemp = 0; indexTemp < dimension; indexTemp++)
	cout << point[indexTemp] << ",";
      cout<< ") point " << iBin << " pointInBin (";
      for(int indexTemp = 0; indexTemp < dimension; indexTemp++)
	cout << fPoint[indexTemp] << ",";      
	cout<< ") indiv "<< x_current << "\n";*/

      double tempValue = fitness(x_current,covarianceMatrix,covarianceInverse,positionMatrix, dimension,mixture);

      if(tempValue >= MIN_FITNESS){
	cov12 += x_current.element(0,0) * x_current.element(1,0) * tempValue;
	mean1 += x_current.element(0,0);// * tempValue;
	mean2 += x_current.element(1,0);// * tempValue;
	//cout << "x_current[0] "<< x_current.element(0,0) << " x_current[1] " << x_current.element(1,0) 
	//<< " cov12 "<< cov12 << " mean1 " << mean1 << " mean2 " << mean2 << "\n";
	sqrBin += tempValue;
	totalHigh++;
      }
	
      binHight += tempValue;
      if(tempValue > maxBin){
	maxBin = tempValue;
      }

      if(tempValue < minBin){
	minBin = tempValue;
      }

      int k = dimension - 1;
      while(k >= 0){
	fPoint[k]++;
	if(fPoint[k] < NrDivBin)
	  break;
	else {
	  fPoint[k] = 0;
	  k--;
	}
      }//end while k

    }// end for iBin


    total += binHight;
    sqrTotal+= sqrBin;
    //write the data in a file
    binHight = binHight/pow((double)NrDivBin,(double)dimension);

    for(int iFile = 0; iFile < dimension; iFile++)
      myFile << x_current.element(iFile,0) << "\t";
    myFile << binHight << "\t" << maxBin << "\t" << minBin << "\t"; 

    kl_distance += log(binHight) * (binHight);

    if(maxBin - minBin > maxDiff) 
      maxDiff = maxBin - minBin;

    //increment the point - start from the last position and go towards the first one
    int j = dimension - 1;
    while(j >= 0){
      point[j]++;
      if(point[j] < NrBins_h){
	myFile << "\n";
	break;
      }
      else {
	point[j] = 0;
	j--;
	if(j == dimension -2) myFile << "\n";
      }
    }//end while j
  } // end for i 

  cout << "maxDiff" << maxDiff << "\n";
  cout << " total " <<  total/ pow((double) NrDivBin * NrBins_h,(double)dimension) << " cov12 " << cov12
    << " mean1 " << mean1 / pow((double) NrDivBin * NrBins_h,(double)dimension)<< " mean2 " << mean2 / pow((double) NrDivBin * NrBins_h,(double)dimension)<< "\n";
    // <<" mean1 " << mean1/sqrTotal << " mean2 " << mean2/sqrTotal << "\n";
  
  double std1 = 0, std2 = 0;
  //for all possible points that you can think of
  //for each of this points compute other points in the bin
  cov12 = 0;
  for(unsigned long long i = 0; i < pow((double)NrBins_h,(double)dimension); i++){
    
    int fPoint[dimension];
    for(int k = 0; k < dimension; k++)
      fPoint[k] = 0;
    
    //the bin value
    double binHight = 0;
    double maxBin = 0, minBin = 32000; //sper ca maximum nu este atit !!!!!!!!!!!!!!!!!!!!

    for(unsigned long long iBin = 0; iBin < pow((double)NrDivBin,(double)dimension); iBin++){

      //compute the fitness of the function 
      for(int index = 0; index < dimension; index++)
	x_current.element(index,0) = -(double)lenght/2.0 + point[index] * (double)lenght/(double)NrBins_h + fPoint[index] * (double) lenght / (double)(NrBins_h * NrDivBin);
      
      double tempValue = fitness(x_current,covarianceMatrix,covarianceInverse,positionMatrix, dimension,mixture);

      /* cout << "bin " << i << " pointBin (";
      for(int indexTemp = 0; indexTemp < dimension; indexTemp++)
	cout << point[indexTemp] << ",";
      cout<< ") point " << iBin << " pointInBin (";
      for(int indexTemp = 0; indexTemp < dimension; indexTemp++)
	cout << fPoint[indexTemp] << ",";      
	cout<< ") indiv "<< x_current << "\n";*/
      if(tempValue >= MIN_FITNESS){
	//cov12 += (x_current.element(0,0) -  mean1 / totalHigh) * (x_current.element(1,0) -  mean2 / totalHigh) * tempValue;
	cov12 += (x_current.element(0,0) -  mean1 / sqrTotal)* (x_current.element(1,0) -  mean2 / sqrTotal) * tempValue;
	std1 += pow(x_current.element(0,0) -  mean1 / sqrTotal,2.0) * tempValue;
	std2 += pow(x_current.element(1,0) -  mean2 / sqrTotal,2.0) * tempValue;
	//cov12 += x_current.element(0,0) * x_current.element(1,0);
	//cout << "x_current[0] "<< x_current.element(0,0) << " x_current[1] " << x_current.element(1,0) << " cov12 "<< cov12 << "\n";
      }

      binHight += tempValue;

      if(tempValue > maxBin){
	maxBin = tempValue;
      }

      if(tempValue < minBin){
	minBin = tempValue;
      }

      int k = dimension - 1;
      while(k >= 0){
	fPoint[k]++;
	if(fPoint[k] < NrDivBin)
	  break;
	else {
	  fPoint[k] = 0;
	  k--;
	}
      }//end while k

    }// end for iBin

     //write the data in a file
    //binHight = binHight/total;

    for(int iFile = 0; iFile < dimension; iFile++)
      histFile << (int)(x_current.element(iFile,0) * NrBins_h/lenght + NrBins_h/2) << "\t";
    histFile << binHight/pow((double)NrDivBin,(double)dimension) << "\t"<< binHight/total << "\t" << maxBin << "\t" << maxBin/total << "\t" 
	     << minBin << "\t" << minBin/total << "\n"; 

    if(maxBin - minBin > maxDiff) 
      maxDiff = maxBin - minBin;

    //increment the point - start from the last position and go towards the first one
    int j = dimension - 1;
    while(j >= 0){
      point[j]++;
      if(point[j] < NrBins_h){
	break;
      }
      else {
	point[j] = 0;
	j--;
      }
    }//end while j
  } // end for i 

  cout << "\n\n# maximum difference " << maxDiff << " total " << sqrTotal / totalHigh;
  cout << " corrected conv " << cov12/sqrTotal << " std1 " << std1/sqrTotal << " std2 " << std2/sqrTotal << " cor12" <<  cov12/sqrt(std2 * std1) << "\n";

  myFile.flush();
  myFile.close();

  histFile.flush();
  histFile.close();


  return kl_distance;
}
