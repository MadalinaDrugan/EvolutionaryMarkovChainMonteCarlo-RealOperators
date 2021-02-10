//generated the matrix with covariance 
#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
                                     // newmatap.h will get include.h
#include<iostream.h>
#include <fstream.h>
#include <stdlib.h>
//#include <math.h>

#include "MatrixUtil\newmatap.h"                // need matrix applications
#include "MatrixUtil\newmatio.h"                // need matrix output routines
//#include "MatrixUtil\newmat.h"                // need matrix applications

#include"Utils\random.cpp"

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif


void RotateCovariance(Matrix* covarianceYes, Matrix* positionMatrix, int dimension, int lenght, int mixture){
    Matrix covarianceNon(dimension,dimension);
    double angles[dimension];
    for(int i = 0; i < mixture; i++)
      positionMatrix[i].element(i,0) = 0;

    //open the file to be wrote
    ofstream myFile("parametri_mixture.file");		
    if(!myFile.is_open()){
	cout << "Error opening file \n";
	exit(1);
    }

    //write general data
    myFile << "1" << "\t" << "rotate\n";
    myFile << "nr_of_mixtures\n";
    myFile << mixture << "\n";
    myFile << "dimension\n";
    myFile << dimension << "\n";
    myFile << "lenght" << "\n" << lenght << "\n";

    //generate at random the non-correlated matrix
    //only the diagonal
    for(int i = 0; i < dimension; i++)
      covarianceNon.element(i,i) = genrand_real1(); // * lenght - lenght/2.0;
    myFile << " matrix without correlation in dimensions \n";
    myFile << covarianceNon;

    //generate the corrltated matrix
    //generate an angle to rotate the non-correlated matrix
    //high angles
    myFile << " rotation\n";
    cout << " rotation angle where pi = " << pi << "\n";
    for(int i = 0; i < (int)(dimension * (dimension - 1) / 2.0); i++){
      //double randomN = genrand_real1(); //gaussian();
      angles[i] = (genrand_int32() % 10 + 39)* (double)pi/180.0; //do not take 90 degrees - > exception in computing the matrices
      //angles[i] = genrand_real1()* (double)pi/2.0;
      cout << " angles[" << i << "] = " << (angles[i]/pi)*180 << " is " << cos(angles[i]) << " \t " << sin(angles[i])<< "\n";
      if(i != dimension * (dimension - 1) / 2.0 -1) myFile << angles[i] << "\t";
      else myFile << angles[i] << "\n";
    }

    //rotate the matrix on the dimensions
    Matrix temp[(int)(dimension * (dimension - 1)/2.0)];
    for(int i = 0; i < dimension*(dimension - 1)/2.0; i++)
      temp[i].ReSize(dimension,dimension);
    
    for(int i = 0; i < dimension*(dimension - 1)/2.0; i++)
      for(int j = 0; j < dimension; j++)
	for(int k = 0; k < dimension; k++)
	  temp[i].element(j,k) = 0;
    
    for(int i = 0; i < dimension*(dimension - 1)/2.0; i++)
      for(int j = 0; j < dimension; j++)
	temp[i].element(j,j) = 1;
    
    //cout << " dimension " << 
    
    if(dimension == 2){
      temp[0].element(0,0) = cos(angles[0]);
      temp[0].element(1,1) = cos(angles[0]);
      temp[0].element(0,1) = sin(angles[0]);
      temp[0].element(1,0) = -sin(angles[0]);
      
      covarianceYes[0] = temp[0];
    }
    else 
	if(dimension >= 3){
	  
	  temp[0].element(1,1) = cos(angles[0]);
	  temp[0].element(2,2) = cos(angles[0]);
	  temp[0].element(1,2) = sin(angles[0]);
	  temp[0].element(2,1) = -sin(angles[0]);
	  
	  temp[1].element(0,0) = cos(angles[1]);
	  temp[1].element(2,2) = cos(angles[1]);
	  temp[1].element(0,2) = -sin(angles[1]);
	  temp[1].element(2,0) = sin(angles[1]);
	  
	  temp[2].element(0,0) = cos(angles[2]);
	  temp[2].element(1,1) = cos(angles[2]);
	  temp[2].element(1,0) = sin(angles[2]);
	  temp[2].element(0,1) = -sin(angles[2]);
	  
	  if(dimension == 3) 
	    covarianceYes[0] = temp[0] * temp[1] * temp[2]; // * temp[1].i() * temp[0].i();
	  else
	    if(dimension >= 4){
	      
	      temp[0].element(0,0) = cos(angles[0]);
	      temp[0].element(1,1) = cos(angles[0]);
	      temp[0].element(0,1) = sin(angles[0]);
	      temp[0].element(1,0) = -sin(angles[0]);
	      
	      temp[1].element(1,1) = cos(angles[1]);
	      temp[1].element(2,2) = cos(angles[1]);
	      temp[1].element(1,2) = sin(angles[1]);
	      temp[1].element(2,1) = -sin(angles[1]);
	      
	      temp[2].element(0,0) = cos(angles[2]);
	      temp[2].element(2,2) = cos(angles[2]);
	      temp[2].element(2,0) = -sin(angles[2]);
	      temp[2].element(0,2) = sin(angles[2]);

	      temp[3].element(0,0) = cos(angles[3]);
	      temp[3].element(3,3) = cos(angles[3]);
	      temp[3].element(0,3) = sin(angles[3]);
	      temp[3].element(3,0) = -sin(angles[3]);
	      
	      temp[4].element(1,1) = cos(angles[4]);
	      temp[4].element(3,3) = cos(angles[4]);
	      temp[4].element(1,3) = -sin(angles[4]);
	      temp[4].element(3,1) = sin(angles[4]);

	      temp[5].element(2,2) = cos(angles[5]);
	      temp[5].element(3,3) = cos(angles[5]);
	      temp[5].element(2,3) = -sin(angles[5]);
	      temp[5].element(3,2) = sin(angles[5]);

	      if(dimension == 4)
		covarianceYes[0] = temp[1] * temp[2] * temp[0] * temp[3] * temp[5] * temp[4] * temp[5].i() * temp[3].i() * temp[0].i() * temp[2].i() * temp[1].i();
	      else 
		if(dimension == 5){
		  
		}
	    }
	} else {
	  //experimental ; not verified yet
	  
	  /* for(int i =0; i < dimension * (dimension - 1)/2.0; i++)
	    for(int i1 = i+1; i1 < dimension * (dimension - 1)/2.0; i1++){
	      //for each element from the rotation matrix
	      for(int j = 0; j < dimension; j++)
		for(int k = 0; k < dimension; k++)
		  if(j == k && (j == 0 || k == i+1))
		    temp.element(j,k) = cos(angles[i]);
		  else 
		    if(j == k) temp.element(j,k) = 1;
		    else 
		      if(k == 0 && j == i + 1) temp.element(j,k) = sin(angles[i]);
		      else 
			if(j == 0 && k == i + 1) temp.element(j,k) = -sin(angles[i]);
			else temp.element(j,k) = 0;
	      
	      //transfer the resuted matrix
	      cout << " rotation for dimension " << i << " is \n" << temp;
	      if(i == 0)	
		covarianceYes = temp;
	      else covarianceYes *= temp;
	      cout << " resulting rotation \n" << covarianceYes; 
	      }*/
	}
    //matrices to compute the rotations
    cout << "covariance_no \n" << covarianceNon << " rotation matrix \n " << covarianceYes[0] << "\n";
    covarianceYes[0]  = covarianceYes[0] * covarianceNon * covarianceYes[0].t();
    //cout << " half covariance \n " << covarianceYes << "\n";
    //covarianceYes *= covarianceYes.t();
    cout << " covariance \n " << covarianceYes[0] << "\n";
    myFile << "covariance matrix - rotated matrix \n";
    myFile << covarianceYes[0];

    cout << " determinant " << covarianceYes[0].Determinant() << " = " << covarianceNon.Determinant() << "\n";

    myFile.flush();
    myFile.close();
}


////////////mixture
//global variables for indicate the position in a cluster
//componenets in a mixture are generated according with a distribution 
//how to generate positions with this clusters
int clusterIndex = 0; //position in the cluster
int clusterCum = 0; //number of generated individuals from a position
//clusters to generate the positions
int clusters = 3;
int* clusterElements;
Matrix* clusterPosition;
Matrix* clusterCovariance;

void mixturePositionsCorrelatedInitialization(int dimension, int lenght, int mixturePar){
  clusterElements = new int[clusters];
    //manual choise summed up to 20
    clusterElements[0] = (int)((mixturePar*10.0)/100.0);
    clusterElements[1] = (int)((mixturePar *30.0)/100.0);
    clusterElements[2] = (int)(mixturePar * (1 - clusterElements[0] - clusterElements[1]));
    
    clusterPosition = new Matrix[clusters];
    clusterCovariance = new Matrix[clusters];
    for(int i = 0; i < clusters; i++){
	clusterPosition[i].ReSize(dimension,1);
	clusterCovariance[i].ReSize(dimension,dimension);
    }
    
    //generate random the positions of the three clusters
    for(int i = 0; i < clusters; i++){
	for(int j = 0; j < dimension; j++)
	    clusterPosition[i].element(j,0) = genrand_real1() * lenght - lenght/2;
	
	//assume that number of multi-variates in a cluter is increasing with i - thus the covariance is decreasing 
	for(int j = 0; j < dimension; j++)
	    for(int k = j; k < dimension; k++){
		//I love Madalina!!! (brown bear)
		if(i == 0) clusterCovariance[i].element(j,k) = gaussian() * 0.2 + 0.8;
		else clusterCovariance[i].element(j,k) = gaussian() * 0.6 + 0.4;
		
		if(clusterCovariance[i].element(j,k) <= 0.01) clusterCovariance[i].element(j,k) = 0.01;
		if(clusterCovariance[i].element(j,k) >= 0.99) clusterCovariance[i].element(j,k) = 0.99;
	    }
	
	for(int j = 0; j < dimension; j++)
	    for(int k = j; k < dimension; k++)
		if(j != k){
		    clusterCovariance[i].element(j,k) = sqrt(clusterCovariance[i].element(j,j)) * sqrt(clusterCovariance[i].element(k,k)) * clusterCovariance[i].element(j,k);
		    clusterCovariance[i].element(k,j) = clusterCovariance[i].element(j,k);
		}
    }

    //initialization of the indexes
    clusterCum = 0;
    clusterIndex = 0;
}

Matrix mixturePositionsCorrelatedIteration(int dimension, int lenght, int mixturePar){
    if(clusterCum >= clusterElements[clusterIndex]){
	clusterIndex++;
	clusterCum = 0;
    } else clusterCum++;

    Matrix dataTemp(dimension,1);
    for(int i = 0; i < dimension; i++)
      dataTemp.element(i,0) = gaussian();

    return clusterPosition[clusterIndex] + dataTemp * clusterCovariance[clusterIndex]; 
}

void mixturePositionsCorrelatedDelete(){
  delete[] clusterElements;
  delete[] clusterPosition;
  delete[] clusterCovariance;
}

void mixtureStatistics(int dimension, int mixturePar, Matrix* positionMatrix, Matrix* covarianceMatrix){
  //write a file with statistics about the mixtures 
  ofstream myFile_S("parametri_mixture.stat");		
  if(!myFile_S.is_open()){
    cout << "Error opening file \n";
    exit(1);
  }

  //write in this file in a shape easier to read the
  myFile_S << "positions in a mixture\n";
  for(int i = 0; i < mixturePar; i++)
    myFile_S << "postion component " << i << " = " << positionMatrix[i].t() << "\n";

  //position mixture the distance computation 
  Matrix euclidianDistance(mixturePar,mixturePar);

  double maxim = 0; //, iMax = -1, jMax = -1;
  //double maximWeight = 0, iMaxW = -1;
  //double minim = 300, iMin = -1, jMin = -1;
  //double minimWeight = 300, iMinW = -1;

  int indexOrdered = 0;
  double orderedDistance[(mixturePar-1) * mixturePar/2];
  int iOrderedDistance[(mixturePar-1) * mixturePar/2];
  int jOrderedDistance[(mixturePar-1) * mixturePar/2];

  for(int i = 0; i < mixturePar; i++)
    for(int j = i; j < mixturePar; j++){
      euclidianDistance.element(i,j) = 0;
      if(i != j){
	for(int k = 0; k < dimension; k++)
	  euclidianDistance.element(i,j) += pow(positionMatrix[i].element(k,0) - positionMatrix[j].element(k,0),2.0);
	euclidianDistance.element(i,j) = sqrt(euclidianDistance.element(i,j));
	euclidianDistance.element(j,i) = euclidianDistance.element(i,j);
	orderedDistance[indexOrdered] = euclidianDistance.element(i,j);
	iOrderedDistance[indexOrdered] = i;
	jOrderedDistance[indexOrdered] = j;
	indexOrdered++;
      }
      //if(euclidianDistance.element(i,j) < minim){
      //	minim = euclidianDistance.element(i,j);
      //iMin = i;
      //jMin = j;
      //}
  }

  myFile_S << " euclidian distance between positions \n";
  myFile_S << euclidianDistance << "\n";

  //sort the distance in the mixture
  for(int i = 0; i < (mixturePar-1) * mixturePar/2; i++)
    for(int j = 0; j < (mixturePar-1) * mixturePar/2; j++)
      if(orderedDistance[i] > orderedDistance[j]){
	maxim = orderedDistance[i];
	orderedDistance[i] = orderedDistance[j];
	orderedDistance[j] = (int)maxim;

	maxim = iOrderedDistance[i];
	iOrderedDistance[i] = iOrderedDistance[j];
	iOrderedDistance[j] = (int)maxim;

	maxim = jOrderedDistance[i];
	jOrderedDistance[i] = jOrderedDistance[j];
	jOrderedDistance[j] = (int)maxim;
      }

  myFile_S << " ordered Euclidian distances \n"; 
  for(int i = 0; i < (mixturePar-1) * mixturePar/2; i++)
    myFile_S << "(" << i << "," << orderedDistance[i] << "," << iOrderedDistance[i] << "," << jOrderedDistance[i] << ")\t";
  myFile_S << "\n\n";

  myFile_S << "ordered weighted Euclidian distances \n"; 

  for(int i = 0; i < mixturePar; i++){
    double weight = 0;
    for(int j = 0; j < mixturePar; j++)
      weight += euclidianDistance.element(i,j);
    weight /= mixturePar;
    orderedDistance[i] = weight;
    iOrderedDistance[i] = i;

  }

  //sort the distance in the mixture
  for(int i = 0; i < mixturePar; i++)
    for(int j = 0; j < mixturePar; j++)
      if(orderedDistance[i] > orderedDistance[j]){
	maxim = orderedDistance[i];
	orderedDistance[i] = orderedDistance[j];
	orderedDistance[j] = (int)maxim;

	maxim = iOrderedDistance[i];
	iOrderedDistance[i] = iOrderedDistance[j];
	iOrderedDistance[j] = (int)maxim;
      }

  for(int i = 0; i < mixturePar; i++)
    myFile_S << "(" << i << "," << orderedDistance[i] << "," << iOrderedDistance[i] << ")\t";
  myFile_S << "\n\n";


  myFile_S << " covariance Matrixes of the mixture \n";
  for(int i = 0; i < mixturePar; i++)
    myFile_S << "covariance component " << i << " \n " << covarianceMatrix[i] << "\n";
  
  myFile_S << " the value of the function in mean for each mixture \n";
  Matrix values(mixturePar,1);
  //Matrix multiM(1,1); 
  for(int i =0; i < mixturePar; i++){
    //multiM = (positionMatrix[i] - positionMatrix[i]).t() * covarianceMatrixInverse[i] * (positionMatrix[i] - positionMatrix[i]);
    values.element(i,0) = pow(2.0 * pi,-((double)dimension)/2.0) * pow(abs(covarianceMatrix[i].Determinant()),-1.0/2.0);
  }
  myFile_S << values.t() << "\n";

  myFile_S << " the difference in hight between picks \n";
  for(int i = 0; i < mixturePar; i++)
    for(int j = i; j < mixturePar; j++){
      euclidianDistance.element(i,j) = values.element(i,0) - values.element(j,0);
  }
  myFile_S << euclidianDistance;  

  //close the file
  myFile_S.flush();
  myFile_S.close();
}

void mixtureCorrelated(Matrix* covarianceMatrix, Matrix* positionMatrix, int dimension, int lenght, int mixturePar){
    //minitialize a matrix like this to compute its determinant and its inverse
    //open the file to be wrote
    ofstream myFile("parametri_mixture.file");		
    if(!myFile.is_open()){
	cout << "Error opening file \n";
	exit(1);
    }
    
    //write general data
    myFile << "1\trotate\n";
    myFile << "nr_of_mixtures\n";
    myFile << mixturePar << "\n";
    myFile << "dimension \n";
    myFile << dimension << "\n";
    myFile << "lenght" << "\n" << lenght << "\n";

    ////initializate the distribution for generating position and covariance
    mixturePositionsCorrelatedInitialization(dimension, lenght, mixturePar);
    
    //position and covariance of a mixture
    for(int i = 0; i < mixturePar; i++){  
	//the number of current componenet of the mixture
	myFile << "mixture " << i << "\n";
	
	//we need correlated and un-correlated components
	int low = 0;
	if(genrand_real1() > 0.5){
	    low = 1;
	}
        
	//generates position following a distribution
	//change the writer of the distribuion 
	positionMatrix[i] = mixturePositionsCorrelatedIteration(dimension,lenght,mixturePar);//genrand_real1() * lenght - lenght/2.0;
	
	myFile << "position \n";
	myFile << positionMatrix[i];
	
	for(int j = 0; j < dimension; j++){
	    //myFile << " covariance matrix \n"; 
	    //variance between (0,1) with the same rules as for the rest 
	    //covarianceMatrix[i].element(j,j) = genrand_real3(); //mixtureCovarianceCorrelatedIteration(dimension,lenght,mixture);
	    if(low == 1){
		covarianceMatrix[i].element(j,j) = gaussian() * 0.2 + 0.2;
	    } else {
		covarianceMatrix[i].element(j,j) = gaussian() * 0.2 + 0.8;	  
	    }

	    if(covarianceMatrix[i].element(j,j) <= 0.01) covarianceMatrix[i].element(j,j) = 0.01;
	    if(covarianceMatrix[i].element(j,j) >= 0.99) covarianceMatrix[i].element(j,j) = 0.99;
	    
	}
	
	for(int j = 0; j < dimension; j++){
	    //covariance: half with low values half with high values
	    for(int k = j+1; k < dimension ; k++){
		//low values
		//covarianceMatrix.element([j][k] = gaussian();
		
	      if(low == 1){
		covarianceMatrix[i].element(j,k) = gaussian() * 0.2 + 0.2;
	      } else {
		covarianceMatrix[i].element(j,k) = gaussian() * 0.2 + 0.8;	  
	      }
	      
	      if(covarianceMatrix[i].element(j,k) <= 0.01) covarianceMatrix[i].element(j,k) = 0.01;
	      if(covarianceMatrix[i].element(j,k) >= 0.99) covarianceMatrix[i].element(j,k) = 0.99;
	

	      covarianceMatrix[i].element(j,k) = sqrt(covarianceMatrix[i].element(j,j)) * sqrt(covarianceMatrix[i].element(k,k)) * covarianceMatrix[i].element(j,k);
	      covarianceMatrix[i].element(k,j) = covarianceMatrix[i].element(j,k);
	      //covarianceMatrix.element([k][j] = covarianceMatrix.element([j][k];
	    }// end k
	} // end j

	myFile << "covariance\n";
	myFile << covarianceMatrix[i];
    }//and i - counter for mixture

    //only for two values
    //test(covarianceMatrix,covarianceMatrixInverse,positionMatrix,dimension,mixturePar);
    //test1(covarianceMatrix,covarianceMatrixInverse,positionMatrix,dimension,mixturePar);
    
    //close the file
    mixturePositionsCorrelatedDelete();

    myFile.flush();
    myFile.close();
}

void mixtureCorrelatedHigh(Matrix* covarianceMatrix, Matrix* positionMatrix, Matrix* correlations, int dimension, int lenght, int mixturePar){
    //minitialize a matrix like this to compute its determinant and its inverse
    //open the file to be wrote
    ofstream myFile("parametri_mixture.file");		
    if(!myFile.is_open()){
	cout << "Error opening file \n";
	exit(1);
    }
    
    //write general data
    if(mixturePar == 1)
      myFile << "3\tcorrelations\n";
    else myFile << "4\tcorrelations\n";
    myFile << "nr_of_mixtures\n";
    myFile << mixturePar << "\n";
    myFile << "dimension \n";
    myFile << dimension << "\n";
    myFile << "lenght" << "\n" << lenght << "\n";

    ////initializate the distribution for generating position and covariance
    //mixturePositionsCorrelatedInitialization(dimension, lenght, mixturePar);
    
    //position and covariance of a mixture
    for(int i = 0; i < mixturePar; i++){  
	//the number of current componenet of the mixture
	//myFile << "mixture " << i << "\n";
	
	//we need correlated and un-correlated components
	int low = 0;
	//if(genrand_real1() > 0.5){
	//  low = 1;
	//}

	if(mixturePar == 1){        
	//generates position following a distribution
	//change the writer of the distribuion 
	positionMatrix[i] = 0;
//mixturePositionsCorrelatedIteration(dimension,lenght,mixturePar);//genrand_real1() * lenght - lenght/2.0;
	
	//myFile << "position \n";
	//myFile << positionMatrix[i];
#define hand_set	
#ifdef hand_set
	if(dimension == 3){
	  covarianceMatrix[i].element(0,0) = 1;
	  correlations[i].element(0,0) = 1;
	  
	  covarianceMatrix[i].element(1,1) = 0.1;
	  correlations[i].element(1,1) = 1;
	  
	  covarianceMatrix[i].element(2,2) = 0.01;
	  correlations[i].element(2,2) = 1;

	  correlations[i].element(0,1) = 0.9;
	  correlations[i].element(1,0) = 0.9;

	  correlations[i].element(0,2) = 0.8;
	  correlations[i].element(2,0) = 0.8;
 
	  correlations[i].element(1,2) = 0.9;
	  correlations[i].element(2,1) = 0.9;
	} else 	if(dimension == 2){
	  covarianceMatrix[i].element(0,0) = 1;
	  correlations[i].element(0,0) = 1;
	   
	  covarianceMatrix[i].element(1,1) = 0.01;
	  correlations[i].element(1,1) = 1;
	  
	  correlations[i].element(0,1) = 0.9;
	  correlations[i].element(1,0) = 0.9;

	} else if(dimension == 5){
	  covarianceMatrix[i].element(0,0) = 1;
	  correlations[i].element(0,0) = 1;
	  
	  covarianceMatrix[i].element(1,1) = 0.1;
	  correlations[i].element(1,1) = 1;
	  
	  covarianceMatrix[i].element(2,2) = 0.05;
	  correlations[i].element(2,2) = 1;

	  covarianceMatrix[i].element(3,3) = 0.01;
	  correlations[i].element(3,3) = 1;

	  covarianceMatrix[i].element(4,4) = 0.005;
	  correlations[i].element(3,3) = 1;

	  correlations[i].element(0,1) = 0.9;
	  correlations[i].element(1,0) = 0.9;

	  correlations[i].element(0,2) = 0.8;
	  correlations[i].element(2,0) = 0.8;

	  correlations[i].element(0,3) = 0.5;
	  correlations[i].element(3,0) = 0.5;

	  correlations[i].element(0,4) = 0.7;
	  correlations[i].element(4,0) = 0.7;
 
	  correlations[i].element(1,2) = 0.9;
	  correlations[i].element(2,1) = 0.9;
 
	  correlations[i].element(1,3) = 0.3;
	  correlations[i].element(3,1) = 0.3;

	  correlations[i].element(1,4) = 0.9;
	  correlations[i].element(4,1) = 0.9;

	  correlations[i].element(2,3) = 0.8;
	  correlations[i].element(3,2) = 0.8;

	  correlations[i].element(2,4) = 0.9;
	  correlations[i].element(4,2) = 0.9;

	  correlations[i].element(3,4) = 0.1;
	  correlations[i].element(4,3) = 0.1;
	} else{
	  cout << " dimension not implemented " << dimension << "\n";
	  exit(1);
	}
#else
	for(int j = 0; j < dimension; j++){
	    //myFile << " covariance matrix \n"; 
	    //variance between (0,1) with the same rules as for the rest 
	    //covarianceMatrix[i].element(j,j) = genrand_real3(); //mixtureCovarianceCorrelatedIteration(dimension,lenght,mixture);
	    if(genrand_real1() > 0.5){
		covarianceMatrix[i].element(j,j) = gaussian() * 0.1 + 0.1;
	    } else {
		covarianceMatrix[i].element(j,j) = gaussian() * 0.1 + 0.9;	  
	    }

	    if(covarianceMatrix[i].element(j,j) <= 0.01) covarianceMatrix[i].element(j,j) = 0.01;
	    if(covarianceMatrix[i].element(j,j) >= 0.99) covarianceMatrix[i].element(j,j) = 0.99;

	    correlations[i].element(j,j) = 1;
	}

	for(int j = 0; j < dimension; j++){
	    //covariance: half with low values half with high values
	    for(int k = j+1; k < dimension ; k++){
		//low values
		//covarianceMatrix.element([j][k] = gaussian();
		
	      if(low == 1){
		correlations[i].element(j,k) = gaussian() * 0.1 + 0.1;
	      } else {
		correlations[i].element(j,k) = gaussian() * 0.1 + 0.9;	  
	      }
	      
	      if(correlations[i].element(j,k) <= 0.01) correlations[i].element(j,k) = 0.01;
	      if(correlations[i].element(j,k) >= 0.99) correlations[i].element(j,k) = 0.99;
	      correlations[i].element(k,j) = correlations[i].element(j,k);
	    }// end k
	} // end j
#endif //hand_set	
	} else {
#define hand_set	
#ifdef hand_set
	if(dimension == 2){
	  if(i == 0){

	    covarianceMatrix[i].element(0,0) = 0.03;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.09;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = 0.1;
	    correlations[i].element(1,0) = 0.1;

	    positionMatrix[i].element(0,0) = 1.9;
	    positionMatrix[i].element(1,0) = -0.2;

	  } else if(i == 8){

	    covarianceMatrix[i].element(0,0) = 0.1;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.03;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = -0.6;
	    correlations[i].element(1,0) = -0.6;

	    positionMatrix[i].element(0,0) = 1.4;
	    positionMatrix[i].element(1,0) = 0.6;

	  } else if(i == 3){

	    covarianceMatrix[i].element(0,0) = 0.02;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.1;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = 0.5;
	    correlations[i].element(1,0) = 0.5;

	    positionMatrix[i].element(0,0) = 2.1;
	    positionMatrix[i].element(1,0) = -0.9;

	  } else if(i == 1){

	    covarianceMatrix[i].element(0,0) = 0.04;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.1;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = 0.5;
	    correlations[i].element(1,0) = 0.5;

	    positionMatrix[i].element(0,0) = -1.7;
	    positionMatrix[i].element(1,0) = -1.0;
	    
	  }else if(i == 10){

	    covarianceMatrix[i].element(0,0) = 0.02;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.13;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = 0.5;
	    correlations[i].element(1,0) = 0.5;

	    positionMatrix[i].element(0,0) = -1.4;
	    positionMatrix[i].element(1,0) = -0.2;
	    
	  } else if(i == 2){

	    covarianceMatrix[i].element(0,0) = 0.12;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.02;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = 0.5;
	    correlations[i].element(1,0) = 0.5;

	    positionMatrix[i].element(0,0) = 0.2;
	    positionMatrix[i].element(1,0) = 1.0;

	  }  else if(i == 4){

	    covarianceMatrix[i].element(0,0) = 0.02;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.13;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = 0.5;
	    correlations[i].element(1,0) = 0.5;

	    positionMatrix[i].element(0,0) = -0.8;
	    positionMatrix[i].element(1,0) = 0.6;
	    
	  }else if(i == 5){

	    covarianceMatrix[i].element(0,0) = 0.11;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.03;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = 0.7;
	    correlations[i].element(1,0) = 0.7;

	    positionMatrix[i].element(0,0) = 0.5;
	    positionMatrix[i].element(1,0) = -2.3;
	    
	  }else if(i == 6){

	    covarianceMatrix[i].element(0,0) = 0.1;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.04;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = -0.5;
	    correlations[i].element(1,0) = -0.5;

	    positionMatrix[i].element(0,0) = -0.5;
	    positionMatrix[i].element(1,0) = -2.3;
	    
	  }else if(i == 7){

	    covarianceMatrix[i].element(0,0) = 0.03;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.1;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = 0.5;
	    correlations[i].element(1,0) = 0.5;

	    positionMatrix[i].element(0,0) = -1.3;
	    positionMatrix[i].element(1,0) = -1.7;
	    
	  }else{

	    covarianceMatrix[i].element(0,0) = 0.04;
	    correlations[i].element(0,0) = 1;
	    
	    covarianceMatrix[i].element(1,1) = 0.1;
	    correlations[i].element(1,1) = 1;
	    
	    correlations[i].element(0,1) = 0.6;
	    correlations[i].element(1,0) = 0.6;

	    positionMatrix[i].element(0,0) = 1.5;
	    positionMatrix[i].element(1,0) = -1.9;

	    }
	}	  
#endif 
	}
 
	if(mixturePar != 1){
	  myFile << "position \n";
	  myFile << positionMatrix[i];
	}

	for(int j = 0; j < dimension; j++){
	  //covariance: half with low values half with high values
	  for(int k = j+1; k < dimension ; k++){
	    //low values
	    covarianceMatrix[i].element(j,k) = sqrt(covarianceMatrix[i].element(j,j)) * sqrt(covarianceMatrix[i].element(k,k)) * correlations[i].element(j,k);
	    covarianceMatrix[i].element(k,j) = covarianceMatrix[i].element(j,k);
	    //covarianceMatrix.element([k][j] = covarianceMatrix.element([j][k];
	  }// end k
	} // end j
	
	myFile << "correlations\n";
	myFile << correlations[i];

	myFile << "covariance\n";
	myFile << covarianceMatrix[i];
    }//and i - counter for mixture

    //only for two values
    //test(covarianceMatrix,covarianceMatrixInverse,positionMatrix,dimension,mixturePar);
    //test1(covarianceMatrix,covarianceMatrixInverse,positionMatrix,dimension,mixturePar);
    
    //close the file
    mixturePositionsCorrelatedDelete();

    myFile.flush();
    myFile.close();
}

