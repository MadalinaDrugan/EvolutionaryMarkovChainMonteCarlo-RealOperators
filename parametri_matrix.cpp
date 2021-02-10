//write the parametri file for the mixture of gaussian
//parametri: bounds
//           dimension of the problem
//           number of mixture
//           position of mixture 
//           covariance matrix for each individual

#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
                                     // newmatap.h will get include.h
#include "MatrixUtil\newmatap.h"                // need matrix applications
#include "MatrixUtil\newmatio.h"                // need matrix output routines

#include<iostream.h>
#include <fstream.h>
#include <stdlib.h>

//#include <math.h>
#include"Utils\random.cpp"
#include "Utils\histogram.cpp"                // need matrix applications

#include "matrix_generation.cpp"                // need matrix applications
#include "test_matrix_generation.cpp"                // need matrix applications

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

void readCovariance(Matrix* covarianceMixture, Matrix* positionMixture, int dimension, int lenght, int mixture){
  char buffer[100];
  int typeLandscape;
  ifstream myFile("parametri_mixture.file");		
  if(!myFile.is_open()){
    cout << "Error opening file \n";
    exit(1);
  }

//type of lanscape
  myFile >> buffer;
  typeLandscape = atoi(buffer);
  cout << buffer << " should be" << typeLandscape << "\n";
  myFile.getline(buffer,100);
  cout << buffer << "\n";
  //cout << " discard " << buffer << "\n";

//number of components
//primel cimp este explicative
  myFile.getline(buffer,100);
  cout << " discard " << buffer << "\n";
  myFile.getline(buffer,100);
  //mixture = atoi(buffer);
  cout << buffer << " should be " << mixture << "\n";

//number of dimensions
  myFile.getline(buffer,100);
  cout << " discard " << buffer << "\n";
  myFile.getline(buffer,100);
  cout << buffer << "should be" << atoi(buffer) << "\n";
  //mixture = atoi(buffer);
  
//the bounds of the landscape
  myFile.getline(buffer,100);
  cout << " discard " << buffer << "\n";
  myFile.getline(buffer,100);
  cout << " discard " << buffer << "should be " <<  atoi(buffer);

//if we have correation in a matrix
  if(typeLandscape == 1 || typeLandscape == 2) {
    //position is always in the middle of the space 0
    cout << " enter correctly on typeLandscape " << typeLandscape << "\n";
      for(int i = 0; i < dimension; i++)
	  positionMixture[0].element(i,0) = 0;
      
      //read the title of the covariance matrix
      myFile.getline(buffer,100);
      cout << "discard last before matrix" << buffer << "\n";

      //read the non-correlated covariance matrix
      for(int i = 0; i < dimension; i++){
	  for(int j = 0; j < dimension; j++){
	      myFile >> buffer;
	      if(typeLandscape == 2){
		  covarianceMixture[0].element(i,j) = atof(buffer);
		  cout << "(" << i << "," << j << "," << covarianceMixture[0].element(i,j) << ")\n";
	      } else 
		cout << "(" << i << "," << j << "," << buffer << ")\n";
	  }
	  cout << "\n";
      }

      //read the angles title and angles
      myFile.getline(buffer,100);
      cout << " discard " << buffer << "\n";
      myFile.getline(buffer,100);
      cout << " discard " << buffer << "\n";

      //read the title of the corelated covariance matrix
      myFile.getline(buffer,100);
      cout << " discard " << buffer << "\n";
      myFile.getline(buffer,100);
      cout << " discard " << buffer << "\n";

      for(int i = 0; i < dimension; i++)
	  for(int j = 0; j < dimension; j++){
	      myFile >> buffer;
	      if(typeLandscape == 1){
		  covarianceMixture[0].element(i,j) = atof(buffer);
		  cout << "(" << i << "," << j << "," << covarianceMixture[0].element(i,j) << ")\n";
	      } else 
		cout << "(" << i << "," << j << "," << buffer << ")\n";
	  }

      cout << " type lanscape " << typeLandscape << " is " << covarianceMixture[0] << "\n";
  } else {
    // extract which dimension
    myFile.getline(buffer,100); 
    for(int i = 0; i < mixture; i++){
      
      myFile.getline(buffer,100); 
      //cout << " numerical position[" << j << "]= "  << buffer << "\n";
      for(int k = 0; k < dimension; k++){
	myFile >> buffer;
	positionMixture[i].element(k,0) = atof(buffer);
      }
      
      myFile.getline(buffer,100); 
      //cout << " numerical covariance[" << j << "][" << j << "]= " << buffer << "\n";
      for(int k = 0; k < dimension; k++)
	for(int t = 0; t < dimension; t++){
	  myFile >> buffer;
	  covarianceMixture[i].element(k,t) = atof(buffer);
	}
      
      }
  }
  
  myFile.close();

}

///// main function
int main (){
    setSeed();
    
//bounds
    int lenght = 8; 
//dimensions
    int dimension = 3;
//number of mixtures
    int mixturePar = 1;
    
    //covariance matrix
    Matrix* covarianceMatrix = new Matrix[mixturePar];
    Matrix* correlations = new Matrix[mixturePar];
    Matrix* positionMatrix = new Matrix[mixturePar];
    for(int i = 0; i < mixturePar; i++){
      covarianceMatrix[i].ReSize(dimension, dimension);
      correlations[i].ReSize(dimension, dimension);
      positionMatrix[i].ReSize(dimension,1);
    }
    
    //first generate a matrix 
    mixtureCorrelatedHigh(covarianceMatrix,positionMatrix,correlations,dimension,lenght,mixturePar);

    //cout << "maxfitness " << maxFitness(covarianceMatrix,dimension,mixturePar) << "\n";

    //read the matrix and generate the histogram
    //readCovariance(covarianceMatrix,positionMatrix,dimension,lenght,mixturePar);
    
    double kl_distance = histogram(covarianceMatrix,positionMatrix,dimension,lenght,mixturePar);

    ofstream myFile("parametri_mixture.file",ios::app);		
    if(!myFile.is_open()){
    	cout << "Error opening file \n";
    	exit(1);
    }

    myFile << "KL - distance \n";
    myFile << kl_distance << "\n";
    
    myFile.flush();
    myFile.close();
}

