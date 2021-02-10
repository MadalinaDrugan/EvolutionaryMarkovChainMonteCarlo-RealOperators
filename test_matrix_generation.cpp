// some test function to see for the generation of the matrix works
#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
                                     // newmatap.h will get include.h
#include<iostream.h>
#include <fstream.h>
#include <stdlib.h>
//#include <math.h>

#include "MatrixUtil\newmatap.h"                // need matrix applications
#include "MatrixUtil\newmatio.h"                // need matrix output routines

#include"Utils\random.cpp"
//#include "matrix_generation.cpp"                // need matrix applications

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

void test(Matrix* covarianceMatrix, Matrix* covarianceMatrixInverse, Matrix* positionMatrix, int dimensionMy, int mixtureMy){
  double value = 0;

  // generate a value of the search space equal with one of the first position
  Matrix x_current(dimensionMy,1);
  for(int i =0; i < dimensionMy; i++)
    if (i < 0) //if(genrand_real1() > 0.5)
      x_current.element(i,0) = positionMatrix[0].element(i,0) - 1.5;
    else x_current.element(i,0) = positionMatrix[0].element(i,0) + 1.5;
  //  for(int i =0; i < dimensionMy; i++)
  //x_current.element(i,0) = positionMatrix[0].element(i,0);
 
  Matrix multiM(1,1);

  //compute the values
  for(int i = 0; i < mixtureMy; i++){
    cout << " mixiture " << i << "\n";

    double tempValue = 0;

    multiM = (x_current - positionMatrix[i]).t() * covarianceMatrixInverse[i] * (x_current - positionMatrix[i]);

    cout << " position mixture " << positionMatrix[i].t() << " difference with current " << (x_current - positionMatrix[i]).t() 
	 << "\n half matrix " << (x_current - positionMatrix[i]).t() * covarianceMatrixInverse[i]
	 << " multiple with " << multiM << " determinant " << abs(covarianceMatrix[i].Determinant()) << "\n"; 

    tempValue = pow(2.0 * pi,-((double)dimensionMy)/2.0) * pow(abs(covarianceMatrix[i].Determinant()),-1.0/2.0) * exp(-0.5 * multiM.element(0,0));

    value += tempValue;

    cout << " exp " << exp(-0.5 * multiM.element(0,0)) << "value " << tempValue << "\n\n";
    
    //only for two dimensions
    double sigma1 = sqrt(covarianceMatrix[i].element(0,0));
    double sigma2 = sqrt(covarianceMatrix[i].element(1,1));
    double rho = covarianceMatrix[i].element(1,0)/(sigma1 * sigma2);
    double z = pow(x_current.element(0,0) - positionMatrix[i].element(0,0),2) / (sigma1* sigma1) - 2 * rho * (x_current.element(0,0) - positionMatrix[i].element(0,0)) * (x_current.element(1,0) - positionMatrix[i].element(1,0)) / (sigma1 * sigma2) + pow(x_current.element(1,0) - positionMatrix[i].element(1,0),2) / (sigma2* sigma2);
    double bivariate = (1.0/(2.0 * pi * sigma1 * sigma2 * sqrt(1.0 - rho*rho))) * exp(-z/(2.0 * (1 - rho * rho)));
    cout << "sigma_1 " << sigma1 << " sigma_2 " << sigma2 << " rho " << rho << " z " << z << " exp " << exp(-z/(2.0 * (1 - rho * rho))) 
	 << " bivariate " << bivariate << "\n"; 
  }

  cout << "fitness value " << value << "\n";

}

void test1(Matrix* covarianceMatrix, Matrix* covarianceMatrixInverse, Matrix* positionMatrix, int dimensionMy, int mixtureMy){
  double value = 0;

  // generate a value of the search space equal with one of the first position
  Matrix x_current(dimensionMy,1);
  for(int i =0; i < dimensionMy; i++)
    if(i >= 2)
      x_current.element(i,0) = positionMatrix[0].element(i,0) - 1.5;
    else x_current.element(i,0) = positionMatrix[0].element(i,0) + 1.5;
  Matrix multiM(1,1);
  cout << "position current " << x_current.t();

  //compute the values
  for(int i = 0; i < mixtureMy; i++){
    cout << " mixiture " << i << "\n";
    double tempValue = 0;
    multiM = (x_current - positionMatrix[i]).t() * covarianceMatrixInverse[i] * (x_current - positionMatrix[i]);

    cout << " position mixture " << positionMatrix[i].t() << " difference with current " << (x_current - positionMatrix[i]).t() 
	 << "\n half matrix " << (x_current - positionMatrix[i]).t() * covarianceMatrixInverse[i]
	 << " multiple with " << multiM << " determinant " << abs(covarianceMatrix[i].Determinant()) << "\n"; 

    tempValue = pow(2.0 * pi,-((double)dimensionMy)/2.0) * pow(abs(covarianceMatrix[i].Determinant()),-1.0/2.0) * exp(-0.5 * abs(multiM.element(0,0)));

    value += tempValue;

    cout << " exp " << exp(-0.5 * multiM.element(0,0)) << "value " << tempValue << "\n\n";
  }

  cout << "fitness value " << value << "\n";

}
