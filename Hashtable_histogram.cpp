#include<iostream.h> 
#include<fstream.h>
#include<stdlib.h>
#include<math.h>

#include"Hashtable_histogram.h"

using namespace std;

Hashtable_histogram :: Hashtable_histogram(unsigned long nr_vars, unsigned long bins)
{
  dimension = nr_vars;
  nrBins = bins;

  histogram = new Hashtable_histogram*[nrBins];  
  for(int i = 0; i < nrBins; i++)
    histogram[i] = NULL;

    mean = 0;
    inverseMean = 0;
    meanSqrt = 0;

    nrElem = 0;
    nrIndividuals = 0;

    //first element from the hashtable
    first = 0;

#ifdef FIRST_ITERATOR
    //irerators
    listCurrent = NULL;
    indexCurrent = NULL;
#endif //iterator

#ifdef SECOND_ITERATOR
    listSecondCurrent = NULL;
    indexSecondCurrent = NULL;
    //first = 1;
#endif //second iterator

    // cout << " dimension " << dimension << " nrBins " << nrBins << "\n";
 
}

Hashtable_histogram :: ~Hashtable_histogram()
{
    //cout << " ajuns aici \n";
    if(first == 0){
#ifdef FIRST_ITERATOR
	delete[] listCurrent;
	delete[] indexCurrent;
#endif //first iterator	 

#ifdef SECOND_ITERATOR
	//cout << "detete the second \n";
	if(listSecondCurrent != NULL)
	    delete[] listSecondCurrent;
	//cout << "delete second index \n";
	if(indexSecondCurrent != NULL)
	    delete[] indexSecondCurrent;
	//cout << "finish to delete \n";
#endif // second iterator	
    }	
}

void Hashtable_histogram :: free_table(Hashtable_histogram *table)
{
  if (table == NULL)
  return;

  for(unsigned long i = 0; i < dimension; i++)
    free_table(table->histogram[i]);
  
  if(first != 0){
      delete table;
      table = NULL;
      //delete listCurrent;
      //if(indexCurrent != NULL) delete indexCurrent;
  } else {
    for(unsigned long i = 0; i < dimension; i++)
      table->histogram[i] = NULL;
    
    nrIndividuals = 0;

    mean = 0;
    inverseMean = 0;
    meanSqrt = 0;

    nrElem = 0;
    dimension = 0;
    nrBins = 0;
  }

  //free(table);
  // free(listCurrent);
  //free(indexCurrent);

  //delete table;
  //table == NULL;
}


void Hashtable_histogram :: reset()
{
    if(this == NULL) {
	cout << "Hashtable_histogram is empty\n";
	return;
    }
#ifdef FIRST_ITERATOR
    if(listCurrent != NULL) 
      for(unsigned long i = 0; i < dimension; i++){
	listCurrent[i] = NULL;
	indexCurrent[i] = 0;

#ifdef SECOND_ITERATOR
	secondListCurrent[i] = NULL;
	secondListCurrent[i] = 0;
#endif //second iterator
      }
#endif // first iterator	

    for(unsigned long i = 0; i < dimension; i++)
      free_table(histogram[i]);
    
    for(unsigned long i = 0; i < dimension; i++)
      histogram[i] = NULL;

    nrElem = 0;
    nrIndividuals = 0;
    mean = 0;
    inverseMean = 0;
    meanSqrt = 0;

    dimension = 0;
    nrBins = 0;
}

unsigned long Hashtable_histogram :: check_individual(unsigned long* indiv, double& value, double& square, double& inverse)
{
  unsigned long i;
  unsigned long output;
  Hashtable_histogram *table;
  unsigned long nrInds;
  
  table = this;

  //cout << " check - dimension " << dimension << " nrBins " << nrBins << "\n";

  //cout << "Ceck ";
  /* for(i = 0; i < dimension; i++)
    cout << indiv[i] << "\t";
    cout << "\n";*/
  
  for(i = 0 ; i < dimension ; i++){
      //check each individual with a fix table value
      int temp = indiv[i];
 
      if(temp < 0 || temp > nrBins - 1){
	cout << " arror in hashtable; attempt to introduce finer order histogram then precizated: " 
	     << temp << "hashtable_histogram.cpp:154\n";
	exit(1);
      }
	
      if (table->histogram[temp] == NULL){
	//cout << " not found " << temp << " in dimension " << i << "\n";
	return 0;
      }
      else
	{
	  //cout << "gasit (" <<indiv[i] << "," << i << ","<< table-> nrIndividuals << "," << table->mean <<")";
	  table = table->histogram[temp];
	  nrInds = table->nrIndividuals;
	  value = table->mean;
	  square = table->meanSqrt;
	  inverse = table->inverseMean;
	}
  }
  
  //  cout << "\n";
  return nrInds;
}

void Hashtable_histogram :: store_individual(unsigned long* indiv, double value)
{
  unsigned long i,j;
  unsigned long output;
  Hashtable_histogram *table,*firsttable;

  //  cout << "store dimension " << dimension << " nrBins "<< nrBins << "\n";
  //for(unsigned long i = 0; i < dimension; i++)
  //cout << indiv[i] << "\t";
  //cout << "\n";
  
  //cout << " before global mean " << firsttable->mean << " nr Indiv " << firsttable->nrIndividuals << " inv Mean " <<
  //  firsttable->inverseMean << "\n";

  firsttable = this;

  firsttable->mean += value;
  firsttable->nrIndividuals++;
  firsttable->inverseMean += 1.0/(double)value;

  table = this;

  for(i = 0 ; i < dimension ; i++)
    {
      int temp = indiv[i];
 
      if(temp < 0 || temp > nrBins - 1){
	cout << " error in hashtable; attempt to introduce finer order histogram then precizated: " << 
	    temp << "hashtable_histogram.cpp:210\n";
	exit(1);
      }

      if (table->histogram[temp] == NULL){
	  table->histogram[temp] = new Hashtable_histogram(dimension,nrBins);
	  table = table->histogram[temp];
	  table->first = 1; 
	  table-> nrIndividuals++;

	  //firsttable-> nrIndividuals++;
	  firsttable->nrElem++;
	  //firsttable->mean += value;
	  //firsttable->nrIndividuals++;
	  //firsttable->inverseMean += 1.0/(double)value;

	  table-> mean += value;
	  table-> meanSqrt += pow(value,2.0);
	  table->inverseMean += 1.0/(double)value;

	  //cout << "nou (" << i << " , val="<< indiv[i] << ", nr=" << table-> nrIndividuals 
	  //   << ", mean=" << table->mean << ") \t";

	  for(j = i+1 ; j < dimension ; j++){
	      int tempJ  = indiv[j];
		
	      if(tempJ < 0 || tempJ > nrBins - 1){
		cout << " arror in hashtable; attempt to introduce finer order histogram then precizated: " << 
		    temp << "hashtable_histogram.cpp:117\n";
		exit(1);
	      }
	    
	      table->histogram[tempJ] = new Hashtable_histogram(dimension,nrBins);
	      table = table->histogram[tempJ];
	      table -> first = 1;
	      table-> nrIndividuals++;
	      table-> mean += value;
	      table->meanSqrt += pow(value,2.0);
	      table->inverseMean += 1.0/(double)value;
	      //  cout << "(" << j << " , val="<<indiv[j] << ", nr="<< table-> nrIndividuals 
	      //   << ", mean=" << table->mean <<") \t";
	    }
	  //cout << "\n";
	  return;
      }else{
	      table = table->histogram[temp];
	      table -> first = 1;
	      table->nrIndividuals++;
	      table -> mean += value; 
	      table->meanSqrt += pow(value,2.0);
	      table->inverseMean += 1.0/(double)value;

	      //add also for the first elelemnt
	      //firsttable -> nrIndividuals++;
	      //firsttable -> mean += value;
	      //firsttable -> inverseMean += 1.0/value;
	      //  cout << "gasit ( val " <<indiv[i] << ", at " << i << ", nr"<< table-> nrIndividuals 
	      //   << ", mean " << table->mean<< ")";
	    }
    }

  //cout << "\n";
}

unsigned long Hashtable_histogram :: delete_individual(unsigned long* indiv){
unsigned long i;
unsigned long output;
Hashtable_histogram *table;
 unsigned long nrInds;

 table = this;

 for(i = 0 ; i < dimension ; i++){
   int temp = indiv[i];

   if(temp < 0 || temp > nrBins - 1){
     cout << " arror in hashtable; attempt to introduce finer order histogram then precizated: " 
	  << temp << "hashtable_histogram.cpp:117\n";
     exit(1);
   }

   if (table->histogram[temp] == NULL)
     return 0;
   else
     {
       table->nrIndividuals --;
       nrInds = table->nrIndividuals;
       table = table->histogram[temp];
     }
 }
 return nrInds;
}

#ifdef FIRST_ITERATOR
long Hashtable_histogram :: iteratorInit(unsigned long* indiv, double& value, double& squared, double& inverseValue){
    //listCurrent = new Hashtable_histogram*[nr_vars];
    //indexCurrent = new unsigned long[nr_vars];

  if(this == NULL)
    return -1;
  
  listCurrent[0] = this;
  
  unsigned long i =0;
  //cout << "indiv=";

  while(i < dimension){
    unsigned long j = 0;
    
    while(j < nrBins){
      if(listCurrent[i]->histogram[j] != NULL){
	indexCurrent[i] = j;
	//cout << "( at " << j << ", nr " << listCurrent[i] ->nrIndividuals << ", val " << listCurrent[i] ->mean << ")";
	indiv[i] = j;
	if(i != dimension) 
	  listCurrent[i+1] = listCurrent[i] -> histogram[j];
	i++;
	break;
      }
      j++;
    }
    
    //if there is nothing in list 
    if(j == nrBins){
      //cout << "notihng in hashtable \n";
      return -1;
    }
  }

  //cout << "indiv  = " << listCurrent[dimension]->nrIndividuals << " mean " << listCurrent[dimension]->mean << "\n";
  value = listCurrent[dimension]->mean;
  squared = listCurrent[dimension]->meanSqrt;
  inverseValue = listCurrent[dimension]->inverseMean;
  //cout << " value " << value << "\n";
  return listCurrent[dimension]->nrIndividuals;
}


void Hashtable_histogram :: iteratorFree()
{
//    if(indexCurrent != NULL)
//	delete[] indexCurrent;
}

long Hashtable_histogram :: iterator(unsigned long* indiv, double& value, double& squared, double& inverseValue)
{
    int i=0;
    //Hashtable_histogram* table = listCurrent[nr_vars-1];
    /*cout << "indiv precedent=";
    for(i = 0; i<dimension; i++){
	indiv[i] = indexCurrent[i];
	cout << indiv[i];
    }
    cout << "\n";
    */
    //if(nrElem == 0 || nrIndividuals == 0)
    // return -1;

    i = dimension-1;

    //inainte daca nu toate dimensiunile posibile au fost vizitate
    //otherwise, inapoi
    while(i>=0){
	//cout << " I am at i = " << i << "\t";
      //go further in histogram
      unsigned long j = 1;
      while(indexCurrent[i] + j < nrBins){
	  //cout << " I am at j = " << j << " + index " << indexCurrent[i] + j << "\n";
	if(listCurrent[i] -> histogram[indexCurrent[i] + j] != NULL){
	  //found other location
	  indiv[i] = indexCurrent[i] + j;
	  indexCurrent[i] = indexCurrent[i] + j;
	  //  cout << "found entrance (at i=" << i << ", j=" << j << ", with " << indiv[i] << ", nr "<< listCurrent[i]->nrIndividuals 
	  //   << ", value " << listCurrent[dimension]->mean << ")\t";
	  if(i != dimension) listCurrent[i+1] = listCurrent[i]->histogram[indexCurrent[i]];
	  i++;
	  goto a12;
	} //else cout << "null j \n;";
	j++;
      }
      //cout << " micsorat " << i << "\n";
      i--;
      if(i<0) {
	//  cout << "trebuie sa te termine\n";
	return -1;
      }
    }
    a12:
    //cout << "la sfirsti i=" <<i <<"\n";

    //for the new values -> inainte
    //daca a gasit un drum nou acesta nu poate fi gol
    //cout << "inaunsigned longe from i=" << i << "\t";
    int gasit = 0;
    while(i < dimension){
      int j = 0;
      gasit = 0;
      //cout << " start counting j at" << indexCurrent[i]+j << " sould be " << listCurrent[i]->histogram[indexCurrent[i] + j] << "\n";
      while(j < nrBins){
	if(listCurrent[i]->histogram[j] != NULL){ 
	  //found another entry
	  gasit = 1;
	  indexCurrent[i] = j;
	  indiv[i] = indexCurrent[i];
	  //cout << "( at i=" << i << ", j=" << j << ", with " << indiv[i] << ",nr " << listCurrent[i]->nrIndividuals 
	  //   << ",val " << listCurrent[dimension]->mean << ") ";
	  if(i != dimension) listCurrent[i+1] = listCurrent[i]->histogram[j];
	  break;
	} 
	j++;
      }
      if(gasit == 0){
	cout << "Greseala de constrictie in hashtable pentru variabla "<<i 
	     << "in Hashtable_histogram.cpp:312\n";
	exit(1);
      }
      i++;
    }

    //  cout << "indiv  = " << listCurrent[dimension]->nrIndividuals << " mean " << listCurrent[dimension]->mean << "\n";
    value = listCurrent[dimension]->mean;
    squared = listCurrent[dimension]->meanSqrt;
    inverseValue = listCurrent[dimension]->inverseMean;
    //cout << " value " << value << "\n";
  //    cout << "\n";
    return listCurrent[dimension]->nrIndividuals;
}
#endif //first iterator

#ifdef SECOND_ITERATOR
long Hashtable_histogram :: iteratorSecondInit(unsigned long* indiv, double& value, double& squared, double& inverseValue){
    //listCurrent = new Hashtable_histogram*[nr_vars];
    //indexSecondCurrent = new unsigned long[nr_vars];

  if(this == NULL)
    return -1;
  
  listSecondCurrent[0] = this;
  
  unsigned long i =0;
  //cout << "indiv=";

  while(i < dimension){
    unsigned long j = 0;
    
    while(j < nrBins){
      if(listSecondCurrent[i]->histogram[j] != NULL){
	indexSecondCurrent[i] = j;
	//cout << "( at " << j << ", nr " << listSecondCurrent[i] ->nrIndividuals << ", val " << listSecondCurrent[i] ->mean << ")";
	indiv[i] = j;
	if(i != dimension) 
	  listSecondCurrent[i+1] = listSecondCurrent[i] -> histogram[j];
	i++;
	break;
      }
      j++;
    }
    
    //if there is nothing in list 
    if(j == nrBins){
      //cout << "notihng in hashtable \n";
      return -1;
    }
  }

  //cout << "indiv  = " << listSecondCurrent[dimension]->nrIndividuals << " mean " << listSecondCurrent[dimension]->mean << "\n";
  value = listSecondCurrent[dimension]->mean;
  squared = listSecondCurrent[dimension]->meanSqrt;
  inverseValue = listSecondCurrent[dimension]->inverseMean;
  // cout << " value " << value << "\n";
  return listSecondCurrent[dimension]->nrIndividuals;
}


long Hashtable_histogram :: iteratorSecond(unsigned long* indiv, double& value, double& squared, double& inverseValue)
{
    int i=0;
    //Hashtable_histogram* table = listSecondCurrent[nr_vars-1];
    /*cout << "indiv precedent=";
    for(i = 0; i<dimension; i++){
	indiv[i] = indexSecondCurrent[i];
	cout << indiv[i];
    }
    cout << "\n";
    */
    //if(nrElem == 0 || nrIndividuals == 0)
    // return -1;

    i = dimension-1;

    //inainte daca nu toate dimensiunile posibile au fost vizitate
    //otherwise, inapoi
    while(i>=0){
	//cout << " I am at i = " << i << "\t";
      //go further in histogram
      unsigned long j = 1;
      while(indexSecondCurrent[i] + j < nrBins){
	  //cout << " I am at j = " << j << " + index " << indexSecondCurrent[i] + j << "\n";
	if(listSecondCurrent[i] -> histogram[indexSecondCurrent[i] + j] != NULL){
	  //found other location
	  indiv[i] = indexSecondCurrent[i] + j;
	  indexSecondCurrent[i] = indexSecondCurrent[i] + j;
	  //  cout << "found entrance (at i=" << i << ", j=" << j << ", with " << indiv[i] << ", nr "<< listSecondCurrent[i]->nrIndividuals 
	  //   << ", value " << listSecondCurrent[dimension]->mean << ")\t";
	  if(i != dimension) listSecondCurrent[i+1] = listSecondCurrent[i]->histogram[indexSecondCurrent[i]];
	  i++;
	  goto a12;
	} //else cout << "null j \n;";
	j++;
      }
      //cout << " micsorat " << i << "\n";
      i--;
      if(i<0) {
	//  cout << "trebuie sa te termine\n";
	return -1;
      }
    }
    a12:
    //cout << "la sfirsti i=" <<i <<"\n";

    //for the new values -> inainte
    //daca a gasit un drum nou acesta nu poate fi gol
    //cout << "inaunsigned longe from i=" << i << "\t";
    int gasit = 0;
    while(i < dimension){
      int j = 0;
      gasit = 0;
      //cout << " start counting j at" << indexSecondCurrent[i]+j << " sould be " << listSecondCurrent[i]->histogram[indexSecondCurrent[i] + j] << "\n";
      while(j < nrBins){
	if(listSecondCurrent[i]->histogram[j] != NULL){ 
	  //found another entry
	  gasit = 1;
	  indexSecondCurrent[i] = j;
	  indiv[i] = indexSecondCurrent[i];
	  //cout << "( at i=" << i << ", j=" << j << ", with " << indiv[i] << ",nr " << listSecondCurrent[i]->nrIndividuals 
	  //   << ",val " << listSecondCurrent[dimension]->mean << ") ";
	  if(i != dimension) listSecondCurrent[i+1] = listSecondCurrent[i]->histogram[j];
	  break;
	} 
	j++;
      }
      if(gasit == 0){
	cout << "Greseala de constrictie in hashtable pentru variabla "<<i 
	     << "in Hashtable_histogram.cpp:312\n";
	exit(1);
      }
      i++;
    }

    //cout << "indiv  = " << listSecondCurrent[dimension]->nrIndividuals << " mean " << listSecondCurrent[dimension]->mean << "\n";
    value = listSecondCurrent[dimension]->mean;
    squared = listSecondCurrent[dimension]->meanSqrt;
    inverseValue = listSecondCurrent[dimension]->inverseMean;
    //cout << " value " << value << "\n";
  //    cout << "\n";
    return listSecondCurrent[dimension]->nrIndividuals;
}

#endif //second iterator

//print each histohram in hashtable in a different file
void Hashtable_histogram::printEachHistogram(Hashtable_histogram** myHashtables, int records){
  double iteratorNumber, iteratorTotalNumber = 0, iteratorTotalValue = 0, iteratorTotalInverse = 0, iteratorSqrtValue = 0, iteratorValue = 0, iteratorInverseValue = 0;
  unsigned long* iteratorTemp = new unsigned long[records]; 
  ofstream myFile2Dstd;
  char nameFile2Dstd[30] = "Hist2DcDim"; 
  for(int j = 0; j < records; j++){
    nameFile2Dstd[10] = (int)(j / 10) + 48;
    nameFile2Dstd[11] = j - 10 * ((int)(j / 10)) + 48;
    nameFile2Dstd[12] = '\0';
    myFile2Dstd.open(nameFile2Dstd);
  
    iteratorNumber = myHashtables[j]->iteratorInit(iteratorTemp,iteratorValue,iteratorSqrtValue,iteratorInverseValue);
    iteratorTotalNumber = myHashtables[j]->nrIndividuals;
    iteratorTotalValue = myHashtables[j]->mean;
    iteratorTotalInverse = myHashtables[j] -> inverseMean;
 
    while(iteratorNumber != -1){
      for(int i = 0; i < dimension; i++)
	myFile2Dstd << iteratorTemp[i] << "\t";
      myFile2Dstd << iteratorValue/iteratorNumber << "\t" << iteratorSqrtValue/iteratorNumber - pow(iteratorValue/iteratorNumber,2.0) << "\t"
		  << iteratorNumber/iteratorTotalNumber << "\n";

      iteratorNumber = myHashtables[j] -> iterator(iteratorTemp,iteratorValue,iteratorSqrtValue,iteratorInverseValue);
      iteratorTotalNumber = myHashtables[j]->nrIndividuals;
      iteratorTotalValue = myHashtables[j]->mean;
      iteratorTotalInverse = myHashtables[j]->inverseMean;
    }

    myFile2Dstd.flush();
    myFile2Dstd.close();    
  }
  delete[] iteratorTemp;
}

//bubble sort
int* Hashtable_histogram::bubble(unsigned long** iteratorTemp, int* orderHistogram, int records, int dimension){
  
    for(int iBin = 0; iBin < records; iBin++)
      for(int jBin = 0; jBin < records; jBin++){
	if(iBin != jBin)
	  for(int j = 0; j < dimension; j++){
	    //cout << " compare iBin = " << iBin << " jBin " << jBin << " j " << j << " (" << orderHistogram[iBin] 
	    //    << "," <<orderHistogram[jBin] << " ," << iteratorTemp[orderHistogram[iBin]][j] << "," << iteratorTemp[orderHistogram[jBin]][j] << ")\n";
	    if(orderHistogram[jBin] != -1 && orderHistogram[iBin] != -1 && iteratorTemp[orderHistogram[iBin]][j] > iteratorTemp[orderHistogram[jBin]][j] && iBin < jBin){
	      //cout << " exchange \n";
	      long temp = orderHistogram[iBin];
	      orderHistogram[iBin] = orderHistogram[jBin];
	      orderHistogram[jBin] = temp;	  
	      break;
	    } else 
	      if(orderHistogram[jBin] != -1 && orderHistogram[iBin] != -1 && iteratorTemp[orderHistogram[iBin]][j] == iteratorTemp[orderHistogram[jBin]][j])
		continue;
	      else break;
	  }     
      }
    
    return orderHistogram;
}

//aditional printing !!! ATENTIE - change with type of function
void Hashtable_histogram::printAdditionalHistograms(ofstream* myFile2Dstd, ofstream* myFile2DContur, ofstream* myFile2DDiag, unsigned long* iteratorTemp, double* tempValue, double* tempStd, int newLine2Dc){
  //void Hashtable_histogram::printAdditionalHistograms(ofstream* myFile2Dstd, ofstream* myFile2DContur, ofstream* myFile3D, unsigned long* iteratorTemp, double* tempValue, double* tempStd){
  for(int k = 0; k < dimension; k++){
    bool stop = true; 
    for(int j = 0; j < dimension; j++)
      if(k != j && iteratorTemp[j] != NrBins/2) stop = false;
    if(stop == true){
      for(int i = 0; i < dimension; i++)
	myFile2Dstd[k] << iteratorTemp[i] << "\t";
      for(int i = 0; i < 9; i++)
	myFile2Dstd[k] << tempValue[i] << "\t" << tempStd[i] << "\t";
      myFile2Dstd[k] << "\n";
    }
  }

  //bool memStop = true;
  for(int k = 0 ; k < dimension; k++)
    for(int j = k+1; j < dimension; j++){
	bool stop = true;
	for(int j1 = 0; j1 < dimension; j1++)
	  if(k != j1 && j != j1 && iteratorTemp[j1] != NrBins/2) stop = false;
	if(stop == true){
	  long tempValueIndex = k * (dimension -1) - k * (k-1)/2 + j - k - 1;
	if(newLine2Dc != iteratorTemp[0]) myFile2DContur[tempValueIndex] << "\n";
	  //cout << "write contru in file " << tempValueIndex << "\n"; 
	for(int i = 0; i < dimension; i++)
	  myFile2DContur[tempValueIndex] << iteratorTemp[i] << "\t";
	for(int i = 0; i < 9; i++)
	  myFile2DContur[tempValueIndex] << tempValue[i] << "\t" << tempStd[i] << "\t";
	myFile2DContur[tempValueIndex] << "\n";
      }
    }

//diag file 
  for(int k = 0 ; k < dimension; k++)
    for(int j = k+1; j < dimension; j++){
	bool stop = true;
	for(int j1 = 0; j1 < dimension; j1++)
	    for(int j2 = j1+1; j2 < dimension; j2++)
		if(iteratorTemp[j1] != iteratorTemp[j2])  stop = false;
	if(stop == true){
	  long tempValueIndex = k * (dimension -1) - k * (k-1)/2 + j - k - 1;
	  //cout << "write contru in file " << tempValueIndex << "\n"; 
	for(int i = 0; i < dimension; i++)
	  myFile2DDiag[tempValueIndex] << iteratorTemp[i] << "\t";
	for(int i = 0; i < 9; i++)
	  myFile2DDiag[tempValueIndex] << tempValue[i] << "\t" << tempStd[i] << "\t";
	myFile2DDiag[tempValueIndex] << "\n";
      }
    }
  
  //3D histogram
  /*if(dimension > 3){
    long tempValueIndex = 0;
    for(int k = 0 ; k < dimension; k++)
      for(int j = k+1; j < dimension; j++)
	for(int l = j+1; l< dimension; l++){
	  bool stop = true;
	    tempValueIndex++;
	  for(int j1 = 0; j1 < dimension; j1++)
	    if(k != j1 && j != j1 && l != j1 && iteratorTemp[j1] != NrBins/2) stop = false;
	  if(stop == true){
	    //cout << "write contru in file " << tempValueIndex << "\n"; 
	    for(int i = 0; i < dimension; i++)
	    myFile3D[tempValueIndex] << iteratorTemp[i] << "\t";
	    for(int i = 0; i < 6; i++)
	      myFile3D[tempValueIndex] << tempValue[i] << "\t" << tempStd[i] << "\t";
	    myFile3D[tempValueIndex] << "\n";
	  }
	}
	}*/
}

// print the mean histogram of all histogrmas
//records - # of separat records that are kept
void Hashtable_histogram::printStd(ofstream& myFile, Hashtable_histogram** myHashtables, int records, double& centerMass, double& populationVariance){
    
  printEachHistogram(myHashtables, records);

  cout << " print dimension " << dimension << " nrBins " << nrBins << " records" << records << "\n";
  double binStd = 0;
  int binNrIndiv = 0;
  double tempValue;
  double meanStd, meanValue;
  double* tempStd = new double[9];
  double* tempMean = new double[9];
  for(int i = 0; i < 9; i++){
    tempStd[i] = 0;
    tempMean[i] = 0;
  }
  //cout << " before " << dimension << " , " << nrBins << "\n";
  int j, iBin;

  double number = pow((double)nrBins,(double)dimension);
//  double* tableStd = new double[(unsigned long)number];
  //double* tableValue = new double[(unsigned long)number];
  //for(unsigned long i = 0; i < number; i++){
      //tableStd[i] = 0;
  //tableValue[i] = 0;
  //}

  //structures for the iterator
#ifdef FIRST_ITERATOR
  //cout << "first elements in iterator " << records << " \n";
  unsigned long** iteratorTemp = new unsigned long*[records];
  double* iteratorValue = new double[records];
  double* iteratorSqrtValue = new double[records];
  double* iteratorInverseValue = new double[records];
  unsigned long* iteratorNumber = new unsigned long[records];
  double* iteratorTotalNumber = new double[records];
  double* iteratorTotalValue = new double[records];
  double* iteratorTotalInverse = new double[records];
  for(iBin = 0; iBin < records; iBin++){
    iteratorTemp[iBin] = new unsigned long[dimension];
    iteratorValue[iBin] = 0;
    iteratorSqrtValue[iBin] = 0;
    iteratorInverseValue[iBin] = 0;
    iteratorNumber[iBin] = myHashtables[iBin]->iteratorInit(iteratorTemp[iBin],iteratorValue[iBin],iteratorSqrtValue[iBin],iteratorInverseValue[iBin]);
    iteratorTotalNumber[iBin] = myHashtables[iBin]->nrIndividuals;
    iteratorTotalValue[iBin] = myHashtables[iBin]->mean;
    iteratorTotalInverse[iBin] = myHashtables[iBin] -> inverseMean;
    /*cout << " indiv = " << iBin << " is (";
    for(int k = 0; k < dimension; k++)
     cout << iteratorTemp[iBin][k] << ",";
    cout << ") = (val = (" << iteratorValue[iBin] << "," << iteratorTotalValue[iBin] << "), numb = (" 
	 << iteratorNumber[iBin] << "," << myHashtables[iBin]->nrElem << " = " << iteratorTotalNumber[iBin] 
	 << ") inv Val = (" << iteratorInverseValue[iBin] << "," << iteratorTotalInverse[iBin]
	 << ")) \n";*/
  }
#endif //first iterator

  /*  ofstream myFile("histogram_hashtable.file");		
  if(!myFile.is_open()){
      cout << "Error opening histogram file \n";
      return;
      }*/

  //aditional histograms
  ofstream* myFile2Dstd = new ofstream[dimension];
  ofstream* myFile2DContur = new ofstream[dimension * (dimension - 1)/2];
  ofstream* myFile2DDiag = new ofstream[dimension * (dimension - 1)/2];
  //if(dimension > 3)
  //ofstream* myFile3D = new ofstream[dimension * (dimension - 1) * (dimension - 2)/6];

#ifdef IsRecombination
  cout << "write in the histogram for recombination \n";
  //#ifdef IsRecombination
  char nameFile2Dstd[30] = "outRecombHist2DsDim"; 
  for(j = 0; j < dimension; j++){
    nameFile2Dstd[19] = (int)(j / 10) + 48;
    nameFile2Dstd[20] = j - 10 * ((int)(j / 10)) + 48;
    nameFile2Dstd[21] = '\0';
    myFile2Dstd[j].open(nameFile2Dstd);
  }

  nameFile2Dstd[15] = 'c'; 
  for(int k = 0; k < dimension -1; k++)
    for(j = k+1; j < dimension; j++){
      int tempValueIndex = k * (dimension -1) - k * (k-1)/2 + j - k - 1;
      //cout << "dimension " << dimension << " k " << k << " j " << j << " tempValueIndex " << tempValueIndex << "\n";
      nameFile2Dstd[19] = (int) (tempValueIndex / 10) + 48;
      nameFile2Dstd[20] = tempValueIndex - 10 * ((int) (tempValueIndex / 10)) + 48;
      nameFile2Dstd[21] = '\0';
      myFile2DContur[tempValueIndex].open(nameFile2Dstd);
    }

  nameFile2Dstd[15] = 'd'; 
  for(int k = 0; k < dimension -1; k++)
    for(j = k+1; j < dimension; j++){
      int tempValueIndex = k * (dimension -1) - k * (k-1)/2 + j - k - 1;
      //cout << "dimension " << dimension << " k " << k << " j " << j << " tempValueIndex " << tempValueIndex << "\n";
      nameFile2Dstd[19] = (int) (tempValueIndex / 10) + 48;
      nameFile2Dstd[20] = tempValueIndex - 10 * ((int) (tempValueIndex / 10)) + 48;
      nameFile2Dstd[21] = '\0';
      myFile2DDiag[tempValueIndex].open(nameFile2Dstd);
    }

  /*nameFile2Dstd[13] = '3'; 
    int tempValueIndex = 0;
  for(int k = 0; k < dimension -1; k++)
    for(j = k+1; j < dimension; j++)
    for(int l = j+1; l < dimension; l++){
    tempValueIndex++;
    //cout << "dimension " << dimension << " k " << k << " j " << j << " tempValueIndex " << tempValueIndex << "\n";
    nameFile2Dstd[19] = tempValueIndex + 48;
    nameFile2Dstd[20] = '\0';
    myFile3D[tempValueIndex].open(nameFile2Dstd);
    }
  */
#else
  char nameFile2Dstd[30] = "outMutHist2DsDim"; 
  for(j = 0; j < dimension; j++){
    nameFile2Dstd[16] = (int)(j / 10) + 48;
    nameFile2Dstd[17] = j - 10 * (int)(j / 10) + 48;
    nameFile2Dstd[18] = '\0';
    myFile2Dstd[j].open(nameFile2Dstd);
  }
  
  nameFile2Dstd[12] = 'c'; 
  for(int k = 0; k < dimension; k++)
    for(j = k+1; j < dimension; j++){
      int tempValueIndex = k * (dimension -1) - k * (k-1)/2 + j - k - 1;
      //cout << "dimension " << dimension << " k " << k << " j " << j << " tempValueIndex " << tempValueIndex << "\n";
      nameFile2Dstd[16] = (int)(tempValueIndex / 10) + 48;
      nameFile2Dstd[17] = tempValueIndex - 10 * (int)(tempValueIndex / 10) + 48;
      nameFile2Dstd[18] = '\0';
      myFile2DContur[tempValueIndex].open(nameFile2Dstd);
    }

  nameFile2Dstd[12] = 'd'; 
  for(int k = 0; k < dimension; k++)
    for(j = k+1; j < dimension; j++){
      int tempValueIndex = k * (dimension -1) - k * (k-1)/2 + j - k - 1;
      //cout << "dimension " << dimension << " k " << k << " j " << j << " tempValueIndex " << tempValueIndex << "\n";
      nameFile2Dstd[16] = (int)(tempValueIndex / 10) + 48;
      nameFile2Dstd[17] = tempValueIndex - 10 * (int)(tempValueIndex / 10) + 48;
      nameFile2Dstd[18] = '\0';
      myFile2DDiag[tempValueIndex].open(nameFile2Dstd);
    }

  /*nameFile2Dstd[11] = '3'; 
    int tempValueIndex = 0;
  for(int k = 0; k < dimension -1; k++)
    for(j = k+1; j < dimension; j++)
    for(int l = j+1; l < dimension; l++){
    tempValueIndex++;
    //cout << "dimension " << dimension << " k " << k << " j " << j << " tempValueIndex " << tempValueIndex << "\n";
    nameFile2Dstd[16] = tempValueIndex + 48;
    nameFile2Dstd[17] = '\0';
    myFile3D[tempValueIndex].open(nameFile2Dstd);
    }
  */
#endif 
  //     myFile2DContur[0].open("output_recomb_histogram_2D_contru_dim1");
  //  myFile2DContur[1].open("output_recomb_histogram_2D_contru_dim2");
  //  myFile2DContur[2].open("output_recomb_histogram_2D_contru_dim3");
  

  //#define WHOLE_HISTOGRAM

#ifdef WHOLE_HISTOGRAM

  double binHight[records];

  unsigned long point[dimension];
  for(j = 0; j < dimension; j++)
    point[j] = 0;

  for(unsigned long i = 0; i < (unsigned long)number; i++){    
    //the bin value
    //binHight = 0;
    /*for(unsigned long j = 0 ; j < dimension; j++)
      cout << point[j] << "\t";
      cout << "\n";*/

    binStd = 0;
    binNrIndiv = 0;

    tempValue = 0;

#ifndef FIRST_ITERATOR
    for(iBin = 0; iBin < records; iBin++){ 
      int temp = myHashtables[iBin]->check_individual(point,tempValue,tempInverseValue);
      if(temp != -1){
	binNrIndiv += temp;

	binHight[iBin]= tempValue;

	binStd += tempValue;
      }
    }// end for iBin
#else // first iterator
    for(iBin = 0; iBin < records; iBin++){
      int count = 0;
      for(j = 0; j < dimension; j++)
	if(point[j] != iteratorTemp[iBin][j])
	  count = 1;
      if(count == 0){
	  binNrIndiv += iteratorNumber[iBin];
	  
	  binHight[iBin]= iteratorValue[iBin];
	  
	  binStd += iteratorValue[iBin];

	  //move to the next level
	  myHashtables[iBin]->iterator(iteratorTemp[iBin],value);
      }
    }
  
#endif //first iterator
    //  cout << " 1 dimension " << dimension << " nrBins " << nrBins << "\n";

    //compute the mean
    if(binNrIndiv != 0){
      tempValue = binStd / binNrIndiv;
         binStd = 0;
      
      for(iBin = 0; iBin < records; iBin++)
	binStd += pow(binHight[iBin] - tempValue,2.0);
      
      binStd = sqrt(binStd/(double)records);
    }
    //cout << " 2 dimension " << dimension << " nrBins " << nrBins << "\n";
 
    //write the point in the file
    for(j = 0; j < dimension; j++)
      myFile << point[j] << "\t";
    myFile << tempValue << "\t" << binStd << "\n";
    //cout << " 3 dimension " << dimension << " nrBins " << nrBins << "\n";

    //increment the point - start from the last position and go towards the first one
    j = dimension - 1;
    while(j >= 0){
      point[j]++;
      if(point[j] < nrBins){
	myFile << "\n";
	j = -1;
      }
      else {
	point[j] = 0;
	j--;
      }
    }//end while j
    //cout << " 4 dimension " << dimension << " nrBins " << nrBins << "\n";

    tableStd[i] = binStd;
    meanStd += binStd * tempValue;
    tableValue[i] = tempValue;
    meanValue += tempValue;
  } // end for i 

  centerMass = meanStd/meanValue;
  
  populationVariance  = 0;
  populationMean  = 0;
  for(unsigned long long i = 0; i < number; i++){    
  populationVariance += pow(tableStd[i] - meanStd,2.0) * tableValue[i];
   populationMean += (tableStd[i] - meanStd) * tableValue[i];
  }
  populationVariance = sqrt(populationVariance/pow((double)nrBins,(double)dimension));

  delete[] tableStd;
  delete[] tableValue;
  
#else //whole histogram

  //table with the order of individuals
  int orderHistogram[records];
  int temp, stop;
  //order 
  for(iBin = 0; iBin < records; iBin++){
    orderHistogram[iBin]= iBin;
  }

  bubble(iteratorTemp,orderHistogram,records,dimension);

  //cout << "orderHistogram (";
  //for(iBin = 0; iBin < records; iBin++)
  //  cout << orderHistogram[iBin] << ",";
  //cout << ")\n";

  //a data for number of records
  unsigned long long number_record = 0;
  //std given different runs
  double binHight = 0;
  double binHightInverse = 0;
  double binHightSqrt = 0;
  //  double tempStd = 0;

  //std for all run
  double sumBinHight = 0;
  double sumBinHightInverse = 0;
  double sumBinHightSqrt = 0;
  double sumNrIndiv = 0;

//means for each separate records
  double meanHight = 0;
  double meanInverse = 0;
  double stdHight = 0;
  double stdInverse = 0;

//frequncies
  double freqHight =0;
  double freqStd = 0;
  double freqMean = 0;
  double freqMeanStd = 0;
  double freqInverse = 0;
  double freqInverseStd = 0;

  populationVariance  = 0;
//  populationMean  = 0;
  int prevK = -1;
  //until the last element from list were visited
  while((temp = orderHistogram[0]) != -1){
    
    number_record++;
    //cout << " temp " << temp << " number_record " << number_record << "\n";
 
    //for the first elelement
    if(iteratorNumber[temp] != -1){
      binNrIndiv = 1 ;
      sumNrIndiv = (double)iteratorNumber[temp]; //iteratorNumber[temp];
    }
    else break;

    if(iteratorValue[temp] > 0){
      binHightInverse = (double)iteratorNumber[temp]/(double)iteratorValue[temp];
      binHight = (double)iteratorValue[temp]/(double)iteratorNumber[temp];
      binHightSqrt = pow((double)iteratorValue[temp]/(double)iteratorNumber[temp],2.0);

      sumBinHight = (double)iteratorValue[temp];
      sumBinHightInverse = (double)iteratorInverseValue[temp];
      sumBinHightSqrt = (double)iteratorSqrtValue[temp];
      
      meanHight = iteratorValue[temp]/iteratorNumber[temp];
      meanInverse = iteratorNumber[temp]/iteratorInverseValue[temp];
      stdHight = pow(iteratorValue[temp]/iteratorNumber[temp],2.0);
      stdInverse = pow(iteratorNumber[temp]/iteratorInverseValue[temp],2.0);

      freqHight = iteratorNumber[temp]/iteratorTotalNumber[temp];
      freqMean = iteratorValue[temp] * iteratorTotalNumber[temp]/ (iteratorNumber[temp] * iteratorTotalValue[temp]);
      freqStd = pow(iteratorNumber[temp]/iteratorTotalNumber[temp],2.0);
      freqMeanStd = 
	 pow(iteratorValue[temp] * iteratorTotalNumber[temp]/(iteratorNumber[temp] * iteratorTotalValue[temp]),2.0);
      freqInverse = 
	 (iteratorNumber[temp]/iteratorTotalNumber[temp]) * (iteratorTotalInverse[temp]/iteratorInverseValue[temp]);
      freqInverseStd = 
  pow((iteratorNumber[temp]/iteratorTotalNumber[temp])*(iteratorTotalInverse[temp]/iteratorInverseValue[temp]),2.0);
    }
    else 
      if(iteratorValue[temp] == 0){
	cout << " 0 value for ";
	for(int k = 0; k < dimension; k++)
	  cout << iteratorTemp[temp][k] << "\t";
	cout << " \n";
	exit(1);
      }
      else break;
    
    //binStd = iteratorValue[temp];
    
    //pop-up the first elements from the orderHistogram
    j = 1;
    //look for equal elements
    stop = 0;
    while(stop == 0 && j < records){
      //cout << " j = " << j << " stop = " << stop << "\n"; 
      stop = 0;
      if(orderHistogram[j] == -1) {
	stop = 1;
	//cout << " notihng to find \n";
      }
      else{
	for(iBin = 0; iBin < dimension; iBin++){
	  //cout << " iBin " << iBin << " orderHistogram[j] " << orderHistogram[j] << " temp " 
	  //   << temp << " iteratorTemp[temp][iBin] " << iteratorTemp[temp][iBin] << " iteratorTemp[orderHistogram[j]][iBin] " 
	  //   << iteratorTemp[orderHistogram[j]][iBin] << "\n"; 
	  if(iteratorTemp[temp][iBin] != iteratorTemp[orderHistogram[j]][iBin]) {
	    stop = iBin + 1;
	    //cout << " stop \n";
	    if(iteratorTemp[temp][iBin] > iteratorTemp[orderHistogram[j]][iBin]){
	      cout << "WRONG order in the print hashtables for ";
	      cout << " iBin " << iBin << " orderHistogram[j] " << orderHistogram[j] << " temp " 
		   << temp << " iteratorTemp[temp][iBin] " << iteratorTemp[temp][iBin] << 
		  " iteratorTemp[orderHistogram[j]][iBin] " << iteratorTemp[orderHistogram[j]][iBin] << "\n"; 
	      for(int k = 0; k < dimension; k++)
		cout <<  iteratorTemp[temp][k] << "\t";
	      cout << "\n";
	      for(int k = 0; k < dimension; k++)
		cout <<  iteratorTemp[orderHistogram[j]][k] << "\t";
	      cout << "\n";
	      
	      exit(1);
	    }
	    break;
	  }
	}
      }
      if(stop == 0){
	//add to the other statistical data
	//cout << " there are equal components \n";
	//if(iteratorValue[orderHistogram[j]] == 0){
	//  cout << " 0 value when (" << j  << "," << orderHistogram[j] << "  for ";
	//  for(int k = 0; k < dimension; k++)
	//	    cout << iteratorTemp[temp][k] << "\t";
	//	  cout << " \n";
	//	  exit(1);
	//	}

	//number of individuals counted
	binNrIndiv++; 
	sumNrIndiv += iteratorNumber[orderHistogram[j]];
	
	//mean given
	binHight += (double)iteratorValue[orderHistogram[j]]/(double)iteratorNumber[orderHistogram[j]];
	binHightSqrt += pow((double)iteratorValue[orderHistogram[j]]/(double)iteratorNumber[orderHistogram[j]],2.0);
	binHightInverse += (double)iteratorNumber[orderHistogram[j]]/(double)iteratorValue[orderHistogram[j]];

	sumBinHight += (double)iteratorValue[orderHistogram[j]];
	sumBinHightInverse += (double)iteratorInverseValue[orderHistogram[j]];
	sumBinHightSqrt += (double)iteratorSqrtValue[orderHistogram[j]];

	meanHight += iteratorValue[temp]/iteratorNumber[temp];
	meanInverse += iteratorNumber[temp]/iteratorInverseValue[temp];

	stdHight += pow(iteratorValue[temp]/iteratorNumber[temp], 2.0);
	stdInverse += pow(iteratorNumber[temp]/iteratorInverseValue[temp], 2.0);

	freqHight += iteratorNumber[temp]/iteratorTotalNumber[temp];
	freqMean += iteratorValue[temp]/iteratorNumber[temp] * iteratorTotalNumber[temp]/iteratorTotalValue[temp];

	freqStd += pow(iteratorNumber[temp]/iteratorTotalNumber[temp],2.0); 
	freqMeanStd += 
	    pow(iteratorValue[temp]/iteratorNumber[temp] * iteratorTotalNumber[temp]/iteratorTotalValue[temp],2.0);

	freqInverse += 
	    iteratorNumber[temp]/iteratorInverseValue[temp] * iteratorTotalInverse[temp]/iteratorTotalNumber[temp];
	freqInverseStd += 
	   pow(iteratorNumber[temp]/iteratorInverseValue[temp]*iteratorTotalInverse[temp]/iteratorTotalNumber[temp],2.0);

	//cout << " value (";
	//for(int k = 0; k < dimension; k++)
	//cout << iteratorTemp[orderHistogram[j]][k] << ",";
	//cout << ") = (" << iteratorValue[orderHistogram[j]] << "," << iteratorNumber[orderHistogram[j]] << ")\n";
	j++;
      }
    }
    
    //compute the mean
    //if(binStd != 0){
    //cout << "first aver val=" << binHight << " nr = (" << records << "," << binNrIndiv << ") sqrt val = " 
    // << binHightSqrt << "\n"; 
    tempMean[0] = binHight/(double)binNrIndiv;
    tempStd[0] = binHightSqrt/(double)binNrIndiv - tempMean[0] * tempMean[0];
    //tempStd[0] = binHightSqrt/(double)records - tempMean[0] * tempMean[0];
    //if(abs(tempStd[0]) < 0.01) tempStd[0] = abs(tempStd[0]);

    tempMean[1] = sumBinHight/(double)sumNrIndiv;
    tempStd[1] = sumBinHightSqrt/(double)sumNrIndiv - pow(tempMean[1],2.0);
    //if(abs(tempStd[1]) < 0.01) tempStd[1] = abs(tempStd[1]);
    
    //tempMean[2] = (double)binNrIndiv/(double)binHightInverse;
    tempMean[2] = (double)binNrIndiv/(double)binHightInverse;
    tempStd[2] = binHight/binHightInverse - pow(tempMean[2],2.0);
    //if(abs(tempStd[2]) < 0.01) tempStd[2] = abs(tempStd[2]);
    
    tempMean[3] = (double)sumNrIndiv/(double)sumBinHightInverse;
    tempStd[3] = sumBinHight/(double)sumBinHightInverse - pow(tempMean[3],2.0);
    //if(abs(tempStd[3]) < 0.01) tempStd[3] = abs(tempStd[3]);
      //} 

    tempMean[4] = meanHight/records;
    tempStd[4] = stdHight/records - pow(tempMean[4],2.0);

    tempMean[5] = meanInverse/records;
    tempStd[5] = stdInverse/records - pow(tempMean[5],2.0);

    //cout << " freq hight =" << freqHight << " records = (" << records << "," << binNrIndiv << " freqStd = " 
    // << freqStd << "\n"; 
    tempMean[6] = freqHight/records;
    tempStd[6] = freqStd/records - pow(tempMean[6],2.0);

    tempMean[7] = freqMean/records;
    tempStd[7] = freqMeanStd/records - pow(tempMean[7],2.0);

    tempMean[8] = freqInverse/records;
    tempStd[8] = freqInverseStd/records - pow(tempMean[8],2.0);
    //write in a file the point found
    //cout << "write in file ";
    for(int k = 0; k < dimension; k++){
      myFile << iteratorTemp[temp][k] << "\t";
      //cout << iteratorTemp[temp][k] << "\t";
    }
    for(int indexMean = 0; indexMean < 9;indexMean++){
      myFile << tempMean[indexMean] << "\t" << tempStd[indexMean] << "\t";
      //cout << tempMean[indexMean] << "\t" << tempStd[indexMean] << "\t";
    }
    //myFile << binNrIndiv << "\t" << binHight << "\t" << binHightSqrt << "\t"<< sumNrIndiv << "\t" << sumBinHight << "\t"<< sumBinHightSqrt;
    myFile << "\n"; // sumNrIndiv << "\t" << binNrIndiv << "\n";
    //cout << "\n";
    // << binNrIndiv << "\t"<< binStd << "\t" << binHight << "\n";
    //if(stop == 1 ) myFile << "\n";
    //cout << tempValue << "\t" << binStd << "\n";
    //printAdditionalHistograms(myFile2Dstd,myFile2DContur,myFile3D,iteratorTemp[temp], tempMean, tempStd);
    
    //cout << "iteratorTemp[0] =" << iteratorTemp[temp][0] << ", prevK = " << prevK << "\n";
    printAdditionalHistograms(myFile2Dstd,myFile2DContur,myFile2DDiag,iteratorTemp[temp], tempMean, tempStd, prevK);

    populationVariance += tempStd[0] * tempStd[0] * tempMean[0];
    meanStd += tempStd[0] * tempMean[0];
    meanValue += tempMean[0];
    
    prevK = iteratorTemp[temp][0];
    //update the value
    //cout << "refill the orderHistogram " << j  << "\n"; 
    for(iBin = 0; iBin < j; iBin++){
      int index = orderHistogram[iBin]; 
      if(index != -1){
	//cout << "the empty iterator " << index << " dimension" << dimension << "\n";
	iteratorNumber[index] = myHashtables[index] -> iterator(iteratorTemp[index],iteratorValue[index],iteratorSqrtValue[index],iteratorInverseValue[index]);
	iteratorTotalNumber[index] = myHashtables[index]->nrIndividuals;
	iteratorTotalValue[index] = myHashtables[index]->mean;
	iteratorTotalInverse[index] = myHashtables[index]->inverseMean;
/*	cout << " indiv " << index << " is = (";
        for(int k = 0; k < dimension; k++)
	    cout << iteratorTemp[index][k] << ",";
	cout << ") = ( val = (" << iteratorValue[index] << "," << iteratorTotalValue[index] << "), numb = (" <<
	    iteratorNumber[index] << "," << myHashtables[index] -> nrElem << " = " << iteratorTotalNumber[index] <<
	    ") inv Val = (" << iteratorInverseValue[index] << ","<< iteratorTotalInverse[index]<< ")) \n";*/
      // for one hashtable we have finish the rearch
      if(iteratorNumber[index] == -1){
	//cout << "finish the records for chain " << index << "\n";
	for(int i = j-1; i < records - 1; i++)
	  orderHistogram[i] = orderHistogram[i+1];
	orderHistogram[records - 1] = -1;
      } /*else {
	//cout << "from that position search another one where to insert it \n";
	//it is a already sorted list, where we add an individual greater than the current one
	int i = j-1; 
	stop = 0;
	while(i < records && stop == 0){
	  if(orderHistogram[i] == -1) stop = 1;
	  else {
	    for(int k = 0; k < dimension; k++)
	      if(iteratorTemp[index][k] == iteratorTemp[orderHistogram[i]][k])
		continue;
	      else{ 
		if(iteratorTemp[index][k] < iteratorTemp[orderHistogram[i]][k] || iteratorTemp[orderHistogram[i]][k] == -1)
		  stop = 1;
		break;
	      }
	  }
	  if(stop == 0) i++;
	}
	//cout << "move all individual for"  << j -1 << " until " << i - 1 << " with \n";
	//if(i - 1 != j){
	for(int k = j-1; k < i - 1; k++)
	  orderHistogram[k] = orderHistogram[k+1];
	//if(i == records){
	orderHistogram[i - 1] = index;
	//}
	}*/
      //cout << " the new order is (";
      // for(int k = 0; k < records; k++)
      //cout << orderHistogram[k] << ",";
      //cout << ")\n";
      } else {
	//cout << " should stop !!!!!!!!!!!!!!!!!!! \n";
	//cout << " the new order is (";
	//for(int k = 0; k < records; k++)
	//cout << orderHistogram[k] << ",";
	//cout << ")\n";
	break;
      }
    }
    //
    bubble(iteratorTemp,orderHistogram,records,dimension);
    //cout << "orderHistogram (";
    //for(iBin = 0; iBin < records; iBin++)
    //  cout << orderHistogram[iBin] << ",";
    //cout << ")\n";
  }

  centerMass = meanStd/meanValue;
  
  populationVariance  = sqrt(populationVariance/meanValue - pow(centerMass, 2.0));

  for(j = 0; j < dimension; j++){
    //cout << "close contur ";
    myFile2Dstd[j].flush();
    //myFile2DContur[j].close();
  }
  for(int k = 0; k < dimension-1; k++)
    for(j = 0; j < dimension - k -1; j++){
      myFile2DContur[k * (dimension-1) - k * (k-1)/2 + j].flush();
      myFile2DContur[k * (dimension-1) - k * (k-1)/2 + j].close();      

      myFile2DDiag[k * (dimension-1) - k * (k-1)/2 + j].flush();
      myFile2DDiag[k * (dimension-1) - k * (k-1)/2 + j].close();      

    }
  /*if(dimension > 3)
    for(int k = 0; k < dimension * (dimension-1) * (dimension - 2) / 6; k++){
      myFile3D[k].flush();
      myFile3D[k].close();      
    }
  */
  delete[] myFile2Dstd;
  delete[] myFile2DContur;
  delete[] myFile2DDiag;

  // delete[] myFile3D;
  delete[] iteratorTotalValue;
  delete[] iteratorTotalNumber;
  delete[] iteratorTotalInverse;
  delete[] iteratorValue;
  delete[] iteratorNumber;
  delete[] iteratorInverseValue;
  delete[] iteratorSqrtValue;
#endif //whole histogram

#ifdef FIRST_ITERATOR
  delete[] iteratorTemp;
#endif //first iterator

  //myFile.flush();
  //myFile.close();
} //end of print

/*
int main (){
  Hashtable_histogram HT(5,20);
  unsigned long proef[5];
  unsigned long trys;
  unsigned long i;
  int it;
  double value;
  HT.listCurrent = new Hashtable_histogram*[5];
  //HT.secondListCurrent = new Hashtable_histogram*[5];
  HT.indexCurrent = new unsigned long[5];
  //HT.secondIndexCurrent = new unsigned long[5];
  
  for(i = 0 ; i < 5 ; i++)
    proef[i] = i;
  //value = 0.5;    
  
  HT.store_individual(proef,0.5);
  trys = HT.check_individual(proef);
  //cout << "answer should be positive check = " << trys << endl;

  //HT.store_individual(proef,value);
  proef[3] = 10;
  
  trys = HT.check_individual(proef);
  //cout << "answer check should be negative = " << trys << endl;
  
  HT.store_individual(proef,0.7);
  proef[3] = 3;
  
  HT.store_individual(proef,0.9);

  proef[1] = 7;
  
  HT.store_individual(proef,0.2);
  
  proef[0] = 7;
  
  HT.store_individual(proef,0.1);
  //iterate the list 
  it = HT.iteratorInit(proef,value);
  cout << " indiv (";
  for(i = 0 ; i < 5 ; i++ )
    cout << proef[i] << ",";
  cout << ") = " << value << " nr Indiv = " << it << "\n";
  while((it = HT.iterator(proef,value)) != -1){
    cout << " indiv (";
    for(i = 0 ; i < 5 ; i++ )
      cout << proef[i] << ",";
    cout << ") = " << value << " nr Indiv = " << it << "\n";
  }
  
  HT.free_table(&HT);
  return 0;
  }*/


