#include<stdio.h>
#include<string>
#include <stdlib.h>
#include <iostream.h>

#include "Hashtable.h"

using namespace std;

Hashtable :: Hashtable(unsigned long nr_vars)
{
    one = NULL;
    zero = NULL;
    nrElem = 0;
    nrIndividuals = 0;
    first = 0;

    listCurrent = NULL;
    indexCurrent = NULL;

    secondListCurrent = NULL;
    secondIndexCurrent = NULL;
    //first = 1;
}

Hashtable :: ~Hashtable()
{
    if(first == 0){
	delete[] listCurrent;
	delete[] indexCurrent;

	delete[] secondListCurrent;
	delete[] secondIndexCurrent;
    }	
}

void Hashtable :: free_table(Hashtable *table)
{
  if (table == NULL)
  return;

  free_table(table->zero);
  free_table(table->one);
  if(first != 0){
      delete table;
      table = NULL;
      //delete listCurrent;
      //if(indexCurrent != NULL) delete indexCurrent;
  } else {
      zero = NULL;
      one = NULL;
      nrIndividuals = 0;
      nrElem = 0;
  }

  //free(table);
  // free(listCurrent);
  //free(indexCurrent);

  //delete table;
  //table == NULL;
}


void Hashtable :: reset(unsigned long nr_vars)
{
    if(this == NULL) {
	cout << "Hashtable is empty\n";
	return;
    }
    if(listCurrent != NULL) 
      for(unsigned long i = 0; i < nr_vars; i++){
	listCurrent[i] = NULL;
	indexCurrent[i] = 0;
	
	secondListCurrent[i] = NULL;
	secondListCurrent[i] = 0;
      }

    free_table(zero);
    free_table(one);

    zero = NULL;
    one = NULL;

    nrElem = 0;
    nrIndividuals = 0;
}

unsigned long Hashtable :: check_individual(unsigned long* indiv, unsigned long nr_vars)
{
  unsigned long i;
  unsigned long output;
  Hashtable *table;
  unsigned long nrInds;
  
  table = this;
  
  //cout << "Cercetat ";
  //for(unsigned long i = 0; i < nr_vars; i++)
  //  cout << indiv[i];
  //cout << "\t";

  for(i = 0 ; i < nr_vars ; i++)
    {
      if (indiv[i] == 0)
	{
	  if (table->zero == NULL)
	    return 0;
	  else
	    {
	      //cout << "(" <<indiv[i] << "= 0," << i << ","<< table-> nrIndividuals <<")";
	      table = table->zero;
	      nrInds = table->nrIndividuals;
	    }
	}
      else
	if (indiv[i] == 1)
	  {
	    if (table->one == NULL)
	      return 0;
	    else
	      {
	        //cout << "(" <<indiv[i] << "= 1," << i << ","<< table-> nrIndividuals <<")";
		table = table->one;
 		nrInds = table->nrIndividuals;
	      }
	  }
    }
  //cout << "\n";
  return nrInds;
}

void Hashtable :: store_individual(unsigned long* indiv, unsigned long nr_vars)
{
  unsigned long i,j;
  unsigned long output;
  Hashtable *table,*firsttable;
  
  /*cout << "Cercetat ";
  for(unsigned long i = 0; i < nr_vars; i++)
    cout << indiv[i];
    cout << "\t";*/

  firsttable = this;
  table = this;
  for(i = 0 ; i < nr_vars ; i++)
    {
      if (indiv[i] == 0)
	{
	  if (table->zero == NULL)
	    {
	      table->zero = new Hashtable(nr_vars);
	      table = table->zero;
	      table->first = 1; 
	      table-> nrIndividuals++;
	      //firsttable-> nrIndividuals++;
	      firsttable->nrElem++;

	      //cout << "nou ("<< indiv[i] << "= 0,"<< i << "," << table-> nrIndividuals << ")";
	      for(j = i+1 ; j < nr_vars ; j++)
		{
		  if (indiv[j] == 0)
		    {
		      table->zero = new Hashtable(nr_vars);
		      table = table->zero;
		      table -> first = 1;
		      table-> nrIndividuals++;
		      //cout << "(" <<indiv[j] << "= 0," << j << ","<< table-> nrIndividuals <<")";
		    }
		  else if (indiv[j] == 1)
		    {
		      table->one = new Hashtable(nr_vars);
		      table = table->one;
		      table -> first = 1;
		      table-> nrIndividuals++;
		      //      cout << "(" <<indiv[j] << "= 1," << j << ","<< table-> nrIndividuals <<")";
		    }
		}
	      //cout << "\n";
	      return;
	    }
	  else
	    {
	      table = table->zero;
	      table -> first = 1;
	      table->nrIndividuals++;
	      // cout << "(" <<indiv[i] << "=0," << i << ","<< table-> nrIndividuals << ")";
	    }
	}
      else
	if (indiv[i] == 1)
	  {
	    if (table->one == NULL)
	      {
		table->one = new Hashtable(nr_vars);
		table = table->one;
		table -> first = 1;
		table->nrIndividuals++;
  	        //firsttable-> nrIndividuals++;
		firsttable->nrElem++;
	        //cout << "nou ("<< indiv[i] << "=1,"<< i << "," << table-> nrIndividuals << ")";
		for(j = i+1 ; j < nr_vars ; j++)
		  {
		    if (indiv[j] == 0)
		      {
			table->zero = new Hashtable(nr_vars);
			table -> first = 1;
			table = table->zero;
			table-> nrIndividuals++;
		        //cout << "(" <<indiv[j] << "=0," << j << ","<< table-> nrIndividuals <<")";
		      }
		    else if (indiv[j] == 1)
		      {
			table->one = new Hashtable(nr_vars);
			table -> first = 1;
			table = table->one;
			table-> nrIndividuals++;
			//      cout << "(" <<indiv[j] << "=1," << j << ","<< table-> nrIndividuals <<")";
		      }
		  }
		//cout << "\n";
		return;
	      }
	    else{
	      table = table->one;
	      table->nrIndividuals++;
	      // cout << "(" <<indiv[i] << "=1," << i << ","<< table-> nrIndividuals << ")";
	    }
	  }
    }
  //cout << "\n";
}

unsigned long Hashtable :: delete_individual(unsigned long* indiv, unsigned long nr_vars)
{
unsigned long i;
unsigned long output;
Hashtable *table;
 unsigned long nrInds;

 table = this;

for(i = 0 ; i < nr_vars ; i++)
 {
   if (indiv[i] == 0)
    {
     if (table->zero == NULL)
	return 0;
     else
       {
         table->nrIndividuals --;
	 nrInds = table->nrIndividuals;
	 table = table->zero;
       }
    }
   else
     if (indiv[i] == 1)
       {
	 if (table->one == NULL)
	   return 0;
	 else
	   {
	     table->nrIndividuals --;
	     nrInds = table->nrIndividuals;
	     table = table->one;
	   }
       }
 }
 return nrInds;
}

long Hashtable :: iteratorInit(unsigned long* indiv, unsigned long nr_vars)
{
    //listCurrent = new Hashtable*[nr_vars];
    //indexCurrent = new unsigned long[nr_vars];

  if(this == NULL)
    return -1;
  
  listCurrent[0] = this;
  
  unsigned long i =0;
  cout << "indiv=";
  while(i < nr_vars){
    if(listCurrent[i]->zero != NULL){
      indexCurrent[i] = 0;
      //cout << "(0" << listCurrent[i] ->nrIndividuals << ")";
      indiv[i] = 0;
      if(i != nr_vars) 
	listCurrent[i+1] = listCurrent[i] -> zero;
      i++;
    } else if(listCurrent[i]->one != NULL){
      indexCurrent[i] = 1;
      indiv[i] = 1;
      //cout << "(1" << listCurrent[i] ->nrIndividuals << ")";
      if(i != nr_vars) 
	listCurrent[i+1] = listCurrent[i] -> one;
      i++;
    } else {
      if(i==0){
	cout << "greaseala in Hashtable.cpp:301 for " << i << " empty list ";
	
	cout << "\n";
	return -1;
      }
      cout << "Mergi in spate la i=" << i << " greseala in Hashtable.cpp:306" << listCurrent[i]->one << "\n";
      exit(1);
      while(i>0){ 
	i--;
	if(indexCurrent[i] == 0 && listCurrent[i] -> one != NULL){
	  indiv[i] = 1;
	  indexCurrent[i] = 1;
	  listCurrent[i+1] = listCurrent[i]->one;
	  cout << "Mergi in spate la i=" << i << " greseala in Hashtable.cpp:239\n";
	  exit(1);
	}
      }
    }
  }
  //cout << " = " << listCurrent[nr_vars]->nrIndividuals << "\n";
  return listCurrent[nr_vars]->nrIndividuals;
}


void Hashtable :: iteratorFree(unsigned long nr_vars)
{
//    if(indexCurrent != NULL)
//	delete[] indexCurrent;
}

long Hashtable :: iterator(unsigned long* indiv, unsigned long nr_vars)
{
    unsigned long i=0;
    //Hashtable* table = listCurrent[nr_vars-1];
    //cout << "indiv precedent=";
    for(i = 0; i<nr_vars; i++){
	indiv[i] = indexCurrent[i];
	//cout << indiv[i];
    }
    //cout << "\n";

    //if(nrElem == 0 || nrIndividuals == 0)
    // return -1;

    i = nr_vars;

    //inapoi
    if (indexCurrent[nr_vars] == 1 || listCurrent[nr_vars]->one == NULL){
      while(i>0){ 
	  i--;
	  if(indexCurrent[i] == 0 && listCurrent[i] -> one != NULL){
	      indiv[i] = 1;
	      indexCurrent[i] = 1;
	      if(i != nr_vars) listCurrent[i+1] = listCurrent[i]->one;
	      //cout << "inapoi i=" << i << "actual 1\n";
	      i++;
	      break;
	  }
	  if(i==0) {
	      //cout << "trebuie sa te termine\n";
	      return -1;
	  }
	  //cout << "i" << i << "\n";
      }
      //cout << "la sfirsti i=" <<i <<"\n";
    } else {
	indexCurrent[i] = 1;
	indiv[i] =1;
	//cout << "change i=" << i << "actual 1";
	return listCurrent[nr_vars]->nrIndividuals;
    }

  //inaunsigned longe
    //cout << "inaunsigned longe from i=" << i << "\t";
  while(i < nr_vars){
      if(listCurrent[i]->zero != NULL){ 
	  indexCurrent[i] = 0;
	  indiv[i] = 0;
	  //cout << "(0" << listCurrent[i]->nrIndividuals << ") ";
	  if(i != nr_vars) listCurrent[i+1] = listCurrent[i]->zero;
	  i++;
      }else if(listCurrent[i]->one != NULL){
	    indexCurrent[i] = 1;
	    indiv[i] = 1;
	    //cout << "1" << listCurrent[i]->nrIndividuals << ") ";
	    if(i != nr_vars) listCurrent[i+1] = listCurrent[i] -> one;
	    ++i;
      } else {
	  cout << "Greseala de constrictie in hashtable pentru variabla "<<i 
	       << "in Hashtable.cpp:312\n";
	  exit(1);
      }
  } 
  //cout << "\n";
  return listCurrent[nr_vars]->nrIndividuals;
}

long Hashtable :: iteratorSecondInit(unsigned long* indiv, unsigned long nr_vars)
{
  unsigned long i = 0;
  
  if(this == NULL)
    return -1;

  secondListCurrent[0] = this;
  
  //cout << "Cercetat ";
  //for(unsigned long i = 0; i < nr_vars; i++)
  //  cout << indiv[i];
  //cout << "\t";

  for(i = 0 ; i < nr_vars ; i++)
    {
      if (indiv[i] == 0)
	{
	  if (secondListCurrent[i]->zero == NULL)
	    return 0;
	  else
	    {
	      //cout << "(" <<indiv[i] << "= 0," << i << ","<< table-> nrIndividuals <<")";
	      if(i != nr_vars)
		secondListCurrent[i+1] = secondListCurrent[i]->zero;
	      secondIndexCurrent[i] = 0;
	    }
	}
      else
	if (indiv[i] == 1)
	  {
	    if (secondListCurrent[i]->one == NULL)
	      return 0;
	    else
	      {
	        //cout << "(" <<indiv[i] << "= 1," << i << ","<< table-> nrIndividuals <<")";
		if(i != nr_vars)
		  secondListCurrent[i+1] = secondListCurrent[i]->one;
 		secondIndexCurrent[i] = 1;
	      }
	  }
    }
  //cout << "\n";
  return listCurrent[nr_vars]->nrIndividuals;

}

long Hashtable :: iteratorSecond(unsigned long* indiv, unsigned long nr_vars)
{
    unsigned long i=0;
    //Hashtable* table = listCurrent[nr_vars-1];
    //cout << "indiv precedent=";
    for(i = 0; i<nr_vars; i++){
	indiv[i] = secondIndexCurrent[i];
	//cout << indiv[i];
    }
    //cout << "\n";

    i = nr_vars;

    //inapoi
    if (secondIndexCurrent[nr_vars] == 1 || secondListCurrent[nr_vars]->one == NULL){
      while(i>0){ 
	  i--;
	  if(secondIndexCurrent[i] == 0 && secondListCurrent[i] -> one != NULL){
	      indiv[i] = 1;
	      secondIndexCurrent[i] = 1;
	      if(i != nr_vars) secondListCurrent[i+1] = secondListCurrent[i]->one;
	      //cout << "inapoi i=" << i << "actual 1\n";
	      i++;
	      break;
	  }
	  if(i==0) {
	      //cout << "trebuie sa te termine\n";
	      return -1;
	  }
	  //cout << "i" << i << "\n";
      }
      //cout << "la sfirsti i=" <<i <<"\n";
    } else {
	secondIndexCurrent[i] = 1;
	indiv[i] =1;
	//cout << "change i=" << i << "actual 1";
	return secondListCurrent[nr_vars]->nrIndividuals;
    }

  //inaunsigned longe
    //cout << "inaunsigned longe from i=" << i << "\t";
  while(i < nr_vars){
      if(secondListCurrent[i]->zero != NULL){ 
	  secondIndexCurrent[i] = 0;
	  indiv[i] = 0;
	  //cout << "(0" << listCurrent[i]->nrIndividuals << ") ";
	  if(i != nr_vars) secondListCurrent[i+1] = secondListCurrent[i]->zero;
	  i++;
      }else if(secondListCurrent[i]->one != NULL){
	    secondIndexCurrent[i] = 1;
	    indiv[i] = 1;
	    //cout << "1" << listCurrent[i]->nrIndividuals << ") ";
	    if(i != nr_vars) secondListCurrent[i+1] = secondListCurrent[i] -> one;
	    ++i;
      } else {
	  cout << "Greseala de constrictie in hashtable pentru variabla "<<i 
	       << "in Hashtable.cpp:312\n";
	  exit(1);
      }
  } 
  //cout << "\n";
  return secondListCurrent[nr_vars]->nrIndividuals;
}

/*
unsigned long main ()
{
Hashtable HT(10);
unsigned long proef[10];
unsigned long trys;
unsigned long i;

  for(i = 0 ; i < 10 ; i++)
    proef[i] = 1;
    

  HT.store_individual(proef);
   proef[9] = 0;
  trys = HT.check_individual(proef);
  cout << "answer = " << trys << endl;

return 0;
}
*/
 

