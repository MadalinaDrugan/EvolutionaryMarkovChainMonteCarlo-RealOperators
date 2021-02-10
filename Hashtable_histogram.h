#include<iostream.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include"PARAMETERS.h"

class Hashtable_histogram {
    Hashtable_histogram **histogram;
    unsigned long dimension;
    unsigned long nrBins;
    //indicates the first element from the hashtable
    int first;
 public:
    Hashtable_histogram(unsigned long,unsigned long);
    ~Hashtable_histogram();
    void free_table(Hashtable_histogram*);
    void store_individual(unsigned long*, double);
    unsigned long check_individual(unsigned long*,double&,double&,double&);
    unsigned long delete_individual(unsigned long*);
    void reset();
    unsigned long nrElem;
    unsigned long nrIndividuals;
    double mean;
    double meanSqrt;
    double inverseMean;

    //one iterator
#define FIRST_ITERATOR
#ifdef FIRST_ITERATOR
    // pointers to histograms
    Hashtable_histogram** listCurrent;
    //the values in the current histogram
    unsigned long* indexCurrent;
    long iterator(unsigned long*, double&,double&,double&);
    long iteratorInit(unsigned long*, double&,double&,double&);
    void iteratorFree();
#endif //first iterator

#ifdef SECOND_ITERATOR
    //second iterator
    Hashtable_histogram** listSecondCurrent;
    unsigned long* indexSecondCurrent;
    long iteratorSecond(unsigned long *, double&,double&,double&);
    long iteratorSecondInit(unsigned long*, double&, double&, double&);
#endif //second iterator

    //statistical functions
    void printEachHistogram(Hashtable_histogram**, int);
    int* bubble(unsigned long**,int*,int,int);
    void printAdditionalHistograms(ofstream*, ofstream*, ofstream*, unsigned long*, double*, double*,int);
    void printStd(ofstream&,Hashtable_histogram**, int, double&, double&);
};

