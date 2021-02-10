#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>

#include "long_doubleList.h"

doubleList::doubleList(){
	object = NULL;
	numbers = 0;
	next = NULL;
	index = 0;
}

doubleList::doubleList(double myObject){
	object = new double(myObject);
	numbers = 1;
	next = NULL;
	index = 0;
}

doubleList::doubleList(const doubleList& mydoubleList){
	object = new double(*(mydoubleList.object));
	numbers = mydoubleList.numbers;
	next = (doubleList*)malloc(sizeof(mydoubleList.next));
	next = mydoubleList.next;
	index = mydoubleList.index;
}

void doubleList::setObject(double myObject){
	object = new double(myObject);
	numbers = 1;
}

int doubleList::setObject(double myObject, unsigned long myIndex){
	if(index == myIndex) {
		object = new double(myObject);
		numbers = 1;
		return 0;
	}	
	if(next == NULL) return -1;
	else return next->setObject(myObject,myIndex);
}

double doubleList::getObject(){
	double myObject = *object;
	return myObject;
}

int doubleList::addElement(double myObject){
	if(object == NULL) {
		object = new double(myObject);
		numbers = 1;
		index = 0;
		return 1;
	}
	
	if(myObject == *object) {
	  numbers++;
	  return index;
	}

	if(next == NULL) 
	{
		doubleList* temp = new doubleList(myObject);
		temp->index = index + 1;
		temp->numbers = 1;
		temp->next = NULL;
		next = temp;
		return temp -> index;
	} 
	  
	return next->addElement(myObject);
}


doubleList::~doubleList(){
	if(next != NULL) delete next;  
	free(object);
	//delete this;
}

void doubleList::print(){
	if(object == NULL) {
		cout << "emply list";
		return;
	}	
	cout << "(" << *object <<"," << index << ") ";
	if(next!= NULL) next->print();	
}

//ostream& doubleList::operator<<(doubleList* mydoubleList){
//	cout << mydoubleList->object <<" ";
//	if(mydoubleList->next!= NULL) cout << mydoubleList->next;	
//}

double doubleList::operator[](int number){
	if(index == number) return *object;	
	else if(index < number && next != NULL) 
		{
			double temp = (*next)[number];
			return temp;
		}
	else return -1; 
}

double doubleList::getObject(unsigned long number){
	if(index == number) return *object;	
	else if(index < number && next != NULL) 
			return next -> getObject(number);
	else return 0; 
}

unsigned long doubleList::operator()(int number){
	if(index == number) return numbers;	
	else if(index < number && next != NULL) 
			return (*next)(number);
	else return -1; 
}

unsigned long doubleList::getNumbers(unsigned long number){
	if(index == number) return numbers;	
	else if(index < number && next != NULL) 
			return next -> getNumbers(number);
	else return 0; 
}

int doubleList::operator==(doubleList& mydoubleList){
	if(next != NULL)
		{if(&object == &(mydoubleList.object) && *next == *(mydoubleList.next)) return 1;}
	else if(*object == *(mydoubleList.object)) return 1;	
	else return 0;
}

int doubleList::operator!=(doubleList& mydoubleList){
	if(next != NULL)
		{if(*object != *(mydoubleList.object) || *next != *(mydoubleList.next)) return 1;}
	else if(*object != *(mydoubleList.object)) return 1;
	else return 0;
}

int doubleList::contains(double myObject){
	if(myObject == *object) return 1;
	if(next == NULL) return 0; 
	return next -> contains(myObject);
}

unsigned long doubleList::getIndex(double myObject){
	if(object == NULL) {
		object = new double(myObject);
		index = 0;
		return index;
	}
	if(myObject != *object) return index;
	if(next == NULL) {
		addElement(myObject);
		return next->index;
	} 
	return next -> getIndex(myObject);	
} 
//int doubleList::operator=(doubleList& mydoubleList)

unsigned long doubleList::getSize(){
	if(object == NULL) return 0;
	if(next != NULL) return next -> getSize() + 1;
	return 1;
}
