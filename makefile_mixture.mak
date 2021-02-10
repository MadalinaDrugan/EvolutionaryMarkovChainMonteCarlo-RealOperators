# Definition

CXX         = g++ -Wall -g  	 # Compiler is gcc
CXXFLAGS = -O2 -Wall

DIFF = ./sdiff
PRE = ./
MAJOR = 1
MINOR = 0

#CC          = g++ -Wall -O3           # Compiler is gc

%.o:           	%.cpp
		$(CXX) $(CXXFLAGS) -c $*.cpp

#------------------------------------------------------------------------------
# Aktion
#------------------------------------------------------------------------------

MCMC:  mcmc_2_binary_trap.o binaryPopulation_trap.o binaryMCMC_trap.o param_letters_read.o utils.o SIMUL.o Hashtable.o Hashtable_histogram.o long_doubleList.o convergence.o random.o libnewmat.a
	$(CXX) -o MCMC mcmc_2_binary_trap.o binaryPopulation_trap.o utils.o param_letters_read.o SIMUL.o Hashtable.o Hashtable_histogram.o binaryMCMC_trap.o long_doubleList.o convergence.o random.o  -L. -lnewmat -lm

mcmc_2_binary_trap.o: mcmc_2_binary_trap.cpp binaryPopulation_trap.h random.cpp 
	$(CXX) -c mcmc_2_binary_trap.cpp 

binaryMCMC_trap.o: binaryMCMC_trap.cpp binaryMCMC_trap.h PARAMETERS.h random.cpp newmatap.h newmatio.h newmat.h include.h boolean.h myexcept.h
	$(CXX) -c binaryMCMC_trap.cpp 

binaryPopulation_trap.o: binaryPopulation_trap.cpp SIMUL.h binaryPopulation_trap.h binaryMCMC_trap.h PARAMETERS.h Hashtable.h Hashtable_histogram.h convergence.cpp random.cpp 
	$(CXX) -c binaryPopulation_trap.cpp

utils.o: utils.c utils.h
	$(CXX) -c utils.c

Hashtable.o: Hashtable.cpp Hashtable.h
	$(CXX) -c Hashtable.cpp

Hashtable_histogram.o: Hashtable_histogram.cpp Hashtable_histogram.h
	$(CXX) -c Hashtable_histogram.cpp

SIMUL.o: SIMUL.c SIMUL.h
	$(CXX) -c SIMUL.c

param_letters_read.o: param_letters_read.c param_letters.h utils.h
	$(CXX) -c param_letters_read.c

convergence.o: convergence.cpp binaryPopulation_trap.h binaryMCMC_trap.h PARAMETERS.h long_doubleList.h random.cpp 
	$(CXX) -c convergence.cpp

long_doubleList.o: long_doubleList.cpp long_doubleList.h
	$(CXX) -c long_doubleList.cpp 

random.o: random.cpp PARAMETERS.h
	$(CXX) -c random.cpp

parametri_matrix:      parametri_matrix.o libnewmat.a
		$(CXX) -o parametri_matrix.o -L. -lnewmat -lm

parametri_matrix.o: parametri_matrix.cpp newmatap.h newmatio.h newmat.h include.h boolean.h myexcept.h random.cpp



