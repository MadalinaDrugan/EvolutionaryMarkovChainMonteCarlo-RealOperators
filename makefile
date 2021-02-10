######################
# Evolutionary Markov chain Monte Carlo @ 2005
# Source code for the article 
# Author: Madalina Drugan
######################### 

CC          = g++ -Wall -g  	 # Compiler is gcc


#------------------------------------------------------------------------------
# Action
#------------------------------------------------------------------------------

MCMC:  mcmc_2_binary_trap.o binaryPopulation_trap.o binaryMCMC_trap.o param_letters_read.o utils.o SIMUL.o Hashtable.o long_doubleList.o convergence.o random.o
	$(CC) -o MCMC mcmc_2_binary_trap.o binaryPopulation_trap.o utils.o param_letters_read.o SIMUL.o Hashtable.o binaryMCMC_trap.o long_doubleList.o convergence.o random.o

mcmc_2_binary_trap.o: mcmc_2_binary_trap.cpp binaryPopulation_trap.h random.cpp
	$(CC) -c mcmc_2_binary_trap.cpp 

binaryMCMC_trap.o: binaryMCMC_trap.cpp binaryMCMC_trap.h PARAMETERS.h random.cpp
	$(CC) -c binaryMCMC_trap.cpp 

binaryPopulation_trap.o: binaryPopulation_trap.cpp SIMUL.h binaryPopulation_trap.h binaryMCMC_trap.h PARAMETERS.h Hashtable.h convergence.cpp random.cpp
	$(CC) -c binaryPopulation_trap.cpp

utils.o: utils.c utils.h
	$(CC) -c utils.c

Hashtable.o: Hashtable.cpp Hashtable.h
	$(CC) -c Hashtable.cpp

SIMUL.o: SIMUL.c SIMUL.h
	$(CC) -c SIMUL.c

param_letters_read.o: param_letters_read.c param_letters.h utils.h
	$(CC) -c param_letters_read.c

convergence.o: convergence.cpp binaryPopulation_trap.h binaryMCMC_trap.h PARAMETERS.h long_doubleList.h random.cpp
	$(CC) -c convergence.cpp

long_doubleList.o: long_doubleList.cpp long_doubleList.h
	$(CC) -c long_doubleList.cpp 

random.o: random.cpp PARAMETERS.h
	$(CC) -c random.cpp


