#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#include "param_letters.h"
 #include "utils.h"
//#include "utils.c"

#define MAX_STRING_LENGTH 256

extern word_list line_to_words(char*);

parameters *read_parameters(char *filename) {

  parameters *par;
  FILE *fp;
  char line[MAX_STRING_LENGTH];
  word_list words;
  int nr_pars;
  int i;
 
  par = (parameters*)malloc(sizeof(parameters));
  if (par == NULL)
     {
    printf("couldn't allocate memory for parameters !\n");
     exit(0);
    }
 
  fp = fopen(filename,"r");
  if (fp == NULL) {
    printf("read_control_parameters: Couldn't find file %s...\n",filename);
    exit(-1);
  }
  while (NULL != fgets(line, MAX_STRING_LENGTH, fp))
    if (line[0] != '#') { /* ignore comments */
      words = line_to_words(line);
      if (words.number_of_words != 0) { /* ignore empty lines */
        if (!strcmp(words.word[0], "nr_simulations:")) 
	  {
	  par->nr_simulations = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "memory_used:")) 
	  {
	  par->memory_used = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "adapt_fitness:")) 
	  {
	  par->adapt_fitness = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "adapt_fitness_map:")) 
	  {
	  par->adapt_fitness_map = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "first_bias:")) 
	  {
	  par->first_bias = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "adapt_map_factor:")) 
	  {
	  par->adapt_map_factor = atof(words.word[1]);
          }
        else if (!strcmp(words.word[0], "optimum:")) 
	  {
	  par->optimum = atof(words.word[1]);
          }
        else if (!strcmp(words.word[0], "pop_size:")) 
	  {
	  par->pop_size = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "local_iter:")) 
	  {
	  par->local_iter = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "function:")) 
	  {
	  par->function = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "local_search:")) 
	  {
	  par->local_search = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "mut_rate:")) 
	  {
	  par->mut_rate = atof(words.word[1]);
          }
        else if (!strcmp(words.word[0], "cross_rate:")) 
	  {
	  par->cross_rate = atof(words.word[1]);
          }
        else if (!strcmp(words.word[0], "elitist:")) 
	  {
	  par->elitist = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "nr_generations:")) 
	  {
	  par->nr_generations = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "tournament_size:")) 
	  {
	  par->tournament_size = atoi(words.word[1]);
          }
        else if (!strcmp(words.word[0], "nr_vars:")) 
	  {
	  par->nr_vars = atoi(words.word[1]);
          }

        else if (!strcmp(words.word[0], "output_fname:"))
	  {
	    par->output_fname = strdup(words.word[1]);
          }
        else if (!strcmp(words.word[0], "output_fname_all:"))
	  {
	    par->output_fname_all = strdup(words.word[1]);
          }
	 else printf("Do not know parameter : %s\n",words.word[0]);
      }
    }
  fclose(fp);
  return(par);

} /* read_control_parameters */



