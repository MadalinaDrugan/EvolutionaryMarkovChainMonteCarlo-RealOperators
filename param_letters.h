typedef struct parameters_struct parameters;
struct parameters_struct {
int nr_simulations;
int nr_vars;
int first_bias;
int adapt_fitness;
int adapt_fitness_map;
double adapt_map_factor;
int memory_used;
double optimum;
int local_iter;
int function;
int local_search;
int pop_size;
int nr_iterations;
double mut_rate;
double cross_rate;
int elitist;
int nr_generations;
int init_method;
int tournament_size;
char *output_fname;
char *output_fname_all;
};


