///PARAMETERS for binaryMCMC and binaryPopulation
//#define Boltzmann 
//#define Normal_Boltzmann

#define bestTEMPERATURE 1
#define worstTEMPERATURE 1

//option for functions with integer values
#ifndef MIXTURE_REALS
//#define FLOOR
#endif
//option for detailed debuging
//#define DEBUG_

//option for debuging
//#define RUN_SEQENCES

//option for scaling the temperature every generation to prevent that all the chains are at same temperature
//#define WITHOUT_MEM
//#define NORMAL_TEMPERING
//#define PARALLEL_SIMULATED_ANNEALING_TEMPERING
//#define COMPLICAT_TEMPERING
//#define FOARTE_COMPLICAT

//cooling types
//#define LINEAR_COOLING
//#define GEOMETRIC_COOLING
//#define LOGARITHMIC_COOLING
#define NO_COOLING
#define DEFAULT_TEMPERATURE 1.0

//no acceptance
//#define RANDOM_WALK

//types of algorithms that we deal with
//#define RECOMB_MUT
#ifdef RECOMB_MUT
#define BOTH_GOOD
#endif

//#define NO_COUPLING
/////#define SINGLE_vs_MULTIPLE

//#define RANDOM_COUPLING

#define PARALLEL_TEMPERING
//#define SWAP_CHAIN
#define SWAP 0.0          

//#define ELITIST_ACCEPTANCE

//#define ELITIST_ACCEPTANCE_WITH_REGENERATION
//the parents are also considered in the selection process
//#define WITH_PARENTS
#define SIMPLE_SELECTION
//#define MAX_SELECTION
//#define GA_SELECTION
//#define DETAILED_BALANCE
#define IRREDUCTIBIL

//#define PARALLEL_SIMULATED_ANNEALING
#define COUPLED_ACCEPTANCE
#define COUPLING_PROBABILITY 1.0
#define JUMP_SINGLE
//#define SUM_DISTRIB
#define WRITE_OUT_FILE

//#define POPULATION_MCMC
#define METROPOLIS_ALGORITHM
//#define METROPOLIS_HASTINGS_ALGORITHM
//#define INDEPENDENT_SAMPLER
//#define EXPLOITATION_INDEPENDENT_SAMPLER
//#define EXPLORATION_INDEPENDENT_SAMPLER
//

#define IsMutation
#define RECOMBINATION_PROB 1.0
#define IsRecombination
#ifdef IsRecombination

//real recombination or binary
//#define BINARY_RECOMB
#ifdef BINARY_RECOMB
#define NO_RECOMB 0
#define ONE_POINT_RECOMB 1
#define TWO_POINTS_RECOMB 2

#define N_POINT_RECOMB 12
//#ifndef MULTIPLE_RUNS
#define PROB_nPOINT 1
//#endif

//#define RECOMBINATION_TEST
#ifdef RECOMBINATION_TEST
#define RECOMBINATION_TEST_PARAM 6
#endif //#ifdef RECOMBINATION_TEST

#define UNIFORM_RECOMB 3
//#ifndef MULTIPLE_RUNS
#define UNIF_RECOMB_PARAM 0.5
//#endif

//#define RECOMBINATION_PROB 1.0
#define UNIFORM_RECOMB_PROBABIL 4
//#define TWO_CHILDREN
//#define MIXTURE 0.5
//#define CYCLE

#else //not binary recomb

#define PARENT_CENTRIC 1
#define MIXTURE 0.7

//#define TWO_PARENTS_SIMPLE

//#ifdef TWO_PARENTS
#define sigma_m 1
//#endif //two parents

//#define CYCLE 0.5
//
#define TRANSL_ROTATE_SCALE

//#define EUCLIDIAN

#ifdef TRANSL_ROTATE_SCALE

#define TRANSLATION
#ifdef TRANSLATION
#define sigma_transl 0.01
//#define TRANSLATION_COMPLEX
#ifdef TRANSLATION_COMPLEX
#define sigma_transl_oriz 0.25
#ifndef TWO_PARENTS
#define SIMPLE_PERP_TRANSL
#endif
#endif //translation complex
//#define PCTX
#ifdef PCTX 
//#define SUM_TRANSL 
//not implemented completely - see the proposal probability
//#define wPCTX
#endif 
//#define SYMMETRICAL_SNOOKER
#ifndef PCTX
//#define SIMPLEX
#ifdef SIMPLEX
//#define 
#endif //simplex

//#define HAFL_PCTX
#ifdef HALF_PCTX
//
#else //half pctx
//#define TerBraak
#ifdef TerBraak
#define constTerBraak 0.5
#endif  //terBraak
#endif //half PCTX
#endif 
#endif 

#define SINGLE_INDIV_ROT
#define sigma_rotate 0.01
#define sigma_scale 0.01

#define ROTATION
#ifdef ROTATION
#define SUM_ROT
//#define ROTATION_PROP 
//#define sigma_rotate 0.750
#endif 

//#define SCALING
#ifdef SCALING
#define sigma_scale 1.0
#define SCALING_MUTATION_SAME
#endif

#define TRANSL_ROT
//#define MIXTURE_TRANSL_ROT 0.75
//#define TRANSL_ROT_TRANSL
//#define ROT_TRANSL_ROT
#define MIXTURE_TRANSL_ROT_SCALE 1.0
//#define CYCLE_TRANSL_ROT_SCALE 0.5
//#define MIXTURE_TRANSL_DIFF_ROT 0.5
//#define SUM_TRANSL_ROT_SCALE
//#define MIXTURE_TRANSL_SCALING 0.8
//#define SYMMETRIC_PROP

#endif //TRANSL_ROTATE_SCALE

#endif // binary recomb
#endif

//#define TRAP_LIKE_FUNCTION
#ifdef TRAP_LIKE_FUNCTION
//#define MULTI_PICK_TRAP_FUNCTION
//#define a 0.75
//#define b 1
//#define BLOCKsize 6
//#define z 1

#define TRAP_FUNCTION
#ifdef TRAP_FUNCTION
#define a 4
#define b 3
#define c 2
#define BLOCKsize 4
#define z (BLOCKsize - 1)

#define LINEAR_TRAP_FUNCTION
#ifndef LINEAR_TRAP_FUNCTION
#define ORDER_BLOCK 1
#endif //LINEAR_TRAP_FUNCTION

//#define THRESHOLD
#ifdef THRESHOLD
//#define HALF_MAX_FITNESS
//#define THRESHOLD_DISTRIBUTION
#endif
#endif //TRAP_FUNCION

//bimomial functions
//#define BINOMIAL
#ifdef BINOMIAL
#define p_binomial 0.1
#define BLOCKsize 4
#define LINEAR_TRAP_FUNCTION
#ifndef LINEAR_TRAP_FUNCTION
#define ORDER_BLOCK 1
#endif //LINEAR_TRAP_FUNCTION
#define THRESHOLD
#ifdef THRESHOLD
//#define HALF_MAX_FITNESS
#define THRESHOLD_DISTRIBUTION
#endif
#endif

#endif //trap like functiona

//#define ONEMAX_FUNCTION

//#define TWO_ATTRACTORS_FUNCTION
#ifdef TWO_ATTRACTORS_FUNCTION
#define wrongBits 4
#define a wrongBits*2
#endif

//bernoulli functions
//#define BERNOULLI

//qudratic bipartitioning problem
//#define QBP
#ifdef QBP
#define FileQBP "qbp.data"
#define generatedQBP
#endif

//mixture real numbers
#define MIXTURE_REALS
#ifdef MIXTURE_REALS

#define REALS
#define BLOCKsize 10 // precizia histogramei approximative

//#define RANDOM_MUTATION
#ifndef RANDOM_MUTATION
//#define OPTIMAL_MUTATION
#define mutation_sigma 0.01
#endif //random mutation

#define pi 3.1415926535
#define e 2.718281828
//#define HISTOGRAM
#ifdef HISTOGRAM
#define hist_points 60
#else //not histogram
#define NrBins 60
#endif //histogram

#define scale 8

//#define UNIFORM_DISTRIBUTION

#define BIVARIATE
#ifdef BIVARIATE

#define TWO_PICKS

#define MIXTURE_BIVARIATE
#ifdef MIXTURE_BIVARIATE
// mixture of 9 bivariate

#ifdef NINE_BIVARIATE

#define sigma_1 0.5
#define sigma_2 1
#define rho1 0.5
#define coef1 1.0

#define sigma_3 1
#define sigma_4 1
#define rho3 0.7
#define coef3 0.9

#define sigma_5 0.6
#define sigma_6 0.8
#define rho5 0.9
#define coef5 0.8

#else
// multivariate is the default
//#define mixture 1
#endif //9_BIVARIATE

#else
#ifndef TWO_PICKS
//one bivariate
#define sigma_1 0.1
#define sigma_2 1
#define rho 0.95

#else
//mixture two bivariate

#define sigma_11 0.5
#define sigma_12 1
#define rho_1 0.9
#define coef1 1.0

#define sigma_21 1
#define sigma_22 0.5
#define rho_2 0.9
#define coef2 1.0

#endif //two picks

#endif //mixture bivariate

// for bounded real spaces
//#define CYCLE_SPACE
#define INIT_HALF
#ifdef INIT_HALF
#define init_coef 2.0
#endif 
//#define INIT_CENTER
#ifdef INIT_CENTER
#define init_coef 2.0
#endif //init_center
#define INIT_SIDE
#ifdef INIT_SIDE
#define init_coef 2.0
#endif //init_center
#define CLOSED_EQ_SPACE
#ifdef CLOSED_EQ_SPACE
#define MIN_FITNESS 0.01
#endif //CLOSED_EQ_SPACE

#else // not bivariate
#define ROSENBROCK
#endif //BIVARIATE

#define KNOWN_MAXIM
#endif

#define mutation_multi 1.0
//#define NEAR_CHANGING

#define sizeCoupling 4
//#define TWO_PARENTS
#define coupling 1

#define SAMPLE

#define sampleSize_Const 80100
//double sampleSize = 2000;
//#define MULTIPLE_RUN_SAMPLE
//#define TRACE_RUN
//#ifdef TRACE_RUN
#define nr_runs 10
#define BURN_IN_GEN
#ifdef BURN_IN_GEN
#define throw_gen 40000
#endif
#define multiple_restarts
//#else
//#define nr_runs 10000
//#endif //trace run

#define SEE_SAMPLE 5000

//#define SIMULATION

//#define FILE_BOOKEEPING

//////////////////
//////convergence methods
/////////////////
//measure distance between two peaks
#define DISTANCE

//#define THRESHOLD

//KL distance to the true distribution
//#define KULLBACK_INFORMATION
#define EXPAND_KULLBACK_INFORMATION
#define TIME_PERFORMANCE
//#define KL_DIFFERENCE
//#define PENALITY

//Aceptance rate
#define ACCEPTANCE_RATIO
#define multiplication 1
#define AwayFromZero 0.001

//plots for external analyses
#define PLOTS_MCMC

//multiple chain
//#define GELMAN_RUBIN
#ifdef GELMAN_RUBIN
//#define HALF_GELMAN
#define GELMAN_Moment 3
//#define corection ()
#endif

//metoda CUSUM
#define CUSUM

//#define BURN_IN
#define constBurn 1

#define BURN_IN_STABLE
#define epsilon 0.05

//#define BURN_IN_NORMAL
#define CREDIBIL 0.95
#define CI 1.96

//#define CONVERGENCE_CUSUM
//#define NORMAL_FUNCTION
#define MEAN_UNKNOWN
//#define STICKY

//metoda regenerari
//#define REGENERATION

//
//#define SPECTRAL
////////////

//fisiere folosite
#ifdef RECOMB_MUT
#define specific "recomb_mut" 
#define CusumFILE "cusum_recomb_mut.file"
#define reductionMultivariateData "tempREDMultiDFILE_recomb_mut.txt"
#define MixingFILE "mixing_recomb_mut.file"
#define tempFILE "temp_recomb_mut.file"
#define gnuplotFILE "gnuFILE_recomb_mut.txt"
#define gnuplotDDFILE "gnuDDFILE_recomb_mut.txt"
#define gnuplotDRSFILE "gnuDRSFILE_recomb_mut.txt"
#define simul_aver "simul_aver_recomb_mut.txt"
#define simul_all_aver "simul_all_aver_recomb_mut.txt"
#define simul_all_proc_right "simul_all_proc_right_recomb_mut.txt"
#define reductionData "gnuREDFILE_recomb_mut.txt"
#define simul_all_bestFitness "simul_all_bestFitness_recomb_mut.txt"
#define histogramFILE "histogram_recomb_mut.file"
#define fileHistogram "histogram_recomb_mut.data"
#define fileAcceptance "acceptance_recomb_mut.data"
#define diversityTempFILE "diverseTemp_recomb_mut.file"
#else
#ifdef NO_COUPLING
#define specific "no_coupling" 
#define CusumFILE "cusum_no_coupling.file"
#define reductionMultivariateData "tempREDMultiDFILE_no_coupling.txt"
#define MixingFILE "mixing_no_coupling.file"
#define tempFILE "temp_no_coupling.file"
#define gnuplotFILE "gnuFILE_no_coupling.txt"
#define gnuplotDDFILE "gnuDDFILE_no_coupling.txt"
#define gnuplotDRSFILE "gnuDRSFILE_no_coupling.txt"
#define simul_aver "simul_aver_no_coupling.txt"
#define simul_all_aver "simul_all_aver_no_coupling.txt"
#define simul_all_proc_right "simul_all_proc_right_no_coupling.txt"
#define reductionData "gnuREDFILE_no_coupling.txt"
#define simul_all_bestFitness "simul_all_bestFitness_no_coupling.txt"
#define histogramFILE "histogram_no_coupling.file"
#define fileHistogram "histogram_no_coupling.data"
#define fileAcceptance "acceptance_no_coupling.data"
#define diversityTempFILE "diverseTemp_no_coupling.file"
#else
#ifdef PARALLEL_TEMPERING
#define specific "par_temp" 
#define CusumFILE "cusum_par_temp.file"
#define reductionMultivariateData "tempREDMultiDFILE_par_temp.txt"
#define MixingFILE "mixing_par_temp.file"
#define tempFILE "temp_par_temp.file"
#define gnuplotFILE "gnuFILE_par_temp.txt"
#define gnuplotDDFILE "gnuDDFILE_par_temp.txt"
#define gnuplotDRSFILE "gnuDRSFILE_par_temp.txt"
#define simul_aver "simul_aver_par_temp.txt"
#define simul_all_aver "simul_all_aver_par_temp.txt"
#define simul_all_proc_right "simul_all_proc_right_par_temp.txt"
#define reductionData "gnuREDFILE_par_temp.txt"
#define simul_all_bestFitness "simul_all_bestFitness_par_temp.txt"
#define histogramFILE "histogram_par_temp.file"
#define fileHistogram "histogram_par_temp.data"
#define fileAcceptance "acceptance_par_temp.data"
#define diversityTempFILE "diverseTemp_par_temp.file"
#else
#ifdef ELITIST_ACCEPTANCE
#define specific "elit"
#define CusumFILE "cusum_elit.file"
#define reductionMultivariateData "tempREDMultiDFILE_elit.txt"
#define MixingFILE "mixing_elit.file"
#define tempFILE "temp_elit.file"
#define gnuplotFILE "gnuFILE_elit.txt"
#define gnuplotDDFILE "gnuDDFILE_elit.txt"
#define gnuplotDRSFILE "gnuDRSFILE_elit.txt"
#define simul_aver "simul_aver_elit.txt"
#define simul_all_aver "simul_all_aver_elit.txt"
#define simul_all_proc_right "simul_all_proc_right_elit.txt"
#define reductionData "gnuREDFILE_elit.txt"
#define simul_all_bestFitness "simul_all_bestFitness_elit.txt"
#else
#ifdef ELITIST_ACCEPTANCE_WITH_REGENERATION
#define specific "elit_regen"
#define CusumFILE "cusum_elit_regen.file"
#define reductionMultivariateData "tempREDMultiDFILE_elit_regen.txt"
#define MixingFILE "mixing_elit_regen.file"
#define tempFILE "temp_elit_regen.file"
#define gnuplotFILE "gnuFILE_elit_regen.txt"
#define gnuplotDDFILE "gnuDDFILE_elit_regen.txt"
#define gnuplotDRSFILE "gnuDRSFILE_elit_regen.txt"
#define simul_aver "simul_aver_elit_regen.txt"
#define simul_all_aver "simul_all_aver_elit_regen.txt"
#define simul_all_proc_right "simul_all_proc_right_elit_regen.txt"
#define reductionData "gnuREDFILE_elit_regen.txt"
#define simul_all_bestFitness "simul_all_bestFitness_elit_regen.txt"
#else
#ifdef RANDOM_COUPLING
#define specific "random_coup"
#define CusumFILE "cusum_rand_coup.file"
#define reductionMultivariateData "tempREDMultiDFILE_random_coup.txt"
#define MixingFILE "mixing_random_coup.file"
#define tempFILE "temp_random_coup.file"
#define gnuplotFILE "gnuFILE_random_coup.txt"
#define gnuplotDDFILE "gnuDDFILE_random_coup.txt"
#define gnuplotDRSFILE "gnuDRSFILE_random_couptxt"
#define simul_aver "simul_aver_random_coup.txt"
#define simul_all_aver "simul_all_aver_random_coup.txt"
#define simul_all_proc_right "simul_all_proc_right_random_coup.txt"
#define reductionData "gnuREDFILE_random_coup.txt"
#define simul_all_bestFitness "simul_all_bestFitness_random_coup.txt"
#else
#ifdef  PARALLEL_SIMULATED_ANNEALING
#define specific "sim_anneal"
#define CusumFILE "cusum_sim_anneal.file"
#define reductionMultivariateData "tempREDMultiDFILE_sim_anneal.txt"
#define MixingFILE "mixing_sim_anneal.file"
#define tempFILE "temp_sim_anneal.file"
#define gnuplotFILE "gnuFILE_sim_anneal.txt"
#define gnuplotDDFILE "gnuDDFILE_sim_anneal.txt"
#define gnuplotDRSFILE "gnuDRSFILE_sim_anneal.txt"
#define simul_aver "simul_aver_sim_anneal.txt"
#define simul_all_aver "simul_all_aver_sim_anneal.txt"
#define simul_all_proc_right "simul_all_proc_right_sim_anneal.txt"
#define reductionData "gnuREDFILE_sim_anneal.txt"
#define simul_all_bestFitness "simul_all_bestFitness_sim_anneal.txt"
#else
#ifdef  POPULATION_MCMC
#define specific "pop_MCMC"
#define CusumFILE "cusum_pop_MCMC.file"
#define reductionMultivariateData "tempREDMultiDFILE_pop_MCMC.txt"
#define MixingFILE "mixing_pop_MCMC.file"
#define tempFILE "temp_pop_MCMC.file"
#define gnuplotFILE "gnuFILE_pop_MCMC.txt"
#define gnuplotDDFILE "gnuDDFILE_pop_MCMC.txt"
#define gnuplotDRSFILE "gnuDRSFILE_pop_MCMC.txt"
#define simul_aver "simul_aver_pop_MCMC.txt"
#define simul_all_aver "simul_all_aver_pop_MCMC.txt"
#define simul_all_proc_right "simul_all_proc_right_pop_MCMC.txt"
#define reductionData "gnuREDFILE_pop_MCMC.txt"
#define simul_all_bestFitness "simul_all_bestFitness_pop_MCMC.txt"
#endif
#endif
#endif
#endif 
#endif
#endif
#endif
#endif

//#define diversityTempFILE "diverseTemp.file"
#define distributionFILE "distributionFILE.file"

#define NORMAL_RUNS

#define SAMPLE_FILE
#define sampleFile "sample.file"

//#define MULTIPLE_RUNS
#ifdef MULTIPLE_RUNS
#define MUTATION_VARIES
//#ifdef TEMPERATURE_CHANGE
#define SIZEGENOME_CT 1
#define TEMPERATURE_CT 1
//#define MUTATION_VARIES

#endif
