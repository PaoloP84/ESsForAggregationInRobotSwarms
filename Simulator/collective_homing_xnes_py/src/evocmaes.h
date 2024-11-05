/*
 * evo.h
 */

#include <QVector>

extern long int maxevaluations;		// max number of evaluations
extern int seed;					// random seed
extern int nreplications;			// number of replications
extern double prange;				// parameters range
extern double *bestggenotype;		// vector of genes
extern int steptimelength;			// step timelength during display
extern char filedir[512];
extern int intNeurons; 			// The number of internal neurons (hidden + output)

extern void initES();
extern void finalizeES();
extern double evaluate(double* genotype, const int glen, int mode, int seed);
extern void createInitialStates(int mode);
extern void saveGenotype(char* filename, double* genotype, int glen);
extern int loadGenotype(char* filename, double* genotype, int glen);
extern void savePop(char* filename, double **pop, int glen, int popsize);
extern void loadPop(char* filename, double **pop, int glen, int popsize);
extern bool readEvoParameters(const char* filename);
extern void evolve(int ngenes, int cseed);
extern void setInternalNeurons(int nneurons);
