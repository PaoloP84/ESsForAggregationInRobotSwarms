/*
 * xNES evolutionary algorithm
 * it require the Eigen and unsupported librieris (costituted only by include files)
 * and the GSL libraries (included in the utilities.h file) for the generation of random numbers
 * it call an eternal evalaate function that is used to evaluate candidate solutions
 */
#include <Python.h>
#include "numpy/arrayobject.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <time.h>
#include <typeinfo>
#include "utilities.h"
#include "evoopenaies.h"

#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./unsupported/Eigen/MatrixFunctions"

//#include "gsl/gsl_rng.h"
//#include "gsl/gsl_randist.h"

#define NUM_PY_FUNCTS 6

int trials=1;
int ttrials=1;
double *bestggenotype;				// the best generalizing genotype
long int maxevaluations = 100;		// max number of evaluations
int seed = 1;						// random seed
int nreplications = 1;				// number of replications
double prange = 1.0;				// parameters range
int steptimelength = 0;				// time length of a step during graphic visualization
int envchangeevery = 20;			// environment changed every n generations
int batchSize = 20;				// Number of samples
bool wdecay = true;
RandomGenerator* RNG;               // randomgenerator class
char filedir[512];
int intNeurons;

#if PY_MAJOR_VERSION >= 3
int
#else
void
#endif
init_numpy()
{
        import_array();
}

/* Python section of the code (to manage python objects and use python wrappers) */

// Python objects
PyObject *pName, *pModule;
PyObject *pFunc1, *pFunc2, *pFunc3, *pFunc4, *pFunc5, *pFunc6;
PyObject *pValue, *pArgs, *pMat, *pVec, *pFirst, *pLast;

void clearPythonObjects()
{
	Py_DECREF(pFunc1);
	Py_DECREF(pFunc2);	
	Py_DECREF(pFunc3);
	Py_DECREF(pFunc4);	
	Py_DECREF(pFunc5);
	Py_DECREF(pFunc6);
}

void fillPythonObjects()
{
	// Here we set numpy as the library to be imported
	// since we know in advance the methods we need to
	// access to... We use a wrapper file of numpy
	const char* pyfile = "numpywrapper";
	const char* pyfuncts[NUM_PY_FUNCTS] = {"seedwrapper", "randomwrapper", "randnwrapper", "normwrapper", "eigwrapper", "dotwrapper"};
	pName = PyString_FromString(pyfile);
	pModule = PyImport_Import(pName);
	Py_DECREF(pName);
	if (pModule == NULL)
	{
		PyErr_Print();
        	fprintf(stderr, "Failed to load \"%s\"\n", pyfile);
		printf("Stopped...\n");
        	exit(-2);
	}	
	/* Get pointers to python functions */
	pFunc1 = PyObject_GetAttrString(pModule, pyfuncts[0]);
	if (pFunc1 == NULL)
	{
		printf("Failed to load %s from numpy/%s... Abort!\n", pyfuncts[0], pyfile);
		exit(-2);
	}
	pFunc2 = PyObject_GetAttrString(pModule, pyfuncts[1]);
	if (pFunc2 == NULL)
	{
		printf("Failed to load %s from numpy/%s... Abort!\n", pyfuncts[1], pyfile);
		exit(-2);
	}
	pFunc3 = PyObject_GetAttrString(pModule, pyfuncts[2]);
	if (pFunc3 == NULL)
	{
		printf("Failed to load %s from numpy/%s... Abort!\n", pyfuncts[2], pyfile);
		exit(-2);
	}
	pFunc4 = PyObject_GetAttrString(pModule, pyfuncts[3]);
	if (pFunc4 == NULL)
	{
		printf("Failed to load %s from numpy/%s... Abort!\n", pyfuncts[3], pyfile);
		exit(-2);
	}
	pFunc5 = PyObject_GetAttrString(pModule, pyfuncts[4]);
	if (pFunc5 == NULL)
	{
		printf("Failed to load %s from numpy/%s... Abort!\n", pyfuncts[4], pyfile);
		exit(-2);
	}
	pFunc6 = PyObject_GetAttrString(pModule, pyfuncts[5]);
	if (pFunc6 == NULL)
	{
		printf("Failed to load %s from numpy/%s... Abort!\n", pyfuncts[5], pyfile);
		exit(-2);
	}
	Py_DECREF(pModule);
	/* Check whether or not all the functions are callable from C++ */
	if (!PyCallable_Check(pFunc1))
	{
		printf("Function %s not callable... Abort!\n", pyfuncts[0]);
		exit(-1);
	}
	if (!PyCallable_Check(pFunc2))
	{
		printf("Function %s not callable... Abort!\n", pyfuncts[1]);
		exit(-1);
	}
	if (!PyCallable_Check(pFunc3))
	{
		printf("Function %s not callable... Abort!\n", pyfuncts[2]);
		exit(-1);
	}
	if (!PyCallable_Check(pFunc4))
	{
		printf("Function %s not callable... Abort!\n", pyfuncts[3]);
		exit(-1);
	}
	if (!PyCallable_Check(pFunc5))
	{
		printf("Function %s not callable... Abort!\n", pyfuncts[4]);
		exit(-1);
	}
	if (!PyCallable_Check(pFunc6))
	{
		printf("Function %s not callable... Abort!\n", pyfuncts[5]);
		exit(-1);
	}
}
/* End Python section */

void initES()
{
	// Initialize python
    	Py_Initialize();
    	PyObject *sys = PyImport_ImportModule("sys");
    	PyObject *path = PyObject_GetAttrString(sys, "path");
    	PyList_Append(path, PyString_FromString("."));
    	if(PyArray_API == NULL)
    	{
		//import_array();
		init_numpy();
    	}
    	fillPythonObjects();
}

void finalizeES()
{
	clearPythonObjects();
    	// Finalize python
    	Py_Finalize();
}

/*
 * save a genotype in a file
 */
void saveGenotype(char* filename, double* genotype, int glen)
{
	FILE *fp;
	
	fp = fopen(filename,"wb!");
	if (!fp)
      		printf("ERROR: unable to open file %s\n", filename);
     	else
      	{
      		fwrite(&glen, sizeof(int), 1, fp);
	  	fwrite(genotype, sizeof(double)*glen, 1, fp);
		fclose(fp);
      	}
}

/*
 * load a genotype from a file
 */
int loadGenotype(char* filename, double* genotype, int glen)
{
	FILE *fp;
	size_t result;
    	int np;
	
	fp = fopen(filename,"rb");
	if (!fp)
	{
	  	printf("ERROR: unable to open file %s\n", filename);
	  	return(0);
	}
     	else
	{
	      	fread(&np, sizeof(int), 1, fp);
	      	if (np != glen)
		{
		 	printf("ERROR: the number of parameters contained in the file %s differ from that expected\n", filename);
		 	//if (!feof(fp) || result != sizeof(double)*glen */)
		 	//  printf("ERROR: the number of parameters contained in the file %s differ from what is expected\n", filename);
		}
	       	else
		{
		 	result = fread(genotype, sizeof(double)*glen, 1, fp);
		}
	      	fclose(fp);
	}
    	return(1);
}

/*
 * save a genotype in a file
 */
void savePop(char* filename, double **pop, int glen, int popsize)
{
	FILE *fp;
	int p;
	
	fp = fopen(filename,"wb!");
	if (!fp)
      		printf("ERROR: unable to open file %s\n", filename);
     	else
      	{
	  	fwrite(&popsize, sizeof(int), 1, fp);
      		fwrite(&glen, sizeof(int), 1, fp);
	  	for(p=0; p < popsize; p++, pop++)
	    		fwrite(*pop, sizeof(double)*glen, 1, fp);
		fclose(fp);
      	}
}

/*
 * load a genotype from a file
 */
void loadPop(char* filename, double **pop, int glen, int popsize)
{
	FILE *fp;
    	int np;
	int ps;
	int p;
	
	fp = fopen(filename,"rb");
	if (!fp)
	{
	  	printf("ERROR: unable to open file %s\n", filename);
	}
     	else
	{
	  	fread(&ps, sizeof(int), 1, fp);
      		fread(&np, sizeof(int), 1, fp);
      		if (np != glen || ps != popsize)
        	{
         		printf("ERROR: the number of parameters contained in the file %s differ from that expected\n", filename);
        	}
       		else
        	{
		 	for(p=0; p < popsize; p++, pop++)
           			fread(*pop, sizeof(double)*glen, 1, fp);
        	}
      		fclose(fp);
	}
}

void setInternalNeurons(int nneurons)
{
	intNeurons = nneurons;
}

/* Sort samples based on fitness. Returns sorted fitnesses and corresponding
   indices. Flag <maximize> specifies the type of sorting (i.e., maximization or minimization) */
static void sortSamples(const double* f, const int size, double* sf, int* si, bool maximize)
{
	int i;
	double bf;
	int bi;
	int coeff;
	double* tmpf = new double[size];
	if (maximize)
		coeff = -1;
	else
		coeff = 1;
	if (sf == NULL)
		sf = new double[size];
	if (si == NULL)
		si = new int[size];
	for (i = 0; i < size; i++)
		tmpf[i] = f[i];
	i = 0;
	while (i < size)
	{
		bf = tmpf[0];
		bi = 0;
		for (int j = 1; j < size; j++)
		{
			if (maximize)
			{
				// Maximization
				if (tmpf[j] > bf)
				{
					bf = tmpf[j];
					bi = j;
				}
			}
			else
			{
				// Minimization
				if (tmpf[j] < bf)
				{
					bf = tmpf[j];
					bi = j;
				}
			}
		}
		sf[i] = bf;
		si[i] = bi;
		tmpf[bi] = coeff * 999999.0;
		i++;
	}
	delete tmpf;
}

//! OpenAI-ES evolution algorithm
void evolve(int _glen, int cseed)
{
	std::vector<double*> _individuals;
	_individuals.resize(1);
	for (int i = 0; i < 1; i++)
	{
		_individuals[i] = new double[_glen];
	}
	Eigen::VectorXd ind(_glen); // The individual from which the offspring will be generated
	Eigen::VectorXd dInd(_glen); // The individual from which the offspring will be generated
	Eigen::VectorXd grad(_glen); // Gradient
	Eigen::VectorXd globalg(_glen); // Actual gradient!!!
	Eigen::VectorXd utilities(batchSize * 2); // Array containing the utilities
	Eigen::VectorXd weights(batchSize); // Array containing the weights
	double* f = new double[batchSize * 2]; // Array containing the fitness of the offspring
	double* ordf = new double[batchSize * 2]; // Array containing the ordered fitness of the offspring
	int* ordoff = new int[batchSize * 2]; // Array containing the ordered indices of the offspring
	double* currInd = new double[_glen]; // Current individual to be evaluated (temporary variable)
	char filen[256];						// filename
	int ceval;								// current evaluation
	int gn;									// current generation
	double* bestind = new double[_glen];	// genotype of the best individual to date
	double* bestgsol = new double[_glen];	// genotype of the best generalizing individual to date
	int i;
	int j;
	int k;
    int s, ss;
	double gfit;							// generalization fitness
	double bestfit;							// best fitness to date
	double bestgfit;						// best generalization fitness to date
	double fitcentroid;						// fitness of the centroid
    double *statistics;                     // performance across generations
    double *stat;                           // performance across generations, local pointer
    long int statn;                         // number of stat data
    FILE *fp;
	char filename[512];
	// Previous individual
	Eigen::VectorXd prevInd(_glen);
	double distInd;
	std::vector<double> dists;
	double noiseStdDev = 0.02;
	// Adam optimizer parameters
	Eigen::VectorXd expAvg(_glen);
	Eigen::VectorXd expAvgSq(_glen);
	double stepSize = 0.01;
	double beta1 = 0.9;
	double beta2 = 0.999;
	double epsilon = pow(10.0,-8.0);

	// Overwrite seed
	seed = cseed;

    printf("OpenAI-ES seed %d MaxEval %ld prange %.1f\n", seed, maxevaluations, prange);
	
    time_t startTimer;
    time_t endTimer;
    time(&startTimer);
    // set the seed and generate the initial states
    RNG->setSeed(seed);
	// Set seed
    // Here we call numpy.random.seed() function
    pArgs = PyTuple_New(1);
	pValue = PyInt_FromLong(seed);
	if (!pValue)
	{
		Py_DECREF(pArgs);
		Py_DECREF(pModule);
		fprintf(stderr, "Cannot convert argument\n");
		exit(-2);
	}
	PyTuple_SetItem(pArgs, 0, pValue);
	pValue = PyObject_CallObject(pFunc1, pArgs);
	Py_DECREF(pArgs);
	if (pValue != NULL) 
	{
		Py_DECREF(pValue);
	}
	else
	{
		Py_DECREF(pFunc1);
		Py_DECREF(pModule);
		PyErr_Print();
		fprintf(stderr,"Call failed\n");
		exit(-2);
	}
    createInitialStates(0);

    // allocate memory for statistics
    statn = (maxevaluations / (trials * batchSize * 2) * 7) + 100;
    statistics = (double *) malloc(statn * sizeof(double));
    stat = statistics;
    statn = 0;

	// Initialize matrices and vectors
	for (i = 0; i < _glen; i++)
	{
		expAvg[i] = 0.0; // Initialize individual
		expAvgSq[i] = 0.0; // Initialize individual
	}
	printf("Seed %d \n", seed);
	fflush(stdout);
	
		bestfit = -999999999.0;
		bestgfit = -999999999.0;
		// Initialize the individual for this replication
		for (i = 0; i < _glen; i++)
        {
            //_individuals[0][i] = - prange + ( double( rand() ) / double( RAND_MAX ) ) * (prange * 2.0);
            //_individuals[0][i] = RNG->getDouble(-prange, prange);
			pArgs = NULL;
		    pValue = PyObject_CallObject(pFunc2, pArgs);
			_individuals[0][i] = -prange + (PyFloat_AsDouble(pValue) * prange * 2.0);
			Py_DECREF(pValue);
			// Reset the current offspring (to avoid numerical issues)
			currInd[i] = 0.0;
		}
		// Initialize fitnesses, ordered offspring and ordered utilities
		for (j = 0; j < batchSize * 2; j++) // Remember that the last slot is reserved for the parent
		{
			f[j] = 0.0;
			ordf[j] = 0.0;
			ordoff[j] = 0.0;
		}

		// Store initial (random) genotype
		strcpy(filename, filedir);
        sprintf(filen,"initCentroidS%d.gen", seed);
        strcat(filename, filen);
        fp = fopen(filename,"w!");
        for(i = 0; i < _glen; i++)
           	fprintf(fp, "%lf ", _individuals[0][i]);
        fclose(fp);

		// While (maxevaluations) loop
		gn = 0;
		ceval = 0;
		while (ceval < maxevaluations)
		{
			if (envchangeevery > 0)
			{
				// change environment every n generations and re-evaluate parents
            			if (gn > 0 && (gn % envchangeevery) == 0)
			     	{
			      		createInitialStates(1); // Change environment only for evolution trials (test trials are fixed!!!)
				  	printf("Environment changed\n");
			     	}
			}
			
			//printf("Gen %d) ", (gn + 1));

			for (i = 0; i < _glen; i++)
			{
				// Copy back individual from previous generation
				ind[i] = _individuals[0][i];
				prevInd[i] = _individuals[0][i];
			}
			
			// Extracting samples (i.e. generate offspring of current individual)
			Eigen::MatrixXd samples(_glen, batchSize);
			Eigen::MatrixXd offspring(_glen, batchSize * 2);
			// Here we call numpy.random.randn() function
			pArgs = PyTuple_New(2);
			pValue = PyInt_FromLong(_glen);
			if (!pValue)
			{
				Py_DECREF(pArgs);
				Py_DECREF(pModule);
				fprintf(stderr, "Cannot convert argument\n");
				exit(-2);
		        }
		   	PyTuple_SetItem(pArgs, 0, pValue);
			pValue = PyInt_FromLong(batchSize);
			if (!pValue)
			{
				Py_DECREF(pArgs);
				Py_DECREF(pModule);
				fprintf(stderr, "Cannot convert argument\n");
				exit(-2);
		        }
		   	PyTuple_SetItem(pArgs, 1, pValue);
			pValue = PyObject_CallObject(pFunc3, pArgs);
			Py_DECREF(pArgs);
			if (pValue != NULL) 
			{
				double* mat = (double*)PyArray_DATA(pValue);
				for (j = 0; j < batchSize; j++)
				{
					for (i = 0; i < _glen; i++)
					{
						samples(i,j) = mat[i + (j * _glen)];
					}
				}
				Py_DECREF(pValue);
			}
			else
			{
				Py_DECREF(pFunc3);
				Py_DECREF(pModule);
				PyErr_Print();
				fprintf(stderr,"Call failed\n");
				exit(-2);
		        }
			// Draw out samples from a Gaussian distribution 
			/*for (i = 0; i < _glen; i++)
			{
				for (j = 0; j < batchSize; j++)
				{
                    			samples(i,j) = RNG->getGaussian(1.0, 0.0);
				}
			}*/
			// Create the offspring
			for (i = 0; i < _glen; i++)
			{
				for (j = 0; j < batchSize; j++)
				{
					// Symmetric samples
					for (k = 0; k < 2; k++)
					{
						if (k == 0)
							offspring(i, 2 * j) = _individuals[0][i] + noiseStdDev * samples(i,j);
						else
							offspring(i, 2 * j + 1) = _individuals[0][i] - noiseStdDev * samples(i,j);
					}
				}
			}
			// Fill and evaluate the current offspring
			for (j = 0; j < batchSize * 2; j++)
			{
				for (i = 0; i < _glen; i++)
				 {
				   currInd[i] = offspring(i,j);
				 }
				f[j] = evaluate(currInd, _glen, 0, -1);
				ceval += trials;
			}
			// Sort samples based on fitness
			sortSamples(f, batchSize * 2, ordf, ordoff, false); // We need ascendent sorting!!!
			// Compute utilities
			for (j = 0; j < batchSize * 2; j++)
			{
				utilities[ordoff[j]] = j;
			}
			utilities /= ((batchSize * 2) - 1);
			for (j = 0; j < batchSize * 2; j++)
			{
				utilities[j] -= 0.5;
			}
			// Compute weights
			for (j = 0; j < batchSize; j++)
			{
				weights[j] = (utilities[2 * j] - utilities[2 * j + 1]);
			}
			// Compute gradient
			/*grad = samples * weights;
			// Normalize for the number of individuals that have been evaluated
			grad /= (batchSize * 2);*/
			//printf("\n");
			double* v = new double[batchSize];
			for (j = 0; j < batchSize; j++)
				v[j] = weights[j];
			long int vdim[] = { batchSize };
			pArgs = PyTuple_New(2);
			pVec = PyArray_SimpleNewFromData( 1, vdim, PyArray_DOUBLE, v );
			if (!pVec)
			{
				Py_DECREF(pArgs);
				fprintf(stderr, "Cannot convert argument\n");
				exit(-2);
		    	}
		   	/* pValue reference stolen here: */
		    	PyTuple_SetItem(pArgs, 0, pVec);
			double** m = new double*[batchSize];
		    	m[0] = new double[batchSize * _glen];
			for (j = 1; j < batchSize; j++)
		    		m[j] = m[j - 1] + _glen;				
		    	// fill in values				
		    	for (j = 0; j < batchSize; j++)
			{
				for (i = 0; i < _glen; i++)
				{
					m[j][i] = samples(i,j);
				}
			}	
		    	long int mdim[] = { batchSize, _glen };
				// Here we call numpy.linalg.eig() function
			pMat = PyArray_SimpleNewFromData( 2, mdim, PyArray_DOUBLE, m[0] );
			if (!pMat)
			{
				Py_DECREF(pArgs);
				fprintf(stderr, "Cannot convert argument\n");
				exit(-2);
		    	}
		    	PyTuple_SetItem(pArgs, 1, pMat);
		    	pValue = PyObject_CallObject(pFunc6, pArgs);
			Py_DECREF(pArgs);
		    	if (pValue != NULL) 
			{
				double* vec = (double*)PyArray_DATA(pValue);
				for (i = 0; i < _glen; i++)
					grad[i] = vec[i];
				Py_DECREF(pValue);
		    	}
		    	else
			{
				Py_DECREF(pFunc6);
				PyErr_Print();
				fprintf(stderr,"Call failed\n");
				exit(-2);
		    	}
			grad /= (batchSize * 2);
			if (wdecay)
				globalg = -grad + 0.005 * ind;
			else
				globalg = -grad;
			// Adam optimizer
			double a = stepSize * sqrt(1.0 - pow(beta2, gn + 1)) / (1.0 - pow(beta1, gn + 1));
			for (i = 0; i < _glen; i++)
			{
				expAvg[i] = beta1 * expAvg[i] + (1.0 - beta1) * globalg[i];
				expAvgSq[i] = beta2 * expAvgSq[i] + (1.0 - beta2) * globalg[i] * globalg[i];
				dInd[i] = -a * expAvg[i] / (sqrt(expAvgSq[i] + epsilon));
				ind[i] += dInd[i];
			}
			// Copy back the individual
			for (i = 0; i < _glen; i++)
			{
				_individuals[0][i] = ind[i];
			}
			// Evaluate the individual (to find out its fitness)
			fitcentroid = evaluate(_individuals[0], _glen, 0, -1);
			ceval += trials;
			// the best of the current generation is the centroid
			if (fitcentroid > ordf[batchSize * 2 - 1])
			{
				if (fitcentroid > bestfit)
				{
					bestfit = fitcentroid;
					for (i = 0; i < _glen; i++)
					{
						bestind[i] = _individuals[0][i];
					}
				}
				// Generalization
				gfit = evaluate(bestind, _glen, 1, -1);
				ceval += ttrials;
				if (gfit > bestgfit)
				{
					bestgfit = gfit;
					for (i = 0; i < _glen; i++)
					{
						bestgsol[i] = _individuals[0][i];
					}
				}
			}
			else
			{
			    // the best of the current generation is one of the offspring
				if (ordf[batchSize * 2 - 1] > bestfit)
				{
					bestfit = ordf[batchSize * 2 - 1];
					for (i = 0; i < _glen; i++)
					{
						bestind[i] = offspring(i,ordoff[batchSize * 2 - 1]);
					}
				}
				// Generalization
				for (i = 0; i < _glen; i++)
				{
					currInd[i] = offspring(i,ordoff[batchSize * 2 - 1]);
				}
				gfit = evaluate(currInd, _glen, 1, -1);
				ceval += ttrials;
				if (gfit > bestgfit)
				{
					bestgfit = gfit;
					for (i = 0; i < _glen; i++)
					{
						bestgsol[i] = offspring(i,ordoff[batchSize * 2 - 1]);
					}
				}
			}
			// Compute distance between previous and current individual
			distInd = 0.0;
			for (i = 0; i < _glen; i++)
			{
				distInd += fabs(_individuals[0][i] - prevInd[i]);
			}
			distInd /= _glen;
			dists.push_back(distInd);
			// Reset weights and fitnesses once used
			for (j = 0; j < batchSize * 2; j++)
			{
				f[j] = 0.0;
			}
			// Reset current offspring
			for (i = 0; i < _glen; i++)
			{
				currInd[i] = 0.0;
			}
            printf("Gen %d Steps %d Centroid %.2f BestOffspring %.2f Bestfit %.2f BestgFit %.2f\n", gn, ceval, fitcentroid, ordf[batchSize * 2 - 1], bestfit, bestgfit);
			fflush(stdout);
			gn++;

            // store stat data
            *stat = ceval; stat++;
            *stat = bestgfit; stat++;
            *stat = bestfit; stat++;
            *stat = ordf[batchSize * 2 - 1]; stat++;
            *stat = fitcentroid; stat++;
            statn++;

			if ((gn % 100) == 0)
			{
				strcpy(filename, filedir);
        		sprintf(filen,"bestgS%d.bgen", seed);
				strcat(filename, filen);
        		saveGenotype(filename, bestgsol, _glen);
			}
			free(v);
			free(m[0]);
			free(m);
        }

	strcpy(filename, filedir);
        sprintf(filen,"bestgS%d.bgen", seed);
	strcat(filename, filen);
        saveGenotype(filename, bestgsol, _glen);
	
	double *gen;
	strcpy(filename, filedir);
        sprintf(filen,"bestgS%d.gen", seed);
        strcat(filename, filen);
        fp = fopen(filename,"w!");
        for(i = 0, gen = bestgsol; i < _glen; i++, gen++)
           	fprintf(fp, "%lf ", *gen);
        fclose(fp);

	strcpy(filename, filedir);
        sprintf(filen, "S%d.fit", seed);
	strcat(filename, filen);
        if  ((fp = fopen(filename,"w")) != NULL)
        {
               fprintf(fp, "Gen %d Centroid %.2f BestOffspring %.2f Bestfit %.2f BestgFit %.2f\n", gn, fitcentroid, ordf[0], bestfit, bestgfit);
               fclose(fp);
        }

	strcpy(filename, filedir);
        sprintf(filen, "statS%d.fit", seed);
        strcat(filename, filen);
        if  ((fp = fopen(filename,"w")) != NULL)
        {
               fprintf(fp, "cevaluation bestgfit bestfit fitbestoffspring fitcentroid \n");
               for(s = 0, stat = statistics; s < statn; s++)
               {
                  	for(ss=0; ss < 5; ss++, stat++)
                    		fprintf(fp, "%lf ", *stat);
                  	fprintf(fp, "\n");
               }
               fclose(fp);
        }

	strcpy(filename, filedir);
        sprintf(filen, "distS%d.txt", seed);
        strcat(filename, filen);
        if  ((fp = fopen(filename,"w")) != NULL)
        {
               for (i = 0; i < dists.size(); i++)
               {
                  	fprintf(fp, "%lf\n", dists[i]);
               }
               fclose(fp);
        }

	time(&endTimer);
	double seconds = difftime(endTimer,startTimer);
	printf("evolve() took %lf seconds\n", seconds);
}

/*
 * read parameters from the configuration file
 */
bool readEvoParameters(const char* filename)
{
    	char *s;
    	char buff[1024];
    	char name[1024];
    	char value[1024];
    	char *ptr;
    	int section;  // 0=before the section 1=in the section 2= after the section

    	section = 0;
    	FILE* fp = fopen(filename, "r");
    	if (fp != NULL)
    	{
        	// Read lines
        	while (fgets(buff, 1024, fp) != NULL)
        	{
            		// check whether the relevant section start or end
            		if (buff[0] == '[')
            		{
              			if (section == 1)
                		{
                    			section = 2;
                    			continue;
                		}
              			if ((section == 0) && (strncmp(buff, "[EVO]",5)==0))
              			{
                			section = 1;
                			continue;
              			}
            		}

            		if (section == 1)
            		{
            			//Skip blank lines and comments
            			if (buff[0] == '\n' || buff[0] == '#' || buff[0] == '/')
            				continue;

            			//Parse name/value pair from line
            			s = strtok(buff, " = ");
            			if (s == NULL)
            				continue;
            			else
            				copyandclear(s, name);

            			s = strtok(NULL, " = ");
            			if (s == NULL)
            				continue;
            			else
            				copyandclear(s, value);

            			// Copy into correct entry in parameters struct
            			if (strcmp(name, "maxevaluations")==0)
            				maxevaluations = (int)strtol(value, &ptr, 10);
            			else if (strcmp(name, "nreplications")==0)
            				nreplications = (int)strtol(value, &ptr, 10);
            			else if (strcmp(name, "prange")==0)
            				prange = strtod(value, &ptr);
            			else if (strcmp(name, "seed")==0)
            				seed = (int)strtol(value, &ptr, 10);
            			else if (strcmp(name, "steptimelength")==0)
            				steptimelength = (int)strtol(value, &ptr, 10);
            			else if (strcmp(name, "envchangeevery")==0)
            				envchangeevery = (int)strtol(value, &ptr, 10);
				else if (strcmp(name, "batchSize")==0)
					batchSize = (int)strtol(value, &ptr, 10);
				else if (strcmp(name, "wdecay")==0)
				{
					int wd = (int)strtol(value, &ptr, 10);
					wdecay = (bool)wd;
				}
				else printf("WARNING: Unknown parameter %s in section [EVO] of file %s \n", name, filename);
            		}
         	}
        	fclose (fp);
        	if (section == 0)
           		printf("WARNING: Missing section [EVO] in file %s \n", filename);
        	return(true);
    	}
    	else
    	{
        	printf("ERROR: unable to open file %s\n", filename);
        	return(false);
    	}
}


void test(int _glen)
{

}


