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
#include "evocmaes.h"

#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./unsupported/Eigen/MatrixFunctions"

//#include "gsl/gsl_rng.h"
//#include "gsl/gsl_randist.h"

int trials=1;
int ttrials=1;
double *bestggenotype;				// the best generalizing genotype
long int maxevaluations = 100;		// max number of evaluations
int seed = 1;						// random seed
int nreplications = 1;				// number of replications
double prange = 1.0;				// parameters range
int steptimelength = 0;				// time length of a step during graphic visualization
int envchangeevery = 20;			// environment changed every n generations
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
PyObject *pFunc1, *pFunc2, *pFunc3, *pFunc4, *pFunc5;
PyObject *pValue, *pArgs, *pMat, *pVec, *pFirst, *pLast;

void clearPythonObjects()
{
	Py_DECREF(pFunc1);
	Py_DECREF(pFunc2);	
	Py_DECREF(pFunc3);
	Py_DECREF(pFunc4);	
	Py_DECREF(pFunc5);
}

void fillPythonObjects()
{
	// Here we set numpy as the library to be imported
	// since we know in advance the methods we need to
	// access to... We use a wrapper file of numpy
	const char* pyfile = "numpywrapper";
	const char* pyfuncts[5] = {"seedwrapper", "randomwrapper", "randnwrapper", "normwrapper", "eigwrapper"};
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

/* Convert an array to a diagonal matrix */
static Eigen::MatrixXd diag(Eigen::VectorXd v, const int size)
{
	Eigen::MatrixXd m(size,size);
	int i;
	int j;
	for (i = 0; i < size; i++)
	{
		m(i,i) = v[i];
		for (j = 0; j < size; j++)
		{
			if (i != j)
			{
				m(i,j) = 0.0;
			}
		}
	}
	return m;
}

//! CMA-ES evolution algorithm
void evolve(int _glen, int cseed)
{

	std::vector<double*> _individuals;
	_individuals.resize(1);
	for (int i = 0; i < 1; i++)
	{
		_individuals[i] = new double[_glen];
	}

	double _stepRate;						// learning rate (usually 1.0)
	double genLen = (double)_glen;
	double sampleSize = 4.0 + floor(3.0 * log(genLen));
    int ssize = (int)sampleSize;            // number of samples
	int mu = (int)floor(sampleSize / 2.0);
	Eigen::VectorXd ind(_glen); // The individual from which the offspring will be generated
	Eigen::VectorXd covPath(_glen);
	Eigen::VectorXd stepPath(_glen);
	Eigen::MatrixXd C(_glen,_glen); // Covariance matrix
	Eigen::MatrixXd B(_glen,_glen); // Eigenvectors
	Eigen::MatrixXd D(_glen,_glen); // Eigenvalues
	double muEff;
	double cc;
	double csigma;
	double c1;
	double cmu;
	double covLearnRate;
	double dsigma;
	double chiN;
	double sigma;
	double tmp;
	Eigen::VectorXd weights(mu); // Array containing the weights
	double wSum; // Sum of weights
	double wSqSum; // Sum of squared weights
	double* f = new double[ssize]; // Array containing the fitness of the offspring
	double* ordf = new double[ssize]; // Array containing the ordered fitness of the offspring
	int* ordoff = new int[ssize]; // Array containing the ordered indices of the offspring
	double* currInd = new double[_glen]; // Current individual to be evaluated (temporary variable)
	double w; // Weight
	char filen[256];						// filename
	int ceval;								// current evaluation
	int gn;									// current generation
	double* bestind = new double[_glen];	// genotype of the best individual to date
	double* bestgsol = new double[_glen];	// genotype of the best generalizing individual to date
	int i;
	int j;
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

	// Overwrite seed
	seed = cseed;

    printf("CMA-ES seed %d MaxEval %ld prange %.1f\n", seed, maxevaluations, prange);
	
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
    statn = (maxevaluations / (trials * ssize) * 7) + 100;
    statistics = (double *) malloc(statn * sizeof(double));
    stat = statistics;
    statn = 0;

	// Calculate the utility values before normalization
	wSum = 0.0;
	for (j = 0; j < mu; j++)
	{
		// w[j] is equal to the uttermost between 0 and (log(lambda / 2 + 1) - log(j))
		tmp = log(mu + 1.0) - log((double)j + 1.0); // Remember that log(0) = -Inf, so 1 is added to index i
		if (tmp < 0.0)
		{
			tmp = 0.0;
		}
		weights[j] = tmp;
		// Update the sum of utility values
		wSum += weights[j];
	}
	// Normalize the utility values
	for (j = 0; j < mu; j++)
	{
		if (wSum != 0.0)
		{
			weights[j] /= wSum;
		}
	}
	// Compute new sum of weights and sum of w^2
	wSum = 0.0;
	wSqSum = 0.0;
	for (j = 0; j < mu; j++)
	{
		wSum += weights[j];
		wSqSum += (weights[j] * weights[j]);
	}
	muEff = (wSum * wSum) / wSqSum;
	cc = 4.0 / (genLen + 4.0);
	csigma = (muEff + 2.0) / (genLen + muEff + 3.0);
	c1 = ((1.0 / muEff) * 2.0) / ((genLen + 1.4) * (genLen + 1.4));
	cmu = ((1.0 - (1.0 / muEff)) * ((2.0 * muEff - 1.0) / ((genLen + 2.0) * (genLen + 2.0) + 2.0 * muEff)));
	covLearnRate = (c1 + cmu);
	dsigma = 1.0 + csigma;
	tmp = sqrt((muEff - 1.0) / (genLen + 1.0)) - 1.0;
	if (tmp > 0.0)
		dsigma += 2.0 * tmp;
	chiN = sqrt(genLen) * (1.0 - 1.0 / (4.0 * genLen) + 1.0 / (21.0 * (genLen * genLen)));
	// Initialize matrices and vectors
	for (i = 0; i < _glen; i++)
	{
		ind[i] = 0.0; // Initialize individual
		B(i,i) = 1.0;
		D(i,i) = 1.0;
		C(i,i) = 1.0;
		for (j = 0; j < _glen; j++)
		{
			if (i != j)
			{
				B(i,j) = 0.0;
				D(i,j) = 0.0;
				C(i,j) = 0.0;
			}
		}
		stepPath[i] = 0.0;
		covPath[i] = 0.0;
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
	    // Here we call numpy.random.random() function
			pArgs = NULL;
		    pValue = PyObject_CallObject(pFunc2, pArgs);
			_individuals[0][i] = -prange + (PyFloat_AsDouble(pValue) * prange * 2.0);
			Py_DECREF(pValue);
			// Reset the current offspring (to avoid numerical issues)
			currInd[i] = 0.0;
		}
		// Initialize fitnesses, ordered offspring and ordered utilities
		for (j = 0; j < ssize; j++) // Remember that the last slot is reserved for the parent
		{
			f[j] = 0.0;
			ordf[j] = 0.0;
			ordoff[j] = 0.0;
		}
		sigma = 0.5;

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
			Eigen::MatrixXd copyCenters(_glen,ssize);
			Eigen::MatrixXd samples(_glen,ssize);
			Eigen::MatrixXd offspring(_glen,ssize);
			Eigen::MatrixXd sortedSamples(_glen,ssize);
			Eigen::MatrixXd sortedOffspring(_glen,ssize);
			Eigen::MatrixXd selectedSamples(_glen,mu);
			Eigen::MatrixXd selectedOffspring(_glen,mu);
			Eigen::MatrixXd mutOffspring(_glen,mu);
			Eigen::MatrixXd muCenters(_glen,mu);
			Eigen::VectorXd meanSamples(_glen);
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
			pValue = PyInt_FromLong(ssize);
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
				for (j = 0; j < ssize; j++)
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
			for (j = 0; j < ssize; j++)
			{
				for (i = 0; i < _glen; i++)
				{
					copyCenters(i,j) = ind[i];
				}
			}
			// Create the offspring
			offspring = copyCenters + sigma * B * D * samples;
			// Fill and evaluate the current offspring
			for (j = 0; j < ssize; j++)
			{
				for (i = 0; i < _glen; i++)
				 {
				   currInd[i] = offspring(i,j);
				 }
				f[j] = evaluate(currInd, _glen, 0, -1);
				ceval += trials;
			}
			// Sort samples based on fitness
			sortSamples(f, ssize, ordf, ordoff, true);
			for (i = 0; i < _glen; i++)
			{
				for (j = 0; j < ssize; j++)
				{
					sortedSamples(i,j) = samples(i,ordoff[j]);
					sortedOffspring(i,j) = offspring(i,ordoff[j]);
				}
			}
			for (i = 0; i < _glen; i++)
			{
				for (j = 0; j < mu; j++)
				{
					selectedSamples(i,j) = sortedSamples(i,j);
					selectedOffspring(i,j) = sortedOffspring(i,j);
					muCenters(i,j) = ind[i];
				}
			}
			mutOffspring = selectedOffspring - muCenters;
			// Update individual (i.e., center)
			meanSamples = selectedSamples * weights;
			ind = selectedOffspring * weights;
			// Cumulation: update evolutionary paths
			stepPath = (1.0 - csigma) * stepPath + sqrt(csigma * (2.0 - csigma) * muEff) * B * meanSamples;
			// Norm of stepPath		
			/*double sumStepPath = 0.0;
			for (i = 0; i < _glen; i++)
			{
				sumStepPath += (stepPath[i] * stepPath[i]);
			}
			double normStepPath = sqrt(sumStepPath);*/
			double normStepPath;
			double* v = new double[_glen];
			for (i = 0; i < _glen; i++)
				v[i] = stepPath[i];
			long int vdim[] = { _glen };
			// Here we call numpy.linalg.norm() function
			pArgs = PyTuple_New(1);
			pVec = PyArray_SimpleNewFromData( 1, vdim, PyArray_DOUBLE, v );
			if (!pVec)
			{
				Py_DECREF(pArgs);
				fprintf(stderr, "Cannot convert argument\n");
				exit(-2);
		    	}
		   	/* pValue reference stolen here: */
		    	PyTuple_SetItem(pArgs, 0, pVec);
		    	pValue = PyObject_CallObject(pFunc4, pArgs);
			Py_DECREF(pArgs);
		    	if (pValue != NULL) 
			{
				normStepPath = PyFloat_AsDouble(pValue);
				Py_DECREF(pValue);
		    	}
		    	else
			{
				Py_DECREF(pFunc4);
				PyErr_Print();
				fprintf(stderr,"Call failed\n");
				exit(-2);
		    	}
			double hsig = ((normStepPath / (sqrt(1.0 - pow((1.0 - csigma), (2 * (ceval / sampleSize))))) / chiN) < (1.4 + 2.0 / (genLen + 1.0))) ? 1.0 : 0.0;
			covPath = (1.0 - cc) * covPath + hsig * sqrt(cc * (2.0 - cc) * muEff) * B * D * meanSamples;
			// Adapt covariance matrix
			C = (1.0 - covLearnRate) * C + covLearnRate * (1.0 / muEff) * ((covPath * covPath.transpose()) + (1.0 - hsig) * cc * (2.0 - cc) * C) + covLearnRate * (1.0 - (1.0 / muEff)) * (mutOffspring * diag(weights, mu) * mutOffspring.transpose());
			// Adapt step size
			sigma *= std::exp((csigma / dsigma) * (normStepPath / chiN - 1.0));
			// Compute eigenvalues and eigenvectors of C
			Eigen::MatrixXd oldC = C;
			Eigen::MatrixXd CT = C.transpose();
			C = (oldC + CT) / 2.0; // Enforce symmetry
			/*Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
			// Use eigenvalues and eigenvectors to update T
			B = es.eigenvectors();
			D = es.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();*/
			double** m = new double*[_glen];
		    	m[0] = new double[_glen * _glen];
			for (i = 1; i < _glen; i++)
		    		m[i] = m[i - 1] + _glen;				
		    	// fill in values				
		    	for (i = 0; i < _glen; i++)
			{
				for (j = 0; j < _glen; j++)
				{
					m[i][j] = C(i,j);
				}
			}	
		    	long int mdim[] = { _glen, _glen };
				// Here we call numpy.linalg.eig() function
			pArgs = PyTuple_New(1);
			pMat = PyArray_SimpleNewFromData( 2, mdim, PyArray_DOUBLE, m[0] );
			if (!pMat)
			{
				Py_DECREF(pArgs);
				fprintf(stderr, "Cannot convert argument\n");
				exit(-2);
		    	}
		    	PyTuple_SetItem(pArgs, 0, pMat);
		    	pValue = PyObject_CallObject(pFunc5, pArgs);
			Py_DECREF(pArgs);
		    	if (pValue != NULL) 
			{
				pFirst = PyTuple_GetItem(pValue, 0);
				pLast = PyTuple_GetItem(pValue, 1);
				// Extract eigenvalues and eigenvectors
				double* ret = (double*)PyArray_DATA(pFirst);
				for (i = 0; i < _glen; i++)
					D(i,i) = ret[i];
				double* ret2 = (double*)PyArray_DATA(pLast);
				for (i = 0; i < _glen; i++)
				{
					for (j = 0; j < _glen; j++)
					{
						B(i,j) = ret2[j + (i * _glen)];
					}
				}
				for (i = 0; i < _glen; i++)
				{
					for (j = 0; j < _glen; j++)
					{
						if (i != j)
							D(i,j) = 0.0;
					}
				}
				Py_DECREF(pValue);
		    	}
		    	else
			{
				Py_DECREF(pFunc5);
				PyErr_Print();
				fprintf(stderr,"Call failed\n");
				exit(-2);
		    	}
			// Copy back the individual
			for (i = 0; i < _glen; i++)
			{
				_individuals[0][i] = (double)ind[i];
			}
			// Evaluate the individual (to find out its fitness)
			fitcentroid = evaluate(_individuals[0], _glen, 0, -1);
			ceval += trials;
			// the best of the current generation is the centroid
			if (fitcentroid > ordf[0])
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
				if (ordf[0] > bestfit)
				{
					bestfit = ordf[0];
					for (i = 0; i < _glen; i++)
					{
						bestind[i] = offspring(i,ordoff[0]);
					}
				}
				// Generalization
				for (i = 0; i < _glen; i++)
				{
					currInd[i] = offspring(i,ordoff[0]);
				}
				gfit = evaluate(currInd, _glen, 1, -1);
				ceval += ttrials;
				if (gfit > bestgfit)
				{
					bestgfit = gfit;
					for (i = 0; i < _glen; i++)
					{
						bestgsol[i] = offspring(i,ordoff[0]);
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
			for (j = 0; j < ssize; j++)
			{
				f[j] = 0.0;
			}
			// Reset current offspring
			for (i = 0; i < _glen; i++)
			{
				currInd[i] = 0.0;
			}
            printf("Gen %d Steps %d Centroid %.2f BestOffspring %.2f Bestfit %.2f BestgFit %.2f\n", gn, ceval, fitcentroid, ordf[0], bestfit, bestgfit);
			fflush(stdout);
			gn++;

            // store stat data
            *stat = ceval; stat++;
            *stat = bestgfit; stat++;
            *stat = bestfit; stat++;
            *stat = ordf[0]; stat++;
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

	strcpy(filename, filedir);
        sprintf(filen,"covMatS%d.txt", seed);
        strcat(filename, filen);
        fp = fopen(filename,"w!");
        for(i = 0; i < _glen; i++)
		{
			for (j = 0; j < _glen; j++)
			{
        		fprintf(fp, "%lf\t", C(i,j));
			}
			fprintf(fp, "\n");
		}
        fclose(fp);

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


