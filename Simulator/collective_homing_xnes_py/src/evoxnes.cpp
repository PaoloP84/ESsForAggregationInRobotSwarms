/*
 * xNES evolutionary algorithm
 * it require the Eigen and unsupported librieris (costituted only by include files)
 * and the GSL libraries (included in the utilities.h file) for the generation of random numbers
 * it call an eternal evalaate function that is used to evaluate candidate solutions
 */
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <time.h>
#include <typeinfo>
#include "utilities.h"
#include "evoxnes.h"

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

//! xNES evolution algorithm
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
	Eigen::VectorXf ind(_glen);				// The individual from which the offspring will be generated
	Eigen::VectorXf dInd(_glen);			// The update of the individual
	Eigen::MatrixXf covMatrix(_glen,_glen); // Covariance Matrix
	Eigen::MatrixXf dCovM(_glen,_glen);		// The update of the covariance matrix
	Eigen::MatrixXf expCovM(_glen,_glen);	// Exponential matrix of CovM
	double etaInd;							// Learning rate for the individual
	double etaCovM;							// Learning rate for the covariance matrix
	Eigen::VectorXf utility(ssize);			// Array containing the utility values
	double utilSum;							// Sum of utility values
	double tmp;								// Temporary variable
	Eigen::VectorXf orderedUtility(ssize);  // Array containing the ordered utilities (i.e. the best individual has the highest utility value)
	double ordUtilSum;						// Sum of ordered utility values (it shall be equal to 1.0)
	double* f = new double[ssize];			// Array containing the fitness of the offspring
	Eigen::VectorXf ordf(ssize);			// Array containing the ordered fitness of the offspring
	int* ordoff = new int[ssize];			// Array containing the ordered indices of the offspring
	double* currInd = new double[_glen];	// Current individual to be evaluated (temporary variable)
	Eigen::MatrixXf eye(_glen,_glen);		// Identity matrix
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

    printf("xNES seed %d MaxEval %ld prange %.1f\n", seed, maxevaluations, prange);
	
    time_t startTimer;
    time_t endTimer;
    time(&startTimer);
    // set the seed and generate the initial states
    RNG->setSeed(seed);
    createInitialStates(0);

    // allocate memory for statistics
    statn = (maxevaluations / (trials * ssize) * 7) + 100;
    statistics = (double *) malloc(statn * sizeof(double));
    stat = statistics;
    statn = 0;

	// Initialize learning rates
	etaInd = 1.0;
	etaCovM = (9.0 + 3.0 * log(genLen)) / (5.0 * genLen * sqrt(genLen));
	_stepRate = 1.0;
	// Calculate the utility values before normalization
	utilSum = 0.0;
	for (j = 0; j < ssize; j++)
	{
		// u[j] is equal to the uttermost between 0 and (log(lambda / 2 + 1) - log(j))
		tmp = log(sampleSize / 2.0 + 1.0) - log((double)j + 1.0); // Remember that log(0) = -Inf, so 1 is added to index i
		if (tmp < 0.0)
		{
			tmp = 0.0;
		}
		utility[j] = tmp;
		// Update the sum of utility values
		utilSum += utility[j];
	}
	// Normalize the utility values
	for (j = 0; j < ssize; j++)
	{
		if (utilSum != 0.0)
		{
			utility[j] /= utilSum;
		}
		utility[j] -= 1.0 / genLen; // In the standard version, utility values might be negative
	}
	for (i = 0; i < _glen; i++)
	{
		// Define the identity matrix
		eye(i,i) = 1.0;
		for (j = 0; j < _glen; j++)
		{
			if (i != j)
			{
				eye(i,j) = 0.0;
			}
			// Initialize covariance matrix, its update and its exponential
			covMatrix(i,j) = 0.0;
			dCovM(i,j) = 0.0;
			//expCovM(i,j) = 0.0;
		}
		// Initialize individual and its update
		ind[i] = 0.0;
		dInd[i] = 0.0;
		// Initialize the individual (temporary variable) used for evaluation
		currInd[i] = 0.0;
	}
	printf("Seed %d \n", seed);
	fflush(stdout);
	
		bestfit = -999999999.0;
		bestgfit = -999999999.0;
		// Initialize the individual for this replication
		for (i = 0; i < _glen; i++)
        {
            //_individuals[0][i] = - prange + ( double( rand() ) / double( RAND_MAX ) ) * (prange * 2.0);
            _individuals[0][i] = RNG->getDouble(-prange, prange);
			// Reset the current offspring (to avoid numerical issues)
			currInd[i] = 0.0;
		}
		// Initialize Covariance Matrix
		for (i = 0; i < _glen; i++)
		{
			for (j = 0; j < _glen; j++)
			{
				covMatrix(i,j) = 0.0;
				//expCovM(i,j) = 0.0;
			}
		}
		// Initialize fitnesses, ordered offspring and ordered utilities
		for (j = 0; j < ssize; j++) // Remember that the last slot is reserved for the parent
		{
			f[j] = 0.0;
			ordf[j] = 0.0;
			ordoff[j] = 0.0;
			orderedUtility[j] = 0.0;
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
			
			// Calculate the matrix exponential
			Eigen::MatrixExponential<Eigen::MatrixXf>(covMatrix).compute(expCovM);
			// Extracting samples (i.e. generate offspring of current individual)
			Eigen::MatrixXf samples(_glen,ssize);
			Eigen::MatrixXf trSamples(_glen,ssize);
			Eigen::MatrixXf offspring(_glen,ssize);
			
			// Draw out samples from a Gaussian distribution 
			for (i = 0; i < _glen; i++)
			{
				for (j = 0; j < ssize; j++)
				{
                    //samples(i,j) = rand_gaussian(0.0,1.0);
                    samples(i,j) = RNG->getGaussian(1.0, 0.0);
				}
			}
			// Create the offspring
			trSamples = expCovM * samples;
			for (j = 0; j < ssize; j++)
			{
				for (i = 0; i < _glen; i++)
				{
					offspring(i,j) = ind[i] + trSamples(i,j);
				}
			}
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
			
			// Now fitnesses have to be ordered in descending order
			// (in the original XNES, fitnesses are ordered in ascending order)
			int startIndex = 0;
			bool* examined = new bool[ssize];
			for (j = 0; j < ssize; j++)
			{
				examined[j] = false;
			}
			int countExamined = 0;
			while (countExamined < ssize)
			{
				bool found = false;
				i = 0;
				double maxf;
				int indMaxf;
				while (i < ssize && !found)
				{
					if (!examined[i])
					{
						maxf = f[i];
						indMaxf = i;
						found = true;
					}
					i++;
				}
				for (j = 0; j < ssize; j++)
				{
					if (f[j] > maxf && !examined[j])
					{
						maxf = f[j];
						indMaxf = j;
					}
				}
				ordf[startIndex] = maxf;
				ordoff[startIndex] = indMaxf;
				examined[indMaxf] = true;
				countExamined++;
				startIndex++;
			}
			delete examined;
			// Initialize the sum of ordered utility values
			ordUtilSum = 0.0;
			// Now fill the array with the utility values
			for (j = 0; j < ssize; j++)
			{
				orderedUtility[ordoff[j]] = utility[j];
				ordUtilSum += orderedUtility[ordoff[j]];
			}
			// Compute the natural gradient for both the individual and the covariance matrix
			// Create a matrix whose rows contain weights * samples
			Eigen::MatrixXf ordUtil(_glen,ssize);
			Eigen::MatrixXf ordUtilSamples(_glen,ssize);
			Eigen::MatrixXf Z(ssize,_glen);
			for (i = 0; i < _glen; i++)
			{
				for (j = 0; j < ssize; j++)
				{
					ordUtil(i,j) = orderedUtility[j];
				}
			}
			for (i = 0; i < _glen; i++)
			{
				for (j = 0; j < ssize; j++)
				{
					ordUtilSamples(i,j) = ordUtil(i,j) * samples(i,j);
					Z(j,i) = samples(i,j);
				}
			}
			Eigen::MatrixXf G = ordUtilSamples * Z - ordUtilSum * eye;
			// Define transpose matrices for both samples and ordered utilities
			// (they are necessary to perform the calculation of gradients)
			Eigen::MatrixXf ordUtilTr(ssize,1);
			Eigen::MatrixXf v(_glen,1);
			for (j = 0; j < ssize; j++)
			{
				ordUtilTr(j,0) = orderedUtility[j];
			}

			v = etaInd * expCovM * samples * ordUtilTr;
			for (i = 0; i < _glen; i++)
			{
				dInd[i] = v(i,0);
			}
			dCovM = etaCovM * G;
			// Update the individual and the covariance matrix
			for (i = 0; i < _glen; i++)
			{
				ind[i] += (_stepRate * dInd[i]);
			}
			covMatrix += (_stepRate * dCovM);
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
				orderedUtility[j] = 0.0;
			}
			// Reset matrix exponential
			/*for (i = 0; i < _glen; i++)
			{
				for (j = 0; j < _glen; j++)
				{
					expCovM(i,j) = 0.0;
				}
			}*/
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
        		fprintf(fp, "%lf\t", covMatrix(i,j));
			}
			fprintf(fp, "\n");
		}
        fclose(fp);

	strcpy(filename, filedir);
        sprintf(filen,"expCovMatS%d.txt", seed);
        strcat(filename, filen);
        fp = fopen(filename,"w!");
        for(i = 0; i < _glen; i++)
		{
			for (j = 0; j < _glen; j++)
			{
        		fprintf(fp, "%lf\t", expCovM(i,j));
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


