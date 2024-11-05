#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "evonet.h"
#include "utilities.h"

#define MAXN 100
#define MAX_BLOCKS 20

// other variables
// display level (0=null, 1=number required parameters & genotye 2= also state of the network every step
int verbose=0;

/*
 * compute the number of free parameters
 * this function is included but it is not used
 */
int computeParameters(struct evonet *net)
{
	int i;
	int t;
	int b;
	int ng;
	int *nbl;
    
	ng  = 0;
	
    // biases
    for(i=0;i < net->nneurons;i++)
       {
        if (net->neuronbias[i] == 1)
            ng++;
        // plogistic neurons require an additional parameter
        if (net->neurontype[i] == 2)
            ng++;

    }

    /* timeconstants
    for(i=0, nt = net->neurontype;i < net->nneurons;i++, nt++) {
        if (*nt == 1) {
            ng++;
        }
    } */

    // blocks
    for (b=0, nbl = net->netblock; b < net->nblocks; b++)
    {
        // connection block
        if (*nbl == 0)
        {
            for(t=0; t < *(nbl + 2);t++)
            {
                for(i=0; i < *(nbl + 4);i++)
                {
                    ng++;
                }
            }
        }
        nbl = (nbl + 5);

    }
    
	return(ng);
}


/*
 * standard logistic
 */
double logistic(double f)
{
	return((double) (1.0 / (1.0 + exp(0.0 - f))));
}

/*
 * parametric logistic
 */
double plogistic(double f, double d)
{
    return((double) (1.0 / (1.0 + exp(0.0 - f * d))));
}

/*
 * hyperbolic tangent
 */
double hyptan(double f)
{
	if (f >= 10.0)
		return 1.0;
	else if (f <= -10.0)
		return -1.0;
	else
		return ((double) ((1.0 - exp(0.0 - 2 * f)) / (1.0 + exp(0.0 - 2 * f))));
}

/*
 * update the network
 */
void updateNet(struct evonet *net, double *inputs, double *freep)
{
    int i;
    int t;
    int b;
    double *p;
    double *a;
    double *ni;
    int *nt;
    double delta;
    int *nbl;


    p = freep;

    // biases
    for(i=0, ni = net->neti;i < net->nneurons;i++, ni++)
    {
        if (net->neuronbias[i] == 1)
        {
            *ni = *p;//((double)*p/net->wrange)*net->brange;
            p++;
        }
        else
        {
            *ni = 0.0;
        }
    }


    // blocks
    for (b=0, nbl = net->netblock; b < net->nblocks; b++)
    {
        // connection block
        if (*nbl == 0)
        {
            for(t=0, ni=(net->neti + *(nbl + 1)); t < *(nbl + 2);t++, ni++)
            {
                for(i=0, a=(net->act + *(nbl + 3)); i < *(nbl + 4);i++, a++)
                {
                    *ni += *a * *p;
                    p++;
                }
            }
        }

        // update block
        if (*nbl == 1)
          {
            for(t=*(nbl + 1), a =(net->act + *(nbl + 1)), ni=(net->neti + *(nbl + 1)),nt = (net->neurontype + *(nbl + 1)); t < (*(nbl + 1) + *(nbl + 2)); t++, a++, ni++, nt++)
                {
                  switch(*nt)
                    {
                      case 0: // input neurons are  simple rely units
                        *a = *(inputs + t); break;
                      case 1: // logistic function
                        *a = logistic(*ni); break;
                      case 2: // parametric logistic function
                        delta = fabs(*p); p++; *a = plogistic(*ni, delta); break;
                      //case 3: // delta function
                        //delta = fabs(*p); p++; *a = (*a * delta)  + (logistic(*ni) * (1.0f - delta)); break;
                      case 3: // binary function
                        if (*ni >= 0.0) *a = 1.0; else *a = 0.0; break;
                      case 4: // rectified linear function
                        if (*ni >= 0.0) *a = *ni; else *a = 0.0; break;
                      case 5: // tanh
                        *a = hyptan(*ni); break;
                    }
                }
            }
        nbl = (nbl + 5);
    }
}



void resetNet(struct evonet *net)
{
    double *a;
    int i;
    
	for (i = 0, a = net->act; i < net->nneurons; i++, a++)
	 {
		*a = 0.0;
	 }
}


double getOutput(struct evonet *net, int i)
{
    double *a;
    
    a = (net->act + (net->ninputs + net->nhiddens + i));
    return *a;
}

/*
 * Initialize the neural network architecture
 */
int initNetArchitecture(struct evonet *net)
{

    int n,nn;
    int *nt;
    int *nbl;
    int *nb;
    
    net->nneurons = net->ninputs + net->nhiddens + net->noutputs;
	
	//initialize
	net->nblocks = 0;
	nbl = net->netblock;
	
    //neurons' type and biases
    for (n=0, nt = net->neurontype, nb = net->neuronbias; n < net->nneurons; n++, nt++, nb++)
      {
          if (n < net->ninputs)
            *nt = 0;
          else
            if (n < (net->ninputs + net->nhiddens))
              *nt = net->afunction;
              else
              *nt = net->aofunction;

          if (net->biases == 1 && n >= net->ninputs)
            *nb = 1;
           else
            *nb = 0;
      }

    /*biases
	for (n=0, nt = net->neurontype; n < net->nneurons; n++, nt++)
      {
          net->neuronbias[n] = false;
          *nt = 0;
      }
     */

	// input update block
	*nbl = 1; nbl++;
	*nbl = 0; nbl++;
	*nbl = net->ninputs; nbl++;
	*nbl = 0; nbl++;
	*nbl = 0; nbl++;
	net->nblocks++;
	
    // WE SHOULD USE A SINGLE BLOCK OF CONNECTION HERE
	if (net->fullRecurrent)
	{
	  // input-neurons connections
	  *nbl = 0; nbl++;
	  *nbl = net->ninputs; nbl++;
	  *nbl = net->nhiddens+net->noutputs; nbl++;
	  *nbl = 0; nbl++;
	  *nbl = net->ninputs; nbl++;
	  net->nblocks++;
		
	  // neurons-neurons connections
	  *nbl = 0; nbl++;
	  *nbl = net->ninputs; nbl++;
	  *nbl = net->nhiddens+net->noutputs; nbl++;
	  *nbl = net->ninputs; nbl++;
	  *nbl = net->nhiddens+net->noutputs; nbl++;
	  net->nblocks++;
		
	  // hidden update block
	  *nbl = 1; nbl++;
	  *nbl = net->ninputs; nbl++;
	  *nbl = net->nhiddens+net->noutputs; nbl++;
	  *nbl = 0; nbl++;
	  *nbl = 0; nbl++;
	  net->nblocks++;
	}
	else
	{
	
      // input-hidden connections
      if (net->nhiddens > 0) {
          *nbl = 0; nbl++;
          *nbl = net->ninputs; nbl++;
          *nbl = net->nhiddens; nbl++;
          *nbl = 0; nbl++;
          *nbl = net->ninputs; nbl++;
          net->nblocks++;
      }
    
      // hidden-hidden connections
      if (net->recurrentHiddens) {
          *nbl = 0; nbl++;
          *nbl = net->ninputs; nbl++;
          *nbl = net->nhiddens; nbl++;
          *nbl = net->ninputs; nbl++;
          *nbl = net->nhiddens; nbl++;
          net->nblocks++;
      }
    
      // hidden update block
      if (net->nhiddens > 0) {
          *nbl = 1; nbl++;
          *nbl = net->ninputs; nbl++;
          *nbl = net->nhiddens; nbl++;
          *nbl = 0; nbl++;
          *nbl = 0; nbl++;
          net->nblocks++;
      }
    
      // input-output connections
      if (net->nhiddens == 0 || net->inputOutputConnections) {
          *nbl = 0; nbl++;
          *nbl = net->ninputs + net->nhiddens; nbl++;
          *nbl = net->noutputs; nbl++;
          *nbl = 0; nbl++;
          *nbl = net->ninputs; nbl++;
          net->nblocks++;
      }
    
      // hidden-output connections
      if (net->nhiddens > 0) {
          *nbl = 0; nbl++;
          *nbl = net->ninputs + net->nhiddens; nbl++;
          *nbl = net->noutputs; nbl++;
          *nbl = net->ninputs; nbl++;
          *nbl = net->nhiddens; nbl++;
          net->nblocks++;
      }
    
      // output-output connections
      if (net->recurrentOutputs) {
          *nbl = 0; nbl++;
          *nbl = net->ninputs + net->nhiddens; nbl++;
          *nbl = net->noutputs; nbl++;
          *nbl = net->ninputs + net->nhiddens; nbl++;
          *nbl = net->noutputs; nbl++;
          net->nblocks++;
      }
    
      // output update block
      *nbl = 1; nbl++;
      *nbl = net->ninputs + net->nhiddens; nbl++;
      *nbl = net->noutputs; nbl++;
      *nbl = 0; nbl++;
      *nbl = 0; nbl++;
      net->nblocks++;
	
	  }
  /* display blocks
  if (verbose >= 2)
  {
	nbl = net->netblock;
	for(n=0; n < net->nblocks; n++)
	{
	    printf("block ");
		for (nn=0; nn < 5; nn++, nbl++)
	      printf("%d ", *nbl);
	    printf("\n");
	}
  }
  */

   return(0);
	
}


void initNet(struct evonet *net, const char* filename, int nid)
{

    char *s;
    char buff[1024];
    char name[1024];
    char value[1024];
    char *ptr;
    int section;  // 0=before the section 1=in the section 2= after the section

    section = 0;
    // default
    net->fullRecurrent = false;
    net->inputOutputConnections = false;
    net->recurrentHiddens = false;
    net->recurrentOutputs = false;
    net->nhiddens = 5;
    net->biases = false;
    net->afunction = 1;  // default: logistic
    net->aofunction = 0; // default: the same of afunction

    verbose = 2;

    // read parameters from .ini file
    FILE* fp = fopen(filename, "r");
    if (fp != NULL)
    {
        // Read lines
        while (fgets(buff, 1024, fp) != NULL && section < 2)
        {

            // check whether the relevant section start or end
            if (buff[0] == '[')
            {
              if (section == 1)
                {
                    section = 2;
                    continue;
                }
              if ((section == 0) && (strncmp(buff, "[NET]",5)==0))
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

            // parse parameters
            if (strcmp(name, "nhiddens")==0)
            net->nhiddens = (int)strtol(value, &ptr, 10);
            else if (strcmp(name, "biases")==0)
            net->biases = (bool)strtol(value, &ptr, 10);
            else if (strcmp(name, "afunction")==0)
            net->afunction = (int)strtol(value, &ptr, 10);
            else if (strcmp(name, "aofunction")==0)
            net->aofunction = (int)strtol(value, &ptr, 10);
            else if (strcmp(name, "fullRecurrent")==0)
            net->fullRecurrent = (bool)strtol(value, &ptr, 10);
            else if (strcmp(name, "inputOutputConnections")==0)
            net->inputOutputConnections = (bool)strtol(value, &ptr, 10);
            else if (strcmp(name, "recurrentHiddens")==0)
            net->recurrentHiddens = (bool)strtol(value, &ptr, 10);
            else if (strcmp(name, "recurrentOutputs")==0)
            net->recurrentOutputs = (bool)strtol(value, &ptr, 10);
            else if (strcmp(name, "prange")==0)
            net->wrange = net->brange = strtod(value, &ptr);
            else if (strcmp(name, "verbose")==0)
            verbose = (int)strtol(value, &ptr, 10);
            else if (nid == 0) printf("WARNING: Unknown parameter %s in section [NET] of file %s \n", name, filename);
            }
        }
        fclose (fp);
        if (section == 0 && nid == 0)
           printf("WARNING: Missing section [NET] in file %s \n", filename);

        // default values
        if (net->aofunction == 0)
            net->aofunction = net->afunction;
        if (verbose >= 1 && nid == 0)
        {
          printf("NET: Neurons: %d->%d->%d, ", net->ninputs, net->nhiddens, net->noutputs);
            if (net->fullRecurrent == 0)
               printf("Topology: IOConn %d RecH %d RecO %d, ", net->inputOutputConnections, net->recurrentHiddens, net->recurrentOutputs);
             else
               printf("Topology: FullRecurrent, ");
            if (net->biases == 0)
               printf("NoBias, ");
             else
               printf("Bias, ");
            switch(net->afunction)
            {
              case 1: printf("aFunct: logistic "); break;
              case 2: printf("aFunct: plogistic "); break;
              case 3: printf("aFunct: binary "); break;
              case 4: printf("aFunction: rectified linear "); break;
              case 5: printf("aFunction: tanh "); break;
            }
            if (net->afunction != net->aofunction)
             switch(net->aofunction)
              {
                case 1: printf("aoFunct: logistic "); break;
                case 2: printf("aoFunct: plogistic "); break;
                case 3: printf("aoFunct: binary "); break;
                case 4: printf("aoFunction: rectified linear "); break;
                case 5: printf("tanh "); break;
              }
           printf("\n");
        }
    }
    else
    {
        printf("ERROR: unable to open file %s\n", filename);
        fflush(stdout);
    }

    // initialize the number of neurons
    net->nneurons = net->ninputs + net->nhiddens + net->noutputs;
    // allocate memory for network vectors and blocks
    net->act=(double*)malloc(net->nneurons*sizeof(double));
    net->neti=(double*)malloc(net->nneurons*sizeof(double));
    net->neurontype=(int*)malloc(net->nneurons*sizeof(int));
    net->neuronbias=(int*)malloc(net->nneurons*sizeof(int));
    net->netblock = (int *)malloc(MAX_BLOCKS * 5 * sizeof(int));
    // define the architecture of the network
    initNetArchitecture(net);
    // initialize neuron activation to 0.0
    resetNet(net);

}







