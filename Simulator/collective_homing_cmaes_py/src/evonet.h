/*
 * evonet.h
 */

#ifndef __cplusplus
#define false 0
#define true 1
typedef int bool;
#endif

#ifndef MYBOOLEAN_H
#define MYBOOLEAN_H


#define MAXN 100
#define MAX_BLOCKS 20

struct evonet
{

	//! whether the network has a full recurrent architecture
	bool fullRecurrent;
	//! whether we have input-output connections
	bool inputOutputConnections;
	//! whether we have recurrent connections on the hiddens
	bool recurrentHiddens;
	//! whether we have recurrent connections on the outputs
	bool recurrentOutputs;
	//! actvations
	double *act;
	//! Net-inputs
	double *neti;
	//! Inputs
	//double *input;
	//! Number of inputs
	int ninputs;
	//! Number of hiddens
	int nhiddens;
	//! Number of outputs
	int noutputs;
	//! Number of neurons
	int nneurons;
	//! Weight range
	double wrange;
	//! Bias range
	double brange;
	//! neurons' activation function (defoult 1 = logistic)
	int afunction;
    //! output neurons' activation function (default = afunction)
    int aofunction;
	//! whether (hidden and output) neurons have biases
	bool biases;
	//! The type of neurons
	int *neurontype;
	//! Whether neurons have bias
	int *neuronbias;
	//! Description of network architecture
	int *netblock;
	//! Number of blocks
	int nblocks;
};

extern struct evonet *net;			// the list of networks publicly available
extern int verbose;                 // verbose level

// public functions
void copyandclear(char *s, char *sc);
int computeParameters(struct evonet *net);
void updateNet(struct evonet *net, double *inputs, double *freep);
void resetNet(struct evonet *net);
double getOutput(struct evonet *net, int i);
int initNetArchitecture(struct evonet *net);
void initNet(struct evonet *net, const char* filename, int nid);

#endif
