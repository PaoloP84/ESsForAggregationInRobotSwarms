void update_world();				// THIS FUNCTION IS INCLUDED IN MAIN.CPP

// public discrim.c functions
bool readEvoConfig(const char* filename);
void initEnvironment();
double calcNestRadius(int nr);
int calcRobots(double radius);
void createInitialStates(int mode);
int initialize(const char* filename);
double evaluate(double* genotype, const int glen, int mode, int seed);
