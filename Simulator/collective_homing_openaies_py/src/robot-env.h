// the robot's structure
struct robot
{
	
	int idn;							// robot's identification number
	int type;                           // robot's type (Khepera, ePuck, MarXBot)
    double x;                           // robot's x position
    double y;                           // robot's y position
    double dir;                         // robot's direction in radiants
    double *sensors;                    // sensors state
    double radius;                      // Robot Radius (mm) - Khepera 27.5m, ePuch 37.5m
    double axleLength;                  // Distance between wheels (mm)
    double maxSpeed;                    // Max linear speed (mm/s)
    double *csensors;                   // pointer to the first sensor to be updated
    double dx;                          // robot's current x offset
    double dy;                          // robot's current y offset
	int color;							// 0=black, 1=red, 2=green, 3=blue, 4=red/green
    bool alive;                         // whether the robot is alive
    double energy;                      // the energy level of the robot
	int sensorinfraredid;				// id of the first infrared sensor neuron
	int sensorinfraredn;                // n of infrared sensory neurons
	int sensorcameraid;				    // id of the first camera sensory neuron
	int sensorcameran;                  // n of infrared sensory neurons
	
	
    // motors
    int motorwheels;                    // number of motors neurons controlling the wheels
    int motorwheelstype;                // type of motor wheel control (0=direct, 1=traslation-rotation)
    int motorleds;                      // number of motors neurons controlling the wheels
	int motorwheelsid;                  // id of the first motorwheels neuron
    int motorledsid;                    // id of the first motorleds neuron

    // variables used by sensors
	double *proprioceptors;				// robots' proprioceptors state
	int nifsensors;                     // number of infrared sensors
    int nlasersensors;                  // number of laser sensors
	int camnsectors;					// the number of camera setctors
	double **camblobs;					// matrix of color blobs delected by the camera, firstlevel (sector) secondlevel (blobs)
	int *camblobsn;                     // the number of color blobs for each sector
};


// the environmental structure (objects contained in the environment)
struct envobjects
{
	int type;
	double x;
	double y;
	double x2;
	double y2;
	double r;
	float color[3];
};

// Robots
#define Khepera 0
#define ePuck   1
#define MarXBot 2


// Environmental variables
#define WALL 0						// wall type (we assume orthogonal with respect to x and y axes)
#define SAMPLEDCYLINDER 1			// sampled cylindrical objects type
#define SAMPLEDSCYLINDER 2			// sampled small cylindrical objects type
#define LIGHTBULB 3					// bulbs emitting lights type
#define RTARGETAREA 4				// rectangular target area type
#define STARGETAREA 5				// spherical target area type

// public variables
extern int trials;					// number of trials
extern int ttrials;					// number of testing trials
extern int steps;					// number of steps x trial
extern int steptimelength;			// time length of a step during graphic visualization
extern int nfreep;					// number of free parameters
extern int renderrobotworld;       // whether we update the graphic renderer

extern struct robot *rob;			// the list of robots is publicly available
extern int nrobots;					// number of robots
extern int nnests;					// number of nests
extern double nestdist;				// minimum nest distance
extern double robotnestdist;        // minimum distance between robot and nest
extern double robotdist;			// minimum robot distance
extern struct envobjects *envobjs;	// list of environemntal objects is publicly available
extern int nenvobjs;				// number of objects is publicly available
extern double worldx;				// world x dimension
extern double worldy;				// world y dimension
extern double *initialstates;		// initial states of the robots and of the environment
extern int steptimelength;			// time length of a step during graphic visualization
extern int cells[100][100];         // virtual cells
extern double cellsize;             // virtual cells size
extern float **wall;				// wall samples
extern int *wallcf;					// description of wall samples
extern float **cylinder;			// cylindrical sample
extern int *cylindercf;				// description of cylindrical sample
extern float **scylinder;			// small cylindrical sample
extern int *scylindercf;			// description of small cylindrical sample
extern float **light;				// light samples
extern int *lightcf;				// description of light samples

extern char globdir[512];

// public functions
double mrand(double range);
double pmrand(double range);
double angv(double x1, double x2, double y1, double y2);
double xvect(double ang, double module);
double yvect(double ang, double module);
double mangrel(double absolute, double orient);
double mangrelr(double absolute, double orient);
void save_obstacle(char *filename, int **object, int  *objectcf);
int ** load_obstacle_old(char *filename, int  *objectcf);
float ** load_obstacle(char *filename, int  *objectcf);
void initRobot(struct robot *cro, int n, int robottype);
void initRobotSensors(struct robot *cro, int nsensors);
int initInfraredSensor(struct robot *rob, int nsensors, int ids);
void updateInfraredSensor(struct robot *rob);
int initCameraSensorRFB(struct robot *cro, int nsectors);
void updateCameraSensorRFB(struct robot *cro, int *rbd);
int initCameraPPSensor(struct robot *cro);
void updateCameraPPSensor(struct robot *cro);
int initProprioceptors(struct robot *ro);
void updateProprioceptors(struct robot *ro);
int initEnergySensor(struct robot *ro);
void updateEnergySensor(struct robot *ro);
int initBiasSensor(struct robot *rob);
void updateBiasSensor(struct robot *rob);
int initTimeSensor(struct robot *cro);
void updateTimeSensor(struct robot *cro, int step, int nsteps);
int initGroundSensor(struct robot *cro);
void updateGroundSensor(struct robot *cro);
int initGroundGradSensor(struct robot *cro);
void updateGroundGradSensor(struct robot *cro);
int initLaserDistanceSensor(struct robot *cro);
void updateLaserDistanceSensor(struct robot *cro);


int robotCheckCollision(struct robot *rob);
int updateRobot(struct robot *cro, struct evonet *ne);
void savegenotype(char* filename, double* genotype, const int glen, int mode);

void update_world();				// THIS IS INCLUDED IN MAIN.CPP

