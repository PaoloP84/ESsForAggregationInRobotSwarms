/*
 * discrim.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <locale.h>
#include "robot-env.h"
#include "evonet.h"
#include "utilities.h"
#include <vector>

#define LOCAL_SEED 12345

int robottype = Khepera;			// the robot we are using
int  nfoods;						// number of food elements
double *robotsdist;                 // matrix containing robot distances;
int *robotsbydistance;              // matrix containing the id of the robots ordered by the inverse of the distance
double wsize = 5000.0; // world size
double asize; // nest size (maximum)
double msize; // minimum size
double nestBonus = 100.0; // bonus awarded for the ability to reach a nest
double kdist = 1.0; // coefficient weighting the distance component of fitness
bool camera = true; // Presence/absence of camera sensor (default is true, i.e. presence)
bool fixedpos = false; // Fixed positions for nests and robots (default is false)

/*
 * read parameters from the configuration file
 */
bool readEvoConfig(const char* filename)
{
    char *s;
    char buff[1024];
    char name[1024];
    char value[1024];
    char *ptr;
    int section;  // 0=before the section 1=in the section 2= after the section

    section = 0;

    // Setting default values for distances (to avoid any numerical issue)
    nestdist = 1000.0;
    robotnestdist = 500.0;
    robotdist = 200.0;

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
              if ((section == 0) && (strncmp(buff, "[EXP]",5)==0))
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
            if (strcmp(name, "nrobots")==0)
            nrobots = (int)strtol(value, &ptr, 10);
            else if (strcmp(name, "wsize")==0)
            wsize = strtod(value,&ptr);
            else if (strcmp(name, "nnests")==0)
            nnests = (int)strtol(value, &ptr, 10);
            else if (strcmp(name, "nestdist")==0)
            nestdist = strtod(value, &ptr);
            else if (strcmp(name, "robotnestdist")==0)
            robotnestdist = strtod(value, &ptr);
            else if (strcmp(name, "robotdist")==0)
            robotdist = strtod(value, &ptr);
            else if (strcmp(name, "nestbonus")==0)
            nestBonus = strtod(value, &ptr);
            else if (strcmp(name, "kdist")==0)
            kdist = strtod(value, &ptr);
            else if (strcmp(name, "trials")==0)
            trials = (int)strtol(value, &ptr, 10);
            else if (strcmp(name, "ttrials")==0)
            ttrials = (int)strtol(value, &ptr, 10);
            else if (strcmp(name, "steps")==0)
            steps = (int)strtol(value, &ptr, 10);
            else if (strcmp(name, "steptimelength")==0)
            steptimelength = (int)strtol(value, &ptr, 10);
			      else if (strcmp(name, "robottype")==0)
            robottype = (int)strtol(value, &ptr, 10);
            else if (strcmp(name, "camera")==0)
            {
              int cam = strtol(value, &ptr, 10);
              camera = (bool)cam;
            }
            else if (strcmp(name, "fixedpos")==0)
            {
              int fixed = strtol(value, &ptr, 10);
              fixedpos = (bool)fixed;
            }
            else printf("WARNING: Unknown parameter %s in section [EXP] of file %s \n", name, filename);
         }
        }
        fclose (fp);
        if (section == 0)
           printf("WARNING: Missing section [NET] in file %s \n", filename);

        printf("Robot %d nrobots %d Trials %d TTrials %d Steps %d\n", robottype, nrobots, trials, ttrials, steps);
        return(true);
    }
    else
    {
        printf("ERROR: unable to open file %s\n", filename);
        fflush(stdout);
        return(false);
    }
}


/*
 * initialize the environment
 */
void
initEnvironment()

{

	int cobj=0;
	//double wsize;
	//double asize;
	int f;
	
    //nfoods = 2;//1;
	nenvobjs = 4 + nnests;//nfoods;	// total number of objects
	switch (robottype)
	 {
	   case (Khepera):
	    //wsize = 2000.0;
		asize = 120.0;
		break;
	   case (ePuck):
	    //wsize = 2500.0;
		asize = 160.0;
		break;
	   case (MarXBot):
	    //wsize = 5000.0; // Arena size(s): it is squared!!!
		asize = 400.0; // Maximum radius of a nest???
		break;
	  }
	
	worldx = wsize;
	worldy = wsize;

	
	envobjs = (struct envobjects *) malloc(nenvobjs * sizeof(struct envobjects));
	
	envobjs[cobj].type = WALL;
	envobjs[cobj].x = 0.0;
	envobjs[cobj].y = 0.0;
	envobjs[cobj].x2 = worldx;
	envobjs[cobj].y2 = 0.0;
	cobj++;
	envobjs[cobj].type = WALL;
	envobjs[cobj].x = 0.0;
	envobjs[cobj].y = 0.0;
	envobjs[cobj].x2 = 0.0;
	envobjs[cobj].y2 = worldy;
	cobj++;
	envobjs[cobj].type = WALL;
	envobjs[cobj].x = worldx;
	envobjs[cobj].y = 0.0;
	envobjs[cobj].x2 = worldx;
	envobjs[cobj].y2 = worldy;
	cobj++;
	envobjs[cobj].type = WALL;
	envobjs[cobj].x = 0.0;
	envobjs[cobj].y = worldy;
	envobjs[cobj].x2 = worldx;
	envobjs[cobj].y2 = worldy;
	cobj++;
	for(f=0; f < nnests; f++)//nfoods; f++)
	{
	  envobjs[cobj].type = STARGETAREA;
	  envobjs[cobj].x = 400.0;
	  envobjs[cobj].y = 400.0;
	  envobjs[cobj].r = asize;
    envobjs[cobj].color[0] = 0.5;
	  cobj++;
	}
	
	if (cobj > nenvobjs)
		{
			printf("ERROR: you should allocate more space for environmental objects");
			fflush(stdout);
		}
	
}

double calcNestRadius(int nr)
{
	double r;
	if (nr <= 0)
	{
		printf("Invalid number of robots %d!!\n", nr);
		exit(-1);
	}
	if (nr == 1)
		r = rob->radius * 2.0;
	else
		r = rob->radius * (1.0 + (1.0 / sin((2.0 * M_PI) / (2.0 * nr))));
	return r;
}

int calcRobots(double radius)
{
	int n;
	double num;
	if (radius <= 0.0)
	{
		printf("Invalid radius %lf\n", radius);
		exit(-1);
	}
	// To compute the number of robots that can occupy the area, we must invert the formula computing the radius of thea area given the robots
	// Compute ratio between radii
	double r_ratio = radius / rob->radius;
	// Now compute argument of arcsin operation...
	// Remember that r = rob->radius * (1.0 + (1.0 / sin((2.0 * M_PI) / (2.0 * nr))))
	// Therefore:
	// 1) r / rob->radius = 1.0 + (1.0 / sin((2.0 * M_PI) / (2.0 * nr)))
	// 2) (r / rob->radius - 1.0) = 1.0 / sin((2.0 * M_PI) / (2.0 * nr))
	// 3) 1.0 / (r / rob->radius - 1.0) = sin((2.0 * M_PI) / (2.0 * nr))
	// 4) arcsin((1.0 / (r / rob->radius - 1.0))) = M_PI / nr
	// 5) nr = M_PI / arcsin((1.0 / (r / rob->radius - 1.0)))
	double arg = asin(1.0 / (r_ratio - 1.0));
	num = M_PI / arg;
	n = ((double)num);
	return n;
}

/*
 * generate the initial states
 * if mode != 0 regenerate only the trial initial conditions
 */
void createInitialStates(int mode)

{
    
  int s;
  double *is;
  int numtrials;
  int attempts;
  double dx, dy;
  struct robot *ro1;
  struct robot *ro2;
  int r1, r2;
  double cdist, mindist;
  double fcx, fcy, fcr;
  double distfromborder;
  double distfromareacenter;
  bool placed;
  double dist;
  int nest;
  int other;
  double* nestPos = (double*)malloc(sizeof(double) * (nnests * 3)); // x and y
  double minsize;
  int nthrobots;
  int nleftrobots;

  if (mode == 0)
    numtrials = trials + ttrials;
  else
    numtrials = trials;

  switch (robottype)
  {
    case (Khepera):
    distfromborder = 300.0;
    distfromareacenter = 150.0;
    break;
    case (ePuck):
    distfromborder = 400.0;
    distfromareacenter = 200.0;
    break;
    case (MarXBot):
    distfromborder = asize + 200.0;//300.0;//1000.0;
    distfromareacenter = 500.0;
    break;
  }

  if (fixedpos)
  {
    printf("#### TEST WITH FIXED INITIAL STATES ####\n");
    printf("ttrials = %d\n", numtrials);
  }

  for (s=0, is = initialstates; s < numtrials; s++)
  {
    if (fixedpos)
      // Set constant seed. This should ensure that the sequence of pseudo-random numbers
      // is the same for all the trials. Consequently, the initial position of robots and 
      // the initial position and radius of the nests are the same!!!
      RNG->setSeed(LOCAL_SEED);
    //printf("trial %d\n", s);
	minsize = msize;
	nthrobots = 0;
	if (nnests == 1)
	{
		// Extract random position and radius
        	fcx = RNG->getDouble(distfromborder, worldx - distfromborder);
        	fcy = RNG->getDouble(distfromborder, worldy - distfromborder);
		// Radius is the maximum possible
		fcr = asize;
		// Set initialstates
	     	*is = fcx;
	     	is++;
	      	*is = fcy;
	     	is++;
	      	*is = fcr;
	      	is++;
	      	nestPos[0] = fcx;
      		nestPos[1] = fcy;
		nestPos[2] = fcr;
	}
	else
	{
	    	for (nest = 0; nest < nnests; nest++)
	    	{
	      placed = false;
	      while (!placed)
	      {
		// Extract random position and radius
		fcx = RNG->getDouble(distfromborder, worldx - distfromborder);
		fcy = RNG->getDouble(distfromborder, worldy - distfromborder);
		//printf("asize %lf\n", asize);
		fcr = RNG->getDouble(minsize, asize);//asize / 4.0, asize);
		//printf("nest %d - radius %lf\n", nest, fcr);
		// Check whether nest position overlaps with already placed nests, if any
		placed = true;
		other = 0;
		while ((other < nest) && placed)
		{
		  //printf("nest %d - other %d\n", nest, other);
		  // Compute distance with already placed nests
		  dx = (fcx - nestPos[3 * other]);
		  dy = (fcy - nestPos[3 * other + 1]);
		  dist = sqrt(dx*dx+dy*dy);
		  if (dist < nestdist)
		    placed = false;
		  other++;
		}
	      }
	      // Set initialstates
	      *is = fcx;
	      is++;
	      *is = fcy;
	      is++;
	      *is = fcr;
	      is++;
	      //printf("placed nest %d\n", nest);
	      // Update nest position array
	      nestPos[3 * nest] = fcx;
	      nestPos[3 * nest + 1] = fcy;
		nestPos[3 * nest + 2] = fcr;
		//printf("Nest %d - radius %lf\n", nest, fcr);
		nthrobots += calcRobots(fcr);
		// Update minimum size
		if (nest == (nnests - 2))
		{
			//printf("nthrobots %d\n", nthrobots);
			// Check the number of robots
			if (nthrobots < nrobots)
			{
				nleftrobots = (nrobots - nthrobots);
				//printf("nleftrobots %d\n", nleftrobots);
				if (nleftrobots > 2)
				{
					minsize = calcNestRadius(nleftrobots);
					//printf("new minimum radius %lf\n", minsize);
				}
			}
		}
	    }										
	}
    // initial positions and orientations of the robots
    // robots are placed far from the nests
    for (r1=0, ro1=rob; r1 < nrobots; r1++, ro1++)
    {
      placed = false;
      //attempts = 0;
      while (!placed)// && attempts < 100)
      {
        ro1->dir = RNG->getDouble(0.0, PI2);
        ro1->x = RNG->getDouble(200.0, wsize - 200.0);
        ro1->y = RNG->getDouble(200.0, wsize - 200.0);
        placed = true;
        // Check whether robot position overlaps with nests
        nest = 0;
        while ((nest < nnests) && placed)
        {
          // Distance from nest
          //printf("robot %d - nest %d\n", r1, nest);
          dx = (ro1->x - nestPos[3 * nest]);
          dy = (ro1->y - nestPos[3 * nest + 1]);
          dist = sqrt(dx * dx + dy * dy);
          if (dist < robotnestdist)
            placed = false;
          nest++;
        }
        if (placed)
        {
          // Check whether robot position overlaps with other robots
          for (r2=0, ro2=rob; r2 < r1; r2++, ro2++)
          {
            dx = (ro1->x - ro2->x);
            dy = (ro1->y - ro2->y);
            dist = sqrt(dx*dx+dy*dy);
            if (dist < robotdist)//(ro1->radius + ro2->radius))
              placed = false;
          }
        }
      }
      // Set initial state of the robot
      *is = ro1->dir; is++;
      *is = ro1->x; is++;
      *is = ro1->y; is++;
    }
  }

}

/*
 *  update the power-distances between robots and the matrix of nearest robots
 */
void updateRobotDistances()

{

    //double *robotsdist;     // matrix containing robot distances;
    //int *robotsbydistance;  // matrix containing the id of the robots ordered by the inverse of the distance

    struct robot *ro1;
    struct robot *ro2;
    int r1, r2, r3;
    double *cdist;
    double *ccdist;
    double smallerdist;
    int smallerid;
    int *rbd;
    double *rbddist;
    int remaining;


    // update the matrix of distances
    cdist = robotsdist;
    for (r1=0, ro1=rob; r1 < nrobots; r1++, ro1++)
    {
      for (r2=0, ro2=rob; r2 < r1; r2++, ro2++)
      {
        *cdist = ((ro1->x - ro2->x)*(ro1->x - ro2->x)) + ((ro1->y - ro2->y)*(ro1->y - ro2->y));
        *(robotsdist + (r1 + (r2 * nrobots))) = *cdist;
        cdist++;
      }
      *cdist = 9999999999.0;; // distance from itself is set to a large number to exclude itself from the nearest
      cdist = (cdist + (nrobots - r1));
    }

    // update the matrix of nearest robots
    rbd = robotsbydistance;
    for (r1=0, cdist = robotsdist; r1 < nrobots; r1++, cdist = (cdist + nrobots))
    {
      remaining = (nrobots - 1);
      for (r2=0; r2 < (nrobots - 1); r2++)
      {
        for (r3=0, ccdist = cdist, smallerid=0, smallerdist = *ccdist, rbddist = ccdist; r3 < nrobots; r3++, ccdist++)
        {
          if (*ccdist <= smallerdist)
          {
            smallerdist = *ccdist;
            smallerid = r3;
            rbddist = ccdist;
          }
        }
        if (smallerdist < 1000000)  // we ignore robots located at a power dstance greater that 750*750
         {
          *rbd = smallerid;
          rbd++;
          *rbddist = 9999999999.0;
          remaining--;
         }
         else
         {
          *rbd = -1;
          rbd = (rbd + remaining);
          break;
         }
      }
      *rbd = -1;  // we use -1 to indicate the end of the list
      rbd++;

    }

    /*
    cdist = robotsdist;
    for (r1=0; r1 < nrobots; r1++)
        {
          for (r2=0; r2 < nrobots; r2++)
             {
               printf("%10.1f ", *cdist);
               cdist++;
             }
             printf("\n");
         }
     printf("\n");
     */

    /*
    rbd = robotsbydistance;
    for (r1=0; r1 < nrobots; r1++)
        {
          for (r2=0; r2 < nrobots; r2++)
             {
               printf("%3d ", *rbd);
               rbd++;
             }
             printf("\n");
         }
    printf("\n");
    */

}



/*
 * initialize the environment, the robots, and the network
 * if mode != 0 regenerate the trial initial conditions
 * return the number of required parameters (-1 in case of failure)
 */
int initialize(const char* filename)

{
    
    struct robot *ro;
    struct evonet *ne;
	int r;


	// default values
    trials = 1;
    ttrials = 1;
    steps = 100;
	
	// set USA local conventions
	setlocale( LC_ALL, "en-US" );
	// read parameters from the .ini file
	if (readEvoConfig(filename) == false)
		return(-1);

    // create and configure the environment
    initEnvironment();

    // create the robots, the networks, and initialize the sensors
    rob = (struct robot *) malloc(nrobots * sizeof(struct robot));
	net = (struct evonet *) malloc(nrobots * sizeof(struct evonet));
    // initialize the robots and the networks
	for (r=0, ro=rob, ne=net; r < nrobots; r++, ro++, ne++)
	 {
	  // initialize the robots
      initRobot(ro, r, robottype);
	  // set maxspeed
	  switch (robottype)
	   {
	   case (Khepera):
	    ro->maxSpeed = 144.0;
		break;
	   case (ePuck):
        ro->maxSpeed = 200.0;
		break;
	   case (MarXBot):
        ro->maxSpeed = 500.0;
		break;
	   }
      ro->color = 4;  // frontal and rear side are red and blue, respectively
      // initialize sensors
      ne->ninputs = 0;
      ne->ninputs += initInfraredSensor(ro, 8, ne->ninputs);
      if (camera)
        ne->ninputs += initCameraSensorRFB(ro, 2);
      ne->ninputs += initGroundSensor(ro);
      //ne->ninputs += initEnergySensor(ro);
      //ne->ninputs += initBiasSensor(ro);
      // initialize motors
      ro->motorwheels = 2;
      ro->motorwheelstype = 1; // encoding of speed and rotation
      ro->motorleds = 0;
      if (camera)
        ro->motorleds = 2;
      ne->noutputs = ro->motorwheels + ro->motorleds;
	    // allocate/initialize the input vector
	    initRobotSensors(ro, ne->ninputs);
      // initialize robot's network
      initNet(ne, filename, r);
      resetNet(ne);
      // set the id number of motors neurons
      ro->motorwheelsid = ne->ninputs + ne->nhiddens;
      ro->motorledsid = ne->ninputs + ne->nhiddens + ro->motorwheels;
	 }

	// Update target areas (nests) size
	asize = calcNestRadius(nrobots);//sqrt(rob->radius * rob->radius * nrobots);//(rob->radius * 2.0 * nrobots) / (2.0 * M_PI); // We update the nest size based on 
  //printf("asize = %lf\n", asize);
  //exit(-1);
	// Compute the minimum nest size
	msize = (rob->radius * 2.0);//rob->radius;//(rob->radius * 2.0) / (2.0 * M_PI); // The minimum size to have a robot in a nest
	//printf("msize %lf asize %lf\n", msize, asize);
	//exit(-1);
	
	// generate initial states
    // (initial positions/orientations of robots and positions and radii of nests)
    initialstates = (double *) malloc(((trials+ttrials)*(3*nrobots+3*nnests)) * sizeof(double));
    createInitialStates(0);

    // allocate distance matrices and initialize the dyagonal distances to 0
    robotsdist = (double *) malloc((nrobots*nrobots) * sizeof(double));
    robotsbydistance = (int *) malloc((nrobots*nrobots) * sizeof(int));

	nfreep = computeParameters(net);
	printf("Network parameters %d\n", nfreep);
	
	return(nfreep);
}

/*
 * evaluate a genotype for capability to perform the discrimination task
 * return the total fitness
 * when mode=1 evaluate a genotype on the test trials
 * seed is used in graphic mode to store data for statistical analysis and plotting
 * when its value is >= 0
 */
double evaluate(double* genotype, const int glen, int mode, int seed)

{
	
	  int r,f,s;
    double *is;
    int trial;
	  int numtrials;
    int step;
    double fitness;
	  struct robot *ro, *ro2;
	  struct evonet *ne;
	  bool stoptrial;
    double oldx, oldy; // updated every 2s
	  double dx, dy;
    int foodarea[100];
	  int collectedfood[100];
	  double cdist;
    int *rbd;
    int x, y;
	
    //double* fit = (double*)malloc(sizeof(double)*nrobots);
    double fit[100];
    double tfit;
    int nest;
    //const double nestBonus = 100.0; // To be encoded
    bool robotOnNest[100];
    bool onNest;
    // Number of robots on a given nest (at the end of each trial)
    int nrobotsOnNest[100];
    // Number of steps each robot is not on a nest
    int nexploresteps[100];
    int rr;
    const double maxDist = sqrt(2.0*wsize*wsize);
    double robotradius;
printf("QUI\n");	
	  cellsize = worldx / 20.0;
	
    // evaluation loop
	  if (glen != nfreep)
	  {
	    printf("WARNING: The length of the parameter vector received by discriminate (%d) does not match the number of parameters required (%d)\n", glen, nfreep);
	  }
    fitness = 0.0;
	  if (mode == 0)
	  {
	    is = initialstates;
		  numtrials = trials;
		}
	  else
	  {
	    is = (initialstates + (trials * ((3 * nrobots) + (3 * nnests))));
		  numtrials = ttrials;
	  }
printf("QUI0\n");
	  for (trial=0; trial < numtrials; trial++)
    {
      // the robots initially do not carry food
		  for(r=0; r < nrobots; r++)
		  {
        foodarea[r] = 0;
		    collectedfood[r] = 0;
        fit[r] = 0.0; // Reset fitness
        robotOnNest[r] = false;
		  }
		  stoptrial = false;

        for(x=0; x < 100; x++)
          for(y=0; y < 100; y++)
            cells[x][y] = 0;
		
	    for (r=0, ro=rob, ne=net; r < nrobots; r++, ro++, ne++)
		  {
           resetNet(ne);
			for(s=0,ro->csensors=ro->sensors; s < ne->ninputs; s++, ro->csensors++)
		      *ro->csensors = 0.0;
		  }
        // initialize the nest position and radius
        for (f=0; f < nnests; f++)
        {
          envobjs[f+4].x = *is; is++;
          envobjs[f+4].y = *is; is++;
          envobjs[f+4].r = *is; is++;
          //printf("%d\t%lf\n", f, envobjs[f+4].r);
          nrobotsOnNest[f] = 0;
        }
printf("QUI1\n");
	// Save initial position of the robots
	char rfilename[256];
    	char rbasename[256];
    	char rsuffix[64];
    	FILE* rfp;
	// Build filename
      if (seed > 0)
      {
      	sprintf(rbasename, "robotPosS%d", seed);
      	strcpy(rfilename, rbasename);
      	sprintf(rsuffix, "T%d.txt", trial + 1);
      	strcat(rfilename, rsuffix);
      	// Open file
      	rfp = fopen(rfilename, "a");
      }
        // initialize robots position, orientation, and energy level
	      for (r=0, ro=rob; r < nrobots; r++, ro++)
	      {
          ro->alive = true;
		      ro->dir = *is; is++;
	        ro->x = *is; is++;
	        ro->y = *is; is++;
		      ro->energy = 1.0;
          nexploresteps[r] = 0;
          if (r == 0)
            robotradius = ro->radius;
            if (seed > 0)
            {
              	if (rfp != NULL)
              	{
			            fprintf(rfp, "%lf\t%lf\t", ro->x, ro->y);
	      	      }
            }
	}
  if (seed > 0)
      {
		//printf("trial %d robot1 %d %.2f %.2f %.2f robot2 %d %.2f %.2f %.2f \n", trial, rob->dir, rob->x, rob->y, b1, b2, b3);
	      fprintf(rfp, "\n");
	      fclose(rfp);
      }
printf("QUI2\n");
        for (step=0; step < steps; step++)
        {
          if (seed > 0)
		        rfp = fopen(rfilename, "a");
          // Reset number of robots on each nest (we compute the number at each step to avoid errors)
          for (nest = 0; nest < nnests; nest++)
            nrobotsOnNest[nest] = 0;
          for (r=0, ro=rob, ne=net, rbd = robotsbydistance; r < nrobots; r++, ro++, ne++, rbd = (rbd + nrobots))
		      {
            onNest = false;
            // update robots distances
            updateRobotDistances();
			      // update sensors
            ro->csensors = ro->sensors;
            updateInfraredSensor(ro);
            if (camera)
              updateCameraSensorRFB(ro, rbd);
			      updateGroundSensor(ro);
            //updateEnergySensor(ro);
		      	//updateBiasSensor(ro);
			      // update network
			      updateNet(ne, ro->sensors, genotype);
			      // update the robot
			      oldx = ro->x;
			      oldy = ro->y;
			      updateRobot(ro, ne);
            for (nest = 0; nest < nnests; nest++)
            {
              dx = ro->x - envobjs[4 + nest].x;
              dy = ro->y - envobjs[4 + nest].y;
              cdist = sqrt((dx*dx)+(dy*dy));
			        if (cdist < envobjs[4 + nest].r)
              {
                onNest = true;
                // Update the number of robots on that nest
                nrobotsOnNest[nest]++;
              }

            }
            if (onNest)
              robotOnNest[r] = true;
            else
            {
              robotOnNest[r] = false;
              // Update the number of steps the robot is "exploring" (i.e., it is not on a nest)
              nexploresteps[r]++;
            }
printf("QUI3\n");
      if (seed > 0)
          {      
		        if (rfp != NULL)
              	{
			      fprintf(rfp, "%lf\t%lf\t", ro->x, ro->y);
	      	      }
          }
			    }
          if (seed > 0)
          {
	          fprintf(rfp, "\n");
	          fclose(rfp);
          }	

          #ifdef __cplusplus
            if (renderrobotworld) update_world();		// update graphic rendering
          #endif

          if (seed >= 0)
          {
            // Store data for graphic purposes
            char filename[256];
            char basename[256];
            char suffix[64];
            FILE* fp;
            if (!fixedpos)
            {
              // Build filename
              sprintf(basename, "metricsS%d", seed);
              strcpy(filename, basename);
              sprintf(suffix, "T%d.txt", trial + 1);
              strcat(filename, suffix);
              // Open file
              fp = fopen(filename, "a");
              if (fp != NULL)
              {
                // Compute cluster and dispersion metrics
                // Cluster metric
                double cluster = 0.0;
                int nclusters = 0;
                for (nest = 0; nest < nnests; nest++)
                {
                  double nmax;
                  double circ;
                  if (nrobotsOnNest[nest] > 0)
                  {
                    circ = 2.0 * M_PI * envobjs[4 + nest].r;
                    nmax = circ / (2.0 * robotradius);//pow(envobjs[4 + nest].r, 2.0) / pow(robotradius, 2.0);// (circ / (2.0 * M_PI * robotradius));
                    cluster += (nrobotsOnNest[nest] / nmax);
                    nclusters++;
                  }
                }
                if (nclusters > 0)
                  cluster /= nclusters;
                // Dispersion metric
                double dispersion = 0.0;
                // Compute COG (Center Of Gravity)
                double cogx = 0.0;
                double cogy = 0.0;
                for (r=0, ro=rob; r < nrobots; r++, ro++)
                {
                  cogx += ro->x;
                  cogy += ro->y;
                }
                cogx /= nrobots;
                cogy /= nrobots;
                // Now compute dispersion metric
                double dist;
                for (r=0, ro=rob; r < nrobots; r++, ro++)
                {
                  dx = (cogx - ro->x);
                  dy = (cogy - ro->y);
                  dist = sqrt(dx * dx + dy * dy);
                  dispersion += (dist * dist);
                }
                dispersion /= (4.0 * robotradius * robotradius);
                fprintf(fp, "%lf\t%lf\n", cluster, dispersion);
                fclose(fp);
              }
            }
            else
            {
              // Now store nest positions and robot positions
              // Build filename
              sprintf(basename, "fixedPositionsS%d", seed);
              strcpy(filename, basename);
              sprintf(suffix, "T%d.txt", trial + 1);
              strcat(filename, suffix);
              // Open file
              fp = fopen(filename, "a");
              if (fp != NULL)
              {
                for (r=0, ro=rob; r < nrobots; r++, ro++)
                {
                  // Robot position
                  fprintf(fp, "%lf\t%lf\t", ro->x, ro->y);
                }
                for (nest = 0; nest < nnests; nest++)
                {
                  // Nest position and radius
                  fprintf(fp, "%lf\t%lf\t%lf\t", envobjs[nest+4].x, envobjs[nest+4].y, envobjs[nest+4].r);
                }
                fprintf(fp, "\n");
                fclose(fp);
              }
            }
          }
        }
printf("QUI4\n");

      // Compute episode fitness
      tfit = 0.0;
      for (r=0, ro=rob; r < nrobots; r++, ro++)
      {
        fit[r] = 0.0;
        if (robotOnNest[r])
        {
          fit[r] += nestBonus; // Bonus
        }
          // Compute average distance with swarm-mates
          cdist = 0.0;
          for (rr = 0, ro2 = rob; rr < nrobots; rr++, ro2++)
          {
            // Distance with others
            dx = (ro->x - ro2->x);
            dy = (ro->y - ro2->y);
            double mydist = sqrt((dx*dx)+(dy*dy));
            //printf("Trial %d: Dist %d - %d: %lf\n", trial, r, rr, mydist);
            cdist += (mydist / maxDist);
          }
          cdist /= (nrobots - 1);
          // Add component rewarding for swarm distance minimization
          fit[r] += kdist * (1.0 - cdist);
        //}
        // Update episode fitness
        tfit += fit[r];
      }
      // Normalize episode fitness over the number of robots
      tfit /= nrobots;
      // Update fitness
      fitness += tfit;
      if (seed >= 0)
      {
        // Store data for graphic purposes
        char filename[256];
        char basename[256];
        char suffix[64];
        FILE* fp;
        double cogx, cogy;
        int nout;
        double outxy[100];
        int idx;
        double avgd;
        int ndists;
        int cwall;
        if (!fixedpos)
        {
          // Build filename
          sprintf(basename, "positionsS%d", seed);
          strcpy(filename, basename);
          sprintf(suffix, "T%d.txt", trial + 1);
          strcat(filename, suffix);
          // Open file
          fp = fopen(filename, "w");
          if (fp != NULL)
          {
            // Data to be written in the following order:
            // 1) position of the walls
            // 2) position and radius of each nest + number of robots on the nest (4 * nnests)
            // 3) position and radius of each robot + number of steps the robot is not on a nest (4 * nrobots)
            // 4) center of gravity of the swarm (2)
            // 5) number and average distance of robots out of the nests
            for (cwall = 0; cwall < 4; cwall++)
            {
              // Wall position
              fprintf(fp, "%lf\t%lf\t%lf\t%lf\t", envobjs[cwall].x, envobjs[cwall].y, envobjs[cwall].x2, envobjs[cwall].y2);
            }
            for (nest = 0; nest < nnests; nest++)
            {
              // Nest position and radius
              fprintf(fp, "%lf\t%lf\t%lf\t%d\t", envobjs[nest+4].x, envobjs[nest+4].y, envobjs[nest+4].r, nrobotsOnNest[nest]);
            }
            cogx = 0.0;
            cogy = 0.0;
            for (r=0, ro=rob; r < nrobots; r++, ro++)
            {
              // Robot position and radius
              fprintf(fp, "%lf\t%lf\t%lf\t%d\t", ro->x, ro->y, ro->radius, nexploresteps[r]);
              // Compute the center of gravity of the swarm
              cogx += ro->x;
              cogy += ro->y;
            }
            cogx /= nrobots;
            cogy /= nrobots;
            // Center of gravity
            fprintf(fp, "%lf\t%lf\t", cogx, cogy);
            // Compute number of robots out of the nests
            nout = 0;
            idx = 0;
            for (r=0, ro=rob; r < nrobots; r++, ro++)
            {
              if (!robotOnNest[r])
              {
                nout++;
                outxy[idx] = ro->x;
                idx++;
                outxy[idx] = ro->y;
                idx++;
              }
            }
            fprintf(fp, "%d\t", nout);
            // Compute average distance of robots out of the nests
            avgd = 0.0;
            if (nout > 1)
            {
              ndists = 0;
              for (r = 0; r < nout; r++)
              {
                for (rr = 0; rr < r; rr++)
                {
                  dx = (outxy[2 * r] - outxy[2 * rr]);
                  dy = (outxy[2 * r + 1] - outxy[2 * rr + 1]);
                  avgd += sqrt((dx*dx)+(dy*dy));
                  ndists++;
                }
              }
              avgd /= ndists;
            }
            fprintf(fp, "%lf\t%lf\n", avgd, tfit);
            fclose(fp);
          }
        }
      }
	}
	
    return(fitness / numtrials);
    
}


