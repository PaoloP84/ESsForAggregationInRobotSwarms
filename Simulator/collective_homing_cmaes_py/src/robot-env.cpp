/*
 * robot-env.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <locale.h>
#include "robot-env.h"
#include "evonet.h"
#include "utilities.h"

extern int trials;						// number of trials
extern int ttrials;					// number of testing trials

// global variables
struct robot *rob;					// list of robots
int nrobots;						// number of robots
int nnests;                 // number of nests
double nestdist;            // minimum nest distance
double robotnestdist;       // minimum distance between robot and nest
double robotdist;           // minimum robot distance
struct evonet *net;					// list of networks
struct envobjects *envobjs;			// list of environemntal objects
int nenvobjs;						// number of objects
double worldx;						// world x dimension
double worldy;						// world y dimension
double *initialstates;				// robots initial positions
int cells[100][100];                // virtual cells
double cellsize=0.0;                // virtual cell size
char globdir[512];

// object samples
float **wall = NULL;				// wall samples
int *wallcf = NULL;					// description of wall samples
float **cylinder = NULL;			// cylindrical sample
int *cylindercf = NULL;				// description of cylindrical sample
float **scylinder = NULL;			// small cylindrical sample
int *scylindercf = NULL;			// description of small cylindrical sample
float **light = NULL;				// light samples
int *lightcf = NULL;				// description of light samples
int **iwall = NULL;                 // wall samples  (old style)
int **iobst = NULL;					// cylindrical sample (old style)
int **isobst = NULL;				// small cylindrical sample (old style)
int **ilight = NULL;				// light samples (old style)

int steps = 100;					// number of steps x trial
int nfreep = 0;						// number of free parameters
int renderrobotworld=0;             // whether we update the graphic renderer


/*
 * load an obstacle from a .sam file
 */
float **
load_obstacle(char *filename, int  *objectcf)

{

   float **object;
   float  **ob;
   float  *o;
   int    *ocf;
   FILE *fp;
   int  i,ii,t,v,s;
   char sbuffer[128];

   if ((fp = fopen(filename,"r")) == NULL)
	 {
	    printf("I cannot open file %s", filename);
		fflush(stdout);
	  }
   // read configuration data: 0-nsensors, 1-nangles, 2-ndistances,
   //                          3-initdist,4-distinterval
   for (i=0, ocf=objectcf; i < 5; i++, ocf++)
     fscanf(fp,"%d ",ocf);
   fscanf(fp,"\n");
	
   // allocate space and initialize
   object = (float **) malloc(objectcf[2] * sizeof(float *));
   ob = object;
   for (i=0; i<objectcf[2]; i++,ob++)
	  {
	   *ob = (float *)malloc((objectcf[0] * objectcf[1]) * sizeof(float));
	   for (ii=0,o=*ob; ii < (objectcf[0] * objectcf[1]); ii++)
		*o = 0;
	  }
   // load data
   for (t=0, ob=object; t < objectcf[2]; t++, ob++)
    {
     fscanf(fp,"TURN %d\n",&v);
     for (i = 0, o = *ob; i < objectcf[1]; i++)
      {
       for (s=0; s < objectcf[0]; s++,o++)
        {
	   fscanf(fp,"%f ",o);
        }
       fscanf(fp,"\n");
      }
    }
   // sanity check
   fscanf(fp, "%s\n",sbuffer);
   if (strcmp(sbuffer,"END") != 0)
     {
      printf("loading file %s: the number of sample is inconsistent", filename);
	  fflush(stdout);
	 }
	
   fclose(fp);
   return object;
}


/*
 *	initialize the infrared sensor
 *  load tha sample files
 */
int initInfraredSensor(struct robot *cro, int nsensors, int ids)
{
   char filename[512];
   cro->nifsensors = nsensors;
   cro->sensorinfraredn = nsensors;
   cro->sensorinfraredid = ids;

   if (cro->idn == 0) printf("Sensor[%d]: sampled infrared sensors (robot %d)\n", cro->nifsensors, cro->type);

   // init the environment and load environmental samples
   switch (cro->type)
    {
	  case Khepera:
		wallcf = (int *) malloc(5 * sizeof(int));
	strcpy(filename, globdir);
	strcat(filename, "khepera-wall.sample");
        wall = load_obstacle(filename, wallcf);
        scylindercf = (int *) malloc(5 * sizeof(int));
	strcpy(filename, globdir);
	strcat(filename, "khepera-scylinder.sample");
        scylinder = load_obstacle(filename, scylindercf);
        cylindercf = (int *) malloc(5 * sizeof(int));
	strcpy(filename, globdir);
	strcat(filename, "khepera-cylinder.sample");
        cylinder = load_obstacle(filename, cylindercf);
	  break;
	  case ePuck: // TO DO: cylider should be resampled (we are using scylinder)
		wallcf = (int *) malloc(5 * sizeof(int));
	strcpy(filename, globdir);
	strcat(filename, "epuck-wall.sample");
        wall = load_obstacle(filename, wallcf);
        scylindercf = (int *) malloc(5 * sizeof(int));
	strcpy(filename, globdir);
	strcat(filename, "epuck-scylinder.sample");
        scylinder = load_obstacle(filename, scylindercf);
        cylindercf = (int *) malloc(5 * sizeof(int));
	strcpy(filename, globdir);
	strcat(filename, "epuck-cylinder.sample");
        cylinder = load_obstacle(filename, cylindercf);
	  break;
	  case MarXBot:
		wallcf = (int *) malloc(5 * sizeof(int));
	strcpy(filename, globdir);
	strcat(filename, "marxbot-wall.sample");
        wall = load_obstacle(filename, wallcf);
        scylindercf = (int *) malloc(5 * sizeof(int));
	strcpy(filename, globdir);
	strcat(filename, "marxbot-scylinder.sample");
        scylinder = load_obstacle(filename, scylindercf);
        cylindercf = (int *) malloc(5 * sizeof(int));
	strcpy(filename, globdir);
	strcat(filename, "marxbot-cylinder.sample");
        cylinder = load_obstacle(filename, cylindercf);
	  break;
   }

   return(cro->nifsensors);
}

/*
 *	update the infrared sensor
 *  that activation produced by multiple objects in the perception range is summed
 */
void updateInfraredSensor(struct robot *cro)

{

	int ob;
	double x, y, dx, dy, dist, angle;
	int relang, idist;
	float *o;
	int s,r;
	struct robot *ro;
	float act[8];
	
	// initialize
	for(s=0; s < cro->nifsensors; s++)
	    act[s] = 0.0;
	// calculate the stimolation generated environmental objects
	for(ob=0; ob < nenvobjs; ob++)
	{
		switch(envobjs[ob].type)
		  {
			// sampled walls (we assume they are orthogonal with respect to the carthesian axes)
			case WALL:
				if (cro->x >= envobjs[ob].x && cro->x <= envobjs[ob].x2)
					x = cro->x;
				 else
					x = envobjs[ob].x;
				if (cro->y >= envobjs[ob].y && cro->y <= envobjs[ob].y2)
					y = cro->y;
				 else
					y = envobjs[ob].y;
				dx = (cro->x - x);
				dy = (cro->y - y);
				dist = sqrt(dx*dx+dy*dy);
				if (dist < (wallcf[3] + cro->radius + (wallcf[2] * wallcf[4])))
				 {
					angle = angv(x, y, cro->x, cro->y);
					angle = mangrelr(angle, cro->dir);
				    if (angle > (M_PI * 2.0)) angle -= (M_PI * 2.0);
					relang  = (int) round(angle / PI2 * 180.0);
					if (relang == 180) relang = 0;
				   //if (x == cro->x)
					//  if (y < cro->y) angle = M_PI + (M_PI * 2.0); else angle = M_PI / 2.0;
				   //else
					// if (x < cro->x) angle = M_PI; else angle = 0.0;
				   //relang  = (int) round((mangrelr(angle,cro->dir) / (M_PI * 2.0)) * 180.0);
				   idist = (int) round((dist - cro->radius - wallcf[3]) / wallcf[4]);
				   if (idist < 0)
					  idist = 0;
				   if (idist >= wallcf[2])
				      idist = wallcf[2]-1;
				   o = *(wall + idist);
				   o = (o + (relang * wallcf[0]));
				   for(s=0;s<wallcf[0];s++,o++)
                      {
						act[s] += *o;
                      }
                 }
			break;
			// sampled cylinder
			case SAMPLEDSCYLINDER:
				x = envobjs[ob].x;
				y = envobjs[ob].y;
				dx = (cro->x - x);
				dy = (cro->y - y);
				dist = sqrt(dx*dx+dy*dy) - (cro->radius + envobjs[ob].r - scylindercf[3]);
				//dot = cro->x*x + cro->y*y;
				//det = cro->x*y - cro->y*x;
				//ang = atan2(det, dot) - cro->dir;
				angle = angv(x,y, cro->x, cro->y);
				angle = mangrelr(angle, cro->dir);
				relang  = (int) round(angle / (M_PI * 2.0) * 180.0);
				if (relang == 180) relang = 0;
				if (dist < (scylindercf[2] * scylindercf[4]))
				 {
				   idist = (int) round(dist / scylindercf[4]);
				   if (idist < 0)
					  idist = 0;
				   if (idist >= scylindercf[2])
				      idist = scylindercf[2]-1;
				   o = *(scylinder + idist);
				   o = (o + (relang * wallcf[0]));
				   for(s=0;s<wallcf[0];s++,o++)
                      {
						act[s] += *o;
                      }
                 }
			break;
			// default
			//default:
			//	printf("ERROR: undefined object type %d \n",envobjs[ob].type);
			break;
		 }
	}
	// calculate the stimulation generated by other croots (treated as cylindrical objects)
	if (nrobots > 1)
	  {
	    for (r=0, ro=rob; r < nrobots; r++, ro++)
	      {
		    if (cro->idn != ro->idn)
			 {
				dx = (cro->x - ro->x);
				dy = (cro->y - ro->y);
				dist = sqrt(dx*dx+dy*dy) - (cro->radius * 2.0);
				angle = angv(ro->x, ro->y, cro->x, cro->y);
				angle = mangrelr(angle, cro->dir);
				if (angle > PI2) angle -= PI2;
				relang  = (int) round(angle / PI2 * 180.0);
				if (relang == 180) relang = 0;
				//printf("dist %.2f soglia %d ang %.2f dir %.2f relang %d \n", dist, cylindercf[2] * cylindercf[4], ang, cro->dir, relang);
				if (dist < (cylindercf[2] * cylindercf[4]))
				 {
				   idist = (int) round(dist / wallcf[4]);
				   if (idist < 0)
					  idist = 0;
				   if (idist >= cylindercf[2])
				      idist = cylindercf[2]-1;
				   o = *(cylinder + idist);
				   o = (o + (relang * wallcf[0]));
				   for(s=0;s<cylindercf[0];s++,o++)
					 act[s] += *o;
                 }
			  }
		   }
		}
	
	// normalize and copy on the imput vector
	for(s=0; s < cro->nifsensors; s++)
		{
		 if (act[s] > 1.0)
             act[s] = 1.0;
         *cro->csensors = act[s];
         cro->csensors++;
		}
}


/*
 *	initialize the PP camera sensor
 */
int initCameraPPSensor(struct robot *cro)
{
   if (cro->idn == 0) printf("Camera[%d]: single color, 8 sectors and distance \n", 9);
   return(9);
}


/*
 *	update the camera sensor
 *  assume that the environment include a single coloured object constituted by the other robot
 *  encode angular centre of the color blob produced by the other robot withon 8 sectors and the average fraction of pixel stimulated
 */
void updateCameraPPSensor(struct robot *cro)

{

    double dx, dy, dist, angle;
    double angleb;
    int n,r,s;
    struct robot *ro;
    double act[8];
    double adist;
    double sectoracenter;
    double pi8;
    double pi16;

    pi8 = PI2 / 8.0;
    pi16 = PI2 / 16.0;
    for (n=0; n < 8; n++)
        act[n] = 0.0;
    for (r=0, ro=rob; r < nrobots; r++, ro++)
       {
         if (cro->idn != ro->idn)
           {
             dx = (cro->x - ro->x);
             dy = (cro->y - ro->y);
             dist = sqrt(dx*dx+dy*dy);
             angle = mangrelr(angv(ro->x, ro->y, cro->x, cro->y), cro->dir);
             for (s=0, sectoracenter = 0.0 + (PI2 / 16.0); s < 8; s++, sectoracenter += pi8)
             {
               angleb = angle;
               if (s == 0 && angle > M_PI)
                 angleb = angle - PI2;
               if (s == 7 && angle < M_PI)
                 angleb = angle + PI2;
               if (fabs(angleb - sectoracenter) < pi8)
                 {
                  act[s] = 1.0 - ((fabs(angleb - sectoracenter) / pi8));
                  //printf("s %d ang %.2f secta %.2f fabs %.2f -> %.2f\n", s, angle, sectoracenter, fabs(angleb - sectoracenter), act[s]);
                 }
             }

           }
       }
    // perceived angular blob (maximum value il 0.93 radiants for two adjacent marxbots)
    adist = ((M_PI / 2.0) - atan(dist/cro->radius)) * 2.0;
    //printf("dist %.2f angular blob %.2f\n", dist, adist);

    // normalize and copy on the imput vector
    for(n=0; n < 8; n++)
        {
         *cro->csensors = act[n];
         //printf("%.2f\n", act[n]);
         cro->csensors++;
        }
    *cro->csensors = adist;
    //printf("%.2f\n", adist);
    cro->csensors++;
    //printf("end\n");
}

/*
 *	initialize the ground gradient sensor
 */
int initGroundGradSensor(struct robot *cro)
{
   return(5);
}


/*
 *	update the marXbot ground gradient sensors
 *  assume that the color of the area is white at the centre and progressively darker in the periphery
 *  the last value encode the average color detected by the four sensors located at the frontal-left/right and at the back-left/right
 *  the first four values variation of the with respect to the average color
 */
void updateGroundGradSensor(struct robot *cro)

{
    double act[5];
    double gcolor[4];
    double x,y;
    double dx, dy;
    double average;
    double wx, wy;
    double maxdist;
    int s, n;
    double sensordist;

    wx = worldx / 2.0;
    wy = worldy / 2.0;
    maxdist = sqrt((wx*wx)+(wy*wy));
    sensordist = cro->radius * 0.9;

    //front left
    x = cro->x + xvect(cro->dir - 0.0872 , sensordist);
    y = cro->y + yvect(cro->dir - 0.0872 , sensordist);
    dx = x - wx;
    dy = y - wy;
    gcolor[0] = (sqrt((dx*dx)+(dy*dy)) / maxdist);
    //front right
    x = cro->x + xvect(cro->dir + 0.0872 , sensordist);
    y = cro->y + yvect(cro->dir + 0.0872 , sensordist);
    dx = x - wx;
    dy = y - wy;
    gcolor[1] = (sqrt((dx*dx)+(dy*dy)) / maxdist);
    //rear left
    x = cro->x + xvect(cro->dir - 0.0872 + M_PI, sensordist);
    y = cro->y + yvect(cro->dir - 0.0872 + M_PI, sensordist);
    dx = x - wx;
    dy = y - wy;
    gcolor[2] = (sqrt((dx*dx)+(dy*dy)) / maxdist);
    //rear right
    x = cro->x + xvect(cro->dir + 0.0872 + M_PI, sensordist);
    y = cro->y + yvect(cro->dir + 0.0872 + M_PI, sensordist);
    dx = x - wx;
    dy = y - wy;
    gcolor[3] = (sqrt((dx*dx)+(dy*dy)) / maxdist);

    // average value
    act[4] = average = (gcolor[0] + gcolor[1] + gcolor[2] + gcolor[3]) / 4.0;
    // variations
    for(s=0; s < 4; s++)
      act[s] = 0.5 + (gcolor[s] - average) * 10.0;

    // normalize and copy on the imput vector
    for(n=0; n < 5; n++)
        {
         *cro->csensors = act[n];
         cro->csensors++;
        }


}


/*
 * verify whether the angle a is in the range [r1,r2]
 * assume that a,r1 and r2 are in the range [0, PI2]
 */
bool anginside(double a, double r1, double r2)

{

    if ((r2 - r1) > 0)
    {
      if ((r2 - r1) < M_PI)
        {
         if (a > r1 && a < r2) // clockwise from r1 to r2
           return(true);
          else
           return(false);
        }
         else
        {
         if (a > r2 || a < r1) // clockwise from r2 to PI2 and from 0 to r1
           return(true);
          else
           return(false);
         }
    }
    else
    {
        if ((r1 - r2) < M_PI)
          {
           if (a > r2 && a < r1) // counter-clockwise from r2 to r1
             return(true);
            else
             return(false);
          }
           else
          {
           if (a > r1 || a < r2) // couner-clockwise from r1 to PI2 and from 0 to r2
             return(true);
            else
             return(false);
           }
    }


}

/*
 *	update the laser sensor
 *  measure the closest distance detected in each perceptual sector
 */
int initLaserDistanceSensor(struct robot *cro)
{
   if (cro->nlasersensors <= 0 || cro->nlasersensors > 8)
     cro->nlasersensors = 8;
   printf("Sensor[%d]: laser sensors\n", cro->nlasersensors);
   return(cro->nlasersensors);
}


void updateLaserDistanceSensor(struct robot *cro)

{

    double maxdist = 100;          // maximum dist at which  objects are detected
    double starta = -M_PI;          // initial angle of the first sectonr
    double inta = M_PI / (double) (cro->nlasersensors / 2);
    int ob;
    double x, y, dx, dy, dist;
    double a1, a2;
    double angle;                   // angle of cylindrical object or of the nearest point of a wall
    double angle1, angle2;          // angle of the initial and final point of a wall
    int s,r;
    struct robot *ro;
    double mdist[8];                // minimum detected distance
    double cdist;
    double ldist;                   // lateral distance
    double angoffset;               // angular offset between the angle of the perceived object and the angular center of the sector
    double maxperceivableaoffset;   // the maximum angular offset of a perceived object with respect to a sector

    // initialize
    for(s=0; s < cro->nlasersensors; s++)
        mdist[s] = maxdist;

    // calculate the stimulation generated environmental objects
    for(ob=0; ob < nenvobjs; ob++)
    {
        switch(envobjs[ob].type)
          {
            // walls (we assume they are orthogonal with respect to the carthesian axes)
            case WALL:
                if (cro->x >= envobjs[ob].x && cro->x <= envobjs[ob].x2)
                    x = cro->x;
                 else
                    x = envobjs[ob].x;
                if (cro->y >= envobjs[ob].y && cro->y <= envobjs[ob].y2)
                    y = cro->y;
                 else
                    y = envobjs[ob].y;
                dx = (cro->x - x);
                dy = (cro->y - y);
                dist = sqrt(dx*dx+dy*dy) - cro->radius;
                angle = angv(x, y, cro->x, cro->y);
                angle = mangrelr(angle, cro->dir);
                if (angle > M_PI) angle -= PI2; // direction of the perceived object normalized between -M_PI and M_PI
                if (angle < -M_PI) angle += PI2;
                if (dist < maxdist)
                 {
                  // we compute the angle of the initial and final point of the wall
                  angle1 = angv(envobjs[ob].x, envobjs[ob].y, cro->x, cro->y);
                  angle1 = mangrelr(angle1, cro->dir);
                  if (angle1 > M_PI) angle1 -= PI2; // direction of the perceived object normalized between -M_PI and M_PI
                  if (angle1 < -M_PI) angle1 += PI2;
                  angle2 = angv(envobjs[ob].x2, envobjs[ob].y2, cro->x, cro->y);
                  angle2 = mangrelr(angle2, cro->dir);
                  if (angle2 > M_PI) angle2 -= PI2; // direction of the perceived object normalized between -M_PI and M_PI
                  if (angle2 < -M_PI) angle2 += PI2;

                  for(s=0, a1 = starta, a2 = starta + inta; s < cro->nlasersensors; s++, a1 += inta, a2 += inta)
                   {
                    // the pointer to the nearest point of wall is inside the sector
                    if (angle > a1 && angle < a2)
                     {
                      if (dist < mdist[s]) mdist[s] = dist;
                     }
                    else
                    {
                     // the first side point is inside the wall
                     if (anginside(a1, angle1, angle2))
                      {
                        // calculates the length of the vector between the sensor and the object providing that the angular offset is smaller than 60 degrees
						// for angular offset between 45 and 60 degrees the distance is increase linearly of a fraction of maxdist
						angoffset = fabs(angdelta(a1, angle));
						if (angoffset < (M_PI / 3.0))
                          {
						   ldist = dist * sin(angoffset);
						   cdist = sqrt((dist*dist)+(ldist*ldist));
						   if (angoffset > (M_PI / 4.0))
							  cdist += maxdist * ((angoffset - (M_PI / 4.0)) / (M_PI / 12.0));
						   if (cdist < mdist[s])
                             mdist[s] = cdist;
						   }
                       }
                     // the second side point is inside the wall
                     if (anginside(a2, angle1, angle2))
                      {
                        // calculates the length of the vector between the sensor and the object providing that the angular offset is smaller than 60 degrees
						// for angular offset between 45 and 60 degrees the distance is increase linearly of a fraction of maxdist
                        angoffset = fabs(angdelta(a2, angle));
                        if (angoffset < (M_PI / 3.0))
                          {
                           ldist = dist * sin(angoffset);
                           cdist = sqrt((dist*dist)+(ldist*ldist));
						   if (angoffset > (M_PI / 4.0))
						     cdist += maxdist * ((angoffset - (M_PI / 4.0)) / (M_PI / 12.0));
                           if (cdist < mdist[s])
                             mdist[s] = cdist;
                          }
                       }
                      }
					  
                     }
					 
                }
            break;
            // sampled cylinder
            case SAMPLEDSCYLINDER:
            break;
            // default
            //default:
            //	printf("ERROR: undefined object type %d \n",envobjs[ob].type);
            break;
         }
    }
    // calculate the stimulation generated by the cylindrical body of other robots
//double absa;
    if (nrobots > 1)
      {
        for (r=0, ro=rob; r < nrobots; r++, ro++)
          {
            if (cro->idn != ro->idn)
             {
                dx = (cro->x - ro->x);
                dy = (cro->y - ro->y);
                dist = sqrt(dx*dx+dy*dy) - (cro->radius * 2.0);
                angle = angv(ro->x, ro->y, cro->x, cro->y);
                angle = mangrelr(angle, cro->dir);
                if (angle > M_PI) angle -= PI2; // direction of the perceived object normalized between -M_PI and M_PI
                if (angle < -M_PI) angle += PI2;
//printf("robot %d absa %.2f robotdir %.2f angle %.2f dist %.2f\n", cro->idn, absa, cro->dir, angle, dist);
                if (dist < maxdist)
                 {
                   for(s=0, a1 = starta, a2 = starta + inta; s < cro->nlasersensors; s++, a1 += inta, a2 += inta)
                      {
                        // the center of the perceived robot is inside the sector
                        if (anginside(angle, a1, a2))
                         {
                          if (dist < mdist[s]) mdist[s] = dist;
//printf("center inside - robot %d angle %.2f sector %.2f %.2f dist %.2f\n", cro->idn, angle, a1, a2, dist);
                         }
                        else
                        {
                         if (dist < (cro->radius * 4.0)) // adjacent sectors can detect only in near proximity
                         {
                          // one of the side of the perceived robot is inside the maximum perceived angular offset of the sector
                          // when the distance is equal to the radius of the robot the angular offset is 45 degrees, longer the distance less the angularoffset
                          angoffset = angdelta(a1+(inta/2), angle);
                          maxperceivableaoffset = cro->radius / dist * (M_PI / 4);
                          if (maxperceivableaoffset > (M_PI / 4.0)) maxperceivableaoffset = M_PI / 4.0;  // cannot exceed 45 degrees
                          if (fabs(angoffset) < maxperceivableaoffset)
                            {
                             ldist = cro->radius * (angoffset / maxperceivableaoffset);
                             cdist = sqrt((dist*dist)+(ldist*ldist));
                             if (cdist < mdist[s])
                                mdist[s] = cdist;
//printf("side inside - robot %d angle %.2f sector %.2f %.2f dist %.2f angoffset %.2f maxperceivableoff %.2f ldist %.2f\n", cro->idn, angle, a1, a2, dist, angoffset, maxperceivableaoffset, ldist);
                            }
                           }
                        }


                      }
                 }
              }
           }
        }

    // normalize and copy on the imput vector
    for(s=0; s < cro->nlasersensors; s++)
        {
         *cro->csensors = (maxdist - mdist[s]) / maxdist;
         *cro->csensors = *cro->csensors * *cro->csensors * *cro->csensors;
         cro->csensors++;
        }
}


/*
 *  omnidirectional camera that perceives the color blobs constituted by the other robots
 *  robots can have their color leds turned on in red or blue
 *  it update four sensors that encode the fraction of red and blue pixels detected in the two frontal sectors of the camera
 */

int initCameraSensorRFB(struct robot *cro, int nsectors)
{

   int i;
   double  **camb;
   int *nb;

   cro->camnsectors = nsectors;

   if (cro->idn == 0)
     printf("Sensor[%d]: camera2, %d sectors, %d colors \n", cro->camnsectors * 2, cro->camnsectors, 2);

   cro->camblobsn = (int *) malloc(cro->camnsectors * sizeof(double));
   // allocate space and initialize
   cro->camblobs = (double **) malloc(cro->camnsectors * sizeof(double *));
   for (i=0, camb = cro->camblobs, nb = cro->camblobsn; i < cro->camnsectors; i++, camb++, nb++)
     {
      *camb = (double *)malloc(nrobots*4*8 * sizeof(double));
      *nb = 0;
     }
   return(cro->camnsectors * 2);
}

/*
 * add a new blob to the blob list of a sector
 * providing that it does not overlap with previously stored blobs
 * it receive as input the pointer to the first blob of the sector-list, the number of existing blobs, the sector range
 * the blob color, and angular range
 */
void updateCameraAddblob(double *cb, int *nb, double rangel, double ranger, double color, double dist, double brangel, double branger)

{

   int b;

   // we ignore small blobs with a negative intervals since they are spurious
   if ((branger - brangel) < 0.0)
   {
      return;
   }

   // if both the start and the end of the blob angle are negative we normalize them in the [0, PI2] range
   if (brangel < 0.0 && branger < 0.0)
    {
      //brangel += PI2; cb++;
      //branger += PI2; cb++;
    }
    // a blob that ends after PI2 is divided in two blobs included in the normal range
    // FOR THE MOMENT WE JUST DISCARD THE SECOND PART. THIS HAS TO BE FIXED ***************
    if (branger > PI2)
    {
     //branger = PI2;
     // second brangel - 0.0
     // second branger - branger - PI2
    }
    // a blob that start after PI2 is divided in two blobs included in the normal range
    // FOR THE MOMENT WE JUST DISCARD THE FIRST PART. THIS HAS TO BE FIXED ***************
    if (brangel < 0.0)
    {
      //brangel = 0.0;
     // second brangel = brangel -PI2
     // second branger - PI2
    }


   if (fabs(branger - brangel) > 1.6)
       printf("wrong blob %.5f %.5f %.5f %6.5f offset %.6f\n", brangel, branger, color, dist, branger - brangel);



   // check whether this blob overlap with preexisting ones
   for (b=0; b < *nb; b++)
   {
       cb++;
       cb++;
       // if fully overlap with previous blobs we do nothing
       if (anginrange(brangel, *cb, *(cb + 1)) && anginrange(branger, *cb, *(cb + 1)))
         {
          return;
         }
         else
         {
         // if start inside an existimg blob but ends after the end of the existing blobs we trim the first part of the blob
         if (anginrange(brangel, *cb, *(cb + 1)))
           {
             //printf("**** start inside (cb %.2f %.2f) (range %.2f %.2f) ", *cb, *(cb + 1), brangel, branger);
//printf("1 %.2f %.2f (%.2f %.2f)-> ", brangel, branger, *cb, *(cb + 1));
             //brangel = *(cb + 1);
             //printf(" -> (newrange %.2f %.2f) offset %.2f\n", brangel, branger, branger - brangel);

           }
           else
            {
              // if end inside an existing blob but starts outside the existing blob we trim the last part of the blob
               if (anginrange(branger, *cb, *(cb + 1)))
                   {
//printf("**** end inside   newblob (%.2f %.2f) (oldblob %.2f %.2f)-> ", brangel, branger, *cb, *(cb + 1));
                     //branger = *cb;
                     //if (branger < brangel) branger += PI2;
//printf("-> new trimmed blob %.2f %.2f - %.2f \n", brangel, branger, branger - brangel);
                   }
            }
         }
         cb++;
         cb++;
   }

   *cb = color; cb++;
   *cb = dist; cb++;
   *cb = brangel; cb++;
   *cb = branger; cb++;
   *nb += 1;


   if ((branger - brangel) < 0 || (branger - brangel) > (M_PI / 2.0))
      printf("wrong blob after processing %.5f %.5f %.5f %6.5f  %.6f\n", brangel, branger, color, dist, branger - brangel);

}


/*
 * add a new blob to the blob list of a sector
 * providing that it does not overlap with previously stored blobs
 * it receive as input the pointer to the first blob of the sector-list, the number of existing blobs, the blob color, and the start and ending angle
 * assume that the starting and ending angles are in the range [0, PI2]
 * blobs smaller than the resolution of the camera (0.1 degrees, 0.00174 radiants) are filtered out
 */
void updateCameraAddBlob(double *cb, int *nb, double color, double dist, double brangel, double branger)

{

   int b;

   // we ignore small blobs with a negative intervals since they are spurious
   if ((branger - brangel) < 0.00174)
   {
      return;
   }

   // check whether this blob overlap with preexisting ones
   for (b=0; b < *nb; b++)
   {
       cb++;
       cb++;
       // if fully overlap with previous blobs we simply filter it out
       if (anginrange(brangel, *cb, *(cb + 1)) && anginrange(branger, *cb, *(cb + 1)))
         {
          return;
         }
         else
         {
         // if start inside an existimg blob but ends after the end of the existing blobs we trim the first part of the blob
         if (anginrange(brangel, *cb, *(cb + 1)))
           {
             brangel = *(cb + 1);
           }
           else
            {
              // if end inside an existing blob but starts outside the existing blob we trim the last part of the blob
               if (anginrange(branger, *cb, *(cb + 1)))
                   {
                     branger = *cb;
                   }
            }
         }
         cb++;
         cb++;
   }

   // we ignore small blobs with a negative intervals since they are spurious
   // the blob could had become too small after being trimmed
   if ((branger - brangel) < 0.00174)
   {
      return;
   }

   *cb = color; cb++;
   *cb = dist; cb++;
   *cb = brangel; cb++;
   *cb = branger; cb++;
   *nb += 1;

}


void updateCameraSensorRFB(struct robot *cro, int *rbd)

{


	int s;						// sector
    struct robot *ro;			// pointer to robot list
	int r;						// robot id
    double v1, v2;				// visible arc of the robot (from angle v1 to angle v2)
	double x1,x2,x3,y1,y2,y3;   // the coordinate of the initial and final points of the two blobs
	double a1,a2,a3;			// angle of the initial and final points of the two adjacent blobs
	double ab1, ab2;			// the angular boundaries between the frontal and rear side
	double ab;					// the angular boundary located within the visible arc
	double d1, d2;				// distance of the two frontal/rear borders
	double ab1x, ab1y, ab2x,ab2y; // coordinates of the two borders
	int ac;						// selected front/rear border
	double ang;					// the angle of the perceiving robot from the perceived robot
	double dist2;               // the power distance between the two robots
	double rangel, ranger;		// initial and final angle of the current sector
	double cblob[3][2];			// color blobs ([Red, Blue, Green][ang-start, ang-end])
	double **camb;              // pointer to blob matrix
	double *cb;					// pointer to blobs of a sectors
	int *nb;					// pointer to the number of blobs x sectors
    double act[10];				// activation of the current visual sector (0=red, 1=blue)
	double secta;				// amplitude of visual sectors
	int c, color;				// color
    double bcolor;              // blob color
    double buf;
	int b;
    struct evonet *ne;
	std::vector<double> out;
	int i;

	
    secta = M_PI / 3.0; //PI2 / (double) cro->camnsectors; // / 3.0;        // angular amplitude of the camera sectors
    for(s=0, nb = cro->camblobsn; s < cro->camnsectors; s++, nb++)
	  *nb = 0;
	// we extract a red or blue color blob for each perceived robot
	// we stored the visible blobs divided by visual sectors and color
	// finally we compute the fraction of pixels for each sector and each color
    for (r=0; r < nrobots; r++)
	    {
          if (rbd[r] < 0)     // if the list of nearby robots ended we exit from the for
             break;
          ro=(rob + rbd[r]);  // the r nearest perceived robot
          ne=(net + rbd[r]);  // the network of the r neareast perceived robot

	out.resize(ne->noutputs);
	// Convert output activations in the proper ranges when tanh is used
	if (net->afunction == 5)
	{
		for (i = 0; i < ne->noutputs; i++)
		{
			out[i] = getOutput(ne, i);
			out[i] = 1.0 + ((out[i] - 1.0) / 2.0);
		}
	}
	else
	{
		for (i = 0; i < ne->noutputs; i++)
			out[i] = getOutput(ne, i);
	}
		

          if (1 > 0 /* cro->idn != ro->idn*/)
		   {
			// angle from perceived to perceiving robot
			ang = angv(cro->x, cro->y, ro->x, ro->y);
			// compute the visibile and coloured angular intervals
			v1 = ang - (M_PI / 2.0);
			v2 = ang + (M_PI / 2.0);
			ab1 = ro->dir - (M_PI / 2.0);
			ab2 = ro->dir + (M_PI / 2.0);
			// identify the relevant boundary (the boundary that is visible from the point of view of the perceiving robot)
			// we do that by checking the point that is nearer from the perceiving robot
			ab1x = ro->x + xvect(ab1, ro->radius);
			ab1y = ro->y + yvect(ab1, ro->radius);
			ab2x = ro->x + xvect(ab2, ro->radius);
			ab2y = ro->y + yvect(ab2, ro->radius);
			d1 =((ab1x - cro->x)*(ab1x - cro->x) + (ab1y - cro->y)*(ab1y - cro->y));
			d2 =((ab2x - cro->x)*(ab2x - cro->x) + (ab2y - cro->y)*(ab2y - cro->y));
			// the left and right border are followed and preceeded by different colors
			if (d1 <= d2)
			  {
			   ab = ab1;
			   ac = 0;
			   }
			  else
			  {
			   ab = ab2;
			   ac = 1;
			  }
			// calculate the xy coordibate of the three points located on the borders of the perceived robot
			x1 = ro->x + xvect(v2, ro->radius);
			y1 = ro->y + yvect(v2, ro->radius);
			x2 = ro->x + xvect(ab, ro->radius);
			y2 = ro->y + yvect(ab, ro->radius);
			x3 = ro->x + xvect(v1, ro->radius);
			y3 = ro->y + yvect(v1, ro->radius);
			// calculate the correspoding three angle from the point of view of the perceiving robot
			a1 = angv(x1, y1, cro->x, cro->y);
			a2 = angv(x2, y2, cro->x, cro->y);
			a3 = angv(x3, y3, cro->x, cro->y);
			// extract the angular intervals of the red and blue subsections
			if (ac == 0)
			 {
			  cblob[0][0] = a1;
			  cblob[0][1] = a1 + angdelta(a1, a2);
			  cblob[1][0] = a3 - angdelta(a2, a3);
			  cblob[1][1] = a3;
			 }
			 else
			 {
			  cblob[1][0] = a1;
			  cblob[1][1] = a1 + angdelta(a1, a2);
			  cblob[0][0] = a3 - angdelta(a2, a3);
			  cblob[0][1] = a3;
			 }

            // angles sanity checks
            for (c=0; c < 2; c++)
            {
              // if the first angle is negative the blog is over the border
              // we make both angle positive
              // it will the be divided in two blobs below because the ending angle will exceed PI2
              if (cblob[c][0] < 0)
                   {
                      cblob[c][0] += PI2;
                      cblob[c][1] += PI2;
                   }
              // if the second angle is smaller than the first and the interval is small, we invert them
              // apparently this is due to limited precision of angle calculation
              if ((cblob[c][1] - cblob[c][0]) < 0)
                 {
                    buf = cblob[c][0];
                    cblob[c][0] = cblob[c][1];
                    cblob[c][1] = buf;
                 }
             }

            /*
            for (c=0; c < 2; c++)
            {
              if ((cblob[c][1] - cblob[c][0]) < 0)
                 {
                    printf("negative %.4f %.4f   %.4f %d ", cblob[c][0], cblob[c][1], cblob[c][1] - cblob[c][0], ac);
                    if (ac == 0 && c == 0) printf("red  (%.4f %.4f %.4f) a1 a2 %.4f %.4f a1 + a1_a2 %.4f \n", a1, a2, a3, a1, a2, angdelta(a1, a2));
                    if (ac == 0 && c == 1) printf("blue (%.4f %.4f %.4f) a3 a2 %.4f %.4f a3 - a2_a3 %.4f\n", a1, a2, a3, a3, a2, angdelta(a2, a3));
                    if (ac == 1 && c == 1) printf("blue (%.4f %.4f %.4f) a1 a2 %.4f %.4f a1 + a1_a2 %.4f \n", a1, a2, a3, a1, a2, angdelta(a1, a2));
                    if (ac == 1 && c == 0) printf("red  (%.4f %.4f %.4f) a3 a2 %.4f %.4f a3 - a2_a3 %.4f \n", a1, a2, a3, a3, a2, angdelta(a2, a3));
                 }

              if ((cblob[c][1] - cblob[c][0]) > 0.8)
                 printf("large %.4f %.4f   %.4f %d\n", cblob[c][0], cblob[c][1], cblob[c][1] - cblob[c][0], ac);
            }
            */

            // we store the two blobs
            // blobs extending over PI2 are divided in two
            dist2 =((ro->x - cro->x)*(ro->x - cro->x) + (ro->y - cro->y)*(ro->y - cro->y));
            camb = cro->camblobs;
            nb = cro->camblobsn;
            cb = *camb;
            // we check whether frontal red leds are turned on or not
            if (ro->motorleds == 0 || out[2] > 0.5)
               bcolor = 1.0;
              else
               bcolor = 0.0;
            if (cblob[0][1] < PI2)
            {
              updateCameraAddBlob(cb, nb, bcolor, dist2, cblob[0][0], cblob[0][1]);
            }
            else
            {
              updateCameraAddBlob(cb, nb, bcolor, dist2, cblob[0][0], PI2);
              updateCameraAddBlob(cb, nb, bcolor, dist2, 0.0, cblob[0][1] - PI2);
            }
            // we check whether rear blue leds are turned on or not
            if (ro->motorleds == 0 || out[3] > 0.5)
               bcolor = 2.0;
              else
               bcolor = 0.0;
            if (cblob[1][1] < PI2)
            {
              updateCameraAddBlob(cb, nb, bcolor, dist2, cblob[1][0], cblob[1][1]);
            }
            else
            {
              updateCameraAddBlob(cb, nb, bcolor, dist2, cblob[1][0], PI2);
              updateCameraAddBlob(cb, nb, bcolor, dist2, 0.0, cblob[1][1] - PI2);
            }

		  }  // end if (cro->idn != ro->idn)
		}  // end for nrobots



    // sum the angular contribution of each relevant blob to each color sector
    double inrange;
    double addsrangel;  // additional sector rangel
    double addsranger;  // additional sector ranger
    int addsid;         // sector to which the additial sector belong
    int ss;             // id of the sector, usually = s, but differ for the additional sector
    double *cbadd;      // pointer to blob list used to add a new sub-blob

    // initialize to 0 neurons actiovation
    for(b=0; b < cro->camnsectors * 2; b++)
       act[b] = 0.0;

    camb = cro->camblobs;
    cb = *camb;
    nb = cro->camblobsn;
    b = 0;
    while (b < *nb)
     {
       inrange=false;
       addsid = -1;  // the id of the additional sensors is initialized to a negative number
       if (*cb == 0.0) color = -1; // black
       if (*cb == 1.0) color = 0; // red
       if (*cb == 2.0) color = 1; // blue
       //if (cro->idn == 0) printf("b %d) %.2f %.2f %.2f %.2f (%.2f) \n", b, *cb, *(cb + 1), *(cb + 2), *(cb + 3), *(cb + 3) - *(cb + 2));
       for(s=0, rangel = cro->dir - secta, ranger = rangel + secta; s < (cro->camnsectors + 1); s++, rangel += secta, ranger += secta)
         {

           if (s < cro->camnsectors)
           {
            ss = s;
            //if (cro->idn == 0) printf("sector %d (ss %d) %.2f %.2f  \n", s, ss, rangel, ranger);
            // we normalize the angles of the sector in the range [0, PI2+sectora]
            if (rangel < 0.0)
            {
             rangel += PI2;
             ranger += PI2;
            }
            // if the current sector extend over PI2 we trim it to PI2 and we initialize the additional sector
            if (ranger > PI2)
            {
              addsrangel = 0.0;
              addsranger = ranger - PI2;
              addsid=s;
            }
           }
           else
           {
            // if an additional sensor has been defined we process is otherwise we exit from the sector for
            if (addsid >= 0)
            {
              ss = addsid;
              // if (cro->idn == 1) printf(" Additional sector s %d ss %d addsid %d range %.2f %.2f\n", s, ss, addsid, addsrangel, addsranger);
              rangel = addsrangel;
              ranger = addsranger;
            }
           else
            {
             break;
            }
           }
           //if (cro->idn == 0) printf("sector %d (ss %d) %.2f %.2f  \n", s, ss, rangel, ranger);
           if (color >= 0)
           {
            if ((*(cb + 2) >= rangel) && (*(cb + 2) < ranger) && (*(cb + 3) >= rangel) && (*(cb + 3) < ranger) ) // coloured blob fully contained in the sector
             {
              act[ss * cro->camnsectors + color] += *(cb + 3) - *(cb + 2);
              //if (cro->idn == 1) printf("fullin rodir %.2f sector %d %.2f %.2f blobcolor %.2f blobang %.2f %.2f (%.2f)\n", cro->dir, s, rangel, ranger, *cb, *(cb + 2), *(cb + 3), *(cb + 3) - *(cb + 2));
              inrange=true;
             }
            else
             {
              if ((*(cb + 2) >= rangel) && (*(cb + 2) < ranger) && (*(cb + 3) >= rangel))  // non-black blob staring inside and ending outside, inside the next sector
                {
                  act[ss * cro->camnsectors + color] += ranger - *(cb + 2);
                  // we use the exceeding part to create a new blob added at the end of the blob list
                  camb = cro->camblobs;
                  cbadd = *camb;
                  cbadd = (cbadd + (*nb * 4));
                  *cbadd = *cb; cbadd++;
                  *cbadd = *(cb + 1); cbadd++;
                  *cbadd = ranger; cbadd++; // the new blob start from the end of the current sector
                  *cbadd = *(cb + 3); cbadd++;
                  *nb = *nb + 1;
                  //printf("added blob %d %.2f %.2f %.6f %.6f  range %.6f %.6f \n", *nb, *cb, *(cb + 1), ranger, *(cb + 3), rangel, ranger);
                  //if (cro->idn == 1) printf("startin rodir %.2f sector %d %.2f %.2f blobcolor %.2f blobang %.2f %.2f (%.2f)\n", cro->dir, s, rangel, ranger, *cb, *(cb + 2), ranger, ranger - *(cb + 2));
                  inrange=true;
                }

             }
           }
         }
        if (!inrange)   // blobs outsiode the view range of all sectors are turned in black for graphical purpose
           *cb = 0.0;

        cb = (cb + 4);
        b++;
     }
    //if (cro->idn == 0) printf("\n\n");


    // we finally store the value in the input units
    //printf("activation s1 red blue s2 red blue robot %d ", cro->idn);
    for(b=0; b < cro->camnsectors * 2; b++, cro->csensors++)
    {
       *cro->csensors = act[b];
       //printf("%.2f ", act[b]);
    }
    //printf("\n");

 
}


/*
 *    update the state of the bias sensor
 */

int initBiasSensor(struct robot *cro)
{

  if (cro->idn == 0) printf("Sensor[%d]: bias\n", 1);
  return(1);
}

void updateBiasSensor(struct robot *cro)
{
    
    *cro->csensors = 1.0;
    cro->csensors++;
	
}

/*
 *    update the state of the proprioceptor sensor
 */
int initProprioceptors(struct robot *cro)
{
  if (cro->idn == 0) printf("Sensor[%d]: Proprioceptor\n", 1);
  cro->proprioceptors = (double *) malloc(1 * sizeof(double));
  return(1);
}

void updateProprioceptors(struct robot *cro)
{
    
    *cro->csensors = *cro->proprioceptors;
    cro->csensors++;
    
}

/*
 *    update the state of the ground sensor on circular target areas
 *    the first neuron is turned on when the robot is on a target area with a color < 0.25
 *    the second neuron is turned on when the robot is on a target area with a color > 0.25 and < 0.75
 */
int initGroundSensor(struct robot *cro)
{
  if (cro->idn == 0) printf("Sensor[%d]: ground color \n", 2);
  return(2);
}

void updateGroundSensor(struct robot *cro)
{

	int o;
	double dx, dy, cdist;
    double act[2];
	
	act[0] = act[1] = 0.0;
	
	for (o=0; o < nenvobjs; o++)
	 {
	  if (envobjs[o].type == STARGETAREA)
	  {
	    dx = cro->x - envobjs[o].x;
		dy = cro->y - envobjs[o].y;
		cdist = sqrt((dx*dx)+(dy*dy));
		if (cdist < envobjs[o].r)
			{
			  if (envobjs[o].color[0] < 0.25)
			    act[0] = 1.0;
			   else
				if (envobjs[o].color[0] < 0.75)
			      act[1] = 1.0;
			 }
		}
	  }
    *cro->csensors = act[0];
    cro->csensors++;
    *cro->csensors = act[1];
    cro->csensors++;
}

/*
 *    update the state of the time sensor
 */
int initTimeSensor(struct robot *cro)
{

  if (cro->idn == 0) printf("Sensor[%d]: trial time \n", 1);
  return(1);
}

void updateTimeSensor(struct robot *cro, int step, int nsteps)
{
    
    *cro->csensors = (double) step / (double)nsteps;
    cro->csensors++;
    
}


/*
 *    update the state of the time sensor
 */
int initEnergySensor(struct robot *cro)
{

  if (cro->idn == 0) printf("Sensor[%d]: energy sensor \n", 1);
  return(1);
}

void updateEnergySensor(struct robot *cro)
{

    *cro->csensors = cro->energy;
    cro->csensors++;

}

/*
 * update robot position on the basis of the speed of the wheels
 * the first two output neurons are used to encode wheel speed
 
    def __init__(self,robotDiameter, axleLength, maxSpeed):
        self.robotDiameter=robotDiameter
        self.axleLength=axleLength
        self.maxSpeed=maxSpeed #mm/s
        self.halfAxleLength=axleLength/2.0
        self.robotRadius=robotDiameter/2.0
        self.x=self.oldx=0
        self.y=self.oldy=0
        self.dir=self.dir=0 #radiants
		
		self.epuck =  WheeledRobot(72,53,144) #epuck data
        self.epuck2= WheeledRobot(150,100,200) #epuck data
		*/
 
 
/*
 * update the robot position and direction
 * if the robot collides compenetrate an obstacle, moves it back of the necessary quantity to eliminate the compenetration and return 1
 * return [0] when the robot does not collide, the id of the collided robot (+1) in case of collision with a robot, -id object number -1 in case of collision with an obstacle
 */
int updateRobot(struct robot *cro, struct evonet *ne)

{
    //double axleLength = 53.0;       // distance betwen the two wheels
    double updateFrequency = 0.1;   // update frequency in s.
    //double maxSpeed = 144.0;        // max linear speed mm/s // diametro 72
    double anglVari;
    double roboPivoDist;
    double roboCentPivoDist;
    double olddir;
    double halfPI = M_PI/2.0;
    double twoPI  = M_PI*2.0;
	double dx, dy, da;
	double x,y, odx, ody;
	double leftSpeed, rightSpeed;
	double dist, smallerd;
    struct envobjects *ob;
    int nob;
    int idcoll;
	struct robot *ro;
	int r;
    double dspeed,dturn;
	int i;
	std::vector<double> out(ne->noutputs);
	
/*
if (cro->idn == 0)
{
  ne->act[ne->ninputs+ne->nhiddens] = 0.0;
  ne->act[ne->ninputs+ne->nhiddens+1] = 0.0;
  ne->act[ne->ninputs+ne->nhiddens+2] = 0.25;
}
else
{
  ne->act[ne->ninputs+ne->nhiddens] = 0.1;
  ne->act[ne->ninputs+ne->nhiddens+1] = 0.0;
  ne->act[ne->ninputs+ne->nhiddens+2] = 0.25;
}
*/

	// Convert output activations in the proper ranges when tanh is used
	if (net->afunction == 5)
	{
		for (i = 0; i < ne->noutputs; i++)
		{
			out[i] = getOutput(ne, i);
			out[i] = 1.0 + ((out[i] - 1.0) / 2.0);
		}
	}
	else
	{
		for (i = 0; i < ne->noutputs; i++)
			out[i] = getOutput(ne, i);
	}

    olddir = cro->dir;

	// calculate the moving and angular speed
    switch (cro->motorwheelstype)
	{
      case 0: // standard direct control
       leftSpeed  = ((out[0]  * 2.0 * cro->maxSpeed) - cro->maxSpeed) * updateFrequency;
       rightSpeed = ((out[1] * 2.0 * cro->maxSpeed) - cro->maxSpeed) * updateFrequency;
       break;
      case 1: // forward movement or straight backward
       dspeed = out[0] * cro->maxSpeed * updateFrequency;  // desired max speed range [0, maxspeed]
       dturn = out[1];			  // desired turn (controlateral wheel speed in the range [-dspeed, dspeed]
       if (dturn < 0.5)
         {
          leftSpeed = dspeed * (1.0 - ((0.5 - dturn) * 4.0));
          rightSpeed = dspeed;
         }
         else
         {
          leftSpeed = dspeed;
          rightSpeed = dspeed * (1.0 - ((dturn - 0.5) * 4.0));
         }
       break;
	}
	//printf("act %.2f %.2f %.2f speed %.2f %.2f\n", ne->act[ne->ninputs+ne->nhiddens],ne->act[ne->ninputs+ne->nhiddens+1], ne->act[ne->ninputs+ne->nhiddens+2], leftSpeed, rightSpeed);
    //leftSpeed = rightSpeed = (maxSpeed * 0.009);
	
    if ((fabs(leftSpeed-rightSpeed)) < 0.00001)
    {
        cro->x += dx = cos(cro->dir) * rightSpeed;
        cro->y += dy = sin(cro->dir) * rightSpeed;
		da = 0.0;
		//printf("straight) leftspeed %.2f rightspeed %.2f dir %.2f dx %.2f dy %.2f \n", leftSpeed, rightSpeed, cro->dir, dx, dy);
    }
    else
	{
        if (leftSpeed < rightSpeed)
        {
          anglVari = (rightSpeed-leftSpeed) / cro->axleLength;
          cro->dir += da = anglVari;
          roboPivoDist = leftSpeed / anglVari;
          roboCentPivoDist = roboPivoDist + (cro->axleLength/2.0);
          cro->x += dx = (cos(cro->dir - halfPI) - cos(olddir - halfPI)) * roboCentPivoDist;
          cro->y += dy = (sin(cro->dir - halfPI) - sin(olddir - halfPI)) * roboCentPivoDist;
		  //printf("l>r) angVari %.2f CentPivoDist %.2f dx %.2f dy %.2f\n", anglVari, roboCentPivoDist, cro->dx, cro->dy);
        }
        else
        {
          anglVari = (leftSpeed - rightSpeed) / cro->axleLength;
		  cro->dir += da = -anglVari;
          roboPivoDist = rightSpeed / anglVari;
          roboCentPivoDist = roboPivoDist + (cro->axleLength/2.0);
		  cro->x += dx = (cos(cro->dir + halfPI) - cos(olddir + halfPI)) * roboCentPivoDist;
		  cro->y += dy = (sin(cro->dir + halfPI) - sin(olddir + halfPI)) * roboCentPivoDist;
		  //printf("r>l) angVari %.2f CentPivoDist %.2f dx %.2f dy %.2f\n", anglVari, roboCentPivoDist, cro->dx, cro->dy);
        }
    }
	
	
	// we now check and handle possible compenetrations
	idcoll = 0;
	smallerd = 99999.0;
    // processing environmental object list
    for(nob=0, ob=envobjs; nob < nenvobjs; ob++, nob++)
    {
        switch(ob->type)
        {
            // walls
            case WALL:
                if (cro->x >= ob->x && cro->x <= ob->x2)
                    x = cro->x;
                else
                    x = ob->x;
                if (cro->y >= ob->y && cro->y <= ob->y2)
                    y = cro->y;
                else
                    y = ob->y;
                odx = (cro->x - x);
                ody = (cro->y - y);
                dist = sqrt(odx*odx+ody*ody) - cro->radius;
                if (dist <= 0.0)
                  {
                    idcoll = -nob - 1;
                    //printf("collision center %.3f %.3f dist %.3f \n", cro->x, cro->y, dist);
                    //cro->x -= dx;
                    //cro->y -= dy;
                    //vectormod = sqrt(dx*dx+dy*dy);
                    //pdist = sqrt(cro->x*cro->x+cro->y*cro->y) - cro->radius;
                    //cro->x += dx * pdist / vectormod;
                    //cro->y += dy * pdist / vectormod;
                  }
				if (dist < smallerd)
				  smallerd = dist;
             break;
            // cylinder
            case SAMPLEDSCYLINDER:
				x = ob->x;
				y = ob->y;
                odx = (cro->x - x);
                ody = (cro->y - y);
                dist = sqrt(odx*odx+ody*ody) - (cro->radius + ob->r);
                if (dist <= 0.0)
                  {
                    idcoll = -nob - 1;
                    //printf("collision center %.3f %.3f dist %.3f \n", cro->x, cro->y, dist);
                    //cro->x -= dx;
                    //cro->y -= dy;
                    //vectormod = sqrt(dx*dx+dy*dy);
                    //pdist = sqrt(cro->dx*cro->dx+dy*dy) - rob->radius;
                    //rob->x += rob->dx * pdist / vectormod;
                    //rob->y += rob->dy * pdist / vectormod;
                    //dist = sqrt(dx*dx+dy*dy) - rob->radius;
                    //printf("dist after collision %.3f \n", dist);
                  }
				if (dist < smallerd)
				  smallerd = dist;
             break;
			 
        }
    }
	// checking collisions with other robots
	if (nrobots > 1)
	  {
	    for (r=0, ro=rob; r < nrobots; r++, ro++)
	      {
		    if (cro->idn != ro->idn)
			 {
				odx = (cro->x - ro->x);
				ody = (cro->y - ro->y);
				dist = sqrt(odx*odx+ody*ody) - (cro->radius * 2.0);
				if (dist < 0.0)
				  idcoll = ro->idn + 1;
				if (dist < smallerd)
				  smallerd = dist;
			 }
		  }
		  
 	  }
	
	// in case the robot compenetrated with an obstacle we move it back
	if (idcoll != 0)
	{
      cro->dir = olddir;
	  cro->x = cro->x - dx;
	  cro->y = cro->y - dy;
	}


	// in case an obstacle at less then 15mm we check that the robot did not pass over it
    //if (smallerd < 15)
    //  {
	    //printf("near obstacle");
    //  }



	// we force robot's direction in the range [0, 2PI]
    if (cro->dir >= twoPI)
        cro->dir -= twoPI;
    if (cro->dir < 0.0)
        cro->dir += twoPI;
	
	return(idcoll);
    
}




/*
 * save a genotype in a file
 */
void savegenotype(char* filename, double* genotype, const int glen, int mode)

{

	double *gene;
	int g;
	FILE *fp;
	
	fp = fopen(filename,"w!");
	for (g=0, gene=genotype; g < glen; g++, gene++)
	  {
	   fprintf(fp, "%lf ", *gene);
	  }
	fprintf(fp, "\n");
	fclose(fp);

}

void initRobot(struct robot *cro, int n, int robottype)

{
	
    switch(robottype)
    {
       case Khepera:
          cro->radius = 27.5;
          cro->axleLength = 53.0;
          cro->maxSpeed = 144.0;
          break;
       case ePuck:
          cro->radius = 35.0;
          cro->axleLength = 55.0;
          cro->maxSpeed = 200.0;
          break;
       case MarXBot:
          cro->radius = 85.0;
          cro->axleLength = 104.0;
          cro->maxSpeed = 360.0; // (270) it was 144, then 240
          break;
		// Khepera bodyradius 0.0275 AxleLenght 0.053  wheelr 0.015
		// ePuck   bodyradius 0.035  AxleLenght 0.055  wheelr 0.0205
	    // Marxbot bodyradius 0.085  AxleLenght 0.104  wheelr 0.027

    }
	
	// general
	cro->idn = n;						// id
	cro->type = robottype;				// type
	cro->color = 0;						// robot's color 0=black, 1=red, 2=green, 3=blue, 4=red/green
	cro->alive = true;					// whether the robot is alive
	
    cro->x = 100;						// x-pos
	cro->y = 100;						// y-pos
    cro->dir = 0.0;                     // direction
    //*cro->sensors = NULL;               // vector of input states
    //*cro->csensors = NULL;              // pointer to the first sensor to be updated
    cro->dx = 0.0;                      // last deltax
    cro->dy = 0.0;                      // last deltay
    cro->energy = 0.0;                  // the energy level of the robot
	
	// sensors
	cro->sensorinfraredid = 0;		    // id of the first infrared sensor neuron
	cro->sensorinfraredn = 0;           // n of infrared sensory neurons
	cro->sensorcameraid = 0;		    // id of the first camera sensory neuron
	cro->sensorcameran = 0;             // n of infrared sensory neurons
    //*cro->proprioceptors = NULL;		// robots' proprioceptors state
	cro->camnsectors = 0;				// the number of camera setctors
	//**cro->camblobs = NULL;				// matrix of color blobs delected by the camera, firstlevel (sector) secondlevel (blobs)
	//*cro->camblobsn = NULL;             // the number of color blobs for each sector
	
	// motors
    cro->motorwheels = 0;               // number of motors neurons controlling the wheels
    cro->motorwheelstype = 0;           // type of motor wheel control (0=direct, 1=traslation-rotation)
    cro->motorleds = 0;                 // number of motors neurons controlling the wheels
	cro->motorwheelsid = 0;             // id of the first motorwheels neuron
    cro->motorledsid = 0;               // id of the first motorleds neuron
}

void initRobotSensors(struct robot *cro, int nsensors)

{

    int s;
	
    cro->sensors = (double *) malloc(nsensors * sizeof(double));
    for(s=0, cro->csensors = cro->sensors; s < nsensors; s++, cro->csensors++)
       *cro->csensors = 0.0;
	cro->csensors = cro->sensors;
}