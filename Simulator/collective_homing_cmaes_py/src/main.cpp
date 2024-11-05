/*
 * Evorobot* - A software for running evolutionary robotics experiments
 * Copyright (C) 2005
 * LARAL, Institute of Cognitive Science and Technologies, CNR, Roma, Italy
 * http://laral.istc.cnr.it

 * This program is free software; you
 ls can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


#include "main.h"
#include <QTime>
#include <QtGui>

#include <math.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utilities.h"
#include "robot-env.h"
#include "evonet.h"
#include "evocmaes.h"

extern double evaluate(double* genotype, const int glen, int mode, int seed);
extern int initialize(const char* filename);

//MainWindow *mainWin;				// The main window
WorldWidget *worldWidget;			// The world graphic widget
ActWidget *actWidget;               // The activation widget

// wait for n milliseconds
void timedelay(int milliseconds)
{
    QTime dieTime = QTime::currentTime().addMSecs( milliseconds );
    while( QTime::currentTime() < dieTime )
    {
        QCoreApplication::processEvents( QEventLoop::AllEvents, 100 );
    }
}

// update the 2D world renderer
void update_world()

{

    worldWidget->update();
    actWidget->update();
    timedelay(steptimelength);
    QCoreApplication::processEvents();

}

// Activation Widget constructor
ActWidget::ActWidget(QWidget *parent)
    : QWidget(parent)
{

    double dx,dy;

    dy = nrobots * 12.0;
    dx = (net->ninputs + net->nhiddens + net->noutputs) * 12.0;

    setBackgroundRole(QPalette::Base);
    resize(dx, dy);
    setMinimumSize(dx, dy);

    setWindowTitle("Activation");

    QPalette pal = palette();
    pal.setColor(QPalette::Background, Qt::white);
    setAutoFillBackground(true);
    setPalette(pal);
}

// Activation Widget painter
void ActWidget::paintEvent(QPaintEvent *)
{

    double sx, sy;
    int r,n;
    struct robot *ro;
    struct evonet *ne;
    double act;

    QPainter painter;
    painter.begin(this);

    // draw neurons activation for all robots
    for (r=0, ro = rob, ne = net, sy=0.0; r < nrobots; r++, ro++, ne++, sy+= 12.0)
     {
       for (n=0, sx = 0.0; n < (ne->ninputs + ne->nhiddens + ne->noutputs); n++, sx += 12.0)
         {
           painter.setPen(Qt::black);
           if (n < ne->ninputs)
             painter.setBrush(QBrush(Qt::red, Qt::SolidPattern));
            else
             if (n < (ne->ninputs + ne->nhiddens))
                painter.setBrush(QBrush(Qt::darkGray, Qt::SolidPattern));
               else
                painter.setBrush(QBrush(Qt::blue, Qt::SolidPattern));
           painter.drawRect(sx, sy, 10.0, 10.0);
           act = (1.0 - ne->act[n]) * 10.0;
           painter.setPen(Qt::white);
           painter.setBrush(QBrush(Qt::white, Qt::SolidPattern));
           painter.drawRect(sx, sy, 10.0, act);
         }
     }
    painter.end();
}


// 2d World widget constructor
WorldWidget::WorldWidget(QWidget *parent)
    : QWidget(parent)
{

    double maxd;
    double maxrate;

    setBackgroundRole(QPalette::Base);
    setMinimumSize(100, 100);

    setBackgroundRole(QPalette::Base);
    move(0,0);

    // we scaledown the size if the x or y axes exceeds 1000
    if (worldx > worldy)
        maxd = worldx;
       else
        maxd = worldy;
    if (maxd > 1000)
       maxrate = maxd / 1000.0;
      else
       maxrate = 1.0;
    resize(worldx/maxrate, worldy/maxrate);

    setWindowTitle("Robots/Environment");

    QPalette pal = palette();
    pal.setColor(QPalette::Background, Qt::white);
    setAutoFillBackground(true);
    setPalette(pal);
}


// 2D world widget painter
void WorldWidget::paintEvent(QPaintEvent *)
{

	int o,r;
	struct robot *ro;
	struct evonet *ne;
	double an, ain,ann;
	double **camb;
	double *cb;
	int *nb;
	double x,y,x2,y2;
    int cx, cy;
    double cxx, cyy;
	int i;
	int b;
	int s;
    double color;
	QColor col;
    double roradius2;
    double angpie2;
 
 
 
    QPainter painter;
    painter.begin(this);
    painter.scale(width()/worldx,height()/worldy);
    painter.setRenderHint(QPainter::Antialiasing);
	QPen pen(Qt::black, 2);
	painter.setPen(pen);

	// paint visited virtual cells
    if (cellsize > 0.0)
     {
       painter.setPen(QColor(0, 0, 255, 20));
       painter.setBrush(QBrush(QColor(0, 0, 255, 20), Qt::SolidPattern));
       for (cx=0, cxx=0.0; cxx < worldx; cx++, cxx += cellsize)
          for (cy=0, cyy=0.0; cyy < worldy; cy++, cyy += cellsize)
             if (cells[cx][cy] != 0)
				{
                painter.drawRect(cxx, cyy, cellsize, cellsize);
				}
     }
	// paint environmental objects
	for (o=0; o < nenvobjs; o++)
	 {
	   if (envobjs[o].type == WALL)
	     painter.drawLine(envobjs[o].x, envobjs[o].y, envobjs[o].x2, envobjs[o].y2);
	   if (envobjs[o].type == SAMPLEDSCYLINDER)
		 painter.drawEllipse(envobjs[o].x - envobjs[o].r, envobjs[o].y - envobjs[o].r, envobjs[o].r*2.0, envobjs[o].r*2.0);
	   if (envobjs[o].type == STARGETAREA)
	    {
		 if (envobjs[o].color[0] == 0.0)
		   painter.setBrush(QBrush(QColor(50, 50, 50, 255), Qt::SolidPattern));
		 if (envobjs[o].color[0] == 0.5)
		   painter.setBrush(QBrush(QColor(150, 150, 150, 255), Qt::SolidPattern));
      //printf("o %d x %lf y %lf r %lf\n", o, envobjs[o].x, envobjs[o].y, envobjs[o].r);
		 painter.drawEllipse(envobjs[o].x - envobjs[o].r, envobjs[o].y - envobjs[o].r, envobjs[o].r*2.0, envobjs[o].r*2.0);
		}
	 }
	
	
	for (r=0, ro = rob; r < nrobots; r++, ro++)
	 {
	 
	   // we TEMPORARILY ASSUME THAT THE NETWORK OF THE ROBOT CORRESPOND TO THE ROBOT ID
	   ne = (net + ro->idn);
	
        // display the color of the robot
        if ((ro->color >= 1) && (ro->color <= 3))
        {
          switch (ro->color)
           {
             case 1:
              painter.setPen(Qt::red);
              painter.setBrush(QBrush(Qt::red, Qt::SolidPattern));
              break;
             case 2:
              painter.setPen(Qt::blue);
              painter.setBrush(QBrush(Qt::blue, Qt::SolidPattern));
              break;
             case 3:
              painter.setPen(Qt::green);
              painter.setBrush(QBrush(Qt::green, Qt::SolidPattern));
              break;
            }
          painter.drawEllipse(ro->x - (ro->radius * 0.75), ro->y - (ro->radius * 0.75), ro->radius*1.5, ro->radius*1.5);
        }
	

        if (ro->color == 4)
        {
          //printf("id %d nneurons %d act %lf\n", ro->motorledsid, ne->nneurons, ne->act[ro->motorledsid]);
          roradius2 = ro->radius * 0.75;
          if (ro->motorleds == 0 || ne->act[ro->motorledsid] > 0.5)
             {
               //printf("COLOR RED\n");
               angpie2 = (360 - (ro->dir / PI2 * 360.0));
               angpie2 -= 90.0; if (angpie2 < 0.0) angpie2 += 360.0;
               painter.setPen(Qt::red);
               painter.setBrush(QBrush(Qt::red, Qt::SolidPattern));
               painter.drawPie(ro->x - roradius2, ro->y - roradius2, roradius2 * 2.0, roradius2 * 2.0, angpie2*16, (180*16));
             }
          if (ro->motorleds == 0 || (ro->motorleds > 1 && ne->act[ro->motorledsid+1] > 0.5))
             {
               //printf("COLOR BLUE\n");
               angpie2 = (360 - (ro->dir / PI2 * 360.0));
               angpie2 += 90.0; if (angpie2 > 360.0) angpie2 -= 360.0;
               painter.setPen(Qt::blue);
               painter.setBrush(QBrush(Qt::blue, Qt::SolidPattern));
               painter.drawPie(ro->x - roradius2, ro->y - roradius2, roradius2 * 2.0, roradius2 * 2.0, angpie2*16, (180*16));
             }
        }
		 

       // paint the boudary circle and the frontal line in black
       painter.setBrush(QBrush(Qt::black, Qt::NoBrush));
       painter.setPen(QColor(0, 0, 0, 255));
	   painter.drawEllipse(ro->x - ro->radius, ro->y - ro->radius, ro->radius*2.0, ro->radius*2.0);
	   painter.drawLine(ro->x, ro->y,ro->x + xvect(ro->dir, ro->radius), ro->y + yvect(ro->dir, ro->radius));

       // paint the infrared sensors
	   // we TEMPORARILY assume they are encoded in the first 8 sensors)
	   if (ro->type != Khepera)
	    {
	     for(s=0, an = ro->dir - M_PI, ain = (M_PI / 4.0); s < ro->nifsensors; s++, an += ain)
	      {
		    x = ro->x + xvect(an, ro->radius);
		    y = ro->y + yvect(an, ro->radius);
		    ann = an + M_PI;
            x2 = x + xvect(ann, ro->radius * ne->act[s]);
            y2 = y + yvect(ann, ro->radius * ne->act[s]);
            painter.setPen(Qt::green);
		    painter.drawLine(x, y, x2, y2);
          }
		}
		else
		{
	     for(s=0, an = ro->dir - (M_PI / 2.0), ain = (M_PI / 5.0); s < ro->nifsensors; s++, an += ain)
	      {
			if (s == 6) an = ro->dir + (M_PI - (M_PI / 6.0));
			if (s == 7) an = ro->dir + (M_PI + (M_PI / 6.0));
		    x = ro->x + xvect(an, ro->radius);
		    y = ro->y + yvect(an, ro->radius);
		    ann = an + M_PI;
            x2 = x + xvect(ann, ro->radius * ne->act[s]);
            y2 = y + yvect(ann, ro->radius * ne->act[s]);
            painter.setPen(Qt::green);
		    painter.drawLine(x, y, x2, y2);
          }
		}
		
	   // paint the laser sensors
	   /* we TEMPORARILY assume they are encoded in the first 8 sensors)
	   for(s=0, an = ro->dir - M_PI, ain = (M_PI / 4.0); s < ro->nlasersensors; s++, an += ain)
	    {
		  x = ro->x + xvect(an, ro->radius);
		  y = ro->y + yvect(an, ro->radius);
		  ann = an + M_PI;
          x2 = x + xvect(ann, ro->radius * ne->act[s]);
          y2 = y + yvect(ann, ro->radius * ne->act[s]);
		  //printf("infrared ang %.2f %.2f act %.3f\n", an, ann, ne->act[s]);
          painter.setPen(Qt::green);
		  painter.drawLine(x, y, x2, y2);
        }
		*/
	

       // paint the angular and distance sensor (predprey experiment) from sensor 8 on
       //painter.setPen(Qt::blue);
       //for (s=0, an = ro->dir + (PI2 / 16.0); s < 8; s++, an += (PI2 / 8.0))
       //  {
       //    x = ro->x + xvect(an, ne->act[8+s] * ro->radius);
       //    y = ro->y + yvect(an, ne->act[8+s] * ro->radius);
       //    painter.drawLine(ro->x, ro->y, x, y);
       //  }

       // paint the color blobs detected by the camera
       if (ro->camnsectors > 0)
	    {
		for(s=0, camb = ro->camblobs, nb = ro->camblobsn; s < ro->camnsectors; s++, camb++, nb++)
	    {
	     for(b=0, cb = *camb; b < *nb; b++)
	      {
//if (ro->idn == 0) printf("*SECTOR %d number blobs %d BLOB %d %.2f %.2f %.2f %.2f\n", s, *nb, b, *cb, *(cb + 1), *(cb + 2), *(cb + 3));

          color = *cb;
          if (color == 0.0)
            painter.setPen(Qt::white);
          if (color == 1.0)
	        painter.setPen(Qt::red);
          if (color == 2.0)
	        painter.setPen(Qt::blue);
		  if (color == 3.0)
	        painter.setPen(Qt::green);
	      cb++; cb++;
          if (color > 0.0)
           for(i=0, an = *cb; (i < 180 && an < *(cb + 1)); i++, an += 0.03644)
		    {
	         painter.drawLine(ro->x, ro->y,ro->x + xvect(an, ro->radius), ro->y + yvect(an, ro->radius));
			}
//if (i > 40) printf("i %d sector %d blob %d %.2f %.2f ***********************\n", i, s, b, *cb, *(cb + 1));
	      cb = (cb + 2);
	     }
	    }
	   }
		 


	   // paint the color blobs detected by the camera
       //if (ro->camnsectors > 0)
       // {
       //    double secta = M_PI / 2.0;
       //    double sectx, secty;
       //    double act;
       //    ne = (net + ro->idn);
			
       //    for(s=0, an = ro->dir - secta + (secta / 2.0); s < ro->camnsectors; s++, an += secta)
       //     {
       //        sectx = ro->x + xvect(an, ro->radius / 3.0 * 2.0);
       //        secty = ro->y + yvect(an, ro->radius / 3.0 * 2.0);
//act = 1.0;
//printf("%.2f \n", ne->act[8+s*2+1]);
               //act = ne->act[1] * 15.0;
//               painter.setPen(QColor(255, 0, 0, 175));
//               painter.drawEllipse(sectx - act, secty - act, (act * 2.0), (act * 2.0));
               //act = ne->act[8+s*2+1] * 15.0;
//               painter.setPen(QColor(0, 0, 255, 175));
//               painter.drawEllipse(sectx - act, secty - act, (act * 2.0), (act * 2.0));
//             }
			

        //}
	
		 
		 
		 
     }
	 painter.end();
	
	 
		 

}


/* create the main window
MainWindow::MainWindow()
{


    Window* window = new Window(this);
    //ActWindow* actWindow = new ActWindow(this);

    //setMinimumSize(worldx+20, worldy+20);
    QScrollArea* scroll = new QScrollArea( this );
    setCentralWidget( scroll );
    scroll->setWidget( window );

    setWindowTitle("2D World");

}


Window::Window(MainWindow* main)
{

    worldWidget = new WorldWidget;
    worldWidget->setMinimumSize(worldx, worldy);

	QPalette pal = palette();
	pal.setColor(QPalette::Background, Qt::white);
    worldWidget->setAutoFillBackground(true);
    worldWidget->setPalette(pal);

    QGridLayout *mainLayout = new QGridLayout;
    mainLayout->addWidget(worldWidget, 0, 0, 1, 2 );
    setLayout(mainLayout);

}
*/

/*
 * copy the content of the ini file into a file SN.log
 */
void copyIniFile(char* filename, int seed)
{

    char fname[1024];
    char buff[1024];
    FILE* fp1;
    FILE* fp2;

    sprintf(fname, "S%d.log", seed);

    fp1 = fopen(filename, "r");
    fp2 = fopen(fname, "w");
    if (fp1 != NULL)
    {
        // Read lines
        while (fgets(buff, 1024, fp1) != NULL)
        {
          fprintf(fp2, "%s", buff);
        }
        fclose (fp1);
        fclose (fp2);
    }
    else
    {
        printf("ERROR: unable to open file %s\n", filename);
        fflush(stdout);
    }
}


/*
 * describe parameters
 */
void usage()

{
   printf("Evorobot 3.0\n");
   printf(".ini parameter is interpreted as the name of the configuration file\n");
   printf(".bgen or .gen parameter interpreted as the name of a file containing a genotype to be tested\n");
   printf("<number> is interpreted as the seed \n");
   printf("test trigger the execution of the test() function\n");
}


extern void test(int ngenes);

/*
 * Main function
 */
int main(int argc, char *argv[])
{

	 int ngenes;
     char fileini[128];
     char filebgen[128];
     char filegen[128];
     int  pseed;
     bool ufileini;
     bool ufilebgen;
     bool ufilegen;
     bool useed;
     bool utest;
     bool ueval;
	bool autotest;
	bool ureplicas;
	int replicas;
	 double f;
     int a;
     bool setGlobDir = false;
     bool setFileDir = false;
	bool setSeed = false;
     char filename[512];

     // initialize the randomgenerator
     RNG = new RandomGenerator();
    
     // graphic on or off
	 renderrobotworld = 0;

	 // set USA local conventions
	 setlocale( LC_ALL, "en-US" );
	 setlocale(LC_NUMERIC, "C");

     // parse arguments
     ufileini = false;
     ufilebgen = false;
     ufilegen = false;
     useed = false;
     utest = false;
     ueval = false;
	autotest = false;
	ureplicas = false;
	replicas = 1;
     bestggenotype = NULL;
     for(a=1; a < argc; a++)
     {
        if (filenendwith(argv[a], ".ini"))
        {
           sprintf(fileini, "%s", argv[a]);
           ufileini = true;
        }
        else if (filenendwith(argv[a], ".bgen"))
        {
           sprintf(filebgen, "%s", argv[a]);
           ufilebgen = true;
        }
        else if (argv[a][0] >= 48 && argv[a][0] <= 57)
        {
		if (!setSeed)
		{
          		pseed = (int)strtol(argv[a], NULL, 10);
          		useed = true;
			setSeed = true;
		}
		else
		{
			replicas = (int)strtol(argv[a], NULL, 10);
			ureplicas = true;
		}
        }
	else if (strchr(argv[a], '/') != NULL)
	{
	    if (!setGlobDir)
	    {
	    	sprintf(globdir, "%s", argv[a]);
	    	if (!filenendwith(globdir, "/"))
		    strcat(globdir, "/");
		setGlobDir = true;
	    }
	    else
	    {
		sprintf(filedir, "%s", argv[a]);
	    	if (!filenendwith(filedir, "/"))
		    strcat(filedir, "/");
		setFileDir = true;
	    }
	}
        else if (strcmp(argv[a], "test")==0)
        {
          utest = true;
        }
        else if (strcmp(argv[a], "eval")==0)
        {
          ueval = true;
        }
        else if (filenendwith(argv[a], ".gen"))
        {
           sprintf(filegen, "%s", argv[a]);
           ufilegen = true;
        }
	else if (strcmp(argv[a], "autotest")==0)
		autotest = true;
     }
	
	 // display usage information (the first parameter is mandatory)
     if (!ufileini)
     {
         printf("ERROR: you need to specify a .ini file as first parameter\n");
         usage();
     }

	if (!setGlobDir)
        {
            printf("You must specify the directory\n");
            exit(-1);
        }

	strcpy(filename, globdir);
	strcat(filename, fileini);

	 // initialize
     ngenes = initialize(filename/*fileini*/);
	 if (ngenes <= 0)
	     printf("ERROR: Failure during the initialization process\n");

	 // read EVO parameters
     if (ufileini)
       readEvoParameters(filename/*fileini*/);

	// Set the number of internal neurons to algorithm
	setInternalNeurons(net->nhiddens + net->noutputs);

     if (useed)
     {
       seed = pseed;
     }

     RNG->setSeed(seed);

	if (ureplicas)
	{
		nreplications = replicas;
	}

     initES();

     if (ufilebgen)
     {
       bestggenotype = (double *)malloc(ngenes * sizeof(double));
       loadGenotype(argv[2], bestggenotype, ngenes);
     }

     if (ufilegen)
     {
         bestggenotype = (double *)malloc(ngenes * sizeof(double));
         double *gen;
         int g;
         FILE *fp;
         fp = fopen(argv[2],"r");
         for(g=0, gen=bestggenotype; g < ngenes; g++, gen++)
             fscanf(fp, "%lf ", gen);
         fclose(fp);
     }

     if (!utest && !ueval && !ufilebgen && !ufilegen && !autotest)
	 {
		printf("nreplications = %d - seed %d\n", nreplications, seed);
	  for(int r=0; r < nreplications; r++, seed++)
	   {
	    RNG->setSeed(seed);
        evolve(ngenes, seed);
		copyIniFile(filename/*fileini*/, seed);
	   }
	 }

	if (autotest)
	{
		printf("nreplications = %d - seed %d\n", nreplications, seed);
		// Automatic test: script must be run from the directory containing .bgen files!!!!
		createInitialStates(0);
		// Allocate space for genotype once only
		if (bestggenotype == NULL)
		{
			bestggenotype = (double *)malloc(ngenes * sizeof(double));
		}
		// Test all best genotypes automatically!!!
		int sseed = seed;
		int eseed = seed + nreplications;
		int cseed;
		for (cseed = sseed; cseed < eseed; cseed++)
		{
			// Build filename bestgS*.bgen
			char genfile[256];
			char seedfile[64];
			sprintf(seedfile, "%d", cseed);
			strcpy(genfile, "bestgS");
			strcat(genfile, seedfile);
			strcat(genfile, ".bgen");
			// Load file
			loadGenotype(genfile, bestggenotype, ngenes);
			printf("Testing %s with seed %d steplength %d \n", genfile, seed, steptimelength);
			f = evaluate(bestggenotype, ngenes, 1, cseed);
			printf("best g genotype retest %.2f\n", f);
		}
	}
	else
	{

     if (utest)
       test(ngenes);

      if (ufilebgen || ufilegen || ueval)
      {
       createInitialStates(0);
       renderrobotworld = 1;

       QApplication app(argc, argv);

       worldWidget = new WorldWidget;
       worldWidget->show();

       actWidget = new ActWidget;
       actWidget->move(worldWidget->width() + 10, 0);
       actWidget->show();

       printf("Testing %s with seed %d steplength %d \n", argv[2], seed, steptimelength);
       f = evaluate(bestggenotype, ngenes, 1, seed);
       printf("best g genotype retest %.2f\n", f);
		 
       return app.exec();
      }
	}
  finalizeES();
	
}
