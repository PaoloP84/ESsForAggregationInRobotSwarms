/********************************************************************************
 *  FARSA Genetic Algorithm Library                                             *
 *  Copyright (C) 2007-2008 Gianluca Massera <emmegian@yahoo.it>                *
 *                                                                              *
 *  This program is free software; you can redistribute it and/or modify        *
 *  it under the terms of the GNU General Public License as published by        *
 *  the Free Software Foundation; either version 2 of the License, or           *
 *  (at your option) any later version.                                         *
 *                                                                              *
 *  This program is distributed in the hope that it will be useful,             *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
 *  GNU General Public License for more details.                                *
 *                                                                              *
 *  You should have received a copy of the GNU General Public License           *
 *  along with this program; if not, write to the Free Software                 *
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
 ********************************************************************************/

// ang constants
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef PI2
#define PI2 (M_PI * 2.0)
#endif

class RandomGenerator;
class RandomGeneratorPrivate;

/*! Global Random Generator 
 *  \warning this may be not-thread safe
 * \ingroup utilities_rng
 */
extern RandomGenerator* RNG;

/*!  \brief Random Generator Utility Class
 *
 *  \par Description
 *    Create, Manage Random Numbers Generators
 *  \par Warnings
 *
 * \ingroup utilities_rng
 */
class RandomGenerator {
public:
	/*! Default Constructor */
	RandomGenerator( int seed = 0 );
	/*! Destructor */
	~RandomGenerator();
	/*! set the seed of this random number generator */
	void setSeed( int seed );
	/*! Return the seed setted */
	int seed();
	/*! return true with probability trueProbability, otherwise false */
	bool getBool( double trueProbability );
	/*! return a random number within range specified (extreme inclusive) with a Flat distribution */
	int getInt( int min, int max );
	/*! return a random number within range specified (extreme inclusive) with a Flat distribution */
	double getDouble( double min, double max );
	/*! return a random number accordlying to a Gaussian distribution
	 *  \param var is the variance of the Gaussian distribution
	 *  \param mean is the centre of the Gaussian distribution
	 */
	double getGaussian( double var, double mean = 0.0 );

private:
	/*! encapsulate all third-party dependent code */
	RandomGeneratorPrivate* prive;
	/*! Seed */
	int seedv;
};

void sortDoubles(int n, double *vector, int *rank);

double xvect(double ang, double module);
double yvect(double ang, double module);
double angv(double x1, double y1, double x2, double y2);
double angdelta(double a1, double a2);
bool anginrange(double a, double r1, double r2);
double mangrelr(double absolute, double orient);

void copyandclear(char *s, char *sc);
bool filenendwith(char *str, char *suffix);
bool parsebool(char *str);
