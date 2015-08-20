/***************************
*  File masplibrary.c      *
*  Lost of useful stuff    *
***************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/********************************************************
* Local memory allocation. Note that for each core I    *
* initialise a large block using the standard malloc    *
* From then on, I manage my own space. Memory is always *
* de-allocated in the opposite order to its allocation  *
* which makes matters very simple. soFar is a memory    *
* stack pointer.                                        *
********************************************************/

int soFar;

char * mem;
void initialiseMemory()
{
	mem = (char *) (malloc(3500000000*sizeof(char)));
	soFar = 0;
}

void * localMalloc(int p)
{
	char * q =mem+soFar;
	soFar +=p;
	return (soFar  >= MEMSIZE)?NULL:(void *)(q);
}



/****************************
* Various utilities follow  *
****************************/
/*********************************
* This function takes the coded  *
* address of a cell and inverts  *
* the order of the bits. It is   *
* reversible. The point is to    *
* make tree building a lot more  *
* efficient (up to 300%)         *
*********************************/

long invert(long p)
{
	long q = 0;
	int j;
	for(j=0; j<64; j++)
	{
		long w = p&1;
		p >>=1;
		q<<=1;
		q+=w;
    }
     return q;
}


/************************************
* This two functions convert a set  *
* of six-dimensional coordinates    *
* into a single unambiguous         *
* unsigned long (and vice versa)    *
************************************/

unsigned long map (short g[6])
{
	unsigned long p = 0;
	int j;
	for (j=0; j<6; j++)
	{
		p *= (dims[j].resolution+1);
        p+= g[j];
    }
    p= invert(p);
    return p;
}

void expand(unsigned long p, short g[6])
{
	int j;
    p=invert(p);
    for (j=5; j>= 0; j--)
    {
    	g[j] = p % (dims[j].resolution+1);
    	p /= (dims[j].resolution+1);
    }
}


/*****************************
* Tests if a character is an *
* upper-case letter          *
*****************************/

int isUpper(char w)
{
	return (w >='A' && w <= 'Z');
}

/*****************************
* Tests if a character is an *
* lower-case letter          *
*****************************/

int isLower (char w)
{
	return (w >='a' && w <= 'z');
}

/*****************************
* Converts lower-case letter *
* to upper-case              *
*****************************/

int toUpper (char w)
{
	if (isLower(w))
	  return (w +'A' - 'a');
	else
	  return w;
}

/*****************************
* Tests if string w is empty *
*****************************/

int empty(char * w)
{
	unsigned int j;
	for ( j = 0; j < strlen(w)-1; j++)
    if (w[j] > ' ')
	  return 0;
	return 1;
}


/***********************************************
* Random number generator taken from Numerical *
* recipes in C (page 210). Said to be very     *
* good indeed!                                 *
* The parameter idum is an integer.If negative *
* it re-initialises the sequence of numbers.   *
* Otherwise it should be 0 for the next number *
*                                              *
***********************************************/





double grandom(int* idum)
{
	/* Returns a uniform random deviate between 0 and 1. Set idum to any negative
	 value to initialse the sequence */
	static long ix1,ix2,ix3;
	static double r[98];
	double temp;
	static int iff=0;
	int j;
	if (*idum < 0 || iff==0)
	{
		iff = 1;
		ix1 = (IC1-(*idum))%M1;
		ix1 = (IA1*ix1+IC1)%M1;
		ix2 = ix2%M2;
		ix1 = (IA1*ix1+IC1)%M1;
		ix3 = ix1 %M3;
		for(j=1; j<=97; j++)
		{
			ix1 = (IA1*ix1+IC1)%M1;
			ix2 = (IA2*ix2+IC2)%M2;
			r[j]=(ix1+ix2*RM2)*RM1;
		}
		*idum = 1;
	}
	ix1 = (IA1*ix1+IC1)%M1;
	ix2 = (IA2*ix2+IC2)%M2;
	ix3 = (IA3*ix3+IC3)%M3;
	j= 1+((97*ix3)/M3);
	if (j>97 || j < 1)
		printf("impossible happened: j = %d\n",j);
	temp = r[j];
	r[j]=(ix1+ix2*RM2)*RM1;
	return temp;
}

//  End of sooper-dooper random number generator


/**************************************
*   This is a crummy Gauss function   *
*   used only for hign values (>15)   *
*   of the Poisson distribution       *
**************************************/

double gauss(double mu, double sigma)
{
	double t = -6;
	int j;
	for (j=0; j<12; j++)
	   t+= grandom(&idum);
	return sigma*t+mu;
}

/************************************************
* Poisson distribution, The parameter is a mean *
* that does not have to be an integer. The      *
* return is an integer from the distribution    *
* that corresponds to the mean.                 *
* This function might have been in the C        *
* library I had - but it wasn't.                *
************************************************/

int poisson (double lambda)
{
	double elambda = exp(-lambda);
	double target = grandom(&idum);
	double term = elambda;
	double sum = term;
	int n = 0;
	if (lambda > 15)
	{    // Use gaussian distribution instead

		return (int)(gauss(lambda, sqrt(lambda)));
	}
	while ( sum < target)
	{
		term *= (lambda/(++n));
		sum+= term;
    }
    return n;
}


/*****************************************************
*   Numerical approximation to erf(x)  (Wikipedia)   *
*   Again, this function might have been in the C    *
*   library I had - but it wasn't.                   *
*****************************************************/


double erf(double x)
{
	 int w;
	 double tau;
	 double z = eparams[9];
	 double t = 1.0/(1+0.5*fabs(x));
	 for (w = 8; w>= 0; w--)
	    z = z*t+eparams[w];
     tau = t * exp(z-x*x);
	 if (x>=0)
	   return (1-tau);
	 else
	   return (tau-1);
}

