/******************************
*  File masplibrary.c         *
*  Lost of useful stuff       *
*  Version 0.7  (21/08/2016)  *
******************************/
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
{ //1
	mem = (char *) (malloc(3500000000*sizeof(char)));
	soFar = 0;
} //1

void * localMalloc(int p)
{ //1
	char * q =mem+soFar;
	soFar +=p;
	return (soFar  >= MEMSIZE)?NULL:(void *)(q);
} //1



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
{ //1
	long q = 0;
	int j;
	for(j=0; j<64; j++)
	{ //2
		long w = p&1;
		p >>=1;
		q<<=1;
		q+=w;
    } //2
     return q;
} //1


/************************************
* This two functions convert a set  *
* of six-dimensional coordinates    *
* into a single unambiguous         *
* unsigned long (and vice versa)    *
* Modified 151080 to ensure that    *
* negative coordinates and co-      *
* ordinates slightly over the res-  *
* olition are correctly handled     *
************************************/

unsigned long map (short g[6])
{ //1
	unsigned long p = 0;
	int j;
	for (j=0; j<6; j++)
	{ //2
		p *= (dims[j].resolution+100);
        p+= (g[j]+16);
    } //2
    p= invert(p);
    return p;
} //1

void expand(unsigned long p, short g[6])
{ //1
	int j;
    p=invert(p);
    for (j=5; j>= 0; j--)
    { //2
    	g[j] = p % (dims[j].resolution+100) -16;
    	p /= (dims[j].resolution+100);
    } //2
} //1


/*****************************
* Tests if a character is an *
* upper-case letter          *
*****************************/

int isUpper(char w)
{ //1
	return (w >='A' && w <= 'Z');
} //1

/*****************************
* Tests if a character is an *
* lower-case letter          *
*****************************/

int isLower (char w)
{ //1
	return (w >='a' && w <= 'z');
} //1

/*****************************
* Converts lower-case letter *
* to upper-case              *
*****************************/

int toUpper (char w)
{ //1
	if (isLower(w))
	  return (w +'A' - 'a');
	else
	  return w;
} //1

/*****************************
* Tests if string w is empty *
*****************************/

int empty(char * w)
{ //1
	unsigned int j;
	for ( j = 0; j < strlen(w)-1; j++)
    if (w[j] > ' ')
	  return 0;
	return 1;
} //1


/***********************************************
* Random number generator taken from Numerical *
* recipes in C (page 210). Said to be very     *
* good indeed!                                 *
* The parameter idum is an integer.If negative *
* it re-initialises the sequence of numbers.   *
* Otherwise it should be 0 for the next number *
*                                              *
***********************************************/


/***********************************************
* Constants for the random number generator    *
***********************************************/

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349
#define PI 3.141592654



double grandom(int* idum)
{ //1
	/* Returns a uniform random deviate between 0 and 1. Set idum to any negative
	 value to initialse the sequence */
	static long ix1,ix2,ix3;
	static double r[98];
	double temp;
	static int iff=0;
	int j;
	if (*idum < 0 || iff==0)
	{ //2
		iff = 1;
		ix1 = (IC1-(*idum))%M1;
		ix1 = (IA1*ix1+IC1)%M1;
		ix2 = ix2%M2;
		ix1 = (IA1*ix1+IC1)%M1;
		ix3 = ix1 %M3;
		for(j=1; j<=97; j++)
		{ //3
			ix1 = (IA1*ix1+IC1)%M1;
			ix2 = (IA2*ix2+IC2)%M2;
			r[j]=(ix1+ix2*RM2)*RM1;
		} //3
		*idum = 1;
	} //2
	ix1 = (IA1*ix1+IC1)%M1;
	ix2 = (IA2*ix2+IC2)%M2;
	ix3 = (IA3*ix3+IC3)%M3;
	j= 1+((97*ix3)/M3);
	if (j>97 || j < 1)
		printf("impossible happened: j = %d\n",j);
	temp = r[j];
	r[j]=(ix1+ix2*RM2)*RM1;
	return temp;
} //1

//  End of sooper-dooper random number generator


/**************************************
*   This is a crummy Gauss function   *
*   used only for high values (>15)   *
*   of the Poisson distribution       *
**************************************/

double gauss(double mu, double sigma)
{ //1
	double t = -6;
	int j;
	for (j=0; j<12; j++)
	   t+= grandom(&idum);
	return sigma*t+mu;
} //1

/************************************************
* Poisson distribution, The parameter is a mean *
* that does not have to be an integer. The      *
* return is an integer from the distribution    *
* that corresponds to the mean.                 *
* This function might have been in the C        *
* library I had - but it wasn't.                *
************************************************/

int poisson (double lambda)
{ //1
	double elambda = exp(-lambda);
	double target = grandom(&idum);
	double term = elambda;
	double sum = term;
	int n = 0;
	if (lambda > 15)
	{ //2    // Use gaussian distribution instead

		return (int)(gauss(lambda, sqrt(lambda)));
	} //2
	while ( sum < target)
	{ //2
		term *= (lambda/(++n));
		sum+= term;
    } //2
    return n;
} //1


/*****************************************************
*   Numerical approximation to erf(x)  (Wikipedia)   *
*   Again, this function might have been in the C    *
*   library I had - but it wasn't.                   *
*****************************************************/

/**********************
* Constants for erf   *
**********************/

double eparams[10] = {-1.26551223, 1.00002368, 0.37409196,0.09678418,
                      -0.18628806, 0.27886807, -1.13520398,1.48851587,
                      -0.82215223,0.17087277};


double erf(double x)
{ //1
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
} //1


/*****************************************************
*  This function gives useful reports on MPI errors  *
*****************************************************/

void handleError(char * message, int r, int kk)
{ //1
   char error_string[100];
   int length_of_error_string;
   printf ("Core #%d %s\n",r,message);

   MPI_Error_string(kk, error_string, &length_of_error_string);
   printf( "Core%3d: %s\n", r, error_string);
   exit(0);
} //1
