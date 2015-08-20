/*****************************************************
*  Program masp.c    for processing the output of a  *
*  pulse generator simulator, to adapt it for input  *
*  to a Free Electron Laser simulator.               *
*  This version is for the ARCHIE_WEST multi-        *
*  processor computer                                *
*  Author:             Andrew Colin                  *
*  User documentation: ElectronBeamConditioning.pdf  *
*  Date :              August 15st  2015             *
*  Version number:     0.04                          *
*  Latest mod:         See notes for 150316,150601,  *
*                      150801                        *
*  Note : all functions with double {{ and }}        *
*  brackets are executed by core 0 only.             *
*****************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <setjmp.h>
// #define DT 1   // To include distribution test
// #define TS 1   // To include tree statistics
#include "mpi.h"
#include "maspstructures.h"
#include "masp.h"
#include "masplibrary.c"
#include "maspparams.c"
#include "maspdatainput.c"
#include "MPIstuff.c"
#include "masptree.c"
#include "maspweight.c"
#include "masphisto.c"
#include "maspbunch.c"
#include "maspparticles.c"

/****************************************************
* This function sets up navigation for plot, which  *
* deals with splitting up a single fat electron.    *
* See "Distributing fat particles" to find out what *
* I am talking about, but basically, here it is:    *
* Every fat electron is aimed at a single cell, but *
* it gets divided into macro particles that appear  *
* in cells close by. The number of slots in any     *
* direction depends on the scatter constant of that *
* dimension.  If it is 1.0, then we assume that     *
* scattering can take place up to SPREAD units in   *
* either direction (like a Gaussian) so the number  *
* of slots is 2*SPREAD+1 (SPREAD on either side     *
* plus the middle one). If the scatter constant is  *
* 2.0, then the  number of slots will be 4*SPREAD+1 *
* The number must always be odd. It is given by the *
* formula                                           *
*   ((int)floor(2*scatter_constant*SPREAD+0.5)) | 1 *
* We also need the index of the middle slot (start- *
* ing from the first slot) to find the actual coor- *
* dinates of a specific dimension.                  *
* This happens for all six dimensions.              *
* When a fat electron comes along, aimed at a spec- *
* ific position within the target cell, the         *
* scattering function is used to calculate the      *
* probability that some of the charge will end up   *
* with a coordinate that corresponds to a slot in a *
* given row.                                        *
* The the probability of it occupying a given cell  *
* is just the product of the probabilities in each  *
* of the slots. Easy, innit??                       *
* Note: in this comment we put SPREAD to show it    *
* is a constant, but actually it's called 'spread'. *
* It is obeyed by all cores.                        *
****************************************************/

void setPlotSizes()
{{
	int j;
	for (j=0; j<DIMENSIONS; j++)
	{
		double q = dims[j].scatter;
		plotSizes[j] = ((int)floor(q*2*spread+0.5)) | 1;
		plotMidPoints[j] = plotSizes[j]>>1;
		plot [j] = (double *) (localMalloc(plotSizes[j] * sizeof(double)));
		if (plot[j] == NULL)
		{
			die("Not enough space for plot in setPlotSizes (very unlikely!)\n",88);
    	}
    }
}}



/**************************************
* This function sets the key axis,    *
* which is the axis with the highest  *
* resolution.                         *
**************************************/

void setKeyAxis()
{{
	int q=0,j;
	key = 0;
	for (j=0; j<6; j++)
	if  (dims[j].resolution > q)
	{
		q = dims[j].resolution;
			key = j;
	}
}}


/************************************************
* This is the initial sequence for core 0 only  *
************************************************/

void startCore0(char * pf)
{{   //1.5  // double brackets identify code exclusive to core 0

    	char nb[100];
    	int j;
    	monitor = fopen("monitor.txt","w");
    	if (monitor == NULL)
    	{  //2
			printf("Could not open monitor file\n"); // Can't use monitor file here!
			MPI_Abort(MPI_COMM_WORLD,4);
		}  //2
		params = fopen(pf,"r");
    	if (params == NULL)
    	{  //2
       		die("Could not open parameter.txt file. \n",3);
    	}  //2
   		fprintf(monitor, "Your parameter file is \"%s\".\n", pf);
   		j = readParameters(params);
        fprintf(monitor,"        -----   \n");

    	displayParameters();   // printf display (includes errors if any)
    	if (j)             // Temporary to test new parameters
    	{  //2
			die("There were errors in the parameters ",9);
    	}  //2

    	//  If no tasks identified, activate them all
    	if (stringSet[1] == 0 && stringSet[2] == 0 && stringSet[3] == 0 && stringSet[4] == 0)
    	{  //2
			stringSet[1]=1;stringSet[2]=1;stringSet[3]=1;stringSet[4]=1;
    	}  //2

    	in = fopen(files[0],"r");
		if (in == NULL)
		{  //2
			sprintf(nb,"Could not open input file %s \n", files[0]);
			die(nb,6);
		}  //2
        printf("Opened %s as input file\n",files[0]);

	    // Open all output files needed

	    for(j=1; j<= 4; j++)
		{
			if (stringSet[j])
		    {
			    out[j-1] = fopen(files[j],"w");
	            if (out[j-1] == NULL)
	            {
			       sprintf(nb,"Could not open output file %s \n", files[j]);
			       die(nb,6);
	            }
	             printf ("Opened %s as an output file. \n",files[j]);
		    }
	    }
        setKeyAxis();

		for (j = 0; j<100000; j++)
		{
			pt[j] = 0;
			pw[j] = 0;
		}

}}  // 1.5


/*********************************************
* This function read in and scales the data. *
*********************************************/

void initialiseData()
{{
		readData(in);
		scale();
                  // Set number of records sent to each core
		if (numberOfRecords % np)
		   localNumberOfRecords = numberOfRecords/np+1;
		else
		   localNumberOfRecords = numberOfRecords/np;
        // Set up coreBoundaries
        printf("Number of Records = %d Number Of Records per core = %d\n",numberOfRecords,localNumberOfRecords);
        for (j = 0; j < np; j++)
			coreBoundaries[j] = v[j*localNumberOfRecords].v[key];
}}

/*******************************************
* This function examines the scaled data   *
* and finds the range with the highest     *
* polulation in each of the six dimensions.*
* This is to get the best sample for the   *
* bunching estimates                       *
*******************************************/



void setHighest()
{ //1
	int * h;
	int k,z,q,a;
	for( k = 0; k < 6; k++)
	{  //2
	   h = (int *)(mem+soFar);   //  Reserve some space for a bit
	   z = dims[k].resolution;
	   for (q = 0; q<= z; q++)
	       h[q] = 0;
	   for (q=0; q < numberOfRecords; q++)
		   h[(int)(v[q].v[k])]++;
	   a=0; highest[k] = 0;
	   for (q = 0; q<=z; q++)
	   if (h[q] > a)
	   {  //3
	      a = h[q];
	      highest[k] = q;
       }  //3
    }//2
} //1









int main(int argc, char* argv[])
{  //1
   // For all cores
    int j, source;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);
    totalTime += (MPI_Wtime() - startTime);
    MPI_Barrier(MPI_COMM_WORLD);
    startTime = MPI_Wtime();
	start = MPI_Wtime();
	startTime = MPI_Wtime();
	totalTime = 0;

	JUMP = 10 * np;
    if (myRank == 0)

		printf("Number of threads  = %d\n",np);

	initialiseMemory();    // Organise local memory

	long qqq = time(&qqq);         // Fix random number generator
	idum = -30000+ (qqq&0x7FFF) +17 * myRank;
	   // depending on the day and different for each core

	dummy = grandom(&idum);        // This will still work after 2038 I think!

	if (myRank == 0)
	  startCore0(argv[1]);

    coreBoundaries = (double*) (localMalloc(np*sizeof(double)));  // This happens for all cores


	transferParameters();  // Transfers parameters to all cores
    setPlotSizes();    // for all cores

    if (myRank == 0)
    {{//2
		initialiseData();
		setHighest();
    }}  //2
   totalTime += (MPI_Wtime() - startTime);
   MPI_Barrier(MPI_COMM_WORLD);
   startTime = MPI_Wtime();
   if(stringSet[1])
   {
	   if (myRank == 0)
	      printf("Doing Macroparticles\n");
	   doMacroparticles();
   }
   if(stringSet[2]&& myRank == 0)
   {
		   printf("Doing Histogram\n");
		   doHistogram();
   }
   if(stringSet[3]&& myRank == 0)
   {
		   printf("Doing Cell weights\n");
		   doWeights();
   }

   if(stringSet[4]&& myRank == 0)
   {
	  printf ("Doing Bunching test\n");
	  doBunching();
   }

         if (myRank ==0)
   {{ //3
        fclose(out[0]);

        printf("Total input weight = %lf\n",totalInWeight);
        printf("Total output weight = %lf\n",totalOutWeight);
        printf("Percentage captured = %4.3lf\n",100*totalOutWeight/totalInWeight);
        fprintf(monitor,"Total input weight = %lf\n",totalInWeight);
        fprintf(monitor,"Total output weight = %lf\n",totalOutWeight);
        fprintf(monitor, "Percentage captured = %4.3lf\n",100*totalOutWeight/totalInWeight);
		printf(" Time with %d cores = %lf seconds\n", np, finish - start);
		fprintf(monitor, " Time with %d cores = %lf seconds\n", np, finish - start);

		fclose(in);
		fclose(monitor);
      }}  //3
   MPI_Finalize();

} //1

