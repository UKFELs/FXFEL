/*****************************************************
*  Program masp.c    for processing the output of a  *
*  pulse generator simulator, to adapt it for input  *
*  to a Free Electron Laser simulator.               *
*  This version is for the ARCHIE_WEST multi-        *
*  processor computer                                *
*  Author:             Andrew Colin                  *
*  User documentation: The FELFX Simulator     .pdf  *
*  Date :              August 21st 2016              *
*  Version number:     0.7                           *
*  Latest mod:         See notes for 150316,150601,  *
*                      150801                        *
**  This version allows control of the number of    **
** significant figures in numbers written to the    **
** output file. Put "Decimals = n" (for 2<=n<=15)   **
** in the parameter file.                           **
*  Note : all functions with double {{ and }}        *
*  brackets are executed by core 0 only. Each double *
*  bracket counts as one-and-a-half brackets.        *
*****************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <setjmp.h>
#include "mpi.h"
#include "maspstructures.h"
#include "masp.h"
#include "masplibrary.c"
#include "maspparams.c"
#include "maspdatainput.c"
#include "MPIstuff.c"
#include "masptree.c"
#include "maspparticles.c"

/****************************************************
* This function sets up navigation for plot, which  *
* deals with splitting up a single fat electron.    *
* See "Distributing fat particles" to find out what *
* I am talking about, but basically, here it is:    *
* Every fat electron is aimed at a single cell, but *
* it gets divided into micro particles that appear  *
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
{{  //1.5
	int j;
	for (j=0; j<DIMENSIONS; j++)
	{  //2
		double q = dims[j].scatter;
		plotSizes[j] = ((int)floor(q*2*spread+0.5)) | 1;
		plotMidPoints[j] = plotSizes[j]>>1;
		plot [j] = (double *) (localMalloc(plotSizes[j] * sizeof(double)));
		if (plot[j] == NULL)
			die("Not enough space for plot in setPlotSizes (very unlikely!)\n",88);
    }  //2
}}  //1.5



/**************************************
* This function sets the key axis,    *
* which is the axis with the highest  *
* resolution.                         *
**************************************/

void setKeyAxis()
{{  //1.5
	int q=0,j;
	key = 0;
	for (j=0; j<6; j++)
	if  (dims[j].resolution > q)
	{ //2
		q = dims[j].resolution;
			key = j;
	}  //2
}}  //1.5


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
    	if (stringSet[1] == 0)
    	{  //2
			stringSet[1]=1;
    	}  //2


        in = fopen(files[0],"r");
		if (in == NULL)
		{  //2
			sprintf(nb,"Could not open input file %s \n", files[0]);
			die(nb,6);
		}  //2
        printf("Opened %s as input file\n",files[0]);

        if (stringSet[1])
        {  //2
			       // Transfer microparticle output file name to all cores
			if (strlen(files[1]) > 38)
			{  //3
				sprintf(nb, "Microparticle file name too long :%s\n",files[1]);
				die(nb,777);
		    }  //3
				strcpy(microfile,files[1]);
	    }  //2
	    else
	      strcpy(microfile, "dummy");
          j=1;
		{  //2
			if (stringSet[j])
		    {  //3
			    out = fopen(files[j],"w");
	            if (out == NULL)
	            {  //4
			       sprintf(nb,"Could not open output file %s \n", files[j]);
			       die(nb,6);
	            } //4
	             printf ("Opened %s as an output file. \n",files[j]);
		    } //3
	    }  //2
        setKeyAxis();
}}  // 1.5


/*********************************************
* This function read in and scales the data. *
*********************************************/

void initialiseData()
{{  //1.5
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
}}  //1.5

/***************************
* Where it all starts ...  *
***************************/


int main(int argc, char* argv[])
{  //1
   // For all cores
    int j, source,jj,kk;
    double grandTotalOut;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

	if (np > MAX_CORES)
	{  //2
		printf("Too many cores. Maximum is%d\n",MAX_CORES);
		exit(1);
    }  //2
    totalTime += (MPI_Wtime() - startTime);
    MPI_Barrier(MPI_COMM_WORLD);
    startTime = MPI_Wtime();
	start = MPI_Wtime();
	startTime = MPI_Wtime();
	totalTime = 0;

	JUMP = 10 * np;
    if (myRank == 0)
		printf("Number of threads  = %d Myrank = %d\n",np,myRank);

	initialiseMemory();    // Organise local memory

	long qqq = time(&qqq);         // Fix random number generator
	idum = -30000+ (qqq&0x7FFF) +17 * myRank;
	   // depending on the day and different for each core

	dummy = grandom(&idum);        // This will still work after 2038 I think!

	if (myRank == 0)
	startCore0(argv[1]);
	                     // Next bit must happen after parameters have been read
    jj = decimals+7;     // Set up format for output
    kk = (decimals%4)+1; // The number of characters must be divisible by 4
    sprintf(format,"%%+%d.%dle,%%+%d.%dle,%%+%d.%dle,%%+%d.%dle,%%+%d.%dle,%%+%d.%dle,%%+%d.%dle%%%ds",
            jj,decimals,jj,decimals,jj,decimals,jj,decimals,jj,decimals,jj,decimals,jj,decimals,kk);
    recordLength = (7*jj+kk+6);
    coreBoundaries = (double*) (localMalloc(np*sizeof(double)));  // This happens for all cores

   transferParameters();  // Transfers parameters to all cores


    setPlotSizes();    // for all cores

    if (myRank == 0)
    {{//2.5
		initialiseData();
    }}  //2.5
   totalTime += (MPI_Wtime() - startTime);
   MPI_Barrier(MPI_COMM_WORLD);
   startTime = MPI_Wtime();
   if(stringSet[1])
   { //2
	   if (myRank == 0)
	      printf("Doing Macroparticles\n");
	   doMacroparticles();
   }  //2
   MPI_Reduce (&totalOutWeight,&grandTotalOut,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
         if (myRank ==0)
   {{ //2.5

        printf("Total input weight = %lf\n",totalInWeight);
        printf("Total output weight = %lf\n",grandTotalOut);
        printf("Percentage captured = %4.3lf\n",100*grandTotalOut/totalInWeight);
        fprintf(monitor,"Total input weight = %lf\n",totalInWeight);
        fprintf(monitor,"Total output weight = %lf\n",grandTotalOut);
        fprintf(monitor, "Percentage captured = %4.3lf\n",100*grandTotalOut/totalInWeight);
		printf(" Time with %d cores = %lf seconds\n", np, finish - start);
		fprintf(monitor, " Time with %d cores = %lf seconds\n", np, finish - start);

		fclose(in);
		fclose(monitor);
      }}  //2.5
   MPI_Finalize();
} //1

