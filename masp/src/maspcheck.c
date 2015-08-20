/*****************************************************
*  Program maspcheck.c                               *
*  This program is used to check the parameter file  *
*  for MASP, and to give estimates of the resources  *
*  needed for a run on ARCHIE-WEST                   *
*  Author:             Andrew Colin                  *
*  User documentation: maspspec2.pdf                 *
*  Date :              August 3rd  2015              *
*  Version number:     0.02                          *
*****************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maspstructures.h"
#include "masp.h"
#include "masplibrary.c"
#include "maspparams.c"
#include "maspcheckdatainput.c"

void die(char * s, int e)
{
	fprintf(monitor, "Fatal fault: %s\n",s);
	printf("Fatal fault: see monitor.txt\n");
	fclose(monitor);
}

int main(int argc, char* argv[])
{  //1
  	int kk = 0;
  	int q, numberOfCores;
  	if (argc!= 2)
  	{
		printf("Wrong number of parameters\n");
		exit(1);
    }
    params = fopen(argv[1],"r");
    JUMP = 1000;
    if (params == NULL)
    {  //2
       printf("Could not find the parameter file.\n");
       exit(1);
    }  //2
    monitor = fopen("monitor.txt","w");
    if (monitor == NULL)
    {  //2
		printf("Could not open monitor file\n");
		exit(1);
	}  //2
    initialiseMemory();
    printf("Please look in file \"monitor.txt\" for a report\n");
    j = readParameters(params);
    displayParameters();
    if (j)
    {  //2
		printf("There were errors in the parameters (see monitor.txt)\n");
		fclose(monitor);
		exit(1);
    }  //2
    //  If no tasks identified, activate them all
    if (stringSet[1] == 0 && stringSet[2] == 0 && stringSet[3] == 0)
    {  //2
		stringSet[1]=1;stringSet[2]=1;stringSet[3]=1;;
    }  //2

    in = fopen(files[0],"r");
	if (in == NULL)
	{  //2
		printf("Could not open input file %s \n", files[0]);
		fprintf(monitor,"Could not open input file %s \n", files[0]);
		fclose(monitor);
		exit(0);
	}  //2
	for(j = 0; j<3; j++)
	{  //2
		if (stringSet[j+1])
		{  //3
		    out[j] = fopen(files[j+1],"w");
	        if (out[j] == NULL)
	        {  //4
		       kk = 1;
		       printf("File %s already open in another program!\n", files[j+1]);
		       fprintf(monitor, "File %s already open in another program!\n", files[j+1]);
	       } //4
	    }  //3
    }  //2
    if (kk)
        fprintf(monitor, "Make sure your output files are usable\n");

    readData (in);// Read data and sort on the k'th value
    scale();
    q = numberOfRecords / 3200;
    numberOfCores = ((q+11)/12)*12;
    fprintf(monitor, "Run MASP using at least %d cores\n", numberOfCores);
    fprintf(monitor, "Approximate run-time estimate = %lf seconds \n",(numberOfRecords/38000.0) * 12*108.0/numberOfCores);
    fprintf(monitor, "This will be more if you lower the Threshold, or increase the number of cells in any dimension\n");
    fprintf(monitor, "If you double the number of cores, there will be some speed gain,\n");
    fprintf(monitor, "but beyond that, don't bother. There will be almost no advantage.\n");
    printf("\n\nNo errors found in the parameter file.\n");
    printf("\nLook in monitor.txt for a full list of \nparameters (including defaults)\n");
    fclose(monitor);
    return 1;

} //1
