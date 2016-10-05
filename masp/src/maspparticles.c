
/*************************************
* Maspparticles.c handles particles  *
* Version 0.7  (21/08/2016)          *
*************************************/


/**************************************
*  This function takes a cell with an *
*  expectation > threshold.           *
*  If noise is needed,the expectation *
*  is replaced by a Poisson variate.  *
*  If this variate is zero the cell   *
*  discarded. Otherwise it is output  *
*  using the MPI parallel-write       *
*  feature. All records must be the   *
*  same length ('cos the file layout  *
*  must decided in advance) and the   *
*  number of characters depends on    *
*  setting of the Decimals parameter. *
**************************************/

void writeRealCell (Cell v)
{   //1
	OutputRecord t;
    double bb[7];     // Final record for printing
    short ww[6];
    char dump [300];
    double s,shift;
    int j,k;
    expand(v.v,ww);         // Decode coordinates
    for (j=0; j<6; j++)     // Set up output record
       t.p.z[j] = ww[j]+0.5;
    t.weight = v.charge;
    if (!nonoise)          // Bring in the noise displacement
                           // Poisson function of weight already
                           // done if necessary
    { //2
	    s = 1.0/sqrt(t.weight);
        for(j=0; j<6; j++)  // Attach rectangular noise
           t.p.z[j] += s*(grandom(&idum)-0.5);
    } //2
    totalOutWeight += t.weight;

    for(j=0; j<6; j++)
    {  //2
       double dd = t.p.z[j];
       double ff = range[j]/dims[j].resolution;
       t.p.z[j] *= ff;
       t.p.z[j] += mins[j];
    }  //2
 		     //Set up the output record in the right order
	for(j=0; j<6; j++)
		  bb[j] = t.p.z[j];
	bb[6] = t.weight;
    sprintf(dump, format,
          bb[0],bb[1],bb[2],bb[3],bb[4],bb[5],bb[6],"\n");
          if (strlen(dump) != recordLength)
    {  //2
		  printf("Format = %s\n",format);
          printf("Dump=%s",dump);
          printf ("Record length = %d \n",(int)(strlen(dump)));
		  exit(1);
	}  //2

	j = MPI_File_write(mf, dump, recordLength,MPI_CHAR,&status);
	if (j != MPI_SUCCESS)
	   handleError("File write fails",myRank,j);
 }   //1


/******************************************
*  This function does the same job for an *
*  arbitrary scatter rule as erf for a    *
*  gaussian. It uses linear interpolation *
*  on cumulative values.                  *
******************************************/


double scatWeight (ScatterRule a, double val)
{  //1
   double t;
   if (val < -3)val = -3;
   if (val > 3) return 1.0;
      // Find boundary x-value
   double v = (3+val) / a.interval;
   int lowx = (int)(floor(v));
   double f = v-lowx;
   return a.values[lowx] * (1-f) + a.values[lowx+1] * f;
}  //1


/**************************************************
*  This function sets up the numbers in PLOT for  *
*  a gaussian or a private distribution. It uses  *
*  the actual position vector of each macro-      *
*  particle, and so must be executed for each one.*
**************************************************/
void setPlot(double positions[])
{  //1
    int j,k;
    double d;
    for (j=0; j<DIMENSIONS; j++)
    {  //2
       if (dims[j].scatter == 0)
		  plot[j][0] = 1.0;
	   else
	   {  //3
          double effectiveCellSize = (2.0*spread)/plotSizes[j];
 		  double effectivePos = positions[j]*effectiveCellSize;
		  k=0;double xx=0;
		  int ww=dims[j].scatterIndex;
		  for (d= -spread; d < spread+0.0001; d+= effectiveCellSize)
		  {  //4
			if (ww == 0)    // Gaussian distribution
			  plot[j][k++] = 0.5*(erf(d+effectiveCellSize-effectivePos) - erf(d+-effectivePos));
			else            // Private distribution
			  plot[j][k++]=  (scatWeight(list[ww], d+ effectiveCellSize-effectivePos)
			                      - scatWeight(list[ww],d-effectivePos));
		    xx+= plot[j][k-1];   // This is just a check!
		  }  //4
		  if (ww != 0 && (xx<0.9999 || xx > 1.0001))
		  {  //4
			 printf("Core #%d:Private scattering duzzent add up| xx = %lf\n", myRank, xx);
			 printf("This is not supposed to happen.\n");
			 exit(0);
		  }  //4

		  if(ww == 0)
		  {  //4
			 d=0.0;
			    // Apply very small correction to make sure area under
			    // the curve is 1.0000...
			 for(k=0; k<plotSizes[j]; k++)
			    d+= plot[j][k];
			 for (k=0; k<plotSizes[j]; k++)
			    plot[j][k] /= d;
		  }  //4
	  }  //3
  } //2
}  //1


/******************************************************
* This function is called for each record in each     *
* core. It splits the macroparticle into many micro-  *
* particles and sends each one to one of three trees: *
*  leftTree if the key coordinate is less than the    *
*  lower bound for that core                          *
*  rightTree if the key coordinate is more than the   *
*  upper bound for that core                          *
*  Otherwise the record goes to middleTree.           *
*  Vot iss dummy record? This is a dummy record with  *
*  impossible values, added to some cores to ensure   *
*  that they all have the same number of records.     *
******************************************************/

void processRecord(int j)
{ //1
	FILE * glob;    // Special
	int k,count,intwt;
	double dweight;
    double totwt = 0;
	int d[6];
	Record m = localv[j];
	Cell cc;
	short ww[6];
	Element * e;
    if (m.v[k] > 100000)
    {  //2
        printf("Dummy record caught\n");
        return;
    }   //2

    for (k = 0; k< 6; k++)
    {  //2
		celads[k] = (short) floor(m.v[k]);
		pos[k]  = m.v[k] - celads[k];
    } //2
	   setPlot(pos);
       //  The order here is staggered to produce a more balanced tree
    for (d[0] =0; d[0]<plotSizes[0]; d[0]++)
        for(d[1] = plotSizes[1]-1; d[1]>= 0; d[1]--)
           for (d[2] =0; d[2]<plotSizes[2]; d[2]++)
             for(d[3] = plotSizes[3]-1; d[3]>= 0; d[3]--)
                for (d[4] =0; d[4]<plotSizes[4]; d[4]++)
                  for(d[5] = plotSizes[5]-1; d[5]>= 0; d[5]--)
              { //2
   				  dweight = (m.v[6] * plot[5][d[5]]*plot[4][d[4]]
   				                    *plot[3][d[3]]*plot[2][d[2]]
   				                    *plot[1][d[1]]*plot[0][d[0]]);
                  if (dweight > threshold)
                  { //3
                     cc.charge = dweight;
                     for (k=0; k<6; k++)
						 ww[k] = celads[k]-plotMidPoints[k]+d[k];
                     cc.v = map(ww);
                     if (ww[key]< lowerBoundary)
                     { //4
						 treeBuild(&cc, &leftTree);
                         left++;
                     }  //4
                     else
                     if ( ww[key] >= upperBoundary)
                     { //4
						 treeBuild (&cc,&rightTree);
                         right++;
				     } //4
                     else
                     {  //4
						 treeBuild(&cc, &middleTree);
						 middle++;
				     } //4
			      } //3
		      } //2
	dun++;
	if (dun == (localNumberOfRecords/3))
	   printf ("Core %d has done one third of its job\n",myRank);
	if(dun == (2*localNumberOfRecords/3))
	   printf("Core %d has completed two thirds of the job\n",myRank);
}


/************************************
* This function splits incoming fat *
* particles into macroparticles.    *
* This is stage 1, where the data   *
* in each core is split among three *
* files.                            *
************************************/

void doMacroParticles1()
{ //1

	int j;
	int marginWidth;
	int soFarDump = soFar; // Where memory allocation starts
    middleTree=NULL;
    leftTree=NULL;
    rightTree = NULL;
    char nb[100];
    // Set boundaries
	if (myRank == 0)
	{  //2
		lowerBoundary = -100;
		upperBoundary = coreBoundaries[1];;
	}  //2
	else if (myRank == np-1)
	{  //2
		lowerBoundary = coreBoundaries[np-1];
		upperBoundary = dims[key].resolution + 100;
	}  //2
	else
	{  //2
		lowerBoundary = coreBoundaries[myRank];
		upperBoundary = coreBoundaries[myRank+1];
	}  //2

	     // Start processing records
    left =0; middle = 0; right = 0;
	for(j = 0; j < localNumberOfRecords; j++)
	if(localv[j].v[0] < 999999)   // Don't pocess dummy records
	   processRecord(j);
	printf("Core #%d left = %d, middle = %d, right = %d\n",myRank,left,middle,right);
}  //1


/********************************************
* This function explores the selected side  *
* tree and moves all the cells to the seam  *
* array in random order                     *
********************************************/

void explore (Element * h)
{ //1
	if (grandom(&idum) > 0.5)
	{  //2
		 if (h->back != NULL)
		    explore (h->back);
		 inSeam[inSeamLength++] = h->cell;
		 if (h->fore != NULL)
		    explore (h->fore);
    } //2
    else
 	{  //2
		 if (h->fore != NULL)
		    explore (h->fore);
		 inSeam[inSeamLength++] = h->cell;
		 if (h->back != NULL)
		    explore (h->back);
    }  //2
}  //1


/***********************************
* This is the second stage, where  *
* the marginal overflows from each *
* pair of cores is sewn together   *
* Special care needed when the size*
* of the seam is zero.             *
***********************************/

void doMacroParticles2()
{ //1
	Build_derived_cell(&MPI_cell);
	inSeam = (Cell *) (mem+soFar);
	inSeamLength = 0;
	int j;
	if(myRank != 0)
	{  //2
		                       // For all cores except core 0:
		if (left > 0)
		    explore(leftTree);
	                           // Initialise outSeam when inSeam is aleady full
                               // Send length of cell table to Core (myRank-1)
        MPI_Send(&inSeamLength,1, MPI_LONG,myRank-1,0, MPI_COMM_WORLD);
    }  //2

    if (myRank < (np-1)) // For all cores except core (np-1)
    {  //2
		MPI_Recv(&outSeamLength,1,MPI_LONG,myRank+1,0,MPI_COMM_WORLD,&status);
		outSeam = &(inSeam[inSeamLength]);
    }  //2

    if(myRank != 0)      // For all cores except 0: send data
    {  //2
		if(inSeamLength > 0)
			MPI_Send(inSeam,inSeamLength,MPI_cell,myRank-1,0, MPI_COMM_WORLD);
    }  //2
	if(myRank < (np-1))   // For all cores except cor(n-1): receive data
	{ // 2                   // and populate the middle tree
	     if(outSeamLength > 0)
	     {  //3
             MPI_Recv((void *)(outSeam),outSeamLength,MPI_cell,myRank+1, 0,MPI_COMM_WORLD,&status);
		     for(j=0; j < outSeamLength; j++)
		        treeBuild(outSeam+j, &middleTree);
	     }  //3
     }  //2
	                        // Now do it all the other way!
	inSeam = (Cell *) (mem+soFar);
	inSeamLength = 0;
	if(myRank < (np-1))
	{  //2
		                   // For all cores except core np-1:
		if (right > 0)
		    explore(rightTree);
	                       // Initialise outSeam when inSeam is aleady full
	    outSeam = (Cell *)(inSeam + inSeamLength);
                           // Send length of cell table to Core (myRank+1)
        MPI_Send(&inSeamLength,1, MPI_LONG,myRank+1,0, MPI_COMM_WORLD);
    } //2

    if (myRank != 0)      // For all cores except core 0)
    {  //2
       MPI_Recv(&outSeamLength,1,MPI_LONG,myRank-1,0,MPI_COMM_WORLD,&status);
       outSeam = &(inSeam[inSeamLength]);
    }  //2
    if(myRank < (np-1))      // For all cores except core np-1: send data
        if(inSeamLength > 0)
            MPI_Send(inSeam,inSeamLength,MPI_cell,myRank+1,0, MPI_COMM_WORLD);

	if(myRank != 0)          // For all cores except core 0 receive data
	 { //2                   // and populate the middle tree
		 if (outSeamLength > 0)
		 {  //3
			 MPI_Recv(outSeam,outSeamLength,MPI_cell, myRank-1,0,MPI_COMM_WORLD,&status);
		     for(j=0; j < outSeamLength; j++)
		        treeBuild(outSeam+j, &middleTree);
	     }  //3
     } //2
}//1

/********************************************
* This function explores the main tree.     *
* If the cell's charge expectation is more  *
* than the THRESHOLD it appends the cell to *
* inSeam.                                   *
* Important modification made on 150904:    *
* If noise is to be used, then the Poisson  *
* adjustment to the charge is made here     *
* and cells with zero charge are discarded. *
* The coordinates are not altered at this   *
* stage.                                    *
********************************************/

void finalOut (Element * h)
{  //1
	 if (h->back != NULL)
		    finalOut (h->back);
	 if (!nonoise)
	 { //2
		 int k = poisson(h->cell.charge);
		 if (k > 0)
		 {  //3
			 h->cell.charge = (double)k;
			 inSeam[inSeamLength++] = h->cell;
		 }  //3
	 }  //2
	 else
	 if (h->cell.charge > threshold)
		 inSeam[inSeamLength++] = h->cell;
	 if (h->fore != NULL)
		finalOut (h->fore);
}  //1


/*************************************
* (Comment entered a bit later)      *
* I THINK this function collects all *
* the records in the middle tree and *
* outputs them.                      *
*************************************/

void collectCells()
{  //1
	inSeam = (Cell *) (mem+soFar);
	inSeamLength = 0;
    if (middle > 0)
    finalOut(middleTree);
}  //1

/***************************************
* This is where all the work gets dun! *
***************************************/

void doMacroparticles()
{  //1
   int source,kk;
   MPI_Datatype  ptype;
   {  //2
    // Do MACROPARTICLES
       {  //3
	   transferData();
	   doMacroParticles1();
	   MPI_Barrier(MPI_COMM_WORLD);  // Synchronise
       doMacroParticles2();
       }  //3
       totalTime += (MPI_Wtime() - startTime);
       startTime = MPI_Wtime();
       collectCells();
       totalTime += (MPI_Wtime() - startTime);
       MPI_Barrier(MPI_COMM_WORLD);
       startTime = MPI_Wtime();
          setRoots();
       if (myRank == 0)
       { //3
           for (j=0; j<np; j++)
              printf ("Root [%d] = %ld\n",j,ww[j]);
       } //3

       MPI_Type_contiguous (recordLength,MPI_CHAR,&ptype);
       MPI_Type_commit (&ptype);
	   kk = MPI_File_open(MPI_COMM_WORLD, microfile,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&mf);
       if (kk != MPI_SUCCESS)
		   handleError("MPI_File_open fails",myRank,kk);
	   if (myRank == 0)
	       printf("Opened microparticle file %s\n",microfile);
	   kk=  MPI_File_set_view(mf,recordLength*root,ptype,ptype,"native",MPI_INFO_NULL);
       if (kk != MPI_SUCCESS)
          handleError("MPI_File_set_view fails",myRank,kk);
       printf("Core #%d  About to start writing %ld records\n",myRank,inSeamLength);
       for(kk = 0; kk<inSeamLength; kk++)
             writeRealCell(inSeam[kk]);
       printf("Core %d finished writing\n",myRank);
       kk= MPI_File_close(&mf);
       if (kk != MPI_SUCCESS)
		   handleError("MPI_File_close fails",myRank,kk);
       totalTime += (MPI_Wtime() - startTime);
       MPI_Barrier(MPI_COMM_WORLD);
       startTime = MPI_Wtime();
	   finish = MPI_Wtime();

   }  //2
} //1

