/**************************************
*  This function takes a cell with an *
*  expectation > threshold.           *
*  If noise is needed,the expectation *
*  is replaced by a Poisson variate.  *
*  If this variate is zero the cell   *
*  discarded
**************************************/

void writeRealCell (Cell v)
{{   //1
	OutputRecord t;
    float bb[7];     // Final record for printing
    short ww[6];
    double s,shift;
    int j,k,OK;
    if (stringSet[4] == 1)  // Decide if this cell is elegible
       OK = elegible (ww);  // for inclusion in the bunching statistics
    else OK = -1;
    expand(v.v,ww);         // Decode coordinates
    for (j=0; j<6; j++)     // Set up output record
       t.p.z[j] = ww[j]+0.5;
    t.weight = v.charge;
    if (!nonoise)          // Bring in the noise
    { //2
		k = poisson(t.weight);
		if(k == 0)         // Ignore cells with zero weight
		   return;
        t.weight = (double) k;
        s = 1.0/sqrt(t.weight);
        for(j=0; j<6; j++)  // Attach rectangular noise
           t.p.z[j] += s*(grandom(&idum)-0.5);
    } //2
    if(OK >=0)    // Keep these values for bunching check
    {  //2
		pt[OK+ww[key]] = t.p.z[key];
		pw[OK+ww[key]] = t.weight;
    }  //2
    totalOutWeight += t.weight;
    for( j=0; j<6; j++)
       t.p.z[j] = mins[j] + t.p.z[j] *(range[j]/dims[j].resolution);
		     //Set up the output record in the right order
	for(j=0; j<7; j++)
	{ //2
		k = output[j];
		if (k == 6)
		  bb[j] = t.weight;
		else
		  bb[j] = t.p.z[k];
	} //2
	if (binary)
	    fwrite((void *) bb,sizeof(float),7,out[0]);
    else
        fprintf(out[0],
          "%10.5lg,%10.5lg,%10.5lg,%10.5lg,%10.5lg,%10.5lg,%10.5lg\n",
          bb[0],bb[1],bb[2],bb[3],bb[4],bb[5],bb[6]);
}}   //1


/******************************************
*  This function does the same job for an *
*  arbitrary scatter rule as erf for a    *
*  gaussian. It uses linear interpolation *
*  on cumulative values.                  *
******************************************/


double scatWeight (ScatterRule a, double val)
{
   double t;
   if (val < -3)val = -3;
   if (val > 3) return 1.0;
      // Find boundary x-value
   double v = (3+val) / a.interval;
   int lowx = (int)(floor(v));
   double f = v-lowx;
   return a.values[lowx] * (1-f) + a.values[lowx+1] * f;
}





/**************************************************
*  This function sets up the numbers in PLOT for  *
*  a gaussian or a private distribution. It uses  *
*  the actual position vector of each macro-      *
*  particle, and so must be executed for each one.*
**************************************************/
void setPlot(double positions[])
{
    int j,k;
    double d;
    for (j=0; j<DIMENSIONS; j++)
    {
           if (dims[j].scatter == 0)
		    plot[j][0] = 1.0;
		else
		{
            double effectiveCellSize = (2.0*spread)/plotSizes[j];

 			double effectivePos = positions[j]*effectiveCellSize;
			k=0;double xx=0;
		    int ww=dims[j].scatterIndex;
			for (d= -spread; d < spread+0.0001; d+= effectiveCellSize)
			{
				if (ww == 0)    // Gaussian distribution
				   plot[j][k++] = 0.5*(erf(d+effectiveCellSize-effectivePos) - erf(d+-effectivePos));
			    else            // Private distribution
			       plot[j][k++]=  (scatWeight(list[ww], d+ effectiveCellSize-effectivePos)
			                      - scatWeight(list[ww],d-effectivePos));
		        xx+= plot[j][k-1];   // This is just a check!
		    }
			if (ww != 0 && (xx<0.9999 || xx > 1.0001))
			{
			    printf("Core #%d:Private scattering duzzent add up| xx = %lf\n", myRank, xx);
				printf("This is not supposed to happen.\n");
				exit(0);
		    }

			// Temporary bit
			if(ww == 0)
			{
				d=0.0;
			// Apply very small correction
			   for(k=0; k<plotSizes[j]; k++)
			     d+= plot[j][k];
			   for (k=0; k<plotSizes[j]; k++)
			     plot[j][k] /= d;
		    }
	  }
  }
}







/******************************************************
* This function is called for each record in each     *
* core. It splits the macroparticle into many micro-  *
* particles and sends each one to one of three trees: *
*  leftTree if the key coordinate is less than the    *
*  lower bound for that core                          *
*  rightTree if the key coordinate is more than the   *
*  upper bound for that core                          *
*  Otherwise the record goes to middleTree.           *
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

              { //3
				  dweight = (m.v[6] * plot[5][d[5]]*plot[4][d[4]]*plot[3][d[3]]*plot[2][d[2]]*plot[1][d[1]]*plot[0][d[0]]);
                  if (dweight > threshold)
                  { //3
                        cc.charge = dweight;
 //                    printf("Core # %d: key variable: %lf\n",myRank,cc.v[key]);
                     for (k=0; k<6; k++)
						 ww[k] = celads[k]-plotMidPoints[k]+d[k];
                     cc.v = map(ww);
 #ifdef DT
                     if (myRank == 1  && j == 100)   // Special
                     { //4
                         fprintf(glob,"%d %d %d %d %d %d, %f\n",ww[0],ww[1],ww[2],ww[3],ww[4],ww[5],cc.charge);
                         totwt += cc.charge;
				     } //4
 #endif
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
	if ((dun%500) == 0)
	   printf ("Core %d has done %d records\n",myRank,dun);
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
    aa=bb=cc=dd=0;

	treeCount = 0;
    char nb[100];
    // Set boundaries
		if (myRank == 0)
		{
			lowerBoundary = -100;
		    upperBoundary = coreBoundaries[1];;
	    }
	    else if (myRank == np-1)
	    {
			lowerBoundary = coreBoundaries[np-1];
			upperBoundary = dims[key].resolution + 100;
	    }
	    else
	    {
			lowerBoundary = coreBoundaries[myRank];
			upperBoundary = coreBoundaries[myRank+1];
		}

	     // Start processing records
    left =0; middle = 0; right = 0;
	for(j = 0; j < localNumberOfRecords; j++)
	if(localv[j].v[0] < 999999)   // Don't pocess dummy records
	   processRecord(j);
}




/********************************************
* This function explores the selected side  *
* tree and moves all the cells to the seam  *
* array in random order                     *
********************************************/

void explore (Element * h)
{

	if (grandom(&idum) > 0.5)
	{
		 if (h->back != NULL)
		    explore (h->back);
		 inSeam[inSeamLength++] = h->cell;
		 if (h->fore != NULL)
		    explore (h->fore);
    }
    else
 	{
		 if (h->fore != NULL)
		    explore (h->fore);
		 inSeam[inSeamLength++] = h->cell;
		 if (h->back != NULL)
		    explore (h->back);
    }

}


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
        MPI_Send(&inSeamLength,1, MPI_INT,myRank-1,0, MPI_COMM_WORLD);
    }  //2

    if (myRank < (np-1)) // For all cores except core (np-1)
    {
		MPI_Recv(&outSeamLength,1,MPI_INT,myRank+1,0,MPI_COMM_WORLD,&status);
		outSeam = &(inSeam[inSeamLength]);
    }

    if(myRank != 0)      // For all cores except 0: send data
    {
		if(inSeamLength > 0)
		{
			MPI_Send(inSeam,inSeamLength,MPI_cell,myRank-1,0, MPI_COMM_WORLD);
	    }
    }
	if(myRank < (np-1))   // For all cores except cor(n-1): receive data
	{ // 2                   // and populate the middle tree
	     if(outSeamLength > 0)
	     {
             MPI_Recv((void *)(outSeam),outSeamLength,MPI_cell,myRank+1, 0,MPI_COMM_WORLD,&status);
		     for(j=0; j < outSeamLength; j++)
		        treeBuild(outSeam+j, &middleTree);
	     }
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
        MPI_Send(&inSeamLength,1, MPI_INT,myRank+1,0, MPI_COMM_WORLD);
    } //2

    if (myRank != 0) // For all cores except core 0)
    {
       MPI_Recv(&outSeamLength,1,MPI_INT,myRank-1,0,MPI_COMM_WORLD,&status);
       outSeam = &(inSeam[inSeamLength]);
    }
    if(myRank < (np-1))      // For all cores except core np-1: send data
        if(inSeamLength > 0)
            MPI_Send(inSeam,inSeamLength,MPI_cell,myRank+1,0, MPI_COMM_WORLD);

	if(myRank != 0)   // For all cores except core 0 receive data
	 { //2                   // and populate the middle tree
		 if (outSeamLength > 0)
		 {
			 MPI_Recv(outSeam,outSeamLength,MPI_cell, myRank-1,0,MPI_COMM_WORLD,&status);
		     for(j=0; j < outSeamLength; j++)
		        treeBuild(outSeam+j, &middleTree);
	     }
     } //2
}//1

/********************************************
* This function explores the main tree.     *
* If the cell's charge expectation is more  *
* than the THRESHOLD it appends the cell to *
* inSeam.                                   *
********************************************/

void finalOut (Element * h)
{
	 if (h->back != NULL)
		    finalOut (h->back);
		 if (h->cell.charge > threshold)
		     inSeam[inSeamLength++] = h->cell;
		 if (h->fore != NULL)
		    finalOut (h->fore);
}

void collectCells()
{
	inSeam = (Cell *) (mem+soFar);
	inSeamLength = 0;
    if (middle > 0)
    finalOut(middleTree);
}

void doMacroparticles()
{
   int source;
   {  //2
    // Do MACROPARTICLES
       {  //3
	   transferData();
	   stage = 1;
	   doMacroParticles1();
	   MPI_Barrier(MPI_COMM_WORLD);  // Synchronise
       doMacroParticles2();
       }  //3
       totalTime += (MPI_Wtime() - startTime);
       MPI_Barrier(MPI_COMM_WORLD);
       startTime = MPI_Wtime();
       collectCells();
       totalTime += (MPI_Wtime() - startTime);
       MPI_Barrier(MPI_COMM_WORLD);
       startTime = MPI_Wtime();
       if (myRank == 0)
       { //3
          for(j = 0; j < inSeamLength; j++)
		     writeRealCell(inSeam[j]);
		  for (source = 1; source < np; source++)
		  { //4
		     MPI_Recv(&inSeamLength,1,MPI_INT, source,0,MPI_COMM_WORLD, &status);
             MPI_Recv(inSeam, inSeamLength, MPI_cell, source,0,MPI_COMM_WORLD, &status);
		     for(j = 0; j < inSeamLength; j++)
		         writeRealCell(inSeam[j]);
	      } //4
	  }  //3
	  else
	  {  //3
		  MPI_Send(&inSeamLength,1,MPI_INT,0,0,MPI_COMM_WORLD);
		  MPI_Send(inSeam,inSeamLength,MPI_cell,0,0,MPI_COMM_WORLD);
      }  //3

      totalTime += (MPI_Wtime() - startTime);
      MPI_Barrier(MPI_COMM_WORLD);
      startTime = MPI_Wtime();
	  finish = MPI_Wtime();

   }  //2
} //1

