/*************************************************
*  MPIstuff.c                                    *
*  This file has functions for transferring data *
*  between processes on ARCHI-WEST. In all cases *
*  process 0 is sending to all the others        *
*************************************************/


/********************************
* These variables are used to   *
* transfer a block of data from *
* core 0 to other cores, aiming *
* at different addresses        *
********************************/

MPI_Status status;
/**************************************
*  General termination message. s is  *
*  descriptive string and e is the    *
*  MPI error number                   *
**************************************/

void die(char * s, int e)
{
	fprintf(monitor, "Fatal fault: %s\n",s);
	printf("Fatal fault: see monitor.txt\n");
	fclose(monitor);
	MPI_Abort(MPI_COMM_WORLD,e);
}

/*********************************
* This function moves a block of *
* n doubles from core 0 to a     *
* memory block newly allocated   *
* in the other cores, and        *
* returns the address of this    *
* block.                         *
*********************************/


double * moveDoublesAcross(int n, double * s)
{
   int j,e;
   double buffer[500];
   if (myRank == 0)
   {
	  for (j = 0; j<n; j++)
	     buffer[j] = s[j];
   }
   e = MPI_Bcast(buffer, n, MPI_DOUBLE,0,MPI_COMM_WORLD);
   if (e != MPI_SUCCESS)
      die("MPI_Bcast fails in moveDoublesAcross",e);
   if(myRank == 0)
       return s;
   else
   {
	   double * v = (double *)(localMalloc(n*sizeof(double)));
	   for(j=0; j<n; j++)
	      v[j] = buffer[j];
       return v;
   }
}


/**************************************
* To broadcast a single int variable  *
**************************************/

void transferInt (int * z)
{
   int e = MPI_Bcast(z,1,MPI_INT,0,MPI_COMM_WORLD);
   if (e != MPI_SUCCESS)
      die("MPI_Bcast fails in transferInt",e);
}

/*************************************************
* To broadcast a an array of 6 double variables  *
*************************************************/

void transfer6DoubleArray (double * z)
{
   int e = MPI_Bcast(z,6,MPI_DOUBLE,0,MPI_COMM_WORLD);
   if (e != MPI_SUCCESS)
      die("MPI_Bcast fails in transferDoubleArray",e);
}

/*************************************************
* To broadcast a an array of 6 int variables     *
*************************************************/

void transfer6IntArray (int * z)
{
   int e = MPI_Bcast(z,6,MPI_INT,0,MPI_COMM_WORLD);
   if (e != MPI_SUCCESS)
      die("MPI_Bcast fails in transferIntArray",e);
}

/****************************************
* To broadcast a single double variable *
****************************************/

void transferDouble (double * z)
{
   int e = MPI_Bcast(z,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   if (e != MPI_SUCCESS)
      die("MPI_Bcast fails in transferDouble",e);
}

/*************************************************
* To broadcast an integer array to pos of length *
* len, where l is not generally known            *
*************************************************/

void transferIntArray( int * len, int * pos)
{
   int e =  MPI_Bcast(len,1,MPI_INT,0,MPI_COMM_WORLD);
   if (e != MPI_SUCCESS)
      die("MPI_Bcast fails in transferIntArray",e);

   totalTime += (MPI_Wtime() - startTime);
   MPI_Barrier(MPI_COMM_WORLD);
   startTime = MPI_Wtime();
   if (myRank != 0 )
   {
      pos = (int *)(localMalloc(sizeof(int)* (*len)));
      if (pos == NULL)
         die("localMalloc fails in transferInt Array",1);
   }
   e = MPI_Bcast(pos,*len,MPI_INT,0,MPI_COMM_WORLD);
      if (e != MPI_SUCCESS)
         die("MPI_Bcast fails in transferIntArray",e);
}

/*****************************************************
*  Next, consider Dimension structures, defined by : *
*                                                    *
*  typedef struct Dimension{                         *
*                     int resolution;                *
*                     double scatter;                *
*                     int scatterIndex;              *
*                     char scatterName[20];          *
*                  } Dimension;                      *
*                                                    *
*    We need to build an object of type MPI_Datatype *
*****************************************************/


MPI_Datatype MPI_dimension;
            // MPI variable for transferring dims structures

/*************************************************
* Call this function once to build an MPI        *
* variable for transferring dims records         *
* The record format is                           *
*                                                *
*  typedef struct Dimension{                     *
*                 int resolution;                *
*                 double scatter;                *
*                 int scatterIndex;              *
*                 char scatterName[20];          *
*                  } Dimension;                  *
*                                                *
*************************************************/

void Build_derived_dimension(MPI_Datatype * zzz)
{
    int block_length [4];
    int j;
    for (j = 0; j<4; j++)
       block_length[j] = 1;
    block_length[3] = 20;
    MPI_Datatype typeList[4];
    typeList[0] = typeList[2] = MPI_INT;
    typeList[1] = MPI_DOUBLE;
    typeList[3] = MPI_CHAR;
    MPI_Aint displacement[4];
    displacement[0] = 0;
    displacement[1] = sizeof(int);
    displacement[2] = displacement[1] + sizeof(double);
    displacement[3] = displacement[2] + sizeof(int);
    MPI_Type_struct(4,block_length,displacement,typeList, zzz);
    MPI_Type_commit(zzz);
}

MPI_Datatype MPI_list;
/*************************************************
* Call this function once to build an MPI        *
* variable for transferring ScatterRule records  *
* The record format is                           *
*                                                *
*  typedef struct ScatterRule{                   *
*	              char scatterName[20];          *
*	              int size;                      *
*	              double interval;               *
*	              double * values;               *
*			  } ScatterRule;                     *
*;                                               *
*  But we don't want to include the last com-    *
*  ponent as it is set by other means            *
*************************************************/

void Build_derived_list(MPI_Datatype * yyy)
{
    int block_length [3];
    block_length[0] = 20;
    block_length[1] = block_length[2] = 1;
    MPI_Datatype typeList[3];
    typeList[0] = MPI_CHAR;
    typeList[1] = MPI_INT;
    typeList[2] = MPI_DOUBLE;
    MPI_Aint displacement[3];
    displacement[0] = 0;
    displacement[1] = 20 * sizeof(char);
    displacement[2] = displacement[1] + sizeof(int);
    MPI_Type_struct(3,block_length,displacement,typeList, yyy);
    MPI_Type_commit(yyy);
}

MPI_Datatype MPI_cell;
            // MPI variable for transferring dims structures
            // I'll leave it here 'cos it's never used elsewhere

/*************************************************
* Call this function once to build an MPI        *
* variable for transferring cells                *
* The record format is                           *
*                                                *
*  typedef struct Cell{                          *
*                 long v;                        *
*                 double charge;                 *
*                  } Cell;                       *
*                                                *
*************************************************/

void Build_derived_cell(MPI_Datatype * xxx)
{
    int block_length [2];
    int j;
    for (j = 0; j<2; j++)
       block_length[j] = 1;
    MPI_Datatype typeList[2];
    typeList[0] = MPI_LONG;
    typeList[1] = MPI_FLOAT;
    MPI_Aint displacement[2];
    displacement[0] = 0;
    displacement[1] = sizeof(long);
    MPI_Type_struct(2,block_length,displacement,typeList, xxx);
    MPI_Type_commit(xxx);
}



void transferAllDimsRecords()
{
   for (j = 0; j<6; j++)
   {
      int e = MPI_Bcast(dims+j,1,MPI_dimension,0,MPI_COMM_WORLD);
         if (e != MPI_SUCCESS)
             die("MPI_BCast fails in transferAllDimsRecords",e);
   }
}


/*******************************************
* This function transfers all the stuff    *
* from the parameter file to all the cores *
* It is obeyed by all cores.               *
*******************************************/


void transferParameters()
{
    int j;
    transferInt( &spread);
    transferDouble (&threshold);
    transferInt (&nonoise);
    transferInt(&key);
    Build_derived_dimension (&MPI_dimension);
    transferAllDimsRecords();
    transferIntArray(&dataLen,data);
    transfer6DoubleArray(mins);
    transfer6DoubleArray(maxs);
    transfer6DoubleArray(range);
    transfer6DoubleArray(multiplier);
    transferInt(&rules);
    transfer6IntArray(stringSet);
    if (rules > 1)
        Build_derived_list(&MPI_list);
    for ( j = 1; j< rules; j++)
    {
	   int e = MPI_Bcast(list+j,1,MPI_list,0,MPI_COMM_WORLD);
	   if (e != MPI_SUCCESS)
	       die("MPI_Bcast fails in transferParameters",e);
       list[j].values = moveDoublesAcross(list[j].size+1,list[j].values);
    }
}

/****************************************
* This function distibutes the actual   *
* data among the cores. It is executed  *
* by all processes                      *
****************************************/

void transferData()
{
	 int j;
	 transferInt(&localNumberOfRecords);
     localv = (Record *) localMalloc(sizeof(Record) * localNumberOfRecords);
	 int e = MPI_Scatter(v,localNumberOfRecords*7,MPI_DOUBLE,
		                 localv, localNumberOfRecords*7,MPI_DOUBLE,
		                 0, MPI_COMM_WORLD);
	 if (e != MPI_SUCCESS)
	       die("MPI_Scatter fails in transferData",e);
     e = MPI_Bcast(coreBoundaries,np,MPI_DOUBLE,0,MPI_COMM_WORLD);
	 if (e != MPI_SUCCESS)
	       die("MPI_BCast fails in transferData",e);
}

