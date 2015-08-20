/*****************************************
*  maspdatainput.c  Input, randomisation *
*  and scaling of raw data               *
*****************************************/


/******************************************
* Comparison function for sorting         *
******************************************/

int comparison(const void * a, const void * b)
{
	double q;
	Record *x = (Record *) a;
	Record *y = (Record *) b;
	q= (x->v[key] -y->v[key]);
	if (q<0) return (-1);
	else if (q>0)  return 1;
    else return 0;
}



/****************************************
*  This auxiliary function scans a line *
*  of characters in *line, starting at  *
*  line[p]. It returns the next number  *
*  in *zz, and an index to the rest of  *
*  the line. If there are no more char- *
*  acters the funtion returns -1.       *
****************************************/

int getNext( int p, char * zz, char * line)
{
	int q = 0, s=p;
	while (line[s] == ' ')
	   s++;
	if (line[s] == 0)
	   return (-1);
	while (line[s] != ' ' && line[s] != 0)
	   zz[q++] = line[s++];
	zz[q] = 0;
	return s;
}

/****************************************
* This function reads data into array v *
* and sorts on the key axis. The key    *
* axis is the one with the highest res- *
* olution, as determined by the input   *
* parameters.                           *
****************************************/

int readData( FILE * fs)
{  //1
	int k,b,r;
	char nb[100];   // Buffer for error reports
	double dummy;
	char  line[400];
	char  next[20];
	char * zz;
	Record w;
    int p;
    int n=0,kk;
    Record dummyRecord;
    int qq;
    v = (Record *)(malloc(sizeof(Record) * JUMP));
    if (v == NULL)
    { //2
		die("malloc failed in readData",7);
	}  //2
	numberOfRecords = 0;
    vSize = JUMP;
    while(1)
    {  //2
      zz = fgets(line,399,fs);
      if (zz == NULL)
      {  //3
		  if (feof(fs))
		    break;
		  else
		  { //4
			  sprintf(nb,"Input failed after %d records\n",numberOfRecords);
			  die(nb,77);
	      } //4
       } //3
       p=0;kk=0;n=0;
       while (1)
       { //3
		  p = getNext(p,next,line);
          if (p < 0)
		  { //4
			  sprintf(nb,"Not enough data in record %d\n",numberOfRecords+1);
			  die(nb,9);

	      } //4
		  r = sscanf(next, "%le", &dummy);
		  if (r != 1)
		  {  //4
			  sprintf(nb,"Bad data in record %d\n",numberOfRecords+1);
              die (nb,10);
	      }  //4
	      n++;
	      for(j=0; j<7; j++)
	      if (n-1 == data[j])
	      { //4
			  w.v[j] = dummy;
	          kk++;
	      }  //4
          if(kk==7)

				break;
       } //3
     if (numberOfRecords == vSize)
     { //3          // Just grab more space (knowing it will be contiguous)
		 Record * qqq = (Record *)(malloc(sizeof(Record) * JUMP));
		 if (qqq == NULL)
		 { //4
			 sprintf(nb,"malloc failed in readData after %d records\n",vSize);
			 die(nb,76);
	     }  //4
	     vSize += JUMP;
     }  //3
     totalInWeight += w.v[6];
     v[numberOfRecords++] = w;
  } //2
  // Now sort on the key axis
  fprintf(monitor,"Number of records = %d\n",numberOfRecords);
  qsort(v,numberOfRecords,sizeof(Record),comparison);

    // Pad out the buffer with irrelevant records
/*  for (qq = 0; qq<7; qq++)
     dummyRecord.v[qq] = 1000000;  // All fields are a million
  qq = numberOfRecords;
  while (qq < vSize)
      v[qq++] = dummyRecord;
*/
      // Now the total number of records is an exact multiple
      // of np, the number or processes.
  return numberOfRecords;
} //1

/****************************************
* This function finds the max and min   *
* of each entry in the Record, and      *
* scales the entries to the range 0 to  *
* dims[j].resolution for each dimension *
****************************************/

  void scale ()
  { //1
	  int j,k;
	  for (j=0; j<6; j++)
      mins[j]=maxs[j] =v[0].v[j];
      for (k=1; k<numberOfRecords; k++)
         for (j=0; j<6; j++)
         {  //2
              if (v[k].v[j] < mins[j]) mins[j] = v[k].v[j];
              if (v[k].v[j] > maxs[j]) maxs[j] = v[k].v[j];
         } //2
      for (j=0; j<6; j++)
   	     range[j] = maxs[j]-mins[j];
      volume = 1.0;
      for (k = 0; k < 6; k+=2)
		  volume *= range[k]/dims[k].resolution;
      fprintf(monitor,"Cell volume = %lg cubic mm.\n", volume*1.0e9);
      for(k=0; k<numberOfRecords; k++)
          for (j=0; j<6; j++)
             v[k].v[j] = (v[k].v[j] - mins[j]) * dims[j].resolution / range[j];
} //1




