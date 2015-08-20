/***************************************************
*  file maspparams.c                               *
*  This module reads the parameter file for masp,  *
*  fills in defaults, and reports any errors it    *
*  can find. See "Particle conditioner (140302)    *
***************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>


/*************************************************
* This function reads in the details of scatter  *
* descriptions like TOPHAT or BLANCMANGE         *
*************************************************/


	void getScatter(char * name, int dim)
	{
		                  // First, see if it's already there
		FILE * w;
		int j,k,m;
        double * mm;
		for (j = 0; j < rules; j++)
		{
		if (strcmp(name,list[j].scatterName) == 0)
		    {
			dims[dim].scatterIndex = j;
			return;
	        }
	    }        // Otherwise read it in
	    w = fopen(name,"r");
	    if (w == NULL)
	    {
		   error = 1;
		   fprintf(monitor,"Scatter file %s missing\n",name);
		   return;
	    }
	    k = fscanf(w, "%d",&m);
	    if (k != 1)
	    {
			error = 1;
			fprintf (monitor,"In scatter file %s, count of scatter values missing\n",name);
            return;
	    }
	    mm = (double *)(localMalloc((m+1)*sizeof(double)));
	    if (mm == NULL)
	    {
			error = 1;
			fprintf(monitor,"Unlikely error: not enough space for scatter array\n");
			return;
	    }
        j=0;
        while (1)
        {
			k =fscanf(w, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  mm+j,mm+j+1,mm+j+2,mm+j+3,mm+j+4,mm+j+5,mm+j+6,mm+j+7,mm+j+8,mm+j+9);
            if (k<0)
               break;
            j+=k;
            if (j == m)
               break;
	    }
	    for (j=1; j<= m; j++)
			mm[j]+= mm[j-1];

	    for (j=1; j<=m; j++)
	       mm[j] /= mm[m];
	    strcpy(list[rules].scatterName,name);
		list[rules].size = m;
		list[rules].values = mm;
		list[rules].interval = 6.0/m;
		dims[dim].scatterIndex = rules++;
   }

/******************************************************
* Reads all the parameters, and sets up the parameter *
* structure. Reports faults throght the monitor file. *
******************************************************/

int readParameters(FILE * f)
{ //1
	char dump[200];
	char word[30];   // To extract parameter name
	char * m;
	int j,k;
	int lineNumber = 0;
	int cul;
	unsigned int p,q;
	while (1)
	{ //2
		m =fgets(dump, 200, f);
		if (strlen(dump) > 1)
		     fprintf(monitor,"   %s",dump);
		if (m==NULL)
		   break;
        lineNumber++;
        if (empty(m))   // Skip empty lines
           continue;
        if (!(isUpper(dump[0]) || isLower(dump[0])))
           continue;    // Skip comment lines that don't start with a letter

	    j=0;     // Read parameter name and convert to upper case
	    while ((isUpper(dump[j]) || isLower(dump[j])))
	    { //3
		   if (isUpper(dump[j]))
		    word[j] = dump[j];
		   else
		   if (isLower(dump[j]))
		      word[j]= (char)(toUpper(dump[j]));
		   j++;

	    }  //3
	    word[j] = 0;  // Plant terminator;
	        // Handle integer parameters
        for (k = 0; k < INTEGERS; k++)
	    if (strcmp(word, intNames[k])==0)
	    {  //  3
	       if(integerSet[k])
	       {  //4
	    	     dump[strlen(dump)-1] = 0;
                 fprintf(monitor,"Line %2d: Error in line \"%s\":- %s set twice\n",lineNumber,dump,word);
	    	     error = 1;
	       }  //4
	       else
	      {  //4
	    	p= sscanf(dump+j,"%ld ",&spread);
		    if (p!=1)
		    {  //5
		        dump[strlen(dump)-1] = 0;
                fprintf (monitor,"Line %2d: Error in line \"%s\":- Value for %s missing\n",lineNumber,dump, word);
		        error = 1;
		    }  //5

		   integerSet[k] = 1;
          }  //4
	      break;  // Because we have recognised an integer
	            // and don't need to go enny further
	    }  //3

	    if (k < INTEGERS)
	       continue;
	    // Dun handling integers. Now look at doubles
	    // Handle double parameters
        for (k = 0; k < DOUBLES; k++)
	    if (strcmp(word, doubleNames[k])==0)
	    {  //  3
	       if(doubleSet[k])
	       {  //4
	    	     dump[strlen(dump)-1] = 0;
                fprintf(monitor,"Line %2d: Error in line \"%s\":- %s set twice\n",lineNumber,dump,word);
	    	    error = 1;
	       }  //4
	       else
	       {  //4
		       p= sscanf(dump+j,"%lf ",&threshold);
		       if (p!=1)
		       {  //5
		          dump[strlen(dump)-1] = 0;
		          fprintf (monitor,"Line %2d: Error in line \"%s\":- Value for %s missing or incorrect\n",lineNumber, dump, word);
		          error = 1;
		       }  //5

		       doubleSet[k] = 1;
          }  //4
	      break;  // Because we have recognised a double
	              // and don't need to go enny further
	   }  //3

	   if (k < DOUBLES)
	      continue;
	      // Dun handling doubles. Now look at strings

       for (k = 0; k < STRINGS; k++)
	   if (strcmp(word, stringNames[k])==0)
	   {  //  3
	      if(stringSet[k])
	      { //4
			dump[strlen(dump)-1] = 0;
            fprintf(monitor,"Line %2d: Error in line \"%s\":- %s set twice\n",lineNumber,dump,word);
		    error = 1;
	      } //4
	      else
	      { //4
	        char ddd[40];
		    p= sscanf(dump+j,"%s",ddd);
		    if (p!=1)
		    {  //5
		        dump[strlen(dump)-1] = 0;
                fprintf (monitor,"Line %2d: Error in line \"%s\":- Value for %s missing\n",lineNumber,dump, word);
		        error = 1;
		    }  //5
		    if (strlen(ddd) > 35)
		    {  //5
		       dump[strlen(dump)-1] = 0;
               fprintf (monitor,"Line %2d: Error in line \"%s\":- \n    %s\n    is too long for a file name\n",lineNumber,dump,ddd);
		       error = 1;
		    }  //5
            strcpy(files[k],ddd);
		    stringSet[k] = 1;
	     } //4
	     break;
	   }  //3
       if (k < STRINGS)
          continue;
               // Dun handling strings. Now look at NONOISE and BINARY

	   for(k=0; k<WORD; k++)
	   if (strcmp(word,wordNames[k])==0)
	   {  //3
		 if (wordSet[k] )
		 {  //4
			 dump[strlen(dump)-1] = 0;
             fprintf(monitor,"Line %2d: Error in line \"%s\":- %s set twice\n",lineNumber,dump,word);
			 error = 1;
		 } //4
		 else
		 {   //4
			 if (k==0)
			    nonoise = 1;
			 else
			    binary = 1;
			 wordSet[k] = 1;
		 } //4
         break;
      } //3
      if (k< WORD)
         continue;
           // Dun handling NONOISE. Now look at DATA

	  for(k=0; k<DATA; k++)
	  if (strcmp(word,dataNames[k])==0)
	  {  //3
		 if (dataSet[k] )
		 {  //4
			 dump[strlen(dump)-1] = 0;
             fprintf(monitor,"Line %2d: Error in line \"%s\":- DATA set twice\n",lineNumber,dump);
			 error = 1;
		 } //4
		 else
		 {   //4
			 // Allow for more values than are required
			 p = sscanf(dump+j, "%d %d %d %d %d %d %d,%d,%d,%d,%d,%d,%d,%d",
			            data,data+1,data+2,data+3,data+4,data+5,data+6,&cul,&cul,&cul,&cul,&cul,&cul,&cul);
			 if (p != 7)
			 {  //5
				 dump[strlen(dump)-1] = 0;
                 fprintf(monitor, "Line %2d: Error in line \"%s\":- Wrong number of values in data command\n",lineNumber,dump);
				 error = 1;
		     }  //5
		     if(p == 7)  // Only do check if the number of values is correct
		     {  //5
				 for (p=0; p<6; p++)
		            for (q = p+1; q<7;q++)
		              if (data[p] == data[q])
		              { //6
					    dump[strlen(dump)-1] = 0;
                        fprintf(monitor, "Line %2d: Error in line \"%s\":-\n Column numbers in DATA command not unique\n",lineNumber,dump);
					    error = 1;
			           } //6
		     }//5
			 dataSet[k] = 1;
		 } //4
         break;
      } //3

     if ( k < DATA )
         continue;
         // Dun dealing with DATA Now look at OUTPUT

	 for(k=0; k<OUTPUT; k++)
	 if (strcmp(word,outputNames[k])==0)
	 {  //3
		 if (outputSet[k] )
		 {  //4
			 dump[strlen(dump)-1] = 0;
             fprintf(monitor,"Line %2d: Error in line \"%s\":- OUTPUT set twice\n",lineNumber,dump);
			 error = 1;
		 } //4
		 else
		 {   //4
			 // Allow for more values than are required
			 p = sscanf(dump+j, "%d %d %d %d %d %d %d,%d,%d,%d,%d,%d,%d",
			            output,output+1,output+2,output+3,output+4,output+5,output+6,&cul,&cul,&cul,&cul,&cul,&cul,&cul);
			 if (p != 7)
			 {  //5
				 dump[strlen(dump)-1] = 0;
                 fprintf(monitor, "Line %2d: Error in line \"%s\":- Wrong number of values in OUTPUT command\n",lineNumber,dump);
				 error = 1;
		     }  //5
		     if(p == 7)  // Only do check if the number of values is correct
		     {  //5
				 for (p=0; p<6; p++)
		            for (q = p+1; q<7;q++)
		              if (output[p] == output[q])
		              { //6
					    dump[strlen(dump)-1] = 0;
                        fprintf(monitor, "Line %2d: Error in line \"%s\":-\n Column numbers in OUTPUT command not unique\n",lineNumber,dump);
					    error = 1;
			           } //6
		     }//5
			 for (p=0; p<7; p++)
			 if ( output[p] < 0 || output[p] > 6)
		     { //5
					dump[strlen(dump)-1] = 0;
                    fprintf(monitor, "Line %2d: Error in line \"%s\":-\n Column numbers in OUTPUT command not in the range 0-6\n",lineNumber,dump);
					    error = 1;
			 } //5

			 outputSet[k] = 1;
		 } //4
         break;
      } //3
      if ( k < OUTPUT )
         continue;
         // Dun dealing with Output Now look at Dimensions

      for (k=0; k< DIMENSIONS; k++)
      if(strcmp(word,dimensionNames[k]) == 0)
      {   //3
           if (dimensionSet[k])
           { //4
		         dump[strlen(dump)-1] = 0;
                 fprintf(monitor,"Line %2d: Dimension \"%s\":%s set twice\n",lineNumber,word);
		         error = 1;
	       } //4
	       else
	       {  //4
	           char ddd[40];
	           int res;
	           double scatter;
	           p = sscanf(dump+j, "%d %lf %s", &res,&scatter,ddd);
		       if (p!=3)
		       {  //5
		          dump[strlen(dump)-1] = 0;
                  fprintf (monitor,"Line %2d: Values for axis %s are missing\n", lineNumber,word);
		          error = 1;
		       }  //5
		       if(res == 1 && scatter >0)
		       {  //5
				   dump[strlen(dump)-1] = 0;
				   fprintf(monitor,
				   "Line %2d:In axis %s the only sensible scatter for one cell is zero\n", lineNumber,word);
				  error = 1;
		       }  //5
		       if (res <= 0)
		       { //5
		            dump[strlen(dump)-1] = 0;
                    fprintf( monitor,"Line %2d: %d is not an acceptable number of cells for %s\n",lineNumber,res,word);
	                error = 1;
	           }  //5
	           if (strcmp(ddd, "GAUSSIAN") != 0)
	           {  //5
				  FILE * testFile = fopen(ddd,"r");
				  if (testFile == NULL)
				  { //6
					 dump[strlen(dump)-1] = 0;
					 fprintf( monitor,"Line %2d: Cannot find file %s\n",lineNumber,ddd);
				     error = 1;
			      }  //6
			      else
			       fclose(testFile);
		      } //5

			 for (p=0; p<strlen(ddd); p++)
	         ddd[p] = toUpper(ddd[p]);
	         dims[k].resolution = res;
	         dims[k].scatter = scatter;
	         strcpy(dims[k].scatterName,ddd);
	         dimensionSet[k] = 1;
	       } // 4
           break;
	   } //3
       if (k == DIMENSIONS)
	   {  //3
		    dump[strlen(dump)-1] = 0;
            fprintf (monitor,"Line %2d: Error in line \"%s\":-This is a bad parameter line: %s \n",lineNumber,dump);
			error = 1;
	   } //3
     } //2
     strcpy(list[0].scatterName, "GAUSSIAN");
     list[0].size = 0;
     list[0].values = NULL;
     for (j=0; j<6; j++)
		getScatter(dims[j].scatterName, j);
         // some global error checks
     if (nonoise && (stringSet[4] == 1))
	 {  //2
         fprintf (monitor,"If NONOISE is set you can't ask for BUNCHING\n");
		 error = 1;
	 } //2
     if ( !nonoise && (stringSet[1] == 0) && (stringSet[4] == 1))
	 {  //2
         fprintf (monitor," If MICROPARTICLES is not set you can't ask for BUNCHING\n",lineNumber);
		 error = 1;
	 } //2
     for (j = 0; j<3; j++)
        for (k = j+1; k<4; k++)
          if (strcmp(files[j], files[k]) == 0)
          {  //3
				fprintf(monitor, "File names for %s and %s clash\n",stringNames[j],stringNames[k]);
				error = 1;
		  } //3
 return error;
} //1

/********************************************
* Sends all parameters (including defaults) *
* to the monitor file.                      *
********************************************/

void displayParameters()
{
	int j;
	fprintf(monitor,"\n\nInteger parameters:");
	for (j = 0; j< INTEGERS; j++)
	   fprintf(monitor,"  %6s=%-6d\n", intNames[j], spread);
	fprintf(monitor,"Double parameters: ");
	for (j = 0; j< DOUBLES; j++)
	   fprintf(monitor,"  %6s=%lf\n", doubleNames[j], threshold);
	fprintf(monitor,"\nString parameters:\n");
		for (j = 0; j< STRINGS; j++)
		if(j == 0 || stringSet[j])
		    fprintf(monitor,"  %14s = %s\n", stringNames[j], files[j]);

	fprintf(monitor,"\nDimension parameters:\n\n");
	for (j = 0; j< DIMENSIONS; j++)
			fprintf(monitor," %2s: Res =  %6d scatter =  %lf  Name = %12s \n",
	           dimensionNames[j],dims[j].resolution,dims[j].scatter, dims[j].scatterName);
    fprintf(monitor,"\nDATA: %d %d %d %d %d %d %d\n", data[0],data[1],data[2],data[3],data[4],data[5],data[6]);

	fprintf(monitor,"\n");
	fprintf(monitor,"Scatter rules\n");
	for (j=0; j<rules; j++)
		fprintf(monitor,"Name = %12s\n",list[j].scatterName);
	if (nonoise)
	    fprintf(monitor, "NONOISE\n");
	else
	    fprintf (monitor, "Noise adjustments included\n");
	if (binary)
	    fprintf(monitor, "BINARY\n");
	else
	    fprintf(monitor, "Output in ASCII (decimal)\n");

	fprintf(monitor, "OUTPUT order is %d %d %d %d %d %d %d\n", output[0],output[1],output[2],output[3],output[4],output[5],output[6]);
}

