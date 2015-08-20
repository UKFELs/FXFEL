/***************************************
* masphisto.c                          *
This function works out and displays *
* the distribution histograms          *
***************************************/



void doHistogram()
{{  //1
   if (myRank == 0)
   {  //2
	  int j,k;
	  char top[100];
	  *top = 0;
	  for(j=0; j<6; j++)
	  {  //3      // Assemble top label
		  strcat(top, labels[output[j]]);
		  if (j < 5)
		      strcat(top,",");
		  else
		      strcat(top,"\n");
          }  //3

	  for(j = 0; j<101; j++)  // Clear histograms
       for(k=0; k<6; k++)
          histograms[j][k] = 0;

       for(j=0; j<numberOfRecords; j++)
	   for (k = 0; k< 6; k++)
	   { //3            // Assemble the data
		double q = v[j].v[k]/dims[k].resolution;
		int s = (int)(100*q+0.5);
	         histograms[s][k] ++;
	   }  //3
           fprintf(out[1],top);
	   for(j=0; j<101; j++)
                fprintf(out[1],"%d,%d,%d,%d,%d,%d\n",
                    histograms[j][output[0]],histograms[j][output[1]],histograms[j][output[2]],
                    histograms[j][output[3]],histograms[j][output[4]],histograms[j][output[5]]);
       fclose(out[1]);
    } //2
}}  //1


