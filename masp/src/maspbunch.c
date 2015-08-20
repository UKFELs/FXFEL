/*********************************************
* maspbunch.c                                *
* All the data is prepared ennyway. All this *
* has to do is the analysis.                 *
*********************************************/

#define PERIOD 10

void doBunching()
{{  //1
   printf("Starting doBunching\n");
   int j,n,cycles,start,x,m;
   int hsl[150], bsl[150];
   double b [10000];
   double br,bi,modb,ex,eb2=0.0,ebs=0.0;
   double dtheta = 2*PI/PERIOD;
            // Zeroise buffers to accumulate results
   for(j=0; j<150; j++)
        hsl[j]=bsl[j] = 0;
   n= 51*dims[4].resolution;  // Total number of readings
   cycles = n/PERIOD;
   printf("Cycles = %d\n",cycles);
   for(j=0; j<cycles;j++)
   {  //2
      b[j] = 0;
      bi = br = 0.0;
      ex=0; x=0;
      start = j*PERIOD;   // Where this period starts
      for(m = 0; m<PERIOD; m++)
      {  //3
         int w = start + m;
         double noise = pt[w]-floor(pt[w]) - 0.5;
         if (pw[w] > 0) x++;    // Count non-zero elements
         ex += pw[w];
         br += pw[w] * cos(m*dtheta+noise);
         bi += pw[w] * sin(m*dtheta+noise);
      } //3
      if (x<3) continue;   // Ignore periods with < 3 valid points
      ex/=PERIOD;
      br *= PERIOD/ex;
      bi *= PERIOD/ex;
      modb = sqrt(br*br+bi*bi);
      b[j] = modb;
      eb2+= modb;
      ebs += modb*modb;
   }  //  2

   eb2 /=(cycles-1);
   ebs /=(cycles-1);
   int maxValue = 0;
   for ( j = 0; j<cycles; j++)
   {   // 2
      if (b[j] > 0)
      {   //3
         int w = (int) (10*b[j]/(eb2+0.5));
         if (w < 150)
             hsl[w]++;
         if(w > maxValue) maxValue = w;
         w = (int)(10*b[j]*b[j]/ebs);
         if (w<150)
            bsl[w]++;
         if(w > maxValue) maxValue = w;
      } //3
   } // 2
   fprintf(out[3],"Rayleigh,Negative-exponential\n");
   for (j=0; j<=maxValue; j++)
  	  fprintf(out[3], "%d,%d\n",hsl[j],bsl[j]);
   fclose(out[3]);
}}

