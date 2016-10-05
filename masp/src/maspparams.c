/****************************************************
*  file maspparams.c                                *
*  This module reads the parameter file for masp,   *
*  fills in defaults, and reports any errors it     *
*  can find. See "Particle conditioner (140302)     *
*  Late amandment: Now reads parameters in NameList *
*  format. Any change to parameters must be coded   *
*  into this module.                                *
*  Version 0.7  (21/08/2016)                        *
****************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define stricmp strcasecmp

/*************************************************
* This function reads in the details of scatter  *
* descriptions like TOPHAT or BLANCMANGE         *
*************************************************/


void getScatter(char * name, int dim)
{  //1
					  // First, see if it's already there
	FILE * w;
	int j,k,m;
	double * mm;
	for (j = 0; j < rules; j++)
	{  //2
	    if (strcmp(name,list[j].scatterName) == 0)
		  {  //3
		  dims[dim].scatterIndex = j;
		  return;
		  }  //3
	}//2        // Otherwise read it in
	w = fopen(name,"r");
	if (w == NULL)
	{  //2
	   error = 1;
	   fprintf(monitor,"Scatter file %s missing\n",name);
	   return;
	} //2
	k = fscanf(w, "%d",&m);
	if (k != 1)
	{  //2
		error = 1;
		fprintf (monitor,"In scatter file %s, count of scatter values missing\n",name);
		return;
	} //2
	mm = (double *)(localMalloc((m+1)*sizeof(double)));
	if (mm == NULL)
	{  //2
		error = 1;
		fprintf(monitor,"Unlikely error: not enough space for scatter array\n");
		return;
	} //2
	j=0;
	while (1)
	{  //2
		k =fscanf(w, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			  mm+j,mm+j+1,mm+j+2,mm+j+3,mm+j+4,mm+j+5,mm+j+6,mm+j+7,mm+j+8,mm+j+9);
		if (k<0)
		   break;
		j+=k;
		if (j == m)
		   break;
	} //2
	for (j=1; j<= m; j++)
		mm[j]+= mm[j-1];

	for (j=1; j<=m; j++)
	   mm[j] /= mm[m];
	strcpy(list[rules].scatterName,name);
	list[rules].size = m;
	list[rules].values = mm;
	list[rules].interval = 6.0/m;
	dims[dim].scatterIndex = rules++;
} //1

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
        if( strstr(dump,"&MASP")!= NULL)
        {
			int error = readNameList(f);
            return error;
	    }
	}  //2
} //1


/**************************************
*  This function is used where a      *
*  boolean value is expected. The     *
*  input is a string. Any string that *
*  starts with "t"  of ".t" counts    *
*  as true, and any other string      *
*  counts as false. There is no       *
*  error return.                      *
**************************************/

int truth (char * x)
{ //1
	return (   strncasecmp(x,"t",1) == 0
	    || strncasecmp(x,".t",2) == 0
	   );
} //1

/*******************************
*  This function checks that a *
*  scatter file name is OK. It *
*  must be either GAUSSIAN or  *
*  the name of a file that can *
*  be opened                   *
********************************/

int checkName(char * v)
{  //1
	if (strcmp(v,"GAUSSIAN")==0)
	  return 1;
	FILE * w;
	w = fopen(v,"r");
	if (w == NULL)
	{ //2
		error = 1;
		return 0;
	} //2
	fclose(w);
	return 1;
} //1


/*******************************************************
*  This function (added 160213) reads the parameters   *
*  in NameList format. It is assumed that a header     *
*  in the form &MASP has already been read.            *
*******************************************************/

int readNameList(FILE * f)
{  //1
	char s[5000] ;
	char dump[200];
	char clump[150];
	char name[40];
	char value[500];
	int p = 0; // Pointer to s
	int j; int k;
	error = 0;
	*s = 0;   // S contains only a terminator
	char * m;
	         // This section reads in the namelist and arranges it
	         // in one long string, with each item terminated by a comma.
	printf("called readNameList\n");
        while (1)
	{  //2
		m =fgets(dump, 200, f);
		if (dump[0] == '|')
		   continue;             // Ignore comments
                if(dump[0] == '/')
                break;
	    j=k=0;
        while (j < strlen(dump)-1)
	    {  //3
			char w = dump[j++];
			if (w == '"' || w == '\'')
			{  //4
				do
				{  //5
					clump[k++] = dump[j++];
			    }   //5
			    while (clump[k-1] != w);
			    k--;
		     }  //4
		     else
		     if (w != ' ')
		       clump[k++] = w;
	    }  //3
	    clump[k++] = ',';
	    clump[k] = 0;
            strcat ( s,clump);
    }//2
    while (p < strlen(s))
    {//2
        int q = 0;
        while (s[p++] !='=')
          name[q++]= s[p-1];
		  if (q > 30)
		  {  //6
			name[k] = 0;
			fprintf(monitor,"Parameter name too long. Starts with %s\n",name);
			error = 1;
			return error;
		  }  //6
        name[q] = '\0';
        q=0;
        while (s[p++] != ',')
        { //3
			if (q>140)
			{  //4
			   value[q] = 0;
			   fprintf(monitor,"Value too long : %s\n",value);
			   error = 1;
			   return error;
			}  //4
			value[q++] = s[p-1];
	    }  //3

        value[q] = '\0';

        if(stricmp(name,"nonoise")== 0)   // Handle Boolean parameters
           nonoise = truth(value);

        else if (stricmp(name,"Input")== 0)    // Handle file names
          strcpy(files[0],value);
        else if (stricmp(name,"Output") == 0)
          strcpy(files[1],value);

        else if(stricmp(name,"X_resolution") == 0)    // Handle integer parameters
             sscanf(value,"%d", &(dims[0].resolution));
        else if(stricmp(name,"Y_resolution") == 0)
             sscanf(value,"%d", &(dims[2].resolution));
        else if(stricmp(name,"Z_resolution") == 0)
             sscanf(value,"%d", &(dims[4].resolution));
        else if(stricmp(name,"Px_resolution") == 0)
             sscanf(value,"%d", &(dims[1].resolution));
        else if(stricmp(name,"Py_resolution") == 0)
             sscanf(value,"%d", &(dims[3].resolution));
        else if(stricmp(name,"Pz_resolution") == 0)
             sscanf(value,"%d", &(dims[5].resolution));
        else if(stricmp(name,"Range") == 0)
             sscanf(value,"%d",&Range);
        else if(stricmp(name,"Decimals") == 0)
             sscanf(value,"%d",&decimals);


        else if(stricmp(name,"X_width")==0)     // Handle real parameters
             sscanf(value,"%lf",&(dims[0].scatter));
        else if(stricmp(name,"Y_width")==0)
             sscanf(value,"%lf",&(dims[2].scatter));
        else if(stricmp(name,"Z_width")==0)
             sscanf(value,"%lf",&(dims[4].scatter));
        else if(stricmp(name,"Px_width")==0)
             sscanf(value,"%lf",&(dims[1].scatter));
        else if(stricmp(name,"Py_width")==0)
             sscanf(value,"%lf",&(dims[3].scatter));
        else if(stricmp(name,"Pz_width")==0)
             sscanf(value,"%lf",&(dims[5].scatter));
        else if(stricmp(name,"Threshold")==0)
             sscanf(value,"%lf",&threshold);

        else if(stricmp(name,"X_scatter") == 0)  // Handle scatter file names
        {  //3
             if(checkName(value))
                strcpy(dims[0].scatterName,value);
             else
                fprintf(monitor,"Scatter name not recognised: %s\n",value);
	    }  //3
        else if(stricmp(name,"Y_scatter") == 0)
        {  //3
             if(checkName(value))
                strcpy(dims[2].scatterName,value);
             else
                fprintf(monitor,"Scatter name not recognised: %s\n",value);
	    }  //3
        else if(stricmp(name,"Z_scatter") == 0)
        {  //3
             if(checkName(value))
                strcpy(dims[4].scatterName,value);
             else
                fprintf(monitor,"Scatter name not recognised: %s\n",value);
	    }  //3
        else if(stricmp(name,"Px_scatter") == 0)
        {  //3
             if(checkName(value))
                strcpy(dims[1].scatterName,value);
             else
                fprintf(monitor,"Scatter name not recognised: %s\n",value);
	    }  //3
        else if(stricmp(name,"Py_scatter") == 0)
        {  //3
             if(checkName(value))
                strcpy(dims[3].scatterName,value);
             else
                fprintf(monitor,"Scatter name not recognised: %s\n",value);
	    }  //3
        else if(stricmp(name,"Pz_scatter") == 0)
        {  //3
             if(checkName(value))
                strcpy(dims[5].scatterName,value);
             else
                fprintf(monitor,"Scatter name not recognised: %s\n",value);
	    }  //3
        else
        {  //3
			fprintf(monitor,"Parameter name not recognised: %s\n",name);
			error = 1;
	    }  //3
    }//2
    strcpy(list[0].scatterName, "GAUSSIAN");
    list[0].size = 0;
    list[0].values = NULL;
    for (j=0; j<6; j++)
	   getScatter(dims[j].scatterName, j);
    return error;

} // 1




/********************************************
* Sends all parameters (including defaults) *
* to the monitor file.                      *
********************************************/

void displayParameters()
{  //1
	int j;
	fprintf(monitor,"\n\nInteger parameters:");
	   fprintf(monitor,"  %8s=%-6d  ", intNames[0], Range);
	   fprintf(monitor,"  %8s=%-6d\n", intNames[1], decimals);
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

	fprintf(monitor,"Scatter rules\n");
	for (j=0; j<rules; j++)
		fprintf(monitor,"Name = %12s\n",list[j].scatterName);
	if (nonoise)
	    fprintf(monitor, "NONOISE\n");
	else
	    fprintf (monitor, "Noise adjustments included\n");
}  //1

