/*********************************************
* masptree.c                                 *
* This module holds (most) of the code       *
* that handles trees                         *
*********************************************/



#ifdef TS

/************************************
* Following module analyses a tree  *
* to determine its statistics       *
************************************/

#define MAXD  (5000)
int fff[MAXD]; int high = 0;
double ttt = 0;   // Counts number of nodes
double sss = 0;   // Counts work
void analyseTree(Element * h, int q)
{
	 if (h->back != NULL)
	    analyseTree(h->back,q+1);
	 if (q > MAXD-1)
	    printf("Core #%d: Tree too deep to analyse\n",myRank);
	 fff[q]++;
     ttt += 1;
     sss += q;
	 if (q > high) high = q;
	 if (h->fore != NULL)
	    analyseTree(h->fore,q+1);
}

void treeDepth(Element * z,char * name)
{
	char ww[100];
	sprintf(ww,"%s%d.csv",name, myRank);
	FILE * f = fopen(ww,"w");
	int j;
	for(j=0; j<MAXD; j++)  fff[j] = 0;
	analyseTree(z,0);
	for (j=0; j<= high; j++)
	  fprintf(f,"  %d, %d %d\n",j,fff[j],j*fff[j]);
    fprintf(f,"Total nodes = %lf, Total work = %lf\n", ttt,sss);
	fclose(f);
}
#endif


/*************************************************
* This function determines if a given vector is  *
* elegible for inclusion in the data set used to *
* assess bunching. A vector is eligible if at    *
* most two of its coordinates differ by 1 from   *
* the 'middle' vector, its other coordinates     *
* being the same. The exact rules are expressed  *
* in array  tab.                                 *
*************************************************/

int elegible (short ww[])
{
     short z[5];
     int j,p=0;
        // Set up a coordinate array that excludes the key coordinate
     for (j=0; j<6; j++)
     if ( j != key)
     {
		 z[p] = ww[j] - highest[j];
		 if(z[p] < (-1) || z[p] > 1)
		    return -1;
		 p++;
     }
//     printf("%d %d %d %d %d accepted\n",z[0],z[1],z[2],z[3],z[4]);
     for (j =0; j<51; j++)
     {
		 if (z[0]== tab[j][0] && z[1]== tab[j][1] && z[2]== tab[j][2] && z[3]== tab[j][3] && z[4]== tab[j][4])
			 return j* dims[key].resolution;
     }
     return -1;
}



Element * makeNode(Cell q)
{
    newElement  = (Element *)(localMalloc(sizeof (Element)));
	if (newElement == NULL)
	{
		char nb[100];
		sprintf(nb,"Core # %d: Not enough memory to build node tree \n",myRank);
		die(nb,99);
    }
    newElement->cell = q;
    newElement->fore = NULL;
    newElement->back = NULL;
    return newElement;
}




/****************************************************
* This function searches and builds a binary tree   *
* of elements. If the item presented is already     *
* there it adds to charge. Otherwise it inserts the *
* element into the selected tree.                   *
****************************************************/

void treeBuild ( Cell  *w, Element ** selectedTree)
{  //1
	long j;
    Position  ww;
    Element * a = *selectedTree;
	while (1)
	{  //2
		if (a == NULL)
		{  //3
			*selectedTree =makeNode(*w);dd++;
			return;
	    }  //3
        else
        j = (w->v-a->cell.v);
        if (j == 0)
        { //3
            a->cell.charge += w->charge;
			cc++;
			return;
	    }  //3
        else
        if ( j < 0)
        {  //3
			  ee++;
			  if (a->back != NULL)
			  {  //4
			      a=a->back;
			      continue;
			  }  //4
			  else
			  {   //4
				  Element* ttt =makeNode(*w);
				  a->back = ttt;

			      bb++;dd++;
			      return;
			  }  //4
	     }  //3
	     else
	     {  //3
			 ff++;
			 if (a->fore != NULL)
			 {  //4
				a=a->fore;
			    continue;
			 }   //4
			 else
			 {   //4
			      a->fore  =makeNode(*w);
			      aa++;dd++;
			      return;
			  } //4
	      }  //3
      }  //2
} //1


