/***************************************
* masptree.c                           *
* This module holds (most) of the code *
* that handles trees                   *
* Version 0.7  (21/08/2016)            *
***************************************/

Element * makeNode(Cell q)
{  //1
    newElement  = (Element *)(localMalloc(sizeof (Element)));
	if (newElement == NULL)
	{  //2
		char nb[100];
		sprintf(nb,"Core # %d: Not enough memory to build node tree \n",myRank);
		die(nb,99);
    }  //2
    newElement->cell = q;
    newElement->fore = NULL;
    newElement->back = NULL;
    return newElement;
}  //2


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
			*selectedTree =makeNode(*w);
			return;
	    }  //3
        else
        j = (w->v-a->cell.v);
        if (j == 0)
        { //3
            a->cell.charge += w->charge;
			return;
	    }  //3
        else
        if ( j < 0)
        {  //3
			  if (a->back != NULL)
			  {  //4
			      a=a->back;
			      continue;
			  }  //4
			  else
			  {   //4
				  Element* ttt =makeNode(*w);
				  a->back = ttt;

			      return;
			  }  //4
	     }  //3
	     else
	     {  //3
			 if (a->fore != NULL)
			 {  //4
				a=a->fore;
			    continue;
			 }   //4
			 else
			 {   //4
			      a->fore  =makeNode(*w);
			      return;
			  } //4
	      }  //3
      }  //2
} //1


