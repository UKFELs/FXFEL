
void treeBuild ( Cell  *w, Element ** selectedTree);
void setPlot(double positions[]);



void presentWeight(Cell q)
{{
	short int g[6];
	double X,Y,Z;
	Binary3 bb;
	expand (q.v, g);
	tw += q.charge;
	X = mins[0]+(g[0]+0.5)*range[0]/dims[0].resolution;
	Y = mins[2]+(g[2]+0.5)*range[2]/dims[2].resolution;
	Z = mins[4]+(g[4]+0.5)*range[4]/dims[4].resolution;
	if (binary)
	{
		bb.v[0] = (float)X;
		bb.v[1] = (float)Y;
		bb.v[2] = (float)Z;
		bb.v[3] = (float) q.charge;
		fwrite((void *)(&bb),sizeof(Binary3),1,out[2]);
    }
    else
	fprintf(out[2], "%10.5lg %10.5lg %10.5lg %lf\n",X,Y,Z, q.charge);
}}

void exploreWeight(Element * h)
{{

	 if (h->back != NULL)
		exploreWeight (h->back);
		presentWeight(h->cell);
	 if (h->fore != NULL)
		exploreWeight (h->fore);
}}


/******************************************************
*  This function is used to process each record for   *
*  the WEIGHT parameter. Only the spatial coordinates *
*  are used, so the tree will be small enough to live *
*  entirely in core 0. No stripes needed!             *
******************************************************/

void processWeightRecord(int j)
{{ //1
    int k;
    double dweight;
    double totwt = 0;
    int d[6];
    Record m = v[j];
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
             // Only three dimension!
    for (d[0] =0; d[0]<plotSizes[0]; d[0]++)
           for (d[2] =0; d[2]<plotSizes[2]; d[2]++)
                for (d[4] =0; d[4]<plotSizes[4]; d[4]++)
               { //2
		           dweight = (m.v[6] *plot[4][d[4]]*plot[2][d[2]]*plot[0][d[0]]);
                   if (dweight > threshold)
                  { //3
                     cc.charge = dweight;
                     for (k=0; k<6; k++)
			         ww[k] = (k&1)? 0:celads[k]-plotMidPoints[k]+d[k];  // Miss out momentum coordinates
                     cc.v = map(ww);
		             treeBuild(&cc, &middleTree);
		          } //3
	          } //2
}} //1



/****************************************
* This function finds the total charge  *
* in each physical block of the output  *
* I guess we don't need to use stripes  *
****************************************/

void doWeights()
{{
    int j;
    printf("Doing Weights\n");
    middleTree = NULL;
    for (j=0; j<numberOfRecords; j++)
	processWeightRecord(j);
    exploreWeight(middleTree);
    printf("Total weight found in WEIGHT option = %lf\n", tw);
    fclose(out[2]);
}}
