/****************************
* maspstructures.h          *
* Various structures        *
* Version 0.7  (21/08/2016) *
****************************/

extern int key;

// A line of input from the electron generator
// Gamma  and various other numbers
// in the input data are omitted. The
// weight is put into the seventh column

typedef struct Record {
	                 double v[7];
                      } Record;


// A true position in six dimensions

typedef struct Position{
	                      double z[6];
		               } Position;


// An output record

typedef struct OutputRecord{
	                            Position p;
	                            double weight;
	             		   }   OutputRecord;


// A cell in internal representation. v is
// six cell coordinates collapsed into a long

typedef struct Cell {
                       long v;
                       float charge;
		            }  Cell;


// What binary trees are made of

typedef struct Element {
	                      Cell cell;       // A cell in internal representation
	                      struct Element * fore;  // Forward and backward
	                      struct Element * back;  // references
		               }  Element;

// Properties of a dimension

typedef struct Dimension{
                          int resolution;
                          double scatter;
                          int scatterIndex;
                          char scatterName[20];
                        } Dimension;


// Properties of a scatter rule

typedef struct ScatterRule{
	                    char scatterName[20];
	                    int size;
	                    double interval;
	                    double * values;
			              } ScatterRule;






