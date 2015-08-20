/**********************
* maspstructures.h    *
* Various structures  *
**********************/

typedef struct Record {           // A line of input from the electron generator
	                 double v[7]; // Gamma  and arious other numbers
	                              //  in the input data are omitted. The
                      } Record;   // weight is put into the seventh column
extern int key;


typedef struct Position{
	                  double z[6];    // A true position in six dimensions
		       } Position;


typedef struct OutputRecord{    // An output record
	                     Position p;
	                     double weight;
			   }   OutputRecord;


typedef struct DensityRecord{   // A record of electron density
                                // for a given cell
                              double x;
                              double y;
                              double z;
                              double charge;
                              double density;
			    } DensityRecord;


typedef struct Cell {                 // A cell in internal representation
                       long v;        // six cell coordinates collapsed into a long
                       float charge;
		    }  Cell;


typedef struct Element {       // What binary trees are made of
	                      Cell cell;       // A cell in internal representation
	                      struct Element * fore;  // Forward and backward
	                      struct Element * back;  // references
		       }  Element;

typedef struct Dimension{       // Properties of a dimension
                          int resolution;
                          double scatter;
                          int scatterIndex;
                          char scatterName[20];
                        } Dimension;


typedef struct ScatterRule{      // Properties of a scatter rule
	                    char scatterName[20];
	                    int size;
	                    double interval;
	                    double * values;
			  } ScatterRule;


typedef struct Binary3{    // Used for binary output of weights
                          float v[4];
				      }  Binary3;




