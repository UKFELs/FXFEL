/*********************************************
* masp.h. This file holds defined constants, *
* prototypes and various global variables.   *
* Version 0.7  (21/08/2016)                  *
*********************************************/


// Number of each type of parameter follows

#define STRINGS 2      // Up to 2 string parameters (file names)
                       // One for INPUT, and the other for output
#define DOUBLES 1      // Up to 1 double parameter (THRESHOLD)
#define DIMENSIONS 6   // Up to six dimension parameters
#define INTEGERS 2     // One for Range, one for Decimals
#define MEMSIZE 2000000000
#define MAX_CORES  100

int JUMP;               // Increment in expanding v by realloc. It must be
                        // an exact multiple of np

/***********************************
* Details for multi-process code   *
***********************************/

int myRank; //  Core allocation figures
int np;     //  Number of cores


int idum = -4563;      // For random number generation
double dummy;          // Use in main to set up random number generator

FILE * params;         // for Parameters
FILE * in;             // for data input
FILE * monitor;        // for monitor file
FILE * out;            // output files for tasks

int j;
int recordLength;      // Length of the standard output record
char format[200];      // Format statement for standard output record
double totalInWeight = 0.0;   // Total weight of the input pulse
double totalOutWeight = 0.0;  // Total weight of the output pulse

long celads[6];  // Expanded Cell address of the current record
double pos[6];    // Position within the cell of the current record

Element * leftTree = NULL;      // Left overflow tree
Element * middleTree = NULL;    // Middle tree
Element * rightTree = NULL;     // Right overflow tree

int error = 0;  // Set to 1 if any error is found.
                // The program can't go on.

short plotSizes[6];   // Used to split up a macroparticle
short plotMidPoints[6];
double * plot [6];


// These variables used to scale the input records to
// fit the grid

double mins[6],maxs[6],range[6],multiplier[6];

double volume;  // Cell volume


Record * v;               // Array to hold all input records in core 0
Record * localv;          // Array to hold records in each core
int numberOfRecords = 0;  // Actual number of records read so far
int vSize;                // Size of input record buffer
int localNumberOfRecords; // Number of records ent to each process
double * coreBoundaries;  // key coordinate where data in each core starts
//   These variables are used to handle analysis by stripes
int key;                  // The axis used as basis

double start,finish;      // Timing variables for total run
double startTime,totalTime;  // Timing variables for individual cores



//  Variables and default values for parameters follow

char files [STRINGS][40] =
   {"electrons.txt","microparticles.txt"};
char microfile[50];    // Name of file for microparticle output (used by all cores)
MPI_File mf;             // Handle for microfile
int spread = 3;
int Range = 3;          // Capital letter is deliberate
int decimals = 6;       // Number of decimals in the output;
double threshold = 0.005;
Dimension dims [DIMENSIONS] = {{20,1.0,0,"GAUSSIAN"}, {20,1.0,0,"GAUSSIAN"},
                               {20,1.0,0,"GAUSSIAN"},  {20,1.0,0,"GAUSSIAN"},
                               {1000,1.0,0,"GAUSSIAN"},{100,1.0,0,"GAUSSIAN"}};
int nonoise = 0;    // Set for NONOISE and BINARY

// Parameter names follow This section to be changed if the parameter
// names prove unpopular


char * stringNames[] = {"INPUT","MICROPARTICLES"};
char * intNames[] = {"RANGE","DECIMALS"};
char * dimensionNames[] = {"X","PX","Y","PY","Z","PZ"};
char * doubleNames[] = {"THRESHOLD"};
char * wordNames[] = {"NONOISE"};

// Markers to catch double settings
int stringSet[] = {0,0};
int integerSet[] ={0,0};
int doubleSet[] = {0};
int dimensionSet[] = {0,0,0,0,0,0};
int wordSet[] = {0};
extern FILE * monitor;

           // Rules for alternative scattering

ScatterRule list[DIMENSIONS];
int rules = 1;  // Number of actual scatter rules

FILE * za;    //  Variables for maspparticles
FILE * zb;
FILE * zc;
int lowerBoundary;
int upperBoundary;
Element * newElement;
int dun = 0;
int left=0, middle = 0, right = 0;
Cell * inSeam = NULL;
Cell * outSeam = NULL;
long inSeamLength,outSeamLength;

double tw = 0.0;    // Total Weight of output
MPI_Offset root;
long ww[MAX_CORES];   // Starting output addresses for all the cores






/****************
*  Prototypes   *
****************/

int isUpper(char w);
int isLower (char w);
int toUpper (char w);
int empty(char * w);
double grandom(int* idum);
double gauss(double mu, double sigma);
int poisson (double lambda);
double erf(double x);
void die(char * s, int e);

