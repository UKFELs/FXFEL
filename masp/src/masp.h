/************************************************
* masp.h. This file holds defined constants,    *
* prototypes and various global variables.      *
************************************************/




/*  constants for the random number generator */

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349
#define PI 3.141592654

// Number of each type of parameter follows

#define STRINGS 5      // Up to 5 string parameters (file names)
#define INTEGERS 1     // One integer parameter (SPREAD)
#define DOUBLES 1      // Up to 1 double parameter (THRESHOLD)
#define DIMENSIONS 6   // Up to six dimension parameters
#define WORD 2         // NONOISE, BINARY
#define DATA 1         // DATA command
#define OUTPUT 1       // OUTPUT command

#define MEMSIZE 2000000000

/**********************
* Constants for erf   *
**********************/

double eparams[10] = {-1.26551223, 1.00002368, 0.37409196,0.09678418,
                      -0.18628806, 0.27886807, -1.13520398,1.48851587,
                      -0.82215223,0.17087277};



int JUMP;               // Increment in expanding v by realloc. It must be
                        // an exact multiple of np





int stage;  // The point reached within each task

int myRank; //  Core allocation figures
int np;     //  Number of cores


int idum = -4563;      // For random number generation
double dummy;

FILE * params;         // for Parameters
FILE * in;             // for data input
FILE * monitor;        // for monitor file
FILE * out[4];         // output files for tasks

int j;

double totalInWeight = 0.0;   // Total weight of the input pulse
double totalOutWeight = 0.0;  // Total weight of the output pulse
double prePoissonWeight = 0.0; // Total weight before applying Poisson rule
long celads[6];  // Expanded Cell address of the current record
double pos[6];    // Position within the cell of the current record

Element * leftTree = NULL;      // Left overflow tree
Element * middleTree = NULL;    // Middle tree
Element * rightTree = NULL;     // Right overflow tree
int treeCount = 0;
int outputRecordCount = 0;

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
int histograms[102][6];   // For data histograms
double * coreBoundaries;  // key coordinate where data in each core starts
//   These variables are used to handle analysis by stripes
int key;                  // The axis used as basis

int aa=0,bb=0,cc=0,dd=0, ee=0, ff=0;   // Performance monitors
double start,finish;      // Timing variables for total run
double startTime,totalTime;  // Timing variables for individual cores
int highest[6];              // Cells with the most entries in each dimension

float pt[100000];          //  Buffer for raw data for noise profile
int pw[100000];
FILE * noiseBuffer;
int tab[51][5] = {{0,0,0,0,0},
              {1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},
              {-1,0,0,0,0},{0,-1,0,0,0},{0,0,-1,0,0},{0,0,0,-1,0},{0,0,0,0,-1},
              {1,1,0,0,0},{1,0,1,0,0},{1,0,0,1,0},{1,0,0,0,1},{0,1,1,0,0},{0,1,0,1,0},{0,1,0,0,1},{0,0,1,1,0},{0,0,1,0,1},{0,0,0,1,1,},
              {1,-1,0,0,0},{1,0,-1,0,0},{1,0,0,-1,0},{1,0,0,0,-1},{0,1,-1,0,0},{0,1,0,-1,0},{0,1,0,0,-1},{0,0,1,-1,0},{0,0,1,0,-1},{0,0,0,1,-1},
              {-1,1,0,0,0},{-1,0,1,0,0},{-1,0,0,1,0},{-1,0,0,0,1},{0,-1,1,0,0},{0,-1,0,1,0},{0,-1,0,0,1},{0,0,-1,1,0},{0,0,-1,0,1},{0,0,0,-1,1},
              {-1,-1,0,0,0},{-1,0,-1,0,0},{-1,0,0,-1,0},{-1,0,0,0,-1},{0,-1,-1,0,0},{0,-1,0,-1,0},{0,-1,0,0,-1},{0,0,-1,-1,0},{0,0,-1,0,-1},{0,0,0,-1,-1}};


//  Variables and default values for parameters follow

char files [STRINGS][40] =
   {"electrons.txt","microparticles.txt","datapicture.csv","weights.txt","bunching.csv"};
int spread = 3;
double threshold = 0.001;
Dimension dims [DIMENSIONS] = {{1,0.0,0,"GAUSSIAN"}, {1,0.0,0,"GAUSSIAN"},
                               {1,0.0,0,"GAUSSIAN"},  {1,0.0,0,"GAUSSIAN"},
                               {1000,0.0,0,"GAUSSIAN"},{100,0.0,0,"GAUSSIAN"}};
int nonoise = 0,binary = 0;    // Set for NONOISE and BINARY
int data[7] = {0,1,2,3,4,5,6};
int dataLen = 7;
int output[7] = {0,1,2,3,4,5,6};
// Parameter names follow

// Parameter names follow This section to be changed if the parameter
// names prove unpopular

char * stringNames[] = {"INPUT","MICROPARTICLES","HISTOGRAM","WEIGHTS","BUNCHING"};
char * intNames[] = {"RANGE"};
char * dimensionNames[] = {"X","PX","Y","PY","Z","PZ"};
char * doubleNames[] = {"THRESHOLD"};
char * dataNames[] = {"DATA"};
char * wordNames[] = {"NONOISE","BINARY"};
char * outputNames[] = {"OUTPUT"};

// Markers to catch double settings
int stringSet[] = {0,0,0,0,0};
int integerSet[] ={0};
int doubleSet[] = {0};
int dimensionSet[] = {0,0,0,0,0,0};
int wordSet[] = {0,0};
int dataSet[] = {0};
int outputSet[] = {0};
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
int inSeamLength,outSeamLength;

// Labels for histogram

char * labels[] = {"X","Px","Y","Py","Z","Pz"};

double tw = 0.0;    // Total Weight of output

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

