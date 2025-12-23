#define DEFAULT		0
#define ADDITION	1

#define NoDirection   0
#define RIGHT  1
#define LEFT   2
#define BOTH   3

#define NoLaser  0
#define Gaussian  1
#define SSTF   2

#define Boundary  1
#define Shot   2

typedef struct _LaserList  {
   int polarity;
   int mode,loadMethod;
   double lambda;
   double omega;
   double amplitude;   //unit is a0.
   double rU;
   double rD;
   double retard;
   double flat;
   int loadPointX; 
   int loadPointY; 
   int loadPointZ; 

   double rayleighLength,beamWaist,elliptic;
   double focus;
   int direction;
   int add;
   double gdd;

   // SSTF
   double lensFocus,sigOmegaSSTF,sigX_t;
   double alphaYSSTF,sigYSSTF,sigY_t;
   double alphaZSSTF,sigZSSTF,sigZ_t;


   struct _LaserList *next;
} LaserList;
