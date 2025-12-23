#define Electron 	1
#define Positron 	10
#define Test	 	2
#define HPlus0 	 	100
#define HPlus1 	 	101
#define HePlus0 	200
#define HePlus1 	201
#define HePlus2 	202
#define CPlus0          600
#define CPlus1          601
#define CPlus2          602
#define CPlus3          603
#define CPlus4          604
#define CPlus5          605
#define CPlus6          606
#define NPlus0          700
#define NPlus1          701
#define NPlus2          702
#define NPlus3          703
#define NPlus4          704
#define NPlus5          705
#define NPlus6          706
#define NPlus7          707
#define AlPlus4         1304
#define userDefined   	9999999

#define Polygon    	1
#define Defined    	2
#define Beam		3
#define Channel    	4
#define BoostFrame	5
#define Circle    	6
#define Exp    		7

#define Constant   	0
#define Gaussian   	1
#define Polynomial   	2

#define byNumber	0
#define byDensity	1


// Token types (use typedef enum for clarity and portability)
typedef enum {
    TOK_NUM,
    TOK_VAR,
    TOK_OP
} TokenType;

// Token structure
typedef struct {
    TokenType type;
    double num;  // For TOK_NUM
    struct op_s *op;  // For TOK_OP
} Token;




typedef struct _LoadList  {
   int type;
   int species;
   int pair;
   double superP;
   double density;
   double numberInCell;
   double criticalDensity;
   double targetW;
   double num;      //exceeded number of particle which is less than 1
   int xnodes, ynodes, znodes;     //point numbers
   double *xn, *yn, *zn;      			//density (1 is P->density)
   double *xpoint, *ypoint, *zpoint;  	// point
   char **expr_x, **expr_y, **expr_z;    //function_XYZ
   Token **rpn_x, **rpn_y, **rpn_z;
   int *rpn_size_x, *rpn_size_y, *rpn_size_z;
   //center offset
   char **expr_x_y0, **expr_x_z0;
   Token **rpn_x_y0, **rpn_x_z0;
   int *rpn_size_x_y0, *rpn_size_x_z0;
   double givenMinPx;	//for saveParticle option

   //defined plasma
   int defineMode;
//   double *xPosition;
//   double *yPosition;
//   double *zPosition;
   int numDefined;
   double createGuard,deleteGuard;
   struct _DefPtcl *def;

//   double xLengDef;
//   double yLengDef;
   double RDef;
   int numDefPtcls;
//   double **define;
   int maxLoadTime;
   int minLoadTime;
   double minX,maxX;
   double minY,maxY;
   double minZ,maxZ;

   int pointPosition;   
   double px,py,pz,gamma;

   double mass;
   int charge;
   
   double temperature;
   
   //applying function
   double centerX;
   double centerY;
   double centerZ;
   double gaussCoefX;
   double polyCoefX;
   double gaussCoefYZ;
   double polyCoefYZ;
   int modeX;
   int modeYZ;

   //ionization
   int levels;
   int ionFinal;
   double givenMinA0;
   double *ionEnergy;
   double *ionW;

   //beam
//   double energy,posX,posY,posZ,spread,focus;
//   double sigY,sigZ,betaY,betaZ,emitY,emitZ;
//   double rjac,ratio;
   int loadingStep,gaussMode,numGam;
   double energy,posY,posZ,spread,peakCurr,eChirp;
   double sigY,sigZ,betaY,betaZ,alphaY,alphaZ,emitY,emitZ;
   double rjac,ratio,rjacComp,focus;
   double focalLy,focalLz;
   double sigX,gaussPower,posX;


   struct _LoadList *next;
} LoadList;

typedef struct _DefPtcl  {
   double xPos,yPos,zPos;
   double *define;
   int numDefPtcls;
   int flag;
   struct _DefPtcl *next;
} DefPtcl;

typedef struct _PlasmaLens  {
   int xnodes;       //longitudinal point number
   double *xn;       //longitudinal density (1 is P->density)
   double *xpoint;      //longitudinal point
   double radius;
   double current;

   struct _PlasmaLens *next;
} PlasmaLens;


