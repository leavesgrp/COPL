/*
#define HPMACHINE
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>

#ifdef HPMACHINE
#include <string.h>
#include <sys/time.h>
#else
#include <ctype.h>
#endif
/*
#define TEST
*/
#define E3FT "%3.1E"
#define EGFT "%+16.9E"
#define IFMT "%d"
#define I4FT "%-4d"
#define I6FT "%6d"
#define PEFT "%.9E"
#define SFMT "%s"
#define TFMT "%.2f"
#define SPC4 "    "
#define SPC5 SPC4 " "
#define SPC6 SPC5 " "
#define SPC8 SPC5 "   "
#define SPCA "            "   // 10
#define SPCE "              " // 14

#define true 1
#define false 0
#define NumberPrimes 30
#define GoldFactor 0.6180339887
#define LineSize 255
#define StringSize 255
#define MaxPositive 1.E+30
#define MinPositive 1.E-30

#define COPLFILE_BIN "coplfile.bin"
#define COPLFILE_MPS "coplfile.mps"
#define COPLFILE_OUT "coplfile.out"
#define COPLFILE_RPT "coplfile.rpt"
#define COPLFILE_STK "coplfile.stk"

#if !defined(min)
#define min(a, b) ((a <= b) ? (a) : (b))
#endif
#if !defined(max)
#define max(a, b) ((a >= b) ? (a) : (b))
#endif

typedef enum
{
  START,
  DATAIN,
  CRUSH,
  PRSLV,
  SYBOL,
  NUMER,
  CROSS,
  ELAPS
} tim_set;

typedef enum
{
  BOUN = 0,
  CACH,
  CGAP,
  CROS,
  DENS,
  DINF,
  INFE,
  ITER,
  LOWE,
  NONZ,
  NUMC,
  NUME,
  OBJS,
  OPTI,
  PCGT,
  PIVT,
  PRES,
  PLOP,
  PGAP,
  PINF,
  RANG,
  RHSS,
  SAVE,
  STPF,
  UPPE,
  PEND = 26
} specid;

typedef enum
{
  SysWarning = 30,
  NOT_ROW_TYPE,
  NULL_ROW_SEC,
  DPL_ROW_NAME,
  DPL_COL_NAME,
  ZERO_ELEM,
  NOT_ROW_NAME,
  NOT_COL_NAME,
  DPL_RNG_ROW,
  NOT_BND_TYPE,

  SysError = 100,
  NOT_MEMSPC,
  NOT_DSKSPC,
  NOT_SCHFIL,
  NULL_FILE,
  NULL_SIZE,
  FAC_ERROR,
  OVER_FLOW,
  MIS_INDEX,
  NOT_NAMECARD,
  NOT_ROWSCARD,
  NOT_COLSCARD,
  NOT_ENDCARD,
  DPL_SYBENTRY,
  EXC_ROWLIM,
  EXC_COLLIM,
  EXC_NNZLIM,
  SZE_ERROR,
  SYB_PROC_ERR
} errmsg;

typedef enum
{
  UN_SOLVD = -1,
  OPT_FOUND,
  INF_PROB,
  UNB_PROB,
  EXC_ITER,
  FEA_PROB,
  UNK_SOL,
  UNK_PROB,
  NUM_PROB,
  CRS_TURN
} solmsg;

typedef enum
{
  NAME = 0,
  ROWS = 1,
  COLUMNS = 2,
  RHS = 3,
  RANGES = 4,
  BOUNDS = 5,
  ENDATA = 6,
  HEADCARD = 10,
  MIDCARD = 11
} mpssec;

typedef enum
{
  VALUE1 = 24,
  VALUE2 = 25
} fields;

typedef enum
{
  LOW = 0,
  UPP = 1,
  FIX = 2,
  FRE = 3,
  MIS = 4,
  PLS = 5,
  BOH = 6
} mpsbnd;

typedef enum
{
  ProcOk = 0,
  OutOfSpc,
  ItrLimit,
  ProcStall,
  NumDiff,
  FacError,
  UnkError
} procsta;

typedef enum
{
  FIRST = 0,
  BASIC,
  SPBAS,
  LOWER,
  UPPER,
  NUFRE,
  FIXED
} vartype;

typedef enum
{
  ChlOk = 0,
  ChlFail
} facsta;

typedef enum
{
  CDMN,
  CDUP,
  CFIX,
  CNUL,
  LBND,
  MATX,
  MINS,
  NMES,
  RDBL,
  RDMN,
  RFRC,
  RNUL,
  RDUP,
  RSNG,
  XSOL,
  YSOL
} stkid;
