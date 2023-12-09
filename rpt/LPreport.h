#define HPMACHINE

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

#define true 1
#define false 0

#ifdef HPMACHINE
#define tfmt "%-13s"
#define dfmt "%+12e"
#else
#define tfmt "%-14s"
#define dfmt "%+13e"
#endif
#if !defined(min)
#define min(a, b) ((a <= b) ? (a) : (b))
#endif
#if !defined(max)
#define max(a, b) ((a >= b) ? (a) : (b))
#endif

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

typedef enum
{
  UNSOLVD = -1,
  OPTFUND,
  INFPROB,
  UNBPROB,
  NUMPROB
} solid;

typedef struct
{
  int nn0;
  int *ja;
  double *an;
} array;

typedef struct
{
  int ncol;
  int non0;
  array *cols;
} matrix;

typedef struct
{
  char sname[20];
  char *rp;

  int resp;
  int m;
  int n;
  int m0;
  int n0;
  int non0;
  int *subj;
  int *subi;

  int cdmn;
  int cdup;
  int cfix;
  int cnul;
  int lbnd;
  int mins;
  int rdbl;
  int rdmn;
  int rdup;
  int rfrc;
  int rnul;
  int rsng;

  double *x;
  double *y;
  double *z;

  double *x0;
  double *y0;
  double *z0;

  double *r;
  double *b;
  double *c;
  double *l;
  double *u;

  double *ai;
  double *aj;

  matrix *a;
} solpt;

extern FILE *fstk, *fbin,
    *fsol, *frpt;

/*
 * functions in coplpost.c
 */
char cGet(FILE *fp);
int iGet(FILE *fp);
double dGet(FILE *fp);
void PostSol(int, solpt *);

/*
 * functions in coplrepo.c
 */
void PrintHead(void);
int LeftDots(char *, int);
clock_t GetTime(void);
double TimeInSec(clock_t, clock_t);
char *cAlloc(int);
void cFree(char **);
int *iAlloc(int);
void iFree(int **);
double *dAlloc(int);
void dFree(double **);
matrix *mAlloc(int, int);
void mFree(matrix **);
