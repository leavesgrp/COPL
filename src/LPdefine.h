#include "LPconst.h"

typedef struct
{
  int mcol;
  int melm;

  int iterlmt;
  int cachesz;

  int save;
  int prlev;
  int prlop;
  int cross;

  int bemin;
  int idden;

  char *rhsset;
  char *rngset;
  char *bndset;
  char *objname;

  double lower;
  double upper;

  double bplus;
  double aijtol;
  double tolpiv;

  double tolgap;
  double tolp;
  double told;

  double tolx;
  double tols;

  double tolmu;
  double tolinf;
  double pcgtol;

  double Gamma;
  double beta;
} optpar;

typedef struct list *plst;

typedef struct list
{
  int addr;
  char *key;
  plst next;
} hashlst;

typedef struct
{
  int lssz;
  plst *list;
} hashtab;

typedef struct
{
  double low;
  double upp;
  double rhs;
  int nn0;
  int *ja;
  double *an;
} optrow;

typedef struct
{
  int nn0;
  int *ja;
  double *an;
} array;

typedef struct
{
  int mnrs;
  int nrow;
  int mems;
  array *ia;
} matrix;

typedef struct
{
  int maxr;
  int nrow;
  int maxn;
  int *iaat;
  int *naat;
  int *jaat;

  int *perm;
  int *invp;

  int nsnd;
  int *isnd;

  int ndnd;
  int ndns;
  int *idn;
  int *jdn;

  int *dmj;
  int plmt;
  int dnsp;

  int nu;
  int nju;
  int *iju;
  int *iu;
  int *jun;
  int *ju;
  double *ud;
  double *un;
} chol;

typedef struct
{

  int ista;
  int solsta;
  int infsta;

  int m;
  int n;
  int nz;

  int nlb;
  int nub;
  int nfr;
  int nfx;

  int nsk;
  int nrg;
  int ineq;
  int nlow;
  int nupp;
  int nfix;
  int nfre;
  int nbmi;
  int nbpl;

  int rnamesz;
  int cnamesz;

  hashtab *rowtab;
  hashtab *coltab;

  optrow *row;
  optrow *col;

  char *mpscard;
  char *prbname;
  char *objname;
  char *rowtype;
  char **rowname;
  char **colname;
  int *csk;

  int *imem;
  double *rmem;

  double *b;
  double *c;
  double *l;
  double *u;
  double *r;
  double obj;
  int *bid;

  char *rflg;
  char *cflg;

  int prtst;
  int nrld;
  int nrnul;
  int nrsg;
  int nrdb;
  int nrfc;
  int nrdm;
  int nrdu;
  int ncnul;
  int ncfx;
  int ncdm;
  int ncdu;

  double *rlb;
  double *rub;
  double *clb;
  double *cub;

  int items;
  int bytes;

  int iter;

  double pobj;
  double dobj;
  double rgap;
  double mu;

  double inf0;
  double rho;
  double mu0;
  double tau0;
  double kap0;
  double cx;
  double by;
  double rg;
  double rdf;

  double dfdt;
  double dfnm;

  double rddt;
  double ddtu;
  double rdnm;
  double rrdnm;

  double rpdt;
  double pdtu;
  double pnmu;
  double rpnm;
  double rrpnm;

  double rgdt;
  double smdt;
  double cgap;

  double infe;
  double nall;

  double *x;
  double *y;
  double *z;
  double *xu;
  double *yu;
  double *zu;
  double tau;
  double kap;
  double *xf;

  matrix *a;
  matrix *at;
  chol *cl;

  int nden;
  int *cid;
  matrix *dac;
  matrix *sac;
  double **ld;
  int *lperm;
  int *linvp;
  double *pcgbf;
} optmod;

extern FILE *fbin,
    *fmps,
    *fout,
    *fspc,
    *fstk,
    *fres;

extern char *nmps,
    *nout,
    *nspc;

extern optmod *ot;
extern optpar *pr;

#include "LPfunct.h"
