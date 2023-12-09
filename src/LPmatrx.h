#include "LPdefine.h"

typedef struct
{
  int idep;
  int last;
  int hprt;
  int cent;
  int bwvt;
  int lwpr;
  int ntot;

  int *pbeg;
  int *pval;
  int *prev;
  int *succ;
} xlist;

typedef struct
{
  int nnod;
  int nn0;
  int raft;
  int head;
  int last;
  int ntot;

  int *adjn;
  int *rbeg;
  int *rexs;
  int *rlen;
  int *rend;
  int *pres;
  int *succ;
} order;

#include "LPmfunc.h"
