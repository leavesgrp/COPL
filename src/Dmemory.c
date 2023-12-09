#include "LPmatrx.h"

chol *ChlAlloc(int maxr, char *info)
{
  chol *r = NULL;

  if (maxr)
  {
    r = (chol *)calloc(1, sizeof(chol));
    if (!r)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }

    r->maxr = maxr;
    r->nrow = maxr;

    r->maxn = 0;
    r->iaat = iAlloc(maxr, info);
    r->naat = iAlloc(maxr, info);
    r->jaat = NULL;
    r->ud = dAlloc(maxr, info);
    r->nu = 0;
    r->nju = 0;
    r->iju = iAlloc(maxr, info);
    r->iu = iAlloc(maxr, info);
    r->jun = iAlloc(maxr, info);
    r->ju = NULL;
    r->un = NULL;
    r->perm = iAlloc(maxr, info);
    r->invp = iAlloc(maxr, info);
    r->nsnd = 0;
    r->isnd = iAlloc(maxr + 1, info);
  }
  return r;
} /* ChlAlloc */

void ChlFree(chol **cl)
{
  chol *r = *cl;

  if (*cl)
  {
    if (r->iaat)
      iFree(&r->iaat);
    if (r->naat)
      iFree(&r->naat);
    if (r->jaat)
      iFree(&r->jaat);
    if (r->ud)
      dFree(&r->ud);
    if (r->iju)
      iFree(&r->iju);
    if (r->iu)
      iFree(&r->iu);
    if (r->jun)
      iFree(&r->jun);
    if (r->ju)
      iFree(&r->ju);
    if (r->un)
      dFree(&r->un);
    if (r->perm)
      iFree(&r->perm);
    if (r->invp)
      iFree(&r->invp);
    if (r->isnd)
      iFree(&r->isnd);
    if (r->idn)
      iFree(&r->idn);
    if (r->dmj)
      iFree(&r->dmj);
    if (r->jdn)
      iFree(&r->jdn);
    free(r);
  }
  *cl = NULL;
} /* ChlFree */

void iPtAlloc(int m,
              int n,
              int *ipt[],
              char *info)
{
  int i;

  if (n)
  {
    for (i = 0; i < m; i++)
    {
      ipt[i] = (int *)calloc(n, sizeof(int));
      if (!ipt[i])
      {
        ErrorProc(NOT_MEMSPC, info);
        ShutDown();
        exit(0);
      }
    }
  }
} /* iPtAlloc */

void iPtFree(int m,
             int *ipt[])
{
  int i;

  for (i = 0; i < m; i++)
    iFree(&ipt[i]);
} /* iPtFree */

order *OdAlloc(int m,
               int nnz,
               char *info)
{
  order *r;

  r = (order *)calloc(1, sizeof(order));
  if (!r)
  {
    ErrorProc(NOT_MEMSPC, info);
    ShutDown();
    exit(0);
  }

  r->nnod = m;
  r->nn0 = nnz;

  r->adjn = iAlloc(nnz, info);
  r->rbeg = iAlloc(m, info);
  r->rexs = iAlloc(m, info);
  r->rlen = iAlloc(m, info);
  r->rend = iAlloc(m, info);
  r->pres = iAlloc(m, info);
  r->succ = iAlloc(m, info);

  return (r);
} /* OdAlloc */

void OdFree(order **od)
{
  order *r;

  if (*od)
  {
    r = *od;
    iFree(&r->adjn);
    iFree(&r->rbeg);
    iFree(&r->rexs);
    iFree(&r->rlen);
    iFree(&r->rend);
    iFree(&r->pres);
    iFree(&r->succ);
    free(*od);
    *od = NULL;
  }
} /* OdFree */
