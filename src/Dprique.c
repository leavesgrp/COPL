#include "LPmatrx.h"

void XtFree(xlist **xt)
{
  xlist *r = *xt;

  if (r)
  {
    if (r->bwvt)
    {
      iFree(&r->pbeg);
      iFree(&r->pval);
      iFree(&r->prev);
      iFree(&r->succ);
    }
    free(r);
  }

  *xt = NULL;
} /* XtFree */

static void XtInit(xlist *xt)
{
  int k, m, n;

  xt->idep = xt->hprt + 1;
  xt->lwpr = xt->idep;
  xt->cent = xt->last;
  xt->ntot = 0;

  m = xt->idep;
  n = xt->last;
  for (k = 0; k < m; k++)
    xt->pbeg[k] = n;
  for (k = 0; k < n; k++)
  {
    xt->pval[k] = m;
    xt->prev[k] = n;
    xt->succ[k] = n;
  }
} /* XtInit */

xlist *XtAlloc(int last,
               int maxp)
{
  xlist *r;

  r = (xlist *)calloc(1, sizeof(xlist));
  if (r)
  {
    r->bwvt = true;
    r->last = last;
    r->hprt = maxp;
    r->ntot = 0;

    r->pbeg = iAlloc(maxp + 1, "PqAlloc");
    r->pval = iAlloc(last, "PqAlloc");
    r->prev = iAlloc(last, "PqAlloc");
    r->succ = iAlloc(last, "PqAlloc");

    if (r->pbeg && r->pval &&
        r->prev && r->succ)
      XtInit(r);
    else
      XtFree(&r);
  }

  return (r);
} /* XtAlloc */

int setXt(xlist *xt,
          int last,
          int maxp,
          int isze,
          int *imem)
{
  xt->hprt = maxp;
  xt->last = last;

  if (isze < (1 + xt->hprt + 3 * xt->last))
    return (false);

  xt->bwvt = false;
  xt->last = last;
  xt->hprt = maxp;

  if (imem == NULL)
    exit(0);

  xt->pbeg = imem;

  xt->pval = xt->pbeg + (maxp + 1);
  xt->prev = xt->pval + last;
  xt->succ = xt->prev + last;

  XtInit(xt);

  return (true);
} /* setXt */

int nextXt(xlist *xt)
{
  int t, last, maxp, *pbegp;

  last = xt->last;
  maxp = xt->hprt;

  if (xt->cent == last)
    return (false);
  else
  {
    if (xt->prev[xt->cent] != last)
      xt->cent = xt->prev[xt->cent];
    else
    {
      pbegp = xt->pbeg;

      for (t = xt->pval[xt->cent] + 1; t <= maxp && pbegp[t] == last; ++t)
        ;

      if (t > maxp)
        xt->cent = last;
      else
        xt->cent = xt->pbeg[t];
    }

    return (true);
  }
} /* nextXt */

void delXt(xlist *xt,
           int e)
{
  int p;

  if (xt->pval[e] != xt->idep)
  {
    if (xt->ntot <= 0)
    {
      printf("\n\n system error, xt->ntot<=0.\n");
      ShutDown();
      exit(0);
    }

    --xt->ntot;

    if (xt->cent == e)
    {
      if (xt->ntot)
        nextXt(xt);
      else
        xt->cent = xt->last;
    }

    p = xt->pval[e];
    xt->pval[e] = xt->idep;

    if (xt->succ[e] != xt->last)
      xt->prev[xt->succ[e]] = xt->prev[e];
    else
      xt->pbeg[p] = xt->prev[e];

    if (xt->prev[e] != xt->last)
      xt->succ[xt->prev[e]] = xt->succ[e];

    if (xt->pbeg[p] == xt->last && xt->lwpr == p)
    {
      xt->lwpr = xt->idep;
      if (xt->ntot)
      {
        for (++p; p <= xt->hprt; ++p)
        {
          if (xt->pbeg[p] != xt->last)
          {
            xt->lwpr = p;
            break;
          }
        }
      }
    }
  }
} /* delXt */

int infXt(xlist *xt)
{
  if (xt->lwpr == xt->idep)
  {
    if (xt->ntot != 0)
    {
      printf("\n\n system error, xt->ntot!=0.\n");
      ShutDown();
      exit(0);
    }

    xt->cent = xt->last;
    return (false);
  }
  else
  {
    if (xt->ntot <= 0)
    {
      printf("\n\n system error, xt->ntot<=0.\n");
      ShutDown();
      exit(0);
    }

    xt->cent = xt->pbeg[xt->lwpr];
    return (true);
  }
} /* infXt */

int supXt(xlist *xt,
          int p)
{
  if (0 <= p && p <= xt->hprt)
  {
    if (xt->pbeg[p] == xt->last)
      return (false);
    else
    {
      xt->cent = xt->pbeg[p];
      return (true);
    }
  }
  else
  {
    printf("\n\n arglast error, p<0 or p>xt->hprt.\n");
    ShutDown();
    exit(0);
  }

  return (false);
} /* supXt */

int succXt(xlist *xt,
           int *e)
{
  int last;

  last = xt->last;
  if (xt->cent != last && xt->prev[xt->cent] != last)
  {
    xt->cent = xt->prev[xt->cent];
    *e = xt->cent;
    return (true);
  }
  return (false);
} /* succXt */

void putXt(xlist *xt,
           int e,
           int p)
{
  if (0 <= e && e < xt->last && 0 <= p && p <= xt->hprt)
  {
    delXt(xt, e);

    ++xt->ntot;

    xt->pval[e] = p;
    xt->prev[e] = xt->pbeg[p];
    xt->succ[e] = xt->last;

    if (xt->pbeg[p] != xt->last)
      xt->succ[xt->pbeg[p]] = e;

    xt->pbeg[p] = e;
    xt->lwpr = min(p, xt->lwpr);
  }

  else
  {
    if (e < 0 || e >= xt->last)
    {
      printf("\n\n argulast error, e<0 or e>=xt->last.\n");
      ShutDown();
      exit(0);
    }
    else
    {
      printf("\n\n argulast error, 0<=e<xt->last.\n");
      ShutDown();
      exit(0);
    }
  }
} /* putXt */

void miXt(xlist *xt,
          int e,
          int dec)
{
  putXt(xt, e, xt->pval[e] - dec);
} /* miXt */

void plXt(xlist *xt,
          int e,
          int inc)
{
  putXt(xt, e, xt->pval[e] + inc);
} /* plXt */

int getXt(xlist *xt,
          int *e,
          int *p)
{
  if (xt->cent > xt->last)
  {
    printf("\n\n internal error, xt->cent>xt->last.\n");
    ShutDown();
    exit(0);
  }
  if (xt->cent == xt->last)
    return (false);
  else
  {
    *e = xt->cent;
    *p = xt->pval[*e];
    return (true);
  }
} /* getXt */
