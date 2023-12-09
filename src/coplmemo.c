#include "LPdefine.h"

static int PrimeSet[NumberPrimes] = {
    29, 229, 883, 1671, 2791,
    4801, 8629, 15289, 25303, 34843,
    65269, 99709, 129403, 147673, 166669,
    201403, 222163, 242729, 261431, 303491,
    320237, 402761, 501131, 602309, 701507,
    800999, 900551, 1000619, 1100837, 1200359};

char *cAlloc(int len, char *info)
{
  char *r = NULL;

  if (len)
  {
    r = (char *)calloc(len, sizeof(char));
    if (!r)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }
  }
  return r;
} /* cAlloc */

void cFree(char **x)
{
  char *r = *x;

  if (r)
  {
    free(r);
    *x = NULL;
  }
} /* cFree */

char *cAllocCopy(char *s, char *info)
{
  int len;
  char *r = NULL;

  if (s)
  {
    len = strlen(s);
    if (len)
    {
      r = cAlloc(len + 1, info);
      strcpy(r, s);
    }
  }
  return r;
} /* cAllocCopy */

int *iAlloc(int len, char *info)
{
  int *r = NULL;

  if (len)
  {
    r = (int *)calloc(len, sizeof(int));
    if (!r)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }
  }
  return r;
} /* iAlloc */

void iFree(int **x)
{
  int *r = *x;

  if (r)
  {
    free(r);
    *x = NULL;
  }
} /* iFree */

double *dAlloc(int len, char *info)
{
  double *r = NULL;

  if (len)
  {
    r = (double *)calloc(len, sizeof(double));
    if (!r)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }
  }
  return r;
} /* dAlloc */

void dFree(double **x)
{
  double *r = *x;

  if (r)
  {
    free(r);
    *x = NULL;
  }
} /* dFree */

double **dPtAlloc(int m, int n, char *info)
{
  int i;
  double **r = NULL;

  if (m)
  {
    r = (double **)calloc(m, sizeof(double *));
    if (!r)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }

    if (n)
    {
      r[0] = dAlloc(m * n, info);
      for (i = 1; i < m; i++)
        r[i] = r[i - 1] + n;
    }
  }

  return r;
} /* dPtAlloc */

void dPtFree(double ***m)
{
  double **r = *m;

  if (r)
  {
    dFree(&r[0]);
    free(r);
    *m = NULL;
  }
} /* dPtFree */

char **cPtAlloc(int nptr, char *info)
{
  char **r = NULL;

  if (nptr)
  {
    r = (char **)calloc(nptr, sizeof(char *));
    if (!r)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }
  }
  return r;
} /* cPtAlloc */

void cPtFree(char ***x, int n)
{
  int i;
  char **r = *x;

  if (r)
  {
    for (i = 0; i < n; i++)
      cFree(&r[i]);
    free(r);
    *x = NULL;
  }
} /* cPtFree */

hashtab *HashAlloc(int lssz, char *info)
{
  int i;
  hashtab *tab;

  tab = (hashtab *)malloc(sizeof(hashtab));
  if (!tab)
  {
    ErrorProc(NOT_MEMSPC, info);
    ShutDown();
    exit(0);
  }

  if (lssz)
  {
    for (i = 0; i < NumberPrimes; i++)
    {
      if (lssz < PrimeSet[i])
      {
        lssz = PrimeSet[i];
        break;
      }
    }

    tab->lssz = lssz;
    tab->list = (plst *)calloc(lssz, sizeof(plst));
    if (!tab->list)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }

    for (i = 0; i < lssz; i++)
      tab->list[i] = NULL;
  }
  return tab;
} /* HashAlloc */

void HashFree(hashtab **tab)
{
  int i;
  plst ptr, next;
  hashtab *r = *tab;

  if (!r)
    return;

  if (r->list)
  {
    for (i = 0; i < r->lssz; i++)
    {
      ptr = r->list[i];
      while (ptr)
      {
        next = ptr->next;
        cFree(&ptr->key);
        cFree((char **)&ptr);
        ptr = next;
      }
    }
    free(r->list);
    r->list = NULL;
  }
  free(r);
  *tab = NULL;
} /* HashFree */

optrow *RowAlloc(int m, int nz, char *info)
{
  optrow *r = NULL;

  if (m)
  {
    r = (optrow *)calloc(m, sizeof(optrow));
    if (!r)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }
    r->ja = iAlloc(nz, info);
    r->an = dAlloc(nz, info);
  }
  return r;
} /* RowAlloc */

void RowFree(optrow **x)
{
  optrow *r = *x;

  if (r)
  {
    iFree(&r->ja);
    dFree(&r->an);
    free(r);
    *x = NULL;
  }
} /* RowFree */

matrix *MtxAlloc(int m, int nz, char *info)
{
  matrix *r;

  r = (matrix *)calloc(1, sizeof(matrix));
  if (!r)
  {
    ErrorProc(NOT_MEMSPC, info);
    ShutDown();
    exit(0);
  }

  if (m)
  {
    r->ia = (array *)calloc(m, sizeof(array));
    if (!r)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }

    if (nz)
    {
      r->ia->ja = iAlloc(nz, info);
      r->ia->an = dAlloc(nz, info);
    }
  }
  r->mnrs = m;
  r->nrow = m;
  r->mems = nz;
  return r;
} /* MtxAlloc */

void MtxFree(matrix **a)
{
  matrix *r = *a;

  if (r)
  {
    if (r->ia)
    {
      iFree(&r->ia->ja);
      dFree(&r->ia->an);
      free(r->ia);
      r->ia = NULL;
    }

    r->mnrs = 0;
    r->nrow = 0;
    free(r);
  }
  *a = NULL;
} /* MtxFree */

int BndAlloc(optmod *opt)
{
  int m = opt->m, n = opt->n;

  opt->b = dAlloc(m, "for opt->b in BndAlloc");
  opt->r = dAlloc(m, "for opt->r in BndAlloc");
  opt->l = dAlloc(n, "for opt->l in BndAlloc");
  opt->u = dAlloc(n, "for opt->u in BndAlloc");
  return true;
} /* BndAlloc */
