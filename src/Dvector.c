#include "LPmatrx.h"

int iSum(int n, int *v)
{
  int i, r = 0;

  for (i = 0; i < n; i++)
    r += v[i];
  return (r);
} /* iSum */

int getPos(int n,
           int i,
           int *v)
{
  int j;

  for (j = 0; j < n && i != v[j]; ++j)
    ;
  return (j);
} /* getPos */

void iZero(int n,
           int *x,
           int *sub)
{
  int i;

  if (sub)
  {
    for (i = 0; i < n; ++i)
      x[sub[i]] = 0;
  }
  else
    memset(x, 0, n * sizeof(n));
} /* iZero */

void iSet(int n,
          int a,
          int *x,
          int *s)
{
  int i;

  if (!s)
    for (i = 0; i < n; ++i)
      x[i] = a;
  else
    for (i = 0; i < n; ++i)
      x[s[i]] = a;
} /* iSet */

void iSwap(int i1,
           int i2,
           int *v)
{
  int temp;

  if (i1 < 0 || i2 < 0)
  {
    printf(" negative index.\n");
    ShutDown();
    exit(0);
  }
  temp = v[i1];
  v[i1] = v[i2];
  v[i2] = temp;
} /* iSwap */

void iCopy(int n,
           int *x,
           int *y)
{
  memcpy(y, x, n * sizeof(int));
} /* iCopy */

void fSort(int num,
           int maxi,
           int *val,
           int *fir,
           int *link)
{
  int k, v;

  iSet(maxi + 1, num, fir, NULL);

  for (k = 0; k < num; ++k)
  {
    v = val[k];
    if (v <= maxi)
    {
      link[k] = fir[v];
      fir[v] = k;
    }
  }
} /* fSort */

void dCopy(int n,
           double *s,
           double *d)
{
  memcpy(d, s, sizeof(double) * n);
} /* dCopy */

double dSums(int nnz,
             double x[],
             int indx[])
{
  double r = 0.0;
  int i;

  if (indx)
  {
    for (i = 0; i < nnz; ++i)
      r += x[indx[i]] * x[indx[i]];
  }
  else
  {
    for (i = 0; i < nnz; ++i)
      r += x[i] * x[i];
  }
  return (r);
} /* dSums */

void dZero(int n,
           double x[],
           int ind[])
{
  int i;

  if (!ind)
  {
    memset(x, 0, sizeof(double) * n);
  }

  else
  {
    for (i = 0; i < n; ++i)
      x[ind[i]] = 0.0;
  }
} /* dZero */

void dSwap(int i,
           int j,
           double *x)
{
  double temp;

  temp = x[i];
  x[i] = x[j];
  x[j] = temp;
} /* dSwap */

double dNorm0(int n,
              double x[],
              int indx[])
{
  int i;
  double r = 0.0;

  if (!indx)
    for (i = 0; i < n; ++i)
      r = max(r, fabs(x[i]));

  else
    for (i = 0; i < n; ++i)
      r = max(r, fabs(x[indx[i]]));

  return (r);
} /* dNorm0 */

double dNorm1(int n,
              double *x)
{
  int i;
  double r = 0.0;

  if (n)
  {
    for (i = 0; i < n; ++i)
      r += fabs(x[i]);
  }
  return (r);
} /* dNorm1 */

double dNorm2(int n,
              double *x)
{
  int i;
  double r = 0.0;

  for (i = 0; i < n; i++)
    r += x[i] * x[i];

  return (sqrt(fabs(r)));
} /* dNorm2 */

void addVect(int n,
             double alf,
             double *x,
             int *s,
             double *y)
{
  int i;

  if (n <= 0 || alf == 0.0)
    return;

  if (!s)
  {
    if (alf == 1.0)
    {
      for (i = 0; i < n; i++)
        y[i] += x[i];
    }
    else if (alf == -1.0)
    {
      for (i = 0; i < n; i++)
        y[i] -= x[i];
    }
    else
    {
      for (i = 0; i < n; i++)
        y[i] += alf * x[i];
    }
  }

  else
  {
    if (alf == 1.0)
    {
      for (i = 0; i < n; i++)
        y[s[i]] += x[i];
    }
    else if (alf == -1.0)
    {
      for (i = 0; i < n; i++)
        y[s[i]] -= x[i];
    }
    else
    {
      for (i = 0; i < n; i++)
        y[s[i]] += alf * x[i];
    }
  }
} /* addVect */

double svDot(array *x,
             double *y)
{
  int i;
  double r = 0.0;

  for (i = 0; i < x->nn0; ++i)
    r += x->an[i] * y[x->ja[i]];

  return (r);
} /* svDot */

void setArray(double alf,
              array *x,
              double *y)
{
  int k;

  if (x->nn0 <= 0 || alf == 0.0)
    return;

  if (alf == 1.0)
  {
    for (k = 0; k < x->nn0; k++)
      y[x->ja[k]] += x->an[k];
  }
  else if (alf == -1.0)
  {
    for (k = 0; k < x->nn0; k++)
      y[x->ja[k]] -= x->an[k];
  }
  else
  {
    for (k = 0; k < x->nn0; k++)
      y[x->ja[k]] += alf * x->an[k];
  }
} /* setArray */
