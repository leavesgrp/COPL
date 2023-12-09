#include "LPmatrx.h"

//#define TIME_COUNT

static void setXYind(int nnz,
                     double *y,
                     double *x,
                     int *s)
{
  int i;

  for (i = 0; i < nnz; ++i)
  {
    x[i] = y[s[i]];
    y[s[i]] = 0.0;
  }
} /* setXYind */

static int setColi(chol *cl,
                   int i,
                   double *ai)
{
  cl->ud[i] += ai[i];
  ai[i] = 0.0;
  setXYind(cl->jun[i], ai, cl->un + cl->iu[i], cl->ju + cl->iju[i]);

  return (ChlOk);
} /* setColi */

static void modDiag(double *v,
                    int n,
                    double *uni)
{
  int j;

  if (*v < 1.0e-13)
  {
    *v = 1.0;
    for (j = 0; j < n; j++)
      uni[j] = 0.0;
  }
} /* modDiag */

static void modSnods(int m,
                     int n,
                     int s,
                     double diaga[],
                     double *a,
                     int fira[],
                     double diagb[],
                     double *b,
                     int firb[])
{
  int i, k, t, sze;
  double rtemp1, rtemp2, rtemp3, rtemp4,
      rtemp5, rtemp6, rtemp7, rtemp8,
      rtemp9, rtemp10, rtemp11, rtemp12,
      *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8,
      *a9, *a10, *a11, *a12, *b0;

  i = 0;
  for (; i < s; ++i)
  {
    b0 = b + firb[i];
    sze = m - i - 1;
    k = 0;

    for (; k + 11 < n; k += 12)
    {
      a1 = a + (fira[k + 0] + i);
      a2 = a + (fira[k + 1] + i);
      a3 = a + (fira[k + 2] + i);
      a4 = a + (fira[k + 3] + i);
      a5 = a + (fira[k + 4] + i);
      a6 = a + (fira[k + 5] + i);
      a7 = a + (fira[k + 6] + i);
      a8 = a + (fira[k + 7] + i);
      a9 = a + (fira[k + 8] + i);
      a10 = a + (fira[k + 9] + i);
      a11 = a + (fira[k + 10] + i);
      a12 = a + (fira[k + 11] + i);

      rtemp1 = *a1 / diaga[k + 0];
      rtemp2 = *a2 / diaga[k + 1];
      rtemp3 = *a3 / diaga[k + 2];
      rtemp4 = *a4 / diaga[k + 3];
      rtemp5 = *a5 / diaga[k + 4];
      rtemp6 = *a6 / diaga[k + 5];
      rtemp7 = *a7 / diaga[k + 6];
      rtemp8 = *a8 / diaga[k + 7];
      rtemp9 = *a9 / diaga[k + 8];
      rtemp10 = *a10 / diaga[k + 9];
      rtemp11 = *a11 / diaga[k + 10];
      rtemp12 = *a12 / diaga[k + 11];

      diagb[i] -= rtemp1 * (*a1) + rtemp2 * (*a2) + rtemp3 * (*a3) + rtemp4 * (*a4) + rtemp5 * (*a5) + rtemp6 * (*a6) + rtemp7 * (*a7) + rtemp8 * (*a8) + rtemp9 * (*a9) + rtemp10 * (*a10) + rtemp11 * (*a11) + rtemp12 * (*a12);

      ++a1;
      ++a2;
      ++a3;
      ++a4;
      ++a5;
      ++a6;
      ++a7;
      ++a8;
      ++a9;
      ++a10;
      ++a11;
      ++a12;

      for (t = 0; t < sze; ++t)
        b0[t] -= rtemp1 * a1[t] + rtemp2 * a2[t] + rtemp3 * a3[t] + rtemp4 * a4[t] + rtemp5 * a5[t] + rtemp6 * a6[t] + rtemp7 * a7[t] + rtemp8 * a8[t] + rtemp9 * a9[t] + rtemp10 * a10[t] + rtemp11 * a11[t] + rtemp12 * a12[t];
    }

    for (; k + 7 < n; k += 8)
    {
      a1 = a + (fira[k + 0] + i);
      a2 = a + (fira[k + 1] + i);
      a3 = a + (fira[k + 2] + i);
      a4 = a + (fira[k + 3] + i);
      a5 = a + (fira[k + 4] + i);
      a6 = a + (fira[k + 5] + i);
      a7 = a + (fira[k + 6] + i);
      a8 = a + (fira[k + 7] + i);

      rtemp1 = *a1 / diaga[k + 0];
      rtemp2 = *a2 / diaga[k + 1];
      rtemp3 = *a3 / diaga[k + 2];
      rtemp4 = *a4 / diaga[k + 3];
      rtemp5 = *a5 / diaga[k + 4];
      rtemp6 = *a6 / diaga[k + 5];
      rtemp7 = *a7 / diaga[k + 6];
      rtemp8 = *a8 / diaga[k + 7];

      diagb[i] -= rtemp1 * (*a1) + rtemp2 * (*a2) + rtemp3 * (*a3) + rtemp4 * (*a4) + rtemp5 * (*a5) + rtemp6 * (*a6) + rtemp7 * (*a7) + rtemp8 * (*a8);

      ++a1;
      ++a2;
      ++a3;
      ++a4;
      ++a5;
      ++a6;
      ++a7;
      ++a8;

      for (t = 0; t < sze; ++t)
        b0[t] -= rtemp1 * a1[t] + rtemp2 * a2[t] + rtemp3 * a3[t] + rtemp4 * a4[t] + rtemp5 * a5[t] + rtemp6 * a6[t] + rtemp7 * a7[t] + rtemp8 * a8[t];
    }

    for (; k + 3 < n; k += 4)
    {
      a1 = a + (fira[k + 0] + i);
      a2 = a + (fira[k + 1] + i);
      a3 = a + (fira[k + 2] + i);
      a4 = a + (fira[k + 3] + i);

      rtemp1 = *a1 / diaga[k + 0];
      rtemp2 = *a2 / diaga[k + 1];
      rtemp3 = *a3 / diaga[k + 2];
      rtemp4 = *a4 / diaga[k + 3];

      diagb[i] -= rtemp1 * (*a1) + rtemp2 * (*a2) + rtemp3 * (*a3) + rtemp4 * (*a4);

      ++a1;
      ++a2;
      ++a3;
      ++a4;

      for (t = 0; t < sze; ++t)
        b0[t] -= rtemp1 * a1[t] + rtemp2 * a2[t] + rtemp3 * a3[t] + rtemp4 * a4[t];
    }

    for (; k + 1 < n; k += 2)
    {
      a1 = a + (fira[k + 0] + i);
      a2 = a + (fira[k + 1] + i);

      rtemp1 = *a1 / diaga[k + 0];
      rtemp2 = *a2 / diaga[k + 1];

      diagb[i] -= rtemp1 * (*a1) + rtemp2 * (*a2);

      ++a1;
      ++a2;

      for (t = 0; t < sze; ++t)
        b0[t] -= rtemp1 * a1[t] + rtemp2 * a2[t];
    }

    for (; k < n; ++k)
    {
      a1 = a + (fira[k + 0] + i);
      rtemp1 = *a1 / diaga[k + 0];
      diagb[i] -= rtemp1 * (*a1);
      ++a1;

      for (t = 0; t < sze; ++t)
        b0[t] -= rtemp1 * a1[t];
    }
  }
} /* modSnods */

static void modInts(chol *cl,
                    int snde,
                    int f,
                    int l,
                    int uf,
                    int ul,
                    int *i1nrow)
{
  int k, *jun, *iu, *isnd;
  double *diag, *un;

  jun = cl->jun;
  iu = cl->iu;
  isnd = cl->isnd;
  diag = cl->ud;
  un = cl->un;

  if (f == l || uf == ul)
    return;

  f += isnd[snde];
  l += isnd[snde];
  uf += isnd[snde];
  ul += isnd[snde];

  for (k = f; k < l; ++k)
    i1nrow[k - f] = iu[k] + uf - k - 1;

  modSnods(1 + jun[uf], l - f, ul - uf, diag + f, un, i1nrow, diag + uf, un, iu + uf);
} /* modInts */

static int decSnode(chol *cl,
                    int snde,
                    int f,
                    int l,
                    int *i1nrow)
{
  int itemp, k;

  if (f == l)
    return (ChlOk);

  itemp = cl->isnd[snde] + f;

  modDiag(&cl->ud[itemp], cl->jun[itemp], cl->un + cl->iu[itemp]);

  if (cl->ud[itemp] <= 1.0e-13)
  {
    printf("\nsingular d[%d]=%e",
           cl->isnd[snde] + f, cl->ud[cl->isnd[snde] + f]);
    return (ChlFail);
  }

  for (k = f + 1; k < l; ++k)
  {
    modInts(cl, snde, f, k, k, k + 1, i1nrow);

    itemp = cl->isnd[snde] + k;

    modDiag(&cl->ud[itemp], cl->jun[itemp], cl->un + cl->iu[itemp]);

    if (cl->ud[itemp] <= 1.0e-13)
    {
      printf("\nsingular d[%d]=%e",
             cl->isnd[snde] + k, cl->ud[cl->isnd[snde] + k]);
      return (ChlFail);
    }
  }

  return (ChlOk);
} /* decSnode */

static void updSnods(int m,
                     int n,
                     int s,
                     double *diaga,
                     double *a,
                     int *fira,
                     double *diagb,
                     double *b,
                     int *firb,
                     int *subb)
{
  int i, j, k, t, u, sze, delay,
      *ls;
  double rtemp1, rtemp2, rtemp3, rtemp4,
      rtemp5, rtemp6, rtemp7, rtemp8,
      rtemp9, rtemp10, rtemp11, rtemp12,
      *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8,
      *a9, *a10, *a11, *a12, *b0;

  if (m < s)
    exit(0);

  if (m == 0 || n == 0)
    return;

  for (i = 0; i < s; ++i)
  {
    j = subb[i];
    u = j - subb[0];

    b0 = b + firb[u];

    delay = j + 1;
    sze = m - i - 1;
    ls = subb + i + 1;

    k = 0;

    for (; k + 11 < n; k += 12)
    {
      a1 = a + (fira[k + 0] + i);
      a2 = a + (fira[k + 1] + i);
      a3 = a + (fira[k + 2] + i);
      a4 = a + (fira[k + 3] + i);
      a5 = a + (fira[k + 4] + i);
      a6 = a + (fira[k + 5] + i);
      a7 = a + (fira[k + 6] + i);
      a8 = a + (fira[k + 7] + i);
      a9 = a + (fira[k + 8] + i);
      a10 = a + (fira[k + 9] + i);
      a11 = a + (fira[k + 10] + i);
      a12 = a + (fira[k + 11] + i);

      rtemp1 = *a1 / diaga[k + 0];
      rtemp2 = *a2 / diaga[k + 1];
      rtemp3 = *a3 / diaga[k + 2];
      rtemp4 = *a4 / diaga[k + 3];
      rtemp5 = *a5 / diaga[k + 4];
      rtemp6 = *a6 / diaga[k + 5];
      rtemp7 = *a7 / diaga[k + 6];
      rtemp8 = *a8 / diaga[k + 7];
      rtemp9 = *a9 / diaga[k + 8];
      rtemp10 = *a10 / diaga[k + 9];
      rtemp11 = *a11 / diaga[k + 10];
      rtemp12 = *a12 / diaga[k + 11];

      diagb[u] -= rtemp1 * (*a1) + rtemp2 * (*a2) + rtemp3 * (*a3) + rtemp4 * (*a4) + rtemp5 * (*a5) + rtemp6 * (*a6) + rtemp7 * (*a7) + rtemp8 * (*a8) + rtemp9 * (*a9) + rtemp10 * (*a10) + rtemp11 * (*a11) + rtemp12 * (*a12);

      ++a1;
      ++a2;
      ++a3;
      ++a4;
      ++a5;
      ++a6;
      ++a7;
      ++a8;
      ++a9;
      ++a10;
      ++a11;
      ++a12;

      for (t = 0; t < sze; ++t)
        b0[ls[t] - delay] -= rtemp1 * a1[t] + rtemp2 * a2[t] + rtemp3 * a3[t] + rtemp4 * a4[t] + rtemp5 * a5[t] + rtemp6 * a6[t] + rtemp7 * a7[t] + rtemp8 * a8[t] + rtemp9 * a9[t] + rtemp10 * a10[t] + rtemp11 * a11[t] + rtemp12 * a12[t];
    }
    for (; k + 7 < n; k += 8)
    {
      a1 = a + (fira[k + 0] + i);
      a2 = a + (fira[k + 1] + i);
      a3 = a + (fira[k + 2] + i);
      a4 = a + (fira[k + 3] + i);
      a5 = a + (fira[k + 4] + i);
      a6 = a + (fira[k + 5] + i);
      a7 = a + (fira[k + 6] + i);
      a8 = a + (fira[k + 7] + i);

      rtemp1 = *a1 / diaga[k + 0];
      rtemp2 = *a2 / diaga[k + 1];
      rtemp3 = *a3 / diaga[k + 2];
      rtemp4 = *a4 / diaga[k + 3];
      rtemp5 = *a5 / diaga[k + 4];
      rtemp6 = *a6 / diaga[k + 5];
      rtemp7 = *a7 / diaga[k + 6];
      rtemp8 = *a8 / diaga[k + 7];

      diagb[u] -= rtemp1 * (*a1) + rtemp2 * (*a2) + rtemp3 * (*a3) + rtemp4 * (*a4) + rtemp5 * (*a5) + rtemp6 * (*a6) + rtemp7 * (*a7) + rtemp8 * (*a8);

      ++a1;
      ++a2;
      ++a3;
      ++a4;
      ++a5;
      ++a6;
      ++a7;
      ++a8;

      for (t = 0; t < sze; ++t)
        b0[ls[t] - delay] -= rtemp1 * a1[t] + rtemp2 * a2[t] + rtemp3 * a3[t] + rtemp4 * a4[t] + rtemp5 * a5[t] + rtemp6 * a6[t] + rtemp7 * a7[t] + rtemp8 * a8[t];
    }

    for (; k + 3 < n; k += 4)
    {
      a1 = a + (fira[k + 0] + i);
      a2 = a + (fira[k + 1] + i);
      a3 = a + (fira[k + 2] + i);
      a4 = a + (fira[k + 3] + i);

      rtemp1 = *a1 / diaga[k + 0];
      rtemp2 = *a2 / diaga[k + 1];
      rtemp3 = *a3 / diaga[k + 2];
      rtemp4 = *a4 / diaga[k + 3];

      diagb[u] -= rtemp1 * (*a1) + rtemp2 * (*a2) + rtemp3 * (*a3) + rtemp4 * (*a4);

      ++a1;
      ++a2;
      ++a3;
      ++a4;

      for (t = 0; t < sze; ++t)
        b0[ls[t] - delay] -= rtemp1 * a1[t] + rtemp2 * a2[t] + rtemp3 * a3[t] + rtemp4 * a4[t];
    }

    for (; k + 1 < n; k += 2)
    {
      a1 = a + (fira[k + 0] + i);
      a2 = a + (fira[k + 1] + i);

      rtemp1 = *a1 / diaga[k + 0];
      rtemp2 = *a2 / diaga[k + 1];

      diagb[u] -= rtemp1 * (*a1) + rtemp2 * (*a2);

      ++a1;
      ++a2;

      for (t = 0; t < sze; ++t)
        b0[ls[t] - delay] -= rtemp1 * a1[t] + rtemp2 * a2[t];
    }

    for (; k < n; ++k)
    {
      a1 = a + (fira[k + 0] + i);
      rtemp1 = *a1 / diaga[k + 0];
      diagb[u] -= rtemp1 * (*a1);
      a1++;

      for (t = 0; t < sze; ++t)
        b0[ls[t] - delay] -= rtemp1 * a1[t];
    }
  }
} /*  updSnods */

static void modExts(chol *cl,
                    int snde,
                    int usnde,
                    int f,
                    int l,
                    int start,
                    int *i1nrow)
{
  int k, sze, *ls, *isnd, *jun, *ju, *iju, *iu;
  double *diag, *un;

  isnd = cl->isnd;
  jun = cl->jun;
  ju = cl->ju;
  iju = cl->iju;
  iu = cl->iu;
  diag = cl->ud;
  un = cl->un;

  f += isnd[snde];
  l += isnd[snde];

  if (usnde == cl->nsnd - 1)
  {
    if (ju[iju[f] + start] < isnd[usnde])
    {
      printf("\n Index error");
      exit(0);
    }

    if (cl->dnsp)
      exit(0);

    ls = ju + iju[f] + start;
    sze = jun[f] - start;

    for (k = f; k < l; ++k)
      i1nrow[k - f] = iu[k] + start - (k - f);

    updSnods(sze, l - f, sze, diag + f, un, i1nrow, diag + ls[0], un, iu + ls[0], ls);
  }
  else
    exit(0);
} /* modExts */

static void ClFwrd(chol *cl,
                   int snde,
                   int f,
                   int l,
                   int *i1nrow)
{
  int j, s, t, u, k, stops, offset, sze, itemp,
      *ls0, *ls1, *map, *isnd, *jun, *iu, *iju, *ju;
  double rtemp1, *l0, *l1, *diag, *un;

  jun = cl->jun;
  iu = cl->iu;
  iju = cl->iju;
  ju = cl->ju;
  diag = cl->ud;
  un = cl->un;
  isnd = cl->isnd;
  map = i1nrow;

  if (f > isnd[snde + 1] - isnd[snde])
  {
    printf("\n ClFwrd");
    exit(0);
  }

  if (f == l)
    return;

  f += isnd[snde];
  l += isnd[snde];

  offset = isnd[snde + 1] - f - 1;
  sze = jun[f] - offset;
  ls1 = ju + iju[f] + offset;

  if (f + 1 == l)
  {
    l1 = un + iu[f] + offset;

    stops = sze;
    for (t = 0; t < sze; ++t)
    {
      j = ls1[0];

      if (j >= isnd[cl->nsnd - 1])
        break;

      rtemp1 = l1[0] / diag[f];
      diag[j] -= rtemp1 * l1[0];
      ++l1;

      l0 = un + iu[j];
      ls0 = ju + iju[j];

      ++ls1;
      --stops;

      if (stops && ls1[stops - 1] == ls0[stops - 1])
      {
        for (s = 0; s < stops; ++s)
          l0[s] -= rtemp1 * l1[s];
      }
      else
      {
        for (s = 0, u = 0; s < stops; ++u)
        {
          if (ls0[u] == ls1[s])
          {
            l0[u] -= rtemp1 * l1[s];
            ++s;
          }
        }
      }
    }

    if (t < sze && !cl->dnsp)
      modExts(cl, snde, cl->nsnd - 1, f - isnd[snde], l - isnd[snde], t, i1nrow);
  }
  else
  {
    stops = sze;
    for (t = 0; t < sze; ++t, ++offset)
    {
      j = ls1[0];
      if (j >= isnd[cl->nsnd - 1])
      {
        if (!cl->dnsp)
          modExts(cl, snde, cl->nsnd - 1, f - isnd[snde], l - isnd[snde], offset, i1nrow);
        break;
      }

      ls0 = ju + iju[j];
      l0 = un + iu[j];

      ++ls1;
      --stops;

      k = f;
      itemp = offset + f;

      if (stops && ls1[stops - 1] == ls0[stops - 1])
      {
        for (k = f; k < l; ++k)
          map[k - f] = iu[k] + itemp - k;

        modSnods(1 + stops, l - f, 1, diag + f, un, map, diag + j, un, iu + j);
      }
      else
      {
        map[l] = 0;
        for (s = 0, u = 0; s < stops; ++u)
        {
          if (ls1[s] == ls0[u])
          {
            map[1 + l + s] = 1 + u;
            ++s;
          }
        }

        for (k = f; k < l; ++k)
          map[k - f] = iu[k] + itemp - k;

        updSnods(1 + stops, l - f, 1, diag + f, un, map, diag + j, un, iu + j, map + l);
      }
    }
  }
} /* ClFwrd */

static int decDense(chol *cl,
                    int *i1nrow,
                    double *r1nrow)
{
  int c, d, j, s, t, sncl, k, k0, m, cacsze, sze, offset,
      *isnd, *jun, *iju, *iu, *ju, *idn, *jdn, *dmj,
      *ls, sresp;
  double *diag, *un;
#ifdef TIME_COUNT
  clock_t Tm;
  int cnt;
  Tm = GetTime();
#endif

  isnd = cl->isnd,
  jun = cl->jun,
  iju = cl->iju,
  iu = cl->iu,
  ju = cl->ju,
  idn = cl->idn,
  jdn = cl->jdn,
  dmj = cl->dmj;
  diag = cl->ud,
  un = cl->un;

  cacsze = pr->cachesz;

  if (cl->dnsp)
  {
    for (d = 0; d < cl->ndnd; ++d)
    {
      c = 0;
      for (k = idn[d]; k < idn[d + 1]; ++k)
      {
        offset = dmj[k];
        s = jdn[k];
        if (ju[iju[isnd[s]] + offset] < isnd[cl->nsnd - 1])
        {
          printf("\nindex error1");
          exit(0);
        }

        for (j = isnd[s]; j < isnd[s + 1]; ++c, ++j)
        {
          r1nrow[c] = diag[j];
          i1nrow[c] = iu[j] + offset - (j - isnd[s]);

          if (ju[iju[j] + offset - (j - isnd[s])] < isnd[cl->nsnd - 1])
          {
            printf("\nindex error");
            exit(0);
          }
        }
      }

      if (c)
      {
        k = idn[d];
        s = jdn[k];
        m = jun[isnd[s]] - dmj[k];
        ls = ju + iju[isnd[s]] + dmj[k];
        if (m)
        {
          for (k = 0; k < c;)
          {
            t = cacsze / (m * sizeof(double));
            t = max(t, 1);
            t = min(c - k, t);

            updSnods(m, t, m, r1nrow + k, un, i1nrow + k,
                     diag + ls[0], un, iu + ls[0], ls);
            k += t;
          }
        }
      }
    }
  }

  s = cl->nsnd - 1;

  sncl = cl->isnd[s + 1] - cl->isnd[s];
#ifdef TIME_COUNT
  printf("decDense-A: %d ms\n", GetTime() - Tm);
  cnt = 0;
#endif
  for (k = 0; k < sncl;)
  {
    k0 = k;
    for (sze = 0; sze < cacsze && k < sncl; ++k)
      sze += jun[isnd[s] + k] * sizeof(double);

    if (k == k0)
      ++k;
    else if (k >= k0 + 2 && sze > cacsze)
      --k;

    if (k > sncl)
      exit(0);

    sresp = decSnode(cl, s, k0, k, i1nrow); // HOTSPOT
    if (sresp != ChlOk)
      return (sresp);

    modInts(cl, s, k0, k, k, sncl, i1nrow);
  }
#ifdef TIME_COUNT
  printf("decDense-B: %d ms\n", GetTime() - Tm);
#endif
  return (ChlOk);
} /* decDense */

int puchChol(chol *cl,
             int *i1nrow,
             double *zr1nrow)
{
  int s, sncl, k, k0, cacsze, sze, *isnd, *jun;
  int sresp;

  isnd = cl->isnd;
  jun = cl->jun;
  cacsze = pr->cachesz;

  for (s = 0; s + 1 < cl->nsnd; ++s)
  {
    sncl = cl->isnd[s + 1] - cl->isnd[s];

    for (k = 0; k < sncl;)
    {
      k0 = k;
      for (sze = 0; sze <= cacsze && k < sncl; ++k)
        sze += jun[isnd[s] + k] * sizeof(double);

      if (k == k0)
        k++;
      else if (k >= k0 + 2 && sze > cacsze)
        k--;

      if (k > sncl)
        exit(0);

      sresp = decSnode(cl, s, k0, k, i1nrow);
      if (sresp != ChlOk)
        return (sresp);

      modInts(cl, s, k0, k, k, sncl, i1nrow);

      ClFwrd(cl, s, k0, k, i1nrow);
    }
  }
  decDense(cl, i1nrow, zr1nrow); // HOTSPOT
  return (ChlOk);
} /* puchChol */

/*
int sfc_numfac(chol   *cl,
               int    i1nrow[],
               double r1nrow[])
{
  int     nrow=cl->nrow;

  dZero(nrow,r1nrow,NULL);
  return (puchChol(cl,i1nrow,r1nrow));
} / sfc_numfac */

int decChol(chol *cl,
            double *diag,
            matrix *at,
            int isze,
            int *imem,
            int rsze,
            double *rmem)
{
  int k, j, i, ii, next, nnzaj, beg, end, nrow, ncol,
      nil, *invp, *subaj, *jfir, *jnex, *jcur;
  double rtemp, *valaj, *bval, *rnex;
  array *aj;
  int ret;

#ifdef TIME_COUNT
  clock_t Tm;
  int COUNT1, COUNT2;
  Tm = GetTime();
  COUNT1 = 0;
  COUNT2 = 0;
#endif

  nrow = cl->nrow;
  ncol = at->nrow;
  nil = ncol;

  if (isze < 2 * ncol + nrow || rsze < ncol + 2 * nrow)
    return (OutOfSpc);

  jcur = imem;
  jfir = jcur + ncol;
  jnex = jfir + nrow;
  invp = cl->invp;
  bval = rmem;
  rnex = bval + nrow;

  iSet(nrow, nil, jfir, NULL);

  for (j = 0; j < ncol; ++j)
  {
    jcur[j] = 0;

    aj = at->ia + j;
    nnzaj = aj->nn0;
    subaj = aj->ja;

    if (nnzaj)
    {
      i = invp[subaj[0]];
      jnex[j] = jfir[i];
      jfir[i] = j;
    }
    else
      jnex[j] = nil;
  }
  dZero(nrow, bval, NULL);

#ifdef TIME_COUNT
  printf("Chol-A: %d ms\n", GetTime() - Tm);
#endif
  for (k = 0; k < nrow; ++k)
  {
    for (j = jfir[k]; j != nil; j = next)
    {
      next = jnex[j];
      aj = at->ia + j;
      subaj = aj->ja;
      valaj = aj->an;

      rtemp = aj->an[jcur[j]] / diag[j];
      beg = jcur[j];
      end = aj->nn0;

#ifdef TIME_COUNT
      COUNT1++;
      COUNT2 += end - beg;
#endif
      // printf("%d %d\n", beg, end);
#pragma omp parallel for num_threads(1)
      for (i = beg; i < end; ++i)
      {
        bval[invp[subaj[i]]] += rtemp * valaj[i]; // DENSE: HOTSPOT
                                                  // if (invp[subaj[i]] != i) printf("%d %d %d\n", i, subaj[i], invp[subaj[i]]);
      }
      /*
      for(i=beg; i<end; ++i)
        for(ii=i+1; ii<end; ++ii)
          if (invp[subaj[i]] == invp[subaj[ii]])
            printf("bug\n");
            */
      jcur[j]++;

      if (jcur[j] < end)
      {
        i = invp[aj->ja[jcur[j]]];
        jnex[j] = jfir[i];
        jfir[i] = j;
      }
    }

    setColi(cl, k, bval);
  }
#ifdef TIME_COUNT
  printf("[%d, %d]\n", COUNT1, COUNT2);
  printf("Chol-B: %d ms\n", GetTime() - Tm);
#endif
  for (k = 0; k < nrow; ++k)
  {
    rmem[k] = 0.0;
    rnex[k] = cl->ud[cl->invp[k]];
  }
  ret = puchChol(cl, imem, rmem); // SPARSE: HOTSPOT
#ifdef TIME_COUNT
  printf("Chol-C: %d ms\n", GetTime() - Tm);
#endif
  return (ret);
} /* decChol */

static void subsFwrd(chol *cl,
                     int snde,
                     int f,
                     int l,
                     double *x)
{
  int i, t, sze, *ls, *isnd, *iju, *iu, *ju;
  double xi, *l1, *diag, *un;

  isnd = cl->isnd,
  iju = cl->iju,
  iu = cl->iu,
  ju = cl->ju;
  diag = cl->ud,
  un = cl->un;
  f += isnd[snde];
  l += isnd[snde];

  for (i = f; i < l; ++i)
  {
    x[i] /= diag[i];
    xi = x[i];
    ls = ju + iju[i];
    l1 = un + iu[i];
    sze = l - i - 1;

    for (t = 0; t < sze; ++t)
      x[ls[t]] -= l1[t] * xi;
  }
} /* subsFwrd */

static void subsBwrd(int nrow,
                     double *diag,
                     double *un,
                     int *fir,
                     double *x)
{
  int i, t, sze;
  double x1, x2, rtemp, *x0, *l1, *l2;

  for (i = nrow; i;)
  {
    for (; i > 1; --i)
    {
      --i;
      l1 = un + fir[i - 1] + 1;
      l2 = un + fir[i] + 0;
      sze = nrow - i - 1;
      x1 = 0.0;
      x2 = 0.0;
      x0 = x + 1 + i;

      for (t = 0; t < sze; ++t)
      {
        rtemp = x0[t];

        x1 += l1[t] * rtemp;
        x2 += l2[t] * rtemp;
      }

      x[i] -= x2 / diag[i];
      x[i - 1] -= (un[fir[i - 1]] * x[i] + x1) / diag[i - 1];
    }

    for (; i;)
    {
      --i;
      l1 = un + fir[i];
      sze = nrow - i - 1;
      x1 = 0.0;
      x0 = x + 1 + i;

      for (t = 0; t < sze; ++t)
        x1 += l1[t] * x0[t];

      x[i] -= x1 / diag[i];
    }
  }
} /* subsBwrd */

void CholSol(chol *cl,
             double *b,
             double *x)
{
  int i, k, s, t, sze, f, l, itemp, *ls, *isnd, *jun, *ju, *iju, *iu;
  double x1, x2, rtemp1, rtemp2, rtemp3, rtemp4, rtemp5, rtemp6,
      rtemp7, rtemp8, *l1, *l3, *l2, *l4, *l5, *l6, *l7, *l8,
      *diag, *un;

  isnd = cl->isnd;
  jun = cl->jun;
  ju = cl->ju;
  iju = cl->iju;
  iu = cl->iu;
  diag = cl->ud;
  un = cl->un;

  for (i = 0; i < cl->nrow; ++i)
    x[i] = b[cl->perm[i]];

  for (s = 0; s < cl->nsnd; ++s)
  {
    f = isnd[s];
    l = isnd[s + 1];

    subsFwrd(cl, s, 0, l - f, x);

    itemp = l - f - 1;
    ls = ju + iju[f] + itemp;
    sze = jun[f] - itemp;
    k = f;

    itemp = l - 1;
    for (; k + 7 < l; k += 8)
    {
      l1 = un + iu[k + 0] + itemp - (k + 0);
      l2 = un + iu[k + 1] + itemp - (k + 1);
      l3 = un + iu[k + 2] + itemp - (k + 2);
      l4 = un + iu[k + 3] + itemp - (k + 3);
      l5 = un + iu[k + 4] + itemp - (k + 4);
      l6 = un + iu[k + 5] + itemp - (k + 5);
      l7 = un + iu[k + 6] + itemp - (k + 6);
      l8 = un + iu[k + 7] + itemp - (k + 7);

      rtemp1 = x[k + 0];
      rtemp2 = x[k + 1];
      rtemp3 = x[k + 2];
      rtemp4 = x[k + 3];
      rtemp5 = x[k + 4];
      rtemp6 = x[k + 5];
      rtemp7 = x[k + 6];
      rtemp8 = x[k + 7];

      for (t = 0; t < sze; ++t)
        x[ls[t]] -= rtemp1 * l1[t] + rtemp2 * l2[t] + rtemp3 * l3[t] + rtemp4 * l4[t] + rtemp5 * l5[t] + rtemp6 * l6[t] + rtemp7 * l7[t] + rtemp8 * l8[t];
    }

    for (; k + 3 < l; k += 4)
    {
      l1 = un + iu[k + 0] + itemp - (k + 0);
      l2 = un + iu[k + 1] + itemp - (k + 1);
      l3 = un + iu[k + 2] + itemp - (k + 2);
      l4 = un + iu[k + 3] + itemp - (k + 3);

      rtemp1 = x[k + 0];
      rtemp2 = x[k + 1];
      rtemp3 = x[k + 2];
      rtemp4 = x[k + 3];

      for (t = 0; t < sze; ++t)
        x[ls[t]] -= rtemp1 * l1[t] + rtemp2 * l2[t] + rtemp3 * l3[t] + rtemp4 * l4[t];
    }

    for (; k + 1 < l; k += 2)
    {
      l1 = un + iu[k + 0] + itemp - (k + 0);
      l2 = un + iu[k + 1] + itemp - (k + 1);
      rtemp1 = x[k + 0];
      rtemp2 = x[k + 1];
      for (t = 0; t < sze; ++t)
        x[ls[t]] -= rtemp1 * l1[t] + rtemp2 * l2[t];
    }

    for (; k < l; ++k)
    {
      l1 = un + iu[k + 0] + itemp - (k + 0);
      rtemp1 = x[k + 0];
      for (t = 0; t < sze; ++t)
        x[ls[t]] -= rtemp1 * l1[t];
    }
  }

  if (cl->nsnd)
  {
    s = cl->nsnd - 1;
    f = isnd[s];
    l = isnd[s + 1];

    dCopy(l - f, x + f, b + f);

    subsBwrd(l - f, diag + f, un, iu + f, b + f);

    s = cl->nsnd - 1;

    for (; s >= 1; --s)
    {
      f = isnd[s - 1];
      l = isnd[s];
      i = l;

      for (; i > 1 + f; --i)
      {
        --i;
        ls = ju + iju[i];
        l1 = un + iu[i - 1] + 1;
        l2 = un + iu[i] + 0;
        sze = jun[i];
        x1 = 0.0;
        x2 = 0.0;

        for (t = 0; t < sze; ++t)
        {
          rtemp1 = b[ls[t]];
          x1 += l1[t] * rtemp1;
          x2 += l2[t] * rtemp1;
        }

        b[i] = x[i] - x2 / diag[i];
        b[i - 1] = x[i - 1] - (x1 + un[iu[i - 1]] * b[i]) / diag[i - 1];
      }

      for (; i > f;)
      {
        --i;
        l1 = un + iu[i];
        ls = ju + iju[i];
        sze = jun[i];
        x1 = 0.0;

        for (t = 0; t < sze; ++t)
          x1 += l1[t] * b[ls[t]];

        b[i] = x[i] - x1 / diag[i];
      }
    }
  }

  for (i = 0; i < cl->nrow; ++i)
    x[i] = b[cl->invp[i]];
} /* CholSol */

void SolArgLES(chol *cl,
               double *d,
               matrix *at,
               double *b,
               double *x)
{
  int j, nrow, ncol;

  nrow = cl->nrow;
  ncol = at->nrow;

  for (j = 0; j < ncol; ++j)
    x[j] = b[j] / d[j];

  if (at)
    mTimesv(false, nrow, ncol, -1.0, at, x, 1.0, b + ncol);

  dCopy(ncol, b, x);

  CholSol(cl, b + ncol, x + ncol);

  if (at)
    mTimesv(true, nrow, ncol, 1.0, at, x + ncol, 1.0, x);

  for (j = 0; j < nrow; ++j)
    x[ncol + j] = -x[ncol + j];

  for (j = 0; j < ncol; ++j)
    x[j] /= d[j];
} /* SolArgLES */
