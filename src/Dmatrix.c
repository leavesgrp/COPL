#include "LPmatrx.h"

static void InitMtx(matrix *a)
{
  int i;

  for (i = 0; i < a->nrow; i++)
    a->ia[i].nn0 = 0;
} /* InitMtx */

static void ScanMtx(matrix *a,
                    matrix *at)
{
  int i, j, k, nrow = a->nrow, nn0;

  for (i = 0; i < nrow; i++)
  {
    nn0 = a->ia[i].nn0;
    for (k = 0; k < nn0; k++)
    {
      j = a->ia[i].ja[k];
      at->ia[j].nn0++;
    }
  }
} /* ScanMtx */

static void makeMtxrw(matrix *a,
                      int mems,
                      char *info)
{
  if (a->ia)
  {
    iFree(&a->ia->ja);
    dFree(&a->ia->an);
    a->ia->ja = iAlloc(mems, info);
    a->ia->an = dAlloc(mems, info);
  }
} /* makeMtxrw */

int MtxTrans(int ncol,
             matrix *a,
             int pid,
             int *p,
             int vid,
             matrix **at)
{
  int i, j, k, t, nn0;
  matrix *r;
  array *ai, *aj;

  nn0 = 0;
  for (i = 0; i < a->nrow; i++)
    nn0 += a->ia[i].nn0;

  r = *at;

  if (r)
  {
    if (!vid && nn0 > r->mems)
      return false;

    if (ncol > r->mnrs)
      return false;

    r->nrow = ncol;
  }

  else
  {
    r = MtxAlloc(ncol, 0, "r, MtxTrans");
    makeMtxrw(r, nn0, "r, MtxTrans");
  }

  if (vid)
    for (j = 0; j < ncol; j++)
      r->ia[j].nn0 = 0;

  else
  {
    InitMtx(r);
    ScanMtx(a, r);

    nn0 = 0;
    for (j = 0; j < ncol; j++)
    {
      r->ia[j].ja = r->ia->ja + nn0;
      r->ia[j].an = r->ia->an + nn0;
      nn0 += r->ia[j].nn0;
      r->ia[j].nn0 = 0;
    }

    r->mems = nn0;
    if (nn0 > r->mems)
    {
      printf("\n\n Exit -- 140: out of internal space.\n\n");
      exit(140);
    }
  }

  for (k = 0; k < a->nrow; k++)
  {
    i = k;
    if (p)
    {
      ai = a->ia + p[k];
      if (!pid)
        i = p[k];
    }
    else
      ai = a->ia + k;

    for (t = 0; t < ai->nn0; t++)
    {
      j = ai->ja[t];
      aj = r->ia + j;
      nn0 = aj->nn0;

      aj->ja[nn0] = i;
      aj->an[nn0] = ai->an[t];
      aj->nn0++;
    }
  }

  *at = r;

  return true;
} /* MtxTrans */

static int SetChol(chol *cl,
                   int *nnzi)
{
  int i, annz, nrow;

  nrow = cl->nrow;
  annz = 0;
  for (i = 0; i < nrow; i++)
    annz += nnzi[i];

  if (!annz)
    return true;

  cl->jaat = iAlloc(annz, "cl->jaat, SetSfact");
  cl->maxn = annz;
#ifdef TEST
  fprintf(fres, "&%8d", nrow + 2 * annz);
#endif
  cl->iaat[0] = 0;
  cl->naat[0] = 0;
  for (i = 1; i < nrow; i++)
  {
    cl->iaat[i] = cl->iaat[i - 1] + nnzi[i - 1];
    cl->naat[i] = 0;
  }

  return true;
} /* SetChol */

static int ChkSub(chol *cl,
                  int i,
                  int j)
{
  int r, c;

  c = min(i, j);
  r = max(i, j);

  if (c < 0 || r >= cl->nrow)
    return false;

  if (r == c)
    return true;

  cl->jaat[cl->iaat[c] + cl->naat[c]] = r;
  cl->naat[c]++;

  return true;
} /* ChkSub */

static void plusXs(int n,
                   int *x,
                   int *s,
                   int *m)
{
  int i;

  if (!s && !m)
  {
    for (i = 0; i < n; i++)
      x[i]++;
  }
  else if (s && !m)
  {
    for (i = 0; i < n; i++)
      x[s[i]]++;
  }
  else
    exit(0);
} /* plusXs */

static void addRow(int i,
                   int ifir,
                   int *ihead,
                   int *ilink)
{
  int temp;

  temp = ihead[ifir];
  ihead[ifir] = i;
  ilink[i] = temp;
} /* addRow */

static void formSymtx(int nrow,
                      int *fir,
                      int *nnz,
                      int *sub,
                      int *p,
                      int rwws,
                      int *firt,
                      int *nnzt,
                      int *subt)
{
  int i, j, s, t, stopt;

  iZero(nrow, nnzt, NULL);

  if (rwws)
  {
    if (p)
    {
      for (s = 0; s < nrow; ++s)
      {
        j = p[s];
        for (t = fir[s], stopt = t + nnz[s]; t < stopt; ++t)
        {
          i = p[sub[t]];
          nnzt[max(i, j)]++;
        }
      }
    }
    else
    {
      for (j = 0; j < nrow; j++)
      {
        for (t = fir[j], stopt = t + nnz[j]; t < stopt; ++t)
        {
          i = sub[t];
          nnzt[max(i, j)]++;
        }
      }
    }
  }

  else
  {
    if (p)
    {
      for (s = 0; s < nrow; ++s)
      {
        j = p[s];
        for (t = fir[s], stopt = t + nnz[s]; t < stopt; ++t)
        {
          i = p[sub[t]];
          nnzt[min(i, j)]++;
        }
      }
    }
    else
    {
      for (j = 0; j < nrow; j++)
      {
        for (t = fir[j], stopt = t + nnz[j]; t < stopt; ++t)
        {
          i = sub[t];
          nnzt[min(i, j)]++;
        }
      }
    }
  }

  firt[0] = 0;
  for (i = 1; i < nrow; i++)
  {
    firt[i] = firt[i - 1] + nnzt[i - 1];
    nnzt[i - 1] = 0;
  }
  nnzt[nrow - 1] = 0;

  if (rwws)
  {
    if (p)
    {
      for (s = 0; s < nrow; ++s)
      {
        j = p[s];
        for (t = fir[s], stopt = t + nnz[s]; t < stopt; ++t)
        {
          i = p[sub[t]];
          if (i > j)
          {
            subt[firt[i] + nnzt[i]] = j;
            nnzt[i]++;
          }
          else
          {
            subt[firt[j] + nnzt[j]] = i;
            nnzt[j]++;
          }
        }
      }
    }
    else
    {
      for (j = 0; j < nrow; j++)
      {
        for (t = fir[j], stopt = t + nnz[j]; t < stopt; ++t)
        {
          i = sub[t];
          if (i > j)
          {
            subt[firt[i] + nnzt[i]] = j;
            nnzt[i]++;
          }
          else
          {
            subt[firt[j] + nnzt[j]] = i;
            nnzt[j]++;
          }
        }
      }
    }
  }

  else
  {
    if (p)
    {
      for (s = 0; s < nrow; ++s)
      {
        j = p[s];
        for (t = fir[s], stopt = t + nnz[s]; t < stopt; ++t)
        {
          i = p[sub[t]];
          if (i < j)
          {
            subt[firt[i] + nnzt[i]] = j;
            nnzt[i]++;
          }
          else
          {
            subt[firt[j] + nnzt[j]] = i;
            nnzt[j]++;
          }
        }
      }
    }
    else
    {
      for (j = 0; j < nrow; j++)
      {
        for (t = fir[j], stopt = t + nnz[j]; t < stopt; ++t)
        {
          i = sub[t];
          if (i < j)
          {
            subt[firt[i] + nnzt[i]] = j;
            nnzt[i]++;
          }
          else
          {
            subt[firt[j] + nnzt[j]] = i;
            nnzt[j]++;
          }
        }
      }
    }
  }
} /* formSymtx */

static void getDProw(int nrow,
                     int ncol,
                     int *fir,
                     int *sze,
                     int *sub,
                     int *map,
                     int *ihead,
                     int *ilink,
                     int *ilist,
                     int *irec,
                     int *iptr1)
{
  int i, new, n, oisze, isze, s, count, temp, k, nexti, *cur;

  n = nrow;
  *irec = 0;
  cur = iptr1;

  for (i = 0; i < nrow; i++)
  {
    cur[i] = 0;
    ilink[i] = n;
    ilist[i] = n;
  }
  for (i = 0; i < ncol; i++)
    ihead[i] = n;

  isze = 0;
  count = 0;
  oisze = isze;
  new = n;
  for (i = 0; i < nrow; i++)
  {
    if (map)
      for (; cur[i] < sze[i] && !map[sub[fir[i] + cur[i]]]; ++cur[i])
        ;

    if (cur[i] < sze[i])
    {
      s = sub[fir[i] + cur[i]];
      if (ihead[s] == n)
        ilist[isze++] = s;

      addRow(i, s, ihead, ilink);

      cur[i]++;
    }

    else
    {
      temp = new;
      new = i;
      ilink[i] = temp;
    }
  }

  for (k = oisze; k < isze; k++)
  {
    temp = ihead[ilist[k]];
    ihead[ilist[k]] = n;
    ilist[k] = temp;
  }

  if (new != n)
  {
    count++;
    ilist[nrow - count] = new;
  }

  while (isze)
  {
    isze--;
    oisze = isze;

    i = ilist[isze];
    ilist[isze] = n;

    if (i == n)
      exit(0);

    new = n;
    if (ilink[i] == n)
      new = i;
    else
    {
      for (; i != n; i = nexti)
      {
        nexti = ilink[i];
        ilink[i] = n;

        if (map)
          for (; cur[i] < sze[i] && !map[sub[fir[i] + cur[i]]]; ++cur[i])
            ;

        if (cur[i] < sze[i])
        {
          s = sub[fir[i] + cur[i]];
          cur[i]++;

          if (ihead[s] == n)
            ilist[isze++] = s;

          temp = ihead[s];
          ihead[s] = i;
          ilink[i] = temp;
        }

        else
        {
          temp = new;
          new = i;
          ilink[i] = temp;
        }
      }
    }

    for (k = oisze; k < isze; k++)
    {
      temp = ihead[ilist[k]];
      ihead[ilist[k]] = n;
      ilist[k] = temp;
    }

    if (new != n)
    {
      count++;
      ilist[nrow - count] = new;
    }
  }

  *irec = count;
  for (k = 0; k < count; k++)
    ilist[k] = ilist[nrow - count + k];
} /* getDProw */

static int iComp(const void *e1,
                 const void *e2)
{
  int *i1, *i2;

  i1 = (int *)e1;
  i2 = (int *)e2;

  if (*i1 < *i2)
    return (-1);
  else if (*i1 > *i2)
    return (1);
  return (0);
} /* iComp */

static void iSort(int n,
                  int *x)
{
  qsort((void *)x, n, sizeof(int), iComp);
} /* iSort */

static void chkDense(chol *cl,
                     int *iptr1,
                     int *iptr2,
                     int *iptr3,
                     int *iptr4,
                     int *iptr5,
                     int *iptr6)
{
  int j, k, l, n, t, ndnd, *isnd, *iju, *jun,
      *ju, *fir, *sze, *ilist, *ilink;

  n = cl->nrow;
  isnd = cl->isnd;
  iju = cl->iju;
  jun = cl->jun;
  ju = cl->ju;

  if (!cl->nsnd ||
      !iptr1 | !iptr2 || !iptr3 ||
      !iptr4 || !iptr5 || !iptr6)
  {
    cl->dnsp = false;
    return;
  }

  cl->dnsp = true;
  fir = iptr1;
  sze = iptr2;
  ilist = iptr3;
  ilink = iptr4;

  cl->ndns = 0;

  l = isnd[cl->nsnd - 1];
  for (k = 0; k + 1 < cl->nsnd; k++)
  {
    j = isnd[k];
    for (t = 0; t < jun[j] && ju[iju[j] + t] < l; ++t)
      ;

    fir[k] = iju[j] + t;
    sze[k] = jun[j] - t;
  }

  getDProw(cl->nsnd - 1, cl->nrow, fir, sze, ju,
           NULL, iptr6, ilink, ilist, &ndnd, iptr5);

  cl->idn = iAlloc(ndnd + 1, "cl->idn, chkDense");
  cl->jdn = iAlloc(cl->nsnd, "cl->jdn, chkDense");
  cl->dmj = iAlloc(cl->nsnd, "cl->dmj, chkDense");

  n = cl->nsnd - 1;
  cl->ndnd = 0;
  cl->ndns = 0;
  cl->idn[0] = 0;

  for (k = 0; k < ndnd; k++)
  {
    j = ilist[k];
    if (sze[j])
    {
      cl->idn[cl->ndnd + 1] = cl->idn[cl->ndnd];
      cl->ndnd++;
      for (; j != n; j = ilink[j])
      {
        cl->ndns += cl->isnd[j + 1] - cl->isnd[j];
        cl->jdn[cl->idn[cl->ndnd]] = j;
        cl->dmj[cl->idn[cl->ndnd]] = fir[j] - iju[isnd[j]];
        cl->idn[cl->ndnd]++;
      }
      iSort(cl->idn[cl->ndnd] - cl->idn[cl->ndnd - 1],
            cl->jdn + cl->idn[cl->ndnd - 1]);

      for (t = cl->idn[cl->ndnd - 1]; t < cl->idn[cl->ndnd]; ++t)
        cl->dmj[t] = fir[cl->jdn[t]] - iju[isnd[cl->jdn[t]]];
    }
  }
} /* chkDense */

static int SymboChol(chol *cl,
                     int usz)
{
  int chksn, i, j, t, stopt, sze, first, cur, k, paft, ipos, nrow, iend,
      *nnz, *fir, *pja, *link, *buf, *mask, *ju, *tju, *iptr1, *iptr2,
      *iptr3, *iptr4, *p, *invp, *iju, *jun, *isnd;

  paft = 0;
  nrow = cl->nrow;
  iend = nrow;
  p = cl->perm;
  invp = cl->invp;
  iju = cl->iju;
  jun = cl->jun;
  isnd = cl->isnd;

  pja = iAlloc(cl->maxn, "jat, SymboChol");

  for (i = 0; i < nrow; i++)
    invp[p[i]] = i;

  nnz = cl->iu;
  fir = cl->isnd;

  formSymtx(nrow, cl->iaat, cl->naat, cl->jaat, invp, true, fir, nnz, pja);
  formSymtx(nrow, fir, nnz, pja, NULL, false, cl->iaat, cl->naat, cl->jaat);

  iFree(&pja);

  k = usz + nrow;
  ju = iAlloc(k, "ju, SymboChol");
  buf = ju + usz;

  mask = cl->iu;

  link = invp;

  for (i = 0; i < nrow; i++)
  {
    mask[i] = 0;
    link[i] = iend;
  }

  paft = 0;
  cl->nsnd = 0;
  isnd[0] = 0;
  for (i = 0; i < nrow; i++)
  {
    sze = cl->naat[i];
    first = iend;
    cur = link[i];
    chksn = false;

    if (cur == iend)
    {

      isnd[cl->nsnd + 1] = 1 + isnd[cl->nsnd];
      jun[i] = sze;
      iju[i] = paft;
      paft += sze;

      iCopy(sze, cl->jaat + cl->iaat[i], ju + iju[i]);
      if (sze)
      {
        first = ju[iju[i]];
        for (cur = first; link[cur] != iend; cur = link[cur])
          ;
        link[cur] = cl->nsnd;
        link[cl->nsnd] = iend;
      }
      cl->nsnd++;
    }

    else
    {
      mask[i] = 1;

      iCopy(sze, cl->jaat + cl->iaat[i], buf);
      iSet(sze, 1, mask, buf);

      for (; cur != iend; cur = link[cur])
      {
        chksn |= (1 + cur == cl->nsnd);
        k = isnd[cur];

        for (t = iju[k], stopt = t + jun[k]; t < stopt; ++t)
        {
          j = ju[t];
          if (j > i && !mask[j])
          {
            buf[sze] = j;
            mask[j] = 1;
            sze++;
          }
        }
      }

      if (chksn)
      {
        k = isnd[cl->nsnd - 1];
        chksn = sze == (jun[k] - (isnd[cl->nsnd] - isnd[cl->nsnd - 1]));
      }

      first = nrow;
      mask[i] = 0;
      for (t = 0; t < sze; ++t)
      {
        j = buf[t];
        mask[j] = 0;
        first = min(j, first);
      }

      if (chksn)
      {
        ipos = getPos(jun[i - 1], i, ju + iju[i - 1]);

        if (ipos == jun[i - 1])
        {
          printf("\n\n Exit -- 133: position error.\n\n");
          exit(0);
        }

        iSwap(iju[i - 1], ipos + iju[i - 1], ju);

        isnd[cl->nsnd]++;
        iju[i] = iju[i - 1] + 1;
        jun[i] = jun[i - 1] - 1;

        if (ju[iju[i] - 1] != i)
        {
          printf("\n\n Exit -- 143: link error lk=" IFMT " li" IFMT "\n\n",
                 iju[k], iju[i]);
          exit(0);
        }

        if (first != iend)
        {
          for (cur = first; link[cur] != iend; cur = link[cur])
            ;
          link[cur] = cl->nsnd - 1;
          link[cl->nsnd - 1] = iend;
        }
      }

      else
      {
        isnd[cl->nsnd + 1] = 1 + isnd[cl->nsnd];
        iju[i] = paft;
        jun[i] = sze;
        paft += sze;

        if (paft > usz)
        {
          printf("\n\n Exit -- 153: space prob paft=" IFMT ".\n\n", paft);
          exit(0);
        }

        iCopy(sze, buf, ju + iju[i]);

        if (first != iend)
        {
          for (cur = first; link[cur] != iend; cur = link[cur])
            ;
          link[cur] = cl->nsnd;
          link[cl->nsnd] = iend;
        }
        cl->nsnd++;
      }
    }

    if (jun[i] + 1 == nrow - i)
      break;
  }

  for (i++; i < nrow; i++)
  {
    jun[i] = jun[i - 1] - 1;
    iju[i] = iju[i - 1] + 1;

    isnd[cl->nsnd]++;
  }

  tju = iAlloc(paft, "tju, SymboChol");

  fir = buf;
  nnz = cl->iu;

  iZero(nrow, nnz, NULL);

  for (k = 0; k < cl->nsnd; k++)
  {
    j = isnd[k];
    plusXs(jun[j], nnz, ju + iju[j], NULL);
  }

  fir[0] = 0;
  for (k = 1; k < nrow; k++)
    fir[k] = fir[k - 1] + nnz[k - 1];

  iZero(nrow, nnz, NULL);

  for (k = 0; k < cl->nsnd; k++)
  {
    j = isnd[k];
    for (t = iju[j], stopt = t + jun[j]; t < stopt; ++t)
    {
      i = ju[t];
      tju[fir[i] + nnz[i]] = j;
      nnz[i]++;
    }
    jun[j] = 0;
  }

  for (i = 0; i < nrow; i++)
  {
    for (t = fir[i], stopt = t + nnz[i]; t < stopt; ++t)
    {
      j = tju[t];
      ju[iju[j] + jun[j]] = i;
      jun[j]++;
    }
  }

  iFree(&tju);

  if (paft <= cl->nju)
  {
    iCopy(paft, ju, cl->ju);
    iFree(&ju);
  }

  else
  {
    cl->nju = 0;
    iFree(&cl->ju);

    cl->ju = iAlloc(paft, "cl->ju, SymboChol");
    iCopy(paft, ju, cl->ju);

    cl->nju = paft;
    iFree(&ju);
  }

  iptr1 = iAlloc(4 * nrow, "iptr1, SymboChol");
  iptr2 = iptr1 + nrow;
  iptr3 = iptr2 + nrow;
  iptr4 = iptr3 + nrow;

  chkDense(cl, cl->iu, cl->invp, iptr1, iptr2, iptr3, iptr4);

  iFree(&iptr1);

  cl->iu[0] = 0;
  for (i = 1; i < nrow; i++)
    cl->iu[i] = cl->iu[i - 1] + cl->jun[i - 1];

  for (i = 0; i < nrow; i++)
    invp[p[i]] = i;

  for (k = 0; k < cl->nsnd; k++)
    if (isnd[k] + 1 != isnd[k + 1])
      break;

  cl->plmt = k;

  return true;
} /* SymboChol */

static int makeClval(chol *cl,
                     char *info)
{
  int nnz;

  nnz = iSum(cl->nrow, cl->jun);

  if (nnz <= cl->nu)
    return true;

  cl->nu = 0;
  if (cl->un)
    dFree(&cl->un);
  cl->un = dAlloc(nnz, info);

  cl->nu = nnz;

  return true;
} /* makeClval */

static int SymboCfac(chol *cl)
{
  int i, t, lnnz, *nnzi, nrow, ista;
  order *od;

  nrow = cl->nrow;

  od = OdAlloc(nrow, 2 * cl->maxn, "od, SymboCfac");

  nnzi = cl->perm;
  iZero(nrow, nnzi, NULL);

  for (i = 0; i < nrow; i++)
  {
    nnzi[i] += cl->naat[i];
    plusXs(cl->naat[i], nnzi, cl->jaat + cl->iaat[i], NULL);
  }

  OdInit(od, nnzi);
  for (i = 0; i < nrow; i++)
    for (t = 0; t < cl->naat[i]; ++t)
      OdIndex(od, i, cl->jaat[cl->iaat[i] + t]);

  GetOrder(od, cl->perm);

#ifdef TEST
  fprintf(fres, "&%8d", od->ntot + nrow);
#endif

  lnnz = od->ntot;
  OdFree(&od);

  ista = SymboChol(cl, lnnz);

  return ista;
} /* SymboCfac */

int SymboProc(chol *cl,
              matrix *a,
              matrix *at)
{
  int i, j, sze, k, t, nnzaj, next, nrow, ncol, iend,
      *bmap, *imem, *jfir, *jcur, *jnex, *nnzi, *subaj, *p;
  array *aj;

  nrow = cl->nrow;
  ncol = at->nrow;
  iend = ncol;

  MtxTrans(nrow, at, false, NULL, false, &a);
  MtxTrans(ncol, a, false, NULL, true, &at);

  sze = 2 * ncol + 4 * nrow;
  ot->imem = iAlloc(sze, "ot->imem, StructKfact");

  jfir = ot->imem;
  jcur = jfir + nrow;
  jnex = jcur + ncol;
  imem = jnex + ncol;
  bmap = imem + nrow;
  nnzi = bmap + nrow;

  iSet(nrow, iend, jfir, NULL);

  for (j = 0; j < ncol; j++)
  {
    jcur[j] = 0;
    aj = at->ia + j;
    nnzaj = aj->nn0;
    subaj = aj->ja;

    if (nnzaj)
    {
      i = subaj[0];
      jnex[j] = jfir[i];
      jfir[i] = j;
    }
    else
      jnex[j] = iend;
  }

  iZero(nrow, bmap, NULL);

  for (i = 0; i < nrow; i++)
  {
    bmap[i] = true;
    sze = 0;

    for (j = jfir[i]; j != iend; j = next)
    {

      next = jnex[j];
      aj = at->ia + j;
      subaj = aj->ja;

      for (k = jcur[j]; k < aj->nn0; k++)
      {
        t = subaj[k];
        if (!bmap[t])
        {
          imem[sze++] = t;
          bmap[t] = true;
        }
      }

      jcur[j]++;

      if (jcur[j] < aj->nn0)
      {
        jnex[j] = jfir[aj->ja[jcur[j]]];
        jfir[aj->ja[jcur[j]]] = j;
      }
    }

    nnzi[i] = sze;
    bmap[i] = false;
    iZero(sze, bmap, imem);
  }

  SetChol(cl, nnzi);
  iSet(nrow, iend, jfir, NULL);

  for (j = 0; j < ncol; j++)
  {
    jcur[j] = 0;
    aj = at->ia + j;
    nnzaj = aj->nn0;
    subaj = aj->ja;

    if (nnzaj)
    {
      i = subaj[0];
      jnex[j] = jfir[i];
      jfir[i] = j;
    }
    else
      jnex[j] = iend;
  }

  iZero(nrow, bmap, NULL);

  for (i = 0; i < nrow; i++)
  {
    bmap[i] = true;
    sze = 0;
    for (j = jfir[i]; j != iend; j = next)
    {
      next = jnex[j];
      aj = at->ia + j;
      subaj = aj->ja;

      for (k = jcur[j]; k < aj->nn0; k++)
      {
        t = subaj[k];
        if (!bmap[t])
        {
          imem[sze++] = t;
          bmap[t] = true;
        }
      }

      jcur[j]++;

      if (jcur[j] < aj->nn0)
      {
        jnex[j] = jfir[aj->ja[jcur[j]]];
        jfir[aj->ja[jcur[j]]] = j;
      }
    }

    for (t = 0; t < sze; ++t)
      if (!ChkSub(cl, i, imem[t]))
        exit(0);

    bmap[i] = false;
    iZero(sze, bmap, imem);
  }

  if (!SymboCfac(cl))
    exit(0);

  p = cl->perm;
  MtxTrans(ncol, a, false, p, true, &at);

  MtxFree(&a);
  iFree(&ot->imem);

  makeClval(cl, "cl, StructKfact");

  return true;
} /* SymboProc */

static void ScalVect(int n,
                     double alfa,
                     double *x)
{
  int i;

  if (n <= 0 || alfa == 1.0)
    return;

  if (alfa == 0.0)
  {
    memset(x, 0, sizeof(double) * n);
    return;
  }

  else if (alfa == -1.0)
  {
    for (i = 0; i < n; i++)
      x[i] = -x[i];
  }

  else
  {
    for (i = 0; i < n; i++)
      x[i] *= alfa;
  }

} /* ScalVect */

void mTimesv(int tran,
             int nrow,
             int ncol,
             double alfa,
             matrix *a,
             double *x,
             double beta,
             double *y)
{
  int xsz, ysz, i;

  if (nrow < 0 || ncol < 0)
  {
    ShutDown();
    exit(0);
  }

  if (!tran)
  {
    xsz = ncol;
    ysz = nrow;

    ScalVect(ysz, beta, y);

    if (alfa == 0.0)
      return;

    for (i = 0; i < xsz; i++)
      if (x[i])
        setArray(alfa * x[i], a->ia + i, y);
  }
  else
  {
    xsz = nrow;
    ysz = ncol;

    ScalVect(ysz, beta, y);

    if (alfa == 0.0)
      return;

    for (i = 0; i < ysz; i++)
      y[i] += alfa * svDot(a->ia + i, x);
  }
} /* mTimesv */
