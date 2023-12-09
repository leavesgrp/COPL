#include "LPmatrx.h"

void OdInit(order *od, int *nnzi)
{
  int i, n = od->nnod;

  if (n)
  {
    od->rexs[0] = nnzi[0];
    od->rlen[0] = nnzi[0];
    od->rbeg[0] = 0;
    od->pres[0] = n;
    od->succ[0] = 1;
    for (i = 1; i < od->nnod; ++i)
    {
      od->pres[i] = i - 1;
      od->succ[i] = i + 1;
      od->rexs[i] = nnzi[i];
      od->rlen[i] = nnzi[i];
      od->rbeg[i] = od->rbeg[i - 1] + od->rlen[i - 1];
    }

    od->succ[n - 1] = n;
    od->last = n - 1;

    od->raft = od->rbeg[n - 1] + nnzi[n - 1];

    if (od->raft > od->nn0)
    {
      ErrorProc(OutOfSpc, "InitMmd");
      ShutDown();
      exit(0);
    }
  }
} /* OdInit */

static void OdSet(order *od,
                  int allow_eli,
                  xlist *elist,
                  int *idrw,
                  int *marker,
                  int *isize,
                  int *ilink,
                  int *oinfo,
                  int *osize,
                  int *e,
                  int *p)
{
  int i, n, deg, *rbeg, *rexs, *rend;

  n = od->nnod;
  rbeg = od->rbeg;
  rexs = od->rexs;
  rend = od->rend;
  *e = 0;

  for (i = 0; i < n; i++)
  {
    isize[i] = 0;
    ilink[i] = n;
    osize[i] = 0;
    oinfo[i] = n;
    marker[i] = 0;
  }

  for (i = 0; i < n; ++i)
  {
    rbeg[i] -= rexs[i];
    rend[i] = 0;
  }

  for (i = 0; i < n; ++i)
  {
    deg = rexs[i];
    if (!allow_eli || deg)
    {
      idrw[i] = 1;
      putXt(elist, i, deg);
    }
    else
    {
      idrw[i] = 0;
      marker[i] = true;
      p[*e] = i;
      (*e)++;
    }
  }
} /* OdSet */

void OdIndex(order *od,
             int i,
             int j)
{
  if (i != j)
  {
    od->adjn[od->rbeg[i]++] = j;
    od->adjn[od->rbeg[j]++] = i;
  }
} /* OdIndex */

static void OdArriv(order *od,
                    int *idrw,
                    int *marker,
                    int *isize,
                    int x,
                    int *xdeg,
                    int *rsze,
                    int *esze,
                    int *rcrw)
{
  int *visited, i, n, y, z, l, s, t, f, stopt,
      stops, *adjn, *rbeg, *rexs, *rend;

  n = od->nnod;
  adjn = od->adjn;
  rbeg = od->rbeg;
  rexs = od->rexs;
  rend = od->rend;
  *rsze = 0;
  *esze = 0;

  if (rexs[x])
  {
    l = n;

    visited = marker;
    visited[x] = true;

    for (t = rbeg[x], stopt = rbeg[x] + rend[x]; t < stopt; ++t)
    {
      y = adjn[t];

      if (idrw[y] != 0)
      {
        l--;
        rcrw[l] = y;
        visited[y] = true;

        for (s = rbeg[y], stops = rbeg[y] + rexs[y]; s < stops; ++s)
        {
          z = adjn[s];

          if (idrw[z] != 0)
          {
            if (!visited[z])
            {
              visited[z] = true;

              rcrw[*rsze] = z;
              (*rsze)++;
            }
          }
        }
      }
    }

    f = rbeg[x] + rend[x];
    for (t = f, stopt = rbeg[x] + rexs[x]; t < stopt; ++t)
    {
      y = adjn[t];
      if (!visited[y])
      {
        adjn[f++] = y;
        visited[y] = true;

        rcrw[*rsze] = y;
        (*rsze)++;
      }
    }

    rexs[x] = f - rbeg[x];

    *esze = n - l;
    visited[x] = false;
    iZero(*rsze, visited, rcrw);
    iZero(n - l, visited, rcrw + l);
  }

  if (xdeg)
  {
    *xdeg = *rsze + isize[x];
    for (i = 0; i < *rsze; ++i)
      *xdeg += isize[rcrw[i]];
  }
} /* OdArriv */

static void OdRenew(order *od,
                    int *ilink,
                    int x,
                    int xdeg,
                    int *e,
                    int *p)
{
  int c, n;

  n = od->nnod;
  od->ntot += xdeg--;
  p[*e] = x;
  (*e)++;
  for (c = x; ilink[c] != n; c = ilink[c])
  {
    od->ntot += xdeg--;
    p[*e] = ilink[c];
    (*e)++;
  }
} /* OdRenew */

static void OdCheck(order *od,
                    int *idrw)
{
  int f, i, t, stopt, new, z, previous, n, *adjn,
      *rbeg, *rexs, *rlen, *rend, *pres, *succ;

  n = od->nnod;
  adjn = od->adjn;
  rbeg = od->rbeg;
  rexs = od->rexs;
  rlen = od->rlen;
  rend = od->rend;
  pres = od->pres;
  succ = od->succ;

  f = 0;
  previous = n;
  for (i = od->head; i != n; i = succ[i])
  {
    if (idrw[i] != 0)
    {
      new = f;

      for (t = rbeg[i], stopt = rbeg[i] + rend[i]; t < stopt; ++t)
      {
        z = adjn[t];
        if (idrw[z] == 3)
          adjn[f++] = z;
      }

      rend[i] = f - new;

      for (stopt = rbeg[i] + rexs[i]; t < stopt; ++t)
      {
        z = adjn[t];
        if (idrw[z] != 0)
          adjn[f++] = z;
      }

      rexs[i] = f - new;
      rlen[i] = rexs[i];

      rbeg[i] = new;

      if (previous == n)
      {
        od->head = i;
        pres[i] = n;
      }

      else
      {
        succ[previous] = i;
        pres[i] = previous;
      }
      previous = i;
    }
  }

  if (previous != n)
  {
    succ[previous] = n;
    od->raft = rbeg[previous] + rexs[previous];
  }

  od->last = previous;
} /* OdCheck */

static void OdAdd(order *od,
                  int *idrw,
                  int x,
                  int newsze)
{
  int n, *adjn, *rbeg, *rexs, *rlen, *pres, *succ;

  n = od->nnod;
  adjn = od->adjn;
  rbeg = od->rbeg;
  rexs = od->rexs;
  rlen = od->rlen;
  pres = od->pres;
  succ = od->succ;

  if (newsze <= rlen[x])
    return;

  if (od->raft + newsze > od->nn0)
    OdCheck(od, idrw);

  if (od->raft + newsze > od->nn0)
  {
    ErrorProc(OutOfSpc, "OdAdd");
    ShutDown();
    exit(0);
  }

  if (pres[x] != n)
    rlen[pres[x]] += rlen[x];

  iCopy(rexs[x], adjn + rbeg[x], adjn + od->raft);
  rbeg[x] = od->raft;
  rlen[x] = newsze;
  od->raft += newsze;

  if (pres[x] == n)
  {
    if (succ[x] == n)
      od->head = x;
    else
      od->head = succ[x];
  }

  else
  {
    if (succ[x] == n)
      succ[pres[x]] = x;
    else
      succ[pres[x]] = succ[x];
  }

  if (succ[x] != n)
    pres[succ[x]] = pres[x];

  if (od->last != x)
  {
    succ[od->last] = x;
    pres[x] = od->last;
  }

  succ[x] = n;
  od->last = x;
} /* OdAdd */

static int OdComb(order *od,
                  int *idrw,
                  int *marker,
                  int *isize,
                  int *ilink,
                  int *osize,
                  int xsize,
                  int *xset)
{
  int i, n, new, rlen, x, icur;

  n = od->nnod;
  rlen = 0;

  if (xsize == 0)
    new = n;
  else if (xsize == 1)
    new = xset[0];
  else
  {
    new = xset[0];
    for (i = 1; i < xsize; ++i)
      rlen += 1 + isize[xset[i]];

    idrw[new] = 1;
    osize[new] = 0;

    for (icur = new; ilink[icur] != n; icur = ilink[icur])
      ;
    isize[new] += rlen;

    for (i = 1; i < xsize; ++i)
    {
      x = xset[i];

      idrw[x] = 0;
      marker[x] = true;

      ilink[icur] = x;

      for (icur = x; ilink[icur] != n; icur = ilink[icur])
        ;

      isize[x] = 0;
    }
  }

  return (new);
} /* OdComb */

static int OdSelect(order *od,
                    xlist *elist,
                    int *idrw,
                    int *marker,
                    int *isize,
                    int *ilink,
                    int *oinfo,
                    int *osize,
                    int x,
                    int *rsze,
                    int *rcrw,
                    int *iptr1,
                    int *iptr2,
                    int *mask2,
                    int *e,
                    int *p)
{
  int absorp, old, i, j, n, esze, y, z, l, f, t, stopt, s,
      o, stops, indsze, xdeg, e0, ssze, *slist, tsze,
      *tlist, sze, *adjn, *rbeg, *rexs, *rlen, *rend;

  adjn = od->adjn;
  rbeg = od->rbeg;
  rexs = od->rexs;
  rlen = od->rlen;
  rend = od->rend;
  n = od->nnod;
  slist = iptr1;

  e0 = *e;
  OdArriv(od, idrw, marker, isize, x, &xdeg, rsze, &esze, rcrw);

  delXt(elist, x);

  OdRenew(od, ilink, x, xdeg, e, p);

  for (i = n - esze; i < n; ++i)
  {
    idrw[rcrw[i]] = 0;
    marker[rcrw[i]] = true;
  }

  marker[x] = true;
  iSet(*rsze, true, marker, rcrw);

  ssze = 0;
  for (i = 0; i < *rsze;)
  {
    y = rcrw[i];

    if (idrw[y] == 0 || idrw[y] == 3)
    {
      ErrorProc(SysError, NULL);
      ShutDown();
      exit(0);
    }

    f = rbeg[y];
    for (t = f, stopt = f + rend[y]; t < stopt; ++t)
    {
      z = adjn[t];
      if (idrw[z] == 3)
      {
        adjn[f++] = z;

        if (!mask2[z])
        {
          slist[ssze++] = z;
          mask2[z] = true;
        }
      }
    }
    rend[y] = f - rbeg[y];

    for (stopt = rbeg[y] + rexs[y]; t < stopt; ++t)
    {
      z = adjn[t];
      if (!marker[z])
        adjn[f++] = z;
    }

    rexs[y] = f - rbeg[y];

    if (rexs[y] == 0)
    {
      OdRenew(od, ilink, y, xdeg - (*e - e0), e, p);
      idrw[y] = 0;
      marker[y] = true;

      (*rsze)--;
      iSwap(i, *rsze, rcrw);
    }

    else
    {
      if (rexs[y] >= rlen[y])
      {
        ErrorProc(SysError, NULL);
        ShutDown();
        exit(0);
      }

      if (rexs[y] > rend[y])
        adjn[rbeg[y] + rexs[y]] = adjn[rbeg[y] + rend[y]];

      rexs[y]++;

      adjn[rbeg[y] + rend[y]] = x;
      rend[y]++;

      i++;
    }
  }

  iSet(ssze, false, mask2, slist);

  if (*rsze == 0)
  {
    idrw[x] = 0;
    marker[x] = true;
  }

  else
  {
    idrw[x] = 3;

    rend[x] = 0;
    rexs[x] = 0;
    if (*rsze > rlen[x])
      OdAdd(od, idrw, x, *rsze);

    rexs[x] = *rsze;
    iCopy(*rsze, rcrw, adjn + rbeg[x]);

    tsze = 0;
    tlist = iptr2;
    for (i = 0; i < ssze; ++i)
    {
      y = slist[i];
      old = marker[y];
      marker[y] = true;

      absorp = true;

      indsze = n;
      l = n;

      f = rbeg[y];
      for (t = f, stopt = f + rexs[y]; t < stopt; ++t)
      {
        z = adjn[t];
        if (idrw[z] != 0)
        {
          adjn[f++] = z;

          if (marker[z])
          {
            l--;
            slist[l] = z;

            if (!mask2[z])
            {
              for (s = rbeg[z], stops = rbeg[z] + rexs[z];
                   s < stops && marker[adjn[s]]; ++s)
                ;

              if (s == stops)
              {
                indsze--;
                iSwap(l, indsze, slist);
              }

              mask2[z] = true;
              tlist[tsze++] = z;
            }
          }
          else
            absorp = false;
        }
      }

      marker[y] = old;
      rexs[y] = f - rbeg[y];

      if (indsze < n)
      {
        z = OdComb(od, idrw, marker,
                   isize, ilink, osize,
                   n - indsze, slist + indsze);

        idrw[z] = 1;

        sze = 0;
        for (j = l; j < indsze; ++j)
        {
          o = slist[j];
          sze += 1 + isize[o];
          idrw[o] = 2;
          oinfo[o] = z;
        }
        osize[z] = max(osize[z], sze);
      }

      if (absorp)
      {
        idrw[y] = 0;
        marker[y] = true;
      }
    }

    iSet(tsze, false, mask2, tlist);
  }

  marker[x] = (idrw[x] == 0);

  for (t = 0; t < *rsze; ++t)
  {
    z = rcrw[t];
    marker[z] = (idrw[z] == 0);
  }

  return (false);
} /* OdSelect */

static int OdOrder(order *od,
                   int *idrw,
                   int *marker,
                   int *isize,
                   int x,
                   int *iptr1)
{
  int rsze, esze, deg;

  OdArriv(od, idrw, marker, isize,
          x, &deg, &rsze, &esze, iptr1);

  return deg;
} /* OdOrder */

static void OdModf(order *od,
                   xlist *elist,
                   int *idrw,
                   int *marker,
                   int *isize,
                   int *oinfo,
                   int rsze,
                   int *rcrw,
                   int *iptr1)
{

  int i, x, deg;

  for (i = 0; i < rsze; ++i)
  {
    x = rcrw[i];
    if (idrw[x] == 2)
      if (idrw[oinfo[x]] == 0 || idrw[oinfo[x]] == 3)
        idrw[x] = 1;

    if (idrw[x] == 1)
    {
      deg = OdOrder(od, idrw, marker, isize, x, iptr1);
      putXt(elist, x, deg - isize[x]);
    }

    else
      delXt(elist, x);
  }
} /* OdModf */

static void OdProc(order *od,
                   xlist *xt,
                   int *bptra,
                   int *bptrb,
                   int *iptr1,
                   int *iptr2,
                   int *iptr3,
                   int *iptr4,
                   int *iptr5,
                   int *iptr6,
                   int *iptr7,
                   int *iptr8,
                   int *iptr9,
                   int *mapi,
                   int *p)
{
  int *mask2, nfilsub, beord, *marker, *idrw, i, n, e, x,
      y, rsze, deg, mindeg, xsize, sze, j, *isize, *ilink,
      *oinfo, *osize, *rcrw, *ldpt, *slist;
  xlist *elist;

  beord = true;
  elist = xt;
  marker = bptra;
  isize = iptr1;
  ilink = iptr2;
  oinfo = iptr3;
  osize = iptr4;
  rcrw = iptr5;
  ldpt = iptr6;
  slist = iptr7;
  idrw = mapi;
  mask2 = bptrb;

  OdSet(od, true, elist, idrw, marker,
        isize, ilink, oinfo, osize, &e, p);

  n = od->nnod;

  iSet(n, 0, ldpt, NULL);

  nfilsub = false;
  for (; e < n && !nfilsub;)
  {

    infXt(elist);

    if (!getXt(elist, &y, &mindeg))
    {
      printf("\n No new nodes e=%d  n=%d", e, n);
      printf(" Node status: ");

      for (i = 0; i < n; ++i)
        if (idrw[i] == 1)
          printf("A\n");
        else if (idrw[i] == 2)
          printf("\n O%d: rlen=%d oinfo=%d\n",
                 i, isize[i], oinfo[i]);

      ErrorProc(SysError, NULL);
      ShutDown();
      exit(0);
    }

    xsize = 0;
    for (; beord;)
    {
      if (!getXt(elist, &x, &deg) || deg > mindeg)
        break;

      if (idrw[x] != 1)
        delXt(elist, x);

      else
      {
        nextXt(elist);
        if (!ldpt[x])
        {

          nfilsub = OdSelect(od, elist, idrw, marker,
                             isize, ilink, oinfo, osize, x,
                             &rsze, rcrw, iptr8, iptr9, mask2,
                             &e, p);

          if (!nfilsub)
          {
            ldpt[x] = 2;
            slist[xsize++] = x;

            for (i = 0; i < rsze; ++i)
            {
              y = rcrw[i];
              if (!ldpt[y])
              {
                ldpt[y] = 1;
                slist[xsize++] = y;
              }
            }
          }

          if (!beord)
            break;
        }
      }
    }

    if (!nfilsub)
    {
      sze = 0;
      for (j = 0; j < xsize; ++j)
      {
        y = slist[j];
        if (ldpt[y] == 1 && idrw[y] != 0)
          slist[sze++] = y;
        ldpt[y] = 0;
      }

      OdModf(od, elist, idrw, marker,
             isize, oinfo, sze, slist, iptr8);
    }
  }

  if (e < n)
  {
    sze = 0;
    for (i = 0; i < n; ++i)
      if (idrw[i] == 2 || idrw[i] == 1)
        iptr8[sze++] = i;

    x = OdComb(od, idrw, marker, isize, ilink, osize,
               sze, iptr8);

    OdRenew(od, ilink, x, n - e - 1, &e, p);
  }

} /* OdProc */

int GetOrder(order *od,
             int *p)
{
  int bbufs = 2, *bbuf[2] = {0},
      ibufs = 9, *ibuf[9] = {0},
      n, *iwmd;
  xlist *xt;

  n = od->nnod;

  xt = XtAlloc(n, n + 1);
  iwmd = iAlloc(n, "iptr21, GetOrder");

  iPtAlloc(ibufs, n, ibuf, "ibuf, GetOrder");
  iPtAlloc(bbufs, n, bbuf, "bbuf, GetOrder");

  OdProc(od, xt, ibuf[0], ibuf[1], ibuf[2], ibuf[3], ibuf[4], ibuf[5],
         ibuf[6], ibuf[7], ibuf[8], iwmd, bbuf[0], bbuf[1], p);

  XtFree(&xt);
  iFree(&iwmd);
  iPtFree(ibufs, ibuf);
  iPtFree(bbufs, bbuf);

  return true;
} /* GetOrder */
