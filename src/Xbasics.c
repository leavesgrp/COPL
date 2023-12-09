#include "LPcross.h"
//#define TIME_COUNT

pivbas *PbsAlloc(int n, char *info)
{
  pivbas *r = NULL;

  if (n)
  {
    r = (pivbas *)calloc(n, sizeof(pivbas));
    if (!r)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }
  }
  return r;
} /* PbsAlloc */

void PbsFree(pivbas **x)
{
  pivbas *r = *x;

  if (r)
  {
    free(r);
    (*x) = NULL;
  }
} /* PbsFree */

optbas *ObsAlloc(optmod *opt,
                 char *info)
{
  optbas *r;
  int m = opt->m,
      n = opt->n;

  r = (optbas *)calloc(1, sizeof(optbas));
  if (!r)
  {
    ErrorProc(NOT_MEMSPC, info);
    ShutDown();
    exit(0);
  }

  r->maxnrow = m;
  r->maxncol = n;
  r->nrow = m;
  r->ncol = n;

  r->basis = PbsAlloc(m, info);
  r->rsk = iAlloc(m, info);
  r->csk = iAlloc(n, info);
  r->rx = dAlloc(m, info);
  r->cx = dAlloc(n, info);
  r->xb = dAlloc(m, info);
  r->y = dAlloc(m, info);
  r->cs = dAlloc(n, info);

  return r;
} /* ObsAlloc */

void ObsFree(optbas **x)
{
  optbas *r = *x;

  if (r)
  {
    PbsFree(&r->basis);
    iFree(&r->rsk);
    iFree(&r->csk);
    dFree(&r->rx);
    dFree(&r->cx);
    dFree(&r->xb);
    dFree(&r->y);
    dFree(&r->cs);
    free(r);
    *x = NULL;
  }
} /* ObsFree */

gsdec *LufAlloc(int m,
                int n,
                int naat,
                int bfac2,
                char *info)
{
  gsdec *r;

  if (m <= 0 || n < m || naat < 0)
  {
    ErrorProc(SZE_ERROR, info);
    ShutDown();
    exit(0);
  }

  r = (gsdec *)calloc(1, sizeof(gsdec));

  if (!r)
  {
    ErrorProc(NOT_MEMSPC, info);
    ShutDown();
    exit(0);
  }

  r->nrow0 = m;
  r->ncol0 = n;
  r->nrow = m;
  r->ncol = n;

  r->wsze = max(2 * naat, 100) + 2 * m;

  r->toldrop = 1.0e-12;
  r->tolfill = 1.0;
  r->tolrpiv = 0.01;
  r->tolapiv = 1.0e-8;
  r->tolupiv = 1.0e-6;
  r->tolstb = 1.0e15;

  r->trhsch = 10;
  r->trhuov = 2;
  r->trhrfq = 100;

  r->dropl = false;
  r->dropu = false;

  r->val = dAlloc(r->wsze, info);
  r->subi = iAlloc(r->wsze, info);
  r->subj = iAlloc(r->wsze, info);

  r->chead = iAlloc(n, info);
  r->ctail = iAlloc(n, info);
  r->ccur = iAlloc(n, info);
  r->cprev = iAlloc(r->wsze, info);
  r->cnext = iAlloc(r->wsze, info);

  r->ufir = iAlloc(m, info);
  r->usze = iAlloc(m, info);
  r->uuse = iAlloc(m, info);
  r->unex = iAlloc(m, info);
  r->upre = iAlloc(m, info);

  if (bfac2)
  {
    r->nzrow = (xlist *)calloc(1, sizeof(xlist));
    r->nzcol = (xlist *)calloc(1, sizeof(xlist));
    if (!r->nzrow || !r->nzcol)
    {
      ErrorProc(NOT_MEMSPC, info);
      ShutDown();
      exit(0);
    }
  }
  else
  {
    r->nzrow = XtAlloc(m, n);
    r->nzcol = XtAlloc(n, m);
  }

  r->p = iAlloc(m, info);
  r->invp = iAlloc(m, info);
  r->invq = iAlloc(n, info);
  r->rowsta = iAlloc(m, info);

  return r;
} /* LufAlloc */

int LufRenew(gsdec *r,
             int wsze,
             char *info)
{
  if (r->wsze)
  {
    dFree(&r->val);
    iFree(&r->subi);
    iFree(&r->subj);
    iFree(&r->cprev);
    iFree(&r->cnext);
  }

  r->wsze = wsze;

  if (wsze)
  {
    r->val = dAlloc(wsze, info);
    r->subi = iAlloc(wsze, info);
    r->subj = iAlloc(wsze, info);
    r->cprev = iAlloc(wsze + r->nrow, info);
    r->cnext = iAlloc(wsze + r->ncol, info);
  }
  return lu_ok;
} /* LufRenew */

void LufFree(gsdec **lfac)
{
  gsdec *r = *lfac;

  if (r)
  {
    dFree(&r->val);
    iFree(&r->subi);
    iFree(&r->subj);
    iFree(&r->chead);
    iFree(&r->ctail);
    iFree(&r->ccur);
    iFree(&r->cprev);
    iFree(&r->cnext);
    iFree(&r->ufir);
    iFree(&r->usze);
    iFree(&r->uuse);
    iFree(&r->unex);
    iFree(&r->upre);
    XtFree(&r->nzrow);
    XtFree(&r->nzcol);
    iFree(&r->p);
    iFree(&r->invp);
    iFree(&r->invq);
    iFree(&r->rowsta);
    free(r);
  }
  *lfac = NULL;
} /* lu_delete */

int blkord_maxtran(matrix *at,
                   int rsze,
                   int rset[],
                   int invrset[],
                   int csze,
                   int cset[],
                   int cmark[],
                   int wisze,
                   int wimem[])
{
  int asgn;
  int k, i, j, sze, s, t, nasgn, numasgn, nil, itemp,
      *stack, *cur, *sfir, *rvis;
  array *aj, *ak;

  /*
   * Find maximum tranversal.
   */
  if (wisze < 4 * rsze)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  stack = wimem;      /* Column stack (implicit via rset).     */
  cur = stack + rsze; /* Pointer to current in column.         */
  sfir = cur + rsze;  /* Pointer to first sub-diagonal element
                         in column.                            */
  rvis = sfir + rsze; /* Mark if a row has been visited.       */

  for (k = 0; k < rsze; ++k)
    rvis[k] = 0;

  /*
   * Simple greedy step.
   */

  for (k = 0, itemp = 0; k < csze && itemp < rsze; ++k)
  {
    /*
     * Try a new column.
     */

    ak = at->ia + cset[k];
    asgn = false;

    for (s = 0; s < ak->nn0 && !asgn; ++s)
    {
      i = invrset[ak->ja[s]];
      asgn = (itemp <= i && i < rsze);
    }

    if (asgn)
    {
      /*
       * Cheap assigment.
       */

      iSwap(itemp, i, rset);
      invrset[rset[itemp]] = itemp;
      invrset[rset[i]] = i;

      iSwap(itemp, k, cset);
      if (cmark)
        iSwap(itemp, k, cmark);

      ++itemp;
    }
  }

  nil = 1;
  nasgn = 0; /* Columns that can not be assigned. */
  for (k = 0; k < csze - nasgn && k < rsze;)
  {
    /*
     * Try a new column.
     */

    ak = at->ia + cset[k];
    asgn = false;

    for (s = 0; s < ak->nn0 && !asgn; ++s)
    {
      i = invrset[ak->ja[s]];
      asgn = (k <= i && i < rsze);
    }

    sfir[k] = s;

    if (asgn)
    {
      iSwap(k, i, rset);
      invrset[rset[k]] = k;
      invrset[rset[i]] = i;

      ++k;
    }
    else
    {
      asgn = false;

      for (cur[k] = 0; cur[k] < ak->nn0 && !asgn; ++cur[k])
      {
        i = invrset[ak->ja[cur[k]]];
        if (i < rsze)
        {
          sze = 1;
          stack[0] = i;
          cur[i] = 0;
          rvis[i] = nil;

          do
          {
            t = stack[sze - 1];

            if (t >= k)
            {
              printf("\n\n system error: t>=k.\n");
              ShutDown();
              exit(0);
            }

            j = cset[t];
            aj = at->ia + j;

            for (s = sfir[t]; s < aj->nn0 && !asgn; ++s)
            {
              i = invrset[aj->ja[s]];
              if (k <= i && i < rsze)
              {
                asgn = true;
                stack[sze++] = i;
              }
            }
            sfir[t] = s;

            if (!asgn)
            {
              for (; cur[t] < aj->nn0; ++cur[t])
              {
                i = invrset[aj->ja[cur[t]]];
                if (i < rsze && i != t && rvis[i] != nil)
                {
                  rvis[i] = nil;
                  cur[i] = 0;
                  stack[sze++] = i;

                  break;
                }
              }

              if (cur[t] == aj->nn0)
              {
                --sze;
              }
            }
          } while (sze && !asgn);
        }
      }

      if (asgn)
      {
        itemp = rset[stack[0]];
        for (t = 1; t < sze; ++t)
          rset[stack[t - 1]] = rset[stack[t]];
        rset[stack[sze - 1]] = itemp;

        for (t = 0; t < sze; ++t)
          invrset[rset[stack[t]]] = stack[t];

        i = stack[sze - 1];

        iSwap(k, i, rset);
        invrset[rset[k]] = k;
        invrset[rset[i]] = i;

        ++k;
      }
      else
      {
        ++nasgn;
        iSwap(csze - nasgn, k, cset);

        if (cmark)
          iSwap(csze - nasgn, k, cmark);

        sfir[k] = 0;
      }
    }
    ++nil;
  }
  numasgn = k;
  return (numasgn);
} /* blkord_maxtran */

int blkord_maxlwrsqblk(matrix *at,
                       int rsze,
                       int rset[],
                       int invrset[],
                       int csze,
                       int cset[],
                       int cmark[],
                       int *maxtrsze,
                       int wisze,
                       int wimem[])
{
  int t, k, i, j, maxtran, blksze, sze,
      *list, *mask;
  array *aj;

  maxtran = blkord_maxtran(at, rsze, rset, invrset, csze, cset, cmark, wisze, wimem);
  *maxtrsze = maxtran;

  list = wimem;
  mask = list + rsze;

  iZero(rsze, mask, NULL);

  sze = 0;
  for (k = maxtran; k < csze; ++k)
  {
    j = cset[k];
    aj = at->ia + j;
    for (t = 0; t < aj->nn0; ++t)
    {
      i = invrset[aj->ja[t]];
      if (i < rsze)
      {
        if (maxtran <= i)
        {
          printf("\nError 1");
          exit(0);
        }

        if (!mask[i])
        {
          mask[i] = 1;
          list[sze++] = i;
        }
      }
    }
  }

  while (sze)
  {
    --sze;

    i = list[sze];
    j = cset[i];
    aj = at->ia + j;
    for (t = 0; t < aj->nn0; ++t)
    {
      i = invrset[aj->ja[t]];
      if (i < rsze && !mask[i])
      {
        mask[i] = 1;
        list[sze++] = i;
      }
    }
  }

  blksze = 0;
  for (k = 0; k < maxtran; ++k)
  {
    if (!mask[k])
    {
      iSwap(blksze, k, cset);
      if (cmark)
        iSwap(blksze, k, cmark);

      iSwap(blksze, k, rset);

      invrset[rset[k]] = k;
      invrset[rset[blksze]] = blksze;

      ++blksze;
    }
  }

  /*
   * Checking code which can be removed.
   */
  for (k = blksze; k < *maxtrsze; ++k)
  {
    j = cset[k];
    aj = at->ia + j;
    for (t = 0; t < aj->nn0; ++t)
      if (aj->ja[t] == rset[k])
        break;

    if (t == aj->nn0)
      exit(0);
  }

  return (blksze);
} /* blkord_maxlwrsqblk */

int blkord_maxuprsqblk2(int ncol,
                        matrix *a,
                        matrix *at,
                        int rsze,
                        int rset[],
                        int invrset[],
                        int csze,
                        int cset[],
                        int cmark[],
                        int wisze,
                        int wimem[],
                        int wzisze,
                        int wzimem[])
{
  int k, s, t, i, j, maxtran, blksze, blksze0, itemp, *invj;
  array *ai, *aj;

  if (wisze < 4 * rsze || wzisze < ncol)
  {
    printf("\n\n system error, work space size.\n");
    ShutDown();
    exit(0);
  }

  maxtran = blkord_maxtran(at, rsze, rset, invrset,
                           csze, cset, cmark, wisze, wimem);
  blksze = maxtran;
  invj = wzimem;

  for (k = 0; k < maxtran; ++k)
    invj[cset[k]] = ncol - k;

  for (k = 0; k < blksze;)
  {
    j = cset[k];
    aj = at->ia + j;
    for (t = 0; t < aj->nn0; ++t)
    {
      i = invrset[aj->ja[t]];
      if (maxtran <= i && i < rsze)
        break;
    }

    if (t < aj->nn0)
    {
      /*
       * Perform symmetric permutation.
       */

      --blksze;

      iSwap(k, blksze, cset);
      if (cmark)
        iSwap(k, blksze, cmark);

      invj[cset[k]] = ncol - k;
      invj[cset[blksze]] = ncol - blksze;

      iSwap(k, blksze, rset);

      invrset[rset[k]] = k;
      invrset[rset[blksze]] = blksze;
    }
    else
      ++k;
  }

  blksze0 = maxtran;
  while (blksze < blksze0)
  {
    itemp = blksze;
    for (t = blksze; t < blksze0; ++t)
    {
      i = rset[t];
      ai = a->ia + i;
      for (s = 0; s < ai->nn0; ++s)
      {
        k = ncol - invj[ai->ja[s]];

        if (k < maxtran && cset[k] != ai->ja[s])
        {
          printf("\n\n system error, invj wrong.\n");
          ShutDown();
          exit(0);
        }

        if (k < blksze)
        {
          /*
           * Perform symmetric permutation.
           */

          --blksze;

          iSwap(k, blksze, cset);
          if (cmark)
            iSwap(k, blksze, cmark);

          invj[cset[k]] = ncol - k;
          invj[cset[blksze]] = ncol - blksze;

          iSwap(k, blksze, rset);

          invrset[rset[k]] = k;
          invrset[rset[blksze]] = blksze;
        }
      }
    }
    blksze0 = itemp;
  }

  for (k = 0; k < maxtran; ++k)
    invj[cset[k]] = 0;

  return (blksze);
} /* blkord_maxuprsqblk2 */

int basord_blklwr(int nrow,
                  int ncol,
                  int nrowc,
                  int inew[],
                  int iinv[],
                  int rsze,
                  int rset[],
                  matrix *a,
                  matrix *at,
                  pivbas *basis,
                  int wisze,
                  int wimem[],
                  int wzisze,
                  int wzimem[])
{
  int i, j, k, t, cszel, blksze, wiszel, nrowc0,
      *csetl, *jsta, *unfin, *cmark, *wimeml;

  if (wisze < 6 * nrowc || wzisze < ncol)
  {
    printf("\n\n system error, work space size.\n");
    ShutDown();
    exit(0);
  }
  nrowc0 = nrowc;

  csetl = wimem;
  unfin = csetl + nrowc;

  wiszel = wisze - 2 * nrowc;
  wimeml = wimem + 2 * nrowc;

  iZero(nrowc, unfin, NULL);

  for (k = 0; k < rsze; ++k)
  {
    if (inew[rset[k]] >= nrowc)
    {
      printf("\n rset[k]=%d", rset[k]);
      printf("  k=%d"
             "  nrowc=%d",
             k, nrowc);
      printf(" inew[rset[k]]=%d", inew[rset[k]]);
      printf("\n\n system error, rset error inew[rset[k]]>=nrowc.\n");
      ShutDown();
      exit(0);
    }

    unfin[inew[rset[k]]] = 1;
  }

  /*
   * Remove finished rows and artificial variables.
   */
  for (k = 0; k < nrowc;)
  {
    i = iinv[k];
    j = basis[i].j;
    if (basis[i].vt == PIV_ART || basis[i].vt == PIV_ROW)
    {
      if (i == j)
        ++k;
      else
      {
        pivlp_swapbasis(i, j, basis);
        iSwap(inew[i], inew[j], unfin);
      }
    }
    else
      ++k;
  }

  for (k = 0; k < nrowc;)
  {
    i = iinv[k];
    if (!unfin[k] && (basis[i].vt == PIV_ART || basis[i].vt == PIV_ROW))
    {
      /*
       * The row/artificial variable is finished.
       */

      j = basis[i].j;
      if (i != j)
        exit(0);

      --nrowc;

      iSwap(k, nrowc, iinv);

      inew[iinv[k]] = k;
      inew[iinv[nrowc]] = nrowc;

      iSwap(k, nrowc, unfin);

      if (unfin[nrowc])
        exit(0);
    }
    else
      ++k;
  }

  cmark = unfin;
  cszel = 0;
  for (k = 0; k < nrowc; ++k)
  {
    if (!unfin[k])
    {
      i = iinv[k];

      if (basis[i].vt == PIV_COL)
      {
        cmark[cszel] = i;
        csetl[cszel++] = basis[i].j;
      }
    }
  }

  /*
   * csetl[0..cszel-1] is the basis finshed columns.
   * All finished rows and artificial variable
   * has been removed.
   */

  blksze = blkord_maxuprsqblk2(ncol,
                               a, at,
                               nrowc, iinv, inew,
                               cszel, csetl, cmark,
                               wiszel, wimeml,
                               wzisze, wzimem);

  /*
   * Move death basic variables to death rows.
   */

  jsta = cmark + nrowc;

  /* jsta[k] : 0 - basis[k].j is not finished.
   *           1 - basis[k].j is finished but can not be removed.
   *           2 - basis[k].j is finished and can be removed.
   */

  iZero(nrowc, jsta, NULL);

  for (k = 0; k < blksze; ++k)
  {
    if (basis[cmark[k]].vt != PIV_COL)
    {
      printf("\n cmark error");
      ShutDown();
      exit(0);
    }

    if (basis[cmark[k]].j != csetl[k])
    {
      printf("\n1. Error cmark");

      printf("\n k=%d", k);
      printf("  cseti=%d", csetl[k]);
      printf(" != basis[%d"
             "].j=%d\n",
             cmark[k], basis[cmark[k]].j);
      ShutDown();
      exit(0);
    }

    jsta[inew[cmark[k]]] = 2;
  }

  for (k = blksze; k < cszel; ++k)
  {
    if (basis[cmark[k]].vt != PIV_COL)
    {
      printf("\n cmark error");
      ShutDown();
      exit(0);
    }

    if (basis[cmark[k]].j != csetl[k])
    {
      printf("\n1. Error cmark");

      printf("\nk=%d" IFMT, k);
      printf("  cseti=%d", csetl[k]);
      printf(" != basis[%d"
             "].j=%d\n",
             cmark[k], basis[cmark[k]].j);
      ShutDown();
      exit(0);
    }

    jsta[inew[cmark[k]]] = 1;
  }

  t = nrowc;
  for (k = 0; k < t;)
  {
    if (jsta[k] == 2)
      ++k;
    else
    {
      --t;
      pivlp_swapbasis(iinv[t], iinv[k], basis);
      iSwap(t, k, jsta);
    }
  }

  if (t != blksze)
  {
    printf("\n internal error, t!=blksze.\n");
    ShutDown();
    exit(0);
  }

  t = 0;
  for (k = blksze; k < nrowc; ++k)
  {
    if (jsta[k] == 0)
      rset[t++] = iinv[k];
  }

  if (t != rsze)
  {
    printf("\n LAst t!=rsze");
    printf("\n nrowc=%d", nrowc);
    printf("\nt=%d"
           " rsze=%d\n",
           t, rsze);
    ShutDown();
    exit(0);
  }

  for (k = 0; k < blksze; ++k)
  {
    if (jsta[k] != 2)
    {
      printf("\n internal error, jsta[k]!=2.\n");
      ShutDown();
      exit(0);
    }

    --nrowc;

    iSwap(k, nrowc, iinv);

    inew[iinv[k]] = k;
    inew[iinv[nrowc]] = nrowc;
  }
  return (nrowc0 - nrowc);
} /* bas_blklwr */

int lu_putncol(gsdec *lf,
               int ncol)
{
  if (ncol > lf->ncol0)
    return (lu_arg2);

  lf->ncol = ncol;

  return (lu_ok);
} /* lu_putncol */

int lu_putnrow(gsdec *lf,
               int nrow)
{
  if (nrow > lf->nrow0)
    return (lu_arg2);

  lf->nrow = nrow;

  return (lu_ok);
} /* lu_putncol */

void pivlp_formxb(int nrow,
                  int ncol,
                  matrix *at,
                  int rbk[],
                  double rbl[],
                  double rbu[],
                  pivbas *basis,
                  double rx[],
                  double cx[],
                  double xb[],
                  gsdec *lu,
                  int wrsze,
                  double wrmem[])
{
  int k;
  double *res, *dx;

  if (wrsze < 2 * nrow)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  res = wrmem;
  dx = wrmem + nrow;

  /*
   * Compute res
   */

  for (k = 0; k < nrow; ++k)
  {
    if (basis[k].vt == PIV_ROW)
      rx[basis[k].j] = xb[k];
    else if (basis[k].vt == PIV_COL)
      cx[basis[k].j] = xb[k];
  }

  if (rx)
  {
    for (k = 0; k < nrow; ++k)
      res[k] = rx[k];
  }
  else
  {
    for (k = 0; k < nrow; ++k)
      res[k] = optdef_getrbli(rbk, rbl, k);
  }

  mTimesv(false, nrow, ncol, -1.0, at, cx, 1.0, res);

  for (k = 0; k < nrow; ++k)
    if (basis[k].vt == PIV_ART)
      res[basis[k].j] -= xb[k];

  lu_ftran(lu, res, dx);

  for (k = 0; k < nrow; ++k)
  {
    xb[k] += dx[k];

    if (basis[k].vt == PIV_ROW)
      rx[basis[k].j] = xb[k];
    else if (basis[k].vt == PIV_COL)
      cx[basis[k].j] = xb[k];
  }

  if (rx)
  {
    for (k = 0; k < nrow; ++k)
      res[k] = rx[k];
  }
  else
  {
    for (k = 0; k < nrow; ++k)
      res[k] = optdef_getrbli(rbk, rbl, k);
  }

  mTimesv(false, nrow, ncol, -1.0, at, cx, 1.0, res);

  for (k = 0; k < nrow; ++k)
    if (basis[k].vt == PIV_ART)
      res[basis[k].j] -= xb[k];
} /* pivlp_formxb */

void pivlp_xb2rxcx(int nrow,
                   pivbas *basis,
                   double rx[],
                   double cx[],
                   double xb[])
{
  int k, j;

  for (k = 0; k < nrow; ++k)
  {
    j = basis[k].j;

    if (basis[k].vt == PIV_ROW)
      rx[j] = xb[k];
    else if (basis[k].vt == PIV_COL)
      cx[j] = xb[k];
  }
} /* pivlp_xb2rxcx */

int pivlp_insart(gsdec *lu,
                 int nrow,
                 int nrowc,
                 int ncol,
                 int iinv[],
                 int inew[],
                 int rsety,
                 pivbas *basis,
                 int rsk[],
                 int csk[],
                 double rx[],
                 double cx[],
                 double xb[],
                 double yb[],
                 int wisze,
                 int wimem[])
{
  int i, j, k, num = 0, *insj;

  if (wisze < 2 * nrowc)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  insj = wimem;

  lu_insart(lu, insj, wisze - nrowc, wimem + nrowc);

  if (iinv)
  {
    for (k = 0; k < nrowc; ++k)
    {
      if (insj[k] != nrowc)
      {
        /*
         * The ith basic variables should be removed.
         */

        ++num;

        i = iinv[k];

        switch (basis[i].vt)
        {
        case PIV_ART:
          printf("\n\n system error, not art var.\n");
          ShutDown();
          exit(0);
        case PIV_ROW:
          j = basis[i].j;
          if (rsk)
            rsk[j] = SPBAS;

          if (xb && cx)
            rx[j] = xb[i];
          break;
        case PIV_COL:
          j = basis[i].j;
          if (csk)
            csk[j] = SPBAS;

          if (xb && cx)
            cx[j] = xb[i];
          break;
        }

        if (insj[k] >= nrowc)
        {
          printf("\n\n internal error, insj wrong.\n");
          ShutDown();
          exit(0);
        }

        basis[i].vt = PIV_ART;
        basis[i].j = iinv[insj[k]];
        if (xb)
          xb[i] = 0.0;

        if (rsety && yb)
          yb[i] = 0.0;
      }
    }
  }
  else
  {
    for (i = 0; i < nrow; ++i)
    {
      if (insj[i] != nrow)
      {
        /*
         * The ith basic variables should be removed.
         */

        ++num;

        switch (basis[i].vt)
        {
        case PIV_ART:
          printf("\n\n system error, not art var.\n");
          ShutDown();
          exit(0);
        case PIV_ROW:
          j = basis[i].j;
          if (rsk)
            rsk[j] = SPBAS;

          if (xb && cx)
            rx[j] = xb[i];
          break;
        case PIV_COL:
          j = basis[i].j;
          if (csk)
            csk[j] = SPBAS;

          if (xb && cx)
            cx[j] = xb[i];
          break;
        }

        basis[i].vt = PIV_ART;
        basis[i].j = insj[i];
        if (xb)
          xb[i] = 0.0;

        if (rsety && yb)
          yb[i] = 0.0;
      }
    }
  }
  return (num);
} /* pivlp_insart */

int pivlp_factor1(int nrow,
                  int nrowc,
                  int iinv[],
                  int inew[],
                  matrix *at,
                  int inca,
                  pivbas *basis,
                  gsdec *lu,
                  int wisze,
                  int wimem[],
                  int wrsze,
                  double wrmem[])
{
  int i, k, s, t, try_, sze;
  gsmsg *lfi = &lu->info;
  int lresp = lu_ok;
  array *aj;
#ifdef TIME_COUNT
  clock_t Tm;
  Tm = GetTime();
  printf("pivlp_factor1: inca = %d\n", inca);
#endif
  lresp = lu_putncol(lu, nrowc);
  if (lresp != lu_ok)
    return (lresp);

  lresp = lu_putnrow(lu, nrowc);
  if (lresp != lu_ok)
    return (lresp);

  if (wisze < 2 + 11 * nrowc || wrsze < 2 * nrowc)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  if (!wimem || !wrmem)
  {
    printf("\n\n null work space.\n");
    ShutDown();
    exit(0);
  }

  if (lu->incwsze)
    lresp = LufRenew(lu, max(lu->wsze / 5, 5 * nrow) + lu->wsze, NULL);
  else if (lfi->ngc > 2)
    lresp = LufRenew(lu, max(lu->wsze / 10, 2 * nrow) + lu->wsze, NULL);
#ifdef TIME_COUNT
  printf("pivlp_factor1-A: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif

  for (try_ = 0; try_ < 5 && lresp == lu_ok; ++try_)
  {
    lresp = lu_init(lu);
#ifdef TIME_COUNT
    printf("pivlp_factor1-B%d: lresp = %d\n", try_, lresp);
    Tm = GetTime();
#endif

    if (lresp == lu_ok)
    {
      if (iinv)
      {
        for (k = 0; k < nrowc && lresp == lu_ok; ++k)
        {
          i = iinv[k];
          if (basis[i].vt == PIV_ART)
          {
            if (inca && inew[basis[i].j] < nrowc)
              lresp = lu_inaij(lu, 1.0, inew[basis[i].j], k);
          }
          else if (basis[i].vt == PIV_ROW)
          {
            if (inew[basis[i].j] < nrowc)
              lresp = lu_inaij(lu, -1.0, inew[basis[i].j], k);
          }
          else
          {
            aj = at->ia + basis[i].j;
            lresp = lu_inaj2(lu, k, aj->nn0, aj->an, aj->ja, nrowc, inew);
          }
        }
      }
      else
      {
        for (i = 0; i < nrow && lresp == lu_ok; ++i)
        {
          if (basis[i].vt == PIV_ART)
          {
            if (inca)
              lresp = lu_inaij(lu, 1.0, basis[i].j, i);
          }
          else if (basis[i].vt == PIV_ROW)
            lresp = lu_inaij(lu, -1.0, basis[i].j, i);
          else
          {
            aj = at->ia + basis[i].j;
            for (t = 0; t < aj->nn0; ++t)
              lresp = lu_inaij(lu, aj->an[t], aj->ja[t], i);
          }
        }
      }

      if (lresp == lu_ok)
      {
#ifdef TIME_COUNT
        printf("pivlp_factor1-B%d: INTO-A - START %d ms\n", try_, GetTime() - Tm);
#endif
        lresp = lu_factor(lu, wisze, wimem, wrsze, wrmem); // HOTSPOT
#ifdef TIME_COUNT
        printf("pivlp_factor1-B%d: INTO-A - AFTER lu_factor %d ms\n", try_, GetTime() - Tm);
#endif
        if (lresp == lu_ok || lresp == lu_insta)
          break;
        else if (lresp == lu_wfull)
          lresp = LufRenew(lu, 5 * nrowc + 2 * lu->wsze, NULL);
        else
        {
          printf("\n\n work space error.\n");
          ShutDown();
          exit(0);
        }
#ifdef TIME_COUNT
        printf("pivlp_factor1-B%d: INTO-A - AFTER LufRenew %d ms\n", try_, GetTime() - Tm);
#endif
      }
      else if (try_ < 1)
      {
#ifdef TIME_COUNT
        printf("pivlp_factor1-B%d: INTO-B\n", try_);
#endif
        /*
         * Increase work area.
         *
         * Count number of non-zeros in the basis.
         */

        sze = 0;

        if (iinv)
        {
          for (k = 0; k < nrowc; ++k)
          {
            i = iinv[k];
            if (basis[i].vt == PIV_ART)
            {
              if (inca && inew[basis[i].j] < nrowc)
                ++sze;
            }
            else if (basis[i].vt == PIV_ROW)
            {
              if (inew[basis[i].j] < nrowc)
                ++sze;
            }
            else
            {
              aj = at->ia + basis[i].j;
              for (t = 0; t < aj->nn0; ++t)
              {
                s = inew[aj->ja[t]];
                if (s < nrowc)
                  ++sze;
              }
            }
          }
        }
        else
        {
          for (i = 0; i < nrow; ++i)
          {
            if (basis[i].vt == PIV_ART)
            {
              if (inca)
                ++sze;
            }
            else if (basis[i].vt == PIV_ROW)
              ++sze;
            else
            {
              aj = at->ia + basis[i].j;
              for (t = 0; t < aj->nn0; ++t)
                ++sze;
            }
          }
        }

        lresp = lu_ok;
        sze = max(3 * nrowc, 5 * sze);
        if (sze > lu->wsze)
        {
          lresp = LufRenew(lu, sze, NULL);
          if (lresp != lu_ok)
            break;
        }
      }
    }
#ifdef TIME_COUNT
    printf("pivlp_factor1-B%d: %d ms\n", try_, GetTime() - Tm);
    Tm = GetTime();
#endif
  }
#ifdef TIME_COUNT
  printf("pivlp_factor1-B%d: %d ms\n", try_, GetTime() - Tm);
  Tm = GetTime();
#endif

  return (lresp);
} /* pivlp_factor1 */

int pivlp_factor2(int prlev,
                  int nrow,
                  int nrowc,
                  int ncol,
                  int iinv[],
                  int inew[],
                  matrix *at,
                  pivbas *basis,
                  int rsk[],
                  int csk[],
                  double rx[],
                  double cx[],
                  double xb[],
                  double yb[],
                  int rsety,
                  gsdec *lu,
                  int wisze,
                  int wimem[],
                  int wrsze,
                  double wrmem[])
{
  int num, rank;
  int lresp;
  int resp = ProcOk;

  lresp = pivlp_factor1(nrow, nrowc,
                        iinv, inew, at, true, basis, lu,
                        wisze, wimem,
                        wrsze, wrmem);

  switch (lresp)
  {
  case lu_insta:
    printf("\n\n unstable lu.\n");
    ShutDown();
    exit(0);
  case lu_ok:
    lu_getrank(lu, &rank);

    if (rank < nrowc)
    {
      num = pivlp_insart(lu, nrow, nrowc, ncol, iinv, inew, rsety, basis,
                         rsk, csk,
                         rx, cx, xb, yb,
                         wisze, wimem);

      if (num)
        printf("\n Insert %d"
               " artificials into the basis.",
               num);
    }
    break;
  case lu_space:
    resp = OutOfSpc;
    break;
  default:
    printf("\n\n system error, unknown case.\n");
    ShutDown();
    exit(0);
  }
  return (resp);
} /* pivlp_factor2 */

int pivlp_prolresp(int prlev,
                   int nrow,
                   int nrowc,
                   int ncol,
                   int iinv[],
                   int inew[],
                   matrix *at,
                   pivbas *basis,
                   int rsk[],
                   int csk[],
                   double rx[],
                   double cx[],
                   double xb[],
                   double yb[],
                   int rsety,
                   gsdec *lu,
                   int lresp,
                   int *refac,
                   int wisze,
                   int wimem[],
                   int wrsze,
                   double wrmem[])
{
  int resp;

  *refac = false;

  switch (lresp)
  {
  case lu_ok:
    return (ProcOk);
  case lu_refac:
  case lu_insta:
    *refac = true;
    resp = pivlp_factor2(prlev,
                         nrow, nrowc, ncol,
                         iinv, inew,
                         at,
                         basis,
                         rsk, csk,
                         rx, cx, xb,
                         yb,
                         rsety,
                         lu,
                         wisze, wimem,
                         wrsze, wrmem);

    break;
  default:
    printf("\n\n system error, unknown case.\n");
    ShutDown();
    exit(0);
  }
  return (resp);
} /* pivlp_prolresp */

void pivlp_swapbasis(int r1,
                     int r2,
                     pivbas *basis)
{
  int jtemp;
  int vttemp;

  vttemp = basis[r1].vt;
  jtemp = basis[r1].j;

  basis[r1].vt = basis[r2].vt;
  basis[r1].j = basis[r2].j;

  basis[r2].vt = vttemp;
  basis[r2].j = jtemp;
} /* pivlp_swapbasis */

void pivlp_mvartrow(int nrowc,
                    int iinv[],
                    pivbas *basis)
{
  int k, i, j;

  for (k = 0; k < nrowc;)
  {
    i = iinv[k];
    j = basis[i].j;
    if (basis[i].vt == PIV_ART || basis[i].vt == PIV_ROW)
    {
      if (i == j)
        ++k;
      else
        pivlp_swapbasis(i, j, basis);
    }
    else
      ++k;
  }
} /* pivlp_mvartrow */

void pivlp_formyb(int nrow,
                  int ncol,
                  matrix *at,
                  double cc[],
                  pivbas *basis,
                  double yb[],
                  gsdec *lu,
                  int wrsze,
                  double wrmem[])
{
  int i;
  double *res, *dy;

  if (wrsze < 2 * nrow)
  {
    printf("\n\n system error. work space size.\n");
    ShutDown();
    exit(0);
  }

  res = wrmem;
  dy = wrmem + nrow;

  for (i = 0; i < nrow; ++i)
  {
    switch (basis[i].vt)
    {
    case PIV_ART:
      res[i] = -yb[basis[i].j];
      break;
    case PIV_ROW:
      res[i] = yb[basis[i].j];
      break;
    case PIV_COL:
      res[i] = cc[basis[i].j] - svDot(at->ia + basis[i].j, yb);
      break;
    }
  }

  lu_btran(lu, res, dy);

  addVect(nrow, 1.0, dy, NULL, yb);
} /* pivlp_formyb */

double GetDualObj(int nrow,
                  int ncol,
                  double cf,
                  int rbk[],
                  double rbl[],
                  double rbu[],
                  int cbk[],
                  double cbl[],
                  double cbu[],
                  int rsk[],
                  int csk[],
                  double yb[],
                  double csb[])
{
  int i, j;
  double dobj = 0.0;

  for (i = 0; i < nrow; ++i)
  {
    if (optdef_getrbki(rbk, i) == FX)
      dobj += yb[i] * optdef_getrbli(rbk, rbl, i);
    else
    {
      switch (rsk[i])
      {
      case LOWER:
        dobj += yb[i] * optdef_getrbli(rbk, rbl, i);
        break;
      case UPPER:
        dobj -= yb[i] * optdef_getrbui(rbk, rbu, i);
        break;
      }
    }
  }

  for (j = 0; j < ncol; ++j)
  {
    switch (csk[j])
    {
    case LOWER:
      dobj -= csb[j] * optdef_getcblj(cbk, cbl, j);
      break;
    case UPPER:
      dobj += csb[j] * optdef_getcbuj(cbk, cbu, j);
      break;
    }
  }

  return (dobj + cf);
} /* GetDualObj */
