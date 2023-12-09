#include "LPcross.h"

static void setYX0(int n,
                   double alfa,
                   double *x,
                   int *s,
                   double *y)
{
  int i;

  if (s)
  {
    if (alfa == 0.0)
    {
      for (i = 0; i < n; ++i)
        x[s[i]] = 0.0;
    }
    else if (alfa == 1.0)
    {
      for (i = 0; i < n; ++i)
      {
        y[s[i]] += x[s[i]];
        x[s[i]] = 0.0;
      }
    }
    else if (alfa == -1.0)
    {
      for (i = 0; i < n; ++i)
      {
        y[s[i]] -= x[s[i]];
        x[s[i]] = 0.0;
      }
    }
    else
    {
      for (i = 0; i < n; ++i)
      {
        y[s[i]] += alfa * x[s[i]];
        x[s[i]] = 0.0;
      }
    }
  }
  else
  {
    ShutDown();
    exit(0);
  }
} /* setYX0 */

void prf_swapbasis(int r1,
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
} /* prf_swapbasis */

void prf_mvartrow(int nrowc,
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
        prf_swapbasis(i, j, basis);
    }
    else
      ++k;
  }
} /* prf_mvartrow */

static void prf_drtestc(prfparamt *param,
                        int nrow,
                        int ncol,
                        int orsk[],
                        int ocsk[],
                        int nnzj[],
                        double yb[],
                        double csb[],
                        int szeu,
                        double valu[],
                        int subu[],
                        int vtj,
                        int j,
                        double sj,
                        int *vti,
                        int *i,
                        double *theta)
{
  int decj, incj, t, q, imaxj, dmaxj, nilj = ncol;
  double dtheta = 0.0, itheta = 0.0, sq, aiq, abs_aiq, rtemp,
         tols, tolapiv;

  tols = param->tols;
  tolapiv = param->tolapiv;

  if (orsk)
  {
    printf("\n\n system error, orsk is defined.\n");
    ShutDown();
    exit(0);
  }

  /*
   * Check the columns.
   */

  decj = nilj;
  incj = nilj;

  for (t = 0; t < szeu; ++t)
  {
    q = subu[t];
    if (ocsk[q] != BASIC && ocsk[q] != FIXED)
    {
      aiq = valu[t];
      abs_aiq = fabs(aiq);

      if (abs_aiq > tolapiv)
      {
        switch (ocsk[q])
        {
        case LOWER:
        case UPPER:
          sq = csb[q];

          if (ocsk[q] == UPPER)
          {
            sq = -sq;
            aiq = -aiq;
          }

          if (sq < 0.0)
            sq = 0.0;

          if (aiq < 0.0)
          {
            rtemp = sq + aiq * itheta;
            if (incj == nilj || rtemp < -tols)
            {
              incj = q;
              itheta = -sq / aiq;
              imaxj = nnzj[q];
            }
            else if (rtemp <= tols && imaxj < nnzj[q])
              imaxj = nnzj[q];
          }
          else if (aiq > 0.0)
          {
            rtemp = sq - aiq * dtheta;
            if (decj == nilj || rtemp < -tols)
            {
              decj = q;
              dtheta = sq / aiq;
              dmaxj = nnzj[q];
            }
            else if (rtemp <= tols && dmaxj < nnzj[q])
              dmaxj = nnzj[q];
          }
          break;
        case SPBAS:
          if (incj == nilj || itheta != 0.0)
          {
            incj = q;
            itheta = 0.0;
          }

          if (decj == nilj || dtheta != 0.0)
          {
            decj = q;
            dtheta = 0.0;
          }
          break;
        }
      }
    }
  }

  switch (vtj)
  {
  case PIV_ART:
    if (incj == nilj && decj == nilj)
    {
      *vti = vtj;
      *i = j;
      *theta = -sj;
    }
    else if (decj == nilj || (incj != nilj && imaxj < dmaxj))
    {
      *vti = PIV_COL;
      *i = incj;
      *theta = itheta;
    }
    else
    {
      *vti = PIV_COL;
      *i = decj;
      *theta = -dtheta;
    }
    break;
  case PIV_ROW:
  case PIV_COL:
    printf("\n\n leaving is not an art var vtj=%d.\n", vtj);
    ShutDown();
    exit(0);
  }
} /* prf_drtestc */

static void prf_dualcrash(prfparamt *param,
                          int nrow,
                          int nrowc,
                          int ncol,
                          int iinv[],
                          int inew[],
                          matrix *a,
                          matrix *at,
                          int irsk[],
                          int icsk[],
                          pivbas *basis,
                          int orsk[],
                          int ocsk[],
                          double yb[],
                          double csb[],
                          prfinfot *prfinf,
                          int wisze,
                          int wimem[])
{
  int i, j, t, szeu, k, inc, itemp, sze, cur, sing, ncrash, nrowc0,

      *nnzi, *nnzj, *list, *fir, *link, *subu;
  int vtj, vtinc;
  double sj, theta,
      tols = param->tols,
      *valu;
  array *ai, *aj;

  /*
  if (param->prlev>=0)
    printf("     dual crash... (");
  */
  if (wisze < 1 + 2 * nrow + 3 * ncol)
  {
    printf("\n\n system error, work space.\n");
    ShutDown();
    exit(0);
  }

  nnzi = wimem;
  nnzj = nnzi + nrow;
  list = nnzj + ncol;
  fir = list + ncol;
  link = fir + (ncol + 1);

  /*
   * Count number of non-zeros in each row basis.
   */
  iZero(nrow, nnzi, NULL);
  for (k = 0; k < nrowc; ++k)
  {
    j = basis[iinv[k]].j;
    switch (basis[iinv[k]].vt)
    {
    case PIV_ROW:
      nnzi[j]++;
      break;
    case PIV_COL:
      aj = at->ia + j;
      for (t = 0; t < aj->nn0; ++t)
      {
        i = aj->ja[t];
        if (inew[i] < nrowc)
          ++nnzi[i];
      }
      break;
    }
  }

  /*
   * Move non-null rows to end.
   */
  for (k = 0; k < nrowc;)
  {
    i = iinv[k];
    if (nnzi[i] > 0 && basis[i].vt == PIV_ART)
    {
      --nrowc;

      iSwap(k, nrowc, iinv);

      inew[iinv[k]] = k;
      inew[iinv[nrowc]] = nrowc;
    }
    else
      ++k;
  }
  /*
  if (param->prlev>=0)
    printf(" nrows=%d",nrowc);
  */

  /*
   * Count the numbers of non-zeros in each row.
   */

  for (k = 0; k < nrowc; ++k)
  {
    i = iinv[k];
    nnzi[i] = 0;
    ai = a->ia + i;
    for (t = 0; t < ai->nn0; ++t)
    {
      j = ai->ja[t];
      if (ocsk[j] != FIXED)
        ++nnzi[i];
    }
  }

  for (k = nrowc; k < nrow; ++k)
    nnzi[iinv[k]] = ncol + 1;

  fSort(nrow, ncol, nnzi, fir, link);

  for (j = 0; j < ncol; ++j)
  {
    nnzj[j] = 0;
    if (ocsk[j] != BASIC && ocsk[j] != FIXED)
    {
      aj = at->ia + j;
      itemp = 0;
      for (t = 0; t < aj->nn0; ++t)
        if (inew[aj->ja[t]] < nrowc)
          ++itemp;

      nnzj[j] = itemp;
    }
  }

  cur = ncol;

  while (true)
  {
    if (fir[cur] != nrow)
    {
      i = fir[cur];
      k = inew[i];

      if (k >= nrowc)
        exit(0);

      i = iinv[k];

      j = basis[i].j;
      vtj = basis[i].vt;
      sj = -yb[i];

      ai = a->ia + i;

      szeu = ai->nn0;
      subu = ai->ja;
      valu = ai->an;

      prf_drtestc(param,
                  nrow, ncol,
                  orsk, ocsk,
                  nnzj,
                  yb, csb,
                  szeu, valu, subu,
                  vtj, j, sj,
                  &vtinc, &inc, &theta);

      yb[i] -= theta;

      sing = 0;
      for (t = 0; t < szeu; ++t)
      {
        csb[subu[t]] += theta * valu[t];
        sing |= (ocsk[j] != FIXED && nnzj[subu[t]] == 1 &&
                 fabs(csb[subu[t]]) <= tols);
      }

      if (sing)
        for (t = 0; t < szeu; ++t)
          --nnzj[subu[t]];

      fir[cur] = link[i];
    }
    else
    {
      if (cur == 0)
        break;
      else
        --cur;
    }
  }

  /*
   * Count the numbers of non-zeros in each row.
   */

  for (k = 0; k < nrowc; ++k)
  {
    i = iinv[k];
    nnzi[i] = 0;
    ai = a->ia + i;
    for (t = 0; t < ai->nn0; ++t)
    {
      j = ai->ja[t];
      if (ocsk[j] != FIXED)
        ++nnzi[i];
    }
  }

  for (k = nrowc; k < nrow; ++k)
    nnzi[iinv[k]] = ncol + 1;

  /*
   * Count number of non-zeros in each column.
   */

  sze = 0;
  for (j = 0; j < ncol; ++j)
  {
    nnzj[j] = 0;
    if (ocsk[j] != BASIC && ocsk[j] != FIXED)
    {
      aj = at->ia + j;
      itemp = 0;
      for (t = 0; t < aj->nn0; ++t)
        if (inew[aj->ja[t]] < nrowc)
          ++itemp;

      nnzj[j] = itemp;
      if (itemp == 1 && fabs(csb[j]) <= param->tols)
        list[sze++] = j;
    }
  }

  fSort(nrow, ncol, nnzi, fir, link);

  ncrash = 0;
  nrowc0 = nrowc;
  cur = ncol;
  do
  {
    while (sze)
    {
      --sze;

      j = list[sze];
      if (nnzj[j] == 1)
      {
        ++ncrash;

        aj = at->ia + j;

        for (t = 0; t < aj->nn0; ++t)
          if (inew[aj->ja[t]] < nrowc)
            break;

        if (t == aj->nn0)
        {
          /*
          printf("  j=%d  csk=%d",j,ocsk[j]);
          */
          for (k = 0; k < nrow; ++k)
            if (iinv[inew[k]] != k)
            {
              printf("\n\n system error, iinv wrong.\n");
              ShutDown();
              exit(0);
            }

          printf("\n\n system error, not found.\n");
          ShutDown();
          exit(0);
        }

        i = aj->ja[t];
        k = inew[i];

        if (basis[i].vt != PIV_ART)
        {
          printf("\n\n system error, leaving var is not art.\n");
          ShutDown();
          exit(0);
        }
        basis[i].vt = PIV_COL;
        basis[i].j = j;
        ocsk[j] = BASIC;

        --nrowc;

        iSwap(k, nrowc, iinv);

        inew[iinv[k]] = k;
        inew[iinv[nrowc]] = nrowc;

        ai = a->ia + i;

        for (t = 0; t < ai->nn0; ++t)
        {
          j = ai->ja[t];
          if (ocsk[j] != BASIC && ocsk[j] != FIXED)
          {
            --nnzj[j];

            if (nnzj[j] == 1 && fabs(csb[j]) <= param->tols)
              list[sze++] = j;
          }
        }
      }
    }

    for (; cur;)
    {
      if (fir[cur] != nrow)
      {
        i = fir[cur];
        k = inew[i];

        if (k < nrowc)
        {
          --nrowc;

          iSwap(k, nrowc, iinv);

          inew[iinv[k]] = k;
          inew[iinv[nrowc]] = nrowc;

          ai = a->ia + i;

          for (t = 0; t < ai->nn0; ++t)
          {
            j = ai->ja[t];
            if (ocsk[j] != BASIC && ocsk[j] != FIXED)
            {
              --nnzj[j];

              if (nnzj[j] == 1 && fabs(csb[j]) <= param->tols)
                list[sze++] = j;
            }
          }
          break;
        }

        fir[cur] = link[i];
      }
      else
        --cur;
    }
  } while (cur);
  /*
  if (param->prlev>=0)
    printf(",ncrash=%d) end.\n",ncrash);
  */
} /* prf_dualcrash */

static void prf_duared(prfparamt *param,
                       int nrow,
                       int *nrowc,
                       int ncol,
                       int iinv[],
                       int inew[],
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
  int blksze;
  /*
  if (param->prlev>=0)
  {
    printf("       Blkred.");
    printf(" nrow0: %-7d",*nrowc);
    printf(" rsze: %-7d",rsze);
  }
  */
  blksze = basord_blklwr(nrow, ncol, *nrowc,
                         inew, iinv,
                         rsze, rset,
                         a, at,
                         basis,
                         wisze, wimem,
                         wzisze, wzimem);

  *nrowc -= blksze;
  /*
  if (param->prlev>=0) {
    printf(" blsze: %-7d",blksze);
    printf(" nrow1: %-7d\n",*nrowc);
  }
  */
} /* prf_duared */

static void prf_insertsing(prfparamt *param,
                           int nrow,
                           int nrowc,
                           int ncol,
                           int iinv[],
                           int inew[],
                           matrix *a,
                           pivbas *basis,
                           int orsk[],
                           int ocsk[],
                           double yb[],
                           double csb[],
                           int wisze,
                           int wimem[])
{
  int i, j, k, t, count = 0, foundj, *nnzj;
  double rtemp, aij;
  array *ai;

  if (wisze < ncol)
  {
    printf("\n\n work space size error.\n");
    ShutDown();
    exit(0);
  }

  nnzj = wimem;

  iZero(ncol, nnzj, NULL);
  for (k = 0; k < nrowc; ++k)
  {
    i = iinv[k];
    ai = a->ia + i;
    for (t = 0; t < ai->nn0; ++t)
      ++nnzj[ai->ja[t]];
  }

  for (k = 0; k < nrowc; ++k)
  {
    if (basis[iinv[k]].vt == PIV_ART)
    {
      i = basis[iinv[k]].j;

      if (inew[i] >= nrowc)
        exit(0);

      foundj = ncol;

      ai = a->ia + i;
      for (t = 0; t < ai->nn0; ++t)
      {
        j = ai->ja[t];
        aij = ai->an[t];
        if (nnzj[j] == 1 &&
            ocsk[j] != FIXED &&
            fabs(aij) >= 1.0e-3)
        {
          if (foundj == ncol || fabs(csb[j] / aij) < rtemp)
          {
            foundj = t;
            rtemp = fabs(csb[j] / aij);
          }
        }
      }

      if (foundj != ncol)
      {
        j = ai->ja[foundj];

        count++;
        if (ocsk[j] == BASIC)
        {
          printf("\n\n system error, singular basis.\n");
          ShutDown();
          exit(0);
        }
        basis[iinv[k]].vt = PIV_COL;
        basis[iinv[k]].j = j;
        ocsk[j] = BASIC;
      }
    }
  }
  /*
  if (param->prlev>=0&&count)
    printf("\n     insert singletons: %-8d\n",count);
  */
} /* prf_insertsing */

static void prf_compress(int ncol,
                         int nrowc,
                         matrix *a,
                         int iinv[],
                         pivbas *basis,
                         int ocsk[],
                         int rsze,
                         int rset[],
                         int offset[],
                         int wzisze,
                         int wzimem[])
{
  int use_compress = true;
  int k, i, s, t, j,
      *finj;
  double rtemp;
  array *ai;

  if (wzisze < ncol)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  if (use_compress)
  {
    finj = wzimem;

    for (k = 0; k < nrowc; ++k)
    {
      i = iinv[k];
      j = basis[i].j;
      if (basis[i].vt == PIV_COL)
        finj[j] = 1;
    }

    for (k = 0; k < rsze; ++k)
    {
      i = rset[k];
      j = basis[i].j;
      if (basis[i].vt == PIV_COL)
        finj[j] = 0;
    }

    for (k = 0; k < nrowc; ++k)
    {
      i = iinv[k];
      ai = a->ia + i;
      s = 0;
      for (t = 0; t < ai->nn0; ++t)
      {
        j = ai->ja[t];
        if (finj[j])
        {
          rtemp = ai->an[t];
          ai->an[t] = ai->an[s];
          ai->an[s] = rtemp;

          ai->ja[t] = ai->ja[s];
          ai->ja[s] = j;
          ++s;
        }
      }
      offset[i] = s;
    }

    /*
     * Reset the buffer.
     */

    for (k = 0; k < nrowc; ++k)
    {
      i = iinv[k];
      j = basis[i].j;
      if (basis[i].vt == PIV_COL)
        finj[j] = 0;
    }
  }
  else
  {
    for (k = 0; k < nrowc; ++k)
      offset[iinv[k]] = 0;
  }
} /* prf_compress */

static void prf_drtest(prfparamt *param,
                       int nrow,
                       int ncol,
                       int irsk[],
                       int icsk[],
                       int orsk[],
                       int ocsk[],
                       double csb[],
                       int szeu,
                       double valu[],
                       int subu[],
                       int vtj,
                       int j,
                       double sj,
                       int *vti,
                       int *i,
                       double *theta)
{
  int decj, incj, t, q, nilj = ncol;
  double dtheta = 0.0, itheta = 0.0, sq, aiq, abs_aiq, abs_iaij, abs_daij, rtemp,
         tols, tolapiv;

  tols = param->tols;
  tolapiv = param->tolapiv;

  if (orsk)
  {
    printf("\n\n system error, orsk is defined.\n");
    ShutDown();
    exit(0);
  }

  /*
   * Check the columns.
   */
  decj = nilj;
  incj = nilj;

  for (t = 0; t < szeu; ++t)
  {
    q = subu[t];
    if (ocsk[q] != BASIC && ocsk[q] != FIXED)
    {
      aiq = valu[q];
      abs_aiq = fabs(aiq);

      if (abs_aiq > tolapiv)
      {
        switch (ocsk[q])
        {
        case LOWER:
        case UPPER:
          sq = csb[q];

          if (ocsk[q] == UPPER)
          {
            sq = -sq;
            aiq = -aiq;
          }

          if (sq < 0.0)
            sq = 0.0;

          if (aiq < 0.0)
          {
            rtemp = sq + aiq * itheta;
            if (incj == nilj || rtemp < -tols)
            {
              incj = q;
              itheta = -sq / aiq;
              abs_iaij = abs_aiq;
            }
            else if (rtemp <= tols && abs_aiq > abs_iaij)
            {
              incj = q;
              abs_iaij = abs_aiq;
            }
          }
          else if (aiq > 0.0)
          {
            rtemp = sq - aiq * dtheta;
            if (decj == nilj || rtemp < -tols)
            {
              decj = q;
              dtheta = sq / aiq;
              abs_daij = abs_aiq;
            }
            else if (rtemp <= tols && abs_aiq > abs_iaij)
            {
              decj = q;
              abs_daij = abs_aiq;
            }
          }
          break;
        case SPBAS:
          if (incj == nilj || itheta * aiq > tols || abs_aiq > abs_iaij)
          {
            incj = q;
            itheta = 0.0;
            abs_iaij = abs_aiq;
          }

          if (decj == nilj || dtheta * aiq < -tols || abs_aiq > abs_daij)
          {
            decj = q;
            dtheta = 0.0;
            abs_daij = abs_aiq;
          }
          break;
        }
      }
    }
  }

  switch (vtj)
  {
  case PIV_ART:
    if (incj == nilj && decj == nilj)
    {
      *vti = vtj;
      *i = j;
      *theta = -sj;
    }
    else if (decj == nilj || (incj != nilj && abs_iaij >= abs_daij))
    {
      *vti = PIV_COL;
      *i = incj;
      *theta = itheta;
    }
    else
    {
      *vti = PIV_COL;
      *i = decj;
      *theta = -dtheta;
    }
    break;
  case PIV_ROW:
    printf("\n\n system error, rows is not allowed.\n");
    ShutDown();
    exit(0);
    break;
  case PIV_COL:
    if (icsk[j] == BASIC || icsk[j] == SPBAS)
    {
      printf("\n\n system error, var type wrong.\n");
      ShutDown();
      exit(0);
    }

    if ((incj == nilj && decj == nilj) ||
        (icsk[j] == LOWER && (decj == nilj || sj - dtheta <= tols)) ||
        (icsk[j] == UPPER && (incj == nilj || sj + itheta >= -tols)))
    {
      *vti = vtj;
      *i = j;
      *theta = -sj;
    }
    else if (icsk[j] == LOWER)
    {
      *vti = PIV_COL;
      *i = decj;
      *theta = -dtheta;
    }
    else if (icsk[j] == UPPER)
    {
      *vti = PIV_COL;
      *i = incj;
      *theta = itheta;
    }
    else
    {
      printf("\n\n system error, ocsk wrong.\n");
      ShutDown();
      exit(0);
    }
    break;
  }

} /* prf_drtest */

static void prf_purdualp(prfparamt *param,
                         int nrow,
                         int *nrowc,
                         int ncol,
                         int iinv[],
                         int inew[],
                         int rsze,
                         int rset[],
                         matrix *a,
                         matrix *at,
                         int rbk[],
                         double rbl[],
                         double rbu[],
                         int cbk[],
                         double cbl[],
                         double cbu[],
                         int irsk[],
                         int icsk[],
                         double y[],
                         double rs[],
                         double cs[],
                         pivbas *basis,
                         int orsk[],
                         int ocsk[],
                         double orx[],
                         double ocx[],
                         double xb[],
                         double yb[],
                         double csb[],
                         gsdec *lu,
                         prfinfot *prfinf,
                         int iszero,
                         int wisze,
                         int wimem[],
                         int wzisze,
                         int wzimem[],
                         int wzrsze,
                         double wzrmem[])
{
  int exch, updyb = false, refac, i, j, t, r, szeu,
            k, l, q, inc, itemp, nnz, nnz1, nnz2, *offset,
            *sub1, *sub2, *subu, *map1, *map2, *mapu,
            *ajsub, *arsub, lresp = lu_ok, vtj, vtinc;
  double sj, rtemp, theta, nrmres, tols = param->tols,
                                   *val1, *val2, *valu, *ajval, *arval;
  array *ar, *aj;

  if (wisze < 2 + 5 * nrow + max(6 * nrow, 1 * ncol) ||
      wzisze < ncol || wzrsze < 2 * nrow + 1 * ncol)
  {
    printf("\n\n system error, work space.\n");
    ShutDown();
    exit(0);
  }

  if (!rsze)
    return;

  map1 = wimem;
  map2 = map1 + nrow;

  offset = map2 + nrow;

  sub1 = offset + nrow;
  sub2 = sub1 + nrow;
  subu = sub2 + nrow;

  val1 = wzrmem;
  val2 = val1 + nrow;
  valu = val2 + nrow;

  mapu = wzimem;

  if (iszero)
  {
    /*
     * Check if zero.
     */
  }
  else
  {
    iZero(*nrowc, map1, NULL);
    iZero(*nrowc, map2, NULL);
  }

  prf_compress(ncol, *nrowc,
               a,
               iinv,
               basis,
               ocsk,
               rsze, rset,
               offset,
               ncol, mapu);

  itemp = rsze;
  for (t = 0; t < itemp;)
  {
    if (basis[rset[t]].vt == PIV_ART)
      ++t;
    else
    {
      --itemp;
      iSwap(itemp, t, rset);
    }
  }

  for (t = rsze; t && prfinf->resp == ProcOk;)
  {
    --t;

    if (updyb)
    {
      for (k = 0; k < *nrowc; ++k)
      {
        r = iinv[k];

        j = basis[r].j;
        switch (basis[r].vt)
        {
        case PIV_ART:
          val1[k] = 0.0;
          break;
        case PIV_ROW:
          val1[k] = 0.0;
          break;
        case PIV_COL:
          val1[k] = svDot(at->ia + j, y) - csb[j] - svDot(at->ia + j, yb);

          switch (icsk[j])
          {
          case LOWER:
            val1[k] += max(0.0, cs[j]);
            break;
          case UPPER:
            val1[k] += min(0.0, cs[j]);
            break;
          }
          break;
        }
      }

      nrmres = dNorm0(*nrowc, val1, NULL);

      /*
      if (param->prlev>=0)
        printf("       Dres=%8.1e\n",nrmres);
       */
      if (nrmres >= 1.0e-3)
      {
        prfinf->resp = NumDiff;
        dZero(*nrowc, val1, NULL);
        break;
      }
      else if (nrmres >= 1.0e-13)
      {
        lu_btran(lu, val1, val2);

        for (k = 0; k < *nrowc; ++k)
          yb[iinv[k]] += val2[k];

        dCopy(nrow, y, val2);
        addVect(nrow, -1.0, yb, NULL, val2);

        for (j = 0; j < ncol; ++j)
        {
          aj = at->ia + j;
          csb[j] = svDot(aj, val2);

          switch (icsk[j])
          {
          case LOWER:
            csb[j] += max(0.0, cs[j]);
            break;
          case UPPER:
            csb[j] += min(0.0, cs[j]);
            break;
          }
        }

        dZero(nrow, val2, NULL);
      }
      else
        dZero(*nrowc, val1, NULL);

      updyb = false;
    }

    i = rset[t];
    j = basis[i].j;
    vtj = basis[i].vt;
    switch (vtj)
    {
    case PIV_ART:
      sj = -yb[j];
      break;
    case PIV_ROW:
      printf("\n\n system error, no row allowed.\n");
      ShutDown();
      exit(0);
    case PIV_COL:
      sj = csb[j];
      break;
    }

    if (fabs(sj) >= tols)
    {
      prfinf->ited++;

      itemp = inew[i];
      nnz2 = 1;
      val2[itemp] = 1.0;
      sub2[0] = itemp;
      map2[itemp] = 1;

      if (itemp >= *nrowc)
      {
        printf("\n\n system error, index itemp>=*nrowc.\n");
        ShutDown();
        exit(0);
      }

      lu_btran2(lu, nnz2, val2, sub2, map2, &nnz1, val1, sub1, map1);

      szeu = 0;
      if (nnz1 == 1)
      {
        r = sub1[0];
        rtemp = val1[r];
        r = iinv[r];

        ar = a->ia + r;

        nnz = ar->nn0;
        arsub = ar->ja;
        arval = ar->an;

        if (rtemp == 1.0)
        {
          for (l = offset[r]; l < nnz; ++l)
          {
            valu[arsub[l]] = arval[l];
            subu[szeu++] = arsub[l];
          }
        }

        else
        {
          for (l = offset[r]; l < nnz; ++l)
          {
            valu[arsub[l]] = rtemp * arval[l];
            subu[szeu++] = arsub[l];
          }
        }
      }
      else
      {
        for (k = 0; k < nnz1; ++k)
        {
          r = sub1[k];
          rtemp = val1[r];
          r = iinv[r];

          ar = a->ia + r;

          nnz = ar->nn0;
          arsub = ar->ja;
          arval = ar->an;

          for (l = offset[r]; l < nnz; ++l)
          {
            q = arsub[l];
            valu[q] += rtemp * arval[l];
            if (!mapu[q])
            {
              mapu[q] = 1;
              subu[szeu++] = q;
            }
          }
        }
        iZero(szeu, mapu, subu);
      }

      prf_drtest(param,
                 nrow, ncol,
                 irsk, icsk,
                 orsk, ocsk,
                 csb,
                 szeu, valu, subu,
                 vtj, j, sj,
                 &vtinc, &inc, &theta);

      exch = false;
      if (vtj != vtinc || j != inc)
      {
        exch = true;
        switch (basis[i].vt)
        {
        case PIV_ART:
          break;
        case PIV_ROW:
          if (irsk[basis[i].j] == SPBAS)
            exit(0);

          orsk[basis[i].j] = irsk[basis[i].j];
          break;
        case PIV_COL:
          if (icsk[basis[i].j] == SPBAS)
          {
            printf("\n\n system error, no superbasic is allowed.\n");
            ShutDown();
            exit(0);
          }
          ocsk[basis[i].j] = icsk[basis[i].j];

          break;
        default:
          exit(0);
        }

        switch (vtinc)
        {
        case PIV_ART:
          printf("\n\n system error, artcol can not be incomming.\n");
          ShutDown();
          exit(0);
        case PIV_ROW:
          basis[i].vt = PIV_ROW;
          basis[i].j = inc;
          orsk[inc] = BASIC;
          break;
        case PIV_COL:
          basis[i].vt = PIV_COL;
          basis[i].j = inc;
          ocsk[inc] = BASIC;

          valu[inc] = 0.0;
          csb[inc] = 0.0;
          break;
        }
      }

      if (theta != 0.0)
      {
        for (k = 0; k < nnz1; ++k)
        {
          r = sub1[k];
          yb[iinv[r]] -= theta * val1[r];

          val1[r] = 0.0;
          map1[r] = 0;
        }

        setYX0(szeu, theta, valu, subu, csb);
      }
      else
      {
        dZero(nnz1, val1, sub1);
        iZero(nnz1, map1, sub1);

        dZero(szeu, valu, subu);
      }

      if (exch)
      {
        switch (vtinc)
        {
        case PIV_ART:
          exit(0);
        case PIV_ROW:
          r = inew[inc];

          nnz1 = 1;
          val1[r] = 1.0;
          map1[r] = 1;
          sub1[0] = itemp;
          break;
        case PIV_COL:
          aj = at->ia + inc;
          ajsub = aj->ja;
          ajval = aj->an;
          nnz1 = 0;
          for (k = 0; k < aj->nn0; ++k)
          {
            r = inew[ajsub[k]];
            if (r < *nrowc)
            {
              val1[r] = ajval[k];
              sub1[nnz1] = r;
              map1[r] = 1;
              ++nnz1;
            }
          }
          break;
        }

        lresp = lu_supd(lu, inew[i], inew[i],
                        nnz1, val1, sub1, map1,
                        sub2, subu);

        if (lresp == lu_refac)
        {
          prf_duared(param,
                     nrow, nrowc, ncol,
                     iinv, inew,
                     t, rset,
                     a, at,
                     basis,
                     wisze, wimem,
                     ncol, mapu);
        }

        prfinf->resp = pivlp_prolresp(param->prlev,
                                      nrow, *nrowc, ncol,
                                      iinv, inew, at,
                                      basis,
                                      orsk, ocsk,
                                      orx, ocx, xb,
                                      yb, false,
                                      lu, lresp,
                                      &refac,
                                      min(wisze, 2 + 11 * (*nrowc)), wimem,
                                      min(wzrsze, 2 * (*nrowc)), wzrmem);

        if (prfinf->resp != ProcOk || refac)
        {
          dZero(min(wzrsze, 2 * (*nrowc)), wzrmem, NULL);
          iZero(*nrowc, map1, NULL);
          iZero(*nrowc, map2, NULL);
        }

        if (refac)
        {
          prf_compress(ncol, *nrowc,
                       a,
                       iinv,
                       basis, ocsk,
                       t, rset,
                       offset,
                       ncol, mapu);
        }
        updyb = refac;
      }
      else
        ++prfinf->mved;
    }
  }
} /* prf_purdualp */

void XcrDualMsg(prfparamt *param,
                int nrow,
                int ncol,
                double cc[],
                double cf,
                matrix *a,
                matrix *at,
                int rbk[],
                double rbl[],
                double rbu[],
                int cbk[],
                double cbl[],
                double cbu[],
                pivbas *basis,
                int orsk[],
                int ocsk[],
                double yb[],
                prfinfot *prfinf,
                int wrsze,
                double wrmem[])
{
  int j, k, ndbi = 0;
  double sj, dbi, sdbi = 0.0, mdbi = 0.0, dfeas = 0.0, dfeasi = 0.0, viol,
                  nrmcb = 0.0, nrmcbi = 0.0,
                  *csb;

  if (wrsze < ncol)
  {
    printf("\n\n system error, work space size.\n");
    ShutDown();
    exit(0);
  }

  csb = wrmem;

  dCopy(ncol, cc, csb);
  mTimesv(false, ncol, nrow, -1.0, a, yb, 1.0, csb);

  for (k = 0; k < nrow; ++k)
  {
    switch (basis[k].vt)
    {
    case PIV_ROW:
    case PIV_ART:
      viol = fabs(yb[basis[k].j]);
      break;
    case PIV_COL:
      viol = fabs(cc[basis[k].j] - svDot(at->ia + basis[k].j, yb));
      nrmcb += fabs(cc[basis[k].j]);
      nrmcbi = max(nrmcbi, fabs(cc[basis[k].j]));
      break;
    default:
      printf("\n\n system error, work space size.\n");
      ShutDown();
      exit(0);
    }
    dfeas += fabs(viol);
    dfeasi = max(viol, dfeasi);
  }

  prfinf->nrmcbi = nrmcbi;
  prfinf->dfeasi = dfeasi;

  for (j = 0; j < ncol; ++j)
  {
    sj = cc[j] - svDot(at->ia + j, yb);
    dbi = 0;
    switch (ocsk[j])
    {
    case LOWER:
      dbi = max(0.0, -sj);
      break;
    case UPPER:
      dbi = max(0.0, sj);
      break;
    case NUFRE:
    case BASIC:
    case SPBAS:
      dbi = fabs(sj);
      break;
    default:
      printf("\n\n system error, unknown var id.\n");
      ShutDown();
      exit(0);
      break;
    }

    if (dbi > param->tols)
      ++ndbi;

    mdbi = max(mdbi, dbi);
    sdbi += dbi;
  }

  prfinf->dobj = GetDualObj(nrow, ncol, cf, rbk,
                            rbl, rbu, cbk, cbl, cbu,
                            orsk, ocsk, yb, csb);

  prfinf->ndbi = ndbi;
  prfinf->mdbi = mdbi;
  prfinf->sdbi = sdbi;
  prfinf->dfeas1 = dfeas;
  prfinf->nrmy = dNorm1(nrow, yb);
  prfinf->nrmcb = nrmcb;
} /* XcrDualMsg */

void XcrDualProc(prfparamt *param,
                 int finf,
                 int nrow,
                 int ncol,
                 double cc[],
                 double cf,
                 matrix *a,
                 matrix *at,
                 int rbk[],
                 double rbl[],
                 double rbu[],
                 int cbk[],
                 double cbl[],
                 double cbu[],
                 int irsk[],
                 int icsk[],
                 double y[],
                 double rs[],
                 double cs[],
                 pivbas *basis,
                 int orsk[],
                 int ocsk[],
                 double orx[],
                 double ocx[],
                 double xb[],
                 double yb[],
                 gsdec *lu,
                 prfinfot *prfinf,
                 int wisze,
                 int wimem[],
                 int wrsze,
                 double wrmem[])
{
  int refac;
  int i, j, k, t, rsze, wiszel, wzisze, wzrszel, nrowc, itemp,
      *iinv, *inew, *rset, *nnzi, *wimeml, *wzimem;
  double tols = param->tols,
         *csb, *wzrmeml;
  array *aj;

  printf("   Begin dual phase...\n");

  wiszel = 3 * nrow + ncol;
  itemp = 2 + max(8 * nrow + max(6 * nrow, 1 * ncol), 2 * nrow + 3 * ncol);

  if (wisze < wiszel + itemp || wrsze < 2 * nrow + 2 * ncol)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  iinv = wimem;
  inew = iinv + nrow;
  rset = inew + nrow;

  wzisze = ncol;
  wzimem = rset + nrow;

  wimeml = wimem + wiszel;
  wiszel = wisze - wiszel;

  iZero(wzisze, wzimem, NULL);

  for (i = 0; i < nrow; ++i)
  {
    iinv[i] = i;
    inew[i] = i;
  }

  csb = wrmem;

  wzrszel = wrsze - ncol;
  wzrmeml = wrmem + ncol;

  dCopy(nrow, y, yb);
  for (j = 0; j < ncol; ++j)
  {
    aj = at->ia + j;
    csb[j] = 0.0;
    switch (icsk[j])
    {
    case LOWER:
      csb[j] += max(0.0, cs[j]);
      break;
    case UPPER:
      csb[j] += min(0.0, cs[j]);
      break;
    }
  }

  dZero(wzrszel, wzrmeml, NULL);

  /*
   * First pivot on NULL rows.
   */

  nnzi = wimeml;

  iZero(nrow, nnzi, NULL);

  for (k = 0; k < nrow; ++k)
  {
    j = basis[k].j;
    switch (basis[k].vt)
    {
    case PIV_ROW:
      ++nnzi[j];
      break;
    case PIV_COL:
      aj = at->ia + j;
      for (t = 0; t < aj->nn0; ++t)
        ++nnzi[aj->ja[t]];
      break;
    }
  }

  nrowc = 0;
  for (k = 0; k < nrow; ++k)
    if (nnzi[k] == 0)
      ++nrowc;

  refac = false;
  if (nrowc >= 98)
  {
    rsze = 0;
    nrowc = 0;
    for (k = 0; k < nrow; ++k)
    {
      if (nnzi[k] == 0)
      {
        rset[rsze++] = k;
        itemp = inew[k];

        iSwap(itemp, nrowc, iinv);

        inew[iinv[nrowc]] = nrowc;
        inew[iinv[itemp]] = itemp;

        ++nrowc;
      }
    }

    prf_mvartrow(nrow, iinv, basis);

    prf_dualcrash(param, nrow, nrowc, ncol,
                  iinv, inew,
                  a, at,
                  irsk, icsk,
                  basis, orsk, ocsk,
                  yb, csb,
                  prfinf,
                  wiszel, wimeml);

    refac = true;
  }

  nrowc = nrow;
  rsze = 0;
  for (i = 0; i < nrow; ++i)
  {
    switch (basis[i].vt)
    {
    case PIV_ART:
      if (fabs(yb[basis[i].j]) > tols)
        rset[rsze++] = i;
      break;
    case PIV_ROW:
      exit(0);
    case PIV_COL:
      j = basis[i].j;
      if (icsk[j] != SPBAS && fabs(csb[j]) > tols)
        rset[rsze++] = i;
    }
  }
  /*
  if (param->prlev>=0) {
    printf("     final problem...\n");
    printf("       Subpro.");
    printf(" nrow0: %-8d",nrowc);
    printf(" rsze: %-8d\n",rsze);
  }
  */
  if (prfinf->resp == ProcOk && rsze > 98)
  {
    prf_duared(param,
               nrow, &nrowc, ncol,
               iinv, inew,
               rsze, rset,
               a, at,
               basis,
               wiszel, wimeml,
               wzisze, wzimem);

    prf_insertsing(param,
                   nrow, nrowc, ncol,
                   iinv, inew,
                   a, basis, orsk, ocsk,
                   yb, csb, wiszel, wimeml);

    refac = true;
  }

  if (prfinf->resp == ProcOk && (nrowc < nrow || refac))
  {
    if (nrowc == nrow)
    {
      for (i = 0; i < nrow; ++i)
      {
        iinv[i] = i;
        inew[i] = i;
      }
    }

    prfinf->resp = pivlp_factor2(param->prlev,
                                 nrow, nrowc, ncol,
                                 iinv, inew,
                                 at,
                                 basis,
                                 orsk, ocsk,
                                 orx, ocx, xb,
                                 yb, false,
                                 lu,
                                 wiszel, wimeml,
                                 min(2 * nrowc, wzrszel), wzrmeml);

    dZero(min(2 * nrowc, wzrszel), wzrmeml, NULL);
  }

  if (prfinf->resp == ProcOk && rsze)
  {
    prf_purdualp(param,
                 nrow, &nrowc, ncol,
                 iinv, inew,
                 rsze, rset,
                 a, at,
                 rbk, rbl, rbu,
                 cbk, cbl, cbu,
                 irsk, icsk,
                 y, rs, cs,
                 basis,
                 orsk, ocsk,
                 orx, ocx, xb,
                 yb, csb,
                 lu, prfinf,
                 false,
                 wiszel, wimeml,
                 wzisze, wzimem,
                 wzrszel, wzrmeml);
  }

  if (prfinf->resp == ProcOk)
  {
    refac = (nrowc < nrow);
    for (k = 0; k < nrow; ++k)
    {
      refac |= (iinv[k] != k);

      iinv[k] = k;
      inew[k] = k;
    }

    if (refac)
    {
      prfinf->resp = pivlp_factor2(param->prlev,
                                   nrow, nrow, ncol,
                                   iinv, inew,
                                   at,
                                   basis,
                                   orsk, ocsk,
                                   orx, ocx, xb,
                                   NULL, false,
                                   lu,
                                   wiszel, wimeml,
                                   wzrszel, wzrmeml);
    }
  }
  /*
  if (param->prlev>=0)
    printf("     end.\n");
  */
  if (prfinf->resp == ProcOk)
  {
    if (finf)
    {
      pivlp_formyb(nrow, ncol, at, cc, basis, yb, lu, wrsze, wrmem);

      XcrDualMsg(param,
                 nrow, ncol,
                 cc, cf,
                 a, at,
                 rbk, rbl, rbu,
                 cbk, cbl, cbu,
                 basis, orsk, ocsk,
                 yb, prfinf,
                 wrsze, wrmem);
    }
  }
  prfinf->tited += prfinf->ited;
  prfinf->tmved += prfinf->mved;

  if (xb)
  {
    for (k = 0; k < nrow; ++k)
    {
      switch (basis[k].vt)
      {
      case PIV_ART:
        xb[k] = 0.0;
        break;
      case PIV_ROW:
        xb[k] = orx[basis[k].j];
        break;
      case PIV_COL:
        xb[k] = ocx[basis[k].j];
        break;
      }
    }
  }
  printf("   End dual phase.\n");
} /* XcrDualProc */
