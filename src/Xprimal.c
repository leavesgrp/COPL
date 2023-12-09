#include "LPcross.h"

//#define TIME_COUNT

static void BasisCrash(prfparamt *param,
                       int nrow,
                       int *nrowc,
                       int ncol,
                       int iinv[],
                       int inew[],
                       int rsze,
                       int rset[],
                       int csze,
                       int cset[],
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
                       gsdec *lu,
                       prfinfot *prfinf,
                       int wisze,
                       int wimem[],
                       int wrsze,
                       double wrmem[])
{
  int done;
  int i, j, s, t, try_, nnz, blksze, mis, maxnnz, k, itemp, sze, maxtrsze, rank, // Tian Xie 2016.04.27
      *bas, *nnzj, *list;
  gsmsg *linfo;
  int lresp;
  xlist *rnnz;
  double uii, temp_tolapiv, rtemp,
      tolapiv = 1.0e-3;
  array *ai, *aj;
#ifdef TIME_COUNT
  clock_t Tm;
  Tm = GetTime();
#endif

  if (wisze < max(2 + 3 * nrow + 1 * ncol, 4 * (*nrowc)))
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  blksze = blkord_maxlwrsqblk(at, *nrowc,
                              iinv, inew,
                              csze, cset, NULL,
                              &maxtrsze,
                              min(wisze, 4 * (*nrowc)), wimem);

  for (k = 0; k < blksze; ++k)
  {
    i = iinv[k];
    j = cset[k];
    basis[i].vt = PIV_COL;
    basis[i].j = j;
    ocsk[j] = BASIC;
  }

  nnzj = wimem;
  list = nnzj + ncol;

  sze = 0;
  for (k = 0; k < csze; ++k)
  {
    j = cset[k];
    if (ocsk[j] == SPBAS)
    {
      aj = at->ia + j;
      itemp = 0;
      for (t = 0; t < aj->nn0; ++t)
      {
        i = inew[aj->ja[t]];
        if (blksze <= i && i < *nrowc)
          ++itemp;
      }
      nnzj[j] = itemp;
      if (itemp == 1)
        list[sze++] = j;
    }
  }

  sze = 0;
  while (sze)
  {
    --sze;
    j = list[sze];
    if (nnzj[j] == 1)
    {
      aj = at->ia + j;

      for (t = 0; t < aj->nn0; ++t)
      {
        i = aj->ja[t];
        if (blksze <= inew[i] && inew[i] < *nrowc)
          break;
      }

      if (t == aj->nn0)
        exit(0);

      if (fabs(aj->an[t]) >= tolapiv)
      {
        if (basis[i].vt != PIV_ART)
        {
          printf("\nNot artif");
          exit(0);
        }

        basis[i].vt = PIV_COL;
        basis[i].j = j;
        ocsk[j] = BASIC;

        k = inew[i];

        iSwap(blksze, k, iinv);

        inew[iinv[k]] = k;
        inew[iinv[blksze]] = blksze;

        ++blksze;

        ai = a->ia + i;
        for (t = 0; t < ai->nn0; ++t)
        {
          j = ai->ja[t];
          if (ocsk[j] == SPBAS)
          {
            --nnzj[j];
            if (nnzj[j] == 1)
              list[sze++] = j;
          }
        }
      }
    }
  }

  rnnz = (xlist *)calloc(1, sizeof(xlist));
  if (!rnnz)
  {
    prfinf->resp = OutOfSpc;
    return;
  }

  if (!setXt(rnnz, nrow, ncol + 1, min(2 + 3 * nrow + ncol, wisze), wimem))
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  for (k = blksze; k < *nrowc; ++k)
  {
    i = iinv[k];

    if (basis[i].vt != PIV_ART)
    {
      printf("\nbasis diff art2");
      exit(0);
    }

    putXt(rnnz, i, 0);

    if (orsk && orsk[i] == SPBAS)
      plXt(rnnz, i, 1);
  }

  for (k = blksze; k < csze; ++k)
  {
    j = cset[k];
    if (ocsk[j] == SPBAS)
    {
      aj = at->ia + j;
      for (t = 0; t < aj->nn0; ++t)
      {
        i = aj->ja[t];
        if (blksze <= inew[i] && inew[i] < *nrowc)
          plXt(rnnz, i, 1);
      }
    }
  }

  while (infXt(rnnz))
  {
    getXt(rnnz, &i, &nnz);
    delXt(rnnz, i);

    if (nnz)
    {
      ai = a->ia + i;
      s = ncol;
      maxnnz = 0;

      for (t = 0; t < ai->nn0; ++t)
      {
        j = ai->ja[t];
        if (ocsk[j] == SPBAS && fabs(ai->an[t]) >= 1.0e-6)
        {
          aj = at->ia + j;
          nnz = aj->nn0;
          if (s == ncol || nnz < maxnnz)
          {
            s = t;
            maxnnz = nnz;
          }
        }
      }

      if (s != ncol)
      {
        j = ai->ja[s];

        basis[i].vt = PIV_COL;
        basis[i].j = j;
        ocsk[j] = BASIC;

        aj = at->ia + j;
        for (t = 0; t < aj->nn0; ++t)
        {
          i = aj->ja[t];
          if (rnnz->pval[i] != rnnz->idep &&
              blksze <= inew[i] &&
              inew[i] < *nrowc)
            miXt(rnnz, i, 1);
        }
      }
    }
  }

  XtFree(&rnnz);

  mis = 0;
  for (k = 0; k < *nrowc; ++k)
  {
    i = iinv[k];
    if (basis[i].vt == PIV_COL)
      ocsk[basis[i].j] = SPBAS;

    if (basis[i].vt == PIV_ART)
      ++mis;
  }

  if (nrow == *nrowc)
  {
    for (i = 0; i < nrow; ++i)
    {
      iinv[i] = i;
      inew[i] = i;
    }
  }

  temp_tolapiv = lu->tolapiv;
  lu->tolapiv = tolapiv;
  lu->dropl = true;
  lu->dropu = true;
#ifdef TIME_COUNT
  printf("BasicCrash-A: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
  lresp = pivlp_factor1(nrow, *nrowc,
                        iinv, inew,
                        at,
                        false,
                        basis, lu,
                        wisze, wimem,
                        wrsze, wrmem); // HOTSPOT
#ifdef TIME_COUNT
  printf("BasicCrash-B: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
  lu->tolapiv = temp_tolapiv;
  lu->dropl = false;
  lu->dropu = false;
#ifdef TIME_COUNT
  printf("lresp: %d\n", lresp);
#endif
  switch (lresp)
  {
  case lu_ok:
  case lu_insta:
    bas = wimem;

    for (k = 0; k < *nrowc; ++k)
    {
      i = iinv[k];
      if (basis[i].vt == PIV_ART)
        bas[i] = nrow + ncol;
      else if (basis[i].vt == PIV_ROW)
        bas[i] = basis[i].j;
      else
        bas[i] = basis[i].j + nrow;
    }

    mis = 0;
    for (k = 0; k < *nrowc; ++k)
    {
      lu_getuii(lu, k, &i, &uii);

      if (i == *nrowc)
      {
        basis[iinv[k]].vt = PIV_ART;
        basis[iinv[k]].j = iinv[k];
      }
      else
      {
        if (bas[iinv[i]] < nrow)
        {
          basis[iinv[k]].vt = PIV_ROW;
          basis[iinv[k]].j = bas[iinv[i]];
        }
        else if (bas[iinv[i]] < nrow + ncol)
        {
          basis[iinv[k]].vt = PIV_COL;
          basis[iinv[k]].j = bas[iinv[i]] - nrow;
        }
      }
    }

    for (k = 0; k < *nrowc; ++k)
    {
      i = iinv[k];
      switch (basis[i].vt)
      {
      case PIV_ROW:
        orsk[basis[i].j] = BASIC;
        break;
      case PIV_COL:
        ocsk[basis[i].j] = BASIC;
        break;
      }
    }

    linfo = &lu->info;
    itemp = linfo->nnzl + linfo->nnzu;
    itemp += 3 * nrow + itemp / 2;

    if (itemp < lu->wsze)
    {
      lresp = LufRenew(lu, itemp, NULL);
      if (lresp != lu_ok)
      {
        itemp = linfo->nnzl + linfo->nnzu + nrow;

        lresp = LufRenew(lu, itemp, NULL);
        if (lresp != lu_ok)
          prfinf->resp = OutOfSpc;
      }
    }

    for (try_ = 0, done = false; try_ < 2 && !done && prfinf->resp == ProcOk; ++try_)
    {
      lresp = pivlp_factor1(nrow, *nrowc,
                            iinv, inew,
                            at, true, basis,
                            lu,
                            wisze, wimem,
                            wrsze, wrmem); // HOTSPOT

      switch (lresp)
      {
      case lu_ok:
      case lu_insta:
        lu_getrank(lu, &rank);

        if (rank < *nrowc || lresp == lu_insta)
        {
          if (rank < *nrowc)
          {
            mis = pivlp_insart(lu, nrow, *nrowc, ncol, iinv, inew, false, basis,
                               orsk, ocsk,
                               NULL, NULL, NULL, NULL, wisze, wimem);
          }

          rtemp = lu->tolrpiv;
          if (rtemp >= 0.9)
            break;
          else
          {
            rtemp = min(0.9, 2 * rtemp);
            lu->tolrpiv = rtemp;
          }

          done = false;
        }
        else
          done = true;

        break;
      case lu_wfull:
        prfinf->resp = OutOfSpc;
        break;
      default:
        printf("\n\n gob_inibas.\n");
        ShutDown();
        exit(0);
      }
    }
    break;
  case lu_space:
  case lu_wfull:
    prfinf->resp = OutOfSpc;
    break;
  default:
    exit(0);
  }
#ifdef TIME_COUNT
  printf("BasicCrash-C: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
} /* BasisCrash */

static void GetPrimResid(int nrow,
                         int ncol,
                         matrix *at,
                         int rbk[],
                         double rbl[],
                         double rbu[],
                         int cbk[],
                         double cbl[],
                         double cbu[],
                         int irsk[],
                         int icsk[],
                         double irx[],
                         double icx[],
                         pivbas *basis,
                         int orsk[],
                         int ocsk[],
                         double orx[],
                         double ocx[],
                         double xb[],
                         double pres[])
{
  int i, j, t, nnz, *subi;
  double xj, *vali;
  array *aj;

  for (i = 0; i < nrow; ++i)
  {
    j = basis[i].j;
    pres[i] = 0.0;
    switch (basis[i].vt)
    {
    default:
      break;
    case PIV_ROW:
      orx[j] = xb[i];
      break;
    case PIV_COL:
      ocx[j] = xb[i];
      break;
    }
  }

  for (i = 0; i < nrow; ++i)
  {
    if (optdef_getrbki(rbk, i) != FX &&
        (irsk[i] != orsk[i] ||
         orsk[i] == BASIC ||
         orsk[i] == SPBAS))
    {
      switch (irsk[i])
      {
      case LOWER:
      case FIXED:
        xj = -optdef_getrbli(rbk, rbl, i);
        break;
      case UPPER:
        xj = -optdef_getrbui(rbk, rbu, i);
        break;
      case SPBAS:
        xj = -irx[i];
        break;
      default:
        exit(0);
      }
      pres[i] += xj + orx[i];
    }
  }

  for (j = 0; j < ncol; ++j)
  {
    if (ocsk[j] != icsk[j] ||
        ocsk[j] == BASIC ||
        ocsk[j] == SPBAS)
    {
      switch (icsk[j])
      {
      case LOWER:
      case FIXED:
        xj = optdef_getcblj(cbk, cbl, j);
        break;
      case UPPER:
        xj = optdef_getcbuj(cbk, cbu, j);
        break;
      case SPBAS:
        xj = icx[j];

        if (cmps_haslo(optdef_getcbkj(cbk, j)))
          xj = max(xj, optdef_getcblj(cbk, cbl, j));

        if (cmps_hasup(optdef_getcbkj(cbk, j)))
          xj = min(xj, optdef_getcbuj(cbk, cbu, j));

        break;
      case NUFRE:
        xj = 0.0;
        break;
      default:
        exit(0);
      }

      xj -= ocx[j];

      if (xj)
      {
        aj = at->ia + j;
        nnz = aj->nn0;
        subi = aj->ja;
        vali = aj->an;

        for (t = 0; t < nnz; ++t)
          pres[subi[t]] += xj * vali[t];
      }
    }
  }
} /* GetPrimResid */

static void GenPrimBasis(prfparamt *param,
                         int nrow,
                         int nrowc,
                         int ncol,
                         int iinv[],
                         int inew[],
                         matrix *at,
                         int rbk[],
                         double rbl[],
                         double rbu[],
                         int cbk[],
                         double cbl[],
                         double cbu[],
                         int irsk[],
                         int icsk[],
                         double irx[],
                         double icx[],
                         pivbas *basis,
                         int orsk[],
                         int ocsk[],
                         double orx[],
                         double ocx[],
                         double xb[],
                         gsdec *lu,
                         prfinfot *prfinf,
                         double r1nrow[],
                         double r2nrow[])
{
  int k;
  double pfeas0, rtemp, umin, umax,
      *val1, *val2;

  val1 = r1nrow;
  val2 = r2nrow;

  GetPrimResid(nrow, ncol,
               at,
               rbk, rbl, rbu,
               cbk, cbl, cbu,
               irsk, icsk,
               irx, icx,
               basis,
               orsk, ocsk,
               orx, ocx,
               xb,
               val1);

  pfeas0 = dNorm0(nrowc, val1, iinv);

  for (k = 0; k < nrowc; ++k)
    val2[k] = val1[iinv[k]];

  lu_ftran(lu, val2, val1);

  for (k = 0; k < nrowc; ++k)
    xb[iinv[k]] += val1[k];

  GetPrimResid(nrow, ncol,
               at,
               rbk, rbl, rbu,
               cbk, cbl, cbu,
               irsk, icsk,
               irx, icx,
               basis,
               orsk, ocsk,
               orx, ocx,
               xb,
               val1);

  rtemp = dNorm0(nrowc, val1, iinv);
  if (rtemp > 2.0 * pfeas0 && rtemp > 1.0e-3)
  {
    lu_getucond(lu, &umin, &umax);
    prfinf->resp = NumDiff;
  }
} /* GenPrimBasis */

static int prf_pchuzr(prfparamt *param,
                      int nrow,
                      int iinv[],
                      int rbk[],
                      double rbl[],
                      double rbu[],
                      int cbk[],
                      double cbl[],
                      double cbu[],
                      pivbas *basis,
                      double orx[],
                      double ocx[],
                      double xb[],
                      int jisrow,
                      int j,
                      int nnzaj,
                      double valaj[],
                      int subaj[],
                      int *r,
                      double *max_absaij,
                      double *dx)
{
  int ichs, dchs;
  int bki, bkj;
  int dr, ir, i, t, itemp;
  double ddeltax, ideltax, abs_daij, abs_iaij, aij, abs_aij,
      xbi, rtemp, bli, bui, tolx, tolapiv;

  tolx = param->tolx;
  tolapiv = param->tolapiv;
  bkj = optdef_getcbkj(cbk, j);

  *max_absaij = 0.0;
  *dx = 0.0;
  t = 0;

  dr = nrow;
  ir = nrow;

  ideltax = 0.0;
  ddeltax = 0.0;

  for (; t < nnzaj; ++t)
  {
    i = iinv[subaj[t]];
    aij = valaj[subaj[t]];
    abs_aij = fabs(aij);

    if (abs_aij > *max_absaij)
      *max_absaij = abs_aij;

    xbi = xb[i];

    switch (basis[i].vt)
    {
    case PIV_COL:
      itemp = basis[i].j;
      bki = optdef_getcbkj(cbk, itemp);
      bli = optdef_getcblj(cbk, cbl, itemp);
      bui = optdef_getcbuj(cbk, cbu, itemp);
      break;
    case PIV_ART:
      bki = FX;
      bli = 0.0;
      bui = 0.0;
      break;
    case PIV_ROW:
      bki = optdef_getrbki(rbk, basis[i].j);
      bli = optdef_getrbli(rbk, rbl, basis[i].j);
      bui = optdef_getrbui(rbk, rbu, basis[i].j);
      break;
    }

    if (aij > tolapiv)
    {
      if (cmps_haslo(bki))
      {
        rtemp = (xbi - aij * ideltax) - bli;
        if (ir == nrow || rtemp < -tolx)
        {
          ir = i;
          ideltax = max((xbi - bli) / aij, 0.0);
          ichs = chuzr_blockl;
          abs_iaij = abs_aij;
        }
        else if (rtemp <= tolx && abs_aij > abs_iaij)
        {
          ir = i;
          ichs = chuzr_blockl;
          abs_iaij = abs_aij;
        }
      }

      if (cmps_hasup(bki))
      {
        rtemp = bui - (xbi + aij * ddeltax);
        if (dr == nrow || rtemp < -tolx)
        {
          dr = i;
          dchs = chuzr_blocku;
          ddeltax = max((bui - xbi) / aij, 0.0);
          abs_daij = abs_aij;
        }
        else if (rtemp <= tolx && abs_aij > abs_daij)
        {
          dr = i;
          dchs = chuzr_blocku;
          abs_daij = abs_aij;
        }
      }
    }
    else if (aij < -tolapiv)
    {
      if (cmps_haslo(bki))
      {
        rtemp = (xbi + aij * ddeltax) - bli;
        if (dr == nrow || rtemp < -tolx)
        {
          dr = i;
          dchs = chuzr_blockl;
          ddeltax = max(-(xbi - bli) / aij, 0.0);
          abs_daij = abs_aij;
        }
        else if (rtemp <= tolx && abs_aij > abs_daij)
        {
          dr = i;
          dchs = chuzr_blockl;
          abs_daij = abs_aij;
        }
      }

      if (cmps_hasup(bki))
      {
        rtemp = bui - (xbi - aij * ideltax);
        if (ir == nrow || rtemp < -tolx)
        {
          ir = i;
          ichs = chuzr_blocku;
          ideltax = max(-(bui - xbi) / aij, 0.0);
          abs_iaij = abs_aij;
        }
        else if (rtemp <= tolx && abs_aij > abs_iaij)
        {
          ir = i;
          ichs = chuzr_blocku;
          abs_iaij = abs_aij;
        }
      }
    }
  }

  *r = nrow;
  if (jisrow)
  {
    printf("\n system error.\n");
    exit(0);
  }
  else
  {
    if (bkj == FR)
    {
      if (ocx[j] >= 0.0 && (dr == nrow || ocx[j] <= ddeltax + tolx))
      {
        *dx = -ocx[j];
        return (chuzr_low);
      }

      if (ocx[j] <= 0.0 && (ir == nrow || -ocx[j] <= ideltax + tolx))
      {
        *dx = -ocx[j];
        return (chuzr_upr);
      }
    }
    else
    {
      if (cmps_haslo(bkj) && dr != nrow)
      {
        if (ocx[j] - optdef_getcblj(cbk, cbl, j) <= ddeltax + tolx)
        {
          *dx = -max(ocx[j] - optdef_getcblj(cbk, cbl, j), 0.0);
          return (chuzr_low);
        }
      }

      if (cmps_hasup(bkj) && ir != nrow)
      {
        if (optdef_getcbuj(cbk, cbu, j) - ocx[j] <= ideltax + tolx)
        {
          *dx = max(optdef_getcbuj(cbk, cbu, j) - ocx[j], 0.0);
          return (chuzr_upr);
        }
      }

      if (ir == nrow && dr == nrow)
      {
        if (cmps_haslo(bkj))
        {
          *dx = -max(ocx[j] - optdef_getcblj(cbk, cbl, j), 0.0);
          return (chuzr_low);
        }
        else
        {
          *dx = max(optdef_getcbuj(cbk, cbu, j) - ocx[j], 0.0);
          return (chuzr_upr);
        }
      }
    }
  }

  if (dr != nrow && (ir == nrow || (abs_daij >= abs_iaij)))
  {
    *r = dr;
    *dx = -ddeltax;
    return (dchs);
  }
  else if (ir != nrow)
  {
    *r = ir;
    *dx = ideltax;
    return (ichs);
  }

  return (chuzr_unb);
} /* prf_pchuzr */

static void BiPrimalLp(prfparamt *param,
                       double safe,
                       int nrow,
                       int *nrowc,
                       int ncol,
                       int iinv[],
                       int inew[],
                       int rsze,
                       int rset[],
                       int csze,
                       int cset[],
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
                       double irx[],
                       double icx[],
                       pivbas *basis,
                       int orsk[],
                       int ocsk[],
                       double orx[],
                       double ocx[],
                       double xb[],
                       gsdec *lu,
                       prfinfot *prfinf,
                       int wisze,
                       int wimem[],
                       int wrsze,
                       double wrmem[])
{
  int bypass, exch, refac, updxb;
  int chs;
  int i, j, k, t, r, nnz1, nnz2, nnz3, itemp, lea,
      itep = 0, mvep = 0,
      *sub1, *sub2, *sub3, *map1, *map3, *ajsub;
  int lresp;
  double rtemp, dx, max_absaij, aij, tolrpiv = 1.0e-12,
                                     *val1, *val2, *val3, *ajval;
  array *aj;
#ifdef TIME_COUNT
  clock_t Tm;
  Tm = GetTime();
#endif

  if (wisze < 2 + 11 * nrow || wrsze < 3 * nrow)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  sub1 = wimem;
  sub2 = sub1 + nrow;
  sub3 = sub2 + nrow;

  map1 = sub3 + nrow;
  map3 = map1 + nrow;

  val1 = wrmem;
  val2 = val1 + nrow;
  val3 = val2 + nrow;

  updxb = false;

  iZero(*nrowc, map1, NULL);
  iZero(*nrowc, map3, NULL);

  if (!updxb)
  {
    dZero(*nrowc, val1, NULL);
    dZero(*nrowc, val2, NULL);
  }

  dZero(*nrowc, val3, NULL);

  for (t = csze; t && prfinf->resp == ProcOk;)
  {
    --t;
    if (t == 217)
      t = t;

    if (updxb)
    {
      GenPrimBasis(param,
                   nrow, *nrowc, ncol,
                   iinv, inew,
                   at,
                   rbk, rbl, rbu,
                   cbk, cbl, cbu,
                   irsk, icsk,
                   irx, icx,
                   basis,
                   orsk, ocsk,
                   orx, ocx,
                   xb,
                   lu,
                   prfinf,
                   val1, val2);

      dZero(*nrowc, val1, NULL);
      dZero(*nrowc, val2, NULL);

      updxb = false;
    }

    if (prfinf->resp != ProcOk)
      break;

    bypass = false;
    j = cset[t];

    if (ocsk[j] == SPBAS)
    {
      itep++;

      aj = at->ia + j;
      ajsub = aj->ja;
      ajval = aj->an;

      nnz1 = 0;
      for (k = 0; k < aj->nn0; ++k)
      {
        i = inew[ajsub[k]];
        if (i < *nrowc)
        {
          val1[i] = ajval[k];
          map1[i] = 1;
          sub1[nnz1++] = i;
        }
      }

      lu_ftranl(lu, &nnz1, val1, sub1, map1);

      nnz3 = nnz1;
      for (k = 0; k < nnz1; ++k)
      {
        itemp = sub1[k];
        rtemp = val1[itemp];

        val3[itemp] = rtemp;
        sub3[k] = itemp;
        map3[itemp] = 1;
      }

      lu_ftranu(lu, &nnz3, val3, sub3, map3, &nnz2, val2, sub2, NULL);

      chs = prf_pchuzr(param,
                       nrow, iinv,
                       rbk, rbl, rbu,
                       cbk, cbl, cbu,
                       basis,
                       orx, ocx, xb,
                       false, j, nnz2, val2, sub2,
                       &r, &max_absaij, &dx);

      exch = false;

      if ((chs == chuzr_low || chs == chuzr_upr) &&
          optdef_getcbkj(cbk, j) == FR)
      {
        ocsk[j] = NUFRE;
      }
      else if (chs == chuzr_low)
      {
        if (!cmps_haslo(optdef_getcbkj(cbk, j)))
        {
          printf("\n\n has no lower.\n");
          ShutDown();
          exit(0);
        }
        ocsk[j] = LOWER;
      }
      else if (chs == chuzr_upr)
      {
        if (!cmps_hasup(optdef_getcbkj(cbk, j)))
        {
          printf("\n\n has no upper.\n");
          ShutDown();
          exit(0);
        }
        ocsk[j] = UPPER;
      }
      else if (chs == chuzr_blockl || chs == chuzr_blocku)
      {
        aij = val2[inew[r]];
        if (fabs(aij) >= tolrpiv * max_absaij &&
            fabs(aij) >= safe * param->tolapiv)
          exch = true;
        else
          bypass = true;
      }
      else
      {
        if (optdef_getcbkj(cbk, j) != FR)
        {
          printf("\n\n is not FR.\n");
          ShutDown();
          exit(0);
        }
        ocsk[j] = NUFRE;
      }

      if (bypass)
      {
        dZero(nnz1, val1, sub1);
        iZero(nnz1, map1, sub1);

        dZero(nnz2, val2, sub2);
      }
      else
      {
        if (!exch)
          ++mvep;

        for (k = 0; k < nnz2; ++k)
        {
          itemp = sub2[k];
          xb[iinv[itemp]] -= dx * val2[itemp];

          val2[itemp] = 0.0;
        }
        ocx[j] += dx;

        if (exch)
        {
          lresp = lu_supdl(lu, inew[r], inew[r], nnz1, val1, sub1, map1,
                           sub2, sub3);

          if (lresp != lu_insta)
          {
            switch (basis[r].vt)
            {
            case PIV_ART:
              break;
            case PIV_ROW:
              printf("\nRow in basis");
              exit(0);
            case PIV_COL:
              lea = basis[r].j;
              ocx[lea] = xb[r];

              if (chs == chuzr_blockl)
                ocsk[lea] = LOWER;
              else if (chs == chuzr_blocku)
                ocsk[lea] = UPPER;
              else
              {
                printf("\n\n system error.\n");
                ShutDown();
                exit(0);
              }
              break;
            }

            basis[r].vt = PIV_COL;
            basis[r].j = j;
            xb[r] = ocx[j];
            ocsk[j] = BASIC;
          }

          prfinf->resp = pivlp_prolresp(param->prlev,
                                        nrow, *nrowc, ncol,
                                        iinv, inew, at,
                                        basis,
                                        orsk, ocsk,
                                        orx, ocx, xb,
                                        NULL, false,
                                        lu, lresp,
                                        &refac,
                                        min(wisze, 2 + 11 * (*nrowc)), wimem,
                                        min(wrsze, 2 * (*nrowc)), wrmem);

          if (refac)
          {
            iZero(*nrowc, map1, NULL);
            iZero(*nrowc, map3, NULL);
            updxb = true;
          }
        }
        else
        {
          dZero(nnz1, val1, sub1);
          iZero(nnz1, map1, sub1);
        }
      }
    }
  }
#ifdef TIME_COUNT
  printf("BiPrimalLp: %d ms, itep = %d\n", GetTime() - Tm, itep);
#endif
  prfinf->itep += itep;
  prfinf->mvep += mvep;

  prfinf->titep += itep;
  prfinf->tmvep += mvep;
} /* BiPrimalLp */

void XcrPrimMsg(prfparamt *param,
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
                double orx[],
                double ocx[],
                double xb[],
                prfinfot *prfinf,
                int wrsze,
                double wrmem[])
{
  int i, j, k, npbi = 0;
  double xj, pbi, mpbi = 0.0, spbi = 0.0, pobj,
                  *bbar;

  if (wrsze < nrow)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }
  bbar = wrmem;

  dZero(nrow, bbar, NULL);

  pobj = cf;

  for (i = 0; i < nrow; ++i)
  {
    if (optdef_getrbki(rbk, i) == FX)
      xj = optdef_getrbli(rbk, rbl, i);
    else
    {
      switch (orsk[i])
      {
      case LOWER:
      case FIXED:
        xj = optdef_getrbli(rbk, rbl, i);
        break;
      case UPPER:
        xj = optdef_getrbui(rbk, rbu, i);
        break;
      case NUFRE:
      case SPBAS:
      case BASIC:
        xj = 0.0;
        break;
      default:
        exit(0);
      }
    }

    bbar[i] += xj;
  }

  for (j = 0; j < ncol; ++j)
  {
    switch (ocsk[j])
    {
    case LOWER:
    case FIXED:
      if (!cmps_haslo(optdef_getcbkj(cbk, j)))
        exit(0);

      xj = optdef_getcblj(cbk, cbl, j);
      break;
    case UPPER:
      if (!cmps_hasup(optdef_getcbkj(cbk, j)))
        exit(0);

      xj = optdef_getcbuj(cbk, cbu, j);
      break;
    case NUFRE:
      xj = 0.0;
      break;
    case SPBAS:
      printf("\n\n no superbasic is allowed.\n");
      ShutDown();
      exit(0);
    case BASIC:
      xj = 0.0;
      break;
    default:
      exit(0);
    }

    if (ocsk[j] != BASIC && fabs(ocx[j] - xj) > 1.0e-6)
    {
      printf("\n\n incorrect ocx compared to (status keys).\n");
      ShutDown();
      exit(0);
    }

    pobj += cc[j] * xj;

    setArray(-xj, at->ia + j, bbar);
  }

  prfinf->nrmb = dNorm1(nrow, bbar);
  prfinf->nrmbi = dNorm0(nrow, bbar, NULL);

  for (k = 0; k < nrow; ++k)
  {
    pbi = 0.0;
    switch (basis[k].vt)
    {
    case PIV_ART:
      pbi = fabs(xb[k]);
      break;
    case PIV_ROW:
      pbi = CompBndInf(xb[k],
                       optdef_getrbki(rbk, basis[k].j),
                       optdef_getrbli(rbk, rbl, basis[k].j),
                       optdef_getrbui(rbk, rbu, basis[k].j));
      break;
    case PIV_COL:
      j = basis[k].j;

      if (j < 0 || j >= ncol)
      {
        printf("\n\n wrong basis.j.\n");
        ShutDown();
        exit(0);
      }

      pobj += cc[j] * xb[k];
      pbi = CompBndInf(xb[k],
                       optdef_getcbkj(cbk, j),
                       optdef_getcblj(cbk, cbl, j),
                       optdef_getcbuj(cbk, cbu, j));

      if (ocsk[j] != BASIC)
      {
        printf("\n\n not basic=%d\n", j);
        ShutDown();
        exit(0);
      }
      break;
    default:
      exit(0);
    }

    if (pbi > param->tolx)
      ++npbi;

    mpbi = max(mpbi, pbi);
    spbi += pbi;

    if (basis[k].vt == PIV_ART)
      bbar[basis[k].j] -= xb[k];
    else if (basis[k].vt == PIV_ROW)
      bbar[basis[k].j] += xb[k];
    else
    {
      setArray(-xb[k], at->ia + basis[k].j, bbar);
      if (fabs(xb[k] - ocx[basis[k].j]) > 1.0e-6)
      {
        printf("\nocx 1");
        exit(0);
      }
    }
  }

  prfinf->pobj = pobj;
  prfinf->npbi = npbi;
  prfinf->mpbi = mpbi;
  prfinf->spbi = spbi;
  prfinf->pfeas1 = dNorm1(nrow, bbar);
  prfinf->pfeasi = dNorm0(nrow, bbar, NULL);
  prfinf->nrmxb = dNorm1(nrow, xb);
} /* XcrPrimMsg */

void XcrPrimProc(prfparamt *param,
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
                 double irx[],
                 double icx[],
                 pivbas *basis,
                 int orsk[],
                 int ocsk[],
                 double orx[],
                 double ocx[],
                 double xb[],
                 gsdec *lu,
                 prfinfot *prfinf,
                 int wisze,
                 int wimem[],
                 int wrsze,
                 double wrmem[])
{
  int refac, init_basis = true;
  int cutoff = 398;
  int i, j, k, rsze, csze, wiszel, nrowc, itemp,
      *iinv, *inew, *rset, *cset, *nnzj, *wimeml;
  double safe = 100.0;
#ifdef TIME_COUNT
  clock_t Tm;
  Tm = GetTime();
#endif

  printf("   Begin primal phase...\n");

  prfinf->itep = 0;
  prfinf->mvep = 0;

  wiszel = 2 * nrow + 2 * ncol;
  if (rbk)
    wiszel += nrow;

  itemp = 2 + max(11 * nrow, 3 * nrow + 3 * ncol);
  if (wisze < wiszel + itemp || wrsze < 3 * nrow)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  iinv = wimem;
  inew = iinv + nrow;

  if (rbk)
    rset = inew + nrow;
  else
    rset = inew;

  cset = rset + nrow;
  nnzj = cset + ncol;

  wiszel = wisze - wiszel;
  wimeml = nnzj + ncol;

  for (i = 0; i < nrow; ++i)
  {
    iinv[i] = i;
    inew[i] = i;
  }

  rsze = 0;
  csze = 0;

  if (rbk)
  {
    for (i = 0; i < nrow; ++i)
      if (orsk[i] == SPBAS)
        rset[rsze++] = i;
  }

  for (j = 0; j < ncol; ++j)
    if (ocsk[j] == SPBAS)
      cset[csze++] = j;

  nrowc = nrow;
#ifdef TIME_COUNT
  printf("  PRIMAL-START: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
  if (init_basis)
  {
    for (i = 0; i < nrow; ++i)
    {
      basis[i].vt = PIV_ART;
      basis[i].j = i;
      xb[i] = 0.0;
    }

    BasisCrash(param,
               nrow, &nrowc, ncol,
               iinv, inew,
               rsze, rset,
               csze, cset,
               a, at,
               rbk, rbl, rbu,
               cbk, cbl, cbu,
               basis, orsk, ocsk,
               lu, prfinf,
               wiszel, wimeml,
               wrsze, wrmem);
  }
#ifdef TIME_COUNT
  printf("  PRIMAL-A: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif

  if (prfinf->resp == ProcOk)
  {
    refac = false;

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

  rsze = 0;
  csze = 0;

  if (rbk)
  {
    for (i = 0; i < nrow; ++i)
      if (orsk[i] == SPBAS)
        rset[rsze++] = i;
  }

  for (j = 0; j < ncol; ++j)
    if (ocsk[j] == SPBAS)
      cset[csze++] = j;

  if (prfinf->resp == ProcOk)
  {
    csze = 0;
    for (j = 0; j < ncol; ++j)
      if (ocsk[j] == SPBAS)
        cset[csze++] = j;

    if (csze > 0)
    {
      if (refac)
      {
        prfinf->resp = pivlp_factor2(param->prlev,
                                     nrow, nrowc, ncol,
                                     iinv, inew,
                                     at,
                                     basis,
                                     orsk, ocsk,
                                     orx, ocx, xb,
                                     NULL, false,
                                     lu,
                                     wiszel, wimeml,
                                     wrsze, wrmem);
      }

      if (prfinf->resp == ProcOk)
      {
        BiPrimalLp(param,
                   safe,
                   nrow, &nrowc, ncol,
                   iinv, inew,
                   rsze, rset,
                   csze, cset,
                   a, at,
                   rbk, rbl, rbu,
                   cbk, cbl, cbu,
                   irsk, icsk,
                   irx, icx,
                   basis, orsk, ocsk,
                   orx, ocx, xb,
                   lu, prfinf,
                   wiszel, wimeml,
                   wrsze, wrmem);
      }
    }
  }
#ifdef TIME_COUNT
  printf("  PRIMAL-B: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
  if (prfinf->resp == ProcOk)
  {
    safe = 0.0;

    csze = 0;
    for (j = 0; j < ncol; ++j)
      if (ocsk[j] == SPBAS)
        cset[csze++] = j;

    if (csze > 0)
    {
      /*
      printf("   primal cleaning up...");
       */
      nrowc = nrow;
      for (k = 0; k < nrow; ++k)
      {
        refac |= (iinv[k] != k);

        iinv[k] = k;
        inew[k] = k;
      }
      prfinf->resp = pivlp_factor2(param->prlev,
                                   nrow, nrowc, ncol,
                                   iinv, inew,
                                   at,
                                   basis,
                                   orsk, ocsk,
                                   orx, ocx, xb,
                                   NULL, false,
                                   lu,
                                   wiszel, wimeml,
                                   wrsze, wrmem);

      if (prfinf->resp == ProcOk)
        BiPrimalLp(param,
                   safe,
                   nrow, &nrowc, ncol,
                   iinv, inew,
                   rsze, rset,
                   csze, cset,
                   a, at,
                   rbk, rbl, rbu,
                   cbk, cbl, cbu,
                   irsk, icsk,
                   irx, icx,
                   basis, orsk, ocsk,
                   orx, ocx, xb,
                   lu, prfinf,
                   wiszel, wimeml,
                   wrsze, wrmem);

      for (j = 0; j < ncol; ++j)
      {
        if (ocsk[j] == SPBAS)
        {
          if (cmps_hasbo(optdef_getcbkj(cbk, j)))
          {
            if (ocx[j] - optdef_getcblj(cbk, cbl, j) <= optdef_getcbuj(cbk, cbu, j) - ocx[j])
            {
              ocsk[j] = LOWER;
              ocx[j] = optdef_getcblj(cbk, cbl, j);
            }
            else
            {
              ocsk[j] = UPPER;
              ocx[j] = optdef_getcbuj(cbk, cbu, j);
            }
          }
          else if (cmps_haslo(optdef_getcbkj(cbk, j)))
          {
            ocsk[j] = LOWER;
            ocx[j] = optdef_getcblj(cbk, cbl, j);
          }
          else if (cmps_hasup(optdef_getcbkj(cbk, j)))
          {
            ocsk[j] = UPPER;
            ocx[j] = optdef_getcbuj(cbk, cbu, j);
          }
          else
          {
            ocsk[j] = NUFRE;
            ocx[j] = 0.0;
          }
        }
      }
    }
  }
#ifdef TIME_COUNT
  printf("  PRIMAL-C: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
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
                                   wrsze, wrmem);
    }
  }
#ifdef TIME_COUNT
  printf("  PRIMAL-D: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
  for (j = 0; j < ncol; ++j)
  {
    if (ocsk[j] == LOWER)
      ocx[j] = optdef_getcblj(cbk, cbl, j);
    else if (ocsk[j] == UPPER)
      ocx[j] = optdef_getcbuj(cbk, cbu, j);
    else if (ocsk[j] == NUFRE)
      ocx[j] = 0.0;
    else if (ocsk[j] == SPBAS)
    {
      if (cmps_hasbo(optdef_getcbkj(cbk, j)))
      {
        if (ocx[j] - optdef_getcblj(cbk, cbl, j) <=
            optdef_getcbuj(cbk, cbu, j) - ocx[j])
        {
          ocsk[j] = LOWER;
          ocx[j] = optdef_getcblj(cbk, cbl, j);
        }
        else
        {
          ocsk[j] = UPPER;
          ocx[j] = optdef_getcbuj(cbk, cbl, j);
        }
      }
      else if (cmps_haslo(optdef_getcbkj(cbk, j)))
      {
        ocsk[j] = LOWER;
        ocx[j] = optdef_getcblj(cbk, cbl, j);
      }
      else if (cmps_hasup(optdef_getcbkj(cbk, j)))
      {
        ocsk[j] = UPPER;
        ocx[j] = optdef_getcbuj(cbk, cbl, j);
      }
      else
      {
        ocsk[j] = NUFRE;
        ocx[j] = 0.0;
      }
    }
  }
#ifdef TIME_COUNT
  printf("  PRIMAL-E: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
  if (prfinf->resp == ProcOk)
  {
    pivlp_xb2rxcx(nrow, basis, orx, ocx, xb);

    if (finf)
    {
      pivlp_formxb(nrow, ncol, at, rbk, rbl, rbu, basis,
                   orx, ocx, xb,
                   lu,
                   wrsze, wrmem);

      XcrPrimMsg(param,
                 nrow, ncol,
                 cc, cf, a, at,
                 rbk, rbl, rbu,
                 cbk, cbl, cbu,
                 basis,
                 orsk, ocsk,
                 orx, ocx, xb,
                 prfinf,
                 wrsze, wrmem);
    }
  }

  printf("   End primal phase.\n");
#ifdef TIME_COUNT
  printf("  PRIMAL-F: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
} /* XcrPrimProc */
