#include "LPcross.h"

static void pivlp_setup(pivlppart *param,
                        pivlpt *pdat,
                        pivlpsolt *sdat,
                        pivlpslvit *slvi,
                        int rsze,
                        int rset[],
                        int invrset[],
                        int csze,
                        int cset[],
                        int invcset[],
                        iteinft *iinf)
{
  int *rsk = sdat->rsk, *csk = sdat->csk;
  int k, t;

  slvi->sresp = ProcOk;

  if (param->wlev == 0)
  {
    sdat->solsta = SOL_INFEAS;
    slvi->iter = 0;
    slvi->mves = 0;
  }

  t = 0;
  for (k = 0; k < rsze; ++k)
  {
    if (rsk[rset[k]] == BASIC)
    {
      iSwap(t, k, rset);

      invrset[rset[k]] = k;
      invrset[rset[t]] = t;
      t++;
    }
  }

  iinf->nar = rsze - t;
  iinf->nrp = rsze - t;

  t = 0;
  for (k = 0; k < csze; ++k)
  {
    if (csk[cset[k]] == BASIC)
    {
      iSwap(t, k, cset);

      invcset[cset[k]] = k;
      invcset[cset[t]] = t;
      t++;
    }
  }

  iinf->nac = csze - t;
  iinf->ncp = 0;

  if (sdat->solsta == SOL_OPTIMAL)
    sdat->solsta = SOL_FEAS;

  if (sdat->solsta == SOL_INFEAS)
    param->w12 = 0.0;
  else
    param->w12 = 1.0;

} /* pivlp_setup */

static void pivlp_formsubi(pivlppart *param,
                           pivlpt *pdat,
                           pivlpsolt *sdat,
                           pivlpslvit *slvi,
                           int *numi,
                           int subi[])
{
  int bk,
      *rbk = pdat->rbk, *cbk = pdat->cbk;
  int i, j, k,
      nrow = pdat->nrow;
  pivbas *basis = sdat->basis;
  double xj, l, u, spbi = 0.0, mpbi = 0.0, pbi,
                   tolx = param->tolx,
                   *rbl = pdat->rbl, *rbu = pdat->rbu,
                   *cbl = pdat->cbl, *cbu = pdat->cbu;

  *numi = 0;
  for (k = 0; k < nrow; ++k)
  {
    xj = sdat->xb[k];
    switch (sdat->basis[k].vt)
    {
    case PIV_ART:
      bk = FX;
      l = 0.0;
      u = 0.0;
      break;
    case PIV_ROW:
      i = basis[k].j;
      bk = optdef_getrbki(rbk, i);
      l = optdef_getrbli(rbk, rbl, i),
      u = optdef_getrbui(rbk, rbu, i);
      break;
    case PIV_COL:
      j = basis[k].j;
      bk = optdef_getcbkj(cbk, j);
      l = optdef_getcblj(cbk, cbl, j);
      u = optdef_getcbuj(cbk, cbu, j);
      break;
    }

    if (cmps_haslo(bk))
    {
      pbi = max(l - xj, 0.0);

      if (pbi > mpbi)
        mpbi = pbi;

      spbi += pbi;

      if (pbi > tolx)
      {
        subi[*numi] = k;
        ++*numi;
      }
    }

    if (cmps_hasup(bk))
    {
      pbi = max(xj - u, 0.0);

      if (pbi > mpbi)
        mpbi = pbi;

      spbi += pbi;

      if (pbi > tolx)
      {
        subi[*numi] = k;
        ++*numi;
      }
    }
  }

  slvi->mpbi = mpbi;
  slvi->spbi = spbi;
  slvi->npbi = *numi;

  if (numi == 0)
    sdat->solsta = SOL_FEAS;
} /* pivlp_formsubi */

static void pivlp_formcbs(pivlppart *param,
                          pivlpt *pdat,
                          pivlpsolt *sdat,
                          pivlpslvit *slvi,
                          int *nnz,
                          double cbs[],
                          int numi,
                          int subi[],
                          int sub[],
                          int map[])
{
  int bk,
      *rbk = pdat->rbk, *cbk = pdat->cbk;
  int i, j, k, npbi = 0;
  pivbas *basis = sdat->basis;
  double xj, l, u, spbi = 0.0, mpbi = 0.0, pbi,
                   tolx = param->tolx,
                   *rbl = pdat->rbl, *rbu = pdat->rbu,
                   *cbl = pdat->cbl, *cbu = pdat->cbu;

  for (i = 0; i < numi; ++i)
  {
    k = subi[i];
    xj = sdat->xb[k];
    switch (sdat->basis[k].vt)
    {
    case PIV_ART:
      bk = FX;
      l = 0.0;
      u = 0.0;
      break;
    case PIV_ROW:
      i = basis[k].j;
      bk = optdef_getrbki(rbk, i);
      l = optdef_getrbli(rbk, rbl, i),
      u = optdef_getrbui(rbk, rbu, i);
      break;
    case PIV_COL:
      j = basis[k].j;
      bk = optdef_getcbkj(cbk, j);
      l = optdef_getcblj(cbk, cbl, j);
      u = optdef_getcbuj(cbk, cbu, j);
      break;
    }

    if (cmps_haslo(bk))
    {
      pbi = max(l - xj, 0.0);

      if (pbi > mpbi)
        mpbi = pbi;

      spbi += pbi;

      if (pbi > tolx)
      {
        cbs[k] = -1.0;
        sub[npbi] = k;
        map[k] = 1;

        npbi++;
      }
    }

    if (cmps_hasup(bk))
    {
      pbi = max(xj - u, 0.0);

      if (pbi > mpbi)
        mpbi = pbi;

      spbi += pbi;

      if (pbi > tolx)
      {
        cbs[k] = 1.0;
        sub[npbi] = k;
        map[k] = 1;

        npbi++;
      }
    }
  }

  *nnz = npbi;

  slvi->mpbi = mpbi;
  slvi->spbi = spbi;
  slvi->npbi = npbi;

  if (npbi == 0)
    sdat->solsta = SOL_FEAS;
} /* pivlp_formcbs */

static void pivlp_formy(pivlppart *param,
                        pivlpt *pdat,
                        pivlpsolt *sdat,
                        pivlpslvit *slvi,
                        int numi,
                        int subi[],
                        gsdec *lu,
                        int wisze,
                        int wimem[],
                        int wiszez,
                        int wimemz[],
                        int wrszez,
                        double wrmemz[])
{
  int nnz1, nnz2,
      nrow = pdat->nrow, ncol = pdat->ncol,
      *sub1, *sub2, *map1, *map2;
  double *val1;

  if (wisze < 2 * nrow || wiszez < 2 * nrow || wrszez < 1 * nrow)
  {
    printf("\n\n system error, work space size.\n");
    ShutDown();
    exit(0);
  }

  if (sdat->solsta == SOL_INFEAS)
  {
    /*
     * Solution is infeasible.
     */

    nnz1 = 0;
    nnz2 = 0;

    sub1 = wimem;
    sub2 = sub1 + nrow;

    map1 = wimemz;
    map2 = map1 + nrow;

    val1 = wrmemz;

    pivlp_formcbs(param, pdat, sdat, slvi,
                  &nnz1, val1,
                  numi, subi,
                  sub1, map1);

    if (sdat->solsta == SOL_INFEAS)
    {
      /*
       * y must be zeroed before use
       * of lu_btran2.
       */

      dZero(nrow, sdat->y, NULL);

      lu_btran2(lu,
                nnz1, val1, sub1, map1,
                &nnz2, sdat->y, sub2, map2);

      iZero(nnz2, map2, sub2);

      return;
    }

    if (nnz1 != 0)
    {
      printf("\n\n system error, nnz1!=0.\n");
      ShutDown();
      exit(0);
    }
  }

  pivlp_formyb(nrow, ncol, pdat->at, pdat->cc, sdat->basis, sdat->y,
               lu,
               2 * nrow, wrmemz);

  dZero(2 * nrow, wrmemz, NULL);
} /* pivlp_fomry */

static void pivlp_price(pivlppart *param,
                        pivlpt *pdat,
                        pivlpsolt *sdat,
                        int rsze,
                        int rset[],
                        int csze,
                        int cset[],
                        iteinft *iinf,
                        int *found,
                        int *isrow,
                        int *jinc,
                        double *sinc)

{
  int sk, *csk = sdat->csk;
  int fsbas, rsbas;
  int k, j, jsbas, nums, ncp;
  double mins, s, w12, tols = param->tols,
                       *cc = pdat->cc, *y = sdat->y;
  matrix *at = pdat->at;

  fsbas = false;
  rsbas = false;

  *found = false;
  *isrow = false;

  mins = -param->tols;

  if (sdat->solsta == SOL_FEAS)
    w12 = 1.0;
  else
    w12 = param->w12;

  if (rsze)
  {
    printf("\n\n internal error, rsze is defined.\n");
    ShutDown();
    exit(0);
  }

  iinf->ncp = 0;
  nums = 0;
  for (k = csze - iinf->ncp; k < csze; ++k)
  {
    j = cset[k];
    sk = csk[j];

    if (j == 1 && sk == BASIC)
      printf(" bas");

    if (sk != BASIC)
    {
      s = w12 * cc[j] - svDot(at->ia + j, y);

      if (j == 1)
        printf(" %e", s);

      if (sk != UPPER && s < tols)
        ++nums;

      if (sk != UPPER && s < mins)
      {
        *found = true;
        *isrow = false;
        *jinc = j;

        mins = s;
      }

      if (sk != LOWER && s > tols)
        ++nums;

      if (sk != LOWER && -s < mins)
      {
        *found = true;
        *isrow = false;
        *jinc = j;

        mins = -s;
      }

      if (sk == SPBAS)
      {
        fsbas = true;
        rsbas = false;
        jsbas = j;
      }
    }
  }

  if (!*found && param->nosbas && fsbas)
  {
    *found = true;
    *isrow = rsbas;
    *jinc = jsbas;
  }

  if (!*found)
  {
    /*
     * Major pricing step.
     */

    ncp = 0;
    nums = 0;
    for (k = csze - iinf->nac; k < csze; ++k)
    {
      j = cset[k];
      sk = csk[j];

      if (sk != BASIC)
      {
        s = w12 * cc[j] - svDot(at->ia + j, y);

        if (sk != UPPER && s < mins)
        {
          *found = true;
          *isrow = false;
          *jinc = j;

          mins = s;
        }

        if (sk != LOWER && -s < mins)
        {
          *found = true;
          *isrow = false;
          *jinc = j;

          mins = -s;
        }

        if (sk == SPBAS)
        {
          fsbas = true;
          rsbas = false;
          jsbas = j;
        }
      }
    }
    iinf->ncp = ncp;
  }

  if (*found)
    *sinc = w12 * cc[*jinc] - svDot(at->ia + (*jinc), y);
} /* pivlp_price */

static int prisim_chkopt(pivlppart *param,
                         pivlpt *pdat,
                         pivlpsolt *sdat,
                         pivlpslvit *slvi,
                         iteinft *iinf,
                         int wrsze,
                         double wrmem[])
{
  iinf->chkite = slvi->iter;

  pivlp_solinfo(param, pdat, sdat, slvi, wrsze, wrmem);

  return (slvi->mpbi <= param->tolx && slvi->mdbi <= param->tols);
} /* prisim_chkopt */

static void pivlp_chooserow(pivlppart *param,
                            pivlpt *pdat,
                            pivlpsolt *sdat,
                            double sj,
                            int nnzaj,
                            double valaj[],
                            int subaj[],
                            double xj,
                            int bkj,
                            double blj,
                            double buj,
                            chuzr_status *csta,
                            int *r,
                            double *alpha)
{
  int found = false;
  int bk;
  int i, k, t;
  pivbas *basis = sdat->basis;
  double l, u, aij, xbi, theta, mintheta, abs_aij, abs_arj, rtemp,
      tolapiv = param->tolapiv, tolx = param->tolx,
      *xb = sdat->xb;

  if (sj > 0.0)
  {
    for (k = 0; k < nnzaj; ++k)
      valaj[subaj[k]] = -valaj[subaj[k]];
  }

  /*
   * Phase 1 of Harris's ratiotest.
   */
  for (k = 0; k < nnzaj; ++k)
  {
    i = subaj[k];
    aij = valaj[i];

    if (fabs(aij) > tolapiv)
    {
      xbi = xb[i];
      switch (basis[i].vt)
      {
      case PIV_ART:
        bk = FX;
        l = 0.0;
        u = 0.0;
        break;
      case PIV_ROW:
        t = sdat->basis[i].j;
        bk = optdef_getrbki(pdat->rbk, t);
        l = optdef_getrbli(pdat->rbk, pdat->rbl, t),
        u = optdef_getrbui(pdat->rbk, pdat->rbu, t);
        break;
      case PIV_COL:
        t = sdat->basis[i].j;
        bk = optdef_getcbkj(pdat->cbk, t);
        l = optdef_getcblj(pdat->cbk, pdat->cbl, t);
        u = optdef_getcbuj(pdat->cbk, pdat->cbu, t);
        break;
      }

      if (aij > 0.0)
      {
        if (cmps_haslo(bk) && xbi > l - tolx)
        {
          theta = (xbi - l + tolx) / aij;
          if (!found || theta < mintheta)
          {
            found = true;
            mintheta = theta;
          }
        }
      }
      else
      {
        /*
         * aij <0.0
         */

        if (cmps_hasup(bk) && xbi < u + tolx)
        {
          theta = -(u + tolx - xbi) / aij;
          if (!found || theta < mintheta)
          {
            found = true;
            mintheta = theta;
          }
        }
      }
    }
  }

  *csta = chuzr_unb;
  if (sj <= 0.0)
  {
    if (cmps_hasup(bkj) &&
        (!found || (buj - xj) <= mintheta))
    {
      found = true;
      *alpha = max(buj - xj, 0.0);
      *csta = chuzr_upr;
    }
  }

  if (*csta == chuzr_unb && sj >= 0.0)
  {
    if (cmps_haslo(bkj) &&
        (!found || (xj - blj) <= mintheta))
    {
      found = true;
      *alpha = max(0.0, xj - blj);
      *csta = chuzr_low;
    }
  }

  if (!found && sdat->solsta == SOL_INFEAS)
  {
    found = false;
    abs_arj = 0.0;
    *alpha = 0.0;
    for (k = 0; k < nnzaj; ++k)
    {
      i = subaj[k];

      aij = valaj[i];
      abs_aij = fabs(aij);

      if (abs_aij > tolapiv)
      {
        xbi = xb[i];
        switch (basis[i].vt)
        {
        case PIV_ART:
          bk = FX;
          l = 0.0;
          u = 0.0;
          break;
        case PIV_ROW:
          t = sdat->basis[i].j;
          bk = optdef_getrbki(pdat->rbk, t);
          l = optdef_getrbli(pdat->rbk, pdat->rbl, t),
          u = optdef_getrbui(pdat->rbk, pdat->rbu, t);
          break;
        case PIV_COL:
          t = sdat->basis[i].j;
          bk = optdef_getcbkj(pdat->cbk, t);
          l = optdef_getcblj(pdat->cbk, pdat->cbl, t);
          u = optdef_getcbuj(pdat->cbk, pdat->cbu, t);
          break;
        }

        if (aij < 0.0)
        {
          if (cmps_haslo(bk) && xbi < l - tolx)
          {
            rtemp = xbi - aij * (*alpha) - l;
            if (rtemp < -tolx || (rtemp <= tolx && (!found || abs_aij > abs_arj)))
            {
              found = true;

              *r = i;
              *csta = chuzr_blockl;
              *alpha = (xbi - l) / aij;
              abs_arj = abs_aij;
            }
          }
        }
        else
        {
          /*
           * aij > 0.0
           */

          if (cmps_hasup(bk) && xbi > u + tolx)
          {
            rtemp = u - (xbi - aij * (*alpha));
            if (rtemp < -tolx || (rtemp <= tolx && (!found || abs_aij > abs_arj)))
            {
              found = true;

              *r = i;
              *csta = chuzr_blocku;
              *alpha = (xbi - u) / aij;
              abs_arj = abs_aij;
            }
          }
        }
      }
    }

    if (!found)
    {
      printf("\n\n no pivot is found.\n");
      ShutDown();
      exit(0);
    }
  }
  else
  {
    if (*csta == chuzr_unb)
    {
      found = false;
      abs_arj = 0.0;
      for (k = 0; k < nnzaj; ++k)
      {
        i = subaj[k];

        aij = valaj[i];
        abs_aij = fabs(aij);

        if (abs_aij > tolapiv)
        {
          xbi = xb[i];
          switch (basis[i].vt)
          {
          case PIV_ART:
            bk = FX;
            l = 0.0;
            u = 0.0;
            break;
          case PIV_ROW:
            t = sdat->basis[i].j;
            bk = optdef_getrbki(pdat->rbk, t);
            l = optdef_getrbli(pdat->rbk, pdat->rbl, t),
            u = optdef_getrbui(pdat->rbk, pdat->rbu, t);
            break;
          case PIV_COL:
            t = sdat->basis[i].j;
            bk = optdef_getcbkj(pdat->cbk, t);
            l = optdef_getcblj(pdat->cbk, pdat->cbl, t);
            u = optdef_getcbuj(pdat->cbk, pdat->cbu, t);
            break;
          }

          if (aij > 0.0)
          {
            if (cmps_haslo(bk) && xbi >= l - tolx)
            {
              theta = max(xbi - l, 0.0) / aij;
              if (theta <= mintheta && abs_aij > abs_arj)
              {
                found = true;

                *r = i;
                *csta = chuzr_blockl;
                *alpha = theta;
                abs_arj = abs_aij;
              }
            }
          }
          else
          {
            /*
             * aij < 0.0
             */

            if (cmps_hasup(bk) && xbi < u + tolx)
            {
              theta = min(xbi - u, 0.0) / aij;
              if (theta < mintheta && abs_aij > abs_arj)
              {
                found = true;

                *r = i;
                *csta = chuzr_blocku;
                *alpha = theta;
                abs_arj = abs_aij;
              }
            }
          }
        }
      }
    }
  }

  if (sj > 0.0)
  {
    for (k = 0; k < nnzaj; ++k)
      valaj[subaj[k]] = -valaj[subaj[k]];

    *alpha = -*alpha;
  }

} /* pivlp_chooserow */

static void pivlp_updsol(pivlppart *param,
                         pivlpt *pdat,
                         pivlpsolt *sdat,
                         pivlpslvit *slvi,
                         int rsze,
                         int rset[],
                         int invrset[],
                         int csze,
                         int cset[],
                         int invcset[],
                         int isrow,
                         int incj,
                         int nnza,
                         double vala[],
                         int suba[],
                         int mapa[],
                         int csta,
                         int r,
                         double alpha)
{
  int *rsk = sdat->rsk, *csk = sdat->csk;
  int i, jr, k, itemp1, itemp2;
  pivbas *basis = sdat->basis;
  int vtr;
  double *rx = sdat->rx, *cx = sdat->cx, *xb = sdat->xb;

  if (isrow)
    rx[incj] += alpha;
  else
    cx[incj] += alpha;

  for (k = 0; k < nnza; ++k)
  {
    i = suba[k];
    xb[i] -= alpha * vala[i];
    vala[i] = 0.0;
    mapa[i] = 0;
  }

  if (csta == chuzr_low || csta == chuzr_upr)
    ++slvi->mves;

  if (isrow)
  {
    printf("\n\n row is not allowed.\n");
    ShutDown();
    exit(0);
  }
  else
  {
    switch (csta)
    {
    case chuzr_unb:
      printf("\n\n system error, unbounded pivot.\n");
      ShutDown();
      exit(0);
    case chuzr_low:
    case chuzr_upr:
      if (optdef_getcbkj(pdat->cbk, incj) == FR)
        csk[incj] = NUFRE;
      else if (csta == chuzr_low)
        csk[incj] = LOWER;
      else
        csk[incj] = UPPER;
      break;
    case chuzr_blockl:
    case chuzr_blocku:
      vtr = basis[r].vt;
      jr = basis[r].j;

      if (vtr == PIV_ROW)
        rx[jr] = xb[r];
      else
        cx[jr] = xb[r];

      if (csta == chuzr_blockl)
      {
        if (vtr == PIV_ROW)
          rsk[jr] = LOWER;
        else if (vtr == PIV_COL)
          csk[jr] = LOWER;
      }
      else
      {
        if (vtr == PIV_ROW)
          rsk[jr] += UPPER;
        else if (vtr == PIV_COL)
          csk[jr] = UPPER;
      }

      basis[r].vt = PIV_COL;
      basis[r].j = incj;
      csk[incj] = BASIC;
      xb[r] = cx[incj];

      if (vtr == PIV_COL)
      {
        itemp1 = invcset[jr];
        itemp2 = invcset[incj];

        cset[itemp1] = incj;
        cset[itemp2] = jr;

        invcset[incj] = itemp1;
        invcset[jr] = itemp2;
      }
      else if (vtr == PIV_ROW)
      {
        printf("\n\n row type is not allowed.\n");
        ShutDown();
        exit(0);
      }
      break;
    default:
      printf("\n\n system error, unknown case.\n");
      ShutDown();
      exit(0);
    }
  }
} /* pivlp_updsol */

static void pivlp_solve(pivlppart *param,
                        pivlpt *pdat,
                        pivlpsolt *sdat,
                        pivlpslvit *slvi,
                        void *prinf,
                        int prlog,
                        int rsze,
                        int rset[],
                        int invrset[],
                        int csze,
                        int cset[],
                        int invcset[],
                        iteinft *iinf,
                        gsdec *lu,
                        int wisze,
                        int wimem[],
                        int wrsze,
                        double wrmem[])
{
  int found, isrow, formyb, refac, formxb = false;
  chuzr_status csta;
  int bki,
      *cbk = pdat->cbk;
  int i, j, k, incj, itemp, r, nnz1, nnz2, nnz3, numi, wiszel, wiszez, wrszez,
      nrow = pdat->nrow, ncol = pdat->ncol,
      *ajsub, *wimeml, *wimemz,
      *sub1, *sub2, *sub3, *subi, *map1, *map2, *map3;
  int lresp;
  double xi, bli, bui, rtemp, incs, alpha,
      *cx = sdat->cx, *cbl = pdat->cbl, *cbu = pdat->cbu,
      *y = sdat->y,
      *ajval, *wrmemz, *val1, *val2, *val3;
  matrix *at = pdat->at;
  array *aj;

  if (wisze < 11 * nrow || wrsze < 3 * nrow)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  val1 = wrmem;
  val2 = val1 + nrow;
  val3 = val2 + nrow;

  sub1 = wimem;
  sub2 = sub1 + nrow;
  sub3 = sub2 + nrow;
  subi = sub3 + nrow;

  map1 = subi + nrow;
  map2 = map1 + nrow;
  map3 = map2 + nrow;

  wiszel = 3 * nrow;
  wimeml = sub1;

  wiszez = 3 * nrow;
  wimemz = map1;

  wrszez = 3 * nrow;
  wrmemz = val1;

  formxb = (param->wlev == 0);
  formyb = (param->wlev == 0);

  dZero(wrszez, wrmemz, NULL);
  iZero(wiszez, wimemz, NULL);

  pivlp_formsubi(param, pdat, sdat, slvi, &numi, subi);

  do
  {
    ++slvi->iter;

    if (slvi->iter >= param->maxiter)
    {
      slvi->sresp = ItrLimit;
      break;
    }

    /*
     * Check cset.
     */
    do
    {
      if (formxb)
      {
        pivlp_formxb(nrow, ncol,
                     pdat->at,
                     pdat->rbk, pdat->rbl, pdat->rbu,
                     sdat->basis,
                     sdat->rx, sdat->cx, sdat->xb, lu,
                     2 * nrow, val1);

        dZero(2 * nrow, val1, NULL);

        formxb = false;

        pivlp_formsubi(param, pdat, sdat, slvi, &numi, subi);
      }

      if (formyb)
      {
        pivlp_formy(param, pdat, sdat, slvi,
                    numi, subi,
                    lu,
                    wiszel, wimeml,
                    wiszez, wimemz,
                    wrszez, wrmemz);

        formyb = false;
      }

      pivlp_price(param, pdat, sdat,
                  rsze, rset, csze, cset,
                  iinf,
                  &found, &isrow, &incj, &incs);

      if (!found)
      {
        if (sdat->solsta == SOL_FEAS)
        {
          if (iinf->chkite == slvi->iter ||
              prisim_chkopt(param, pdat, sdat, slvi, iinf, nrow, wrmemz))
          {
            sdat->solsta = SOL_OPTIMAL;
            return;
          }
          else
          {
            formyb = true;
            dZero(nrow, wrmem, NULL);
          }
        }
        else
        {
          if (param->w12 == 0.0)
            return;
          else
          {
            param->w12 /= 100.0;
            if (param->w12 <= 1.0e-10)
              param->w12 = 0.0;

            formyb = true;
          }
        }
      }

    } while (!found);

    if (param->prfrq && !(slvi->iter % param->prfrq))
    {
      pivlp_xb2rxcx(pdat->nrow, sdat->basis, sdat->rx, sdat->cx, sdat->xb);
      pivlp_solinfo(param, pdat, sdat, slvi, wrszez, wrmemz);
      dZero(wrszez, wrmemz, NULL);
      prf_prilog(prinf, slvi);
    }

    if (!found)
    {
      printf("\n\n not found.\n");
      ShutDown();
      exit(0);
    }

    /*
     * Update incoming column.
     */
    if (isrow)
    {
      i = incj;
      printf("\n\n isrow is defined.\n");
      ShutDown();
      exit(0);
    }
    else
    {
      j = incj;

      xi = cx[j];

      bki = optdef_getcbkj(cbk, j);
      bli = optdef_getcblj(cbk, cbl, j);
      bui = optdef_getcbuj(cbk, cbu, j);

      aj = at->ia + j;
      ajsub = aj->ja;
      ajval = aj->an;

      nnz1 = 0;
      for (k = 0; k < aj->nn0; ++k)
      {
        i = ajsub[k];
        val1[i] = ajval[k];
        map1[i] = 1;
        sub1[nnz1++] = i;
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

      lu_ftranu(lu, &nnz3, val3, sub3, map3, &nnz2, val2, sub2, map2);
    }

    /*
     * Choose leaving variable.
     */
    pivlp_chooserow(param, pdat, sdat,
                    incs, nnz2, val2, sub2,
                    xi, bki, bli, bui,
                    &csta, &r, &alpha);

    /*
     * Update the solution.
     */
    pivlp_updsol(param, pdat, sdat, slvi,
                 rsze, rset, invrset,
                 csze, cset, invcset,
                 isrow, incj,
                 nnz2, val2, sub2, map2,
                 csta, r, alpha);

    if (csta == chuzr_blockl || csta == chuzr_blocku)
    {
      /*
       * Update the LU factorization and yb.
       */
      lresp = lu_supdl(lu, r, r, nnz1, val1, sub1, map1, sub2, sub3);

      slvi->sresp = pivlp_prolresp(param->prlev,
                                   nrow, nrow, ncol,
                                   NULL, NULL,
                                   at,
                                   sdat->basis, sdat->rsk, sdat->csk,
                                   sdat->rx, sdat->cx, sdat->xb, sdat->y,
                                   false,
                                   lu, lresp, &refac,
                                   min(2 + 11 * nrow, wisze), wimem,
                                   min(2 * nrow, wrsze), wrmem);

      if (refac)
      {
        formxb = true;
        iZero(wiszez, wimemz, NULL);
        dZero(min(wrszez, 2 * nrow), wrmemz, NULL);
      }

      if (slvi->sresp == ProcOk)
      {
        if (sdat->solsta == SOL_INFEAS)
          formyb = true;
        else if (sdat->solsta == SOL_FEAS)
        {
          /*
           * Update yb.
           */

          formyb = false;

          nnz1 = 1;
          val1[r] = 1.0;
          sub1[0] = r;
          map1[r] = 1;

          lu_btran2(lu, nnz1, val1, sub1, map1, &nnz2, val2, sub2, map2);

          for (k = 0; k < nnz2; ++k)
          {
            y[sub2[k]] += incs * val2[sub2[k]];
            val2[sub2[k]] = 0.0;
            map2[sub2[k]] = 0;
          }
        }
      }
    }
    else
    {
      dZero(nnz1, val1, sub1);
      iZero(nnz1, map1, sub1);
    }

  } while (slvi->sresp == ProcOk);
} /* pivlp_solve */

void pivlp_primalsimplex(pivlppart *param,
                         pivlpt *pdat,
                         pivlpsolt *sdat,
                         pivlpslvit *slvi,
                         void *priinf,
                         int prlog,
                         int rsze,
                         int rset[],
                         int invrset[],
                         int csze,
                         int cset[],
                         int invcset[],
                         gsdec *lu,
                         int wisze,
                         int wimem[],
                         int wrsze,
                         double wrmem[])
{
  int j, k, nrow = pdat->nrow, ncol = pdat->ncol;
  iteinft iinf = {0};
  pivbas *basis = sdat->basis;
  double tolx, tols;

  pivlp_setup(param, pdat, sdat, slvi,
              rsze, rset, invrset,
              csze, cset, invcset,
              &iinf);

  tolx = param->tolx;
  tols = param->tols;

  param->tolx = tolx * 10.0;
  param->tols = tols * 10.0;

  pivlp_solve(param, pdat, sdat, slvi, priinf, prlog,
              rsze, rset, invrset,
              csze, cset, invcset,
              &iinf,
              lu,
              wisze, wimem,
              wrsze, wrmem);

  for (k = 0; k < nrow; ++k)
  {
    j = basis[k].j;
    switch (basis[k].vt)
    {

    case PIV_ROW:
      sdat->rx[basis[k].j] = sdat->xb[k];
      break;
    case PIV_COL:
      if (sdat->csk[j] == LOWER &&
          fabs(sdat->cx[j] - optdef_getcblj(pdat->cbk, pdat->cbl, j)) > 1.0e-6)
        exit(0);

      if (sdat->csk[j] == UPPER &&
          fabs(sdat->cx[j] - optdef_getcbuj(pdat->cbk, pdat->cbu, j)) > 1.0e-6)
        exit(0);

      sdat->cx[basis[k].j] = sdat->xb[k];
      break;
    }
  }

  for (j = 0; j < ncol; ++j)
  {
    switch (sdat->csk[j])
    {
    case LOWER:
    case FIXED:
      sdat->cx[j] = optdef_getcblj(pdat->cbk, pdat->cbl, j);
      break;
    case UPPER:
      sdat->cx[j] = optdef_getcbuj(pdat->cbk, pdat->cbu, j);
      break;
    case NUFRE:
      sdat->cx[j] = 0.0;
      break;
    }
  }

  sdat->solsta = SOL_INFEAS;

  param->tolx = tolx;
  param->tols = tols;

  pivlp_formxb(nrow, ncol,
               pdat->at,
               pdat->rbk, pdat->rbl, pdat->rbu,
               sdat->basis,
               sdat->rx, sdat->cx, sdat->xb, lu,
               wrsze, wrmem);

  pivlp_solve(param, pdat, sdat, slvi, priinf, prlog,
              rsze, rset, invrset,
              csze, cset, invcset,
              &iinf,
              lu,
              wisze, wimem,
              wrsze, wrmem);

} /* pivlp_primalsimplex */
