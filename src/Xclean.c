#include "LPcross.h"

double CompBndInf(double x,
                  int bk,
                  double bl,
                  double bu)
{
  double r = 0.0;

  if (cmps_haslo(bk))
    r = max(r, bl - x);

  if (cmps_hasup(bk))
    r = max(r, x - bu);

  return (r);
} /* CompBndInf */

int pivlp_solinfo(pivlppart *param,
                  pivlpt *pdat,
                  pivlpsolt *sdat,
                  pivlpslvit *slvi,
                  int wrsze,
                  double wrmem[])
{
  int i, j, k, ndbi;
  double viol, sj, dbi, mdbi, sdbi, cbi, dfeasi, pobj = pdat->cf;

  if (wrsze < pdat->nrow)
    return (OutOfSpc);

  /*
   * Primal check.
   */
  slvi->mpbi = 0.0;
  slvi->spbi = 0.0;
  slvi->npbi = 0;
  viol = 0;

  dZero(pdat->nrow, wrmem, NULL);
  for (k = 0; k < pdat->nrow; ++k)
  {
    wrmem[k] = 0.0;
    switch (sdat->basis[k].vt)
    {
    case PIV_ART:
      viol = fabs(sdat->xb[k]);
      wrmem[k] = -sdat->xb[k];
      break;
    case PIV_ROW:
      if (sdat->rx)
      {
        i = sdat->basis[k].j;
        viol = CompBndInf(sdat->rx[i],
                          optdef_getrbki(pdat->rbk, i),
                          optdef_getrbli(pdat->rbk, pdat->rbl, i),
                          optdef_getrbui(pdat->rbk, pdat->rbu, i));
      }
      break;
    case PIV_COL:
      j = sdat->basis[k].j;
      viol = CompBndInf(sdat->cx[j],
                        optdef_getcbkj(pdat->cbk, j),
                        optdef_getcblj(pdat->cbk, pdat->cbl, j),
                        optdef_getcbuj(pdat->cbk, pdat->cbu, j));

      break;
    }

    slvi->mpbi = max(slvi->mpbi, viol);
    slvi->spbi += viol;
    if (viol > param->tolx)
      ++slvi->npbi;
  }

  mTimesv(false, pdat->nrow, pdat->ncol, -1.0, pdat->at, sdat->cx, 1.0, wrmem);

  if (sdat->rx)
    addVect(pdat->nrow, 1.0, sdat->rx, NULL, wrmem);
  else
  {
    for (k = 0; k < pdat->nrow; ++k)
      wrmem[k] += optdef_getrbli(pdat->rbk, pdat->rbl, k);
  }

  slvi->pfeasi = dNorm0(pdat->nrow, wrmem, NULL);
  slvi->pfeas1 = dNorm1(pdat->nrow, wrmem);

  dfeasi = 0.0;
  cbi = 0.0;
  mdbi = 0.0;
  ndbi = 0;
  sdbi = 0.0;

  for (j = 0; j < pdat->ncol; ++j)
  {
    pobj += pdat->cc[j] * sdat->cx[j];
    sj = pdat->cc[j] - svDot(pdat->at->ia + j, sdat->y);
    switch (sdat->csk[j])
    {
    case BASIC:
      dbi = 0.0;
      dfeasi = max(dfeasi, fabs(sj));
      cbi = max(cbi, fabs(pdat->cc[j]));
      break;
    case SPBAS:
    case NUFRE:
      dbi = fabs(sj);
      break;
    case LOWER:
      dbi = max(0.0, -sj);
      break;
    case UPPER:
      dbi = max(0.0, sj);
      break;
    default:
      dbi = 0.0;
    }

    if (dbi >= param->tols)
      ++ndbi;

    mdbi = max(dbi, mdbi);
    sdbi += dbi;
  }

  slvi->pobj = pobj;

  slvi->dfeasi = dfeasi;
  slvi->mdbi = mdbi;
  slvi->sdbi = sdbi;
  slvi->cbi = cbi;
  slvi->ndbi = ndbi;

  return (ProcOk);
} /* pivlp_solinfo */
/*
static void PrintXcleanHead(void)
{
  printf("\n     %-8s  %-8s  %-12s  %-12s  %-16s",
           "ITER","NPBI","MPBI","SPBI","POBJ");

}  PrintXcleanHead */

int prf_prilog(void *info,
               pivlpslvit *slvi)
{
  prfparamt *param = (prfparamt *)info;

  if (param->prlev >= 0)
  {
    printf("\n     %-8d  %-8d  %-12.3e  %-12.3e  %-16.6e",
           slvi->iter, slvi->npbi, slvi->mpbi, slvi->spbi, slvi->pobj);
  }
  return (true);
} /* prf_prilog */

void XcrCleaProc(prfparamt *param,
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
                 double rx[],
                 double cx[],
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
  int incall;
  int j, rsze, csze, csze0, wiszel, wrszel, itemp, itr,
      *rset, *invrset, *cset, *invcset, *wimeml;
  pivlppart lppar;
  pivlpslvit lpsinf = {0};
  pivlpsolt lpsol;
  pivlpt lp;
  double theta, *wrmeml;

  if (param->prlev >= 0)
    printf("   Begin clean phase...\n");

  lp.nrow = nrow;
  lp.ncol = ncol;

  lp.cc = cc;
  lp.cf = cf;

  lp.a = a;
  lp.at = at;

  lp.rbk = rbk;
  lp.rbl = rbl;
  lp.rbu = rbu;

  lp.cbk = cbk;
  lp.cbl = cbl;
  lp.cbu = cbu;

  lpsol.basis = basis;
  lpsol.rsk = orsk;
  lpsol.csk = ocsk;

  lpsol.xb = xb;
  lpsol.rx = orx;
  lpsol.cx = ocx;

  lpsol.y = yb;
  lpsol.cs = cs;

  lppar.wlev = 0;
  lppar.nosbas = true;
  lppar.maxiter = param->cmaxite;
  lppar.prfrq = param->prfrq;
  lppar.tolx = param->tolx;
  lppar.tols = param->tols;
  lppar.tolapiv = 1.0e-8;

  wrszel = wrsze;
  wrmeml = wrmem;

  pivlp_solinfo(&lppar, &lp, &lpsol, &lpsinf, wrszel, wrmeml);

  /*
   * Setup cset.
   */
  itemp = 2 * ncol;
  if (irsk)
    itemp = 2 * nrow;

  if (wisze < itemp)
  {
    printf("\n\n system error, work space size.\n");
    ShutDown();
    exit(0);
  }

  wiszel = wisze - itemp;
  wimeml = wimem + itemp;

  /***
  if (param->prlev>=0)
    PrintXcleanHead();
  ***/

  rset = NULL;
  // if (irsk)  { // Tian Xie, 2016.04.23
  rset = wimem;
  invrset = wimem + nrow;
  //} // Tian Xie, 2016.04.23

  cset = wimem + (itemp - 2 * ncol);
  invcset = cset + ncol;

  if (irsk)
    exit(0);

  if (lppar.maxiter)
  {

    itr = 0;
    theta = 10.0;
    incall = false;
    csze = 0;

    do
    {
      csze0 = csze;

      ++itr;
      rsze = 0;
      csze = 0;
      for (j = 0; j < ncol; ++j)
      {
        invcset[j] = ncol;
        if ((incall && ocsk[j] != FIXED) ||
            ocsk[j] == BASIC ||
            icsk[j] == SPBAS ||
            (icsk[j] == LOWER && theta * fabs(cx[j] - optdef_getcblj(cbk, cbl, j)) >= fabs(cs[j])) ||
            (icsk[j] == UPPER && theta * fabs(optdef_getcbuj(cbk, cbu, j) - cx[j]) >= fabs(cs[j])))
        {
          cset[csze] = j;
          invcset[j] = csze++;
        }
      }

      if (csze != csze0)
      {

        pivlp_primalsimplex(&lppar,
                            &lp,
                            &lpsol,
                            &lpsinf,
                            (void *)param,
                            j,
                            rsze, rset, invrset,
                            csze, cset, invcset,
                            lu,
                            wiszel, wimeml,
                            wrszel, wrmeml);

        lppar.wlev = 1;
      }

      theta *= theta;

      if (incall)
        break;

      incall = (itr >= 3 || csze == ncol);
    } while (lpsinf.sresp == ProcOk);
  }

  pivlp_xb2rxcx(lp.nrow, lpsol.basis, lpsol.rx, lpsol.cx, lpsol.xb);
  pivlp_solinfo(&lppar, &lp, &lpsol, &lpsinf, wrszel, wrmeml);

  prfinf->itec = lpsinf.iter;
  prfinf->mvec = lpsinf.mves;

  prfinf->titec += lpsinf.iter;
  prfinf->tmvec += lpsinf.mves;

  if (param->prlev >= 0)
    printf("   End clean.\n");
} /* XcrCleaProc */
