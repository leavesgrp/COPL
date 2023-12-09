#include "LPcross.h"

//#define TIME_COUNT

static void XcrMsg(prfinfot *prfinf, gsdec *lu)
{
  gsmsg *linfo;
  char ss[80];
#ifdef TEST
  int bopt;
#endif

  fprintf(fout, " Cross-over Partition message\n");

  printf("\n   Partition message\n");
  sprintf(ss, " " IFMT, prfinf->psze);
  LeftDots(ss, 46);
  printf("     number of SB set " SFMT "\n", ss);
  LeftDots(ss, 48);
  fprintf(fout, "   number of SB set " SFMT "\n", ss);

  sprintf(ss, " " IFMT, prfinf->pbsze);
  LeftDots(ss, 41);
  printf("     number of basis in SB " SFMT "\n", ss);
  LeftDots(ss, 43);
  fprintf(fout, "   number of basis in SB " SFMT "\n", ss);

  sprintf(ss, " " IFMT, prfinf->zsze);
  LeftDots(ss, 46);
  printf("     number of NB set " SFMT "\n", ss);
  LeftDots(ss, 48);
  fprintf(fout, "   number of NB set " SFMT "\n", ss);

  sprintf(ss, " " IFMT, prfinf->zbsze);
  LeftDots(ss, 41);
  printf("     number of basis in NB " SFMT "\n", ss);
  LeftDots(ss, 43);
  fprintf(fout, "   number of basis in NB " SFMT "\n", ss);

  sprintf(ss, " " IFMT, prfinf->asze);
  LeftDots(ss, 37);
  printf("     number of basis in AT set " SFMT "\n", ss);
  LeftDots(ss, 39);
  fprintf(fout, "   number of basis in AT set " SFMT "\n", ss);

  printf("   Iteration message\n");
  fprintf(fout, " Cross-over Iteration Message\n");

  sprintf(ss, " " IFMT, prfinf->itep);
  LeftDots(ss, 35);
  printf("     number of primal iterations " SFMT "\n", ss);
  LeftDots(ss, 37);
  fprintf(fout, "   number of primal iterations " SFMT "\n", ss);

  sprintf(ss, " " IFMT, prfinf->ited);
  LeftDots(ss, 37);
  printf("     number of dual iterations " SFMT "\n", ss);
  LeftDots(ss, 39);
  fprintf(fout, "   number of dual iterations " SFMT "\n", ss);

  sprintf(ss, " %d", prfinf->itec);
  LeftDots(ss, 36);
  printf("     number of clean iterations " SFMT "\n", ss);
  LeftDots(ss, 38);
  fprintf(fout, "   number of clean iterations " SFMT "\n", ss);

  linfo = &lu->info;

  printf("   Final results\n");
  fprintf(fout, " Cross-over Optimization Results\n");

  sprintf(ss, " %3.1e", prfinf->pfeas1 + prfinf->mpbi);
  LeftDots(ss, 42);
  printf("     primal infeasibility " SFMT "\n", ss);
  LeftDots(ss, 44);
  fprintf(fout, "   primal infeasibility " SFMT "\n", ss);

  sprintf(ss, " %3.1e", prfinf->dfeas1 + prfinf->mdbi);
  LeftDots(ss, 44);
  printf("     dual infeasibility " SFMT "\n", ss);
  LeftDots(ss, 46);
  fprintf(fout, "   dual infeasibility " SFMT "\n", ss);

  sprintf(ss, " %.9e", prfinf->pobj);
  LeftDots(ss, 40);
  printf("     primal objective value " SFMT "\n", ss);
  LeftDots(ss, 42);
  fprintf(fout, "   primal objective value " SFMT "\n", ss);

  sprintf(ss, " %.9e", prfinf->dobj);
  LeftDots(ss, 42);
  printf("     dual objective value " SFMT "\n", ss);
  LeftDots(ss, 44);
  fprintf(fout, "   dual objective value " SFMT "\n", ss);

#ifdef TEST
  bopt = ((fabs(prfinf->pobj - prfinf->dobj) <=
           (1.0 + fabs(prfinf->dobj)) * pr->tolgap) &&
          (prfinf->pfeasi <= pr->tolp) &&
          (prfinf->dfeasi <= pr->told) &&
          (prfinf->mpbi <= pr->tolp) &&
          (prfinf->mdbi <= pr->told));

  if (bopt)
    fprintf(fres, "COPT&");
  else
    fprintf(fres, "CUNK&");
#endif
} /* XcrMsg */

void XcrProc(prfparamt *param,
             int nrow,
             int ncol,
             double cc[],
             double cf,
             matrix *a,
             matrix *at,
             int rbk[], /* NULL */
             double rbl[],
             double rbu[],
             int cbk[],
             double cbl[], /* NULL */
             double cbu[],
             int irsk[], /* NULL */
             int icsk[],
             double rx[], /* NULL */
             double cx[],
             double y[],
             double rs[], /* NULL */
             double cs[],
             pivbas *basis,
             int orsk[], /* NULL */
             int ocsk[],
             double orx[], /* NULL */
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
  int chksol, solok, sbas;
  int i, j, k, wiszel, wrszel, iter,
      *wimeml;
  double xj, rtemp, *wrmeml;
#ifdef TIME_COUNT
  clock_t Tm;
#endif

  printf(" BEGIN cross-over...\n");

  prfinf->resp = ProcOk;
  prfinf->solsta = SOL_UNK;
  prfinf->psze = 0;
  prfinf->pbsze = 0;
  prfinf->zsze = 0;
  prfinf->zbsze = 0;
  prfinf->asze = 0;

  prfinf->titep = 0;
  prfinf->tited = 0;
  prfinf->titec = 0;

  prfinf->tmvep = 0;
  prfinf->tmved = 0;
  prfinf->tmvec = 0;

  lu->tolrpiv = param->tollurpiv;

  for (j = 0; j < ncol; ++j)
    ocsk[j] = icsk[j];

  for (j = 0; j < ncol; ++j)
  {
    switch (icsk[j])
    {
    case SPBAS:
      xj = cx[j];
      if (cmps_haslo(optdef_getcbkj(cbk, j)))
        xj = max(xj, optdef_getcblj(cbk, cbl, j));

      if (cmps_hasup(optdef_getcbkj(cbk, j)))
        xj = min(xj, optdef_getcbuj(cbk, cbu, j));

      ocx[j] = xj;
      break;
    case LOWER:
    case FIXED:
      ocx[j] = optdef_getcblj(cbk, cbl, j);
      break;
    case UPPER:
      ocx[j] = optdef_getcbuj(cbk, cbu, j);
      break;
    case NUFRE:
      ocx[j] = 0.0;
    default:
      exit(0);
    }
  }
#ifdef TIME_COUNT
  Tm = GetTime();
#endif
  if (prfinf->resp == ProcOk)
  {
    /*
     * Primal phase.
     */
    wiszel = 2 + 2 * nrow + 2 * ncol + max(11 * nrow, 3 * nrow + 3 * ncol);

    if (wiszel > wisze)
      wimeml = iAlloc(wiszel, NULL);
    else
      wimeml = wimem;

    wrszel = 3 * nrow;
    if (wrszel > wrsze)
      wrmeml = dAlloc(wrszel, NULL);
    else
      wrmeml = wrmem;

    if ((!wiszel || wimeml) && (!wrszel || wrmeml))
      XcrPrimProc(param,
                  false,
                  nrow, ncol,
                  cc, cf, a, at,
                  rbk, rbl, rbu,
                  cbk, cbl, cbu,
                  irsk, icsk,
                  rx, cx,
                  basis,
                  orsk, ocsk,
                  orx, ocx, xb,
                  lu, prfinf,
                  wiszel, wimeml,
                  wrszel, wrmeml);
    else
      prfinf->resp = OutOfSpc;

    if (wrszel > wrsze)
      dFree(&wrmeml);

    if (wiszel > wisze)
      iFree(&wimeml);
  }
#ifdef TIME_COUNT
  printf("PRIMAL: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
  if (prfinf->resp == ProcOk)
  {
    /*
     * Dual phase.
     */
    wiszel = 2 + 3 * nrow + 1 * ncol + max(8 * nrow + max(6 * nrow, 1 * ncol), 2 * nrow + 3 * ncol);

    if (wiszel > wisze)
      wimeml = iAlloc(wiszel, NULL);
    else
      wimeml = wimem;

    wrszel = 2 * nrow + 2 * ncol;
    if (wrszel > wrsze)
      wrmeml = dAlloc(wrszel, NULL);
    else
      wrmeml = wrmem;

    if ((!wiszel || wimeml) && (!wrszel || wrmeml))
    {
      XcrDualProc(param,
                  false,
                  nrow, ncol,
                  cc, cf, a, at,
                  rbk, rbl, rbu,
                  cbk, cbl, cbu,
                  irsk, icsk,
                  y, rs, cs,
                  basis,
                  orsk, ocsk,
                  orx, ocx, xb,
                  yb,
                  lu,
                  prfinf,
                  wiszel, wimeml,
                  wrszel, wrmeml);
    }
    else
      prfinf->resp = OutOfSpc;

    if (wrszel > wrsze)
      free(wrmeml);

    if (wiszel > wisze)
      free(wimeml);
  }
#ifdef TIME_COUNT
  printf("DUAL: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
  sbas = false;

  chksol = true;
  for (iter = 0; prfinf->resp == ProcOk && iter < 2; ++iter)
  {
    solok = false;
    if (chksol)
    {
      wrszel = max(2 * nrow, ncol);
      if (wrszel > wrsze)
        wrmeml = dAlloc(wrszel, NULL);
      else
        wrmeml = wrmem;

      if (!wrszel || wrmeml)
      {
        pivlp_formxb(nrow, ncol, at,
                     rbk, rbl, rbu,
                     basis,
                     orx, ocx, xb,
                     lu,
                     wrszel, wrmeml);

        pivlp_formyb(nrow, ncol, at, cc, basis, yb, lu, wrszel, wrmeml);

        sbas = false;
        if (orsk)
        {
          for (k = 0; k < nrow && orsk[k] != SPBAS; ++k)
            ;
          sbas = (k < nrow);
        }

        if (!sbas)
        {
          for (k = 0; k < ncol && ocsk[k] != SPBAS; ++k)
            ;
          sbas = (k < ncol);
        }

        if (!sbas)
        {
          XcrPrimMsg(param,
                     nrow, ncol,
                     cc, cf, a, at,
                     rbk, rbl, rbu,
                     cbk, cbl, cbu,
                     basis, orsk, ocsk,
                     orx, ocx, xb, prfinf,
                     wrszel, wrmeml);

          XcrDualMsg(param,
                     nrow, ncol,
                     cc, cf, a, at,
                     rbk, rbl, rbu,
                     cbk, cbl, cbu,
                     basis, orsk, ocsk,
                     yb, prfinf,
                     wrszel, wrmeml);

          solok = ((prfinf->mpbi <= max(param->tolx,
                                        param->rtolx * max(1.0, prfinf->nrmbi))) &&
                   (prfinf->mdbi <= max(param->tols,
                                        param->rtols * max(1.0, prfinf->nrmcbi))));
        }
        else
          solok = false;
      }
      else
        prfinf->resp = OutOfSpc;

      if (wrszel > wrsze)
        free(wrmeml);
    }

    if (solok || iter == 1)
    {
      prfinf->solsta = SOL_FEAS;
      break;
    }
#ifdef TIME_COUNT
    printf("CHECK: %d ms\n", GetTime() - Tm);
    Tm = GetTime();
#endif

    /*
     * Clean phase.
     */
    wiszel = 2 * ncol + 14 * nrow;
    if (irsk)
      wiszel += 2 * nrow;

    if (wiszel > wisze)
      wimeml = iAlloc(wiszel, "wimenl, XcrProc.");
    else
      wimeml = wimem;

    wrszel = 3 * nrow;
    if (wrszel > wrsze)
      wrmeml = dAlloc(wrszel, NULL);
    else
      wrmeml = wrmem;

    if ((!wiszel || wimeml) &&
        (!wrszel || wrmeml))
    {
      rtemp = lu->tolrpiv;
      if (rtemp <= 0.9)
      {
        rtemp = min(0.9, 2 * rtemp);
        lu->tolrpiv = rtemp;
      }

      XcrCleaProc(param,
                  nrow, ncol,
                  cc, cf,
                  a, at,
                  rbk, rbl, rbu,
                  cbk, cbl, cbu,
                  irsk, icsk,
                  rx, cx,
                  y, rs, cs,
                  basis,
                  orsk, ocsk,
                  orx, ocx,
                  xb, yb,
                  lu,
                  prfinf,
                  wiszel, wimeml,
                  wrszel, wrmeml);

      if (wiszel > wisze)
        free(wimeml);

      if (wrszel > wrsze)
        free(wrmeml);
    }
  }
#ifdef TIME_COUNT
  printf("CLEAN: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif

  if (prfinf->resp == ProcOk)
  {
    /*
     * Setup partition information.
     */
    for (j = 0; j < ncol; ++j)
    {
      if (icsk[j] == SPBAS)
      {
        ++prfinf->psze;
        if (ocsk[j] == BASIC)
          ++prfinf->pbsze;
      }
      else if (icsk[j] == LOWER || icsk[j] == UPPER)
      {
        ++prfinf->zsze;
        if (ocsk[j] == BASIC)
          ++prfinf->zbsze;
      }
    }

    for (i = 0; i < nrow; ++i)
      if (basis[i].vt == PIV_ART)
        ++prfinf->asze;
  }

  if (param->prlev >= 0)
    XcrMsg(prfinf, lu);

  printf(" END cross-over\n\n");

} /* XcrProc */

void XcrMain(int nrow, int ncol,
             int nblo, int nbup,
             int nbfr, int nbfx,
             int bset[], double cf,
             matrix *a, matrix *at,
             double b[], double c[],
             double u[], int csk[],
             double x[], double y[],
             double cs[], pivbas *basis,
             int cskb[], double cxb[],
             double xb[], double yb[])
{
  int *cbk;
  int k, pri, dua, wisze, wrsze,
      *bsetu = bset + (nblo - nbup), *wimem;
  gsdec *lu;
  double *cbu, *wrmem;
  matrix *a0 = a;
  prfinfot prfinf = {0};
  prfparamt prfpar = {0};

  prfpar.tolrpiv = 1e-12;
  prfpar.tolapiv = 1e-08;
  prfpar.tollurpiv = 0.01;
  prfpar.tolx = 1e-08;
  prfpar.tols = 1e-08;
  prfpar.rtolx = 1e-12;
  prfpar.rtols = 1e-12;
  prfpar.prlev = 6;
  prfpar.cmaxite = 10000;
  prfpar.prfrq = 100;

  prfinf.resp = ProcOk;

  if (!a)
    MtxTrans(nrow, at, false, NULL, false, &a);

  pri = 2 + 2 * nrow + 2 * ncol + max(11 * nrow, 3 * nrow + 3 * ncol);
  dua = 2 + 3 * nrow + 1 * ncol + max(8 * nrow + max(6 * nrow, 1 * ncol), 2 * nrow + 3 * ncol);

  wisze = max(pri, dua);
  wrsze = 2 * nrow + 2 * max(nrow, ncol);

  cbu = dAlloc(ncol, "ot->cbu, CrossOver");
  if (nbup || nbfr || nbfx)
    cbk = iAlloc(ncol, "ot->cbk, CrossOver");
  else
    cbk = NULL;

  wimem = iAlloc(wisze, "imem, CrossOver");
  wrmem = dAlloc(wrsze, "rmem, CrossOver");
  lu = LufAlloc(nrow, nrow, 0, true, "lu, CrossOver");

  if (nbup || nbfx || nbfr)
  {
    for (k = 0; k < nblo - nbup; ++k)
      cbk[bset[k]] = LO;

    for (k = 0; k < nbup; ++k)
    {
      cbk[bsetu[k]] = BO;
      cbu[bsetu[k]] = u[k];
    }

    for (k = nblo; k < nblo + nbfr; ++k)
      cbk[bset[k]] = FR;

    for (k = nblo + nbfr; k < nblo + nbfr + nbfx; ++k)
      cbk[bset[k]] = FX;
  }

  XcrProc(&prfpar,
          nrow, ncol,
          c, cf,
          a, at,
          NULL, b, b,
          cbk, NULL, cbu,
          NULL, csk,
          NULL,
          x, y, NULL, cs,
          basis,
          NULL, cskb,
          NULL, cxb, xb, yb,
          lu, &prfinf,
          wisze, wimem,
          wrsze, wrmem);

  iFree(&wimem);
  dFree(&wrmem);
  LufFree(&lu);

  if (!a0)
    MtxFree(&a);
} /* XcrMain */

void CrossOverProc(optmod *opt)
{
  int *cskb, i, j;
  double *cxb, *xb, *yb;
  pivbas *basis;

  basis = PbsAlloc(opt->m, "opt->bass, HsdProc");
  cskb = iAlloc(opt->n, "opt->cskb, HsdProc");
  cxb = dAlloc(opt->n, "opt->cxb, HsdProc");
  xb = dAlloc(opt->m, "opt->xb, HsdProc");
  yb = dAlloc(opt->m, "opt->yb, HsdProc");

  opt->nfx = 0;

  XcrMain(opt->m, opt->n,
          opt->nlb, opt->nub,
          opt->nfr, opt->nfx,
          opt->bid, opt->obj,
          opt->a, opt->at,
          opt->b, opt->c, opt->u,
          opt->csk, opt->x,
          opt->y, opt->z,
          basis, cskb, cxb, xb, yb);

  for (j = 0; j < opt->n; j++)
  {
    opt->x[j] = 0.0;
    opt->z[j] = 0.0;
  }

  for (j = 0; j < opt->n; j++)
    opt->x[j] = cxb[j];

  for (i = 0; i < opt->m; i++)
    opt->y[i] = yb[i];

  mTimesv(true, opt->m, opt->n, 1.0, opt->at, opt->y, 0.0, opt->z);

  for (j = 0; j < opt->n; j++)
  {
    opt->z[j] = opt->c[j] - opt->z[j];
    if (opt->z[j] < opt->x[j])
      opt->z[j] = 0.0;
  }

  PbsFree(&basis);
  iFree(&cskb);
  dFree(&cxb);
  dFree(&xb);
  dFree(&yb);

} /* CrossOverProc */
