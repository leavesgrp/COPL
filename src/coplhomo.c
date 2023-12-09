#include "LPdefine.h"

////#define TIME_COUNT

static void InitOptSol(optmod *opt)
{
  int i, j, k, m, n, nub, nfr, nlb, *sub;
  array *ak, *rows;

  m = opt->m;
  n = opt->n;
  nlb = opt->nlb;
  nub = opt->nub;
  nfr = opt->nfr;
  sub = opt->bid;
  rows = opt->at->ia;

  opt->x = dAlloc(n, "opt->x, InitOptSol");
  opt->xu = dAlloc(nub, "opt->xu, InitOptSol");
  opt->y = dAlloc(m, "opt->y, InitOptSol");
  opt->yu = dAlloc(nub, "opt->yu, InitOptSol");
  opt->z = dAlloc(n, "opt->z, InitOptSol");
  opt->zu = dAlloc(nub, "opt->zu, InitOptSol");
  opt->xf = dAlloc(nfr, "opt->xf, InitOptSol");
  opt->csk = iAlloc(n, "opt->csk, InitOptSol");

  opt->lperm = iAlloc(opt->nden, "opt->lperm, InitOptSol");
  opt->linvp = iAlloc(opt->nden, "opt->linvp, InitOptSol");

  for (k = 0; k < nlb; k++)
  {
    j = sub[k];
    opt->x[j] = 1.0;
    opt->z[j] = 1.0;
  }

  for (k = 0; k < nub; k++)
  {
    opt->xu[k] = 1.0;
    opt->yu[k] = 0.0;
    opt->zu[k] = 1.0;
  }

  for (k = 0; k < m; k++)
    opt->y[k] = 0.0;

  for (k = 0; k < nfr; k++)
  {
    ak = rows + sub[nlb + k];
    opt->xf[k] = 1.0 + dNorm0(ak->nn0, ak->an, NULL);
    j = sub[nlb + k];
    opt->x[j] = 0.0;
    opt->z[j] = 0.0;
  }

  for (i = 0; i < n; i++)
    opt->csk[i] = 0;

  opt->nfx = 0;
  opt->tau = 1.0;
  opt->kap = 1.0;
} /* InitOptSol */

static void ClClean(chol *cl)
{
  int i, n;

  n = cl->nrow;
  for (i = 0; i < n; i++)
  {
    cl->ud[i] = 0.0;
    cl->naat[i] = 0;
  }
} /* ClClean */

static double dDot(int n, double *x, double *y)
{
  int i;
  double r;

  r = 0.0;
  if (n <= 0)
    return r;

  for (i = 0; i < n; i++)
    r += x[i] * y[i];
  return r;
} /* dDot */

static double dDots(int n, double *x, double *y, int *s)
{
  int i;
  double r;

  r = 0.0;
  if (n <= 0)
    return r;

  if (!s)
    ErrorProc(SysError, "null index set.");

  for (i = 0; i < n; ++i)
    r += x[s[i]] * y[s[i]];
  return r;
} /* dDots */

static void SetVect(int n, double alfa, double *x, double *y)
{
  int i;

  if (alfa == 1.0)
  {
    for (i = 0; i < n; ++i)
      y[i] = x[i];
  }

  else if (alfa == -1.0)
  {
    for (i = 0; i < n; ++i)
      y[i] = -x[i];
  }

  else
  {
    for (i = 0; i < n; ++i)
      y[i] = alfa * x[i];
  }
} /* SetVect */

void ShowHead(void)
{
  printf("\n ITER" SPC8 "POBJ" SPCE "DOBJ" SPCA "GAP" SPC6 "PINF" SPC5 "DINF\n");
  fprintf(fout, " ITER" SPC8 "POBJ" SPCE
                "DOBJ" SPCA "GAP" SPC6 "PINF" SPC5 "DINF\n");
} /* ShowHead */

void ShowInfor(optmod *ot)
{
  printf(" " I4FT "  " EGFT "  " EGFT "  " E3FT "  " E3FT "  " E3FT "\n", ot->iter, ot->pobj,
         ot->dobj, ot->rgap, ot->rrpnm, ot->rrdnm);

  fprintf(fout, " " I4FT "  " EGFT "  " EGFT "  " E3FT "  " E3FT "  " E3FT "\n", ot->iter,
          ot->pobj, ot->dobj, ot->rgap, ot->rrpnm,
          ot->rrdnm);
} /* ShowInfor */

static void CompValues(optmod *opt,
                       double *tht,
                       double *thtu,
                       double *thtg,
                       double *rp,
                       double *rd,
                       double *rg,
                       double *mu,
                       double *rv,
                       int isze,
                       int *imem,
                       int rsze,
                       double *rmem)
{
  int j, k, m, n, nlb, nub, nfr, nfx, ntt, unu, *bid, *bidu, *bidf;
  double cgap, xzmn, pnmu, rtmp, tau, *c, *b, *u, *x, *xu, *y, *yu, *z, *zu;

  m = opt->m;
  n = opt->n;
  nlb = opt->nlb;
  nub = opt->nub;
  nfr = opt->nfr;
  nfx = opt->nfx;
  ntt = nlb + nfr;
  unu = nlb - nub;
  bid = opt->bid;
  bidu = bid + unu;
  bidf = bid + nlb;

  tau = opt->tau,
  c = opt->c,
  b = opt->b,
  u = opt->u,
  x = opt->x,
  xu = opt->xu,
  y = opt->y,
  yu = opt->yu,
  z = opt->z,
  zu = opt->zu;

  xzmn = MaxPositive;
  cgap = 0.0;

  if (opt->ista == ProcOk)
  {
    for (k = 0; k < nlb; ++k)
    {
      j = bid[k];
      tht[j] = z[j] / x[j];

      if (tht[j] <= 0.0)
      {
        opt->ista = ProcStall;
        /* ErrorProc(SysError,"negative iterate.");*/
        return;
      }
      rtmp = x[j] * z[j];
      cgap += rtmp;

      if (rtmp < xzmn)
        xzmn = rtmp;
    }

    for (k = 0; k < nub; ++k)
    {
      j = bid[unu + k];
      thtu[k] = zu[k] / xu[k];
      tht[j] += thtu[k];

      if (thtu[k] <= 0.0)
      {
        opt->ista = ProcStall;
        return;
        /* ErrorProc(SysError,"negative iterate.");*/
      }

      rtmp = xu[k] * zu[k];
      cgap += rtmp;

      if (rtmp < xzmn)
        xzmn = rtmp;
    }

    rtmp = opt->tau * opt->kap;
    cgap += rtmp; // Tian Xie, 2016.04.24

    if (rtmp < xzmn)
      xzmn = rtmp;

    *rv = opt->nall * xzmn / cgap;
    *mu = cgap / opt->nall;

    opt->rdf = 0.0;
    for (k = 0; k < nfr; ++k)
    {
      j = bidf[k];
      tht[j] = opt->rho * opt->xf[k];
      opt->rdf += opt->xf[k] * x[j] * x[j];
    }

    opt->rdf *= opt->rho;
    *thtg = opt->kap / opt->tau + opt->rdf / (tau * tau);

    opt->cx = dDots(nlb, c, x, bid) + dDots(nfr, c, x, bidf);
    opt->by = dDot(m, b, y) + dDot(nub, u, yu);

    SetVect(m, tau, b, rp);
    mTimesv(false, m, n, -1.0, opt->at, x, 1.0, rp);

    for (k = 0; k < nlb; ++k)
      rd[bid[k]] = z[bid[k]] - c[bid[k]] * tau;
    for (k = 0; k < nfr; ++k)
      rd[bidf[k]] = -c[bidf[k]] * tau;

    mTimesv(true, m, n, 1.0, opt->at, y, 1.0, rd);

    for (k = 0; k < nub; ++k)
      rd[bidu[k]] += yu[k];

    for (k = 0; k < nfx; ++k)
      rd[bid[ntt + k]] = 0.0;

    *rg = opt->kap + opt->cx - opt->by;
    opt->dfnm = dSums(nfr, rd, bidf);
    opt->rdnm = dSums(nlb, rd, bid) + opt->dfnm;

    opt->ddtu = 0.0;
    for (k = 0; k < nub; ++k)
    {
      rtmp = yu[k] + zu[k];
      opt->ddtu += rtmp * rtmp;
    }

    opt->rdnm += opt->ddtu;

    opt->rpdt = dSums(m, rp, NULL);
    pnmu = 0.0;

    for (k = 0; k < nub; ++k)
    {
      rtmp = u[k] * tau - xu[k] - x[bidu[k]];
      pnmu += rtmp * rtmp;
    }

    opt->pdtu = pnmu;
    opt->rpdt += pnmu;
    opt->pnmu = sqrt(pnmu);
    opt->rpnm = sqrt(opt->rpdt);

    opt->rdnm = sqrt(opt->rdnm);
    opt->dfdt = dSums(nfr, rd, bidf);
    opt->rddt = dSums(nlb, rd, bid) + opt->dfdt;

    for (k = 0; k < nub; ++k)
    {
      rtmp = zu[k] + yu[k];
      opt->rddt += rtmp * rtmp;
    }

    opt->rgdt = (*rg) * (*rg);
    opt->smdt = opt->rpdt + opt->rddt + opt->rgdt;

    opt->pobj = opt->cx / tau + opt->obj;
    opt->dobj = opt->by / tau + opt->obj;

    opt->cgap = cgap;

    opt->rgap = fabs(opt->pobj - opt->dobj) / (1.0 + fabs(opt->dobj));
    opt->rrpnm = opt->rpnm / (opt->tau + dNorm2(n, x));
    opt->rrdnm = opt->rdnm / (opt->tau + sqrt(dSums(nlb, z, bid) + dSums(nub, zu, NULL) + dSums(m, y, NULL) + dSums(nub, yu, NULL)));
    opt->infe = opt->tau / (opt->inf0 * opt->kap);
  }
} /* CompValues */

static void SolLES(optmod *opt,
                   double *dx,
                   double *du,
                   double *sx,
                   double *su,
                   double *tx,
                   double *tu,
                   int rsze,
                   double *rmem)
{
  int k, m, n, nlb, nub, nfr,
      *bid, *bidu, *bidf;
  double rtmp, *rhs, *sol;

  if (rsze < 2 * (opt->m + opt->n))
    ErrorProc(SysError, "system error.");

  m = opt->m;
  n = opt->n;
  nlb = opt->nlb;
  nub = opt->nub;
  nfr = opt->nfr;
  bid = opt->bid;
  bidu = bid + (nlb - nub);
  bidf = bid + nlb;
  rhs = rmem;
  sol = rhs + (m + n);

  for (k = 0; k < nlb; ++k)
    rhs[bid[k]] = sx[bid[k]];
  for (k = 0; k < nub; ++k)
    rhs[bidu[k]] += du[k] * tu[k] - su[k];
  for (k = 0; k < nfr; ++k)
    rhs[bidf[k]] = sx[bidf[k]];

  dCopy(m, tx, rhs + n);

  SolArgLES(opt->cl, dx, opt->at, rhs, sol);

  dCopy(n, sol, sx);
  SetVect(m, -1.0, sol + n, tx);

  for (k = 0; k < nub; ++k)
  {
    rtmp = su[k];
    su[k] = tu[k] - sx[bidu[k]];
    tu[k] = du[k] * (tu[k] - sx[bidu[k]]) - rtmp;
  }
} /* SolLES */

static double MixDots(optmod *ot,
                      double *wc,
                      double *wb,
                      double *wu)
{
  int k, *bid, *bidf;
  double rtmp, tarh;

  bid = ot->bid;
  bidf = bid + ot->nlb;

  rtmp = 0.0;
  for (k = 0; k < ot->nlb; ++k)
    rtmp -= ot->c[bid[k]] * wc[bid[k]];
  tarh = ot->rho / ot->tau;
  for (k = 0; k < ot->nfr; ++k)
    rtmp -= (ot->c[bidf[k]] + tarh * ot->xf[k] * ot->x[bidf[k]]) * wc[bidf[k]];
  rtmp += dDot(ot->m, ot->b, wb) + dDot(ot->nub, ot->u, wu);
  return (rtmp);
} /* MixDots */

static void GetConst(optmod *opt,
                     double *d,
                     double *du,
                     double *pi,
                     double *p,
                     double *pu,
                     double *q,
                     double *qu,
                     int rsze,
                     double *rmem)
{
  int m, n, nub, k, *bidf;
  double rtmp;

  m = opt->m;
  n = opt->n;
  nub = opt->nub;
  bidf = opt->bid + opt->nlb;

  dCopy(n, opt->c, p);

  rtmp = opt->rho / opt->tau;
  for (k = 0; k < opt->nfr; ++k)
    p[bidf[k]] -= rtmp * opt->xf[k] * opt->x[bidf[k]];

  SetVect(m, -1.0, opt->b, q);
  for (k = 0; k < nub; k++)
  {
    pu[k] = 0.0;
    qu[k] = -opt->u[k];
  }

  SolLES(opt, d, du, p, pu, q, qu, rsze, rmem);

  *pi = MixDots(opt, p, q, qu);
} /* GetConst */

static int GenSearchDir(optmod *opt, double eta,
                        double *dl, double *du,
                        double dg, double *rp,
                        double *rd, double rg,
                        double cnst, double *s,
                        double *v, double *t,
                        double *w, double *dx,
                        double *dxu, double *dtau,
                        double *dy, double *dyu,
                        double *dz, double *dzu,
                        double *dkap, double *rpns,
                        double *rdns, double *rgns,
                        int rsze, double *rmem)
{
  int j, k, m, n, itmp, rszel, nlb, nub, nfr, *bl, *bu, *bf;
  double rtmp, tau, rr, rdx, rdz, smdt, smdt0, *rmeml, *fl,
      *fu, *gl, *gu, *pl, *pu, *ql, *qu, *c, *u, *x, *xu,
      *z, *zu, *yu;
  array *aj;

  m = opt->m;
  n = opt->n;
  nlb = opt->nlb;
  nub = opt->nub;
  nfr = opt->nfr;

  bl = opt->bid;
  bu = bl + (nlb - nub);
  bf = bl + nlb;

  itmp = (n + 2 * nub + m);
  if (rsze < 2 * itmp)
  {
    ErrorProc(NOT_MEMSPC, NULL);
    ShutDown();
    exit(0);
  }

  c = opt->c;
  u = opt->u;
  x = opt->x;
  xu = opt->xu;
  z = opt->z;
  zu = opt->zu;
  yu = opt->yu;

  rszel = rsze - itmp;
  rmeml = rmem + itmp;

  fl = rmem;
  fu = fl + n;
  gl = fu + nub;
  gu = gl + m;

  tau = opt->tau;

  for (k = 0; k < nlb; ++k)
  {
    j = bl[k];
    if (x[j] >= z[j])
      dz[j] /= x[j];
    else
      dz[j] /= z[j];
  }

  for (k = 0; k < nub; ++k)
  {
    if (xu[k] >= zu[k])
      dzu[k] /= xu[k];
    else
      dzu[k] /= zu[k];
  }

  SetVect(n, eta, rd, fl);
  SetVect(m, eta, rp, gl);

  for (k = 0; k < nub; ++k)
  {
    fu[k] = eta * (zu[k] + yu[k]);
    gu[k] = eta * (u[k] * tau - x[bu[k]] - xu[k]);
  }

  rr = eta * rg;

  for (k = 0; k < nlb; ++k)
  {
    j = bl[k];
    if (x[j] >= z[j])
      fl[j] += dz[j];
  }

  for (k = 0; k < nub; ++k)
    if (xu[k] >= zu[k])
      fu[k] += dzu[k];

  rr += (*dkap) / tau;

  for (k = 0; k < nlb; ++k)
  {
    j = bl[k];
    if (x[j] < z[j])
    {
      rtmp = dz[j];
      rr += c[j] * rtmp;
      aj = opt->at->ia + j;
      setArray(-rtmp, aj, gl);
    }
  }

  for (k = 0; k < nub; ++k)
  {
    j = bu[k];
    if (x[j] < z[j])
      gu[k] -= dz[j];

    if (xu[k] < zu[k])
      gu[k] -= dzu[k];
  }

  pl = dx;
  pu = dxu;
  ql = dy;
  qu = dyu;

  dCopy(n, fl, pl);
  dCopy(m, gl, ql);

  for (k = 0; k < nub; k++)
  {
    pu[k] = fu[k];
    qu[k] = gu[k];
  }

  SolLES(opt, dl, du, pl, pu, ql, qu, rszel, rmeml);

  if (dg - cnst <= 1.0e-13)
    *dtau = 0.0;
  else
    *dtau = (rr - MixDots(opt, pl, ql, qu)) / max(opt->kap / opt->tau, (dg - cnst));

  rtmp = *dtau;

  for (k = 0; k < nlb; ++k)
  {
    j = bl[k];

    rdx = pl[j] - s[j] * rtmp;
    if (x[j] >= z[j])
    {
      dz[j] = dz[j] - (z[j] * rdx) / x[j];
      dx[j] = rdx;
    }
    else
    {
      dx[j] = dz[j] + rdx;
      dz[j] = -(z[j] * rdx) / x[j];
    }
  }

  for (k = 0; k < nfr; ++k)
  {
    j = bf[k];
    dx[j] = pl[j] - s[j] * rtmp;
  }

  for (k = 0; k < nub; ++k)
  {
    rdz = pu[k] - v[k] * rtmp;
    if (xu[k] >= zu[k])
    {
      dzu[k] = dzu[k] - (zu[k] * rdz) / xu[k];
      dxu[k] = rdz;
    }
    else
    {
      dxu[k] = dzu[k] + rdz;
      dzu[k] = -(zu[k] * rdz) / xu[k];
    }
  }

  for (k = 0; k < m; ++k)
    dy[k] = ql[k] - t[k] * rtmp;

  for (k = 0; k < nub; ++k)
    dyu[k] = qu[k] - w[k] * rtmp;

  *dkap = ((*dkap) - opt->kap * rtmp) / tau;

  smdt0 = opt->smdt;

  SetVect(n, eta, rd, fl);
  SetVect(m, eta, rp, gl);

  for (k = 0; k < nub; ++k)
  {
    fu[k] = eta * (zu[k] + yu[k]);
    gu[k] = eta * (u[k] * tau - x[bu[k]] - xu[k]);
  }

  for (k = 0; k < nlb; ++k)
    fl[bl[k]] += -opt->c[bl[k]] * (*dtau) + dz[bl[k]];

  rtmp = opt->rho / opt->tau;
  for (k = 0; k < nfr; ++k)
    fl[bf[k]] -= (opt->c[bf[k]] - rtmp * opt->xf[k] * x[bf[k]]) * (*dtau) + opt->rho * opt->xf[k] * dx[bf[k]];

  mTimesv(true, m, n, 1.0, opt->at, dy, 1.0, fl);
  addVect(m, (*dtau), opt->b, NULL, gl);
  mTimesv(false, m, n, -1.0, opt->at, dx, 1.0, gl);

  (*rdns) = 0.0;
  (*rpns) = 0.0;

  for (k = 0; k < nub; ++k)
  {
    fl[bu[k]] += dyu[k];
    fu[k] = eta * (zu[k] + yu[k]) + dyu[k] + dzu[k];
    gu[k] += u[k] * (*dtau) - dxu[k] - dx[bu[k]];
    (*rdns) += fu[k] * fu[k];
    (*rpns) += gu[k] * gu[k];
  }

  (*rdns) += dSums(n, fl, NULL);
  (*rpns) += dSums(m, gl, NULL);

  rtmp = MixDots(opt, dx, dy, dyu);
  rr = (*dkap) + eta * (rg)-rtmp - opt->rdf * (*dtau) / (tau * tau);
  *rgns = rr * rr;
  smdt = *rdns + *rpns + *rgns;

  if (sqrt(*rpns + *rdns) <= 1.0e-10 ||
      *rpns + *rdns + *rgns <= opt->smdt / 10.0)
    return (true);

  return (false);
} /* GenSearchDir */

static void setPartmsg(optmod *opt,
                       double *dx,
                       double *dxu,
                       double dtau,
                       double *dz,
                       double *dzu,
                       double dkap,
                       int *cks,
                       int *nchg,
                       int *nlow)
{
  int ocks, k, j, sbas, low, lsta, nlb,
      nub, nfr, *bid, *bidu, *bidf;
  double rtmp, delta, *x, *xu, *z, *zu;

  sbas = 0;
  low = 0;
  lsta = 0;

  nlb = opt->nlb;
  nub = opt->nub;
  nfr = opt->nfr;
  bid = opt->bid;
  bidu = bid + (nlb - nub);
  bidf = bid + nlb;

  delta = -0.475;
  x = opt->x;
  xu = opt->xu;
  z = opt->z;
  zu = opt->zu;

  if (cks)
  {
    *nchg = 0;

    rtmp = MaxPositive;

    if (dtau <= delta * opt->tau && dkap <= delta * opt->kap)
    {
      rtmp = min(rtmp, min(opt->tau, opt->kap));
      lsta++;
    }

    for (k = 0; k < nlb - nub; ++k)
    {
      j = bid[k];
      ocks = cks[j];

      if (dx[j] <= delta * x[j] && dz[j] <= delta * z[j])
      {
        rtmp = min(rtmp, min(x[j], z[j]));
        lsta++;
      }

      if (fabs(dx[j]) * z[j] <= fabs(dz[j]) * x[j])
      {
        cks[j] = SPBAS;
        sbas++;
      }
      else
      {
        cks[j] = LOWER;
        low++;
      }

      if (cks[j] != ocks)
        (*nchg)++;
    }

    for (k = 0; k < nub; ++k)
    {
      ocks = cks[k];
      j = bidu[k];

      if (dxu[k] <= delta * xu[k] && dzu[k] <= delta * zu[k])
      {
        rtmp = min(rtmp, min(zu[k], xu[k]));
        lsta++;
      }

      if (fabs(dx[j]) * z[j] <= fabs(dz[j]) * x[j] &&
          fabs(dxu[k]) * zu[k] <= fabs(dzu[k]) * xu[k])
      {
        sbas++;
        cks[j] = SPBAS;
      }
      else if (fabs(dx[j]) * z[j] >= fabs(dz[j]) * x[j])
      {
        low++;
        cks[j] = LOWER;
      }
      else
      {
        low++;
        cks[j] = UPPER;
      }

      if (cks[j] != ocks)
        (*nchg)++;
    }

    sbas += nfr;
    for (k = 0; k < nfr; ++k)
      cks[bidf[k]] = SPBAS;

    sbas += nfr;
    for (k = 0; k < nfr; ++k)
      cks[bidf[k]] = SPBAS;
    *nlow = low;
  }
} /* setPartmsg */

static double RatioStep(optmod *opt,
                        double x[],
                        double xu[],
                        double dx[],
                        double dxu[])
{
  int i, j, k, n;
  double step;

  n = opt->nlb;
  step = 1.0;

  for (i = 0; i < n && dx[opt->bid[i]] > -1.0e-13; ++i)
    ;

  if (i < n)
  {
    for (j = i; j < n; j++)
    {
      k = opt->bid[j];
      if (x[k] + step * dx[k] < 0.0)
        step = -x[k] / dx[k];
    }
  }

  n = opt->nub;

  for (i = 0; i < n && dxu[i] > -1.0e-13; ++i)
    ;

  if (i < n)
  {
    for (j = i; j < n; j++)
    {
      if (xu[j] + step * dxu[j] < 0.0)
        step = -xu[j] / dxu[j];
    }
  }

  return (step);
} /* RatioStep */

static int ChkStep(optmod *opt,
                   double rg,
                   double eta,
                   double dx[],
                   double dxu[],
                   double dtau,
                   double dz[],
                   double dzu[],
                   double dkap,
                   double alfp,
                   double alft,
                   double alfd,
                   double alfk,
                   double beta,
                   double *lpval)
{
  int k, j, nlb, nub, *bid;
  double rtmp, mxz, xz, *x, *xu, *z, *zu;

  nlb = opt->nlb;
  nub = opt->nub;
  bid = opt->bid;

  x = opt->x;
  xu = opt->xu;
  z = opt->z;
  zu = opt->zu;

  rtmp = 0.0;
  mxz = MaxPositive;
  for (k = 0; k < nlb; ++k)
  {
    j = bid[k];
    xz = (x[j] + alfp * dx[j]) * (z[j] + alfd * dz[j]);
    rtmp += xz;
    mxz = min(mxz, xz);
  }

  for (k = 0; k < nub; ++k)
  {
    xz = (xu[k] + alfp * dxu[k]) * (zu[k] + alfd * dzu[k]);
    rtmp += xz;
    mxz = min(mxz, xz);
  }

  xz = (opt->tau + alft * dtau) * (opt->kap + alfk * dkap);
  rtmp += xz;
  mxz = min(mxz, xz);
  rtmp /= opt->nall;

  *lpval = mxz / rtmp;

  return (*lpval >= beta);
} /* ChkStep */

static void getInitPoint(optmod *opt,
                         optpar *par,
                         int isze,
                         int *imem,
                         int rsze,
                         double *rmem)
{
  matrix *at;
  int h, j, k, m, n, nlb, nub, nfr, afine, end, idch, boka,
      pid, *bid, *bidu, *bidf, *iwk;
  double mu, pi, rg, eta, gamma, lpv, phi, alfa, alfd, alfp, alfk, alft,
      dkap, dtau, thtg, rpnorm, rdnorm, rgnorm, rfval, dstp, pstp, dfac,
      pfac, taup, taud, btaup, rtmp, rstep, *s, *v, *t, *w, *x, *y, *xu, *yu,
      *z, *zu, *tht, *thtu, *thtd, *dx, *dy, *dz, *dxu, *dyu, *dzu,
      *rp, *rd, *rwk;

  m = opt->m;
  n = opt->n;
  nlb = opt->nlb;
  nub = opt->nub;
  nfr = opt->nfr;
  bid = opt->bid;
  bidu = bid + (nlb - nub);
  bidf = bid + nlb;

  dtau = 0.0;
  x = opt->x;
  xu = opt->xu;
  y = opt->y;
  yu = opt->yu;
  z = opt->z;
  zu = opt->zu;
  at = opt->at;

  rp = rmem;
  rd = rp + m;
  tht = rd + n;
  thtu = tht + n;
  thtd = thtu + nub;

  s = thtd + n;
  v = s + n;
  t = v + nub;
  w = t + m;

  dx = w + nub;
  dxu = dx + n;
  dy = dxu + nub;
  dyu = dy + m;
  dz = dyu + nub;
  dzu = dz + n;

  iwk = imem;
  rwk = dzu + nub;

  opt->nall = 1.0 + (double)(nlb + nub);
  opt->mu0 = (dDots(nlb, x, z, bid) + dDot(nub, xu, zu) + opt->tau * opt->kap) / opt->nall;
  opt->tau0 = opt->tau;
  opt->kap0 = opt->kap;
  mu = opt->mu0;
  opt->inf0 = opt->tau0 / opt->kap0;
  opt->rho = 1.0;
  opt->iter = 0;

  for (k = 0; k < nfr; k++)
    opt->xf[k] = 1.0 / (1.0 + 10.0 * fabs(opt->x[bidf[k]]));

  CompValues(opt, tht, thtu, &thtg, rp, rd, &rg,
             &mu, &rfval, isze, iwk, rsze, rwk);

  if ((opt->ista != ProcOk) ||
      (mu <= par->tolmu * opt->mu0 && opt->infe < par->tolinf) ||
      (opt->rgap <= par->tolgap && opt->rrpnm <= par->tolp && opt->rrdnm <= par->told))
    return;

  ClClean(opt->cl);
  idch = decChol(opt->cl, tht, at, isze, iwk, rsze, rwk);

  if (idch != ChlOk)
  {
    ShutDown();
    exit(0);
  }

  GetConst(opt, tht, thtu, &pi, s, v, t, w, rsze, rwk);
  phi = 1.0;
  end = 2;

  for (h = 0; h < end; h++)
  {
    eta = 1.0;

    for (k = 0; k < nlb; ++k)
      dz[bid[k]] = -x[bid[k]] * z[bid[k]] + mu;
    for (k = 0; k < nfr; ++k)
      dz[bidf[k]] = 0.0;
    for (k = 0; k < nub; ++k)
      dzu[k] = -xu[k] * zu[k] + mu;

    dkap = -opt->tau * opt->kap + mu;

    afine = GenSearchDir(opt, eta, tht, thtu, thtg, rp, rd, rg, pi, s, v, t, w,
                         dx, dxu, &dtau, dy, dyu, dz, dzu, &dkap, &rpnorm,
                         &rdnorm, &rgnorm, rsze, rwk);

    afine = afine && (rpnorm + rdnorm + rgnorm <= opt->smdt / 10.0 || sqrt(rpnorm + rdnorm + rgnorm) <= 1.0e-8);

    if (afine)
    {
      pstp = RatioStep(opt, x, xu, dx, dxu);
      dstp = RatioStep(opt, z, zu, dz, dzu);

      if (dtau < -par->tolx)
      {
        pstp = min(pstp, -opt->tau / dtau);
        dstp = min(dstp, -opt->tau / dtau);
      }

      if (dkap < -par->tolx)
      {
        pstp = min(pstp, -opt->kap / dkap);
        dstp = min(dstp, -opt->kap / dkap);
      }
    }

    else
    {
      fprintf(fout, " inaccurate affine search direction.\n");
      phi = 0.0;
    }

    if (!h)
    {
      gamma = 10.0;
      eta = 1.0;

      rtmp = (dDots(nlb, dx, dz, bid) + dDot(nub, dxu, dzu) + dkap * dtau) / opt->nall;
      rtmp = max(rtmp, 0.0);

      alfp = min(pstp, dstp);
      alfd = max(0.0, 1.0 - alfp);
      alfp = alfd * gamma * mu;

      for (k = 0; k < nlb; ++k)
        dz[bid[k]] = -x[bid[k]] * z[bid[k]] + alfp - phi * (dx[bid[k]] * dz[bid[k]] - rtmp);
      for (k = 0; k < nfr; ++k)
        dz[bidf[k]] = 0.0;
      for (k = 0; k < nub; ++k)
        dzu[k] = -xu[k] * zu[k] + alfp - phi * (dxu[k] * dzu[k] - rtmp);

      dkap = -opt->tau * opt->kap + alfp - phi * (dtau * dkap - rtmp);

      pid = GenSearchDir(opt, eta, tht, thtu, thtg, rp, rd, rg, pi, s, v, t, w,
                         dx, dxu, &dtau, dy, dyu, dz, dzu, &dkap, &rpnorm,
                         &rdnorm, &rgnorm, rsze, rwk);

      pid = pid && (rpnorm + rdnorm + rgnorm <= 10.0 * opt->smdt);
      if (!pid)
      {
        if (afine)
          continue;
        else
          return;
      }
    }

    if (pid)
    {
      pstp = RatioStep(opt, x, xu, dx, dxu);
      dstp = RatioStep(opt, z, zu, dz, dzu);

      if (dtau < -par->tolx)
      {
        pstp = min(pstp, -opt->tau / dtau);
        dstp = min(dstp, -opt->tau / dtau);
      }

      if (dkap < -par->tolx)
      {
        pstp = min(pstp, -opt->kap / dkap);
        dstp = min(dstp, -opt->kap / dkap);
      }
    }

    alfa = min(pstp, dstp);
    taup = opt->tau + pstp * dtau;
    taud = opt->tau + dstp * dtau;

    if (taup > taud)
    {
      btaup = true;
      alft = dstp;
      alfk = pstp;
    }
    else
    {
      btaup = false;
      alft = pstp;
      alfk = dstp;
    }

    if (alfa < 1.0e-13)
      return;

    else
    {
      rtmp = max(min(sqrt(par->beta), rfval / 100.0), par->beta);

      for (k = 0, boka = false; !boka; ++k)
      {
        rstep = par->Gamma * pow(0.95, (double)k);

        boka = ChkStep(opt, rg, eta, dx, dxu, dtau, dz, dzu, dkap, rstep * pstp,
                       rstep * alft, rstep * dstp, rstep * alfk, rtmp, &lpv);
      }
    }
    end = 1;
  }

  alfp = rstep * pstp;
  alft = rstep * alft;

  alfd = rstep * dstp;
  alfk = rstep * alfk;

  taup = opt->tau + alfp * dtau;
  taud = opt->tau + alfd * dtau;

  if (btaup)
  {
    pfac = taud / taup;
    dfac = 1.0;
  }

  else
  {
    pfac = 1.0;
    dfac = taup / taud;
  }

  for (k = 0; k < nlb; k++)
  {
    j = bid[k];
    x[j] += alfp * dx[j];
    z[j] += alfd * dz[j];
    x[j] *= pfac;
    z[j] *= dfac;
  }

  for (j = 0; j < m; j++)
  {
    y[j] += alfd * dy[j];
    y[j] *= dfac;
  }

  for (k = 0; k < nub; k++)
  {
    xu[k] += alfp * dxu[k];
    yu[k] += alfd * dyu[k];
    zu[k] += alfd * dzu[k];
    xu[k] *= pfac;
    yu[k] *= dfac;
    zu[k] *= dfac;
  }

  for (k = 0; k < nfr; k++)
  {
    j = bidf[k];
    x[j] += alfp * dx[j];
    x[j] *= pfac;
  }

  opt->tau += alft * dtau;
  opt->kap += alfk * dkap;
  opt->kap *= pfac * dfac;

} /* getInitPoint */

void HomOptSol(optmod *opt,
               optpar *par)
{
  int j, k, m, n, idir, nlb, nub, nfr, pid, boka, idch, isze, nchg,
      nlow, term, rsze, afine, btaup, rfac, rbac,
      xproc, *bid, *bidu, *bidf, *imem, *iwk;
  double r, mu, pi, rg, eta, lpv, mua, phi, alfa, alfd, alfk, alfp,
      alft, dfac, dkap, dstp, dtau, pfac, pstp, rtmp, taup, taud,
      thtg, gamma, rfval, rstep, rsum, rdnorm, rgnorm, rpnorm,
      *s, *v, *t, *w, *x, *y, *z, *xu, *yu, *zu, *tht, *thtu, *thtd,
      *dx, *dxu, *dz, *dzu, *dy, *dyu, *rp, *rd, *rmem, *rwk;
  matrix *at;
#ifdef TIME_COUNT
  clock_t Tm;
#endif
  term = false;
  m = opt->m;
  n = opt->n;
  nlb = opt->nlb;
  nub = opt->nub;
  nfr = opt->nfr;
  bid = opt->bid;
  bidu = bid + (nlb - nub);
  bidf = bid + nlb;

  dtau = 0.0;
  mu = 0.0;
  x = opt->x;
  xu = opt->xu;
  y = opt->y;
  yu = opt->yu;
  z = opt->z;
  zu = opt->zu;
  at = opt->at;

  isze = m + n + n;
  rsze = 8 * m + 13 * n + 12 * nub;

  imem = iAlloc(isze, "imem, HomOptSol.");
  rmem = dAlloc(rsze, "rmem, HomOptSol.");

  rp = rmem;
  rd = rp + m;
  tht = rd + n;
  thtu = tht + n;
  thtd = thtu + nub;

  s = thtd + n;
  v = s + n;
  t = v + nub;
  w = t + m;

  dx = w + nub;
  dxu = dx + n;
  dy = dxu + nub;
  dyu = dy + m;
  dz = dyu + nub;
  dzu = dz + n;

  iwk = imem;
  rwk = dzu + nub;

  ShowHead();
#ifdef TIME_COUNT
  Tm = GetTime();
#endif
  getInitPoint(opt, par, isze, imem, rsze, rmem);
#ifdef TIME_COUNT
  printf("getInitPoint: %d ms\n", GetTime() - Tm);
#endif

  opt->mu0 = (dDots(nlb, x, z, bid) + dDot(nub, xu, zu) + opt->tau * opt->kap) / opt->nall;
  opt->tau0 = opt->tau;
  opt->kap0 = opt->kap;
  mu = opt->mu0;
  opt->inf0 = opt->tau0 / opt->kap0;
  xproc = false;

  for (opt->iter = 0; opt->iter < par->iterlmt; opt->iter++)
  {
#ifdef TIME_COUNT
    Tm = GetTime();
#endif
    rfac = false;
    afine = false;

    opt->rho = max(mu / opt->mu0, 1.0e-12);

    for (k = 0; k < nfr; ++k)
      opt->xf[k] = 1.0 / (1.0 + 10.0 * fabs(opt->x[bidf[k]]));

    CompValues(opt, tht, thtu, &thtg, rp, rd, &rg,
               &mu, &rfval, isze, iwk, rsze, rwk);
#ifdef TIME_COUNT
    printf("Part1-A: %d ms\n", GetTime() - Tm);
#endif
    ShowInfor(opt);

    if ((opt->ista != ProcOk) ||
        (mu <= par->tolmu * opt->mu0 && opt->infe < par->tolinf) ||
        (opt->rgap <= par->tolgap && opt->rrpnm <= par->tolp && opt->rrdnm <= par->told))
      break;
#ifdef TIME_COUNT
    printf("Part1-B: %d ms\n", GetTime() - Tm);
#endif
    ClClean(opt->cl);
#ifdef TIME_COUNT
    printf("Part1-C: %d ms\n", GetTime() - Tm);
#endif
    idch = decChol(opt->cl, tht, at, isze, iwk, rsze, rwk); // *****
#ifdef TIME_COUNT
    printf("Part1-D: %d ms\n", GetTime() - Tm);
#endif

    if (idch != ChlOk)
    {
      ShutDown();
      exit(0);
    }
#ifdef TIME_COUNT
    printf("Part1: %d ms\n", GetTime() - Tm);
#endif
#ifdef TIME_COUNT
    Tm = GetTime();
#endif
    GetConst(opt, tht, thtu, &pi, s, v, t, w, rsze, rwk);
    phi = 1.0;
    rbac = true;

    for (idir = 0; idir < 4 && opt->ista == ProcOk && rbac && !rfac; ++idir)
    {
      rbac = false;
      if (phi != 0.0)
      {
        eta = 1.0;

        for (k = 0; k < nlb; ++k)
          dz[bid[k]] = -x[bid[k]] * z[bid[k]];
        for (k = 0; k < nfr; ++k)
          dz[bidf[k]] = 0.0;
        for (k = 0; k < nub; ++k)
          dzu[k] = -xu[k] * zu[k];

        dkap = -opt->tau * opt->kap;

        pid = GenSearchDir(opt, eta, tht, thtu, thtg, rp, rd, rg, pi, s, v, t, w,
                           dx, dxu, &dtau, dy, dyu, dz, dzu, &dkap, &rpnorm,
                           &rdnorm, &rgnorm, rsze, rwk);

        afine = pid;
        rsum = rpnorm * rdnorm + rgnorm;

        if (pid && (rpnorm + rdnorm + rgnorm <= opt->smdt / 10.0 || sqrt(rpnorm + rdnorm + rgnorm) <= 1.0e-8))
        {
          if (idir == 0)
          {
            setPartmsg(opt, dx, dxu, dtau, dz, dzu, dkap, opt->csk, &nchg, &nlow);

            pstp = RatioStep(opt, x, xu, dx, dxu);
            dstp = RatioStep(opt, z, zu, dz, dzu);

            if (dtau < -par->tolx)
            {
              pstp = min(pstp, -opt->tau / dtau);
              dstp = min(dstp, -opt->tau / dtau);
            }

            if (dkap < -par->tolx)
            {
              pstp = min(pstp, -opt->kap / dkap);
              dstp = min(dstp, -opt->kap / dkap);
            }

            rtmp = 1.0 - min(pstp, dstp);
            xproc = (par->cross && rtmp <= 0.1);

            mua = 0.0;
            for (k = 0; k < nlb; ++k)
              mua += (x[bid[k]] + pstp * dx[bid[k]]) * (z[bid[k]] + dstp * dz[bid[k]]);
            for (k = 0; k < nub; ++k)
              mua += (xu[k] + pstp * dxu[k]) * (zu[k] + dstp * dzu[k]);
            mua += (opt->tau + pstp * dtau) * (opt->kap + dstp * dkap);
            mua /= opt->nall;
            r = mua / mu;
            r = min(r * r * r, 0.2 * min(r, 1.0));
          }

          if (mu < 0.001)
            r = max(mu, r);
          r = max(1.0e-6, r);
        }

        else
        {
          fprintf(fout, " inaccurate search direction.\n");
          opt->ista = ProcStall;
          r = 0.1;
          phi = 0.0;
        }
      }

      gamma = r;
      eta = 1.0 - gamma;
      rtmp = (dDots(nlb, dx, dz, bid) + dDot(nub, dxu, dzu) + dkap * dtau) / opt->nall;
      rtmp = max(rtmp, 0.0);

      for (k = 0; k < nlb; ++k)
        dz[bid[k]] = -x[bid[k]] * z[bid[k]] + gamma * mu - phi * (dx[bid[k]] * dz[bid[k]] - rtmp);
      for (k = 0; k < nfr; ++k)
        dz[bidf[k]] = 0.0;
      for (k = 0; k < nub; ++k)
        dzu[k] = -xu[k] * zu[k] + gamma * mu - phi * (dxu[k] * dzu[k] - rtmp);

      dkap = -opt->tau * opt->kap + gamma * mu - phi * (dtau * dkap - rtmp);

      pid = GenSearchDir(opt, eta, tht, thtu, thtg, rp, rd, rg, pi, s, v, t, w,
                         dx, dxu, &dtau, dy, dyu, dz, dzu, &dkap, &rpnorm,
                         &rdnorm, &rgnorm, rsze, rwk);

      if (!pid || rpnorm + rdnorm + rgnorm > 10.0 * opt->smdt)
      {
        if (afine && rsum < opt->smdt / 10.0 && idir < 3)
        {
          r = 0.0;
          phi = 0.0;
          rbac = true;
        }

        else
        {
          if ((opt->rgap < par->tolgap) &&
                  (opt->rrpnm <= par->tolp * 100.0) &&
                  (opt->rrdnm <= par->told * 100.0) ||
              term)
          {
            fprintf(fout, "\n The solver is stalled\n");
            opt->ista = ProcStall;
            break;
          }

          term = true;
          rfac = true;

          rtmp = 100.0 * mu;
          for (k = 0; k < nlb; ++k)
          {
            j = bid[k];
            if (x[j] >= z[j])
              z[j] = rtmp / x[j];
            else
              x[j] = rtmp / z[j];
          }

          for (k = 0; k < nub; ++k)
          {
            if (xu[k] >= zu[k])
              zu[k] = rtmp / xu[k];
            else
              xu[k] = rtmp / zu[k];
          }
        }
      }

      if (opt->ista == ProcOk && !rfac && !rbac)
      {
        pstp = RatioStep(opt, x, xu, dx, dxu);
        dstp = RatioStep(opt, z, zu, dz, dzu);

        if (dtau < -par->tolx)
        {
          pstp = min(pstp, -opt->tau / dtau);
          dstp = min(dstp, -opt->tau / dtau);
        }

        if (dkap < -par->tolx)
        {
          pstp = min(pstp, -opt->kap / dkap);
          dstp = min(dstp, -opt->kap / dkap);
        }

        alfa = min(pstp, dstp);

        taup = opt->tau + pstp * dtau;
        taud = opt->tau + dstp * dtau;

        if (taup > taud)
        {
          btaup = true;
          alft = dstp;
          alfk = pstp;
        }
        else
        {
          btaup = false;
          alft = pstp;
          alfk = dstp;
        }

        if (alfa < 1.0e-13)
          opt->ista = UnkError;
        else
        {

          rtmp = max(min(sqrt(par->beta), rfval / 100.0), par->beta);

          for (k = 0, boka = false; !boka; ++k)
          {

            rstep = par->Gamma * pow(0.95, (double)k);

            boka = ChkStep(opt, rg, eta, dx, dxu, dtau,
                           dz, dzu, dkap, rstep * pstp,
                           rstep * alft, rstep * dstp,
                           rstep * alfk, rtmp, &lpv);
          }
        }
      }
    }

    if (opt->ista == ProcOk && !rfac && !rbac)
    {

      alfp = rstep * pstp;
      alft = rstep * alft;

      alfd = rstep * dstp;
      alfk = rstep * alfk;

      taup = opt->tau + alfp * dtau;
      taud = opt->tau + alfd * dtau;

      if (btaup)
      {
        pfac = taud / taup;
        dfac = 1.0;
      }

      else
      {
        pfac = 1.0;
        dfac = taup / taud;
      }

      for (k = 0; k < nlb; k++)
      {
        j = bid[k];
        x[j] += alfp * dx[j];
        z[j] += alfd * dz[j];
        x[j] *= pfac;
        z[j] *= dfac;
      }

      for (j = 0; j < m; j++)
      {
        y[j] += alfd * dy[j];
        y[j] *= dfac;
      }

      for (k = 0; k < nub; k++)
      {
        xu[k] += alfp * dxu[k];
        yu[k] += alfd * dyu[k];
        zu[k] += alfd * dzu[k];
        xu[k] *= pfac;
        yu[k] *= dfac;
        zu[k] *= dfac;
      }

      for (k = 0; k < nfr; k++)
      {
        j = bidf[k];
        x[j] += alfp * dx[j];
        x[j] *= pfac;
      }

      opt->tau += alft * dtau;
      opt->kap += alfk * dkap;

      opt->kap *= pfac * dfac;
    }
#ifdef TIME_COUNT
    printf("Part2: %d ms\n", GetTime() - Tm);
#endif
  } /* main loop */

  if (opt->iter >= par->iterlmt)
    opt->ista = ItrLimit;

  if (opt->ista != ProcOk)
  {
    for (k = 0; k < nlb - nub; ++k)
    {
      j = bid[k];

      if (z[j] <= x[j])
        opt->csk[j] = SPBAS;
      else
        opt->csk[j] = LOWER;
    }

    for (k = 0; k < nub; ++k)
    {
      j = bidu[k];
      if (z[j] <= x[j] && zu[k] <= xu[k])
        opt->csk[j] = SPBAS;
      else if (z[j] > x[j])
        opt->csk[j] = LOWER;
      else
        opt->csk[j] = UPPER;
    }

    for (k = 0; k < nfr; ++k)
      opt->csk[bidf[k]] = SPBAS;
  }

  free(imem);
  free(rmem);
} /* HomOptSol */

static int SolResult(optmod *opt,
                     optpar *par)
{
  int i, j, k, oka = true, n = opt->n,
               nlb = opt->nlb, nub = opt->nub,
               *bidu = opt->bid + (nlb - nub);
  double tau;

  if (opt->rrpnm <= par->tolp &&
      opt->rrdnm <= par->told &&
      opt->rgap <= par->tolgap &&
      opt->infe > par->tolinf)
  {
    opt->solsta = OPT_FOUND;
    opt->infsta = FEA_PROB;
#ifdef TEST
    fprintf(fres, "IOPT&");
#endif
  }
  else if (opt->infe <= par->tolinf)
  {
    opt->solsta = INF_PROB;
    opt->infsta = INF_PROB;
    oka = false;
#ifdef TEST
    fprintf(fres, "IINF&");
#endif
  }

  else if (opt->ista == ProcStall)
  {
    opt->solsta = NUM_PROB;
    opt->infsta = UNK_PROB;

#ifdef TEST
    fprintf(fres, "IUNK&");
#endif
  }

  else
  {
    opt->solsta = UNK_SOL;
    opt->infsta = UNK_PROB;
#ifdef TEST
    fprintf(fres, "IUNK&");
#endif
  }

  ErrorProc(opt->solsta, NULL);

  if (!oka)
    return false;

  tau = opt->tau;

  for (j = 0; j < n; ++j)
  {
    opt->x[j] /= tau;
    opt->z[j] /= tau;
  }

  for (i = 0; i < opt->m; ++i)
    opt->y[i] /= tau;

  for (k = 0; k < nub; ++k)
    opt->z[bidu[k]] += opt->yu[k] / tau;

  dFree(&opt->xu);
  dFree(&opt->yu);
  dFree(&opt->zu);

  iFree(&opt->cid);
  iFree(&opt->lperm);
  iFree(&opt->linvp);

  dFree(&opt->pcgbf);
  dPtFree(&opt->ld);

  if (opt->dac)
  {
    if (opt->dac->ia)
    {
      opt->dac->ia->ja = NULL;
      opt->dac->ia->an = NULL;
    }
    MtxFree(&opt->dac);
  }
  if (opt->sac)
  {
    if (opt->sac->ia)
    {
      opt->sac->ia->ja = NULL;
      opt->sac->ia->an = NULL;
    }
    MtxFree(&opt->sac);
  }

  return oka;
} /* SolResult */

void HsdProc(clock_t tim[],
             optmod *opt,
             optpar *par)
{
  int ok = true, m = opt->m,
      n = opt->n, nub = opt->nub;
  char ss[LineSize];
  FILE *fp;
  int i;
  clock_t tm;

  opt->nden = ChkDenseCol(opt, par);
  sprintf(ss, " " IFMT, opt->nden);
  LeftDots(ss, 51);
  printf("   dense columns " SFMT "\n\n", ss);

  opt->cl = ChlAlloc(m, "cl, HsdProc");
  ok = SymboProc(opt->cl, opt->a, opt->at);
  if (!ok)
    ErrorProc(SYB_PROC_ERR, "symbolic process err.");
#ifdef TEST1
  return;
#endif

  InitOptSol(opt);

#ifdef TIME_COUNT
  tm = GetTime();
#endif

  HomOptSol(opt, par);
#ifdef TIME_COUNT
  printf("HomOptSol: %d ms", GetTime() - tm);
#endif
  SetTime(tim, NUMER);

  ResList(true, par, opt);

  iFree(&opt->imem);
  dFree(&opt->rmem);

  ChlFree(&opt->cl);

  if (SolResult(opt, par) && par->cross)
  {
    CrossOverProc(opt);
    SetTime(tim, CROSS);
  }
} /* HsdProc */
