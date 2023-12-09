#include "LPdefine.h"

int ChkDenseCol(optmod *opt, optpar *par)
{
  int j, nden, idden = par->idden,
               n = opt->n, *cid;
  array *rows = opt->at->ia;

  if (!idden)
  {
    idden = opt->m;
    if (opt->m > 50000)
      idden = opt->m / 80;
    else if (opt->m > 30000)
      idden = opt->m / 50;
    else if (opt->m > 10000)
      idden = opt->m / 40;
    else if (opt->m > 1000)
      idden = opt->m / 10;
    else if (opt->m > 100)
      idden = opt->m / 5;
  }
  cid = iAlloc(n, "cid, ChkDenseCol");
  nden = 0;
  for (j = 0; j < n; j++)
  {
    cid[j] = false;
    if (rows[j].nn0 > idden)
    {
      /* Tian Xie 20160423
        nden++;
        cid[j]=true;
      */
    }
  }

  if (nden)
    opt->cid = cid;
  else
    iFree(&cid);

#ifdef TEST
  fprintf(fres, "&%4d", nden);
#endif

  return nden;
} /* ChkDenseCol */

void StripMatrix(optmod *opt)
{
  int j, sj = 0, dj = 0, n = opt->n,
         nden = opt->nden, *cid = opt->cid;
  array *rows = opt->at->ia,
        *drow = opt->dac->ia,
        *srow = opt->sac->ia;

  for (j = 0; j < n; j++)
  {
    if (cid[j])
    {
      drow[dj].nn0 = rows[j].nn0;
      drow[dj].ja = rows[j].ja;
      drow[dj].an = rows[j].an;
      dj++;
    }
    else
    {
      srow[sj].nn0 = rows[j].nn0;
      srow[sj].ja = rows[j].ja;
      srow[sj].an = rows[j].an;
      sj++;
    }
  }
} /* StripMatrix */

void StripDiag(optmod *opt,
               double *d,
               double *dd,
               double *sd)
{
  int j, n = opt->n,
         jd, js, *cid = opt->cid;

  jd = 0;
  js = 0;

  for (j = 0; j < n; j++)
  {
    if (!cid[j])
    {
      sd[js] = d[j];
      js++;
    }
    else
    {
      dd[jd] = d[j];
      jd++;
    }
  }
} /* StripDiag */

int LuPivot(int m,
            double **a,
            double *w,
            int *itmp,
            int *p,
            int *q)
{
  int g, h, i, j, k, mm;
  double temp, maxelm;

  for (i = 0; i < m; i++)
    p[i] = i;
  for (j = 0; j < m; j++)
    itmp[j] = j;

  mm = m - 1;
  for (k = 0; k < mm; k++)
  {

    /*
     * pivoting selection
     */
    maxelm = 0;
    for (i = k; i < m; i++)
    {
      for (j = k; j < m; j++)
      {
        temp = fabs(a[i][j]);
        if (temp > maxelm)
        {
          maxelm = temp;
          g = i;
          h = j;
        }
      }
    }

    if (maxelm < 1.0E-15)
    {
      printf("\n too small pivoting %3.1e.\n", maxelm);
      return false;
    }

    /*
     * interchange row k and row g
     * interchange col k and col h
     */
    if (g != k)
    {
      for (j = 0; j < m; j++)
      {
        temp = a[g][j];
        a[g][j] = a[k][j];
        a[k][j] = temp;
      }
      i = p[k];
      p[k] = p[g];
      p[g] = i;
    }

    if (h != k)
    {
      for (i = 0; i < m; i++)
      {
        temp = a[i][h];
        a[i][h] = a[i][k];
        a[i][k] = temp;
      }
      j = itmp[k];
      itmp[k] = itmp[h];
      itmp[h] = j;
    }

    /*
     * lu decomposition
     */
    for (j = k + 1; j < m; j++)
      w[j] = a[k][j];
    for (i = k + 1; i < m; i++)
    {
      temp = a[i][k] / a[k][k];
      a[i][k] = temp;
      for (j = k + 1; j < m; j++)
        a[i][j] = a[i][j] - temp * w[j];
    }
  }

  for (j = 0; j < m; j++)
  {
    k = itmp[j];
    q[k] = j;
  }

  return true;
} /* LuPivot */

void DenseMatrix(optmod *opt,
                 double *dd,
                 double *rmem,
                 int *imem)
{
  int n = opt->n, nden = opt->nden,
      nspa = n - nden, m = opt->m,
      i, j, k, nn0;
  double rtmp, **r = NULL, *rhs, *sol;
  matrix *as = opt->sac, *ad = opt->dac;

  opt->ld = dPtAlloc(nden, nden, "opt->ld, DenseMatrix");

  r = opt->ld;
  rhs = rmem;
  sol = rhs + m;

  for (j = 0; j < nden; j++)
  {
    for (k = 0; k < m; k++)
      rhs[k] = 0.0;

    nn0 = ad->ia[j].nn0;
    for (k = 0; k < nn0; k++)
      rhs[ad->ia[j].ja[k]] = ad->ia[j].an[k];

    CholSol(opt->cl, rhs, sol);

    rtmp = dd[j];
    for (k = 0; k < nn0; k++)
      rtmp += ad->ia[j].an[k] * sol[ad->ia[j].ja[k]];
    r[j][j] = rtmp;

    if (j > 0)
    {
      for (i = 0; i < j; i++)
      {
        nn0 = ad->ia[i].nn0;
        rtmp = 0.0;
        for (k = 0; k < nn0; k++)
          rtmp += ad->ia[i].an[k] * sol[ad->ia[i].ja[k]];
        r[j][i] = rtmp;
        r[i][j] = rtmp;
      }
    }
  }

  if (!LuPivot(nden, r, rhs, imem,
               opt->lperm, opt->linvp))
  {
    ErrorProc(FAC_ERROR, NULL);
    ShutDown();
    exit(0);
  }
} /* DenseMatrix */

static int LesSolver(int size,
                     double **a,
                     double *b,
                     double *x,
                     int *p,
                     int *q)
{
  int i, j;

  for (i = 0; i < size; i++)
  {
    if (fabs(a[i][i]) < 1.E-15)
    {
      printf("\n zero diagonal found\n");
      ShutDown();
      exit(0);
    }
    x[i] = b[p[i]];
  }

  /*
   * forward substitution
   */
  for (i = 0; i < size; i++)
    for (j = 0; j < i; j++)
      x[i] = x[i] - a[i][j] * x[j];

  /*
   * backward substitution
   */
  for (i = size - 1; i >= 0; i--)
  {
    for (j = i + 1; j < size; j++)
      x[i] = x[i] - a[i][j] * x[j];
    x[i] = x[i] / a[i][i];
  }

  for (i = 0; i < size; i++)
    b[i] = x[q[i]];

  return true;
} /* LesSolver */

static void GetResidual(optmod *opt,
                        double *diag,
                        double *b,
                        double *x,
                        double *tmpx,
                        double *tmpy)
{
  int i, m = opt->m,
         n = opt->n;
  matrix *at = opt->at;

  mTimesv(true, m, n, 1.0, at, x, 0.0, tmpx);

  for (i = 0; i < n; i++)
    tmpx[i] /= diag[i];

  for (i = 0; i < m; i++)
    tmpy[i] = b[i];

  mTimesv(false, m, n, -1.0, at, tmpx, 1.0, tmpy);
} /* GetResidual */

static int PreConjGrad(optmod *opt,
                       double *rsd,
                       double *ptol,
                       double nrhs,
                       int iter,
                       double *sol,
                       double *diag)
{
  int i, k, m = opt->m,
            n = opt->n;
  double alpha, beta, ztrk, ztr0, papk,
      *r, *z, *p, *t, *y, *x, norm;
  matrix *at = opt->at;

  opt->pcgbf = dAlloc(3 * m + n, "opt->pcgbf, PreConjGrad");
  r = opt->pcgbf;
  z = r + m;
  p = z + m;
  t = p + m;

  y = rsd;
  x = sol;

  norm = 0.0;
  for (i = 0; i < opt->m; i++)
  {
    r[i] = y[i];
    norm = norm + r[i] * r[i];
  }
  norm = sqrt(norm);

  for (k = 0; k < iter && norm >= (*ptol); k++)
  {

    /*
     * preconditioning residual
     */
    CholSol(opt->cl, y, z);

    /*
     * compute beta^k and p^k
     */
    if (!k)
    {
      ztrk = 0.0;
      for (i = 0; i < m; i++)
      {
        p[i] = z[i];
        ztrk = ztrk + r[i] * z[i];
      }
    }

    else
    {
      ztr0 = ztrk;
      ztrk = 0.0;
      for (i = 0; i < opt->m; i++)
        ztrk = ztrk + r[i] * z[i];

      if (fabs(ztr0) < 1.0E-15 * ztrk)
      {
        fprintf(fout, "   preconditioner "
                      "conjugate grandient"
                      " failed\n");
        dFree(&opt->pcgbf);
        return false;
      }

      beta = ztrk / ztr0;
      for (i = 0; i < m; i++)
        p[i] = z[i] + beta * p[i];
    }

    /*
     * compute y^k=ADA^T*p^k
     */

    mTimesv(true, m, n, 1.0, at, p, 0.0, t);

    for (i = 0; i < n; i++)
      t[i] /= diag[i];

    mTimesv(false, m, n, 1.0, at, t, 0.0, y);

    /*
     * compute alpha=z^T*r/p^T*y;
     */
    papk = 0.0;
    for (i = 0; i < m; i++)
      papk = papk + p[i] * y[i];

    if (fabs(papk) < 1.0E-15 * ztrk)
    {
      fprintf(fout, "   preconditioner "
                    "conjugate grandient"
                    " failed\n");
      dFree(&opt->pcgbf);
      return false;
    }

    alpha = ztrk / papk;

    /*
     * update x^k and r^k
     */
    norm = 0.0;
    for (i = 0; i < opt->m; i++)
    {
      x[i] = x[i] + alpha * p[i];
      r[i] = r[i] - alpha * y[i];
      y[i] = r[i];
      norm = norm + r[i] * r[i];
    }
    norm = sqrt(norm);
  }

  fprintf(fout, " PCG reduced it to %5.1e in %d iterations\n",
          norm / nrhs, k + 1);
  if (k >= iter)
    fprintf(fout,
            "     (Max PCG iterations exceeded,"
            " use the latest answer ...)\n");
  *ptol = norm;

  dFree(&opt->pcgbf);
  return true;

} /* PreConjGrad */

static void LuSolve(optmod *opt,
                    double *diag,
                    double **ld,
                    double *b,
                    double *x)
{
  int i, j, k, nn0, nall, nden, nrnk,
      m = opt->m;
  double *dtmp, normb, normr, *z,
      normr0, pcgrf, *p, *q, rtmp;
  matrix *dac = opt->dac, *sac = opt->sac;
  chol *kf = opt->cl;

  nrnk = 0;
  nden = opt->nden;

  nall = nden + nrnk;
  dtmp = dAlloc(nall, "dtmp, LuSolve");
  if (m > nall)
    z = dAlloc(m, "z, LuSolve");
  else
    z = dAlloc(nall, "z, LuSolve");
  p = dtmp;
  q = p + nden;

  /*
   * normalize b
   */
  normb = 0.0;
  for (i = 0; i < m; i++)
  {
    normb += fabs(b[i]);
    z[i] = b[i];
  }
  if (normb < 1.0E-15)
  {
    for (i = 1; i < m; i++)
      x[i] = 0.0;
    normb = 1.0;
  }
  else
  {
    normb = log(normb) / log(2.0);
    k = (int)(normb + 0.5);
    normb = pow(2.0, (double)k);
    for (i = 0; i < m; i++)
      z[i] = z[i] / normb;

    /*
     * sol=(l*d*lt)^-1*rhs
     */
    CholSol(kf, z, x);
  }

  /*
   * p=Ad^T*sol
   */
  for (i = 0; i < nden; i++)
  {
    rtmp = 0.0;
    nn0 = dac->ia[i].nn0;
    for (k = 0; k < nn0; k++)
      rtmp += dac->ia[i].an[k] * x[dac->ia[i].ja[k]];
    p[i] = rtmp;
  }

  /*
   * normalize p
   */
  normr = 0.0;
  for (i = 0; i < nden; i++)
    normr += fabs(p[i]);

  rtmp = log(normr) / log(2.0);
  k = (int)(rtmp + 0.5);
  normr = pow(2.0, (double)k);

  for (i = 0; i < nall; i++)
    p[i] /= normr;

  /*
   * solve dense equation system
   */
  LesSolver(nden, opt->ld, p, z,
            opt->lperm, opt->linvp);

  /*
   * undo normalization
   */
  rtmp = normb * normr;
  for (i = 0; i < nden; i++)
    p[i] *= rtmp;

  /*
   * z=b-Ad*p
   */
  for (i = 0; i < m; i++)
    z[i] = b[i];

  for (j = 0; j < nden; j++)
  {
    nn0 = dac->ia[j].nn0;
    for (k = 0; k < nn0; k++)
      z[dac->ia[j].ja[k]] -= dac->ia[j].an[k] * p[j];
  }

  /*
   * normalize z, sol=(l*d*lt)^-1*z
   */
  normb = 0.0;
  for (i = 0; i < m; i++)
    normb += fabs(z[i]);
  if (normb < 1.0E-15)
  {
    for (i = 0; i < m; i++)
      x[i] = 0.0;
    normb = 1.0;
  }
  else
  {
    rtmp = log(normb) / log(2.0);
    k = (int)(rtmp + 0.5);
    normb = pow(2.0, (double)k);

    for (i = 0; i < m; i++)
      z[i] /= normb;

    CholSol(kf, z, x);
  }
  for (i = 0; i < m; i++)
    x[i] *= normb;

  dFree(&dtmp);

  /*
   * improve accuracy if necessary
   */
  k = max(opt->n, m);
  dtmp = dAlloc(k, "dtmp, LuSolve");
  GetResidual(opt, diag, b, x, dtmp, z);

  normb = 0.0;
  normr = 0.0;
  for (i = 0; i < m; i++)
  {
    normb += b[i] * b[i];
    normr += z[i] * z[i];
  }

  normb = sqrt(normb);
  normr = sqrt(normr);

  pcgrf = normb * pr->pcgtol;

  if (normr > pcgrf)
  {
    fprintf(fout, " relative residual = " E3FT "\n", normr / normb);
    normr0 = normr;
    for (i = 0; i < m; i++)
      dtmp[i] = x[i];

    PreConjGrad(opt, z, &pcgrf, normb, 10 * opt->nden, x, diag);

    if (pcgrf > normr0)
    {
      fprintf(fout, " no improvement achieved.\n");
      for (i = 0; i < m; i++)
        x[i] = dtmp[i];
    }
  }
  dFree(&z);
  dFree(&dtmp);
} /* LuSolve */

void DensNormSol(optmod *opt,
                 chol *kf,
                 double *diag,
                 matrix *at,
                 double *b,
                 double *x)
{
  int j, m = opt->m, n = opt->n;

  for (j = 0; j < n; ++j)
    x[j] = b[j] / diag[j];

  if (at)
    mTimesv(false, m, n, -1.0, at, x, 1.0, b + n);

  dCopy(n, b, x);

  LuSolve(opt, diag, opt->ld, b + n, x + n);

  if (at)
    mTimesv(true, m, n, 1.0, at, x + n, 1.0, x);

  for (j = 0; j < m; ++j)
    x[n + j] = -x[n + j];

  for (j = 0; j < n; ++j)
    x[j] /= diag[j];
} /* DensNormSol */
