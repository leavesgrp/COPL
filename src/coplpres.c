#include "LPdefine.h"

static void PrsInit(optmod *opt, optpar *par)
{
  int m = opt->m, n = opt->n,
      nz = opt->nz, i, j, k, nn0;
  optrow *col = opt->col;
  matrix *at;

  dFree(&opt->c);
  dFree(&opt->l);
  dFree(&opt->u);

  opt->c = dAlloc(n, "opt->c, PrsInit");
  opt->l = dAlloc(n, "opt->l, PrsInit");
  opt->u = dAlloc(n, "opt->u, PrsInit");
  opt->rflg = cAlloc(m, "opt->rflg, PrsInit");
  opt->cflg = cAlloc(n, "opt->cflg, PrsInit");
  opt->at = MtxAlloc(n, nz, "opt->at, PrsInit");

  for (i = 0; i < m; i++)
    opt->rflg[i] = 'T';

  for (j = 0; j < n; j++)
    opt->cflg[j] = 'T';

  k = 0;
  at = opt->at;
  for (j = 0; j < n; j++)
  {
    nn0 = col[j].nn0;
    at->ia[j].ja = at->ia->ja + k;
    at->ia[j].an = at->ia->an + k;
    at->ia[j].nn0 = nn0;
    for (i = 0; i < nn0; i++)
    {
      at->ia[j].ja[i] = col[j].ja[i];
      at->ia[j].an[i] = col[j].an[i];
    }
    k += nn0;

    opt->c[j] = col[j].rhs;
    opt->l[j] = col[j].low;
    opt->u[j] = col[j].upp;
  }

  RowFree(&opt->col);

  MtxTrans(m, at, false, NULL, false, &opt->a);

} /* PrsInit */

static void FinalModel(optmod *opt,
                       optpar *par)
{
  int h, i, j, k, nlb, nub, nfr, m = opt->m, nn,
                                 n = opt->n, nz = opt->nz, nn0, *bid, *ja;
  double *c, *b, *l = opt->l, *u = opt->u, *an, xj,
                 cj;
  array *cols, *rows;
  char *rflg = opt->rflg, *cflg = opt->cflg;

  nlb = 0;
  nub = 0;
  nfr = 0;

  h = 0;
  nz = 0;
  cols = opt->at->ia;

  ja = cols->ja;
  an = cols->an;

  for (j = 0; j < n; j++)
  {

    if (cflg[j] != 'T')
      continue;

    nn0 = cols[j].nn0;
    nn = 0;
    for (k = 0; k < nn0; k++)
    {
      i = cols[j].ja[k];
      if (rflg[i] == 'T')
      {
        cols[j].ja[nn] = cols[j].ja[k];
        cols[j].an[nn] = cols[j].an[k];
        nn++;
      }
    }

    if (!nn)
    {

      cflg[j] = 'F';
      xj = 0.0;
      cj = opt->c[j];

      if (cj > 0.0)
        xj = opt->l[j];
      else if (cj < 0.0)
        xj = opt->u[j];

      if (fabs(xj) > par->bplus)
      {
        ErrorProc(UNB_PROB, NULL);
        ShutDown();
        exit(0);
      }

      opt->obj += cj * xj;

      AddStack();
      Int2Dbl2(CNUL, j, xj, cj);
    }
    else
    {
      cols[h].nn0 = nn;
      cols[h].ja = cols[j].ja;
      cols[h].an = cols[j].an;
      h++;
      nz += nn;
    }
  }

  AddStack();
  Int2(XSOL, n);

  h = 0;
  for (j = 0; j < n; j++)
  {
    if (cflg[j] != 'T')
      i = -1;

    else
    {
      i = h;
      h++;
    }

    Int1(i);
  }

  for (j = 0; j < n; j++)
  {

    if (cflg[j] != 'T')
      continue;

    if (u[j] <= par->bplus)
    {
      nub++;
      nlb++;
    }
    else if (l[j] >= -par->bplus)
      nlb++;
    else
      nfr++;
  }

  nn = nlb + nfr;
  opt->bid = iAlloc(nn, "opt->bid, FinalModel");
  l = dAlloc(nub, "opt->c, FinalModel");
  bid = opt->bid;

  h = 0;
  k = nlb - nub;
  i = nlb;
  nn = 0;
  nub = 0;

  for (j = 0; j < n; j++)
  {
    if (cflg[j] != 'T')
      continue;

    if (u[j] <= par->bplus)
    {
      bid[k] = nn;
      l[nub] = u[j];
      k++;
      nub++;
      nn++;
    }
    else if (opt->l[j] >= -par->bplus)
    {
      bid[h] = nn;
      h++;
      nn++;
    }
    else
    {
      bid[i] = nn;
      i++;
      nn++;
    }
  }

  dFree(&opt->l);
  dFree(&opt->u);
  dFree(&opt->r);

  opt->u = l;
  c = dAlloc(nn, "opt->c, FinalModel");
  nn = 0;
  for (j = 0; j < n; j++)
  {
    if (cflg[j] == 'T')
    {
      c[nn] = opt->c[j];
      nn++;
    }
  }

  dFree(&opt->c);
  opt->c = c;
  opt->n = nn;
  opt->at->nrow = nn;

  MtxTrans(m, opt->at, false, NULL, false, &opt->a);

  cols->ja = ja;
  cols->an = an;
  MtxFree(&opt->at);

  nn = 0;
  nn0 = 0;
  rows = opt->a->ia;
  ja = rows->ja;
  an = rows->an;

  for (i = 0; i < m; i++)
  {
    nn0 += rows[i].nn0;
    if (rows[i].nn0)
    {
      rows[nn].nn0 = rows[i].nn0;
      rows[nn].ja = rows[i].ja;
      rows[nn].an = rows[i].an;
      opt->b[nn] = opt->b[i];
      nn++;
    }
    else
    {
      if (opt->rflg[i] == 'T')
      {
        AddStack();
        Int2(RNUL, i);
      }
    }
  }

  AddStack();
  Int2(YSOL, m);

  nn = 0;
  for (i = 0; i < m; i++)
  {
    if (opt->rflg[i] == 'T')
    {
      j = nn;
      nn++;
    }
    else
      j = -1;

    Int1(j);
  }

  opt->a->nrow = nn;

  b = dAlloc(nn, "opt->b, FinalModel");
  for (i = 0; i < nn; i++)
    b[i] = opt->b[i];
  dFree(&opt->b);
  opt->b = b;

  MtxTrans(opt->n, opt->a, false, NULL, false, &opt->at);

  rows->ja = ja;
  rows->an = an;

  MtxFree(&opt->a);

  cFree(&opt->rflg);
  cFree(&opt->cflg);

  opt->m = nn;
  opt->nlb = nlb;
  opt->nub = nub;
  opt->nfr = nfr;
  opt->nz = nn0;

} /* FinalModel */

int PreSolve(optpar *par,
             optmod *opt)
{
  char ss[LineSize];
  int m, n, idproc, front, rear, isze, isze1,
      rsze, *nnzi, *nnzj, *iw1, *iw, *ix,
      iloop;
  double *rw;

  m = opt->m;
  n = opt->n;

  PrsInit(opt, par);

  if (par->prlev)
  {
    printf(" BEGIN presolving (level=" IFMT ")\n", par->prlev);

    isze = max(6 * m, n);
    isze = max(2 + 6 * m + n + isze, 1 + 2 * m + 3 * n);
    isze += 4 * m + 3 * n;
    rsze = 4 * m + n + n;

    opt->imem = iAlloc(isze, "opt->imem, prsreduce");
    opt->rmem = dAlloc(rsze, "opt->rmem, prsreduce");

    iw = opt->imem;
    rw = opt->rmem;
    ix = iw + 2 * (n + m);

    opt->rlb = rw;
    opt->rub = opt->rlb + m;
    opt->clb = opt->rub + m;
    opt->cub = opt->clb + n;

    nnzi = iw;
    nnzj = nnzi + opt->m;

    isze1 = isze - (opt->m + opt->n);
    iw1 = iw + (opt->m + opt->n);

    ChkRowLdpt(par, opt, isze1, iw1, rsze, rw);

    idproc = NullFix(par, opt);
    if (!idproc)
      return false;

    iloop = 0;

    do
    {
      opt->prtst = false;

      if (par->prlev <= 1)
        continue;

      idproc = ChkRowSing(m, par, opt, iw, ix);
      if (!idproc)
        return false;

      SetRowBnds(m, par, opt, iw, ix);
      idproc = ChkRowForc(m, par, opt, iw, ix);
      if (!idproc)
        return false;

      SetColList(n, &front, &rear, opt, iw, ix);
      idproc = ChkRowDomn(front, par, opt, ix);
      if (!idproc)
        return false;

      SetColList(n, &front, &rear, opt, iw, ix);
      idproc = ChkRowDblt(front, par, opt, ix);
      if (!idproc)
        return false;

      if (par->prlev <= 2)
        continue;

      SetColList(n, &front, &rear, opt, iw, ix);
      SetColBnds(m, front, par, opt, ix);
      idproc = ChkColDomn(n, rear, par, opt, ix);
      if (!idproc)
        return false;

      if (par->prlev <= 3)
        continue;

      SetColList(n, &front, &rear, opt, iw, ix);
      idproc = ChkRowDupl(m, n, rear, par, opt, iw, ix);
      if (!idproc)
        return false;

      if (par->prlev <= 4)
        continue;

      idproc = ChkColDupl(m, n, par, opt, iw, ix);
      if (!idproc)
        return false;

      iloop++;
    } while (opt->prtst && iloop < par->prlop);

    iFree(&opt->imem);
    dFree(&opt->rmem);

    sprintf(ss, " %d", opt->nrld);
    LeftDots(ss, 43);
    printf("   linear dependent rows " SFMT "\n", ss);

    sprintf(ss, " %d", opt->nrnul);
    LeftDots(ss, 55);
    printf("   null rows " SFMT "\n", ss);

    sprintf(ss, " %d", opt->nrsg);
    LeftDots(ss, 50);
    printf("   singleton rows " SFMT "\n", ss);

    sprintf(ss, " %d", opt->nrdb);
    LeftDots(ss, 50);
    printf("   doubleton rows " SFMT "\n", ss);

    sprintf(ss, " %d", opt->nrfc);
    LeftDots(ss, 52);
    printf("   forcing rows " SFMT "\n", ss);

    sprintf(ss, " %d", opt->nrdm);
    LeftDots(ss, 50);
    printf("   dominated rows " SFMT "\n", ss);

    sprintf(ss, " %d", opt->nrdu);
    LeftDots(ss, 50);
    printf("   duplicate rows " SFMT "\n", ss);

    sprintf(ss, " %d", opt->ncnul);
    LeftDots(ss, 52);
    printf("   null columns " SFMT "\n", ss);

    sprintf(ss, " %d", opt->ncfx);
    LeftDots(ss, 51);
    printf("   fixed columns " SFMT "\n", ss);

    sprintf(ss, " %d", opt->ncdm);
    LeftDots(ss, 47);
    printf("   domimated columns " SFMT "\n", ss);

    sprintf(ss, " %d", opt->ncdu);
    LeftDots(ss, 47);
    printf("   duplicate columns " SFMT "\n", ss);

    printf(" END presolving\n");
  }

  FinalModel(opt, par);

  printf("\n FINAL MODEL\n");

  if (opt->n > 0 && opt->m > 0)
  {
    sprintf(ss, " %.3f%s",
            100.0 * (((double)opt->nz /
                      (double)opt->m) /
                     (double)(opt->n)),
            "%");
    LeftDots(ss, 57);
    printf("   density " SFMT "\n", ss);

    sprintf(ss, " %.3f%s : " IFMT "/(" IFMT "*" IFMT ")",
            100.0 * (((double)opt->nz /
                      (double)opt->m) /
                     (double)opt->n),
            "%", opt->nz, opt->m, opt->n);
    LeftDots(ss, 53);
    fprintf(fout, "   final model " SFMT "\n\n", ss);
  }
  else
  {
    opt->pobj = opt->obj;
    opt->dobj = opt->obj;
    sprintf(ss, " " IFMT "/(" IFMT "*" IFMT ")",
            opt->nz, opt->m, opt->n);
    LeftDots(ss, 53);
    fprintf(fout, "   final model " SFMT "\n\n", ss);
  }

  sprintf(ss, " " IFMT, opt->m);
  LeftDots(ss, 60);
  printf("   rows " SFMT "\n", ss);

  sprintf(ss, " " IFMT, opt->n);
  LeftDots(ss, 57);
  printf("   columns " SFMT "\n", ss);

  sprintf(ss, " " IFMT, opt->nz);
  LeftDots(ss, 56);
  printf("   nonzeros " SFMT "\n", ss);

  sprintf(ss, " " IFMT, opt->nub);
  LeftDots(ss, 52);
  printf("   upper bounds " SFMT "\n", ss);

  sprintf(ss, " " IFMT, opt->nfr);
  LeftDots(ss, 52);
  printf("   free columns " SFMT "\n", ss);

#ifdef TEST1
  {
    int nc, i;

    nc = 0;
    for (i = 0; i < opt->n; i++)
      if (fabs(opt->c[i]) > par->aijtol)
        nc++;

    fprintf(fres, "&%6d&%6d&%7d&%6d",
            opt->m + 1, opt->n, opt->nz + nc, opt->nub);
  }
#endif

  return true;
} /* PreSolve */
