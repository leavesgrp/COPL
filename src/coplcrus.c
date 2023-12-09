#include "LPdefine.h"

int DataTransf(optmod *opt, optpar *par)
{
  int i, j, k, beg, end, nz, *ja,
      n = opt->n, *ind = opt->imem;
  double *l = opt->l, *u = opt->u,
         *c = opt->rmem, *an;
  optrow *col = opt->col;

  ja = ind + par->mcol + 1;
  an = c + par->mcol;

  AddStack();
  Int3(MATX, opt->nz, n);
  fwrite(c, 8, n, fbin);
  opt->bytes += 8 * n;

  for (j = 0; j < n; j++)
  {
    beg = ind[j];
    end = ind[j + 1];
    nz = end - beg;

    col[j].nn0 = nz;
    col[j].rhs = c[j];
    col[j].low = l[j];
    col[j].upp = u[j];

    col[j].ja = col->ja + beg;
    col[j].an = col->an + beg;

    for (k = beg; k < end; k++)
    {
      i = k - beg;
      col[j].ja[i] = ja[k];
      col[j].an[i] = an[k];
    }

    Int1(nz);
    fwrite(col[j].ja, 4, nz, fbin);
    fwrite(col[j].an, 8, nz, fbin);
    opt->bytes += 12 * nz;
  }

  iFree(&opt->imem);
  dFree(&opt->rmem);
  return true;
} /* DataTransf */

int TreatBound(double bound, int indx, optpar *par, optmod *opt, int type)
{
  int i, j, k;
  double rtmp;

  switch (type)
  {
  case LOW:
    AddStack();
    Int2Dbl1(LBND, indx, bound);

    k = opt->col[indx].nn0;
    for (j = 0; j < k; j++)
    {
      i = opt->col[indx].ja[j];
      rtmp = -opt->col[indx].an[j] * bound;
      opt->b[i] += rtmp;
      opt->r[i] += rtmp;
    }

    opt->obj += opt->col[indx].rhs * bound;
    if (opt->col[indx].upp <= par->bplus)
      opt->col[indx].upp -= bound;
    opt->col[indx].low = 0.0;
    break;

  case MIS:
    AddStack();
    Int2Dbl1(MINS, indx, bound);

    k = opt->col[indx].nn0;
    opt->col[indx].low = 0;
    opt->col[indx].upp = MaxPositive;
    opt->col[indx].rhs = -opt->col[indx].rhs;

    if (fabs(bound) >= par->aijtol)
    {
      for (j = 0; j < k; j++)
      {
        i = opt->col[indx].ja[j];
        opt->col[indx].an[j] = -opt->col[indx].an[j];
        rtmp = opt->col[indx].an[j] * bound;
        opt->b[i] += rtmp;
        opt->r[i] += rtmp;
      }
      opt->obj += opt->c[indx] * bound;
    }
    else
      for (j = 0; j < k; j++)
        opt->col[indx].an[j] = -opt->col[indx].an[j];
    break;

  default:
    break;
  }
  return true;
} /* TreatBound */

int CrushProc(optpar *par, optmod *opt)
{
  int i, j, k, n, nz;
  double lj, uj;
  char ss[LineSize];
  optrow *col;

  printf(" BEGIN crushing...\n");

  n = opt->n + opt->nsk;
  nz = opt->nz + opt->nsk;
  opt->col = RowAlloc(n, nz, "col, CrushProc");

  DataTransf(opt, par);
  col = opt->col;

  /*
   *  standardize bounds, i.e. 0<=x<=u.
   */
  n = opt->n;
  for (j = 0; j < n; j++)
  {

    lj = col[j].low;
    uj = col[j].upp;

    if (lj > uj)
    {
      sprintf(ss, "inconsistent bounds on var %s.", opt->colname[j]);
      return ErrorProc(INF_PROB, ss);
    }

    else if (fabs(lj) < par->aijtol)
      continue;

    else if (lj >= -par->bplus)
      TreatBound(lj, j, par, opt, LOW);

    else if (uj <= par->bplus)
      TreatBound(uj, j, par, opt, MIS);
  }

  /*
   * standardize constriants, i.e. Ax=b.
   */
  nz = opt->nz;
  for (i = 0; i < opt->m; i++)
  {

    if (opt->rowtype[i] == 'E')
      continue;

    j = n;
    k = j - 1;
    col[j].nn0 = 1;
    col[j].rhs = 0.0;
    col[j].low = 0.0;
    col[j].upp = MaxPositive;

    col[j].ja = col->ja + nz;
    col[j].an = col->an + nz;

    col[j].ja[0] = i;
    n++;
    nz++;

    switch (opt->rowtype[i])
    {
    case 'L':
      col[j].an[0] = 1.0;
      break;

    case 'G':
      col[j].an[0] = -1.0;
      break;

    case 'R':
      col[j].an[0] = -1.0;
      col[j].upp = opt->b[i] - opt->r[i];
      opt->b[i] = opt->r[i];
      break;

    default:
      break;
    }
  }

  opt->nz = nz;
  opt->n = n;

  printf(" END crushing\n");
  printf(" DENSITY\n");
  sprintf(ss, " %.3f%s : " IFMT "/(" IFMT "*" IFMT ")",
          100.0 * (((double)opt->nz /
                    (double)opt->m) /
                   (double)opt->n),
          "%", opt->nz, opt->m, opt->n);
  LeftDots(ss, 50);
  printf("   after crushing " SFMT "\n\n", ss);
  fprintf(fout, "   after crushing " SFMT "\n", ss);

  return true;
} /* CrushProc */
