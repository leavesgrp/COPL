#include "LPdefine.h"

int NullFix(optpar *par, optmod *opt)
{
  int m, n, j;
  array *aj;
  char *rflg, *cflg;

  m = opt->m;
  n = opt->n;
  rflg = opt->rflg;
  cflg = opt->cflg;

  for (j = 0; j < n; j++)
  {

    if (opt->u[j] < opt->l[j])
    {
      ErrorProc(INF_PROB, " inconsistant bounds.\n");
      return false;
    }

    if (opt->u[j] - opt->l[j] <= par->aijtol)
    {
      aj = opt->at->ia + j;

      AddStack();
      Int2Dbl2(CFIX, j, opt->u[j], opt->c[j]);
      RecVect(aj, false);

      opt->cflg[j] = 'F';
      opt->ncfx++;
    }

    else
    {
      aj = opt->at->ia + j;
      if (!aj->nn0)
      {
        if (opt->c[j] < 0.0)
        {
          if (opt->u[j] > par->bplus)
          {
            ErrorProc(UNB_PROB, NULL);
            return false;
          }
          else
          {
            AddStack();
            Int2Dbl2(CNUL, j, opt->u[j], opt->c[j]);
            opt->obj += opt->u[j] * opt->c[j];
          }
        }

        else
        {
          if (opt->l[j] < -par->bplus)
          {
            ErrorProc(UNB_PROB, NULL);
            return false;
          }
          else
          {
            AddStack();
            Int2Dbl2(CNUL, j, opt->l[j], opt->c[j]);
            opt->obj += opt->c[j] * opt->l[j];
          }
        }

        opt->cflg[j] = 'F';
        opt->ncnul++;
      }
    }
  }

  return true;
} /* NullFix */

int ChkRowSing(int m,
               optpar *par,
               optmod *opt,
               int *iw,
               int *ix)
{
  int h, i, j, k, front, rear;
  double xx, aij, lj, uj;
  array *ai, *aj;

  ix = iw + m;
  front = -1;
  rear = -1;

  for (i = 0; i < m; i++)
  {
    if (opt->rflg[i] != 'T')
      continue;

    iw[i] = 0;
    ai = opt->a->ia + i;
    for (k = 0; k < ai->nn0; k++)
    {
      j = ai->ja[k];
      if (opt->cflg[j] == 'T')
        iw[i]++;
    }
  }

  for (i = 0; i < m; i++)
  {
    if (opt->rflg[i] != 'T')
      continue;
    ix[i] = -1;

    if (iw[i] <= 1)
    {
      if (front < 0)
      {
        front = i;
        rear = i;
      }
      else
      {
        ix[rear] = i;
        rear = i;
      }
    }
  }

  for (i = front; i >= 0; i = ix[i])
  {

    opt->rflg[i] = 'F';
    opt->nrsg++;

    ai = opt->a->ia + i;
    for (k = 0; k < ai->nn0; k++)
    {
      j = ai->ja[k];
      if (opt->cflg[j] == 'T')
        break;
    }

    if (k >= ai->nn0)
    {
      if (fabs(opt->b[i]) >= par->aijtol)
      {
        ErrorProc(INF_PROB, NULL);
        return false;
      }
      AddStack();
      Int2(RNUL, i);
    }

    else
    {
      aij = ai->an[k];

      if (fabs(aij) < par->aijtol)
      {
        if (fabs(opt->b[i]) >= par->aijtol)
        {
          ErrorProc(INF_PROB, NULL);
          return false;
        }
        else
        {
          AddStack();
          Int2(RNUL, i);
        }
      }

      else
      {
        xx = opt->b[i] / aij;
        lj = opt->l[j];
        uj = opt->u[j];

        opt->cflg[j] = 'F';
        opt->ncfx++;

        opt->prtst = true;

        if (xx < lj || xx > uj)
        {
          ErrorProc(INF_PROB, NULL);
          return false;
        }

        aj = opt->at->ia + j;

        AddStack();
        Int3Dbl3(RSNG, i, j, xx, aij, opt->c[j]);
        RecVect(aj, false);

        for (k = 0; k < aj->nn0; k++)
        {
          h = aj->ja[k];
          if (opt->rflg[h] != 'T')
            continue;

          iw[h]--;
          opt->b[h] -= xx * aj->an[k];

          if (iw[h] == 1)
          {
            ix[rear] = h;
            rear = h;
          }
        }

        opt->obj += opt->c[j] * xx;
      }
    }
  }

  return true;
} /* ChkRowSing */

int SetRowBnds(int m,
               optpar *par,
               optmod *opt,
               int *iw,
               int *ix)
{
  int i, j, k, ig, ih;
  double xx, li, lj, ui, uj;
  array *ai;

  for (i = 0; i < m; i++)
  {

    if (opt->rflg[i] != 'T')
      continue;

    li = 0.0;
    ui = 0.0;
    ig = 0;
    ih = 0;

    ai = opt->a->ia + i;

    for (k = 0; k < ai->nn0; k++)
    {
      j = ai->ja[k];
      if (opt->cflg[j] != 'T')
        continue;

      xx = ai->an[k];
      uj = opt->u[j];
      lj = opt->l[j];
      if (xx > 0)
      {
        if (uj > par->bplus)
          ih++;
        else
          ui += xx * uj;
        if (lj < -par->bplus)
          ig++;
        else
          li += xx * lj;
      }
      else
      {
        if (uj > par->bplus)
          ig++;
        else
          li += xx * uj;
        if (lj < -par->bplus)
          ih++;
        else
          ui += xx * lj;
      }
    }

    opt->rlb[i] = li;
    opt->rub[i] = ui;
    ix[i] = ig;
    iw[i] = ih;
  }

  return true;
} /* SetRowBnds */

int ChkRowForc(int m,
               optpar *par,
               optmod *opt,
               int *iw,
               int *ix)
{
  int i, j, k, p, q;
  double xx, li, ui, lj, aij;
  char ch = '&';
  array *ai, *aj;

  for (i = 0; i < m; i++)
  {

    if (opt->rflg[i] != 'T')
      continue;

    if (ix[i])
      li = -MaxPositive;
    else
      li = opt->rlb[i];

    if (iw[i])
      ui = MaxPositive;
    else
      ui = opt->rub[i];

    xx = opt->b[i];
    if (li > xx || ui < xx)
    {
      ErrorProc(INF_PROB, NULL);
      return false;
    }

    if (xx - li >= par->aijtol && ui - xx >= par->aijtol)
      continue;

    ai = opt->a->ia + i;

    opt->rflg[i] = 'F';
    opt->prtst = true;
    opt->nrfc++;

    AddStack();
    Int2(RFRC, i);

    for (k = 0; k < ai->nn0; k++)
    {
      j = ai->ja[k];
      if (opt->cflg[j] != 'T')
        continue;

      lj = opt->l[j];
      aij = ai->an[k];

      if ((aij > 0 && (ui - xx) < par->aijtol) || (aij < 0 && (xx - li) < par->aijtol))
      {
        lj = opt->u[j];
        Char1('u');
      }
      else
        Char1('l');

      aj = opt->at->ia + j;
      Int1Dbl3(j, lj, aij, opt->c[j]);
      RecVect(aj, false);

      for (p = 0; p < aj->nn0; p++)
      {
        q = aj->ja[p];
        if (opt->rflg[q] != 'T' || q == i)
          continue;

        aij = aj->an[p];
        opt->b[q] -= aij * lj;
      }
      opt->obj += opt->c[j] * lj;
      opt->cflg[j] = 'F';
      opt->ncfx++;
    }

    Char1('e');
  }
  return true;
} /* ChkRowForc */

int SetColList(int n,
               int *front,
               int *rear,
               optmod *opt,
               int *iw,
               int *ix)
{
  int i, j, k;
  array *aj;

  *front = -1;
  *rear = -1;

  for (j = n - 1; j >= 0; j--)
  {

    iw[j] = 0;
    ix[j] = -1;

    if (opt->cflg[j] != 'T')
      continue;

    aj = opt->at->ia + j;
    for (k = 0; k < aj->nn0; k++)
    {
      i = aj->ja[k];
      if (opt->rflg[i] != 'T')
        continue;

      iw[j]++;
      if (iw[j] > 2)
        break;
    }

    if (iw[j] != 1)
      continue;

    if ((*front) < 0)
    {
      *front = j;
      *rear = j;
    }

    else
    {
      ix[*rear] = j;
      *rear = j;
    }
  }

  return true;
} /* SetCList */

int ChkRowDomn(int front,
               optpar *par,
               optmod *opt,
               int *ix)
{
  int i, j, k, ig, p, q, r, ih;
  double aij, lj, uj, gi, bi,
      cj, li, ui, hi, aik;
  array *aj, *ai;

  for (j = front; j >= 0; j = ix[j])
  {
    aj = opt->at->ia + j;
    for (k = 0; k < aj->nn0; k++)
    {
      i = aj->ja[k];
      if (opt->rflg[i] == 'T')
        break;
    }

    if (k >= aj->nn0)
      continue;

    aij = aj->an[k];
    if (fabs(aij) < par->aijtol)
      return ErrorProc(OVER_FLOW, NULL);

    lj = opt->l[j];
    uj = opt->u[j];
    cj = opt->c[j];
    bi = opt->b[i];
    ig = 0;
    ih = 0;
    ui = bi;
    li = bi;

    ai = opt->a->ia + i;

    for (p = 0; p < ai->nn0; p++)
    {
      q = ai->ja[p];
      if (opt->cflg[q] != 'T' || j == q)
        continue;

      hi = ai->an[p];
      if (hi > 0)
      {
        if (!ig)
        {
          if (opt->u[q] > par->bplus)
          {
            li = -MaxPositive;
            ig = 1;
          }
          else
            li -= hi * opt->u[q];
        }
        if (!ih)
        {
          if (opt->l[q] < -par->bplus)
          {
            ui = MaxPositive;
            ih = 1;
          }
          else
            ui -= hi * opt->l[q];
        }
      }

      else
      {
        if (!ig)
        {
          if (opt->l[q] < -par->bplus)
          {
            li = -MaxPositive;
            ig = 1;
          }
          else
            li -= hi * opt->l[q];
        }
        if (!ih)
        {
          if (opt->u[q] > par->bplus)
          {
            ui = MaxPositive;
            ih = 1;
          }
          else
            ui -= hi * opt->u[q];
        }
      }
    }

    if (li < -par->bplus)
    {
      li = -MaxPositive;
      if (aij < 0.0)
        li = -li;
    }
    else
      li = li / aij;

    if (ui > par->bplus)
    {
      ui = MaxPositive;
      if (aij < 0.0)
        ui = -ui;
    }
    else
      ui = ui / aij;

    gi = li;
    li = min(li, ui);
    ui = max(ui, gi);

    if (li < lj || ui > uj)
    {
      opt->u[j] = min(ui, uj);
      if (li > lj)
      {
        for (; k < aj->nn0; k++)
        {
          r = aj->ja[k];
          if (opt->rflg[r] != 'T')
            continue;
          opt->b[r] -= aj->an[k] * li;
        }
        opt->obj += cj * li;
        if (opt->u[j] < par->bplus)
          opt->u[j] -= li;
        opt->l[j] = 0.0;

        AddStack();
        Int2Dbl1(LBND, j, li);
      }
      continue;
    }

    opt->prtst = true;
    opt->cflg[j] = 'F';
    opt->rflg[i] = 'F';
    opt->ncfx++;
    opt->nrdm++;

    AddStack();
    Int3Dbl3(RDMN, i, j, opt->b[i], aij, cj);
    RecVect(ai, true);

    aij = cj / aij;

    for (p = 0; p < ai->nn0; p++)
    {
      q = ai->ja[p];
      if (opt->cflg[q] != 'T' || j == q)
        continue;

      aik = ai->an[p];
      opt->c[q] -= aij * aik;
    }

    opt->obj += aij * bi;
  }
  return true;
} /* ChkRowDomn */

int ChkRowDblt(int front,
               optpar *par,
               optmod *opt,
               int *ix)
{
  int i, j, k, jmin, jmax, p, q, r;
  double aij, lj, uj, bi, aik, lk, uk,
      apk, cj;
  array *aj, *ai, *aq;

  for (j = front; j >= 0; j = ix[j])
  {

    if (opt->cflg[j] != 'T')
      continue;

    aj = opt->at->ia + j;
    for (k = 0; k < aj->nn0; k++)
    {
      i = aj->ja[k];
      if (opt->rflg[i] == 'T')
        break;
    }

    if (k >= aj->nn0)
      continue;

    jmin = k;
    ai = opt->a->ia + i;

    r = 0;
    q = 0;

    for (k = 0; k < ai->nn0; k++)
    {
      p = ai->ja[k];
      if (opt->cflg[p] != 'T')
        continue;

      r++;
      if (!q && p != j)
      {
        q = p;
        jmax = k;
      }
      if (r > 2)
        break;
    }

    if (r != 2)
      continue;

    aij = aj->an[jmin];
    aik = ai->an[jmax];

    if ((fabs(aij) < par->aijtol) || (fabs(aik) < par->aijtol))
      return ErrorProc(OVER_FLOW, NULL);

    lj = opt->l[j];
    uj = opt->u[j];
    bi = opt->b[i];
    cj = opt->c[j];

    opt->rflg[i] = 'F';
    opt->cflg[j] = 'F';
    opt->prtst = true;
    opt->nrdb++;
    opt->ncfx++;

    AddStack();

    Int4(RDBL, i, j, q);
    Dbl3(bi, aij, aik);
    Dbl2(cj, opt->c[q]);

    if (fabs(cj) >= par->aijtol)
    {
      lk = cj / aij;
      opt->c[q] -= lk * aik;
      opt->obj += bi * lk;
    }

    if (aij * aik > 0)
    {
      if (lj < -par->bplus)
        uk = MaxPositive;
      else
        uk = (bi - aij * lj) / aik;
      if (uj > par->bplus)
        lk = -MaxPositive;
      else
        lk = (bi - aij * uj) / aik;
    }

    else
    {
      if (lj < -par->bplus)
        lk = -MaxPositive;
      else
        lk = (bi - aij * lj) / aik;
      if (uj > par->bplus)
        uk = MaxPositive;
      else
        uk = (bi - aij * uj) / aik;
    }

    k = 0;
    if (lk > opt->l[q])
    {
      k = 1;
      opt->l[q] = lk;
    }
    Int1Dbl1(k, lk);

    k = 0;
    if (opt->u[q] > uk)
    {
      k = 1;
      opt->u[q] = uk;
    }
    Int1Dbl1(k, uk);

    lk = opt->l[q];
    k = 0;
    if (fabs(lk) >= par->aijtol)
      k = 1;

    aq = opt->at->ia + q;
    RecVect(aq, false);

    for (p = 0; p < aq->nn0; p++)
    {
      r = aq->ja[p];
      if (r == i || opt->rflg[r] != 'T')
        continue;

      apk = aq->an[p];
      if (k)
        opt->b[r] -= apk * lk;
    }

    if (k)
    {
      opt->obj += opt->c[q] * lk;
      if (opt->u[q] < par->bplus)
        opt->u[q] -= lk;
      opt->l[q] = 0.0;
      AddStack();
      Int2Dbl1(LBND, q, lk);
    }
  }

  return true;
} /* ChkRowDblt */

int SetColBnds(int m,
               int front,
               optpar *par,
               optmod *opt,
               int *ix)
{
  int i, j, k, p;
  double aij, uk, lk;
  array *aj;

  for (i = 0; i < m; i++)
  {
    opt->rlb[i] = -MaxPositive;
    opt->rub[i] = MaxPositive;
  }

  for (j = front; j >= 0; j = ix[j])
  {

    if (opt->cflg[j] != 'T')
      continue;

    aj = opt->at->ia + j;
    for (k = 0; k < aj->nn0; k++)
    {
      i = aj->ja[k];
      if (opt->rflg[i] == 'T')
        break;
    }

    if (k >= aj->nn0)
      continue;

    p = k;

    if (opt->u[j] > par->bplus)
    {

      aij = aj->an[p];

      if (aij > 0)
      {
        uk = opt->c[j] / aij;
        if (opt->rub[i] > uk)
          opt->rub[i] = uk;
      }

      else
      {
        lk = opt->c[j] / aij;
        if (opt->rlb[i] < lk)
          opt->rlb[i] = lk;
      }
    }
  }
  return true;
} /* SetCBnds */

int ChkColDomn(int n,
               int rear,
               optpar *par,
               optmod *opt,
               int *ix)
{
  int i, j, k, ilow, iupp;
  double dj, ej, cdj, cej, aij,
      cj, li, ui;
  array *aj;

  for (j = 0; j < n; j++)
  {

    if (opt->cflg[j] != 'T')
      continue;

    dj = 0.0;
    ej = 0.0;
    cj = opt->c[j];
    ilow = 0;
    iupp = 0;

    aj = opt->at->ia + j;

    for (k = 0; k < aj->nn0; k++)
    {
      i = aj->ja[k];
      if (opt->rflg[i] != 'T')
        continue;

      aij = aj->an[k];
      li = opt->rlb[i];
      ui = opt->rub[i];

      if (aij > 0)
      {
        if (li < -par->bplus)
          ilow++;
        else
          ej += aij * li;
        if (ui > par->bplus)
          iupp++;
        else
          dj += aij * ui;
      }

      else
      {
        if (li < -par->bplus)
          iupp++;
        else
          dj += aij * li;
        if (ui > par->bplus)
          ilow++;
        else
          ej += aij * ui;
      }
    }

    if (ilow)
      cej = MaxPositive;
    else
      cej = cj - ej;
    if (iupp)
      cdj = -MaxPositive;
    else
      cdj = cj - dj;

    if (cej > 0 && cdj < 0)
    {

      if (opt->u[j] > par->bplus && ilow < 2)
      {

        for (k = 0; k < aj->nn0; k++)
        {
          i = aj->ja[k];
          if (opt->rflg[i] != 'T')
            continue;

          aij = aj->an[k];
          if (fabs(aij) < par->aijtol)
            return ErrorProc(OVER_FLOW, NULL);

          if (aij > 0)
          {
            ui = opt->rlb[i];
            if (ui < -par->bplus)
              ui = (cj - ej) / aij;
            else if (ilow)
              ui = MaxPositive;
            else
              ui += cej / aij;
            opt->rub[i] = min(opt->rub[i], ui);
          }

          else
          {
            li = opt->rub[i];
            if (li > par->bplus)
              li = (cj - ej) / aij;
            else if (ilow)
              li = -MaxPositive;
            else
              li += cej / aij;
            opt->rlb[i] = max(opt->rlb[i], li);
          }
        }
      }
      continue;
    }

    else if (fabs(cej) < par->aijtol)
    {

      if (ix[j] >= 0 || j == rear)
        continue;

      ui = opt->u[j];
      if (ui < par->bplus)
      {
        opt->cflg[j] = 'F';
        opt->prtst = true;
        opt->ncdm++;

        AddStack();
        Int2Dbl2(CDMN, j, ui, opt->c[j]);
        RecVect(aj, false);

        for (k = 0; k < aj->nn0; k++)
        {

          i = aj->ja[k];
          if (opt->rflg[i] != 'T')
            continue;

          aij = aj->an[k];
          opt->b[i] -= aij * ui;
        }
        opt->obj += cj * ui;
      }
    }

    else if (fabs(cdj) < par->aijtol)
    {

      if (ix[j] >= 0 || j == rear)
        continue;

      if (opt->l[j] > -par->bplus)
      {
        opt->cflg[j] = 'F';
        opt->prtst = true;
        opt->ncdm++;
        AddStack();
        Int2Dbl2(CDMN, j, opt->l[j], opt->c[j]);
        RecVect(aj, false);
        for (k = 0; k < aj->nn0; k++)
        {
          i = aj->ja[k];
          if (opt->rflg[i] == 'T')
            opt->b[i] -= aij * opt->l[j];
        }
        opt->obj += opt->c[j] * opt->l[j];
      }
    }

    else if (cdj > 0)
    {
      if (opt->l[j] < -pr->aijtol)
      {
        ErrorProc(UNB_PROB, NULL);
        return false;
      }

      opt->cflg[j] = 'F';
      opt->prtst = true;
      opt->ncdm++;

      AddStack();
      Int2Dbl2(CDMN, j, opt->l[j], opt->c[j]);
      RecVect(aj, false);

      for (k = 0; k < aj->nn0; k++)
      {
        i = aj->ja[k];
        if (opt->rflg[i] == 'T')
          opt->b[i] -= aij * opt->l[j];
      }
      opt->obj += opt->c[j] * opt->l[j];
    }

    else
    {

      ui = opt->u[j];
      if (ui > par->bplus)
      {
        ErrorProc(INF_PROB, NULL);
        return false;
      }

      opt->cflg[j] = 'F';
      opt->prtst = true;
      opt->ncdm++;

      AddStack();
      Int2Dbl2(CDMN, j, ui, opt->c[j]);
      RecVect(aj, false);

      for (k = 0; k < aj->nn0; k++)
      {

        i = aj->ja[k];
        if (opt->rflg[i] != 'T')
          continue;

        aij = aj->an[k];
        opt->b[i] -= aij * ui;
      }
      opt->obj += cj * ui;
    }
  }
  return true;
} /* ChkCDomn */

int ChkRowDupl(int m,
               int n,
               int rear,
               optpar *par,
               optmod *opt,
               int *iw,
               int *ix)
{
  int i, j, k, o, p, q, f, g, h, x, y, z, *relm,
      *idsg, iii, *lfnt, *lrer, ppp, ik;
  double v, cj, *psg, xx, bk, *isg, aij, akl,
      lj, uj, lg, ug, cg;
  array *ai, *aj, *ap;

  psg = opt->clb;
  isg = opt->cub;
  relm = iw + n;
  lfnt = ix + n;
  idsg = relm + m;
  lrer = lfnt + m;

  for (i = 0; i < m; i++)
  {

    lfnt[i] = -1;
    lrer[i] = -1;

    if (opt->rflg[i] != 'T')
      continue;

    o = 0;
    p = 0;
    q = 0;

    ai = opt->a->ia + i;
    for (k = 0; k < ai->nn0; k++)
    {
      j = ai->ja[k];
      if (opt->cflg[j] != 'T')
        continue;

      if (ix[j] >= 0 || j == rear)
      {
        q++;
        continue;
      }

      if (!o)
      {
        aij = ai->an[k];
        if (fabs(aij) >= par->aijtol)
        {
          opt->rlb[i] = ai->an[k];
          o = 1;
        }
      }

      p++;
    }
    relm[i] = p;
    opt->rub[i] = q;
  }

  for (j = 0; j < n; j++)
  {

    iw[j] = -1;
    if (opt->cflg[j] != 'T')
      continue;

    o = 0;
    aj = opt->at->ia + j;

    for (k = 0; k < aj->nn0; k++)
    {
      i = aj->ja[k];
      if (opt->rflg[i] != 'T')
        continue;
      o++;
    }

    if (lfnt[o] < 0)
    {
      lfnt[o] = j;
      lrer[o] = j;
    }

    else
    {
      iw[lrer[o]] = j;
      lrer[o] = j;
    }
  }

  for (h = 2; h <= m; h++)
  {

    for (j = lfnt[h]; j >= 0; j = iw[j])
    {

      aj = opt->at->ia + j;

      for (k = 0; k < aj->nn0; k++)
      {

        i = aj->ja[k];
        if (opt->rflg[i] != 'T')
          continue;

        o = relm[i];
        v = opt->rlb[i];

        ai = opt->a->ia + i;

        opt->rflg[i] = '+';

        for (q = k + 1; q < aj->nn0; q++)
        {

          p = aj->ja[q];
          if (opt->rflg[p] != 'T')
            continue;

          if (o != relm[p])
            continue;

          ap = opt->a->ia + p;

          iii = 0;
          ppp = 0;
          v = opt->rlb[p] / v;

          x = 0;
          y = 0;

          while (x < ai->nn0 || y < ap->nn0)
          {
            if (x < ai->nn0)
            {
              f = ai->ja[x];
              if (opt->cflg[f] != 'T')
              {
                x++;
                continue;
              }
              else if (ix[f] >= 0 || f == rear)
              {
                idsg[iii] = f;
                isg[iii] = ai->an[x];
                x++;
                iii++;
                continue;
              }
            }

            if (y < ap->nn0)
            {
              g = ap->ja[y];
              if (opt->cflg[g] != 'T')
              {
                y++;
                continue;
              }
              else if (ix[g] >= 0 || g == rear)
              {
                lrer[ppp] = g;
                psg[ppp] = ap->an[y];
                y++;
                ppp++;
                continue;
              }
            }
            if (f != g || fabs(v * ai->an[x] - ap->an[y]) >= par->aijtol)
              break;
            x++;
            y++;
          }

          if (x < ai->nn0 || y < ap->nn0)
            continue;

          if (ppp)
          {

            z = iii + ppp;

            if (z < 2)
            {

              opt->prtst = true;
              opt->rflg[p] = 'F';
              opt->nrdu++;

              AddStack();
              Int3Dbl1(RDUP, i, p, v);

              g = lrer[0];
              bk = opt->b[p] - v * opt->b[i];
              bk /= psg[0];
              cg = opt->c[g];

              if (bk < opt->l[g] || bk > opt->u[g])
              {
                ErrorProc(INF_PROB, NULL);
                return false;
              }

              opt->cflg[g] = 'F';
              opt->obj += cg * bk;
              opt->ncfx++;

              AddStack();
              Int3Dbl3(RSNG, p, g, bk, psg[0], cg);
            } /*z<2*/

            else if (z == 2)
            {

              opt->prtst = true;
              opt->rflg[p] = 'F';
              opt->nrdu++;

              if (iii)
              {

                f = idsg[0];
                g = lrer[0];
                bk = opt->b[p] - v * opt->b[i];
                aij = v * isg[0];
                akl = psg[0];
                lg = opt->l[g];
                ug = opt->u[g];

                if (aij * akl > 0)
                {
                  if (ug > par->bplus)
                    uj = MaxPositive;
                  else
                    uj = (akl * ug - bk) / aij;
                  lj = (akl * lg - bk) / aij;
                }

                else
                {
                  if (uj > par->bplus)
                    lj = -MaxPositive;
                  else
                    lj = (akl * ug - bk) / aij;
                  uj = (akl * lg - bk) / aij;
                }

                lj = max(lj, opt->l[f]);
                uj = min(uj, opt->u[f]);

                if (lj > uj)
                {
                  ErrorProc(INF_PROB, NULL);
                  return false;
                }

                opt->cflg[g] = 'F';
                opt->ncfx++;

                cj = opt->c[f];
                cg = opt->c[g];

                if (fabs(cg) >= par->aijtol)
                {
                  xx = cg / akl;
                  opt->c[f] += xx * aij;
                  opt->obj += xx * bk;
                }
                aij = -aij;

                AddStack();
                Int3Dbl1(RDUP, i, p, v);

                AddStack();
                Int4(RDBL, p, g, f);
                Dbl3(bk, akl, aij);
                Dbl2(cg, cj);

                ik = 0;
                if (lj > opt->l[f])
                  ik = 1;
                Int1Dbl1(ik, lj);

                ik = 0;
                if (uj < opt->u[f])
                  ik = 1;
                Int1Dbl1(ik, uj);
                Int1Dbl1(i, isg[0]);
                ik = -1;
                Int1(ik);

                if (lj > opt->l[f])
                {

                  opt->obj += opt->c[f] * lj;
                  opt->b[i] -= isg[0] * lj;

                  if (uj <= par->bplus)
                    uj -= lj;

                  AddStack();
                  Int2Dbl1(LBND, f, lj);

                  opt->l[f] = 0.0;
                }

                opt->u[f] = uj;
              }

              else
              {

                f = lrer[0];
                g = lrer[1];
                bk = opt->b[p] - v * opt->b[i];
                aij = psg[0];
                akl = psg[1];
                lj = opt->l[f];
                uj = opt->u[f];

                if (aij * akl > 0)
                {
                  if (uj > par->bplus)
                    lg = -MaxPositive;
                  else
                    lg = (bk - aij * uj) / akl;
                  if (lj < -par->bplus)
                    ug = MaxPositive;
                  else
                    ug = (bk - aij * lj) / akl;
                }

                else
                {
                  if (lj < -par->bplus)
                    lg = -MaxPositive;
                  else
                    lg = (bk - aij * lj) / akl;
                  if (uj > par->bplus)
                    ug = MaxPositive;
                  else
                    ug = (bk - aij * uj) / akl;
                }

                ug = min(ug, opt->u[g]);
                lg = max(lg, opt->l[g]);

                if (ug < lg)
                {
                  ErrorProc(INF_PROB, NULL);
                  return false;
                }

                opt->cflg[f] = 'F';
                cj = opt->c[f];
                cg = opt->c[g];

                if (fabs(cj) >= par->aijtol)
                {
                  xx = cj / aij;
                  opt->c[g] -= xx * akl;
                  opt->obj += xx * bk;
                }

                AddStack();
                Int3Dbl1(RDUP, i, p, v);

                AddStack();
                Int4(RDBL, p, f, g);
                Dbl3(bk, aij, akl);
                Dbl2(cj, cg);

                ik = 0;
                if (lg > opt->l[g])
                  ik = 1;
                Int1Dbl1(ik, lg);

                ik = 0;
                if (ug < opt->u[g])
                  ik = 1;
                Int1Dbl1(ik, ug);
                ik = -1;
                Int1(ik);

                if (lg > opt->l[g])
                {
                  opt->obj += cg * lg;
                  if (ug < par->bplus)
                    ug -= lg;
                  AddStack();
                  Int2Dbl1(LBND, g, lg);
                  opt->l[g] = 0.0;
                }

                opt->u[g] = ug;
              }
            } /*z=2*/
          }   /*if(ppp)*/

          else
          {
            if (!iii)
            {
              bk = opt->b[p] - v * opt->b[i];
              if (fabs(bk) >= par->aijtol)
              {
                ErrorProc(INF_PROB, NULL);
                return false;
              }

              AddStack();
              Int3Dbl1(RDUP, i, p, v);

              AddStack();
              Int2(RNUL, p);

              opt->rflg[p] = 'F';
              opt->prtst = true;
              opt->nrdu++;
            }

            else if (iii == 1)
            {

              f = idsg[0];
              bk = opt->b[i] - opt->b[p] / v;
              bk /= isg[0];
              if (bk < opt->l[f] || bk > opt->u[f])
              {
                ErrorProc(INF_PROB, NULL);
                return false;
              }

              opt->cflg[f] = 'F';
              opt->rflg[i] = 'F';
              opt->prtst = true;
              opt->nrdu++;
              opt->ncfx++;

              cj = opt->c[f];
              opt->obj += cj * bk;
              v = 1.0 / v;

              AddStack();
              Int3Dbl1(RDUP, p, i, v);

              AddStack();
              Int3Dbl3(RSNG, i, f, bk, isg[0], cj);
              break;
            }

            else if (iii == 2)
            {

              f = idsg[0];
              g = idsg[1];
              bk = opt->b[i] - opt->b[p] / v;
              aij = isg[0];
              akl = isg[1];
              lj = opt->l[f];
              uj = opt->u[f];
              cj = opt->c[f];
              cg = opt->c[g];

              if (aij * akl > 0)
              {
                if (uj > par->bplus)
                  lg = -MaxPositive;
                else
                  lg = (bk - aij * uj) / akl;
                if (lj < -par->bplus)
                  ug = MaxPositive;
                else
                  ug = (bk - aij * lj) / akl;
              }

              else
              {
                if (lj < -par->bplus)
                  lg = -MaxPositive;
                else
                  lg = (bk - aij * lj) / akl;
                if (uj > par->bplus)
                  ug = MaxPositive;
                else
                  ug = (bk - aij * uj) / akl;
              }

              ug = min(ug, opt->u[g]);
              lg = max(lg, opt->l[g]);

              if (ug < lg)
              {
                ErrorProc(INF_PROB, NULL);
                return false;
              }

              opt->cflg[f] = 'F';
              opt->ncfx++;
              if (fabs(cj) >= par->aijtol)
              {
                xx = cj / aij;
                opt->c[g] -= xx * akl;
                opt->obj += xx * bk;
              }
              v = 1.0 / v;

              AddStack();
              Int3Dbl1(RDUP, p, i, v);

              AddStack();
              Int4(RDBL, i, f, g);
              Dbl3(bk, aij, akl);
              Dbl2(cj, cg);

              ik = 0;
              if (lg > opt->l[g])
                ik = 1;
              Int1Dbl1(ik, lg);

              ik = 0;
              if (ug < opt->u[g])
                ik = 1;
              Int1Dbl1(ik, ug);
              ik = -1;
              Int1(ik);

              if (lg > opt->l[g])
              {
                opt->obj += opt->c[g] * lg;
                if (ug < par->bplus)
                  ug -= lg;

                AddStack();
                Int2Dbl1(LBND, g, lg);
                opt->l[g] = 0.0;
              }

              opt->u[g] = ug;
              opt->rflg[i] = 'F';
              opt->prtst = true;
              opt->nrdu++;
              break;
            }
          }
        }
      }
    }
  }

  for (i = 0; i < m; i++)
    if (opt->rflg[i] == '+')
      opt->rflg[i] = 'T';

  return true;
} /* ChkRowDupl */

int ChkColDupl(int m,
               int n,
               optpar *par,
               optmod *opt,
               int *iw,
               int *ix)
{
  int i, j, k, o, p, q, f, g, h, x, y, z,
      *ft, *re;
  double cj, ck, dj, ej, uj, lk, uk, v,
      aij, lj;
  array *ai, *aj, *ak, *af;

  ft = iw + n;
  re = ix + m;

  for (j = 0; j < n; j++)
  {

    ft[j] = -1;
    re[j] = -1;

    if (opt->cflg[j] != 'T')
      continue;
    o = 0;
    p = 0;

    aj = opt->at->ia + j;

    for (k = 0; k < aj->nn0; k++)
    {

      i = aj->ja[k];
      if (opt->rflg[i] != 'T')
        continue;

      if (!o)
      {
        aij = aj->an[k];
        if (fabs(aij) >= par->aijtol)
        {
          opt->rlb[j] = aij;
          o = 1;
        }
      }
      p++;
    }
    iw[j] = p;
  }

  ft[n] = -1;
  re[n] = -1;

  for (i = 0; i < m; i++)
  {

    if (opt->rflg[i] != 'T')
      continue;

    o = 0;
    ix[i] = -1;
    ai = opt->a->ia + i;

    for (k = 0; k < ai->nn0; k++)
    {
      j = ai->ja[k];
      if (opt->cflg[j] != 'T')
        continue;
      o++;
    }

    if (ft[o] < 0)
    {
      ft[o] = i;
      re[o] = i;
    }

    else
    {
      ix[re[o]] = i;
      re[o] = i;
    }
  }

  for (h = 2; h <= n; h++)
  {

    for (i = ft[h]; i >= 0; i = ix[i])
    {

      ai = opt->a->ia + i;

      for (p = 0; p < ai->nn0; p++)
      {

        j = ai->ja[p];
        if (opt->cflg[j] != 'T')
          continue;

        opt->cflg[j] = '+';
        o = iw[j];
        v = opt->rlb[j];
        cj = opt->c[j];
        uj = opt->u[j];
        lj = opt->l[j];

        aj = opt->at->ia + j;

        for (q = p + 1; q < ai->nn0; q++)
        {

          k = ai->ja[q];
          if (opt->cflg[k] != 'T')
            continue;
          if (o != iw[k])
            continue;

          ak = opt->at->ia + k;

          v = opt->rlb[k] / v;
          x = 0;
          y = 0;

          while (x < aj->nn0 && y < ak->nn0)
          {

            f = aj->ja[x];
            if (opt->rflg[f] != 'T')
            {
              x++;
              continue;
            }

            g = ak->ja[y];
            if (opt->rflg[g] != 'T')
            {
              y++;
              continue;
            }

            if (f != g || fabs(ak->an[y] - v * aj->an[x]) >= par->aijtol)
              break;

            x++;
            y++;
          }

          if (x < aj->nn0 || y < ak->nn0)
            continue;

          ck = opt->c[k] - v * cj;
          uk = opt->u[k];
          lk = opt->l[k];

          if (fabs(ck) < par->aijtol)
          {

            if (v > 0.0)
            {
              if (lj < -par->bplus || lk < -par->bplus)
                dj = -MaxPositive;
              else
                dj = lj + v * lk;
              if (uj > par->bplus || uk > par->bplus)
                ej = MaxPositive;
              else
                ej = uj + v * uk;
            }

            else
            {
              if (uk > par->bplus || lj < -par->bplus)
                dj = -MaxPositive;
              else
                dj = lj + v * uk;
              if (uj > par->bplus || lk < -par->bplus)
                ej = MaxPositive;
              else
                ej = uj + v * lk;
            }

            if (ej > par->bplus && dj < -par->bplus)
            {

              if (o != 1)
                continue;
              aij = opt->rlb[j];

              AddStack();
              Int3Dbl3(CDUP, k, j, v, lk, uk);
              Dbl2(lj, uj);

              AddStack();
              Int3Dbl3(RDMN, f, j, opt->b[f], aij, cj);

              cj /= aij;

              af = opt->a->ia + f;
              RecVect(af, true);

              for (g = 0; g < af->nn0; g++)
              {
                z = af->ja[g];
                if (opt->cflg[z] == 'F')
                  continue;
                if (z == j || z == k)
                  continue;
                opt->c[z] -= cj * af->an[g];
              }

              opt->obj += opt->b[f] * cj;
              opt->prtst = true;
              opt->cflg[j] = 'F';
              opt->cflg[k] = 'F';
              opt->rflg[f] = 'F';
              opt->ncdu += 2;
              opt->nrdm++;
              break;
            }

            else
            {
              opt->prtst = true;
              opt->cflg[k] = 'F';
              opt->ncdu++;

              AddStack();
              Int3Dbl3(CDUP, k, j, v, lk, uk);
              Dbl2(lj, uj);

              if (dj < -par->bplus)
              {

                for (g = 0; g < aj->nn0; g++)
                {
                  x = aj->ja[g];
                  if (opt->rflg[x] == 'F')
                    continue;

                  aij = aj->an[g];
                  aj->an[g] = -aij;
                  opt->b[x] -= aij * ej;

                  af = opt->a->ia + x;
                  for (z = 0; z < af->nn0; z++)
                    if (af->ja[z] == j)
                      af->an[z] = -af->an[z];
                }

                AddStack();
                Int2Dbl1(MINS, j, ej);

                opt->obj += cj * ej;
                opt->c[j] = -cj;
                opt->u[j] = MaxPositive;
                opt->l[j] = 0.0;
              }

              else
              {
                if (fabs(dj) >= par->aijtol)
                {
                  if (ej <= par->bplus)
                    ej -= dj;
                  opt->obj += cj * dj;
                  for (g = 0; g < aj->nn0; g++)
                  {
                    x = aj->ja[g];
                    if (opt->rflg[x] != 'T')
                      continue;
                    opt->b[x] -= aj->an[g] * dj;
                  }
                  AddStack();
                  Int2Dbl1(LBND, j, dj);
                  opt->l[j] = 0.0;
                }
                opt->u[j] = ej;
              }
            } /*else*/
          }   /*if ck==0*/

          else if (ck > 0.0)
          {
            if (v > 0.0 && uj > par->bplus)
            {

              if (lk < -par->bplus)
              {
                ErrorProc(UNB_PROB, NULL);
                return false;
              }

              opt->prtst = true;
              opt->cflg[k] = 'F';
              opt->ncdu++;

              AddStack();
              Int2Dbl2(CFIX, k, lk, opt->c[k]);
              RecVect(ak, false);

              for (g = 0; g < ak->nn0; g++)
              {
                z = ak->ja[g];
                if (opt->rflg[z] != 'T')
                  continue;
                opt->b[z] -= ak->an[g] * lk;
              }
              opt->obj += opt->c[k] * lk;

              if (lj < -par->bplus)
                opt->l[j] = -MaxPositive;

              else
              {
                dj = lj + v * lk;
                if (fabs(dj) >= par->aijtol)
                {
                  opt->obj += opt->c[j] * dj;
                  for (g = 0; g < aj->nn0; g++)
                  {
                    x = aj->ja[g];
                    if (opt->rflg[x] == 'F')
                      continue;
                    opt->b[x] -= aj->an[g] * dj;
                  }
                }
                AddStack();
                Int2Dbl1(LBND, j, dj);
                opt->l[j] = 0.0;
              }
              opt->u[j] = MaxPositive;
            }
          } /*ck>0*/

          else
          {

            if (v < 0.0 && uj > par->bplus)
            {

              if (uk > par->bplus)
              {
                ErrorProc(UNB_PROB, NULL);
                return false;
              }

              opt->prtst = true;
              opt->cflg[k] = 'F';
              opt->ncdu++;

              AddStack();
              Int2Dbl2(CFIX, k, uk, opt->c[k]);
              RecVect(ak, false);

              for (g = 0; g < ak->nn0; g++)
              {
                z = ak->ja[g];
                if (opt->rflg[z] == 'F')
                  continue;
                opt->b[z] -= ak->an[g] * uk;
              }

              opt->obj += opt->c[k] * uk;

              if (lj < -par->bplus)
                opt->l[j] = -MaxPositive;

              else
              {
                dj = lj + v * uk;
                if (fabs(dj) >= par->aijtol)
                {
                  opt->obj += cj * dj;
                  for (g = 0; g < aj->nn0; g++)
                  {
                    x = aj->ja[g];
                    if (opt->rflg[x] == 'F')
                      continue;
                    opt->b[x] -= aj->an[g] * dj;
                  }

                  AddStack();
                  Int2Dbl1(LBND, j, dj);
                }
                opt->l[j] = 0.0;
              }

              opt->u[j] = MaxPositive;
            }

            else if (v < 0.0 && uk > par->bplus)
            {

              if (uj > par->bplus)
              {
                ErrorProc(UNB_PROB, NULL);
                return false;
              }

              opt->prtst = true;
              opt->cflg[j] = 'F';
              opt->ncdu++;

              AddStack();
              Int2Dbl2(CFIX, j, uj, cj);
              RecVect(aj, false);

              for (g = 0; g < aj->nn0; g++)
              {
                z = aj->ja[g];
                if (opt->rflg[z] != 'T')
                  continue;
                opt->b[z] -= aj->an[g] * uj;
              }

              opt->obj += cj * uj;

              if (lk < -par->bplus)
                opt->l[k] = -MaxPositive;

              else
              {
                dj = lk + uj / v;
                if (fabs(dj) >= par->aijtol)
                {
                  AddStack();
                  Int2Dbl1(LBND, k, dj);

                  opt->obj += opt->c[k] * dj;
                  for (g = 0; g < ak->nn0; g++)
                  {
                    x = ak->ja[g];
                    if (opt->rflg[x] == 'F')
                      continue;
                    opt->b[x] -= ak->an[g] * dj;
                  }
                }
                opt->l[k] = 0.0;
              }
              opt->u[k] = MaxPositive;
              break;
            } /*6: v<0&&uk=+inf*/

            else if (v > 0.0 && uk > par->bplus)
            {

              if (lj < -par->bplus)
              {
                ErrorProc(UNB_PROB, NULL);
                return false;
              }

              opt->prtst = true;
              opt->cflg[j] = 'F';
              opt->ncdu++;

              AddStack();
              Int2Dbl2(CFIX, j, lj, cj);
              RecVect(aj, false);

              for (g = 0; g < aj->nn0; g++)
              {
                z = aj->ja[g];
                if (opt->rflg[z] != 'T')
                  continue;
                opt->b[z] -= aj->an[g] * lj;
              }
              opt->obj += opt->c[j] * lj;

              if (lk < -par->bplus)
                opt->l[k] = -MaxPositive;

              else
              {
                dj = lk + lj / v;
                if (fabs(dj) >= par->aijtol)
                {
                  opt->obj += opt->c[k] * dj;
                  for (g = 0; g < ak->nn0; g++)
                  {
                    x = ak->ja[g];
                    if (opt->rflg[x] == 'F')
                      continue;
                    opt->b[x] -= ak->an[g] * dj;
                  }
                  AddStack();
                  Int2Dbl1(LBND, k, dj);
                }
                opt->l[k] = 0.0;
              }
              opt->u[k] = MaxPositive;
              break;
            } /* 6: else v>0&&uk=+inf*/
          }   /*5: else ck<0*/
        }     /*4: for q*/
      }       /*3: for p*/
    }         /*2: for i*/
  }           /*1: for h*/

  for (j = 0; j < n; j++)
    if (opt->cflg[j] == '+')
      opt->cflg[j] = 'T';

  return true;
} /* ChkCDupl */
