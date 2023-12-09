#include "LPreport.h"

char cGet(FILE *fp)
{
  char cs;

  fread(&cs, 1, 1, fbin);
  return cs;
} /* cGet */

int iGet(FILE *fp)
{
  int i;

  fread(&i, 4, 1, fp);
  return (i);
} /* iGet */

double dGet(FILE *fp)
{
  double r;

  fread(&r, 8, 1, fp);
  return (r);
} /* dGet */

void PostSol(int key,
             solpt *pt)
{
  int i, j, k, m, n, nn0, il, iu, itmp,
      *subi, *subj;
  double *x, *y, *z, bi, cj, ck, lj, uj, lk, uk,
      aij, aik, ymin, ymax, rtmp;
  char ckey;
  array *aj;

  switch (key)
  {
  case CDMN:
    pt->cdmn++;
    j = iGet(fbin);
    pt->x0[j] = dGet(fbin);
    cj = dGet(fbin);
    i = iGet(fbin);
    while (i >= 0)
    {
      rtmp = dGet(fbin);
      cj -= pt->y0[i] * rtmp;
      i = iGet(fbin);
    }
    pt->z0[j] = cj;
    break;

  case CDUP:
    pt->cdup++;
    k = iGet(fbin);
    j = iGet(fbin);
    rtmp = dGet(fbin);
    lk = dGet(fbin);
    uk = dGet(fbin);
    lj = dGet(fbin);
    uj = dGet(fbin);
    cj = pt->x0[j];
    if (fabs(rtmp) < 1e-15)
    {
      printf(" err data for duplicate column.\n");
      exit(0);
    }
    if (rtmp > 0.0)
    {
      if (uj <= 1e29)
        ymin = (cj - uj) / rtmp;
      else
        ymin = -1e30;
      if (lk >= -1e29)
        ymax = (cj - lj) / rtmp;
      else
        ymax = 1e30;
      ymin = max(ymin, lk);
      ymax = min(ymax, uk);
    }
    else
    {
      if (lj >= -1e29)
        ymin = (cj - lj) / rtmp;
      else
        ymin = -1e30;
      if (uk <= 1e29)
        ymax = (cj - uj) / rtmp;
      else
        ymax = 1e30;
      ymin = max(ymin, lk);
      ymax = min(ymax, uk);
    }
    if (ymin < -1e29 && ymax > 1e29)
      cj = 0.0;
    else if (ymax > 1e29)
      cj = ymin;
    else if (ymin < -1e29)
      cj = ymax;
    else
      cj = 0.5 * (ymax + ymin);

    pt->x0[j] = pt->x0[j] - rtmp * cj;
    pt->x0[k] = cj;
    pt->z0[k] = rtmp * pt->z0[j];
    break;

  case CFIX:
    pt->cfix++;
    j = iGet(fbin);
    pt->x0[j] = dGet(fbin);
    pt->z0[j] = dGet(fbin);

    i = iGet(fbin);
    while (i >= 0)
    {
      rtmp = dGet(fbin);
      pt->z0[j] -= pt->y0[i] * rtmp;
      i = iGet(fbin);
    }
    break;

  case CNUL:
    pt->cnul++;
    j = iGet(fbin);
    pt->x0[j] = dGet(fbin);
    pt->z0[j] = dGet(fbin);
    break;

  case LBND:
    pt->lbnd++;
    j = iGet(fbin);
    lj = dGet(fbin);
    pt->x0[j] += lj;
    break;

  case MATX:
    pt->non0 = iGet(fbin);
    pt->n0 = iGet(fbin);

    pt->c = dAlloc(pt->n0);
    pt->a = mAlloc(pt->n0, pt->non0);

    fread(pt->c, 8, pt->n0, fbin);

    nn0 = 0;
    for (j = 0; j < pt->n0; j++)
    {
      aj = pt->a->cols + j;

      k = iGet(fbin);
      aj->nn0 = k;
      aj->ja = pt->a->cols->ja + nn0;
      aj->an = pt->a->cols->an + nn0;

      fread(aj->ja, 4, k, fbin);
      fread(aj->an, 8, k, fbin);

      nn0 += k;
    }
    break;

  case MINS:
    pt->mins++;
    j = iGet(fbin);
    rtmp = dGet(fbin);
    pt->x0[j] = rtmp - pt->x0[j];
    pt->z0[j] = -pt->z0[j];
    break;

  case NMES:
    pt->m0 = iGet(fbin);
    pt->n0 = iGet(fbin);

    k = iGet(fbin);
    fread(pt->sname, 1, k, fbin);
    pt->sname[k] = '\0';

    pt->r = dAlloc(pt->m0);
    pt->b = dAlloc(pt->m0);
    pt->rp = cAlloc(pt->m0);
    pt->l = dAlloc(pt->n0);
    pt->u = dAlloc(pt->n0);

    fread(pt->r, 8, pt->m0, fbin);
    fread(pt->b, 8, pt->m0, fbin);
    fread(pt->rp, 1, pt->m0, fbin);
    fread(pt->l, 8, pt->n0, fbin);
    fread(pt->u, 8, pt->n0, fbin);

    break;
  case RDBL:
    pt->rdbl++;
    i = iGet(fbin);
    j = iGet(fbin);
    k = iGet(fbin);
    bi = dGet(fbin);
    aij = dGet(fbin);
    aik = dGet(fbin);
    cj = dGet(fbin);
    ck = dGet(fbin);
    il = iGet(fbin);
    lk = dGet(fbin);
    iu = iGet(fbin);
    uk = dGet(fbin);

    pt->x0[j] = (bi - aik * pt->x0[k]) / aij;
    if ((il && pt->x0[k] == lk) || (iu && pt->x0[k] == uk))
    {
      pt->z0[k] = 0.0;
      itmp = iGet(fbin);
      while (itmp >= 0)
      {
        rtmp = dGet(fbin);
        ck -= rtmp * pt->y0[itmp];
        itmp = iGet(fbin);
      }
      pt->y0[i] = ck / aik;
    }
    else
      pt->y0[i] = cj / aij;
    pt->z0[j] = cj - aij * pt->y0[i];
    break;

  case RDMN:
    pt->rdmn++;
    i = iGet(fbin);
    j = iGet(fbin);
    bi = dGet(fbin);
    aij = dGet(fbin);
    cj = dGet(fbin);
    k = iGet(fbin);
    while (k >= 0)
    {
      rtmp = dGet(fbin);
      if (k != j)
        bi -= rtmp * pt->x0[k];
      k = iGet(fbin);
    }
    pt->x0[j] = bi / aij;
    pt->y0[i] = cj / aij;
    pt->z0[j] = 0.0;
    break;

  case RDUP:
    pt->rdup++;
    i = iGet(fbin);
    k = iGet(fbin);
    rtmp = dGet(fbin);
    pt->y0[i] -= rtmp * pt->y0[k];
    break;

  case RFRC:
    pt->rfrc++;
    i = iGet(fbin);
    ymin = 1.0e40;
    ymax = -ymin;
    itmp = 0;

    if (!pt->subi)
      pt->subi = iAlloc(pt->n);
    if (!pt->ai)
      pt->ai = dAlloc(pt->n);

    ckey = cGet(fbin);
    while (ckey != 'e')
    {
      j = iGet(fbin);
      lj = dGet(fbin);
      aij = dGet(fbin);
      cj = dGet(fbin);
      k = iGet(fbin);
      while (k >= 0)
      {
        rtmp = dGet(fbin);
        cj -= rtmp * pt->y0[k];
        k = iGet(fbin);
      }
      if (fabs(aij) < 1.0e-13)
        exit(0);
      pt->x0[j] = lj;
      pt->z0[j] = cj;
      pt->subi[itmp] = j;
      pt->ai[itmp] = aij;
      itmp++;

      rtmp = cj / aij;
      if (aij > 0.0)
      {
        if (ckey == 'l')
        {
          if (ymin > rtmp)
            ymin = rtmp;
        }
        else
        {
          if (ymax < rtmp)
            ymax = rtmp;
        }
      }
      else
      {
        if (ckey == 'l')
        {
          if (ymax < rtmp)
            ymax = rtmp;
        }
        else
        {
          if (ymin > rtmp)
            ymin = rtmp;
        }
      }

      ckey = cGet(fbin);
    }

    if (ymax < -1.0e30)
    {
      if (ymin > 1.0e30)
        pt->y0[i] = 0.0;
      else
        pt->y0[i] = ymin;
    }
    else if (ymin > 1.0e30)
      pt->y0[i] = ymax;
    else if (ymin < ymax)
      exit(0);
    else
      pt->y0[i] = 0.5 * (ymin + ymax);

    for (k = 0; k < itmp; k++)
    {
      j = pt->subi[k];
      pt->z0[j] -= pt->ai[k] * pt->y0[i];
    }
    break;

  case RNUL:
    pt->rnul++;
    i = iGet(fbin);
    pt->y0[i] = 0.0;
    break;

  case RSNG:
    pt->rsng++;
    i = iGet(fbin);
    j = iGet(fbin);
    pt->x0[j] = dGet(fbin);
    aij = dGet(fbin);
    pt->y0[i] = dGet(fbin);
    k = iGet(fbin);
    while (k >= 0)
    {
      rtmp = dGet(fbin);
      if (k != i)
        pt->y0[i] -= rtmp * pt->y0[k];
      k = iGet(fbin);
    }
    pt->y0[i] /= aij;
    pt->z0[j] = 0.0;
    break;

  case XSOL:
    n = iGet(fbin);
    subi = iAlloc(n);
    x = dAlloc(n);
    z = dAlloc(n);

    fread(subi, 4, n, fbin);
    for (k = 0; k < n; k++)
    {
      j = subi[k];
      if (j < 0)
      {
        x[k] = 0.0;
        z[k] = 0.0;
      }
      else
      {
        x[k] = pt->x[j];
        z[k] = pt->z[j];
      }
    }

    pt->x0 = x;
    pt->z0 = z;
    pt->subi = subi;
    pt->n = n;
    break;

  case YSOL:
    m = iGet(fbin);
    subj = iAlloc(m);
    y = dAlloc(m);

    fread(subj, 4, m, fbin);
    for (k = 0; k < m; k++)
    {
      i = subj[k];
      if (i < 0)
        y[k] = 0.0;
      else
        y[k] = pt->y[i];
    }

    pt->y0 = y;
    pt->subj = subj;
    pt->m = m;
    break;
  default:
    break;
  }
} /* PostSol */
