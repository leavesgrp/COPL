#include "LPcross.h"
//#define TIME_COUNT

void lu_delcur(gsdec *lf,
               int j)
{
  int iend = lf->wsze, cpre, cnex;

  if (lf->ccur[j] == lf->chead[j])
  {
    cnex = lf->cnext[lf->chead[j]];
    if (cnex == iend)
      lf->ctail[j] = iend;
    else
      lf->cprev[cnex] = iend;

    lf->chead[j] = cnex;
  }
  else
  {
    cpre = lf->cprev[lf->ccur[j]];
    cnex = lf->cnext[lf->ccur[j]];

    if (cpre == iend)
      lf->chead[j] = iend;
    else
      lf->cnext[cpre] = cnex;

    if (cnex == iend)
      lf->ctail[j] = cpre;
    else
      lf->cprev[cnex] = cpre;
  }

  lf->ccur[j] = cnex;
} /* lu_delcur */

void lu_delpos(gsdec *lf,
               int j,
               int pos)
{
  int temp;

  if (pos == lf->ccur[j])
    temp = lf->cnext[lf->ccur[j]];
  else
    temp = lf->ccur[j];

  lf->ccur[j] = pos;

  lu_delcur(lf, j);

  lf->ccur[j] = temp;
} /* lf_delpos */

void lu_inscur(gsdec *lf,
               int j,
               int pos)
{
  int cpre, iend = lf->wsze;

  if (lf->chead[j] == lf->wsze)
  {
    lf->cprev[pos] = iend;
    lf->cnext[pos] = iend;

    lf->chead[j] = pos;
    lf->ctail[j] = pos;
  }
  else if (lf->ccur[j] == lf->chead[j])
  {
    lf->cprev[pos] = iend;
    lf->cnext[pos] = lf->chead[j];

    lf->cprev[lf->chead[j]] = pos;

    lf->chead[j] = pos;
  }
  else if (lf->ccur[j] == lf->wsze)
  {

    lf->cprev[pos] = lf->ctail[j];
    lf->cnext[pos] = iend;

    lf->cnext[lf->ctail[j]] = pos;

    lf->ctail[j] = pos;
  }
  else
  {
    cpre = lf->cprev[lf->ccur[j]];

    lf->cprev[pos] = cpre;
    lf->cnext[pos] = lf->ccur[j];

    lf->cprev[lf->ccur[j]] = pos;
    lf->cnext[cpre] = pos;
  }

  lf->ccur[j] = pos;
} /* lu_inscur */

void lu_buildcol(gsdec *lf,
                 int incall,
                 double colmax[])
{
  int i, j, t, itemp, stop,
      iend = lf->wsze, nrow = lf->nrow, ncol = lf->ncol,
      *chead = lf->chead, *ctail = lf->ctail,
      *cprev = lf->cprev, *cnext = lf->cnext,
      *ccur = lf->ccur,
      *ufir = lf->ufir, *uuse = lf->uuse,
      *subj = lf->subj;
  double abs_aij, *val = lf->val;

  for (i = 0; i < ncol; i++)
  {
    chead[i] = iend;
    ctail[i] = iend;
    ccur[i] = iend;
  }

  if (colmax)
  {
    dZero(ncol, colmax, NULL);

    for (i = 0; i < nrow; ++i)
    {
      if ((incall || lf->rowsta[i] == lu_rinit))
      {
        for (t = ufir[i], stop = ufir[i] + uuse[i]; t < stop; ++t)
        {
          j = subj[t];

          itemp = ctail[j];
          ctail[j] = t;

          if (itemp == iend)
          {
            chead[j] = t;
            cprev[t] = iend;
          }
          else
          {
            cprev[t] = itemp;
            cnext[itemp] = t;
          }

          abs_aij = fabs(val[t]);
          if (abs_aij > colmax[j])
            colmax[j] = abs_aij;
        }
      }
    }
  }
  else
  {
    for (i = 0; i < nrow; ++i)
    {
      if ((incall || lf->rowsta[i] == lu_rinit))
      {
        for (t = ufir[i], stop = ufir[i] + uuse[i]; t < stop; ++t)
        {
          j = subj[t];

          itemp = ctail[j];
          ctail[j] = t;

          if (itemp == iend)
          {
            chead[j] = t;
            cprev[t] = iend;
          }
          else
          {
            cprev[t] = itemp;
            cnext[itemp] = t;
          }
        }
      }
    }
  }

  for (j = 0; j < ncol; ++j)
    if (ctail[j] != iend)
    {
      ccur[j] = chead[j];
      cnext[ctail[j]] = iend;
    }
} /* lu_buildcol */

static void lu_updcolmax(gsdec *lf,
                         double colmax[],
                         int cmiu[],
                         int j)
{
  double rtemp, cmax = 0.0;

  for (lf->ccur[j] = lf->chead[j]; lf->ccur[j] != lf->wsze;
       lf->ccur[j] = lf->cnext[lf->ccur[j]])
  {
    rtemp = fabs(lf->val[lf->ccur[j]]);
    if (rtemp > cmax)
      cmax = rtemp;
  }

  colmax[j] = cmax;
  cmiu[j] = 1;
} /* lu_updcolmax */

void lu_choosepiv(gsdec *lf,
                  int step,
                  double colmax[],
                  int cmiu[],
                  int *pivi,
                  int *pivj)
{
  int found = false, found_best = false;
  int t, i, j, k, nnzi, nnzj, nnzk, best_i, best_j, stop,
      nsea, curnnz, maxnnz, trhnnz,
      nrow = lf->nrow, trhsch = lf->trhsch,
      iendi = lf->nrow, iendj = lf->ncol,
      *subj = lf->subj,
      *ufir = lf->ufir, *uuse = lf->uuse;
  long trh_mcost, mcost, best_mcost = 0;
  xlist *nzrow = lf->nzrow, *nzcol = lf->nzcol;
  double rtemp, coltrh,
      tolrpiv = lf->tolrpiv, tolapiv = lf->tolapiv,
      *val = lf->val;

  if (infXt(nzcol))
    while (getXt(nzcol, &j, &nnzj) && nnzj == 0)
      delXt(nzcol, j);

  if (getXt(nzcol, &j, &nnzj))
  {
    if (nnzj == 1)
    {
      do
      {
        lf->ccur[j] = lf->chead[j];
        if (fabs(lf->val[lf->ccur[j]]) >= tolapiv)
        {
          *pivi = lf->subi[lf->ccur[j]];
          *pivj = j;
          return;
        }
      } while (succXt(nzcol, &j));

      nnzj = 2;
    }
  }
  else
  {
    *pivi = iendi;
    *pivj = iendj;

    return;
  }

  if (infXt(nzrow))
    while (getXt(nzrow, &i, &nnzi) && nnzi == 0)
      delXt(nzrow, i);

  if (getXt(nzrow, &i, &nnzi))
  {
    if (nnzi == 1)
    {
      {
        j = subj[ufir[i]];
        rtemp = fabs(val[ufir[i]]);

        if (!cmiu[j])
          lu_updcolmax(lf, colmax, cmiu, j);

        if (rtemp >= tolapiv && rtemp >= tolrpiv * colmax[j])
        {
          *pivi = i;
          *pivj = j;
          return;
        }
      }
      while (succXt(nzrow, &i))
        ;

      nnzi = 2;
    }
  }
  else
  {
    *pivi = iendi;
    *pivj = iendj;

    return;
  }

  nsea = 0;
  curnnz = min(nnzi, nnzj);
  maxnnz = 1 + nrow - step;

  trh_mcost = (long)(nnzj - 1) * (long)(nnzi - 1);

  for (found = false, found_best = false; !found && curnnz < maxnnz; ++curnnz)
  {
    if (supXt(nzcol, curnnz))
    {
      getXt(nzcol, &j, &nnzj);
      {
        ++nsea;

        if (!cmiu[j])
          lu_updcolmax(lf, colmax, cmiu, j);

        coltrh = tolrpiv * colmax[j];
        if (coltrh < tolapiv)
          coltrh = tolapiv;

        for (lf->ccur[j] = lf->chead[j];
             lf->ccur[j] != lf->wsze && !found;
             lf->ccur[j] = lf->cnext[lf->ccur[j]])
        {
          if (fabs(lf->val[lf->ccur[j]]) >= coltrh)
          {
            k = lf->subi[lf->ccur[j]];
            nnzk = uuse[k];

            if (curnnz <= nnzk)
            {
              mcost = (long)(curnnz - 1) * (long)(nnzk - 1);
              if (!found_best || mcost < best_mcost)
              {
                found = (mcost <= trh_mcost);
                found_best = true;
                best_mcost = mcost;
                best_i = k;
                best_j = j;
              }
            }
          }
        }

        found |= (found_best && nsea >= trhsch);
      }
      while (!found && succXt(nzcol, &j))
        ;

      nnzj = curnnz + 1;
    }

    trh_mcost = (long)(nnzj - 1) * (long)(nnzi - 1);

    found |= (found_best && best_mcost <= trh_mcost);

    if (!found && supXt(nzrow, curnnz))
    {
      getXt(nzrow, &i, &nnzi);
      if (found_best)
        trhnnz = (int)(best_mcost / ((long)(nnzi - 1)) + 1);
      else
        trhnnz = nrow;

      do
      {
        ++nsea;

        for (t = ufir[i], stop = ufir[i] + uuse[i]; t < stop && !found; ++t)
        {
          k = subj[t];

          if (nzcol->pval[k] == nzcol->idep)
          {
            printf("\n\n system error, column is removed.\n");
            ShutDown();
            exit(0);
          }

          nnzk = nzcol->pval[k];
          if (nnzk > curnnz && nnzk <= trhnnz)
          {
            if (!cmiu[k])
              lu_updcolmax(lf, colmax, cmiu, k);

            rtemp = fabs(val[t]);
            if (rtemp >= tolapiv && rtemp >= tolrpiv * colmax[k])
            {
              mcost = (long)(curnnz - 1) * (long)(nnzk - 1);
              if (!found_best || mcost < best_mcost)
              {
                found = (mcost <= trh_mcost);
                trhnnz = (int)(mcost / ((long)nnzi - 1) + 1);
                found_best = true;
                best_mcost = mcost;
                best_i = i;
                best_j = k;
              }
            }
          }
        }

        found |= (found_best && nsea >= trhsch);

      } while (!found && succXt(nzrow, &i));

      nnzi = curnnz + 1;
    }

    trh_mcost = (long)(nnzj - 1) * (long)(nnzi - 1);

    found |= (found_best && best_mcost <= trh_mcost);
  }

  if (found)
  {
    *pivi = best_i;
    *pivj = best_j;
  }
  else
  {
    *pivi = iendi;
    *pivj = iendj;
  }
} /* lu_choosepiv */

int lu_dopivot(gsdec *lf,
               int step,
               double colmax[],
               int cmiu[],
               int pivi,
               int pivj,
               int zmap[],
               double rmemn1[])
{
  int nnzi, nnzj, jpos, r, t, stop, itemp, mark,
      isec, fill, jtemp, newpos, lcur,
      iend = lf->wsze, iendj = lf->ncol,
      *subi = lf->subi, *subj = lf->subj,
      *cprev = lf->cprev, *cnext = lf->cnext,
      *ufir = lf->ufir, *uuse = lf->uuse, *usze = lf->usze;
  int lresp;
  xlist *nzrow = lf->nzrow, *nzcol = lf->nzcol;
  double uii, mul, vtemp, avtemp,
      toldrop = lf->toldrop,
      *ai, *val = lf->val;
#ifdef TIME_COUNT
//  clock_t  Tm;
//  Tm = GetTime();
#endif

  ai = rmemn1;
  lf->p[step] = pivi;
  nnzi = uuse[pivi];

  if (pivj != iendj)
    nnzj = lf->nzcol->pval[pivj];
  else
    nnzj = 1;

  lf->info.nnzl += nnzj - 1;
  lf->info.nnzu += nnzi;

  /*
   * Expand pivot row.
   */
  if (pivj == iendj)
  {
    for (t = ufir[pivi], stop = ufir[pivi] + uuse[pivi]; t < stop; ++t)
    {
      jtemp = subj[t];

      cmiu[jtemp] = 0;

      lu_delpos(lf, jtemp, t);
      miXt(nzcol, jtemp, 1);
    }
  }
  else
  {
    if (nnzj == 1)
    {
      for (t = ufir[pivi], stop = ufir[pivi] + uuse[pivi]; t < stop; ++t)
      {
        jtemp = subj[t];

        lu_delpos(lf, jtemp, t);
        miXt(nzcol, jtemp, 1);

        if (jtemp == pivj)
          jpos = t;

        cmiu[jtemp] = 0;
      }
    }
    else
    {
      for (t = ufir[pivi], stop = t + uuse[pivi]; t < stop; ++t)
      {
        jtemp = subj[t];

        lu_delpos(lf, jtemp, t);
        miXt(nzcol, jtemp, 1);

        if (jtemp == pivj)
          jpos = t;

        zmap[jtemp] = 1;
        cmiu[jtemp] = 0;
        ai[jtemp] = val[t];
      }
    }

    /* Move pivot element to front. */
    dSwap(ufir[pivi], jpos, val);
    iSwap(ufir[pivi], jpos, subj);

    zmap[pivj] = 0;
    uii = val[ufir[pivi]];
  }
#ifdef TIME_COUNT
//  printf("*** dopivot-A: %d ms\n", GetTime() - Tm);
//  Tm = GetTime();
#endif
  /* Set of row status. */
  lf->rowsta[pivi] = lu_rfin;

  if (pivj == iendj)
    uuse[pivi] = 0;
  else
  {
    if (nnzj > 1)
    {
      if (!lf->dropl)
      {
        /*
         * Setup information for the L.
         */

        if (lf->uend + nnzj - 1 >= lf->lbeg)
        {
          lresp = lu_garcol(lf, step, true, false);
          if (lresp != lu_ok)
            return (lresp);

          if (lf->uend + nnzj - 1 >= lf->lbeg)
            return (lu_wfull);
        }

        lcur = lf->lbeg - (nnzj - 1);
        lf->lbeg = lcur;

        ++lf->lupds;
        cprev[lf->wsze - lf->lupds] = lcur;
        cnext[lf->wsze - lf->lupds] = 0;
        subi[lcur] = pivi;
      }

      /*
       * Perform submatrix update.
       */
      mark = 1;

      for (lf->ccur[pivj] = lf->chead[pivj]; lf->ccur[pivj] != lf->wsze;)
      {
        ++mark;
        r = lf->subi[lf->ccur[pivj]];

        /*
         * Pivot element position.
         */
        jpos = lf->ccur[pivj];

        lu_delcur(lf, pivj);

        /*
         * Multiplier.
         */

        mul = -val[jpos] / uii;
        if (!lf->dropl)
        {
          /*
           * Store multiplier.
           */

          val[lcur] = mul;
          subj[lcur] = r;
          ++lcur;
        }

        /*
         * Delete pivot column element.
         */

        --uuse[r];
        itemp = ufir[r] + uuse[r];
        if (itemp != jpos)
        {
          jtemp = subj[itemp];

          lf->ccur[jtemp] = itemp;
          lu_delcur(lf, jtemp);

          val[jpos] = val[itemp];
          subj[jpos] = jtemp;

          lu_inscur(lf, jtemp, jpos);
        }

        /*
         * Add (mul) times of pivot row to non-pivot row.
         */

        isec = 0;
        for (t = ufir[r], stop = t + uuse[r]; t < stop;)
        {
          jtemp = subj[t];
          if (zmap[jtemp])
          {
            zmap[jtemp] = mark;
            ++isec;
            val[t] += mul * ai[jtemp];
            vtemp = fabs(val[t]);
            if (vtemp < toldrop)
            {
              lf->ccur[jtemp] = t;
              lu_delcur(lf, jtemp);
              miXt(nzcol, jtemp, 1);

              --uuse[r];
              --stop;

              if (stop != t)
              {
                jtemp = subj[stop];

                val[t] = val[stop];
                subj[t] = jtemp;

                lf->ccur[jtemp] = stop;
                lu_delcur(lf, jtemp);

                lu_inscur(lf, jtemp, t);
              }
            }
            else
              ++t;
          }
          else
            ++t;
        }

        fill = nnzi - 1 - isec;
        if (fill == 1)
        {
          for (t = ufir[pivi] + 1, stop = ufir[pivi] + uuse[pivi];
               zmap[subj[t]] == mark && t < stop;
               ++t)
            ;

          newpos = ufir[r] + uuse[r];
          jtemp = subj[t];
          vtemp = mul * ai[jtemp];
          avtemp = fabs(vtemp);
          if (avtemp > toldrop)
          {
            val[newpos] = vtemp;
            subi[newpos] = r;
            subj[newpos] = jtemp;

            lu_inscur(lf, jtemp, newpos);
            plXt(nzcol, jtemp, 1);
            ++uuse[r];
          }
        }
        else if (fill)
        {
          if (uuse[r] + fill > usze[r])
          {
            lresp = lu_incrowsze(lf, step, r, uuse[r] + fill, true, false);
            if (lresp != lu_ok)
              return (lresp);
          }

          newpos = ufir[r] + uuse[r];
          for (t = ufir[pivi] + 1, stop = ufir[pivi] + uuse[pivi]; t < stop; ++t)
          {
            jtemp = subj[t];
            if (zmap[jtemp] != mark)
            {
              vtemp = mul * ai[jtemp];
              avtemp = fabs(vtemp);
              if (avtemp > toldrop)
              {
                val[newpos] = vtemp;
                subi[newpos] = r;
                subj[newpos] = jtemp;

                lu_inscur(lf, jtemp, newpos);
                plXt(nzcol, jtemp, 1);

                ++newpos;
              }
            }
          }
          uuse[r] = newpos - ufir[r];
        }

        /*
         * Update row count.
         */

        putXt(nzrow, r, uuse[r]);
      }

      for (t = ufir[pivi], stop = t + uuse[pivi]; t < stop; ++t)
        zmap[subj[t]] = 0;
    }
  }
#ifdef TIME_COUNT
//  printf("*** dopivot-B: %d ms\n", GetTime() - Tm);
//  Tm = GetTime();
#endif
  /*
   * Remove all information about pivot row and column.
   */

  delXt(nzrow, pivi);

  if (uuse[pivi] && lf->dropu)
    uuse[pivi] = 1;

  if (pivj != iendj)
  {
    delXt(nzcol, pivj);

    lf->chead[pivj] = iend;
    lf->ctail[pivj] = iend;
    lf->ccur[pivj] = iend;
  }
#ifdef TIME_COUNT
//  printf("*** dopivot-C: %d ms\n", GetTime() - Tm);
//  Tm = GetTime();
#endif
  return (lu_ok);
} /* lu_dopivot */

int lu_garcol(gsdec *lf,
              int step,
              int buildcol,
              int incall)
{
  int i, t, k, c, tsze, first,
      itemp, jtemp, cpre, cnex, rnew,
      nrow = lf->nrow,
      iend = lf->wsze, iendi = lf->nrow, iendj = lf->ncol,
      *rord, *rfir, *rmar, *rhead,
      *p = lf->p,
      *subi = lf->subi, *subj = lf->subj,
      *upre = lf->upre, *unex = lf->unex,
      *ufir = lf->ufir, *uuse = lf->uuse, *usze = lf->usze,
      *cprev = lf->cprev, *cnext = lf->cnext;
  double vtemp,
      *val = lf->val;

  ++lf->info.tngc;
  ++lf->info.ngc;

  rord = lf->chead; /* Must not be changed. */
  rmar = upre;
  rfir = unex;

  iZero(nrow, rmar, NULL);

  for (i = 0; i < step; ++i)
  {
    rord[i] = p[i];
    rmar[rord[i]] = 1;
    rfir[rord[i]] = 0;
  }

  c = step;
  for (i = 0; i < nrow; ++i)
  {
    if (!rmar[i])
    {
      rord[c] = i;
      rfir[i] = 0;
      ++c;
    }
  }

  if (c != nrow)
  {
    printf("\n\n internal error, wrong garcol c.\n");
    ShutDown();
    exit(0);
  }

  c = 0;
  for (k = 0; k < nrow; ++k)
  {
    i = rord[k];
    if (ufir[i] == c && uuse[i] + lf->trhuov == usze[i])
      c += usze[i];
    else
      break;
  }

  first = k;

  if (first == nrow)
    return (lu_ok);

  if (first)
    tsze = ufir[rord[first - 1]] + usze[rord[first - 1]];
  else
    tsze = 0;

  rhead = rmar;
  for (k = first; k < nrow; ++k)
  {
    i = rord[k];

    if (uuse[i])
    {
      rhead[i] = ufir[i];
      cprev[ufir[i]] = iend;
      subi[ufir[i]] = i;
      for (t = ufir[i] + 1; t < ufir[i] + uuse[i]; ++t)
      {
        subi[t] = i;

        cprev[t] = t - 1;
        cnext[t - 1] = t;
      }
      cnext[ufir[i] + uuse[i] - 1] = iend;
    }
    else
      rhead[i] = iend;

    for (t = ufir[i] + uuse[i]; t < ufir[i] + usze[i]; ++t)
      subj[t] = iendj;

    usze[i] = uuse[i] + lf->trhuov;
    tsze += usze[i];
  }

  if (tsze > lf->lbeg)
    return (lu_wfull);

  for (t = 0; t < lf->ufir[lf->uhead]; ++t)
    subj[t] = iendj;

  for (t = lf->uend; t < tsze; ++t)
    subj[t] = iendj;

  if (first)
    rfir[rord[first]] = ufir[rord[first - 1]] + usze[rord[first - 1]];
  else
    rfir[rord[first]] = 0;

  for (k = first + 1; k < nrow; ++k)
    rfir[rord[k]] = rfir[rord[k - 1]] + usze[rord[k - 1]];

  for (k = first; k < nrow; ++k)
  {
    uuse[rord[k]] = 0;
    ufir[rord[k]] = rfir[rord[k]];
  }

  if (ufir[rord[nrow - 1]] + usze[rord[nrow - 1]] != tsze)
  {
    printf("\nSize error 1 first=%" IFMT "  nrow=%" IFMT, first, nrow);
    exit(0);
  }

  for (k = first; k < nrow; ++k)
  {
    i = rord[k];
    for (; rhead[i] != iend;)
    {
      rnew = cnext[rhead[i]];

      if (rhead[i] != rfir[i] + uuse[i])
      {
        jtemp = subj[rfir[i] + uuse[i]];
        if (jtemp != iendj)
        {
          vtemp = val[rfir[i] + uuse[i]];
          itemp = subi[rfir[i] + uuse[i]];
          cpre = cprev[rfir[i] + uuse[i]];
          cnex = cnext[rfir[i] + uuse[i]];
        }

        val[rfir[i] + uuse[i]] = val[rhead[i]];
        subi[rfir[i] + uuse[i]] = i;
        subj[rfir[i] + uuse[i]] = subj[rhead[i]];

        if (jtemp == iendj)
          subj[rhead[i]] = iendj;
        else
        {
          val[rhead[i]] = vtemp;
          subi[rhead[i]] = itemp;
          subj[rhead[i]] = jtemp;
          cprev[rhead[i]] = cpre;
          cnext[rhead[i]] = cnex;

          if (cpre == rhead[i])
          {
            rnew = rhead[i];
          }
          else
          {
            if (cpre == iend)
            {
              rhead[itemp] = rhead[i];
            }
            else
              cnext[cpre] = rhead[i];
          }

          if (cnex != iend)
            cprev[cnex] = rhead[i];
        }
      }

      ++uuse[i];

      rhead[i] = rnew;
      if (rnew != iend)
        cprev[rnew] = iend;
    }
  }

  lf->uend = tsze;
  lf->uhead = rord[0];
  lf->utail = rord[nrow - 1];

  lf->unex[rord[0]] = iendi;
  lf->upre[rord[0]] = iendi;

  for (k = 1; k + 1 < nrow; ++k)
  {
    lf->upre[rord[k]] = rord[k - 1];
    lf->unex[rord[k]] = rord[k + 1];
  }

  lf->unex[rord[nrow - 1]] = iendi;
  lf->upre[rord[nrow - 1]] = iendi;

  if (nrow > 1)
  {
    lf->unex[rord[0]] = rord[1];
    lf->upre[rord[nrow - 1]] = rord[nrow - 2];
  }

  if (buildcol)
    lu_buildcol(lf, incall, NULL);

  return (lu_ok);
} /* lu_garcol */

static void lu_moverowtoend(gsdec *lf,
                            int i,
                            int newsze,
                            int updcol)
{
  int t, c, nfir, pre, nex, stop,
      iend = lf->wsze, iendi = lf->nrow, iendj = lf->ncol,
      *subi = lf->subi, *subj = lf->subj,
      *ufir = lf->ufir, *uuse = lf->uuse, *usze = lf->usze,
      *upre = lf->upre, *unex = lf->unex;
  double *val = lf->val;

  /* Compress and transfer row. */
  nfir = lf->uend;
  for (t = ufir[i], c = 0; t < ufir[i] + uuse[i]; ++t)
  {
    val[nfir + c] = val[t];
    subi[nfir + c] = i;
    subj[nfir + c] = subj[t];
    ++c;
  }

  pre = upre[i];
  nex = unex[i];
  if (pre == iendi)
  {
    if (nex == iendi)
      lf->uhead = i;
    else
      lf->uhead = nex;
  }
  else
  {
    usze[pre] += usze[i];

    if (nex == iendi)
      unex[pre] = i;
    else
      unex[pre] = nex;
  }

  if (nex != iendi)
    upre[nex] = pre;

  if (lf->utail == i)
  {
    /* upre[i] and unex[i] is unchanged. */
  }
  else
  {
    upre[i] = lf->utail;
    unex[lf->utail] = i;
  }

  unex[i] = iendi;
  lf->utail = i;

  if (updcol)
  {
    for (t = nfir, stop = nfir + uuse[i]; t < stop; ++t)
    {
      lf->ccur[subj[t]] = ufir[i] + t - nfir;
      lu_delcur(lf, subj[t]);
      lu_inscur(lf, subj[t], t);
    }
  }

  lf->uend += newsze;
  ufir[i] = nfir;
  usze[i] = newsze;
  uuse[i] = c;
} /* lu_moverowtoend */

int lu_incrowsze(gsdec *lf,
                 int step,
                 int i,
                 int newsze,
                 int updcol,
                 int incall)
{
  int lresp;

  if (newsze < 0)
  {
    printf("\n\n system error, newsze=%d<0.\n", newsze);
    ShutDown();
    exit(0);
  }
  if (i < 0 || i >= lf->nrow)
  {
    printf("\n\n system error, %d is negative or too large.\n", i);
    ShutDown();
    exit(0);
  }

  newsze += lf->trhuov;
  if (newsze > lf->uuse[i])
  {
    if (lf->uend + newsze > lf->lbeg)
    {
      lresp = lu_garcol(lf, step, updcol, incall);

      if (lresp != lu_ok)
        return (lresp);

      if (lf->uend + newsze > lf->lbeg)
        return (lu_wfull);
    }
    lu_moverowtoend(lf, i, newsze, updcol);
  }
  else
  {
    printf("\n\n system error, smaller newsze.\n");
    ShutDown();
    exit(0);
  }

  return (lu_ok);
} /* lu_incrowsze */
