#include "LPcross.h"

static double vDots(int n,
                    double *x,
                    int *s,
                    double *y)
{
  int i;
  double r;

  r = 0.0;
  if (n <= 0)
    return r;

  if (!s)
  {
    for (i = 0; i < n; i++)
      r += x[i] * y[i];
  }
  else
  {
    for (i = 0; i < n; i++)
      r += x[i] * y[s[i]];
  }
  return r;
} /* vDots */

void lu_btranu(gsdec *lf,
               int bsze,
               double b[],
               int bsub[],
               int bmap[],
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[])
{
  int i, j, t, k, stop, first, itemp,
      nrow = lf->nrow,
      *subj = lf->subj,
      *ufir = lf->ufir, *uuse = lf->uuse,
      *p = lf->p, *invp = lf->invp, *invq = lf->invq;
  double mul, bj,
      *val = lf->val;

  if (xsub)
    *xnnz = 0;

  if (bsub)
  {
    if (bsze)
    {
      i = invp[invq[bsub[0]]];
      for (t = 1; t < bsze; ++t)
      {
        itemp = invp[invq[bsub[t]]];
        if (i > itemp)
          i = itemp;
      }

      i = p[i];

      if (bsze == 1)
      {
        j = subj[ufir[i]];

        bj = b[j];
        b[j] = 0.0;
        bmap[j] = 0;

        if (xsub)
        {
          while (uuse[i] == 2)
          {
            x[i] = bj / val[ufir[i]];

            bj = -val[ufir[i] + 1] * x[i];

            xsub[*xnnz] = i;
            ++(*xnnz);

            xmap[i] = 1;

            j = subj[ufir[i] + 1];
            i = invq[j];
          }
        }
        else
        {
          while (uuse[i] == 2)
          {
            x[i] = bj / val[ufir[i]];
            bj = -val[ufir[i] + 1] * x[i];
            xmap[i] = 1;

            j = subj[ufir[i] + 1];
            i = invq[j];
          }
        }

        if (uuse[i] == 1)
        {
          x[i] = bj / val[ufir[i]];
          if (xsub)
          {
            xsub[*xnnz] = i;
            ++(*xnnz);
          }
          xmap[i] = 1;

          return;
        }
        else
        {
          b[j] = bj;
          bmap[j] = 1;
        }
      }

      first = invp[i];
    }
    else
      first = nrow;
  }
  else
    first = 0;

  for (t = first; t < nrow; ++t)
  {
    i = p[t];
    if (uuse[i])
    {
      j = subj[ufir[i]];

      if (bmap[j])
      {
        if (xsub)
        {
          xsub[*xnnz] = i;
          ++*xnnz;
        }

        xmap[i] = 1;
        x[i] = b[j] / val[ufir[i]];

        mul = -x[i];

        for (k = ufir[i] + 1, stop = ufir[i] + uuse[i]; k < stop; ++k)
        {
          b[subj[k]] += mul * val[k];
          bmap[subj[k]] = 1;
        }
        b[j] = 0.0;
        bmap[j] = 0;
      }
    }
  }
} /* lu_btranu */

void lu_btranl(gsdec *lf,
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[])
{
  int i, j, t, start, end, itemp,
      *cnext = lf->cnext, *cprev = lf->cprev,
      *subi = lf->subi, *subj = lf->subj;
  double rtemp, *val = lf->val;

  for (i = 0; i < lf->lupds; ++i)
  {
    start = cprev[lf->wsze + i - lf->lupds];
    if (i + 1 < lf->lupds)
      end = cprev[lf->wsze + i + 1 - lf->lupds];
    else
      end = lf->wsze;

    if (cnext[lf->wsze + i - lf->lupds])
    {
      /*
       * A row pivot.
       */

      itemp = subi[start];

      if (xmap[itemp])
      {
        rtemp = x[itemp];
        for (t = start; t < end; ++t)
          x[subj[t]] += rtemp * val[t];
        if (xsub)
        {
          for (t = start; t < end; ++t)
          {
            j = subj[t];
            if (!xmap[j])
            {
              xmap[j] = 1;
              xsub[*xnnz] = j;
              ++*xnnz;
            }
          }
        }
        else
          iSet(end - start, 1, xmap, subj + start);
      }
    }
    else
    {
      /*
       * A column pivot.
       */

      for (t = start; t < end && !xmap[subj[t]]; ++t)
        ;

      if (t < end)
      {
        itemp = subi[start];

        if (!xmap[itemp])
        {
          if (xsub)
          {
            xsub[*xnnz] = itemp;
            ++*xnnz;
          }
          xmap[itemp] = 1;
        }

        /*
         * Change start to t.
         */

        rtemp = 0.0;

        for (; t < end; ++t)
          rtemp += x[subj[t]] * val[t];

        x[itemp] += rtemp;
      }
    }
  }
} /* lu_btranl */

void lu_btran2(gsdec *lf,
               int bsze,
               double b[],
               int bsub[],
               int bmap[],
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[])
{
  if (xsub)
    *xnnz = 0;

  lu_btranu(lf, bsze, b, bsub, bmap, xnnz, x, xsub, xmap);
  lu_btranl(lf, xnnz, x, xsub, xmap);
} /* lu_btran2 */

void lu_btran(gsdec *lf,
              double b[],
              double x[])
{
  int i, t, start, end,
      nrow = lf->nrow,
      *cprev = lf->cprev, *cnext = lf->cnext,
      *subi = lf->subi, *subj = lf->subj,
      *ufir = lf->ufir, *uuse = lf->uuse;
  double *val = lf->val;

  for (t = 0; t < nrow; ++t)
  {
    i = lf->p[t];
    x[i] = b[subj[ufir[i]]] / val[ufir[i]];
    addVect(uuse[i] - 1, -x[i], val + ufir[i] + 1, subj + ufir[i] + 1, b);

    b[subj[ufir[i]]] = 0.0;
  }

  for (i = 0; i < lf->lupds; ++i)
  {
    start = cprev[lf->wsze - lf->lupds + i];
    if (i + 1 < lf->lupds)
      end = cprev[lf->wsze - lf->lupds + i + 1];
    else
      end = lf->wsze;

    if (cnext[lf->wsze - lf->lupds + i])
      addVect(end - start, x[subi[start]], val + start, subj + start, x);
    else
      x[subi[start]] += vDots(end - start, val + start, subj + start, x);
  }
} /* lu_btran */

void lu_ftranl(gsdec *lf,
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[])
{
  int t, i, itemp, jtemp, start, end,
      *subi = lf->subi, *subj = lf->subj,
      *cnext = lf->cnext, *cprev = lf->cprev;
  double rtemp, *val = lf->val;

  end = lf->wsze;
  for (i = 0; i < lf->lupds; ++i)
  {
    start = cprev[lf->wsze - 1 - i];
    if (cnext[lf->wsze - 1 - i])
    {
      /* A row pivot. */
      for (t = start; t < end && !xmap[subj[t]]; ++t)
        ;

      if (t < end)
      {
        itemp = subi[start];
        rtemp = 0.0;
        for (; t < end; ++t)
          if (xmap[subj[t]])
            rtemp += val[t] * x[subj[t]];

        x[itemp] += rtemp;
        if (xsub)
        {
          if (!xmap[itemp])
          {
            xmap[itemp] = 1;
            xsub[*xnnz] = itemp;
            ++*xnnz;
          }
        }
        else
          xmap[itemp] = 1;
      }
    }
    else
    {
      /*
       * A column pivot
       */

      itemp = subi[start];
      if (xmap[itemp])
      {
        if (xsub)
        {
          for (t = start; t < end; ++t)
          {
            jtemp = subj[t];
            if (!xmap[jtemp])
            {
              xmap[jtemp] = 1;
              xsub[*xnnz] = jtemp;
              ++*xnnz;
            }
          }
        }
        else
          iSet(end - start, 1, xmap, subj + start);

        rtemp = x[itemp];
        for (t = start; t < end; ++t)
          x[subj[t]] += rtemp * val[t];
      }
    }
    end = start;
  }
} /* lu_ftranl */

void lu_ftranu(gsdec *lf,
               int *bnnz,
               double b[],
               int bsub[],
               int bmap[],
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[])
{
  int s, t, i, j, itemp, cur, wsze,
      nrow = lf->nrow,
      *subi = lf->subi, *subj = lf->subj,
      *ufir = lf->ufir, *cnext = lf->cnext,
      *p = lf->p, *invp = lf->invp;
  double rtemp, *val = lf->val;

  if (!xsub)
  {
    printf("\n\n argument error, xsub==0.\n");
    ShutDown();
    exit(0);
  }

  if (xsub)
    *xnnz = 0;

  if (bnnz)
  {
    s = 0;
    for (t = 0; t < *bnnz; ++t)
    {
      itemp = invp[bsub[t]];
      if (itemp > s)
        s = itemp;
    }
    ++s;
  }
  else
    s = nrow;

  wsze = lf->wsze;

  for (t = s; t > 0;)
  {
    --t;

    i = p[t];
    j = subj[ufir[i]];
    if (bmap[i])
    {
      xsub[*xnnz] = j;
      ++*xnnz;

      rtemp = b[i] / val[ufir[i]];
      x[j] = rtemp;

      for (cur = lf->chead[j]; cur != wsze; cur = cnext[cur])
      {
        itemp = subi[cur];
        bmap[itemp] = 1;
        b[itemp] -= rtemp * val[cur];
      }

      b[i] = 0.0;
      bmap[i] = 0;
    }
    else
      x[j] = 0.0;
  }

  if (xmap)
    iSet(*xnnz, 1, xmap, xsub);
} /* lu_ftranu */

void lu_ftran(gsdec *lf,
              double b[],
              double x[])
{
  int s, t, i, start, end, itemp, stop,
      nrow = lf->nrow,
      *subi = lf->subi, *subj = lf->subj,
      *cprev = lf->cprev, *cnext = lf->cnext,
      *ufir = lf->ufir, *uuse = lf->uuse;
  double rtemp, *val = lf->val;

  end = lf->wsze;
  for (i = 0; i < lf->lupds; ++i)
  {
    start = cprev[lf->wsze - 1 - i];
    itemp = subi[start];
    rtemp = b[itemp];
    if (cnext[lf->wsze - 1 - i])
    {
      for (t = start; t < end; ++t)
        rtemp += b[subj[t]] * val[t];

      b[itemp] = rtemp;
    }
    else
    {
      for (t = start; t < end; ++t)
        b[subj[t]] += rtemp * val[t];
    }
    end = start;
  }

  for (s = nrow; s;)
  {
    --s;
    i = lf->p[s];

    rtemp = b[i];
    b[i] = 0.0;

    itemp = ufir[i];

    for (t = itemp + 1, stop = itemp + uuse[i]; t < stop; ++t)
      rtemp -= x[subj[t]] * val[t];

    x[subj[itemp]] = rtemp / val[itemp];
  }
} /* lu_ftran */

int lu_initupdate(gsdec *lf)
{
  lf->updates = 0;

  return (lu_ok);
} /* lu_iniupdate */

void perm_rotateleft(int start,
                     int end,
                     int p[])
{
  int itemp, t;

  itemp = p[start];
  for (t = start; t < end; ++t)
    p[t] = p[t + 1];
  p[end] = itemp;
} /* perm_rotateleft */

int lu_supdl(gsdec *lf,
             int i,
             int j,
             int nnzla,
             double la[],
             int subla[],
             int mapla[],
             int inrow1[],
             int incol1[])
{
  int update_of_p, *spike_map;
  int t, jpos, k, l, c, lnnz, lcur, start, stop, end, itemp, jtemp, r, pivj, end0, f,
      pos,
      *spike_sub, *newp,
      nrow = lf->nrow,
      *p = lf->p, *invp = lf->invp, *invq = lf->invq,
      *cprev = lf->cprev, *cnext = lf->cnext,
      *subi = lf->subi, *subj = lf->subj,
      *ufir = lf->ufir, *uuse = lf->uuse, *usze = lf->usze;
  int lresp = lu_ok;
  double mul,
      *val = lf->val, *spike_row;

  ++lf->updates;
  if (lf->updates >= lf->trhrfq)
    lresp = lu_refac;
  else
    ++lf->info.tnnu;

  if (lresp == lu_ok)
  {
    start = invp[invq[i]];

    for (lf->ccur[i] = lf->chead[i]; lf->ccur[i] != lf->wsze;)
    {
      itemp = lf->subi[lf->ccur[i]];
      jpos = lf->ccur[i];

      if (subj[jpos] != i)
      {
        printf("\n\n system error, subj wrong.\n");
        ShutDown();
        exit(0);
      }

      lu_delcur(lf, i);

      --uuse[itemp];
      pos = ufir[itemp] + uuse[itemp];
      if (jpos != pos)
      {
        val[jpos] = val[pos];
        subj[jpos] = subj[pos];

        lf->ccur[subj[jpos]] = pos;
        lu_delcur(lf, subj[jpos]);

        lu_inscur(lf, subj[jpos], jpos);
      }
    }

    end = 0;
    for (t = 0; t < nnzla; ++t)
    {
      k = subla[t];
      if (fabs(la[k]) > lf->toldrop)
      {
        if (uuse[k] == usze[k])
        {
          lresp = lu_incrowsze(lf, nrow, k, usze[k] + 1, true, true);
          if (lresp != lu_ok)
          {
            c = nnzla;
            lresp = lu_refac;
            break;
          }
        }

        val[ufir[k] + uuse[k]] = la[subla[t]];
        subi[ufir[k] + uuse[k]] = k;
        subj[ufir[k] + uuse[k]] = j;

        lu_inscur(lf, j, ufir[k] + uuse[k]);

        ++uuse[k];

        if (end < invp[k])
          end = invp[k];
      }

      la[k] = 0.0;
      mapla[k] = false;
    }

    for (; t < nnzla; ++t)
    {
      k = subla[t];
      la[k] = 0.0;
      mapla[k] = false;
    }

    if (lresp == lu_ok && start > end)
      lresp = lu_insta;

    c = 0;
    lnnz = 0;
    spike_row = la;
    spike_sub = subla;
    spike_map = mapla;

    update_of_p = false;

    if (lresp == lu_ok)
    {
      end0 = end;
      newp = inrow1;
      r = p[start];

      if (start == end)
        f = end + 1;
      else
      {
        f = nrow;
        for (t = ufir[r], stop = ufir[r] + uuse[r]; t < stop; ++t)
          if (subj[t] != j && invp[invq[subj[t]]] < f)
            f = invp[invq[subj[t]]];
      }

      if (f <= start)
      {
        printf("\n f is smaller than start");

        printf("\n r=%d"
               " invp[r]=%d"
               "  f=%d"
               "  start=%d",
               r, lf->invp[r], f, start);

        printf("\n end=%d", end);
        printf("  nrow=%d", nrow);

        for (t = 0; t < nrow; ++t)
        {
          if (uuse[p[t]] && invq[subj[ufir[p[t]]]] != p[t])
          {
            printf("\n invq(t=%d"
                   ",p[t]=%d"
                   ",j=%d"
                   ")",
                   t, p[t], subj[ufir[p[t]]]);
          }
          else
            printf(" Null column ");
        }

        printf("\n\n system error, f<=start.\n");
        ShutDown();
        exit(0);
      }

      if (f > end)
      {
        t = getPos(uuse[r], j, subj + ufir[r]);
        itemp = ufir[r] + t;

        if (t == uuse[r] || fabs(val[itemp]) <= lf->tolapiv)
          lresp = lu_insta;
        else
        {
          invq[j] = r;
          if (t > 0)
          {
            lf->ccur[j] = itemp;
            lu_delcur(lf, j);

            jtemp = subj[ufir[r]];
            lf->ccur[jtemp] = ufir[r];
            lu_delcur(lf, jtemp);

            dSwap(ufir[r], itemp, val);
            iSwap(ufir[r], itemp, subj);

            lu_inscur(lf, j, ufir[r]);
            lu_inscur(lf, jtemp, itemp);
          }
        }
      }
      else
      {
        update_of_p = true;

        --f;
      }

      perm_rotateleft(start, end, p);

      if (update_of_p)
        for (t = start; t < f; ++t)
          newp[t] = p[t];

      /*
       * Main iteration loop.
       */

      for (k = f; k < end0 && lresp == lu_ok;)
      {
        /*
         * Unpacking of pivot row.
         */

        r = p[end0];
        c = 0;
        for (t = ufir[r], stop = ufir[r] + uuse[r]; t < stop; ++t)
        {
          jtemp = subj[t];
          spike_map[jtemp] = 1;
          spike_row[jtemp] = val[t];
          spike_sub[c] = jtemp;
          ++c;
        }

        /*
         * Is index of column at the diagonal
         * of the spike row.
         */

        pivj = j;

        lnnz = 0;

        for (; k < end0; ++k)
        {
          itemp = p[k];
          jtemp = subj[ufir[itemp]]; /* First column in row k. */

          if (!uuse[itemp])
          {
            printf("\n\n system error, uuse wrong.\n");
            ShutDown();
            exit(0);
          }

          if (itemp == r)
          {
            printf("\n\n system error, r wrong.\n");
            ShutDown();
            exit(0);
          }

          newp[k - (end0 - end)] = itemp;
          if (spike_map[jtemp])
          {
            if (fabs(spike_row[jtemp]) > lf->toldrop)
            {
              if (uuse[itemp] == 1 ||
                  (uuse[itemp] == 2 &&
                   invp[invq[subj[ufir[itemp] + 1]]] > end0))
              {
                newp[end] = itemp;
                end--;
              }
              else
              {
                if (fabs(val[ufir[itemp]]) >=
                    lf->tolupiv * fabs(spike_row[jtemp]))
                {
                  spike_map[jtemp] = false;
                  incol1[lnnz++] = k - (end0 - end);

                  mul = -spike_row[jtemp] / val[ufir[itemp]];
                  spike_row[jtemp] = mul;

                  for (t = ufir[itemp] + 1, stop = ufir[itemp] + uuse[itemp]; t < stop; ++t)
                  {
                    l = subj[t];
                    spike_row[l] += mul * val[t];

                    if (!spike_map[l])
                    {
                      spike_map[l] = true;
                      spike_sub[c++] = l;
                    }
                  }
                }
                else
                {
                  pivj = jtemp;
                  iSwap(k, end0, p);

                  break;
                }
              }
            }
            else
            {
              spike_row[jtemp] = 0.0;
              spike_map[jtemp] = false;
            }
          }
        }
        if (lresp == lu_ok)
        {
          if (lf->uend + lnnz > lf->lbeg)
            lresp = lu_garcol(lf, nrow, true, true);

          if (lresp != lu_ok || lf->uend + lnnz > lf->lbeg)
          {
            lf->incwsze = true;
            lresp = lu_refac;
          }
          else
          {
            if (lnnz)
            {
              lcur = lf->lbeg - lnnz;

              ++lf->lupds;
              cprev[lf->wsze - lf->lupds] = lcur;
              cnext[lf->wsze - lf->lupds] = 1;

              lf->lbeg = lcur;
              subi[lcur] = r;
            }
          }

          if (lresp == lu_ok && c - lnnz > usze[r])
          {
            lresp = lu_incrowsze(lf, nrow, r, c - lnnz, true, true);

            if (lresp != lu_ok)
            {
              lf->incwsze = true;
              lresp = lu_refac;
            }
          }

          if (lresp == lu_ok)
          {
            if (fabs(spike_row[pivj]) <= lf->tolapiv)
              lresp = lu_insta;

            newp[k - (end0 - end)] = r;

            invq[pivj] = r;

            if (lnnz)
            {
              for (t = ufir[r], stop = ufir[r] + uuse[r]; t < stop; ++t)
              {
                jtemp = subj[t];
                lf->ccur[jtemp] = t;
                lu_delcur(lf, jtemp);
              }

              val[ufir[r]] = spike_row[pivj];
              subi[ufir[r]] = r;
              subj[ufir[r]] = pivj;

              spike_map[pivj] = false;
              spike_row[pivj] = 0.0;

              lu_inscur(lf, pivj, ufir[r]);

              itemp = ufir[r] + 1;
              for (t = 0; t < c; ++t)
              {
                jtemp = spike_sub[t];

                if (spike_map[jtemp])
                {
                  spike_map[jtemp] = false;
                  if (fabs(spike_row[jtemp]) > lf->toldrop)
                  {
                    val[itemp] = spike_row[jtemp];
                    subi[itemp] = r;
                    subj[itemp] = jtemp;

                    lu_inscur(lf, jtemp, itemp); /* Append. */
                    ++itemp;
                  }
                  spike_row[jtemp] = 0.0;
                }
              }
              uuse[r] = itemp - ufir[r];

              for (t = 0; t < lnnz; ++t, ++lcur)
              {
                itemp = newp[incol1[t]];

                jtemp = subj[ufir[itemp]];

                val[lcur] = spike_row[jtemp];
                subj[lcur] = itemp;

                spike_row[jtemp] = 0.0;
              }
            }
            else
            {
              itemp = ufir[r] + uuse[r];
              for (t = ufir[r], stop = ufir[r] + uuse[r]; t < stop; ++t)
              {
                jtemp = subj[t];
                spike_row[jtemp] = 0.0;
                spike_map[jtemp] = false;

                if (jtemp == pivj)
                  itemp = t;
              }

              if (itemp == ufir[r] + uuse[r])
              {
                printf("\n\n system error, end wrong.\n");
                ShutDown();
                exit(0);
              }

              t = itemp;

              if (t > ufir[r])
              {
                lf->ccur[pivj] = t;
                lu_delcur(lf, pivj);

                jtemp = subj[ufir[r]];
                lf->ccur[jtemp] = ufir[r];
                lu_delcur(lf, jtemp);

                dSwap(ufir[r], t, val);
                iSwap(ufir[r], t, subj);

                lu_inscur(lf, pivj, ufir[r]);
                lu_inscur(lf, jtemp, t);
              }
            }
          }
        }
      }
    }

    if (lresp == lu_ok)
    {
      /*
       * Update invp.
       */

      if (update_of_p)
      {
        for (t = start; t <= end0; ++t)
        {
          k = newp[t];
          p[t] = k;
          invp[k] = t;
        }
      }
      else
      {
        for (t = start; t <= end0; ++t)
        {
          k = p[t];
          invp[k] = t;
        }
      }
    }
    else
    {
      dZero(c, spike_row, spike_sub);
      iZero(c, spike_map, spike_sub);

      for (t = 0; t < lnnz; ++t)
      {
        itemp = newp[incol1[t]];
        jtemp = subj[ufir[itemp]];

        spike_row[jtemp] = 0.0;
      }
    }
  }
  else
  {
    dZero(nnzla, la, subla);
    iZero(nnzla, mapla, subla);
  }

  return (lresp);
} /* lu_supdl */

lu_resp lu_supd(gsdec *lf,
                int i,
                int j,
                int nnzaj,
                double aj[],
                int subaj[],
                int mapaj[],
                int imemm1[],
                int imemn1[])
{
  if (lf->updates + 1 < lf->trhrfq)
  {
    lu_ftranl(lf, &nnzaj, aj, subaj, mapaj);
    return (lu_supdl(lf, i, j, nnzaj, aj, subaj, mapaj, imemm1, imemn1));
  }
  else
  {
    dZero(nnzaj, aj, subaj);
    iZero(nnzaj, mapaj, subaj);

    return (lu_refac);
  }
} /* lu_supd */
