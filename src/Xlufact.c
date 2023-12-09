#include "LPcross.h"
//#define TIME_COUNT

static int lu_cinit(gsdec *lf)
{
  int i, nrow = lf->nrow;

  lf->incwsze = false;
  lf->lbeg = lf->wsze;
  lf->uend = 0;

  lf->lupds = 0;

  for (i = 0; i < nrow; ++i)
    lf->rowsta[i] = lu_rinit;

  lf->info.nnzl = 0;
  lf->info.nnzu = 0;
  lf->info.ngc = 0;

  return (lu_ok);
} /* lu_cinit */

int lu_init(gsdec *lf)
{
  int nrow = lf->nrow, ncol = lf->ncol, nil = lf->wsze,
      *firi = lf->upre, *firj = lf->chead,
      *nnzi = lf->uuse, *nnzj = lf->ctail;

  lu_cinit(lf);

  iSet(nrow, nil, firi, NULL);
  iSet(ncol, nil, firj, NULL);

  iZero(nrow, nnzi, NULL);
  iZero(ncol, nnzj, NULL);

  return (lu_ok);
} /* lu_init */

int lu_inaij(gsdec *lf,
             double aij,
             int i,
             int j)
{
  int t,
      *subi = lf->subi, *subj = lf->subj,
      *lnki = lf->cprev, *lnkj = lf->cnext,
      *firi = lf->upre, *firj = lf->chead,
      *nnzi = lf->uuse, *nnzj = lf->ctail;
  double *val = lf->val;

  if (fabs(aij) > lf->toldrop)
  {
    if (lf->uend < lf->wsze)
    {
      t = lf->uend;

      lnki[t] = firi[i];
      firi[i] = t;
      nnzi[i]++;

      lnkj[t] = firj[j];
      firj[j] = t;
      nnzj[j]++;

      val[t] = aij;
      subi[t] = i;
      subj[t] = j;

      lf->uend++;
    }
    else
      return (lu_wfull);
  }
  return (lu_ok);
} /* lu_inaij */

int lu_inaj(gsdec *lf,
            int j,
            int sze,
            double valj[],
            int subj[],
            int mapj[])
{
  int i, t;
  int lresp = lu_ok;

  if (lf->uend + sze < lf->wsze)
  {
    if (mapj)
    {
      for (t = 0; t < sze && lresp == lu_ok; ++t)
      {
        i = subj[t];
        if (mapj[i])
          lresp = lu_inaij(lf, valj[t], i, j);
      }
    }
    else
      for (t = 0; t < sze && lresp == lu_ok; ++t)
        lresp = lu_inaij(lf, valj[t], subj[t], j);
  }
  else
    lresp = lu_wfull;

  return (lresp);
} /* lu_inaj */

int lu_inaj2(gsdec *lf,
             int j,
             int sze,
             double valj[],
             int subj[],
             int cutoff,
             int inew[])
{
  int i, t, itemp,
      *lsubi = lf->subi, *lsubj = lf->subj,
      *lnki = lf->cprev, *lnkj = lf->cnext,
      *firi = lf->upre, *firj = lf->chead,
      *nnzi = lf->uuse, *nnzj = lf->ctail;
  int lresp = lu_ok;
  double *lval = lf->val;

  itemp = lf->uend;
  if (itemp + sze < lf->wsze)
  {
    for (t = 0; t < sze; ++t)
    {
      i = inew[subj[t]];
      if (fabs(valj[t]) > lf->toldrop && i < cutoff)
      {
        lnki[itemp] = firi[i];
        firi[i] = itemp;
        nnzi[i]++;

        lnkj[itemp] = firj[j];
        firj[j] = itemp;
        nnzj[j]++;

        lval[itemp] = valj[t];
        lsubi[itemp] = i;
        lsubj[itemp] = j;

        itemp++;
      }
    }
    lf->uend = itemp;
  }
  else
    lresp = lu_wfull;

  return (lresp);
} /* lu_inaj2 */

void lu_restore(gsdec *lf)
{
  int i;

  for (i = 0; i < lf->nrow; ++i)
  {
    lf->invp[lf->p[i]] = i;

    if (lf->uuse[i])
      lf->invq[lf->subj[lf->ufir[i]]] = i;
  }

  lu_buildcol(lf, true, NULL);

} /* lu_restore */

int lu_insart(gsdec *lf,
              int insj[],
              int wisze,
              int wimem[])
{
  int i, j, t, r, pos, itemp,
      nrow = lf->nrow, ncol = lf->ncol, nilj = ncol,
      *subi = lf->subi, *subj = lf->subj,
      *ufir = lf->ufir, *uuse = lf->uuse, *usze = lf->usze,
      *p = lf->p,
      *keepj;
  double uii,
      *val = lf->val;

  for (t = 0; t < lf->ncol; ++t)
    for (lf->ccur[t] = lf->chead[t]; lf->ccur[t] != lf->wsze;
         lf->ccur[t] = lf->cnext[lf->ccur[t]])
      if (lf->subj[lf->ccur[t]] != t)
      {
        printf("\n\n system error, incorrect column link.\n");
        ShutDown();
        exit(0);
      }

  if (wisze < ncol)
    return (lu_wispc);

  keepj = wimem;

  iZero(ncol, keepj, NULL);

  for (i = 0; i < nrow; ++i)
  {
    lu_getuii(lf, i, &j, &uii);
    if (j != nilj)
    {
      if (keepj[j])
      {
        printf("\n\n system error, double keep.\n");
        ShutDown();
        exit(0);
      }
      keepj[j] = 1;
    }
  }

  iSet(ncol, ncol, insj, NULL);

  t = 0;
  for (i = 0; i < nrow; ++i)
  {

    lu_getuii(lf, i, &j, &uii);
    if (j == nilj)
    {
      for (; keepj[t] && t < ncol; ++t)
        ;

      if (t == ncol)
      {
        printf("\n\n system error, t==ncol.\n");
        ShutDown();
        exit(0);
      }

      /*
       * Delete column t.
       */

      for (lf->ccur[t] = lf->chead[t]; lf->ccur[t] != lf->wsze;)
      {
        r = lf->subi[lf->ccur[t]];

        if (lf->subj[lf->ccur[t]] != t)
        {
          printf("\n\n system error, incorrect column link.\n");
          ShutDown();
          exit(0);
        }

        pos = lf->ccur[t];

        lu_delcur(lf, t);

        if (uuse[r] == 0)
        {
          printf("\n\n system error, uuse wrong.\n");
          ShutDown();
          exit(0);
        }

        --uuse[r];
        itemp = ufir[r] + uuse[r];

        if (pos != itemp)
        {
          /*
           * Remove element.
           */

          lf->ccur[subj[itemp]] = itemp;
          lu_delcur(lf, subj[itemp]);

          subj[pos] = subj[itemp];
          val[pos] = val[itemp];

          lu_inscur(lf, subj[pos], pos);
        }
      }

      /*
       * Insert an artificial variable
       * on column t.
       */

      if (!usze[i])
      {
        printf("\n\n system error, no space.\n");
        ShutDown();
        exit(0);
      }

      if (uuse[i])
      {
        printf("\n\n system error, empty row.\n");
        ShutDown();
        exit(0);
      }

      uuse[i] = 1;
      subj[ufir[i]] = t;
      subi[ufir[i]] = i;
      val[ufir[i]] = 1.0;

      lf->ccur[t] = lf->chead[t];
      lu_inscur(lf, t, ufir[i]);

      insj[t] = i;
      ++t;
    }
  }

  lu_restore(lf);

  for (t = 0; t < nrow; ++t)
    if (lf->p[lf->invp[t]] != t)
    {
      printf("\n\n system error, wrong p or invp.\n");
      ShutDown();
      exit(0);
    }

  for (t = 0; t < nrow; ++t)
    if (lf->invq[lf->subj[lf->ufir[lf->p[t]]]] != lf->p[t])
    {
      if (lf->uuse[lf->p[t]] == 0)
      {
        printf("\n\n internal error.\n");
        ShutDown();
        exit(0);
      }

      printf("\nt=%" IFMT "  p[t]=%" IFMT, t, p[t]);
      printf("\ninvq[subj[ufir[p[t]]]]=%" IFMT,
             lf->invq[subj[ufir[p[t]]]]);
      printf("  subj=%" IFMT, lf->subj[lf->ufir[p[t]]]);

      printf("\n\n system error, wrong invp.\n");
      ShutDown();
      exit(0);
    }

  return (lu_ok);
} /* lu_insart */

static int lu_setup(gsdec *lf,
                    int *step,
                    int nnzj[],
                    int wrsze,
                    double wrmem[])
{
  int i, j, k, s, t, sze0, sze, pivi, pivj, cur, pos, utrisze, itemp,
      stop, jtemp, inew, jnew, left, lbeg, uend,
      trhuov = lf->trhuov, nrow = lf->nrow, ncol = lf->ncol,
      nil = lf->wsze, nili = nrow,
      *firi, *firj, *lnki, *lnkj, *nnzi, *list, *q, *p,
      *lcur, *lfir, *ucur,
      *marki, *markj,
      *subi = lf->subi, *subj = lf->subj,
      *cprev = lf->cprev, *cnext = lf->cnext,
      *ufir = lf->ufir, *uuse = lf->uuse, *usze = lf->usze;
  int *rowsta = lf->rowsta;
  double rtemp, cmax, absaij, vnew, v,
      tolapiv = lf->tolapiv, tolrpiv = lf->tolrpiv,
      *val = lf->val,
      *uij;

  lf->info.nnza = lf->uend;
  lf->info.tnnza += lf->uend;

  lf->info.nnzl = 0;
  lf->info.nnzu = 0;

  lf->lupds = 0;
  lf->lbeg = lf->wsze;

  if (!nrow)
    return (lu_ok);

  if (lf->uend + trhuov * nrow > lf->wsze)
    return (lu_wfull);

  if (wrsze < ncol)
    return (lu_wrspc);

  firi = lf->upre;
  firj = lf->chead;
  nnzi = lf->unex;
  list = lf->ctail;

  iCopy(nrow, lf->uuse, nnzi);
  iCopy(ncol, lf->ctail, nnzj);

  p = lf->p;
  q = lf->invq;

  lnki = cprev;
  lnkj = cnext;

  /*
   * Columns triangulazation.
   */
  sze = 0;
  for (j = 0; j < ncol; ++j)
    if (nnzj[j] == 1)
      list[sze++] = j;

  while (sze)
  {
    --sze;
    pivj = list[sze];

    if (nnzj[pivj] == 1)
    {
      for (cur = firj[pivj]; cur != nil &&
                             rowsta[subi[cur]] != lu_rinit;
           cur = lnkj[cur])
        ;

      if (cur == nil)
      {
        printf("\n\n system error, not a sing step.\n");
        ShutDown();
        exit(0);
      }

      pos = cur;
      pivi = subi[pos];

      sze0 = sze;
      for (cur = firi[pivi]; cur != nil; cur = lnki[cur])
      {
        j = subj[cur];
        --nnzj[j];

        if (nnzj[j] == 1)
          list[sze++] = j;
      }

      rtemp = val[pos];
      if (fabs(rtemp) >= tolapiv)
      {
        rowsta[pivi] = lu_rfin;

        firj[pivj] = nil;

        p[*step] = pivi;
        q[*step] = pivj;

        ++*step;
      }
      else
      {
        for (cur = firi[pivi]; cur != nil; cur = lnki[cur])
          ++nnzj[subj[cur]];

        sze = sze0;
      }
    }
  }

  utrisze = *step;

  uij = wrmem;

  if (*step < nrow)
  {
    /*
     * Row triangulazation.
     */

    sze = 0;
    for (i = 0; i < nrow; ++i)
      if (rowsta[i] == lu_rinit && nnzi[i] == 1)
        list[sze++] = i;

    while (sze)
    {
      --sze;
      pivi = list[sze];

      if (nnzi[pivi] == 1)
      {
        for (cur = firi[pivi]; cur != nil &&
                               firj[subj[cur]] == nil;
             cur = lnki[cur])
          ;

        if (cur == nil)
        {
          printf("\nnrow=%" IFMT, nrow);
          printf("  cur=%" IFMT, cur);
          printf("\npivi=%" IFMT "  subi[cur]=%" IFMT, pivi, subi[cur]);
          printf("\nnil=%" IFMT "  firi[pivi]=%" IFMT,
                 nil, firi[pivi]);
          printf(" nnzi[pivi]=%" IFMT, nnzi[pivi]);

          for (cur = firi[pivi]; cur != nil &&
                                 nnzj[subj[cur]] == nrow + 1;
               cur = lnki[cur])
            printf(" j=%" IFMT, subj[cur]);

          printf("\n\n system error, cur==nil.\n");
          ShutDown();
          exit(0);
        }

        if (subi[cur] != pivi)
        {
          printf("\nsubi[cur]=%" IFMT, subi[cur]);
          printf("\nnrow=%" IFMT, nrow);
          printf("\nnil=%" IFMT "  firi[pivi]=%" IFMT "  pivi=%" IFMT,
                 nil, firi[pivi], pivi);
          printf(" nnzi[pivi]=%" IFMT, nnzi[pivi]);

          for (cur = firi[pivi]; cur != nil; cur = lnki[cur])
          {
            if (firj[subj[cur]] != nil)
              printf(" j=%" IFMT, subj[cur]);
          }
          printf("\n\n system error, subi[cur]!=pivi.\n");
          ShutDown();
          exit(0);
        }

        pos = cur;
        pivj = subj[pos];
        rtemp = fabs(val[pos]);

        rowsta[pivi] = lu_rfin;
        itemp = 0;
        cmax = 0.0;
        sze0 = sze;
        for (cur = firj[pivj]; cur != nil; cur = lnkj[cur])
        {
          i = subi[cur];
          if (rowsta[i] == lu_rinit)
          {
            ++itemp;

            absaij = fabs(val[cur]);
            if (absaij > cmax)
              cmax = absaij;

            --nnzi[i];
            if (nnzi[i] == 1)
              list[sze++] = i;
          }
        }

        if (rtemp >= tolapiv && rtemp >= tolrpiv * cmax && itemp > 0)
        {
          uij[pivj] = val[pos];

          firj[pivj] = nil;
          nnzj[pivj] = itemp;

          p[*step] = pivi;
          q[*step] = pivj;

          ++*step;
        }
        else
        {
          /*
           * Update nnzi
           */

          for (cur = firj[pivj]; cur != nil; cur = lnkj[cur])
          {
            i = subi[cur];
            if (rowsta[i] == lu_rinit)
              ++nnzi[i];
          }

          rowsta[pivi] = lu_rinit;
          sze = sze0;
        }
      }
    }
  }

  /*
   * Setup data-structures.
   */
  uend = lf->uend;
  lbeg = lf->lbeg;

  marki = firi;
  markj = firj;

  itemp = *step;
  if (*step < nrow)
  {
    for (t = 0; t < nrow; ++t)
      if (rowsta[t] == lu_rinit)
        p[itemp++] = t;
  }

  if (utrisze < *step)
  {
    /*
     * Separate L and U.
     */

    for (t = 0; t < nrow; ++t)
      marki[p[t]] = t;

    iSet(ncol, nrow, markj, NULL);

    for (t = 0; t < *step; ++t)
      markj[q[t]] = t;

    for (t = 0; t < uend;)
    {
      i = subi[t];
      j = subj[t];

      if (marki[i] > markj[j])
      {
        /*
         * L element.
         */

        v = val[t];

        for (--uend; t < uend;)
        {
          inew = subi[uend];
          jnew = subj[uend];

          if (marki[inew] > markj[jnew])
          {
            /*
             * L element.
             */

            --lbeg;

            subi[lbeg] = inew;
            subj[lbeg] = jnew;
            val[lbeg] = -val[uend] / uij[jnew];

            --uend;
          }
          else
            break;
        }

        subi[t] = subi[uend];
        subj[t] = subj[uend];
        val[t] = val[uend];

        --lbeg;
        subi[lbeg] = i;
        subj[lbeg] = j;
        val[lbeg] = -v / uij[j];
      }
      else
        ++t;
    }

    lf->uend = uend;
    lf->lbeg = lbeg;
  }

  if (lf->info.nnza != uend + lf->wsze - lf->lbeg)
  {
    printf("\nuend=%" IFMT "  wsze=%" IFMT,
           uend, lf->wsze);
    printf("\nlbeg=%" IFMT, lbeg);
    printf("\nannz=%" IFMT, lf->info.nnza);
    printf("\n\n system error, lu size wrong.\n");
    ShutDown();
    exit(0);
  }

  lf->info.nnzl = lf->wsze - lf->lbeg;
  lf->info.nnzu = lf->info.nnza - lf->info.nnzl;

  /*
   * Setup of L.
   */

  if (utrisze < *step && !lf->dropl)
  {
    lfir = markj;
    lcur = list;

    j = q[utrisze];
    lfir[j] = lf->wsze - nnzj[j];
    lcur[j] = 0;

    for (t = utrisze + 1; t < *step; ++t)
    {
      j = q[t];
      lfir[j] = lfir[q[t - 1]] - nnzj[j];
      lcur[j] = 0;
    }

    for (t = utrisze; t < *step; ++t)
    {
      j = q[t];

      if (nnzj[j] == 0)
        exit(0);

      cur = lfir[j];
      stop = cur + nnzj[j];

      ++lf->lupds;
      cprev[lf->wsze - lf->lupds] = cur;
      cnext[lf->wsze - lf->lupds] = 0;

      for (cur += lcur[j]; cur < stop; ++cur)
      {
        itemp = subi[cur];
        jtemp = subj[cur];
        rtemp = val[cur];

        while (jtemp != j)
        {
          s = lfir[jtemp] + lcur[jtemp];
          ++lcur[jtemp];

          vnew = val[s];
          val[s] = rtemp;
          rtemp = vnew;

          inew = subi[s];
          subi[s] = jtemp;

          jnew = subj[s];
          subj[s] = itemp;

          itemp = inew;
          jtemp = jnew;
        }

        val[cur] = rtemp;
        subi[cur] = jtemp;
        subj[cur] = itemp;
      }

      lcur[j] = nnzj[j];
    }

    /*
     * Setup af subi.
     */
    for (t = utrisze; t < *step; ++t)
      subi[lfir[q[t]]] = p[t];
  }

  /*
   * Setup of U.
   */

  ucur = lf->invp;

  /*
   * marki[i] : Pivot column on ith row (otherwise = ncol).
   */
  for (t = 0; t < *step; ++t)
    marki[p[t]] = q[t];

  for (; t < nrow; ++t)
    marki[p[t]] = ncol;

  left = lf->uend;
  lf->uend = left + nrow * trhuov;

  i = p[0];
  ufir[i] = 0;
  ucur[i] = 0;
  uuse[i] = nnzi[i];
  usze[i] = nnzi[i] + trhuov;
  nnzi[i] = min(usze[i], left);
  left -= nnzi[i];
  j = i;

  for (t = 1; t < nrow; ++t)
  {
    i = p[t];
    ufir[i] = ufir[j] + usze[j];
    ucur[i] = 0;
    uuse[i] = nnzi[i];
    usze[i] = nnzi[i] + trhuov;
    nnzi[i] = min(usze[i], left);
    left -= nnzi[i];
    j = i;
  }

  /*
   * Setup row structure of U.
   */
  for (t = 0; t < nrow; ++t)
  {
    i = p[t];

    cur = ufir[i] + ucur[i];
    stop = ufir[i] + nnzi[i];

    for (; cur < stop;)
    {
      itemp = subi[cur];

      if (itemp == i)
      {
        if (subj[cur] == marki[i] && ucur[i])
        {
          /*
           * Located a pivot element. Move it to front.
           */

          iSwap(ufir[i], cur, subj);
          dSwap(ufir[i], cur, val);
        }
        ++ucur[i];
        ++cur;
      }
      else
      {
        rtemp = val[cur];
        jtemp = subj[cur];

        while (itemp != i)
        {
          s = ufir[itemp] + ucur[itemp];

          if (ucur[itemp] < nnzi[itemp])
          {
            if (jtemp == marki[itemp] && ucur[itemp])
            {
              /*
               * Located a pivot element.
               * Move it to the front of the row.
               */

              k = ufir[itemp];

              vnew = val[k];
              jnew = subj[k];

              val[k] = rtemp;
              subj[k] = jtemp;

              rtemp = vnew;
              jtemp = jnew;

              /*
               * itemp is unchanged.
               */
            }

            ++ucur[itemp];

            vnew = val[s];
            val[s] = rtemp;
            rtemp = vnew;

            inew = subi[s];
            subi[s] = itemp;
            itemp = inew;

            jnew = subj[s];
            subj[s] = jtemp;
            jtemp = jnew;
          }
          else
          {
            /*
             * Empty position has been located.
             */
            if (jtemp == marki[itemp] && ucur[itemp])
            {
              /*
               * Located a pivot element.
               */

              k = ufir[itemp];

              val[s] = val[k];
              subi[s] = subi[k];
              subj[s] = subj[k];

              s = k;
            }

            val[s] = rtemp;
            subi[s] = itemp;
            subj[s] = jtemp;

            ++ucur[itemp];

            --nnzi[i];
            --stop;
            rtemp = val[stop];
            itemp = subi[stop];
            jtemp = subj[stop];

            break;
          }
        }

        val[cur] = rtemp;
        subi[cur] = itemp;
        subj[cur] = jtemp;
      }
    }
  }

  if (lf->uend >= lf->lbeg ||
      lf->uend != lf->ufir[p[nrow - 1]] + lf->usze[p[nrow - 1]])
  {
    printf("\n\n uend error uend=%" IFMT, lf->uend);
    printf("\n test=%" IFMT, lf->ufir[p[nrow - 1]] + lf->usze[p[nrow - 1]]);
    ShutDown();
    exit(0);
  }

  if (lf->dropu)
  {
    for (t = 0; t < *step; ++t)
      uuse[p[t]] = 1;
  }

  /*
   * Setup links.
   */

  lf->uhead = p[0];
  lf->utail = p[nrow - 1];

  lf->unex[p[0]] = nili;
  lf->upre[p[0]] = nili;

  for (i = 1; i + 1 < nrow; ++i)
  {
    lf->upre[p[i]] = p[i - 1];
    lf->unex[p[i]] = p[i + 1];
  }

  lf->unex[p[nrow - 1]] = nili;
  lf->upre[p[nrow - 1]] = nili;

  if (nrow > 1)
  {
    lf->unex[p[0]] = p[1];
    lf->upre[p[nrow - 1]] = p[nrow - 2];
  }

  return (lu_ok);
} /* lu_setup */

int lu_factor(gsdec *lu,
              int wisze,
              int wimem[],
              int wrsze,
              double wrmem[])
{
  int step, pivi, pivj, i, itemp, j, wiszel,
      nrow = lu->nrow, ncol = lu->ncol,
      nili = nrow,
      *zmap, *nnzaj, *cmiu, *wimeml;
  int lresp;
  double *colmax, *rbn1;
#ifdef TIME_COUNT
  clock_t Tm;
  int CntA = 0, CntB = 0, CntPivot = 0;
  Tm = GetTime();
#endif

  if (wisze < 2 + 5 * nrow + 6 * ncol)
    return (lu_wispc);

  if (wrsze < 2 * nrow)
    return (lu_wrspc);

  if (ncol != nrow)
  {
    printf("\n\n system error, nrow!=ncol.\n");
    ShutDown();
    exit(0);
  }

  lu->info.nnzl = 0;
  lu->info.nnzu = 0;

  step = 0;

  wiszel = wisze;
  wimeml = wimem;

  nnzaj = wimeml;
  cmiu = nnzaj + ncol;

  wiszel -= 2 * ncol;
  wimeml += 2 * ncol;

  lresp = lu_setup(lu, &step, nnzaj, wrsze, wrmem);

  if (step < nrow)
  {
    for (i = 0; i < step; ++i)
      if (lu->uuse[lu->p[i]])
        nnzaj[lu->subj[lu->ufir[lu->p[i]]]] = nrow + 1;
  }

  if (!setXt(lu->nzrow, nrow, ncol, wiszel, wimeml))
  {
    printf("\n\n system error, pque setup.\n");
    ShutDown();
    exit(0);
  }

  itemp = 1 + 3 * nrow + 1 * ncol;

  wiszel -= itemp;
  wimeml += itemp;

  if (!setXt(lu->nzcol, ncol, nrow, wiszel, wimeml))
  {
    printf("\n\n system error, pque setup.\n");
    ShutDown();
    exit(0);
  }

  colmax = wrmem;
  rbn1 = colmax + ncol;
#ifdef TIME_COUNT
  printf("lu_factor-A: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif
  if (step < nrow && lresp == lu_ok)
  {
    for (i = 0; i < nrow; ++i)
      if (lu->rowsta[i] == lu_rinit)
        putXt(lu->nzrow, i, lu->uuse[i]);

    for (j = 0; j < ncol; ++j)
      if (nnzaj[j] != nrow + 1)
        putXt(lu->nzcol, j, nnzaj[j]);

    lu_buildcol(lu, false, colmax);

    zmap = nnzaj;

    iZero(ncol, zmap, NULL);
    iSet(ncol, 1, cmiu, NULL);
#ifdef TIME_COUNT
    printf("lu_factor-B1: %d ms\n", GetTime() - Tm);
    Tm = GetTime();
#endif

    for (; step < nrow && lresp == lu_ok; ++step)
    {
      lu_choosepiv(lu, step, colmax, cmiu, &pivi, &pivj);
#ifdef TIME_COUNT
      CntA += GetTime() - Tm;
      Tm = GetTime();
#endif
      if (pivi == nili)
        break;

      lresp = lu_dopivot(lu, step, colmax, cmiu, pivi, pivj, zmap, rbn1); // HOTSPOT
#ifdef TIME_COUNT
      CntB += GetTime() - Tm;
      Tm = GetTime();
      CntPivot++;
#endif
    }
#ifdef TIME_COUNT
    printf("lu_factor-B2: choosepiv = %d ms, dopivot = %d ms, CntPivot = %d\n", CntA, CntB, CntPivot);
    Tm = GetTime();
#endif
    lu->rank = step;

    if (step < lu->nrow)
    {
      for (i = 0; i < nrow; ++i)
      {
        if (lu->rowsta[i] == lu_rinit)
        {
          lu->p[step++] = i;
          lu->rowsta[i] = lu_rfin;
          lu->uuse[i] = 0;
        }
      }
    }
  }
  else
    lu->rank = nrow;

  if (lresp == lu_ok)
  {
    for (i = 0; i < nrow; ++i)
    {
      lu->invp[lu->p[i]] = i;

      if (lu->uuse[i] > 0)
        lu->invq[lu->subj[lu->ufir[i]]] = i;
    }

    lu_buildcol(lu, true, NULL);

    lu_initupdate(lu);
  }

  ++lu->info.tnnfa;
#ifdef TIME_COUNT
  printf("lu_factor-C: %d ms\n", GetTime() - Tm);
  Tm = GetTime();
#endif

  return (lresp);
} /* lu_factor */

int lu_getrank(gsdec *lf,
               int *rank)
{
  *rank = lf->rank;
  return (lu_ok);
} /* lu_getrank */

int lu_getucond(gsdec *lf,
                double *abs_min_uii,
                double *abs_max_uii)
{
  int i;
  double rtemp;

  *abs_min_uii = MaxPositive;
  *abs_max_uii = 0.0;
  for (i = 0; i < lf->nrow; ++i)
  {
    if (lf->uuse[i])
      rtemp = fabs(lf->val[lf->ufir[i]]);
    else
      rtemp = 0.0;

    *abs_min_uii = min(*abs_min_uii, rtemp);
    *abs_max_uii = max(*abs_max_uii, rtemp);
  }

  return (lu_ok);
} /* lf_getucond */

int lu_getuii(gsdec *lf,
              int i,
              int *j,
              double *uii)
{
  if (lf->uuse[i])
  {
    *j = lf->subj[lf->ufir[i]];
    *uii = lf->val[lf->ufir[i]];
  }
  else
  {
    *j = lf->ncol;
    *uii = 0.0;
  }

  return (lu_ok);
} /* lf_getuii */
