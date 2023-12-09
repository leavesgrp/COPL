#include "LPcross.h"

static int lindep_crash(double tolapiv,
                        int m,
                        int n,
                        int *mc,
                        int *nc,
                        int inew[],
                        int iinv[],
                        int jnew[],
                        int jinv[],
                        matrix *a,
                        matrix *at,
                        pivbas *basis,
                        int *nol,
                        int lbuf[],
                        int wisze,
                        int wimem[],
                        int wrsze,
                        double wrmem[])
{
  int i, j, k, t, sze, pos, ncrash, mc0, cur,
      *nnzi, *nnzj, *list, *fir, *link,
      sresp = ProcOk;
  array *ai, *aj;

  if (wisze < 1 + 2 * m + 3 * n)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  nnzi = wimem;
  nnzj = nnzi + m;
  list = nnzj + n;
  fir = list + n;
  link = fir + (n + 1);

  /*
   * Count the numbers of non-zeros in each row and column.
   */

  iZero(n, nnzj, NULL);

  for (k = 0; k < *mc; ++k)
  {
    i = iinv[k];
    nnzi[i] = 0;
    ai = a->ia + i;
    for (t = 0; t < ai->nn0; ++t)
    {
      j = ai->ja[t];
      if (jnew[j] < *nc)
      {
        ++nnzi[i];
        ++nnzj[j];
      }
    }
  }

  for (k = *mc; k < m; ++k)
    nnzi[iinv[k]] = n + 1;

  for (k = 0; k < *mc;)
  {
    i = iinv[k];
    if (nnzi[i] == 0)
    {
      lbuf[*nol] = i;
      ++*nol;

      --*mc;
      iSwap(k, *mc, iinv);

      inew[iinv[k]] = k;
      inew[iinv[*mc]] = *mc;
    }
    else
      ++k;
  }

  sze = 0;
  for (k = 0; k < *nc; ++k)
  {
    j = jinv[k];
    if (nnzj[j] == 1)
      list[sze++] = j;
  }

  fSort(m, n, nnzi, fir, link);

  ncrash = 0;
  mc0 = *mc;
  cur = n;
  do
  {
    while (sze)
    {
      --sze;

      j = list[sze];
      if (nnzj[j] == 1)
      {
        ++ncrash;

        aj = at->ia + j;

        for (t = 0; t < aj->nn0; ++t)
          if (inew[aj->ja[t]] < *mc)
            break;

        if (t == aj->nn0)
        {
          printf("\n\n system error, not found.\n");
          ShutDown();
          exit(0);
        }

        if (fabs(aj->an[t]) > tolapiv)
        {
          i = aj->ja[t];
          k = inew[i];

          if (basis[i].vt != PIV_ART)
          {
            printf("\n\n system error, not art var.\n");
            ShutDown();
            exit(0);
          }

          basis[i].vt = PIV_COL;
          basis[i].j = j;

          --*mc;

          iSwap(k, *mc, iinv);

          inew[iinv[k]] = k;
          inew[iinv[*mc]] = *mc;

          --*nc;
          pos = jnew[j];

          iSwap(pos, *nc, jinv);

          jnew[jinv[pos]] = pos;
          jnew[jinv[*nc]] = *nc;
          ai = a->ia + i;

          for (t = 0; t < ai->nn0; ++t)
          {
            j = ai->ja[t];
            if (nnzj[j])
            {
              --nnzj[j];

              if (nnzj[j] == 1)
                list[sze++] = j;
            }
          }
        }
      }
    }

    for (; cur;)
    {
      i = fir[cur];
      if (i != m)
      {
        k = inew[i];

        if (k < *mc)
        {
          --*mc;

          iSwap(k, *mc, iinv);

          inew[iinv[k]] = k;
          inew[iinv[*mc]] = *mc;

          ai = a->ia + i;

          for (t = 0; t < ai->nn0; ++t)
          {
            j = ai->ja[t];
            if (nnzj[j])
            {
              --nnzj[j];

              if (nnzj[j] == 1)
                list[sze++] = j;
            }
          }
          break;
        }

        fir[cur] = link[i];
      }
      else
        --cur;
    }
  } while (cur);

  *mc = mc0;

  return (sresp);
} /* lindep_crash */

static int lindep_crash2(double tolapiv,
                         int m,
                         int n,
                         int mc,
                         int *nc,
                         int inew[],
                         int iinv[],
                         int jnew[],
                         int jinv[],
                         matrix *a,
                         matrix *at,
                         pivbas *basis,
                         gsdec *lu,
                         int wisze,
                         int wimem[],
                         int wrsze,
                         double wrmem[])
{
  int i, j, k, t, try_, pos, minnz, minj, rsze,
      *nnzj, *bas, *wimeml, lresp;
  double loc_tolapiv = 1.0e-3, temp_tolapiv, uii;
  int sresp = ProcOk;
  array *ai;

  if (wisze < m + n)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  nnzj = wimem;
  wimeml = wimem + n;

  /*
   * Count the numbers of non-zeros in each row and column.
   */

  iZero(n, nnzj, NULL);

  for (k = 0; k < mc; ++k)
  {
    i = iinv[k];
    ai = a->ia + i;
    for (t = 0; t < ai->nn0; ++t)
    {
      j = ai->ja[t];
      if (jinv[j] < *nc)
        ++nnzj[j];
    }
  }

  for (try_ = 0; try_ < 3 && sresp == ProcOk; ++try_)
  {
    rsze = 0;
    for (k = 0; k < mc; ++k)
    {
      if (basis[iinv[k]].vt == PIV_ART)
        ++rsze;
    }

    if (rsze <= 200)
      break;

    pivlp_mvartrow(mc, iinv, basis);

    for (k = 0; k < mc; ++k)
    {
      i = iinv[k];

      if (basis[i].vt == PIV_ART)
      {
        ai = a->ia + i;
        minnz = m;

        for (t = 0; t < ai->nn0; ++t)
        {
          j = ai->ja[t];
          if (jnew[j] < *nc &&
              nnzj[j] < minnz &&
              fabs(ai->an[t]) >= 1.0e-3)
          {
            minj = j;
            minnz = nnzj[j];
          }
        }

        if (minnz < m)
        {
          j = minj;
          basis[i].vt = PIV_COL;
          basis[i].j = j;

          --*nc;
          pos = jnew[j];

          iSwap(*nc, pos, jinv);

          jnew[jinv[pos]] = pos;
          jnew[jinv[*nc]] = *nc;
        }
      }
    }

    temp_tolapiv = lu->tolapiv;
    lu->tolapiv = loc_tolapiv;
    lu->dropl = true;
    lu->dropu = true;

    lresp = pivlp_factor1(m, mc,
                          iinv, inew,
                          at,
                          false,
                          basis, lu,
                          wisze, wimem,
                          wrsze, wrmem);

    lu->tolapiv = temp_tolapiv;
    lu->dropl = false;
    lu->dropu = false;

    switch (lresp)
    {
    case lu_ok:
    case lu_insta:
      bas = wimeml;

      /*
       * Store old basis.
       */

      for (k = 0; k < mc; ++k)
      {
        i = iinv[k];
        if (basis[i].vt == PIV_ART)
          bas[i] = m + n;
        else
        {
          bas[i] = basis[i].j + m;

          j = basis[i].j;

          pos = jnew[j];

          if (pos < *nc)
            exit(0);

          iSwap(*nc, pos, jinv);

          jnew[jinv[pos]] = pos;
          jnew[jinv[*nc]] = *nc;
          ++*nc;
        }
      }

      for (k = 0; k < mc; ++k)
      {
        lu_getuii(lu, k, &i, &uii);

        if (i == mc)
        {
          basis[iinv[k]].vt = PIV_ART;
          basis[iinv[k]].j = iinv[k];
        }
        else
        {
          if (bas[iinv[i]] < m)
          {
            basis[iinv[k]].vt = PIV_ROW;
            basis[iinv[k]].j = bas[iinv[i]];
          }
          else if (bas[iinv[i]] < m + n)
          {
            basis[iinv[k]].vt = PIV_COL;
            basis[iinv[k]].j = bas[iinv[i]] - m;

            j = bas[iinv[i]] - m;

            --*nc;
            pos = jnew[j];

            iSwap(pos, *nc, jinv);

            jnew[jinv[pos]] = pos;
            jnew[jinv[*nc]] = *nc;
          }
        }
      }

      break;
    default:
      sresp = FacError;
    }
  }

  return (sresp);
} /* lindep_crash2 */

static void lindep_reorder(int m,
                           int n,
                           int mc,
                           int nc,
                           int inew[],
                           int iinv[],
                           int jnew[],
                           int jinv[],
                           matrix *a,
                           matrix *at,
                           pivbas *basis,
                           int wisze,
                           int wimem[])
{
  int i, j, k, s, cur,
      *nnzi, *fir, *link, *new_iinv;
  array *ai;

  if (wisze < 1 + 1 * m + 2 * n)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  nnzi = wimem;
  fir = nnzi + m;
  link = fir + (n + 1);

  /*
   * Count the numbers of non-zeros in each row and column.
   */
  pivlp_mvartrow(mc, iinv, basis);

  for (k = 0; k < mc; ++k)
  {
    i = iinv[k];
    if (basis[i].vt == PIV_ART)
    {
      nnzi[i] = 0;
      ai = a->ia + i;
      for (s = 0; s < ai->nn0; ++s)
      {
        j = ai->ja[s];
        if (jnew[j] < nc)
          ++nnzi[i];
      }
    }
    else
      nnzi[i] = n + 1;
  }

  for (k = mc; k < m; ++k)
    nnzi[iinv[k]] = n + 1;

  fSort(m, n, nnzi, fir, link);

  new_iinv = nnzi;

  s = 0;
  for (k = 0; k < mc; ++k)
  {
    i = iinv[k];
    if (basis[i].vt != PIV_ART)
      new_iinv[s++] = i;
  }

  for (cur = 0; cur <= n;)
  {
    i = fir[cur];
    if (i != m)
    {
      new_iinv[s++] = i;

      fir[cur] = link[i];
    }
    else
      ++cur;
  }

  if (s != mc)
    exit(0);

  for (k = 0; k < mc; ++k)
  {
    iinv[k] = new_iinv[k];
    inew[iinv[k]] = k;
  }
} /* lindep_reorder */

static void lindep_compress(int m,
                            int n,
                            int mc,
                            int nc,
                            matrix *a,
                            int iinv[],
                            int jnew[],
                            pivbas *basis,
                            int rsze,
                            int rset[],
                            int offset[],
                            int wzisze,
                            int wzimem[])
{
  int k, i, s, t, j,
      *finj;
  double rtemp;
  array *ai;

  if (wzisze < n)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  finj = wzimem;

  for (k = 0; k < mc; ++k)
  {
    i = iinv[k];
    j = basis[i].j;
    if (basis[i].vt == PIV_COL)
      finj[j] = 1;
  }

  for (k = 0; k < rsze; ++k)
  {
    i = rset[k];
    j = basis[i].j;
    if (basis[i].vt == PIV_COL)
      finj[j] = 0;
  }

  for (k = 0; k < mc; ++k)
  {
    i = iinv[k];
    ai = a->ia + i;
    s = 0;
    for (t = 0; t < ai->nn0; ++t)
    {
      j = ai->ja[t];
      if (jnew[j] >= nc && finj[j])
      {
        rtemp = ai->an[t];

        ai->an[t] = ai->an[s];
        ai->an[s] = rtemp;

        ai->ja[t] = ai->ja[s];
        ai->ja[s] = j;
        ++s;
      }
    }
    offset[i] = s;
  }

  /*
   * Reset the buffer.
   */

  for (k = 0; k < mc; ++k)
  {
    i = iinv[k];
    j = basis[i].j;
    if (basis[i].vt == PIV_COL)
      finj[j] = 0;
  }
} /* lindep_compress */

static int lindep_pivot(double tolapiv,
                        int m,
                        int n,
                        int *mc,
                        int *nc,
                        int iinv[],
                        int inew[],
                        int jinv[],
                        int jnew[],
                        int rsze,
                        int rset[],
                        matrix *a,
                        matrix *at,
                        pivbas *basis,
                        int *nol,
                        int lbuf[],
                        gsdec *lu,
                        int iszero,
                        int wisze,
                        int wimem[],
                        int wzisze,
                        int wzimem[],
                        int wzrsze,
                        double wzrmem[])
{
  int refac;
  int i, j, t, r, szeu, k, l, q, inc, itemp, maxj, blksze,
      nnz1, nnz2, nnz, nilj,
      *offset, *sub1, *sub2, *subu,
      *map1, *map2, *mapu, *ajsub, *arsub;
  int lresp = lu_ok;
  double maxa, rtemp,
      *val1, *val2, *valu, *ajval, *arval;
  int resp = ProcOk;
  array *ar, *aj;

  if (wisze < 2 + 5 * m + max(6 * m, 1 * n) ||
      wzisze < n || wzrsze < 2 * m + 1 * n)
  {
    printf("\n\n work space error.\n");
    ShutDown();
    exit(0);
  }

  /*
   * Storage layout is important.
   */

  if (!rsze)
    return (resp);

  map1 = wimem;
  map2 = map1 + m;

  offset = map2 + m;

  sub1 = offset + m;
  sub2 = sub1 + m;
  subu = sub2 + m;

  val1 = wzrmem;
  val2 = val1 + m;
  valu = val2 + m;

  mapu = wzimem;

  nilj = n;

  if (iszero)
  {
    /*
     * Check if zero.
     */
  }
  else
  {
    iZero(*mc, map1, NULL);
    iZero(*mc, map2, NULL);
  }

  lindep_compress(m, n,
                  *mc, *nc,
                  a,
                  iinv,
                  jnew,
                  basis,
                  rsze, rset,
                  offset,
                  n, mapu);

  for (t = rsze; t && resp == ProcOk;)
  {
    t--;
    i = rset[t];
    j = basis[i].j;

    itemp = inew[i];
    nnz2 = 1;
    val2[itemp] = 1.0;
    sub2[0] = itemp;
    map2[itemp] = 1;

    lu_btran2(lu, nnz2, val2, sub2, map2, &nnz1, val1, sub1, map1);

    if (nnz1 == 1)
    {
      r = sub1[0];
      rtemp = val1[r];
      r = iinv[r];

      ar = a->ia + r;

      nnz = ar->nn0;
      arsub = ar->ja;
      arval = ar->an;

      szeu = 0;
      if (rtemp == 1.0)
      {
        for (l = offset[r]; l < nnz; ++l)
        {
          valu[arsub[l]] = arval[l];
          subu[szeu++] = arsub[l];
        }
      }
      else
      {
        for (l = offset[r]; l < nnz; ++l)
        {
          valu[arsub[l]] = rtemp * arval[l];
          subu[szeu++] = arsub[l];
        }
      }
    }
    else
    {
      szeu = 0;
      for (k = 0; k < nnz1; ++k)
      {
        r = sub1[k];

        rtemp = val1[r];
        r = iinv[r];

        ar = a->ia + r;

        nnz = ar->nn0;
        arsub = ar->ja;
        arval = ar->an;

        for (l = offset[r]; l < nnz; ++l)
        {
          q = arsub[l];
          valu[q] += rtemp * arval[l];
          if (!mapu[q])
          {
            mapu[q] = 1;
            subu[szeu++] = q;
          }
        }
      }

      /*
       * mapu is not needed any more.
       */

      iZero(szeu, mapu, subu);
    }

    dZero(nnz1, val1, sub1);
    iZero(nnz1, map1, sub1);

    /*
     * Find pivot element.
     */

    maxj = nilj;
    maxa = tolapiv;
    for (k = 0; k < szeu; ++k)
    {
      q = subu[k];
      if (jnew[q] < *nc)
      {
        if ((maxj == nilj && fabs(valu[q]) > tolapiv) || fabs(valu[q]) > maxa)
        {
          maxj = q;
          maxa = fabs(valu[q]);
        }
      }
      valu[q] = 0.0;
    }

    if (maxj == nilj)
    {
      lbuf[*nol] = j;
      ++*nol;
    }
    else
    {
      /*
       * Check for numerical stability of pivot.
       */

      inc = maxj;

      /*
       * Update the basis record.
       */

      basis[i].vt = PIV_COL;
      basis[i].j = inc;

      --*nc;
      itemp = jnew[inc];

      if (jinv[itemp] != inc)
        exit(0);

      iSwap(*nc, itemp, jinv);

      jnew[jinv[*nc]] = *nc;
      jnew[jinv[itemp]] = itemp;

      /*
       * Unpack the incoming column.
       */

      aj = at->ia + inc;

      nnz = aj->nn0;
      ajsub = aj->ja;
      ajval = aj->an;

      nnz1 = 0;
      for (k = 0; k < nnz; ++k)
      {
        r = inew[ajsub[k]];
        if (r < *mc)
        {
          val1[r] = ajval[k];
          sub1[nnz1] = r;
          map1[r] = 1;
          ++nnz1;
        }
      }

      lresp = lu_supd(lu, inew[i], inew[i],
                      nnz1, val1, sub1, map1,
                      sub2, subu);

      if (lresp == lu_refac)
      {
        /*
         * t: Is number of unfinished rows.
         */

        blksze = basord_blklwr(m, n, *mc,
                               inew, iinv,
                               t, rset,
                               a, at,
                               basis,
                               wisze, wimem,
                               wzisze, wzimem);
        *mc -= blksze;
      }

      resp = pivlp_prolresp(10,
                            m, *mc, n,
                            iinv, inew,
                            at,
                            basis,
                            NULL, NULL,
                            NULL, NULL, NULL,
                            NULL, false,
                            lu, lresp,
                            &refac,
                            min(wisze, 2 + 11 * (*mc)), wimem,
                            min(wzrsze, 2 * (*mc)), wzrmem);

      if (resp != ProcOk || refac)
      {
        dZero(min(wzrsze, 2 * (*mc)), wzrmem, NULL);

        iZero(*mc, map1, NULL);
        iZero(*mc, map2, NULL);
      }

      if (refac)
      {
        /*
         * Note: t is the current rsze.
         */

        lindep_compress(m, n,
                        *mc, *nc,
                        a,
                        iinv,
                        jnew,
                        basis,
                        rsze, rset,
                        offset,
                        n, mapu);
      }
    }
  }

  return (resp);
} /* lindep_pivot */

static int DetectLdpRow(double tolapiv,
                        int m,
                        int n,
                        int mc,
                        int nc,
                        int inew[],
                        int iinv[],
                        int jnew[],
                        int jinv[],
                        matrix *a,
                        matrix *at,
                        pivbas *basis,
                        gsdec *lu,
                        int *nol,
                        int lbuf[],
                        int wisze,
                        int wimem[],
                        int wrsze,
                        double wrmem[])
{
  int i, k, wiszel, wziszel, wzrszel, rank, rsze, blksze, mc0,
      *rset,
      *wimeml, *wzimeml;
  int lresp;
  double *wzrmeml;
  int sresp = ProcOk;

  if (wisze < 2 + 6 * m + 1 * n + max(6 * m, 1 * n) || wrsze < 2 * m + 1 * n)
    printf("\n\n work space error.\n");

  for (k = 0; k < mc; ++k)
  {
    basis[iinv[k]].vt = PIV_ART;
    basis[iinv[k]].j = iinv[k];
  }

  sresp = lindep_crash(tolapiv,
                       m, n,
                       &mc, &nc,
                       inew, iinv,
                       jnew, jinv,
                       a, at,
                       basis,
                       nol, lbuf,
                       wisze, wimem,
                       wrsze, wrmem);

  if (sresp == ProcOk)
  {
    sresp = lindep_crash2(tolapiv,
                          m, n,
                          mc, &nc,
                          inew, iinv,
                          jnew, jinv,
                          a, at,
                          basis, lu,
                          wisze, wimem,
                          wrsze, wrmem);

    if (sresp == ProcOk)
    {
      wziszel = n;
      wiszel = wisze - wziszel - m;

      rset = wimem;
      wzimeml = rset + m;
      wimeml = wzimeml + wziszel;

      wzrszel = min(2 * m + 1 * n, wrsze);
      wzrmeml = wrmem;

      iZero(wziszel, wzimeml, NULL);
      dZero(wzrszel, wzrmeml, NULL);

      rsze = 0;
      for (k = 0; k < mc; ++k)
      {
        i = iinv[k];
        if (basis[i].vt == PIV_ART)
          rset[rsze++] = i;
      }

      blksze = basord_blklwr(m, n, mc,
                             inew, iinv,
                             rsze, rset,
                             a, at,
                             basis,
                             wiszel, wimeml,
                             wziszel, wzimeml);

      mc -= blksze;

      lindep_reorder(m, n,
                     mc, nc,
                     inew, iinv,
                     jnew, jinv,
                     a, at,
                     basis,
                     wiszel, wimeml);

      rsze = 0;
      for (k = 0; k < mc && sresp == ProcOk; ++k)
      {
        i = iinv[k];
        if (basis[i].vt == PIV_ART)
          rset[rsze++] = i;

        if (rsze >= mc || k + 1 == mc)
        {
          mc0 = mc;
          mc = k + 1;

          blksze = basord_blklwr(m, n, mc,
                                 inew, iinv,
                                 rsze, rset,
                                 a, at,
                                 basis,
                                 wiszel, wimeml,
                                 wziszel, wzimeml);

          mc -= blksze;

          lresp = pivlp_factor1(m, mc,
                                iinv, inew,
                                at,
                                true,
                                basis,
                                lu,
                                wiszel, wimeml,
                                wzrszel, wzrmeml);

          dZero(wzrszel, wzrmeml, NULL);

          lu_getrank(lu, &rank);
          if (lresp == lu_ok && rank == mc)
          {
            sresp = lindep_pivot(tolapiv,
                                 m, n,
                                 &mc, &nc,
                                 iinv, inew,
                                 jinv, jnew,
                                 rsze, rset,
                                 a, at,
                                 basis,
                                 nol, lbuf,
                                 lu,
                                 false,
                                 wiszel, wimeml,
                                 wziszel, wzimeml,
                                 wzrszel, wzrmeml);
          }
          else
            sresp = FacError;

          mc = mc0;
          rsze = 0;
        }
      }
    }
  }
  return (sresp);
} /* DetectLdpRow */

int LdpRowProc(double tolaij,
               double tolapiv,
               int m,
               int n,
               int mc,
               int nc,
               int inew[],
               int iinv[],
               int jnew[],
               int jinv[],
               matrix *a,
               matrix *at,
               double b[],
               pivbas *basis,
               double xb[],
               int *nol,
               int lbuf[],
               int wisze,
               int wimem[],
               int wrsze,
               double wrmem[])
{
  int k, rank;
  gsdec *lu;
  int lresp;
  int resp = ProcOk;

  *nol = 0;

  if (!mc)
    return (resp);

  lu = LufAlloc(m, m, 0, true, "lu, lindep_check1");

  resp = DetectLdpRow(tolapiv,
                      m, n,
                      mc, nc,
                      inew, iinv,
                      jnew, jinv,
                      a, at,
                      basis, lu,
                      nol, lbuf,
                      wisze, wimem,
                      wrsze, wrmem);

  pivlp_mvartrow(mc, iinv, basis);

  if (b && xb)
  {
    lresp = pivlp_factor1(m, mc,
                          iinv, inew,
                          at,
                          true,
                          basis,
                          lu,
                          wisze, wimem,
                          wrsze, wrmem);

    lu_getrank(lu, &rank);
    if (lresp == lu_ok && rank == mc)
    {

      for (k = 0; k < mc; ++k)
        xb[k] = b[iinv[k]];

      lu_ftran(lu, xb, b);

      for (k = 0; k < mc; ++k)
        xb[iinv[k]] = b[k];
    }
    else
      resp = FacError;
  }

  LufFree(&lu);

  return (resp);
} /* LdpRowProc */

void ChkRowLdpt(optpar *param,
                optmod *pdat,
                int wisze,
                int wimem[],
                int wrsze,
                double wrmem[])
{
  int i, j, k, t, wiszel, wrszel, m = pdat->m,
                                  n = pdat->n, mc = 0, nc = 0, nld = 0,
                                  *inew, *iinv, *jnew, *jinv, *lbuf, *wimeml;
  pivbas *basis;
  double *b, *xb, *wrmeml;
  int resp;

  if (wisze < 3 * m + 2 * n || wrsze < 2 * m)
  {
    printf("\n\n Out of work space.\n");
    ShutDown();
    exit(0);
  }

  inew = wimem;
  iinv = inew + m;
  jnew = iinv + m;
  jinv = jnew + n;
  lbuf = jinv + n;

  wiszel = wisze - (3 * m + 2 * n);
  wimeml = wimem + (3 * m + 2 * n);

  /*
   * include all linear rows
   */
  b = wrmem;
  xb = wrmem + m;

  wrszel = wrsze - 2 * m;
  wrmeml = xb + m;

  t = m;
  for (i = 0; i < m; ++i)
  {
    if (pdat->rflg[i] == 'T')
    {
      iinv[mc] = i;
      inew[i] = mc;
      b[i] = pdat->b[i];
      mc++;
    }
    else
    {
      t--;
      iinv[t] = i;
      inew[i] = t;
    }
  }

  t = n;
  for (j = 0; j < n; ++j)
  {
    if (pdat->cflg[j] == 'T')
    {
      jinv[nc] = j;
      jnew[j] = nc;
      nc++;
    }
    else
    {
      t--;
      jinv[t] = j;
      jnew[j] = t;
    }
  }

  basis = PbsAlloc(m, "basis, ChkRowLdpt");

  for (i = 0; i < m; ++i)
  {
    basis[i].vt = PIV_ART;
    basis[i].j = i;
  }

  resp = LdpRowProc(param->aijtol, param->tolpiv,
                    pdat->m, pdat->n,
                    mc, nc,
                    inew, iinv, jnew, jinv,
                    pdat->a, pdat->at,
                    b, basis, xb,
                    &nld, lbuf,
                    wiszel, wimeml,
                    wrszel, wrmeml);

  for (k = 0; k < nld; k++)
  {
    i = lbuf[k];
    if (basis[i].vt != PIV_ART || basis[i].j != i)
    {
      printf("\n\n internal error, basis error.\n");
      ShutDown();
      exit(0);
    }

    if (fabs(xb[i]) > param->aijtol)
    {
      printf("\n\n Exit -- 2: primal infeasible"
             " in row %d, rhs=%G!=0.0\n",
             i + 1, xb[i]);
      ShutDown();
      exit(0);
    }

    AddStack();
    Int2(RNUL, i);
    pdat->rflg[i] = 'F';
  }

  pdat->nrld = nld;

  PbsFree(&basis);
} /* ChkRowLdpt */
