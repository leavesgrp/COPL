#include "LPdefine.h"

static int ParseStr(char *str, char **key, char **val)
{
  int i, slen, klen, vlen;

  slen = strlen(str);

  for (i = 0, klen = 0; i < slen; i++)
  {
    if (isgraph(str[i]) && !klen)
      *key = str + i;

    if (!isalpha(str[i]))
    {
      if (str[i] == '_')
      {
        klen++;
        continue;
      }
      else if (klen)
        break;
    }
    else
    {
      str[i] = tolower(str[i]);
      klen++;
    }
  }

  if (!klen || i >= slen)
    return false;

  for (i = slen - 1, vlen = 0; i > 0; i--)
  {
    if (!isalnum(str[i]))
    {
      if (str[i] == '.' ||
          str[i] == '+' ||
          str[i] == '-')
        vlen++;
      else if (vlen)
      {
        *val = str + (i + 1);
        break;
      }
    }
    else
      vlen++;
  }

  if (!vlen || i <= 0)
    *val = NULL;

  if (str[slen - 1] == '\n')
    str[slen - 1] = '\0';
  return true;

} /* ParseStr */

int ReadParms(optpar *par,
              optmod *opt)
{
  FILE *fp;
  char set[][30] = {
      "bound_set_name", "cache_memory_size",
      "complementarity_tolerance", "cross_over_switch",
      "dense_column_threshold", "dual_feasible_tolerance",
      "infeasible_tolerance", "iteration_limit",
      "lower_bound_default", "nonzero_tolerance",
      "number_of_columns", "number_of_elements",
      "objective_set_name", "optimization_type",
      "pcg_starting_tolerance", "pivoting_tolerance",
      "presolve_level", "presolve_loop",
      "primal_dual_gap_tolerance", "primal_feasible_tolerance",
      "range_set_name", "right_hand_side_set_name",
      "save_solution", "step_factor",
      "upper_bound_default", ""},
       ss[LineSize], *key, *val;
  specid k;
  int i, size, ncard = 0;

  par->mcol = 100000;
  par->melm = 500000;

  par->iterlmt = 500;
  par->cachesz = 256;

  par->save = false;
  par->prlev = 5;
  par->prlop = 10;
  par->cross = false;

  par->bemin = true;
  par->idden = 0;

  par->rhsset = NULL;
  par->rngset = NULL;
  par->bndset = NULL;

  par->lower = 0.0;
  par->upper = MaxPositive;

  par->aijtol = 1.0E-12;
  par->tolpiv = 1.0E-14;
  par->tolgap = 1.0E-08;

  par->tolp = 1.0E-08;
  par->told = 1.0E-08;

  par->tolx = 1.0E-13;
  par->tols = 1.0E-08;

  par->tolmu = 1.0E-12;
  par->tolinf = 1.0E-08;
  par->pcgtol = 1.0E-10;

  par->Gamma = 0.99995;
  par->beta = 1.0E-08;
  par->bplus = par->Gamma * MaxPositive;

  opt->nub = 0;
  opt->nfx = 0;
  opt->nfr = 0;
  opt->nsk = 0;
  opt->nrg = 0;

  if (nspc)
  {
    strcpy(ss, nspc);
    fp = fopen(ss, "r");
    if (!fp)
    {
      sprintf(ss, SFMT ".spc", nspc);
      fp = fopen(ss, "r");
      if (!fp)
        return (ErrorProc(NOT_SCHFIL, nspc));
      sprintf(ss, SFMT, nspc);
    }
  }
  else
  {
    sprintf(ss, "COPLLP.spc");
    fp = fopen(ss, "r");
  }

  if (fp)
  {
    printf(" BEGIN SPC file ...\n");
    printf("   file name                   : " SFMT "\n", ss);
    fprintf(fout, "SPC FILE\n");
    fprintf(fout, "--------\n");

    while (!feof(fp))
    {

      ss[0] = '\0';
      fgets(ss, LineSize, fp);
      ncard++;
      fprintf(fout, I6FT "" SPC4, ncard);
      fputs(ss, fout);
      if (*ss == '*' || *ss == '\n' || *ss == '\0')
        continue;

      if (!ParseStr(ss, &key, &val))
        continue;

      for (k = BOUN; k < PEND; k++)
        if (strstr(key, set[k]))
          break;

      switch (k)
      {
      case BOUN:
        size = strlen(val);
        for (i = 0; i < size; i++)
          val[i] = toupper(val[i]);

        if (strncmp(val, "NULL", 4))
        {
          par->bndset = cAlloc(size + 1,
                               "for par->bndset in ReadParms");
          strcpy(par->bndset, val);
          printf("   bound set name              : " SFMT "\n", val);
        }
        break;
      case CACH:
        par->cachesz = atoi(val);
        if (par->cachesz < 0)
          par->cachesz = -par->cachesz;
        printf("   cache memory size           : " IFMT "\n", par->cachesz);
        break;
      case CGAP:
        par->tolmu = atof(val);
        if (par->tolmu < 1.E-15)
          par->tolmu = 1.0E-15;
        printf("   complementarity tolerance   : "
               "%G\n",
               par->tolmu);
        break;
      case CROS:
        printf("   cross-over                  : ");
        size = strlen(val);
        for (i = 0; i < size; i++)
          val[i] = toupper(val[i]);
        if (!strncmp(val, "ON", 2))
        {
          par->cross = 1;
          printf(SFMT "\n", val);
        }
        else if (!strncmp(val, "OFF", 3))
        {
          par->cross = 0;
          printf(SFMT "\n", val);
        }
        else
          fprintf(fout, " this setting ignored\n");
        break;
      case DENS:
        par->idden = atoi(val);
        if (par->idden < 1)
          par->idden = 1;
        printf("   dense column threshold      : " IFMT "\n", par->idden);
        break;
      case DINF:
        par->told = atof(val);
        if (par->told < 1.E-15)
          par->told = 1.0E-15;
        printf("   dual feasibility tolerance  : "
               "%G\n",
               par->told);
        break;
      case INFE:
        par->tolinf = atof(val);
        if (par->tolinf < 1.E-15)
          par->tolinf = 1.0E-15;
        printf("   infeasibility tolerance     : "
               "%G\n",
               par->tolinf);
        break;
      case ITER:
        par->iterlmt = atoi(val);
        if (par->iterlmt < 0)
          par->iterlmt = -par->iterlmt;
        printf("   iteration limit             : " IFMT "\n", par->iterlmt);
        break;
      case LOWE:
        par->lower = atof(val);
        printf("   defualt for lower bound     : "
               "%G\n",
               par->lower);
        break;
      case NONZ:
        par->aijtol = atof(val);
        if (par->aijtol < 0)
          par->aijtol = -par->aijtol;
        printf("   nonzero tolerance           : "
               "%G\n",
               par->aijtol);
        break;
      case NUMC:
        par->mcol = atoi(val);
        if (par->mcol < 0)
          par->mcol = -par->mcol;
        printf("   maximum number of columns   : " IFMT "\n", par->mcol);
        break;
      case NUME:
        par->melm = atoi(val);
        if (par->melm < 0)
          par->mcol = -par->melm;
        printf("   maximum number of elements  : " IFMT "\n", par->melm);
        break;
      case OBJS:
        size = strlen(val);
        for (i = 0; i < size; i++)
          val[i] = toupper(val[i]);

        if (strncmp(val, "NULL", 4))
        {
          par->objname = cAlloc(size + 1,
                                "for par->objname in ReadParms");
          strcpy(par->objname, val);
          printf("   objective row name          : " SFMT "\n", par->objname);
        }
        break;
      case OPTI:
        size = strlen(val);
        for (i = 0; i < size; i++)
          val[i] = toupper(val[i]);

        if (!strncmp(val, "MAX", 3))
        {
          printf("   optimization type           : maximize\n");
          par->bemin = false;
        }
        else if (!strncmp(val, "MIN", 3))
        {
          printf("   optimization type           : minimize\n");
          par->bemin = true;
        }
        else
        {
          printf("   the setting is ignored\n");
          par->bemin = true;
        }
        break;
      case PCGT:
        par->pcgtol = atof(val);
        if (par->pcgtol < 1.E-15)
          par->pcgtol = 1.0E-15;
        printf("   pcg starting tolerance      : "
               "%G\n",
               par->pcgtol);
        break;
      case PIVT:
        par->tolpiv = atof(val);
        if (par->tolpiv < 1.E-15)
          par->tolpiv = 1.0E-15;
        printf("   cholesky pivoting tolerance : "
               "%G\n",
               par->tolpiv);
        break;
      case PRES:
        par->prlev = atoi(val);
        if (par->prlev < 0)
          par->prlev = 0;
        if (par->prlev > 5)
          par->prlev = 5;
        printf("   presolving level            : " IFMT "\n", par->prlev);
        break;
      case PLOP:
        par->prlop = atoi(val);
        if (par->prlop <= 0)
          par->prlop = 1;
        printf("   presolving loop             : " IFMT "\n", par->prlop);
        break;
      case PGAP:
        par->tolgap = atof(val);
        if (par->tolgap < 1.E-15)
          par->tolgap = 1.0E-15;
        printf("   primal dual gap tolerance   : "
               "%G\n",
               par->tolgap);
        break;
      case PINF:
        par->tolp = atof(val);
        if (par->tolp < 1.E-15)
          par->tolp = 1.0E-15;
        printf("   primal feasibility tolerance: "
               "%G\n",
               par->tolp);
        break;
      case RANG:
        size = strlen(val);
        for (i = 0; i < size; i++)
          val[i] = toupper(val[i]);

        if (strncmp(val, "NULL", 4))
        {
          par->rngset = cAlloc(size + 1,
                               "par->rngset in ReadParms");
          strcpy(par->rngset, val);
          printf("   range set name              : " SFMT "\n", par->rngset);
        }
        break;
      case RHSS:
        size = strlen(val);
        for (i = 0; i < size; i++)
          val[i] = toupper(val[i]);

        if (strncmp(val, "NULL", 4))
        {
          par->rhsset = cAlloc(size + 1,
                               "par->rhsset in ReadParms");
          strcpy(par->rhsset, val);
          printf("   rhs set name                : " SFMT "\n", par->rhsset);
        }
        break;
      case SAVE:
        printf("   save solution                 : ");
        size = strlen(val);
        for (i = 0; i < size; i++)
          val[i] = toupper(val[i]);
        if (!strncmp(val, "ON", 2))
        {
          par->save = 1;
          printf(SFMT "\n", val);
        }
        else if (!strncmp(val, "OFF", 3))
        {
          par->save = 0;
          printf(SFMT "\n", val);
        }
        else
          fprintf(fout, " this setting ignored\n");
        break;
      case STPF:
        par->Gamma = atof(val);
        if (par->Gamma > 0.99995)
          par->Gamma = 0.99995;
        if (par->Gamma < 0.1)
          par->Gamma = 0.1;
        printf("   step size factor            : "
               "%G\n",
               par->Gamma);
        break;
      case UPPE:
        par->upper = atof(val);
        printf("   default for upper bound     : "
               "%G\n",
               par->upper);
        break;
      default:
        fprintf(fout, "   no such setting\n");
        break;
      }
    }

    FileClose(&fp);
    printf(" END SPC file\n\n");
    fprintf(fout, "\n");
  }

  par->cachesz *= 1000;
  return true;
} /* ReadParms */
