#include "LPdefine.h"

#ifdef TEST
void CleanOptdat(optmod *opt, optpar *par, clock_t tim[])
{
  opt->ista = ProcOk;

  opt->m = 0;
  opt->n = 0;
  opt->nz = 0;

  opt->nlb = 0;
  opt->nub = 0;
  opt->nfr = 0;
  opt->nfx = 0;

  opt->nsk = 0;
  opt->nrg = 0;
  opt->ineq = 0;
  opt->nlow = 0;
  opt->nupp = 0;
  opt->nfix = 0;
  opt->nfre = 0;
  opt->nbmi = 0;
  opt->nbpl = 0;

  opt->rnamesz = 0;
  opt->cnamesz = 0;

  opt->obj = 0.0;

  opt->items = 0;
  opt->bytes = 0;

  opt->prtst = false;
  opt->nrld = 0;
  opt->nrnul = 0;
  opt->nrsg = 0;
  opt->nrdb = 0;
  opt->nrdm = 0;
  opt->ncfx = 0;
  opt->ncnul = 0;
  opt->nrfc = 0;
  opt->nrdu = 0;
  opt->ncdm = 0;
  opt->ncdu = 0;

  opt->iter = 0;

  opt->pobj = 0.0;
  opt->dobj = 0.0;
  opt->rgap = 0.0;
  opt->mu = 0.0;

  opt->inf0 = 0.0;
  opt->rho = 0.0;
  opt->tau0 = 0.0;
  opt->kap0 = 0.0;

  opt->cx = 0.0;
  opt->by = 0.0;
  opt->rg = 0.0;

  opt->rdf = 0.0;
  opt->dfdt = 0.0;
  opt->dfnm = 0.0;
  opt->rddt = 0.0;
  opt->ddtu = 0.0;
  opt->rdnm = 0.0;
  opt->rrdnm = 0.0;

  opt->rpdt = 0.0;
  opt->pdtu = 0.0;
  opt->pnmu = 0.0;
  opt->rpnm = 0.0;
  opt->rrpnm = 0.0;

  opt->rgdt = 0.0;
  opt->smdt = 0.0;
  opt->cgap = 0.0;

  opt->infe = 0.0;
  opt->nall = 0.0;
  opt->mu0 = 0.0;
}
#endif

clock_t GetTime(void)
{
  clock_t t;

#ifdef HPMACHINE
  struct tms
  {
    clock_t tms_utime;
    clock_t tms_stime;
    clock_t tms_cutime;
    clock_t tms_cstime;
  } tmrec;

  t = times(&tmrec);
  t = tmrec.tms_utime + tmrec.tms_stime;
#else
  t = clock();
#endif

  return t;
} /* GetTime */

int SetTime(clock_t tim[],
            int phase)
{
  int i;

  tim[phase] = GetTime();
  for (i = phase + 1; i < ELAPS; i++)
    tim[i] = tim[phase];

  return true;
}

int CheckArgs(int argc, char *argv[], char *ss)
{
  int i, id;

  id = 0;
  for (i = 1; i < argc; i++)
  {
    if (strstr(argv[i], ss))
    {
      argv[i] += 2;
      if (argv[i][0] == ':')
      {
        argv[i]++;
        id = i;
        break;
      }
      if (argv[i][0])
      {
        id = i;
        break;
      }
      id = i + 1;
      if (argv[id][0] == '-')
        id = 0;
      break;
    }
  }
  return id;
} /* Checkargs */

int FilesOpen(void)
{
  fout = fopen(nout, "w");
  if (!fout)
    return ErrorProc(NOT_DSKSPC, COPLFILE_BIN);

  return true;
} /* FilesOpen */

void FileClose(FILE **fp)
{
  if ((*fp))
  {
    fclose(*fp);
    *fp = NULL;
  }
}

void PrintHead(void)
{
  char Title[8][60] = {
      " ==========================================\n",
      " Computational Optimization Program Library\n",
      "                                           \n",
      "     (  C  O  P  L  )     Version 1.0      \n",
      "                                           \n",
      " ==========================================\n"};
  int i;

  printf("\n\n\n");
  for (i = 0; i < 8; i++)
    printf("%s", Title[i]);
  printf("\n *************************** COPL STARTS ");
  printf("***************************\n\n");

  return;
} /* PrintHead */

int FprintHead(void)
{
  if (!FilesOpen())
    return false;
  fprintf(fout, "====================\n");
  fprintf(fout, "C O P L  1.0\n");
  fprintf(fout, "====================\n\n");
  return true;
} /* FprintHead */

int GetFilen(char flag,
             char *ss)
{
  int k;

  k = 1;
  switch (flag)
  {
  case 'f':
    printf(" Please type the name of your data file: ");
    fflush(stdin);
    gets(ss);
    break;
  case 'o':
    sprintf(ss, SFMT, COPLFILE_OUT);
    break;

  default:
    k = 0;
    break;
  }
  return true;
}

int CompString(const char *s1,
               const char *s2)
{
  if (!s1 || !s2)
    return true;
  return strcmp(s1, s2);
} /* CompString */

int LeftDots(char *ss, int len)
{
  char sz[60];
  int i, j, k;

  k = strlen(ss);
  j = len - k;

  strcpy(sz, ss);
  for (i = 0; i < j; i++)
    ss[i] = '.';
  for (i = j; i < len; i++)
    ss[i] = sz[i - j];
  ss[len] = '\0';
  return true;
}

void PrintEnd(clock_t tim[])
{
  double tmproc, tmscal;
  char ss[StringSize];

#ifdef HPMACHINE
  tmscal = 0.01;
#else
  tmscal = 0.001;
#endif

  fprintf(fout, "\nTime distribution");
  fprintf(fout, "\n-----------------\n");

  tmproc = tmscal * (double)(tim[DATAIN] - tim[START]);
  sprintf(ss, " " TFMT " sec", tmproc);
  LeftDots(ss, 39);
  printf(" time for initial data input " SFMT "\n", ss);
  fprintf(fout, " time for initial data input " SFMT "\n", ss);

  tmproc = (double)(tim[PRSLV] - tim[DATAIN]) * tmscal;
#ifdef TEST
  fprintf(fres, "%8.2f&", tmproc);
#endif
  sprintf(ss, " " TFMT " sec", tmproc);
  LeftDots(ss, 41);
  printf(" time for presolve process " SFMT "\n", ss);
  fprintf(fout, " time for presolve process " SFMT "\n", ss);

  tmproc = (double)(tim[SYBOL] - tim[PRSLV]) * tmscal;
  sprintf(ss, " " TFMT " sec", tmproc);
  LeftDots(ss, 37);
  printf(" time for symbolic computation " SFMT "\n", ss);
  fprintf(fout, " time for symbolic computation " SFMT "\n", ss);

  tmproc = (double)(tim[NUMER] - tim[SYBOL]) * tmscal;
  sprintf(ss, " " TFMT " sec", tmproc);
  LeftDots(ss, 36);
  printf(" time for numerical computation " SFMT "\n", ss);
  fprintf(fout, " time for numerical computation " SFMT "\n", ss);

  tmproc = (double)(tim[CROSS] - tim[NUMER]) * tmscal;
#ifdef TEST
  fprintf(fres, "%8.2f&", tmproc);
#endif
  sprintf(ss, " " TFMT " sec", tmproc);
  LeftDots(ss, 35);
  printf(" time for cross-over computation " SFMT "\n", ss);
  fprintf(fout, " time for cross-over computation " SFMT "\n", ss);

  tmproc = (double)(tim[ELAPS] - tim[START]) * tmscal;
#ifdef TEST
  fprintf(fres, "%8.2f\\\\\n", tmproc);
#endif
  sprintf(ss, " " TFMT " sec", tmproc);
  LeftDots(ss, 40);
  printf(" time for the whole process " SFMT "\n", ss);
  fprintf(fout, " time for the whole process " SFMT "\n", ss);

  printf("\n **************************** COPL ENDS ");
  printf("****************************\n\n");

  fprintf(fout, "\n------------------------------  EOF  ");
  fprintf(fout, "-------------------------------\n\n");
  FileClose(&fout);

  // system("pause");
  return;
} /* PrintEnd */

int ErrorProc(int code,
              char *str)
{
  printf("\n Exit -- " IFMT ": ", code);

  switch (code)
  {
  case OPT_FOUND:
    printf("optimal solution found.\n\n");
    return code;
    break;
  case INF_PROB:
    printf("infeasible problem.\n\n");
    return code;
    break;
  case UNB_PROB:
    printf("unbounded problem.\n\n");
    return code;
    break;
  case EXC_ITER:
    printf("too many iterations.\n\n");
    return code;
    break;
  case NUM_PROB:
    printf("numerical difficulty.\n\n");
    return code;
    break;
  case NOT_MEMSPC:
    printf("not enough memory space ");
    break;
  case NOT_DSKSPC:
    printf("can't create file ");
    break;
  case NOT_SCHFIL:
    printf("can't open file ");
    break;
  case NULL_FILE:
    printf("empty file ");
    break;
  case NULL_SIZE:
    printf("null size.");
    break;
  case NOT_NAMECARD:
    printf("invalid name card: ");
    break;
  case NOT_ROWSCARD:
    printf("no ROWS section in MPS file.");
    break;
  case NOT_COLSCARD:
    printf("no COLUMNS section in MPS file.");
    break;
  case NOT_ENDCARD:
    printf("no ENDATA card in MPS file.");
    break;
  default:
    break;
  }
  if (str)
    printf(SFMT "\n\n", str);
  else
    printf("\n\n");

  ShutDown();

  exit(0);
} /* ErrorProc */
