#include "LPdefine.h"

/* global variables */
FILE *fbin = NULL,
     *fmps = NULL,
     *fout = NULL,
     *fstk = NULL,
     *fres = NULL;
char *nmps = NULL,
     *nout = NULL,
     *nspc = NULL;

optmod *ot = NULL;
optpar *pr = NULL;

int main(int argc, char *argv[])
{
  /* local variables */
  int idproc = 0, idxspc = 0,
      idxmps = 0, idxout = 0,
      ntst = 1, i, k, start = 0;
  clock_t otim[8];
  optpar parm = {0};
  optmod optt = {0};

  char name[][12] = {
      "25FV47", "80BAU3B", "ADLTTLE", "AFIRO", "AGG",
      "AGG2", "AGG3", "BANDM", "BEACONFD", "BLEND",
      "BNL1", "BNL2", "BOEING1", "BOEING2", "BORE3D",
      "BRANDY", "CAPRI", "CYCLE", "CZPROB", "D2Q06C",
      "D6CUBE", "DEGEN2", "DEGEN3", "DFL001", "E226",
      "ETAMACRO", "FFFFF800", "FINNIS", "FIT1D", "FIT1P",
      "FIT2D", "FIT2P", "FORPLAN", "GANGES",
      "GFRD-PNC", "GREENBEA", "GREENBEB", "GROW15", "GROW22",
      "GROW7", "ISRAEL", "KB2", "LOTFI", "MAROS",
      "MAROS-R7", "NESM", "PEROLD", "PILOT", "PILOT.JA",
      "PILOT.WE", "PILOT4", "PILOT87", "PILOTNOV", "RECIPE",
      "SC105", "SC205", "SC50A", "SC50B", "SCAGR25",
      "SCAGR7", "SCFXM1", "SCFXM2", "SCFXM3", "SCORPION",
      "SCRS8", "SCSD1", "SCSD6", "SCSD8", "SCTAP1",
      "SCTAP2", "SCTAP3", "SEBA", "SHARE1B", "SHARE2B",
      "SHELL", "SHIP04L", "SHIP04S", "SHIP08L", "SHIP08S",
      "SHIP12L", "SHIP12S", "SIERRA", "STAIR", "STANDATA",
      "STANDGUB", "STANDMPS", "STOCFOR1", "STOCFOR2", "STOCFOR3",
      "TRUSS", "TUFF", "VTP.BASE", "WOOD1P", "WOODW",

      "bas1lp", "dbir1", "dbir2", "dbic1", "lpl1",
      "lpl2", "lpl3", "nemsafm", "nemscem", "nemsemm1",
      "nemsemm2", "nemspmm1", "nemspmm2", "nemswrld", "nsct1",
      "nsct2", "nsic1", "nsic2", "nsir1", "nsir2",
      "cre-a", "cre-b", "cre-c", "cre-d", "ken-07",
      "ken-11", "ken-13", "ken-18", "osa-07", "osa-14",
      "osa-30", "osa-60", "pds-02", "pds-06", "pds-10",
      "pds-20", "ch", "co5", "co9", "cq5",
      "cq9", "ge", "mod2", "nl", "world",

      "bgdbg1", "bgetam", "bgindy", "bgprtr", "box1",
      "chemcom", "cplex1", "cplex2", "ex72a", "ex73a",
      "forest6", "galenet", "gosh", "gran", "greenbea",
      "itest2", "itest6", "klein1", "klein2", "klein3",
      "mondou2", "pang", "pilot4i", "qual", "reactor",
      "refinery", "vol1", "woodinfe"},
       *cs;

  ot = &optt;
  pr = &parm;

  idxout = CheckArgs(argc, argv, "-t");
  if (idxout)
  {
    cs = argv[idxout];
    k = strlen(cs);
    for (i = 0; i < k; i++)
      cs[i] = tolower(cs[i]);

    if (!strcmp(cs, "on"))
      ntst = 95;
  }

#ifdef TEST
  ntst = 94;
  fres = fopen("coplfile.res", "w");
  start = 0;
#endif

  for (k = start; k < ntst; k++)
  {

#ifdef TEST
    CleanOptdat(&optt, &parm, otim);
#endif

    SetTime(otim, START);
    PrintHead();

    idxout = CheckArgs(argc, argv, "-o");
    if (!idxout)
    {
      nout = cAlloc(LineSize, "for nout in main");
      idproc = GetFilen('o', nout);
    }
    else
      nout = argv[idxout];

    idproc = FprintHead();
    if (!idproc)
      return false;

    idxspc = CheckArgs(argc, argv, "-s");
    if (idxspc)
      nspc = argv[idxspc];

    idproc = ReadParms(&parm, &optt);

    if (idproc)
    {

#ifdef TEST
      nmps = name[k];
      fprintf(fres, "%-12s&", nmps);
      idxspc = strlen(nmps);
      for (i = 0; i < idxspc; i++)
        nmps[i] = tolower(nmps[i]);
#else
      idxmps = CheckArgs(argc, argv, "-f");
      if (!idxmps)
      {
        nmps = cAlloc(LineSize, "for nmps in main");
        idproc = GetFilen('f', nmps);
      }
      else
        nmps = argv[idxmps];
#endif

      idproc = OptProc(otim, &optt, &parm);
    }

#ifdef TEST
    idxmps = true;
    fprintf(fres, "%3d&%+10e&%+10e&%7.1e&%7.1e&%7.1e&",
            optt.iter, optt.pobj, optt.dobj,
            optt.rgap, optt.rrpnm, optt.rrdnm);
#endif

    if (!idxmps)
      cFree(&nmps);
    if (!idxout)
      cFree(&nout);
    ShutDown();

    SetTime(otim, ELAPS);
    PrintEnd(otim);
  }

  FileClose(&fres);
  // system("pause");
  return idproc;
} /* main */
