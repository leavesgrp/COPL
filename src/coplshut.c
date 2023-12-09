#include "LPdefine.h"

void ShutDown(void)
{
  cFree(&pr->rhsset);
  cFree(&pr->rngset);
  cFree(&pr->bndset);
  cFree(&pr->objname);

  HashFree(&ot->rowtab);
  HashFree(&ot->coltab);

  cFree(&ot->mpscard);
  cFree(&ot->prbname);
  cFree(&ot->objname);
  cFree(&ot->rowtype);

  cPtFree(&ot->rowname, ot->rnamesz);
  cPtFree(&ot->colname, ot->cnamesz);

  iFree(&ot->bid);
  iFree(&ot->imem);

  dFree(&ot->rmem);
  dFree(&ot->b);
  dFree(&ot->c);

  dFree(&ot->l);
  dFree(&ot->u);
  dFree(&ot->r);

  dFree(&ot->x);
  dFree(&ot->y);
  dFree(&ot->z);

  dFree(&ot->xu);
  dFree(&ot->yu);
  dFree(&ot->zu);
  dFree(&ot->xf);

  RowFree(&ot->col);
  RowFree(&ot->row);

  MtxFree(&ot->a);
  MtxFree(&ot->at);

  ChlFree(&ot->cl);

  iFree(&ot->cid);
  iFree(&ot->lperm);
  iFree(&ot->linvp);

  dFree(&ot->pcgbf);
  dPtFree(&ot->ld);
  if (ot->dac)
  {
    if (ot->dac->ia)
    {
      ot->dac->ia->ja = NULL;
      ot->dac->ia->an = NULL;
    }
    MtxFree(&ot->dac);
  }
  if (ot->sac)
  {
    if (ot->sac->ia)
    {
      ot->sac->ia->ja = NULL;
      ot->sac->ia->an = NULL;
    }
    MtxFree(&ot->sac);
  }
} /* ShutDown */
