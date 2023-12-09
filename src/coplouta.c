#include "LPdefine.h"

void OutPutA(matrix *at)
{
  int j;
  array *aj;
  FILE *fp;

  fp = fopen("matrix.dat", "wb");
  if (!fp)
    ErrorProc(NOT_DSKSPC, "matrix.dat in OutPutA");

  fwrite(&at->nrow, sizeof(int), 1, fp);
  fwrite(&at->mems, sizeof(int), 1, fp);

  for (j = 0; j < at->nrow; j++)
  {
    aj = at->ia + j;
    fwrite(&aj->nn0, sizeof(int), 1, fp);
    fwrite(aj->ja, sizeof(int), aj->nn0, fp);
    fwrite(aj->an, sizeof(double), aj->nn0, fp);
  }

  fclose(fp);
} /* OutPutA */
