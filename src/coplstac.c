#include "LPdefine.h"

void CreateStack(void)
{
  fstk = fopen("coplfile.stk", "wb");
  if (!fstk)
    ErrorProc(NOT_DSKSPC, "coplfile.stk");
  fbin = fopen("coplfile.bin", "wb");
  if (!fbin)
    ErrorProc(NOT_DSKSPC, "coplfile.bin");
} /* CreateStack */

void Char1(char sc)
{
  fwrite(&sc, 1, 1, fbin);
  ot->bytes++;
} /* Char1 */

void Int1(int i1)
{
  fwrite(&i1, 4, 1, fbin);
  ot->bytes += 4;
} /* Int1 */

void AddStack(void)
{
  fwrite(&ot->bytes, 4, 1, fstk);
  ot->items++;
} /* AddStack */

void Int2(int i1,
          int i2)
{
  Int1(i1);
  Int1(i2);
} /* Int2 */

void Int3(int i1,
          int i2,
          int i3)
{
  Int1(i1);
  Int1(i2);
  Int1(i3);
} /* Int3 */

void Int4(int i1,
          int i2,
          int i3,
          int i4)
{
  Int2(i1, i2);
  Int2(i3, i4);
} /* Int4 */

void Dbl1(double v1)
{
  fwrite(&v1, 8, 1, fbin);
  ot->bytes += 8;
} /* Dbl1 */

void Dbl2(double v1,
          double v2)
{
  Dbl1(v1);
  Dbl1(v2);
} /* Dbl2 */

void Dbl3(double v1,
          double v2,
          double v3)
{
  Dbl1(v1);
  Dbl1(v2);
  Dbl1(v3);
} /* Dbl3 */

void Int1Dbl1(int i1,
              double v1)
{
  Int1(i1);
  Dbl1(v1);
} /* Int1Dbl1 */

void Int1Dbl2(int i1,
              double v1,
              double v2)
{
  Int1(i1);
  Dbl2(v1, v2);
} /* Int1Dbl2 */

void Int1Dbl3(int i1,
              double v1,
              double v2,
              double v3)
{
  Int1(i1);
  Dbl3(v1, v2, v3);
} /* Int1Dbl3 */

void Int2Dbl1(int i1,
              int i2,
              double v1)
{
  Int2(i1, i2);
  Dbl1(v1);
} /* Int2Dbl1 */

void Int2Dbl2(int i1,
              int i2,
              double v1,
              double v2)
{
  Int2(i1, i2);
  Dbl2(v1, v2);
} /* Int2Dbl2 */

void Int3Dbl1(int i1,
              int i2,
              int i3,
              double v1)
{
  Int3(i1, i2, i3);
  Dbl1(v1);
} /* Int3Dbl1 */

void Int3Dbl2(int i1,
              int i2,
              int i3,
              double v1,
              double v2)
{
  Int3(i1, i2, i3);
  Dbl2(v1, v2);
} /* Int3Dbl2 */

void Int3Dbl3(int i1,
              int i2,
              int i3,
              double v1,
              double v2,
              double v3)
{
  Int3(i1, i2, i3);
  Dbl3(v1, v2, v3);
} /* Int3Dbl3 */

void Str1(char *s)
{
  int k;

  k = strlen(s);
  Int1(k);
  fwrite(s, 1, k, fbin);
  ot->bytes += k;
} /* Str1 */

void RecNames(char *nmps)
{
  int i;
  char ss[60];

  AddStack();

  if (!(*nmps))
    sprintf(ss, "coplfile");
  else
    sprintf(ss, "%s", nmps);

  Int3(NMES, ot->m, ot->n);
  Str1(ss);

  fwrite(ot->r, 8, ot->m, fbin);
  fwrite(ot->b, 8, ot->m, fbin);
  fwrite(ot->rowtype, 1, ot->m, fbin);
  fwrite(ot->l, 8, ot->n, fbin);
  fwrite(ot->u, 8, ot->n, fbin);
  ot->bytes += 17 * ot->m + 16 * ot->n;

  for (i = 0; i < ot->m; i++)
    Str1(ot->rowname[i]);

  for (i = 0; i < ot->n; i++)
    Str1(ot->colname[i]);

} /* RecNames */

void RecVect(array *a,
             int t)
{
  int i, j, k;

  if (t)
  { /* row vector */
    for (k = 0; k < a->nn0; k++)
    {
      j = a->ja[k];
      if (ot->cflg[j] == 'T')
        Int1Dbl1(j, a->an[k]);
    }
  }

  else
  { /* col vector */
    for (k = 0; k < a->nn0; k++)
    {
      i = a->ja[k];
      if (ot->rflg[i] == 'T')
        Int1Dbl1(i, a->an[k]);
    }
  }
  i = -1;
  Int1(i);

} /* RecVect */

void CloseStack(void)
{
  fwrite(&ot->bytes, 4, 1, fstk);
  fwrite(&ot->items, 4, 1, fstk);
  fclose(fstk);
  fclose(fbin);
} /* CloseStack */
