
int MtxTrans(int, matrix *, int, int *, int, matrix **);
int SymboProc(chol *, matrix *, matrix *);

chol *ChlAlloc(int, char *);
void ChlFree(chol **);
void iPtAlloc(int, int, int **, char *);
void iPtFree(int, int **);
order *OdAlloc(int, int, char *);
void OdFree(order **);
xlist *XtAlloc(int, int);
void XtFree(xlist **);

void OdInit(order *, int *);
int GetOrder(order *, int *);
void OdIndex(order *, int, int);

int setXt(xlist *xt,
          int last,
          int maxp,
          int isze,
          int *imem);
int nextXt(xlist *xt);
void delXt(xlist *xt,
           int e);
int infXt(xlist *xt);
int supXt(xlist *xt,
          int p);
int succXt(xlist *xt,
           int *e);
void putXt(xlist *xt,
           int e,
           int p);
void miXt(xlist *xt,
          int e,
          int dec);
void plXt(xlist *xt,
          int e,
          int inc);
int getXt(xlist *xt,
          int *e,
          int *p);

void dZero(int, double *, int *);
void dcopy(int, double *, int *, int *, double *);
void dSwap(int, int, double *);
double dNorm1(int, double *);
int iSum(int, int *);
int getPos(int, int, int *);
void iZero(int, int *, int *);
void iSet(int, int, int *, int *);
void iSwap(int, int, int *);
void iCopy(int, int *, int *);
void fSort(int, int, int *, int *, int *);
double svDot(array *, double *);
