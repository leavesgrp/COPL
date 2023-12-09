/*
 * functions in coplbasi.c
 */
void CleanOptdat(optmod *, optpar *, clock_t *);
clock_t GetTime(void);
int SetTime(clock_t *, int);
int CheckArgs(int, char *argv[], char *);
int FilesOpen(void);
void FileClose(FILE **);
void PrintHead(void);
int FprintHead(void);
int GetFilen(char, char *);
int CompString(const char *, const char *);
int LeftDots(char *ss, int);
void PrintEnd(clock_t *);
int ErrorProc(int, char *);

/*
 * functions in coplcrus.c
 */
int DataTransf(optmod *, optpar *);
int TreatBound(double, int, optpar *, optmod *, int);
int CrushProc(optpar *, optmod *);

/*
 * functions in copldens.c
 */
int ChkDenseCol(optmod *, optpar *);
void StripMatrix(optmod *);
void StripDiag(optmod *, double *, double *, double *);
int LuPivot(int, double **, double *, int *, int *, int *);
void DenseMatrix(optmod *, double *, double *, int *);
void DensNormSol(optmod *, chol *, double *, matrix *, double *, double *);

/*
 * functions in coplhash.c
 */
int HashValue(hashtab *, char *, unsigned);
int InsertHash(hashtab *, char *, int, unsigned);
int FoundHash(hashtab *, char *, unsigned);

/*
 * functions in coplhomo.c
 */
void ShowHead(void);
void ShowInfor(optmod *);
void HomOptSol(optmod *, optpar *);
static int SolResult(optmod *, optpar *);
void HsdProc(clock_t *, optmod *, optpar *);

/*
 * functions in coplmemo.c
 */
char *cAlloc(int, char *);
void cFree(char **);
char *cAllocCopy(char *, char *);
char *cAllocCopy(char *, char *);
int *iAlloc(int, char *);
void iFree(int **);
double *dAlloc(int, char *);
void dFree(double **);
double **dPtAlloc(int, int, char *);
void dPtFree(double ***);
char **cPtAlloc(int, char *);
void cPtFree(char ***, int);
hashtab *HashAlloc(int, char *);
void HashFree(hashtab **);
optrow *RowAlloc(int, int, char *);
void RowFree(optrow **);
matrix *MtxAlloc(int, int, char *);
void MtxFree(matrix **);
int BndAlloc(optmod *);
chol *ChlAlloc(int, char *);
void ChlFree(chol **);

/*
 * functions in coplmpsf.c
 */
int GetCard(char *, int *);
int ParseCard(char *, char **);
int CopyName(char *, char *);
int ReadSection(mpssec, optpar *, optmod *, int *, int *);
int PrintStatis(optpar *, optmod *);
int ReadMpsf(optpar *, optmod *, char *);

/*
 * functions in coplopti.c
 */
void ResList(int, optpar *, optmod *);
int RecSolution(optmod *);
int OptProc(clock_t *, optmod *, optpar *);

/*
 * functions in coplpara.c
 */
int ReadParms(optpar *, optmod *);

/*
 * functions in coplpres.c
 */
int PreSolve(optpar *, optmod *);

/*
 * functions in coplpsub.c
 */
int NullFix(optpar *, optmod *);
int ChkRowSing(int, optpar *, optmod *, int *, int *);
int SetRowBnds(int, optpar *, optmod *, int *, int *);
int ChkRowForc(int, optpar *, optmod *, int *, int *);
int SetColList(int, int *, int *, optmod *, int *, int *);
int ChkRowDomn(int, optpar *, optmod *, int *);
int ChkRowDblt(int, optpar *, optmod *, int *);
int SetColBnds(int, int, optpar *, optmod *, int *);
int ChkColDomn(int, int, optpar *, optmod *, int *);
int ChkRowDupl(int, int, int, optpar *, optmod *, int *, int *);
int ChkColDupl(int, int, optpar *, optmod *, int *, int *);

/*
 * function in coplshut.c
 */
void ShutDown(void);

/*
 * function in matxlib.a
 */
void dCopy(int, double *, double *);
void dScal(int, double, double *, int *, int *);
double dSums(int, double *, int *);
void addVect(int, double, double *, int *, double *);
double dNorm0(int, double *, int *);
double dNorm2(int, double *);
void dZero(int, double *, int *);
void mTimesv(int, int, int, double, matrix *, double *, double, double *);
void SetSvec(double, array *, double *);

int MtxTrans(int, matrix *, int, int *, int, matrix **);
void CholSol(chol *, double *, double *);
void SolArgLES(chol *, double *, matrix *, double *, double *);
int SymboProc(chol *, matrix *, matrix *);
int decChol(chol *, double *, matrix *, int, int *, int, double *);
void CrossOverProc(optmod *);
void ChkRowLdpt(optpar *, optmod *, int, int *, int, double *);
void setArray(double, array *, double *);

/*
 * function in coplstac.c
 */
void CreateStack(void);
void Char1(char);
void Int1(int);
void AddStack(void);
void Int2(int, int);
void Int3(int, int, int);
void Int4(int, int, int, int);
void Dbl1(double);
void Dbl2(double, double);
void Dbl3(double, double, double);
void Int1Dbl1(int, double);
void Int1Dbl2(int, double, double);
void Int1Dbl3(int, double, double, double);
void Int2Dbl1(int, int, double);
void Int2Dbl2(int, int, double, double);
void Int3Dbl1(int, int, int, double);
void Int3Dbl2(int, int, int, double, double);
void Int3Dbl3(int, int, int, double, double, double);
void Str1(char *);
void RecNames(char *);
void RecVect(array *, int);
void CloseStack(void);
