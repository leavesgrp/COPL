#define optdef_getrbki(_rbk, _i) (_rbk ? _rbk[_i] : FX)
#define optdef_getrbli(_rbk, _rbl, _i) (_rbl ? _rbl[_i] : 0.0)
#define optdef_getrbui(_rbk, _rbu, _i) (_rbu ? _rbu[_i] : 0.0)

#define optdef_getcbkj(_cbk, _j) (_cbk ? _cbk[_j] : LO)
#define optdef_getcblj(_cbk, _cbl, _j) (_cbl ? _cbl[_j] : 0.0)
#define optdef_getcbuj(_cbk, _cbu, _j) (_cbu ? _cbu[_j] : MaxPositive)

#define cmps_haslo(_bk) ((_bk) == LO || (_bk) == BO || \
                         (_bk) == FX || (_bk) == BV || (_bk) == UI)
#define cmps_hasup(_bk) ((_bk) == UP || (_bk) == BO || \
                         (_bk) == FX || (_bk) == BV || (_bk) == UI)
#define cmps_hasbo(_bk) ((_bk) == BO || (_bk) == FX || \
                         (_bk) == BV || (_bk) == UI)
#define cmps_isfx(_bk) ((_bk) == FX)
#define cmps_isfr(_bk) ((_bk) == FR)

/*
 * functions in Xbasics.c
 */
pivbas *PbsAlloc(int, char *);
void PbsFree(pivbas **);
optbas *ObsAlloc(optmod *, char *);
void ObsFree(optbas **);
gsdec *LufAlloc(int, int, int, int, char *);
int LufRenew(gsdec *, int, char *);
void LufFree(gsdec **);
int blkord_maxtran(matrix *at,
                   int rsze,
                   int rset[],
                   int invrset[],
                   int csze,
                   int cset[],
                   int cmark[],
                   int wisze,
                   int wimem[]);
int blkord_maxlwrsqblk(matrix *at,
                       int rsze,
                       int rset[],
                       int invrset[],
                       int csze,
                       int cset[],
                       int cmark[],
                       int *maxtrsze,
                       int wisze,
                       int wimem[]);
int blkord_maxuprsqblk2(int ncol,
                        matrix *a,
                        matrix *at,
                        int rsze,
                        int rset[],
                        int invrset[],
                        int csze,
                        int cset[],
                        int cmark[],
                        int wisze,
                        int wimem[],
                        int wzisze,
                        int wzimem[]);
int basord_blklwr(int nrow,
                  int ncol,
                  int nrowc,
                  int inew[],
                  int iinv[],
                  int rsze,
                  int rset[],
                  matrix *a,
                  matrix *at,
                  pivbas *basis,
                  int wisze,
                  int wimem[],
                  int wzisze,
                  int wzimem[]);
int lu_putncol(gsdec *lf,
               int ncol);
int lu_putnrow(gsdec *lf,
               int nrow);
void pivlp_formxb(int nrow,
                  int ncol,
                  matrix *at,
                  int rbk[],
                  double rbl[],
                  double rbu[],
                  pivbas *basis,
                  double rx[],
                  double cx[],
                  double xb[],
                  gsdec *lu,
                  int wrsze,
                  double wrmem[]);
void pivlp_xb2rxcx(int nrow,
                   pivbas *basis,
                   double rx[],
                   double cx[],
                   double xb[]);
int pivlp_insart(gsdec *lu,
                 int nrow,
                 int nrowc,
                 int ncol,
                 int iinv[],
                 int inew[],
                 int rsety,
                 pivbas *basis,
                 int rsk[],
                 int csk[],
                 double rx[],
                 double cx[],
                 double xb[],
                 double yb[],
                 int wisze,
                 int wimem[]);
int pivlp_factor1(int nrow,
                  int nrowc,
                  int iinv[],
                  int inew[],
                  matrix *at,
                  int inca,
                  pivbas *basis,
                  gsdec *lu,
                  int wisze,
                  int wimem[],
                  int wrsze,
                  double wrmem[]);
int pivlp_factor2(int prlev,
                  int nrow,
                  int nrowc,
                  int ncol,
                  int iinv[],
                  int inew[],
                  matrix *at,
                  pivbas *basis,
                  int rsk[],
                  int csk[],
                  double rx[],
                  double cx[],
                  double xb[],
                  double yb[],
                  int rsety,
                  gsdec *lu,
                  int wisze,
                  int wimem[],
                  int wrsze,
                  double wrmem[]);
int pivlp_prolresp(int prlev,
                   int nrow,
                   int nrowc,
                   int ncol,
                   int iinv[],
                   int inew[],
                   matrix *at,
                   pivbas *basis,
                   int rsk[],
                   int csk[],
                   double rx[],
                   double cx[],
                   double xb[],
                   double yb[],
                   int rsety,
                   gsdec *lu,
                   int lresp,
                   int *refac,
                   int wisze,
                   int wimem[],
                   int wrsze,
                   double wrmem[]);
void pivlp_swapbasis(int, int, pivbas *);
void pivlp_mvartrow(int nrowc,
                    int iinv[],
                    pivbas *basis);
void pivlp_formyb(int nrow,
                  int ncol,
                  matrix *at,
                  double cc[],
                  pivbas *basis,
                  double yb[],
                  gsdec *lu,
                  int wrsze,
                  double wrmem[]);
double GetDualObj(int nrow,
                  int ncol,
                  double cf,
                  int rbk[],
                  double rbl[],
                  double rbu[],
                  int cbk[],
                  double cbl[],
                  double cbu[],
                  int rsk[],
                  int csk[],
                  double yb[],
                  double csb[]);

/*
 * functions in Xclean.c
 */
void XcrCleaProc(prfparamt *, int, int, double *, double, matrix *, matrix *, int *,
                 double *, double *, int *, double *, double *, int *, int *, double *,
                 double *, double *, double *, double *, pivbas *, int *, int *,
                 double *, double *, double *, double *, gsdec *, prfinfot *,
                 int, int *, int, double *);
double CompBndInf(double, int, double, double);
int pivlp_solinfo(pivlppart *, pivlpt *, pivlpsolt *, pivlpslvit *, int, double *);
int prf_prilog(void *, pivlpslvit *);

/*
 * functions in Xcross.c
 */
static void XcrMsg(prfinfot *, gsdec *);
void XcrProc(prfparamt *, int, int, double *, double, matrix *, matrix *, int *, double *,
             double *, int *, double *, double *, int *, int *, double *, double *, double *,
             double *, double *, pivbas *, int *, int *, double *, double *, double *,
             double *, gsdec *, prfinfot *, int, int *, int, double *);
void XcrMain(int, int, int, int, int, int, int *, double, matrix *, matrix *,
             double *, double *, double *, int *, double *, double *, double *,
             pivbas *, int *, double *, double *, double *);
void CrossOverProc(optmod *);

/*
 * functions in Xdual.c
 */
void XcrDualProc(prfparamt *, int, int, int, double *, double, matrix *, matrix *,
                 int *, double *, double *, int *, double *, double *, int *, int *,
                 double *, double *, double *, pivbas *, int *, int *, double *,
                 double *, double *, double *, gsdec *, prfinfot *, int, int *,
                 int, double *);
void XcrDualMsg(prfparamt *, int, int, double *, double, matrix *, matrix *, int *,
                double *, double *, int *, double *, double *, pivbas *, int *, int *,
                double *, prfinfot *, int, double *);

/*
 * functions in Xlubasi.c
 */
void lu_delcur(gsdec *lf,
               int j);
void lu_delpos(gsdec *lf,
               int j,
               int pos);
void lu_inscur(gsdec *lf,
               int j,
               int pos);
void lu_buildcol(gsdec *lf,
                 int incall,
                 double colmax[]);
void lu_choosepiv(gsdec *lf,
                  int step,
                  double colmax[],
                  int cmiu[],
                  int *pivi,
                  int *pivj);
int lu_dopivot(gsdec *lf,
               int step,
               double colmax[],
               int cmiu[],
               int pivi,
               int pivj,
               int zmap[],
               double rmemn1[]);
int lu_garcol(gsdec *lf,
              int step,
              int buildcol,
              int incall);
static void lu_moverowtoend(gsdec *lf,
                            int i,
                            int newsze,
                            int updcol);
int lu_incrowsze(gsdec *lf,
                 int step,
                 int i,
                 int newsze,
                 int updcol,
                 int incall);

/*
 * functions in Xgsdec.c
 */
int lu_init(gsdec *lf);
int lu_inaij(gsdec *lf,
             double aij,
             int i,
             int j);
int lu_inaj(gsdec *lf,
            int j,
            int sze,
            double valj[],
            int subj[],
            int mapj[]);
int lu_inaj2(gsdec *lf,
             int j,
             int sze,
             double valj[],
             int subj[],
             int cutoff,
             int inew[]);
void lu_restore(gsdec *lf);
int lu_insart(gsdec *lf,
              int insj[],
              int wisze,
              int wimem[]);
int lu_factor(gsdec *lu,
              int wisze,
              int wimem[],
              int wrsze,
              double wrmem[]);
int lu_getrank(gsdec *lf,
               int *rank);
int lu_getucond(gsdec *lf,
                double *abs_min_uii,
                double *abs_max_uii);
int lu_getuii(gsdec *lf,
              int i,
              int *j,
              double *uii);

/*
 * functions in Xlusolv.c
 */
void lu_btranu(gsdec *lf,
               int bsze,
               double b[],
               int bsub[],
               int bmap[],
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[]);
void lu_btranl(gsdec *lf,
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[]);
void lu_btran2(gsdec *lf,
               int bsze,
               double b[],
               int bsub[],
               int bmap[],
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[]);
void lu_btran(gsdec *lf,
              double b[],
              double x[]);
void lu_ftranl(gsdec *lf,
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[]);
void lu_ftranu(gsdec *lf,
               int *bnnz,
               double b[],
               int bsub[],
               int bmap[],
               int *xnnz,
               double x[],
               int xsub[],
               int xmap[]);
void lu_ftran(gsdec *lf,
              double b[],
              double x[]);
int lu_initupdate(gsdec *lf);
void perm_rotateleft(int start,
                     int end,
                     int p[]);
int lu_supdl(gsdec *lf,
             int i,
             int j,
             int nnzla,
             double la[],
             int subla[],
             int mapla[],
             int inrow1[],
             int incol1[]);
lu_resp lu_supd(gsdec *lf,
                int i,
                int j,
                int nnzaj,
                double aj[],
                int subaj[],
                int mapaj[],
                int imemm1[],
                int imemn1[]);

/*
 * functions in Xprimal.c
 */
double GetDualObj(int, int, double, int *, double *, double *, int *, double *, double *,
                  int *, int *, double *, double *);

/*
 * functions in Xprimal.c
 */
void XcrPrimProc(prfparamt *, int, int, int, double *, double, matrix *, matrix *,
                 int *, double *, double *, int *, double *, double *, int *, int *,
                 double *, double *, pivbas *, int *, int *, double *, double *,
                 double *, gsdec *, prfinfot *, int, int *, int, double *);
void XcrPrimMsg(prfparamt *, int, int, double *, double, matrix *, matrix *, int *,
                double *, double *, int *, double *, double *, pivbas *, int *, int *,
                double *, double *, double *, prfinfot *, int, double *);

/*
 * functions in Xsimplx.c
 */
void pivlp_primalsimplex(pivlppart *, pivlpt *, pivlpsolt *, pivlpslvit *,
                         void *, int, int, int *, int *, int, int *, int *,
                         gsdec *, int, int *, int, double *);
