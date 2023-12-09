#include "LPmatrx.h"

typedef enum
{
  SOL_OPTIMAL,
  SOL_FEAS,
  SOL_INFEAS,
  SOL_UNK
} sol_sta;

typedef enum
{
  PIV_ROW = 0,
  PIV_COL,
  PIV_ART
} bastype;

typedef enum
{
  LO = 0, /*   lj<= xj <=INF       */
  UP = 1, /* -INF<= xj <=uj        */
  FX = 2, /*   lj = xj  =uj        */
  FR = 3, /* -INF<= xj <=INF       */
  BO = 4, /*   lj<= xj <=uj        */
  UI = 5, /*   xj \in {0,1, .. uj} */
  BV = 6, /*   xj \in {0,1}        */
  MI = 7, /* Can never appear.     */
  PL = 8  /* Can never appear.     */
} bndkey;

typedef enum
{
  lu_ok = 1, /* Everything is OK. */

  /* Argument errors. */
  lu_arg1,
  lu_arg2,
  lu_arg3,
  lu_arg4,
  lu_arg5,
  lu_arg6,
  lu_arg7,
  lu_arg8,

  lu_space, /* Could not allocate enough space.       */
  lu_wispc, /* Not enough integer work space.         */
  lu_wrspc, /* Not enough real work space.            */
  lu_wfull, /* Out of work space.                     */
  lu_insta, /* Terminated due to numerical stability. */
  lu_refac  /* The matrix should be refactored.       */
} lu_resp;

typedef enum
{
  lu_rinit = 0,
  lu_rfin,
  lu_rnull
} rsta;

typedef enum
{
  chuzr_unb,    /* Problem is unbounded.    */
  chuzr_blockl, /* Lower bound is blocking. */
  chuzr_blocku,
  chuzr_low,
  chuzr_upr
} chuzr_status;

typedef struct
{
  int vt;
  int j;
} pivbas;

typedef struct
{
  int maxnrow;
  int maxncol;
  int nrow;
  int ncol;
  pivbas *basis; /* Optimal basis [maxnrow]      */
  int *rsk;      /* Row status keys [maxnrow]    */
  int *csk;      /* Column status keys [maxncol] */
  double *rx;    /* rx  [maxnrow]        */
  double *cx;    /* cx  [maxncol]        */
  double *xb;    /* xb  [maxncol]        */
  double *y;     /* y   [maxnrow]        */
  double *cs;    /* csb [maxncol]        */
} optbas;

typedef struct
{
  int nnza;   /* Number f non-zeros in A. */
  int nnzl;   /* Non-zeros in L after last factor. */
  int nnzu;   /* Non-zeros in U after last factor. */
  int tnnfa;  /* Total number of factorizations. */
  int tnnft;  /* Total number of ftrans.         */
  int tnnbt;  /* Total number of factorizations. */
  int tnnu;   /* Total number of updates.        */
  int tnns;   /* Total number of solves.         */
  long tnnza; /* Total number of non-zeros in A. */
  long tnnzl; /* Total non-zeros in L for all factors. */
  long tnnzu; /* Total non-zeros in U for all factors. */
  int tngc;   /* Total number of garbage collections. */
  int ngc;    /* Number of garbage collections since
                 last factorization. */
} gsmsg;

typedef struct
{
  int nrow0;   /* Max number of rows. */
  int ncol0;   /* Max number of columns. */
  int nrow;    /* Number of rows. */
  int ncol;    /* Number of columns. */
  int wsze0;   /* Previous size of work arrays. */
  int wsze;    /* Size of work arrays. */
  int lbeg;    /* Position of the first element in L. */
  int uend;    /* Position of the first element in U. */
  double *val; /* Store the values [wsze]. */
  int *subi;   /* Store the i subscripts [wsze]. */
  int *subj;   /* Store the j subscripts [wsze]. */
  int *chead;  /* Pointer to the first row [ncol,nil=wsze]. */
  int *ctail;
  int *ccur;    /* Pointer to the current row [ncol]. */
  int *cprev;   /* Pointer to previous element [wsze]. */
  int *cnext;   /* Pointer to the next element [wsze]. */
  int *ufir;    /* ufir[i] is position first
                   element in the i'th row of U [nrow]. */
  int *usze;    /* Total size of the i'th
                   row of U incl. overhead [nrow]. */
  int *uuse;    /* Actual size of the i'th
                   row of U [nrow]. */
  int uhead;    /* Pointer to first stored row. */
  int utail;    /* Pointer to last stored row. */
  int *unex;    /* unex[i]: Index of row stored after
                   row i [len=nrow,nil=nrow]. */
  int *upre;    /* upre[i]: Index of row stored previous
                   to row i [len=nrow,nil=nrow]. */
  xlist *nzrow; /* Non-zeros in rows. */
  xlist *nzcol; /* Non-zeros in cols. */
  int *p;       /* Row ordering (implicit) [nrow]. */
  int *invp;    /* Inverse row ordering [nrow].    */
  int *invq;    /* Inverse column ordering [ncol]. */
  int *rowsta;
  int updates;
  int lupds; /* Number of L updates. */

  double toldrop;
  double tolfill;
  double tolrpiv;
  double tolapiv;
  double tolupiv;
  double tolstb;

  int trhsch;
  int trhuov;
  int trhrfq;
  int dropl;
  int dropu;
  int incwsze; /* Increase storage. */
  int rank;    /* Numerical rank after last refactorization. */
  gsmsg info;
} gsdec;

typedef struct
{
  int tries; /*  Number of purification tries. */
  int psze;
  int pbsze;
  int zsze;
  int zbsze;
  int asze; /* Number of artificial in the
                       optimal basis. */
  int itep;
  int mvep;
  int titep;
  int tmvep;

  double pobj;
  int npbi;      /* Number of primal bound infeasibilities. */
  double mpbi;   /* Max primal bound infeasibility.         */
  double spbi;   /* Sum of primal bound infeasiblity.       */
  double pfeas1; /* || Ax - b ||_1                          */
  double pfeasi; /* || Ax - b ||_inf                        */
  double nrmxb;  /* || x_B ||_1                             */
  double nrmb;   /* || b - N x_n ||_1                       */
  double nrmbi;  /* || b - N x_n ||_inf                     */

  int ited;
  int mved;
  int tited;
  int tmved;

  double dobj;
  int ndbi;      /* Number of dual bound infeasibilities. */
  double mdbi;   /* Max dual bound infeasibility.         */
  double sdbi;   /* Sum of dual bound infeasibility.      */
  double dfeas1; /* || B^T y - c_B ||_1                   */
  double dfeasi; /* || B^T y - c_B ||_inf                 */
  double nrmy;   /* || y ||_1                             */
  double nrmcb;  /* || c_B ||_1                           */
  double nrmcbi; /* || c_B ||_inf                         */

  int itec;
  int mvec;
  int titec;
  int tmvec;
  int resp;
  int solsta;

  gsmsg luinf;
} prfinfot;

typedef struct
{
  double tolrpiv;
  double tolapiv;
  double tollurpiv;

  double tolx;
  double tols;

  double rtolx;
  double rtols;

  int prlev;   /* Print level. */
  int cmaxite; /* Max iterations allowed in clean.
                */
  int prfrq;   /* Print freq. for the simplex.      */

  int rdf; /* Temporary option. */
  int wrf;
} prfparamt;

typedef struct
{
  int nrow;
  int ncol;
  double *cc;  /* [ncol] */
  double cf;   /* Fixed term in objective. */
  matrix *a;   /* Row wise representation of A. */
  matrix *at;  /* Column wise representation of A. */
  int *rbk;    /* [nrow] */
  double *rbl; /* [nrow] */
  double *rbu; /* [nrow] */
  int *cbk;    /* [ncol] */
  double *cbl; /* [ncol] */
  double *cbu; /* [ncol] */
} pivlpt;

typedef struct
{
  int solsta;    /* Solution status. */
  pivbas *basis; /* [nrow] */
  int *rsk;      /* [nrow] */
  int *csk;      /* [ncol] */
  double *rx;    /* [nrow] */
  double *cx;    /* [ncol] */
  double *xb;    /* [nrow] */
  double *y;     /* [nrow] */
  double *cs;    /* [ncol] */
} pivlpsolt;

typedef struct
{
  int maxiter;    /* Max iterations allowed.             */
  double tolx;    /* x tolerence.                        */
  double tols;    /* s tolerence.                        */
  int nosbas;     /* No super-basics is allowed.         */
  double w12;     /* Weight between phase 1 and phase 2 [ 0<= w12 <= 1]. */
  double tolapiv; /* Absolut pivot tolerence in ratio test. */
  int wlev;       /* Warm starts level:
                       0 - Cold start
                       1 - Warm start
                           (assumes xb and yb are correct
                            and sol->status are correct). */
  int prfrq;      /* Number of iterations between each printing. */
  int prlev;      /* Level of printing. */
} pivlppart;

typedef struct
{
  int sresp; /* Solver resp                        */
  int iter;  /* Iterations.                        */
  int mves;  /* Moves.                             */
  double pobj;
  double dobj;
  double mpbi;   /* Max primal bound violations.       */
  double mdbi;   /* Max dual bound volations.          */
  double spbi;   /* Sum of primal bound infeasibility. */
  double sdbi;   /* Sum of dual bound infeasibility.   */
  int npbi;      /* Number of primal bound violations. */
  int ndbi;      /* Number of dual bound violations.   */
  double pfeas1; /* || rx - Ax ||_1                    */
  double dfeas1; /* || A^T y  + s - c ||_1             */
  double pfeasi; /* || rx - Ax     ||_inf              */
  double dfeasi; /* || B^T y - c_B ||_inf              */
  double bi;     /* || b - N x_N ||_inf                */
  double cbi;    /* || c_B ||_inf                      */
} pivlpslvit;

typedef struct
{
  int nar;    /* Number of active rows for pricing.        */
  int nac;    /* Number of active columns for pricing.     */
  int nrp;    /* Potential number of rows to be priced.    */
  int ncp;    /* Potential number of columns to be priced. */
  int chkite; /* Iteration number of last optimality check. */
} iteinft;

#include "LPxfunc.h"
