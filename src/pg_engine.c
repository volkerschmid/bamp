/* ===========================================================================
 * Native-C port of the Polya-Gamma age-period-cohort MCMC engine.
 *
 * This is a compiled implementation of R/pg_engine.R's .bamp_pg_chain inner
 * loop. The R function is kept as the reference engine (engine="R"); this code
 * is engine="C". It uses the R C API + the LAPACK/BLAS that the package already
 * links (src/Makevars), with NO Rcpp / RcppEigen dependency.
 *
 * Two entry points (registered in init.c):
 *   pg_assemble_c : the deterministic joint precision assembler, exposed on its
 *                   own for the parity test against the R assemble_prec().
 *   pg_chain_c    : one full MCMC chain (the full inner-loop port).
 *
 * REPRODUCIBILITY -- the C code MUST consume R's RNG stream in EXACTLY the same
 * order and count as R/pg_engine.R, so set.seed(seed) + this code reproduce the
 * R engine. Per sweep the locked draw order is:
 *   1. omega: I*J normals, COLUMN-MAJOR (i fastest), via .pg_rpg
 *   2. block draw: P normals
 *   3. overdisp (if on): I*J normals (delta), then 1 gamma (zeta)
 *   4. precisions: gammas kappa,lambda,ny (present), then het kappa2,lambda2,ny2
 *   5. ASIS: per present effect in order a,p,c -> 1 normal (propk) + 1 unif (acc)
 *   6. Laplace-MH: (P-k) normals (gstar) + 1 unif (accept)
 * NOTE: R's rgamma(shape, rate) == C rgamma(shape, scale=1/rate).
 *       R's rnorm(n, m_vec, s_vec) consumes one norm_rand() per element in order.
 * Any change to the R draw order MUST be mirrored here (and vice versa).
 * ===========================================================================*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Random.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ---- small helpers -------------------------------------------------------*/

/* solve Q g = rhs in place on g, given upper Cholesky factor U (Q=U'U),
 * column-major P x P. Two triangular solves: U' y = rhs, then U g = y. */
static void chol_solve(const double *U, int P, double *g) {
  int one = 1;
  F77_CALL(dtrsv)("U", "T", "N", &P, U, &P, g, &one FCONE FCONE FCONE);
  F77_CALL(dtrsv)("U", "N", "N", &P, U, &P, g, &one FCONE FCONE FCONE);
}

/* numerically stable log(1+exp(e)) */
static inline double softplus(double e) {
  return (e > 0.0) ? e + log1p(exp(-e)) : log1p(exp(e));
}

/* quadratic form x' M x for a column-major n x n M and length-n x */
static double quad_form(const double *M, const double *x, int n) {
  double s = 0.0;
  for (int c = 0; c < n; c++) {
    double mc = 0.0;
    for (int r = 0; r < n; r++) mc += M[r + (size_t)n * c] * x[r];
    s += mc * x[c];
  }
  return s;
}

/* ===========================================================================
 * Joint precision assembler (deterministic; mirrors R assemble_prec()).
 *
 * Layout of the P-vector beta: [mu] [theta(I)] [phi(J)] [psi(K)]
 *                              [theta2(I?)] [phi2(J?)] [psi2(K?)]   (het blocks)
 * Block presence is given by has_a/has_p/has_c and het_a/het_p/het_c.
 * w is the I x J per-cell weight; saK/spK/scK the (scaled) structure matrices;
 * cov_p/cov_c are NULL or the mean-1 covariate vectors. coh is the 1-based
 * cohort index per cell (I x J, column-major).
 *
 * Smooth and het indices of the same effect SHARE the likelihood coupling and
 * differ only in the prior block (GMRF K for smooth, kap2*I for het). We build
 * the smooth blocks then add the het duplicate rows/cols by copying.
 * ===========================================================================*/
typedef struct {
  int I, J, K, P;
  int has_a, has_p, has_c, het_a, het_p, het_c;
  int ia, ith, iph, ips, ith2, iph2, ips2;   /* 0-based block start offsets */
  const int *coh;                            /* I*J, 1-based */
  const double *saK, *spK, *scK;
  const double *cov_p, *cov_c;               /* NULL if absent */
  int has_pcov, has_ccov;
} apc_dims;

/* fill an already-zeroed P x P column-major Q */
static void assemble_prec_c(double *Q, const apc_dims *D, const double *w,
                            double kap, double lam, double nyv,
                            double kap2, double lam2, double ny2v) {
  int I = D->I, J = D->J, K = D->K, P = D->P;
  int ia = D->ia, ith = D->ith, iph = D->iph, ips = D->ips;
  int ith2 = D->ith2, iph2 = D->iph2, ips2 = D->ips2;
  const int *coh = D->coh;
  #define QQ(r,c) Q[(r) + (size_t)P * (c)]

  /* marginal weight sums */
  double sw = 0.0;
  double *rs = (double*)R_alloc(I, sizeof(double));
  double *cs = (double*)R_alloc(J, sizeof(double));
  double *cm = (double*)R_alloc(K, sizeof(double));
  for (int i = 0; i < I; i++) rs[i] = 0.0;
  for (int j = 0; j < J; j++) cs[j] = 0.0;
  for (int k = 0; k < K; k++) cm[k] = 0.0;
  for (int j = 0; j < J; j++)
    for (int i = 0; i < I; i++) {
      double wij = w[i + (size_t)I * j];
      sw += wij; rs[i] += wij; cs[j] += wij;
      cm[ coh[i + (size_t)I * j] - 1 ] += wij;
    }

  /* covariate-scaled couplings (no-op when absent): period col j *cov_p[j],
   * cohort col k *cov_c[k]. cs_l = cs*cov_p (x^1), cs_d = cs*cov_p^2 (x^2). */
  const double *cp = D->cov_p, *cc = D->cov_c;
  int hp = D->has_pcov, hc = D->has_ccov;

  /* intercept diagonal + couplings to each present block */
  QQ(ia, ia) = sw;
  if (D->has_a) for (int i=0;i<I;i++){ double v=rs[i]; QQ(ia,ith+i)=v; QQ(ith+i,ia)=v; }
  if (D->has_p) for (int j=0;j<J;j++){ double v=hp?cs[j]*cp[j]:cs[j]; QQ(ia,iph+j)=v; QQ(iph+j,ia)=v; }
  if (D->has_c) for (int k=0;k<K;k++){ double v=hc?cm[k]*cc[k]:cm[k]; QQ(ia,ips+k)=v; QQ(ips+k,ia)=v; }
  if (D->het_a) for (int i=0;i<I;i++){ double v=rs[i]; QQ(ia,ith2+i)=v; QQ(ith2+i,ia)=v; }
  if (D->het_p) for (int j=0;j<J;j++){ double v=hp?cs[j]*cp[j]:cs[j]; QQ(ia,iph2+j)=v; QQ(iph2+j,ia)=v; }
  if (D->het_c) for (int k=0;k<K;k++){ double v=hc?cm[k]*cc[k]:cm[k]; QQ(ia,ips2+k)=v; QQ(ips2+k,ia)=v; }

  /* diagonal weight blocks: age diag=rs, period diag=cs*cov_p^2, cohort diag=cm*cov_c^2.
   * The smooth and het indices of the SAME effect share the likelihood coupling
   * (theta2_i enters cell (i,j) exactly like theta_i), so when both are present
   * the (smooth,het) and (het,smooth) cross-diagonal blocks get the SAME diagonal
   * weight as each self-block -- mirroring R's nested `for g1,g2 in {ith,ith2}`
   * which fills all four sub-blocks. They differ only in the prior added later. */
  for (int i=0;i<I;i++) { double v=rs[i];
    if (D->has_a) QQ(ith+i,ith+i)   += v;
    if (D->het_a) QQ(ith2+i,ith2+i) += v;
    if (D->has_a && D->het_a) { QQ(ith+i,ith2+i) += v; QQ(ith2+i,ith+i) += v; }
  }
  for (int j=0;j<J;j++) { double v=hp?cs[j]*cp[j]*cp[j]:cs[j];
    if (D->has_p) QQ(iph+j,iph+j)   += v;
    if (D->het_p) QQ(iph2+j,iph2+j) += v;
    if (D->has_p && D->het_p) { QQ(iph+j,iph2+j) += v; QQ(iph2+j,iph+j) += v; }
  }
  for (int k=0;k<K;k++) { double v=hc?cm[k]*cc[k]*cc[k]:cm[k];
    if (D->has_c) QQ(ips+k,ips+k)   += v;
    if (D->het_c) QQ(ips2+k,ips2+k) += v;
    if (D->has_c && D->het_c) { QQ(ips+k,ips2+k) += v; QQ(ips2+k,ips+k) += v; }
  }

  /* age<->period coupling w_ap = w * cov_p[col]; and age<->cohort / period<->cohort
   * via the cohort index. Build for smooth blocks, then duplicate to het blocks. */
  for (int j = 0; j < J; j++)
    for (int i = 0; i < I; i++) {
      double wij = w[i + (size_t)I * j];
      double wap = hp ? wij * cp[j] : wij;
      int k = coh[i + (size_t)I * j] - 1;
      /* age-period (and het duplicates) */
      if (D->has_a && D->has_p) { QQ(ith+i,iph+j)+=wap; QQ(iph+j,ith+i)+=wap; }
      if (D->het_a && D->has_p) { QQ(ith2+i,iph+j)+=wap; QQ(iph+j,ith2+i)+=wap; }
      if (D->has_a && D->het_p) { QQ(ith+i,iph2+j)+=wap; QQ(iph2+j,ith+i)+=wap; }
      if (D->het_a && D->het_p) { QQ(ith2+i,iph2+j)+=wap; QQ(iph2+j,ith2+i)+=wap; }
      /* age-cohort: TP col k *cov_c[k] */
      if (D->has_c) {
        double wtc = hc ? wij * cc[k] : wij;
        if (D->has_a) { QQ(ith+i,ips+k)+=wtc; QQ(ips+k,ith+i)+=wtc; }
        if (D->het_a) { QQ(ith2+i,ips+k)+=wtc; QQ(ips+k,ith2+i)+=wtc; }
        if (D->het_c) {
          if (D->has_a){ QQ(ith+i,ips2+k)+=wtc; QQ(ips2+k,ith+i)+=wtc; }
          if (D->het_a){ QQ(ith2+i,ips2+k)+=wtc; QQ(ips2+k,ith2+i)+=wtc; }
        }
        /* period-cohort: PP row j *cov_p[j], col k *cov_c[k] */
        double wpc = wij;
        if (hp) wpc *= cp[j];
        if (hc) wpc *= cc[k];
        if (D->has_p) { QQ(iph+j,ips+k)+=wpc; QQ(ips+k,iph+j)+=wpc; }
        if (D->het_p) { QQ(iph2+j,ips+k)+=wpc; QQ(ips+k,iph2+j)+=wpc; }
        if (D->het_c) {
          if (D->has_p){ QQ(iph+j,ips2+k)+=wpc; QQ(ips2+k,iph+j)+=wpc; }
          if (D->het_p){ QQ(iph2+j,ips2+k)+=wpc; QQ(ips2+k,iph2+j)+=wpc; }
        }
      }
    }

  /* prior blocks: smooth get prior*structure-matrix; het get prec2*I */
  if (D->has_a) for (int c=0;c<I;c++) for (int r=0;r<I;r++) QQ(ith+r,ith+c) += kap*D->saK[r+(size_t)I*c];
  if (D->has_p) for (int c=0;c<J;c++) for (int r=0;r<J;r++) QQ(iph+r,iph+c) += lam*D->spK[r+(size_t)J*c];
  if (D->has_c) for (int c=0;c<K;c++) for (int r=0;r<K;r++) QQ(ips+r,ips+c) += nyv*D->scK[r+(size_t)K*c];
  if (D->het_a) for (int i=0;i<I;i++) QQ(ith2+i,ith2+i) += kap2;
  if (D->het_p) for (int j=0;j<J;j++) QQ(iph2+j,iph2+j) += lam2;
  if (D->het_c) for (int k=0;k<K;k++) QQ(ips2+k,ips2+k) += ny2v;
  #undef QQ
}

/* read a named integer/logical flag list element */
static int geti(SEXP list, const char *nm) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), nm) == 0) {
      SEXP e = VECTOR_ELT(list, i);
      return (TYPEOF(e) == LGLSXP) ? LOGICAL(e)[0] : INTEGER(coerceVector(e, INTSXP))[0];
    }
  return 0;
}
static double getd(SEXP list, const char *nm) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), nm) == 0)
      return REAL(coerceVector(VECTOR_ELT(list, i), REALSXP))[0];
  return 0.0;
}
static SEXP getv(SEXP list, const char *nm) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), nm) == 0) return VECTOR_ELT(list, i);
  return R_NilValue;
}

/* build apc_dims from a config list (keys: I,J,K, has_a/p/c, het_a/p/c, coh,
 * saK,spK,scK, cov_p,cov_c) */
static void fill_dims(apc_dims *D, SEXP cfg) {
  D->I = geti(cfg,"I"); D->J = geti(cfg,"J"); D->K = geti(cfg,"K");
  D->has_a = geti(cfg,"has_a"); D->has_p = geti(cfg,"has_p"); D->has_c = geti(cfg,"has_c");
  D->het_a = geti(cfg,"het_a"); D->het_p = geti(cfg,"het_p"); D->het_c = geti(cfg,"het_c");
  int I=D->I,J=D->J,K=D->K;
  int off = 1;
  D->ia = 0;
  D->ith  = off; if (D->has_a) off += I;
  D->iph  = off; if (D->has_p) off += J;
  D->ips  = off; if (D->has_c) off += K;
  D->ith2 = off; if (D->het_a) off += I;
  D->iph2 = off; if (D->het_p) off += J;
  D->ips2 = off; if (D->het_c) off += K;
  D->P = off;
  D->coh = INTEGER(getv(cfg,"coh"));
  D->saK = D->has_a ? REAL(getv(cfg,"saK")) : NULL;
  D->spK = D->has_p ? REAL(getv(cfg,"spK")) : NULL;
  D->scK = D->has_c ? REAL(getv(cfg,"scK")) : NULL;
  SEXP cp = getv(cfg,"cov_p"), cc = getv(cfg,"cov_c");
  D->cov_p = (cp==R_NilValue) ? NULL : REAL(cp);
  D->cov_c = (cc==R_NilValue) ? NULL : REAL(cc);
  D->has_pcov = (D->cov_p != NULL && D->has_p);
  D->has_ccov = (D->cov_c != NULL && D->has_c);
}

/* .Call entry: assemble Q for a given weight matrix + precisions (no ridge,
 * to match exactly what we want to compare; the R test adds the ridge itself
 * or we add it here on request). */
SEXP pg_assemble_c(SEXP cfg, SEXP w_, SEXP prec_, SEXP ridge_) {
  apc_dims D; fill_dims(&D, cfg);
  int P = D.P;
  const double *prec = REAL(prec_);  /* kap,lam,ny,kap2,lam2,ny2 */
  SEXP Q_ = PROTECT(allocMatrix(REALSXP, P, P));
  double *Q = REAL(Q_);
  memset(Q, 0, (size_t)P * P * sizeof(double));
  assemble_prec_c(Q, &D, REAL(w_), prec[0],prec[1],prec[2],prec[3],prec[4],prec[5]);
  double ridge = REAL(ridge_)[0];
  if (ridge != 0.0) {
    double dsum = 0.0; for (int p=0;p<P;p++) dsum += Q[p+(size_t)P*p];
    double r = ridge * dsum / P;
    for (int p=0;p<P;p++) Q[p+(size_t)P*p] += r;
  }
  UNPROTECT(1);
  return Q_;
}

/* ===========================================================================
 * .pg_rpg : vectorised Polya-Gamma draw, exact mean/variance normal approx.
 * Fills out[0..n-1] consuming n norm_rand() draws IN ORDER (matches R's
 * rnorm(n, b*m, sqrt(b*v)) element order). b,c are length n. ----------------*/
static void pg_rpg_c(int n, const double *b, const double *c, double *out) {
  for (int i = 0; i < n; i++) {
    double ac = fabs(c[i]), m, v;
    if (ac < 1e-4) { m = 0.25; v = 1.0/24.0; }
    else {
      double th = tanh(ac/2.0), sech2 = 1.0 - th*th;
      m = th/(2.0*ac);
      v = ((2.0/ac)*th - sech2)/(4.0*ac*ac);
    }
    double mean = b[i]*m, sd = sqrt(fmax(b[i]*v, DBL_EPSILON));
    double z = norm_rand();                      /* one draw per cell, in order */
    double val = mean + sd*z;
    out[i] = (val > 1e-9) ? val : 1e-9;
  }
}

/* scratch buffers for the Laplace-MH state() evaluation, allocated once */
typedef struct {
  double *e, *Wt, *ymnp, *g, *H, *HZ, *rsy, *csy, *cgy, *Kx, *gg;
} mh_scratch;

/* Laplace-MH state(): mirrors R state(bv). Fills Rout (Pf x Pf upper Cholesky of
 * the Zbasis-projected Hessian), mean_g (Pf), gamma = Z'bv (Pf), and *lp_out.
 * Reuses the verified assemble_prec_c for the Hessian. */
static int mh_state(const apc_dims *D, const double *Y, const double *N,
                    const double *Z, int Pf, const double *bv,
                    double kappa, double lambda, double ny,
                    double kappa2, double lambda2, double ny2,
                    const double *delta, int has_od,
                    mh_scratch *S,
                    double *Rout, double *mean_g, double *gamma, double *lp_out) {
  int I=D->I,J=D->J,K=D->K,P=D->P,IJ=I*J;
  int ia=D->ia,ith=D->ith,iph=D->iph,ips=D->ips,ith2=D->ith2,iph2=D->iph2,ips2=D->ips2;
  int hp=D->has_pcov, hc=D->has_ccov; const double *cp=D->cov_p, *cc=D->cov_c;
  const int *coh=D->coh;
  double *e=S->e,*Wt=S->Wt,*ymnp=S->ymnp,*g=S->g,*H=S->H,*HZ=S->HZ;
  double *rsy=S->rsy,*csy=S->csy,*cgy=S->cgy,*Kx=S->Kx;
  double one=1.0, zero=0.0; int ione=1, info;

  for (int i=0;i<I;i++) rsy[i]=0; for (int j=0;j<J;j++) csy[j]=0; for (int k=0;k<K;k++) cgy[k]=0;
  double lp_lik=0.0;
  for (int j=0;j<J;j++) for (int i=0;i<I;i++) {
    int t=i+I*j; int k=coh[t]-1;
    double ev = bv[ia];
    if (D->has_a) ev += bv[ith+i];
    if (D->has_p) ev += (hp?bv[iph+j]*cp[j]:bv[iph+j]);
    if (D->has_c) ev += (hc?bv[ips+k]*cc[k]:bv[ips+k]);
    if (D->het_a) ev += bv[ith2+i];
    if (D->het_p) ev += bv[iph2+j];
    if (D->het_c) ev += bv[ips2+k];
    if (has_od)   ev += delta[t];
    e[t]=ev;
    double pp = 1.0/(1.0+exp(-ev));
    Wt[t] = N[t]*pp*(1.0-pp);
    double r = Y[t]-N[t]*pp; ymnp[t]=r;
    rsy[i]+=r; csy[j]+=r; cgy[k]+=r;
    lp_lik += Y[t]*ev - N[t]*softplus(ev);
  }
  /* gradient g */
  for (int p=0;p<P;p++) g[p]=0.0;
  { double s=0; for(int t=0;t<IJ;t++) s+=ymnp[t]; g[ia]=s; }
  if (D->has_a){ F77_CALL(dgemv)("N",&I,&I,&one,D->saK,&I,bv+ith,&ione,&zero,Kx,&ione FCONE);
                 for(int i=0;i<I;i++) g[ith+i]=rsy[i]-kappa*Kx[i]; }
  if (D->has_p){ F77_CALL(dgemv)("N",&J,&J,&one,D->spK,&J,bv+iph,&ione,&zero,Kx,&ione FCONE);
                 for(int j=0;j<J;j++) g[iph+j]=(hp?csy[j]*cp[j]:csy[j])-lambda*Kx[j]; }
  if (D->has_c){ F77_CALL(dgemv)("N",&K,&K,&one,D->scK,&K,bv+ips,&ione,&zero,Kx,&ione FCONE);
                 for(int k=0;k<K;k++) g[ips+k]=(hc?cgy[k]*cc[k]:cgy[k])-ny*Kx[k]; }
  if (D->het_a) for(int i=0;i<I;i++) g[ith2+i]=rsy[i]-kappa2*bv[ith2+i];
  if (D->het_p) for(int j=0;j<J;j++) g[iph2+j]=(hp?csy[j]*cp[j]:csy[j])-lambda2*bv[iph2+j];
  if (D->het_c) for(int k=0;k<K;k++) g[ips2+k]=(hc?cgy[k]*cc[k]:cgy[k])-ny2*bv[ips2+k];
  /* Hessian H = assemble_prec(Wt) + ridge */
  memset(H,0,(size_t)P*P*sizeof(double));
  assemble_prec_c(H,D,Wt,kappa,lambda,ny,kappa2,lambda2,ny2);
  { double dsum=0; for(int p=0;p<P;p++) dsum+=H[p+(size_t)P*p];
    double r=1e-6*dsum/P; for(int p=0;p<P;p++) H[p+(size_t)P*p]+=r; }
  /* HZ = H Z (P x Pf); Hg = Z' HZ (Pf x Pf) into Rout */
  F77_CALL(dgemm)("N","N",&P,&Pf,&P,&one,H,&P,Z,&P,&zero,HZ,&P FCONE FCONE);
  F77_CALL(dgemm)("T","N",&Pf,&Pf,&P,&one,Z,&P,HZ,&P,&zero,Rout,&Pf FCONE FCONE);
  /* gg = Z' g (Pf); gamma = Z' bv (Pf) */
  double *gg = S->gg;
  F77_CALL(dgemv)("T",&P,&Pf,&one,Z,&P,g,&ione,&zero,gg,&ione FCONE);
  F77_CALL(dgemv)("T",&P,&Pf,&one,Z,&P,bv,&ione,&zero,gamma,&ione FCONE);
  /* chol(Hg) -> Rout (upper) */
  F77_CALL(dpotrf)("U",&Pf,Rout,&Pf,&info FCONE);
  if (info!=0) return info;
  /* mean_g = gamma + Hg^{-1} gg */
  for (int f=0;f<Pf;f++) mean_g[f]=gg[f];
  { int oneb=1;
    F77_CALL(dtrsv)("U","T","N",&Pf,Rout,&Pf,mean_g,&oneb FCONE FCONE FCONE);
    F77_CALL(dtrsv)("U","N","N",&Pf,Rout,&Pf,mean_g,&oneb FCONE FCONE FCONE); }
  for (int f=0;f<Pf;f++) mean_g[f]+=gamma[f];
  /* lp = lik - 0.5*(prior quadratic forms) */
  double pq=0.0;
  if (D->has_a) pq += kappa  * quad_form(D->saK, bv+ith, I);
  if (D->has_p) pq += lambda * quad_form(D->spK, bv+iph, J);
  if (D->has_c) pq += ny     * quad_form(D->scK, bv+ips, K);
  if (D->het_a){ double q=0; for(int i=0;i<I;i++) q+=bv[ith2+i]*bv[ith2+i]; pq+=kappa2*q; }
  if (D->het_p){ double q=0; for(int j=0;j<J;j++) q+=bv[iph2+j]*bv[iph2+j]; pq+=lambda2*q; }
  if (D->het_c){ double q=0; for(int k=0;k<K;k++) q+=bv[ips2+k]*bv[ips2+k]; pq+=ny2*q; }
  *lp_out = lp_lik - 0.5*pq;
  return 0;
}

/* logq(R,x,m) = sum(log diag R) - 0.5 ||R (x-m)||^2 ; R upper-tri Pf x Pf */
static double logq_c(const double *R, const double *x, const double *m, int Pf, double *tmp) {
  double ld=0.0; for (int f=0;f<Pf;f++) ld += log(R[f+(size_t)Pf*f]);
  for (int f=0;f<Pf;f++) tmp[f]=x[f]-m[f];
  int ione=1; F77_CALL(dtrmv)("U","N","N",&Pf,R,&Pf,tmp,&ione FCONE FCONE FCONE);
  double ss=0.0; for (int f=0;f<Pf;f++) ss+=tmp[f]*tmp[f];
  return ld - 0.5*ss;
}

/* ===========================================================================
 * pg_chain_c : full MCMC chain. Single SEXP `args` list carries everything
 * (mirrors the .bamp_pg_chain signature). Returns the same named list as the
 * R engine. The whole-chain RNG order is locked in the file header.
 * ===========================================================================*/
SEXP pg_chain_c(SEXP args) {
  apc_dims D; fill_dims(&D, args);
  int I=D.I, J=D.J, K=D.K, P=D.P, IJ=I*J;
  int ia=D.ia, ith=D.ith, iph=D.iph, ips=D.ips, ith2=D.ith2, iph2=D.iph2, ips2=D.ips2;
  int has_a=D.has_a, has_p=D.has_p, has_c=D.has_c;
  int het_a=D.het_a, het_p=D.het_p, het_c=D.het_c;
  int hp=D.has_pcov, hc=D.has_ccov;
  const double *cp=D.cov_p, *cc=D.cov_c;
  const int *coh = D.coh;                        /* 1-based cohort idx per cell */

  const double *Y = REAL(getv(args,"Y"));        /* I x J */
  const double *N = REAL(getv(args,"N"));
  const double *ymh = REAL(getv(args,"ymh"));    /* Y - N/2, precomputed in R */
  int n_iter = geti(args,"n_iter"), burn_in = geti(args,"burn_in"), thin = geti(args,"thin");
  int has_od = geti(args,"overdisp");
  double zha = getd(args,"z_hyper_a"), zhb = getd(args,"z_hyper_b");
  /* structure-matrix hyperparameters a,b and ranks (precomputed in R, scaled K) */
  double sa_a=getd(args,"sa_a"), sa_b=getd(args,"sa_b"); int sa_rank=geti(args,"sa_rank");
  double sp_a=getd(args,"sp_a"), sp_b=getd(args,"sp_b"); int sp_rank=geti(args,"sp_rank");
  double sc_a=getd(args,"sc_a"), sc_b=getd(args,"sc_b"); int sc_rank=geti(args,"sc_rank");
  double k2a_a=getd(args,"k2a_a"), k2a_b=getd(args,"k2a_b");
  double k2p_a=getd(args,"k2p_a"), k2p_b=getd(args,"k2p_b");
  double k2c_a=getd(args,"k2c_a"), k2c_b=getd(args,"k2c_b");
  const double *saK=D.saK, *spK=D.spK, *scK=D.scK;
  /* constraint matrix A (k x P col-major) and Zbasis (P x (P-k) col-major) */
  SEXP A_ = getv(args,"A"); SEXP Z_ = getv(args,"Zbasis");
  const double *A = (A_==R_NilValue)?NULL:REAL(A_);
  int kcon = (A_==R_NilValue)?0:INTEGER(getAttrib(A_,R_DimSymbol))[0];
  const double *Z = REAL(Z_);
  int Pf = INTEGER(getAttrib(Z_,R_DimSymbol))[1];   /* free dim = P - kcon */
  /* initial values (precomputed in R: mu, theta, phi from empirical init) */
  double mu = getd(args,"mu0");
  double *theta=(double*)R_alloc(I,sizeof(double));
  double *phi  =(double*)R_alloc(J,sizeof(double));
  double *psi  =(double*)R_alloc(K,sizeof(double));
  { SEXP t0=getv(args,"theta0"),p0=getv(args,"phi0");
    for(int i=0;i<I;i++) theta[i]= has_a? REAL(t0)[i]:0.0;
    for(int j=0;j<J;j++) phi[j]  = has_p? REAL(p0)[j]:0.0;
    for(int k=0;k<K;k++) psi[k]  = 0.0; }
  double kappa=1,lambda=1,ny=1, kappa2=1,lambda2=1,ny2=1;
  double zeta = zha/zhb;
  double *theta2=(double*)R_alloc(I,sizeof(double));
  double *phi2  =(double*)R_alloc(J,sizeof(double));
  double *psi2  =(double*)R_alloc(K,sizeof(double));
  double *delta =(double*)R_alloc(IJ,sizeof(double));    /* I x J cell effect */
  for(int i=0;i<I;i++) theta2[i]=0; for(int j=0;j<J;j++) phi2[j]=0;
  for(int k=0;k<K;k++) psi2[k]=0;   for(int t=0;t<IJ;t++) delta[t]=0;

  /* output storage: nkeep rows. nkeep MUST equal the number of stored sweeps,
   * i.e. #{it in (burn_in, n_iter] : (it-burn_in) % thin == 0} = floor((n_iter-burn_in)/thin),
   * matching the store grid at the `it>burn_in && (it-burn_in)%thin==0` test below.
   * (The old seq(burn_in+1,..,by=thin) count was one too large when thin did not divide
   * n_iter-burn_in, leaving the final row unwritten.) */
  int nkeep = (n_iter - burn_in) / thin;
  SEXP out_theta=PROTECT(allocMatrix(REALSXP,nkeep,I));
  SEXP out_phi  =PROTECT(allocMatrix(REALSXP,nkeep,J));
  SEXP out_psi  =PROTECT(allocMatrix(REALSXP,nkeep,K));
  SEXP out_t2=PROTECT(allocMatrix(REALSXP,nkeep,I));
  SEXP out_p2=PROTECT(allocMatrix(REALSXP,nkeep,J));
  SEXP out_s2=PROTECT(allocMatrix(REALSXP,nkeep,K));
  SEXP out_kap=PROTECT(allocVector(REALSXP,nkeep)), out_lam=PROTECT(allocVector(REALSXP,nkeep));
  SEXP out_ny=PROTECT(allocVector(REALSXP,nkeep)), out_mu=PROTECT(allocVector(REALSXP,nkeep));
  SEXP out_zeta=PROTECT(allocVector(REALSXP,nkeep)), out_dev=PROTECT(allocVector(REALSXP,nkeep));
  SEXP out_k2=PROTECT(allocVector(REALSXP,nkeep)), out_l2=PROTECT(allocVector(REALSXP,nkeep));
  SEXP out_n2=PROTECT(allocVector(REALSXP,nkeep));
  int nprot = 15;
  double *oth=REAL(out_theta),*oph=REAL(out_phi),*ops=REAL(out_psi);
  double *ot2=REAL(out_t2),*op2=REAL(out_p2),*os2=REAL(out_s2);
  double *okap=REAL(out_kap),*olam=REAL(out_lam),*ony=REAL(out_ny),*omu=REAL(out_mu);
  double *ozeta=REAL(out_zeta),*odev=REAL(out_dev);
  double *ok2=REAL(out_k2),*ol2=REAL(out_l2),*on2=REAL(out_n2);
  /* zero the output storage (allocMatrix/allocVector do NOT initialise REALSXP),
   * matching the R reference's matrix(0,...)/numeric(): defence-in-depth so any
   * row the store loop does not reach is a defined 0 rather than heap garbage. */
  memset(oth,0,(size_t)nkeep*I*sizeof(double)); memset(ot2,0,(size_t)nkeep*I*sizeof(double));
  memset(oph,0,(size_t)nkeep*J*sizeof(double)); memset(op2,0,(size_t)nkeep*J*sizeof(double));
  memset(ops,0,(size_t)nkeep*K*sizeof(double)); memset(os2,0,(size_t)nkeep*K*sizeof(double));
  memset(okap,0,(size_t)nkeep*sizeof(double));  memset(olam,0,(size_t)nkeep*sizeof(double));
  memset(ony,0,(size_t)nkeep*sizeof(double));   memset(omu,0,(size_t)nkeep*sizeof(double));
  memset(ozeta,0,(size_t)nkeep*sizeof(double)); memset(odev,0,(size_t)nkeep*sizeof(double));
  memset(ok2,0,(size_t)nkeep*sizeof(double));   memset(ol2,0,(size_t)nkeep*sizeof(double));
  memset(on2,0,(size_t)nkeep*sizeof(double));
  double *ksi_sum=(double*)R_alloc(IJ,sizeof(double)); for(int t=0;t<IJ;t++) ksi_sum[t]=0; int ksi_n=0;

  /* scratch */
  double *eta=(double*)R_alloc(IJ,sizeof(double));
  double *omega=(double*)R_alloc(IJ,sizeof(double));
  double *wgt=(double*)R_alloc(IJ,sizeof(double));
  double *bvec=(double*)R_alloc(P,sizeof(double));
  double *Q=(double*)R_alloc((size_t)P*P,sizeof(double));
  double *U=(double*)R_alloc((size_t)P*P,sizeof(double));
  double *beta=(double*)R_alloc(P,sizeof(double));
  double *pgb=(double*)R_alloc(IJ,sizeof(double));   /* b arg to pg_rpg (=N) */
  int one=1; int info;
  /* ASIS scratch */
  double *zwork=(double*)R_alloc(IJ,sizeof(double));
  double *etaA=(double*)R_alloc(IJ,sizeof(double));  /* eta used by ASIS, updated in place */
  /* Laplace-MH scratch (allocated once) */
  mh_scratch S;
  S.e=(double*)R_alloc(IJ,sizeof(double)); S.Wt=(double*)R_alloc(IJ,sizeof(double));
  S.ymnp=(double*)R_alloc(IJ,sizeof(double)); S.g=(double*)R_alloc(P,sizeof(double));
  S.H=(double*)R_alloc((size_t)P*P,sizeof(double)); S.HZ=(double*)R_alloc((size_t)P*Pf,sizeof(double));
  S.rsy=(double*)R_alloc(I,sizeof(double)); S.csy=(double*)R_alloc(J,sizeof(double));
  S.cgy=(double*)R_alloc(K,sizeof(double)); S.Kx=(double*)R_alloc(P,sizeof(double));
  S.gg=(double*)R_alloc(Pf,sizeof(double));
  double *curR=(double*)R_alloc((size_t)Pf*Pf,sizeof(double));
  double *propR=(double*)R_alloc((size_t)Pf*Pf,sizeof(double));
  double *cur_mg=(double*)R_alloc(Pf,sizeof(double)), *cur_gam=(double*)R_alloc(Pf,sizeof(double));
  double *prop_mg=(double*)R_alloc(Pf,sizeof(double)), *prop_gam=(double*)R_alloc(Pf,sizeof(double));
  double *bcur=(double*)R_alloc(P,sizeof(double)), *bstar=(double*)R_alloc(P,sizeof(double));
  double *gstar=(double*)R_alloc(Pf,sizeof(double)), *qtmp=(double*)R_alloc(Pf,sizeof(double));
  double *zP=(double*)R_alloc(P,sizeof(double));

  /* helper: contribution of period coef j -> pe, cohort -> ce */
  #define PE(j,val) (hp ? (val)*cp[j] : (val))
  #define CE(k,val) (hc ? (val)*cc[k] : (val))

  GetRNGstate();
  int keep=0;
  for (int it=1; it<=n_iter; it++) {
    /* ---- eta (col-major i fastest) and omega ---- */
    for (int j=0;j<J;j++) for (int i=0;i<I;i++) {
      int t=i+I*j; int k=coh[t]-1;
      double e = mu + theta[i] + PE(j,phi[j]) + CE(k,psi[k]);
      if (het_a) e += theta2[i];
      if (het_p) e += phi2[j];
      if (het_c) e += psi2[k];
      if (has_od) e += delta[t];
      eta[t]=e;
    }
    /* omega ~ PG(N, eta): IJ norm_rand in col-major order (== R as.numeric) */
    for (int t=0;t<IJ;t++) pgb[t]=N[t];
    pg_rpg_c(IJ, pgb, eta, omega);

    /* ---- build b and working weight wgt (collapsed-delta if overdisp) ---- */
    for (int p=0;p<P;p++) bvec[p]=0.0;
    if (has_od) {
      for (int t=0;t<IJ;t++) wgt[t]= omega[t]*zeta/(omega[t]+zeta);
      /* bm = ymh*zeta/(omega+zeta) */
      double bmu=0;
      for (int j=0;j<J;j++) for (int i=0;i<I;i++){ int t=i+I*j; double bm=ymh[t]*zeta/(omega[t]+zeta);
        bmu+=bm;
        if(has_a) bvec[ith+i]+=bm;
        if(has_p) bvec[iph+j]+=bm;
        if(has_c) bvec[ips+coh[t]-1]+=bm; }
      bvec[ia]=bmu;
    } else {
      for (int t=0;t<IJ;t++) wgt[t]=omega[t];
      double bmu=0;
      for (int j=0;j<J;j++) for (int i=0;i<I;i++){ int t=i+I*j; double bm=ymh[t];
        bmu+=bm;
        if(has_a) bvec[ith+i]+=bm;
        if(has_p) bvec[iph+j]+=bm;
        if(has_c) bvec[ips+coh[t]-1]+=bm; }
      bvec[ia]=bmu;
    }
    if (hp) for(int j=0;j<J;j++) bvec[iph+j]*=cp[j];
    if (hc) for(int k=0;k<K;k++) bvec[ips+k]*=cc[k];
    if (het_a) for(int i=0;i<I;i++) bvec[ith2+i]=bvec[ith+i];
    if (het_p) for(int j=0;j<J;j++) bvec[iph2+j]=bvec[iph+j];
    if (het_c) for(int k=0;k<K;k++) bvec[ips2+k]=bvec[ips+k];

    /* ---- assemble Q + ridge, Cholesky, constrained draw ---- */
    memset(Q,0,(size_t)P*P*sizeof(double));
    assemble_prec_c(Q,&D,wgt,kappa,lambda,ny,kappa2,lambda2,ny2);
    { double dsum=0; for(int p=0;p<P;p++) dsum+=Q[p+(size_t)P*p];
      double r=1e-6*dsum/P; for(int p=0;p<P;p++) Q[p+(size_t)P*p]+=r; }
    memcpy(U,Q,(size_t)P*P*sizeof(double));
    F77_CALL(dpotrf)("U",&P,U,&P,&info FCONE);
    if(info!=0){ PutRNGstate(); UNPROTECT(nprot); error("dpotrf (Gibbs) info=%d",info); }
    /* mean = Q^{-1} b */
    for(int p=0;p<P;p++) beta[p]=bvec[p];
    chol_solve(U,P,beta);
    /* + U^{-1} z : P norm_rand (zP preallocated) */
    for(int p=0;p<P;p++) zP[p]=norm_rand();
    F77_CALL(dtrsv)("U","N","N",&P,U,&P,zP,&one FCONE FCONE FCONE);
    for(int p=0;p<P;p++) beta[p]+=zP[p];
    /* constraint projection: beta -= Q^{-1}A'(A Q^{-1}A')^{-1} A beta */
    if (kcon>0) {
      double *W=(double*)R_alloc((size_t)P*kcon,sizeof(double));
      for(int c=0;c<kcon;c++){ double *wc=W+(size_t)P*c;
        for(int p=0;p<P;p++) wc[p]=A[c+(size_t)kcon*p];
        chol_solve(U,P,wc); }
      double *Ax=(double*)R_alloc(kcon,sizeof(double));
      double *AW=(double*)R_alloc((size_t)kcon*kcon,sizeof(double));
      for(int c=0;c<kcon;c++){ double s=0; for(int p=0;p<P;p++) s+=A[c+(size_t)kcon*p]*beta[p]; Ax[c]=s;
        for(int d=0;d<kcon;d++){ double sd=0; const double *wd=W+(size_t)P*d;
          for(int p=0;p<P;p++) sd+=A[c+(size_t)kcon*p]*wd[p]; AW[c+(size_t)kcon*d]=sd; } }
      F77_CALL(dposv)("U",&kcon,&one,AW,&kcon,Ax,&kcon,&info FCONE);
      if(info!=0){ PutRNGstate(); UNPROTECT(nprot); error("dposv info=%d",info); }
      for(int p=0;p<P;p++){ double s=0; for(int c=0;c<kcon;c++) s+=W[p+(size_t)P*c]*Ax[c]; beta[p]-=s; }
    }
    mu=beta[ia];
    if(has_a) for(int i=0;i<I;i++) theta[i]=beta[ith+i];
    if(has_p) for(int j=0;j<J;j++) phi[j]=beta[iph+j];
    if(has_c) for(int k=0;k<K;k++) psi[k]=beta[ips+k];
    if(het_a) for(int i=0;i<I;i++) theta2[i]=beta[ith2+i];
    if(het_p) for(int j=0;j<J;j++) phi2[j]=beta[iph2+j];
    if(het_c) for(int k=0;k<K;k++) psi2[k]=beta[ips2+k];

    /* ---- overdispersion: delta_ij | rest ~ N, then zeta ~ Gamma ---- */
    if (has_od) {
      /* eta0 = smooth + het (no delta) */
      double ss=0;
      for (int j=0;j<J;j++) for (int i=0;i<I;i++){ int t=i+I*j; int k=coh[t]-1;
        double e0=mu+theta[i]+PE(j,phi[j])+CE(k,psi[k]);
        if(het_a)e0+=theta2[i]; if(het_p)e0+=phi2[j]; if(het_c)e0+=psi2[k];
        double precd=omega[t]+zeta;
        double mden=(ymh[t]-omega[t]*e0)/precd;
        delta[t]=mden + norm_rand()/sqrt(precd);   /* IJ norm_rand, col-major */
        ss+=delta[t]*delta[t];
      }
      zeta = rgamma(zha+0.5*IJ, 1.0/(zhb+0.5*ss));  /* R rate -> C scale */
    }

    /* ---- precisions: Gamma full conditionals (quadratic forms via saK etc) ---- */
    if (has_a){ double q=quad_form(saK,theta,I); kappa =rgamma(sa_a+sa_rank/2.0, 1.0/(sa_b+0.5*q)); }
    if (has_p){ double q=quad_form(spK,phi,J);   lambda=rgamma(sp_a+sp_rank/2.0, 1.0/(sp_b+0.5*q)); }
    if (has_c){ double q=quad_form(scK,psi,K);   ny    =rgamma(sc_a+sc_rank/2.0, 1.0/(sc_b+0.5*q)); }
    if (het_a){ double q=0; for(int i=0;i<I;i++) q+=theta2[i]*theta2[i]; kappa2 =rgamma(k2a_a+0.5*(I-1), 1.0/(k2a_b+0.5*q)); }
    if (het_p){ double q=0; for(int j=0;j<J;j++) q+=phi2[j]*phi2[j];     lambda2=rgamma(k2p_a+0.5*(J-1), 1.0/(k2p_b+0.5*q)); }
    if (het_c){ double q=0; for(int k=0;k<K;k++) q+=psi2[k]*psi2[k];     ny2    =rgamma(k2c_a+0.5*(K-1), 1.0/(k2c_b+0.5*q)); }

    /* ---- ASIS interweaving: non-centred re-draw of each precision ----
     * Per present effect in order a,p,c: propk = exp(log(prec)+rnorm(0,0.4)),
     * accept via runif. ll(xn) = -0.5 sum omega*(eta + dEta(xn-x) - zwork)^2.
     * etaA is the current full predictor and is updated in place on accept. */
    for (int j=0;j<J;j++) for (int i=0;i<I;i++){ int t=i+I*j; int k=coh[t]-1;
      double e=mu+theta[i]+PE(j,phi[j])+CE(k,psi[k]);
      if(het_a)e+=theta2[i]; if(het_p)e+=phi2[j]; if(het_c)e+=psi2[k];
      if(has_od)e+=delta[t];
      etaA[t]=e; zwork[t]=ymh[t]/omega[t];
    }
    /* helper macro: current ll given etaA */
    #define ASIS_EFFECT(present, LEN, XVEC, PREC, AA, BB, DETA_CODE)               \
      if (present) {                                                              \
        int L=LEN; double prec=PREC;                                              \
        double propk=exp(log(prec)+norm_rand()*0.4);                              \
        double sp_=sqrt(prec), spp=sqrt(propk);                                   \
        /* xp = xt/sqrt(propk) = x*sqrt(prec)/sqrt(propk); d = xp - x */          \
        double scale=sp_/spp;                                                     \
        double ll_cur=0.0, ll_prop=0.0;                                           \
        /* compute ll at current (d=0) and at proposal (d=(scale-1)*x) */         \
        for (int t2=0;t2<IJ;t2++){ double diff=etaA[t2]-zwork[t2]; ll_cur+=omega[t2]*diff*diff; } \
        ll_cur*=-0.5;                                                             \
        { DETA_CODE /* fills tmpdeta[t] = make_deta(xp-x) into etaA-add space */ }\
        for (int t2=0;t2<IJ;t2++){ double diff=etaA[t2]+tmpdeta[t2]-zwork[t2]; ll_prop+=omega[t2]*diff*diff; } \
        ll_prop*=-0.5;                                                            \
        double lp_cur = ll_cur + AA*log(prec) - BB*prec;                          \
        double lp_prop= ll_prop+ AA*log(propk)- BB*propk;                         \
        if (log(unif_rand()) < lp_prop - lp_cur) {                               \
          for (int q=0;q<L;q++) XVEC[q]*=scale;                                   \
          for (int t2=0;t2<IJ;t2++) etaA[t2]+=tmpdeta[t2];                        \
          PREC=propk;                                                            \
        }                                                                         \
      }
    {
      double *tmpdeta=zP;  /* reuse: length>=IJ? zP is length P; need IJ. use S.e */
      tmpdeta=S.e;         /* S.e is length IJ, free here (state() not yet called) */
      /* age: dEta(d_i) added to row i of every column -> tmpdeta[i,j]=d_i */
      ASIS_EFFECT(has_a, I, theta, kappa, sa_a, sa_b, {
        for(int j=0;j<J;j++) for(int i=0;i<I;i++){ tmpdeta[i+I*j]=(scale-1.0)*theta[i]; }
      })
      ASIS_EFFECT(has_p, J, phi, lambda, sp_a, sp_b, {
        for(int j=0;j<J;j++){ double d=PE(j,(scale-1.0)*phi[j]); for(int i=0;i<I;i++) tmpdeta[i+I*j]=d; }
      })
      ASIS_EFFECT(has_c, K, psi, ny, sc_a, sc_b, {
        for(int j=0;j<J;j++) for(int i=0;i<I;i++){ int k=coh[i+I*j]-1; tmpdeta[i+I*j]=CE(k,(scale-1.0)*psi[k]); }
      })
    }
    #undef ASIS_EFFECT

    /* ---- Laplace (Newton) Metropolis-Hastings refinement ---- */
    for (int p=0;p<P;p++) bcur[p]=0.0;
    bcur[ia]=mu;
    if(has_a) for(int i=0;i<I;i++) bcur[ith+i]=theta[i];
    if(has_p) for(int j=0;j<J;j++) bcur[iph+j]=phi[j];
    if(has_c) for(int k=0;k<K;k++) bcur[ips+k]=psi[k];
    if(het_a) for(int i=0;i<I;i++) bcur[ith2+i]=theta2[i];
    if(het_p) for(int j=0;j<J;j++) bcur[iph2+j]=phi2[j];
    if(het_c) for(int k=0;k<K;k++) bcur[ips2+k]=psi2[k];
    double cur_lp, prop_lp;
    /* cur state (no RNG); chol failure errors, matching R's chol(Hg) */
    if (mh_state(&D,Y,N,Z,Pf,bcur,kappa,lambda,ny,kappa2,lambda2,ny2,delta,has_od,
                 &S,curR,cur_mg,cur_gam,&cur_lp)!=0) {
      PutRNGstate(); UNPROTECT(nprot); error("MH cur-state Cholesky failed"); }
    /* gstar = cur_mg + curR^{-1} z   (Pf norm_rand, ALWAYS); bstar = Z gstar */
    for (int f=0;f<Pf;f++) gstar[f]=norm_rand();
    F77_CALL(dtrsv)("U","N","N",&Pf,curR,&Pf,gstar,&one FCONE FCONE FCONE);
    for (int f=0;f<Pf;f++) gstar[f]+=cur_mg[f];
    { double zero=0.0, oned=1.0; int ione=1;
      F77_CALL(dgemv)("N",&P,&Pf,&oned,Z,&P,gstar,&ione,&zero,bstar,&ione FCONE); }
    if (mh_state(&D,Y,N,Z,Pf,bstar,kappa,lambda,ny,kappa2,lambda2,ny2,delta,has_od,
                 &S,propR,prop_mg,prop_gam,&prop_lp)!=0) {
      PutRNGstate(); UNPROTECT(nprot); error("MH prop-state Cholesky failed"); }
    double la = prop_lp - cur_lp
              + logq_c(propR, cur_gam, prop_mg, Pf, qtmp)
              - logq_c(curR,  gstar,   cur_mg,  Pf, qtmp);
    /* runif drawn ONLY when la is finite (R's short-circuit &&) */
    if (R_FINITE(la) && log(unif_rand()) < la) {
      mu=bstar[ia];
      if(has_a) for(int i=0;i<I;i++) theta[i]=bstar[ith+i];
      if(has_p) for(int j=0;j<J;j++) phi[j]=bstar[iph+j];
      if(has_c) for(int k=0;k<K;k++) psi[k]=bstar[ips+k];
      if(het_a) for(int i=0;i<I;i++) theta2[i]=bstar[ith2+i];
      if(het_p) for(int j=0;j<J;j++) phi2[j]=bstar[iph2+j];
      if(het_c) for(int k=0;k<K;k++) psi2[k]=bstar[ips2+k];
    }

    /* ---- store ---- */
    if (it>burn_in && ((it-burn_in)%thin==0)) {
      for(int i=0;i<I;i++) oth[keep+nkeep*i]=theta[i];
      for(int j=0;j<J;j++) oph[keep+nkeep*j]=phi[j]*( hp?cp[j]:1.0 );
      for(int k=0;k<K;k++) ops[keep+nkeep*k]=psi[k]*( hc?cc[k]:1.0 );
      for(int i=0;i<I;i++) ot2[keep+nkeep*i]=theta2[i];
      for(int j=0;j<J;j++) op2[keep+nkeep*j]=phi2[j];
      for(int k=0;k<K;k++) os2[keep+nkeep*k]=psi2[k];
      okap[keep]=kappa; olam[keep]=lambda; ony[keep]=ny; omu[keep]=mu;
      ozeta[keep]=zeta; ok2[keep]=kappa2; ol2[keep]=lambda2; on2[keep]=ny2;
      /* fitted predictor + deviance (absolute/covariate-scaled eta) */
      double dev=0;
      for (int j=0;j<J;j++) for (int i=0;i<I;i++){ int t=i+I*j; int k=coh[t]-1;
        double e=mu+theta[i]+PE(j,phi[j])+CE(k,psi[k]);
        if(het_a)e+=theta2[i]; if(het_p)e+=phi2[j]; if(het_c)e+=psi2[k];
        if(has_od)e+=delta[t];
        ksi_sum[t]+=e;
        double pr=1.0/(1.0+exp(-e)), yhat=N[t]*pr;
        double d2;
        if (Y[t]>0) d2=2.0*(Y[t]*log(Y[t]/yhat)+(N[t]-Y[t])*log((N[t]-Y[t])/(N[t]-yhat)));
        else        d2=2.0*((N[t]-Y[t])*log((N[t]-Y[t])/(N[t]-yhat)));
        if (!R_FINITE(d2)) d2=0.0;
        dev+=d2;
      }
      ksi_n++;
      odev[keep]=dev;
      keep++;
    }
  }
  PutRNGstate();

  /* ksi = ksi_sum / ksi_n, returned as I x J */
  SEXP ksi_=PROTECT(allocMatrix(REALSXP,I,J)); nprot++;
  for(int t=0;t<IJ;t++) REAL(ksi_)[t]=ksi_sum[t]/ksi_n;

  /* assemble return list */
  const char *nms[]={"theta","phi","psi","kappa","lambda","ny","my","zeta",
                     "theta2","phi2","psi2","kappa2","lambda2","ny2","deviance","ksi",""};
  SEXP res=PROTECT(mkNamed(VECSXP,nms)); nprot++;
  SET_VECTOR_ELT(res,0,out_theta); SET_VECTOR_ELT(res,1,out_phi); SET_VECTOR_ELT(res,2,out_psi);
  SET_VECTOR_ELT(res,3,out_kap); SET_VECTOR_ELT(res,4,out_lam); SET_VECTOR_ELT(res,5,out_ny);
  SET_VECTOR_ELT(res,6,out_mu); SET_VECTOR_ELT(res,7,out_zeta);
  SET_VECTOR_ELT(res,8,out_t2); SET_VECTOR_ELT(res,9,out_p2); SET_VECTOR_ELT(res,10,out_s2);
  SET_VECTOR_ELT(res,11,out_k2); SET_VECTOR_ELT(res,12,out_l2); SET_VECTOR_ELT(res,13,out_n2);
  SET_VECTOR_ELT(res,14,out_dev); SET_VECTOR_ELT(res,15,ksi_);
  UNPROTECT(nprot);
  return res;
}
