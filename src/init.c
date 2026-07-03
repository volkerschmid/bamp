#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* .C calls */
extern void bamp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls (native Polya-Gamma engine, src/pg_engine.c) */
extern SEXP pg_assemble_c(SEXP, SEXP, SEXP, SEXP);
extern SEXP pg_chain_c(SEXP);

R_CMethodDef CEntries[] = {
  {"bamp", (DL_FUNC) &bamp, 25},
  {NULL}
};

R_CallMethodDef CallEntries[] = {
  {"pg_assemble_c", (DL_FUNC) &pg_assemble_c, 4},
  {"pg_chain_c",    (DL_FUNC) &pg_chain_c,    1},
  {NULL}
};

void R_init_bamp(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}
