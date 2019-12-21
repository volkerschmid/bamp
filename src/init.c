#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void bamp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

R_CMethodDef CEntries[] = {
  {"bamp", (DL_FUNC) &bamp, 25},
  {NULL}
};

void R_init_bamp(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}
