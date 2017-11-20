#include <R.h>
#include <R_ext/Rdynload.h>
#include "bamp.h"

R_CMethodDef cMethods[] = {
  {"bampc", (DL_FUNC)&bamp, 8},
  NULL
};

void R_init_bioimagetools(DllInfo *info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

