#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void fittedllar(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void llarinc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"fittedllar", (DL_FUNC) &fittedllar, 10},
    {"llarinc",       (DL_FUNC) &llarinc,       12},
    {NULL, NULL, 0}
};

void R_init_tsDyn(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
