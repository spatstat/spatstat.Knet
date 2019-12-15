/* 
   Native symbol registration table for spatstat.Knet package

*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  Prototype declarations

*/

void createGraphNet(int *, int *, int *, int *, double *, int *,
		    int *, int *, double *, double *, double *,
		    int *, double *, int *, double *);

void I_createGraphNet(int *, int *, int *, int *, double *, int *,
		      double *,
		      int *, int *, double *, double *, double *,
		      int *, double *, int *, double *);

/*
  declared symbol table
*/

static const R_CMethodDef CEntries[] = {
    {"createGraphNet",           (DL_FUNC) &createGraphNet,   15},
    {"I_createGraphNet",         (DL_FUNC) &I_createGraphNet, 16},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {NULL, NULL, 0}
};

/* register */ 

void R_init_spatstat_Knet(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
