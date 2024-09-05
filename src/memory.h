/*
  memory.h

  Macros for memory allocation

*/

#define CALLOC(LEN, TYP)      R_Calloc(LEN, TYP)
#define REALLOC(PTR,LEN,TYP)  R_Realloc(PTR, LEN, TYP)
#define FREE(PTR)             R_Free(PTR)
