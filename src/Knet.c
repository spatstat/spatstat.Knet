#include <R.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* 
         Knet.c

	 This file defines the 'homogeneous' version of the K function
	 
	 The actual C code for the functions is
	 given in 'netKcode.h', which is #included here.

*/

/* interface to memory allocation */
#include "memory.h"

/* type definitions */
#include "networkdef.h"

/* prototypes of functions from netbase.c */
#include "netbase.h"

/* function definitions */
#undef INHOM
#include "netKcode.h"



