#include <R.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* 
         KnetInhom.c

	 This file defines the 'inhomogeneous' version of the K function
	 
	 The actual C code for the functions is
	 given in 'netKcode.h', which is #included here.

	 The macro 'INHOM' is #defined in this file
	 so that 'netKcode.h' defines functions for the
	 inhomogeneous case, with names beginning 'I_'

*/


/* type definitions */
#include "networkdef.h"
/* prototypes of functions from netbase.c */
#include "netbase.h"

/* function definitions */
#define INHOM
#include "netKcode.h"



