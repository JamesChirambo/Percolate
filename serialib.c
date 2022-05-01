#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>

#include "percolate.h"



void mpstart(int *size, int *rank) 
{
  *size   = 1;
  *rank   = 0;
}


double gettime(void) 
{
    struct timeval tp;
    gettimeofday (&tp, NULL);
    return tp.tv_sec + tp.tv_usec/(double)1.0e6;
}

void mpstop(void) 
{
    
}

void fnd2dnbrs(int size, int *rank, int *left, int *right, int *top, int *bottom) 
{
  *left = *right = *top = *bottom = -1;
}

void fnd2ddecomp(int l, int *istart, int *istop, int *jstart, int *jstop ) 
{
  *istart = *jstart = 1;
  *istop = *jstop = l;
}


void bcast(int *ival) 
{
}

void globalsum(int *ival) 
{

}


void startbcswap(int **old, int ml, int nl, int left, int right, int top, int bottom, int size) 
{
      /*
    *  Implement periodic boundaries in the first dimension "i",
    *  i.e. map i=0 to i=M and i=M+1 to i=1 (for all values of j)
    */
  int j;

  for (j=1; j <= nl; j++)
	{
	  old[0][j]   = old[ml][j];
	  old[ml+1][j] = old[1][j];
	}

  /*
    *  Map is mostly periodic across the first dimension "i". In
    *  serial, we can do this with a simple copy as above.
    * 
    *  However, it is only periodic for certain elements: in
    *  blocks of "jpbcblock" cells separated by "jpbcstride". For
    *  example, if jpcblock=4 and jpbcstride=6 then boundary
    *  elements 1 and 2 will be non-periodic, 3, 4, 5 and 6
    *  periodic, 7 and 8 non-periodic, 9, 10, 11 and 12 periodic,
    *  etc. etc.
    * 
    *  The simplest approach is to copy across the whole boundary
    *  then discard the unwanted (non-periodic) elements.  It may
    *  seem wasteful to copy the halos then zero some of them, but
    *  this approach is probably the simplest to implement in
    *  parallel.
    */

  int jpbcstride = 10;
  int jpbcblock  =  8;
  int jremain;

  for (j=1; j <= nl; j++)
	{
	  // Is this one of the non-periodic elements?

	  jremain = (j-1)%jpbcstride;

	  if (jremain < (jpbcstride-jpbcblock))
	    {
	      // Periodic boundaries do not operate in this region

	      old[0][j]   = 0;
	      old[ml+1][j] = 0;
	    }
	}


}

void endbcswap(void)
{
}


void gathermlist(int *inbuff, int **outbuff, int ndata) 
{

  //printf("Gather Success!\n");
}

void split2dmap(int **inbuff, int **outbuff, int **alldimstore, int istart, int istop, int jstart, int jstop, int mdata, int ndata, int rank, int size)
{
  initold(inbuff, outbuff, istart, istop, jstart, jstop);

}

void merge2dmap(int **inbuff, int **outbuff, int **alldimstore, int istart, int istop, int jstart, int jstop, int mdata, int ndata, int rank, int size)
{
  updatemap(outbuff, inbuff, istart, istop, jstart, jstop);

}

