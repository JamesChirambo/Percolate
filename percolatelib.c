#include <stdio.h>
#include <string.h>
#include <math.h>
#include "percolate.h"


void initold(int **map, int **old, int istart, int istop, int jstart, int jstop)
{
    int i, j, myi, myj;
    myi = 0;
    for(i = istart; i <= istop; i++)
    {
    myi++;
    myj = 0;
    for(j = jstart; j <= jstop; j++)
    {
        myj++;
        old[myi][myj] = map[i-1][j-1];
    }
    }

}

void updatestep(int **old, int **new, int ml, int nl, int *nchangelocal)
{
  int i, j, oldval, newval;
        for (i=1; i<=ml; i++)
      {
        //#pragma ivdep
        for (j=1; j<=nl; j++)
        {
          oldval = old[i][j];
          newval = oldval;

          /*
          * Set new[i][j] to be the maximum value of old[i][j]
          * and its four nearest neighbours
          */

          if (oldval != 0)
          {
            if (old[i-1][j] > newval) newval = old[i-1][j];
            if (old[i+1][j] > newval) newval = old[i+1][j];
            if (old[i][j-1] > newval) newval = old[i][j-1];
            if (old[i][j+1] > newval) newval = old[i][j+1];

            if (newval != oldval)
              {
                ++*nchangelocal;
              }
          }
          new[i][j] = newval;
        }
      }
}

void updatemap(int **map, int **old, int istart, int istop, int jstart, int jstop)
{
    int i, j, myi, myj;
    myi = 0;
    for(i = istart; i <= istop; i++)
    {
      myi++;
      myj = 0;
      for(j = jstart; j <= jstop; j++)
      {
        myj++;
        map[i-1][j-1] = old[myi][myj];
      }
    }

}

void inithalos(int **old, int ml, int nl)
{
  int i, j;
  for (i=0; i <= ml+1; i++)  // zero the bottom and top halos
  {
    old[i][0]   = 0;
    old[i][nl+1] = 0;
  }

  for (j=0; j <= nl+1; j++)  // zero the left and right halos
  {
    old[0][j]   = 0;
    old[ml+1][j] = 0;
  }
}

void nextrank(int **dstore, int i, int *rank, int *istart, int *jstart, int *m, int *n)
{
      *rank   = dstore[i][0];
      *istart = dstore[i][1];
      *jstart = dstore[i][2];
      *m      = dstore[i][3];
      *n      = dstore[i][4];
}

