/**
 * @file percolate.c
 * @author Exam Number B144511
 * @brief Parallel program to test for percolation of a cluster - 2D Decomposition
 * @version 1.1
 * @date 2022-05-01
 * 
 * @copyright Copyright (c) 2022
 * 
 ***************************************************************************
 * 
 * This version uses dynamic array allocation via the arralloc routine.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "percolate.h"
#include "arralloc.h"
#include "mpilib.h"


int main(int argc, char *argv[])
{

   /* Define the main local arrays for the simulation */

  int **old, **new;

  /* Additional array WITHOUT halos (of size LxL) for initialisation and IO. */

  int **map;

  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */

  int i, j, nhole, step, maxstep; 
  int nchange, tchanges, printfreq;
  int itop, ibot, perc;
  int ncluster = 2;
  double r;

  /*
   * First and second dimension of the array (Simulation local size) 
   */
  int ml, nl; 


  /*
   * the array dimensions (optionally supplied at runtime from the command line arguments)
   */

  int l;

  int checkchange; //whether to quite simulation immediately after there are no change made
  int testmode = FALSE; // help me do sensible performance tests using reasonable amounts of CPU time.

  /*
   * collect timing data
   */

  double ostart, ostop, tot, sstart, sstop, tst;
  tst = 0.0;
  

  /*
   * Parallel variables
   */

  int rank, size;
  int left, right, top, bottom; //neighbours variables


  /*
   *  dimensional boundaries
   */

  int istart, istop, jstart, jstop;

  /*
   * keep decomposition in an array
   */
  int count = 5;
  int *dimstore;
  int **alldimstore; //master list of the decomposition


  /*
   * initialise message passing
   */

  mpstart(&size, &rank);
  
  //log start time
  ostart = gettime();

  /** check number of arguments; quit early if not between 2 and 6.... */
  if (argc < 2 || argc > 7) 
  {  
    if (rank == ROOT) 
    {
            printf("ERROR! Usage: percolate <seed> [scale factor] [checkchange] [ncluster] [rho] [testmode].\n");
    }
    mpstop();
    return 0;
  }


  /*
   * Check command line parameters and parse arguments
   * only on the controller process
   */
  if (rank == ROOT) 
  {
    printf("percolate: running on %d process(es)\n", size);


    /**  Set the randum number seed and initialise the generator */
    seed = 7777;
    if(argc > 1){
      if (atoi(argv[1]) > 0 && atoi(argv[1]) <= 10000){
        seed = atoi(argv[1]);
      } 
      else {
        if(ROOT == rank)
        {
          printf("Using Default seed value %d instead of %d.\n", seed, atoi(argv[2]));
        }
      }
    }


    /** Accept new global simulation size l or use default L */
    l = L; 
    if(argc > 2){
      if (atoi(argv[2]) > 0 && atoi(argv[2]) <= 100000){
        l = atoi(argv[2]);
      } 
      else {
        if(ROOT == rank)
        {
          printf("Using Default L value %d instead of %d.\n", l, atoi(argv[2]));
        }
      }
    }

    /**  Set most important value: the rock density rho (between 0 and 1) */
    rho = 0.4040;
    if(argc > 3){ 
      if(atof(argv[3]) > 0.0 && atof(argv[3]) <= 1.0){
          rho = atof(argv[3]);
      }
    }
    
    /** Set number of clusters to be printed */
    if(argc > 4){
      if(atoi(argv[4]) > 0 && atoi(argv[4]) <= 9)
      {
        ncluster = atoi(argv[4]);
      }
    }

    /** Decide whether to stop loop when there's no change */
    checkchange = FALSE;
    if(argc > 5){
      if(atoi(argv[5]) >= 0 && atoi(argv[5]) <= 1)
      {
        checkchange = atoi(argv[5]);
      }
    }

    /** Set test mode */
    if(argc > 6){ 
      if(atof(argv[6]) >= 0 && atof(argv[6]) <= 1){
          testmode = atof(argv[6]);
      }
    }

  }

  // broadcast runtime parameters to other processors
  bcast(&l);
  bcast(&checkchange);
  bcast(&testmode);

  //consistency check and report
  if(l%size != 0){
    if(rank == ROOT){
      printf("WARNING! m=%d does not divide onto %d processes\n", l, size);
    }
    /** 
     * we could terminate here if desired.
     * However, our decomposition allows for inexact multiples of L and size
    */
    //mpstop();
    //return 0;    
  }

  // arrays to keep a master list of ranks and their indices on the big map
  alldimstore = (int **) arralloc(sizeof(int), 2, size, count);
  dimstore = (int *) arralloc(sizeof(int), 1, count);
  if (NULL == alldimstore || NULL == dimstore)
  {
    printf("percolate: array allocation failed\n");
    return 1;
  }

  for(i=0;i<size;i++){
    for(j=0;j<count;j++){
      alldimstore[i][j]=0;
      dimstore[j]=0;
    }
  }


  /*
   *  Find 2D Decomposition:
   *      my neighbours
   *      my start and stop positions in the map
   *      my m x n values (to determine the size) of my section of the map
   */

  fnd2dnbrs(size, &rank, &left, &right, &top, &bottom);  //my neighbours
  //printf("@rank %d: left %d, right %d; top %d, bottom %d\n", rank, left, right, top, bottom);

  /* Compute the decomposition */
  fnd2ddecomp( l, &istart, &istop, &jstart, &jstop);
  //printf("@rank %d: istart %d, istop %d; jstart %d, jstop %d\n", rank, istart, istop, jstart, jstop);

  ml = (istop - istart + 1);
  nl = (jstop - jstart + 1);


  // gather a master list of ranks and their indices on the big map
  dimstore[0]  = rank; 
  dimstore[1]  = istart;  dimstore[2]  = jstart;
  dimstore[3]  = ml;      dimstore[4]  = nl;

  // only gather the list of ranks and their decomposition on the
  // controller process
  gathermlist(dimstore, alldimstore, count);


  /*
   *  Update for a number of steps (depending on the problem size)
   *    the calculation below only helps to get a maxtep such that  
   *    scale = (l/l)^3 
   *    scale * l * maxstep is approx (L*L*5) = 2949120 
   *    not a perfect solution but works
   */

  if (testmode)
  {
    checkchange = FALSE;
    double scalefactor = (double) l/L;
    scalefactor = scalefactor*scalefactor*scalefactor;
    maxstep = (int) ((l/scalefactor)* (5/scalefactor));
    printfreq = (int) 100/scalefactor;
    if(printfreq<=0) printfreq = 1;
  }

  /*
   *  Update for a fixed number of steps, periodically report progress
   */
  else
  {
    maxstep = 5*l; 
    printfreq = 100;
  }


  if(rank == ROOT){
    printf("percolate: L = %d, rho = %f, seed = %d, maxstep = %d, checkchange = %d , testmode = %d", l, rho, seed, maxstep, checkchange, testmode);
  }


    /*
   * Allocate 2D mlxnl / lxl integer arrays dynamically
   */

  old = (int **) arralloc(sizeof(int), 2, ml+2, nl+2);
  new = (int **) arralloc(sizeof(int), 2, ml+2, nl+2);
  map = (int **) arralloc(sizeof(int), 2, l, l);
  if (NULL == old || NULL == new || NULL == map) 
  {
    printf("percolate: array allocation failed\n");
    return 1;
  }

  /*
   *  Initialise map on the controller process only
   */
  if (rank == ROOT) 
  {
    rinit(seed);

    /*
    *  Initialise map with density rho. Zero indicates rock, a positive
    *  value indicates a hole. For the algorithm to work, all the holes
    *  must be initialised with a unique integer
    */

    nhole = 0;

    for (i=0; i < l; i++)
      {
        for (j=0; j < l; j++)
        {
          r=uni();
          if(r < rho)
          {
            map[i][j] = 0;
          }
          else
          {
            nhole++;
            map[i][j] = nhole;
          }
        }
      }
    printf(", actual density = %f\n", 1.0 - ((double) nhole)/((double) l*l) );
 }


  /*
   * Distribute the map from the controller to all other processes 
   *  using 
   *      MPI_Ssend with sendbuf = map and
   *      MPI_Recv with recvbuf = old
   */

  //copy my section of the map
  split2dmap(map, old, alldimstore, istart, istop, jstart, jstop, ml, nl, rank, size);

  /*
   * Zero the halo values.
   */

  inithalos(old, ml, nl);

  step = 1;
  nchange = 1;
  tchanges = 0;

  sstart = gettime(); // log calculation time


  while (step <= maxstep)
  {

    /* 
      * start boundary swaps to the old array (non blocking)
      */
    startbcswap(old, ml, nl, left, right, top, bottom, size);

    
    //start computation
    nchange = 0;

    updatestep(old, new, ml, nl, &nchange);


    /*
      * Compute global number of change on rank 0
      */

    globalsum(&nchange);
    tchanges += nchange;

    /*
      *  Report progress every now and then
      */

    if (step % printfreq == 0)
    {
      if (rank == ROOT) {
        printf("percolate: changes on step %d is %d\n", step, nchange);
      }
    }

    /*
      *  Copy back in preparation for next step, omitting halos
      */


    for (i=1; i<=ml; i++)
    {
      for (j=1; j<=nl; j++)
      {
        old[i][j] = new[i][j];
      }
    }

    /*
      * finalise boundary swaps to the old array
      */
    endbcswap();

    //Quit early if there are zero changes

    if (checkchange) {
      if (nchange == 0) {
        if (rank == ROOT) {
          printf("percolate: changes on step %d is %d\n", step, nchange);
        }
        break;
      }
    }


    step++;
  }

  sstop = gettime();
  tst = sstop - sstart;

  /*
   *  We set a maximum number of steps to ensure the algorithm always
   *  terminates. However, if we hit this limit before the algorithm
   *  has finished then there must have been a problem (e.g. the value
   *  of maxstep is too small)
   */

  if (rank == ROOT) {
    if (nchange != 0)
      {
          printf("percolate: WARNING max steps = %d reached; Calculation time = %f seconds; time per step = %g; nchange != 0; total changes = %d\n", step, tst, tst/step, tchanges);
      }
      else if (maxstep >= step-1)
      {
        printf("percolate: Update: steps = %d reached; Calculation time  = %f seconds; time per step = %g; nchange = 0; total changes = %d.\n", step, tst, tst/step, tchanges);
      }
      
  }

  ostop = gettime();
  tot = ostop - ostart;

  /*
   * h) gather the map back to the controller from all other processes 
   *    using MPI_Send and MPI_Recv with sendbuf = old and recvbuf = map
   *    NOTE: Only send the centre of old, excluding the halos, into map
   */
  merge2dmap(old, map, alldimstore, istart, istop, jstart, jstop, ml, nl, rank, size);

  /*
   *  Controller Only!
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */

  if (rank == ROOT) 
  {

    perc = 0;

    for (itop=0; itop < l; itop++)
      {
        if (map[itop][l-1] > 0)
        {
          for (ibot=0; ibot < l; ibot++)
            {
              if (map[ibot][0] == map[itop][l-1])
              {
                perc = 1;
              }
            }
        }
      }

    if (perc != 0)
    {
      printf("percolate: cluster DOES percolate\n");
    }
    else
    {
      printf("percolate: cluster DOES NOT percolate\n");
    }

    printf("Total Execution time:  %f seconds.\n", tot);


    /*
    *  Write the map to the file "map.pgm", displaying the two
    *  largest clusters. If the last argument here was 3, it would
    *  display the three largest clusters etc. The picture looks
    *  cleanest with only a single cluster, but multiple clusters
    *  are useful for debugging.
    */

    mapwrite("map.pgm", map, l, ncluster);

    printf("\nDone!!\n");


  }

  /*
   * Free the arrays
   */

  free(old);
  free(new);
  free(map);


  mpstop();
  
  return 0;
}
