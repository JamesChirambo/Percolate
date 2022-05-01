#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "percolate.h"

/** Parallel variables */
static MPI_Comm comm; // to store MPI_COMM_WORLD

/** Cartesian Topology Variables */
static MPI_Comm comm2d;   // Communicator for a decomposition of the domain. 
static int dims[ndims]; 
static int periods[ndims];
static int coords[ndims];


/** For swaping halos */
static MPI_Status status[4];
static MPI_Request request[4];

/**
 * @brief Initialise MPI
 * 
 * @param size  number of processes running the program
 * @param rank  this process's ID
 */
void mpstart(int *size, int *rank) 
{

  comm = MPI_COMM_WORLD;
   
  /** Initialise MPI */
  MPI_Init(NULL, NULL);

  /** compute the size and rank  */
  MPI_Comm_size(comm, size);
  MPI_Comm_rank(comm, rank);
  
}

/**
 * @brief Terminate MPI execution environment
 *
 */
void mpstop(void) 
{
  MPI_Finalize();
}

/**
 * @brief Return time on the calling processor
 *
 * @return double
 */
double gettime(void) 
{
  MPI_Barrier(comm);
  return MPI_Wtime();
  
}


/**
 * @brief Neighbours are computed with a Cartesian topology
 * Note that the special rank of MPI_PROC_NULL is a "black hole" for
 * communications. Using this value for processes off the edges of the
 * image means there is no additional logic needed to ensure processes
 * at the edges do not attempt to send to or receive from invalid
 * ranks (i.e. rank = -1 and rank = NPROC).
 */
void fnd2dnbrs(int size, int *rank, int *left, int *right, int *top, int *bottom) 
{

  int reorder;   //Cartesian Topology Variables

/**  initialise variables  */
  dims[0] = 0;
  periods[0] = TRUE; // the map periodic across the first dimension "i"
  dims[1] = 0;  periods[1] = FALSE;
  reorder = FALSE;

  /** Let MPI find a "good" decomposition */
  MPI_Dims_create(size, ndims, dims);

  /** Create a new communicator to which topology information has been attached  */
  MPI_Cart_create(comm, ndims, dims, periods, reorder, &comm2d);

  /** Get my position in this communicator */
  MPI_Comm_rank(comm2d, rank);

  /** Get this rank's coordinates on the grid  */
  MPI_Cart_coords(comm2d, *rank, ndims, coords);

  /*
   * Each process/rank needs to know the ranks of its four neighbours in the grid.
   * MPI_PROC_NULL is assigned automatically.
   */
    MPI_Cart_shift( comm2d, 0,  1, left,   right );
    MPI_Cart_shift( comm2d, 1,  1, bottom, top );
    

}

/**
 * @brief The function below is a 1D decomposition that does not
 * requires problem size L to be the exact multiple of the number 
 * of processes for the program to run correctly.  This has no limitation on matching
 * problem size and number of processes.

 * Solution from:
 * https://ftp.mcs.anl.gov/pub/mpi/usingmpi-1st/examples/intermediate/decomp.f
 **/
void decomp1d( int n, int size, int position, int *mystart, int *mystop )
{

    int nlocal;
    int deficit;  // remainder after dividing simulation size (n) by number of processes (size)

    nlocal   = ceil(n/size); // get the largest interger value from n/size
    deficit  = (n%size); // keep the remainder

    *mystart	 = position * nlocal + 1;  // my start initial value
    if(position < deficit) {*mystart += position;} else { *mystart += deficit;} // add to my start, my rank or the remainder
    
    if (position < deficit)  nlocal += 1;  // 
    *mystop   = *mystart + nlocal - 1;  //
    if (*mystop > n || position == size-1)  *mystop = n;  //
    
}

/**
 * @brief 2D Decomposition; uses the 1D decomposition to compute the parameters below 
 *        
 * @param istart 
 * @param istop 
 * @param jstart 
 * @param jstop 
 */
void fnd2ddecomp(int l, int *istart, int *istop, int *jstart, int *jstop ) 
{
  /* Compute the decomposition */
   
    decomp1d( l, dims[0], coords[0], istart, istop );
    decomp1d( l, dims[1], coords[1], jstart, jstop );

}

/**
 * @brief broadcast an integer value from the controller processor
 * 
 * @param ival 
 */
void bcast(int *ival) 
{
    MPI_Bcast(ival, 1, MPI_INT, ROOT, comm);
}

void gathermlist(int *localist, int **masterlist, int nl) 
{
  MPI_Gather(localist, nl, MPI_INT, &(masterlist[0][0]), nl, MPI_INT, ROOT, comm);
  //printf("Gather Success!\n");
}

/**
 * @brief The controller process splits the map into small map (old)
 *        1. it copies its section of the map into old (with halos)
 *        2. then sends them to the desgnated process
 *        The rest of the processes receives their allocated map into old (with halos)
 *          Uses blocking point-to-point comm
 * 
 */
void split2dmap(int **map, int **old, int **alldimstore, int istart, int istop, int jstart, int jstop, int ml, int nl, int rank, int size)
{
  
  // Create new datatypes maptype and oldtype to send and receive M columns, respectively.
  MPI_Datatype oldtype;  //stride type for old
  MPI_Datatype maptype;  //stride type for map


  MPI_Type_vector(ml, nl, nl+2, MPI_INT, &oldtype);
  MPI_Type_commit(&oldtype);

  MPI_Status status;

  int tag = 1;
  int i;

  if(rank == ROOT)
  {
    // copy my section of the map
    initold(map, old, istart, istop, jstart, jstop);

    // send map sections to their corresponding processes on the grid
    for (i = 0; i < size; i++) {
      int dest, istartd, jstartd, mld, nld;
      nextrank(alldimstore, i, &dest, &istartd, &jstartd, &mld, &nld);


      MPI_Type_vector(mld, nld, nld*dims[1], MPI_INT, &maptype);
      MPI_Type_commit(&maptype);

      if(alldimstore[i][0] != ROOT){ // don't send it to myself
          MPI_Ssend(&map[istartd-1][jstartd-1], 1, maptype, dest, tag, comm);
      }
    }
  }
  else
  {
    // receive my section of the large map
    MPI_Recv(&old[1][1], 1, oldtype, ROOT, tag, comm, &status);
  }
}


/**
 * @brief Initialise Non-Blocking Boundary swaping
 *  the map is mostly periodic across the first dimension "i" 
 *  except certain cells remain non-periodic 
 *  (i.e. the halo values are set to zero for some values of "j" 
 *  as in the non-periodic case).
 */
void startbcswap(int **old, int ml, int nl, int left, int right, int top, int bottom, int size) 
{

    int tag = 1;
  // If single process, MPI_PROC_NULL will ensure no halo swap :
  
    //send right boundaries and receive left ones
    MPI_Isend(&old[ml][1], nl, MPI_INT, right, tag, comm, &request[0]);
    MPI_Recv(&old[0][1], nl, MPI_INT, left, tag, comm, &status[0]);

    //send left boundary and receive right
    MPI_Isend(&old[1][1], nl, MPI_INT, left, tag, comm, &request[1]);
    MPI_Recv(&old[ml+1][1], nl, MPI_INT, right, tag, comm, &status[1]);

    // This uses the vector datatype stridetype
    MPI_Datatype stridetype;
    MPI_Type_vector(ml, 1, nl+2, MPI_INT, &stridetype);
    MPI_Type_commit(&stridetype);

    //send top boundaries and receive bottom ones
    MPI_Isend(&old[1][nl], 1, stridetype, top, tag, comm, &request[2]);
    MPI_Recv(&old[1][0], 1, stridetype, bottom, tag, comm, &status[2]);

    // //send bottom boundary and receive top
    MPI_Isend(&old[1][1], 1, stridetype, bottom, tag, comm, &request[3]);
    MPI_Recv(&old[1][nl+1], 1, stridetype, top, tag, comm, &status[3]);


}

/**
 * @brief Finalise Non-blocking Boundary Swap
 * 
 */
void endbcswap(void)
{
  // wait for the send request to complete
  MPI_Waitall(4, request, status);
}


/**
 * @brief Sum up ival from all processes
 * 
 * @param ival 
 */
void globalsum(int *ival) 
{
  int tmpval;
  MPI_Allreduce(ival, &tmpval, 1, MPI_INT, MPI_SUM, comm);
  *ival = tmpval;
}


/**
 * @brief Merge small maps from processes (excluding halos) into 
 *        the desgnated space of the original map
 * 
 * Since this is a 2D decomposition a simple send/receive was used for simplicity
 * This uses Blocking point-to-point send/receive instead of MPI_Gather
 * 
 */
void merge2dmap(int **old, int **map, int **alldimstore, int istart, int istop, int jstart, int jstop, int ml, int nl, int rank, int size)
{

  // Create new mpi vector datatypes oldtype and maptype to send M columns.
  MPI_Datatype oldtype;
  MPI_Datatype maptype;

  MPI_Type_vector(ml, nl, nl+2, MPI_INT, &oldtype);
  MPI_Type_commit(&oldtype);

  MPI_Status status;

  int tag = rank;
  int i;

  if(rank != ROOT)
  {
    // send old (without halos) to the controller process
    MPI_Ssend(&old[1][1], 1, oldtype, ROOT, tag, comm);
  }
  else if(rank == ROOT)
  {
    // copy old into map on the controller process
    updatemap(map, old, istart, istop, jstart, jstop);

    //receive small map from proccesses into their designated space
    for (i = 0; i < size; i++) { 
      int source, istartd, jstartd, mld, nld;
      nextrank(alldimstore, i, &source, &istartd, &jstartd, &mld, &nld);

      tag = source;

      MPI_Type_vector(mld, nld, nld*dims[1], MPI_INT, &maptype);
      MPI_Type_commit(&maptype);

      if(alldimstore[i][0] != ROOT){
          //receiving data from source
          MPI_Recv(&map[istartd-1][jstartd-1], 1, maptype, source, tag, comm, &status);
      }
    }
  }
}



