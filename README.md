> # MPP Coursework Assessment: Parallel Percolation Simulation
> #### **Author: Exam Number B144511**
--------------------------------------------------------------------------

>> This is a parallel program (**MPI**) to test for percolation of a cluster - 2D Decomposition.  
When run on multiple processes, the code initiates simulation variables, including the simulation size (L) and number of processes, on the controller process and broadcasts the requied value to all processes.  
Then, using MPI cartesian topology, each process is arranged into a 2D grid and the grid information is stored into a new communicator.  
This information is then used to identify each processor's 4 neighbours (left/right; up/down), determine processor's section of the global map, and determine the local map size. The simulation arrays are then initialised on each process with appropriate size.  The controller process then splits the global map, sends them to the appropriate processes, while the all the other processes receive their part of the map. On each process, the percolation loop starts by initialising a halo swap then the update process takes place, followed by gathering the number of changes that occured, completing the halo swap before checking if the there were no changes in the update (to exit the loop if necessary). This is repeated for a defined number of steps.  Once the loop completes, all processes send their part of the updated map (without halos) to the controller process, and the controller processer receives them into the desgnated space of the map.


> ## Source and Header Files: (all at the top level of the source code)
--------------------------------------------------------------------------

>> **Source file and Functionality**
>> - **percolate.c:** The solution to the percolation problem.  
>> - **percolate.h:** Header for the solution that links to all other libraries required for the program.
>> - **percolatelib.c:** User defined functions update steps, array initialisation, copying map to old with halos, copying old to map
>> - **percio.c:**  I/O code to write the final image
>> - **mpilib.c:**  Communications library for MPI.  Includes all MPI functions to be called from the main program. Such as: 2d decomposition, Split map, Merge Map, Boundary swaps, etc.
>> - **mpilib.h:**  Header that links percolate.c with the mpilib.c (or serialib.c when running in serial)
>> - **serialib.c:**  Serial dummy for the percolate file, to run the program in serial.
>> - **arrayalloc.c:**   Dynamic memory allocation and dynamic store allocators for C.
>> - **arrayalloc.h:**   Header for the arraylloc.c
>> - **unirand.c:**      Random value generation library


> ## Building the code on Cirrus. 
--------------------------------------------------------------------------

>> - On Cirrus, make sure mpt and Intel compilers are loaded before compiling the code
>>>>>> module load mpt/2.22
>>>>>> module load intel-compilers-19

>> - Compile using "make" and the compiler flags specified in the Make file.
>> For example: 
>>>>>> make


> ## Running the code on P processes, including any command-line arguments.
--------------------------------------------------------------------------

>> Run the executable "percolate" with at least one argument.  The first argument must be the random number seed 7777. To run the program on Cirrus login node on 16 processes, use the following command:
>>> mpirun -n 16 ./percolate 7777 

>> ### Note:
>>> The default L = 768, rho = 0:4040
>>> This is adjustable with 6 additional acceptable arguments in the following order:
>>>> mpirun -n 16 *<program> *<seed> [L] [rho] [nclusters] [checkchange] [testmode] 

> #### You MUST run the program with the following required parameters:
----------------------------------------------------
>> - <program>  ./percolate
>> - <seed>      7777 (can be changed)

> #### Optional Parameter (in specific order):
----------------------------------------------------
>> L       { >0 - <= 100000 }
>> rho         { 0 - 1 }
>> nclusters   { 1 - 9 }
>> checkchange { 1=true | 0=false }
>> testmode    { 1=true |0=false }

> ##### Description of the optional parameters
>> - [L] is the problem size (L). 
>> - [rho] is the density.
>> - [ncluster] is the number of clusters to be printed on the map.pmg.
>> - [checkchange] 
>>> 1 = exit the program when there are no more changes to the map.
>>> 0 = exit the program when it reaches the set maximum iteration value.
>> - [testmode] 
>>> 1 = change the maxiter relative to the problem size; set [checkchange] to 0 used to limit how long the program should to run to give a reasonable perfomance assessment.

>> **NOTE:** For the optional params, if the last option is desired (i.e. testmode), then all the preceding params must be specified to the desired spec. For example:
>>> mpirun -n 1 ./percolate 7777 768 0.4040 2 1

>> The program should run on any L and nproc combinations supplied as runtime arguments.

>> The program is submitted to run using MPI (mpilib) by default. To run a it in serial mode, open the Makefile, replace the mpilib.c with serialib.c.


>> The code produces an output file called “map.pgm” using a call to “mapwrite”.

