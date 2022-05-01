/*
 *  Main header file for percolation exercise.
 */

/*
 *  Default System/Array size L = 768
 */

#define L 768

/*
 * 1)
 * It is easiest to program the exercise if you make ROOT, TRUE, FALSE compile time
 * constants (a #define
 */
#define ROOT 0

#define TRUE 1
#define FALSE 0
#define ndims 2


/*
 *  Prototypes for supplied functions
 */

void initold(int **map, int **old, int istart, int istop, int jstart, int jstop);

void updatestep(int **old, int **new, int ml, int nl, int *nchangelocal);

void updatemap(int **map, int **old, int istart, int istop, int jstart, int jstop);

void inithalos(int **old, int ml, int nl);

void nextrank(int **dstore, int i, int *rank, int *istart, int *jstart, int *m, int *n);

/*
 *  Visualisation
 */

void mapwrite(char *percfile, int **map, int l, int ncluster);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);

