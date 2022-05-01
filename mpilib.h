
void mpstart(int *size, int *rank);

void mpstop(void);

double gettime(void);


void fnd2dnbrs(int size, int *rank, int *left, int *right, int *top, int *bottom);

void decomp1d( int n, int size, int position, int *mystart, int *mystop );

void fnd2ddecomp(int l, int *istart, int *istop, int *jstart, int *jstop );

void bcast(int *ival);

void globalsum(int *ival);

void startbcswap(int **old, int ml, int nl, int left, int right, int top, int bottom, int size);

void endbcswap(void);


void gathermlist(int *inbuff, int **outbuff, int nl);

void split2dmap(int **inbuff, int **outbuff, int **alldimstore, int istart, int istop, int jstart, int jstop, int ml, int nl, int rank, int size);

void merge2dmap(int **inbuff, int **outbuff, int **alldimstore, int istart, int istop, int jstart, int jstop, int ml, int nl, int rank, int size);

