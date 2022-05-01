MF=	Makefile

# For Local (MPICH)
CC=	mpicc
CFLAGS=	-O3 -lm -Wall

# For Cirrus
#CC=	mpicc
#CFLAGS=	-cc=icc

# For ARCHER2
#CC=	cc
#CFLAGS=	-O3 -Wall

LFLAGS= $(CFLAGS)

EXE=	percolate

INC= \
	percolate.h \
	arralloc.h \
	mpilib.h

SRC= \
	percolate.c \
	percolatelib.c \
	percio.c \
	unirand.c \
	arralloc.c \
	mpilib.c 

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(OBJ):	$(INC)

$(EXE):	$(OBJ)
	$(CC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ):	$(MF)

clean:
	rm -f $(EXE) $(OBJ) core
