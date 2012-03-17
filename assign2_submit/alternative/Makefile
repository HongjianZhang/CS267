#
# Franklin and Hopper
#
CC = CC
MPCC = CC
OPENMP = -mp
CFLAGS = -O3
LIBS =


TARGETS = serial mpi

all:	$(TARGETS)

#Serial
serial: serial.o common.o prune_sweep.o
	$(CC) -o $@ $(LIBS) serial.o common.o prune_sweep.o

#MPI
mpi: mpi.o common.o mpi_migrants.o ghost.o prune_sweep.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o mpi_migrants.o ghost.o prune_sweep.o

prune_sweep.o: prune_sweep.cpp common.h
	$(CC) -c $(CFLAGS) prune_sweep.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp

mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
mpi_migrants.o: mpi_migrants.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi_migrants.cpp
ghost.o: ghost.cpp common.h
	$(MPCC) -c $(CFLAGS) ghost.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS)
