#
# Franklin and Hopper
#
CC = CC
MPCC = CC
OPENMP = -mp
CFLAGS = -O3
LIBS =


TARGETS = serial pthreads openmp mpi serial_mb

all:	$(TARGETS)

serial: serial.o prune_sweep.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o prune_sweep.o
serial_mb: serial_mb.o common.o
	$(CC) -o $@ $(LIBS) serial_mb.o common.o prune_sweep.o
pthreads: pthreads.o prune_sweep.o common.o
	$(CC) -o $@ $(LIBS) -lpthread pthreads.o prune_sweep.o common.o
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common.o mpi_migrants.o ghost.o prune_sweep.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o mpi_migrants.o ghost.o prune_sweep.o

openmp.o: openmp.cpp common.h microblock.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
serial_mb.o: serial_mb.cpp microblock.h common.h
	$(CC) -c $(CFLAGS) serial_mb.cpp
prune_sweep.o: prune_sweep.cpp common.h
	$(CC) -c $(CFLAGS) prune_sweep.cpp
pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp
mpi.o: mpi.cpp common.h mpi_particles.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
mpi_migrants.o: mpi_migrants.cpp common.h mpi_particles.h
	$(MPCC) -c $(CFLAGS) mpi_migrants.cpp
ghost.o: ghost.cpp common.h mpi_particles.h
	$(MPCC) -c $(CFLAGS) ghost.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp
dumb.o: dumb.cpp common.h
	$(CC) -c $(CFLAGS) dumb.cpp

clean:
	rm -f *.o $(TARGETS)
