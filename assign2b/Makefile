# Load CUDA using the following command
# module load cuda

#
# Dirac
#
CC = nvcc
MPCC = nvcc
OPENMP = 
CFLAGS = -O3 -arch sm_20
NVCCFLAGS = -O3 -arch sm_20
LIBS = -lm

TARGETS = serial gpu

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
gpu: gpu.o common.o gpu_main.o
	$(CC) -o $@ $(NVCCLIBS) gpu.o common.o gpu_main.o

serial.o: serial.cu common.h
	$(CC) -c $(CFLAGS) serial.cu
gpu.o: gpu.cu common.h
	$(CC) -c $(NVCCFLAGS) gpu.cu
gpu_main.o: gpu_main.cu common.h
	$(CC) -c $(NVCCFLAGS) gpu_main.cu
common.o: common.cu common.h
	$(CC) -c $(CFLAGS) common.cu

clean:
	rm -f *.o $(TARGETS)
