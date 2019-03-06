CC=/usr/local/cuda-5.0/bin/nvcc

all: q2.out

q2.out: q2.o
	$(CC) -o q2.out q2.o

q2.o: q2.cu
	$(CC) -o q2.o -c q2.cu

clean:
	rm -f q2.o q2.out
