CC=/usr/local/cuda-5.0/bin/nvcc

all: bout

bout: csrmv.o
	$(CC) -o bout csrmv.o

csrmv.o: csrmv.cu
	$(CC) -o csrmv.o -c csrmv.cu

clean:
	rm -f csrmv.o bout
