CC=/opt/mvapich2-1.8-r5423/gcc/bin/mpicc

all: q1.out

q1.out: q1.o
	$(CC) -mpe=mpilog -o q1.out q1.o -lm

q1.o: q1.c
	$(CC) -mpe=mpilog -o q1.o -c q1.c -lm

clean:
	rm -f q1.o q1.out
