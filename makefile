# basic makefile. Barnett 4/16/18

FC = gfortran
CC = gcc
FLAGS = -O3 -march=native

.SUFFIXES: .c .f90

.f90.o:
	$(FC) $(FLAGS) -c $<
.c.o:
	$(CC) $(FLAGS) -c $<

SRC = bessel_eval.o cdriver.o

default: cdriver

# note gfortran is linker
cdriver: $(SRC)
	$(FC) $(SRC) -o cdriver
	./cdriver

clean:
	rm -f *.o cdriver
