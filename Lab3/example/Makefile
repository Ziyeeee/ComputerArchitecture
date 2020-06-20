# Makefile Performance Lab
CC = gcc
CFLAGS = -Wall -O2 -m32
LIBS = -lm

OBJS = vec.o benchmark.o combine.o cpe.o clock.o fcyc.o lsquare.o

all: benchmark

benchmark: $(OBJS) vec.h combine.h cpe.h fcyc.h lsquare.h clock.h
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o benchmark

clean: 
	-rm -f $(OBJS) benchmark *~ *.o


