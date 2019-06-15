CC = gcc
CFLAGS = -g

default: test

test: test.o linalg.o
	$(CC) $(CFLAGS) -o test test.o linalg.o

test.o: test.c
	$(CC) $(CFLAGS) -c test.c

linalg.o: linalg.c
	$(CC) $(CFLAGS) -c linalg.c


.PHONY: clean
clean:
	rm -f test *.o *.gch
