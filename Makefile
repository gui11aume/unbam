CC= gcc
CFLAGS= -std=gnu99 -g -Wall -O0

prog: *.c
	$(CC) $(CFLAGS) -o $@ $^ -lz

clean:
	rm -f prog *.o
