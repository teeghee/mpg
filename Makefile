CFLAGS = -Wall -g

LDFLAGS = -lgmp

all: mpg

mpg: mpz_extra.o

.PHONY: clean

clean:
	rm *.o
	rm mpg
