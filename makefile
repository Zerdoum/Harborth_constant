# makefile pour somme nulle
#

all: sommenulle.o
	gcc -o sommenulle sommenulle.o -march=native

sommenulle.o: sommenulle.c
	gcc -c -Wall -Ofast -m64 sommenulle.c
clean:
	rm -f sommenulle.o sommenulle
exec: sommenulle
	./sommenulle
