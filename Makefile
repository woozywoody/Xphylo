Xphylo : GTree.o Xphylo.o
	gcc -D BUFFERLEN=100000 -g -Wall -o Xphylo GTree.o Xphylo.o

GTree.o : GTree.c GTree.h
	gcc -Wall -g -c GTree.c

Xphylo.o : Xphylo.c
	gcc -Wall -g -c Xphylo.c

clean:
	rm *.o Xphylo
