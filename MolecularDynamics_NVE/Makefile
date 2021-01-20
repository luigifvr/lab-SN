CC = g++
CFLAGS = -Wall  --std=c++11

main : MD.o
	$(CC) MD.o -o MolecularDynamics
MD.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MD.o $(CFLAGS)
clean :
	rm *.o main