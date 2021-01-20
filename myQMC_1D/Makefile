CC = g++
CFLAGS = -Wall -O3 --std=c++11

myQMC1D.exe : myQMC1D.o random.o
	$(CC) random.o myQMC1D.o -o myQMC1D.exe
myQMC1D.o : myQMC1D.cpp
	$(CC) -c myQMC1D.cpp -o myQMC1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o myQMC1D.exe
