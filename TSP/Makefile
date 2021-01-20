CC = g++
CFLAGS = -Wall -O3 --std=c++11

TSP : functions.o tsp.o random.o
	$(CC) functions.o tsp.o random.o -o TSP
functions.o : functions.cpp
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
tsp.o : tsp.cpp
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean : 
	rm *.o TSP
