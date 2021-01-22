#include <vector>
#include "random.h"

#define PI 3.14159265

extern int seed[4];
extern Random rnd;

extern std::vector<std::vector<int>> population, newpopulation, send_migr, rec_migr;
extern std::vector<int> start, gen1, gen2, son1, son2, conn;
extern std::vector<double> positions_c;
extern std::vector<std::vector<double>> positions_s;

extern int size, rankMPI;
extern int nmigr, best_indiv, continente;
extern int npop, ncities, ngener, state;
extern double pMutations, pCrossover;

/////////////////////////////
//GA

void GeneratePopulation(void);
void Selection(void);
void Crossover(void);
void Mutation1(std::vector<int>&);
void Mutation2(std::vector<int>&);
void Mutation3(std::vector<int>&);
void Mutation4(std::vector<int>&);
void SortPop(void);
int GetIndex(std::vector<int>, int);

///////////////

void Input();
void GeneratePositions(void);
void checkTrack(std::vector<int>);
double LossFunction(std::vector<int>);
double Diff(int, int);

void PrintPositions(void);
void PrintTrack(void);
void PrintBestLength(void);
void PrintAverageLength(void);
