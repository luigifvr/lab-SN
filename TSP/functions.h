#include <vector>
#include "random.h"

#define PI 3.14159265

extern int seed[4];
extern Random rnd;

extern std::vector<std::vector<int>> population, newpopulation;
extern std::vector<int> start, gen1, gen2, son1, son2;
extern std::vector<double> positions_c;
extern std::vector<std::vector<double>> positions_s;

extern int npop, ncities, ngener, state;
extern double pMutations, pCrossover;

/////////////////////////////

extern int algorithm, nsteps, niter;
extern double beta_max, beta_step, beta;
extern double accepted, attempted;
extern std::vector<int> track;

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

///////////////
//SA

void InitialTrack(void);
void DoStep(void);
void MetroVerb(void);

///////////////

void Input();
void GeneratePositions(double);
void checkTrack(std::vector<int>);
double LossFunction(std::vector<int>);
double Diff(int, int);

void PrintPositions(void);
void PrintTrack(void);
void PrintBestLength(void);
void PrintAverageLength(void);
