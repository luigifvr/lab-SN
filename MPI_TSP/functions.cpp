#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "functions.h"
#include "random.h"

using namespace std;

int seed[4];
Random rnd;

int nmigr=75;       //Scambi ogni nmigr generazioni
int best_indiv=50;   // primi miglior individui da scambiare
int continente;
int size, rankMPI;
int npop, ncities, ngener, state=0;
double pMutations, pCrossover;

std::vector<std::vector<int>> population, newpopulation, send_migr, rec_migr;
std::vector<int> start, gen1, gen2, son1, son2, conn;
std::vector<double> positions_c;
std::vector<std::vector<double>> positions_s;

/////////////////////////////

int algorithm, nsteps, niter;
double beta_max, beta_step, beta=0;
double accepted, attempted;
std::vector<int> trackSA;

/////////////////////////////

//Genetic algorithm
void GeneratePopulation(void){
    std::vector<int> track, shuffle;
    
    shuffle = start;
    
    for (int indiv=0; indiv<npop; ++indiv) {
        random_shuffle(shuffle.begin()+1, shuffle.end());
        
        track.push_back(1);
        track.insert(std::end(track), std::begin(shuffle)+1, std::end(shuffle));
        
        checkTrack(track);
        
        population.push_back(track);
        
        track.clear();
    }
    return;
}

void SortPop(void){
    
    sort(population.begin(), population.end(), [](const std::vector<int>& first, const std::vector<int>& second){
        return LossFunction(first) < LossFunction(second);}
         );
    return;
}

void Selection(void){
    double sumL=0, partL=0, r, s;
    std::vector<std::vector<int>>::const_iterator j = population.begin();
    
    for (std::vector<std::vector<int>>::const_iterator i = population.begin(); i != population.end(); ++i) {
        sumL += LossFunction(*i);
    }
    r = rnd.Rannyu(0., sumL);
    s = rnd.Rannyu(0., sumL);
    
    do {
        partL += LossFunction(*j);
        ++j;
    } while (partL < r);
    gen1 = population[j - population.begin()-1];

    j = population.begin();
    do {
        partL += LossFunction(*j);
        ++j;
    } while (partL < s);
    gen2 = population[j - population.begin()-1];
    
    checkTrack(gen1);
    checkTrack(gen2);
    return;
}

void Crossover(void){
    int n;
    double r;
    std::vector<int> gen1old = gen1, gen2old = gen2;
    
    r = rnd.Rannyu();
    
    if (r<pCrossover) {
        n = rand()%(ncities-2) + 2;
        
        for (int i=0; i<ncities-n; ++i) {
            son1.push_back(gen1[i]);
            son2.push_back(gen2[i]);
        }
        
        for (int i=0; i<ncities-n; ++i) {
            gen1.erase(std::remove(gen1.begin(), gen1.end(), son2[i]), gen1.end());
            gen2.erase(std::remove(gen2.begin(), gen2.end(), son1[i]), gen2.end());
        }
        son1.insert(son1.end(), gen2.begin(), gen2.end());
        son2.insert(son2.end(), gen1.begin(), gen1.end());
        
    } else{
        son1 = gen1old;
        son2 = gen2old;
    }
    
    checkTrack(son1);
    checkTrack(son2);
    return;
}

void Mutation1(std::vector<int> &itrk){
    int p1, p2;
    double r;
    
    r = rnd.Rannyu();
    if (r<pMutations) {
        do {
            p1 = rand()%(ncities-1) + 1;
            p2 = rand()%(ncities-1) + 1;

        } while (p1 == p2);
        
        std::swap(itrk[p1], itrk[p2]);
    }
    
    checkTrack(itrk);
    return;
}

void Mutation2(std::vector<int> &itrk){
    int n, m, pos;
    double r;
    
    r = rnd.Rannyu();
    
    if (r < pMutations) {
        m = rand()%(ncities-2) + 1;     // number of contiguos cities
        pos = rand()%(ncities - (m+1) ) + 1;         // position of the first city
        n = rand()%(ncities-(pos+m)) + 1;             //shift length
        
        std::rotate(itrk.begin()+pos, itrk.begin()+pos+m, itrk.begin()+pos+m+n);
    }
    
    checkTrack(itrk);
    return;
}

void Mutation3(std::vector<int> &itrk){
    int m, lpos, fpos;
    double r;
    
    r = rnd.Rannyu();
    
    if (r<pMutations) {
        m = rand()%((int) floor((double(ncities)-1)/2)) + 1;          // Number pf contiguos cities
        
        fpos = rand()%((int) floor((double(ncities)-1)/2 - m+1)) + 1;    // position of first city
        lpos = fpos + m + rand()%(ncities-(fpos+2*m-1));           // position of first city to swap
        
        std::swap_ranges(itrk.begin()+fpos, itrk.begin()+fpos+m, itrk.begin()+lpos);
    }
    
    checkTrack(itrk);
    return;
}

void Mutation4(std::vector<int> &itrk){
    double r;
    int m, pos;
    
    r = rnd.Rannyu();
    
    if (r<pMutations) {
        m = rand()%(ncities-2) + 2;
        pos = rand()%(ncities - m) + 1;
        
        std::reverse(itrk.begin()+pos, itrk.begin()+pos+m);
    }
    
    checkTrack(itrk);
    return;
}

int GetIndex(std::vector<int> vec, int val){
    std::vector<int>::const_iterator it;
    
    it = std::find(vec.begin(), vec.end(), val);
    return it-vec.begin();
}

void PrintTrack(void){
    ofstream Track;
    
    Track.open("best_tracks"+std::to_string(rankMPI)+".dat", ios::app);
    
    for (std::vector<int>::const_iterator i = population[0].begin(); i != population[0].end(); ++i) {
        Track << *i << "    ";
    }
    Track << endl;
    
    Track.close();
    return;
}

void PrintBestLength(void){
    ofstream Len;
    
    Len.open("best_lengths"+std::to_string(rankMPI)+".dat", ios::app);
    
    Len << LossFunction(population[0]) << endl;
    
    Len.close();
    return;
}

void PrintAverageLength(void){
    ofstream AvLen;
    double len=0;
    
    AvLen.open("ave_length"+std::to_string(rankMPI)+".dat", ios::app);
    
    for (int i = 0; i<npop/2; ++i) {
        len += LossFunction(population[i]);
    }
    
    AvLen << len/(npop/2) << endl;
    
    AvLen.close();
    return;
}
