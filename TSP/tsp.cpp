#include <iostream>
#include <cmath>
#include <assert.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "functions.h"

using namespace std;

int main(){
    srand(10);
    Input();
    GeneratePositions(state);
    if (algorithm==0) {
        GeneratePopulation();           //Genero la popolazione per GA
        SortPop();                      // Ordino in base alla loss function
        for (int i=0; i<ngener; ++i) {
            if (i%10==0) {
                cout << "Generation number = " << i << endl;
            }
            for (int j=0; j<int(npop/2); ++j) {
                Selection();            // Selezione dei genitori
                Crossover();            // Crossover dei geni
                
                // Mutazioni
                Mutation1(son1);
                Mutation1(son2);
                Mutation2(son1);
                Mutation2(son2);
                Mutation3(son1);
                Mutation3(son2);
                Mutation4(son1);
                Mutation4(son2);

                newpopulation.push_back(son1);
                newpopulation.push_back(son2);
                
                son1.clear();
                son2.clear();
            }
            population = newpopulation;         // Sostituisco la vecchia popolazione con la nuova
            newpopulation.clear();
            
            SortPop();
            PrintTrack();
            PrintBestLength();
            PrintAverageLength();
        }
    } else {
        // Simulated Annealing
        InitialTrack();
        for (int i=0; i<niter; ++i) {
            for (int j=0; j<nsteps; ++j) {
                DoStep();                   //Metropolis step
            }
            if (i%10==0) {
                MetroVerb();
            }
            beta += beta_step;
            PrintTrack();
            PrintBestLength();
        }
        
    }

    PrintPositions();
}

void Input(){
    ifstream ReadInput;
    
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    ReadInput.open("input.dat");
    
    ReadInput >> ncities;
    ReadInput >> pMutations;
    ReadInput >> pCrossover;
    
    ReadInput >> ngener;
    ReadInput >> npop;
    
    ReadInput >> nsteps;
    ReadInput >> beta_max;
    ReadInput >> beta_step;
    
    ReadInput >> state;
    ReadInput >> algorithm;
    
    cout << "The Traveling Salesman Problem... " << endl;
    if (state==0) {
        cout <<"... on a circle" << endl;
    } else{
        cout << "... inside a square" << endl;
    }
    cout << "Number of iterations = " << ngener << endl;
    cout << "Number of cities = " << ncities << endl;
    cout << "Number of individuals for the GA = " << npop << endl;
    cout << "Number of MC steps for the SA = " << nsteps << endl << endl;
    cout << "The algorithm used is..." << endl;
    if (algorithm==0) {
        cout << "...a Genetic Algorithm." << endl << endl;
    } else{
        cout << "...Simulated Annealing." << endl << endl;
    }
    
    ReadInput.close();
    
    for (int i=1; i<=ncities; ++i) {        // crea un vettore di partenza con città ordinate da 1 a ncities
        start.push_back(i);
    }
    niter = int(beta_max/beta_step);        // definisco il numero di iterazione per il SA
    
    return;
}

void GeneratePositions(double){
    double angle, x, y;
    std::vector<double> v;
    
    if (state == 0) {
        for (int i=0; i<ncities; ++i) {
            angle = rnd.Rannyu(0, 2*PI);        //genero gli angoli sulla circonferenza
            positions_c.push_back(angle);
        }
    } else{
        for (int i=0; i<ncities; ++i) {         // Genero le coordinate x ed y in un quadrato tra 0 e 1
            x = rnd.Rannyu();
            y = rnd.Rannyu();
            
            v.push_back(x);
            v.push_back(y);
            positions_s.push_back(v);
            
            v.clear();
        }
    }
    
    return;
}

void checkTrack(std::vector<int> track){
    
    assert(track[0]==1);            // Accerto che inizio dalla città 1
    
    sort(track.begin(), track.end());
    
    assert(track==start);           // accerto che il percorso è una permutazione del vettore di start (passo per le città una ed una sola volta)
    
    return;
}

double LossFunction(std::vector<int> track){            // Calcolo L1 come loss function
    double L=0;
    std::vector<int> nextpos = track;
    std::vector<double> distance;
    
    nextpos.push_back(1);
    std::transform(track.begin(), track.end(), nextpos.begin()+1, back_inserter(distance), Diff);
    
    L = accumulate(distance.begin(), distance.end(), 0.);
    
    return L;
}

double Diff(int x, int y){
    double diff_x, diff_y;
    
    if (state==0) {
        return abs(2*sin( (positions_c[x-1]-positions_c[y-1])/2 ));
    } else{
        diff_x = positions_s[x-1][0] - positions_s[y-1][0];
        diff_y = positions_s[x-1][1] - positions_s[y-1][1];
        return sqrt(diff_x*diff_x + diff_y*diff_y);
    }
}

void PrintPositions(void){
    ofstream Pos;
    
    Pos.open("positions.dat", ios::app);
    
    if (state==0) {
        for (std::vector<double>::const_iterator i = positions_c.begin(); i != positions_c.end(); i++) {
            Pos << *i << endl;
        }
    } else{
        for (std::vector<std::vector<double>>::const_iterator i = positions_s.begin(); i != positions_s.end(); i++) {
            for (std::vector<double>::const_iterator j = i->begin(); j != i->end(); j++) {
                Pos << *j << "   ";
            }
            Pos << endl;
        }
    }
    
    Pos.close();
    return;
}
