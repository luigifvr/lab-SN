#include <iostream>
#include <cmath>
#include <assert.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "functions.h"
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]){
    MPI_Status status;
    MPI_Request req;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI);

    srand(rankMPI);
    Input();
    GeneratePopulation();
    SortPop();
    
    for (int i=0; i<ngener; ++i) {
        if (i%10==0 && rankMPI==0) {
            cout << "Generation number = " << i << endl;
        }
        if (i%nmigr==0) {
            if (rankMPI==0) {
                random_shuffle(conn.begin()+1, conn.end()-1);
            }
            
            MPI_Bcast(&conn[0], size+1, MPI_INTEGER, 0, MPI_COMM_WORLD);    // faccio il broadcast delle connessioni fra continenti
            rec_migr.resize(best_indiv, std::vector<int>(ncities));         // vettore che riceve gli individui
            send_migr.insert(send_migr.end(), population.begin(), population.begin()+best_indiv);       //vettore che invia gli individui
            population.erase(population.begin(), population.begin()+best_indiv);
            
            // Scambio gli individui tra i continenti       conn[GetIndex(conn, rankMPI)+1]
            for (int indiv=0; indiv<best_indiv; ++indiv) {
                MPI_Isend(&send_migr[indiv][0], ncities, MPI_INTEGER, conn[GetIndex(conn, rankMPI)+1], 0, MPI_COMM_WORLD, &req);
                MPI_Recv(&rec_migr[indiv][0], ncities, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            }
            population.insert(population.end(), rec_migr.begin(), rec_migr.end());
            SortPop();
            send_migr.clear();
            rec_migr.clear();
        }
        for (int j=0; j<int(npop/2); ++j) {
            Selection();
            Crossover();

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
        population = newpopulation;
        newpopulation.clear();
        
        SortPop();
        PrintTrack();
        PrintBestLength();
        PrintAverageLength();
    }
    MPI_Finalize();
}

void Input(){
    ifstream ReadInput, ReadPos;
    
    int p1, p2;
    ifstream Primes("Primes");
    for (int prim=0; prim<=rankMPI; ++prim) {
        Primes >> p1 >> p2 ;
    }
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

    ReadInput.close();
    
    for (int i=1; i<=ncities; ++i) {        // crea un vettore di partenza con città ordinate da 1 a ncities
        start.push_back(i);
    }
    for (int i=0; i<size; ++i) {        // crea il vettore che definirà le connessioni tra i continenti
        conn.push_back(i);
    }
    conn.push_back(0);
    
    ReadPos.open("positions.dat");
    positions_s.resize(ncities, std::vector<double>(2));
    for (std::vector<std::vector<double>>::iterator i = positions_s.begin(); i != positions_s.end(); i++) {
        for (std::vector<double>::iterator j = i->begin(); j != i->end(); j++) {
            ReadPos >> *j;
        }
    }
    return;
}

void GeneratePositions(void){
    double x, y;
    std::vector<double> v;
    
    for (int i=0; i<ncities; ++i) {
        x = rnd.Rannyu();
        y = rnd.Rannyu();
            
        v.push_back(x);
        v.push_back(y);
        positions_s.push_back(v);
            
        v.clear();
    }
    
    return;
}

void checkTrack(std::vector<int> track){
    
    assert(track[0]==1);
    
    sort(track.begin(), track.end());
    
    assert(track==start);
    
    return;
}

double LossFunction(std::vector<int> track){
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
    
    diff_x = positions_s[x-1][0] - positions_s[y-1][0];
    diff_y = positions_s[x-1][1] - positions_s[y-1][1];
    return sqrt(diff_x*diff_x + diff_y*diff_y);
}

void PrintPositions(void){
    ofstream Pos;
    
    Pos.open("positions.dat"+std::to_string(rankMPI), ios::app);

    for (std::vector<std::vector<double>>::const_iterator i = positions_s.begin(); i != positions_s.end(); i++) {
        for (std::vector<double>::const_iterator j = i->begin(); j != i->end(); j++) {
            Pos << *j << "   ";
        }
        Pos << endl;
    }
    
    Pos.close();
    return;
}
