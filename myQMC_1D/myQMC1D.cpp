/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include "myQMC1D.h"

using namespace std;

int main(int argc, char* argv[])
{ 
  Input(argv[1], argv[2]); //Inizialization
    if (var_flag) {
        cout << "Finding best parameters..." << endl;
        for (int i=0; i<mu_trials; ++i) {
            for (int j=0; j<sigma_trials; ++j) {
                SetParam(i, j);
                for (int eqstep=0; eqstep< eq_steps; eqstep++) {    // Equilibration
                    Move();
                }
                for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
                {
                    Reset(iblk);   //Reset block averages
                    for(int istep=1; istep <= nstep; ++istep)
                    {
                        Move();
                        Measure();
                        Accumulate(); //Update block averages
                    }
                    Averages(iblk);   //Print results for current block
                }
                CheckParam();
            }
        }
        PrintParam();
    } else {
        for (int eqstep=0; eqstep< eq_steps; eqstep++) {    // Equilibration
            Move();
        }
        for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
        {
            Reset(iblk);   //Reset block averages
            for(int istep=1; istep <= nstep; ++istep)
            {
                Move();
                Measure();
                Accumulate(); //Update block averages
            }
            Averages(iblk);   //Print results for current block
        }
    }
}


void Input(string inputFile, string varF)
{
  ifstream ReadInput,ReadParam;

  cout << "1D Quantum Particle        " << endl;
  cout << "Variational Monte Carlo simulation             " << endl << endl;
  cout << "Potential v(x) = x^4 - 2.5*x^2" << endl;
  cout << "Trial Wave Function psi(x; mu, sigma) = exp(-0.5*( (x-mu)/sigma )^2 + exp(-0.5*( (x+mu)/sigma )^2" << endl << endl;
    
//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
    ReadInput.open(inputFile);

    ReadInput >> mu_min;
    ReadInput >> mu_max;
  
    cout << "Minimum mean used is = " << mu_min << endl;
    cout << "Maximum mean used is = " << mu_max << endl;
    
    ReadInput >> sigma_min;
    ReadInput >> sigma_max;
    
    cout << "Minimum sigma used is = " << sigma_min << endl;
    cout << "Maximum sigma used is = " << sigma_max << endl;

    ReadInput >> length;
    cout << "Length of the box is = " << length << endl;

    ReadInput >> nblk;

    ReadInput >> nstep;

    ReadInput >> eq_steps;
    
    ReadInput >> delta;
    cout << "The program perform Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl;
    cout << "Number of equilibration steps = " << eq_steps << endl << endl;
    ReadInput.close();

    var_flag = atoi(varF.c_str());
    mu_trials = int((mu_max - mu_min)/var_step);
    sigma_trials = int((sigma_max - sigma_min)/var_step);
    
//Prepare arrays for measurements
    ie = 0; //Energy
 
    n_props = 1; //Number of observables
    
//histogram
    ihist = 1;
    nbins = 200;
    bin_size = length/(double)nbins;
    
    x = 0;
    ene_old = 10;
    
    if (!var_flag) {
        cout << "Reading best parameters..." << endl << endl;
        ReadParam.open("output.minvar.0");
        
        ReadParam >> mu;
        ReadParam >> sigma;
        cout << " mu = " << mu << " sigma = " << sigma << endl << endl;
    }
    return;
}

double Potential(double x){
    return pow(x,4) - 2.5*x*x;
}

double TrialWaveF(double x, double mu, double sigma){
    return exp(-0.5*pow( (x-mu)/sigma , 2)) + exp(-0.5*pow( (x+mu)/sigma , 2));
}

double TrialWaveF_secondDer(double x, double mu, double sigma){
    double first = exp(-0.5*pow( (x-mu)/sigma , 2));
    double second = exp(-0.5*pow( (x+mu)/sigma , 2));
    
    double secDer = pow(sigma, -2)*( -(first + second) + pow(sigma, -2)*(pow( x-mu, 2)*first + pow( x+mu, 2)*second ));
    return secDer;
}

double Distribution(double x, double mu, double sigma){
    double wf = TrialWaveF(x, mu, sigma);
    
    return wf*wf;
}

double Energy(double x, double mu, double sigma){
    double energy = 0;
    energy = (-0.5*TrialWaveF_secondDer(x, mu, sigma) + Potential(x)*TrialWaveF(x, mu, sigma))/TrialWaveF(x, mu, sigma);
    
    return energy;
}

void SetParam(int mu_t, int sigma_t){
    mu = mu_min + mu_t*var_step;
    sigma = sigma_min + sigma_t*var_step;
    x = 0;
    
    return;
}

void CheckParam(){
    double ene = glob_av[ie]/(double)nblk;
    if (ene < ene_old) {
        final_mu = mu;
        final_sigma = sigma;
        
        cout << mu << " " << sigma << " " << ene << endl;
        ene_old = ene;
    }
    return;
}

void PrintParam(void){
    ofstream Min;
    const int wd=12;
    
    Min.open("output.minvar.0");
    
    cout << "Best Parameters founds are: mu = " << final_mu << " sigma = " << final_sigma << endl;
    Min << final_mu << setw(wd) << final_sigma << setw(wd) << length << endl;
    Min << "mu" << setw(wd) << "sigma" << setw(wd) << "length" << endl;
    
    return;
}

void PrintPoint(double x){
    ofstream Points;
    
    Points.open("output.points.0", ios::app);

    Points << x << endl;
    
    Points.close();
    
    return;
}

void Move(void)
{
    double p, xold, xnew;
    
    xold = x;
    
  //New
    xnew = Pbc( x + rnd.Rannyu(-delta, delta) );

  //Metropolis test
    p = Distribution(xnew, mu, sigma)/Distribution(xold, mu, sigma);
    if(p >= rnd.Rannyu())  
    {
    //Update
       x = xnew;
    
       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
    if (!var_flag) {
        PrintPoint(x);
    }
  }

void Measure(void)
{
    int bin;
    double e = 0.0;
    e = Energy(x, mu, sigma);
    bin = int((length/2 + x)/bin_size);
  walker[ie] = e;
    walker[bin+ihist] += 1;
    return;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
   ofstream Ene;
   const int wd=12;
    
    stima_ene = blk_av[ie]/blk_norm; //Potential energy
    glob_av[ie] += stima_ene;
    glob_av2[ie] += stima_ene*stima_ene;
    err_ene=Error(glob_av[ie],glob_av2[ie],iblk);
    
    if (!var_flag) {
        cout << "Block number " << iblk << endl;
        cout << "Acceptance rate " << accepted/attempted << endl << endl;
        cout << "----------------------------" << endl << endl;
    
        Ene.open("output.energy.0",ios::app);
  
//Energy
        Ene << setw(wd) << iblk <<  setw(wd) << stima_ene << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_ene << endl;
    }
    Ene.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - length * rint(r/length);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
