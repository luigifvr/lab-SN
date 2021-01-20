/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __NVT__
#define __NVT__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, ie, ihist, nbins;
double walker[m_props], bin_size;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_ene,err_ene, ene_old;

// variational parameters
double mu_min, mu_max, sigma_min, sigma_max, length, mu, sigma;
double final_mu, final_sigma;

// simulation
int nblk, eq_steps, nstep, mu_trials, sigma_trials, var_flag;
double delta, var_step=0.01;
double x;

//pigreco
const double pi=3.1415927;

//functions
void Input(std::string, std::string);
double Potential(double);
double TrialWaveF(double, double, double);
double TrialWaveF_secDer(double, double, double);
double Distribution(double, double, double);
double Energy(double, double, double);
void SetParam(int, int);
void CheckParam(void);
void PrintParam(void);
void PrintPoint(double);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
double Pbc(double);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
