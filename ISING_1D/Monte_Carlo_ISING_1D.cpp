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
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  Input(); //Inizialization
    for (int t=0; t<int(sizeof(temp_list)/sizeof(temp_list[0])); ++t)
    {
        beta = 1/temp_list[t];
        cout << "Equilibration period..." << endl;
        for (int eq_step=0; eq_step<eq_steps; ++eq_step) {
            Move(metro);
        }
        for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
        {
            Reset(iblk);   //Reset block averages
            for(int istep=1; istep <= nstep; ++istep)
            {
                Move(metro);
                Measure();
                Accumulate(); //Update block averages
            }
            Averages(iblk);   //Print results for current block
        }
        ConfFinal(); //Write final configuration
    }

  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConfig;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

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
  ReadInput.open("input.dat");

  ReadInput >> temp_i >> temp_f;
    
  cout << "Initial temperature = " << temp_i << endl;
    cout << "Final temperature = " << temp_f << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
    
  ReadInput >> eq_steps;
  cout << "Number of equilibration steps = " << eq_steps << endl;
    
  ReadInput >> old;  // Start from old configuration
    
  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
    if (old == 1) {
        cout << "Reading configuration in file old.0" << endl;
        ReadConfig.open("old.0");
        for (int i=0; i<nspin; ++i) {
            ReadConfig >> s[i];
        }
    } else {
        for (int i=0; i<nspin; ++i)
        {
            if(rnd.Rannyu() >= 0.5) s[i] = 1;
            else s[i] = -1;
        }
    }
    for (int i=0; i<int(sizeof(temp_list)/sizeof(temp_list[0])); ++i) {
        temp_list[i] = temp_i + (temp_f-temp_i)/temp_meas*i;
    }
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm, r;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
        sm = s[o];
        energy_old = Boltzmann(sm, o);
        
        if (s[o] == +1) {              // Cambio lo stato di uno spin
            s[o] = -1;
        } else{
            s[o] = +1;
        }
        energy_new = Boltzmann(s[o], o);
        if (energy_new - energy_old < 0.) {     // probabilità di effettuare un flip dello spin
            p = 1;
        } else{
            p = exp(-beta*(energy_new-energy_old));
        }
        
        r = rnd.Rannyu();
        if (r<=p) {
            accepted += 1;
        } else{
            s[o] = sm;
        }
        attempted += 1;
    }
    else //Gibbs sampling
    {
        p = 1/(1+exp(-2*beta*J*(s[Pbc(o-1)]+s[Pbc(o+1)])-2*h*beta));
        r = rnd.Rannyu();
        if (r<p) {
            s[o] = +1;
        } else{
            s[o] = -1;
        }
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0, chi=0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
      u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
      chi += s[i];
      m += s[i];
  }
    walker[iu] = u;
    walker[ic] = u*u;
    walker[im] = m;
    walker[ix] = pow(chi,2);
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
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output."+std::to_string(1/beta)+"_"+std::to_string(h)+".ene.0",ios::app);
    Heat.open("output."+std::to_string(1/beta)+"_"+std::to_string(h)+".heat.0",ios::app);
    Chi.open("output."+std::to_string(1/beta)+"_"+std::to_string(h)+".chi.0",ios::app);
    Mag.open("output."+std::to_string(1/beta)+"_"+std::to_string(h)+".mag.0",ios::app);
    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    stima_c = pow(beta,2)*(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm, 2))/(double)nspin; //Heat capacity
    stima_chi = beta*blk_av[ix]/blk_norm/(double)nspin;
    stima_m = blk_av[im]/blk_norm/(double)nspin;
    
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    glob_av[ix]  += stima_chi;
    glob_av2[ix] += stima_chi*stima_chi;
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    err_chi=Error(glob_av[ix],glob_av2[ix],iblk);
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Chi << setw(wd) << iblk << setw(wd) << stima_chi << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_chi << endl;
    Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    
    Ene.close();
    Heat.close();
    Chi.close();
    Mag.close();
    
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config_"+std::to_string(1/beta)+"_"+std::to_string(h)+".final " << endl << endl;
    WriteConf.open("config_"+std::to_string(1/beta)+"_"+std::to_string(h)+".final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
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
