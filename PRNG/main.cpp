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
#include <string>
#include "random.h"
#include <math.h>

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   int lambda, mu, gamma;
   int x, y, z;
   double sx=0., sy=0., sz=0.;
   int N[4] = {1, 2, 10, 100};

   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
   for(int i=0; i<4; i++){
      ofstream unif("Numbers/Unif" + std::to_string(N[i]) + ".txt");
      ofstream expo("Numbers/Expo" + std::to_string(N[i]) + ".txt");
      ofstream lorentz("Numbers/Lorentz" + std::to_string(N[i]) + ".txt");

   }   
   for(int i=0; i<4; i++){

      for(int j=0; j<10000; j++){
	    
	 for(int k=0; k<N[i]; k++){
            x = round(rnd.Rannyu(0.5, 6.5));  
	    y = ceil(rnd.Expo(lambda=1));
	    if (y>6){
	       y=6; }
            do{
	    z = ceil(abs(rnd.Lorentz(mu=0, gamma=1)));
	    }
	    while (z>=7.);

            sz += z;	    
	    sx += x;
 	    sy += y;
  
         }
	 ofstream unif("Numbers/Unif" + std::to_string(N[i]) + ".txt", ios_base::app);
      	 ofstream expo("Numbers/Expo" + std::to_string(N[i]) + ".txt", ios_base::app);
     	 ofstream lorentz("Numbers/Lorentz" + std::to_string(N[i]) + ".txt", ios_base::app);
         unif << sx/N[i] << endl;   
	 expo << sy/N[i] << endl;
	 lorentz << sz/N[i] << endl;  
         sx = 0.;
	 sy = 0.;
	 sz = 0.;
      } 
   }

   rnd.SaveSeed();
   return 0;
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
