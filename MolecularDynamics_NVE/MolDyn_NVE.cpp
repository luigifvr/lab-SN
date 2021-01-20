/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include <iomanip>

using namespace std;

int main(int argc, char* argv[]){
    int nconf = 1;
    Input(argv[1]);             //Inizialization
    cout << "Starting equilibration..." << endl;
    for (int eqstep=0; eqstep<eq_steps; ++eqstep) {
        if (eqstep%10==0) {
            Scale();
        }
        Move();     // Equilibration steps
        Measure();
    }
    cout << "End of equilibration..." << endl << endl;
    for(int iblock=1; iblock <= nblock; ++iblock){
        Reset(nblock);
        if(iblock%iprint == 0) cout << "Number of block: " << iblock << endl;
        for (int istep=1; istep<=nstep/nblock; ++istep) {
            Move();            //Move particles with Verlet algorithm
            Measure();
            Accumulate();
            if(istep%10 == 0){
                //ConfXYZ(nconf, iblock);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf += 1;
            }
        }
        MeanValues(iblock);
    }
    ConfFinal();         //Write final configuration to restart
    
  return 0;
}


void Input(string InputFile){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf,ReadActualConf, ReadOldConf;
    
    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;
        
    seed = 1;    //Set seed for random numbers
    srand(seed); //Initialize random number generator
        
    ReadInput.open(InputFile); //Read input
        
    ReadInput >> res;
    cout << "Restart option from old configuration =  " << res << endl;
        
    ReadInput >> temp;
        
    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;
        
    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol,1.0/3.0);
    cout << "Edge of the simulation box = " << box << endl;
        
    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;
    ReadInput >> eq_steps;
    
    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl;
    cout << "Number of blocks = " << nblock << endl;
    cout << "Number of equilibration steps = " << eq_steps << endl << endl;
        
    ReadInput.close();
    
    //Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    ip = 4; //Pressure
    n_props = 5; //Number of observables
    
    //g(r)
    igofr = 5;
    bin_size = box/2/(double)nbins;
    n_props += nbins;
    
    if (res==1) {
        cout << "Reading actual conf from old.final and old conf from old.0" << endl << endl;
        
        ReadActualConf.open("old.final");
        ReadOldConf.open("old.0");
        
        for (int i=0; i<npart; ++i){
            ReadActualConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        for (int i=0; i<npart; ++i){
            ReadOldConf >> xold[i] >> yold[i] >> zold[i];
            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
        }
        ReadActualConf.close();
        ReadOldConf.close();
    } else{
        //Read initial configuration
        cout << "Read initial configuration from file config.0 " << endl << endl;
        ReadConf.open("config.0");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();
        
        //Prepare initial velocities
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i=0; i<npart; ++i){
            vx[i] = rand()/double(RAND_MAX) - 0.5;
            vy[i] = rand()/double(RAND_MAX) - 0.5;
            vz[i] = rand()/double(RAND_MAX) - 0.5;
            
            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
        double sumv2 = 0.0, fs;
        for (int i=0; i<npart; ++i){
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];
            
            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;
        
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;
            
            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
    }
    return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
    
    double v, t, vij, vir, p;
    double dx, dy, dz, dr;
    int bin;
    ofstream Epot, Ekin, Temp, Etot, Press;
    
  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Press.open("output_pres.dat",ios::app);

    //reset the hystogram of g(r)
    for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;
    
  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;
    
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
         vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
         vir = 48*((1/pow(dr,12) - 0.5*1/pow(dr,6)));
//Potential energy & pressure
         v += vij;
         p += vir;
     }
// histogram g(r)
        if (dr < box/2) {
            bin = int(dr/bin_size);
            walker[bin+igofr] += 2;
        }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_pres = p/(3*vol) + rho*stima_temp;
    
    walker[iv] = stima_pot;
    walker[ik] = stima_kin;
    walker[ie] = stima_etot;
    walker[it] = stima_temp;
    walker[ip] = stima_pres;
    
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_pres << endl;
    
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();
    
    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf, WriteActualConf, WriteOldConf;
    
    if (res == -1) {
        cout << "Print actual configuration to file old.final and old configuration to file old.0 " << endl << endl;
        WriteActualConf.open("old.final");
        WriteOldConf.open("old.0");
        
        for (int i=0; i<npart; ++i){
            WriteActualConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
            WriteOldConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
        }
        WriteActualConf.close();
        WriteOldConf.close();
    } else {
        cout << "Print final configuration to file config.final " << endl << endl;
        WriteConf.open("config.final");
        
        for (int i=0; i<npart; ++i){
            WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
        }
        WriteConf.close();
    }
  return;
}


void ConfXYZ(int nconf, int block){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + "_" + to_string(block) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void Scale(void){
    double xnew[m_part], ynew[m_part], znew[m_part], fx[m_part], fy[m_part], fz[m_part];
    double scale, kin=0.0;       //Scale factor & Kinetic energy
    double v1[m_part], v2[m_part], v3[m_part];
    
    for(int i=0; i<npart; ++i){ //Force acting on particle i
        fx[i] = Force(i,0);
        fy[i] = Force(i,1);
        fz[i] = Force(i,2);
    }
    
    for(int i=0; i<npart; ++i){ //Verlet integration scheme
        
        xnew[i] = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
        ynew[i] = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
        znew[i] = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
        
        v1[i] = Pbc(xnew[i] - xold[i])/(2.0 * delta);
        v2[i] = Pbc(ynew[i] - yold[i])/(2.0 * delta);
        v3[i] = Pbc(znew[i] - zold[i])/(2.0 * delta);
    }
    //Kinetic Energy
    for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_temp = (2.0 / 3.0) * kin/(double)npart; //Temperature
    
    scale = sqrt(temp/stima_temp);
    
    for (int i=0; i<npart; ++i){
        v1[i] *= scale;
        v2[i] *= scale;
        v3[i] *= scale;
        
        xold[i] = Pbc(x[i] - v1[i] * delta);
        yold[i] = Pbc(y[i] - v2[i] * delta);
        zold[i] = Pbc(z[i] - v3[i] * delta);
    }
    return;
}

void MeanValues(int nblocks){
    ofstream AveKin, AvePot, AveTot, AveTemp, AvePres, AveGofr;
    double err_k, err_p, err_t, err_temp, err_pres, err_gdir;
    double r, gdir;
    const int wd=12;
    
    AveKin.open("ave_ekin.out", ios::app);
    AvePot.open("ave_epot.out", ios::app);
    AveTot.open("ave_etot.out", ios::app);
    AveTemp.open("ave_temp.out", ios::app);
    AvePres.open("ave_pres.out", ios::app);
    AveGofr.open("ave_gofr.out", ios::app);
    
    for (int i=0; i<igofr; ++i) {
        glob_av[i] += blk_av[i]/blk_norm;
        glob_av2[i] += pow(blk_av[i]/blk_norm,2);
    }
    err_p = Error(glob_av[iv], glob_av2[iv], nblocks);
    err_k = Error(glob_av[ik], glob_av2[ik], nblocks);
    err_t = Error(glob_av[ie], glob_av2[ie], nblocks);
    err_temp = Error(glob_av[it], glob_av2[it], nblocks);
    err_pres = Error(glob_av[ip], glob_av2[ip], nblocks);
    
    if (print_gofr) {
        for (int i=0; i<nbins; ++i) {
            r = i*bin_size;
            
            AveGofr << r << setw(wd) << r << setw(wd);
        }
        AveGofr << endl;
        
        print_gofr = 0;
    }
    
    for (int i=0; i<nbins; ++i) {
        r = i*bin_size;
        
        gdir = blk_av[i+igofr]/blk_norm/rho/npart/(4*pi/3*(pow(r+bin_size,3)-pow(r,3)));
        glob_av[i+igofr] += gdir;
        glob_av2[i+igofr] += gdir*gdir;
        err_gdir = Error(glob_av[i+igofr], glob_av2[i+igofr], nblocks);
        
        AveGofr << setw(wd) << glob_av[i+igofr]/(double)nblocks << setw(wd) << err_gdir;
    }
    AveGofr << endl;
    
    AveKin << glob_av[ik]/(double)nblocks << setw(wd) << err_k << endl;
    AvePot << glob_av[iv]/(double)nblocks << setw(wd) << err_p << endl;
    AveTot << glob_av[ie]/(double)nblocks << setw(wd) << err_t << endl;
    AveTemp << glob_av[it]/(double)nblocks << setw(wd) << err_temp << endl;
    AvePres << glob_av[ip]/(double)nblocks << setw(wd) << err_pres << endl;
    
    AveKin.close();
    AvePot.close();
    AveTot.close();
    AveTemp.close();
    AvePres.close();
    
    return;
}

void Reset(int iblk){
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
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void Accumulate(void) //Update block averages
{
    
    for(int i=0; i<n_props; ++i)
    {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
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
