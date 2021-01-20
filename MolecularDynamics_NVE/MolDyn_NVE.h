/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//restart option
int res;

//parameters, observables
const int m_props=200;
int n_props;
int iv,ik,it,ie, ip, igofr;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
double walker[m_props];

// averages
double acc,att;
double blk_av[m_props], glob_av[m_props], glob_av2[m_props], blk_norm;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
int nbins=100, print_gofr=1;
double bin_size;

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nblock=100, eq_steps;
double delta;
const double pi=3.1415927;

//functions
void Input(std::string);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int, int);
void Measure(void);
double Force(int, int);
double Pbc(double);

void Scale(void);
void MeanValues(int);
void Reset(int);
double Error(double, double, int);
void Accumulate(void);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
