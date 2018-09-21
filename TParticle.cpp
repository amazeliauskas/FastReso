/*
 * Copyright (c) 2018 Aleksas Mazeliauskas
 * All rights reserved.
 *
 * FastReso is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/amazeliauskas/FastReso/
 */
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include "TParticle.h"

using namespace std;

//! Class contructor, which takes particle properties and the size of momentum grid
TParticle::TParticle(std::string name, double mass, double gamma, double spin,
    double isospin, double i3, double nq, double ns, double naq, double nas, double nc, double nac,
    int pdgcode): fParticleName(name), fMass(mass), fGamma(gamma), fSpin(spin),
  fIsospin(isospin), fI3(i3), fNq(nq), fNs(ns), fNaq(naq), fNas(nas), fNc(nc), fNac(nac),
  fPDGCode(pdgcode) {
//! Check that mass is and calculated number of spins/polarization
    if (fMass>0){ fNu = int(2*fSpin+1); } else if (fMass==0) { fNu = int(2*fSpin+1)-1; } else {
      cerr << "\033[1mTParticle.cpp\033[0m : \033[1;31merror\033[0m : Negative mass particle: M="<< fMass  << endl;
      exit(EXIT_FAILURE);
    }
    cout << fParticleName << " nu = " << fNu << endl;
    //! Decide on particle type (Boltzmann if spin is negative)
    if (fSpin>=0) { 
      if ( int(2*fSpin+1) % 2==1) {
        fParticleType = EParticleType::kBoson;
      } else {
        fParticleType =  EParticleType::kFermion; 
      }
    } else {
      fParticleType =  EParticleType::kBoltzman; 
    }
// Calculate Baryon, Strangeness and Charm charges
    fQB =round((fNq+fNs+fNc-fNaq-fNas-fNac)/3.0);
    fQS = round(fNas-fNs);
    fQC = round(fNc-fNac);
    // Initialize the integrated yield
    fNtotal=0;
    // Initialize Fj grid
    fSpline_Fj = gsl_spline_alloc(fSpline_type, grid_params::fNpbar);
    fPbar_acc = gsl_interp_accel_alloc();

    //! Momentum grid initialization (use Tan sampling)
    //double dp = grid_params::fPbarMax/grid_params::fNpbar;
    for (int i = 0; i < grid_params::fNpbar; i++) {
      //fPbar_arr[i]=dp*(i+0.5); 
      fPbar_arr[i]=grid_params::fMref*tan(atan(grid_params::fPbarMax/grid_params::fMref)*(i+0.5)/grid_params::fNpbar); 
      for (int jj=0; jj<grid_params::fNf; jj++){
        fFj_arr[jj][i]= 0.0;
        fFj_arr_buffer[jj][i]= 0.0;
        fIsModified[jj]=true;
        fIsLoaded[jj]=false;
      }
    }

    for (int jj=0; jj<grid_params::fNf; jj++){

      fIsLoaded[jj]=false;
      fIsModified[jj]=true;
    }
    // printf("# \033[1mNew particle\033[0m : \033[1;32m%s\033[0m ;\tm = %9f GeV;\tnu = %d ;\n", fParticleName.c_str(), fM, fNu);
  }


//! Particle destructor
TParticle::~TParticle(){
  gsl_interp_accel_free(fPbar_acc);
  gsl_spline_free(fSpline_Fj);
  // printf("# \033[1mDelete particle\033[0m : \033[1;32m%s\033[0m ;\tm = %9f GeV;\tnu = %d ;\n", fParticleName.c_str(), fM, fNu);
}
//! Initialize grid interpolators for function Fj
void TParticle::init(int j){
  gsl_spline_init(fSpline_Fj, fPbar_arr, fFj_arr[j], grid_params::fNpbar);
  fIsModified[j]=false;
  for (int jj=0; jj<grid_params::fNf; jj++){
    fIsLoaded[jj]=false;
  }
  fIsLoaded[j]=true;
}
//! return value of function Fj at Ebar
double TParticle::get_Fj(int j, double Ebar) 
{
  //! if not loaded or is modified, re-initialize the interpolator
  if (!fIsLoaded[j] or fIsModified[j]  ){
    init(j);
  }
  double pbar = sqrt(Ebar*Ebar-fMass*fMass);
  // If pbar below lowest bin, return the first value of the grid
  if ( pbar <fPbar_arr[0] ){
    //std::cerr << "\033[1mTParticle.h\033[0m : \033[1;31merror\033[0m : out of range pbar = " << pbar << std::endl;
    return fFj_arr[j][0];
    //  exit(EXIT_FAILURE);
  }
  // if pbar above the highest bin, return the logarithmically extrapolated value
  else if (pbar >fPbar_arr[grid_params::fNpbar-1]){

    int ip=grid_params::fNpbar-1;
    int im=grid_params::fNpbar-2;
    double diffp= pbar-fPbar_arr[ip];
    return fFj_arr[j][ip]*exp(log(fFj_arr[j][ip]/fFj_arr[j][im])/(fPbar_arr[ip]-fPbar_arr[im])*diffp);
    //  return 0.0; //fFj_arr[j][ip]*pow(fFj_arr[j][ip]/fFj_arr[j][im], diffp/dp);
  } 
  // else return interpolated value
  else {
    return gsl_spline_eval(fSpline_Fj, pbar, fPbar_acc);
  }
} 


////////////////////////////////////////////////////////////////////////////////
// Printing proceedures

// print Fj to the terminal
void TParticle::print(){
  printf("#%9s\t%9s\t%9s\t%5s\t%4s\t%5s\t%s\t%s\t%s\t%s\t%s\t%s\t%9s\t%9s\n",
      "Name", "Mass", "Gamma","Spin","Isospin","I3","Nq","Ns","Naq","Nas","Nc","Nac","MC", "Nyield");
  printf("# %s\t%9f\t%9.5e\t%5.1f\t%5.1f\t%5.1f\t%g\t%g\t%g\t%g\t%g\t%g\t%9d\t%9f\n",
      fParticleName.c_str(), fMass, fGamma, fSpin, fIsospin, fI3, fNq, fNs, fNaq, fNas, fNc, fNac, fPDGCode, fNtotal);
  printf("#%15s\t%s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n",
      "1:pbar [GeV]","2:m [GeV]", 
      "3:feq 1", "4:feq 2",
      "5:fshear 1", "6:fshear 2", "7:fshear 3", 
      "8:fbulk 1", "9:fbulk 2", 
      "10:ftemperature 1", "11:ftemperature 2",
      "12:fvelocity 1", "13:fvelocity 2", "14:fvelocity 3");
  for(int i = 0; i <grid_params::fNpbar; i++){
    printf("%15.5e\t%15.5e", fPbar_arr[i], fMass);
    for(int j = 0; j <grid_params::fNf; j++){
      printf("\t%15.5e", fFj_arr[j][i]);
    }
    printf("\n");
  }
}
// print buffered values of Fj to the terminal (for particular decay chain)
void TParticle::print_buffer(){

  printf("#%9s\t%9s\t%9s\t%5s\t%4s\t%5s\t%s\t%s\t%s\t%s\t%s\t%s\t%9s\t%9s\n",
      "Name", "Mass", "Gamma","Spin","Isospin","I3","Nq","Ns","Naq","Nas","Nc","Nac","MC", "Nyield");
  printf("# %s\t%9f\t%9.5e\t%5.1f\t%5.1f\t%5.1f\t%g\t%g\t%g\t%g\t%g\t%g\t%9d\t%9f\n",
      fParticleName.c_str(), fMass, fGamma, fSpin, fIsospin, fI3, fNq, fNs, fNaq, fNas, fNc, fNac, fPDGCode, fNtotal_buffer);
  printf("#%15s\t%s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n",
      "1:pbar [GeV]","2:m [GeV]", 
      "3:feq 1", "4:feq 2",
      "5:fshear 1", "6:fshear 2", "7:fshear 3", 
      "8:fbulk 1", "9:fbulk 2", 
      "10:ftemperature 1", "11:ftemperature 2",
      "12:fvelocity 1", "13:fvelocity 2", "14:fvelocity 3");  for(int i = 0; i <grid_params::fNpbar; i++){
    printf("%15.5e\t%15.5e", fPbar_arr[i], fMass);
    for(int j = 0; j <grid_params::fNf; j++){
      printf("\t%15.5e", fFj_arr_buffer[j][i]);
    }
    printf("\n");
  }
}
void TParticle::print(string tag){
  FILE * pFile;
  string fname = tag+"_Fj.out";
  pFile = fopen (fname.c_str(),"w");
  fprintf(pFile,"#%9s\t%9s\t%9s\t%5s\t%4s\t%5s\t%s\t%s\t%s\t%s\t%s\t%s\t%9s\t%9s\n",
      "Name", "Mass", "Gamma","Spin","Isospin","I3","Nq","Ns","Naq","Nas","Nc","Nac","MC", "Nyield");
  fprintf(pFile,"# %s\t%9f\t%9.5e\t%5.1f\t%5.1f\t%5.1f\t%g\t%g\t%g\t%g\t%g\t%g\t%9d\t%9f\n",
      fParticleName.c_str(), fMass, fGamma, fSpin, fIsospin, fI3, fNq, fNs, fNaq, fNas, fNc, fNac, fPDGCode, fNtotal);
  fprintf(pFile,"#%15s\t%s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n",
      "1:pbar [GeV]","2:m [GeV]", 
      "3:feq 1", "4:feq 2",
      "5:fshear 1", "6:fshear 2", "7:fshear 3", 
      "8:fbulk 1", "9:fbulk 2", 
      "10:ftemperature 1", "11:ftemperature 2",
      "12:fvelocity 1", "13:fvelocity 2", "14:fvelocity 3"); 
  for(int i = 0; i <grid_params::fNpbar; i++){
    fprintf(pFile,"%15.5e\t%15.5e", fPbar_arr[i], fMass);
    for(int j = 0; j <grid_params::fNf; j++){
      fprintf(pFile,"\t%15.5e", fFj_arr[j][i]);
    }
    fprintf(pFile,"\n");
  }
}

void TParticle::print_buffer(string tag){
  FILE * pFile;
  string fname = tag+"_Fj.out";
  pFile = fopen (fname.c_str(),"w");
  fprintf(pFile,"#%9s\t%9s\t%9s\t%5s\t%4s\t%5s\t%s\t%s\t%s\t%s\t%s\t%s\t%9s\t%9s\n",
      "Name", "Mass", "Gamma","Spin","Isospin","I3","Nq","Ns","Naq","Nas","Nc","Nac","MC", "Nyield");
  fprintf(pFile,"# %s\t%9f\t%9.5e\t%5.1f\t%5.1f\t%5.1f\t%g\t%g\t%g\t%g\t%g\t%g\t%9d\t%9f\n",
      fParticleName.c_str(), fMass, fGamma, fSpin, fIsospin, fI3, fNq, fNs, fNaq, fNas, fNc, fNac, fPDGCode, fNtotal_buffer);
  fprintf(pFile,"#%15s\t%s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n",
      "1:pbar [GeV]","2:m [GeV]", 
      "3:feq 1", "4:feq 2",
      "5:fshear 1", "6:fshear 2", "7:fshear 3", 
      "8:fbulk 1", "9:fbulk 2", 
      "10:ftemperature 1", "11:ftemperature 2",
      "12:fvelocity 1", "13:fvelocity 2", "14:fvelocity 3");
  for(int i = 0; i <grid_params::fNpbar; i++){
    fprintf(pFile,"%15.5e\t%15.5e", fPbar_arr[i], fMass);
    for(int j = 0; j <grid_params::fNf; j++){
      fprintf(pFile,"\t%15.5e", fFj_arr_buffer[j][i]);
    }
    fprintf(pFile,"\n");
  }
}




