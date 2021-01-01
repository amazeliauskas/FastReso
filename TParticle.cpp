/*
 * Copyright (c) 2018-2021 Aleksas Mazeliauskas, Stefan Floerchinger, 
 *                    Eduardo Grossi, and Derek Teaney
 * All rights reserved.
 *
 * FastReso is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/amazeliauskas/FastReso/
 */
#include <iostream>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <sstream>
#include "TParticle.h"
using namespace std;

//! Particle destructor
TParticle::~TParticle(){
  gsl_interp_accel_free(fPbar_acc);

  for (int jj=0; jj<grid_params::fNf; jj++){
  gsl_spline_free(fSpline_Fj[jj]);
  }
//   printf("# \033[1mDelete particle\033[0m : \033[1;32m%s\033[0m ;\tm = %9f GeV;\tnu = %d ;\n", fParticleName.c_str(), fMass, fNu);
}

void TParticle::init_grid(string comps){
  stringstream sstr(comps);
  string col;
  while (getline(sstr, col, ' ')){
    if (col=="Feq")         { 
      
fComponents.push_back(EFjIndex::kFeq1);
fComponents.push_back(EFjIndex::kFeq2);
      }
    else if (col=="Fshear") { 
    
fComponents.push_back(EFjIndex::kFshear1);
fComponents.push_back(EFjIndex::kFshear2);
fComponents.push_back(EFjIndex::kFshear3);
    
    }
    else if (col=="Fbulk")  { 
fComponents.push_back(EFjIndex::kFbulk1);
fComponents.push_back(EFjIndex::kFbulk2);
    }
    else if (col=="Ftemp")  {
fComponents.push_back(EFjIndex::kFtemp1);
fComponents.push_back(EFjIndex::kFtemp2);
    }
    else if (col=="Fvel")   {
fComponents.push_back(EFjIndex::kFvel1);
fComponents.push_back(EFjIndex::kFvel2);
fComponents.push_back(EFjIndex::kFvel3);
    }
    else { std::cerr << "\033[1mTFastReso.cpp\033[0m : \033[1;31mwarning\033[0m : unkwnown label " << col  <<std::endl; }
  }

  // Initialize Fj grid
  for (int jj=0; jj<grid_params::fNf; jj++){
  fSpline_Fj[jj] = gsl_spline_alloc(fSpline_type, grid_params::fNpbar);
  }
  fPbar_acc = gsl_interp_accel_alloc();
  //! Momentum grid initialization (use Tan sampling)
  //double dp = grid_params::fPbarMax/grid_params::fNpbar;
  for (int i = 0; i < grid_params::fNpbar; i++) {
    //fPbar_arr[i]=dp*(i+0.5); 
    fPbar_arr[i]=grid_params::fMref*tan(atan(grid_params::fPbarMax/grid_params::fMref)*(i+1)/(grid_params::fNpbar)); 
    for (int jj=0; jj<grid_params::fNf; jj++){
      fFj_arr[jj][i]= 0.0;
      fFj_arr_buffer[jj][i]= 0.0;
    }
  }

  for (int jj=0; jj<grid_params::fNf; jj++){
    fIsModified[jj]=true;
  }
}

//! Initialize grid interpolators for function Fj
void TParticle::init(int j){
  gsl_spline_init(fSpline_Fj[j], fPbar_arr, fFj_arr[j], grid_params::fNpbar);
  fIsModified[j]=false;
}
//! return value of function Fj at Ebar
double TParticle::get_Fj(int j, double Ebar) 
{
  //! if  modified, re-initialize the interpolator
  if ( fIsModified[j]  ){
    init(j);
  }
  double pbar = sqrt(Ebar*Ebar-fMass*fMass);
  // If pbar below lowest bin, return the first value of the grid
  if ( pbar >=fPbar_arr[0] and pbar <=fPbar_arr[grid_params::fNpbar-1] ){
    //  return interpolated value
    return gsl_spline_eval(fSpline_Fj[j], pbar, fPbar_acc);
  } else if ( pbar <fPbar_arr[0] ){
    //std::cerr << "\033[1mTParticle.h\033[0m : \033[1;31merror\033[0m : out of range pbar = " << pbar << std::endl;
    return fFj_arr[j][0];
    //  exit(EXIT_FAILURE);
  }
  // if pbar above the highest bin, return the logarithmically extrapolated value
  else if (pbar >fPbar_arr[grid_params::fNpbar-1]){

    int ip=grid_params::fNpbar-1;
    int im=grid_params::fNpbar-2;
    if (fFj_arr[j][ip]/fFj_arr[j][im]<=1 ) {
      //cout << pbar  <<" " << fPbar_arr[grid_params::fNpbar-1] << endl;
      double diffp= pbar-fPbar_arr[ip];
      //    exit(EXIT_FAILURE);
      return fFj_arr[j][ip]*exp(log(fFj_arr[j][ip]/fFj_arr[j][im])/(fPbar_arr[ip]-fPbar_arr[im])*diffp);
      //   return 0.0; //fFj_arr[j][ip]*pow(fFj_arr[j][ip]/fFj_arr[j][im], diffp/dp);
    }
    else if (fFj_arr[j][im]==0) { return 0.0; }
    else if (fFj_arr[j][ip]/fFj_arr[j][im]>1 and abs(fFj_arr[j][im]) < 10*DBL_EPSILON ) { return 0.0;}
    else if (fFj_arr[j][ip]/fFj_arr[j][im]>1 and abs(fFj_arr[j][im]) >= 10*DBL_EPSILON ) { 
      std::cerr << "\033[1mTParticle.h\033[0m : \033[1;31merror\033[0m : values increasing " << std::endl;
      exit(EXIT_FAILURE);
    }

  } 
  else {
    std::cerr << "\033[1mTParticle.h\033[0m : \033[1;31merror\033[0m : case not found " << std::endl;

    exit(EXIT_FAILURE);
  }
    exit(EXIT_FAILURE);
    return 0;
} 


////////////////////////////////////////////////////////////////////////////////
// Printing proceedures

// print Fj to the terminal
void TParticle::print(){
  printf("# %s\t%9s\n", fDescriptionHeader.c_str(), "Nyield");
  printf("# %s\t%9g\n", fDescription.c_str(), fNtotal);
  printf("# %s\n", fHeader.c_str());
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


  printf("# %s\t%9s\n", fDescriptionHeader.c_str(), "Nyield");
  printf("# %s\t%9g\n", fDescription.c_str(), fNtotal_buffer);
  printf("# %s\n", fHeader.c_str());
  for(int i = 0; i <grid_params::fNpbar; i++){
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
  fprintf(pFile,"# %s\t%9s\n", fDescriptionHeader.c_str(), "Nyield");
  fprintf(pFile,"# %s\t%9g\n", fDescription.c_str(), fNtotal);
  fprintf(pFile,"# %s", fHeader.c_str());
  int inx = 3;
  for(auto j: fComponents) {
    fprintf(pFile,"\t%d:%s",inx, fComponentNames[(int) j].c_str());
    inx++;
  }
  fprintf(pFile,"\n");
  for(int i = 0; i <grid_params::fNpbar; i++){
    fprintf(pFile,"%15.5e\t%15.5e", fPbar_arr[i], fMass);
    //for(int j = 0; j <grid_params::fNf; j++){
    for(auto j: fComponents) {
      fprintf(pFile,"\t%15.5e", fFj_arr[ (int)j][i]);
    }
    fprintf(pFile,"\n");
  }
  fclose(pFile);
}

void TParticle::print_buffer(string tag){
  FILE * pFile;
  string fname = tag+"_Fj.out";
  pFile = fopen (fname.c_str(),"w");
  fprintf(pFile,"# %s\t%9s\n", fDescriptionHeader.c_str(), "Nyield");
  fprintf(pFile,"# %s\t%9g\n", fDescription.c_str(), fNtotal_buffer);
  fprintf(pFile,"# %s\n", fHeader.c_str());
  for(int i = 0; i <grid_params::fNpbar; i++){
    fprintf(pFile,"%15.5e\t%15.5e", fPbar_arr[i], fMass);
    //for(int j = 0; j <grid_params::fNf; j++){
    for(auto j: fComponents) {
      fprintf(pFile,"\t%15.5e", fFj_arr_buffer[ (int) j][i]);
    }
    fprintf(pFile,"\n");
  }
  fclose(pFile);
}




