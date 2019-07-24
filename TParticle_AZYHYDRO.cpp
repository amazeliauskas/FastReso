/*
 * Copyright (c) 2019 Aleksas Mazeliauskas, Stefan Floerchinger, 
 *                    Eduardo Grossi, and Derek Teaney
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
#include "TParticle_AZYHYDRO.h"
using namespace std;
//! Class contructor, which takes particle properties and the size of momentum grid
TParticle_AZYHYDRO::TParticle_AZYHYDRO(int pdgcode, std::string name, double mass, double gamma, int nu,
    double qb, double qs, double qc, double qbot, double isospin, double charge, int ndecays)
{
 fPDGCode = pdgcode;
 fParticleName = name;
 fMass = mass;
 fGamma = gamma;
 fNu = nu;
 fQB =qb;
 fQS = qs;
 fQC = qc;
 fQBot = qbot;
 fIsospin = isospin;
 fCharge = charge;
 fNdecays = ndecays;
//! Check that mass is  positive 
    if (fMass< 0){
      cerr << "\033[1mTParticle_AZYHYDRO.cpp\033[0m : \033[1;31merror\033[0m : Negative mass particle: M="<< fMass  << endl;
      exit(EXIT_FAILURE);
    }
    //! Decide on particle type
    if ( fNu % 2==0 and fMass > 0.001) {
      fParticleType =  EParticleType::kFermion; 
    } else if ( fNu % 2==1 and fMass > 0.001) {
      fParticleType = EParticleType::kBoson;
    } else if ( fNu % 2==0 and fMass <= 0.001) {
      //photon
      fParticleType = EParticleType::kBoson;
    }
    else 
    {
      fParticleType =  EParticleType::kBoltzman; 
    }
    // Initialize the integrated yield
    fNtotal=0;
  //! Momentum grid initialization
  init_grid();
    // printf("# \033[1mNew particle\033[0m : \033[1;32m%s\033[0m ;\tm = %9f GeV;\tnu = %d ;\n", fParticleName.c_str(), fM, fNu);
  char buffer[500];
 sprintf(buffer,"%9s\t%17s\t%9s\t%9s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
      "mc-number",  "name", "mass", "width", "deg","Qb","Qs","Qc","Qbot","isospin","charge","decays");
  fDescriptionHeader =  buffer;
      sprintf(buffer,"%9d\t%17s\t%9.5f\t%9.5e\t%d\t%g\t%3.g\t%g\t%g\t%g\t%3.g\t%3d",
          fPDGCode, fParticleName.c_str(), fMass, fGamma, fNu, fQB, fQS, fQC, fQBot, fIsospin, fCharge, fNdecays);
  fDescription =  buffer;
  }
TParticle_AZYHYDRO::TParticle_AZYHYDRO(string tag){
  FILE * pFile;
  string fname = tag+"_Fj.out";
  char buffer[200];
  pFile = fopen (fname.c_str(),"r");
  fscanf(pFile,"# %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n");
if( fscanf(pFile, "# %d  %s %lf %le %d %lg %lg %lg %lg %lg %lg %d %lf\n",
          &fPDGCode, buffer, &fMass, &fGamma, &fNu, &fQB, &fQS, &fQC, &fQBot, &fIsospin, &fCharge, &fNdecays,&fNtotal)!=13){
   cerr  << "\033[1mTParticle_AZYHYDRO.cpp\033[0m : \033[1;31merror\033[0m :  reading error " << endl;
    exit(EXIT_FAILURE);
}
  fParticleName = buffer;
//! Check that mass is  positive 
    if (fMass< 0){
      cerr << "\033[1mTParticle_AZYHYDRO.cpp\033[0m : \033[1;31merror\033[0m : Negative mass particle: M="<< fMass  << endl;
      exit(EXIT_FAILURE);
    }
    //! Decide on particle type
    if ( fNu % 2==0 and fMass > 0.001) {
      fParticleType =  EParticleType::kFermion; 
    } else if ( fNu % 2==1 and fMass > 0.001) {
      fParticleType = EParticleType::kBoson;
    } else if ( fNu % 2==0 and fMass <= 0.001) {
      //photon
      fParticleType = EParticleType::kBoson;
    }
    else 
    {
      fParticleType =  EParticleType::kBoltzman; 
    }
    // Initialize the integrated yield
    fNtotal=0;
  //! Momentum grid initialization
  init_grid();
  // Initialize the integrated yield
  fscanf(pFile,"#%*s [GeV]\t%*s [GeV]\t%*s 1\t%*s 2\t%*s 1\t%*s 2\t%*s 3\t%*s 1\t%*s 2\t%*s 1\t%*s 2\t%*s 1\t%*s 2\t%*s 3\n");
  for(int i = 0; i <grid_params::fNpbar; i++){
    if(fscanf(pFile,"%le\t%le", &fPbar_arr[i], &fMass)!=2){
    cerr << "\033[1mTParticle_AZYHYDRO.cpp\033[0m : \033[1;31merror\033[0m : fscan error"  << endl;
    exit(EXIT_FAILURE);
      }
    for(int j = 0; j <grid_params::fNf; j++){
     if(fscanf(pFile,"\t%le", &fFj_arr[j][i])!=1){
       
    cerr << "\033[1mTParticle_AZYHYDRO.cpp\033[0m : \033[1;31merror\033[0m : fscan error"  << endl;
    exit(EXIT_FAILURE);
       
     }
    }
    fscanf(pFile,"\n");
  }
  fclose(pFile);
}
