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
#include <math.h>
#include <fstream>
#include <sstream>
#include "TParticle_THERMINATOR.h"

using namespace std;

//! Class contructor, which takes particle properties and the size of momentum grid
TParticle_THERMINATOR::TParticle_THERMINATOR(std::string name, double mass, double gamma, double spin,
    double isospin, double i3, double nq, double ns, double naq, double nas, double nc, double nac,
    int pdgcode, string comps)
{
  fParticleName =name;
  fMass = mass;
  fGamma = gamma;
  fSpin = spin;
  fIsospin = isospin;
  fI3 = i3;
  fNq = nq;
  fNs = ns;
  fNaq = naq;
  fNas = nas;
  fNc = nc;
  fNac = nac;

  fPDGCode =pdgcode;
  //! Check that mass is and calculated number of spins/polarization (photons are treated as 1 MeV particles)
  if (fMass>0.001){ fNu = int(2*fSpin+1); } else if (fMass==0 or fMass <= 0.001) { fNu = int(2*fSpin+1)-1; } else {
    cerr << "\033[1mTParticle_THERMINATOR.cpp\033[0m : \033[1;31merror\033[0m : Negative mass particle: M="<< fMass  << endl;
    exit(EXIT_FAILURE);
  }
  //    cout << fParticleName << " nu = " << fNu << endl;
  //! Decide on particle type (Boltzmann if spin is negative)
  if (fSpin>=0) 
  { 
    if ( int(2*fSpin+1) % 2==1) 
    {
      fParticleType = EParticleType::kBoson;
    } 
    else 
    {
      fParticleType =  EParticleType::kFermion; 
    }
  } 
  else 
  {
    fParticleType =  EParticleType::kBoltzman; 
  }

  // Calculate Baryon, Strangeness and Charm charges
  fQB =round((fNq+fNs+fNc-fNaq-fNas-fNac)/3.0);
  fQS = round(fNas-fNs);
  fQC = round(fNc-fNac);
  // Initialize the integrated yield
  fNtotal=0;


  

  //! Momentum grid initialization
  init_grid(comps);



  // printf("# \033[1mNew particle\033[0m : \033[1;32m%s\033[0m ;\tm = %9f GeV;\tnu = %d ;\n", fParticleName.c_str(), fM, fNu);
  char buffer[500];
  sprintf(buffer,"%9s\t%9s\t%9s\t%5s\t%4s\t%5s\t%s\t%s\t%s\t%s\t%s\t%s\t%9s",
      "Name", "Mass", "Gamma","Spin","Isospin","I3","Nq","Ns","Naq","Nas","Nc","Nac","MC");
  fDescriptionHeader =  buffer;
  sprintf(buffer,"%s\t%9f\t%9.5e\t%5.1f\t%5.1f\t%5.1f\t%g\t%g\t%g\t%g\t%g\t%g\t%9d",
      fParticleName.c_str(), fMass, fGamma, fSpin, fIsospin, fI3, fNq, fNs, fNaq, fNas, fNc, fNac, fPDGCode);
  fDescription =  buffer;

  
}

TParticle_THERMINATOR::TParticle_THERMINATOR(string tag){

  FILE * pFile;
  string fname = tag+"_Fj.out";
  char buffer[200];
  pFile = fopen (fname.c_str(),"r");
  fscanf(pFile,"# %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n");
if( fscanf(pFile, "# %s %lf %le %lf %lf %lf %lg %lg %lg %lg %lg %lg %d %lf\n",
      buffer, &fMass, &fGamma, &fSpin, &fIsospin, &fI3, &fNq, &fNs, &fNaq, &fNas, &fNc, &fNac, &fPDGCode, &fNtotal)!=14){
   cerr  << "\033[1mTParticle_THERMINATOR.cpp\033[0m : \033[1;31merror\033[0m :  reading error " << endl;
    exit(EXIT_FAILURE);
}

  fParticleName = buffer;

  //! Check that mass is and calculated number of spins/polarization (photons are treated as 0.01 MeV particles)
  if (fMass>0.01){ fNu = int(2*fSpin+1); } else if (fMass==0 or fMass < 0.01) { fNu = int(2*fSpin+1)-1; } else {
    cerr << "\033[1mTParticle_THERMINATOR.cpp\033[0m : \033[1;31merror\033[0m : Negative mass particle: M="<< fMass  << endl;
    exit(EXIT_FAILURE);
  }
  //    cout << fParticleName << " nu = " << fNu << endl;
  //! Decide on particle type (Boltzmann if spin is negative)
  if (fSpin>=0) 
  { 
    if ( int(2*fSpin+1) % 2==1) 
    {
      fParticleType = EParticleType::kBoson;
    } 
    else 
    {
      fParticleType =  EParticleType::kFermion; 
    }
  } 
  else 
  {
    fParticleType =  EParticleType::kBoltzman; 
  }

  // Calculate Baryon, Strangeness and Charm charges
  fQB =round((fNq+fNs+fNc-fNaq-fNas-fNac)/3.0);
  fQS = round(fNas-fNs);
  fQC = round(fNc-fNac);

  //! Momentum grid initialization
  init_grid();
  // Initialize the integrated yield
  fscanf(pFile,"#%*s [GeV]\t%*s [GeV]\t%*s 1\t%*s 2\t%*s 1\t%*s 2\t%*s 3\t%*s 1\t%*s 2\t%*s 1\t%*s 2\t%*s 1\t%*s 2\t%*s 3\n");
  for(int i = 0; i <grid_params::fNpbar; i++){
    if(fscanf(pFile,"%le\t%le", &fPbar_arr[i], &fMass)!=2){

    cerr << "\033[1mTParticle_THERMINATOR.cpp\033[0m : \033[1;31merror\033[0m : fscan error"  << endl;
    exit(EXIT_FAILURE);
      }
    for(int j = 0; j <grid_params::fNf; j++){
     if(fscanf(pFile,"\t%le", &fFj_arr[j][i])!=1){
       
    cerr << "\033[1mTParticle_THERMINATOR.cpp\033[0m : \033[1;31merror\033[0m : fscan error"  << endl;
    exit(EXIT_FAILURE);
       
     }
    }
    fscanf(pFile,"\n");
  }
  fclose(pFile);



}


