/*
 * Copyright (c) 2018 Aleksas Mazeliauskas, Stefan Floerchinger, 
 *                    Eduardo Grossi, and Derek Teaney
 * All rights reserved.
 *
 * FastReso is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/amazeliauskas/FastReso/
 */
#include <iostream>
#include "TFastReso_AZYHYDRO.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include "TFastReso_formulas.h"


using namespace std;



//! Initialization of thermal particle distribution functions
void TFastReso_AZYHYDRO::do_thermal(double Tfo, double MuB, double MuI3, double MuS, double MuC, double Cs2){

  using namespace qag_params;
  gsl_integration_workspace * fWorkspacethermal ;
  fWorkspacethermal = gsl_integration_workspace_alloc (fLimit);
  for (std::vector<unique_ptr<TParticle_AZYHYDRO>>::iterator it = fParticleData.begin() ; it != fParticleData.end(); ++it){
    // don't initialize photon
    if ((*it)->getName()=="Gamma") continue;
    double m = (*it)->getM();
    double nu = (*it)->getNu();
    double QB = (*it)->getQB();
    double QS = (*it)->getQS();
    double QC = (*it)->getQC();
    double QMu=QB*MuB+QS*MuS+QC*MuC;
    gsl_function functhermal;
    functhermal.function = get_F_p_dvbar;
    thermal_params params = {m, Tfo, QMu,(*it)->getType() };
    functhermal.params = &params;
    using namespace qag_params;
    double result, error;
    // calculated yield
    gsl_integration_qag (&functhermal, 0, 1, fEpsAbs,fEpsRel, fLimit,fKey ,fWorkspacethermal, &result, &error); 
    // updated yield
    (*it)->addN(result*nu);
    for (int i=0; i < (*it)->getNpbar(); i++){
      double Ebar = (*it)->getEbar(i);
      for (int j=0; j < grid_params::fNf ; j++){
        (*it)->addFj(j,i, get_initial_Fj(j, Ebar,Tfo,QMu, (*it)->getType(), m, Cs2));
      }
    }

  }

  gsl_integration_workspace_free (fWorkspacethermal);
}

//! Read in the particle data file and create TParticle_AZYHYDRO class object, which is stored in fParticleData vector
void TFastReso_AZYHYDRO::read_particles_data(string inputname, bool verbose){
  // Open the input file
  ifstream fInputFile(inputname);
  if (fInputFile.fail()) {
    cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::read_particles_data\033[0m : \033[1;31merror\033[0m :  input file does not exist "  << inputname << endl;
    exit(EXIT_FAILURE);
  }
if (verbose) {
  printf("#%9s\t%17s\t%9s\t%9s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
      "mc-number",  "name", "mass", "width", "deg","Qb","Qs","Qc","Qbot","isospin","charge","decays");
}

  int PDGCode;
  string Name;
  double Mass, Gamma; 
  int Nu;
  double  Qb, Qs,Qc, Qbot, Isospin, Charge;
  int Ndecays;

  string line;
  while (getline(fInputFile, line)){
    //! check that there is the correct number of elements
    stringstream sstrc(line); int count=0; while (sstrc >> Name){count++;} if (count!=12) continue;
    stringstream sstr(line);
    sstr >> Name; 
    // check that the first non-whitespace character is not #, i.e. commented line
    if(Name.front()!='#'){
      PDGCode = atoi(Name.c_str());
      sstr >> Name >> Mass >> Gamma >> Nu >> Qb >> Qs >> Qc >> Qbot >> Isospin >> Charge >> Ndecays ; 

if (verbose) {
      printf("%9d\t%17s\t%9.5f\t%9.5f\t%d\t%g\t%3.g\t%g\t%g\t%g\t%3.g\t%3d\n",
          PDGCode, Name.c_str(), Mass, Gamma, Nu, Qb, Qs, Qc, Qbot, Isospin, Charge, Ndecays);
}
      fIndexTable[PDGCode] = fParticleData.size();
      fParticleData.push_back(move(unique_ptr<TParticle_AZYHYDRO>(new TParticle_AZYHYDRO(PDGCode, Name, Mass,Gamma, Nu, Qb, Qs, Qc, Qbot, Isospin, Charge, Ndecays))));
      // add anti-baryions
      if (Qb>0){
        fIndexTable[-PDGCode] = fParticleData.size();
        fParticleData.push_back(move(unique_ptr<TParticle_AZYHYDRO>(new TParticle_AZYHYDRO(-PDGCode, Name+"bar", Mass,Gamma, Nu, -Qb, -Qs, -Qc, -Qbot, Isospin, -Charge, Ndecays))));

if (verbose) {
        printf("%9d\t%17s\t%9.5f\t%9.5f\t%d\t%g\t%3.g\t%g\t%g\t%g\t%3.g\t%3d\n",
            -PDGCode, (Name+"bar").c_str(), Mass, Gamma, Nu, -Qb,-Qs, -Qc, -Qbot, Isospin, -Charge, Ndecays);
}
      }
    }
  }
if (verbose) {
  clog << "Number of particles read " << fParticleData.size() << endl;
}
  fParticleNumber = fParticleData.size();
}

//! Read the decay table and perform 2-body and 3-body decays
void TFastReso_AZYHYDRO::do_decays(string inputname, bool verbose){

  string line;

  // Open the input file
  ifstream fInputFile(inputname);
  if (fInputFile.fail()) {
    cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m :  input file does not exist "  << inputname<< endl;
    exit(EXIT_FAILURE);
  }
  if (verbose){
  // Print the header
  printf("#%9s\t%9s\t%9s\t%9s\t%9s\n",
      "Father", "Child1", "Child2","Child3","BranchingRatio");
  }

  int Father, Child1, Child2, Child3, Child4, Child5;
  int ndecayproducts;
  string Name;
  double BranchingRatio;

  int count_decays=0;
  int count_2decays=0;
  int count_3decays=0;

  while (getline(fInputFile, line)){
    //! check that there is the correct number of elements per line (5 for 2-body decay, 6 for three body decay), continue otherwise.
    stringstream sstrc(line); int count=0; while (sstrc >> Name){count++;} if (not (count==8)) continue;
    stringstream sstr(line);
    sstr >> Name; 
    // check that the first non-whitespace character is not #, i.e. commented line
    if(Name.front()!='#'){

      Father = atoi(Name.c_str());
      sstr>> ndecayproducts >> BranchingRatio >>  Child1 >> Child2 >> Child3 >> Child4 >> Child5; 
      ndecayproducts = abs(ndecayproducts);
      if (ndecayproducts==1) { continue;}
      else if (ndecayproducts==2) { count_2decays++;
if (verbose) {printf("%9d\t%9d\t%9d\t%f\n", Father, Child1, Child2, BranchingRatio);}
      }
      else if (ndecayproducts==3) {count_3decays++;
      if (verbose) { printf("%9d\t%9d\t%9d\t%9d\t%f\n", Father, Child1, Child2, Child3, BranchingRatio);}
      }
      else if (ndecayproducts==4) {
 if (verbose) {       cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : skipping 4 particle decay! " << endl;
        printf("%9d\t%9d\t%9d\t%9d\t%9d\t%f\n", Father, Child1, Child2, Child3, Child4, BranchingRatio);}
        continue;
      }
      else { 
        cerr << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : unrecognined number of decays! " << ndecayproducts << endl;
        cout << line << endl;
        exit(EXIT_FAILURE);}
      count_decays++;

      if (ndecayproducts==2)
      {
        double Ma = getParticleByPDG(Father)->getM();
        double QBa = getParticleByPDG(Father)->getQB();
        double GammaA = getParticleByPDG(Father)->getGamma();
        double Mb = getParticleByPDG(Child1)->getM();
        double Mc = getParticleByPDG(Child2)->getM();
        // Skip decays violating mass inequality (only possible by including resonance widths)
        if (Ma < Mb+Mc) {
     if(verbose){     cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc << endl;}
          continue;
        }
        //// Skip decays whose parents have decay width less than 10KeV 
        //if (GammaA <  10e-6) {
        //  cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle width smaller than 10 KeV ! Gamma = " << GammaA << " GeV" << endl;
        //  continue;
        //}

        double QS = getParticleByPDG(Father)->getQS();
        double QS1 = getParticleByPDG(Child1)->getQS();
        double QS2 = getParticleByPDG(Child2)->getQS();
        if (QS != QS1+QS2){
          cout << getParticleByPDG(Father)->getName() << endl;
          exit(EXIT_FAILURE);
        }
        do_2bodydecay(getParticleByPDG(Father), getParticleByPDG(Child1), getParticleByPDG(Child2), BranchingRatio);
        if (Child1!=Child2){ do_2bodydecay(getParticleByPDG(Father),getParticleByPDG( Child2),getParticleByPDG( Child1), BranchingRatio); }

        // decay anti-baryions
        if (QBa>0){
          count_2decays++;
      count_decays++;
 if (verbose){         printf("%9d\t%9d\t%9d\t%f\n", -Father, -Child1, -Child2, BranchingRatio);}
          do_2bodydecay(getParticleByPDG(-Father), getParticleByPDG(-Child1), getParticleByPDG(-Child2), BranchingRatio);
          if (Child1!=Child2){ do_2bodydecay(getParticleByPDG(-Father), getParticleByPDG(-Child2), getParticleByPDG(-Child1), BranchingRatio); }
        }




        //if (Child1=="pi0139plu" or Child2=="pi0139plu") { getParticleByPDG("pi0139plu")->print("pion_"+to_string(fIndexTable[Father]) +"_"+to_string(fIndexTable[Child1])+"_"+to_string(fIndexTable[Child2])); }
      }
      else if (ndecayproducts==3)
      {
        double Ma = getParticleByPDG(Father)->getM();
        double QBa = getParticleByPDG(Father)->getQB();
        double Mb = getParticleByPDG(Child1)->getM();
        double Mc = getParticleByPDG(Child2)->getM();
        double Md = getParticleByPDG(Child3)->getM();
        if (Ma < Mb+Mc+Md) {
       if(verbose) {   cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc+Md << endl;}
          continue;
        }

        double QS = getParticleByPDG(Father)->getQS();
        double QS1 = getParticleByPDG(Child1)->getQS();
        double QS2 = getParticleByPDG(Child2)->getQS();
        double QS3 = getParticleByPDG(Child3)->getQS();
        if (QS != QS1+QS2+QS3){
          cout << getParticleByPDG(Father)->getName() << endl;
          exit(EXIT_FAILURE);
        }
        do_3bodydecay(getParticleByPDG(Father), getParticleByPDG(Child1), getParticleByPDG(Child2),getParticleByPDG(Child3), BranchingRatio);
        if (Child1!=Child2){ do_3bodydecay(getParticleByPDG(Father), getParticleByPDG(Child2), getParticleByPDG(Child1),getParticleByPDG(Child3), BranchingRatio); }
        if (Child1!=Child3 and Child2!=Child3){ do_3bodydecay(getParticleByPDG(Father), getParticleByPDG(Child3),getParticleByPDG( Child1),getParticleByPDG(Child2), BranchingRatio); }

        // if (Child1=="pi0139plu" or Child2=="pi0139plu" or Child3=="pi0139plu") {getParticleByPDG("pi0139plu")->print("pion_"+to_string(fIndexTable[Father]) +"_"+to_string(fIndexTable[Child1])+"_"+to_string(fIndexTable[Child2])+"_"+to_string(fIndexTable[Child3])); }
        // decay anti-baryions
        if (QBa>0){
          count_3decays++;
      count_decays++;
 if (verbose){         printf("%9d\t%9d\t%9d\t%9d\t%f\n", -Father, -Child1, -Child2,-Child3, BranchingRatio);}

          do_3bodydecay(getParticleByPDG(-Father), getParticleByPDG(-Child1), getParticleByPDG(-Child2),getParticleByPDG(-Child3), BranchingRatio);
          if (Child1!=Child2){ do_3bodydecay(getParticleByPDG(-Father),getParticleByPDG( -Child2),getParticleByPDG( -Child1),getParticleByPDG(-Child3), BranchingRatio); }
          if (Child1!=Child3 and Child2!=Child3){ do_3bodydecay(getParticleByPDG(-Father), getParticleByPDG(-Child3), getParticleByPDG(-Child1),getParticleByPDG(-Child2), BranchingRatio); }
        }
      }


    }
  }


if (verbose) {
  clog << "Tota number of decays read " << count_decays << endl;
  clog << "Number of 2-body decays read " << count_2decays << endl;
  clog << "Number of 3-body decays read " << count_3decays << endl;
}

}



