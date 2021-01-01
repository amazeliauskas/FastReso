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
#include "TFastReso_THERMINATOR.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include "TFastReso_formulas.h"

using namespace std;
//! Initialization of thermal particle distribution functions
void TFastReso_THERMINATOR::do_thermal(double Tfo, double MuB, double MuI3, double MuS, double MuC, double Cs2){

  Cs2 = cs2FluiduM(Tfo);
  using namespace qag_params;
  gsl_integration_workspace * fWorkspacethermal ;
  fWorkspacethermal = gsl_integration_workspace_alloc (fLimit);
  for (std::vector<unique_ptr<TParticle_THERMINATOR>>::iterator it = fParticleData.begin() ; it != fParticleData.end(); ++it){
    // don't initialize photon
    if ((*it)->getM() < 0.001) continue;
    double m = (*it)->getM();
    double nu = (*it)->getNu();
    double QB = (*it)->getQB();
    double I3 = (*it)->getI3();
    double QS = (*it)->getQS();
    double QC = (*it)->getQC();
    double QMu=QB*MuB+I3*MuI3+QS*MuS+QC*MuC;
    gsl_function functhermal;
    functhermal.function = get_F_p_dvbar;
    thermal_params params = {m, Tfo, QMu,(*it)->getType() };
    functhermal.params = &params;
    using namespace qag_params;
    double result, error;
    // calculated yield
int status = gsl_integration_qag (&functhermal, 0, 1, fEpsAbs,fEpsRel, fLimit,fKey ,fWorkspacethermal, &result, &error); 
  // Check status of integration and reduce accuracy if failing
  if (status) { 
status = gsl_integration_qag (&functhermal, 0, 1, fEpsAbs,1e-4, fLimit,GSL_INTEG_GAUSS15 ,fWorkspacethermal, &result, &error); 
    if (status) { 
status = gsl_integration_qag (&functhermal, 0, 1, 1e-4,1e-4, fLimit,GSL_INTEG_GAUSS15 ,fWorkspacethermal, &result, &error); 
      if (status) {      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }


    // updated yield
    (*it)->addN(result*nu);
    for (int i=0; i < (*it)->getNpbar(); i++){
      double Ebar = (*it)->getEbar(i);
    //  for (int j=0; j < grid_params::fNf ; j++){
    for(auto j: (*it)->fComponents) {
        (*it)->addFj( (int)j,i, get_initial_Fj(j, Ebar,Tfo,QMu, (*it)->getType(), m, Cs2));
      }
    }

  }

  gsl_integration_workspace_free (fWorkspacethermal);
}
void TFastReso_THERMINATOR::do_thermal(int pi, double Tfo, double QMu, double Cs2){
Cs2 = cs2FluiduM(Tfo); 
  using namespace qag_params;
  gsl_integration_workspace * fWorkspacethermal ;
  fWorkspacethermal = gsl_integration_workspace_alloc (fLimit);
    // don't initialize photon
    double m = getParticle(pi)->getM();
    double nu = getParticle(pi)->getNu();
    gsl_function functhermal;
    functhermal.function = get_F_p_dvbar;
    thermal_params params = {m, Tfo, QMu,getParticle(pi)->getType() };
    functhermal.params = &params;
    using namespace qag_params;
    double result, error;
    // calculated yield
    int status = gsl_integration_qag (&functhermal, 0, 1, fEpsAbs,fEpsRel, fLimit,fKey ,fWorkspacethermal, &result, &error); 
  // Check status of integration and reduce accuracy if failing
  if (status) { 
status = gsl_integration_qag (&functhermal, 0, 1, fEpsAbs,1e-4, fLimit,GSL_INTEG_GAUSS15 ,fWorkspacethermal, &result, &error); 
    if (status) { 
status = gsl_integration_qag (&functhermal, 0, 1, 1e-4,1e-4, fLimit,GSL_INTEG_GAUSS15 ,fWorkspacethermal, &result, &error); 
      if (status) {      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }


    // updated yield
    getParticle(pi)->addN(result*nu);
    for (int i=0; i < getParticle(pi)->getNpbar(); i++){
      double Ebar = getParticle(pi)->getEbar(i);
    //  for (int j=0; j < grid_params::fNf ; j++){
    for(auto j: getParticle(pi)->fComponents) {
        getParticle(pi)->addFj((int) j,i, get_initial_Fj(j, Ebar,Tfo,QMu, getParticle(pi)->getType(), m, Cs2));
      }
    }
  gsl_integration_workspace_free (fWorkspacethermal);
}

//! Read in the particle data file and create TParticle_THERMINATOR class object, which is stored in fParticleData vector
void TFastReso_THERMINATOR::read_particles_data(string inputname, string comps, bool verbose){
  // Open the input file
  ifstream fInputFile(inputname);
  if (fInputFile.fail()) {
    cerr  << "\033[1mTFastReso_THERMINATOR.cpp::read_particles_data\033[0m : \033[1;31merror\033[0m :  input file does not exist "  << inputname << endl;
    exit(EXIT_FAILURE);
  }

  if (verbose){
  printf("#%9s\t%9s\t%9s\t%5s\t%4s\t%5s\t%s\t%s\t%s\t%s\t%s\t%s\t%9s\n",
      "Name", "Mass", "Gamma","Spin","Isospin","I3","Nq","Ns","Naq","Nas","Nc","Nac","MC");
  }

  string Name;
  double Mass, Gamma, Spin, Isospin, I3;
  double Nq, Ns, Naq, Nas,  Nc,Nac;
  int PDGCode;

  string line;
  while (getline(fInputFile, line)){
    //! check that there is the correct number of elements
    stringstream sstrc(line); int count=0; while (sstrc >> Name){count++;} if (count!=13) continue;
    stringstream sstr(line);
    sstr >> Name; 
    // check that the first non-whitespace character is not #, i.e. commented line
    if(Name.front()!='#'){
      sstr >> Mass >> Gamma >> Spin >> Isospin >> I3 >> Nq >> Ns >> Naq >> Nas >> Nc >> Nac >> PDGCode; 

  if (verbose){
      printf("%s\t%9f\t%9.5e\t%5.1f\t%5.1f\t%5.1f\t%g\t%g\t%g\t%g\t%g\t%g\t%9d\n",
          Name.c_str(), Mass, Gamma, Spin, Isospin, I3, Nq, Ns, Naq, Nas, Nc, Nac, PDGCode);
  }

      fIndexTable[Name] = (int) fParticleData.size();
      fIndexTablePDG[PDGCode] = (int) fParticleData.size();
      fParticleData.push_back(move(unique_ptr<TParticle_THERMINATOR>(new TParticle_THERMINATOR(Name, Mass,Gamma,Spin, Isospin, I3, Nq,Ns,Naq,Nas,Nc,Nac,PDGCode, comps))));
    }
  }


  if (verbose) {
    clog << "Number of particles read " << fParticleData.size() << endl;

  }
  fParticleNumber = (int) fParticleData.size();

}
//! Read the decay table and perform 2-body and 3-body decays
void TFastReso_THERMINATOR::do_decays(string inputname, bool verbose){

  // Precalcullated Clebsch Gordan coefficients for reweighting branching ratios
  std::map <std::string, double> fClebschGordanTable;
  string line;
  // Read a table of precalculated Clebsch Gordon coefficient
  // |<j1, m1, j2, m2 |J,M>|^2
  ifstream fCGCInputFile("CGC.dat");
  if (fCGCInputFile.fail()) {
    cerr  << "\033[1mTFastReso_THERMINATOR.cpp::do_decays\033[0m : \033[1;31merror\033[0m :  CGC.dat file does not exist "<< endl;
    exit(EXIT_FAILURE);
  }
  getline(fCGCInputFile, line);
  double j1,m1,j2,m2,J,M, CGC;
  char CGCname[100];
  while (getline(fCGCInputFile, line)){
    stringstream sstr(line);
    sstr >> j1 >> m1 >> j2 >> m2 >> J >> M >> CGC;
    sprintf(CGCname,"j1%2.1fm2%+3.1fj2%+2.1fm2%+3.1fJ%2.1fM%+3.1f", j1,m1,j2,m2,J,m1+m2);
    fClebschGordanTable[CGCname] = CGC;
  }
  fCGCInputFile.close();

  // Open the input file
  ifstream fInputFile(inputname);
  if (fInputFile.fail()) {
    cerr  << "\033[1mTFastReso_THERMINATOR.cpp::do_decays\033[0m : \033[1;31merror\033[0m :  input file does not exist "  << inputname<< endl;
    exit(EXIT_FAILURE);
  }
  // Print the header
  //
  if (verbose){
  printf("#%9s\t%9s\t%9s\t%9s\t%9s\t%9s\n",
      "Parent", "Child1", "Child2","Child3","BranchingRatio","ClebschGordanCoeff");
  }
  string Parent, Child1, Child2, Child3;
  double BranchingRatio;
  int ClebschGordanCoeff;

  int count_decays=0;
  int count_2decays=0;
  int count_3decays=0;

  while (getline(fInputFile, line)){
    //! check that there is the correct number of elements per line (5 for 2-body decay, 6 for three body decay), continue otherwise.
    stringstream sstrc(line); int count=0; while (sstrc >> Parent){count++;} if (not (count==5 or count==6)) continue;
    stringstream sstr(line);
    sstr >> Parent; 
    // check that the first non-whitespace character is not #, i.e. commented line
    if(Parent.front()!='#'){
      Child3="";
      if (count==5) {sstr >>  Child1 >> Child2 >> BranchingRatio >> ClebschGordanCoeff; }
      else if (count==6) {sstr >>  Child1 >> Child2 >> Child3>> BranchingRatio >> ClebschGordanCoeff; }
      //if (Child1==Child2 or Child1==Child3 or Child2==Child3) {
      // Print out the read-in coefficients
      //
  if (verbose){
      printf("%9s\t%9s\t%9s\t%9s\t%f\t%d\n",
          Parent.c_str(), Child1.c_str(), Child2.c_str(), Child3.c_str(), BranchingRatio, ClebschGordanCoeff);
  }
      //}
      count_decays++;
      (count==5) ? count_2decays++ : count_3decays++;

      if (count==5)
      {
        double Ma = getParticle(Parent)->getM();
        double Mb = getParticle(Child1)->getM();
        double Mc = getParticle(Child2)->getM();
        // Skip decays violating mass inequality (only possible by including resonance widths)
        if (Ma < Mb+Mc) {
          //cerr  << "\033[1mTFastReso_THERMINATOR.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc << endl;
          continue;
        }
        // Skip decays whose parents have decay width less than 10KeV 
        // If ClebschGodranCoeff flag is non-zero, rescaling the branching ratio.
        if (ClebschGordanCoeff){
          double J = getParticle(Parent)->getIsospin();
          double M = getParticle(Parent)->getI3();
          double j1 = getParticle(Child1)->getIsospin();
          double m1 = getParticle(Child1)->getI3();
          double j2 = getParticle(Child2)->getIsospin();
          double m2 = getParticle(Child2)->getI3();
          double spin1 = getParticle(Child1)->getSpin();
          double spin2 = getParticle(Child2)->getSpin();
          double QS1 = getParticle(Child1)->getQS();
          double QS2 = getParticle(Child2)->getQS();
          double QC1 = getParticle(Child1)->getQC();
          double QC2 = getParticle(Child2)->getQC();
          if (M != m1+m2) {
            cerr  << "\033[1mTFastReso_THERMINATOR.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle I3 doesnt' match the that of childern! " << M << " != " << m1+m2 << endl;
            continue;
          }


          sprintf(CGCname,"j1%2.1fm2%+3.1fj2%+2.1fm2%+3.1fJ%2.1fM%+3.1f", j1,m1,j2,m2,J,m1+m2);
          BranchingRatio *= fClebschGordanTable[CGCname];
          //cout << " Multiplying branching ratio by Clebsch-Gordan coefficient " <<  fClebschGordanTable[CGCname]<<endl;
          // If decay products belong to the same isospin multiplet, bultiply branching ration by two, to account for the process
          // of interchanging particles, c.f. pi^0+pi^+ and pi^+ + pi^0 vs pi^0 + K^+ and pi^+ K^0
          if (abs(Mb-Mc) <0.01 and abs(spin1-spin2) <0.01 and abs(QC1-QC2) < 1e-6 and abs(QS1-QS2) < 1e-6 and abs(QC1-QC2)<1e-6 and abs(m1-m2)>0.01)
          {
            //cout << " Multiplying branching ratio by 2 "<<endl;
            BranchingRatio *= 2;

          }
        }    
        //double GammaA = getParticle(Parent)->getGamma();
        //double QS = getParticle(Parent)->getQS();
        //double QS1 = getParticle(Child1)->getQS();
        //double QS2 = getParticle(Child2)->getQS();
        //if (QS!=QS1+QS2) {
        //  cerr  << "\033[1mTFastReso_THERMINATOR.cpp::do_decays\033[0m : \033[1;31merror\033[0m : strangeness changing weak decay" << endl;
        //  cerr << Parent << " " << Child1 << " " <<  Child2 << endl;
        //  cerr << QS << " " << QS1 << " " <<  QS2 << endl;
        //  continue;
        //}
        do_2bodydecay(getParticle(Parent), getParticle(Child1), getParticle(Child2), BranchingRatio);
        if (Child1!=Child2){ do_2bodydecay(getParticle(Parent), getParticle(Child2), getParticle(Child1), BranchingRatio); }

        //if (Child1=="pi0139plu" or Child2=="pi0139plu") { getParticle("pi0139plu")->print("pion_"+to_string(fIndexTable[Parent]) +"_"+to_string(fIndexTable[Child1])+"_"+to_string(fIndexTable[Child2])); }
      }
      else if (count==6)
      {
        double Ma = getParticle(Parent)->getM();
        double Mb = getParticle(Child1)->getM();
        double Mc = getParticle(Child2)->getM();
        double Md = getParticle(Child3)->getM();
        if (Ma < Mb+Mc+Md) {
          cerr  << "\033[1mTFastReso_THERMINATOR.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc+Md << endl;
          continue;
        }
        //double GammaA = getParticle(Parent)->getGamma();
        //double QS = getParticle(Parent)->getQS();
        //double QS1 = getParticle(Child1)->getQS();
        //double QS2 = getParticle(Child2)->getQS();
        //double QS3 = getParticle(Child3)->getQS();
        //if (QS!=QS1+QS2+QS3) {
        //  cerr  << "\033[1mTFastReso_THERMINATOR.cpp::do_decays\033[0m : \033[1;31merror\033[0m : strangeness changing weak decay" << endl;
        //  cerr << Parent << " " << Child1 << " " <<  Child2 << " " << Child3 << endl;
        //  cerr << QS << " " << QS1 << " " <<  QS2  << " " << QS3<< endl;
        //  continue;
        //}
        do_3bodydecay(getParticle(Parent), getParticle(Child1), getParticle(Child2),getParticle(Child3), BranchingRatio);
        if (Child1!=Child2){ do_3bodydecay(getParticle(Parent), getParticle(Child2), getParticle(Child1),getParticle(Child3), BranchingRatio); }
        if (Child1!=Child3 and Child2!=Child3){ do_3bodydecay(getParticle(Parent), getParticle(Child3), getParticle(Child1),getParticle(Child2), BranchingRatio); }

        // if (Child1=="pi0139plu" or Child2=="pi0139plu" or Child3=="pi0139plu") {getParticle("pi0139plu")->print("pion_"+to_string(fIndexTable[Parent]) +"_"+to_string(fIndexTable[Child1])+"_"+to_string(fIndexTable[Child2])+"_"+to_string(fIndexTable[Child3])); }

      }


    }
  }


  if (verbose) {
    clog << "Tota number of decays read " << count_decays << endl;
    clog << "Number of 2-body decays read " << count_2decays << endl;
    clog << "Number of 3-body decays read " << count_3decays << endl;
  }
}


