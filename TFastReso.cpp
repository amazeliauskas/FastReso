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
#include "TFastReso.h"
#include <fstream>
#include <sstream>
#include <math.h>


using namespace std;

//! Eabc -- energy of b in restframe of a in a->b+c decay
double get_Eabc(double ma,double mb,double mc) { return (ma*ma+mb*mb-mc*mc)/(2*ma);}
//! particle momenta
double get_p(double E, double m) { return sqrt(E*E-m*m);}
//! pabc -- momentum of b or c in restframe of a in a->b+c decay
double get_pabc(double ma,double mb,double mc) { return get_p(get_Eabc(ma,mb,mc),mb);}

//! returns thermal Bose-Einstein, Fermi-Dirac or Boltzmann distribution
double get_thermal_F(double Ebar, double T, double QMu, EParticleType type);
//! returns derivative of thermal Bose-Einstein, Fermi-Dirac or Boltzmann distribution
double get_thermal_dF(double Ebar, double T, double QMu, EParticleType type);


//! Initialize f_j on the freeze-out surface
//! Here you can add additional fj's
double get_initial_Fj(int j, double Ebar, double T, double QMu, EParticleType type, double m, double Cs2)
{
  switch(j){
    case 0: //feq 1
      return get_thermal_F(Ebar,T, QMu, type);
      break;
    case 1: //feq 2
      return get_thermal_F(Ebar,T, QMu, type);
      break;
    case 2: //fshear 1
      return get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 3: //fshear 2
      return get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 4: //fshear 3
      return get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 5: //fbulk 1
      return get_thermal_dF(Ebar,T, QMu, type)*(1./3.*m*m/T/Ebar-Ebar/T*(1./3.-Cs2));
      break;
    case 6: //fbulk 2
      return get_thermal_dF(Ebar,T, QMu, type)*(1./3.*m*m/T/Ebar-Ebar/T*(1./3.-Cs2));
      break;
    case 7: //ftemperature. 1
      return get_thermal_dF(Ebar,T, QMu, type)*Ebar/T;
      break;
    case 8: //ftemperature 2
      return get_thermal_dF(Ebar,T, QMu, type)*Ebar/T;
      break;
    case 9: //fvelocity  1
      return get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 10: //fvelocity 2
      return get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 11: //fvelocity 3
      return get_thermal_dF(Ebar,T, QMu, type);
      break;
    default:
      cerr  << "\033[1mTFastReso.cpp::get_initial_Fj\033[0m : \033[1;31merror\033[0m :  wrong index "  << j << endl;
      exit(EXIT_FAILURE);
  }
}
//! Transformation rule factor A for Fj in a 2-body decay
//! Here you can add additional rules for fj's
double get_factor_Aj(int j, double w, double Eabc, double pabc, double Ebar, double pbar, double Ma, double Mb)
{

  double A= (Ma*Eabc/Mb/Mb-w*Ma*pabc/Mb/Mb*Ebar/pbar);
  double Ew= (Ma*Eabc*Ebar/Mb/Mb-w*Ma*pabc*pbar/Mb/Mb);
  switch(j){
    case 0: // feq 1
      return A;
      break;
    case 1: // feq 2
      return Ew/Ebar;
      break;
    case 2: // fshear 1
      return A*(2.5*A*A-1.5*(Ew*Ew-Ma*Ma)/pbar/pbar);
      break;
    case 3: // fshear 2
      return A*(Ew*Ew-Ma*Ma)/pbar/pbar;
      break;
    case 4: // fshear 3
      return (1.5*A*A-0.5*(Ew*Ew-Ma*Ma)/pbar/pbar)*Ew/Ebar;
      break;
    case 5: // fbulk 1
      return A;
      break;
    case 6: // fbulk 2
      return Ew/Ebar;
      break;
    case 7: // ftemperature 1
      return A;
      break;
    case 8: // ftemperature 2
      return Ew/Ebar;
      break;
    case 9: // fvelocity 1
      return 1.5*A*A-0.5*(Ew*Ew-Ma*Ma)/pbar/pbar;
      break;
    case 10: // fvelocity 2
      return (Ew*Ew-Ma*Ma)/pbar/pbar;
      break;
    case 11: // fvelocity 3
      return A*Ew/Ebar;
      break;
    default:
      cerr  << "\033[1mTFastReso.cpp::get_factor_Aj\033[0m : \033[1;31merror\033[0m :  wrong index "  << j << endl;
      exit(EXIT_FAILURE);
  }
}

//! Parameters for integration of thermal distribution
struct thermal_params{
  double m;
  double Tfo;
  double QMu; 
  EParticleType type;
};
//! Thermal distribution function for integration over vbar=[0,1]
double get_F_p_dvbar(double vbar, void * p) {
  double &m           = ((struct thermal_params *) p)->m;
  double &Tfo           = ((struct thermal_params *) p)->Tfo;
  double &QMu           = ((struct thermal_params *) p)->QMu;
  EParticleType &type           = ((struct thermal_params *) p)->type;
  double Ebar=m/sqrt(1-vbar*vbar);
  double pbar=vbar*Ebar;
  double F_p = get_thermal_F(Ebar,Tfo,QMu, type);
  return F_p*pbar*pbar*4*M_PI*m/pow(2*M_PI*sqrt(1-vbar*vbar),3);
}
//! Parameters for integration of two body decay function
struct twobody_params{
  double Ebar;
  double Ma;
  double Mb;
  double Mc;
  TParticle * Father;
  int index;     
};               
//! Integrand of two body decay formula for f_j components
double get_F_pu_dw(double w, void * p) {
  double &Ebar           = ((struct twobody_params *) p)->Ebar;
  double &Ma             = ((struct twobody_params *) p)->Ma;
  double &Mb             = ((struct twobody_params *) p)->Mb;
  double &Mc             = ((struct twobody_params *) p)->Mc;
  TParticle *Father      = ((struct twobody_params *) p)->Father;
  int &index             = ((struct twobody_params *) p)->index;

  double Eabc = get_Eabc(Ma,Mb,Mc);
  double pabc = get_p(Eabc, Mb);  
  //double pbar = GSL_MAX(get_p(Ebar, Mb),DBL_EPSILON);  
  double pbar = get_p(Ebar, Mb);// ,DBL_EPSILON);  
  double Ebar_old = Ebar*Eabc*Ma/Mb/Mb-w*pbar*pabc*Ma/Mb/Mb;
  double F_old=Father->get_Fj(index, Ebar_old);
  double A = get_factor_Aj(index, w, Eabc, pabc, Ebar, pbar, Ma, Mb);
  return F_old*A;
}

//! Parameters for three body decay mass integral
struct threebody_params{
  gsl_function  *twobody_function;
  double Mc;
  double Md;
  gsl_integration_workspace *twobody_workspace;
};

//! Integrand for three body decay normalization
double get_B_dm(double mct, void * p) {
  gsl_function *twobody_function = ((struct threebody_params *) p)->twobody_function;
  double &Ma = ((struct twobody_params *) twobody_function->params)->Ma;
  double &Mb = ((struct twobody_params *) twobody_function->params)->Mb;
  double &Mc             = ((struct threebody_params *) p)->Mc;
  double &Md             = ((struct threebody_params *) p)->Md;
  return get_pabc(Ma,Mb,mct)*get_pabc(mct, Mc, Md);
}
//! Integrand for three body decay function
double get_F_pu_dm(double mct, void * p) {
  gsl_function *twobody_function = ((struct threebody_params *) p)->twobody_function;
  ((struct twobody_params *) twobody_function->params)->Mc = mct;
  gsl_integration_workspace *twobody_workspace = ((struct threebody_params *) p)->twobody_workspace;
  double &Ma = ((struct twobody_params *) twobody_function->params)->Ma;
  double &Mb = ((struct twobody_params *) twobody_function->params)->Mb;
  double &Mc             = ((struct threebody_params *) p)->Mc;
  double &Md             = ((struct threebody_params *) p)->Md;
  double result,error;

  using namespace qag_params;
  gsl_integration_qag (twobody_function, -1, 1, fEpsAbs,fEpsRel, fLimit,fKey ,twobody_workspace, &result, &error); 
  return result*get_pabc(Ma,Mb,mct)*get_pabc(mct, Mc,Md);
}



//! Initialization of thermal particle distribution functions
void TFastReso::do_thermal(double Tfo, double MuB, double MuI3, double MuS, double MuC, double Cs2){

  using namespace qag_params;
  gsl_integration_workspace * fWorkspacethermal ;
  fWorkspacethermal = gsl_integration_workspace_alloc (fLimit);
  for (std::vector<unique_ptr<TParticle>>::iterator it = fParticleData.begin() ; it != fParticleData.end(); ++it){
    // don't initialize photon
    if ((*it)->getName()=="gam000zer") continue;
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

//! Read in the particle data file and create TParticle class object, which is stored in fParticleData vector
void TFastReso::read_particles_data(string inputname){
  // Open the input file
  ifstream fInputFile(inputname);
  if (fInputFile.fail()) {
    cerr  << "\033[1mTFastReso.cpp::read_particles_data\033[0m : \033[1;31merror\033[0m :  input file does not exist "  << inputname << endl;
    exit(EXIT_FAILURE);
  }

  // printf("#%9s\t%9s\t%9s\t%5s\t%4s\t%5s\t%s\t%s\t%s\t%s\t%s\t%s\t%9s\n",
  //      "Name", "Mass", "Gamma","Spin","Isospin","I3","Nq","Ns","Naq","Nas","Nc","Nac","MC");

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
      //   printf("%s\t%9f\t%9.5e\t%5.1f\t%5.1f\t%5.1f\t%g\t%g\t%g\t%g\t%g\t%g\t%9d\n",
      //       Name.c_str(), Mass, Gamma, Spin, Isospin, I3, Nq, Ns, Naq, Nas, Nc, Nac, MC);

      fIndexTable[Name] = fParticleData.size();
      fParticleData.push_back(move(unique_ptr<TParticle>(new TParticle(Name, Mass,Gamma,Spin, Isospin, I3, Nq,Ns,Naq,Nas,Nc,Nac,PDGCode))));
    }
  }
  clog << "Number of particles read " << fParticleData.size() << endl;
}

//! Read the decay table and perform 2-body and 3-body decays
void TFastReso::do_decays(string inputname){

  string line;
  // Read a table of precalculated Clebsch Gordon coefficient
  // |<j1, m1, j2, m2 |J,M>|^2
  ifstream fCGCInputFile("CGC.dat");
  if (fCGCInputFile.fail()) {
    cerr  << "\033[1mTFastReso.cpp::do_decays\033[0m : \033[1;31merror\033[0m :  CGC.dat file does not exist "<< endl;
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
    cerr  << "\033[1mTFastReso.cpp::do_decays\033[0m : \033[1;31merror\033[0m :  input file does not exist "  << inputname<< endl;
    exit(EXIT_FAILURE);
  }
 // Print the header
  printf("#%9s\t%9s\t%9s\t%9s\t%9s\t%9s\n",
      "Father", "Child1", "Child2","Child3","BranchingRatio","ClebschGordanCoeff");

  string Father, Child1, Child2, Child3;
  double BranchingRatio;
  int ClebschGordanCoeff;

  int count_decays=0;
  int count_2decays=0;
  int count_3decays=0;

  while (getline(fInputFile, line)){
    //! check that there is the correct number of elements per line (5 for 2-body decay, 6 for three body decay), continue otherwise.
    stringstream sstrc(line); int count=0; while (sstrc >> Father){count++;} if (not (count==5 or count==6)) continue;
    stringstream sstr(line);
    sstr >> Father; 
    // check that the first non-whitespace character is not #, i.e. commented line
    if(Father.front()!='#'){
      Child3="";
      if (count==5) {sstr >>  Child1 >> Child2 >> BranchingRatio >> ClebschGordanCoeff; }
      else if (count==6) {sstr >>  Child1 >> Child2 >> Child3>> BranchingRatio >> ClebschGordanCoeff; }
      //if (Child1==Child2 or Child1==Child3 or Child2==Child3) {
      // Print out the read-in coefficients
      printf("%9s\t%9s\t%9s\t%9s\t%f\t%d\n",
          Father.c_str(), Child1.c_str(), Child2.c_str(), Child3.c_str(), BranchingRatio, ClebschGordanCoeff);
      //}
      count_decays++;
      (count==5) ? count_2decays++ : count_3decays++;
      if (count==5)
      {
        double Ma = getParticle(Father)->getM();
        double GammaA = getParticle(Father)->getGamma();
        double Mb = getParticle(Child1)->getM();
        double Mc = getParticle(Child2)->getM();
        // Skip decays violating mass inequality (only possible by including resonance widths)
        if (Ma < Mb+Mc) {
          cerr  << "\033[1mTFastReso.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc << endl;
          continue;
        }
        //// Skip decays whose parents have decay width less than 10KeV 
        //if (GammaA <  10e-6) {
        //  cerr  << "\033[1mTFastReso.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle width smaller than 10 KeV ! Gamma = " << GammaA << " GeV" << endl;
        //  continue;
        //}
        // If ClebschGodranCoeff flag is non-zero, rescaling the branching ratio.
        if (ClebschGordanCoeff){
          double J = getParticle(Father)->getIsospin();
          double M = getParticle(Father)->getI3();
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
            cerr  << "\033[1mTFastReso.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle I3 doesnt' match the that of childern! " << M << " != " << m1+m2 << endl;
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
        do_2bodydecay(Father, Child1, Child2, BranchingRatio);
        if (Child1!=Child2){ do_2bodydecay(Father, Child2, Child1, BranchingRatio); }

        //if (Child1=="pi0139plu" or Child2=="pi0139plu") { getParticle("pi0139plu")->print("pion_"+to_string(fIndexTable[Father]) +"_"+to_string(fIndexTable[Child1])+"_"+to_string(fIndexTable[Child2])); }
      }
      else if (count==6)
      {
        double Ma = getParticle(Father)->getM();
        double Mb = getParticle(Child1)->getM();
        double Mc = getParticle(Child2)->getM();
        double Md = getParticle(Child3)->getM();
        if (Ma < Mb+Mc+Md) {
          cerr  << "\033[1mTFastReso.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc+Md << endl;
          continue;
        }

        do_3bodydecay(Father, Child1, Child2,Child3, BranchingRatio);
        if (Child1!=Child2){ do_3bodydecay(Father, Child2, Child1,Child3, BranchingRatio); }
        if (Child1!=Child3 and Child2!=Child3){ do_3bodydecay(Father, Child3, Child1,Child2, BranchingRatio); }

        // if (Child1=="pi0139plu" or Child2=="pi0139plu" or Child3=="pi0139plu") {getParticle("pi0139plu")->print("pion_"+to_string(fIndexTable[Father]) +"_"+to_string(fIndexTable[Child1])+"_"+to_string(fIndexTable[Child2])+"_"+to_string(fIndexTable[Child3])); }

      }


    }
  }


  clog << "Tota number of decays read " << count_decays << endl;
  clog << "Number of 2-body decays read " << count_2decays << endl;
  clog << "Number of 3-body decays read " << count_3decays << endl;

}

//! Initialize thermal distribution
double get_thermal_F(double Ebar, double T, double QMu, EParticleType type)
{
  switch(type){
    case EParticleType::kBoson:
      return 1.0/(exp((Ebar-QMu)/T)-1);
      break;
    case EParticleType::kFermion: 
      return 1.0/(exp((Ebar-QMu)/T)+1);
      break;
    case EParticleType::kBoltzman: 
      return exp(-(Ebar-QMu)/T);
      break;
    default: 
      cerr << "\033[1mTFastReso.cpp::get_thermal_F\033[0m : \033[1;31merror\033[0m : unrecognined particle type!" <<endl;
      exit(EXIT_FAILURE);}
}
//! Initialize  derivative of t hermal distribution (With a minus sign)
double get_thermal_dF(double Ebar, double T, double QMu, EParticleType type)
{
  switch(type){
    case EParticleType::kBoson:
      return 1.0/(exp((Ebar-QMu)/T)-1)*(1+1.0/(exp((Ebar-QMu)/T)-1));
      break;
    case EParticleType::kFermion: 
      return 1.0/(exp((Ebar-QMu)/T)+1)*(1-1.0/(exp((Ebar-QMu)/T)+1));
      break;
    case EParticleType::kBoltzman: 
      return exp(-(Ebar-QMu)/T);
      break;
    default: 
      cerr << "\033[1mTFastReso.cpp::get_thermal_F_pu\033[0m : \033[1;31merror\033[0m : unrecognined particle type!" <<endl;
      exit(EXIT_FAILURE);}
}


//! Perform two body decay integral: Father->Child1 + Child2
void TFastReso::do_2bodydecay(string Father, string Child1, string Child2, double BranchingRatio){
  using namespace qag_params;

  gsl_integration_workspace * fWorkspace2body ;
  fWorkspace2body = gsl_integration_workspace_alloc (fLimit);
  getParticle(Father)->lock();
  double Ma = getParticle(Father)->getM();
  double Mb = getParticle(Child1)->getM();
  double Mc = getParticle(Child2)->getM();
  double nua = getParticle(Father)->getNu();
  double nub = getParticle(Child1)->getNu();
  double Ebar=0.0;
  getParticle(Child1)->clean_buffer();
  if (Ma < Mb+Mc) {
    cerr  << "\033[1mTFastReso.cpp::do_2bodydecays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc << endl;
    exit(EXIT_FAILURE);
  }
  gsl_function func2body;
  func2body.function = get_F_pu_dw;
  twobody_params params = {Ebar, Ma, Mb, Mc, getParticle(Father), 0};
  func2body.params = &params;
  double error, result ;
  int sym=1 ;
  if (Child1==Child2 ){ sym=2; } else { sym=1; }
  double Nfather = getParticle(Father)->getN();
  //! add fraction of the total yield to child
  getParticle(Child1)->addN(BranchingRatio*Nfather*sym);
  for (int i=0; i < getParticle(Child1)->getNpbar(); i++){
    double Ebar =getParticle(Child1)->getEbar(i);
    params.Ebar =Ebar;

    for (int j=0; j < grid_params::fNf ; j++){
      params.index =j;
      gsl_integration_qag (&func2body, -1, 1, fEpsAbs,fEpsRel, fLimit,fKey ,fWorkspace2body, &result, &error); 
      double fac = sym*BranchingRatio*nua/nub*pow(Ma/Mb,2)/2;
      getParticle(Child1)->addFj(j,i, fac*result);

    }
  }
  //if (Child1=="pi0139plu"  ) { getParticle(Child1)->print_buffer(Child1+"_"+Father+"_"+Child2);}
  gsl_integration_workspace_free (fWorkspace2body);

}

//! Perform three body decay integral Father->Child1 + Child2 + Child3
void TFastReso::do_3bodydecay(string Father, string Child1, string Child2,string Child3, double BranchingRatio){
  using namespace qag_params;
  gsl_integration_workspace * fWorkspace2body ;
  gsl_integration_workspace * fWorkspace3body ;
  fWorkspace2body = gsl_integration_workspace_alloc (fLimit);
  fWorkspace3body = gsl_integration_workspace_alloc (fLimit);
  getParticle(Father)->lock();
  double Ma = getParticle(Father)->getM();
  double Mb = getParticle(Child1)->getM();
  double Mc = getParticle(Child2)->getM();
  double Md = getParticle(Child3)->getM();
  double nua = getParticle(Father)->getNu();
  double nub = getParticle(Child1)->getNu();
  double Ebar=0.0;

  getParticle(Child1)->clean_buffer();
  if (Ma < Mb+Mc+Md) {
    cerr  << "\033[1mTFastReso.cpp::do_3bdoydecays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc+Md << endl;
    exit(EXIT_FAILURE);
  }
  gsl_function func2body;
  gsl_function func3body;
  func2body.function = get_F_pu_dw;
  func3body.function = get_B_dm;
  twobody_params params2 = {Ebar, Ma, Mb, Mc+Md, getParticle(Father), 0};
  func2body.params = &params2;
  threebody_params params3 = {&func2body, Mc,Md, fWorkspace2body};
  func3body.params = &params3;
  double error, result ;
  gsl_integration_qag (&func3body, Mc+Md, Ma-Mb, fEpsAbs,fEpsRel, fLimit,fKey ,fWorkspace3body, &result, &error); 
  double B3body = result;
  func3body.function = get_F_pu_dm;

  int sym = 1 ;
  if (Child1==Child2 and Child2==Child3 ){ sym=3; } else if (Child1==Child2 or Child1==Child3 ){ sym=2;} else { sym=1; }
  double Nfather = getParticle(Father)->getN();
  // add fraction of yield to child
  getParticle(Child1)->addN(BranchingRatio*Nfather*sym);
  for (int i=0; i < getParticle(Child1)->getNpbar(); i++){
    double Ebar =getParticle(Child1)->getEbar(i);
    params2.Ebar =Ebar;
    for (int j=0; j < grid_params::fNf ; j++){
      params2.index =j;
      gsl_integration_qag (&func3body, Mc+Md, Ma-Mb, fEpsAbs,fEpsRel, fLimit,fKey ,fWorkspace3body, &result, &error); 
      double fac = sym*BranchingRatio*nua/nub*pow(Ma/Mb,2)/2/B3body;
      getParticle(Child1)->addFj(j,i, fac*result);

    }
  }


  //if (Child1=="pi0139plu" ) { getParticle(Child1)->print_buffer(Child1+"_"+Father+"_"+Child2+"_"+Child3);}
  gsl_integration_workspace_free (fWorkspace2body);
  gsl_integration_workspace_free (fWorkspace3body);
}


