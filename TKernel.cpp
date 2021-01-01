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
#include <sstream>
#include <gsl/gsl_integration.h>
#include "TKernel.h"

using namespace std;
// struct for eta and phi numerical integrals
struct dK_dphideta_params{
  double phi,pT,ur;
  int index;
  TKernel * kernel_class;
};
struct dK_dphi_params{
  gsl_function  *Fdphideta;
  gsl_integration_workspace *weta ;
  double etamax;
};
double dK_dphideta(double eta, void * p) {
  double &phi          = ((struct dK_dphideta_params *) p)->phi;
  double &pT           = ((struct dK_dphideta_params *) p)->pT;
  double &ur           = ((struct dK_dphideta_params *) p)->ur;
  int &index           = ((struct dK_dphideta_params *) p)->index;
  TKernel * kernel_class  = ((struct dK_dphideta_params *) p)->kernel_class;
  //double dK_mu_detadphi[8];
  //kernel_class->get_K_mu_detadphi(eta, phi, pT, ur, dK_mu_detadphi);
  double dK_mu_detadphi = kernel_class->get_K_mu_detadphi_i(eta, phi, pT, ur,  index);
  return 2*dK_mu_detadphi;
}

double dK_dphi(double phi, void * p) {
  gsl_function *Fdphideta          = ((struct dK_dphi_params *) p)->Fdphideta;
  gsl_integration_workspace *weta  = ((struct dK_dphi_params *) p)->weta;
  double &etamax                   = ((struct dK_dphi_params *) p)->etamax;
  int &index           = ((struct dK_dphideta_params *) ((struct dK_dphideta_params *) Fdphideta->params))->index;
  ((struct dK_dphideta_params *) Fdphideta->params)->phi = phi;
  double result, error;
  using namespace qag_params;
  int status =  gsl_integration_qag (Fdphideta, 0.0, etamax,  fEpsAbs,fEpsRel,fLimit,fKey,weta, &result, &error); 
  // Check status of integration and reduce accuracy if failing
  if (status) { 
    status =  gsl_integration_qag (Fdphideta, 0.0, etamax,  fEpsAbs,1e-4,fLimit,GSL_INTEG_GAUSS15,weta, &result, &error); 
    if (status) { 
      status =  gsl_integration_qag (Fdphideta, 0.0, etamax,  1e-4,1e-4,fLimit,GSL_INTEG_GAUSS15,weta, &result, &error); 
      //if (status) {      std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;

      if (status) {      std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << "index = " << index  << " " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }


  return 2*result;
}

TKernel::TKernel(TParticle *particle) : fParticle(particle) {
  initialize();
}
TKernel::~TKernel(){
  gsl_integration_workspace_free (fWeta);
  gsl_integration_workspace_free (fWphi);
}
void TKernel::initialize(){
  fWphi = gsl_integration_workspace_alloc (qag_params::fLimit);
  fWeta = gsl_integration_workspace_alloc (qag_params::fLimit);
}

// pT = pT ur =sinh(chi)
double TKernel::get_K_mu_detadphi_i(double eta, double phi,double pT, double ur, int index) {

  double fMb = fParticle->getM();
  int fNu = fParticle->getNu();
  double mT = sqrt(fMb*fMb+pT*pT);
  double utau = sqrt(1+ur*ur); 
  double Ebar=mT*utau*cosh(eta)-pT*ur*cos(phi);
  double pbar = sqrt(Ebar*Ebar-fMb*fMb);
  switch (index) {
    case 0: 
      {
        double feq1 = fNu*fParticle->get_Fj(0, Ebar)/pbar;
        double feq2 = fNu*fParticle->get_Fj(1, Ebar)/pbar;
        return feq1*mT*cosh(eta) + (feq2-feq1)*Ebar*utau;
      } 
      break;
    case 1:
      {
        double feq1 = fNu*fParticle->get_Fj(0, Ebar)/pbar;
        double feq2 = fNu*fParticle->get_Fj(1, Ebar)/pbar;
        return feq1*cos(phi)*pT + (feq2-feq1)*Ebar*ur;
      }
      break;
    case 2:
      {
        double fshear1 = fNu*fParticle->get_Fj(2, Ebar)/pbar/pbar/pbar;
        double fshear2 = fNu*fParticle->get_Fj(3, Ebar)/pbar/pbar/pbar;
        double fshear3 = fNu*fParticle->get_Fj(4, Ebar)/pbar/pbar/pbar;
        return (fshear1*cosh(eta)*mT + (fshear3-fshear1)*Ebar*utau)*
          (mT*mT*sinh(eta)*sinh(eta)-pow(mT*ur*cosh(eta)-pT*utau*cos(phi), 2)  )
          +(fshear2-fshear1)*2.0/5.0*pbar*pbar*ur*(mT*ur*cosh(eta)-pT*utau*cos(phi)) ;
      } 
      break;
    case 3:
      {
        double fshear1 = fNu*fParticle->get_Fj(2, Ebar)/pbar/pbar/pbar;
        double fshear2 = fNu*fParticle->get_Fj(3, Ebar)/pbar/pbar/pbar;
        double fshear3 = fNu*fParticle->get_Fj(4, Ebar)/pbar/pbar/pbar;
        return (fshear1*cos(phi)*pT + (fshear3-fshear1)*Ebar*ur)*
          (mT*mT*sinh(eta)*sinh(eta)-pow(mT*ur*cosh(eta)-pT*utau*cos(phi), 2)  )
          +(fshear2-fshear1)*2.0/5.0*pbar*pbar*utau*(mT*ur*cosh(eta)-pT*utau*cos(phi));
      } 
      break;
    case 4:
      {
        double fshear1 = fNu*fParticle->get_Fj(2, Ebar)/pbar/pbar/pbar;
        double fshear2 = fNu*fParticle->get_Fj(3, Ebar)/pbar/pbar/pbar;
        double fshear3 = fNu*fParticle->get_Fj(4, Ebar)/pbar/pbar/pbar;
        return (fshear1*cosh(eta)*mT+ (fshear3-fshear1)*Ebar*utau)*
          (pT*pT*sin(phi)*sin(phi)-pow(mT*ur*cosh(eta)-pT*utau*cos(phi), 2)  )
          +(fshear2-fshear1)*2.0/5.0*pbar*pbar*ur*(mT*ur*cosh(eta)-pT*utau*cos(phi));
      } 
      break;
    case 5:
      {
        double fshear1 = fNu*fParticle->get_Fj(2, Ebar)/pbar/pbar/pbar;
        double fshear2 = fNu*fParticle->get_Fj(3, Ebar)/pbar/pbar/pbar;
        double fshear3 = fNu*fParticle->get_Fj(4, Ebar)/pbar/pbar/pbar;
        return (fshear1*cos(phi)*pT + (fshear3-fshear1)*Ebar*ur)*
          (pT*pT*sin(phi)*sin(phi)-pow(mT*ur*cosh(eta)-pT*utau*cos(phi), 2)  )
          +(fshear2-fshear1)*2.0/5.0*pbar*pbar*utau*(mT*ur*cosh(eta)-pT*utau*cos(phi));
      } 
      break;
    case 6:
      {
        double fbulk1 = fNu*fParticle->get_Fj(5, Ebar)/pbar;
        double fbulk2 = fNu*fParticle->get_Fj(6, Ebar)/pbar;
        return fbulk1*cosh(eta)*mT + (fbulk2-fbulk1)*Ebar*utau;
      } 
      break;
    case 7:
      {
        double fbulk1 = fNu*fParticle->get_Fj(5, Ebar)/pbar;
        double fbulk2 = fNu*fParticle->get_Fj(6, Ebar)/pbar;
        return fbulk1*cos(phi)*pT + (fbulk2-fbulk1)*Ebar*ur;
      } 
      break;
    default:
      cerr  << "\033[1mTKernel.cpp::get_K_mu_detadphi_i\033[0m : \033[1;31merror\033[0m :  index out of bounds i =  "  << index << endl;
      exit(EXIT_FAILURE);

  }

}

double TKernel::get_K_mu_qag(double pT, double ur, int index) {
  using namespace qag_params;
  double K_mu =0;
  gsl_function Fdphi;
  Fdphi.function = dK_dphi;
  gsl_function Fdphideta;
  Fdphideta.function = dK_dphideta;
  dK_dphideta_params dphideta_params = {0.0, pT,ur,index, this};
  Fdphideta.params = &dphideta_params;
  dK_dphi_params dphi_params = {&Fdphideta,fWeta,fEtaMax};
  Fdphi.params = &dphi_params;
  double error;
  int status;

  status = gsl_integration_qag (&Fdphi, 0.0, M_PI,fEpsAbs,fEpsRel,fLimit,fKey ,fWphi, &K_mu, &error); 

  if (status) { 
    status = gsl_integration_qag (&Fdphi, 0.0, M_PI,fEpsAbs,1e-4,fLimit,GSL_INTEG_GAUSS15 ,fWphi, &K_mu, &error); 
    if (status) { 
      status = gsl_integration_qag (&Fdphi, 0.0, M_PI,1e-4,1e-4,fLimit,GSL_INTEG_GAUSS15 ,fWphi, &K_mu, &error); 
      if (status) {      std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << "index = " << index<< " " << K_mu << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << K_mu<< "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << K_mu << "+-" << error <<std::endl;
    }
  }
  return K_mu;
}

void TKernel::print(string tag, string columns){

  stringstream sstr(columns);
  string col;
  bool IsKeq=false;
  bool IsKshear=false;
  bool IsKbulk=false;
  while (getline(sstr, col, ' ')){
    if (col=="Keq")         { IsKeq=true; }
    else if (col=="Kshear") { IsKshear=true; }
    else if (col=="Kbulk")  { IsKbulk=true; }
    else { std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : unkwnown label " << col  <<std::endl; }

  }

  FILE * pFile;
  string fname = tag+"_Kj.out";
  pFile = fopen (fname.c_str(),"w");
  fprintf(pFile,"#%15s\t%s", "1:pT [GeV]","2:ur"); 
  int index = 3;
  if (IsKeq) {
    fprintf(pFile,"\t%d:%s\t%d:%s", index, "Keq 1", index+1, "Keq 2");
    index +=2;
  }
  if (IsKshear) {
    fprintf(pFile,"\t%d:%s\t%d:%s\t%d:%s\t%d:%s",
        index, "Kshear 1", index+1, "Kshear 2" , index+2,  "Kshear 3", index+3, "Kshear 4");
    index +=4;
  }
  if (IsKbulk) {
    fprintf(pFile,"\t%d:%s\t%d:%s", index,  "Kbulk 1", index+1, "Kbulk 2");
    index +=2;
  }

  fprintf(pFile,"\n");
  for(int j = 0; j <grid_params::fNur; j++){
    double ur = tan(atan(grid_params::fUrMax)*(j)/(grid_params::fNur-1)); 
    for(int i = 0; i <grid_params::fNpT; i++){
      double pT = grid_params::fMref*tan(atan(grid_params::fPTMax/grid_params::fMref)*(i+0.5)/(grid_params::fNpT-1)); 
      fprintf(pFile,"%15.5e\t%15.5e", pT, ur);

      if (IsKeq) {
        double Keq1 = get_K_mu_qag(pT,  ur, 0) ;
        double Keq2 = get_K_mu_qag(pT,  ur, 1) ;
        fprintf(pFile,"\t%15.5e\t%15.5e", Keq1, Keq2);
      } 
      if (IsKshear) {
        double Kshear1 = get_K_mu_qag(pT,  ur, 2) ;
        double Kshear2 = get_K_mu_qag(pT,  ur, 3) ;
        double Kshear3 = get_K_mu_qag(pT,  ur, 4) ;
        double Kshear4 = get_K_mu_qag(pT,  ur, 5) ;
        fprintf(pFile,"\t%15.5e\t%15.5e\t%15.5e\t%15.5e", Kshear1,Kshear2,Kshear3,Kshear4);
      }
      if (IsKbulk) {
        double Kbulk1 = get_K_mu_qag(pT,  ur, 6) ;
        double Kbulk2 = get_K_mu_qag(pT,  ur, 7) ;
        fprintf(pFile,"\t%15.5e\t%15.5e", Kbulk1, Kbulk2);
      } 
      fprintf(pFile,"\n");
    }
    fprintf(pFile,"\n");
  }
  fclose(pFile);
}
