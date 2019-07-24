#ifndef FASTRESO_TFastReso_formulas_h
#define FASTRESO_TFastReso_formulas_h
#include <gsl/gsl_integration.h>
#include "TParticle.h"
// List of parameters for GSL integration procedure.
#include "qag_params.h"
#include "grid_params.h"

//using namespace std;
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
  double pbar = sqrt(Ebar*Ebar-m*m);
  switch(j){
    case 0: //feq 1
      return pbar*get_thermal_F(Ebar,T, QMu, type);
      break;
    case 1: //feq 2
      return pbar*get_thermal_F(Ebar,T, QMu, type);
      break;
    case 2: //fshear 1
      return pbar*pbar*pbar*get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 3: //fshear 2
      return pbar*pbar*pbar*get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 4: //fshear 3
      return pbar*pbar*pbar*get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 5: //fbulk 1
      if (Cs2<0 or Cs2>1.0){
      std::cerr  << "\033[1mTFastReso_formulas.h::get_initial_Fj\033[0m : \033[1;31merror\033[0m :  please provide physical speed of sound for bulk initialization" << std::endl;
      exit(EXIT_FAILURE);}
      else{
      return pbar*get_thermal_dF(Ebar,T, QMu, type)*(1./3.*m*m/T/Ebar-Ebar/T*(1./3.-Cs2));
      }
      break;
    case 6: //fbulk 2
      if (Cs2<0 or Cs2>1.0){
      std::cerr  << "\033[1mTFastReso_formulas.h::get_initial_Fj\033[0m : \033[1;31merror\033[0m :  please provide physical speed of sound for bulk initialization" << std::endl;
      exit(EXIT_FAILURE);}
      else{
      return pbar*get_thermal_dF(Ebar,T, QMu, type)*(1./3.*m*m/T/Ebar-Ebar/T*(1./3.-Cs2));
      }
      break;
    case 7: //ftemperature. 1
      return pbar*get_thermal_dF(Ebar,T, QMu, type)*Ebar/T;
      break;
    case 8: //ftemperature 2
      return pbar*get_thermal_dF(Ebar,T, QMu, type)*Ebar/T;
      break;
    case 9: //fvelocity  1
      return pbar*pbar*get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 10: //fvelocity 2
      return pbar*pbar*get_thermal_dF(Ebar,T, QMu, type);
      break;
    case 11: //fvelocity 3
      return pbar*pbar*get_thermal_dF(Ebar,T, QMu, type);
      break;
    default:
      std::cerr  << "\033[1mTFastReso_formulas.h::get_initial_Fj\033[0m : \033[1;31merror\033[0m :  wrong index "  << j << std::endl;
      exit(EXIT_FAILURE);
  }
}
//! Transformation rule factor A for Fj in a 2-body decay
//! Here you can add additional rules for fj's
double get_factor_Aj(int j, double Qw, double Ew, double Ebar, double pbar, double Ma, double pw)
{

//  double A= (Ma*Eabc/Mb/Mb-w*Ma*pabc/Mb/Mb*Ebar/pbar);
//  double Ew= (Ma*Eabc*Ebar/Mb/Mb-w*Ma*pabc*pbar/Mb/Mb);
 
  switch(j){
    case 0: // feq 1
      return Qw/pw;
      break;
    case 1: // feq 2
      return pbar/pw*Ew/Ebar;
      break;
    case 2: // fshear 1
      return Qw/pw*(2.5*Qw/pw*Qw/pw-1.5);
      break;
    case 3: // fshear 2
      return Qw/pw;
      break;
    case 4: // fshear 3
      return (1.5*Qw/pw*Qw/pw-0.5)*Ew/Ebar*pbar/pw;
      break;
    case 5: // fbulk 1
      return Qw/pw;
      break;
    case 6: // fbulk 2
      return pbar/pw*Ew/Ebar;
      break;
    case 7: // ftemperature 1
      return Qw/pw;
      break;
    case 8: // ftemperature 2
      return pbar/pw*Ew/Ebar;
      break;
    case 9: // fvelocity 1
      return 1.5*Qw/pw*Qw/pw-0.5;
      break;
    case 10: // fvelocity 2
      return 1;
      break;
    case 11: // fvelocity 3
      return pbar/pw*Qw/pw*Ew/Ebar;
      break;
    default:
      std::cerr  << "\033[1mTFastReso_formulas.h::get_factor_Aj\033[0m : \033[1;31merror\033[0m :  wrong index "  << j << std::endl;
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
  //double Ebar_old = pbar*pabc*Ma/Mb/Mb*(1-w);
  //double Ebar_old = pbar*pabc/Ma*u;
  double F_old=Father->get_Fj(index, Ebar_old);
  double Qw= (Ma*Eabc*pbar/Mb/Mb-w*Ma*pabc/Mb/Mb*Ebar);
  //double Qw = pabc/Ma*u;
//  double Ew= (Ma*Eabc*Ebar/Mb/Mb-w*Ma*pabc*pbar/Mb/Mb);
  double A = get_factor_Aj(index, Qw, Ebar_old, Ebar, pbar, Ma,sqrt(Qw*Qw+(1-w*w)*Ma*Ma/Mb/Mb*pabc*pabc));
  return F_old*A;
}
double get_F_pu_du(double u, void * p) {
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
  //double Ebar_old = Ebar*Eabc*Ma/Mb/Mb-w*pbar*pabc*Ma/Mb/Mb;
  //double Ebar_old = pbar*pabc*Ma/Mb/Mb*(1-w);
  if (pbar==0) { return 0;};
  double Ebar_old = Ma/2*(pbar/pabc+pabc/pbar)+u*Ma;
  double F_old=Father->get_Fj(index, Ebar_old);
  //double Qw= (Ma*Eabc/Mb/Mb-w*Ma*pabc/Mb/Mb*Ebar/pbar);
  double Qw = Ma/2/pbar*(pbar*pbar/pabc-pabc)+Ma*u;
//  double Ew= (Ma*Eabc*Ebar/Mb/Mb-w*Ma*pabc*pbar/Mb/Mb);
  double A = get_factor_Aj(index, Qw, Ebar_old, Ebar, pbar, Ma,sqrt(Qw*Qw+2*u*pabc/pbar));
  return F_old*A/(pbar*pabc/Ma/Ma);
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
  int status = gsl_integration_qag (twobody_function, -1, 1, fEpsAbs,fEpsRel, fLimit,fKey ,twobody_workspace, &result, &error); 
  // Check status of integration and reduce accuracy if failing
  if (status) { 
    status = gsl_integration_qag (twobody_function, -1, 1, fEpsAbs,1e-4, fLimit,GSL_INTEG_GAUSS15,twobody_workspace, &result, &error); 
    if (status) { 
      status = gsl_integration_qag (twobody_function, -1, 1, 1e-4,1e-4, fLimit,GSL_INTEG_GAUSS15,twobody_workspace, &result, &error); 
      if (status) {      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }
  return result*get_pabc(Ma,Mb,mct)*get_pabc(mct, Mc,Md);
}
double get_F_pu_dm_u(double mct, void * p) {
  gsl_function *twobody_function = ((struct threebody_params *) p)->twobody_function;
  ((struct twobody_params *) twobody_function->params)->Mc = mct;
  gsl_integration_workspace *twobody_workspace = ((struct threebody_params *) p)->twobody_workspace;
  double &Ma = ((struct twobody_params *) twobody_function->params)->Ma;
  double &Mb = ((struct twobody_params *) twobody_function->params)->Mb;
  double &Mc             = ((struct threebody_params *) p)->Mc;
  double &Md             = ((struct threebody_params *) p)->Md;
  double result,error;

  using namespace qag_params;
  int status = gsl_integration_qagiu (twobody_function, 0.0, fEpsAbs,fEpsRel, fLimit,twobody_workspace, &result, &error); 
  // Check status of integration and reduce accuracy if failing
  if (status) { 
  status = gsl_integration_qagiu (twobody_function, 0.0, fEpsAbs,1e-4, fLimit,twobody_workspace, &result, &error); 
    if (status) { 
  status = gsl_integration_qagiu (twobody_function, 0.0, 1e-4,1e-4, fLimit,twobody_workspace, &result, &error); 
      if (status) {      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }

  return result*get_pabc(Ma,Mb,mct)*get_pabc(mct, Mc,Md);
}

void do_2bodydecay(TParticle * Father, TParticle * Child1, TParticle * Child2, double BranchingRatio);
void do_3bodydecay(TParticle * Father, TParticle * Child1, TParticle * Child2, TParticle * Child3, double BranchingRatio);

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
      std::cerr << "\033[1mTFastReso_formulas.h::get_thermal_F\033[0m : \033[1;31merror\033[0m : unrecognined particle type!" <<std::endl;
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
      std::cerr << "\033[1mTFastReso_formulas.h::get_thermal_F_pu\033[0m : \033[1;31merror\033[0m : unrecognined particle type!" <<std::endl;
      exit(EXIT_FAILURE);}
}


//! Perform two body decay integral: Father->Child1 + Child2
void do_2bodydecay(TParticle * Father, TParticle * Child1, TParticle * Child2, double BranchingRatio){
  using namespace qag_params;

  gsl_integration_workspace * fWorkspace2body ;
  fWorkspace2body = gsl_integration_workspace_alloc (fLimit);
  Father->lock();
  double Ma = Father->getM();
  double Mb = Child1->getM();
  double Mc = Child2->getM();
  double nua = Father->getNu();
  double nub = Child1->getNu();
  double Ebar=0.0;
  Child1->clean_buffer();
  if (Ma < Mb+Mc) {
    std::cerr  << "\033[1mTFastReso_formulas.h::do_2bodydecays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc << std::endl;
    exit(EXIT_FAILURE);
  }
  gsl_function func2body;
      if (Mb==0) {
  func2body.function = get_F_pu_du;
      } else {
  func2body.function = get_F_pu_dw;
      }
  twobody_params params = {Ebar, Ma, Mb, Mc, Father, 0};
  func2body.params = &params;
  double error, result ;
  int sym=1 ;
  if (Child1==Child2 ){ sym=2; } else { sym=1; }

  double Nfather = Father->getN();
  //! add fraction of the total yield to child
  Child1->addN(BranchingRatio*Nfather*sym);
  for (int i=0; i < Child1->getNpbar(); i++){
    double Ebar =Child1->getEbar(i);
    params.Ebar =Ebar;

    for (int j=0; j < grid_params::fNf ; j++){
      params.index =j;
      double fac = 0;
      if (Mb==0) {
int status = gsl_integration_qagiu (&func2body, 0.0, fEpsAbs,fEpsRel, fLimit ,fWorkspace2body, &result, &error); 

  // Check status of integration and reduce accuracy if failing
  if (status) { 
status = gsl_integration_qagiu (&func2body, 0.0, fEpsAbs,1e-4, fLimit ,fWorkspace2body, &result, &error); 
    if (status) { 
status = gsl_integration_qagiu (&func2body, 0.0, 1e-4,1e-4, fLimit ,fWorkspace2body, &result, &error); 
      if (status) {      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }





      fac = sym*BranchingRatio*nua/nub/2;
      }
      else {
  double Eabc = get_Eabc(Ma,Mb,Mc);
  double pabc = get_p(Eabc, Mb);  
  //double pbar = GSL_MAX(get_p(Ebar, Mb),DBL_EPSILON);  
  double pbar = get_p(Ebar, Mb);// ,DBL_EPSILON);  
int status = gsl_integration_qag (&func2body, -1, 1, fEpsAbs,fEpsRel, fLimit,fKey ,fWorkspace2body, &result, &error); 

  // Check status of integration and reduce accuracy if failing
  if (status) { 
status = gsl_integration_qag (&func2body, -1, 1, fEpsAbs,1e-4, fLimit, GSL_INTEG_GAUSS15,fWorkspace2body, &result, &error); 
    if (status) { 
status = gsl_integration_qag (&func2body, -1, 1, 1e-4,1e-4, fLimit, GSL_INTEG_GAUSS15 ,fWorkspace2body, &result, &error); 
      if (status) {      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }








      fac = sym*BranchingRatio*nua/nub*pow(Ma/Mb,2)/2;

      }
      Child1->addFj(j,i, fac*result);

    }
  }
  //if (Child1=="pi0139plu"  ) { getParticle(Child1)->print_buffer(Child1+"_"+Father+"_"+Child2);}
  gsl_integration_workspace_free (fWorkspace2body);

}

//! Perform three body decay integral Father->Child1 + Child2 + Child3
void do_3bodydecay(TParticle * Father, TParticle * Child1, TParticle * Child2, TParticle * Child3, double BranchingRatio){
  using namespace qag_params;
  gsl_integration_workspace * fWorkspace2body ;
  gsl_integration_workspace * fWorkspace3body ;
  fWorkspace2body = gsl_integration_workspace_alloc (fLimit);
  fWorkspace3body = gsl_integration_workspace_alloc (fLimit);
  Father->lock();
  double Ma = Father->getM();
  double Mb = Child1->getM();
  double Mc = Child2->getM();
  double Md = Child3->getM();
  double nua = Father->getNu();
  double nub = Child1->getNu();
  double Ebar=0.0;

  Child1->clean_buffer();
  if (Ma < Mb+Mc+Md) {
    std::cerr  << "\033[1mTFastReso_formulas.h::do_3bdoydecays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc+Md << std::endl;
    exit(EXIT_FAILURE);
  }
  gsl_function func2body;
  gsl_function func3body;
      if (Mb==0) {
  func2body.function = get_F_pu_du;
      } else {
  func2body.function = get_F_pu_dw;
      }
  func3body.function = get_B_dm;
  twobody_params params2 = {Ebar, Ma, Mb, Mc+Md, Father, 0};
  func2body.params = &params2;
  threebody_params params3 = {&func2body, Mc,Md, fWorkspace2body};
  func3body.params = &params3;
  double error, result ;
  int status =gsl_integration_qag (&func3body, Mc+Md, Ma-Mb, fEpsAbs,fEpsRel, fLimit,fKey ,fWorkspace3body, &result, &error); 
  // Check status of integration and reduce accuracy if failing
  if (status) { 
  status =gsl_integration_qag (&func3body, Mc+Md, Ma-Mb, fEpsAbs,1e-4, fLimit,GSL_INTEG_GAUSS15 ,fWorkspace3body, &result, &error); 
    if (status) { 
  status =gsl_integration_qag (&func3body, Mc+Md, Ma-Mb, 1e-4,1e-4, fLimit,GSL_INTEG_GAUSS15 ,fWorkspace3body, &result, &error); 
      if (status) {      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }


  double B3body = result;
      if (Mb==0) {
  func3body.function = get_F_pu_dm_u;
      } else {
  func3body.function = get_F_pu_dm;
      }

  int sym = 1 ;
  if (Child1==Child2 and Child2==Child3 ){ sym=3; } else if (Child1==Child2 or Child1==Child3 ){ sym=2;} else { sym=1; }
  double Nfather = Father->getN();
  // add fraction of yield to child
  Child1->addN(BranchingRatio*Nfather*sym);
  for (int i=0; i < Child1->getNpbar(); i++){
    double Ebar =Child1->getEbar(i);
    params2.Ebar =Ebar;
    for (int j=0; j < grid_params::fNf ; j++){
      params2.index =j;
      int status = gsl_integration_qag (&func3body, Mc+Md, Ma-Mb, fEpsAbs,fEpsRel, fLimit,fKey ,fWorkspace3body, &result, &error); 
  // Check status of integration and reduce accuracy if failing
  if (status) { 
      status = gsl_integration_qag (&func3body, Mc+Md, Ma-Mb, fEpsAbs,1e-4, fLimit,GSL_INTEG_GAUSS15 ,fWorkspace3body, &result, &error); 
    if (status) { 
      status = gsl_integration_qag (&func3body, Mc+Md, Ma-Mb, 1e-4,1e-4, fLimit,GSL_INTEG_GAUSS15 ,fWorkspace3body, &result, &error); 
      if (status) {      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE); 
      }else {
        std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTFastReso_formulas.h\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }


      double fac =0;
      if (Mb==0) {
        fac = sym*BranchingRatio*nua/nub/2/B3body;
      }
      else {
        fac = sym*BranchingRatio*nua/nub*pow(Ma/Mb,2)/2/B3body;
      }
      Child1->addFj(j,i, fac*result);

    }
  }


  //if (Child1=="pi0139plu" ) { getParticle(Child1)->print_buffer(Child1+"_"+Father+"_"+Child2+"_"+Child3);}
  gsl_integration_workspace_free (fWorkspace2body);
  gsl_integration_workspace_free (fWorkspace3body);
}


#endif
