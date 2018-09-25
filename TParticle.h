/*
 * Copyright (c) 2018 Aleksas Mazeliauskas
 * All rights reserved.
 *
 * FastReso is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/amazeliauskas/FastReso/
 */
#ifndef FASTFO_TParticle_h
#define FASTFO_TParticle_h

#include <string>
#include <math.h>
#include <grid_params.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
//! Enum class for deciding on Bose-Einstein, Fermi-Dirac or Botlzmann particle dsitribution
enum class EParticleType : int { kBoson = 1, kFermion = -1, kBoltzman = 0 };

//! Particle class with basic properties and universal decay spectra components
class TParticle {
  private:
    //! Particle properties
    std::string fParticleName;  // particle name
    double fMass;   // [GeV] resonance mass
    double fGamma;  // [GeV] resonance width (currently not used)
    double fSpin;   // spin
    double fIsospin; // isospin J |J,m>
    double fI3;     // isospin projection m: |J,m>
    double fNq;     // number of up or down quarks
    double fNs;     // number of strange quarks
    double fNaq;    // number of up or down antiquarks
    double fNas;    // number of strange antiquarks
    double fNc;     // number of charm quarks
    double fNac;    // number of charm antiquarks
    int fPDGCode;   // number code according to PDG
    int fNu;        // degeneracy of a particle
    double fQB;     // baryon charge
    double fQS;     // strange charge
    double fQC;     // charme charge
    EParticleType fParticleType; // particle statistic (bose-einstein, fermi-dirac, boltzmann)

    double fTfo;
    double fNtotal; //integrated particle number
    double fNtotal_buffer;

    //! Grid and interpolator for decay spectra components
    const gsl_interp_type   *fSpline_type = gsl_interp_cspline;
    double fPbar_arr[grid_params::fNpbar]; // array of fluid restframe momentum
    double fFj_arr[grid_params::fNf][grid_params::fNpbar]; // array of scalar functions Fj
    double fFj_arr_buffer[grid_params::fNf][grid_params::fNpbar]; // buffer array of scalar functions Fj 
    gsl_spline *fSpline_Fj ;
    gsl_interp_accel *fPbar_acc;
    //! Bool variables to control data read in/out-put
    bool fIsLocked = false; //if true, prevent modifying fFj_arr data
    bool fIsModified[grid_params::fNf]; //if true, require initilization of interpolator
    bool fIsLoaded[grid_params::fNf] ; //if false, load new data to interpolator
    //! After any update of the grid, don't forget to reinitialized the interpolator!
    //! initiliazes the interpolator
    void init(int j);
  public:
    //! Contructor with particle properties
    TParticle(std::string name, double mass, double gamma, double spin,
    double isospin, double i3, double nq, double ns, double naq, double nas, double nc, double nac,
    int pdgcode);
    ~TParticle() ;
    //! return parameters about the particle
    std::string getName() {return fParticleName;};
    double getM() {return fMass;};
    double getGamma() {return fGamma;};
    double getQB() {return fQB;};
    double getQS() {return fQS;};
    double getQC() {return fQC;};
    double getI3() {return fI3;};
    double getIsospin() {return fIsospin;};
    double getSpin() {return fSpin;};
    double getN() {return fNtotal;};
    int getNu() {return fNu;};
    EParticleType getType() {return fParticleType;};

    //! increment integrated particle number
    void addN(double n) {
            if (fIsLocked){
        std::cerr << "\033[1mTParticle.h\033[0m : \033[1;31merror\033[0m : attemped to modify locked "<< fParticleName <<" data" << std::endl;
        exit(EXIT_FAILURE);}
      else { fNtotal+=n; fNtotal_buffer+=n; }};
    //! set initial thermal temperature
    void setTfo(double Tfo) {fTfo=Tfo;};
    //! set lock on data modification
    void lock() {fIsLocked=true;};
    //! clean the buffer
    void clean_buffer() { 
      for(int i = 0; i <grid_params::fNpbar; i++){
        for (int jj=0; jj<grid_params::fNf; jj++){
          fFj_arr_buffer[jj][i]=0.0;
        }   
        fNtotal_buffer=0.0; }
    };

   //! return the energy and momentum of grid point i
   double getEbar(int i) {return sqrt(fMass*fMass+getPbar(i)*getPbar(i));};
   double getPbar(int i) {return fPbar_arr[i];};
   //! Return the number of momentum grid points
    int getNpbar() {return grid_params::fNpbar;};
    //! Return of the initialization/freeze-out temperature
    double getTfo() {return fTfo;};

    //! increment Fj array at cite i by F
    void addFj(int j, int i, double F) { 
      if (fIsLocked){
        std::cerr << "\033[1mTParticle.h\033[0m : \033[1;31merror\033[0m : attemped to modify locked "<< fParticleName <<" data" << std::endl;
        exit(EXIT_FAILURE);}
      else {
        fFj_arr[j][i]+=F;
        fFj_arr_buffer[j][i]+=F;
        fIsModified[j]=true;
      }
    };

    //! Returns the universal scalar functions Fj at Ebar
    double get_Fj(int j, double Ebar);

    //Routines to print the data to screen or file
    void print();
    void print(std::string tag);
    void print_buffer();
    void print_buffer(std::string tag);
};

#endif
