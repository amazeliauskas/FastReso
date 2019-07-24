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
  protected:
    //! Particle properties
    std::string fParticleName;  // particle name
    std::string fDescription;  // particle name
    std::string fDescriptionHeader;  // particle name
    std::string fHeader = "1:pbar [GeV]\t2:m [GeV]\t3:pbar*feq_1\t4:pbar*feq_2\t5:pbar^3*fshear_1\t6:pbar^3*fshear_2\t7:pbar^3*fshear_3\t8:pbar*fbulk_1\t9:pbar*fbulk_2\t10:pbar*ftemperature_1\t11:pbar*ftemperature_2\t12:pbar^2*fvelocity 1\t13:pbar^2*fvelocity_2\t14:pbar^2*fvelocity_3";
      
        // particle name
    double fMass;   // [GeV] resonance mass
    double fGamma;  // [GeV] resonance width (currently not used)
    double fIsospin; // isospin J |J,m>
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
    gsl_spline *fSpline_Fj[grid_params::fNf] ;
    gsl_interp_accel *fPbar_acc;
    //! Bool variables to control data read in/out-put
    bool fIsLocked = false; //if true, prevent modifying fFj_arr data
    bool fIsModified[grid_params::fNf]; //if true, require initilization of interpolator
    //! After any update of the grid, don't forget to reinitialized the interpolator!
    //! initiliazes the interpolator
    void init_grid();
    void init(int j);
  public:
    //! Contructor with particle properties
    ~TParticle() ;
    //! return parameters about the particle
    std::string getName() {return fParticleName;};
    double getM() {return fMass;};
    double getGamma() {return fGamma;};
    double getQB() {return fQB;};
    double getQS() {return fQS;};
    double getQC() {return fQC;};
    double getIsospin() {return fIsospin;};
    double getN() {return fNtotal;};
    int getNu() {return fNu;};
    int getPDG() {return fPDGCode;};
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
    bool hasdecayed() {return fIsLocked;};
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
    std::string getDescriptionHeader(){ return fDescriptionHeader;};
    std::string getDescription(){return fDescription;};
    void print();
    void print(std::string tag);
    void print_buffer();
    void print_buffer(std::string tag);
};
#endif
