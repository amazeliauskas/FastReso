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
#ifndef FASTRESO_TFastReso_THERMINATOR_h
#define FASTRESO_TFastReso_THERMINATOR_h

#include "TFastReso.h"
#include <fstream>
#include <memory>
#include <map>
#include <vector>
#include <string>
#include <gsl/gsl_integration.h>

// The class storing the properties of particles and their irreducible distrubtion
// function components f_i.  Printing out routines are also defined there.
#include "TParticle_THERMINATOR.h"


// List of parameters controlling the momentum discretization of irreducible components
// See also TParticle_THERMINATOR class definition
#include "grid_params.h"

//! TFastReso_THERMINATOR class is the top wrapper class of the fast freeze-out implementation.
class TFastReso_THERMINATOR: public TFastReso {
  private:
    //! Data vector of TParticle_THERMINATOR class of all read-in particles
    //! See TParticle_THERMINATOR class definition.
    std::vector <std::unique_ptr<TParticle_THERMINATOR>> fParticleData;

    //! Map between the particle name and its position in fParticleData vector;
    std::map <std::string, int> fIndexTable;
    std::map <int, int> fIndexTablePDG;
    
  public:

    //TFastReso_THERMINATOR(std::string comps) : TFastReso(comps){};
    //! reads the inputfile of resonances and creates TParticle_THERMINATOR entry in fParticleData vector.
    //! This routine should be edited by the user to match their input file format.
    //! Currently accepted format is (sstr is a string stream from the input file) is the one used
    //! in particle.data file in THERMINATOR 2 	[arXiv:1102.0273]
    //!  Name   Mass   Gamma   Sppin  Isospin  I3    Nq    Ns    Naq    Nas    Nc    Nac    PDGCode; 
    //!  Lines starting with '#' are ignored.
    void read_particles_data(std::string inputname, std::string comps="Feq Fshear Fbulk Ftemp Fvel", bool verbose=false) override;

    //! The function call to initialized the particle distribution functions in fParticleData. Input parameters
    //! are the freeze-out temperature, chemical potential and (optionally), speed of sound (only needed for bulk perturbations).
    void do_thermal(double Tfo, double MuB=0, double MuI3=0, double MuS=0, double MuC=0, double Cs2=0.14) override;
    void do_thermal(int pi, double Tfo, double QMu=0, double Cs2 = 0.14 );

    //! Read in the input file of resonance decays and perfom the decays on the fly. It is important that decays 
    //! are mass ordered to included all feeddown contributions to lower mass resonances.
    //! The table format is the same as decays.data in THERMINATOR 2 	[arXiv:1102.0273]
    //! Father Child1 Child2 (optional Child3) BranchingRatio ClebschGordanCoeff
    //! Note that ClebschGordanCoeff term is 0 or 1 identifying if the branching ration needs
    //! to be reweighted by actual Clebsch Gordon coefficients.
    void do_decays(std::string inputname, bool verbose=false) override;
  private:
    //! The sub-routines applying 2-body and 3-body transformation rules to the irreducible
    //! components. 
  //  void do_2bodydecay(TParticle_THERMINATOR * Father, TParticle_THERMINATOR  * Child1, TParticle_THERMINATOR * Child2, double BranchingRatio);
//    void do_3bodydecay(TParticle_THERMINATOR * Father, TParticle_THERMINATOR * Child1, TParticle_THERMINATOR * Child2, TParticle_THERMINATOR * Child3, double BranchingRatio);
  public: 
    //TFastReso_THERMINATOR() : TFastReso() {} ;
//#~TFastReso_THERMINATOR(): ~TFastReso() {};
    //! routine to get the pointer to TParticle_THERMINATOR class for a particle particle name"
    TParticle_THERMINATOR * getParticleByPDG(int pdgid){return fParticleData[fIndexTablePDG[pdgid]].get();};
    TParticle_THERMINATOR * getParticle(std::string name){return fParticleData[fIndexTable[name]].get();};
    TParticle_THERMINATOR * getParticle(int i){return (i<0) ? fParticleData[fParticleData.size()+i].get()  : fParticleData[i].get();};
    int  getParticleIndex(std::string name){return fIndexTable[name];};
    //! The printout of the calculated components is done by printing procedures in TParticle_THERMINATOR class
};


#endif
