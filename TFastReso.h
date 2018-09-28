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
#ifndef FASTRESO_TFastReso_h
#define FASTRESO_TFastReso_h

#include <fstream>
#include <memory>
#include <map>
#include <vector>
#include <string>
#include <gsl/gsl_integration.h>

// The class storing the properties of particles and their irreducible distrubtion
// function components f_i.  Printing out routines are also defined there.
#include "TParticle.h"

// List of parameters for GSL integration procedure.
#include "qag_params.h"

// List of parameters controlling the momentum discretization of irreducible components
// See also TParticle class definition
#include "grid_params.h"

//! TFastReso class is the top wrapper class of the fast freeze-out implementation.
class TFastReso {
  private:
    //! Data vector of TParticle class of all read-in particles
    //! See TParticle class definition.
    std::vector <std::unique_ptr<TParticle>> fParticleData;

    //! Map between the particle name and its position in fParticleData vector;
    std::map <std::string, size_t> fIndexTable;
    
    // Precalcullated Clebsch Gordan coefficients for reweighting branching ratios
    std::map <std::string, double> fClebschGordanTable;
  public:
    //! reads the inputfile of resonances and creates TParticle entry in fParticleData vector.
    //! This routine should be edited by the user to match their input file format.
    //! Currently accepted format is (sstr is a string stream from the input file) is the one used
    //! in particle.data file in THERMINATOR 2 	[arXiv:1102.0273]
    //!  Name   Mass   Gamma   Sppin  Isospin  I3    Nq    Ns    Naq    Nas    Nc    Nac    PDGCode; 
    //!  Lines starting with '#' are ignored.
    void read_particles_data(std::string inputname);

    //! The function call to initialized the particle distribution functions in fParticleData. Input parameters
    //! are the freeze-out temperature, chemical potential and (optionally), speed of sound (only needed for bulk perturbations).
    void do_thermal(double Tfo, double MuB=0, double MuI3=0, double MuS=0, double MuC=0, double Cs2=0.14);

    //! Read in the input file of resonance decays and perfom the decays on the fly. It is important that decays 
    //! are mass ordered to included all feeddown contributions to lower mass resonances.
    //! The table format is the same as decays.data in THERMINATOR 2 	[arXiv:1102.0273]
    //! Father Child1 Child2 (optional Child3) BranchingRatio ClebschGordanCoeff
    //! Note that ClebschGordanCoeff term is 0 or 1 identifying if the branching ration needs
    //! to be reweighted by actual Clebsch Gordon coefficients.
    void do_decays(std::string inputname);
  private:
    //! The sub-routines applying 2-body and 3-body transformation rules to the irreducible
    //! components. 
    void do_2bodydecay(std::string Father, std::string Child1, std::string Child2, double BranchingRatio);
    void do_3bodydecay(std::string Father, std::string Child1, std::string Child2, std::string Child3, double BranchingRatio);
  public: 
    //! routine to get the pointer to TParticle class for a particle particle name"
    TParticle * getParticle(std::string name){return fParticleData[fIndexTable[name]].get();};
    TParticle * getParticle(int i){return (i<0) ? fParticleData[fParticleData.size()+i].get()  : fParticleData[i].get();};
    
    //! The printout of the calculated components is done by printing procedures in TParticle class
};

#endif
