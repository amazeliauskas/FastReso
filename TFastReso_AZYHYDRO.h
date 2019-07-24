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
#ifndef FASTRESO_TFastReso_AZYHYDRO_h
#define FASTRESO_TFastReso_AZYHYDRO_h

#include "TFastReso.h"
#include <fstream>
#include <memory>
#include <map>
#include <vector>
#include <string>
#include <gsl/gsl_integration.h>

// The class storing the properties of particles and their irreducible distrubtion
// function components f_i.  Printing out routines are also defined there.
#include "TParticle_AZYHYDRO.h"

// List of parameters for GSL integration procedure.
#include "qag_params.h"

// List of parameters controlling the momentum discretization of irreducible components
// See also TParticle_AZYHYDRO class definition
#include "grid_params.h"

//! TFastReso_AZYHYDRO class is the top wrapper class of the fast freeze-out implementation.
class TFastReso_AZYHYDRO: public TFastReso {
  private:
    //! Data vector of TParticle_AZYHYDRO class of all read-in particles
    //! See TParticle_AZYHYDRO class definition.
    std::vector <std::unique_ptr<TParticle_AZYHYDRO>> fParticleData;

    //! Map between the particle name and its position in fParticleData vector;
    std::map <int, size_t> fIndexTable;
    
  public:
    //! reads the inputfile of resonances and creates TParticle_AZYHYDRO entry in fParticleData vector.
    //! This routine should be edited by the user to match their input file format.

    void read_particles_data(std::string inputname, bool verbose=false) override;

    //! The function call to initialized the particle distribution functions in fParticleData. Input parameters
    //! are the freeze-out temperature, chemical potential and (optionally), speed of sound (only needed for bulk perturbations).
    void do_thermal(double Tfo, double MuB=0, double MuI3=0, double MuS=0, double MuC=0, double Cs2=-1.0) override;

    //! Read in the input file of resonance decays and perfom the decays on the fly. It is important that decays 
    //! are mass ordered to included all feeddown contributions to lower mass resonances.
    //! The table format is pdg_weak.dat found in many derivatives of AZYHYDRO code, e.g. see MUSIC or VISHNU
    void do_decays(std::string inputname, bool verbose=false) override;
  public: 
    //! routine to get the pointer to TParticle_AZYHYDRO class for a particle particle name"
    TParticle_AZYHYDRO * getParticleByPDG(int pdgid){return fParticleData[fIndexTable[pdgid]].get();};
    TParticle_AZYHYDRO * getParticle(int i){return (i<0) ? fParticleData[fParticleData.size()+i].get()  : fParticleData[i].get();};
    int getParticleIndexByPDG(int pdgid){return fIndexTable[pdgid];};
    
    //! The printout of the calculated components is done by printing procedures in TParticle_AZYHYDRO class
};

#endif
