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


// List of parameters controlling the momentum discretization of irreducible components
// See also TParticle class definition
#include "grid_params.h"

//! TFastReso class is the top wrapper class of the fast freeze-out implementation.
class TFastReso {
  public:
    int fParticleNumber;
    //! reads the inputfile of resonances and creates TParticle entry in fParticleData vector.
    //! This routine should be edited by the user to match their input file format.
    //! Currently accepted format is (sstr is a string stream from the input file) is the one used
    //! in particle.data file in THERMINATOR 2 	[arXiv:1102.0273]
    //!  Name   Mass   Gamma   Sppin  Isospin  I3    Nq    Ns    Naq    Nas    Nc    Nac    PDGCode; 
    //!  Lines starting with '#' are ignored.
    virtual void read_particles_data(std::string inputname, std::string comps, bool verbose=false) = 0;

    //! The function call to initialized the particle distribution functions in fParticleData. Input parameters
    //! are the freeze-out temperature, chemical potential and (optionally), speed of sound (only needed for bulk perturbations).
    virtual void do_thermal(double Tfo, double MuB=0, double MuI3=0, double MuS=0, double MuC=0, double Cs2=0.14) = 0;

    //! Read in the input file of resonance decays and perfom the decays on the fly. It is important that decays 
    //! are mass ordered to included all feeddown contributions to lower mass resonances.
    //! The table format is the same as decays.data in THERMINATOR 2 	[arXiv:1102.0273]
    //! Father Child1 Child2 (optional Child3) BranchingRatio ClebschGordanCoeff
    //! Note that ClebschGordanCoeff term is 0 or 1 identifying if the branching ration needs
    //! to be reweighted by actual Clebsch Gordon coefficients.
    virtual void do_decays(std::string inputname, bool verbose=false) = 0;
  public: 
    //TFastReso(std::string comps);
    //! routine to get the pointer to TParticle class for a particle particle name"
    int getNumberOfParticles(){return fParticleNumber;};
    //virtual TParticle * getParticle(std::string name);
    //virtual TParticle * getParticle(int i);
};

// Speed of sound squared as a function of temperature T=x in GeV
inline double cs2FluiduM(double x) {
  return (x*x*(0.0010461910330715387 + x*(-0.016997351850612453 + x*(0.12595528994069613 + (-0.510477729039857 + x)*x)))*
     (3.831255349484765e-8 + x*(0.00002140058802480486 + x*(-0.00008693192041426982 + 
             x*(-0.018444802527704172 + x*(0.5368632648594895 + x*(-7.4522565817870365 + x*(62.30776783074574 + x*(-336.8744183413884 + x*(1191.026405027395 + x*(-2590.1225767072065 + 2712.1934046206975*x)))))))))
       ))/(7.313548422127797e-15 + x*(8.356306148115294e-12 + x*(2.3387986355395423e-9 + 
           x*(-6.605395132516945e-9 + x*(-1.790148526211301e-6 + x*(8.958041058260661e-6 + 
                    x*(0.0014381954623440515 + x*(-0.045458133957877116 + x*(0.726263758547687 + 
                             x*(-7.564793448616153 + x*(56.12131032495715 + x*(-308.4534327185232 + 
                                      x*(1271.3766801470479 + x*(-3884.3701789306033 + x*(8436.311867208919 + x*(-11805.768858200583 + 8136.580213862093*x))))))))))))))));
}


#endif
