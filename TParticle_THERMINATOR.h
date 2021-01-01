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
#ifndef FASTFO_TParticle_THERMINATOR_h
#define FASTFO_TParticle_THERMINATOR_h

#include "TParticle.h"
#include <string>
#include <math.h>
#include <grid_params.h>

//! Particle class with basic properties and universal decay spectra components
class TParticle_THERMINATOR: public TParticle {
  private:
    //! Particle properties
    double fSpin;   // spin
    double fI3;     // isospin projection m: |J,m>
    double fNq;     // number of up or down quarks
    double fNs;     // number of strange quarks
    double fNaq;    // number of up or down antiquarks
    double fNas;    // number of strange antiquarks
    double fNc;     // number of charm quarks
    double fNac;    // number of charm antiquarks
  public:
    //! Contructor with particle properties
    TParticle_THERMINATOR(std::string name, double mass, double gamma, double spin,
    double isospin, double i3, double nq, double ns, double naq, double nas, double nc, double nac,
    int pdgcode, std::string comps);

    TParticle_THERMINATOR(std::string nametag);

    //! return parameters about the particle
    double getI3() {return fI3;};
    double getSpin() {return fSpin;};

};

#endif
