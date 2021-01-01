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
#ifndef FASTFO_TParticle_AZYHYDRO_h
#define FASTFO_TParticle_AZYHYDRO_h

#include "TParticle.h"
#include <string>
#include <math.h>

//! Particle class with basic properties and universal decay spectra components
class TParticle_AZYHYDRO : public TParticle {
  private:
    double fQBot;     // bottom charge
    double fCharge;     // bottom charge
    int fNdecays;        // number of decays
  public:
    //! Contructor with particle properties
    TParticle_AZYHYDRO(int pdgcode, std::string name, double mass, double gamma, int nu,
    double qb, double qs, double qc, double qbot, double isospin, double charge, int ndecays, std::string comps);
    TParticle_AZYHYDRO(std::string tag);
};

#endif
