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
#ifndef FASTFO_TKernel_h
#define FASTFO_TKernel_h

#include "TParticle.h"
#include "qag_params.h"
#include "grid_params.h"


//! Particle decay kernel integrate in eta and phi
class TKernel {
  private:
    TParticle * fParticle;
    //! Internal variables
    const double fEtaMax = 5.0;
    gsl_integration_workspace * fWphi;
    gsl_integration_workspace * fWeta;
    void initialize();
  public:
    TKernel(TParticle *particle);
    ~TKernel() ;
    double get_K_mu_detadphi_i(double eta, double phi,double pT, double ur,  int index) ;
    double get_K_mu_qag(double pT, double ur, int index);
    void print(std::string tag, std::string columns );
};


#endif
