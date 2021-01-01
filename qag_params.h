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
#ifndef FASTFO_qag_params_h
#define FASTFO_qag_params_h
#include <gsl/gsl_integration.h>
namespace qag_params{
  const double fEpsRel = 1e-6;
  const double fEpsAbs = 1e-10;
  const int fKey = GSL_INTEG_GAUSS51;
  const int fLimit = 1000;
}
#endif
