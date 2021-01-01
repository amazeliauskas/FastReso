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
#ifndef FASTFO_grid_params_h
#define FASTFO_grid_params_h
// number of points in interpolation grid of rest-frame momentum fPbar_arr[i]
// As defined in TParticle class, the momentum is discretized on  bin centered grid
// pbar = M*tan(xmax*(i+0.5)/N), where M is some momentum scale (chosen to be 1GeV).
// Explicitly:
////    for (int i = 0; i < fNpbar; i++) {
////      fPbar_arr[i]=dp*(i+0.5); 
////      fPbar_arr[i]=grid_params::fMref*tan(atan(grid_params::fPbarMax/grid_params::fMref)*(i+0.5)/grid_params::fNpbar); 
////      }

namespace grid_params{
 static const int fNpbar = 401; // number of points fluid restfram pbar grid
 static const int fNf = 12; //  the maximum number of f_i functions to compute
 static const int fNpT = 41; // number of points fluid restfram pbar grid
 static const int fNur = 41; // number of points fluid restfram pbar grid
 static const double fPbarMax = 8.5; //  the max pbar value
 static const double fPTMax = 3.5; //  the max pbar value
 static const double fUrMax = 3.5; //  the max Ur value
 static const double fMref = 1.0; //  [GeV] momentum scale for tan discretization.
}

#endif
