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

#include <iostream>
#include <string>
#include <TKernel.h>
#include <TParticle_AZYHYDRO.h>
#include <gsl/gsl_errno.h>
using namespace std ;


int main(int argc, const char * argv[])
{
  if (argc<2) {
    cerr << argv[0] << " not enough input arguments " << endl;
    cerr << argv[0] << " ./inputfolder/ ./outputfolder/" << endl;
    exit(EXIT_FAILURE) ;
  }

  string tag1(argv[1]);
  string tag2(argv[2]);

    const int ns=6;

    const int pdg[ns]={
        211,
        321,
        2212,
        3122,
        3312,
        3334
    };

//3334	Omega	1.67245	0	4	1	-3	0	0	0	-1	1
//3122	Lambda	1.11568	0	2	1	-1	0	0	0	0	1
//3312	Ksi-	1.32131	0	2	1	-2	0	0	0.5	-1	1
//333	phi(1020)	1.019	0.00426	3	0	0	0	0	0	0	8
//331	eta'(958)	0.958	0.0002	1	0	0	0	0	0	0	4
//2112	n	0.9396	0	2	1	0	0	0	0.5	0	1
//2212	p	0.9383	0	2	1	0	0	0	0.5	1	1
//313	K*(892)0	0.896	0.0507	3	0	1	0	0	0.5	0	2
//323	K*(892)+	0.892	0.0508	3	0	1	0	0	0.5	1	2
//223	omega(782)	0.783	0.00844	3	0	0	0	0	0	0	3
//213	rho+	0.775	0.1492	3	0	0	0	0	1	1	1
//113	rho0	0.775	0.1492	3	0	0	0	0	1	0	1
//221	eta	0.547	0.00118	1	0	0	0	0	0	0	3
//311	K0	0.498	0	1	0	1	0	0	0.5	0	1
//321	K+	0.494	0	1	0	1	0	0	0.5	1	1
//211	pi+	0.14	0	1	0	0	0	0	1	1	1
//111	pi0	0.135	0	1	0	0	0	0	1	0	1
//22	Gamma	0	0	2	0	0	0	0	0	0	1
  /////////////////// DANGER ///////////////////////////////////////////////////
  // If integration fails, you can turn the handler off. Integration functions
  // in TFastReso_AZYHYDRO.cpp and TFastReso_formulas.h will try reducing accuracy
  // for the failing cases.  Use with caution!
  gsl_set_error_handler_off ();
  /////////////////// DANGER ///////////////////////////////////////////////////



  // freeze-out temperature

//  const int NT = 13;
//  double Tmin=0.120;
//  double Tmax=0.180;
//  double dT=(Tmax-Tmin)/(NT-1);
//
//#pragma omp parallel for
//  for (int i =0; i <NT; i++)
  {
//    double Tfo=Tmin+dT*i; //GeV
   double Tfo=0.145; //GeV
    char buffer [50];
    sprintf (buffer, "%.4f", Tfo);
    string Ttag(buffer);    
    cout << " Do kernels = " <<Ttag  <<" GeV" <<endl;
    for (int k =0; k<ns; k++){
        sprintf (buffer, "PDGid_%d", pdg[k]);
        string name(buffer);
      cout << " Do particle  " <<tag1+name+"_thermal_T"+Ttag  <<endl;
      TParticle_AZYHYDRO particle(tag1+name+"_thermal_T"+Ttag);
      TKernel kernel(&particle);
      kernel.print(tag2+name+"_thermal_T"+Ttag,"Keq Kshear Kbulk" );
    }
    for (int k =0; k<ns; k++){
        sprintf (buffer, "PDGid_%d", pdg[k]);
        string name(buffer);
      cout << " Do particle  " <<tag1+name+"_total_T"+Ttag  <<endl;
      TParticle_AZYHYDRO particle(tag1+name+"_total_T"+Ttag);
      TKernel kernel(&particle);
      kernel.print(tag2+name+"_total_T"+Ttag,"Keq Kshear Kbulk" );
    }

  }
}
