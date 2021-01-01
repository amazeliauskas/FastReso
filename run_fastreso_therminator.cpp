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
#include <memory>
#include <string>
#include <TFastReso_THERMINATOR.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_errno.h>
using namespace std ;


int main(int argc, const char * argv[])
{

  if (argc<4) {
    cerr << argv[0] << " not enough input arguments " << endl;
    cerr << argv[0] << " particles.data decays.data ./outputfolder/" << endl;
    exit(EXIT_FAILURE) ;
  }
/////////////////// DANGER ///////////////////////////////////////////////////
  // If integration fails, you can turn the handler off. Integration functions
  // in TFastReso_THERMINATOR.cpp and TFastReso_formulas.h will try reducing accuracy
  // for the failing cases.  Use with caution!
  gsl_set_error_handler_off ();
/////////////////// DANGER ///////////////////////////////////////////////////


  // freeze-out temperature
  string pdata(argv[1]);// = "particles.data";
  string ddata(argv[2]);// = "decays_without_weak.data";
  string tag(argv[3]);

    bool verbose = true ; // print what is being done
    const int ns=6;
    // which particles to print out
    const int pdg[ns]={
        211,
        321,
        2212,
        3122,
        3312,
        3334
    };

    // freeze-out temperature
//    const int NT = 6; // number of temperatures
//    double Tmin=0.120; // min temp
//    double Tmax=0.170; // max temp
//    double dT=(Tmax-Tmin)/(NT-1); // temp step
//#pragma omp parallel for
//    for (int i =0; i <NT; i++){
    {
//        double Tfo=Tmin+dT*i; //GeV
        double Tfo=0.145; //GeV freeze-out temperature
        char buffer [50];
        sprintf (buffer, "%.4f", Tfo);
        string Ttag(buffer);
        string inarg= "Feq Fshear Fbulk Ftemp Fvel"; // which components to print out
        
        // create the TFastReso class
        TFastReso_THERMINATOR fastreso;
        // read all particles
        fastreso.read_particles_data(pdata,inarg, verbose);

        // initialize components with thermal distributions at freeze-out temperature
        cout << " Do thermal T = " <<Tfo  <<" GeV" <<endl;
        fastreso.do_thermal(Tfo);

        // print out particles
        for (int k =0; k<ns; k++){ 
        sprintf (buffer, "PDGid_%d", pdg[k]);
        string name(buffer);
          fastreso.getParticleByPDG(pdg[k])->print(tag+name+"_thermal_T"+Ttag); }


        // perform decays
        cout << " Do decays T = " <<Tfo  <<" GeV" <<endl;
        fastreso.do_decays(ddata, verbose);

        // print out components
        for (int k =0; k<ns; k++){ 
        sprintf (buffer, "PDGid_%d", pdg[k]);
        string name(buffer);
          fastreso.getParticleByPDG(pdg[k])->print(tag+name+"_total_T"+Ttag); }
    }
}
