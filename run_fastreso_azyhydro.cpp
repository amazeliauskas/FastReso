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

#include <iostream>
#include <memory>
#include <string>
#include <TFastReso_AZYHYDRO.h>

using namespace std ;


int main(int argc, const char * argv[])
{

    if (argc<4) {
        cerr << argv[0] << " not enough input arguments " << endl;
        cerr << argv[0] << " particles.data decays.data ./outputfolder/" << endl;
        exit(EXIT_FAILURE) ;
    }

    // freeze-out temperature
    string pdata(argv[1]);// = "pdg_weak.dat";
    string ddata(argv[2]);// = "reverse_pdg_weak.dat";
    string tag(argv[3]);
    const int ns=22;
    const string names[ns]={
        "gam000zer",
        "pi0135zer",
        "pi0139plu",
        "pi0139min",
        "Ka0492plu",
        "Ka0492min",
        "Ka0492zer",
        "Ka0492zrb",
        "pr0938plu", 
        "pr0938plb", 
        "Lm1115zer",
        "Lm1115zrb",
        "Sg1189plu",
        "Sg1189plb",
        "Sg1189min",
        "Sg1189mnb",
        "Xi1321zer",
        "Xi1321zrb",
        "Xi1321min",
        "Xi1321mnb",
        "UM1672min",
        "UM1672mnb"
    };
    const int pdg[ns]={
        22,
        111,
        211,
        -211,
        321,
        -321,
        311,
        -311,
        2212,
        -2212,
        3122,
        -3122,
        3222,
        -3222,
        3112,
        -3112,
        3322,
        -3322,
        3312,
        -3312,
        3334,
        -3334
    };
    // freeze-out temperature
    const int NT = 21;
    double Tmin=0.140;
    double Tmax=0.160;
    double dT=(Tmax-Tmin)/(NT-1);


#pragma omp parallel for
    for (int i =0; i <NT; i++){
        double Tfo=Tmin+dT*i; //GeV
        char buffer [50];
        sprintf (buffer, "%.4f", Tfo);
        string Ttag(buffer);
        TFastReso_AZYHYDRO fastreso;
        // read all particles
        fastreso.read_particles_data(pdata);
        cout << " Do thermal T = " <<Tfo  <<" GeV" <<endl;
        fastreso.do_thermal(Tfo);

        for (int k =0; k<ns; k++){ fastreso.getParticleByPDG(pdg[k])->print(tag+names[k]+"_thermal_T"+Ttag); }
        cout << " Do decays T = " <<Tfo  <<" GeV" <<endl;
        fastreso.do_decays(ddata);
        cout << " Finished decays T = " <<Tfo  <<" GeV" <<endl;

        for (int k =0; k<ns; k++){ fastreso.getParticleByPDG(pdg[k])->print(tag+names[k]+"_total_T"+Ttag); }
    }

}
