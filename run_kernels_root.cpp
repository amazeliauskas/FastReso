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
#include "TTree.h"
#include "TThread.h"
#include "TDirectory.h"
#include "TFile.h"
#include "omp.h"
#include "TH3D.h"
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




    //ROOT::EnableThreadSafety();
    TThread::Initialize();
    // freeze-out temperature
    const int NT = 1;
    double Tmin=0.145;
    double Tmax=0.145;
    double betamin = 0.0;
    double betamax = 1.0;
    int Nbeta=100;
    double dbeta=(betamax-betamin)/(Nbeta);
    //double dT=(Tmax-Tmin)/(NT-1);
    double dT = 0.010; // replace for NT>1
    double pTmin = 0.0175;
    double pTmax = 3.5175;
    int NpT =101;
    double dpT=(pTmax-pTmin)/(NpT-1);
    const int Nk=8;
    const string ktag[Nk]={
        "Keq1",
        "Keq2",
        "Kshear1",
        "Kshear2",
        "Kshear3",
        "Kshear4",
        "Kbulk1",
        "Kbulk2"
    };


#pragma omp parallel for
    for (int k =0; k <NT; k++){
        char buffer [150];
        double Tfo=Tmin+dT*k; //GeV
        sprintf (buffer, "%.4f", Tfo);
        string Ttag(buffer);    
        cout << " Do kernels = " <<Ttag  <<" GeV" <<endl;
        TFile * f = TFile::Open((tag2+"K_total_full_T"+Ttag+".root").c_str(),"RECREATE");
        TH3D *hist[ns][Nk];
        for (int ki =0; ki<Nk; ki++){ 
            for (int si =0; si<ns; si++){ 
        sprintf (buffer, "PDGid_%d", pdg[si]);
        string name(buffer);
                hist[si][ki] = new TH3D((name+"_" + ktag[ki]).c_str(),(name+"_" +ktag[ki]).c_str(), NpT, pTmin-dpT/2, pTmax+dpT/2,Nbeta, betamin-dbeta/2, betamax-dbeta/2,NT, Tmin-dT/2,Tmax+dT/2);
                hist[si][ki]->Write();
            }
        }
        for (int si =0; si<ns; si++){

        sprintf (buffer, "PDGid_%d", pdg[si]);
        string name(buffer);

            cout << " Do particle  " <<tag1+name+"_total_T"+Ttag  <<endl;
            TParticle_AZYHYDRO particle(tag1+name+"_total_T"+Ttag);
            TKernel kernel(&particle);
            for (int i =1; i<=hist[0][0]->GetXaxis()->GetNbins(); i++){ 
                double pT = hist[0][0] ->GetXaxis()->GetBinCenter(i);
                for (int j =1; j<=hist[0][0]->GetYaxis()->GetNbins(); j++){ 
                    double beta = hist[0][0] ->GetYaxis()->GetBinCenter(j);
                    double ur = beta/sqrt(1-beta*beta);
                    // increase lambdas by two because neutral (particle and antiparticle)
                    if (pdg[si]==3312){ 
                        for (int ki =0; ki<Nk; ki++){ 
                    double K_mu =kernel.get_K_mu_qag( pT, ur, ki);
                          hist[si][ki]->Fill(pT,beta,Tfo,2*K_mu); }
                    } else {
                        for (int ki =0; ki<Nk; ki++){
                    double K_mu =kernel.get_K_mu_qag( pT, ur, ki);
                          hist[si][ki]->Fill(pT,beta,Tfo,K_mu); }
                    }
                }
            }
        }

        for (int ki =0; ki<Nk; ki++){ 
            for (int si =0; si<ns; si++){
                hist[si][ki]->Write("",TObject::kOverwrite);
            }
        }

        f->Close();
        //    f->Write();
    }
    TH3D *tothist[ns][Nk];
    for (int ki =0; ki<Nk; ki++){ 
        char buffer [150];
        for (int si =0; si<ns; si++){ 

        sprintf (buffer, "PDGid_%d", pdg[si]);
        string name(buffer);
            tothist[si][ki] = new TH3D((name+"_"+ktag[ki]).c_str(),(name+"_"+ktag[ki]).c_str(), NpT, pTmin-dpT/2, pTmax+dpT/2,Nbeta, betamin-dbeta/2, betamax-dbeta/2,NT, Tmin-dT/2,Tmax+dT/2);
            tothist[si][ki]->SetDirectory(0);
        }
    }

    for (int k =0; k <NT; k++){
        //        double Tfo = h3[0] ->GetZaxis()->GetBinCenter(k);
        double Tfo=Tmin+dT*k; //GeV
        char buffer [150];
        sprintf (buffer, "%.4f", Tfo);
        string Ttag(buffer);    
        cout << " Read kernels = " <<Ttag  <<" GeV" <<endl;
        TFile * f = TFile::Open((tag2+"K_total_full_T"+Ttag+".root").c_str());
        TH3D *hist[ns][Nk];
        for (int ki =0; ki<Nk; ki++){ 
            for (int si =0; si<ns; si++){ 

        sprintf (buffer, "PDGid_%d", pdg[si]);
        string name(buffer);
                hist[si][ki] = (TH3D *) f -> Get((name+"_" +ktag[ki]).c_str());
                tothist[si][ki]->Add(hist[si][ki]);
            }
        }

        f->Close();
    }

    TFile * totf = TFile::Open((tag2+"K_total_full.root").c_str(),"RECREATE");
    for (int ki =0; ki<Nk; ki++){ 
        for (int si =0; si<ns; si++){ tothist[si][ki]->Write(); }}
    totf->Close();

}
