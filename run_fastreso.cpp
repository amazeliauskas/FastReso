/*
 * Copyright (c) 2018 Aleksas Mazeliauskas
 * All rights reserved.
 *
 * FastReso is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/amazeliauskas/FastReso/
 */

#include <iostream>
#include <TFastReso.h>

using namespace std ;


int main()
{

  // freeze-out temperature
  double Tfo=0.145; //GeV
  TFastReso fastreso;
  // read all particles
  fastreso.read_particles_data("particles.data");
  // initialize irreducible functions with freeze-out parameters
  // void do_thermal(double Tfo, double MuB=0, double MuI3=0, double MuS=0, double MuC=0, double Cs2=0.14);
  fastreso.do_thermal(Tfo);
  cout << " Done Thermal " << endl;
  // print out thermal irreducible functions of some particles
  //fastreso.getParticle("pi0139plu")->print("pi0139plu_thermal");
  //fastreso.getParticle("om0782zer")->print("om0782zer_thermal");
  //fastreso.getParticle("rho770zer")->print("rho770zer_thermal");
  // read and do decays (if you want information from partial decay chains, put print
  // out inside the decay loop
  fastreso.do_decays("decays.data");
  // print out final rreducible functions
  fastreso.getParticle("pr0938plu")->print("pr0938plu_total");
  fastreso.getParticle("pr0938plb")->print("pr0938plb_total");
  fastreso.getParticle("Ka0492plu")->print("Ka0492plu_total");
  fastreso.getParticle("Ka0492min")->print("Ka0492min_total");
  fastreso.getParticle("pi0139plu")->print("pi0139plu_total");
  fastreso.getParticle("pi0139min")->print("pi0139min_total");
}
