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
  fastreso.getParticle("pi0139plu")->print("pi0139plu_thermal");
  fastreso.getParticle("om0782zer")->print("om0782zer_thermal");
  fastreso.getParticle("rho770zer")->print("rho770zer_thermal");
  // read and do decays (if you want information from partial decay chains, put print
  // out inside the decay loop
  fastreso.do_decays("decays.data");
  // print out final rreducible functions
  fastreso.getParticle("pi0139plu")->print("pi0139plu_total");
}
