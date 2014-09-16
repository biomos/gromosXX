//GP AT WORK 22
/**
 * @file in_rdc.cc
 * implements methods of In_RDC
 */

#include <stdheader.h>
#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <interaction/interaction_types.h>
#include <configuration/configuration.h>
#include <io/instream.h>
#include <io/blockinput.h>
#include <vector>
#include <iosfwd>
#include <math/random.h>
#include "in_rdc.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

using namespace std;

/**
 * @section conversion CONVERSION block
 * The CONVERSION block is read from the RDC restraint specification file.
 *
 * - Unit conversion factor from Hz to ps-1.
 * - Unit conversion factor from 10^7 (rad /T s) to (e/u).
 *
 * @verbatim

CONVERSION
 # factors
 # to convert the frequency from [RDC]=(s-1) to (ps-1)
 # and to convert gyromagnetic ratios from [gamma]=10^6*(rad/T s)=10^6(C/kg) to (e/u)
 0.000000000001
 0.0103643533
END
@endverbatim

 * @section magfieldc MAGFIELDC block
 * The MAGFIELDC block is read from the RDC restraint specification file. You need
 * this section if you choose NTRDCT 0 (cartesian representation of magnetic field
 * vectors) in your gromos configuration file.
 *
 * - Variable \c NMF is the number of magnetic field direction vectors
 * - Variables \c x, \c y, \z are the cartesian coordinates of the moveable pseudo-atoms representing
 *   the magnetic field vectors. They specify the direction of the magnetic
 *   field vector, as the other coordinate is always (0, 0, 0). The magnetic field vectors
 *   are normalised when read in.
 * - Variable \c mass is the mass of each moveable atom and can be used as a damping
 *   factor.
 *
 * @verbatim
MAGFIELDC
# NMF
    2
#      x      y      z    mass
     1.0    0.0    0.0     1.0
     0.0    1.0    0.0     1.0
END
 @endverbatim

 * @section alignt ALIGNT block
 * The ALIGNT block is read from the RDC restraint specification file. You need
 * this section if you choose NTRDCT 1 (tensor representation of the alignment of
 * the molecule in the magnetic field) in your gromos configuration file.
 *
 * - Variables \c Axx, \c Ayy, \c Axy, \c Axz, \c Ayz, are the five independent tensor
 *   components of the 3x3 alignment tensor representing the alignment of the molecule
 *   in the magnetic field.
 * - Variable \c mass is the pseudo mass of each component and can be used as a damping
 *   factor.
 * - Variable \c vel is the velocity of each component
 *
 * @verbatim
ALIGNT
#    Axx mass1 Ayy mass2 Axy mass3 Axz mass4 Ayz mass5
     0.5   1.0 0.5   1.0 0.5   1.0 0.5   1.0 0.5   1.0
END
 @endverbatim

"# For each RDC restraint the following should be specified:\n"
"# IPRDCR, JPRDCR, KPRDCR, LPRDCR, atom numbers (defining the vector that forms the angle with the magnetic field)\n"
"# WRDCR                           weight factor (for weighting some RDCs higher than others)\n"
"# PRDCR0                          RDC restraint value (i.e. the experimental data)\n"
"# RDCGI, RDCGJ                    gyromagnetic ratios for atoms i and j\n"
"# RDCRIJ, RDCRIK                  distance between atoms i and j or i and k (RIJ = RCH for CA:HA)\n"
"# TYPE                            code to define type of RDC\n"
"# IPRDCR JPRDCR KPRDCR LPRDCR   WRDCR    PRDCR0      RDCGI      RDCGJ     RDCRIJ     RDCRIK    RDCTYPE\n";

 * @section rdcresspec RDCRESSPEC block
 * The RDCRESSPEC block is read from the RDC restraint specification file.
 *
 * - Variables \c IPRDCR, \c JPRDCR, \c KPRDCR, \c LPRDCR are atom sequence numbers
 *   defining the vector that forms the angle with the magnetic field. Only IPRDCR
 *   and JPRDCR are read in (KPRDCR and LPRDCR are used by the gromos++ program svd_fit)
 * - Variable \c WRDCR is an individual RDC restraint weight factor by which
 *   the RDC restraining term for a specific RDC may be multiplied
 * - Variable \c PRDCR0 is the experimental or reference RDC value in Hz.
 * - Variables \c RDCGI, \c RDCGJ are the gyromagnetic ratios of atoms i and j in
 *   10^-7 (rad/T s)
 * - Variables \c RDCRIJ, \c RDCIK are user-defined distances between atoms i and j or i and k
 *   (these are only for use in the gromos++ program svd_fit)
 * - Variable \c TYPE is a code to specify the type of RDC. Only codes 1 (N:H),
 *   2 (CA:C) and 3 (C:N) are valid for RDC restraining.
 *
 * @verbatim
RDCRESSPEC
# For each RDC restraint the following is to be specified:
# IPRDCR, JPRDCR, KPRDCR, LPRDCR    atom sequence numbers
# WRDCR                             weight factor
# PRDCR0                            RDC value
# RDCGI, RDCGJ                      gyromagnetic ratios
# RDCRIJ, RDCRIK                    inter-nuclear distances
# TYPE                              type of RDC
# IPRDCR JPRDCR KPRDCR LPRDCR   WRDCR    PRDCR0      RDCGI      RDCGJ     RDCRIJ     RDCRIK    RDCTYPE
      16     17      0      0       1     -1.78      2.712     267.52      0.104      0.104          1
      24     25      0      0       1      7.52      2.712     267.52      0.104      0.104          1
      41     42      0      0       1     -6.92      2.712     267.52      0.104      0.104          1
      46     47      0      0       1     -6.47      2.712     267.52      0.104      0.104          1
 * END
@endverbatim
 */

double generate_boltzmann_velocities(math::RandomGenerator* rng, const double temp, const double mass) {
  const double sd = sqrt(math::k_Boltzmann * temp / mass);
  rng->stddev(sd);
  return rng->get_gauss();
}

// FIXME
// the following validity test is highly inelegant and misses some invalid a,
// but at least all the rejected ones *are* invalid.
bool a_is_valid(double a1, double a2, double a3, double a4, double a5){

  a1 = 2.0/3.0*(a1+.5);
  a2 = 2.0/3.0*(a2+.5);
  a3 = 2.0/3.0*a3;
  a4 = 2.0/3.0*a4;
  a5 = 2.0/3.0*a5;

  if(a1 <  0.0 || a1 > 1.5) return false;
  if(a2 <  0.0 || a2 > 1.5) return false;
  if(a3 < -0.5 || a3 > 0.5) return false;
  if(a4 < -0.5 || a4 > 0.5) return false;
  if(a5 < -0.5 || a5 > 0.5) return false;

  if(a3 > sqrt(a1)*sqrt(a2) || a3 < -sqrt(a1)*sqrt(a2)) return false;

//  if(pow(Axx + 0.5 ,2) + pow(Ayy + 0.5 ,2) + 2* pow(Axy ,2) + 2* pow(Axz ,2) + 2* pow(Ayz ,2) > 2.25) return false;
  if(pow(a1 ,2) + pow(a2 ,2) + 2* pow(a3 ,2) + 2* pow(a4 ,2) + 2* pow(a5 ,2) > 1.0) return false;

  return true;
}

math::VArray create_points_on_sphere(const unsigned int N){
  if (N<2) {
    io::messages.add("Please choose a number of at least 3 magnetic field vectors", "In_RDC", io::message::critical);
  }
  math::VArray coordinates;
  switch(N){
    case 3:{
      coordinates.push_back(math::Vec(1, 0, 0));
      coordinates.push_back(math::Vec(0, 1, 0));
      coordinates.push_back(math::Vec(0, 0, 1));
    break;
    }
    case 4:{
      const double a = 1.0/sqrt(3.0);
      coordinates.push_back(math::Vec( a, a, a));
      coordinates.push_back(math::Vec( a,-a,-a));
      coordinates.push_back(math::Vec(-a, a,-a));
      coordinates.push_back(math::Vec(-a,-a, a));
    break;
    }
    case 5:{ //FIXME, this is not a good distribution
      coordinates.push_back(math::Vec( 0,   0,             1));
      coordinates.push_back(math::Vec( 1,   0,             0));
      coordinates.push_back(math::Vec(-0.5, cos(M_PI/6.0), 0));
      coordinates.push_back(math::Vec(-0.5,-cos(M_PI/6.0), 0));
      coordinates.push_back(math::Vec( 0,   0,            -1));
    break;
    }
    case 6:{
      coordinates.push_back(math::Vec( 1, 0, 0));
      coordinates.push_back(math::Vec( 0, 1, 0));
      coordinates.push_back(math::Vec( 0, 0, 1));
      coordinates.push_back(math::Vec(-1, 0, 0));
      coordinates.push_back(math::Vec( 0,-1, 0));
      coordinates.push_back(math::Vec( 0, 0,-1));
    break;
    }
    default:{
      const double s = 3.6/sqrt(N);
      const double dz = 2.0/N;
      double longitude = 0.0;
      double z = 1.0 - dz/2.0;
      for (unsigned int i=0; i<N; ++i){
        double r = sqrt(1.0 - z*z);
        coordinates.push_back(math::Vec(cos(longitude)*r, sin(longitude)*r, z));
        DEBUG(10, cos(longitude)*r << ", " <<  sin(longitude)*r << ", " << z)
        z -= dz;
        longitude += s/r;
      }
    }
  }
  return coordinates;
}


void io::In_RDC::read(topology::Topology& topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        ostream & os) {

  DEBUG(7, "reading in a rdc restraints specification file")

#ifdef DEBUG
      // Create random number generator, set seed to 0 (to make sure it's reproducible)
      math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), "0");
#else
      // Create random number generator, choose a remotely random random-seed
      ostringstream ss;
      ss << time(NULL);
      math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), ss.str());
#endif // DEBUG

  const double magic_value = -27.1;

  //////////////////////
  // CONVERSION BLOCK //
  //////////////////////

  {
    // get two conversion factors from (input) to ps-1 (internal)  typically input is Hz and the factor is 10^10
    //                        and from (input) to e/u (internal)   typically input is 10^7*rad/T*s and the factor is 0.10375

    DEBUG(10, "RDC CONVERSION block")
    vector<string> bufferC = m_block["CONVERSION"];
    if (!bufferC.size()) {
      io::messages.add("no CONVERSION block in RDC restraints file", "In_RDC", io::message::critical);
    }
    double factor;
    int swtch = 0;
    os.precision(12);
    os.setf(ios_base::fixed, ios_base::floatfield);

//FIXME this is insanely complicated

    vector<string>::const_iterator itC = bufferC.begin() + 1, toC = bufferC.end() - 1;
    for (; itC != toC; ++itC) {
      _lineStream.clear();
      _lineStream.str(*itC);
      _lineStream >> factor;
      if (!swtch) {
        DEBUG(10, setw(14) << "factorFreq: " << factor)
        conf.special().rdc.factorFreq = factor;
        swtch++;
      } else {
        DEBUG(10, setprecision(4) << setw(14) << "factorGyr: " << factor)
        conf.special().rdc.factorGyr = factor;
      }
    }
    DEBUG(10, "END")
  }
  // convert DeltaRDC to internal GROMOS units
  //GP TODO check DeltaRDC
  sim.param().rdc.delta = sim.param().rdc.delta * conf.special().rdc.factorFreq;

  DEBUG(10, "DeltaRDC (in Gromos units) " << sim.param().rdc.delta)


  //////////////////
  //  INPUTMODE  //
  //////////////////

  // there is a 'basic' input mode ('0') in which most settings are chosen for
  // the user and an expert mode ('1') in which more choices can be made by the
  // user

  DEBUG(10, "RDC INPUTMODE block")
  vector<string> bufferC = m_block["INPUTMODE"];
  if (!bufferC.size()) {
    io::messages.add("no INPUTMODE block in RDC restraints file", "In_RDC", io::message::critical);
  }
  unsigned int mode = 0;
  os.precision(12);
  os.setf(ios_base::fixed, ios_base::floatfield);

  vector<string>::const_iterator itC = bufferC.begin() + 1, toC = bufferC.end() - 1;

  _lineStream.clear();
  _lineStream.str(*itC);
  _lineStream >> mode;

  if(_lineStream.fail()){
    io::messages.add("bad line in INPUTMODE block: failed to read mode", "In_RDC", io::message::error);
  }

  DEBUG(10, "chosen inputmode: " << mode)

  DEBUG(10, "END")


  switch(sim.param().rdc.type){

    /////////////////////
    // MAGFIELDC BLOCK //
    /////////////////////

    case simulation::rdc_mf: {
      DEBUG(10, "RDC MAGFIELDC block")
      vector<string> bufferMF = m_block["MAGFIELDC"];

      if(mode==0){ // basic
        if (!bufferMF.size()) { // one line for nmf but no coordinates
          io::messages.add("no or empty MAGFIELDC block in RDC restraints file", "In_RDC", io::message::error);
        }
        else {
          vector<string>::const_iterator itMF = bufferMF.begin() + 1,
                  toMF = bufferMF.end() - 1;

          unsigned int nmf;
          unsigned int mass;

          _lineStream.clear();
          _lineStream.str(*itMF);
          _lineStream >> nmf >> mass;

          if (_lineStream.fail()) {
            io::messages.add("bad line in MAGFIELDC block: failed to read in NMF", "In_RDC", io::message::error);
          }

          os.setf(ios_base::fixed, ios_base::floatfield);
          DEBUG(10, setprecision(4) << setw(8) << "NMF" << setw(8) << "mass")
          DEBUG(10, setprecision(4) << setw(8) << nmf   << setw(8) << mass  )

          conf.special().rdc.MFpoint = create_points_on_sphere(nmf);

          DEBUG(10, setw(13) << "x" << setw(13) << "vx" << setw(13) << "y" << setw(13) << "vy" << setw(13) << "z" << setw(13) << "vz")
          for (unsigned int i=0; i<nmf; ++i) {
            const double randomvel1 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass); // [nm/ps]
            const double randomvel2 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass); // [nm/ps]
            const double randomvel3 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass); // [nm/ps]
            const math::Vec tmpVecV(randomvel1, randomvel2, randomvel3);

            conf.special().rdc.MFpointVel.push_back(tmpVecV);
            conf.special().rdc.MFpointMass.push_back(mass);

            DEBUG(10, setprecision(8)
                 << setw(13) << conf.special().rdc.MFpoint[i][0]
                 << setw(13) << tmpVecV[0]
                 << setw(13) << conf.special().rdc.MFpoint[i][1]
                 << setw(13) << tmpVecV[1]
                 << setw(13) << conf.special().rdc.MFpoint[i][2]
                 << setw(13) << tmpVecV[2])
          }
        }
      }else if (mode==1){ // advanced

        if (bufferMF.size() <3) { // one line for nmf but no coordinates
            io::messages.add("no, empty or incomplete (at least 3 vectors required) MAGFIELDC block in RDC restraints file", "In_RDC", io::message::error);
        } else {
          DEBUG(10, setw(13) << "x" << setw(13) << "vx" << setw(13) << "y" << setw(13) << "vy" << setw(13) << "z" << setw(13) << "vz")
          vector<string>::const_iterator itMF = bufferMF.begin() + 1,
                  toMF = bufferMF.end() - 1;
          for (; itMF != toMF; ++itMF){
            double MFx = 0.0, MFy = 0.0, MFz = 0.0, MFm = 0.0; // in units of nm, nm, nm and u
            _lineStream.clear();
            _lineStream.str(*itMF);
            _lineStream >> MFx >> MFy >> MFz >> MFm;

            if(_lineStream.fail()){
                 io::messages.add("bad line in ALIGNT block", "In_RDC", io::message::error);
            }

            math::Vec tmpVecP(MFx, MFy, MFz);
            tmpVecP = tmpVecP.norm(); // normalise

            const double randomvel1 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, MFm); // [nm/ps]
            const double randomvel2 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, MFm); // [nm/ps]
            const double randomvel3 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, MFm); // [nm/ps]

            math::Vec tmpVecV(randomvel1, randomvel2, randomvel3);
            conf.special().rdc.MFpoint.push_back(tmpVecP);
            conf.special().rdc.MFpointVel.push_back(tmpVecV);
            conf.special().rdc.MFpointMass.push_back(MFm);

            DEBUG(10, setw(13) << tmpVecP[0] << setw(13) << tmpVecV[0] << setw(13) << tmpVecP[1] << setw(13) << tmpVecV[1] << setw(13) << tmpVecP[2] << setw(13) << tmpVecV[2] << setw(13) << MFm)
          }
        }
      }
      conf.special().rdc.stochastic_integral_mf.resize(conf.special().rdc.MFpoint.size());

      DEBUG(10, "END")
  // FIXME check for remaining chars in the buffer
      break;
    }

  ///////////////////
  // ALIGNT BLOCK //
  ///////////////////

    case simulation::rdc_t: { // only parse the tensor block if we use it
      DEBUG(10, "RDC ALIGNT block")
      vector<string> bufferA = m_block["ALIGNT"];
      if(mode==0){ // simple input

        if (!bufferA.size()) {
          io::messages.add("no or empty ALIGNT block in RDC restraints file", "In_RDC", io::message::warning);
        }

        vector<string>::const_iterator itA = bufferA.begin() + 1, toA = bufferA.end() - 1;

        DEBUG(10, "reading in ALIGNT data")

        os.setf(ios_base::fixed, ios_base::floatfield);

        double mass;
        double Axx=0, vAxx, Ayy=0, vAyy, Axy=0, vAxy, Axz=0, vAxz, Ayz=0, vAyz;

        _lineStream.clear();
        _lineStream.str(*itA);
        _lineStream >> mass;

        if(_lineStream.fail()){
             io::messages.add("bad line in ALIGNT block", "In_RDC", io::message::error);
        }

        // Generate Boltzmann distributed velocities
        vAxx = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass);
        vAyy = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass);
        vAxy = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass);
        vAxz = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass);
        vAyz = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass);

        DEBUG(10, setw(13) << "Axx" << setw(13) << "vAxx" << setw(13) << "mAxx"
              << setw(13) << "Ayy" << setw(13) << "vAyy" << setw(13) << "mAyy"
              << setw(13) << "Axy" << setw(13) << "vAxy" << setw(13) << "mAxy"
              << setw(13) << "Axz" << setw(13) << "vAxz" << setw(13) << "mAxz"
	          << setw(13) << "Ayz" << setw(13) << "vAyz" << setw(13) << "mAyz")

        DEBUG(10, setprecision(4)
              << setw(13) << Axx << setw(13) << vAxx << setw(13) << mass
              << setw(13) << Ayy << setw(13) << vAyy << setw(13) << mass
              << setw(13) << Axy << setw(13) << vAxy << setw(13) << mass
              << setw(13) << Axz << setw(13) << vAxz << setw(13) << mass
              << setw(13) << Ayz << setw(13) << vAyz << setw(13) << mass)

        conf.special().rdc.Tensor.push_back(Axx);
        conf.special().rdc.Tensor.push_back(Ayy);
        conf.special().rdc.Tensor.push_back(Axy);
        conf.special().rdc.Tensor.push_back(Axz);
        conf.special().rdc.Tensor.push_back(Ayz);

        conf.special().rdc.TensorVel.push_back(vAxx);
        conf.special().rdc.TensorVel.push_back(vAyy);
        conf.special().rdc.TensorVel.push_back(vAxy);
        conf.special().rdc.TensorVel.push_back(vAxz);
        conf.special().rdc.TensorVel.push_back(vAyz);

        conf.special().rdc.TensorMass.push_back(mass);
        conf.special().rdc.TensorMass.push_back(mass);
        conf.special().rdc.TensorMass.push_back(mass);
        conf.special().rdc.TensorMass.push_back(mass);
        conf.special().rdc.TensorMass.push_back(mass);

        // set first c value to -2 (which is an invalid value) so signal to the
        // MD or SD routine, that correct values should be calculated.  In case
        // of EM it is overwritten anyway.
        conf.special().rdc.Tensor[0] = -2;

      }else if(mode==1){ // advanced input

        if (!bufferA.size()) {
          io::messages.add("no or empty ALIGNT block in RDC restraints file", "In_RDC", io::message::warning);
        }

        vector<string>::const_iterator itA = bufferA.begin() + 1, toA = bufferA.end() - 1;

        DEBUG(10, "reading in ALIGNT data")

        os.setf(ios_base::fixed, ios_base::floatfield);

        double Axx, mAxx, vAxx, Ayy, mAyy, vAyy, Axy, mAxy, vAxy, Axz, mAxz, vAxz, Ayz, mAyz=magic_value, vAyz; // there is no straight forward way to check how many numbers were read on a line.  This is ugly, but short

        _lineStream.clear();
        _lineStream.str(*itA);
        _lineStream >> Axx >> mAxx >> Ayy >> mAyy >> Axy >> mAxy >> Axz >> mAxz >> Ayz >> mAyz;

        if(_lineStream.fail()){
             io::messages.add("bad line in ALIGNT block", "In_RDC", io::message::error);
        }
        if(abs(mAyz-magic_value) < 1.e-10) io::messages.add("incomplete ALIGN_T block in RDC restraints file (less than 5 values and masses read.)", "In_RDC", io::message::error);
        if(mAxx<0.0 || mAyy<0.0 || mAxy<0.0 || mAxz<0.0 || mAyz<0.0 ) io::messages.add("masses can only be positive", "In_RDC", io::message::error);

	    // the tensor components a_1 to a_5 can be expressed as averages over
	    // products of cosines and are, hence, limited to a range of values
        if(!a_is_valid(Axx, Ayy, Axy, Axz, Ayz)){
          io::messages.add("some or all a_h have insensible values.  Setting all components to 0.0 (i.e. isotropic alignment) ...", "In_RDC", io::message::warning);
          Axx = 0.0;
          Ayy = 0.0;
          Axy = 0.0;
          Axz = 0.0;
          Ayz = 0.0;
        }

        // Generate Boltzmann distributed velocities
        vAxx = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAxx);
        vAyy = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAyy);
        vAxy = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAxy);
        vAxz = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAxz);
        vAyz = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAyz);

        DEBUG(10, setw(13) << "Axx" << setw(13) << "vAxx" << setw(13) << "mAxx"
              << setw(13) << "Ayy" << setw(13) << "vAyy" << setw(13) << "mAyy"
              << setw(13) << "Axy" << setw(13) << "vAxy" << setw(13) << "mAxy"
              << setw(13) << "Axz" << setw(13) << "vAxz" << setw(13) << "mAxz"
	          << setw(13) << "Ayz" << setw(13) << "vAyz" << setw(13) << "mAyz")

        DEBUG(10, setprecision(4)
              << setw(13) << Axx << setw(13) << vAxx << setw(13) << mAxx
              << setw(13) << Ayy << setw(13) << vAyy << setw(13) << mAyy
              << setw(13) << Axy << setw(13) << vAxy << setw(13) << mAxy
              << setw(13) << Axz << setw(13) << vAxz << setw(13) << mAxz
              << setw(13) << Ayz << setw(13) << vAyz << setw(13) << mAyz)

        conf.special().rdc.Tensor.push_back(Axx);
        conf.special().rdc.Tensor.push_back(Ayy);
        conf.special().rdc.Tensor.push_back(Axy);
        conf.special().rdc.Tensor.push_back(Axz);
        conf.special().rdc.Tensor.push_back(Ayz);

        conf.special().rdc.TensorVel.push_back(vAxx);
        conf.special().rdc.TensorVel.push_back(vAyy);
        conf.special().rdc.TensorVel.push_back(vAxy);
        conf.special().rdc.TensorVel.push_back(vAxz);
        conf.special().rdc.TensorVel.push_back(vAyz);

        conf.special().rdc.TensorMass.push_back(mAxx);
        conf.special().rdc.TensorMass.push_back(mAyy);
        conf.special().rdc.TensorMass.push_back(mAxy);
        conf.special().rdc.TensorMass.push_back(mAxz);
        conf.special().rdc.TensorMass.push_back(mAyz);
      }

      conf.special().rdc.stochastic_integral_t.resize(conf.special().rdc.Tensor.size());

      DEBUG(10, "END")
  //   FIXME check for remaining chars in the buffer
      break;
    }   // ALIGNT


  ///////////////////////////////
  // SPHERICAL HARMONICS BLOCK //
  ///////////////////////////////

    case simulation::rdc_sh: { // only parse the spherical-harmonics block if we use it
      DEBUG(10, "RDC SPHERICAL HARMONICS block")
      vector<string> bufferSH = m_block["SPHERICALHARMONICS"];
      if(mode==0){ // basic input
        if (!bufferSH.size()) {
          io::messages.add("no or empty SPHERICALHARMONICS block in RDC restraints file", "In_RDC", io::message::warning);
        }

        vector<string>::const_iterator itSH = bufferSH.begin() + 1, toA = bufferSH.end() - 1;

        DEBUG(10, "reading in SPHERICALHARMONICS data")

        os.setf(ios_base::fixed, ios_base::floatfield);

        const int n_clm = 5; // (-2,2), (-1,2), (0,2), (1,2), (2,2)
        double mass;
        vector<double> c;

        _lineStream.clear();
        _lineStream.str(*itSH);

        _lineStream >> mass;

        if(_lineStream.fail()){
          io::messages.add("bad line in SPHERICALHARMONICS block: failed to read 'm(l)'", "In_RDC", io::message::error);
        }

        // resize clm
        conf.special().rdc.clm.resize(n_clm, 0.0);

        // resize + init masses
        conf.special().rdc.clmMass.resize(n_clm,mass);
        // resize + init velocities
        conf.special().rdc.clmVel.resize(n_clm);
        DEBUG(10, setw(13) << "c" << setw(13) << "mass" <<  setw(13) << "velocity")
        for(int i=0; i<n_clm; i++){
          conf.special().rdc.clmVel[i]=generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass);
          DEBUG(10, setprecision(8) << setw(13) << conf.special().rdc.clm[i] << setw(13) << conf.special().rdc.clmMass[i] << setw(13) << conf.special().rdc.clmVel[i])
        }

        // set first c value to -2 (which is an invalid value) so signal to the
        // MD or SD routine, that correct values should be calculated.  In case
        // of EM it is overwritten anyway.
        conf.special().rdc.clm[0] = -2;

      }else if (mode==1){ // advanced input

        if (!bufferSH.size()) {
          io::messages.add("no or empty SPHERICALHARMONICS block in RDC restraints file", "In_RDC", io::message::error);
        }

        vector<string>::const_iterator itSH = bufferSH.begin() + 1, toA = bufferSH.end() - 1;

        DEBUG(10, "reading in SPERICAL HARMONICS data")

        const int n_clm = 5; // (-2,2), (-1,2), (0,2), (1,2), (2,2)
        vector<double> c(5,0.0), mass(5,0.0);

        _lineStream.clear();
        _lineStream.str(*itSH);

        mass[4] = magic_value;
        _lineStream >> c[0] >> mass[0] >> c[1] >> mass[1] >> c[2] >> mass[2] >> c[3] >> mass[3] >> c[4] >> mass[4];

        if(_lineStream.fail()){
          io::messages.add("bad line in SPHERICALHARMONICS block: failed to read in five values c and five masses", "In_RDC", io::message::error);
        }
        if(abs(mass[4]-magic_value) < 1.e-10) io::messages.add("incomplete SPHERICALHARMONICS block in RDC restraints file (less than 5 values and masses read.)", "In_RDC", io::message::error);
        if(mass[0]<0.0 || mass[1]<0.0 || mass[2]<0.0 || mass[3]<0.0 || mass[4]<0.0 ) io::messages.add("masses can only be positive", "In_RDC", io::message::error);

        conf.special().rdc.clm = c;
        conf.special().rdc.clmMass = mass;

        // resize + init velocities
        conf.special().rdc.clmVel.resize(n_clm);
        DEBUG(10, setw(13) << "c" << setw(13) << "mass" << setw(13) << "velocity")
        for(int i=0; i<n_clm;i++){
          conf.special().rdc.clmVel[i]=generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass[i]);
          DEBUG(10, setprecision(8) << setw(13) << conf.special().rdc.clm[i] << setw(13) << conf.special().rdc.clmMass[i] <<  setw(13) << conf.special().rdc.clmVel[i])
        }

        conf.special().rdc.stochastic_integral_sh.resize(n_clm);
      }

      DEBUG(10, "END")
      break;
    } // SH
    default:
      io::messages.add("no valid method chosen", "In_RDC", io::message::critical);
  }


  //////////////////////
  // RDCRESSPEC BLOCK //
  //////////////////////

  vector<string> buffer;
  DEBUG(10, "RDC RESTRAINTS")
  { // RDCRESSPEC
    DEBUG(10, "RDCRESSPEC block")
    buffer = m_block["RDCRESSPEC"];

    if (!buffer.size()){
      io::messages.add("no RDCRESSPEC block in RDC restraints file", "In_RDC", io::message::critical);
    }

    DEBUG(10, "reading in RDCRESSPEC data")

    int i, j, k, l, type;
    double weight, R, gyri, gyrj, rij, rik;

    conf.special().rdc.av.clear();

    os.setf(ios_base::fixed, ios_base::floatfield);

    DEBUG(10, setw(6) << "i"
	       << setw(6) << "j"
	       << setw(8) << "weight"
           << setw(19) << "R"
           << setw(19) << "RDCav"
           << setw(9) << "GYRi"
           << setw(9) << "GYRj")

    vector<string>::const_iterator it = buffer.begin()+1,
    to = buffer.end()-1;

    for(int n=0; it != to; ++it, ++n){

      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> i >> j >> k >> l >> weight >> R >> gyri >> gyrj >> rij >> rik >> type;

      if(_lineStream.fail()){
           io::messages.add("bad line in RDCRESSPEC block", "In_RDC", io::message::error);
      }

      // CONVERT RDC from <input> to 1/ps
      R *= conf.special().rdc.factorFreq;

      // CONVERT GYR from <input> to (e/u)
      gyri *= conf.special().rdc.factorGyr;
      gyrj *= conf.special().rdc.factorGyr;

      // check for sensible choice of RDC
      if(abs(-( math::eps0_i * math::h_bar * gyri * gyrj )/( pow(math::spd_l,2) * 4.0 * pow(math::Pi,2)) * 1000) < abs(R)){
        io::messages.add("The chosen RDC is larger in magnitude than RDC_max.  This is probably a mistake and may result in strange behaviour.",
         "In_RDC", io::message::warning);
      }

      topo.rdc_restraints().push_back(topology::rdc_restraint_struct(i-1, j-1, weight, R, gyri, gyrj)); // in units of 1, 1, 1, 1/ps, e/u and e/u
      conf.special().rdc.av.push_back(R); // init history as experimental values

      DEBUG(10, setprecision(4) << setw(6) << i
                                << setw(6) << j
                                << setw(8) << weight
            << setprecision(14) << setw(19) << R
                                << setw(19) << conf.special().rdc.av[n]
            << setprecision(4)  << setw(9) << gyri
                                << setw(9) << gyrj)
    }

    if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
       sim.param().rdc.mode == simulation::rdc_restr_av ||
       sim.param().rdc.mode == simulation::rdc_restr_biq){
      // as we don't weight any RDCs, we set all weights to one
      // no further distinction is required then
      io::messages.add("No weighting selected: Setting all weights to 1.0", "In_RDC", io::message::notice);
      vector<topology::rdc_restraint_struct>::iterator
        it = topo.rdc_restraints().begin(),
        to = topo.rdc_restraints().end();
      for(; it!=to; ++it){ it->weight = 1.0; }
    }

    conf.special().rdc.curr.resize(topo.rdc_restraints().size());

    DEBUG(10, "END")

  } // RDCRESSPEC
} // io::In_RDC::read

