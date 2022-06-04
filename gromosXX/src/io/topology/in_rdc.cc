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

// three structs, local to this file
struct mf_struct{
  math::VArray cart_coords;
  math::VArray cart_vels;
  vector<double> masses;
};

struct t_struct{
  vector<double> a;
  vector<double> vels;
  vector<double> masses;
};

struct sh_struct{
  vector<double> clm;
  vector<double> vels;
  vector<double> masses;
};

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
  if(pow(a1, 2) + pow(a2, 2) + 2* pow(a3, 2) + 2* pow(a4, 2) + 2* pow(a5, 2) > 1.0) return false;

  return true;
}

math::VArray create_points_on_sphere(const unsigned int N){
  if (N<5) {
    io::messages.add("Please choose a number of at least 5 magnetic field vectors, to describe general alignment", "In_RDC", io::message::critical);
  }
  math::VArray coordinates;
  switch(N){
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


/**
 * @section conversion CONVERSION block
 * The CONVERSION block is read from the RDC restraint specification file.
 *
 * - Unit conversion factor from Hz to ps-1.
 * - Unit conversion factor from 10^6 (rad /T s) to (e/u).
 *
 * @verbatim
CONVERSION
 # factors
 # to convert the frequency from [RDC]=(s^-1) to (ps^-1)
 # and to convert gyromagnetic ratios from [gamma]=10^6*(rad/T s)=10^6(C/kg) to (e/u)
 0.000000000001 0.010364272
END
@endverbatim

 * @section magfieldc MAGFIELDC block
 * The MAGFIELDC block is read from the RDC restraint specification file. You need
 * this section if you choose NTRDCT 0 (cartesian representation of magnetic field
 * vectors) in your gromos configuration file.
 *
 * - Variable \c NMF is the number of magnetic field direction vectors
 * - Variables \c x, \c y, \c z are the cartesian coordinates of the moveable pseudo-atoms representing
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
 *   10^6 (rad/T s)
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
      16     17      0      0       1     -1.78      -27.12     267.52      0.104      0.104          1
      24     25      0      0       1      7.52      -27.12     267.52      0.104      0.104          1
      41     42      0      0       1     -6.92      -27.12     267.52      0.104      0.104          1
      46     47      0      0       1     -6.47      -27.12     267.52      0.104      0.104          1
END
@endverbatim
 */

void io::In_RDC::read(topology::Topology& topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        ostream & os) {

  DEBUG(7, "reading in an RDC restraints specification file")

#ifdef DEBUG
      // Create random number generator, set seed to 0 (to make sure it's reproducible)
      math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), "0");
#else
      // Create random number generator, choose a remotely random random-seed
      ostringstream ss;
      ss << time(NULL);
      math::RandomGenerator* rng = math::RandomGenerator::create(sim.param(), ss.str());
#endif // DEBUG

//  const double magic_value = -27.1;


////////////////////////////////////////
//                                    //
//  reading values to temp-variables  //
//                                    //
////////////////////////////////////////

  // buffer for all input blocks
  vector<string> buffer;

  //////////////////////
  // CONVERSION BLOCK //
  //////////////////////
  double factorfreq = 0.0, factorgyr = 0.0;
  {
    // get two conversion factors from (input) to ps-1 (internal)  typically input is Hz and the factor is 10^-12
    //                        and from (input) to e/u (internal)   typically input is 10^7*rad/T*s and the factor is 0.10375

    DEBUG(10, "RDC CONVERSION")

    vector<string> buffer = m_block["CONVERSION"];
    DEBUG(10, "CONVERSION block: " << buffer.size())
    if (!buffer.size()) {
      io::messages.add("no CONVERSION block in RDC restraints file", "In_RDC", io::message::critical);
    }

    os.precision(12);
    os.setf(ios_base::fixed, ios_base::floatfield);

    _lineStream.clear();
    _lineStream.str(*(buffer.begin()+1));
    _lineStream >> factorfreq >> factorgyr;

    DEBUG(10, scientific << setprecision(6) << setw(14) << "factorFreq: " << factorfreq)
    DEBUG(10, scientific << setprecision(6) << setw(14) << "factorGyr: " << factorgyr)

    DEBUG(10, "END")
  }


  //////////////////
  //  INPUTMODE  //
  //////////////////
  // there is a 'basic' input mode ('0') in which most settings are chosen for
  // the user and an expert mode ('1') in which more choices can be made by the
  // user
  unsigned int input_mode = 0;
  {
    DEBUG(10, "RDC INPUTMODE")
    buffer = m_block["INPUTMODE"];
    DEBUG(10, "INPUTMODE block: " << buffer.size())
    if (!buffer.size()) {
      io::messages.add("no INPUTMODE block in RDC restraints file", "In_RDC", io::message::critical);
    }
    os.precision(12);
    os.setf(ios_base::fixed, ios_base::floatfield);

    _lineStream.clear();
    _lineStream.str(*(buffer.begin()+1));
    _lineStream >> input_mode;

    if(_lineStream.fail()){
      io::messages.add("bad line in INPUTMODE block: failed to read input mode", "In_RDC", io::message::error);
    }
    if(input_mode!=0 && input_mode!=1){
      io::messages.add("the only valid input modes are '0' (basic) and '1' (advanced)", "In_RDC", io::message::error);
    }

    DEBUG(10, "chosen inputmode: " << input_mode)

    DEBUG(10, "END")
  }

// FIXME having these lines here is a bit annoying
  mf_struct temp_mf;
  t_struct temp_t;
  sh_struct temp_sh;
  switch(sim.param().rdc.type){

    /////////////////////
    // MAGFIELDC BLOCK //
    /////////////////////

    case simulation::rdc_mf: {
      DEBUG(10, "RDC MAGFIELDC")
      buffer = m_block["MAGFIELDC"];
      DEBUG(10, "MAGFIELDC block: " << buffer.size())

      if(input_mode==0){ // basic
        if (!buffer.size()) { // one line for nmf but no coordinates
          io::messages.add("no or empty MAGFIELDC block in RDC restraints file", "In_RDC", io::message::error);
        }
        else {
          vector<string>::const_iterator itMF = buffer.begin() + 1;

          unsigned int n_mf = 0;
          unsigned int mf_mass = 0;

          _lineStream.clear();
          _lineStream.str(*itMF);
          _lineStream >> n_mf >> mf_mass;

          if (_lineStream.fail()) {
            io::messages.add("bad line in MAGFIELDC block: failed to read in NMF", "In_RDC", io::message::error);
          }
          if (mf_mass < 1e-5) {
            io::messages.add("You chose a mass of < 1e-5 for the mf particles.  This is either an error or a bad idea.", "In_RDC", io::message::error);
          }

          os.setf(ios_base::fixed, ios_base::floatfield);
          DEBUG(10, setprecision(4) << setw(8) << "NMF" << setw(8) << "mass")
          DEBUG(10, setprecision(4) << setw(8) << n_mf   << setw(8) << mf_mass  )

          temp_mf.cart_coords = create_points_on_sphere(n_mf);

//          DEBUG(10, setw(13) << "x" << setw(13) << "vx" << setw(13) << "y" << setw(13) << "vy" << setw(13) << "z" << setw(13) << "vz")
          for (unsigned int i=0; i<n_mf; ++i) {
            const double randomvel1 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]
            const double randomvel2 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]
            const double randomvel3 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]

            temp_mf.cart_vels.push_back(math::Vec(randomvel1, randomvel2, randomvel3));
            temp_mf.masses.push_back(mf_mass);

//            DEBUG(10, setprecision(8)
//                 << setw(13) << conf.special().rdc.MFpoint[i][0]
//                 << setw(13) << tmpVecV[0]
//                 << setw(13) << conf.special().rdc.MFpoint[i][1]
//                 << setw(13) << tmpVecV[1]
//                 << setw(13) << conf.special().rdc.MFpoint[i][2]
//                 << setw(13) << tmpVecV[2])
          }
        }
      }else if (input_mode==1){ // advanced

        if (buffer.size()<3) { // one line for nmf but no coordinates
            io::messages.add("no, empty or incomplete (at least 3 vectors required) MAGFIELDC block in RDC restraints file", "In_RDC", io::message::error);
        } else {
          DEBUG(10, setw(13) << "x" << setw(13) << "vx" << setw(13) << "y" << setw(13) << "vy" << setw(13) << "z" << setw(13) << "vz")
          vector<string>::const_iterator itMF = buffer.begin() + 1, toMF = buffer.end() - 1;
          for (; itMF != toMF; ++itMF){
            double MFx = 0.0, MFy = 0.0, MFz = 0.0, mf_mass = 0.0; // in units of nm, nm, nm and u
            _lineStream.clear();
            _lineStream.str(*itMF);
            _lineStream >> MFx >> MFy >> MFz >> mf_mass;

            if(_lineStream.fail()){
                 io::messages.add("bad line in ALIGNT block", "In_RDC", io::message::error);
            }
            if (mf_mass < 1e-5) {
              io::messages.add("You chose a mass of < 1e-5 for an mf particle.  This is either an error or a bad idea.", "In_RDC", io::message::error);
            }

            const double randomvel1 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]
            const double randomvel2 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]
            const double randomvel3 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]

            temp_mf.cart_coords.push_back(math::Vec(MFx, MFy, MFz).norm());
            temp_mf.cart_vels.push_back(math::Vec(randomvel1, randomvel2, randomvel3));
            temp_mf.masses.push_back(mf_mass);

//            DEBUG(10, setw(13) << tmpVecP[0] << setw(13) << tmpVecV[0] << setw(13) << tmpVecP[1] << setw(13) << tmpVecV[1] << setw(13) << tmpVecP[2] << setw(13) << tmpVecV[2] << setw(13) << mf_mass)
          }
        }
      }

      DEBUG(10, "END")
  // FIXME check for remaining chars in the buffer
      break;
    }

  ///////////////////
  // ALIGNT BLOCK //
  ///////////////////

    case simulation::rdc_t: { // only parse the tensor block if we use it
      DEBUG(10, "RDC ALIGNT block")
      buffer = m_block["ALIGNT"];
      DEBUG(10, "ALIGNT block: " << buffer.size())

      if(input_mode==0){ // simple input

        if (!buffer.size()) {
          io::messages.add("no or empty ALIGNT block in RDC restraints file", "In_RDC", io::message::warning);
        }

        vector<string>::const_iterator itA = buffer.begin() + 1;

        DEBUG(10, "reading in ALIGNT data")

        os.setf(ios_base::fixed, ios_base::floatfield);

        double mass = 0.0;
        double Axx=0, vAxx = 0.0, Ayy=0, vAyy = 0.0, Axy=0, vAxy = 0.0, Axz=0, vAxz = 0.0, Ayz=0, vAyz = 0.0;

        _lineStream.clear();
        _lineStream.str(*itA);
        _lineStream >> mass;

        if(_lineStream.fail()){
             io::messages.add("bad line in ALIGNT block", "In_RDC", io::message::error);
        }
        if (mass < 1e-5) {
          io::messages.add("You chose a mass of < 1e-5 for the tensor component particles.  This is either an error or a bad idea.", "In_RDC", io::message::error);
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

        temp_t.a.push_back(Axx);
        temp_t.a.push_back(Ayy);
        temp_t.a.push_back(Axy);
        temp_t.a.push_back(Axz);
        temp_t.a.push_back(Ayz);

        temp_t.vels.push_back(vAxx);
        temp_t.vels.push_back(vAyy);
        temp_t.vels.push_back(vAxy);
        temp_t.vels.push_back(vAxz);
        temp_t.vels.push_back(vAyz);

        temp_t.masses.push_back(mass);
        temp_t.masses.push_back(mass);
        temp_t.masses.push_back(mass);
        temp_t.masses.push_back(mass);
        temp_t.masses.push_back(mass);

//        // set first c value to -2 (which is an invalid value) so signal to the
//        // MD or SD routine, that correct values should be calculated.  In case
//        // of EM it is overwritten anyway.
//        conf.special().rdc.Tensor[0] = -2;

      }else if(input_mode==1){ // advanced input

        if (!buffer.size()) {
          io::messages.add("no or empty ALIGNT block in RDC restraints file", "In_RDC", io::message::warning);
        }

        vector<string>::const_iterator itA = buffer.begin() + 1;

        DEBUG(10, "reading in ALIGNT data")

        os.setf(ios_base::fixed, ios_base::floatfield);

        double Axx = 0.0, mAxx = 0.0, vAxx = 0.0, Ayy = 0.0, mAyy = 0.0, vAyy = 0.0, Axy = 0.0, mAxy = 0.0, vAxy = 0.0, Axz = 0.0, mAxz = 0.0, vAxz = 0.0, Ayz = 0.0, mAyz = 0.0, vAyz = 0.0;

        _lineStream.clear();
        _lineStream.str(*itA);
        _lineStream >> Axx >> mAxx >> Ayy >> mAyy >> Axy >> mAxy >> Axz >> mAxz >> Ayz >> mAyz;

        if(_lineStream.fail()){
             io::messages.add("bad line in ALIGNT block", "In_RDC", io::message::error);
        }
        if(mAxx<1e-5 || mAyy<1e-5 || mAxy<1e-5 || mAxz<1e-5 || mAyz<1e-5 ) io::messages.add("masses can only be positive and shouldn't be too small", "In_RDC", io::message::error);

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

        temp_t.a.push_back(Axx);
        temp_t.a.push_back(Ayy);
        temp_t.a.push_back(Axy);
        temp_t.a.push_back(Axz);
        temp_t.a.push_back(Ayz);

        temp_t.vels.push_back(vAxx);
        temp_t.vels.push_back(vAyy);
        temp_t.vels.push_back(vAxy);
        temp_t.vels.push_back(vAxz);
        temp_t.vels.push_back(vAyz);

        temp_t.masses.push_back(mAxx);
        temp_t.masses.push_back(mAyy);
        temp_t.masses.push_back(mAxy);
        temp_t.masses.push_back(mAxz);
        temp_t.masses.push_back(mAyz);
      }

      DEBUG(10, "END")
  //   FIXME check for remaining chars in the buffer
      break;
    }   // ALIGNT


  ///////////////////////////////
  // SPHERICAL HARMONICS BLOCK //
  ///////////////////////////////

    case simulation::rdc_sh: { // only parse the spherical-harmonics block if we use it
      DEBUG(10, "RDC SPHERICAL HARMONICS")
      vector<string> buffer = m_block["SPHERICALHARMONICS"];
      DEBUG(10, "SPHERICALHARMONICS block: " << buffer.size())

      if(input_mode==0){ // basic input
        if (!buffer.size()) {
          io::messages.add("no or empty SPHERICALHARMONICS block in RDC restraints file", "In_RDC", io::message::warning);
        }

        vector<string>::const_iterator itSH = buffer.begin() + 1;

        DEBUG(10, "reading in SPHERICALHARMONICS data")

        os.setf(ios_base::fixed, ios_base::floatfield);

        const int n_clm = 5; // (-2,2), (-1,2), (0,2), (1,2), (2,2)
        double mass = 0.0;
        vector<double> c;

        _lineStream.clear();
        _lineStream.str(*itSH);

        _lineStream >> mass;

        if(_lineStream.fail()){
          io::messages.add("bad line in SPHERICALHARMONICS block: failed to read 'm(l)'", "In_RDC", io::message::error);
        }
        if (mass < 1e-5) {
          io::messages.add("You chose a mass of < 1e-5 for the sh particles.  This is either an error or a bad idea.", "In_RDC", io::message::error);
        }


//        // set first c value to -2 (which is an invalid value) so signal to the
//        // MD or SD routine, that correct values should be calculated.  In case
//        // of EM it is overwritten anyway.
//        conf.special().rdc.clm[0] = -2;

        temp_sh.clm.resize(n_clm, 0.0);
        temp_sh.masses.resize(n_clm, mass);
        temp_sh.vels.resize(n_clm);
        for(int i=0; i<n_clm; i++){
          temp_sh.vels[i] = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass);
        }


      }else if (input_mode==1){ // advanced input

        if (!buffer.size()) {
          io::messages.add("no or empty SPHERICALHARMONICS block in RDC restraints file", "In_RDC", io::message::error);
        }

        vector<string>::const_iterator itSH = buffer.begin() + 1;

        DEBUG(10, "reading in SPERICAL HARMONICS data")

        const int n_clm = 5; // (-2,2), (-1,2), (0,2), (1,2), (2,2)
        vector<double> c(5,0.0), mass(5,0.0);

        _lineStream.clear();
        _lineStream.str(*itSH);

        _lineStream >> c[0] >> mass[0] >> c[1] >> mass[1] >> c[2] >> mass[2] >> c[3] >> mass[3] >> c[4] >> mass[4];

        if(_lineStream.fail()){
          io::messages.add("bad line in SPHERICALHARMONICS block: failed to read in five values c and five masses", "In_RDC", io::message::error);
        }
        if(mass[0]<1e-5 || mass[1]<1e-5 || mass[2]<1e-5 || mass[3]<1e-5 || mass[4]<1e-5 ) io::messages.add("masses can only be positive and shouldn't be too small", "In_RDC", io::message::error);


        temp_sh.clm = c;
        temp_sh.masses = mass;
        temp_sh.vels.resize(n_clm);
        for(int i=0; i<n_clm; i++){
          temp_sh.vels[i] = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mass[i]);
        }

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

  DEBUG(10, "RDC RESTRAINTS")
  vector<topology::rdc_restraint_struct>  tmp_rdc_rest_strct;
  { // RDCRESSPEC
    DEBUG(10, "RDCRESSPEC")
    buffer = m_block["RDCRESSPEC"];

    if (!buffer.size()){
      io::messages.add("no RDCRESSPEC block in RDC restraints file", "In_RDC", io::message::critical);
    }

    DEBUG(10, "reading in RDCRESSPEC data")

    int atom_i = 0, atom_j = 0, k = 0, l = 0, type = 0;
    double weight = 0.0, R = 0.0, gyri = 0.0, gyrj = 0.0, rij = 0.0, rik = 0.0;

    os.setf(ios_base::fixed, ios_base::floatfield);

    DEBUG(10, setw(6) << "i"
           << setw(6) << "j"
           << setw(8) << "weight"
           << setw(19) << "R [1/ps]"
      //     << setw(19) << "RDCav [1/ps]"
           << setw(12) << "GYRi [e/u]"
           << setw(12) << "GYRj [e/u]")

    vector<string>::const_iterator it = buffer.begin()+1, to = buffer.end()-1;

    for(int n=0; it != to; ++it, ++n){

      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> atom_i >> atom_j >> k >> l >> weight >> R >> gyri >> gyrj >> rij >> rik >> type;

      if(_lineStream.fail()){
           io::messages.add("bad line in RDCRESSPEC block", "In_RDC", io::message::error);
      }

      R *= factorfreq;
      gyri *= factorgyr;
      gyrj *= factorgyr;

      // check for sensible choice of RDC
      if(abs(-( math::eps0_i * math::h_bar * gyri * gyrj )/( pow(math::spd_l,2) * 4.0 * pow(math::Pi,2)) * 1000) < abs(R)){
        io::messages.add("The chosen RDC is larger in magnitude than RDC_max.  This is probably a mistake and may result in strange behaviour.",
         "In_RDC", io::message::warning);
      }

      tmp_rdc_rest_strct.push_back(topology::rdc_restraint_struct(atom_i-1, atom_j-1, weight, R, gyri, gyrj)); // in units of 1, 1, 1, 1/ps, e/u and e/u

      DEBUG(10,                    setw(6) << atom_i
                                << setw(6) << atom_j
            << setprecision(2)  << setw(8) << weight
            << setprecision(14) << setw(19) << R
                                //<< setw(19) << conf.special().rdc.av[n]
            << setprecision(4)  << setw(12) << gyri
                                << setw(12) << gyrj)
    }

    if(sim.param().rdc.mode == simulation::rdc_restr_inst ||
       sim.param().rdc.mode == simulation::rdc_restr_av ||
       sim.param().rdc.mode == simulation::rdc_restr_biq){
      // as we don't weight any RDCs, we set all weights to one
      // no further distinction is required then
      io::messages.add("No weighting selected: Setting all weights to 1.0", "In_RDC", io::message::notice);
      vector<topology::rdc_restraint_struct>::iterator
        topo_it = tmp_rdc_rest_strct.begin(),
        topo_to = tmp_rdc_rest_strct.end();
      for(; topo_it!=topo_to; ++topo_it){ topo_it->weight = 1.0; }
    }

    DEBUG(10, "END")

  } // RDCRESSPEC


  //////////////////////
  // RDCGROUPS BLOCK //
  //////////////////////
  vector<vector<unsigned int> > rdc_groups;
  { // RDCGROUPS
    DEBUG(10, "RDCGROUPS")
    buffer = m_block["RDCGROUPS"];

    if (!buffer.size()){
      io::messages.add("no RDCGROUPS block in RDC restraints file", "In_RDC", io::message::critical);
    }

    unsigned int int_buf = 0;

    os.setf(ios_base::fixed, ios_base::floatfield);

    vector<string>::const_iterator
      it = buffer.begin()+1,
      to = buffer.end()-1;
    for(; it!=to; ++it){

      _lineStream.clear();
      _lineStream.str(*it);

      rdc_groups.push_back(vector<unsigned int>());
      while(_lineStream >> int_buf){
        rdc_groups.back().push_back(int_buf);
      }
      // cannot check state of stream because it will always be eof/fail
    }

    DEBUG(10, "number of RDCs: " << tmp_rdc_rest_strct.size())
    DEBUG(10, "number of RDC groups: " << rdc_groups.size())

    //sort RDCs in each group
    for(unsigned int i=0; i< rdc_groups.size(); i++){
       std::sort(rdc_groups[i].begin(), rdc_groups[i].end());
       vector<unsigned int>::iterator it = std::unique(rdc_groups[i].begin(), rdc_groups[i].end());
       if(it != rdc_groups[i].end()){
         io::messages.add("Removing duplicate RDC from group.  It might be a typo to add an RDC to a group twice.", "In_RDC", io::message::warning);
       }
       rdc_groups[i].resize( std::distance(rdc_groups[i].begin(),it) );
    }

#ifdef DEBUG
    for(unsigned int i=0; i< rdc_groups.size(); i++){
      cout << "{";
      for(unsigned int j=0; j< rdc_groups[i].size(); j++){
         cout << rdc_groups[i][j];
         (j<rdc_groups[i].size()-1) ? cout << ", " : cout << "}" << endl;
      }
    }
#endif

    DEBUG(10, "END")
  } // RDCGROUPS


/////////////////////////////////////////////////
//                                             //
//       further checking for validity         //
//                                             //
/////////////////////////////////////////////////


  DEBUG(10, "checking validity of RDC groups ...")

// check if blocks are larger than 5 each
  for(vector<vector<unsigned int> >::iterator it = rdc_groups.begin(), to = rdc_groups.end(); it!=to; it++){
    if(it->size() < 5){
      io::messages.add("less than five RDCs in one block lead to undefined situations ...", "In_RDC", io::message::error);
    } else if(it->size() == 5){
      io::messages.add("exactly five RDCs in one block lead to *no* restraint ... you probably don't want this", "In_RDC", io::message::warning);
    }
  }

// check if non-existing rdcs are in a groups
  for(unsigned int i=0; i!=rdc_groups.size(); i++){
    for(unsigned int j=0; j!=rdc_groups[i].size(); j++){
      if(rdc_groups[i][j] > tmp_rdc_rest_strct.size()){
        DEBUG(10, "RDC #" << rdc_groups[i][j] << " is part of RDC group #" << i << " but does not exist.")
        io::messages.add("non-existing RDCs included in RDC groups", "In_RDC", io::message::critical);
      }
    }
  }

// count in how many rdc groups every rdc appears
  const int n_rdc = tmp_rdc_rest_strct.size();
  vector<int> occurrence_count(n_rdc, 0);
  for (unsigned int i=0; i<rdc_groups.size(); ++i) {
    for (unsigned int j=0; j<rdc_groups[i].size(); ++j) {
      occurrence_count[rdc_groups[i][j]-1]++;
    }
  }
#ifdef DEBUG
  cout << "{";
  for(unsigned int i=0; i<occurrence_count.size(); i++){
    cout << occurrence_count[i];
    (i<occurrence_count.size()-1) ? (cout << ", ") : (cout << "}" << endl);
  }
#endif
 
// check if rdc groups contain each rdc at least once
  bool ignored_rdcs_exist = false;
  for (unsigned int i=0; i<occurrence_count.size(); ++i) {
    if (occurrence_count[i] == 0) {
      ignored_rdcs_exist = true;
      DEBUG(10, "RDC #" << i << " is not part of any RDC group.")
    }
  }
  if(ignored_rdcs_exist){
    io::messages.add("One or more RDCs are not part of any RDC group and will be ignored entirely.  You probably don't want this.  In particular, if you run svd-fit to analyse the trajectory, *all* RDCs will be taken into account, even the ones that are not restrained.",
      "In_RDC", io::message::warning);
  }

  DEBUG(10, "setting RDC weights according to occurrence in rdc groups")
  for (unsigned int i=0; i<occurrence_count.size(); ++i) {
    if(occurrence_count[i] != 0) tmp_rdc_rest_strct[i].weight /= occurrence_count[i];
  }


  DEBUG(10, "RDC group checking done")


/////////////////////////////////////////////////
//                                             //
//  writing temp-variables to the right places //
//                                             //
/////////////////////////////////////////////////

  sim.param().rdc.delta *= factorfreq;
  DEBUG(12, "delta (flat-bottom-pot) set to " << sim.param().rdc.delta )


  DEBUG(10, "writing parsed data to topo.rdc_restraints() and conf.special().rdc ...")

  const int n_clm = 5,  n_ah = 5;

// write stuff to topo.rdc_restraints() which has the type std::vector<std::vector<topology::rdc_restraint_struct> >
// reshuffle temporary vector<topology::rdc_restraint_struct> into vector<vector<topology::rdc_restraint_struct> >
    topo.rdc_restraints().resize(rdc_groups.size());
    unsigned int i=0;
    vector<vector<unsigned int> >::iterator it = rdc_groups.begin(), to = rdc_groups.end();
    for(; it!=to; it++, i++){
      topo.rdc_restraints()[i].resize(rdc_groups[i].size());
      unsigned int j=0;
      vector<unsigned int>::iterator jt = it->begin(), tj = it->end();
      for(; jt!=tj; jt++, j++){
        topo.rdc_restraints()[i][j] = tmp_rdc_rest_strct[rdc_groups[i][j]-1]; // in units of 1, 1, 1, 1/ps, e/u and e/u
      }
    }


// write stuff to conf.special().rdc which has the type std::vector<configuration::Configuration::special_struct::rdc_struct>

   conf.special().rdc.resize(rdc_groups.size());
   for(unsigned int i=0; i!=rdc_groups.size(); i++){
     const int group_size = rdc_groups[i].size();

     conf.special().rdc[i].av.resize(group_size);
     conf.special().rdc[i].curr.resize(group_size);

     conf.special().rdc[i].factorFreq = factorfreq;
     conf.special().rdc[i].factorGyr = factorgyr;

     switch(sim.param().rdc.type){
       case simulation::rdc_mf: {
         conf.special().rdc[i].MFpoint.resize(group_size);
         conf.special().rdc[i].MFpointVel.resize(group_size);
         conf.special().rdc[i].MFpointMass.resize(group_size);
         conf.special().rdc[i].stochastic_integral_mf.resize(temp_mf.cart_coords.size());
         break;
       }
       case simulation::rdc_t: {
         conf.special().rdc[i].Tensor.resize(group_size);
         conf.special().rdc[i].TensorVel.resize(group_size);
         conf.special().rdc[i].TensorMass.resize(group_size);
         conf.special().rdc[i].stochastic_integral_t.resize(n_ah);
         break;
       }
       case simulation::rdc_sh: {
         conf.special().rdc[i].clm.resize(group_size);
         conf.special().rdc[i].clmMass.resize(group_size);
         conf.special().rdc[i].clmVel.resize(group_size);
         conf.special().rdc[i].stochastic_integral_sh.resize(n_clm);
         break;
       }
       default: assert(false);
     }

     for(unsigned int j=0; j!=rdc_groups[i].size(); j++){
       conf.special().rdc[i].av[j] = topo.rdc_restraints()[i][j].R0; // init history as experimental values

       switch(sim.param().rdc.type){
         case simulation::rdc_mf: {
           conf.special().rdc[i].MFpoint = temp_mf.cart_coords;
           conf.special().rdc[i].MFpointVel = temp_mf.cart_vels;
           conf.special().rdc[i].MFpointMass = temp_mf.masses;
           break;
         }
         case simulation::rdc_t: {
           conf.special().rdc[i].Tensor = temp_t.a;
           conf.special().rdc[i].TensorVel = temp_t.vels;
           conf.special().rdc[i].TensorMass = temp_t.masses;
           break;
         }
         case simulation::rdc_sh: {
           conf.special().rdc[i].clm = temp_sh.clm;
           conf.special().rdc[i].clmVel = temp_sh.vels;
           conf.special().rdc[i].clmMass = temp_sh.masses;
           break;
         }
         default: assert(false);
       }
     }
   }

  // destroy the rng that was only created for this file
  delete rng;

  DEBUG(10, "done")

} // io::In_RDC::read

