/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
#include <vector>
#include <iosfwd>
#include <math/random.h>

#include "../../io/instream.h"
#include "../../io/blockinput.h"
#include "in_rdc.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

using namespace std;

static std::set<std::string> block_read;

// three structs, local to this file
struct mf_struct
{
  math::VArray cart_coords;
  math::VArray cart_vels;
  vector<double> masses;
};

struct t_struct
{
  vector<double> a;
  vector<double> vels;
  vector<double> masses;
};

struct sh_struct
{
  vector<double> clm;
  vector<double> vels;
  vector<double> masses;
};
mf_struct temp_mf;
t_struct temp_t;
sh_struct temp_sh;
std::vector<topology::rdc_restraint_struct> tmp_rdc_rest_strct;
std::vector<std::vector<unsigned int>> rdc_groups;

double generate_boltzmann_velocities(math::RandomGenerator *rng, const double temp, const double mass)
{
  const double sd = sqrt(math::k_Boltzmann * temp / mass);
  rng->stddev(sd);
  return rng->get_gauss();
}

// FIXME
// the following validity test is highly inelegant and misses some invalid a,
// but at least all the rejected ones *are* invalid.
bool a_is_valid(double a1, double a2, double a3, double a4, double a5)
{

  a1 = 2.0 / 3.0 * (a1 + .5);
  a2 = 2.0 / 3.0 * (a2 + .5);
  a3 = 2.0 / 3.0 * a3;
  a4 = 2.0 / 3.0 * a4;
  a5 = 2.0 / 3.0 * a5;

  if (a1 < 0.0 || a1 > 1.5)
    return false;
  if (a2 < 0.0 || a2 > 1.5)
    return false;
  if (a3 < -0.5 || a3 > 0.5)
    return false;
  if (a4 < -0.5 || a4 > 0.5)
    return false;
  if (a5 < -0.5 || a5 > 0.5)
    return false;

  if (a3 > sqrt(a1) * sqrt(a2) || a3 < -sqrt(a1) * sqrt(a2))
    return false;

  if (pow(a1, 2) + pow(a2, 2) + 2 * pow(a3, 2) + 2 * pow(a4, 2) + 2 * pow(a5, 2) > 1.0)
    return false;

  return true;
}

math::VArray create_points_on_sphere(const unsigned int N)
{
  if (N < 5)
  {
    io::messages.add("Please choose a number of at least 5 magnetic field vectors, to describe general alignment", "In_RDC", io::message::critical);
  }
  math::VArray coordinates;
  switch (N)
  {
  case 5:
  { // FIXME, this is not a good distribution
    coordinates.push_back(math::Vec(0, 0, 1));
    coordinates.push_back(math::Vec(1, 0, 0));
    coordinates.push_back(math::Vec(-0.5, cos(M_PI / 6.0), 0));
    coordinates.push_back(math::Vec(-0.5, -cos(M_PI / 6.0), 0));
    coordinates.push_back(math::Vec(0, 0, -1));
    break;
  }
  case 6:
  {
    coordinates.push_back(math::Vec(1, 0, 0));
    coordinates.push_back(math::Vec(0, 1, 0));
    coordinates.push_back(math::Vec(0, 0, 1));
    coordinates.push_back(math::Vec(-1, 0, 0));
    coordinates.push_back(math::Vec(0, -1, 0));
    coordinates.push_back(math::Vec(0, 0, -1));
    break;
  }
  default:
  {
    const double s = 3.6 / sqrt(N);
    const double dz = 2.0 / N;
    double longitude = 0.0;
    double z = 1.0 - dz / 2.0;
    for (unsigned int i = 0; i < N; ++i)
    {
      double r = sqrt(1.0 - z * z);
      coordinates.push_back(math::Vec(cos(longitude) * r, sin(longitude) * r, z));
      DEBUG(10, cos(longitude) * r << ", " << sin(longitude) * r << ", " << z)
      z -= dz;
      longitude += s / r;
    }
  }
  }
  return coordinates;
}

void io::In_RDC::read(topology::Topology &topo,
                      configuration::Configuration &conf,
                      simulation::Simulation &sim,
                      ostream &os)
{
  DEBUG(7, "reading in an RDC restraints specification file")

  os << "RDC RESTRAINTS\n";

#ifndef NDEBUG
  // Create random number generator, set seed to 0 (to make sure it's reproducible)
  math::RandomGenerator *rng = math::RandomGenerator::create(sim.param(), "0");
#else
  // Create random number generator, choose a remotely random random-seed
  ostringstream ss;
  ss << time(NULL);
  math::RandomGenerator *rng = math::RandomGenerator::create(sim.param(), ss.str());
#endif // DEBUG

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

  read_CONVERSION(topo, sim, os);
  read_INPUTMODE(topo, sim, os);

  switch (sim.param().rdc.type)
  {
  case simulation::rdc_mf:
  {

    read_MAGFIELDC(topo, sim, os, rng);
    break;
  }
  case simulation::rdc_t:
  { // only parse the tensor block if we use it
    read_ALIGNT(topo, sim, os, rng);
    break;
  }

  case simulation::rdc_sh:
  { // only parse the spherical-harmonics block if we use it
    read_SPHERICALHARMONICS(topo, sim, os, rng);
    break;
  }
  default:
    io::messages.add("no valid method chosen", "In_RDC", io::message::critical);
  }

  DEBUG(10, "RDC RESTRAINTS")
  read_RDCRESSPEC(topo, sim, os);
  read_RDCMOLAXIS(topo, sim, os);

  read_RDCGROUPS(topo, sim, os);

  /////////////////////////////////////////////////
  //                                             //
  //       further checking for validity         //
  //                                             //
  /////////////////////////////////////////////////

  DEBUG(10, "checking validity of RDC groups ...")

  // check if blocks are larger than 5 each
  for (vector<vector<unsigned int>>::iterator it = rdc_groups.begin(), to = rdc_groups.end(); it != to; it++)
  {
    if (it->size() < 5)
    {
      io::messages.add("less than five RDCs in one block lead to undefined situations ...", "In_RDC", io::message::error);
    }
    else if (it->size() == 5)
    {
      io::messages.add("exactly five RDCs in one block lead to *no* restraint ... you probably don't want this", "In_RDC", io::message::warning);
    }
  }

  // check if non-existing rdcs are in a group
  for (unsigned int i = 0; i != rdc_groups.size(); i++)
  {
    for (unsigned int j = 0; j != rdc_groups[i].size(); j++)
    {
      if (rdc_groups[i][j] > tmp_rdc_rest_strct.size())
      {
        DEBUG(10, "RDC #" << rdc_groups[i][j] << " is part of RDC group #" << i << " but does not exist.")
        io::messages.add("non-existing RDCs included in RDC groups", "In_RDC", io::message::critical);
        return;
      }
    }
  }

  // count in how many rdc groups every rdc appears
  const int n_rdc = tmp_rdc_rest_strct.size();
  vector<int> occurrence_count(n_rdc, 0);
  for (unsigned int i = 0; i < rdc_groups.size(); ++i)
  {
    for (unsigned int j = 0; j < rdc_groups[i].size(); ++j)
    {
      occurrence_count[rdc_groups[i][j] - 1]++;
    }
  }
#ifndef NDEBUG
  cout << "{";
  for (unsigned int i = 0; i < occurrence_count.size(); i++)
  {
    cout << occurrence_count[i];
    (i < occurrence_count.size() - 1) ? (cout << ", ") : (cout << "}" << endl);
  }
#endif

  // check if rdc groups contain each rdc at least once
  bool ignored_rdcs_exist = false;
  for (unsigned int i = 0; i < occurrence_count.size(); ++i)
  {
    if (occurrence_count[i] == 0)
    {
      ignored_rdcs_exist = true;
      DEBUG(10, "RDC #" << i << " is not part of any RDC group.")
    }
  }
  if (ignored_rdcs_exist)
  {
    io::messages.add("One or more RDCs are not part of any RDC group and will be ignored entirely.  You probably don't want this.  In particular, if you run svd-fit to analyse the trajectory, *all* RDCs will be taken into account, even the ones that are not restrained.",
                     "In_RDC", io::message::warning);
  }

  DEBUG(10, "setting RDC weights according to occurrence in rdc groups")
  for (unsigned int i = 0; i < occurrence_count.size(); ++i)
  {
    if (occurrence_count[i] != 0)
      tmp_rdc_rest_strct[i].weight /= occurrence_count[i];
  }

  DEBUG(10, "RDC group checking done")

  /////////////////////////////////////////////////
  //                                             //
  //  writing temp-variables to the right places //
  //                                             //
  /////////////////////////////////////////////////

  DEBUG(10, "writing parsed data to topo.rdc_restraints() and conf.special().rdc ...")

  const int n_clm = 5, n_ah = 5;

  // write stuff to topo.rdc_restraints() which has the type std::vector<std::vector<topology::rdc_restraint_struct> >
  // reshuffle temporary vector<topology::rdc_restraint_struct> into vector<vector<topology::rdc_restraint_struct> >
  topo.rdc_restraints().resize(rdc_groups.size());
  unsigned int i = 0;
  vector<vector<unsigned int>>::iterator it = rdc_groups.begin(), to = rdc_groups.end();
  for (; it != to; it++, i++)
  {
    topo.rdc_restraints()[i].resize(rdc_groups[i].size());
    unsigned int j = 0;
    vector<unsigned int>::iterator jt = it->begin(), tj = it->end();
    for (; jt != tj; jt++, j++)
    {
      topo.rdc_restraints()[i][j] = tmp_rdc_rest_strct[rdc_groups[i][j] - 1]; // in units of 1, 1, 1, 1/ps, e/u and e/u
    }
  }

  // write stuff to conf.special().rdc_groups which has the type std::vector<configuration::Configuration::special_struct::rdc_struct>

  conf.special().rdc_groups.resize(rdc_groups.size());
  for (unsigned int i = 0; i != rdc_groups.size(); i++)
  {
    const int group_size = rdc_groups[i].size();

    conf.special().rdc_groups[i].av.resize(group_size);
    conf.special().rdc_groups[i].curr.resize(group_size);
    conf.special().rdc_groups[i].RDC_cumavg.resize(group_size, 0);
    conf.special().rdc_groups[i].num_averaged = 0;

    switch (sim.param().rdc.type)
    {
    case simulation::rdc_mf:
    {
      conf.special().rdc_groups[i].MFpoint.resize(group_size);
      conf.special().rdc_groups[i].MFpointVel.resize(group_size);
      conf.special().rdc_groups[i].MFpointMass.resize(group_size);
      conf.special().rdc_groups[i].stochastic_integral_mf.resize(temp_mf.cart_coords.size());
      break;
    }
    case simulation::rdc_t:
    {
      conf.special().rdc_groups[i].Tensor.resize(group_size);
      conf.special().rdc_groups[i].Tensor_av.resize(group_size, 0);
      conf.special().rdc_groups[i].TensorVel.resize(group_size);
      conf.special().rdc_groups[i].TensorMass.resize(group_size);
      conf.special().rdc_groups[i].stochastic_integral_t.resize(n_ah);
      break;
    }
    case simulation::rdc_sh:
    {
      conf.special().rdc_groups[i].clm.resize(group_size);
      conf.special().rdc_groups[i].clmMass.resize(group_size);
      conf.special().rdc_groups[i].clmVel.resize(group_size);
      conf.special().rdc_groups[i].stochastic_integral_sh.resize(n_clm);
      break;
    }
    default:
      assert(false);
    }

    for (unsigned int j = 0; j != group_size; j++)
    {
      conf.special().rdc_groups[i].av[j] = topo.rdc_restraints()[i][j].D0; // init history as experimental values

      switch (sim.param().rdc.type)
      {
      case simulation::rdc_mf:
      {
        conf.special().rdc_groups[i].MFpoint = temp_mf.cart_coords;
        conf.special().rdc_groups[i].MFpointVel = temp_mf.cart_vels;
        conf.special().rdc_groups[i].MFpointMass = temp_mf.masses;
        break;
      }
      case simulation::rdc_t:
      {
        conf.special().rdc_groups[i].Tensor = temp_t.a;
        conf.special().rdc_groups[i].TensorVel = temp_t.vels;
        conf.special().rdc_groups[i].TensorMass = temp_t.masses;
        break;
      }
      case simulation::rdc_sh:
      {
        conf.special().rdc_groups[i].clm = temp_sh.clm;
        conf.special().rdc_groups[i].clmVel = temp_sh.vels;
        conf.special().rdc_groups[i].clmMass = temp_sh.masses;
        break;
      }
      default:
        assert(false);
      }
    }
  }

  // destroy the rng that was only created for this file
  delete rng;

  DEBUG(10, "done")

  os << "END\n";

} // io::In_RDC::read

/**
 * @section conversion CONVERSION block
 * The CONVERSION block is read from the RDC restraint specification file.
 *
 * - Frequency unit conversion factor to ps-1.
 * - Gyr.magn. ratio unit conversion factor to (e/u).
 * @snippet snippets/snippets.cc CONVERSION
 */
void io::In_RDC::read_CONVERSION(topology::Topology &topo,
                                 simulation::Simulation &sim,
                                 std::ostream &os)
{ // CONVERSION

  DEBUG(10, "CONVERSION block");
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "CONVERSION\n";
  exampleblock << "# FACFREQ   frequency conversion factor from input units to ps-1\n";
  exampleblock << "#           typically input is Hz and the factor is 10^-12\n";
  exampleblock << "# FACGYR    rgyr conversion factor from input units to e/u\n";
  exampleblock << "#           typically input is 10^6*rad/T*s and the factor is 0.010375\n";
  exampleblock << "#  FACFREQ        FACGYR\n";
  exampleblock << "   0.000000000001 0.010364272\n";
  exampleblock << "END\n";

  std::string blockname = "CONVERSION";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], true) == 0)
  {
    block_read.insert(blockname);

    DEBUG(10, "reading in RDC CONVERSION block")
    block.get_next_parameter("FACFREQ", sim.param().rdc.factorFreq, ">0", "");
    block.get_next_parameter("FACGYR", sim.param().rdc.factorGyr, ">0", "");

    os.precision(12);
    os.setf(ios_base::fixed, ios_base::floatfield);

    DEBUG(10, scientific << setprecision(6) << setw(14) << "factorFreq: " << sim.param().rdc.factorFreq)
    DEBUG(10, scientific << setprecision(6) << setw(14) << "factorGyr: " << sim.param().rdc.factorGyr)

    block.get_final_messages();
  } // if block empty or not there
} // CONVERSION

/**
 * @section rdcinputmode INPUTMODE block
 * This block is to set a 'basic' input mode ('0') in which most settings are chosen for
 * the user or an expert mode ('1') in which more choices can be made by the
 * user
 * @snippet snippets/snippets.cc INPUTMODE
 */
void io::In_RDC::read_INPUTMODE(topology::Topology &topo,
                                simulation::Simulation &sim,
                                std::ostream &os)
{
  DEBUG(10, "INPUTMODE block")
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "INPUTMODE\n";
  exampleblock << "# IM   0,1  input mode\n";
  exampleblock << "#      0:  basic mode, most settings are chosen for the user\n";
  exampleblock << "#      1:  expert mode, more choices can be made\n";
  exampleblock << "#  IM\n";
  exampleblock << "    0\n";
  exampleblock << "END\n";
  std::string blockname = "INPUTMODE";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0)
  {
    block_read.insert(blockname);

    DEBUG(10, "reading in INPUTMODE data");
    block.get_next_parameter("IM", input_mode, "", "0,1");

    DEBUG(10, "chosen inputmode: " << input_mode)

    block.get_final_messages();
  } // if block empty or not there
} // INPUTMODE

/**
 * @section magfieldc MAGFIELDC block
 * The MAGFIELDC block is read from the RDC restraint specification file. You need
 * this section if you choose NTRDCT 0 (cartesian representation of magnetic field
 * vectors) in your gromos configuration file.
 * Depending on the value in the INPUTMODE block the basic or advanced mode
 * parameters are expected.
 * @snippet snippets/snippets.cc MAGFIELDC
 */
void io::In_RDC::read_MAGFIELDC(topology::Topology &topo,
                                simulation::Simulation &sim,
                                std::ostream &os,
                                math::RandomGenerator *rng)
{
  DEBUG(10, "MAGFIELDC block")
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "MAGFIELDC\n";
  exampleblock << "# NMF   number of magnetic field direction vectors\n";
  exampleblock << "# x,y,z     cartesian coordinates of the moveable pseudo-atoms representing\n";
  exampleblock << "#   the magnetic field vectors. They specify the direction of the magnetic\n";
  exampleblock << "#   field vector, as the other coordinate is always (0, 0, 0). The magnetic field \n";
  exampleblock << "#   vectors are normalised when read in.\n";
  exampleblock << "# mass   the mass of each moveable atom and can be used as a damping\n";
  exampleblock << "#        factor \n";
  exampleblock << "# -- basic mode:\n";
  exampleblock << "# NMF  mass\n";
  exampleblock << "    5   1.0\n";
  exampleblock << "# -- advanced mode:\n";
  exampleblock << "#      x      y      z    mass\n";
  exampleblock << "     1.0    0.0    0.0     1.0\n";
  exampleblock << "     0.0    1.0    0.0     1.0\n";
  exampleblock << "END\n";

  std::string blockname = "MAGFIELDC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], true) == 0)
  {
    block_read.insert(blockname);

    DEBUG(10, "reading in MAGFIELDC data");

    if (input_mode == 0)
    { // basic
      unsigned int n_mf;
      double mf_mass;
      block.get_next_parameter("NMF", n_mf, ">=5", "");
      block.get_next_parameter("mass", mf_mass, ">1e-5", "");

      os.setf(ios_base::fixed, ios_base::floatfield);
      DEBUG(10, setprecision(4) << setw(8) << "NMF" << setw(8) << "mass")
      DEBUG(10, setprecision(4) << setw(8) << n_mf << setw(8) << mf_mass)

      temp_mf.cart_coords = create_points_on_sphere(n_mf);

      for (unsigned int i = 0; i < n_mf; ++i)
      {
        const double randomvel1 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]
        const double randomvel2 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]
        const double randomvel3 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]

        temp_mf.cart_vels.push_back(math::Vec(randomvel1, randomvel2, randomvel3));
        temp_mf.masses.push_back(mf_mass);
      }
    }
    else if (input_mode == 1)
    { // advanced
      DEBUG(10, setw(13) << "x" << setw(13) << "vx" << setw(13) << "y" << setw(13) << "vy" << setw(13) << "z" << setw(13) << "vz")
      unsigned int num = block.numlines() - 2;
      for (unsigned int line_number = 0; line_number < num; ++line_number)
      {
        DEBUG(11, "\tnr " << line_number);
        double MFx = 0.0, MFy = 0.0, MFz = 0.0, mf_mass = 0.0; // in units of nm, nm, nm and u
        block.get_next_parameter("x", MFx, "", "");
        block.get_next_parameter("y", MFy, "", "");
        block.get_next_parameter("z", MFz, "", "");
        block.get_next_parameter("mass", mf_mass, ">1e-5", "");

        const double randomvel1 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]
        const double randomvel2 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]
        const double randomvel3 = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mf_mass); // [nm/ps]

        temp_mf.cart_coords.push_back(math::Vec(MFx, MFy, MFz).norm());
        temp_mf.cart_vels.push_back(math::Vec(randomvel1, randomvel2, randomvel3));
        temp_mf.masses.push_back(mf_mass);
      }
    }
    block.get_final_messages();
  } // if block empty or not there

} // MAGFIELDC

/**
 * @section alignt ALIGNT block
 * The ALIGNT block is read from the RDC restraint specification file. You need
 * this section if you choose NTRDCT 1 (tensor representation of the alignment of
 * the molecule in the magnetic field) in your gromos configuration file.
 * Depending on the value in the INPUTMODE block the basic or advanced mode
 * parameters are expected.
 * @snippet snippets/snippets.cc ALIGNT
 */
void io::In_RDC::read_ALIGNT(topology::Topology &topo,
                             simulation::Simulation &sim,
                             std::ostream &os,
                             math::RandomGenerator *rng)
{
  DEBUG(10, "ALIGNT block")
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "ALIGNT\n";
  exampleblock << "# -- basic mode:\n";
  exampleblock << "# mass   the pseudo mass of each moveable atom, can be used as a\n";
  exampleblock << "#        damping factor \n";
  exampleblock << "#  mass\n";
  exampleblock << "    1.0\n";
  exampleblock << "# -- advanced mode:\n";
  exampleblock << "# Axx, Ayy, Axy, Axz, Ayz     the five independent tensor\n";
  exampleblock << "#    components of the 3x3 alignment tensor representing the alignment of\n";
  exampleblock << "#    the molecule in the magnetic field.\n";
  exampleblock << "# mAxx, mAyy, mAxy, mAxz, mAyz     masses\n";
  exampleblock << "#    Axx  mAxx Ayy  mAyy Axy  mAxy Axz  mAxz Ayz  mAyz\n";
  exampleblock << "     0.5   1.0 0.5   1.0 0.5   1.0 0.5   1.0 0.5   1.0\n";
  exampleblock << "END\n";
  std::string blockname = "ALIGNT";
  Block block(blockname, exampleblock.str());

  double mAxx, mAyy, mAxy, mAxz, mAyz;
  double Axx = 0, vAxx, Ayy = 0, vAyy, Axy = 0, vAxy, Axz = 0, vAxz, Ayz = 0, vAyz;

  bool block_required = false;
  if (sim.param().rdc.type != simulation::rdc_t)
    block_required = true;
  if (block.read_buffer(m_block[blockname], block_required) == 0)
  {
    block_read.insert(blockname);

    DEBUG(10, "reading in ALIGNT data");

    if (input_mode == 0)
    { // simple input

      block.get_next_parameter("mass", mAxx, ">1e-5", "");
      mAyy = mAxy = mAxz = mAyz = mAxx;

      //        // set first c value to -2 (which is an invalid value) so signal to the
      //        // MD or SD routine, that correct values should be calculated.  In case
      //        // of EM it is overwritten anyway.
      //        conf.special().rdc.Tensor[0] = -2;
    }
    else if (input_mode == 1)
    { // advanced input
      os.setf(ios_base::fixed, ios_base::floatfield);

      block.get_next_parameter("Axx", Axx, "", "");
      block.get_next_parameter("mAxx", mAxx, ">1e-5", "");
      block.get_next_parameter("Ayy", Ayy, "", "");
      block.get_next_parameter("mAyy", mAyy, ">1e-5", "");
      block.get_next_parameter("Axy", Axy, "", "");
      block.get_next_parameter("mAxy", mAxy, ">1e-5", "");
      block.get_next_parameter("Axz", Axz, "", "");
      block.get_next_parameter("mAxz", mAxz, ">1e-5", "");
      block.get_next_parameter("Ayz", Ayz, "", "");
      block.get_next_parameter("mAyz", mAyz, ">1e-5", "");

      // the tensor components a_1 to a_5 can be expressed as averages over
      // products of cosines and are, hence, limited to a range of values
      if (!a_is_valid(Axx, Ayy, Axy, Axz, Ayz))
      {
        io::messages.add("some or all of Axx, Ayy, Axy, Axz, Ayz have insensible values.  Setting all components to 0.0 (i.e. isotropic alignment) ...", "In_RDC", io::message::error);
        Axx = 0.0;
        Ayy = 0.0;
        Axy = 0.0;
        Axz = 0.0;
        Ayz = 0.0;
      }
    }
    // Generate Boltzmann distributed velocities
    vAxx = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAxx);
    vAyy = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAyy);
    vAxy = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAxy);
    vAxz = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAxz);
    vAyz = generate_boltzmann_velocities(rng, sim.param().rdc.temp, mAyz);

    os.setf(ios_base::fixed, ios_base::floatfield);
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
    block.get_final_messages();
  } // if block empty or not there
} // ALIGNT

/**
 * @section sphericalharmonics SPHERICALHARMONICS block
 * The SPHERICALHARMONICS block is read from the RDC restraint specification file. You need
 * this section if you choose NTRDCT 2 in your gromos configuration file.
 * Depending on the value in the INPUTMODE block the basic or advanced mode
 * parameters are expected.
 * @snippet snippets/snippets.cc SPHERICALHARMONICS
 */
void io::In_RDC::read_SPHERICALHARMONICS(topology::Topology &topo,
                                         simulation::Simulation &sim,
                                         std::ostream &os,
                                         math::RandomGenerator *rng)
{
  DEBUG(10, "SPHERICALHARMONICS block")
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "SPHERICALHARMONICS\n";
  exampleblock << "# -- basic mode:\n";
  exampleblock << "# mass   the pseudo mass, can be used as a\n";
  exampleblock << "#        damping factor \n";
  exampleblock << "#  mass\n";
  exampleblock << "    1.0\n";
  exampleblock << "# -- advanced mode:\n";
  exampleblock << "# c1-5     Spherical harmonics coefficients\n";
  exampleblock << "# mass1-5     Spherical harmonics coefficient masses\n";
  exampleblock << "#    c1    mass1   c2    mass2 c3    mass3   c4    mass4  c5   mass5\n";
  exampleblock << "     0.5   1.0     0.5   1.0 0.5   1.0 0.5   1.0 0.5   1.0\n";
  exampleblock << "END\n";
  std::string blockname = "SPHERICALHARMONICS";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], true) == 1)
  {
    block_read.insert(blockname);

    DEBUG(10, "reading in SPHERICALHARMONICS data");

    const int n_clm = 5; // (-2,2), (-1,2), (0,2), (1,2), (2,2)
    if (input_mode == 0)
    { // basic input

      double mass;
      block.get_next_parameter("mass", mass, ">1e-5", "");

      //        // set first c value to -2 (which is an invalid value) so signal to the
      //        // MD or SD routine, that correct values should be calculated.  In case
      //        // of EM it is overwritten anyway.
      //        conf.special().rdc.clm[0] = -2;

      temp_sh.clm.resize(n_clm, 0.0);
      temp_sh.masses.resize(n_clm, mass);
    }
    else if (input_mode == 1)
    { // advanced input

      vector<double> c(5, 0.0), mass(5, 0.0);

      block.get_next_parameter("c1", c[0], "", "");
      block.get_next_parameter("mass1", mass[0], ">1e-5", "");
      block.get_next_parameter("c2", c[1], "", "");
      block.get_next_parameter("mass2", mass[1], ">1e-5", "");
      block.get_next_parameter("c3", c[2], "", "");
      block.get_next_parameter("mass3", mass[2], ">1e-5", "");
      block.get_next_parameter("c4", c[3], "", "");
      block.get_next_parameter("mass4", mass[3], ">1e-5", "");
      block.get_next_parameter("c5", c[4], "", "");
      block.get_next_parameter("mass5", mass[4], ">1e-5", "");

      temp_sh.clm = c;
      temp_sh.masses = mass;
    }
    temp_sh.vels.resize(n_clm);
    for (int i = 0; i < n_clm; i++)
    {
      temp_sh.vels[i] = generate_boltzmann_velocities(rng, sim.param().rdc.temp, temp_sh.masses[i]);
    }
    block.get_final_messages();
  } // if block empty or not there
} // SPHERICALHARMONICS

/**
 * @section rdcresspec RDCRESSPEC block
 * The RDCRESSPEC block is read from the RDC restraint specification file.
 * @snippet snippets/snippets.cc RDCRESSPEC
 */
void io::In_RDC::read_RDCRESSPEC(topology::Topology &topo,
                                 simulation::Simulation &sim,
                                 std::ostream &os)
{ // RDCRESSPEC
  DEBUG(10, "RDCRESSPEC")
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "RDCRESSPEC\n";
  exampleblock << "# DISH, DISC carbon-hydrogen/carbon-carbon distance\n";
  exampleblock << "# i,j,k,l  atoms comprising the virtual atom (put 0 if less than four atoms in use)\n";
  exampleblock << "# type   virtual atom type\n";
  exampleblock << "# R0                inter-nuclear distance\n";
  exampleblock << "# G1, G2            gyromagnetic ratios\n";
  exampleblock << "# D0, WRDC  target RDC and force constant weighting factor\n";
  exampleblock << "# DD0 half the width of the flatbottom potential\n";
  exampleblock << "# DISH  DISC\n";
  exampleblock << "  0.1   0.153\n";
  exampleblock << "#  i   j   k   l   type  i    j    k   l  type  R0      G1      G2      D0      DD0   WRDC\n";
  exampleblock << "  16   0   0   0   0    17    0    0   0  0    0.104   -27.12   267.52  -1.78   1.0    1\n";
  exampleblock << "  24   0   0   0   0    25    0    0   0  0    0.104   -27.12   267.52  8.78   1.0    1\n";
  exampleblock << "END\n";

  std::string blockname = "RDCRESSPEC";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], true) == 0)
  {
    block_read.insert(blockname);

    block.get_next_parameter("DISH", dish, ">0", "");
    block.get_next_parameter("DISC", disc, ">0", "");

    DEBUG(10, "reading in RDCRESSPEC data")

    os.setf(ios_base::fixed, ios_base::floatfield);

    DEBUG(10, setw(6) << "i"
                      << setw(6) << "j"
                      << setw(6) << "k"
                      << setw(6) << "l"
                      << setw(6) << "type"
                      << setw(6) << "i"
                      << setw(6) << "j"
                      << setw(6) << "k"
                      << setw(6) << "l"
                      << setw(6) << "type"
                      << setw(19) << "R0"
                      << setw(12) << "G1 [e/u]"
                      << setw(12) << "G2 [e/u]"
                      << setw(8) << "D0 [1/ps]"
                      << setw(8) << "DD0 [1/ps]"
                      << setw(8) << "WRDC");

    unsigned int num = block.numlines() - 3;
    for (unsigned int line_number = 0; line_number < num; line_number++)
    {
      std::vector<int> atom1, atom2;
      int type1, type2;
      double weight, D, DD, gyr1, gyr2, r0;

      DEBUG(11, "\trdc line " << line_number);

      for (unsigned int i = 0; i < 4; i++)
      {
        unsigned int atom;
        std::string str_i = io::to_string(i);
        block.get_next_parameter("ATOM[" + str_i + "]", atom, ">=0", "");
        DEBUG(11, "\tat " << i << " " << atom);
        if (atom > topo.num_atoms())
        {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is " << topo.num_atoms();
          io::messages.add(msg.str(), "In_RDC", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0)
          atom1.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type1, "", "-2,-1,0,1,2,3,4,5,6,7");

      for (unsigned int i = 0; i < 4; i++)
      {
        unsigned int atom;
        std::string str_i = io::to_string(i);
        block.get_next_parameter("ATOM[" + str_i + "]", atom, ">=0", "");
        if (atom > topo.num_atoms())
        {
          std::ostringstream msg;
          msg << blockname << " block: atom number out of range: " << atom << ", last atom is " << topo.num_atoms();
          io::messages.add(msg.str(), "In_RDC", io::message::error);
        }

        // -1 because we directly convert to array indices
        if (atom > 0)
          atom2.push_back(atom - 1);
      }
      block.get_next_parameter("TYPE", type2, "", "-2,-1,0,1,2,3,4,5,6,7");

      block.get_next_parameter("R0", r0, ">0", "");
      block.get_next_parameter("GYR1", gyr1, "", "");
      block.get_next_parameter("GYR2", gyr2, "", "");
      block.get_next_parameter("D0", D, "", "");
      block.get_next_parameter("DD0", DD, ">=0", "");
      block.get_next_parameter("W0", weight, ">=0", "");

      if (!block.error())
      {

        D *= sim.param().rdc.factorFreq;
        DD *= sim.param().rdc.factorFreq;
        gyr1 *= sim.param().rdc.factorGyr;
        gyr2 *= sim.param().rdc.factorGyr;

        // check for sensible choice of RDC
        if (abs(-(math::eps0_i * math::h_bar * gyr1 * gyr2) / (pow(math::spd_l, 2) * 4.0 * pow(math::Pi, 2)) * 1000) < abs(D))
        {
          io::messages.add("The chosen RDC is larger in magnitude than RDC_max.  This is probably a mistake and may result in strange behaviour.",
                           "In_RDC", io::message::warning);
        }

        util::virtual_type t1 = util::virtual_type(type1);
        util::virtual_type t2 = util::virtual_type(type2);

        util::Virtual_Atom v1(t1, atom1, dish, disc);
        util::Virtual_Atom v2(t2, atom2, dish, disc);

        tmp_rdc_rest_strct.push_back(topology::rdc_restraint_struct(v1, v2, weight, D, DD, gyr1, gyr2, r0)); // in units of 1, 1, 1, 1/ps, e/u and e/u

        DEBUG(10, setw(6) << atom1[0]
                          << setw(6) << atom1[1]
                          << setw(6) << atom1[2]
                          << setw(6) << atom1[3]
                          << setw(6) << type1
                          << setw(6) << atom2[0]
                          << setw(6) << atom2[1]
                          << setw(6) << atom2[2]
                          << setw(6) << atom2[3]
                          << setw(6) << type2
                          << setw(8) << r0
                          << setprecision(4) << setw(12) << gyr1
                          << setw(12) << gyr2
                          << setprecision(14) << setw(19) << D
                          << setprecision(14) << setw(19) << DD
                          << setprecision(2) << setw(8) << weight);
      }
    } // for restraint-lines
    block.get_final_messages();

    if (sim.param().rdc.mode == simulation::rdc_restr_inst ||
        sim.param().rdc.mode == simulation::rdc_restr_av ||
        sim.param().rdc.mode == simulation::rdc_restr_biq)
    {
      // as we don't weight any RDCs, we set all weights to one
      // no further distinction is required then
      io::messages.add("No weighting selected: Setting all weights to 1.0", "In_RDC", io::message::notice);
      vector<topology::rdc_restraint_struct>::iterator
          topo_it = tmp_rdc_rest_strct.begin(),
          topo_to = tmp_rdc_rest_strct.end();
      for (; topo_it != topo_to; ++topo_it)
      {
        topo_it->weight = 1.0;
      }
    }
  } // if block content
} // RDCRESSPEC

/**

 * @section rdcmolaxis RDCMOLAXIS block
 * The RDCMOLAXIS block is read from the RDC restraint specification file.
 * It defines a molecular axis, the angle of which with the magn. field vector is calculated to write out its distribution at the end of a run.
 * @snippet snippets/snippets.cc RDCMOLAXIS
 */
void io::In_RDC::read_RDCMOLAXIS(topology::Topology &topo,
                                 simulation::Simulation &sim,
                                 std::ostream &os)
{ // RDCMOLAXIS
  DEBUG(10, "RDCMOLAXIS block")
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "RDCMOLAXIS\n";
  exampleblock << "# i,j,k,l  atoms comprising the virtual atom (put 0 if less than four atoms in use)\n";
  exampleblock << "# type   virtual atom type \n";
  exampleblock << "# i   j   k   l   type  i    j    k   l  type\n";
  exampleblock << "  16  0   0   0    0    17   0    0   0  0\n";
  exampleblock << "END\n";
  std::string blockname = "RDCMOLAXIS";
  Block block(blockname, exampleblock.str());

  bool block_required = false;
  if (sim.param().rdc.type == simulation::rdc_t)
    block_required = true;

  if (block.read_buffer(m_block[blockname], block_required) == 0)
  {
    if (sim.param().rdc.type != simulation::rdc_t) {
        io::messages.add("Ignoring RDCMOLAXIS block, it is only used with the alignm. tensor representation ", "In_RDC", io::message::warning);
        return;
    }
    block_read.insert(blockname);

    DEBUG(10, "reading in RDCMOLAXIS data");

    int type1, type2;
    std::vector<int> atom1, atom2;
    // double dist_norm, gyri, gyrj, D0, dD0, w;

    for (unsigned int i = 0; i < io::In_RDC::MAX_ATOMS; i++)
    {
      unsigned int atom;
      std::string str_i = io::to_string(i);
      std::string condition = i == 0 ? ">0" : ">=0";
      block.get_next_parameter("ATOM[" + str_i + "]", atom, condition, "");
      if (atom > topo.num_atoms())
      {
        std::ostringstream msg;
        msg << blockname << " block: atom number out of range: " << atom << ", last atom is " << topo.num_atoms();
        io::messages.add(msg.str(), "In_RDC", io::message::error);
      }

      // -1 because we directly convert to array indices
      if (atom > 0)
        atom1.push_back(atom - 1);
    }
    block.get_next_parameter("type1", type1, "", "");

    for (unsigned int i = 0; i < io::In_RDC::MAX_ATOMS; i++)
    {
      unsigned int atom;
      std::string str_i = io::to_string(i);
      std::string condition = i == 0 ? ">0" : ">=0";
      block.get_next_parameter("ATOM[" + str_i + "]", atom, condition, "");
      if (atom > topo.num_atoms())
      {
        std::ostringstream msg;
        msg << blockname << " block: atom number out of range: " << atom << ", last atom is " << topo.num_atoms();
        io::messages.add(msg.str(), "In_RDC", io::message::error);
      }

      // -1 because we directly convert to array indices
      if (atom > 0)
        atom2.push_back(atom - 1);
    }
    block.get_next_parameter("type2", type2, "", "");

    util::virtual_type t1 = util::virtual_type(type1);
    util::virtual_type t2 = util::virtual_type(type2);

    util::Virtual_Atom v1(t1, atom1, dish, disc);
    util::Virtual_Atom v2(t2, atom2, dish, disc);

    topo.rdc_molaxis().push_back(v1);
    topo.rdc_molaxis().push_back(v2);
    block.get_final_messages();
  } // if block empty or not there
} // RDCMOLAXIS

/**
 * @section rdcgroups RDCGROUPS block
 * The RDCGROUPS block is read from the RDC restraint specification file.
 * Each line contains the indices of RDCs belonging to a group.
 * RDCs in one group share the same alignment, their atoms to not move much with respect to each other.
 * @snippet snippets/snippets.cc RDCGROUPS
 */
void io::In_RDC::read_RDCGROUPS(topology::Topology &topo,
                                simulation::Simulation &sim,
                                std::ostream &os)
{ // RDCGROUPS
  DEBUG(10, "RDCGROUPS block")
  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line is the tag
  exampleblock << "RDCGROUPS\n";
  exampleblock << "#  RDC indices\n";
  exampleblock << "    1 2 3 4 5\n";
  exampleblock << "    6 7 8\n";
  exampleblock << "END\n";
  std::string blockname = "RDCGROUPS";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], true) == 0)
  {
    block_read.insert(blockname);

    DEBUG(10, "reading in RDCGROUPS data");

    os.setf(ios_base::fixed, ios_base::floatfield);
    block.get_as_array("RG", rdc_groups, "", "");
    block.get_final_messages();

    os << "  "
       << "number of RDCs: " << tmp_rdc_rest_strct.size() << endl;
    os << "  "
       << "number of RDC groups: " << rdc_groups.size() << endl;
    DEBUG(10, "number of RDCs: " << tmp_rdc_rest_strct.size())
    DEBUG(10, "number of RDC groups: " << rdc_groups.size())

    // sort RDCs in each group
    for (unsigned int i = 0; i < rdc_groups.size(); i++)
    {
      std::sort(rdc_groups[i].begin(), rdc_groups[i].end());
      vector<unsigned int>::iterator it = std::unique(rdc_groups[i].begin(), rdc_groups[i].end());
      if (it != rdc_groups[i].end())
      {
        io::messages.add("Removing duplicate RDC from group.  It might be a typo to add an RDC to a group twice.", "In_RDC", io::message::warning);
      }
      rdc_groups[i].resize(std::distance(rdc_groups[i].begin(), it));
    }

#ifndef NDEBUG
    for (unsigned int i = 0; i < rdc_groups.size(); i++)
    {
      cout << "{";
      for (unsigned int j = 0; j < rdc_groups[i].size(); j++)
      {
        cout << rdc_groups[i][j];
        (j < rdc_groups[i].size() - 1) ? cout << ", " : cout << "}" << endl;
      }
    }
#endif
  }
} // RDCGROUPS