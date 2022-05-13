/**
 * @file in_topology.cc
 * implements methods of In_Topology.
 */


#include "../../stdheader.h"

#include "../../topology/topology.h"
#include "../../simulation/multibath.h"
#include "../../simulation/parameter.h"
#include "../../interaction/interaction_types.h"
#include "../../io/instream.h"
#include "../../util/parse_tcouple.h"

#include "../../io/blockinput.h"
#include <vector>

#include "in_topology.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

template<typename T>
bool check_type(std::vector<std::string> const & buffer, std::vector<T> term) {
  if (buffer.size()) {
    std::istringstream is(buffer[1]);
    int num = 0;
    if (!(is >> num) || num < 0) {
      std::cout << "ERROR:\tcould not read num or smaller 0!"
              << std::endl;
      return false;
    }

    for (typename std::vector<T>::const_iterator
      it = term.begin(),
            to = term.end();
            it != to;
            ++it) {

      if (int(it->type) >= num) {
        std::cout << "ERROR:\tused type " << it->type + 1
                << " larger than max type (" << num << ")"
                << std::endl;

        return false;
      }
    }
    return true;
  }
  if (term.size() == 0) return true;

  std::cout << "ERROR:\tblock not found!" << std::endl;
  return false;
}

template<typename T>
bool check_type(int num, std::vector<T> term) {

    for (typename std::vector<T>::const_iterator
      it = term.begin(),
            to = term.end();
            it != to;
            ++it) {

      if (int(it->type) >= num) {
        std::cout << "ERROR:\tused type " << it->type + 1
                << " larger than max type (" << num << ")"
                << std::endl;

        return false;
      }
    }
    return true;
}

void io::In_Topology::read(topology::Topology& topo,
        simulation::Parameter &param,
        std::ostream & os) {

  DEBUG(7, "reading in topology");

  if (!quiet) {
    os << "TOPOLOGY\n";
    os << title << "\n";
  }

  block_read.clear();

  if (m_block["TOPVERSION"].size()) {
    block_read.insert("TOPVERSION");
  }

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  // allow only BONDSTRETCHTYPE or one or both of the other two
  if (m_block.count("HARMBONDTYPE") + m_block.count("BONDTYPE") != 0 && m_block.count("BONDSTRETCHTYPE") != 0)
          io::messages.add("Topology should contain either a BONDSTRETCHTYPE block or BONDTYPE/HARMBONDTYPE",
                  "In_Topology", io::message::error);

  // allow only BONDANGLEBENDTYPE or one or both of the other two
  if (m_block.count("HARMBONDANGLETYPE") + m_block.count("BONDANGLETYPE") != 0 && m_block.count("BONDANGLEBENDTYPE") != 0)
          io::messages.add("Topology should contain either a BONDANGLEBENDTYPE block or BONDANGLETYPE/HARMBONDANGLETYPE",
                  "In_Topology", io::message::error);

  read_bond_types(topo, param, os);
  read_bondangle_types(topo, param, os);
  read_dihedral_types(topo, param, os);

  read_block_IMPDIHEDRALTYPE(topo, param, os);

  read_block_TYPE(topo, param, os);
  read_block_PHYSICALCONSTANTS(topo, param, os);

  if (param.system.npm) {

  read_block_RESNAME(topo, param, os);
  read_block_ATOMTYPENAME(topo, param, os);
  read_block_SOLUTEATOM(topo, param, os);

  // os << "time after SOLUTEATOM: " << util::now() - start << std::endl;

    { //scale mass for adiabatic decoupling
      if (param.addecouple.adgr > 0) {
        double sm = 1;
        int adc_index = 0;
        int num = topo.num_solute_atoms();
        for (int n = 0; n < num; ++n) {
          sm = 1;
          adc_index = param.addecouple.check_index_adc(n);
          if (adc_index != -1) {
            sm = param.addecouple.adc_index()[adc_index].sm;
          }
          topo.mass()(n) *= sm;
        }
      }
    }//adiabatic decoupling

    read_block_SOLUTEPOLARISATION(topo, param, os);
    read_block_CGSOLUTE(topo, param, os);
    read_block_LJEXCEPTIONS(topo, param, os);
    read_block_BONDH(topo, param, os);
    read_block_BOND(topo, param, os);
    read_block_BONDDP(topo, param, os);

    { // check the bonds
      if (param.force.bond) {
        if (!check_type(num_solute_bondtypes, topo.solute().bonds()) ||
                !check_type(num_solute_bondtypes, topo.solute().cgbonds())) {
          io::messages.add("Illegal bond type in BOND(H) or BONDDP block"
                           "(type not defined in topology file)",
                  "In_Topology", io::message::error);
        }
      }
    }

    read_block_CONSTRAINT(topo, param, os);

    // check the bonds in constraints
    if (!check_type(num_solute_bondtypes, topo.solute().distance_constraints())) {
      io::messages.add("Illegal bond type in CONSTRAINT block",
              "In_Topology", io::message::error);
    }
    read_block_BONDANGLEH(topo, param, os);
    read_block_BONDANGLE(topo, param, os);

    // os << "time after BONDANGLE: " << util::now() - start << std::endl;

    // check the angles
    if (param.force.angle) {
      if (!check_type(num_solute_angletypes, topo.solute().angles())) {
        io::messages.add("Illegal bond angle type in BONDANGLE(H) block",
              "In_Topology", io::message::error);
      }
    }

    read_block_IMPDIHEDRAL(topo, param, os);
    read_block_IMPDIHEDRALH(topo, param, os);

    // check the impropers
    if (!check_type(num_solute_impropertypes, topo.solute().improper_dihedrals())) {
      io::messages.add("Illegal improper dihedral type in IMPDIHEDRAL(H) block",
              "In_Topology", io::message::error);
    }

    read_block_DIHEDRAL(topo, param, os);
    read_block_DIHEDRALH(topo, param, os);

    // check the dihedrals
    if (!check_type(num_solute_dihedraltypes, topo.solute().dihedrals())) {
      io::messages.add("Illegal dihedral type in DIHEDRAL(H) block",
              "In_Topology", io::message::error);
    }

    read_block_CROSSDIHEDRAL(topo, param, os);
    read_block_CROSSDIHEDRALH(topo, param, os);

    // check the crossdihedrals
    if (m_block["CROSSDIHEDRALH"].size() && !m_block["TORSDIHEDRALTYPE"].size()) {
      io::messages.add("TORSDIHEDRALTYPE block must be specified to use crossdihedrals",
              "In_Topology", io::message::error);
    }
    
    if (!check_type(num_solute_dihedraltypes, topo.solute().crossdihedrals())) {
      io::messages.add("Illegal crossdihedral type in CROSSDIHEDRAL(H) block",
              "In_Topology", io::message::error);
    }

    read_block_VIRTUALGRAIN(topo, param, os);

    // add the solute molecules (should be done before solvate ;-)

    read_block_SOLUTEMOLECULES(topo, param, os);

    if (topo.molecules().size() == 0) {
      topo.molecules().push_back(0);
      topo.molecules().push_back(topo.num_solute_atoms());
    }

    // solutemolecules check
    if (topo.molecules().back()
            != topo.num_solute_atoms()) {

      std::cout << "ERROR: Solute molecules wrong\n"
              << "\tlast atom in SOLUTEMOLECULES block = " << topo.molecules().back()
              << "\n\tlast atom in topology = " << topo.num_solute_atoms()
              << std::endl;

      io::messages.add("Error in SOLUTEMOLECULES block: "
              "last solute molecule has to end with last solute atom",
              "In_Topology", io::message::error);
    }

    topo.num_solute_molecules() = topo.molecules().size() - 1;


    read_block_TEMPERATUREGROUPS(topo, param, os);

    if (topo.temperature_groups().size() == 0) {
      topo.temperature_groups().push_back(0);
      topo.temperature_groups().push_back(topo.num_solute_atoms());
    }

    // temperature groups check
    if (topo.temperature_groups().back()
            != topo.num_solute_atoms()) {

      std::cout << "ERROR: temperature groups wrong\n"
              << "\tlast atom in TEMPERATUREGROUPS block = " << topo.temperature_groups().back()
              << "\n\tlast atom in topology = " << topo.num_solute_atoms()
              << std::endl;

      io::messages.add("Error in TEMPERATUREGROUPS block: "
              "last temperature group has to end with last solute atom",
              "In_Topology", io::message::error);
    }

    topo.num_solute_temperature_groups() = topo.temperature_groups().size() - 1;


    read_block_PRESSUREGROUPS(topo, param, os);

    if (topo.pressure_groups().size() == 0) {
      topo.pressure_groups().push_back(0);
      topo.pressure_groups().push_back(topo.num_solute_atoms());
    }

    // pressure groups check
    if (topo.pressure_groups().back()
            != topo.num_solute_atoms()) {

      std::cout << "ERROR: pressure groups wrong\n"
              << "\tlast atom in PRESSUREGROUPS block = " << topo.pressure_groups().back()
              << "\n\tlast atom in topology = " << topo.num_solute_atoms()
              << std::endl;

      io::messages.add("Error in PRESSUREGROUPS block: "
              "last pressure group has to end with last solute atom",
              "In_Topology", io::message::error);
    }

    topo.num_solute_pressure_groups() = topo.pressure_groups().size() - 1;

    { // PATHINTSPEC
      buffer = m_block["PATHINTSPEC"];
      if (buffer.size()) {
        block_read.insert("PATHINTSPEC");
        io::messages.add("md++ does not support path-integral simulations (PATHINTSPEC block).",
              "InTopology", io::message::warning);
      }
    }
  } // npm != 0

  { // SOLVENTATOM and SOLVENTCONSTR
    // give it a number (SOLVENTATOM1, SOLVENTATOM2) for multiple
    // solvents...
      read_SOLVENT_blocks(topo, param, os);
      // solvent atoms have been read into s
  }

  // add the solvent to the topology
  if (!quiet)
    os << "\n\t\tadding " << param.system.nsm
          << " solvents.";

  if (param.system.nsm)
    topo.solvate(0, param.system.nsm);

  if (!quiet)
    os << "\n\tEND\n";

  {// SASAPARAMETER
    if (param.sasa.switch_sasa) { // if SASA is switched on

      read_sasa_parameter(topo, topo.sasa_parameter());

      if (!quiet) {
        os << "\tSASAPARAMETERS\n\t\tSASA is switched on\n\t\tN_sasa_atoms :\t"
        << topo.sasa_parameter().size() << "\n\t\tp_12 :\t\t" << param.sasa.p_12
        << "\n\t\tp_13 :\t\t" << param.sasa.p_13 << "\n\t\tp_1x :\t\t" << param.sasa.p_1x
        << "\n\t\tR_solv :\t" << param.sasa.r_solv << std::endl;

        if (param.sasa.switch_volume) {
          os << "\t\tVOL is switched on\n\t\tsigma_vol :\t" << param.sasa.sigma_v << std::endl;
        }
        os << "\tEND\n";
      }
    }
    // SASAPARAMETER
  }

  // set lambda (new and old one, yes it looks strange...)
  topo.lambda(param.perturbation.lambda);
  topo.lambda(param.perturbation.lambda);

  topo.lambda_exp(param.perturbation.lambda_exponent);

  //==================================================
  // CHECKING
  //==================================================

  // solute molecule check
  if (topo.molecules().back()
          != topo.num_atoms()) {

    io::messages.add("Error in SOLUTEMOLECULE / solvation block: "
            "last solute molecule has to end with last atom",
            "In_Topology", io::message::error);
  }

  // temperature group check
  if (topo.temperature_groups().back()
          != topo.num_atoms()) {

    io::messages.add("Error in TEMPERATUREGROUPS / solvation block: "
            "last temperature group has to end with last atom",
            "In_Topology", io::message::error);
  }

  // pressure group check
  if (topo.pressure_groups().back()
          != topo.num_atoms()) {

    io::messages.add("Error in PRESSUREGROUPS / solvation block: "
            "last pressure group has to end with last atom",
            "In_Topology", io::message::error);
  }

  // chargegroup check (starts with 0)
  if (topo.chargegroups()[topo.num_solute_chargegroups()] != int(topo.num_solute_atoms())) {
    io::messages.add("Error: last solute atom has to be end of chargegroup",
            "In_Topology",
            io::message::error);
    os << "ERROR:"
            << "\tsolute cg    : " << topo.num_solute_chargegroups() << "\n"
            << "\tsolute atoms : " << topo.num_solute_atoms() << "\n"
            << "\tlast cg      : " << topo.chargegroups()[topo.num_solute_chargegroups()] << "\n";
  }

  if (!quiet)
    os << "\tSOLUTE [sub]molecules: "
          << unsigned(topo.molecules().size()) - param.system.nsm - 1 << "\n";

  DEBUG(10, "molecules().size: " << unsigned(topo.molecules().size())
          << " nsm : " << param.system.nsm);

  if (!quiet)
    os << "\tSOLUTE temperature groups: "
          << unsigned(topo.temperature_groups().size()) - param.system.nsm - 1 << "\n";

  DEBUG(10, "temperature_groups().size: " << unsigned(topo.temperature_groups().size()));

  if (!quiet)
    os << "\tSOLUTE pressure groups: "
          << unsigned(topo.pressure_groups().size()) - param.system.nsm - 1 << "\n";

  DEBUG(10, "pressure_groups().size: " << unsigned(topo.pressure_groups().size()));

  // energy group check
  if (param.force.energy_group.size() == 0) {
    param.force.energy_group.push_back(topo.num_atoms() - 1);
  }

  if (param.force.energy_group.back() > topo.num_atoms() - 1) {
    io::messages.add("Error in FORCE block: "
            "last energy group has to end with last atom",
            "In_Topology", io::message::error);
  } else if (param.force.energy_group.back() < topo.num_atoms() - 1) {
    param.force.energy_group.push_back(topo.num_atoms() - 1);
    io::messages.add("FORCE block: "
            "added an additional energy group",
            "In_Topology", io::message::warning);
  }

  // and add them
  DEBUG(10, "adding energy groups : " << param.force.energy_group.size());
  unsigned int atom = 0;
  for (unsigned int i = 0; i < param.force.energy_group.size(); ++i) {
    assert(param.force.energy_group.size() > i);
    topo.energy_groups().push_back(param.force.energy_group[i]);
    DEBUG(10, "energy group " << i << " start = " << atom << " end = " << param.force.energy_group[i]);
    for (; atom <= param.force.energy_group[i]; ++atom) {
      topo.atom_energy_group().push_back(i);
      // DEBUG(11, "atom " << atom << ": " << i);
    }
  }

  // Now that we have the energy groups, we initialize the
  // LAMBDAS parameters that depend on them.
  int maxnilg = param.force.energy_group.size();
  std::vector< double > one(maxnilg, 1.0);
  std::vector< double > zero(maxnilg, 0.0);
  for (unsigned int i = 0; i < param.lambdas.a.size(); i++) {
    // check whether we have to resize the lambdas arrays
    int lam_size = param.lambdas.a[i].size();
    if (lam_size < maxnilg) {
      param.lambdas.a[i].insert(param.lambdas.a[i].end(), maxnilg - lam_size, zero);
      lam_size = param.lambdas.b[i].size();
      param.lambdas.b[i].insert(param.lambdas.b[i].end(), maxnilg - lam_size, zero);
      lam_size = param.lambdas.c[i].size();
      param.lambdas.c[i].insert(param.lambdas.c[i].end(), maxnilg - lam_size, zero);
      lam_size = param.lambdas.d[i].size();
      param.lambdas.d[i].insert(param.lambdas.d[i].end(), maxnilg - lam_size, one);
      lam_size = param.lambdas.e[i].size();
      param.lambdas.e[i].insert(param.lambdas.e[i].end(), maxnilg - lam_size, zero);
      for (unsigned int j = 0; j < param.lambdas.a[i].size(); j++) {
        lam_size = param.lambdas.a[i][j].size();
        param.lambdas.a[i][j].insert(param.lambdas.a[i][j].end(), maxnilg - lam_size, 0.0);
        lam_size = param.lambdas.b[i][j].size();
        param.lambdas.b[i][j].insert(param.lambdas.b[i][j].end(), maxnilg - lam_size, 0.0);
        lam_size = param.lambdas.c[i][j].size();
        param.lambdas.c[i][j].insert(param.lambdas.c[i][j].end(), maxnilg - lam_size, 0.0);
        lam_size = param.lambdas.d[i][j].size();
        param.lambdas.d[i][j].insert(param.lambdas.d[i][j].end(), maxnilg - lam_size, 1.0);
        lam_size = param.lambdas.e[i][j].size();
        param.lambdas.e[i][j].insert(param.lambdas.e[i][j].end(), maxnilg - lam_size, 0.0);
      }
    } // do resize lambdas arrays
  }

  DEBUG(10, "multibath?");
  if (!param.multibath.found_multibath && param.multibath.found_tcouple) {
    if (!quiet)
      os << "\tparsing a (deprecated) TCOUPLE block into the new "
            << "MULTIBATH format.\n";

    util::parse_TCOUPLE(param, topo);
  }

  // additional known blocks not read in with this routine
  if (m_block["LJPARAMETERS"].size()) {
    block_read.insert("LJPARAMETERS");
  }
  if (m_block["CGPARAMETERS"].size()) {
    block_read.insert("CGPARAMETERS");
  }

  // warn for unread input data
  for(std::map<std::string, std::vector<std::string> >::const_iterator
      it = m_block.begin(), to = m_block.end(); it != to; ++it){
    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add("block " + it->first + " unknown and not read in!",
                       "In_Topology", io::message::warning);
    }
  }

  if (!quiet)
    os << "END\n";

  DEBUG(10, "topology read");
}


void io::In_Topology::read_bond_types(topology::Topology& topo,
        simulation::Parameter &param,
        std::ostream & os) {
  
  DEBUG(10, "reading bond types");
  if (m_block["BONDSTRETCHTYPE"].size()) {
      io::messages.add("Reading bond types from BONDSTRETCHTYPE block",
          "In_Topology::bondtype", io::message::notice);
      read_block_BONDSTRETCHTYPE(topo, param, os);
  } else if (m_block["HARMBONDTYPE"].size() && m_block["BONDTYPE"].size()) {
      io::messages.add("Reading bond types from HARMBONDTYPE and BONDTYPE block",
          "In_Topology::bondtype", io::message::notice);
      read_block_BONDTYPE(topo, param, os);
      read_block_HARMBONDTYPE(topo, param, os);

      // require that the blocks have the same length
      if (topo.bond_types_harm().size() != topo.bond_types_quart().size())
        io::messages.add("BONDTYPE and HARMBONDTYPE blocks contain different numbers of bond types",
                  "In_Topology", io::message::error);

      //check consistency of bond lengths
      for (unsigned int i=0; i< topo.bond_types_harm().size(); i++) {
        double rh=topo.bond_types_harm()[i].r0, r=topo.bond_types_quart()[i].r0;
        if (rh != r) {
          std::ostringstream msg;
          msg << "Bond lengths in BONDTYPE and HARMBONDTYPE differ (" << r << " vs. " << rh << ").";
          io::messages.add(msg.str(), "InTopology", io::message::warning);
        }
      }
  } else if (m_block["HARMBONDTYPE"].size()) {
      io::messages.add("Reading bond types from HARMBONDTYPE block",
          "In_Topology::bondtype", io::message::notice);
      if (param.force.bond == 1)
        io::messages.add("Found only bond types with harmonic force constants (HARMBONDTYPE) but NTBBH=0\n"
                  "\tCalculating quartic force constants from the harmonic ones.",
                  "In_Topology", io::message::warning);

      read_block_HARMBONDTYPE(topo, param, os);

      // calculate the quartic ones
      for (unsigned int i=0; i< topo.bond_types_harm().size(); i++) {
          double kh=topo.bond_types_harm()[i].K;
          double r=topo.bond_types_harm()[i].r0;
          double k=kh/(2*r*r);
          topo.bond_types_quart().push_back(interaction::bond_type_struct(k, r));
      }

  } else if (m_block["BONDTYPE"].size()) {
      io::messages.add("Reading bond types from BONDTYPE block",
          "In_Topology::bondtype", io::message::notice);
      if (param.force.bond == 2)
        io::messages.add("Found only bond types with quartic force constants (BONDTYPE) but NTBBH=1\n"
                  "\tCalculating harmonic force constants from the quartic ones.",
                  "In_Topology", io::message::warning);

      read_block_BONDTYPE(topo, param, os);

      // calculate the harmonic ones
      for (unsigned int i=0; i< topo.bond_types_quart().size(); i++) {
          double k=topo.bond_types_quart()[i].K;
          double r=topo.bond_types_quart()[i].r0;
          double kh=k*2*r*r;
          topo.bond_types_harm().push_back(interaction::bond_type_struct(kh, r));
      }

  } else {
        io::messages.add("Found no bond type parameters (BONDSTRETCHTYPE, BONDTYPE or HARMBONDTYPE block) in the topology.",
                  "In_Topology", io::message::error);
  }

  // check consistency: bond_types_quart and bond_types_harm have to be the same
  for (unsigned int i=0; i< topo.bond_types_harm().size(); i++) {
      double kh=topo.bond_types_harm()[i].K, k=topo.bond_types_quart()[i].K;
      double r=topo.bond_types_harm()[i].r0;
      double calc_kh = k * 2.0 * r*r;
      if (fabs(kh - calc_kh) / kh > 1.0e-4) {
        std::ostringstream msg;
        msg << "harmonic and quartic force constant do not match (CHB!=CB*2*B0*B0): " << std::endl
                << "\tCHB = " << kh << ", CB*2*B0*B0 = " << calc_kh
                << " |CHB-CB*2*B0*B0| = " << fabs(kh - calc_kh) << std::endl;
        io::messages.add(msg.str(), "InTopology", io::message::warning);
      }
  }

  num_solute_bondtypes = topo.bond_types_quart().size();

  // add solvent constraints as bond types
  std::string blockname = "SOLVENTCONSTR";
  Block block(blockname);
  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);
    unsigned int num = 0;
    block.get_next_parameter("NCONS", num, ">=0", "");

    int ic = 0, jc = 0;
    double b0 = 0.0;
    for (unsigned int i=0; i<num; i++) {
      block.get_next_parameter("ICONS", ic, ">=1", "");
      block.get_next_parameter("JCONS", jc, ">=1", "");
      block.get_next_parameter("CONS", b0, ">=0", "");

      if (block.error()) break;

      topo.bond_types_harm().push_back(interaction::bond_type_struct(0, b0));
      topo.bond_types_quart().push_back(interaction::bond_type_struct(0, b0));
    }

    block.get_final_messages(false);
  }
}

void io::In_Topology::read_bondangle_types(topology::Topology& topo,
    simulation::Parameter &param,
    std::ostream & os) {

  // avoid "block not read" messages
  block_read.insert("BONDANGLEBENDTYPE");
  block_read.insert("BONDANGLETYPE");
  block_read.insert("HARMBONDANGLETYPE");

  // read BONDANGLEBENDTYPE, or, if it does not exist, any other bond angle
  // type blocks that exist
  if (m_block["BONDANGLEBENDTYPE"].size()) {
    io::messages.add("Reading angle types from BONDANGLEBENDTYPE block",
        "In_Topology::angletype", io::message::notice);
    read_block_BONDANGLEBENDTYPE(topo, param, os);
  } else {
    if (m_block["BONDANGLETYPE"].size()) {
        io::messages.add("Reading angle types from BONDANGLETYPE block",
          "In_Topology::angletype", io::message::notice);
        read_block_BONDANGLETYPE(topo, param, os);
    }
    if (m_block["HARMBONDANGLETYPE"].size()) {
        io::messages.add("Reading angle types from HARMBONDANGLETYPE block",
          "In_Topology::angletype", io::message::notice);
        read_block_HARMBONDANGLETYPE(topo, param, os);
    }
  }

  // now check if we have the ones that we need
  if (param.force.angle == 2) {
      // harmonic potential
      if (!m_block["HARMBONDANGLETYPE"].size() && !m_block["BONDANGLEBENDTYPE"].size()) {
        io::messages.add("Harmonic angle potential requested (NTBAH=1) but "
            "no BONDANGLEBENDTYPE or HARMBONDANGLETYPE block found.",
            "In_Topology", io::message::error);
      }
  } else if (param.force.angle == 1) {
      // cosine harmonic potential
      if (!m_block["BONDANGLETYPE"].size() && !m_block["BONDANGLEBENDTYPE"].size()) {
        io::messages.add("Cosine armonic angle potential requested (NTBAH=0) but "
            "no BONDANGLEBENDTYPE or BONDANGLETYPE block found.",
            "In_Topology", io::message::error);
      }
  } 

  // either both lists have the same length or one of them has size 0:
  num_solute_angletypes = topo.angle_types_cosharm().size();
  if (num_solute_angletypes == 0) num_solute_angletypes = topo.angle_types_harm().size();

}


void io::In_Topology::read_dihedral_types(topology::Topology& topo,
    simulation::Parameter &param,
    std::ostream & os) {

  // avoid "block not read" messages
  block_read.insert("TORSDIHEDRALTYPE");
  block_read.insert("DIHEDRALTYPE");

  if (m_block["TORSDIHEDRALTYPE"].size()) {
    io::messages.add("Reading angle types from TORSDIHEDRALTYPE block",
        "In_Topology::dihedraltype", io::message::notice);
    read_block_TORSDIHEDRALTYPE(topo, param, os);
  } else if (m_block["DIHEDRALTYPE"].size()) {
    if (param.force.dihedral == 1) {
        io::messages.add("DIHEDRALTYPE block can not be used with NTBDN=0 (COVALENTFORM block).",
            "In_Topology", io::message::error);
    }
    io::messages.add("Reading angle types from DIHEDRALTYPE block",
        "In_Topology::dihedraltype", io::message::notice);
    read_block_DIHEDRALTYPE(topo, param, os);
  } else {
    io::messages.add("Either TORSDIHEDRALTYPE or DIHEDRALTYPE block must be present.",
        "In_Topology", io::message::error);
  }
  num_solute_dihedraltypes=topo.dihedral_types().size();
}

void io::In_Topology::read_block_BONDSTRETCHTYPE(topology::Topology &topo,
        simulation::Parameter &param, std::ostream & os) {

  std::string blockname = "BONDSTRETCHTYPE";
  Block block(blockname);

  if (block.read_buffer(m_block[blockname], false) == 0) {
    DEBUG(7, "reading in "+blockname+" block");

    block_read.insert(blockname);

    int num = 0;
    double k = 0.0, kh = 0.0, r = 0.0;
    block.get_next_parameter("NBTY", num, ">=0", "");

    if (num<=0 || block.error() ) return;

    for (int i=0; i<num; i++) {
      block.get_next_parameter("k", k, ">=0", "");
      block.get_next_parameter("kh", kh, ">=0", "");
      block.get_next_parameter("r", r, ">=0", "");
      if (block.error()) {
        std::string linenumber=io::to_string(num+1);
        io::messages.add("Bad values in line "+linenumber+" in "+blockname+" block",
                "In_Topology",
                io::message::error);
        break; // TODO: should we break here or rather just set k and r to 0 and loop through everything?
      } else {
        topo.bond_types_harm().push_back(interaction::bond_type_struct(kh, r));
        topo.bond_types_quart().push_back(interaction::bond_type_struct(k, r));
      }
    }

    block.get_final_messages(false);
  }

}

void io::In_Topology::read_block_BONDTYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) {

  std::string blockname = "BONDTYPE";
  Block block(blockname);

  if (block.read_buffer(m_block[blockname], false) == 0) {
    DEBUG(7, "reading in "+blockname+" block");

    block_read.insert(blockname);

    int num = 0;
    double k = 0.0, r = 0.0;
    block.get_next_parameter("NBTY", num, ">=0", "");

    if (num<=0 || block.error() ) return;

    for (int i=0; i<num; i++) {
      block.get_next_parameter("k", k, ">=0", "");
      block.get_next_parameter("r", r, ">=0", "");
      if (block.error()) {
        std::string linenumber=io::to_string(num+1);
        io::messages.add("Bad values in line "+linenumber+" in "+blockname+" block",
                "In_Topology",
                io::message::error);
        break; // TODO: should we break here or rather just set k and r to 0 and loop through everything?
      } else {
        topo.bond_types_quart().push_back(interaction::bond_type_struct(k, r));
      }
    }

    block.get_final_messages(false);
  }
}

void io::In_Topology::read_block_HARMBONDTYPE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) {

  std::string blockname = "HARMBONDTYPE";
  Block block(blockname);

  if (block.read_buffer(m_block[blockname], false) == 0) {
    DEBUG(7, "reading in "+blockname+" block");

    block_read.insert(blockname);

    int num = 0;
    double kh = 0.0, r = 0.0;
    block.get_next_parameter("NBTY", num, ">=0", "");

    if (num<=0 || block.error() ) return;

    for (int i=0; i<num; i++) {
      block.get_next_parameter("kh", kh, ">=0", "");
      block.get_next_parameter("r", r, ">=0", "");
      if (block.error()) {
        std::string linenumber=io::to_string(num+1);
        io::messages.add("Bad values in line "+linenumber+" in "+blockname+" block",
                "In_Topology",
                io::message::error);
        break; // TODO: should we break here or rather just set k and r to 0 and loop through everything?
      } else {
        topo.bond_types_harm().push_back(interaction::bond_type_struct(kh, r));
      }
    }

    block.get_final_messages(false);
  }
}

void io::In_Topology::read_block_BONDANGLEBENDTYPE(topology::Topology &topo,
        simulation::Parameter &param, std::ostream & os) {

  std::string blockname = "BONDANGLEBENDTYPE";
  Block block(blockname);

  if (block.read_buffer(m_block[blockname], false) == 0) {
    DEBUG(7, "reading in "+blockname+" block");

    block_read.insert(blockname);

    int num = 0;
    double kcosh = 0.0, kh = 0.0, theta = 0.0;
    block.get_next_parameter("NTTY", num, ">=0", "");

    if (num<=0 || block.error() ) return;

    for (int i=0; i<num; i++) {
      block.get_next_parameter("CT", kcosh, ">=0", "");
      block.get_next_parameter("CHT", kh, ">=0", "");
      block.get_next_parameter("T0", theta, "", "");
      if (block.error()) {
        std::string linenumber=io::to_string(num+1);
        io::messages.add("Bad values in line "+linenumber+" in "+blockname+" block",
                "In_Topology",
                io::message::error);
        break;
      } else {
        topo.angle_types_harm().push_back(interaction::angle_type_struct(kh * (180.0 * 180.0 / (math::Pi * math::Pi)),
               theta *  math::Pi / 180.0));
        topo.angle_types_cosharm().push_back(interaction::angle_type_struct(kcosh, cos(theta * math::Pi / 180.0)));
      }
    }

    block.get_final_messages(false);
  }

}

void io::In_Topology::read_block_BONDANGLETYPE(topology::Topology &topo,
        simulation::Parameter &param, std::ostream & os) {

  std::string blockname = "BONDANGLETYPE";
  Block block(blockname);

  if (block.read_buffer(m_block[blockname], false) == 0) {
    DEBUG(7, "reading in "+blockname+" block");

    block_read.insert(blockname);

    int num = 0;
    double kcosh = 0.0, theta = 0.0;
    block.get_next_parameter("NTTY", num, ">=0", "");

    if (num<=0 || block.error() ) return;

    for (int i=0; i<num; i++) {
      block.get_next_parameter("CT", kcosh, ">=0", "");
      block.get_next_parameter("T0", theta, "", "");
      if (block.error()) {
        std::string linenumber=io::to_string(num+1);
        io::messages.add("Bad values in line "+linenumber+" in "+blockname+" block",
                "In_Topology",
                io::message::error);
        break;
      } else {
        topo.angle_types_cosharm().push_back(interaction::angle_type_struct(kcosh, cos(theta * math::Pi / 180.0)));
      }
    }

    block.get_final_messages(false);
  }

}

void io::In_Topology::read_block_HARMBONDANGLETYPE(topology::Topology &topo,
        simulation::Parameter &param, std::ostream & os) {

  std::string blockname = "HARMBONDANGLETYPE";
  Block block(blockname);

  if (block.read_buffer(m_block[blockname], false) == 0) {
    DEBUG(7, "reading in "+blockname+" block");

    block_read.insert(blockname);

    int num = 0;
    double kh = 0.0, theta = 0.0;
    block.get_next_parameter("NTTY", num, ">=0", "");

    if (num<=0 || block.error() ) return;

    for (int i=0; i<num; i++) {
      block.get_next_parameter("CHT", kh, ">=0", "");
      block.get_next_parameter("T0", theta, "", "");
      if (block.error()) {
        std::string linenumber=io::to_string(num+1);
        io::messages.add("Bad values in line "+linenumber+" in "+blockname+" block",
                "In_Topology",
                io::message::error);
        break;
      } else {
        topo.angle_types_harm().push_back(interaction::angle_type_struct(kh * (180.0 * 180.0 / (math::Pi * math::Pi)),
               theta *  math::Pi / 180.0));
      }
    }

    block.get_final_messages(false);
  }

}

void io::In_Topology::read_block_TORSDIHEDRALTYPE(topology::Topology &topo,
        simulation::Parameter &param, std::ostream & os) {
  // dihedrals in degrees
  std::string blockname = "TORSDIHEDRALTYPE";
  Block block(blockname);

  if (block.read_buffer(m_block[blockname], false) == 0) {
    DEBUG(7, "reading in "+blockname+" block");

    block_read.insert(blockname);

    int num = 0, m = 0;
    double k = 0.0, pdl = 0.0;
    bool asymmetric=false;
    block.get_next_parameter("NPTY", num, ">=0", "");

    if (num<=0 || block.error() ) return;

    for (int i=0; i<num; i++) {
      block.get_next_parameter("CP", k, ">=0", "");
      block.get_next_parameter("PD", pdl, "", "");
      block.get_next_parameter("NP", m, ">0", "");
      if (param.force.dihedral == 2) {
        if (pdl != 0 && pdl != 180)
          io::messages.add("TORSDIHEDRALTYPE block: dihedral neither 0 nor 180, this is not allowed if NTBDN=1 (COVALENTFORM block)", "In_Topology",
                io::message::error);
      } else if (param.force.dihedral == 1) {
          if (pdl!=0 || pdl!=180) {
            asymmetric=true;
          }
      }
      if (block.error()) {
        std::string linenumber=io::to_string(num+1);
        io::messages.add("Bad values in line "+linenumber+" in "+blockname+" block",
                "In_Topology",
                io::message::error);
        break;
      } else {
        topo.dihedral_types().push_back(interaction::dihedral_type_struct(k, cos(pdl * math::Pi / 180.0), pdl * math::Pi / 180.0, m));
      }
    }
    if (!asymmetric)
        io::messages.add("TORSDIHEDRALTYPE block contains no values different from 0 or 180, you could also use NTBDN=1 (COVALENTFORM block)", "In_Topology",
              io::message::notice);

    block.get_final_messages(false);
  }

}

void io::In_Topology::read_block_DIHEDRALTYPE(topology::Topology &topo,
        simulation::Parameter &param, std::ostream & os) {
  // dihedrals as cos(theta)
  std::string blockname = "DIHEDRALTYPE";
  Block block(blockname);

  if (block.read_buffer(m_block[blockname], false) == 0) {
    DEBUG(7, "reading in "+blockname+" block");

    block_read.insert(blockname);

    int num = 0, m = 0;
    double k = 0.0, pd = 0.0;
    block.get_next_parameter("NPTY", num, ">=0", "");

    if (num<=0 || block.error() ) return;

    for (int i=0; i<num; i++) {
      block.get_next_parameter("CP", k, ">=0", "");
      block.get_next_parameter("PD", pd, "", "1, -1");
      block.get_next_parameter("NP", m, ">0", "");
      if (block.error()) {
        std::string linenumber=io::to_string(num+1);
        io::messages.add("Bad values in line "+linenumber+" in "+blockname+" block",
                "In_Topology",
                io::message::error);
        break;
      } else {
        topo.dihedral_types().push_back(interaction::dihedral_type_struct(k, pd, acos(pd), m));
      }
    }

    block.get_final_messages(false);
  }

}
void io::In_Topology::read_block_IMPDIHEDRALTYPE(topology::Topology &topo,
        simulation::Parameter &param, std::ostream & os) {
  std::string blockname = "IMPDIHEDRALTYPE";
  Block block(blockname);

  if (block.read_buffer(m_block[blockname], false) == 0) {
    DEBUG(7, "reading in "+blockname+" block");

    block_read.insert(blockname);

    int num = 0;
    double k = 0.0, q0 = 0.0;
    block.get_next_parameter("NQTY", num, ">=0", "");

    if (num<=0 || block.error() ) return;

    for (int i=0; i<num; i++) {
      block.get_next_parameter("CQ", k, ">=0", "");
      block.get_next_parameter("Q0", q0, "", "");
      if (block.error()) {
        std::string linenumber=io::to_string(num+1);
        io::messages.add("Bad values in line "+linenumber+" in "+blockname+" block",
                "In_Topology",
                io::message::error);
        break;
      } else {
        topo.impdihedral_types().push_back(interaction::improper_dihedral_type_struct(k * 180 * 180 / math::Pi / math::Pi,
                q0 * math::Pi / 180.0));
      }
    }

    num_solute_impropertypes=topo.impdihedral_types().size();
    block.get_final_messages(false);
  }

}

void io::In_Topology
::read_lj_parameter(std::vector<std::vector
        <interaction::lj_parameter_struct> >
        & lj_parameter,
        std::ostream & os) {
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  { // LJPARAMETERS

    DEBUG(10, "LJPARAMETERS block");

    buffer = m_block["LJPARAMETERS"];
    if (!buffer.size()) {
      io::messages.add("No LJPARAMETERS block found in topology!",
              "In_Topology",
              io::message::error);
      return;
    }

    int num = 0, n = 0;

    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num;

    // calculate the matrix size from: x = n*(n+1)/2
    unsigned int sz = unsigned(sqrt(double((8 * num + 1) - 1)) / 2);

    lj_parameter.resize(sz);
    std::vector< std::vector<interaction::lj_parameter_struct> >::iterator
    lj_it = lj_parameter.begin(),
            lj_to = lj_parameter.end();

    for (; lj_it != lj_to; ++lj_it)
      lj_it->resize(sz);

    ++it;

    for (n = 0; it != buffer.end() - 1; ++it, ++n) {

      interaction::lj_parameter_struct s;
      int i = 0, j = 0;

      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> i >> j >> s.c12 >> s.c6 >> s.cs12 >> s.cs6;

      --i;
      --j;

      if (_lineStream.fail() )
        io::messages.add("bad line in LJPARAMETERS block", "In_Topology",
              io::message::error);

      if (i >= int(sz) || j >= int(sz)) {
        DEBUG(7, "wrong iac in LJPARAMETERS: i=" << i << " j=" << j
                << " sz=" << sz);
        io::messages.add("wrong integer atom code in LJPARAMETERS block",
                "In_Topology",
                io::message::error);
      }

      lj_parameter[i][j] = s;
      lj_parameter[j][i] = s;

    }

    if (num != n) {
      io::messages.add("Reading the LJPARAMETERS failed (n != num)",
              "InTopology",
              io::message::error);
    }
  } // LJPARAMETER

}

void io::In_Topology
::read_cg_parameter(std::vector<std::vector
        <interaction::lj_parameter_struct> >
        & cg_parameter,
        std::ostream & os) {
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  { // CGPARAMETERS

    DEBUG(10, "CGPARAMETERS block");

    buffer = m_block["CGPARAMETERS"];
    if (!buffer.size()) {
      io::messages.add("No CGPARAMETERS block found in topology!",
              "In_Topology",
              io::message::error);
      return;
    }

    int num = 0, n = 0;

    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num;

    // calculate the matrix size from: x = n*(n+1)/2
    unsigned int sz = unsigned(sqrt(double((8 * num + 1) - 1)) / 2);

    cg_parameter.resize(sz);
    std::vector< std::vector<interaction::lj_parameter_struct> >::iterator
    cg_it = cg_parameter.begin(),
            cg_to = cg_parameter.end();

    for (; cg_it != cg_to; ++cg_it)
      cg_it->resize(sz);

    ++it;

    for (n = 0; it != buffer.end() - 1; ++it, ++n) {

      interaction::lj_parameter_struct s;
      int i = 0, j = 0;

      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> i >> j >> s.c12 >> s.c6;

      --i;
      --j;

      if (_lineStream.fail() || !_lineStream.eof())
        io::messages.add("bad line in CGPARAMETERS block", "In_Topology",
              io::message::error);

      if (i >= int(sz) || j >= int(sz)) {
        DEBUG(7, "wrong iac in CGPARAMETERS: i=" << i << " j=" << j
                << " sz=" << sz);
        io::messages.add("wrong integer atom code in CGPARAMETERS block",
                "In_Topology",
                io::message::error);
        break;
      }

      // no different 1,4 interactions
      s.cs6 = s.c6;
      s.cs12 = s.c12;

      cg_parameter[i][j] = s;
      cg_parameter[j][i] = s;

    }

    if (num != n) {
      io::messages.add("Reading the CGPARAMETERS failed (n != num)",
              "InTopology",
              io::message::error);
    }

  } // CGPARAMETER
}

void io::In_Topology::
read_sasa_parameter(topology::Topology & topo, std::vector<topology::sasa_parameter_struct>
    & sasa_parameter) {
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  { // SASAPARAMETERS

    buffer = m_block["SASAPARAMETERS"];
    block_read.insert("SASAPARAMETERS");
    // if no SASA block is present
    if (!buffer.size()) {
      io::messages.add("No SASAPARAMETERS block found in topology!",
                      "In_Topology", io::message::error);
      return;
    }

    unsigned int num = 0; // number of sasa atoms;

    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num; // reads in number of sasa atoms
    ++it;
    DEBUG(1, "Number of sasa atoms : " << num);

    for (; it != buffer.end() - 1; ++it) {

      topology::sasa_parameter_struct s;

      _lineStream.clear();
      _lineStream.str(*it);

      int atom = 0;
      _lineStream >> atom >> s.r >> s.p >> s.sigma;

      if (_lineStream.fail() || !_lineStream.eof())
        io::messages.add("bad line in SASAPARAMETERS block",
              "In_Topology", io::message::error);

      if (atom < 0 || atom > int(topo.num_solute_atoms()))
        io::messages.add("atom out of range in SASAPARAMETERS block",
          "In_Topology", io::message::error);
      s.atom = atom - 1; // convert human atom numbers to internal atom numbers

      if (s.r < 0) {
        DEBUG(7, "negative radius in SASAPARAMETERS: r = " << s.r);
        io::messages.add("negative radius in SASAPARAMETERS block",
                "In_Topology", io::message::error);
      }

      if (s.p < 0) {
        DEBUG(7, "negative probability in SASAPARAMETERS: r = " << s.p);
        io::messages.add("negative probability in SASAPARAMETERS block",
                "In_Topology", io::message::error);
      }

      sasa_parameter.push_back(s);
    }

    if (num != sasa_parameter.size()) {
      io::messages.add("Number of SASAPARAMETERS not equal to number of sasa atoms",
              "InTopology", io::message::error);
    }
  } // SASAPARAMETERS
}


void io::In_Topology::read_block_TYPE(topology::Topology& topo,
        simulation::Parameter &param,
        std::ostream & os)   { //TYPE

    std::vector<std::string> buffer = m_block["TYPE"];
    if (buffer.size()) {
      block_read.insert("TYPE");
      if (buffer.size() != 3) {
        io::messages.add("Bad line in TYPE block",
                "InTopology", io::message::error);
      }

      _lineStream.clear();
      _lineStream.str(buffer[1]);
      std::string s;
      _lineStream >> s;

      if (_lineStream.fail())
        io::messages.add("Bad line in TYPE block",
              "InTopology", io::message::error);

      if (s.length() > MAX_NAME) {
        std::ostringstream msg;
        msg << "Error in TYPE block: type " << s
                << " is too long (> " << MAX_NAME << " characters).";
        io::messages.add(msg.str(), "InTopology", io::message::error);
      }


      std::transform(s.begin(), s.end(), s.begin(), tolower);

      if (s == "atomistic") {
        if (!quiet)
          os << "\tatomistic topology\n";

        if (param.force.interaction_function != simulation::lj_crf_func &&
                param.force.interaction_function != simulation::pol_lj_crf_func) {
          io::messages.add("wrong interaction function selected for atomistic topology",
                  "InTopology", io::message::error);
        }
      } else if (s == "coarse-grained") {
        if (!quiet)
          os << "\tcoarse-grained topology\n";

        if (param.force.interaction_function != simulation::cgrain_func) {
          io::messages.add("wrong interaction function selected for coarse-grained topology",
                  "InTopology", io::message::error);
        }
      } else {
        io::messages.add("TYPE block: unknown topology type",
                "InTopology", io::message::error);
      }
    } else {
      if (!quiet) {
        os << "\tunknown topology type (atomistic / coarse-grained)\n";
        if (param.force.interaction_function == simulation::lj_crf_func ||
                param.force.interaction_function == simulation::pol_lj_crf_func ||
                param.force.interaction_function == simulation::pol_off_lj_crf_func)
          os << "\tusing atomistic parameters\n";
        else if (param.force.interaction_function == simulation::cgrain_func)
          os << "\tusing coarse-grained parameters\n";
        else
          os << "\tunknown interaction function selected!\n";
      }
    }
} //TYPE



void io::In_Topology::read_block_PHYSICALCONSTANTS(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)  { // PHYSICALCONSTANTS
    std::vector<std::string> buffer = m_block["PHYSICALCONSTANTS"];
    if (buffer.size()) {
      block_read.insert("PHYSICALCONSTANTS");
      std::string s;
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin() + 1,
              buffer.end() - 1, s));
      double four_pi_eps0_i = 0.0;

      _lineStream >> four_pi_eps0_i >> math::h_bar >> math::spd_l >> math::k_Boltzmann;
      math::eps0_i = four_pi_eps0_i * 4.0 * math::Pi;
      math::four_pi_eps_i = four_pi_eps0_i / param.nonbonded.epsilon;

      if (_lineStream.fail())
        io::messages.add("Bad line in PHYSICALCONSTANTS block",
              "InTopology", io::message::error);
    } else {
      io::messages.add("no PHYSICALCONSTANTS block in topology",
              "InTopology", io::message::error);
    }

} // PHYSICALCONSTANTS

void io::In_Topology::read_block_RESNAME(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) { // RESNAME
      if (!quiet)
        os << "\tRESNAME\n\t";

      DEBUG(10, "RESNAME block");
      std::vector<std::string> buffer = m_block["RESNAME"];
      block_read.insert("RESNAME");
      std::vector<std::string>::const_iterator it = buffer.begin() + 1;
      int n = 0, num = 0;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> num;
      ++it;

      if (!quiet && num > 10) {
        for (n = 0; n < 10; ++n)
          os << std::setw(8) << n + 1;
        os << "\n\t";
      }

      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        std::string s;
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> s;

        if (!quiet) {
          if (n && ((n) % 10) == 0) os << std::setw(10) << n << "\n\t";
          os << std::setw(8) << s;
        }

        if (s.length() > MAX_NAME) {
          std::ostringstream msg;
          msg << "Error in RESNAME block: residue name " << s
                  << " is too long (> " << MAX_NAME << " characters).";
          io::messages.add(msg.str(), "InTopology", io::message::error);
        }

        topo.residue_names().push_back(s);
      }

      if (n != num) {
        io::messages.add("Error in RESNAME block: n!=num.",
                "InTopology", io::message::error);
      }

      if (!quiet)
        os << "\n\tEND\n";

} // RESNAME

void io::In_Topology::read_block_ATOMTYPENAME(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)  { // ATOMTYPENAME
      if (!quiet)
        os << "\tATOMTYPENAME\n\t";

      DEBUG(10, "ATOMTYPENAME block");
      std::vector<std::string> buffer = m_block["ATOMTYPENAME"];
      block_read.insert("ATOMTYPENAME");
      std::vector<std::string>::const_iterator it = buffer.begin() + 1;
      int n = 0, num = 0;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> num;
      ++it;

      topo.set_num_atomtype(num);

      if (!quiet)
        os << "\t" << num << " atom types\n";

      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        std::string s;
        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> s;

        if (s.length() > MAX_NAME) {
          std::ostringstream msg;
          msg << "Error in ATOMTYPENAME block: type " << s
                  << " is too long (> " << MAX_NAME << " characters).";
          io::messages.add(msg.str(), "InTopology", io::message::error);
        }

        if (topo.atom_names().find(s) == topo.atom_names().end()) {  // not found
          topo.atom_names()[s] = n;
        } else {  // found
          std::ostringstream msg;
          msg << "Error in ATOMTYPENAME block: atom type name " << s
                  << " used more than once.";
          io::messages.add(msg.str(), "InTopology", io::message::error);
        }
      }

      if (n != num) {
        io::messages.add("Error in ATOMTYPENAME block: n!=num.",
                "InTopology", io::message::error);
      }

      if (!quiet)
        os << "\tEND\n";

} // ATOMTYPENAME

void io::In_Topology::read_block_SOLUTEATOM(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)   { // SOLUTEATOM
      DEBUG(10, "SOLUTEATOM block");
      std::vector<std::string> buffer = m_block["SOLUTEATOM"];
      block_read.insert("SOLUTEATOM");

      std::vector<std::string>::const_iterator it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num = 0, n = 0;
      _lineStream >> num;
      topo.resize(num);

      if (!quiet)
        os << "\tSOLUTEATOM\n\t"
              << "\tnumber of atoms : " << num;

      // put the rest of the block into a single stream
      ++it;

      // std::string soluteAtoms;
      // concatenate(it, buffer.end()-1, soluteAtoms);

      // _lineStream.clear();
      // _lineStream.str(soluteAtoms);

      // os << "\ntime after concatenate: " << util::now() - start << std::endl;

      int a_nr = 0, r_nr = 0, t = 0, cg = 0, n_ex = 0, a_ex = 0;
      double m = 0.0, q = 0.0;
      std::string s;
      topology::excl_cont_t::value_type ex;
      topology::excl_cont_t::value_type ex14;

      for (n = 0; n < num; ++n) {

        _lineStream.clear();
        _lineStream.str(*it);
        ++it;

        _lineStream >> a_nr >> r_nr >> s >> t >> m >> q >> cg;

        if (a_nr != n + 1) {
          io::messages.add("Error in SOLUTEATOM block: atom number not sequential.",
                  "InTopology", io::message::error);
        }

        if (r_nr > int(topo.residue_names().size()) || r_nr < 1) {
          io::messages.add("Error in SOLUTEATOM block: residue number out of range.",
                  "InTopology", io::message::error);
        }

        if (s.length() > MAX_NAME) {
          std::ostringstream msg;
          msg << "Error in SOLUTEATOM block: atom name " << s
                  << " is too long (> " << MAX_NAME << " characters).";
          io::messages.add(msg.str(), "InTopology", io::message::error);
        }


        if (t < 1) {
          io::messages.add("Error in SOLUTEATOM block: iac < 1.",
                  "InTopology", io::message::error);
        }

        if (t > topo.num_atomtype()) {
          io::messages.add("Error in SOLUTEATOM block: iac > number of atom types.",
                  "InTopology", io::message::error);
        }

        if (m <= 0) {
          io::messages.add("Error in SOLUTEATOM block: mass < 0.",
                  "InTopology", io::message::error);
        }

        if (cg != 0 && cg != 1) {
          io::messages.add("Error in SOLUTEATOM block: cg = 0 / 1.",
                  "InTopology", io::message::error);
        }

        if (!(_lineStream >> n_ex)) {
          if (_lineStream.eof()) {
            _lineStream.clear();
            _lineStream.str(*it);
            ++it;
            _lineStream >> n_ex;
          } else {
            io::messages.add("Error in SOLUTEATOM block: number of exclusions "
                    "could not be read.",
                    "InTopology", io::message::error);
          }
        }

        if (n_ex < 0) {
          io::messages.add("Error in SOLUTEATOM block: number of exclusions < 0.",
                  "InTopology", io::message::error);
        }

        // exclusions
        ex.clear();
        for (int i = 0; i < n_ex; ++i) {
          if (!(_lineStream >> a_ex)) {
            if (_lineStream.eof()) {
              _lineStream.clear();
              _lineStream.str(*it);
              ++it;
              _lineStream >> a_ex;
            } else {
              io::messages.add("Error in SOLUTEATOM block: exclusion "
                      "could not be read.",
                      "InTopology", io::message::error);
            }
          }

          if (a_ex <= a_nr)
            io::messages.add("Error in SOLUTEATOM block: exclusions only to "
                  "larger atom numbers.",
                  "InTopology", io::message::error);

          ex.insert(a_ex - 1);
        }

        // 1,4 - pairs
        if (!(_lineStream >> n_ex)) {
          if (_lineStream.eof()) {
            _lineStream.clear();
            _lineStream.str(*it);
            ++it;
            _lineStream >> n_ex;
          } else {
            io::messages.add("Error in SOLUTEATOM block: number of 1,4 - exclusion "
                    "could not be read.",
                    "InTopology", io::message::error);
          }
        }

        if (n_ex < 0) {
          io::messages.add("Error in SOLUTEATOM block: number of 1,4 exclusions < 0.",
                  "InTopology", io::message::error);
        }

        ex14.clear();
        for (int i = 0; i < n_ex; ++i) {
          if (!(_lineStream >> a_ex)) {
            if (_lineStream.eof()) {
              _lineStream.clear();
              _lineStream.str(*it);
              ++it;
              _lineStream >> a_ex;
            } else {
              io::messages.add("Error in SOLUTEATOM block: 1,4 - exclusion "
                      "could not be read.",
                      "InTopology", io::message::error);
            }
          }

          if (a_ex <= a_nr)
            io::messages.add("Error in SOLUTEATOM block: 1,4 - exclusions only to "
                  "larger atom numbers.",
                  "InTopology", io::message::error);

          ex14.insert(a_ex - 1);
        }

        if (_lineStream.fail())
          io::messages.add("bad line in SOLUTEATOM block",
                "In_Topology",
                io::message::critical);

        topo.add_solute_atom(s, r_nr - 1, t - 1, m, q, cg, ex, ex14);
      }
      if (!quiet)
        os << "\n\tEND\n";

} // SOLUTEATOM


void io::In_Topology::read_block_SOLUTEPOLARISATION(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)
    { // SOLUTEPOLARISATION
      DEBUG(10, "SOLUTEPOLARISATION block");
      std::vector<std::string> buffer = m_block["SOLUTEPOLARISATION"];

      if (buffer.size()) {
        block_read.insert("SOLUTEPOLARISATION");
        if (!quiet)
          os << "\tSOLUTEPOLARISATION\n"
                << "\t\tblock present\n"
                << "\tEND\n";
        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, k = 0;
          double polarisability = 0.0, coscharge = 0.0, damping_level = 0.0, damping_power = 0.0, gamma = 0.0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> polarisability >> coscharge >> damping_level
                  >> damping_power >> gamma >> j >> k;
          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in SOLUTEPOLARISATION block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || i < 1) {
            io::messages.add("Atom number out of range in SOLUTEPOLARISATION block",
                    "In_Topology", io::message::error);
          } else if (gamma != 0.0 && (j < 1 || k < 1 || j > int(topo.num_solute_atoms())
                  || k > int(topo.num_solute_atoms()) || i == k || i == j || k == j)) {
            io::messages.add("Atom number for off atom out of range in SOLUTEPOLARISATION block",
                    "In_Topology", io::message::error);
          } else {
            DEBUG(10, "\tpolarisable atom: " << i);
            topo.polarisability()[i - 1] = polarisability / math::four_pi_eps_i;
            topo.coscharge()[i - 1] = coscharge;
            topo.damping_level()[i - 1] = damping_level * sqrt(math::four_pi_eps_i);
            topo.damping_power()[i - 1] = damping_power;
            topo.is_polarisable()[i - 1] = bool(polarisability > 0.0);
            if (param.polarise.cos == 2) {
              topo.gamma()[i - 1] = 2 * gamma;
              topo.gamma_j()[i - 1] = j - 1;
              topo.gamma_k()[i - 1] = k - 1;
            }

          }
        }

        if (n != num) {
          io::messages.add("Wrong number of polarisable atoms in SOLUTEPOLARISATION block",
                  "In_Topology", io::message::error);
        }
      }
} // SOLUTEPOLARISATION

void io::In_Topology::read_block_CGSOLUTE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)
    { // CGSOLUTE
      DEBUG(10, "CGSOLUTE block");

      if (!quiet)
        os << "\tCGSOLUTE\n";

      std::vector<std::string> buffer = m_block["CGSOLUTE"];
      if (buffer.size()) {
        block_read.insert("CGSOLUTE");
        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;
        os << "\t\tnumber of ranges: " << num << "\n";

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int cg_begin = 0, cg_end = 0, cg_fac = 0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> cg_begin >> cg_end >> cg_fac;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in CGSOLUTE block",
                    "In_Topology", io::message::error);
          }

          if ((cg_begin > int(topo.num_solute_atoms()) || cg_begin < 1)
                  || (cg_end > int(topo.num_solute_atoms()) || cg_end < 1)
                  || (cg_begin > cg_end)) {
            io::messages.add("Sequence number out of range in CGSOLUTE block",
                    "In_Topology", io::message::error);
          } else if (cg_fac < 1) {
            io::messages.add("CG factor out of range in CGSOLUTE block",
                    "In_Topology", io::message::error);
          } else {
            DEBUG(10, "\tcoarse grained range: " << cg_begin << " " << cg_end);
            DEBUG(10, "\tcoarse grained factor: " << cg_fac);
            os << "\t\trange: " << cg_begin << " " << cg_end << "\n"
               << "\t\tcoarse grained factor: " << cg_fac << "\n";
            for (unsigned int i = (cg_begin - 1); i < unsigned(cg_end); ++i) {
              topo.is_coarse_grained()[i] = true;
              topo.cg_factor()[i] = cg_fac;
            }
          }
        }
        os << "\tEND\n";

        if (n != num) {
          io::messages.add("Wrong number of ranges in CGSOLUTE block",
                  "In_Topology", io::message::error);
        }
      }
} // CGSOLUTE


void io::In_Topology::read_block_LJEXCEPTIONS(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)
     { // LJEXCEPTIONS
      DEBUG(10, "LJEXCEPTIONS block");
      std::vector<std::string> buffer = m_block["LJEXCEPTIONS"];

      if (buffer.size()) {
        block_read.insert("LJEXCEPTIONS");
        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;
        if (_lineStream.fail() || !_lineStream.eof()) {
          io::messages.add("Bad line in LJEXCEPTIONS block",
                  "In_Topology", io::message::error);
        } else {
          if (!quiet) {
            os << "\tLJEXCEPTIONS\n"
                    << "\t\t" << num << " Lennard-Jones exceptions.\n"
                    << "\tEND\n";
          }
        }

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0;
          double c6 = 0.0, c12 = 0.0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> c12 >> c6;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in LJEXCEPTIONS block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || i < 1 ||
                  j > int(topo.num_solute_atoms()) || j < 1) {
            io::messages.add("Atom number out of range in LJEXCEPTIONS block",
                    "In_Topology", io::message::error);
          } else {
            i--;
            j--;
            if (i > j)
              std::swap(i, j);
            topo.lj_exceptions().push_back(topology::lj_exception_struct(i, j, c6, c12));
            topo.all_exclusion(i).insert(j);
          }
        }

        if (n != num) {
          io::messages.add("Wrong number of lines in LJEXCEPTIONS block",
                  "In_Topology", io::message::error);
        }
      }
} // LJEXCEPTION

void io::In_Topology::read_block_BONDH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)
      { // BONDH
      DEBUG(10, "BONDH block");

      if (!quiet)
        os << "\tBOND";

      std::vector<std::string> buffer = m_block["BONDH"];
      if (buffer.size()) {
        block_read.insert("BONDH");
        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet) {
          if (param.constraint.ntc == 2 || param.constraint.ntc == 3) {
            os << "\n\t\t"
                    << num
                    << " bonds from BONDH block added to CONSTRAINT";
          } else
            os << "\n\t\tbonds containing hydrogens : "
                  << num;
        }

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, t = 0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in BONDH block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1) {
            io::messages.add("Atom number out of range in BONDH block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in BONDH block: bond type < 1",
                    "In_Topology", io::message::error);
          }

          if (param.constraint.ntc == 2 || param.constraint.ntc == 3) {
            topo.solute().distance_constraints().
                    push_back(topology::two_body_term_struct(i - 1, j - 1, t - 1));
          } else
            topo.solute().bonds().
            push_back(topology::two_body_term_struct(i - 1, j - 1, t - 1));

        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in BONDH block",
                  "In_Topology", io::message::error);
        }
      }

    } // BONDH

void io::In_Topology::read_block_BOND(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)
       { // BOND
      DEBUG(10, "BOND block");
      std::vector<std::string> buffer = m_block["BOND"];

      if (buffer.size()) {
        block_read.insert("BOND");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet) {
          if (param.constraint.ntc == 3) {
            os << "\n\t\t"
                    << num
                    << " bonds from BOND block added to CONSTRAINT";
          } else
            os << "\n\t\tbonds not containing hydrogens : "
                  << num;
        }

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, t = 0;

          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in BOND block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1) {
            io::messages.add("Atom number out of range in BOND block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in BOND block: bond type < 1",
                    "In_Topology", io::message::error);
          }

          if (param.constraint.ntc == 3) {
            topo.solute().distance_constraints().
                    push_back(topology::two_body_term_struct(i - 1, j - 1, t - 1));
          } else
            topo.solute().bonds().
            push_back(topology::two_body_term_struct(i - 1, j - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in BOND block",
                  "In_Topology", io::message::error);
        }
      }

      if (!quiet)
        os << "\n\tEND\n";

} // BOND

void io::In_Topology::read_block_BONDDP(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) { // BONDDP
      DEBUG(10, "BONDDP block");

      if (!quiet)
        os << "\tBONDDP";

      std::vector<std::string> buffer = m_block["BONDDP"];
      if (buffer.size()) {
        block_read.insert("BONDDP");
        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet) {
            os << "\n\t\tspecial bonds to dipole particles : "
                  << num;
        }

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, t = 0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in BONDDP block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1) {
            io::messages.add("Atom number out of range in BONDDP block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in BONDDP block: bond type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().cgbonds().
                  push_back(topology::two_body_term_struct(i - 1, j - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in BONDDP block",
                  "In_Topology", io::message::error);
        }
      }

      if (!quiet)
        os << "\n\tEND\n";

} // BONDDP

void io::In_Topology::read_block_CONSTRAINT(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)  { // CONSTRAINT
      DEBUG(10, "CONSTRAINT block");
      std::vector<std::string> buffer = m_block["CONSTRAINT"];

      if (buffer.size() && param.constraint.ntc > 1) {
        block_read.insert("CONSTRAINT");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\tCONSTRAINT\n\t\t"
                << num
                << " bonds in CONSTRAINT block."
                << "\n\t\ttotal of constraint bonds : "
                << num + unsigned(topo.solute().distance_constraints().size())
          << "\n\tEND\n";

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, t = 0;

          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in CONSTRAINT block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1) {
            io::messages.add("Atom number out of range in CONSTRAINT block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in CONSTRAINT block: bond type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().distance_constraints().
                  push_back(topology::two_body_term_struct(i - 1, j - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in CONSTRAINT block",
                  "In_Topology", io::message::error);
        }
      } else if (buffer.size()) {
        block_read.insert("CONSTRAINT");
      }

} // CONSTRAINT

void io::In_Topology::read_block_BONDANGLEH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)   { // BONDANGLEH

      if (!quiet)
        os << "\tBONDANGLE";

      DEBUG(10, "BONDANGLEH block");
      std::vector<std::string> buffer = m_block["BONDANGLEH"];

      if (buffer.size()) {
        block_read.insert("BONDANGLEH");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\n\t\tbondangles containing hydrogens : " << num;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, k = 0, t = 0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> k >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in BONDANGLEH block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  k > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1 || k < 1) {
            io::messages.add("Atom number out of range in BONDANGLEH block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in BONDANGLE block: bond angle type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().angles().
                  push_back(topology::three_body_term_struct(i - 1, j - 1, k - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in BONDANGLEH block",
                  "In_Topology", io::message::error);
        }
      }

} // BONDANGLEH

void io::In_Topology::read_block_BONDANGLE(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) { // BONDANGLE
      DEBUG(10, "BONDANGLE block");
      std::vector<std::string> buffer = m_block["BONDANGLE"];

      if (buffer.size()) {
        block_read.insert("BONDANGLE");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\n\t\tbondangles not containing hydrogens : " << num;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, k = 0, t = 0;

          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> k >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in BONDANGLE block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  k > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1 || k < 1) {
            io::messages.add("Atom number out of range in BONDANGLE block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in BONDANGLEH block: bond angle type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().angles().
                  push_back(topology::three_body_term_struct(i - 1, j - 1, k - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in BONDANGLE block",
                  "In_Topology", io::message::error);
        }
      }

      if (!quiet)
        os << "\n\tEND\n";

} // BONDANGLE


void io::In_Topology::read_block_IMPDIHEDRAL(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) { // IMPDIHEDRAL
      DEBUG(10, "IMPDIHEDRAL block");
      std::vector<std::string> buffer = m_block["IMPDIHEDRAL"];

      if (!quiet)
        os << "\tIMPDIHEDRAL";

      if (buffer.size()) {
        block_read.insert("IMPDIHEDRAL");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\n\t\timproper dihedrals not containing hydrogens : "
                << num;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, k = 0, l = 0, t = 0;

          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> k >> l >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in IMPDIHEDRAL block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  k > int(topo.num_solute_atoms()) || l > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1 || k < 1 || l < 1) {
            io::messages.add("Atom number out of range in IMPDIHEDRAL block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in IMPDIHEDRAL block: improper dihedral type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().improper_dihedrals().
                  push_back(topology::four_body_term_struct(i - 1, j - 1, k - 1, l - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in IMPDIHEDRAL block",
                  "In_Topology", io::message::error);
        }
      }

    } // IMPDIHEDRAL
void io::In_Topology::read_block_IMPDIHEDRALH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) { // IMPDIHEDRALH
      DEBUG(10, "IMPDIHEDRALH block");
      std::vector<std::string> buffer = m_block["IMPDIHEDRALH"];

      if (buffer.size()) {
        block_read.insert("IMPDIHEDRALH");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\n\t\timproper dihedrals containing hydrogens : "
                << num;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, k = 0, l = 0, t = 0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> k >> l >> t;


          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in IMPDIHEDRALH block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  k > int(topo.num_solute_atoms()) || l > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1 || k < 1 || l < 1) {
            io::messages.add("Atom number out of range in IMPDIHEDRALH block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in IMPDIHEDRALH block: improper dihedral type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().improper_dihedrals().
                  push_back(topology::four_body_term_struct(i - 1, j - 1, k - 1, l - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in IMPDIHEDRALH block",
                  "In_Topology", io::message::error);
        }
      }

      if (!quiet)
        os << "\n\tEND\n";

    } // IMPDIHEDRALH
void io::In_Topology::read_block_DIHEDRAL(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) { // DIHEDRAL
      DEBUG(10, "DIHEDRAL block");
      std::vector<std::string> buffer = m_block["DIHEDRAL"];

      if (!quiet)
        os << "\tDIHEDRAL";

      if (buffer.size()) {
        block_read.insert("DIHEDRAL");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\n\t\tdihedrals not containing hydrogens : "
                << num;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, k = 0, l = 0, t = 0;

          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> k >> l >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in DIHEDRAL block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  k > int(topo.num_solute_atoms()) || l > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1 || k < 1 || l < 1) {
            io::messages.add("Atom number out of range in DIHEDRAL block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in DIHEDRAL block: dihedral type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().dihedrals().
                  push_back(topology::four_body_term_struct(i - 1, j - 1, k - 1, l - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in DIHEDRAL block",
                  "In_Topology", io::message::error);
        }
      }

    } // DIHEDRAL
void io::In_Topology::read_block_DIHEDRALH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)  { // DIHEDRALH
      DEBUG(10, "DIHEDRALH block");

      std::vector<std::string> buffer = m_block["DIHEDRALH"];
      if (buffer.size()) {
        block_read.insert("DIHEDRALH");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\n\t\tdihedrals containing hydrogens : "
                << num;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, k = 0, l = 0, t = 0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> k >> l >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in DIHEDRALH block",
                    "In_Topology", io::message::error);
          }

          if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
                  k > int(topo.num_solute_atoms()) || l > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1 || k < 1 || l < 1) {
            io::messages.add("Atom number out of range in DIHEDRALH block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in DIHEDRALH block: dihedral type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().dihedrals().
                  push_back(topology::four_body_term_struct(i - 1, j - 1, k - 1, l - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in DIHEDRALH block",
                  "In_Topology", io::message::error);
        }
      }

      if (!quiet)
        os << "\n\tEND\n";

    } // DIHEDRALH
void io::In_Topology::read_block_CROSSDIHEDRAL(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)  { // CROSSDIHEDRAL
      DEBUG(10, "CROSSDIHEDRAL block");
      std::vector<std::string> buffer = m_block["CROSSDIHEDRAL"];

      if (!quiet)
        os << "\tCROSSDIHEDRAL";

      if (buffer.size()) {
        block_read.insert("CROSSDIHEDRAL");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\n\t\tcrossdihedrals not containing hydrogens : "
                << num;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, t = 0;

          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> a >> b >> c >> d >> e >> f >> g >> h >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in CROSSDIHEDRAL block",
                    "In_Topology", io::message::error);
          }

          if (a > int(topo.num_solute_atoms()) || b > int(topo.num_solute_atoms()) ||
                  c > int(topo.num_solute_atoms()) || d > int(topo.num_solute_atoms()) ||
                  e > int(topo.num_solute_atoms()) || f > int(topo.num_solute_atoms()) ||
                  g > int(topo.num_solute_atoms()) || h > int(topo.num_solute_atoms()) ||
                  a < 1 || b < 1 || c < 1 || d < 1 || e < 1 || f < 1 || g < 1 || h < 1) {
            io::messages.add("Atom number out of range in CROSSDIHEDRAL block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in CROSSDIHEDRAL block: cross dihedral type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().crossdihedrals().
                  push_back(topology::eight_body_term_struct(a - 1, b - 1, c - 1, d - 1, e - 1,
                  f - 1, g - 1, h - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in CROSSDIHEDRAL block",
                  "In_Topology", io::message::error);
        }
      }

    } // CROSSDIHEDRAL
void io::In_Topology::read_block_CROSSDIHEDRALH(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) { // CROSSDIHEDRALH
      DEBUG(10, "CROSSDIHEDRALH block");

      std::vector<std::string> buffer = m_block["CROSSDIHEDRALH"];
      if (buffer.size()) {
        block_read.insert("CROSSDIHEDRALH");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\n\t\tcrossdihedrals containing hydrogens : "
                << num;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, t = 0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> a >> b >> c >> d >> e >> f >> g >> h >> t;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in CROSSDIHEDRALH block",
                    "In_Topology", io::message::error);
          }

          if (a > int(topo.num_solute_atoms()) || b > int(topo.num_solute_atoms()) ||
                  c > int(topo.num_solute_atoms()) || d > int(topo.num_solute_atoms()) ||
                  e > int(topo.num_solute_atoms()) || f > int(topo.num_solute_atoms()) ||
                  g > int(topo.num_solute_atoms()) || h > int(topo.num_solute_atoms()) ||
                  a < 1 || b < 1 || c < 1 || d < 1 || e < 1 || f < 1 || g < 1 || h < 1) {
            io::messages.add("Atom number out of range in CROSSDIHEDRALH block",
                    "In_Topology", io::message::error);
          }

          if (t < 1) {
            io::messages.add("Error in CROSSDIHEDRALH block: cross dihedral type < 1",
                    "In_Topology", io::message::error);
          }

          topo.solute().crossdihedrals().
                  push_back(topology::eight_body_term_struct(a - 1, b - 1, c - 1, d - 1, e - 1,
                  f - 1, g - 1, h - 1, t - 1));
        }

        if (n != num) {
          io::messages.add("Wrong number of bonds in CROSSDIHEDRALH block",
                  "In_Topology", io::message::error);
        }
      }

      if (!quiet)
        os << "\n\tEND\n";

} // CROSSDIHEDRALH

void io::In_Topology::read_block_VIRTUALGRAIN(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) { // VIRTUALGRAIN
      DEBUG(10, "VIRTUALGRAIN block");

      std::vector<std::string> buffer = m_block["VIRTUALGRAIN"];
      if (buffer.size()) {
        block_read.insert("VIRTUALGRAIN");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\tVIRTUALGRAIN\n\t\tVirtual Grains : "
                << num;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          _lineStream.clear();
          _lineStream.str(*it);

          // number of real atoms to define virtual atom
          int index = 0, i = 0, q = 0;
          _lineStream >> index >> i;

          std::vector<int> cog;

          for (int j = 0; j < i; ++j) {
            _lineStream >> q;
            cog.push_back(q - 1);
          }

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in VIRTUALGRAIN block",
                    "In_Topology", io::message::error);
          }

          topo.virtual_grains().
                  push_back(topology::virtual_grain_struct(index - 1,
                  util::Virtual_Atom(util::va_cog, cog)));
        }

        if (n != num) {
          io::messages.add("Wrong number of elements in VIRTUALGRAIN block",
                  "In_Topology", io::message::error);
        }

        if (!quiet)
          os << "\n\tEND\n";
      }

} // VIRTUALGRAIN

void io::In_Topology::read_block_SOLUTEMOLECULES(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) { // SOLUTEMOLECULES
      DEBUG(10, "read SOLUTEMOLECULES");

      std::string s;

      std::vector<std::string> buffer = m_block["SOLUTEMOLECULES"];

      if (!buffer.size()) {
        io::messages.add("no SOLUTEMOLECULES block in topology",
                "In_Topology", io::message::error);
      } else {
        block_read.insert("SOLUTEMOLECULES");
        _lineStream.clear();
        _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

        int num = 0;

        _lineStream >> num;

        if (num < 0) {
          io::messages.add("negative number of SOLUTEMOLECULES is not allowed",
                  "In_Topology", io::message::error);
        }

        unsigned int m = 0;
        unsigned int old_m = 0;

        for (int i = 0; i < num; ++i) {
          _lineStream >> m;
          topo.molecules().push_back(m);
          DEBUG(11, "add submol " << m);
          if (m < old_m) {
            io::messages.add("wrong order in SOLUTEMOLECULES block",
                    "In_Topology", io::message::error);
            break;
          }
          old_m = m;
        }

        if (_lineStream.fail())
          io::messages.add("bad line in SOLUTEMOLECULES block",
                "In_Topology", io::message::error);
      }
    } // SOLUTEMOLECULES
void io::In_Topology::read_block_TEMPERATUREGROUPS(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)  { // TEMPERATUREGROUPS
      DEBUG(10, "read TEMPERATUREGROUPS");

      std::string s;

      std::vector<std::string> buffer = m_block["TEMPERATUREGROUPS"];

      if (!buffer.size()) {
        io::messages.add("no TEMPERATUREGROUPS block in topology",
                "In_Topology", io::message::error);
      } else {
        block_read.insert("TEMPERATUREGROUPS");
        _lineStream.clear();
        _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

        int num = 0;

        _lineStream >> num;

        if (num < 0) {
          io::messages.add("negative number of TEMPERATUREGROUPS is not allowed",
                  "In_Topology", io::message::error);
        }

        unsigned int m = 0;
        unsigned int old_m = 0;

        for (int i = 0; i < num; ++i) {
          _lineStream >> m;
          topo.temperature_groups().push_back(m);
          DEBUG(11, "add temperature group " << m);
          if (m < old_m) {
            io::messages.add("wrong order in TEMPERATUREGROUPS block",
                    "In_Topology", io::message::error);
            break;
          }
          old_m = m;
        }

        if (_lineStream.fail())
          io::messages.add("bad line in TEMPERATUREGROUPS block",
                "In_Topology", io::message::error);
      }
    } // TEMPERATUREGROUPS

void io::In_Topology::read_block_PRESSUREGROUPS(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os)  { // PRESSUREGROUPS
      DEBUG(10, "read PRESSUREGROUPS");

      std::string s;

      std::vector<std::string> buffer = m_block["PRESSUREGROUPS"];

      if (!buffer.size()) {
        io::messages.add("no PRESSUREGROUPS block in topology",
                "In_Topology", io::message::error);
      } else {
        block_read.insert("PRESSUREGROUPS");
        _lineStream.clear();
        _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

        int num = 0;

        _lineStream >> num;

        if (num < 0) {
          io::messages.add("negative number of PRESSUREGROUPS is not allowed",
                  "In_Topology", io::message::error);
        }

        unsigned int m = 0;
        unsigned int old_m = 0;

        for (int i = 0; i < num; ++i) {
          _lineStream >> m;
          topo.pressure_groups().push_back(m);
          DEBUG(11, "add pressure group " << m);
          if (m < old_m) {
            io::messages.add("wrong order in PRESSUREGROUPS block",
                    "In_Topology", io::message::error);
            break;
          }
          old_m = m;
        }

        if (_lineStream.fail())
          io::messages.add("bad line in PRESSUREGROUPS block",
                "In_Topology", io::message::error);
      }
} // PRESSUREGROUPS



void io::In_Topology::read_SOLVENT_blocks(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os) {
    DEBUG(10, "SOLVENTATOM block");
    std::vector<std::string> buffer = m_block["SOLVENTATOM"];

    if (!quiet)
      os << "\tSOLVENT";

    if (buffer.size()) {
        block_read.insert("SOLVENTATOM");


      unsigned int res_nr = unsigned(topo.residue_names().size());

      topo.residue_names().push_back("SOLV");

      std::vector<std::string>::const_iterator it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num = 0, n = 0;
      _lineStream >> num;
      ++it;

      if (!quiet)
        os << "\n\t\tatoms : " << num;

      topology::Solvent s;

      std::string name;
      int i = 0, iac = 0;
      double mass = 0.0, charge = 0.0;
      //adiabatic decoupling
      int adc_index = 0;
      double sm = 1;
      int total_nr = topo.num_solute_atoms() - 1;


      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        _lineStream.clear();
        _lineStream.str(*it);
        sm = 1;
        total_nr++;
        adc_index = param.addecouple.check_index_adc(total_nr);
        if (adc_index != -1) {
          sm = param.addecouple.adc_index()[adc_index].sm;
        }
        _lineStream >> i >> name >> iac >> mass >> charge;
        mass *= sm;

        if (_lineStream.fail() || !_lineStream.eof()) {
          io::messages.add("Bad line in SOLVENTATOM block",
                  "In_Topology", io::message::error);
        }

        if (name.length() > MAX_NAME) {
          std::ostringstream msg;
          msg << "Error in SOLVENTATOM block: name " << name
                  << " is too long (> " << MAX_NAME << " characters).";
          io::messages.add(msg.str(), "InTopology", io::message::error);
        }
        // IAC indices in topology file start from "1" to match documentation
        // IAC indices in Gromos C++ standard start from "0" to match C++ data structure layouts
        s.add_atom(name, res_nr, iac - 1, mass, charge);
      }

      if (n != num) {
        io::messages.add("Error in SOLVENTATOM block (num != n)",
                "In_Topology", io::message::error);
      }
      read_block_SOLVENTPOLARISATION(topo, param, os, s);
      read_block_SOLVENTCONSTR(topo, param, os, s);

      topo.add_solvent(s);

    } else {
      io::messages.add("no solvent topology specified",
              "In_Topology", io::message::warning);
    }

} // end SOLVENTATOM

void io::In_Topology::read_block_SOLVENTPOLARISATION(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os, topology::Solvent &s) {
      DEBUG(10, "SOLVENTPOLARISATION block");

      std::vector<std::string> buffer = m_block["SOLVENTPOLARISATION"];
      if (buffer.size()) {
        block_read.insert("SOLVENTPOLARISATION");
        if (!quiet)
          os << "\n\t\tpolarisation parameters present";
        std::vector<std::string>::const_iterator it = buffer.begin() + 1;

        _lineStream.clear();
        _lineStream.str(*it);

        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          int i = 0, j = 0, k = 0;
          double polarisability = 0.0, coscharge = 0.0, damping_level = 0.0, damping_power = 0.0, gamma = 0.0;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> polarisability >> coscharge >> damping_level
                  >> damping_power >> gamma >> j >> k;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in SOLVENTPOLARISATION block",
                    "In_Topology", io::message::error);
          }

          if (i > int(s.num_atoms()) || i < 1) {
            io::messages.add("Atom number out of range in SOLVENTPOLARISATION block",
                    "In_Topology", io::message::error);
          } else if (gamma != 0.0 && (j < 1 || k < 1 || j > int(s.num_atoms())
                  || k > int(s.num_atoms()) || i == k || i == j || k == j)) {
            io::messages.add("Atom number for off atom out of range in SOLVENTPOLARISATION block",
                    "In_Topology", io::message::error);
          } else {
            DEBUG(10, "\tpolarisable atom: " << i);

            s.atoms()[i - 1].polarisability = polarisability / math::four_pi_eps_i;
            s.atoms()[i - 1].coscharge = coscharge;
            s.atoms()[i - 1].damping_level = damping_level * sqrt(math::four_pi_eps_i);
            s.atoms()[i - 1].damping_power = damping_power;
            if (param.polarise.cos == 2) {
              s.atoms()[i - 1].gamma = 2 * gamma;
              s.atoms()[i - 1].gamma_j = j - 1;
              s.atoms()[i - 1].gamma_k = k - 1;
            }

          }
        }

        if (n != num) {
          io::messages.add("Wrong number of polarisable atoms in SOLVENTPOLARISATION block",
                  "In_Topology", io::message::error);
        }
      }

} // SOLVENTPOLARISATION

void io::In_Topology::read_block_SOLVENTCONSTR(topology::Topology& topo,
        simulation::Parameter &param, std::ostream & os, topology::Solvent &s) {
      std::vector<std::string> buffer = m_block["SOLVENTCONSTR"];
      if (!buffer.size()) {
        io::messages.add("no SOLVENTCONST (block missing).",
                "In_Topology", io::message::notice);
      } else {
        block_read.insert("SOLVENTCONSTR");

        std::vector<std::string>::const_iterator it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);

        int i = 0, j = 0, n = 0, num = 0;
        _lineStream >> num;
        ++it;

        if (!quiet)
          os << "\n\t\tconstraints : " << num;

        double b0 = 0.0;

        for (n = 0; it != buffer.end() - 1; ++it, ++n) {
          _lineStream.clear();
          _lineStream.str(*it);

          _lineStream >> i >> j >> b0;

          if (_lineStream.fail() || !_lineStream.eof()) {
            io::messages.add("Bad line in SOLVENTCONSTR block",
                    "In_Topology", io::message::error);
          }

          // the solvent (distance constraints) bond types
          s.add_distance_constraint
                  (topology::two_body_term_struct(i - 1, j - 1, num_solute_bondtypes + n));
        }

        if (n != num) {
          io::messages.add("Error in SOLVENTCONSTR block (num != n)",
                  "In_Topology", io::message::error);
        }
      }
} // SOLVENTCONSTR
