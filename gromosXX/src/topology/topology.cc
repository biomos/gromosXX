/**
 * @file topology.cc
 * methods definition
 */

#include <set>

#include "../stdheader.h"

double topology_ver = 0.30;

namespace simulation {
  class Multibath;
}


#include "../util/debug.h"
#include "../topology/core/core.h"

#include "../topology/solute.h"
#include "../topology/solvent.h"
#include "../topology/perturbed_atom.h"
#include "../topology/perturbed_solute.h"

#include "../topology/topology.h"
#include "../simulation/multibath.h"
#include "../simulation/simulation.h"

#include "../interaction/interaction_types.h"
#include "../util/le_coordinate.h"

#include "../simulation/parameter.h"
#include <limits>

#undef MODULE
#undef SUBMODULE
#define MODULE topology
#define SUBMODULE topology

/**
 * Constructor
 */
topology::Topology::Topology()
: m_mass(0),
m_inverse_mass(0),
m_charge(0),
m_num_solute_chargegroups(0),
m_num_solute_molecules(0),
m_num_solute_temperature_groups(0),
m_num_solute_pressure_groups(0),
m_rottrans_last_atom(0),
m_multicell_topo(NULL),
m_polarisability(0),
m_coscharge(0),
m_damping_level(0),
m_damping_power(0),
m_gamma(0),
m_gamma_j(0),
m_gamma_k(0),
m_cg_factor(1),
m_tot_cg_factor(1.0),
m_le_coordinates(),
m_sasa_parameter(),
m_sasa_volume_tot(0.0),
m_sasa_first_neighbour(),
m_sasa_second_neighbour(),
m_sasa_third_neighbour(),
m_sasa_higher_neighbour(),
m_is_qm(0),
m_is_qm_buffer(0),
m_qm_delta_charge(0),
m_qm_atomic_number(0) {
  m_chargegroup.push_back(0);
  m_molecule.push_back(0);
  m_temperature_group.push_back(0);
  m_pressure_group.push_back(0);
}

topology::Topology::~Topology() {
  if (m_multicell_topo != NULL)
    delete m_multicell_topo;

  for (unsigned int i = 0; i < m_le_coordinates.size(); ++i)
    delete m_le_coordinates[i];
}

/**
 * copy (multiply) constructor
 */
topology::Topology::Topology(topology::Topology const & topo, int mul_solute, int mul_solvent)
: m_multicell_topo(NULL) {

  DEBUG(7, "multiplying topology");

  if (mul_solvent == -1) mul_solvent = mul_solute;

  const int num_solute = topo.num_solute_atoms();
  //assert(num_solute >= 0 && topo.m_is_perturbed.size() == unsigned(num_solute));

  m_bond_types_harm = topo.bond_types_harm();
  m_bond_types_quart = topo.bond_types_quart();
  m_angle_types_harm = topo.angle_types_harm();
  m_angle_types_cosharm = topo.angle_types_cosharm();
  m_dihedral_types = topo.dihedral_types();
  m_impdihedral_types = topo.impdihedral_types();
  m_virtual_atom_types = topo.virtual_atom_types();

// CHRIS: this seems to be copy-pasted piece of code?
//    std::vector<interaction::bond_type_struct> m_bond_types_quart;

    /**
     * store all available angle types with harmonic/cosine harmonic force constant
     */
//    std::vector<interaction::angle_type_struct> m_angle_types_harm;
//    std::vector<interaction::angle_type_struct> m_angle_types_cosharm;

    /**
     * store all available dihedral types
     */
//    std::vector<interaction::dihedral_type_struct> m_dihedral_types;

    /**
     * store all available improper dihedral types
     */
//    std::vector<interaction::improper_dihedral_type_struct> m_impdihedral_types;
  m_is_perturbed.clear();
  m_is_eds_perturbed.clear();
  m_gamd_accel_group.clear();
  m_gamd_interaction_pairs.clear();
  m_is_polarisable.clear();
  m_is_coarse_grained.clear();
  m_iac.clear();
  m_mass.clear();
  m_inverse_mass.clear();
  m_charge.clear();
  m_polarisability.clear();
  m_coscharge.clear();
  m_damping_level.clear();
  m_damping_power.clear();
  m_gamma.clear();
  m_gamma_j.clear();
  m_gamma_k.clear();
  m_cg_factor.clear();
  m_exclusion.clear();
  m_one_four_pair.clear();
  m_all_exclusion.clear();
  m_atom_energy_group.clear();
  m_chargegroup.clear();
  m_chargegroup.push_back(0);
  m_stochastic.clear();
  m_lj_exceptions.clear();
  // QMMM TEST
  m_qm_exclusion.clear();;
  m_qm_one_four_pair.clear();;
  m_qm_all_exclusion.clear();;
  m_qm_lj_exceptions.clear();;
  // END QMMM TEST
  m_is_qm.clear();
  m_is_qm_buffer.clear();
  m_qm_delta_charge.clear();
  m_qm_atomic_number.clear();
  m_qmmm_link.clear();

  DEBUG(10, "solute chargegroups = " << topo.num_solute_chargegroups());

  m_num_solute_chargegroups = topo.num_solute_chargegroups() * mul_solute;
  m_num_solute_molecules = topo.num_solute_molecules() * mul_solute;
  m_num_solute_temperature_groups = topo.num_solute_temperature_groups() * mul_solute;
  m_num_solute_pressure_groups = topo.num_solute_pressure_groups() * mul_solute;

  m_molecule.clear();
  m_molecule.push_back(0);

  m_temperature_group.clear();
  m_temperature_group.push_back(0);

  m_pressure_group.clear();
  m_pressure_group.push_back(0);

  solute().bonds().clear();
  solute().angles().clear();
  solute().improper_dihedrals().clear();
  solute().dihedrals().clear();
  solute().crossdihedrals().clear();

  DEBUG(8, "\tmultiplying solute");

  for (int m = 0; m < mul_solute; ++m) {
    DEBUG(10, "\tmul " << m);

    m_is_perturbed.insert(m_is_perturbed.end(), topo.m_is_perturbed.begin(), topo.m_is_perturbed.end());
    m_is_eds_perturbed.insert(m_is_eds_perturbed.end(), topo.m_is_eds_perturbed.begin(), topo.m_is_eds_perturbed.end());
    m_gamd_accel_group.insert(m_gamd_accel_group.end(), topo.m_gamd_accel_group.begin(), topo.m_gamd_accel_group.end());
    m_is_polarisable.insert(m_is_polarisable.end(), topo.m_is_polarisable.begin(), topo.m_is_polarisable.end());
    m_is_coarse_grained.insert(m_is_coarse_grained.end(), topo.m_is_coarse_grained.begin(), topo.m_is_coarse_grained.end());
    m_qmmm_link.insert(topo.m_qmmm_link.begin(), topo.m_qmmm_link.end());

    for (int i = 0; i < num_solute; ++i) {
      m_iac.push_back(topo.m_iac[i]);
      m_mass.push_back(topo.m_mass[i]);
      m_inverse_mass.push_back(topo.m_inverse_mass[i]);
      m_charge.push_back(topo.m_charge[i]);
      m_polarisability.push_back(topo.m_polarisability[i]);
      m_coscharge.push_back(topo.m_coscharge[i]);
      m_damping_level.push_back(topo.m_damping_level[i]);
      m_damping_power.push_back(topo.m_damping_power[i]);
      m_gamma.push_back(topo.m_gamma[i]);
      m_gamma_j.push_back(topo.m_gamma_j[i]);
      m_gamma_k.push_back(topo.m_gamma_k[i]);
      m_cg_factor.push_back(topo.m_cg_factor[i]);
      m_is_qm.push_back(topo.m_is_qm[i]);
      m_is_qm_buffer.push_back(topo.m_is_qm_buffer[i]);
      m_qm_delta_charge.push_back(topo.m_qm_delta_charge[i]);
      m_qm_atomic_number.push_back(topo.m_qm_atomic_number[i]);

      topology::excl_cont_t::value_type ex;
      topology::excl_cont_t::value_type::const_iterator it = topo.m_exclusion[i].begin(),
              to = topo.m_exclusion[i].end();
      for (; it != to; ++it)
        ex.insert(*it + num_solute * m);
      m_exclusion.push_back(ex);

      ex.clear();
      it = topo.m_one_four_pair[i].begin();
      to = topo.m_one_four_pair[i].end();
      for (; it != to; ++it)
        ex.insert(*it + num_solute * m);
      m_one_four_pair.push_back(ex);

      ex.clear();
      it = topo.m_all_exclusion[i].begin();
      to = topo.m_all_exclusion[i].end();
      for (; it != to; ++it)
        ex.insert(*it + num_solute * m);
      m_all_exclusion.push_back(ex);

      std::vector<topology::lj_exception_struct>::const_iterator ljex_it =
              topo.lj_exceptions().begin(), ljex_to =
              topo.lj_exceptions().end();
      for (; ljex_it != ljex_to; ++ljex_it) {
        m_lj_exceptions.push_back(
                topology::lj_exception_struct(ljex_it->i + num_solute * m, ljex_it->j + num_solute * m,
                ljex_it->c6, ljex_it->c12));
      }

      // QMMM TEST
      {
        topology::excl_cont_t::value_type ex;
        topology::excl_cont_t::value_type::const_iterator it = topo.m_qm_exclusion[i].begin(),
                to = topo.m_qm_exclusion[i].end();
        for (; it != to; ++it)
          ex.insert(*it + num_solute * m);
        m_qm_exclusion.push_back(ex);

        ex.clear();
        it = topo.m_qm_one_four_pair[i].begin();
        to = topo.m_qm_one_four_pair[i].end();
        for (; it != to; ++it)
          ex.insert(*it + num_solute * m);
        m_qm_one_four_pair.push_back(ex);

        ex.clear();
        it = topo.m_qm_all_exclusion[i].begin();
        to = topo.m_qm_all_exclusion[i].end();
        for (; it != to; ++it)
          ex.insert(*it + num_solute * m);
        m_qm_all_exclusion.push_back(ex);

        std::vector<topology::lj_exception_struct>::const_iterator ljex_it =
                topo.qm_lj_exceptions().begin(), ljex_to =
                topo.qm_lj_exceptions().end();
        for (; ljex_it != ljex_to; ++ljex_it) {
          m_qm_lj_exceptions.push_back(
                  topology::lj_exception_struct(ljex_it->i + num_solute * m, ljex_it->j + num_solute * m,
                  ljex_it->c6, ljex_it->c12));
        }
      }
      // END QMMM TEST





      m_atom_energy_group.push_back(topo.m_atom_energy_group[i]);

      m_stochastic.gamma.push_back(topo.stochastic().gamma[i]);
      m_stochastic.c1.push_back(topo.stochastic().c1[i]);
      m_stochastic.c2.push_back(topo.stochastic().c2[i]);
      m_stochastic.c3.push_back(topo.stochastic().c3[i]);
      m_stochastic.c4.push_back(topo.stochastic().c4[i]);
      m_stochastic.c5.push_back(topo.stochastic().c5[i]);
      m_stochastic.c6.push_back(topo.stochastic().c6[i]);
      m_stochastic.c7.push_back(topo.stochastic().c7[i]);
      m_stochastic.c8.push_back(topo.stochastic().c8[i]);
      m_stochastic.c9.push_back(topo.stochastic().c9[i]);

      solute().add_atom(topo.solute().atom(i).name, topo.solute().atom(i).residue_nr);

    }

    DEBUG(10, "\tcg");
    for (unsigned int i = 1; i <= topo.num_solute_chargegroups(); ++i) {
      m_chargegroup.push_back(topo.m_chargegroup[i] + m * num_solute);
    }

    DEBUG(10, "\tmol");
    for (unsigned int i = 1; i < topo.molecules().size() - topo.num_solvent_molecules(0); ++i) {
      m_molecule.push_back(topo.molecules()[i] + m * num_solute);
    }

    DEBUG(10, "\ttemperature groups");
    for (unsigned int i = 1; i < topo.temperature_groups().size() - topo.num_solvent_molecules(0); ++i) {
      m_temperature_group.push_back(topo.temperature_groups()[i] + m * num_solute);
    }

    DEBUG(10, "\tpressure groups");
    for (unsigned int i = 1; i < topo.pressure_groups().size() - topo.num_solvent_molecules(0); ++i) {
      m_pressure_group.push_back(topo.pressure_groups()[i] + m * num_solute);
    }

    DEBUG(10, "\tbonds");
    for (unsigned int i = 0; i < topo.solute().bonds().size(); ++i) {
      solute().bonds().push_back
              (two_body_term_struct(topo.solute().bonds()[i].i + m * num_solute,
              topo.solute().bonds()[i].j + m * num_solute,
              topo.solute().bonds()[i].type));
    }

    DEBUG(10, "\tdistance constraints");
    for (unsigned int i = 0; i < topo.solute().distance_constraints().size(); ++i) {
      solute().add_distance_constraint
              (two_body_term_struct(topo.solute().distance_constraints()[i].i + m * num_solute,
              topo.solute().distance_constraints()[i].j + m * num_solute,
              topo.solute().distance_constraints()[i].type));
    }

    DEBUG(10, "\tangles");
    for (unsigned int i = 0; i < topo.solute().angles().size(); ++i) {
      solute().angles().push_back
              (three_body_term_struct(topo.solute().angles()[i].i + m * num_solute,
              topo.solute().angles()[i].j + m * num_solute,
              topo.solute().angles()[i].k + m * num_solute,
              topo.solute().angles()[i].type));
    }

    DEBUG(10, "\timps");
    for (unsigned int i = 0; i < topo.solute().improper_dihedrals().size(); ++i) {
      solute().improper_dihedrals().push_back
              (four_body_term_struct(topo.solute().improper_dihedrals()[i].i + m * num_solute,
              topo.solute().improper_dihedrals()[i].j + m * num_solute,
              topo.solute().improper_dihedrals()[i].k + m * num_solute,
              topo.solute().improper_dihedrals()[i].l + m * num_solute,
              topo.solute().improper_dihedrals()[i].type));
    }

    DEBUG(10, "\tdihedrals");
    for (unsigned int i = 0; i < topo.solute().dihedrals().size(); ++i) {
      solute().dihedrals().push_back
              (four_body_term_struct(topo.solute().dihedrals()[i].i + m * num_solute,
              topo.solute().dihedrals()[i].j + m * num_solute,
              topo.solute().dihedrals()[i].k + m * num_solute,
              topo.solute().dihedrals()[i].l + m * num_solute,
              topo.solute().dihedrals()[i].type));
    }
    DEBUG(10, "\tcrossdihedrals");
    for (unsigned int i = 0; i < topo.solute().crossdihedrals().size(); ++i) {
      solute().crossdihedrals().push_back
              (eight_body_term_struct(topo.solute().crossdihedrals()[i].a + m * num_solute,
              topo.solute().crossdihedrals()[i].b + m * num_solute,
              topo.solute().crossdihedrals()[i].c + m * num_solute,
              topo.solute().crossdihedrals()[i].d + m * num_solute,
              topo.solute().crossdihedrals()[i].e + m * num_solute,
              topo.solute().crossdihedrals()[i].f + m * num_solute,
              topo.solute().crossdihedrals()[i].g + m * num_solute,
              topo.solute().crossdihedrals()[i].h + m * num_solute,
              topo.solute().crossdihedrals()[i].type));
    }

    // perturbed solute
    DEBUG(10, "\tmultiplying perturbed solute...");
    DEBUG(10, "\tperturbed solute: atoms");
    std::map<unsigned int, Perturbed_Atom>::const_iterator it, to;
    it = topo.perturbed_solute().atoms().begin();
    to = topo.perturbed_solute().atoms().end();
    for (; it != to; ++it) {
      unsigned int seq = it->first + m * num_solute;
      perturbed_solute().atoms()[seq] = it->second;
      perturbed_solute().atoms()[seq].sequence_number(seq);

      // and the exlusions and one-four pairs... again!
      perturbed_solute().atoms()[seq].exclusion().clear();
      perturbed_solute().atoms()[seq].one_four_pair().clear();

      topology::excl_cont_t::value_type::const_iterator ex_it, ex_to;
      ex_it = it->second.exclusion().begin();
      ex_to = it->second.exclusion().end();
      for (; ex_it != ex_to; ++ex_it)
        perturbed_solute().atoms()[seq].exclusion().insert(*ex_it + m * num_solute);

      ex_it = it->second.one_four_pair().begin();
      ex_to = it->second.one_four_pair().end();
      for (; ex_it != ex_to; ++ex_it)
        perturbed_solute().atoms()[seq].one_four_pair().insert(*ex_it + m * num_solute);

    }

    DEBUG(10, "\tperturbed solute: bonds");
    for (unsigned int i = 0; i < topo.perturbed_solute().bonds().size(); ++i) {
      perturbed_solute().bonds().push_back
              (perturbed_two_body_term_struct(
              topo.perturbed_solute().bonds()[i].i + m * num_solute,
              topo.perturbed_solute().bonds()[i].j + m * num_solute,
              topo.perturbed_solute().bonds()[i].A_type,
              topo.perturbed_solute().bonds()[i].B_type));
    }

    DEBUG(10, "\tperturbed solute: atom pairs");
    for (unsigned int i = 0; i < topo.perturbed_solute().atompairs().size(); ++i) {
      perturbed_solute().atompairs().push_back
              (perturbed_two_body_term_struct(
              topo.perturbed_solute().atompairs()[i].i + m * num_solute,
              topo.perturbed_solute().atompairs()[i].j + m * num_solute,
              topo.perturbed_solute().atompairs()[i].A_type,
              topo.perturbed_solute().atompairs()[i].B_type));
    }

    DEBUG(10, "\tperturbed solute: distance constraints");
    for (unsigned int i = 0; i < topo.perturbed_solute().distance_constraints().size(); ++i) {
      perturbed_solute().distance_constraints().push_back
              (perturbed_two_body_term_struct(
              topo.perturbed_solute().distance_constraints()[i].i + m * num_solute,
              topo.perturbed_solute().distance_constraints()[i].j + m * num_solute,
              topo.perturbed_solute().distance_constraints()[i].A_type,
              topo.perturbed_solute().distance_constraints()[i].B_type));
    }

    DEBUG(10, "\tperturbed solute: angles");
    for (unsigned int i = 0; i < topo.perturbed_solute().angles().size(); ++i) {
      perturbed_solute().angles().push_back
              (perturbed_three_body_term_struct(topo.perturbed_solute().angles()[i].i + m * num_solute,
              topo.perturbed_solute().angles()[i].j + m * num_solute,
              topo.perturbed_solute().angles()[i].k + m * num_solute,
              topo.perturbed_solute().angles()[i].A_type,
              topo.perturbed_solute().angles()[i].B_type));
    }

    DEBUG(10, "\tperturbed solute: imps");
    for (unsigned int i = 0; i < topo.perturbed_solute().improper_dihedrals().size(); ++i) {
      perturbed_solute().improper_dihedrals().push_back
              (perturbed_four_body_term_struct(topo.perturbed_solute().improper_dihedrals()[i].i + m * num_solute,
              topo.perturbed_solute().improper_dihedrals()[i].j + m * num_solute,
              topo.perturbed_solute().improper_dihedrals()[i].k + m * num_solute,
              topo.perturbed_solute().improper_dihedrals()[i].l + m * num_solute,
              topo.perturbed_solute().improper_dihedrals()[i].A_type,
              topo.perturbed_solute().improper_dihedrals()[i].B_type));
    }

    DEBUG(10, "\tperturbed solute: dihedrals");
    for (unsigned int i = 0; i < topo.perturbed_solute().dihedrals().size(); ++i) {
      perturbed_solute().dihedrals().push_back
              (perturbed_four_body_term_struct(topo.perturbed_solute().dihedrals()[i].i + m * num_solute,
              topo.perturbed_solute().dihedrals()[i].j + m * num_solute,
              topo.perturbed_solute().dihedrals()[i].k + m * num_solute,
              topo.perturbed_solute().dihedrals()[i].l + m * num_solute,
              topo.perturbed_solute().dihedrals()[i].A_type,
              topo.perturbed_solute().dihedrals()[i].B_type));
    }

  }

  DEBUG(10, "\tegroups");
  m_energy_group = topo.m_energy_group;

  DEBUG(10, "\tperturbation param");
  m_lambda = topo.m_lambda;
  m_old_lambda = topo.m_old_lambda;
  m_lambda_exp = topo.m_lambda_exp;
  m_energy_group_scaling = topo.m_energy_group_scaling;
  //m_energy_group_lambdadep = topo.m_energy_group_lambdadep;
  m_individual_lambda_parameters = topo.m_individual_lambda_parameters;
  m_individual_lambda = topo.m_individual_lambda;
  m_individual_lambda_derivative = topo.m_individual_lambda_derivative;
  //m_lambda_prime = topo.m_lambda_prime;
  //m_lambda_prime_derivative = topo.m_lambda_prime_derivative;
  //m_perturbed_energy_derivative_alpha =   topo.m_perturbed_energy_derivative_alpha;

  DEBUG(10, "\tspecial");
  m_position_restraint = topo.m_position_restraint;
  m_jvalue_restraint = topo.m_jvalue_restraint;
  m_rdc_restraint = topo.m_rdc_restraint; 
  m_rottrans_last_atom = topo.m_rottrans_last_atom;
  for (unsigned int i = 0; i < topo.m_le_coordinates.size(); ++i) {
    m_le_coordinates.push_back(topo.m_le_coordinates[i]->clone());
  }

  DEBUG(8, "\tnum_solute_atoms()=" << num_solute_atoms());

  if (topo.num_solvent_atoms()) {
    // solvent
    DEBUG(8, "\tmultiplying solvent");

    assert(topo.num_solvents() == 1);
    add_solvent(topo.solvent(0));
    solvate(0, topo.num_solvent_atoms() * mul_solvent / topo.solvent(0).num_atoms());


    for (int m = 0; m < mul_solvent; ++m)
      for (unsigned int s = topo.num_solute_atoms(); s < topo.num_atoms(); ++s) {
        m_atom_energy_group.push_back(topo.m_atom_energy_group[s]);
      }


    DEBUG(8, "\tnum_atoms()=" << num_atoms());
  }

  m_atom_name = topo.m_atom_name;

  m_tot_cg_factor = topo.m_tot_cg_factor;
}

/**
 * set the capacity of solute atoms by resizeing
 * the apropriate arrays.
 */
void topology::Topology::resize(unsigned int const atoms) {
  DEBUG(8, "resizing topology for " << atoms << " atoms!");

  // if you want to shrink, first change num_atoms...
  assert(atoms >= num_solute_atoms() + num_solvent_atoms());

  // m_mass.resizeAndPreserve(atoms);
  m_mass.resize(atoms);
  m_inverse_mass.resize(atoms);
  // m_charge.resizeAndPreserve(atoms);
  m_charge.resize(atoms);
  m_polarisability.resize(atoms);
  m_coscharge.resize(atoms);
  m_damping_level.resize(atoms);
  m_damping_power.resize(atoms);
  m_gamma.resize(atoms);
  m_gamma_j.resize(atoms);
  m_gamma_k.resize(atoms);
  m_cg_factor.resize(atoms, 1);
  m_exclusion.resize(atoms);
  m_one_four_pair.resize(atoms);
  m_all_exclusion.resize(atoms);
  m_stochastic.resize(atoms);
  //QMMM TEST
  m_qm_exclusion.resize(atoms);
  m_qm_one_four_pair.resize(atoms);
  m_qm_all_exclusion.resize(atoms);
  // END QMMM TEST
  m_qm_atomic_number.resize(atoms, 0);
  m_is_qm.resize(atoms, 0);
  m_is_qm_buffer.resize(atoms, 0);
  m_qm_delta_charge.resize(atoms, 0);

  m_iac.resize(atoms);
  // chargegroups???
  m_is_perturbed.resize(atoms, false);
  m_is_eds_perturbed.resize(atoms, false);
  m_gamd_accel_group.resize(atoms, 0);
  m_is_polarisable.resize(atoms, false);
  m_is_coarse_grained.resize(atoms, false);

}

void topology::Topology::init(simulation::Simulation const & sim,
        std::ostream & os, bool quiet) {
  DEBUG(6, "Topology INIT");

  if (!m_molecule.size()) {
    m_molecule.push_back(0);
    m_molecule.push_back(num_solute_atoms());
  }

  if (!m_temperature_group.size()) {
    m_temperature_group.push_back(0);
    m_temperature_group.push_back(num_solute_atoms());
  }

  if (!m_pressure_group.size()) {
    m_pressure_group.push_back(0);
    m_pressure_group.push_back(num_solute_atoms());
  }

  if (!m_energy_group.size()) {
    m_energy_group.push_back(num_solute_atoms() - 1);
    for (unsigned int i = 0; i < num_solute_atoms(); ++i)
      m_atom_energy_group.push_back(0);
  }

  if (sim.param().multicell.multicell) {

    const int num = sim.param().multicell.x * sim.param().multicell.y * sim.param().multicell.z;
    simulation::Simulation multi_sim = sim;
    multi_sim.param().multicell.multicell = false;
    m_multicell_topo = new Topology(*this, num);
    m_multicell_topo->init(multi_sim, os, true);
    m_multicell_topo->check_state();

    if (!quiet) {
      os << "\tMULTICELL\n" << "\t\t"
              << std::setw(4) << "X"
              << std::setw(6) << sim.param().multicell.x
              << std::setw(4) << "Y"
              << std::setw(6) << sim.param().multicell.y
              << std::setw(4) << "Z"
              << std::setw(6) << sim.param().multicell.z
              << "\n\t\t"
              << "molecules : " << m_multicell_topo->molecules().size() - 1
              << "\n\t\t"
              << "atoms     : " << m_multicell_topo->num_atoms()
              << "\n\tEND\n";
    }
  }

  // add chargegroup exclusions (a clever technique to improve pairlisting...)
  update_chargegroup_exclusion();

  if (!sim.param().force.nonbonded_crf) {
    for (unsigned int i = 0; i < m_charge.size(); ++i) {
      m_charge(i) = 0.0;
    }
    for (std::map<unsigned int, Perturbed_Atom>::iterator it=perturbed_solute().atoms().begin(); it!=perturbed_solute().atoms().end(); ++it) {
      it->second.A_charge(0.0);
      it->second.B_charge(0.0);
    }
  }

  if (!sim.param().force.nonbonded_vdw) {
    if (m_atom_name.find("DUM") == m_atom_name.end()) {
      io::messages.add("no dummy atomtype (DUM) found in topology",
              "topology",
              io::message::error);
    } else {
      int dum = m_atom_name["DUM"];
      for (unsigned int i = 0; i < m_iac.size(); ++i)
        m_iac[i] = dum;
      for (std::map<unsigned int, Perturbed_Atom>::iterator it=perturbed_solute().atoms().begin(); it!=perturbed_solute().atoms().end(); ++it) {
        it->second.A_IAC(dum);
        it->second.B_IAC(dum);
      }
    }
  }

  // initialize the individual lambdas
  m_individual_lambda_parameters.a.resize(simulation::last_interaction_lambda);
  m_individual_lambda_parameters.b.resize(simulation::last_interaction_lambda);
  m_individual_lambda_parameters.c.resize(simulation::last_interaction_lambda);
  m_individual_lambda_parameters.d.resize(simulation::last_interaction_lambda);
  m_individual_lambda_parameters.e.resize(simulation::last_interaction_lambda);
  m_individual_lambda.resize(simulation::last_interaction_lambda);
  m_individual_lambda_derivative.resize(simulation::last_interaction_lambda);

  for (int i = 0; i < simulation::last_interaction_lambda; i++) {
    m_individual_lambda[i].resize(m_energy_group.size());
    m_individual_lambda_derivative[i].resize(m_energy_group.size());
    for (unsigned int n1 = 0; n1 < m_energy_group.size(); n1++) {
      m_individual_lambda[i][n1].resize(m_energy_group.size());
      m_individual_lambda_derivative[i][n1].resize(m_energy_group.size());

      for (unsigned int n2 = 0; n2 < m_energy_group.size(); n2++) {
        std::pair<int, int> p(n1, n2);

        assert(i < int(sim.param().lambdas.a.size()));
        assert(n1 < sim.param().lambdas.a[i].size());
        assert(n2 < sim.param().lambdas.a[i][n1].size());
        assert(i < int(m_individual_lambda_parameters.a.size()));

        m_individual_lambda_parameters.a[i][p] = sim.param().lambdas.a[i][n1][n2];
        m_individual_lambda_parameters.b[i][p] = sim.param().lambdas.b[i][n1][n2];
        m_individual_lambda_parameters.c[i][p] = sim.param().lambdas.c[i][n1][n2];
        m_individual_lambda_parameters.d[i][p] = sim.param().lambdas.d[i][n1][n2];
        m_individual_lambda_parameters.e[i][p] = sim.param().lambdas.e[i][n1][n2];

      }
    }
  }
  // And calculate the values for all individual lambdas and their derivatives
  update_for_lambda();

  // build neighbour lists for SASA implicit solvent model using POPS method
  if (sim.param().sasa.switch_sasa) {

    const unsigned int num_sasa_atoms = m_sasa_parameter.size();

    // resize sasa neighbours vectors
    m_sasa_first_neighbour.resize(num_sasa_atoms); // resize sasa neighbours vector
    m_sasa_second_neighbour.resize(num_sasa_atoms); // resize sasa neighbours vector
    m_sasa_third_neighbour.resize(num_sasa_atoms); // resize sasa neighbours vector
    m_sasa_higher_neighbour.resize(num_sasa_atoms); // resize sasa neighbours vector

    // define a bondnumber storage list
    std::vector<unsigned int> pathlength(num_sasa_atoms * num_sasa_atoms);
    // and fill with large numbers ("infinity")
    const unsigned int infinity = std::numeric_limits<unsigned int>::max();
    for (unsigned int i = 0; i < pathlength.size(); ++i) {
      pathlength[i] = infinity;
    }

    // loop over bonds
    std::vector<topology::two_body_term_struct>::const_iterator b_it, b_to;
    DEBUG(6, "building SASA neighbour lists from bond list");

    b_it = solute().bonds().begin();
    b_to = solute().bonds().end();
    for (; b_it != b_to; ++b_it) {
      DEBUG(6, "bonded atoms " << b_it->i << "\t" << b_it->j);
      // we have to do it this way because we want to store things by their
      // index in the list of sasa atoms
      for (unsigned int l = 0; l < num_sasa_atoms; ++l) {
        if (b_it->i == m_sasa_parameter[l].atom) { // l is our index for atom i
          for (unsigned int m = 0; m < num_sasa_atoms; ++m) {
            if (b_it->j == m_sasa_parameter[m].atom) { // m is our index for atom j
              // store first neighbour
              m_sasa_first_neighbour[l].insert(m);
              // store path between atoms i and j in both directions
              pathlength[l * num_sasa_atoms + m] = 1;
              pathlength[m * num_sasa_atoms + l] = 1;
            }
          }
        }
      }
    } // end loop over bonds (empty)

    // now consider constraints
    DEBUG(6, "building SASA first neighbour list from constraints");
    b_it = solute().distance_constraints().begin();
    b_to = solute().distance_constraints().end();
    for (; b_it != b_to; ++b_it) {
      DEBUG(6, "constrained atoms " << b_it->i << "\t" << b_it->j);
      // we have to do it this way because we want to store things by their
      // index in the list of sasa atoms
      for (unsigned int l = 0; l < num_sasa_atoms; ++l) {
        if (b_it->i == m_sasa_parameter[l].atom) { // l is our index for atom i
          for (unsigned int m = 0; m < num_sasa_atoms; ++m) {
            if (b_it->j == m_sasa_parameter[m].atom) { // m is our index for atom j
              // store first neighbour
              m_sasa_first_neighbour[l].insert(m);
              // store path between atoms i and j in both directions
              pathlength[l * num_sasa_atoms + m] = 1;
              pathlength[m * num_sasa_atoms + l] = 1;
            }
          }
        }
      }
    } // end loop over constraints


    // new loop over num_sasa_atoms (Floyds algorithm)
    for (unsigned int k = 0; k < num_sasa_atoms; ++k) {
      for (unsigned int i = 0; i < num_sasa_atoms - 1; ++i) {
        // stored path between i and k
        const unsigned int pathlength_ik = pathlength[i * num_sasa_atoms + k];
        if (pathlength_ik < infinity) {
          // stored path between k and j
          for (unsigned int j = i + 1; j < num_sasa_atoms; ++j) {
            const unsigned int pathlength_kj = pathlength[k * num_sasa_atoms + j];
            if (pathlength_kj < infinity) {
              // path between i and j via k
              const unsigned int complength = pathlength_ik + pathlength_kj;
              // compare to stored path between i and j
              if (complength < pathlength[i * num_sasa_atoms + j]) {
                pathlength[i * num_sasa_atoms + j] = complength;
              }
            }
          }
        }
      }
    } // end floyd's

    // now directly make the third and higher lists
    // c.f. POPS which makes exclusion and 1-4 lists and then neighbour lists
    for (unsigned int i = 0; i < num_sasa_atoms - 1; ++i) {
      for (unsigned int j = i + 1; j < num_sasa_atoms; ++j) {
        const unsigned int pathlength_ij = pathlength[i * num_sasa_atoms + j];
        if (pathlength_ij < infinity) {
          if (pathlength_ij == 3) {
            // store in 3rd neighbours
            // note we are storing the indexes not the actual atom numbers
            m_sasa_third_neighbour[i].insert(j);
          } else if (pathlength_ij > 3) {
            // store in higher neighbours
            // POPS takes only those closer than nonbonded cutoff but I do not
            // have access to the coords and thus inter-atomic distances here
            // also, the distances could change during the simulation
            // note we are storing the indexes not the actual atom numbers
            m_sasa_higher_neighbour[i].insert(j);
          }
        }
      }
    } // end third and higher

    // second neighbours
    std::vector<topology::three_body_term_struct>::const_iterator a_it, a_to;
    DEBUG(6, "building sasa second neighbour list from angles");
    a_it = solute().angles().begin();
    a_to = solute().angles().end();
    for (; a_it != a_to; ++a_it) {
      DEBUG(6, "second neighbours " << a_it->i << "\t" << a_it->k);
      for (unsigned int l = 0; l < num_sasa_atoms; ++l) {
        if (a_it->i == m_sasa_parameter[l].atom) { // l is our index for atom i
          for (unsigned int m = 0; m < num_sasa_atoms; ++m) {
            if (a_it->k == m_sasa_parameter[m].atom) { // m is our index for atom k
              m_sasa_second_neighbour[l].insert(m);
            }
          }
        }
      }
    } // end angles

    // initialise (total) surface areas, ri+rh2o and, if needed, volumes
    m_sasa_volume_tot = 0.0;
    for (unsigned int i = 0; i < num_sasa_atoms; ++i) {
      // total surface area
      const double surface = 4.0 * math::Pi * (m_sasa_parameter[i].r + sim.param().sasa.r_solv) *
              (m_sasa_parameter[i].r + sim.param().sasa.r_solv);
      DEBUG(15, "\tTotal surface of atom " << m_sasa_parameter[i].atom << " is " << surface);
      m_sasa_parameter[i].surface = surface;
      // ri+rh2o
      const double ri_rh2o = m_sasa_parameter[i].r + sim.param().sasa.r_solv;
      m_sasa_parameter[i].r_rh2o = ri_rh2o;
      // volumes
      if (sim.param().sasa.switch_volume) {
        const double volume = (4.0 / 3.0) * math::Pi * m_sasa_parameter[i].r *
                m_sasa_parameter[i].r * m_sasa_parameter[i].r;
        DEBUG(15, "\tTotal volume of atom " << m_sasa_parameter[i].atom << " is " << volume);
        m_sasa_parameter[i].vol = volume;
        m_sasa_volume_tot += volume;
      } // end volume
    } // end sasa atom loop
  } // end sasa

    // initialize the roto-translational constraints
    //std::cout << "@@ sim.param().rottrans.rottrans " << sim.param().rottrans.rottrans << "\n";
    if ( sim.param().rottrans.rottrans ) {
    //std::cout << "@@ sim.param().rottrans.last " << sim.param().rottrans.last << "\n";
    m_rottrans_last_atom = sim.param().rottrans.last;
    }

} // end init

/**
 * add a solute atom to the topology.
 * if the arrays are too small they will be increased.
 * if adding multiple solute atoms, first call solute_atoms_capacity...
 */
void topology::Topology
::add_solute_atom(std::string name, int residue_nr,
        int iac, double mass,
        double charge, bool chargegroup,
        topology::excl_cont_t::value_type exclusions,
        topology::excl_cont_t::value_type one_four_pairs) {

  if (unsigned(m_mass.size()) < num_solute_atoms() + 1) {
    resize(num_solute_atoms() + 1);
  }

  Topology::mass()(num_solute_atoms()) = mass;
  Topology::inverse_mass()(num_solute_atoms()) = 1.0 / mass;
  Topology::charge()(num_solute_atoms()) = charge;

  Topology::polarisability()(num_solute_atoms()) = 0.0;
  Topology::coscharge()(num_solute_atoms()) = 0.0;
  Topology::damping_level()(num_solute_atoms()) = 0.0;
  Topology::damping_power()(num_solute_atoms()) = 0.0;
  Topology::cg_factor()[num_solute_atoms()] = 1;

  if (chargegroup) {
    m_chargegroup.push_back(num_solute_atoms() + 1);
    ++m_num_solute_chargegroups;
  }

  DEBUG(15, "iac[" << num_solute_atoms() << "] = " << iac);

  m_iac[num_solute_atoms()] = iac;

  m_exclusion[num_solute_atoms()] = exclusions;
  m_one_four_pair[num_solute_atoms()] = one_four_pairs;

  std::set_union(exclusions.begin(), exclusions.end(),
          one_four_pairs.begin(), one_four_pairs.end(),
          std::inserter(m_all_exclusion[num_solute_atoms()],
          m_all_exclusion[num_solute_atoms()].end())
          );

  // this increases num_solute_atoms()
  solute().add_atom(name, residue_nr);
}

/**
 * add solvent molecules to the simulation (system).
 */
void topology::Topology::solvate(unsigned int solv, unsigned int num_molecules) {
  // only add in the correct order!
  assert(solv == m_num_solvent_atoms.size());
  assert(solv < m_solvent.size());

  int n = num_solute_atoms() + num_solvent_atoms();

  m_num_solvent_molecules.push_back(num_molecules);
  m_num_solvent_atoms.push_back(num_molecules * m_solvent[solv].num_atoms());

  DEBUG(5, "solvate: solvent atoms: " << num_solvent_atoms());
  DEBUG(10, "solvate: total atoms: " << num_solute_atoms() + num_solvent_atoms());

  resize(num_solute_atoms() + num_solvent_atoms());

  // add to iac, mass, charge
  for (unsigned int i = 0; i < num_molecules; ++i) {
    for (unsigned int j = 0; j < m_solvent[solv].num_atoms(); ++j, ++n) {

      DEBUG(15, "iac[" << n << "]=" << m_solvent[solv].atom(j).iac);
      DEBUG(15, "charge[" << n << "]=" << m_solvent[solv].atom(j).charge);

      m_iac[n] = m_solvent[solv].atom(j).iac;
      m_mass(n) = m_solvent[solv].atom(j).mass;
      m_inverse_mass(n) = 1.0 / m_solvent[solv].atom(j).mass;
      m_charge(n) = m_solvent[solv].atom(j).charge;

      m_polarisability[n] = m_solvent[solv].atom(j).polarisability;
      m_coscharge(n) = m_solvent[solv].atom(j).coscharge;
      m_damping_level(n) = m_solvent[solv].atom(j).damping_level;
      m_damping_power(n) = m_solvent[solv].atom(j).damping_power;
      m_is_polarisable[n] = bool(m_solvent[solv].atom(j).polarisability > 0.0);
      m_gamma[n] = m_solvent[solv].atom(j).gamma;
      m_gamma_j[n] = n + m_solvent[solv].atom(j).gamma_j;
      m_gamma_k[n] = n + m_solvent[solv].atom(j).gamma_k;
      m_cg_factor[n] = 1;
      // no exclusions or 1-4 interactions for solvent ?!
    }

    // add to the chargegroups
    DEBUG(8, "solvent cg: " << n);
    m_chargegroup.push_back(n);

    // and to the molecules
    m_molecule.push_back(n);
    // and temperature and pressure groups
    m_temperature_group.push_back(n);
    m_pressure_group.push_back(n);
  }

}

void topology::Topology::resolvate(unsigned int solv, unsigned int num_molecules) {
  if (solv != 0) {
    io::messages.add("resolvation for multiple solvents not implemented",
            "topology",
            io::message::error);
    return;
  }

  assert(m_num_solvent_atoms.size() == 1);
  assert(m_solvent.size() == 1);

  int n = num_solute_atoms();

  m_num_solvent_molecules[m_num_solvent_molecules.size() - 1] = num_molecules;
  m_num_solvent_atoms[m_num_solvent_atoms.size() - 1] = num_molecules * m_solvent[solv].num_atoms();

  DEBUG(5, "solvate: solvent atoms: " << num_solvent_atoms());
  DEBUG(10, "solvate: total atoms: " << num_solute_atoms() + num_solvent_atoms());

  DEBUG(8, "resizing for solute: " << num_solute_atoms() << " and solvent: "
          << num_solvent_atoms() << " (total: " << num_solute_atoms() + num_solvent_atoms());
  resize(num_solute_atoms() + num_solvent_atoms());

  m_chargegroup.erase(m_chargegroup.begin() + m_num_solute_chargegroups + 1, m_chargegroup.end());
  m_molecule.erase(m_molecule.begin() + m_num_solute_molecules + 1, m_molecule.end());
  m_temperature_group.erase(m_temperature_group.begin() + m_num_solute_temperature_groups + 1, m_temperature_group.end());
  m_pressure_group.erase(m_pressure_group.begin() + m_num_solute_pressure_groups + 1, m_pressure_group.end());

  DEBUG(9, "solute chargegroups: " << m_num_solute_chargegroups
          << " size: " << m_chargegroup.size());
  DEBUG(9, "solute molecules: " << m_num_solute_molecules
          << " size: " << m_molecule.size());
  DEBUG(9, "solute temperature groups: " << m_num_solute_temperature_groups
          << " size: " << m_temperature_group.size());
  DEBUG(9, "solute pressure groups: " << m_num_solute_pressure_groups
          << " size: " << m_pressure_group.size());

  // add to iac, mass, charge
  for (unsigned int i = 0; i < num_molecules; ++i) {
    for (unsigned int j = 0; j < m_solvent[solv].num_atoms(); ++j, ++n) {

      DEBUG(15, "iac[" << n << "]=" << m_solvent[solv].atom(j).iac);
      DEBUG(15, "charge[" << n << "]=" << m_solvent[solv].atom(j).charge);

      m_iac[n] = m_solvent[solv].atom(j).iac;
      m_mass(n) = m_solvent[solv].atom(j).mass;
      m_inverse_mass(n) = 1.0 / m_solvent[solv].atom(j).mass;
      m_charge(n) = m_solvent[solv].atom(j).charge;

      m_polarisability(n) = m_solvent[solv].atom(j).polarisability;
      m_coscharge(n) = m_solvent[solv].atom(j).coscharge;
      m_damping_level(n) = m_solvent[solv].atom(j).damping_level;
      m_damping_power(n) = m_solvent[solv].atom(j).damping_power;
      m_is_polarisable[n] = bool(m_solvent[solv].atom(j).polarisability > 0.0);
      m_gamma(n) = m_solvent[solv].atom(j).gamma;
      m_gamma_j[n] = m_solvent[solv].atom(j).gamma_j;
      m_gamma_k[n] = m_solvent[solv].atom(j).gamma_k;
      m_cg_factor[n] = 1;

      // no exclusions or 1-4 interactions for solvent
    }

    // add to the chargegroups
    DEBUG(8, "solvent cg: " << n);
    m_chargegroup.push_back(n);

    // and to the molecules
    DEBUG(8, "solvent mol: " << n);
    m_molecule.push_back(n);
    // and temperature and pressure groups
    m_temperature_group.push_back(n);
    m_pressure_group.push_back(n);
  }
}

/**
 * total number of solvent atoms.
 */
unsigned int topology::Topology::num_solvent_atoms()const {
  unsigned int n = 0;
  for (std::vector<unsigned int>::const_iterator it = m_num_solvent_atoms.begin(),
          to = m_num_solvent_atoms.end();
          it != to; ++it)
    n += *it;
  return n;
}

/**
 * calculate constraint degrees of freedom.
 */
void
topology::Topology::
calculate_constraint_dof(simulation::Multibath &multibath,
        bool rottrans_constraints, bool position_constraints,
        bool dih_constraints,
        bool ang_constraints)const {

  // save the positionally constrained atom inidices for fast checking
  // only make sure the distance constraints dof are not removed for
  // positionally constrained atoms.
  std::set<unsigned int> pos_cons_atom;
  if (position_constraints) {
    std::vector<topology::position_restraint_struct>::const_iterator
    it = position_restraints().begin(),
            to = position_restraints().end();

    for (; it != to; ++it) {
      pos_cons_atom.insert(it->seq);
    }
  }
  DEBUG(6, "Number of pos. constraint atoms: " << pos_cons_atom.size());
  DEBUG(6, "Number of non pos. constraint atoms: " << num_atoms() - pos_cons_atom.size());

  // substract constraints
  {
    std::vector<two_body_term_struct>::const_iterator
    c_it = solute().distance_constraints().begin(),
            c_to = solute().distance_constraints().end();

    unsigned int com_bath_i = 0, ir_bath_i = 0, com_bath_j = 0, ir_bath_j = 0;

    for (; c_it != c_to; ++c_it) {

      DEBUG(10, "Constraint: " << c_it->i << " - " << c_it->j);
      if (pos_cons_atom.find(c_it->i) != pos_cons_atom.end() &&
              pos_cons_atom.find(c_it->j) != pos_cons_atom.end()) {
        continue;
      } else if (pos_cons_atom.find(c_it->i) != pos_cons_atom.end()) {
        multibath.in_bath(c_it->j, com_bath_j, ir_bath_j);
        multibath[ir_bath_j].dof -= 1.0;
        multibath[ir_bath_j].ir_dof -= 1.0;
        multibath[ir_bath_j].solute_constr_dof += 1.0;
        continue;
      } else if (pos_cons_atom.find(c_it->j) != pos_cons_atom.end()) {
        multibath.in_bath(c_it->i, com_bath_i, ir_bath_i);
        multibath[ir_bath_i].dof -= 1.0;
        multibath[ir_bath_i].ir_dof -= 1.0;
        multibath[ir_bath_i].solute_constr_dof += 1.0;
        continue;
      }

      multibath.in_bath(c_it->i, com_bath_i, ir_bath_i);
      multibath.in_bath(c_it->j, com_bath_j, ir_bath_j);

      multibath[ir_bath_i].dof -= 0.5;
      multibath[ir_bath_j].dof -= 0.5;

      multibath[ir_bath_i].ir_dof -= 0.5;
      multibath[ir_bath_j].ir_dof -= 0.5;

      multibath[ir_bath_i].solute_constr_dof += 0.5;
      multibath[ir_bath_j].solute_constr_dof += 0.5;

    }

    for (unsigned int i = 0; i < multibath.size(); ++i) {
      DEBUG(7, "dof           " << multibath[i].dof);
      DEBUG(7, "solute constr " << multibath[i].solute_constr_dof);
    }

    // solvent constraints
    int index = num_solute_atoms();
    for (unsigned int s = 0; s < num_solvents(); ++s) {

      for (unsigned int m = 0; m < num_solvent_molecules(s); ++m) {

        c_it = solvent(s).distance_constraints().begin();
        c_to = solvent(s).distance_constraints().end();

        for (; c_it != c_to; ++c_it) {
          unsigned int c_i = index + c_it->i;
          unsigned int c_j = index + c_it->j;
          if (pos_cons_atom.find(c_i) != pos_cons_atom.end() &&
                  pos_cons_atom.find(c_j) != pos_cons_atom.end()) {
            continue;
          } else
            if (pos_cons_atom.find(c_i) != pos_cons_atom.end()) {
            multibath.in_bath(c_j, com_bath_j, ir_bath_j);
            multibath[ir_bath_j].dof -= 1.0;
            multibath[ir_bath_j].ir_dof -= 1.0;
            multibath[ir_bath_j].solvent_constr_dof += 1.0;
            continue;
          } else if (pos_cons_atom.find(c_j) != pos_cons_atom.end()) {
            multibath.in_bath(c_i, com_bath_i, ir_bath_i);
            multibath[ir_bath_i].dof -= 1.0;
            multibath[ir_bath_i].ir_dof -= 1.0;
            multibath[ir_bath_i].solvent_constr_dof += 1.0;
            continue;
          }

          multibath.in_bath(c_i, com_bath_i, ir_bath_i);
          multibath.in_bath(c_j, com_bath_j, ir_bath_j);

          multibath[ir_bath_i].dof -= 0.5;
          multibath[ir_bath_j].dof -= 0.5;

          multibath[ir_bath_i].ir_dof -= 0.5;
          multibath[ir_bath_j].ir_dof -= 0.5;

          multibath[ir_bath_i].solvent_constr_dof += 0.5;
          multibath[ir_bath_j].solvent_constr_dof += 0.5;

        }

        index += solvent(s).num_atoms();

      }
    }
  }

  DEBUG(7, "and the perturbed distance constraints (DOF calc)");

  {
    // substract perturbed constraints
    std::vector<perturbed_two_body_term_struct>
            ::const_iterator
            c_it = perturbed_solute().distance_constraints().begin(),
            c_to = perturbed_solute().distance_constraints().end();

    unsigned int com_bath_i = 0, ir_bath_i = 0, com_bath_j = 0, ir_bath_j = 0;

    for (; c_it != c_to; ++c_it) {

      DEBUG(10, "Constraint: " << c_it->i << " - " << c_it->j);
      if (pos_cons_atom.find(c_it->i) != pos_cons_atom.end() &&
              pos_cons_atom.find(c_it->j) != pos_cons_atom.end()) {
        // both atoms are position constrained
        continue;
      } else if (pos_cons_atom.find(c_it->i) != pos_cons_atom.end()) {
        // atom i position constrained
        multibath.in_bath(c_it->j, com_bath_j, ir_bath_j);
        multibath[ir_bath_j].dof -= 1.0;
        multibath[ir_bath_j].ir_dof -= 1.0;
        multibath[ir_bath_j].solute_constr_dof += 1.0;
        continue;
      } else if (pos_cons_atom.find(c_it->j) != pos_cons_atom.end()) {
        // atom j position constrained
        multibath.in_bath(c_it->i, com_bath_i, ir_bath_i);
        multibath[ir_bath_i].dof -= 1.0;
        multibath[ir_bath_i].ir_dof -= 1.0;
        multibath[ir_bath_i].solute_constr_dof += 1.0;
        continue;
      }

      // both atoms are not position-constrained
      multibath.in_bath(c_it->i, com_bath_i, ir_bath_i);
      multibath.in_bath(c_it->j, com_bath_j, ir_bath_j);

      multibath[ir_bath_i].dof -= 0.5;
      multibath[ir_bath_j].dof -= 0.5;

      multibath[ir_bath_i].ir_dof -= 0.5;
      multibath[ir_bath_j].ir_dof -= 0.5;

      multibath[ir_bath_i].solute_constr_dof += 0.5;
      multibath[ir_bath_j].solute_constr_dof += 0.5;

    }

    for (unsigned int i = 0; i < multibath.size(); ++i) {
      DEBUG(7, "dof           " << multibath[i].dof);
      DEBUG(7, "solute constr " << multibath[i].solute_constr_dof);
    }
  }

  DEBUG(7, "and the angle constraints (DOF calc)");
  if (ang_constraints) {
    std::vector<angle_restraint_struct>::const_iterator angit = angle_restraints().begin(), 
                                                   angto = angle_restraints().end();
    for (; angit != angto; ++angit){
      std::vector<int> not_pos_constrained;
      if (pos_cons_atom.find(angit->i) == pos_cons_atom.end()) not_pos_constrained.push_back(angit->i);
      if (pos_cons_atom.find(angit->j) == pos_cons_atom.end()) not_pos_constrained.push_back(angit->j);
      if (pos_cons_atom.find(angit->k) == pos_cons_atom.end()) not_pos_constrained.push_back(angit->k);

      double num_not_pos_const = not_pos_constrained.size();
      unsigned int ir_bath = 0, com_bath = 0;

      for (unsigned int i=0; i < num_not_pos_const; i++) {
        double part=1/num_not_pos_const;
        multibath.in_bath(not_pos_constrained[i], com_bath, ir_bath);
        multibath[ir_bath].dof-=part;
        multibath[ir_bath].ir_dof-=part;
        multibath[ir_bath].solute_constr_dof+=part;
      }
    }

    std::vector<perturbed_angle_restraint_struct>::const_iterator pangit = perturbed_angle_restraints().begin(), 
                                                   pangto = perturbed_angle_restraints().end();
    for (; pangit != pangto; ++pangit){
      std::vector<int> not_pos_constrained;
      if (pos_cons_atom.find(pangit->i) == pos_cons_atom.end()) not_pos_constrained.push_back(pangit->i);
      if (pos_cons_atom.find(pangit->j) == pos_cons_atom.end()) not_pos_constrained.push_back(pangit->j);
      if (pos_cons_atom.find(pangit->k) == pos_cons_atom.end()) not_pos_constrained.push_back(pangit->k);

      double num_not_pos_const = not_pos_constrained.size();
      unsigned int ir_bath = 0, com_bath = 0;

      for (unsigned int i=0; i < num_not_pos_const; i++) {
        double part=1/num_not_pos_const;
        multibath.in_bath(not_pos_constrained[i], com_bath, ir_bath);
        multibath[ir_bath].dof-=part;
        multibath[ir_bath].ir_dof-=part;
        multibath[ir_bath].solute_constr_dof+=part;
      }
    }

    for (unsigned int i = 0; i < multibath.size(); ++i) {
      DEBUG(7, "dof           " << multibath[i].dof);
      DEBUG(7, "solute constr " << multibath[i].solute_constr_dof);
    }
  }

  DEBUG(7, "and the dihedral constraints (DOF calc)");
  if (dih_constraints) {
    std::vector<dihedral_restraint_struct>::const_iterator dihit = dihedral_restraints().begin(), 
                                                   dihto = dihedral_restraints().end();
    for (; dihit != dihto; ++dihit){
      std::vector<int> not_pos_constrained;
      if (pos_cons_atom.find(dihit->i) == pos_cons_atom.end()) not_pos_constrained.push_back(dihit->i);
      if (pos_cons_atom.find(dihit->j) == pos_cons_atom.end()) not_pos_constrained.push_back(dihit->j);
      if (pos_cons_atom.find(dihit->k) == pos_cons_atom.end()) not_pos_constrained.push_back(dihit->k);
      if (pos_cons_atom.find(dihit->l) == pos_cons_atom.end()) not_pos_constrained.push_back(dihit->l);

      double num_not_pos_const = not_pos_constrained.size();
      unsigned int ir_bath = 0, com_bath = 0;

      for (unsigned int i=0; i < num_not_pos_const; i++) {
        double part=1/num_not_pos_const;
        multibath.in_bath(not_pos_constrained[i], com_bath, ir_bath);
        multibath[ir_bath].dof-=part;
        multibath[ir_bath].ir_dof-=part;
        multibath[ir_bath].solute_constr_dof+=part;
      }
    }

    std::vector<perturbed_dihedral_restraint_struct>::const_iterator pdihit = perturbed_dihedral_restraints().begin(), 
                                                   pdihto = perturbed_dihedral_restraints().end();
    for (; pdihit != pdihto; ++pdihit){
      std::vector<int> not_pos_constrained;
      if (pos_cons_atom.find(pdihit->i) == pos_cons_atom.end()) not_pos_constrained.push_back(pdihit->i);
      if (pos_cons_atom.find(pdihit->j) == pos_cons_atom.end()) not_pos_constrained.push_back(pdihit->j);
      if (pos_cons_atom.find(pdihit->k) == pos_cons_atom.end()) not_pos_constrained.push_back(pdihit->k);
      if (pos_cons_atom.find(pdihit->l) == pos_cons_atom.end()) not_pos_constrained.push_back(pdihit->l);

      double num_not_pos_const = not_pos_constrained.size();
      unsigned int ir_bath = 0, com_bath = 0;

      for (unsigned int i=0; i < num_not_pos_const; i++) {
        double part=1/num_not_pos_const;
        multibath.in_bath(not_pos_constrained[i], com_bath, ir_bath);
        multibath[ir_bath].dof-=part;
        multibath[ir_bath].ir_dof-=part;
        multibath[ir_bath].solute_constr_dof+=part;
      }
    }

    for (unsigned int i = 0; i < multibath.size(); ++i) {
      DEBUG(7, "dof           " << multibath[i].dof);
      DEBUG(7, "solute constr " << multibath[i].solute_constr_dof);
    }
  }

  if (position_constraints) {
    DEBUG(6, "position contraints dof");

    topology::Temperaturegroup_Iterator tmpit = temperature_group_begin(),
            tmpto = temperature_group_end();
    for (; tmpit != tmpto; ++tmpit) {
      unsigned int num_cons = 0;
      topology::Atom_Iterator ait = tmpit.begin(), ato = tmpit.end();
      unsigned int com_bath_i = 0, ir_bath_i = 0;
      unsigned int first_atom = *ait;
      multibath.in_bath(*ait, com_bath_i, ir_bath_i);

      // loop over atoms in temperature group
      for (; ait != ato; ++ait) {
        if (pos_cons_atom.find(*ait) != pos_cons_atom.end()) {
          ++num_cons;
        }
      }

      if (num_cons >= 1) {
        multibath[com_bath_i].com_dof -= 3.0;
        multibath[com_bath_i].dof -= 3.0;

        if (first_atom < num_solute_atoms())
          multibath[com_bath_i].solute_constr_dof += 3.0;
        else
          multibath[com_bath_i].solvent_constr_dof += 3.0;

        if (num_cons > 1) {
          multibath[ir_bath_i].ir_dof -= (num_cons - 1) * 3.0;
          multibath[ir_bath_i].dof -= (num_cons - 1) * 3.0;

          if (first_atom < num_solute_atoms())
            multibath[ir_bath_i].solute_constr_dof += (num_cons - 1) * 3.0;
          else
            multibath[ir_bath_i].solvent_constr_dof += (num_cons - 1) * 3.0;
        }
      }
    }

    for (unsigned int i = 0; i < multibath.size(); ++i) {
      DEBUG(7, "dof           " << multibath[i].dof);
      DEBUG(7, "solute constr " << multibath[i].solute_constr_dof);
    }
  }

  if (rottrans_constraints) {
    // check whether all solute is in one bath
    if (num_solute_atoms()) {

      unsigned int ir_bath = 0, ir_bath_0 = 0, com_bath = 0, com_bath_0 = 0;
      bool ok = true;

      multibath.in_bath(0, com_bath_0, ir_bath_0);
// only loop over the roto-translationally-constrained solute part
//      for (unsigned int i = 1; i < num_solute_atoms(); ++i) {
       // std::cout << "@@ m_rottrans_last_atom " << m_rottrans_last_atom << "\n";
      for (unsigned int i = 1; i < m_rottrans_last_atom; ++i) {
        multibath.in_bath(i, com_bath, ir_bath);
       // std::cout << "@@ now i " << i << " com_bath " << com_bath << " ir_bath " << ir_bath << "\n";
        if (com_bath != com_bath_0 || ir_bath != ir_bath_0) {
//          io::messages.add("roto-translational constraints: all solute has to be coupled "
          io::messages.add("roto-translational constraints: roto-translationally-constrained solute part has to be coupled "
                  "to one ir and one com bath", "calc_dof",
                  io::message::error);
          ok = false;
          break;
        }
      }
      if (ok) {
        std::cout << "roto-translational constraints: removing 3 dof from com bath "
                << com_bath_0 + 1 << " and from ir bath " << ir_bath_0 + 1 << "\n";
        multibath[com_bath_0].dof -= 3.0;
        multibath[com_bath_0].com_dof -= 3.0;
        multibath[ir_bath_0].dof -= 3.0;
        multibath[ir_bath_0].ir_dof -= 3.0;
      }
    }
  }

  // check whether we have dof in every (coupled) bath
  for (unsigned int i = 0; i < multibath.size(); ++i) {
    if (multibath[i].dof <= 0 && multibath[i].tau != -1) {
      io::messages.add("removing coupling of bath with 0 dof",
              "calc_dof",
              io::message::notice);
      multibath[i].tau = -1;
    }
  }

  DEBUG(10, "end dof calc");

}

void
topology::Topology::update_for_lambda() {
  DEBUG(10, "update for lambda");

  // update the individual lambdas
  const unsigned int size = m_individual_lambda_parameters.a.size();
  double lam = lambda();
  assert(size == m_individual_lambda_parameters.b.size() &&
          size == m_individual_lambda_parameters.c.size() &&
          size == m_individual_lambda_parameters.d.size() &&
          size == m_individual_lambda_parameters.e.size());

  assert(m_individual_lambda.size() == size &&
          m_individual_lambda_derivative.size() == size);

  const unsigned int num_energy_groups = m_energy_group.size();
  for (unsigned int i = 0; i < size; ++i) {
    assert(m_individual_lambda[i].size() == num_energy_groups);
    assert(m_individual_lambda_derivative[i].size() == num_energy_groups);
    for (unsigned int n1 = 0; n1 < num_energy_groups; ++n1) {
      assert(m_individual_lambda[i][n1].size() == num_energy_groups);
      assert(m_individual_lambda_derivative[i][n1].size() == num_energy_groups);
      for (unsigned int n2 = 0; n2 < num_energy_groups; ++n2) {
        std::pair<int, int> p(n1, n2);
        m_individual_lambda[i][n1][n2] =
                m_individual_lambda_parameters.a[i][p] * lam * lam * lam * lam +
                m_individual_lambda_parameters.b[i][p] * lam * lam * lam +
                m_individual_lambda_parameters.c[i][p] * lam * lam +
                m_individual_lambda_parameters.d[i][p] * lam +
                m_individual_lambda_parameters.e[i][p];
        m_individual_lambda_derivative[i][n1][n2] =
                4.0 * m_individual_lambda_parameters.a[i][p] * lam * lam * lam +
                3.0 * m_individual_lambda_parameters.b[i][p] * lam * lam +
                2.0 * m_individual_lambda_parameters.c[i][p] * lam +
                m_individual_lambda_parameters.d[i][p];
        m_individual_lambda[i][n2][n1] = m_individual_lambda[i][n1][n2];
        m_individual_lambda_derivative[i][n2][n1] =
                m_individual_lambda_derivative[i][n1][n2];

      }
    }
  }
  // update the masses using the correct lambda
  for (std::map<unsigned int, topology::Perturbed_Atom>::const_iterator
    it = perturbed_solute().atoms().begin(),
          to = perturbed_solute().atoms().end();
          it != to; ++it) {
    const int n1 = atom_energy_group()[it->second.sequence_number()];
    const double lambda = m_individual_lambda[simulation::mass_lambda][n1][n1];

    mass()(it->second.sequence_number()) =
            (1 - lambda) * it->second.A_mass() +
            lambda * it->second.B_mass();
    inverse_mass()(it->second.sequence_number()) = 1.0 / mass()(it->second.sequence_number());

    DEBUG(8, "mass A : " << it->second.A_mass() << " B : "
            << it->second.B_mass());
    DEBUG(8, "mass(" << it->second.sequence_number()
            << ") = " << mass()(it->second.sequence_number()));
  }
}

void topology::Topology::update_chargegroup_exclusion() {
  m_chargegroup_exclusion.clear();
  m_chargegroup_exclusion.resize(num_solute_chargegroups());

  for (size_t cg1 = 0; cg1 < num_solute_chargegroups(); ++cg1) {

    m_chargegroup_exclusion[cg1].insert(cg1);

    for (size_t cg2 = cg1 + 1; cg2 < num_solute_chargegroups(); ++cg2) {

      // std::cerr << "\tchecking cg1=" << cg1 << " cg2=" << cg2 << std::endl;

      for (int at1 = m_chargegroup[cg1]; at1 < m_chargegroup[cg1 + 1]; ++at1) {

              topology::excl_cont_t::value_type::const_iterator ex = m_all_exclusion[at1].begin(),
                ex_to = m_all_exclusion[at1].end();
        for (; ex != ex_to; ++ex) {

          // std::cerr << "cg2: " << m_chargegroup[cg2] << " - " << m_chargegroup[cg2+1]
          // << " ex=" << *ex << std::endl;

          if (m_chargegroup[cg2] <= *ex &&
                  m_chargegroup[cg2 + 1] > *ex) {

            // std::cerr << "cg1=" << cg1 << " cg2=" << cg2 << " : excluded!" << std::endl;

            m_chargegroup_exclusion[cg1].insert(cg2);
            m_chargegroup_exclusion[cg2].insert(cg1);

          } // exclude cg
        } // exclusions of at1
      } // at1 in cg1
    } // cg2
  } // cg1
}

void topology::Topology::update_all_exclusion() {
  const unsigned size = m_all_exclusion.size();
  const unsigned n_atoms = num_atoms();
  assert(size == n_atoms
      && size == m_exclusion.size()
      && size == m_one_four_pair.size());
  m_all_exclusion.clear();
  m_all_exclusion.resize(n_atoms);

  DEBUG(15, "Merging exclusions and one_four_pairs into all_exclusions");
  // Insert exclusions and 1,4-pairs
  for (unsigned i = 0; i < n_atoms; ++i)
  {
    std::set_union(m_exclusion[i].begin(), m_exclusion[i].end(),
          m_one_four_pair[i].begin(), m_one_four_pair[i].end(),
          std::inserter(m_all_exclusion[i], m_all_exclusion[i].end()));
  }
  // Insert LJ exceptions
  DEBUG(15, "Inserting LJ exceptions");
  for (std::vector<topology::lj_exception_struct>::const_iterator
        it = m_lj_exceptions.begin(); it != m_lj_exceptions.end(); ++it) {
    DEBUG(15, "LJ exception: " << it->i << " - " << it->j);
    assert((unsigned) it->i < size && it->i < it->j);
    m_all_exclusion[it->i].insert(it->j);
  }
  /// Chargegroup exclusions need to be updated too
  update_chargegroup_exclusion();
}

/**
 * check state
 */
int
topology::Topology::check_state()const {
  int result = 0;

  // check that we have an energy group for every atom
  if (m_atom_energy_group.size() != num_atoms()) {
    io::messages.add("not every atom has an energy group index",
            "Topology::check_state", io::message::error);
    ++result;
  }
  for (std::vector<unsigned int>::const_iterator it = m_atom_energy_group.begin(),
          to = m_atom_energy_group.end(); it != to; ++it) {
    if (*it >= m_energy_group.size()) {
      io::messages.add("energy group index of atom too large",
              "Topology::check_state", io::message::error);
      ++result;
    }
  }

  return result;
}

double topology::Topology::squared_sum_charges() const {
  double sum = 0.0;
  const unsigned int natoms = num_atoms();
  for (unsigned int i = 0; i < natoms; i++)
    sum += charge(i);
  return sum*sum;
}

double topology::Topology::sum_squared_charges() const {
  double sum = 0.0;
  const unsigned int natoms = num_atoms();
  for (unsigned int i = 0; i < natoms; i++)
    sum += charge(i) * charge(i);
  return sum;
}

namespace topology {

  /**
   * output information about the topology.
   */
  std::ostream & operator<<(std::ostream &os, Topology &topo) {
    os << "a topology";
    return os;
  }
}

