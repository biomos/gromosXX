/**
 * @file multibath.cc
 * the multibath parameter class.
 */

#include "../stdheader.h"

#include "../simulation/multibath.h"
#include "../simulation/simulation.h"

#include "../configuration/energy.h"
#include "../topology/topology.h"

#include "../util/error.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE simulation
#define SUBMODULE simulation

void simulation::Multibath
::calculate_degrees_of_freedom(topology::Topology &topo,
        bool rottrans_constraints,
        bool position_constraints,
        double dof_to_subtract,
        bool dih_constraints,
        bool ang_constraints) {
  // check whether we have at least one bath
  if (size() == 0) {
    io::messages.add("Adding a bath, no temperature coupling",
            "Multibath::calculate_degrees_of_freedom",
            io::message::notice);
    add_bath(0.0);
  }

  // check whether the last last is really the last_atom
  if ((!(m_bath_index.size() == 0))
          && (m_bath_index.end() - 1)->last_atom != topo.num_atoms() - 1) {
    io::messages.add("Last atom of last bath is not the last atom in the sytem!",
            "Multibath::calculate_degrees_of_freedom",
            io::message::error);
    return;
    /*
    io::messages.add("Adding atoms to the last bath!",
         "Multibath::calculate_degrees_of_freedom",
         io::message::notice);
    add_bath_index(topo.num_atoms()-1, unsigned(topo.temperature_groups().size())-1, 
       unsigned(size())-1, unsigned(size())-1);
     */
  }
  if (m_bath_index.size() == 0) {
    io::messages.add("Adding atoms to the last bath!",
            "Multibath::calculate_degrees_of_freedom",
            io::message::notice);
    add_bath_index(topo.num_atoms() - 1, unsigned(topo.temperature_groups().size()) - 1,
            unsigned(size()) - 1, unsigned(size()) - 1);
  }

  // loop over the ranges
  std::vector<bath_index_struct>::iterator it = m_bath_index.begin(),
          to = m_bath_index.end();

  DEBUG(8, "number of baths: " << unsigned(size()));
  DEBUG(8, "temperature groups: " << unsigned(topo.temperature_groups().size()));

  for (int last = -1; it != to; ++it) {
    // get the number of temperature groups in the range
    int num_tg = 0;
    int tg = 0;

    DEBUG(8, "last atom: " << it->last_atom);
    DEBUG(8, "end of last group: " << last);

    for (topology::Temperaturegroup_Iterator tg_it = topo.temperature_group_begin(),
            tg_to = topo.temperature_group_end();
            tg_it != tg_to;
            ++tg_it, ++tg) {

      DEBUG(8, "current temperature group begins: " << (*tg_it.begin()));

      if ((*(tg_it.begin())) > it->last_atom) {
        break;
      }

      if (int(*(tg_it.begin())) > last)
        ++num_tg;
    }

    // add the last temperature group
    DEBUG(8, "last temperature group is " << tg - 1);
    it->last_temperature_group = tg - 1;

    DEBUG(8, "num temperature groups is " << num_tg);
    // add the molecular translational dof
    (*this)[it->com_bath].dof += num_tg * 3;
    (*this)[it->com_bath].com_dof += num_tg * 3;

    // and the internal and molecular rotational dof
    (*this)[it->ir_bath].dof += (it->last_atom - last) * 3 - num_tg * 3;
    (*this)[it->ir_bath].ir_dof += (it->last_atom - last) * 3 - num_tg * 3;

    last = it->last_atom;
  }

  // subtract constraints
  topo.calculate_constraint_dof(*this, rottrans_constraints, position_constraints, dih_constraints, ang_constraints);

  // subtract user dof (NDFMIN in BOUNDCOND block)
  //
  // Martin Stroet edit:
  //     loop over baths and remove bath_dof/total_dof from every bath
  int dof_total = 0;
  for (simulation::Multibath::iterator it = begin(), to = end();
            it != to; ++it) {
    dof_total += it->dof;
  }
  for (simulation::Multibath::iterator it = begin(), to = end();
       it != to; ++it) {
    it->dof -= dof_to_subtract*(it->dof/dof_total);
  }

  //dof_to_subtract /= double(size());
  //for (simulation::Multibath::iterator it = begin(), to = end();
  //        it != to; ++it) {
  //  it->dof -= dof_to_subtract;
  //}

  // check whether all temperature groups are attached to only one bath
  for (topology::Temperaturegroup_Iterator tg_it = topo.temperature_group_begin(),
          tg_to = topo.temperature_group_end(); tg_it != tg_to; ++tg_it) {
    unsigned int first_com = 0, first_ir = 0;
    topology::Atom_Iterator a_it = tg_it.begin(),
            a_to = tg_it.end();
    in_bath(*a_it, first_com, first_ir);
    ++a_it;
    for (; a_it != a_to; ++a_it) {
      unsigned int com_i = 0, ir_i = 0;
      in_bath(*a_it, com_i, ir_i);
      if (com_i != first_com || ir_i != first_ir) {
        io::messages.add("Multibath: Temperature group distributed over multiple baths.",
                "Multibath::check_state",
                io::message::error);
        return;
      }
    }
  }


}

int simulation::Multibath::check_state(unsigned int num_atoms)const {
  DEBUG(8, "checking multibath state");
  int result = 0;
  unsigned int last_atom = 0;
  std::vector<bath_index_struct>::const_iterator it = m_bath_index.begin(),
          to = m_bath_index.end();
  for (; it != to; ++it) {
    if (it->last_atom < last_atom) {
      io::messages.add("Multibath not sorted", "Multibath::check_state",
              io::message::error);
      ++result;
    }
    if (it->last_atom > num_atoms) {
      io::messages.add("Multibath last atom index too large",
              "Multibath::check_state",
              io::message::error);
      ++result;
    }
    if (it->com_bath >= size()) {
      io::messages.add("Multibath: com bath index out of range",
              "Multibath::check_state",
              io::message::error);
      return E_INPUT_ERROR;
    }
    if (it->ir_bath >= size()) {
      io::messages.add("Multibath: ir bath index out of range",
              "Multibath::check_state",
              io::message::error);
      return E_INPUT_ERROR;
    }
    if ((*this)[it->com_bath].dof == 0)
      io::messages.add("Multibath: bath with 0 degrees of freedom?",
            "Multibath::check_state",
            io::message::warning);
    if ((*this)[it->ir_bath].dof == 0)
      io::messages.add("Multibath: bath with 0 degrees of freedom?",
            "Multibath::check_state",
            io::message::warning);
    if ((*this)[it->ir_bath].solute_constr_dof < 0
            || (*this)[it->ir_bath].solvent_constr_dof < 0) {
      io::messages.add("Multibath: constrained degrees of freedom negative",
              "Multibath::check_state",
              io::message::error);
      ++result;
    }
    if ((*this)[it->com_bath].solute_constr_dof < 0
            || (*this)[it->com_bath].solvent_constr_dof < 0) {
      io::messages.add("Multibath: constrained degrees of freedom negative",
              "Multibath::check_state",
              io::message::error);
      ++result;
    }
    if ((*this)[it->ir_bath].tau < 0 && (*this)[it->ir_bath].tau != -1) {
      io::messages.add("Multibath: tau < 0 && tau != -1",
              "Multibath::check_state",
              io::message::error);
      ++result;
    }
    if ((*this)[it->com_bath].tau < 0 && (*this)[it->com_bath].tau != -1) {
      io::messages.add("Multibath: tau < 0 && tau != -1",
              "Multibath::check_state",
              io::message::error);
      ++result;
    }
    if ((*this)[it->ir_bath].temperature < 0 ||
            (*this)[it->com_bath].temperature < 0) {
      io::messages.add("Multibath: temperature < 0",
              "Multibath::check_state",
              io::message::error);
      ++result;
    }

  }

  return result;

}

void simulation::Multibath::calc_totals(configuration::Energy const &energy,
        double & ekin, double & ekin_mol,
        double & ekin_ir,
        double & temp, double & temp_mol,
        double & temp_ir,
        double & scale)const {
  ekin = 0.0;
  ekin_mol = 0.0;
  ekin_ir = 0.0;
  scale = 0.0;

  double sum_dof = 0.0;
  double sum_dof_com = 0.0;
  double sum_dof_ir = 0.0;
  double tau_dof = 0.0;

  std::vector<bath_struct>::const_iterator
  it = begin(),
          to = end();

  for (unsigned int i = 0; it != to; ++it, ++i) {

    const double e_kin = energy.kinetic_energy[i];
    const double e_kin_com = energy.com_kinetic_energy[i];
    const double e_kin_ir = energy.ir_kinetic_energy[i];

    if (it->tau != -1) {
      tau_dof += it->dof;
      scale += it->scale * it->dof;
    }

    sum_dof += it->dof;
    sum_dof_ir += it->ir_dof;
    sum_dof_com += it->com_dof;

    ekin += e_kin;
    ekin_mol += e_kin_com;
    ekin_ir += e_kin_ir;
  }

  if (sum_dof)
    temp = 2 * ekin / (math::k_Boltzmann * sum_dof);
  else temp = 0.0;

  if (sum_dof_com)
    temp_mol = 2 * ekin_mol / (math::k_Boltzmann * sum_dof_com);
  else temp_mol = 0.0;

  if (sum_dof_ir)
    temp_ir = 2 * ekin_ir / (math::k_Boltzmann * sum_dof_ir);
  else temp_ir = 0.0;

  if (tau_dof)
    scale /= tau_dof;
  else scale = 0.0;

}
