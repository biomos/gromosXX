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
 * @file create_bonded.cc
 * create the bonded terms.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../interaction/forcefield/forcefield.h"

// interactions
#include "../../interaction/interaction_types.h"
#include "../../interaction/bonded/quartic_bond_interaction.h"
#include "../../interaction/bonded/harmonic_bond_interaction.h"
#include "../../interaction/bonded/cg_bond_interaction.h"
#include "../../interaction/bonded/angle_interaction.h"
#include "../../interaction/bonded/harm_angle_interaction.h"
#include "../../interaction/bonded/dihedral_interaction.h"
#include "../../interaction/bonded/dihedral_new_interaction.h"
#include "../../interaction/bonded/crossdihedral_interaction.h"
#include "../../interaction/bonded/improper_dihedral_interaction.h"

// perturbed interactions
#include "../../interaction/bonded/perturbed_quartic_bond_interaction.h"
#include "../../interaction/bonded/perturbed_harmonic_bond_interaction.h"
#include "../../interaction/bonded/perturbed_soft_bond_interaction.h"
#include "../../interaction/bonded/perturbed_cg_bond_interaction.h"
#include "../../interaction/bonded/perturbed_angle_interaction.h"
#include "../../interaction/bonded/perturbed_soft_angle_interaction.h"
#include "../../interaction/bonded/perturbed_improper_dihedral_interaction.h"
#include "../../interaction/bonded/perturbed_soft_improper_interaction.h"
#include "../../interaction/bonded/perturbed_dihedral_interaction.h"
#include "../../interaction/bonded/perturbed_dihedral_new_interaction.h"

// #include "../../io/instream.h"
#include "../../io/ifp.h"

#include "create_bonded.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

int interaction::create_g96_bonded(interaction::Forcefield & ff,
				   topology::Topology const & topo,
				   simulation::Simulation const & sim,
				   io::IFP & it,
				   std::ostream & os,
				   bool quiet)
{
  const simulation::Parameter & param = sim.param();
  DEBUG(8, "creating g96 bonded");
  if (param.force.bond == 1){
    if (!quiet)
      os << "\tquartic bond interaction\n";

    interaction::Quartic_Bond_Interaction *b =
              new interaction::Quartic_Bond_Interaction();

    ff.push_back(b);

    if (param.perturbation.perturbation) {
      if (!quiet)
        os << "\tperturbed quartic bond interaction\n";

      interaction::Perturbed_Quartic_Bond_Interaction * pb =
        new interaction::Perturbed_Quartic_Bond_Interaction(*b);
      ff.push_back(pb);
    }
  } else if (param.force.bond == 2) {
    if (!quiet)
      os << "\tharmonic bond interaction\n";

    interaction::Harmonic_Bond_Interaction *b =
      new interaction::Harmonic_Bond_Interaction();

    ff.push_back(b);

    io::messages.add("using harmonic bond potential",
                     "create bonded", io::message::notice);

    if (param.perturbation.perturbation) {
      if (!quiet)
        os << "\tperturbed harmonic bond interaction\n";

      interaction::Perturbed_Harmonic_Bond_Interaction * pb =
        new interaction::Perturbed_Harmonic_Bond_Interaction(*b);
      ff.push_back(pb);
    }
  } 
  
  // do the perturbed soft bonds regardless of constraints
  // the affected bonds will have been removed from bonds and constraints
  if (param.perturbation.perturbation) { 
    if (param.force.bond) {
      if (!quiet)
        os << "\tperturbed soft harmonic bond interaction\n";   
      interaction::Perturbed_Soft_Bond_Interaction * sb =
        new interaction::Perturbed_Soft_Bond_Interaction();
      ff.push_back(sb);
    }
    if (param.force.angle) {
      if (!quiet)
        os << "\tperturbed soft harmonic angle interaction\n";   
      interaction::Perturbed_Soft_Angle_Interaction * sa =
        new interaction::Perturbed_Soft_Angle_Interaction();
      ff.push_back(sa);
    }
    if (param.force.improper) {
      if (!quiet)
        os << "\tperturbed soft improper dihedral interaction\n";   
      interaction::Perturbed_Soft_Improper_Interaction * si =
        new interaction::Perturbed_Soft_Improper_Interaction();
      ff.push_back(si);
    }
  }

  if (param.cgrain.level > 1) {
    if (!quiet)
      os << "\tdipole-particle bond interaction\n";

    interaction::DP_Bond_Interaction * bcg =
      new interaction::DP_Bond_Interaction();

    ff.push_back(bcg);

    if (param.perturbation.perturbation) {
      if (!quiet)
        os << "\tperturbed dipole-particle bond interaction\n";

      interaction::Perturbed_DP_Bond_Interaction * pbcg =
              new interaction::Perturbed_DP_Bond_Interaction(*bcg);
      ff.push_back(pbcg);
    }

  }

  if (param.force.angle == 1){
    if (!quiet)
      os <<"\tbond angle (cosine) interaction\n";
    interaction::Angle_Interaction *a =
      new interaction::Angle_Interaction();

    ff.push_back(a);

    if (param.perturbation.perturbation){
      if (!quiet)
	os <<"\tperturbed bond angle interaction\n";
      interaction::Perturbed_Angle_Interaction * pa =
	new interaction::Perturbed_Angle_Interaction(*a);
      ff.push_back(pa);
    }
  }

  if (param.force.angle == 2){
    if (!quiet)
      os <<"\tharmonic bond angle interaction\n";
    interaction::Harm_Angle_Interaction *a =
      new interaction::Harm_Angle_Interaction();

    ff.push_back(a);

    if (param.perturbation.perturbation){
      io::messages.add("perturbed harmonic (g87) angle potential not implemented",
		       "create bonded", io::message::error);
      /*
      if (!quiet)
	os <<"\tperturbed bond angle interaction\n";
      interaction::Perturbed_Angle_Interaction * pa =
	new interaction::Perturbed_Angle_Interaction(*a);
      ff.push_back(pa);
      */
    }
  }

  if (param.force.improper == 1){
    if (!quiet)
      os << "\timproper dihedral interaction\n";

    interaction::Improper_Dihedral_Interaction * i =
      new interaction::Improper_Dihedral_Interaction();
    ff.push_back(i);

    if (param.perturbation.perturbation){
      if(!quiet)
	os << "\tperturbed improper dihedral interaction\n";
      interaction::Perturbed_Improper_Dihedral_Interaction * pi =
	new interaction::Perturbed_Improper_Dihedral_Interaction(*i);
      ff.push_back(pi);
    }

  }

  if (param.force.dihedral == 1){
    if (!quiet)
      os <<"\tdihedral interaction\n";

    interaction::Dihedral_new_Interaction * d =
      new interaction::Dihedral_new_Interaction();
    ff.push_back(d);

    if (param.perturbation.perturbation){
      if(!quiet)
	os <<"\tperurbed dihedral interaction\n";
      interaction::Perturbed_Dihedral_new_Interaction * pd =
	new interaction::Perturbed_Dihedral_new_Interaction(*d);
      ff.push_back(pd);
    }
  }
   if (param.force.dihedral == 2){
    if (!quiet)
      os <<"\tdihedral interaction\n";

    interaction::Dihedral_Interaction * d =
      new interaction::Dihedral_Interaction();
    ff.push_back(d);

    if (param.perturbation.perturbation){
      if(!quiet)
	os <<"\tperurbed dihedral interaction\n";
      interaction::Perturbed_Dihedral_Interaction * pd =
	new interaction::Perturbed_Dihedral_Interaction(*d);
      ff.push_back(pd);
    }
  }
  if (param.force.crossdihedral == 1){
    if (!quiet)
      os <<"\tcrossdihedral interaction\n";

    interaction::Crossdihedral_Interaction * c =
      new interaction::Crossdihedral_Interaction();
    ff.push_back(c);

    /*if (param.perturbation.perturbation){
      if(!quiet)
	os <<"\tperurbed crossdihedral interaction\n";
      interaction::Perturbed_Crossdihedral_Interaction * pc =
	new interaction::Perturbed_Crossdihedral_Interaction(*c);
      ff.push_back(pc);
    }*/
  }

  return 0;

}
