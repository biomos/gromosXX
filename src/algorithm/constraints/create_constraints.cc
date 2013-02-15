/**
 * @file create_constraints.cc
 * create the constraint algorithm
 * choosing from SHAKE, LINCS
 * and enabling or disabling perturbation
 */

#include "../../stdheader.h"
#include <fstream>

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../configuration/state_properties.h"

#include "../../algorithm/algorithm/algorithm_sequence.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../interaction/forcefield/forcefield.h"

#include "../../math/periodicity.h"
#include "../../math/volume.h"

#include "../../io/argument.h"
#include "../../io/blockinput.h"
#include "../../io/instream.h"
#include "../../io/topology/in_topology.h"

#include "../../util/error.h"

#include "../../algorithm/constraints/shake.h"
#include "../../algorithm/constraints/perturbed_shake.h"
#include "../../algorithm/constraints/lincs.h"
#include "../../algorithm/constraints/flexible_constraint.h"
#include "../../algorithm/constraints/perturbed_flexible_constraint.h"
#include "../../algorithm/constraints/settle.h"
//#include "../../algorithm/constraints/gpu_settle.h"
#include "../../algorithm/constraints/m_shake.h"
#include "../../algorithm/constraints/gpu_shake.h"

#include "../../algorithm/constraints/rottrans.h"

#include "../../io/print_block.h"

#include "create_constraints.h"


#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

int algorithm::create_constraints(algorithm::Algorithm_Sequence &md_seq,
        topology::Topology &topo,
        simulation::Simulation & sim,
        io::In_Topology &it,
        bool quiet)
 {
  DEBUG(7, "solute:  " << sim.param().constraint.solute.algorithm);
  DEBUG(7, "solvent: " << sim.param().constraint.solvent.algorithm);
  DEBUG(7, "\tNTC: " << sim.param().constraint.ntc);
  DEBUG(7, "rottrans: " << sim.param().rottrans.rottrans);


  // CONSTRAINTS
  DEBUG(7, "Constrain solute?");

  // solute constraints have to be set to SHAKE, but NTC may be 1
  // for dihedral constraints to be used.
  if (sim.param().dihrest.dihrest == simulation::dihedral_constr &&
          sim.param().constraint.solute.algorithm != simulation::constr_shake) {
    io::messages.add("dihedral constraints require SHAKE as solute algorithm.",
            "create_constraints", io::message::error);
    return 1;
  }

  switch (sim.param().constraint.solute.algorithm) {
    case simulation::constr_shake :
    {
      if (!sim.param().perturbation.perturbation) {
        // SHAKE
        algorithm::Shake * s =
                new algorithm::Shake
                (sim.param().constraint.solute.shake_tolerance);
        it.read_harmonic_bonds(s->parameter());
        md_seq.push_back(s);

      } else {
        // perturbed shake also calls normal shake...
        algorithm::Perturbed_Shake * ps =
                new algorithm::Perturbed_Shake
                (sim.param().constraint.solute.shake_tolerance);
        it.read_harmonic_bonds(ps->parameter());
        md_seq.push_back(ps);

      }
      break;
    }
    case simulation::constr_lincs :
    {
      algorithm::Lincs * s =
              new algorithm::Lincs;
      it.read_harmonic_bonds(s->parameter());
      md_seq.push_back(s);

      if (sim.param().perturbation.perturbation) {
        io::messages.add("no free energy derivatives for LINCS, so you better don't "
                "change constrained bond lengths", "create_constraints",
                io::message::warning);
      }
      break;
    }
    case simulation::constr_flexshake :
    {
      // let's try to get the forcefield
      interaction::Forcefield * ff = NULL;

      for (size_t i = 0; i < md_seq.size(); ++i) {
        if (md_seq[i]->name == "Forcefield") {
          DEBUG(8, "flexible shake: forcefield found");
          ff = dynamic_cast<interaction::Forcefield *> (md_seq[i]);
          break;
        }
      }

      if (!ff) {
        io::messages.add("no forcefield found", "create_constraints",
                io::message::error);
      }

      if (!sim.param().perturbation.perturbation) {

        algorithm::Flexible_Constraint * fs =
                new algorithm::Flexible_Constraint
                (sim.param().constraint.solute.shake_tolerance, 1000, ff);

        it.read_harmonic_bonds(fs->parameter());

        md_seq.push_back(fs);

      } else {

        algorithm::Perturbed_Flexible_Constraint * pfc =
                new algorithm::Perturbed_Flexible_Constraint
                (sim.param().constraint.solute.shake_tolerance, 1000, ff);

        it.read_harmonic_bonds(pfc->parameter());

        md_seq.push_back(pfc);

      }

      break;

    }
    case simulation::constr_settle :
    {
      io::messages.add("SETTLE is only available for solvent.", "create_constraints",
              io::message::error);
    }
    case simulation::constr_m_shake :
    {
      io::messages.add("M_SHAKE is only available for solvent.", "create_constraints",
              io::message::error);
    }
    case simulation::constr_gpu_settle:
    {
      io::messages.add("SETTLE is only available for solvent.", "create_constraints",
              io::message::error);
    }
     case simulation::constr_gpu_shake :
    {
      io::messages.add("GPU_SHAKE is only available for solvent.", "create_constraints",
              io::message::error);
    }
    default:
    {
      // no constraints...
    }
  }

  // solvent (if not the same as solute)
  if (topo.num_solvent_atoms() > 0 &&
          sim.param().constraint.solute.algorithm !=
          sim.param().constraint.solvent.algorithm) {

    switch (sim.param().constraint.solvent.algorithm) {
      case simulation::constr_shake :
      {
        // SHAKE
        algorithm::Shake * s =
                new algorithm::Shake
                (sim.param().constraint.solvent.shake_tolerance);
        it.read_harmonic_bonds(s->parameter());
        md_seq.push_back(s);

        break;
      }
      case simulation::constr_lincs :
      {
        algorithm::Lincs * s =
                new algorithm::Lincs;
        it.read_harmonic_bonds(s->parameter());
        md_seq.push_back(s);

        break;
      }
      case simulation::constr_flexshake :
      {
        io::messages.add("Flexible Shake not implemented for solvent",
                "create_constraints", io::message::error);
        break;
      }
      case simulation::constr_settle :
      {
        algorithm::Settle * s =
                new algorithm::Settle;
        it.read_harmonic_bonds(s->parameter());
        md_seq.push_back(s);

        break;
      }
      case simulation::constr_m_shake :
      {
        // first check that the solvent(s) has only three constraints!
        for (unsigned int i = 0; i < topo.num_solvents(); ++i) {
          if (topo.solvent(i).distance_constraints().size() != 3)
            io::messages.add("M_Shake only implemented for solvents with three constraints",
                "create_constraints", io::message::error);
        }
        algorithm::M_Shake * s =
                new algorithm::M_Shake(sim.param().constraint.solvent.shake_tolerance);
        it.read_harmonic_bonds(s->parameter());
        md_seq.push_back(s);

        break;
      }
      case simulation::constr_gpu_settle :
      {
        /*algorithm::GPU_Settle * s = new algorithm::GPU_Settle;
        it.read_harmonic_bonds(s->parameter());
        md_seq.push_back(s); */

        break;
      }
      case simulation::constr_gpu_shake :
      {
        // first check that the solvent(s) has only three constraints!
        for (unsigned int i = 0; i < topo.num_solvents(); ++i) {
          if (topo.solvent(i).distance_constraints().size() != 3)
            io::messages.add("GPU_Shake only implemented for solvents with three constraints",
                "create_constraints", io::message::error);
        }
        algorithm::GPU_Shake * s =
                new algorithm::GPU_Shake(sim.param().constraint.solvent.shake_tolerance);
        it.read_harmonic_bonds(s->parameter());
        md_seq.push_back(s);

        break;
      }
      default:
      {
        // no constraints
        // should already be warned from In_Parameter
      }
    }
  }

  // roto-translational constraints
  if (sim.param().rottrans.rottrans) {
    DEBUG(8, "creating roto-translational constraints");

    algorithm::Rottrans_Constraints * rtc =
            new algorithm::Rottrans_Constraints();
    md_seq.push_back(rtc);
  }

  return 0;
}

