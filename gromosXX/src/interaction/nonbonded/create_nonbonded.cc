/**
 * @file create_nonbonded.cc
 * create the nonbonded interaction.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>

#include <interaction/interaction.h>
#include <interaction/interaction_types.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/forcefield/forcefield.h>

#include <configuration/energy.h>
#include <interaction/nonbonded/interaction/storage.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/standard_pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/grid_pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/vgrid_pairlist_algorithm.h>

#include <interaction/nonbonded/interaction/nonbonded_outerloop.h>
#include <interaction/nonbonded/interaction/nonbonded_set.h>

#include <interaction/nonbonded/interaction/nonbonded_interaction.h>
#include <interaction/nonbonded/interaction/omp_nonbonded_interaction.h>
#include <interaction/nonbonded/interaction/mpi_nonbonded_master.h>
#include <interaction/nonbonded/interaction/mpi_nonbonded_slave.h>

#include <io/ifp.h>

#include <interaction/nonbonded/create_nonbonded.h>

#include <util/debug.h>

#ifdef XXMPI
#include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

using namespace std;

int interaction::create_g96_nonbonded
(
 interaction::Forcefield & ff,
 topology::Topology const & topo,
 simulation::Simulation const & sim,
 io::IFP & it,
 std::ostream & os,
 bool quiet
 )
{
  if (sim.param().force.nonbonded_vdw == 0 &&
      sim.param().force.nonbonded_crf == 0) return 0;
  
  if(!quiet){
    
    os << "\t" << setw(20) << left << "nonbonded force";
    if (sim.param().force.nonbonded_vdw)
      os << setw(30) << left << "van-der-Waals";
    if (sim.param().force.nonbonded_crf)
      os << setw(30) << left << "Coulomb-reaction-field";    
    os << "\n";
    
    if (sim.param().pairlist.grid == 0)
      os << "\t" << setw(20) << left << "Pairlist Algorithm" << setw(30) 
	   << left << "Standard Pairlist Algorithm" << right << "\n";
    else if (sim.param().pairlist.grid == 1){
      os << "\t" << setw(20) << left << "Pairlist Algorithm" << setw(30) 
	   << left << "Grid Pairlist Algorithm" << right << "\n";
    }
    else if (sim.param().pairlist.grid == 2){
      os << "\t" << setw(20) << left << "Pairlist Algorithm" << setw(30) 
	   << left << "Vector Grid Pairlist Algorithm" << right << "\n";
    }
    
    switch(sim.param().boundary.boundary){
      case math::vacuum:
	os << "\t" << setw(20) << left << "boundary" << setw(30) 
	     << left << "vacuum" << right << "\n";
	break;
      case math::rectangular:
	os << "\t" << setw(20) << left << "boundary" << setw(30) 
	     << left << "rectangular" << right << "\n";
	break;
      case math::truncoct:
	os << "\t" << setw(20) << left << "boundary" << setw(30) 
	     << left << "truncoct" << right << "\n";
	break;
      case math::triclinic:
	os << "\t" << setw(20) << left << "boundary" << setw(30) 
	     << left << "triclinic" << right << "\n";
	break;
      default:
	os << "\t" << setw(20) << left << "boundary" << setw(30) 
	     << left << "unknown" << right << "\n";
    }

    switch(sim.param().pcouple.virial){
      case math::no_virial:
	os << "\t" << setw(20) << left << "virial" << setw(30) << left << "none" << right << "\n";
	break;
      case math::molecular_virial:
	os << "\t" << setw(20) << left << "virial" << setw(30) << left << "molecular" << right << "\n";
	break;
      case math::atomic_virial:
	os << "\t" << setw(20) << left << "virial" << setw(30) << left << "atomic" << right << "\n";
	break;
      default:
	os << "\t" << setw(20) << left << "virial" << setw(30) << left << "unknown" << right << "\n";
    }

    if (sim.param().pairlist.atomic_cutoff)
      os << "\t" << setw(20) << left << "cutoff" << setw(30) << left << "atomic" << right << "\n";
    else
      os << "\t" << setw(20) << left << "cutoff" << setw(30) << left << "chargegroup" << right << "\n";

  }
  
  Pairlist_Algorithm * pa;
  if (sim.param().pairlist.grid == 0){
    pa = new Standard_Pairlist_Algorithm();
  }
  else if (sim.param().pairlist.grid == 1){
    pa = new Grid_Pairlist_Algorithm();
  }
  else if (sim.param().pairlist.grid == 2){
    pa = new VGrid_Pairlist_Algorithm();
  }
  
#if defined(OMP)

  Nonbonded_Interaction * ni = new OMP_Nonbonded_Interaction(pa);

#elif defined(XXMPI)

  Nonbonded_Interaction * ni;

  if (sim.mpi){
    int rank = MPI::COMM_WORLD.Get_rank();
    if (rank == 0)
      ni = new MPI_Nonbonded_Master(pa);
    else
      ni = new MPI_Nonbonded_Slave(pa);
  }
  else
    ni = new Nonbonded_Interaction(pa);

#else

  Nonbonded_Interaction * ni = new Nonbonded_Interaction(pa);

#endif

  // standard LJ parameter
  if (sim.param().force.interaction_function ==
      simulation::lj_crf_func)
    it.read_lj_parameter(ni->parameter().lj_parameter());
  // and coarse-grained parameter
  if (sim.param().force.interaction_function ==
      simulation::cgrain_func)
    it.read_cg_parameter(ni->parameter().cg_parameter());
  
  // check if DUM really has no LJ interactons.
  pa->set_parameter(&ni->parameter());
  if (!sim.param().force.nonbonded_vdw && sim.param().force.interaction_function ==
      simulation::lj_crf_func){
    unsigned int dum = topo.iac(0); // has been previously set to DUM.

    bool error = false;
    for(unsigned int i = 0; i < dum; i++) {
      interaction::lj_parameter_struct& s = ni->parameter().lj_parameter()[i][dum];
      if (s.c6 != 0.0 || s.c12 != 0.0 || s.cs6 != 0.0 || s.cs12 != 0.0 )
        error = true;
    }
    for(unsigned int i = dum + 1; i < topo.atom_names().size(); i++)  {
      interaction::lj_parameter_struct& s = ni->parameter().lj_parameter()[dum][i];
      if (s.c6 != 0.0 || s.c12 != 0.0 || s.cs6 != 0.0 || s.cs12 != 0.0 )
        error = true;
    }
      
    if (error) {
      io::messages.add("Dummy atomtype (DUM) in topology has Lennard-Jones "
                       "interactions.", "topology", io::message::error);
    }
  }
  ff.push_back(ni);

  if (!quiet){
    os
      << "\tshortrange cutoff      : "
      << sim.param().pairlist.cutoff_short << "\n"
      << "\tlongrange cutoff       : "
      << sim.param().pairlist.cutoff_long << "\n"
      << "\treactionfield cutoff   : "
      << sim.param().longrange.rf_cutoff << "\n"
      << "\tepsilon                : "
      << sim.param().longrange.epsilon << "\n"
      << "\treactionfield epsilon  : "
      << sim.param().longrange.rf_epsilon << "\n"
      << "\tkappa                  : "
      << sim.param().longrange.rf_kappa << "\n"
      << "\tpairlist creation every "
      << sim.param().pairlist.skip_step
      << " steps\n\n";
  }
  
  return 0;
}
