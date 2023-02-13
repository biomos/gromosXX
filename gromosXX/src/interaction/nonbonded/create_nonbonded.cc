/**
 * @file create_nonbonded.cc
 * create the nonbonded interaction.
 */
#ifdef XXMPI
#include <mpi.h>
#endif

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"
#include "../../interaction/nonbonded/interaction/nonbonded_parameter.h"

#include "../../interaction/forcefield/forcefield.h"

#include "../../configuration/energy.h"
#include "../../interaction/nonbonded/interaction/storage.h"

#include "../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../../interaction/nonbonded/pairlist/standard_pairlist_algorithm.h"
#include "../../interaction/nonbonded/pairlist/extended_grid_pairlist_algorithm.h"
#include "../../interaction/nonbonded/pairlist/grid_cell_pairlist.h"

#include "../../interaction/nonbonded/interaction/nonbonded_outerloop.h"
#include "../../interaction/nonbonded/interaction/nonbonded_set.h"

#include "../../interaction/nonbonded/interaction/nonbonded_interaction.h"
#include "../../interaction/nonbonded/interaction/omp_nonbonded_interaction.h"
#include "../../interaction/nonbonded/interaction/mpi_nonbonded_master.h"
#include "../../interaction/nonbonded/interaction/mpi_nonbonded_slave.h"

#include "../../io/ifp.h"

#include "../../interaction/nonbonded/create_nonbonded.h"

#include "../../util/debug.h"

#include "../../simulation/parameter.h"

#ifdef HAVE_HOOMD
#include <HOOMD_GROMOSXX_interface.h>
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
  //if (sim.param().force.nonbonded_vdw == 0 &&
  //    sim.param().force.nonbonded_crf == 0) return 0;
  
  if(!quiet){
    
    os << "\t" << setw(20) << left << "nonbonded force";
    if (sim.param().force.nonbonded_vdw)
      os << setw(30) << left << "van-der-Waals";
    if (sim.param().force.nonbonded_crf)
      os << setw(30) << left << "Coulomb-reaction-field";
    os << "\n";

    if (sim.param().sasa.switch_sasa) {
      os << "\t" << "SASA Interaction" << "\n";
      if (sim.param().sasa.switch_volume)
        os << "\t" << "VOL Interaction" << "\n";
    }
   
    if (sim.param().pairlist.grid == 0)
      os << "\t" << setw(20) << left << "Pairlist Algorithm" << setw(30) 
	   << left << "Standard Pairlist Algorithm" << right << "\n";
    else if (sim.param().pairlist.grid == 1){
      os << "\t" << setw(20) << left << "Pairlist Algorithm" << setw(30) 
	   << left << "Grid Pairlist Algorithm" << right << "\n";
    }
    else if (sim.param().pairlist.grid == 2){
      os << "\t" << setw(20) << left << "Pairlist Algorithm" << setw(30) 
	   << left << "Grid-Cell Pairlist Algorithm" << right << "\n";
    }
    
    if (sim.param().polarise.cos) {
      os << "\t" << setw(20) << left << "Polarisation enabled" << right << "\n";
      os << "\t" << setw(20) << left << "Electric Field";
      switch(sim.param().polarise.efield_site) {
        case simulation::ef_atom :
          os << setw(30) << left << "calculated at atom position" << right << "\n";
          break;
        case simulation::ef_cos :
          os << setw(30) << left << "calculated at charge-on-spring position" << right << "\n";
          break;
        default :
          os << setw(30) << left << "calculated at unknown site" << right << "\n";
      }
      if (sim.param().polarise.damp) {
        os << "\t" << setw(20) << left
           << "Polarisation is damped using parameters from topology" << right << "\n";
      }
        
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

  Pairlist_Algorithm * pa = nullptr;
#ifdef HAVE_HOOMD
  // some sanity checks
  #ifdef OMP
  #error "HOOMD yet not ready for OMP. Might break!"
  #endif
  #ifdef XXMPI
  #error "HOOMD not yet ready for MPI. Might break!"
  #endif
  if (sim.param().hoomd.processor != simulation::unknown) { // hoomd pairlist
    pa = new HOOMD_Pairlist_Algorithm(sim); 
  } else { // gromosxx pairlist
#endif
    if (sim.param().pairlist.grid == 0) {
      pa = new Standard_Pairlist_Algorithm();
    } else if (sim.param().pairlist.grid == 1) {
      pa = new Extended_Grid_Pairlist_Algorithm();
    } else if (sim.param().pairlist.grid == 2) {
      pa = new Grid_Cell_Pairlist(topo, sim);
    }  else {
      io::messages.add("unkown pairlist algorithm.", "create_nonbonded",
             io::message::error);
      return 1;
	}
#ifdef HAVE_HOOMD
  }
#endif
   
#if defined(OMP)
  Nonbonded_Interaction * ni = new OMP_Nonbonded_Interaction(pa);
#elif defined(XXMPI)
  Nonbonded_Interaction * ni;

  if (sim.mpi){
    if (sim.mpiControl().threadID == sim.mpiControl().masterID)
      ni = new MPI_Nonbonded_Master(pa);
    else
      ni = new MPI_Nonbonded_Slave(pa);
  }
  else{
    ni = new Nonbonded_Interaction(pa);
  }
#else
  Nonbonded_Interaction * ni = new Nonbonded_Interaction(pa);
#endif

  // standard LJ parameter
  if (sim.param().force.interaction_function ==
      simulation::lj_crf_func ||
      sim.param().force.interaction_function ==
      simulation::lj_shifted_crf_corr_func ||
      sim.param().force.interaction_function ==
      simulation::cggromos_func ||
      sim.param().force.interaction_function ==
      simulation::pol_lj_crf_func ||
      sim.param().force.interaction_function ==
      simulation::pol_off_lj_crf_func ||
      sim.param().force.interaction_function ==
      simulation::lj_ls_func || 
      sim.param().force.interaction_function ==
      simulation::lj_func
          )
    it.read_lj_parameter(ni->parameter().lj_parameter());
  // and coarse-grained parameter (MARTINI model)
  if (sim.param().force.interaction_function ==
      simulation::cgrain_func)
    it.read_cg_parameter(ni->parameter().cg_parameter());
  
  // check if DUM really has no LJ interactons.
  pa->set_parameter(&ni->parameter());
  if ((!sim.param().force.nonbonded_vdw && sim.param().force.interaction_function ==
          simulation::lj_crf_func) ||
          (!sim.param().force.nonbonded_vdw && sim.param().force.interaction_function ==
          simulation::lj_shifted_crf_corr_func) ||
          (!sim.param().force.nonbonded_vdw && sim.param().force.interaction_function ==
          simulation::cggromos_func) ||
          (!sim.param().force.nonbonded_vdw && sim.param().force.interaction_function ==
          simulation::pol_lj_crf_func ) || (!sim.param().force.nonbonded_vdw
          && sim.param().force.interaction_function == simulation::pol_off_lj_crf_func )) {
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

  // scaling factor for electrostatic 1,4-interactions
  if (sim.param().amber.amber) {
    ni->parameter().set_coulomb_scaling(sim.param().amber.coulomb_scaling);
  } else {
    ni->parameter().set_coulomb_scaling(1.0);
  }

  ff.push_back(ni);

  if (!quiet){
    os
      << "\tshortrange cutoff      : "
      << sim.param().pairlist.cutoff_short << "\n"
      << "\tlongrange cutoff       : "
      << sim.param().pairlist.cutoff_long << "\n"
      << "\tpairlist creation every "
      << sim.param().pairlist.skip_step
      << " steps\n\n";
    
    if (sim.param().force.interaction_function == simulation::lj_crf_func ||
        sim.param().force.interaction_function == simulation::lj_shifted_crf_corr_func ||
        sim.param().force.interaction_function == simulation::cgrain_func || 
        sim.param().force.interaction_function == simulation::pol_lj_crf_func ||
        sim.param().force.interaction_function == simulation::pol_off_lj_crf_func) {
      os << "\tREACTION FIELD PARAMETERS\n"
              << "\treactionfield cutoff   : "
              << sim.param().nonbonded.rf_cutoff << "\n"
              << "\tepsilon                : "
              << sim.param().nonbonded.epsilon << "\n"
              << "\treactionfield epsilon  : "
              << sim.param().nonbonded.rf_epsilon << "\n"
              << "\tkappa                  : "
              << sim.param().nonbonded.rf_kappa << "\n";
    }
    
     if (sim.param().force.interaction_function == simulation::lj_ls_func) {
      os << "\tLATTICE SUM PARAMETERS\n"
              << "\tcharge shaping function         : "
              << sim.param().nonbonded.ls_charge_shape << "\n"
              << "\tcharge shaping function width   : "
              << sim.param().nonbonded.ls_charge_shape_width << "\n"
              << "\tA2 calculation method           : "
              << sim.param().nonbonded.ls_calculate_a2 << "\n"
              << "\tA2 relative tolerance           : "
              << sim.param().nonbonded.ls_a2_tolerance << "\n"
              << "\tLS permittivity                 : "
              << sim.param().nonbonded.ls_epsilon << "\n";
    }
    if (sim.param().nonbonded.method == simulation::el_ewald) {
      os << "\tMax. absolute Ewald k component : \n"
              << "\t" << std::setw(5) << sim.param().nonbonded.ewald_max_k_x
              << std::setw(5) << sim.param().nonbonded.ewald_max_k_y
              << std::setw(5) << sim.param().nonbonded.ewald_max_k_z << "\n"
              << "\tEwald k space cutoff            : "
              << sim.param().nonbonded.ewald_kspace_cutoff << "\n";
    }
    if (sim.param().nonbonded.method == simulation::el_p3m) {
      os << "\tP3M number of grid point along axes  : \n"
              << "\t" << std::setw(5) << sim.param().nonbonded.p3m_grid_points_x
              << std::setw(5) << sim.param().nonbonded.p3m_grid_points_y
              << std::setw(5) << sim.param().nonbonded.p3m_grid_points_z << "\n";
      if (sim.param().multicell.multicell) {
        os << "\tgrid dimensions are going to be scaled according to MULTICELL information.\n";
      }
      os << "\tP3M assignment function order   : "
              << sim.param().nonbonded.p3m_charge_assignment << "\n"
              << "\tP3M finite differences order    : "
              << sim.param().nonbonded.p3m_finite_differences_operator << "\n"
              << "\tP3M number of mesh alias vecs   : "
              << sim.param().nonbonded.p3m_mesh_alias << "\n"
              << "\tP3M accuracy evaluation         : "
              << sim.param().nonbonded.accuracy_evaluation << "\n"
              << "\tP3M RMS force error threshold   : "
              << sim.param().nonbonded.influence_function_rms_force_error << "\n";
    }
  }
  
  return 0;
}
