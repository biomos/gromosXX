/**
 * @file configuration.cc
 * methods definition
 */

#ifdef XXMPI
#include <mpi.h>
#endif

#include "../stdheader.h"

#include "../configuration/configuration_global.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../configuration/configuration.h"
#include "../configuration/mesh.h"
#include "../configuration/influence_function.h"
#include "../simulation/simulation.h"
#include "../simulation/multibath.h"
#include "../simulation/parameter.h"

#include "../math/periodicity.h"
#include "../math/boundary_checks.h"
#include "../util/template_split.h"

#include "configuration.h"

#ifdef OMP
#include <omp.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE configuration

double configuration_ver = 0.10;

/**
 * Constructor
 */
configuration::Configuration::Configuration() {
  m_current = &m_state1;
  m_old = &m_state2;

  current().virial_tensor = 0.0;
  old().virial_tensor = 0.0;

  current().kinetic_energy_tensor = 0.0;
  old().kinetic_energy_tensor = 0.0;

  current().pressure_tensor = 0.0;
  old().pressure_tensor = 0.0;

  current().sasa_tot = 0.0;
  old().sasa_tot = 0.0;
  current().sasa_buriedvol_tot = 0.0;
  old().sasa_buriedvol_tot = 0.0;

  for (unsigned int k = 0; k < special().eds.virial_tensor_endstates.size(); ++k) {
    special().eds.virial_tensor_endstates[k] = 0.0;
  }
  
  special().shake_failure_occurred = false;
}

/**
 * copy constructor
 */
configuration::Configuration::Configuration
(
 configuration::Configuration const & conf
 )
{
  m_current = &m_state1;
  m_old = &m_state2;

  current().virial_tensor = conf.current().virial_tensor;
  old().virial_tensor = conf.old().virial_tensor;

  current().kinetic_energy_tensor = conf.current().kinetic_energy_tensor;
  old().kinetic_energy_tensor = conf.old().kinetic_energy_tensor;

  current().pressure_tensor = conf.current().pressure_tensor;
  old().pressure_tensor = conf.old().pressure_tensor;
  
  for (unsigned int k = 0; k < special().eds.virial_tensor_endstates.size(); ++k) {
    special().eds.virial_tensor_endstates[k] =
            conf.special().eds.virial_tensor_endstates[k];
  }


  current().pos = conf.current().pos;
  old().pos = conf.old().pos;
  current().posV = conf.current().posV;
  old().posV = conf.old().posV;
  current().vel = conf.current().vel;
  old().vel = conf.old().vel;
  current().force = conf.current().force;
  current().constraint_force = conf.current().constraint_force;
  current().cgrad = conf.current().cgrad;
  old().force = conf.old().force;
  old().constraint_force = conf.current().constraint_force;
  old().cgrad = conf.old().cgrad;
  current().stochastic_integral = conf.current().stochastic_integral;
  old().stochastic_integral = conf.old().stochastic_integral;
  current().stochastic_seed = conf.current().stochastic_seed;
  old().stochastic_seed = conf.old().stochastic_seed;
  
  current().box = conf.current().box;
  old().box = conf.old().box;
 
  current().energies = conf.current().energies;
  old().energies = conf.old().energies;
  current().averages = conf.current().averages;
  old().averages = conf.old().averages;
  
  current().perturbed_energy_derivatives =
    conf.current().perturbed_energy_derivatives;
  old().perturbed_energy_derivatives =
    conf.old().perturbed_energy_derivatives;

  current().sasa_area = conf.current().sasa_area;
  old().sasa_area = conf.old().sasa_area;
  current().sasa_buriedvol = conf.current().sasa_buriedvol;
  old().sasa_buriedvol = conf.old().sasa_buriedvol;
  current().gvol = conf.current().gvol;
  old().gvol = conf.old().gvol;
  current().dgvol = conf.current().dgvol;
  old().dgvol = conf.old().dgvol;

  current().sasa_tot = conf.current().sasa_tot;
  old().sasa_tot = conf.old().sasa_tot;
  current().sasa_buriedvol_tot = conf.current().sasa_buriedvol_tot;
  old().sasa_buriedvol_tot = conf.old().sasa_buriedvol_tot;

  special().dihangle_trans.dihedral_angle_minimum = conf.special().dihangle_trans.dihedral_angle_minimum;
  special().dihangle_trans.old_minimum = conf.special().dihangle_trans.old_minimum;
  special().dihangle_trans.resid = conf.special().dihangle_trans.resid;
  special().dihangle_trans.i = conf.special().dihangle_trans.i;
  special().dihangle_trans.j = conf.special().dihangle_trans.j;
  special().dihangle_trans.k = conf.special().dihangle_trans.k;
  special().dihangle_trans.l = conf.special().dihangle_trans.l;
  
  special().umbrellas = conf.special().umbrellas;
  special().flexible_constraint = conf.special().flexible_constraint;
  
  special().jvalue_av = conf.special().jvalue_av;
  special().jvalue_curr = conf.special().jvalue_curr;
  special().jvalue_epsilon = conf.special().jvalue_epsilon;
  
  special().distanceres.av = conf.special().distanceres.av;
  special().distanceres.energy = conf.special().distanceres.energy;
  special().distanceres.d = conf.special().distanceres.d;
  
  special().distancefield.ngrid = conf.special().distancefield.ngrid;
  special().distancefield.distance = conf.special().distancefield.distance;
  special().distancefield.dist = conf.special().distancefield.dist;
  special().distancefield.energy = conf.special().distancefield.energy;
  special().distancefield.energy_deriv = conf.special().distancefield.energy_deriv;
  
  special().angleres.energy = conf.special().angleres.energy;
  special().angleres.d = conf.special().angleres.d;

  special().dihedralres.energy = conf.special().dihedralres.energy;
  special().dihedralres.d = conf.special().dihedralres.d;
  
  special().pscale = conf.special().pscale;

  special().orderparamres.S2_avg = conf.special().orderparamres.S2_avg;
  special().orderparamres.energy = conf.special().orderparamres.energy;
  special().orderparamres.Q_avg = conf.special().orderparamres.Q_avg;
  special().orderparamres.D_avg = conf.special().orderparamres.D_avg;
  special().orderparamres.Q_winavg = conf.special().orderparamres.Q_winavg;
  special().orderparamres.D_winavg = conf.special().orderparamres.D_winavg;

  special().rdc = conf.special().rdc;
 
  special().rottrans_constr = conf.special().rottrans_constr;

  // if this works just like this, why do we need to explicitly copy the virial tensor?
  special().eds = conf.special().eds;
  
  special().lattice_shifts = conf.special().lattice_shifts;
  
  special().shake_failure_occurred = conf.special().shake_failure_occurred;
  
  special().force_groups = conf.special().force_groups;
  
  boundary_type = conf.boundary_type;
}

/**
 * operator equal
 */
configuration::Configuration & configuration::Configuration::operator=
(
 configuration::Configuration const & conf
 )
{
  m_current = &m_state1;
  m_old = &m_state2;

  current().virial_tensor = conf.current().virial_tensor;
  old().virial_tensor = conf.old().virial_tensor;

  current().kinetic_energy_tensor = conf.current().kinetic_energy_tensor;
  old().kinetic_energy_tensor = conf.old().kinetic_energy_tensor;

  current().pressure_tensor = conf.current().pressure_tensor;
  old().pressure_tensor = conf.old().pressure_tensor;

  for (unsigned int k = 0; k < special().eds.virial_tensor_endstates.size(); ++k) {
    special().eds.virial_tensor_endstates[k] =
            conf.special().eds.virial_tensor_endstates[k];
  }

  
  current().pos = conf.current().pos;
  old().pos = conf.old().pos;
  current().posV = conf.current().posV;
  old().posV = conf.old().posV;
  current().vel = conf.current().vel;
  old().vel = conf.old().vel;
  current().force = conf.current().force;
  old().force = conf.old().force;
  current().cgrad = conf.current().cgrad;
  old().cgrad = conf.old().cgrad;
  current().stochastic_integral = conf.current().stochastic_integral;
  old().stochastic_integral = conf.old().stochastic_integral;
  current().stochastic_seed = conf.current().stochastic_seed;
  old().stochastic_seed = conf.old().stochastic_seed;
  
  current().box = conf.current().box;
  old().box = conf.old().box;
  
  current().energies = conf.current().energies;
  old().energies = conf.old().energies;
  current().averages = conf.current().averages;
  old().averages = conf.old().averages;
  
  current().perturbed_energy_derivatives =
    conf.current().perturbed_energy_derivatives;
  old().perturbed_energy_derivatives =
    conf.old().perturbed_energy_derivatives;

  current().sasa_area = conf.current().sasa_area;
  old().sasa_area = conf.old().sasa_area;
  current().sasa_buriedvol = conf.current().sasa_buriedvol;
  old().sasa_buriedvol = conf.old().sasa_buriedvol;
  current().gvol = conf.current().gvol;
  old().gvol = conf.old().gvol;
  current().dgvol = conf.current().dgvol;
  old().dgvol = conf.old().dgvol;

  current().sasa_tot = conf.current().sasa_tot;
  old().sasa_tot = conf.old().sasa_tot;
  current().sasa_buriedvol_tot = conf.current().sasa_buriedvol_tot;
  old().sasa_buriedvol_tot = conf.old().sasa_buriedvol_tot;

  special().dihangle_trans.dihedral_angle_minimum = conf.special().dihangle_trans.dihedral_angle_minimum;
  special().dihangle_trans.old_minimum = conf.special().dihangle_trans.old_minimum;
  special().dihangle_trans.resid = conf.special().dihangle_trans.resid;
  special().dihangle_trans.i = conf.special().dihangle_trans.i;
  special().dihangle_trans.j = conf.special().dihangle_trans.j;
  special().dihangle_trans.k = conf.special().dihangle_trans.k;
  special().dihangle_trans.l = conf.special().dihangle_trans.l;
  
  special().umbrellas = conf.special().umbrellas;
  special().flexible_constraint = conf.special().flexible_constraint;
  
  special().jvalue_av = conf.special().jvalue_av;
  special().jvalue_curr = conf.special().jvalue_curr;
  special().jvalue_epsilon = conf.special().jvalue_epsilon;
  
  special().distanceres.av = conf.special().distanceres.av;
  special().distanceres.energy = conf.special().distanceres.energy;
  special().distanceres.d = conf.special().distanceres.d;
 
  special().distancefield.ngrid = conf.special().distancefield.ngrid;
  special().distancefield.distance = conf.special().distancefield.distance;
  special().distancefield.dist = conf.special().distancefield.dist;
  special().distancefield.energy = conf.special().distancefield.energy;
  special().distancefield.energy_deriv = conf.special().distancefield.energy_deriv;
  
  special().angleres.energy = conf.special().angleres.energy;
  special().angleres.d = conf.special().angleres.d;
  
  special().dihedralres.energy = conf.special().dihedralres.energy;
  special().dihedralres.d = conf.special().dihedralres.d;
 
  special().rdc = conf.special().rdc;

  special().pscale = conf.special().pscale;
  
  special().rottrans_constr = conf.special().rottrans_constr;

  special().eds = conf.special().eds;
  
  special().lattice_shifts = conf.special().lattice_shifts;
  
  special().shake_failure_occurred = conf.special().shake_failure_occurred;
  
  special().force_groups = conf.special().force_groups;
  
  boundary_type = conf.boundary_type;

  return *this;
}

void configuration::Configuration::init(topology::Topology const & topo,
					simulation::Parameter & param,
					bool gather)
{
  // resize the energy arrays
  const unsigned int num = unsigned(topo.energy_groups().size());
  const unsigned int numb = unsigned(param.multibath.multibath.size());
  // ANITA
  const unsigned int numl = unsigned(param.precalclam.nr_lambdas); //
  
  DEBUG(5, "number of energy groups: " << num 
	<< "\nnumber of baths: " << numb
        << "\nnumber of lambdas: " << numl); // ANITA

  DEBUG(5, "ANITA resizing energies, configuration::init"); 
  current().energies.resize(num, numb,numl);
  old().energies.resize(num, numb,numl);

//  current().energies.resize(num, numb);
//  old().energies.resize(num, numb); // ANITA
  if (param.force.force_groups) {
    special().force_groups.resize(num, 
            std::vector<math::VArray>(num, math::VArray(
            topo.num_atoms(), math::Vec(0.0, 0.0, 0.0))));
  }

  // resize sasa vectors
  const unsigned int num_sasa_atoms = topo.sasa_parameter().size();
  DEBUG(5, "Number of sasa atoms: " << num_sasa_atoms);

  current().sasa_area.resize(num_sasa_atoms, 0.0);
  old().sasa_area.resize(num_sasa_atoms, 0.0);
  current().sasa_buriedvol.resize(num_sasa_atoms, 0.0);
  old().sasa_buriedvol.resize(num_sasa_atoms, 0.0);
  current().gvol.resize(num_sasa_atoms, 0.0);
  old().gvol.resize(num_sasa_atoms, 0.0);
  current().dgvol.resize(num_sasa_atoms, 0.0);
  old().dgvol.resize(num_sasa_atoms, 0.0);

  // check whether this can really stay here! see resize function below
  special().eds.force_endstates.resize(param.eds.numstates);
  for (unsigned int i = 0; i < special().eds.force_endstates.size(); i++){
    special().eds.force_endstates[i].resize(topo.num_atoms());
  }  
  special().eds.virial_tensor_endstates.resize(param.eds.numstates);
  current().energies.eds_vi.resize(param.eds.numstates);
  current().energies.eds_eir.resize(param.eds.numstates);
  current().energies.eds_vi_special.resize(param.eds.numstates);
  old().energies.eds_vi.resize(param.eds.numstates);
  old().energies.eds_eir.resize(param.eds.numstates);
  old().energies.eds_vi_special.resize(param.eds.numstates);
  
  current().energies.ewarn(param.ewarn.limit);
  old().energies.ewarn(param.ewarn.limit);

  //ANITA
//  current().perturbed_energy_derivatives.resize(num, numb);
//  old().perturbed_energy_derivatives.resize(num, numb);
  DEBUG(5, "ANITA resizing perturbed energies, configuration::init");
  current().perturbed_energy_derivatives.resize(num, numb,numl);
  old().perturbed_energy_derivatives.resize(num, numb,numl); //

  DEBUG(5, "ANITA resizing averages , configuration::init");
  current().averages.resize(topo, *this, param);
  old().averages.resize(topo, *this, param);

  // possibly resize the dihedral angle monitoring array
  // initialize or set to such a value that it is recalculated in
  // the first step, initialization would require an additional function
  // which can be done only after reading of the coordinates. Would be a
  // nicer solution, but also requires the parameters...
  if(param.print.monitor_dihedrals){
    special().dihangle_trans.dihedral_angle_minimum.resize
      (topo.solute().dihedrals().size(), 4*math::Pi);
    special().dihangle_trans.old_minimum.resize
      (topo.solute().dihedrals().size(), 0.0);
    special().dihangle_trans.resid.resize
      (topo.solute().dihedrals().size(), 0);
    special().dihangle_trans.i.resize
      (topo.solute().dihedrals().size(), 0);
    special().dihangle_trans.j.resize
      (topo.solute().dihedrals().size(), 0);
    special().dihangle_trans.k.resize
      (topo.solute().dihedrals().size(), 0);
    special().dihangle_trans.l.resize
      (topo.solute().dihedrals().size(), 0);
  }
  
  if (param.constraint.solute.algorithm == simulation::constr_flexshake &&
      special().flexible_constraint.flexible_vel.size() == 0){

    special().flexible_constraint.flexible_vel.resize(topo.solute().distance_constraints().size() +
				  topo.perturbed_solute().distance_constraints().size());

    special().flexible_constraint.flexible_ekin.resize(numb);
  }

  if(param.nemd.nemd != simulation::nemd_off){
    special().nemd_conf.Px = 0.0; //To accumulate the momemtum
    special().nemd_conf.counter = 0; 
    unsigned int num_grid = 2 * param.nemd.slabnum;
    special().nemd_conf.stored_data_per_bin.clear();
    for(unsigned int i = 0; i < num_grid; ++i) {
      special().nemd_conf.stored_data_per_bin.push_back(0.0);
    }
    for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
      special().nemd_conf.vel_per_atom.push_back(0.0);
      special().nemd_conf.dvx_per_atom.push_back(0.0);
    }
  }
  
  
  // resize the arrays
  // to make scripting easier...
  resize(topo.num_atoms());

  // gather the molecules!
  // check box size

  // mc: bugfix: chargegroups should be gathered
  //             problem if submolecules are set to 1 atom
  //             (for whatever esoteric reasons)

  if (gather){
    switch(boundary_type){
      case math::vacuum:
	break;
      case math::rectangular:
	{
	  math::Periodicity<math::rectangular> periodicity(current().box);
	  // periodicity.gather_molecules_into_box(*this, topo);
	  periodicity.gather_chargegroups(*this, topo);
	  
	  break;
	}
      case math::truncoct:
      case math::triclinic:
	{
	  // NO CUTOFF CHECK -- IMPLEMENT!!!
	  math::Periodicity<math::triclinic> periodicity(current().box);
	  // periodicity.gather_molecules_into_box(*this, topo);
	  periodicity.gather_chargegroups(*this, topo);
	  
	  break;
	}
      default:
	std::cout << "wrong periodic boundary conditions!";
	io::messages.add("wrong PBC!", "In_Configuration", io::message::error);
    }
  }

  // check periodicity
  if (!param.multicell.multicell) {
    if (!math::boundary_check_cutoff(current().box, boundary_type, param.pairlist.cutoff_long)) {
      io::messages.add("box is too small: not twice the cutoff!",
              "configuration",
              io::message::error);
    }
  }

  if (boundary_type != math::vacuum){
    if (param.centreofmass.remove_rot){
      io::messages.add("disabling removing of centre of mass rotation (PBC)",
		       "configuration",
		       io::message::notice);
      param.centreofmass.remove_rot = false;
    }
  }
}



/**
 * set the number of atoms.
 */
void configuration::Configuration::resize(unsigned int s)
{
  DEBUG(7, "Configuration resize: " << s);
  
  current().resize(s);
  old().resize(s);
  
  special().lattice_shifts.resize(s);
}

/**
 * set the number of atoms.
 * using resizeAndPreserve. Therefore
 * you can enlarge the system (or shrink it)
 * while keeping all existing positions/velocities/...
 * a faster version would be just resize, but then
 * the arrays contain garbage...
 * the energies have to be sized separately!
 */
void configuration::Configuration::state_struct::resize(unsigned int s)
{
  DEBUG(7, "state struct resize: " << s);

  pos.resize(s, math::Vec(0.0, 0.0, 0.0));
  posV.resize(s, math::Vec(0.0, 0.0, 0.0));
  vel.resize(s, math::Vec(0.0, 0.0, 0.0));
  force.resize(s, math::Vec(0.0, 0.0, 0.0));
  cgrad.resize(s, math::Vec(0.0, 0.0, 0.0));
  constraint_force.resize(s, math::Vec(0.0, 0.0, 0.0));
  stochastic_integral.resize(s, math::Vec(0.0, 0.0, 0.0));
}

void configuration::Configuration::lattice_sum_struct::init(topology::Topology const & topo,
        simulation::Simulation & sim) {
  DEBUG(1,"Lattice Sum initialization.");
  simulation::Parameter & param = sim.param();
#ifdef OMP
  int tid = 0, size = 0;
#pragma omp parallel private(tid)
  {
    tid = omp_get_thread_num();
    if (tid == 0) {
      size = omp_get_num_threads();
    }
  }
  FFTW3(init_threads());
  FFTW3(plan_with_nthreads(size));
  sim.openmp = true;
#endif
  // get the k space
  if (param.nonbonded.method == simulation::el_ewald) {
    kspace.reserve(param.nonbonded.ewald_max_k_x *
            param.nonbonded.ewald_max_k_y *
            param.nonbonded.ewald_max_k_z);
  }
  
  if (param.nonbonded.method == simulation::el_p3m) {
    const unsigned int Nx = param.nonbonded.p3m_grid_points_x;
    const unsigned int Ny = param.nonbonded.p3m_grid_points_y;
    const unsigned int Nz = param.nonbonded.p3m_grid_points_z;

    const bool do_a2t =
          param.nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact ||
          param.nonbonded.ls_calculate_a2 == simulation::ls_a2t_exact_a2_numerical ||
          param.nonbonded.ls_calculate_a2 == simulation::ls_a2t_ave_a2_numerical;

#ifdef XXMPI
    if (sim.mpi) {
      int rank = sim.mpiControl().threadID;
      int num_threads = sim.mpiControl().numberOfThreads;
      MPI_Comm com = sim.mpiControl().comm;
      
      const int cache_size = std::max(param.nonbonded.p3m_charge_assignment - 1,
              param.nonbonded.p3m_finite_differences_operator);
      
      charge_density = new configuration::ParallelMesh(num_threads, rank, cache_size, com);

      potential = new configuration::ParallelMesh(num_threads, rank, cache_size, com);
      electric_field.x = new configuration::ParallelMesh(num_threads, rank, cache_size, com);
      electric_field.y = new configuration::ParallelMesh(num_threads, rank, cache_size, com);
      electric_field.z = new configuration::ParallelMesh(num_threads, rank, cache_size, com);

      if (do_a2t)
        squared_charge = new configuration::ParallelMesh(num_threads, rank, cache_size, com);
      
      ((configuration::ParallelMesh*)charge_density)->resize(Nx, Ny, Nz);
      ((configuration::ParallelMesh*)potential)->resize(Nx, Ny, Nz);
      ((configuration::ParallelMesh*)electric_field.x)->resize(Nx, Ny, Nz);
      ((configuration::ParallelMesh*)electric_field.y)->resize(Nx, Ny, Nz);
      ((configuration::ParallelMesh*)electric_field.z)->resize(Nx, Ny, Nz);
      if (do_a2t)
        ((configuration::ParallelMesh*)squared_charge)->resize(Nx, Ny, Nz);
    } else {
#endif
      charge_density = new configuration::Mesh();
      potential = new configuration::Mesh();
      electric_field.x = new configuration::Mesh();
      electric_field.y = new configuration::Mesh();
      electric_field.z = new configuration::Mesh();
      if (do_a2t)
        squared_charge = new configuration::Mesh();
      charge_density->resize(Nx, Ny, Nz);
      potential->resize(Nx, Ny, Nz);
      electric_field.x->resize(Nx, Ny, Nz);
      electric_field.y->resize(Nx, Ny, Nz);
      electric_field.z->resize(Nx, Ny, Nz);
      if (do_a2t)
        squared_charge->resize(Nx, Ny, Nz);
#ifdef XXMPI
    }
#endif
    
    influence_function.init(param);
  }

  // reset the A term
  a2_tilde = 0.0;


}
namespace configuration
{
  std::ostream &operator<<(std::ostream &os, Configuration &conf)
  {
    os << "a configuration";
    return os;
  }
}


bool configuration::Configuration::check(topology::Topology const & topo, simulation::Simulation & sim) {
  int error = 0;
  
  // check the positions if nonbonded forces are computed
  if (sim.param().force.nonbonded_crf || sim.param().force.nonbonded_vdw) {
    SPLIT_MY_BOUNDARY(boundary_type, check_positions, topo, error);
    SPLIT_MY_BOUNDARY(boundary_type, check_excluded_positions, topo, sim);
  }
  // check the positions if (improper) dihedrals are computed
  if(sim.param().force.improper){
    SPLIT_MY_BOUNDARY(boundary_type, check_dihedrals, topo.solute().improper_dihedrals(), error);
  }
  if(sim.param().force.dihedral){
    SPLIT_MY_BOUNDARY(boundary_type, check_dihedrals, topo.solute().dihedrals(), error);
  }   
  
  return error == 0;
}

template<math::boundary_enum B> 
void configuration::Configuration::check_positions(topology::Topology const & topo, int & error) const {  
  math::Periodicity<B> periodicity(current().box);
  const unsigned int num_pos = current().pos.size();
  const math::VArray & pos = current().pos;
  math::Vec r;
  for(unsigned int i = 0; i < num_pos; ++i) {
    for(unsigned int j = i + 1; j < num_pos; ++j) {
      periodicity.nearest_image(pos(i), pos(j), r);
      if (math::abs2(r) < math::epsilon) {
        // if they are excluded, it is a warning, if not it is an error
	if (topo.exclusion(i).is_excluded(j)){
	    std::ostringstream msg;
            msg << "Singularity: Atoms " << i+1 << " and " << j+1 << " are at the same "
		<< "position. They are excluded from nonbonded interactions.";
            io::messages.add(msg.str(), "Configuration", io::message::warning);
	} else {
	  ++error;
	  std::ostringstream msg;
	  msg << "Singularity: Atoms " << i+1 << " and " << j+1 << " at same position.";
	  io::messages.add(msg.str(), "Configuration", io::message::error);
	}
      }
      
    }
  }
}
template<math::boundary_enum B> 
void configuration::Configuration::check_dihedrals(std::vector<topology::four_body_term_struct> const & dihedrals, int & error) const { 
  math::Periodicity<B> periodicity(current().box);
  const math::VArray & pos = current().pos;
  math::Vec r;
  std::vector<topology::four_body_term_struct>::const_iterator d_it = dihedrals.begin(), d_to = dihedrals.end();
  for( ; d_it != d_to; ++d_it){
    periodicity.nearest_image(pos(d_it->i), pos(d_it->j), r);
    if (math::abs2(r) < math::epsilon) {
      ++error;
      std::ostringstream msg;
      msg << "Singularity: Atoms " << d_it->i+1 << " and " << d_it->j+1 << " at same position. Cannot compute (improper) dihedral.";
      io::messages.add(msg.str(), "Configuration", io::message::error);
    }
    periodicity.nearest_image(pos(d_it->i), pos(d_it->k), r);
    if (math::abs2(r) < math::epsilon) {
      ++error;
      std::ostringstream msg;
      msg << "Singularity: Atoms " << d_it->i+1 << " and " << d_it->k+1 << " at same position. Cannot compute (improper) dihedral.";
      io::messages.add(msg.str(), "Configuration", io::message::error);
    }
    periodicity.nearest_image(pos(d_it->i), pos(d_it->l), r);
    if (math::abs2(r) < math::epsilon) {
      ++error;
      std::ostringstream msg;
      msg << "Singularity: Atoms " << d_it->i+1 << " and " << d_it->l+1 << " at same position. Cannot compute (improper) dihedral.";
      io::messages.add(msg.str(), "Configuration", io::message::error);
    }
    periodicity.nearest_image(pos(d_it->j), pos(d_it->k), r);
    if (math::abs2(r) < math::epsilon) {
      ++error;
      std::ostringstream msg;
      msg << "Singularity: Atoms " << d_it->j+1 << " and " << d_it->k+1 << " at same position. Cannot compute (improper) dihedral.";
      io::messages.add(msg.str(), "Configuration", io::message::error);
    }
    periodicity.nearest_image(pos(d_it->j), pos(d_it->l), r);
    if (math::abs2(r) < math::epsilon) {
      ++error;
      std::ostringstream msg;
      msg << "Singularity: Atoms " << d_it->j+1 << " and " << d_it->l+1 << " at same position. Cannot compute (improper) dihedral.";
      io::messages.add(msg.str(), "Configuration", io::message::error);
    }
    periodicity.nearest_image(pos(d_it->k), pos(d_it->l), r);
    if (math::abs2(r) < math::epsilon) {
      ++error;
      std::ostringstream msg;
      msg << "Singularity: Atoms " << d_it->k+1 << " and " << d_it->l+1 << " at same position. Cannot compute (improper) dihedral.";
      io::messages.add(msg.str(), "Configuration", io::message::error);
    }

  }
}


template<math::boundary_enum B> 
void configuration::Configuration::check_excluded_positions(topology::Topology const & topo, simulation::Simulation & sim) {
  math::Periodicity<B> periodicity(current().box);
  const int num_solute = topo.num_solute_atoms();
  const math::VArray & pos = current().pos;
  math::Vec r;
  const double cutoff_2 = sim.param().pairlist.cutoff_short * sim.param().pairlist.cutoff_short;
  
  if (sim.param().pairlist.atomic_cutoff) {
    // loop over the solute atoms
    for(int a1 = 0 ; a1 < num_solute; ++a1) {
      for(int a2 = a1 + 1; a2 < num_solute; ++a2){
        // check if they are outside the inner cut off
        periodicity.nearest_image(pos(a1), pos(a2), r);
        const double d2 = math::abs2(r);
        if (d2 > cutoff_2) {
          // if yes, check if they are excluded
          if (topo.all_exclusion(a1).is_excluded(a2)) {
            //check if reeds is on and both atoms are perturbed - then subpress the warning
            if(sim.param().reeds.reeds){
                continue;
            }
            // if yes, issue warning!
            std::ostringstream msg;
            msg << "Warning: Atoms " << a1 << " and " << a2
                << " are excluded, but they are further apart than the"
                << " short cut-off radius distance. Any two atoms further"
                << " apart than this distance interact fully.";
            io::messages.add(msg.str(), "Configuration", io::message::warning);
          }
        }
      } 
    }
      
  } else {
    // first put the chargegroups into the box
    periodicity.put_chargegroups_into_box(*this, topo);

    // calculate cg cog's
    DEBUG(10, "calculating cg cog (" << topo.num_solute_chargegroups() << ")");
    math::VArray cg_cog;
    cg_cog.resize(topo.num_solute_chargegroups());

    // calculate solute center of geometries
    topology::Chargegroup_Iterator
      cg1 =   topo.chargegroup_begin();
    unsigned int i = 0, num_cg = topo.num_solute_chargegroups();
    
    for(i=0; i < num_cg; ++cg1, ++i){
      cg1.cog(pos, cg_cog(i));
    }

    // loop over the solute charge groups
    unsigned int idx_cg1 = 0, idx_cg2 = 0;
    for (idx_cg1 = 0; idx_cg1 < num_cg; idx_cg1++) {
      for (idx_cg2 = idx_cg1 + 1; idx_cg2 < num_cg; idx_cg2++) {
        // check if they are outside of inner cut off
        periodicity.nearest_image(cg_cog(idx_cg1), cg_cog(idx_cg2), r);
        const double d2 = math::abs2(r);
        if (d2 > cutoff_2) {
          // if yes, check if any of the atoms are excluded
            for (int a1 = topo.chargegroup(idx_cg1), a1_to = topo.chargegroup(idx_cg1 + 1);
                   a1 != a1_to; ++a1) {
            for (int a2 = topo.chargegroup(idx_cg2), a2_to = topo.chargegroup(idx_cg2 + 1);
                     a2 != a2_to; ++a2) {

              bool not_excluded_via_eds= !((sim.param().reeds.reeds || sim.param().eds.eds) && (topo.eds_perturbed_solute().atoms().count(a1)>=1 && topo.eds_perturbed_solute().atoms().count(a2)>=1));
              if (topo.all_exclusion(a1).is_excluded(a2) && not_excluded_via_eds) {
                // if yes, issue warning!
                std::ostringstream msg;
                msg << "Warning: Atoms " << a1 << " and " << a2
                    << " are excluded, but their respective charge groups "
                    << idx_cg1 << " and " << idx_cg2 << " are further apart than the"
                    << " short cut-off radius distance. Any two atoms whose charge"
                    << " groups are further apart than this distance interact fully.";
                io::messages.add(msg.str(), "Configuration", io::message::warning);
              }
            } // loop over atom of cg2
          } // loop over atom of cg1
        }
      }
    }
  }
}
