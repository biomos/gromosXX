/**
 * @file out_configuration.cc
 * definition of the Out_Configuration methods.
 */
#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <simulation/parameter.h>
#include <configuration/configuration.h>
#include <configuration/energy.h>

#include <math/periodicity.h>
#include <math/volume.h>
#include <math/transformation.h>

#include <io/print_block.h>
#include <io/argument.h>

#include "out_configuration.h"

#include <util/template_split.h>
#include <util/debug.h>
#include <limits>
#include <util/le_coordinate.h>
#include <util/umbrella_weight.h>
#include <util/bs_potentials.h>

#undef MODULE
#undef SUBMODULE

#define MODULE io
#define SUBMODULE configuration

// Energy trajectory version
// For details, see definition in out_configuration.cc
const std::string io::Out_Configuration::ene_version = "2021-04-12";

// declarations
static void _print_energyred_helper(std::ostream & os, configuration::Energy const &e);

static void _print_volumepressurered_helper(std::ostream &os,
        double mass,
        double const & phi, double const & theta, double const &psi,
        simulation::Multibath const & m,
        std::vector<double> const & s,
        configuration::Energy const & e,
        math::Box const & b,
        math::boundary_enum t,
        math::Matrix const & p,
        math::Matrix const & v,
        math::Matrix const & k);

io::Out_Configuration::Out_Configuration(std::string title,
        std::ostream & os)
: m_has_replica_traj(false),
m_output(os),
m_final(false),
m_replica(false),
m_write_special(false),
m_every_pos(0),
m_every_vel(0),
m_every_force(0),
m_every_energy(0),
m_every_free_energy(0),
m_every_blockaverage(0),
m_every_cos_pos(0),
m_every_jvalue(0),
m_every_xray(0),
m_every_disres(0),
m_every_disfieldres(0),
m_every_angres(0),
m_every_dihres(0),
m_every_dat(0),
m_every_leus(0),
m_every_bsleus(0),
m_every_dipole(0),
m_every_current(0),
m_every_adde(0),
m_every_nemd(0),
m_every_oparam(0),
m_every_rdc(0),
m_write_blockaverage_energy(false),
m_write_blockaverage_free_energy(false),
m_precision(9),
m_force_precision(9),
m_distance_restraint_precision(7),
m_disfield_restraint_precision(7),
m_angle_restraint_precision(7),
m_dihedral_restraint_precision(7),
m_width(15),
m_force_width(18),
m_title(title),
minimum_energy(std::numeric_limits<double>::max()) {
  _print_title(m_title, "output file", os);
}

/**
 * destructor
 */
io::Out_Configuration::~Out_Configuration() 
{
  // std::cout << "out_configuration destructor" << std::endl;

  if (m_every_pos) {
    m_pos_traj.flush();
    m_pos_traj.close();
  }

  if (m_final) {
    m_final_conf.flush();
    m_final_conf.close();
  }

  if (m_every_vel) {
    m_vel_traj.flush();
    m_vel_traj.close();
  }

  if (m_every_force) {
    m_force_traj.flush();
    m_force_traj.close();
  }

  if (m_every_energy) {
    m_energy_traj.flush();
    m_energy_traj.close();
  }

  if (m_has_replica_traj) {
    m_replica_traj.flush();
    m_replica_traj.close();
  }

  if (m_every_free_energy) {
    m_free_energy_traj.flush();
    m_free_energy_traj.close();
  }

  if (m_write_blockaverage_energy) {
    m_blockaveraged_energy.flush();
    m_blockaveraged_energy.close();
  }

  if (m_write_blockaverage_free_energy) {
    m_blockaveraged_free_energy.flush();
    m_blockaveraged_free_energy.close();
  }

  if (m_write_special) { 
    m_special_traj.flush();
    m_special_traj.close();
  }
}

void io::Out_Configuration::_print_title(std::string title,
        std::string name,
        std::ostream &os) {
  os << "TITLE\n\t"
          << title << "\n"
          << "\t" << name
          << "\nEND\n";
}

/*
     * Prints the ENEVERSION block for the (free) energy trajectories
     */
void io::Out_Configuration::_print_ene_version(std::ostream &os) {
  os << "ENEVERSION\n\t"
     << ene_version
     << "\nEND\n";
}

void io::Out_Configuration::init(io::Argument & args,
        simulation::Parameter const & param) {
  if (args.count(argname_fin) > 0) {
    if (!param.analyze.analyze) final_configuration(args[argname_fin]);
  } else io::messages.add("argument " + argname_fin + " for final configuration required!",
          "Out_Configuration",
          io::message::error);

  if (args.count(argname_trj) > 0)
    trajectory(args[argname_trj], param.write.position);
  else if (param.write.position)
    io::messages.add("write trajectory but no " + argname_trj + " argument",
          "Out_Configuration",
          io::message::error);

  if (args.count(argname_trv) > 0)
    velocity_trajectory(args[argname_trv], param.write.velocity);
  else if (param.write.velocity)
    io::messages.add("write velocity trajectory but no trv argument",
          "Out_Configuration",
          io::message::error);

  if (args.count(argname_trf) > 0)
    force_trajectory(args[argname_trf], param.write.force);
  else if (param.write.force)
    io::messages.add("write force trajectory but no trf argument",
          "Out_Configuration",
          io::message::error);

  m_write_special = param.polarise.write || param.jvalue.write || param.xrayrest.write 
     || param.distanceres.write || param.distancefield.write || param.angrest.write
     || param.dihrest.write || param.print.monitor_dihedrals
     || param.localelev.write || param.electric.dip_write || param.electric.cur_write 
     || param.addecouple.write || param.nemd.write || param.orderparamrest.write || param.rdc.write
     || param.bsleus.write;    // add others if there are any

  if (args.count(argname_trs) > 0)
    special_trajectory(args[argname_trs], param.polarise.write, 
            param.jvalue.write, param.xrayrest.write, param.distanceres.write, 
            param.distancefield.write, param.angrest.write, param.dihrest.write,
            param.print.monitor_dihedrals,param.localelev.write, 
            param.electric.dip_write, param.electric.cur_write, param.addecouple.write,
            param.nemd.write, param.orderparamrest.write, param.rdc.write,
            param.bsleus.write);
  else if (m_write_special)
    io::messages.add("write special trajectory but no trs argument",
          "Out_Configuration",
          io::message::error);

  if (args.count(argname_tre) > 0)
    energy_trajectory(args[argname_tre], param.write.energy);
  else if (param.write.energy)
    io::messages.add("write energy trajectory but no " + argname_tre + " argument",
          "Out_Configuration",
          io::message::error);

  if (args.count(argname_trg) > 0)
    free_energy_trajectory(args[argname_trg], param.write.free_energy);
  else if (param.write.free_energy)
    io::messages.add("write free energy trajectory but no trg argument",
          "Out_Configuration",
          io::message::error);

  if (args.count(argname_bae) > 0)
    block_averaged_energy(args[argname_bae], param.write.block_average);
  else if (param.write.block_average && param.write.energy)
    io::messages.add("write block averaged energy but no bae argument",
          "Out_Configuration",
          io::message::error);

  if (param.perturbation.perturbation) {
    if (args.count(argname_bag) > 0)
      block_averaged_free_energy(args[argname_bag],
            param.write.block_average);
    else if (param.write.block_average && param.write.free_energy)
      io::messages.add("write block averaged free energy "
            "but no bag argument",
            "Out_Configuration",
            io::message::error);
  }

  if (param.replica.num_T && param.replica.num_l) {
    m_replica = true;
  }

}

void io::Out_Configuration::write(configuration::Configuration &conf,
        topology::Topology const &topo,
        simulation::Simulation const &sim,
        output_format const form) {
  // standard trajectories

  bool constraint_force = sim.param().constraint.solute.algorithm == simulation::constr_shake ||
          sim.param().constraint.solvent.algorithm == simulation::constr_shake;

  // check whether a new energy minimum was found
  bool minimum_found = false;
  if (sim.param().write.energy_index > 0) {
    double current_energy = conf.old().energies.get_energy_by_index(sim.param().write.energy_index);

    // found a new minimum?
    if (current_energy < minimum_energy) {
      minimum_found = true;
      minimum_energy = current_energy;
    }
  }

  if (form == reduced) {
    /**
     * set this to true when you print the timestep to the special traj.
     * make sure you don't print it twice. 
     */

    if (m_every_pos && ((sim.steps() % m_every_pos) == 0 || minimum_found)) {
        _print_timestep(sim, m_pos_traj);

        if (sim.param().write.position_solute_only)
          _print_positionred(conf, topo, topo.num_solute_atoms(), m_pos_traj);
        else
          _print_positionred(conf, topo, topo.num_atoms(), m_pos_traj);

        if (conf.boundary_type != math::vacuum)
          _print_box(conf, m_pos_traj);

        m_pos_traj.flush();
      // a new block begins. let's reset the minimum
      minimum_energy = conf.old().energies.get_energy_by_index(sim.param().write.energy_index);
    }

    if (m_every_vel && (sim.steps() % m_every_vel) == 0) {
      _print_timestep(sim, m_vel_traj);
      if (sim.param().write.velocity_solute_only)
        _print_velocityred(conf, topo.num_solute_atoms(), m_vel_traj);
      else
        _print_velocityred(conf, topo.num_atoms(), m_vel_traj);
      m_vel_traj.flush();
    }

    if (m_every_force && ((sim.steps()-sim.param().analyze.stride) % m_every_force) == 0) {
      if (sim.steps()) {
        _print_old_timestep(sim, m_force_traj);
        if (sim.param().write.force_solute_only)
          _print_forcered(conf, topo.num_solute_atoms(), m_force_traj, constraint_force);
        else
          _print_forcered(conf, topo.num_atoms(), m_force_traj, constraint_force);
        m_force_traj.flush();
      }
    }
    
    // special trajectory blocks
    // IMPORTANT: because this writing is done before the md_sequence is applied
    // most special terms for the current configuration are written at the 
    // beginning of the next step; they therefore have to be printed before 
    // _print_special_timestep is executed -- MariaP 2016/04/29

    if (m_every_disres && sim.steps() && ((sim.steps()-sim.param().analyze.stride) % m_every_disres) == 0) {
      _print_distance_restraints(conf, topo, m_special_traj);
      m_special_traj.flush();
    }

    if (m_every_disfieldres && sim.steps() && ((sim.steps()-sim.param().analyze.stride) % m_every_disfieldres) == 0) {
      _print_disfield_restraints(conf, topo, m_special_traj);
      m_special_traj.flush();
    }

    if (m_every_angres && sim.steps() && ((sim.steps()-sim.param().analyze.stride) % m_every_angres) == 0) {
      _print_angle_restraints(conf, topo, m_special_traj);
      m_special_traj.flush();
    }

    if (m_every_dihres && sim.steps() && ((sim.steps()-sim.param().analyze.stride) % m_every_dihres) == 0) {
      _print_dihedral_restraints(conf, topo, m_special_traj);
      m_special_traj.flush();
    }

    if (m_every_jvalue  && sim.steps() && ((sim.steps()-sim.param().analyze.stride) % m_every_jvalue) == 0) {
      _print_jvalue(sim.param(), conf, topo, m_special_traj, true);
      m_special_traj.flush();
    }

    if (m_every_xray && sim.steps() && ((sim.steps()-1) % m_every_xray) == 0) {
      _print_xray_rvalue(sim.param(), conf, m_special_traj);
      _print_xray_umbrellaweightthresholds(sim.param(), topo, m_special_traj);
      _print_xray_bfactors(sim.param(), conf, m_special_traj);
      m_special_traj.flush();
    }

    if (m_every_oparam && sim.steps() && ((sim.steps()-1) % m_every_oparam) == 0) {
      _print_order_parameter_restraints(conf, topo, m_special_traj);
      m_special_traj.flush();
    }

    if (m_every_rdc && sim.steps() && ((sim.steps()-1) % m_every_rdc) == 0) {
      if(sim.param().rdc.mode != simulation::rdc_restr_off) {
        io::Out_Configuration::_print_rdc_values(sim.param(), conf, topo, m_special_traj, /*formatted*/ true);
      }
      if(sim.param().rdc.mode == simulation::rdc_restr_av ||
         sim.param().rdc.mode == simulation::rdc_restr_av_weighted ||
         sim.param().rdc.mode == simulation::rdc_restr_biq ||
         sim.param().rdc.mode == simulation::rdc_restr_biq_weighted) {
        io::Out_Configuration::_print_rdc_averages(sim.param(), conf, topo, m_special_traj, /*formatted*/ true);
      }
      if(sim.param().rdc.mode != simulation::rdc_restr_off) {
        io::Out_Configuration::_print_rdc_representation(sim.param(), conf, topo, m_special_traj, /*formatted*/ true);
      }
      if (sim.param().rdc.method == simulation::rdc_sd) {
        io::Out_Configuration::_print_rdc_stochastic_integrals(sim.param(), conf, topo, m_special_traj, /*formatted*/ true);
      }
      m_special_traj.flush();
    }

    if (m_every_adde && sim.steps() && ((sim.steps()-1) % m_every_adde) == 0) {
      _print_adde(sim, topo, conf, m_special_traj);
      m_special_traj.flush();
    }
    
    if (m_every_nemd && sim.steps() && ((sim.steps()-1) % m_every_nemd) == 0) {
      _print_nemd(sim, topo, conf, m_special_traj);
      m_special_traj.flush();
    }
    
    if (m_every_leus && sim.steps() && ((sim.steps()-sim.param().analyze.stride) % m_every_leus) == 0) {
      _print_umbrellas(conf, m_special_traj);
      m_special_traj.flush();
    }
    
    if ((sim.param().bsleus.bsleus == simulation::bsleus_on) && 
            m_every_bsleus && sim.steps() && ((sim.steps()-sim.param().analyze.stride) % m_every_bsleus) == 0) {
      _print_bsleusmem(conf, m_special_traj);
      _print_bsleus(conf, m_special_traj);
      m_special_traj.flush();
    }  

    if (m_every_cos_pos && sim.steps() && ((sim.steps()-sim.param().analyze.stride) % m_every_cos_pos) == 0) {
      _print_cos_position(conf, topo, m_special_traj);
      m_special_traj.flush();
    }
    
    // print timestep to the special trajectory if any special terms 
    // are to be written for the current timestep (even if most printing is only 
    // done at the start of the next timestep
    _print_special_timestep(sim);  

    if (m_every_dat) {
      _print_dihangle_trans(conf, topo, m_special_traj);
      m_special_traj.flush();
    }
    
    // WARNING: with polarization it seems like _print_dipole is using 
    // positions and cos-displacement that don't belong to the same configuration
    // this should be looked into -- MariaP 16/05/12
    if (m_every_dipole && (sim.steps() % m_every_dipole) == 0) {
      SPLIT_BOUNDARY(_print_dipole, sim, topo, conf, m_special_traj);
      m_special_traj.flush();
    }

    if (m_every_current && (sim.steps() % m_every_current) == 0) {
      _print_current(sim, topo, conf, m_special_traj);
      m_special_traj.flush();
    }    
    // end of special trajectory printing statements
    
    if (m_every_energy && (((sim.steps()-sim.param().analyze.stride) % m_every_energy) == 0 || minimum_found)) {
      if (sim.steps()) {
        _print_old_timestep(sim, m_energy_traj);
        _print_energyred(conf, m_energy_traj);
        _print_volumepressurered(topo, conf, sim, m_energy_traj);
        m_energy_traj.flush();
      }
    }

    if (m_every_free_energy && ((sim.steps()-sim.param().analyze.stride) % m_every_free_energy) == 0) {
      if (sim.steps()) {
        _print_old_timestep(sim, m_free_energy_traj);
        _print_free_energyred(conf, topo, m_free_energy_traj);
        m_free_energy_traj.flush();
      }
    }

    if (m_every_blockaverage && ((sim.steps()-sim.param().analyze.stride) % m_every_blockaverage) == 0) {

      if (m_write_blockaverage_energy) {
        if (sim.steps()) {
          _print_old_timestep(sim, m_blockaveraged_energy);
          _print_blockaveraged_energyred(conf, m_blockaveraged_energy);
          _print_blockaveraged_volumepressurered(conf, sim, m_blockaveraged_energy);
          m_blockaveraged_energy.flush();
        }
      }

      if (m_write_blockaverage_free_energy) {
        if (sim.steps()) {
          _print_old_timestep(sim, m_blockaveraged_free_energy);
          _print_blockaveraged_free_energyred(conf, sim.param().perturbation.dlamt,
                  m_blockaveraged_free_energy);
          m_blockaveraged_free_energy.flush();
        }
      }
      conf.current().averages.block().zero();
    }

  } else if (form == final) {
    if (m_final) {
    _print_timestep(sim, m_final_conf);
    _print_position(conf, topo, m_final_conf);
    _print_lattice_shifts(conf, topo, m_final_conf);

    if (sim.param().polarise.cos)
      _print_cos_position(conf, topo, m_final_conf);

    if (sim.param().minimise.ntem == 0)
      _print_velocity(conf, topo, m_final_conf);

    _print_box(conf, m_final_conf);

    if (sim.param().constraint.solute.algorithm
            == simulation::constr_flexshake) {
      _print_flexv(conf, topo, m_final_conf);
    }

    if (sim.param().stochastic.sd) {
      _print_stochastic_integral(conf, topo, m_final_conf);
    }

    if (sim.param().perturbation.perturbation) {
      _print_pertdata(topo, m_final_conf);
    }

    if (sim.param().distanceres.distanceres < 0) {
      _print_distance_restraint_averages(conf, topo, m_final_conf);
    }

    if(sim.param().distancefield.printgrid) {
      _print_disfield_grid(sim.param(), conf, topo, m_final_conf);
    }
    
    if (sim.param().posrest.posrest != simulation::posrest_off) {
      _print_position_restraints(sim, topo, conf, m_final_conf);
    }

    if (sim.param().jvalue.mode != simulation::jvalue_restr_off) {
      _print_jvalue(sim.param(), conf, topo, m_final_conf, false);
    }

    if (sim.param().xrayrest.xrayrest != simulation::xrayrest_off) {
      _print_xray(sim.param(), conf, topo, m_final_conf, /*final=*/ true);
      _print_xray_umbrellaweightthresholds(sim.param(), topo, m_final_conf);
      _print_xray_bfactors(sim.param(), conf, m_final_conf);
    }

    if (sim.param().localelev.localelev != simulation::localelev_off) {
      _print_umbrellas(conf, m_final_conf);
    }
    
    if (sim.param().bsleus.bsleus != simulation::bsleus_off){
      _print_bsleusmem(conf, m_final_conf);
      _print_bsleuspos(conf, m_final_conf);
    }

    if (sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_av ||
        sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_av_weighted) {
      _print_order_parameter_restraint_averages(conf, topo, m_final_conf);
    } else if (sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_winav ||
        sim.param().orderparamrest.orderparamrest == simulation::oparam_restr_winav_weighted) {
      _print_order_parameter_restraint_average_window(conf, topo, m_final_conf);
    }

    if(sim.param().rdc.mode == simulation::rdc_restr_av ||
       sim.param().rdc.mode == simulation::rdc_restr_av_weighted ||
       sim.param().rdc.mode == simulation::rdc_restr_biq ||
       sim.param().rdc.mode == simulation::rdc_restr_biq_weighted) {
      io::Out_Configuration::_print_rdc_averages(sim.param(), conf, topo, m_final_conf, /*formatted*/ false);
    }
  
    if(sim.param().rdc.mode != simulation::rdc_restr_off) {
      io::Out_Configuration::_print_rdc_representation(sim.param(), conf, topo, m_final_conf, /*formatted*/ false);
    }
  
    if (sim.param().rdc.method == simulation::rdc_sd) {
      io::Out_Configuration::_print_rdc_stochastic_integrals(sim.param(), conf, topo, m_final_conf, /*formatted*/ false);
    }

    if (sim.param().rottrans.rottrans) {
      _print_rottrans(conf, sim, m_final_conf);
    }

    if (sim.param().pscale.jrest) {
      _print_pscale_jrest(conf, topo, m_final_conf);
    }

    if (sim.param().multibath.algorithm > 1) {
      _print_nose_hoover_chain_variables(sim.multibath(), m_final_conf);
    }
    if (conf.special().shake_failure_occurred) {
      _print_shake_failure(conf, topo, m_final_conf);
    }
    }

    if (sim.param().eds.form == simulation::aeds_search_eir || sim.param().eds.form == simulation::aeds_search_emax_emin || sim.param().eds.form == simulation::aeds_search_all) {
      _print_aedssearch(conf, sim, m_final_conf);
    }

    // forces and energies still go to their trajectories
    if (m_every_force && ((sim.steps()-sim.param().analyze.stride) % m_every_force) == 0) {
      _print_old_timestep(sim, m_force_traj);
      if (sim.param().write.force_solute_only)
        _print_forcered(conf, topo.num_solute_atoms(), m_force_traj, constraint_force);
      else
        _print_forcered(conf, topo.num_atoms(), m_force_traj, constraint_force);
    }

    if (m_every_energy && ((sim.steps()-sim.param().analyze.stride) % m_every_energy) == 0) {
      _print_old_timestep(sim, m_energy_traj);
      _print_energyred(conf, m_energy_traj);
      _print_volumepressurered(topo, conf, sim, m_energy_traj);
    }

    if (m_every_free_energy && ((sim.steps()-sim.param().analyze.stride) % m_every_free_energy) == 0) {
      _print_old_timestep(sim, m_free_energy_traj);
      _print_free_energyred(conf, topo, m_free_energy_traj);
    }
    if (m_every_blockaverage && ((sim.steps()-sim.param().analyze.stride) % m_every_blockaverage) == 0) {

      if (m_write_blockaverage_energy) {
        if (sim.steps()) {
          _print_old_timestep(sim, m_blockaveraged_energy);
          _print_blockaveraged_energyred(conf, m_blockaveraged_energy);
          _print_blockaveraged_volumepressurered(conf, sim, m_blockaveraged_energy);
        }
      }

      if (m_write_blockaverage_free_energy) {
        if (sim.steps()) {
          _print_old_timestep(sim, m_blockaveraged_free_energy);
          _print_blockaveraged_free_energyred(conf, sim.param().perturbation.dlamt,
                  m_blockaveraged_free_energy);
        }
      }
      conf.current().averages.block().zero();
    }
   
    // the last step of special terms which are printed after the md sequence
    // still has to be written to the special trajectory
    // the timestep has already been printed
 
    if (m_every_disres && ((sim.steps()-sim.param().analyze.stride) % m_every_disres) == 0) {   
      _print_distance_restraints(conf, topo, m_special_traj);
    }
    
    if (m_every_disfieldres && ((sim.steps()-sim.param().analyze.stride) % m_every_disfieldres) == 0) {    
      _print_disfield_restraints(conf, topo, m_special_traj);
    }
    
    if (m_every_angres && ((sim.steps()-sim.param().analyze.stride) % m_every_angres) == 0) {   
      _print_angle_restraints(conf, topo, m_special_traj);
    }
    
    if (m_every_dihres && ((sim.steps()-sim.param().analyze.stride) % m_every_dihres) == 0) {   
      _print_dihedral_restraints(conf, topo, m_special_traj);
    }
    
    if (m_every_jvalue && ((sim.steps()-sim.param().analyze.stride) % m_every_jvalue) == 0) {    
      _print_jvalue(sim.param(), conf, topo, m_special_traj, true);
    }
    
    if (m_every_xray && ((sim.steps()-sim.param().analyze.stride) % m_every_xray) == 0) {    
      _print_xray_rvalue(sim.param(), conf, m_special_traj);
      _print_xray_umbrellaweightthresholds(sim.param(), topo, m_special_traj);
      _print_xray_bfactors(sim.param(), conf, m_special_traj);
    }

    if (m_every_oparam &&  ((sim.steps()-1) % m_every_oparam) == 0) {
      _print_order_parameter_restraints(conf, topo, m_special_traj);
    }

    if (m_every_rdc && ((sim.steps()-1) % m_every_rdc) == 0) {
      if(sim.param().rdc.mode != simulation::rdc_restr_off) {
        io::Out_Configuration::_print_rdc_values(sim.param(), conf, topo, m_special_traj, /*formatted*/ true);
      }
      if(sim.param().rdc.mode == simulation::rdc_restr_av ||
         sim.param().rdc.mode == simulation::rdc_restr_av_weighted ||
         sim.param().rdc.mode == simulation::rdc_restr_biq ||
         sim.param().rdc.mode == simulation::rdc_restr_biq_weighted) {
        io::Out_Configuration::_print_rdc_averages(sim.param(), conf, topo, m_special_traj, /*formatted*/ true);
      }
      if(sim.param().rdc.mode != simulation::rdc_restr_off) {
        io::Out_Configuration::_print_rdc_representation(sim.param(), conf, topo, m_special_traj, /*formatted*/ true);
      }
      if (sim.param().rdc.method == simulation::rdc_sd) {
        io::Out_Configuration::_print_rdc_stochastic_integrals(sim.param(), conf, topo, m_special_traj, /*formatted*/ true);
      }
    }

    if (m_every_adde  && ((sim.steps()-1) % m_every_adde) == 0) {
      _print_adde(sim, topo, conf, m_special_traj);
    }
    
    if (m_every_nemd && ((sim.steps()-1) % m_every_nemd) == 0) {
      _print_nemd(sim, topo, conf, m_special_traj);
    }
    
    if (m_every_leus && ((sim.steps()-sim.param().analyze.stride) % m_every_leus) == 0) {
      _print_umbrellas(conf, m_special_traj);
    }

    if ((sim.param().bsleus.bsleus == simulation::bsleus_on) &&
              m_every_bsleus && ((sim.steps()-sim.param().analyze.stride) % m_every_bsleus) == 0) {
        _print_bsleusmem(conf, m_special_traj);
        _print_bsleus(conf, m_special_traj);
    }
    
    if (m_every_cos_pos && ((sim.steps()-sim.param().analyze.stride) % m_every_cos_pos) == 0) {
      _print_cos_position(conf, topo, m_special_traj);
    }
    
    
    //end of special traj printing
  } else {

    // not reduced or final (so: decorated)

    if (m_every_pos && (sim.steps() % m_every_pos) == 0) {
      _print_timestep(sim, m_pos_traj);
      _print_position(conf, topo, m_pos_traj);
      if (conf.boundary_type != math::vacuum)
        _print_box(conf, m_pos_traj);
    }

    if (m_every_vel && (sim.steps() % m_every_vel) == 0) {
      _print_timestep(sim, m_vel_traj);
      _print_velocity(conf, topo, m_vel_traj);
    }

    if (m_every_force && (sim.steps() % m_every_force) == 0) {
      if (sim.steps()) {
        _print_timestep(sim, m_force_traj);
        _print_force(conf, topo, m_force_traj, constraint_force);
      }
    }
  }

  // done writing!

}

void io::Out_Configuration
::final_configuration(std::string name) {
  m_final_conf.open(name.c_str());

  _print_title(m_title, "final configuration", m_final_conf);
  m_final = true;
}

void io::Out_Configuration
::trajectory(std::string name, int every) {
  m_pos_traj.open(name.c_str());

  m_every_pos = every;
  _print_title(m_title, "position trajectory", m_pos_traj);
}

void io::Out_Configuration
::velocity_trajectory(std::string name, int every) {
  m_vel_traj.open(name.c_str());

  m_every_vel = every;
  _print_title(m_title, "velocity trajectory", m_vel_traj);
}

void io::Out_Configuration
::force_trajectory(std::string name, int every) {
  m_force_traj.open(name.c_str());

  m_every_force = every;
  _print_title(m_title, "force trajectory", m_force_traj);
}

void io::Out_Configuration
::special_trajectory(std::string name, int every_cos, int every_jvalue, 
                     int every_xray, int every_disres, int every_disfieldres, 
                     int every_angres, int every_dihres, int every_dat, 
                     int every_leus, int every_dipole, int every_current,
                     int every_adde, int every_nemd, int every_oparam, int every_rdc,
                     int every_bsleus) {

  m_special_traj.open(name.c_str());

  m_every_cos_pos = every_cos;
  m_every_jvalue = every_jvalue;
  m_every_xray = every_xray;
  m_every_disres = every_disres;
  m_every_disfieldres = every_disfieldres;
  m_every_angres = every_angres;
  m_every_dihres = every_dihres;
  m_every_dat = every_dat;
  m_every_leus = every_leus;
  m_every_dipole = every_dipole;
  m_every_current = every_current;
  m_every_adde = every_adde;
  m_every_nemd = every_nemd;
  m_every_oparam = every_oparam;
  m_every_rdc = every_rdc; 
  m_every_bsleus = every_bsleus;
  _print_title(m_title, "special trajectory", m_special_traj);
}

void io::Out_Configuration
::energy_trajectory(std::string name, int every) {
  m_energy_traj.open(name.c_str());

  m_every_energy = every;
  _print_title(m_title, "energy trajectory", m_energy_traj);
  _print_ene_version(m_energy_traj);
}

void io::Out_Configuration
::free_energy_trajectory(std::string name, int every) {
  m_free_energy_traj.open(name.c_str());

  m_every_free_energy = every;
  _print_title(m_title, "free energy trajectory", m_free_energy_traj);
  _print_ene_version(m_free_energy_traj);
}

void io::Out_Configuration
::block_averaged_energy(std::string name, int every) {
  m_blockaveraged_energy.open(name.c_str());

  if (m_every_blockaverage && m_every_blockaverage != every) {
    io::messages.add("overwriting how often block averages are written out illegal",
            "Out_Configuration",
            io::message::error);
  }
  m_every_blockaverage = every;
  m_write_blockaverage_energy = true;
  _print_title(m_title, "block averaged energies", m_blockaveraged_energy);
}

void io::Out_Configuration
::block_averaged_free_energy(std::string name, int every) {
  m_blockaveraged_free_energy.open(name.c_str());

  if (m_every_blockaverage && m_every_blockaverage != every) {
    io::messages.add("overwriting how often block averages are written out illegal",
            "Out_Configuration",
            io::message::error);
  }
  m_every_blockaverage = every;
  m_write_blockaverage_free_energy = true;
  _print_title(m_title, "block averaged free energies", m_blockaveraged_free_energy);
}

void io::Out_Configuration
::_print_timestep(simulation::Simulation const &sim,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "TIMESTEP\n"
          << std::setw(m_width) << sim.steps()
          << " "
          << std::setw(m_width - 1) << sim.time()
          << "\nEND\n";
}

void io::Out_Configuration
::_print_old_timestep(simulation::Simulation const &sim,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "TIMESTEP\n"
          << std::setw(m_width) << sim.steps()-sim.param().analyze.stride
          << " "
          << std::setw(m_width - 1) << sim.time() - sim.time_step_size()*sim.param().analyze.stride
          << "\nEND\n";

}

void io::Out_Configuration
::_print_special_timestep(simulation::Simulation const &sim) {
   if (m_write_special  && ( (m_every_cos_pos && ((sim.steps() % m_every_cos_pos) == 0)) ||
                             (m_every_jvalue && ((sim.steps() % m_every_jvalue) == 0)) ||
                             (m_every_xray && ((sim.steps() % m_every_xray) == 0)) ||
                             (m_every_disres && ((sim.steps() % m_every_disres) == 0)) ||
                             (m_every_disfieldres && ((sim.steps() % m_every_disfieldres) == 0)) ||
                             (m_every_angres && ((sim.steps() % m_every_angres) == 0)) ||
                             (m_every_dihres && ((sim.steps() % m_every_dihres) == 0)) ||
                             (m_every_dat && ((sim.steps() % m_every_dat) == 0)) ||
                             (m_every_leus && ((sim.steps() % m_every_leus) == 0)) ||
                             (m_every_bsleus && ((sim.steps() % m_every_bsleus) == 0)) ||
                             (m_every_dipole && ((sim.steps() % m_every_dipole) == 0)) ||
                             (m_every_current && ((sim.steps() % m_every_current) == 0)) ||
                             (m_every_adde && ((sim.steps() % m_every_adde) == 0)) ||
                             (m_every_nemd && ((sim.steps() % m_every_nemd) == 0)) ||
                             (m_every_oparam && ((sim.steps() % m_every_oparam) == 0)) ||
                             (m_every_rdc && ((sim.steps() % m_every_rdc) == 0)) ))
        _print_timestep(sim, m_special_traj);
}

template<math::boundary_enum b>
void _print_g96_position_bound(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os, int width,
        bool old_conf = false) {
  const configuration::Configuration::state_struct * state = nullptr;
  if (old_conf)
    state = & conf.old();
  else
    state = & conf.current();

  math::Periodicity<b> periodicity(state->box);
  math::VArray const &pos = state->pos;
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();


  math::Vec v, v_box, trans, r;
  //matrix to rotate back into orignial Cartesian Coordinat system
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));
  if (conf.boundary_type == math::truncoct) {
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));
  }

  os << "# first 24 chars ignored\n";

  // put chargegroups into the box (on the fly)
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
          cg_to = topo.chargegroup_end();

  // solute chargegroups...
  unsigned int i = 0;
  for (; i < topo.num_solute_chargegroups(); ++cg_it, ++i) {
    // gather on first atom...
    v = pos(*cg_it.begin());
    v_box = v;
    periodicity.put_into_positive_box(v_box);
    trans = v_box - v;

    // atoms in a chargegroup
    topology::Atom_Iterator at_it = cg_it.begin(),
            at_to = cg_it.end();

    for (; at_it != at_to; ++at_it) {
      r = pos(*at_it) + trans;
      //rotate to original Cartesian coordinates
      r = math::Vec(math::product(Rmat, r));
      os << std::setw(5) << solute.atom(*at_it).residue_nr + 1 << " "
              << std::setw(5) << std::left
              << residue_name[solute.atom(*at_it).residue_nr] << " "
              << std::setw(6) << std::left << solute.atom(*at_it).name
              << std::right
              << std::setw(6) << *at_it + 1
              << std::setw(width) << r(0)
              << std::setw(width) << r(1)
              << std::setw(width) << r(2)
              << "\n";
    }
  }

  // solvent chargegroups
  unsigned int s = 0;
  unsigned int mol = 0;

  for (; cg_it != cg_to; ++cg_it, ++mol) {
    v = pos(**cg_it);
    v_box = v;
    periodicity.put_into_positive_box(v_box);
    trans = v_box - v;

    if (mol >= topo.num_solvent_molecules(s))++s;

    // loop over the atoms
    topology::Atom_Iterator at_it = cg_it.begin(),
            at_to = cg_it.end();
    // one chargegroup per solvent
    for (unsigned int atom = 0; at_it != at_to; ++at_it, ++atom) {
      r = pos(*at_it) + trans;
      //rotate to original Cartesian coordinates
      r = math::Vec(math::product(Rmat, r));
      os << std::setw(5) << mol + 1
              << ' ' << std::setw(5) << std::left
              << residue_name[topo.solvent(s).atom(atom).residue_nr] << " "
              << std::setw(6) << std::left << topo.solvent(s).atom(atom).name
              << std::right
              << std::setw(6) << *at_it + 1
              << std::setw(width) << r(0)
              << std::setw(width) << r(1)
              << std::setw(width) << r(2)
              << "\n";
    }
  }

}

/**
 * i need a specialized function to put the particles into the box.
 */
template<math::boundary_enum b>
void _print_position_bound(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os, int width) {
  math::Periodicity<b> periodicity(conf.current().box);
  math::VArray const &pos = conf.current().pos;
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();

  math::Vec v;
  //matrix to rotate back into orignial Cartesian Coordinat system
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));
  if (conf.boundary_type == math::truncoct) {
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));
  }

  os << "# first 24 chars ignored\n";

  for (int i = 0, to = topo.num_solute_atoms(); i < to; ++i) {

    v = pos(i);
    periodicity.put_into_box(v);
    //rotate to original Cartesian coordinates
    v = math::Vec(math::product(Rmat, v));
    os << std::setw(5) << solute.atom(i).residue_nr + 1 << " "
            << std::setw(5) << std::left
            << residue_name[solute.atom(i).residue_nr] << " "
            << std::setw(6) << std::left << solute.atom(i).name << std::right
            << std::setw(6) << i + 1
            << std::setw(width) << v(0)
            << std::setw(width) << v(1)
            << std::setw(width) << v(2)
            << "\n";
  }

  int index = topo.num_solute_atoms();
  int res_nr = 1;

  for (unsigned int s = 0; s < topo.num_solvents(); ++s) {

    for (unsigned int m = 0; m < topo.num_solvent_molecules(s); ++m, ++res_nr) {

      for (unsigned int a = 0; a < topo.solvent(s).num_atoms(); ++a, ++index) {

        v = pos(index);
        periodicity.put_into_positive_box(v);
        //rotate to original Cartesian coordinates
        v = math::Vec(math::product(Rmat, v));
        os << std::setw(5) << res_nr
                << ' ' << std::setw(5) << std::left
                << residue_name[topo.solvent(s).atom(a).residue_nr] << " "
                << std::setw(6) << std::left
                << topo.solvent(s).atom(a).name << std::right
                << std::setw(6) << index + 1
                << std::setw(width) << v(0)
                << std::setw(width) << v(1)
                << std::setw(width) << v(2)
                << "\n";
      }
    }
  }
}

void io::Out_Configuration
::_print_position(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "POSITION\n";

  SPLIT_BOUNDARY(_print_g96_position_bound,
          conf, topo, os, m_width);

  os << "END\n";

}

void io::Out_Configuration
::_print_shake_failure(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "SHAKEFAILPOSITION\n";
  SPLIT_BOUNDARY(_print_g96_position_bound,
          conf, topo, os, m_width);
  os << "END\n";
  os << "SHAKEFAILPREVPOSITION\n";
  SPLIT_BOUNDARY(_print_g96_position_bound,
          conf, topo, os, m_width, true /* old coords */);
  os << "END\n";

}

template<math::boundary_enum b>
void _print_g96_positionred_bound(configuration::Configuration const &conf,
        topology::Topology const &topo,
        int num,
        std::ostream &os, int width) {
  DEBUG(10, "g96 positionred");

  math::Periodicity<b> periodicity(conf.current().box);
  math::VArray pos = conf.current().pos;

  math::Vec v, v_box, trans, r;
  //matrix to rotate back into original Cartesian Coordinate system
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));

  if (conf.boundary_type == math::truncoct)
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));

  assert(num >= 0);

  // put chargegroups into the box (on the fly)
  topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
          cg_to = topo.chargegroup_end();
  DEBUG(10, "cg to : " << **cg_to << std::endl);

  // solute chargegroups...
  unsigned int i = 0, count = 0;
  for (; i < topo.num_solute_chargegroups(); ++cg_it, ++i) {
    DEBUG(10, "solute cg: " << i);
    // gather on first atom...
    v = pos(*cg_it.begin());
    v_box = v;
    periodicity.put_into_positive_box(v_box);
    trans = v_box - v;

    // atoms in a chargegroup
    topology::Atom_Iterator at_it = cg_it.begin(),
            at_to = cg_it.end();

    for (; at_it != at_to; ++at_it, ++count) {

      if (*at_it >= unsigned(num)) return;

      DEBUG(10, "atom: " << count);
      r = pos(*at_it) + trans;
      //rotate back to original Carthesian Coord
      r = math::Vec(math::product(Rmat, r));
      os << std::setw(width) << r(0)
              << std::setw(width) << r(1)
              << std::setw(width) << r(2)
              << "\n";

      if ((count + 1) % 10 == 0) os << '#' << std::setw(10) << count + 1 << "\n";

    }
  }

  DEBUG(10, "solvent");

  // solvent chargegroups
  unsigned int mol = 0;

  for (; cg_it != cg_to; ++cg_it, ++mol) {
    DEBUG(10, "solvent " << mol);

    v = pos(**cg_it);
    v_box = v;
    periodicity.put_into_positive_box(v_box);
    trans = v_box - v;

    // loop over the atoms
    topology::Atom_Iterator at_it = cg_it.begin(),
            at_to = cg_it.end();
    // one chargegroup per solvent
    for (; at_it != at_to; ++at_it, ++count) {
      DEBUG(10, "\tatom " << count);

      if (*at_it >= unsigned(num)) return;

      r = pos(*at_it) + trans;
      //rotate back to original Carthesian Coord
      r = math::Vec(math::product(Rmat, r));
      os << std::setw(width) << r(0)
              << std::setw(width) << r(1)
              << std::setw(width) << r(2)
              << "\n";

      if ((count + 1) % 10 == 0) os << '#' << std::setw(10) << count + 1 << "\n";

    }
  }
}

template<math::boundary_enum b>
void
_print_positionred_bound(configuration::Configuration const &conf,
        int num,
        std::ostream &os, int width) {
  math::Periodicity<b> periodicity(conf.current().box);

  math::VArray const &pos = conf.current().pos;
  math::Vec v;
  //matrix to rotate back into original Cartesian Coordinate system
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));

  if (conf.boundary_type == math::truncoct)
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));

  DEBUG(10, "writing POSITIONRED " << pos.size());

  for (int i = 0; i < num; ++i) {

    v = pos(i);
    periodicity.put_into_box(v);
    //rotate back into original Cartesian coordinates
    v = math::Vec(math::product(Rmat, v));
    os << std::setw(width) << v(0)
            << std::setw(width) << v(1)
            << std::setw(width) << v(2)
            << "\n";

    if ((i + 1) % 10 == 0) os << '#' << std::setw(10) << i + 1 << "\n";
  }

}

inline void io::Out_Configuration
::_print_positionred(configuration::Configuration const &conf,
        topology::Topology const &topo,
        int num,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "POSITIONRED\n";
  DEBUG(7, "configuration boundary type :" << conf.boundary_type);

  SPLIT_BOUNDARY(_print_g96_positionred_bound, conf, topo, num, os, m_width);

  os << "END\n";

}

inline void io::Out_Configuration
::_print_cos_position(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "COSDISPLACEMENTS\n";

  math::VArray const &posV = conf.current().posV;
  //rotate back to original carthesian coordinates
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));
  if (conf.boundary_type == math::truncoct) {
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));
  }
  math::Vec posV_rot;
  for (unsigned int i = 0; i < posV.size(); ++i) {
    posV_rot = math::Vec(math::product(Rmat, posV(i)));
    os << std::setw(m_width) << posV_rot(0)
            << std::setw(m_width) << posV_rot(1)
            << std::setw(m_width) << posV_rot(2)
            << "\n";

    if ((i + 1) % 10 == 0) os << '#' << std::setw(10) << i + 1 << "\n";
  }

  os << "END\n";
}

void io::Out_Configuration
::_print_velocity(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "VELOCITY\n";

  math::VArray const &vel = conf.current().vel;
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();

  os << "# first 24 chars ignored\n";
  //matrix to rotate back into orignial Cartesian Coordinat system
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));
  if (conf.boundary_type == math::truncoct)
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));

  math::Vec vel_rot;
  for (int i = 0, to = topo.num_solute_atoms(); i < to; ++i) {
    //rotate back to original Carthesian coordinates
    vel_rot = math::Vec(math::product(Rmat, vel(i)));
    os << std::setw(5) << solute.atom(i).residue_nr + 1 << " "
            << std::setw(5) << std::left << residue_name[solute.atom(i).residue_nr] << " "
            << std::setw(6) << std::left << solute.atom(i).name << std::right
            << std::setw(6) << i + 1
            << std::setw(m_width) << vel_rot(0)
            << std::setw(m_width) << vel_rot(1)
            << std::setw(m_width) << vel_rot(2)
            << "\n";
  }

  int index = topo.num_solute_atoms();
  int res_num = 1;

  for (unsigned int s = 0; s < topo.num_solvents(); ++s) {

    for (unsigned int m = 0; m < topo.num_solvent_molecules(s); ++m, ++res_num) {

      for (unsigned int a = 0; a < topo.solvent(s).num_atoms(); ++a, ++index) {
        vel_rot = math::Vec(math::product(Rmat, vel(index)));
        os << std::setw(5) << res_num << " "
                << std::setw(5) << std::left
                << residue_name[topo.solvent(s).atom(a).residue_nr] << " "
                << std::setw(6) << std::left << topo.solvent(s).atom(a).name << std::right
                << std::setw(6) << index + 1
                << std::setw(m_width) << vel_rot(0)
                << std::setw(m_width) << vel_rot(1)
                << std::setw(m_width) << vel_rot(2)
                << "\n";
      }
    }
  }

  os << "END\n";

}

void io::Out_Configuration
::_print_lattice_shifts(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  os << "LATTICESHIFTS\n";
  math::VArray const &shift = conf.special().lattice_shifts;
  for (int i = 0, to = topo.num_atoms(); i < to; ++i) {
    os << std::setw(10) << int(rint(shift(i)(0)))
            << std::setw(10) << int(rint(shift(i)(1)))
            << std::setw(10) << int(rint(shift(i)(2)))
            << "\n";
  }
  os << "END\n";
}

void io::Out_Configuration
::_print_velocityred(configuration::Configuration const &conf,
        int num, std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  os << "VELOCITYRED\n";

  math::VArray const &vel = conf.current().vel;
  //matrix to rotate back into orignial Cartesian Coordinat system
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));
  if (conf.boundary_type == math::truncoct)
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));
  math::Vec vel_rot;

  assert(num <= int(vel.size()));
  for (int i = 0; i < num; ++i) {
    vel_rot = math::Vec(math::product(Rmat, vel(i)));
    os << std::setw(m_width) << vel_rot(0)
            << std::setw(m_width) << vel_rot(1)
            << std::setw(m_width) << vel_rot(2)

            << "\n";
    if ((i + 1) % 10 == 0) os << '#' << std::setw(10) << i + 1 << "\n";
  }

  os << "END\n";
}

void io::Out_Configuration
::_print_force(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os, bool constraint_force) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_force_precision);

  os << "FREEFORCE\n";

  math::VArray const & force = conf.current().force;
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();
  //matrix to rotate back into orignial Cartesian Coordinat system
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));
  if (conf.boundary_type == math::truncoct)
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));
  math::Vec force_rot;
  os << "# first 24 chars ignored\n";
  for (int i = 0, to = topo.num_solute_atoms(); i < to; ++i) {
    force_rot = math::Vec(math::product(Rmat, force(i)));
    os << std::setw(6) << solute.atom(i).residue_nr + 1
            << std::setw(5) << residue_name[solute.atom(i).residue_nr]
            << std::setw(6) << solute.atom(i).name
            << std::setw(8) << i + 1
            << std::setw(m_force_width) << force_rot(0)
            << std::setw(m_force_width) << force_rot(1)
            << std::setw(m_force_width) << force_rot(2)
            << "\n";
  }
  int index = topo.num_solute_atoms();

  for (unsigned int s = 0; s < topo.num_solvents(); ++s) {
    for (unsigned int m = 0; m < topo.num_solvent_molecules(s); ++m) {
      for (unsigned int a = 0; a < topo.solvent(s).num_atoms(); ++a, ++index) {
        force_rot = math::Vec(math::product(Rmat, force(index)));
        os << std::setw(6) << topo.solvent(s).atom(a).residue_nr + 1
                << std::setw(5) << residue_name[topo.solvent(s).atom(a).residue_nr]
                << std::setw(6) << topo.solvent(s).atom(a).name
                << std::setw(8) << index + 1
                << std::setw(m_force_width) << force_rot(0)
                << std::setw(m_force_width) << force_rot(1)
                << std::setw(m_force_width) << force_rot(2)
                << "\n";
      }
    }
  }
  os << "END\n";

  if (constraint_force) {
    os << "CONSFORCE\n";

    math::VArray const & cons_force = conf.current().constraint_force;

    os << "# first 24 chars ignored\n";
    for (int i = 0, to = topo.num_solute_atoms(); i < to; ++i) {
      force_rot = math::Vec(math::product(Rmat, cons_force(i)));
      os << std::setw(6) << solute.atom(i).residue_nr + 1
              << std::setw(5) << residue_name[solute.atom(i).residue_nr]
              << std::setw(6) << solute.atom(i).name
              << std::setw(8) << i + 1
              << std::setw(m_force_width) << force_rot(0)
              << std::setw(m_force_width) << force_rot(1)
              << std::setw(m_force_width) << force_rot(2)
              << "\n";
    }
    index = topo.num_solute_atoms();

    for (unsigned int s = 0; s < topo.num_solvents(); ++s) {
      for (unsigned int m = 0; m < topo.num_solvent_molecules(s); ++m) {
        for (unsigned int a = 0; a < topo.solvent(s).num_atoms(); ++a, ++index) {
           force_rot = math::Vec(math::product(Rmat, cons_force(index)));
          os << std::setw(6) << topo.solvent(s).atom(a).residue_nr + 1
                  << std::setw(5) << residue_name[topo.solvent(s).atom(a).residue_nr]
                  << std::setw(6) << topo.solvent(s).atom(a).name
                  << std::setw(8) << index + 1
                  << std::setw(m_force_width) << force_rot(0)
                  << std::setw(m_force_width) << force_rot(1)
                  << std::setw(m_force_width) << force_rot(2)
                  << "\n";
        }
      }
    }
    os << "END\n";
  }
}

void io::Out_Configuration
::_print_forcered(configuration::Configuration const &conf,
        int num,
        std::ostream &os,
        bool constraint_force) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_force_precision);

  os << "FREEFORCERED\n";

  const math::VArray &force = conf.old().force;
//matrix to rotate back into orignial Cartesian Coordinat system
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));
  if (conf.boundary_type == math::truncoct)
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));
  math::Vec force_rot;
  
  assert(num <= int(force.size()));

  for (int i = 0; i < num; ++i) {
    force_rot = math::Vec(math::product(Rmat, force(i)));
    os << std::setw(m_force_width) << force_rot(0)
            << std::setw(m_force_width) << force_rot(1)
            << std::setw(m_force_width) << force_rot(2)
            << "\n";

    if ((i + 1) % 10 == 0) os << '#' << std::setw(10) << i + 1 << "\n";
  }

  os << "END\n";
  // group wise forces
  if (conf.special().force_groups.size()) {
    os << "FREEFORCEGROUPRED\n";
    for(unsigned int i = 0; i < conf.special().force_groups.size(); ++i) {
      for (unsigned int j = i; j < conf.special().force_groups.size(); ++j) {
        os << "# group " << i+1 << "-" << j+1 << "\n";
        for (int k = 0; k < num; ++k) {
          force_rot = math::Vec(math::product(Rmat, conf.special().force_groups[j][i](k)));
          os << std::setw(m_force_width) << force_rot(0)
                  << std::setw(m_force_width) << force_rot(1)
                  << std::setw(m_force_width) << force_rot(2)
                  << "\n";

          if ((k + 1) % 10 == 0) os << '#' << std::setw(10) << k + 1 << "\n";
        }
      }
    }
    os << "END\n";
  }

  if (constraint_force) {
    os << "CONSFORCERED\n";

    const math::VArray & cons_force = conf.old().constraint_force;
    assert(num <= int(cons_force.size()));
    math::Vec cons_force_rot;
    for (int i = 0; i < num; ++i) {
      cons_force_rot = math::Vec(math::product(Rmat, cons_force(i)));
      os << std::setw(m_force_width) << cons_force_rot(0)
              << std::setw(m_force_width) << cons_force_rot(1)
              << std::setw(m_force_width) << cons_force_rot(2)
              << "\n";

      if ((i + 1) % 10 == 0) os << '#' << std::setw(10) << i + 1 << "\n";
    }

    os << "END\n";
  }
}

/**
 * 
 * @section energyredhelper ENERGY03
 *
@verbatim
ENERGY03
# totals
   1.598443703e+02 # total
   1.709320077e+02 # kinetic
  -4.205278511e+01 # potential total
   3.980457903e+01 # covalent total
   1.935607149e+01 # bonds total
   1.267041319e+01 # angles total
   1.485470503e+00 # impropers total
   6.292623846e+00 # dihedrals total
   0.000000000e+00 # crossdihedrals total
  -8.185736413e+01 # nonbonded total
  -1.442773701e+00 # Lennard-Jones total
  -8.041459043e+01 # Coulomb/Reaction-Field total
   0.000000000e+00 # lattice total
   0.000000000e+00 # lattice sum pair total
   0.000000000e+00 # lattice sum real space total
   0.000000000e+00 # lattice sum k (reciprocal) space total
   0.000000000e+00 # lattice sum A term total
   0.000000000e+00 # lattice sum self total
   0.000000000e+00 # lattice sum surface total
   0.000000000e+00 # polarisation self total
   0.000000000e+00 # special total
   0.000000000e+00 # SASA total
   0.000000000e+00 # SASA volume total
   0.000000000e+00 # constraints total
   0.000000000e+00 # distance restraints total
   0.000000000e+00 # distancefield restraints total
   0.000000000e+00 # dihedral restraints total
   0.000000000e+00 # position restraints total
   0.000000000e+00 # J-value restraints total
   0.000000000e+00 # X-ray restraints total
   0.000000000e+00 # Local elevation total
   0.000000000e+00 # order parameter restraints total
   0.000000000e+00 # symmetry restraints total
   3.096514777e+01 # EDS: energy of reference state
   0.000000000e+00 # Entropy
   0.000000000e+00 # QM/MM total
   0.000000000e+00 # B&S-LEUS energy
   0.000000000e+00 # RDC-value total
   0.000000000e+00 # angle restraints total
# baths
# number of baths
2
#  kinetic total     centre of mass    internal/rotational
   2.579417475e+01   8.149832229e-01   2.497919152e+01 # 1-st bath
   1.451378329e+02   7.403204119e+01   7.110579174e+01 # 2-nd bath
# bonded
# number of energy groups
2
#  bond              angle             improper          dihedral        crossdihedral
   1.935607149e+01   1.267041319e+01   1.485470503e+00   6.292623846e+00 0.0000000e+00  # energy group 1
   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00 0.0000000e+00  # energy group 2
# nonbonded
#  Lennard-Jones     Coulomb/RF        lattice sum real  lattice sum reciproc.
  -4.896399521e-02  -6.781355366e+01   0.000000000e+00   0.000000000e+00  # 1 - 1
  -6.774710271e-01   4.920740888e-01   0.000000000e+00   0.000000000e+00  # 1 - 2
  -7.163386790e-01  -1.309311086e+01   0.000000000e+00   0.000000000e+00  # 2 - 2
# special
#  constraints       pos. restraints   dist. restraints  disfield res      dihe. restr.      SASA              SASA volume       jvalue            rdc               local elevation   path integral   angle restraint
   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00 0.000000000e+00 # group 1
   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00   0.000000000e+00 0.000000000e+00 # group 2
# eds (enveloping distribution sampling)
# numstates
2
           # total         nonbonded          special
   3.096514777e+01   3.096514777e+01   0.000000000e+00
   3.096514777e+01   3.096514777e+01   0.000000000e+00
END
@endverbatim
 *
 */
void io::Out_Configuration
::_print_energyred(configuration::Configuration const &conf,
        std::ostream &os) {
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  os << "ENERGY03\n";
  _print_energyred_helper(os, conf.old().energies);
  os << "END\n";

}

void io::Out_Configuration
::_print_free_energyred(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  // assert(m_free_energy_traj.is_open());

  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  os << "FREEENERDERIVS03\n"
          << "# lambda\n"
          << std::setw(18) << topo.old_lambda() << "\n";

  _print_energyred_helper(os, conf.old().perturbed_energy_derivatives);

  os << "END\n";

}

/**
 * 
 * @section volumepressurered VOLUMEPRESSURE03
 *
@verbatim
VOLUMEPRESSURE03
# mass
   5.044822000e+02
# temperature
# number of temperature coupling baths
2
#  total             com               ir                scaling factor
   1.723525431e+02   6.534704791e+01   1.820803154e+02   1.128031828e+00 # 1st bath
   2.909363241e+02   2.968021431e+02   2.850705051e+02   1.002070048e+00 # 2nd bath
# volume
   5.345718909e+01
# box
   3.767055681e+00   0.000000000e+00   0.000000000e+00 # K
   0.000000000e+00   3.767055681e+00   0.000000000e+00 # L
   0.000000000e+00   0.000000000e+00   3.767055681e+00 # M
# pressure
   7.348382177e-01
   5.261908890e+00
   2.490310167e+01
#  pressure tensor
   3.260036299e-01  -3.892307215e-01   3.337546495e-01
  -4.705767676e-02   8.540193624e-01   1.012276325e-01
   1.604604267e-01   2.852680401e-01   1.024491661e+00
#  virial tensor
   1.508936821e+01   1.501353732e+01   1.232971839e+00
   5.867732745e+00   4.634376722e+00   1.731944968e+00
   5.864882856e+00  -3.187196467e+00  -3.938018259e+00
#  molecular kinetic energy tensor
   2.380298705e+01   4.609947183e+00   1.015376454e+01
   4.609947183e+00   2.746111399e+01   4.437617314e+00
   1.015376454e+01   4.437617314e+00   2.344520396e+01
@endverbatim
 *
 */
void io::Out_Configuration
::_print_volumepressurered(topology::Topology const & topo,
        configuration::Configuration const &conf,
        simulation::Simulation const &sim,
        std::ostream &os) {
  std::vector<double> const s;

  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  os << "VOLUMEPRESSURE03\n";

  _print_volumepressurered_helper(os,
          math::sum(topo.mass()),
          conf.old().phi,conf.old().theta, conf.old().psi,
          sim.multibath(),
          s,
          conf.old().energies,
          conf.old().box,
          conf.boundary_type,
          conf.old().pressure_tensor,
          conf.old().virial_tensor,
          conf.old().kinetic_energy_tensor);

  os << "END\n";

}

void io::Out_Configuration
::_print_box(configuration::Configuration const &conf,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  //change to GENBOX
  os << "GENBOX\n";

  math::Box box = conf.current().box;

  os << std::setw(5) << conf.boundary_type << "\n";

  long double a = 0.0, b = 0.0, c = 0.0, alpha = 0.0, beta = 0.0, gamma = 0.0, phi = 0.0, theta = 0.0, psi = 0.0;
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));
  
  math::Box m(0.0);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        m(j)(i) += Rmat(k, i) * box(j)(k);
  box = m;

  // convert it back to truncoct
  if (conf.boundary_type == math::truncoct)
    math::truncoct_triclinic_box(box, false);
  
  a = math::abs(box(0));
  b = math::abs(box(1));
  c = math::abs(box(2));

  os << std::setw(m_width) << a
          << std::setw(m_width) << b
          << std::setw(m_width) << c
          << "\n";


  if (a != 0.0 && b != 0.0 && c != 0.0) {
    alpha = acos(math::costest(dot(box(1), box(2)) / (b * c)));
    beta = acos(math::costest(dot(box(0), box(2)) / (a * c)));
    gamma = acos(math::costest(dot(box(0), box(1)) / (a * b)));

    os << std::setw(m_width) << alpha * 180 / math::Pi
            << std::setw(m_width) << beta * 180 / math::Pi
            << std::setw(m_width) << gamma * 180 / math::Pi
            << "\n";

    /* we are already in the frame of the box!
        math::Matrixl Rmat = (math::rmat(box));

        long double R11R21 = sqrtl(Rmat(0, 0) * Rmat(0, 0) + Rmat(0, 1) * Rmat(0, 1));
        if (R11R21 == 0.0) {
          theta = -math::sign(Rmat(0, 2)) * M_PI / 2;
          psi = 0.0;
          phi = -math::sign(Rmat(1, 0)) * acosl(math::costest(Rmat(1, 1)));
        } else {
          theta = -math::sign(Rmat(0, 2)) * acosl(math::costest(R11R21));
          long double costheta = cosl(theta);
          psi = math::sign(Rmat(1, 2) / costheta) * acosl(math::costest(Rmat(2, 2) / costheta));
          phi = math::sign(Rmat(0, 1) / costheta) * acosl(math::costest(Rmat(0, 0) / costheta));

        }

     */
  
    if(fabs(conf.current().phi)<math::epsilon)
      phi=0.0;
    else phi=conf.current().phi;
    if(fabs(conf.current().theta)<math::epsilon)
      theta=0.0;
    else theta=conf.current().theta;
    if(fabs(conf.current().psi)<math::epsilon)
      psi=0.0;
    else psi=conf.current().psi;
    os << std::setw(m_width) << phi * 180 / math::Pi
            << std::setw(m_width) << theta * 180 / math::Pi
            << std::setw(m_width) << psi * 180 / math::Pi
            << "\n";
  } else {
    os << std::setw(m_width) << 0.0
            << std::setw(m_width) << 0.0
            << std::setw(m_width) << 0.0
            << "\n";
    os << std::setw(m_width) << 0.0
            << std::setw(m_width) << 0.0
            << std::setw(m_width) << 0.0
            << "\n";
  }
  double origin = 0.0;

  os << std::setw(m_width) << origin
          << std::setw(m_width) << origin
          << std::setw(m_width) << origin
          << "\n";


  os << "END\n";

}

void io::Out_Configuration
::precision(int prec, int add) {
  m_precision = prec;
  m_width = prec + add;
}

void io::Out_Configuration
::force_precision(int prec, int add) {
  m_force_precision = prec;
  m_force_width = prec + add;
}

int io::Out_Configuration
::precision() {
  return m_precision;
}

int io::Out_Configuration
::force_precision() {
  return m_force_precision;
}

void io::Out_Configuration
::print(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation const & sim) {
  if (sim.param().print.stepblock && (sim.steps() % sim.param().print.stepblock) == 0) {

    m_output << "\n---------------------------------------------------"
            << "-----------------------------\n";

    _print_timestep(sim, m_output);

    //print_MINIMISATION(m_output, rmsd_force, max_force, total_iterations);

    print_ENERGY(m_output, conf.old().energies, topo.energy_groups());

    if (sim.param().sasa.switch_sasa)
      print_sasa(m_output, topo, conf, sim,
              "SOLVENT ACCESSIBLE SURFACE AREAS AND VOLUMES OF BURIED ATOMS");

    if (sim.param().perturbation.perturbation) {

      m_output << "lambda: " << topo.old_lambda() << "\n";

      print_ENERGY(m_output, conf.old().perturbed_energy_derivatives,
              topo.energy_groups(), "dE/dLAMBDA", "dE_");
    }
    if (!sim.param().minimise.ntem) {
      print_MULTIBATH(m_output, sim.multibath(), conf.old().energies);
    }

    // flexible shake kinetic energy
    if (sim.param().constraint.solute.algorithm == simulation::constr_flexshake) {
      m_output << "FLEXSHAKE\n";
      m_output << "\tflex_ekin";
      for (unsigned int i = 0; i < conf.special().flexible_constraint.flexible_ekin.size(); ++i)
        m_output << std::setw(12) << std::setprecision(4) << std::scientific
              << conf.special().flexible_constraint.flexible_ekin[i];
      m_output << "\nEND\n";
    }

    if (sim.param().pcouple.calculate)
      print_PRESSURE(m_output, conf);

    m_output.flush();

  }
}

void io::Out_Configuration
::print_final(topology::Topology const & topo,
        configuration::Configuration & conf,
        simulation::Simulation const & sim) {
  m_output << "\n============================================================\n";
  m_output << "FINAL DATA\n";
  m_output << "============================================================\n\n\n";

  m_output << "\tsimulation time  :" << std::setw(10) << sim.time() << "\n"
          << "\tsimulation steps :" << std::setw(10) << sim.steps() << "\n\n";

  configuration::Energy e, ef;
  math::Matrix p, pf, v, vf, et, etf;

  std::vector<double> sasa_a, sasa_af, sasa_vol, sasa_volf;
  double sasa_tot = 0.0, sasa_totf = 0.0, sasa_voltot = 0.0, sasa_voltotf = 0.0;

  if (sim.param().minimise.ntem) {
    //print_MINIMISATION(m_output, rmsd_force, max_force, sim.minimisation);
    print_ENERGY(m_output, conf.current().energies, topo.energy_groups(), "MINIMIZED ENERGY",
            "<EMIN>_");
  }

  // new averages
  conf.current().averages.simulation().energy_average(e, ef);
  conf.current().averages.simulation().pressure_average(p, pf, v, vf, et, etf);

  // averages and fluctuation for sasa and volume calculation
  if (sim.param().sasa.switch_sasa) {
    conf.current().averages.simulation().sasa_average(sasa_a, sasa_af, sasa_tot, sasa_totf);
    if (sim.param().sasa.switch_volume)
      conf.current().averages.simulation().sasavol_average(sasa_vol, sasa_volf, sasa_voltot,
              sasa_voltotf);
  }

  print_ENERGY(m_output, e, topo.energy_groups(), "ENERGY AVERAGES", "<E>_");
  print_ENERGY(m_output, ef, topo.energy_groups(), "ENERGY FLUCTUATIONS", "<<E>>_");
  if (!sim.param().minimise.ntem) {
    print_MULTIBATH(m_output, sim.multibath(), e, "TEMPERATURE AVERAGES");
    print_MULTIBATH(m_output, sim.multibath(), ef, "TEMPERATURE FLUCTUATIONS");
  }

  if (sim.param().pcouple.calculate) {
    print_MATRIX(m_output, p, "PRESSURE AVERAGE");
    print_MATRIX(m_output, pf, "PRESSURE FLUCTUATION");
  }

  if (sim.param().perturbation.perturbation) {

    double lambda = 0.0, lambda_fluct = 0.0;
    conf.current().averages.simulation().
            energy_derivative_average(e, ef, lambda, lambda_fluct, sim.param().perturbation.dlamt);

    if (sim.param().perturbation.dlamt) {

      print_ENERGY(m_output, e, topo.energy_groups(), "CUMULATIVE DG", "DG_");

      // what's that anyway...
      //print_ENERGY(m_output, ef, topo.energy_groups(), "DG FLUCTUATIONS", "<<DG>>_");
    }
    else {

      std::ostringstream ss, pre;
      ss << "dE/dLAMBDA ";
      pre << "dE/dl";

      print_ENERGY(m_output, e, topo.energy_groups(),
              ss.str() + "AVERAGES", "<" + pre.str() + ">_");


      print_ENERGY(m_output, ef, topo.energy_groups(),
              ss.str() + "FLUCTUATIONS", "<<" + pre.str() + ">>_");

    }
  }

  // print sasa and volume averages, fluctuations
  if (sim.param().sasa.switch_sasa){
    int volume = sim.param().sasa.switch_volume;
    print_sasa_avg(m_output, topo, sasa_a, sasa_vol, sasa_tot, sasa_voltot, 
            "SASA AND BURIED VOLUME AVERAGES", volume);
    print_sasa_fluct(m_output, topo, sasa_af, sasa_volf, sasa_totf, sasa_voltotf, 
            "SASA AND BURIED VOLUME FLUCTUATIONS", volume);
  }

}

void io::Out_Configuration
::_print_blockaveraged_energyred(configuration::Configuration const &conf,
        std::ostream &os) {
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  configuration::Energy e, ef;
  // energies are in old(), but averages in current()!
  conf.current().averages.block().energy_average(e, ef);

  os << "BAENERGY03\n";
  _print_energyred_helper(os, e);
  os << "END\n";

  os << "BAEFLUCT03\n";
  _print_energyred_helper(os, ef);
  os << "END\n";

}

void io::Out_Configuration
::_print_blockaveraged_volumepressurered(configuration::Configuration const & conf,
        simulation::Simulation const & sim,
        std::ostream &os) {
  double mass = 0.0, massf = 0.0, vol = 0.0, volf = 0.0;
  std::vector<double> s, sf;
  configuration::Energy e, ef;
  math::Box b, bf;
  math::Matrix p, pf, v, vf, k, kf;

  // needed again. do in 1 function (energyred && volumepressurered)
  // energies in old(), but averages in current()!
  conf.current().averages.block().energy_average(e, ef);
  conf.current().averages.block().pressure_average(p, pf, v, vf, k, kf);
  conf.current().averages.block().mbs_average(mass, massf, vol, volf, b, bf, s, sf);

  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  os << "BAVOLUMEPRESSURE03\n";

  _print_volumepressurered_helper(os,
          mass,
          conf.old().phi,conf.old().theta, conf.old().psi,
          sim.multibath(),
          s,
          e,
          b,
          conf.boundary_type,
          p,
          v,
          k);

  os << "END\n";

  os << "BAVPFLUCT03\n";

  _print_volumepressurered_helper(os,
          massf,
          conf.old().phi,conf.old().theta, conf.old().psi,
          sim.multibath(),
          sf,
          ef,
          bf,
          conf.boundary_type,
          pf,
          vf,
          kf);

  os << "END\n";

}

void io::Out_Configuration
::_print_blockaveraged_free_energyred(configuration::Configuration const &conf,
        double dlamt,
        std::ostream &os) {
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);

  configuration::Energy e, ef;
  double lambda = 0.0, lambda_fluct = 0.0;

  // energies in old(), but averages in current()!
  conf.current().averages.block().energy_derivative_average(e, ef, lambda, lambda_fluct, dlamt);

  os << "BAFREEENERDERIVS03\n"
          << "# lambda\n"
          << std::setw(18) << lambda << "\n";

  _print_energyred_helper(os, e);

  os << "END\n";

  os << "BAFREEEFLUCTS03\n"
          << "# lambda\n"
          << std::setw(18) << lambda_fluct << "\n";

  _print_energyred_helper(os, ef);

  os << "END\n";

}

void io::Out_Configuration::_print_flexv(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  DEBUG(10, "FLEXV");

  unsigned int k = 0;

  std::vector<topology::two_body_term_struct>::const_iterator
  constr_it = topo.solute().distance_constraints().begin(),
          constr_to = topo.solute().distance_constraints().end();

  os << "FLEXV\n";
  os << "#\tflexible constraints ("
          << topo.solute().distance_constraints().size()
          << ")\n";

  for (; constr_it != constr_to; ++constr_it, ++k) {

    assert(conf.special().flexible_constraint.flex_len.size() > k);
    assert(conf.special().flexible_constraint.flexible_vel.size() > k);

    os << std::setw(15) << constr_it->i + 1
            << std::setw(10) << constr_it->j + 1
            << std::setw(20) << conf.special().flexible_constraint.flex_len[k]
            << std::setw(20) << conf.special().flexible_constraint.flexible_vel[k]
            << "\n";
  }

  std::vector<topology::perturbed_two_body_term_struct>::const_iterator
  pconstr_it = topo.perturbed_solute().distance_constraints().begin(),
          pconstr_to = topo.perturbed_solute().distance_constraints().end();

  os << "#\tperturbed flexible constraints ("
          << topo.perturbed_solute().distance_constraints().size()
          << " of "
          << conf.special().flexible_constraint.flex_len.size()
          << ")\n";

  for (; pconstr_it != pconstr_to; ++pconstr_it, ++k) {

    assert(conf.special().flexible_constraint.flex_len.size() > k);
    assert(conf.special().flexible_constraint.flexible_vel.size() > k);

    os << std::setw(15) << pconstr_it->i + 1
            << std::setw(10) << pconstr_it->j + 1
            << std::setw(20) << conf.special().flexible_constraint.flex_len[k]
            << std::setw(20) << conf.special().flexible_constraint.flexible_vel[k]
            << "\n";
  }

  os << "END\n";

}

void io::Out_Configuration::_print_stochastic_integral(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();

  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision - 2); // because it's scientific

  os << "STOCHINT\n";
  os << "# first 24 chars ignored\n";

  for (int i = 0, to = topo.num_solute_atoms(); i < to; ++i) {
    os.setf(std::ios::fixed, std::ios::floatfield);
    os << std::setw(5) << solute.atom(i).residue_nr + 1 << " "
            << std::setw(5) << std::left
            << residue_name[solute.atom(i).residue_nr] << " "
            << std::setw(6) << std::left << solute.atom(i).name << std::right
            << std::setw(6) << i + 1;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os << std::setw(m_width) << conf.current().stochastic_integral(i)(0)
            << std::setw(m_width) << conf.current().stochastic_integral(i)(1)
            << std::setw(m_width) << conf.current().stochastic_integral(i)(2)
            << "\n";
  }

  int index = topo.num_solute_atoms();
  int res_nr = 1;

  for (unsigned int s = 0; s < topo.num_solvents(); ++s) {

    for (unsigned int m = 0; m < topo.num_solvent_molecules(s); ++m, ++res_nr) {

      for (unsigned int a = 0; a < topo.solvent(s).num_atoms(); ++a, ++index) {
        os.setf(std::ios::fixed, std::ios::floatfield);
        os << std::setw(5) << res_nr
                << ' ' << std::setw(5) << std::left
                << residue_name[topo.solvent(s).atom(a).residue_nr] << " "
                << std::setw(6) << std::left
                << topo.solvent(s).atom(a).name << std::right
                << std::setw(6) << index + 1;
        os.setf(std::ios::scientific, std::ios::floatfield);
        os << std::setw(m_width) << conf.current().stochastic_integral(index)(0)
                << std::setw(m_width) << conf.current().stochastic_integral(index)(1)
                << std::setw(m_width) << conf.current().stochastic_integral(index)(2)
                << "\n";
      }
    }
  }
  os.setf(std::ios::fixed, std::ios::floatfield);
  os << "# seed\n" << std::setw(10) << std::right
          << conf.current().stochastic_seed << "\n";
  os << "END\n";
}

void io::Out_Configuration::_print_pertdata(topology::Topology const &topo,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(7);
  os << "PERTDATA" << std::endl
          << std::setw(15) << topo.lambda() << std::endl
          << "END" << std::endl;
}

void io::Out_Configuration::_print_jvalue(simulation::Parameter const & param,
        configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os, bool formated) {
  DEBUG(10, "JVALUE Averages and LE data");

  if (param.jvalue.mode != simulation::jvalue_restr_inst &&
      param.jvalue.mode != simulation::jvalue_restr_inst_weighted) {
    os << "JVALUERESEXPAVE\n";
    if (formated) {
      os << "# I J K L <J>\n";
    }
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(7);
    std::vector<double>::const_iterator av_it = conf.special().jvalue_av.begin(),
            av_to = conf.special().jvalue_av.end();
    std::vector<topology::jvalue_restraint_struct>::const_iterator jv_it =
            topo.jvalue_restraints().begin();
    for (; av_it != av_to; ++av_it, ++jv_it) {
      if (formated) {
        os << std::setw(5) << jv_it->i+1
           << std::setw(5) << jv_it->j+1
           << std::setw(5) << jv_it->k+1
           << std::setw(5) << jv_it->l+1;
      }
      os << std::setw(15) << *av_it << "\n";
    }
    os << "END\n";
  }

  if (param.jvalue.le && param.jvalue.mode != simulation::jvalue_restr_off) {
    os << "JVALUERESEPS\n";
    if (formated) {
      os << "# I J K L\n# GRID[1.." << param.jvalue.ngrid << "]\n";
    }
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(7);
    std::vector<std::vector<double> >::const_iterator
    le_it = conf.special().jvalue_epsilon.begin(),
            le_to = conf.special().jvalue_epsilon.end();
    std::vector<topology::jvalue_restraint_struct>::const_iterator jv_it =
            topo.jvalue_restraints().begin();

    for (; le_it != le_to; ++le_it, ++jv_it) {
      if (formated) {
        os << std::setw(5) << jv_it->i+1
           << std::setw(5) << jv_it->j+1
           << std::setw(5) << jv_it->k+1
           << std::setw(5) << jv_it->l+1 << "\n";
      }
      for (unsigned int i = 0; i < le_it->size(); ++i)
        os << std::setw(15) << (*le_it)[i];
      os << "\n";
    }
    os << "END\n";
  }
}

void io::Out_Configuration::_print_xray(simulation::Parameter const & param,
        configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os, bool final) {
  DEBUG(10, "XRAY Averages");

  if (param.xrayrest.xrayrest != simulation::xrayrest_inst) {
    os << "XRAYRESEXPAVE\n";
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(7);
    for (unsigned int i=0; i<conf.special().xray_rest.size(); ++i) {
      os << std::setw(15) << conf.special().xray_rest[i].sf_av
         << std::setw(15) << conf.special().xray_rest[i].phase_av << "\n";
    }
    os << "END\n";
  }
}

void io::Out_Configuration::_print_xray_rvalue(simulation::Parameter const & param,
        configuration::Configuration const &conf,
        std::ostream &os) {
  DEBUG(10, "XRAY scaling constants and R-values");

  double k_inst = conf.special().xray.k_inst;
  double k_free_inst = conf.special().xray.k_free_inst;
  double k_avg  = conf.special().xray.k_avg;
  double k_free_avg  = conf.special().xray.k_free_avg;
  double R_inst = conf.special().xray.R_inst;
  double R_avg  = conf.special().xray.R_avg;
  double R_free_inst = conf.special().xray.R_free_inst;
  double R_free_avg  = conf.special().xray.R_free_avg;

  // make sure no rubbish is written
  switch(param.xrayrest.xrayrest) {
    case simulation::xrayrest_off : return;
    case simulation::xrayrest_inst :
      k_avg = k_free_avg = R_avg = R_free_avg = 0.0; break;
    case simulation::xrayrest_avg :
      k_inst = k_free_inst = R_inst = R_free_inst = 0.0; break;
    default: ;// value are OK. do nothing
  }
  os << "XRAYRVALUE\n";
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  os << std::setw(15) << k_inst << std::endl
          << std::setw(15) << R_inst << std::endl
          << std::setw(15) << k_free_inst << std::endl
          << std::setw(15) << R_free_inst << std::endl
          << std::setw(15) << k_avg << std::endl
          << std::setw(15) << R_avg << std::endl
          << std::setw(15) << k_free_avg << std::endl
          << std::setw(15) << R_free_avg << std::endl;
  os << "END\n";
}

void io::Out_Configuration::_print_xray_umbrellaweightthresholds(simulation::Parameter const & param,
        topology::Topology const & topo,
        std::ostream &os) {
  DEBUG(10, "XRAY umbrella weight thresholds");

  if (param.xrayrest.local_elevation == false) return;

  os << "XRAYUMBRELLAWEIGHTTHRESHOLDS\n";
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  std::vector<topology::xray_umbrella_weight_struct>::const_iterator it =
            topo.xray_umbrella_weights().begin(), to =
            topo.xray_umbrella_weights().end();
  for (; it != to; ++it) {
    os << std::setw(15) << it->threshold
            << std::setw(15) << it->threshold_growth_rate
            << std::setw(15) << it->threshold_overshoot
            << std::setw(3) << (it->threshold_freeze ? 1 : 0) << std::endl;
  }
  os << "END\n";
}

void io::Out_Configuration::_print_xray_bfactors(simulation::Parameter const & param,
        configuration::Configuration const & conf,
        std::ostream &os) {
  DEBUG(10, "XRAY b factors");

  if (param.xrayrest.bfactor.step == 0) return;

  os << "XRAYBFOCCSPEC\n";
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  std::vector<configuration::Configuration::special_struct::xray_bfoc_struct>::const_iterator
  it = conf.special().xray_bfoc.begin(), to = conf.special().xray_bfoc.end();
  for (; it != to; ++it) {
    os << std::setw(15) << it->b_factor
            << std::setw(15) << it->occupancy << std::endl;
  }
  os << "END\n";
}

void io::Out_Configuration::_print_distance_restraints(
        configuration::Configuration const &conf,
        topology::Topology const & topo,
        std::ostream &os) {
  DEBUG(10, "distance restraints");

  std::vector<double>::const_iterator av_it = conf.special().distanceres.av.begin(),
          av_to = conf.special().distanceres.av.end();
  std::vector<double>::const_iterator ene_it = conf.special().distanceres.energy.begin();
  std::vector<double>::const_iterator d_it = conf.special().distanceres.d.begin();

  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_distance_restraint_precision);
  
  if (conf.special().distanceres.d.size() > 0) {
    os << "DISRESDATA" << std::endl;
    os << conf.special().distanceres.d.size() << "\n";
    int i = 0;
    for (i = 1; av_it != av_to; ++av_it, ++ene_it, ++d_it, ++i) {
       os << std::setw(m_width) << *d_it
       << std::setw(m_width) << *ene_it;
       if (*av_it != 0) {
         os << std::setw(m_width) << pow(*av_it, -1.0 / 3.0);
       } else {
         os << std::setw(m_width) << 0.0;
       }
       os << std::endl;
    }
    os << "END" << std::endl;
  }
  
  std::vector<double>::const_iterator pav_it = conf.special().pertdistanceres.av.begin(),
          pav_to = conf.special().pertdistanceres.av.end();
  std::vector<double>::const_iterator pene_it = conf.special().pertdistanceres.energy.begin();
  std::vector<double>::const_iterator pd_it = conf.special().pertdistanceres.d.begin();
  
  if (conf.special().pertdistanceres.d.size() > 0) {
    os << "PERTDISRESDATA" << std::endl;
    os << conf.special().pertdistanceres.d.size() << "\n";
    int i = 0;
    for (i = 1; pav_it != pav_to; ++pav_it, ++pene_it, ++pd_it, ++i) {
       os << std::setw(m_width) << *pd_it
       << std::setw(m_width) << *pene_it;
       if (*pav_it != 0) {
         os << std::setw(m_width) << pow(*pav_it, -1.0 / 3.0);
       } else {
         os << std::setw(m_width) << 0.0;
       }
       os << std::endl;
    }
    os << "END" << std::endl;
  }
}

void io::Out_Configuration::_print_disfield_restraints(
        configuration::Configuration const &conf,
        topology::Topology const & topo,
        std::ostream &os) {
  DEBUG(10, "disfield restraints");

  double distance = conf.special().distancefield.dist;
  double energy = conf.special().distancefield.energy;
  
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_disfield_restraint_precision);

  os << "DISFIELDRESDATA" << std::endl;
  os << std::setw(m_width) << distance
     << std::setw(m_width) << energy
     << std::endl;
  os << "END" << std::endl;
}

void io::Out_Configuration::_print_distance_restraint_averages(
        configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  DEBUG(10, "distance restraint averages");

  std::vector<double>::const_iterator it = conf.special().distanceres.av.begin(),
          to = conf.special().distanceres.av.end();

  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_distance_restraint_precision);

  os << "DISRESEXPAVE" << std::endl;
  int i = 0;
  for (i = 1; it != to; ++it, ++i) {
    os << std::setw(m_width) << *it;
    if (i % 5 == 0)
      os << std::endl;
  }
  if ((i-1) % 5 != 0)
    os << std::endl;

  os << "END" << std::endl;
}

void io::Out_Configuration::_print_disfield_grid(
	simulation::Parameter const &param,
        configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  DEBUG(10, "distancefield grid");
  
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_distance_restraint_precision);

  os << "DISTANCEFIELDGRID" << std::endl;
 
  const std::vector<int> &ngrid = conf.special().distancefield.ngrid;
  const double grid = param.distancefield.grid;
  math::Box box = conf.current().box;
  
  for(unsigned int j=0; j< conf.special().distancefield.distance.size(); j++){ 
    int nz = int(double(j) / double(ngrid[0] * ngrid[1]));
    int ny = int(double( j % (ngrid[0] * ngrid[1]) ) / double(ngrid[0]));
    int nx = (j % (ngrid[0] * ngrid[1])) % ngrid[0];
    // the grid is defined for -half the box lenght to +half the box length
    math::Vec gpos(nx*grid - box(0,0)/2, ny*grid - box(1,1)/2, nz*grid - box(2,2)/2);
    
    // we put it in the positive box in an ugly way because:
    // 1. distancefield only works for rectangular boxes anyway
    // 2. otherwise this function has to be templated for periodicity
    if(gpos[0] < 0.0) gpos[0] += box(0,0);
    if(gpos[1] < 0.0) gpos[1] += box(1,1);
    if(gpos[2] < 0.0) gpos[2] += box(2,2);
    
    os << std::setw(m_width) << gpos[0]
       << std::setw(m_width) << gpos[1]
       << std::setw(m_width) << gpos[2]
       << std::setw(m_width) << conf.special().distancefield.distance[j]
       << std::endl;
    
  }
  os << "END" << std::endl;
}

void io::Out_Configuration::_print_angle_restraints(
        configuration::Configuration const &conf,
        topology::Topology const & topo,
        std::ostream &os) {
  DEBUG(10, "angle restraints");

  std::vector<double>::const_iterator ene_it = conf.special().angleres.energy.begin();
  std::vector<double>::const_iterator d_it = conf.special().angleres.d.begin(), 
                                      d_to = conf.special().angleres.d.end();

  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_angle_restraint_precision);
  
  if (conf.special().angleres.d.size() > 0) {
    os << "ANGRESDATA" << std::endl;
    os << std::setw(m_width) << conf.special().angleres.d.size() << "\n";
    int i = 0;
    for (i = 1; d_it != d_to; ++d_it, ++ene_it, ++i) {
       double theta = *d_it * 180 /math::Pi;
       os << std::setw(m_width) << theta << " "
       << std::setw(m_width) << *ene_it;
       os << std::endl;
    }
    os << "END" << std::endl;
  }
  
  std::vector<double>::const_iterator pene_it = conf.special().pertangleres.energy.begin();
  std::vector<double>::const_iterator pd_it = conf.special().pertangleres.d.begin(), 
                                      pd_to = conf.special().pertangleres.d.end();
  
  if (conf.special().pertangleres.d.size() > 0) {
    os << "PERTANGRESDATA" << std::endl;
    os << conf.special().pertangleres.d.size() << "\n";
    int i = 0;
    for (i = 1; pd_it != pd_to; ++pd_it, ++pene_it, ++i) {
       double theta = *pd_it * 180 /math::Pi;
       os << std::setw(m_width) << theta
       << std::setw(m_width) << *pene_it;
       os << std::endl;
    }
    os << "END" << std::endl;
  }
}

void io::Out_Configuration::_print_dihedral_restraints(
        configuration::Configuration const &conf,
        topology::Topology const & topo,
        std::ostream &os) {
  DEBUG(10, "dihedral restraints");

  std::vector<double>::const_iterator ene_it = conf.special().dihedralres.energy.begin();
  std::vector<double>::const_iterator d_it = conf.special().dihedralres.d.begin(), 
                                      d_to = conf.special().dihedralres.d.end();

  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_dihedral_restraint_precision);
  
  if (conf.special().dihedralres.d.size() > 0) {
    os << "DIHRESDATA" << std::endl;
    os << std::setw(m_width) 
       << conf.special().dihedralres.d.size() << "\n";
    int i = 0;
    for (i = 1; d_it != d_to; ++d_it, ++ene_it, ++i) {
       double phi = *d_it * 360 /(2 * math::Pi);
       os << std::setw(m_width) << phi
       << std::setw(m_width) << *ene_it;
       os << std::endl;
    }
    os << "END" << std::endl;
  }
  
  std::vector<double>::const_iterator pene_it = conf.special().pertdihedralres.energy.begin();
  std::vector<double>::const_iterator pd_it = conf.special().pertdihedralres.d.begin(), 
                                      pd_to = conf.special().pertdihedralres.d.end();
  
  if (conf.special().pertdihedralres.d.size() > 0) {
    os << "PERTDIHRESDATA" << std::endl;
    os << conf.special().pertdihedralres.d.size() << "\n";
    int i = 0;
    for (i = 1; pd_it != pd_to; ++pd_it, ++pene_it, ++i) {
       double phi = *pd_it * 360 /(2 * math::Pi);
       os << std::setw(m_width) << phi
       << std::setw(m_width) << *pene_it;
       os << std::endl;
    }
    os << "END" << std::endl;
  }
}

void io::Out_Configuration::_print_order_parameter_restraints(
        configuration::Configuration const &conf,
        topology::Topology const & topo,
        std::ostream &os) {
  DEBUG(10, "order parameter restraints");

  std::vector<double>::const_iterator S2_avg_it = conf.special().orderparamres.S2_avg.begin(),
          S2_avg_to = conf.special().orderparamres.S2_avg.end();
  std::vector<double>::const_iterator e_it = conf.special().orderparamres.energy.begin();

  os << "OPARAMRESDATA" << std::endl;
  for (int i = 1; S2_avg_it != S2_avg_to; ++S2_avg_it, ++e_it, ++i) {
    os << std::setw(6) << i;
    os.precision(m_precision);
    os.setf(std::ios::fixed, std::ios::floatfield);
    os << std::setw(m_width) << *S2_avg_it;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(m_distance_restraint_precision);
    os << std::setw(m_width) << *e_it << std::endl;
  }
  os << "END" << std::endl;
}

void io::Out_Configuration::_print_order_parameter_restraint_averages(
        configuration::Configuration const & conf,
        topology::Topology const & topo,
        std::ostream & os) {
  DEBUG(10, "order parameter restraint averages");

  std::vector<math::Matrix>::const_iterator it = conf.special().orderparamres.Q_avg.begin(),
          to = conf.special().orderparamres.Q_avg.end();
  std::vector<double>::const_iterator d_it = conf.special().orderparamres.D_avg.begin();

  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_distance_restraint_precision); // use a lower precision due to scientific formats

  os << "ORDERPARAMRESEXPAVE" << std::endl;
  int l = 0;
  for (l = 0; it != to; ++it, ++d_it) {
    for(unsigned int i = 0; i < 3; ++i) {
      for(unsigned int j = 0; j < 3; ++j) {
        os << std::setw(m_width) << std::right << (*it)(i,j);
        if (++l % 5 == 0)
          os << std::endl;
      } 
    }
    os << std::setw(m_width) << std::right << *d_it;
    if (++l % 5 == 0)
      os << std::endl;
  }

  os << "END" << std::endl;
}

void io::Out_Configuration::_print_order_parameter_restraint_average_window(
        configuration::Configuration const & conf,
        topology::Topology const & topo,
        std::ostream & os) {
  DEBUG(10, "order parameter restraint averages");

  std::vector<std::list<math::Matrix> >::const_iterator Qwin_it = conf.special().orderparamres.Q_winavg.begin(),
          Qwin_to = conf.special().orderparamres.Q_winavg.end();
  std::vector<std::list<double> >::const_iterator Dwin_it = conf.special().orderparamres.D_winavg.begin();

  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_distance_restraint_precision); // use a lower precision due to scientific formats

  os << "ORDERPARAMRESWINAVE" << std::endl;

  for (unsigned int l = 0; Qwin_it != Qwin_to; ++Qwin_it, ++Dwin_it) {
    std::list<math::Matrix>::const_iterator it = Qwin_it->begin(),
      to = Qwin_it->end();
    std::list<double>::const_iterator D_it = Dwin_it->begin();
    // loop over window
    for(; it != to; ++it, ++D_it) {
      for(unsigned int i = 0; i < 3; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
          os << std::setw(m_width) << std::right << (*it)(i,j);
          if (++l % 5 == 0)
            os << std::endl;
        } 
      }
      os << std::setw(m_width) << std::right << *D_it;
      if (++l % 5 == 0)
        os << std::endl;
    }
  }

  os << "END" << std::endl;
}

void io::Out_Configuration::_print_rdc_values(simulation::Parameter const &param,
                             configuration::Configuration const &conf,
                             topology::Topology const &topo,
                             std::ostream &os,
                             bool formatted) {
  os << "RDCVALUES" << std::endl;
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);
  for (unsigned int i=0; i<conf.special().rdc.size(); ++i){
    for (unsigned int j=0; j<conf.special().rdc[i].curr.size(); ++j){
      os << std::setw(m_width) << conf.special().rdc[i].curr[j] << std::endl;
    }
  }
  os << "END" << std::endl;
}

void io::Out_Configuration::_print_rdc_averages(simulation::Parameter const &param,
                             configuration::Configuration const &conf,
                             topology::Topology const &topo,
                             std::ostream &os,
                             bool formatted) {
  os << "RDCAVERAGES" << std::endl;
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);
  for (unsigned int i=0; i<conf.special().rdc.size(); ++i){
    for (unsigned int j=0; j<conf.special().rdc[i].av.size(); ++j){
      os << std::setw(m_width) << conf.special().rdc[i].av[j] << std::endl;
    }
  }
  os << "END" << std::endl;
}

void io::Out_Configuration::_print_rdc_representation(simulation::Parameter const &param,
                             configuration::Configuration const &conf,
                             topology::Topology const &topo,
                             std::ostream &os,
                             bool formatted) {
  switch(param.rdc.type){
    case(simulation::rdc_mf): {
      os << "RDCMF" << std::endl;
      switch(param.rdc.method) {
        case simulation::rdc_em: {
          if (formatted) {
            os << "#   x              y              z" << std::endl;
          }
          os.setf(std::ios::fixed, std::ios::floatfield);
          os.precision(m_precision);
          for (unsigned int i=0; i<conf.special().rdc.size(); ++i){
            for (unsigned int j=0; j<conf.special().rdc[i].MFpoint.size(); ++j) {
              os << std::setw(m_width) << conf.special().rdc[i].MFpoint[j](0)
                 << std::setw(m_width) << conf.special().rdc[i].MFpoint[j](1)
                 << std::setw(m_width) << conf.special().rdc[i].MFpoint[j](2) << std::endl;
            }
          }
          
          break;
        }
        case simulation::rdc_md:
        case simulation::rdc_sd: {
          if (formatted) {
            os << "#   x              y              z              vx             vy             vz             mass" << std::endl;
          }
          os.setf(std::ios::fixed, std::ios::floatfield);
          os.precision(m_precision);
          for (unsigned int i=0; i<conf.special().rdc.size(); ++i){
            for (unsigned int j=0; j<conf.special().rdc[i].MFpoint.size(); ++j) {
              os << std::setw(m_width) << conf.special().rdc[i].MFpoint[j](0)
                 << std::setw(m_width) << conf.special().rdc[i].MFpoint[j](1)
                 << std::setw(m_width) << conf.special().rdc[i].MFpoint[j](2)
                 << std::setw(m_width) << conf.special().rdc[i].MFpointVel[j](0)
                 << std::setw(m_width) << conf.special().rdc[i].MFpointVel[j](1)
                 << std::setw(m_width) << conf.special().rdc[i].MFpointVel[j](2)
                 << std::setw(m_width) << conf.special().rdc[i].MFpointMass[j] << std::endl;
            }
          }
          break;
        }
      }
      os << "END" << std::endl;
      break;
    }
    case(simulation::rdc_t): {
      os << "RDCT\n";
      switch(param.rdc.method) {
        case simulation::rdc_em: {
          const int n_ah = 5;
          if (formatted) {
            os << "#   Axx            Ayy            Axy            Axz            Ayz\n";
          }
          os.setf(std::ios::fixed, std::ios::floatfield);
          os.precision(m_precision);
          for (unsigned int i=0; i<conf.special().rdc.size(); ++i){
            for (unsigned int j=0; j<n_ah; ++j) {
              os << std::setw(m_width) << conf.special().rdc[i].Tensor[j];
            }
            os << std::endl;
          }
          break;
        }
        case simulation::rdc_md:
        case simulation::rdc_sd: {
          const int n_ah = 5;
          if (formatted) {
            os <<      "#   Axx            vel1           mass1"
              << "          Ayy            vel2           mass2"
              << "          Axy            vel3           mass3"
              << "          Axz            vel4           mass4"
              << "          Ayz            vel5           mass5\n";
          }
          os.setf(std::ios::fixed, std::ios::floatfield);
          os.precision(m_precision);
          for (unsigned int i=0; i<conf.special().rdc.size(); ++i){
            for (unsigned int j=0; j<n_ah; ++j) {
              os << std::setw(m_width) << conf.special().rdc[i].Tensor[j]
                 << std::setw(m_width) << conf.special().rdc[i].TensorVel[j]
                 << std::setw(m_width) << conf.special().rdc[i].TensorMass[j];
            }
            os << std::endl;
          }
          break;
        }
      }
      os << "END" << std::endl;
      break;
    }
    case(simulation::rdc_sh): {
      os << "RDCSH" << std::endl;
      switch(param.rdc.method) {
        case simulation::rdc_em: {
          const int n_clm = 5;
          if (formatted) {
            os << "#   c2,-2          c2,-1          c2,0           c2,1           c2,2\n";
          }
          os.setf(std::ios::fixed, std::ios::floatfield);
          os.precision(m_precision);
          for (unsigned int i=0; i<conf.special().rdc.size(); ++i){
            for (unsigned int j=0; j<n_clm; ++j) {
              os << std::setw(m_width) << conf.special().rdc[i].clm[j];
            }
            os << std::endl;
          }
          break;
        }
        case simulation::rdc_md:
        case simulation::rdc_sd: {
          const int n_clm = 5;
          if (formatted) {
            os <<       "#   c2,-2          vel1           mass1"
               << "          c2,-1          vel2           mass2"
               << "          c2,0           vel3           mass3"
               << "          c2,1           vel4           mass4"
               << "          c2,2           vel5           mass5\n";
          }
          os.setf(std::ios::fixed, std::ios::floatfield);
          os.precision(m_precision);
          for (unsigned int i=0; i<conf.special().rdc.size(); ++i){
            for (unsigned int j=0; j<n_clm; ++j) {
              os << std::setw(m_width) << conf.special().rdc[i].clm[j]
                 << std::setw(m_width) << conf.special().rdc[i].clmVel[j]
                 << std::setw(m_width) << conf.special().rdc[i].clmMass[j];
            }
            os << std::endl;
          }
          break;
        }
      }
      os << "END" << std::endl;
      break;
    }
  }
}


void io::Out_Configuration::_print_rdc_stochastic_integrals(simulation::Parameter const &param,
                             configuration::Configuration const &conf,
                             topology::Topology const &topo,
                             std::ostream &os,
                             bool formatted) {
  os << "RDCSTOCHINT" << std::endl;
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(m_precision);
  for (unsigned int i=0; i<conf.special().rdc.size(); ++i){
    switch(param.rdc.type){
      case(simulation::rdc_mf): {
        for (unsigned int j=0; j<conf.special().rdc[i].stochastic_integral_mf.size(); ++j){
          os << std::setw(m_width) << conf.special().rdc[i].stochastic_integral_mf[j][0];
          os << std::setw(m_width) << conf.special().rdc[i].stochastic_integral_mf[j][1];
          os << std::setw(m_width) << conf.special().rdc[i].stochastic_integral_mf[j][2] << std::endl;
        }
        break;
      }
      case(simulation::rdc_t): {
        for (unsigned int j=0; j<conf.special().rdc[i].stochastic_integral_t.size(); ++j){
          os << std::setw(m_width) << conf.special().rdc[i].stochastic_integral_t[j] << std::endl;
        }
        break;
      }
      case(simulation::rdc_sh): {
        for (unsigned int j=0; j<conf.special().rdc[i].stochastic_integral_sh.size(); ++j){
          os << std::setw(m_width) << conf.special().rdc[i].stochastic_integral_sh[j] << std::endl;
        }
        break;
      }
    }
  }
  os << "END" << std::endl;
}


void io::Out_Configuration::_print_dihangle_trans(
        configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  DEBUG(10, "dihedral angle transitions");

  std::vector<double>::const_iterator it = conf.special().dihangle_trans.dihedral_angle_minimum.begin(),
          to = conf.special().dihangle_trans.dihedral_angle_minimum.end();

  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(1);
  unsigned int atom_i = 0, atom_j = 0, atom_k = 0, atom_l = 0;

  os << "D-A-T" << std::endl;
  os << "# Dih. No.    Resid  Atoms                                        Old min. -> New min." << std::endl;
  int i = 0;
  for (i = 1; it != to; ++it, ++i) {
    if (conf.special().dihangle_trans.old_minimum[i] > math::epsilon) {
      atom_i = conf.special().dihangle_trans.i[i];
      atom_j = conf.special().dihangle_trans.j[i];
      atom_k = conf.special().dihangle_trans.k[i];
      atom_l = conf.special().dihangle_trans.l[i];
      os << i << std::setw(14) << conf.special().dihangle_trans.resid[i] + 1
         << std::setw(5) << topo.residue_names()[topo.solute().atom(atom_i).residue_nr]
         << std::setw(4) << topo.solute().atom(atom_i).name << " - "
         << topo.solute().atom(atom_j).name << " - "
         << topo.solute().atom(atom_k).name << " - " << topo.solute().atom(atom_l).name
         << std::setw(5) << atom_i + 1 << " - " << atom_j + 1 << " - "
         << atom_k + 1 << " - " << atom_l + 1
         << std::setw(m_width) << 180.0 * conf.special().dihangle_trans.old_minimum[i] / math::Pi << " -> "
         << 180.0 * conf.special().dihangle_trans.dihedral_angle_minimum[i] / math::Pi
         << std::endl;
    }
  }

  os << "END" << std::endl;
}

void io::Out_Configuration::_print_pscale_jrest(configuration::Configuration const &conf,
        topology::Topology const &topo,
        std::ostream &os) {
  DEBUG(10, "JVALUEPERSCALE data");

  std::vector<topology::jvalue_restraint_struct>::const_iterator
  jval_it = topo.jvalue_restraints().begin(),
          jval_to = topo.jvalue_restraints().end();

  os << "JVALUEPERSCALE\n";
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  for (int i = 0; jval_it != jval_to; ++jval_it, ++i) {
    os << std::setw(10) << conf.special().pscale.scaling[i]
            << std::setw(15) << conf.special().pscale.t[i]
            << "\n";
  }
  os << "END\n";
}

static void _print_energyred_helper(std::ostream & os, configuration::Energy const &e) {

  const int numenergygroups = unsigned(e.bond_energy.size());
  const int numbaths = unsigned(e.kinetic_energy.size());
  // const int energy_group_size = numenergygroups * (numenergygroups + 1) /2;

  os << "# totals\n";
  
  os << std::setw(18) << e.total << "\n"// 1
          << std::setw(18) << e.kinetic_total << "\n" // 2
          << std::setw(18) << e.potential_total << "\n" // 3
          << std::setw(18) << e.bonded_total << "\n" // 4
          << std::setw(18) << e.bond_total << "\n" // 5
          << std::setw(18) << e.angle_total << "\n" // 6
          << std::setw(18) << e.improper_total << "\n" // 7
          << std::setw(18) << e.dihedral_total << "\n" // 8
          << std::setw(18) << e.crossdihedral_total << "\n" // 9
          << std::setw(18) << e.nonbonded_total << "\n" // 10
          << std::setw(18) << e.lj_total << "\n" // 11
          << std::setw(18) << e.crf_total << "\n" // 12
          << std::setw(18) << e.ls_total << "\n" // 13
          << std::setw(18) << e.ls_pair_total << "\n" // 14
          << std::setw(18) << e.ls_realspace_total << "\n" // 15
          << std::setw(18) << e.ls_kspace_total << "\n" // 16
          << std::setw(18) << e.ls_a_term_total << "\n" // 17
          << std::setw(18) << e.ls_self_total << "\n" // 18
          << std::setw(18) << e.ls_surface_total << "\n" // 19
          << std::setw(18) << e.self_total << "\n" // 20
          << std::setw(18) << e.special_total << "\n" // 21
          << std::setw(18) << e.sasa_total << "\n" // 22
          << std::setw(18) << e.sasa_volume_total << "\n" // 23
          << std::setw(18) << e.constraints_total << "\n" // 24
          << std::setw(18) << e.distanceres_total << "\n" // 25
          << std::setw(18) << e.disfieldres_total << "\n" // 26
          << std::setw(18) << e.dihrest_total << "\n" // 27
          << std::setw(18) << e.posrest_total << "\n" // 28
          << std::setw(18) << e.jvalue_total << "\n" // 29
          << std::setw(18) << e.xray_total << "\n" // 30
          << std::setw(18) << e.leus_total << "\n" // 31
          << std::setw(18) << e.oparam_total << "\n" // 32
          << std::setw(18) << e.symrest_total << "\n" // 33
          << std::setw(18) << e.eds_vmix << "\n" // 34
          << std::setw(18) << e.eds_vr << "\n" // 35
          << std::setw(18) << e.eds_emax << "\n" // 36
          << std::setw(18) << e.eds_emin << "\n" // 37
          << std::setw(18) << e.eds_globmin << "\n" // 38
          << std::setw(18) << e.eds_globminfluc << "\n" // 39
          << std::setw(18) << e.entropy_term << "\n" // 40
          << std::setw(18) << e.qm_total << "\n" // 41
          << std::setw(18) << e.bsleus_total << "\n" // 42
          << std::setw(18) << e.rdc_total << "\n" // 43
          << std::setw(18) << e.angrest_total << "\n" // 44
          << std::setw(18) << e.nn_valid << "\n"; // 45

  os << "# baths\n";
  os << numbaths << "\n";

  for (int i = 0; i < numbaths; ++i) {
    os << std::setw(18) << e.kinetic_energy[i]
            << std::setw(18) << e.com_kinetic_energy[i]
            << std::setw(18) << e.ir_kinetic_energy[i] << "\n";
  }

  os << "# bonded\n";
  os << numenergygroups << "\n";
  for (int i = 0; i < numenergygroups; i++) {
    os << std::setw(18) << e.bond_energy[i]
            << std::setw(18) << e.angle_energy[i]
            << std::setw(18) << e.improper_energy[i]
            << std::setw(18) << e.dihedral_energy[i]
            << std::setw(18) << e.crossdihedral_energy[i] << "\n";
  }
  
  os << "# nonbonded\n";
  for (int i = 0; i < numenergygroups; i++) {
    for (int j = i; j < numenergygroups; j++) {
      os << std::setw(18) << e.lj_energy[j][i]
              << std::setw(18) << e.crf_energy[j][i]
      //      << std::setw(18) << e.ls_real_energy[j][i]
      // currently only one energy group is permitted for LS calculations. 
      // Therefore the total LS energy is written out.        
      // As soon as multiple energy groups are possible this has to be revised      
              << std::setw(18) << e.ls_total
              << std::setw(18) << e.ls_k_energy[j][i] << "\n";
    }
  }

  os << "# special\n";
  for (int i = 0; i < numenergygroups; i++) {
    os << std::setw(18) << e.constraints_energy[i]
            << std::setw(18) << e.posrest_energy[i]
            << std::setw(18) << e.distanceres_energy[i] // disres
            << std::setw(18) << e.disfieldres_energy[i] // disfieldres
            << std::setw(18) << e.angrest_energy[i]// angle res
            << std::setw(18) << e.dihrest_energy[i] // dihedral res
            << std::setw(18) << e.sasa_energy[i]
            << std::setw(18) << e.sasa_volume_energy[i]
            << std::setw(18) << 0.0 // jval
            << std::setw(18) << e.rdc_energy[i] // rdc
            << std::setw(18) << 0.0 // local elevation
            << std::setw(18) << 0.0 << "\n"; //path integral
  }

  // eds energy of end states
  os << "# eds\n";
  os << "# numstates\n";
  const unsigned int numstates = e.eds_vi.size();
  os << numstates << "\n";
  os << std::setw(18) << "# total"
          << std::setw(18) << "nonbonded"
          << std::setw(18) << "special"
          << std::setw(18) << "offset\n";
  for (unsigned i = 0; i < e.eds_vi.size(); i++) {
    os << std::setw(18) << e.eds_vi[i]
            << std::setw(18) << e.eds_vi[i] - e.eds_vi_special[i]
            << std::setw(18) << e.eds_vi_special[i] 
            << std::setw(18) << e.eds_eir[i] << "\n";
  }

  // write eds energies (vr,{V_i}) here

  // ANITA
  // write precalculate energies for all lambda values
  DEBUG(5,"ANITA writing precalc data");
  DEBUG(5,"ANITA nr lambdas: " << e.A_lj_total.size());
  os << "# precalclam\n";
  os << "# nr_lambdas\n";
  const unsigned int nr_lambdas = e.A_lj_total.size();
  os << nr_lambdas << "\n";
  os << std::setw(18) << "# A_e_lj"
     << std::setw(18) << "B_e_lj"
     << std::setw(18) << "A_e_crf" 
     << std::setw(18) << "B_e_crf"
     << std::setw(18) << "AB_kinetic"
     << std::setw(18) << "AB_bond" 
     << std::setw(18) << "AB_angle"
     << std::setw(18) << "AB_improper"
     // special interactions - Betty
     << std::setw(18) << "AB_disres"
     << std::setw(18) << "AB_angres"
     << std::setw(18) << "AB_dihres"
     << std::setw(18) << "AB_disfld\n";
     /*<< std::setw(18) << /"AB_dihedral\n";*/
  
  for (unsigned i = 0; i < nr_lambdas; i++) {
    os << std::setw(18) << e.A_lj_total[i]
       << std::setw(18) << e.B_lj_total[i]
       << std::setw(18) << e.A_crf_total[i]
       << std::setw(18) << e.B_crf_total[i]
       << std::setw(18) << e.AB_kinetic[i]
       << std::setw(18) << e.AB_bond[i]
       << std::setw(18) << e.AB_angle[i]
       << std::setw(18) << e.AB_improper[i]
       // special interations - Betty
       << std::setw(18) << e.AB_disres[i]
       << std::setw(18) << e.AB_angres[i]
       << std::setw(18) << e.AB_dihres[i]
       << std::setw(18) << e.AB_disfld[i]
       //<< std::setw(18) << e.AB_dihedral[i]
       << "\n";
  }
  os << "# ABdih\n";
  os << std::setw(18) << "# A_dihedral"
     << std::setw(18) << "B_dihedral\n";
  os << std::setw(18) << e.A_dihedral
     << std::setw(18) << e.B_dihedral << "\n";
 // ANITA
}

static void _print_volumepressurered_helper(std::ostream &os,
        double mass,
        double const & phi,
        double const & theta,
        double const & psi,
        simulation::Multibath const & m,
        std::vector<double> const & s,
        configuration::Energy const & e,
        math::Box const & b,
        math::boundary_enum t,
        math::Matrix const & p,
        math::Matrix const & v,
        math::Matrix const & k) {
  const int numbaths = int(m.size());

  os << "# mass\n";
  os << std::setw(18) << mass << "\n";

  os << "# temperature\n";
  os << numbaths << "\n";

  for (int i = 0; i < numbaths; ++i) {
    if (m[i].dof)
      os << std::setw(18) << 2 * e.kinetic_energy[i] / math::k_Boltzmann / m[i].dof;
    else
      os << std::setw(18) << 0.0;
    if (m[i].com_dof)
      os << std::setw(18) << 2 * e.com_kinetic_energy[i] / math::k_Boltzmann / m[i].com_dof;
    else
      os << std::setw(18) << 0.0;
    if (m[i].ir_dof)
      os << std::setw(18) << 2 * e.ir_kinetic_energy[i] / math::k_Boltzmann / m[i].ir_dof;
    else
      os << std::setw(18) << 0.0;

    if (s.size())
      os << std::setw(18) << s[i] << "\n";
    else
      os << std::setw(18) << m[i].scale << "\n";

  }

  os << "# volume\n";
  os << std::setw(18) << math::volume(b, t) << "\n";
  //rotate the volume and pressure tensors into the original Cartesian Coordinate system
  math::Matrixl Rmat(math::rmat(phi, theta, psi));
  math::Box b2(math::product(Rmat, b));
  os << std::setw(18) << b2(0)(0) << std::setw(18) << b2(0)(1) << std::setw(18) << b2(0)(2) << "\n"
          << std::setw(18) << b2(1)(0) << std::setw(18) << b2(1)(1) << std::setw(18) << b2(1)(2) << "\n"
          << std::setw(18) << b2(2)(0) << std::setw(18) << b2(2)(1) << std::setw(18) << b2(2)(2) << "\n";

  os << "# pressure\n";
  math::Matrixl auxp(math::product(Rmat, p));
  os << std::setw(18) << (auxp(0, 0) + auxp(1, 1) + auxp(2, 2)) / 3.0 << "\n";

  math::Matrixl auxv(math::product(Rmat, v));
  // the virial is stored internally as just the outer product of positions and forces
  // so without the -0.5 prefactor
  auxv *= -0.5;

  os << std::setw(18) << (auxv(0, 0) + auxv(1, 1) + auxv(2, 2)) / 3.0 << "\n";
  math::Matrixl auxk(math::product(Rmat, k));
  os << std::setw(18) << (auxk(0, 0) + auxk(1, 1) + auxk(2, 2)) / 3.0 << "\n";

  os << std::setw(18) << auxp(0, 0) << std::setw(18) << auxp(0, 1) << std::setw(18) << auxp(0, 2) << "\n"
          << std::setw(18) << auxp(1, 0) << std::setw(18) << auxp(1, 1) << std::setw(18) << auxp(1, 2) << "\n"
          << std::setw(18) << auxp(2, 0) << std::setw(18) << auxp(2, 1) << std::setw(18) << auxp(2, 2) << "\n";

  os << std::setw(18) << auxv(0, 0) << std::setw(18) << auxv(0, 1) << std::setw(18) << auxv(0, 2) << "\n"
          << std::setw(18) << auxv(1, 0) << std::setw(18) << auxv(1, 1) << std::setw(18) << auxv(1, 2) << "\n"
          << std::setw(18) << auxv(2, 0) << std::setw(18) << auxv(2, 1) << std::setw(18) << auxv(2, 2) << "\n";

  os << std::setw(18) << k(0, 0) << std::setw(18) << k(0, 1) << std::setw(18) << k(0, 2) << "\n"
          << std::setw(18) << k(1, 0) << std::setw(18) << k(1, 1) << std::setw(18) << k(1, 2) << "\n"
          << std::setw(18) << k(2, 0) << std::setw(18) << k(2, 1) << std::setw(18) << k(2, 2) << "\n";

}

void io::Out_Configuration
::_print_position_restraints(simulation::Simulation const & sim,
        topology::Topology const &topo,
        configuration::Configuration const &conf,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "REFPOSITION\n";

  topology::Solute const &solute = topo.solute();
  std::vector<std::string> const &residue_name = topo.residue_names();
  topology::Solvent const &solvent = topo.solvent(0);

  const math::VArray & ref = conf.special().reference_positions;
  math::Matrixl Rmat(math::rmat(conf.current().phi,
          conf.current().theta, conf.current().psi));
  if (conf.boundary_type == math::truncoct) {
    Rmat = math::product(Rmat, math::truncoct_triclinic_rotmat(false));
  }
  const math::SArray & b = conf.special().bfactors;

  for (unsigned int i = 0; i < topo.num_atoms(); ++i) {
    math::Vec r(math::product(Rmat, ref(i)));
    if (i < topo.num_solute_atoms()) {
      os << std::setw(5) << solute.atom(i).residue_nr + 1 << " "
              << std::setw(5) << std::left << residue_name[solute.atom(i).residue_nr] << " "
              << std::setw(6) << std::left << solute.atom(i).name << std::right
              << std::setw(6) << i + 1
              << std::setw(m_width) << r(0)
              << std::setw(m_width) << r(1)
              << std::setw(m_width) << r(2)
              << "\n";
    } else { // just writing out dummy values for first 17 chars
      os << std::setw(5) << "0" << " "
              << std::setw(5) << std::left
              << "SOLV" << " "
              << std::setw(6) << std::left << solvent.atom((i - topo.num_solute_atoms())
              % solvent.atoms().size()).name << std::right
              << std::setw(6) << i + 1
              << std::setw(m_width) << r(0)
              << std::setw(m_width) << r(1)
              << std::setw(m_width) << r(2)
              << "\n";
    }
  }

  os << "END\n";

  if (sim.param().posrest.posrest == simulation::posrest_bfactor) {
    os << "BFACTOR\n";

    for (unsigned int i = 0; i < topo.num_atoms(); ++i) {
      if (i < topo.num_solute_atoms()) {
        os << std::setw(5) << solute.atom(i).residue_nr + 1 << " "
                << std::setw(5) << std::left << residue_name[solute.atom(i).residue_nr] << " "
                << std::setw(6) << std::left << solute.atom(i).name << std::right
                << std::setw(6) << i + 1
                << std::setw(m_width) << b(i)
                << "\n";
      } else {
        os << std::setw(5) << "0" << " "
                << std::setw(5) << std::left
                << "SOLV" << " "
                << std::setw(6) << std::left << solvent.atom((i - topo.num_solute_atoms())
              % solvent.atoms().size()).name << std::right
                << std::setw(6) << i + 1
                << std::setw(m_width) << b(i)
                << "\n";
      }
    }

    os << "END\n";
  }

}

void io::Out_Configuration::
_print_nose_hoover_chain_variables(const simulation::Multibath & multibath,
        std::ostream & os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "NHCVARIABLES\n";
  for (simulation::Multibath::const_iterator it = multibath.begin(), to = multibath.end();
          it != to; ++it) {
    const std::vector<double> & zeta = it->zeta;
    for (std::vector<double>::const_iterator z_it = zeta.begin(), z_to = zeta.end();
            z_it != z_to; ++z_it) {
      os << std::setw(m_width) << *z_it;
    }
    os << "\n";
  }
  os << "END\n";
}

void io::Out_Configuration::
_print_rottrans(configuration::Configuration const &conf,
        simulation::Simulation const &sim,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  os << "ROTTRANSREFPOS\n";
  os << "# translational matrix\n";
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      os << std::setw(m_width) << conf.special().rottrans_constr.theta_inv_trans(i, j);
    }
    os << "\n";
  }
  os << "# rotational matrix\n";
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      os << std::setw(m_width) << conf.special().rottrans_constr.theta_inv_rot(i, j);
    }
    os << "\n";
  }
  os << "# reference positions\n";
  unsigned int last = sim.param().rottrans.last;
  const math::VArray & refpos = conf.special().rottrans_constr.pos;
  for (unsigned int i = 0; i < last; ++i) {
    os << std::setw(m_width) << refpos(i)(0)
            << std::setw(m_width) << refpos(i)(1)
            << std::setw(m_width) << refpos(i)(2) << "\n";
  }
  os << "END\n";
}

void io::Out_Configuration::
_print_umbrellas(configuration::Configuration const & conf, std::ostream & os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  const std::vector<util::Umbrella> & umb = conf.special().umbrellas;
  os << "LEUSBIAS\n# NUMUMBRELLAS\n"
     << umb.size()
     << "\n";
  for(unsigned int i = 0; i < umb.size(); ++i) {
    os << "# NLEPID NDIM CLES\n";
    os << std::setw(10) << umb[i].id << std::setw(10) << umb[i].dim() << std::setw(15)
            << umb[i].force_constant << "\n";
    os << "# VARTYPE(N) NTLEFU(N) WLES(N) RLES(N) NGRID(N) GRIDMIN(N) GRIDMAX(N)\n";
    for(unsigned int j = 0; j < umb[i].dim(); ++j) {
      os << std::setw(10) << umb[i].variable_type[j]
              << std::setw(10) << umb[i].functional_form[j]
              << std::setw(15) << umb[i].width[j]
              << std::setw(15) << umb[i].cutoff[j]
              << std::setw(10) << umb[i].num_grid_points[j]
              << std::setw(15) << umb[i].grid_min[j]
              << std::setw(15) << umb[i].grid_max[j] << "\n";
    }
    os << "# NCONLE\n";
    os << std::setw(10) << umb[i].configurations.size() << "\n";
    os << "# NVISLE ICONF(1..NDIM)\n";
    std::map<util::Umbrella::leus_conf,util::Umbrella_Weight*>::const_iterator conf_it =
    umb[i].configurations.begin(), conf_to = umb[i].configurations.end();
    for(; conf_it != conf_to; ++conf_it) {
      os << *(conf_it->second);
      for(unsigned int j = 0; j < umb[i].dim(); j++) {
        // here we have to convert to fortran and add 1 because the grid
        // goes form 1 to NDIM in the output format
        os << std::setw(10) << conf_it->first.pos[j]+1;
      }
      os << "\n";
    }
  }
  os << "END\n";
  os << "LEUSPOS\n# NUMUMBRELLAS\n"
     << umb.size()
     << "\n";
  for(unsigned int i = 0; i < umb.size(); ++i) {
    os << "# NLEPID NDIM\n";
    os << std::setw(10) << umb[i].id << std::setw(10) << umb[i].dim() << "\n";
    os << "# POS(1..NDIM)\n";
    std::vector<util::LE_Coordinate*>::const_iterator coord_it =
    umb[i].coordinates.begin(), coord_to = umb[i].coordinates.end();
    for(; coord_it != coord_to; ++coord_it) {
      for(unsigned int j = 0; j < umb[i].dim(); j++) {
        os << std::setw(10) << (*coord_it)->get_value(umb[i].grid_min[j], umb[i].grid_max[j]);
      }
      os << "\n";
    }
  }
  os << "END\n";
}

/**
 * @section bsleusmem BSLEUSMEM block
 * 
 * Defines the state of the memory and the auxiliary memory as well as the 
 * auxiliary and the reduction counter of the B&S-LEUS scheme.
 * 
 * @verbatim
 BSMEM
#
# The current state of the Memory in the BS&LEUS algorithm.
# NUMSPH:   The Number of Spheres
# NUMSTK:   The Number of Sticks
# AUXMEM:   Is the auxiliary memory in use (BSLEUS:FRED != 1)?
#   0:          no auxiliary memory; set to zero, if needed.
#   1:          yes, there is auxiliary memory
#
# NUMSPH    NUMSTK  AUXMEM
  4         2       1
#
# SPHID:    The ID of the Sphere
# SUBSP:    In which subspace the Sphere/Stick is
# NUMGP:    The number of Grid Points
# MEM:      The Memory
#
# SPHID SUBSP   NUMGP   MEM[1..NUMGP]
  1     1       10      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  2     1       10      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  3     1       10      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  4     1       10      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#
# STKID:    The ID of the Stick
# SUBSP:    In which subspace the Sphere/Stick is
# NUMGP:    The Number of Grid Points
# MEM:      The Memory
#
# STKID SUBSP   NUMGP   MEM[1..NUMGP]
  1     1       20      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  2     1       20      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#
# The Auxiliary Memory & the subspace
#
# NUMSUB:   The number of subspaces
# SUBID:    The ID of the subspace
# AUXC:     The Auxiliary Counter
# REDC:     The Reduction Counter
#
# NUMSUB
  1
# SUBID AUXC    REDC
  1     0       0
#
# SPHID:    The ID of the Sphere
# STKID:    The ID of the Stick
# SUBSP:    In which subspace the Sphere/Stick is
# NUMGP:    The Number of Grid Points
# AUXMEM:   The Memory
#
# SPHID SUBSP   NUMGP   AUXMEM[1..NUMGP]
  1     1       10      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  2     1       10      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  3     1       10      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  4     1       10      0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#
# STKID SUBSP   NUMGP   AUXMEM[1..NUMGP]
  1     1       20      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  2     1       20      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
END
@endverbatim
 * 
 */
void io::Out_Configuration::_print_bsleusmem(const configuration::Configuration &conf, 
                                          std::ostream& os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  const util::BS_Umbrella *umb = &conf.special().bs_umbrella;
  if (umb == 0){
    io::messages.add("There seems to be no umbrella!", "Out_Configuration",
            io::message::error);
  }
  os << "BSLEUSMEM\n"
     << "#\n"
     << "# The current state of the Memory in the BS&LEUS algorithm.\n"
     << "# NUMPOT:   The Number of Potentials\n"
     << "# AUXMEM:   Is the auxiliary memory set? (0: No; 1: Yes)\n"
     << "#\n"
     << "# NUMPOT AUXMEM\n";
  int numPotentials = 0;
  umb->getNumPotentials(numPotentials);
  os << std::setw(3) << numPotentials << std::setw(7) << umb->printAuxMem() << "\n";
  os << "#\n"
          << "# ID:       The ID of the Sphere\n"
          << "# SUBID:    The ID of the Subspace\n"
          << "# NUMGP:    The number of Grid Points\n"
          << "# MEM:      The Memory\n"
          << "#\n"
          << "# ID    SUBSP   NUMGP   MEM[1..NUMGP]\n";
  for (int i = 0; i < numPotentials; i++) {
    std::vector<double> memory;
    unsigned int subid = 0;
    if (!umb->getMemory(i + 1, subid, memory)) {
      std::ostringstream msg;
      msg << "Could not find the memory of potential " << i + 1 << "!\n";
      io::messages.add(msg.str(), "Out_Configuration", io::message::error);
      return;
    }
    int num_gp = memory.size();
    os << std::setw(3) << (i + 1)
            << std::setw(6) << subid + 1
            << std::setw(8) << num_gp;
    std::vector<double>::iterator it = memory.begin(),
            to = memory.end();
    for (; it != to; it++) {
      os << " " << *it;
    }
    os << "\n";
  }
  
  if (umb->printAuxMem()) {
    unsigned int numSubspaces = umb->getNumSubspaces();
    os << "#\n"
            << "# The Auxiliary Memory & the subspaces\n"
            << "#\n"
            << "# NUMSUB:   The number of subspaces\n"
            << "# SUBID:    The ID of the subspace\n"
            << "# AUXC:     The Auxiliary Counter\n"
            << "# REDC:     The Reduction Counter\n"
            << "# AUXMEM:   The Auxiliary Memory\n"
            << "#\n"
            << "# NUMSUB\n  "
            << numSubspaces << "\n"
            << "# SUBID AUXC    REDC\n";
    for (unsigned int i = 0; i < numSubspaces; i++) {
      unsigned int auxc = 0, redc = 0;
      umb->getCounter(i, auxc, redc);
      os << std::setw(3) << i + 1
              << std::setw(6) << auxc
              << std::setw(8) << redc << "\n";
    }

    os << "#\n"
       << "# ID    SUBSP   NUMGP   AUXMEM[1..NUMGP]\n";
    for (int i = 0; i < numPotentials; i++) {
      std::vector<double> memory;
      unsigned int subid = 0;
      if (!umb->getAuxMemory(i + 1, subid, memory)) {
        std::ostringstream msg;
        msg << "Could not find the memory of potential " << i + 1 << "!\n";
        io::messages.add(msg.str(), "Out_Configuration", io::message::error);
        return;
      }
      int num_gp = memory.size();
      os << std::setw(3) << (i + 1)
              << std::setw(6) << subid + 1
              << std::setw(8) << num_gp;
      std::vector<double>::iterator it = memory.begin(),
              to = memory.end();
      for (; it != to; it++) {
        os << " " << *it;
      }
      os << "\n";
    }
  }
  os << "END\n";
}

/**
 * @section bsleuspos BSLEUSPOS block
 * 
 * Defines the position in the B&S-LEUS subspaces.
 * 
 * @verbatim
BLEUSPOS
#
# The current position in the BSLEUS subspaces
#
# NUMSUB:   The Number of Subspaces
# SUBID:    The Number of the Subspace
# NUMDIM:   The number of dimensions of the subspace
# POS:      The position
#
# NUMSUB
  1
#
# SUBID NUMDIM  POS[1..NUMDIM]
  1     2       300 379
END
@endverbatim
 * 
 */
void io::Out_Configuration::_print_bsleuspos(const configuration::Configuration &conf, 
                                          std::ostream& os)
{
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  const util::BS_Umbrella *umb = &conf.special().bs_umbrella;
  if (umb == 0){
    io::messages.add("There seems to be no umbrella!", "Out_Configuration",
            io::message::error);
  }
  os << "BSLEUSPOS\n"
      << "#\n"
      << "# The current position in the BSLEUS subspaces\n"
      << "#\n"
      << "# NUMSUB:   The Number of Subspaces\n"
      << "# SUBID:    The Number of the Subspace\n"
      << "# NUMDIM:   The number of dimensions of the subspace\n"
      << "# POS:      The position\n"
      << "#\n"
      << "# NUMSUB\n"
      << std::setw(3) << umb->getNumSubspaces() << "\n"
      << "#\n"
      << "# SUBID NUMDIM  POS[1..NUMDIM]\n";
  for (unsigned int i = 0; i < umb->getNumSubspaces(); i++){
      os << std::setw(3) << i + 1;
      std::vector<double> position;
      umb->getPosition(i, position);
      os << std::setw(5) << position.size();
      for (unsigned int j = 0; j < position.size(); j++){
          os << " " << position[j];
      }
      os << "\n";
  }
  os << "END\n";
}

void io::Out_Configuration::
_print_bsleus(const configuration::Configuration &conf, std::ostream& os)
{
  os << "BSLEUS\n" 
          << conf.special().bs_umbrella.traj_str() 
          << "END\n";
}

template<math::boundary_enum b>
void io::Out_Configuration::
_print_dipole(simulation::Simulation const & sim,
              topology::Topology const &topo,
              configuration::Configuration const & conf, std::ostream & os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  //topology::Solute const &solute = conf.init(topology::Topology.solute());
  //topology::Solvent const &solvent = conf.init(topology::Topology.solvent(0));



  //Calculate dipole
  

  math::Vec box_dipole_moment(0.0);
  math::Vec box_centre = conf.current().box(0) / 2.0 +
                         conf.current().box(1) / 2.0 +
                         conf.current().box(2) / 2.0;
  math::Periodicity<b> periodicity(conf.current().box);
 

  // Only solute
  if(sim.param().electric.dip_groups == 0) 
  {
    os << "DIPOLE\n#Solute atoms\n";
    for(unsigned int i = 0; i < topo.num_solute_atoms(); ++i) 
    {
      math::Vec r = conf.current().pos(i);
      if (topo.is_polarisable(i)) 
      {
        //offset position
	math::Vec rm=r;
        
	//cos dipol contribution
        box_dipole_moment += topo.coscharge(i) * conf.current().posV(i);
        
	if(sim.param().polarise.cos == 2 && topo.gamma(i)!=0.0)
        {
            math::Vec rij, rik, rim;
            periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_j(i)), rij);
            periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(topo.gamma_k(i)), rik);
            rim=topo.gamma(i)*(rij+rik)/2;
            rm-=rim;
            box_dipole_moment += topo.charge(i) * (rm- box_centre);
        }
        else
        {
            box_dipole_moment += topo.charge(i) * (r - box_centre);
        }
      }
      else
      {
          box_dipole_moment += topo.charge(i) * (r - box_centre);
      }
    }
  }

  // Only solvent
  if(sim.param().electric.dip_groups == 1) {
    os << "DIPOLE\n#Solvent atoms\n";
    for(unsigned int i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i) 
    {
     math::Vec r = conf.current().pos(i);
     if (topo.is_polarisable(i)) 
     {
       //offset position
       math::Vec rm=r;
        
       //cos dipol contribution
       box_dipole_moment += topo.coscharge(i) * conf.current().posV(i);
        
       if(sim.param().polarise.cos == 2 && topo.gamma(i)!=0.0)
       {
	 math::Vec rij, rik, rim;
	 periodicity.nearest_image(conf.current().pos(i),
				   conf.current().pos(topo.gamma_j(i)), rij);
	 periodicity.nearest_image(conf.current().pos(i),
				   conf.current().pos(topo.gamma_k(i)), rik);
	 rim=topo.gamma(i)*(rij+rik)/2;
	 rm-=rim;
	 box_dipole_moment += topo.charge(i) * (rm- box_centre);
       }
       else
       {
	 box_dipole_moment += topo.charge(i) * (r - box_centre);
       }
     }
     else
     {
       box_dipole_moment += topo.charge(i) * (r - box_centre);
     }
    }
  }

  // All atoms in the box
  if(sim.param().electric.dip_groups == 2) {
    os << "DIPOLE\n#All atoms\n";
    for(unsigned int i = 0; i < topo.num_atoms(); ++i) 
    {
      math::Vec r = conf.current().pos(i);
      if (topo.is_polarisable(i)) 
      {
	//offset position
	math::Vec rm=r;
        
	//cos dipol contribution
	box_dipole_moment += topo.coscharge(i) * conf.current().posV(i);
        
	if(sim.param().polarise.cos == 2 && topo.gamma(i)!=0.0)
	{
	  math::Vec rij, rik, rim;
	  periodicity.nearest_image(conf.current().pos(i),
				    conf.current().pos(topo.gamma_j(i)), rij);
	  periodicity.nearest_image(conf.current().pos(i),
				    conf.current().pos(topo.gamma_k(i)), rik);
	  rim=topo.gamma(i)*(rij+rik)/2;
	  rm-=rim;
          box_dipole_moment += topo.charge(i) * (rm- box_centre);
	}
	else
        {
	  box_dipole_moment += topo.charge(i) * (r - box_centre);
	}
      }
      else
      {
	box_dipole_moment += topo.charge(i) * (r - box_centre);
      }
    }
  }

  // Print dipole to special topology

  //os << sim.param().polarise.cos <<"\t" << topo.gamma(0) << "\n";

  os << std::setw(15) << box_dipole_moment(0)
     << std::setw(15) << box_dipole_moment(1)
     << std::setw(15) << box_dipole_moment(2)
     << std::setw(15) << math::volume(conf.current().box, conf.boundary_type)
     << "\n";
  
  os << "END\n";
  
}

void io::Out_Configuration::
_print_current(simulation::Simulation const & sim,
              topology::Topology const &topo,
              configuration::Configuration const & conf, std::ostream & os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  unsigned int ngroups = sim.param().electric.cur_groups;
  int first = 0;
  math::Vec cur_current(0.0);
  //double scale = math::four_pi_eps_i;
  double scale = 1.0;

  os << "CURRENT\n";
  os << "#GROUP        COMPONENTS\n";

  for (unsigned  int i = 0; i < ngroups; ++i){
    
    for (unsigned int j = first; j < sim.param().electric.current_group[i]; ++j){
      math::Vec v = conf.current().vel(j);
      cur_current += scale*topo.charge(j) * v;
    }
    int grp = i+1;
    os << std::setw(10) << grp << std::setw(15) << cur_current(0)
                               << std::setw(15) << cur_current(1)
                               << std::setw(15) << cur_current(2) << "\n";
    cur_current(0.0);
  }
  os << "END\n";
}

void io::Out_Configuration::
_print_adde(simulation::Simulation const & sim,
        topology::Topology const &topo,
        configuration::Configuration const & conf, std::ostream & os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  int not_adde = 0;
  double vhh = conf.special().adde.vhh;

  for (unsigned int i = 0; i < sim.param().multibath.multibath.size(); ++i) {
    if (i != sim.param().addecouple.adc_index()[0].tg)
      not_adde = i;
  }
  double betal = 1 /
          (math::k_Boltzmann * sim.param().multibath.multibath[not_adde].temperature);
  double betah = sim.param().addecouple.adc_index()[0].sv /
          (math::k_Boltzmann *
          sim.param().multibath.multibath[sim.param().addecouple.adc_index()[0].tg].temperature);
  double lnevhl = conf.special().adde.evhl;
  //if(conf.special().adde.evhl==0.0)
  //  lnevhl = 0;
  //else
  //  lnevhl = std::log(conf.special().adde.evhl)+betal*conf.special().adde.vhl0;

  //double scale = math::four_pi_eps_i;

  os << "ADDEREWEIGHTING\n";
  os << "#VHH LNEVHL BETAL BETAH  \n";

  os << std::setw(15) << vhh << " " << std::setw(15) << lnevhl << " " 
          << std::setw(15) << betal << std::setw(15) << betah << "\n";
  os << "END\n";
    
}

void io::Out_Configuration::
_print_aedssearch(configuration::Configuration const &conf,
        simulation::Simulation const &sim,
        std::ostream &os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);
  os << "AEDSSEARCH\n";
  os << std::setw(m_width) << sim.param().eds.emax << "\n";
  os << std::setw(m_width) << sim.param().eds.emin << "\n";
  os << std::setw(m_width) << sim.param().eds.searchemax << "\n";
  os << std::setw(m_width) << sim.param().eds.emaxcounts << "\n";
  os << std::setw(m_width) << sim.param().eds.oldstate << "\n";
  os << std::setw(m_width) << sim.param().eds.fullemin << "\n";
  for (unsigned int i = 0; i < sim.param().eds.numstates; i++) {
    os << std::setw(m_width) << sim.param().eds.eir[i] << " "
      << std::setw(m_width) << sim.param().eds.lnexpde[i] << " "
      << std::setw(m_width) << sim.param().eds.statefren[i] << " "
      << std::setw(m_width) << sim.param().eds.visitedstates[i] << " "
      << std::setw(m_width) << sim.param().eds.visitcounts[i] << " "
      << std::setw(m_width) << sim.param().eds.avgenergy[i] << " "
      << std::setw(m_width) << sim.param().eds.eiravgenergy[i] << " "
      << std::setw(m_width) << sim.param().eds.bigs[i] << " "
      << std::setw(m_width) << sim.param().eds.stdevenergy[i] << "\n";
  }
  os << "END\n";
}

void io::Out_Configuration::
_print_nemd(simulation::Simulation const & sim,
        topology::Topology const &topo,
        configuration::Configuration const & conf, std::ostream & os) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(m_precision);

  
  /*
     * METHOD SELECTION 
     * 
     */

    
    int prop_vs_method = 0;
    if (sim.param().nemd.property == 0){
      if (sim.param().nemd.method == 0)
        prop_vs_method = 0;
      if (sim.param().nemd.method == 1)
        prop_vs_method = 1;
      //enter selection for other methods for this property here
    }
    //enter selection for other properties with respective methods here
  
  
  

  switch (prop_vs_method){
    
    
    case 0:
    {
      double k = 2 * math::Pi / conf.current().box(2)(2);
      double k2 = k*k;
      double lzdiv2pisq = 1/k2;
      //double tau = sim.param().nemd.pertfrq * sim.time_step_size();
      double rho=0.0;

      for (unsigned int i=0; i<topo.num_atoms(); ++i)
        rho += topo.mass(i);

      double volume = conf.current().box(0)(0)*conf.current().box(1)(1)*conf.current().box(2)(2);

      rho /= volume;
      os << "NEMD\n";
      os << "# Amplitude of applied perturbation: " << std::setw(15) << sim.param().nemd.ampbath <<"\n";
      os << "# (Lz/2*pi)^2                      : " << std::setw(15) << lzdiv2pisq << "\n";
      //os << "# tau:             " << std::setw(15) << tau << "\n";
      os << "# rho (mass density)               : " << std::setw(15) << rho << "\n";
      //os << "# k^2 * tau / rho: " << std::setw(15) << k2*tau/rho;

      // print amplitudes
      double vel = 0.0;
      math::Periodicity<math::rectangular> periodicity(conf.current().box);

      for (unsigned int i=0; i<topo.num_atoms(); ++i){
        math::Vec pos = conf.current().pos(i);
        periodicity.put_into_positive_box(pos);
        double zbath = pos(2);
        double a = conf.current().vel(i)(0) * cos(zbath*k);
        vel += a;
      }
      vel *= 2;
      vel /= topo.num_atoms();

      // vel = 2/N * sum_i_N(v_ix * cos(2*pi/L_z))

      os << "\nVELOCITY AMPLITUDE" << std::setw(15) << vel << "\n";


      /*
       * Only printing the grids after stdyaft steps
       *
       */
      if(sim.steps() > sim.param().nemd.stdyaft) {
        os << "#AVERAGE X-COMPONENT OF VELOCITY PER SLAB\n";
        for(unsigned int i = 0; i < conf.special().nemd_conf.stored_data_per_bin.size(); ++i) {
          //per steps
          double averaged = conf.special().nemd_conf.stored_data_per_bin[i] / (sim.steps() - sim.param().nemd.stdyaft);
          //not averaged, but using the variable name
          //double averaged = conf.special().nemd_conf.stored_data_per_bin[i];
          os << std::setw(15) << i + 1 << std::setw(15) << averaged << "\n";
        }
      }
      
      os << "END\n";
      break;
    }
    
    case 1:
    {
      os << "NEMD\n";
      //math::Periodicity<math::rectangular> periodicity(conf.current().box);
      double flux = conf.special().nemd_conf.Px/(2*conf.current().box(0)(0)
                    *conf.current().box(1)(1) * sim.steps() * sim.time_step_size());
      os << "#MOMENTUM Px     FLUX\n" << std::setw(15) << conf.special().nemd_conf.Px
              << std::setw(15) << flux;
      //conf.special().remd_conf.Px = 0;
      os << "\n#AVERAGE X-COMPONENT OF VELOCITY PER SLAB\n";

      //int nslab = 2*sim.param().nemd.slabnum;

      //std::vector<std::vector<unsigned int> > grid(nslab);

        /*for(unsigned int i = 0; i < grid.size(); ++i) {
          grid[i].clear();
        }*/
        
        //Put atom indexes into grid
        /*for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
          math::Vec pos = conf.current().pos(i);
          periodicity.put_into_positive_box(pos);
          double z = pos(2);
          int bin = z / conf.current().box(2)(2) * nslab;
          grid[bin].push_back(i);
        }*/

        /*
        for(unsigned int slice = 0; slice < nslab; ++slice) {
          double vel_i = 0.0;
          int counter = 0;
          for(unsigned int i = 0; i < grid[slice].size(); ++i){
            vel_i += conf.current().vel(grid[slice][i])(0);
            counter++;
          }
          double av_vel=vel_i/counter;
          os <<  std::setw(15) << slice+1 << std::setw(15) << av_vel << "\n";

        }*/
      // CHANGE THIS!!!! This is not considering the relaxation process
      if(sim.steps() > sim.param().nemd.stdyaft) {
      for (unsigned int  i=0; i<conf.special().nemd_conf.stored_data_per_bin.size(); ++i){
        //per steps
        double averaged = conf.special().nemd_conf.stored_data_per_bin[i]/(sim.steps()-sim.param().nemd.stdyaft);
        //not averaged, but using the variable name
        //double averaged = conf.special().nemd_conf.stored_data_per_bin[i];
        os <<  std::setw(15) << i+1 << std::setw(15) << averaged << "\n";
      }
      }

      os << "END\n";
      break;
    }
    
    default: break;
  }
}
