/**
 * @file mg.cc
 * md using multi-graining
 * test program, later to be incorporated into repex
 *
 * @todo document it
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>
#include <interaction/special/external_interaction.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/virtual_grain.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>

#include <io/configuration/out_configuration.h>

#ifdef OMP
#include <omp.h>
#endif

int main(int argc, char *argv[]){

  const double start = util::now();

  util::Known knowns;
  knowns << "topo" << "cg_topo" << "conf" << "cg_conf" << "input" << "cg_input" 
	 << "verb" << "pttopo" << "cg_pttopo"
	 << "trc" << "cg_trc" << "fin" << "cg_fin" << "trv" << "trf" << "trs" << "tre" << "cg_tre" << "trg"
	 << "bae" << "bag" << "posresspec" << "refpos" <<"distrest" << "jval" << "xray" << "sym" << "order" << "lud" << "led"
	 << "anatrj" << "print" << "friction" << "qmmm"
	 << "version";
  
  
  std::string usage;
  util::get_usage(knowns, usage, argv[0]);

  io::Argument args;

  if (args.parse(argc, argv, knowns)){
    std::cerr << usage << std::endl;
    return 1;
  }

  util::print_title(false);
  if (args.count("version") >= 0){
    return 0;
  }
    
  // parse the verbosity flag and set debug levels
  if (util::parse_verbosity(args)){
    std::cerr << "could not parse verbosity argument" << std::endl;
    return 1;
  }

  // create the simulation classes
  topology::Topology topo;
  configuration::Configuration conf;
  algorithm::Algorithm_Sequence md;
  simulation::Simulation sim;

  topology::Topology cg_topo;
  configuration::Configuration cg_conf;
  algorithm::Algorithm_Sequence cg_md;
  simulation::Simulation cg_sim;

  io::Out_Configuration traj(GROMOSXX "\n\tfine grained\n");
  io::Out_Configuration cg_traj(GROMOSXX "\n\tcoarse grained\n");
  
  // add an external interaction
  sim.param().force.external_interaction = 1;
  io::read_input(args, topo, conf, sim,  md);

  interaction::Forcefield * ff = 
    dynamic_cast<interaction::Forcefield *>(md.algorithm("Forcefield"));
  if (ff == NULL){
    std::cout << "Error: no forcefield in MD" << std::endl;
    return 1;
  }
  interaction::External_Interaction * ei = 
    dynamic_cast<interaction::External_Interaction *>(ff->interaction("External"));
  if (ei == NULL){
    std::cout << "Error: no external interaction in forcefield" << std::endl;
    return 1;
  }

  ei->set_coarsegraining(cg_topo, cg_conf, cg_sim);

  traj.title(GROMOSXX "\n\tfine-grained\n" + sim.param().title);
  traj.init(args, sim.param());

  io::argname_conf = "cg_conf";
  io::argname_topo = "cg_topo";
  io::argname_pttopo = "cg_pttopo";
  io::argname_input = "cg_input";
  io::argname_trj = "cg_trj";
  io::argname_fin = "cg_fin";
  io::argname_tre = "cg_tre";
  
  io::read_input(args, cg_topo, cg_conf, cg_sim, cg_md);
  interaction::Forcefield * cg_ff =
    dynamic_cast<interaction::Forcefield *>(cg_md.algorithm("Forcefield"));
  if (cg_ff == NULL){
    std::cout << "Error: no forcefield in cg_MD" << std::endl;
    return 1;
  }
  
  cg_traj.title(GROMOSXX "\n\tcoarse-grained\n" + sim.param().title);
  cg_traj.init(args, cg_sim.param());

  // initialises all algorithms (and therefore also the forcefield)
  md.init(topo, conf, sim);
  cg_md.init(cg_topo, cg_conf, cg_sim);

  std::cout << "\nMESSAGES FROM INITIALISATION\n";
  if (io::messages.display(std::cout) >= io::message::error){
    std::cout << "\nErrors during initialisation!\n" << std::endl;
    return 1;
  }
    
  io::messages.clear();

  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
  std::cout << "\nenter the next level of molecular "
	    << "dynamics simulations\n" << std::endl;
    
  int msteps = sim.param().multistep.steps;
  if (msteps < 1) msteps = 1;
  
  
  std::cout << "==================================================\n"
	    << " MAIN MD LOOP\n"
	    << "==================================================\n"
	    << std::endl;

  int error;
  int percent = 0;
  const double init_time = util::now() - start;
    
  while(int(sim.steps()) < sim.param().step.number_of_steps){  
      
    traj.write(conf, topo, sim, io::reduced);
    cg_traj.write(cg_conf, cg_topo, cg_sim, io::reduced);

    // coarse grained atom positions are based upon
    // real atom positions
    // make sure the lambdas are identical
    cg_sim.param().perturbation.lambda = topo.lambda();
    cg_topo.lambda(topo.lambda());
    cg_topo.lambda(topo.lambda());
    cg_topo.update_for_lambda();

    // if ((sim.steps() % msteps) == 0){
      
    // std::cerr << "----- cg md --------------------" << std::endl;
    // std::cerr << "\tupdate_virtual_pos " << sim.steps() << std::endl;
    util::update_virtual_pos(cg_topo, cg_conf, topo, conf);

    // calculate the cg forces first!
    if ((error = cg_md.run(cg_topo, cg_conf, cg_sim))){

      if (error == E_MINIMUM_REACHED) // ignore this...
	error = 0;
      else{
	io::print_ENERGY(traj.output(), cg_conf.current().energies,
			 cg_topo.energy_groups(),
			 "CGERROR", "CGERR_");
	
	io::print_ENERGY(traj.output(), cg_conf.old().energies, 
			 cg_topo.energy_groups(),
			 "CGOLDERROR", "CGOLDERR_");
	
	std::cout << "\nError during CG MD run!\n" << std::endl;
	break;
      }
    }
    // std::cerr << "----- cg md done ---------------" << std::endl;
    // std::cerr << "----- at md --------------------" << std::endl;

    // run a step
    if ((error = md.run(topo, conf, sim))){
      
      if (error == E_MINIMUM_REACHED){
	conf.old().energies.calculate_totals();
	traj.print_timestep(sim, traj.output());
	io::print_ENERGY(traj.output(), conf.old().energies, 
			 topo.energy_groups(),
			 "MINIMUM ENERGY", "EMIN_");
	  
	error = 0; // clear error condition
	break;
      }
      else { 
	// try to print energies anyway
	// if (error == E_NAN){
	io::print_ENERGY(traj.output(), conf.current().energies,
			 topo.energy_groups(),
			 "OLDERROR", "OLDERR_");
	
	io::print_ENERGY(traj.output(), conf.old().energies, 
			 topo.energy_groups(),
			 "ERROR", "ERR_");
      }

      std::cout << "\nError during MD run!\n" << std::endl;
      // try to save the final structures...
      break;
    }
    // std::cerr << "----- at md done ---------------" << std::endl;

    traj.print(topo, conf, sim);

    ++sim.steps();
    sim.time() = sim.param().step.t0 + sim.steps()*sim.time_step_size();

    ++cg_sim.steps();
    cg_sim.time() = cg_sim.param().step.t0 + cg_sim.steps()*cg_sim.time_step_size();

    if ((sim.steps() % (sim.param().step.number_of_steps * msteps / 10)) == 0){
      ++percent;
      const double spent = util::now() - start;
      const int hh = int(spent / 3600);
      const int mm = int((spent - hh * 3600) / 60);
      const int ss = int(spent - hh * 3600 - mm * 60);

      std::cerr << "MD:       " << std::setw(3) << percent * 10 << "% done..." << std::endl;
      std::cout << "MD:       " << std::setw(3) << percent * 10 << "% done..." << std::endl;
      std::cerr << "MD: spent " << std::setw(3) << hh << ":" 
		<< std::setw(2) << mm << ":" 
		<< std::setw(2) << ss << std::endl;
      std::cout << "MD: spent " << std::setw(3) << hh << ":" 
		<< std::setw(2) << mm << ":" 
		<< std::setw(2) << ss << std::endl;

      const double eta_spent = spent / sim.steps() * (msteps * sim.param().step.number_of_steps) - spent;
      const int eta_hh = int(eta_spent / 3600);
      const int eta_mm = int((eta_spent - eta_hh * 3600) / 60);
      const int eta_ss = int(eta_spent - eta_hh * 3600 - eta_mm * 60);
      
      std::cerr << "MD: ETA   " << std::setw(3) << eta_hh << ":" 
		<< std::setw(2) << eta_mm << ":"
		<< std::setw(2) << eta_ss << std::endl;
      std::cout << "MD: ETA   " << std::setw(3) << eta_hh << ":"
		<< std::setw(2) << eta_mm << ":" 
		<< std::setw(2) << eta_ss << std::endl;
    }
  }
    
  std::cout << "writing final configuration" << std::endl;
    
  traj.write(conf, topo, sim, io::final);
  traj.print_final(topo, conf, sim);

  util::update_virtual_pos(cg_topo, cg_conf, topo, conf);
  cg_traj.write(cg_conf, cg_topo, cg_sim, io::final);
    
  std::cout << "\nMESSAGES FROM SIMULATION\n";
  io::message::severity_enum err_msg = io::messages.display(std::cout);

  std::cout << "\n\n";
    
  md.print_timing(std::cout);

  std::cout << "Overall time used:\t" << util::now() - start << "\n"
	    << "(initialisation took " << init_time << ")\n\n";

  const time_t time_now = time_t(util::now());
  std::cout << ctime(&time_now) << "\n\n";
    
  if (error)
    std::cout << "\nErrors encountered during run - check above!\n" << std::endl;
  else if(err_msg > io::message::notice){
    std::cout << "\n" GROMOSXX " finished. "
	      << "Check the messages for possible problems during the run."
	      << std::endl;
  } else {
    std::cout << "\n" GROMOSXX "finished successfully\n" << std::endl;
  }
  
  return 0;
}