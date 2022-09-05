/**
 * @file md.cc
 * the main md program
 */

/**
 * @page programs Program Documentation
 *
 * @anchor md
 * @section md molecular dynamics
 * @date 28.10.2008
 *
 * Program md is used to run molecular dynamics simulations. The command line
 * arguments are summarized in the following table:
 *
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@pttopo</td><td>&lt;@ref pttopo "molecular perturbation topology file"&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@conf</td><td>&lt;coordinates and restart data&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@input</td><td>&lt;@ref input "input parameters"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@fin</td><td>&lt;final configuration&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trc</td><td>&lt;coordinate trajectory&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trv</td><td>&lt;velocity trajectory&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trf</td><td>&lt;force trajectory&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trs</td><td>&lt;special trajectory&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@tre</td><td>&lt;@ref energy_trajectory "energy trajectory"&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@bae</td><td>&lt;@ref block_averaged_energy_trajectory "block averaged energy trajectory"&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trg</td><td>&lt;@ref free_energy_trajectory "free energy trajectory"&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@bag</td><td>&lt;@ref block_averaged_free_energy_trajectory "block averaged free energy trajectory"&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@posresspec</td><td>&lt;@ref posres "position restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@refpos</td><td>&lt;@ref posres "position restraints"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@distrest</td><td>&lt;@ref disres "distance restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@angrest</td><td>&lt;@ref angrest "angle restraints specification"&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@dihrest</td><td>&lt;@ref dihrest "dihedral restraints specification"&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@jval</td><td>&lt;@ref jvalue "J-value restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@xray</td><td>&lt;@ref xrayresfile "X-ray restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@sym</td><td>&lt;@ref symrest "Symmetry restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@order</td><td>&lt;@ref orderparamresspec "Order-parameter restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@lud</td><td>&lt;@ref leusdb "local elevation umbrella database"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@led</td><td>&lt;@ref leus "local elevation coordinate specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@bsleus</td><td>&lt;@ref bsleus "Ball & Stick Local Elevation topology"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@friction</td><td>&lt;@ref friction "atomic friction coefficients"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@qmmm</td><td>&lt;@ref qmmm "QM/MM specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@print</td><td>&lt;print additional information&gt; </td><td></td></tr>
 * <tr><td> \@anatrj</td><td>&lt;re-analyze trajectory&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@verb</td><td>&lt;@ref debug "control verbosity"&gt;</td><td></td></tr>
 * <tr><td> \@version</td><td>&lt;print version information&gt; </td><td></td></tr>
 * </table>
 * 
 * @sa md_mpi
*/

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>

#include <io/configuration/out_configuration.h>
#include <signal.h>

bool exit_md;

void signal_handler(int sig) {
  std::cerr << "\nThe program will exit as soon as the MD step has finished. Press CTRL-C to force quit." << std::endl;
  exit_md = true;
  signal(SIGINT, SIG_DFL);
}

int main(int argc, char *argv[]){

  const double start = util::now();
  exit_md = false;
  signal(SIGINT, signal_handler);

  util::Known knowns;
  knowns << "topo" << "conf" << "input" << "verb" << "pttopo"
	 << "trc" << "fin" << "trv" << "trf" << "trs" << "tre" << "trg"
	 << "bae" << "bag" << "posresspec" << "refpos" <<"distrest" 
	 << "angrest" << "dihrest"
         << "jval" << "xray" << "sym" << "order" << "rdc" << "lud" << "led" << "bsleus" 
         << "anatrj" << "print" << "friction" << "qmmm" << "version" << "develop";
  
  
  std::string usage;
  util::get_usage(knowns, usage, argv[0]);
  usage += "#\n\n";

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
  simulation::Simulation sim;
  algorithm::Algorithm_Sequence md;

  io::Out_Configuration traj(GROMOSXX "\n");

  if (io::read_input(args, topo, conf, sim,  md)){
    io::messages.display(std::cout);
    std::cout << "\nErrors during initialization!\n" << std::endl;
    return 1;
  }

  // check for development 
  if (sim.param().develop.develop==true && args.count("develop") < 0) { 
    io::messages.add(sim.param().develop.msg, io::message::develop); 
  } 
  
  traj.title(GROMOSXX "\n" + sim.param().title);

  // create output files...
  traj.init(args, sim.param());

  // initialises all algorithms (and therefore also the forcefield)
  md.init(topo, conf, sim);

  std::cout << "\nMESSAGES FROM INITIALISATION\n";
  {
    int iom = io::messages.display(std::cout);
    if (iom >= io::message::error) {
      std::cout << "\nErrors during initialisation!\n" << std::endl;
      return 1;
    } else if (iom == io::message::develop) {
      std::cout << "\nUse @develop to run untested code.\n" << std::endl;
      return 1;
    }
  }
      
  io::messages.clear();

  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
  std::cout << "\nenter the next level of molecular "
	    << "dynamics simulations\n" << std::endl;

    
  int percent = 0;
    
  std::cout << "==================================================\n"
	    << " MAIN MD LOOP\n"
	    << "==================================================\n"
	    << std::endl;

  int error = 0;

  const double init_time = util::now() - start;
  while(int(sim.steps()) < sim.param().step.number_of_steps && !exit_md){
      
    traj.write(conf, topo, sim, io::reduced);

    // run a step
    if ((error = md.run(topo, conf, sim))){

      if (error == E_MINIMUM_REACHED) {
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
      std::cout << "\tat step " << sim.steps() << " (time " << sim.time() << ")\n" << std::endl;
      // try to save the final structures...
      break;
    }

    if (exit_md)
      std::cout << "\nMD run terminated by SIGINT (CTRL-C)." << std::endl;

    traj.print(topo, conf, sim);

    sim.steps()=sim.steps()+sim.param().analyze.stride;
    sim.time() = sim.param().step.t0 + sim.steps()*sim.time_step_size();
    
    if ((sim.param().step.number_of_steps / 10 > 0) &&
	(sim.steps() % (sim.param().step.number_of_steps / 10) == 0)){
      percent=int(sim.steps())*10/sim.param().step.number_of_steps;
      const double spent = util::now() - start;
      const int hh = int(spent / 3600);
      const int mm = int((spent - hh * 3600) / 60);
      const int ss = int(spent - hh * 3600 - mm * 60);

      std::cerr << "MD:       " << std::setw(4) << percent * 10 << "% done..." << std::endl;
      std::cout << "MD:       " << std::setw(4) << percent * 10 << "% done..." << std::endl;
      std::cerr << "MD: spent " << hh << ":" << mm << ":" << ss << std::endl;
      std::cout << "MD: spent " << hh << ":" << mm << ":" << ss << std::endl;

      const double eta_spent = spent / sim.steps() * sim.param().step.number_of_steps - spent;
      const int eta_hh = int(eta_spent / 3600);
      const int eta_mm = int((eta_spent - eta_hh * 3600) / 60);
      const int eta_ss = int(eta_spent - eta_hh * 3600 - eta_mm * 60);
      
      std::cerr << "MD: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
      std::cout << "MD: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
    }
  } // main md loop
    
  std::cout << "writing final configuration" << std::endl;
    
  traj.write(conf, topo, sim, io::final);
  traj.print_final(topo, conf, sim);
    
  std::cout << "\nMESSAGES FROM SIMULATION\n";
  io::message::severity_enum err_msg = io::messages.display(std::cout);

  std::cout << "\n\n";
    
  md.print_timing(std::cout);

  std::cout << "Overall time used:\t" << util::now() - start << "\n"
	    << "(initialisation took " << init_time << ")\n\n";

  const time_t time_now = time_t(util::now());
  std::cout << ctime(&time_now) << "\n\n";
    
  if (error){
    std::cout << "\nErrors encountered during run - check above!\n" << std::endl;
    return 1;
  } else if (exit_md) {
     std::cout <<"\n" GROMOSXX " finished due to SIGINT. Returning non-zero value." <<std::endl;
     return 2;
  } else if(err_msg > io::message::notice){
    std::cout << "\n" GROMOSXX " finished. "
	      << "Check the messages for possible problems during the run."
	      << std::endl;
    return 0;
  } else{
    std::cout << "\n" GROMOSXX " finished successfully\n" << std::endl;
  }
  
  return 0;
}


