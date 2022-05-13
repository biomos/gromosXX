/**
 * @file eds_2box.cc
 * the main md program for 2-box eds simulations
 */

/**
 * @page programs Program Documentation
 *
 * @anchor eds_2box
 * @section eds_2box molecular dynamics
 * @date 17.08.2011
 *
 * Program eds_2box is used to run eds simulations with 2 boxes. The command line
 * arguments are summarized in the following table:
 *
 * <table border=0 cellpadding=0>
 * <tr><td> \@topoX</td><td>&lt;molecular topology file for box X&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@pttopoX</td><td>&lt;@ref pttopo "molecular perturbation topology file for box X"&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@confX</td><td>&lt;coordinates and restart data for box X&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@inputX</td><td>&lt;@ref input "input parameters for box X"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@finX</td><td>&lt;final configuration for box X&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trcX</td><td>&lt;coordinate trajectory for box X&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trvX</td><td>&lt;velocity trajectory for box X&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trfX</td><td>&lt;force trajectory for box X&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trsX</td><td>&lt;special trajectory for box X&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@treX</td><td>&lt;@ref energy_trajectory "energy trajectory"&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@bae</td><td>&lt;@ref block_averaged_energy_trajectory "block averaged energy trajectory"&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@trg</td><td>&lt;@ref free_energy_trajectory "free energy trajectory"&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@bag</td><td>&lt;@ref block_averaged_free_energy_trajectory "block averaged free energy trajectory"&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@posresspec</td><td>&lt;@ref posres "position restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@refpos</td><td>&lt;@ref posres "position restraints"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@distrest</td><td>&lt;@ref disres "distance restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@dihtrest</td><td>&lt;@ref dihrest "dihedral restraints specification"&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@jval</td><td>&lt;@ref jvalue "J-value restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@xray</td><td>&lt;@ref xrayresfile "X-ray restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@order</td><td>&lt;@ref orderparamresspec "Order-parameter restraints specification"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@lud</td><td>&lt;@ref leusdb "local elevation umbrella database"&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@led</td><td>&lt;@ref leus "local elevation coordinate specification"&gt;</td><td style="color:#088A08">in</td></tr>
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
#include <algorithm/integration/eds.h>

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
  knowns << "topo1" << "topo2" << "conf1" << "conf2" << "input1" << "input2" << "verb" 
         << "pttopo1" << "pttopo2" << "trc1" << "trc2" << "fin1" << "fin2"
	 << "trv1" << "trv2" << "trf1" << "trf2" << "trs1" << "trs2"  
         << "tre1" << "tre2" << "trg1" << "trg2" 
	 << "bae" << "bag" << "posresspec" << "refpos" << "distrest1" << "distrest2" << "dihrest"
         << "jval" << "xray" << "order" << "rdc" << "lud" << "led" << "anatrj" << "print" << "friction"
         << "qmmm" << "version";
  
  
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
  topology::Topology topo1, topo2;
  configuration::Configuration conf1, conf2;
  simulation::Simulation sim1, sim2;
  algorithm::Algorithm_Sequence md1, md2;
  algorithm::Algorithm_Sequence md1_first(false), md2_first(false),
          md1_second(false), md2_second(false);

  io::Out_Configuration traj1(GROMOSXX "\n"), traj2(GROMOSXX "\n");
  io::Argument args1, args2;
  for(io::Argument::const_iterator it = args.begin(), to = args.end();
          it != to; ++it) {
    const std::string & argname = it->first;
    if (argname[argname.size() - 1] == '1') {
      std::string newargname(argname.substr(0, argname.size() - 1));
      args1.insert(std::pair<std::string, std::string>(newargname, it->second));
    }
  }
  
  for(io::Argument::const_iterator it = args.begin(), to = args.end();
          it != to; ++it) {
    const std::string & argname = it->first;
    if (argname[argname.size() - 1] == '2') {
      std::string newargname(argname.substr(0, argname.size() - 1));
      args2.insert(std::pair<std::string, std::string>(newargname, it->second));
    }
  }
  
  
  if (io::read_input(args1, topo1, conf1, sim1,  md1)){
    io::messages.display(std::cout);
    std::cout << "\nErrors during initialization of box 1!\n" << std::endl;
    return 1;
  }

  if (io::read_input(args2, topo2, conf2, sim2,  md2)){
   io::messages.display(std::cout);
   std::cout << "\nErrors during initialization of box 2!\n" << std::endl;
   return 1;
  }
  
  traj1.title(GROMOSXX "\n" + sim1.param().title);
  traj2.title(GROMOSXX "\n" + sim2.param().title);

  // create output files...
  traj1.init(args1, sim1.param());
  traj2.init(args2, sim2.param());
  
  algorithm::EDS * eds1 = (algorithm::EDS *) md1.algorithm("EDS");
  if (eds1 != NULL) {
      eds1->set_conf2(conf2);
  }
  algorithm::EDS * eds2 = (algorithm::EDS *) md2.algorithm("EDS");
  if (eds2 != NULL) {
      eds2->set_conf2(conf1);
  }
  
  // initialises all algorithms (and therefore also the forcefield)
  md1.init(topo1, conf1, sim1);
  md2.init(topo2, conf2, sim2);
  
  bool second = false;
  for(algorithm::Algorithm_Sequence::const_iterator it = md1.begin(),
          to = md1.end(); it != to; ++it) {
    if (*it == eds1) {
      second = true;
      continue;
    }
    
    if (second)
      md1_second.push_back(*it);
    else
      md1_first.push_back(*it);
  }
  second = false;
  for(algorithm::Algorithm_Sequence::const_iterator it = md2.begin(),
          to = md2.end(); it != to; ++it) {
    if (*it == eds2) {
      second = true;
      continue;
    }
    
    if (second)
      md2_second.push_back(*it);
    else
      md2_first.push_back(*it);
  }

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

    
  int percent = 0;
    
  std::cout << "==================================================\n"
	    << " MAIN MD LOOP\n"
	    << "==================================================\n"
	    << std::endl;

  int error1 = 0, error2 = 0;

  const double init_time = util::now() - start;
  while(int(sim1.steps()) < sim1.param().step.number_of_steps &&
          int(sim2.steps()) < sim2.param().step.number_of_steps  && !exit_md){
      
    traj1.write(conf1, topo1, sim1, io::reduced);
    traj2.write(conf2, topo2, sim2, io::reduced);

    // run a step
    error1 = md1_first.run(topo1, conf1, sim1);
    error2 = md2_first.run(topo2, conf2, sim2);
    error1 += eds1->apply(topo1, conf1, sim1);
    error2 += eds2->apply(topo2, conf2, sim2);
    error1 += md1_second.run(topo1, conf1, sim1);
    error2 += md2_second.run(topo2, conf2, sim2);
    if ((error1 || error2)){
    //if ((error = md.run(topo, conf, sim))){

      if (error1 == E_MINIMUM_REACHED){
	conf1.old().energies.calculate_totals();
	traj1.print_timestep(sim1, traj1.output());
	io::print_ENERGY(traj1.output(), conf1.old().energies, 
			 topo1.energy_groups(),
			 "MINIMUM ENERGY", "EMIN_");
	  
	error1 = 0; // clear error condition
	break;
      }
      else { 
	// try to print energies anyway
	// if (error == E_NAN){
	io::print_ENERGY(traj1.output(), conf1.current().energies,
			 topo1.energy_groups(),
			 "OLDERROR", "OLDERR_");
	
        io::print_ENERGY(traj2.output(), conf2.current().energies,
			 topo2.energy_groups(),
			 "OLDERROR", "OLDERR_");

	io::print_ENERGY(traj1.output(), conf1.old().energies, 
			 topo1.energy_groups(),
			 "ERROR", "ERR_");
        
	io::print_ENERGY(traj2.output(), conf2.old().energies, 
			 topo2.energy_groups(),
			 "ERROR", "ERR_");

      }

      std::cout << "\nError during MD run!\n" << std::endl;
      std::cout << "\tat step " << sim1.steps() << " (time " << sim1.time() << ")\n" << std::endl;
      // try to save the final structures...
      break;
    }

    if (exit_md)
      std::cout << "\nMD run terminated by SIGINT (CTRL-C)." << std::endl;

    traj1.print(topo1, conf1, sim1);
    traj2.print(topo2, conf2, sim2);

    ++sim1.steps();
    ++sim2.steps();
    sim1.time() = sim1.param().step.t0 + sim1.steps()*sim1.time_step_size();
    sim2.time() = sim2.param().step.t0 + sim2.steps()*sim2.time_step_size();
    
    if ((sim1.param().step.number_of_steps / 10 > 0) &&
	(sim1.steps() % (sim1.param().step.number_of_steps / 10) == 0)){
      ++percent;
      const double spent = util::now() - start;
      const int hh = int(spent / 3600);
      const int mm = int((spent - hh * 3600) / 60);
      const int ss = int(spent - hh * 3600 - mm * 60);

      std::cerr << "MD:       " << std::setw(4) << percent * 10 << "% done..." << std::endl;
      std::cout << "MD:       " << std::setw(4) << percent * 10 << "% done..." << std::endl;
      std::cerr << "MD: spent " << hh << ":" << mm << ":" << ss << std::endl;
      std::cout << "MD: spent " << hh << ":" << mm << ":" << ss << std::endl;

      const double eta_spent = spent / sim1.steps() * sim1.param().step.number_of_steps - spent;
      const int eta_hh = int(eta_spent / 3600);
      const int eta_mm = int((eta_spent - eta_hh * 3600) / 60);
      const int eta_ss = int(eta_spent - eta_hh * 3600 - eta_mm * 60);
      
      std::cerr << "MD: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
      std::cout << "MD: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
    }
  } // main md loop
    
  std::cout << "writing final configuration" << std::endl;
    
  traj1.write(conf1, topo1, sim1, io::final);
  traj1.print_final(topo1, conf1, sim1);
  
  traj2.write(conf2, topo2, sim2, io::final);
  traj2.print_final(topo2, conf2, sim2);

    
  std::cout << "\nMESSAGES FROM SIMULATION\n";
  io::message::severity_enum err_msg = io::messages.display(std::cout);

  std::cout << "\n\n";
    
  md1.print_timing(std::cout);

  std::cout << "Overall time used:\t" << util::now() - start << "\n"
	    << "(initialisation took " << init_time << ")\n\n";

  const time_t time_now = time_t(util::now());
  std::cout << ctime(&time_now) << "\n\n";
    
  if (error1){
    std::cout << "\nErrors encountered in box 1 during run - check above!\n" << std::endl;
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


