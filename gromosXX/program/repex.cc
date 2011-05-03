/**
 * @file repex.cc
 * the main md program for replica exchange simulations
 */
/**
 * @page programs Program Documentation
 *
 * @anchor repex
 * @section repex replica exchange
 * @date 28.10.2008
 *
 * Program repex is used to run replica exchange simulations.
 *
 * See @ref md for the documentation of the command line arguments.
 * Addition command line arguments are:
 * <table border=0 cellpadding=0>
 * <tr><td> \@cg_topo</td><td>&lt;coarse grain topology data&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@cg_pttopo</td><td>&lt;coarse grain perturbation topology file&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@cg_conf</td><td>&lt;coarse grain coordinates and restart data&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@cg_input</td><td>&lt;coarse grain input parameters&gt;</td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@cg_fin</td><td>&lt;coarse grain final configuration&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@cg_trc</td><td>&lt;coordinate trajectory&gt; </td><td style="color:#FF0000">out</td></tr>
 * <tr><td> \@master</td><td>&lt;master process&gt; </td><td></td></tr>
 * <tr><td> \@slave</td><td>&lt;slave process&gt; </td><td></td></tr>
 * </table>
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
#include <unistd.h>

#include <io/configuration/out_configuration.h>

#include <gsl/gsl_rng.h>
#include <util/replica_exchange.h>

int main(int argc, char *argv[]){
#ifdef REPEX

  util::Known knowns;
  knowns << "topo" << "cg_topo" << "conf" << "cg_conf"
     << "input" << "cg_input" << "verb" << "pttopo" << "cg_pttopo"
	 << "trc" << "cg_trc" << "fin" << "cg_fin" << "trv" << "trf" << "trs" << "tre" << "re" << "trg"
	 << "bae" << "bag" << "posresspec" << "refpos" << "distrest" << "jval" << "xray" << "lud" << "led"
         << "friction" << "rep" << "master" << "slave" << "control" << "qmmm" << "version";
  
  
  std::string usage;
  util::get_usage(knowns, usage, argv[0]);
  usage += "#\n\n";

  io::Argument args;

  if (args.parse(argc, argv, knowns)){
    std::cerr << usage << std::endl;
    return 1;
  }

  util::print_title(false, std::cout, true, args.count("master") >= 0);
  if (args.count("version") >= 0){
    return 0;
  }
    
  // parse the verbosity flag and set debug levels
  if (util::parse_verbosity(args)){
    std::cerr << "could not parse verbosity argument" << std::endl;
    return 1;
  }
  
  int error = 0;

  if (args.count("master") >= 0){

    // the master
    util::Replica_Exchange_Master rep_master;
    std::cout << "starting master thread" << std::endl;
    error = rep_master.run(args);
  }
  else if (args.count("slave") >= 0){
    
    // and the slave
    std::cout << "repex: starting slave" << std::endl;
    util::Replica_Exchange_Slave rep_slave;

    error = rep_slave.run(args);

  }
  else if (args.count("control") >= 0){
    std::cout << "repex: starting control!" << std::endl;
    util::Replica_Exchange_Control rep_control;
    error = rep_control.run(args);
  }
  else{
    std::cout << "repex: either @master @slave or @interactive required" << std::endl;
  }

  if (error){
    std::cerr << "errors during repex run! (error code " << error << ")\n"
	      << std::endl;

    io::messages.display();
    
    return 1;
  }
  
  return 0;

#else
  std::cout << argv[0] << " needs REPEX to run\n\tuse --enable-repex at configure\n" << std::endl;
  return 1;
#endif
  
}
