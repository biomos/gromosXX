/**
 * @file script.cc
 * example for a simple md script
 */

#include <util/stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction.h>
#include <interaction/interaction_types.h>
#include <interaction/forcefield/forcefield.h>
#include <interaction/forcefield/create_forcefield.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/error.h>

#include <time.h>

#include <io/forcefield/in_ifp.h>
#include <io/configuration/out_configuration.h>

int main(int argc, char *argv[]){

  try{
    
    char *knowns[] = 
      {
        "ifp", "verb"
      };
    
    int nknowns = 2;
    
    std::string usage = argv[0];
    usage += "\n\t@ifp    <interaction function parameter file>\n";
    usage += "\t[@verb   <[module:][submodule:]level>]\n";

    io::Argument args(argc, argv, nknowns, knowns, usage);
    
    // parse the verbosity flag and set debug levels
    util::parse_verbosity(args);

    // create the simulation classes
    topology::Topology topo;
    configuration::Configuration conf;
    simulation::Simulation sim;

    // add two atoms
    topo.resize(2);
    topo.residue_names().push_back("AR");
    topo.residue_names().push_back("AR");
    
    topo.add_solute_atom("AR", 0, 27, 39, 0, 1, std::set<int>(), std::set<int>());
    topo.add_solute_atom("AR", 1, 27, 39, 0, 1, std::set<int>(), std::set<int>());
    
    topo.initialise();

    // perturb one of the atoms
    std::cout <<"perturb an atom" << std::endl;

    topology::Perturbed_Atom atom(0, 
				  27, 39, 0,
				  27, 39, 0,
				  0.1, 0.005);
    
    topo.perturbed_solute().atoms()[0] = atom;
    topo.is_perturbed()[0] = true;

    topo.lambda(0);
    topo.lambda_exp(1);
    
    conf.current().box(0) = math::Vec(5, 0, 0);
    conf.current().box(1) = math::Vec(0, 5, 0);
    conf.current().box(2) = math::Vec(0, 0, 5);
    conf.boundary_type = math::rectangular;
    sim.param().boundary.boundary = math::rectangular;

    // no gather
    conf.initialise(topo, sim.param(), false);
    
    // create a forcefield
    std::ifstream fifp(args["ifp"].c_str());
    if (!fifp.is_open())
      throw std::string("could not open interaction parameter file " + args["ifp"]);
    
    io::In_IFP ifp(fifp);

    sim.param().perturbation.perturbation = true;

    interaction::Forcefield ff;
    interaction::create_g96_forcefield(ff,
				       topo,
				       sim,
				       conf,
				       ifp);

    // start the thing...
    conf.current().pos(0) = math::Vec(0, 0, 0);
    conf.current().pos(1) = math::Vec(1, 0, 0);
    
    io::messages.display();
    
    for(double i=0; i<1.01; i+=0.01){

      conf.current().energies.zero();
      conf.current().force = 0.0;
      
      topo.lambda(i);

      ff.apply(topo, conf, sim);
      
      conf.current().energies.calculate_totals();
      std::cout << "lambda: " 
		<< i
		<< "\tetot: "
		<< conf.current().energies.total << "\n";
      
    }
    
  }
  catch (std::string s){
    io::messages.display();
    std::cerr << "error encountered:\n" << s << std::endl;
    return 1;
  }
  
  return 0;
}
