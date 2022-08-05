/**
 * @file aladip.t.cc
 * tests using aladip
 */


#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../interaction/interaction.h"
#include "../interaction/forcefield/forcefield.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../util/error.h"

#include "../interaction/interaction_types.h"
#include "../io/instream.h"
#include "../util/parse_tcouple.h"
#include "../io/blockinput.h"
#include "../io/topology/in_topology.h"
#include "../io/message.h"

#include "../algorithm/integration/leap_frog.h"
#include "../algorithm/temperature/temperature_calculation.h"
#include "../algorithm/temperature/berendsen_thermostat.h"
#include "../algorithm/pressure/pressure_calculation.h"
#include "../algorithm/pressure/berendsen_barostat.h"

#include "../interaction/forcefield/create_forcefield.h"

#include "../util/create_simulation.h"
#include "../algorithm/create_md_sequence.h"

#include <time.h>

#include "check.h"

#include "check_forcefield.h"
#include "check_state.h"

void hard_coded_values(std::map<std::string, double> & m){
  m["QuarticBond"] = 18.055276;
  m["PerturbedQuarticBond"] = 1.271962;
  m["Angle"] = 12.170290;
  m["PerturbedAngle"] =  0.500162;
  m["ImproperDihedral"] = 0.965060;
  m["PerturbedImproperDihedral"] =  0.560988;
  m["Dihedral"] = 2.255206;
  m["PerturbedDihedral"] = 4.626301;
  m["NonBonded_cg"] = -7.352312;
  m["NonBonded"] =  -52.832175;
  m["NonBonded_newRF"] = -70.855468;
  m["NonBonded_atomic"] =  -49.912;
  m["DistanceRestraint"] = 257.189539;
  m["PerturbedDistanceRestraint"] = 195.899012;
  m["DihedralRestraint"] = 2127.910749;
  m["PerturbedDihedralRestraint"] = 279.207857;
}


#ifdef OMP
  #include <omp.h>
#endif

int main(int argc, char* argv[]) {

#ifdef OMP
  omp_set_num_threads(1);
#endif

  int total = 0;
  
  util::Known knowns;
  knowns << "topo" << "pttopo" << "conf" << "input" << "verb";

  std::string usage = argv[0];
  usage += "\n\t[@topo       <topology>]\n";
  usage += "\t[@pttopo     <perturbation topology>]\n";
  usage += "\t[@conf       <starting configuration>]\n";
  usage += "\t[@input_off  <input without individual lambdas>]\n";
  usage += "\t[@input_on   <input with individual lambdas>]\n";
  usage += "\t[@verb       <[module:][submodule:]level>]\n";

  io::Argument args;
  if (args.parse(argc, argv, knowns, true)){
    std::cerr << usage << std::endl;
    return 1;
  }

  // parse the verbosity flag and set debug levels
  util::parse_verbosity(args);
      
  std::string stopo, spttopo, sconf, sinputon, sinputoff;
  bool quiet = true;

  if (args.count("verb") != -1) quiet = false;
      
  if(args.count("topo") == 1)
    stopo = args["topo"];
  else
    GETFILEPATH(stopo, "aladip.topo", "src/check/data/");

  if(args.count("pttopo") == 1)
    spttopo = args["pttopo"];
  else
    GETFILEPATH(spttopo, "lambdas.pttopo", "src/check/data/");
      
  if(args.count("conf") == 1)
    sconf = args["conf"];
  else
    GETFILEPATH(sconf, "aladip.conf", "src/check/data/");

  if(args.count("input_off") == 1)
    sinputoff = args["input_off"];
  else
    GETFILEPATH(sinputoff, "lambdas.off.in", "src/check/data/");
  
  if(args.count("input_on") == 1)
    sinputon = args["input_on"];
  else
    GETFILEPATH(sinputon, "lambdas.on.in", "src/check/data/");
  
  if (!quiet)
    std::cout << "\n\n"
	      << "topology :      " << stopo << "\n"
	      << "perturbation :  " << spttopo << "\n"
	      << "input (off):    " << sinputoff << "\n"
	      << "input (on):     " << sinputon << "\n"
	      << "configuration : " << sconf << "\n"
	      << std::endl;

  // set hard coded values to compare to
  std::map<std::string, double> ref_values;
  hard_coded_values(ref_values);
  
  util::simulation_struct aladip_sim_on;
  util::simulation_struct aladip_sim_off;
  
  io::In_Topology in_topo_on, in_topo_off;

  in_topo_on.quiet = quiet;
  in_topo_off.quiet = quiet;
      
  if (util::create_simulation(stopo,
			      spttopo,
			      sconf,
			      sinputon,
			      aladip_sim_on,
			      in_topo_on,
			      "", "", "", "", "", "", "", "",
			      quiet
			      )
      != 0){
    std::cerr << "creating simulation (on) failed!" << std::endl;
    return 1;
  }
  io::messages.display(std::cout);
  io::messages.clear();

  if (util::create_simulation(stopo,
			      spttopo,
			      sconf,
			      sinputoff,
			      aladip_sim_off,
			      in_topo_off,
			      "", "", "", "", "", "", "", "",
			      quiet
			      )
      != 0){
    std::cerr << "creating simulation (off) failed!" << std::endl;
    return 1;
  }     
  io::messages.display(std::cout);
  io::messages.clear();
  
  // create a forcefield
  interaction::Forcefield *ff_on = new interaction::Forcefield;
  
  if (interaction::create_g96_forcefield(*ff_on, 
					 aladip_sim_on.topo,
					 aladip_sim_on.sim,
					 in_topo_on,
					 std::cout,
					 quiet)
      != 0){
    std::cerr << "creating forcefield failed!" << std::endl;
    return 1;
  }
  
  ff_on->init(aladip_sim_on.topo, aladip_sim_on.conf, aladip_sim_on.sim, 
	   std::cout,  quiet);
  
  // first check the forcefield for the on simulation
  total += check::check_forcefield(aladip_sim_on.topo, aladip_sim_on.conf, 
				   aladip_sim_on.sim, *ff_on, ref_values);
  
  // check virial, ...
  total += check::check_state(aladip_sim_on.topo, aladip_sim_on.conf, 
			      aladip_sim_on.sim, *ff_on);

  // create and apply the force field for the off-simulation
  interaction::Forcefield *ff_off = new interaction::Forcefield;
  
  if (interaction::create_g96_forcefield(*ff_off, 
					 aladip_sim_off.topo,
					 aladip_sim_off.sim,
					 in_topo_off,
					 std::cout,
					 quiet)
      != 0){
    std::cerr << "creating forcefield failed!" << std::endl;
    return 1;
  }
  ff_off->init(aladip_sim_off.topo, aladip_sim_off.conf, aladip_sim_off.sim, 
	      std::cout,  quiet);
      
  // calculate all dl_int / dl derivatives directly from the input parameters
  // we can only check the bonded and special interactions this way
  int n_groups = aladip_sim_off.topo.energy_groups().size();
  double lam = aladip_sim_on.topo.lambda();
  if(lam != aladip_sim_off.topo.lambda()){
    std::cerr << "on and of simulations are at different lambdas" << std::endl;
    return 1;
  }

  double a=aladip_sim_on.sim.param().lambdas.a[simulation::mass_lambda][0][0];
  double b=aladip_sim_on.sim.param().lambdas.b[simulation::mass_lambda][0][0];
  double c=aladip_sim_on.sim.param().lambdas.c[simulation::mass_lambda][0][0];
  double d=aladip_sim_on.sim.param().lambdas.d[simulation::mass_lambda][0][0];
  double e=aladip_sim_on.sim.param().lambdas.e[simulation::mass_lambda][0][0];
  bool do_kin=true;

  for(int i=0; i< simulation::last_interaction_lambda; i++){
    std::string nm;
    double E_on = 0.0, E_off = 0.0, dE_on = 0.0, dE_off = 0.0;

    for(int n1=0; n1 < n_groups; n1++){
      int res=0;
      
      if(aladip_sim_on.sim.param().lambdas.a[i][n1][n1] != 0 ||
	 aladip_sim_on.sim.param().lambdas.b[i][n1][n1] != 0 ||
	 aladip_sim_on.sim.param().lambdas.c[i][n1][n1] != 0 ||
	 aladip_sim_on.sim.param().lambdas.d[i][n1][n1] != 1 ||
	 aladip_sim_on.sim.param().lambdas.e[i][n1][n1] != 0){
	// we have an individual lambda interaction
	double Lint=
	  aladip_sim_on.sim.param().lambdas.a[i][n1][n1] * lam * lam * lam * lam +
	  aladip_sim_on.sim.param().lambdas.b[i][n1][n1] * lam * lam * lam +
	  aladip_sim_on.sim.param().lambdas.c[i][n1][n1] * lam * lam +
	  aladip_sim_on.sim.param().lambdas.d[i][n1][n1] * lam +
	  aladip_sim_on.sim.param().lambdas.e[i][n1][n1];
	double dLint = 
	  4.0 * aladip_sim_on.sim.param().lambdas.a[i][n1][n1] * lam * lam * lam +
	  3.0 * aladip_sim_on.sim.param().lambdas.b[i][n1][n1] * lam * lam +
	  2.0 * aladip_sim_on.sim.param().lambdas.c[i][n1][n1] * lam +
	  aladip_sim_on.sim.param().lambdas.d[i][n1][n1];
	
	aladip_sim_off.topo.lambda(Lint);
	aladip_sim_off.topo.update_for_lambda();
	ff_off->apply(aladip_sim_off.topo, 
		      aladip_sim_off.conf,
		      aladip_sim_off.sim);
	aladip_sim_off.conf.current().perturbed_energy_derivatives.
	  calculate_totals();
	if(i==simulation::bond_lambda){
	  E_on=aladip_sim_on.conf.current().energies.
	    bond_energy[n1];
	  E_off=aladip_sim_off.conf.current().energies.
	    bond_energy[n1];
	  dE_on=aladip_sim_on.conf.current().perturbed_energy_derivatives.
	    bond_energy[n1];
	  dE_off=aladip_sim_off.conf.current().perturbed_energy_derivatives.
	    bond_energy[n1];
	  nm="bond";
	}
	else if(i==simulation::angle_lambda){
	  E_on=aladip_sim_on.conf.current().energies.
	    angle_energy[n1];
	  E_off=aladip_sim_off.conf.current().energies.
	    angle_energy[n1];
	  dE_on=aladip_sim_on.conf.current().perturbed_energy_derivatives.
	    angle_energy[n1];
	  dE_off=aladip_sim_off.conf.current().perturbed_energy_derivatives.
	    angle_energy[n1];
	  nm="angle";
	}
	else if(i==simulation::improper_lambda){
	  E_on=aladip_sim_on.conf.current().energies.
	    improper_energy[n1];
	  E_off=aladip_sim_off.conf.current().energies.
	    improper_energy[n1];
	  dE_on=aladip_sim_on.conf.current().perturbed_energy_derivatives.
	    improper_energy[n1];
	  dE_off=aladip_sim_off.conf.current().perturbed_energy_derivatives.
	    improper_energy[n1];
	  nm="improper";
	}
	else if(i==simulation::dihedral_lambda){
	  E_on=aladip_sim_on.conf.current().energies.
	    dihedral_energy[n1];
	  E_off=aladip_sim_off.conf.current().energies.
	    dihedral_energy[n1];
	  dE_on=aladip_sim_on.conf.current().perturbed_energy_derivatives.
	    dihedral_energy[n1];
	  dE_off=aladip_sim_off.conf.current().perturbed_energy_derivatives.
	    dihedral_energy[n1];
	  nm="dihedral";
	}
	else if(i==simulation::angres_lambda){
	  E_on=aladip_sim_on.conf.current().energies.
	    angrest_energy[n1];
	  E_off=aladip_sim_off.conf.current().energies.
	    angrest_energy[n1];
	  dE_on=aladip_sim_on.conf.current().perturbed_energy_derivatives.
	    angrest_energy[n1];
	  dE_off=aladip_sim_off.conf.current().perturbed_energy_derivatives.
	    angrest_energy[n1];
	  nm="angle restraint";
	}
	else if(i==simulation::dihres_lambda){
	  E_on=aladip_sim_on.conf.current().energies.
	    dihrest_energy[n1];
	  E_off=aladip_sim_off.conf.current().energies.
	    dihrest_energy[n1];
	  dE_on=aladip_sim_on.conf.current().perturbed_energy_derivatives.
	    dihrest_energy[n1];
	  dE_off=aladip_sim_off.conf.current().perturbed_energy_derivatives.
	    dihrest_energy[n1];
	  nm="dihedral restraint";
	}
	else if(i==simulation::disres_lambda){
	  E_on=aladip_sim_on.conf.current().energies.
	    distanceres_energy[n1];
	  E_off=aladip_sim_off.conf.current().energies.
	    distanceres_energy[n1];
	  dE_on=aladip_sim_on.conf.current().perturbed_energy_derivatives.
	    distanceres_energy[n1];
	  dE_off=aladip_sim_off.conf.current().perturbed_energy_derivatives.
	    distanceres_energy[n1];
	  nm="distance restraint";
	}
	// let's not claim we checked something that was zero to begin with
	if(dE_on !=0){
	    
	  CHECKING("individual lambdas ("+nm+")", res);
	  CHECK_APPROX_EQUAL(E_on, E_off, 0.0000001, res);
	  CHECK_APPROX_EQUAL(dE_on, dE_off*dLint, 0.0000001, res);
	  RESULT(res, total);
	}
      }
      
      // check if we can do the kinetic energy
      if(i==simulation::mass_lambda &&
	 (aladip_sim_on.sim.param().lambdas.a[i][n1][n1] != a ||
	  aladip_sim_on.sim.param().lambdas.b[i][n1][n1] != b ||
	  aladip_sim_on.sim.param().lambdas.c[i][n1][n1] != c ||
	  aladip_sim_on.sim.param().lambdas.d[i][n1][n1] != d ||
	  aladip_sim_on.sim.param().lambdas.e[i][n1][n1] != e)){
	do_kin=false;
      }
    } // loop over energy groups
  } // loop over interactions
  
  // kinetic (we can't really do this energy group dependent)
  if(do_kin){
    int res = 0;
    
    if(a != 0 || b != 0 || c != 0 || d != 1 || e != 0){

      // we have an individual lambda interaction
      double Lint = a * lam * lam * lam * lam +
	b * lam * lam * lam +
	c * lam * lam +
	d * lam +
	e;
      double dLint = 
	4.0 * a * lam * lam * lam +
	3.0 * b * lam * lam +
	2.0 * c * lam +
	d;
      aladip_sim_off.topo.lambda(Lint);
      aladip_sim_off.topo.update_for_lambda();
      ff_off->apply(aladip_sim_off.topo, 
		    aladip_sim_off.conf,
		    aladip_sim_off.sim);
      aladip_sim_off.conf.current().perturbed_energy_derivatives.
	calculate_totals();
      // we probably need to do some temperature calculation as well?
      CHECKING("indivual lambdas (kinetic)", res);
      double dE_on=aladip_sim_on.conf.current().perturbed_energy_derivatives.
	kinetic_total;
      double dE_off=aladip_sim_off.conf.current().perturbed_energy_derivatives.
	kinetic_total;
      CHECK_APPROX_EQUAL(dE_on, dE_off*dLint, 0.0000001, res);
      RESULT(res, total);
    }
  }
  io::messages.display(std::cout);
  io::messages.clear();

  return total;
}
