/**
 * @file create_simulation.cc
 * implementation of function create_simulation
 * to easily create a simulation from (non-complete) data.
 */

#include "../stdheader.h"
#include <fstream>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../interaction/interaction.h"
#include "../interaction/interaction_types.h"

#include "../io/argument.h"
#include "../io/blockinput.h"
#include "../io/instream.h"
#include "../io/configuration/inframe.h"
#include "../io/configuration/in_configuration.h"
#include "../io/topology/in_topology.h"
#include "../io/topology/in_perturbation.h"
#include "../io/topology/in_distanceres.h"
#include "../io/topology/in_dihrest.h"
#include "../io/topology/in_xray.h"
#include "../io/topology/in_leus.h"
#include "../io/topology/in_order.h"
#include "../io/parameter/in_parameter.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../algorithm/create_md_sequence.h"

#include "../interaction/forcefield/forcefield.h"

#include "create_simulation.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util

int util::create_simulation(std::string topo,
			    std::string pttopo,
			    std::string conf,
			    std::string param,
			    util::simulation_struct & sim,
			    io::In_Topology & in_topo,
			    std::string distanceres,
			    std::string dihrest,
                            std::string xray,
			    std::string led,
			    std::string lud,
                            std::string order,
			    bool quiet)
{

  // a topology is needed
  if (topo == ""){
    io::messages.add("No topology specified",
		     "create_simulation",
		     io::message::error);
    return 1;
  }

  io::igzstream input_file, topo_file, pttopo_file, conf_file, 
    distanceres_file, dihrest_file, xray_file, led_file, lud_file, order_file;
  
  // if we got a parameter file, try to read it...
  if (param != ""){

    input_file.open(param.c_str());
    
    if (!input_file.is_open()){
      std::cout << "\n\ncould not open " << param << "!\n" << std::endl;
      io::messages.add("opening input failed", "read_input",
		       io::message::error);
      return -1;
    }
    
    io::In_Parameter ip(input_file);
    ip.quiet = quiet;
    ip.read(sim.sim.param());

    sim.sim.time_step_size() = sim.sim.param().step.dt;
    sim.sim.time() = sim.sim.param().step.t0;

  }
  
  
  topo_file.open(topo.c_str());

  if (!topo_file.is_open()){
    std::cout << "\n\ncould not open " << topo << "!\n" << std::endl;
    io::messages.add("opening topology failed", "read_input",
		     io::message::error);
    return -1;
  }

  in_topo.stream(topo_file);
  in_topo.read(sim.topo, sim.sim.param());

  if(pttopo != ""){
    
    pttopo_file.open(pttopo.c_str());

    if (!pttopo_file.is_open()){
      std::cout << "\n\ncould not open " << pttopo << "!\n" << std::endl;
      io::messages.add("opening perturbation topology failed", "read_input",
		       io::message::error);
      return -1;
    }
    
    io::In_Perturbation ipt(pttopo_file);
    ipt.quiet = quiet;
    ipt.read(sim.topo, sim.sim.param());

  }

  sim.topo.init(sim.sim);

  // read special
  if (distanceres != ""){
    
    distanceres_file.open(distanceres.c_str());

    if (!distanceres_file.is_open()){
      std::cout << "\n\ncould not open " << distanceres << "!\n" << std::endl;
      io::messages.add("opening distance restraints failed", "read_input",
		       io::message::error);
      return -1;
    }

    io::In_Distanceres idr(distanceres_file);
    idr.quiet = quiet;
    idr.read(sim.topo, sim.sim);
  }
  if (dihrest != ""){
    
    dihrest_file.open(dihrest.c_str());
    
    if(!dihrest_file.is_open()){
      std::cout << "\n\ncould not open " << dihrest << "!\n" << std::endl;
      io::messages.add("opening dihedral restraints failed", "read_input",
		       io::message::error);
      return -1;
    }

    io::In_Dihrest idr(dihrest_file);
    idr.quiet = quiet;
    idr.read(sim.topo, sim.sim);
  }

  
  if (xray != "") {
    xray_file.open(xray.c_str());
    if (!xray_file.is_open()) {
      std::cout << "\n\ncould not open " << xray << "!\n" << std::endl;
      io::messages.add("opening xray restraints failed", "read_input",
              io::message::error);
      return -1;
    }
    io::In_Xrayresspec ixr(xray_file);
    ixr.quiet = true;
    ixr.read(sim.topo, sim.sim);
  }

  // do this after reading in a perturbation topology
  sim.sim.multibath().calculate_degrees_of_freedom(sim.topo,
          sim.sim.param().rottrans.rottrans,
          sim.sim.param().posrest.posrest == simulation::posrest_const,
          sim.sim.param().boundary.dof_to_subtract);

  if (conf != ""){

    conf_file.open(conf.c_str());

    if (!conf_file.is_open()){
      std::cout << "\n\ncould not open " << conf << "!\n" << std::endl;
      io::messages.add("opening configuration failed", "read_input",
		       io::message::error);
      return -1;
    }

    io::In_Configuration ic(conf_file);
    ic.quiet = quiet;
    ic.read(sim.conf, sim.topo, sim.sim);
    sim.conf.init(sim.topo, sim.sim.param());
    
  }

  if (lud != "" && led != "") {
    
    lud_file.open(lud.c_str());
    
    if(!lud_file.is_open()){
      std::cout << "\n\ncould not open " << lud << "!\n" << std::endl;
      io::messages.add("opening local elevation umbrella file (lud) failed", "read_input",
		       io::message::error);
      return -1;
    }      
    
    io::In_LEUSBias ilud(lud_file);
    ilud.quiet = quiet;
    ilud.read(sim.topo, sim.conf, sim.sim);


    led_file.open(led.c_str());
    
    if(!led_file.is_open()){
      std::cout << "\n\ncould not open " << led << "!\n" << std::endl;
      io::messages.add("opening local elevation definition file (led) failed", "read_input",
		       io::message::error);
      return -1;
    }

    io::In_Localelevspec iled(led_file);
    iled.quiet = quiet;
    iled.read(sim.topo, sim.sim);
  }
  
    if (order != "") {
    order_file.open(order.c_str());
    if (!order_file.is_open()) {
      std::cout << "\n\ncould not open " << order << "!\n" << std::endl;
      io::messages.add("opening order parameter specification file failed", "read_input",
              io::message::error);
      return -1;
    }
    io::In_Orderparamresspec iorder(order_file);
    iorder.quiet = true;
    iorder.read(sim.topo, sim.sim);
    
    //io::In_Configuration::read_order_parameter_restraint_averages  order_res(sim.conf);
    //io::In_Configuration::read_order_parameter_restraint_averages(sim.topo, sim.conf, sim.sim, std::cout);
    
  }


  return 0;
}



