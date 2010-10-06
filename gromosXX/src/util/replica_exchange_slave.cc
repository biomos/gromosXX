/**
 * @file replica_exchange_slave.cc
 * replica exchange
 */

#include <stdheader.h>

#ifdef REPEX

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>
#include <interaction/special/external_interaction.h>

#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen_thermostat.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/error.h>
#include <util/virtual_grain.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>

#include <io/configuration/out_configuration.h>

#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
// linux includes?
#include <sys/types.h>
#include <sys/socket.h>
// end linux includes
#include <netdb.h>

#include "replica_exchange.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica

////////////////////////////////////////////////////////////////////////////////
// replica slave    ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

util::Replica_Exchange_Slave::Replica_Exchange_Slave() {
}

int util::Replica_Exchange_Slave::run
(
        io::Argument & args) {
  multigraining = false;
  if (args.count("cg_topo") >= 0) {
    multigraining = true;
  }

  // create the simulation classes
  topology::Topology topo;
  configuration::Configuration conf;
  simulation::Simulation sim;
  algorithm::Algorithm_Sequence md;

  io::Out_Configuration traj(GROMOSXX "\n\tfine-grained\n", std::cout);

  topology::Topology cg_topo;
  configuration::Configuration cg_conf;
  algorithm::Algorithm_Sequence cg_md;
  simulation::Simulation cg_sim;

  io::Out_Configuration cg_traj(GROMOSXX "\n\tcoarse-grained\n", std::cout);

  init(args,
          topo, conf, sim, md, traj,
          cg_topo, cg_conf, cg_sim, cg_md, cg_traj);

  // aliases
  std::vector<double> const & T = sim.param().replica.temperature;
  std::vector<double> const & l = sim.param().replica.lambda;
  std::vector<double> const & dt = sim.param().replica.dt;

  struct addrinfo *addrinfo_p;
  struct addrinfo hints;
  std::string server_name;

  try {
    addrinfo_p = get_server(args, hints, server_name);
  }  catch (std::runtime_error e) {
    std::cout << "Exception: " << e.what() << std::endl;
    std::cerr << "Exception: " << e.what() << std::endl;
    return 1;
  }

  sockaddr * s_addr_p = addrinfo_p->ai_addr;
  int len = addrinfo_p->ai_addrlen;

  ////////////////////////////////////////////////////////////////////////////////
  // MAIN SLAVE LOOP
  ////////////////////////////////////////////////////////////////////////////////


  for (int run = 0; run < sim.param().replica.slave_runs; ++run) {

    cl_socket = socket(addrinfo_p->ai_family, addrinfo_p->ai_socktype,
            addrinfo_p->ai_protocol);
    if (cl_socket < 0) {
      std::cerr << "could not create client socket" << std::endl;
      return 1;
    }

    DEBUG(8, "slave: connecting..");
    int result = connect(cl_socket, s_addr_p, len);

    if (result == -1) {
      std::cout << "could not connect to master. master finished?"
              << std::endl;
      return 1;
    }

    if (!magic_cookie(false)) {
      close(cl_socket);
      continue;
    }

    DEBUG(9, "slave: connected");

    // request a job
    DEBUG(8, "slave: requesting job");
    char ch = 1;
    if (write(cl_socket, &ch, 1) != 1) {
      std::cerr << "could not write to socket" << std::endl;
      close(cl_socket);
      return 1;
    }

    char server_response;

    if (read(cl_socket, &server_response, 1) != 1) {
      std::cerr << "could not read" << std::endl;
      close(cl_socket);
      return 1;
    }

    if (server_response == 0) {
      int ID;
      if (get_replica_data(replica_data, ID)) {
        close(cl_socket);
        return 1;
      }

      if (!quiet) {
        // std::cout << "slave: got a job! (replica "
        // << replica_data.ID << ")" << std::endl;

        std::cout << "slave:  running replica " << replica_data.ID
                << " at T=" << T[replica_data.Ti]
                << " and l=" << l[replica_data.li]
                << " (run = " << replica_data.run
                << " and dt = " << dt[replica_data.li] << ")"
                << std::endl;
      }

      // init replica parameters (t, conf, T, lambda)
      if (init_replica(topo, conf, sim, cg_topo, cg_conf, cg_sim)) {
        std::cerr << "init_replica returned error!" << std::endl;
        std::cerr << "slave: disconnecting..." << std::endl;
        close(cl_socket);
        return 1;
      }

      // close connection
      close(cl_socket);

      ////////////////////////////////////////////////////////////////////////////////
      // SLAVE: run the MD
      ////////////////////////////////////////////////////////////////////////////////

      if (!quiet)
        std::cerr << "running md" << std::endl;
      int error = run_md(topo, conf, sim, md, traj,
              cg_topo, cg_conf, cg_sim, cg_md, cg_traj);
      if (!quiet)
        std::cerr << "run finished" << std::endl;

      ////////////////////////////////////////////////////////////////////////////////
      // POSTPROCESSING
      ////////////////////////////////////////////////////////////////////////////////

      if (!error) {

        recalc_energy(topo, conf, sim, md, traj,
                cg_topo, cg_conf, cg_sim, cg_md, cg_traj);

        ++replica_data.run;
        replica_data.time = sim.time();
        replica_data.steps = sim.steps();
        replica_data.state = waiting;
      } else {
        replica_data.state = st_error;
      }

      DEBUG(8, "slave: connecting (after run)...");
      cl_socket = socket(addrinfo_p->ai_family, addrinfo_p->ai_socktype,
              addrinfo_p->ai_protocol);

      result = connect(cl_socket, s_addr_p, len);
      if (result == -1) {
        std::cout << "could not (re-)connect to master. master finished?"
                << std::endl;
        return 1;
      }

      DEBUG(9, "slave: connected");

      if (!magic_cookie(false)) {
        close(cl_socket);
        continue;
      }

      DEBUG(8, "slave: finished job " << replica_data.ID);
      char ch = 2;
      if (write(cl_socket, &ch, 1) != 1) {
        std::cerr << "could not write to socket" << std::endl;
        close(cl_socket);
        return 1;
      }

      // get ACK
      if (read(cl_socket, &ch, 1) != 1) {
        std::cerr << "could not read" << std::endl;
        close(cl_socket);
        return 1;
      }
      if (ch != 0) {
        io::messages.add("server reported error",
                "replica_exchange",
                io::message::error);
        close(cl_socket);
        return 1;
      }

      ////////////////////////////////////////////////////////////////////////////////
      // store the final configuration on server
      ////////////////////////////////////////////////////////////////////////////////

      if (put_replica_data(replica_data)) {
        close(cl_socket);
        return 1;
      }

      if (put_configuration(conf)) {
        std::cout << "could not update configuration!" << std::endl;
        return 1;
      }
      if (multigraining) {
        if (put_configuration(cg_conf)) {
          std::cout << "could not update coarse-grained configuration!" << std::endl;
          return 1;
        }
      }

      ////////////////////////////////////////////////////////////////////////////////
      // that's it
      ////////////////////////////////////////////////////////////////////////////////
      close(cl_socket);
    } else if (server_response == 9) {
      std::cout << "server has finished!\n"
              << "exiting...\n"
              << std::endl;

      std::cerr << "slave: disconnecting..." << std::endl;
      close(cl_socket);

      return 0;
    } else {
      close(cl_socket);
      if (!quiet)
        std::cout << "slave sleeping ..." << std::endl;
      // wait apropriate time for master to prepare
      sleep(timeout);
    }
  } // for slave_runs

  std::cerr << "slave: finished runs. terminating..." << std::endl;

  freeaddrinfo(addrinfo_p);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// init slave
////////////////////////////////////////////////////////////////////////////////

int util::Replica_Exchange_Slave::init
(
        io::Argument & args,
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        algorithm::Algorithm_Sequence & md,
        io::Out_Configuration & traj,
        topology::Topology & cg_topo,
        configuration::Configuration & cg_conf,
        simulation::Simulation & cg_sim,
        algorithm::Algorithm_Sequence & cg_md,
        io::Out_Configuration & cg_traj
        ) {
  // add an external interaction
  if (multigraining)
    sim.param().force.external_interaction = 1;

  std::cout << "reading (standard) input" << std::endl;
  io::read_input(args, topo, conf, sim, md, std::cout);
  md.init(topo, conf, sim, std::cout);

  traj.title(GROMOSXX "\n\tfine-grained\n" + sim.param().title);
  traj.init(args, sim.param());

  ////////////////////////////////////////////////////////////////////////////////
  // MULTIGRAINING
  ////////////////////////////////////////////////////////////////////////////////
  cg_traj.title(GROMOSXX "\n\tcoarse-grained\n" + cg_sim.param().title);

  if (multigraining) {

    std::cout << "\n\n"
            << "==================================================\n"
            << "MULTIGRAINING\n"
            << "==================================================\n"
            << std::endl;

    interaction::Forcefield * ff =
            dynamic_cast<interaction::Forcefield *> (md.algorithm("Forcefield"));
    if (ff == NULL) {
      std::cout << "Error: no forcefield in MD" << std::endl;
      return 1;
    }
    interaction::External_Interaction * ei =
            dynamic_cast<interaction::External_Interaction *> (ff->interaction("External"));
    if (ei == NULL) {
      std::cout << "Error: no external interaction in forcefield" << std::endl;
      return 1;
    }

    ei->set_coarsegraining(cg_topo, cg_conf, cg_sim);

    io::argname_conf = "cg_conf";
    io::argname_topo = "cg_topo";
    io::argname_pttopo = "cg_pttopo";
    io::argname_input = "cg_input";
    io::argname_trj = "cg_trc";
    io::argname_fin = "cg_fin";
    io::argname_tre = "cg_tre";

    io::read_input(args, cg_topo, cg_conf, cg_sim, cg_md);

    cg_md.init(cg_topo, cg_conf, cg_sim);

    cg_traj.init(args, sim.param());

  }

  // should be the same as from master... (more or less, multigraining)
  std::cout << "\nMESSAGES FROM (SLAVE) INITIALIZATION\n";
  if (io::messages.display(std::cout) >= io::message::error) {
    std::cout << "\nErrors during initialization!\n" << std::endl;
    return 1;
  }
  io::messages.clear();

  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

  std::cout << "\n\nslave initialised\n" << std::endl;
  if (multigraining) {
    std::cout << "\tmultigrained simulation" << std::endl;
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// RUN
////////////////////////////////////////////////////////////////////////////////

int util::Replica_Exchange_Slave::run_md
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        algorithm::Algorithm_Sequence & md,
        io::Out_Configuration & traj,
        topology::Topology & cg_topo,
        configuration::Configuration & cg_conf,
        simulation::Simulation & cg_sim,
        algorithm::Algorithm_Sequence & cg_md,
        io::Out_Configuration & cg_traj
        ) {
  // do we scale the initial temperatures?
  if (sim.param().replica.scale) {
    algorithm::Temperature_Calculation tcalc;
    tcalc.apply(topo, conf, sim);

    algorithm::Berendsen_Thermostat tcoup;
    tcoup.calc_scaling(topo, conf, sim, true);
    tcoup.scale(topo, conf, sim);

    if (multigraining) {
      tcalc.apply(cg_topo, cg_conf, cg_sim);
      tcoup.calc_scaling(cg_topo, cg_conf, cg_sim, true);
      tcoup.scale(cg_topo, cg_conf, cg_sim);
    }

  }

  int msteps = sim.param().multistep.steps;
  if (msteps < 1) msteps = 1;

  /*
  const double end_time = sim.time() + 
    sim.time_step_size() *
    (sim.param().step.number_of_steps * msteps - 1);
   */

  int error;

  const unsigned int end_steps =
          sim.steps() + sim.param().step.number_of_steps * msteps;

  // while(sim.time() < end_time + math::epsilon){
  while (sim.steps() < end_steps) {

    ////////////////////////////////////////////////////
    // multigraining
    // and
    // multiple time stepping
    ////////////////////////////////////////////////////
    if (multigraining) {

      // std::cout << "MULTISTEP: doing coarse-grained calculation\n";
      // coarse grained atom positions are based upon
      // real atom positions
      util::update_virtual_pos(cg_topo, cg_conf, topo, conf);

      // calculate the cg forces first!
      // if ((error = cg_ff->apply(cg_topo, cg_conf, cg_sim))){
      if ((error = cg_md.run(cg_topo, cg_conf, cg_sim))) {
        io::print_ENERGY(cg_traj.output(), cg_conf.current().energies,
                cg_topo.energy_groups(),
                "CGERROR", "CGERR_");

        io::print_ENERGY(cg_traj.output(), cg_conf.old().energies,
                cg_topo.energy_groups(),
                "CGOLDERROR", "CGOLDERR_");

        std::cout << "\nError during CG MD run!\n" << std::endl;
        break;
      }
    }

    // run a step
    if ((error = md.run(topo, conf, sim))) {

      io::print_ENERGY(traj.output(), conf.current().energies,
              topo.energy_groups(),
              "ERROR", "ERR_");

      io::print_ENERGY(traj.output(), conf.old().energies,
              topo.energy_groups(),
              "OLDERROR", "OLDERR_");

      std::cout << "\nError during MD run!\n" << std::endl;
      break;
    }

    traj.print(topo, conf, sim);

    sim.time() += sim.time_step_size();
    ++sim.steps();

    // the very first position will be missing...
    // still, this is much easier than any of the other
    // variants.
    // otherwise, the very last position is missing ;-)

    traj.write_replica_step(sim, replica_data);
    traj.write(conf, topo, sim, io::reduced);

    if (multigraining) {
      cg_sim.time() += cg_sim.time_step_size();
      ++cg_sim.steps();

      cg_traj.write_replica_step(cg_sim, replica_data);
      cg_traj.write(cg_conf, cg_topo, cg_sim, io::reduced);
    }

  }

  util::update_virtual_pos(cg_topo, cg_conf, topo, conf);

  if (!quiet) {
    std::cout << "\nslave: run finished!\n\n";
    std::cout << "\tlast energies were: potential=" << conf.old().energies.potential_total
            << "  special=" << conf.old().energies.special_total
            << "  kinetic=" << conf.old().energies.kinetic_total << "\n\n";
  }

  return error;
}

////////////////////////////////////////////////////////////////////////////////
// initialisation
////////////////////////////////////////////////////////////////////////////////

int util::Replica_Exchange_Slave::init_replica
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        topology::Topology & cg_topo,
        configuration::Configuration & cg_conf,
        simulation::Simulation & cg_sim
        ) {
  // get configuration from master
  if (get_configuration(conf)) {
    std::cerr << "get configuration failed" << std::endl;
    return 1;
  }

  if (multigraining) {
    if (get_configuration(cg_conf)) {
      std::cerr << "get configuration failed" << std::endl;
      return 1;
    }
  }

  std::vector<double> const & T = sim.param().replica.temperature;
  std::vector<double> const & l = sim.param().replica.lambda;

  // change all the temperature coupling temperatures
  for (unsigned int i = 0; i < sim.multibath().size(); ++i) {
    sim.multibath()[i].temperature = T[replica_data.Ti];
    // temperature constraining
    sim.multibath()[i].tau = sim.param().replica.dt[replica_data.li];

    if (multigraining) {
      cg_sim.multibath()[i].temperature = T[replica_data.Ti];
      // temperature constraining
      cg_sim.multibath()[i].tau = sim.param().replica.dt[replica_data.li];
    }
  }
  sim.param().stochastic.temp = T[replica_data.Ti];
  if (multigraining)
    cg_sim.param().stochastic.temp = T[replica_data.Ti];

  // change the lambda value
  sim.param().perturbation.lambda = l[replica_data.li];
  topo.lambda(l[replica_data.li]);
  // twice, to set old_lambda...
  topo.lambda(l[replica_data.li]);
  topo.update_for_lambda();

  if (!quiet)
    std::cout << "\tlambda = " << topo.lambda() << "\n";

  if (multigraining) {
    cg_sim.param().perturbation.lambda = l[replica_data.li];
    cg_topo.lambda(l[replica_data.li]);
    cg_topo.lambda(l[replica_data.li]);
    cg_topo.update_for_lambda();
  }

  // change simulation time
  sim.time() = replica_data.time;
  sim.steps() = replica_data.steps;

  // don't change time_step_size
  // sim.time_step_size() = sim.param().replica.dt[replica_data.li];
  // but enable multistepping
  sim.param().multistep.steps =
          int(rint(sim.param().replica.dt[replica_data.li] / sim.time_step_size()));

  if (!quiet)
    std::cout << "\tmultistepping steps = "
          << int(rint(sim.param().replica.dt[replica_data.li] / sim.time_step_size()))
    << "\n";

  if (multigraining) {
    cg_sim.time() = replica_data.time;
    cg_sim.param().multistep.steps =
            int(rint(sim.param().replica.dt[replica_data.li] / sim.time_step_size()));

    cg_sim.steps() = replica_data.steps;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// recalculate the energy after run
////////////////////////////////////////////////////////////////////////////////

int util::Replica_Exchange_Slave::recalc_energy
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        algorithm::Algorithm_Sequence & md,
        io::Out_Configuration & traj,
        topology::Topology & cg_topo,
        configuration::Configuration & cg_conf,
        simulation::Simulation & cg_sim,
        algorithm::Algorithm_Sequence & cg_md,
        io::Out_Configuration & cg_traj
        ) {
  // aliases
  std::vector<double> const & T = sim.param().replica.temperature;
  std::vector<double> const & l = sim.param().replica.lambda;
  // std::vector<double> const & dt = sim.param().replica.dt;

  int error = 0;

  if (!quiet)
    std::cout << "\tcalculating potential energies of last configuration\n\n";

  // do we need to reevaluate the potential energy ?
  // yes! 'cause otherwise it's just the energy for
  // the previous configuration...

  algorithm::Algorithm * cg_ff;

  if (multigraining) {
    // coarse grained atom positions are based upon
    // real atom positions
    util::update_virtual_pos(cg_topo, cg_conf, topo, conf);

    cg_ff = cg_md.algorithm("Forcefield");
    if (cg_ff == NULL) {
      std::cerr << "forcefield not found in MD algorithm sequence"
              << std::endl;
      return 1;
    }

    // calculate the cg forces first!
    if ((error = cg_ff->apply(cg_topo, cg_conf, cg_sim))) {
      cg_conf.current().energies.calculate_totals();
      io::print_ENERGY(traj.output(), cg_conf.current().energies,
              cg_topo.energy_groups(),
              "CGOLDERROR", "CGOLDERR_");

      io::print_ENERGY(traj.output(), cg_conf.old().energies,
              cg_topo.energy_groups(),
              "CGERROR", "CGERR_");

      std::cout << "\nError during CG final energy calc (" << error << ")!\n" << std::endl;
      return 1;
    }

    cg_conf.current().energies.calculate_totals();

    if (!quiet) {
      io::print_ENERGY(traj.output(), cg_conf.current().energies,
              cg_topo.energy_groups(),
              "CGENERGY_Li", "CGELi_");
    }

  } // multigraining

  algorithm::Algorithm * ff = md.algorithm("Forcefield");

  if (ff == NULL) {
    std::cerr << "forcefield not found in MD algorithm sequence"
            << std::endl;
    return 1;
  }

  if (ff->apply(topo, conf, sim)) {
    io::print_ENERGY(traj.output(), conf.current().energies,
            topo.energy_groups(),
            "ERROR", "ERR_");

    io::print_ENERGY(traj.output(), conf.old().energies,
            topo.energy_groups(),
            "OLDERROR", "OLDERR_");

    std::cout << "\nError during final energy calc!\n" << std::endl;
    return 1;
  }

  conf.current().energies.calculate_totals();
  replica_data.epot_i = conf.current().energies.potential_total +
          conf.current().energies.special_total;

  if (!quiet) {
    io::print_ENERGY(traj.output(), conf.current().energies,
            topo.energy_groups(),
            "ENERGY_Li", "ELi_");
  }

  traj.write_replica_energy(replica_data, sim,
          conf.current().energies,
          0);

  if (!quiet) {
    if (replica_data.Ti != replica_data.Tj)
      std::cout << "replica_energy " << replica_data.epot_i
            << " @ " << T[replica_data.Ti] << "K" << std::endl;
    else if (replica_data.li != replica_data.lj)
      std::cout << "replica_energy " << replica_data.epot_i
            << " @ " << l[replica_data.li] << std::endl;
  }

  if (replica_data.li != replica_data.lj) {

    if (!quiet)
      std::cout << "energies at switched lambda:\n\n";

    // change the lambda value
    sim.param().perturbation.lambda = l[replica_data.lj];
    topo.lambda(l[replica_data.lj]);
    // twice: set old lambda as well
    topo.lambda(l[replica_data.lj]);
    topo.update_for_lambda();
    if (!quiet)
      std::cout << "\tlambda = " << topo.lambda() << "\n";

    if (multigraining) {
      cg_sim.param().perturbation.lambda = l[replica_data.lj];
      cg_topo.lambda(l[replica_data.lj]);
      cg_topo.lambda(l[replica_data.lj]);
      cg_topo.update_for_lambda();

      if ((error = cg_ff->apply(cg_topo, cg_conf, cg_sim))) {
        cg_conf.current().energies.calculate_totals();
        io::print_ENERGY(traj.output(), cg_conf.current().energies,
                cg_topo.energy_groups(),
                "CGOLDERROR", "CGOLDERR_");

        io::print_ENERGY(traj.output(), cg_conf.old().energies,
                cg_topo.energy_groups(),
                "CGERROR", "CGERR_");

        std::cout << "\nError during CG force recalc at different lambda!\n"
                << std::endl;
        return 1;
      }

      cg_conf.current().energies.calculate_totals();
      if (!quiet) {
        io::print_ENERGY(traj.output(), cg_conf.current().energies,
                cg_topo.energy_groups(),
                "CGENERGY_Lj", "CGELj_");
      }
    }

    // recalc energy
    if (ff->apply(topo, conf, sim)) {
      io::print_ENERGY(traj.output(), conf.current().energies,
              topo.energy_groups(),
              "OLDERROR", "OLDERR_");

      io::print_ENERGY(traj.output(), conf.old().energies,
              topo.energy_groups(),
              "ERROR", "ERR_");

      std::cout << "\nError during force recalc at different lambda!\n" << std::endl;
    }

    conf.current().energies.calculate_totals();
    replica_data.epot_j = conf.current().energies.potential_total +
            conf.current().energies.special_total;

    if (!quiet) {
      io::print_ENERGY(traj.output(), conf.current().energies,
              topo.energy_groups(),
              "ENERGY_Lj", "ELj_");
    }

    traj.write_replica_energy(replica_data, sim,
            conf.current().energies,
            1);

  }
  else {
    replica_data.epot_j = replica_data.epot_i;
  }
  return 0;
}

#endif
