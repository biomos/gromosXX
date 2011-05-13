/**
 * @file replica_exchange_master.cc
 * replica exchange
 */

#include "../stdheader.h"

#ifdef REPEX

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../interaction/interaction.h"
#include "../interaction/forcefield/forcefield.h"

#include "../algorithm/temperature/temperature_calculation.h"
#include "../algorithm/temperature/berendsen_thermostat.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/error.h"

#include "../io/read_input.h"
#include "../io/print_block.h"

#include <time.h>
#include <unistd.h>

#include "../io/configuration/out_configuration.h"

#include "../math/volume.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

#include "replica_exchange.h"


#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica

////////////////////////////////////////////////////////////////////////////////
// replica master   ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

util::Replica_Exchange_Master::Replica_Exchange_Master()
: switch_T(0),
switch_l(0) {
  DEBUG(7, "creating replica master");

  // enable control via environment variables
  gsl_rng_env_setup();
  const gsl_rng_type * rng_type = gsl_rng_default;

  // get the random number generator
  m_rng = gsl_rng_alloc(rng_type);
}

int util::Replica_Exchange_Master::run
(
        io::Argument & args) {
  ////////////////////////////////////////////////////////////
  // socket
  sockaddr_in server_address;
  sockaddr_in client_address;
  serv_socket = socket(AF_INET, SOCK_STREAM, 0);

  server_address.sin_family = AF_INET;
  server_address.sin_addr.s_addr = htonl(INADDR_ANY);

  int port_nr = 29375;
  if (args.count("master") > 0) {
    std::istringstream is(args["master"]);
    if (!(is >> port_nr)) {
      io::messages.add("master [port number]",
              "replica_exchange",
              io::message::error);
      return 1;
    }
  }
  server_address.sin_port = htons(port_nr);
  socklen_t server_len = sizeof (server_address);

  if (bind(serv_socket, (sockaddr *) & server_address, server_len)) {

    io::messages.add("could not bind address",
            "replica_exchange",
            io::message::error);

    std::cout << "replica exchange: could not bind address (error = "
            << errno << ")" << std::endl;

    close(serv_socket);
    return 1;
  }

  multigraining = false;
  if (args.count("cg_topo") >= 0) {
    multigraining = true;
  }

  // create the simulation classes (necessary to store the configurations)
  topology::Topology topo;
  simulation::Simulation sim;
  algorithm::Algorithm_Sequence md;

  // read the files
  if (io::read_replica_input(args, topo, m_conf, sim,
          md, replica_data, std::cout)) {
    std::cerr << "could not read input!!!" << std::endl;
    io::messages.add("replica exchange: could not read input!",
            "replica exchange",
            io::message::critical);
    close(serv_socket);
    return 2;
  }

  // save the initial run information
  init_run.resize(replica_data.size(), 0);
  for(unsigned int i = 0; i < replica_data.size(); ++i) {
    init_run[i] = replica_data[i].run;
  }


  // write whenever we want!
  sim.param().write.position = 1;
  if (sim.param().xrayrest.xrayrest != simulation::xrayrest_off &&
          sim.param().xrayrest.write != 0) {
    sim.param().xrayrest.write = 1;
  } else {
    sim.param().xrayrest.write = 0;
  }

  // initialises everything
  for (unsigned int i = 0; i < m_conf.size(); ++i) {
    if (md.init(topo, m_conf[i], sim, std::cout)) {
      std::cerr << "md init failed!" << std::endl;
      close(serv_socket);
      return 3;
    }
  }

  io::Out_Configuration traj(GROMOSXX "\n\treplica master\n\tfine-grained\n", std::cout);
  traj.title(GROMOSXX "\n\treplica master\n\tfine-grained\n" + sim.param().title);
  traj.init(args, sim.param());

  ////////////////////////////////////////////////////////////////////////////////
  // multigraining
  ////////////////////////////////////////////////////////////////////////////////
  topology::Topology cg_topo;
  simulation::Simulation cg_sim;
  algorithm::Algorithm_Sequence cg_md;
  io::Out_Configuration cg_traj(GROMOSXX "\n\treplica master\n\tcoarse-grained\n", std::cout);
  cg_traj.title(GROMOSXX "\n\treplica master\n\tcoarse-grained\n" + sim.param().title);

  if (multigraining) {

    std::cout << "MULTIGRAINING\n\n";

    io::argname_conf = "cg_conf";
    io::argname_topo = "cg_topo";
    io::argname_pttopo = "cg_pttopo";
    io::argname_input = "cg_input";
    io::argname_trj = "cg_trc";
    io::argname_fin = "cg_fin";
    io::argname_tre = "cg_tre";

    m_cg_conf.resize(m_conf.size());

    std::vector<Replica_Data> rep_data;
    if (io::read_replica_input(args, cg_topo, m_cg_conf, cg_sim, cg_md,
            rep_data, std::cout)) {
      std::cerr << "could not read coarse-grained input!!!" << std::endl;
      io::messages.add("replica exchange: could not read coarse-grained input!",
              "replica exchange",
              io::message::critical);
      close(serv_socket);
      return 2;
    }


    for (unsigned int i = 0; i < m_cg_conf.size(); ++i) {
      if (cg_md.init(cg_topo, m_cg_conf[i], cg_sim, std::cout)) {
        std::cerr << "md init failed (coarse-grained)!" << std::endl;
        close(serv_socket);
        return 3;
      }
    }

    cg_sim.param().write.position = 1;

    cg_traj.init(args, cg_sim.param());

    std::cout << "END\n";
  }
  ////////////////////////////////////////////////////////////////////////////////
  // end multigraining input
  ////////////////////////////////////////////////////////////////////////////////

  // check input
  if (sim.param().replica.num_T * sim.param().replica.num_l < 1) {
    io::messages.add("replica exchange with no replicas?!",
            "Replica_Exchange",
            io::message::error);
    std::cerr << "RE: not enough replicas" << std::endl;
    close(serv_socket);
    return 4;
  }

  if (sim.param().replica.trials < 1) {
    io::messages.add("replica exchange with less than 1 trials?!",
            "Replica_Exchange",
            io::message::error);
  }

  switch_T = sim.param().replica.num_T;
  switch_l = sim.param().replica.num_l;

  const int rep_num = switch_T * switch_l;

  std::cout << "\nMESSAGES FROM (MASTER) INITIALIZATION\n";
  if (io::messages.display(std::cout) >= io::message::error) {
    std::cout << "\nErrors during initialization!\n" << std::endl;
    std::cerr << "RE: not enough trials" << std::endl;
    return 4;
  }

  io::messages.display();
  io::messages.clear();

  // set seed if not set by environment variable
  if (gsl_rng_default_seed == 0)
    gsl_rng_set(m_rng, sim.param().start.ig);

  std::cout << "master thread initialised" << std::endl;

  std::cout << "Replica Exchange\n"
          << "\treplicas (temperature) :\t" << switch_T << "\n"
          << "\treplicas (lambda)      :\t" << switch_l << "\n"
          << "\treplicas (total)       :\t" << rep_num << "\n\t"
          << std::setw(10) << "ID"
          << std::setw(20) << "Temp"
          << std::setw(20) << "lambda\n";

  {
    int i = 0;
    for (int l = 0; l < switch_l; ++l) {
      for (int t = 0; t < switch_T; ++t, ++i) {

        std::cout << "\t"
                << std::setw(10) << i + 1
                << std::setw(20)
                << sim.param().replica.temperature[replica_data[i].Ti]
                << std::setw(20)
                << sim.param().replica.lambda[replica_data[i].li]
                << "\n";
      }
    }
  }

  std::cout << "\n\ttrials       :\t" << sim.param().replica.trials
          << "\n\truns (slave) :\t" << sim.param().replica.slave_runs
          << "\n\nEND\n" << std::endl;

  trials = 1;
  runs = 0;

  rep_out.open("replica.dat");
  rep_out << "num_T\t" << switch_T << "\n"
          << "num_l\t" << switch_l << "\n";

  rep_out.precision(4);
  rep_out.setf(std::ios::fixed, std::ios::floatfield);

  rep_out << "T    \t";
  for (int t = 0; t < switch_T; ++t)
    rep_out << std::setw(12) << sim.param().replica.temperature[t];

  rep_out << "\nl    \t";
  for (int l = 0; l < switch_l; ++l)
    rep_out << std::setw(12) << sim.param().replica.lambda[l];

  rep_out << "\n\n";

  rep_out << "#"
          << std::setw(5) << "ID"
          << std::setw(7) << "run"
          << std::setw(13) << "Ti"
          << std::setw(13) << "li"
          << std::setw(14) << "Epoti"
          << std::setw(13) << "Tj"
          << std::setw(13) << "lj"
          << std::setw(14) << "Epotj"
          << std::setw(13) << "p"
          << std::setw(4) << "s"
          << "\n";

  if (!quiet) {
    std::cout << std::setw(6) << "ID"
            << std::setw(7) << "run"
            << std::setw(13) << "Ti"
            << std::setw(13) << "li"
            << std::setw(14) << "Epoti"
            << std::setw(13) << "Tj"
            << std::setw(13) << "lj"
            << std::setw(14) << "Epotj"
            << std::setw(13) << "p"
            << std::setw(4) << "s"
            << std::endl;
  }

  bool quit = false;

  // listen, keep a queue for max all replicas
  // that should be enough, i guess
  if (listen(serv_socket, switch_l * switch_T)) {

    io::messages.add("could not listen on address", "replica_exchange",
            io::message::error);

    std::cout << "replica exchange: could not listen" << std::endl;
    close(serv_socket);
    return 1;
  }

  while (quit == false) {

    if (runs == rep_num) {
      // replicas run depth-first. so some replicas are further than others
      // the average finished the trial...

      // write recovery trajectory
      if (sim.param().replica.write &&
              ((trials % sim.param().replica.write) == 0)) {
        std::cout << "writing trajectory..." << std::endl;
        traj.write_replica(replica_data, m_conf, topo, sim);
        if (multigraining)
          cg_traj.write_replica(replica_data, m_cg_conf, cg_topo, cg_sim);
      }

      ++trials;
      ++sim.steps();
      sim.time() += sim.param().step.number_of_steps * sim.param().step.dt;
      runs = 0;

      if (!quiet) {
        std::cout
                << "\n======================================================="
                << "==============================================================\n"
                << std::setw(6) << "ID"
                << std::setw(7) << "run"
                << std::setw(13) << "T"
                << std::setw(13) << "l"
                << std::setw(19) << "Epot"
                << std::setw(13) << "sT"
                << std::setw(13) << "sl"
                << std::setw(19) << "sEpot"
                << std::setw(13) << "p"
                << std::setw(4) << "s"
                << "\n-------------------------------------------------------"
                << "--------------------------------------------------------------"
                << std::endl;
      }

      rep_out.flush();
    }

    if (trials > sim.param().replica.trials) {
      DEBUG(8, "master: finished all trials...");
      rep_out.flush();
      std::cout << "master: finished with all trials...\n";
      // keep around 'till slaves have finished
      // break;
    }

    // wait for a thread to connect
    DEBUG(9, "master: accepting connection...");
    socklen_t client_len = sizeof (client_address);
    cl_socket = accept(serv_socket, (struct sockaddr *) &client_address, &client_len);
    if (cl_socket == -1) {
      std::cout << "accept failed!" << std::endl;
      io::messages.add("accept failed",
              "replica_exchange",
              io::message::error);
      close(serv_socket);
      return 4;
    }

    if (!magic_cookie(true)) {
      close(cl_socket);
      continue;
    }

    char ch;
    if (read(cl_socket, &ch, 1) != 1) {
      std::cerr << "could not read from socket!" << std::endl;
      close(cl_socket);
      continue;
    }

    DEBUG(9, "client connected!");

    switch (ch) {
      case 0:
      {
        std::cout << "master: got a 'hello'\n";
        break;
      }

      case 1:
      {
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        // client asks for job
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        select_job(sim);

        close(cl_socket);
        break;
      }

      case 2:
      {
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        // client finished job
        //////////////////////////////////////////////////////////////////////////////////////////////////////

        finish_job(sim);

        close(cl_socket);
        break;
      }

      case 3:
      {
        ch = 0;
        if (write(cl_socket, &ch, 1) != 1) {
          std::cerr << "could not write to socket" << std::endl;
          close(cl_socket);
          return 1;
        }

        std::cout << "process " << ch << " has aborted run\n";
        break;
      }

      case 4: // interactive session
      {

        if (interactive(sim) == 2)
          quit = true;
        close(cl_socket);
        break;

        default:
        std::cout << "message not understood\n";
      }
    }

    DEBUG(9, "disconnecting");
    close(cl_socket);

  } // while trials to do

  // write out final configurations
  traj.write_replica(replica_data, m_conf, topo, sim, io::final);

  if (multigraining)
    cg_traj.write_replica(replica_data, m_cg_conf, cg_topo, cg_sim, io::final);

  // simulation done
  DEBUG(9, "master: done");

  rep_out.flush();
  rep_out.close();

  std::cout << "exiting..." << std::endl;
  close(serv_socket);

  return 0;
}

int util::Replica_Exchange_Master::select_job(simulation::Simulation & sim) {

  char ch;

  if (trials > sim.param().replica.trials) { // terminate...
    std::cout << "master: sending quit (9) signal..." << std::endl;
    // quit signal
    ch = 9;
    //int res = write(cl_socket, &ch, 1);
    return 1;
  }

  const int rep_num = switch_T * switch_l;

  // select a replica to run
  DEBUG(9, "request a job");
  // std::cout << "requesting a job" << std::endl;

  
  int select_replica = -1;

  for (int r = 0; r < rep_num; ++r) {
    int min_run = sim.param().replica.trials + init_run[r];
    if (replica_data[r].state == waiting) {
      // try a switch
      switch_replica(r, sim.param());
    }

    if (replica_data[r].state == ready &&
            replica_data[r].run < min_run) {
      select_replica = r;
      break;
    }
  }

  if (select_replica == -1) {
    // no replica available, wait...
    std::cout << "could not select replica!!!" << std::endl;
    ch = 1;
    if (write(cl_socket, &ch, 1) != 1) {
      std::cerr << "could not set job to waiting" << std::endl;
      return 1;
    }
  } else {

    // assign it!
    replica_data[select_replica].state = running;

    // all ok, sending replica
    ch = 0;
    if (write(cl_socket, &ch, 1) != 1) {
      std::cerr << "could not write to socket" << std::endl;
      return 1;
    }

    if (put_replica_data(replica_data[select_replica])) {
      std::cerr << "could not write replica information" << std::endl;
      return 1;
    }

    if (put_configuration(m_conf[select_replica])) {
      std::cerr << "could not write configuration" << std::endl;
      // try again
      replica_data[select_replica].state = ready;
      return 1;
    }

    if (multigraining) {

      if (put_configuration(m_cg_conf[select_replica])) {
        std::cerr << "could not write configuration" << std::endl;
        // try again
        replica_data[select_replica].state = ready;
        return 1;
      }
    } // multigraining
  }

  return 0;
}

int util::Replica_Exchange_Master::finish_job(simulation::Simulation & sim) {
  char ch = 0;

  if (write(cl_socket, &ch, 1) != 1) {
    std::cerr << "could not wirte to socket" << std::endl;
    return 1;
  }

  // make a backup copy
  // (if transaction fails)
  std::vector<Replica_Data> tmp_rep = replica_data;

  int i;
  try {
    get_replica_data(replica_data, i);
  }  catch (std::runtime_error e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    std::cout << "Exception: " << e.what() << std::endl;
    replica_data = tmp_rep;
    if (i >= 0) replica_data[i].state = ready;
    return 1;
  }

  if (replica_data[i].state != st_error) {

    // backup pos / vel / box
    math::VArray tmp_pos = m_conf[i].current().pos;
    math::VArray tmp_vel = m_conf[i].current().vel;
    math::Box tmp_box = m_conf[i].current().box;

    if (get_configuration(m_conf[i])) {
      // restore original replica data
      replica_data = tmp_rep;
      // restore pos / vel / box
      m_conf[i].current().pos = tmp_pos;
      m_conf[i].current().vel = tmp_vel;
      m_conf[i].current().box = tmp_box;

      replica_data[i].state = ready;
      return 1;
    }

    if (multigraining) {

      // backup pos / vel / box
      math::VArray cg_tmp_pos = m_cg_conf[i].current().pos;
      math::VArray cg_tmp_vel = m_cg_conf[i].current().vel;
      math::Box cg_tmp_box = m_cg_conf[i].current().box;

      if (get_configuration(m_cg_conf[i])) {

        // restore original replica data
        replica_data = tmp_rep;
        // restore pos / vel / box
        m_conf[i].current().pos = tmp_pos;
        m_conf[i].current().vel = tmp_vel;
        m_conf[i].current().box = tmp_box;

        m_cg_conf[i].current().pos = cg_tmp_pos;
        m_cg_conf[i].current().vel = cg_tmp_vel;
        m_cg_conf[i].current().box = cg_tmp_box;

        replica_data[i].state = ready;

        return 1;
      }
    } // multigraining

    DEBUG(9, "master: got replica " << replica_data[i].ID
            << " temperature=" << replica_data[i].Ti
            << " lambda=" << replica_data[i].li);

    if (replica_data[i].state == waiting) { // should always be true (after run) ?!
      // try a switch
      switch_replica(i, sim.param());
    }

    ++runs;
  } else { // state error
    std::cout << "received replica " << i << " with state error!" << std::endl;
    std::cout << "ID = " << replica_data[i].ID
            << "\ntemperature = " << replica_data[i].Ti
            << "\nlambda = " << replica_data[i].li
            << "\n" << std::endl;
  }

  return 0;
}

int util::Replica_Exchange_Master::interactive(simulation::Simulation & sim) {
  char ch = 0;

  if (write(cl_socket, &ch, 1) != 1) {
    std::cerr << "could not write to socket" << std::endl;
    close(cl_socket);
    return 1;
  }

  if (read(cl_socket, &ch, 1) != 1) {
    std::cerr << "could not read from socket" << std::endl;
    close(cl_socket);
    return 1;
  }

  switch (ch) {
    case 1: // replica information
    {
      int sz = replica_data.size();
      if (write(cl_socket, (char *) &sz, sizeof (int)) != sizeof (int)) {
        std::cerr << "could not write replica id" << std::endl;
        return 1;
      }

      if (write(cl_socket, (char *) &replica_data[0],
              sz * sizeof (Replica_Data)) != sz * int(sizeof (Replica_Data))) {
        std::cerr << "could not write replica data" << std::endl;
        return 1;
      }
      return 0;
    }
    case 2: // replica change
    {
      int sz = replica_data.size();
      if (write(cl_socket, (char *) &sz, sizeof (int)) != sizeof (int)) {
        std::cerr << "could not write replica ID" << std::endl;
        return 1;
      }

      // std::cerr << "sending data" << std::endl;
      if (write(cl_socket, (char *) &replica_data[0],
              sz * sizeof (Replica_Data)) != sz * int(sizeof (Replica_Data))) {
        std::cerr << "could not write replica data" << std::endl;
        return 1;
      }

      int r;
      // std::cerr << "reading ID" << std::endl;
      if (read(cl_socket, (char *) &r, sizeof (int)) != sizeof (int)) {
        std::cerr << "could not read replica ID" << std::endl;
        return 1;
      }

      if (r < 0 || r >= int(replica_data.size())) {
        io::messages.add("replica ID out of range",
                "replica_exchange",
                io::message::error);
        return 1;
      }

      if (read(cl_socket, (char *) &replica_data[r], sizeof (Replica_Data))
              != sizeof (Replica_Data)) {
        std::cerr << "could not read replica data!" << std::endl;
      }
      return 0;
    }
    case 3: // quit
    {
      std::cout << "master: stopping" << std::endl;
      return 2;
    }
  }
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// switching methods
////////////////////////////////////////////////////////////////////////////////

int util::Replica_Exchange_Master::switch_replica(int i, simulation::Parameter const & param) {
  if (i < 0 || unsigned(i) >= replica_data.size()) {
    std::cerr << "trying to switch replica: index out of range (" << i << ")" << std::endl;
    return 0;
  }

  assert(i >= 0 && unsigned(i) < replica_data.size());
  DEBUG(8, "switch replica: " << i);

  // aliases
  std::vector<double> const & T = param.replica.temperature;
  std::vector<double> const & l = param.replica.lambda;

  if (replica_data[i].state != waiting) {
    assert(false);
    return i;
  }

  const int j = find_switch_partner(i);
  if (j == -1) {
    // partner not yet ready...
    return 0;
  }
  if (j < 0 || unsigned(j) >= replica_data.size()) {
    std::cerr << "trying to switch replica: index j out of range (" << j << ")" << std::endl;
    return 0;
  }

  if (j == i) { // no switch this time...
    replica_data[i].probability = 0.0;
    replica_data[i].switched = false;

    if (!quiet)
      print_replica(i, param, std::cout);
    print_replica(i, param, rep_out);

    set_next_switch(i);
    replica_data[i].state = ready;

    return 0;
  }

  // try switch...
  double probability = switch_probability(i, j, param);

  // equilibrate starting structures w/o switching
  if (replica_data[i].run < param.replica.equilibrate)
    probability = 0.0;

  replica_data[i].probability = probability;
  replica_data[j].probability = probability;
  replica_data[i].switched = false;
  replica_data[j].switched = false;

  const double r = gsl_rng_uniform(m_rng);

  if (r < probability) {
    // SUCCEEDED!!!

    if (replica_data[i].Ti != replica_data[j].Ti) {
      if (!quiet) {
        std::cout << "-----> switch: "
                << T[replica_data[i].Ti]
                << " <-> "
                << T[replica_data[j].Ti]
                << "\n";
      }
    } else {
      if (!quiet) {
        std::cout << "-----> switch: "
                << l[replica_data[i].li]
                << " <-> "
                << l[replica_data[j].li]
                << "\n";
      }
    }

    replica_data[i].switched = true;
    replica_data[j].switched = true;
  }

  if (!quiet)
    print_replica(i, param, std::cout);
  print_replica(i, param, rep_out);

  if (!quiet)
    print_replica(j, param, std::cout);
  print_replica(j, param, rep_out);

  if (replica_data[i].switched) {

    replica_data[i].li = replica_data[i].lj;
    replica_data[i].Ti = replica_data[i].Tj;

    replica_data[j].li = replica_data[j].lj;
    replica_data[j].Ti = replica_data[j].Tj;
  }

  set_next_switch(i);
  set_next_switch(j);

  replica_data[i].state = ready;
  replica_data[j].state = ready;

  return 0;
}

double util::Replica_Exchange_Master::switch_probability
(
        int i, int j,
        simulation::Parameter const & param
        ) {
  // aliases
  std::vector<double> const & T = param.replica.temperature;
  // std::vector<double> const & l = param.replica.lambda;

  double delta = 0;
  const double bi = 1.0 / (math::k_Boltzmann * T[replica_data[i].Ti]);
  const double bj = 1.0 / (math::k_Boltzmann * T[replica_data[j].Ti]);

  if (replica_data[i].li != replica_data[j].li) {
    // 2D formula
    delta =
            bi * (replica_data[j].epot_j - replica_data[i].epot_i) -
            bj * (replica_data[j].epot_i - replica_data[i].epot_j);
    /*
      std::cout << "i=" << i << " j=" << j << " bi=" << bi << "  bj=" << bj
      << "  [i]epot_i=" << replica_data[i].epot_i
      << "  [j]epot_i=" << replica_data[j].epot_i
      << "  [i]epot_j=" << replica_data[i].epot_j
      << "  [j]epot_j=" << replica_data[j].epot_j
      << "\n";*/
  } else {
    // standard formula
    delta =
            (bi - bj) *
            (replica_data[j].epot_i - replica_data[i].epot_i);
  }

  // and pressure coupling
  if (param.pcouple.scale != math::pcouple_off) {
    delta += (bi - bj) *
            (param.pcouple.pres0(0, 0) + param.pcouple.pres0(1, 1) + param.pcouple.pres0(2, 2)) / 3.0 *
            (math::volume(m_conf[j].current().box, m_conf[j].boundary_type) -
            math::volume(m_conf[i].current().box, m_conf[i].boundary_type));
  }

  if (!quiet)
    std::cout << "\tswitching: delta = " << std::setw(18) << delta << "\n";

  double probability = 1.0;
  if (delta > 0.0)
    probability = exp(-delta);

  return probability;
}

int util::Replica_Exchange_Master::find_switch_partner(int i) {
  assert(i >= 0 && unsigned(i) < replica_data.size());
  DEBUG(8, "find switch partner of replica: " << i);

  for (unsigned int j = 0; j < replica_data.size(); ++j) {

    if (replica_data[j].state != waiting) continue;
    if (replica_data[j].run != replica_data[i].run) continue;

    if (replica_data[j].Ti == replica_data[i].Tj &&
            replica_data[j].li == replica_data[i].lj)
      return j;
  }

  return -1;
}

void util::Replica_Exchange_Master::set_next_switch(int i) {
  int l_change = 0, T_change = 0;

  if (switch_l > 1 && switch_T > 1) {
    const int c = replica_data[i].run % 4;
    switch (c) {
      case 0: T_change = 1;
        break;
      case 1: l_change = 1;
        break;
      case 2: T_change = -1;
        break;
      case 3: l_change = -1;
        break;
    }
  } else if (switch_T > 1) {
    const int c = replica_data[i].run % 2;
    switch (c) {
      case 0: T_change = 1;
        break;
      case 1: T_change = -1;
        break;
    }
  } else if (switch_l > 1) {
    const int c = replica_data[i].run % 2;
    switch (c) {
      case 0: l_change = 1;
        break;
      case 1: l_change = -1;
        break;
    }
  } else {
    std::cout << "Only one replica: why are you running replica exchange???" << std::endl;
    /*
      io::messages.add("No exchanges in replica exchange?",
      "Replica Exchange",
      io::message::critical);
     */
  }

  // and the modifiers
  if ((replica_data[i].Ti % 2) == 1) T_change = -T_change;
  if ((replica_data[i].li % 2) == 1) l_change = -l_change;

  replica_data[i].Tj = replica_data[i].Ti + T_change;
  replica_data[i].lj = replica_data[i].li + l_change;

  // check if at the edge...
  if (replica_data[i].Tj < 0 || replica_data[i].Tj >= switch_T)
    replica_data[i].Tj = replica_data[i].Ti;
  if (replica_data[i].lj < 0 || replica_data[i].lj >= switch_l)
    replica_data[i].lj = replica_data[i].li;
}


////////////////////////////////////////////////////////////////////////////////
// output
////////////////////////////////////////////////////////////////////////////////

void util::Replica_Exchange_Master::print_replica(int r,
        simulation::Parameter const & param,
        std::ostream & os) {
  os << std::setw(6) << r + 1
          << " "
          << std::setw(6) << replica_data[r].run
          << std::setw(13) << param.replica.temperature[replica_data[r].Ti]
          << std::setw(13) << param.replica.lambda[replica_data[r].li]
          << " "
          << std::setw(18) << replica_data[r].epot_i
          << std::setw(13) << param.replica.temperature[replica_data[r].Tj]
          << std::setw(13) << param.replica.lambda[replica_data[r].lj]
          << " "
          << std::setw(18) << replica_data[r].epot_j
          << std::setw(13) << replica_data[r].probability
          << std::setw(4) << replica_data[r].switched
          << std::endl;
}

#endif
