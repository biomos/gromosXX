/**
 * @file replica.h
 * Holds all the information for a single replica
 */
#include <stdheader.h>
#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>
#include <io/configuration/out_configuration.h>
#include <io/argument.h>
#include <io/read_input.h>
#include <math/random.h>
#include <util/repex_mpi.h>
#include <util/usage.h>

#ifndef REPLICA_H
#define REPLICA_H


namespace util {

  /**
   * replica Class
   * All the data members are held public, so Master/Slave classes have full 
   * control over it.
   */
  class replica {
  private:
    /**
     * copy constructor
     * don't allow copies because there's a bug in the copy constructor of configuration::configurtiaon
     */
    replica(const replica & r);
  public:
    /**
     * Constructor
     * Reads in all input parameters and assigns rank and ID of the replica
     * @param _args io::Argument, copy of the one initialized in repex_mpi.cc 
     * @param _ID integer, unique identifier of the replica
     * @param _rank integer, for knowing where to send the data
     */
    replica(io::Argument _args, int cont, int _ID, int _rank);

    /**
     * Deconstructor
     */
    ~replica();

    /**
     * runs an MD simulation
     */
    void run_MD();
    /**
     * init an MD simulation
     */
    void init();
    /**
     * Finds partner for current switch
     * @return ID of partner, own ID if no switching in current trial
     */
    int find_partner() const;
    /**
     * calculates probability of switch with current partner, may involve MPI communication
     * @param partnerID
     * @param partnerRank
     * @return probability
     */
    double calc_probability(const int partnerID, const int partnerRank);
    /**
     * calculates potential energy of current configuration with lambda(Hamiltonian) of partner
     * @param partner ID of partner
     * @return potential energy of configuration with lambda(Hamiltonian) of partner
     */
    double calculate_energy(const int partner);
    /**
     * calculates potential energy for current configuration with own lambda
     */
    double calculate_energy();
    /**
     * finds out if configurations are to be switched; sets switched to true if so
     */
    void swap(const unsigned int partner, const unsigned int partnerRank);
    
    /**
     * switch back the averages
     */
     void exchange_averages();
     
    /**
     * write final cnf
     */
     void write_final_conf();
     
    /**
     * Initiates MPI communication to receive new configuration information
     * @param senderID
     * @param senderRank
     */
    void receive_new_coord(const int senderID, const int senderRank);
    /**
     * Initiates MPI communication to send new configuration information
     * @param receiverID
     * @param receiverRank
     */
    void send_coord(const int receiverID, const int senderRank);
    /**
     * Prints information of replica to std::cout for debugging purposes
     */
    void print_info(std::string bla)const;
    /**
     * Scales the velocities after an exchange in temperature replica exchange
     */
    void velscale(int i);

  private:
    /**
     * Sets lambda parameter to original value of replica
     */
    void set_lambda();
    /**
     * Sets temperature to original value of replica
     */
    void set_temp();
    /**
     * Sets lambda parameter to value of partner (for energy calculations)
     */
    void change_lambda(const unsigned int partner);
    /**
     * Sets temperature to value of partner (for energy calculations)
     */
    void change_temp(const unsigned int partner);

  public:
    /**
     * ID, unique identifier
     */
    const unsigned int ID;
    /**
     * rank of node the replica is on
     */
    const unsigned int rank;

    /* NOTE:
     * These cannot be references, because every node has only one copy, so
     * the repicas would then reference the same conf/topo
     */
    /**
     * io::Argument object to read in all simulation parameters
     */
    io::Argument args;
    /**
     * trajectories
     */
    io::Out_Configuration * traj;
    /**
     * output file stream for out_configuration
     */
    std::ofstream * os;
    /**
     * configuration of replica
     */
    configuration::Configuration conf;
    /**
     * a temporary VArray for switching the lattice shifts
     */
    math::VArray * latticeTMP;
    /**
     * algorithm sequence
     */
    algorithm::Algorithm_Sequence md;
    /**
     * simulation
     */
    simulation::Simulation sim;
    /**
     * topology
     */
    topology::Topology topo;
    /**
     * random number generator
     */
    math::RandomGeneratorGSL * rng;
    /**
     * Temperature of replica
     */
    double T;
    /**
     * Lambda of replica
     */
    double l;
    /**
     * Timestep for current run
     */
    double dt;
    /**
     * potential energy using current Hamiltonian
     */
    double epot;
    /**
     * potential energy of the partner Hamiltonian (for bookkeeping)
     */
    double epot_partner;
    /**
     * ID of partner
     */
    unsigned int partner;

    /**
     * run number of replica
     */
    unsigned int run;
    
    /**
     * total number of runs per replica
     */
    unsigned int total_runs;

    /**
     * total number of steps simulated
     */
    unsigned int steps;
    /**
     * maximal number of steps to simulate per run
     */
    unsigned int maxSteps;

    /**
     * probability of last switch
     */
    double probability;
    /**
     * switched last time?
     */
    bool switched;
    /**
     * current simulation time
     */
    double time;

  };
}

#endif  /* REPLICA_H */
