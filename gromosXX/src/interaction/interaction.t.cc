/**
 * @file interaction.t.cc
 * test routines for the interaction library
 */

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "../debug.h"

#include "../math/gmath.h"
#include "../simulation/simulation.h"

#include "interaction.h"

#ifndef NDEBUG
int debug_level = 0;
#endif

using namespace math;

/**
 * provide comparision operators for the blitz TinyVector.
 * they should be implemented by blitz, but i cannot get
 * them to work?!
 */
bool operator==(math::Vec &t1, math::Vec &t2)
{
  bool b = true;
  for(int i=0; i<3; ++i)
    if (t1(i) != t2(i)) b = false;
  return b;
}
/**
 * != operator
 */
bool operator!=(math::Vec &t1, math::Vec &t2)
{
  return !(t1 == t2);
}

/**
 * Testing the forcefield
 */
int forcefield_test()
{
  try{
    
    int result = 0;
  
    const int SIZE = 10;
  
    simulation::system the_system;
    simulation::Topology the_topology;
  
    the_system.resize(SIZE);
    the_topology.resize(SIZE);

    // have to add real atoms...
    std::set<int> ex;
    for(int i=0; i<SIZE; ++i)
      the_topology.add_solute_atom("HI", 0, 0, 1.0, 0.0, false, ex, ex);
  

    the_system.pos() = tensor::i;
    the_system.vel() = 0.0;
    the_system.force() = 0.0;
  
    typedef simulation::Simulation<simulation::Topology, simulation::system>
      simulation_type;  

    simulation_type the_simulation(the_topology, the_system);

    // set a boundary condition
    the_system.boundary_condition(math::vacuum);

    // increase the cutoffs
    the_simulation.nonbonded().cutoff_short(2.0);
    the_simulation.nonbonded().cutoff_long(8.0);

    interaction::twin_range_pairlist<simulation_type> the_pairlist;
    the_pairlist.update(the_simulation);
  
    cout << the_pairlist.short_range();
    cout << the_pairlist.long_range();
  
    // let's create a forcefield...
    interaction::forcefield<simulation_type> the_forcefield;

    std::cout << "forcefield created" << std::endl;
    
    interaction::Nonbonded_Interaction<simulation_type,
      interaction::twin_range_pairlist<simulation_type> >
      *the_nb_interaction =
      new interaction::Nonbonded_Interaction<simulation_type,
      interaction::twin_range_pairlist<simulation_type> >;
  
    // add a lennard jones parameter
    interaction::lj_parameter_struct s;
    s.c6   = 2.261954e-03;
    s.c12  = 7.414932e-07;
    s.cs6  = 2.261954e-03;
    s.cs12 = 7.414932e-07;

    the_nb_interaction->resize(1);
    the_nb_interaction->add_lj_parameter(0, 0, s);

    std::cout << "nonbonded" << std::endl;
    
    the_forcefield.add_interaction(the_nb_interaction);

    std::cout << "added" << std::endl;
    
    // and calculate the forces...
    the_forcefield.calculate_interactions(the_simulation);

    std::cout << "calculated" << std::endl;
    
    // total force should be 0
    Vec t(math::epsilon, math::epsilon, math::epsilon);
    Vec v = sum(the_simulation.system().force());

    std::cout << "sum of the forces: " << v << std::endl;
    
    bool b = true;
    for(int i=0; i<3; ++i)
      if (v(i) < t(i)) b=true;
      else b = false;

    if (!b) ++result;
  
    return result;

  }
  catch(std::runtime_error e){
    std::cout << "exception in forcefield_test\n"
	      << e.what() << std::endl;
    throw;
  }

}

int main(int argc, char*argv[])
{
  if (argc >= 2){
    ::debug_level = atoi(argv[1]);
  }
  if (argc >= 3){
    interaction::debug_level = atoi(argv[2]);
  }
  if (argc >= 4){
    interaction::interaction_debug_level = atoi(argv[3]);
  }
  
  int r1;
  if ((r1 = forcefield_test()))
    std::cout << "forcefield_test failed" << std::endl;
  
  return r1;
}

