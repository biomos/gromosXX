/**
 * @file interaction.t.cc
 * test routines for the interaction library
 */

#include "../math/gmath.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "../simulation/simulation.h"

#include "interaction/interaction.h"
#include "pairlist/simple_pairlist.h"
#include "interaction/nonbonded_interaction.h"


#include "forcefield/forcefield.h"

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
  int result = 0;
  
  const int SIZE = 10;
  
  simulation::system the_system;
  simulation::topology the_topology;
  
  the_system.resize(SIZE);
  the_topology.resize(SIZE);
  the_topology.num_solute_atoms(SIZE);

  // initialize everything with zero
  Vec t(0.0, 0.0, 0.0), t2(1, 1, 1);

  the_system.pos() = tensor::i;
  the_system.vel() = 0.0;
  the_system.force() = 0.0;
  
  the_topology.mass() = 1.0;

  typedef simulation::simulation<simulation::topology, simulation::system>
    simulation_type;  

  simulation_type the_simulation(the_topology, the_system);

  interaction::simple_pairlist<simulation_type> the_pairlist;
  the_pairlist.make_pairlist(the_simulation);
  
  // the_pairlist.print_pairlist(cout);
  
  interaction::simple_pairlist<simulation_type>::iterator it =
    the_pairlist.begin();
  
  // let's create a forcefield...

  interaction::forcefield<simulation_type> the_forcefield;

  interaction::nonbonded_interaction<simulation_type> *the_nb_interaction =
    new interaction::nonbonded_interaction<simulation_type>;
  
  the_forcefield.add_interaction(the_nb_interaction);

  // and calculate the forces...
  the_forcefield.calculate_interactions(the_simulation);

  // total force should be 0
  Vec v = sum(the_simulation.system().force());

  bool b = true;
  if (v == t) b=true;
  else b = false;

  if (!b) ++result;

  return result;
}

int main()
{
  int r1;
  if ((r1 = forcefield_test()))
    std::cout << "forcefield_test failed" << std::endl;
  
  return r1;
}

