/**
 * @file read_input.h
 * input parameters
 */

#ifndef INCLUDED_READ_INPUT_H
#define INCLUDED_READ_INPUT_H

class Simulation;

/**
 * read in data and setup an md simulation.
 */
int read_input(Argument const &args,
	       Simulation & sim,
	       std::ostream & os = std::cout);

#endif
