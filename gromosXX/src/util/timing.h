/**
 * @file timing.h
 * timing routines.
 */

/**
 * @struct timing_struct
 * timings.
 */

struct timing_struct
{
  double total;

  double startup;

  double integration;
  int count_integration;

  double constraints;
  int count_constraints;

  double force;
  int count_force;

  double pairlist;
  int count_pairlist;

  double shortrange;
  int count_shortrange;

  double output;

  timing_struct()
    : total(0),
      startup(0),
      integration(0),
      count_integration(0),
      constraints(0),
      count_constraints(0),
      force(0),
      count_force(0),
      pairlist(0),
      count_pairlist(0),
      shortrange(0),
      count_shortrange(0),
      output(0)
  {
  }
};

/**
 * return the time now.
 */
double now();
void print_timing(std::ostream &os);

extern timing_struct timing;

  
