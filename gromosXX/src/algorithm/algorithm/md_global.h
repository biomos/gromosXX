/**
 * @file md_global.h
 * global functions to get an md simulation started.
 */

#ifndef INCLUDED_MD_GLOBAL_H
#define INCLUDED_MD_GLOBAL_H

namespace algorithm
{
  /**
   * perform an MD simulation.
   */
  int do_md(io::Argument &args);
  
} // algorithm

#endif
