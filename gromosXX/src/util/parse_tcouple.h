/**
 * @file parse_tcouple.h
 * parse tcouple into multibath
 */

#ifndef INCLUDED_PARSE_TCOUPLE_H
#define INCLUDED_PARSE_TCOUPLE_H

namespace util
{
  /**
   * parse TCOUPLE block into MULTIBATH block.
   */
  void parse_TCOUPLE(simulation::Parameter &param,
		     topology::Topology const & topo);
  
}

#endif
