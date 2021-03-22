/* 
 * File:   eds.h
 * Author: haniels
 *
 * Created on August 2, 2011, 1:41 PM
 */

#ifndef EDS_H
#define	EDS_H
namespace algorithm
{
  /**
   * @class EDS
   * implements EDS.
   */
  class EDS : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    EDS() : Algorithm("EDS"), conf2(NULL) {}
    
    void set_conf2(configuration::Configuration & conf) {
      conf2 = &conf;
    }
    

    /**
     * Destructor.
     */
    virtual ~EDS(){}
    
    /**
     * 
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    
    {
      if (!quiet)
	os << "\tEDS\nEND\n";
      return 0;
    };

   private:
     configuration::Configuration * conf2;
  
  };
   
} // algorithm

#endif



