/**
 * @file average.h
 * storage of the various averages
 */

#ifndef INCLUDED_AVERAGE_H
#define INCLUDED_AVERAGE_H

namespace topology{
	class Topology;
}
namespace configuration{
	class Configuration;
}
namespace simulation
{
  class Parameter;
  class Simulation;
}

namespace configuration
{
  
  /**
   * @class Average
   * storage of the various average properties.
   * also does block averaging.
   */
  class Average : public algorithm::Algorithm
  {
  public:
    /**
     * @class Block_Average
     * contains block averages.
     */
    class Block_Average;

    /**
     * Constructor.
     */
    Average();

    /**
     * set to zero.
     */
    void zero();

    /**
     * update the averages.
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * resize the arrays.
     */
    void resize(topology::Topology const & topo,
		configuration::Configuration const & conf,
		simulation::Parameter const &param);
    
    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false)
    {
      os << "Averages\n";
      return 0;
    };

    /**
     * simulation averages accessor.
     */
    Block_Average const & simulation()const 
    {
      return m_sim_avg;
    }
    
    /**
     * const block averages accessor.
     */
    Block_Average const & block()const
    {
      return m_block_avg;
    }

    /**
     * block averages accessor.
     */
    Block_Average & block()
    {
      return m_block_avg;
    }
    
    class Block_Average 
    {
    public:
      /**
       * Constructor.
       */
      Block_Average() { zero(); }
      
      /**
       * set to zero.
       */
      void zero();

      /**
       * resize the arrays.
       */
      void resize(topology::Topology const & topo,
		  configuration::Configuration const & conf,
		  simulation::Parameter const &param);

      /**
       * update the averages.
       */
      void update(Block_Average const & old,
		  topology::Topology &topo, 
		  configuration::Configuration &conf,
		  simulation::Simulation &sim);

      /**
       * get the average energy and fluctuations
       * @param energy average energy
       * @param fluctuation energy fluctuations
       */
      void energy_average(Energy &energy, Energy &fluctuation)const;
      
      /**
       * get the average energy lambda derivative and its fluctuations
       * @param energy average energy lambda derivative
       * @param fluctuation energy lambda derivative fluctuations
       * @param lambda the current lambda value
       * @param lambda_fluct lambda value for the fluctuations
       * @param dlamt slow growth change of lambda per step
       * if dlamt != 0 the integrated energy lambda derivative is returned.
       */
      void energy_derivative_average(Energy &energy, Energy &fluctuation,
				     double & lambda, double & lambda_fluct,
				     double const dlamt)const;
      
      /**
       * get the average energy and fluctuations
       * @param pressure average
       * @param perssure_fluctuations pressure fluctuations
       * @param virial average
       * @param virial_fluctuations virial fluctuations
       * @param kinetic_energy tensor average
       * @param kinetic_energy_fluctuations average
       */
      void pressure_average(math::Matrix &pressure, 
			    math::Matrix &pressure_fluctuations,
			    math::Matrix &virial,
			    math::Matrix &virial_fluctuations,
			    math::Matrix &kinetic_energy,
			    math::Matrix &kinetic_energy_fluctuations)const;
      
      /**
       * get the averages and fluctuations of various other properties.
       * @param mass average
       * @param mass_fluctuations mass fluctuations
       * @param volume volume average
       * @param volume_fluctuations volume fluctuations
       * @param box average
       * @param box_fluctuations box fluctuations
       * @param scaling average of the scaling factors for temperature scaling
       * @param scaling_fluctuations fluctuations of the scaling factors
       */
      void mbs_average(double & mass, double & mass_fluctuations,
		       double & volume, double & volume_fluctuations,
		       math::Box & box, math::Box & box_fluctuations,
		       std::vector<double> & scaling,
		       std::vector<double> & scaling_fluctuations)const;

      /**
       * get the average surface areas and fluctuation
       */
      void sasa_average(std::vector<double> & sasa,
                        std::vector<double> & sasa_fluctuations,
                        double & sasatot, double & sasatot_fluct)const;
      /**
       * get the average volumes and fluctuation
       */
      void sasavol_average(std::vector<double> & sasavol,
                           std::vector<double> & sasavol_fluctuations,
                           double & sasavol_tot, double & sasavol_totfluct)const;

    private:
      /**
       * the time.
       */
      double time;
      /**
       * the average energies.
       */
      Energy energy_avg;
      /**
       * the squared averages.
       */
      Energy energy_fluct;
      /**
       * the average energy lambda derivatives.
       */
      Energy energy_derivative_avg;
      /**
       * the squared averages of the energy lambda derivatives.
       */
      Energy energy_derivative_fluct;
      /**
       * the average pressure.
       */
      math::Matrix pressure_avg;
      /**
       * the squared average pressure.
       */
      math::Matrix pressure_fluct;
      /**
       * average virial.
       */
      math::Matrix virial_avg;
      /**
       * virial fluctuations.
       */
      math::Matrix virial_fluct;
      /**
       * average kinetic energy tensor.
       */
      math::Matrix ekin_tensor_avg;
      /**
       * kinetic energy tensor fluctuations.
       */
      math::Matrix ekin_tensor_fluct;
      /**
       * average mass.
       */
      double mass_avg;
      /**
       * mass fluctuations.
       */
      double mass_fluct;
      /**
       * average temperature scaling factors.
       */
      std::vector<double> temp_scaling_avg;
      /**
       * temperature scaling factors fluctuation.
       */
      std::vector<double> temp_scaling_fluct;
      /**
       * volume average.
       */
      double volume_avg;
      /**
       * volume fluctuations.
       */
      double volume_fluct;
      /**
       * box average.
       */
      math::Box box_avg;
      /**
       * box fluctuations.
       */
      math::Box box_fluct;
      /**
       * lambda average
       */
      double lambda_avg;
      /**
       * lambda fluctutaitons.
       */
      double lambda_fluct;
      /**
       * the average total sasa.
       */
      double sasatot_avg;
      /**
       * the total sasa fluctuation.
       */
      double sasatot_fluct;
      /**
       * the average total volume.
       */
      double sasa_buriedvol_tot_avg;
      /**
       * the total volume fuctuation.
       */
      double sasa_buriedvol_tot_fluct;
      /**
       * the average sasa.
       */
      std::vector<double> sasa_avg;
      /**
       * the sasa fluctuations.
       */
      std::vector<double> sasa_fluct;
      /**
       * the average volume.
       */
      std::vector<double> sasa_buriedvol_avg;
      /**
       * the volume fuctuations.
       */
      std::vector<double> sasa_buriedvol_fluct;

      ////////////////////////////////////////////////////

      /**
       * update energies.
       */
      void update_energies(Energy & avg, Energy & fluct,
			   Energy const & e,
			   Energy const & old_avg,
			   Energy const & old_fluct,
			   double dt, double dlamt = 0.0);
      
      /**
       * update the pressure.
       */
      void update_volumepressure(topology::Topology const & topo,
				 configuration::Configuration const & conf,
				 simulation::Simulation const & sim,
				 Block_Average const & old);

      /**
       * update the sasa.
       */
      void update_sasa(topology::Topology const & topo,
                       configuration::Configuration const & conf,
                       simulation::Simulation const & sim,
                       Block_Average const & old);

      /**
       * update the volume.
       */
      void update_sasavol(topology::Topology const & topo,
                          configuration::Configuration const & conf,
                          simulation::Simulation const & sim,
                          Block_Average const & old);
      
      /**
       * calculate average energies and fluctuations.
       */
      void energy_average_helper(configuration::Energy const & avg,
				 configuration::Energy const & fluct,
				 configuration::Energy &energy, 
				 configuration::Energy &fluctuation,
				 double const dlamt)const;
      
      // friend class Average;

    };
    
  private:

    /**
     * averages over the simulation.
     */
    Block_Average m_sim_avg;
    /**
     * block averages.
     */
    Block_Average m_block_avg;

    ////////////////////////////////////////////////////


  }; // Average
  
} // configuration

#endif
