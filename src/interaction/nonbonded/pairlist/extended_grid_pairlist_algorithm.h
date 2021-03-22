/**
 * @file grid_pairlist_algorithm.h
 * extended grid pairlist algorithm (beautiful implementation)
 */

#ifndef INCLUDED_EXTENDED_GRID_PAIRLIST_ALGORITHM_H
#define INCLUDED_EXTENDED_GRID_PAIRLIST_ALGORITHM_H

namespace math
{
  template<math::boundary_enum>
  class Periodicity;
}

namespace interaction
{
  class Pairlist;
  class Nonbonded_Parameter;
  
  template<typename t_interaction_spec>
  class Nonbonded_Innerloop;
  template<typename t_interaction_spec, typename t_perturbation_details>
  class Perturbed_Nonbonded_Innerloop;
  
  /**
   * @class Extended_Grid_Pairlist_Algorithm
   * create an atomic pairlist with a
   * chargegroup based or atom based
   * cut-off criterion.
   * Frist implementation only rectangular,
   * non-perturbed, chargegroup based
   */
  class Extended_Grid_Pairlist_Algorithm :
    public Failing_Pairlist_Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Extended_Grid_Pairlist_Algorithm();

    /**
     * Destructor.
     */
    virtual ~Extended_Grid_Pairlist_Algorithm();

    /**
     * init grid
     * get some information
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

    /**
     * prepare the pairlists
     */    
    virtual int prepare
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation &sim
     );
    
    /**
     * update the pairlist
     */
    virtual void update
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::PairlistContainer & pairlist,
     unsigned int begin,
     unsigned int end,
     unsigned int stride
     );

    /**
     * update the perturbed pairlist
     */
    virtual void update_perturbed
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::PairlistContainer & pairlist,
     interaction::PairlistContainer & perturbed_pairlist,
     unsigned int begin,
     unsigned int end, 
     unsigned int stride
     );
    
  protected:

    void set_cutoff(double const cutoff_short, double const cutoff_long)
    {
      m_cutoff_long = cutoff_long;
      m_cutoff_short = cutoff_short;
      m_cutoff_short_2 = cutoff_short * cutoff_short;
      m_cutoff_long_2  = cutoff_long * cutoff_long;
    }

    /**
     * timing information.
     */
    virtual void print_timing(std::ostream & os)
    {
      /*
      os << "            "
	 << std::setw(32) << std::left << "solv - solv"
	 << std::setw(20) << m_solvent_solvent_timing << "\n"
	 << "            "
	 << std::setw(32) << std::left << "spc - spc"
	 << std::setw(20) << m_spc_timing << "\n";
      */
    }
      
  private:
    /**
     * @class Grid
     * grid information
     */
    class Grid
    {
    public:
      Grid()
	: shift_space(0)
      {
      }
      
      struct Particle
      {
	Particle() : i(0), x(0.0), y(0.0), z(0.0), shift_index(0) 
	{
	}
	Particle(int i, double x, double y, double z, int shift_index=0) 
	  : i(i), x(x), y(y), z(z), shift_index(shift_index)
	{
	}
	Particle(int i, math::Vec const & v, int shift_index=0) 
	  : i(i), x(v(0)), y(v(1)), z(v(2)), shift_index(shift_index)
	{
	}
	Particle(Particle const & p) 
	  : i(p.i), x(p.x), y(p.y), z(p.z), shift_index(p.shift_index)
	{
	}
	void shift(Particle const & p, int shift_ind, std::vector<math::Vec> const & shift_vec)
	{
	  this->shift_index = shift_ind;
	  x = p.x + shift_vec[shift_index](0);
	  y = p.y + shift_vec[shift_index](1);
	  z = p.z + shift_vec[shift_index](2);
	  i = p.i;
	}
	
	int i;
	double x, y, z;
	int shift_index;
      };
      /**
       * cell dimensions
       */
      double a, b, c;
      /**
       * number of cells
       */
      int Na, Nb, Nc;
      /**
       * number of cells for extended (plane)
       */
      int Na_ex, Nb_ex, Nc_ex;
      /**
       * particles per cell (int()+1)
       */
      int Pcell;
      
      /**
       * additional space to put in particles
       */
      int shift_space;
      
      /**
       * particles in cells
       */
      std::vector<std::vector<Particle> > p_cell;
      /**
       * cell start in p_cell list
       */
      std::vector<std::vector<int> > cell_start;
      /**
       * count of particles per cell
       */
      std::vector<std::vector<int> > count;
      /**
       * cell index from particle (start) index
       */
      std::vector<std::vector<int> > cell_index;
      /**
       * planes inside cutoff
       */
      int mask_z;
      /**
       * the mask
       */
      std::vector<std::vector<int> > mask;
    };

    /**
     * calculate properties of the grid
     */
    void grid_properties
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     );

    /**
     * calculate the mask
     */
    void calculate_mask();
    
    /**
     * print the mask
     */
    void print_mask();

    /**
     * prepare the grid
     */
    int prepare_grid
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim
     );

    /**
     * collapse the grid
     * (remove all empty space)
     */
    void collapse_grid();

    /**
     * check the grid and print it...
     */
    void print_grid();
    
    /**
     * prepare a plane...
     */
    void prepare_plane
    (
     int z,
     std::vector<Grid::Particle> & p_plane, 
     std::vector<int> & cell_start
     );

    void print_plane
    (
     int z,
     std::vector<Grid::Particle> & p_plane, 
     std::vector<int> & cell_start
     );

    /**
     * update the pairlist
     */
    void _update
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::PairlistContainer & pairlist,
     unsigned int begin,
     unsigned int end,
     unsigned int stride
     );
     
    void _update_perturbed
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     interaction::PairlistContainer & pairlist,
     interaction::PairlistContainer & perturbed_pairlist,
     unsigned int begin,
     unsigned int end, 
     unsigned int stride
     );

    bool insert_pair
    (
     topology::Topology & topo,
     interaction::Pairlist & pairlist,
     interaction::Pairlist & perturbed_pairlist,
     int a1, int a2,
     bool scaled_only
    );
    
    bool excluded_solute_pair
    (
     topology::Topology & topo,
     unsigned int i,
     unsigned int j
     );

    /** 
     * the grid
     */
    Grid m_grid;

    /**
     * squared shortrange cutoff.
     */
    double m_cutoff_short_2;
    /**
     * squared longrange cutoff.
     */
    double m_cutoff_long_2;
    /**
     * longrange cutoff.
     */
    double m_cutoff_long;
    /**
     * shortrange cutoff.
     */
    double m_cutoff_short;

    /**
     * periodic shifts
     */
    std::vector<math::Vec> m_shift_vector;

    /**
     * periodic shifts for interchanged atoms (second perturbed)
     */
    std::vector<math::Vec> m_reverse_shift_vector;
  };
} // interaction

#endif
