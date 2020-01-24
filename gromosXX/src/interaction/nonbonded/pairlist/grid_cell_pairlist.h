/**
 * @file grid_cell_pairlist.h
 * Implementation of a fast grid-cell algorithm
 */

#ifndef _GRID_CELL_PARILIST_H
#define	_GRID_CELL_PARILIST_H

#include "pairlist_algorithm.h"

namespace interaction {

  /**@class Reg_Mask
   * A trait class used, so make_mask_pointer can be implemented as template,
   * once for regulare shapes and and once for irregular shapes
   */
  class Reg_Mask {};

  /**@class Irr_Mask
   * A trait class used, so make_mask_pointer can be implemented as template,
   * once for regulare shapes and and once for irregular shapes
   */
  class Irr_Mask {};

  /**
   * @class Grid_Cell_Pairlist
   * Implementation of a fast grid-cell algorithm from
   * Tim N. Heinz, Philippe H. Hueneberger: A fast pairlist-construction
   * algorithm for molecular simulations under periodic boundary conditions.
   * J Comput Chem 25: 1474-1486, 2004
   */
  class Grid_Cell_Pairlist : public Failing_Pairlist_Algorithm {
  public:
    /**
     * Constructor
     */
    Grid_Cell_Pairlist(const topology::Topology & topo,
            const simulation::Simulation &sim);
    /**
     * Default Destructor
     */
    virtual ~Grid_Cell_Pairlist();
    /**
     * Initialize the pairlist
     */
    virtual int init(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation &sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    /**
     * prepare the pairlist(s).
     */
    virtual int prepare(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation &sim);
    /**
     * update the pairlist
     */
    virtual void update(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation &sim,
            interaction::PairlistContainer &pairlist,
            unsigned int begin, unsigned int end,
            unsigned int stride);

    /**
     * update the pairlist, separating perturbed and nonperturbed interactions
     */
    virtual void update_perturbed(
            topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            interaction::PairlistContainer & pairlist,
            interaction::PairlistContainer & perturbed_pairlist,
            unsigned int begin, unsigned int end,
            unsigned int stride);

  private:
    /**
     * Calculate parameters like dim_m, length_i etc.
     */
    int calc_par();
    /**
     * Creates the mask for regular shapes (as described in the appendix A)
     */
    bool make_mask(unsigned int delta_m, Reg_Mask &t);
    /**
     * creates the mask for irregular shapes (as described in the appendix B)
     */
    bool make_mask(int delta_m, Irr_Mask &t);
    /**
     * Create the mask pointer
     */
    template<typename trait_type>
    int make_mask_pointer();
    /*int make_mask_pointer(const std::vector<bool> & mask,
            const unsigned int m,
            const unsigned int size);*/
    /**
     * make the cell array and the cell pointer
     */
    template<math::boundary_enum b>
    int make_cell(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);
    /**
     * CG based cutoff trait
     */
    class cg_cutoff {};
    /**
     * atomic based cutoff trait
     */
    class atomic_cutoff {};
    /**
     * no perturbation trait
     */
    class no_perturbation {};
    /**
     * no perturbation trait
     */
    class do_perturbation {};
    /**
     * make the pairlist
     */
    template<math::boundary_enum b, class cutoff_trait, class perturbation_trait>
    int _pairlist(interaction::PairlistContainer & pairlist,
            interaction::PairlistContainer & perturbed_pairlist,
            unsigned int offset, unsigned int stride,
            const cutoff_trait & cutoff,
            const perturbation_trait & perturbation);
    /**
     * put the  charge groups inside the cell into the pairlist
     */
    template<math::boundary_enum b, class cutoff_trait, class perturbation_trait>
    inline int pair(unsigned int n1, unsigned int n2,
            interaction::PairlistContainer & pairlist,
            interaction::PairlistContainer & perturbed_pairlist,
            const math::Periodicity<b> &periodicity,
            const cutoff_trait & cutoff,
            const perturbation_trait & perturbation);
    /**
     * pair two solvent molecules
     */
    inline int pair_solvent(const unsigned int first, const unsigned int second,
            interaction::Pairlist &pairlist, const cg_cutoff & cutoff);
    /**
     * pair two solvent atoms
     */
    inline int pair_solvent(const unsigned int first, const unsigned int second,
            interaction::Pairlist &pairlist, const atomic_cutoff & cutoff);
    /**
     * put the atoms of a chargegroupe into the pairlist
     */
    template<class perturbation_trait>
    inline int pair_solute(const unsigned int first, const unsigned int second,
            interaction::Pairlist &pairlist, interaction::Pairlist &perturbed_pairlist, const cg_cutoff & cutoff,
            const perturbation_trait & perturbation);
    /**
     * put the atoms into the pairlist
     */
    template<class perturbation_trait>
    inline int pair_solute(const unsigned int first, const unsigned int second,
            interaction::Pairlist &pairlist, interaction::Pairlist &perturbed_pairlist, const atomic_cutoff & cutoff,
            const perturbation_trait & perturbation);
    /**
     * put the atoms of a chargegroupe into the pairlist
     */
    template<class perturbation_trait>
    inline int pair_solute_solvent(const unsigned int first, const unsigned int second,
            interaction::Pairlist &pairlist, interaction::Pairlist &perturbed_pairlist, const cg_cutoff & cutoff,
            const perturbation_trait & perturbation);
    /**
     * put the atoms into the pairlist
     */
    template<class perturbation_trait>
    inline int pair_solute_solvent(const unsigned int first, const unsigned int second,
            interaction::Pairlist &pairlist, interaction::Pairlist &perturbed_pairlist, const atomic_cutoff & cutoff,
            const perturbation_trait & perturbation);
    /**
     * Check, if two atoms are exclude from the pairlist
     */
    inline bool excluded_solute_pair(
            topology::Topology & topo,
            unsigned int i, unsigned int j);

    /**
     * insert a pair to the pairlist
     */
    inline void insert_pair(
            interaction::Pairlist & pairlist,
            interaction::Pairlist & perturbed_pairlist,
            int first, int second, const no_perturbation & perturbation);
    /**
     * insert a pair to the pairlist (perturbation code)
     */
    inline void insert_pair(
            interaction::Pairlist & pairlist,
            interaction::Pairlist & perturbed_pairlist,
            int first, int second, const do_perturbation & perturbation);

    /**
     * put the chargegroups into the box
     */
    template<math::boundary_enum b>
    void _prepare_cog(configuration::Configuration & conf,
            topology::Topology & topo);
    /**
     * Makes the combination table of the codes for irregular shapes aswell
     * as the corresponding code table for the parameters
     * Table 1 & 2, p. 1484
     */
    inline void make_code_table();
    /**
     * returns the minimum image value of its first argument based on the
     * periodicity defined by its second argument.
     */
    inline int minIm(int n, int num);
    /**
     * Returns the sign
     */
    template <typename T> inline int signum(T arg){
      if (arg > 0)
        return 1;
      else if (arg < 0)
        return -1;
      else
        return 0;
    }
    /**
     * returns the maximum of n or 1
     */
    template<typename t>
    inline t max(t n);
    void put_into_brickwall(math::Vec &v);

    // values
    /**
     * pointer to the topology
     */
    topology::Topology * mytopo;
    /**
     * pointer to the configuration
     */
    configuration::Configuration * myconf;
    /**
     * pointer to the simulation
     */
    simulation::Simulation * mysim;
    /**
     * Number of boxes in x direction (in the paper: N_x)
     */
    int num_x;
    /**
     * Number of boxes in y direction (in the paper: N_y)
     */
    int num_y;
    /**
     * Number of boxes in z direction (in the paper: N_z)
     */
    int num_z;
    /**
     * Number of boxes in y direction lowered by one (needed for calculations
     * with irregular shape)
     */
    int num_y_minus_1;
    /**
     * Number of boxes in z direction lowered by one (needed for calculations
     * with irregular shape)
     */
    int num_z_minus_1;
    /**
     * Number of boxes in x direction divided by two (needed for calculations
     * with irregular shape)
     */
    double num_x_half;
    /**
     * Number of boxes in y direction divided by two (needed for calculations
     * with irregular shape)
     */
    double num_y_half;
    /**
     * Number of boxes in z direction divided by two (needed for calculations
     * with irregular shape)
     */
    double num_z_half;
    /**
     * N_x * N_y + P = num_x * num_y + padding
     */
    int nxnyp;
    /**
     * length of the box in x direction (in the paper: L_x)
     */
    double length_x;
    /**
     * length of the box in y direction (in the paper: L_y)
     */
    double length_y;
    /**
     * length of the box in z direction (in the paper: L_z)
     */
    double length_z;
    /**
     * longrange cuttoff
     */
    double cutoff_lr;
    /**
     * shortrange cutoff squared
     */
    double cutoff_sr2;
    /**
     * longrange cuttoff squared
     */
    double cutoff_lr2;
    /**
     * shortrange cutoff
     */
    double cutoff_sr;
    /**
     * Length of the mask (in the paper: M = N_x * N_y * N_z)
     */
    int dim_m;
    /**
     * The length per cell (in the paper: l_x^2 = (L_x / N_x)^2 )
     */
    double l_x2;
    /**
     * The length per cell (in the paper: l_y^2 = (L_x / N_x)^2 )
     */
    double l_y2;
    /**
     * The length per cell (in the paper: l_z^2 = (L_x / N_x)^2 )
     */
    double l_z2;

    // For the irregularly shaped boxes
    /**
     * Are there irregular shapes?
     */
    bool irregular_shape;
    /**
     * The padding (P)
     */
    int padding;
    double delta_l_xy;
    double delta_l_xz;
    double delta_l_yz;
    /**
     * The offset lambda_xy
     */
    double lambda_xy;
    /**
     * The offset lambda_xz
     */
    double lambda_xz;
    /**
     * The offset lambda_yz
     */
    double lambda_yz;
    /**
     * the cases
     */
    int cases[10];
    /**
     * the tables (Table 2)
     */
    int a_code[7];
    int b_code[7];
    int c_code[7];
    /**
     * Mask_Pointer
     */
    std::vector<unsigned int> mask_pointer;
    /**
     * The number of stripes (french fries)
     */
    int stripes;
    /**
     * element of the cell vector
     */
    struct cell_element{
      /**
       * pointer to the charge group
       */
      //topology::Chargegroup_Iterator cg_it;
      /**
       * index of the grid box
       */
      int m;
      /**
       * index of the position list
       */
      unsigned int i;

      /**
       * Compare two cell elements
       */
      //bool operator<(const cell_element & lhs, const cell_element & rhs) {
      bool operator<(const cell_element & rhs) const {
        return (this->m < rhs.m);
      }
    };
    /**
     * The cell array 
     */
    std::vector<cell_element> cell;
    /**
     * Cell Pointer: points to the first element of a cell
     */
    std::vector<unsigned int> cell_pointer;
    /**
     * the chargegroup center of geometries.
     */
    math::VArray m_cg_cog;
    /**
     * Does cuda look for the solvent-solvent pairlist?
     */
    bool cuda;

    /**
     * the saved box for vacuum
     */
    math::Box box_backup;
    /**
     * is it vacuum?
     */
    bool is_vacuum;
    /**
     * create a box for vacuum
     */
    void create_vacuum_box();
     /**
     * restore the original box
     */
    void restore_vacuum_box();

    /**
     * index of the first solvent atom or cg
     */
    unsigned int first_solvent;

    /**
     * the number of atoms per solvent molecule
     */
    unsigned int num_atoms_per_solvent;

  };
}


#endif	/* _GRID_CELL_PARILIST_H */

