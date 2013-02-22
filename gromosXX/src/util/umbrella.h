/**
 * @file umbrella.h
 * Umbrellas for LE-US
 */

#ifndef INCLUDED_UMBRELLA_H
#define	INCLUDED_UMBRELLA_H

namespace configuration {
  class Configuration;
}

namespace util {
  class LE_Coordinate;
  class Umbrella_Weight;
  class Umbrella_Weight_Factory;

  /**
   * @class Umbrella
   * @ingroup util
   * @short Umbrella potential class
   *
   * The class Umbrella is used to combine multiple collective variables
   * (see @ref util::LE_Coordinate) into one umbrella potential. It holds
   * the visited configurations and applies biasing potentials to this
   * configurations.
   */
  struct Umbrella {
    /**
     * constructor
     */
    Umbrella(int id, unsigned int dim,
            Umbrella_Weight_Factory * factory = NULL);
    /**
     * copy constructor
     * this will pass the factory to the new instance
     */
    Umbrella(const Umbrella & u);
    /**
     * assignment operator
     */
    Umbrella & operator=(const Umbrella & u);
    /**
     * destructor
     */
    ~Umbrella();
    /**
     * the ID of the umbrella potential
     */
    int id;
    /**
     * a vector of the LE coordinates attached to the umbrella
     * @sa util::LE_Coordinate
     */
    std::vector<LE_Coordinate*> coordinates;

    /**
     * dimension of the potential
     */
    unsigned int dim() const {
      return functional_form.size();
    };
    /**
     * force constant
     */
    double force_constant;

    /**
     * @enum functional_form_enum
     * functional form of the biasing potential
     */
    enum functional_form_enum {
      /** polynomial (truncated) */ ff_polynomial,
      /** Gaussian */ ff_gaussian
    };
    /**
     * functional form of the biasing potential per dimension
     */
    std::vector<functional_form_enum> functional_form;
    /**
     * @enum variable_type_enum
     * variable type
     */
    enum variable_type_enum {
      /** unkown */ vt_unkown = 0,
      /** dihedral */ vt_dihedral = 1,
      /** distance */ vt_distance = 2,
      /** distancefield */ vt_distancefield = 6
    };
    /**
     * type of the variables
     */
    std::vector<variable_type_enum> variable_type;
    /**
     * width of function in units of grid spacing
     */
    std::vector<double> width;
    /**
     * width of function in units of grid spacing (in transformed units)
     */
    std::vector<double> width_rel;
    /**
     * cutoff applied to the range of action of the local functions
     * in units of grid spacing
     */
    std::vector<double> cutoff;
    /**
     * cutoff applied to the range of action of the local functions
     * in units of grid spacing (in tranformed units)
     */
    std::vector<double> cutoff_rel;
    /**
     * number of grid points along dimensions
     */
    std::vector<unsigned int> num_grid_points;
    /**
     * minimum grid point used along each dimension
     */
    std::vector<double> grid_min;
    /**
     * minimum grid point used along each dimension (in transformed units)
     */
    std::vector<double> grid_min_rel;
    /**
     * maximum grid point used along each dimension
     */
    std::vector<double> grid_max;
    /**
     * maximum grid point used along each dimension (in transformed units
     */
    std::vector<double> grid_max_rel;
    /**
     * grid spacing (in transformed units)
     */
    std::vector<double> grid_spacing_rel;

    /**
     * @struct leus_conf
     * @short a visted LEUS configuration
     *
     * Struct leus_conf holds the integer codes of the visited LEUS configurations
     * in its pos field. The struct is less-than comparable for efficient tree
     * searching.
     */
    struct leus_conf {
      /**
       * constructor
       * @param[in] dim the dimension of the configuration
       */
      leus_conf(unsigned int dim) {
        pos.resize(dim);
      }
      /**
       * constructor from a position vector
       * @param[in] pos the position vector to copy from
       */
      leus_conf(const std::vector<int> & pos) :
      pos(pos) {
      }
      /**
       * the coded index (along dimension) of the visted grid point
       */
      std::vector<int> pos;
      /**
       * operator to sort and compare (less-than comparable)
       */
      bool operator<(const leus_conf & c) const;
    };
    /**
     * a map to hold the visited configurations and the number of times it
     * was visited. A @ref leus_conf is used as key type. This allows for
     * efficient searching of the visited grid points. The value is the weight of
     * a configuration. If the configuration element
     * is not found in this map, the configuration wasn't visited so far.
     */
    std::map<leus_conf, Umbrella_Weight *> configurations;
    /**
     * is the umbrella still building up?
     */
    bool building;
    /**
     * is the umbrella enabled (i.e. applied, energy and forces are calculated)
     */
    bool enabled;
    /**
     * the factory for the umbrella weights
     */
    Umbrella_Weight_Factory * umbrella_weight_factory;
    /**
     * a function to read the configuration after setting the umbrella weights
     */
    void read_configuration();
    /**
     * transform the units of the grid. This is carried out for user friendlyness
     * For example we calculate internally in radients but it is more convenient
     * to have a input/output format in degree.
     */
    void transform_units();
    /**
     * calculate the attached coordinates.
     * This will invoke a loop over the attached coordinates and calculate
     * their values.
     * @param[inout] the configuration that is used to calculate the energy
     *               and force and wherin the force is stored.
     */
    void calculate_coordinates(configuration::Configuration & conf);
    /**
     * built-up the umbrella
     * This will find the grid point and add it to the visited configurations
     * @param[inout] conf the configuration
     */
    void build(configuration::Configuration & conf);
    /**
     * apply the umbrella
     * @param[inout] conf the configuration
     */
    void apply(configuration::Configuration & conf);
    /**
     * convert to a nice formated string
     */
    std::string str() const;
    /**
     * the block from the configuration file
     */
    std::string configuration_block;
    /**
     * the position of the stream in this block
     */
    std::istringstream::pos_type configuration_block_pos;
  };

}

#endif	/* INCLUDED_UMBRELLA_H */

