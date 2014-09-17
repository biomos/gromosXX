/**
 * @file bs_potentials.h
 * Implements the B&S-LEUS algorithm of Halvor Hansen et al.
 * Hansen et al., JCTC (2010), 6, 2622-2646. DOI: 10.1021/ct1003065
 */

#ifndef BS_POTENTIALS_H
#define	BS_POTENTIALS_H
#include "bs_vector.h"

namespace util{
  
  /**
   @class BS_Potential
   @ingroup util
   The base class for the ball and stick potentials
   */
  class BS_Potential{
  public:
    /**
     Constructor
     */
    BS_Potential(int id, unsigned int num_gp, double force_const):
            id(id), num_gp(num_gp), m_force_const(force_const) {
      m_potential = 0;
      m_weight = 0;
      m_BoltzmannFactor = 0;
      m_activeGridPoint = -1;
      m_memory.assign(num_gp, 0);
      m_auxiliaryMemory.assign(num_gp, 0);
      m_half_force_const = 0.5 * m_force_const;
    }
    /**
     * Destructor
     */
    virtual ~BS_Potential(){};
            
    /* Members */
    /**
     * The id of the potential
     */
    int id;
    /**
     * The number of grid points
     */
    unsigned int num_gp;
    /**
     * An intelligent version of eq. (1) & (2) for calculating the potential
     * inside the sphere / stick
     * @param gridPointCoordinate[in] = (num_gp - 1) / Characteristic_Length * distance
     * @param potential[out]
     * @param force[out]
     * @return the potential
     */
    void calcPotAndForce(double gridPointCoordinate, 
                         double &potential, 
                         double &force);
    /**
     * Update the memory
     * @param updateAuxMem  Should we update the Auxiliary Memory aswell?
     * @return  Are all auxillary memory points bigger than the local cutoff?
     */
    bool updateMemory(bool updateAuxMem);
    /**
     * Set the memory to newMemory
     * @param newMemory the new Memory
     */
    void setMemory(std::vector<double> &newMemory); // {memory = newMemory;}
    /**
     * Get the memory of the potential
     * @param[out] newMemory
     */
    void getMemory(std::vector<double> &newMemory) {newMemory = m_memory;}
    /**
     * set the memory to zero
     */
    void setMemoryToZero();
    /**
     * Set the auxiliary memory to newMemory
     * @param[in] newMemory the new Memory
     * @param[in] auxCounter    The auxiliary memory Counter
     * @param[in] redCounter    The reduction Counter
     */
    void setAuxMemory(std::vector<double> &newMemory); 
    /**
     * Get the auxiliary memory of the potential
     * @param[out] newMemory
     * @param[out] auxCounter    The auxiliary memory Counter
     * @param[out] redCounter    The reduction Counter
     */
    void getAuxMemory(std::vector<double> &newMemory);
    /**
     * set the auxiliary memory to zero
     */
    void setAuxMemoryToZero();
    /**
     * Set the paramter for the memory update scheme
     * @param[in] forceIncrement    basis force constant increment k_LE
     * @param[in] localCutoff       local visiting cutoff gamma_LE
     */
    void setMemoryParameters(double forceIncrement, int localCutoff);
    /**
     * Calculate the potential and the derivatives
     * @return The Potential
     */
    virtual double calcPotential(BS_Vector &bs_pos) = 0;
    /**
     * Return the potential (B_m)
     */
    double getPotential() {return m_potential;}
    /**
     * get the derivatives of the potential dB_m / dQ
     * @param[inout] target the derivatives will be saved in target
     */
    void getDerivatives(BS_Vector &target) {target = m_potDerivatives;}
    /**
     * Calculate the weight:
     * w_m = BoltzmannFactor / totalPartitionFct
     * w_m = exp[ -beta B_m] / sum_n( exp[-beta * (B_n - B_m)])
     * @param[in] {B_n}
     * @param[in] beta
     * @return the weight w_m
     */
    double calcWeight(std::vector<double> &potentials, double beta);
    /**
     * Set the wight of the potential
     * @param weight
     */
    void setWeight(double weight){m_weight = weight;}
    /**
     * Return the weight
     * @return  the weight
     */
    //double getWeight(){return m_weight;}
    int getActiveGridPoint(){return m_activeGridPoint;}
    /**
     * The Type of the Potential (Sphere / Stick)
     */
    enum potential_enum {
      bs_sphere,
      bs_stick,
      bs_snake,
      bs_pipe
    } m_potentialType;
    // Return a string for the type. Maybe not the nicest solution....
    static const char* potentialType(int i){
      const char *type[] = {"Sphere", "Stick", "Snake", "Pipe"};
      return type[i];
    }
    /**
     * Returns a string with informations about the potential
     * @return 
     */
    virtual std::string str(){return "Base Potential";}
    /**
     * Returns a string for writing the special trajectory 
     */
    std::string traj_str();
  protected:
    /**
     * The weight of a potential. Eq. (17)
     */
    double m_weight;
    /**
     * The current active grid point
     */
    int m_activeGridPoint;
    /**
     * used to calculate coordinate in "grid point coordinates"
     * Gamma_k / r_k (or u_l)
     */
    double m_gridPointScaling;
    /**
     * The memory
     */
    std::vector<double> m_memory;
    /**
     * The force constant increment of the memory. k_LE
     */
    double m_forceIncrement;
    /**
     * The auxiliary memory: A(t)
     */
    std::vector<double> m_auxiliaryMemory;
    /**
     * local cutoff: gamma_LE
     */
    double m_localCutoff;
    /**
     * The force constant of the potential c_LE
     */
    double m_force_const;
    double m_half_force_const;
    /**
     * The derivative of the potential with respect to the internal coordinates
     * dB_k / dQ
     */
    BS_Vector m_potDerivatives;
    /**
     * The value of the potential at the current place. Should be in kJ / mol.
     */
    double m_potential;
    /**
     * The Boltzmann factor exp [ - beta B_m)]
     */
    double m_BoltzmannFactor;
  };
  
  /**
   * @class BS_Sphere
   * @ingroup util
   * Contains all the informations about a B&S-LEUS sphere
   */
  class BS_Sphere : public BS_Potential {
  public:
    BS_Sphere(int id, int num_gp, double force_const,
              BS_Vector center, double radius) : 
                BS_Potential(id, num_gp, force_const),
                m_center(center), m_radius(radius)
                {
                  m_gridPointScaling = ((num_gp - 1) / radius);
                  m_potentialType = bs_sphere;
                }
    virtual ~BS_Sphere(){}
    virtual double calcPotential(BS_Vector &bs_pos);
    virtual std::string str();
    BS_Vector getCenter() {return m_center;}
  private:
    /**
     The center of the sphere
     */
    BS_Vector m_center;
    /**
     * The radius of the sphere (in reduced untis)
     */
    double m_radius;
  };
  
  /** 
   * @class BS_Stick
   * @ingroup util
   * Contains all the information about the sticks
   */
  class BS_Stick : public util::BS_Potential {
  public:
    BS_Stick(int id, int num_gp, double force_const, 
             BS_Vector startPoint, BS_Vector endPoint, 
             double half_width) :
              BS_Potential(id, num_gp, force_const), 
              m_startPoint(startPoint), m_endPoint(endPoint), 
              m_half_width(half_width)
              {
                endPoint.minus(startPoint, m_unitLongitudinal);
                m_length = m_unitLongitudinal.normalize();
                m_gridPointScaling = (num_gp - 1) / m_length;
                m_potentialType = bs_stick;
              }
    virtual ~BS_Stick(){}
    virtual double calcPotential(BS_Vector &bs_pos);
    virtual std::string str();
    double distance() {return m_distance;}
  protected:
    /**
     * The start point  of the stick in reduced coordinates
     */
    BS_Vector m_startPoint;
    /**
     * The end point  of the stick in reduced coordinates
     */
    BS_Vector m_endPoint;
    /**
     * the length between the start and the end point
     */
    double m_length;
    /**
     * The unit vector pointing from the start to the end
     */
    BS_Vector m_unitLongitudinal;
    /**
     * The half width of the line
     */
    double m_half_width;
    double m_distance;
  };
  
  class BS_Distorted_Stick : public util::BS_Stick {
  public:
    BS_Distorted_Stick (int id, int num_gp, double force_const, 
             BS_Vector startPoint, BS_Vector endPoint, 
             double half_width, int is_where);
    BS_Vector direction() { return m_unitLongitudinal; }
    /**
     * Initilaize all vectors
     */
    void init_vectors(BS_Vector prev_dir, BS_Vector next_dir);
    /**
     * Is there a possibility for the Distorted stick to interact with the system?
     * Calculates already the distance.
     */
    bool interacts(BS_Vector &bs_pos);
    /**
     * Calculate the potential
     */
    virtual double calcPotential(BS_Vector &bs_pos);
    virtual std::string str();
    
    enum {
      in_between = 0,
      start,
      end
    } where;
  protected:
    //BS_Vector m_start_norm;
    //BS_Vector m_end_norm;
    
    // $\vec{r}$ or $\vec{\rho}$
    BS_Vector radialVec;
    
    // $\vec{\tilde{p}}$
    //BS_Vector m_unit_tilted_perp;
    //BS_Vector m_tilted_perp_deriv;
    //double m_tilted_perp_dist;
    BS_Vector m_unit_perp;
    double m_perp_dist;
    
    BS_Vector m_h_base;
    BS_Vector m_H_base;
    
    BS_Vector m_deriv_h_over_H;
    double m_h;
    double m_H;
    double m_h_over_H;
    bool m_at_end;
  };

  /**
   * 
   * @class BS_Snake
   * 
   * A snake like potential that behave a bit like connected sticks.
   */
  class BS_Snake : public util::BS_Potential {
  public:
    BS_Snake(int id, int num_gp, double force_const,
            std::vector<BS_Vector> points, double halfwidth);
    virtual double calcPotential(BS_Vector &bs_pos);
    virtual std::string str();
  private:
    //std::vector<BS_Vector> m_points;
    std::vector<util::BS_Distorted_Stick> m_sticks;
  };
  
  struct BS_Pipe_Param {
    BS_Pipe_Param(BS_Vector point, double inner_width, double outer_width) :
    point(point), inner_width(inner_width), outer_width(outer_width) {}
    BS_Pipe_Param() :
    point(), inner_width(-1), outer_width(-1) {}
    BS_Vector point;
    double inner_width;
    double outer_width;
  };
  
  class BS_Pipe : public util::BS_Potential {
  public:
    BS_Pipe(int id, int num_gp_long, int num_gp_perp, double force_const,
            util::BS_Pipe_Param start, util::BS_Pipe_Param end);
    virtual double calcPotential(BS_Vector &bs_pos);
    virtual std::string str();
  protected:
    BS_Pipe_Param m_start;
    BS_Pipe_Param m_end;
    
    BS_Vector m_unitLongitudinal;
    
    double m_long_conversion;
    double m_length;
    double m_inner_slope;
    double m_outer_slope;
    
    const int m_num_long_gp;
    const int m_num_perp_gp;
    int m_long_gp;
    int m_perp_gp;
    
    enum where_in_potential {
      inside = 0,
      below_lower,
      above_upper
    };
  };

}


#endif	/* BS_POTENTIALS_H */

