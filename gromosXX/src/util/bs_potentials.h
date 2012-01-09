/**
 * @file bs_umbrella.h
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
    BS_Potential(int id, int num_gp, double force_const):
            id(id), num_gp(num_gp), m_force_const(force_const) {
      m_potential = 0;
      m_weight = 0;
      m_BoltzmannFactor = 0;
      m_memory.assign(num_gp, 0);
      m_auxiliaryMemory.assign(num_gp, 0);
      m_half_force_const = 0.5 * m_force_const;
      m_auxiliaryCounter = 0;
      m_reductionCounter = 0;
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
    int num_gp;
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
     */
    void updateMemory();
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
    void setAuxMemory(std::vector<double> &newMemory, int auxCounter, int redCounter); 
    /**
     * Get the auxiliary memory of the potential
     * @param[out] newMemory
     * @param[out] auxCounter    The auxiliary memory Counter
     * @param[out] redCounter    The reduction Counter
     */
    void getAuxMemory(std::vector<double> &newMemory, int &auxCounter, int &redCounter);
    /**
     * set the auxiliary memory to zero
     */
    void setAuxMemoryToZero();
    /**
     * Set the paramter for the memory update scheme
     * @param[in] forceIncrement    basis force constant increment k_LE
     * @param[in] reductionFactor   reduction factor f_LE
     * @param[in] localCutoff       local visiting cutoff gamma_LE
     * @param[in] globalCutoff      global visiting cutoff n_LE
     */
    void setMemoryParameters(double forceIncrement, double reductionFactor,
                             int localCutoff, int globalCutoff);
    /**
     * Calculate the potential and the derivatives
     * @return are we inside of the boundaries and do therefore have to update
     *         the memory?
     */
    virtual bool calcPotential(BS_Vector &bs_pos) = 0;
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
     * Calculate the Boltzmann Factor and return  it
     * @param[in] beta the current value of beta = 1/k_B * T
     * @return The Boltzmann Factor
     */
    double calcBoltzmann(const double beta);
    /**
     * Calculate the weight:
     * w_m = BoltzmannFactor / totalPartitionFct
     * w_m = exp[ -beta B_m] / sum( exp[-beta * B_m])
     * @param[in] totalPartitionFct
     * @return the weight w_m
     */
    double calcWeight(double totalPartitionFct);
    /**
     * Assign the current position the a current position and report 
     * back, wheter sucessfull or not.
     * @param bs_pos
     * @return 
     */
    //virtual bool assignGridPoint(BS_Vector &bs_pos) = 0;
    /**
     * Set the wight of the potential
     * @param weight
     */
    void setWeight(double weight){m_weight = weight;}
    /**
     * Return the weight
     * @return  the weight
     */
    double getWeight(){return m_weight;}
    /**
     * The Type of the Potential (Sphere / Stick)
     */
    enum potential_enum {
      bs_sphere,
      bs_stick
    } m_potentialType;
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
     * The force constant reduction factor. f_LE
     */
    double m_reductionFactor;
    /**
     * The auxiliary counter: N_C
     */
    unsigned int m_auxiliaryCounter;
    /**
     * the reduction counter: I_R
     */
    unsigned int m_reductionCounter;
    /**
     * local cutoff: gamma_LE
     */
    double m_localCutoff;
    /**
     * global cutoff: n_LE
     */
    double m_globalCutoff;
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
    BS_Sphere(int id, int num_gp, double force_const, //double forceConstIncr,
              BS_Vector center, double radius) : 
                BS_Potential(id, num_gp, force_const),
                m_center(center), m_radius(radius)
                {
                  m_gridPointScaling = ((num_gp - 1) / radius);
                  m_potentialType = bs_sphere;
                }
    virtual ~BS_Sphere(){}
    virtual bool calcPotential(BS_Vector &bs_pos);
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
    BS_Stick(int id, int num_gp, double force_const, //double forceConstIncr, 
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
    virtual bool calcPotential(BS_Vector &bs_pos);
    virtual std::string str();
  private:
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
  };

}


#endif	/* BS_POTENTIALS_H */

