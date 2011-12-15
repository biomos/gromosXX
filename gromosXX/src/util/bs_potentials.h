/**
 * @file bs_umbrella.h
 * Implements the B&S-LEUS algorithm of Halvor Hansen et al.
 * Hansen et al., JCTC (2010), 6, 2622-2646. DOI: 10.1021/ct1003065
 */

#ifndef BS_POTENTIALS_H
#define	BS_POTENTIALS_H
#include "bs_coordinate.h"

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
    BS_Potential(int id, int num_gp, double force_const, double forceConstIncr):
            id(id), num_gp(num_gp), memForce(forceConstIncr), 
            force_const(force_const) {
      m_potential = 0;
      m_weight = 0;
      m_BoltzmannFactor = 0;
      memory.assign(num_gp, 0);
      auxillaryMemory.assign(num_gp, 0);
      half_force_const = 0.5 * force_const;
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
    void getMemory(std::vector<double> &newMemory) {newMemory = memory;}
    /**
     * set the memory to zero
     */
    void setMemoryToZero();
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
    void getDerivatives(BS_Vector &target) {target = potDerivatives;}
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
    } potentialType;
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
    int activeGridPoint;
    /**
     * used to calculate coordinate in "grid point coordinates"
     * Gamma_k / r_k (or u_l)
     */
    double gridPointScaling;
    /**
     * The memory
     */
    std::vector<double> memory;
    /**
     * The auxillary memory
     */
    std::vector<double> auxillaryMemory;
    /**
     * The force constant increment of the memory. k_LE
     */
    double memForce;
    /**
     * The force constant reduction factor. f_LE
     */
    //double memForceReduction;
    /**
     * the reduction counter
     */
    //unsigned int reductionCounter;
    /**
     * The force constant of the potential c_LE
     */
    double force_const;
    double half_force_const;
    /**
     * The derivative of the potential with respect to the internal coordinates
     * dB_k / dQ
     */
    BS_Vector potDerivatives;
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
    BS_Sphere(int id, int num_gp, double force_const, double forceConstIncr,
              BS_Vector center, double radius) : 
                BS_Potential(id, num_gp, force_const, forceConstIncr),
                center(center), radius(radius)
                {
                  gridPointScaling = ((num_gp - 1) / radius);
                  potentialType = bs_sphere;
                }
    virtual ~BS_Sphere(){}
    virtual bool calcPotential(BS_Vector &bs_pos);
    virtual std::string str();
    BS_Vector getCenter() {return center;}
  private:
    /**
     The center of the sphere
     */
    BS_Vector center;
    /**
     * The radius of the sphere (in reduced untis)
     */
    double radius;
  };
  
  /** 
   * @class BS_Stick
   * @ingroup util
   * Contains all the information about the sticks
   */
  class BS_Stick : public util::BS_Potential {
  public:
    BS_Stick(int id, int num_gp, double force_const, double forceConstIncr, 
             BS_Vector startPoint, BS_Vector endPoint, 
             double half_width) :
              BS_Potential(id, num_gp, force_const, forceConstIncr), 
              startPoint(startPoint), endPoint(endPoint), half_width(half_width)
              {
                endPoint.minus(startPoint, unitLongitudinal);
                length = unitLongitudinal.normalize();
                gridPointScaling = (num_gp - 1) / length;
                potentialType = bs_stick;
              }
    virtual ~BS_Stick(){}
    virtual bool calcPotential(BS_Vector &bs_pos);
    virtual std::string str();
  private:
    /**
     * The start point  of the stick in reduced coordinates
     */
    BS_Vector startPoint;
    /**
     * The end point  of the stick in reduced coordinates
     */
    BS_Vector endPoint;
    /**
     * the length between the start and the end point
     */
    double length;
    /**
     * The unit vector pointing from the start to the end
     */
    BS_Vector unitLongitudinal;
    /**
     * The half width of the line
     */
    double half_width;
  };

}


#endif	/* BS_POTENTIALS_H */

