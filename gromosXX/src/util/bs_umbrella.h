/** 
 * @file   bs_umbrella.h
 * Contains all the B&S-LEUS potentials and unites them.
 *
 */

#ifndef BS_UMBRELLA_H
#define	BS_UMBRELLA_H
#include "bs_subspace.h"

namespace configuration {
  class Configuration;
}
namespace util {
  class BS_Potential;
  class BS_Coordinate;
  
  class BS_Umbrella {
  public:
    BS_Umbrella(){
      m_smoothing = 1.0;
      m_bsleus_total = 0;
    }
    ~BS_Umbrella();
    /**
     * Add a subspace to the umbrella
     * @param subspace
     */
    void addSubspace(BS_Subspace *subspace);
    /**
     * Apply the potentials and the forces to the configuration.
     * @param conf
     */
    void apply(configuration::Configuration &conf,
            simulation::Simulation &sim);
    /**
     * Set the memory of Potential
     * @param id    The id of the potential
     * @param type  The Type (Sphere / Stick)
     * @param memory    The memory
     */
    void setMemory(int id, BS_Potential::potential_enum type, 
            std::vector<double> &memory);
    /**
     * Set all memories to zero
     */
    void setMemoryToZero();
    /**
     * Get the Memory of a sphere or stick
     * @param[in] id
     * @param[in] type
     * @param[out] memory
     * @return wheter memory with specified id was found or not
     */
    bool getMemory(int id, BS_Potential::potential_enum type, 
            std::vector<double> &memory) const;
    /**
     * Get Number of Spheres and Sticks
     * @param[out] numSpheres
     * @param[out] numSticks
     */
    void getNumPotentials(int &numSpheres, int &numSticks) const;
    /**
     * Return the total B&S-LEUS Potential
     * @return 
     */
    double getTotalPotential() const;
    /**
     * Get the force in reduced coordinate space
     * @param[inout] force
     */
    void getForce(std::vector<double> &force) const;
    /**
     * Print out the informations about the potentials for a trajectory
     * @return 
     */
    std::string traj_str() const;
  private:
    /**
     * The smoothing factor for calculating the weight
     */
    double m_smoothing;
    /**
     * The active subspaces of the B&S-LEUS algorithm
     */
    std::vector<BS_Subspace *> m_subspaces;
    /**
     * Total bsleus energy
     */
    double m_bsleus_total;
  };
}

#endif	/* BS_UMBRELLA_H */

