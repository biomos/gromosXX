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
     * @param[in] subspaces Add the subspaces to the umbrella
     */
    void addSubspace(std::vector<BS_Subspace  *> &subspaces);
    /**
     * Apply the potentials and the forces to the configuration.
     * @param conf
     */
    void apply(configuration::Configuration &conf,
            simulation::Simulation &sim);
    /**
     * Set the memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[in] memory    The memory
     */
    void setMemory(int id, BS_Potential::potential_enum type, 
            std::vector<double> &memory);
    /**
     * Set the auxiliary memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[in] memory    The auxiliary memory
     * @param[in] auxCounter    The auxiliary memory Counter
     * @param[in] redCounter    The reduction Counter
     */
    void setAuxMemory(int id, BS_Potential::potential_enum type, 
            std::vector<double> &memory, int auxCounter, int redCounter);
    /**
     * Set all memories to zero
     */
    void setMemoryToZero();
    /**
     * Set all auxiliary memories to zero
     */
    void setAuxMemoryToZero();
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
     * Get the Auxiliary Memory of a sphere or stick
     * @param[in] id
     * @param[in] type
     * @param[out] memory
     * @param[out] auxCounter    The auxiliary memory Counter
     * @param[out] redCounter    The reduction Counter
     * @return whether memory with specified id was found or not
     */
    bool getAuxMemory(int id, BS_Potential::potential_enum type, 
            std::vector<double> &memory, int &auxCounter, int &redCounter) const;
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
    /**
     * The force constant reduction factor: f_LE
     */
    double m_reductionFactor;
    /**
     * The local visiting cutoff: gamma_LE
     */
    int m_localCutoff;
    /**
     * The global visting cutoff: n_LE
     */
    int m_globalCutoff;
  };
}

#endif	/* BS_UMBRELLA_H */

