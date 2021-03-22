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
  
  /**
   * @class BS_Umbrella
   * 
   * Holds the subspaces and is basically 'vase', in which the whole B&S-LEUS
   * scheme sits.
   * 
   * @sa BS_Subspace
   */
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
    void addSubspaces(std::vector<BS_Subspace  *> &subspaces);
    /**
     * Apply the potentials and the forces to the configuration.
     * @param conf
     */
    void apply(configuration::Configuration &conf,
            simulation::Simulation &sim);
    /**
     * Set the memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] subid The id of the subspace
     * @param[in] type  The Type (Sphere / Stick)
     * @param[in] memory    The memory
     */
    void setMemory(int id, int subid, std::vector<double> &memory);
    /**
     * Set the auxiliary memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] subid The id of the subspace
     * @param[in] type  The Type (Sphere / Stick)
     * @param[in] memory    The auxiliary memory
     */
    void setAuxMemory(int id, int subid, std::vector<double> &memory);
    /**
     * Set the counters for a subspace
     * @param subid         The id of the subspace
     * @param auxCounter    The auxiliary Counter
     * @param redCounter    The reduction Counter
     */
    void setCounter(int subid, int auxCounter, int redCounter);
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
    bool getMemory(int id, unsigned int &subid, std::vector<double> &memory) const;
    /**
     * Get the Auxiliary Memory of a sphere or stick
     * @param[in] id
     * @param[in] type
     * @param[out] memory
     * @return whether memory with specified id was found or not
     */
    bool getAuxMemory(int id, unsigned int &subid, std::vector<double> &memory) const;
    /**
     * Get the counters for a subspace
     * @param[in]  subid         The id of the subspace
     * @param[out] auxCounter    The auxiliary Counter
     * @param[out] redCounter    The reduction Counter
     */
    void getCounter(int subid, 
                    unsigned int &auxCounter, unsigned int &redCounter) const;
    /**
     * Set the position of the subspace <subid>.
     */
    void setPosition(unsigned int subid, std::vector<double> &position);
    /**
     * Get the position of the subspace <subid>.
     */
    void getPosition(unsigned int subid, std::vector<double> &position) const;
    /**
     * Get Number of Spheres and Sticks
     * @param[out] numSpheres
     * @param[out] numSticks
     */
    void getNumPotentials(int &numPotentials) const;
    /**
     * return the number of subspaces
     */
    unsigned int getNumSubspaces() const {return m_subspaces.size();}
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
     * Do we need to print the auxiliary memory?
     * @return 
     */
    bool printAuxMem() const;
    /**
     * Print out the informations about the potentials for a trajectory
     */
    std::string traj_str() const;
    /**
     * Return a string about the topology of the umbrella
     */
    std::string str();
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

