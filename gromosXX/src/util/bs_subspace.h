/**
 * @file   bs_subspace.h
 */

#ifndef BS_SUBSPACE_H
#define	BS_SUBSPACE_H
#include "bs_coordinate.h"
#include "bs_potentials.h"

namespace util {

  class BS_Coordinate;
  class BS_Vector;
  
  class BS_Subspace {
  public:
    BS_Subspace(int id, double forceIncrement, double reductionFactor, 
                int localCutoff, int globalCutoff) :
            id(id), m_numSpheres(0), m_numSticks(0),
            m_forceIncrement(forceIncrement), m_reductionFactor(reductionFactor), 
            m_localCutoff(localCutoff), m_globalCutoff(globalCutoff){}

    ~BS_Subspace(){}
    /**
     * Free the potentials and coordinates from the memory.
     */
    void freeMemory();
    /**
     * transform the current configuration into internal Configurations
     * @param[in]    conf   The current configuration
     * @param[inout] result The vector containing the internal coordinates
     *                      (Q)
     */
    void transformCoordinates(configuration::Configuration &conf,
            BS_Vector &result);
    /**
     * Update the current position in internal coordinates
     * @param conf
     */
    void transformCurrent(configuration::Configuration &conf);
    /**
     * Add the Forces of one Potential (in internal Coordinates) to the 
     * atoms.
     * @param[inout] conf               The current configuration
     * @param[in]    totalPartitionFct  The total Partition Function for
     *                                  calculating the weight
     */
    void addForces(configuration::Configuration &conf,
            double totalPartitionFct);
    /**
     * Calculate the potentials and update the current position in the active
     * subspace and return the total Partition Function of the potentials.
     * @param conf The current configuration
     * @return The total Partition Function of all the potentials
     */
    double calculatePotential(configuration::Configuration &conf,
            simulation::Simulation &sim);
    /**
     * Update the memory of the potentials
     */
    void updateMemory();
    /**
     * Set the memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[in] memory    The memory
     */
    void setMemory(int id, BS_Potential::potential_enum type, 
            std::vector<double> &memory);
    /**
     * Get the memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[out] memory    The memory
     * @return whether Potential was found
     */
    bool getMemory (int id, BS_Potential::potential_enum type,
            std::vector<double> &memory) const;
    /**
     * Set all the memories to zero.
     */
    void setMemoryToZero();
    /**
     * Set the auxiliary memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[in] memory    The memory
     * @param[in] auxCounter    The auxiliary memory Counter
     * @param[in] redCounter    The reduction Counter
     */
    void setAuxMemory(int id, BS_Potential::potential_enum type, 
            std::vector<double> &memory, int auxCounter, int redCounter);
    /**
     * Get the auxiliary memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[out] memory    The auxiliary memory
     * @param[out] auxCounter    The auxiliary memory Counter
     * @param[out] redCounter    The reduction Counter
     * @return whether Potential was found
     */
    bool getAuxMemory (int id, BS_Potential::potential_enum type,
            std::vector<double> &memory, int &auxCounter, int &redCounter) const;
    /**
     * Set all the auxiliary memories to zero.
     */
    void setAuxMemoryToZero();
    /**
     * Return the number of spheres
     */
    int getNumSpheres() const {return m_numSpheres;}
    /**
     * Return the number of sticks.
     */
    int getNumSticks() const {return m_numSticks;}
    /**
     * Add a new coordinate to the definition of the subspace
     * @param newCoordinate
     */
    void addCoordinate(BS_Coordinate *newCoordinate);
    /**
     * Add a new potential to define the active subspace
     * @param newPotential
     */
    void addPotential(BS_Potential *newPotential);
    /**
     * Return the force in internal coordinates
     * @param[inout] force
     */
    void getForce(std::vector<double> &force);
    /**
     * Return the center of the sphere id
     * @param id the id of the sphere
     * @return  the center
     */
    BS_Vector getCenter(int id);
    /**
     * A string describing the subspace
     * @return 
     */
    std::string str() {return debug_str();}
    /**
     * Return a string with all informations of the potential for the
     * trajectory
     * @return string with the potentials
     */
    std::string traj_str();
    /**
     * The id of the subspace
     */
    int id;
  private:
    /**
     * The definition of the subspace
     */
    std::vector<BS_Coordinate *> m_definition;
    /**
     * The current position expressed in internal coordinates
     */
    BS_Vector bs_pos;
    /**
     * The potentials of this subspace and the ones with a weight > 0
     */
    std::vector<BS_Potential *> potentials;
    /**
     * The number of Spheres
     */
    int m_numSpheres;
    /**
     * The number of Sticks
     */
    int m_numSticks;
    /**
     * The force in the subspace
     */
    BS_Vector m_force;
    /**
     * Basis force constant increment: k_LE
     */
    double m_forceIncrement;
    /**
     * The force constant reduction factor: f_LE
     */
    double m_reductionFactor;
    /**
     * The local visiting cutoff: gamma_LE
     */
    int m_localCutoff;
    /**
     * The global visiting cutoff: n_LE
     */
    int m_globalCutoff;
    
    /** 
     * Returns a string with informations of state of the subspace
     */
    std::string debug_str();
  };
}
#endif	/* BS_SUBSPACE_H */

