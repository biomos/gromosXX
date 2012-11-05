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
  
  /**
   * @class BS_Subspace
   * @ingroup util
   * 
   * Holds its definition via BS_Coordinates and the Potentials.
   */
  class BS_Subspace {
  public:
    /**
     * The constructor of a subspace
     * 
     * @param id                The ID of the subspace
     * @param forceIncrement    The force increment factor for the subspace
     * @param reductionFactor   Reduce forceIncrement by this factor, when the
     *                          auxiliary counter is larger than the global cutoff
     * @param localCutoff       The local cutoff
     * @param globalCutoff      The global cutoff
     */
    BS_Subspace(int id, double forceIncrement, double reductionFactor, 
                int localCutoff, int globalCutoff) :
            id(id), m_numSpheres(0), m_numSticks(0), m_numSnakes(0), m_numPipes(0),
            m_forceIncrement(forceIncrement), m_reductionFactor(reductionFactor), 
            m_localCutoff(localCutoff), m_globalCutoff(globalCutoff),
            m_auxilliaryCounter(0), m_reductionCounter(0) {}

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
     * @param[in]    potentials         A vector with all the potentials
     * @param[in]    beta               1 / kB T
     */
    void addForces(configuration::Configuration &conf,
            //double totalPartitionFct,
            std::vector<double> &potentials,
            double beta);
    /**
     * Calculate the potentials and update the current position in the active
     * subspace and return the total Partition Function of the potentials.
     * @param conf The current configuration
     * @return The total Partition Function of all the potentials
     */
    void calculatePotential(configuration::Configuration &conf,
            std::vector<double> &potentials);
    /**
     * Update the memory of the potentials
     */
    void updateMemory();
    /**
     * Set the memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[in] memory    The memory
     * @return whether Potential was found
     */
    bool setMemory(int id, std::vector<double> &memory);
    /**
     * Get the memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[out] memory    The memory
     * @return whether Potential was found
     */
    bool getMemory (int id, std::vector<double> &memory) const;
    /**
     * Set all the memories to zero.
     */
    void setMemoryToZero();
    /**
     * Set the auxiliary memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[in] memory    The memory
     * @return whether Potential was found
     */
    bool setAuxMemory(int id, std::vector<double> &memory);
    /**
     * Get the auxiliary memory of Potential
     * @param[in] id    The id of the potential
     * @param[in] type  The Type (Sphere / Stick)
     * @param[out] memory    The auxiliary memory
     * @return whether Potential was found
     */
    bool getAuxMemory (int id, std::vector<double> &memory) const;
    /**
     * Set all the auxiliary memories to zero.
     */
    void setAuxMemoryToZero();
    /**
     * Set the auxilliary and reduction counter
     * @param auxillaryCounter
     * @param reductionCounter
     */
    void setCounter(unsigned int auxilliaryCounter, 
                     unsigned int reductionCounter);
    /**
     * Set the auxilliary and reduction counter
     * @param[inout] auxillaryCounter
     * @param[inout] reductionCounter
     */
    void getCounter(unsigned int &auxilliaryCounter, 
                     unsigned int &reductionCounter);
    /**
     * Set the position in the subspace.
     */
    void setPosition(std::vector<double> &position);
    /**
     * Get the position in the subspace.
     */
    void getPosition(std::vector<double> &position);
    /**
     * Return the number of spheres
     */
    int getNumSpheres() const {return m_numSpheres;}
    /**
     * Return the number of sticks.
     */
    int getNumSticks() const {return m_numSticks;}
    /**
     * Return the number of snakes.
     */
    int getNumSnakes() const {return m_numSnakes;}
    /**
     * Return the number of snakes.
     */
    int getNumPipes() const {return m_numPipes;}
    /**
     * Return the number of potentials.
     */
    int getNumPotentials() const {
      return m_numSnakes + m_numSticks + m_numSpheres + m_numPipes;
    }
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
     * Return the number of dimensions of each coordinate  in the subspace.
     */
    std::vector<int> getDimensionality();
    /**
     * Return the Total Number of Dimensions of the subspace
     * (sum of all elements in getDimensionality)
     */
    unsigned int getNumDimensions();
    /**
     * Test wheter coordinate "id" has type "type".
     */
    bool testType(int id, BS_Coordinate::Coord_type type);
    /**
     * Do we need to print the auxiliary memory for this subspace?
     */
    bool printAuxMem(){ return !(m_reductionFactor == 1.0);}
    /**
     * A string describing the subspace
     * @return 
     */
    std::string str();
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
    std::vector<BS_Potential *> m_potentials;
    /**
     * The number of Spheres
     */
    int m_numSpheres;
    /**
     * The number of Sticks
     */
    int m_numSticks;
    /**
     * The number of Snakes
     */
    int m_numSnakes;
    /**
     * The number of Pipes
     */
    int m_numPipes;
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
    unsigned int m_globalCutoff;
    /**
     * Auxilliary Counter: N_c(t)
     */
    unsigned int m_auxilliaryCounter;
    /**
     * The reduction Counter: I_R 
     */
    unsigned int m_reductionCounter;
    /** 
     * Returns a string with informations of state of the subspace
     */
    std::string debug_str();
  };
}
#endif	/* BS_SUBSPACE_H */

