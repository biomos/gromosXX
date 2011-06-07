/**
 * @file qm_storage.h
 * storage for the QM workers
 */

#ifndef QM_STORAGE_H
#define	QM_STORAGE_H

namespace interaction {
  /**
   * @class QM_Storage
   * storage for the results of the QM calculation
   */
  class QM_Storage {
  public:
    /**
     * Constructor
     */
    QM_Storage() : energy(0), force(0), charge(0) {}
    /**
     * the energy of the system
     */
    double energy;
    /**
     * the force acting on every atom
     */
    math::VArray force;
    /**
     * the force acting on every charge-on-spring
     */
    math::VArray cos_force;
    /**
     * the charge of every atom
     */
    math::SArray charge;
    /**
     * zero the storage
     */
    void zero() {
      energy = 0.0; force = 0.0; cos_force = 0.0; charge = 0.0; 
    }
    /**
     * resize the storage
     * @param num_atoms number of atoms
     * @param num_qm_atoms number of QM atoms
     */
    void resize(unsigned int num_atoms, unsigned int num_qm_atoms) {
      force.resize(num_atoms, math::Vec(0.0, 0.0, 0.0));
      cos_force.resize(num_atoms, math::Vec(0.0, 0.0, 0.0));
      charge.resize(num_qm_atoms, 0.0);
    }
  };
}

#endif	/* QM_STORAGE_H */

