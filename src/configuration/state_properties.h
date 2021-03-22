/**
 * @file state_properties.h
 * calculates state properties.
 */

#ifndef INCLUDED_STATE_PROPERTIES_H
#define INCLUDED_STATE_PROPERTIES_H

namespace configuration
{
  /**
   * @class State_Properties
   * calculates state properties.
   */
  class State_Properties
  {
  public:
    State_Properties(Configuration const & conf) : m_configuration(conf) {};
    
    /**
     * calculate the center of mass and the
     * translational kinetic energy of a group of
     * atoms.
     * @param start begin of a group of atoms.
     * @param end of a group of atoms.
     * @param mass the masses of (all) atoms.
     * @param com_pos returns the center of mass.
     * @param com_e_kin returns the tranlational kinetic energy tensor.
     * 
     * @todo the gathering of the molecule is hardcoded in here.
     * Maybe this should be changed to a generic implementation.
     * Gathering is done in respect to the previous atom. An idea would
     * be to gather as default with respect to the previous atom but
     * letting the user override this (GATHER block).
     * This does not yield the same answer as Phils approach for all cases
     * but maybe for the practical ones???
     */
    void center_of_mass(topology::Atom_Iterator start, 
			topology::Atom_Iterator end,
			math::SArray const &mass, 
			math::Vec &com_pos, math::Matrix &com_e_kin);
    
    /**
     * calculate molecular translational kinetic energy.
     * @param start begin of molecule.
     * @param end end of molecule.
     * @param mass the atomic masses.
     * @param mol_v delivered with the average 
     *        molecular center of mass velocity.
     * @param com_e_kin delivered with molecular center of mass 
     *        kinetic energy (-dt/2).
     * @param e_kin average kinetic energy
     * @param new_mol_v delivered with molecular center 
     *        of mass velocity (+dt/2).
     * @param new_com_e_kin delivered with molecular center of mass kinetic energy (+dt/2).
     * @param new_e_kin the kinetic energy at (+dt/2).
     */
    void molecular_translational_ekin(simulation::Simulation &sim,
                                          topology::Atom_Iterator start, 
				      topology::Atom_Iterator end,
				      math::SArray const &mass, 
				      math::Vec &mol_v, 
				      double &com_e_kin,
				      double &e_kin,
				      math::Vec &new_mol_v,
				      double &new_com_e_kin,
				      double &new_e_kin);

  protected:
    /**
     * reference to the configuration.
     */
    Configuration const & m_configuration;
    

  };
  
}

#endif
