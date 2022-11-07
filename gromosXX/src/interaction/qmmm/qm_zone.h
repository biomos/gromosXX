/**
 * @file qm_zone.h
 * QM zone class - implements QM zone for QM/MM calculations
 */

#ifndef INCLUDED_QM_ZONE_H
#define	INCLUDED_QM_ZONE_H

namespace math
{
  template<math::boundary_enum B>
  class Periodicity;
}

namespace interaction {
  struct QM_Atom;
  struct MM_Atom;
  class QM_Link;
  struct PairlistContainer;
  /**
   * @class QM_Zone
   * Holds all information on QM Zone - this will be passed between QM worker and QMMM interaction
   * merges function of QM storage and QM zone and eventually allows multiple QM zones
   */
  class QM_Zone {
  public:
    /**
     * Constructor
     */
    explicit QM_Zone(int net_charge = 0
                   , int spin_mult = 1
                  );

    /**
     * Copy constructor
     */
    //explicit QM_Zone(const QM_Zone& qmz);

    /**
     * Assignment operator
     */
    //QM_Zone& operator= (const QM_Zone& qmz);

    /**
     * Destructor
     */
    ~QM_Zone();
    
    /**
     * QM atoms
     */
    std::set<interaction::QM_Atom> qm;

    /**
     * MM atoms and point charges (including COS) - these are expected to be rebuilt between steps (due to cutoff)
     */
    std::set<interaction::MM_Atom> mm;

    /**
     * QMMM links
     */
    std::set<interaction::QM_Link> link;
    
    /**
     * Initialize QM zone
     */
    int init(topology::Topology& topo, 
             const configuration::Configuration& conf, 
             const simulation::Simulation& sim);
    
    /**
     * Update positions of QM and re-gather MM atoms
     */
    int update(topology::Topology& topo, 
               const configuration::Configuration& conf, 
               const simulation::Simulation& sim);

    /**
     * Write data to the configuration and topology
     */
    void write(topology::Topology& topo, 
               configuration::Configuration& conf, 
               const simulation::Simulation& sim);
    
    /**
     * Update QM-MM pairlist
     */
    void update_pairlist(const topology::Topology& topo
                       , const simulation::Simulation& sim
                       , PairlistContainer& pairlist
                       , unsigned begin, unsigned end
                       , unsigned stride) const;

    /**
     * Update charges-on-spring
     */
    void update_cos(const topology::Topology& topo
                  , const configuration::Configuration& conf
                  , const simulation::Simulation& sim);

    /**
     * Provide electric field on MM atoms
     */
    void electric_field(const simulation::Simulation& sim
                      , math::VArray& electric_field);

    /**
     * Zero energies
     */
    void inline zero();

    /**
     * Clear the QM zone
     */
    void inline clear();

    /**
     * QM energy accessor
     */
    double QM_energy() const {
      return m_qm_energy;
    }

    /**
     * QM energy mutator
     */
    double & QM_energy() {
      return m_qm_energy;
    }

    /**
     * net charge accessor
     */
    int charge() const {
      return m_charge;
    }

    /**
     * net charge mutator
     */
    int & charge() {
      return m_charge;
    }

    /**
     * spin multiplicity accessor
     */
    int spin_mult() const {
      return m_spin_mult;
    }

    /**
     * spin multiplicity mutator
     */
    int & spin_mult() {
      return m_spin_mult;
    }

    /**
     * Creates QM buffer zone
     */
    QM_Zone* create_buffer_zone(const topology::Topology& topo
                              , const simulation::Simulation& sim);

  protected:
    /**
     * the QM energy of the zone
     */
    double m_qm_energy;

    /**
     * the QM zone net charge
     */
    int m_charge;

    /**
     * the QM zone net spin multiplicity
     */
    int m_spin_mult;

    /**
     * Scale charges with distance to the nearest QM atom
     */
    void scale_charges(const simulation::Simulation& sim);

    /**
     * Gather QM atoms
     */
    int get_qm_atoms(const topology::Topology& topo, 
                     const configuration::Configuration& conf, 
                     const simulation::Simulation& sim);
    
    /**
     * Gather buffer atoms in adaptive scheme
     */
    void get_buffer_atoms(topology::Topology& topo, 
                          const configuration::Configuration& conf, 
                          const simulation::Simulation& sim);

    /**
     * Update positions of QM atoms - wrapper
     */
    void update_qm_pos(const topology::Topology& topo, 
                       const configuration::Configuration& conf, 
                       const simulation::Simulation& sim);

    /**
     * Gather MM atoms - wrapper
     */
    void get_mm_atoms(const topology::Topology& topo, 
                      const configuration::Configuration& conf, 
                      const simulation::Simulation& sim);
    
    /**
     * Get QM-MM bonds
     */
    void get_links(const topology::Topology& topo, 
                   const simulation::Simulation& sim);
    
    /**
     * Update capping atoms positions
     */
    void update_links(const simulation::Simulation& sim);

  private:  
    /**
     * Gather chargegroups within cutoff from QM atoms
     * QM buffer atoms (QM_Atom) or MM atoms (MM_Atom)
     */
    template <math::boundary_enum B, class AtomType>
    int gather_chargegroups(const topology::Topology& topo, 
                            const configuration::Configuration& conf, 
                            const simulation::Simulation& sim,
                            std::set<AtomType>& atom_set,
                            const double cutoff2);
    
    /**
     * Gather buffer atoms in adaptive scheme - internal function
     */
    template<math::boundary_enum B>
    int _get_buffer_atoms(topology::Topology& topo, 
                          const configuration::Configuration& conf, 
                          const simulation::Simulation& sim);

    /**
     * Update positions of QM atoms - internal function
     */
    template<math::boundary_enum B>
    int _update_qm_pos(const topology::Topology& topo, 
                       const configuration::Configuration& conf, 
                       const simulation::Simulation& sim);

    /**
     * Gather MM atoms (chargegroup-based cutoff) - internal function
     */
    template<math::boundary_enum B>
    int _get_mm_atoms(const topology::Topology& topo, 
                      const configuration::Configuration& conf, 
                      const simulation::Simulation& sim);

    /**
     * Gather MM atoms (atom-based cutoff) - internal function
     */
    template<math::boundary_enum B>
    int _get_mm_atoms_atomic(const topology::Topology& topo, 
                             const configuration::Configuration& conf, 
                             const simulation::Simulation& sim);

    /**
     * Gather linked MM atoms - internal function
     */
    template<math::boundary_enum B>
    void get_linked_mm_atoms(const topology::Topology& topo, 
                            const configuration::Configuration& conf,
                            const math::Periodicity<B>& periodicity);

    /**
     * Emplace the atom of proper type to the provided set - internal function
     */
    template<class AtomType>
    inline void emplace_atom(std::set<AtomType>& set,
                             typename std::set<AtomType>::const_iterator& it,
                             const unsigned index,
                             const math::Vec& pos,
                             const unsigned atomic_number,
                             const double charge,
                             const bool is_qm_buffer);

    /**
     * Skip the chargegroup search - internal function
     */
    template<class AtomType>
    inline bool skip_cg(const topology::Topology& topo,
                        const unsigned index);
  };

  /**
   * Emplace the QM atom
   */
  template<>
  inline void QM_Zone::emplace_atom<interaction::QM_Atom>
                     (std::set<interaction::QM_Atom>& set,
                      std::set<interaction::QM_Atom>::const_iterator& it,
                      const unsigned index,
                      const math::Vec& pos,
                      const unsigned atomic_number,
                      const double charge,
                      const bool is_qm_buffer);

  /**
   * Emplace the MM atom
   */
  template<>
  inline void QM_Zone::emplace_atom<interaction::MM_Atom>
                     (std::set<interaction::MM_Atom>& set,
                      std::set<interaction::MM_Atom>::const_iterator& it,
                      const unsigned index,
                      const math::Vec& pos,
                      const unsigned atomic_number,
                      const double charge,
                      const bool is_qm_buffer);

  /**
   * Skip the chargegroup search - QM_Atom version
   */
  template<>
  inline bool QM_Zone::skip_cg<interaction::QM_Atom>
                              (const topology::Topology& topo,
                               const unsigned index);

  /**
   * Skip the chargegroup search - MM_Atom version
   */
  template<>
  inline bool QM_Zone::skip_cg<interaction::MM_Atom>
                              (const topology::Topology& topo,
                               const unsigned index);
}
#endif	/* QM_ZONE_H */

