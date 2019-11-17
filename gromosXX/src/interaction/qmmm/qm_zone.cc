/**
 * @file qm_zone.cc
 * Implements methods of QM_Zone
 */
#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <interaction/nonbonded/pairlist/pairlist.h>

#include <math/periodicity.h>
#include <util/template_split.h>

#include <util/debug.h>

#include <interaction/qmmm/qm_atom.h>
#include <interaction/qmmm/mm_atom.h>
#include <interaction/qmmm/qm_link.h>
#include <interaction/qmmm/qm_zone.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::QM_Zone::QM_Zone(int net_charge
                            , int spin_multiplicity)
                              : qm_energy(0.0)
                              , mm_energy(0.0)
{}

interaction::QM_Zone::~QM_Zone()
{}

void interaction::QM_Zone::zero() {
  this->qm_energy = 0.0; this->mm_energy = 0.0; 
}

void interaction::QM_Zone::clear() {
  this->zero();
  this->qm.clear();
  this->mm.clear();
  this->link.clear();
}

int interaction::QM_Zone::init(const topology::Topology& topo, 
                               const configuration::Configuration& conf, 
                               const simulation::Simulation& sim) {
  int err;
  
  DEBUG(15,"Getting QM atoms");
  err = this->get_qm_atoms(topo,conf,sim);
  if (err) return err;

  DEBUG(15,"Getting MM atoms");
  err = this->get_mm_atoms(topo,conf,sim);
  if (err) return err;

  DEBUG(15,"Getting QM-MM links");
  this->get_links(topo, sim);
  return 0;
}

int interaction::QM_Zone::update(const topology::Topology& topo, 
                                 const configuration::Configuration& conf, 
                                 const simulation::Simulation& sim) {
  this->zero();
  DEBUG(15,"Updating QM atoms positions");
  int err = this->update_qm_pos(topo, conf, sim);
  if (err) return err;
  
  DEBUG(15,"Getting MM atoms positions");
  err = this->get_mm_atoms(topo, conf, sim);
  if (err) return err;
  
  DEBUG(15,"Updating links");
  this->update_links(sim);
  return 0;
}

void interaction::QM_Zone::write_pos(math::VArray& pos) {
  /**
   * Check, if changed (if optimization took place),
   * otherwise we may add some noise from truncated precision
  */
  for (std::set<QM_Atom>::const_iterator
      it = this->qm.begin(), to = this->qm.end(); it != to; ++it)
    {
    pos(it->index) = it->pos;
  }
}

void interaction::QM_Zone::write_force(math::VArray& force) {
  // First distribute link atom forces
  for (std::set<QM_Link>::const_iterator
      it = this->link.begin(), to = this->link.end(); it !=to; ++it)
    {
    it->distribute_force(*this);
  }

  // Then write QM atoms forces
  for (std::set<QM_Atom>::const_iterator
      it = this->qm.begin(), to = this->qm.end(); it != to; ++it)
    {
    force(it->index) += it->force;
  }
  // And MM atoms forces
  for (std::set<MM_Atom>::const_iterator
      it = this->mm.begin(), to = this->mm.end(); it != to; ++it)
    {
    force(it->index) += it->force;
  }
}

void interaction::QM_Zone::write_charge(math::SArray& charge) {
  for (std::set<QM_Atom>::const_iterator
    it = this->qm.begin(), to = this->qm.end(); it != to; ++it)
    {
    charge(it->index) = it->qm_charge;
  }
}

int interaction::QM_Zone::get_qm_atoms(const topology::Topology& topo, 
                                       const configuration::Configuration& conf, 
                                       const simulation::Simulation& sim) {
  this->qm.clear();
  DEBUG(15,"QM_Zone::get_qm_atoms: Collecting atoms");
  // Check if one CG is not split between QM and MM
  // Also collect QM chargegroups
  const unsigned num_solute_atoms = topo.num_solute_atoms();
  const unsigned num_solute_cg = topo.num_solute_chargegroups();
  DEBUG(15, "num_solute_atoms: " << num_solute_atoms);
  DEBUG(15, "num_solute_cg: " << num_solute_cg);
  for (unsigned cg = 0; cg < num_solute_cg; ++cg) {
    unsigned a = topo.chargegroup(cg);
    const unsigned a_to = topo.chargegroup(cg + 1);
    DEBUG(15, "cg: " << cg);
    DEBUG(15, "a_to: " << a_to);
    assert(a_to <= num_solute_atoms);
    const bool cg_is_qm = topo.is_qm(a);
    for (; a < a_to; ++a) {
      const bool a_is_qm = topo.is_qm(a);
      if (cg_is_qm != a_is_qm) {
        std::ostringstream msg;
        msg << "Chargegroup " << cg << " is split between QM and MM zone - "
            << "atoms " << a + 1 << " and " << a;
        io::messages.add(msg.str(), "QM Zone", io::message::error);
        return 1;
      }
      if (a_is_qm) {
        if (topo.is_polarisable(a)
                || topo.is_coarse_grained(a)
                || topo.is_perturbed(a)
                || topo.is_eds_perturbed(a)) {
          std::ostringstream msg;
          msg << "QM atom should not be polarisable, coarse-grained or perturbed";
          io::messages.add(msg.str(), "QM Zone", io::message::error);
          return 1;
        }
        DEBUG(7,"Adding QM atom " << a << ", atomic number: " << topo.qm_atomic_number(a));
        this->qm.emplace(a, topo.qm_atomic_number(a));
      }
    }
  }
  if (this->qm.size() == 0) {
    io::messages.add("QMMM requested but no QM atoms found", "QM Zone", io::message::error);
    return 1;
  }
  DEBUG(15,"QM_Zone::get_qm_atoms: Updating QM atoms positions");
  int err = this->update_qm_pos(topo, conf, sim);
  if (err) return err;
  return 0;
}

int interaction::QM_Zone::update_qm_pos(const topology::Topology& topo, 
                                        const configuration::Configuration& conf, 
                                        const simulation::Simulation& sim) {
  DEBUG(15,"QM_Zone::update_qm_pos: Splitting boundary");
  SPLIT_BOUNDARY(this->_update_qm_pos, topo, conf, sim);
  return 0;
}

template<math::boundary_enum B>
int interaction::QM_Zone::_update_qm_pos(const topology::Topology& topo, 
                                         const configuration::Configuration& conf, 
                                         const simulation::Simulation& sim) {
  math::Periodicity<B> periodicity(conf.current().box);
  const math::VArray& pos = conf.current().pos;

  // Get position of first QM atom, will be reference for gathering
  const std::set<QM_Atom>::iterator first_it = this->qm.begin();

  assert(first_it != this->qm.end());

  DEBUG(15, "First QM index: " << first_it->index);
  const math::Vec& ref_pos = pos(first_it->index);

  DEBUG(15, "First QM as ref_pos: " << math::v2s(ref_pos));
  first_it->pos = ref_pos;
  math::Vec nim;

  // Gather the rest of QM atoms
  if (this->qm.size() > 1) {
    for (std::set<QM_Atom>::iterator qm_it = std::next(first_it),
        qm_to = this->qm.end(); qm_it != qm_to; ++qm_it)
      {
      periodicity.nearest_image(ref_pos, pos(qm_it->index), nim);
      qm_it->pos = ref_pos - nim;
      DEBUG(15, "QM atom " << qm_it->index << " : " << math::v2s(qm_it->pos));
    }
  }
  // Check if QM zone sees its periodic image
  if (conf.boundary_type != math::vacuum && this->qm.size() > 1) {
    const std::set<QM_Atom>::const_iterator to1 = std::prev(this->qm.end())
                                          , to2 = this->qm.end();
    std::set<QM_Atom>::const_iterator it1, it2;
    DEBUG(15, "QM zone gathering check");
    for (it1 = this->qm.begin(); it1 != to1; ++it1) {
      const math::Vec& i_pos = it1->pos;
      DEBUG(15, "Atom " << it1->index << ": " << math::v2s(i_pos));
      for (it2 = std::next(it1); it2 != to2; ++it2) {
        const math::Vec& j_pos = it2->pos;
        math::Vec nim;
        periodicity.nearest_image(i_pos, j_pos, nim);
        DEBUG(15, "nim to " << it2->index << " : " << math::v2s(nim));
        const math::Vec j_pos_2 = i_pos - nim;
        DEBUG(15, "j_pos:   " << math::v2s(j_pos));
        DEBUG(15, "j_pos_2: " << math::v2s(j_pos_2));
        const double delta = math::abs2(j_pos_2 - j_pos);
        DEBUG(15, "delta:   " << delta);
        if (delta > math::epsilon) {
          std::ostringstream msg;
          msg << "QM zone sees own periodic image (atoms "
              << (it1->index + 1) << " and " << (it2->index + 1) << ")";
          io::messages.add(msg.str(), "QM Zone", io::message::error);
        }
      }
    }
  }
  return 0;
}

void interaction::QM_Zone::update_pairlist(
                  const topology::Topology& topo
                , const simulation::Simulation& sim
                , PairlistContainer& pairlist
                , unsigned begin, unsigned end
                , unsigned stride) const {
  pairlist.clear();
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    // Generate QM-MM pairs
    unsigned counter = begin;
    for (std::set<QM_Atom>::const_iterator
          qm_it = this->qm.begin(), qm_to = this->qm.end();
          qm_it != qm_to; ++qm_it) {
      for (std::set<MM_Atom>::const_iterator
            mm_it = this->mm.begin(), mm_to = this->mm.end();
            mm_it != mm_to; ++mm_it) {
        if (counter == 0) {
          unsigned i = qm_it->index;
          unsigned j = mm_it->index;
          if (i > j) std::swap(i,j);
          // There should not be any QM pair in standard exclusions
          assert(!topo.all_exclusion(i).is_excluded(j));
          if (!topo.qm_all_exclusion(i).is_excluded(j)) {
            pairlist.solute_short[i].push_back(j);
          }
          counter = stride;
        }
        else {
          --counter;
        }
      }
    }
  }
  if (sim.param().qmmm.qm_lj && this->qm.size() > 1) {
    // Generate QM-QM pairs
    unsigned counter = begin;
    DEBUG(15, "Updating QM-QM pairlist");
    const std::set<QM_Atom>::const_iterator qm_to1 = std::prev(this->qm.end())
                                          , qm_to2 = this->qm.end();
    std::set<QM_Atom>::const_iterator qm_it1, qm_it2;
    for (qm_it1 = this->qm.begin(); qm_it1 != qm_to1; ++qm_it1) {
      for (qm_it2 = std::next(qm_it1); qm_it2 != qm_to2; ++qm_it2) {
        if (counter == 0) {
          unsigned i = qm_it1->index;
          unsigned j = qm_it2->index;
          DEBUG(15, "Adding QM-QM pair: " << i << " - " << j);
          if (!topo.all_exclusion(i).is_excluded(j)
              && !topo.qm_all_exclusion(i).is_excluded(j)) {
            pairlist.solute_short[i].push_back(j);
          }
          counter = stride;
        }
        else {
          --counter;
        }
      }
    }
    DEBUG(15, "QM-QM pairlist done");
  }
}

int interaction::QM_Zone::get_mm_atoms(const topology::Topology& topo, 
                                       const configuration::Configuration& conf, 
                                       const simulation::Simulation& sim) {
  /** Could be probably done employing standard pairlist, but we are
   * using different cutoff. We also do gathering at the same time and perform
   * some checks on cutoff and periodic copy interactions so it becomes
   * rather specific for QM and useless for other purposes, thus could
   * be done here
   */
  this->mm.clear();
  if (sim.param().boundary.boundary == math::vacuum
      && sim.param().qmmm.qmmm == simulation::qmmm_mechanical) {
    // Include only link atoms
    DEBUG(9, "Gathering only linked MM atoms");
    for (std::set< std::pair<unsigned,unsigned> >::const_iterator
        it = topo.qmmm_link().begin(), to = topo.qmmm_link().end();
        it != to; ++it) {
      const unsigned i = it->second;
      DEBUG(9, "Adding MM link atom " << i);
      this->mm.emplace(i, topo.charge(i), conf.current().pos(i));
    }
  }
  else if (sim.param().boundary.boundary == math::vacuum
            && sim.param().qmmm.qmmm != simulation::qmmm_mechanical
            && sim.param().qmmm.cutoff == 0.0) {
    // Include all - allowed only with vacuum PBC
    DEBUG(9, "Gathering all MM atoms");
    for (unsigned int i = 0; i < topo.num_atoms(); ++i) {
      if (!topo.is_qm(i)) {
        DEBUG(9, "Adding atom " << i);
        this->mm.emplace(i, topo.charge(i), conf.current().pos(i));
      }
    }
  }
  else {
    // Apply chargegroup- or atom-based cutoff
    if (sim.param().qmmm.atomic_cutoff) {
      SPLIT_BOUNDARY(this->_get_mm_atoms_atomic, topo, conf, sim);
    }
    else {
      SPLIT_BOUNDARY(this->_get_mm_atoms, topo, conf, sim);
    }
  }
  if (sim.param().qmmm.qmmm == simulation::qmmm_polarisable)
    this->update_cos(topo, conf, sim);
  if (sim.param().qmmm.mm_scale >= 0.0)
    this->scale_charges(sim);
  return 0;
}

template<math::boundary_enum B>
int interaction::QM_Zone::_get_mm_atoms(const topology::Topology& topo, 
                                        const configuration::Configuration& conf, 
                                        const simulation::Simulation& sim)
  {
  math::Periodicity<B> periodicity(conf.current().box);
          DEBUG(15, "Boundary type: " << conf.boundary_type);
  // MM link atoms need to be always gathered
  this->get_linked_mm_atoms<B>(topo, conf, periodicity);
  // Gather MM atoms only for electrostatic or polarisable embedding
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    DEBUG(9, "Gathering MM atoms using chargegroup-based cutoff");
    /** We are running twice almost the same loop for solute and solvent.
     * Lambdas cannot be used because passing capturing lambda to lambda
     * is not allowed. We can solve this by passing local functors to lambda.
     */
      class cgs_funct {
      public:
        virtual ~cgs_funct() {};
        virtual const math::Vec& operator()(unsigned i) const = 0;
      };

      class solute_cgs_funct : public cgs_funct {
        const math::VArray& m_cogs;
      public:
        solute_cgs_funct(const math::VArray& cogs)
                  : m_cogs(cogs) {};
        //~solute_cgs_funct();
        const math::Vec& operator()(unsigned i) const {
          return m_cogs(i);
        }
      };

      class solvent_cgs_funct : public cgs_funct {
        const topology::Topology& m_topo;
        const math::VArray& m_pos;
      public:
        solvent_cgs_funct(const topology::Topology& topo, const math::VArray& pos)
                  :m_topo(topo), m_pos(pos) {};
        const math::Vec& operator()(unsigned i) const {
          return m_pos(m_topo.chargegroup(i));
        }
      };
    // End of functors definition

    // Initialize constants
    const math::VArray& pos = conf.current().pos;
    const double cutoff2 = sim.param().qmmm.cutoff * sim.param().qmmm.cutoff;
    DEBUG(15, "cutoff2 = " << cutoff2);

    const unsigned num_cg = topo.num_chargegroups()
                 , num_solute_cg = topo.num_solute_chargegroups();
    DEBUG(15, "num_cg = " << num_cg);
    DEBUG(15, "num_solute_cg = " << num_solute_cg);

    // Calculate solute CG COGs
    math::VArray cogs(num_solute_cg);
    {
      topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin();
      for (unsigned cg = 0; cg < num_solute_cg; ++cg_it, ++cg) {
        cg_it.cog(pos, cogs(cg));
      }
    }
    // Functor for solute CG COGs
    cgs_funct* solute_cgs = new solute_cgs_funct(cogs);
    // Functor for solvent CG COGs
    cgs_funct* solvent_cgs = new solvent_cgs_funct(topo, pos);

    // Collect MM chargegroups within every QM atom cutoff
    for (std::set<QM_Atom>::const_iterator
        qm_it = this->qm.begin(), qm_to = this->qm.end(); qm_it != qm_to; ++qm_it)
      {
      DEBUG(15, "Gathering around QM atom " << qm_it->index);
      const math::Vec& qm_pos = qm_it->pos;
      math::Vec r_qm_cg;
      std::set<MM_Atom>::iterator mm_it = this->mm.begin();
      // We run almost the same loop twice
      // With functors and lambda, one definition is enough
      auto cg_cogs_loop = [&](unsigned cg_from, unsigned cg_to, cgs_funct* cgs)-> int {
        for(unsigned cg = cg_from; cg < cg_to; ++cg) {
          // If first atom is QM, then whole chargegroup is, skip
          unsigned a = topo.chargegroup(cg);
          if (topo.is_qm(a)) continue;
          // get distance
          periodicity.nearest_image(qm_pos, (*cgs)(cg), r_qm_cg);
          DEBUG(15, "qm_pos: " << math::v2s(qm_pos));
          DEBUG(15, "(*cgs)(cg): " << math::v2s((*cgs)(cg)));
          DEBUG(15, "Nearest image qm_pos - cg " << cg << " : " << math::v2s(r_qm_cg));
          DEBUG(15, "math::abs2(r_qm_cg) = " << math::abs2(r_qm_cg));
          
          if (math::abs2(r_qm_cg) < cutoff2) {
            // Iterate over atoms in cg
            for (unsigned a_to = topo.chargegroup(cg + 1); a < a_to; ++a) {
              DEBUG(15, "Adding MM atom " << a);
              // nearest image of atom to CG COG
              math::Vec r_mm_cg;
              periodicity.nearest_image(pos(a), (*cgs)(cg), r_mm_cg);
              DEBUG(15, "nearest image mm_pos-cg_cog = " << math::v2s(r_mm_cg));
              DEBUG(15, "math::abs2(r_mm_cg) = " << math::abs2(r_mm_cg));
              // MM atom position
              const math::Vec mm_pos = qm_pos - r_qm_cg + r_mm_cg;
              /** this will add every MM atom only once. We have to check, if we already
               * have the atom and if its position is the same. Otherwise we see its
               * periodic copy and that is bad.
               */
              DEBUG(15, "mm_pos: " << math::v2s(mm_pos));
              const size_t size = this->mm.size();
              if (mm_it != this->mm.end())
                std::advance(mm_it,1); // hint on set insertion, improves performance
              mm_it = this->mm.emplace_hint(mm_it, a, topo.charge(a), mm_pos);

              if (size == this->mm.size()) {
                DEBUG(15, "Atom " << a << " already in the list");
                const double delta = math::abs2(mm_it->pos - mm_pos);
                DEBUG(15, "mm_it->pos: " << math::v2s(mm_it->pos));
                DEBUG(15, "mm_pos:     " << math::v2s(mm_pos));
                DEBUG(15, "delta: " << delta);
                if (delta > math::epsilon) {
                  // Atom already in the list but as different periodic image, exit
                  std::ostringstream msg;
                  msg << "QM atom " << (qm_it->index + 1)
                      << " sees different periodic copy of MM atom "
                      << (mm_it->index + 1);
                  io::messages.add(msg.str(), "QM Zone", io::message::error);
                  //return 1;
                } // if not in the same position
              } // if already in the set
            } // for atoms in cg
          } // if cg within cutoff
        } // for cgs in the range
        return 0;
      }; // lambda cg_cogs_loop
      int err;

      DEBUG(15, "Iterating over solute cgs");
      err = cg_cogs_loop(0, num_solute_cg, solute_cgs);
      if (err) return err;

      DEBUG(15, "Iterating over solvent cgs");
      err = cg_cogs_loop(num_solute_cg, num_cg, solvent_cgs);
      if (err) return err;
    }
    delete solute_cgs;
    delete solvent_cgs;
  }
  return 0;
}

template<math::boundary_enum B>
int interaction::QM_Zone::_get_mm_atoms_atomic(const topology::Topology& topo, 
                                               const configuration::Configuration& conf, 
                                               const simulation::Simulation& sim)
  {
  math::Periodicity<B> periodicity(conf.current().box);

  // MM link atoms need to be always gathered
  this->get_linked_mm_atoms<B>(topo, conf, periodicity);

  // Gather MM atoms only for electrostatic or polarisable embedding
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    DEBUG(9, "Gathering MM atoms using atomic cutoff");
    const math::VArray& pos = conf.current().pos;

    // Initialize constants
    const double cutoff2 = sim.param().qmmm.cutoff * sim.param().qmmm.cutoff;
    DEBUG(15, "cutoff2 = " << cutoff2);

    const unsigned num_atoms = topo.num_atoms();

    // Collect MM atoms within every QM atom cutoff
    for (std::set<QM_Atom>::const_iterator
        qm_it = this->qm.begin(), qm_to = this->qm.end(); qm_it != qm_to; ++qm_it)
      {
      const math::Vec& qm_pos = qm_it->pos;
      math::Vec r_qm_mm;
      std::set<MM_Atom>::iterator mm_it = this->mm.begin();
      for(unsigned i = 0; i < num_atoms; ++i) {
        if (topo.is_qm(i)) continue;
        // get distance
        periodicity.nearest_image(qm_pos, pos(i), r_qm_mm);
        DEBUG(15, "Nearest image of atom " << i << " : " << math::v2s(r_qm_mm));
        DEBUG(15, "math::abs2(r_qm_mm) = " << math::abs2(r_qm_mm));
        if (math::abs2(r_qm_mm) < cutoff2) {
          DEBUG(9, "Adding MM atom " << i);
          // MM atom position
          const math::Vec mm_pos = qm_pos - r_qm_mm;
          /** this will add every MM atom only once. We have to check, if we already
           * have the atom and if its position is the same. Otherwise we see its
           * periodic copy and that is bad.
           */
          const size_t size = this->mm.size();
          if (mm_it != this->mm.end())
            std::advance(mm_it,1); // hint on set insertion, improves performance
          mm_it = this->mm.emplace_hint(mm_it, i, topo.charge(i), mm_pos);
          if (size == this->mm.size()) {
            DEBUG(9, "Atom " << i << " already in the list");
            const double delta = math::abs2(mm_it->pos - mm_pos);
            if (delta > math::epsilon) {
              // Atom already in the list but as different periodic image, exit
              std::ostringstream msg;
              msg << "QM atom " << (qm_it->index + 1)
                  << " sees multiple periodic copies of MM atom "
                  << (mm_it->index + 1);
              io::messages.add(msg.str(), "QM Zone", io::message::error);
              return 1;
            } // if not in the same position
          } // if already in the set
        } // if atom within cutoff
      } // for all MM atoms
    } // for all QM atoms
  } // if electrostatic embedding
  return 0;
}

template<math::boundary_enum B>
void interaction::QM_Zone::get_linked_mm_atoms(const topology::Topology& topo
                                             , const configuration::Configuration& conf
                                             , const math::Periodicity<B>& periodicity)
{
  for (std::set< std::pair<unsigned,unsigned> >::const_iterator
        it = topo.qmmm_link().begin(), to = topo.qmmm_link().end();
        it != to; ++it) {
    const unsigned qm_i = it->first,
                    mm_i = it->second;
    DEBUG(15, "Gathering MM atom " << mm_i << " linked to QM atom " << qm_i);
    assert(topo.is_qm(qm_i));
    assert(!topo.is_qm(mm_i));
    std::set<QM_Atom>::const_iterator qm_it = this->qm.find(QM_Atom(qm_i));
    assert(qm_it != this->qm.end());
    const math::Vec& qm_pos = qm_it->pos;
    math::Vec r_qm_mm;
    periodicity.nearest_image(qm_pos, conf.current().pos(mm_i), r_qm_mm);

    DEBUG(15, "Linked QM position: " << math::v2s(qm_pos));
    DEBUG(15, "Nearest image vector: " << math::v2s(r_qm_mm));
    const math::Vec mm_pos = qm_pos - r_qm_mm;

    DEBUG(15, "Linked MM position: " << math::v2s(mm_pos));
    this->mm.emplace(mm_i, topo.charge(mm_i), mm_pos);
    DEBUG(15, "Adding MM atom: " << mm_i);
  }
}

void interaction::QM_Zone::scale_charges(const simulation::Simulation& sim) {
  for (std::set<MM_Atom>::iterator
        mm_it = this->mm.begin(), mm_to = this->mm.end();
        mm_it != mm_to; ++mm_it)
    {
    if (mm_it->charge != 0.0) {
      // Find closest QM atom and scale accordingly
      double min_d2;
      for (std::set<QM_Atom>::const_iterator
            qm_it = this->qm.begin(), qm_to = this->qm.end();
            qm_it != qm_to; ++qm_it)
        {
          min_d2 = std::min(math::abs2(qm_it->pos - mm_it->pos), min_d2);
      }
      const double min_d = sqrt(min_d2);
      mm_it->charge = 0.6366197723675814/* = (2/pi) */
                    * mm_it->charge * (atan(sim.param().qmmm.mm_scale * min_d));
    }
  }
}

void interaction::QM_Zone::get_links(const topology::Topology& topo,
                                     const simulation::Simulation& sim)
  {
  this->link.clear();
  typedef std::set<std::pair<unsigned,unsigned>> linkset;
  const linkset& links = topo.qmmm_link();
  for (linkset::const_iterator it = links.begin(), to = links.end(); it != to; ++it)
    {
    // assert that atoms are gathered in the sets
    assert(this->qm.find(QM_Atom(it->first)) != this->qm.end());
    assert(this->mm.find(MM_Atom(it->second)) != this->mm.end());

    DEBUG(9,"Adding QM-MM link: " << it->first << " - " << it->second);
    this->link.emplace(it->first, it->second, 0);
  }
  this->update_links(sim);
}

void interaction::QM_Zone::update_links(const simulation::Simulation& sim) {
  for (std::set<QM_Link>::iterator
    it = this->link.begin(), to = this->link.end(); it != to; ++it)
    {
    DEBUG(9, "Updating position of capping atom: link " << it->qm_index << " - " << it->mm_index);
    std::set<QM_Atom>::const_iterator
      qm_it = this->qm.find(QM_Atom(it->qm_index));
    assert(qm_it != this->qm.end());
    std::set<MM_Atom>::const_iterator
      mm_it = this->mm.find(MM_Atom(it->mm_index));
    assert(mm_it != this->mm.end());
  
    const math::Vec& qm_pos = qm_it->pos;
    const math::Vec& mm_pos = mm_it->pos;

    it->update_cap_position(qm_pos, mm_pos, sim.param().qmmm.cap_length);
  }
}

void interaction::QM_Zone::update_cos(const topology::Topology& topo
                                    , const configuration::Configuration& conf
                                    , const simulation::Simulation& sim) {
  for (std::set<MM_Atom>::iterator
      it = this->mm.begin(), to = this->mm.end(); it != to; ++it)
    {
    const unsigned i = it->index;
    if (topo.is_polarisable(i)) {
      DEBUG(9, "Updating charge-on-spring of MM atom: " << i);
      it->is_polarisable = true;
      it->cos_charge = topo.coscharge(i);
      it->cosV = conf.current().posV(i);
    }
  }
}

void interaction::QM_Zone::electric_field(const simulation::Simulation& sim
                                        , math::VArray& electric_field) {
  // get electric field at either the charge or the COS site.
  for (std::set<MM_Atom>::iterator
      it = this->mm.begin(), to = this->mm.end(); it != to; ++it)
    {
    const unsigned i = it->index;
    math::Vec e;
    switch (sim.param().polarise.efield_site) {
      case simulation::ef_atom:
        e = it->force / (it->charge - it->cos_charge);
        break;
      case simulation::ef_cos:
        e = it->cos_force / it->cos_charge;
        break;
      default:
        io::messages.add("Electric field site not implemented.",
                "QMMM_Interaction", io::message::critical);
    }
    electric_field(i) += e;
  }
}