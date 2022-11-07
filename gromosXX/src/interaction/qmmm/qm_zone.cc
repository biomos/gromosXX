/**
 * @file qm_zone.cc
 * Implements methods of QM_Zone
 */
#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/interaction.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"

#include "../../../math/periodicity.h"
#include "../../../util/template_split.h"

#include "../../../util/debug.h"
#include "../../../util/error.h"

#include "qm_atom.h"
#include "mm_atom.h"
#include "qm_link.h"
#include "qm_zone.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::QM_Zone::QM_Zone(int net_charge
                            , int spin_mult)
                              : m_qm_energy(0.0)
                              , m_charge(net_charge)
                              , m_spin_mult(spin_mult)
{}

interaction::QM_Zone::~QM_Zone() = default;

void interaction::QM_Zone::zero() {
  this->m_qm_energy = 0.0; 
}

void interaction::QM_Zone::clear() {
  this->zero();
  this->qm.clear();
  this->mm.clear();
  this->link.clear();
}

int interaction::QM_Zone::init(topology::Topology& topo, 
                               const configuration::Configuration& conf, 
                               const simulation::Simulation& sim) {
  int err = 0;
  
  DEBUG(15,"Getting QM atoms");
  if ((err = this->get_qm_atoms(topo, conf, sim)))
    return err;

  if (sim.param().qmmm.buffer_zone.cutoff) {
    DEBUG(15,"Getting buffer atoms");
    this->get_buffer_atoms(topo, conf, sim);
  }

  DEBUG(15,"Getting MM atoms");
  this->get_mm_atoms(topo, conf, sim);

  DEBUG(15,"Getting QM-MM links");
  this->get_links(topo, sim);
  return 0;

}

int interaction::QM_Zone::update(topology::Topology& topo, 
                                 const configuration::Configuration& conf, 
                                 const simulation::Simulation& sim) {
  this->zero();

  DEBUG(15,"Updating QM atoms positions");
  update_qm_pos(topo, conf, sim);

  if (sim.param().qmmm.buffer_zone.cutoff) {
    DEBUG(15,"Getting buffer atoms");
    this->get_buffer_atoms(topo, conf, sim);
  }
  
  DEBUG(15,"Getting MM atoms positions");
  this->get_mm_atoms(topo, conf, sim);
  
  DEBUG(15,"Updating links");
  this->update_links(sim);
  return 0;
}

void interaction::QM_Zone::write(topology::Topology& topo, 
                                 configuration::Configuration& conf, 
                                 const simulation::Simulation& sim) {
  // Write positions
  // If geometry minimisation within QM is requested
  /*if (minimise) {
  for (std::set<QM_Atom>::const_iterator
      it = this->qm.begin(), to = this->qm.end(); it != to; ++it)
      {
      conf.current().pos(it->index) = it->pos;
    }
  }*/
  
  // Write forces
  // First distribute capping atom forces
  for (std::set<QM_Link>::const_iterator
      it = this->link.begin(), to = this->link.end(); it !=to; ++it)
    {
    DEBUG(15, "Redistributing force of capping atom");
    DEBUG(15, "QM-MM link: " << it->qm_index << " - " << it->mm_index);
    std::set<QM_Atom>::iterator qm_it = this->qm.find(it->qm_index);
    std::set<MM_Atom>::iterator mm_it = this->mm.find(it->mm_index);
    assert(qm_it != this->qm.end());
    assert(mm_it != this->mm.end());
    it->distribute_force(qm_it->pos, mm_it->pos, qm_it->force, mm_it->force);
  }

  // Then write QM atoms forces and calculate virial
  math::Matrix virial_tensor(0.0);
  for (std::set<QM_Atom>::const_iterator
      it = this->qm.begin(), to = this->qm.end(); it != to; ++it)
    {
    DEBUG(15, "Atom " << it->index << ", force: " << math::v2s(it->force));
    conf.current().force(it->index) += it->force;
    math::Vec& pos = conf.current().pos(it->index);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        virial_tensor(j,i) += pos(j) * it->force(i);
      }
    }
  }

  // And MM atoms forces
  for (std::set<MM_Atom>::const_iterator
      it = this->mm.begin(), to = this->mm.end(); it != to; ++it)
    {
    DEBUG(15, "Atom " << it->index << ", force: " << math::v2s(it->force));
    conf.current().force(it->index) += it->force;
    math::Vec& pos = conf.current().pos(it->index);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        virial_tensor(j,i) += pos(j) * it->force(i);
      }
    }
  }
  
  if (sim.param().pcouple.virial) {
    conf.current().virial_tensor += virial_tensor;
  }
  
  // Write charges
  if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical
      && sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
    for (std::set<QM_Atom>::const_iterator
          it = this->qm.begin(), to = this->qm.end(); it != to; ++it)
      {
      if (topo.is_qm(it->index)) {
        DEBUG(15, "Atom " << it->index << ", new charge: " << it->qm_charge);
        topo.charge()(it->index) = it->qm_charge;
      }
    }
    // Add capping atom charge to QM link atom
    for (std::set<QM_Link>::const_iterator
          it = this->link.begin(), to = this->link.end(); it != to; ++it)
      {
      DEBUG(15, "Capping atom " << it->qm_index << "-" << it->mm_index << ", charge: " << it->qm_charge);
      topo.charge()(it->qm_index) += it->qm_charge;
      DEBUG(15, "Charge added to QM atom " << it->qm_index << ", new charge: " << topo.charge(it->qm_index));
    }
  }

  // write separate charges for interaction between the buffer region and MM region
  if (sim.param().qmmm.use_qm_buffer
      && sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic
      && (sim.param().qmmm.software != simulation::qm_nn
          || sim.steps() % sim.param().qmmm.nn.charge_steps == 0))
    {
    // First reset delta charges
    std::fill(topo.qm_delta_charge().begin(), topo.qm_delta_charge().end(), 0.0);
    // Fill with delta-charges
    for (std::set<QM_Atom>::const_iterator
        it = this->qm.begin(), to = this->qm.end(); it != to; ++it) {
      if (topo.is_adaptive_qm_buffer(it->index)) {
        const double delta_charge = it->qm_charge - topo.charge(it->index);
        DEBUG(15, "Atom " << it->index << ", new delta charge: " << delta_charge);
        topo.qm_delta_charge(it->index) = delta_charge;
      }
    }
  }

  // Write energies
  conf.current().energies.qm_total = this->m_qm_energy;
}

int interaction::QM_Zone::get_qm_atoms(const topology::Topology& topo, 
                                       const configuration::Configuration& conf, 
                                       const simulation::Simulation& sim) {
  this->qm.clear();
  DEBUG(15,"QM_Zone::get_qm_atoms: Collecting atoms");
  const bool adaptive_buffer = sim.param().qmmm.buffer_zone.cutoff;
  // Check if one CG is not split between QM and MM
  // Also collect QM chargegroups
  for (unsigned cg = 0; cg < topo.num_chargegroups(); ++cg) {
    unsigned a = topo.chargegroup(cg);
    const unsigned a_to = topo.chargegroup(cg + 1);
    DEBUG(15, "cg: " << cg);
    DEBUG(15, "a_to: " << a_to);
    assert(a_to <= topo.num_atoms());
    /** QM_Zone will also contain QM_Buffer atoms, since the full QM calculation will be run
     * on this system, but we are not excluding the QM_Buffer atoms in MM pairlist generation.
     */
    const bool cg_is_qm = topo.is_qm(a);
    const bool cg_is_buf = topo.is_qm_buffer(a);
    assert(!(cg_is_qm && cg_is_buf));
    for (; a < a_to; ++a) {
      const bool a_is_qm = topo.is_qm(a);
      const bool a_is_buf = topo.is_qm_buffer(a);
      assert(!(a_is_qm && a_is_buf));
      if ((cg_is_qm != a_is_qm)
          || (cg_is_buf != a_is_buf)
        ) {
        std::ostringstream msg;
        msg << "Chargegroup " << cg << " is split between QM and MM zone - "
            << "atoms " << a + 1 << " and " << a;
        io::messages.add(msg.str(), "QM_Zone", io::message::error);
        return E_INPUT_ERROR;
      }
      if (a_is_qm || a_is_buf
        ) {
        if (topo.is_polarisable(a)
                || topo.is_coarse_grained(a)
                || topo.is_perturbed(a)
                || topo.is_eds_perturbed(a)) {
          std::ostringstream msg;
          msg << "QM or buffer atom should not be polarisable, coarse-grained or perturbed";
          io::messages.add(msg.str(), "QM_Zone", io::message::error);
          return E_INPUT_ERROR;
        }
        // Skip adaptive buffer zone
        if (adaptive_buffer && a_is_buf) continue;

        DEBUG(7,"Adding QM atom " << a << ", atomic number: " << topo.qm_atomic_number(a));
        this->qm.emplace(a, math::Vec(0.0), topo.qm_atomic_number(a));
      }
    }
  }
  if (this->qm.empty()) {
    io::messages.add("QMMM requested but no QM atoms found", "QM_Zone", io::message::error);
    return E_INPUT_ERROR;
  }
  DEBUG(15,"QM_Zone::get_qm_atoms: Updating QM atoms positions");
  this->update_qm_pos(topo, conf, sim);
  return 0;
}

void interaction::QM_Zone::update_qm_pos(const topology::Topology& topo, 
                                         const configuration::Configuration& conf, 
                                         const simulation::Simulation& sim) {
  DEBUG(15,"QM_Zone::update_qm_pos: Splitting boundary");
  SPLIT_BOUNDARY(this->_update_qm_pos, topo, conf, sim);
}

template<math::boundary_enum B>
int interaction::QM_Zone::_update_qm_pos(const topology::Topology& topo, 
                                         const configuration::Configuration& conf, 
                                         const simulation::Simulation& sim) {
  int err = 0;
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
          io::messages.add(msg.str(), "QM_Zone", io::message::error);
          err = 1;
        }
      }
    }
  }
  return err;
}

void interaction::QM_Zone::update_pairlist(
                  const topology::Topology& topo
                , const simulation::Simulation& sim
                , PairlistContainer& pairlist
                , unsigned begin, unsigned end
                , unsigned stride) const {
  /** This pairlist generator is not adapted for striding.
   *  To implement parallelism, uncomment the if statement.
   *  Since we never use more than 1 set, this is disabled
   *  due to the performance reasons.
   *  So we also assert...
  **/
  assert(begin == 0 && stride == 1);
  pairlist.clear();
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    // Generate QM-MM pairs
    //unsigned counter = 0; // Not implemented
    for (std::set<QM_Atom>::const_iterator
          qm_it = this->qm.begin(), qm_to = this->qm.end();
          qm_it != qm_to; ++qm_it) {
      for (std::set<MM_Atom>::const_iterator
            mm_it = this->mm.begin(), mm_to = this->mm.end();
            mm_it != mm_to; ++mm_it) {
        //if ((counter++ - begin) % stride == 0) { // Not implemented
          unsigned i = qm_it->index;
          unsigned j = mm_it->index;
          if (!topo.is_qm_buffer(i)) {
            if (i > j) std::swap(i,j);
            // There should not be any QM pair in standard exclusions
            assert(!topo.all_exclusion(i).is_excluded(j));
            if (!topo.qm_all_exclusion(i).is_excluded(j)) {
              pairlist.solute_short[i].push_back(j);
            }
          }
        //}
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

template <math::boundary_enum B, class AtomType>
int interaction::QM_Zone::gather_chargegroups(const topology::Topology& topo, 
                                              const configuration::Configuration& conf, 
                                              const simulation::Simulation& sim,
                                              std::set<AtomType>& atom_set,
                                              const double cutoff2) {
  math::Periodicity<B> periodicity(conf.current().box);
  /** We are running twice almost the same loop for solute and solvent.
   * COG of solute is calculated normally, but in solvent, first atom is
   * considered the COG
   * Lambdas cannot be used because passing capturing lambda to lambda
   * is not allowed. We can solve this by passing local functors to lambda.
   */
  class cog_calculator {
    public:
    virtual ~cog_calculator() {};
      virtual const math::Vec& operator()(unsigned i) const = 0;
    };

  class solute_cog_calculator : public cog_calculator {
      const math::VArray& m_cogs;
    public:
    solute_cog_calculator(const math::VArray& cogs)
                : m_cogs(cogs) {};
      const math::Vec& operator()(unsigned i) const {
        return m_cogs(i);
      }
    };

  class solvent_cog_calculator : public cog_calculator {
      const topology::Topology& m_topo;
      const math::VArray& m_pos;
    public:
    solvent_cog_calculator(const topology::Topology& topo, const math::VArray& pos)
                :m_topo(topo), m_pos(pos) {};
      const math::Vec& operator()(unsigned i) const {
        return m_pos(m_topo.chargegroup(i));
      }
    };
  // End of functors definition

  // Initialize constants
  const math::VArray& pos = conf.current().pos;
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
  cog_calculator* solute_cgs = new solute_cog_calculator(cogs);
  // Functor for solvent CG COGs
  cog_calculator* solvent_cgs = new solvent_cog_calculator(topo, pos);

  // Collect MM chargegroups within every QM atom cutoff
  for (std::set<QM_Atom>::const_iterator
      qm_it = this->qm.begin(), qm_to = this->qm.end(); qm_it != qm_to; ++qm_it)
    {
    // Initial dumb implementation - gather only around the first atom, if
    // AtomType == QM_Atom (true only for QM buffer atoms gathering)
    if (std::is_same<AtomType,QM_Atom>::value && qm_it != this->qm.begin()) {
      break;
    }
    DEBUG(15, "Gathering around QM atom " << qm_it->index);
    const math::Vec& qm_pos = qm_it->pos;
    math::Vec r_qm_cg;
    typename std::set<AtomType>::iterator set_it = atom_set.begin();
    // We run almost the same loop twice
    // With functors and lambda, one definition is enough
    auto cg_cogs_loop = [&](unsigned cg_from, unsigned cg_to, cog_calculator* cog)-> int {
      for(unsigned cg = cg_from; cg < cg_to; ++cg) {
        unsigned a = topo.chargegroup(cg);
        // The whole chargegroup is the same type, so skip based on the type
        if (this->skip_cg<AtomType>(topo, a)) continue;
        // get distance
        periodicity.nearest_image(qm_pos, (*cog)(cg), r_qm_cg);
        DEBUG(15, "qm_pos: " << math::v2s(qm_pos));
        DEBUG(15, "(*cog)(cg): " << math::v2s((*cog)(cg)));
        DEBUG(15, "Nearest image qm_pos - cg " << cg << " : " << math::v2s(r_qm_cg));
        DEBUG(15, "math::abs2(r_qm_cg) = " << math::abs2(r_qm_cg));
        
        if (math::abs2(r_qm_cg) < cutoff2) {
          // Iterate over atoms in cg
          for (unsigned a_to = topo.chargegroup(cg + 1); a < a_to; ++a) {
            // nearest image of atom to CG COG
            math::Vec r_mm_cg;
            periodicity.nearest_image(pos(a), (*cog)(cg), r_mm_cg);
            DEBUG(15, "nearest image mm_pos-cg_cog = " << math::v2s(r_mm_cg));
            DEBUG(15, "math::abs2(r_mm_cg) = " << math::abs2(r_mm_cg));
            // MM atom position
            const math::Vec mm_pos = qm_pos - r_qm_cg + r_mm_cg;
            /** this will add every MM atom only once. We have to check, if we already
             * have the atom and if its position is the same. Otherwise we see its
             * periodic copy and that is bad.
             */
            DEBUG(15, "mm_pos: " << math::v2s(mm_pos));
            const size_t size = atom_set.size();
            if (set_it != atom_set.end())
              std::advance(set_it,1); // hint on set insertion, improves performance
            this->emplace_atom<AtomType>(atom_set, set_it, a, mm_pos,
                                          topo.qm_atomic_number(a), topo.charge(a),
                                          topo.is_qm_buffer(a));

            if (size == atom_set.size()) {
              DEBUG(15, "Atom " << a << " already in the list");
              const double delta = math::abs2(set_it->pos - mm_pos);
              DEBUG(15, "set_it->pos: " << math::v2s(set_it->pos));
              DEBUG(15, "mm_pos:     " << math::v2s(mm_pos));
              DEBUG(15, "delta: " << delta);
              if (delta > math::epsilon) {
                // Atom already in the list but as different periodic image, exit
                std::ostringstream msg;
                msg << "QM atom " << (qm_it->index + 1)
                    << " sees different periodic copy of MM atom "
                    << (set_it->index + 1);
                io::messages.add(msg.str(), "QM_Zone", io::message::error);
                return E_BOUNDARY_ERROR;
              } // if not in the same position
            } // if already in the set
          } // for atoms in cg
        } // if cg within cutoff
      } // for cgs in the range
      return 0;
    }; // lambda cg_cogs_loop
    int err = 0;
    DEBUG(15, "Iterating over solute cgs");
    err = cg_cogs_loop(0, num_solute_cg, solute_cgs);
    if (err) return err;

    DEBUG(15, "Iterating over solvent cgs");
    err = cg_cogs_loop(num_solute_cg, num_cg, solvent_cgs);
    if (err) return err;
  }
  delete solute_cgs;
  delete solvent_cgs;
  return 0;
}

// emplace is very similar for QM and MM atom, except the buffer logic is flipped, so we use template here
template<>
void interaction::QM_Zone::emplace_atom(std::set<interaction::QM_Atom>& set,
                                        std::set<interaction::QM_Atom>::const_iterator& it,
                                        const unsigned index,
                                        const math::Vec& pos,
                                        const unsigned atomic_number,
                                        const double charge,
                                        const bool is_qm_buffer)
  {
  if (is_qm_buffer) {
    DEBUG(15, "Adding QM buffer atom " << index);
    it = set.emplace_hint(it, index, pos, atomic_number);
  }
}

template<>
void interaction::QM_Zone::emplace_atom(std::set<interaction::MM_Atom>& set,
                                        std::set<interaction::MM_Atom>::const_iterator& it,
                                        const unsigned index,
                                        const math::Vec& pos,
                                        const unsigned atomic_number,
                                        const double charge,
                                        const bool is_qm_buffer) 
{
  if (!is_qm_buffer) {
    DEBUG(15, "Adding MM atom " << index);
    it = set.emplace_hint(it, index, pos, atomic_number, charge);
  }
}

template<>
bool interaction::QM_Zone::skip_cg<interaction::QM_Atom>
                                (const topology::Topology& topo,
                                 const unsigned index)
  {
  return (topo.is_qm(index) || !topo.is_qm_buffer(index));
}

template<>
bool interaction::QM_Zone::skip_cg<interaction::MM_Atom>
                                (const topology::Topology& topo,
                                 const unsigned index)
  {
  return (topo.is_qm(index) || topo.is_qm_buffer(index));
}

void interaction::QM_Zone::get_buffer_atoms(topology::Topology& topo, 
                                            const configuration::Configuration& conf, 
                                            const simulation::Simulation& sim) {
  SPLIT_BOUNDARY(this->_get_buffer_atoms, topo, conf, sim);
}

template<math::boundary_enum B>
int interaction::QM_Zone::_get_buffer_atoms(topology::Topology& topo, 
                                            const configuration::Configuration& conf, 
                                            const simulation::Simulation& sim) {
  /** Here we need to iterate over all buffer atoms and see, if they are within the
   * adaptive buffer cutoff
   */
  int err = 0;
  // Firstly remove buffer atoms from the QM set
  DEBUG(15, "Firstly removing the old buffer atoms");
  for (std::set<QM_Atom>::const_iterator qm_it = this->qm.begin();
        qm_it != this->qm.end();) {
    const unsigned i = qm_it->index;
    if (topo.is_qm_buffer(i)) {
      DEBUG(12, "Removing QM buffer atom " << qm_it->index);
      qm_it = this->qm.erase(qm_it);
    }
    else {
      ++qm_it;
    }
  }
  const double cutoff2 = sim.param().qmmm.buffer_zone.cutoff * sim.param().qmmm.buffer_zone.cutoff;
  std::set<QM_Atom> buffer_atoms;
  if ((err = this->gather_chargegroups<B>(topo, conf, sim, buffer_atoms, cutoff2)))
    return err;
  
  // update the QM buffer topology so the pairlist algorithm can skip them
  for (unsigned i = 0; i < topo.num_atoms(); ++i) {
    if (topo.is_qm_buffer(i)) {
      if (buffer_atoms.count(i)) {
        topo.is_qm_buffer(i) = 1;
        DEBUG(9, "Atom " << i << " in adaptive buffer");
      } else {
        topo.is_qm_buffer(i) = -1; // temporarily disabled buffer atom
        DEBUG(15, "Atom " << i << " not in adaptive buffer");
      }
    }
  }

  DEBUG(15, "Buffer atoms:");
  for (std::set<QM_Atom>::const_iterator it = buffer_atoms.begin(); it != buffer_atoms.end(); ++it) {
    DEBUG(15, it->index);
    DEBUG(15, "Distance to the 1st atom: " << math::abs(it->pos - this->qm.begin()->pos));
  }

  // And merge with QM
  this->qm.insert(buffer_atoms.begin(), buffer_atoms.end());

  // Also update their positions
  if ((err = this->_update_qm_pos<B>(topo, conf, sim)))
    return err;
  return 0;
}

void interaction::QM_Zone::get_mm_atoms(const topology::Topology& topo, 
                                        const configuration::Configuration& conf, 
                                        const simulation::Simulation& sim) {
  /** The default pairlist algorithms cannot be applied here since we are
   * using different cutoff. We also do gathering at the same time and perform
   * some checks on cutoff and periodic copy interactions so it becomes
   * rather specific for QM
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
      this->mm.emplace(i, conf.current().pos(i), topo.qm_atomic_number(i), topo.charge(i));
    }
  }
  else if (sim.param().boundary.boundary == math::vacuum
            && sim.param().qmmm.qmmm != simulation::qmmm_mechanical
            && sim.param().qmmm.cutoff == 0.0) {
    // Include all - allowed only with vacuum PBC
    DEBUG(9, "Gathering all MM atoms");
    for (unsigned int i = 0; i < topo.num_atoms(); ++i) {
      if ( !topo.is_qm(i) && !topo.is_qm_buffer(i) ) {
        DEBUG(9, "Adding atom " << i);
        this->mm.emplace(i, conf.current().pos(i), topo.qm_atomic_number(i), topo.charge(i));
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
}

template<math::boundary_enum B>
int interaction::QM_Zone::_get_mm_atoms(const topology::Topology& topo, 
                                        const configuration::Configuration& conf, 
                                        const simulation::Simulation& sim)
  {
  int err = 0;
  math::Periodicity<B> periodicity(conf.current().box);
  DEBUG(15, "Boundary type: " << conf.boundary_type);
  // MM link atoms need to be always gathered
  this->get_linked_mm_atoms<B>(topo, conf, periodicity);
  // Gather MM atoms only for electrostatic or polarisable embedding
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    DEBUG(9, "Gathering MM atoms using chargegroup-based cutoff");
    const double cutoff2 = sim.param().qmmm.cutoff * sim.param().qmmm.cutoff;
    if ((err = this->gather_chargegroups<B>(topo, conf, sim, this->mm, cutoff2)))
      return err;
    
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
        if (topo.is_qm(i) || topo.is_qm_buffer(i)) continue;
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
          mm_it = this->mm.emplace_hint(mm_it, i, mm_pos, topo.qm_atomic_number(i), topo.charge(i));
          if (size == this->mm.size()) {
            DEBUG(9, "Atom " << i << " already in the list");
            const double delta = math::abs2(mm_it->pos - mm_pos);
            if (delta > math::epsilon) {
              // Atom already in the list but as different periodic image, exit
              std::ostringstream msg;
              msg << "QM atom " << (qm_it->index + 1)
                  << " sees multiple periodic copies of MM atom "
                  << (mm_it->index + 1);
              io::messages.add(msg.str(), "QM_Zone", io::message::error);
              return E_BOUNDARY_ERROR;
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
    assert(topo.is_qm(qm_i) || topo.is_qm_buffer(qm_i));
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
    this->mm.emplace(mm_i, mm_pos, topo.qm_atomic_number(mm_i), topo.charge(mm_i));
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
      // Initialize to large value
      double min_d2 = std::numeric_limits<double>::max();
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
    assert(this->qm.find(it->first) != this->qm.end());
    assert(this->mm.find(it->second) != this->mm.end());

    DEBUG(9,"Adding QM-MM link: " << it->first << " - " << it->second);
    this->link.emplace(QM_Atom(0), it->first, it->second);
    this->qm.find(it->first)->is_linked = true;
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
        io::messages.add("Electric field site not implemented",
                "QMMM_Interaction", io::message::critical);
    }
    electric_field(i) += e;
  }
}

interaction::QM_Zone* interaction::QM_Zone::create_buffer_zone(
                                          const topology::Topology& topo
                                        , const simulation::Simulation& sim
                                          ) {
  DEBUG(15,"Generating QM buffer zone");
  // We make the buffer zone by copying the QM zone and deleting the QM atoms
  QM_Zone* qm_buffer = new interaction::QM_Zone(*this);
  qm_buffer->zero();
  qm_buffer->charge() = sim.param().qmmm.buffer_zone.charge;
  qm_buffer->spin_mult() = sim.param().qmmm.buffer_zone.spin_mult;

  for (std::set<interaction::QM_Atom>::const_iterator
        it = this->qm.begin(); it != this->qm.end(); ++it) {
    if (topo.is_qm(it->index)) {
      DEBUG(10,"Erasing QM atom " << it->index << " from buffer zone");
      // 'it' remains valid, since we are iterating over 'this->qm' but erasing in 'qm_buffer->qm'
      qm_buffer->qm.erase(*it);
    }
  }
  return qm_buffer;
}
