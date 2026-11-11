/**
 * @file colvar.cc
 * @brief Implementation of base Colvar class
 */

#include "../../stdheader.h"
#include "colvar.h"
#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../configuration/configuration.h"
#include "../../../math/periodicity.h"
#include "../../../util/template_split.h"
#include "../../../util/virtual_atom.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

namespace {

template<math::boundary_enum B>
void apply_colvar_forces(topology::Topology &topo,
                         configuration::Configuration &conf,
                         const std::vector<util::Virtual_Atom*> &atoms,
                         const math::VArray &derivatives,
                         double dE_dinstant, bool include_virial,
                         math::VArray &forces, math::Matrix &virial,
                         bool molecular_virial)
{
  math::Periodicity<B> periodicity(conf.current().box);
  // Forcefield's atomic-to-molecular conversion precedes special interactions.
  // Convert only this CV contribution when it is added afterwards.
  std::map<int,math::Vec> before;
  if (include_virial && molecular_virial && B!=math::vacuum)
    for (const auto *site : atoms) for (int i=0;i<site->size();++i)
      before.emplace(site->atom(i),forces(site->atom(i)));
  for (size_t endpoint = 0; endpoint < atoms.size(); ++endpoint) {
    const util::Virtual_Atom &atom = *atoms[endpoint];
    const math::Vec force = -dE_dinstant * derivatives[endpoint];
    DEBUG(15, "COLVARRES endpoint " << endpoint + 1
          << " derivative " << math::v2s(derivatives[endpoint])
          << " force " << math::v2s(force));
    // Explicit atoms and linear centers already deform affinely.
    if (!include_virial || atom.type() == util::va_explicit
        || atom.type() == util::va_cog || atom.type() == util::va_com) {
      atom.force(conf, topo, force, forces);
      continue;
    }

    // Snapshot only this site's real atoms. Subtracting their previous forces
    // isolates this endpoint, including when different endpoints share atoms.
    // A map also avoids counting repeated defining atoms more than once.
    std::map<int, math::Vec> previous_forces;
    for (int i = 0; i < atom.size(); ++i)
      previous_forces[atom.atom(i)] = forces(atom.atom(i));

    const math::Vec anchor = conf.current().pos(atom.atom(0));
    const math::Vec site_relative = atom.pos(conf, topo) - anchor;
    atom.force(conf, topo, force, forces);

    // The geometry tensor treats virtual endpoints as affine points. Replace
    // that assumption by the actual redistributed atomic forces:
    // correction = sum_i (r_i - r_site) outer F_i.
    // Gather defining atoms relative to the same anchor used by Virtual_Atom.
    for (const auto &entry : previous_forces) {
      math::Vec relative;
      periodicity.nearest_image(conf.current().pos(entry.first), anchor, relative);
      relative -= site_relative;
      const math::Vec atomic_force = forces(entry.first) - entry.second;
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          virial(a, b) += relative(a) * atomic_force(b);
        }
      }
    }
  }
  if (!before.empty()) {
    const auto &groups=topo.pressure_groups();
    const auto &positions=conf.current().pos;
    for (size_t group=0;group+1<groups.size();++group) {
      auto affected=before.lower_bound(groups[group]);
      if (affected==before.end()) break;
      if (static_cast<unsigned int>(affected->first)>=groups[group+1]) continue;
      math::Vec previous=positions(groups[group]),center(0);
      double mass=0;
      for (unsigned int atom=groups[group];atom<groups[group+1];++atom) {
        math::Vec delta;
        periodicity.nearest_image(positions(atom),previous,delta);
        previous+=delta;
        center+=topo.mass(atom)*previous;
        mass+=topo.mass(atom);
      }
      center/=mass;
      // Match util::atomic_to_molecular_virial's pressure-group convention.
      for (;affected!=before.end() && static_cast<unsigned int>(affected->first)<groups[group+1];++affected) {
        math::Vec relative;
        periodicity.nearest_image(positions(affected->first),center,relative);
        const auto force=forces(affected->first)-affected->second;
        for (int a=0;a<3;++a) for (int b=0;b<3;++b)
          virial(a,b)-=relative(a)*force(b);
      }
    }
  }
}

} // anonymous namespace

namespace interaction {

void Colvar::apply_forces(topology::Topology &topo, configuration::Configuration &conf,
                          double dE_dinstant, bool include_virial,bool molecular_virial) const
{
  apply_forces(topo, conf, dE_dinstant, include_virial,
               conf.current().force, conf.current().virial_tensor,molecular_virial);
}

void Colvar::apply_forces(topology::Topology &topo, configuration::Configuration &conf,
                          double dE_dinstant, bool include_virial,
                          math::VArray &forces, math::Matrix &virial,bool molecular_virial) const
{
  SPLIT_BOUNDARY(apply_colvar_forces, topo, conf, atoms, derivatives,
                 dE_dinstant, include_virial, forces, virial,molecular_virial);
  if (include_virial) add_virial(virial, dE_dinstant);
}

void Colvar::add_virial(math::Matrix &virial, double dE_dinstant) const
{
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      virial(a, b) -= dE_dinstant * strain_derivative(a, b);
    }
  }
}

} // namespace interaction
