/** Adapt LE kernel derivatives to the shared CV geometry/force layer. */
#include "../stdheader.h"
#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../configuration/configuration.h"
#include "../simulation/simulation.h"
#include "../interaction/special/colvar/colvar.h"
#include "../interaction/special/colvar/distance.h"
#include "../interaction/special/colvar/dihedral.h"
#include "../interaction/special/colvar/angle.h"
#include "../interaction/special/colvar/distance_difference.h"
#include "../interaction/special/colvar/coordnum.h"
#include "le_colvar_coordinate.h"

namespace {
std::vector<util::Virtual_Atom> explicit_sites(const std::vector<unsigned int> &atoms) {
  std::vector<util::Virtual_Atom> result;
  for (auto atom : atoms) result.emplace_back(util::va_explicit,std::vector<int>(1,atom));
  return result;
}
}

struct util::LE_Colvar_Coordinate::Impl {
  int type;
  std::vector<unsigned int> indices;
  std::vector<util::Virtual_Atom> sites;
  bool site_definition;
  std::unique_ptr<topology::distance_restraint_struct> distance;
  std::unique_ptr<topology::dihedral_restraint_struct> dihedral;
  std::unique_ptr<topology::angle_restraint_struct> angle;
  std::unique_ptr<topology::distance_difference_restraint_struct> difference;
  std::unique_ptr<topology::coordnum_restraint_struct> coordnum;
  size_t group1_size;
  double rcut;
  int nn, mm;
  std::unique_ptr<interaction::Colvar> cv;
  topology::Topology *topo;
  configuration::Configuration *conf;
  bool initialised, valid, molecular_virial;
  Impl(int type, const std::vector<unsigned int> &atoms, size_t group1_size=0,
       double rcut=0, int nn=0, int mm=0)
      : Impl(type,explicit_sites(atoms),group1_size,rcut,nn,mm,false) {}
  Impl(int type,const std::vector<util::Virtual_Atom> &sites,size_t group1_size,
       double rcut,int nn,int mm,bool site_definition=true)
      : type(type), sites(sites), site_definition(site_definition),
        group1_size(group1_size), rcut(rcut), nn(nn), mm(mm),
        topo(NULL), conf(NULL), initialised(false), valid(false), molecular_virial(false) {
    for (const auto &site : sites) {
      if (site.size()==0) throw std::invalid_argument("Empty LE virtual endpoint");
      for (int i=0;i<site.size();++i) indices.push_back(site.atom(i));
    }
    if (type == 2 && sites.size() == 2) {
      distance.reset(new topology::distance_restraint_struct(sites[0],sites[1],0,1,0));
      auto *geometry = new interaction::Distance_Colvar();
      geometry->params = distance.get();
      cv.reset(geometry);
    } else if (type == 1 && sites.size() == 4) {
      dihedral.reset(new topology::dihedral_restraint_struct(sites[0].atom(0),sites[1].atom(0),sites[2].atom(0),sites[3].atom(0),0,0,1));
      auto *geometry = new interaction::Dihedral_Colvar();
      geometry->params = dihedral.get();
      geometry->endpoints = &this->sites;
      cv.reset(geometry);
    } else if (type == 3 && sites.size() == 3) {
      angle.reset(new topology::angle_restraint_struct(sites[0].atom(0),sites[1].atom(0),sites[2].atom(0),0,1));
      auto *geometry = new interaction::Angle_Colvar();
      geometry->params = angle.get();
      geometry->endpoints = &this->sites;
      cv.reset(geometry);
    } else if (type == 4 && sites.size() == 4) {
      difference.reset(new topology::distance_difference_restraint_struct(
          sites[0],sites[1],sites[2],sites[3]));
      auto *geometry = new interaction::Distance_Difference_Colvar();
      geometry->params = difference.get();
      cv.reset(geometry);
    } else if (type == 5 && group1_size > 0 && group1_size < sites.size()
               && std::isfinite(rcut) && rcut > 0 && nn > 0
               && nn <= std::numeric_limits<int>::max()/2 && mm > nn) {
      const std::vector<util::Virtual_Atom> first(sites.begin(),sites.begin()+group1_size);
      const std::vector<util::Virtual_Atom> second(sites.begin()+group1_size,sites.end());
      coordnum.reset(new topology::coordnum_restraint_struct(true,first,second,0,1,rcut,nn,mm));
      coordnum->masscut=0;
      auto *geometry = new interaction::Coordnum_Colvar();
      geometry->params = coordnum.get();
      cv.reset(geometry);
    } else {
      throw std::invalid_argument("Invalid LE CV geometry or CN switching parameters");
    }
  }
};

util::LE_Colvar_Coordinate::LE_Colvar_Coordinate(int umbrella, int type,
    const std::vector<unsigned int> &atoms) : m_impl(new Impl(type,atoms))
{ umbrella_id(umbrella); }
util::LE_Colvar_Coordinate::LE_Colvar_Coordinate(int umbrella,int type,
    const std::vector<util::Virtual_Atom> &sites,size_t group1_size,double rcut,int nn,int mm)
    : m_impl(new Impl(type,sites,group1_size,rcut,nn,mm)) { umbrella_id(umbrella); }
util::LE_Colvar_Coordinate::LE_Colvar_Coordinate(int umbrella,
    const std::vector<unsigned int> &group1,const std::vector<unsigned int> &group2,
    double rcut,int nn,int mm) {
  auto atoms=group1;
  atoms.insert(atoms.end(),group2.begin(),group2.end());
  m_impl.reset(new Impl(5,atoms,group1.size(),rcut,nn,mm));
  umbrella_id(umbrella);
}
util::LE_Colvar_Coordinate::~LE_Colvar_Coordinate() {}
util::LE_Coordinate *util::LE_Colvar_Coordinate::clone() const
{
  if (m_impl->site_definition)
    return new LE_Colvar_Coordinate(umbrella_id(),m_impl->type,m_impl->sites,
        m_impl->group1_size,m_impl->rcut,m_impl->nn,m_impl->mm);
  if (m_impl->type == 5) {
    const auto split=m_impl->indices.begin()+m_impl->group1_size;
    return new LE_Colvar_Coordinate(umbrella_id(),
        std::vector<unsigned int>(m_impl->indices.begin(),split),
        std::vector<unsigned int>(split,m_impl->indices.end()),m_impl->rcut,m_impl->nn,m_impl->mm);
  }
  return new LE_Colvar_Coordinate(umbrella_id(),m_impl->type,m_impl->indices);
}
int util::LE_Colvar_Coordinate::get_type() const { return m_impl->type; }
void util::LE_Colvar_Coordinate::calculate(configuration::Configuration &)
{
  m_impl->valid=false;
  io::messages.add("LE CV coordinates require calculate_checked with topology and simulation context.","LE_Colvar",io::message::error);
}
int util::LE_Colvar_Coordinate::calculate_checked(topology::Topology &topo,
    configuration::Configuration &conf, simulation::Simulation &sim)
{
  m_impl->valid=false;
  m_impl->molecular_virial=sim.param().pcouple.virial==math::molecular_virial;
  for (unsigned int atom : m_impl->indices)
    if (atom >= topo.num_atoms()) {
      io::messages.add("LE CV atom index out of range.","LE_Colvar",io::message::error);
      return 1;
    }
  if (!m_impl->initialised) {
    if (m_impl->cv->init(topo,conf,sim,std::cout,true)) return 1;
    m_impl->initialised=true;
  }
  if (m_impl->cv->calculate(topo,conf,sim)) return 1;
  double norm=0;
  for (const auto &derivative : m_impl->cv->derivatives) norm+=math::abs2(derivative);
  if (!std::isfinite(m_impl->cv->ct) || !std::isfinite(norm)
      || ((m_impl->type == 1 || m_impl->type == 3) && norm == 0)) {
    io::messages.add("LE CV has a non-finite value/derivative or singular dihedral/angle geometry.","LE_Colvar",io::message::error);
    return 1;
  }
  m_impl->topo=&topo; m_impl->conf=&conf; m_impl->valid=true;
  return 0;
}
double util::LE_Colvar_Coordinate::get_value(double grid_min, double grid_max) const
{
  if (!m_impl->valid) return std::numeric_limits<double>::quiet_NaN();
  double value=m_impl->cv->ct;
  if (m_impl->type == 1) {
    const double period=2*math::Pi;
    // Preserve LE's existing image convention, including partial grid windows.
    if (value < grid_min) value+=period*std::ceil((grid_min-value)/period);
    if (value >= grid_max) value-=period*(1+std::floor((value-grid_max)/period));
  }
  return value;
}
double util::LE_Colvar_Coordinate::get_deviation(const double &grid_value) const
{
  if (!m_impl->valid) return std::numeric_limits<double>::quiet_NaN();
  double delta=m_impl->cv->ct-grid_value;
  if (m_impl->type == 1) delta-=2*math::Pi*std::floor(delta/(2*math::Pi)+.5);
  return delta;
}
void util::LE_Colvar_Coordinate::apply(double derivative)
{ apply(derivative,false); }
void util::LE_Colvar_Coordinate::apply(double derivative,bool include_virial)
{
  if (!m_impl->valid || !std::isfinite(derivative)) {
    io::messages.add("Cannot apply LE CV force without valid evaluated geometry and finite derivative.","LE_Colvar",io::message::error);
    return;
  }
  m_impl->cv->apply_forces(*m_impl->topo,*m_impl->conf,derivative,include_virial,m_impl->molecular_virial);
}
std::string util::LE_Colvar_Coordinate::str() const
{
  std::ostringstream out;
  const char *names[]={"unknown","dihedral","distance","angle","distance difference","coordination number"};
  out << "LE CV " << names[m_impl->type] << " (";
  for (size_t i=0;i<m_impl->indices.size();++i) out << (i ? "-" : "") << m_impl->indices[i]+1;
  return out.str()+")";
}

std::string util::LE_Colvar_Coordinate::definition() const
{
  std::ostringstream out;
  out.setf(std::ios::scientific,std::ios::floatfield);
  out.precision(std::numeric_limits<double>::max_digits10);
  if (m_impl->site_definition) {
    out << "CV " << m_impl->type << " SITES " << m_impl->sites.size();
    for (const auto &site : m_impl->sites) {
      out << " | " << site.type() << ' ' << site.size();
      for (int i=0;i<site.size();++i) out << ' ' << site.atom(i)+1;
      out << ' ' << site.dish() << ' ' << site.disc() << ' ' << site.orientation();
    }
  } else {
  out << "CV " << m_impl->type << ' ' << m_impl->indices.size();
  for (unsigned int atom : m_impl->indices) out << ' ' << atom+1;
  }
  if (m_impl->type == 5) {
    out.setf(std::ios::scientific,std::ios::floatfield);
    out.precision(std::numeric_limits<double>::max_digits10);
    out << " CN " << m_impl->group1_size << ' ' << m_impl->rcut
        << ' ' << m_impl->nn << ' ' << m_impl->mm;
  }
  return out.str();
}
