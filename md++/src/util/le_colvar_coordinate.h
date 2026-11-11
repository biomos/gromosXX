/** LE coordinate adapter for the existing geometric CV implementations. */
#ifndef INCLUDED_LE_COLVAR_COORDINATE_H
#define INCLUDED_LE_COLVAR_COORDINATE_H
#include <memory>
#include <vector>
#include "le_coordinate.h"
namespace util {
class Virtual_Atom;
class LE_Colvar_Coordinate : public LE_Coordinate {
public:
  // LE types: 1 dihedral, 2 distance, 3 angle, 4 distance difference, 5 CN.
  // These are LE identifiers, not COLVARRES TYPE numbers.
  LE_Colvar_Coordinate(int umbrella, int type, const std::vector<unsigned int> &atoms);
  LE_Colvar_Coordinate(int umbrella, const std::vector<unsigned int> &group1,
      const std::vector<unsigned int> &group2, double rcut, int nn, int mm);
  LE_Colvar_Coordinate(int umbrella,int type,const std::vector<Virtual_Atom> &sites,
      size_t group1_size=0,double rcut=0,int nn=0,int mm=0);
  ~LE_Colvar_Coordinate();
  LE_Coordinate *clone() const;
  int get_type() const;
  void calculate(configuration::Configuration &conf);
  int calculate_checked(topology::Topology &topo, configuration::Configuration &conf,
                        simulation::Simulation &sim);
  double get_value(double grid_min, double grid_max) const;
  double get_deviation(const double &grid_value) const;
  void apply(double deriv);
  void apply(double deriv,bool include_virial);
  bool supports_virial() const { return true; }
  std::string str() const;
  std::string definition() const;
private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
  LE_Colvar_Coordinate(const LE_Colvar_Coordinate &) = delete;
  LE_Colvar_Coordinate &operator=(const LE_Colvar_Coordinate &) = delete;
};
}
#endif
