/**
 * @file perturbed_atom.tcc
 * template methods of Perturbed_Atom
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

inline size_t simulation::Perturbed_Atom::sequence_number(){
  return m_sequence_number;
}
inline size_t const simulation::Perturbed_Atom::sequence_number()const{
  return m_sequence_number;
}
inline size_t simulation::Perturbed_Atom::A_IAC(){
  return m_A_IAC;
}
inline size_t const simulation::Perturbed_Atom::A_IAC()const{
  return m_A_IAC;
}
inline double simulation::Perturbed_Atom::A_mass(){
  return m_A_mass; 
}
inline double const simulation::Perturbed_Atom::A_mass()const{
  return m_A_mass; 
}
inline double simulation::Perturbed_Atom::A_charge(){
  return m_A_charge; 
}
inline double const simulation::Perturbed_Atom::A_charge()const{
  return m_A_charge; 
}
inline size_t simulation::Perturbed_Atom::B_IAC(){
  return m_B_IAC; 
}
inline size_t const simulation::Perturbed_Atom::B_IAC()const{
  return m_B_IAC;
}
inline double simulation::Perturbed_Atom::B_mass(){
  return m_B_mass; 
}
inline double const simulation::Perturbed_Atom::B_mass()const{
  return m_B_mass; 
}
inline double simulation::Perturbed_Atom::B_charge(){
  return m_B_charge; 
}
inline double const simulation::Perturbed_Atom::B_charge()const{
  return m_B_charge; 
}
inline double simulation::Perturbed_Atom::LJ_softcore(){
  return m_LJ_softcore; 
}
inline double const simulation::Perturbed_Atom::LJ_softcore()const{
  return m_LJ_softcore; 
}
inline double simulation::Perturbed_Atom::CRF_softcore(){
  return m_crf_softcore; 
}
inline double const simulation::Perturbed_Atom::CRF_softcore()const{
  return m_crf_softcore; 
}

inline void simulation::Perturbed_Atom::sequence_number(const size_t a){
  m_sequence_number = a;
}
inline void simulation::Perturbed_Atom::A_IAC(const size_t a){
  m_A_IAC = a;
}
inline void simulation::Perturbed_Atom::A_mass(const double a){ 
  m_A_mass = a;
}
inline void simulation::Perturbed_Atom::A_charge(const double a){
  m_A_charge = a;
}
inline void simulation::Perturbed_Atom::B_IAC(const size_t a){
  m_B_IAC = a;
}
inline void simulation::Perturbed_Atom::B_mass(const double a){
  m_B_mass = a;
}
inline void simulation::Perturbed_Atom::B_charge(const double a){
  m_B_charge = a;
}
inline void simulation::Perturbed_Atom::LJ_softcore(const double a){
  m_LJ_softcore = a;
}
inline void simulation::Perturbed_Atom::CRF_softcore(const double a){
  m_crf_softcore = a;
}
inline std::set<int> & simulation::Perturbed_Atom::exclusion()
{
  return m_exclusion;
}
inline std::set<int> const & simulation::Perturbed_Atom::exclusion()const
{
  return m_exclusion;
}

inline std::set<int> & simulation::Perturbed_Atom::one_four_pair()
{
  return m_one_four_pair;
}

inline std::set<int> const & simulation::Perturbed_Atom::one_four_pair()const
{
  return m_one_four_pair;
}

