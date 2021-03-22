/**
 * @file perturbed_atom.h
 * a perturbed atom.
 */

#ifndef INCLUDED_PERTURBED_ATOM_H
#define INCLUDED_PERTURBED_ATOM_H

namespace topology
{
  /**
   * @class Perturbed_Atom
   * holds the perturbed atom information.
   */
  class Perturbed_Atom
  {
  public:
    /**
     * Default Constructor.
     */
    Perturbed_Atom()
      : m_sequence_number(0),
	m_A_IAC(0),
	m_A_mass(0),
	m_A_charge(0),
        m_A_polarisability(0.0),
        m_A_damping_level(0.0),
  	m_B_IAC(0),
	m_B_mass(0),
	m_B_charge(0),
        m_B_polarisability(0.0),
        m_B_damping_level(0.0),
       	m_LJ_softcore(0),
        m_crf_softcore(0),
        m_exclusion(),
        m_one_four_pair()
    {};

    /**
     * Constructor.
     * @param JLA atom sequence number.
     * @param IACA integer atom code of state A.
     * @param WMA mass of state A.
     * @param CGA charge of state A.
     * @param IACB integer atom code of state B.
     * @param WMB mass of state B.
     * @param CGB charge of state B.
     * @param SCLJ soft core van der Waals parameter.
     * @param SCC soft core electrostatic parameter.
     */
    Perturbed_Atom(unsigned int JLA, unsigned int IACA,
		   double WMA, double CGA,
                   double pol_A, double lev_A,
		   unsigned int IACB, double WMB,
		   double CGB,
                   double pol_B, double lev_B,
		   double SCLJ, double SCC)
      : m_sequence_number(JLA),
	m_A_IAC(IACA),
	m_A_mass(WMA),
	m_A_charge(CGA),
        m_A_polarisability(pol_A),
        m_A_damping_level(lev_A),
	m_B_IAC(IACB),
	m_B_mass(WMB),
	m_B_charge(CGB),
        m_B_polarisability(pol_B),
        m_B_damping_level(lev_B),
	m_LJ_softcore(SCLJ),
        m_crf_softcore(SCC),
        m_exclusion(),
        m_one_four_pair()
    {};
  public:
    /**
     * @name accessors
     * @{
     */
    unsigned int sequence_number();
    unsigned int sequence_number()const;
    void sequence_number(unsigned int);
    unsigned int A_IAC()const;
    void A_IAC( unsigned int);
    
    double A_mass()const;
    void A_mass(double);
    
    double A_charge()const;
    void A_charge(double);
    
    double A_polarisability()const;
    void A_polarisability(double);
    
    double A_damping_level()const;
    void A_damping_level(double);
   
    unsigned int B_IAC()const;
    void B_IAC(unsigned int);
    
    double B_mass()const;
    void B_mass(double);
    
    double B_charge()const;
    void B_charge(double);
 
    double B_polarisability()const;
    void B_polarisability(double);
    
    double B_damping_level()const;
    void B_damping_level(double);
   
    double LJ_softcore()const;
    void LJ_softcore(double);
    
    double CRF_softcore()const;
    void CRF_softcore(double);
    
    std::set<int> & exclusion();
    std::set<int> const & exclusion()const;
    
    std::set<int> & one_four_pair();
    std::set<int> const & one_four_pair()const;
    /**
     * @}
     */

  private:
    unsigned int m_sequence_number;
    unsigned int m_A_IAC;
    double m_A_mass;
    double m_A_charge;
    double m_A_polarisability;
    double m_A_damping_level;
    unsigned int m_B_IAC;
    double m_B_mass;
    double m_B_charge;
    double m_B_polarisability;
    double m_B_damping_level;
    double m_LJ_softcore;
    double m_crf_softcore;
    std::set<int> m_exclusion;
    std::set<int> m_one_four_pair;
  };
  
} // topology

inline unsigned int topology::Perturbed_Atom::sequence_number()const{
  return m_sequence_number;
}

inline unsigned int topology::Perturbed_Atom::A_IAC()const{
  return m_A_IAC;
}

inline double topology::Perturbed_Atom::A_mass()const{
  return m_A_mass; 
}

inline double topology::Perturbed_Atom::A_charge()const{
  return m_A_charge; 
}

inline double topology::Perturbed_Atom::A_polarisability()const{
  return m_A_polarisability; 
}

inline double topology::Perturbed_Atom::A_damping_level()const{
  return m_A_damping_level; 
}

inline unsigned int topology::Perturbed_Atom::B_IAC()const{
  return m_B_IAC; 
}

inline double topology::Perturbed_Atom::B_mass()const{
  return m_B_mass; 
}

inline double topology::Perturbed_Atom::B_charge()const{
  return m_B_charge; 
}

inline double topology::Perturbed_Atom::B_polarisability()const{
  return m_B_polarisability; 
}

inline double topology::Perturbed_Atom::B_damping_level()const{
  return m_B_damping_level; 
}

inline double topology::Perturbed_Atom::LJ_softcore()const{
  return m_LJ_softcore; 
}

inline double topology::Perturbed_Atom::CRF_softcore()const{
  return m_crf_softcore; 
}

inline void topology::Perturbed_Atom::sequence_number(unsigned int a){
  m_sequence_number = a;
}
inline void topology::Perturbed_Atom::A_IAC(unsigned int a){
  m_A_IAC = a;
}
inline void topology::Perturbed_Atom::A_mass(double a){ 
  m_A_mass = a;
}
inline void topology::Perturbed_Atom::A_charge(double a){
  m_A_charge = a;
}
inline void topology::Perturbed_Atom::A_polarisability(double a){ 
  m_A_polarisability = a;
}
inline void topology::Perturbed_Atom::A_damping_level(double a){ 
  m_A_damping_level = a;
}
inline void topology::Perturbed_Atom::B_IAC(unsigned int a){
  m_B_IAC = a;
}
inline void topology::Perturbed_Atom::B_mass(double a){
  m_B_mass = a;
}
inline void topology::Perturbed_Atom::B_charge(double a){
  m_B_charge = a;
}
inline void topology::Perturbed_Atom::B_polarisability(double a){ 
  m_B_polarisability = a;
}
inline void topology::Perturbed_Atom::B_damping_level(double a){ 
  m_B_damping_level = a;
}
inline void topology::Perturbed_Atom::LJ_softcore(double a){
  m_LJ_softcore = a;
}
inline void topology::Perturbed_Atom::CRF_softcore(double a){
  m_crf_softcore = a;
}
inline std::set<int> & topology::Perturbed_Atom::exclusion()
{
  return m_exclusion;
}
inline std::set<int> const & topology::Perturbed_Atom::exclusion()const
{
  return m_exclusion;
}

inline std::set<int> & topology::Perturbed_Atom::one_four_pair()
{
  return m_one_four_pair;
}

inline std::set<int> const & topology::Perturbed_Atom::one_four_pair()const
{
  return m_one_four_pair;
}


namespace topology
{
  /**
   * @class EDS_Perturbed_Atom
   * holds the EDS-perturbed atom information.
   */
  class EDS_Perturbed_Atom
  {
  public:
    /**
     * Default Constructor.
     */
    EDS_Perturbed_Atom()
          : m_sequence_number(0),
            m_M_IAC(),
            m_M_charge(),
            m_LJ_softcore(0),
            m_crf_softcore(0),
            m_exclusion(),
            m_one_four_pair()
    {};

    /**
     * Constructor.
     * @param JLA atom sequence number.
     * @param IACM integer atom code of states 0...N.
     * @param CGM  charge of states 0...N.
     * @param SCLJ soft core LJ parameter
     * @param SCC soft core CRF parameter
     */
    EDS_Perturbed_Atom(unsigned int JLA,
            std::vector<unsigned int> IACM,
            std::vector<double> CGM,
            double SCLJ, double SCC)
          : m_sequence_number(JLA),
            m_M_IAC(IACM),
            m_M_charge(CGM),
            m_LJ_softcore(SCLJ),
            m_crf_softcore(SCC),
            m_exclusion(),
            m_one_four_pair()
    {};

  public:
    /**
     * @name accessors
     * @{
     */
    unsigned int sequence_number();
    unsigned int sequence_number()const;
    void sequence_number(unsigned int);
    
    std::vector<unsigned int> M_IAC()const;
    void M_IAC(std::vector<unsigned int>);

    std::vector<double> M_charge()const;
    void M_charge(std::vector<double>);
    
    double LJ_softcore()const;
    void LJ_softcore(double);
    
    double CRF_softcore()const;
    void CRF_softcore(double);
    
    std::set<int> & exclusion();
    std::set<int> const & exclusion()const;
    
    std::set<int> & one_four_pair();
    std::set<int> const & one_four_pair()const;
    
    /**
     * @}
     */

  private:
    unsigned int m_sequence_number;
    std::vector<unsigned int> m_M_IAC;
    std::vector<double> m_M_charge;
    double m_LJ_softcore;
    double m_crf_softcore;
    std::set<int> m_exclusion;
    std::set<int> m_one_four_pair;
  };
  
} // topology

inline unsigned int topology::EDS_Perturbed_Atom::sequence_number()const{
  return m_sequence_number;
}
inline std::vector<unsigned int> topology::EDS_Perturbed_Atom::M_IAC()const{
  return m_M_IAC; 
}
inline std::vector<double> topology::EDS_Perturbed_Atom::M_charge()const{
  return m_M_charge; 
}
inline double topology::EDS_Perturbed_Atom::LJ_softcore()const{
  return m_LJ_softcore; 
}
inline double topology::EDS_Perturbed_Atom::CRF_softcore()const{
  return m_crf_softcore; 
}
inline std::set<int> const & topology::EDS_Perturbed_Atom::exclusion()const{
  return m_exclusion;
}
inline std::set<int> const & topology::EDS_Perturbed_Atom::one_four_pair()const{
  return m_one_four_pair;
}

inline void topology::EDS_Perturbed_Atom::sequence_number(unsigned int a){
  m_sequence_number = a;
}
inline void topology::EDS_Perturbed_Atom::M_IAC(std::vector<unsigned int> a){
  m_M_IAC = a;
}
inline void topology::EDS_Perturbed_Atom::M_charge(std::vector<double> a){
  m_M_charge = a;
}
inline void topology::EDS_Perturbed_Atom::LJ_softcore(double a){
  m_LJ_softcore = a;
}
inline void topology::EDS_Perturbed_Atom::CRF_softcore(double a){
  m_crf_softcore = a;
}
inline std::set<int> & topology::EDS_Perturbed_Atom::exclusion()
{
  return m_exclusion;
}
inline std::set<int> & topology::EDS_Perturbed_Atom::one_four_pair()
{
  return m_one_four_pair;
}

#endif
  
