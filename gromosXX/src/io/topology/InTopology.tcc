/**
 * @file InTopology.tcc
 * implements methods of InTopology.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

#include "../../debug.h"

inline io::InTopology &io::InTopology::operator>>(simulation::Topology& topo){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  
  { // BONDH
    buffer.clear();
    buffer = m_block["BONDH"];
  
    it = buffer.begin() + 1;

    _lineStream.clear();
    _lineStream.str(*it);

    int num, n;
    _lineStream >> num;
    ++it;

    for(n=0; it != buffer.end() - 1; ++it, ++n){
      int i, j, t;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i >> j >> t;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in BONDH block");

      topo.solute().bonds().push_back(simulation::Bond(i-1, j-1, t-1));
    }
    
    if(n != num){
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("error in BONDH block (n != num)");
    }
  } // BONDH

  { // BOND
    buffer = m_block["BOND"];
  
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    int num, n;
    _lineStream >> num;
    ++it;
    
    for(n=0; it != buffer.end() - 1; ++it, ++n){
      int i, j, t;
      
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i >> j >> t;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in BOND block");
      
      topo.solute().bonds().push_back(simulation::Bond(i-1, j-1, t-1));
    }
    
    if(n != num){
      if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in BOND block (n != num)");
    }
  } // BOND

  { // BONDANGLE
    buffer = m_block["BONDANGLE"];
  
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    int num, n;
    _lineStream >> num;
    ++it;
    
    for(n=0; it != buffer.end() - 1; ++it, ++n){
      int i, j, k, t;
      
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i >> j >> k >> t;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in BONDANGLE block");
      
      topo.solute().angles().push_back(simulation::Angle(i-1, j-1, k-1, t-1));
    }
    
    if(n != num){
      if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in BONDANGLE block (n != num)");
    }
  } // BONDANGLE

  { // BONDANGLEH
    buffer.clear();
    buffer = m_block["BONDANGLEH"];
  
    it = buffer.begin() + 1;

    _lineStream.clear();
    _lineStream.str(*it);

    int num, n;
    _lineStream >> num;
    ++it;

    for(n=0; it != buffer.end() - 1; ++it, ++n){
      int i, j, k, t;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i >> j >> k >> t;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in BONDANGLEH block");

      topo.solute().angles().push_back(simulation::Angle(i-1, j-1, k-1, t-1));
    }
    
    if(n != num){
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("error in BONDANGLEH block (n != num)");
    }
  } // BONDANGLEH

  { // IMPDIHEDRAL
    buffer = m_block["IMPDIHEDRAL"];
  
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    int num, n;
    _lineStream >> num;
    ++it;
    
    for(n=0; it != buffer.end() - 1; ++it, ++n){
      int i, j, k, l, t;
      
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i >> j >> k >> l >> t;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in IMPDIHEDRAL block");
      
      topo.solute().improper_dihedrals().
	push_back(simulation::Improper_Dihedral(i-1, j-1, k-1, l-1, t-1));
    }
    
    if(n != num){
      if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in IMPDIHEDRAL block (n != num)");
    }
  } // IMPDIHEDRAL

  { // IMPDIHEDRALH
    buffer.clear();
    buffer = m_block["IMPDIHEDRALH"];
  
    it = buffer.begin() + 1;

    _lineStream.clear();
    _lineStream.str(*it);

    int num, n;
    _lineStream >> num;
    ++it;

    for(n=0; it != buffer.end() - 1; ++it, ++n){
      int i, j, k, l, t;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i >> j >> k >> l >> t;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in IMPDIHEDRALH block");

      topo.solute().improper_dihedrals().
	push_back(simulation::Improper_Dihedral(i-1, j-1, k-1, l-1, t-1));
    }
    
    if(n != num){
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("error in IMPDIHEDRALH block (n != num)");
    }
  } // IMPDIHEDRALH

  { // DIHEDRAL
    buffer = m_block["DIHEDRAL"];
  
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    int num, n;
    _lineStream >> num;
    ++it;
    
    for(n=0; it != buffer.end() - 1; ++it, ++n){
      int i, j, k, l, t;
      
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i >> j >> k >> l >> t;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in DIHEDRAL block");
      
      topo.solute().dihedrals().
	push_back(simulation::Dihedral(i-1, j-1, k-1, l-1, t-1));
    }
    
    if(n != num){
      if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in DIHEDRAL block (n != num)");
    }
  } // DIHEDRAL

  { // DIHEDRALH
    buffer.clear();
    buffer = m_block["DIHEDRALH"];
  
    it = buffer.begin() + 1;

    _lineStream.clear();
    _lineStream.str(*it);

    int num, n;
    _lineStream >> num;
    ++it;

    for(n=0; it != buffer.end() - 1; ++it, ++n){
      int i, j, k, l, t;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> i >> j >> k >> l >> t;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	io::messages.add("bad line in DIHEDRALH block",
			 "InTopology", io::message::error);

      topo.solute().dihedrals().
	push_back(simulation::Dihedral(i-1, j-1, k-1, l-1, t-1));
    }
    
    if(n != num){
      if (_lineStream.fail() || ! _lineStream.eof())
	io::messages.add("error in DIHEDRALH block (n != num)",
			 "InTopology", io::message::error);
    }
  } // DIHEDRALH
  
  { // RESNAME
    buffer = m_block["RESNAME"];
    it = buffer.begin()+1;
    int n, num;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num;
    ++it;
    
    for(n=0; it != buffer.end() - 1; ++it, ++n){
      std::string s;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> s;
      
      topo.residue_name().push_back(s);
      DEBUG(10, "RESNAME: " << s);
    }

    if (n != num){
      io::messages.add("Error in RESNAME block: n!=num.",
		       "InTopology", io::message::critical);
      throw std::runtime_error("error in RESNAME block (n != num)");
    }
    
  } // RESNAME
  
  { // SOLUTEATOM
    buffer = m_block["SOLUTEATOM"];
  
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    int num, n;
    _lineStream >> num;
    topo.resize(num);

    // put the rest of the block into a single stream
    ++it;
    std::string soluteAtoms;
    std::vector<std::string>::const_iterator bla = it;
    std::vector<std::string>::const_iterator fasel = buffer.end() - 1;
    concatenate(bla, fasel, soluteAtoms);
    _lineStream.clear();
    _lineStream.str(soluteAtoms);

    int a_nr, r_nr, t, cg, n_ex, a_ex;
    double m, q;
    std::string s;
    std::set<int> ex;
    std::set<int> ex14;
      
    for(n=0; n < num; ++n){

      _lineStream >> a_nr >> r_nr >> s >> t >> m >> q >> cg >> n_ex;

      // exclusions
      ex.clear();
      for(int i=0; i<n_ex; ++i){
	_lineStream >> a_ex;
	ex.insert(a_ex-1);
      }
      
      // 1,4 - pairs
      _lineStream >> n_ex;
      ex14.clear();
      for(int i=0; i<n_ex; ++i){
	_lineStream >> a_ex;
	ex14.insert(a_ex-1);
      }
      
      if (_lineStream.fail())
	throw std::runtime_error("bad line in SOLUTEATOM block");

      topo.add_solute_atom(s, r_nr-1, t-1, m, q, cg, ex, ex14);
    }
  } // SOLUTEATOM
  
  { // SOLVENTATOM and SOLVENTCONSTR
    // give it a number (SOLVENTATOM1, SOLVENTATOM2) for multiple
    // solvents...
    buffer = m_block["SOLVENTATOM"];
    
    int res_nr = topo.residue_name().size();
    DEBUG(10, "res name: " << topo.residue_name().size());
    for(int i=0; i<res_nr; ++i){
      DEBUG(10, topo.residue_name()[i]);
    }
    
    topo.residue_name().push_back("SOLV");

    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    int num, n;
    _lineStream >> num;
    ++it;
    
    simulation::Solvent s;
    
    std::string name;
    int i, iac;
    double mass, charge;
    
    for(n=0; it != buffer.end()-1; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> i >> name >> iac >> mass >> charge;

      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in SOLVENTATOM block");

      s.add_atom(name, res_nr, iac-1, mass, charge);
    }
    
    if (n!=num)
      throw std::runtime_error("error in SOLVENTATOM block (n != num)");

    buffer = m_block["SOLVENTCONSTR"];
    
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> num;
    ++it;
    
    int j;
    double b0;
    
    for(n=0; it != buffer.end()-1; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> i >> j >> b0;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in SOLVENTCONSTR block");
 
      s.add_distance_constraint(i-1, j-1, b0);
    }

    if (n!=num)
      throw std::runtime_error("error in SOLVENTATOM block (n != num)");

    topo.add_solvent(s);
  }

  return *this;
}

template<typename t_simulation>
io::InTopology &io::InTopology
::operator>>(interaction::harmonic_bond_interaction<t_simulation> &hbi){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["BONDTYPE"];

  io::messages.add("converting bond force constants from quartic to harmonic form", "InTopology::bondtype", io::message::notice);

  // 1. BONDTYPE 2. number of types
  for (it = buffer.begin() + 2; 
   it != buffer.end() - 1; ++it) {

    double k, r;
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> k >> r;

    if (_lineStream.fail() || ! _lineStream.eof())
      throw std::runtime_error("bad line in BONDTYPE block");

    // we are reading into harmonic bond term, so convert k
    k *= 2 * r * r;

    // and add...
    hbi.add(k, r);

  }
  

  return *this;
}
template<typename t_simulation>
io::InTopology &io::InTopology
::operator>>(interaction::Quartic_bond_interaction<t_simulation> &qbi){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["BONDTYPE"];

  // 1. BONDTYPE 2. number of types
  for (it = buffer.begin() + 2; 
   it != buffer.end() - 1; ++it) {

    double k, r;
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> k >> r;

    if (_lineStream.fail() || ! _lineStream.eof())
      throw std::runtime_error("bad line in BONDTYPE block");

    // and add...
    qbi.add(k, r);

  }
  

  return *this;
}

template<typename t_simulation>
io::InTopology &io::InTopology
::operator>>(interaction::angle_interaction<t_simulation> &ai){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["BONDANGLETYPE"];

  // 1. BONDTYPE 2. number of types
  for (it = buffer.begin() + 2; 
   it != buffer.end() - 1; ++it) {

    double k, t;
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> k >> t;

    if (_lineStream.fail() || ! _lineStream.eof())
      throw std::runtime_error("bad line in BONDANGLETYPE block");

    // and add...
    ai.add(k, cos( t * 2 * math::Pi / 360.0 ) );

  }
  return *this;
}

template<typename t_simulation>
io::InTopology &io::InTopology
::operator>>(interaction::Improper_dihedral_interaction<t_simulation> &ii){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["IMPDIHEDRALTYPE"];

  // 1. IMPDIHEDRALTYPE 2. number of types
  for (it = buffer.begin() + 2; 
   it != buffer.end() - 1; ++it) {

    double k, q;
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> k >> q;

    if (_lineStream.fail() || ! _lineStream.eof())
      throw std::runtime_error("bad line in IMPDIHEDRALTYPE block");

    // and add...
    ii.add(k*180*180/math::Pi/math::Pi , q * math::Pi / 180);

  }
  return *this;
}
template<typename t_simulation>
io::InTopology &io::InTopology
::operator>>(interaction::Dihedral_interaction<t_simulation> &di){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["DIHEDRALTYPE"];

  // 1. DIHEDRALTYPE 2. number of types
  for (it = buffer.begin() + 2; 
   it != buffer.end() - 1; ++it) {

    double k, pd;
    int m;
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> k >> pd >> m;

    if (_lineStream.fail() || ! _lineStream.eof())
      io::messages.add("bad line in DIHEDRALTYPE block", 
		       "InTopology", io::message::error);

    // and add...
    di.add(k, pd, m);
    DEBUG(10, "DIHEDRALTYPE " << k << " " << pd << " " << m);
  }

  return *this;
}

template<typename t_simulation, typename t_pairlist, typename t_innerloop>
io::InTopology &io::InTopology
::operator>>(interaction
	     ::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop> &nbi){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  { // TOPPHYSCON
    buffer = m_block["TOPPHYSCON"];
    
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    double fpepsi, hbar;
    _lineStream >> fpepsi;
    if (_lineStream.fail())
      io::messages.add("Bad line in TOPPHYSCON block",
			"InTopology", io::message::error);
    nbi.coulomb_constant(fpepsi);
  } // TOPPHYSCON

  { // LJPARAMETERS
    
    buffer = m_block["LJPARAMETERS"];
    
    int num, n;
    
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num;
    
    // calculate the matrix size from: x = n*(n+1)/2
    size_t sz = size_t(sqrt(double((8*num+1)-1))/2);
    nbi.resize(sz);
    
    ++it;
    
    for (n=0; it != buffer.end() - 1; ++it, ++n) {
      
      typename interaction::lj_parameter_struct s;
      int i, j;
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> i >> j >> s.c12 >> s.c6 >> s.cs12 >> s.cs6;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in LJPARAMETERS block");
      
      // and add...
      nbi.add_lj_parameter(i-1, j-1, s);
    } 
    
  

    if (num != n){
      io::messages.add("Reading the LJPARAMETERS failed (n != num)",
		       "InTopology",
		       io::message::critical);
      throw std::runtime_error("bad LJPARAMETERS block (n != num)");
    }
  } // LJPARAMETER
   
  return *this;  
}

