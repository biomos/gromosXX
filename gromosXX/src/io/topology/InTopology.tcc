/**
 * @file InTopology.tcc
 * implements methods of InTopology.
 */

inline io::InTopology &io::InTopology::operator>>(simulation::topology& topo){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
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
      
      topo.bonds().add(i-1, j-1, t-1);
    }
    
    if(n != num){
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("error in BOND block (n != num)");
    }
  } // BOND

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

      topo.bonds().add(i-1, j-1, t-1);
    }
    
    if(n != num){
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("error in HBOND block (n != num)");
    }
  } // BONDH
  
  { // RESNAME
    buffer = m_block["RESNAME"];
    buffer.erase(buffer.begin());
    buffer.erase(buffer.begin());
    buffer.erase(buffer.end()-1);
    
    topo.residue_name() = buffer;
    
  } // RESNAME
  
  { // SOLUTEATOM
    buffer = m_block["SOLUTEATOM"];
  
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    int num, n;
    _lineStream >> num;
    ++it;
    
    topo.solute_atoms_capacity(num);

    for(n=0; it != buffer.end() - 1; ++it, ++n){
      int a_nr, r_nr, t, cg, n_ex, a_ex;
      double m, q;
      std::string s;
      std::set<int> ex;
      std::set<int> ex14;
      
      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> a_nr >> r_nr >> s >> t >> m >> q >> cg >> n_ex;
      for(int i=0; i<n_ex; ++i){
	_lineStream >> a_ex;
	ex.insert(a_ex);
      }
      ++it;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> n_ex;
      for(int i=0; i<n_ex; ++i){
	_lineStream >> a_ex;
	ex14.insert(a_ex);
      }
      
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line SOLUTEATOM block");

      topo.add_solute_atom(s, r_nr-1, t-1, m, q, cg, ex, ex14);

    }
    
    if(n != num){
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("error in SOLUTEATOM block (n != num)");
    }
  } // SOLUTEATOM

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
      throw std::runtime_error("bad line in bond type block");

    // we are reading into harmonic bond term, so convert k
    k *= 2 * r * r;

    // and add...
    hbi.add(k, r);

  }
  

  return *this;
}

inline void io::InTopology::read_stream()
{
  std::vector<std::string> buffer;
  
  while(!stream().eof()){

    try{
      io::getblock(stream(), buffer);
    }
    catch(std::runtime_error e){
      break;
    }
    
    m_block[buffer[0]] = buffer;    
    buffer.clear();
    
  }
}
