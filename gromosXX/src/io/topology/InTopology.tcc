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
      
      topo.solute().bonds().add(i-1, j-1, t-1);
    }
    
    if(n != num){
      if (_lineStream.fail()|| ! _lineStream.eof())
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

      topo.solute().bonds().add(i-1, j-1, t-1);
    }
    
    if(n != num){
      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("error in HBOND block (n != num)");
    }
  } // BONDH
  
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

      for(int i=0; i<n_ex; ++i){
	_lineStream >> a_ex;
	ex.insert(a_ex);
      }
      
      _lineStream >> n_ex;
      for(int i=0; i<n_ex; ++i){
	_lineStream >> a_ex;
	ex14.insert(a_ex);
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

template<typename t_simulation, typename t_pairlist>
io::InTopology &io::InTopology
::operator>>(interaction
	     ::Nonbonded_Interaction<t_simulation, t_pairlist> &nbi){

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

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
