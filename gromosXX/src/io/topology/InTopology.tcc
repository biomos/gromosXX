/**
 * @file InTopology.tcc
 * implements methods of InTopology.
 */

inline io::InTopology &io::InTopology::operator>>(simulation::topology& topo){

  std::vector<std::string> soluteatomblock;
  std::vector<std::string> bondHblock;
  std::vector<std::string> bondblock;

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  while (!soluteatomblock.size()
    && ! bondHblock.size()
    && ! bondblock.size()
  ) {

    io::getblock(stream(), buffer);

    if (buffer[0] == "SOLUTEATOM")
      soluteatomblock = buffer;
    else if (buffer[0] == "BONDH")
      bondHblock = buffer;
    else if (buffer[0] == "BOND")
      bondblock = buffer;
  }  

  return *this;
}

template<typename t_simulation>
io::InTopology &io::InTopology
::operator>>(interaction::harmonic_bond_interaction<t_simulation> &hbi){

  std::vector<std::string> bondTypeblock;

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  while (! bondTypeblock.size()) {

    io::getblock(stream(), buffer);

    if (buffer[0] == "BONDTYPE")
      bondTypeblock = buffer;
  }

  for (it = bondTypeblock.begin() + 1; 
   it != bondTypeblock.end() - 1; it++) {

    double k, r;
    _lineStream.clear();
    _lineStream.str(*it);

    _lineStream >> k >> r;

    if (_lineStream.fail() || ! _lineStream.eof())
      throw std::runtime_error("bad line in bond type block");

    // and add...
    hbi.add(k, r);

  }
  

  return *this;
}
