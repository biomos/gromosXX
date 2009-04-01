/**
 * @file in_topology.cc
 * implements methods of In_Topology.
 */


#include <stdheader.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <interaction/interaction_types.h>
#include <io/instream.h>
#include <util/parse_tcouple.h>

#include <io/blockinput.h>

#include "in_topology.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

template<typename T>
bool check_type(std::vector<std::string> const & buffer, std::vector<T> term)
{
  if (buffer.size()){
    std::istringstream is(buffer[1]);
    int num;
    if (!(is >> num) || num < 0){
      std::cout << "ERROR:\tcould not read num or smaller 0!"
		<< std::endl;
      return false;
    }
    
    for(typename std::vector<T>::const_iterator
	  it = term.begin(),
	  to = term.end();
	it != to;
	++it){
      
      if (int(it->type) >= num){
	std::cout << "ERROR:\tused type " << it->type + 1
		  << " larger than max type (" << num << ")" 
		  << std::endl;
	
	return false;
      }
    }
    return true;
  }
  if (term.size() == 0) return true;

  std::cout << "ERROR:\tblock not found!" << std::endl;
  return false;
}

void 
io::In_Topology::read(topology::Topology& topo,
		      simulation::Parameter &param,
		      std::ostream & os){

  DEBUG(7, "reading in topology");

  if (!quiet){
    os << "TOPOLOGY\n";
    os << title << "\n";
  }
  
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  std::string bond_bname = "BONDTYPE";
  if (param.force.bond == 2){
    if(m_block["HARMBONDTYPE"].size())
      bond_bname = "HARMBONDTYPE";
  }
  if (param.force.bond == 0 &&
          m_block["BONDTYPE"].size() == 0 &&
          m_block["HARMBONDTYPE"].size() > 0)
    bond_bname = "HARMBONDTYPE";
 
  if(m_block["BONDSTRETCHTYPE"].size())
    bond_bname = "BONDSTRETCHTYPE";
  
  std::string angle_bname = "BONDANGLETYPE";
  // no automatic conversion for angles!
  if (param.force.angle == 2)
    angle_bname = "HARMBONDANGLETYPE";
  if(m_block["BONDANGLEBENDTYPE"].size())
    angle_bname = "BONDANGLEBENDTYPE";
  
  std::string dihedral_bname = "DIHEDRALTYPE";
  if(m_block["TORSDIHEDRALTYPE"].size())
    dihedral_bname = "TORSDIHEDRALTYPE";
  
  {
    buffer = m_block["TYPE"];
    if (buffer.size()){
      if (buffer.size() != 3){
	io::messages.add("Bad line in TYPE block",
			 "InTopology", io::message::error);
      }

      _lineStream.clear();
      _lineStream.str(buffer[1]);
      std::string s;
      _lineStream >> s;
      
      if (_lineStream.fail())
	io::messages.add("Bad line in TYPE block",
			 "InTopology", io::message::error);

      if (s.length() > MAX_NAME) {
        std::ostringstream msg;
        msg << "Error in TYPE block: type " << s 
            << " is too long (> " << MAX_NAME << " characters).";
        io::messages.add(msg.str(), "InTopology", io::message::error);
      }

      
      std::transform(s.begin(), s.end(), s.begin(), tolower);
      
      if (s == "atomistic"){
	if (!quiet)
	  os << "\tatomistic topology\n";
	
	if (param.force.interaction_function != simulation::lj_crf_func &&
            param.force.interaction_function != simulation::pol_lj_crf_func){
	  io::messages.add("wrong interaction function selected for atomistic topology",
			   "InTopology", io::message::error);
	}
      }
      else if (s == "coarse-grained"){
	if (!quiet)
	  os << "\tcoarse-grained topology\n";
	
	if (param.force.interaction_function != simulation::cgrain_func){
	  io::messages.add("wrong interaction function selected for coarse-grained topology",
			   "InTopology", io::message::error);
	}
      }
      else{
	io::messages.add("TYPE block: unknown topology type",
			 "InTopology", io::message::error);
      }
    }
    else{
      if (!quiet){
	os << "\tunknown topology type (atomistic / coarse-grained)\n";
	if (param.force.interaction_function == simulation::lj_crf_func ||
            param.force.interaction_function == simulation::pol_lj_crf_func)
	  os << "\tusing atomistic parameters\n";
	else if (param.force.interaction_function == simulation::cgrain_func)
	  os << "\tusing coarse-grained parameters\n";
	else
	  os << "\tunknown interaction function selected!\n";
      }
    }
  }

  { // PHYSICALCONSTANTS
    buffer = m_block["PHYSICALCONSTANTS"];
    if(buffer.size()){
      std::string s;
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin()+1,
              buffer.end()-1, s));
      double four_pi_eps0_i;
      
      _lineStream >> four_pi_eps0_i >> math::h_bar >> math::k_Boltzmann;
      math::eps0_i = four_pi_eps0_i * 4.0 * math::Pi;
      math::four_pi_eps_i = four_pi_eps0_i / param.nonbonded.epsilon;
      
      if (_lineStream.fail())
        io::messages.add("Bad line in PHYSICALCONSTANTS block",
                "InTopology", io::message::error);
    }
    else{
      io::messages.add("no PHYSICALCONSTANTS block in topology",
              "InTopology", io::message::error);
    }
    
  }

  if (param.system.npm){
    
    { // RESNAME
      if (!quiet)
	os << "\tRESNAME\n\t";
    
      DEBUG(10, "RESNAME block");
      buffer = m_block["RESNAME"];
      it = buffer.begin()+1;
      int n, num;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> num;
      ++it;
    
      if (!quiet && num > 10){
	for(n=0; n<10; ++n)
	  os << std::setw(8) << n+1;
	os << "\n\t";
      }
    
      for(n=0; it != buffer.end() - 1; ++it, ++n){
	std::string s;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> s;
      
	if (!quiet){
	  if (n && ((n) % 10) == 0) os << std::setw(10) << n << "\n\t";
	  os << std::setw(8) << s;      
	}

        if (s.length() > MAX_NAME) {
          std::ostringstream msg;
          msg << "Error in RESNAME block: residue name " << s 
              << " is too long (> " << MAX_NAME << " characters).";
          io::messages.add(msg.str(), "InTopology", io::message::error);
        }
      
	topo.residue_names().push_back(s);
      }

      if (n != num){
	io::messages.add("Error in RESNAME block: n!=num.",
			 "InTopology", io::message::error);
      }

      if (!quiet)
	os << "\n\tEND\n";
    
    
    } // RESNAME

    { // ATOMTYPENAME
      if (!quiet)
	os << "\tATOMTYPENAME\n\t";
    
      DEBUG(10, "ATOMTYPENAME block");
      buffer = m_block["ATOMTYPENAME"];
      it = buffer.begin()+1;
      int n, num;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> num;
      ++it;
    
      if (!quiet)
	os << "\t" << num << " atom types\n";
      
      for(n=0; it != buffer.end() - 1; ++it, ++n){
	std::string s;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> s;

        if (s.length() > MAX_NAME) {
          std::ostringstream msg;
          msg << "Error in ATOMTYPENAME block: type " << s 
              << " is too long (> " << MAX_NAME << " characters).";
          io::messages.add(msg.str(), "InTopology", io::message::error);
        }
      
	topo.atom_names()[s] = n;
      }

      if (n != num){
	io::messages.add("Error in ATOMTYPENAME block: n!=num.",
			 "InTopology", io::message::error);
      }

      if (!quiet)
	os << "\tEND\n";
    
    } // ATOMTYPENAME


    // os << "time after RESNAME: " << util::now() - start << std::endl;

    { // SOLUTEATOM
      DEBUG(10, "SOLUTEATOM block");
      buffer = m_block["SOLUTEATOM"];
  
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      topo.resize(num);

      if (!quiet)
	os << "\tSOLUTEATOM\n\t"
		  << "\tnumber of atoms : " << num;
    
      // put the rest of the block into a single stream
      ++it;

      // std::string soluteAtoms;
      // concatenate(it, buffer.end()-1, soluteAtoms);

      // _lineStream.clear();
      // _lineStream.str(soluteAtoms);

      // os << "\ntime after concatenate: " << util::now() - start << std::endl;

      int a_nr, r_nr, t, cg, n_ex, a_ex;
      double m, q;
      std::string s;
      std::set<int> ex;
      std::set<int> ex14;
      
      for(n=0; n < num; ++n){

	_lineStream.clear();
	_lineStream.str(*it);
	++it;
      
	_lineStream >> a_nr >> r_nr >> s >> t >> m >> q >> cg;

	if (a_nr != n+1){
	  io::messages.add("Error in SOLUTEATOM block: atom number not sequential.",
			   "InTopology", io::message::error);
	}

	if (r_nr > int(topo.residue_names().size()) || r_nr < 1){
	  io::messages.add("Error in SOLUTEATOM block: residue number out of range.",
			   "InTopology", io::message::error);
	}

        if (s.length() > MAX_NAME) {
          std::ostringstream msg;
          msg << "Error in SOLUTEATOM block: atom name " << s
              << " is too long (> " << MAX_NAME << " characters).";
          io::messages.add(msg.str(), "InTopology", io::message::error);
        }

      
	if (t < 1){
	  io::messages.add("Error in SOLUTEATOM block: iac < 1.",
			   "InTopology", io::message::error);
	}

	if (m <= 0){
	  io::messages.add("Error in SOLUTEATOM block: mass < 0.",
			   "InTopology", io::message::error);
	}

	if (cg != 0 && cg != 1){
	  io::messages.add("Error in SOLUTEATOM block: cg = 0 / 1.",
			   "InTopology", io::message::error);
	}

	if (!(_lineStream >> n_ex)){
	  if (_lineStream.eof()){
	    _lineStream.clear();
	    _lineStream.str(*it);
	    ++it;
	    _lineStream >> n_ex;
	  }
	  else{
	    io::messages.add("Error in SOLUTEATOM block: number of exclusions "
			     "could not be read.",
			     "InTopology", io::message::error);	
	  }
	}

	if (n_ex < 0){
	  io::messages.add("Error in SOLUTEATOM block: number of exclusions < 0.",
			   "InTopology", io::message::error);
	}
      
	// exclusions
	ex.clear();
	for(int i=0; i<n_ex; ++i){
	  if (!(_lineStream >> a_ex)){
	    if (_lineStream.eof()){
	      _lineStream.clear();
	      _lineStream.str(*it);
	      ++it;
	      _lineStream >> a_ex;
	    }
	    else{
	      io::messages.add("Error in SOLUTEATOM block: exclusion "
			       "could not be read.",
			       "InTopology", io::message::error);	
	    }
	  }

	  if (a_ex <= a_nr)
	    io::messages.add("Error in SOLUTEATOM block: exclusions only to "
			     "larger atom numbers.",
			     "InTopology", io::message::error);	

	  ex.insert(a_ex-1);
	}
      
	// 1,4 - pairs
	if (!(_lineStream >> n_ex)){
	  if (_lineStream.eof()){
	    _lineStream.clear();
	    _lineStream.str(*it);
	    ++it;
	    _lineStream >> n_ex;
	  }
	  else{
	    io::messages.add("Error in SOLUTEATOM block: number of 1,4 - exclusion "
			     "could not be read.",
			     "InTopology", io::message::error);	
	  }
	}

	if (n_ex < 0){
	  io::messages.add("Error in SOLUTEATOM block: number of 1,4 exclusions < 0.",
			   "InTopology", io::message::error);
	}

	ex14.clear();
	for(int i=0; i<n_ex; ++i){
	  if (!(_lineStream >> a_ex)){
	    if (_lineStream.eof()){
	      _lineStream.clear();
	      _lineStream.str(*it);
	      ++it;
	      _lineStream >> a_ex;
	    }
	    else{
	      io::messages.add("Error in SOLUTEATOM block: 1,4 - exclusion "
			       "could not be read.",
			       "InTopology", io::message::error);	
	    }
	  }

	  if (a_ex <= a_nr)
	    io::messages.add("Error in SOLUTEATOM block: 1,4 - exclusions only to "
			     "larger atom numbers.",
			     "InTopology", io::message::error);	
	
	  ex14.insert(a_ex-1);
	}
      
	if (_lineStream.fail())
	  io::messages.add("bad line in SOLUTEATOM block",
			   "In_Topology",
			   io::message::critical);

	topo.add_solute_atom(s, r_nr-1, t-1, m, q, cg, ex, ex14);
      }
      if (!quiet)
	os << "\n\tEND\n";
    
    } // SOLUTEATOM
    
    // os << "time after SOLUTEATOM: " << util::now() - start << std::endl;
    
    { // SOLUTEPOLARISATION
      DEBUG(10, "SOLUTEPOLARISATION block");

      if (!quiet)
	os << "\tSOLUTEPOLARISATION";
    
      buffer.clear();
      buffer = m_block["SOLUTEPOLARISATION"];
      if (buffer.size()){
	it = buffer.begin() + 1;
      
	_lineStream.clear();
	_lineStream.str(*it);
      
	int num, n;
	_lineStream >> num;
	++it;
        
	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i;
          double polarisability, coscharge, damping_level, damping_power;
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> polarisability >> coscharge >> damping_level
                      >> damping_power;
	
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in SOLUTEPOLARISATION block",
			     "In_Topology", io::message::error);
	  } 
          
          if (i > int(topo.num_solute_atoms()) || i < 1){
	    io::messages.add("Atom number out of range in SOLUTEPOLARISATION block",
			     "In_Topology", io::message::error);
          } else {
	    DEBUG(10, "\tpolarisable atom: " << i);
	    topo.polarisability()[i-1] = polarisability / math::four_pi_eps_i;
            topo.coscharge()[i-1] = coscharge;
            topo.damping_level()[i-1] = damping_level;
            topo.damping_power()[i-1] = damping_power;
            topo.is_polarisable()[i-1] = bool(polarisability > 0.0);
          }
	}
      
	if(n != num){
	  io::messages.add("Wrong number of polarisable atoms in SOLUTEPOLARISATION block",
			   "In_Topology", io::message::error);
	}
      }
    } // SOLUTEPOLARISATION
  
    { // BONDH
      DEBUG(10, "BONDH block");

      if (!quiet)
	os << "\tBOND";
    
      buffer.clear();
      buffer = m_block["BONDH"];
      if (buffer.size()){
	it = buffer.begin() + 1;
      
	_lineStream.clear();
	_lineStream.str(*it);
      
	int num, n;
	_lineStream >> num;
	++it;

        if (!quiet) {
          if (param.constraint.ntc == 2 || param.constraint.ntc == 3) {
            os << "\n\t\t"
                    << num
                    << " bonds from BONDH block added to CONSTRAINT";
          } else
            os << "\n\t\tbonds containing hydrogens : "
                  << num;
        }
      
	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i, j, t;
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> t;
	
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in BONDH block",
			     "In_Topology", io::message::error);
	  }
	
	  if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	      i < 1 || j < 1){
	    io::messages.add("Atom number out of range in BONDH block",
			     "In_Topology", io::message::error);
	  }
	
	  if (param.constraint.ntc == 2 || param.constraint.ntc == 3){
	    topo.solute().distance_constraints().
	      push_back(topology::two_body_term_struct(i-1, j-1, t-1));
	  }
	  else
	    topo.solute().bonds().
	      push_back(topology::two_body_term_struct(i-1, j-1, t-1));
	
	}
      
	if(n != num){
	  io::messages.add("Wrong number of bonds in BONDH block",
			   "In_Topology", io::message::error);
	}
      }

    } // BONDH
  
    { // BOND
      DEBUG(10, "BOND block");
      buffer = m_block["BOND"];

      if (buffer.size()){
      
	it = buffer.begin() + 1;
	_lineStream.clear();
	_lineStream.str(*it);
	int num, n;
	_lineStream >> num;
	++it;

	if (!quiet){
	  if (param.constraint.ntc == 3){
	    os << "\n\t\t"
		      << num
		      << " bonds from BOND block added to CONSTRAINT";
	  }
	  else
	    os << "\n\t\tbonds not containing hydrogens : "
		      << num;
	}
      
	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i, j, t;
	
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> t;
	
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in BOND block",
			     "In_Topology", io::message::error);
	  }
      
	  if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	      i < 1 || j < 1){
	    io::messages.add("Atom number out of range in BOND block",
			     "In_Topology", io::message::error);
	  }
      
	  if (param.constraint.ntc == 3){
	    topo.solute().distance_constraints().
	      push_back(topology::two_body_term_struct(i-1, j-1, t-1));
	  }
	  else
	    topo.solute().bonds().
	      push_back(topology::two_body_term_struct(i-1, j-1, t-1));
	}
    
	if(n != num){
	  io::messages.add("Wrong number of bonds in BOND block",
			   "In_Topology", io::message::error);
	}
      }
  
      if (!quiet)
	os << "\n\tEND\n";

    } // BOND

    { // check the bonds

      if (param.force.bond){
	if (!check_type(m_block[bond_bname], topo.solute().bonds())){
	  io::messages.add("Illegal bond type in BOND(H) block (type not in "
			   + bond_bname + ")",
			   "In_Topology", io::message::error);
	}
      }
    }
    
    { // CONSTRAINT
      DEBUG(10, "CONSTRAINT block");
      buffer = m_block["CONSTRAINT"];
  
      if (buffer.size() && param.constraint.ntc != 1){
      
	it = buffer.begin() + 1;
	_lineStream.clear();
	_lineStream.str(*it);
	int num, n;
	_lineStream >> num;
	++it;

	if (!quiet)
	  os << "\tCONSTRAINT\n\t\t"
		    << num
		    << " bonds in CONSTRAINT block."
		    << "\n\t\ttotal of constraint bonds : " 
		    << num + unsigned(topo.solute().distance_constraints().size())
		    << "\n\tEND\n";
      
	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i, j, t;
	
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> t;
	
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in CONSTRAINT block",
			     "In_Topology", io::message::error);
	  }
      
	  if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	      i < 1 || j < 1){
	    io::messages.add("Atom number out of range in CONSTRAINT block",
			     "In_Topology", io::message::error);
	  }
      
	  topo.solute().distance_constraints().
	    push_back(topology::two_body_term_struct(i-1, j-1, t-1));
	}
    
	if(n != num){
	  io::messages.add("Wrong number of bonds in CONSTRAINT block",
			   "In_Topology", io::message::error);
	}
      }
    
    } // CONSTRAINT

    // check the bonds in constraints
    if (!check_type(m_block[bond_bname], topo.solute().distance_constraints())){
      io::messages.add("Illegal bond type in " + bond_bname + " block",
		       "In_Topology", io::message::error);
    }

    { // BONDANGLEH

      if (!quiet)
	os << "\tBONDANGLE";

      DEBUG(10, "BONDANGLEH block");
      buffer.clear();
      buffer = m_block["BONDANGLEH"];
  
      if(buffer.size()){
      
	it = buffer.begin() + 1;

	_lineStream.clear();
	_lineStream.str(*it);

	int num, n;
	_lineStream >> num;
	++it;

	if (!quiet)
	  os << "\n\t\tbondangles containing hydrogens : " << num;

	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i, j, k, t;
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> k >> t;
      
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in BONDANGLEH block",
			     "In_Topology", io::message::error);
	  }
      
	  if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	      k > int(topo.num_solute_atoms()) ||
	      i < 1 || j < 1 || k < 1){
	    io::messages.add("Atom number out of range in BONDANGLEH block",
			     "In_Topology", io::message::error);
	  }

	  topo.solute().angles().
	    push_back(topology::three_body_term_struct(i-1, j-1, k-1, t-1));
	}
    
	if(n != num){
	  io::messages.add("Wrong number of bonds in BONDANGLEH block",
			   "In_Topology", io::message::error);
	}
      }
        
    } // BONDANGLEH
  
    { // BONDANGLE
      DEBUG(10, "BONDANGLE block");
      buffer = m_block["BONDANGLE"];
  
      if (buffer.size()){
      
	it = buffer.begin() + 1;
	_lineStream.clear();
	_lineStream.str(*it);
	int num, n;
	_lineStream >> num;
	++it;

	if (!quiet)
	  os << "\n\t\tbondangles not containing hydrogens : " << num;
    
	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i, j, k, t;
      
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> k >> t;
      
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in BONDANGLE block",
			     "In_Topology", io::message::error);
	  }
      
	  if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	      k > int(topo.num_solute_atoms()) ||
	      i < 1 || j < 1 || k < 1){
	    io::messages.add("Atom number out of range in BONDANGLE block",
			     "In_Topology", io::message::error);
	  }
      
	  topo.solute().angles().
	    push_back(topology::three_body_term_struct(i-1, j-1, k-1, t-1));
	}
    
	if(n != num){
	  io::messages.add("Wrong number of bonds in BONDANGLE block",
			   "In_Topology", io::message::error);
	}
      }

      if (!quiet)
	os << "\n\tEND\n";

    } // BONDANGLE

    // os << "time after BONDANGLE: " << util::now() - start << std::endl;

    // check the angles
    if (!check_type(m_block[angle_bname], topo.solute().angles())){
      io::messages.add("Illegal bond angle type in BONDANGLE(H) block",
		       "In_Topology", io::message::error);
    }

    { // IMPDIHEDRAL
      DEBUG(10, "IMPDIHEDRAL block");
      buffer = m_block["IMPDIHEDRAL"];
  
      if (!quiet)
	os << "\tIMPDIHEDRAL";

      if(buffer.size()){
      
	it = buffer.begin() + 1;
	_lineStream.clear();
	_lineStream.str(*it);
	int num, n;
	_lineStream >> num;
	++it;

	if (!quiet)
	  os << "\n\t\timproper dihedrals not containing hydrogens : "
		    << num;
    
	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i, j, k, l, t;
      
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> k >> l >> t;
      
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in IMPDIHEDRAL block",
			     "In_Topology", io::message::error);
	  }
      
	  if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	      k > int(topo.num_solute_atoms()) || l > int(topo.num_solute_atoms()) ||
	      i < 1 || j < 1 || k < 1 || l < 1){
	    io::messages.add("Atom number out of range in IMPDIHEDRAL block",
			     "In_Topology", io::message::error);
	  }
      
	  topo.solute().improper_dihedrals().
	    push_back(topology::four_body_term_struct(i-1, j-1, k-1, l-1, t-1));
	}
    
	if(n != num){
	  io::messages.add("Wrong number of bonds in IMPDIHEDRAL block",
			   "In_Topology", io::message::error);
	}
      }
    
    } // IMPDIHEDRAL

    { // IMPDIHEDRALH
      DEBUG(10, "IMPDIHEDRALH block");
      buffer.clear();
      buffer = m_block["IMPDIHEDRALH"];
  
      if(buffer.size()){
      
	it = buffer.begin() + 1;

	_lineStream.clear();
	_lineStream.str(*it);

	int num, n;
	_lineStream >> num;
	++it;

	if (!quiet)
	  os << "\n\t\timproper dihedrals containing hydrogens : "
		    << num;

	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i, j, k, l, t;
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> k >> l >> t;
      
      
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in IMPDIHEDRALH block",
			     "In_Topology", io::message::error);
	  }
      
	  if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	      k > int(topo.num_solute_atoms()) || l > int(topo.num_solute_atoms()) ||
	      i < 1 || j < 1 || k < 1 || l < 1){
	    io::messages.add("Atom number out of range in IMPDIHEDRALH block",
			     "In_Topology", io::message::error);
	  }
      
	  topo.solute().improper_dihedrals().
	    push_back(topology::four_body_term_struct(i-1, j-1, k-1, l-1, t-1));
	}
    
	if(n != num){
	  io::messages.add("Wrong number of bonds in IMPDIHEDRALH block",
			   "In_Topology", io::message::error);
	}
      }

      if (!quiet)
	os << "\n\tEND\n";
    
    } // IMPDIHEDRALH

    // check the imporopers
    if (!check_type(m_block["IMPDIHEDRALTYPE"], topo.solute().improper_dihedrals())){
      io::messages.add("Illegal improper dihedral type in IMPDIHEDRAL(H) block",
		       "In_Topology", io::message::error);
    }  

    { // DIHEDRAL
      DEBUG(10, "DIHEDRAL block");    
      buffer = m_block["DIHEDRAL"];

      if (!quiet)
	os << "\tDIHEDRAL";
    
      if(buffer.size()){
      
	it = buffer.begin() + 1;
	_lineStream.clear();
	_lineStream.str(*it);
	int num, n;
	_lineStream >> num;
	++it;
    
	if (!quiet)
	  os << "\n\t\tdihedrals not containing hydrogens : "
		    << num;

	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i, j, k, l, t;
      
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> k >> l >> t;
      
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in DIHEDRAL block",
			     "In_Topology", io::message::error);
	  }
      
	  if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	      k > int(topo.num_solute_atoms()) || l > int(topo.num_solute_atoms()) ||
	      i < 1 || j < 1 || k < 1 || l < 1){
	    io::messages.add("Atom number out of range in DIHEDRAL block",
			     "In_Topology", io::message::error);
	  }
      
	  topo.solute().dihedrals().
	    push_back(topology::four_body_term_struct(i-1, j-1, k-1, l-1, t-1));
	}
    
	if(n != num){
	  io::messages.add("Wrong number of bonds in DIHEDRAL block",
			   "In_Topology", io::message::error);
	}
      }
    
    } // DIHEDRAL

    { // DIHEDRALH
      DEBUG(10, "DIHEDRALH block");
      buffer.clear();
      buffer = m_block["DIHEDRALH"];
      if(buffer.size()){
      
	it = buffer.begin() + 1;

	_lineStream.clear();
	_lineStream.str(*it);

	int num, n;
	_lineStream >> num;
	++it;

	if (!quiet)
	  os << "\n\t\tdihedrals containing hydrogens : "
		    << num;

	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  int i, j, k, l, t;
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> k >> l >> t;
      
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in DIHEDRALH block",
			     "In_Topology", io::message::error);
	  }
      
	  if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	      k > int(topo.num_solute_atoms()) || l > int(topo.num_solute_atoms()) ||
	      i < 1 || j < 1 || k < 1 || l < 1){
	    io::messages.add("Atom number out of range in DIHEDRALH block",
			     "In_Topology", io::message::error);
	  }
      
	  topo.solute().dihedrals().
	    push_back(topology::four_body_term_struct(i-1, j-1, k-1, l-1, t-1));
	}
    
	if(n != num){
	  io::messages.add("Wrong number of bonds in DIHEDRALH block",
			   "In_Topology", io::message::error);
	}
      }
    
      if (!quiet)
	os << "\n\tEND\n";
    
    } // DIHEDRALH

    // check the dihedrals
    if (!check_type(m_block[dihedral_bname], topo.solute().dihedrals())){
      io::messages.add("Illegal dihedral type in DIHEDRAL(H) block",
		       "In_Topology", io::message::error);
    }

    { // VIRTUALGRAIN
      DEBUG(10, "VIRTUALGRAIN block");
      buffer.clear();
      buffer = m_block["VIRTUALGRAIN"];
      if(buffer.size()){
      
	it = buffer.begin() + 1;

	_lineStream.clear();
	_lineStream.str(*it);

	int num, n;
	_lineStream >> num;
	++it;

	if (!quiet)
	  os << "\tVIRTUALGRAIN\n\t\tVirtual Grains : "
	     << num;

	for(n=0; it != buffer.end() - 1; ++it, ++n){
	  _lineStream.clear();
	  _lineStream.str(*it);

	  // number of real atoms to define virtual atom
	  int index, i, q;
	  _lineStream >> index >> i;

	  std::vector<int> cog;
	  
	  for(int j=0; j<i; ++j){
	    _lineStream >> q;
	    cog.push_back(q-1);
	  }
      
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in VIRTUALGRAIN block",
			     "In_Topology", io::message::error);
	  }
            
	  topo.virtual_grains().
	    push_back(topology::virtual_grain_struct(index-1, 
						     util::Virtual_Atom(util::va_cog, cog)));
	}
    
	if(n != num){
	  io::messages.add("Wrong number of elements in VIRTUALGRAIN block",
			   "In_Topology", io::message::error);
	}

	if (!quiet)
	  os << "\n\tEND\n";
      }
    
    } // VIRTUALGRAIN

    // add the solute molecules (should be done before solvate ;-)
    
    { // SOLUTEMOLECULES
      DEBUG(10, "read SOLUTEMOLECULES");

      buffer.clear();
      std::string s;
  
      buffer = m_block["SOLUTEMOLECULES"];
  
      if (!buffer.size()){
        io::messages.add("no SOLUTEMOLECULES block in topology",
                         "In_Topology", io::message::error);
      } else {
        _lineStream.clear();
        _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
        int num;
  
        _lineStream >> num;

        if(num<0){
          io::messages.add("negative number of SOLUTEMOLECULES is not allowed",
		     "In_Topology", io::message::error);
        }
  
        unsigned int m;
        unsigned int old_m=0;

        for(int i=0; i<num; ++i){
          _lineStream >> m;
          topo.molecules().push_back(m);
          DEBUG(11, "add submol " << m);
          if(m<old_m){
            io::messages.add("wrong order in SOLUTEMOLECULES block",
                             "In_Topology", io::message::error);
            break;
          }
          old_m=m;
        }
  
        if (_lineStream.fail())
          io::messages.add("bad line in SOLUTEMOLECULES block",
                           "In_Topology", io::message::error);
      }
    } // SOLUTEMOLECULES

    if (topo.molecules().size() == 0){
      topo.molecules().push_back(0);
      topo.molecules().push_back(topo.num_solute_atoms());
    }

    // solutemolecules check
    if (topo.molecules().back()
	!= topo.num_solute_atoms()){
    
      std::cout << "ERROR: Solute molecules wrong\n"
		<< "\tlast atom in SOLUTEMOLECULES block = " << topo.molecules().back()
		<< "\n\tlast atom in topology = " << topo.num_solute_atoms()
		<< std::endl;
      
      io::messages.add("Error in SOLUTEMOLECULES block: "
		       "last solute molecule has to end with last solute atom",
		       "In_Topology", io::message::error);
    }
    
    topo.num_solute_molecules() = topo.molecules().size() - 1;

    { // TEMPERATUREGROUPS
      DEBUG(10, "read TEMPERATUREGROUPS");

      buffer.clear();
      std::string s;
  
      buffer = m_block["TEMPERATUREGROUPS"];
  
      if (!buffer.size()){
        io::messages.add("no TEMPERATUREGROUPS block in topology",
                         "In_Topology", io::message::error);
      } else {
        _lineStream.clear();
        _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
        int num;
  
        _lineStream >> num;

        if(num<0){
          io::messages.add("negative number of TEMPERATUREGROUPS is not allowed",
		     "In_Topology", io::message::error);
        }
  
        unsigned int m;
        unsigned int old_m=0;

        for(int i=0; i<num; ++i){
          _lineStream >> m;
          topo.temperature_groups().push_back(m);
          DEBUG(11, "add temperature group " << m);
          if(m<old_m){
            io::messages.add("wrong order in TEMPERATUREGROUPS block",
                             "In_Topology", io::message::error);
            break;
          }
          old_m=m;
        }
  
        if (_lineStream.fail())
          io::messages.add("bad line in TEMPERATUREGROUPS block",
                           "In_Topology", io::message::error);
      }
    } // TEMPERATUREGROUPS

    if (topo.temperature_groups().size() == 0){
      topo.temperature_groups().push_back(0);
      topo.temperature_groups().push_back(topo.num_solute_atoms());
    }

    // temperature groups check
    if (topo.temperature_groups().back()
	!= topo.num_solute_atoms()){
    
      std::cout << "ERROR: temperature groups wrong\n"
		<< "\tlast atom in TEMPERATUREGROUPS block = " << topo.temperature_groups().back()
		<< "\n\tlast atom in topology = " << topo.num_solute_atoms()
		<< std::endl;
      
      io::messages.add("Error in TEMPERATUREGROUPS block: "
		       "last temperature group has to end with last solute atom",
		       "In_Topology", io::message::error);
    }
    
    topo.num_solute_temperature_groups() = topo.temperature_groups().size() - 1;
    
    { // PRESSUREGROUPS
      DEBUG(10, "read PRESSUREGROUPS");

      buffer.clear();
      std::string s;
  
      buffer = m_block["PRESSUREGROUPS"];
  
      if (!buffer.size()){
        io::messages.add("no PRESSUREGROUPS block in topology",
                         "In_Topology", io::message::error);
      } else {
        _lineStream.clear();
        _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
        int num;
  
        _lineStream >> num;

        if(num<0){
          io::messages.add("negative number of PRESSUREGROUPS is not allowed",
		     "In_Topology", io::message::error);
        }
  
        unsigned int m;
        unsigned int old_m=0;

        for(int i=0; i<num; ++i){
          _lineStream >> m;
          topo.pressure_groups().push_back(m);
          DEBUG(11, "add pressure group " << m);
          if(m<old_m){
            io::messages.add("wrong order in PRESSUREGROUPS block",
                             "In_Topology", io::message::error);
            break;
          }
          old_m=m;
        }
  
        if (_lineStream.fail())
          io::messages.add("bad line in PRESSUREGROUPS block",
                           "In_Topology", io::message::error);
      }
    } // PRESSUREGROUPS

    if (topo.pressure_groups().size() == 0){
      topo.pressure_groups().push_back(0);
      topo.pressure_groups().push_back(topo.num_solute_atoms());
    }

    // pressure groups check
    if (topo.pressure_groups().back()
	!= topo.num_solute_atoms()){
    
      std::cout << "ERROR: pressure groups wrong\n"
		<< "\tlast atom in PRESSUREGROUPS block = " << topo.pressure_groups().back()
		<< "\n\tlast atom in topology = " << topo.num_solute_atoms()
		<< std::endl;
      
      io::messages.add("Error in PRESSUREGROUPS block: "
		       "last pressure group has to end with last solute atom",
		       "In_Topology", io::message::error);
    }
    
    topo.num_solute_pressure_groups() = topo.pressure_groups().size() - 1;
     
    { // PATHINTSPEC
      buffer = m_block["PATHINTSPEC"];
      if (buffer.size())
        io::messages.add("md++ does not support path-integral simulations (PATHINTSPEC block).",
                "InTopology", io::message::warning);
     }
  } // npm != 0
    
  { // SOLVENTATOM and SOLVENTCONSTR
    // give it a number (SOLVENTATOM1, SOLVENTATOM2) for multiple
    // solvents...
    DEBUG(10, "SOLVENTATOM block");
    buffer = m_block["SOLVENTATOM"];

    if (!quiet)
      os << "\tSOLVENT";
    
    if (buffer.size()){
      
      unsigned int res_nr = unsigned(topo.residue_names().size());
    
      topo.residue_names().push_back("SOLV");

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
    
      if (!quiet)
	os << "\n\t\tatoms : " << num;

      topology::Solvent s;
    
      std::string name;
      int i, iac;
      double mass, charge;
    
      for(n=0; it != buffer.end()-1; ++it, ++n){
	_lineStream.clear();
	_lineStream.str(*it);
      
	_lineStream >> i >> name >> iac >> mass >> charge;

	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in SOLVENTATOM block",
			   "In_Topology", io::message::error);	
	}

        if (name.length() > MAX_NAME) {
          std::ostringstream msg;
          msg << "Error in SOLVENTATOM block: name " << name 
              << " is too long (> " << MAX_NAME << " characters).";
          io::messages.add(msg.str(), "InTopology", io::message::error);
        }

	s.add_atom(name, res_nr, iac-1, mass, charge);
      }
    
      if (n!=num){
	io::messages.add("Error in SOLVENTATOM block (num != n)",
			 "In_Topology", io::message::error);
      }

      // end SOLVENTATOM
      
      DEBUG(10, "SOLVENTEPOLARISATION block");
      
      if (!quiet)
        os << "\tSOLVENTEPOLARISATION";
      
      buffer.clear();
      buffer = m_block["SOLVENTPOLARISATION"];
      if (buffer.size()){
        it = buffer.begin() + 1;
        
        _lineStream.clear();
        _lineStream.str(*it);
        
        int num, n;
        _lineStream >> num;
        ++it;
        
        for(n=0; it != buffer.end() - 1; ++it, ++n){
          int i;
          double polarisability, coscharge, damping_level, damping_power;
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> polarisability >> coscharge >> damping_level
                  >> damping_power;
          
          if (_lineStream.fail() || ! _lineStream.eof()){
            io::messages.add("Bad line in SOLVENTPOLARISATION block",
                    "In_Topology", io::message::error);
          }
          
          if (i > int(s.num_atoms()) || i < 1){
            io::messages.add("Atom number out of range in SOLVENTPOLARISATION block",
                    "In_Topology", io::message::error);
          } else {  
            DEBUG(10, "\tpolarisable atom: " << i);
          
            s.atoms()[i-1].polarisability = polarisability / math::four_pi_eps_i;
            s.atoms()[i-1].coscharge = coscharge;
            s.atoms()[i-1].damping_level = damping_level;
            s.atoms()[i-1].damping_power = damping_power;
          }
        }
        
        if(n != num){
          io::messages.add("Wrong number of polarisable atoms in SOLVENTPOLARISATION block",
                  "In_Topology", io::message::error);
        }
      }
      // solvent atoms have been read into s
    
      //--------------------------------------------------
      // lookup the number of bond types
      // add additional ones for the solvent constraints
      {
	std::vector<interaction::bond_type_struct> b;
	std::ostringstream os;
	read_harmonic_bonds(b, os);
      }

      //--------------------------------------------------
    
      // now read the solvent constraints
      buffer = m_block["SOLVENTCONSTR"];
      if (!buffer.size()){
	io::messages.add("no SOLVENTCONST (block missing).",
			 "In_Topology", io::message::notice);	
      }
      else{

	it = buffer.begin() + 1;
	_lineStream.clear();
	_lineStream.str(*it);

	_lineStream >> num;
	++it;
    
	if (!quiet)
	  os << "\n\t\tconstraints : " << num;
	
	int j;
	double b0;
	
	for(n=0; it != buffer.end()-1; ++it, ++n){
	  _lineStream.clear();
	  _lineStream.str(*it);
	  
	  _lineStream >> i >> j >> b0;
	  
	  if (_lineStream.fail() || ! _lineStream.eof()){
	    io::messages.add("Bad line in SOLVENTCONSTR block",
			     "In_Topology", io::message::error);	
	  }
	  
	  // the solvent (distance constraints) bond types
	  s.add_distance_constraint
	    (topology::two_body_term_struct(i-1, j-1, num_solute_bondtypes + n));
	}

	if (n!=num){
	  io::messages.add("Error in SOLVENTCONSTR block (num != n)",
			   "In_Topology", io::message::error);
	}
      }
      topo.add_solvent(s);
    }
    else{
      io::messages.add("no solvent topology specified",
		       "In_Topology", io::message::warning);
    }
  } // SOLVENTCONSTR

  // add the solvent to the topology
  if (!quiet)
    os << "\n\t\tadding " << param.system.nsm 
	      << " solvents.";
  
  // if (param.system.nsm) 
  topo.solvate(0, param.system.nsm);  
  
  if (!quiet)
    os << "\n\tEND\n";

  // set lambda (new and old one, yes it looks strange...)
  topo.lambda(param.perturbation.lambda);
  topo.lambda(param.perturbation.lambda);
  
  topo.lambda_exp(param.perturbation.lambda_exponent);

  //==================================================
  // CHECKING
  //==================================================
    
  // solute molecule check
  if (topo.molecules().back()
      != topo.num_atoms()){
    
    io::messages.add("Error in SOLUTEMOLECULE / solvation block: "
		     "last solute molecule has to end with last atom",
		     "In_Topology", io::message::error);
  }
  
  // temperature group check
  if (topo.temperature_groups().back()
      != topo.num_atoms()){
    
    io::messages.add("Error in TEMPERATUREGROUPS / solvation block: "
		     "last temperature group has to end with last atom",
		     "In_Topology", io::message::error);
  }  
  
  // pressure group check
  if (topo.pressure_groups().back()
      != topo.num_atoms()){
    
    io::messages.add("Error in PRESSUREGROUPS / solvation block: "
		     "last pressure group has to end with last atom",
		     "In_Topology", io::message::error);
  }  

  // chargegroup check (starts with 0)
  if (topo.chargegroups()[topo.num_solute_chargegroups()] != int(topo.num_solute_atoms())){
    io::messages.add("Error: last solute atom has to be end of chargegroup",
		     "In_Topology",
		     io::message::error);
    os << "ERROR:"
	      << "\tsolute cg    : " << topo.num_solute_chargegroups() << "\n"
	      << "\tsolute atoms : " << topo.num_solute_atoms() << "\n"
	      << "\tlast cg      : " << topo.chargegroups()[topo.num_solute_chargegroups()] << "\n";
  }

  if (!quiet)
    os << "\n\tSOLUTE [sub]molecules: " 
	      << unsigned(topo.molecules().size()) - param.system.nsm - 1 << "\n";

  DEBUG(10, "molecules().size: " << unsigned(topo.molecules().size())
	<< " nsm : " << param.system.nsm);
  
  if (!quiet)
    os << "\n\tSOLUTE temperature groups: " 
	      << unsigned(topo.temperature_groups().size()) - param.system.nsm - 1 << "\n";

  DEBUG(10, "temperature_groups().size: " << unsigned(topo.temperature_groups().size()));
  
  if (!quiet)
    os << "\n\tSOLUTE pressure groups: " 
	      << unsigned(topo.pressure_groups().size()) - param.system.nsm - 1 << "\n";

  DEBUG(10, "pressure_groups().size: " << unsigned(topo.pressure_groups().size()));

  // energy group check
  if (param.force.energy_group.size() == 0){
    param.force.energy_group.push_back(topo.num_atoms() -1);
  }

  if (param.force.energy_group.back() > topo.num_atoms()-1){
    io::messages.add("Error in FORCE block: "
		     "last energy group has to end with last atom",
		     "In_Topology", io::message::error);
  }
  else if (param.force.energy_group.back() < topo.num_atoms()-1){
    param.force.energy_group.push_back(topo.num_atoms() -1);
    io::messages.add("FORCE block: "
            "added an additional energy group",
            "In_Topology", io::message::warning);
  }

  // and add them
  DEBUG(10, "adding energy groups : " << param.force.energy_group.size());
  unsigned int atom = 0;
  for(unsigned int i=0; i<param.force.energy_group.size(); ++i){
    assert(param.force.energy_group.size() > i);
    topo.energy_groups().push_back(param.force.energy_group[i]);
    DEBUG(10, "energy group " << i << " start = " << atom << " end = " << param.force.energy_group[i]);
    for( ; atom <= param.force.energy_group[i]; ++atom){
      topo.atom_energy_group().push_back(i);
      // DEBUG(11, "atom " << atom << ": " << i);
    }
  }

  // Now that we have the energy groups, we initialize the 
  // LAMBDAS parameters that depend on them. 
  int maxnilg = param.force.energy_group.size();
  std::vector< double > one(maxnilg, 1.0);
  std::vector< double > zero(maxnilg, 0.0);
  for (unsigned int i = 0; i < param.lambdas.a.size(); i++) {
    // check whether we have to resize the lambdas arrays
    int lam_size = param.lambdas.a[i].size();
    if (lam_size < maxnilg) {
      param.lambdas.a[i].insert(param.lambdas.a[i].end(), maxnilg - lam_size, zero);
      lam_size = param.lambdas.b[i].size();
      param.lambdas.b[i].insert(param.lambdas.b[i].end(), maxnilg - lam_size, zero);
      lam_size = param.lambdas.c[i].size();
      param.lambdas.c[i].insert(param.lambdas.c[i].end(), maxnilg - lam_size, zero);
      lam_size = param.lambdas.d[i].size();
      param.lambdas.d[i].insert(param.lambdas.d[i].end(), maxnilg - lam_size, one);
      lam_size = param.lambdas.e[i].size();
      param.lambdas.e[i].insert(param.lambdas.e[i].end(), maxnilg - lam_size, zero);
      for (unsigned int j = 0; j < param.lambdas.a[i].size(); j++) {
        lam_size = param.lambdas.a[i][j].size();
        param.lambdas.a[i][j].insert(param.lambdas.a[i][j].end(),maxnilg - lam_size,0.0);
        lam_size = param.lambdas.b[i][j].size();
        param.lambdas.b[i][j].insert(param.lambdas.b[i][j].end(),maxnilg - lam_size,0.0);
        lam_size = param.lambdas.c[i][j].size();
        param.lambdas.c[i][j].insert(param.lambdas.c[i][j].end(),maxnilg - lam_size,0.0);
        lam_size = param.lambdas.d[i][j].size();
        param.lambdas.d[i][j].insert(param.lambdas.d[i][j].end(),maxnilg - lam_size,1.0);
        lam_size = param.lambdas.e[i][j].size();
        param.lambdas.e[i][j].insert(param.lambdas.e[i][j].end(),maxnilg - lam_size,0.0);
      }
    } // do resize lambdas arrays
  }

  DEBUG(10, "multibath?");
  if(!param.multibath.found_multibath && param.multibath.found_tcouple){
    if (!quiet)
      os << "\tparsing a (deprecated) TCOUPLE block into the new "
		<< "MULTIBATH format.\n";
    
    util::parse_TCOUPLE(param, topo);
  }

  if (!quiet)
    os << "END\n";

  DEBUG(10, "topology read");
}

void io::In_Topology
::read_harmonic_bonds(std::vector<interaction::bond_type_struct> &b,
		      std::ostream & os)
{
  DEBUG(8, "read_harmonic_bonds");
  DEBUG(10, "(HARM)BONDTYPE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["HARMBONDTYPE"];
  if (buffer.size()){
    DEBUG(7, "reading in a DIRK (HARMBONDTYPE) block)");
    io::messages.add("harmonic bond force constants from HARMBONDTYPE block",
		     "In_Topology::bondtype", io::message::notice);
    
    int num, n=0;
    it = buffer.begin()+1;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num;
    ++it;
    for(; it!=buffer.end()-1; ++it, ++n){
      double k, r;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> k >> r;

      if (_lineStream.fail() || ! _lineStream.eof()){
	io::messages.add("bad line in HARMBONDTYPE block",
			 "In_Topology",
			 io::message::error);
	k = 0;
	r = 0;
      }

      // and add...
      b.push_back(interaction::bond_type_struct(k, r));
    }
  
    if (num != n)
      io::messages.add("not enough bond types in HARMBONDTYPE block",
		       "In_Topology",
		       io::message::error);
  }
  else if (m_block["BONDTYPE"].size()){
    buffer = m_block["BONDTYPE"];
    _lineStream.clear();
    _lineStream.str(buffer[1]);
    unsigned int num_types;
    _lineStream >> num_types;
    if (_lineStream.fail()) {
      io::messages.add("bad line in BONDSTRETCHTYPE block: number of types",
          "In_Topology",  io::message::error);
    }
    
    // 1. BONDTYPE 2. number of types
    for (it = buffer.begin() + 2; 
	 it != buffer.end() - 1; ++it) {

      double k, r;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> r;
      
      if (_lineStream.fail()){
	os << *it << std::endl;
	io::messages.add("bad line in BONDTYPE block!",
			 "In_Topology",
			 io::message::error);
      }
      if (! _lineStream.eof()){
	os << *it << std::endl;
	io::messages.add("eof not reached in BONDTYPE block",
			 "InTopology", io::message::warning);
      }

      // we are reading into harmonic bond term, so convert k
      k *= 2 * r * r;
      DEBUG(10, "\t\tbond type k: " << k << ", r: " << r << "\n");
      // and add...
      b.push_back(interaction::bond_type_struct(k, r));
    }
    if (b.size() != num_types) {
      io::messages.add("BONDTYPE block: number of types does not "
          "correspond with number of lines",
          "In_Topology",  io::message::error);
    }
  }
  else if (m_block["BONDSTRETCHTYPE"].size()){
    // read in the new block
    buffer = m_block["BONDSTRETCHTYPE"];
    _lineStream.clear();
    _lineStream.str(buffer[1]);
    unsigned int num_types;
    _lineStream >> num_types;
    if (_lineStream.fail()) {
      io::messages.add("bad line in BONDSTRETCHTYPE block: number of types",
          "In_Topology",  io::message::error);
    }
    
    // we are reading into harmonic bond term, so only read the harmonic force constant (kh)
    // 1. BONDTYPE 2. number of types
    for (it = buffer.begin() + 2; 
	 it != buffer.end() - 1; ++it) {
      double k, kh, r;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> kh >> r;
      
      if (_lineStream.fail()){
	os << *it << std::endl;
	io::messages.add("bad line in BONDSTRETCHTYPE block!",
			 "In_Topology",
			 io::message::error);
      }
      if (! _lineStream.eof()){
	os << *it << std::endl;
	io::messages.add("eof not reached in BONDSTRETCHTYPE block",
			 "InTopology", io::message::warning);
      }
      // check for consistency, allow 0.01% error
      double calc_kh = k*2.0*r*r;
      if (fabs(kh-calc_kh) / kh > 1.0e-4){
        std::ostringstream msg;
        msg << "harmonic and quartic force constant do not match (CHB!=CB*2*B0*B0): " << std::endl
            << "CHB = " << kh << ", CB*2*B0*B0 = " << calc_kh
                << " |CHB-CB*2*B0*B0| = " << fabs(kh-calc_kh) << std::endl;
        io::messages.add(msg.str(), "InTopology", io::message::warning);
      }
      // and add r and the harmonic force constant
      b.push_back(interaction::bond_type_struct(kh, r));
    }
    
    if (b.size() != num_types) {
      io::messages.add("BONDSTRETCHTYPE block: number of types does not "
          "correspond with number of lines",
          "In_Topology",  io::message::error);
    }
      
  }
  else{
    io::messages.add("either BONDTYPE, BONDSTRETCHTYPE or HARMBONDTYPE block must be present!",
		       "In_Topology",
		       io::message::error);
  }

  // also add the solvent constraints to the bond types...
  // (if there is one)
  
  num_solute_bondtypes = b.size();

  buffer = m_block["SOLVENTCONSTR"];
  if (buffer.size()){

    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);

    int num;
    _lineStream >> num;
    ++it;
    
    int i, j, n;
    double b0;
    
    for(n=0; it != buffer.end()-1; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> i >> j >> b0;
      
      if (_lineStream.fail() || ! _lineStream.eof()){
	io::messages.add("Bad line in SOLVENTCONSTR block",
			 "In_Topology", io::message::error);	
      }
      
      // the solvent (distance constraints) bond types
      b.push_back(interaction::bond_type_struct(0, b0));
      // (K is set to 0.0)
    }

    if (n!=num){
      io::messages.add("Error in SOLVENTCONSTR block (num != n)",
		       "In_Topology", io::message::error);
      
    }
  }
  
}

void io::In_Topology
::read_g96_bonds(std::vector<interaction::bond_type_struct> &b,
		 std::ostream & os)
{
  DEBUG(8, "read_g96_bonds");
  DEBUG(10, "BONDTYPE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  if (m_block["BONDTYPE"].size()){
    buffer = m_block["BONDTYPE"];

    // 1. BONDTYPE 2. number of types
    for (it = buffer.begin() + 2; 
	 it != buffer.end() - 1; ++it) {

      double k, r;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> r;
      
      if (_lineStream.fail()){
	os << *it << std::endl;
	io::messages.add("bad line in BONDTYPE block!",
			 "In_Topology",
			 io::message::error);
      }
      if (! _lineStream.eof()){
	os << *it << std::endl;
	io::messages.add("eof not reached in BONDTYPE block",
			 "InTopology", io::message::warning);
      }

      // and add...
      b.push_back(interaction::bond_type_struct(k, r));
    }
  }
  else if (m_block["BONDSTRETCHTYPE"].size()){
    // read in the new block
    buffer = m_block["BONDSTRETCHTYPE"];
    // we are reading into quartic bond term, so only read the quartic force constant (k)
    // 1. BONDTYPE 2. number of types
    for (it = buffer.begin() + 2; 
	 it != buffer.end() - 1; ++it) {
      double k, kh, r;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> kh >> r;
      
      if (_lineStream.fail()){
	os << *it << std::endl;
	io::messages.add("bad line in BONDSTRETCHTYPE block!",
			 "In_Topology",
			 io::message::error);
      }
      if (! _lineStream.eof()){
	os << *it << std::endl;
	io::messages.add("eof not reached in BONDSTRETCHTYPE block",
			 "InTopology", io::message::warning);
      }
      // check for consistency
      if (kh!=k*2*r*r){
        os << *it << std::endl;
        io::messages.add("harmonic and quartic force constant do not match (CHB!=CB*2*B0*B0)",
			 "InTopology", io::message::warning);
      }
      // and add r and the quartic force constant
      b.push_back(interaction::bond_type_struct(k, r));
    }
  }
  else{
    io::messages.add("either BONDTYPE or BONDSTRETCHTYPE block must be present!",
		       "In_Topology",
		       io::message::error);
  }

  num_solute_bondtypes = b.size();

  // also add the solvent constraints to the bond types...
  // (if there is one)
  buffer = m_block["SOLVENTCONSTR"];
  if (buffer.size()){
    
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    
    int num;
    _lineStream >> num;
    ++it;
    
    int i,j, n;
    double b0;
    
    for(n=0; it != buffer.end()-1; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> i >> j >> b0;
      
      if (_lineStream.fail() || ! _lineStream.eof()){
	io::messages.add("Bad line in SOLVENTCONSTR block",
			 "In_Topology", io::message::error);	
      }
      
      // the solvent (distance constraints) bond types
      // (K is set to 0.0)
      b.push_back(interaction::bond_type_struct(0, b0));
    }

    if (n!=num){
      io::messages.add("Error in SOLVENTCONSTR block (num != n)",
		       "In_Topology", io::message::error);
    }
  }
  
}

void io::In_Topology
        ::read_angles(std::vector<interaction::angle_type_struct> &b,
        std::ostream & os) {
    
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  if(m_block["BONDANGLEBENDTYPE"].size()){
    DEBUG(10, "BONDANGLEBENDTYPE block");
    buffer = m_block["BONDANGLEBENDTYPE"];
    // 1. BONDTYPE 2. number of types
    for (it = buffer.begin() + 2;
    it != buffer.end() - 1; ++it) {
      
      double k, kh, cos0;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> kh >> cos0;
      
      if (_lineStream.fail()){
        os << *it << std::endl;
        io::messages.add("bad line in BONDANGLEBENDTYPE block", "In_Topology",
                io::message::error);
      }
      if (! _lineStream.eof()){
        os << *it << std::endl;
        io::messages.add("eof not reached in BONDANGLEBENDTYPE block",
                "InTopology", io::message::warning);
      }
      
      // and add (force constant based on a potential harmonic in the angle cosine)
      b.push_back(interaction::angle_type_struct(k, cos(cos0 * 2 * math::Pi / 360.0)));
    }
  }
  else if(m_block["BONDANGLETYPE"].size()){
    DEBUG(10, "BONDANGLETYPE block");
    buffer = m_block["BONDANGLETYPE"];
      
    // 1. BONDTYPE 2. number of types
    for (it = buffer.begin() + 2;
    it != buffer.end() - 1; ++it) {
      
      double k, cos0;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> cos0;
      
      if (_lineStream.fail()){
        os << *it << std::endl;
        io::messages.add("bad line in BONDANGLETYPE block", "In_Topology",
                io::message::error);
      }
      if (! _lineStream.eof()){
        os << *it << std::endl;
        io::messages.add("eof not reached in BONDANGLETYPE block",
                "InTopology", io::message::warning);
      }
      
      // and add...
      b.push_back(interaction::angle_type_struct(k, cos(cos0 * 2 * math::Pi / 360.0)));
    }
  }
  else{
    io::messages.add("either BONDANGLEBENDTYPE or BONDANGLETYPE block must be present",
            "In_Topology", io::message::error);
  }
  
}

void io::In_Topology
::read_harm_angles(std::vector<interaction::angle_type_struct> &b,
		   std::ostream & os)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  if(m_block["BONDANGLEBENDTYPE"].size()){
    DEBUG(10, "BONDANGLEBENDTYPE block");
    buffer = m_block["BONDANGLEBENDTYPE"];
    // 1. BONDTYPE 2. number of types
    for (it = buffer.begin() + 2;
    it != buffer.end() - 1; ++it) {
      
      double k, kh, cos0;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> kh >> cos0;
      
      if (_lineStream.fail()){
        os << *it << std::endl;
        io::messages.add("bad line in BONDANGLEBENDTYPE block", "In_Topology",
                io::message::error);
      }
      if (! _lineStream.eof()){
        os << *it << std::endl;
        io::messages.add("eof not reached in BONDANGLEBENDTYPE block",
                "InTopology", io::message::warning);
      }
      
      // and add (force constant based on a potential harmonic in the angle cosine)
      b.push_back(interaction::angle_type_struct(kh, cos0 * 2 * math::Pi / 360.0));
    }
  }
  else if (m_block["HARMBONDANGLETYPE"].size()){
    DEBUG(10, "HARMBONDANGLETYPE block");
     
    buffer = m_block["HARMBONDANGLETYPE"];
    
    if (buffer.size()==0){
      io::messages.add("HARMBONDANGLETYPE block not found!", "In_Topology",
              io::message::error);
      return;
    }
    
    // 1. BONDANGLETYPE 2. number of types
    for (it = buffer.begin() + 2;
    it != buffer.end() - 1; ++it) {
      
      double k, theta;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> theta;
      
      if (_lineStream.fail()){
        os << *it << std::endl;
        io::messages.add("bad line in HARMBONDANGLETYPE block", "In_Topology",
                io::message::error);
      }
      if (! _lineStream.eof()){
        os << *it << std::endl;
        io::messages.add("eof not reached in HARMBONDANGLETYPE block",
                "InTopology", io::message::warning);
      }
      
      // and add...
      b.push_back(interaction::angle_type_struct(k, theta * 2 * math::Pi / 360.0));
    }
  }
  else{
    io::messages.add("either BONDANGLEBENDTYPE or HARMBONDANGLETYPE block must be present",
            "In_Topology", io::message::error);
  }

}

void io::In_Topology
::read_improper_dihedrals(std::vector<interaction::improper_dihedral_type_struct> &i,
			  std::ostream & os)
{
  
  DEBUG(10, "IMPDIHEDRALTYPE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["IMPDIHEDRALTYPE"];

  if (buffer.size()==0){
    io::messages.add("IMPDIHEDRALTYPE block not found!", "In_Topology",
		     io::message::error);
    return;
  }
  
  // 1. IMPDIHEDRALTYPE 2. number of types
  for (it = buffer.begin() + 2; 
       it != buffer.end() - 1; ++it) {

    double k, q0;
    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> k >> q0;
      
    if (_lineStream.fail()){
      os << *it << std::endl;
      io::messages.add("bad line in IMPDIHEDRALTYPE block", "In_Topology",
		       io::message::error);
    }
    if (! _lineStream.eof()){
      os << *it << std::endl;
      io::messages.add("eof not reached in IMPDIHEDRALTYPE block",
		       "InTopology", io::message::warning);
    }

    // and add...
    i.push_back(interaction::improper_dihedral_type_struct(k*180*180/math::Pi/math::Pi,
							   q0 * math::Pi / 180.0));
  }

}

void io::In_Topology
::read_dihedrals(std::vector<interaction::dihedral_type_struct> &d,
		 std::ostream & os)
{
  if (m_block["TORSDIHEDRALTYPE"].size()){
    DEBUG(10, "TORSDIHEDRALTYPE block");
    
    io::messages.add("parsing TORSDIHEDRALTYPE to DIHEDRALTYPE format i.e. PD=cos(PDL)", "In_Topology",
                io::message::warning);
    
    DEBUG(10, "TORSDIHEDRALTYPE block");
    
    std::vector<std::string> buffer;
    std::vector<std::string>::const_iterator it;
    
    buffer = m_block["TORSDIHEDRALTYPE"];
    
    if (buffer.size()==0){
      io::messages.add("TORSDIHEDRALTYPE block not found!", "In_Topology",
              io::message::error);
      return;
    }
    
    // 1. DIHEDRALTYPE 2. number of types
    for (it = buffer.begin() + 2;
    it != buffer.end() - 1; ++it) {
      
      double k, pdl;
      int m;
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> pdl >> m;
      
      if (_lineStream.fail()){
        os << *it << std::endl;
        io::messages.add("bad line in TORSDIHEDRALTYPE block", "In_Topology",
                io::message::error);
      }
      if (! _lineStream.eof()){
        os << *it << std::endl;
        io::messages.add("eof not reached in TORSDIHEDRALTYPE block",
                "InTopology", io::message::warning);
      }
      
      // and add...
      d.push_back(interaction::dihedral_type_struct(k, cos(pdl*math::Pi/180.0), m));
    }
  }
  else if(m_block["DIHEDRALTYPE"].size()){
    DEBUG(10, "DIHEDRALTYPE block");
    
    std::vector<std::string> buffer;
    std::vector<std::string>::const_iterator it;
    
    buffer = m_block["DIHEDRALTYPE"];
    
    if (buffer.size()==0){
      io::messages.add("DIHEDRALTYPE block not found!", "In_Topology",
              io::message::error);
      return;
    }
    
    // 1. DIHEDRALTYPE 2. number of types
    for (it = buffer.begin() + 2;
    it != buffer.end() - 1; ++it) {
      
      double k, pd;
      int m;
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> pd >> m;
      
      if (_lineStream.fail()){
        os << *it << std::endl;
        io::messages.add("bad line in DIHEDRALTYPE block", "In_Topology",
                io::message::error);
      }
      if (! _lineStream.eof()){
        os << *it << std::endl;
        io::messages.add("eof not reached in DIHEDRALTYPE block",
                "InTopology", io::message::warning);
      }
      
      // and add...
      d.push_back(interaction::dihedral_type_struct(k, pd, m));
    }
  }
  else{
    // complain that we need a block
    io::messages.add("either TORSDIHEDRALTYPE or DIHEDRALTYPE block must be present", "In_Topology",
                io::message::error);
  }
}


void io::In_Topology
::read_lj_parameter(std::vector<std::vector
		    <interaction::lj_parameter_struct> > 
		    & lj_parameter,
		    std::ostream & os)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  { // LJPARAMETERS
    
    DEBUG(10, "LJPARAMETERS block");
    
    buffer = m_block["LJPARAMETERS"];
    if (!buffer.size()){
      io::messages.add("No LJPARAMETERS block found in topology!",
		       "In_Topology",
		       io::message::error);
      return;
    }
    
    int num, n;
    
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num;
    
    // calculate the matrix size from: x = n*(n+1)/2
    unsigned int sz = unsigned(sqrt(double((8*num+1)-1))/2);

    lj_parameter.resize(sz);
    std::vector< std::vector<interaction::lj_parameter_struct> >::iterator
      lj_it = lj_parameter.begin(),
      lj_to = lj_parameter.end();
  
    for(; lj_it!=lj_to; ++lj_it)
      lj_it->resize(sz);
    
    ++it;
    
    for (n=0; it != buffer.end() - 1; ++it, ++n) {
      
      interaction::lj_parameter_struct s;
      int i, j;
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> i >> j >> s.c12 >> s.c6 >> s.cs12 >> s.cs6;

      --i;
      --j;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	io::messages.add("bad line in LJPARAMETERS block", "In_Topology",
			 io::message::error);
      
      if (i >= int(sz) || j >= int(sz)){
	DEBUG(7, "wrong iac in LJPARAMETERS: i=" << i << " j=" << j
	      << " sz=" << sz);
	io::messages.add("wrong integer atom code in LJPARAMETERS block",
			 "In_Topology", 
			 io::message::error);
      }

      lj_parameter[i][j] = s;
      lj_parameter[j][i] = s;
      
    }

    if (num != n){
      io::messages.add("Reading the LJPARAMETERS failed (n != num)",
		       "InTopology",
		       io::message::error);
    }
  } // LJPARAMETER

}


void io::In_Topology
::read_cg_parameter(std::vector<std::vector
		    <interaction::lj_parameter_struct> > 
		    & cg_parameter,
		    std::ostream & os)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  { // CGPARAMETERS
    
    DEBUG(10, "CGPARAMETERS block");
    
    buffer = m_block["CGPARAMETERS"];
    if (!buffer.size()){
      io::messages.add("No CGPARAMETERS block found in topology!",
		       "In_Topology",
		       io::message::error);
      return;
    }
    
    int num, n;
    
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num;
    
    // calculate the matrix size from: x = n*(n+1)/2
    unsigned int sz = unsigned(sqrt(double((8*num+1)-1))/2);

    cg_parameter.resize(sz);
    std::vector< std::vector<interaction::lj_parameter_struct> >::iterator
      cg_it = cg_parameter.begin(),
      cg_to = cg_parameter.end();
  
    for(; cg_it!=cg_to; ++cg_it)
      cg_it->resize(sz);
    
    ++it;
    
    for (n=0; it != buffer.end() - 1; ++it, ++n) {
      
      interaction::lj_parameter_struct s;
      int i, j;
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> i >> j >> s.c12 >> s.c6;

      --i;
      --j;
      
      if (_lineStream.fail() || ! _lineStream.eof())
	io::messages.add("bad line in CGPARAMETERS block", "In_Topology",
			 io::message::error);
      
      if (i >= int(sz) || j >= int(sz)){
	DEBUG(7, "wrong iac in CGPARAMETERS: i=" << i << " j=" << j
	      << " sz=" << sz);
	io::messages.add("wrong integer atom code in CGPARAMETERS block",
			 "In_Topology", 
			 io::message::error);
	break;
      }

      // no different 1,4 interactions
      s.cs6 = s.c6;
      s.cs12 = s.c12;
      
      cg_parameter[i][j] = s;
      cg_parameter[j][i] = s;
      
    }

    if (num != n){
      io::messages.add("Reading the CGPARAMETERS failed (n != num)",
		       "InTopology",
		       io::message::error);
    }
   
  } // CGPARAMETER
}

