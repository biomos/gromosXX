/**
 * @file in_topology.tcc
 * implements methods of In_Topology.
 */


#include <util/stdheader.h>

#include <topology/core/core.h>
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
static bool check_type(std::vector<std::string> const & buffer, std::vector<T> term)
{
  if (buffer.size()){
    std::istringstream is(buffer[1]);
    int num;
    if (!(is >> num) || num < 0)
      return false;
      
    for(typename std::vector<T>::const_iterator
	  it = term.begin(),
	  to = term.end();
	it != to;
	++it){
      
      if (int(it->type) >= num){
	return false;
      }
    }
  }
  else return false;
  return true;
}

void 
io::In_Topology::read(topology::Topology& topo,
		      simulation::Parameter &param){

  DEBUG(7, "reading in topology");

  std::cout << "TOPOLOGY\n";

  std::cout << title << "\n";
  
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  { // TOPPHYSCON
    buffer = m_block["TOPPHYSCON"];
    
    std::string s;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1,
				buffer.end()-1, s));
    double four_pi_eps0_i;
    
    _lineStream >> four_pi_eps0_i >> math::h_bar;
    math::four_pi_eps_i = four_pi_eps0_i / param.longrange.epsilon;
    
    if (_lineStream.fail())
      io::messages.add("Bad line in TOPPHYSCON block",
		       "InTopology", io::message::error);
  }

  { // RESNAME
    std::cout << "\tRESNAME\n\t";
    
    DEBUG(10, "RESNAME block");
    buffer = m_block["RESNAME"];
    it = buffer.begin()+1;
    int n, num;
    _lineStream.clear();
    _lineStream.str(*it);
    _lineStream >> num;
    ++it;
    
    for(n=0; n<10; ++n)
      std::cout << std::setw(8) << n+1;
    std::cout << "\n\t";

    for(n=0; it != buffer.end() - 1; ++it, ++n){
      std::string s;
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> s;
      
      if (n && ((n) % 10) == 0) std::cout << std::setw(10) << n << "\n\t";
      std::cout << std::setw(8) << s;      

      topo.residue_names().push_back(s);
    }

    if (n != num){
      io::messages.add("Error in RESNAME block: n!=num.",
		       "InTopology", io::message::error);
      throw std::runtime_error("error in RESNAME block (n != num)");
    }

    std::cout << "\n\tEND\n";
    
    
  } // RESNAME

  { // SOLUTEATOM
    DEBUG(10, "SOLUTEATOM block");
    buffer = m_block["SOLUTEATOM"];
  
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);
    int num, n;
    _lineStream >> num;
    topo.resize(num);

    std::cout << "\tSOLUTEATOM\n\t";
    std::cout << "\tnumber of atoms : " << num;
    
    // put the rest of the block into a single stream
    ++it;
    std::string soluteAtoms;
    concatenate(it, buffer.end()-1, soluteAtoms);

    _lineStream.clear();
    _lineStream.str(soluteAtoms);

    int a_nr, r_nr, t, cg, n_ex, a_ex;
    double m, q;
    std::string s;
    std::set<int> ex;
    std::set<int> ex14;
      
    for(n=0; n < num; ++n){

      _lineStream >> a_nr >> r_nr >> s >> t >> m >> q >> cg >> n_ex;
      
      if (a_nr != n+1){
	io::messages.add("Error in SOLUTEATOM block: atom number not sequential.",
			 "InTopology", io::message::error);
      }

      if (r_nr > int(topo.residue_names().size()) || r_nr < 1){
	io::messages.add("Error in SOLUTEATOM block: residue number out of range.",
			 "InTopology", io::message::error);
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

      if (n_ex < 0){
	io::messages.add("Error in SOLUTEATOM block: number of exclusions < 0.",
			 "InTopology", io::message::error);
      }
      
      // exclusions
      ex.clear();
      for(int i=0; i<n_ex; ++i){
	_lineStream >> a_ex;

	if (a_ex <= a_nr)
	  io::messages.add("Error in SOLUTEATOM block: exclusions only to "
			   "larger atom numbers.",
			   "InTopology", io::message::error);	

	ex.insert(a_ex-1);
      }
      
      // 1,4 - pairs
      _lineStream >> n_ex;

      if (n_ex < 0){
	io::messages.add("Error in SOLUTEATOM block: number of 1,4 exclusions < 0.",
			 "InTopology", io::message::error);
      }

      ex14.clear();
      for(int i=0; i<n_ex; ++i){
	_lineStream >> a_ex;
	if (a_ex <= a_nr)
	  io::messages.add("Error in SOLUTEATOM block: 1,4 - exclusions only to "
			   "larger atom numbers.",
			   "InTopology", io::message::error);	
	
	ex14.insert(a_ex-1);
      }
      
      if (_lineStream.fail())
	throw std::runtime_error("bad line in SOLUTEATOM block");

      topo.add_solute_atom(s, r_nr-1, t-1, m, q, cg, ex, ex14);
    }
    std::cout << "\n\tEND\n";
    
  } // SOLUTEATOM
    
  
  { // BOND
    DEBUG(10, "BOND block");
    buffer = m_block["BOND"];

    std::cout << "\tBOND";
  
    if (buffer.size()){
      
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;

      if (param.shake.ntc == 3){
	std::cout << "\n\t\t"
		  << num
		  << " bonds from BOND block added to CONSTRAINT";
      }
      else
	std::cout << "\n\t\tbonds not containing hydrogens : "
		  << num;
      
      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, t;
	
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> t;
	
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in BOND block",
			   "In_Topology", io::message::error);
	  throw std::runtime_error("bad line in BOND block");
	}
      
	if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	    i < 1 || j < 1){
	  io::messages.add("Atom number out of range in BOND block",
			   "In_Topology", io::message::error);
	}
      
	if (param.shake.ntc == 3){
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
	throw std::runtime_error("error in BOND block (n != num)");
      }
    }
  
  } // BOND

  { // BONDH
    DEBUG(10, "BONDH block");
    
    buffer.clear();
    buffer = m_block["BONDH"];
    if (buffer.size()){
      it = buffer.begin() + 1;
      
      _lineStream.clear();
      _lineStream.str(*it);
      
      int num, n;
      _lineStream >> num;
      ++it;
      
      if (param.shake.ntc == 2 || param.shake.ntc == 3){
	std::cout << "\n\t\t"
		  << num
		  << " bonds from BONDH block added to CONSTRAINT";
      }
      else
	std::cout << "\n\t\tbonds containing hydrogens : "
		  << num;
      
      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, t;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> t;
	
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in BONDH block",
			   "In_Topology", io::message::error);
	  throw std::runtime_error("bad line in BONDH block");
	}
	
	if (i > int(topo.num_solute_atoms()) || j > int(topo.num_solute_atoms()) ||
	    i < 1 || j < 1){
	  io::messages.add("Atom number out of range in BONDH block",
			   "In_Topology", io::message::error);
	}
	
	if (param.shake.ntc == 2 || param.shake.ntc == 3){
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
	throw std::runtime_error("error in BONDH block (n != num)");
      }
    }

    std::cout << "\n\tEND\n";
    
  } // BONDH

  // check the bonds
  if (!check_type(m_block["BONDTYPE"], topo.solute().bonds())){
    io::messages.add("Illegal bond type in BOND(H) block",
		     "In_Topology", io::message::error);
  }

  { // CONSTRAINT
    DEBUG(10, "CONSTRAINT block");
    buffer = m_block["CONSTRAINT"];
  
    if (buffer.size() && param.shake.ntc != 1){
      
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;

      std::cout << "\tCONSTRAINT\n\t\t"
		<< num
		<< " bonds in CONSTRAINT block."
		<< "\n\t\ttotal of constraint bonds : " 
		<< num + topo.solute().distance_constraints().size()
		<< "\n\tEND\n";
      
      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, t;
	
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> t;
	
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in CONSTRAINT block",
			   "In_Topology", io::message::error);
	  throw std::runtime_error("bad line in CONSTRAINT block");
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
	throw std::runtime_error("error in CONSTRAINT block (n != num)");
      }
    }
    
  } // CONSTRAINT

  // check the bonds in constraints
  if (!check_type(m_block["BONDTYPE"], topo.solute().distance_constraints())){
    io::messages.add("Illegal bond type in CONSTRAINT (or BOND(H)) block",
		     "In_Topology", io::message::error);
  }

  { // BONDANGLEH

    std::cout << "\tBONDANGLE";

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

      std::cout << "\n\t\tbondangles not containing hydrogens : " << num;

      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, k, t;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> t;
      
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in BONDANGLEH block",
			   "In_Topology", io::message::error);
	  throw std::runtime_error("bad line in BONDANGLEH block");
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
	// if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in BONDANGLEH block (n != num)");
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

      std::cout << "\n\t\tbondangles containing hydrogens : " << num;
    
      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, k, t;
      
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> t;
      
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in BONDANGLE block",
			   "In_Topology", io::message::error);
	  throw std::runtime_error("bad line in BONDANGLE block");
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
	// if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in BONDANGLE block (n != num)");
      }
    }

    std::cout << "\n\tEND\n";

  } // BONDANGLE


  // check the angles
  if (!check_type(m_block["BONDANGLETYPE"], topo.solute().angles())){
    io::messages.add("Illegal bond angle type in BONDANGLE(H) block",
		     "In_Topology", io::message::error);
  }

  { // IMPDIHEDRAL
    DEBUG(10, "IMPDIHEDRAL block");
    buffer = m_block["IMPDIHEDRAL"];
  
    std::cout << "\tIMPDIHEDRAL";

    if(buffer.size()){
      
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;

      std::cout << "\n\t\timproper dihedrals not containing hydrogens : "
		<< num;
    
      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, k, l, t;
      
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> l >> t;
      
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in IMPDIHEDRAL block",
			   "In_Topology", io::message::error);
	  throw std::runtime_error("bad line in IMPDIHEDRAL block");
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
	// if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in IMPDIHEDRAL block (n != num)");
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

      std::cout << "\n\t\timproper dihedrals containing hydrogens : "
		<< num;

      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, k, l, t;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> l >> t;
      
      
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in IMPDIHEDRALH block",
			   "In_Topology", io::message::error);
	  throw std::runtime_error("bad line in IMPDIHEDRALH block");
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
	// if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in IMPDIHEDRALH block (n != num)");
      }
    }

    std::cout << "\n\tEND\n";
    
  } // IMPDIHEDRALH

  // check the imporopers
  if (!check_type(m_block["IMPDIHEDRALTYPE"], topo.solute().improper_dihedrals())){
    io::messages.add("Illegal improper dihedral type in IMPDIHEDRAL(H) block",
		     "In_Topology", io::message::error);
  }  

  { // DIHEDRAL
    DEBUG(10, "DIHEDRAL block");    
    buffer = m_block["DIHEDRAL"];

    std::cout << "\tDIHEDRAL";
    
    if(buffer.size()){
      
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
    
      std::cout << "\n\t\tdihedrals not containing hydrogens : "
		<< num;

      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, k, l, t;
      
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> l >> t;
      
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in DIHEDRAL block",
			   "In_Topology", io::message::error);
	  throw std::runtime_error("bad line in DIHEDRAL block");
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
	// if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in DIHEDRAL block (n != num)");
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

      std::cout << "\n\t\tdihedrals containing hydrogens : "
		<< num;

      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, k, l, t;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> l >> t;
      
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in DIHEDRALH block",
			   "In_Topology", io::message::error);
	  throw std::runtime_error("bad line in DIHEDRALH block");
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
	// if (_lineStream.fail()|| ! _lineStream.eof())
	throw std::runtime_error("error in DIHEDRALH block (n != num)");
      }
    }
    
    std::cout << "\n\tEND\n";
    
  } // DIHEDRALH

  // check the dihedrals
  if (!check_type(m_block["DIHEDRALTYPE"], topo.solute().dihedrals())){
    io::messages.add("Illegal dihedral type in DIHEDRAL(H) block",
		     "In_Topology", io::message::error);
  }
      
  { // SOLVENTATOM and SOLVENTCONSTR
    // give it a number (SOLVENTATOM1, SOLVENTATOM2) for multiple
    // solvents...
    DEBUG(10, "SOLVENTATOM block");
    buffer = m_block["SOLVENTATOM"];

    std::cout << "\tSOLVENT";
    
    if (buffer.size()){
      
      int res_nr = topo.residue_names().size();
    
      topo.residue_names().push_back("SOLV");

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
    
      std::cout << "\n\t\tatoms : " << num;

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
	  throw std::runtime_error("bad line in SOLVENTATOM block");
	}
      

	s.add_atom(name, res_nr, iac-1, mass, charge);
      }
    
      if (n!=num){
	io::messages.add("Error in SOLVENTATOM block (num != n)",
			 "In_Topology", io::message::error);
	throw std::runtime_error("error in SOLVENTATOM block (n != num)");
      }
    
      // lookup the number of bond types
      // add additional ones for the solvent constraints
      int num_bond_types;
    
      buffer = m_block["BONDTYPE"];
      if (!buffer.size())
	num_bond_types = 0;
      else{
	it = buffer.begin() + 1;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> num_bond_types;
      
	if (num_bond_types < 0){
	  io::messages.add("Illegal value for number of bond types in BONDTYPE block",
			   "In_Topology", io::message::error);	
	  num_bond_types = 0;
	}
      }
    
      buffer = m_block["SOLVENTCONSTR"];
      if (!buffer.size()){
	io::messages.add("SOLVENTCONSTR block missing.",
			 "In_Topology", io::message::error);	

	throw std::runtime_error("SOLVENTCONSTR block missing.");
      }
    
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> num;
      ++it;
    
      std::cout << "\n\t\tconstraints : " << num;
      
      int j;
      double b0;
    
      for(n=0; it != buffer.end()-1; ++it, ++n){
	_lineStream.clear();
	_lineStream.str(*it);
      
	_lineStream >> i >> j >> b0;
      
	if (_lineStream.fail() || ! _lineStream.eof()){
	  io::messages.add("Bad line in SOLVENTCONSTR block",
			   "In_Topology", io::message::error);	
	  throw std::runtime_error("bad line in SOLVENTCONSTR block");
	}
      
	// the solvent (distance constraints) bond types
	s.add_distance_constraint
	  (topology::two_body_term_struct(i-1, j-1, num_bond_types + n));
      }

      if (n!=num){
	io::messages.add("Error in SOLVENTCONSTR block (num != n)",
			 "In_Topology", io::message::error);
      
	throw std::runtime_error("error in SOLVENTCONSTR block (n != num)");

      }
      topo.add_solvent(s);
    }

  }

  // add the solvent to the topology
  std::cout << "\n\t\tadding " << param.system.nsm 
	    << " solvents.";
  
  if (param.system.nsm) topo.solvate(0, param.system.nsm);  

  std::cout << "\n\tEND\n";

  // set lambda
  topo.lambda(param.perturbation.lambda);
  topo.lambda_exp(param.perturbation.lambda_exponent);

  //==================================================
  // CHECKING
  //==================================================
    
  // submolecules check
  if (param.submolecules.submolecules[param.submolecules.submolecules.size()-1]
      != topo.num_solute_atoms()){
    
    io::messages.add("Error in SUBMOLECULE block: "
		     "last submolecule has to end with last solute atom",
		     "In_Topology", io::message::error);
  }

  // add the submolecules
  topo.molecules() = param.submolecules.submolecules;

  // energy group check
  if (param.force.energy_group[param.force.energy_group.size()-1]
      != topo.num_atoms()-1){
    io::messages.add("Error in FORCE block: "
		     "last energy group has to end with last atom",
		     "In_Topology", io::message::error);
  }
  // and add them
  size_t atom = 0;
  for(size_t i=0; i<param.force.energy_group.size(); ++i){
    topo.energy_groups().push_back(param.force.energy_group[i]);
    for( ; atom <= param.force.energy_group[i]; ++atom){
      topo.atom_energy_group().push_back(i);
      // DEBUG(11, "atom " << atom << ": " << i);
    }
  }

  if(!param.multibath.found_multibath && param.multibath.found_tcouple){
    std::cout << "\tparsing a (deprecated) TCOUPLE block into the new "
	      << "MULTIBATH format.\n";
    
    util::parse_TCOUPLE(param, topo);
  }

  std::cout << "END\n";
  
}

void io::In_Topology
::read_harmonic_bonds(std::vector<interaction::bond_type_struct> &b)
{
  
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

      if (_lineStream.fail() || ! _lineStream.eof())
	throw std::runtime_error("bad line in HARMBONDTYPE block");
      
      // and add...
      b.push_back(interaction::bond_type_struct(k, r));
    }
  
    if (num != n)
      throw std::runtime_error("not enough bond types in HARMBONDTYPE block");
  
  }
  else{
    buffer = m_block["BONDTYPE"];

    io::messages.add("converting bond force constants from quartic "
		     "to harmonic form", "InTopology::bondtype",
		     io::message::notice);

    if (buffer.size()==0)
      throw std::runtime_error("BONDTYPE block not found!");

    // 1. BONDTYPE 2. number of types
    for (it = buffer.begin() + 2; 
	 it != buffer.end() - 1; ++it) {

      double k, r;
      _lineStream.clear();
      _lineStream.str(*it);
      
      _lineStream >> k >> r;
      
      if (_lineStream.fail()){
	std::cout << *it << std::endl;
	throw std::runtime_error("bad line in BONDTYPE block");
      }
      if (! _lineStream.eof()){
	std::cout << *it << std::endl;
	io::messages.add("eof not reached in BONDTYPE block",
			 "InTopology", io::message::warning);
      }

      // we are reading into harmonic bond term, so convert k
      k *= 2 * r * r;
      
      // and add...
      b.push_back(interaction::bond_type_struct(k, r));
    }
  }

  // also add the solent constraints to the bond types...
  // (if there is one)
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
	throw std::runtime_error("bad line in SOLVENTCONSTR block");
      }
      
      // the solvent (distance constraints) bond types
      b.push_back(interaction::bond_type_struct(0, b0));
      // (K is set to 0.0)
    }

    if (n!=num){
      io::messages.add("Error in SOLVENTCONSTR block (num != n)",
		       "In_Topology", io::message::error);
      
      throw std::runtime_error("error in SOLVENTCONSTR block (n != num)");
      
    }
  }
  
}

void io::In_Topology
::read_g96_bonds(std::vector<interaction::bond_type_struct> &b)
{
  
  DEBUG(10, "BONDTYPE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["BONDTYPE"];

  if (buffer.size()==0)
    throw std::runtime_error("BONDTYPE block not found!");

  // 1. BONDTYPE 2. number of types
  for (it = buffer.begin() + 2; 
       it != buffer.end() - 1; ++it) {

    double k, r;
    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> k >> r;
      
    if (_lineStream.fail()){
      std::cout << *it << std::endl;
      throw std::runtime_error("bad line in BONDTYPE block");
    }
    if (! _lineStream.eof()){
      std::cout << *it << std::endl;
      io::messages.add("eof not reached in BONDTYPE block",
		       "InTopology", io::message::warning);
    }

    // and add...
    b.push_back(interaction::bond_type_struct(k, r));
  }

  // also add the solent constraints to the bond types...
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
	throw std::runtime_error("bad line in SOLVENTCONSTR block");
      }
      
      // the solvent (distance constraints) bond types
      // (K is set to 0.0)
      b.push_back(interaction::bond_type_struct(0, b0));
    }

    if (n!=num){
      io::messages.add("Error in SOLVENTCONSTR block (num != n)",
		       "In_Topology", io::message::error);
      
      throw std::runtime_error("error in SOLVENTCONSTR block (n != num)");
      
    }
  }
  
}

void io::In_Topology
::read_angles(std::vector<interaction::angle_type_struct> &b)
{
  
  DEBUG(10, "BONDANGLETYPE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["BONDANGLETYPE"];

  if (buffer.size()==0)
    throw std::runtime_error("BONDANGLETYPE block not found!");

  // 1. BONDTYPE 2. number of types
  for (it = buffer.begin() + 2; 
       it != buffer.end() - 1; ++it) {

    double k, cos0;
    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> k >> cos0;
      
    if (_lineStream.fail()){
      std::cout << *it << std::endl;
      throw std::runtime_error("bad line in BONDANGLETYPE block");
    }
    if (! _lineStream.eof()){
      std::cout << *it << std::endl;
      io::messages.add("eof not reached in BONDANGLETYPE block",
		       "InTopology", io::message::warning);
    }

    // and add...
    b.push_back(interaction::angle_type_struct(k, cos(cos0 * 2 * math::Pi / 360.0)));
  }

}

void io::In_Topology
::read_improper_dihedrals(std::vector<interaction::improper_dihedral_type_struct> &i)
{
  
  DEBUG(10, "IMPDIHEDRALTYPE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["IMPDIHEDRALTYPE"];

  if (buffer.size()==0)
    throw std::runtime_error("IMPDIHEDRALTYPE block not found!");

  // 1. IMPDIHEDRALTYPE 2. number of types
  for (it = buffer.begin() + 2; 
       it != buffer.end() - 1; ++it) {

    double k, q0;
    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> k >> q0;
      
    if (_lineStream.fail()){
      std::cout << *it << std::endl;
      throw std::runtime_error("bad line in IMPDIHEDRALTYPE block");
    }
    if (! _lineStream.eof()){
      std::cout << *it << std::endl;
      io::messages.add("eof not reached in IMPDIHEDRALTYPE block",
		       "InTopology", io::message::warning);
    }

    // and add...
    i.push_back(interaction::improper_dihedral_type_struct(k*180*180/math::Pi/math::Pi,
							   q0 * math::Pi / 180.0));
  }

}

void io::In_Topology
::read_dihedrals(std::vector<interaction::dihedral_type_struct> &d)
{
  
  DEBUG(10, "DIHEDRALTYPE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["DIHEDRALTYPE"];

  if (buffer.size()==0)
    throw std::runtime_error("DIHEDRALTYPE block not found!");

  // 1. DIHEDRALTYPE 2. number of types
  for (it = buffer.begin() + 2; 
       it != buffer.end() - 1; ++it) {

    double k, pd;
    int m;

    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> k >> pd >> m;
      
    if (_lineStream.fail()){
      std::cout << *it << std::endl;
      throw std::runtime_error("bad line in DIHEDRALTYPE block");
    }
    if (! _lineStream.eof()){
      std::cout << *it << std::endl;
      io::messages.add("eof not reached in DIHEDRALTYPE block",
		       "InTopology", io::message::warning);
    }

    // and add...
    d.push_back(interaction::dihedral_type_struct(k, pd, m));
  }

}


void io::In_Topology
::read_lj_parameter(std::vector<std::vector
		    <interaction::lj_parameter_struct> > 
		    & lj_parameter)
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
    size_t sz = size_t(sqrt(double((8*num+1)-1))/2);

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
	throw std::runtime_error("bad line in LJPARAMETERS block");
      
      if (i >= int(sz) || j >= int(sz)){
	DEBUG(7, "wrong iac in LJPARAMETERS: i=" << i << " j=" << j
	      << " sz=" << sz);
	throw std::string("wrong integer atom code in LJPARAMETERS block");
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

