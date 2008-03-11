/**
 * @file in_perturbation.cc
 * implements methods of In_Perturbation.
 */

#include <stdheader.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <io/instream.h>

#include <io/blockinput.h>

#include "in_perturbation.h"


#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

/**
 * Constructor.
 */

io::In_Perturbation::In_Perturbation(std::istream &is) 
  : GInStream(is) 
{
  // read the whole file at beginning
  readStream();
};

template<typename T>
bool check_type(std::vector<std::string> const & buffer, std::vector<T> term)
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
      
      if (int(it->A_type) >= num && int(it->B_type) >= num)
	return false;
    }
  }
  else return false;
  return true;
}

void
io::In_Perturbation::read(topology::Topology &topo,
			  simulation::Parameter &param)
{

  if (!param.perturbation.perturbation){
    io::messages.add("Ignoring perturbation topology because perturbation is not enabled.",
		     "In_Perturbation",
		     io::message::warning);
    return;
  }

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  DEBUG(7, "Reading PERTURBATION");
  if (!quiet){
    std::cout << "PERTURBATION TOPOLOGY\n";
    std::cout << title << "\n";
  }

  // prepare arrays
  topo.is_perturbed().resize(topo.num_solute_atoms(), false);

  { // PERTBONDSTRETCH(H)
      
    if (buffer.size()){
        
    block_read.insert("PERTBOND03");
    io::messages.add("The PERTBOND03 block was renamed to PERTBONDSTRETCH.",
             "In_Perturbation", io::message::error); 
    }
    std::vector<std::string> pertbondstretch;  
    pertbondstretch.push_back("PERTBONDSTRETCHH");
    pertbondstretch.push_back("PERTBONDSTRETCH");
    for (unsigned int hh=0; hh < pertbondstretch.size(); hh++) {
    buffer = m_block[pertbondstretch.at(hh)];
    if (buffer.size()){
            
      block_read.insert(pertbondstretch.at(hh));
      
      if (!quiet)
	std::cout << "\t" << pertbondstretch.at(hh) << "\n";

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
      
      if (param.constraint.ntc == 2){
	io::messages.add("No perturbed distance constraints for "
			 "NTC = 2 from perturbed bonds",
			 "in_perturbation",
			 io::message::warning);
      }
      else if (param.constraint.ntc == 3){
	if (!quiet)
	  std::cout << "\n\t\t"
		    << num
		    << " perturbed bonds from " << pertbondstretch.at(hh) << " block added to "
		    << "perturbed distance constraints.";
      }

      if (!quiet)
	std::cout << "\n\t"
		  << std::setw(10) << "atom i"
		  << std::setw(10) << "atom j"
		  << std::setw(10) << "type A"
		  << std::setw(10) << "type B"
		  << "\n";
      
      for(n=0; it != buffer.end() -1; ++it, ++n){
	int i, j, t_A, t_B;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> t_A >> t_B;
	
	if (_lineStream.fail()){
	  io::messages.add("Bad line in " + pertbondstretch.at(hh) + " block.",
			   "In_Perturbation", io::message::error);
	}
	
	topology::two_body_term_struct b(i-1, j-1, t_A-1);

	if (param.constraint.ntc != 3){
	  std::vector<topology::two_body_term_struct>::iterator b_it
	    = std::find(topo.solute().bonds().begin(), 
			topo.solute().bonds().end(), 
			b);
	
	  if (b_it == topo.solute().bonds().end()){
		// or B-type
		b.type = t_B-1;
		b_it = std::find(topo.solute().bonds().begin(),
						 topo.solute().bonds().end(),
						 b);
		if (b_it == topo.solute().bonds().end()){
		  // should this be a warning?
		  std::cout << "perturbing " << i << " - " << j << " type " << t_A << std::endl;
		  std::cout << "\tnot in topology!" << std::endl;

	      io::messages.add("Perturbation of a non-existing bond "
		          	       "in " + pertbondstretch.at(hh) + " block.",
			     		   "In_Perturbation", io::message::error);
	    }
	  }

	  topo.solute().bonds().erase(b_it);
	  topology::perturbed_two_body_term_struct 
	    pb(i-1, j-1, t_A-1, t_B-1);

	  if (!quiet)
	    std::cout << "\t" 
		      << std::setw(10) << pb.i+1 
		      << std::setw(10) << pb.j+1
		      << std::setw(10) << pb.A_type+1 
		      << std::setw(10) << pb.B_type+1 
		      << "\n";
	  
	  topo.perturbed_solute().bonds().push_back(pb);
	}
	else{
	  std::vector<topology::two_body_term_struct>::iterator b_it
	    = std::find(topo.solute().distance_constraints().begin(), 
			topo.solute().distance_constraints().end(), 
			b);
	  
	  if (b_it == topo.solute().distance_constraints().end()){
	    io::messages.add("Perturbation of a non-existing distance "
			     "constraint in " + pertbondstretch.at(hh) + " block.",
			     "In_Perturbation", io::message::error);
	  }
	
	  topo.solute().distance_constraints().erase(b_it);
	  topology::perturbed_two_body_term_struct 
	    pb(i-1, j-1, t_A-1, t_B-1);

	  if (!quiet)
	    std::cout << "\t" 
		      << std::setw(10) << pb.i+1 
		      << std::setw(10) << pb.j+1
		      << std::setw(10) << pb.A_type+1 
		      << std::setw(10) << pb.B_type+1 
		      << "\n";
	  
	  topo.perturbed_solute().distance_constraints().push_back(pb);

	}

      }
      
      if (n != num){
	  io::messages.add("Wrong number of bonds in " + pertbondstretch.at(hh) + " block.",
			   "In_Perturbation", io::message::error);	
      }
      else if (_lineStream.fail()){
	io::messages.add("Bad line in " + pertbondstretch.at(hh) + " block.",
			 "In_Perturbation", io::message::error);
      }

      if (!quiet)
	std::cout << "\n\t\tbonds :                          " 
		  << unsigned(topo.solute().bonds().size())
		  << "\n\t\tperturbed bonds :                "
		  << unsigned(topo.perturbed_solute().bonds().size())
		  << "\n\t\tdistance constraints :           "
		  << unsigned(topo.solute().distance_constraints().size())
		  << "\n\t\tperturbed distance constraints : "
		  << unsigned(topo.perturbed_solute().distance_constraints().size())
		  << "\n\n"
		  << "\tEND\n";

    } // if block present
   }  // loop over H/non H blocks
  } // PERTBONDSTRETCH(H)

  { // PERTCONSTRAINT03
    DEBUG(10, "PERTCONSTRAINT03 block");
    buffer = m_block["PERTCONSTRAINT03"];
  
    if (buffer.size() && param.constraint.ntc != 1){
      block_read.insert("PERTCONSTRAINT03");

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;

      if (!quiet)
	std::cout << "\tPERTCONSTRAINT03\n\t\t"
		  << num
		  << " bonds in PERTCONSTRAINT03 block."
		  << "\n\t\ttotal of perturbed constraint bonds : " 
		  << unsigned(num + topo.perturbed_solute().distance_constraints().size())
		  << "\n"  
		  << "\t"
		  << std::setw(10) << "atom i"
		  << std::setw(10) << "atom j"
		  << std::setw(10) << "type A"
		  << std::setw(10) << "type B"
		  << "\n";
    
      for(n=0; it != buffer.end() - 1; ++it, ++n){
	int i, j, t_A, t_B;
	
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> t_A >> t_B;
	
	if (_lineStream.fail()){
	  io::messages.add("Bad line in PERTCONSTRAINT03 block",
			   "In_Topology", io::message::error);
	}
      
	if (i > int(topo.num_solute_atoms()) || 
	    j > int(topo.num_solute_atoms()) ||
	    i < 1 || j < 1){
	  io::messages.add("Atom number out of range in PERTCONSTRAINT03 "
			   " block", "In_Topology", io::message::error);
	}
      
	topology::two_body_term_struct b(i-1, j-1, t_A-1);
	
	std::vector<topology::two_body_term_struct>::iterator b_it
	  = std::find(topo.solute().distance_constraints().begin(), 
		      topo.solute().distance_constraints().end(), 
		      b);
	  
	if (b_it == topo.solute().distance_constraints().end()){
	  io::messages.add("Perturbation of a non-existing distance "
			   "constraint in PERTCONSTRAINT03 block.",
			   "In_Perturbation", io::message::error);
	  
	}
	
	topo.solute().distance_constraints().erase(b_it);
	topology::perturbed_two_body_term_struct 
	  pb(i-1, j-1, t_A-1, t_B-1);

	topo.perturbed_solute().distance_constraints().push_back(pb);

	if (!quiet)
	  std::cout << "\t" 
		    << std::setw(10) << pb.i+1 
		    << std::setw(10) << pb.j+1
		    << std::setw(10) << pb.A_type+1 
		    << std::setw(10) << pb.B_type+1 
		    << "\n";
	
      }
    
      if(n != num){
	io::messages.add("Wrong number of bonds in PERTCONSTRAINT03 block",
			 "In_Perturbation", io::message::error);
      }
    }
    
  } // PERTCONSTRAINT03

  { // PERTBONDANGLE(H)
      
    if (buffer.size()){
        
    block_read.insert("PERTBANGLE03");
    io::messages.add("The PERTBANGLE03 block was renamed to PERTBONDANGLE.",
             "In_Perturbation", io::message::error); 
    }
      
    std::vector<std::string> pertbondangle;  
    pertbondangle.push_back("PERTBONDANGLEH");
    pertbondangle.push_back("PERTBONDANGLE");
    for (unsigned int hh=0; hh < pertbondangle.size(); hh++) {
    buffer = m_block[pertbondangle.at(hh)];
    if (buffer.size()){
      block_read.insert(pertbondangle.at(hh));
      
      if (!quiet)
	std::cout << "\t" << pertbondangle.at(hh) << "\n";
   
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
      
      if (!quiet)
	std::cout << "\t"
		  << std::setw(10) << "atom i"
		  << std::setw(10) << "atom j"
		  << std::setw(10) << "atom k"
		  << std::setw(10) << "type A"
		  << std::setw(10) << "type B"
		  << "\n";

      for(n=0; it != buffer.end() -1; ++it, ++n){
	int i, j, k, t_A, t_B;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> t_A >> t_B;
	
	if (_lineStream.fail()){
	  io::messages.add("Bad line in " + pertbondangle.at(hh) + " block.",
			   "In_Perturbation", io::message::error);
	}
	
	topology::three_body_term_struct a(i-1, j-1, k-1, t_A-1);
	std::vector<topology::three_body_term_struct>::iterator a_it
	  = std::find(topo.solute().angles().begin(), 
		      topo.solute().angles().end(), 
		      a);
	
	if (a_it == topo.solute().angles().end()){
	  
	  // try for B
	  topology::three_body_term_struct a2(i-1, j-1, k-1, t_B-1);
	  a_it = std::find(topo.solute().angles().begin(), 
			   topo.solute().angles().end(), 
			   a2);
	  if(a_it == topo.solute().angles().end()){
	    
	    io::messages.add("Perturbation of a non-existing angle in "
			     + pertbondangle.at(hh) + " block.",
			     "In_Perturbation", io::message::error);	
	  }
	}
	
	
	if (a_it != topo.solute().angles().end())
	  topo.solute().angles().erase(a_it);
	//	if (a_it2 != topo.solute().angles().end())
	// topo.solute().angles().erase(a_it2);

	topology::perturbed_three_body_term_struct pa(i-1, j-1, k-1, t_A-1, t_B-1);
	topo.perturbed_solute().angles().push_back(pa);

	if (!quiet)
	  std::cout << "\t"
		    << std::setw(10) << pa.i+1 
		    << std::setw(10) << pa.j+1
		    << std::setw(10) << pa.k+1
		    << std::setw(10) << pa.A_type+1 
		    << std::setw(10) << pa.B_type+1 
		    << "\n";

      }
      if (n != num){
	  io::messages.add("Wrong number of angles in " + pertbondangle.at(hh) + " block.",
			   "In_Perturbation", io::message::error);	
      }
      else if (_lineStream.fail()){
	io::messages.add("Bad line in " + pertbondangle.at(hh) + " block.",
			 "In_Perturbation", io::message::error);
      }
      
      if (!quiet)
	std::cout << "\tEND\n";
      
    } // if block present
   } // loop over H/non H blocks
  } // PERTBONDANGLE(H)
  
  { // PERTIMPROPERDIH(H)
      
    if (buffer.size()){
        
    block_read.insert("PERTIMPDIHEDRAL03");
    io::messages.add("The PERTIMPDIHEDRAL03 block was renamed to PERTIMPROPERDIH.",
             "In_Perturbation", io::message::error); 
    }
      
    std::vector<std::string> pertimproperdih;  
    pertimproperdih.push_back("PERTIMPROPERDIHH");
    pertimproperdih.push_back("PERTIMPROPERDIH");
    for (unsigned int hh=0; hh < pertimproperdih.size(); hh++) {
    buffer = m_block[pertimproperdih.at(hh)];
    if (buffer.size()){
      block_read.insert(pertimproperdih.at(hh));
      
      if (!quiet)
	std::cout << "\t" << pertimproperdih.at(hh) << "\n";
   
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
      
      if (!quiet)
	std::cout << "\t"
		  << std::setw(10) << "atom i"
		  << std::setw(10) << "atom j"
		  << std::setw(10) << "atom k"
		  << std::setw(10) << "atom l"
		  << std::setw(10) << "type A"
		  << std::setw(10) << "type B"
		  << "\n";

      for(n=0; it != buffer.end() -1; ++it, ++n){
	int i, j, k, l, t_A, t_B;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> l >> t_A >> t_B;
	
	if (_lineStream.fail()){
	  io::messages.add("Bad line in " + pertimproperdih.at(hh) + " block.",
			   "In_Perturbation", io::message::error);
	}
	
	topology::four_body_term_struct id(i-1, j-1, k-1,l-1, t_A-1);
	std::vector<topology::four_body_term_struct>::iterator id_it
	  = std::find(topo.solute().improper_dihedrals().begin(), 
		      topo.solute().improper_dihedrals().end(), 
		      id);
	
	if (id_it == topo.solute().improper_dihedrals().end()){
	  io::messages.add("Perturbation of a non-existing improper dihedral in "
			   + pertimproperdih.at(hh) + " block.",
			   "In_Perturbation", io::message::error);
	}
	
	topo.solute().improper_dihedrals().erase(id_it);
	topology::perturbed_four_body_term_struct pid(i-1, j-1, k-1, l-1, 
						      t_A-1, t_B-1);
	
	  if (!quiet)
	    std::cout << "\t"
		      << std::setw(10) << pid.i+1 
		      << std::setw(10) << pid.j+1
		      << std::setw(10) << pid.k+1
		      << std::setw(10) << pid.l+1
		      << std::setw(10) << pid.A_type+1 
		      << std::setw(10) << pid.B_type+1 
		      << "\n";
	
	topo.perturbed_solute().improper_dihedrals().push_back(pid);
      }
      
      if (n != num){
	  io::messages.add("Wrong number of bonds in " + pertimproperdih.at(hh) + " block.",
			   "In_Perturbation", io::message::error);	
      }
      else if (_lineStream.fail()){
	io::messages.add("Bad line in " + pertimproperdih.at(hh) + " block.",
			 "In_Perturbation", io::message::error);
      }
      if (!quiet)
	std::cout << "\tEND\n";

    } // if block present
   } // loop over H/non H blocks
  } // PERTIMPROPERDIH(H)

  { // PERTPROPERDIH(H)
      
    std::vector<std::string> pertproperdih;  
    pertproperdih.push_back("PERTPROPERDIHH");
    pertproperdih.push_back("PERTPROPERDIH");
    for (unsigned int hh=0; hh < pertproperdih.size(); hh++) {
    buffer = m_block[pertproperdih.at(hh)];
    if (buffer.size()){
      block_read.insert(pertproperdih.at(hh));
      
      if (!quiet)
	std::cout << "\t" << pertproperdih.at(hh) << "\n"; 
    
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;

      if (!quiet)
	std::cout << "\t"
		  << std::setw(10) << "atom i"
		  << std::setw(10) << "atom j"
		  << std::setw(10) << "atom k"
		  << std::setw(10) << "atom l"
		  << std::setw(10) << "type A"
		  << std::setw(10) << "type B"
		  << "\n";
      
      for(n=0; it != buffer.end() -1; ++it, ++n){
	int i, j, k, l, t_A, t_B;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> l >> t_A >> t_B;
	
	if (_lineStream.fail()){
	  io::messages.add("Bad line in " + pertproperdih.at(hh) + " block.",
			   "In_Perturbation", io::message::error);	  
	}
	
	topology::four_body_term_struct id(i-1, j-1, k-1,l-1, t_A-1);
	std::vector<topology::four_body_term_struct>::iterator id_it
	  = std::find(topo.solute().dihedrals().begin(), 
		      topo.solute().dihedrals().end(), 
		      id);
	
	if (id_it == topo.solute().dihedrals().end()){
	  // try B type
	  topology::four_body_term_struct idB(i-1, j-1, k-1,l-1, t_B-1);
	  id_it =  std::find(topo.solute().dihedrals().begin(), 
			     topo.solute().dihedrals().end(), 
			     idB);

	  if (id_it == topo.solute().dihedrals().end()){
	    io::messages.add("Perturbation of a non-existing dihedral in "
			     + pertproperdih.at(hh) + " block.",
			     "In_Perturbation", io::message::error);
	    return;
	  }
	}
	
	topo.solute().dihedrals().erase(id_it);
	topology::perturbed_four_body_term_struct pid(i-1, j-1, k-1, l-1, 
						      t_A-1, t_B-1);

	if (!quiet)
	  std::cout << "\t"
		    << std::setw(10) << pid.i+1 
		    << std::setw(10) << pid.j+1
		    << std::setw(10) << pid.k+1
		    << std::setw(10) << pid.l+1
		    << std::setw(10) << pid.A_type+1 
		    << std::setw(10) << pid.B_type+1 
		    << "\n";
	
	topo.perturbed_solute().dihedrals().push_back(pid);
      }
      
      if (n != num){
	  io::messages.add("Wrong number of bonds in " + pertproperdih.at(hh) + " block.",
			   "In_Perturbation", io::message::error);	
      }
      else if (_lineStream.fail()){
	io::messages.add("Bad line in " + pertproperdih.at(hh) + " block.",
			 "In_Perturbation", io::message::error);
      }

      if (!quiet)
	std::cout << "\tEND\n";

    } // if block present
   } // loop over H/non H blocks
  } // PERTPROPERDIH(H)

  { // PERTATOMPAIR03
    // has to be read in before(!!) PERTATOMPARAM
    // because the exclusions and 1,4 exclusions have to be adapted...

    buffer = m_block["PERTATOMPAIR03"];
    if (buffer.size()){
      block_read.insert("PERTATOMPAIR03");

      if (!quiet)
	std::cout << "\tPERTATOMPAIR\n";

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;

      int i, j, A, B;

      if (!quiet)
	std::cout << "\t"
		  << std::setw(10) << "atom i"
		  << std::setw(10) << "atom j"
		  << std::setw(10) << "type A"
		  << std::setw(10) << "type B"
		  << "\n";

      for(n = 0; it != buffer.end() - 1; ++it, ++n){
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> A >> B;
	
	if (_lineStream.fail()){
	  io::messages.add("Bad line in PERTATOMPAIR03 block.",
			   "In_Perturbation", io::message::error);
	}
	
	if(j<i) { int t=j; j=i; i=t; }

	topology::perturbed_two_body_term_struct ap(i-1,j-1,A,B);
	
	if (!quiet)
	  std::cout << "\t"
		    << std::setw(10) << ap.i+1
		    << std::setw(10) << ap.j+1
		    << std::setw(10) << ap.A_type
		    << std::setw(10) << ap.B_type
		    << std::endl;
	
	topo.perturbed_solute().atompairs().push_back(ap);

	// make sure it's excluded
	if (topo.all_exclusion(ap.i).count(ap.j) != 1){
	  topo.all_exclusion(ap.i).insert(ap.j);
	  DEBUG(7, "excluding perturbed pair " << ap.i << " and " << ap.j);
	  
	}
	else{
	  // it was already excluded, let's remove it from the
	  // exclusions or 1,4 pairs...
	  
	  // is it in the exclusions
	  if (topo.exclusion(ap.i).count(ap.j)){
	    DEBUG(7, "removing perturbed pair from exclusion " 
		  << ap.i << " and " << ap.j);
	    topo.exclusion(ap.i).erase(ap.j);
	  }
	  if (topo.one_four_pair(ap.i).count(ap.j)){
	    DEBUG(7, "removing perturbed pair from one four " 
		  << ap.i << " and " << ap.j);
	    topo.one_four_pair(ap.i).erase(ap.j);
	  }
	  
	}
      }
      
      if (n != num){
	  io::messages.add("Wrong number of bonds in PERTATOMPAIR03 block.",
			   "In_Perturbation", io::message::error);	
      }
      else if (_lineStream.fail()){
	io::messages.add("Bad line in PERTATOMPAIR03 block.",
			 "In_Perturbation", io::message::error);
      }

      if (!quiet)
	std::cout << "\tEND\n";
      
    } // if block present
  } // PERTATOMPAIR03
  
  { // PERTATOMPARAM
    
    if (buffer.size()){
        
    block_read.insert("PERTATOM03");
    io::messages.add("The PERTATOM03 block was renamed to PERTATOMPARAM.",
             "In_Perturbation", io::message::error); 
    }
    
    buffer = m_block["PERTATOMPARAM"];
    if (buffer.size()){
      block_read.insert("PERTATOMPARAM");

      if (!quiet)
	std::cout << "\tPERTATOMPARAM\n";
      DEBUG(7, "PERTATOMPARAM block");
      
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
      
      int seq, res, a_iac, b_iac;
      double a_mass, b_mass, a_charge, b_charge;
      double lj_soft, crf_soft;
      std::string name;

      if (!quiet)
	std::cout << "\t"
		  << std::setw(5) << "seq"
		  << std::setw(8) << "IAC(A)"
		  << std::setw(10) << "mass(A)"
		  << std::setw(10) << "charge(A)"
		  << std::setw(8) << "IAC(B)"
		  << std::setw(10) << "mass(B)"
		  << std::setw(10) << "charge(B)"
		  << std::setw(10) << "LJ(soft)"
		<< std::setw(10) << "CRF(soft)"
		  << "\n";
      
      for(n = 0; it != buffer.end() - 1; ++it, ++n){
	DEBUG(10, "\treading a line: " << n);
	
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> seq >> res >> name >> a_iac >> a_mass >> a_charge
		    >> b_iac >> b_mass >> b_charge
		    >> lj_soft >> crf_soft;
	
	if (_lineStream.fail()){
	  io::messages.add("Bad line in PERTATOMPARAM block.",
			   "In_Perturbation", io::message::error);
	}
	
	--seq;
	--a_iac;
	--b_iac;

	// weight by input
	lj_soft *= param.perturbation.soft_vdw;
	crf_soft *= param.perturbation.soft_crf;

	if (seq < 0 || seq >= int(topo.num_solute_atoms())){
	  io::messages.add("atom sequence number wrong in PERTATOMPARAM block",
			   "In_Perturbation", io::message::critical);
	  return;
	}

	if (a_iac < 0 || b_iac < 0){
	  std::cout << "n=" << n << "a_iac=" << a_iac << " b_iac=" << b_iac << std::endl;
	  io::messages.add("integer atom code wrong in PERTATOMPARAM block",
			   "In_Perturbation", io::message::critical);
	  return;
	}

	topology::Perturbed_Atom atom(seq, a_iac, a_mass, a_charge,
					b_iac, b_mass, b_charge,
					lj_soft, crf_soft);

	DEBUG(10, "\tcreated an atom");
	
	if (!quiet)
	  std::cout << "\t"
		    << std::setw(5) << seq + 1
		    << std::setw(5) << a_iac + 1
		    << "   "
		    << std::setw(10) << a_mass
		    << std::setw(10) << a_charge
		    << std::setw(5) << b_iac + 1
		    << "   "
		    << std::setw(10) << b_mass
		    << std::setw(10) << b_charge
		    << std::setw(10) << lj_soft
		    << std::setw(10) << crf_soft
		    << "\n";
	
	atom.exclusion() = topo.exclusion(seq);
	topo.exclusion(seq).clear();
	DEBUG(10, "\treplace the exclusions to perturbation");

	std::vector<std::set<int> > & ex = topo.exclusion();
	int seq2=0;
	
	for(std::vector<std::set<int> >::iterator eit=ex.begin(),
	      eto=ex.end(); eit!=eto; ++eit, ++seq2){
	  if(eit->count(seq)){
	    atom.exclusion().insert(seq2);
	    eit->erase(seq);
	  }
	}
	DEBUG(10, "\tadapted perturbed exclusions");
	
	atom.one_four_pair() = topo.one_four_pair(seq);
	topo.one_four_pair(seq).clear();
	DEBUG(10, "\treplaced the 14 interactions");
	
	std::vector<std::set<int> > & ofp = topo.one_four_pair();
	seq2=0;
	
	for(std::vector<std::set<int> >::iterator pit=ofp.begin(), 
	      pito= ofp.end(); pit!=pito; ++pit, ++seq2){
	  if(pit->count(seq)){
	    atom.one_four_pair().insert(seq2);
	    pit->erase(seq);
	  }
	}
	DEBUG(10, "\tadapted 14 interactions");
	
	
	topo.perturbed_solute().atoms()[seq] = atom;

	assert(seq<int(topo.is_perturbed().size()));
	topo.is_perturbed()[seq] = true;
	
      }
      if (n != num){
	  io::messages.add("Wrong number of bonds in PERTATOMPARAM block.",
			   "In_Perturbation", io::message::error);	
      }
      else if (_lineStream.fail()){
	io::messages.add("Bad line in PERTATOMPARAM block.",
			 "In_Perturbation", io::message::error);
      }


      if (!quiet)
	std::cout << "\tEND\n";
      
    } // if block present
    
  } // PERTATOMPARAM
    
  { // SCALEDINTERACTIONS

    buffer = m_block["SCALEDINTERACTIONS"];
    if (buffer.size()){
      block_read.insert("SCALEDINTERACTIONS");

      if(!param.perturbation.scaling){
	io::messages.add("Scaled interactions not turned on, ignoring SCALEDINTERACTIONS block.",
			 "InPerturbation", io::message::warning);
      }
      else{
	if (!quiet)
	  std::cout << "\tSCALED INTERACTIONS\n";
	
	it = buffer.begin() + 1;
	_lineStream.clear();
	_lineStream.str(*it);
	int num, n;
	_lineStream >> num;
	++it;
	
	int i, j;
	double A, B;
	
	if (!quiet)
	  std::cout << "\t"
		    << std::setw(10) << "group i"
		    << std::setw(10) << "group j"
		    << std::setw(10) << "scale A"
		    << std::setw(10) << "scale B"
		    << "\n";
	
	for(n = 0; it != buffer.end() - 1; ++it, ++n){
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> A >> B;
	  
	  if (_lineStream.fail()){
	    io::messages.add("Bad line in SCALEDINTERACTIONS block.",
			     "In_Perturbation", io::message::error);
	  }
	  
	  --i;
	  --j;
	  
	  std::pair<int, int> energy_pair(i,j);
	  std::pair<int, int> energy_pair2(j,i);
	  
	  std::pair<double, double> scale_pair(A,B);
	  
	  topo.energy_group_scaling()[energy_pair]=scale_pair;
	  topo.energy_group_scaling()[energy_pair2]=scale_pair;
	  
	  if (!quiet)
	    std::cout << "\t"
		      << std::setw(10) << i+1
		      << std::setw(10) << j+1
		      << std::setw(10) << A
		      << std::setw(10) << B
		      << std::endl;
	  
	}
	
	if (n != num){
	  io::messages.add("Wrong number of pairs in SCALEDINTERACTIONS block.",
			   "In_Perturbation", io::message::error);	
	}
	else if (_lineStream.fail()){
	  io::messages.add("Bad line in SCALEDINTERACTIONS block.",
			   "In_Perturbation", io::message::error);
	}
	if (!quiet)
	  std::cout << "\tEND\n";
      } // if scaling turned on
      
    } // if block present
    else{
      if(param.perturbation.scaling){
	io::messages.add("Scaling turned on but no SCALEDINTERACTIONS block.",
			 "In_Perturbation", io::message::error);
      }
    }
  } // SCALEDINTERACTIONS

  { // LAMBDADEP
    
    // the standard lambda derivatives
    topo.perturbed_energy_derivative_alpha().push_back(0.0);

    buffer = m_block["LAMBDADEP"];
    if (buffer.size()){
      block_read.insert("LAMBDADEP");

      if(!param.perturbation.scaling){
	io::messages.add("Changed lambda dependence for interactions requires scaling to"
			 " be turned on. Ignoring LAMBDADEP block.",
			 "InPerturbation", io::message::warning);
      }
      else{
	if (!quiet)
	  std::cout << "\tINTERACTIONS with changed lambda dependence\n";
	
	it = buffer.begin() + 1;
	_lineStream.clear();
	_lineStream.str(*it);
	int num, n;
	_lineStream >> num;
	++it;
	
	int i, j;
	double a;
	
	if (!quiet)
	  std::cout << "\t"
		    << std::setw(10) << "index"
		    << std::setw(10) << "group i"
		    << std::setw(10) << "group j"
		    << std::setw(10) << "parameter"
		    << "\n";
	
	for(n = 0; it != buffer.end() - 1; ++it, ++n){
	  _lineStream.clear();
	  _lineStream.str(*it);
	  _lineStream >> i >> j >> a;
	  
	  if (_lineStream.fail()){
	    io::messages.add("Bad line in LAMBDADEP block.",
			     "In_Perturbation", io::message::error);
	  }
	  
	  --i;
	  --j;
	  
	  std::pair<int, int> energy_pair(i,j);
	  std::pair<int, int> energy_pair2(j,i);

	  // check a
	  // this is limited now
	  if (a < -1 || a > 1){
	    /*
	    io::messages.add("LAMBDADEP: parameter a results in l' >1 or <0 "
			     "in the range [0,1] for l (recommended: -1<=a<=1).",
			     "In_Perturbation",
			     io::message::warning);
	    */
	    io::messages.add("LAMBDADEP: limiting l' because -1.0 <= a <= 1.0 not fulfilled",
			    "In_Perturbation",
			    io::message::notice);
	  }
	  
	  std::pair<int, double> lambdadep_pair(n, a);
	  
	  topo.energy_group_lambdadep()[energy_pair]  = lambdadep_pair;
	  topo.energy_group_lambdadep()[energy_pair2] = lambdadep_pair;

	  topo.perturbed_energy_derivative_alpha().push_back(a);
	  
	  if (!quiet)
	    std::cout << "\t"
		      << std::setw(10) << n+1
		      << std::setw(10) << i+1
		      << std::setw(10) << j+1
		      << std::setw(10) << a
		      << std::endl;
	  
	}

	// resize the arrays to cache the lambda primes
	// and the lambda prime / lambda derivative
	topo.lambda_prime().resize(n);
	topo.lambda_prime_derivative().resize(n);

	if (n != num){
	  io::messages.add("Wrong number of pairs in LAMBDADEP block.",
			   "In_Perturbation", io::message::error);	
	}
	else if (_lineStream.fail()){
	  io::messages.add("Bad line in LAMBDADEP block.",
			   "In_Perturbation", io::message::error);
	}
	if (!quiet)
	  std::cout << "\tEND\n";
      } // if scaling turned on
      
    } // if block present
    else{
      // if(param.perturbation.scaling){
      // io::messages.add("Scaling turned on but no LAMBDADEP block.",
      // "In_Perturbation", io::message::warning);
    }
  } // LAMBDADEP

  for(std::map<std::string, std::vector<std::string> >::const_iterator
	it = m_block.begin(),
	to = m_block.end();
      it != to;
      ++it){
    
    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add("block " + it->first + " not supported!",
		       "In_Perturbation (Topology)",
		       io::message::warning);
    }
  }

  if (!quiet)
    std::cout << "END\n";

  // and update the properties for lambda
  topo.update_for_lambda();
  
}

