/**
 * @file InPerturbationTopology.tcc
 * implements methods of InPerturbationTopology.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

#include "../../debug.h"

/**
 * Constructor.
 */
io::InPerturbationTopology::InPerturbationTopology(std::istream &is) 
  : GInStream(is) 
{
  // read the whole file at beginning
  readStream();
};

inline io::InPerturbationTopology &
io::InPerturbationTopology::operator>>(simulation::Perturbation_Topology &topo)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  std::cout << "PERTURBATION\n";

  { // PERTBOND03
    buffer = m_block["PERTBOND03"];
    if (buffer.size()){
      std::cout << "\tPERTBONDS\n";

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
      
      for(n=0; it != buffer.end() -1; ++it, ++n){
	int i, j, t_A, t_B;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> t_A >> t_B;
	
	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in PERTBOND03 block");
	
	simulation::Bond b(i-1, j-1, t_A-1);
	std::vector<simulation::Bond>::iterator b_it
	  = std::find(topo.solute().bonds().begin(), 
		      topo.solute().bonds().end(), 
		      b);
	
	if (b_it == topo.solute().bonds().end())
	  {
	    throw std::runtime_error("trying to perturb non-existing bond");
	  }
	
	topo.solute().bonds().erase(b_it);
	simulation::Perturbed_Bond pb(b, t_B-1);

	std::cout << std::setw(10) << pb.i+1 
		  << std::setw(10) << pb.j+1
		  << std::setw(10) << pb.type+1 
		  << std::setw(10) << pb.B_type+1 
		  << "\n";
	
	topo.perturbed_solute().bonds().push_back(pb);
      }
      
      if (n != num)
	throw std::runtime_error("error in PERTBOND03 block (n != num)");
      else if (_lineStream.fail())
    	throw std::runtime_error("error in PERTBOND03 block (fail)");

    } // if block present
    
  } // PERTBOND03
  { // PERtBANGLE03
    buffer = m_block["PERTBANGLE03"];
    if (buffer.size()){
      std::cout << "\tPERTANGLES\n";

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
      
      for(n=0; it != buffer.end() -1; ++it, ++n){
	int i, j, k, t_A, t_B;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> t_A >> t_B;
	
	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in PERTBANGLE03 block");
	
	simulation::Angle a(i-1, j-1, k-1, t_A-1);
	std::vector<simulation::Angle>::iterator a_it
	  = std::find(topo.solute().angles().begin(), 
		      topo.solute().angles().end(), 
		      a);
	
	if (a_it == topo.solute().angles().end())
	  {
	    throw std::runtime_error("trying to perturb non-existing angle");
	  }
	
	topo.solute().angles().erase(a_it);
	simulation::Perturbed_Angle pa(a, t_B-1);

	std::cout << std::setw(10) << pa.i+1 
		  << std::setw(10) << pa.j+1
		  << std::setw(10) << pa.k+1
		  << std::setw(10) << pa.type+1 
		  << std::setw(10) << pa.B_type+1 
		  << "\n";
	
	topo.perturbed_solute().angles().push_back(pa);
      }
      
      if (n != num)
	throw std::runtime_error("error in PERTBANGLE03 block (n != num)");
      else if (_lineStream.fail())
    	throw std::runtime_error("error in PERTBANGLE03 block (fail)");

    } // if block present
    
  } // PERTANGLE03
  
  { // PERTATOM03
    
    buffer = m_block["PERTATOM03"];
    if (buffer.size()){
      std::cout << "\tPERTATOMS\n";
      
      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
      
      size_t seq, res, a_iac, b_iac;
      double a_mass, b_mass, a_charge, b_charge;
      double lj_soft, crf_soft;
      std::string name;
      
      for(n = 0; it != buffer.end() - 1; ++it, ++n){
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> seq >> res >> name >> a_iac >> a_mass >> a_charge
		    >> b_iac >> b_mass >> b_charge
		    >> lj_soft >> crf_soft;
	
	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in PERTATOM03 block");
	
	--seq;
	--a_iac;
	--b_iac;
	simulation::Perturbed_Atom atom(seq, a_iac, a_mass, a_charge,
					b_iac, b_mass, b_charge,
					lj_soft, crf_soft);

	topo.perturbed_solute().atoms()[seq] = atom;
	topo.perturbed_atom()[seq] = true;
      }
      
      if (n != num)
	throw std::runtime_error("error in PERTATOM03 block (n != num)");
      else if (_lineStream.fail())
	throw std::runtime_error("error in PERTATOM03 block (fail)");
      
    } // if block present
    
  } // PERTATOM03
  
  std::cout << "END\n";
  
  return *this;
}

