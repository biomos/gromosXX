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
inline
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

  // std::cout << "PERTURBATION\n";
  DEBUG(7, "Reading PERTURBATION");
  
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

      std::cout << "\tEND\n";

    } // if block present
    
  } // PERTBOND03

  { // PERTBANGLE03
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

      std::cout << "\tEND\n";

    } // if block present
  } // PERTANGLE03
  
  { // PERTIMPDIHEDRAL03
    buffer = m_block["PERTIMPDIHEDRAL03"];
    if (buffer.size()){
      std::cout << "\tPERTIMPDIHEDRALS\n";

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
      
      for(n=0; it != buffer.end() -1; ++it, ++n){
	int i, j, k, l, t_A, t_B;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> l >> t_A >> t_B;
	
	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in PERTIMPDIHEDRAL03 block");
	
	simulation::Improper_Dihedral id(i-1, j-1, k-1,l-1, t_A-1);
	std::vector<simulation::Improper_Dihedral>::iterator id_it
	  = std::find(topo.solute().improper_dihedrals().begin(), 
		      topo.solute().improper_dihedrals().end(), 
		      id);
	
	if (id_it == topo.solute().improper_dihedrals().end()){
	  throw std::runtime_error(
		     "trying to perturb non-existing improper dihedral");
	}
	
	topo.solute().improper_dihedrals().erase(id_it);
	simulation::Perturbed_Improper_Dihedral pid(id, t_B-1);

	std::cout << std::setw(10) << pid.i+1 
		  << std::setw(10) << pid.j+1
		  << std::setw(10) << pid.k+1
		  << std::setw(10) << pid.l+1
		  << std::setw(10) << pid.type+1 
		  << std::setw(10) << pid.B_type+1 
		  << "\n";
	
	topo.perturbed_solute().improper_dihedrals().push_back(pid);
      }
      
      if (n != num)
	throw std::runtime_error("error in PERTIMPDIHEDRAL03 block (n != num)");
      else if (_lineStream.fail())
    	throw std::runtime_error("error in PERTIMPDIHEDRAL03 block (fail)");

      std::cout << "\tEND\n";

    } // if block present
   
  } // PERTIMPDIHEDRAL03

  { // PERTDIHEDRAL03
    buffer = m_block["PERTDIHEDRAL03"];
    if (buffer.size()){
      std::cout << "\tPERTDIHEDRALS\n";

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;
      
      for(n=0; it != buffer.end() -1; ++it, ++n){
	int i, j, k, l, t_A, t_B;
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> k >> l >> t_A >> t_B;
	
	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in PERTDIHEDRAL03 block");
	
	simulation::Dihedral id(i-1, j-1, k-1,l-1, t_A-1);
	std::vector<simulation::Dihedral>::iterator id_it
	  = std::find(topo.solute().dihedrals().begin(), 
		      topo.solute().dihedrals().end(), 
		      id);
	
	if (id_it == topo.solute().dihedrals().end()){
	  throw std::runtime_error(
		     "trying to perturb non-existing dihedral dihedral");
	}
	
	topo.solute().dihedrals().erase(id_it);
	simulation::Perturbed_Dihedral pid(id, t_B-1);

	std::cout << std::setw(10) << pid.i+1 
		  << std::setw(10) << pid.j+1
		  << std::setw(10) << pid.k+1
		  << std::setw(10) << pid.l+1
		  << std::setw(10) << pid.type+1 
		  << std::setw(10) << pid.B_type+1 
		  << "\n";
	
	topo.perturbed_solute().dihedrals().push_back(pid);
      }
      
      if (n != num)
	throw std::runtime_error("error in PERTDIHEDRAL03 block (n != num)");
      else if (_lineStream.fail())
    	throw std::runtime_error("error in PERTDIHEDRAL03 block (fail)");

      std::cout << "\tEND\n";

    } // if block present
   
  } // PERTDIHEDRAL03

  { // PERTATOMPAIR03
    // has to be read in before(!!) PERTATOM03
    // because the exclusions and 1,4 exclusions have to be adapted...

    buffer = m_block["PERTATOMPAIR03"];
    if (buffer.size()){
      std::cout << "\tPERTATOMPAIRS\n";

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;

      int i, j, A, B;
      for(n = 0; it != buffer.end() - 1; ++it, ++n){
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> A >> B;
	
	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in PERTATOMPAIR03 block\n\t"+ *it);
	
	--i;
	--j;

	simulation::Perturbed_Atompair ap(i,j,A,B);

	std::cout << std::setw(10) << ap.i+1
		  << std::setw(10) << ap.j+1
		  << std::setw(10) << ap.A_interaction
		  << std::setw(10) << ap.B_interaction
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
	    DEBUG(7, "removing perturbed pair from exclusion " << ap.i << " and " << ap.j);
	    topo.exclusion(ap.i).erase(ap.j);
	  }
	  if (topo.one_four_pair(ap.i).count(ap.j)){
	    DEBUG(7, "removing perturbed pair from one four " << ap.i << " and " << ap.j);
	    topo.one_four_pair(ap.i).erase(ap.j);
	  }
	  
	}
      }
      
      if (n != num)
	throw std::runtime_error("error in PERTATOM03 block (n != num)");
      else if (_lineStream.fail())
	throw std::runtime_error("error in PERTATOM03 block (fail)");
    
      std::cout << "\tEND\n";
      
    } // if block present
  } // PERTATOMPAIR03
  
  { // PERTATOM03
    
    buffer = m_block["PERTATOM03"];
    if (buffer.size()){
      std::cout << "\tPERTATOMS\n";
      DEBUG(7, "PERTATOM03 block");
      
      
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

      std::cout << std::setw(5) << "seq"
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
	
	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in PERTATOM03 block\n"+*it);
	
	--seq;
	--a_iac;
	--b_iac;
	simulation::Perturbed_Atom atom(seq, a_iac, a_mass, a_charge,
					b_iac, b_mass, b_charge,
					lj_soft, crf_soft);

	DEBUG(10, "\tcreated an atom");
	
	std::cout << std::setw(5) << seq + 1
		  << std::setw(8) << a_iac + 1
		  << std::setw(10) << a_mass
		  << std::setw(10) << a_charge
		  << std::setw(8) << b_iac + 1
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

	assert(seq<topo.perturbed_atom().size());
	
	topo.perturbed_atom()[seq] = true;
	
      }
      
      if (n != num)
	throw std::runtime_error("error in PERTATOM03 block (n != num)");
      else if (_lineStream.fail())
	throw std::runtime_error("error in PERTATOM03 block (fail)");
      
      std::cout << "\tEND\n";
      
    } // if block present
    
  } // PERTATOM03
    
  { // SCALEDINTERACTIONS

    buffer = m_block["SCALEDINTERACTIONS"];
    if (buffer.size()){
      std::cout << "\tSCALED INTERACTIONS\n";

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      int num, n;
      _lineStream >> num;
      ++it;

      int i, j;
      double A, B;
      for(n = 0; it != buffer.end() - 1; ++it, ++n){
	_lineStream.clear();
	_lineStream.str(*it);
	_lineStream >> i >> j >> A >> B;
	
	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in SCALEDINTERACTIONS block\n\t"
				   + *it);
	
	--i;
	--j;

	std::pair<int, int> energy_pair(i,j);
	std::pair<int, int> energy_pair2(j,i);
	
	std::pair<double, double> scale_pair(A,B);
	
	topo.energy_group_scaling()[energy_pair]=scale_pair;
	topo.energy_group_scaling()[energy_pair2]=scale_pair;
	
	std::cout << std::setw(10) << i+1
		  << std::setw(10) << j+1
		  << std::setw(10) << A
		  << std::setw(10) << B
		  << std::endl;
	
      }
      
      if (n != num)
	throw std::runtime_error("error in SCALEDINTERACTIONS block (n != num)");
      else if (_lineStream.fail())
	throw std::runtime_error("error in SCALEDINTERACTIONS block (fail)");
    
      std::cout << "\tEND\n";
      
    } // if block present
  } // PERTATOMPAIR03
  
  return *this;
}

