/**
 * @file OutFlexibleConstraints.tcc
 * implements output for flexible constraints.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE trajectory

#include "../../debug.h"

/**
 * Constructor.
 */
inline
io::OutFlexibleConstraints::OutFlexibleConstraints(std::ostream &os) 
  : m_os(os)
{
};

inline
void io::OutFlexibleConstraints::write_title(std::string title)
{
  m_os << "TITLE\n"
       << title
       << "\nEND\n";
}

/**
 * write out flexible constraints information from a topology.
 */
inline void
io::OutFlexibleConstraints::write_FLEXCON(std::vector<double> &vel,
					  simulation::Topology &topo)
{
  std::vector<simulation::compound::distance_constraint_struct>
    & constr = topo.solute().distance_constraints();
  
  { // FLEXCON
    DEBUG(8, "writing in FLEXCON block");
    m_os << "FLEXCON\n";
    
    std::vector<simulation::compound::distance_constraint_struct>::iterator
      cit = constr.begin(),
      cto = constr.end();

    size_t k = 0;
    m_os.precision(9);
    m_os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
    for( ; cit != cto; ++cit, ++k){
	
      m_os << std::setw(10) << cit->i 
	   << std::setw(10) << cit->j 
	   << std::setw(20) << cit->b0 
	   << std::setw(20) << vel[k]
	   << "\n";
    }

    m_os << "END\n";
    
  } // FLEXCON
  
}

/**
 * write out flexible constraints information from a perturbation topology.
 */
inline void
io::OutFlexibleConstraints
::write_FLEXCON(std::vector<double> &vel,
		simulation::Perturbation_Topology &topo)
{
  std::vector<simulation::compound::distance_constraint_struct>
    & constr = topo.solute().distance_constraints();

  std::vector<simulation::Perturbed_Solute
    ::perturbed_distance_constraint_struct>
    & pert_constr = topo.perturbed_solute().distance_constraints();
  
  { // FLEXCON
    DEBUG(8, "writing in FLEXCON block");
    m_os << "FLEXCON\n";
    
    std::vector<simulation::compound::distance_constraint_struct>::iterator
      cit = constr.begin(),
      cto = constr.end();

    size_t k = 0;
    m_os.precision(9);
    m_os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
    for( ; cit != cto; ++cit, ++k){
	
      m_os << std::setw(10) << cit->i 
	   << std::setw(10) << cit->j 
	   << std::setw(20) << cit->b0 
	   << std::setw(20) << vel[k]
	   << "\n";
    }

    std::vector<simulation::Perturbed_Solute
      ::perturbed_distance_constraint_struct>::iterator
      pcit = pert_constr.begin(),
      pcto = pert_constr.end();

    m_os.precision(9);
    m_os.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
    for( ; pcit != pcto; ++pcit, ++k){
	
      m_os << std::setw(10) << pcit->i 
	   << std::setw(10) << pcit->j 
	   << std::setw(20) << pcit->b0 
	   << std::setw(20) << vel[k]
	   << "\n";
    }

    m_os << "END\n";
    
  } // FLEXCON
  
}

