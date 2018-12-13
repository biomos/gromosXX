/**
 * @file multibath.h
 * the multibath parameter class.
 */

#ifndef INCLUDED_MULTIBATH_H
#define INCLUDED_MULTIBATH_H

namespace topology{
  class Topology;
}

namespace configuration
{
  class Energy;
}

namespace simulation
{
  /**
   * @struct bath_struct
   * holds the bath / degree of freedom information
   */
  struct bath_struct
  {
    bath_struct(double t, double tau, double dof, double ir_dof, 
		double com_dof, double solu_constr_dof, 
		double solv_constr_dof, double scale=0, 
		std::vector<double> zeta=std::vector<double>(1, 0.0), double ekin=0)
      : temperature(t),tau(tau), dof(dof), ir_dof(ir_dof),com_dof(com_dof),
	solute_constr_dof(solu_constr_dof), 
	solvent_constr_dof(solv_constr_dof), scale(scale), zeta(zeta), ekin(ekin)
    {}
    
    double temperature;
    double tau;
    double dof;
    double ir_dof;
    double com_dof;
    double solute_constr_dof;
    double solvent_constr_dof;
    double scale;
    std::vector<double> zeta;
    double ekin;
  };

  /**
   * @struct bath_index_struct
   * holds bath index for a range of atoms.
   */
  struct bath_index_struct
  {
    bath_index_struct(unsigned int last_atom, unsigned int last_temp_group, 
		      unsigned int com_bath, unsigned int ir_bath)
      : last_atom(last_atom), last_temperature_group(last_temp_group), 
	com_bath(com_bath), ir_bath(ir_bath){}
    
    unsigned int last_atom;
    unsigned int last_temperature_group;
    unsigned int com_bath;
    unsigned int ir_bath;
  };
  
  /**
   * @class Multibath
   * holds multibath and degree of freedom information.
   */
  class Multibath : public std::vector<bath_struct>
  {
  public:
    
    /**
     * Constructor.
     */
    Multibath(){};

    /**
     * add a bath.
     */
    void add_bath(bath_struct s) { push_back(s);}

    /**
     * add a bath.
     */
    void add_bath(double temperature,
		  double tau = -1, double dof = 0, 
		  double com_dof = 0, double ir_dof = 0,
		  double solute_constr_dof = 0, double solvent_constr_dof = 0)
    {
      push_back(bath_struct(temperature, tau, dof, com_dof, ir_dof, 
			    solute_constr_dof, solvent_constr_dof));
    }
      
    /**
     * add the bath indices for a range of atoms.
     */
    void add_bath_index(bath_index_struct s){ m_bath_index.push_back(s);}

    /**
     * add the bath indices for a range of atoms.
     */
    void add_bath_index(unsigned int last, unsigned int last_m, 
			unsigned int com_bath, unsigned int ir_bath){
      
      m_bath_index.push_back(bath_index_struct(last, last_m, 
					       com_bath, ir_bath));
    }
    
    /**
     * get bath i.
     */
    bath_struct & bath(unsigned int i) {
      assert(i < size());
      return (*this)[i];  
    }

    /**
     * get const bath i.
     */
    bath_struct const & bath(unsigned int i)const{
      assert(i < size());
      return (*this)[i];  
    }

    /**
     * bath indices accessor.
     */
    std::vector<bath_index_struct> & bath_index(){return m_bath_index;}
    
    /**
     * const bath indices accessor.
     */
    std::vector<bath_index_struct> const & bath_index()const{
      return m_bath_index;}

    /**
     * get the bath number of particle number i.
     */
    void in_bath(unsigned int const atom,
		 unsigned int &com, unsigned int &ir)const{
      std::vector<bath_index_struct>::const_iterator 
	it = m_bath_index.begin(),
	to = m_bath_index.end();
  
      for(; it != to; ++it){
	if (it->last_atom >= atom){
	  com = it->com_bath;
	  ir = it->ir_bath;
	  return;
	}
      }
      
      // if no bath read in, the 0 bath is not yet added (change that!)
      // calculate degrees of freedom does the two necessary checks.
      // they should probably be done much earlier
      assert(false);
      
      com = 0;
      ir = 0;  
    }
    
    /**
     * calculate degrees of freedom.
     */
    void calculate_degrees_of_freedom(topology::Topology & topo,
				      bool rottrans_constraints,
                                      bool position_constraints,
                                      double dof_to_subtract);
    
    /**
     * calculate total kinetic energies and temperatures.
     */
    void calc_totals(configuration::Energy const &energy,
		     double & ekin, double & ekin_mol, double & ekin_ir,
		     double & temp, double & temp_mol, double & temp_ir,
		     double & scale)const;

    /**
     * check the state.
     */
    int check_state(unsigned int num_atoms)const;

  private:
    /**
     * the bath index for a range of atoms.
     */
    std::vector<bath_index_struct> m_bath_index;
    
  };
  
} // simulation

#endif
