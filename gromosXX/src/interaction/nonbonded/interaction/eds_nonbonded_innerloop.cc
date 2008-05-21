/**
 * @file eds_nonbonded_innerloop.cc
 * template methods of Eds_Nonbonded_Innerloop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

template<typename t_interaction_spec,
	 typename t_perturbation_details>
void interaction::Eds_Nonbonded_Innerloop<
  t_interaction_spec, t_perturbation_details>
::eds_lj_crf_innerloop
( 
 topology::Topology & topo, configuration::Configuration & conf,
 unsigned int i, unsigned int j, Storage &storage,
 Periodicity_type const & periodicity
 )
{
  DEBUG(8, "\teds-perturbed pair\t" << i << "\t" << j << " (inner loop)");
  
  math::Vec r; 
  const unsigned int numstates = conf.special().eds.force_endstates.size();
  assert(storage.energies.eds_vi.size() == numstates);
  assert(storage.force_endstates.size() == numstates);
  std::vector<double> e_nb(numstates);
  std::vector<double> f(numstates);
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);

  std::vector<double> c6(numstates);
  std::vector<double> c12(numstates);
  std::vector<double> q(numstates);
  
  // speed!
  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
   // the following increases topo.eds_perturbed_solute().atoms().size() by 1
  // PROBLEM: j can also be not perturbed!!!
  /*
  const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
  const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
   */
  const std::vector<std::vector<lj_parameter_struct> > &m_param_lj_parameter=m_param->lj_parameter();
  

  
  unsigned int iac_j=topo.iac(j);
  double charge_j=topo.charge()(j);
  // --speed
  
  int both_perturbed=0;
  if(j < topo.num_solute_atoms() &&
  topo.is_eds_perturbed(j) == true){
    both_perturbed=1;
  }
  
  switch(t_interaction_spec::interaction_func){
    case simulation::lj_crf_func :{
      switch(both_perturbed){
        case 0 : {
          for (unsigned int state=0;state<numstates;state++){
            const lj_parameter_struct &lj =m_param_lj_parameter[(pert_i_M_IAC[state])][iac_j];
            c6[state] = lj.c6;
            c12[state] = lj.c12;
            q[state]=pert_i_M_charge[state]*charge_j;
          }
          break;
        }
        case 1 : {
          const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
          const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
          for (unsigned int state=0;state<numstates;state++){
            const lj_parameter_struct &lj =
            m_param_lj_parameter[(pert_i_M_IAC[state])][(pert_j_M_IAC[state])];
            c6[state] = lj.c6;
            c12[state] = lj.c12;
            q[state] = pert_i_M_charge[state]* (pert_j_M_charge[state]);
          }
          break;
        }
      }
      // give numstates as reference to const int argument to avoid .size()
      eds_lj_crf_interaction(r, c6, c12, q, f, e_nb, numstates);
      
      DEBUG(10, "\t\tatomic virial");
      for (unsigned int state=0;state<numstates;state++){
        for (int a=0; a<3; ++a){
          const double term = f[state] * r(a);
          storage.force_endstates[state](i)(a) += term;
          storage.force_endstates[state](j)(a) -= term;
          
          for(int b=0; b<3; ++b)
            storage.virial_tensor_endstates[state](b, a) += r(b) * term;
        }
        // energy
        assert(storage.energies.eds_vi.size() == numstates);
        storage.energies.eds_vi[state] += e_nb[state];
      }
      break;
    }
    case simulation::pol_lj_crf_func :
    case simulation::cgrain_func :
    default:
      io::messages.add("EDS-Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }
  
  DEBUG(7, "\teds-perturbed lj_crf_innerloop " << i << " - " << j << " done!");
}


template<typename t_interaction_spec,
	 typename t_perturbation_details>
void interaction::Eds_Nonbonded_Innerloop<
  t_interaction_spec, t_perturbation_details>
::eds_one_four_interaction_innerloop
( topology::Topology & topo,
  configuration::Configuration & conf,
  unsigned int i, unsigned int j,
  Periodicity_type const & periodicity)
{
  DEBUG(8, "\teds one four pair\t" << i << "\t" << j);
  
  math::Vec r; 
  const unsigned int numstates = conf.special().eds.force_endstates.size();
  assert(conf.current().energies.eds_vi.size() == numstates);
  assert(conf.special().eds.force_endstates.size() == numstates);
  std::vector<double> e_nb(numstates);
  std::vector<double> f(numstates);
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);

  std::vector<double> c6(numstates);
  std::vector<double> c12(numstates);
  std::vector<double> q(numstates);
  
  // speed!
  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
  /*
  const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
  const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
   */
  const std::vector<std::vector<lj_parameter_struct> > &m_param_lj_parameter=m_param->lj_parameter();
  
  unsigned int iac_j=topo.iac(j);
  double charge_j=topo.charge()(j);
  // --speed
  
  int both_perturbed=0;
  if(j < topo.num_solute_atoms() &&
  topo.is_eds_perturbed(j) == true){
    both_perturbed=1;
  }
  
  switch(t_interaction_spec::interaction_func){
    case simulation::lj_crf_func :{
      switch(both_perturbed){
        case 0 : {
          for (unsigned int state=0;state<numstates;state++){
            const lj_parameter_struct &lj =m_param_lj_parameter[(pert_i_M_IAC[state])][iac_j];
            c6[state] = lj.cs6;
            c12[state] = lj.cs12;
            q[state]=pert_i_M_charge[state]*charge_j;
          }
          break;
        }
        case 1 : {
          const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
          const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
          for (unsigned int state=0;state<numstates;state++){
            const lj_parameter_struct &lj =
            m_param_lj_parameter[(pert_i_M_IAC[state])][(pert_j_M_IAC[state])];
            c6[state] = lj.cs6;
            c12[state] = lj.cs12;
            q[state] = pert_i_M_charge[state]* (pert_j_M_charge[state]);
          }
          break;
        }
      }
      // give numstates as reference to const int argument to avoid .size()
      eds_lj_crf_interaction(r, c6, c12, q, f, e_nb, numstates);
      
      DEBUG(10, "\t\tatomic virial");
      for (unsigned int state=0;state<numstates;state++){
        for (int a=0; a<3; ++a){
          const double term = f[state] * r(a);
          conf.special().eds.force_endstates[state](i)(a) += term;
          conf.special().eds.force_endstates[state](j)(a) -= term;
          
          for(int b=0; b<3; ++b)
            conf.special().eds.virial_tensor_endstates[state](b, a) += r(b) * term;
        }
        // energy
        assert(conf.current().energies.eds_vi.size() == numstates);
        conf.current().energies.eds_vi[state] += e_nb[state];
      }
      break;
    }
    case simulation::pol_lj_crf_func :
    case simulation::cgrain_func :
    default:
      io::messages.add("EDS-Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }
      
  DEBUG(9, "\teds-perturbed one-four lj_crf_innerloop " << i << " - " << j << " done!");
}


template<typename t_interaction_spec,
	 typename t_perturbation_details>
inline void 
interaction::Eds_Nonbonded_Innerloop<
  t_interaction_spec, t_perturbation_details>
::eds_RF_excluded_interaction_innerloop
(topology::Topology & topo, configuration::Configuration & conf,
 std::map<unsigned int, topology::EDS_Perturbed_Atom>::const_iterator const & mit,
 Periodicity_type const & periodicity)
{

  math::VArray &pos   = conf.current().pos;
  std::vector<math::VArray> &force = conf.special().eds.force_endstates;
  math::Vec r;
  double e_rf;

  std::set<int>::const_iterator it, to;
  
  const int i=mit->second.sequence_number();
  // self term has already been calculated for state A, 
  // correct for that and 
  // calculate it for this lambda
  // only a distance independent part
  math::Vec f_rf;
  r=0.0;
  eds_rf_interaction(r, topo.charge()(i) * topo.charge()(i), f_rf, e_rf);
      conf.current().energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] -= 0.5 * e_rf;
  
 
  unsigned int numstates = conf.special().eds.force_endstates.size();
  DEBUG(7,"mit->second.M_charge().size() = " << mit->second.M_charge().size());
  assert(mit->second.M_charge().size() == numstates);
  
  for(unsigned int state = 0; state < numstates; state++){
    r=0.0;
    double q_i = mit->second.M_charge()[state];
    // now calculate everything
    switch(t_interaction_spec::interaction_func){
      case simulation::lj_crf_func :{
        
        eds_rf_interaction(r, q_i*q_i, f_rf, e_rf);
        break;
      }
      default:
        io::messages.add("Nonbonded_Innerloop",
                "interaction function not implemented",
                io::message::critical);
    }
    DEBUG(7, "Self term for atom " << i << " in state " << state << " = "  << e_rf << " (q*q = " << q_i * q_i << ")");
    
    conf.current().energies.eds_vi[state] += 0.5 * e_rf;
    
    // now loop over the exclusions
    // those are fortunately not in the normal exclusions!
    it = mit->second.exclusion().begin();
    to = mit->second.exclusion().end();
    //it = topo.exclusion(i).begin();
    //to = topo.exclusion(i).end();
    
    for( ;it!= to; ++it){
      periodicity.nearest_image(pos(i), pos(*it), r);
      
      DEBUG(8, "r2 i(" << i << "-" << *it << ") " << abs2(r));
      
      double q_j;
            
      if(unsigned(*it) < topo.num_solute_atoms() && topo.is_eds_perturbed(*it)){
        // j perturbed
        q_j = topo.eds_perturbed_solute().atoms()[*it].M_charge()[state];
      }
      else{
        // only i perturbed
        q_j = topo.charge()(*it);
      }

      switch(t_interaction_spec::interaction_func){
        case simulation::lj_crf_func :{
          math::Vec f_rf;
          eds_rf_interaction(r, q_i*q_j, f_rf, e_rf);
          
          DEBUG(7, "excluded atoms " << i << " & " << *it << ": " << e_rf);
          
          // and add everything to the correct arrays
          conf.current().energies.eds_vi[state] += e_rf;
          
          force[state](i) += f_rf;
          force[state](*it) -=f_rf;
          
          // if (t_interaction_spec::do_virial != math::no_virial){
          for(int a=0; a<3; ++a)
            for(int b=0; b<3; ++b)
              conf.special().eds.virial_tensor_endstates[state](a, b) +=
              r(a) * f_rf(b);
          
          DEBUG(7, "\tatomic virial done");
          // }
          break;
        }
        default:
          io::messages.add("Nonbonded_Innerloop",
                  "interaction function not implemented",
                  io::message::critical);
      }
    }
  } // loopover states
}

