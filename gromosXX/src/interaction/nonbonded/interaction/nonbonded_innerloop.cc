/**
 * @file nonbonded_innerloop.cc
 * template methods of Nonbonded_Innerloop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::lj_crf_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 unsigned int j,
 Storage & storage,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  DEBUG(8, "\tpair\t" << i << "\t" << j);
  
  math::Vec r;
  double f;
  double e_lj, e_crf;
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);
  DEBUG(10, "\tni r " << r(0) << " / " << r(1) << " / " << r(2));
  
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
      {
	const lj_parameter_struct & lj = 
	  m_param->lj_parameter(topo.iac(i),
				topo.iac(j));
	
	DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
	
	lj_crf_interaction(r, lj.c6, lj.c12,
			   topo.charge(i) * 
			   topo.charge(j),
			   f, e_lj, e_crf);
        
        DEBUG(10, "\t\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
            
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
        }
	break;
      }
    
    case simulation::cgrain_func :
      {
	DEBUG(11, "\tiac(i) = " << topo.iac(i) << " iac(j) = " << topo.iac(j));
	DEBUG(11, "\tcg_parameter.size() = " << m_param->cg_parameter().size());
	DEBUG(11, "\tcg_parameter()[i].size() = " << m_param->cg_parameter()[i].size());

	const lj_parameter_struct & cg = 
	  m_param->cg_parameter(topo.iac(i),
				topo.iac(j));
	
	DEBUG(11, "\tcg-parameter c6=" << cg.c6 << " c12=" << cg.c12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
	
	cgrain_interaction(r, cg.c6, cg.c12,
			   topo.charge(i) * 
			   topo.charge(j),
			   f, e_lj, e_crf);
        
        DEBUG(10, "\t\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
            
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
        }
	break;
      }
    case simulation::pol_lj_crf_func :
      {
        math::Vec rp1, rp2, rpp;
        std::vector<double> f_pol(4, 0.0);
        
        rp1 = r - conf.current().posV(j);
        rp2 = r + conf.current().posV(i);
        rpp = r + conf.current().posV(i) - conf.current().posV(j);
        
        DEBUG(10, "\t rpp " << rpp(0) << " / " << rpp(1) << " / " << rpp(2));
        
        const lj_parameter_struct & lj = m_param->lj_parameter(topo.iac(i), topo.iac(j));
	
	DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        DEBUG(11, "\tcoscharge i=" << topo.coscharge()(i) << " j=" << topo.coscharge()(j));
        
        pol_lj_crf_interaction(r, rp1, rp2, rpp, lj.c6, lj.c12,
                               topo.charge(i), topo.charge(j),
                               topo.coscharge(i), topo.coscharge(j), 
                               f_pol, e_lj, e_crf);
                
        DEBUG(10, "\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f_pol[0]*r(a) + f_pol[1]*rp1(a) + f_pol[2]*rp2(a) + f_pol[3]*rpp(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
          
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b)*r(a)*f_pol[0] + rp1(b)*f_pol[1]*rp1(a) + rp2(b)*f_pol[2]*rp2(a) + rpp(b)*f_pol[3]*rpp(a);
        }

	break;
      }
    
    default:
      io::messages.add("Nonbonded_Innerloop",
		       "interaction function not implemented",
		       io::message::critical);
  }
  
  //////////////////////////////////////////////////
  // molecular virial is calculated by a
  // correction to the atomic virial
  //////////////////////////////////////////////////
  
  /*
  if (fabs(f) > 10000){
    std::cerr << "pair " << i << " - " << j << " force = " << f << std::endl;
    std::cerr << "\tiac " << topo.iac(i) << " " << topo.iac(j)
	      << "\tq " << topo.charge(i) << " " << topo.charge(j)
	      << "\n\tr " << math::v2s(r)
	      << "\n\tr2 " << abs2(r) << std::endl;
  }
  */
  
  // energy
  assert(storage.energies.lj_energy.size() > 
	 topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
	 topo.atom_energy_group(j));

  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));
  
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
  
  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
}
/*
template<typename t_nonbonded_spec>
inline void
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::lj_crf_innerloop_central
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 unsigned int j,
 Storage & storage
 )
{
  DEBUG(8, "\tpair\t" << i << "\t" << j);
  
  double f;
  double e_lj, e_crf;
  
  const math::Vec r = conf.current().pos(i) - conf.current().pos(j);
  
  const lj_parameter_struct &lj = 
    m_param->lj_parameter(topo.iac(i),
			  topo.iac(j));
  
  DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
  DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));

  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
      {
	const lj_parameter_struct & lj = 
	  m_param->lj_parameter(topo.iac(i),
				topo.iac(j));
	
	DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
	
	lj_crf_interaction(r, lj.c6, lj.c12,
			   topo.charge(i) * 
			   topo.charge(j),
			   f, e_lj, e_crf);
	break;
      }
    
    case simulation::cgrain_func :
      {
	const lj_parameter_struct & cg = 
	  m_param->cg_parameter(topo.iac(i),
				topo.iac(j));
	
	DEBUG(11, "\tcg-parameter c6=" << cg.c6 << " c12=" << cg.c12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
	
	cgrain_interaction(r, cg.c6, cg.c12,
			   topo.charge(i) * 
			   topo.charge(j),
			   f, e_lj, e_crf);
	break;
      }
    
    default:
      io::messages.add("Nonbonded_Innerloop",
		       "interaction function not implemented",
		       io::message::critical);
  }

  for (int a=0; a<3; ++a){
    
    const double term = f * r(a);
    storage.force(i)(a) += term;
    storage.force(j)(a) -= term;
    
    // if (t_nonbonded_spec::do_virial != math::no_virial){
    for(int b=0; b<3; ++b){
      storage.virial_tensor(b, a) += 
	r(b) * term;
    }
    // }
  }

  // energy
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));

  assert(storage.energies.lj_energy.size() > 
	 topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
	 topo.atom_energy_group(j));
  
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
  
  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
}
*/
template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::lj_crf_innerloop_shift
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 unsigned int j,
 Storage & storage,
 math::Vec const & shift
 )
{
  DEBUG(8, "\tpair\t" << i << "\t" << j);
  
  double f;
  double e_lj, e_crf;
  
  math::Vec r = conf.current().pos(i) - conf.current().pos(j);
  r -= shift;

  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
      {
	const lj_parameter_struct & lj = 
	  m_param->lj_parameter(topo.iac(i),
				topo.iac(j));
	
	DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
	
	lj_crf_interaction(r, lj.c6, lj.c12,
			   topo.charge(i) * 
			   topo.charge(j),
			   f, e_lj, e_crf);
        
        DEBUG(10, "\t\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
            
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
        }
	break;
      }
    
    case simulation::cgrain_func :
      {
	DEBUG(11, "\tiac(i) = " << topo.iac(i) << " iac(j) = " << topo.iac(j));
	DEBUG(11, "\tcg_parameter.size() = " << m_param->cg_parameter().size());
	DEBUG(11, "\tcg_parameter()[i].size() = " << m_param->cg_parameter()[i].size());

	const lj_parameter_struct & cg = 
	  m_param->cg_parameter(topo.iac(i),
				topo.iac(j));
	
	DEBUG(11, "\tcg-parameter c6=" << cg.c6 << " c12=" << cg.c12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
	
	cgrain_interaction(r, cg.c6, cg.c12,
			   topo.charge(i) * 
			   topo.charge(j),
			   f, e_lj, e_crf);
        
        DEBUG(10, "\t\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f * r(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
            
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b) * term;
        }
	break;
      }
    case simulation::pol_lj_crf_func :
      {
        math::Vec rp1, rp2, rpp;
        std::vector<double> f_pol(4, 0.0);
        
        rp1 = r - conf.current().posV(j);
        rp2 = r + conf.current().posV(i);
        rpp = r + conf.current().posV(i) - conf.current().posV(j);
        
        DEBUG(10, "\t rpp " << rpp(0) << " / " << rpp(1) << " / " << rpp(2));
        
        const lj_parameter_struct & lj = m_param->lj_parameter(topo.iac(i), topo.iac(j));
	
	DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        DEBUG(11, "\tcoscharge i=" << topo.coscharge()(i) << " j=" << topo.coscharge()(j));
        
        pol_lj_crf_interaction(r, rp1, rp2, rpp, lj.c6, lj.c12,
                               topo.charge(i), topo.charge(j),
                               topo.coscharge(i), topo.coscharge(j), 
                               f_pol, e_lj, e_crf);
                
        DEBUG(10, "\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f_pol[0]*r(a) + f_pol[1]*rp1(a) + f_pol[2]*rp2(a) + f_pol[3]*rpp(a);
          storage.force(i)(a) += term;
          storage.force(j)(a) -= term;
          
          for(int b=0; b<3; ++b)
            storage.virial_tensor(b, a) += r(b)*r(a)*f_pol[0] + rp1(b)*f_pol[1]*rp1(a) + rp2(b)*f_pol[2]*rp2(a) + rpp(b)*f_pol[3]*rpp(a);
        }

	break;
      }
    
    default:
      io::messages.add("Nonbonded_Innerloop",
		       "interaction function not implemented",
		       io::message::critical);
  }
  
  // energy
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));

  assert(storage.energies.lj_energy.size() > 
	 topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
	 topo.atom_energy_group(j));
  
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
  
  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
  
}

template<typename t_nonbonded_spec>
void interaction::Nonbonded_Innerloop<t_nonbonded_spec>::one_four_interaction_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 int i,
 int j,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  DEBUG(8, "\t1,4-pair\t" << i << "\t" << j);
  
  math::Vec r;
  double f, e_lj, e_crf;
  
  periodicity.nearest_image(conf.current().pos(i), 
			    conf.current().pos(j), r);

  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
      {
	const lj_parameter_struct & lj = 
	  m_param->lj_parameter(topo.iac(i),
				topo.iac(j));
	
	DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
	
        lj_crf_interaction(r, lj.cs6, lj.cs12,
			 topo.charge()(i) * 
			 topo.charge()(j),
			 f, e_lj, e_crf);
        
        DEBUG(10, "\t\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f * r(a);
          conf.current().force(i)(a) += term;
          conf.current().force(j)(a) -= term;
            
          for(int b=0; b<3; ++b)
            conf.current().virial_tensor(b, a) += r(b) * term;
        }
	break;
      }
    
    case simulation::cgrain_func :
      {
        io::messages.add("Nonbonded_Innerloop",
                         "no 1,4 interactions for coarse-grained simulations!",
		         io::message::critical);
        break;
      }
    case simulation::pol_lj_crf_func :
      {
        math::Vec rp1, rp2, rpp;
        std::vector<double> f_pol(4, 0.0);
        
        rp1 = r - conf.current().posV(j);
        rp2 = r + conf.current().posV(i);
        rpp = r + conf.current().posV(i) - conf.current().posV(j);
        
        DEBUG(10, "\t rpp " << rpp(0) << " / " << rpp(1) << " / " << rpp(2));
        
        const lj_parameter_struct & lj = m_param->lj_parameter(topo.iac(i), topo.iac(j));
	
	DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);
	DEBUG(11, "\tcharge i=" << topo.charge()(i) << " j=" << topo.charge()(j));
        DEBUG(11, "\tcoscharge i=" << topo.coscharge()(i) << " j=" << topo.coscharge()(j));
        
        pol_lj_crf_interaction(r, rp1, rp2, rpp, lj.cs6, lj.cs12,
                               topo.charge(i), topo.charge(j),
                               topo.coscharge(i), topo.coscharge(j), 
                               f_pol, e_lj, e_crf);
                
        DEBUG(10, "\tatomic virial");
        for (int a=0; a<3; ++a){
          const double term = f_pol[0]*r(a) + f_pol[1]*rp1(a) + f_pol[2]*rp2(a) + f_pol[3]*rpp(a);
          conf.current().force(i)(a) += term;
          conf.current().force(j)(a) -= term;
          
          for(int b=0; b<3; ++b)
            conf.current().virial_tensor(b, a) += r(b)*r(a)*f_pol[0] + rp1(b)*f_pol[1]*rp1(a) + rp2(b)*f_pol[2]*rp2(a) + rpp(b)*f_pol[3]*rpp(a);
        }

	break;
      }
    
    default:
      io::messages.add("Nonbonded_Innerloop",
		       "interaction function not implemented",
		       io::message::critical);
  }
  
  // energy
  conf.current().energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
    
  conf.current().energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
    
  DEBUG(11, "\tenergy group i " << topo.atom_energy_group(i)
	<< " j " << topo.atom_energy_group(j));
    
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::RF_excluded_interaction_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 int i,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity)
{
  math::Vec r, f;
  double e_crf;
  
  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;

  std::set<int>::const_iterator it, to;
  it = topo.exclusion(i).begin();
  to = topo.exclusion(i).end();
  
  DEBUG(8, "\tself-term " << i );
  r=0;
  
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
    {
      // this will only contribute in the energy, the force should be zero.
      rf_interaction(r, topo.charge()(i) * topo.charge()(i), f, e_crf);
      conf.current().energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] += 0.5 * e_crf;
      DEBUG(11, "\tcontribution " << 0.5*e_crf);
      
      for( ; it != to; ++it){
        
        DEBUG(11, "\texcluded pair " << i << " - " << *it);
        
        periodicity.nearest_image(pos(i), pos(*it), r);
        
        rf_interaction(r, topo.charge()(i) * topo.charge()(*it), f, e_crf);
        
        force(i) += f;
        force(*it) -= f;
        
        for(int a=0; a<3; ++a)
          for(int b=0; b<3; ++b)
            conf.current().virial_tensor(a, b) += r(a) * f(b);
        DEBUG(11, "\tatomic virial done");
             
        // energy
        conf.current().energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;
        
        DEBUG(11, "\tcontribution " << e_crf);
        
      } // loop over excluded pairs
      break;
    }
    case simulation::cgrain_func :{
      io::messages.add("Nonbonded_Innerloop",
                       "no RF excluded interactions for coarse-grained simulations!",
		       io::message::critical);
      break;
    }
    case simulation::pol_lj_crf_func : {
      math::Vec rp1, rp2, rpp;
      math::VArray f_pol(4);
      f_pol = 0.0;
  
      DEBUG(8, "\tself-term " << i );
      rp1 = -conf.current().posV(*it);
      rp2 = conf.current().posV(i);
      rpp = 0.0;
      
      // this will only contribute in the energy, the force should be zero.
      pol_rf_interaction(r, rp1, rp2, rpp, topo.charge(i), topo.charge(i),
              topo.coscharge(i), topo.coscharge(i), f_pol, e_crf);
      conf.current().energies.crf_energy[topo.atom_energy_group(i)]
              [topo.atom_energy_group(i)] += 0.5 * e_crf;
      DEBUG(11, "\tcontribution " << 0.5*e_crf);
      
      for( ; it != to; ++it){
        
        DEBUG(11, "\texcluded pair " << i << " - " << *it);
        
        periodicity.nearest_image(pos(i), pos(*it), r);
        rp1 = r - conf.current().posV(*it);
        rp2 = r + conf.current().posV(i);
        rpp = r + conf.current().posV(i) - conf.current().posV(*it);
        
        pol_rf_interaction(r, rp1, rp2, rpp, topo.charge()(i),
                topo.charge()(*it), topo.coscharge(i),
                topo.coscharge(*it), f_pol, e_crf);
        
        force(i) += f_pol(0) + f_pol(1) + f_pol(2) + f_pol(3);
        force(*it) -= f_pol(0) + f_pol(1) + f_pol(2) + f_pol(3);
        
        for(int a=0; a<3; ++a)
          for(int b=0; b<3; ++b)
            conf.current().virial_tensor(a, b) +=
                    r(a)*f_pol(0)(b) + rp1(a)*f_pol(1)(b) + rp2(a)*f_pol(2)(b) + rpp(a)*f_pol(3)(b);
        DEBUG(11, "\tatomic virial done");
     
        // energy
        conf.current().energies.crf_energy[topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_crf;
        
        DEBUG(11, "\tcontribution " << e_crf);
        
      } // loop over excluded pairs
      
      break;
    }
    default:
      io::messages.add("Nonbonded_Innerloop",
		       "interaction function not implemented",
		       io::message::critical);
  }
  

}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::RF_solvent_interaction_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 topology::Chargegroup_Iterator const & cg_it,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  math::Vec r;
  double e_crf;
  
  math::VArray &pos   = conf.current().pos;

  // loop over the atoms
  topology::Atom_Iterator at_it = cg_it.begin(),
    at_to = cg_it.end();
  
  switch(t_nonbonded_spec::interaction_func){
    case simulation::lj_crf_func :
    {
      for ( ; at_it != at_to; ++at_it){
        DEBUG(11, "\tsolvent self term " << *at_it);
        // no solvent self term. The distance dependent part and the forces
        // are zero. The distance independent part should add up to zero
        // for the energies and is left out.
        
        for(topology::Atom_Iterator at2_it=at_it+1; at2_it!=at_to; ++at2_it){
          
          DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
          periodicity.nearest_image(pos(*at_it),
                  pos(*at2_it), r);
          
          // for solvent, we don't calculate internal forces (rigid molecules)
          // and the distance independent parts should go to zero
          e_crf = -topo.charge()(*at_it) *
                  topo.charge()(*at2_it) *
                  math::four_pi_eps_i *
                  crf_2cut3i() * abs2(r);
          
          // energy
          conf.current().energies.crf_energy
                  [topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_crf;
        } // loop over at2_it
      } // loop over at_it
      break;
    }
    case simulation::cgrain_func :{
      io::messages.add("Nonbonded_Innerloop",
                       "no RF solvent interaction innerloop for coarse-grained simulations!",
		       io::message::critical);
      break;
    }
    case simulation::pol_lj_crf_func: {
      math::Vec rp1, rp2, rpp;
      
      for ( ; at_it != at_to; ++at_it){
        DEBUG(11, "\tsolvent self term " << *at_it);
        // no solvent self term. The distance dependent part and the forces
        // are zero. The distance independent part should add up to zero
        // for the energies and is left out.
        
        for(topology::Atom_Iterator at2_it=at_it+1; at2_it!=at_to; ++at2_it){
          
          DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
          periodicity.nearest_image(pos(*at_it),
                  pos(*at2_it), r);
          rp1 = r - conf.current().posV(*at2_it);
          rp2 = r + conf.current().posV(*at_it);
          rpp = r + conf.current().posV(*at_it) - conf.current().posV(*at2_it);
          
          // for solvent, we don't calculate internal forces (rigid molecules)
          // and the distance independent parts should go to zero
          const double cqi = topo.coscharge(*at_it);
          const double cqj = topo.coscharge(*at2_it);
          const double qeps = -topo.charge()(*at_it) * topo.charge()(*at2_it) * math::four_pi_eps_i;
          const double qepsp1 = -topo.charge()(*at_it)*(topo.charge()(*at2_it)-cqj)*math::four_pi_eps_i;
          const double qepsp2 = (-topo.charge()(*at_it)-cqi)*topo.charge()(*at2_it)*math::four_pi_eps_i;
          const double qepspp = cqi*cqj*math::four_pi_eps_i;
          e_crf = qeps*crf_2cut3i()*abs2(r) + qepsp1*crf_2cut3i()*abs2(rp1)
          + qepsp2*crf_2cut3i()*abs2(rp2) + qepspp*crf_2cut3i()*abs2(rpp);
          
          // energy
          conf.current().energies.crf_energy
                  [topo.atom_energy_group(*at_it) ]
                  [topo.atom_energy_group(*at2_it)] += e_crf;
        } // loop over at2_it
      } // loop over at_it
      
      break;
    }
    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }
  
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::spc_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 int i,
 int j,
 Storage & storage,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  math::Vec r;
  double x[4], y[4], z[4], r2[4], r2i[4], ri[4], ff[4], tx, ty, tz, fx, fy, fz;
  // , rx, ry, rz;

  // only one energy group
  const int egroup = topo.atom_energy_group(topo.num_solute_atoms());
  DEBUG(8, "\tspc pair\t" << i << "\t" << j << " egroup " << egroup);
  
  double dist6i, e_lj, e_crf, f;
  
  const int ii = i;
  math::Vec const * const pos_i = &conf.current().pos(ii);
  math::Vec * const force_i = &storage.force(ii);

  const int jj = j;
  math::Vec const * const pos_j = &conf.current().pos(jj);
  math::Vec * const force_j = &storage.force(jj);
  DEBUG(9, "ii = " << ii << " jj = " << jj);    

  // O - O
  periodicity.nearest_image(*pos_i, *pos_j, r);
      
  tx = r(0) - (*pos_i)(0) + (*pos_j)(0);
  ty = r(1) - (*pos_i)(1) + (*pos_j)(1);
  tz = r(2) - (*pos_i)(2) + (*pos_j)(2);
      
  assert(abs2(r) != 0);
      
  r2[0] = abs2(r);
  r2i[0] = 1.0 / r2[0];
  dist6i = r2i[0] * r2i[0] * r2i[0];
  ri[0] = sqrt(r2i[0]);
      
  e_lj = (2.634129E-6 * dist6i - 2.617346E-3) * dist6i;
  e_crf = 0.82 * 0.82 * 138.9354 * (ri[0] - m_crf_2cut3i * r2[0] - m_crf_cut);
      
  f = (12 * 2.634129E-6 * dist6i - 6 * 2.617346E-3) * dist6i * r2i[0] +
    0.82 * 0.82 * 138.9354 * (ri[0] * r2i[0] + m_crf_cut3i);
      
  fx = f * r(0);
  fy = f * r(1);
  fz = f * r(2);
      
  (*force_i)(0) += fx;
  (*force_j)(0) -= fx;
  (*force_i)(1) += fy;
  (*force_j)(1) -= fy;
  (*force_i)(2) += fz;
  (*force_j)(2) -= fz;
      
  /*
  rx = r(0) -
    conf.special().rel_mol_com_pos(ii)(0) +
    conf.special().rel_mol_com_pos(jj)(0);
  ry = r(1) -
    conf.special().rel_mol_com_pos(ii)(1) +
    conf.special().rel_mol_com_pos(jj)(1);
  rz = r(2) -
    conf.special().rel_mol_com_pos(ii)(2) +
    conf.special().rel_mol_com_pos(jj)(2);
  */
  
  storage.virial_tensor(0, 0) += r(0) * fx;
  storage.virial_tensor(0, 1) += r(0) * fy;
  storage.virial_tensor(0, 2) += r(0) * fz;
  storage.virial_tensor(1, 0) += r(1) * fx;
  storage.virial_tensor(1, 1) += r(1) * fy;
  storage.virial_tensor(1, 2) += r(1) * fz;
  storage.virial_tensor(2, 0) += r(2) * fx;
  storage.virial_tensor(2, 1) += r(2) * fy;
  storage.virial_tensor(2, 2) += r(2) * fz;
      
  storage.energies.lj_energy[egroup][egroup] += e_lj;
      
  // O - H interactions...
      
  x[0] = (*pos_i)(0) - (*(pos_j+1))(0) + tx;
  y[0] = (*pos_i)(1) - (*(pos_j+1))(1) + ty;
  z[0] = (*pos_i)(2) - (*(pos_j+1))(2) + tz;
      
  x[1] = (*pos_i)(0) - (*(pos_j+2))(0) + tx;
  y[1] = (*pos_i)(1) - (*(pos_j+2))(1) + ty;
  z[1] = (*pos_i)(2) - (*(pos_j+2))(2) + tz;
      
  x[2] = (*(pos_i+1))(0) - (*(pos_j))(0) + tx;
  y[2] = (*(pos_i+1))(1) - (*(pos_j))(1) + ty;
  z[2] = (*(pos_i+1))(2) - (*(pos_j))(2) + tz;
      
  x[3] = (*(pos_i+2))(0) - (*(pos_j))(0) + tx;
  y[3] = (*(pos_i+2))(1) - (*(pos_j))(1) + ty;
  z[3] = (*(pos_i+2))(2) - (*(pos_j))(2) + tz;
      
  r2[0] = x[0]*x[0] + y[0]*y[0] + z[0]*z[0];
  r2[1] = x[1]*x[1] + y[1]*y[1] + z[1]*z[1];
  r2[2] = x[2]*x[2] + y[2]*y[2] + z[2]*z[2];
  r2[3] = x[3]*x[3] + y[3]*y[3] + z[3]*z[3];
      
  r2i[0] = 1.0 / r2[0];
  r2i[1] = 1.0 / r2[1];
  r2i[2] = 1.0 / r2[2];
  r2i[3] = 1.0 / r2[3];
      
  ri[0] = sqrt(r2i[0]);
  ri[1] = sqrt(r2i[1]);
  ri[2] = sqrt(r2i[2]);
  ri[3] = sqrt(r2i[3]);
      
  e_crf -= 0.82 * 0.41 * 138.9354 * (ri[0] + ri[1] + ri[2] + ri[3] -
				    m_crf_2cut3i * (r2[0] + r2[1] + r2[2] + r2[3]) - 4 * m_crf_cut);
      
  ff[0] = -0.82 * 0.41 * 138.9354 * (ri[0] * r2i[0] + m_crf_cut3i);
  ff[1] = -0.82 * 0.41 * 138.9354 * (ri[1] * r2i[1] + m_crf_cut3i);
  ff[2] = -0.82 * 0.41 * 138.9354 * (ri[2] * r2i[2] + m_crf_cut3i);
  ff[3] = -0.82 * 0.41 * 138.9354 * (ri[3] * r2i[3] + m_crf_cut3i);
      
  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
      
  (*force_i)(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+1))(2) -= fz;
      
  /*
  rx = x[0] -
    conf.special().rel_mol_com_pos(ii)(0) +
    conf.special().rel_mol_com_pos(jj+1)(0);
  ry = y[0] - 
    conf.special().rel_mol_com_pos(ii)(1) +
    conf.special().rel_mol_com_pos(jj+1)(1);
  rz = z[0] -
    conf.special().rel_mol_com_pos(ii)(2) +
    conf.special().rel_mol_com_pos(jj+1)(2);
  */    
  
  storage.virial_tensor(0, 0) += x[0] * fx;
  storage.virial_tensor(0, 1) += x[0] * fy;
  storage.virial_tensor(0, 2) += x[0] * fz;
  storage.virial_tensor(1, 0) += y[0] * fx;
  storage.virial_tensor(1, 1) += y[0] * fy;
  storage.virial_tensor(1, 2) += y[0] * fz;
  storage.virial_tensor(2, 0) += z[0] * fx;
  storage.virial_tensor(2, 1) += z[0] * fy;
  storage.virial_tensor(2, 2) += z[0] * fz;
      
  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
      
  (*force_i)(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*force_i)(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*force_i)(2) += fz;
  (*(force_j+2))(2) -= fz;
    
  /*
  rx = x[1] -
    conf.special().rel_mol_com_pos(ii)(0) +
    conf.special().rel_mol_com_pos(jj+2)(0);
  ry = y[1] - 
    conf.special().rel_mol_com_pos(ii)(1) +
    conf.special().rel_mol_com_pos(jj+2)(1);
  rz = z[1] -
    conf.special().rel_mol_com_pos(ii)(2) +
    conf.special().rel_mol_com_pos(jj+2)(2);
  */
  
  storage.virial_tensor(0, 0) += x[1] * fx;
  storage.virial_tensor(0, 1) += x[1] * fy;
  storage.virial_tensor(0, 2) += x[1] * fz;
  storage.virial_tensor(1, 0) += y[1] * fx;
  storage.virial_tensor(1, 1) += y[1] * fy;
  storage.virial_tensor(1, 2) += y[1] * fz;
  storage.virial_tensor(2, 0) += z[1] * fx;
  storage.virial_tensor(2, 1) += z[1] * fy;
  storage.virial_tensor(2, 2) += z[1] * fz;
      
  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
      
  (*(force_i+1))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j))(2) -= fz;
      
  /*
  rx = x[2] -
    conf.special().rel_mol_com_pos(ii+1)(0) +
    conf.special().rel_mol_com_pos(jj)(0);
  ry = y[2] - 
    conf.special().rel_mol_com_pos(ii+1)(1) +
    conf.special().rel_mol_com_pos(jj)(1);
  rz = z[2] -
    conf.special().rel_mol_com_pos(ii+1)(2) +
    conf.special().rel_mol_com_pos(jj)(2);
  */
  
  storage.virial_tensor(0, 0) += x[2] * fx;
  storage.virial_tensor(0, 1) += x[2] * fy;
  storage.virial_tensor(0, 2) += x[2] * fz;
  storage.virial_tensor(1, 0) += y[2] * fx;
  storage.virial_tensor(1, 1) += y[2] * fy;
  storage.virial_tensor(1, 2) += y[2] * fz;
  storage.virial_tensor(2, 0) += z[2] * fx;
  storage.virial_tensor(2, 1) += z[2] * fy;
  storage.virial_tensor(2, 2) += z[2] * fz;
      
  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
      
  (*(force_i+2))(0) += fx;
  (*(force_j))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j))(2) -= fz;

  /*
  rx = x[3] -
    conf.special().rel_mol_com_pos(ii+2)(0) +
    conf.special().rel_mol_com_pos(jj)(0);
  ry = y[3] - 
    conf.special().rel_mol_com_pos(ii+2)(1) +
    conf.special().rel_mol_com_pos(jj)(1);
  rz = z[3] -
    conf.special().rel_mol_com_pos(ii+2)(2) +
    conf.special().rel_mol_com_pos(jj)(2);
  */
  
  storage.virial_tensor(0, 0) += x[3] * fx;
  storage.virial_tensor(0, 1) += x[3] * fy;
  storage.virial_tensor(0, 2) += x[3] * fz;
  storage.virial_tensor(1, 0) += y[3] * fx;
  storage.virial_tensor(1, 1) += y[3] * fy;
  storage.virial_tensor(1, 2) += y[3] * fz;
  storage.virial_tensor(2, 0) += z[3] * fx;
  storage.virial_tensor(2, 1) += z[3] * fy;
  storage.virial_tensor(2, 2) += z[3] * fz;
      
  // H - H interactions...
      
  x[0] = (*(pos_i+1))(0) - (*(pos_j+1))(0) + tx;
  y[0] = (*(pos_i+1))(1) - (*(pos_j+1))(1) + ty;
  z[0] = (*(pos_i+1))(2) - (*(pos_j+1))(2) + tz;

  x[1] = (*(pos_i+1))(0) - (*(pos_j+2))(0) + tx;
  y[1] = (*(pos_i+1))(1) - (*(pos_j+2))(1) + ty;
  z[1] = (*(pos_i+1))(2) - (*(pos_j+2))(2) + tz;

  x[2] = (*(pos_i+2))(0) - (*(pos_j+1))(0) + tx;
  y[2] = (*(pos_i+2))(1) - (*(pos_j+1))(1) + ty;
  z[2] = (*(pos_i+2))(2) - (*(pos_j+1))(2) + tz;

  x[3] = (*(pos_i+2))(0) - (*(pos_j+2))(0) + tx;
  y[3] = (*(pos_i+2))(1) - (*(pos_j+2))(1) + ty;
  z[3] = (*(pos_i+2))(2) - (*(pos_j+2))(2) + tz;

  r2[0] = x[0]*x[0] + y[0]*y[0] + z[0]*z[0];
  r2[1] = x[1]*x[1] + y[1]*y[1] + z[1]*z[1];
  r2[2] = x[2]*x[2] + y[2]*y[2] + z[2]*z[2];
  r2[3] = x[3]*x[3] + y[3]*y[3] + z[3]*z[3];
  
  r2i[0] = 1.0 / r2[0];
  r2i[1] = 1.0 / r2[1];
  r2i[2] = 1.0 / r2[2];
  r2i[3] = 1.0 / r2[3];

  ri[0] = sqrt(r2i[0]);
  ri[1] = sqrt(r2i[1]);
  ri[2] = sqrt(r2i[2]);
  ri[3] = sqrt(r2i[3]);

  e_crf += 0.41 * 0.41 * 138.9354 * (ri[0] + ri[1] + ri[2] + ri[3] -
				    m_crf_2cut3i * (r2[0] + r2[1] + r2[2] + r2[3]) - 4 * m_crf_cut);

  ff[0] = 0.41 * 0.41 * 138.9354 * (ri[0] * r2i[0] + m_crf_cut3i);
  ff[1] = 0.41 * 0.41 * 138.9354 * (ri[1] * r2i[1] + m_crf_cut3i);
  ff[2] = 0.41 * 0.41 * 138.9354 * (ri[2] * r2i[2] + m_crf_cut3i);
  ff[3] = 0.41 * 0.41 * 138.9354 * (ri[3] * r2i[3] + m_crf_cut3i);
  
  fx = ff[0] * x[0];
  fy = ff[0] * y[0];
  fz = ff[0] * z[0];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+1))(2) -= fz;

  /*
  rx = x[0] -
    conf.special().rel_mol_com_pos(ii+1)(0) +
    conf.special().rel_mol_com_pos(jj+1)(0);
  ry = y[0] - 
    conf.special().rel_mol_com_pos(ii+1)(1) +
    conf.special().rel_mol_com_pos(jj+1)(1);
  rz = z[0] -
    conf.special().rel_mol_com_pos(ii+1)(2) +
    conf.special().rel_mol_com_pos(jj+1)(2);
  */

  storage.virial_tensor(0, 0) += x[0] * fx;
  storage.virial_tensor(0, 1) += x[0] * fy;
  storage.virial_tensor(0, 2) += x[0] * fz;
  storage.virial_tensor(1, 0) += y[0] * fx;
  storage.virial_tensor(1, 1) += y[0] * fy;
  storage.virial_tensor(1, 2) += y[0] * fz;
  storage.virial_tensor(2, 0) += z[0] * fx;
  storage.virial_tensor(2, 1) += z[0] * fy;
  storage.virial_tensor(2, 2) += z[0] * fz;

  fx = ff[1] * x[1];
  fy = ff[1] * y[1];
  fz = ff[1] * z[1];
  
  (*(force_i+1))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+1))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+1))(2) += fz;
  (*(force_j+2))(2) -= fz;

  /*
  rx = x[1] -
    conf.special().rel_mol_com_pos(ii+1)(0) +
    conf.special().rel_mol_com_pos(jj+2)(0);
  ry = y[1] - 
    conf.special().rel_mol_com_pos(ii+1)(1) +
    conf.special().rel_mol_com_pos(jj+2)(1);
  rz = z[1] -
    conf.special().rel_mol_com_pos(ii+1)(2) +
    conf.special().rel_mol_com_pos(jj+2)(2);
  */

  storage.virial_tensor(0, 0) += x[1] * fx;
  storage.virial_tensor(0, 1) += x[1] * fy;
  storage.virial_tensor(0, 2) += x[1] * fz;
  storage.virial_tensor(1, 0) += y[1] * fx;
  storage.virial_tensor(1, 1) += y[1] * fy;
  storage.virial_tensor(1, 2) += y[1] * fz;
  storage.virial_tensor(2, 0) += z[1] * fx;
  storage.virial_tensor(2, 1) += z[1] * fy;
  storage.virial_tensor(2, 2) += z[1] * fz;

  fx = ff[2] * x[2];
  fy = ff[2] * y[2];
  fz = ff[2] * z[2];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+1))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+1))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+1))(2) -= fz;

  /*
  rx = x[2] -
    conf.special().rel_mol_com_pos(ii+2)(0) +
    conf.special().rel_mol_com_pos(jj+1)(0);
  ry = y[2] - 
    conf.special().rel_mol_com_pos(ii+2)(1) +
    conf.special().rel_mol_com_pos(jj+1)(1);
  rz = z[2] -
    conf.special().rel_mol_com_pos(ii+2)(2) +
    conf.special().rel_mol_com_pos(jj+1)(2);
  */

  storage.virial_tensor(0, 0) += x[2] * fx;
  storage.virial_tensor(0, 1) += x[2] * fy;
  storage.virial_tensor(0, 2) += x[2] * fz;
  storage.virial_tensor(1, 0) += y[2] * fx;
  storage.virial_tensor(1, 1) += y[2] * fy;
  storage.virial_tensor(1, 2) += y[2] * fz;
  storage.virial_tensor(2, 0) += z[2] * fx;
  storage.virial_tensor(2, 1) += z[2] * fy;
  storage.virial_tensor(2, 2) += z[2] * fz;

  fx = ff[3] * x[3];
  fy = ff[3] * y[3];
  fz = ff[3] * z[3];
  
  (*(force_i+2))(0) += fx;
  (*(force_j+2))(0) -= fx;
  (*(force_i+2))(1) += fy;
  (*(force_j+2))(1) -= fy;
  (*(force_i+2))(2) += fz;
  (*(force_j+2))(2) -= fz;

  /*
  rx = x[3] -
    conf.special().rel_mol_com_pos(ii+2)(0) +
    conf.special().rel_mol_com_pos(jj+2)(0);
  ry = y[3] - 
    conf.special().rel_mol_com_pos(ii+2)(1) +
    conf.special().rel_mol_com_pos(jj+2)(1);
  rz = z[3] -
    conf.special().rel_mol_com_pos(ii+2)(2) +
    conf.special().rel_mol_com_pos(jj+2)(2);
  */

  storage.virial_tensor(0, 0) += x[3] * fx;
  storage.virial_tensor(0, 1) += x[3] * fy;
  storage.virial_tensor(0, 2) += x[3] * fz;
  storage.virial_tensor(1, 0) += y[3] * fx;
  storage.virial_tensor(1, 1) += y[3] * fy;
  storage.virial_tensor(1, 2) += y[3] * fz;
  storage.virial_tensor(2, 0) += z[3] * fx;
  storage.virial_tensor(2, 1) += z[3] * fy;
  storage.virial_tensor(2, 2) += z[3] * fz;
    
  storage.energies.crf_energy[egroup][egroup] += e_crf;
      
}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::electric_field_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i, unsigned int j, math::Vec &e_eli, math::Vec &e_elj,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{

  math::Vec r, rp1, rp2, e_el1, e_el2;

  // energy field term at position i and j
  DEBUG(11, "\tenergy field calculation i: "<<i<<" j: "<<j);

  periodicity.nearest_image(conf.current().pos(i),
                            conf.current().pos(j), r);
  
  switch(t_nonbonded_spec::efield_site) {
    case simulation::ef_atom : {
      rp1 = r - conf.current().posV(j);
      rp2 = -r - conf.current().posV(i);
      
      electric_field_interaction(r, rp1, topo.charge(j),
              topo.coscharge(j), e_el1);
      electric_field_interaction(-r, rp2, topo.charge(i),
              topo.coscharge(i), e_el2); 
      break;
    }
    case simulation::ef_cos : {
      const math::Vec r_cos1 = conf.current().posV(i) + r,
                      r_cos2 = conf.current().posV(j) - r;
      rp1 = r_cos1 - conf.current().posV(j);
      rp2 = r_cos2 - conf.current().posV(i);
      
      electric_field_interaction(r_cos1, rp1, topo.charge(j),
              topo.coscharge(j), e_el1);
      electric_field_interaction(r_cos2, rp2, topo.charge(i),
              topo.coscharge(i), e_el2);
      break;
    }
    default :
      io::messages.add("electric_field_innerloop",
                       "this site for electric field calculation was not implemented.",
		       io::message::critical);
      
  }

  e_eli = e_el1;
  e_elj = e_el2;

}

template<typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_nonbonded_spec>::self_energy_innerloop
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 unsigned int i,
 Storage & storage,
 math::Periodicity<t_nonbonded_spec::boundary_type> const & periodicity
 )
{
  DEBUG(8, "\tself energy of molecule i " << i);

  double self_e;
  const double e_i2 = math::abs2(storage.electric_field(i));
  
  if (t_nonbonded_spec::pol_damping)
    self_energy_interaction(topo.polarizability(i), e_i2, 
                            topo.damping_level(i), topo.damping_power(i), self_e);  
  else 
    self_energy_interaction(topo.polarizability(i), e_i2, self_e);

    

  // self energy term
  DEBUG(10, "\tself energy storage:\t" << self_e);

  storage.energies.self_energy[topo.atom_energy_group(i)] += self_e;

}

