/**
 * @file eds_nonbonded_innerloop.cc
 * template methods of Eds_Nonbonded_Innerloop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

template<typename t_interaction_spec,
typename t_perturbation_details>
void interaction::Eds_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>
::eds_lj_crf_innerloop
(
        topology::Topology & topo, configuration::Configuration & conf,
        unsigned int i, unsigned int j, Storage &storage,
        Periodicity_type const & periodicity
        ) {
  DEBUG(8, "\teds-perturbed pair\t" << i << "\t" << j << " (inner loop)");

  DEBUG(10, "\t\t        EDS   pert");
  DEBUG(10, "\t\tatom i " << topo.is_eds_perturbed(i) 
	<< " " << topo.is_perturbed(i));
  DEBUG(10, "\t\tatom j " << topo.is_eds_perturbed(j) 
	<< " " << topo.is_perturbed(j));

  // First determine what kind of interaction we have
  // atom i is always eds or perturbed
  // atom j is eds, perturbed or normal
  // gives 6 combinations

  int both_perturbed = 0;
  if(i< topo.num_solute_atoms() && topo.is_perturbed(i)){
    both_perturbed +=3;
  }
  if (j < topo.num_solute_atoms() &&
          topo.is_eds_perturbed(j) == true) {
    both_perturbed += 1;
  }
  if (j < topo.num_solute_atoms() &&
      topo.is_perturbed(j) == true){
    both_perturbed += 2;
  }

  math::Vec r;
  const unsigned int numstates = conf.special().eds.force_endstates.size();
  assert(storage.energies.eds_vi.size() == numstates);
  assert(storage.force_endstates.size() == numstates);

  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);

  const std::vector<std::vector<lj_parameter_struct> > &m_param_lj_parameter = m_param->lj_parameter();

  unsigned int iac_j = topo.iac(j);
  double charge_j = topo.charge()(j);


  double alpha_lj = 0, alpha_crf = 0;


  switch (t_interaction_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      std::vector<math::VArray> &force_endstates = storage.force_endstates;
      double c6 = 0.0, c12 = 0.0, q = 0.0, A_q = 0.0, B_q = 0.0, e_nb = 0.0, e_lj = 0.0,
       e_crf = 0.0, f = 0.0, de_lj = 0.0, de_crf = 0.0;
      assert(abs2(r) != 0);
      const double dist2 = abs2(r);
      const double disti = 1 / abs(r);
      const double dist6 = dist2 * dist2 * dist2;
      switch (both_perturbed) {
        case 0:
	  // EDS - normal
        {
	  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
          std::vector<double> & storage_energies_eds_vi = storage.energies.eds_vi;
          std::vector<math::Matrix> & storage_virial_tensor_endstates = storage.virial_tensor_endstates;
          for (unsigned int state = 0; state < numstates; state++) {
            math::Vec & force_endstates_state_i = force_endstates[state](i);
            math::Vec & force_endstates_state_j = force_endstates[state](j);
            math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
            const lj_parameter_struct &lj = m_param_lj_parameter[(pert_i_M_IAC[state])][iac_j];
            c6 = lj.c6;
            c12 = lj.c12;
            q = pert_i_M_charge[state] * charge_j;

            eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb);

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              force_endstates_state_i(a) += term;
              force_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                virial_tensor_endstates_state(b, a) += r(b) * term;
            }
            // energy
            assert(storage.energies.eds_vi.size() == numstates);
            storage_energies_eds_vi[state] += e_nb;
          }
          break;
        }
        case 1:
	  // EDS - EDS
	{
	  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
          const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
          const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
          std::vector<double> & storage_energies_eds_vi = storage.energies.eds_vi;
          std::vector<math::Matrix> & storage_virial_tensor_endstates = storage.virial_tensor_endstates;
          for (unsigned int state = 0; state < numstates; state++) {
            math::Vec & force_endstates_state_i = force_endstates[state](i);
            math::Vec & force_endstates_state_j = force_endstates[state](j);
            math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
            const lj_parameter_struct &lj =
                    m_param_lj_parameter[(pert_i_M_IAC[state])][(pert_j_M_IAC[state])];
            c6 = lj.c6;
            c12 = lj.c12;
            q = pert_i_M_charge[state]* (pert_j_M_charge[state]);

            // give numstates as reference to const int argument to avoid .size()
            eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb);

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              force_endstates_state_i(a) += term;
              force_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                virial_tensor_endstates_state(b, a) += r(b) * term;
            }
            // energy
            assert(storage.energies.eds_vi.size() == numstates);
            storage_energies_eds_vi[state] += e_nb;
            // }

          }
          break;
        }
        case 2:
	  // EDS - perturbed
	{
	  alpha_lj = topo.perturbed_solute().atoms()[j].LJ_softcore();
	  alpha_crf = topo.perturbed_solute().atoms()[j].CRF_softcore();

	  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
 
          std::vector<double> & storage_energies_eds_vi = storage.energies.eds_vi;
	  std::vector<double> & storage_energies_eds_dvi = storage.perturbed_energy_derivatives.eds_vi;
          std::vector<math::Matrix> & storage_virial_tensor_endstates = storage.virial_tensor_endstates;
          for (unsigned int state = 0; state < numstates; state++) {
            math::Vec & force_endstates_state_i = force_endstates[state](i);
            math::Vec & force_endstates_state_j = force_endstates[state](j);
            math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
            const lj_parameter_struct &A_lj = 
	      m_param_lj_parameter[(pert_i_M_IAC[state])][topo.perturbed_solute().atoms()[j].A_IAC()];
	    const lj_parameter_struct &B_lj = 
	      m_param_lj_parameter[(pert_i_M_IAC[state])][topo.perturbed_solute().atoms()[j].B_IAC()];

            A_q = pert_i_M_charge[state] * topo.perturbed_solute().atoms()[j].A_charge();
	    B_q = pert_i_M_charge[state] * topo.perturbed_solute().atoms()[j].B_charge();
	    

	    int n1 = topo.atom_energy_group(i);
	    int n2 = topo.atom_energy_group(j);

	    set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
		       topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
		       topo.individual_lambda(simulation::crf_lambda)[n1][n2],
		       topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
		       topo.lambda_exp());

            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
					alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              force_endstates_state_i(a) += term;
              force_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                virial_tensor_endstates_state(b, a) += r(b) * term;
            }
            // energy
            assert(storage.energies.eds_vi.size() == numstates);
            storage_energies_eds_vi[state] += e_lj + e_crf;
            assert(storage.perturbed_energy_derivatives.eds_vi.size() == numstates);
	    storage_energies_eds_dvi[state] += de_lj + de_crf;
	    
          }
          break;

        }
       case 3:
	  // perturbed - normal
	{
	  DEBUG(10, "perturbed - normal");
	  

          alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
          alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();
 

          //std::vector<double> & storage_energies_eds_vi = storage.energies.eds_vi;
	  //std::vector<double> & storage_energies_eds_dvi = storage.perturbed_energy_derivatives.eds_vi;
          //std::vector<math::Matrix> & storage_virial_tensor_endstates = storage.virial_tensor_endstates;
          //for (unsigned int state = 0; state < numstates; state++) {
	  //math::Vec & force_endstates_state_i = force_endstates[state](i);
          //  math::Vec & force_endstates_state_j = force_endstates[state](j);
          //  math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
	  const lj_parameter_struct &A_lj = 
	    m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.iac(j)];
	  const lj_parameter_struct &B_lj = 
	    m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.iac(j)];
	  
	  A_q = topo.perturbed_solute().atoms()[i].A_charge() * topo.charge(j);
	  B_q = topo.perturbed_solute().atoms()[i].B_charge() * topo.charge(j);
	    
	  
	  int n1 = topo.atom_energy_group(i);
	  int n2 = topo.atom_energy_group(j);

	  DEBUG(14, "energy groups " << n1 << " " << n2);
	  
	  
	  DEBUG(10, "parameters gathered");
	  
	  set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
		     topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
		     topo.individual_lambda(simulation::crf_lambda)[n1][n2],
		     topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
		     topo.lambda_exp());
	  
	  DEBUG (10, "lambdas set");
	  
	  eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
				      alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
	  
	  DEBUG(10, "interactions computed");
	  
	  // In this case, we can store everything immediately and do not need the endstates
	  
	  DEBUG(10, "\t\tatomic virial");
	  storage.force(i) += f*r;
	  storage.force(j) -= f*r;
	  
	  for (int a = 0; a < 3; ++a) {
	    const double term = f * r(a);
	    for (int b = 0; b < 3; ++b)
	      storage.virial_tensor(b,a) += r(b) * term;
	  }
	  
	  // energy 
	  assert(storage.energies.lj_energy.size() > n1);
	  assert(storage.energies.lj_energy.size() > n2);
	  
	  storage.energies.lj_energy[n1][n2] += e_lj;
	  
	  storage.energies.crf_energy[n1][n2] += e_crf;
	  
	  assert(storage.perturbed_energy_derivatives.lj_energy.size() > n1 &&
		 storage.perturbed_energy_derivatives.lj_energy.size() > n2);
	  
	  assert(storage.perturbed_energy_derivatives.lj_energy[n1].size() > n1 &&
		 storage.perturbed_energy_derivatives.lj_energy[n2].size() > n2);
	  
	  storage.perturbed_energy_derivatives.lj_energy[n1][n2] += de_lj;
	  
	  storage.perturbed_energy_derivatives.crf_energy[n1][n2] += de_crf;
	  
	  DEBUG(10, "forces, energies and energy derivatives stored in storage");
	  
	  break;
	}
       case 4:
	  // perturbed - EDS
	{
	  alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
	  alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

	  const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
	  const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
 
          std::vector<double> & storage_energies_eds_vi = storage.energies.eds_vi;
	  std::vector<double> & storage_energies_eds_dvi = storage.perturbed_energy_derivatives.eds_vi;
          std::vector<math::Matrix> & storage_virial_tensor_endstates = storage.virial_tensor_endstates;
          for (unsigned int state = 0; state < numstates; state++) {
            math::Vec & force_endstates_state_i = force_endstates[state](i);
            math::Vec & force_endstates_state_j = force_endstates[state](j);
            math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
            const lj_parameter_struct &A_lj = 
	      m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][(pert_j_M_IAC[state])];
	    const lj_parameter_struct &B_lj = 
	      m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][(pert_j_M_IAC[state])];

            A_q = topo.perturbed_solute().atoms()[i].A_charge() * pert_j_M_charge[state] ;
	    B_q = topo.perturbed_solute().atoms()[i].B_charge() * pert_j_M_charge[state];
	    

	    int n1 = topo.atom_energy_group(i);
	    int n2 = topo.atom_energy_group(j);

	    set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
		       topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
		       topo.individual_lambda(simulation::crf_lambda)[n1][n2],
		       topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
		       topo.lambda_exp());

            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
					alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              force_endstates_state_i(a) += term;
              force_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                virial_tensor_endstates_state(b, a) += r(b) * term;
            }
            // energy
            assert(storage.energies.eds_vi.size() == numstates);
            storage_energies_eds_vi[state] += e_lj + e_crf;
            assert(storage.perturbed_energy_derivatives.eds_vi.size() == numstates);
	    storage_energies_eds_dvi[state] += de_lj + de_crf;
	    
          }

	  break;
	}
       case 5:
	  // perturbed - perturbed
	{
	  alpha_lj = (topo.perturbed_solute().atoms()[i].LJ_softcore() +
		      topo.perturbed_solute().atoms()[j].LJ_softcore()) /
                  2.0;
          alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
		       topo.perturbed_solute().atoms()[j].CRF_softcore()) /
                  2.0;

	  const lj_parameter_struct &A_lj = 
	    m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.perturbed_solute().atoms()[j].A_IAC()];
	  const lj_parameter_struct &B_lj = 
	    m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.perturbed_solute().atoms()[j].B_IAC()];
	  
	  A_q = topo.perturbed_solute().atoms()[i].A_charge() * topo.perturbed_solute().atoms()[j].A_charge();
	  B_q = topo.perturbed_solute().atoms()[i].B_charge() * topo.perturbed_solute().atoms()[j].B_charge();
	    
	  
	  int n1 = topo.atom_energy_group(i);
	  int n2 = topo.atom_energy_group(j);
	  
	  set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
		     topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
		     topo.individual_lambda(simulation::crf_lambda)[n1][n2],
		     topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
		     topo.lambda_exp());
	  
	  eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
				      alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
	  
	  // In this case, we can store everything immediately and do not need the endstates
	  
	  DEBUG(10, "\t\tatomic virial");
	  storage.force(i) += f*r;
	  storage.force(j) -= f*r;
	  
	  for (int a = 0; a < 3; ++a) {
	    const double term = f * r(a);
	    for (int b = 0; b < 3; ++b)
	      storage.virial_tensor(b,a) += r(b) * term;
	  }
	  
	  // energy 
	  assert(storage.energies.lj_energy.size() > n1);
	  assert(storage.energies.lj_energy.size() > n2);
	  DEBUG(10, "energies lj " << e_lj << " electro " << e_crf);
	  storage.energies.lj_energy[n1][n2] += e_lj;
	  
	  storage.energies.crf_energy[n1][n2] += e_crf;
	  
	  assert(storage.perturbed_energy_derivatives.
		 lj_energy.size() > n1 &&
		 storage.perturbed_energy_derivatives.
		 lj_energy.size() > n2);
	  
	  assert(storage.perturbed_energy_derivatives.
		 lj_energy[n1].size() > n1 &&
		 storage.perturbed_energy_derivatives.
		 lj_energy[n2].size() > n2);
	  
	  storage.perturbed_energy_derivatives.lj_energy[n1][n2] += de_lj;
	  
	  storage.perturbed_energy_derivatives.crf_energy[n1][n2] += de_crf;
	  

	  break;
	}
	
      }
      break;
    }
    case simulation::cggromos_func:
      if (topo.is_coarse_grained(i)) { // CG-CG
        io::messages.add("EDS_Nonbonded_Innerloop",
                "interaction function not implemented for Gromos coarse-grained simulations!",
                io::message::critical);
      } else if (topo.is_coarse_grained(j)) { // FG-CG
        std::vector<math::VArray> &force_endstates = storage.force_endstates;
        double c6 = 0.0, c12 = 0.0, q = 0.0, e_nb = 0.0, f = 0.0;
        assert(abs2(r) != 0);
        const double dist2 = abs2(r);
        const double disti = 1 / abs(r);
        const double dist6 = dist2 * dist2 * dist2;
        switch (both_perturbed) {
          case 0:
          {
	    const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	    const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
            std::vector<double> & storage_energies_eds_vi = storage.energies.eds_vi;
            std::vector<math::Matrix> & storage_virial_tensor_endstates = storage.virial_tensor_endstates;
            for (unsigned int state = 0; state < numstates; state++) {
              math::Vec & force_endstates_state_i = force_endstates[state](i);
              math::Vec & force_endstates_state_j = force_endstates[state](j);
              math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
              const lj_parameter_struct &lj = m_param_lj_parameter[(pert_i_M_IAC[state])][iac_j];
              c6 = lj.c6;
              c12 = lj.c12;
              q = pert_i_M_charge[state] * charge_j;

              eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q / cgrain_eps[1], f, e_nb, 1);

              DEBUG(10, "\t\tatomic virial");
              for (int a = 0; a < 3; ++a) {
                const double term = f * r(a);
                force_endstates_state_i(a) += term;
                force_endstates_state_j(a) -= term;

                for (int b = 0; b < 3; ++b)
                  virial_tensor_endstates_state(b, a) += r(b) * term;
              }
              // energy
              assert(storage.energies.eds_vi.size() == numstates);
              storage_energies_eds_vi[state] += e_nb;
            }
            break;
          }
          case 1:
          {
            io::messages.add("EDS_Nonbonded_Innerloop",
                    "interaction function not implemented for perturbed CG atoms!",
                    io::message::critical);
            break;
          }
        }
      } else { // FG-FG interaction
        std::vector<math::VArray> &force_endstates = storage.force_endstates;
        double c6 = 0.0, c12 = 0.0, q = 0.0, e_nb = 0.0, f = 0.0;
        assert(abs2(r) != 0);
        const double dist2 = abs2(r);
        const double disti = 1 / abs(r);
        const double dist6 = dist2 * dist2 * dist2;
        switch (both_perturbed) {
          case 0:
          {
	    const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	    const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
            std::vector<double> & storage_energies_eds_vi = storage.energies.eds_vi;
            std::vector<math::Matrix> & storage_virial_tensor_endstates = storage.virial_tensor_endstates;
            for (unsigned int state = 0; state < numstates; state++) {
              math::Vec & force_endstates_state_i = force_endstates[state](i);
              math::Vec & force_endstates_state_j = force_endstates[state](j);
              math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
              const lj_parameter_struct &lj = m_param_lj_parameter[(pert_i_M_IAC[state])][iac_j];
              c6 = lj.c6;
              c12 = lj.c12;
              q = pert_i_M_charge[state] * charge_j;

              eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb, 2);

              DEBUG(10, "\t\tatomic virial");
              for (int a = 0; a < 3; ++a) {
                const double term = f * r(a);
                force_endstates_state_i(a) += term;
                force_endstates_state_j(a) -= term;

                for (int b = 0; b < 3; ++b)
                  virial_tensor_endstates_state(b, a) += r(b) * term;
              }
              // energy
              assert(storage.energies.eds_vi.size() == numstates);
              storage_energies_eds_vi[state] += e_nb;
            }
            break;
          }
          case 1:
          {
	    const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	    const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
            const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
            const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
            std::vector<double> & storage_energies_eds_vi = storage.energies.eds_vi;
            std::vector<math::Matrix> & storage_virial_tensor_endstates = storage.virial_tensor_endstates;
            for (unsigned int state = 0; state < numstates; state++) {
              math::Vec & force_endstates_state_i = force_endstates[state](i);
              math::Vec & force_endstates_state_j = force_endstates[state](j);
              math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
              const lj_parameter_struct &lj =
                      m_param_lj_parameter[(pert_i_M_IAC[state])][(pert_j_M_IAC[state])];
              c6 = lj.c6;
              c12 = lj.c12;
              q = pert_i_M_charge[state]* (pert_j_M_charge[state]);

              // give numstates as reference to const int argument to avoid .size()
              eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb, 2);

              DEBUG(10, "\t\tatomic virial");
              for (int a = 0; a < 3; ++a) {
                const double term = f * r(a);
                force_endstates_state_i(a) += term;
                force_endstates_state_j(a) -= term;

                for (int b = 0; b < 3; ++b)
                  virial_tensor_endstates_state(b, a) += r(b) * term;
              }
              // energy
              assert(storage.energies.eds_vi.size() == numstates);
              storage_energies_eds_vi[state] += e_nb;
              // }

            }
            break;
          }
        }
      }
      break;
    case simulation::pol_lj_crf_func:
    case simulation::pol_off_lj_crf_func:
    case simulation::cgrain_func:
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
(topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i, unsigned int j,
        Periodicity_type const & periodicity) {
  DEBUG(8, "\teds one four pair\t" << i << "\t" << j);

  DEBUG(8, "\teds-perturbed pair\t" << i << "\t" << j << " (inner loop)");


  DEBUG(10, "\t\t        EDS   pert");
  DEBUG(10, "\t\tatom i " << topo.is_eds_perturbed(i) 
	<< " " << topo.is_perturbed(i));
  DEBUG(10, "\t\tatom j " << topo.is_eds_perturbed(j) 
	<< " " << topo.is_perturbed(j));

  // First determine what kind of interaction we have
  // atom i is always eds or perturbed
  // atom j is eds, perturbed or normal
  // gives 6 combinations

  int both_perturbed = 0;
  if(i< topo.num_solute_atoms() && topo.is_perturbed(i)){
    both_perturbed +=3;
  }
  if (j < topo.num_solute_atoms() &&
          topo.is_eds_perturbed(j) == true) {
    both_perturbed += 1;
  }
  if (j < topo.num_solute_atoms() &&
      topo.is_perturbed(j) == true){
    both_perturbed += 2;
  }
  DEBUG(10, "\t\tboth_perturbed " << both_perturbed);

  math::Vec r;
  const unsigned int numstates = conf.special().eds.force_endstates.size();
  assert(conf.current().energies.eds_vi.size() == numstates);
  assert(conf.special().eds.force_endstates.size() == numstates);


  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);

  // speed!
  // the following increases topo.eds_perturbed_solute().atoms().size() by 1
  // PROBLEM: j can also be not perturbed!!!
  /*
  const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
  const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
   */

  const std::vector<std::vector<lj_parameter_struct> > &m_param_lj_parameter = m_param->lj_parameter();

  unsigned int iac_j = topo.iac(j);
  double charge_j = topo.charge()(j);
  // --speed

  double alpha_lj = 0, alpha_crf = 0;

  switch (t_interaction_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      double c6 = 0.0, c12 = 0.0, q = 0.0, A_q = 0.0, B_q = 0.0, e_nb = 0.0 , 
      e_lj = 0.0, e_crf = 0.0, f = 0.0, de_lj = 0.0, de_crf = 0.0;
      assert(abs2(r) != 0);
      const double dist2 = abs2(r);
      const double disti = 1 / abs(r);
      const double dist6 = dist2 * dist2 * dist2;
      switch (both_perturbed) {
        case 0:
        {
	  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
          std::vector<math::VArray> & conf_special_eds_force_endstates = conf.special().eds.force_endstates;
          std::vector<math::Matrix> & conf_special_eds_virial_tensor_endstates = conf.special().eds.virial_tensor_endstates;
          std::vector<double> & conf_current_energies_eds_vi = conf.current().energies.eds_vi;
          for (unsigned int state = 0; state < numstates; state++) {
            math::Vec & conf_special_eds_force_endstates_state_i = conf_special_eds_force_endstates[state](i);
            math::Vec & conf_special_eds_force_endstates_state_j = conf_special_eds_force_endstates[state](j);
            math::Matrix & conf_special_eds_virial_tensor_endstates_state = conf_special_eds_virial_tensor_endstates[state];
            const lj_parameter_struct &lj = m_param_lj_parameter[(pert_i_M_IAC[state])][iac_j];
            c6 = lj.cs6;
            c12 = lj.cs12;
            q = pert_i_M_charge[state] * charge_j;

            eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb, 0, m_param->get_coulomb_scaling());

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              conf_special_eds_force_endstates_state_i(a) += term;
              conf_special_eds_force_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                conf_special_eds_virial_tensor_endstates_state(b, a) += r(b) * term;
            }
            // energy
            assert(conf_current_energies_eds_vi.size() == numstates);
            conf_current_energies_eds_vi[state] += e_nb;
          }
          break;
        }
        case 1:
        {
	  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
          const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
          const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
          std::vector<math::VArray> & conf_special_eds_force_endstates = conf.special().eds.force_endstates;
          std::vector<math::Matrix> & conf_special_eds_virial_tensor_endstates = conf.special().eds.virial_tensor_endstates;
          std::vector<double> & conf_current_energies_eds_vi = conf.current().energies.eds_vi;
          for (unsigned int state = 0; state < numstates; state++) {
            math::Vec & conf_special_eds_force_endstates_state_i = conf_special_eds_force_endstates[state](i);
            math::Vec & conf_special_eds_force_endstates_state_j = conf_special_eds_force_endstates[state](j);
            math::Matrix & conf_special_eds_virial_tensor_endstates_state = conf_special_eds_virial_tensor_endstates[state];
            const lj_parameter_struct &lj =
                    m_param_lj_parameter[(pert_i_M_IAC[state])][(pert_j_M_IAC[state])];
            c6 = lj.cs6;
            c12 = lj.cs12;
            q = pert_i_M_charge[state]* (pert_j_M_charge[state]);

            // give numstates as reference to const int argument to avoid .size()
            eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb, 0, m_param->get_coulomb_scaling());

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              conf_special_eds_force_endstates_state_i(a) += term;
              conf_special_eds_force_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                conf_special_eds_virial_tensor_endstates_state(b, a) += r(b) * term;
            }
            // energy
            assert(conf_current_energies_eds_vi.size() == numstates);
            conf_current_energies_eds_vi[state] += e_nb;
            // }

          }
          break;
        }
        case 2:
	  // EDS - perturbed
	{
	  alpha_lj = topo.perturbed_solute().atoms()[j].LJ_softcore();
	  alpha_crf = topo.perturbed_solute().atoms()[j].CRF_softcore();
 
	  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();

	  std::vector<math::VArray> & force_endstates = conf.special().eds.force_endstates;
          std::vector<double> & storage_energies_eds_vi = conf.current().energies.eds_vi;
	  std::vector<double> & storage_energies_eds_dvi = conf.current().perturbed_energy_derivatives.eds_vi;
          std::vector<math::Matrix> & storage_virial_tensor_endstates = conf.special().eds.virial_tensor_endstates;
          for (unsigned int state = 0; state < numstates; state++) {
            math::Vec & force_endstates_state_i = force_endstates[state](i);
            math::Vec & force_endstates_state_j = force_endstates[state](j);
            math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
            const lj_parameter_struct &A_lj = 
	      m_param_lj_parameter[(pert_i_M_IAC[state])][topo.perturbed_solute().atoms()[j].A_IAC()];
	    const lj_parameter_struct &B_lj = 
	      m_param_lj_parameter[(pert_i_M_IAC[state])][topo.perturbed_solute().atoms()[j].B_IAC()];

            A_q = pert_i_M_charge[state] * topo.perturbed_solute().atoms()[j].A_charge();
	    B_q = pert_i_M_charge[state] * topo.perturbed_solute().atoms()[j].B_charge();
	    

	    int n1 = topo.atom_energy_group(i);
	    int n2 = topo.atom_energy_group(j);

	    set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
		       topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
		       topo.individual_lambda(simulation::crf_lambda)[n1][n2],
		       topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
		       topo.lambda_exp());

            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.cs6, A_lj.cs12, B_lj.cs6, B_lj.cs12, A_q, B_q,
					alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              force_endstates_state_i(a) += term;
              force_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                virial_tensor_endstates_state(b, a) += r(b) * term;
            }
            // energy
            assert(conf.current().energies.eds_vi.size() == numstates);
            storage_energies_eds_vi[state] += e_lj+e_crf;
            assert(conf.current().perturbed_energy_derivatives.eds_vi.size() == numstates);
	    storage_energies_eds_dvi[state] += de_lj+de_crf;
	    
          }
          break;

        }

       case 3:
	  // perturbed - normal
	{
          alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
          alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();
 
 
          //std::vector<double> & storage_energies_eds_vi = storage.energies.eds_vi;
	  //std::vector<double> & storage_energies_eds_dvi = storage.perturbed_energy_derivatives.eds_vi;
          //std::vector<math::Matrix> & storage_virial_tensor_endstates = storage.virial_tensor_endstates;
          //for (unsigned int state = 0; state < numstates; state++) {
	  //math::Vec & force_endstates_state_i = force_endstates[state](i);
          //  math::Vec & force_endstates_state_j = force_endstates[state](j);
          //  math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
	  const lj_parameter_struct &A_lj = 
	    m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.iac(j)];
	  const lj_parameter_struct &B_lj = 
	    m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.iac(j)];
	  
	  A_q = topo.perturbed_solute().atoms()[i].A_charge() * topo.charge(j);
	  B_q = topo.perturbed_solute().atoms()[i].B_charge() * topo.charge(j);
	    
	  
	  int n1 = topo.atom_energy_group(i);
	  int n2 = topo.atom_energy_group(j);
	  
	  set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
		     topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
		     topo.individual_lambda(simulation::crf_lambda)[n1][n2],
		     topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
		     topo.lambda_exp());
	  
	  eds_pert_lj_crf_interaction(dist2, dist6, A_lj.cs6, A_lj.cs12, B_lj.cs6, B_lj.cs12, A_q, B_q,
				      alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
	  
	  // In this case, we can store everything immediately and do not need the endstates
	  
	  DEBUG(10, "\t\tatomic virial");
	  conf.current().force(i) += f*r;
	  conf.current().force(j) -= f*r;
	  
	  for (int a = 0; a < 3; ++a) {
	    const double term = f * r(a);
	    for (int b = 0; b < 3; ++b)
	      conf.current().virial_tensor(b,a) += r(b) * term;
	  }
	  
	  // energy 
	  assert(conf.current().energies.lj_energy.size() > n1);
	  assert(conf.current().energies.lj_energy.size() > n2);
	  
	  conf.current().energies.lj_energy[n1][n2] += e_lj;
	  
	  conf.current().energies.crf_energy[n1][n2] += e_crf;
	  
	  assert(conf.current().perturbed_energy_derivatives.
		 lj_energy.size() > n1 &&
		 conf.current().perturbed_energy_derivatives.
		 lj_energy.size() > n2);
	  
	  assert(conf.current().perturbed_energy_derivatives.
		 lj_energy[n1].size() > n1 &&
		 conf.current().perturbed_energy_derivatives.
		 lj_energy[n2].size() > n2);
	  
	  conf.current().perturbed_energy_derivatives.lj_energy[n1][n2] += de_lj;
	  
	  conf.current().perturbed_energy_derivatives.crf_energy[n1][n2] += de_crf;
	  
	  break;
	}
       case 4:
	  // perturbed - EDS
	{
	  alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
	  alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();
 
	  const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
	  const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();

	  std::vector<math::VArray> & force_endstates = conf.special().eds.force_endstates;
          std::vector<double> & storage_energies_eds_vi = conf.current().energies.eds_vi;
	  std::vector<double> & storage_energies_eds_dvi = conf.current().perturbed_energy_derivatives.eds_vi;
          std::vector<math::Matrix> & storage_virial_tensor_endstates = conf.special().eds.virial_tensor_endstates;
          for (unsigned int state = 0; state < numstates; state++) {
            math::Vec & force_endstates_state_i = force_endstates[state](i);
            math::Vec & force_endstates_state_j = force_endstates[state](j);
            math::Matrix & virial_tensor_endstates_state = storage_virial_tensor_endstates[state];
            const lj_parameter_struct &A_lj = 
	      m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][(pert_j_M_IAC[state])];
	    const lj_parameter_struct &B_lj = 
	      m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][(pert_j_M_IAC[state])];

            A_q = topo.perturbed_solute().atoms()[i].A_charge() * pert_j_M_charge[state] ;
	    B_q = topo.perturbed_solute().atoms()[i].B_charge() * pert_j_M_charge[state];
	    

	    int n1 = topo.atom_energy_group(i);
	    int n2 = topo.atom_energy_group(j);

	    set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
		       topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
		       topo.individual_lambda(simulation::crf_lambda)[n1][n2],
		       topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
		       topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
		       topo.lambda_exp());

            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.cs6, A_lj.cs12, B_lj.cs6, B_lj.cs12, A_q, B_q,
					alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              force_endstates_state_i(a) += term;
              force_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                virial_tensor_endstates_state(b, a) += r(b) * term;
            }
            // energy
            assert(conf.current().energies.eds_vi.size() == numstates);
            storage_energies_eds_vi[state] += e_lj+e_crf;
            assert(conf.current().perturbed_energy_derivatives.eds_vi.size() == numstates);
	    storage_energies_eds_dvi[state] += de_lj+de_crf;
	    
          }

	  break;
	}
       case 5:
	  // perturbed - perturbed
	{
	  alpha_lj = (topo.perturbed_solute().atoms()[i].LJ_softcore() +
		      topo.perturbed_solute().atoms()[j].LJ_softcore()) /
                  2.0;
          alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
		       topo.perturbed_solute().atoms()[j].CRF_softcore()) /
                  2.0;

	  const lj_parameter_struct &A_lj = 
	    m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.perturbed_solute().atoms()[j].A_IAC()];
	  const lj_parameter_struct &B_lj = 
	    m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.perturbed_solute().atoms()[j].B_IAC()];
	  
	  A_q = topo.perturbed_solute().atoms()[i].A_charge() * topo.perturbed_solute().atoms()[j].A_charge();
	  B_q = topo.perturbed_solute().atoms()[i].B_charge() * topo.perturbed_solute().atoms()[j].B_charge();
	    
	  
	  int n1 = topo.atom_energy_group(i);
	  int n2 = topo.atom_energy_group(j);
	  
	  set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
		     topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
		     topo.individual_lambda(simulation::crf_lambda)[n1][n2],
		     topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
		     topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
		     topo.lambda_exp());
	  
	  eds_pert_lj_crf_interaction(dist2, dist6, A_lj.cs6, A_lj.cs12, B_lj.cs6, B_lj.cs12, A_q, B_q,
				      alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
	  
	  // In this case, we can store everything immediately and do not need the endstates
	  
	  DEBUG(10, "\t\tatomic virial");
	  conf.current().force(i) += f*r;
	  conf.current().force(j) -= f*r;
	  
	  for (int a = 0; a < 3; ++a) {
	    const double term = f * r(a);
	    for (int b = 0; b < 3; ++b)
	      conf.current().virial_tensor(b,a) += r(b) * term;
	  }
	  
	  // energy 
	  assert(conf.current().energies.lj_energy.size() > n1);
	  assert(conf.current().energies.lj_energy.size() > n2);
	  
	  conf.current().energies.lj_energy[n1][n2] += e_lj;
	  
	  conf.current().energies.crf_energy[n1][n2] += e_crf;
	  
	  assert(conf.current().perturbed_energy_derivatives.
		 lj_energy.size() > n1 &&
		 conf.current().perturbed_energy_derivatives.
		 lj_energy.size() > n2);
	  
	  assert(conf.current().perturbed_energy_derivatives.
		 lj_energy[n1].size() > n1 &&
		 conf.current().perturbed_energy_derivatives.
		 lj_energy[n2].size() > n2);
	  
	  conf.current().perturbed_energy_derivatives.lj_energy[n1][n2] += de_lj;
	  
	  conf.current() .perturbed_energy_derivatives.crf_energy[n1][n2] += de_crf;
	  

	  break;
	}

      }
      break;
    }
    case simulation::cggromos_func:
    {
      if (topo.is_coarse_grained(i)) { // CG-FG or CG-CG
        io::messages.add("EDS_Nonbonded_Innerloop",
                "1,4-interaction function not implemented for Gromos coarse-grained simulations!",
                io::message::critical);
      } else { // FG-FG
        double c6 = 0.0, c12 = 0.0, q = 0.0, e_nb = 0.0, f = 0.0;
        assert(abs2(r) != 0);
        const double dist2 = abs2(r);
        const double disti = 1 / abs(r);
        const double dist6 = dist2 * dist2 * dist2;
        switch (both_perturbed) {
          case 0:
          {
	    const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	    const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
            std::vector<math::VArray> & conf_special_eds_force_endstates = conf.special().eds.force_endstates;
            std::vector<math::Matrix> & conf_special_eds_virial_tensor_endstates = conf.special().eds.virial_tensor_endstates;
            std::vector<double> & conf_current_energies_eds_vi = conf.current().energies.eds_vi;
            for (unsigned int state = 0; state < numstates; state++) {
              math::Vec & conf_special_eds_force_endstates_state_i = conf_special_eds_force_endstates[state](i);
              math::Vec & conf_special_eds_force_endstates_state_j = conf_special_eds_force_endstates[state](j);
              math::Matrix & conf_special_eds_virial_tensor_endstates_state = conf_special_eds_virial_tensor_endstates[state];
              const lj_parameter_struct &lj = m_param_lj_parameter[(pert_i_M_IAC[state])][iac_j];
              c6 = lj.cs6;
              c12 = lj.cs12;
              q = pert_i_M_charge[state] * charge_j;

              eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb, 2);

              DEBUG(10, "\t\tatomic virial");
              for (int a = 0; a < 3; ++a) {
                const double term = f * r(a);
                conf_special_eds_force_endstates_state_i(a) += term;
                conf_special_eds_force_endstates_state_j(a) -= term;

                for (int b = 0; b < 3; ++b)
                  conf_special_eds_virial_tensor_endstates_state(b, a) += r(b) * term;
              }
              // energy
              assert(conf_current_energies_eds_vi.size() == numstates);
              conf_current_energies_eds_vi[state] += e_nb;
            }
            break;
          }
          case 1:
          {
	    const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	    const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
            const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
            const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
            std::vector<math::VArray> & conf_special_eds_force_endstates = conf.special().eds.force_endstates;
            std::vector<math::Matrix> & conf_special_eds_virial_tensor_endstates = conf.special().eds.virial_tensor_endstates;
            std::vector<double> & conf_current_energies_eds_vi = conf.current().energies.eds_vi;
            for (unsigned int state = 0; state < numstates; state++) {
              math::Vec & conf_special_eds_force_endstates_state_i = conf_special_eds_force_endstates[state](i);
              math::Vec & conf_special_eds_force_endstates_state_j = conf_special_eds_force_endstates[state](j);
              math::Matrix & conf_special_eds_virial_tensor_endstates_state = conf_special_eds_virial_tensor_endstates[state];
              const lj_parameter_struct &lj =
                      m_param_lj_parameter[(pert_i_M_IAC[state])][(pert_j_M_IAC[state])];
              c6 = lj.cs6;
              c12 = lj.cs12;
              q = pert_i_M_charge[state]* (pert_j_M_charge[state]);

              // give numstates as reference to const int argument to avoid .size()
              eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb, 2);

              DEBUG(10, "\t\tatomic virial");
              for (int a = 0; a < 3; ++a) {
                const double term = f * r(a);
                conf_special_eds_force_endstates_state_i(a) += term;
                conf_special_eds_force_endstates_state_j(a) -= term;

                for (int b = 0; b < 3; ++b)
                  conf_special_eds_virial_tensor_endstates_state(b, a) += r(b) * term;
              }
              // energy
              assert(conf_current_energies_eds_vi.size() == numstates);
              conf_current_energies_eds_vi[state] += e_nb;
              // }
            }
            break;
          }
        }
      }
      break;
    }
    case simulation::pol_lj_crf_func:
    case simulation::pol_off_lj_crf_func:
    case simulation::cgrain_func:
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
        Periodicity_type const & periodicity) {

  math::VArray &pos = conf.current().pos;
  std::vector<math::VArray> &force = conf.special().eds.force_endstates;
  math::Vec r;
  double e_rf = 0.0;

  //std::set<int>::const_iterator it, to;

  const int i = mit->second.sequence_number();
  // self term has already been calculated for state A, 
  // correct for that and 
  // calculate it for this lambda
  // only a distance independent part
  math::Vec f_rf;
  r = 0.0;

  switch (t_interaction_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      eds_rf_interaction(r, topo.charge()(i) * topo.charge()(i), f_rf, e_rf);
      break;
    }
    case simulation::cggromos_func:
    {
      if (topo.is_coarse_grained(i)) { // CG-CG or CG-FG
        eds_rf_interaction(r, topo.charge()(i) * topo.charge()(i), f_rf, e_rf, 0);
      } else { // FG-FG
        eds_rf_interaction(r, topo.charge()(i) * topo.charge()(i), f_rf, e_rf, 2);
      }
      break;
    }
    default:
      io::messages.add("EDS_Nonbonded_Innerloop",
              "rf excluded interaction function not implemented",
              io::message::critical);
  }
  conf.current().energies.crf_energy[topo.atom_energy_group(i)]
            [topo.atom_energy_group(i)] -= 0.5 * e_rf;

  unsigned int numstates = conf.special().eds.force_endstates.size();
  DEBUG(7, "mit->second.M_charge().size() = " << mit->second.M_charge().size());
  assert(mit->second.M_charge().size() == numstates);

  for (unsigned int state = 0; state < numstates; state++) {
    r = 0.0;
    double q_i = mit->second.M_charge()[state];
    // now calculate everything
    switch (t_interaction_spec::interaction_func) {
      case simulation::lj_crf_func:
      {
        eds_rf_interaction(r, q_i*q_i, f_rf, e_rf);
        break;
      }
      case simulation::cggromos_func:
      {
        if (topo.is_coarse_grained(i)) { // CG-CG or CG-FG
          eds_rf_interaction(r, q_i*q_i / cgrain_eps[0], f_rf, e_rf, 0);
        } else { // FG-FG
          eds_rf_interaction(r, q_i*q_i, f_rf, e_rf, 2);
        }
        break;
      }
      default:
        io::messages.add("EDS_Nonbonded_Innerloop",
                "rf excluded interaction function not implemented",
                io::message::critical);
    }
    DEBUG(7, "Self term for atom " << i << " in state " << state << " = " << e_rf << " (q*q = " << q_i * q_i << ")");
    conf.current().energies.eds_vi[state] += 0.5 * e_rf;
  }
  eds_perturbed_RF_exclusions_loop(topo, conf, i, mit->second.exclusion(), periodicity);
}


template<typename t_interaction_spec,
typename t_perturbation_details>
inline void
interaction::Eds_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>
::perturbed_RF_excluded_interaction_innerloop
(topology::Topology & topo, configuration::Configuration & conf,
        std::map<unsigned int, topology::Perturbed_Atom>::const_iterator const & mit,
        Periodicity_type const & periodicity
) {

  math::VArray &pos = conf.current().pos;
  math::VArray &force = conf.current().force;

  math::Vec r;
  double e_rf, de_rf;

  // self term has already been calculated for state A, 
  // correct for that and 
  // calculate it for this lambda
  // only a distance independent part
  r = 0.0;
  const int i = mit->second.sequence_number();
  const double q_i_a = mit->second.A_charge();
  const double q_i_b = mit->second.B_charge();

  int n1 = topo.atom_energy_group(i);

  // calculate all the lambda values for lj and crf interactions, softness and A and B
  set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n1],
          topo.individual_lambda(simulation::lj_softness_lambda)[n1][n1],
          topo.individual_lambda(simulation::crf_lambda)[n1][n1],
          topo.individual_lambda(simulation::crf_softness_lambda)[n1][n1],
          topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n1],
          topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n1],
          topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n1],
          topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n1],
          topo.lambda_exp());

  // now calculate everything
  // rf_interaction(r, q_i_a*q_i_a, f_old_A, e_crf_old_A);
  switch (t_interaction_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      math::Vec f_rf;
      eds_perturbed_rf_interaction(r, q_i_a*q_i_a, q_i_b * q_i_b,
              mit->second.CRF_softcore(),
              f_rf, e_rf, de_rf);
      break; 
    }
    case simulation::cggromos_func:
    {
      math::Vec f_rf;
      // check if...
      if (topo.is_coarse_grained(i)) { // CG-CG interaction
        eds_perturbed_rf_interaction(r, q_i_a * q_i_a / cgrain_eps[0],
                q_i_b * q_i_b / cgrain_eps[0],
                mit->second.CRF_softcore(),
                f_rf, e_rf, de_rf, 0);
      } else { // FG-FG interaction
        eds_perturbed_rf_interaction(r, q_i_a*q_i_a, q_i_b * q_i_b,
                mit->second.CRF_softcore(),
                f_rf, e_rf, de_rf, 2);
      }
      break;
    }
    default:
      io::messages.add("Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
  }
  DEBUG(7, "Self term for atom " << i << " = "
          << e_rf);

  conf.current().energies.crf_energy[topo.atom_energy_group(i)][topo.atom_energy_group(i)] += 0.5 * e_rf;
  conf.current().perturbed_energy_derivatives.crf_energy[topo.atom_energy_group(i)][topo.atom_energy_group(i)] += 0.5 * de_rf;

  // now loop over the exclusions
  // those are fortunately not in the normal exclusions!
  // this handled by an extra function
  eds_perturbed_RF_exclusions_loop(topo, conf, i, mit->second.exclusion(), periodicity);
}


template<typename t_interaction_spec,
typename t_perturbation_details>
inline void
interaction::Eds_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>
::eds_perturbed_RF_exclusions_loop(topology::Topology & topo, configuration::Configuration & conf, int atom_i,
                                   const topology::excl_cont_t::value_type &exclusions, Periodicity_type const & periodicity){
  // Self term has already been accounted for
  // now loop over all the exclusions
  // those are fortunately not in the normal exclusions!
  math::VArray &pos = conf.current().pos;
  std::vector<math::VArray> &force_states = conf.special().eds.force_endstates;
  math::VArray &force = conf.current().force;
  math::Vec r;

  std::vector<int>::const_iterator it, to;
  it = exclusions.begin();
  to = exclusions.end();

  for (; it != to; ++it) {
      periodicity.nearest_image(pos(atom_i), pos(*it), r);
      DEBUG(8, "r2 i(" << atom_i << "-" << *it << ") " << abs2(r));

      int both_perturbed = 0;
      if(atom_i< topo.num_solute_atoms() && topo.is_perturbed(atom_i)){
        both_perturbed +=3;
      }
      if (*it < topo.num_solute_atoms() &&
              topo.is_eds_perturbed(*it) == true) {
        both_perturbed += 1;
      }
      if (*it < topo.num_solute_atoms() &&
          topo.is_perturbed(*it) == true){
        both_perturbed += 2;
      }

      switch (t_interaction_spec::interaction_func) {
        case simulation::lj_crf_func:
        {
          switch (both_perturbed) {

            case 0:
            {
              // EDS - normal
              unsigned int numstates = conf.special().eds.force_endstates.size();
              math::Vec f_rf;
              double q_i = 0.0, q_j = 0.0, e_rf = 0.0;

              for (unsigned int state = 0; state < numstates; state++) {

                q_i = topo.eds_perturbed_solute().atoms()[atom_i].M_charge()[state];
                q_j = topo.charge()(*it);
                eds_rf_interaction(r, q_i*q_j, f_rf, e_rf);

                DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

                // and add everything to the correct arrays
                conf.current().energies.eds_vi[state] += e_rf;

                force_states[state](atom_i) += f_rf;
                force_states[state](*it) -= f_rf;

                // if (t_interaction_spec::do_virial != math::no_virial){
                for (int a = 0; a < 3; ++a)
                  for (int b = 0; b < 3; ++b)
                    conf.special().eds.virial_tensor_endstates[state](a, b) += r(a) * f_rf(b);

                DEBUG(7, "\tatomic virial done");
                // }
              } // loop over all states
              break;
            }
            case 1:
            {
              // EDS - EDS
              unsigned int numstates = conf.special().eds.force_endstates.size();
              math::Vec f_rf;
              double q_i = 0.0, q_j = 0.0, e_rf = 0.0;
              for (unsigned int state = 0; state < numstates; state++) {

                q_i = topo.eds_perturbed_solute().atoms()[atom_i].M_charge()[state];
                q_j = topo.eds_perturbed_solute().atoms()[*it].M_charge()[state];
                eds_rf_interaction(r, q_i*q_j, f_rf, e_rf);

                DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

                // and add everything to the correct arrays
                conf.current().energies.eds_vi[state] += e_rf;

                force_states[state](atom_i) += f_rf;
                force_states[state](*it) -= f_rf;

                // if (t_interaction_spec::do_virial != math::no_virial){
                for (int a = 0; a < 3; ++a)
                  for (int b = 0; b < 3; ++b)
                    conf.special().eds.virial_tensor_endstates[state](a, b) += r(a) * f_rf(b);

                DEBUG(7, "\tatomic virial done");
                // }
              } // loop over all states
              break;
            }
            case 2:
            {
            // EDS - perturbed
              unsigned int numstates = conf.special().eds.force_endstates.size();
              math::Vec f_rf;
              double A_q = 0.0, B_q = 0.0, e_rf = 0.0, de_rf = 0.0, alpha_crf = 0.0;
	            alpha_crf = topo.perturbed_solute().atoms()[*it].CRF_softcore();

              int n1 = topo.atom_energy_group(atom_i);
	            int n2 = topo.atom_energy_group(*it);

              set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
                    topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
                    topo.individual_lambda(simulation::crf_lambda)[n1][n2],
                    topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
                    topo.lambda_exp());
 
              for (unsigned int state = 0; state < numstates; state++) {

                A_q = topo.eds_perturbed_solute().atoms()[atom_i].M_charge()[state] * topo.perturbed_solute().atoms()[*it].A_charge();
                B_q = topo.eds_perturbed_solute().atoms()[atom_i].M_charge()[state] * topo.perturbed_solute().atoms()[*it].B_charge();

                eds_perturbed_rf_interaction(r, A_q, B_q, alpha_crf, f_rf, e_rf, de_rf);

                DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

                // and add everything to the correct arrays
                conf.current().energies.eds_vi[state] += e_rf;

                force_states[state](atom_i) += f_rf;
                force_states[state](*it) -= f_rf;
                // energy derivatives
                conf.current().perturbed_energy_derivatives.eds_vi[state] += de_rf;
                
                for (int a = 0; a < 3; ++a)
                  for (int b = 0; b < 3; ++b)
                    conf.special().eds.virial_tensor_endstates[state](a, b) += r(a) * f_rf(b);

                DEBUG(7, "\tatomic virial done");
                // }
              } // loop over all states
              break;
            }
            case 3:
            {
              // Perturbed - Normal
              math::Vec f_rf;
              double A_q = 0.0, B_q = 0.0, e_rf = 0.0, de_rf = 0.0, alpha_crf = 0.0;

              alpha_crf = topo.perturbed_solute().atoms()[atom_i].CRF_softcore();

              int n1 = topo.atom_energy_group(atom_i);
	            int n2 = topo.atom_energy_group(*it);

              set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
                    topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
                    topo.individual_lambda(simulation::crf_lambda)[n1][n2],
                    topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
                    topo.lambda_exp());

              A_q = topo.perturbed_solute().atoms()[atom_i].A_charge() * topo.charge()(*it);
              B_q = topo.perturbed_solute().atoms()[atom_i].B_charge() * topo.charge()(*it);

              eds_perturbed_rf_interaction(r, A_q, B_q, alpha_crf, f_rf, e_rf, de_rf);

              DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

              // and add everything to the correct arrays

              conf.current().energies.crf_energy[topo.atom_energy_group(atom_i)][topo.atom_energy_group(*it)] += e_rf;
              conf.current().perturbed_energy_derivatives.crf_energy[topo.atom_energy_group(atom_i)][topo.atom_energy_group(*it)] += de_rf;

              force(atom_i) += f_rf;
              force(*it) -= f_rf;
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    conf.current().virial_tensor(a, b) += r(a) * f_rf(b);

              DEBUG(7, "\tatomic virial done");

              break;
            }
            case 4:
            {
              // Perturbed - EDS
              unsigned int numstates = conf.special().eds.force_endstates.size();
              math::Vec f_rf;
              double A_q = 0.0, B_q = 0.0, e_rf = 0.0, de_rf = 0.0, alpha_crf = 0.0;
	            alpha_crf = topo.perturbed_solute().atoms()[atom_i].CRF_softcore();


              int n1 = topo.atom_energy_group(atom_i);
	            int n2 = topo.atom_energy_group(*it);

              set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
                    topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
                    topo.individual_lambda(simulation::crf_lambda)[n1][n2],
                    topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
                    topo.lambda_exp());
 
              for (unsigned int state = 0; state < numstates; state++) {

                A_q = topo.eds_perturbed_solute().atoms()[*it].M_charge()[state] * topo.perturbed_solute().atoms()[atom_i].A_charge();
                B_q = topo.eds_perturbed_solute().atoms()[*it].M_charge()[state] * topo.perturbed_solute().atoms()[atom_i].B_charge();

                eds_perturbed_rf_interaction(r, A_q, B_q, alpha_crf, f_rf, e_rf, de_rf);

                DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

                // and add everything to the correct arrays
                conf.current().energies.eds_vi[state] += e_rf;

                force_states[state](atom_i) += f_rf;
                force_states[state](*it) -= f_rf;
                // energy derivatives
                conf.current().perturbed_energy_derivatives.eds_vi[state] += de_rf;
                
                for (int a = 0; a < 3; ++a)
                  for (int b = 0; b < 3; ++b)
                    conf.special().eds.virial_tensor_endstates[state](a, b) += r(a) * f_rf(b);

                DEBUG(7, "\tatomic virial done");
                // }
              } // loop over all states
              break;
            }
            case 5:
            {
              // Perturbed - Normal
              math::Vec f_rf;
              double A_q = 0.0, B_q = 0.0, e_rf = 0.0, de_rf = 0.0, alpha_crf = 0.0;

              alpha_crf = (topo.perturbed_solute().atoms()[atom_i].CRF_softcore() + topo.perturbed_solute().atoms()[*it].CRF_softcore())/2.0;

              int n1 = topo.atom_energy_group(atom_i);
	            int n2 = topo.atom_energy_group(*it);

              set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
                    topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
                    topo.individual_lambda(simulation::crf_lambda)[n1][n2],
                    topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
                    topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
                    topo.lambda_exp());

              A_q = topo.perturbed_solute().atoms()[atom_i].A_charge() * topo.perturbed_solute().atoms()[*it].A_charge();
              B_q = topo.perturbed_solute().atoms()[atom_i].B_charge() * topo.perturbed_solute().atoms()[*it].B_charge();

              eds_perturbed_rf_interaction(r, A_q, B_q, alpha_crf, f_rf, e_rf, de_rf);

              DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

              // and add everything to the correct arrays

              conf.current().energies.crf_energy[topo.atom_energy_group(atom_i)][topo.atom_energy_group(*it)] += e_rf;
              conf.current().perturbed_energy_derivatives.crf_energy[topo.atom_energy_group(atom_i)][topo.atom_energy_group(*it)] += de_rf;

              force(atom_i) += f_rf;
              force(*it) -= f_rf;
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    conf.current().virial_tensor(a, b) += r(a) * f_rf(b);

              DEBUG(7, "\tatomic virial done");

              break;
            }
          }
         break;
        }

        case simulation::cggromos_func:
        {
          switch (both_perturbed){
            case 0:
            {
              // EDS - normal
              unsigned int numstates = conf.special().eds.force_endstates.size();
              math::Vec f_rf;
              double q_i = 0.0, q_j = 0.0, e_rf = 0.0;

              for (unsigned int state = 0; state < numstates; state++) {

                q_i = topo.eds_perturbed_solute().atoms()[atom_i].M_charge()[state];
                q_j = topo.charge()(*it);

                if (topo.is_coarse_grained(atom_i) && topo.is_coarse_grained(*it) ) { // CG-CG
                  eds_rf_interaction(r, q_i*q_j / cgrain_eps[0], f_rf, e_rf, 0);
                } else if (topo.is_coarse_grained(atom_i) || topo.is_coarse_grained(*it)) { // FG-CG
                  eds_rf_interaction(r, q_i*q_j / cgrain_eps[1],f_rf, e_rf, 1);
                } else { // FG-FG
                  eds_rf_interaction(r, q_i*q_j, f_rf, e_rf, 2);
                }

                DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

                // and add everything to the correct arrays
                conf.current().energies.eds_vi[state] += e_rf;

                force_states[state](atom_i) += f_rf;
                force_states[state](*it) -= f_rf;

                // if (t_interaction_spec::do_virial != math::no_virial){
                for (int a = 0; a < 3; ++a)
                  for (int b = 0; b < 3; ++b)
                    conf.special().eds.virial_tensor_endstates[state](a, b) += r(a) * f_rf(b);

                DEBUG(7, "\tatomic virial done");
                // }
              } // loop over all states
              break;
            }
            case 1:
            {
              // EDS - EDS
              unsigned int numstates = conf.special().eds.force_endstates.size();
              math::Vec f_rf;
              double q_i = 0.0, q_j = 0.0, e_rf = 0.0;
              for (unsigned int state = 0; state < numstates; state++) {

                q_i = topo.eds_perturbed_solute().atoms()[atom_i].M_charge()[state];
                q_j = topo.eds_perturbed_solute().atoms()[*it].M_charge()[state];


                if (topo.is_coarse_grained(atom_i) && topo.is_coarse_grained(*it) ) { // CG-CG
                  eds_rf_interaction(r, q_i*q_j / cgrain_eps[0], f_rf, e_rf, 0);
                } else if (topo.is_coarse_grained(atom_i) || topo.is_coarse_grained(*it)) { // FG-CG
                  eds_rf_interaction(r, q_i*q_j / cgrain_eps[1], f_rf, e_rf, 1);
                } else { // FG-FG
                  eds_rf_interaction(r, q_i*q_j, f_rf, e_rf, 2);
                }

                DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

                // and add everything to the correct arrays
                conf.current().energies.eds_vi[state] += e_rf;

                force_states[state](atom_i) += f_rf;
                force_states[state](*it) -= f_rf;

                // if (t_interaction_spec::do_virial != math::no_virial){
                for (int a = 0; a < 3; ++a)
                  for (int b = 0; b < 3; ++b)
                    conf.special().eds.virial_tensor_endstates[state](a, b) += r(a) * f_rf(b);

                DEBUG(7, "\tatomic virial done");
                // }
              } // loop over all states
              break;
            }
            default:
                  io::messages.add("EDS_Nonbonded_Innerloop",
                  "rf excluded interaction function not implemented",
                  io::message::critical);
          }
          break;
        }
      default:
      io::messages.add("EDS_Nonbonded_Innerloop",
      "rf excluded interaction function not implemented",
      io::message::critical);
    }
  } // loop over all atoms
}