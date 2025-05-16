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
  Periodicity_type const & periodicity,  simulation::Simulation & sim
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

  //MULTIAEDS: todo: if multiaeds; check size accordingly
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
/*    // MULTIAEDS
    case simulation::lj_crf_eds_mult_func:
    {
      //todo: forces!!
      std::map<std::vector<int>, math::VArray> &force_mult_endstates = storage.force_mult_endstates;
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
	      const int site_i = topo.eds_perturbed_solute().atoms()[i].site_number();
	      const int numstates_i = sim.param().eds.multnumstates[site_i]; 
	      const int site_j = sim.param().eds.numsites; // non eds atom, belongs to site R (rest), last site
	      const int state_j = 0; // non eds atom, state is set to 0
	      std::vector<unsigned int> states_i_j {site_i, state_i, site_j, state_j};

          const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
          const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
	        std::map<std::vector<int>, double> & storage_energies_eds_mult_vi = storage.energies.eds_mult_vi;
          //std::vector<double> & storage_energies_eds_mult_vi = storage.energies.eds_mult_vi;
	        std::map<std::vector<int>, math::Matrix > & storage_virial_tensor_mult_endstates = storage.virial_tensor_mult_endstates;
          for (unsigned int state_i = 0; state_i < numstates_i; state_i++) {
            math::Vec & force_mult_endstates_state_i = force_mult_endstates[states_i_j](i);//??
            math::Vec & force_mult_endstates_state_j = force_mult_endstates[states_i_j](j);//??
            math::Matrix & virial_tensor_mult_endstates_state = storage_virial_tensor_mult_endstates[states_i_j];//??
            const lj_parameter_struct &lj = m_param_lj_parameter[(pert_i_M_IAC[state_i])][iac_j];
            c6 = lj.c6;
            c12 = lj.c12;
            q = pert_i_M_charge[state_i] * charge_j;

            eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb);

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              force_mult_endstates_state_i(a) += term;//??
              force_mult_endstates_state_j(a) -= term;//??

              for (int b = 0; b < 3; ++b)
                virial_tensor_mult_endstates_state(b, a) += r(b) * term;//??
            }
            // energy
            //assert(storage.energies.eds_vi.size() == numstates);//remove?
	          storage_energies_eds_mult_vi[states_i_j] += e_nb;
          }
          break;
        }
        case 1:
        // EDS - EDS
        {
          const int site_i = topo.eds_perturbed_solute().atoms()[i].site_number();
          const int numstates_i = sim.param().eds.multnumstates[site_i];
          const int site_j =topo.eds_perturbed_solute().atoms()[j].site_number(); 
	        const int numstates_j = sim.param().eds.multnumstates[site_j];

          const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
          const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
          const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
          const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
	        std::map<std::vector<int>, double> & storage_energies_eds_mult_vi = storage.energies.eds_mult_vi;
          //std::vector<double> & storage_energies_eds_mult_vi = storage.energies.eds_mult_vi;
          std::map<std::vector<int>, math::Matrix> & storage_virial_tensor_mult_endstates = storage.virial_tensor_mult_endstates;

      	  if (site_j >= site_i) {
            for (unsigned int state_i = 0; state_i < numstates_i; state_i++) {
              for (unsigned int state_j = 0; state_j < numstates_j; state_j++) {
		            if ((site_j == site_i) && (state_i != state_j)){
		              continue;
		            }
   		          std::vector<unsigned int> states_i_j {site_i, state_i, site_j, state_j}; 
                math::Vec & force_mult_endstates_state_i = force_mult_endstates[states_i_j](i);//??
                math::Vec & force_mult_endstates_state_j = force_mult_endstates[states_i_j](j);//??
                math::Matrix & virial_tensor_mult_endstates_state = storage_virial_tensor_mult_endstates[states_i_j];//??
                const lj_parameter_struct &lj =
                    m_param_lj_parameter[(pert_i_M_IAC[state_i])][(pert_j_M_IAC[state_j])];
                c6 = lj.c6;
                c12 = lj.c12;
                q = pert_i_M_charge[state_i]* (pert_j_M_charge[state_j]);

                // give numstates as reference to const int argument to avoid .size()
                eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb);

                DEBUG(10, "\t\tatomic virial");
                for (int a = 0; a < 3; ++a) {
                  const double term = f * r(a);
                  force_mult_endstates_state_i(a) += term;
                  force_mult_endstates_state_j(a) -= term;
    
                  for (int b = 0; b < 3; ++b)
                    virial_tensor_mult_endstates_state(b, a) += r(b) * term;
                }
                // energy
                //assert(storage.energies.eds_vi.size() == numstates);
                storage_energies_eds_mult_vi[states_i_j] += e_nb;
              }
            }
          }
          break;
        }
        case 2:
        // EDS - perturbed
        {
          alpha_lj = topo.perturbed_solute().atoms()[j].LJ_softcore();
          alpha_crf = topo.perturbed_solute().atoms()[j].CRF_softcore();

 	  const int site_i = topo.eds_perturbed_solute().atoms()[i].site_number();
          const int numstates_i = sim.param().eds.multnumstates[site_i];
          const int site_j = sim.param().eds.numsites; // non eds atom, belongs to site R (rest), last site
          const int state_j = 0; // non eds atom, state is set to 0

          const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
          const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();

          std::map<std::vector<int>, double> & storage_energies_eds_mult_vi = storage.energies.eds_mult_vi;
          std::map<std::vector<int>, double> & storage_energies_eds_mult_dvi = storage.perturbed_energy_derivatives.eds_mult_vi;//?
          std::map<std::vector<int>, math::Matrix> & storage_virial_tensor_mult_endstates = storage.virial_tensor_mult_endstates;

          for (unsigned int state_i = 0; state_i < numstates_i; state_i++) {
            std::vector<unsigned int> states_i_j {site_i, state_i, site_j, state_j};
            math::Vec & force_mult_endstates_state_i = force_mult_endstates[states_i_j](i);
            math::Vec & force_mult_endstates_state_j = force_mult_endstates[states_i_j](j);
            math::Matrix & virial_tensor_mult_endstates_state = storage_virial_tensor_mult_endstates[states_i_j];
            const lj_parameter_struct &A_lj =
                    m_param_lj_parameter[(pert_i_M_IAC[state_i])][topo.perturbed_solute().atoms()[j].A_IAC()];
            const lj_parameter_struct &B_lj =
                    m_param_lj_parameter[(pert_i_M_IAC[state_i])][topo.perturbed_solute().atoms()[j].B_IAC()];

            A_q = pert_i_M_charge[state_i] * topo.perturbed_solute().atoms()[j].A_charge();
                  B_q = pert_i_M_charge[state_i] * topo.perturbed_solute().atoms()[j].B_charge();

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


            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                                            alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            // Extended TI is not compatible with EDS yet
            if (sim.param().precalclam.nr_lambdas){
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Extended TI-EDS function not implemented",
              io::message::critical);
            }

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              force_mult_endstates_state_i(a) += term;
              force_mult_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                virial_tensor_mult_endstates_state(b, a) += r(b) * term;
            }
            // energy
            //assert(storage.energies.eds_vi.size() == numstates);
            storage_energies_eds_mult_vi[states_i_j] += e_lj + e_crf;
            //assert(storage.perturbed_energy_derivatives.eds_vi.size() == numstates);
            storage_energies_eds_mult_dvi[states_i_j] += de_lj + de_crf;

          }
          break;          if (multiAEDS){
            io::messages.add("EDS-Nonbonded_Innerloop",
            "multiAEDS cannot yet be combined with a perturbation",
            io::message::critical);
          }
        }
        case 3:
	// perturbed - normal
	{
	  DEBUG(10, "perturbed - normal");
          alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
          alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();
 
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

          if (t_perturbation_details::do_scaling) {
            // SCALING ON
            // check whether we need to do scaling
            // based on energy groups
            // Check if we want to scale interactions
            std::pair<int, int> energy_group_pair(topo.atom_energy_group(i),
              topo.atom_energy_group(j));
      
            if (topo.energy_group_scaling().count(energy_group_pair)) {

              // YES, we do scale the interactions!
              lj_crf_scaled_interaction(dist2, dist6, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12,
                A_q, B_q,
                alpha_lj, alpha_crf,
                topo.energy_group_scaling()[energy_group_pair].first,
                topo.energy_group_scaling()[energy_group_pair].second,
                f, e_lj, e_crf, de_lj, de_crf);
            }
            else{  // No scaling
                eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                  alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
            }
          } //END Scaling
          else { // No scaling 
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            //---------------------------------------------------------
            //                     ANITA
            // extended TI: calculate A_e_lj, B_e_lj, A_e_crf, B_e_crf,
            //              A_de_LJ, B_de_lj, A_de_crf, B_de_crf
            //---------------------------------------------------------

            // TODO: could add another parameter, to only calculate every x steps
            // if nr_lambdas > 1, we apply extended TI 
            if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
              DEBUG(8, "precalculate lj_crf_soft");
              double A_e_lj = 0.0, B_e_lj = 0.0, A_e_crf = 0.0, B_e_crf = 0.0,
                A_de_lj = 0.0, B_de_lj = 0.0, A_de_crf = 0.0, B_de_crf = 0.0;

              // determine lambda stepsize from min,max and nr of lambdas
              double lambda_step = (sim.param().precalclam.max_lam - 
                   sim.param().precalclam.min_lam) / 
                   (sim.param().precalclam.nr_lambdas-1);

              //loop over nr_lambdas
              for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){ 

              // determine current lambda for this index
              double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

              // start the calculations
              lj_crf_soft_interaction_ext(r, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12, A_q, B_q, alpha_lj, alpha_crf,
                A_e_lj,  B_e_lj, A_e_crf, B_e_crf,
                A_de_lj, B_de_lj, A_de_crf, B_de_crf,
                lam);

              DEBUG(8, "ANITA: precalculated energies for lambda " << lam
                   << "\n now starting storage");
              DEBUG(8, "\n  A_e_lj " << A_e_lj << "\n  lambda index " << lam_index <<
                   "\n  storage.energies.A_lj_energy.size() " << storage.energies.A_lj_energy.size()
                   << "\n  energy group1 " << topo.atom_energy_group(i) << " energy group2 " 
                   << topo.atom_energy_group(j));
              //            assert(storage.energies.A_lj_energy.size() > lam_index);
              //            assert(storage.energies.A_lj_energy[lam_index].size() > topo.atom_energy_group(i));
              //            assert(storage.energies.A_lj_energy[lam_index][topo.atom_energy_group(i)].size() 
              //                     > topo.atom_energy_group(j));

              storage.energies.A_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_lj;
              storage.energies.B_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_lj;

              storage.energies.A_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_crf;
              storage.energies.B_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_crf;

              storage.perturbed_energy_derivatives.A_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_lj;
              storage.perturbed_energy_derivatives.B_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_lj;

              storage.perturbed_energy_derivatives.A_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_crf;
              storage.perturbed_energy_derivatives.B_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_crf;
              DEBUG(8, "\ndone with storing energies ");
              } //all 101 lambda points done
            } // done with extended TI
          }
	  
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

          const int site_i = sim.param().eds.numsites; // non eds atom, belongs to site R (rest), last site
          const int state_i = 0; // non eds atom, state is set to 0
          const int site_j = topo.eds_perturbed_solute().atoms()[j].site_number();
          const int numstates_j = sim.param().eds.multnumstates[site_j];

          const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
          const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();

          std::map<std::vector<int>, double> & storage_energies_eds_mult_vi = storage.energies.eds_mult_vi;
          std::map<std::vector<int>, double> & storage_energies_eds_mult_dvi = storage.perturbed_energy_derivatives.eds_mult_vi;
          std::map<std::vector<int>, math::Matrix> & storage_virial_tensor_mult_endstates = storage.virial_tensor_mult_endstates;

          for (unsigned int state_j = 0; state_j < numstates_j; state_j++) {
            math::Vec & force_mult_endstates_state_i = force_endstates[states_i_j](i);
            math::Vec & force_mult_endstates_state_j = force_endstates[states_i_j](j);
            math::Matrix & virial_tensor_mult_endstates_state = storage_virial_tensor_mult_endstates[states_i_j];
            const lj_parameter_struct &A_lj =
              m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][(pert_j_M_IAC[state])];
            const lj_parameter_struct &B_lj =
              m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][(pert_j_M_IAC[state])];

            A_q = topo.perturbed_solute().atoms()[i].A_charge() * pert_j_M_charge[state_j] ;
            B_q = topo.perturbed_solute().atoms()[i].B_charge() * pert_j_M_charge[state_j];

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

	    // Scaling is not compatible with EDS
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }

            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
            alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            // Extended TI is not compatible with EDS yet
            if (sim.param().precalclam.nr_lambdas){
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Extended TI-EDS function not implemented",
              io::message::critical);
            }

            DEBUG(10, "\t\tatomic virial");
            for (int a = 0; a < 3; ++a) {
              const double term = f * r(a);
              force_mult_endstates_state_i(a) += term;
              force_mult_endstates_state_j(a) -= term;

              for (int b = 0; b < 3; ++b)
                virial_tensor_mult_endstates_state(b, a) += r(b) * term;
            }
            // energy
            //assert(storage.energies.eds_vi.size() == numstates);
            storage_energies_eds_mult_vi[states_i_j] += e_lj + e_crf;
            //assert(storage.perturbed_energy_derivatives.eds_vi.size() == numstates);
            storage_energies_eds_mult_dvi[states_i_j] += de_lj + de_crf;

          }
          break;
        }

        case 5:
	// perturbed - perturbed
	{
	  alpha_lj = (topo.perturbed_solute().atoms()[i].LJ_softcore() +
		topo.perturbed_solute().atoms()[j].LJ_softcore()) / 2.0;
          alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
            	topo.perturbed_solute().atoms()[j].CRF_softcore()) / 2.0;

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
	  
          if (t_perturbation_details::do_scaling) {
            // SCALING ON
            // check whether we need to do scaling
            // based on energy groups
            // Check if we want to scale interactions
            std::pair<int, int> energy_group_pair(topo.atom_energy_group(i),
              topo.atom_energy_group(j));
      
            if (topo.energy_group_scaling().count(energy_group_pair)) {

              // YES, we do scale the interactions!
              lj_crf_scaled_interaction(dist2, dist6, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12,
                A_q, B_q,
                alpha_lj, alpha_crf,
                topo.energy_group_scaling()[energy_group_pair].first,
                topo.energy_group_scaling()[energy_group_pair].second,
                f, e_lj, e_crf, de_lj, de_crf);

            }
            else{  // No scaling
      
              eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
            }
          } //END Scaling
          else { // No scaling 
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            //---------------------------------------------------------
            //                     ANITA
            // extended TI: calculate A_e_lj, B_e_lj, A_e_crf, B_e_crf,
            //              A_de_LJ, B_de_lj, A_de_crf, B_de_crf
            //---------------------------------------------------------

            // TODO: could add another parameter, to only calculate every x steps
            // if nr_lambdas > 1, we apply extended TI 
            if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
              DEBUG(8, "precalculate lj_crf_soft");
              double A_e_lj = 0.0, B_e_lj = 0.0, A_e_crf = 0.0, B_e_crf = 0.0,
              A_de_lj = 0.0, B_de_lj = 0.0, A_de_crf = 0.0, B_de_crf = 0.0;

              // determine lambda stepsize from min,max and nr of lambdas
              double lambda_step = (sim.param().precalclam.max_lam - 
                   sim.param().precalclam.min_lam) / 
                   (sim.param().precalclam.nr_lambdas-1);

              //loop over nr_lambdas
              for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){ 

                // determine current lambda for this index
                double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

                // start the calculations
                lj_crf_soft_interaction_ext(r, A_lj.c6, A_lj.c12,
                  B_lj.c6, B_lj.c12, A_q, B_q, alpha_lj, alpha_crf,
                  A_e_lj,  B_e_lj, A_e_crf, B_e_crf,
                  A_de_lj, B_de_lj, A_de_crf, B_de_crf,
                  lam);

                DEBUG(8, "ANITA: precalculated energies for lambda " << lam
                   << "\n now starting storage");
                DEBUG(8, "\n  A_e_lj " << A_e_lj << "\n  lambda index " << lam_index <<
                   "\n  storage.energies.A_lj_energy.size() " << storage.energies.A_lj_energy.size()
                   << "\n  energy group1 " << topo.atom_energy_group(i) << " energy group2 " 
                   << topo.atom_energy_group(j));
                //            assert(storage.energies.A_lj_energy.size() > lam_index);
                //            assert(storage.energies.A_lj_energy[lam_index].size() > topo.atom_energy_group(i));
                //            assert(storage.energies.A_lj_energy[lam_index][topo.atom_energy_group(i)].size() 
                //                     > topo.atom_energy_group(j));

                storage.energies.A_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_lj;
                storage.energies.B_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_lj;

                storage.energies.A_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_crf;
                storage.energies.B_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_crf;

                storage.perturbed_energy_derivatives.A_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_lj;
                storage.perturbed_energy_derivatives.B_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_lj;

                storage.perturbed_energy_derivatives.A_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_crf;
                storage.perturbed_energy_derivatives.B_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_crf;
                DEBUG(8, "\ndone with storing energies ");
                } //all 101 lambda points done
              } // done with extended TI
            }

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
	break;
      }
      break;
    } */
    case simulation::lj_crf_func:
    {
      DEBUG(9, "case: lj_crf_func");
      double c6 = 0.0, c12 = 0.0, q = 0.0, A_q = 0.0, B_q = 0.0, e_nb = 0.0, e_lj = 0.0,
       e_crf = 0.0, f = 0.0, de_lj = 0.0, de_crf = 0.0;
      assert(abs2(r) != 0);
      const double dist2 = abs2(r);
      const double disti = 1 / abs(r);
      const double dist6 = dist2 * dist2 * dist2;

      bool multiAEDS = false;
      if (sim.param().eds.eds == 3){
        multiAEDS = true;
      }
      /*if (sim.param().eds.numsites > 1){
        multiAEDS = true;
        std::map<std::vector<int>, math::VArray> &force_mult_endstates = storage.force_mult_endstates;
      } else {
        std::vector<math::VArray> &force_endstates = storage.force_endstates;
      } */

      // can we do this????
      std::map<std::vector<int>, math::VArray> &force_mult_endstates = storage.force_mult_endstates;
      std::vector<math::VArray> &force_endstates = storage.force_endstates;

      switch (both_perturbed) {
        case 0:
	      // EDS - normal
        {
          DEBUG(9,"\t case 0: eds-normal");
          const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
          const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();

          if (multiAEDS) {
            DEBUG(9,"\t multiAEDS");
            const int &site_i = topo.eds_perturbed_solute().atoms()[i].site_number();
            const int numstates_i = sim.param().eds.multnumstates[site_i];
            const int site_j = sim.param().eds.numsites; // non eds atom, belongs to site R (rest), last site
            const int state_j = 0; // non eds atom, just one state
            DEBUG(9,"\t\t site_i: " << site_i);
            DEBUG(9,"\t\t numstates_i: " << numstates_i);
            DEBUG(9,"\t\t site_j: " << site_j);
            DEBUG(9,"\t\t state_j: " << state_j);
            std::map<std::vector<int>, double> & storage_energies_eds_mult_vi = storage.energies.eds_mult_vi;
            std::map<std::vector<int>, math::Matrix > & storage_virial_tensor_mult_endstates = storage.virial_tensor_mult_endstates;
          
            for (unsigned int state_i = 0; state_i < numstates_i; state_i++) {
              DEBUG(9,"looping over states, now state: " << state_i);
              std::vector<int> states_i_j {site_i, state_i, site_j, state_j};
              DEBUG(9, "force_mult_endstates.size(): " << force_mult_endstates.size());
              DEBUG(9, "[states_i_j] present?: " << force_mult_endstates.count(states_i_j));
              DEBUG(9, "force_atom_i_x: " << force_mult_endstates[states_i_j](i)[0]);
              DEBUG(9, "force_atom_i_y: " << force_mult_endstates[states_i_j](i)(1));
              //DEBUG(9, "force_atom_i_z: " << force_mult_endstates_states_i[2]);
              math::Vec & force_mult_endstates_state_i = force_mult_endstates[states_i_j](i);
              math::Vec & force_mult_endstates_state_j = force_mult_endstates[states_i_j](j);
              math::Matrix & virial_tensor_mult_endstates_state = storage_virial_tensor_mult_endstates[states_i_j];//??
              const lj_parameter_struct &lj = m_param_lj_parameter[(pert_i_M_IAC[state_i])][iac_j];
              c6 = lj.c6;
              c12 = lj.c12;
              q = pert_i_M_charge[state_i] * charge_j;
              DEBUG(9, "starting with eds_lj_crf_interaction");
              eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb);

              DEBUG(9, "\t\tatomic virial");
              for (int a = 0; a < 3; ++a) {
                const double term = f * r(a);
                force_mult_endstates_state_i(a) += term;//??
                force_mult_endstates_state_j(a) -= term;//??

                for (int b = 0; b < 3; ++b)
                  virial_tensor_mult_endstates_state(b, a) += r(b) * term;
              }
              DEBUG(9, "\t\tdone with virial");
              // energy
              //assert(storage.energies.eds_vi.size() == numstates);//remove?
              storage_energies_eds_mult_vi[states_i_j] += e_nb;
              DEBUG(9, "\t\tenergies stored");
            }
          } else { //only single EDS site

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

          if (multiAEDS) {
            const int &site_i = topo.eds_perturbed_solute().atoms()[i].site_number();
            const int numstates_i = sim.param().eds.multnumstates[site_i];
            const int &site_j =topo.eds_perturbed_solute().atoms()[j].site_number();
            const int numstates_j = sim.param().eds.multnumstates[site_j];

            std::map<std::vector<int>, double> & storage_energies_eds_mult_vi = storage.energies.eds_mult_vi;
            std::map<std::vector<int>, math::Matrix> & storage_virial_tensor_mult_endstates = storage.virial_tensor_mult_endstates;
            if (site_j >= site_i) {
              for (int state_i = 0; state_i < numstates_i; state_i++) {
                for (int state_j = 0; state_j < numstates_j; state_j++) {
                  if ((site_j == site_i) && (state_i != state_j)){
                    continue; //when both atoms from same site, only include same states for both sites
                  }
                  std::vector<int> states_i_j {site_i, state_i, site_j, state_j};
                  math::Vec & force_mult_endstates_state_i = force_mult_endstates[states_i_j][i];
                  math::Vec & force_mult_endstates_state_j = force_mult_endstates[states_i_j][j];
                  math::Matrix & virial_tensor_mult_endstates_state = storage_virial_tensor_mult_endstates[states_i_j];//??
                  const lj_parameter_struct &lj =
                      m_param_lj_parameter[(pert_i_M_IAC[state_i])][(pert_j_M_IAC[state_j])];
                  c6 = lj.c6;
                  c12 = lj.c12;
                  q = pert_i_M_charge[state_i]* (pert_j_M_charge[state_j]);

                  eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb);

                  DEBUG(10, "\t\tatomic virial");
                  for (int a = 0; a < 3; ++a) {
                    const double term = f * r(a);
                    force_mult_endstates_state_i(a) += term;
                    force_mult_endstates_state_j(a) -= term;

                    for (int b = 0; b < 3; ++b)
                      virial_tensor_mult_endstates_state(b, a) += r(b) * term;
                  }
                  // energy
                  //assert(storage.energies.eds_vi.size() == numstates);
                  storage_energies_eds_mult_vi[states_i_j] += e_nb;
                }
              }
            }

          } else {
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
              // 
            }
          }
          break;
        }
        case 2:
	      // EDS - perturbed
	      {
          if (multiAEDS){
            io::messages.add("EDS-Nonbonded_Innerloop",
            "multiAEDS cannot yet be combined with a perturbation",
            io::message::critical);
          }
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


            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }

            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
					    alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            // Extended TI is not compatible with EDS yet
            if (sim.param().precalclam.nr_lambdas){
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Extended TI-EDS function not implemented",
              io::message::critical);
            }

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
          //for (unsigned int state = 0; state < numstates; state++) 
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

          if (t_perturbation_details::do_scaling) {
            // SCALING ON
            // check whether we need to do scaling
            // based on energy groups
            // Check if we want to scale interactions
            std::pair<int, int> energy_group_pair(topo.atom_energy_group(i),
              topo.atom_energy_group(j));
      
            if (topo.energy_group_scaling().count(energy_group_pair)) {

              // YES, we do scale the interactions!
              lj_crf_scaled_interaction(dist2, dist6, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12,
                A_q, B_q,
                alpha_lj, alpha_crf,
                topo.energy_group_scaling()[energy_group_pair].first,
                topo.energy_group_scaling()[energy_group_pair].second,
                f, e_lj, e_crf, de_lj, de_crf);
            }
            else{  // No scaling
                eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                  alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
            }
          } //END Scaling
          else { // No scaling 
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            //---------------------------------------------------------
            //                     ANITA
            // extended TI: calculate A_e_lj, B_e_lj, A_e_crf, B_e_crf,
            //              A_de_LJ, B_de_lj, A_de_crf, B_de_crf
            //---------------------------------------------------------

            // TODO: could add another parameter, to only calculate every x steps
            // if nr_lambdas > 1, we apply extended TI 
            if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
              DEBUG(8, "precalculate lj_crf_soft");
              //        if ( sim.param().precalclam.nr_lambdas )  
              double A_e_lj = 0.0, B_e_lj = 0.0, A_e_crf = 0.0, B_e_crf = 0.0,
                A_de_lj = 0.0, B_de_lj = 0.0, A_de_crf = 0.0, B_de_crf = 0.0;

              // determine lambda stepsize from min,max and nr of lambdas
              double lambda_step = (sim.param().precalclam.max_lam - 
                   sim.param().precalclam.min_lam) / 
                   (sim.param().precalclam.nr_lambdas-1);

              //loop over nr_lambdas
              for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){ 

              // determine current lambda for this index
              double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

              // start the calculations
              lj_crf_soft_interaction_ext(r, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12, A_q, B_q, alpha_lj, alpha_crf,
                A_e_lj,  B_e_lj, A_e_crf, B_e_crf,
                A_de_lj, B_de_lj, A_de_crf, B_de_crf,
                lam);

              DEBUG(8, "ANITA: precalculated energies for lambda " << lam
                   << "\n now starting storage");
              DEBUG(8, "\n  A_e_lj " << A_e_lj << "\n  lambda index " << lam_index <<
                   "\n  storage.energies.A_lj_energy.size() " << storage.energies.A_lj_energy.size()
                   << "\n  energy group1 " << topo.atom_energy_group(i) << " energy group2 " 
                   << topo.atom_energy_group(j));
              //            assert(storage.energies.A_lj_energy.size() > lam_index);
              //            assert(storage.energies.A_lj_energy[lam_index].size() > topo.atom_energy_group(i));
              //            assert(storage.energies.A_lj_energy[lam_index][topo.atom_energy_group(i)].size() 
              //                     > topo.atom_energy_group(j));

              storage.energies.A_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_lj;
              storage.energies.B_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_lj;

              storage.energies.A_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_crf;
              storage.energies.B_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_crf;

              storage.perturbed_energy_derivatives.A_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_lj;
              storage.perturbed_energy_derivatives.B_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_lj;

              storage.perturbed_energy_derivatives.A_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_crf;
              storage.perturbed_energy_derivatives.B_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_crf;
              DEBUG(8, "\ndone with storing energies ");
              } //all 101 lambda points done
            } // done with extended TI
          }
	  
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
          if (multiAEDS){
            io::messages.add("EDS-Nonbonded_Innerloop",
            "multiAEDS cannot yet be combined with a perturbation",
            io::message::critical);
          }

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

            // Scaling is not compatible with EDS
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }

            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
            alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            // Extended TI is not compatible with EDS yet
            if (sim.param().precalclam.nr_lambdas){
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Extended TI-EDS function not implemented",
              io::message::critical);
            }

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
		        topo.perturbed_solute().atoms()[j].LJ_softcore()) / 2.0;
          alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
            topo.perturbed_solute().atoms()[j].CRF_softcore()) / 2.0;

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
	  
          if (t_perturbation_details::do_scaling) {
            // SCALING ON
            // check whether we need to do scaling
            // based on energy groups
            // Check if we want to scale interactions
            std::pair<int, int> energy_group_pair(topo.atom_energy_group(i),
              topo.atom_energy_group(j));
      
            if (topo.energy_group_scaling().count(energy_group_pair)) {

              // YES, we do scale the interactions!
              lj_crf_scaled_interaction(dist2, dist6, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12,
                A_q, B_q,
                alpha_lj, alpha_crf,
                topo.energy_group_scaling()[energy_group_pair].first,
                topo.energy_group_scaling()[energy_group_pair].second,
                f, e_lj, e_crf, de_lj, de_crf);

            }
            else{  // No scaling
      
              eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
            }
          } //END Scaling
          else { // No scaling 
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);

            //---------------------------------------------------------
            //                     ANITA
            // extended TI: calculate A_e_lj, B_e_lj, A_e_crf, B_e_crf,
            //              A_de_LJ, B_de_lj, A_de_crf, B_de_crf
            //---------------------------------------------------------

            // TODO: could add another parameter, to only calculate every x steps
            // if nr_lambdas > 1, we apply extended TI 
            if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
              DEBUG(8, "precalculate lj_crf_soft");
              //        if ( sim.param().precalclam.nr_lambdas )  
              double A_e_lj = 0.0, B_e_lj = 0.0, A_e_crf = 0.0, B_e_crf = 0.0,
              A_de_lj = 0.0, B_de_lj = 0.0, A_de_crf = 0.0, B_de_crf = 0.0;

              // determine lambda stepsize from min,max and nr of lambdas
              double lambda_step = (sim.param().precalclam.max_lam - 
                   sim.param().precalclam.min_lam) / 
                   (sim.param().precalclam.nr_lambdas-1);

              //loop over nr_lambdas
              for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){ 

                // determine current lambda for this index
                double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;

                // start the calculations
                lj_crf_soft_interaction_ext(r, A_lj.c6, A_lj.c12,
                  B_lj.c6, B_lj.c12, A_q, B_q, alpha_lj, alpha_crf,
                  A_e_lj,  B_e_lj, A_e_crf, B_e_crf,
                  A_de_lj, B_de_lj, A_de_crf, B_de_crf,
                  lam);

                DEBUG(8, "ANITA: precalculated energies for lambda " << lam
                   << "\n now starting storage");
                DEBUG(8, "\n  A_e_lj " << A_e_lj << "\n  lambda index " << lam_index <<
                   "\n  storage.energies.A_lj_energy.size() " << storage.energies.A_lj_energy.size()
                   << "\n  energy group1 " << topo.atom_energy_group(i) << " energy group2 " 
                   << topo.atom_energy_group(j));
                //            assert(storage.energies.A_lj_energy.size() > lam_index);
                //            assert(storage.energies.A_lj_energy[lam_index].size() > topo.atom_energy_group(i));
                //            assert(storage.energies.A_lj_energy[lam_index][topo.atom_energy_group(i)].size() 
                //                     > topo.atom_energy_group(j));

                storage.energies.A_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_lj;
                storage.energies.B_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_lj;

                storage.energies.A_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_crf;
                storage.energies.B_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_crf;

                storage.perturbed_energy_derivatives.A_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_lj;
                storage.perturbed_energy_derivatives.B_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_lj;

                storage.perturbed_energy_derivatives.A_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_crf;
                storage.perturbed_energy_derivatives.B_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_crf;
                DEBUG(8, "\ndone with storing energies ");
                } //all 101 lambda points done
              } // done with extended TI
            }

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
        break;
      }
      break;
    }
    case simulation::cggromos_func:
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
          if (topo.is_coarse_grained(i)) { // CG-CG
              io::messages.add("EDS_Nonbonded_Innerloop",
              "interaction function not implemented for Gromos coarse-grained simulations!",
              io::message::critical);
          } 
          else if (topo.is_coarse_grained(j)) {
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
          }
          else{ // FG-FG interaction
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
          }
          break;
        }
        case 1:
        // EDS - EDS
        {
          if (topo.is_coarse_grained(i) || topo.is_coarse_grained(j)) { // CG-CG or FG-CG
            io::messages.add("EDS_Nonbonded_Innerloop",
                "interaction function not implemented for Gromos coarse-grained simulations!",
                io::message::critical);
          } 
          else { // FG -FG
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
          }
          break;
        }
        case 2:
        // EDS - perturbed
        {
          io::messages.add("EDS_Nonbonded_Innerloop",
          "EDS - TI not implemented for Gromos coarse-grained simulations!",
          io::message::critical);
          break;
        }
        case 3:
        //perturbed - normal
        {

          alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
          alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();
 
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

          // Scaling is not compatible
          if (t_perturbation_details::do_scaling) {
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Scaling function not implemented",
            io::message::critical);
          }

          if (topo.is_coarse_grained(i) && topo.is_coarse_grained(j)) { // CG-CG interaction
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, 
              A_q / cgrain_eps[0], B_q / cgrain_eps[0],
              alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
          } else if (topo.is_coarse_grained(i) || topo.is_coarse_grained(j)) { // FG-CG interaction
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, 
              A_q / cgrain_eps[1], B_q / cgrain_eps[1],
              alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
          } else { // FG-FG interaction
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, 
              A_q, B_q, alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
          }


          // Extended TI is not compatible
          if (sim.param().precalclam.nr_lambdas){
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Extended TI function not implemented",
            io::message::critical);
          }

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
          io::messages.add("EDS_Nonbonded_Innerloop",
          "EDS - TI not implemented for Gromos coarse-grained simulations!",
          io::message::critical);
          break;
        }
        case 5:
        // perturbed - perturbed
        {
          alpha_lj = (topo.perturbed_solute().atoms()[i].LJ_softcore() +
		        topo.perturbed_solute().atoms()[j].LJ_softcore()) / 2.0;
          alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
            topo.perturbed_solute().atoms()[j].CRF_softcore()) / 2.0;

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

          DEBUG (10, "lambdas set");

          // Scaling is not compatible
          if (t_perturbation_details::do_scaling) {
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Scaling function not implemented",
            io::message::critical);
          }

          if (topo.is_coarse_grained(i) && topo.is_coarse_grained(j)) { // CG-CG interaction
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, 
              A_q / cgrain_eps[0], B_q / cgrain_eps[0],
              alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
          } else if (topo.is_coarse_grained(i) || topo.is_coarse_grained(j)) { // FG-CG interaction
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, 
              A_q / cgrain_eps[1], B_q / cgrain_eps[1],
              alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
          } else { // FG-FG interaction
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, 
              A_q, B_q, alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf);
          }

          // Extended TI is not compatible
          if (sim.param().precalclam.nr_lambdas){
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Extended TI function not implemented",
            io::message::critical);
          }

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
        break;
      }
      break;
    }
    case simulation::pol_lj_crf_func:
    {
      DEBUG(7, "\tpol_lj_crf_func");
      math::Vec rp1, rp2, rpp;
      double f1[4];
      math::VArray f(4);
      f = 0.0;
      double f6 = 0.0, f12 = 0.0;
      double A_qi = 0.0, A_qj = 0.0, B_qi = 0.0, B_qj = 0.0, A_q = 0.0, B_q = 0.0;
      double e_lj = 0.0, e_crf = 0.0, de_lj = 0.0, de_crf = 0.0;

      rp1 = r - conf.current().posV(j);
      rp2 = r + conf.current().posV(i);
      rpp = r + conf.current().posV(i) - conf.current().posV(j);
      switch (both_perturbed) {
        case 0:
        // EDS - normal
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
            break;
        }
        case 1:
        // EDS - EDS
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
            break;
        }
        case 2:
        // EDS - perturbed
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
            break;
        }
        case 3:
        // perturbed - normal
        {
          DEBUG(10, "perturbed - normal");
          alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
          alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

          const lj_parameter_struct &A_lj = 
            m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.iac(j)];
          const lj_parameter_struct &B_lj = 
            m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.iac(j)];

          A_qi = topo.perturbed_solute().atoms()[i].A_charge();
          A_qj = topo.charge()(j);
          B_qi = topo.perturbed_solute().atoms()[i].B_charge();
          B_qj = topo.charge()(j);
          A_q = A_qi * A_qj;
          B_q = B_qi * B_qj;               

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

          // Scaling is not compatible
          if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
          }
          else { // No scaling 
            // In this case, we can store everything immediately and do not need the endstates
            pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
                    A_lj.cs6, A_lj.cs12,
                    B_lj.cs6, B_lj.cs12,
                    A_qi, B_qi, A_qj, B_qj,
                    topo.coscharge(i),
                    topo.coscharge(j),
                    alpha_lj, alpha_crf,
                    f1, f6, f12,
                    e_lj, e_crf, de_lj, de_crf);
          }

          // Extended TI is not compatible
          if (sim.param().precalclam.nr_lambdas){
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Extended TI function not implemented",
              io::message::critical);
          }

          //--------------------------------------------------
          // interactions have been calculated
          //--------------------------------------------------

          // now combine everything
          f(0) = (f1[0] + f6 + f12) * r;
          f(1) = f1[1] * rp1;
          f(2) = f1[2] * rp2;
          f(3) = f1[3] * rpp;

          conf.current().force(i) += f(0) + f(1) + f(2) + f(3);
          conf.current().force(j) -= f(0) + f(1) + f(2) + f(3);

          DEBUG(7, "\tforces stored");

          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
              conf.current().virial_tensor(a, b) += r(a)*(f(0)(b) +
                      f(1)(b) + f(2)(b) + f(3)(b));

          DEBUG(7, "\tatomic virial done");

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
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
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
          
          A_qi = topo.perturbed_solute().atoms()[i].A_charge();
          A_qj = topo.perturbed_solute().atoms()[j].A_charge();
          B_qi = topo.perturbed_solute().atoms()[i].B_charge();
          B_qj = topo.perturbed_solute().atoms()[j].B_charge();
          A_q = A_qi * A_qj;
          B_q = B_qi * B_qj;   

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
          

          // Scaling is not compatible
          if (t_perturbation_details::do_scaling) {
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Scaling function not implemented",
            io::message::critical);
          }
          else { // No scaling 
            // In this case, we can store everything immediately and do not need the endstates
            pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
                  A_lj.cs6, A_lj.cs12,
                  B_lj.cs6, B_lj.cs12,
                  A_qi, B_qi, A_qj, B_qj,
                  topo.coscharge(i),
                  topo.coscharge(j),
                  alpha_lj, alpha_crf,
                  f1, f6, f12,
                  e_lj, e_crf, de_lj, de_crf);
          }

          // Extended TI is not compatible
          if (sim.param().precalclam.nr_lambdas){
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Extended TI function not implemented",
            io::message::critical);
          }

          //--------------------------------------------------
          // interactions have been calculated
          //--------------------------------------------------

          // now combine everything
          f(0) = (f1[0] + f6 + f12) * r;
          f(1) = f1[1] * rp1;
          f(2) = f1[2] * rp2;
          f(3) = f1[3] * rpp;

          conf.current().force(i) += f(0) + f(1) + f(2) + f(3);
          conf.current().force(j) -= f(0) + f(1) + f(2) + f(3);

          DEBUG(7, "\tforces stored");

          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
              conf.current().virial_tensor(a, b) += r(a)*(f(0)(b) +
                    f(1)(b) + f(2)(b) + f(3)(b));

          DEBUG(7, "\tatomic virial done");
                    
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
        break;   
      }
      break;
    }
    case simulation::pol_off_lj_crf_func:
    {
      DEBUG(7, "\tpol_off_lj_crf_func");
      math::Vec rm = r;
      if (topo.gamma(i)!=0.0) {
        math::Vec rij, rik;
        periodicity.nearest_image(conf.current().pos(i),
                conf.current().pos(topo.gamma_j(i)), rij);
        periodicity.nearest_image(conf.current().pos(i),
                conf.current().pos(topo.gamma_k(i)), rik);
        rm -= topo.gamma(i)*(rij + rik) / 2;
      }
      if (topo.gamma(j)!=0.0) {
        math::Vec rjj, rjk;
        periodicity.nearest_image(conf.current().pos(j),
                conf.current().pos(topo.gamma_k(j)), rjk);
        periodicity.nearest_image(conf.current().pos(j),
                conf.current().pos(topo.gamma_j(j)), rjj);

        rm += topo.gamma(j)*(rjj + rjk) / 2;
      }
      math::Vec rp1, rp2, rpp;
      double f1[4];
      math::VArray f(4);
      f = 0.0;
      double f6 = 0.0, f12 = 0.0;
      double A_qi = 0.0, A_qj = 0.0, B_qi = 0.0, B_qj = 0.0, A_q = 0.0, B_q = 0.0;
      double e_lj = 0.0, e_crf = 0.0, de_lj = 0.0, de_crf = 0.0;

      rp1 = r - conf.current().posV(j);
      rp2 = r + conf.current().posV(i);
      rpp = r + conf.current().posV(i) - conf.current().posV(j);
      switch (both_perturbed) {
        case 0:
        // EDS - normal
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
            break;
        }
        case 1:
        // EDS - EDS
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
            break;
        }
        case 2:
        // EDS - perturbed
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
            break;
        }
        case 3:
        // perturbed - normal
        {
          DEBUG(10, "perturbed - normal");
          alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
          alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

          const lj_parameter_struct &A_lj = 
            m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.iac(j)];
          const lj_parameter_struct &B_lj = 
            m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.iac(j)];

          A_qi = topo.perturbed_solute().atoms()[i].A_charge();
          A_qj = topo.charge()(j);
          B_qi = topo.perturbed_solute().atoms()[i].B_charge();
          B_qj = topo.charge()(j);
          A_q = A_qi * A_qj;
          B_q = B_qi * B_qj;         

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

          // Scaling is not compatible
          if (t_perturbation_details::do_scaling) {
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Scaling function not implemented",
            io::message::critical);
          }

          pol_off_lj_crf_soft_interaction(r, rm, rp1, rp2, rpp,
                A_lj.cs6, A_lj.cs12,
                B_lj.cs6, B_lj.cs12,
                A_qi, B_qi, A_qj, B_qj,
                topo.coscharge(i),
                topo.coscharge(j),
                alpha_lj, alpha_crf,
                f1, f6, f12,
                e_lj, e_crf, de_lj, de_crf); 

          // Extended TI is not compatible
          if (sim.param().precalclam.nr_lambdas){
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Extended TI function not implemented",
            io::message::critical);
          }    
          //--------------------------------------------------
          // interactions have been calculated
          //--------------------------------------------------

          // now combine everything
          f(0) = f1[0] * rm;
          f(1) = f1[1] * rp1;
          f(2) = f1[2] * rp2;
          f(3) = f1[3] * rpp;
          for (int a = 0; a < 3; ++a) {
            const double term = f(0)(a) + f(1)(a)
                    + f(2)(a) + f(3)(a);

            storage.force(i)(a) += (1 - topo.gamma(i)) * term + (f6 + f12) * r(a);
            storage.force(j)(a) -= (1 - topo.gamma(j)) * term + (f6 + f12) * r(a);

            storage.force(topo.gamma_j(i))(a) += topo.gamma(i) / 2 * term;
            storage.force(topo.gamma_j(j))(a) -= topo.gamma(j) / 2 * term;

            storage.force(topo.gamma_k(i))(a) += topo.gamma(i) / 2 * term;
            storage.force(topo.gamma_k(j))(a) -= topo.gamma(j) / 2 * term;

            DEBUG(7, "\tforces stored");

            for (int b = 0; b < 3; ++b)
              storage.virial_tensor(b, a) += r(b) * term + r(b)*(f6 + f12) * r(a);
          }

          DEBUG(7, "\tatomic virial done");

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
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
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
          
          A_qi = topo.perturbed_solute().atoms()[i].A_charge();
          A_qj = topo.perturbed_solute().atoms()[j].A_charge();
          B_qi = topo.perturbed_solute().atoms()[i].B_charge();
          B_qj = topo.perturbed_solute().atoms()[j].B_charge();
          A_q = A_qi * A_qj;
          B_q = B_qi * B_qj;   

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
          

          // Scaling is not compatible
          if (t_perturbation_details::do_scaling) {
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Scaling function not implemented",
            io::message::critical);
          }
          else { // No scaling 
            pol_off_lj_crf_soft_interaction(r, rm, rp1, rp2, rpp,
                  A_lj.cs6, A_lj.cs12,
                  B_lj.cs6, B_lj.cs12,
                  A_qi, B_qi, A_qj, B_qj,
                  topo.coscharge(i),
                  topo.coscharge(j),
                  alpha_lj, alpha_crf,
                  f1, f6, f12,
                  e_lj, e_crf, de_lj, de_crf); 
          }

          // Extended TI is not compatible
          if (sim.param().precalclam.nr_lambdas){
            io::messages.add("EDS-Nonbonded_Innerloop",
            "Extended TI function not implemented",
            io::message::critical);
          }

          //--------------------------------------------------
          // interactions have been calculated
          //--------------------------------------------------

          // now combine everything
          f(0) = f1[0] * rm;
          f(1) = f1[1] * rp1;
          f(2) = f1[2] * rp2;
          f(3) = f1[3] * rpp;
          for (int a = 0; a < 3; ++a) {
            const double term = f(0)(a) + f(1)(a)
                    + f(2)(a) + f(3)(a);

            storage.force(i)(a) += (1 - topo.gamma(i)) * term + (f6 + f12) * r(a);
            storage.force(j)(a) -= (1 - topo.gamma(j)) * term + (f6 + f12) * r(a);

            storage.force(topo.gamma_j(i))(a) += topo.gamma(i) / 2 * term;
            storage.force(topo.gamma_j(j))(a) -= topo.gamma(j) / 2 * term;

            storage.force(topo.gamma_k(i))(a) += topo.gamma(i) / 2 * term;
            storage.force(topo.gamma_k(j))(a) -= topo.gamma(j) / 2 * term;

            DEBUG(7, "\tforces stored");

            for (int b = 0; b < 3; ++b)
              storage.virial_tensor(b, a) += r(b) * term + r(b)*(f6 + f12) * r(a);
          }


          DEBUG(7, "\tatomic virial done");

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
        break;  
      }
      break;
    }
    case simulation::cgrain_func:
    {
      double c6 = 0.0, c12 = 0.0, q = 0.0, A_q = 0.0, B_q = 0.0, e_nb = 0.0, e_lj = 0.0,
      e_crf = 0.0, f = 0.0, de_lj = 0.0, de_crf = 0.0;
      const std::vector<std::vector<lj_parameter_struct> > &m_param_cg_parameter = m_param->cg_parameter();

      switch (both_perturbed) {
        case 0:
        // EDS - normal
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
            break;
        }
        case 1:
        // EDS - EDS
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
            break;
        }
        case 2:
        // EDS - perturbed
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
            break;
        }
        case 3:
        // perturbed - normal
        {
          DEBUG(10, "perturbed - normal");
          alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
          alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

          const lj_parameter_struct &A_lj = 
            m_param_cg_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.iac(j)];
          const lj_parameter_struct &B_lj = 
            m_param_cg_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.iac(j)];
          
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

            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }

            cgrain_soft_interaction(r, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12, A_q, B_q,alpha_lj, alpha_crf,
                f, e_lj, e_crf, de_lj, de_crf);

            // Extended TI is not compatible
            if (sim.param().precalclam.nr_lambdas){
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Extended TI function not implemented",
              io::message::critical);
            }

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
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
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
              m_param_cg_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.perturbed_solute().atoms()[j].A_IAC()];
            const lj_parameter_struct &B_lj = 
              m_param_cg_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.perturbed_solute().atoms()[j].B_IAC()];
            
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

            DEBUG (10, "lambdas set");

            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }

            cgrain_soft_interaction(r, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12, A_q, B_q,alpha_lj, alpha_crf,
                f, e_lj, e_crf, de_lj, de_crf);

            // Extended TI is not compatible
            if (sim.param().precalclam.nr_lambdas){
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Extended TI function not implemented",
              io::message::critical);
            }

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
        break;   
      }
      break;
    }
    default:
      io::messages.add("EDS-Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
    break;
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
        Periodicity_type const & periodicity,
        simulation::Simulation & sim) {
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

  bool multiAEDS = false;
  if (sim.param().eds.eds == 3){
    multiAEDS = true;
  }

  switch (t_interaction_spec::interaction_func) {
    case simulation::cggromos_func:
    {
      if (topo.is_coarse_grained(i)) { // CG-FG or CG-CG
        io::messages.add("EDS_Nonbonded_Innerloop",
        "1,4-interaction function not implemented for Gromos coarse-grained simulations!",
        io::message::critical);
      }
      // No break statement makes next block execute and thus the interaction is treated as FG-FG 
    }
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
        // EDS - normal
        {
	        const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
	        const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();

          if (multiAEDS){
            std::map<std::vector<int>,math::Matrix> & conf_special_eds_virial_tensor_mult_endstates = conf.special().eds.virial_tensor_mult_endstates;
            std::map<std::vector<int>,double> & conf_current_energies_eds_mult_vi = conf.current().energies.eds_mult_vi;
            std::map<std::vector<int>, math::VArray> & conf_special_eds_force_mult_endstates = conf.special().eds.force_mult_endstates;

            const int &site_i = topo.eds_perturbed_solute().atoms()[i].site_number();
            const int numstates_i = sim.param().eds.multnumstates[site_i];
            const int site_j = sim.param().eds.numsites; // non eds atom, belongs to site R (rest), last site
            const int state_j = 0; // non eds atom, just one state

            for (unsigned int state_i = 0; state_i < numstates_i; state_i++) {
              std::vector<int> states_i_j {site_i, state_i, site_j, state_j};
              math::Vec & conf_special_eds_force_mult_endstates_state_i = conf_special_eds_force_mult_endstates[states_i_j][i];
              math::Vec & conf_special_eds_force_mult_endstates_state_j = conf_special_eds_force_mult_endstates[states_i_j][j];
              math::Matrix & conf_special_eds_virial_tensor_mult_endstates_state = conf_special_eds_virial_tensor_mult_endstates[states_i_j];//??
              const lj_parameter_struct &lj = m_param_lj_parameter[(pert_i_M_IAC[state_i])][iac_j];
              c6 = lj.cs6;
              c12 = lj.cs12;
              q = pert_i_M_charge[state_i] * charge_j;

              eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb, 0, m_param->get_coulomb_scaling());

              DEBUG(10, "\t\tatomic virial");
              for (int a = 0; a < 3; ++a) {
                const double term = f * r(a);
                conf_special_eds_force_mult_endstates_state_i(a) += term;
                conf_special_eds_force_mult_endstates_state_j(a) -= term;

                for (int b = 0; b < 3; ++b)
                  conf_special_eds_virial_tensor_mult_endstates_state(b, a) += r(b) * term;
              }
              // energy
              //assert(conf_current_energies_eds_vi.size() == numstates);
              conf_current_energies_eds_mult_vi[states_i_j] += e_nb;
            }

          } else {
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

          if (multiAEDS){
            const int &site_i = topo.eds_perturbed_solute().atoms()[i].site_number();
            const int numstates_i = sim.param().eds.multnumstates[site_i];
            const int &site_j =topo.eds_perturbed_solute().atoms()[j].site_number();
            const int numstates_j = sim.param().eds.multnumstates[site_j];

            std::map<std::vector<int>, double> & conf_current_energies_eds_mult_vi = conf.current().energies.eds_mult_vi;
            std::map<std::vector<int>, math::Matrix> & conf_special_eds_virial_tensor_mult_endstates = conf.special().eds.virial_tensor_mult_endstates;
            std::map<std::vector<int>, math::VArray> & conf_special_eds_force_mult_endstates = conf.special().eds.force_mult_endstates;
            
            if (site_j >= site_i) {
              for (int state_i = 0; state_i < numstates_i; state_i++) {
                for (int state_j = 0; state_j < numstates_j; state_j++) {
                  if ((site_j == site_i) && (state_i != state_j)){
                    continue; //when both atoms from same site, only include same states for both sites
                  }
                  std::vector<int> states_i_j {site_i, state_i, site_j, state_j};
                  math::Vec & conf_special_eds_force_mult_endstates_state_i = conf_special_eds_force_mult_endstates[states_i_j][i];
                  math::Vec & conf_special_eds_force_mult_endstates_state_j = conf_special_eds_force_mult_endstates[states_i_j][j];
                  math::Matrix & conf_special_eds_virial_tensor_mult_endstates_state = conf_special_eds_virial_tensor_mult_endstates[states_i_j];//??
                  const lj_parameter_struct &lj =
                      m_param_lj_parameter[(pert_i_M_IAC[state_i])][(pert_j_M_IAC[state_j])];
                  c6 = lj.cs6;
                  c12 = lj.cs12;
                  q = pert_i_M_charge[state_i]* (pert_j_M_charge[state_j]);

                  eds_lj_crf_interaction(dist2, dist6, disti, c6, c12, q, f, e_nb, 0, m_param->get_coulomb_scaling());

                  DEBUG(10, "\t\tatomic virial");
                  for (int a = 0; a < 3; ++a) {
                    const double term = f * r(a);
                    conf_special_eds_force_mult_endstates_state_i(a) += term;
                    conf_special_eds_force_mult_endstates_state_j(a) -= term;

                    for (int b = 0; b < 3; ++b)
                      conf_special_eds_virial_tensor_mult_endstates_state(b, a) += r(b) * term;
                  }
                  // energy
                  //assert(storage.energies.eds_vi.size() == numstates);
                  conf_current_energies_eds_mult_vi[states_i_j] += e_nb;
                }
              }
            }

          } else {
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

            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }

            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.cs6, A_lj.cs12, B_lj.cs6, B_lj.cs12, A_q, B_q,
					    alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf, 0, m_param->get_coulomb_scaling());

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

          if (t_perturbation_details::do_scaling) {
            // SCALING ON
            // check whether we need to do scaling
            // based on energy groups
            // Check if we want to scale interactions
            std::pair<int, int> energy_group_pair(topo.atom_energy_group(i),
              topo.atom_energy_group(j));
      
            if (topo.energy_group_scaling().count(energy_group_pair)) {

              // YES, we do scale the interactions!
              lj_crf_scaled_interaction(dist2, dist6, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12,
                A_q, B_q,
                alpha_lj, alpha_crf,
                topo.energy_group_scaling()[energy_group_pair].first,
                topo.energy_group_scaling()[energy_group_pair].second,
                f, e_lj, e_crf, de_lj, de_crf, 0, m_param->get_coulomb_scaling());
            }
            else{  // No scaling
                eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                  alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf, 0, m_param->get_coulomb_scaling());
            }
          } //END Scaling
          else {	  
	          eds_pert_lj_crf_interaction(dist2, dist6, A_lj.cs6, A_lj.cs12, B_lj.cs6, B_lj.cs12, A_q, B_q,
				      alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf, 0, m_param->get_coulomb_scaling());
          }
	  
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

            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }

            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.cs6, A_lj.cs12, B_lj.cs6, B_lj.cs12, A_q, B_q,
					  alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf, 0, m_param->get_coulomb_scaling());

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
		        topo.perturbed_solute().atoms()[j].LJ_softcore()) / 2.0;
          alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
		       topo.perturbed_solute().atoms()[j].CRF_softcore()) / 2.0;

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
	  
          if (t_perturbation_details::do_scaling) {
            // SCALING ON
            // check whether we need to do scaling
            // based on energy groups
            // Check if we want to scale interactions
            std::pair<int, int> energy_group_pair(topo.atom_energy_group(i),
              topo.atom_energy_group(j));
      
            if (topo.energy_group_scaling().count(energy_group_pair)) {

              // YES, we do scale the interactions!
              lj_crf_scaled_interaction(dist2, dist6, A_lj.c6, A_lj.c12,
                B_lj.c6, B_lj.c12,
                A_q, B_q,
                alpha_lj, alpha_crf,
                topo.energy_group_scaling()[energy_group_pair].first,
                topo.energy_group_scaling()[energy_group_pair].second,
                f, e_lj, e_crf, de_lj, de_crf, 0, m_param->get_coulomb_scaling());

            }
            else{  // No scaling
      
              eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf, 0, m_param->get_coulomb_scaling());
            }
          } //END Scaling
          else { // No scaling 
            eds_pert_lj_crf_interaction(dist2, dist6, A_lj.c6, A_lj.c12, B_lj.c6, B_lj.c12, A_q, B_q,
                alpha_lj, alpha_crf, f, e_lj, e_crf, de_lj, de_crf, 0, m_param->get_coulomb_scaling());
          } 
          

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
        break;
      }
      break;
    }
    case simulation::pol_lj_crf_func:
    {
      DEBUG(7, "\tpol_lj_crf_func");
      math::Vec rp1, rp2, rpp;
      double f1[4];
      math::VArray f(4);
      f = 0.0;
      double f6 = 0.0, f12 = 0.0;
      double A_qi = 0.0, A_qj = 0.0, B_qi = 0.0, B_qj = 0.0, A_q = 0.0, B_q = 0.0;
      double e_lj = 0.0, e_crf = 0.0, de_lj = 0.0, de_crf = 0.0;

      rp1 = r - conf.current().posV(j);
      rp2 = r + conf.current().posV(i);
      rpp = r + conf.current().posV(i) - conf.current().posV(j);
      switch (both_perturbed)
      {
        case 0:
        // EDS - Normal
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
          break;
        }
        case 1:
        // EDS - EDS
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
          break;
        }
        case 2:
        // EDS - perturbed
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
          break;
        }
        case 3:
        // perturbed - normal
        {
            DEBUG(10, "perturbed - normal");
            alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
            alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

            const lj_parameter_struct &A_lj = 
              m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.iac(j)];
            const lj_parameter_struct &B_lj = 
              m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.iac(j)];

            A_qi = topo.perturbed_solute().atoms()[i].A_charge();
            A_qj = topo.charge()(j);
            B_qi = topo.perturbed_solute().atoms()[i].B_charge();
            B_qj = topo.charge()(j);
            A_q = A_qi * A_qj;
            B_q = B_qi * B_qj;               

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

            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
                io::messages.add("EDS-Nonbonded_Innerloop",
                "Scaling function not implemented",
                io::message::critical);
            }
            else { // No scaling 
              // In this case, we can store everything immediately and do not need the endstates
              pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
                      A_lj.cs6, A_lj.cs12,
                      B_lj.cs6, B_lj.cs12,
                      A_qi, B_qi, A_qj, B_qj,
                      topo.coscharge(i),
                      topo.coscharge(j),
                      alpha_lj, alpha_crf,
                      f1, f6, f12,
                      e_lj, e_crf, de_lj, de_crf);
            }

            // Extended TI is not compatible
            if (sim.param().precalclam.nr_lambdas){
                io::messages.add("EDS-Nonbonded_Innerloop",
                "Extended TI function not implemented",
                io::message::critical);
            }

            //--------------------------------------------------
            // interactions have been calculated
            //--------------------------------------------------

            // now combine everything
            f(0) = (f1[0] + f6 + f12) * r;
            f(1) = f1[1] * rp1;
            f(2) = f1[2] * rp2;
            f(3) = f1[3] * rpp;

            conf.current().force(i) += f(0) + f(1) + f(2) + f(3);
            conf.current().force(j) -= f(0) + f(1) + f(2) + f(3);

            DEBUG(7, "\tforces stored");

            for (int a = 0; a < 3; ++a)
              for (int b = 0; b < 3; ++b)
                conf.current().virial_tensor(a, b) += r(a)*(f(0)(b) +
                        f(1)(b) + f(2)(b) + f(3)(b));

            DEBUG(7, "\tatomic virial done");

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
              
            DEBUG(10, "forces, energies and energy derivatives stored in storage");
          break;
        }
        case 4:
        // perturbed - EDS
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
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
            
            A_qi = topo.perturbed_solute().atoms()[i].A_charge();
            A_qj = topo.perturbed_solute().atoms()[j].A_charge();
            B_qi = topo.perturbed_solute().atoms()[i].B_charge();
            B_qj = topo.perturbed_solute().atoms()[j].B_charge();
            A_q = A_qi * A_qj;
            B_q = B_qi * B_qj;   

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
            

            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }
            else { // No scaling 
              // In this case, we can store everything immediately and do not need the endstates
              pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
                    A_lj.cs6, A_lj.cs12,
                    B_lj.cs6, B_lj.cs12,
                    A_qi, B_qi, A_qj, B_qj,
                    topo.coscharge(i),
                    topo.coscharge(j),
                    alpha_lj, alpha_crf,
                    f1, f6, f12,
                    e_lj, e_crf, de_lj, de_crf);
            }

            // Extended TI is not compatible
            if (sim.param().precalclam.nr_lambdas){
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Extended TI function not implemented",
              io::message::critical);
            }

            //--------------------------------------------------
            // interactions have been calculated
            //--------------------------------------------------

            // now combine everything
            f(0) = (f1[0] + f6 + f12) * r;
            f(1) = f1[1] * rp1;
            f(2) = f1[2] * rp2;
            f(3) = f1[3] * rpp;

            conf.current().force(i) += f(0) + f(1) + f(2) + f(3);
            conf.current().force(j) -= f(0) + f(1) + f(2) + f(3);

            DEBUG(7, "\tforces stored");

            for (int a = 0; a < 3; ++a)
              for (int b = 0; b < 3; ++b)
                conf.current().virial_tensor(a, b) += r(a)*(f(0)(b) +
                      f(1)(b) + f(2)(b) + f(3)(b));

            DEBUG(7, "\tatomic virial done");
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
        break;            
      }
      break;
    }
    case simulation::pol_off_lj_crf_func:
    {
      DEBUG(7, "\tpol_off_lj_crf_func");
      math::Vec rm = r;
      if (topo.gamma(i)!=0.0) {
        math::Vec rij, rik;
        periodicity.nearest_image(conf.current().pos(i),
                conf.current().pos(topo.gamma_j(i)), rij);
        periodicity.nearest_image(conf.current().pos(i),
                conf.current().pos(topo.gamma_k(i)), rik);
        rm -= topo.gamma(i)*(rij + rik) / 2;
      }
      if (topo.gamma(j)!=0.0) {
        math::Vec rjj, rjk;
        periodicity.nearest_image(conf.current().pos(j),
                conf.current().pos(topo.gamma_k(j)), rjk);
        periodicity.nearest_image(conf.current().pos(j),
                conf.current().pos(topo.gamma_j(j)), rjj);

        rm += topo.gamma(j)*(rjj + rjk) / 2;
      }
      math::Vec rp1, rp2, rpp;
      double f1[4];
      math::VArray f(4);
      f = 0.0;
      double f6 = 0.0, f12 = 0.0;
      double A_qi = 0.0, A_qj = 0.0, B_qi = 0.0, B_qj = 0.0, A_q = 0.0, B_q = 0.0;
      double e_lj = 0.0, e_crf = 0.0, de_lj = 0.0, de_crf = 0.0;

      rp1 = r - conf.current().posV(j);
      rp2 = r + conf.current().posV(i);
      rpp = r + conf.current().posV(i) - conf.current().posV(j);
      switch (both_perturbed)
      {
        case 0:
        // EDS - Normal
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
          break;
        }
        case 1:
        // EDS - EDS
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
          break;
        }
        case 2:
        // EDS - perturbed
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
          break;
        }
        case 3:
        // perturbed - normal
        {
            DEBUG(10, "perturbed - normal");
            alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
            alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

            const lj_parameter_struct &A_lj = 
              m_param_lj_parameter[topo.perturbed_solute().atoms()[i].A_IAC()][topo.iac(j)];
            const lj_parameter_struct &B_lj = 
              m_param_lj_parameter[topo.perturbed_solute().atoms()[i].B_IAC()][topo.iac(j)];

            A_qi = topo.perturbed_solute().atoms()[i].A_charge();
            A_qj = topo.charge()(j);
            B_qi = topo.perturbed_solute().atoms()[i].B_charge();
            B_qj = topo.charge()(j);
            A_q = A_qi * A_qj;
            B_q = B_qi * B_qj;               

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

            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
                io::messages.add("EDS-Nonbonded_Innerloop",
                "Scaling function not implemented",
                io::message::critical);
            }
            else { // No scaling 
              // In this case, we can store everything immediately and do not need the endstates
              pol_off_lj_crf_soft_interaction(r, rm, rp1, rp2, rpp,
                      A_lj.cs6, A_lj.cs12,
                      B_lj.cs6, B_lj.cs12,
                      A_qi, B_qi, A_qj, B_qj,
                      topo.coscharge(i),
                      topo.coscharge(j),
                      alpha_lj, alpha_crf,
                      f1, f6, f12,
                      e_lj, e_crf, de_lj, de_crf);
            }

            // Extended TI is not compatible
            if (sim.param().precalclam.nr_lambdas){
                io::messages.add("EDS-Nonbonded_Innerloop",
                "Extended TI function not implemented",
                io::message::critical);
            }

            //--------------------------------------------------
            // interactions have been calculated
            //--------------------------------------------------

            // now combine everything
            f(0) = f1[0] * rm;
            f(1) = f1[1] * rp1;
            f(2) = f1[2] * rp2;
            f(3) = f1[3] * rpp;
            for (int a = 0; a < 3; ++a) {
              const double term = f(0)(a) + f(1)(a)
                      + f(2)(a) + f(3)(a);

              conf.current().force(i)(a) += (1 - topo.gamma(i)) * term + (f6 + f12) * r(a);
              conf.current().force(j)(a) -= (1 - topo.gamma(j)) * term + (f6 + f12) * r(a);

              conf.current().force(topo.gamma_j(i))(a) += topo.gamma(i) / 2 * term;
              conf.current().force(topo.gamma_j(j))(a) -= topo.gamma(j) / 2 * term;

              conf.current().force(topo.gamma_k(i))(a) += topo.gamma(i) / 2 * term;
              conf.current().force(topo.gamma_k(j))(a) -= topo.gamma(j) / 2 * term;

              for (int b = 0; b < 3; ++b)
                conf.current().virial_tensor(b, a) += rm(b) * term + (f6 + f12) * r(a) * r(b);
            }
            DEBUG(7, "\tforces stored");
            DEBUG(7, "\tatomic virial done");

            conf.current().force(i) += f(0) + f(1) + f(2) + f(3);
            conf.current().force(j) -= f(0) + f(1) + f(2) + f(3);

            DEBUG(7, "\tforces stored");

            for (int a = 0; a < 3; ++a)
              for (int b = 0; b < 3; ++b)
                conf.current().virial_tensor(a, b) += r(a)*(f(0)(b) +
                      f(1)(b) + f(2)(b) + f(3)(b));

            DEBUG(7, "\tatomic virial done");

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
              
            DEBUG(10, "forces, energies and energy derivatives stored in storage");
          break;
        }
        case 4:
        // perturbed - EDS
        {
          io::messages.add("EDS-Nonbonded_Innerloop",
            "interaction function not implemented",
            io::message::critical);
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
            
            A_qi = topo.perturbed_solute().atoms()[i].A_charge();
            A_qj = topo.perturbed_solute().atoms()[j].A_charge();
            B_qi = topo.perturbed_solute().atoms()[i].B_charge();
            B_qj = topo.perturbed_solute().atoms()[j].B_charge();
            A_q = A_qi * A_qj;
            B_q = B_qi * B_qj;   

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
            

            // Scaling is not compatible
            if (t_perturbation_details::do_scaling) {
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Scaling function not implemented",
              io::message::critical);
            }
            else { // No scaling 
              // In this case, we can store everything immediately and do not need the endstates
              pol_off_lj_crf_soft_interaction(r, rm,rp1, rp2, rpp,
                    A_lj.cs6, A_lj.cs12,
                    B_lj.cs6, B_lj.cs12,
                    A_qi, B_qi, A_qj, B_qj,
                    topo.coscharge(i),
                    topo.coscharge(j),
                    alpha_lj, alpha_crf,
                    f1, f6, f12,
                    e_lj, e_crf, de_lj, de_crf);
            }

            // Extended TI is not compatible
            if (sim.param().precalclam.nr_lambdas){
              io::messages.add("EDS-Nonbonded_Innerloop",
              "Extended TI function not implemented",
              io::message::critical);
            }

            //--------------------------------------------------
            // interactions have been calculated
            //--------------------------------------------------

            // now combine everything
            f(0) = f1[0] * rm;
            f(1) = f1[1] * rp1;
            f(2) = f1[2] * rp2;
            f(3) = f1[3] * rpp;
            for (int a = 0; a < 3; ++a) {
              const double term = f(0)(a) + f(1)(a)
                      + f(2)(a) + f(3)(a);

              conf.current().force(i)(a) += (1 - topo.gamma(i)) * term + (f6 + f12) * r(a);
              conf.current().force(j)(a) -= (1 - topo.gamma(j)) * term + (f6 + f12) * r(a);

              conf.current().force(topo.gamma_j(i))(a) += topo.gamma(i) / 2 * term;
              conf.current().force(topo.gamma_j(j))(a) -= topo.gamma(j) / 2 * term;

              conf.current().force(topo.gamma_k(i))(a) += topo.gamma(i) / 2 * term;
              conf.current().force(topo.gamma_k(j))(a) -= topo.gamma(j) / 2 * term;

              for (int b = 0; b < 3; ++b)
                conf.current().virial_tensor(b, a) += rm(b) * term + (f6 + f12) * r(a) * r(b);
            }
            DEBUG(7, "\tforces stored");
            DEBUG(7, "\tatomic virial done");

            conf.current().force(i) += f(0) + f(1) + f(2) + f(3);
            conf.current().force(j) -= f(0) + f(1) + f(2) + f(3);

            DEBUG(7, "\tforces stored");

            for (int a = 0; a < 3; ++a)
              for (int b = 0; b < 3; ++b)
                conf.current().virial_tensor(a, b) += r(a)*(f(0)(b) +
                      f(1)(b) + f(2)(b) + f(3)(b));

            DEBUG(7, "\tatomic virial done");
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
        break;            
      }
      break;
    }
    case simulation::cgrain_func:
    {
      io::messages.add("Nonbonded_Innerloop",
          "no perturbed 1,4 interactions for Martini coarse-grained simulations!",
          io::message::critical);
      break;
    }
    default:
      io::messages.add("EDS-Nonbonded_Innerloop",
              "interaction function not implemented",
              io::message::critical);
      break;
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
        Periodicity_type const & periodicity, simulation::Simulation & sim) {

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

  bool multiAEDS = false;
  if (sim.param().eds.eds == 3){
    multiAEDS = true;
  }

  if (multiAEDS){
    //add self term to the <state_i,state_i, site_i, state_i> terms 
    for (unsigned int site_i = 0; site_i < sim.param().eds.numsites; site_i++){
      for (unsigned int state_i = 0; state_i < sim.param().eds.multnumstates[site_i]; state_i++){
        std::vector<int> states_i_i {site_i, state_i, site_i, state_i};
        r =0.0;
        double q_i = mit->second.M_charge()[state_i];
        eds_rf_interaction(r, q_i*q_i, f_rf, e_rf);
        DEBUG(7, "Self term for multiAEDS atom " << i << " in state " << state_i << " = " << e_rf << " (q*q = " << q_i * q_i << ")");
        conf.current().energies.eds_mult_vi[states_i_i] += 0.5 *e_rf;
      }
    }
  } else {
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
  }
  eds_perturbed_RF_exclusions_loop(topo, conf, i, mit->second.exclusion(), periodicity, sim);
}


template<typename t_interaction_spec,
typename t_perturbation_details>
inline void
interaction::Eds_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>
::perturbed_RF_excluded_interaction_innerloop
(topology::Topology & topo, configuration::Configuration & conf,
        std::map<unsigned int, topology::Perturbed_Atom>::const_iterator const & mit,
        Periodicity_type const & periodicity,  simulation::Simulation & sim
) {

  math::VArray &pos = conf.current().pos;
  math::VArray &force = conf.current().force;

  math::Vec r;
  double e_rf = 0.0, de_rf = 0.0;

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
              f_rf, e_rf, de_rf, true);
      // ANITA
      if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
        double A_e_rf = 0.0, B_e_rf = 0.0, A_de_rf = 0.0, B_de_rf = 0.0;

        // determine lambda stepsize from min,max and nr of lambdas
        double lambda_step = (sim.param().precalclam.max_lam -
                 sim.param().precalclam.min_lam) /
                 (sim.param().precalclam.nr_lambdas-1);

        //loop over nr_lambdas
        for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){

          // determine current lambda for this index
          double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;
          rf_soft_interaction_ext(r, q_i_a*q_i_a, q_i_b * q_i_b,
                  mit->second.CRF_softcore(), A_e_rf, B_e_rf, 
                  A_de_rf, B_de_rf, lam);
          conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(i)]
            [topo.atom_energy_group(i)] += 0.5 * A_e_rf;
          conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(i)]
            [topo.atom_energy_group(i)] += 0.5 * B_e_rf;
          conf.current().perturbed_energy_derivatives.A_crf_energy[lam_index]
            [topo.atom_energy_group(i)]
            [topo.atom_energy_group(i)] += 0.5 * A_de_rf;
          conf.current().perturbed_energy_derivatives.B_crf_energy[lam_index]
            [topo.atom_energy_group(i)]
            [topo.atom_energy_group(i)] += 0.5 * B_de_rf;

        }
      } // ANITA 
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
    case simulation::pol_lj_crf_func:
    case simulation::pol_off_lj_crf_func:
    {
      math::Vec rp1, rp2, rpp;

      //new introduced function to calculate the perturbed self reaction field term
      pol_rf_self_soft_interaction(q_i_a,q_i_b, e_rf, de_rf,true);

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
  eds_perturbed_RF_exclusions_loop(topo, conf, i, mit->second.exclusion(), periodicity, sim);
}


template<typename t_interaction_spec,
typename t_perturbation_details>
inline void
interaction::Eds_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>
::eds_perturbed_RF_exclusions_loop(topology::Topology & topo, configuration::Configuration & conf, int atom_i,
                                   const topology::excl_cont_t::value_type &exclusions, Periodicity_type const & periodicity,
                                   simulation::Simulation & sim){
  // Self term has already been accounted for
  // now loop over all the exclusions
  // those are fortunately not in the normal exclusions!
  math::VArray &pos = conf.current().pos;
  std::vector<math::VArray> &force_states = conf.special().eds.force_endstates;
  std::map<std::vector<int>, math::VArray> &force_mult_states = conf.special().eds.force_mult_endstates;

  math::VArray &force = conf.current().force;
  math::Vec r;

  std::vector<int>::const_iterator it, to;
  it = exclusions.begin();
  to = exclusions.end();

  r = 0.0;

  //do we use multiAEDS?
  bool multiAEDS = (sim.param().eds.eds == 3);

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
              math::Vec f_rf;
              double q_i = 0.0, q_j = 0.0, e_rf = 0.0;

              if (multiAEDS){
                const int &site_i = topo.eds_perturbed_solute().atoms()[atom_i].site_number();
                const int numstates_i = sim.param().eds.multnumstates[site_i];
                const int site_j = sim.param().eds.numsites; // non eds atom, belongs to site R (rest), last site
                const int state_j = 0; // non eds atom, just one state
                for (unsigned int state_i = 0; state_i < numstates_i; state_i++){
                  std::vector<int> states_i_j {site_i, state_i, site_j, state_j};
                  q_i = topo.eds_perturbed_solute().atoms()[atom_i].M_charge()[state_i];
                  q_j = topo.charge()(*it);
                  eds_rf_interaction(r,q_i*q_j, f_rf, e_rf);

                  DEBUG(7, "multiAEDS: excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

                  conf.current().energies.eds_mult_vi[states_i_j] += e_rf;
                  force_mult_states[states_i_j][atom_i] += f_rf;
                  force_mult_states[states_i_j][*it] -= f_rf;

                  for (int a = 0; a < 3; ++a){
                    for (int b = 0; b < 3; ++b){
                      conf.special().eds.virial_tensor_mult_endstates[states_i_j](a, b) += r(a) * f_rf(b);
                    }
                  }
                }
              } else {
                unsigned int numstates = conf.special().eds.force_endstates.size();

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
              }
              break;
            }
            case 1:
            {
              // EDS - EDS
              math::Vec f_rf;
              double q_i = 0.0, q_j = 0.0, e_rf = 0.0;

              if (multiAEDS){
                const int &site_i = topo.eds_perturbed_solute().atoms()[atom_i].site_number();
                const int numstates_i = sim.param().eds.multnumstates[site_i];
                const int &site_j =topo.eds_perturbed_solute().atoms()[*it].site_number();
                const int numstates_j = sim.param().eds.multnumstates[site_j];

                if (site_j >= site_i){
                  for (int state_i = 0; state_i < numstates_i; state_i++){
                    for (int state_j = 0; state_j < numstates_j; state_j++){
                      if ((site_j == site_i) && (state_i != state_j)){
                        continue; //when both atoms from same site, only include same states for both sites
                      }
                      std::vector<int> states_i_j {site_i, state_i, site_j, state_j};
                      q_i = topo.eds_perturbed_solute().atoms()[atom_i].M_charge()[state_i];
                      q_j = topo.eds_perturbed_solute().atoms()[*it].M_charge()[state_j];
                      eds_rf_interaction(r, q_i*q_j, f_rf, e_rf);
                      
                      DEBUG(7, "multiAEDS excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

                      // and add everything to the correct arrays
                      conf.current().energies.eds_mult_vi[states_i_j] += e_rf;

                      force_mult_states[states_i_j](atom_i) += f_rf;
                      force_mult_states[states_i_j](*it) -= f_rf;

                      // if (t_interaction_spec::do_virial != math::no_virial){
                      for (int a = 0; a < 3; ++a){
                        for (int b = 0; b < 3; ++b){
                          conf.special().eds.virial_tensor_mult_endstates[states_i_j](a, b) += r(a) * f_rf(b);
                        }
                      }
                      DEBUG(7, "\tatomic virial done");
                    }
                  }
                }
              } else{
                unsigned int numstates = conf.special().eds.force_endstates.size();
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
              }
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
              double A_q = 0.0, B_q = 0.0, e_rf = 0.0, de_rf = 0.0, alpha_crf = 0.0, q_ij_a = 0.0, q_ij_b = 0.0, q_i_a = 0.0,
              q_j_a = 0.0, q_i_b = 0.0, q_j_b = 0.0;

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
              q_j_a = topo.charge()(*it);
              q_ij_a = q_i_a * q_j_a;
              q_j_b = topo.charge()(*it);
              q_ij_b = q_i_b * q_j_b;

              eds_perturbed_rf_interaction(r, A_q, B_q, alpha_crf, f_rf, e_rf, de_rf);

              DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);

              // and add everything to the correct arrays

              conf.current().energies.crf_energy[topo.atom_energy_group(atom_i)][topo.atom_energy_group(*it)] += e_rf;
              conf.current().perturbed_energy_derivatives.crf_energy[topo.atom_energy_group(atom_i)][topo.atom_energy_group(*it)] += de_rf;

              force(atom_i) += f_rf;
              force(*it) -= f_rf;

              // ANITA
              if (sim.param().precalclam.nr_lambdas && ((sim.steps() % sim.param().write.free_energy) == 0)){
                double A_e_rf = 0.0, B_e_rf = 0.0, A_de_rf = 0.0, B_de_rf = 0.0;
      
                // determine lambda stepsize from min,max and nr of lambdas
                double lambda_step = (sim.param().precalclam.max_lam -
                        sim.param().precalclam.min_lam) /
                        (sim.param().precalclam.nr_lambdas-1);
      
                //loop over nr_lambdas
                for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
      
                  // determine current lambda for this index
                  double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;
                  rf_soft_interaction_ext(r, q_ij_a, q_ij_b, alpha_crf, A_e_rf, B_e_rf,
                          A_de_rf, B_de_rf, lam);

                  conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(atom_i)]
                    [topo.atom_energy_group(*it)] += A_e_rf;
                  conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(atom_i)]
                    [topo.atom_energy_group(*it)] += B_e_rf;
                  conf.current().perturbed_energy_derivatives.A_crf_energy[lam_index]
                    [topo.atom_energy_group(atom_i)]
                    [topo.atom_energy_group(*it)] += A_de_rf;
                  conf.current().perturbed_energy_derivatives.B_crf_energy[lam_index]
                    [topo.atom_energy_group(atom_i)]
                    [topo.atom_energy_group(*it)] += B_de_rf;
      
                }
              } // ANITA
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
              // Perturbed - Pertubed
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

              // ANITA
              if (sim.param().precalclam.nr_lambdas && ((sim.steps() % sim.param().write.free_energy) == 0)){
                double A_e_rf = 0.0, B_e_rf = 0.0, A_de_rf = 0.0, B_de_rf = 0.0;
      
                // determine lambda stepsize from min,max and nr of lambdas
                double lambda_step = (sim.param().precalclam.max_lam -
                        sim.param().precalclam.min_lam) /
                        (sim.param().precalclam.nr_lambdas-1);
      
                //loop over nr_lambdas
                for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
      
                  // determine current lambda for this index
                  double lam=(lam_index * lambda_step) + sim.param().precalclam.min_lam;
                  rf_soft_interaction_ext(r, A_q, B_q, alpha_crf, A_e_rf, B_e_rf,
                          A_de_rf, B_de_rf, lam);

                  conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(atom_i)]
                    [topo.atom_energy_group(*it)] += A_e_rf;
                  conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(atom_i)]
                    [topo.atom_energy_group(*it)] += B_e_rf;
                  conf.current().perturbed_energy_derivatives.A_crf_energy[lam_index]
                    [topo.atom_energy_group(atom_i)]
                    [topo.atom_energy_group(*it)] += A_de_rf;
                  conf.current().perturbed_energy_derivatives.B_crf_energy[lam_index]
                    [topo.atom_energy_group(atom_i)]
                    [topo.atom_energy_group(*it)] += B_de_rf;
      
                }
              } // ANITA
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
            case 3:
            // perturbed - normal
            {
              
              double q_ij_a = 0.0, q_i_a = 0.0, q_j_a = 0.0;
              double q_ij_b = 0.0, q_i_b = 0.0, q_j_b = 0.0;
              double e_rf = 0.0, de_rf = 0.0;
              double alpha_crf = 0.0;
              
              alpha_crf = topo.perturbed_solute().atoms()[atom_i].CRF_softcore();
              q_j_a = topo.charge()(*it);
              q_ij_a = q_i_a * q_j_a;
              q_j_b = topo.charge()(*it);
              q_ij_b = q_i_b * q_j_b;
              math::Vec f_rf;
              // check if...
              if (topo.is_coarse_grained(atom_i) && topo.is_coarse_grained(*it)) { // CG-CG interaction
                rf_soft_interaction(r, q_ij_a / cgrain_eps[0], q_ij_b / cgrain_eps[0],
                        alpha_crf, f_rf, e_rf, de_rf, 0);
              } else if (topo.is_coarse_grained(atom_i) || topo.is_coarse_grained(*it)) { // FG-CG interaction
                rf_soft_interaction(r, q_ij_a / cgrain_eps[1], q_ij_b / cgrain_eps[1],
                        alpha_crf, f_rf, e_rf, de_rf, 1);
              } else { // FG-FG interaction
                rf_soft_interaction(r, q_ij_a, q_ij_b, alpha_crf, f_rf, e_rf, de_rf, 2);
              }

              DEBUG(8, "alpha_crf : " << alpha_crf);
              DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);
              DEBUG(7, "\tde_rf: " << de_rf);

              // and add everything to the correct arrays
              conf.current().energies.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += e_rf;
              conf.current().perturbed_energy_derivatives.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += de_rf;
              force(atom_i) += f_rf;
              force(*it) -= f_rf;

              // if (t_interaction_spec::do_virial != math::no_virial){
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                  conf.current().virial_tensor(a, b) +=
                        r(a) * f_rf(b);

              DEBUG(7, "\tatomic virial done");
              break;
            }
            case 5:
            // perturbed - perturbed
            {
              
              double q_ij_a = 0.0, q_i_a = 0.0, q_j_a = 0.0;
              double q_ij_b = 0.0, q_i_b = 0.0, q_j_b = 0.0;
              double e_rf = 0.0, de_rf = 0.0;
              double alpha_crf = 0.0;

              alpha_crf = (topo.perturbed_solute().atoms()[atom_i].CRF_softcore() +
                            topo.perturbed_solute().atoms()[*it].CRF_softcore()) / 2.0;
              q_j_a = topo.perturbed_solute().atoms()[*it].A_charge();
              q_ij_a = q_i_a * q_j_a;
              q_j_b = topo.perturbed_solute().atoms()[*it].B_charge();
              q_ij_b = q_i_b * q_j_b;
              math::Vec f_rf;
              // check if...
              if (topo.is_coarse_grained(atom_i) && topo.is_coarse_grained(*it)) { // CG-CG interaction
                rf_soft_interaction(r, q_ij_a / cgrain_eps[0], q_ij_b / cgrain_eps[0],
                        alpha_crf, f_rf, e_rf, de_rf, 0);
              } else if (topo.is_coarse_grained(atom_i) || topo.is_coarse_grained(*it)) { // FG-CG interaction
                rf_soft_interaction(r, q_ij_a / cgrain_eps[1], q_ij_b / cgrain_eps[1],
                        alpha_crf, f_rf, e_rf, de_rf, 1);
              } else { // FG-FG interaction
                rf_soft_interaction(r, q_ij_a, q_ij_b, alpha_crf, f_rf, e_rf, de_rf, 2);
              }

              DEBUG(8, "alpha_crf : " << alpha_crf);
              DEBUG(7, "excluded atoms " << atom_i << " & " << *it << ": " << e_rf);
              DEBUG(7, "\tde_rf: " << de_rf);

              // and add everything to the correct arrays
              conf.current().energies.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += e_rf;
              conf.current().perturbed_energy_derivatives.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += de_rf;
              force(atom_i) += f_rf;
              force(*it) -= f_rf;

              // if (t_interaction_spec::do_virial != math::no_virial){
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                  conf.current().virial_tensor(a, b) +=
                        r(a) * f_rf(b);

              DEBUG(7, "\tatomic virial done");
              break;
            }
            default:
                  io::messages.add("EDS_Nonbonded_Innerloop",
                  "rf excluded interaction function not implemented",
                  io::message::critical);
                  break;
          }
          break;
        }
        case simulation::pol_lj_crf_func:
        {
          switch (both_perturbed){
            case 3:
            // Perturbed - Normal
            {
              double q_ij_a = 0.0, q_i_a = 0.0, q_j_a = 0.0;
              double q_ij_b = 0.0, q_i_b = 0.0, q_j_b = 0.0;
              double e_rf = 0.0, de_rf = 0.0;
              double alpha_crf = 0.0;
              
              alpha_crf = topo.perturbed_solute().atoms()[atom_i].CRF_softcore();
              q_j_a = topo.charge()(*it);
              q_ij_a = q_i_a * q_j_a;
              q_j_b = topo.charge()(*it);
              q_ij_b = q_i_b * q_j_b;

              math::Vec rp1, rp2, rpp;
              double f_rf[4];
              rp1 = r - conf.current().posV(*it);
              rp2 = r + conf.current().posV(atom_i);
              rpp = r + conf.current().posV(atom_i) - conf.current().posV(*it);

              pol_rf_soft_interaction(r, rp1, rp2, rpp,
                      q_i_a, q_j_a, q_i_b, q_j_b,
                      topo.coscharge(atom_i), topo.coscharge(*it),
                      alpha_crf, f_rf, e_rf, de_rf);

              // and add everything to the correct arrays
              conf.current().energies.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += e_rf;
              conf.current().perturbed_energy_derivatives.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += de_rf;
              force(atom_i) += f_rf[0] * r + f_rf[1] * rp1 + f_rf[2] * rp2 + f_rf[3] * rpp;
              force(*it) -= f_rf[0] * r + f_rf[1] * rp1 + f_rf[2] * rp2 + f_rf[3] * rpp;

              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                  conf.current().virial_tensor(a, b) +=
                        r(a)*(f_rf[0] * r(b) + f_rf[1] * rp1(b) +
                        f_rf[2] * rp2(b) + f_rf[3] * rpp(b));
              break;
            }
            case 5:
            // Perturbed - Perturbed
            {
              double q_ij_a = 0.0, q_i_a = 0.0, q_j_a = 0.0;
              double q_ij_b = 0.0, q_i_b = 0.0, q_j_b = 0.0;
              double e_rf = 0.0, de_rf = 0.0;
              double alpha_crf = 0.0;

              alpha_crf = (topo.perturbed_solute().atoms()[atom_i].CRF_softcore() +
                            topo.perturbed_solute().atoms()[*it].CRF_softcore()) / 2.0;
              q_j_a = topo.perturbed_solute().atoms()[*it].A_charge();
              q_ij_a = q_i_a * q_j_a;
              q_j_b = topo.perturbed_solute().atoms()[*it].B_charge();
              q_ij_b = q_i_b * q_j_b;

              math::Vec rp1, rp2, rpp;
              double f_rf[4];
              rp1 = r - conf.current().posV(*it);
              rp2 = r + conf.current().posV(atom_i);
              rpp = r + conf.current().posV(atom_i) - conf.current().posV(*it);

              pol_rf_soft_interaction(r, rp1, rp2, rpp,
                      q_i_a, q_j_a, q_i_b, q_j_b,
                      topo.coscharge(atom_i), topo.coscharge(*it),
                      alpha_crf, f_rf, e_rf, de_rf);

              // and add everything to the correct arrays
              conf.current().energies.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += e_rf;
              conf.current().perturbed_energy_derivatives.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += de_rf;
              force(atom_i) += f_rf[0] * r + f_rf[1] * rp1 + f_rf[2] * rp2 + f_rf[3] * rpp;
              force(*it) -= f_rf[0] * r + f_rf[1] * rp1 + f_rf[2] * rp2 + f_rf[3] * rpp;

              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                  conf.current().virial_tensor(a, b) +=
                        r(a)*(f_rf[0] * r(b) + f_rf[1] * rp1(b) +
                        f_rf[2] * rp2(b) + f_rf[3] * rpp(b));
              break;
            }
            default:
              io::messages.add("EDS_Nonbonded_Innerloop",
              "rf excluded interaction function not implemented",
              io::message::critical);
              break;
          }
          break;
        }
        case simulation::pol_off_lj_crf_func:
        {
          switch (both_perturbed){
            case 3:
            {
              double q_ij_a = 0.0, q_i_a = 0.0, q_j_a = 0.0;
              double q_ij_b = 0.0, q_i_b = 0.0, q_j_b = 0.0;
              double e_rf = 0.0, de_rf = 0.0;
              double alpha_crf = 0.0;
              
              alpha_crf = topo.perturbed_solute().atoms()[atom_i].CRF_softcore();
              q_j_a = topo.charge()(*it);
              q_ij_a = q_i_a * q_j_a;
              q_j_b = topo.charge()(*it);
              q_ij_b = q_i_b * q_j_b;

              math::Vec rm = r;
              if (topo.gamma(atom_i)!=0.0) {
                math::Vec rij, rik;
                periodicity.nearest_image(conf.current().pos(atom_i),
                        conf.current().pos(topo.gamma_j(atom_i)), rij);
                periodicity.nearest_image(conf.current().pos(atom_i),
                        conf.current().pos(topo.gamma_k(atom_i)), rik);
                rm -= topo.gamma(atom_i)*(rij + rik) / 2;
              }
              if (topo.gamma(*it)!=0.0) {
                math::Vec rjj, rjk;
                periodicity.nearest_image(conf.current().pos(*it),
                        conf.current().pos(topo.gamma_k(*it)), rjk);
                periodicity.nearest_image(conf.current().pos(*it),
                        conf.current().pos(topo.gamma_j(*it)), rjj);

                rm += topo.gamma(*it)*(rjj + rjk) / 2;
              }

              math::Vec rp1, rp2, rpp;
              double f_rf[4];
              rp1 = rm - conf.current().posV(*it);
              rp2 = rm + conf.current().posV(atom_i);
              rpp = rm + conf.current().posV(atom_i) - conf.current().posV(*it);

              pol_rf_soft_interaction(rm, rp1, rp2, rpp,
                      q_i_a, q_j_a, q_i_b, q_j_b,
                      topo.coscharge(atom_i), topo.coscharge(*it),
                      alpha_crf, f_rf, e_rf, de_rf);

              // and add everything to the correct arrays
              conf.current().energies.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += e_rf;
              conf.current().perturbed_energy_derivatives.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += de_rf;
              for (int a = 0; a < 3; ++a) {
                const double term = f_rf[0] * rm(a) + f_rf[1] * rp1(a)
                        + f_rf[2] * rp2(a) + f_rf[3] * rpp(a);

                force(atom_i)(a) += (1 - topo.gamma(atom_i)) * term;
                force(*it)(a) -= (1 - topo.gamma(*it)) * term;

                force(topo.gamma_j(atom_i))(a) += topo.gamma(atom_i) / 2 * term;
                force(topo.gamma_j(*it))(a) -= topo.gamma(*it) / 2 * term;

                force(topo.gamma_k(atom_i))(a) += topo.gamma(atom_i) / 2 * term;
                force(topo.gamma_k(*it))(a) -= topo.gamma(*it) / 2 * term;

                for (int b = 0; b < 3; ++b)
                  conf.current().virial_tensor(b, a) += r(b) * term;
              }
              break;
            }
            case 5:
            {
              double q_ij_a = 0.0, q_i_a = 0.0, q_j_a = 0.0;
              double q_ij_b = 0.0, q_i_b = 0.0, q_j_b = 0.0;
              double e_rf = 0.0, de_rf = 0.0;
              double alpha_crf = 0.0;

              alpha_crf = (topo.perturbed_solute().atoms()[atom_i].CRF_softcore() +
                            topo.perturbed_solute().atoms()[*it].CRF_softcore()) / 2.0;
              q_j_a = topo.perturbed_solute().atoms()[*it].A_charge();
              q_ij_a = q_i_a * q_j_a;
              q_j_b = topo.perturbed_solute().atoms()[*it].B_charge();
              q_ij_b = q_i_b * q_j_b;

              math::Vec rm = r;
              if (topo.gamma(atom_i)!=0.0) {
                math::Vec rij, rik;
                periodicity.nearest_image(conf.current().pos(atom_i),
                        conf.current().pos(topo.gamma_j(atom_i)), rij);
                periodicity.nearest_image(conf.current().pos(atom_i),
                        conf.current().pos(topo.gamma_k(atom_i)), rik);
                rm -= topo.gamma(atom_i)*(rij + rik) / 2;
              }
              if (topo.gamma(*it)!=0.0) {
                math::Vec rjj, rjk;
                periodicity.nearest_image(conf.current().pos(*it),
                        conf.current().pos(topo.gamma_k(*it)), rjk);
                periodicity.nearest_image(conf.current().pos(*it),
                        conf.current().pos(topo.gamma_j(*it)), rjj);

                rm += topo.gamma(*it)*(rjj + rjk) / 2;
              }

              math::Vec rp1, rp2, rpp;
              double f_rf[4];
              rp1 = rm - conf.current().posV(*it);
              rp2 = rm + conf.current().posV(atom_i);
              rpp = rm + conf.current().posV(atom_i) - conf.current().posV(*it);

              pol_rf_soft_interaction(rm, rp1, rp2, rpp,
                      q_i_a, q_j_a, q_i_b, q_j_b,
                      topo.coscharge(atom_i), topo.coscharge(*it),
                      alpha_crf, f_rf, e_rf, de_rf);

              // and add everything to the correct arrays
              conf.current().energies.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += e_rf;
              conf.current().perturbed_energy_derivatives.crf_energy
                      [topo.atom_energy_group(atom_i)]
                      [topo.atom_energy_group(*it)] += de_rf;
              for (int a = 0; a < 3; ++a) {
                const double term = f_rf[0] * rm(a) + f_rf[1] * rp1(a)
                        + f_rf[2] * rp2(a) + f_rf[3] * rpp(a);

                force(atom_i)(a) += (1 - topo.gamma(atom_i)) * term;
                force(*it)(a) -= (1 - topo.gamma(*it)) * term;

                force(topo.gamma_j(atom_i))(a) += topo.gamma(atom_i) / 2 * term;
                force(topo.gamma_j(*it))(a) -= topo.gamma(*it) / 2 * term;

                force(topo.gamma_k(atom_i))(a) += topo.gamma(atom_i) / 2 * term;
                force(topo.gamma_k(*it))(a) -= topo.gamma(*it) / 2 * term;

                for (int b = 0; b < 3; ++b)
                  conf.current().virial_tensor(b, a) += r(b) * term;
              }
              break;
            }
            break;
          default:
              io::messages.add("EDS_Nonbonded_Innerloop",
              "rf excluded interaction function not implemented",
              io::message::critical);
              break;
          }
          break;
        }
        default:
        io::messages.add("EDS_Nonbonded_Innerloop",
        "rf excluded interaction function not implemented",
        io::message::critical);
        break;
    }
  } // loop over all atoms
}

template<typename t_interaction_spec, typename t_perturbation_details >
inline void interaction::Eds_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>::perturbed_electric_field_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i, unsigned int j, math::Vec &e_eli, math::Vec &e_elj,
        Periodicity_type const & periodicity
        ) {
  math::Vec r, rp1, rp2, e_el1, e_el2;
  double alpha_crf = 0.0, A_qi = 0.0, A_qj = 0.0, B_qi = 0.0, B_qj = 0.0;

  if (topo.is_perturbed(j) == true) {
    // both i and j are perturbed
    A_qi = topo.perturbed_solute().atoms()[i].A_charge();
    A_qj = topo.perturbed_solute().atoms()[j].A_charge();
    B_qi = topo.perturbed_solute().atoms()[i].B_charge();
    B_qj = topo.perturbed_solute().atoms()[j].B_charge();

    alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
            topo.perturbed_solute().atoms()[j].CRF_softcore()) / 2.0;
  } else {
    // only i is perturbed
    A_qi = topo.perturbed_solute().atoms()[i].A_charge();
    A_qj = topo.charge()(j);
    B_qi = topo.perturbed_solute().atoms()[i].B_charge();
    B_qj = topo.charge()(j);

    alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();
  }

  // energy field term at position i and j
  DEBUG(11, "\tenergy field calculation i: " << i << " j: " << j);

  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);

  int n1 = topo.atom_energy_group(i);
  int n2 = topo.atom_energy_group(j);

  // calculate all the lambda values for lj and crf interactions, softness and A and B
  set_lambda(topo.individual_lambda(simulation::lj_lambda)[n1][n2],
          topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2],
          topo.individual_lambda(simulation::crf_lambda)[n1][n2],
          topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2],
          topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2],
          topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2],
          topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2],
          topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2],
          topo.lambda_exp());

  if (topo.gamma(i) != 0.0 && simulation::pol_off_lj_crf_func != 0)
  {
     math::Vec rij, rik, rim;
     periodicity.nearest_image(conf.current().pos(i),
       conf.current().pos(topo.gamma_j(i)), rij); 
     periodicity.nearest_image(conf.current().pos(i),
       conf.current().pos(topo.gamma_k(i)), rik);
     rim = topo.gamma(i)*(rij + rik) / 2;
     r -= rim;
  }
  if (topo.gamma(j) != 0.0 && simulation::pol_off_lj_crf_func != 0) 
  {
     math::Vec rjj, rjk, rjm;
     periodicity.nearest_image(conf.current().pos(j),
       conf.current().pos(topo.gamma_j(j)), rjj); 
     periodicity.nearest_image(conf.current().pos(j),
       conf.current().pos(topo.gamma_k(j)), rjk);   
     rjm = topo.gamma(j)*(rjj + rjk) / 2; 
     r += rjm;
  }
  switch (t_interaction_spec::efield_site) {
    case simulation::ef_atom:
    {

      rp1 = r - conf.current().posV(j);
      rp2 = -r - conf.current().posV(i);

      electric_field_soft_interaction(r, rp1, alpha_crf, A_qj, B_qj,
              topo.coscharge(j), e_el1);
      electric_field_soft_interaction(-r, rp2, alpha_crf, A_qi, B_qi,
              topo.coscharge(i), e_el2);
      break;
    }
    case simulation::ef_cos:
    {
      const math::Vec r_cos1 = conf.current().posV(i) + r,
              r_cos2 = conf.current().posV(j) - r;
      rp1 = r_cos1 - conf.current().posV(j);
      rp2 = r_cos2 - conf.current().posV(i);

      electric_field_soft_interaction(r_cos1, rp1, alpha_crf, A_qj, B_qj,
              topo.coscharge(j), e_el1);
      electric_field_soft_interaction(r_cos2, rp2, alpha_crf, A_qi, B_qi,
              topo.coscharge(i), e_el2);
      break;
    }
    default:
      io::messages.add("electric_field_innerloop",
              "this site for electric field calculation was not implemented.",
              io::message::critical);

  }

  e_eli = e_el1;
  e_elj = e_el2;

}

template<typename t_interaction_spec, typename t_perturbation_details>
inline void
interaction::Eds_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>::perturbed_self_energy_innerloop
(
        topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i,
        Storage & storage,
        Periodicity_type const & periodicity
        ) {
  DEBUG(8, "\tself energy of molecule i " << i);

  double self_e = 0.0, de_self = 0.0;
  const double e_i2 = math::abs2(storage.electric_field(i));
  const double alpha_A = topo.perturbed_solute().atom(i).A_polarisability(),
          alpha_B = topo.perturbed_solute().atom(i).B_polarisability(),
          e0_A = topo.perturbed_solute().atom(i).A_damping_level(),
          e0_B = topo.perturbed_solute().atom(i).B_damping_level();

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

  if (t_interaction_spec::pol_damping)
    self_energy_soft_interaction(alpha_A, alpha_B, e_i2,
          e0_A, e0_B, topo.damping_power(i),
          self_e, de_self);
  else
    self_energy_soft_interaction(alpha_A, alpha_B, e_i2, self_e, de_self);


  // self energy term
  DEBUG(10, "\tself energy storage:\t" << self_e);
  storage.energies.self_energy[topo.atom_energy_group(i)] += self_e;

  DEBUG(11, "\tself energy derivative storage:\t" << de_self);
  storage.perturbed_energy_derivatives.self_energy[topo.atom_energy_group(i)]
          += de_self;
}
