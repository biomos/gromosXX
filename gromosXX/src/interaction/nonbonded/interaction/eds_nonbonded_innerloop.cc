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
        ) {
  DEBUG(8, "\teds-perturbed pair\t" << i << "\t" << j << " (inner loop)");

  math::Vec r;
  const unsigned int numstates = conf.special().eds.force_endstates.size();
  assert(storage.energies.eds_vi.size() == numstates);
  assert(storage.force_endstates.size() == numstates);

  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);

  // speed!
  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
  // the following increases topo.eds_perturbed_solute().atoms().size() by 1
  // PROBLEM: j can also be not perturbed!!!
  /*
  const std::vector<unsigned int> &pert_j_M_IAC = topo.eds_perturbed_solute().atoms()[j].M_IAC();
  const std::vector<double> &pert_j_M_charge = topo.eds_perturbed_solute().atoms()[j].M_charge();
   */
  assert(pert_i_M_charge.size() == numstates);
  const std::vector<std::vector<lj_parameter_struct> > &m_param_lj_parameter = m_param->lj_parameter();

  unsigned int iac_j = topo.iac(j);
  double charge_j = topo.charge()(j);
  // --speed

  double alpha_lj = 0, alpha_crf = 0;

  int both_perturbed = 0;
  if (j < topo.num_solute_atoms() &&
          topo.is_eds_perturbed(j) == true) {
    both_perturbed = 1;
  }

  switch (t_interaction_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      std::vector<math::VArray> &force_endstates = storage.force_endstates;
      double c6 = 0.0, c12 = 0.0, q = 0.0, e_nb = 0.0, f = 0.0;
      assert(abs2(r) != 0);
      const double dist2 = abs2(r);
      const double dist6 = dist2 * dist2 * dist2;
      switch (both_perturbed) {
        case 0:
        {
          alpha_lj = topo.eds_perturbed_solute().atoms()[i].LJ_softcore();
          alpha_crf = topo.eds_perturbed_solute().atoms()[i].CRF_softcore();
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

            eds_lj_crf_interaction(dist2, dist6, c6, c12, q, alpha_lj,
                    alpha_crf, f, e_nb);

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
          alpha_lj = (topo.eds_perturbed_solute().atoms()[i].LJ_softcore() +
                  topo.eds_perturbed_solute().atoms()[j].LJ_softcore()) /
                  2.0;
          alpha_crf = (topo.eds_perturbed_solute().atoms()[i].CRF_softcore() +
                  topo.eds_perturbed_solute().atoms()[j].CRF_softcore()) /
                  2.0;
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
            eds_lj_crf_interaction(dist2, dist6, c6, c12, q,
                    alpha_lj, alpha_crf, f, e_nb);

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
        const double dist6 = dist2 * dist2 * dist2;
        switch (both_perturbed) {
          case 0:
          {
            alpha_lj = topo.eds_perturbed_solute().atoms()[i].LJ_softcore();
            alpha_crf = topo.eds_perturbed_solute().atoms()[i].CRF_softcore();
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

              eds_lj_crf_interaction(dist2, dist6, c6, c12, q / cgrain_eps[1],
                      alpha_lj, alpha_crf, f, e_nb, 1);

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
        const double dist6 = dist2 * dist2 * dist2;
        switch (both_perturbed) {
          case 0:
          {
            alpha_lj = topo.eds_perturbed_solute().atoms()[i].LJ_softcore();
            alpha_crf = topo.eds_perturbed_solute().atoms()[i].CRF_softcore();
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

              eds_lj_crf_interaction(dist2, dist6, c6, c12, q,
                      alpha_lj, alpha_crf, f, e_nb, 2);

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
            alpha_lj = (topo.eds_perturbed_solute().atoms()[i].LJ_softcore() +
                    topo.eds_perturbed_solute().atoms()[j].LJ_softcore()) /
                    2.0;
            alpha_crf = (topo.eds_perturbed_solute().atoms()[i].CRF_softcore() +
                    topo.eds_perturbed_solute().atoms()[j].CRF_softcore()) /
                    2.0;
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
              eds_lj_crf_interaction(dist2, dist6, c6, c12, q,
                      alpha_lj, alpha_crf, f, e_nb, 2);

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

  math::Vec r;
  const unsigned int numstates = conf.special().eds.force_endstates.size();
  assert(conf.current().energies.eds_vi.size() == numstates);
  assert(conf.special().eds.force_endstates.size() == numstates);


  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);

  // speed!
  const std::vector<unsigned int> &pert_i_M_IAC = topo.eds_perturbed_solute().atoms()[i].M_IAC();
  const std::vector<double> &pert_i_M_charge = topo.eds_perturbed_solute().atoms()[i].M_charge();
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

  int both_perturbed = 0;
  if (j < topo.num_solute_atoms() &&
          topo.is_eds_perturbed(j) == true) {
    both_perturbed = 1;
  }

  switch (t_interaction_spec::interaction_func) {
    case simulation::lj_crf_func:
    {
      double c6 = 0.0, c12 = 0.0, q = 0.0, e_nb = 0.0, f = 0.0;
      assert(abs2(r) != 0);
      const double dist2 = abs2(r);
      const double dist6 = dist2 * dist2 * dist2;
      switch (both_perturbed) {
        case 0:
        {
          alpha_lj = topo.eds_perturbed_solute().atoms()[i].LJ_softcore();
          alpha_crf = topo.eds_perturbed_solute().atoms()[i].CRF_softcore();
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

            eds_lj_crf_interaction(dist2, dist6, c6, c12, q, alpha_lj,
                    alpha_crf, f, e_nb, 0, m_param->get_coulomb_scaling());

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
          alpha_lj = (topo.eds_perturbed_solute().atoms()[i].LJ_softcore() +
                  topo.eds_perturbed_solute().atoms()[j].LJ_softcore()) /
                  2.0;
          alpha_crf = (topo.eds_perturbed_solute().atoms()[i].CRF_softcore() +
                  topo.eds_perturbed_solute().atoms()[j].CRF_softcore()) /
                  2.0;
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
            eds_lj_crf_interaction(dist2, dist6, c6, c12, q, alpha_lj,
                    alpha_crf, f, e_nb, 0, m_param->get_coulomb_scaling());

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
        const double dist6 = dist2 * dist2 * dist2;
        switch (both_perturbed) {
          case 0:
          {
            alpha_lj = topo.eds_perturbed_solute().atoms()[i].LJ_softcore();
            alpha_crf = topo.eds_perturbed_solute().atoms()[i].CRF_softcore();
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

              eds_lj_crf_interaction(dist2, dist6, c6, c12, q, alpha_lj,
                      alpha_crf, f, e_nb, 2);

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
            alpha_lj = (topo.eds_perturbed_solute().atoms()[i].LJ_softcore() +
                    topo.eds_perturbed_solute().atoms()[j].LJ_softcore()) /
                    2.0;
            alpha_crf = (topo.eds_perturbed_solute().atoms()[i].CRF_softcore() +
                    topo.eds_perturbed_solute().atoms()[j].CRF_softcore()) /
                    2.0;
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
              eds_lj_crf_interaction(dist2, dist6, c6, c12, q, alpha_lj,
                      alpha_crf, f, e_nb, 2);

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

  double alpha_crf = mit->second.CRF_softcore();

  topology::excl_cont_t::value_type::const_iterator it, to;

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
      eds_rf_interaction(r, topo.charge()(i) * topo.charge()(i), alpha_crf, f_rf, e_rf);
      break;
    }
    case simulation::cggromos_func:
    {
      if (topo.is_coarse_grained(i)) { // CG-CG or CG-FG
        eds_rf_interaction(r, topo.charge()(i) * topo.charge()(i), alpha_crf, f_rf, e_rf, 0);
      } else { // FG-FG
        eds_rf_interaction(r, topo.charge()(i) * topo.charge()(i), alpha_crf, f_rf, e_rf, 2);
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
        eds_rf_interaction(r, q_i*q_i, alpha_crf, f_rf, e_rf);
        break;
      }
      case simulation::cggromos_func:
      {
        if (topo.is_coarse_grained(i)) { // CG-CG or CG-FG
          eds_rf_interaction(r, q_i*q_i / cgrain_eps[0], alpha_crf, f_rf, e_rf, 0);
        } else { // FG-FG
          eds_rf_interaction(r, q_i*q_i, alpha_crf, f_rf, e_rf, 2);
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

    // now loop over the exclusions
    // those are fortunately not in the normal exclusions!
    it = mit->second.exclusion().begin();
    to = mit->second.exclusion().end();
    //it = topo.exclusion(i).begin();
    //to = topo.exclusion(i).end();

    for (; it != to; ++it) {
      periodicity.nearest_image(pos(i), pos(*it), r);

      DEBUG(8, "r2 i(" << i << "-" << *it << ") " << abs2(r));

      double q_j = 0.0;

      if (unsigned(*it) < topo.num_solute_atoms() && topo.is_eds_perturbed(*it)) {
        // j perturbed
        q_j = topo.eds_perturbed_solute().atoms()[*it].M_charge()[state];
      } else {
        // only i perturbed
        q_j = topo.charge()(*it);
      }

      switch (t_interaction_spec::interaction_func) {
        case simulation::lj_crf_func:
        {
          math::Vec f_rf;
          eds_rf_interaction(r, q_i*q_j, alpha_crf, f_rf, e_rf);

          DEBUG(7, "excluded atoms " << i << " & " << *it << ": " << e_rf);

          // and add everything to the correct arrays
          conf.current().energies.eds_vi[state] += e_rf;

          force[state](i) += f_rf;
          force[state](*it) -= f_rf;

          // if (t_interaction_spec::do_virial != math::no_virial){
          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
              conf.special().eds.virial_tensor_endstates[state](a, b) +=
                    r(a) * f_rf(b);

          DEBUG(7, "\tatomic virial done");
          // }
          break;
        }
        case simulation::cggromos_func:
        {
          math::Vec f_rf;
          if (topo.is_coarse_grained(i) && topo.is_coarse_grained(*it) ) { // CG-CG
            eds_rf_interaction(r, q_i*q_j / cgrain_eps[0], alpha_crf, f_rf, e_rf, 0);
          } else if (topo.is_coarse_grained(i) || topo.is_coarse_grained(*it)) { // FG-CG
            eds_rf_interaction(r, q_i*q_j / cgrain_eps[1], alpha_crf, f_rf, e_rf, 1);
          } else { // FG-FG
            eds_rf_interaction(r, q_i*q_j, alpha_crf, f_rf, e_rf, 2);
          }

          DEBUG(7, "excluded atoms " << i << " & " << *it << ": " << e_rf);

          // and add everything to the correct arrays
          conf.current().energies.eds_vi[state] += e_rf;

          force[state](i) += f_rf;
          force[state](*it) -= f_rf;

          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
              conf.special().eds.virial_tensor_endstates[state](a, b) +=
                    r(a) * f_rf(b);

          DEBUG(7, "\tatomic virial done");
          
          break;
        }
        default:
          io::messages.add("EDS_Nonbonded_Innerloop",
                  "rf excluded interaction function not implemented",
                  io::message::critical);
      }
    }
  } // loopover states
}

