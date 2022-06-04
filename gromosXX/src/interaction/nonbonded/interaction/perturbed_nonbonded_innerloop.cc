/**
 * @file perturbed_nonbonded_innerloop.cc
 * template methods of Perturbed_Nonbonded_Innerloop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

template<typename t_interaction_spec,
typename t_perturbation_details>
void interaction::Perturbed_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>
::perturbed_lj_crf_innerloop
(
        topology::Topology & topo, configuration::Configuration & conf,
        unsigned int i, unsigned int j, Storage &storage,
        Periodicity_type const & periodicity,
        // ANITA 
        simulation::Simulation & sim //
        ) {
  DEBUG(8, "\tperturbed pair\t" << i << "\t" << j << " (inner loop)");
  math::Vec r;
  double e_lj = 0.0, e_crf = 0.0, de_lj = 0.0, de_crf = 0.0;

  int energy_derivative_index = -1;

  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);

  lj_parameter_struct const *A_lj = nullptr;
  lj_parameter_struct const *B_lj = nullptr;
  double A_q = 0.0, B_q = 0.0, A_qi = 0.0, A_qj = 0.0, B_qi = 0.0, B_qj = 0.0;

  double alpha_lj = 0, alpha_crf = 0;

  // const double l = topo.lambda();

  if (j < topo.num_solute_atoms() && topo.is_perturbed(j) == true) {

    switch (t_interaction_spec::interaction_func) {
      case simulation::lj_crf_func:
      case simulation::cggromos_func:
      case simulation::pol_lj_crf_func:
      case simulation::pol_off_lj_crf_func:
        A_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
                topo.perturbed_solute().atoms()[j].A_IAC());
        B_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
                topo.perturbed_solute().atoms()[j].B_IAC());
        break;
      case simulation::cgrain_func:
        A_lj = &m_param->cg_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
                topo.perturbed_solute().atoms()[j].A_IAC());
        B_lj = &m_param->cg_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
                topo.perturbed_solute().atoms()[j].B_IAC());
        break;
      default:
        io::messages.add("Nonbonded_Innerloop",
                "interaction function not implemented",
                io::message::critical);
    }

    A_qi = topo.perturbed_solute().atoms()[i].A_charge();
    A_qj = topo.perturbed_solute().atoms()[j].A_charge();
    B_qi = topo.perturbed_solute().atoms()[i].B_charge();
    B_qj = topo.perturbed_solute().atoms()[j].B_charge();
    A_q = A_qi * A_qj;
    B_q = B_qi * B_qj;

    alpha_lj = (topo.perturbed_solute().atoms()[i].LJ_softcore() +
            topo.perturbed_solute().atoms()[j].LJ_softcore()) /
            2.0;
    alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
            topo.perturbed_solute().atoms()[j].CRF_softcore()) /
            2.0;

  } else {
    switch (t_interaction_spec::interaction_func) {
      case simulation::lj_crf_func:
      case simulation::cggromos_func:
      case simulation::pol_lj_crf_func:
      case simulation::pol_off_lj_crf_func:
        A_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
                topo.iac(j));
        B_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
                topo.iac(j));
        break;
      case simulation::cgrain_func:
        A_lj = &m_param->cg_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
                topo.iac(j));
        B_lj = &m_param->cg_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
                topo.iac(j));
        break;
      default:
        io::messages.add("Nonbonded_Innerloop",
                "interaction function not implemented",
                io::message::critical);
        return;
    }
    DEBUG(10, "\tiac-i (A) : " << topo.perturbed_solute().atoms()[i].A_IAC()
            << " iac-i (B) : " << topo.perturbed_solute().atoms()[i].B_IAC()
            << " iac-j : " << topo.iac(j));

    A_qi = topo.perturbed_solute().atoms()[i].A_charge();
    A_qj = topo.charge()(j);
    B_qi = topo.perturbed_solute().atoms()[i].B_charge();
    B_qj = topo.charge()(j);
    A_q = A_qi * A_qj;
    B_q = B_qi * B_qj;

    alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
    alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

  }

  DEBUG(8, "\tlj-parameter state A c6=" << A_lj->c6 << " c12=" << A_lj->c12);
  DEBUG(8, "\tlj-parameter state B c6=" << B_lj->c6 << " c12=" << B_lj->c12);
  DEBUG(8, "\tcharges state A i*j = " << A_q);
  DEBUG(8, "\tcharges state B i*j = " << B_q);
  DEBUG(8, "\talpha lj = " << alpha_lj);
  DEBUG(8, "\talpha crf = " << alpha_crf);

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

  DEBUG(8, "individual lambdas\nlj_lambda\t"
          << topo.individual_lambda(simulation::lj_lambda)[n1][n2]
          << "\nlj_softness_lambda\t"
          << topo.individual_lambda(simulation::lj_softness_lambda)[n1][n2]
          << "\ncrf_lambda\t"
          << topo.individual_lambda(simulation::crf_lambda)[n1][n2]
          << "\ncrf_softness_lambda\t"
          << topo.individual_lambda(simulation::crf_softness_lambda)[n1][n2]
          << "\nlj_lambda_derivative\t"
          << topo.individual_lambda_derivative(simulation::lj_lambda)[n1][n2]
          << "\nlj_softness_lambda_derivative\t"
          << topo.individual_lambda_derivative(simulation::lj_softness_lambda)[n1][n2]
          << "\ncrf_lambda_derivative\t"
          << topo.individual_lambda_derivative(simulation::crf_lambda)[n1][n2]
          << "\ncrf_softness_lambda_derivative\t"
          << topo.individual_lambda_derivative(simulation::crf_softness_lambda)[n1][n2]);

  if (t_perturbation_details::do_scaling) {
    // SCALING ON
    math::Vec f;
    double f1 = 0.0, f6 = 0.0, f12 = 0.0;
    if (t_interaction_spec::interaction_func == simulation::cgrain_func ||
            t_interaction_spec::interaction_func == simulation::cggromos_func) {
      io::messages.add("Nonbonded Innerloop",
              "scaling not implemented for coarse-grained simulations!",
              io::message::critical);
      return;
    }
    if (t_interaction_spec::interaction_func == simulation::pol_lj_crf_func
            || t_interaction_spec::interaction_func == simulation::pol_off_lj_crf_func) {
      io::messages.add("Nonbonded Innerloop",
              "scaling not implemented for COS polarisation simulations!",
              io::message::critical);
      return;
    }

    // check whether we need to do scaling
    // based on energy groups
    std::pair<int, int> energy_group_pair(topo.atom_energy_group(i),
            topo.atom_energy_group(j));
    
    if (topo.energy_group_scaling().count(energy_group_pair)) {

      // YES, we do scale the interactions!
      lj_crf_scaled_interaction(r, A_lj->c6, A_lj->c12,
              B_lj->c6, B_lj->c12,
              A_q, B_q,
              alpha_lj, alpha_crf,
              topo.energy_group_scaling()[energy_group_pair].first,
              topo.energy_group_scaling()[energy_group_pair].second,
              f1, f6, f12,
              e_lj, e_crf, de_lj, de_crf);

    } else {
      // no scaling
      lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
              B_lj->c6, B_lj->c12,
              A_q, B_q,
              alpha_lj, alpha_crf,
              f1, f6, f12,
              e_lj, e_crf, de_lj, de_crf);
    }

    //--------------------------------------------------
    // interactions have been calculated
    //--------------------------------------------------

    DEBUG(8, "\tdoscaling: calculated interaction state A:\n\t\tf: "
            << f1 << " f6: " << f6 << " f12: " << f12
            << "\n\t\te_lj: " << e_lj
            << " e_crf: " << e_crf
            << " de_lj: " << de_lj << " de_crf: " << de_crf
            << "\n\t\tr: " << sqrt(math::abs2(r)));

    // now combine everything
    f = (f1 + f6 + f12) * r;

    storage.force(i) += f;
    storage.force(j) -= f;

    DEBUG(8, "\tforces stored");

    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        storage.virial_tensor(a, b) +=
              r(a) * f(b);

    DEBUG(8, "\tatomic virial done");

  }// END OF SCALING ON ---
    //
  else {

    switch (t_interaction_spec::interaction_func) {
      case simulation::lj_crf_func:
      {
        math::Vec f;
        double f1 = 0.0, f6 = 0.0, f12 = 0.0;

        lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
                B_lj->c6, B_lj->c12,
                A_q, B_q,
                alpha_lj, alpha_crf,
                f1, f6, f12,
                e_lj, e_crf, de_lj, de_crf);
        //--------------------------------------------------
        // interactions have been calculated
        //--------------------------------------------------

        DEBUG(8, "\tnoscaling: calculated interaction state A:\n\t\tf: "
                << f1 << " f6: " << f6 << " f12: " << f12
                << "\n\t\te_lj: " << e_lj
                << " e_crf: " << e_crf
                << " de_lj: " << de_lj << " de_crf: " << de_crf
                << "\n\t\tr: " << sqrt(math::abs2(r)));

        // now combine everything
        f = (f1 + f6 + f12) * r;

        storage.force(i) += f;
        storage.force(j) -= f;

        DEBUG(8, "\tforces stored");

        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            storage.virial_tensor(a, b) +=
                  r(a) * f(b);

        DEBUG(8, "\tatomic virial done");

        //---------------------------------------------------------
        //                     ANITA
        // extended TI: calculate A_e_lj, B_e_lj, A_e_crf, B_e_crf,
        //              A_de_LJ, B_de_lj, A_de_crf, B_de_crf
        //---------------------------------------------------------

        // TODO: could add another parameter, to only calculate every x steps
        // if nr_lambdas > 1, we apply extended TI 
        if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
            DEBUG(8, "precalculate lj_crf_soft");
//        if ( sim.param().precalclam.nr_lambdas ) { 
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
            lj_crf_soft_interaction_ext(r, A_lj->c6, A_lj->c12,
                B_lj->c6, B_lj->c12, A_q, B_q, alpha_lj, alpha_crf,
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
        DEBUG(8, "\ndone with lj_crf_func ");

        break;
      }
      case simulation::cgrain_func:
      {
        DEBUG(7, "\tcgrain_func");
        math::Vec f;
        double f1 = 0.0, f6 = 0.0, f12 = 0.0;

        cgrain_soft_interaction(r, A_lj->c6, A_lj->c12,
                B_lj->c6, B_lj->c12,
                A_q, B_q,
                alpha_lj, alpha_crf,
                f1, f6, f12,
                e_lj, e_crf, de_lj, de_crf);
        //--------------------------------------------------
        // interactions have been calculated
        //--------------------------------------------------

        DEBUG(8, "\tcalculated interaction state A:\n\t\tf: "
                << f1 << " f6: " << f6 << " f12: " << f12
                << "\n\t\te_lj: " << e_lj
                << " e_crf: " << e_crf
                << " de_lj: " << de_lj << " de_crf: " << de_crf
                << "\n\t\tr: " << sqrt(math::abs2(r)));

        // now combine everything
        f = (f1 + f6 + f12) * r;

        storage.force(i) += f;
        storage.force(j) -= f;

        DEBUG(8, "\tforces stored");

        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            storage.virial_tensor(a, b) +=
                  r(a) * f(b);

        DEBUG(8, "\tatomic virial done");
        break;
      }
      case simulation::cggromos_func:
      {
        DEBUG(7, "\tcggromos_func");
        math::Vec f;
        double f1 = 0.0, f6 = 0.0, f12 = 0.0;

        // check if...
        if (topo.is_coarse_grained(i) && topo.is_coarse_grained(j)) { // CG-CG interaction
          lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
                  B_lj->c6, B_lj->c12,
                  A_q / cgrain_eps[0], B_q / cgrain_eps[0],
                  alpha_lj, alpha_crf,
                  f1, f6, f12,
                  e_lj, e_crf, de_lj, de_crf, 0);
        } else if (topo.is_coarse_grained(i) || topo.is_coarse_grained(j)) { // FG-CG interaction
          lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
                  B_lj->c6, B_lj->c12,
                  A_q / cgrain_eps[1], B_q / cgrain_eps[1],
                  alpha_lj, alpha_crf,
                  f1, f6, f12,
                  e_lj, e_crf, de_lj, de_crf, 1);
        } else { // FG-FG interaction
          lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
                  B_lj->c6, B_lj->c12,
                  A_q, B_q,
                  alpha_lj, alpha_crf,
                  f1, f6, f12,
                  e_lj, e_crf, de_lj, de_crf, 2);
        }

        //--------------------------------------------------
        // interactions have been calculated
        //--------------------------------------------------

        DEBUG(8, "\tcg: calculated interaction state A:\n\t\tf: "
                << f1 << " f6: " << f6 << " f12: " << f12
                << "\n\t\te_lj: " << e_lj
                << " e_crf: " << e_crf
                << " de_lj: " << de_lj << " de_crf: " << de_crf
                << "\n\t\tr: " << sqrt(math::abs2(r)));

        // now combine everything
        f = (f1 + f6 + f12) * r;

        storage.force(i) += f;
        storage.force(j) -= f;

        DEBUG(8, "\tforces stored");

        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            storage.virial_tensor(a, b) +=
                  r(a) * f(b);

        DEBUG(8, "\tatomic virial done");
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

        rp1 = r - conf.current().posV(j);
        rp2 = r + conf.current().posV(i);
        rpp = r + conf.current().posV(i) - conf.current().posV(j);

        pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
                A_lj->c6, A_lj->c12,
                B_lj->c6, B_lj->c12,
                A_qi, B_qi, A_qj, B_qj,
                topo.coscharge(i),
                topo.coscharge(j),
                alpha_lj, alpha_crf,
                f1, f6, f12,
                e_lj, e_crf, de_lj, de_crf);

        //--------------------------------------------------
        // interactions have been calculated
        //--------------------------------------------------

        // now combine everything
        f(0) = (f1[0] + f6 + f12) * r;
        f(1) = f1[1] * rp1;
        f(2) = f1[2] * rp2;
        f(3) = f1[3] * rpp;

        storage.force(i) += f(0) + f(1) + f(2) + f(3);
        storage.force(j) -= f(0) + f(1) + f(2) + f(3);

        DEBUG(7, "\tforces stored");

        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            storage.virial_tensor(a, b) += r(a)*(f(0)(b) +
                  f(1)(b) + f(2)(b) + f(3)(b));

        DEBUG(7, "\tatomic virial done");
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

        rp1 = rm - conf.current().posV(j);
        rp2 = rm + conf.current().posV(i);
        rpp = rm + conf.current().posV(i) - conf.current().posV(j);

        pol_off_lj_crf_soft_interaction(r, rm, rp1, rp2, rpp,
                A_lj->c6, A_lj->c12,
                B_lj->c6, B_lj->c12,
                A_qi, B_qi, A_qj, B_qj,
                topo.coscharge(i),
                topo.coscharge(j),
                alpha_lj, alpha_crf,
                f1, f6, f12,
                e_lj, e_crf, de_lj, de_crf);

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
        break;
      }
      default:
        io::messages.add("Nonbonded_Innerloop",
                "interaction function not implemented",
                io::message::critical);
    }
  }
DEBUG(8, "\tenergies perturbed lj_crf_innerloop;");
  // energy
  assert(storage.energies.lj_energy.size() >
          topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
          topo.atom_energy_group(j));

  storage.energies.lj_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_lj;

  storage.energies.crf_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_crf;

  DEBUG(8, "\tenergy group: i and j " << topo.atom_energy_group(i)
          << " " << topo.atom_energy_group(j)
          << " pert der index = " << energy_derivative_index);

  assert(storage.perturbed_energy_derivatives.
          lj_energy.size() > topo.atom_energy_group(i) &&
          storage.perturbed_energy_derivatives.
          lj_energy.size() > topo.atom_energy_group(j));

  assert(storage.perturbed_energy_derivatives.
          lj_energy[topo.atom_energy_group(i)].size() > topo.atom_energy_group(i) &&
          storage.perturbed_energy_derivatives.
          lj_energy[topo.atom_energy_group(i)].size() > topo.atom_energy_group(j));

  storage.perturbed_energy_derivatives.lj_energy
          [topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += de_lj;

  storage.perturbed_energy_derivatives.crf_energy
          [topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += de_crf;

  DEBUG(7, "\tperturbed lj_crf_innerloop " << i << " - " << j << " done!");
}

template<typename t_interaction_spec,
typename t_perturbation_details>
void interaction::Perturbed_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>
::perturbed_one_four_interaction_innerloop
(topology::Topology & topo,
        configuration::Configuration & conf,
        unsigned int i, unsigned int j,
        Periodicity_type const & periodicity,
       // ANITA 
        simulation::Simulation & sim //
) {
  DEBUG(7, "\tone four pair\t" << i << "\t" << j);

  math::Vec r;
  double e_lj = 0.0, e_crf = 0.0, de_lj = 0.0, de_crf = 0.0;

  periodicity.nearest_image(conf.current().pos(i),
          conf.current().pos(j), r);

  lj_parameter_struct const * A_lj = nullptr;
  lj_parameter_struct const * B_lj = nullptr;
  double A_q = 0.0, B_q = 0.0, A_qi = 0.0, A_qj = 0.0, B_qi = 0.0, B_qj = 0.0;
  double alpha_lj = 0, alpha_crf = 0;

  // const double l = topo.lambda();

  if (topo.is_perturbed(j) == true) {
    A_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
            topo.perturbed_solute().atoms()[j].A_IAC());
    B_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
            topo.perturbed_solute().atoms()[j].B_IAC());

    A_qi = topo.perturbed_solute().atoms()[i].A_charge();
    A_qj = topo.perturbed_solute().atoms()[j].A_charge();
    B_qi = topo.perturbed_solute().atoms()[i].B_charge();
    B_qj = topo.perturbed_solute().atoms()[j].B_charge();
    A_q = A_qi * A_qj;
    B_q = B_qi * B_qj;

    alpha_lj = (topo.perturbed_solute().atoms()[i].LJ_softcore() +
            topo.perturbed_solute().atoms()[j].LJ_softcore()) /
            2.0;
    alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
            topo.perturbed_solute().atoms()[j].CRF_softcore()) /
            2.0;

  } else {
    A_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
            topo.iac(j));
    B_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
            topo.iac(j));

    A_qi = topo.perturbed_solute().atoms()[i].A_charge();
    A_qj = topo.charge()(j);
    B_qi = topo.perturbed_solute().atoms()[i].B_charge();
    B_qj = topo.charge()(j);
    A_q = A_qi * A_qj;
    B_q = B_qi * B_qj;

    alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
    alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

  }

  DEBUG(7, "\tlj-parameter state A c6=" << A_lj->cs6
          << " c12=" << A_lj->cs12);
  DEBUG(7, "\tlj-parameter state B c6=" << B_lj->cs6
          << " c12=" << B_lj->cs12);
  DEBUG(7, "\tcharges state A i*j = " << A_q);
  DEBUG(7, "\tcharges state B i*j = " << B_q);

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

  if (t_perturbation_details::do_scaling) {
    // SCALING ON
    math::Vec f;
    double f1 = 0.0, f6 = 0.0, f12 = 0.0;
    if (t_interaction_spec::interaction_func == simulation::pol_lj_crf_func
            || t_interaction_spec::interaction_func == simulation::pol_lj_crf_func) {
      io::messages.add("Nonbonded Innerloop",
              "scaling not implemented for COS polarisation simulations!",
              io::message::critical);
      return;
    }
    // check whether we need to do scaling
    std::pair<int, int> energy_group_pair(topo.atom_energy_group(i),
            topo.atom_energy_group(j));
    
    // The coulomb scaling factor corresponds to the scaling of 
    // 1-4 interactions with the amber force field. It is not the 
    // same as the "scaled" interaction this if checks for.
    // note: if amber FF is used, both scaling are used.

    if (topo.energy_group_scaling().count(energy_group_pair)) {

      // YES, we do scale the interactions!
      lj_crf_scaled_interaction
              (r, A_lj->cs6, A_lj->cs12,
              B_lj->cs6, B_lj->cs12,
              A_q, B_q,
              alpha_lj, alpha_crf,
              topo.energy_group_scaling()[energy_group_pair].first,
              topo.energy_group_scaling()[energy_group_pair].second,
              f1, f6, f12,
              e_lj, e_crf, de_lj, de_crf,
              0, m_param->get_coulomb_scaling());
    } else {
      // no scaling
      lj_crf_soft_interaction
              (r, A_lj->cs6, A_lj->cs12,
              B_lj->cs6, B_lj->cs12,
              A_q, B_q,
              alpha_lj, alpha_crf,
              f1, f6, f12,
              e_lj, e_crf, de_lj, de_crf, 
              0, m_param->get_coulomb_scaling());
    }

    //--------------------------------------------------
    // interactions have been calculated
    //--------------------------------------------------

    DEBUG(7, "\tcalculated interaction state A:\n\t\tf: "
            << f1 << " " << f6 << " " << f12 << " e_lj: " << e_lj
            << " e_crf: " << e_crf
            << " de_lj: " << de_lj << " de_crf: " << de_crf);

    // now combine everything
    f = (f1 + f6 + f12) * r;

    conf.current().force(i) += f;
    conf.current().force(j) -= f;

    DEBUG(7, "\tforces stored");

    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        conf.current().virial_tensor(a, b) +=
              r(a) * f(b);

    DEBUG(7, "\tatomic virial done");

  }// END OF SCALING ON ---
    //
  else {
    switch (t_interaction_spec::interaction_func) {
      case simulation::lj_crf_func:
      {
        math::Vec f;
        double f1 = 0.0, f6 = 0.0, f12 = 0.0;

        lj_crf_soft_interaction
                (r, A_lj->cs6, A_lj->cs12,
                B_lj->cs6, B_lj->cs12,
                A_q, B_q,
                alpha_lj, alpha_crf,
                f1, f6, f12,
                e_lj, e_crf, de_lj, de_crf,
                0, m_param->get_coulomb_scaling());

        //---------------------------------------------------------
        //                     ANITA
        // extended TI: calculate A_e_lj, B_e_lj, A_e_crf, B_e_crf,
        //              A_de_LJ, B_de_lj, A_de_crf, B_de_crf
        //---------------------------------------------------------

        // TODO: could add another parameter, to only calculate every x steps
        // if nr_lambdas > 1, we apply extended TI
        if (sim.param().precalclam.nr_lambdas && ((sim.steps()  % sim.param().write.free_energy) == 0)){
 
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
            lj_crf_soft_interaction_ext(r, A_lj->cs6, A_lj->cs12,
                B_lj->cs6, B_lj->cs12, A_q, B_q, alpha_lj, alpha_crf,
                A_e_lj,  B_e_lj, A_e_crf, B_e_crf,
                A_de_lj, B_de_lj, A_de_crf, B_de_crf,
                lam, 0, m_param->get_coulomb_scaling());

            DEBUG(8, "ANITA: perturbed one_four: precalculated energies for lambda " << lam
                   << "\n now starting storage");
            DEBUG(8, "\n  A_e_lj " << A_e_lj << "\n  lambda index " << lam_index <<
                   "\n  conf.current().energies.A_lj_energy.size() " << conf.current().energies.A_lj_energy.size()
                   << "\n  energy group1 " << topo.atom_energy_group(i) << " energy group2 " 
                   << topo.atom_energy_group(j));
//            assert(storage.energies.A_lj_energy.size() > lam_index);
//            assert(storage.energies.A_lj_energy[lam_index].size() > topo.atom_energy_group(i));
//            assert(storage.energies.A_lj_energy[lam_index][topo.atom_energy_group(i)].size() 
//                     > topo.atom_energy_group(j));

            conf.current().energies.A_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_lj;
            conf.current().energies.B_lj_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_lj;

            conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_e_crf;
            conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_e_crf;

            conf.current().perturbed_energy_derivatives.A_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_lj;
            conf.current().perturbed_energy_derivatives.B_lj_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_lj;

            conf.current().perturbed_energy_derivatives.A_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += A_de_crf;
            conf.current().perturbed_energy_derivatives.B_crf_energy
                    [lam_index][topo.atom_energy_group(i)]
                    [topo.atom_energy_group(j)] += B_de_crf;
          }
        }
        //--------------------------------------------------
        // interactions have been calculated
        //--------------------------------------------------

        DEBUG(7, "\tcalculated interaction state A:\n\t\tf: "
                << f1 << " " << f6 << " " << f12 << " e_lj: " << e_lj
                << " e_crf: " << e_crf
                << " de_lj: " << de_lj << " de_crf: " << de_crf);

        // now combine everything
        f = (f1 + f6 + f12) * r;

        conf.current().force(i) += f;
        conf.current().force(j) -= f;

        DEBUG(7, "\tforces stored");

        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            conf.current().virial_tensor(a, b) += r(a) * f(b);

        DEBUG(7, "\tatomic virial done");

        break;
      }
      case simulation::cgrain_func:
      {
        io::messages.add("Nonbonded_Innerloop",
                "no perturbed 1,4 interactions for Martini coarse-grained simulations!",
                io::message::critical);
        break;
      }
      case simulation::cggromos_func:
      {
        // check if...
        if (topo.is_coarse_grained(i)) { // CG-CG or CG-FG interaction
          io::messages.add("Nonbonded_Innerloop",
                "no perturbed 1,4 interactions for Gromos coarse-grained simulations!",
                io::message::critical);
        } else { // FG-FG interaction
          math::Vec f;
          double f1 = 0.0, f6 = 0.0, f12 = 0.0;
          lj_crf_soft_interaction
                  (r, A_lj->cs6, A_lj->cs12,
                  B_lj->cs6, B_lj->cs12,
                  A_q, B_q,
                  alpha_lj, alpha_crf,
                  f1, f6, f12,
                  e_lj, e_crf, de_lj, de_crf, 2);

          DEBUG(7, "\tcalculated interaction state A:\n\t\tf: "
                  << f1 << " " << f6 << " " << f12 << " e_lj: " << e_lj
                  << " e_crf: " << e_crf
                  << " de_lj: " << de_lj << " de_crf: " << de_crf);

          // now combine everything
          f = (f1 + f6 + f12) * r;

          conf.current().force(i) += f;
          conf.current().force(j) -= f;

          DEBUG(7, "\tforces stored");

          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
              conf.current().virial_tensor(a, b) += r(a) * f(b);

          DEBUG(7, "\tatomic virial done");
        }
        break;
      }
      case simulation::pol_lj_crf_func:
      {
        math::Vec rp1, rp2, rpp;
        double f1[4];
        math::VArray f(4);
        f = 0.0;
        double f6 = 0.0, f12 = 0.0;

        rp1 = r - conf.current().posV(j);
        rp2 = r + conf.current().posV(i);
        rpp = r + conf.current().posV(i) - conf.current().posV(j);

        pol_lj_crf_soft_interaction(r, rp1, rp2, rpp,
                A_lj->cs6, A_lj->cs12,
                B_lj->cs6, B_lj->cs12,
                A_qi, B_qi, A_qj, B_qj,
                topo.coscharge(i),
                topo.coscharge(j),
                alpha_lj, alpha_crf,
                f1, f6, f12,
                e_lj, e_crf, de_lj, de_crf);

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
        break;
      }
      case simulation::pol_off_lj_crf_func:
      {
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

        rp1 = rm - conf.current().posV(j);
        rp2 = rm + conf.current().posV(i);
        rpp = rm + conf.current().posV(i) - conf.current().posV(j);

        pol_off_lj_crf_soft_interaction(r, rm, rp1, rp2, rpp,
                A_lj->cs6, A_lj->cs12,
                B_lj->cs6, B_lj->cs12,
                A_qi, B_qi, A_qj, B_qj,
                topo.coscharge(i),
                topo.coscharge(j),
                alpha_lj, alpha_crf,
                f1, f6, f12,
                e_lj, e_crf, de_lj, de_crf);

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
        break;
      }
      default:
        io::messages.add("Nonbonded_Innerloop",
                "interaction function not implemented",
                io::message::critical);
    }
  }


  // energy
  assert(conf.current().energies.lj_energy.size() >
          topo.atom_energy_group(i));
  assert(conf.current().energies.lj_energy.size() >
          topo.atom_energy_group(j));

  conf.current().energies.lj_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_lj;

  conf.current().energies.crf_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += e_crf;

  conf.current().perturbed_energy_derivatives.lj_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += de_lj;
  conf.current().perturbed_energy_derivatives.crf_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(j)] += de_crf;


}

template<typename t_interaction_spec,
typename t_perturbation_details>
inline void
interaction::Perturbed_Nonbonded_Innerloop<
t_interaction_spec, t_perturbation_details>
::perturbed_RF_excluded_interaction_innerloop
(topology::Topology & topo, configuration::Configuration & conf,
        std::map<unsigned int, topology::Perturbed_Atom>::const_iterator const & mit,
        Periodicity_type const & periodicity,
        // ANITA 
        simulation::Simulation & sim //
) {

  math::VArray &pos = conf.current().pos;
  math::VArray &force = conf.current().force;

  math::Vec r;
  double e_rf = 0.0, de_rf = 0.0;

  topology::excl_cont_t::value_type::const_iterator it, to;

  // self term has already been calculated for state A, 
  // correct for that and 
  // calculate it for this lambda
  // only a distance independent part
  r = 0.0;
  const int i = mit->second.sequence_number();
  const double q_i_a = mit->second.A_charge();
  const double q_i_b = mit->second.B_charge();
  double q_j_a = 0.0, q_j_b = 0.0;

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
      rf_soft_interaction(r, q_i_a*q_i_a, q_i_b * q_i_b,
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
        rf_soft_interaction(r, q_i_a * q_i_a / cgrain_eps[0],
                q_i_b * q_i_b / cgrain_eps[0],
                mit->second.CRF_softcore(),
                f_rf, e_rf, de_rf, true, 0);
      } else { // FG-FG interaction
        rf_soft_interaction(r, q_i_a*q_i_a, q_i_b * q_i_b,
                mit->second.CRF_softcore(),
                f_rf, e_rf, de_rf, true, 2);
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

  conf.current().energies.crf_energy[topo.atom_energy_group(i)]
          [topo.atom_energy_group(i)] += 0.5 * e_rf;
  conf.current().perturbed_energy_derivatives.crf_energy
          [topo.atom_energy_group(i)]
          [topo.atom_energy_group(i)] += 0.5 * de_rf;

  // now loop over the exclusions
  // those are fortunately not in the normal exclusions!
  it = mit->second.exclusion().begin();
  to = mit->second.exclusion().end();

  for (; it != to; ++it) {
    periodicity.nearest_image(pos(i),
            pos(*it), r);

    DEBUG(8, "r2 i(" << i << "-" << *it << ") " << abs2(r));

    double q_ij_a = 0.0;
    double q_ij_b = 0.0;

    double alpha_crf = 0;

    if (unsigned(*it) < topo.num_solute_atoms() && topo.is_perturbed(*it)) {
      // j perturbed
      q_j_a = topo.perturbed_solute().atoms()[*it].A_charge();
      q_ij_a = q_i_a * q_j_a;
      q_j_b = topo.perturbed_solute().atoms()[*it].B_charge();
      q_ij_b = q_i_b * q_j_b;

      alpha_crf = (mit->second.CRF_softcore() +
              topo.perturbed_solute().
              atoms()[*it].CRF_softcore()) * 0.5;
    } else {
      // only i perturbed
      q_j_a = topo.charge()(*it);
      q_ij_a = q_i_a * q_j_a;
      q_j_b = topo.charge()(*it);
      q_ij_b = q_i_b * q_j_b;

      alpha_crf = mit->second.CRF_softcore();

    }
    DEBUG(8, "q_i_a " << q_i_a << " q_i_b " << q_i_b
            << " q_ij_a " << q_ij_a << " q_ij_b " << q_ij_b
            << " A_l " << m_A_crf_lambda << " B_l " << m_B_crf_lambda);

    int n1 = topo.atom_energy_group(i);
    int n2 = topo.atom_energy_group(*it);

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

    switch (t_interaction_spec::interaction_func) {
      case simulation::lj_crf_func:
      {
        math::Vec f_rf;

        rf_soft_interaction(r, q_ij_a, q_ij_b, alpha_crf, f_rf, e_rf, de_rf);

        DEBUG(8, "alpha_crf : " << alpha_crf);
        DEBUG(7, "excluded atoms " << i << " & " << *it << ": " << e_rf);
        DEBUG(7, "\tde_rf: " << de_rf);

        // and add everything to the correct arrays
        conf.current().energies.crf_energy
                [topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_rf;
        conf.current().perturbed_energy_derivatives.crf_energy
                [topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += de_rf;
        force(i) += f_rf;
        force(*it) -= f_rf;

        // ANITA
        if (sim.param().precalclam.nr_lambdas && ((sim.steps() % sim.param().write.free_energy) == 0)){
//        if ( sim.param().precalclam.nr_lambdas ) {
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

            conf.current().energies.A_crf_energy[lam_index][topo.atom_energy_group(i)]
              [topo.atom_energy_group(*it)] += A_e_rf;
            conf.current().energies.B_crf_energy[lam_index][topo.atom_energy_group(i)]
              [topo.atom_energy_group(*it)] += B_e_rf;
            conf.current().perturbed_energy_derivatives.A_crf_energy[lam_index]
              [topo.atom_energy_group(i)]
              [topo.atom_energy_group(*it)] += A_de_rf;
            conf.current().perturbed_energy_derivatives.B_crf_energy[lam_index]
              [topo.atom_energy_group(i)]
              [topo.atom_energy_group(*it)] += B_de_rf;
 
          }
        } // ANITA

        // if (t_interaction_spec::do_virial != math::no_virial){
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            conf.current().virial_tensor(a, b) +=
                  r(a) * f_rf(b);

        DEBUG(7, "\tatomic virial done");
        // }
        break;
      }
      case simulation::cggromos_func:
      {
        math::Vec f_rf;
        // check if...
        if (topo.is_coarse_grained(i) && topo.is_coarse_grained(*it)) { // CG-CG interaction
          rf_soft_interaction(r, q_ij_a / cgrain_eps[0], q_ij_b / cgrain_eps[0],
                  alpha_crf, f_rf, e_rf, de_rf, 0);
        } else if (topo.is_coarse_grained(i) || topo.is_coarse_grained(*it)) { // FG-CG interaction
          rf_soft_interaction(r, q_ij_a / cgrain_eps[1], q_ij_b / cgrain_eps[1],
                  alpha_crf, f_rf, e_rf, de_rf, 1);
        } else { // FG-FG interaction
          rf_soft_interaction(r, q_ij_a, q_ij_b, alpha_crf, f_rf, e_rf, de_rf, 2);
        }

        DEBUG(8, "alpha_crf : " << alpha_crf);
        DEBUG(7, "excluded atoms " << i << " & " << *it << ": " << e_rf);
        DEBUG(7, "\tde_rf: " << de_rf);

        // and add everything to the correct arrays
        conf.current().energies.crf_energy
                [topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_rf;
        conf.current().perturbed_energy_derivatives.crf_energy
                [topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += de_rf;
        force(i) += f_rf;
        force(*it) -= f_rf;

        // if (t_interaction_spec::do_virial != math::no_virial){
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            conf.current().virial_tensor(a, b) +=
                  r(a) * f_rf(b);

        DEBUG(7, "\tatomic virial done");
        // }
        break;
      }
      case simulation::pol_lj_crf_func:
      {
        math::Vec rp1, rp2, rpp;
        double f_rf[4];
        rp1 = r - conf.current().posV(*it);
        rp2 = r + conf.current().posV(i);
        rpp = r + conf.current().posV(i) - conf.current().posV(*it);

        pol_rf_soft_interaction(r, rp1, rp2, rpp,
                q_i_a, q_j_a, q_i_b, q_j_b,
                topo.coscharge(i), topo.coscharge(*it),
                alpha_crf, f_rf, e_rf, de_rf);

        // and add everything to the correct arrays
        conf.current().energies.crf_energy
                [topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_rf;
        conf.current().perturbed_energy_derivatives.crf_energy
                [topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += de_rf;
        force(i) += f_rf[0] * r + f_rf[1] * rp1 + f_rf[2] * rp2 + f_rf[3] * rpp;
        force(*it) -= f_rf[0] * r + f_rf[1] * rp1 + f_rf[2] * rp2 + f_rf[3] * rpp;

        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            conf.current().virial_tensor(a, b) +=
                  r(a)*(f_rf[0] * r(b) + f_rf[1] * rp1(b) +
                  f_rf[2] * rp2(b) + f_rf[3] * rpp(b));

        break;
      }
      case simulation::pol_off_lj_crf_func:
      {
        math::Vec rm = r;
        if (topo.gamma(i)!=0.0) {
          math::Vec rij, rik;
          periodicity.nearest_image(conf.current().pos(i),
                  conf.current().pos(topo.gamma_j(i)), rij);
          periodicity.nearest_image(conf.current().pos(i),
                  conf.current().pos(topo.gamma_k(i)), rik);
          rm -= topo.gamma(i)*(rij + rik) / 2;
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
        rp2 = rm + conf.current().posV(i);
        rpp = rm + conf.current().posV(i) - conf.current().posV(*it);

        pol_rf_soft_interaction(rm, rp1, rp2, rpp,
                q_i_a, q_j_a, q_i_b, q_j_b,
                topo.coscharge(i), topo.coscharge(*it),
                alpha_crf, f_rf, e_rf, de_rf);

        // and add everything to the correct arrays
        conf.current().energies.crf_energy
                [topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += e_rf;
        conf.current().perturbed_energy_derivatives.crf_energy
                [topo.atom_energy_group(i)]
                [topo.atom_energy_group(*it)] += de_rf;
        for (int a = 0; a < 3; ++a) {
          const double term = f_rf[0] * rm(a) + f_rf[1] * rp1(a)
                  + f_rf[2] * rp2(a) + f_rf[3] * rpp(a);

          force(i)(a) += (1 - topo.gamma(i)) * term;
          force(*it)(a) -= (1 - topo.gamma(*it)) * term;

          force(topo.gamma_j(i))(a) += topo.gamma(i) / 2 * term;
          force(topo.gamma_j(*it))(a) -= topo.gamma(*it) / 2 * term;

          force(topo.gamma_k(i))(a) += topo.gamma(i) / 2 * term;
          force(topo.gamma_k(*it))(a) -= topo.gamma(*it) / 2 * term;

          for (int b = 0; b < 3; ++b)
            conf.current().virial_tensor(b, a) += r(b) * term;
        }
        break;
      }
      default:
        io::messages.add("Nonbonded_Innerloop",
                "interaction function not implemented",
                io::message::critical);
    }
  }
}

template<typename t_interaction_spec, typename t_perturbation_details >
inline void interaction::Perturbed_Nonbonded_Innerloop<
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
interaction::Perturbed_Nonbonded_Innerloop<
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



