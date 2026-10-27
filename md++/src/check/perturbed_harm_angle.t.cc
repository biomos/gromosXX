/*
 * This file is part of GROMOS.
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 *
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file perturbed_harm_angle.t.cc
 * harmonic angle perturbation checks.
 */

#include "../stdheader.h"
#include "../algorithm/algorithm.h"
#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../interaction/interaction.h"
#include "../interaction/bonded/harm_angle_interaction.h"
#include "../interaction/bonded/perturbed_harm_angle_interaction.h"
#include "../interaction/bonded/create_bonded.h"
#include "../interaction/forcefield/forcefield.h"
#include "../interaction/forcefield/create_forcefield.h"
#include "../io/topology/in_topology.h"
#include "../io/topology/in_perturbation.h"
#include "../io/parameter/check_parameter.h"
#include "../util/create_simulation.h"
#include "check.h"

#ifdef XXMPI
#include <mpi.h>
#endif

namespace {
  int failures = 0;
  double max_force_error = 0.0, max_lambda_error = 0.0;

  void require(bool condition, std::string const &name)
  {
    if (!condition){
      std::cerr << "Failed: " << name << '\n';
      ++failures;
    }
  }

  void close(double actual, double expected, std::string const &name,
             double tolerance = 1.e-9)
  {
    require(std::isfinite(actual) && std::isfinite(expected) &&
            fabs(actual - expected) <= tolerance * (1 + fabs(expected)), name);
  }

  void reset(configuration::Configuration &conf)
  {
    conf.current().force = 0.0;
    conf.current().virial_tensor = 0.0;
    conf.current().energies.zero();
    conf.current().perturbed_energy_derivatives.zero();
    for(unsigned int g=0; g<conf.special().gamd.total_force.size(); ++g){
      conf.special().gamd.total_force[g] = 0.0;
      conf.special().gamd.virial_tensor[g] = 0.0;
    }
  }

  double evaluate(interaction::Interaction &term, util::simulation_struct &s)
  {
    reset(s.conf);
    require(term.calculate_interactions(s.topo, s.conf, s.sim) == 0, term.name + " calculation");
    double energy = 0.0;
    for(unsigned int g=0; g<s.conf.current().energies.angle_energy.size(); ++g)
      energy += s.conf.current().energies.angle_energy[g];
    return energy;
  }

  void compare_forces(math::VArray const &actual, math::VArray const &expected,
                      std::string const &name, double scale = 1.0)
  {
    for(unsigned int i=0; i<actual.size(); ++i)
      for(int d=0; d<3; ++d)
        close(actual(i)(d), scale * expected(i)(d), name);
  }

  void compare_virial(math::Matrix const &actual, math::Matrix const &expected,
                      std::string const &name, double scale = 1.0)
  {
    for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
        close(actual(a,b), scale * expected(a,b), name);
  }

  void check_forces(interaction::Interaction &term, util::simulation_struct &s)
  {
    evaluate(term, s);
    const math::VArray force = s.conf.current().force;
    const double h = 1.e-6;
    for(unsigned int i=0; i<3; ++i){
      for(int d=0; d<3; ++d){
        const double original = s.conf.current().pos(i)(d);
        s.conf.current().pos(i)(d) = original + h;
        const double plus = evaluate(term, s);
        s.conf.current().pos(i)(d) = original - h;
        const double minus = evaluate(term, s);
        s.conf.current().pos(i)(d) = original;
        const double numerical = (minus - plus) / (2 * h);
        max_force_error = std::max(max_force_error, fabs(force(i)(d) - numerical));
        close(force(i)(d), numerical, term.name + " finite-difference force", 2.e-6);
      }
    }
  }

  void check_lambda(interaction::Interaction &term, util::simulation_struct &s, double lambda)
  {
    const double h = 1.e-6;
    s.topo.lambda(lambda);
    s.topo.update_for_lambda();
    evaluate(term, s);
    const double derivative = s.conf.current().perturbed_energy_derivatives.angle_energy[1];
    s.topo.lambda(lambda + h);
    s.topo.update_for_lambda();
    const double plus = evaluate(term, s);
    s.topo.lambda(lambda - h);
    s.topo.update_for_lambda();
    const double minus = evaluate(term, s);
    const double numerical = (plus - minus) / (2 * h);
    max_lambda_error = std::max(max_lambda_error, fabs(derivative - numerical));
    close(derivative, numerical, "finite-difference lambda derivative", 2.e-7);
    s.topo.lambda(lambda);
    s.topo.update_for_lambda();
  }

  int load(util::simulation_struct &s, io::In_Topology &reader)
  {
    std::string topo, conf, input;
    GETFILEPATH(topo, "aladip.topo", "src/check/data/");
    GETFILEPATH(conf, "aladip.conf", "src/check/data/");
    GETFILEPATH(input, "aladip.in", "src/check/data/");
    reader.quiet = true;
    const int error = util::create_simulation(topo, "", conf, input, s, reader);
    require(error == 0 && !io::messages.contains(io::message::error), "load alanine dipeptide");
    io::messages.clear();
    return error;
  }

  void check_numerics()
  {
    util::simulation_struct s;
    io::In_Topology reader;
    if (load(s, reader)) return;
    topology::Topology &topo = s.topo;
    configuration::Configuration &conf = s.conf;
    simulation::Simulation &sim = s.sim;
    interaction::Harm_Angle_Interaction harmonic;
    interaction::Perturbed_Harm_Angle_Interaction perturbed;

    topo.solute().angles().clear();
    topo.solute().angles().push_back(topology::three_body_term_struct(0, 1, 2, 0));
    topo.perturbed_solute().angles().push_back(topology::perturbed_three_body_term_struct(0, 1, 2, 0, 1));
    topo.angle_types_harm().clear();
    topo.angle_types_harm().push_back(interaction::angle_type_struct(300.0, 1.6));
    topo.angle_types_harm().push_back(interaction::angle_type_struct(550.0, 2.2));
    topo.atom_energy_group()[0] = 1;
    conf.current().pos(0) = math::Vec(0.13, -0.04, 0.06);
    conf.current().pos(1) = math::Vec(0.0, 0.0, 0.0);
    conf.current().pos(2) = math::Vec(-0.05, 0.11, 0.03);
    conf.boundary_type = math::vacuum;
    sim.param().pcouple.virial = math::atomic_virial;
    const math::Vec u = conf.current().pos(0), v = conf.current().pos(2);
    const double theta = atan2(sqrt(abs2(cross(u,v))), dot(u,v));
    const double ordinary = evaluate(harmonic, s);
    close(ordinary, 150.0 * (theta - 1.6) * (theta - 1.6), "ordinary harmonic energy");
    check_forces(harmonic, s);

    const double lambdas[] = {0.0, 0.2, 0.5, 0.8, 1.0};
    for(unsigned int l=0; l<5; ++l){
      const double lambda = lambdas[l];
      topo.lambda(lambda);
      topo.update_for_lambda();
      const double energy = evaluate(perturbed, s);
      const double K = 300.0 + 250.0 * lambda;
      const double diff = theta - (1.6 + 0.6 * lambda);
      close(energy, 0.5 * K * diff * diff, "interpolated harmonic energy");
      close(conf.current().perturbed_energy_derivatives.angle_energy[1],
            125.0 * diff * diff - 0.6 * K * diff, "analytical lambda derivative");
      close(conf.current().energies.angle_energy[0], 0.0, "energy group selection");
      const double derivative = conf.current().perturbed_energy_derivatives.angle_energy[1];
      const math::VArray forces = conf.current().force;
      const math::Matrix virial = conf.current().virial_tensor;
      math::Matrix expected_virial(0.0);
      for(int a=0; a<3; ++a)
        for(int b=0; b<3; ++b)
          expected_virial(a,b) = u(a) * forces(0)(b) + v(a) * forces(2)(b);
      compare_virial(virial, expected_virial, "atomic virial");
      for(int d=0; d<3; ++d)
        close(forces(0)(d) + forces(1)(d) + forces(2)(d), 0.0, "total force");
      require(perturbed.calculate_interactions(topo, conf, sim) == 0, "additive calculation");
      close(conf.current().energies.angle_energy[1], 2 * energy, "energy accumulation");
      close(conf.current().perturbed_energy_derivatives.angle_energy[1], 2 * derivative,
            "lambda derivative accumulation");
      compare_forces(conf.current().force, forces, "force accumulation", 2.0);
      compare_virial(conf.current().virial_tensor, virial, "virial accumulation", 2.0);
      if (l == 2){
        const math::VArray positions = conf.current().pos;
        const double h = 1.e-6;
        for(int a=0; a<3; ++a){
          for(int b=0; b<3; ++b){
            for(unsigned int i=0; i<3; ++i)
              conf.current().pos(i)(b) = positions(i)(b) + h * positions(i)(a);
            const double plus = evaluate(perturbed, s);
            conf.current().pos = positions;
            for(unsigned int i=0; i<3; ++i)
              conf.current().pos(i)(b) = positions(i)(b) - h * positions(i)(a);
            const double minus = evaluate(perturbed, s);
            conf.current().pos = positions;
            close(virial(a,b), (minus - plus) / (2 * h), "finite-difference virial", 2.e-6);
          }
        }
      }
      if (l == 0 || l == 4){
        topo.solute().angles()[0].type = l == 0 ? 0 : 1;
        close(evaluate(harmonic, s), energy, "endpoint energy");
        compare_forces(conf.current().force, forces, "endpoint forces");
        compare_virial(conf.current().virial_tensor, virial, "endpoint virial");
      }
      check_forces(perturbed, s);
      check_lambda(perturbed, s, lambda);
    }
    topo.perturbed_solute().angles()[0].B_type = 0;
    topo.solute().angles()[0].type = 0;
    close(evaluate(perturbed, s), ordinary, "identical-state energy");
    close(conf.current().perturbed_energy_derivatives.angle_energy[1], 0.0, "identical-state derivative");
    topo.perturbed_solute().angles()[0].B_type = 1;

    sim.param().lambdas.individual_lambdas = true;
    sim.param().lambdas.c[simulation::angle_lambda][1][1] = 1.0;
    sim.param().lambdas.d[simulation::angle_lambda][1][1] = 0.0;
    topo.init(sim);
    topo.lambda(0.37);
    topo.update_for_lambda();
    close(topo.individual_lambda(simulation::angle_lambda)[1][1], 0.37 * 0.37, "individual lambda schedule");
    check_lambda(perturbed, s, 0.37);

    sim.param().precalclam.nr_lambdas = 5;
    sim.param().precalclam.min_lam = 0.1;
    sim.param().precalclam.max_lam = 0.9;
    sim.param().write.free_energy = 2;
    conf.current().energies.resize(2, 2, 5);
    conf.current().perturbed_energy_derivatives.resize(2, 2, 5);
    sim.steps() = 0;
    evaluate(perturbed, s);
    for(unsigned int l=0; l<5; ++l){
      const double lam = 0.1 + l * 0.2;
      const double K = 300 + 250 * lam, diff = theta - 1.6 - 0.6 * lam;
      close(conf.current().energies.AB_angle[l], 0.5 * K * diff * diff, "extended-TI energy");
      close(conf.current().perturbed_energy_derivatives.AB_angle[l],
            125 * diff * diff - 0.6 * K * diff, "extended-TI derivative");
    }
    sim.steps() = 1;
    evaluate(perturbed, s);
    for(unsigned int l=0; l<5; ++l){
      close(conf.current().energies.AB_angle[l], 0.0, "extended-TI output stride");
      close(conf.current().perturbed_energy_derivatives.AB_angle[l], 0.0, "extended-TI derivative stride");
    }
    sim.steps() = 0;
    sim.param().precalclam.nr_lambdas = 1;
    evaluate(perturbed, s);
    close(conf.current().energies.AB_angle[0], 0.5 * 325 * pow(theta - 1.66, 2), "single lambda sample");
    sim.param().precalclam.nr_lambdas = 0;

    const double energy = evaluate(perturbed, s);
    const math::VArray force = conf.current().force;
    const math::Matrix virial = conf.current().virial_tensor;
    const math::boundary_enum boundaries[] = {math::rectangular, math::triclinic, math::truncoct};
    const math::virial_enum virials[] = {math::no_virial, math::atomic_virial, math::molecular_virial};
    conf.current().box = math::Box(math::Vec(3,0,0), math::Vec(0,3,0), math::Vec(0,0,3));
    for(unsigned int b=0; b<3; ++b){
      conf.boundary_type = boundaries[b];
      conf.current().pos(0) += conf.current().box(0);
      for(unsigned int w=0; w<3; ++w){
        sim.param().pcouple.virial = virials[w];
        close(evaluate(perturbed, s), energy, "periodic-image energy");
        compare_forces(conf.current().force, force, "periodic-image forces");
        compare_virial(conf.current().virial_tensor, virial, "periodic-image virial");
      }
      conf.current().pos(0) -= conf.current().box(0);
    }
    conf.boundary_type = math::vacuum;
    sim.param().gamd.gamd = 1;
    topo.gamd_accel_group().assign(topo.num_atoms(), 1);
    const std::vector<unsigned int> key = {1, 1};
    topo.get_gamd_interaction_pairs()[key] = 1;
    conf.special().gamd.total_force.resize(2);
    for(unsigned int g=0; g<2; ++g)
      conf.special().gamd.total_force[g].resize(topo.num_atoms());
    conf.special().gamd.virial_tensor.resize(2);
    conf.current().energies.gamd_potential_total.resize(2);
    close(evaluate(perturbed, s), energy, "GaMD physical energy");
    close(conf.current().energies.gamd_potential_total[0], 0.0, "GaMD group selection");
    close(conf.current().energies.gamd_potential_total[1], energy, "GaMD energy");
    compare_forces(conf.special().gamd.total_force[1], force, "GaMD forces");
    compare_virial(conf.special().gamd.virial_tensor[1], virial, "GaMD virial");
    require(perturbed.calculate_interactions(topo, conf, sim) == 0, "GaMD accumulation calculation");
    close(conf.current().energies.gamd_potential_total[1], 2 * energy, "GaMD energy accumulation");
    compare_forces(conf.special().gamd.total_force[1], force, "GaMD force accumulation", 2.0);
    compare_virial(conf.special().gamd.virial_tensor[1], virial, "GaMD virial accumulation", 2.0);
    sim.param().gamd.gamd = 0;

    topo.angle_types_harm()[0].cos0 = math::Pi;
    topo.angle_types_harm()[1].cos0 = math::Pi;
    conf.current().pos(0) = math::Vec(1,0,0);
    conf.current().pos(2) = math::Vec(-1,0,0);
    close(evaluate(perturbed, s), 0.0, "linear equilibrium energy");
    close(conf.current().perturbed_energy_derivatives.angle_energy[1], 0.0, "linear equilibrium derivative");
    for(unsigned int i=0; i<3; ++i)
      for(int d=0; d<3; ++d)
        close(conf.current().force(i)(d), 0.0, "linear equilibrium forces");
    const double offsets[] = {1.e-3, 1.e-6, 1.e-9};
    for(unsigned int i=0; i<3; ++i){
      conf.current().pos(2) = math::Vec(-cos(offsets[i]), sin(offsets[i]), 0);
      const double near_energy = evaluate(perturbed, s);
      const math::VArray near_force = conf.current().force;
      const math::Matrix near_virial = conf.current().virial_tensor;
      const double lam = topo.individual_lambda(simulation::angle_lambda)[1][1];
      const double K = 300 + 250 * lam;
      topo.angle_types_harm()[0].K = K;
      topo.solute().angles()[0].type = 0;
      close(evaluate(harmonic, s), near_energy, "near-linear ordinary energy");
      compare_forces(conf.current().force, near_force, "near-linear ordinary forces");
      compare_virial(conf.current().virial_tensor, near_virial, "near-linear ordinary virial");
      topo.angle_types_harm()[0].K = 300;
      if (i < 2) check_forces(perturbed, s);
      check_lambda(perturbed, s, 0.37);
    }
    conf.current().pos(2) = math::Vec(-1,0,0);
    topo.angle_types_harm()[0].cos0 = topo.angle_types_harm()[1].cos0 = 2.0;
    reset(conf);
    require(perturbed.calculate_interactions(topo, conf, sim) != 0 &&
            io::messages.contains(io::message::critical), "singular nonequilibrium diagnostic");
    io::messages.clear();
    conf.current().pos(0)(0) = std::numeric_limits<double>::quiet_NaN();
    reset(conf);
    require(perturbed.calculate_interactions(topo, conf, sim) != 0 &&
            io::messages.contains(io::message::critical), "nonfinite-geometry diagnostic");
    io::messages.clear();
    conf.current().pos(0) = math::Vec(1,0,0);
    conf.current().pos(2) = math::Vec(1,0,0);
    reset(conf);
    require(perturbed.calculate_interactions(topo, conf, sim) != 0 &&
            io::messages.contains(io::message::critical), "zero-angle diagnostic");
    io::messages.clear();
    conf.current().pos(0) = conf.current().pos(1);
    reset(conf);
    require(perturbed.calculate_interactions(topo, conf, sim) != 0 &&
            io::messages.contains(io::message::critical), "coincident-atoms diagnostic");
    io::messages.clear();
  }

  void check_linear_perturbation()
  {
    util::simulation_struct s;
    io::In_Topology reader;
    if (load(s, reader)) return;
    topology::Topology &topo = s.topo;
    configuration::Configuration &conf = s.conf;
    simulation::Simulation &sim = s.sim;
    interaction::Perturbed_Harm_Angle_Interaction perturbed;
    topo.perturbed_solute().angles().push_back(topology::perturbed_three_body_term_struct(0,1,2,0,1));
    topo.angle_types_harm().clear();
    topo.angle_types_harm().push_back(interaction::angle_type_struct(300.0, math::Pi));
    topo.angle_types_harm().push_back(interaction::angle_type_struct(550.0, 170.0 * math::Pi / 180.0));
    topo.atom_energy_group()[0] = 1;
    conf.current().pos(0) = math::Vec(1,0,0);
    conf.current().pos(1) = math::Vec(0,0,0);
    conf.boundary_type = math::vacuum;
    conf.current().energies.resize(2, 2, 5);
    conf.current().perturbed_energy_derivatives.resize(2, 2, 5);
    const double theta_diff = topo.angle_types_harm()[1].cos0 - math::Pi;
    const double eps = std::numeric_limits<double>::epsilon();
    const double offsets[] = {0.0, 0.5e-12, 2.e-12, 1.e-9, 1.e-8, 2.e-8, 1.e-6, 1.e-4};
    const double lambdas[] = {0.0, 1.e-13, 1.e-12, 2.e-12, 4.e-12, 6.e-12, 1.e-8, 1.e-4};
    double limit_error = 0.0, regular_error = 0.0, grid_error = 0.0;
    unsigned int limit_count = 0, regular_count = 0, rejected_count = 0;
    for(unsigned int schedule=0; schedule<2; ++schedule){
      sim.param().lambdas.individual_lambdas = true;
      sim.param().lambdas.c[simulation::angle_lambda][1][1] = schedule ? 0.5 : 0.0;
      sim.param().lambdas.d[simulation::angle_lambda][1][1] = schedule ? 0.5 : 1.0;
      topo.init(sim);
      for(unsigned int o=0; o<8; ++o){
        conf.current().pos(2) = math::Vec(-cos(offsets[o]), sin(offsets[o]), 0);
        const double cost = conf.current().pos(2)(0) / sqrt(abs2(conf.current().pos(2)));
        const double theta = acos(std::max(-1.0, std::min(1.0, cost)));
        const bool limit = sin(theta) < math::epsilon;
        for(unsigned int l=0; l<8; ++l){
          const double lambda = lambdas[l];
          topo.lambda(lambda);
          topo.update_for_lambda();
          const double lam = topo.individual_lambda(simulation::angle_lambda)[1][1];
          const double prime = topo.individual_lambda_derivative(simulation::angle_lambda)[1][1];
          const double K = (1 - lam) * 300.0 + lam * 550.0;
          const double theta0 = (1 - lam) * math::Pi + lam * topo.angle_types_harm()[1].cos0;
          sim.param().precalclam.nr_lambdas = 0;
          if (limit && fabs(theta0 - math::Pi) > math::epsilon){
            reset(conf);
            require(perturbed.calculate_interactions(topo, conf, sim) != 0 &&
                    io::messages.contains(io::message::critical), "changing linear equilibrium diagnostic");
            io::messages.clear();
            ++rejected_count;
            continue;
          }
          const double diff = theta - theta0;
          const double expected_energy = 0.5 * K * diff * diff;
          const double expected_derivative = prime * (125.0 * diff * diff - K * diff * theta_diff);
          const double energy = evaluate(perturbed, s);
          const double derivative = conf.current().perturbed_energy_derivatives.angle_energy[1];
          close(energy, expected_energy, "changing linear equilibrium energy",
                32 * eps * std::max(fabs(expected_energy), 1.e-300));
          close(derivative, expected_derivative, "changing linear equilibrium derivative",
                32 * eps * std::max(fabs(expected_derivative), 1.e-300));
          const double h = limit ? 5.e-13 : 1.e-7;
          topo.lambda(lambda + h);
          topo.update_for_lambda();
          const double plus = evaluate(perturbed, s);
          topo.lambda(lambda < h ? lambda + 2 * h : lambda - h);
          topo.update_for_lambda();
          const double other = evaluate(perturbed, s);
          const double numerical = lambda < h ? (-3 * energy + 4 * plus - other) / (2 * h) :
            (plus - other) / (2 * h);
          const double error = fabs(derivative - numerical);
          const double tolerance = 8 * eps * K * (fabs(diff) / h + fabs(theta_diff)) +
            (limit ? 0.0 : 2.e-12);
          require(std::isfinite(numerical) && error <= tolerance, "near-linear finite-difference derivative");
          if (limit){
            limit_error = std::max(limit_error, error);
            ++limit_count;
          }
          else{
            regular_error = std::max(regular_error, error);
            ++regular_count;
          }
          topo.lambda(lambda);
          topo.update_for_lambda();
          sim.param().precalclam.nr_lambdas = 5;
          sim.param().precalclam.min_lam = lam;
          sim.param().precalclam.max_lam = 1.0;
          sim.param().write.free_energy = 2;
          sim.steps() = 0;
          evaluate(perturbed, s);
          const std::vector<double> &e = conf.current().energies.AB_angle;
          const std::vector<double> &d = conf.current().perturbed_energy_derivatives.AB_angle;
          close(e[0], energy, "near-linear PRECALCLAM energy consistency",
                32 * eps * std::max(fabs(e[0]), 1.e-300));
          close(prime * d[0], derivative, "near-linear PRECALCLAM chain rule",
                32 * eps * std::max(fabs(prime * d[0]), 1.e-300));
          const double step = (1.0 - lam) / 4;
          const double finite[] = {(-25 * e[0] + 48 * e[1] - 36 * e[2] + 16 * e[3] - 3 * e[4]) / (12 * step),
                                  (e[0] - 8 * e[1] + 8 * e[3] - e[4]) / (12 * step),
                                  (25 * e[4] - 48 * e[3] + 36 * e[2] - 16 * e[1] + 3 * e[0]) / (12 * step)};
          for(unsigned int g=0; g<3; ++g){
            close(d[2*g], finite[g], "near-linear PRECALCLAM finite difference", 2.e-12);
            grid_error = std::max(grid_error, fabs(d[2*g] - finite[g]));
          }
        }
      }
    }
    sim.param().write.free_energy = 0;
    for(unsigned int step=0; step<2; ++step){
      sim.steps() = step;
      evaluate(perturbed, s);
      for(unsigned int g=0; g<5; ++g){
        close(conf.current().energies.AB_angle[g], 0.0, "PRECALCLAM zero output interval");
        close(conf.current().perturbed_energy_derivatives.AB_angle[g], 0.0,
              "PRECALCLAM zero derivative output interval");
      }
    }
    io::simple_crosschecks(sim);
    std::ostringstream messages;
    io::messages.display(messages);
    require(messages.str().find("PRECALCLAM requires NTWE=NTWG > 0") != std::string::npos,
            "PRECALCLAM input validation rejects zero interval");
    io::messages.clear();
    require(limit_count && regular_count && rejected_count, "both linear-angle calculation branches exercised");
    std::cout << std::scientific << std::setprecision(10)
              << "Changing-equilibrium limit cases: " << limit_count << ", regular cases: " << regular_count
              << ", singular cases rejected: " << rejected_count << '\n'
              << "Maximum limit-branch lambda error: " << limit_error << '\n'
              << "Maximum near-linear lambda error: " << regular_error << '\n'
              << "Maximum near-linear PRECALCLAM derivative error: " << grid_error << '\n';
  }

  void check_dispatch()
  {
    util::simulation_struct s;
    io::In_Topology reader;
    if (load(s, reader)) return;
    s.sim.param().force.angle = 2;
    std::istringstream input("TITLE\nNonbonded-only perturbation\nEND\nPERTATOMPARAM\n1\n"
                            "1 1 CB 14 15.035 0.0 11 15.035 1.0 0.3 0.6\nEND\n");
    io::In_Perturbation perturbation(input);
    perturbation.quiet = true;
    perturbation.read(s.topo, s.sim.param());
    s.topo.update_for_lambda();
    require(!s.topo.perturbed_solute().atoms().empty() && s.topo.perturbed_solute().angles().empty(),
            "nonbonded-only perturbation topology");
    interaction::Harm_Angle_Interaction harmonic;
    s.sim.param().perturbation.perturbation = false;
    const double reference = evaluate(harmonic, s);
    s.sim.param().perturbation.perturbation = true;
    const math::VArray force = s.conf.current().force;
    const math::Matrix virial = s.conf.current().virial_tensor;
    interaction::Forcefield ff;
    require(interaction::create_g96_forcefield(ff, s.topo, s.sim, reader, std::cout, true) == 0 &&
            !io::messages.contains(io::message::error), "harmonic/nonbonded perturbation dispatch");
    interaction::Interaction *angle = ff.interaction("HarmAngle");
    interaction::Interaction *perturbed = ff.interaction("PerturbedHarmAngle");
    require(angle && perturbed, "harmonic interactions registered");
    if (!angle || !perturbed) return;
    close(evaluate(*angle, s), reference, "nonbonded-only angle energy");
    require(perturbed->calculate_interactions(s.topo, s.conf, s.sim) == 0, "empty perturbed-angle list");
    compare_forces(s.conf.current().force, force, "nonbonded-only angle forces");
    compare_virial(s.conf.current().virial_tensor, virial, "nonbonded-only angle virial");
    require(ff.init(s.topo, s.conf, s.sim, std::cout, true) == 0, "full forcefield initialization");
    require(ff.calculate_interactions(s.topo, s.conf, s.sim) == 0, "full forcefield calculation");
    double total_angle = 0.0;
    for(unsigned int g=0; g<s.conf.current().energies.angle_energy.size(); ++g){
      total_angle += s.conf.current().energies.angle_energy[g];
      close(s.conf.current().perturbed_energy_derivatives.angle_energy[g], 0.0, "nonbonded-only angle derivative");
    }
    close(total_angle, reference, "full forcefield angle energy");

    const unsigned int count = s.topo.solute().angles().size();
    std::istringstream angles("TITLE\nHarmonic angle perturbation\nEND\n"
                              "PERTBONDANGLE\n1\n4 6 7 12 7\nEND\n"
                              "PERTBONDANGLEH\n1\n5 4 6 17 7\nEND\n");
    io::In_Perturbation angle_topology(angles);
    angle_topology.quiet = true;
    angle_topology.read(s.topo, s.sim.param());
    require(s.topo.perturbed_solute().angles().size() == 2 && s.topo.solute().angles().size() == count - 2,
            "explicit angle topology transfer");
    require(!io::messages.contains(io::message::error), "explicit angle topology input");
    require(std::isfinite(evaluate(*perturbed, s)), "explicit angle calculation");
    s.topo.perturbed_solute().softangles().push_back(topology::perturbed_three_body_term_struct(0,1,2,0,1));
    interaction::Forcefield unsupported;
    require(interaction::create_g96_bonded(unsupported, s.topo, s.sim.param(), reader, std::cout, true) != 0 &&
            io::messages.contains(io::message::error), "harmonic soft-angle rejection");
    io::messages.clear();
    s.sim.param().force.angle = 1;
    interaction::Forcefield cosine;
    require(interaction::create_g96_bonded(cosine, s.topo, s.sim.param(), reader, std::cout, true) == 0 &&
            cosine.interaction("PerturbedAngle"), "cosine soft-angle dispatch preserved");
    require(!io::messages.contains(io::message::error), "cosine dispatch diagnostics");
    io::messages.clear();
  }
}

int main(int argc, char **argv)
{
#ifdef XXMPI
  MPI_Init(&argc, &argv);
#endif
  check_numerics();
  check_linear_perturbation();
  check_dispatch();
  std::cout << "Maximum Cartesian force error: " << max_force_error << '\n'
            << "Maximum lambda derivative error: " << max_lambda_error << '\n'
            << "Failed checks: " << failures << '\n';
#ifdef XXMPI
  MPI_Finalize();
#endif
  return failures ? 1 : 0;
}
