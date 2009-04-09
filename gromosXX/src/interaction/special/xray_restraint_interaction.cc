/**
 * @file xray_restraint_interaction.cc
 * template methods of Xray_Restraint_Interaction
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <math/periodicity.h>

// special interactions
#include <interaction/interaction_types.h>

#include <interaction/special/xray_restraint_interaction.h>

#include <util/template_split.h>
#include <util/debug.h>
#include <vector>

// gsl headers
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <string>
#include <ios>


#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

interaction::Xray_Restraint_Interaction::Xray_Restraint_Interaction() : Interaction("XrayRestraint") {
}

interaction::Xray_Restraint_Interaction::~Xray_Restraint_Interaction() {
}

/**
 * calculate xray restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
void interaction::Xray_Restraint_Interaction::_calculate_xray_restraint_interactions
(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim) {

  m_timer.start();
  // get number of atoms in simulation
  const unsigned int atoms_size = topo.num_atoms();
  //update clipper atomvec
  for (unsigned int i = 0; i < atoms_size; i++) {
    atoms[i].set_coord_orth(clipper::Coord_orth(conf.current().pos(i)(0)*10.0,
            conf.current().pos(i)(1)*10.0,
            conf.current().pos(i)(2)*10.0));
  }
  // Calculate structure factors
  clipper::SFcalc_iso_fft<float> sfc;
  // run it
  sfc(fphi, atoms);

  // sqr_calc:       sum of squared Fcalc
  // obs:            sum of Fobs
  // calc:           sum of Fcalc
  // obs_calc:       sum of Fobs*Fcalc
  // obs_calcavg:    sum of Fobs*Fcalc(averaged)
  // obs_k_calcavg:  sum of Fobs-k_avg*Fcalc(averaged)
  // obs_k_calc:     sum of Fobs-k*Fcalc
  // sqr_calcavg:    sum of squared time-averaged Fcalc
  // calcavg:        sum of time-averaged Fcalc
  double sqr_calc = 0.0, obs = 0.0, calc = 0.0, obs_calc = 0.0, obs_k_calc = 0.0,
          sqr_calcavg = 0.0, calcavg = 0.0, obs_calcavg = 0.0, obs_k_calcavg = 0.0;
  // Number of reflections
  const unsigned int num_xray_rest = topo.xray_restraints().size();
  // e-term for time-average
  const double eterm = exp(-sim.time_step_size() / sim.param().xrayrest.tau);

  for (unsigned int i = 0; i < num_xray_rest; i++) {
    //filter calculated sf's
    clipper::HKL hkl(topo.xray_restraints()[i].h, topo.xray_restraints()[i].k, topo.xray_restraints()[i].l);
    conf.special().xray_rest[i].sf_curr = fabs(fphi[hkl].f());
    if (!sim.param().xrayrest.readavg && sim.steps() == 0) {
      // reset the averages at the beginning if requested
      conf.special().xray_rest[i].sf_av = conf.special().xray_rest[i].sf_curr;
    }
    conf.special().xray_rest[i].phase_curr = fphi[hkl].phi();
    conf.special().xray_rest[i].sf_av = fabs((1.0 - eterm) * conf.special().xray_rest[i].sf_curr + eterm * conf.special().xray_rest[i].sf_av);

    // calc sums
    obs_calc += conf.special().xray_rest[i].sf_curr * topo.xray_restraints()[i].sf;
    obs_calcavg += conf.special().xray_rest[i].sf_av * topo.xray_restraints()[i].sf;
    sqr_calc += conf.special().xray_rest[i].sf_curr * conf.special().xray_rest[i].sf_curr;
    obs += topo.xray_restraints()[i].sf;
    calc += conf.special().xray_rest[i].sf_curr;
    sqr_calcavg += conf.special().xray_rest[i].sf_av * conf.special().xray_rest[i].sf_av;
    calcavg += conf.special().xray_rest[i].sf_av;
  }


  //calc k_inst and k_avg
  double & k_inst = conf.special().xray.k_inst;
  k_inst = obs_calc / sqr_calc;
  double & k_avg = conf.special().xray.k_avg;
  k_avg = obs_calcavg / sqr_calcavg;
  DEBUG(10, "k_inst value: " << k_inst);
  DEBUG(10, "k_avg  value: " << k_avg);

  for (unsigned int i = 0; i < num_xray_rest; i++) {
    obs_k_calc += fabs(topo.xray_restraints()[i].sf - k_inst * conf.special().xray_rest[i].sf_curr);
    obs_k_calcavg += fabs(topo.xray_restraints()[i].sf - k_avg * conf.special().xray_rest[i].sf_av);
  }

  // calculate R_inst and R_avg
  double & R_inst = conf.special().xray.R_inst;
  R_inst = obs_k_calc / obs;
  double & R_avg = conf.special().xray.R_avg;
  R_avg = obs_k_calcavg / obs;
  DEBUG(10, "R_inst value: " << std::setw(15) << std::setprecision(8) << R_inst);
  DEBUG(10, "R_avg  value: " << std::setw(15) << std::setprecision(8) << R_avg);

  // calculate gradient
  D_k = std::complex<float> (0.0f, 0.0f); // zero it

  double energy_sum = 0.0;
  for (unsigned int i = 0; i < num_xray_rest; i++) {

    // SWITCH FOR DIFFERENT METHODS
    switch (sim.param().xrayrest.xrayrest) {
      case simulation::xrayrest_off :
      {
        break;
      }
      case simulation::xrayrest_inst :
      {
        // INSTANTANEOUS
        // calculate energy-sum
        const double term = topo.xray_restraints()[i].sf - k_inst * conf.special().xray_rest[i].sf_curr;
        energy_sum += term * term;
        // calculate derivatives of target function
        const double k_d = (topo.xray_restraints()[i].sf * sqr_calc - 2.0 * conf.special().xray_rest[i].sf_curr * obs_calc) / (sqr_calc * sqr_calc);
        const double dterm = (k_inst * conf.special().xray_rest[i].sf_curr - topo.xray_restraints()[i].sf)*(k_inst + conf.special().xray_rest[i].sf_curr * k_d);
        clipper::HKL hkl(topo.xray_restraints()[i].h, topo.xray_restraints()[i].k, topo.xray_restraints()[i].l);
        D_k.set_data(hkl, clipper::data32::F_phi(sim.param().xrayrest.force_constant * (dterm), conf.special().xray_rest[i].phase_curr));

        DEBUG(15, "SF inst: "
                << std::setw(5) << topo.xray_restraints()[i].h
                << std::setw(5) << topo.xray_restraints()[i].k
                << std::setw(5) << topo.xray_restraints()[i].l
                << std::setw(18) << topo.xray_restraints()[i].sf
                << std::setw(18) << k_inst * conf.special().xray_rest[i].sf_curr
                << std::setw(18) << conf.special().xray_rest[i].sf_curr
                << std::setw(18) << D_k[hkl].f() << std::setw(13) << D_k[hkl].phi());
        break;
      }
      case simulation::xrayrest_avg :
      {
        // TIMEAVERAGED
        // calculate energy-sum
        const double term = topo.xray_restraints()[i].sf - k_avg * conf.special().xray_rest[i].sf_av;
        energy_sum += term * term;
        // calculate derivatives of target function
        const double k_d = (topo.xray_restraints()[i].sf * sqr_calcavg - 2.0 * conf.special().xray_rest[i].sf_av * obs_calcavg) / (sqr_calcavg * sqr_calcavg);
        const double dterm = (k_avg * conf.special().xray_rest[i].sf_av - topo.xray_restraints()[i].sf)*((1.0 - eterm)*(k_avg + conf.special().xray_rest[i].sf_av * k_d));
        clipper::HKL hkl(topo.xray_restraints()[i].h, topo.xray_restraints()[i].k, topo.xray_restraints()[i].l);
        D_k.set_data(hkl, clipper::data32::F_phi(sim.param().xrayrest.force_constant * (dterm), conf.special().xray_rest[i].phase_curr));
        DEBUG(15, "SF avg: "
                << std::setw(5) << topo.xray_restraints()[i].h
                << std::setw(5) << topo.xray_restraints()[i].k
                << std::setw(5) << topo.xray_restraints()[i].l
                << std::setw(18) << topo.xray_restraints()[i].sf
                << std::setw(18) << k_avg * conf.special().xray_rest[i].sf_av
                << std::setw(18) << conf.special().xray_rest[i].sf_av
                << std::setw(18) << D_k[hkl].f() << std::setw(13) << D_k[hkl].phi());
        break;
      }
      case simulation::xrayrest_biq :
      {
        // BIQUADRATIC TIME-AVERAGED/INSTANTANEOUS
        // calculate energy-sum
        const double inst_term = topo.xray_restraints()[i].sf - k_inst * conf.special().xray_rest[i].sf_curr;
        const double av_term = topo.xray_restraints()[i].sf - k_avg * conf.special().xray_rest[i].sf_av;
        energy_sum += (inst_term * inst_term)*(av_term * av_term);
        // calculate derivatives of target function
        const double k_d_inst = (topo.xray_restraints()[i].sf * sqr_calc - 2.0 * conf.special().xray_rest[i].sf_curr * obs_calc) / (sqr_calc * sqr_calc);
        const double k_d_avg = (topo.xray_restraints()[i].sf * sqr_calcavg - 2.0 * conf.special().xray_rest[i].sf_av * obs_calcavg) / (sqr_calcavg * sqr_calcavg);
        const double dterm = ((k_inst * conf.special().xray_rest[i].sf_curr - topo.xray_restraints()[i].sf)*(topo.xray_restraints()[i].sf - k_avg * conf.special().xray_rest[i].sf_av)
                *(topo.xray_restraints()[i].sf - k_avg * conf.special().xray_rest[i].sf_av)*(k_inst + conf.special().xray_rest[i].sf_curr * k_d_inst))
                +((k_avg * conf.special().xray_rest[i].sf_av - topo.xray_restraints()[i].sf)*(topo.xray_restraints()[i].sf - k_inst * conf.special().xray_rest[i].sf_curr)
                *(topo.xray_restraints()[i].sf - k_inst * conf.special().xray_rest[i].sf_curr)*((1.0 - eterm)*(k_avg + conf.special().xray_rest[i].sf_curr * k_d_avg)));
        clipper::HKL hkl(topo.xray_restraints()[i].h, topo.xray_restraints()[i].k, topo.xray_restraints()[i].l);
        D_k.set_data(hkl, clipper::data32::F_phi(sim.param().xrayrest.force_constant * (dterm), conf.special().xray_rest[i].phase_curr));
        /*DEBUG(15, "SF avg: "
                << std::setw(5) << topo.xray_restraints()[i].h
                << std::setw(5) << topo.xray_restraints()[i].k
                << std::setw(5) << topo.xray_restraints()[i].l
                << std::setw(18) << topo.xray_restraints()[i].sf
                << std::setw(18) << k_avg * conf.special().xray_rest[i].sf_av
                << std::setw(18) << conf.special().xray_rest[i].sf_av
                << std::setw(18) << D_k[hkl].f() << std::setw(13) << D_k[hkl].phi());
         */
        break;
      }
      case simulation::xrayrest_loel :
      {
        break;
      }
    }
  }
  // finally calc energy
  conf.current().energies.xray_total = 0.5 * sim.param().xrayrest.force_constant * energy_sum;
  DEBUG(10, "energy: " << conf.current().energies.xray_total);

  // perform FFT
  d_r.fft_from(D_k);
  clipper::Coord_grid g0, g1;

  // 2.5 is hardcoded atom radius for grid sampling.
  clipper::Grid_range gd(fphi.base_cell(), d_r.grid_sampling(), 2.5);
  clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;
  math::Vec gradient(0.0, 0.0, 0.0);
  // calculate gradients of structure factors
  for (unsigned int i = 0; i < atoms_size; i++) {
    if (!atoms[i].is_null()) {
      gradient = math::Vec(0.0, 0.0, 0.0);
      clipper::AtomShapeFn sf(atoms[i].coord_orth(), atoms[i].element(),
              atoms[i].u_iso(), atoms[i].occupancy());
      sf.agarwal_params().resize(3);
      // conversion for clipper coordinate enum
      sf.agarwal_params()[0] = clipper::AtomShapeFn::X;
      sf.agarwal_params()[1] = clipper::AtomShapeFn::Y;
      sf.agarwal_params()[2] = clipper::AtomShapeFn::Z;
      // determine grid-ranges
      clipper::Coord_frac uvw = atoms[i].coord_orth().coord_frac(fphi.base_cell());
      g0 = uvw.coord_grid(d_r.grid_sampling()) + gd.min();
      g1 = uvw.coord_grid(d_r.grid_sampling()) + gd.max();
      i0 = clipper::Xmap<clipper::ftype32>::Map_reference_coord(d_r, g0);
      std::vector<clipper::ftype> rho_grad(3, 0.0f);
      clipper::ftype temp_rho = 0.0f;
      // loop over grid
      for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
        for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
          for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
            // get gradient from clipper
            sf.rho_grad(iw.coord_orth(), temp_rho, rho_grad);
            gradient(0) += d_r[iw] * rho_grad[0];
            gradient(1) += d_r[iw] * rho_grad[1];
            gradient(2) += d_r[iw] * rho_grad[2];
          }
        }
      } // loop over map
      // convert from Angstrom to nm
      gradient *= 10.0;
      // add to force
      conf.current().force(i) -= gradient;
      DEBUG(10, "grad(" << i << "): " << math::v2s(gradient));
    } // if atom not null
  } // for atoms

  m_timer.stop();

  // write xmap to external file
  if (sim.param().xrayrest.writexmap != 0 && sim.steps() % sim.param().xrayrest.writexmap == 0) {
    const clipper::Grid_sampling grid(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), fphi.base_hkl_info().resolution(), 1.5);
    clipper::Xmap<clipper::ftype32> density(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), grid);
    density.fft_from(fphi);
    clipper::CCP4MAPfile mapfile;
    std::ostringstream file_name, asu_file_name;
    file_name << "density_frame_" << std::setw(int(log10(sim.param().step.number_of_steps)))
            << std::setfill('0') << sim.steps() << ".ccp4";
    if (sim.param().xrayrest.writedensity == 1 || sim.param().xrayrest.writedensity == 3) {
      mapfile.open_write(file_name.str());
      mapfile.export_xmap(density);
      mapfile.close_write();
    }
    // Non Cristallographic Map
    clipper::NXmap<float> asu(density.grid_asu(), density.operator_orth_grid());
    asu_file_name << "density_asu_frame_" << std::setw(int(log10(sim.param().step.number_of_steps)))
            << std::setfill('0') << sim.steps() << ".ccp4";
    if (sim.param().xrayrest.writedensity == 2 || sim.param().xrayrest.writedensity == 3) {
      mapfile.open_write(asu_file_name.str());
      mapfile.export_nxmap(asu);
      mapfile.close_write();
    }
  }
}

int interaction::Xray_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim) {
  SPLIT_VIRIAL_BOUNDARY(_calculate_xray_restraint_interactions,
          topo, conf, sim);

  return 0;
}

int interaction::Xray_Restraint_Interaction::init(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os,
        bool quiet) {
  DEBUG(15, "Xray_Restraint_Interaction: init")
          const float sqpi2 = (math::Pi * math::Pi * 8.0 / 3.0);

  // Redirect clipper errors
  clipper::Message message;
  message.set_stream(os);

  // Construct clipper objects
  clipper::Spacegroup spacegr;
  try {
    clipper::Spgr_descr spgrinit(clipper::String(sim.param().xrayrest.spacegroup), clipper::Spgr_descr::HM);
    spacegr.init(spgrinit);
  } catch (const clipper::Message_fatal & msg) {
    io::messages.add("Xray_restraint_interaction", msg.text(), io::message::error);
    return 1;
  }

  clipper::Cell_descr cellinit(float(sim.param().xrayrest.cell_a * 10.0),
          float(sim.param().xrayrest.cell_b * 10.0),
          float(sim.param().xrayrest.cell_c * 10.0),
          float(sim.param().xrayrest.cell_alpha),
          float(sim.param().xrayrest.cell_beta),
          float(sim.param().xrayrest.cell_gamma));
  clipper::Cell cell(cellinit);

  clipper::Resolution reso(float(sim.param().xrayrest.resolution * 10.0));

  hkls.init(spacegr, cell, reso, true);
  fphi.init(hkls, hkls.cell());
  D_k.init(hkls, hkls.cell());
  // 1.5 is shannon-rate for oversampled FFT
  const clipper::Grid_sampling grid(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), fphi.base_hkl_info().resolution(), 1.5);
  d_r.init(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), grid);

  // Fill clipper atom-vector
  std::vector<clipper::Atom> atomvec;
  // Fill solute
  for (unsigned int i = 0; i < topo.num_solute_atoms(); i++) {
    clipper::Atom atm;
    atm.set_coord_orth(clipper::Coord_orth(0.0, 0.0, 0.0));
    atm.set_occupancy(1.0);
    atm.set_u_iso(sim.param().xrayrest.bfactor * 100 / sqpi2);
    DEBUG(1, "i: " << i << " size: " <<  topo.xray_elements().size());
    assert(i < topo.xray_elements().size());
    atm.set_element(topo.xray_elements()[i]);
    atomvec.push_back(atm);
  }
  // Fill solvent
  for (unsigned int i = 0; i < topo.num_solvent_atoms(); i++) {
    clipper::Atom atm;
    atm.set_coord_orth(clipper::Coord_orth(0.0, 0.0, 0.0));
    atm.set_occupancy(1.0);
    atm.set_u_iso(sim.param().xrayrest.bfactor * 100 / sqpi2);
    atm.set_element(topo.xray_solvelements()[i % topo.solvent(0).num_atoms()]);
    atomvec.push_back(atm);
  }
  atoms = clipper::Atom_list(atomvec);

  conf.special().xray_rest.resize(topo.xray_restraints().size());

  // Scale Fobs (needed for constant force-constant) -> scaled to sfscale
  const double sfscale = 100.0;
  double maxsf = 0.0;
  // Get max structure factor
  for (unsigned int i = 0; i < topo.xray_restraints().size(); i++) {
    if (maxsf < topo.xray_restraints()[i].sf)
      maxsf = topo.xray_restraints()[i].sf;
  }
  const double scalefactor = maxsf / sfscale;
  // scale
  for (unsigned int i = 0; i < topo.xray_restraints().size(); i++) {
    topo.xray_restraints()[i].sf /= scalefactor;
  }

  if (!quiet) {
    os.precision(2);
    os << "\nXRAYRESTINIT\n";
    os << "Restraint type              : ";
    switch (sim.param().xrayrest.xrayrest) {
      case simulation::xrayrest_off :
      {
        os << "No xray restraining\n";
        break;
      }
      case simulation::xrayrest_inst :
      {
        os << "Instantaneous xray restraining\n";
        break;
      }
      case simulation::xrayrest_avg :
      {
        os << "Time-averaged xray restraining\n";
        break;
      }
      case simulation::xrayrest_biq :
      {
        os << "Biquadratic time-averaged/instantaneous xray restraining\n";
        break;
      }
      case simulation::xrayrest_loel :
      {
        os << "Local-elevation xray restraining\n";
        break;
      }
      if (sim.param().xrayrest.readavg)
      os << "\treading xray averages from file\n";
    }
    os << "Restraint force-constant    : " << sim.param().xrayrest.force_constant << std::endl;
    os << "Spacegroup                  : " << sim.param().xrayrest.spacegroup << std::endl;
    os << "Cell                        : " << sim.param().xrayrest.cell_a
            << std::setw(8) << sim.param().xrayrest.cell_b
            << std::setw(8) << sim.param().xrayrest.cell_c
            << std::setw(8) << sim.param().xrayrest.cell_alpha
            << std::setw(8) << sim.param().xrayrest.cell_beta
            << std::setw(8) << sim.param().xrayrest.cell_gamma << std::endl;
    os.precision(4);
    os << "Resolution                  : " << sim.param().xrayrest.resolution << std::endl;
    os << "Bfactor                     : " << sim.param().xrayrest.bfactor << std::endl;
    os << "Num experimental reflections: " << topo.xray_restraints().size() << std::endl;
    os << "Num expected reflections    : " << hkls.num_reflections() << std::endl;
    os << "Writeing electron density   : " << sim.param().xrayrest.writexmap << std::endl;
    os << "END\n\n";
  }
  // Check if too low resolution
  double expnorm = 0.0, calcnorm = 0.0, tempnorm = 0.0;
  for (unsigned int i = 0; i < topo.xray_restraints().size(); i++) {
    // calc max. experimental-index norm
    tempnorm = sqrt(double((topo.xray_restraints()[i].h * topo.xray_restraints()[i].h)+
            (topo.xray_restraints()[i].k * topo.xray_restraints()[i].k)+
            (topo.xray_restraints()[i].l * topo.xray_restraints()[i].l)));
    if (tempnorm > expnorm)
      expnorm = tempnorm;
  }
  for (int i = 0; i < hkls.num_reflections(); i++) {
    // calc max. calculation-index norm
    tempnorm = sqrt(double((hkls.hkl_of(i).h() * hkls.hkl_of(i).h())+
            (hkls.hkl_of(i).k() * hkls.hkl_of(i).k())+
            (hkls.hkl_of(i).l() * hkls.hkl_of(i).l())));
    if (tempnorm > calcnorm)
      calcnorm = tempnorm;
  }

  if (expnorm > calcnorm) {
    io::messages.add("Xray_restraint_interaction", "Too little reflections. Set higher resolution!", io::message::error);
    return 1;
  }

  return 0;
}

