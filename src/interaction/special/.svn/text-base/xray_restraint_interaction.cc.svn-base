/**
 * @file xray_restraint_interaction.cc
 * template methods of Xray_Restraint_Interaction
 */
#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

// special interactions
#include "../../interaction/interaction_types.h"
#include "../../util/umbrella.h"
#include "../../util/umbrella_weight.h"
#include "../../interaction/special/xray_restraint_interaction.h"

#include "../../math/periodicity.h"
#include "../../util/template_split.h"
#include "../../util/debug.h"
#include <set>
#include <vector>

#ifdef HAVE_CLIPPER
#include "../../interaction/special/xray/dens.h"
#include "../../interaction/special/xray/sf.h"
#include "../../interaction/special/xray/bfactor.h"
using namespace interaction::xray;
#include <clipper/clipper.h>
#include <limits>
#endif




#ifdef OMP
#include <omp.h>
#endif

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
int interaction::Xray_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim) {
#ifdef HAVE_CLIPPER
    const double & to_ang = sim.param().xrayrest.to_angstrom;
    clipper::Atom_list & atoms = conf.special().xray_conf.atoms;
    clipper::Atom_list & atoms_sf = conf.special().xray_conf.atoms_sf;
    clipper::Xmap<clipper::ftype32> & rho_calc = conf.special().xray_conf.rho_calc;
    clipper::Xmap<clipper::ftype32> & rho_obs = conf.special().xray_conf.rho_obs;
    clipper::HKL_data<clipper::data32::F_phi> & fphi_calc = conf.special().xray_conf.fphi_calc;
    clipper::HKL_data<clipper::data32::F_phi> & fphi = conf.special().xray_conf.fphi;
    clipper::HKL_data<clipper::data32::F_phi> & fphi_obs = conf.special().xray_conf.fphi_obs;
    clipper::FFTmap_p1 & D_k = conf.special().xray_conf.D_k;
    clipper::Xmap<clipper::ftype32> & d_r = conf.special().xray_conf.d_r;
    clipper::Spacegroup & sym_spacegroup = conf.special().xray_conf.sym_spacegroup;
    const clipper::Cell & cell = rho_calc.cell();
  m_timer.start();
  // get number of atoms in simulation
  const int atoms_size = topo.num_atoms();
  // update clipper atomvec: convert the position to Angstrom
  math::Periodicity<math::triclinic> periodicity(conf.current().box);

  bool update = false;
  unsigned int update_step = sim.param().xrayrest.structure_factor_calculation.steps_nb_constant;
  if (!sim.steps() || (update_step != 0 && sim.steps() % update_step == 0))
    update = true;

  if (conf.special().change_on_slave == 1){
    update = true;
    conf.special().change_on_slave = 0;
    std::cout << "lambda value changed on slave" << std::endl;
  }
  if (conf.special().change_on_slave == 2){
    update = true;
    conf.special().change_on_slave = 0;
    std::cout << "energy recalculated on slave" << std::endl;
  }
  const double tol2 = sim.param().xrayrest.structure_factor_calculation.atom_move_tolerance *
          sim.param().xrayrest.structure_factor_calculation.atom_move_tolerance;
  for (int i = 0; i < atoms_size && update == false; i++) {
    const math::Vec & new_pos = conf.current().pos(i);
    const math::Vec old_pos(atoms_sf[i].coord_orth().x() / to_ang,
            atoms_sf[i].coord_orth().y() / to_ang,
            atoms_sf[i].coord_orth().z() / to_ang);
    math::Vec ni;
    periodicity.nearest_image(new_pos, old_pos, ni);
    if (math::abs2(ni) > tol2)
      update = true;
  }

  for (int i = 0; i < atoms_size; i++) {
    math::Vec in_box = conf.current().pos(i);
    periodicity.put_into_positive_box(in_box);
    in_box *= to_ang;
    atoms[i].set_coord_orth(clipper::Coord_orth(in_box(0),
            in_box(1), in_box(2)));
  }

  // compute the new force constant from the lambda value
  if(sim.param().xrayrest.replica_exchange_parameters.switcher == simulation::replica_exchange_force){
    double min = sim.param().xrayrest.replica_exchange_parameters.lambda_dependant_min;
    double max = sim.param().xrayrest.replica_exchange_parameters.lambda_dependant_max;
    sim.param().xrayrest.force_constant = topo.lambda() * (max - min) + min;
  }

  // compute the resolution from the lambda value
  if (sim.param().xrayrest.replica_exchange_parameters.switcher == simulation::replica_exchange_resolution){
    double min = sim.param().xrayrest.replica_exchange_parameters.lambda_dependant_min;
    double max = sim.param().xrayrest.replica_exchange_parameters.lambda_dependant_max;
    sim.param().xrayrest.resolution = topo.lambda() * (max - min) + min;
  }

  // lets start with sym stuff.
  if (sim.param().xrayrest.symrest != simulation::xray_symrest_off) {
    DEBUG(6, "symmetry restraints");
    m_timer.start("symmetry restraints");
    std::vector<unsigned int>::const_iterator
            atom_it = topo.xray_sym_restraints().begin(),
            atom_to = topo.xray_sym_restraints().end();

    for (; atom_it != atom_to; ++atom_it) {
      const unsigned int atom_p = *atom_it - topo.xray_asu()[0];
      DEBUG(6, "atom: " << *atom_it);

      switch (sim.param().xrayrest.symrest) {
        case simulation::xray_symrest_off : break;
        case simulation::xray_symrest_ind :
        {
          DEBUG(6, "symmetry individual atoms");
          for (unsigned int i = 0; i < topo.xray_asu().size()-1; ++i) {
            const unsigned int atom_i = topo.xray_asu()[i] + atom_p;
            const clipper::Coord_orth pos_i_ang(sym_spacegroup.symop(i).inverse().rtop_orth(cell) * atoms[atom_i].coord_orth());
            math::Vec pos_i(pos_i_ang.x(), pos_i_ang.y(), pos_i_ang.z());
            pos_i /= to_ang;
            for (unsigned int j = i + 1; j < topo.xray_asu().size(); ++j) {
              const unsigned int atom_j = topo.xray_asu()[j] + atom_p;
              const clipper::Coord_orth pos_j_ang(sym_spacegroup.symop(j).inverse().rtop_orth(cell) * atoms[atom_j].coord_orth());
              math::Vec pos_j(pos_j_ang.x(), pos_j_ang.y(), pos_j_ang.z());
              pos_j /= to_ang;

              DEBUG(8, "i: " << i << " j: " << j);
              DEBUG(8, "pos i   : " << math::v2s(pos_i));
              DEBUG(8, "pos j   : " << math::v2s(pos_j));

              math::Vec dist;
              periodicity.nearest_image(pos_i, pos_j, dist);
              DEBUG(9, "dist    : " << math::v2s(dist));
              const clipper::Coord_orth dist_cl(dist(0), dist(1), dist(2));
              // do the rotation
              const clipper::Coord_orth r_i(sym_spacegroup.symop(i).rtop_orth(cell).rot() * dist_cl);
              const clipper::Coord_orth r_j(sym_spacegroup.symop(j).rtop_orth(cell).rot() * dist_cl);

              const math::Vec f_i(-sim.param().xrayrest.sym_force_constant * math::Vec(r_i.x(), r_i.y(), r_i.z()));
              const math::Vec f_j(sim.param().xrayrest.sym_force_constant * math::Vec(r_j.x(), r_j.y(), r_j.z()));
              
              DEBUG(8, "f i     : " << math::v2s(f_i));
              DEBUG(8, "f j     : " << math::v2s(f_j));

              conf.current().force(atom_i) += f_i;
              conf.current().force(atom_j) += f_j;
              const double V = 0.5 * sim.param().xrayrest.sym_force_constant * math::abs2(dist);
              conf.current().energies.xray_total += V;
              DEBUG(7, "energy : " << V);
            }
          } // loop over images
          break;
        }
        case simulation::xray_symrest_constr :
        {
          DEBUG(6, "symmetry constrain atoms");
          // loop over images
          for (unsigned int i = 1; i < topo.xray_asu().size(); ++i) {
            const unsigned int atom_img = topo.xray_asu()[i] + atom_p;
            // optain the image position
            const clipper::Coord_orth pos_img_ang(sym_spacegroup.symop(i).rtop_orth(cell) * atoms[*atom_it].coord_orth());
            math::Vec pos_img(pos_img_ang.x(), pos_img_ang.y(), pos_img_ang.z());
            pos_img /= to_ang;
            DEBUG(8, "pos     : " << math::v2s(conf.current().pos(atom_img)));
            periodicity.put_into_positive_box(pos_img);
            DEBUG(8, "new pos : " << math::v2s(pos_img));
            conf.current().pos(atom_img) = pos_img;
            pos_img *= to_ang;
            atoms[atom_img].set_coord_orth(clipper::Coord_orth(pos_img(0),
                    pos_img(1), pos_img(2)));
          } // loop over images
          break;
        }
      } // method switch
    } // loop over restrained atoms
    m_timer.stop("symmetry restraints");
  }
  
  // fit the b factor
  if (sim.param().xrayrest.bfactor.step && (!sim.steps() || sim.steps() % sim.param().xrayrest.bfactor.step == 0)) {
    m_timer.start("B factor fitting");
    fit_bfactor(topo, conf, sim, cell, atoms, fphi_calc, fphi, fphi_obs, rho_calc, D_k, d_r);
    update = true;
    m_timer.stop("B factor fitting");
  }

  if (update) {
    atoms_sf = atoms;
    // Calculate structure factors
    m_timer.start("structure factor");
    calculate_electron_density(rho_calc, atoms);
    // FFT the electron density to obtain the structure factors
    rho_calc.fft_to(fphi_calc);
    m_timer.stop("structure factor");

    if (sim.param().xrayrest.overall_bfactor.B_overall_switcher == simulation::B_overall_on) {
      m_timer.start("overall B factor fitting");
      fit_overall_bfactor(topo, conf, sim, cell, atoms, fphi_calc, fphi, fphi_obs, rho_calc, D_k, d_r);
      m_timer.stop("overall B factor fitting");
    }

    m_timer.start("scaling");
    scale_sf(topo, conf, sim, fphi_calc, fphi, fphi_obs);
    m_timer.stop("scaling");
  }

  math::VArray force(atoms_size);
  force = 0.0;

  if (sim.param().xrayrest.local_elevation) {
    ///////////////////////////////////////////////
    // LOCAL ELEVATION
    ///////////////////////////////////////////////
    if (update) {
      m_timer.start("obs. electron density");
      rho_obs.fft_from(fphi_obs);
      m_timer.stop("obs. electron density");
    }
  } else if (sim.param().xrayrest.mode != simulation::xrayrest_mode_electron_density) {
    ///////////////////////////////////////////////
    // STRUCTURE FACTOR RESTRAINING
    ///////////////////////////////////////////////
    if (update) {
      m_timer.start("energy");
      calculate_energy_sf(sim, fphi,
              topo.xray_restraints(),
              conf.special().xray_rest,
              sim.param().xrayrest.xrayrest,
              conf.special().xray.k_inst, conf.special().xray.k_avg,
              D_k, sim.param().xrayrest.force_constant,
              conf.current().energies.xray_total);
      m_timer.stop("energy");
    } else {
      conf.current().energies.xray_total = conf.old().energies.xray_total;
    }

    // start to calculate the forces
    m_timer.start("force");
    math::SArray b_deriv(atoms.size());
    calculate_force_sf(update, D_k, d_r, atoms, force, b_deriv, to_ang);
    for(unsigned int i = 0; i < atoms.size(); ++i) {
      conf.special().xray_bfoc[i].b_factor_gradient = b_deriv(i);
    }
    m_timer.stop("force");
  } else {
    ///////////////////////////////////////////////
    // ELECTRON DENSITY RESTRAINING
    ///////////////////////////////////////////////
    m_timer.start("energy & force");
    if (update)
      rho_obs.fft_from(fphi_obs);
    calculate_energy_rho(atoms, rho_obs, rho_calc, sim.param().xrayrest.force_constant,
            conf.current().energies.xray_total, force, to_ang);
    m_timer.stop("energy & force");
  }
  DEBUG(10, "energy: " << conf.current().energies.xray_total);

  // add the force and calculate some statistics
  double sum_force_phys = 0.0, sum_force_xray = 0.0;
  double sum_force_phys2 = 0.0, sum_force_xray2 = 0.0;

  for(int i = 0; i < atoms_size; ++i) {
    const double f_phys = math::abs(conf.current().force(i));
    sum_force_phys += f_phys;
    sum_force_phys2 += f_phys * f_phys;
    const double f_xray = math::abs(force(i));
    sum_force_xray += f_xray;
    sum_force_xray2 += f_xray * f_xray;
    conf.current().force(i) += force(i);
  }
  const double mean_force_phys = sum_force_phys / atoms_size;
  const double stddev_force_phys = sqrt((sum_force_phys2 - mean_force_phys * mean_force_phys) / atoms_size);
  const double mean_force_xray = sum_force_xray / atoms_size;
  const double stddev_force_xray = sqrt((sum_force_xray2 - mean_force_xray * mean_force_xray) / atoms_size);

  if (sim.param().print.stepblock && sim.steps() % sim.param().print.stepblock == 0) {
    std::cout << "XRAY RESTRAINING\n"
            << "Lambda               : " << std::setw(4) << topo.lambda() << std::endl
            << "Resolution (in use)  : " << std::setw(4) << sim.param().xrayrest.resolution << std::endl
            << "Reflections (in use) : " << std::setw(4) << topo.xray_restraints().size() << std::endl
            << "Energy               : " << std::setw(13) << conf.current().energies.xray_total << std::endl
            << "R value (inst.)      : " << std::setw(4) << conf.special().xray.R_inst << std::endl
            << "free R value (inst.) : " << std::setw(4) << conf.special().xray.R_free_inst << std::endl
            << "R value (avg.)       : " << std::setw(4) << conf.special().xray.R_avg << std::endl
            << "free R value (avg.)  : " << std::setw(4) << conf.special().xray.R_free_avg << std::endl
            << std::endl
            << "physical mean force   : " << std::setw(13) << mean_force_phys
            << " +/- " << std::setw(13) << stddev_force_phys << std::endl
            << "restraints mean force : " << std::setw(13) << mean_force_xray
            << " +/- " << std::setw(13) << stddev_force_xray << std::endl
            << "ratio phys/rest       : " << std::setw(8) << mean_force_phys / mean_force_xray << std::endl
            << "END" << std::endl;
  }

  m_timer.stop();

  // write xmap to external file
  if (sim.param().xrayrest.writexmap != 0 && sim.steps() % sim.param().xrayrest.writexmap == 0) {
    if (sim.param().xrayrest.mode != simulation::xrayrest_mode_electron_density)
      rho_obs.fft_from(fphi_obs);
    clipper::CCP4MAPfile mapfile;
    std::ostringstream file_name, asu_file_name;
    file_name << "density_frame_" << std::setw(int(log10(sim.param().step.number_of_steps)))
            << std::setfill('0') << sim.steps() << ".ccp4";
    if (sim.param().xrayrest.writedensity == 1 || sim.param().xrayrest.writedensity == 3) {
      mapfile.open_write(file_name.str());
      mapfile.export_xmap(rho_obs);
      mapfile.close_write();
    }
    // Non Cristallographic Map
    clipper::NXmap<double> asu(rho_obs.grid_asu(), rho_obs.operator_orth_grid());
    asu_file_name << "density_asu_frame_" << std::setw(int(log10(sim.param().step.number_of_steps)))
            << std::setfill('0') << sim.steps() << ".ccp4";
    if (sim.param().xrayrest.writedensity == 2 || sim.param().xrayrest.writedensity == 3) {
      mapfile.open_write(asu_file_name.str());
      mapfile.export_nxmap(asu);
      mapfile.close_write();
    }
  }

  // increase the local elevation thresholds here if required
  if (sim.steps() > 0) {
    std::vector<topology::xray_umbrella_weight_struct>::iterator it =
            topo.xray_umbrella_weights().begin(), to =
            topo.xray_umbrella_weights().end();
    for (; it != to; ++it) {
      if (!it->threshold_freeze) {
        // increase it
        it->threshold += it->threshold_growth_rate;
        if (it->threshold < 0.0) {
          it->threshold = 0.0;
          it->threshold_freeze = true;
          std::ostringstream msg;
          msg << "Threshold of umbrella " << it->id << " was zero and thus frozen.";
          io::messages.add(msg.str(), "Xray_Restraint_Interaction", io::message::warning);
        }
      }

      // now let's take care of signals
      switch(it->signal) { 
        case 1: 
        {
          // raise the threshold a little more and freeze
          it->threshold += it->threshold_overshoot;
          it->threshold_freeze = true;
          // set all umbrella weights to zero again.
          std::vector<util::Umbrella>::iterator umb_it = conf.special().umbrellas.begin(),
             umb_to = conf.special().umbrellas.end();
          for (; umb_it != umb_to; ++umb_it) {
            if (umb_it->id == it->id) {
              umb_it->configurations.clear();
              break;
            }
          }
          // reset signal
          it->signal = 0;
          break;
        }
        default: break;
      } // switch signal
    } // for umbrellas
  } // steps?

  // restore the restraints data in the topology from the backup
#endif
  return 0;
}

int interaction::Xray_Restraint_Interaction::init(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os,
        bool quiet) {

#ifdef HAVE_CLIPPER
  DEBUG(15, "Xray_Restraint_Interaction: init")
  const double sqpi2 = (math::Pi * math::Pi * 8.0);
  clipper::Atom_list & atoms = conf.special().xray_conf.atoms;
  clipper::Atom_list & atoms_sf = conf.special().xray_conf.atoms_sf;
  clipper::Xmap<clipper::ftype32> & rho_calc = conf.special().xray_conf.rho_calc;
  clipper::Xmap<clipper::ftype32> & rho_obs = conf.special().xray_conf.rho_obs;
  clipper::HKL_info & hkls = conf.special().xray_conf.hkls;
  clipper::HKL_data<clipper::data32::F_phi> & fphi = conf.special().xray_conf.fphi;
  clipper::HKL_data<clipper::data32::F_phi> & fphi_calc = conf.special().xray_conf.fphi_calc;
  clipper::HKL_data<clipper::data32::F_phi> & fphi_obs = conf.special().xray_conf.fphi_obs;
  clipper::FFTmap_p1 & D_k = conf.special().xray_conf.D_k;
  clipper::Xmap<clipper::ftype32> & d_r = conf.special().xray_conf.d_r;
  clipper::Spacegroup & sym_spacegroup = conf.special().xray_conf.sym_spacegroup;

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

  const double & to_ang = sim.param().xrayrest.to_angstrom;

  math::Box & box = conf.current().box;
  const double a = math::abs(box(0));
  const double b = math::abs(box(1));
  const double c = math::abs(box(2));
  const double alpha = acos(math::costest(dot(box(1), box(2)) / (b * c)));
  const double beta = acos(math::costest(dot(box(0), box(2)) / (a * c)));
  const double gamma = acos(math::costest(dot(box(0), box(1)) / (a * b)));
  DEBUG(6, "a: " << a << " b: " << b << " c: " << c << " al.: " << alpha <<
          " be.: " << beta << " ga.: " << gamma);
  clipper::Cell_descr cellinit(a * to_ang, b * to_ang, c * to_ang,
          alpha, beta, gamma);
  clipper::Cell cell(cellinit);

  clipper::Resolution reso(sim.param().xrayrest.resolution * to_ang);
  DEBUG(6, "resolution: " << sim.param().xrayrest.resolution );


    // create a grid and a crystalographic map
  const double shannon_rate = 1.5;
  const clipper::Grid_sampling grid_rho_calc(spacegr, cell, reso, shannon_rate);
  rho_calc.init(spacegr, cell, grid_rho_calc);
  const clipper::Grid_sampling grid_rho_obs(spacegr, cell, reso, shannon_rate);
  rho_obs.init(spacegr, cell, grid_rho_obs);

  hkls.init(spacegr, cell, reso, true);
  fphi.init(hkls, hkls.cell());
  fphi_calc.init(hkls, hkls.cell());
  fphi_obs.init(hkls, hkls.cell());

  // The difference map has to be a P 1 map in order to get agreement with
  // the finite difference results. However, the reasons for this are
  // not 100% clear. In principle (from theory) a spacegroup depdendent
  // Xmap should do the job.
  clipper::Spgr_descr spgrinit(clipper::String("P 1"), clipper::Spgr_descr::HM);
  clipper::Spacegroup p1_spacegr;
  p1_spacegr.init(spgrinit);
  // create a grid and a P1 FFT map. Here we can use the FFTmap_p1 which was
  // designed for fast P1.
  // 1.5 is shannon-rate for oversampled FFT
  const clipper::Grid_sampling fftgrid(p1_spacegr, cell, reso, 1.5);
  D_k.init(fftgrid);
  // we still need an Xmap for convenient looping over the data in order
  // to do the convolution.
  const clipper::Grid_sampling grid_d_r(spacegr, cell, reso, 1.5);
  d_r.init(spacegr, cell, grid_d_r);

  // Fill clipper atom-vector
  std::vector<clipper::Atom> atomvec;

  if (sim.param().xrayrest.bfactor.init)
     conf.special().xray_bfoc.resize(topo.num_atoms());

  // Fill solute
  unsigned int b_index = 0;
  for (unsigned int i = 0; i < topo.num_solute_atoms(); i++, ++b_index) {
    clipper::Atom atm;
    assert(i < topo.xray_occupancies().size());
    if (sim.param().xrayrest.bfactor.init) {
      conf.special().xray_bfoc[b_index].occupancy = topo.xray_occupancies()[i];
    }
    atm.set_occupancy(conf.special().xray_bfoc[b_index].occupancy);
    atm.set_coord_orth(clipper::Coord_orth(0.0, 0.0, 0.0));
    assert(i < topo.xray_b_factors().size());

    conf.special().xray_bfoc[b_index].b_factor_gradient = 0.0;
    if (sim.param().xrayrest.bfactor.init) {
      conf.special().xray_bfoc[b_index].b_factor = topo.xray_b_factors()[i];
    }
    atm.set_u_iso(conf.special().xray_bfoc[b_index].b_factor * to_ang * to_ang / sqpi2);
    assert(i < topo.xray_elements().size());
    atm.set_element(topo.xray_elements()[i]);
    atomvec.push_back(atm);
  }
  // Fill solvent
  for (unsigned int i = 0; i < topo.num_solvent_atoms(); i++, ++b_index) {
    clipper::Atom atm;
    unsigned int index = i % topo.solvent(0).num_atoms();
    assert(index < topo.xray_solv_occupancies().size());
    if (sim.param().xrayrest.bfactor.init) {
      conf.special().xray_bfoc[b_index].occupancy = topo.xray_solv_occupancies()[index];
    }
    atm.set_occupancy(conf.special().xray_bfoc[b_index].occupancy);
    atm.set_coord_orth(clipper::Coord_orth(0.0, 0.0, 0.0));
    assert(index < topo.xray_solv_b_factors().size());
    
    if (sim.param().xrayrest.bfactor.init) {
      conf.special().xray_bfoc[b_index].b_factor = topo.xray_solv_b_factors()[index];
    }
    atm.set_u_iso(conf.special().xray_bfoc[b_index].b_factor * to_ang * to_ang / sqpi2);

    conf.special().xray_bfoc[b_index].b_factor_gradient = 0.0;
    assert(index < topo.xray_solvelements().size());
    atm.set_element(topo.xray_solvelements()[index]);
    atomvec.push_back(atm);
  }
  atoms = clipper::Atom_list(atomvec);
  atoms_sf = atoms;

  if (sim.param().xrayrest.local_elevation) {
    std::vector<topology::xray_umbrella_weight_struct>::iterator it =
            topo.xray_umbrella_weights().begin(), to =
            topo.xray_umbrella_weights().end();
    for (; it != to; ++it) {
      std::vector<util::Umbrella>::iterator umb_it = conf.special().umbrellas.begin(),
              umb_to = conf.special().umbrellas.end();
      bool found = false;
      for (; umb_it != umb_to; ++umb_it) {
        if (umb_it->id == it->id) {
          found = true;
          umb_it->umbrella_weight_factory = new interaction::Electron_Density_Umbrella_Weight_Factory(*it, conf, atoms, rho_calc, rho_obs, to_ang);
        }

      }
      if (!found) {
        std::ostringstream msg;
        msg << "Cannot find umbrella " << it->id;
        io::messages.add(msg.str(), "Xray_Restraint_Interaction", io::message::error);
        return 1;
      }
    } // for weights
  } // if local elev

  if (sim.param().xrayrest.symrest != simulation::xray_symrest_off) {
    try {
      clipper::Spgr_descr spgrinit(clipper::String(sim.param().xrayrest.sym_spacegroup), clipper::Spgr_descr::HM);
      sym_spacegroup.init(spgrinit);
    } catch (const clipper::Message_fatal & msg) {
      io::messages.add("Xray_Restraint_Interaction", "symmetry spacegroup: " + msg.text(), io::message::error);
      return 1;
    }
  } // if symmetry rest

  if (sim.param().xrayrest.bfactor.step) {
    conf.special().xray_bfoc.resize(atoms.size());
  }

  // remove HKLs of too high resolution
  double min_reso = std::numeric_limits<double>::min();
  double max_reso = std::numeric_limits<double>::max();
  std::vector<topology::xray_restraint_struct> refl;
  unsigned int filtered = 0;

  for (unsigned int i = 0; i < topo.xray_restraints().size(); i++) {
    // calc max. experimental-index norm
    const clipper::HKL hkl(topo.xray_restraints()[i].h, topo.xray_restraints()[i].k, topo.xray_restraints()[i].l);
    double hkl_resolution = sqrt(1.0 / hkl.invresolsq(cell)) / to_ang;
    if (hkl_resolution >= sim.param().xrayrest.resolution) {
      refl.push_back(topo.xray_restraints()[i]);
      min_reso = std::max(min_reso, hkl_resolution);
      max_reso = std::min(max_reso, hkl_resolution);
    } else {
      ++filtered;
    }
  }
  topo.xray_restraints() = refl;
  std::vector<topology::xray_restraint_struct> refl_free;
  unsigned int filtered_free = 0;
  for (unsigned int i = 0; i < topo.xray_rfree().size(); i++) {
    // calc max. experimental-index norm
    const clipper::HKL hkl(topo.xray_rfree()[i].h, topo.xray_rfree()[i].k, topo.xray_rfree()[i].l);
    double hkl_resolution = sqrt(1.0 / hkl.invresolsq(cell)) / to_ang;
    if (hkl_resolution >= sim.param().xrayrest.resolution) {
      refl_free.push_back(topo.xray_rfree()[i]);
      min_reso = std::max(min_reso, hkl_resolution);
      max_reso = std::min(max_reso, hkl_resolution);
    } else {
      filtered_free++;
    }
  }
  topo.xray_rfree() = refl_free;

  double calc_max_reso = std::numeric_limits<double>::max();
  for (int i = 0; i < hkls.num_reflections(); i++) {
    // calc max. calculation-index norm
    double hkl_resolution = sqrt(1.0 / hkls.hkl_of(i).invresolsq(cell)) / to_ang;
    calc_max_reso = std::min(calc_max_reso, hkl_resolution);
  }

  conf.special().xray_rest.resize(topo.xray_restraints().size() + topo.xray_rfree().size());

  conf.special().xray.B_overall = sim.param().xrayrest.overall_bfactor.init;

  if (!quiet) {
    os.precision(2);
    os << "\nXRAYREST\n";
    os << "Restraint type              : ";
    switch (sim.param().xrayrest.xrayrest) {
      case simulation::xrayrest_off :
        os << "No xray restraining";
        break;
      case simulation::xrayrest_inst :
        os << "Instantaneous restraining";
        break;
      case simulation::xrayrest_avg :
        os << "Time-averaged restraining";
        break;
      case simulation::xrayrest_biq :
        os << "Biquadratic time-averaged/instantaneous restraining";
        break;
    }
    if (sim.param().xrayrest.mode == simulation::xrayrest_mode_electron_density) {
      os << " on the electron density";
      if (sim.param().xrayrest.local_elevation)
        os << " using local elevation";
      os << ".\n";
    } else {
      os << " on the structure factors.\n";
    }

    if (sim.param().xrayrest.readavg)
      os << "\treading xray averages from file\n";

    os << "Restraint force-constant    : " << sim.param().xrayrest.force_constant << std::endl;
    os << "Spacegroup                  : " << sim.param().xrayrest.spacegroup << std::endl;
    os.precision(4);
    os << "Resolution                  : " << sim.param().xrayrest.resolution << std::endl;
    os << "Num experimental reflections: " << topo.xray_restraints().size() << std::endl
       << "                              " << filtered << " were filtered away due to requested resolution." << std::endl;
    os << "Num R-free reflections      : " << topo.xray_rfree().size() << std::endl
       << "                              " << filtered_free << " were filtered away due to requested resolution." << std::endl;
    os << "Num expected reflections    : " << hkls.num_reflections() << std::endl;
    os.precision(8);
    os << "Max. resolution of data     : " << max_reso << std::endl;
    os << "Min. resolution of data     : " << min_reso << std::endl;
    os << "Max. resolution calculated  : " << calc_max_reso << std::endl;
    if (sim.param().xrayrest.overall_bfactor.B_overall_switcher == simulation::B_overall_on)
      os << "Overall B factor            : " << conf.special().xray.B_overall << std::endl;
    else
      os << "Overall B factor            : not used." << std::endl;
    os << "Writeing electron density   : " << sim.param().xrayrest.writexmap << std::endl << std::endl;
    if (sim.param().xrayrest.local_elevation) {
      os << "The following local elevation umbrellas are weighted by the electron density:" << std::endl;
      std::vector<topology::xray_umbrella_weight_struct>::const_iterator it =
              topo.xray_umbrella_weights().begin(), to =
              topo.xray_umbrella_weights().end();
      for (; it != to; ++it) {
        os.precision(5);
        os << "  - " << std::setw(5) << it->id << ": Threshold " << std::setw(15) << it->threshold << " atoms: ";
        for (std::vector<unsigned int>::const_iterator a_it = it->atoms.begin(),
                a_to = it->atoms.end(); a_it != a_to; ++a_it) {
          os << std::setw(5) << (*a_it) + 1;
        }
        os << std::endl;
      }
    }

    if (sim.param().xrayrest.symrest != simulation::xray_symrest_off) {
      os << std::endl;
      switch (sim.param().xrayrest.symrest) {
        case simulation::xray_symrest_off : break;
        case simulation::xray_symrest_ind:
          os << "symmetry restraints on individual atom positions" << std::endl;
          break;
        case simulation::xray_symrest_constr:
          os << "symmetry constraints on individual atom positions" << std::endl;
          break;
      }
      os << " - force constant : " << sim.param().xrayrest.sym_force_constant << std::endl;
      os << " - Spacegroup     : " << sim.param().xrayrest.sym_spacegroup << std::endl;
      os << " - " << topo.xray_asu().size() << " asymmetric units. First atoms: " << std::endl;
      for (unsigned int i = 0; i < topo.xray_asu().size(); ++i) {
        os << "       - " << topo.xray_asu()[i]+1 << std::endl;
      }
      os << " - the restraint atoms: " << std::endl;
      for(unsigned int i = 0; i < topo.xray_sym_restraints().size(); ++i) {
        os << std::setw(6) << (topo.xray_sym_restraints()[i]+1);
        if ((i+1) % 10 == 0) os << std::endl;
      }
      os << std::endl;
    }

    if (sim.param().xrayrest.bfactor.step) {
      os << "B-factor optimisation: " << std::endl
              << "  - every " << sim.param().xrayrest.bfactor.step << " steps." << std::endl
              << "  - terminate after  " << sim.param().xrayrest.bfactor.terminate_iterations
              << " iterations or if |dR/dB| < " << sim.param().xrayrest.bfactor.terminate_gradient << std::endl
              << "  - boundaries for B factor: [" << sim.param().xrayrest.bfactor.min << " - "
              << sim.param().xrayrest.bfactor.max << "]" << std::endl
              << "  - number of groups: " << sim.param().xrayrest.bfactor.groups.size() << std::endl;
      for(unsigned int i = 0; i < sim.param().xrayrest.bfactor.groups.size(); ++i) {
        std::set<unsigned int>::const_iterator it = sim.param().xrayrest.bfactor.groups[i].begin(),
                to = sim.param().xrayrest.bfactor.groups[i].end();
        os << "        - ";
        unsigned int size = sim.param().xrayrest.bfactor.groups[i].size();
        for(unsigned int num = 1; it != to; ++it, ++num) {
          os << std::setw(6) << *it+1;
          if (num % 10 == 0 && num != size)
            os << std::endl << "          ";
        }
        os << std::endl;
      }
    }

    os << "END" << std::endl;
  }

  if (max_reso < calc_max_reso && fabs(max_reso - calc_max_reso) > math::epsilon) {
    io::messages.add("Xray_restraint_interaction", "Too little reflections. Set higher resolution!", io::message::error);
    return 1;
  }

#endif
  return 0;
}

void interaction::Electron_Density_Umbrella_Weight::increment_weight() {
#ifdef HAVE_CLIPPER
  // convert to angstrom
  const float cutoff2 = param.cutoff * param.cutoff * to_ang * to_ang;
  // create the range (size of atom)
  const clipper::Cell & cell = rho_calc.cell();
  const clipper::Grid_sampling & grid = rho_calc.grid_sampling();
  clipper::Grid_range gd(cell, grid, param.cutoff * to_ang);

  const int atoms_size = param.atoms.size();

  // loop over atoms to find covered grid points
  std::set<int> grid_points;

  for (int i = 0; i < atoms_size; i++) {
    DEBUG(10, "Atom index: " << param.atoms[i]);
    const clipper::Coord_orth & atom_pos = atoms[param.atoms[i]].coord_orth();

    // determine grad range of atom
    const clipper::Coord_frac uvw = atom_pos.coord_frac(cell);
    const clipper::Coord_grid g0 = uvw.coord_grid(grid) + gd.min();
    const clipper::Coord_grid g1 = uvw.coord_grid(grid) + gd.max();

    // loop over atom's grid
    clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;
    i0 = clipper::Xmap<clipper::ftype32>::Map_reference_coord(rho_obs, g0);
    for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
      for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
        for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
          // check if this cell is within the cutoff
          if ((iw.coord_orth() - atom_pos).lengthsq() < cutoff2) {
            grid_points.insert(iw.index());
          } // if within cutoff
        }
      }
    } // loop over grid
  } // loop over atoms
  // now we have a list of grid points without duplicates. Calculate density
  // differences (weight)

   // the fitting parameters
  double a, b;
  fit_rho(rho_calc, rho_obs, grid_points, a, b);
  double sum_obs_m_calc = 0.0;
  double sum_obs_p_calc = 0.0;
  DEBUG(10, "old weight: " << weight);
  // loop over grid points
  for (std::set<int>::const_iterator it = grid_points.begin(), to = grid_points.end();
          it != to; ++it) {
    // bring obs on the same scale as calc
    const double obs = a*rho_obs.get_data(*it)+b;
    DEBUG(12, "obs: " << obs);
    const double calc = rho_calc.get_data(*it);
    DEBUG(12, "calc: " << calc);

    // for R-real
    sum_obs_m_calc += fabs(obs - calc);
    sum_obs_p_calc += fabs(obs + calc);
  }
  const double r_real = sum_obs_m_calc / sum_obs_p_calc;
  DEBUG(10, "R real space: " << r_real);

  if (r_real > param.threshold) {
    const double term = r_real - param.threshold;
    weight += 0.5 * term * term;
  } 
  if (!param.threshold_freeze && r_real <= param.threshold) {
    // signal that the optimal value was found.
    param.signal = 1;
  }

  DEBUG(10, "new weight: " << weight);
#endif
}

void interaction::Electron_Density_Umbrella_Weight::write(std::ostream& os) const {
  os.precision(8);
  os.flags(std::ios::scientific);
  os << std::setw(15) << weight;
}

util::Umbrella_Weight * interaction::Electron_Density_Umbrella_Weight_Factory::get_instance() {
#ifdef HAVE_CLIPPER
  util::Umbrella_Weight * instance =
          new interaction::Electron_Density_Umbrella_Weight(param, conf, atoms, rho_calc, rho_obs, to_ang);
  // make sure the weight is increased as soon as the instance is created
  instance->increment_weight();
  instances.push_back(instance);
  return instance;
#endif
  return NULL;
}

