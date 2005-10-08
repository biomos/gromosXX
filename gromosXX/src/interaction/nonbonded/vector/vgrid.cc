/**
 * @file vgrid.cc
 * contains the grid nonbonded interaction implementation
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/nonbonded/interaction/nonbonded_parameter.h>
#include <interaction/nonbonded/pairlist/pairlist.h>

#include "vgrid.h"

#include <algorithm>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

namespace interaction
{

  typedef float cg_real_t;
  typedef double at_real_t;

  // grid properties
  cg_real_t box_x, box_y, box_z;
  at_real_t ex_box_x, ex_box_y, ex_box_z;

  cg_real_t dx, dy, dz;
  int Nx, Ny, Nz;
  int Cx, Cy, Cz;
  std::vector<int> cell_start;
  std::vector<int> cell_num;
  std::vector<int> mask;

  std::vector<int> pl;
  std::vector<std::vector<double> > lj_energy;
  std::vector<std::vector<double> > crf_energy;
  std::vector<std::vector<double> > lr_lj_energy;
  std::vector<std::vector<double> > lr_crf_energy;

  std::vector<cg_real_t> cg_shift_x, cg_shift_y, cg_shift_z;
  std::vector<at_real_t> at_shift_x, at_shift_y, at_shift_z;
  
  // chargegroup properties
  std::vector<int> cg_index, cg_tmp_index;
  std::vector<int> cg_ex_index;
  
  std::vector<cg_real_t> cg_x, cg_y, cg_z;
  std::vector<cg_real_t> cg_tmp_x, cg_tmp_y, cg_tmp_z;
  
  std::vector<int> cg_shift;
  std::vector<int> cg_tmp_shift;
  
  std::vector<int> cg_at_start;
  std::vector<int> cg_cell;
  std::vector<int> cg_tmp_cell;

  // atom properties
  std::vector<at_real_t> at_x, at_y, at_z;
  std::vector<int> at_iac;
  std::vector<at_real_t> at_q;
  std::vector<int> at_egroup;
  std::vector<at_real_t> at_fx, at_fy, at_fz;
  std::vector<at_real_t> at_lr_fx, at_lr_fy, at_lr_fz;

  // constants
  at_real_t cut3i;
  at_real_t crf;
  at_real_t crf_cut3i;
  at_real_t crf_2cut3i;
  at_real_t crf_cut;

  Nonbonded_Parameter * parameter;

  // calculation vectors have to be on the stack => OMP parallelization!

  // forward declarations
  void grid_properties(topology::Topology const & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation const & sim);

  void calc_mask(cg_real_t cutoff, std::vector<int> & mask);

  void grid_extend(topology::Topology const & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation const & sim);
  
  void grid_cg(topology::Topology const & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation const & sim);

  void grid_pl(topology::Topology const & topo,
	       std::vector<int> & pl, cg_real_t cutl, cg_real_t cuts,
	       bool print = false);

  void grid_print_pl(std::vector<int> const & pl);
  void grid_print_pl_atomic(topology::Topology const & topo,
			    std::vector<int> const & pl);
  
  void grid_interaction(std::vector<int> const & pl);
  
  void grid_init(simulation::Simulation const & sim);
  
  void grid_store_lr();
  
  void grid_restore_lr();
  
  void grid_store(topology::Topology const & topo,
		  configuration::Configuration & conf);

  // the fun begins
  void grid_prepare(topology::Topology const & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation const & sim,
		    Nonbonded_Parameter & param)
  {
    parameter = &param;
  
    if(!(sim.steps() % sim.param().pairlist.skip_step)){
      grid_init(sim);
      grid_properties(topo, conf, sim);
      grid_extend(topo, conf, sim);
      grid_cg(topo, conf, sim);
    }
  }
  
  void grid_update(topology::Topology const & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation const & sim)
  {
    if(!(sim.steps() % sim.param().pairlist.skip_step)){
      grid_pl(topo, 
	      pl, sim.param().pairlist.cutoff_long,
	      sim.param().pairlist.cutoff_short,
	      sim.param().pairlist.print);
      
      grid_store_lr();

      if (sim.param().pairlist.print) 
	grid_print_pl_atomic(topo, pl);
    }
    else{
      grid_restore_lr();
    }
    
    // short-range
    grid_interaction(pl);

    grid_store(topo, conf);

  }
  
  void grid_properties(topology::Topology const & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation const & sim)
  {
    cg_real_t cutoff = sim.param().pairlist.cutoff_long;

    math::Vec const & K = conf.current().box(0);
    math::Vec const & L = conf.current().box(1);
    math::Vec const & M = conf.current().box(2);

    box_x = cg_real_t(K(0));
    box_y = cg_real_t(L(1));
    box_z = cg_real_t(M(2));

    ex_box_x = K(0) + 2.0 * cutoff;
    ex_box_y = L(1) + 2.0 * cutoff;
    ex_box_z = M(2) + 2.0 * cutoff;

    // std::cerr << "ex box: (" << ex_box_x << ", " << ex_box_y
    // << ", " << ex_box_z << ")" << std::endl;

    cg_real_t ideal_length = sim.param().pairlist.grid_size;
    
    dx = box_x / rint(box_x / ideal_length);
    dy = box_y / rint(box_y / ideal_length);
    dz = box_z / rint(box_z / ideal_length);
  
    Cx = int((box_x) / dx + 1);
    Cy = int((box_y) / dy + 1);
    Cz = int((box_z) / dz + 1);

    Nx = int((box_x + 2.0 * cutoff) / dx + 1);
    Ny = int((box_y + 2.0 * cutoff) / dy + 1);
    Nz = int((box_z + 2.0 * cutoff) / dz + 1);

    cell_num.assign(Nx * Ny * Nz, 0);
    cell_start.resize(Nx * Ny * Nz + 1);

    // std::cerr << "grid: C(" << Cx << ", " << Cy << ", " << Cz << ")"
    // << " N(" << Nx << ", " << Ny << ", " << Nz << ")"
    // << " d(" << dx << ", " << dy << ", " << dz << ")"
    // << std::endl;

    cg_shift_x.resize(27);
    cg_shift_y.resize(27);
    cg_shift_z.resize(27);
    at_shift_x.resize(27);
    at_shift_y.resize(27);
    at_shift_z.resize(27);
    
    int c = 0;
    
    for(int k=-1; k<2; ++k){
      for(int l=-1; l<2; ++l){
	for(int m=-1; m<2; ++m){
	    
	  cg_shift_x[c] = cg_real_t(k*K(0) + l*L(0) + m*M(0));
	  cg_shift_y[c] = cg_real_t(k*K(1) + l*L(1) + m*M(1));
	  cg_shift_z[c] = cg_real_t(k*K(2) + l*L(2) + m*M(2));

	  at_shift_x[c] = at_real_t(k*K(0) + l*L(0) + m*M(0));
	  at_shift_y[c] = at_real_t(k*K(1) + l*L(1) + m*M(1));
	  at_shift_z[c] = at_real_t(k*K(2) + l*L(2) + m*M(2));
	    
	  ++c;
	}
      }
    }

    calc_mask(cutoff, mask);
  }

  void calc_mask(cg_real_t cutoff, std::vector<int> & mask)
  {
    cg_real_t distx2, disty2, distz2;
    cg_real_t cutoff2 = cutoff * cutoff;
    cg_real_t dx2 = dx*dx, dy2 = dy*dy, dz2 = dz*dz;
    int p = Ny * Nz;
    int min_y, max_y;
    int ay;
    
    mask.clear();

    for(int x=0; true; ++x){

      if (x>1)
	distx2 = (x-1)*(x-1) * dx2;
      else distx2 = 0.0;
    
      if (distx2 > cutoff2) break;
    
      // calculate max y
      max_y = int(sqrt((cutoff2 - distx2) / dy2)) + 1;
      if (x)
	min_y = -max_y;
      else
	min_y = 0;
      
      for(int y=min_y; y<=max_y; ++y){
	
	ay = abs(y);
	
	if (ay>1)
	  disty2 = (ay-1) * (ay-1) * dy2;
	else disty2 = 0.0;
	
	if (distx2 + disty2 > cutoff2) break;
	
	for(int z=1; true; ++z){
	  
	  if (z>1)
	    distz2 = (z-1)*(z-1) * dz2;
	  else distz2 = 0.0;
	  
	  if (distx2 + disty2 + distz2 > cutoff2){
	    
	    int beg, end;
	    if (x || y){
	      beg = x*p + y*Nz - z + 1;
	      end = x*p + y*Nz + z;
	    }
	    else{
	      beg = 0;
	      end = z;
	    }
	    mask.push_back(beg);
	    mask.push_back(end);
	    
	    /*
	      cout << "mask: " 
	      << setw(3) << x << " | " << setw(3) << y << " | " << setw(3) << z
	      << "    ( "
	      << setw(5) << beg << " -> " << setw(5) << end 
	      << " )" << endl;
	    */
	    
	    break;
	  }
	  
	} // z
      } // y
    } // x
  }

  // put into brick-wall rectangular box around
  // 0.5 * (box_x+2*cutoff, box_y+2*cutoff, box_z+2*cutoff)
  void grid_box(math::Vec & v, math::Box const & box)
  {
    math::Vec const & K = box(0);
    math::Vec const & L = box(1);
    math::Vec const & M = box(2);

    double center_x = 0.5 * ex_box_x,
      center_y = 0.5 * ex_box_y,
      center_z = 0.5 * ex_box_z;

    // std::cerr << "vx=" << v(0) << " vy=" << v(1) << " vz=" << v(2) << std::endl;
    // std::cerr << "box center " << center_x << " " << center_y << " " << center_z << std::endl;
    
    double dd = center_z - v(2);
    // std::cerr << "dz=" << dd << " Mz=" << M(2) << std::endl;
    if (fabs(dd) > 0.5 * M(2)){
      v(0) += rint(dd / M(2)) * M(0);
      v(1) += rint(dd / M(2)) * M(1);
      v(2) += rint(dd / M(2)) * M(2);
    }

    dd = center_y - v(1);
    // std::cerr << "dy=" << dd << " Ly=" << L(1) << std::endl;
    if (fabs(dd) > 0.5 * L(1)){
      v(0) += rint(dd / L(1)) * L(0);
      v(1) += rint(dd / L(1)) * L(1);
    }

    dd = center_x - v(0);
    // std::cerr << "dx=" << dd << " Kx=" << K(0) << std::endl;
    if (fabs(dd) > 0.5 * K(0)){
      v(0) += rint(dd / K(0)) * K(0);
    }
  }

  void grid_extend(topology::Topology const & topo,
		   configuration::Configuration & conf,
		   simulation::Simulation const & sim)
  {
    math::VArray &pos = conf.current().pos;
    math::Vec v, v_box, trans;

    cg_ex_index.clear();
    cg_tmp_index.clear();
    cg_tmp_x.clear();
    cg_tmp_y.clear();
    cg_tmp_z.clear();
    cg_tmp_shift.clear();
    cg_tmp_cell.clear();

    topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
      cg_to = topo.chargegroup_end();

    // solute chargegroups...
    int cg_index = 0;
    int new_cg_index = 0;
    int solute_cg = topo.num_solute_chargegroups();
    
    for( ; cg_it != cg_to; ++cg_it, ++cg_index){

      // std::cerr << "cg " << cg_index << std::endl;
      
      if (cg_index < solute_cg)
	cg_it.cog(pos, v);
      else
	v = pos(**cg_it);
      
      v_box = v;
      grid_box(v_box, conf.current().box);
      trans = v_box - v;
      
      // atoms in a chargegroup
      // std::cerr << "moving atoms" << std::endl;
      topology::Atom_Iterator at_it = cg_it.begin(),
	at_to = cg_it.end();
      for( ; at_it != at_to; ++at_it){
	assert(pos.size() > *at_it);
	pos(*at_it) += trans;
      } // loop over atoms

      // try extensions (bot not for k=-1)
      for(int s=0; s<27; ++s){
      
	cg_real_t sx = v_box(0) + cg_shift_x[s];
	cg_real_t sy = v_box(1) + cg_shift_y[s];
	cg_real_t sz = v_box(2) + cg_shift_z[s];

	// std::cerr << "cg[" << cg_index << "] x=" << v_box(0)
	// << " y=" << v_box(1) << " z=" << v_box(2)
	// << " sx=" << sx << " sy=" << sy << " sz=" << sz << "\t";

	// check if still inside
	if (sx >= 0 && sx <= ex_box_x &&
	    sy >= 0 && sy <= ex_box_y &&
	    sz >= 0 && sz <= ex_box_z){

	  // std::cerr << "adding s=" << s << std::endl;
	
	  // std::cerr << "\tshift " << s << " inside!" << std::endl;

	  cg_tmp_index.push_back(cg_index);
	  cg_ex_index.push_back(new_cg_index);
	  cg_tmp_x.push_back(sx);
	  cg_tmp_y.push_back(sy);
	  cg_tmp_z.push_back(sz);
	  cg_tmp_shift.push_back(s);
	  
	  // calculate grid cell
	  int grid_index =
	    int(sx / dx) * Ny * Nz +
	    int(sy / dy) * Nz +
	    int(sz / dz);
	  cg_tmp_cell.push_back(grid_index);
	  ++cell_num[grid_index];

	  // std::cerr << "\tgrid " << grid_index << std::endl;

	  ++new_cg_index;
	}
	else{
	  // std::cerr << "outside" << std::endl;
	}
      } // shifts
    } // loop over cg's

    int cs = 0;
    for(int i=0; i<Nx*Ny*Nz; ++i){
      cell_start[i] = cs;
      cs += cell_num[i];
    }
    cell_start[Nx*Ny*Nz] = cs;
  }

  struct cg_grid_less : std::binary_function<int, int, bool>
  {
    bool operator()(int const & cg1,
		    int const & cg2)const
    {
      if (cg_tmp_cell[cg1] < cg_tmp_cell[cg2])
	return true;
      if (cg_tmp_cell[cg1] > cg_tmp_cell[cg2])
	return false;
      
      return cg1 < cg2;
    }
  };

  void grid_cg(topology::Topology const & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation const & sim)
  {
    // std::cerr << "grid cg" << std::endl;
    std::sort(cg_ex_index.begin(), cg_ex_index.end(), cg_grid_less());

    size_t s = cg_ex_index.size();

    cg_index.resize(s);
    cg_x.resize(s);
    cg_y.resize(s);
    cg_z.resize(s);
    cg_shift.resize(s);
    cg_cell.resize(s);
    cg_at_start.resize(s+1);
    
    at_x.clear();
    at_y.clear();
    at_z.clear();
    at_iac.clear();
    at_q.clear();
    at_egroup.clear();

    int at_start = 0;

    for(size_t i=0; i<s; ++i){
      
      const int ind = cg_ex_index[i];
      
      cg_index[i] = cg_tmp_index[ind];

      cg_x[i] = cg_tmp_x[ind];
      cg_y[i] = cg_tmp_y[ind];
      cg_z[i] = cg_tmp_z[ind];

      cg_shift[i] = cg_tmp_shift[ind];
      cg_cell[i] = cg_tmp_cell[ind];

      cg_at_start[i] = at_start;
      
      // and atom lists
      int at_index = topo.chargegroup(cg_index[i]),
	at_to = topo.chargegroup(cg_index[i]+1);

      // std::cerr << "cg " << cg_index[i] 
      // << " atoms " << at_index << " - " << at_to - 1 
      // << "\tcell " << cg_cell[i] << std::endl;
      
      for( ; at_index < at_to; ++at_index, ++at_start){
	at_x.push_back(conf.current().pos(at_index)(0) 
		       + at_shift_x[cg_shift[i]]);
	at_y.push_back(conf.current().pos(at_index)(1)
		       + at_shift_y[cg_shift[i]]);
	at_z.push_back(conf.current().pos(at_index)(2)
		       + at_shift_z[cg_shift[i]]);
	
	at_iac.push_back(topo.iac(at_index));
	at_q.push_back(topo.charge(at_index));
	at_egroup.push_back(topo.atom_energy_group(at_index));
      }
    }

    cg_at_start[s] = at_start;

    size_t as = at_x.size();
    at_fx.assign(as, 0.0);
    at_fy.assign(as, 0.0);
    at_fz.assign(as, 0.0);

  }

  void grid_pl(topology::Topology const & topo, 
	       std::vector<int> & pl, 
	       cg_real_t cutl, cg_real_t cuts,
	       bool print)
  {
    // std::cerr << "grid pl" << std::endl;
    
    const int solute_cg = topo.num_solute_chargegroups();
    
    cg_real_t cutl2 = cutl * cutl;
    cg_real_t cuts2 = cuts * cuts;

    pl.clear();
    
    int mask_size = mask.size();

    std::vector<int> lr_pl;

    const int num_cg = cg_index.size();
    
    std::cout << "longrange pairlist" << std::endl;
    
    // loop over chargegroups in central computational box
    for(int cg1 = 0; cg1 < num_cg; ++cg1){

      if (cg_shift[cg1] != 13) continue;
      const int cell = cg_cell[cg1];
      // std::cerr << "central cg " << cg1 << std::endl;
      
      pl.push_back(cg1);
      int range_ind = pl.size();
      pl.push_back(0); // placeholder

      lr_pl.clear();
      lr_pl.push_back(cg1);
      lr_pl.push_back(0); // placeholder
	  
      bool sr_range = false, lr_range = false;
	    
      for(int m=0; m<mask_size; m+=2){
	
	int cg2 = cell_start[cell + mask[m]];
	int cg2_end = cell_start[cell + mask[m+1]];
	
	for( ; cg2 < cg2_end; ++cg2){
	  
	  cg_real_t dist2 = 0.0;
	  cg_real_t dd = cg_x[cg1] - cg_x[cg2];
	  dist2 += dd * dd;
	  dd = cg_y[cg1] - cg_y[cg2];
	  dist2 += dd * dd;
	  dd = cg_z[cg1] - cg_z[cg2];
	  dist2 += dd * dd;
	  
	  if (dist2 > cutl2){
	    if (lr_range){
	      lr_pl.push_back(cg_at_start[cg2]);
	      lr_range = false;
	    }
	    if (sr_range){
	      pl.push_back(cg_at_start[cg2]);
	      sr_range = false;
	    }
	    continue;
	  }

	  // check exclusions: SOLVENT ?
	  // std::cerr << "pl: cg1=" << cg1 << " (real " << cg_index[cg1] << ")"
	  // << " cg2="  << cg2 << " (real " << cg_index[cg2] << ")"
	  // << std::endl;

	  if (cg_index[cg1] > solute_cg || cg_index[cg2] > solute_cg){
	    if (cg_index[cg1] == cg_index[cg2]){
	      if (lr_range){
		lr_pl.push_back(cg_at_start[cg2]);
		lr_range = false;
	      }
	      if (sr_range){
		pl.push_back(cg_at_start[cg2]);
		sr_range = false;
	      }
	      continue;
	    }
	  }
	  else if (std::find(topo.chargegroup_exclusion(cg_index[cg1]).begin(), 
			     topo.chargegroup_exclusion(cg_index[cg1]).end(), cg_index[cg2])
		   != topo.chargegroup_exclusion(cg_index[cg1]).end()) {
	    // excluded
	    // std::cerr << "\texcluded" << std::endl;
	    
	    if (lr_range){
	      lr_pl.push_back(cg_at_start[cg2]);
	      lr_range = false;
	    }
	    if (sr_range){
	      pl.push_back(cg_at_start[cg2]);
	      sr_range = false;
	    }
	    continue;
	  }

	  // std::cerr << "\tnot excluded" << std::endl;
	  
	  if (dist2 > cuts2){
	    if (lr_range) continue;
	    if (sr_range){
	      pl.push_back(cg_at_start[cg2]);
	      sr_range = false;
	    }
	    lr_pl.push_back(cg_at_start[cg2]);
	    lr_range = true;
	  }
	  else{
	    if (sr_range) continue;
	    if (lr_range){
	      lr_pl.push_back(cg_at_start[cg2]);
	      lr_range = false;
	    }
	    pl.push_back(cg_at_start[cg2]);
	    sr_range = true;
	  }
	  
	} // cg2 range

	// don't forget to close the ranges
	if (sr_range){
	  pl.push_back(cg_at_start[cg2_end]);
	  sr_range = false;
	}
	if (lr_range){
	  lr_pl.push_back(cg_at_start[cg2_end]);
	  lr_range = false;
	}
      } // mask

      lr_pl[1] = (lr_pl.size() - 2) / 2;
      // grid_print_pl(lr_pl);
      grid_interaction(lr_pl);
      
      pl[range_ind] = (pl.size() - range_ind - 1) / 2;
      
    } // central cg's (cg1)
  } // cg_pl()

  
  void grid_print_pl(std::vector<int> const & pl)
  {
    size_t i = 0;
    while(i < pl.size() - 1){

      std::cout << std::setw(4) << pl[i] << " : ";
      size_t i_to = i + pl[i+1] * 2 + 2;
      i += 2;

      for( ; i < i_to; i += 2){
	std::cout << std::setw(4) << pl[i] << " - " << std::setw(4) << pl[i+1] << "  ";
      }
      std::cout << "\n";
    }
  }

  void grid_init(simulation::Simulation const & sim)
  {
    cut3i = 1.0 / ( sim.param().longrange.rf_cutoff
		    * sim.param().longrange.rf_cutoff
		    * sim.param().longrange.rf_cutoff);
    
    crf = 2*(sim.param().longrange.epsilon - sim.param().longrange.rf_epsilon) * 
      (1.0 + sim.param().longrange.rf_kappa * sim.param().longrange.rf_cutoff) -
      sim.param().longrange.rf_epsilon * (sim.param().longrange.rf_kappa  * 
					  sim.param().longrange.rf_cutoff *
					  sim.param().longrange.rf_kappa  *
					  sim.param().longrange.rf_cutoff);
    
    crf /= (sim.param().longrange.epsilon +2* sim.param().longrange.rf_epsilon) *
      (1.0 + sim.param().longrange.rf_kappa * sim.param().longrange.rf_cutoff) +
      sim.param().longrange.rf_epsilon * (sim.param().longrange.rf_kappa  * 
					  sim.param().longrange.rf_cutoff *
					  sim.param().longrange.rf_kappa  *
					  sim.param().longrange.rf_cutoff);
    crf_cut3i = crf * cut3i;
    
    crf_2cut3i = crf_cut3i / 2.0;
    
    crf_cut = (1 - crf / 2.0) / sim.param().longrange.rf_cutoff;

    const int egroups = sim.param().force.energy_group.size();
    
    lj_energy.resize(egroups);
    crf_energy.resize(egroups);
    lr_lj_energy.resize(egroups);
    lr_crf_energy.resize(egroups);
    
    for(int i=0; i < egroups; ++i){
      lj_energy[i].assign(egroups, 0.0);
      crf_energy[i].assign(egroups, 0.0);
      lr_lj_energy[i].assign(egroups, 0.0);
      lr_crf_energy[i].assign(egroups, 0.0);
    }
  }

  at_real_t c6[100], c12[100], q[100];
  at_real_t rx[100], ry[100], rz[100];
  at_real_t r2[100], ir2[100], ir[100], ir6[100];
  at_real_t q_eps[100], c12_ir6[100];
  at_real_t e_lj[100], e_crf[100], f[100];
  
  void grid_interaction(std::vector<int> const & pl)
  {
    std::vector<std::vector<lj_parameter_struct> > const & lj_param
      = parameter->lj_parameter();

    size_t i = 0;
    while(i < pl.size() - 1){
      
      int cg1 = pl[i];
      size_t i_to = i + pl[i+1] * 2 + 2;
      i += 2;

      // ranges
      for( ; i < i_to; i += 2){
	// at1
	int at1 = cg_at_start[cg1], at1_to = cg_at_start[cg1+1];
	for( ; at1 != at1_to; ++at1){
	  // at2
	  int at2 = pl[i], at2_to = pl[i+1];
	  int num = at2_to - at2;

	  // prepare c6, c12, rx, ry, rz

	  // not vectorized: dependencies (nr?), still?
	  for(int nr = 0; nr < num; ++nr){
	    c6[nr] = lj_param[at_iac[at1]][at_iac[nr + at2]].c6;
	    c12[nr] = lj_param[at_iac[at1]][at_iac[nr + at2]].c12;
	    q[nr] = at_q[at1] * at_q[nr + at2];
	    rx[nr] = at_x[at1] - at_x[nr + at2];
	    ry[nr] = at_y[at1] - at_y[nr + at2];
	    rz[nr] = at_z[at1] - at_z[nr + at2];
	  } // at2
	  // calc r2
	  for(int nr = 0; nr < num; ++nr){
	    r2[nr] = rx[nr] * rx[nr] + ry[nr] * ry[nr] + rz[nr]*rz[nr];
	  }
	  // calc 1/sqrt(r2)
	  for(int nr = 0; nr < num; ++nr){
	    ir2[nr] = 1.0 / r2[nr];
	  }
	  for(int nr = 0; nr < num; ++nr){
	    ir[nr] = sqrt(ir2[nr]);
	  }
	  for(int nr = 0; nr < num; ++nr){
	    ir6[nr] = ir2[nr] * ir2[nr] * ir2[nr];
	  }
	  for(int nr = 0; nr < num; ++nr){
	    c12_ir6[nr] = c12[nr] * ir6[nr];
	  }
	  for(int nr = 0; nr < num; ++nr){
	    e_lj[nr] = (c12_ir6[nr] - c6[nr]) * ir6[nr];
	  }
	  for(int nr = 0; nr < num; ++nr){
	    q_eps[nr] = math::four_pi_eps_i * q[nr];
	  }
	  for(int nr = 0; nr < num; ++nr){
	    e_crf[nr] = q_eps[nr] * (ir[nr] - crf_2cut3i * r2[nr] - crf_cut);
	  }
	  for(int nr = 0; nr < num; ++nr){
	    f[nr] = (c12_ir6[nr] + c12_ir6[nr] - c6[nr]) * 6.0 * ir6[nr] * ir2[nr] + 
	      q_eps[nr] * (ir[nr] * ir2[nr] + crf_cut3i);
	  }

	  // dependencies ? (nr?), still?
	  for(int nr = 0; nr < num; ++nr){
	    at_fx[at1] += f[nr] * rx[nr];
	    at_fx[nr + at2] -= f[nr] * rx[nr];
	    at_fy[at1] += f[nr] * ry[nr];
	    at_fy[nr + at2] -= f[nr] * ry[nr];
	    at_fz[at1] += f[nr] * rz[nr];
	    at_fz[nr + at2] -= f[nr] * rz[nr];
	    
	    lj_energy[at_egroup[at1]][at_egroup[nr+at2]] += e_lj[nr];
	    crf_energy[at_egroup[at1]][at_egroup[nr+at2]] += e_crf[nr];
	    
	  } // at2

	  // store the energies...

	} // at1
      } // range
    } // pairlist cg1

  }

  void grid_store_lr()
  {
    // store long-range forces
    at_lr_fx = at_fx;
    at_lr_fy = at_fy;
    at_lr_fz = at_fz;
    
    // and energies
    const int egroups = lj_energy.size();
    for(int i=0; i<egroups; ++i){
      for(int j=0; j<egroups; ++j){
	lr_lj_energy[i][j] = lj_energy[i][j];
	lr_crf_energy[i][j] = crf_energy[i][j];
      }
    }
  }

  void grid_restore_lr()
  {
    // set forces to the long-range forces
    at_fx = at_lr_fx;
    at_fy = at_lr_fy;
    at_fz = at_lr_fz;

    // and energies
    const int egroups = lj_energy.size();
    for(int i=0; i<egroups; ++i){
      for(int j=0; j<egroups; ++j){
	lj_energy[i][j] = lr_lj_energy[i][j];
	crf_energy[i][j] = lr_crf_energy[i][j];
      }
    }

  }

  void grid_store(topology::Topology const & topo,
		  configuration::Configuration & conf)
  {
    // loop over all cg's
    for(size_t cg = 0; cg < cg_index.size(); ++cg){

      int sat = topo.chargegroup(cg_index[cg]);
      
      // unsupported loop structure
      for(int at=cg_at_start[cg]; at<cg_at_start[cg+1]; ++at, ++sat){

	conf.current().force(sat)(0) += at_fx[at];
	conf.current().force(sat)(1) += at_fy[at];
	conf.current().force(sat)(2) += at_fz[at];
	
	

      }
    }

    const int ljs = conf.current().energies.lj_energy.size();
    configuration::Energy & e = conf.current().energies;

    for(int i = 0; i < ljs; ++i){
      for(int j = 0; j < ljs; ++j){
	
	e.lj_energy[i][j] += lj_energy[i][j];
	e.crf_energy[i][j] += crf_energy[i][j];
      }
    }
    
  }

  void grid_print_pl_atomic(topology::Topology const & topo,
			    std::vector<int> const & pl)
  {
    std::map<int, int> at_map;

    int grid_at = 0;
    
    for(size_t cg = 0; cg < cg_index.size(); ++cg){
      for(int at = topo.chargegroup(cg_index[cg]); at < topo.chargegroup(cg_index[cg]+1);
	  ++at, ++grid_at){
	at_map[grid_at] = at;
      }
    }

    Pairlist tmp_pl;
    tmp_pl.resize(topo.num_atoms());
    
    { // create temporary atomic pairlist
      size_t i = 0;
      while(i < pl.size() - 1){
	
	int cg1 = pl[i];
	
	size_t i_to = i + pl[i+1] * 2 + 2;
	i += 2;
	
	for( ; i < i_to; i += 2){
	  
	  int at2_beg = pl[i], at2_end = pl[i+1];
	  
	  for(int at1 = cg_at_start[cg1]; at1 < cg_at_start[cg1+1]; ++at1){
	    for(int at2 = at2_beg; at2 < at2_end; ++at2){
	      
	      int at1_real = at_map[at1];
	      int at2_real = at_map[at2];
	      
	      if (at1_real < at2_real) tmp_pl[at1_real].push_back(at2_real);
	      else tmp_pl[at2_real].push_back(at1_real);
	    }
	  }
	}
      }
    }
    
    for(size_t i=0; i<tmp_pl.size(); ++i)
      std::sort(tmp_pl[i].begin(), tmp_pl[i].end());

    std::cout << tmp_pl << std::endl;

  }
  
} // interaction
