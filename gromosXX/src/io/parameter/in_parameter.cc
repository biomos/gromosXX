
/**
 * @file in_parameter.cc
 * implemes
 * nts methods of In_Parameter
 */

#include <stdheader.h>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>

#include <simulation/multibath.h>
#include <simulation/parameter.h>

#include <io/instream.h>
#include <io/blockinput.h>

#include <math/random.h>

#include "in_parameter.h"

#ifdef OMP
#include <omp.h>
#endif

#undef MODULE
#define MODULE io
#undef SUBMODULE
#define SUBMODULE parameter

static std::set<std::string> block_read;

/**
 * Store standard parameters in the Parameter
 */
void io::In_Parameter::read(simulation::Parameter &param,
			    std::ostream & os)
{
  DEBUG(7, "reading input");

  if (!quiet)
    os << "\nINPUT\n"
	      << title << "\n";

  // store the title...
  param.title = title;

  read_ENERGYMIN(param);
  read_SYSTEM(param);
  read_COMTRANSROT(param); // has to be read before INITIALISE
  read_INITIALISE(param);
  read_STEP(param);
  read_BOUNDCOND(param);
  read_REPLICA(param); // has to be read in before MULTIBATH
  read_MULTIBATH(param);
  read_PRESSURESCALE(param);
  read_PRINTOUT(param);
  read_WRITETRAJ(param);
  read_CONSTRAINT(param);
  read_FORCE(param); 
  read_COVALENTFORM(param);
  read_CGRAIN(param);
  read_PAIRLIST(param);
  read_LONGRANGE(param);
  read_POSITIONRES(param);
  read_DISTANCERES(param);
  read_DIHEDRALRES(param); // needs to be called after CONSTRAINT!
  read_PERTURBATION(param);
  read_JVALUERES(param);
  read_PERSCALE(param);
  read_ROTTRANS(param);
  read_INNERLOOP(param);
  read_MULTICELL(param);
  read_READTRAJ(param);
  read_INTEGRATE(param);
  read_STOCHDYN(param);
  read_EWARN(param);
  read_MULTISTEP(param);
  read_MONTECARLO(param);
  read_RAMD(param);
  read_POLARIZE(param);
  read_RANDOMNUMBERS(param);
  read_EDS(param);
  read_LAMBDAS(param);
  read_known_unsupported_blocks();
  
  DEBUG(7, "input read...");

  for(std::map<std::string, std::vector<std::string> >::const_iterator
	it = m_block.begin(),
	to = m_block.end();
      it != to;
      ++it){

    if (block_read.count(it->first) == 0 && it->second.size()){
      io::messages.add("block " + it->first + " not supported!",
		       "In_Parameter",
		       io::message::error);
    }
  }

  if (!quiet)
    os << "END\n";

}

/**
 * read the SYSTEM block.
 */
void io::In_Parameter::read_SYSTEM(simulation::Parameter &param,
				   std::ostream & os)
{
  DEBUG(8, "reading SYSTEM");

  std::vector<std::string> buffer;
  buffer = m_block["SYSTEM"];
  std::string s;
  
  if (!buffer.size()){
    io::messages.add("no SYSTEM block in input", "In_Parameter", io::message::error);
    param.system.nsm = 0;
    param.system.npm = 0;
    return;
  }

  block_read.insert("SYSTEM");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  _lineStream >> param.system.npm 
	      >> param.system.nsm;
 
  if (_lineStream.fail())
    io::messages.add("bad line in SYSTEM block",
		       "In_Parameter", io::message::error);

  /**
  if (!_lineStream.eof())
    io::messages.add("End of line not reached in SYSTEM block, but should have been: \n" + _lineStream.str() +  "\n",
		       "In_Parameter", io::message::warning);
  */
  
  // we might need to also allow for 0...
  if (param.system.npm != 1 && param.system.npm != 0)
    io::messages.add("SYSTEM: currently only NPM=1 allowed (NPM=0 experimental)",
		     "io::In_Parameter::read_SYSTEM",
		     io::message::error);
  if(param.system.nsm < 0)
    io::messages.add("SYSTEM: NSM should be >0",
		     "io::In_Parameter::read_SYSTEM",
		     io::message::error);

} 

/**
 * read the ENERGYMIN block.
 */
void io::In_Parameter::read_ENERGYMIN(simulation::Parameter &param,
				     std::ostream & os)
{
  DEBUG(8, "reading ENERGYMIN");

  std::vector<std::string> buffer;
  buffer = m_block["ENERGYMIN"];
  std::string s;
  
  if (!buffer.size()){
    // no minimisation
    return;
  }

  buffer = m_block["ENERGYMIN"];
  
  block_read.insert("ENERGYMIN");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  _lineStream >> param.minimise.ntem 
	      >> param.minimise.ncyc
	      >> param.minimise.dele
	      >> param.minimise.dx0
	      >> param.minimise.dxm
              >> param.minimise.nmin
              >> param.minimise.flim;
  
  if (_lineStream.fail())
    io::messages.add("bad line in ENERGYMIN block",
		       "In_Parameter", io::message::error);
 
  // allow 0 to disable feature...
  if (param.minimise.nmin == 0)
    param.minimise.nmin = 1;

  if (param.minimise.ntem != 0 && param.minimise.ntem != 1)
    io::messages.add("ENERGYMIN block: currently only steepest descent implemented",
		     "io::In_Parameter",
		     io::message::error);
  if(param.minimise.ntem == 1 && param.minimise.ncyc > 0)
    io::messages.add("ENERGYMIN block: NCYC > 0 has no effect for steepest descent",
		     "io::In_Parameter",
		     io::message::warning);
    
  if(param.minimise.ncyc < 0)
    io::messages.add("ENERGYMIN block: NCYC should be >0",
		     "io::In_Parameter",
		     io::message::error);
    
  if(param.minimise.dele < 0)
    io::messages.add("ENERGYMIN block: DELE should be >0",
		     "io::In_Parameter",
		     io::message::error);
    
  if(param.minimise.dx0 < 0)
    io::messages.add("ENERGYMIN block: DX0 should be >0",
		     "io::In_Parameter", io::message::error);
    
  if(param.minimise.dxm < param.minimise.dx0)
    io::messages.add("ENERGYMIN block: DXM should be > DX0",
		     "io::In_Parameter", io::message::error);

  if(param.minimise.nmin <= 0)
    io::messages.add("ENERGYMIN block: NMIN should be >= 0",
		     "io::In_Parameter", io::message::error);
  
  if(param.minimise.flim < 0)
    io::messages.add("ENERGYMIN: FLIM should be >= 0",
		     "io::In_Parameter",
		     io::message::error);

} 

/**
 * read the STEP block.
 */
void io::In_Parameter::read_STEP(simulation::Parameter &param,
				 std::ostream & os)
{
  DEBUG(8, "read STEP");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["STEP"];
  
  if (!buffer.size()){
    io::messages.add("no STEP block in input", "In_Parameter", io::message::error);
    return;
  }

  block_read.insert("STEP");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  _lineStream >> param.step.number_of_steps
	      >> param.step.t0
	      >> param.step.dt;
  
  if (_lineStream.fail())
    io::messages.add("bad line in STEP block", "In_Parameter", io::message::error);

  if(param.step.t0 < 0 && param.step.t0 != -1.0)
    io::messages.add("STEP block: Negative time is not supported",
		     "In_Parameter", io::message::error);
  if(param.step.number_of_steps <= 0)
    io::messages.add("STEP block: We want to do at least one step...",
		     "In_Parameter", io::message::error);
}

/**
 * the CONSTRAINT block.
 */
void io::In_Parameter::read_CONSTRAINT(simulation::Parameter &param,
					std::ostream & os)
{
  DEBUG(8, "read CONSTRAINT");

  std::vector<std::string> buffer;
  std::string s, salg;
  
  buffer = m_block["CONSTRAINT"];

  if (!buffer.size()){
    param.constraint.ntc = 1;
    param.constraint.solute.algorithm = simulation::constr_off;
    param.constraint.solvent.algorithm = simulation::constr_shake;

    io::messages.add("no CONSTRAINT block", "In_Parameter",
		     io::message::warning);
    return;
  }
  
  block_read.insert("CONSTRAINT");
  
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  std::string sntc;
  _lineStream >> sntc;

  if(sntc=="solvent") param.constraint.ntc=1;
  else if(sntc=="hydrogen") param.constraint.ntc=2;
  else if(sntc=="all") param.constraint.ntc=3;
  else if(sntc=="specified") param.constraint.ntc=4;
  else {
    std::stringstream ss(sntc);
    if (!(ss >> param.constraint.ntc))
      io::messages.add("CONSTRAINT block: NTC not understood",
		       "In_Parameter", io::message::error);
  }
  
  if(param.constraint.ntc<1 || param.constraint.ntc > 4){
    io::messages.add("CONSTRAINT block: NTC out of rangek",
		     "In_Parameter", io::message::error);
    
    param.constraint.solute.algorithm = simulation::constr_off;
    param.constraint.solvent.algorithm = simulation::constr_off;
    
    return;
  }
  
  // SOLUTE
  _lineStream >> salg;
  
  DEBUG(8, "Constraints (solute): " << salg);

  if (_lineStream.fail())
    io::messages.add("bad line in CONSTRAINT block",
		     "In_Parameter", io::message::error);

  std::transform(salg.begin(), salg.end(), salg.begin(), tolower);

  if (salg == "shake"){
    DEBUG(9, "constraints solute shake");
    
    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = simulation::constr_shake;
    else param.constraint.solute.algorithm = simulation::constr_off;
    
    _lineStream >> param.constraint.solute.shake_tolerance;
    
    if(param.constraint.solute.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: shake tolerance should be > 0",
		       "In_Parameter", io::message::error);
  } else if (salg == "flexshake"){
    DEBUG(9, "constraints solute flexshake");
    
    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = simulation::constr_flexshake;
    else param.constraint.solute.algorithm = simulation::constr_off;
    
    _lineStream >> param.constraint.solute.shake_tolerance
		>> param.constraint.solute.flexshake_readin
		>> param.constraint.solute.flexshake_mode;
    
    if(param.constraint.solute.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: shake tolerance should be > 0.0",
		       "In_Parameter", io::message::error);

    if(param.constraint.solute.flexshake_mode < 0 ||
       param.constraint.solute.flexshake_mode > 3)
      io::messages.add("CONSTRAINT block: flexshake mode should be >= 0 and <= 3",
		       "In_Parameter", io::message::error);
    
  } else if (salg == "lincs"){
    DEBUG(9, "constraints solute lincs");

    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = simulation::constr_lincs;
    else 
      param.constraint.solute.algorithm = simulation::constr_off;
    
    _lineStream >> param.constraint.solute.lincs_order;
    
    if(param.constraint.solute.lincs_order < 1)
      io::messages.add("CONSTRAINT block: lincs order should be > 1",
		       "In_Parameter", io::message::error);

  } else if (salg == "off"){
    DEBUG(9, "constraints solute off");
    param.constraint.solute.algorithm = simulation::constr_off;
  } else{
    DEBUG(9, "constraints solute error");
    
    io::messages.add("CONSTRAINT block: unknown algorithm (solute)",
		     "In_Parameter", io::message::error);
    param.constraint.solute.algorithm = simulation::constr_off;
  }

  // SOLVENT
  _lineStream >> salg;
  
  DEBUG(8, "constraints solvent: " << salg);

  if (_lineStream.fail())
    io::messages.add("bad line in CONSTRAINT block",
		     "In_Parameter", io::message::error);

  std::transform(salg.begin(), salg.end(), salg.begin(), tolower);

  if (salg == "shake") {
    DEBUG(9, "constraints solvent shake");

    param.constraint.solvent.algorithm = simulation::constr_shake;
    _lineStream >> param.constraint.solvent.shake_tolerance;
    
    if(param.constraint.solvent.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: shake tolerance should be > 0.0",
		       "In_Parameter", io::message::error);
  } else if (salg == "flexshake"){
    DEBUG(9, "constraints solvent flexshake");

    param.constraint.solvent.algorithm = simulation::constr_flexshake;
    _lineStream >> param.constraint.solvent.shake_tolerance;
    
    if(param.constraint.solvent.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: shake tolerance should be > 0.0",
		       "In_Parameter", io::message::error);
  } else if (salg == "lincs"){
    DEBUG(9, "constraints solvent lincs");

    param.constraint.solvent.algorithm = simulation::constr_lincs;
    _lineStream >> param.constraint.solvent.lincs_order;
    
    if(param.constraint.solvent.lincs_order < 1)
      io::messages.add("CONSTRAINT block: lincs order should be >1",
		       "In_Parameter", io::message::error);

  } else if (salg == "off") {
    DEBUG(9, "constraints solvent off");

    param.constraint.solvent.algorithm = simulation::constr_off;
    io::messages.add("CONSTRAINT block: no constraints for SOLVENT",
		     "In_Parameter", io::message::warning);
  } else{

    DEBUG(9, "constraints solvent error");
    io::messages.add("CONSTRAINT block: unknown algorithm (solvent)",
		     "In_Parameter", io::message::error);
    
    param.constraint.solvent.algorithm = simulation::constr_off;
  }
}

/**
 * read the PRINTOUT
 */
void io::In_Parameter::read_PRINTOUT(simulation::Parameter &param,
				  std::ostream & os)
{
  DEBUG(8, "read PRINTOUT");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["PRINTOUT"];

  if (!buffer.size()){
    io::messages.add("no PRINTOUT block", "In_Parameter", io::message::notice);
    return;
  }

  block_read.insert("PRINTOUT");

  int ntpp;
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  _lineStream >> param.print.stepblock
	      >> ntpp;
  
  if (_lineStream.fail())
    io::messages.add("bad line in PRINTOUT block",
		     "In_Parameter", io::message::error);
  
  if(param.print.stepblock < 0)
    io::messages.add("PRINTOUT block: NTPR should be >=0.",
		     "In_Parameter", io::message::error);
  
  switch(ntpp) {
    case 0 : 
      param.print.monitor_dihedrals = false;
      break;
    case 1 :
      param.print.monitor_dihedrals = true;
      break;
    default:
    io::messages.add("PRINTOUT block: NTPP should be 0 or 1.",
		     "In_Parameter", io::message::error);      
  }
}

/**
 * read the WRITETRAJ
 */
void io::In_Parameter::read_WRITETRAJ(simulation::Parameter &param,
				  std::ostream & os)
{
  DEBUG(8, "read WRITETRAJ");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["WRITETRAJ"];

  if (!buffer.size()){
    io::messages.add("no WRITETRAJ block", "In_Parameter", io::message::notice);
    return;
  }

  block_read.insert("WRITETRAJ");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  int ntwse;
  
  _lineStream >> param.write.position
	      >> ntwse
	      >> param.write.velocity
              >> param.write.force
	      >> param.write.energy
	      >> param.write.free_energy
	      >> param.write.block_average;
  
  if (_lineStream.fail())
    io::messages.add("bad line in WRITETRAJ block",
		     "In_Parameter", io::message::error);

  if(ntwse!=0)
    io::messages.add("WRITETRAJ block: NTWSE != 0 not supported",
		     "In_Parameter", io::message::error);

  if(param.write.position < 0){
    param.write.position_solute_only = true;
    param.write.position = -param.write.position;
    io::messages.add("writing solute only position trajectory",
		     "In_Parameter", io::message::notice);
  }
  
  if(param.write.velocity < 0){
    param.write.velocity_solute_only = true;
    param.write.velocity = -param.write.velocity;
    io::messages.add("writing solute only velocity trajectory",
		     "In_Parameter", io::message::notice);
  }
  
  if(param.write.force < 0){
    param.write.force_solute_only = true;
    param.write.force = -param.write.force;
    io::messages.add("writing solute only force trajectory",
		     "In_Parameter", io::message::notice);
  }
  
  if(param.write.energy < 0)
    io::messages.add("WRITETRAJ block: NTWE should be >= 0",
		     "In_Parameter", io::message::error);
  if(param.write.free_energy < 0)
    io::messages.add("WRITETRAJ block: NTWG should be >= 0",
		     "In_Parameter", io::message::error);
  if(param.write.block_average < 0)
    io::messages.add("WRITETRAJ block: NTWB should be >= 0",
		     "In_Parameter", io::message::error);
}

/**
 * read the PRESSURESCALE block.
 */
void io::In_Parameter::read_PRESSURESCALE(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read PRESSURESCALE");

  std::vector<std::string> buffer;
  std::string s;
  

  // first try for a PRESSURESCALE block
  buffer = m_block["PRESSURESCALE"];

  if (buffer.size()){

    block_read.insert("PRESSURESCALE");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
    std::string s1, s2, s3;
    
    _lineStream >> s1 >> s2 
		>> param.pcouple.compressibility 
		>> param.pcouple.tau 
		>> s3;

    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j){
	_lineStream >> param.pcouple.pres0(i, j);
      }
    }

    if (_lineStream.fail())
      io::messages.add("bad line in PRESSURESCALE block",
		       "In_Parameter", io::message::error);
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    std::transform(s3.begin(), s3.end(), s3.begin(), tolower);

    if (s1 == "off"){
      param.pcouple.calculate = false;
      param.pcouple.scale = math::pcouple_off;
    } 
    else if (s1 == "calc"){
      param.pcouple.calculate = true;
      param.pcouple.scale = math::pcouple_off;
    }
    else if (s1 == "scale"){
      param.pcouple.calculate = true;

      if (s2 == "off"){
	io::messages.add("PRESSURESCALE block: requesting scaling but SCALE set to OFF",
			 "In_Parameter", io::message::error);
	param.pcouple.scale = math::pcouple_off;
      }
      else if (s2 == "iso")
	param.pcouple.scale = math::pcouple_isotropic;
      else if (s2 == "aniso")
	param.pcouple.scale = math::pcouple_anisotropic;
      else if (s2 == "full")
	param.pcouple.scale = math::pcouple_full_anisotropic;
      else{
	io::messages.add("PRESSURESCALE block: bad value for SCALE switch "
			 "(off,iso,aniso,full)",
			 "In_Parameter", io::message::error);
	param.pcouple.scale = math::pcouple_off;
      }

    }
    else{
      io::messages.add("bad value for calc switch in PRESSURESCALE block\n"
		       "(off,calc,scale)",
		       "In_Parameter", io::message::error);
      param.pcouple.calculate = false;
    }
  
    if (param.pcouple.calculate){
      if (s3 == "none"){
	io::messages.add("requesting pressure calculation but "
			 "no virial specified",
			 "In_Parameter", io::message::error);
	param.pcouple.virial = math::no_virial;
      }
      else if (s3 == "atomic")
	param.pcouple.virial = math::atomic_virial;
      else if (s3 == "molecular")
	param.pcouple.virial = math::molecular_virial;
      else{
	io::messages.add("bad value for virial switch in PRESSURESCALE block\n"
			 "(none,atomic,molecular)",
			 "In_Parameter", io::message::error);
	param.pcouple.virial = math::no_virial;
      }
    }
    else
      param.pcouple.virial = math::no_virial;
    
  } // PRESSURESCALE block
  
  if (param.pcouple.calculate==false && param.pcouple.scale!=math::pcouple_off)
    io::messages.add("PRESSURESCALE block: pressure coupling activated but "
		     "not calculating pressure",
		     "In_Parameter",
		     io::message::error);
  if (param.pcouple.calculate == true && param.pcouple.virial == math::no_virial)
    io::messages.add("PRESSURESCALE block: pressure calculation requested but"
		     " no virial specified!", "In_Parameter",
		     io::message::error);
  if(param.pcouple.compressibility <=0)
    io::messages.add("PRESSURESCALE block: compressibility should be >0 ",
		     "In_Parameter", io::message::error);
  if(param.pcouple.tau <=0)
    io::messages.add("PRESSURESCALE block: tau should be >0 ",
		     "In_Parameter", io::message::error);
}

/**
 * read the BOUNDCOND block.
 */
void io::In_Parameter::read_BOUNDCOND(simulation::Parameter &param,
				     std::ostream & os)
{
  DEBUG(8, "read BOUNDCOND");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["BOUNDCOND"];
  
  if (!buffer.size()){
    io::messages.add("no BOUNDCOND block", "In_Parameter", io::message::error);
    param.boundary.boundary = math::vacuum;
    return;
  }

  block_read.insert("BOUNDCOND");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  int ntb;
  _lineStream >> ntb >> param.boundary.dof_to_subtract;
  
  if (_lineStream.fail())
    io::messages.add("bad line in BOUNDCOND block",
		     "In_Parameter", io::message::error);

  if(ntb==0) param.boundary.boundary=math::vacuum;
  else if(ntb==1) param.boundary.boundary=math::rectangular;
  else if(ntb==2) param.boundary.boundary=math::triclinic;
  else if(ntb==-1) param.boundary.boundary=math::truncoct;
  else {
    std::ostringstream msg;
    msg << "BOUNDCOND block: wrong value for NTB "
        << ntb << "\nvacuum (0), rectangular (1), triclinic (2), truncoct (-1)";
    io::messages.add(msg.str(), "In_Parameter", io::message::error);
    param.boundary.boundary=math::vacuum;
  }
  
  if (param.boundary.dof_to_subtract < 0) {
    io::messages.add("BOUNDCOND block: NDFMIN must be >= 0.",
		     "In_Parameter", io::message::error);
    param.boundary.dof_to_subtract = 0;    
  }
  
  if (param.boundary.dof_to_subtract > 0) {
    io::messages.add("BOUNDCOND block: NDFMIN > 0 not implemented",
		     "In_Parameter", io::message::warning);
    param.boundary.dof_to_subtract = 0;    
  }
}

/**
 * read the PERTURBATION block.
 */
void io::In_Parameter::read_PERTURBATION(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read PERTURBATION");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["PERTURBATION"];
  if (!buffer.size()) 
    return;

  block_read.insert("PERTURBATION");


  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  std::string b, s1, s2;
  int ntg, scale;
  int nrdgl;

  _lineStream >> ntg 
              >> nrdgl
              >> param.perturbation.lambda
              >> param.perturbation.dlamt
              >> param.perturbation.soft_vdw
              >> param.perturbation.soft_crf
              >> param.perturbation.lambda_exponent
              >> scale;
    
  if (_lineStream.fail())
    io::messages.add("bad line in PERTURBATION block",
                     "In_Parameter", io::message::error);
    
  switch(ntg) {
    case 0 :
      param.perturbation.perturbation = false;
      break;
    case 1 :
      param.perturbation.perturbation = true;
      break;
    default:
      io::messages.add("PERTURBATION block: NTG must be 0 or 1.",
                       "In_Parameter", io::message::error);
  }
  
  switch(nrdgl) {
    case 0 : // use from input file
      param.perturbation.read_initial = false;
      break;
    case 1 : // use from configuration
      param.perturbation.read_initial = true;
      break;
    default :
      io::messages.add("PERTURBATION block: NRDGL must be 0 or 1.",
                       "In_Parameter", io::message::error);
  }
  
  if (param.perturbation.read_initial) {
    io::messages.add("PERTURBATION block: NRDGL != 0 not implemented.",
                     "In_Parameter", io::message::error);    
  }
  
  if (param.perturbation.lambda < 0.0 ||
      param.perturbation.lambda > 1.0) {
    io::messages.add("PERTURBATION block: RLAM must be 0.0 to 1.0.",
                     "In_Parameter", io::message::error);        
  }
  
  if (param.perturbation.dlamt < 0.0){
    io::messages.add("PERTURBATION block: DLAMT must be >= 0.",
                     "In_Parameter", io::message::error);
  }
  
  if (param.perturbation.lambda_exponent<=0){
    io::messages.add("PERTURBATION block: NLAM must be > 0.",
                     "In_Parameter", io::message::error);
  }
  
  if (param.perturbation.soft_vdw < 0.0){
    io::messages.add("PERTURBATION block: ALPHLJ must be >= 0.",
                     "In_Parameter", io::message::error);
  }
  
  if (param.perturbation.soft_vdw < 0.0){
    io::messages.add("PERTURBATION block: ALPHC must be >= 0.",
                     "In_Parameter", io::message::error);
  }
    
  switch(scale) {
    case 0 : // no scaling
      param.perturbation.scaling = false;
      param.perturbation.scaled_only = false;
      break;
    case 1 : // scaling on
      param.perturbation.scaling = true;
      param.perturbation.scaled_only = false;
      break;
    case 2 : // scaled only
      param.perturbation.scaling = true;
      param.perturbation.scaled_only = true;
      break;
    default :
      io::messages.add("PERTURBATION block: NSCALE must be 0 to 2.",
                       "In_Parameter", io::message::error);
  }
}

/**
 * read FORCE block.
 */
void io::In_Parameter::read_FORCE(simulation::Parameter &param,
				  std::ostream & os)
{
  DEBUG(8, "read FORCE");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["FORCE"];

  if (!buffer.size()){
    DEBUG(8, "no force block found???");
    io::messages.add("no FORCE block", "In_Parameter", io::message::error);
    return;
  }

  block_read.insert("FORCE");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  int bondH, angleH, impH, dihedralH;
  unsigned int num, e, old_e=0;
  
  _lineStream >> bondH >> param.force.bond >> angleH >> param.force.angle
	      >> impH >> param.force.improper >> dihedralH >> param.force.dihedral
	      >> param.force.nonbonded_crf >> param.force.nonbonded_vdw;
  _lineStream >> num;
  if(num<=0){
    DEBUG(10, "number of energy group < 0?");
    io::messages.add("FORCE block: number of energy groups should be > 0",
		     "In_Parameter", io::message::error);
    return;
  }
  
  for(unsigned int i=0; i<num; ++i){
    _lineStream >> e;
    DEBUG(10, "\tadding energy group " << e-1);
    param.force.energy_group.push_back(e-1);
    if(e<=old_e){
      DEBUG(10, "energy groups not in order...");
      io::messages.add("FORCE block: energy groups are not in order",
		       "In_Parameter", io::message::error);
      return;
    }
    old_e = e;
  }
  
  DEBUG(10, "number of energy groups: " << param.force.energy_group.size());

  if (_lineStream.fail())
    io::messages.add("FORCE block: bad line in ENERGYGROUP",
		     "In_Parameter", io::message::error);
  
  
  if (bondH ^ param.force.bond)
    io::messages.add("FORCE block: switch for bond and bond H has to be equal",
		     "In_Parameter", io::message::error);

  if (angleH ^ param.force.angle)
    io::messages.add("FORCE block: switch for angle and angle H has to be equal",
		     "In_Parameter", io::message::error);

  if (impH ^ param.force.improper)
    io::messages.add("FORCE block: switch for improper and improper H has to be equal",
		     "In_Parameter", io::message::error);

  if (dihedralH ^ param.force.dihedral)
    io::messages.add("FORCE block: switch for dihedral and dihedral H has to be equal",
		     "In_Parameter", io::message::error);

  if ((!param.force.nonbonded_crf) && param.force.nonbonded_vdw)
    io::messages.add("FORCE block: setting charges to zero",
		     "In_Parameter", io::message::notice);
  
  if (param.force.nonbonded_crf && (!param.force.nonbonded_vdw))
    io::messages.add("FORCE block: setting atom types to dummy",
		     "In_Parameter", io::message::notice);
  
  if (_lineStream.fail())
    io::messages.add("bad line in FORCE block", "In_Parameter", io::message::error);

  /*
  if (!_lineStream.eof())
    io::messages.add("End of line not reached in FORCE block, but should have been: \n" + s +  "\n",
		     "In_Parameter", io::message::warning);
  */
  
  if(param.force.bond < 0 || param.force.bond > 1)
    io::messages.add("FORCE block: Illegal value for force switch for bond",
		     "In_Parameter", io::message::error);
  if(param.force.angle < 0 || param.force.angle > 1)
    io::messages.add("FORCE block: Illegal value for force switch for angle",
		     "In_Parameter", io::message::error);
  if(param.force.improper < 0 || param.force.improper > 1)
    io::messages.add("FORCE block: Illegal value for force switch for improper dihedral",
		     "In_Parameter", io::message::error);
  if(param.force.dihedral < 0 || param.force.dihedral > 1)
    io::messages.add("FORCE block: Illegal value for force switch for dihedral",
		     "In_Parameter", io::message::error);
  if(param.force.nonbonded_vdw < 0 || param.force.nonbonded_vdw > 1)
    io::messages.add("FORCE block: Illegal value for force switch for nonbonded (vdw)",
		     "In_Parameter", io::message::error);
  if(param.force.nonbonded_crf < 0 || param.force.nonbonded_crf > 1)
    io::messages.add("FORCE block: Illegal value for force switch for nonbonded (crf)",
		     "In_Parameter", io::message::error);
}

/**
 * read COVALENTFORM block.
 */
void io::In_Parameter::read_COVALENTFORM(simulation::Parameter &param,
				       std::ostream & os)
{
  DEBUG(8, "read COVALENTFORM");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["COVALENTFORM"];

  if (!buffer.size()){
    return;
  }
  
  block_read.insert("COVALENTFORM");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  int bond, angle, dihedral;
  
  _lineStream >> bond 
	      >> angle
              >> dihedral;
  
  if (bond != 0 && bond != 1) {
    io::messages.add("COVALENTFORM block: NTBBH must be 0 (quartic) "
                     "or 1 (harmonic).",
		     "In_Parameter", io::message::error);
  } else {
    if (param.force.bond != 0) {
      switch (bond) {
        case 1: param.force.bond = 2; break;
        case 0:
        default: param.force.bond = 1;
      }
    }
  }
  
  if (angle != 0 && angle != 1) {
    io::messages.add("COVALENTFORM block: NTBAH must be 0 (quartic) "
                     "or 1 (harmonic).",
		     "In_Parameter", io::message::error);
  } else {
    if (param.force.angle != 0) {
      switch (angle) {
        case 1: param.force.angle = 2; break;
        case 0:
        default: param.force.angle = 1;
      }
    }
  }
  
  if (dihedral != 0 && dihedral != 1) {
    io::messages.add("COVALENTFORM block: NTBDN must be 0 (arbitray "
                     "phase shifts) or 1 (phase shifts limited).",
		     "In_Parameter", io::message::error);
  } else {
    if (dihedral != 0) {
      switch (angle) {
        case 1: 
         io::messages.add("COVALENTFORM block: NTBDN 1 not implemented.",
		          "In_Parameter", io::message::error);
          break;
        case 0:
        default: ;
      }
    }
  }
}

/**
 * read INITIALISE block.
 */
void io::In_Parameter::read_INITIALISE(simulation::Parameter &param,
				  std::ostream & os)
{
  DEBUG(8, "read INITIALISE");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["INITIALISE"];
  
  if (!buffer.size()){
    io::messages.add("no INITIALISE block", "In_Parameter", io::message::error);
    return;
  }

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  block_read.insert("INITIALISE");

  int ntivel, ntishk, ntinhc, ntishi, ntirtc, nticom, ntisti;
  _lineStream >> ntivel >> ntishk >> ntinhc 
              >> ntishi >> ntirtc >> nticom >> ntisti
	      >> param.start.ig 
	      >> param.start.tempi;
  
  if (_lineStream.fail())
    io::messages.add("bad line in INITIALISE block",
		     "In_Parameter", io::message::error);
  
  // generation of initial velocities
  switch(ntivel) {
    case 0 : param.start.generate_velocities = false; break;
    case 1 : param.start.generate_velocities = true; break;
    default : io::messages.add("INITIALISE block: NTIVEL must be 0 or 1",
		     "In_Parameter", io::message::error);
  }
  
  // controls initial SHAKE
  switch(ntishk) {
    case 0 : // no initial SHAKE
      param.start.shake_pos=false;
      param.start.shake_vel=false;
      break;
    case 1 : // SHAKE coordinates
      param.start.shake_pos=true;
      param.start.shake_vel=false;
      break;
    case 2 : // SHAKE velocities
      param.start.shake_pos=false;
      param.start.shake_vel=true;
      break;
    case 3 : // SHAKE coordinates & velocities
      param.start.shake_pos=true;
      param.start.shake_vel=true;
      break;
    default: io::messages.add("INITIALISE block: NTISHK must be 0 to 3",
		     "In_Parameter", io::message::error);
  }
  
  // controls reading of Nose-Hoover chain variables: not implemented.
  switch(ntinhc) {
    case 0 : param.start.read_nosehoover_chains = true; break;
    case 1 : param.start.read_nosehoover_chains = false; break; // reset them
    default : io::messages.add("INITIALISE block: NTINHC must be 0 or 1",
		     "In_Parameter", io::message::error);
  }
  
  io::messages.add("INITIALISE block: NTINHC is ignored.",
		   "In_Parameter", io::message::notice);
  
  // only for lattice sum. ignored without warning
  switch(ntishi) { 
    case 0: break;
    case 1: break;
    default : io::messages.add("INITIALISE block: NTISHI must be 0 or 1",
		     "In_Parameter", io::message::error);
  }
  
  // controls reading of restart data for roto-translational constraints:
  // not implemented.
  switch(ntirtc) {
    case 0: param.start.read_rottrans = true; break;
    case 1: param.start.read_rottrans = false; break;
    default : io::messages.add("INITIALISE block: NTIRTC must be 0 or 1",
		     "In_Parameter", io::message::error);
  }
  
  io::messages.add("INITIALISE block: NTIRTC is ignored.",
		     "In_Parameter", io::message::notice);
  
  // controls removal of COM translation and rotation.
  switch(nticom) {
    case 0:
      param.start.remove_com_rotation = false; 
      param.start.remove_com_translation = false;
      break;
    case 1: 
      param.start.remove_com_rotation = false; 
      param.start.remove_com_translation = true;
      break;
    case 2:
      param.start.remove_com_rotation = true; 
      param.start.remove_com_translation = true;
      break;
    default : io::messages.add("INITIALISE block: NTICOM must be 0 to 2",
		     "In_Parameter", io::message::error);
  }
  
  // controls reading of stochastic integrals: not implemented.
  switch(ntisti) {
    case 0: 
      param.stochastic.generate_integral = false; 
      break;
    case 1: 
      param.stochastic.generate_integral = true;
      break;
    default : io::messages.add("INITIALISE block: NTISTI must be 0 or 1",
		     "In_Parameter", io::message::error);
  }
  
  io::messages.add("INITIALISE block: NTISTI is ignored.",
		     "In_Parameter", io::message::notice);
  
  if(param.start.tempi <0)
    io::messages.add("Illegal value for TEMPI in INITIALISE block (>=0)",
		     "In_Parameter", io::message::error);
}
/**
 * read COMTRANSROT block.
 */
void io::In_Parameter::read_COMTRANSROT(simulation::Parameter &param,
					 std::ostream & os)
{
  DEBUG(8, "read COMTRANSROT");

  std::vector<std::string> buffer;
  std::string s;
  
  DEBUG(10, "reading COMTRANSROT block");
  buffer = m_block["COMTRANSROT"];
  
  if (!buffer.size()){
    return;
  }

  block_read.insert("COMTRANSROT");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  int nscm;
  
  _lineStream >> nscm;

  if (_lineStream.fail())
    io::messages.add("bad line in COMTRANSROT block",
		     "In_Parameter", io::message::error);
  
  if (nscm > 0){
    param.centreofmass.skip_step = nscm;
    param.centreofmass.remove_rot = false;
    param.centreofmass.remove_trans = true;
  } else if (nscm < 0) {
    param.centreofmass.skip_step = -nscm;
    param.centreofmass.remove_rot = true;
    param.centreofmass.remove_trans = true;    
  } else { // nscm == 0;
    param.centreofmass.skip_step = 0;
    param.centreofmass.remove_rot = false;
    param.centreofmass.remove_trans = false;
  }
}

/**
 * read LONGRANGE block.
 */
void io::In_Parameter::read_LONGRANGE(simulation::Parameter &param,
				      std::ostream & os)
{
  DEBUG(8, "read LONGRANGE");

  std::vector<std::string> buffer;
  std::string s;
  
  DEBUG(10, "longrange block");
  buffer = m_block["LONGRANGE"];
  
  if (!buffer.size()){
    io::messages.add("no LONGRANGE block in input",
		     "In_Parameter",io::message::error);
    return;
  }
  
  block_read.insert("LONGRANGE");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  _lineStream >> param.longrange.rf_epsilon 
	      >> param.longrange.rf_kappa 
	      >> param.longrange.rf_cutoff;
  param.longrange.epsilon = 1.0;
  
  if (_lineStream.fail())
    io::messages.add("bad line in LONGRANGE block",
		     "In_Parameter", io::message::error);
  

  if(param.longrange.rf_cutoff < 0) {
    param.longrange.rf_excluded=true;
    param.longrange.rf_cutoff = -param.longrange.rf_cutoff;
  }
  else{
    param.longrange.rf_excluded=false;
  }
  if(param.longrange.rf_epsilon!=0 && param.longrange.rf_epsilon<1)
    io::messages.add("LONGRANGE block: Illegal value for EPSRF (0  / >=1)", 
		     "In_Parameter", io::message::error);
  if(param.longrange.rf_kappa <0)
    io::messages.add("LONGRANGE block: Illegal value for APPAK (>=0)",
		     "In_Parameter", io::message::error);
} // LONGRANGE

/**
 * read PAIRLIST block.
 */
void io::In_Parameter::read_PAIRLIST(simulation::Parameter &param,
				  std::ostream & os)
{
  DEBUG(8, "read PAIRLIST");

  std::vector<std::string> buffer;
  std::string s;

  DEBUG(10, "pairlist block");
  
  // try a PAIRLIST
  buffer = m_block["PAIRLIST"];
  if (buffer.size()){
    block_read.insert("PAIRLIST");

    std::string s1, s2, s3;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> s1 
		>> param.pairlist.skip_step 
		>> param.pairlist.cutoff_short
		>> param.pairlist.cutoff_long 
		>> s2
		>> s3;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in PAIRLIST block",
		       "In_Parameter", io::message::error);
    }
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    std::transform(s3.begin(), s3.end(), s3.begin(), tolower);

    if (s1 == "grid") param.pairlist.grid = 1;    
    else if (s1 == "vgrid") param.pairlist.grid = 2;
    else if (s1 == "standard") param.pairlist.grid = 0;
    else{
      io::messages.add("PAIRLIST block: wrong pairlist algorithm chosen",
		       "In_Parameter", io::message::error);
      param.pairlist.grid = false;
    }
    
    if (param.pairlist.grid){
      if (s2 == "auto") 
	param.pairlist.grid_size = 0.5 * param.pairlist.cutoff_short;
      else{
	std::istringstream css;
	css.str(s2);
	css >> param.pairlist.grid_size;
	// param.pairlist.grid_size = atof(s2.c_str());
	if (css.fail()){
	  io::messages.add("PAIRLIST block: wrong pairlist grid size chosen",
			   "In_Parameter", io::message::error);
	  param.pairlist.grid_size = 0.5 * param.pairlist.cutoff_short;
	}
      }
    }
    else param.pairlist.grid_size = 0;

    if (s3 == "atomic") param.pairlist.atomic_cutoff = true;
    else if (s3 == "chargegroup") param.pairlist.atomic_cutoff = false;
    else {
      std::istringstream css;
      css.str(s3);
      css >> param.pairlist.atomic_cutoff;
      
      // param.pairlist.atomic_cutoff = (atoi(s3.c_str()) != 0);
      if (css.fail()){
	io::messages.add("PAIRLIST block: wrong cutoff type chosen "
                         "(allowed: atomic, chargegroup)",
			 "In_Parameter", io::message::error);
	param.pairlist.atomic_cutoff = false;
      }
    }
  }
  
  if(param.pairlist.grid && param.pairlist.grid_size <=0)
    io::messages.add("PAIRLIST block: Illegal value for grid size (>0)",
		     "In_Parameter", io::message::error);
  if(param.pairlist.cutoff_short < 0){
    io::messages.add("PAIRLIST block: Illegal value for short range cutoff (>0)",
		     "In_Parameter", io::message::error);
  }
  if(param.pairlist.cutoff_long < param.pairlist.cutoff_short){
    io::messages.add("PAIRLIST block: Illegal value for long range cutoff (>=RCUTP)",
		     "In_Parameter", io::message::error);
  }
}

/**
 * read the CGRAIN block.
 */
void io::In_Parameter::read_CGRAIN(simulation::Parameter &param,
				   std::ostream & os)
{
  DEBUG(8, "read CGRAIN");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["CGRAIN"];
  
  if (buffer.size()){
    
    block_read.insert("CGRAIN");
    
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    int ntcgran;
    _lineStream >> ntcgran >> param.cgrain.EPS;
    
    if (_lineStream.fail())
      io::messages.add("bad line in CGRAIN block",
                       "In_Parameter", io::message::error);

    switch(ntcgran) {
      case 0 :
        param.cgrain.level = 0;
        break;
      case 1 :
        param.cgrain.level = 1;
        param.force.interaction_function = simulation::cgrain_func;
        break;
      case 2 :
        param.cgrain.level = 2;
        param.force.interaction_function = simulation::cgrain_func;
        break;
      default :
        param.cgrain.level = 0;
        io::messages.add("CGRAIN block: NTCGRAN must be 0 to 2.",
                         "In_Parameter", io::message::error);       
    }
    
    if (param.cgrain.EPS < 0)
      io::messages.add("CGRAIN block: EPS must be >= 0.0.",
                       "In_Parameter", io::message::error);
    
  }
}

/**
 * read MULTIBATH block.
 */
void io::In_Parameter::read_MULTIBATH(simulation::Parameter &param,
				      std::ostream & os)
{
  DEBUG(8, "read MULTIBATH");

  std::vector<std::string> buffer;
  std::string s;
  
  DEBUG(10, "TEMPERATURE COUPLING block");
  
  param.multibath.couple = false;

  // is there a MULTIBATH block
  buffer = m_block["MULTIBATH"];
  
  if (buffer.size()){

    block_read.insert("MULTIBATH");

    param.multibath.found_multibath=true;
    param.multibath.found_tcouple=false;
    
    DEBUG(11, "MULTIBATH present");
    /*io::messages.add("using MULTIBATH block",
		     "In_Parameter", io::message::notice);*/

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    std::string alg;

    // the algorithm
    _lineStream >> alg;
    std::transform(alg.begin(), alg.end(), alg.begin(), tolower);

    if (alg == "weak-coupling")
      param.multibath.nosehoover = 0;
    else if (alg == "nose-hoover")
      param.multibath.nosehoover = 1;
    else if (alg == "nose-hoover-chains")
      param.multibath.nosehoover = 2;
    else{
      std::stringstream s(alg);
      if (!(s >> param.multibath.nosehoover) ||
	  param.multibath.nosehoover < 0 || param.multibath.nosehoover > 2){
	io::messages.add("MUTLIBATH block: algorithm not understood",
			 "In_Parameter", io::message::error);

	param.multibath.nosehoover = 0;
	return;
      }
    }

    if (param.multibath.nosehoover == 2){
      // read in the number of chains
      int num;
      _lineStream >> num;

      if (num < 2){
	io::messages.add("MULTIBATH block: wrong number of Nose-Hoover chains",
			 "In_Parameter", io::message::error);
	param.multibath.nosehoover = 0;
	return;
      }
      param.multibath.nosehoover = num;
    }
    
    int num;
    unsigned int last;
    unsigned int com_bath, ir_bath;
    double temp, tau;

    // the baths
    _lineStream >> num;
    
    for(int i=0; i<num; ++i){
      _lineStream >> temp >> tau;

      if (temp < 0.0 || (tau <= 0.0 && tau != -1)){
	io::messages.add("MULTIBATH block: illegal value for temp or tau",
			 "In_Parameter", io::message::error);
      }

      param.multibath.multibath.add_bath(temp, tau);
      if (tau != -1) param.multibath.couple = true;
    }
    
    if (_lineStream.fail()){
      io::messages.add("bad line in MULTIBATH block",
		       "In_Parameter", io::message::error);
    }
    
    // now the ranges
    _lineStream >> num;

    if (param.multibath.multibath.size() == 0 &&
	num > 0){
      
      io::messages.add("MULTIBATH block: no baths but coupling groups specified",
		       "In_Parameter", io::message::error);
      num = 0;
    }

    for(int i=0; i<num; ++i){
      _lineStream >> last >> com_bath >> ir_bath;
      // let it figure out the last molecule on its own
      
      if (last < 1 || com_bath < 1 || ir_bath < 1){
	io::messages.add("MULTIBATH block: range parameter < 1",
			 "In_Parameter", io::message::error);
	if (last < 1) last = 1;
	if (com_bath < 1) com_bath = 1;
	if (ir_bath < 1) ir_bath = 1;
      }

      if (com_bath > param.multibath.multibath.size() ||
	  ir_bath > param.multibath.multibath.size()){
	io::messages.add("MULTIBATH block: ir bath or com bath index too large",
			 "In_Parameter", io::message::error);
	if (com_bath > param.multibath.multibath.size()) com_bath = param.multibath.multibath.size();
	if (ir_bath > param.multibath.multibath.size()) ir_bath = param.multibath.multibath.size();
      }

      param.multibath.multibath.add_bath_index(last - 1, 0, com_bath - 1, ir_bath - 1);
    }
    
    if (_lineStream.fail()){
      io::messages.add("bad line in MULTIBATH block",
		       "In_Parameter", io::message::error);
    }
    
  }
} // MULTIBATH

/**
 * read POSITIONRES block.
 */
void io::In_Parameter::read_POSITIONRES(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read POSITIONRES");

  std::vector<std::string> buffer;
  std::string s;
  
  DEBUG(10, "positionres block");
  buffer = m_block["POSITIONRES"];
  
  if (!buffer.size()){
    return;
  }
  
  block_read.insert("POSITIONRES");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  int ntpor, ntpors, read;
  _lineStream >> ntpor
              >> read
              >> ntpors
	      >> param.posrest.force_constant;
  
  if (_lineStream.fail())
    io::messages.add("bad line in POSITIONRES block",
		     "In_Parameter", io::message::error);
  
    switch(ntpor) {
    case 0 :
      param.posrest.posrest = simulation::posrest_off;
      break;
    case 1 :
      param.posrest.posrest = simulation::posrest_on;
      break;
    case 2 :
      param.posrest.posrest = simulation::posrest_bfactor;
      break;
    case 3 :
      param.posrest.posrest = simulation::posrest_const;
      break;
    default:
      io::messages.add("POSITIONRES block: NTPOR must be 0 to 3.",
		       "In_Parameter", io::message::error);             
  }

  switch(read) {
    case 0 :
      param.posrest.read = false;
      break;
    case 1 :
      param.posrest.read = true;
      break;
    default:
      param.posrest.read = false;
      io::messages.add("POSITIONRES block: NTPORS must be 0 or 1.",
		       "In_Parameter", io::message::error);          
  }

  switch(ntpors) {
    case 0 :
      param.posrest.scale_reference_positions = false;
      break;
    case 1 :
      param.posrest.scale_reference_positions = true;
      break;
    default:
      param.posrest.scale_reference_positions = false;
      io::messages.add("POSITIONRES block: NTPORS must be 0 or 1.",
		       "In_Parameter", io::message::error);          
  }
  
  if(param.posrest.force_constant < 0.0)
    io::messages.add("POSITIONRES block: Illegal value for CPOR.",
		     "In_Parameter", io::message::error);

} // POSITIONRES

/**
 * read DISTANCERES block.
 */
void io::In_Parameter::read_DISTANCERES(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read DISTANCERES");

  std::vector<std::string> buffer;
  std::string s;
  
  DEBUG(10, "distenceres block");
  buffer = m_block["DISTANCERES"];
  
  if (!buffer.size()){
    return;
  }
  
  block_read.insert("DISTANCERES");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  int ntdira;
  _lineStream >> param.distanceres.distanceres
  	      >> ntdira
	      >> param.distanceres.K
	      >> param.distanceres.r_linear
	      >> param.distanceres.tau;
  
  if (_lineStream.fail())
    io::messages.add("bad line in DISTANCERES block",
		     "In_Parameter", io::message::error);
  
  
  if(param.distanceres.distanceres < -2 || param.distanceres.distanceres > 2) {
    io::messages.add("DISTANCERES block: NTDIR must 0 to 2.",
                     "In_Parameter", io::message::error);
  }
  
  switch(ntdira) {
    case 0 : param.distanceres.read = false; break;
    case 1 : param.distanceres.read = true; break;
    default:  param.distanceres.read = false; 
      io::messages.add("DISTANCERES block: NTDIRA must be 0 or 1.",
                       "In_Parameter", io::message::error);
  }
  
  if(param.distanceres.tau < 0.0) {
    io::messages.add("DISTANCERES block: TAUDIR must be >= 0.0.",
                     "In_Parameter", io::message::error);
  }

  if(param.distanceres.K <0) {
    io::messages.add("DISTANCERES block: CDIR must be >= 0.0.",
		     "In_Parameter", io::message::error);
  }

} // DISTANCERES

/**
 * read DIHEDRALRES block.
 */
void io::In_Parameter::read_DIHEDRALRES(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read DIHEDRALRES");

  std::vector<std::string> buffer;
  std::string s;
  double phi_lin;
  
  DEBUG(10, "DIHEDRALRES block");
  buffer = m_block["DIHEDRALRES"];
  
  if (!buffer.size()){
    return;
  }
  
  block_read.insert("DIHEDRALRES");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  /*
   * #         0:    off [default]
   * #         1:    instantaneous dihedral restraining using CDLR
   * #         2:    time-averaged dihedral restraining using CDLR
   * #         3:    biquadratic dihedral restraining using CDLR
   * #         4:    dihedral constraining
   * # CDLR    >=0.0 force constant for dihedral restraining (multiplied by WDLR)
   */
  
   int ntdlr;
  _lineStream >> ntdlr // param.dihrest.dihrest
	      >> param.dihrest.K
	      >> phi_lin;

  param.dihrest.phi_lin = phi_lin * 2 * math::Pi / 360;
  
  if (_lineStream.fail())
    io::messages.add("bad line in DIHEDRALRES block",
		     "In_Parameter", io::message::error);
  
  switch(ntdlr){
    case 0: param.dihrest.dihrest = ntdlr; break;
    case 1: param.dihrest.dihrest = 2; break;
    case 2: {
      io::messages.add("DIHEDRALRES block: NTDLR = 2 not implemented.",
                       "In_Parameter", io::message::error);
      break;
    }
    case 3: {
      io::messages.add("DIHEDRALRES block: NTDLR = 3 not implemented.",
                       "In_Parameter", io::message::error);
      break;
    }
    case 4: param.dihrest.dihrest = 3; break;
    default: {
       io::messages.add("DIHEDRALRES block: NTDLR must be 0...4.",
                        "In_Parameter", io::message::error);     
    }
  }
  
  if(param.dihrest.K < 0)
    io::messages.add("DIHEDRALRES block: Illegal value for force constant (>=0)",
		     "In_Parameter", io::message::error);

  if (param.dihrest.dihrest == 3){
    if (param.constraint.ntc == 1 && param.constraint.solute.algorithm == simulation::constr_off)
      param.constraint.solute.algorithm = simulation::constr_shake;

    if (param.constraint.solute.algorithm != simulation::constr_shake){
      io::messages.add("DIHEDRALRES block: needs SHAKE as (solute) constraints algorithm",
		       "In_Parameter",
		       io::message::error);
    }
  }
  
} // DIHREST

/**
 * read the JVALUERES block.
 */
void io::In_Parameter::read_JVALUERES(simulation::Parameter &param,
				   std::ostream & os)
{
  DEBUG(8, "read JVALUERES");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["JVALUERES"];
  if (buffer.size()){

    block_read.insert("JVALUERES");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    std::string s1;
    _lineStream >> s1 // NTJVR
                >> param.jvalue.read_av // NTJVRA
                >> param.jvalue.K // CJVR
                >> param.jvalue.tau // TAUJVR
                >> param.jvalue.le // LE
                >> param.jvalue.ngrid // NGRID
                >> param.jvalue.delta; // DELTA
      
    if (_lineStream.fail())
      io::messages.add("bad line in JVALUERES block",
		       "In_Parameter", io::message::error);
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    
    // NTJVR 0...3
    if (s1 == "off") param.jvalue.mode = simulation::restr_off;
    else if (s1 == "instantaneous") param.jvalue.mode = simulation::restr_inst;
    else if (s1 == "averaged") param.jvalue.mode = simulation::restr_av;
    else if (s1 == "biquadratic") param.jvalue.mode = simulation::restr_biq;
    else{
      std::istringstream css;
      css.str(s1);
      
      int i;
      css >> i;
      if(css.fail() || i < 0 || i > 3){
	io::messages.add("JVALURES block: bad value for MODE:"+s1+"\n"
			 "off, instantaneous, averaged, biquadratic (0-3)",
			 "In_Parameter", io::message::error);
	param.jvalue.mode = simulation::restr_off;
	return;
      }
      param.jvalue.mode = simulation::restr_enum(i);

      DEBUG(10, "setting jvalue mode to " << param.jvalue.mode);

    }
    if (param.jvalue.tau < 0 ||
	(param.jvalue.tau == 0 && (param.jvalue.mode != simulation::restr_off ||
				   param.jvalue.mode != simulation::restr_inst))){
      io::messages.add("JVALUERES block: bad value for TAU, should be > 0.0",
		       "In_Parameter", io::message::error);
    }
    if (param.jvalue.mode != simulation::restr_off && param.jvalue.K < 0.0){
      io::messages.add("JVALUERES block: bad value for K in JVALUERES block,"
		       "should be > 0.0",
		       "In_Parameter", io::message::error);
    }
    if (param.jvalue.le > 0){
      
      if (param.jvalue.ngrid < 1){
	io::messages.add("JVALUERES block: bad value for NGRID in JVALUERES block, "
			 "should be > 1",
			 "In_Parameter", io::message::error);
      }
    }
    if (param.jvalue.read_av && (param.jvalue.mode != simulation::restr_av
        && !param.jvalue.le)){
      io::messages.add("JVALUERES block: Continuation only needed "
                       "with averaging or LE.",
		       "In_Parameter", io::message::error);
    }
  } // JVALUERES
} // JVALUE

/**
 * read the PERSCALE block.
 */
void io::In_Parameter::read_PERSCALE(simulation::Parameter &param,
				   std::ostream & os)
{
  DEBUG(8, "read PERSCALE");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["PERSCALE"];
  if (buffer.size()){

    block_read.insert("PERSCALE");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    std::string s1;
    int read;
    _lineStream >> s1
                >> param.pscale.KDIH >> param.pscale.KJ 
                >> param.pscale.T >> param.pscale.difference
		>> param.pscale.ratio >> read;
    if (_lineStream.fail()){
      io::messages.add("bad line in PERSCALE block",
		       "In_Parameter", io::message::error);
      return;
    }
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    
    if (s1 == "jrest" || s1 == "1"){
      param.pscale.jrest = true;
    } else if(s1 == "off" || s1 == "0") {
      param.pscale.jrest = false;
    } else {
      io::messages.add("PERSCALE block: RESTYPE must be jrest or off.",
		       "In_Parameter", io::message::error);     
    }

    if (param.pscale.KDIH < 0.0)
      io::messages.add("PERSCALE block: KDIH must be >= 0.0.",
                       "In_Parameter", io::message::error);
    if (param.pscale.KJ < 0.0)
      io::messages.add("PERSCALE block: KJ must be >= 0.0.",
                       "In_Parameter", io::message::error);
    if (param.pscale.T < 0.0)
      io::messages.add("PERSCALE block: T must be >= 0.0.",
                       "In_Parameter", io::message::error);
    if (param.pscale.difference < 0.0)
      io::messages.add("PERSCALE block: DIFF must be >= 0.0.",
                       "In_Parameter", io::message::error);
    
    switch(read) {
      case 0 :
         param.pscale.read_data = false;
         break;
      case 1 :
         param.pscale.read_data = true;
         break;
      default :
      io::messages.add("PERSCALE block: READ must be 0 or 1.",
                       "In_Parameter", io::message::error);        
    }
  }
} // PERSCALE


/**
 * read the ROTTRANS block.
 */
void io::In_Parameter::read_ROTTRANS(simulation::Parameter &param,
				     std::ostream & os)
{
  DEBUG(8, "read ROTTRANS");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["ROTTRANS"];

  if (buffer.size()){
    block_read.insert("ROTTRANS");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    int rtc;
    _lineStream >> rtc >> param.rottrans.last;
    
    if (_lineStream.fail())
      io::messages.add("bad line in ROTTRANS block",
		       "In_Parameter", io::message::error);

    switch(rtc) {
      case 0 :
        param.rottrans.rottrans = false;
        break;
      case 1 :
        param.rottrans.rottrans = true;
        break;
      default :
      io::messages.add("ROTTRANS block: RTC must be 0 or 1",
                       "In_Parameter", io::message::error);
    }

    if (param.rottrans.last <= 0)
      io::messages.add("ROTTRANS block: RTCLAST must be > 0.",
		       "In_Parameter", io::message::error);
  }
}

/**
 * read the INNERLOOP block.
 */
void io::In_Parameter::read_INNERLOOP(simulation::Parameter &param,
				      std::ostream & os)
{
  DEBUG(8, "read INNERLOOP");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["INNERLOOP"];

  if (buffer.size()){
 
    block_read.insert("INNERLOOP");
    
    int spc;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
 
    _lineStream >> spc;
 
    if (_lineStream.fail())
      io::messages.add("bad line in INNERLOOP block",
                       "In_Parameter", io::message::error);
 
    switch(spc){
      case 0: {
        // standard solvent loops
        param.force.spc_loop = -1;
        break;
      }
      case 1: {
        // fast spc loops
        param.force.spc_loop = 0;
        // spc_loop will be set to 1 in check_spc_loop (if all tests are ok).
        break;
      }
      default: {
        param.force.spc_loop = -1;
        io::messages.add("INNERLOOP block: bad value for SPC, allowed : 1 (on), 0 (off)",
                "In_Parameter",
                io::message::error);
      }
    }
  }
}

/**
 * read the REPLICA block.
 */
void io::In_Parameter::read_REPLICA(simulation::Parameter &param,
				      std::ostream & os)
{
  DEBUG(8, "read REPLICA");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["REPLICA"];
  if (buffer.size()){
    block_read.insert("REPLICA");

    bool error = false;
    
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.replica.num_T;
    
    if (_lineStream.fail() || param.replica.num_T < 0) {
      io::messages.add("REPLICA block: NRET must be >= 0.",
		       "In_Parameter", io::message::error);
      error = true;
    }
    
    if (!error)
      param.replica.temperature.resize(param.replica.num_T, 0.0);

    for(int i=0; i<param.replica.num_T; ++i){
      _lineStream >> param.replica.temperature[i];
      if (_lineStream.fail() || param.replica.temperature[i] < 0.0) {
        std::ostringstream msg;
        msg << "REPLICA block: RET(" << i+1 << ") must be >= 0.0";
        io::messages.add(msg.str(), "In_Parameter", io::message::error);
        error = true;
      }
    }
    
    int scale;
    _lineStream >> scale;
    switch(scale) {
      case 0 : 
        param.replica.scale = false;
        break;
      case 1 :
        param.replica.scale = true;
        break;
      default :
       io::messages.add("REPLICA block: LRSCALE must be 0 or 1",
		       "In_Parameter", io::message::error); 
       error = true;
    }

    if (_lineStream.fail() || error){
      param.replica.num_T = 0;
      param.replica.num_l = 0;

      param.replica.temperature.clear();
      param.replica.lambda.clear();
      param.replica.dt.clear();
      error = false;
    }
    
    _lineStream >> param.replica.num_l;
    
    if (_lineStream.fail() || param.replica.num_l < 0) {
      io::messages.add("REPLICA block: NRELAM must be >= 0.",
		       "In_Parameter", io::message::error);
      error = true;
    }
    
    if (!error) {
      param.replica.lambda.resize(param.replica.num_l, 0.0);
      param.replica.dt.resize(param.replica.num_l, 0.0);
    }
    
    for(int i=0; i<param.replica.num_l; ++i){
      _lineStream >> param.replica.lambda[i];
      if (_lineStream.fail() || param.replica.lambda[i] < 0.0) {
        std::ostringstream msg;
        msg << "REPLICA block: RELAM(" << i+1 << ") must be >= 0.0";
        io::messages.add(msg.str(), "In_Parameter", io::message::error);
        error = true;
      }
    }
    for(int i=0; i<param.replica.num_l; ++i){
      _lineStream >> param.replica.dt[i];
      if (_lineStream.fail() || param.replica.dt[i] < 0.0) {
        std::ostringstream msg;
        msg << "REPLICA block: RETS(" << i+1 << ") must be >= 0.0";
        io::messages.add(msg.str(), "In_Parameter", io::message::error);
        error = true;
      }
    }

    if (_lineStream.fail() || error){
      param.replica.num_T = 0;
      param.replica.num_l = 0;

      param.replica.temperature.clear();
      param.replica.lambda.clear();
      param.replica.dt.clear();
      error = false;
    }
    
    _lineStream >> param.replica.trials;
    if (_lineStream.fail() || param.replica.trials < 0) {
      io::messages.add("REPLICA block: NRETRIAL must be >= 0.",
		       "In_Parameter", io::message::error);
      error = true;
    }
    _lineStream >> param.replica.equilibrate;
    if (_lineStream.fail() || param.replica.equilibrate < 0) {
      io::messages.add("REPLICA block: NREQUIL must be >= 0.",
		       "In_Parameter", io::message::error);
      error = true;
    }
    _lineStream >> param.replica.slave_runs;
    if (_lineStream.fail() || param.replica.slave_runs < 0) {
      io::messages.add("REPLICA block: NREJOB must be >= 0.",
		       "In_Parameter", io::message::error);
      error = true;
    }
    _lineStream >> param.replica.write;
    if (_lineStream.fail() || param.replica.write < 0) {
      io::messages.add("REPLICA block: NREWRT must be >= 0.",
		       "In_Parameter", io::message::error);
      error = true;
    }
    
    if (_lineStream.fail() || error){
      io::messages.add("bad line in REPLICA block (trials, equi, slave or write)",
		       "In_Parameter", io::message::error);

      param.replica.num_T = 0;
      param.replica.num_l = 0;

      param.replica.temperature.clear();
      param.replica.lambda.clear();
      param.replica.dt.clear();
    }
  }
}

void io::In_Parameter::read_MULTICELL(simulation::Parameter & param,
				      std::ostream & os)
{
  DEBUG(8, "read MULTICELL");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["MULTICELL"];

  if (buffer.size()){

    block_read.insert("MULTICELL");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    int ntm;
    double tolpx, tolpv, tolpf, tolpfw;
    _lineStream >> ntm 
		>> param.multicell.x >> param.multicell.y >> param.multicell.z
                >> tolpx >> tolpv >> tolpf >> tolpfw;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in MULTICELL block",
		       "In_Parameter", io::message::error);

      param.multicell.multicell = false;
    }
    
    switch(ntm) {
      case 1 :
        param.multicell.multicell = true;
        break;
      case 0 :
        param.multicell.multicell = false;
        break;
      default :
        param.multicell.multicell = false;
        io::messages.add("MULTICELL block: NTM must be 0 or 1.",
                         "In_Parameter", io::message::error);
    }
    
    if (param.multicell.multicell) {
      // disable because broken
      param.multicell.multicell = false;
      io::messages.add("MULTICELL simulations are broken in MD++",
                         "In_Parameter", io::message::error);      
      
      
      // do these checks only if mutlicell is really used.
      if (param.multicell.x < 1 || param.multicell.y < 1 ||
          param.multicell.z < 1) {
        io::messages.add("MULTICELL block: NCELLA, NCELLB and NCELLC "
                         "must be >= 1.", "In_Parameter", io::message::error);
      }
    
      if (param.multicell.x == 1 && param.multicell.y == 1 && 
          param.multicell.z == 1) {
        io::messages.add("MULTICELL block: NCELLA, NCELLB and NCELLC are all 1.\n"
                         "disabling MULTICELL simulation.", "In_Parameter",
                         io::message::warning);    
        param.multicell.multicell = false;
      }
    
      if (tolpx || tolpv || tolpf || tolpfw) {
        io::messages.add("MULTICELL block: Periodicity checks not available in "
                         "this version. Disabling MULTICELL simulation.", 
                         "In_Parameter", io::message::warning);    
        param.multicell.multicell = false;
      }  
    }
  }
}

void io::In_Parameter::read_READTRAJ(simulation::Parameter & param,
				    std::ostream & os)
{
  DEBUG(8, "read READTRAJ");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["READTRAJ"];

  if (buffer.size()){

    block_read.insert("READTRAJ");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    int ntrd, ntrn, ntrb, ntshk;
    _lineStream >> ntrd >> ntrn >> ntrb >> ntshk;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in READTRAJ block",
		       "In_Parameter", io::message::error);
      
      param.analyze.analyze = false;
      return;
    }
    
    switch(ntrd) {
      case 1 :
        param.analyze.analyze = true;
        break;
      case 0 :
        param.analyze.analyze = false;
        break;
      default :
        param.analyze.analyze = false;
        io::messages.add("READTRAJ block: NTRD must be 0 or 1",
                         "In_Parameter", io::message::error);
    }
    
    if (ntrn)
      io::messages.add("READTRAJ block: NTRN was ignored", "In_Parameter",
                       io::message::warning);
    
    if (ntrb != 1)
      io::messages.add("READTRAJ block: NTRB must be 1.", "In_Parameter",
                       io::message::error);
    
    switch(ntshk) {
      case 1 :
        param.analyze.copy_pos = true;
        break;
      case 0 :
        param.analyze.copy_pos = false;
        break;
      default :
        param.analyze.copy_pos = false;
        io::messages.add("READTRAJ block: NTSHK must be 0 or 1",
                         "In_Parameter", io::message::error);
    }
    
  }
}

void io::In_Parameter::read_INTEGRATE(simulation::Parameter & param,
				      std::ostream & os)
{
  DEBUG(8, "read INTEGRATE");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["INTEGRATE"];
  
  if (buffer.size()){

    block_read.insert("INTEGRATE");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    int nint;
    _lineStream >> nint;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in INTEGRATE block",
		       "In_Parameter", io::message::error);
      
    }
    
    switch(nint) {
      case 0 :
        param.integrate.method = simulation::integrate_off;
        break;
      case 1 :
        param.integrate.method = simulation::integrate_leap_frog;
        break;
      default :
        io::messages.add("INTEGRATE block: NINT must be 0 or 1",
                         "In_Parameter", io::message::error);        
    }
  }
}

void io::In_Parameter::read_STOCHDYN(simulation::Parameter & param,
				       std::ostream & os)
{
  DEBUG(8, "read STOCHDYN");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["STOCHDYN"];
  
  if (buffer.size()){

    block_read.insert("STOCHDYN");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.stochastic.sd >> param.stochastic.ntfr
                >> param.stochastic.nsfr >> param.stochastic.nbref
		>> param.stochastic.rcutf >> param.stochastic.cfric
                >> param.stochastic.temp;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in STOCHDYN block",
		       "In_Parameter", io::message::error);
      return;
    }
    
    if(param.stochastic.sd < 0 || param.stochastic.sd > 1)
      io::messages.add("STOCHDYN block: NTSD must be 0 or 1",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.ntfr < 0 || param.stochastic.ntfr > 3)
      io::messages.add("STOCHDYN block: NTFR must be 0 to 3",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.nsfr <= 0)
      io::messages.add("STOCHDYN block: NSFR must be > 0",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.nbref <= 0)
      io::messages.add("STOCHDYN block: NBREF must be > 0",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.rcutf < 0)
      io::messages.add("STOCHDYN block: RCUTF must be >= 0",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.cfric < 0)
      io::messages.add("STOCHDYN block: CFRIC must be >= 0",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.temp < 0)
      io::messages.add("STOCHDYN block: TEMPSD must be >= 0",
                       "In_Parameter", io::message::error);
  }
}

void io::In_Parameter::read_EWARN(simulation::Parameter & param,
				  std::ostream & os)
{
  DEBUG(8, "read EWARN");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["EWARN"];
  
  if (buffer.size()){

    block_read.insert("EWARN");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.ewarn.limit;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in EWARN block",
		       "In_Parameter", io::message::error);
      return;
    }
  } else {
#ifndef HAVE_ISNAN
      io::messages.add("std::isnan() is not available for your compilation. "
                       "Consider using the EWARN block.",
		       "In_Parameter", io::message::warning);    
#endif
  }
}

void io::In_Parameter::read_MULTISTEP(simulation::Parameter & param,
				      std::ostream & os)
{
  DEBUG(8, "read MULTISTEP");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["MULTISTEP"];
  
  if (buffer.size()){
    block_read.insert("MULTISTEP");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    int boost; 
    _lineStream >> param.multistep.steps >> boost;

    if (_lineStream.fail()){
      io::messages.add("bad line in MULTISTEP block",
		       "In_Parameter", io::message::error);
      return;
    }

    if (param.multistep.steps < 0) {
      io::messages.add("MULTISTEP block: STEPS must be >= 0.",
                       "In_Parameter", io::message::error);
    }
    switch(boost) {
      case 0 :
        param.multistep.boost = false;
        break;
      case 1 :
        param.multistep.boost = true;
        break;
      default :
        io::messages.add("MULTISTEP block: BOOST must be 0 or 1",
                         "In_Parameter", io::message::error);
    }
  }
}

void io::In_Parameter::read_MONTECARLO(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read MONTECARLO");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["MONTECARLO"];
  
  if (buffer.size()){
    block_read.insert("MONTECARLO");
    
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    int mc;
    _lineStream >> mc >> param.montecarlo.steps >> param.montecarlo.dlambda;

    if (_lineStream.fail()){
      io::messages.add("bad line in MONTECARLO block",
              "In_Parameter", io::message::error);
      return;
    }
    
    switch(mc){
      case 0: param.montecarlo.mc = mc; break;
      case 1: param.montecarlo.mc = mc; break;
      default : {
        io::messages.add("MONTECARLO block: MC must be 0 or 1",
                "In_Parameter", io::message::error);
        break;
      }
    }
          
    // parameters are all positive
    if(param.montecarlo.steps < 0 || param.montecarlo.dlambda < 0){
      io::messages.add("MONTECARLO block: Negative parameter",
              "In_Parameter", io::message::error);
      return;
    }
    // perturbation
    if(param.montecarlo.mc && !param.perturbation.perturbation){
      io::messages.add("Chemical MONTECARLO only possible if perturbation is on."
      " Set NTG to 1.",
              "In_Parameter", io::message::error);
      return;
    }
  }
}

void io::In_Parameter::read_RAMD(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read RAMD");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["RAMD"];
  
  if (buffer.size()){
    block_read.insert("RAMD");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.ramd.fc 
		>> param.ramd.steps
		>> param.ramd.r_min
		>> param.ramd.every
		>> param.ramd.tau
		>> param.ramd.ta_min;
    

    if (_lineStream.fail()){
      io::messages.add("bad line in RAMD block",
		       "In_Parameter", io::message::error);
      return;
    }
    int numAtoms=0, atom;
    _lineStream >> numAtoms;
    for(int i=0; i< numAtoms; i++){
      _lineStream >> atom;
      param.ramd.atom.insert(atom-1);
    }
    if (_lineStream.fail()){
      io::messages.add("bad line in RAMD block",
		       "In_Parameter", io::message::error);
    }

    if(param.ramd.fc < 0)
      io::messages.add("RAMD block: FC should be >=0",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);
    if(param.ramd.steps < 0)
      io::messages.add("RAMD block: STEPS should be >0",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);
    if(param.ramd.r_min < 0)
      io::messages.add("RAMD block: R_MIN should be >0",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);
    if(param.ramd.every < 0)
      io::messages.add("RAMD block: NWRITE should be >=0",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);
    if(param.ramd.tau ==0)
      param.ramd.do_ta=false;
    else if(param.ramd.tau < 0)
      io::messages.add("RAMD block: TAU should be >=0",
		      "io::In_Parameter::read_RAMD",
		      io::message::error);
    else
      param.ramd.do_ta=true;
    if(param.ramd.do_ta && param.ramd.ta_min <= 0)
      io::messages.add("RAMD block: TA_MIN should be >0 if time averaging is included",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);

    if(param.ramd.fc!=0.0 && param.ramd.atom.size()==0)
      io::messages.add("RAMD block: no atoms read in to apply random force",
		       "io::In_Parameter::read_RAMD",
		       io::message::warning);


  }
}

void io::In_Parameter::read_POLARIZE(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read POLARIZE");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["POLARIZE"];
  
  if (buffer.size()) {
    block_read.insert("POLARIZE");
    
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    int cos, damp, efield;
    _lineStream >> cos >> efield >> param.polarize.minfield >> damp
                >> param.polarize.write;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in POLARIZE block",
		       "In_Parameter", io::message::error);
      return;
    }
    
    switch(cos) {
      case 0 : param.polarize.cos = false; break;
      case 1 : {
        param.polarize.cos = true; 
        param.force.interaction_function = simulation::pol_lj_crf_func;
        break;
      }
      default:
        io::messages.add("POLARIZE block: COS must be 0 or 1.",
                         "In_Parameter", io::message::error);
    }
    
    switch(efield) {
      case 0 : param.polarize.efield_site = simulation::ef_atom; break;
      case 1 : param.polarize.efield_site = simulation::ef_cos; break;
      default:
        io::messages.add("POLARIZE block: EFIELD must be 0 or 1.",
                         "In_Parameter", io::message::error);
    }
    
    if (param.polarize.minfield <= 0.0) {
      io::messages.add("POLARIZE block: MINFIELD must be > 0.0",
                         "In_Parameter", io::message::error); 
    }
    
    switch(damp) {
      case 0 : param.polarize.damp = false; break;
      case 1 : param.polarize.damp = true; break;
      default:
        io::messages.add("POLARIZE block: DAMP must be 0 or 1.",
                         "In_Parameter", io::message::error);
    }
    
    if (!param.polarize.cos)
      param.polarize.write = 0;
    
    if (param.polarize.write < 0) {
      io::messages.add("POLARIZE block: WRITE must be >= 0",
                         "In_Parameter", io::message::error);
    }
    
    if (param.polarize.damp && !param.polarize.cos) {
      io::messages.add("POLARIZE block: DAMP is ignored if no polarization is used",
                       "In_Parameter", io::message::warning);
    }   
  }
}

void io::In_Parameter::read_RANDOMNUMBERS(simulation::Parameter & param, 
        std::ostream & os)
{
  DEBUG(8, "read RANDOMNUMBERS");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["RANDOMNUMBERS"];
  
  if (buffer.size()) {
    block_read.insert("RANDOMNUMBERS");
    
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    int rng;
    _lineStream >> rng >> param.rng.gsl_rng;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in RANDOMNUMBERS block",
		       "In_Parameter", io::message::error);
      return;
    }
    
    switch(rng) {
      case 0 : 
        param.rng.rng = simulation::random_g96;
        break;
      case 1 :
        param.rng.rng = simulation::random_gsl;
        break;
      default :
       io::messages.add("RANDOMNUMBERS block: NTRNG must be 0 (G96) "
                        "or 1 (GSL)", "In_Parameter", io::message::error);       
    }  
    
    math::RandomGenerator::check(param);
  }
}  

void io::In_Parameter::read_EDS(simulation::Parameter & param, 
        std::ostream & os)
{
  DEBUG(8, "read EDS");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["EDS"];
  
  if (buffer.size()) {
    block_read.insert("EDS");
    
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    int eds;
    _lineStream >> eds >> param.eds.s >> param.eds.numstates;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in EDS block",
		       "In_Parameter", io::message::error);
      return;
    }
    param.eds.eir.resize(param.eds.numstates,0.0);
    for (unsigned int i = 0; i < param.eds.numstates; i++){
      _lineStream >> param.eds.eir[i];
    }
    if (_lineStream.fail()){
      io::messages.add("Error when reading EIR from EDS block",
                       "In_Parameter", io::message::error);
      return;
    }
    
    switch(eds) {
      case 0 : 
        param.eds.eds = false;
        param.eds.numstates = 0;
        break;
      case 1 :
        param.eds.eds = true;
        break;
      default :
       io::messages.add("Error in EDS block: EDS must be 0 (no EDS) "
       "or 1 (EDS)", "In_Parameter", io::message::error);
    }
    if(param.eds.eds){
      if(param.eds.s <= 0){
        io::messages.add("Error in EDS block: S must be >0",
                "In_Parameter", io::message::error);
      }
      
      if(param.eds.numstates < 1){
        io::messages.add("Error in EDS block: NUMSTATES must be >=1.",
                "In_Parameter", io::message::error);
      }
      // make sure we simulate at a given temperature (unambiguous kT)
      if(!param.multibath.couple){
        io::messages.add("Error in EDS block: EDS requires temperature coupling.",
                "In_Parameter", io::message::error);
      }
      // check whether all baths have the same temperature (unambiguous kT)
      for(unsigned int i = 1; i < param.multibath.multibath.size(); i++){
        if(param.multibath.multibath.bath(i).temperature !=
                param.multibath.multibath.bath(0).temperature){
          io::messages.add("Error in EDS block: all baths must have the same temperature.",
                  "In_Parameter", io::message::error);
        }
      }
    }

    
    
  }
}

void io::In_Parameter::read_LAMBDAS(simulation::Parameter & param, 
        std::ostream & os)
{
  DEBUG(8, "read LAMBDAS");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["LAMBDAS"];
  
  if (buffer.size()) {
    block_read.insert("LAMBDAS");
    
    io::messages.add("LAMBDAS block: not implemented",
		     "In_Parameter", io::message::error);
  }
}  

// two helper data types to simply unsupported block handling
enum unsupported_block_type { 
  ub_unknown, // I know it and know that I don't use it but I have no idea why
  ub_renamed, // it was renamed. e.g. from previous versions
  ub_promd, // it is a PROMD block. Tell alternative if there is any
  ub_g96 // it is a G96 block. Tell alternative if there is any
};

// give a valid block name as an alternative and it will tell the user to
// use it.
struct unsupported_block {
  unsupported_block() : 
        alternative (""), type(ub_unknown) {}
  unsupported_block(std::string a, unsupported_block_type t) :
        alternative(a), type(t) {}
        
        std::string alternative;
        unsupported_block_type type;
};

void io::In_Parameter::read_known_unsupported_blocks() {
  std::map<std::string,unsupported_block> ub;
  // add all those unknown blocks
  ub["ANATRAJ"] = unsupported_block("READTRAJ", ub_renamed);
  ub["MINIMISE"] = unsupported_block("ENERGYMIN", ub_renamed);
  ub["STOCHASTIC"] = unsupported_block("STOCHDYN", ub_renamed);
  ub["BOUNDARY"] = unsupported_block("BOUNDCOND", ub_renamed);
  ub["THERMOSTAT"] = unsupported_block("MULTIBATH", ub_promd);
  ub["TCOUPLE"] = unsupported_block("MULTIBATH", ub_g96);
  ub["BAROSTAT"] = unsupported_block("PRESSURESCALE", ub_promd);
  ub["VIRIAL"] = unsupported_block("PRESSURESCALE", ub_promd);
  ub["PCOUPLE"] = unsupported_block("PRESSURESCALE", ub_g96);
  ub["PCOUPLE03"] = unsupported_block("PRESSURESCALE", ub_renamed);
  ub["GEOMCONSTRAINT"] = unsupported_block("CONSTRAINT", ub_promd);
  ub["SHAKE"] = unsupported_block("CONSTRAINT", ub_g96);
  ub["GROMOS96COMPAT"] = unsupported_block("", ub_promd);
  ub["PATHINT"] = unsupported_block("", ub_promd);
  ub["NEIGHBOURLIST"] = unsupported_block("PAIRLIST", ub_promd);
  ub["PLIST"] = unsupported_block("PAIRLIST", ub_g96);
  ub["PLIST03"] = unsupported_block("PAIRLIST", ub_renamed);
  ub["NONBONDED"] = unsupported_block("LONGRANGE", ub_promd);
  ub["START"] = unsupported_block("INITIALISE", ub_g96);
  ub["OVERALLTRANSROT"] = unsupported_block("COMTRANSROT", ub_promd);
  ub["CENTREOFMASS"] = unsupported_block("COMTRANSROT", ub_g96);
  ub["POSREST"] = unsupported_block("POSITIONRES", ub_g96);
  ub["DISTREST"] = unsupported_block("DISTANCERES", ub_g96);
  ub["DIHEREST"] = unsupported_block("DIHEDRALRES", ub_g96);
  ub["J-VAL"] = unsupported_block("JVALUERES", ub_g96);
  ub["J-VAL03"] = unsupported_block("JVALUERES", ub_renamed);
  ub["PERTURB"] = unsupported_block("PERTURBATION", ub_g96);
  ub["PERTURB03"] = unsupported_block("PERTURBATION", ub_renamed);
  ub["UMBRELLA"] = unsupported_block("", ub_promd);
  ub["PRINT"] = unsupported_block("PRINTOUT", ub_g96);
  ub["WRITE"] = unsupported_block("WRITETRAJ", ub_g96);
#ifdef NDEBUG
  ub["DEBUG"] = unsupported_block("--enable-debug at compile time and "
                                  "the @verb argument", ub_promd);
#else
  ub["DEBUG"] = unsupported_block("the @verb argument", ub_promd);
#endif
  ub["FOURDIM"] = unsupported_block("", ub_g96);
  ub["LOCALELEV"] = unsupported_block("", ub_promd);
  ub["LOCALELEVATION"] = unsupported_block("LOCALELEV in PROMD", ub_g96);
  ub["SUBMOLECULES"] = unsupported_block("SOLUTEMOLECULES and moved to "
                                         "the topology", ub_renamed);
  ub["FORCEFIELD"] = unsupported_block("COVALENTFORM", ub_renamed);
  ub["PSCALE"] = unsupported_block("PERSCALE", ub_renamed);
  ub["REPLICA03"] = unsupported_block("REPLICA", ub_renamed);
  
  std::map<std::string,unsupported_block>::const_iterator
          it = ub.begin(),
          to = ub.end();
  
  // loop over unsupported blocks;
  for(; it != to; ++it) {
    // if it is present
    if (m_block[it->first].size()) {
      block_read.insert(it->first);
      
      std::ostringstream msg;
      msg << it->first << " block";
      
      switch(it->second.type) {
        case ub_renamed :
          msg << " was renamed to " << it->second.alternative;
          break;
        case ub_promd :
          msg << " is PROMD specific.";
          if (it->second.alternative != "") 
            msg << " Use " << it->second.alternative << " instead.";
          break;
        case ub_g96 :
          msg << " is GROMOS96 specific.";
          if (it->second.alternative != "") 
            msg << " Use " << it->second.alternative << " instead.";
          break;
        default : // don't know what to do.
          msg << " is known to be not supported.";
      }
      
      io::messages.add(msg.str(), "In_Parameter", io::message::error);       
    }
  }
}



