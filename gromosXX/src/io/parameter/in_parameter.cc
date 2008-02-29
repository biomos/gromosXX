
/**
 * @file in_parameter.cc
 * implements methods of In_Parameter
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

#include "in_parameter.h"

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

  read_MINIMISE(param);
  read_SYSTEM(param);
  read_CENTREOFMASS(param); // has to be read before START/INITIALISE
  read_START(param); 
  read_INITIALISE(param);
  read_STEP(param);
  read_BOUNDARY(param);
  read_REPLICA03(param); // has to be read in before MULTIBATH
  read_MULTIBATH(param);
  read_PCOUPLE(param);
  read_PRINT(param);
  read_WRITE(param);
  read_CONSTRAINT(param);
  read_FORCE(param); // and FORCEFIELD
  read_COVALENTFORM(param); // and FORCEFIELD
  read_CGRAIN(param);
  read_PAIRLIST(param);
  read_LONGRANGE(param);
  read_POSREST(param);
  read_DISTREST(param);
  read_DIHREST(param); // needs to be called after CONSTRAINT!
  read_PERTURB(param);
  read_JVALUE(param);
  read_PERSCALE(param);
  read_ROTTRANS(param);
  read_INNERLOOP(param);
  read_MULTICELL(param);
  read_READTRAJ(param);
  read_INTEGRATE(param);
  read_STOCHASTIC(param);
  read_EWARN(param);
  read_MULTISTEP(param);
  read_MONTECARLO(param);
  read_RAMD(param);
  read_CONSISTENCYCHECK(param);
  read_THERMOSTAT(param);
  read_BAROSTAT(param);
  read_VIRIAL(param);
  read_GROMOS96COMPAT(param);
  read_PATHINT(param);
  read_NEIGHBOURLIST(param);
  read_NONBONDED(param);
  read_LOCALELEVATION(param);
  read_UMBRELLA(param);
  read_FORCEFIELD(param);
  
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
 * read the MINIMISE block.
 */
void io::In_Parameter::read_MINIMISE(simulation::Parameter &param,
				     std::ostream & os)
{
  DEBUG(8, "reading MINIMISE");

  std::vector<std::string> buffer;
  buffer = m_block["MINIMISE"];
  std::string s;
  
  if (!buffer.size()){
    // no minimisation
    return;
  }

  block_read.insert("MINIMISE");

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
    io::messages.add("bad line in MINIMISE block",
		       "In_Parameter", io::message::error);
 
  // allow 0 to disable feature...
  if (param.minimise.nmin == 0)
    param.minimise.nmin = 1;

  if (param.minimise.ntem != 0 && param.minimise.ntem != 1)
    io::messages.add("MINIMISE: currently only steepest descent implemented",
		     "io::In_Parameter::read_MINIMISE",
		     io::message::error);
  if(param.minimise.ntem == 1 && param.minimise.ncyc > 0)
    io::messages.add("MINIMISE: NCYC > 0 has no effect for steepest descent",
		     "io::In_Parameter::read_MINIMISE",
		     io::message::warning);
    
  if(param.minimise.ncyc < 0)
    io::messages.add("MINIMISE: NCYC should be >0",
		     "io::In_Parameter::read_MINIMISE",
		     io::message::error);
    
  if(param.minimise.dele < 0)
    io::messages.add("MINIMISE: DELE should be >0",
		     "io::In_Parameter::read_MINIMISE",
		     io::message::error);
    
  if(param.minimise.dx0 < 0)
    io::messages.add("MINIMISE: DX0 should be >0",
		     "io::In_Parameter::read_MINIMISE",
		     io::message::error);
    
  if(param.minimise.dxm < param.minimise.dx0)
    io::messages.add("MINIMISE: DXM should be > DX0",
		     "io::In_Parameter::read_MINIMISE",
		     io::message::error);

  if(param.minimise.nmin <= 0)
    io::messages.add("MINIMISE: NMIN should be >= 0",
		     "io::In_Parameter::read_MINIMISE",
		     io::message::error);
  
  if(param.minimise.flim < 0)
    io::messages.add("MINIMISE: FLIM should be >= 0",
		     "io::In_Parameter::read_MINIMISE",
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
    io::messages.add("bad line in STEP block",
		       "In_Parameter", io::message::error);

  /*
  if (!_lineStream.eof())
    io::messages.add("End of line not reached in STEP block, but should have been: \n" + s +  "\n",
		       "In_Parameter", io::message::warning);
  */

  if(param.step.t0 < 0 && param.step.t0 != -1.0)
    io::messages.add("Negative time in STEP block is not supported",
		     "In_Parameter", io::message::error);
  if(param.step.number_of_steps <= 0)
    io::messages.add("We want to do at least one step...",
		     "In_Parameter", io::message::error);
}

/**
 * the SHAKE block.
 */
void io::In_Parameter::read_SHAKE(simulation::Parameter &param,
				  std::ostream & os)
{
  DEBUG(8, "read SHAKE");

  std::vector<std::string> buffer;
  std::string s, sntc;
  
  buffer = m_block["SHAKE"];
  
  if (buffer.size()){
    block_read.insert("SHAKE");

    io::messages.add("There is no SHAKE block in this version. Use the "
                     "CONSTRAINT block instead.", "In_Parameter",
		     io::message::error);

  }
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
      io::messages.add("NTC not understood in CONSTRAINT block",
		       "In_Parameter", io::message::error);
  }
  
  if(param.constraint.ntc<1 || param.constraint.ntc > 4){
    io::messages.add("NTC out of range in CONSTRAINT block",
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
      io::messages.add("shake tolerance in CONSTRAINT block should be > 0",
		       "In_Parameter", io::message::error);
  }
  else if (salg == "flexshake"){

    DEBUG(9, "constraints solute flexshake");
    
    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = simulation::constr_flexshake;
    else param.constraint.solute.algorithm = simulation::constr_off;
    
    _lineStream >> param.constraint.solute.shake_tolerance
		>> param.constraint.solute.flexshake_readin
		>> param.constraint.solute.flexshake_mode;
    
    if(param.constraint.solute.shake_tolerance <= 0.0)
      io::messages.add("shake tolerance in CONSTRAINT block should be > 0",
		       "In_Parameter", io::message::error);

    if(param.centreofmass.remove_rot || param.centreofmass.remove_trans)
      io::messages.add("flexible shake and removal of centre of mass motion "
		       "needs extra care!", "In_Parameter", io::message::warning);

    if(param.constraint.solute.flexshake_mode < 0 ||
       param.constraint.solute.flexshake_mode > 3)
      io::messages.add("flexshake mode in CONSTRAINT block should be >= 0 and <= 3",
		       "In_Parameter", io::message::error);
    
  }
  else if (salg == "lincs"){

    DEBUG(9, "constraints solute lincs");

    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = simulation::constr_lincs;
    else 
      param.constraint.solute.algorithm = simulation::constr_off;
    
    _lineStream >> param.constraint.solute.lincs_order;
    
    if(param.constraint.solute.lincs_order < 1)
      io::messages.add("lincs order should be >1 in CONSTRAINT block",
		       "In_Parameter", io::message::error);

  }
  else if (salg == "off"){

    DEBUG(9, "constraints solute off");
    param.constraint.solute.algorithm = simulation::constr_off;
  }
  else{

    DEBUG(9, "constraints solute error");
    
    io::messages.add("unknown algorithm in CONSTRAINT block (solute)",
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

  if (salg == "shake"){

    DEBUG(9, "constraints solvent shake");

    param.constraint.solvent.algorithm = simulation::constr_shake;
    _lineStream >> param.constraint.solvent.shake_tolerance;
    
    if(param.constraint.solvent.shake_tolerance <= 0.0)
      io::messages.add("shake tolerance in CONSTRAINT block should be > 0",
		       "In_Parameter", io::message::error);
  }
  else if (salg == "flexshake"){

    DEBUG(9, "constraints solvent flexshake");

    param.constraint.solvent.algorithm = simulation::constr_flexshake;
    _lineStream >> param.constraint.solvent.shake_tolerance;
    
    if(param.constraint.solvent.shake_tolerance <= 0.0)
      io::messages.add("shake tolerance in CONSTRAINT block should be > 0",
		       "In_Parameter", io::message::error);
  }
  else if (salg == "lincs"){

    DEBUG(9, "constraints solvent lincs");

    param.constraint.solvent.algorithm = simulation::constr_lincs;
    _lineStream >> param.constraint.solvent.lincs_order;
    
    if(param.constraint.solvent.lincs_order < 1)
      io::messages.add("lincs order should be >1 in CONSTRAINT block",
		       "In_Parameter", io::message::error);

  }
  else if (salg == "off"){

    DEBUG(9, "constraints solvent off");

    param.constraint.solvent.algorithm = simulation::constr_off;
    io::messages.add("no constraints for SOLVENT: are you sure???",
		     "In_Parameter", io::message::warning);
  }
  else{

    DEBUG(9, "constraints solvent error");
    io::messages.add("unknown algorithm in CONSTRAINT block (solvent)",
		     "In_Parameter", io::message::error);
    
    param.constraint.solvent.algorithm = simulation::constr_off;

  }
  
}

/**
 * read the PRINT
 */
void io::In_Parameter::read_PRINT(simulation::Parameter &param,
				  std::ostream & os)
{
  DEBUG(8, "read PRINT");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["PRINT"];

  if (!buffer.size()){
    io::messages.add("no PRINT block", "In_Parameter", io::message::notice);
    return;
  }

  block_read.insert("PRINT");

  int ntpp;
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  _lineStream >> param.print.stepblock
	      >> ntpp;
  
  if (_lineStream.fail())
    io::messages.add("bad line in PRINT block",
		     "In_Parameter", io::message::error);
  
  if(param.print.stepblock < 0)
    io::messages.add("Error in PRINT block: NTPR should be >=0.",
		     "In_Parameter", io::message::error);
  
  switch(ntpp) {
    case 0 : 
      param.print.monitor_dihedrals = false;
      break;
    case 1 :
      param.print.monitor_dihedrals = true;
      break;
    default:
    io::messages.add("Error in PRINT block: NTPP should be 0 or 1.",
		     "In_Parameter", io::message::error);      
  }
}

/**
 * read the WRITE
 */
void io::In_Parameter::read_WRITE(simulation::Parameter &param,
				  std::ostream & os)
{
  DEBUG(8, "read WRITE");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["WRITE"];

  if (!buffer.size()){
    io::messages.add("no WRITE block", "In_Parameter", io::message::notice);
    return;
  }

  block_read.insert("WRITE");

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
    io::messages.add("bad line in WRITE block",
		     "In_Parameter", io::message::error);
  /*
  if (!_lineStream.eof())
    io::messages.add("End of line not reached in WRITE block, "
		     "but should have been: \n" + s +  "\n",
		     "In_Parameter", io::message::warning);
  */
  if(ntwse!=0)
    io::messages.add("NTWSE != 0 not supported",
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
    io::messages.add("WRITE block: NTWE should be >= 0",
		     "In_Parameter", io::message::error);
  if(param.write.free_energy < 0)
    io::messages.add("WRITE block: NTWG should be >= 0",
		     "In_Parameter", io::message::error);
  if(param.write.block_average < 0)
    io::messages.add("WRITE block: NTWB should be >= 0",
		     "In_Parameter", io::message::error);
}

/**
 * read the PCOUPLE block.
 */
void io::In_Parameter::read_PCOUPLE(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read PRESSURESCALE");

  std::vector<std::string> buffer;
  std::string s;
  

  // first try for a PCOUPLESCALE block
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
	io::messages.add("requesting scaling but SCALE set to OFF\n",
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
	io::messages.add("bad value for SCALE switch in PRESSURESCALE block\n"
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
			 "no virial specified\n",
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
    
  } // PCOUPLE03 block
    
  buffer = m_block["PCOUPLE"];
  if (buffer.size()){
    io::messages.add("The PCOUPLE block is not available in this version. "
                     "Use the PRESSURESCALE block instead.",
                     "In_Parameter",
                     io::message::error);

    block_read.insert("PCOUPLE");
    return;
  }
    
  buffer = m_block["PCOUPLE03"];
  if (buffer.size()){
    io::messages.add("The PCOUPLE03 block was renamed to PRESSURESCALE",
                     "In_Parameter",
                     io::message::error);

    block_read.insert("PCOUPLE03");
    return;
  }
  
  if (param.pcouple.calculate==false && param.pcouple.scale!=math::pcouple_off)
    io::messages.add("pressure coupling activated but "
		     "not calculating pressure",
		     "In_Parameter",
		     io::message::error);
  if (param.pcouple.calculate == true && param.pcouple.virial == math::no_virial)
    io::messages.add("PCOUPLE03 block: pressure calculation requested but"
		     " no virial specified!", "In_Parameter",
		     io::message::error);
  if(param.pcouple.compressibility <=0)
    io::messages.add("PCOUPLE block: compressibility should be >0 ",
		     "In_Parameter", io::message::error);
  if(param.pcouple.tau <=0)
    io::messages.add("PCOUPLE block: tau should be >0 ",
		     "In_Parameter", io::message::error);
}

/**
 * read the BOUNDARY block.
 */
void io::In_Parameter::read_BOUNDARY(simulation::Parameter &param,
				     std::ostream & os)
{
  DEBUG(8, "read BOUNDARY");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["BOUNDARY"];
  
  if (!buffer.size()){
    io::messages.add("no BOUNDARY block", "In_Parameter", io::message::error);
    param.boundary.boundary = math::vacuum;
    return;
  }

  block_read.insert("BOUNDARY");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  int ntb;
  _lineStream >> ntb >> param.boundary.dof_to_subtract;
  
  if (_lineStream.fail())
    io::messages.add("bad line in BOUNDARY block",
		     "In_Parameter", io::message::error);

  if(ntb==0) param.boundary.boundary=math::vacuum;
  else if(ntb==1) param.boundary.boundary=math::rectangular;
  else if(ntb==2) param.boundary.boundary=math::triclinic;
  else if(ntb==-1) param.boundary.boundary=math::truncoct;
  else {
    std::ostringstream msg;
    msg << "wrong value for NTB in BOUNDARY block: "
        << ntb << "\nvacuum (0), rectangular (1), triclinic (2), truncoct (-1)";
    io::messages.add(msg.str(), "In_Parameter", io::message::error);
    param.boundary.boundary=math::vacuum;
  }
  
  if (param.boundary.dof_to_subtract < 0) {
    io::messages.add("Error in BOUNDARY block: NDFMIN must be >= 0.",
		     "In_Parameter", io::message::error);
    param.boundary.dof_to_subtract = 0;    
  }
  
  if (param.boundary.dof_to_subtract > 0) {
    io::messages.add("NDFMIN > 0 not implemented",
		     "In_Parameter", io::message::warning);
    param.boundary.dof_to_subtract = 0;    
  }
}

/**
 * read the PERTURB block.
 */
void io::In_Parameter::read_PERTURB(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read PERTURB");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["PERTURB03"];
  if (buffer.size()){
    block_read.insert("PERTURB03");
    
    io::messages.add("The PERTURB03 block was renamed to PERTURB.",
		     "In_Parameter", io::message::error);
  }
  
  buffer = m_block["PERTURB"];
  if (!buffer.size()) 
    return;

  block_read.insert("PERTURB");


  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  std::string b, s1, s2;
  int ntg, scale;
  int nrdgl;
  double mu, dmut, mexp;
  _lineStream >> ntg 
              >> nrdgl
              >> param.perturbation.lambda
              >> param.perturbation.dlamt
              >> mu >> dmut
              >> param.perturbation.soft_vdw
              >> param.perturbation.soft_crf
              >> param.perturbation.lambda_exponent
              >> mexp
              >> scale;
    
  if (_lineStream.fail())
    io::messages.add("bad line in PERTURB block",
                     "In_Parameter", io::message::error);
    
  switch(ntg) {
    case 0 :
      param.perturbation.perturbation = false;
      break;
    case 1 :
      param.perturbation.perturbation = true;
      break;
    case 2 :
    case 3 :
      io::messages.add("Error in PERTURB block: NTG 2/3 not implemented.",
                       "In_Parameter", io::message::error);
      break;
    default:
      io::messages.add("Error in PERTURB block: NTG must be 0 or 1.",
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
      io::messages.add("Error in PERTURB block: NRDGL must be 0 or 1.",
                       "In_Parameter", io::message::error);
  }
  
  if (param.perturbation.read_initial) {
    io::messages.add("Error in PERTURB block: NRDGL != 0 not implemented.",
                     "In_Parameter", io::message::error);    
  }
  
  if (param.perturbation.lambda < 0.0 ||
      param.perturbation.lambda > 1.0) {
    io::messages.add("Error in PERTURB block: RLAM must be 0.0 to 1.0.",
                     "In_Parameter", io::message::error);        
  }
  
  if (param.perturbation.dlamt < 0.0){
    io::messages.add("Error in PERTURB block: DLAMT must be >= 0.",
                     "In_Parameter", io::message::error);
  }
  
  if (param.perturbation.lambda_exponent<=0){
    io::messages.add("Error in PERTURB block: NLAM must be > 0.",
                     "In_Parameter", io::message::error);
  }
  
  if (mu != 0.0 || dmut != 0.0 || dmut != 0.0) {
    io::messages.add("PERTURB block: RMU and related values are ignored.",
                     "In_Parameter", io::message::warning);    
  }
  
  if (param.perturbation.soft_vdw < 0.0){
    io::messages.add("Error in PERTURB block: ALPHLJ must be >= 0.",
                     "In_Parameter", io::message::error);
  }
  
  if (param.perturbation.soft_vdw < 0.0){
    io::messages.add("Error in PERTURB block: ALPHC must be >= 0.",
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
      io::messages.add("Error in PERTURB block: NSCALE must be 0 to 2.",
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
    io::messages.add("number of energy groups in FORCE block should be > 0",
		     "In_Parameter", io::message::error);
    return;
  }
  
  for(unsigned int i=0; i<num; ++i){
    _lineStream >> e;
    DEBUG(10, "\tadding energy group " << e-1);
    param.force.energy_group.push_back(e-1);
    if(e<=old_e){
      DEBUG(10, "energy groups not in order...");
      io::messages.add("energy groups are not in order in FORCE block",
		       "In_Parameter", io::message::error);
      return;
    }
    old_e = e;
  }
  
  DEBUG(10, "number of energy groups: " << param.force.energy_group.size());

  if (_lineStream.fail())
    io::messages.add("bad line in ENERGYGROUP (FORCE) block",
		     "In_Parameter", io::message::error);
  
  
  if (bondH ^ param.force.bond)
    io::messages.add("Force switch for bond and bond H has to be equal",
		     "In_Parameter", io::message::error);

  if (angleH ^ param.force.angle)
    io::messages.add("Force switch for angle and angle H has to be equal",
		     "In_Parameter", io::message::error);

  if (impH ^ param.force.improper)
    io::messages.add("Force switch for improper and improper H has to be equal",
		     "In_Parameter", io::message::error);

  if (dihedralH ^ param.force.dihedral)
    io::messages.add("Force switch for dihedral and dihedral H has to be equal",
		     "In_Parameter", io::message::error);

  if ((!param.force.nonbonded_crf) && param.force.nonbonded_vdw)
    io::messages.add("Force: setting charges to zero",
		     "In_Parameter", io::message::notice);
  
  if (param.force.nonbonded_crf && (!param.force.nonbonded_vdw))
    io::messages.add("Force: setting atom types to dummy",
		     "In_Parameter", io::message::notice);
  
  if (_lineStream.fail())
    io::messages.add("bad line in FORCE block",
		       "In_Parameter", io::message::error);

  /*
  if (!_lineStream.eof())
    io::messages.add("End of line not reached in FORCE block, but should have been: \n" + s +  "\n",
		     "In_Parameter", io::message::warning);
  */
  
  if(param.force.bond < 0 || param.force.bond > 1)
    io::messages.add("Illegal value for force switch for bond",
		     "In_Parameter", io::message::error);
  if(param.force.angle < 0 || param.force.angle > 1)
    io::messages.add("Illegal value for force switch for angle",
		     "In_Parameter", io::message::error);
  if(param.force.improper < 0 || param.force.improper > 1)
    io::messages.add("Illegal value for force switch for improper dihedral",
		     "In_Parameter", io::message::error);
  if(param.force.dihedral < 0 || param.force.dihedral > 1)
    io::messages.add("Illegal value for force switch for dihedral",
		     "In_Parameter", io::message::error);
  if(param.force.nonbonded_vdw < 0 || param.force.nonbonded_vdw > 1)
    io::messages.add("Illegal value for force switch for nonbonded (vdw)",
		     "In_Parameter", io::message::error);
  if(param.force.nonbonded_crf < 0 || param.force.nonbonded_crf > 1)
    io::messages.add("Illegal value for force switch for nonbonded (crf)",
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
    io::messages.add("Error in COVALENTFORM block: NTBBH must be 0 (quartic) "
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
    io::messages.add("Error in COVALENTFORM block: NTBAH must be 0 (quartic) "
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
    io::messages.add("Error in COVALENTFORM block: NTBDN must be 0 (arbitray "
                     "phase shifts) or 1 (phase shifts limited).",
		     "In_Parameter", io::message::error);
  } else {
    if (dihedral != 0) {
      switch (angle) {
        case 1: 
         io::messages.add("Error in COVALENTFORM block: NTBDN 1 not implemented.",
		          "In_Parameter", io::message::error);
          break;
        case 0:
        default: ;
      }
    }
  }
}

/**
 * read START block.
 */
void io::In_Parameter::read_START(simulation::Parameter &param,
				  std::ostream & os)
{
  DEBUG(8, "read START");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["START"];
  
  if (buffer.size()){
    io::messages.add("There is no START block in this version. Use the "
                     "INITIALISE block instead.", "In_Parameter",
                     io::message::error);
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
    default : io::messages.add("Error in INITIALISE block: NTIVEL must be 0 or 1",
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
    default: io::messages.add("Error in INITIALISE block: NTISHK must be 0 to 3",
		     "In_Parameter", io::message::error);
  }
  
  // controls reading of Nose-Hoover chain variables: not implemented.
  switch(ntinhc) {
    case 0 : param.start.read_nosehoover_chains = true; break;
    case 1 : param.start.read_nosehoover_chains = false; break; // reset them
    default : io::messages.add("Error in INITIALISE block: NTINHC must be 0 or 1",
		     "In_Parameter", io::message::error);
  }
  
  io::messages.add("INITIALISE block: NTINHC is ignored.",
		   "In_Parameter", io::message::notice);
  
  // only for lattice sum. ignored without warning
  switch(ntishi) { 
    case 0: break;
    case 1: break;
    default : io::messages.add("Error in INITIALISE block: NTISHI must be 0 or 1",
		     "In_Parameter", io::message::error);
  }
  
  // controls reading of restart data for roto-translational constraints:
  // not implemented.
  switch(ntirtc) {
    case 0: param.start.read_rottrans = true; break;
    case 1: param.start.read_rottrans = false; break;
    default : io::messages.add("Error in INITIALISE block: NTIRTC must be 0 or 1",
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
    default : io::messages.add("Error in INITIALISE block: NTICOM must be 0 to 2",
		     "In_Parameter", io::message::error);
  }
  
  // controls reading of stochastic integrals: not implemented.
  switch(ntisti) {
    case 0: param.start.read_stochastic = true; break;
    case 1: param.start.read_stochastic = false; break;
    default : io::messages.add("Error in INITIALISE block: NTISTI must be 0 or 1",
		     "In_Parameter", io::message::error);
  }
  
  io::messages.add("INITIALISE block: NTISTI is ignored.",
		     "In_Parameter", io::message::notice);
  
  if(param.start.tempi <0)
    io::messages.add("Illegal value for TEMPI in START block (>=0)",
		     "In_Parameter", io::message::error);
}
/**
 * read CENTREOFMASS block.
 */
void io::In_Parameter::read_CENTREOFMASS(simulation::Parameter &param,
					 std::ostream & os)
{
  DEBUG(8, "read CENTREOFMASS");

  std::vector<std::string> buffer;
  std::string s;
  
  DEBUG(10, "reading CENTREOFMASS block");
  buffer = m_block["CENTREOFMASS"];
  
  if (!buffer.size()){
    return;
  }

  block_read.insert("CENTREOFMASS");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  int nscm, npcm, npcon;
  
  _lineStream >> nscm 
	      >> npcm 
	      >> npcon;

  if (_lineStream.fail())
    io::messages.add("bad line in CENTREOFMASS block",
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

  if (npcm != 0) {
    io::messages.add("Error in CENTREOFMASS block: NPCM != 0 not implemented.",
		     "In_Parameter", io::message::error);
  }
  
  if (npcon != 0) {
    io::messages.add("Error in CENTREOFMASS block: NPCON != 0 not implemented.",
		     "In_Parameter", io::message::error);
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
  
  /*
  if (!_lineStream.eof())
    io::messages.add("End of line not reached, but should have been: \n" + s +  "\n",
		     "In_Parameter", io::message::warning);
  */

  if(param.longrange.rf_cutoff < 0) {
    param.longrange.rf_excluded=true;
    param.longrange.rf_cutoff = -param.longrange.rf_cutoff;
    
    /**
    io::messages.add("Reaction field contribution of excluded atoms is " 
		     "taken into account", "In_Parameter", 
		     io::message::notice);
    */
  }
  else{
    param.longrange.rf_excluded=false;
    /**
    io::messages.add("Reaction field contribution of excluded atoms is " 
		     "NOT taken into account", "In_Parameter", 
		     io::message::notice);
    */
  }
  if(param.longrange.rf_epsilon!=0 && param.longrange.rf_epsilon<1)
    io::messages.add("Illegal value for EPSRF in LONGRANGE block (0  / >=1)", 
		     "In_Parameter", io::message::error);
  if(param.longrange.rf_kappa <0)
    io::messages.add("Illegal value for APPAK (who came up with this name?)"
		     " in LONGRANGE block (>=0)",
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
		       "In_Parameter",
		       io::message::error);
    }
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    std::transform(s3.begin(), s3.end(), s3.begin(), tolower);

    if (s1 == "grid") param.pairlist.grid = 1;    
    else if (s1 == "vgrid") param.pairlist.grid = 2;
    else if (s1 == "standard") param.pairlist.grid = 0;
    else{
      io::messages.add("wrong pairlist algorithm chosen (allowed: standard, "
                       "grid) in PAIRLIST block",
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
	  io::messages.add("wrong pairlist grid size chosen (allowed: auto, "
                           "[size]) in PAIRLIST block",
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
	io::messages.add("wrong cutoff type chosen (allowed: atomic, chargegroup)"
			 " in PAIRLIST block",
			 "In_Parameter", io::message::error);
	param.pairlist.atomic_cutoff = false;
      }
    }
  }

  buffer = m_block["PLIST"]; 
  if (buffer.size()){
    io::messages.add("There is no PLIST block in this version. Use the "
                     "PAIRLIST block instead.","In_Parameter",
                     io::message::error);
    block_read.insert("PLIST");
  }
  
  buffer = m_block["PLIST03"]; 
  if (buffer.size()){
    io::messages.add("The PLIST03 block was renamed to PAIRLIST","In_Parameter",
                     io::message::error);
    block_read.insert("PLIST03");
  }
  
  if(param.pairlist.grid && param.pairlist.grid_size <=0)
    io::messages.add("Illegal value for grid size in PAIRLIST block (>0)",
		     "In_Parameter", io::message::error);
  if(param.pairlist.cutoff_short < 0){
    io::messages.add("Illegal value for short range cutoff in PAIRLIST block (>0)",
		     "In_Parameter", io::message::error);
  }
  if(param.pairlist.cutoff_long < param.pairlist.cutoff_short){
    io::messages.add("Illegal value for long range cutoff in PAIRLIST block (>=RCUTP)",
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
        io::messages.add("Error in CGRAIN block: NTCGRAN must be 0 to 2.",
                         "In_Parameter", io::message::error);       
    }
    
    if (param.cgrain.EPS < 0)
      io::messages.add("Error in CGRAIN block: EPS must be >= 0.0.",
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
	io::messages.add("algorithm not understood in multibath block",
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
	io::messages.add("wrong number of Nose-Hoover chains in multibath block",
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
	io::messages.add("illegal value for temp or tau in MULTIBATH block",
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
      
      io::messages.add("Multibath: no baths but coupling groups specified",
		       "In_Parameter", io::message::error);
      num = 0;
    }

    for(int i=0; i<num; ++i){
      _lineStream >> last >> com_bath >> ir_bath;
      // let it figure out the last molecule on its own
      
      if (last < 1 || com_bath < 1 || ir_bath < 1){
	io::messages.add("bad line in MULTIBATH block: range parameter < 1",
			 "In_Parameter", io::message::error);
	if (last < 1) last = 1;
	if (com_bath < 1) com_bath = 1;
	if (ir_bath < 1) ir_bath = 1;
      }

      if (com_bath > param.multibath.multibath.size() ||
	  ir_bath > param.multibath.multibath.size()){
	io::messages.add("bad line in MULTIBATH block: ir bath or com bath index too large",
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
  else{
    // try a TCOUPLE block
    
    buffer = m_block["TCOUPLE"];
    if (buffer.size()){
      block_read.insert("TCOUPLE");

      param.multibath.found_multibath=false;
      param.multibath.found_tcouple=true;
      DEBUG(11, "TCOUPLE present");
      
      io::messages.add("The TCOUPLE block is replaced by the MULTIBATH block.",
			 "In_Parameter", io::message::warning);
      
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

      
      for(int i=0; i<3; ++i){
	_lineStream >> param.multibath.tcouple.ntt[i] 
		    >> param.multibath.tcouple.temp0[i]
		    >> param.multibath.tcouple.tau[i];
      }
      
      if (_lineStream.fail()){
	io::messages.add("bad line in TCOUPLE block",
			 "In_Parameter", io::message::error);
	return;
      }	
    }
    else{
      param.multibath.found_multibath=false;
      param.multibath.found_tcouple=false;
      DEBUG(11, "no TEMPERATURE COUPLING block");
      // no TCOUPLE block
      // that's fine, same as 0,0,0
    }
    
  }
  
} // TEMPERATURE coupling

/**
 * read POSREST block.
 */
void io::In_Parameter::read_POSREST(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read POSREST");

  std::vector<std::string> buffer;
  std::string s;
  
  DEBUG(10, "posrest block");
  buffer = m_block["POSREST"];
  
  if (!buffer.size()){
    return;
  }
  
  block_read.insert("POSREST");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  int ntpors;
  _lineStream >> param.posrest.posrest
              >> param.posrest.nrdrx
              >> ntpors
	      >> param.posrest.force_constant;
  
  if (_lineStream.fail())
    io::messages.add("bad line in POSREST block",
		     "In_Parameter", io::message::error);
  
  if (param.posrest.posrest < 0 || param.posrest.posrest > 3) {
    io::messages.add("Error in POSREST block: NTPOR must be 0 to 3.",
		     "In_Parameter", io::message::error);    
  }
  
  if (param.posrest.nrdrx < 0 || param.posrest.nrdrx > 3) {
    io::messages.add("Error in POSREST block: NTPORB must be 0 or 1.",
		     "In_Parameter", io::message::error);    
  }

  switch(ntpors) {
    case 0 :
      param.posrest.scale_reference_positions=false;
      break;
    case 1 :
      param.posrest.scale_reference_positions=true;
      break;
    default:
      param.posrest.scale_reference_positions=false;
      io::messages.add("Error in POSREST block: NTPORS must be 0 or 1.",
		       "In_Parameter", io::message::error);          
  }
  
  if(param.posrest.posrest == 3 && param.pcouple.scale != math::pcouple_off &&
     param.posrest.scale_reference_positions != true) {
    io::messages.add("Error in POSREST block: Position constraining together "
                     "with pressure coupling only allowed if reference positions"
                     " are scaled (NTPORS = 1).",
                     "In_Parameter", io::message::error);
  }
  
  if(param.posrest.posrest == 2 && param.posrest.nrdrx == false) {
    io::messages.add("Error in POSREST block: NTPOR = 2 requires NTPORB = 1",
                     "In_Parameter", io::message::error);
  }
  
  if(param.posrest.force_constant <0)
    io::messages.add("Error in POSREST block: Illegal value for CPOR.",
		     "In_Parameter", io::message::error);

} // POSREST

/**
 * read DISTREST block.
 */
void io::In_Parameter::read_DISTREST(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read DISTREST");

  std::vector<std::string> buffer;
  std::string s;
  
  DEBUG(10, "distrest block");
  buffer = m_block["DISTREST"];
  
  if (!buffer.size()){
    return;
  }
  
  block_read.insert("DISTREST");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  int ntdira;
  _lineStream >> param.distrest.distrest
  	      >> ntdira
	      >> param.distrest.K
	      >> param.distrest.r_linear
	      >> param.distrest.tau;
  
  if (_lineStream.fail())
    io::messages.add("bad line in DISTREST block",
		     "In_Parameter", io::message::error);
  
  if(param.distrest.distrest <0 || ntdira != 0 || param.distrest.tau != 0.0) {
    io::messages.add("DISTREST block: Distance restraint averaging "
                     "not implemented", "In_Parameter", io::message::notice);
  }
  
  if(param.distrest.distrest <0 || param.distrest.distrest > 2) {
    io::messages.add("Error in DISTREST block: NTDIR must 0 to 2.",
                     "In_Parameter", io::message::error);
  }
  
  if(ntdira != 0) {
    io::messages.add("Error in DISTREST block: NTDIRA must be 0.",
                     "In_Parameter", io::message::error);
  }
  
  if(param.distrest.tau != 0.0) {
    io::messages.add("Error in DISTREST block: TAUDIR must be 0.0.",
                     "In_Parameter", io::message::error);
  }

  if(param.distrest.K <0) {
    io::messages.add("Error in DISTREST block: CDIR must be >= 0.0.",
		     "In_Parameter", io::message::error);
  }

} // DISTREST

/**
 * read DIHREST block.
 */
void io::In_Parameter::read_DIHREST(simulation::Parameter &param,
				    std::ostream & os)
{
  DEBUG(8, "read DIHREST");

  std::vector<std::string> buffer;
  std::string s;
  double phi_lin;
  
  DEBUG(10, "dihrest block");
  buffer = m_block["DIHREST"];
  
  if (!buffer.size()){
    return;
  }
  
  block_read.insert("DIHREST");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  _lineStream >> param.dihrest.dihrest
	      >> param.dihrest.K
	      >> phi_lin;

  param.dihrest.phi_lin = phi_lin * 2 * math::Pi / 360;
  
  if (_lineStream.fail())
    io::messages.add("bad line in DIHREST block",
		     "In_Parameter", io::message::error);
  
  if(param.dihrest.dihrest < 0){
    io::messages.add("Dihedral restraint averaging not implemented",
		     "In_Parameter", 
		     io::message::error);
  }
  
  if(param.dihrest.dihrest > 3) {
    io::messages.add("bad input in DIHREST block, NTDHR must be <=3",
		     "In_Parameter", 
		     io::message::error);
    
  }

  if(param.distrest.K < 0)
    io::messages.add("Illegal value for force constant"
		     " in DIHREST block (>=0)",
		     "In_Parameter", io::message::error);

  if (param.dihrest.dihrest == 3){
    if (param.constraint.ntc == 1 && param.constraint.solute.algorithm == simulation::constr_off)
      param.constraint.solute.algorithm = simulation::constr_shake;

    if (param.constraint.solute.algorithm != simulation::constr_shake){
      io::messages.add("Dihedral angle constraints: needs SHAKE as (solute) constraints algorithm",
		       "In_Parameter",
		       io::message::error);
    }
  }
  
} // DIHREST

/**
 * read the JVALUE block.
 */
void io::In_Parameter::read_JVALUE(simulation::Parameter &param,
				   std::ostream & os)
{
  DEBUG(8, "read J-VAL");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["J-VAL03"];
  if (buffer.size()){

    block_read.insert("J-VAL03");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    std::string s1;
    _lineStream >> s1
		>> param.jvalue.le
		>> param.jvalue.ngrid
		>> param.jvalue.K
		>> param.jvalue.tau
		>> param.jvalue.delta
		>> param.jvalue.read_av;
    
    if (_lineStream.fail())
      io::messages.add("bad line in J-VAL03 block",
		       "In_Parameter", io::message::error);
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    
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
	io::messages.add("bad value for MODE in J-VAL03 block:"+s1+"\n"
			 "off, instantaneous, averaged, biquadratic, scaled (0-4)",
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
      io::messages.add("bad value for TAU in J-VAL03 block\n"
		       "should be > 0.0",
		       "In_Parameter", io::message::error);
    }
    if (param.jvalue.mode != simulation::restr_off && param.jvalue.K < 0.0){
      io::messages.add("bad value for K in J-VAL03 block\n"
		       "should be > 0.0",
		       "In_Parameter", io::message::error);
    }
    if (param.jvalue.le > 0){
      
      if (param.jvalue.ngrid < 1){
	io::messages.add("bad value for NGRID in J-VAL03 block\n"
			 "should be > 1",
			 "In_Parameter", io::message::error);
      }
    }

  } // J-VAL03

  else{

    buffer = m_block["J-VAL"];
    if (buffer.size()){
      
      block_read.insert("J-VAL");
      
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
      
      std::string s1;
      _lineStream >> s1 
		  >> param.jvalue.tau
		  >> param.jvalue.read_av;
      
      if (_lineStream.fail())
	io::messages.add("bad line in J-VAL block",
			 "In_Parameter", io::message::error);
      
      std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
      
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
	  io::messages.add("bad value for MODE in J-VAL block:"+s1+"\n"
			   "off, instantaneous, averaged, biquadratic, scaled (0-4)",
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
	io::messages.add("bad value for TAU in J-VAL block\n"
			 "should be > 0.0",
			 "In_Parameter", io::message::error);
      }
    }
  } // J-VAL
  
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
  
  buffer = m_block["PSCALE"];
  if (buffer.size()) {
    block_read.insert("PSCALE"); 
    io::messages.add("The PSCALE block was renamed to PERSCALE.", "In_Parameter",
                     io::message::error);
  }
  
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
      io::messages.add("Error in PERSCALE block: RESTYPE must be jrest or off.",
		       "In_Parameter", io::message::error);     
    }

    if (param.pscale.KDIH < 0.0)
      io::messages.add("Error in PERSCALE block: KDIH must be >= 0.0.",
                       "In_Parameter", io::message::error);
    if (param.pscale.KJ < 0.0)
      io::messages.add("Error in PERSCALE block: KJ must be >= 0.0.",
                       "In_Parameter", io::message::error);
    if (param.pscale.T < 0.0)
      io::messages.add("Error in PERSCALE block: T must be >= 0.0.",
                       "In_Parameter", io::message::error);
    if (param.pscale.difference < 0.0)
      io::messages.add("Error in PERSCALE block: DIFF must be >= 0.0.",
                       "In_Parameter", io::message::error);
    
    switch(read) {
      case 0 :
         param.pscale.read_data = false;
         break;
      case 1 :
         param.pscale.read_data = true;
         break;
      default :
      io::messages.add("Error in PERSCALE block: READ must be 0 or 1.",
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
      io::messages.add("Error in ROTTRANS block: RTC must be 0 or 1",
                       "In_Parameter", io::message::error);
    }

    if (param.rottrans.last <= 0)
      io::messages.add("Error in ROTTRANS block: RTCLAST must be > 0.",
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
    // to be taken out once spc loops are fully implemented and tested
    io::messages.add("special solvent loops not implemented (INNERLOOP block)",
            "In_Parameter", io::message::error);
        
    block_read.insert("INNERLOOP");
    // uncomment the following once spc loops are fully implemented and tested
    /*
    int spc;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
 
    _lineStream >> spc;
 
    if (_lineStream.fail())
      io::messages.add("bad line in INNERLOOP block",
                       "In_Parameter", io::message::error);
 
    if (spc != -1 && spc != 0){
      io::messages.add("bad value for SPCL in INNERLOOP: allowed : -1, 0",
                       "In_Parameter",
                       io::message::error);
      spc = -1;
    }
    param.force.spc_loop = spc;
     */
  }
}

/**
 * read the REPLICA03 block.
 */
void io::In_Parameter::read_REPLICA03(simulation::Parameter &param,
				      std::ostream & os)
{
  DEBUG(8, "read REPLICA03");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["REPLICA03"];

  if (buffer.size()){

    block_read.insert("REPLICA03");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.replica.num_T;
    param.replica.temperature.resize(param.replica.num_T, 0.0);

    for(int i=0; i<param.replica.num_T; ++i){
      _lineStream >> param.replica.temperature[i];
    }
    
    _lineStream >> param.replica.scale;

    if (_lineStream.fail()){
      io::messages.add("bad line in REPLICA03 block (numT, T or scale)",
		       "In_Parameter", io::message::error);
      param.replica.num_T = 0;
      param.replica.num_l = 0;

      param.replica.temperature.clear();
      param.replica.lambda.clear();
	  param.replica.dt.clear();
    }
    
    _lineStream >> param.replica.num_l;
    param.replica.lambda.resize(param.replica.num_l, 0.0);
    param.replica.dt.resize(param.replica.num_l, 0.0);
    
    for(int i=0; i<param.replica.num_l; ++i){
      _lineStream >> param.replica.lambda[i];
    }
    for(int i=0; i<param.replica.num_l; ++i){
      _lineStream >> param.replica.dt[i];
    }

    if (_lineStream.fail()){
      io::messages.add("bad line in REPLICA03 block (numl, l or dt)",
		       "In_Parameter", io::message::error);
      param.replica.num_T = 0;
      param.replica.num_l = 0;

      param.replica.temperature.clear();
      param.replica.lambda.clear();
	  param.replica.dt.clear();
    }
    _lineStream >> param.replica.trials;
    _lineStream >> param.replica.equilibrate;
    _lineStream >> param.replica.slave_runs;
    _lineStream >> param.replica.write;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in REPLICA03 block (trials, equi, slave or write)",
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
        io::messages.add("Error in MULTICELL block: NTM must be 0 or 1.",
                         "In_Parameter", io::message::error);
    }
    
    if (param.multicell.multicell) {
      // disable because broken
      param.multicell.multicell = false;
      io::messages.add("MULTICELL simulations are broken in md++",
                         "In_Parameter", io::message::error);      
      
      
      // do these checks only if mutlicell is really used.
      if (param.multicell.x < 1 || param.multicell.y < 1 ||
          param.multicell.z < 1) {
        io::messages.add("Error in MULTICELL block: NCELLA, NCELLB and NCELLC "
                         "must be >= 1.", "In_Parameter", io::message::error);
      }
    
      if (param.multicell.x == 1 && param.multicell.y == 1 && 
          param.multicell.z == 1) {
        io::messages.add("MULTICELL block: NCELLA, NCELLB and NCELLC are all 1.\n"
                         "disabling MULTICELL simulation.", "In_Parameter",
                         io::message::warning);    
        param.multicell.multicell = false;
      }
    
      if (param.boundary.boundary != math::rectangular &&
          param.boundary.boundary != math::triclinic) {
        io::messages.add("MULTICELL is only available for rectangular or "
                         "triclinic periodic boundary conditions.", "In_Parameter",
                         io::message::error);   
      }
    
      if (tolpx || tolpv || tolpf || tolpfw) {
        io::messages.add("MULTICELL block: Periodicity checks not available in "
                         "this version.\n"
                         "disabling MULTICELL simulation.", "In_Parameter",
                         io::message::warning);    
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
        io::messages.add("Error in READTRAJ block: NTRD must be 0 or 1",
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
        io::messages.add("Error in READTRAJ block: NTSHK must be 0 or 1",
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
        io::messages.add("Error in INTEGRATE block: NINT must be 0 or 1",
                         "In_Parameter", io::message::error);        
    }
  }
}

void io::In_Parameter::read_STOCHASTIC(simulation::Parameter & param,
				       std::ostream & os)
{
  DEBUG(8, "read STOCHASTIC");

  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["STOCHASTIC"];
  
  if (buffer.size()){

    block_read.insert("STOCHASTIC");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.stochastic.sd >> param.stochastic.ntfr
                >> param.stochastic.nsfr >> param.stochastic.nbref
		>> param.stochastic.rcutf >> param.stochastic.cfric
                >> param.stochastic.temp;
    
    if (_lineStream.fail()){
      io::messages.add("bad line in STOCHASTIC block",
		       "In_Parameter", io::message::error);
      return;
    }
    
    if(param.stochastic.sd < 0 || param.stochastic.sd > 1)
      io::messages.add("Error in STOCHASTIC block: NTSD must be 0 or 1",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.ntfr < 0 || param.stochastic.ntfr > 3)
      io::messages.add("Error in STOCHASTIC block: NTFR must be 0 to 3",
                       "In_Parameter", io::message::error);
    
    io::messages.add("Behaviour of NTFR not clear.", "In_Parameter", 
                     io::message::warning);
    
    if(param.stochastic.nsfr <= 0)
      io::messages.add("Error in STOCHASTIC block: NSFR must be > 0",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.nbref <= 0)
      io::messages.add("Error in STOCHASTIC block: NBREF must be > 0",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.rcutf < 0)
      io::messages.add("Error in STOCHASTIC block: RCUTF must be >= 0",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.cfric < 0)
      io::messages.add("Error in STOCHASTIC block: CFRIC must be >= 0",
                       "In_Parameter", io::message::error);
    
    if(param.stochastic.temp < 0)
      io::messages.add("Error in STOCHASTIC block: TEMPSD must be >= 0",
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
      io::messages.add("Error in MULTISTEP block: STEPS must be >= 0.",
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
        io::messages.add("Error in MULTISTEP block: BOOST must be 0 or 1",
                         "In_Parameter", io::message::error);
    }
  }
}

void io::In_Parameter::read_MONTECARLO(simulation::Parameter & param,
				       std::ostream & os)
{
  DEBUG(8, "read MONTECARLO");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["MONTECARLO"];
  
  if (buffer.size()){
    block_read.insert("MONTECARLO");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.montecarlo.mc >> param.montecarlo.steps;

    if (_lineStream.fail()){
      io::messages.add("bad line in MONTECARLO block",
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
      io::messages.add("RAMD: FC should be >=0",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);
    if(param.ramd.steps < 0)
      io::messages.add("RAMD: STEPS should be >0",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);
    if(param.ramd.r_min < 0)
      io::messages.add("RAMD: R_MIN should be >0",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);
    if(param.ramd.every < 0)
      io::messages.add("RAMD: NWRITE should be >=0",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);
    if(param.ramd.tau ==0)
      param.ramd.do_ta=false;
    else if(param.ramd.tau < 0)
      io::messages.add("RAMD: TAU should be >=0",
		      "io::In_Parameter::read_RAMD",
		      io::message::error);
    else
      param.ramd.do_ta=true;
    if(param.ramd.do_ta && param.ramd.ta_min <= 0)
      io::messages.add("RAMD: TA_MIN should be >0 if time averaging is included",
		       "io::In_Parameter::read_RAMD",
		       io::message::error);

    if(param.ramd.fc!=0.0 && param.ramd.atom.size()==0)
      io::messages.add("RAMD: no atoms read in to apply random force",
		       "io::In_Parameter::read_RAMD",
		       io::message::warning);


  }
}

void io::In_Parameter::read_CONSISTENCYCHECK(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read CONSISTENCYCHECK");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["CONSISTENCYCHECK"];
  
  if (buffer.size()) {
    block_read.insert("CONSISTENCYCHECK");
    io::messages.add("CONSISTENCYCHECK is not supported in this version. "
                     "Use \"make check\" instead.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_THERMOSTAT(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read THERMOSTAT");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["THERMOSTAT"];
  
  if (buffer.size()) {
    block_read.insert("THERMOSTAT");
    io::messages.add("The THERMOSTAT block is not supported in this version. "
                     "Use the MULTIBATH block instead.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_BAROSTAT(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read BAROSTAT");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["BAROSTAT"];
  
  if (buffer.size()) {
    block_read.insert("BAROSTAT");
    io::messages.add("The BAROSTAT block is not supported in this version. "
                     "Use the PRESSURESCALE block instead.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_VIRIAL(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read VIRIAL");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["VIRIAL"];
  
  if (buffer.size()) {
    block_read.insert("VIRIAL");
    io::messages.add("The VIRIAL block is not supported in this version. "
                     "Use the PRESSURESCALE block instead.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_NEIGHBOURLIST(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read NEIGHBOURLIST");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["NEIGHBOURLIST"];
  
  if (buffer.size()) {
    block_read.insert("NEIGHBOURLIST");
    io::messages.add("The NEIGHBOURLIST block is not supported in this version. "
                     "Use the PAIRLIST block instead.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_NONBONDED(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read NONBONDED");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["NONBONDED"];
  
  if (buffer.size()) {
    block_read.insert("NONBONDED");
    io::messages.add("The NONBONDED block is not supported in this version. "
                     "Use the LONGRANGE block instead.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_GROMOS96COMPAT(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read GROMOS96COMPAT");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["GROMOS96COMPAT"];
  
  if (buffer.size()) {
    block_read.insert("GROMOS96COMPAT");
    io::messages.add("The GROMOS96COMPAT block is not supported in this version. "
                     "Use promd for this feature.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_PATHINT(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read PATHINT");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["PATHINT"];
  
  if (buffer.size()) {
    block_read.insert("PATHINT");
    io::messages.add("The PATHINT block is not supported in this version. "
                     "Use promd for this feature.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_LOCALELEVATION(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read LOCALELEVATION");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["LOCALELEVATION"];
  
  if (buffer.size()) {
    block_read.insert("LOCALELEVATION");
    io::messages.add("The LOCALELEVATION block is not supported in this version. "
                     "Use promd for this feature.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_UMBRELLA(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read UMBRELLA");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["UMBRELLA"];
  
  if (buffer.size()) {
    block_read.insert("UMBRELLA");
    io::messages.add("The UMBRELLA block is not supported in this version. "
                     "Use promd for this feature.", "In_Parameter",
                     io::message::error);
  }
}

void io::In_Parameter::read_FORCEFIELD(simulation::Parameter & param,
				 std::ostream & os)
{
  DEBUG(8, "read FORCEFIELD");
  
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["FORCEFIELD"];
  
  if (buffer.size()) {
    block_read.insert("FORCEFIELD");
    io::messages.add("The FORCEFIELD block is not supported in this version. "
                     "Use the COVALENTFORM block instead.", "In_Parameter",
                     io::message::error);
  }
}
