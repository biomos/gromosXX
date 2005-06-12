
/**
 * @file in_parameter.cc
 * implements methods of In_Parameter
 */

#include "../stdheader.h"

#include "../simulation/parameter.h"

#include "instream.h"
#include "blockinput.h"

#include "in_parameter.h"

static std::set<std::string> block_read;

/**
 * Store standard parameters in the Parameter
 */
void In_Parameter::read(Parameter &param,
			std::ostream & os)
{
  if (!quiet)
    os << "\nINPUT\n"
	      << title << "\n";

  // store the title...
  param.title = title;

  read_MINIMISE(param);
  read_SYSTEM(param);
  read_START(param); // and CENTREOFMASS
  read_STEP(param);
  read_BOUNDARY(param);
  read_SUBMOLECULES(param);
  read_REPLICA(param); // has to be read in before MULTIBATH (to overwrite temps)
  read_REPLICA03(param);
  read_MULTIBATH(param);
  read_PCOUPLE(param);
  read_PRINT(param);
  read_WRITE(param);
  read_CONSTRAINTS(param); // read_SHAKE if no CONSTRAINTS
  read_FORCE(param); // and FORCEFIELD
  read_PLIST(param);
  read_LONGRANGE(param);
  read_POSREST(param);

  read_DISTREST(param);
  
  read_PERTURB(param);
  read_JVALUE(param);
  read_PSCALE(param);
  read_ROTTRANS(param);
  read_INNERLOOP(param);
  read_MULTICELL(param);
  read_ANALYZE(param);
  
  for(std::map<std::string, std::vector<std::string> >::const_iterator
	it = m_block.begin(),
	to = m_block.end();
      it != to;
      ++it){

    if (block_read.count(it->first) == 0 && it->second.size()){
      messages.add("block " + it->first + " not supported!",
		   "In_Parameter",
		   Message::warning);
    }
  }

  if (!quiet)
    os << "END\n";

}

/**
 * read the SYSTEM block.
 */
void In_Parameter::read_SYSTEM(Parameter &param,
				   std::ostream & os)
{
  std::vector<std::string> buffer;
  buffer = m_block["SYSTEM"];
  std::string s;
  
  if (!buffer.size()){
    messages.add("no SYSTEM block in input", "In_Parameter", Message::error);
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
    messages.add("bad line in SYSTEM block",
		       "In_Parameter", Message::error);

  // we might need to also allow for 0...
  if (param.system.npm != 1 && param.system.npm != 0)
    messages.add("SYSTEM: currently only NPM=1 allowed (NPM=0 experimental)",
		     "In_Parameter::read_SYSTEM",
		     Message::error);
  if(param.system.nsm < 0)
    messages.add("SYSTEM: NSM should be >0",
		     "In_Parameter::read_SYSTEM",
		     Message::error);

} 

/**
 * read the MINIMISE block.
 */
void In_Parameter::read_MINIMISE(Parameter &param,
				     std::ostream & os)
{
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
	      >> param.minimise.dxm;
  
  if (_lineStream.fail())
    messages.add("bad line in MINIMISE block",
		       "In_Parameter", Message::error);

  if (!(_lineStream >> param.minimise.nmin)){
    // support standard input format...
    param.minimise.nmin = 1;
    param.minimise.flim = 0.0;
    _lineStream.clear();
  }
  else if (!(_lineStream >> param.minimise.flim)){
    // or only number of steps...
    param.minimise.flim = 0.0;
    _lineStream.clear();
  }
 
  // allow 0 to disable feature...
  if (param.minimise.nmin == 0)
    param.minimise.nmin = 1;

  if (param.minimise.ntem != 0 && param.minimise.ntem != 1)
    messages.add("MINIMISE: currently only steepest descent implemented",
		     "In_Parameter::read_MINIMISE",
		     Message::error);
  if(param.minimise.ncyc < 0)
    messages.add("MINIMISE: NCYC should be >0",
		     "In_Parameter::read_MINIMISE",
		     Message::error);
  if(param.minimise.dele < 0)
    messages.add("MINIMISE: DELE should be >0",
		     "In_Parameter::read_MINIMISE",
		     Message::error);
  if(param.minimise.dx0 < 0)
    messages.add("MINIMISE: DX0 should be >0",
		     "In_Parameter::read_MINIMISE",
		     Message::error);
  if(param.minimise.dxm < param.minimise.dx0)
    messages.add("MINIMISE: DXM should be > DX0",
		     "In_Parameter::read_MINIMISE",
		     Message::error);

  if(param.minimise.nmin <= 0)
    messages.add("MINIMISE: NMIN should be >= 0",
		     "In_Parameter::read_MINIMISE",
		     Message::error);
  
  if(param.minimise.flim < 0)
    messages.add("MINIMISE: FLIM should be >= 0",
		     "In_Parameter::read_MINIMISE",
		     Message::error);

} 

/**
 * read the STEP block.
 */
void In_Parameter::read_STEP(Parameter &param,
				 std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["STEP"];
  
  if (!buffer.size()){
    messages.add("no STEP block in input", "In_Parameter", Message::error);
    return;
  }

  block_read.insert("STEP");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  _lineStream >> param.step.number_of_steps
	      >> param.step.t0
	      >> param.step.dt;
  
  if (_lineStream.fail())
    messages.add("bad line in STEP block",
		       "In_Parameter", Message::error);

  if(param.step.t0 < 0 && param.step.t0 != -1.0)
    messages.add("Negative time in STEP block is not supported",
		     "In_Parameter", Message::error);
  if(param.step.number_of_steps <= 0)
    messages.add("We want to do at least one step...",
		     "In_Parameter", Message::error);
}

/**
 * the SHAKE block.
 */
void In_Parameter::read_SHAKE(Parameter &param,
				  std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s, sntc;
  
  buffer = m_block["SHAKE"];

  if (!buffer.size()){
    param.constraint.ntc = 1;
    param.constraint.solute.algorithm = constr_off;
    param.constraint.solvent.algorithm = constr_shake;

    messages.add("no SHAKE / CONSTRAINTS block", "In_Parameter",
		     Message::warning);

    return;
  }
  
  block_read.insert("SHAKE");
  
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  _lineStream >> sntc
	      >> param.constraint.solvent.shake_tolerance;
  
  if (_lineStream.fail())
    messages.add("bad line in SHAKE block",
		       "In_Parameter", Message::error);

  std::transform(sntc.begin(), sntc.end(), sntc.begin(), tolower);

  if(sntc=="solvent") param.constraint.ntc=1;
  else if(sntc=="hydrogen") param.constraint.ntc=2;
  else if(sntc=="all") param.constraint.ntc=3;
  else if(sntc=="specified") param.constraint.ntc=4;
  else {
    std::stringstream ss(sntc);
    if (!(ss >> param.constraint.ntc)){
      messages.add("NTC not understood in CONSTRAINTS block",
		       "In_Parameter", Message::error);

      param.constraint.solute.algorithm = constr_off;
      param.constraint.solvent.algorithm = constr_off;
      
      return;
    }
  }

  if(param.constraint.ntc<1 || param.constraint.ntc > 4){
    messages.add("NTC out of range in CONSTRAINTS block",
		     "In_Parameter", Message::error);
    
    param.constraint.solute.algorithm = constr_off;
    param.constraint.solvent.algorithm = constr_off;
    
    return;
  }
  

  if(param.constraint.solvent.shake_tolerance<=0.0)
    messages.add("tolerance in SHAKE block should be > 0",
		       "In_Parameter", Message::error);

  if (param.constraint.ntc > 1){

    param.constraint.solute.algorithm = constr_shake;
    param.constraint.solute.shake_tolerance = param.constraint.solvent.shake_tolerance;
  }
  else
    param.constraint.solute.algorithm = constr_off;
  
  param.constraint.solvent.algorithm = constr_shake;

}

/**
 * the CONSTRAINTS block.
 */
void In_Parameter::read_CONSTRAINTS(Parameter &param,
					std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s, salg;
  
  buffer = m_block["CONSTRAINTS"];

  if (!buffer.size()){
    // try reading a shake block
    read_SHAKE(param);
    return;
  }
  
  block_read.insert("CONSTRAINTS");
  
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
      messages.add("NTC not understood in CONSTRAINTS block",
		       "In_Parameter", Message::error);
  }
  
  if(param.constraint.ntc<1 || param.constraint.ntc > 4){
    messages.add("NTC out of range in CONSTRAINTS block",
		     "In_Parameter", Message::error);
    
    param.constraint.solute.algorithm = constr_off;
    param.constraint.solvent.algorithm = constr_off;
    
    return;
  }
  
  // SOLUTE
  _lineStream >> salg;
  
  if (_lineStream.fail())
    messages.add("bad line in CONSTRAINTS block",
		     "In_Parameter", Message::error);

  std::transform(salg.begin(), salg.end(), salg.begin(), tolower);

  if (salg == "shake"){

    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = constr_shake;
    else param.constraint.solute.algorithm = constr_off;
    
    _lineStream >> param.constraint.solute.shake_tolerance;
    
    if(param.constraint.solute.shake_tolerance <= 0.0)
      messages.add("shake tolerance in CONSTRAINTS block should be > 0",
		       "In_Parameter", Message::error);
  }
  else if (salg == "flexshake"){

    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = constr_flexshake;
    else param.constraint.solute.algorithm = constr_off;
    
    _lineStream >> param.constraint.solute.shake_tolerance
		>> param.constraint.solute.flexshake_readin
		>> param.constraint.solute.flexshake_mode;
    
    if(param.constraint.solute.shake_tolerance <= 0.0)
      messages.add("shake tolerance in CONSTRAINTS block should be > 0",
		       "In_Parameter", Message::error);

    if(param.centreofmass.remove_rot || param.centreofmass.remove_trans)
      messages.add("flexible shake and removal of centre of mass motion "
		       "needs extra care!", "In_Parameter", Message::warning);

    if(param.constraint.solute.flexshake_mode < 0 ||
       param.constraint.solute.flexshake_mode > 3)
      messages.add("flexshake mode in CONSTRAINTS block should be >= 0 and <= 3",
		       "In_Parameter", Message::error);
    
  }
  else if (salg == "lincs"){

    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = constr_lincs;
    else 
      param.constraint.solute.algorithm = constr_off;
    
    _lineStream >> param.constraint.solute.lincs_order;
    
    if(param.constraint.solute.lincs_order < 1)
      messages.add("lincs order should be >1 in CONSTRAINTS block",
		       "In_Parameter", Message::error);

  }
  else if (salg == "off"){

    param.constraint.solute.algorithm = constr_off;
  }
  else{

    messages.add("unknown algorithm in CONSTRAINTS block (solute)",
		     "In_Parameter", Message::error);
    
    param.constraint.solute.algorithm = constr_off;

  }

  // SOLVENT
  _lineStream >> salg;
  
  if (_lineStream.fail())
    messages.add("bad line in CONSTRAINTS block",
		     "In_Parameter", Message::error);

  std::transform(salg.begin(), salg.end(), salg.begin(), tolower);

  if (salg == "shake"){

    param.constraint.solvent.algorithm = constr_shake;
    _lineStream >> param.constraint.solvent.shake_tolerance;
    
    if(param.constraint.solvent.shake_tolerance <= 0.0)
      messages.add("shake tolerance in CONSTRAINTS block should be > 0",
		       "In_Parameter", Message::error);
  }
  else if (salg == "flexshake"){

    param.constraint.solvent.algorithm = constr_flexshake;
    _lineStream >> param.constraint.solvent.shake_tolerance;
    
    if(param.constraint.solvent.shake_tolerance <= 0.0)
      messages.add("shake tolerance in CONSTRAINTS block should be > 0",
		       "In_Parameter", Message::error);
  }
  else if (salg == "lincs"){

    param.constraint.solvent.algorithm = constr_lincs;
    _lineStream >> param.constraint.solvent.lincs_order;
    
    if(param.constraint.solvent.lincs_order < 1)
      messages.add("lincs order should be >1 in CONSTRAINTS block",
		       "In_Parameter", Message::error);

  }
  else if (salg == "off"){

    param.constraint.solvent.algorithm = constr_off;
    messages.add("no constraints for SOLVENT: are you sure???",
		     "In_Parameter", Message::warning);
  }
  else{

    messages.add("unknown algorithm in CONSTRAINTS block (solvent)",
		     "In_Parameter", Message::error);
    
    param.constraint.solvent.algorithm = constr_off;

  }
  
}

/**
 * read the PRINT
 */
void In_Parameter::read_PRINT(Parameter &param,
				  std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["PRINT"];

  if (!buffer.size()){
    messages.add("no PRINT block", "In_Parameter", Message::notice);
    return;
  }

  block_read.insert("PRINT");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  _lineStream >> param.print.stepblock
	      >> param.print.centreofmass
	      >> param.print.monitor_dihedrals;
  
  if (_lineStream.fail())
    messages.add("bad line in PRINT block",
		     "In_Parameter", Message::error);

  if(param.print.stepblock<=0)
    messages.add("PRINT block: print stepblock should be >0",
		     "In_Parameter", Message::error);
  if(param.print.centreofmass < 0)
    messages.add("PRINT block: print centre of mass should be >0",
		     "In_Parameter", Message::error);
}

/**
 * read the WRITE
 */
void In_Parameter::read_WRITE(Parameter &param,
				  std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["WRITE"];

  if (!buffer.size()){
    messages.add("no WRITE block", "In_Parameter", Message::notice);
    return;
  }

  block_read.insert("WRITE");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  int ntwse, ntpw;
  
  _lineStream >> param.write.position
	      >> ntwse
	      >> param.write.velocity
	      >> param.write.energy
	      >> param.write.free_energy
	      >> param.write.block_average
	      >> ntpw;
  
  if (_lineStream.fail())
    messages.add("bad line in WRITE block",
		     "In_Parameter", Message::error);
  if(ntwse!=0)
    messages.add("NTWSE != 0 not supported",
		     "In_Parameter", Message::error);
  if(ntpw!=1)
    messages.add("NTPW != 1 not supported",
		     "In_Parameter", Message::error);

  if(param.write.position < 0){
    param.write.solute_only = true;
    param.write.position = -param.write.position;
    messages.add("writing solute only trajectory",
		     "In_Parameter", Message::notice);
  }
  
  if(param.write.velocity < 0)
    messages.add("WRITE block: NTWV should be >= 0",
		     "In_Parameter", Message::error);
  if(param.write.energy < 0)
    messages.add("WRITE block: NTWE should be >= 0",
		     "In_Parameter", Message::error);
  if(param.write.free_energy < 0)
    messages.add("WRITE block: NTWG should be >= 0",
		     "In_Parameter", Message::error);
  if(param.write.block_average < 0)
    messages.add("WRITE block: NTWB should be >= 0",
		     "In_Parameter", Message::error);
}

/**
 * read the PCOUPLE block.
 */
void In_Parameter::read_PCOUPLE(Parameter &param,
				    std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  

  // first try for a PCOUPLE03 block
  buffer = m_block["PCOUPLE03"];

  if (buffer.size()){

    block_read.insert("PCOUPLE03");

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
      messages.add("bad line in PCOUPLE03 block",
		       "In_Parameter", Message::error);
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    std::transform(s3.begin(), s3.end(), s3.begin(), tolower);

    if (s1 == "off"){
      param.pcouple.calculate = false; param.pcouple.scale = pcouple_off;
    } 
    else if (s1 == "calc"){
      param.pcouple.calculate = true;
      param.pcouple.scale = pcouple_off;
    }
    else if (s1 == "scale"){
      param.pcouple.calculate = true;

      if (s2 == "off"){
	messages.add("requesting scaling but SCALE set to OFF\n",
			 "In_Parameter", Message::error);
	param.pcouple.scale = pcouple_off;
      }
      else if (s2 == "iso")
	param.pcouple.scale = pcouple_isotropic;
      else if (s2 == "aniso")
	param.pcouple.scale = pcouple_anisotropic;
      else if (s2 == "full")
	param.pcouple.scale = pcouple_full_anisotropic;
      else{
	messages.add("bad value for SCALE switch in PCOUPLE03 block\n"
			 "(off,iso,aniso,full)",
			 "In_Parameter", Message::error);
	param.pcouple.scale = pcouple_off;
      }

    }
    else{
      messages.add("bad value for calc switch in PCOUPLE03 block\n"
		       "(off,calc,scale)",
		       "In_Parameter", Message::error);
      param.pcouple.calculate = false;
    }
  
    if (param.pcouple.calculate){
      if (s3 == "none"){
	messages.add("requesting pressure calculation but "
			 "no virial specified\n",
			 "In_Parameter", Message::error);
	param.pcouple.virial = no_virial;
      }
      else if (s3 == "atomic")
	param.pcouple.virial = atomic_virial;
      else if (s3 == "molecular")
	param.pcouple.virial = molecular_virial;
      else{
	messages.add("bad value for virial switch in PCOUPLE03 block\n"
			 "(none,atomic,molecular)",
			 "In_Parameter", Message::error);
	param.pcouple.virial = no_virial;
      }
    }
    else
      param.pcouple.virial = no_virial;
    
  } // PCOUPLE03 block
  else{
    
    buffer = m_block["PCOUPLE"];
    if (!buffer.size()){
      // everything should be initialized to nothing
      return;
    }

    block_read.insert("PCOUPLE");
 
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    int ntp;    
    double p0;
    
    _lineStream >> ntp >> p0 
		>> param.pcouple.compressibility 
		>> param.pcouple.tau;
    
    if (_lineStream.fail())
      messages.add("bad line in PCOUPLE block",
		       "In_Parameter", Message::error);

    for(int i=0; i<3; ++i){
      for(int j=0; j<3; ++j){
	if (i == j)
	  param.pcouple.pres0(i,j) = p0;
	else param.pcouple.pres0(i,j) = 0.0;
      }
    }
    if(ntp==1){
      param.pcouple.scale = pcouple_isotropic;
      param.pcouple.virial = molecular_virial;
    }
    else if(ntp==2){
      param.pcouple.scale = pcouple_anisotropic;
      param.pcouple.virial = molecular_virial;
    }
    else if(ntp>0 || ntp <0)
      messages.add("PCOUPLE block: illegal value for ntp (0,1,2)",
		       "In_Parameter",
		       Message::error);
  }
  if (param.pcouple.calculate==false && param.pcouple.scale!=pcouple_off)
    messages.add("pressure coupling activated but "
		     "not calculating pressure",
		     "In_Parameter",
		     Message::error);
  if (param.pcouple.calculate == true && param.pcouple.virial == no_virial)
    messages.add("PCOUPLE03 block: pressure calculation requested but"
		     " no virial specified!", "In_Parameter",
		     Message::error);
  if(param.pcouple.compressibility <=0)
    messages.add("PCOUPLE block: compressibility should be >0 ",
		     "In_Parameter", Message::error);
  if(param.pcouple.tau <=0)
    messages.add("PCOUPLE block: tau should be >0 ",
		     "In_Parameter", Message::error);
}

/**
 * read the BOUNDARY block.
 */
void In_Parameter::read_BOUNDARY(Parameter &param,
				     std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["BOUNDARY"];
  
  if (!buffer.size()){
    messages.add("no BOUNDARY block", "In_Parameter", Message::error);
    param.boundary.boundary = vacuum;
    return;
  }

  block_read.insert("BOUNDARY");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  std::string ntb;
  int n, nrdbox;
  double b1, b2, b3, beta;
  std::istringstream cs;
  
  _lineStream >> ntb >> b1 >> b2 >> b3 >> beta >> nrdbox;
  
  if (_lineStream.fail())
    messages.add("bad line in BOUNDARY block",
		     "In_Parameter", Message::error);

  std::transform(ntb.begin(), ntb.end(), ntb.begin(), tolower);
  if(ntb=="vacuum") param.boundary.boundary=vacuum;
  else if(ntb=="rectangular") param.boundary.boundary=rectangular;
  else if(ntb=="triclinic") param.boundary.boundary=triclinic;
  else if(ntb=="truncoct") param.boundary.boundary=truncoct;
  else {
    cs.str(ntb);
    // n=atoi(ntb.c_str());
    cs >> n;
    if(cs.fail()){
      std::cerr << "boundary error number = " << n << std::endl;
      messages.add("wrong value for NTB in BOUNDARY block: "+ntb+"\n"
		       "vacuum, rectangular, triclinic, 0, +/-1, +/-2",
		       "In_Parameter", Message::error);
      param.boundary.boundary=vacuum;
      return;
    }
    if(n==0) param.boundary.boundary=vacuum;
    else if(n>0) param.boundary.boundary=rectangular;
    else param.boundary.boundary=truncoct;
    
    if(abs(n)==2){
      param.pcouple.calculate=true;
      if(param.pcouple.virial==no_virial)
	param.pcouple.virial=molecular_virial;
    }
  }

  if(!nrdbox && param.boundary.boundary != vacuum){
    messages.add("Illegal value for NRDBOX in BOUNDARY block (should be 1)",
		     "In_Parameter", Message::warning);
  }
}

/**
 * read the PERTURB block.
 */
void In_Parameter::read_PERTURB(Parameter &param,
				    std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  // try the new PERTURB03 block
  buffer = m_block["PERTURB03"];
  if (buffer.size()){

    block_read.insert("PERTURB03");


    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    std::string b, s1, s2;
    _lineStream >> s1 
		>> param.perturbation.lambda
		>> param.perturbation.dlamt
		>> param.perturbation.lambda_exponent
		>> s2;
    
    if (_lineStream.fail())
      messages.add("bad line in PERTURB block",
		       "In_Parameter", Message::error);
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    
    if (s1 == "on"){
      param.perturbation.perturbation = true;
      param.perturbation.scaled_only = false;
    }
    else if (s1 == "off") param.perturbation.perturbation = false;
    else if (s1 == "scaled") {
      param.perturbation.perturbation = true;
      param.perturbation.scaled_only = true;
    }
    else {
      std::istringstream css;
      css.str(s1);
      
      // param.perturbation.perturbation = (atoi(s1.c_str())==1);
      param.perturbation.scaled_only = false;
      css >> param.perturbation.perturbation;
      if(css.fail()){
	messages.add("bad value for NTG in PERTURB block:"+s1+"\n"
			 "on, scaled, off, 0, 1",
			 "In_Parameter", Message::error);
	param.perturbation.perturbation=false;
	return;
      }
    }
    
    if (s2 == "on") param.perturbation.scaling = true;
    else if (s2 == "off") param.perturbation.scaling = false;
    else {
      // param.perturbation.scaling = (atoi(s2.c_str())==1);
      std::istringstream css(s2);
      css >> param.perturbation.scaling;
      
      if(css.fail()){
	std::cerr << "got " << s2 << std::endl;
	messages.add("bad value for SCALING in PERTURB block\n"
			 "on,off,0,1",
			 "In_Parameter", Message::error);
	param.perturbation.scaling=false;
	return;
      }
    }
    if(param.perturbation.scaled_only && !param.perturbation.scaling){
      messages.add("inconsistent input: perturbing only scaled interactions, but scaling not turned on",
		       "In_Parameter", Message::error);
    }
  }
  else{
    // a common PERTURB block...
    buffer = m_block["PERTURB"];
    if (!buffer.size()){
      return;
    }

    block_read.insert("PERTURB");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    int ntg, nrdgl;
    double dmu, dmut;
    double mmu;
    double alpha_lj, alpha_crf;
    
    _lineStream >> ntg >> nrdgl 
		>> param.perturbation.lambda 
		>> param.perturbation.dlamt 
		>> dmu >> dmut
		>> alpha_lj >> alpha_crf 
		>> param.perturbation.lambda_exponent
		>> mmu;
    
    param.perturbation.scaling = false;
    param.perturbation.scaled_only = false;
    
    if (_lineStream.fail())
      messages.add("bad line in PERTURB block",
		       "In_Parameter", Message::error);
    
    if (nrdgl)
      messages.add("PERTURB: nrdgl != 0 not allowed",
		       "In_Parameter", Message::error);
    
    if (alpha_lj || alpha_crf){
      messages.add("PERTURB: softness constants taken from topology!",
		       "In_Parameter", Message::notice);
    }
    
    param.perturbation.perturbation=(ntg!=0);
    
    if (ntg != 0 && ntg != 1)
      messages.add("PERTURB: only ntg = 0 or ntg = 1 allowed",
		       "In_Parameter", Message::error);
  }
  
  if (param.perturbation.lambda_exponent<=0){
    messages.add("PERTURB: nlam > 0",
		     "In_Parameter", Message::error);
  }
}

/**
 * read FORCE block.
 */
void In_Parameter::read_FORCE(Parameter &param,
				  std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["FORCE"];

  if (!buffer.size()){
    messages.add("no FORCE block", "In_Parameter", Message::error);
    return;
  }

  block_read.insert("FORCE");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  int bondH, angleH, impH, dihedralH, charge;
  unsigned int num, e, old_e=0;
  
  _lineStream >> bondH >> param.force.bond >> angleH >> param.force.angle
	      >> impH >> param.force.improper >> dihedralH >> param.force.dihedral
	      >> charge >> param.force.nonbonded;
  _lineStream >> num;
  if(num<=0){
    messages.add("number of energy groups in FORCE block should be > 0",
		     "In_Parameter", Message::error);
    return;
  }
  
  for(unsigned int i=0; i<num; ++i){
    _lineStream >> e;
    param.force.energy_group.push_back(e-1);
    if(e<=old_e){
      messages.add("energy groups are not in order in FORCE block",
		       "In_Parameter", Message::error);
      return;
    }
    old_e = e;
  }
  
  if (_lineStream.fail())
    messages.add("bad line in ENERGYGROUP (FORCE) block",
		     "In_Parameter", Message::error);
  
  
  if (bondH ^ param.force.bond)
    messages.add("Force switch for bond and bond H has to be equal",
		     "In_Parameter", Message::error);

  if (angleH ^ param.force.angle)
    messages.add("Force switch for angle and angle H has to be equal",
		     "In_Parameter", Message::error);

  if (impH ^ param.force.improper)
    messages.add("Force switch for improper and improper H has to be equal",
		     "In_Parameter", Message::error);

  if (dihedralH ^ param.force.dihedral)
    messages.add("Force switch for dihedral and dihedral H has to be equal",
		     "In_Parameter", Message::error);

  if (charge ^ param.force.nonbonded)
    messages.add("Force switch for lj and charge has to be equal",
		     "In_Parameter", Message::error);

  if (_lineStream.fail())
    messages.add("bad line in FORCE block",
		       "In_Parameter", Message::error);

  //from here read the force field block
  read_FORCEFIELD(param);
  
  if(param.force.bond < 0 || param.force.bond > 2)
    messages.add("Illegal value for force switch for bond",
		     "In_Parameter", Message::error);
  if(param.force.bond < 0 || param.force.angle > 2)
    messages.add("Illegal value for force switch for angle",
		     "In_Parameter", Message::error);
  if(param.force.bond < 0 || param.force.angle == 2)
    messages.add("Force switch for angle = 2 currently not implemented",
		     "In_Parameter", Message::error);
  if(param.force.bond < 0 || param.force.improper > 1)
    messages.add("Illegal value for force switch for improper dihedral",
		     "In_Parameter", Message::error);
  if(param.force.bond < 0 || param.force.dihedral > 1)
    messages.add("Illegal value for force switch for dihedral",
		     "In_Parameter", Message::error);
  if(param.force.bond < 0 || param.force.nonbonded > 1)
    messages.add("Illegal value for force switch for nonbonded",
		     "In_Parameter", Message::error);
}

/**
 * read FORCEFIELD block.
 */
void In_Parameter::read_FORCEFIELD(Parameter &param,
				       std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["FORCEFIELD"];

  if (!buffer.size()){
    return;
  }
  
  block_read.insert("FORCEFIELD");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

  int bond, angle;
  
  _lineStream >> bond 
	      >> angle;

  if (bond == 2 && param.force.bond != 0){
    if (param.force.bond != 2){
      messages.add("using FORCEFIELD block to determine bond term",
		       "In_Parameter", Message::notice);
      param.force.bond = 2;
    }
  }
  if (bond == 0 && param.force.bond != 0){
    if (param.force.bond != 1){
      messages.add("using FORCEFIELD block to determine bond term",
		       "In_Parameter", Message::notice);
      param.force.bond = 1;
    }
  }
  
  if (angle != 0){
    messages.add("FORCEFIELD: only Gromos96 functional form for angle "
		     "bending allowed.",
		     "In_Parameter", Message::error);
  }
  
  if (_lineStream.fail())
    messages.add("bad line in FORCEFIELD block",
		     "In_Parameter", Message::error);
  
}

/**
 * read START block.
 */
void In_Parameter::read_START(Parameter &param,
				  std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["START"];
  
  if (!buffer.size()){
    messages.add("no START block", "In_Parameter", Message::error);
    return;
  }

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  block_read.insert("START");

  int ntx, init, ntx0;
  double heat;
  
  _lineStream >> ntx >> init 
	      >> param.start.ig 
	      >> param.start.tempi >> heat >> ntx0 
	      >> k_Boltzmann;
  
  if (_lineStream.fail())
    messages.add("bad line in START block",
		     "In_Parameter", Message::error);
  
  read_CENTREOFMASS(param);
  
  if(ntx==1 || param.start.tempi)
    param.start.generate_velocities = true;
  else
    param.start.generate_velocities = false;
  
  switch(init){
    case 1:
      param.start.shake_pos=true;
      param.start.shake_vel=true;
      break;
    case 2:
      param.start.shake_pos=false;
      param.start.shake_vel=true;
      break;
    case 3:
      param.start.shake_pos=false;
      param.start.shake_vel=false;
      break;
    case 4:
      param.start.shake_pos=false;
      param.start.shake_vel=false;
      param.start.remove_com=false;
      break;
    default:
      messages.add("Illegal option for init in START block",
		       "In_Parameter", Message::error); 
  }
  if(param.start.tempi <0)
    messages.add("Illegal value for TEMPI in START block (>=0)",
		     "In_Parameter", Message::error);
  if(heat)
    messages.add("HEAT != 0 is not supported in START block",
		     "In_Parameter", Message::error);
  if(ntx0!=1)
    messages.add("NTX0 != 1 is not supported in START block",
		     "In_Parameter", Message::error);
  if(k_Boltzmann <=0)
    messages.add("BOLTZ <=0 is not appreciated in START block",
		     "In_Parameter", Message::error);
}
/**
 * read CENTREOFMASS block.
 */
void In_Parameter::read_CENTREOFMASS(Parameter &param,
					 std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["CENTREOFMASS"];
  
  if (!buffer.size()){
    return;
  }

  block_read.insert("CENTREOFMASS");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  int ntcm;
  
  _lineStream >> param.centreofmass.ndfmin 
	      >> ntcm 
	      >> param.centreofmass.skip_step;

  if (param.centreofmass.skip_step){
    param.centreofmass.remove_rot = true;
    param.centreofmass.remove_trans = true;
  }
  else{
    param.centreofmass.remove_rot = false;
    param.centreofmass.remove_trans = false;
  }

  if (_lineStream.fail())
    messages.add("bad line in CENTREOFMASS block",
		     "In_Parameter", Message::error);
  if(ntcm!=0) 
    param.start.remove_com=true;
  
  if(param.centreofmass.ndfmin < 0)
    messages.add("Illegal value for NDFMIN in CENTREOFMASS block (>=0)",
		     "In_Parameter", Message::error);
  if(param.centreofmass.skip_step <0)
    messages.add("Illegal value for NSCM in CENTREOFMASS block (>=0)",
		     "In_Parameter", Message::error);
}

/**
 * read LONGRANGE block.
 */
void In_Parameter::read_LONGRANGE(Parameter &param,
				      std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["LONGRANGE"];
  
  if (!buffer.size()){
    messages.add("no LONGRANGE block in input",
		     "In_Parameter",Message::error);
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
    messages.add("bad line in LONGRANGE block",
		     "In_Parameter", Message::error);
  
  if(param.longrange.rf_cutoff < 0) {
    param.longrange.rf_excluded=true;
    param.longrange.rf_cutoff = -param.longrange.rf_cutoff;
    
  }
  else{
    param.longrange.rf_excluded=false;
  }
  if(param.longrange.rf_epsilon!=0 && param.longrange.rf_epsilon<1)
    messages.add("Illegal value for EPSRF in LONGRANGE block (0  / >=1)", 
		     "In_Parameter", Message::error);
  if(param.longrange.rf_kappa <0)
    messages.add("Illegal value for APPAK (who came up with this name?)"
		     " in LONGRANGE block (>=0)",
		     "In_Parameter", Message::error);
} // LONGRANGE

/**
 * read SUBMOLECULES block.
 */
void In_Parameter::read_SUBMOLECULES(Parameter &param,
					 std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["SUBMOLECULES"];
  
  if (!buffer.size()){
    messages.add("no SUBMOLECULES block in input",
		     "In_Parameter", Message::error);
    return;
  }

  block_read.insert("SUBMOLECULES");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  int num;
  
  _lineStream >> num;

  if(num<0){
    messages.add("negative number of SUBMOLECULES is not allowed",
		     "In_Parameter", Message::error);
    return;
  }
  
  unsigned int m;
  unsigned int old_m=0;
  
  param.submolecules.submolecules.push_back(0);

  for(int i=0; i<num; ++i){
    _lineStream >> m;
    param.submolecules.submolecules.push_back(m);
    if(m<old_m){
      messages.add("wrong order in SUBMOLECULES block",
		       "In_Parameter", Message::error);
      return;
    }
    old_m=m;
  }
  
  if (_lineStream.fail())
    messages.add("bad line in SUBMOLECULES block",
		     "In_Parameter", Message::error);
  
  
} // SUBMOLECULES


/**
 * read PLIST block.
 */
void In_Parameter::read_PLIST(Parameter &param,
				  std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;

  // try a PLIST03
  buffer = m_block["PLIST03"];
  if (buffer.size()){
    block_read.insert("PLIST03");

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
      messages.add("bad line in PLIST03 block",
		       "In_Parameter",
		       Message::error);
    }
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    std::transform(s3.begin(), s3.end(), s3.begin(), tolower);
    
    if (s1 == "grid") param.pairlist.grid = true;
    else if (s1 == "standard") param.pairlist.grid = false;
    else{
      messages.add("wrong pairlist algorithm chosen (allowed: standard, grid) in PLIST03 block",
		       "In_Parameter", Message::error);
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
	  messages.add("wrong pairlist grid size chosen (allowed: auto, [size]) in PLIST03 block",
			   "In_Parameter", Message::error);
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
	messages.add("wrong cutoff type chosen (allowed: atomic, chargegroup)"
			 " in PLIST03 block",
			 "In_Parameter", Message::error);
	param.pairlist.atomic_cutoff = false;
      }
    }
  }
  else{
    buffer = m_block["PLIST"];
    
    if (!buffer.size()){
      messages.add("no PLIST block in input","In_Parameter",Message::error);
      return;
    }

    block_read.insert("PLIST");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    int i;
    
    _lineStream >> i 
		>> param.pairlist.skip_step  
		>> param.pairlist.cutoff_short 
		>> param.pairlist.cutoff_long;

    param.pairlist.grid_size = 0;
    param.pairlist.grid = false;
    param.pairlist.atomic_cutoff = false;

    if (_lineStream.fail()){
      messages.add("bad line in PLIST block",
		       "In_Parameter", Message::error);
    }
    
  }
  if(param.pairlist.grid && param.pairlist.grid_size <=0)
    messages.add("Illegal value for grid size in PLIST03 block (>0)",
		     "In_Parameter", Message::error);
  if(param.pairlist.cutoff_short < 0){
    messages.add("Illegal value for short range cutoff in PLIST block (>0)",
		     "In_Parameter", Message::error);
  }
  if(param.pairlist.cutoff_long < param.pairlist.cutoff_short){
    messages.add("Illegal value for long range cutoff in PLIST block (>=RCUTP)",
		     "In_Parameter", Message::error);
  }
  
}


/**
 * read MULTIBATH block.
 */
void In_Parameter::read_MULTIBATH(Parameter &param,
				      std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  param.multibath.couple = false;

  // is there a MULTIBATH block
  buffer = m_block["MULTIBATH"];
  
  if (buffer.size()){

    block_read.insert("MULTIBATH");

    param.multibath.found_multibath=true;
    param.multibath.found_tcouple=false;
    
    messages.add("using MULTIBATH block",
		     "In_Parameter", Message::notice);

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
      std::stringstream ss(alg);
      if (!(ss >> param.multibath.nosehoover) ||
	  param.multibath.nosehoover < 0 || param.multibath.nosehoover > 2){
	messages.add("algorithm not understood in multibath block",
			 "In_Parameter", Message::error);

	param.multibath.nosehoover = 0;
	return;
      }
    }

    if (param.multibath.nosehoover == 2){
      // read in the number of chains
      int num;
      _lineStream >> num;

      if (num < 2){
	messages.add("wrong number of Nose-Hoover chains in multibath block",
			 "In_Parameter", Message::error);
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
	messages.add("illegal value for temp or tau in MULTIBATH block",
			 "In_Parameter", Message::error);
      }

      if (param.replica.T){
	os << "\tsetting temperature to " << param.replica.T << "K "
		  << "for replica exchange\n";
	temp = param.replica.T;
      }
      
      param.multibath.multibath.add_bath(temp, tau);
      if (tau != -1) param.multibath.couple = true;
    }
    
    if (_lineStream.fail()){
      messages.add("bad line in MULTIBATH block",
		       "In_Parameter", Message::error);
    }
    
    // now the ranges
    _lineStream >> num;
    
    for(int i=0; i<num; ++i){
      _lineStream >> last >> com_bath >> ir_bath;
      param.multibath.multibath.add_bath_index(last - 1, 0, com_bath - 1, ir_bath - 1);
    }
    
    if (_lineStream.fail()){
      messages.add("bad line in MULTIBATH block",
		       "In_Parameter", Message::error);
    }
    
  }
  else{
    // try a TCOUPLE block
    
    buffer = m_block["TCOUPLE"];
    if (buffer.size()){
      block_read.insert("TCOUPLE");

      param.multibath.found_multibath=false;
      param.multibath.found_tcouple=true;
      
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

      
      for(int i=0; i<3; ++i){
	_lineStream >> param.multibath.tcouple.ntt[i] 
		    >> param.multibath.tcouple.temp0[i]
		    >> param.multibath.tcouple.tau[i];
      }
      
      if (_lineStream.fail()){
	messages.add("bad line in TCOUPLE block",
			 "In_Parameter", Message::error);
	return;
      }	

      if (param.replica.T){
	os << "\tsetting temperature to " << param.replica.T << "K "
		  << "for replica exchange\n";
	param.multibath.tcouple.temp0[0] = param.replica.T;
	param.multibath.tcouple.temp0[1] = param.replica.T;
	param.multibath.tcouple.temp0[2] = param.replica.T;
      }
    }
    else{
      param.multibath.found_multibath=false;
      param.multibath.found_tcouple=false;
      // no TCOUPLE block
      // that's fine, same as 0,0,0
    }
    
  }
  
} // TEMPERATURE coupling

/**
 * read POSREST block.
 */
void In_Parameter::read_POSREST(Parameter &param,
				    std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["POSREST"];
  
  if (!buffer.size()){
    return;
  }
  
  block_read.insert("POSREST");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  _lineStream >> param.posrest.posrest
	      >> param.posrest.force_constant
	      >> param.posrest.nrdrx;
  
  if (_lineStream.fail())
    messages.add("bad line in POSREST block",
		     "In_Parameter", Message::error);
  
  if(param.posrest.posrest == 3) {
    messages.add("Position constraining is experimental",
		     "In_Parameter", 
		     Message::warning);

    if (param.pcouple.scale != pcouple_off){
      messages.add("Position constraining together with pressure coupling not allowed",
		       "In_Parameter",
		       Message::error);
    }
  }

  if(param.posrest.force_constant <0)
    messages.add("Illegal value for CHO"
		     " in POSREST block (>=0)",
		     "In_Parameter", Message::error);

} // POSREST

/**
 * read DISTREST block.
 */
void In_Parameter::read_DISTREST(Parameter &param,
				    std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["DISTREST"];
  
  if (!buffer.size()){
    return;
  }
  
  block_read.insert("DISTREST");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
  
  _lineStream >> param.distrest.distrest
	      >> param.distrest.K
	      >> param.distrest.r_linear
	      >> param.distrest.tau
	      >> param.distrest.read;
  
  if (_lineStream.fail())
    messages.add("bad line in DISTREST block",
		     "In_Parameter", Message::error);
  
  if(param.distrest.distrest <0) {
    messages.add("Distance restrain averaging not implemented",
		     "In_Parameter", 
		     Message::warning);
  }
  
  if(param.distrest.distrest >2) {
    messages.add("bad input in DISTREST block, NTDR must be <=2",
		     "In_Parameter", 
		     Message::warning);
    
  }
  


  if(param.distrest.K <0)
    messages.add("Illegal value for force constant"
		     " in DISTREST block (>=0)",
		     "In_Parameter", Message::error);

} // DISTREST

/**
 * read the JVALUE block.
 */
void In_Parameter::read_JVALUE(Parameter &param,
				   std::ostream & os)
{

  std::vector<std::string> buffer;
  std::string s;
  
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
      messages.add("bad line in J-VAL block",
		       "In_Parameter", Message::error);
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    
    if (s1 == "off") param.jvalue.mode = restr_off;
    else if (s1 == "instantaneous") param.jvalue.mode = restr_inst;
    else if (s1 == "averaged") param.jvalue.mode = restr_av;
    else if (s1 == "biquadratic") param.jvalue.mode = restr_biq;
    else{
      std::istringstream css;
      css.str(s1);
      
      int i;
      css >> i;
      if(css.fail() || i < 0 || i > 3){
	messages.add("bad value for MODE in J-VAL block:"+s1+"\n"
			 "off, instantaneous, averaged, biquadratic, scaled (0-4)",
			 "In_Parameter", Message::error);
	param.jvalue.mode = restr_off;
	return;
      }
      param.jvalue.mode = restr_enum(i);

    }
    if (param.jvalue.tau < 0 ||
	(param.jvalue.tau == 0 && (param.jvalue.mode != restr_off ||
				   param.jvalue.mode != restr_inst))){
      messages.add("bad value for TAU in J-VAL block\n"
		       "should be > 0.0",
		       "In_Parameter", Message::error);
    }
  }
} // JVALUE

/**
 * read the PSCALE block.
 */
void In_Parameter::read_PSCALE(Parameter &param,
				   std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["PSCALE"];
  if (buffer.size()){

    block_read.insert("PSCALE");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    std::string s1;
    _lineStream >> s1;
    
    if (_lineStream.fail()){
      messages.add("bad line in PSCALE block",
		       "In_Parameter", Message::error);
      return;
    }
    
    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    
    if (s1 == "jrest"){

      param.pscale.jrest = true;

      _lineStream >> param.pscale.KDIH >> param.pscale.KJ 
		  >> param.pscale.T >> param.pscale.difference
		  >> param.pscale.ratio >> param.pscale.read_data;

      if (_lineStream.fail())
	messages.add("bad line in PSCALE block",
			 "In_Parameter", Message::error);
      if (param.pscale.KDIH < 0.0)
	messages.add("bad value for KDIH in PSCALE block (negative)",
			 "In_Parameter", Message::error);
      if (param.pscale.KJ < 0.0)
	messages.add("bad value for KJ in PSCALE block (negative)",
			 "In_Parameter", Message::error);
      if (param.pscale.T < 0.0)
	messages.add("bad value for T in PSCALE block (negative)",
			 "In_Parameter", Message::error);
      if (param.pscale.difference < 0.0)
	messages.add("bad value for difference in PSCALE block (negative)",
			 "In_Parameter", Message::error);
    }
    else{
      messages.add("bad value for periodic scaling mode", "In_Parameter", Message::error);
      return;
    }
  }
} // PSCALE


/**
 * read the ROTTRANS block.
 */
void In_Parameter::read_ROTTRANS(Parameter &param,
				     std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["ROTTRANS"];

  if (buffer.size()){
    int i;

    block_read.insert("ROTTRANS");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> i;
    
    if (_lineStream.fail())
      messages.add("bad line in ROTTRANS block",
		       "In_Parameter", Message::error);

    param.rottrans.rottrans = (i != 0);

  }
  
}

/**
 * read the INNERLOOP block.
 */
void In_Parameter::read_INNERLOOP(Parameter &param,
				      std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["INNERLOOP"];

  if (buffer.size()){
    int spc;

    block_read.insert("INNERLOOP");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> spc;
    
    if (_lineStream.fail())
      messages.add("bad line in INNERLOOP block",
		       "In_Parameter", Message::error);

    if (spc != -1 && spc != 0){
      messages.add("bad value for SPCL in INNERLOOP: allowed : -1, 0",
		       "In_Parameter",
		       Message::error);
      spc = -1;
    }
    param.force.spc_loop = spc;

  }
  
}

/**
 * read the REPLICA block.
 */
void In_Parameter::read_REPLICA(Parameter &param,
				    std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["REPLICA"];

  if (buffer.size()){

    block_read.insert("REPLICA");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.replica.ID >> param.replica.T >> param.replica.scale;
    
    if (_lineStream.fail()){
      messages.add("bad line in REPLICA block",
		       "In_Parameter", Message::error);

      param.replica.ID = 0;
      param.replica.T = 0.0;
      param.replica.scale = 0;
    }

  }
  
}

/**
 * read the REPLICA03 block.
 */
void In_Parameter::read_REPLICA03(Parameter &param,
				      std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["REPLICA03"];

  if (buffer.size()){

    block_read.insert("REPLICA03");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.replica.number;
    param.replica.temperature.resize(param.replica.number, 0.0);
    param.replica.lambda.resize(param.replica.number, 0.0);

    for(int i=0; i<param.replica.number; ++i){
      _lineStream >> param.replica.temperature[i];
    }
    for(int i=0; i<param.replica.number; ++i){
      _lineStream >> param.replica.lambda[i];
    }

    _lineStream >> param.replica.trials;
    
    if (_lineStream.fail()){
      messages.add("bad line in REPLICA03 block",
		       "In_Parameter", Message::error);

      param.replica.number = 0;
      param.replica.temperature.clear();
      param.replica.lambda.clear();
    }
  }
  
}

void In_Parameter::read_MULTICELL(Parameter & param,
				      std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["MULTICELL"];

  if (buffer.size()){

    block_read.insert("MULTICELL");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.multicell.multicell 
		>> param.multicell.x >> param.multicell.y >> param.multicell.z;
    
    if (_lineStream.fail()){
      messages.add("bad line in MULTICELL block",
		       "In_Parameter", Message::error);

      param.multicell.multicell = false;
    }
  }
  
}

void In_Parameter::read_ANALYZE(Parameter & param,
				    std::ostream & os)
{
  std::vector<std::string> buffer;
  std::string s;
  
  buffer = m_block["ANALYZE"];

  if (buffer.size()){

    block_read.insert("ANALYZE");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));
    
    _lineStream >> param.analyze.analyze;
    
    if (_lineStream.fail()){
      messages.add("bad line in ANALYZE block",
		       "In_Parameter", Message::error);
      
      param.analyze.analyze = false;
    }
  }
  
}

