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

// mk_script.h
#include "../src/gcore/Box.h"
#include "../src/gio/Ginstream.h"

void printIO(std::string b, std::string var, std::string val, std::string allow);
void printErrMsg(std::string block, std::string variable, std::string message);

void printWarning(std::string s);
void printError(std::string s);

// reads a value/variable from a stringstraem and saves it in a variable
// returns true if conversion (input variable -> target variable) was possible
// returns false if conversion failed or end up with loss of data (double -> integer)
bool readValue(std::string BLOCK, std::string VAR, std::istringstream &is, int &variable, std::string allowed, bool allow_missing=false);
bool readValue(std::string BLOCK, std::string VAR, std::istringstream &is, double &variable, std::string allowed, bool allow_missing=false);

// We define two global variables. It makes the printing and bookkeeping 
// of the errors and warnings a lot cleaner.
int numWarnings = 0;
int numErrors = 0;
int numTotWarnings = 0;
int numTotErrors = 0;

enum filetype {
  unknownfile, inputfile, topofile, coordfile, refposfile, anatrxfile,
  posresspecfile, xrayfile, disresfile, pttopofile, gamdfile, dihresfile, angresfile, jvaluefile, orderfile,
  symfile, colvarresfile,
  ledihfile, leumbfile, bsleusfile, qmmmfile, frictionfile, outputfile, outtrxfile, outtrvfile,
  outtrffile, outtrefile, outtrgfile,
  jinfile, joutfile, jtrjfile,
  scriptfile, outbaefile, outbagfile,
  outtrsfile, repoutfile, repdatfile
};

typedef std::map<std::string, filetype>::value_type FT;
const FT filetypes[] = {FT("", unknownfile),
  FT("input", inputfile),
  FT("topo", topofile),
  FT("coord", coordfile),
  FT("refpos", refposfile),
  FT("anatrj", anatrxfile),
  FT("posresspec", posresspecfile),
  FT("xray", xrayfile),
  FT("disres", disresfile),
  FT("colvarres", colvarresfile),
  FT("pttopo", pttopofile),
  FT("gamd",   gamdfile),
  FT("dihres", dihresfile),
  FT("angres", angresfile),
  FT("jvalue", jvaluefile),
  FT("order", orderfile),
  FT("sym", symfile),
  FT("ledih", ledihfile),
  FT("leumb", leumbfile),
  FT("bsleus", bsleusfile),
  FT("qmmm", qmmmfile),
  FT("jin", jinfile),
  FT("jout", joutfile),
  FT("jtrj", jtrjfile),
  FT("friction", frictionfile),
  FT("output", outputfile),
  FT("outtrx", outtrxfile),
  FT("outtrv", outtrvfile),
  FT("outtrf", outtrffile),
  FT("outtre", outtrefile),
  FT("outtrg", outtrgfile),
  FT("outbae", outbaefile),
  FT("outbag", outbagfile),
  FT("outtrs", outtrsfile),
  FT("script", scriptfile),
  FT("repout", repoutfile),
  FT("repdat", repdatfile)};
const int numFiletypes = sizeof(filetypes)/sizeof(FT);

static std::map<std::string, filetype> FILETYPE(filetypes, filetypes + numFiletypes);

enum blocktype {
  unknown, addecoupleblock, aedsblock, barostatblock, boundcondblock, bsleusblock,
  cgrainblock, comtransrotblock, consistencycheckblock,
  constraintblock, covalentformblock, debugblock,
  dihedralresblock, angleresblock, distancefieldblock, distanceresblock, edsblock, 
  energyminblock,
  ewarnblock, forceblock, gamdblock, geomconstraintsblock,
  gromos96compatblock, initialiseblock, innerloopblock,
  integrateblock, jvalueresblock, lambdasblock,
  localelevblock, electricblock, multibathblock, multicellblock, 
  multigradientblock, multistepblock,
  neighbourlistblock, nemdblock, nonbondedblock, 
  orderparamresblock, overalltransrotblock,
  pairlistblock, pathintblock, perscaleblock,
  perturbationblock, polariseblock, positionresblock,
  pressurescaleblock, precalclamblock, printoutblock, qmmmblock,
  randomnumbersblock, readtrajblock, replicablock, reedsblock, rottransblock,
  sasablock, stepblock, stochdynblock, symresblock, systemblock,
  thermostatblock, umbrellablock, virialblock, virtualatomblock,
  writetrajblock, xrayresblock, colvarresblock
};

typedef std::map<std::string, blocktype>::value_type BT;
const BT blocktypes[] = {BT("", unknown),
  BT("ADDECOUPLE", addecoupleblock),
  BT("AEDS", aedsblock),
  BT("BAROSTAT", barostatblock),
  BT("BOUNDCOND", boundcondblock),
  BT("BSLEUS", bsleusblock),
  BT("CGRAIN", cgrainblock),
  BT("COMTRANSROT", comtransrotblock),
  BT("CONSISTENCYCHECK", consistencycheckblock),
  BT("CONSTRAINT", constraintblock),
  BT("COVALENTFORM", covalentformblock),
  BT("DEBUG", debugblock),
  BT("DIHEDRALRES", dihedralresblock),
  BT("ANGLERES", angleresblock),
  BT("DISTANCERES", distanceresblock),
  BT("DISTANCEFIELD", distancefieldblock),
  BT("EDS", edsblock),
  BT("ELECTRIC", electricblock),
  BT("ENERGYMIN", energyminblock),
  BT("EWARN", ewarnblock),
  BT("FORCE", forceblock),
  BT("GAMD", gamdblock),
  BT("GEOMCONSTRAINTS", geomconstraintsblock),
  BT("GROMOS96COMPAT", gromos96compatblock),
  BT("INITIALISE", initialiseblock),
  BT("INNERLOOP", innerloopblock),
  BT("INTEGRATE", integrateblock),
  BT("JVALUERES", jvalueresblock),
  BT("LAMBDAS", lambdasblock),
  BT("LOCALELEV", localelevblock),
  BT("MULTIBATH", multibathblock),
  BT("MULTICELL", multicellblock),
  BT("MULTIGRADIENT", multigradientblock),
  BT("MULTISTEP", multistepblock),
  BT("NEIGHBOURLIST", neighbourlistblock),
  BT("NEMD", nemdblock),
  BT("NONBONDED", nonbondedblock),
  BT("ORDERPARAMRES", orderparamresblock),
  BT("OVERALLTRANSROT", overalltransrotblock),
  BT("PAIRLIST", pairlistblock),
  BT("PATHINT", pathintblock),
  BT("PERSCALE", perscaleblock),
  BT("PERTURBATION", perturbationblock),
  BT("POLARISE", polariseblock),
  BT("POSITIONRES", positionresblock),
  BT("PRECALCLAM", precalclamblock),
  BT("PRESSURESCALE", pressurescaleblock),
  BT("PRINTOUT", printoutblock),
  BT("QMMM", qmmmblock),
  BT("RANDOMNUMBERS", randomnumbersblock),
  BT("READTRAJ", readtrajblock),
  BT("REPLICA", replicablock),
  BT("REPLICA_EDS", reedsblock),
  BT("ROTTRANS", rottransblock),
  BT("SASA", sasablock),
  BT("STEP", stepblock),
  BT("STOCHDYN", stochdynblock),
  BT("SYMRES", symresblock),
  BT("SYSTEM", systemblock),
  BT("THERMOSTAT", thermostatblock),
  BT("UMBRELLA", umbrellablock),
  BT("VIRIAL", virialblock),
  BT("VIRTUALATOM", virtualatomblock),
  BT("WRITETRAJ", writetrajblock),
  BT("XRAYRES", xrayresblock),
  BT("COLVARRES", colvarresblock)};
const int numBlocktypes = sizeof(blocktypes)/sizeof(BT);
static std::map<std::string, blocktype> BLOCKTYPE(blocktypes, blocktypes + numBlocktypes);

enum templateelement {
  unknowntemplate, systemtemplate, numbertemplate,
  oldnumbertemplate,
  start_timetemplate, end_timetemplate, queuetemplate
};
typedef std::map<std::string, templateelement>::value_type TE;
const TE templateelements[] = {TE("", unknowntemplate),
  TE("system", systemtemplate),
  TE("number", numbertemplate),
  TE("oldnumber", oldnumbertemplate),
  TE("start_time", start_timetemplate),
  TE("end_time", end_timetemplate),
  TE("queue", queuetemplate)};
static std::map<std::string, templateelement> TEMPLATE(templateelements,
        templateelements + 7);

//BLOCKDEFINITIONS
class iaddecouple {
public:
  int found, adgr,  write;
  std::vector<int> adstart, adend;
  std::vector<double> sm, sv, st, tir;
  double tmf;

  iaddecouple() {
    adgr = 0; 
    found = 0;
    write = 0;
  }
};

class iaeds {
public:
  int found, aeds, form, numstates, ntiaedss, restremin, bmaxtype, asteps, bsteps;
  double emax, emin, bmax;
  std::vector<double> eir;
  
  iaeds() {
    found = 0;
  }
};

class ibarostat {
public:
  int found, ntp, npvar, npcpl[6];
  double comp;

  class pbath {
  public:
    double prsbth;
    std::vector<double> taubba;
  };
  std::vector<pbath> pbaths;

  ibarostat() {
    found = 0;
  }
};

class iboundcond {
public:
  int found, ntb, ndfmin;

  iboundcond() {
    found = 0;
  }
};

class ibsleus {
public:
  int found, bsleus, build, write;
  
  ibsleus() {
    found = 0;
    write = 0;
  }
};

class icgrain {
public:
  int found, ntcgran;
  double eps, epsm;

  icgrain() {
    found = 0;
  }
};

class icolvarres {
public:
  int found, cvr, vcvr, ntwcv;
  double cvk, taucvr;

  icolvarres() {
    found = 0;
    ntwcv = 0;
  }
};

class icomtransrot {
public:
  int found, nscm;

  icomtransrot() {
    found = 0;
  }
};

class iconsistencycheck {
public:
  int found, ntchk, ntckf, ntckv, ntckt;
  int ntcke, ntckr, ntckl;
  double fdckf, fdckv, fdckl;
  std::vector<int> nckf;

  iconsistencycheck() {
    found = 0;
  }
};

class iconstraint {
public:
  int found, ntc, ntcp, ntcs, ntcg;
  std::vector<int> ntcd;
  double ntcp0[3], ntcs0[3];

  iconstraint() {
    found = 0;
    for (int i = 0; i < 0; ++i) {
      ntcp0[i] = ntcs0[i] = -1.0;
    }
  }
};

class icovalentform {
public:
  int found, ntbbh, ntbah, ntbdn;

  icovalentform() {
    found = 0;
  }
};

class idebug {
public:
  int found;

  class routine {
  public:
    std::string piider;
    int iiideo;
  };
  std::vector<routine> routines;

  idebug() {
    found = 0;
  }
};

class idihedralres {
public:
  int found, ntdlr, ntwdlr, vdih;
  double cdlr, philin, toldac;

  idihedralres() {
    found = 0;
    ntwdlr = 0;
  }
};

class iangleres {
public:
  int found, ntalr, ntwalr, vares;
  double calr, tolbac;

  iangleres() {
    found = 0;
    ntwalr = 0;
  }
};

class idistanceres {
public:
  int found, ntdir, ntdira, ntwdir, vdir, forcescale;
  double cdir, dir0, taudir;

  idistanceres() {
    found = 0;
    ntwdir = 0;
  }
};

class idistancefield {
public:
  int found, ntdfr, update, smooth, ntwdf, printgrid;
  double grid, proteinoffset, proteincutoff, rl, protect;
  
  idistancefield() {
    found = 0;
    ntwdf = 0;
  }
};

class ieds {
public:
  int found, eds, form, numstates;
  std::vector<double> eir, smooth;
  std::vector<std::vector<int> > tree;
  
  ieds() {
    found = 0;
  }
};

class ienergymin {
public:
  int found, ntem, ncyc, nmin, cgim;
  double dele, dx0, dxm, flim, cgic;

  ienergymin() {
    found = 0;
    cgim = 1;
    cgic = 0.0;
  }
};

class iewarn {
public:
  int found;
  double maxener;

  iewarn() {
    found = 0;
  }
};

class iforce {
public:
  int found, ntf[6];
  bool force_groups;
  std::vector<int> nre;

  iforce() {
    found = 0;
    force_groups = false;
  }
};

class igamd {
public:
  int found, gamd, search, form, thresh, ntigamds, eqsteps, window, agroups, igroups;  
  double dihstd, totstd;
  std::vector<double> ed;
  std::vector<double> et;
  std::vector<double> kd;
  std::vector<double> kt;
  igamd() {
    found = 0;
  }
};

class igeomconstraints {
public:
  int found, ntcph, ntcpn, ntcs;
  double shktol;

  igeomconstraints() {
    found = 0;
  }
};

class igromos96compat {
public:
  int found, ntnb96, ntr96, ntp96, ntg96;

  igromos96compat() {
    found = 0;
  }
};

class iinitialise {
public:
  int found, ntivel, ntishk, ntinht, ntinhb;
  int ntishi, ntirtc, nticom, ntisti;
  double ig, tempi;

  iinitialise() {
    found = 0;
  }
};

class iinnerloop {
public:
  int found, ntilm, ntils, ngpus;
  std::vector<int> ndevg;

  iinnerloop() {
    found = 0;
    ngpus = 0;
  }
};

class iintegrate {
public:
  int found, nint;

  iintegrate() {
    found = 0;
  }
};

class ijvalueres {
public:
  int found, ntjvr, ntjvra, le, ngrid, write, njvrtars, njvrbiqw;
  double cjvr, taujvr, delta;

  ijvalueres() {
    found = 0;
    write = 0;
  }
};

class ilambdas {
public:
  int found, ntil;

  class lambint {
  public:
    int ntli, nilg1, nilg2;
    double ali, bli, cli, dli, eli;
  };
  std::vector<lambint> lambints;

  ilambdas() {
    found = 0;
  }
};

class ilocalelev {
public:
  int found, ntles, nlepot, ntlesa, ntwle;
  std::map<int, int> nlepid_ntlepfr;

  ilocalelev() {
    found = 0;
    ntwle = 0;
  }
};

class ielectric {
public:
  int found, field, dipole, current;
  double ef_x, ef_y, ef_z;
  int dipgrp, ntwdip, ntwcur, ncurgrp;
  std::vector<int> curgrp;

  ielectric() {
    found = 0;
    current = 0;
    dipole = 0;
  }
};

class imultibath {
public:
  int found, ntbtyp, num, nbaths, dofset;
  std::vector<double> temp0, tau;
  std::vector<int> last, combath, irbath;

  imultibath() {
    found = 0;
    num = -1;
  }
};

class imulticell {
public:
  int found, ntm, ncella, ncellb, ncellc;
  double tolpx, tolpv, tolpf, tolpfw;

  imulticell() {
    found = 0;
  }
};

class imultigradient {
public:
  int found, ntmgre, ntmgrp, ntmgrn;
  std::vector<std::string> mgrvar;
  std::vector<int> mgrfrm;
  std::vector<std::vector<std::pair<double, double> > > curves;

  imultigradient() {
    found = 0;
    ntmgre = 0;
  }
};

class imultistep {
public:
  int found, steps, boost;

  imultistep() {
    found = 0;
  }
};

class ineighbourlist {
public:
  int found, plalgo, nupdpl, nupdis, nupdii, type, ncgcen;
  double rcuts, rcuti, gridszx, gridszy, gridszz;

  ineighbourlist() {
    found = 0;
  }
};

class inemd {
public:
 int found, nemd, property, method, slabnum, pertfrq, stdyaft, write;
 double ampli;
 
 inemd(){
  nemd =0;
  found = 0;
  write = 0;
 }
};


class inonbonded {
public:
  int found, nlrele, nslfexcl, nshape, na2clc, nkx, nky, nkz, ngx, ngy, ngz;
  int nasord, nfdord, nalias, nqeval, nrdgrd, nwrgrd, nlrlj;
  double appak, rcrf, epsrf, ashape, tola2, epsls, kcut, nspord, faccur, slvdns;

  inonbonded() {
    found = 0;
  }
};

class iorderparamres {
public:
  int found, ntopr, ntopra, ntwop, updopr;
  double copr, tauopr;
  
  iorderparamres() {
    found = 0;
    ntwop = 0;
  }
};

class ioveralltransrot {
public:
  int found, ncmtr, ncmro;
  double cmamx, cmamy, cmamz;

  ioveralltransrot() {
    found = 0;
  }
};

class ipairlist {
public:
  int found, algorithm, nsnb, type;
  double rcutp, rcutl, size;

  ipairlist() {
    found = 0;
  }
};

class ipathint {
public:
  int found, ntpi;

  ipathint() {
    found = 0;
  }
};

class iperscale {
public:
  int found, read, restype;
  double t, kdih, kj, diff, ratio;

  iperscale() {
    found = 0;
  }
};

class iperturbation {
public:
  int found, ntg, nrdgl, nlam, nscale;
  double rlam, dlamt, alphlj, alphc;

  iperturbation() {
    found = 0;
  }
};

class ipolarise {
public:
  int found, cos, efield, damp, write;
  double minfield;

  ipolarise() {
    found = 0;
    write = 0;
  }
};

class ipositionres {
public:
  int found, ntpor, ntporb, ntpors;
  double cpor;

  ipositionres() {
    found = 0;
  }
};

class iprecalclam {
public:
  int found, nrlam;
  double minlam, maxlam;

  iprecalclam() {
    found = 0;
  }
};

class ipressurescale {
public:
  int found, couple, scale, virial;
  int x_semi, y_semi, z_semi;
  double comp, taup, pres0[3][3];

  ipressurescale() {
    found = 0;
  }
};

class iprintout {
public:
  int found, ntpr, ntpp;

  iprintout() {
    found = 0;
  }
};

struct iqmmm {
  int found, ntqmmm, ntqmsw, ntwqmmm, qmlj, qmcon;
  double rcutqm, mmscale;

  iqmmm() {
    found = 0;
    mmscale = -1.0;
  }
};

class irandomnumbers {
public:
  int found, ntrng, ntgsl;

  irandomnumbers() {
    found = 0;
  }
};

class ireadtraj {
public:
  int found, ntrd, ntstr, ntrb, ntshk;

  ireadtraj() {
    found = 0;
  }
};

class ireplica {
public:
  int found, nrel, lrescale, nretrial, nrequil, cont;
  std::vector<double> ret, relam, rets;

  ireplica() {
    found = 0;
    cont = 0;
  }
};

class ireeds {
public:
  int found, reeds, nres, numstates, neoff, nretrial, nrequil, cont, eds_stat_out, periodic;
  std::vector<double> res, eir;

  ireeds() {
    found = 0;
    cont = 0;
  }
};

class irottrans {
public:
  int found, rtc, rtclast;

  irottrans() {
    found = 0;
  }
};

class isasa {
public:
  int found, ntsasa, ntvol;
  double p12, p13, p1x, sigmav, rsolv, as1, as2;

  isasa() {
    found = 0;
  }
};

class istep {
public:
  int found, nstlim;
  double t, dt;

  istep() {
    found = 0;
  }
};

class istochdyn {
public:
  int found, ntsd, ntfr, nsfr, nbref;
  double rcutf, cfric, tempsd;

  istochdyn() {
    found = 0;
  }
};

class isymres {
public:
  int found, ntsym;
  double csym;
  
  isymres() {
    found = 0;
  }
};

class isystem {
public:
  int found, npm, nsm;

  isystem() {
    found = 0;
  }
};

class ithermostat {
public:
  int found, ntt, ntbth, ntset;

  class tbath {
  public:
    int index, ntbtyp, ntbvar;
    double tembth;
    std::vector<double> taubth;
  };
  std::vector<tbath> baths;

  class tdofgroup {
  public:
    int ntscpl, ntstyp, ntscns, ntsgt;
    std::vector<int> ntsgtg;
  };
  std::vector<tdofgroup> dofgroups;

  ithermostat() {
    found = 0;
  }
};

class iumbrella {
public:
  int found, ntus;
  double uscst1, uscst2, usref1, usref2;

  iumbrella() {
    found = 0;
  }
};

class ivirial {
public:
  int found, ntv, ntvg;

  ivirial() {
    found = 0;
  }
};

class ivirtualatom {
public:
  int found, virt, numvirt, lastvirt;

  ivirtualatom() {
    found = 0;
  }
};

class iwritetraj {
public:
  int found, ntwx, ntwse, ntwv, ntwf, ntwe, ntwg, ntwb;

  iwritetraj() {
    found = 0;
    ntwx = 0;
    ntwse = 0;
    ntwv = 0;
    ntwf = 0;
    ntwe = 0;
    ntwg = 0;
    ntwb = 0;
  }
};

class ixrayres {
public:
  int found, ntxr, ntxle, ntwxr, ntwde, ntwxm, rdavg;
  double cxr, cxtau;

  ixrayres() {
    found = 0;
    ntwxr = 0;
  }
};

class iunknown {
public:
  std::string name;
  std::string content;

  iunknown(std::string n) : name(n) {
  }
};

class input {
public:

  iaddecouple addecouple;
  iaeds aeds;
  ibarostat barostat;
  iboundcond boundcond;
  ibsleus bsleus;
  icgrain cgrain;
  icolvarres colvarres;
  icomtransrot comtransrot;
  iconsistencycheck consistencycheck;
  iconstraint constraint;
  icovalentform covalentform;
  idebug debug;
  idihedralres dihedralres;
  iangleres angleres;
  idistancefield distancefield;
  idistanceres distanceres;
  ieds eds;
  ielectric electric;
  ienergymin energymin;
  iewarn ewarn;
  iforce force;
  igamd  gamd;
  igeomconstraints geomconstraints;
  igromos96compat gromos96compat;
  iinitialise initialise;
  iinnerloop innerloop;
  iintegrate integrate;
  ijvalueres jvalueres;
  ilambdas lambdas;
  ilocalelev localelev;
  imultibath multibath;
  imulticell multicell;
  imultigradient multigradient;
  imultistep multistep;
  ineighbourlist neighbourlist;
  inemd nemd; 
  inonbonded nonbonded;
  iorderparamres orderparamres;
  ioveralltransrot overalltransrot;
  ipairlist pairlist;
  ipathint pathint;
  iperscale perscale;
  iperturbation perturbation;
  ipolarise polarise;
  ipositionres positionres;
  iprecalclam precalclam;
  ipressurescale pressurescale;
  iprintout printout;
  irandomnumbers randomnumbers;
  ireadtraj readtraj;
  ireplica replica;
  ireeds reeds;
  irottrans rottrans;
  isasa sasa;
  istep step;
  istochdyn stochdyn;
  isymres symres;
  isystem system;
  ithermostat thermostat;
  iumbrella umbrella;
  ivirial virial;
  ivirtualatom virtualatom;
  iwritetraj writetraj;
  ixrayres xrayres;
  iqmmm qmmm;
  std::vector<iunknown> unknown;
};

class fileInfo {
public:
  gcore::Box box;
  std::vector<std::string> blocks;
  std::vector<int> blockslength;
};

// INSTREAM

// Andreas:
// no checks (except the block size), just read in the data/parameters and check
// if the read in was ok
// more complicated checks are made later, together with the cross checks



std::istringstream & operator>>(std::istringstream &is, iaddecouple &s) {
  s.found = 1;
  readValue("ADDECOUPLE", "ADGR", is, s.adgr, ">=0");
  if(s.adgr>0){
    for (int i=0; i < s.adgr; ++i){
      double sm, sv, st, tir;
      int  adstart,adend;
      std::stringstream blockName;
      blockName << "ADSTART[" << i + 1 << "]";
      readValue("ADDECOUPLE",  blockName.str(), is, adstart, ">0");
      blockName << "ADEND[" << i + 1 << "]";
      readValue("ADDECOUPLE",  blockName.str(), is, adend, ">0");
      blockName << "SM[" << i + 1 << "]";
      readValue("ADDECOUPLE",  blockName.str(), is, sm, ">0");
      blockName << "SV[" << i + 1 << "]";
      readValue("ADDECOUPLE",  blockName.str(), is, sv, ">0");
      blockName << "ST[" << i + 1 << "]";
      readValue("ADDECOUPLE",  blockName.str(), is, st, ">0");
      blockName << "TIR[" << i + 1 << "]";
      readValue("ADDECOUPLE",  blockName.str(), is, tir, ">0");
      s.adstart.push_back(adstart);
      s.adend.push_back(adend);
      s.sm.push_back(sm);
      s.sv.push_back(sv);
      s.st.push_back(st);
      s.tir.push_back(tir);

    }
    readValue("ADDECOUPLE", "TMF", is, s.tmf, ">=0");
    readValue("ADDECOUPLE", "WRITE", is, s.write, "<=0");
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iaeds &s) {
  s.found = 1;
  readValue("AEDS", "AEDS", is, s.aeds, "0,1");
  readValue("AEDS", "FORM", is, s.form, "1..4");
  readValue("AEDS", "NUMSTATES", is, s.numstates, ">1");
  if (s.numstates <= 1) {
    std::stringstream ss;
    ss << s.numstates;
    printIO("AEDS", "NUMSTATES", ss.str(), ">1");
  }
  readValue("AEDS", "EMAX", is, s.emax, ">0.0");
  readValue("AEDS", "EMIN", is, s.emin, ">0.0");
  s.eir.resize(s.numstates);
  for (int N = 0; N < s.numstates; N++) {
    std::stringstream blockName;
    blockName << "EIR[" << N + 1 << "]";
    readValue("EDS", blockName.str(), is, s.eir[N], ">0.0");
  }
  readValue("AEDS", "NTIAEDSS", is, s.ntiaedss, "0,1");
  readValue("AEDS", "RESTREMIN", is, s.restremin, "0,1");
  readValue("AEDS", "BMAXTYPE", is, s.bmaxtype, "1,2");
  readValue("AEDS", "BMAX", is, s.bmax, ">0.0");
  readValue("AEDS", "ASTEPS", is, s.asteps, ">0");
  readValue("AEDS", "BSTEPS", is, s.bsteps, ">0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of EDS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ibarostat &s) {
  int dum, npbth;
  double taup;
  s.found = 1;
  readValue("BAROSTAT", "NTP", is, s.ntp, "0..3");
  readValue("BAROSTAT", "NPVAR", is, s.npvar, ">=0");
  readValue("BAROSTAT", "NPBTH", is, npbth, ">=0");
  readValue("BAROSTAT", "COMP", is, s.comp, ">=0");
  if (npbth < 0) {
    std::stringstream ss;
    ss << npbth;
    printIO("BAROSTAT", "NPBTH", ss.str(), ">= 0");
  }
  for (int i = 0; i < npbth; ++i) {
    class ibarostat::pbath p;
    is >> dum;
    std::stringstream blockName;
    blockName << "PRSBTH[" << dum << "]";
    readValue("BAROSTAT", blockName.str(), is, p.prsbth, "0,1,2,3");
    if (s.npvar < 0) {
      std::stringstream ss;
      ss << npbth;
      printIO("BAROSTAT", "NPBTH", ss.str(), ">= 0");
    }
    for (int j = 0; j < s.npvar; ++j) {
      std::stringstream blockName;
      blockName << "TAUBBA[" << dum << "," << j + 1 << "]";
      std::stringstream ss;
      ss << taup;
      readValue("BAROSTAT", blockName.str(), is, taup, "> 0.0");
      p.taubba.push_back(taup);
    }
    s.pbaths.push_back(p);
  }
  for (int i = 0; i < 6; i++) {
    std::stringstream blockName;
    blockName << "NPCPL[" << i + 1 << "]";
    readValue("BAROSTAT", blockName.str(), is, s.npcpl[i], "> 0.0");
  }
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of BAROSTAT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iboundcond &s) {
  s.found = 1;
  readValue("BOUNDCOND", "NTB", is, s.ntb, "-1..2");
  readValue("BOUNDCOND", "NDFMIN", is, s.ndfmin, ">= 0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of BOUNDCOND block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ibsleus &s) {
  s.found = 1;
  readValue("BSLEUS", "BSLEUS", is, s.bsleus, "0 or 1");
  readValue("BSLEUS", "BUILD", is, s.build, "0 or 1");
  readValue("BSLEUS", "WRITE", is, s.write, ">= 0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of BSLEUS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}


std::istringstream & operator>>(std::istringstream &is, icolvarres &s) {
  s.found = 1;
  readValue("COLVARRES", "CVR", is, s.cvr, "0,1");
  readValue("COLVARRES", "CVK", is, s.cvk, ">=0");
  readValue("COLVARRES", "TAUCVR", is, s.taucvr, ">=0");
  readValue("COLVARRES", "VCVR", is, s.vcvr, "0,1");
  readValue("COLVARRES", "NTWCV", is, s.ntwcv, ">=0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of COLVARRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, icgrain &s) {
  s.found = 1;
  readValue("CGRAIN", "NTCGRAN", is, s.ntcgran, "0,1,2");
  readValue("CGRAIN", "EPS", is, s.eps, ">= 0.0");
  readValue("CGRAIN", "EPSM", is, s.epsm, ">= 0.0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of CGRAIN block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, icomtransrot &s) {
  s.found = 1;
  readValue("COMTRANSROT", "NSCM", is, s.nscm, "integers");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of COMTRANSROT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iconsistencycheck &s) {
  int nackf, nckf;
  s.found = 1;
  readValue("CONSISTENCYCHECK", "NTCHK", is, s.ntchk, "0, 1");
  readValue("CONSISTENCYCHECK", "NTCKF", is, s.ntckf, "0, 1");
  readValue("CONSISTENCYCHECK", "FDCKF", is, s.fdckf, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKV", is, s.ntckv, ">0.0");
  readValue("CONSISTENCYCHECK", "FDCKV", is, s.fdckv, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKT", is, s.ntckt, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKE", is, s.ntcke, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKR", is, s.ntckr, ">0.0");
  readValue("CONSISTENCYCHECK", "NTCKL", is, s.ntckl, ">0.0");
  readValue("CONSISTENCYCHECK", "FDCKL", is, s.fdckl, ">0.0");
  readValue("CONSISTENCYCHECK", "NACKF", is, nackf, ">0.0");
  if (nackf < 0) {
    std::stringstream ss;
    ss << nackf;
    printIO("CONSISTENCYCHECK", "NACKF", ss.str(), ">= 0");
  }
  for (int i = 0; i < nackf; i++) {
    readValue("CONSISTENCYCHECK", "NCKF", is, nckf, ">= 1");
    s.nckf.push_back(nckf);
  }
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of CONSISTENCYCHECK block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iconstraint &s) {
  s.found = 1;
  readValue("CONSTRAINT", "NTC", is, s.ntc, "0..4");
  readValue("CONSTRAINT", "NTCP", is, s.ntcp, "1..3");
  readValue("CONSTRAINT", "NTCP0(1)", is, s.ntcp0[0], ">=0");
  if(s.ntcp == 3) {
    readValue("CONSTRAINT", "NTCP0(2)", is, s.ntcp0[1], ">=0");
    readValue("CONSTRAINT", "NTCP0(3)", is, s.ntcp0[2], ">=0");
  }
  readValue("CONSTRAINT", "NTCS", is, s.ntcs, "1..6");
  if (s.ntcs != 4) readValue("CONSTRAINT", "NTCS0(1)", is, s.ntcs0[0], ">=0");
  if (s.ntcs == 3) {
    readValue("CONSTRAINT", "NTCS0(2)", is, s.ntcs0[1], ">=0");
    readValue("CONSTRAINT", "NTCS0(3)", is, s.ntcs0[2], ">=0");
  }
  if (s.ntcs == 6){
    readValue("CONSTRAINT","NTCG", is, s.ntcg, ">0");
    s.ntcd.resize(s.ntcg, -1);
    for (int g = 0; g < s.ntcg; g++) {
        readValue("CONSTRAINT", "NTCD", is, s.ntcd[g], ">=-1", true);
    }
  }
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of CONSTRAINT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, icovalentform &s) {
  s.found = 1;
  readValue("COVALENTFORM", "NTBBH", is, s.ntbbh, "0,1");
  readValue("COVALENTFORM", "NTBAH", is, s.ntbah, "0,1");
  readValue("COVALENTFORM", "NTBDN", is, s.ntbdn, "0,1");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of COVALENTFORM block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, idebug &s) {
  int nrd;
  s.found = 1;
  readValue("DEBUG", "NRDEO", is, nrd, ">=0");
  if (nrd < 0) {
    std::stringstream ss;
    ss << nrd;
    printIO("DEBUG", "NRD", ss.str(), ">= 0");
  }
  for (int i = 0; i < nrd; ++i) {
    class idebug::routine r;
    std::stringstream ss;
    if (!(is >> r.piider)) {
      std::stringstream msg;
      msg << "PIIDER(" << i + 1 << ")";
      printIO("DEBUG", msg.str(), r.piider, " a string");
    }
    ss.clear();
    ss.str("");
    std::stringstream blockName;
    blockName << "IIIDEO(" << i + 1 << ")";
    readValue("DEBUG", blockName.str(), is, r.iiideo, ">=0");
    s.routines.push_back(r);
  }
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of DEBUG block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, idihedralres &s) {
  s.found = 1;
  readValue("DIHEDRALRES", "NTDLR", is, s.ntdlr, "0..3");
  readValue("DIHEDRALRES", "CDLR", is, s.cdlr, ">=0.0");
  readValue("DIHEDRALRES", "PHILIN", is, s.philin, "-1..1");
  readValue("DIHEDRALRES", "VDIH", is, s.vdih, "0,1");
  readValue("DIHEDRALRES", "NTWDLR", is, s.ntwdlr, ">=0");
  readValue("DIHEDRALRES", "TOLDAC", is, s.toldac, ">=0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of DIHEDRALRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iangleres &s) {
  s.found = 1;
  readValue("ANGLERES", "NTALR", is, s.ntalr, "0..3");
  readValue("ANGLERES", "CALR", is, s.calr, ">=0.0");
  readValue("ANGLERES", "VARES", is, s.vares, "0,1");
  readValue("ANGLERES", "NTWALR", is, s.ntwalr, ">=0");
  readValue("ANGLERES", "TOLBAC", is, s.tolbac, ">=0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of ANGLERES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, idistancefield &s){
  s.found = 1;
  readValue("DISTANCEFIELD", "NTDFR", is, s.ntdfr, "0,1");
  readValue("DISTANCEFIELD", "GRID", is, s.grid, ">0.0");
  readValue("DISTANCEFIELD", "PROTEINOFFSET", is, s.proteinoffset, ">0.0");
  readValue("DISTANCEFIELD", "PROTEINCUTOFF", is, s.proteincutoff, ">0.0");
  readValue("DISTANCEFIELD", "PROTECT", is, s.protect, ">=0");
  readValue("DISTANCEFIELD", "UPDATE", is, s.update, ">0");
  readValue("DISTANCEFIELD", "SMOOTH", is, s.smooth, ">=0");
  readValue("DISTANCEFIELD", "RL", is, s.rl, ">=0");
  readValue("DISTANCEFIELD", "NTWDF", is, s.ntwdf, ">=0");
  readValue("DISTANCEFIELD", "PRINTGRID", is, s.printgrid, "0,1");
  
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of DISTANCEFIELD block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, idistanceres &s) {
  s.found = 1;
  readValue("DISTANCERES", "NTDIR", is, s.ntdir, "-2..3");
  readValue("DISTANCERES", "NTDIRA", is, s.ntdira, "0,1");
  readValue("DISTANCERES", "CDIR", is, s.cdir, ">=0.0");
  readValue("DISTANCERES", "DIR0", is, s.dir0, ">=0.0");
  readValue("DISTANCERES", "TAUDIR", is, s.taudir, ">=0.0");
  readValue("DISTANCERES", "FORCESCALE", is, s.forcescale, "0..2");
  readValue("DISTANCERES", "VDIR", is, s.vdir, "0,1");
  readValue("DISTANCERES", "NTWDIR", is, s.ntwdir, ">=0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of DISTANCERES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ieds &s) {
  s.found = 1;
  readValue("EDS", "EDS", is, s.eds, "0,1");
  readValue("EDS", "FORM", is, s.form, "1..3");
  readValue("EDS", "NUMSTATES", is, s.numstates, ">1");
  if (s.numstates <= 1) {
    std::stringstream ss;
    ss << s.numstates;
    printIO("EDS", "NUMSTATES", ss.str(), ">1");
  }
  switch (s.form) {
    case 1:
      s.smooth.resize(1);
      readValue("EDS", "S", is, s.smooth[0], ">0.0");
      break;
    case 2:
      s.smooth.resize(s.numstates*(s.numstates-1)/2);
      for(int i = 0; i < s.numstates*(s.numstates-1)/2; i++) {
        std::stringstream blockName;
        blockName << "S[" << i + 1 << "]";
        readValue("EDS", blockName.str(), is, s.smooth[i], ">0.0");
      }
      break;
    case 3:
      s.smooth.resize(s.numstates-1);
      s.tree.resize(s.numstates-1);
      for(int N = 0; N < s.numstates-1; N++) {
        std::stringstream blockName;
        blockName << "i[" << N + 1 << "]";
        std::vector<int> pair(2);
        readValue("EDS", blockName.str(), is, pair[0], ">0");
        blockName.str("");
        blockName.clear();
        blockName << "j[" << N + 1 << "]";
        readValue("EDS", blockName.str(), is, pair[1], ">0");
        s.tree[N] = pair;
        blockName.str("");
        blockName.clear();
        blockName << "S[" << N + 1 << "]";
        readValue("EDS", blockName.str(), is, s.smooth[N], ">0.0");
      }
      break;
    default:
      std::stringstream ss;
      ss >> s.numstates;
      printIO("EDS", "NUMSTATES", ss.str(), ">1");
      break;
  }
  s.eir.resize(s.numstates);
  for (int N = 0; N < s.numstates; N++) {
    std::stringstream blockName;
    blockName << "EIR[" << N + 1 << "]";
    readValue("EDS", blockName.str(), is, s.eir[N], ">0.0");
  }
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of EDS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ienergymin &s) {
  s.found = 1;
  readValue("ENERGYMIN", "NTEM", is, s.ntem, "0..3");
  readValue("ENERGYMIN", "NCYC", is, s.ncyc, ">0");
  readValue("ENERGYMIN", "DELE", is, s.dele, ">0.0");
  readValue("ENERGYMIN", "DX0", is, s.dx0, ">0.0");
  readValue("ENERGYMIN", "DXM", is, s.dxm, ">0.0");
  readValue("ENERGYMIN", "NMIN", is, s.nmin, ">0");
  readValue("ENERGYMIN", "FLIM", is, s.flim, ">0.0");
  if (s.ntem > 1) {
    readValue("ENERGYMIN", "CGIM", is, s.cgim, ">0");
    readValue("ENERGYMIN", "CGIC", is, s.cgic, ">=0.0");
  }
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of ENERGYMIN block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iewarn &s) {
  s.found = 1;
  readValue("EWARN", "MAXENER", is, s.maxener, "a double");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of EWARN block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iforce &s) {
  s.found = 1;
  int negr, nre;
  
  // new FORCE block
  for (int i = 0; i < 6; i++) {
    readValue("FORCE", "NTF", is, s.ntf[i], "0,1");
  }
  
  readValue("FORCE", "NEGR", is, negr, "!=0");
  if (negr < 0) {
    s.force_groups = true;
    negr = -negr;
  }

  for (int i = 0; i < negr; i++) {
    std::stringstream blockName;
    blockName << "NRE(" << i + 1 << ")";
    readValue("FORCE", blockName.str(), is, nre, ">0");
    s.nre.push_back(nre);
  }
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of FORCE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, igamd &s) {
  s.found = 1;
  readValue("GAMD", "GAMD", is, s.gamd, "0,1");
  readValue("GAMD", "SEARCH", is, s.search, "0..2");
  readValue("GAMD", "FORM", is, s.form, "1..3");
  readValue("GAMD", "THRESH", is, s.thresh, "1,2");
  readValue("GAMD", "NTIGAMDS", is, s.ntigamds, "0,1");
  readValue("GAMD", "AGROUPS", is, s.agroups, ">1");
  readValue("GAMD", "IGROUPS", is, s.igroups, ">1");
  readValue("GAMD", "DIHSTD", is, s.dihstd, ">0.0");
  readValue("GAMD", "TOTSTD", is, s.totstd, ">0.0");

  s.ed.resize(s.igroups);
  for (int N = 0; N < s.igroups; N++) {
    std::stringstream blockName;
    blockName << "ED[" << N + 1 << "]";
    readValue("GAMD", blockName.str(), is, s.ed[N], ">0.0");
  }
  s.et.resize(s.igroups);
  for (int N = 0; N < s.igroups; N++) {
    std::stringstream blockName;
    blockName << "ET[" << N + 1 << "]";
    readValue("GAMD", blockName.str(), is, s.et[N], ">0.0");
  }
  s.kd.resize(s.igroups);
  for (int N = 0; N < s.igroups; N++) {
    std::stringstream blockName;
    blockName << "KD[" << N + 1 << "]";
    readValue("GAMD", blockName.str(), is, s.kd[N], ">0.0");
  }
  s.kt.resize(s.igroups);
  for (int N = 0; N < s.igroups; N++) {
    std::stringstream blockName;
    blockName << "KT[" << N + 1 << "]";
    readValue("GAMD", blockName.str(), is, s.kt[N], ">0.0");
  }
  readValue("GAMD", "EQSTEPS", is, s.eqsteps, ">0");
  readValue("GAMD", "WINDOW", is, s.window, ">0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of EDS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, igeomconstraints &s) {
  s.found = 1;
  readValue("GEOMCONSTRAINTS", "NTCPH", is, s.ntcph, "0,1");
  readValue("GEOMCONSTRAINTS", "NTCPN", is, s.ntcpn, "0,1");
  readValue("GEOMCONSTRAINTS", "NTCS", is, s.ntcs, "0,1");
  readValue("GEOMCONSTRAINTS", "SHKTOL", is, s.shktol, ">0.0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of GEOMCONSTRAINTS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, igromos96compat &s) {
  s.found = 1;
  readValue("GROMOS96COMPAT", "NTNB96", is, s.ntnb96, "0,1");
  readValue("GROMOS96COMPAT", "NTR96", is, s.ntr96, "0,1");
  readValue("GROMOS96COMPAT", "NTP96", is, s.ntp96, "0,1");
  readValue("GROMOS96COMPAT", "NTG96", is, s.ntg96, "0,1");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of GROMOS96COMPAT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iinitialise &s) {
  s.found = 1;
  readValue("INITIALISE", "NTIVEL", is, s.ntivel, "0,1");
  readValue("INITIALISE", "NTISHK", is, s.ntishk, "0..3");
  readValue("INITIALISE", "NTINHT", is, s.ntinht, "0,1");
  readValue("INITIALISE", "NTINHB", is, s.ntinhb, "0,1");
  readValue("INITIALISE", "NTISHI", is, s.ntishi, "0,1");
  readValue("INITIALISE", "NTIRTC", is, s.ntirtc, "0,1");
  readValue("INITIALISE", "NTICOM", is, s.nticom, "0..3");
  readValue("INITIALISE", "NTISTI", is, s.ntisti, ">0");
  readValue("INITIALISE", "IG", is, s.ig, "0,1");
  readValue("INITIALISE", "TEMPI", is, s.tempi, ">0.0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of INITIALISE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iinnerloop &s) {
  s.found = 1;
  readValue("INNERLOOP", "NTILM", is, s.ntilm, "0..4");
  readValue("INNERLOOP", "NTILS", is, s.ntils, "0,1");
  // if CUDA then read number gpus and device ids
  // if no device ids are given == -1 (CUDA driver will decide)
  if (s.ntilm == 4) {
    readValue("INNERLOOP", "NGPUS", is, s.ngpus, ">0");
    s.ndevg.resize(s.ngpus, -1);
    for (int g = 0; g < s.ngpus; g++) {
        readValue("INNERLOOP", "NDEVG", is, s.ndevg[g], ">=-1");
    }
  }

  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of INNERLOOP block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iintegrate &s) {
  s.found = 1;
  readValue("INTEGRATE", "NINT", is, s.nint, "0,1");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of INTEGRATE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ijvalueres &s) {
  s.found = 1;
  readValue("JVALRES", "NTJVR", is, s.ntjvr, "-3..2");
  readValue("JVALRES", "NTJVRA", is, s.ntjvra, "0..1");
  readValue("JVALRES", "CJVR", is, s.cjvr, ">=0");
  readValue("JVALRES", "TAUJVR", is, s.taujvr, ">=0");
  readValue("JVALRES", "NJVRTARS", is, s.njvrtars, "0,1");
  readValue("JVALRES", "NJVRBIQW", is, s.njvrbiqw, "0..2");
  readValue("JVALRES", "LE", is, s.le, "0..1");
  readValue("JVALRES", "NGRID", is, s.ngrid, ">0");
  readValue("JVALRES", "DELTA", is, s.delta, ">0.0");
  readValue("JVALRES", "NTWJV", is, s.write, ">=0");
  std::string st;
  if(is.eof() == false){
    is >> st;
    if(st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of JVALRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ilambdas &s) {
  s.found = 1;
  readValue("LAMBDAS", "NTIL", is, s.ntil, "0,1");
  int i = 0;
  std::string dum;
  while ((is >> dum)) {
    std::istringstream ss;
    ss.str(dum);
    i++;
    class ilambdas::lambint l;
    std::stringstream blockName;
    blockName << "NTLI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), ss, l.ntli, "1..13");
    blockName.str("");
    blockName << "NILG1[" << i << "]";
    readValue("LAMBDAS", blockName.str(), is, l.nilg1, ">0");
    blockName.str("");
    blockName << "NILG2[" << i << "]";
    readValue("LAMBDAS", blockName.str(), is, l.nilg2, ">0");
    blockName.str("");
    blockName << "ALI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), is, l.ali, "a double");
    blockName.str("");
    blockName << "BLI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), is, l.bli, "a double");
    blockName.str("");
    blockName << "CLI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), is, l.cli, "a double");
    blockName.str("");
    blockName << "DLI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), is, l.dli, "a double");
    blockName.str("");
    blockName << "ELI[" << i << "]";
    readValue("LAMBDAS", blockName.str(), is, l.eli, "a double");

    s.lambints.push_back(l);
  }
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of LAMBDAS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ilocalelev &s) {
  s.found = 1;
  readValue("LOCALELEV", "NTLES", is, s.ntles, "0..5");
  readValue("LOCALELEV", "NLEPOT", is, s.nlepot, ">=0");
  readValue("LOCALELEV", "NTLESA", is, s.ntlesa, "0..2");
  readValue("LOCALELEV", "NTWLE", is, s.ntwle, ">=0");
  if (s.nlepot < 0) {
    std::stringstream ss;
    ss << s.nlepot;
    printIO("LOCALELEV", "NLEPOT", ss.str(), ">=0");
  }
  int nlepid, ntlepfr;
  std::string s_nlepid, s_ntlepfr;
  for (int i = 0; i < s.nlepot; i++) {
    std::stringstream blockName;
    blockName << "NLEPID(" << i + 1 << ")";
    readValue("LOCALELEV", blockName.str(), is, nlepid, "1..NLEPOT");
    blockName.str("");
    blockName << "NTLEPFR(" << i + 1 << ")";
    readValue("LOCALELEV", blockName.str(), is, ntlepfr, "0,1");
    s.nlepid_ntlepfr.insert(std::pair<int, int> (nlepid, ntlepfr));
  }
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of LOCALELEV block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ielectric &s) {
  s.found = 1;
  readValue("ELECTRIC", "FIELD", is, s.field, "0..1");
  readValue("ELECTRIC", "DIPOLE", is, s.dipole, "0..1");
  readValue("ELECTRIC", "CURRENT", is, s.current, "0..1");
  readValue("ELECTRIC", "EF_x", is, s.ef_x, "double");
  readValue("ELECTRIC", "EF_y", is, s.ef_y, "double");
  readValue("ELECTRIC", "EF_z", is, s.ef_z, "double");
  if (s.field < 0 || s.field > 1) {
    std::stringstream ss;
    ss << s.field;
    printIO("ELECTRIC", "FIELD", ss.str(), "0..1");
  }
  if (s.dipole < 0 || s.dipole > 1) {
    std::stringstream ss;
    ss << s.dipole;
    printIO("ELECTRIC", "DIPOLE", ss.str(), "0..1");
  }
  if (s.current < 0 || s.current > 1) {
    std::stringstream ss;
    ss << s.current;
    printIO("ELECTRIC", "CURRENT", ss.str(), "0..1");
  }
  readValue("ELECTRIC", "DIPGRP", is, s.dipgrp, "0..2");
  readValue("ELECTRIC", "NTWDIP", is, s.ntwdip, ">=0");
  readValue("ELECTRIC", "NTWCUR", is, s.ntwcur, ">=0");
  readValue("ELECTRIC", "NCURGRP", is, s.ncurgrp, ">=0");

  for (int i = 0; i < s.ncurgrp; ++i) {
    int grp;
    std::stringstream blockName;
    blockName << "CURGRP[" << i + 1 << "]";
    readValue("ELECTRIC", blockName.str(), is, grp, ">=0");
    blockName.str("");
    s.curgrp.push_back(grp);
  }
  
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of ELECTRIC block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, imultibath &s) {
  s.found = 1;
  readValue("MULTIBATH", "NTBTYP", is, s.ntbtyp, "0..2");
  if(s.ntbtyp == 2) {
    readValue("MULTIBATH", "NUM", is, s.num, ">=0");
  }
  readValue("MULTIBATH", "NBATHS", is, s.nbaths, ">=0");
  if (s.nbaths < 0) {
    std::stringstream ss;
    ss << s.nbaths;
    printIO("MULTIBATH", "NBATHS", ss.str(), ">0");
  }
  for (int i = 0; i < s.nbaths; ++i) {
    double temp0, tau;
    std::stringstream blockName;
    blockName << "TEMP[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, temp0, ">=0.0");
    blockName.str("");
    blockName << "TAU[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, tau, ">=0.0");
    s.temp0.push_back(temp0);
    s.tau.push_back(tau);
  }
  readValue("MULTIBATH", "DOFSET", is, s.dofset, ">=0");
  if (s.dofset < 0) {
    std::stringstream ss;
    ss << s.dofset;
    printIO("MULTIBATH", "DOFSET", ss.str(), ">=0");
  }
  for (int i = 0; i < s.dofset; ++i) {
    int last, combath, irbath;
    std::stringstream blockName;
    blockName << "LAST[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, last, ">=0");
    blockName.str("");
    blockName << "COM-BATH[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, combath, ">=1");
    blockName.str("");
    blockName << "IR-BATH[" << i + 1 << "]";
    readValue("MULTIBATH", blockName.str(), is, irbath, ">=1");
    s.last.push_back(last);
    s.combath.push_back(combath);
    s.irbath.push_back(irbath);
  }
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of MULTIBATH block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, imulticell &s) {
  s.found = 1;
  readValue("MULTICELL", "NTM", is, s.ntm, "0,1");
  readValue("MULTICELL", "NCELLA", is, s.ncella, ">=1");
  readValue("MULTICELL", "NCELLB", is, s.ncellb, ">=1");
  readValue("MULTICELL", "NCELLC", is, s.ncellc, ">=1");
  readValue("MULTICELL", "TOLPX", is, s.tolpx, ">=0.0");
  readValue("MULTICELL", "TOLPV", is, s.tolpv, ">=0.0");
  readValue("MULTICELL", "TOLPF", is, s.tolpf, ">=0.0");
  readValue("MULTICELL", "TOLPFW", is, s.tolpfw, ">=0.0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of MULTICELL block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, imultigradient &s) {
  s.found = 1;
  readValue("MULTIGRADIENT", "NTMGRE", is, s.ntmgre, "0,1");
  readValue("MULTIGRADIENT", "NTMGRP", is, s.ntmgrp, "0..3");
  readValue("MULTIGRADIENT", "NTMGRN", is, s.ntmgrn, ">=0");
  if (s.ntmgrn < 0) {
    std::ostringstream os;
    os << "MULTIGRADIENT block: negative number of gradients.";
    printError(os.str());
    return is;
  }

  s.mgrvar.resize(s.ntmgrn, "");
  s.mgrfrm.resize(s.ntmgrn);
  s.curves.resize(s.ntmgrn);
  
  for(int i = 0; i < s.ntmgrn; ++i) {
    int num;
    is >> s.mgrvar[i] >> s.mgrfrm[i] >> num;
    if (is.fail()) {
      printError("MULTIGRADIENT block: cannot read MGRVAR, MGRFRM or MGRNCP.");
      return is;
    }
    if (num < 2) {
      printError("MULTIGRADIENT block: MGRNCP has to be >= 2.");
      return is;
    }
    s.curves[i].resize(num, std::pair<double, double>(0.0, 0.0));
    for(int j = 0; j < num; ++j) {
      is >> s.curves[i][j].first >> s.curves[i][j].second;
      if (is.fail()) {
        printError("MULTIGRADIENT block: cannot read MGRCPT or MGRCPV.");
        return is;
      }
    }
  }
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of MULTIGRADIENT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}


std::istringstream & operator>>(std::istringstream &is, imultistep &s) {
  s.found = 1;
  readValue("MULTISTEP", "STEPS", is, s.steps, ">0");
  readValue("MULTISTEP", "BOOST", is, s.boost, "0,1");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of MULTISTEP block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ineighbourlist &s) {
  s.found = 1;
  readValue("NEIGHBOURLIST", "PLALGO", is, s.plalgo, "0..2");
  readValue("NEIGHBOURLIST", "NUPDPL", is, s.nupdpl, ">0");
  readValue("NEIGHBOURLIST", "NUPDIS", is, s.nupdis, ">0");
  readValue("NEIGHBOURLIST", "NUPDII", is, s.nupdii, ">0");
  readValue("NEIGHBOURLIST", "RCUTS", is, s.rcuts, ">=0");
  readValue("NEIGHBOURLIST", "RCUTI", is, s.rcuti, ">=RCUTS");
  readValue("NEIGHBOURLIST", "GRIDSZX", is, s.gridszx, ">=0");
  readValue("NEIGHBOURLIST", "GRIDSZY", is, s.gridszy, ">=0");
  readValue("NEIGHBOURLIST", "GRIDSZZ", is, s.gridszz, ">=0");
  readValue("NEIGHBOURLIST", "TYPE", is, s.type, "0..1");
  readValue("NEIGHBOURLIST", "NCGCEN", is, s.ncgcen, ">=-2");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of NEIGHBOURLIST block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, inemd &s) {
  s.found = 1;
  readValue("NEMD", "NEMD", is, s.nemd, "0, 1");
  readValue("NEMD", "PROPERTY", is, s.property, "0");
  readValue("NEMD", "METHOD", is, s.method, "0, 1");
  readValue("NEMD", "SLABNUM", is, s.slabnum, ">=1");
  readValue("NEMD", "PERTFRQ", is, s.pertfrq, ">=1");
  readValue("NEMD", "AMPLI", is, s.ampli, ">=0");
  readValue("NEMD", "STDYAFT", is, s.stdyaft, ">=0");
  readValue("NEMD", "WRITE", is, s.write, ">=0");
  
  return is;
}

std::istringstream & operator>>(std::istringstream &is, inonbonded &s) {
  s.found = 1;
  readValue("NONBONDED", "NLRELE", is, s.nlrele, "-4..4");
  readValue("NONBONDED", "APPAK", is, s.appak, ">=0.0");
  readValue("NONBONDED", "RCRF", is, s.rcrf, ">=0.0");
  readValue("NONBONDED", "EPSRF", is, s.epsrf, ">=0.0");
  readValue("NONBONDED", "NSLFEXCL", is, s.nslfexcl, "(0,1)");
  readValue("NONBONDED", "NSHAPE", is, s.nshape, "-1..10");
  readValue("NONBONDED", "ASHAPE", is, s.ashape, ">0.0");
  readValue("NONBONDED", "NA2CLC", is, s.na2clc, "0..4");
  readValue("NONBONDED", "TOLA2", is, s.tola2, ">0.0");
  readValue("NONBONDED", "EPSLS", is, s.epsls, ">0.0");
  readValue("NONBONDED", "NKX", is, s.nkx, ">0");
  readValue("NONBONDED", "NKY", is, s.nky, ">0");
  readValue("NONBONDED", "NKZ", is, s.nkz, ">0");
  readValue("NONBONDED", "KCUT", is, s.kcut, ">0.0");
  readValue("NONBONDED", "NGX", is, s.ngx, ">0");
  readValue("NONBONDED", "NGY", is, s.ngy, ">0");
  readValue("NONBONDED", "NGZ", is, s.ngz, ">0");
  readValue("NONBONDED", "NASORD", is, s.nasord, "1..5");
  readValue("NONBONDED", "NFDORD", is, s.nfdord, "0..5");
  readValue("NONBONDED", "NALIAS", is, s.nalias, ">0");
  readValue("NONBONDED", "NSPORD", is, s.nspord, ">0");
  readValue("NONBONDED", "NQEVAL", is, s.nqeval, ">=0");
  readValue("NONBONDED", "FACCUR", is, s.faccur, ">0.0");
  readValue("NONBONDED", "NRDGRD", is, s.nrdgrd, "0,1");
  readValue("NONBONDED", "NWRGRD", is, s.nwrgrd, "0,1");
  readValue("NONBONDED", "NLRLJ", is, s.nlrlj, "0,1");
  readValue("NONBONDED", "SLVDNS", is, s.slvdns, ">0.0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of NONBONDED block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iorderparamres &s) {
  s.found = 1;
  readValue("ORDERPARAMRES", "NTOPR", is, s.ntopr, "-2..2");
  readValue("ORDERPARAMRES", "NTOPRA", is, s.ntopra, "0,1");
  readValue("ORDERPARAMRES", "COPR", is, s.copr, ">=0.0");
  readValue("ORDERPARAMRES", "TAUOPR", is, s.tauopr, ">=0.0");
  readValue("ORDERPARAMRES", "UPDOPR", is, s.updopr, "> 0");
  readValue("ORDERPARAMRES", "NTWOP", is, s.ntwop, ">=0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of ORDERPARAMRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ioveralltransrot &s) {
  s.found = 1;
  readValue("OVERALLTRANSROT", "NCMTR", is, s.ncmtr, "0,1");
  readValue("OVERALLTRANSROT", "NCMRO", is, s.ncmro, "0,1");
  readValue("OVERALLTRANSROT", "CMAMX", is, s.cmamx, ">=0.0");
  readValue("OVERALLTRANSROT", "CMAMY", is, s.cmamy, ">=0.0");
  readValue("OVERALLTRANSROT", "CMAMZ", is, s.cmamz, ">=0.0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of OVERALLTRANSROT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ipairlist &s) {
  s.found = 1;
  readValue("PAIRLIST", "ALGORITHM", is, s.algorithm, "0,1,2");
  readValue("PAIRLIST", "NSNB", is, s.nsnb, ">0");
  readValue("PAIRLIST", "RCUTP", is, s.rcutp, ">0.0");
  readValue("PAIRLIST", "RCUTL", is, s.rcutl, ">0.0");
  readValue("PAIRLIST", "SIZE", is, s.size, ">0.0");
  readValue("PAIRLIST", "TYPE", is, s.type, "0,1");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of PAIRLIST block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ipathint &s) {
  s.found = 1;
  readValue("PATHINT", "NTPI", is, s.ntpi, "0,1");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of PATHINT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iperscale &s) {
  s.found = 1;
  readValue("PERSCALE", "RESTYPE", is, s.restype, "0,1");
  readValue("PERSCALE", "KDIH", is, s.kdih, ">=0.0");
  readValue("PERSCALE", "KJ", is, s.kj, ">=0.0");
  readValue("PERSCALE", "T", is, s.t, ">0.0");
  readValue("PERSCALE", "DIFF", is, s.diff, ">=0.0");
  readValue("PERSCALE", "RATIO", is, s.ratio, ">0.0");
  readValue("PERSCALE", "READ", is, s.read, "0,1");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of PERSCALE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iperturbation &s) {
  s.found = 1;
  readValue("PERTURBATION", "NTG", is, s.ntg, "0,1");
  readValue("PERTURBATION", "NRDGL", is, s.nrdgl, "0,1");
  readValue("PERTURBATION", "RLAM", is, s.rlam, "0.0,..,1.0");
  readValue("PERTURBATION", "DLAMT", is, s.dlamt, ">=0.0");
  readValue("PERTURBATION", "ALPHLJ", is, s.alphlj, ">=0.0");
  readValue("PERTURBATION", "ALPHC", is, s.alphc, ">=0.0");
  readValue("PERTURBATION", "NLAM", is, s.nlam, ">0");
  readValue("PERTURBATION", "NSCALE", is, s.nscale, "0..2");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of PERTURBATION block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ipolarise &s) {
  s.found = 1;
  readValue("POLARISE", "COS", is, s.cos, "0,2");
  readValue("POLARISE", "EFIELD", is, s.efield, "0,1");
  readValue("POLARISE", "MINFIELD", is, s.minfield, ">=0.0");
  readValue("POLARISE", "DAMP", is, s.damp, "0,1");
  readValue("POLARISE", "WRITE", is, s.write, ">=0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of POLARISE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ipositionres &s) {
  s.found = 1;
  readValue("POSITIONRES", "NTPOR", is, s.ntpor, "0..3");
  readValue("POSITIONRES", "NTPORB", is, s.ntporb, "0,1");
  readValue("POSITIONRES", "NTPORS", is, s.ntpors, "0,1");
  readValue("POSITIONRES", "CPOR", is, s.cpor, ">=0.0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of POSITIONRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iprecalclam &s) {
  s.found = 1;
  readValue("PRECALCLAM", "NRLAM", is, s.nrlam, ">=0");
  readValue("PRECALCLAM", "MINLAM", is, s.minlam, ">=0.0,<MAXLAM");
  readValue("PRECALCLAM", "MAXLAM", is, s.maxlam, ">MINLAM,<=1.0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of PRECALCLAM block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ipressurescale &s) {
  s.found = 1;
  readValue("PRESSURESCALE", "COUPLE", is, s.couple, "0..2");
  readValue("PRESSURESCALE", "SCALE", is, s.scale, "0..4");
  readValue("PRESSURESCALE", "COMP", is, s.comp, ">0.0");
  readValue("PRESSURESCALE", "TAUP", is, s.taup, ">0.0");
  readValue("PRESSURESCALE", "VIRIAL", is, s.virial, "0..2");
  readValue("PRESSURESCALE", "X_SEMI", is, s.x_semi, "0..2");
  readValue("PRESSURESCALE", "Y_SEMI", is, s.y_semi, "0..2");
  readValue("PRESSURESCALE", "Z_SEMI", is, s.z_semi, "0..2");
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      std::stringstream blockName;
      blockName << "PRES0[" << i + 1 << "," << j + 1 << "]";
      readValue("PRESSURESCALE", blockName.str(), is, s.pres0[i][j], "a double");
      blockName.str("");
      blockName.clear();
    }
  }
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of PRESSURESCALE block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iprintout &s) {
  s.found = 1;
  readValue("PRINTOUT", "NTPR", is, s.ntpr, ">=0");
  readValue("PRINTOUT", "NTPP", is, s.ntpp, "0,1");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of PRINTOUT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iqmmm &s) {
  s.found = 1;
  readValue("QMMM", "NTQMMM", is, s.ntqmmm, "-1..3");
  readValue("QMMM", "NTQMSW", is, s.ntqmsw, "0..4");
  readValue("QMMM", "RCUTQM", is, s.rcutqm, "<==0.0 or !=0.0>");
  readValue("QMMM", "NTWQMMM", is, s.ntwqmmm, ">=0");
  readValue("QMMM", "QMLJ", is, s.qmlj, "0,1");
  readValue("QMMM", "QMCON", is, s.qmcon, "0,1");
  readValue("QMMM", "MMSCALE", is, s.mmscale, "<0.0 or >0.0", true);
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of QMMM block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, irandomnumbers &s) {
  s.found = 1;
  readValue("RANDOMNUMBERS", "NTRNG", is, s.ntrng, "0,1");
  readValue("RANDOMNUMBERS", "NTGSL", is, s.ntgsl, ">=-1");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of RANDOMNUMBERS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ireadtraj &s) {
  s.found = 1;
  readValue("READTRAJ", "NTRD", is, s.ntrd, "0,1");
  readValue("READTRAJ", "NTSTR", is, s.ntstr, "1..18");
  readValue("READTRAJ", "NTRB", is, s.ntrb, "0,1");
  readValue("READTRAJ", "NTSHK", is, s.ntshk, "0,1");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of READTRAJ block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ireplica &s) {
  s.found = 1;
  readValue("REPLICA", "RETL", is, s.nrel, "0,1");
  int nret;
  readValue("REPLICA", "NRET", is, nret, ">=0");
  if (nret < 0) {
    std::stringstream ss;
    ss << nret;
    printIO("REPLICA", "NRET", ss.str(), ">= 0");
  }
  for (int i = 0; i < nret; ++i) {
    std::stringstream blockName;
    blockName << "RET(" << i + 1 << ")";
    double ret;
    readValue("REPLICA", blockName.str(), is, ret, ">=0.0");
    s.ret.push_back(ret);
  }
  readValue("REPLICA", "LRESCALE", is, s.lrescale, "0,1");
  int nrelam;
  readValue("REPLICA", "NRELAM", is, nrelam, ">=0");
  if (nrelam < 0) {
    std::stringstream ss;
    ss << nrelam;
    printIO("REPLICA", "NRELAM", ss.str(), ">= 0");
  }
  for (int i = 0; i < nrelam; ++i) {
    std::stringstream blockName;
    blockName << "RELAM(" << i + 1 << ")";
    double relam;
    readValue("REPLICA", blockName.str(), is, relam, ">=0.0");
    s.relam.push_back(relam);
  //  blockName.str("");
  //  blockName.clear();
  }
  for (int i = 0; i < nrelam; ++i) {
    std::stringstream blockName;
    blockName << "RETS(" << i + 1 << ")";
    double rets;
    readValue("REPLICA", blockName.str(), is, rets, ">=0.0");
    s.rets.push_back(rets);
  }
  readValue("REPLICA", "NRETRIAL", is, s.nretrial, ">=0");
  readValue("REPLICA", "NREQUIL", is, s.nrequil, ">=0");
  readValue("REPLICA", "CONT", is, s.cont, "0,1");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of REPLICA block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ireeds &s) {
  s.found = 1;
  readValue("REPLICA_EDS", "REEDS", is, s.reeds, "0,1");
  int nres;
  readValue("REPLICA_EDS", "NRES", is, nres, ">=0");
  if (nres < 0) {
    std::stringstream ss;
    ss << nres;
    printIO("REPLICA_EDS", "NRES", ss.str(), ">= 0");
  }
  s.nres = nres;

  int numstates;
  readValue("REPLICA_EDS", "NUMSTATES", is, numstates, ">=0");
  if (numstates < 0) {
    std::stringstream ss;
    ss << numstates;
    printIO("REPLICA_EDS", "NUMSTATES", ss.str(), ">= 0");
  }
  s.numstates = numstates;

  int neoff;
  readValue("REPLICA_EDS", "NEOFF", is, neoff, ">=0");
  if (neoff < 0) {
    std::stringstream ss;
    ss << neoff;
    printIO("REPLICA_EDS", "NEOFF", ss.str(), ">= 0");
  }
  s.neoff = neoff;

  // Here add s-values
  for (int i = 0; i < nres; ++i) {
    std::stringstream blockName;
    blockName << "RES(1 ... NRES)";
    double res;
    readValue("REPLICA_EDS", blockName.str(), is, res, ">=0.0");
    s.res.push_back(res);
  }

  // here add the eoffs
  for (int i = 0; i < numstates*nres; ++i) {
    std::stringstream blockName;
    blockName << "EIR(NUMSTATES x NRES)";
    double eir;
    readValue("REPLICA_EDS", blockName.str(), is, eir, "");
    s.eir.push_back(eir);
  }

  // then the last line....
  readValue("REPLICA_EDS", "NRETRIAL", is, s.nretrial, ">=0");
  readValue("REPLICA_EDS", "NREQUIL", is, s.nrequil, "");
  readValue("REPLICA_EDS", "CONT", is, s.cont, ">=0");
  readValue("REPLICA_EDS", "EDS_STAT_OUT", is, s.cont, ">=0");
  readValue("REPLICA_EDS", "PERIODIC", is, s.cont, ">=0");

  // then to check end of block
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of REPLICA_EDS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }

  return is;
}

std::istringstream & operator>>(std::istringstream &is, irottrans &s) {
  s.found = 1;
  readValue("ROTTRANS", "RTC", is, s.rtc, "0,1");
  readValue("ROTTRANS", "RTCLAST", is, s.rtclast, ">0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of ROTTRANS block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, istep &s) {
  s.found = 1;
  readValue("STEP", "NSTLIM", is, s.nstlim, ">=0");
  readValue("STEP", "T", is, s.t, ">=0.0 or -1");
  readValue("STEP", "DT", is, s.dt, ">0.0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of STEP block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream &operator>>(std::istringstream &is, isasa &s) {
  s.found = 1;
  readValue("SASA", "NTSASA", is, s.ntsasa, "0,1");
  readValue("SASA", "NTVOL", is, s.ntvol, "0,1");
  readValue("SASA", "P_12", is, s.p12, ">0 and <1");
  readValue("SASA", "P_13", is, s.p13, ">0 and <1");
  readValue("SASA", "P_1X", is, s.p1x, ">0 and <1");
  readValue("SASA", "SIGMAV", is, s.sigmav, "a double");
  readValue("SASA", "RSOLV", is, s.rsolv, ">0");
  readValue("SASA", "AS1", is, s.as1, ">0");
  readValue("SASA", "AS2", is, s.as2, ">0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of SASA block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, istochdyn &s) {
  s.found = 1;
  readValue("STOCHDYN", "NTSD", is, s.ntsd, "0,1");
  readValue("STOCHDYN", "NTFR", is, s.ntfr, "0..3");
  readValue("STOCHDYN", "NSFR", is, s.nsfr, ">0");
  readValue("STOCHDYN", "NBREF", is, s.nbref, ">0");
  readValue("STOCHDYN", "RCUTF", is, s.rcutf, ">=0.0");
  readValue("STOCHDYN", "CFRIC", is, s.cfric, ">=0.0");
  readValue("STOCHDYN", "TEMPSD", is, s.tempsd, ">=0.0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of STOCHDYN block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, isymres &s) {
  s.found = 1;
  readValue("SYMRES", "NTSYM", is, s.ntsym, ">=0");
  readValue("SYMRES", "CSYM", is, s.csym, ">=0.0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of SYMRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, isystem &s) {
  s.found = 1;
  readValue("SYSTEM", "NPM", is, s.npm, ">=0");
  readValue("SYSTEM", "NSM", is, s.nsm, ">=0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of SYSTEM block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ithermostat &s) {
  s.found = 1;
  readValue("THERMOSTAT", "NTT", is, s.ntt, "0,1");
  readValue("THERMOSTAT", "NTBTH", is, s.ntbth, ">=0");
  readValue("THERMOSTAT", "NTSET", is, s.ntset, ">=0");
  if (s.ntbth < 0) {
    std::stringstream ss;
    ss << s.ntbth;
    printIO("THERMOSTAT", "NTBTH", ss.str(), ">=0");
  }
  for (int i = 0; i < s.ntbth; ++i) {
    class ithermostat::tbath bath;
    std::stringstream blockName;
    blockName << "I(" << i + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, bath.index, "1..NTBATH");
    blockName.str("");
    blockName.clear();
    blockName << "NTBTYP(" << i + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, bath.ntbtyp, "0..3");
    blockName.str("");
    blockName.clear();
    blockName << "TEMBTH(" << i + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, bath.tembth, ">=0");
    blockName.str("");
    blockName.clear();
    blockName << "NTBVAR(" << i + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, bath.ntbvar, ">=0 and <=MAX_NTBVAR");
    if (bath.ntbvar < 0) {
      std::stringstream ss;
      ss << bath.ntbvar;
      printIO("THERMOSTAT", blockName.str(), ss.str(), ">=0 and <=MAX_NTBVAR");
    }
    for (int j = 0; j < bath.ntbvar; j++) {
      std::stringstream blockName;
      blockName << "TAUBTH[" << i + 1 << "]";// << j + 1 << "]";
      double tau;
      readValue("THERMOSTAT", blockName.str(), is, tau, ">=0.0");
      bath.taubth.push_back(tau);
    }
    s.baths.push_back(bath);
  }
  if (s.ntset < 0) {
    std::stringstream ss;
    ss << s.ntbth;
    printIO("THERMOSTAT", "NTSET", ss.str(), ">=0");
  }
  for (int j = 0; j < s.ntset; ++j) {
    class ithermostat::tdofgroup dofgroup;
    std::stringstream blockName;
    blockName << "NTSCPL(" << j + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, dofgroup.ntscpl, "1..NTSET");
    blockName.str("");
    blockName.clear();
    blockName << "NTSTYP(" << j + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, dofgroup.ntstyp, "0..2");
    blockName.str("");
    blockName.clear();
    blockName << "NTSCNS(" << j + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, dofgroup.ntscns, "0..1");
    blockName.str("");
    blockName.clear();
    blockName << "NTSGT(" << j + 1 << ")";
    readValue("THERMOSTAT", blockName.str(), is, dofgroup.ntsgt, "-2..0 or >=1");
    
    int upper_k = 0;
    if (dofgroup.ntsgt == -2) {
      upper_k = 2;
    } else if (dofgroup.ntsgt == -1 || dofgroup.ntsgt == 0) {
      upper_k = 0;
    } else { // NTSGT > 0
      upper_k = dofgroup.ntsgt;
    }

    for (int k = 0; k < upper_k; k++) {
      std::stringstream blockName;
      blockName << "NTSGTG(" << j + 1 << "," << k + 1 << ")";
      int ntsgtg;
      readValue("THERMOSTAT", blockName.str(), is, ntsgtg, "-2..0 or >=1");
      dofgroup.ntsgtg.push_back(ntsgtg);
    }
    s.dofgroups.push_back(dofgroup);
  }
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of THERMOSTAT block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iumbrella &s) {
  s.found = 1;
  readValue("UMBRELLA", "NTUS", is, s.ntus, "0,1");
  readValue("UMBRELLA", "USCST1", is, s.uscst1, ">=0");
  readValue("UMBRELLA", "USCST2", is, s.uscst2, ">=0");
  readValue("UMBRELLA", "USREF1", is, s.usref1, ">=0");
  readValue("UMBRELLA", "USREF2", is, s.usref2, ">=0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of UMBRELLA block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ivirial &s) {
  s.found = 1;
  readValue("VIRIAL", "NTV", is, s.ntv, "0,1");
  readValue("VIRIAL", "NTVG", is, s.ntvg, "0..3");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of VIRIAL block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ivirtualatom &s) {
  s.found = 1;
  readValue("VIRTUALATOM", "VIRT", is, s.virt, "0,1");
  readValue("VIRTUALATOM", "NUMVIRT", is, s.numvirt, ">=0");
  readValue("VIRTUALATOM", "LASTVIRT", is, s.lastvirt, ">=0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of VIRTUALATOM block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iwritetraj &s) {
  s.found = 1;
  readValue("WRITETRAJ", "NTWX", is, s.ntwx, "an integer");
  readValue("WRITETRAJ", "NTWSE", is, s.ntwse, ">=0");
  readValue("WRITETRAJ", "NTWV", is, s.ntwv, "an integer");
  readValue("WRITETRAJ", "NTWF", is, s.ntwf, "an integer");
  readValue("WRITETRAJ", "NTWE", is, s.ntwe, ">=0");
  readValue("WRITETRAJ", "NTWG", is, s.ntwg, ">=0");
  readValue("WRITETRAJ", "NTWB", is, s.ntwb, ">=0");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of WRITETRAJ block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, ixrayres &s) {
  s.found = 1;
  readValue("XRAYRES", "NTXR", is, s.ntxr, "-2..3");
  readValue("XRAYRES", "NTXLE", is, s.ntxle, "0,1");
  readValue("XRAYRES", "CXR", is, s.cxr, ">=0.0");
  readValue("XRAYRES", "NTWXR", is, s.ntwxr, ">=0");
  readValue("XRAYRES", "NTWDE", is, s.ntwde, "0..3");
  readValue("XRAYRES", "NTWXM", is, s.ntwxm, ">=0");
  readValue("XRAYRES", "CXTAU", is, s.cxtau, ">=0.0");
  readValue("XRAYRES", "RDAVG", is, s.rdavg, "0,1");
  std::string st;
  if (is.eof() == false) {
    is >> st;
    if (st != "" || is.eof() == false) {
      std::stringstream ss;
      ss << "unexpected end of XRAYRES block, read \"" << st << "\" instead of \"END\"";
      printError(ss.str());
    }
  }
  return is;
}

std::istringstream & operator>>(std::istringstream &is, iunknown &s) {
  std::string e = is.str();
  s.content = e.substr(0, e.find("END"));
  return is;
}

gio::Ginstream & operator>>(gio::Ginstream &is, input &gin) {
  std::vector<std::string> buffer;
  while (!is.stream().eof()) {
    is.getblock(buffer);

    if (buffer.size() == 1) { // for example if there is twice an "END" in a row
      std::stringstream msg;
      msg << buffer[0] << " instead of a block name was ignored.";
      printWarning(msg.str());
      buffer.pop_back();
    } else if (buffer.size() > 1) {
      std::string bufferstring;
      gio::concatenate(buffer.begin() + 1, buffer.end() - 1, bufferstring);
      std::istringstream bfstream(bufferstring);
      switch (BLOCKTYPE[buffer[0]]) {
        case addecoupleblock: bfstream >> gin.addecouple;
          break;
        case aedsblock: bfstream >> gin.aeds;
          break;
        case barostatblock: bfstream >> gin.barostat;
          break;
        case boundcondblock: bfstream >> gin.boundcond;
          break;
        case bsleusblock: bfstream >> gin.bsleus;
          break;
        case colvarresblock: bfstream >> gin.colvarres;
          break;
        case cgrainblock: bfstream >> gin.cgrain;
          break;
        case comtransrotblock: bfstream >> gin.comtransrot;
          break;
        case consistencycheckblock: bfstream >> gin.consistencycheck;
          break;
        case constraintblock: bfstream >> gin.constraint;
          break;
        case covalentformblock: bfstream >> gin.covalentform;
          break;
        case debugblock: bfstream >> gin.debug;
          break;
        case dihedralresblock: bfstream >> gin.dihedralres;
          break;
        case angleresblock: bfstream >> gin.angleres;
          break;
        case distancefieldblock: bfstream >> gin.distancefield;
          break;
        case distanceresblock: bfstream >> gin.distanceres;
          break;
        case energyminblock: bfstream >> gin.energymin;
          break;
        case edsblock: bfstream >> gin.eds;
          break;
        case ewarnblock: bfstream >> gin.ewarn;
          break;
        case forceblock: bfstream >> gin.force;
          break;
        case gamdblock: bfstream >> gin.gamd;
          break;
        case geomconstraintsblock: bfstream >> gin.geomconstraints;
          break;
        case gromos96compatblock: bfstream >> gin.gromos96compat;
          break;
        case initialiseblock: bfstream >> gin.initialise;
          break;
        case innerloopblock: bfstream >> gin.innerloop;
          break;
        case integrateblock: bfstream >> gin.integrate;
          break;
        case jvalueresblock: bfstream >> gin.jvalueres;
          break;
        case lambdasblock: bfstream >> gin.lambdas;
          break;
        case localelevblock: bfstream >> gin.localelev;
          break;
        case electricblock: bfstream >> gin.electric;
          break;
        case multibathblock: bfstream >> gin.multibath;
          break;
        case multicellblock: bfstream >> gin.multicell;
          break;
        case multigradientblock : bfstream >> gin.multigradient;
          break;
        case multistepblock: bfstream >> gin.multistep;
          break;
        case neighbourlistblock: bfstream >> gin.neighbourlist;
          break;
        case nemdblock: bfstream >> gin.nemd;
          break;
        case nonbondedblock: bfstream >> gin.nonbonded;
          break;
        case orderparamresblock: bfstream >> gin.orderparamres;
          break;
        case overalltransrotblock: bfstream >> gin.overalltransrot;
          break;
        case pairlistblock: bfstream >> gin.pairlist;
          break;
        case pathintblock: bfstream >> gin.pathint;
          break;
        case perscaleblock: bfstream >> gin.perscale;
          break;
        case perturbationblock: bfstream >> gin.perturbation;
          break;
        case polariseblock: bfstream >> gin.polarise;
          break;
        case positionresblock: bfstream >> gin.positionres;
          break;
        case precalclamblock: bfstream >> gin.precalclam;
          break;
        case pressurescaleblock: bfstream >> gin.pressurescale;
          break;
        case printoutblock: bfstream >> gin.printout;
          break;
        case qmmmblock: bfstream >> gin.qmmm;
          break;
        case randomnumbersblock: bfstream >> gin.randomnumbers;
          break;
        case readtrajblock: bfstream >> gin.readtraj;
          break;
        case replicablock: bfstream >> gin.replica;
          break;
        case reedsblock: bfstream >> gin.reeds;
          break;
        case rottransblock: bfstream >> gin.rottrans;
          break;
        case sasablock: bfstream >> gin.sasa;
          break;
        case stepblock: bfstream >> gin.step;
          break;
        case stochdynblock: bfstream >> gin.stochdyn;
          break;
        case symresblock: bfstream >> gin.symres;
          break;
        case systemblock: bfstream >> gin.system;
          break;
        case thermostatblock: bfstream >> gin.thermostat;
          break;
        case umbrellablock: bfstream >> gin.umbrella;
          break;
        case virialblock: bfstream >> gin.virial;
          break;
        case virtualatomblock: bfstream >> gin.virtualatom;
          break;
        case writetrajblock: bfstream >> gin.writetraj;
          break;
        case xrayresblock: bfstream >> gin.xrayres;
          break;
        case unknown:
          iunknown newblock(buffer[0]);
          bfstream >> newblock;
          gin.unknown.push_back(newblock);
          std::stringstream msg;
          msg << "Don't know anything about block " << buffer[0]
                  << ". Just storing data.";
          printWarning(msg.str());
      }
    }
  }
  return is;
}

gio::Ginstream & operator>>(gio::Ginstream &is, fileInfo &s) {

  std::string e;
  std::string first;
  std::vector<std::string> buffer;
  is.getline(first);

  while (!is.stream().eof()) {
    is.getblock(buffer);
    s.blocks.push_back(first);
    s.blockslength.push_back(buffer.size() - 1);
    is.getline(first);
  }
  return is;
}

// TEMPLATE handling of (output) filenames

class filename {
  std::vector<std::string> d_parts;
  double d_time, d_dt;
  int d_start;
  std::string d_system;
  std::string d_queue;
  std::string d_template;

public:
  filename();

  filename(std::string s, double t, double dt, int start = 1, std::string q = "") {
    d_system = s;
    d_time = t;
    d_dt = dt;
    d_start = start;
    d_queue = q;
  };

  void setInfo(std::string s, double t, double dt, int start = 1, std::string q = "") {
    d_system = s;
    d_time = t;
    d_dt = dt;
    d_start = start;
    d_queue = q;
  };

  void setTemplate(std::string s);

  std::string temp() {
    return d_template;
  };
  std::string name(int number);
};

void filename::setTemplate(std::string s) {
  d_template = s;
  d_parts.clear();
  std::string::size_type iter;

  std::string sub;
  iter = s.find('%');
  while (iter != std::string::npos) {
    sub = s.substr(0, iter);
    s = s.substr(iter + 1, s.size() - iter - 1);
    iter = s.find('%');
    d_parts.push_back(sub);

  }
  d_parts.push_back(s);
}

std::string filename::name(int number) {
  std::ostringstream os;
  for (unsigned int i = 0; i < d_parts.size(); i++) {
    if (i % 2) {
      switch (TEMPLATE[d_parts[i]]) {
        case systemtemplate: os << d_system;
          break;
        case numbertemplate: os << d_start + number;
          break;
        case oldnumbertemplate: os << d_start + number - 1;
          break;
        case start_timetemplate: os << d_time + number*d_dt;
          break;
        case end_timetemplate: os << d_time + (number + 1) * d_dt;
          break;
        case queuetemplate: os << d_queue;
          break;
        case unknowntemplate:
          std::cout << "Do not know how to handle " << d_parts[i]
                  << " in template. Just printing the words." << std::endl;
          os << d_parts[i];
          break;
      }
    } else os << d_parts[i];
  }
  return os.str();
}

// Jobinfo

class jobinfo {
public:
  std::map<std::string, std::string> param;
  std::string dir;
  int prev_id;
};

// Job submission directive

class directive : public filename {
public:
  directive(std::string system, double t, double dt, int start = 1, std::string q = "") :
      filename(system, t, dt, start, q) {}
};

// Writing out of an input file

std::ostream & operator<<(std::ostream &os, input &gin) {
  // MOLECULAR SYSTEM

  // SYSTEM (g96, promd, md++)
  if (gin.system.found)
    os << "SYSTEM\n"
          << "#      NPM      NSM\n"
          << std::setw(10) << gin.system.npm
          << std::setw(9) << gin.system.nsm
          << "\nEND\n";

  // METHOD EMPLOYED

  // CONSISTENCYCHECK (promd)
  if (gin.consistencycheck.found) {
    os << "CONSISTENCYCHECK\n"
            << "#    NTCHK     NTCKF     FDCKF     NTCKV     FDCKV\n"
            << std::setw(10) << gin.consistencycheck.ntchk
            << std::setw(10) << gin.consistencycheck.ntckf
            << std::setw(10) << gin.consistencycheck.fdckf
            << std::setw(10) << gin.consistencycheck.ntckv
            << std::setw(10) << gin.consistencycheck.fdckv
            << "\n"
            << "#    NTCKT     NTCKE     NTCKR\n"
            << std::setw(10) << gin.consistencycheck.ntckt
            << std::setw(10) << gin.consistencycheck.ntcke
            << std::setw(10) << gin.consistencycheck.ntckr
            << "\n"
            << "#    NTCKL     FDCKL\n"
            << std::setw(10) << gin.consistencycheck.ntckl
            << std::setw(10) << gin.consistencycheck.fdckl
            << "\n"
            << "#    NACKF\n"
            << std::setw(10) << gin.consistencycheck.nckf.size()
            << "\n"
            << "# NCKF(1...NACKF)\n";
    for (unsigned int i = 0; i < gin.consistencycheck.nckf.size(); ++i) {
      os << std::setw(10) << gin.consistencycheck.nckf[i];
    }
    os << "\nEND\n";
  }

  // ENERGYMIN (promd, md++): only write if NTEM != 0
  if (gin.energymin.found && gin.energymin.ntem) {
    os << "ENERGYMIN\n"
            << "#    NTEM    NCYC    DELE    DX0    DXM   NMIN   FLIM\n"
            << std::setw(10) << gin.energymin.ntem
            << std::setw(10) << gin.energymin.ncyc
            << std::setw(10) << gin.energymin.dele
            << std::setw(10) << gin.energymin.dx0
            << std::setw(10) << gin.energymin.dxm
            << std::setw(10) << gin.energymin.nmin
            << std::setw(10) << gin.energymin.flim;
            if (gin.energymin.ntem > 1) {
              os  << std::setw(10) << gin.energymin.cgim
                  << std::setw(10) << gin.energymin.cgic;
            }
            os << "\nEND\n";
  }
  // STOCHDYN (promd, md++)
  if (gin.stochdyn.found) {
    os << "STOCHDYN\n"
            << "#     NTSD      NTFR      NSFR     NBREF     RCUTF     CFRIC    TEMPSD\n"
            << std::setw(10) << gin.stochdyn.ntsd
            << std::setw(10) << gin.stochdyn.ntfr
            << std::setw(10) << gin.stochdyn.nsfr
            << std::setw(10) << gin.stochdyn.nbref
            << std::setw(10) << gin.stochdyn.rcutf
            << std::setw(10) << gin.stochdyn.cfric
            << std::setw(10) << gin.stochdyn.tempsd
            << "\nEND\n";
  }
  // READTRAJ (promd, md++)
  if (gin.readtraj.found) {
    os << "READTRAJ\n"
            << "#     NTRD      NTSTR     NTRB     NTSHK\n"
            << std::setw(10) << gin.readtraj.ntrd
            << std::setw(10) << gin.readtraj.ntstr
            << std::setw(10) << gin.readtraj.ntrb
            << std::setw(10) << gin.readtraj.ntshk
            << "\nEND\n";
  }
  // STEP (promd, md++, g96)
  if (gin.step.found) {
    os << "STEP\n"
            << "#   NSTLIM         T        DT\n"
            << std::setw(10) << gin.step.nstlim
            << std::setw(10) << gin.step.t
            << std::setw(10) << gin.step.dt
            << "\nEND\n";
  }
  // REPLICA (md++)
  if (gin.replica.found) {
    os << "REPLICA\n"
            << "#     RETL      NRET\n"
            << std::setw(10) << gin.replica.nrel
            << std::setw(10) << gin.replica.ret.size()
            << "\n#  RET(1 ... NRET)\n";
    for (unsigned int i = 0; i < gin.replica.ret.size(); ++i) {
      os << std::setw(10) << gin.replica.ret[i];
    }
    os << "\n# LRESCALE\n"
            << gin.replica.lrescale
            << "\n#   NRELAM\n"
            << gin.replica.relam.size()
            << "\n#  RELAM(1 ... NRELAM)\n";
    for (unsigned int i = 0; i < gin.replica.relam.size(); ++i) {
      os << std::setw(10) << gin.replica.relam[i];
    }
    os << "\n#   RETS(1 ... NRELAM)\n";
    for (unsigned int i = 0; i < gin.replica.rets.size(); ++i) {
      os << std::setw(10) << gin.replica.rets[i];
    }
    os << "\n# NRETRIAL   NREQUIL    CONT\n"
            << std::setw(10) << gin.replica.nretrial
            << std::setw(10) << gin.replica.nrequil
            << std::setw(10) << gin.replica.cont
            << "\nEND\n";
  }

  // SPACIAL BOUNDARY CONDITIONS

  // BOUNDCOND (promd, md++)
  if (gin.boundcond.found) {
    os << "BOUNDCOND\n"
            << "#      NTB    NDFMIN\n"
            << std::setw(10) << gin.boundcond.ntb
            << std::setw(10) << gin.boundcond.ndfmin
            << "\nEND\n";
  }
  // MULTICELL (promd, md++)
  if (gin.multicell.found) {
    os << "MULTICELL\n"
            << "#      NTM    NCELLA    NCELLB    NCELLC\n"
            << std::setw(10) << gin.multicell.ntm
            << std::setw(10) << gin.multicell.ncella
            << std::setw(10) << gin.multicell.ncellb
            << std::setw(10) << gin.multicell.ncellc
            << "\n"
            << "#     TOLPX    TOLPV     TOLPF    TOLPFW\n"
            << std::setw(10) << gin.multicell.tolpx
            << std::setw(10) << gin.multicell.tolpv
            << std::setw(10) << gin.multicell.tolpf
            << std::setw(10) << gin.multicell.tolpfw
            << "\nEND\n";
  }

  // MULTIGRADIENT (md++)
  if (gin.multigradient.found) {
    os << "MULTIGRADIENT\n"
            << "#   NTMGRE     NTMGRP     NTMGRN\n"
            << std::setw(10) << gin.multigradient.ntmgre
            << std::setw(10) << gin.multigradient.ntmgrp
            << std::setw(10) << gin.multigradient.ntmgrn
            << "\n";
    for(int i = 0; i < gin.multigradient.ntmgrn; ++i) {
      os << "#   MGRVAR    MGRFRM    MGRNCP\n"
              << std::setw(10) << gin.multigradient.mgrvar[i]
              << std::setw(10) << gin.multigradient.mgrfrm[i]
              << std::setw(10) << gin.multigradient.curves[i].size() << "\n"
              << "#        MGRCPT         MGRCPV\n";
      for(unsigned int j = 0; j < gin.multigradient.curves[i].size(); ++j) {
        os << std::setw(15) << gin.multigradient.curves[i][j].first
                << std::setw(15) << gin.multigradient.curves[i][j].second << "\n";
      }
    }
    os << "END\n";
  }

  // THERMODYNAMIC BOUNDARY CONDITIONS
  // THERMOSTAT (promd)
  if (gin.thermostat.found) {
    os << "THERMOSTAT\n"
            << "#       NTT     NTBTH     NTSET\n"
            << std::setw(11) << gin.thermostat.ntt
            << std::setw(10) << gin.thermostat.ntbth
            << std::setw(10) << gin.thermostat.ntset
            << "\n"
            << "# I = 1 ... NTBTH\n"
            << "#         I NTBTYP(I) TEMBTH(I) NTBVAR(I) TAUBTH(I,1...NTVAR)\n";
    for (unsigned int i = 0; i < gin.thermostat.baths.size(); ++i) {
      os << std::setw(11) << gin.thermostat.baths[i].index
              << std::setw(10) << gin.thermostat.baths[i].ntbtyp
              << std::setw(10) << gin.thermostat.baths[i].tembth
              << std::setw(10) << gin.thermostat.baths[i].ntbvar;
      for (int j = 0; j < gin.thermostat.baths[i].ntbvar; j++) {
        os << std::setw(10) << gin.thermostat.baths[i].taubth[j];
      }
      os << "\n";
    }
    os << "# NTSCPL(J) NTSTYP(J) NTSCNS(J)  NTSGT(J) NTSGTG(J,1...NTGT(J))\n";
    for (unsigned int j = 0; j < gin.thermostat.dofgroups.size(); ++j) {
      os << std::setw(11) << gin.thermostat.dofgroups[j].ntscpl
              << std::setw(10) << gin.thermostat.dofgroups[j].ntstyp
              << std::setw(10) << gin.thermostat.dofgroups[j].ntscns
              << std::setw(10) << gin.thermostat.dofgroups[j].ntsgt;
      int upper_k = 0;
      if (gin.thermostat.dofgroups[j].ntsgt == -2) {
        upper_k = 2;
      } else if (gin.thermostat.dofgroups[j].ntsgt == -1 || gin.thermostat.dofgroups[j].ntsgt == 0) {
        upper_k = 0;
      //} else if (gin.thermostat.dofgroups[j].ntsgt == 0) {
      //  stringstream msg;
      //  msg << "NTSGT(" << j + 1 << ")";
      //  printIO("THERMOSTAT", msg.str(), "0", "not implemented");
      } else { // NTSGT > 0
        upper_k = gin.thermostat.dofgroups[j].ntsgt;
      }
      for (int k = 0; k < upper_k; k++) {
        os << std::setw(10) << gin.thermostat.dofgroups[j].ntsgtg[k];
      }
      os << "\n";
    }
    os << "END\n";
  }
  // MULTIBATH (md++)
  if (gin.multibath.found) {
    os << "MULTIBATH\n"
            << "# NTBTYP:\n"
            << "#      weak-coupling:      use weak-coupling scheme\n"
            << "#      nose-hoover:        use Nose Hoover scheme\n"
            << "#      nose-hoover-chains: use Nose Hoover chains scheme\n"
            << "# NUM: number of chains in Nose Hoover chains scheme\n"
            << "#      !! only specify NUM when needed !!\n"
            << "# NBATHS: number of temperature baths to couple to\n";
    if (gin.multibath.ntbtyp == 2) {
      os << "#          NTBTYP     NUM\n"
              << std::setw(20) << gin.multibath.ntbtyp
              << std::setw(8) << gin.multibath.num;
    } else {
      os << "#          NTBTYP\n"
              << std::setw(20) << gin.multibath.ntbtyp;
    }
    os << "\n#  NBATHS\n"
            << std::setw(10) << gin.multibath.nbaths
            << "\n";
    os << "# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)\n";
    for (int i = 0; i < gin.multibath.nbaths; ++i) {
      os << std::setw(10) << gin.multibath.temp0[i]
              << std::setw(10) << gin.multibath.tau[i] << std::endl;
    }
    os << "\n";
    os << "#   DOFSET: number of distinguishable sets of d.o.f.\n";
    os << std::setw(10) << gin.multibath.dofset << "\n";
    os << "# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)\n";
    for (int i = 0; i < gin.multibath.dofset; ++i) {
      os << std::setw(10) << gin.multibath.last[i]
              << std::setw(10) << gin.multibath.combath[i]
              << std::setw(10) << gin.multibath.irbath[i];
    }
    os << "\nEND\n";
  }
  // BAROSTAT (promd)
  if (gin.barostat.found) {
    os << "BAROSTAT\n"
            << "#      NTP      NPVAR     NPBTH      COMP\n"
            << std::setw(10) << gin.barostat.ntp
            << std::setw(10) << gin.barostat.npvar
            << std::setw(10) << gin.barostat.pbaths.size()
            << std::setw(10) << gin.barostat.comp
            << "\n"
            << "# I = 1 ... NPBTH\n"
            << "# I  PRSBTH(I)  TAUBBA(I,1...NPVAR)\n";
    for (unsigned int i = 0; i < gin.barostat.pbaths.size(); ++i) {
      os << std::setw(3) << i + 1
              << std::setw(11) << gin.barostat.pbaths[i].prsbth;
      for (unsigned int j = 0; j < gin.barostat.pbaths[i].taubba.size(); ++j) {
        os << std::setw(7) << gin.barostat.pbaths[i].taubba[j];
      }
      os << "\n";
    }
    os << "#  NPCPL(1...6)\n"
            << std::setw(6) << gin.barostat.npcpl[0]
            << std::setw(6) << gin.barostat.npcpl[1]
            << std::setw(6) << gin.barostat.npcpl[2]
            << std::setw(6) << gin.barostat.npcpl[3]
            << std::setw(6) << gin.barostat.npcpl[4]
            << std::setw(6) << gin.barostat.npcpl[5]
            << "\nEND\n";
  }
  // VIRIAL (promd)
  if (gin.virial.found) {
    os << "VIRIAL\n"
            << "#      NTV       NTVG\n"
            << std::setw(10) << gin.virial.ntv
            << std::setw(10) << gin.virial.ntvg
            << "\nEND\n";
  }
  // PRESSURESCALE
  if (gin.pressurescale.found) {
    os << "PRESSURESCALE\n"
            << "# COUPLE   SCALE    COMP    TAUP  VIRIAL\n"
            << std::setw(8) << gin.pressurescale.couple
            << std::setw(8) << gin.pressurescale.scale << " "
            << std::setw(8) << gin.pressurescale.comp << " "
            << std::setw(8) << gin.pressurescale.taup << " "
            << std::setw(8) << gin.pressurescale.virial
            << "\n# SEMIANISOTROPIC COUPLINGS(X, Y, Z)\n"
            << std::setw(8) << gin.pressurescale.x_semi << " "
            << std::setw(8) << gin.pressurescale.y_semi << " "
            << std::setw(8) << gin.pressurescale.z_semi
            << "\n# PRES0(1...3,1...3)\n"
            << std::setw(8) << gin.pressurescale.pres0[0][0]
            << std::setw(8) << gin.pressurescale.pres0[0][1]
            << std::setw(8) << gin.pressurescale.pres0[0][2]
            << "\n"
            << std::setw(8) << gin.pressurescale.pres0[1][0]
            << std::setw(8) << gin.pressurescale.pres0[1][1]
            << std::setw(8) << gin.pressurescale.pres0[1][2]
            << "\n"
            << std::setw(8) << gin.pressurescale.pres0[2][0]
            << std::setw(8) << gin.pressurescale.pres0[2][1]
            << std::setw(8) << gin.pressurescale.pres0[2][2]
            << "\nEND\n";
  }

  // INTERACTION EVALUATION

  // FORCE (promd, md++, g96)
  if (gin.force.found) {
    int size = gin.force.nre.size();
    if (gin.force.force_groups)
      size = -size;
    os << "FORCE\n"
            << "#      NTF array\n"
            << "# bonds    angles   imp.     dihe     charge nonbonded\n"
            << std::setw(3) << gin.force.ntf[0] << std::setw(9) << gin.force.ntf[1]
            << std::setw(9) << gin.force.ntf[2] << std::setw(9) << gin.force.ntf[3]
            << std::setw(9) << gin.force.ntf[4] << std::setw(9) << gin.force.ntf[5]
            << "\n# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)\n"
            << std::setw(6) << size << "\n";
    int countnre = 0;
    for (unsigned int i = 0; i < gin.force.nre.size(); i++) {
      os << std::setw(9) << gin.force.nre[i];
      countnre++;
      if (countnre % 8 == 0) os << std::endl;
    }
    os << "\nEND\n";
  }
  // MULTISTEP
  if (gin.multistep.found) {
    os << "MULTISTEP\n"
            << "# STEPS calculate non-bonded every STEPSth step.\n"
            << "# BOOST 0,1\n"
            << "#       0: stored forces of STEPSth step are added every step\n"
            << "#       1: stored forces of STEPSth setp are multiplied by STEPS\n"
            << "#       and added every STEPSth step.\n"
            << "#\n"
            << "#" << std::setw(14) << "STEPS" << std::setw(15) << "BOOST\n"
            << std::setw(15) << gin.multistep.steps
            << std::setw(15) << gin.multistep.boost << std::endl
            << "END\n";
  }
  // COVALENTFORM (promd, md++)
  if (gin.covalentform.found) {
    os << "COVALENTFORM\n"
            << "#    NTBBH    NTBAH     NTBDN\n"
            << std::setw(10) << gin.covalentform.ntbbh
            << std::setw(10) << gin.covalentform.ntbah
            << std::setw(10) << gin.covalentform.ntbdn
            << "\nEND\n";
  }
  // VIRTUALATOM (md++)
  if (gin.virtualatom.found) {
    os << "VIRTUALATOM\n"
       << "#     VIRT   NUMVIRT  LASTVIRT\n"
       << std::setw(10) << gin.virtualatom.virt
       << std::setw(10) << gin.virtualatom.numvirt
       << std::setw(10) << gin.virtualatom.lastvirt
       << "\nEND\n";
  }
  // GEOMCONSTRAINTS (promd)
  if (gin.geomconstraints.found) {
    os << "GEOMCONSTRAINTS\n"
            << "#    NTCPH     NTCPN      NTCS    SHKTOL\n"
            << std::setw(10) << gin.geomconstraints.ntcph
            << std::setw(10) << gin.geomconstraints.ntcpn
            << std::setw(10) << gin.geomconstraints.ntcs
            << std::setw(10) << gin.geomconstraints.shktol
            << "\nEND\n";
  }
  // CONSTRAINT (md++)
  if (gin.constraint.found) {
    os << "CONSTRAINT\n"
            << "# NTC\n"
            << std::setw(5) << gin.constraint.ntc
            << "\n";
    if (gin.constraint.ntcp == 3) {
      os << "#      NTCP  NTCP0(1 ... 3)\n"
              << std::setw(11) << gin.constraint.ntcp
              << std::setw(10) << gin.constraint.ntcp0[0]
              << std::setw(10) << gin.constraint.ntcp0[1]
              << std::setw(10) << gin.constraint.ntcp0[2];
    } else {
      os << "#      NTCP  NTCP0(1)\n"
              << std::setw(11) << gin.constraint.ntcp
              << std::setw(10) << gin.constraint.ntcp0[0];
    }
    os << "\n";
    if (gin.constraint.ntcs == 1 || gin.constraint.ntcs == 2 || gin.constraint.ntcs == 5) {
      os << "#      NTCS  NTCS0(1)\n"
              << std::setw(11) << gin.constraint.ntcs
              << std::setw(10) << gin.constraint.ntcs0[0];
    } else if (gin.constraint.ntcs == 3) {
      os << "#      NTCS  NTCS0(1 ... 3)\n"
              << std::setw(11) << gin.constraint.ntcs
              << std::setw(10) << gin.constraint.ntcs0[0]
              << std::setw(10) << gin.constraint.ntcs0[1]
              << std::setw(10) << gin.constraint.ntcs0[2];
    } else if (gin.constraint.ntcs == 6) {
      os << "#      NTCS  NTCS0(1)      NTCG      NTCD\n"
              << std::setw(11) << gin.constraint.ntcs
              << std::setw(10) << gin.constraint.ntcs0[0]
              << std::setw(10) << gin.constraint.ntcg;
          for (int g = 0; g < gin.constraint.ntcg; g++) {
              os << std::setw(10) << gin.constraint.ntcd[g];
          }
    } else {
      os << "#      NTCS\n"
              << std::setw(11) << gin.constraint.ntcs;
    }
    os << "\nEND\n";
  }
  // GROMOS96COMPAT (promd)
  if (gin.gromos96compat.found) {
    os << "GROMOS96COMPAT\n"
            << "#   NTNB96    NTR96     NTP96     NTG96\n"
            << std::setw(10) << gin.gromos96compat.ntnb96
            << std::setw(10) << gin.gromos96compat.ntr96
            << std::setw(10) << gin.gromos96compat.ntp96
            << std::setw(10) << gin.gromos96compat.ntg96
            << "\nEND\n";
  }
  // SASA (md++)
  if (gin.sasa.found) {
    os << "SASA\n"
            << "#   NTSASA     NTVOL      P_12      P_13      P_1X    SIGMAV     RSOlV       AS1       AS2\n"
            << std::setw(10) << gin.sasa.ntsasa
            << std::setw(10) << gin.sasa.ntvol
            << std::setw(10) << gin.sasa.p12
            << std::setw(10) << gin.sasa.p13
            << std::setw(10) << gin.sasa.p1x
            << std::setw(10) << gin.sasa.sigmav
            << std::setw(10) << gin.sasa.rsolv
            << std::setw(10) << gin.sasa.as1
            << std::setw(10) << gin.sasa.as2
            << "\nEND\n";
  }
  // PATHINT (promd)
  if (gin.pathint.found) {
    os << "PATHINT\n"
            << "#  NTPI\n"
            << std::setw(7) << gin.pathint.ntpi
            << "\nEND\n";
  }
  // POLARISE 
  if (gin.polarise.found) {
    os << "POLARISE\n"
            << "#       COS    EFIELD       MINFIELD      DAMP     WRITE\n"
            << std::setw(11) << gin.polarise.cos
            << std::setw(10) << gin.polarise.efield
            << std::setw(15) << gin.polarise.minfield
            << std::setw(10) << gin.polarise.damp
            << std::setw(10) << gin.polarise.write
            << "\nEND\n";
  }
  // INTEGRATE (md++)
  if (gin.integrate.found) {
    os << "INTEGRATE\n"
            << "#   NINT\n"
            << std::setw(8) << gin.integrate.nint
            << "\nEND\n";
  }
  // CGRAIN (md++)
  if (gin.cgrain.found) {
    os << "CGRAIN\n"
            << "#  NTCGRAN       EPS\n"
            << std::setw(10) << gin.cgrain.ntcgran
            << std::setw(10) << gin.cgrain.eps
            << std::setw(10) << gin.cgrain.epsm
            << "\nEND\n";
  }
  // ROTTRANS (md++)
  if (gin.rottrans.found) {
    os << "ROTTRANS\n"
            << "#      RTC   RTCLAST\n"
            << std::setw(10) << gin.rottrans.rtc
            << std::setw(10) << gin.rottrans.rtclast
            << "\nEND\n";
  }
  // INNERLOOP (md++)
  if (gin.innerloop.found) {
    os << "INNERLOOP\n"
            << "#     NTILM      NTILS      NGPUS      NDEVG\n"
            << std::setw(10) << gin.innerloop.ntilm
            << std::setw(10) << gin.innerloop.ntils;
        if(gin.innerloop.ntilm==4){
          os << std::setw(10) << gin.innerloop.ngpus;
          for (int g = 0; g < gin.innerloop.ngpus; g++) {
            os << std::setw(10) << gin.innerloop.ndevg[g];
          }
	}
        os << "\nEND\n";
  }

  // PAIRLIST GENERATION

  // NEIGHBOURLIST (promd)
  if (gin.neighbourlist.found) {
    os << "NEIGHBOURLIST\n"
            << "#      ALGO  NUPDPL  NUPDIS  NUPDII\n"
            << std::setw(11) << gin.neighbourlist.plalgo
            << std::setw(8) << gin.neighbourlist.nupdpl
            << std::setw(8) << gin.neighbourlist.nupdis
            << std::setw(8) << gin.neighbourlist.nupdii << std::endl
            << "#     RCUTS   RCUTI  GRDSZX  GRDSZY  GRDSZZ\n"
            << std::setw(11) << gin.neighbourlist.rcuts
            << std::setw(8) << gin.neighbourlist.rcuti
            << std::setw(8) << gin.neighbourlist.gridszx
            << std::setw(8) << gin.neighbourlist.gridszy
            << std::setw(8) << gin.neighbourlist.gridszz << std::endl
            << "#      TYPE  NCGCEN\n"
            << std::setw(11) << gin.neighbourlist.type
            << std::setw(8) << gin.neighbourlist.ncgcen
            << "\nEND\n";
  }
  // PAIRLIST (md++)
  if (gin.pairlist.found) {
    os << "PAIRLIST\n"
            << "# algorithm    NSNB   RCUTP   RCUTL    SIZE    TYPE\n"
            << std::setw(11) << gin.pairlist.algorithm
            << std::setw(8) << gin.pairlist.nsnb
            << std::setw(8) << gin.pairlist.rcutp
            << std::setw(8) << gin.pairlist.rcutl
            << std::setw(8) << gin.pairlist.size
            << std::setw(8) << gin.pairlist.type
            << "\nEND\n";
  }

  // LONGRANGE INTERACTIONS

  // NONBONDED (promd)
  if (gin.nonbonded.found) {
    os << "NONBONDED\n"
            << "# NLRELE\n"
            << std::setw(10) << gin.nonbonded.nlrele
            << "\n"
            << "#  APPAK    RCRF   EPSRF    NSLFEXCL\n"
            << std::setw(10) << gin.nonbonded.appak
            << std::setw(10) << gin.nonbonded.rcrf
            << std::setw(10) << gin.nonbonded.epsrf
            << std::setw(10) << gin.nonbonded.nslfexcl
            << "\n"
            << "# NSHAPE  ASHAPE  NA2CLC   TOLA2   EPSLS\n"
            << std::setw(10) << gin.nonbonded.nshape
            << std::setw(10) << gin.nonbonded.ashape
            << std::setw(10) << gin.nonbonded.na2clc
            << std::setw(10) << gin.nonbonded.tola2
            << std::setw(10) << gin.nonbonded.epsls
            << "\n"
            << "#    NKX     NKY     NKZ   KCUT\n"
            << std::setw(10) << gin.nonbonded.nkx
            << std::setw(10) << gin.nonbonded.nky
            << std::setw(10) << gin.nonbonded.nkz
            << std::setw(10) << gin.nonbonded.kcut
            << "\n"
            << "#    NGX     NGY     NGZ  NASORD  NFDORD  NALIAS  NSPORD\n"
            << std::setw(10) << gin.nonbonded.ngx
            << std::setw(10) << gin.nonbonded.ngy
            << std::setw(10) << gin.nonbonded.ngz
            << std::setw(10) << gin.nonbonded.nasord
            << std::setw(10) << gin.nonbonded.nfdord
            << std::setw(10) << gin.nonbonded.nalias
            << std::setw(10) << gin.nonbonded.nspord
            << "\n"
            << "# NQEVAL  FACCUR  NRDGRD  NWRGRD\n"
            << std::setw(10) << gin.nonbonded.nqeval
            << std::setw(10) << gin.nonbonded.faccur
            << std::setw(10) << gin.nonbonded.nrdgrd
            << std::setw(10) << gin.nonbonded.nwrgrd
            << "\n"
            << "#  NLRLJ  SLVDNS\n"
            << std::setw(10) << gin.nonbonded.nlrlj
            << std::setw(10) << gin.nonbonded.slvdns
            << "\nEND\n";
  }

  // INITIALISATION OF THE RUN

  // INITIALISE (promd, md++)
  if (gin.initialise.found) {
    os << "INITIALISE\n"
            << "# Default values for NTI values: 0\n"
            << "#   NTIVEL    NTISHK    NTINHT    NTINHB\n"
            << std::setw(10) << gin.initialise.ntivel
            << std::setw(10) << gin.initialise.ntishk
            << std::setw(10) << gin.initialise.ntinht
            << std::setw(10) << gin.initialise.ntinhb
            << "\n"
            << "#   NTISHI    NTIRTC    NTICOM\n"
            << std::setw(10) << gin.initialise.ntishi
            << std::setw(10) << gin.initialise.ntirtc
            << std::setw(10) << gin.initialise.nticom
            << "\n"
            << "#   NTISTI\n"
            << std::setw(10) << gin.initialise.ntisti
            << "\n"
            << "#       IG     TEMPI\n"
            << std::setw(10) << gin.initialise.ig
            << std::setw(10) << gin.initialise.tempi
            << "\nEND\n";
  }
  // RANDOMNUMBERS (md++)
  if (gin.randomnumbers.found) {
    os << "RANDOMNUMBERS\n"
            << "#   NTRNG   NTGSL\n"
            << std::setw(8) << gin.randomnumbers.ntrng
            << std::setw(8) << gin.randomnumbers.ntgsl
            << "\nEND\n";
  }

  // CENTRE-OF-MASS MOTION

  // OVERALLTRANSROT (promd)
  if (gin.overalltransrot.found) {
    os << "OVERALLTRANSROT\n"
            << "#    NCMTR     NCMRO     CMAMX     CMAMY     CMAMZ\n"
            << std::setw(10) << gin.overalltransrot.ncmtr
            << std::setw(10) << gin.overalltransrot.ncmro
            << std::setw(10) << gin.overalltransrot.cmamx
            << std::setw(10) << gin.overalltransrot.cmamy
            << std::setw(10) << gin.overalltransrot.cmamz
            << "\nEND\n";
  };
  // COMTRANSROT (md++)
  if (gin.comtransrot.found) {
    os << "COMTRANSROT\n"
            << "#     NSCM\n"
            << std::setw(10) << gin.comtransrot.nscm
            << "\nEND\n";
  }

  // SPECIAL FORCES

  // POSITIONRES (promd, md++)
  if (gin.positionres.found) {
    os << "POSITIONRES\n"
            << "#     values for NTPOR\n"
            << "#     0: no position re(con)straining\n"
            << "#     1: use CPOR\n"
            << "#     2: use CPOR/ ATOMIC B-FACTORS\n"
            << "#     3: position constraining\n"
            << "#    NTPOR    NTPORB  NTPORS      CPOR\n"
            << std::setw(10) << gin.positionres.ntpor
            << std::setw(10) << gin.positionres.ntporb
            << std::setw(10) << gin.positionres.ntpors
            << std::setw(10) << gin.positionres.cpor
            << "\nEND\n";
  }
  // XRAYRES (promd, md++)
  if (gin.xrayres.found) {
    os << "XRAYRES\n"
            << "#    NTXR   NTXLE   CXR   NTWXR   NTWDE   NTWXM   CXTAU  RDAVG\n"
            << std::setw(10) << gin.xrayres.ntxr
            << std::setw(10) << gin.xrayres.ntxle
            << std::setw(10) << gin.xrayres.cxr
            << std::setw(10) << gin.xrayres.ntwxr
            << std::setw(10) << gin.xrayres.ntwde
            << std::setw(10) << gin.xrayres.ntwxm
            << std::setw(10) << gin.xrayres.cxtau
            << std::setw(10) << gin.xrayres.rdavg
            << "\nEND\n";
  }
  // COLVARRES (md++)
  if (gin.colvarres.found) {
    os << "COLVARRES\n"
            << "#      CVR       CVK    TAUCVR      VCVR     NTWCV\n"
            << std::setw(10) << gin.colvarres.cvr
            << std::setw(10) << gin.colvarres.cvk
            << std::setw(10) << gin.colvarres.taucvr
            << std::setw(10) << gin.colvarres.vcvr
            << std::setw(10) << gin.colvarres.ntwcv
            << "\nEND\n";
  }
  // DISTANCERES (promd, md++)
  if (gin.distanceres.found) {
    os << "DISTANCERES\n"
            << "# NTDIR\n"
            << "#   0 : no distance restraining\n"
            << "#   -1,1 : use CDIS\n"
            << "#   -2,2: use W0*CDIS\n"
            << "#   NTDIR < 0 : time averaging\n"
            << "#   NTDIR > 0 : no time averaging\n"
            << "# NTDIRA = 1: read in time averaged distances (for continuation run)\n"
            << "# NTDIRA = 0: don't read them in, recalc from scratch\n"
            << "# NTWDIR >= 0 write every NTWDIRth step dist. restr. information to external file\n"
            << "#     NTDIR  NTDIRA    CDIR    DIR0  TAUDIR  FORCESCALE VDIR NTWDIR\n"
            << std::setw(11) << gin.distanceres.ntdir
            << std::setw(8) << gin.distanceres.ntdira
            << std::setw(8) << gin.distanceres.cdir
            << std::setw(8) << gin.distanceres.dir0
            << std::setw(8) << gin.distanceres.taudir
            << std::setw(8) << gin.distanceres.forcescale
            << std::setw(8) << gin.distanceres.vdir
            << std::setw(8) << gin.distanceres.ntwdir
            << "\nEND\n";
  }
  // DISTANCEFIELD (md++)
  if (gin.distancefield.found){
    os << "DISTANCEFIELD\n"
       << "# NTDFR 0,1        controls distance field restraining\n"
       << "#       0: no distance field restraining\n"
       << "#       1: apply distance field restraining\n"
       << "# GRID > 0.0        grid size for distance field\n"
       << "# PROTEINOFFSET > 0 penalty for distances through the host\n"
       << "# PROTEINCUTOFF > 0 distance to protein atoms to be considered inside\n"
       << "# PROTECT >= 0      protect grid points within this radius around the zero-distance\n"
       << "#                   point from being flagged as protein\n"
       << "# UPDATE > 0        update frequency for grid\n"
       << "# RL >= 0           linearize forces for distances larger than RL\n"
       << "# SMOOTH >= 0       smoothen the protein boundary after grid construction\n"
       << "#                   by SMOOTH layers\n"
       << "# NTWDF >= 0        write every NTWDF step disfield information to external file\n"
       << "# PRINTGRID = 0,1   write grid to final configuration file\n"
       << "#\n"
       << "#  NTDFR\n"
       << std::setw(8) << gin.distancefield.ntdfr << "\n"
       << "#   GRID PROTEINOFFSET PROTEINCUTOFF  PROTECT\n"
       << std::setw(8) << gin.distancefield.grid
       << std::setw(14) << gin.distancefield.proteinoffset
       << std::setw(14) << gin.distancefield.proteincutoff
       << std::setw(8) << gin.distancefield.protect << "\n"
       << "# UPDATE  SMOOTH      RL   NTWDF   PRINTGRID\n"
       << std::setw(8) << gin.distancefield.update
       << std::setw(8) << gin.distancefield.smooth
       << std::setw(8) << gin.distancefield.rl
       << std::setw(8) << gin.distancefield.ntwdf 
       << std::setw(8) << gin.distancefield.printgrid << "\n"
       << "END\n";
  }
  
  // DIHEDRALRES (promd,md++)
  if (gin.dihedralres.found) {
    os << "DIHEDRALRES\n"
            << "# NTDLR 0...3 controls dihedral-angle restraining and constraining\n"
            << "#       0:    off [default]\n"
            << "#       1:    dihedral restraining using CDLR\n"
            << "#       2:    dihedral restraining using CDLR * WDLR\n"
            << "#       3:    dihedral constraining\n"
            << "# CDLR    >=0.0 force constant for dihedral restraining\n"
            << "# PHILIN  >0.0  deviation after which the potential energy function is linearized\n"
            << "# VDIH 0,1 controls contribution to virial\n"
            << "#         0: no contribution\n"
            << "#         1: dihedral restraints contribute to virial\n"
            << "# NTWDLR >= 0 write every NTWDLRth step dist. restr. information to external file\n"
            << "# TOLDAC >= 0 tolerance for dihedral constraint\n"
            << "#          NTDLR      CDLR    PHILIN        VDIH    NTWDLR   TOLDAC\n"
            << std::setw(16) << gin.dihedralres.ntdlr
            << std::setw(10) << gin.dihedralres.cdlr
            << std::setw(10) << gin.dihedralres.philin
            << std::setw(10) << gin.dihedralres.vdih
            << std::setw(10) << gin.dihedralres.ntwdlr
            << std::setw(10) << gin.dihedralres.toldac
            << "\nEND\n";
  }
  
  // ANGLERES (promd,md++)
  if (gin.angleres.found) {
    os << "ANGLERES\n"
            << "# NTALR 0...3 controls bond-angle restraining and constraining\n"
            << "#       0:    off [default]\n"
            << "#       1:    bond-angle restraining using CALR\n"
            << "#       2:    bond-angle restraining using CALR * WALR\n"
            << "#       3:    bond-angle constraining\n"
            << "# CALR    >=0.0 force constant for bond-angle restraining\n"
            << "# VARES 0,1 controls contribution to virial\n"
            << "#         0: no contribution\n"
            << "#         1: angle restraints contribute to virial\n"
            << "# NTWALR >= 0 write every NTWDLRth step dist. restr. information to external file\n"
            << "# TOLBAC >= 0 tolerance for bond-angle constraint\n"
            << "#          NTALR      CALR       VARES      NTWALR   TOLBAC\n"
            << std::setw(16) << gin.angleres.ntalr
            << std::setw(10) << gin.angleres.calr
            << std::setw(10) << gin.angleres.vares
            << std::setw(10) << gin.angleres.ntwalr
            << std::setw(10) << gin.angleres.tolbac
            << "\nEND\n";
  }
  // JVALUERES (promd, md++)
  if (gin.jvalueres.found) {
    os << "JVALUERES\n"
            << "#        NTJVR  NTJVRA    CJVR  TAUJVR  NJVRTARS  NJVRBIQW\n"
            << std::setw(16) << gin.jvalueres.ntjvr
            << std::setw(10) << gin.jvalueres.ntjvra
            << std::setw(10) << gin.jvalueres.cjvr
            << std::setw(10) << gin.jvalueres.taujvr
            << std::setw(10) << gin.jvalueres.njvrtars
            << std::setw(10) << gin.jvalueres.njvrbiqw
            << "\n"
            << "#     LE   NGRID    DELTA     NTWJV\n"
            << std::setw(10) << gin.jvalueres.le
            << std::setw(10) << gin.jvalueres.ngrid
            << std::setw(10) << gin.jvalueres.delta
            << std::setw(10) << gin.jvalueres.write
            << "\nEND\n";
  }
  // ORDERPARAMRES (md++)
  if (gin.orderparamres.found) {
    os << "ORDERPARAMRES\n"
            << "#           NTOPR  NTOPRA  COPR   TAUOPR    UPDOPR    NTWOP\n"
            << std::setw(10) << gin.orderparamres.ntopr
            << std::setw(10) << gin.orderparamres.ntopra
            << std::setw(10) << gin.orderparamres.copr
            << std::setw(10) << gin.orderparamres.tauopr
            << std::setw(10) << gin.orderparamres.updopr
            << std::setw(10) << gin.orderparamres.ntwop
            << "\nEND\n";
  }
  // SYMRES (md++)
  if (gin.symres.found) {
    os << "SYMRES\n"
            << "#           NTSYM        CSYM\n"
            << std::setw(10) << gin.symres.ntsym
            << std::setw(10) << gin.symres.csym
            << "\nEND\n";
  }
  // LOCALELEV (promd, md++)
  if (gin.localelev.found) {
    os << "LOCALELEV\n"
            << "#     NTLES  NLEPOT  NTLESA    NTWS\n"
            << std::setw(11) << gin.localelev.ntles
            << std::setw(8) << gin.localelev.nlepot
            << std::setw(8) << gin.localelev.ntlesa
            << std::setw(8) << gin.localelev.ntwle << std::endl
            << "#    NLEPID       NTLEPFR\n";
    for (std::map<int, int>::iterator it = gin.localelev.nlepid_ntlepfr.begin();
            it != gin.localelev.nlepid_ntlepfr.end(); ++it) {
      os << std::setw(10) << it->first << std::setw(10) << it->second << std::endl;
    }
    os << "\nEND\n";
  }
  // BSLEUS (md++)
  if (gin.bsleus.found) {
    os << "BSLEUS\n"
            << "# BSLEUS   BUILD   WRITE\n"
            << std::setw(8) << gin.bsleus.bsleus
            << std::setw(8) << gin.bsleus.build
            << std::setw(8) << gin.bsleus.write;
    os << "\nEND\n";
  }
  // PERSCALE (md++)
  if (gin.perscale.found) {
    os << "PERSCALE\n"
            << "#  RESTYPE\n"
            << std::setw(10) << gin.perscale.restype
            << "\n#     KDIH       KJ         T      DIFF     RATIO      READ\n"
            << std::setw(10) << gin.perscale.kdih
            << std::setw(10) << gin.perscale.kj
            << std::setw(10) << gin.perscale.t
            << std::setw(10) << gin.perscale.diff
            << std::setw(10) << gin.perscale.ratio
            << std::setw(10) << gin.perscale.read
            << "\nEND\n";
  }

  // FREE-ENERGY CALCULATION

  // PERTURBATION (promd, md++)
  if (gin.perturbation.found) {
    os << "PERTURBATION\n"
            << "# NTG: 0 no perturbation is applied\n"
            << "#    : 1 calculate dV/dRLAM perturbation\n"
            << "#      NTG     NRDGL     RLAM     DLAMT\n"
            << std::setw(10) << gin.perturbation.ntg
            << std::setw(10) << gin.perturbation.nrdgl
            << std::setw(9) << gin.perturbation.rlam
            << std::setw(10) << gin.perturbation.dlamt
            << "\n#   ALPHLJ     ALPHC     NLAM\n"
            << std::setw(10) << gin.perturbation.alphlj
            << std::setw(10) << gin.perturbation.alphc
            << std::setw(9) << gin.perturbation.nlam
            << "\n#   NSCALE\n"
            << std::setw(10) << gin.perturbation.nscale
            << "\nEND\n";
  }

  // LAMBDAS (md++)
  if (gin.lambdas.found) {
    os << "LAMBDAS\n"
            << "#       NTIL\n"
            << std::setw(13) << gin.lambdas.ntil
            << "\n# NTLI(1...)  NILG1  NILG2    ALI    BLI    CLI    DLI    ELI\n";
    for (unsigned int i = 0; i < gin.lambdas.lambints.size(); ++i) {
      os << std::setw(13) << gin.lambdas.lambints[i].ntli
              << std::setw(7) << gin.lambdas.lambints[i].nilg1
              << std::setw(7) << gin.lambdas.lambints[i].nilg2
              << std::setw(7) << gin.lambdas.lambints[i].ali
              << std::setw(7) << gin.lambdas.lambints[i].bli
              << std::setw(7) << gin.lambdas.lambints[i].cli
              << std::setw(7) << gin.lambdas.lambints[i].dli
              << std::setw(7) << gin.lambdas.lambints[i].eli
              << "\n";
    }
    os << "END\n";
  }
  
  // PRECALCLAM (md++)
  if (gin.precalclam.found) {
    os << "PRECALCLAM\n"
            << "# NRLAM MINLAM MAXLAM\n"
            << std::setw(7) << gin.precalclam.nrlam
            << std::setw(7) << gin.precalclam.minlam
            << std::setw(7) << gin.precalclam.maxlam
            << "\nEND\n";
  }

  // UMBRELLA (promd)
  if (gin.umbrella.found) {
    os << "UMBRELLA\n"
            << "#  NTUS  USCST1  USCST2 USREF1 USREF2\n"
            << std::setw(7) << gin.umbrella.ntus
            << std::setw(7) << gin.umbrella.uscst1
            << std::setw(7) << gin.umbrella.uscst2
            << std::setw(7) << gin.umbrella.usref1
            << std::setw(7) << gin.umbrella.usref2
            << "\nEND\n";
  }

  // INPUT-OUTPUT

  // PRINTOUT (promd, md++)
  if (gin.printout.found) {
    os << "PRINTOUT\n"
            << "#NTPR: print out energies, etc. every NTPR steps\n"
            << "#NTPP: =1 perform dihedral angle transition monitoring\n"
            << "#     NTPR      NTPP\n"
            << std::setw(10) << gin.printout.ntpr
            << std::setw(10) << gin.printout.ntpp
            << "\nEND\n";
  }
  // WRITETRAJ (g96)
  if (gin.writetraj.found) {
    os << "WRITETRAJ\n"
            << "#    NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB\n"
            << std::setw(9) << gin.writetraj.ntwx
            << std::setw(10) << gin.writetraj.ntwse
            << std::setw(10) << gin.writetraj.ntwv
            << std::setw(10) << gin.writetraj.ntwf
            << std::setw(10) << gin.writetraj.ntwe
            << std::setw(10) << gin.writetraj.ntwg
            << std::setw(10) << gin.writetraj.ntwb
            << "\nEND\n";
  }
  // DEBUG (promd)
  if (gin.debug.found) {
    os << "DEBUG\n"
            << "#    NRDEO"
            << std::setw(10) << gin.debug.routines.size()
            << "\n#  PIIDER(1...NRDEO)  IIIDEO(1...NRDEO)\n";
    for (unsigned int i = 0; i < gin.debug.routines.size(); ++i) {
      os << std::setw(25) << gin.debug.routines[i].piider
              << std::setw(8) << gin.debug.routines[i].iiideo
              << "\n";
    }
    os << "END\n";
  }
  // EWARN (md++)
  if (gin.ewarn.found) {
    os << "EWARN\n"
            << "#  MAXENER\n"
            << std::setw(10) << gin.ewarn.maxener
            << "\nEND\n";
  }
  // EDS
  if (gin.eds.found && gin.eds.eds) {
    os << "EDS\n"
            << "#      EDS\n"
            << std::setw(10) << gin.eds.eds << std::endl
            << "#     FORM\n"
            << std::setw(10) << gin.eds.form << std::endl
            << "# NUMSTATES\n"
            << std::setw(10) << gin.eds.numstates << std::endl;
    switch (gin.eds.form) {
      case 1:
        os << "#        S\n"
                << std::setw(15) << gin.eds.smooth[0] << std::endl;
        break;
      case 2:
        os << "# S [1..NUMSTATES-1]\n";
        for (int N = 0; N < gin.eds.numstates*(gin.eds.numstates-1)/2; N++) {
          os << std::setw(15) << gin.eds.smooth[N];
        }
        os << std::endl;
        break;
      case 3:
        os << "# i [1..NUMSTATES-1]   j [1..NUMSTATES-1]   S [1..NUMSTATES-1]\n";
        for (int N = 0; N < (gin.eds.numstates - 1); N++) {
          os << std::setw(15) << gin.eds.tree[N][0]
                  << std::setw(15) << gin.eds.tree[N][1]
                  << std::setw(15) << gin.eds.smooth[N] << std::endl;
        }
        break;
    }
    os << "# EIR [1..NUMSTATES]\n";
    for(int N = 0; N < gin.eds.numstates; N++) {
      os << std::setw(15) << gin.eds.eir[N];
    }
    os << "\nEND\n";
  }

  // AEDS
  if (gin.aeds.found && gin.aeds.aeds) {
    os << "AEDS\n"
            << "#     AEDS\n"
            << std::setw(10) << gin.aeds.aeds << std::endl
            << "#   FORM      NUMSTATES\n"
            << std::setw(10) << gin.aeds.form
            << std::setw(15) << gin.aeds.numstates << std::endl
            << "#     EMAX      EMIN\n"
            << std::setw(10) << gin.aeds.emax
            << std::setw(10) << gin.aeds.emin << std::endl;
    os << "# EIR [1..NUMSTATES]\n";
    for(int N = 0; N < gin.aeds.numstates; N++) {
      os << std::setw(10) << gin.aeds.eir[N];
    }
         os << "\n# NTIAEDSS  RESTREMIN  BMAXTYPE      BMAX    ASTEPS    BSTEPS\n"
            << std::setw(10) << gin.aeds.ntiaedss
            << std::setw(11) << gin.aeds.restremin
            << std::setw(10) << gin.aeds.bmaxtype
            << std::setw(10) << gin.aeds.bmax
            << std::setw(10) << gin.aeds.asteps
            << std::setw(10) << gin.aeds.bsteps
            << "\nEND\n";
  }
  if (gin.reeds.found) {
    os << "REPLICA_EDS\n"
       << "#     REEDS\n"
       << std::setw(10) << gin.reeds.reeds << std::endl
       << "#   NRES    NUMSTATES       NEOFF\n"
       << std::setw(10) << gin.reeds.nres
       << std::setw(15) << gin.reeds.numstates
       << std::setw(15) << gin.reeds.neoff
       << std::endl
       << "#     RES(1 ... NRES)\n";
       for(int N = 0; N < gin.reeds.nres; N++) {
          os << std::setw(10) << gin.reeds.res[N];
       }
       os << "\n#        EIR(NUMSTATES x NRES)\n";
       int counter = 0;
       for(int N = 0; N < gin.reeds.nres * gin.reeds.numstates; N++) {
          os << std::setw(10) << gin.reeds.eir[N];
          counter += 1;
          if (counter == gin.reeds.nres){
            counter = 0;
            os << std::endl;
          }
       }
    os << "\n#        NRETRIAL        NREQUIL         CONT    EDS_STAT_OUT    PERIODIC" << std::endl
       << std::setw(10) << gin.reeds.nretrial
       << std::setw(10) << gin.reeds.nrequil
       << std::setw(10) << gin.reeds.cont
       << std::setw(10) << gin.reeds.eds_stat_out
       << std::setw(10) << gin.reeds.periodic
       << "\nEND\n";
  }

  // GAMD
  if (gin.gamd.found && gin.gamd.gamd) {
    os << "GAMD\n"
            << "#     GAMD\n"
            << std::setw(10) << gin.gamd.gamd << std::endl
            << "#     SEARCH   FORM    THRESH   NTIGAMDS\n"
            << std::setw(10) << gin.gamd.search
            << std::setw(10) << gin.gamd.form
            << std::setw(10) << gin.gamd.thresh
            << std::setw(15) << gin.gamd.ntigamds << std::endl
            << "#   AGROUPS     IGROUPS\n"
            << std::setw(10) << gin.gamd.agroups
            << std::setw(15) << gin.gamd.igroups << std::endl
            << "#     DIHSTD      TOTSTD\n"
            << std::setw(10) << gin.gamd.dihstd
            << std::setw(10) << gin.gamd.totstd << std::endl;
    os << "# ED [1..IGROUPS]\n";
    for(int N = 0; N < gin.gamd.igroups; N++) {
      os << std::setw(10) << gin.gamd.ed[N];
    }
    os << "\n";
    os << "# ET [1..IGROUPS]\n";
    for(int N = 0; N < gin.gamd.igroups; N++) {
      os << std::setw(10) << gin.gamd.et[N];
    }
    os << "\n";
    os << "# KD [1..IGROUPS]\n";
    for(int N = 0; N < gin.gamd.igroups; N++) {
      os << std::setw(10) << gin.gamd.kd[N];
    }
    os << "\n";
    os << "# KT [1..IGROUPS]\n";
    for(int N = 0; N < gin.gamd.igroups; N++) {
      os << std::setw(10) << gin.gamd.kt[N];
    }
    os << "\n# EQSTEPS  WINDOW\n"
            << std::setw(10) << gin.gamd.eqsteps
            << std::setw(10) << gin.gamd.window
            << "\nEND\n";
  }

  // ELECTRIC (md++)
  if (gin.electric.found) {
    os << "ELECTRIC\n"
            << "#     FIELD     DIPOLE     CURRENT\n"
            << std::setw(10) << gin.electric.field
            << std::setw(10) << gin.electric.dipole
            << std::setw(10) << gin.electric.current << std::endl
            << "#      EF_x       EF_y        EF_z\n"
            << std::setw(10) << gin.electric.ef_x
            << std::setw(10) << gin.electric.ef_y
            << std::setw(10) << gin.electric.ef_z << std::endl
            << "#    DIPGRP     NTWDIP\n"
            << std::setw(10) << gin.electric.dipgrp
            << std::setw(10) << gin.electric.ntwdip << std::endl
            << "#    NTWCUR    NCURGRP   CURGRP[1..NCURGRP]\n"
            << std::setw(10) << gin.electric.ntwcur
            << std::setw(10) << gin.electric.ncurgrp;
    for(int N = 0; N < gin.electric.ncurgrp; N++) {
      os << std::setw(10) << gin.electric.curgrp[N];
    }
            os << "\nEND\n";
  }

   //ADDECOUPLE (md++)
  if (gin.addecouple.found) {
    os << "ADDECOUPLE\n"
       <<"# ADGR\n"
       << std::setw(10) << gin.addecouple.adgr << "\n";
    if(gin.addecouple.adgr > 0){
      os << "# ADSTART ADEND SM SV  ST TIR\n"; 
      for(int i=0; i<gin.addecouple.adgr; ++i){
         os << std::setw(6) << gin.addecouple.adstart[i]
            << std::setw(6) << gin.addecouple.adend[i]
            << std::setw(6) << gin.addecouple.sm[i]
            << std::setw(6) << gin.addecouple.sv[i]
            << std::setw(6) << gin.addecouple.st[i]
            << std::setw(6) << gin.addecouple.tir[i]
            << "\n";
      }
      os << "# TMF STAD\n"
         << std::setw(8) << gin.addecouple.tmf
         << std::setw(8) << gin.addecouple.write
         << "\n";

    }
    os << "END\n";
  }

  //NEMD (md++)
  if (gin.nemd.found) {
  os << "NEMD\n"
     << "# NEMD     PROPERTY  METHOD\n"
     << std::setw(10) << gin.nemd.nemd 
     << std::setw(10) << gin.nemd.property 
     << std::setw(10) << gin.nemd.method << "\n"
     << "# SLABNUM  PERTFRQ    AMPLI   STDYAFT   WRITE \n"
     << std::setw(10) << gin.nemd.slabnum 
     << std::setw(10) << gin.nemd.pertfrq
     << std::setw(10) << gin.nemd.ampli 
     << std::setw(10) << gin.nemd.stdyaft
     << std::setw(10) << gin.nemd.write << "\n";   

    os << "END\n";
  }

  //QMMM (md++)
  if (gin.qmmm.found) {
  os << "QMMM\n"
     << "#   NTQMMM    NTQMSW    RCUTQM   NTWQMMM      QMLJ     QMCON   MMSCALE\n"
     << std::setw(10) << gin.qmmm.ntqmmm
     << std::setw(10) << gin.qmmm.ntqmsw
     << std::setw(10) << gin.qmmm.rcutqm
     << std::setw(10) << gin.qmmm.ntwqmmm
     << std::setw(10) << gin.qmmm.qmlj
     << std::setw(10) << gin.qmmm.qmcon
     << std::setw(10) << gin.qmmm.mmscale;
  os << "\n";  

    os << "END\n";
  }
  // EXTRA

  // Unknown blocks
  for (unsigned int i = 0; i < gin.unknown.size(); i++) {
    os << gin.unknown[i].name << "\n"
            << gin.unknown[i].content
            << "END\n";
  }


  /*
    Maybe we want to add a string-output to the following blocks?
    LAMBDAS CONSTRAINT MULTIBATH PAIRLIST PRESSURESCALE
   */


  return os;

}


void printIO(std::string block, std::string variable, std::string value, std::string allowed) {
  numErrors++;
  numTotErrors++;
  std::cout << "INPUT ERROR\n";
  std::cout << numWarnings + numErrors << ". ERROR [IO check] (" << numErrors
          << ")\n";
  std::cout << "Error in block " << block << "\n";
  std::cout << "Read " << value << " for " << variable << "\n";
  std::cout << "Accepted values are " << allowed << std::endl;
}

void printErrMsg(std::string block, std::string variable, std::string message) {
  numErrors++;
  numTotErrors++;
  std::cout << "INPUT ERROR\n";
  std::cout << numWarnings + numErrors << ". ERROR [IO check] (" << numErrors
          << ")\n";
  std::cout << "Error in block " << block << " when reading " << variable << std::endl;
  std::cout << message << std::endl;
}

void printWarning(std::string s) {
  numWarnings++;
  numTotWarnings++;
  std::cout << numWarnings + numErrors << ". WARNING (" << numWarnings << ")\n";
  std::cout << s;
  std::cout << std::endl;
}

void printError(std::string s) {
  numErrors++;
  numTotErrors++;
  std::cout << numWarnings + numErrors << ". ERROR (" << numErrors << ")\n";
  std::cout << s;
  std::cout << std::endl;
}

bool readValue(std::string BLOCK, std::string VAR, std::istringstream &is, int &variable, std::string allowed, bool allow_missing) {
  std::string s = "";
  std::stringstream ss;
  if (is.eof() == true) {
    if (allow_missing) return true;
    printIO(BLOCK, VAR, "nothing", allowed);
    return false;
  } else {
    is >> s;
  }
  if (s == "") {
    if (allow_missing) return true;
    printIO(BLOCK, VAR, "nothing", allowed);
    return false;
  } else {
    ss << s;
    ss >> variable;
    if (ss.fail() == true || ss.bad() == true) {
      printIO(BLOCK, VAR, s, allowed);
      return false;
    } else if (ss.eof() == false) {
      std::stringstream msg;
      msg << "Could not convert " << s << " to an integer";
      printErrMsg(BLOCK, VAR, msg.str());
      return false;
    }
  }
  return true;
}

bool readValue(std::string BLOCK, std::string VAR, std::istringstream &is, double &variable, std::string allowed, bool allow_missing) {
  std::string s = "";
  std::stringstream ss;
  if(is.eof() == true){
    if (allow_missing) return true;
    printIO(BLOCK, VAR, "nothing", allowed);
    return false;
  }
  else {
    is >> s;
  }
  if(s == ""){
    if (allow_missing) return true;
    printIO(BLOCK, VAR, "nothing", allowed);
    return false;
  }
  else {
    ss << s;
    ss >> variable;
    if (ss.fail() == true || ss.bad() == true) {
      printIO(BLOCK, VAR, s, allowed);
      return false;
    }
  }
  return true;
}
