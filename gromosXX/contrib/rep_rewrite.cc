/**
 * @file rep_rewrite.cc
 * 
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rep_rewrite
 * @section rep_rewrite re-sort replica exchange trajectories
 * @date 29. 10. 2008
 * 
 * rep_rewrite re-sorts replica exchange trajectories according to lambda or 
 * temperature.
 * <table border=0 cellpadding=0>
 * <tr><td> \@input</td><td>&lt;input file&gt; </td></tr>
 * <tr><td> \@trj</td><td>&lt;cordinate trajectories&gt; </td></tr>
 * <tr><td> \@name</td><td>&lt;prefix and postfix of output trajectories&gt; </td></tr>
 * </table>
 */
#include <stdheader.h>

#include <simulation/multibath.h>
#include <simulation/parameter.h>

#include <io/argument.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/parameter/in_parameter.h>

using namespace std;

struct repframe
{
  int id;
  int run;
  std::string T;
  std::string l;
  int Ti;
  int li;
  int Tj;
  int lj;
  bool reeval;

  bool operator==(repframe const & r)
  {
    return (id == r.id && run == r.run && T == r.T &&
	    l == r.l && Ti == r.Ti && li == r.li &&
	    Tj == r.Tj && lj == r.lj && reeval == r.reeval);
  }
};

std::string get_name(std::string prefix, std::string suffix, repframe const & rep, simulation::Parameter const & param);

int read_repframe(ifstream & traj, repframe & stream_state);

int init(vector<ifstream*> & traj,
	 vector<repframe> & stream_state);

int write_frame(ifstream & ifile, 
		ofstream & ofile,
		repframe & stream_state,
		simulation::Parameter const & param);

void finish(vector<ifstream * > & ifile,
	    map<string, ofstream *> & ofile);

int find_ind(vector<string> const &v, string s);

int main(int argc, char *argv[])
{

  util::Known knowns;
  knowns << "input"
	 << "name"
	 << "trj";
  
  std::string usage;
  util::get_usage(knowns, usage, argv[0]);
  usage += "#\n\n";

  usage += "\t# prefix and postfix of output trajectories\n";
  usage += "\t@name       prefix  postfix\n";

  io::Argument args;

  if (args.parse(argc, argv, knowns)){
    std::cerr << usage << std::endl;
    return 1;
  }
  
  simulation::Parameter param;
  {
    std::ifstream input_file;
    input_file.open(args["input"].c_str());
    
    if (!input_file.is_open()){
      std::cerr << "\n\ncould not open " << args["input"] << "!\n" << std::endl;
      return -1;
    }
    io::In_Parameter ip(input_file);
    ip.read(param, std::cout);
  }
  
  std::string prefix, suffix;
  {
    io::Argument::const_iterator iter=args.lower_bound("name");
    io::Argument::const_iterator to=args.upper_bound("name");
    if (iter == to){
      cerr << "argument @name is required" << endl;
      return 1;
    }
    
    {
      std::istringstream is(iter->second);
      is >> prefix;
    }
    
    if (++iter == to){
      cerr << "argument @name prefix suffix is required" << endl;
      return 1;
    }
    
    {
      std::istringstream is(iter->second);
      is >> suffix;
    }
  }
  
  vector<ifstream*> ifile;
  map<string, ofstream *> ofile;
  
  {
    io::Argument::const_iterator iter=args.lower_bound("trj");
    io::Argument::const_iterator to=args.upper_bound("trj");
    
    for( ; iter != to; ++iter){
      ifstream * ifp = new ifstream(iter->second.c_str());
      if (!(ifp->is_open())){
	cerr << "could not open trajectory " + iter->second << endl;
	return 1;
      }
      
      ifile.push_back(ifp);
    }
  }

  vector<repframe> stream_state(ifile.size());
  map<string, int> out_run;

  std::cout << "initialising..." << std::endl;
  
  if (init(ifile, stream_state)){
    cout << "error in init!" << endl;
    return 1;
  }

  std::cout << "reading trajectories..." << std::endl;

  // select any frame that is ready to be written
  bool done = false;
  
  while(!done){

    done = true;
    for(unsigned int i=0; i<stream_state.size(); ++i){
      
      if (stream_state[i].id == -1) continue;
      
      std::string filename = get_name(prefix, suffix, stream_state[i], param);
      
      // check whether filename is open
      if (ofile.count(filename) == 0){
	ofstream * ofp = new ofstream(filename.c_str());
	out_run[filename] = -1;
	ofile[filename] = ofp;
	
	(*ofp) << "TITLE\n";
	
	(*ofp) << "\ttrajectory for\n";
	
	if (param.replica.temperature.size() > 1)
	  (*ofp) << "\tTi " << param.replica.temperature[stream_state[i].Ti] << "\n";

	if (param.replica.lambda.size() > 1)
	  (*ofp) << "\tlambda " << param.replica.lambda[stream_state[i].li] << "\n";
	
	if (stream_state[i].reeval){
	  (*ofp) << "\treevaluating at";
	  if (stream_state[i].Ti != stream_state[i].Tj)
	    (*ofp) << " T " << param.replica.temperature[stream_state[i].Tj];
	  if (stream_state[i].li != stream_state[i].lj)
	    (*ofp) << " lambda " << param.replica.lambda[stream_state[i].lj];
	  (*ofp) << "\n";
	}
	
	(*ofp) << "END\n";

	cout << "created " << filename << endl;
      }
      /*
      else{
	cout << "using " << filename 
	     << " run=" << stream_state[i].run
	     << " out=" << out_run[filename] + 1
	     << endl;
      }
      */

      // reeval is allowed to skip one
      if ((stream_state[i].run == out_run[filename] + 1) ||
	  (stream_state[i].reeval && stream_state[i].run == out_run[filename] + 2)){
	
	// ready to write!
	done = false;

	repframe old_state = stream_state[i];
	
	if (write_frame(*ifile[i],
			*ofile[filename],
			stream_state[i], param)){
	  cout << "error while writing frame!" << endl;
	  done = true;
	  break;
	}
	if (stream_state[i] == old_state){
	  cout << " (multiframe)" << endl;
	}
	else{
	  out_run[filename] = old_state.run;
	  cout << endl;
	}

      } // selecte a trajectory that is ready to write

    } // for all input trajectories
  } // frames
  
  finish(ifile, ofile);
  return 0;
}

int read_repframe(ifstream & traj, repframe & stream_state)
{
  std::string line;
  { // ID RUN T l
    getline(traj, line);
    istringstream is(line);
    is >> stream_state.id >> stream_state.run
       >> stream_state.T >> stream_state.l;
    if (is.fail()){
      cout << "error while reading file on line " << line << endl;
      return 1;
    }
  }
  { // Ti li Tj lj reeval
    getline(traj, line);
    istringstream is(line);
    is >> stream_state.Ti >> stream_state.li
       >> stream_state.Tj >> stream_state.lj
       >> stream_state.reeval;
    if (is.fail()){
      cout << "error while reading file on line " << line << endl;
      return 1;
    }
  }
  --stream_state.id;
  --stream_state.Ti;
  --stream_state.Tj;
  --stream_state.li;
  --stream_state.lj;
  
  // and the END
  getline(traj, line);
  return 0;
}

int init(vector<ifstream*> &traj, vector<repframe> &stream_state)
{
  string line;
  for(unsigned int i=0; i<traj.size(); ++i){

    cout << "reading first frame of trajectory " << i+1 << std::endl;
    
    do{
      getline(*traj[i], line);
    } while (line != "REMD" && !traj[i]->eof());

    if(traj[i]->eof()){
      stream_state[i].id = -1;
      cout << "done: empty trajectory file" << endl;
    }
    else if (read_repframe(*traj[i], stream_state[i]))
      return 1;
  }

  return 0;
}

  
int read_conf(string fn,
	      vector<string> & T, vector<string> & l,
	      vector<ifstream *> & ifile,
	      vector<vector<ofstream *> > & ofile)
{
  ifstream conf(fn.c_str());
  string line;
  
  do{
    getline(conf, line);
  } while (line.substr(0, 1) == "#");
  
  {
    string q;
    istringstream is(line);
    while(true){
      is >> q;
      if (is.fail()) break;
      T.push_back(q);
    }
  }
  
  do{
    getline(conf, line);
  } while (line.substr(0, 1) == "#");
  
  {
    string q;
    istringstream is(line);
    while(true){
      is >> q;
      if (is.fail()) break;
      l.push_back(q);
    }
  }
  
  const int numrep = T.size() * l.size();
  cout << "number of replicas = " << numrep << endl;

  // now the output files
  do{
    getline(conf, line);
  } while (line.substr(0, 1) == "#");

  string prefix, suffix;
  {
    istringstream is(line);
    is >> prefix >> suffix;
  }

  ofile.resize(T.size());

  for(unsigned int tt=0; tt<T.size(); ++tt){

    ofile[tt].resize(l.size());
    
    for(unsigned int ll=0; ll<l.size(); ++ll){
      
      string fn = prefix + "_" + T[tt] + "_" + l[ll] + "." + suffix;
      ofile[tt][ll] = new ofstream(fn.c_str());
      cout << "output file " << fn << " opened" << endl;
      *ofile[tt][ll] << "TITLE\n\tT = " << T[tt] << "\tl = " << l[ll]
		     << "\nEND\n";
    }
  }
  
  // finally the input files!
  while(true){
    do{
      getline(conf, line);
      if (conf.eof()) break;
    } while (line.substr(0, 1) == "#");
    if (conf.eof()) break;
    
    string fn;
    istringstream is(line);
    is >> fn;
    if (is.fail()) break;
    
    ifile.push_back(new ifstream(fn.c_str()));
    cout << "input file " << fn << " opened" << endl;
  }
  
  return 0;
}

void finish(vector<ifstream * > & ifile,
	    map<string, ofstream *> & ofile)
{
  for(unsigned int i=0; i<ifile.size(); ++i){
    ifile[i]->close();
    delete ifile[i];
  }

  map<string, ofstream *>::iterator
    it = ofile.begin(),
    to = ofile.end();
  
  for( ; it != to; ++it){
    it->second->close();
    delete(it->second);
  }
}


int write_frame(ifstream & ifile, 
		ofstream & ofile,
		repframe & stream_state,
		simulation::Parameter const & param)
{
  cout << "writing frame "
       << setw(4) << stream_state.id << " " 
       << setw(4) << stream_state.run;
  
  cout.precision(3);
  cout.setf(ios::fixed, ios::floatfield);
  
  if (param.replica.temperature.size() > 1)
    cout << " " << setw(10) << param.replica.temperature[stream_state.Ti];
  
  if (param.replica.lambda.size() > 1)
    cout << " " << setw(10) << param.replica.lambda[stream_state.li];
  
  if (stream_state.reeval)
    cout << " (reeval)";

  cout << " ... ";
  cout.flush();
  
  string line;
  while(true){
    getline(ifile, line);
    if (line == "REMD") break;

    if (ifile.eof()){
      stream_state.id = -1;
      ++stream_state.run;
      
      cout << "done: end of trajectory" << endl;
      return 0;
    }

    ofile << line << "\n";
  }
  
  read_repframe(ifile, stream_state);
  cout << "done";

  return 0;
}

std::string get_name(std::string prefix, std::string suffix, repframe const & rep, simulation::Parameter const & param)
{
  std::ostringstream os;
  os << prefix;

  if (param.replica.temperature.size() > 1)
    os << "_" << param.replica.temperature[rep.Ti];
  
  if (param.replica.lambda.size() > 1)
    os << "_" << param.replica.lambda[rep.li];
  
  if (rep.reeval){
    os << "_reeval";
    
    if (rep.Ti != rep.Tj)
      os << "_" << param.replica.temperature[rep.Tj];
    if (rep.li != rep.lj)
      os << "_" << param.replica.lambda[rep.lj];
  }
  
  os << "." << suffix;

  return os.str();
}

