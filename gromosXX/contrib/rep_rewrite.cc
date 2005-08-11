
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

int read_conf(string fn,
	      vector<string> & T, vector<string> & l,
	      vector<ifstream *> & ifile,
	      vector<vector<ofstream *> > & ofile);

int init(vector<ifstream*> &traj, 
	 vector<int> &repid, vector<int> &reprun,
	 vector<string> &repT, vector<string> &repl);

int write_frame(int i, 
		vector<ifstream *> & ifile, 
		vector<vector<ofstream *> > & ofile,
		vector<string> & T,
		vector<string> & l,
		vector<int> & repid,
		vector<int> & reprun, 
		vector<string> & repT,
		vector<string> &repl,
		vector<vector<int> > & orun);

void finish(vector<ifstream * > & ifile,
	    vector<vector<ofstream *> > & ofile);

int find_ind(vector<string> const &v, string s);

int main(int argc, char *argv[])
{

  if (argc < 2){
    cout << "usage: " << argv[0]
	 << " filename"
	 << endl;
    return 1;
  }
  
  vector<string> T;
  vector<string> l;
  
  vector<vector<ofstream* > > ofile;
  vector<ifstream*> ifile;
  
  if (read_conf(argv[1], T, l, ifile, ofile)){
    cout << "error while reading configuration file!" << endl;
    return 1;
  }
  
  vector<int> repid(ifile.size(), -1);
  vector<int> reprun(ifile.size(), -1);
  vector<string> repT(ifile.size(), "");
  vector<string> repl(ifile.size(), "");

  vector<vector<int> > orun(T.size(), vector<int>(l.size(), -1));
  
  if (init(ifile, repid, reprun, repT, repl)){
    cout << "error in init!" << endl;
    return 1;
  }

  for(unsigned int i=0; i<repid.size(); ++i){
    cout << "\t" << repid[i] << " (" << reprun[i] << ")";
  }
  cout << endl;

  // select lowest run number, write to file
  bool done = false;
  
  while(!done){

    done = true;
    for(unsigned int i=0; i<reprun.size(); ++i){

      if (repid[i] == -1) continue;
      
      const int T_ind = find_ind(T, repT[i]);
      const int l_ind = find_ind(l, repl[i]);

      if (T_ind < 0 || l_ind < 0){
	cout << "\ncould not lookup indices!\n"
	     << "T_ind = " << T_ind << " ("
	     << repT[i] << ")\n"
	     << "l_ind = " << l_ind << " ("
	     << repl[i] << ")\n" << endl;
	break;
      }

      if ((reprun[i] == orun[T_ind][l_ind] + 1)){
	// ready to write!
	done = false;

	if (write_frame(i, ifile, ofile, T, l, repid, reprun, repT, repl, orun)){
	  cout << "error while writing frame!" << endl;
	  done = true;
	  break;
	}

      } // selecte a trajectory that is ready to write

    } // for all input trajectories
  } // frames
  
  finish(ifile, ofile);
  return 0;
}

int init(vector<ifstream*> &traj, vector<int> &repid, vector<int> &reprun,
	 vector<string> & repT, vector<string> & repl)
{
  string line;
  for(unsigned int i=0; i<traj.size(); ++i){

    do{
      getline(*traj[i], line);
    } while (line != "REMD");
    
    getline(*traj[i], line);
    istringstream is(line);
    is >> repid[i] >> reprun[i] >> repT[i] >> repl[i];
    
    if (is.fail()){
      cout << "error while reading file: " << i << ", line " << line << endl;
      return 1;
    }

    // and the END
    getline(*traj[i], line);
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
	    vector<vector<ofstream *> > & ofile)
{
  for(unsigned int i=0; i<ifile.size(); ++i)
    delete ifile[i];

  for(unsigned int i=0; i<ofile.size(); ++i)
    for(unsigned int j=0; j<ofile[i].size(); ++j)
      delete ofile[i][j];
}


int write_frame(int i, 
		vector<ifstream *> & ifile, 
		vector<vector<ofstream *> > & ofile,
		vector<string> & T,
		vector<string> & l,
		vector<int> & repid,
		vector<int> & reprun, 
		vector<string> & repT,
		vector<string> &repl,
		vector<vector<int> > & orun)
{
  const int T_ind = find_ind(T, repT[i]);
  const int l_ind = find_ind(l, repl[i]);

  cout << "writing frame " 
       << repid[i] << " " 
       << reprun[i] << " "
       << repT[i] << " " 
       << repl[i] << " ... ";
  cout.flush();

  if (T_ind < 0 || l_ind < 0){
    cout << "\ncould not lookup indices!\n"
	 << "T_ind = " << T_ind << " ("
	 << repT[i] << ")\n"
	 << "l_ind = " << l_ind << " ("
	 << repl[i] << ")\n" << endl;
    return 1;
  }

  if ((reprun[i] != orun[T_ind][l_ind] + 1)){
    cout << "\nwrong run number!\n"
	 << "last run " << orun[T_ind][l_ind]
	 << " current run: " << reprun[i] << endl;
    cout << "all last runs:\n";
    for(unsigned int i=0; i<orun.size(); ++i)
      for(unsigned int j=0; j<orun[i].size(); ++j)
	cout << setw(8) << orun[i][j];
    cout << endl;
    return 2;
  }
  
  string line;
  while(true){
    getline(*ifile[i], line);
    if (line == "REMD") break;

    if (ifile[i]->eof()){
      repid[i] = -1;
      ++reprun[i];
      repT[i] = "";
      repl[i] = "";
      cout << "done: end of trajectory" << endl;

      ++orun[T_ind][l_ind];

      return 0;
    }

    *ofile[T_ind][l_ind] << line << "\n";
  }
  
  getline(*ifile[i], line);
  istringstream is(line);
  int old_id = repid[i], old_run = reprun[i];
  string old_T = repT[i], old_l = repl[i];
  
  is >> repid[i] >> reprun[i] >> repT[i] >> repl[i];

  if (old_id == repid[i] && old_run == reprun[i] && old_T == repT[i] && old_l == repl[i]){
    cout << "(multiframe) ";
  }
  else
    ++orun[T_ind][l_ind];
  
  if (is.fail()){
    cout << "error while reading file: " << i << ", line " << line << endl;
    return 1;
  }
  
  // and the END
  getline(*ifile[i], line);

  cout << "done" << endl;

  return 0;
}

int find_ind(vector<string> const &v, string s)
{
  for(unsigned int i=0; i<v.size(); ++i)
    if (v[i] == s) return i;
  return -1;
}
