
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

using namespace std;

struct Replica_Data
{
  int ID;
  int run;
  double Ti;
  double li;
  double epot_i;
  double p;
  int s;
};

int find_replica(vector<Replica_Data> const & v, double T, double l);

int main(int argc, char *argv[])
{
  if (argc < 2){
    cerr << "usage: " << argv[0] << " filename\n\n";
    return 1;
  }
  
  ifstream rep(argv[1]);
  if (!rep.is_open()){
    cerr << "could not open file!" << endl;
    return 1;
  }
  
  string l;
  getline(rep, l);

  int num_T, num_l;
  string s;
  
  istringstream is(l);
  is >> s >> num_T;
  
  getline(rep, l);
  is.clear();
  is.str(l);
  is >> s >> num_l;
  
  const int num_rep = num_T * num_l;
  
  vector<Replica_Data> entry(num_rep);
  vector<vector<Replica_Data> > rep_data;

  vector<double> T;
  vector<double> lam;
  
  getline(rep, l); // temperatures
  is.clear();
  is.str(l);
  double d;
  is >> s;
  for(int i=0; i<num_T; ++i){
    is >> d;
    T.push_back(d);
  }
    
  getline(rep, l); // lambdas
  is.clear();
  is.str(l);
  is >> s;
  for(int i=0; i<num_l; ++i){
    is >> d;
    lam.push_back(d);
  }

  getline(rep, l); // empty
  getline(rep, l); // header
  
  Replica_Data r;
  double Tj, lj, epot_j;
  
  while(true){
    getline(rep, l);
    if (rep.eof()) break;
    
    is.clear();
    is.str(l);
    
    is >> r.ID >> r.run >> r.Ti >> r.li >> r.epot_i
       >> Tj >> lj >> epot_j >> r.p >> r.s;
  
    if ((unsigned)r.run > rep_data.size())
      rep_data.push_back(entry);
    
    rep_data[r.run-1][r.ID-1] = r;
  }
  
  std::cerr << "read " << rep_data.size() << " runs.\n";

  // write the files
  {
    ofstream temp("temperature.dat");
    temp.precision(4);
    temp.setf(ios::fixed, ios::floatfield);
    
    for(unsigned int i=0; i<rep_data.size(); ++i){
      temp << std::setw(10) << i;
      for(int rr=0; rr<num_rep; ++rr)
	temp << std::setw(12) << rep_data[i][rr].Ti;
      temp << "\n";
    }
    temp.close();
  }
  
  {
    ofstream lambda("lambda.dat");
    lambda.precision(4);
    lambda.setf(ios::fixed, ios::floatfield);
    
    for(unsigned int i=0; i<rep_data.size(); ++i){
      lambda << std::setw(10) << i;
      for(int rr=0; rr<num_rep; ++rr)
	lambda << std::setw(12) << rep_data[i][rr].li;
      lambda << "\n";
    }
    lambda.close();
  }
  
  {
    vector<double> prob(num_rep, 0.0);
    vector<int> swi(num_rep, 0);

    ofstream epot("epot.dat");
    epot.precision(4);
    epot.setf(ios::fixed, ios::floatfield);

    ofstream probfile("probability.dat");
    probfile.precision(4);
    probfile.setf(ios::fixed, ios::floatfield);
    
    ofstream switchedfile("switches.dat");

    for(unsigned int i=0; i<rep_data.size(); ++i){
      
      epot << std::setw(10) << i;
      probfile << std::setw(10) << i;
      switchedfile << std::setw(10) << i;

      int index = 0;
      for(unsigned int tt=0; tt<T.size(); ++tt){
	for(unsigned int ll=0; ll<lam.size(); ++ll, ++index){

	  int rr = find_replica(rep_data[i], T[tt], lam[ll]);
	  if (rr == -1){
	    cerr << "could not find replica for T=" << T[tt]
		 << " and l=" << lam[ll] << endl;
	    return 1;
	  }
	  
	  prob[rr] += rep_data[i][rr].p;
	  swi[rr] += rep_data[i][rr].s;
	  
	  epot << std::setw(18) << rep_data[i][rr].epot_i;
	  probfile << std::setw(12) << rep_data[i][rr].p;
	  switchedfile << std::setw(4) << rep_data[i][rr].s;

	}
      }

      epot << "\n";
      probfile << "\n";
      switchedfile << "\n";

    }
    epot.close();
    probfile.close();
    switchedfile.close();

    cout << "average switching probabilities:\n";
    double ppp = 0;
    for(unsigned int rr=0; rr<prob.size(); ++rr){
      cout << setw(12) << prob[rr] / num_rep;
      ppp += prob[rr];
    }
    cout << "\n\ntotal average switching probability = " 
	 << ppp / num_rep / num_rep;

    cout << "\n\nnumber of switches:\n";
    int sss = 0;
    for(unsigned int rr=0; rr<swi.size(); ++rr){
      cout << setw(12) << swi[rr];
      sss += swi[rr];
    }
    cout << "\n\ntotal number of switches = " << sss / 2;
  }

  cout << "\n\n" << endl;
  
  return 0;
}


int find_replica(vector<Replica_Data> const & v, double T, double l)
{
  for(unsigned int i=0; i<v.size(); ++i){
    if (v[i].Ti == T &&
	v[i].li == l)
      return i;
  }
  return -1;
}

