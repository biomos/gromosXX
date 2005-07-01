
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
  double Ti, Tj;
  double li, lj;
  double epot_i, epot_j;
  double p;
  int s;
};

int find_replica(vector<Replica_Data> const & v, double T, double l);
int find_lambda(vector<double> l, double lam);
int find_temperature(vector<double> T, double temp);

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

  vector<int> lastT(num_rep);
  vector<int> lastl(num_rep);

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
  
  while(true){
    getline(rep, l);
    if (rep.eof()) break;
    
    is.clear();
    is.str(l);
    
    is >> r.ID >> r.run >> r.Ti >> r.li >> r.epot_i
       >> r.Tj >> r.lj >> r.epot_j >> r.p >> r.s;

    if (r.run == 0){
      lastT[r.ID-1] = find_temperature(T, r.Ti);
      lastl[r.ID-1] = find_lambda(lam, r.li);
    }
    else{

      int newT = find_temperature(T, r.Ti);
      int newl = find_lambda(lam, r.li);

      if (abs(lastT[r.ID-1] - newT) + abs(lastl[r.ID-1] - newl) > 1){
	cout << "file corrupt: change T "
	     << lastT[r.ID-1] << " -> " << newT << " and l "
	     << lastl[r.ID-1] << " -> " << newl << "\n";
      }
      lastT[r.ID-1] = newT;
      lastl[r.ID-1] = newl;
    }
  
    if ((unsigned)r.run >= rep_data.size())
      rep_data.push_back(entry);
    
    rep_data[r.run][r.ID-1] = r;
    
  }
  
      cerr << "read " << rep_data.size() << " runs.\n";

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
    vector<double> probT(num_T, 0.0);
    vector<vector<double> > probTv(num_T);
    vector<int> swiT(num_T, 0);
    vector<int> trT(num_T, 0);

    vector<double> probl(num_l, 0.0);
    vector<vector<double> > problv(num_l);
    vector<int> swil(num_l, 0);
    vector<int> trl(num_l, 0);

    ofstream epot("epot.dat");
    epot.precision(4);
    epot.setf(ios::fixed, ios::floatfield);

    ofstream probfile("probability.dat");
    probfile.precision(4);
    probfile.setf(ios::fixed, ios::floatfield);
    
    ofstream switchedfile("switches.dat");

    ofstream prob_l("prob_l.dat");
    prob_l.precision(4);
    prob_l.setf(ios::fixed, ios::floatfield);

    ofstream prob_T("prob_T.dat");
    prob_T.precision(4);
    prob_T.setf(ios::fixed, ios::floatfield);

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

	  if (rep_data[i][rr].li < rep_data[i][rr].lj){
	    // lambda switch (up)

	    const int lind = find_lambda(lam, rep_data[i][rr].li);
	    if (lind < 0){
	      cerr << "could not find lambda!" << endl;
	      return 1;
	    }
	    
	    probl[lind] += rep_data[i][rr].p;
	    problv[lind].push_back(rep_data[i][rr].p);

	    swil[lind] += rep_data[i][rr].s;
	    ++trl[lind];
	  }

	  if (rep_data[i][rr].Ti < rep_data[i][rr].Tj){
	    // temperature switch (up)
	    const int Tind = find_temperature(T, rep_data[i][rr].Ti);
	    if (Tind < 0){
	      cerr << "could not find temperature!" << endl;
	      return 1;
	    }
	    
	    probT[Tind] += rep_data[i][rr].p;
	    probTv[Tind].push_back(rep_data[i][rr].p);

	    swiT[Tind] += rep_data[i][rr].s;
	    ++trT[Tind];
	  }
	  
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

    prob_l << "# switching probabilities per lambda\n\n";
    unsigned int q = 0;
    while(true){
      bool br = true;
      for(int i=0; i<num_l-1; ++i){
	if (problv[i].size() > q){
	  br = false;
	  prob_l << std::setw(12) << problv[i][q];
	}
	else{
	  prob_l << std::setw(12) << "-";
	}
      }
      prob_l << "\n";
      if (br) break;
      ++q;
    }

    prob_T << "# switching probabilities per temperature\n\n";
    q = 0;
    while(true){
      bool br = true;
      for(int i=0; i<num_T-1; ++i){
	if (probTv[i].size() > q){
	  br = false;
	  prob_T << std::setw(12) << probTv[i][q];
	}
	else{
	  prob_T << std::setw(12) << "-";
	}
      }
      prob_T << "\n";
      if (br) break;
      ++q;
    }
    
    prob_l.close();
    prob_T.close();

    cout << "average switching probabilities:\n"
	 << "\tTemperature:\n";
    
    double ppp = 0;
    for(unsigned int rr=0; rr<probT.size() - 1; ++rr){
      if (trT[rr] > 0){
	cout << setw(12) << probT[rr] / trT[rr];
	ppp += probT[rr] / trT[rr];
      }
      else
	cout << setw(12) << 0.0;
    }
    cout << "\n\ntotal average switching probability = ";
    if (probT.size() > 1)
      cout << ppp / (probT.size() - 1);
    else
      cout << 0;

    cout << "\n\nnumber of switches:\n";
    int sss = 0;
    for(unsigned int rr=0; rr<swiT.size() - 1; ++rr){
      cout << setw(12) << swiT[rr];
      sss += swiT[rr];
    }
    cout << "\n\ntotal number of switches = " << sss;

    cout << "\n\naverage switching probabilities:\n"
	 << "\tlambda:\n";
    
    ppp = 0;
    for(unsigned int rr=0; rr<probl.size() - 1; ++rr){
      if (trl[rr] > 0){
	cout << setw(12) << probl[rr] / trl[rr];
	ppp += probl[rr] / trl[rr];	
      }
      else
	cout << setw(12) << 0.0;
    }
    cout << "\n\ntotal average switching probability = ";
    if (probl.size() > 1)
      cout << ppp / (probl.size() - 1);
    else
      cout << 0.0;

    cout << "\n\nnumber of switches:\n";
    sss = 0;
    for(unsigned int rr=0; rr<swil.size() - 1; ++rr){
      cout << setw(12) << swil[rr];
      sss += swil[rr];
    }
    cout << "\n\ntotal number of switches = " << sss;
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

int find_lambda(vector<double> l, double lam)
{
  for(unsigned int i=0; i<l.size(); ++i){
    if (lam == l[i]) return i;
  }
  return -1;
}

int find_temperature(vector<double> T, double temp)
{
  for(unsigned int i=0; i<T.size(); ++i){
    if (temp == T[i]) return i;
  }
  return -1;
}
