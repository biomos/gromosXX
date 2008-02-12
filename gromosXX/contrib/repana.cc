/**
 * @file repana.cc
 * analyses the REMD output file replica.dat
 * and suggests an optimized temperature
 * and lambda set.
 */

/**
 * @page programs Program Documentation
 *
 * @anchor repana
 * @section repana analyze replica.dat
 * @author @ref mc cc
 * @date 3. 12. 2007
 *
 * repana extracts the information stored in the replica
 * exchange molecular dynamics (REMD) output file replica.dat.
 * It produces seven multi column output files:
 * temperature.dat: run number vs. temperature of replica 
 * (column 2 corresponds to replica 1, column 3 to replica 2 etc.)
 * lambda.dat: run number vs. lambda value of replica
 * epot.dat: run number vs. potential energy of replica
 * probability.dat: run number vs. switching probability of replica
 * switches.dat: run number vs. switching data of replica 
 * (0 = no switch in this run, 1 = switch in this run)
 * prob_T.dat: switching probabilities per temperature
 * prob_l.dat: switching probabilities per lambda
 *
 * Furthermore it calculates an optimized temperature or lambda set
 * based on the fraction of replicas diffusing from the lowest
 * to the highest temperature (lambda value). This algorithm is 
 * based on:
 * Katzgraber, H. G.; Trebst, S.; Huse, D. A. & Troyer, M.
 * Feedback-optimized parallel tempering Monte Carlo
 * Journal of Statistical Mechanics-Theory And Experiment,
 * Iop Publishing Ltd, 2006, P03018"
 * 
 * The 'opt' function is based on a perl script by H.G. Katzgraber.
 *
 * <b>usage:</b>
 * <table border=0 cellpadding=0>
 * <tr><td>./repana replica.dat</td><td>&lt;REMD output file&gt; </td></tr>
 * </table>
 *
 */

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

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

struct hist
{
  hist(): nup(0.0),ndown(0.0),f(0.0){}
  double nup;
  double ndown;
  double f;
};

int find_replica(vector<Replica_Data> const & v, double T, double l);
int find_lambda(vector<double> l, double lam);
int find_temperature(vector<double> T, double temp);
void fraction(vector<vector<Replica_Data> > rep_data, vector<double> T, vector<double> lam,ostream &os = cout);
void opt(vector<hist> f_T, vector<double>T,ostream &os = cout);

int main(int argc, char *argv[])
{
  if (argc < 2){
    cerr << "usage: " << argv[0] << " replica.dat\n\n";
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
  
      cout << "# read " << rep_data.size() << " runs.\n";

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
		 << " and l=" << lam[ll] << " run=" << i+1 <<  endl;
	    cerr << "there are only " << rep_data[i].size() << " replicas this far" << endl;
	    cerr << "(T/l) :\n";
	    vector<string> sval, sind;
	    ostringstream os;
	    
	    os.precision(4);
	    os.setf(ios::fixed, ios::floatfield);
	    
	    for(unsigned int rrr=0; rrr<rep_data[i].size(); ++rrr){
	      os.clear();
	      os.str("");
	      os << "(" << std::setw(10) << rep_data[i][rrr].Ti << "," << std::setw(10) << rep_data[i][rrr].li << ")";
	      sval.push_back(os.str());
	      os.clear();
	      os.str("");
	      os << "(" << std::setw(3) << find_temperature(T, rep_data[i][rrr].Ti) 
		 << "," << std::setw(3) << find_lambda(lam, rep_data[i][rrr].li) << ")";
	      sind.push_back(os.str());
	    }

	    sort(sval.begin(), sval.end());
	    sort(sind.begin(), sind.end());

	    for(unsigned int rrr=0; rrr<sval.size(); ++rrr){
	      cerr << std::setw(25) << sval[rrr];
	      if (((rrr + 1) % 5) == 0) cerr << endl;
	    }
	    cerr << endl;
	    for(unsigned int rrr=0; rrr<rep_data[i].size(); ++rrr){
	      cerr << std::setw(10)  << sind[rrr];
	      if (((rrr + 1) % 5) == 0) cerr << endl;
	    }
	    cerr << endl;

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

    cout << "# average switching probabilities:" << endl;
        
    cout << setw(12) << "# T" << setw(12) << "P(T)" << endl;
    double ppp = 0;
    for(unsigned int rr=0; rr<probT.size() - 1; ++rr){
      cout << setw(12) << T[rr];
      if (trT[rr] > 0){
	cout << setw(12) << probT[rr] / trT[rr];
	ppp += probT[rr] / trT[rr];
      }
      else
	cout << setw(12) << 0.0;
      cout << endl;
    }
    cout << "\n\n# total average switching probability = ";
    if (probT.size() > 1)
      cout << ppp / (probT.size() - 1);
    else
      cout << 0;

    cout << "\n\n# number of switches:\n";
    cout << setw(12) << "# T" << setw(12) << "switches" << endl;
    int sss = 0;
    for(unsigned int rr=0; rr<swiT.size() - 1; ++rr){
      cout << setw(12) << T[rr]
	   << setw(12) << swiT[rr] << endl;
      sss += swiT[rr];
    }
    cout << "\n\n# total number of switches = " << sss;

    cout << "\n\n# average switching probabilities:" << endl;

    cout << setw(12) << "# l" << setw(12) << "P(l)" << endl;
    ppp = 0;
    for(unsigned int rr=0; rr<probl.size() - 1; ++rr){
      if (trl[rr] > 0){
	cout << setw(12) << lam[rr];
	cout << setw(12) << probl[rr] / trl[rr] << endl;
	ppp += probl[rr] / trl[rr];	
      }
      else
	cout << setw(12) << 0.0;
    }
    cout << "\n\n# total average switching probability = ";
    if (probl.size() > 1)
      cout << ppp / (probl.size() - 1);
    else
      cout << 0.0;

    cout << "\n\n# number of switches:\n";
    cout << setw(12) << "# l" << setw(12) << "switches" << endl;
    sss = 0;
    for(unsigned int rr=0; rr<swil.size() - 1; ++rr){
      cout << setw(12) << lam[rr]; 
      cout << setw(12) << swil[rr] << endl;
      sss += swil[rr];
    }
    cout << "\n\n# total number of switches = " << sss;
  }

  cout << "\n\n" << endl;
  
  fraction(rep_data, T, lam);
  
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
void fraction(vector<vector<Replica_Data> > rep_data, vector<double> T, vector<double> lam,ostream &os)
{
  os << "# Optimization of temperature and lambda set\n"
     << "# according to:\n"
     << "# Katzgraber, H. G.; Trebst, S.; Huse, D. A. & Troyer, M. \n"
     << "# Feedback-optimized parallel tempering Monte Carlo \n"
     << "# Journal of Statistical Mechanics-Theory And Experiment, \n"
     << "# Iop Publishing Ltd, 2006, P03018" << endl;

  // get min and max temperature and lambda values
  double minT=*(T.begin());
  double maxT=*(T.end()-1);
  double minlam=*(lam.begin());
  double maxlam=*(lam.end()-1);
    
  // vector of labels
  vector<std::string> T_label(rep_data[0].size(),"none");
  vector<std::string> lam_label(rep_data[0].size(),"none");

  // map temperature and lambda values to index
  map<double,unsigned int> map_T;
  for(unsigned int i=0; i<T.size();i++){
    map_T[T[i]]=i;
  }
  map<double,unsigned int> map_lam;
  for(unsigned int i=0; i<lam.size();i++){
    map_lam[lam[i]]=i;
  }
  
  // record the histograms
  hist temp;
  vector<hist> f_T(rep_data[0].size(),temp);
  vector<hist> f_lam(rep_data[0].size(),temp);

  // loop over all runs
  for(unsigned int i=0; i<rep_data.size(); i++){
    // loop over the number of replicas
    for(unsigned int rr=0; rr<rep_data[i].size(); rr++){
      // temperature 
      int index =map_T[rep_data[i][rr].Ti];
      if(T_label[rr]=="up"){
	if(rep_data[i][rr].Ti==maxT){
	  T_label[rr]="down";
	  f_T[index].ndown++;
	}
	else
	  f_T[index].nup++;
      }
      else if (T_label[rr]=="down"){
	if(rep_data[i][rr].Ti==minT){
	  T_label[rr]="up";
	  f_T[index].nup++;
	}
	else
	  f_T[index].ndown++;
      }
      else if (T_label[rr]=="none"){
	if(rep_data[i][rr].Ti==minT)
	  T_label[rr]="up";
	else if(rep_data[i][rr].Ti==maxT)
	  T_label[rr]="down";
      }
      // lambda 
      index=map_lam[rep_data[i][rr].li];
      if(lam_label[rr]=="up"){
	if(rep_data[i][rr].li==maxlam){
	  lam_label[rr]="down";
	  f_lam[index].ndown++;
	}
	else
	  f_lam[index].nup++;
      }
      else if (lam_label[rr]=="down"){
	if(rep_data[i][rr].li==minlam){
	  lam_label[rr]="up";
	  f_lam[index].nup++;
	}
	else
	  f_lam[index].ndown++;
      }
      else if (lam_label[rr]=="none"){
	if(rep_data[i][rr].li==minlam)
	  lam_label[rr]="up";
	else if(rep_data[i][rr].li==maxlam)
	  lam_label[rr]="down";
      }
    } // loop over all replicas
  } // loop over all runs
  
  // output results
  os.precision(4);
  // temperature
  if(T.size()>1){
    os << "\n# fraction T" << endl;
    os << setw(12) << "# ID" 
	 << setw(12) << "f(T)"
	 << setw(12) << "T" << endl;
    
    for(unsigned int i=0;i<T.size();i++){
      double sum = f_T[i].nup+f_T[i].ndown;
      if(sum!=0)
	f_T[i].f=f_T[i].nup/sum;
      os << setw(12) << i+1 
	   << setw(12) << f_T[i].f
	   << setw(12) << T[i] << endl;
    }
    os << "\n# optimized temperature set" << endl;
    opt(f_T,T,os);
    os << endl;
  }
  
  // lambda
  if(lam.size()>1){
    os << "\n# fraction lambda" << endl;
    os << setw(12) << "# ID" 
	 << setw(12) << "f(lam)"
	 << setw(12) << "lam" << endl;
    for(unsigned int i=0;i<lam.size();i++){
      double sum = f_lam[i].nup+f_lam[i].ndown;
      if(sum!=0)
	f_lam[i].f=f_lam[i].nup/sum;
      os << setw(12) << i+1 
	   << setw(12) << f_lam[i].f
	   << setw(12) << lam[i] << endl;            
    }
    // optimize the lambda set
    os << "\n# optimized lambda set" << endl;
    opt(f_lam,lam,os);
    os << endl;
  }
 
}

void opt(vector<hist> f_T, vector<double>T,ostream &os)
{
  os.precision(5);
  // optimize the temperature (or lambda) set
  vector<double> fnew;
  vector<double> Tnew;
  for(unsigned int i=0;i<T.size();i++){
    if(i==(T.size()-1) || f_T[i+1].f!=f_T[i].f){
      fnew.push_back(f_T[i].f);
      Tnew.push_back(T[i]);
    }    
  }
  vector<double> dTnew(fnew.size(),0.0);
  vector<double> derivative(fnew.size(),0.0);
  vector<double> dTnewprime(fnew.size(),0.0);
  for(unsigned int i=0;i<(fnew.size()-1);i++){
    dTnew[i]=Tnew[i+1]-Tnew[i];
    derivative[i]=(fnew[i+1]-fnew[i])/dTnew[i];
    dTnewprime[i]=1.0/sqrt(fabs((1.0/dTnew[i])*derivative[i]));        
  }

  double constant=0.0;
  for(unsigned int i=0;i<fnew.size()-1;i++){
    constant+=dTnew[i]/dTnewprime[i];
  }
  constant = (T.size()-1)/constant;
  unsigned int n = 0;
  double n_prime = 0.0;
  unsigned int t = 0;
  double this_dT = 0.0;
  int dtc = 0;
  vector<double> Tnewprime(T.size(),0.0);
  Tnewprime[0]=T[0];
  do{
    if(n_prime+constant*dTnew[t]/dTnewprime[t] >= n+1){
      double tau=dTnewprime[t]/constant*(n-n_prime+1);
      if(tau>dTnew[t]){
	cerr << "ERROR: Inconsistent tau " << tau << endl;
      }
      this_dT+=tau;
      Tnewprime[dtc+1] = Tnewprime[dtc] + this_dT;
      os << setw(12) << Tnewprime[dtc] << endl;
      dtc++;
      n++;
      n_prime+=constant*tau/dTnewprime[t];
      dTnew[t]-=tau;
      if(dTnew[t] < 0){
	os << "ERROR: Inconsistent deltaT" << dTnew[t] << endl;
      }
      this_dT = 0;
    }
    else{
      n_prime += constant * dTnew[t] / dTnewprime[t];
      this_dT += dTnew[t];
      t++;
    }             
  }
  while(n<(T.size()-1) && t < (fnew.size()-1));
  // check whether the following two lines are really necessary      
  os << "# WARNING: printing fishy line for second to last T value" << endl;
  os << setw(12) << Tnewprime[T.size()-2] << endl;
  os << setw(12) << T[T.size()-1] << endl;
}

  
