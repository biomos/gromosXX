
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 2){
    cout << "\n" << argv[0] << " filename\n" << endl;
    return 1;
  }
  
  ifstream fin(argv[1]);
  if (!fin.is_open()){
    cout << "\ncould not open file " << argv[1] << "\n" << std::endl;
    return 1;
  }

  string s;
  do{
    fin >> s;
  } while (s != "REPLICAFRAME");
  
  int r;
  do{

    fin >> r; // replica number
    ostringstream os;
    os << "rep_" << r << ".crd";
    ofstream fout(os.str().c_str());
    getline(fin, s);
    getline(fin, s); // END
    getline(fin, s); // next block

    cout << "writing replica " << r << " ... ";
    cout.flush();
    
    fout << "TITLE\n\treplica " << r << "\nEND\n";

    while(s != "REPLICAFRAME"){
      fout << s << "\n";
      getline(fin, s);
      // fin >> s;
      if (fin.eof()) break;
    }
    fout.close();
    
    cout << "done" << std::endl;

  } while(!fin.eof());

  fin.close();
  return 0;
}


  
