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

/**
 * @file split_frame.cc
 * 
 */

/**
 * @page programs Program Documentation
 *
 * @anchor split_frame
 * @section split_frame 
 * @date 29. 10. 2008
 *
 */
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


  
