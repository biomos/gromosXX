#include "blockinput.h"

#ifndef INCLUDED_IOSTREAM
#include <iostream>
#define INCLUDED_IOSTREAM
#endif

using namespace std;
using namespace io;

int main(){

  vector<string> b;

  while (getblock(cin, b).good()) {
    for (unsigned int ii = 0; ii < b.size(); ii++) 
      cout << b[ii] << "\n";
    cout << flush;
  }

  return 0;
}
