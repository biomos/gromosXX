#ifndef INCLUDED_IO_BLOCKINPUT
#define INCLUDED_IO_BLOCKINPUT

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

#ifndef INCLUDED_ISTREAM
#include <istream>
#define INCLUDED_ISTREAM
#endif

#ifndef INCLUDED_SSTREAM
#include <sstream>
#define INCLUDED_SSTREAM
#endif

#ifndef INCLUDED_STDEXCEPT
#include <stdexcept>
#define INCLUDED_STDEXCEPT
#endif

using namespace std;

namespace io {

  /*
   * The function io::getline provides an override of 
   * std::getline in that it retrieves the next line (separated 
   * by const char& sep) from the input stream, which, after 
   * having been stripped of comments (indicated by 
   * const char& comm), is not empty.
   */
  inline istream& getline(
    istream& is, 
    string& s, 
    const char& sep = '\n',
    const char& comm = '#'
  );

  /*
   * The function io::getblock retrieves the next block  
   * (separated by const string& sep) using io::getline.
   *
   * It throws a runtime_error if the stream's good bit is unset 
   * before it's finished reading a block.
   *
   * Finally, the vector<string> it writes to is resized to the 
   * number of strings read.
   */
  inline istream& getblock(
    istream& is, 
    vector<string>& b, 
    const string& sep = "END"
  );

  /*
   * The io::concatenate utility allows the concatenation
   * of entries in a vector<string> into just one string. The
   * resulting entries are separated by const char& sep 
   * (defaulting to a newline character).
   */
  inline
  string& concatenate(
    vector<string>::iterator begin,
    vector<string>::iterator end,
    string& s 
  );

}

#include "blockinput.tcc"

#endif
