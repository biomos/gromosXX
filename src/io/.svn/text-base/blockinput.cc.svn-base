/**
 * @file blockinput.cc
 * defines blockinput functions.
 */

#include "../stdheader.h"
#include "../io/message.h"
#include "../ios"
#include "blockinput.h"

template<class size_type>
inline std::basic_string<size_type>&
trim_right( std::basic_string<size_type>& str )
{
    return( str = str.substr( 0, str.find_last_not_of( ' ' ) + 1 ) );
}

template<class size_type>
inline std::basic_string<size_type>&
trim( std::basic_string<size_type>& str )
{
  if (str.find_first_not_of( ' ' ) == std::string::npos) return (str = "");
  return( trim_right( str ) );
}


std::istream& 
io::getline(
	    std::istream& is, 
	    std::string& s, 
	    const char& sep,
	    const char& comm
	    )
{
  std::string::size_type ii;

  while (is.good()) {
    std::getline(is, s, sep);
    std::string bak(s);
    trim(bak);
    ii = s.find(comm, 0);

    if (!bak.size()) continue; // empty/whitespace only line
    else if (ii == std::string::npos) break; // no comment
    else if (!ii) continue; // comment on first position
    else {
      s.erase(s.begin() + ii, s.end());
      if (!trim_right(s).size()) continue; // line with comment only
      break;
    }
  }

  return is;
}

bool 
io::getblock(
	     std::istream& is, 
	     std::vector<std::string>& b, 
	     const std::string& sep
  )
{

  if (!b.size())
    b.push_back("");
  std::vector<std::string>::iterator dest = b.begin();

  bool first = true;
  
  while (true) {

    if (dest == b.end()) {
      b.push_back("");
      dest = b.end() - 1;
    }       
    
    getline(is, *dest);

    if (dest->find(sep) == 0)
      break;

    if (!is.good()){
      return false;
    }

    if (first){
      // first has to be a valid blockname
      // otherwise try next line
      if (trim(*dest) == "") continue;
      first = false;
    }

    ++dest;
  }

  ++dest;
  b.erase(dest, b.end());

  return true;
}

std::string& 
io::concatenate(
		std::vector<std::string>::const_iterator begin,
		std::vector<std::string>::const_iterator end,
		std::string& s,
		const char& sep
		)
{
  s.clear();
  while (begin != end) {
    s += *begin;
    s += sep;
    begin++;
  }
  
  return s;
}

void
io::trimblock(std::vector<std::string> &block)
{
  std::string s;
  std::vector<std::string>::iterator it = block.begin();
  while(true){
    std::istringstream is(*it);
    is >> s;
    trim(s);

    if (is.fail() || s == ""){
      block.erase(it);
      it = block.begin();
      continue;
    }
    else break;
  }  
}

std::string io::replace_string(std::string str, const std::string &search,
const std::string &replace) {
  std::string::size_type pos = str.find(search, 0);
  std::string::size_type length = search.length();

  while(pos != std::string::npos) {
    str.replace(pos, length, replace);
    pos = str.find(search, 0);
  }
  return str;
}


