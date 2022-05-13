/**
 * @file blockinput.cc
 * defines blockinput functions.
 */

#include <cstring>
#include "../stdheader.h"
#include "../io/message.h"
#include <ios>
#include "blockinput.h"

#undef MODULE
#define MODULE io
#undef SUBMODULE
#define SUBMODULE parameter

std::istream& 
io::getline(
	    std::istream& is, 
	    std::string& s, 
	    const char& sep,
	    const char& comm
	    )
{
  std::string::size_type ii = 0;

  while (is.good()) {
    std::getline(is, s, sep);
    rtrim(s);
    ii = s.find(comm, 0);

    if (!s.size()) continue; // empty/whitespace only line
    else if (ii == std::string::npos) break; // no comment
    else if (!ii) continue; // comment on first position
    else {
      s.erase(s.begin() + ii, s.end());
      rtrim(s);
      if (!s.size()) continue; // line with comment only
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
      rtrim(*dest);
      if (*dest == "") continue;
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

std::string
io::concatenate(
                std::vector<std::string>::const_iterator begin,
                std::vector<std::string>::const_iterator end,
                const char& sep
                )
{
  std::string s;
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


/**
 * split a string at delimiter
 * @returns vector of substrings
 */
std::vector<std::string> io::split(const std::string &s, std::string delim) {

    std::vector<std::string> elems;
    size_t start=0, pos=0; 
    std::string item;

    while (( pos=s.find(delim,start)) != std::string::npos) {
      item = s.substr(start, pos-start);
      elems.push_back(item);
      start=pos+delim.length();
    } 
    // last element 
    item = s.substr(start, s.length()-start);
    elems.push_back(item); 

    return elems;
}

bool io::is_valid_int(const char* x, bool sign)
{
	bool checked = true;

	int i = 0;
	do
	{   if (sign && i==0 && x[i] == '-') i++;
        else if (isdigit(x[i])) i++;
		else {
			i++;
			checked = false;
			break;
		}
	} while (x[i] != '\0');

	return checked;
}


template <>
bool io::is_valid_type <int> (int & var, const char* x) {
  return is_valid_int(x, true);
}

template<>
bool io::is_valid_type<unsigned int>(unsigned int & var, const char* x) {
  return is_valid_int(x, false);
}

template <>
bool io::is_valid_type <bool> (bool & var, const char* x) {
  if (std::strcmp(x,"1")==0 ||std::strcmp(x,"0")==0 ) return true;
  else return false;
}

int io::Block::read_buffer(std::vector<std::string> &buffer, bool required) {
  //check if buffer exists and there is something between blockname and END
  _numlines=buffer.size();
  if (buffer.size() == 0 && !required) return 1; 

  if (buffer.size() <= 2){
    if (required) {
      std::string separator="#---------------------------------------\n";
      io::messages.add(_blockname+" block: missing or empty, but required\n"+separator+"# Example block:\n"+separator+_exampleblock+separator,
		     "blockinput", io::message::error);
      _block_error=1;
      return 1;
    } else {
      io::messages.add("empty block found: "+_blockname,
		     "blockinput", io::message::warning);
      return 1;
    }
  } else {
    // we have something in the block
    //std::cerr << "# "<<buffer[0] << std::endl;
    
    if (_blockname != buffer[0])
      io::messages.add(_blockname+" block: wrong block format, the first entry has to be the blockname but is: "+buffer[0],
		     "blockinput", io::message::error);
    //_lineStream.clear();  
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1));
    //std::cerr << concatenate(buffer.begin() , buffer.end()) << std::endl;
    DEBUG(11, _blockname+" present");
    return 0;
  }  
}

void io::Block::get_final_messages(bool print_read) {
    if (_block_error && print_read) {
      if (print_read) {
          std::string readparameters = get_readparameters();
          std::string separator="#---------------------------------------\n";
          io::messages.add(_blockname +" block: Parameters read:\n"+readparameters+"\n"+separator+"# Example block:\n"+separator+_exampleblock+separator,
             "In_Parameter", io::message::error);
      } else {
             io::messages.add(_blockname +" block: not enough values or wrong value type.", 
             "In_Parameter", io::message::error);
      }
    } else {
      std::string leftover;
      _lineStream >> leftover;
      if (!_lineStream.eof())
        io::messages.add("Left-over parameters in "+ _blockname +" block: "+leftover+"\n",
             "In_Parameter", io::message::error);
    }
    // temporary for debugging:
    //io::messages.add(get_readparameters()+"\n", "In_Parameter", io::message::warning);
}

std::string io::Block::get_readparameters() {
  std::stringstream oss;
  //oss << blockname << "\n";
  
  unsigned int entries_per_line=5;
  bool end=false;
  unsigned int j=0;
  while (!end){
    oss << "    # ";
    for (unsigned int i=j*entries_per_line; i < (j+1)*entries_per_line; i++) {
      if (i<_par_names.size()) {
         oss << std::setw(14) << _par_names[i];
      } else {
        end=true;
        break;
      }
    }
    oss << std::endl;
    oss << "      ";
    for (unsigned int i=j*entries_per_line; i < (j+1)*entries_per_line; i++) {
      if (i<_par_values.size()) {
         oss << std::setw(14) << _par_values[i];
      } else{
        break;
      }
    }
    oss << std::endl;
    j+=1;
  }
  return oss.str();
}
