/**
 * @file blockinput.tcc
 * defines blockinput functions.
 */

inline 
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
    // ii = std::find(s.begin(), s.end(), comm) - s.begin();
    ii = s.find(comm, 0);

    // if (ii == s.size()) break; // no comment
    if (ii == std::string::npos) break; // no comment
    else if (!ii) continue;    // comment on first position
    else s.erase(s.begin() + ii, s.end());
  }
  
  return is;
}

inline 
std::istream& 
io::getblock(
	     std::istream& is, 
	     std::vector<std::string>& b, 
	     const std::string& sep
  )
{

  if (!b.size())
    b.push_back("");
  std::vector<std::string>::iterator dest = b.begin();

  while (1) {

    if (dest == b.end()) {
      b.push_back("");
      dest = b.end() - 1;
    }       
    
    getline(is, *dest);

    // if (*dest == sep)
    if (dest->find(sep) == 0)
      break;

    if (!is.good()) 
      throw std::runtime_error("error reading block: " + *(b.begin()));

    ++dest;
  }

  ++dest;
  b.erase(dest, b.end());

  return is;
}

inline
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

inline
void
io::trimblock(std::vector<std::string> &block)
{
  std::string s;
  std::vector<std::string>::iterator it = block.begin();
  while(true){
    std::istringstream is(*it);
    if (!(is >> s)){
      block.erase(it);
      it = block.begin();
      continue;
    }
    else break;
  }  
}
