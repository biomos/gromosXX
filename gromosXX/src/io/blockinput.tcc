// included by BlockInput.h

inline 
std::istream& 
io::getline(
	    std::istream& is, 
	    std::string& s, 
	    const char& sep,
	    const char& comm
	    )
{
  unsigned short int ii;

  while (is.good()) {
    std::getline(is, s, sep);
    ii = std::find(s.begin(), s.end(), comm) - s.begin();
    if (ii == s.size()) break; // no comment
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

    if (*dest == sep)
      break;

    if (!is.good()) 
      throw std::runtime_error("error reading block.");

    dest++;
  }

  // b.erase(dest, b.end());

  return is;
}

inline
std::string& 
io::concatenate(
		std::vector<std::string>::iterator begin,
		std::vector<std::string>::iterator end,
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
