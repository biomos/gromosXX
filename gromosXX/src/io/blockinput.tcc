// included by BlockInput.h

inline 
istream& 
io::getline(
  istream& is, 
  string& s, 
  const char& sep,
  const char& comm
  ) {

  unsigned short int ii;

  while (is.good()) {
    std::getline(is, s, sep);
    ii = find(s.begin(), s.end(), comm) - s.begin();
    if (ii == s.size()) break; // no comment
    else if (!ii) continue;    // comment on first position
    else s.erase(s.begin() + ii, s.end());
  }
  
  return is;
}

inline 
istream& 
io::getblock(
  istream& is, 
  vector<string>& b, 
  const string& sep
  ) {

  if (!b.size())
    b.push_back("");
  vector<string>::iterator dest = b.begin();

  while (1) {

    if (dest == b.end()) {
      b.push_back("");
      dest = b.end() - 1;
    }       

    getline(is, *dest);

    if (*dest == sep)
      break;

    if (!is.good()) 
      throw runtime_error("error reading block.");

    dest++;
  }

  b.erase(dest, b.end());

  return is;
}

inline
string& 
io::concatenate(
  vector<string>::iterator begin,
  vector<string>::iterator end,
  string& s,
  const char& sep = '\n'
) {

  s.clear();
  while (begin != end) {
    s += *begin;
    s += sep;
    begin++;
  }

  return s;
}
