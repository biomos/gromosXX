/**
 * @file blockinput.h
 * read in blocks.
 */

#ifndef INCLUDED_BLOCKINPUT_H
#define INCLUDED_BLOCKINPUT_H

namespace io {

  /*
   * The function io::getline provides an override of 
   * std::getline in that it retrieves the next line (separated 
   * by const char& sep) from the input stream, which, after 
   * having been stripped of comments (indicated by 
   * const char& comm) and leading and trailing whitespace, is not empty.
   */
  std::istream& getline(
			       std::istream& is, 
			       std::string& s, 
			       const char& sep = '\n',
			       const char& comm = '#'
			       );

  /*
   * The function io::getblock retrieves the next block  
   * (separated by const string& sep) using io::getline.
   *
   * Finally, the vector<string> it writes to is resized to the 
   * number of strings read.
   */
  bool getblock(
				std::istream& is, 
				std::vector<std::string>& b, 
				const std::string& sep = "END"
				);

  /*
   * The io::concatenate utility allows the concatenation
   * of entries in a vector<string> into just one string. The
   * resulting entries are separated by const char& sep 
   * (defaulting to a newline character).
   */
  std::string& concatenate(
			   std::vector<std::string>::const_iterator begin,
			   std::vector<std::string>::const_iterator end,
			   std::string& s,
			   const char& sep = '\n'
			   );
  // MP: overload, why do we need to give s anyways if it is always first
  // cleared AND returned
  std::string concatenate(
                           std::vector<std::string>::const_iterator begin,
                           std::vector<std::string>::const_iterator end,
                           const char& sep = '\n'
                           );


  /**
   * do replacements in a string
   * @param str the string
   * @param search the search string
   * @param replace the replace string
   * @return the string with "search" being replaced with "replace"
   */
  std::string replace_string(std::string str, const std::string &search,
          const std::string &replace);


  /**
   * trim away leading empty lines (only whitespace)
   */
  void
  trimblock(std::vector<std::string> &block);
  
  /**
   * split string
   * @param s the string
   * @param delim delimiter
   * @returns vector of substrings
   */
  std::vector<std::string> split(const std::string &s, std::string delim);

  bool is_valid_int(const char* x, bool sign=true);

  template <typename T>
  inline std::string to_string ( T Number )
  {
     std::ostringstream ss;
     ss << Number;
     return ss.str();
  }
  

  template< class T>
  bool is_valid_type(T & var, const char* x) {
    std::istringstream is(x);
    T tmp_dbl;
    is >> tmp_dbl;
    if (is.fail()) return false;
    else return true;
  }

  template <>
  bool is_valid_type <int> (int & var, const char* x);

  template<>
  bool is_valid_type<unsigned int>(unsigned int & var, const char* x) ;

  /**
   * true if 0 or 1
   */
  template <>
  bool is_valid_type <bool> (bool & var, const char* x);

  // trim from start (in place)
   static inline void ltrim(std::string &s) {
     s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
   }

  // trim from end (in place)
   static inline void rtrim(std::string &s) {
     s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
   }

  // trim from both ends (in place)
   static inline void trim(std::string &s) {
     ltrim(s);
     rtrim(s);
   }
   
   // trim from start (copying)
   static inline std::string ltrimmed(std::string s) {
     ltrim(s);
     return s;
   }

   // trim from end (copying)
   static inline std::string rtrimmed(std::string s) {
     rtrim(s);
     return s;
   }

   // trim from both ends (copying)
   static inline std::string trimmed(std::string s) {
     trim(s);
     return s;
   }

  /**
   * simple logical expression evaluator
   * @param var the value to evaluate
   * @param expr terms composed of one of <=, >=, <, > and a number, several 
   * of them can be connected by || or &&, && has precedence, no brackets!
   */
  template < class T >
  bool evaluate_comparison(const T &var, std::string expr) {
    bool result = 0;
    
    std::size_t found = expr.find("<=");
    if (found != std::string::npos) {
      std::vector<std::string> parts = io::split(expr, "<=");
      try {
        T value;
        std::stringstream stream(parts[1]);
        stream >> value;
        result = var <= value;
        return result;
      } catch (...) {
        io::messages.add("evaluate_comparison: could not evaluate "+expr+"!",
		     "In_Parameter", io::message::error);
        return false;        
      }
    }
    
    found = expr.find(">=");
    if (found != std::string::npos) {
      std::vector<std::string> parts = io::split(expr, ">=");
      try {
        T value;
        std::stringstream stream(parts[1]);
        stream >> value;
        result = var >= value;
        return result;
      } catch (...) {
        io::messages.add("evaluate_comparison: could not evaluate "+expr+"!",
		     "In_Parameter", io::message::error);
        return false;       
      }        
    }
    
    found = expr.find("<");
    if (found != std::string::npos) {
      std::vector<std::string> parts = io::split(expr, "<");
      try {
        T value;
        std::stringstream stream(parts[1]);
        stream >> value;
        result = var < value;
        return result;
      } catch (...) {
        io::messages.add("evaluate_comparison: could not evaluate "+expr+"!",
		     "In_Parameter", io::message::error);
        return false;       
      }        
    }
    
    found = expr.find(">");
    if (found != std::string::npos) {
      std::vector<std::string> parts = io::split(expr, ">");
      try {
        T value;
        std::stringstream stream(parts[1]);
        stream >> value;
        result = var > value;
        return result;
      } catch (...) {
        io::messages.add("evaluate_comparison: could not evaluate "+expr+"!",
		     "In_Parameter", io::message::error);
        return false;        
      }       
    }
    if (expr != "") {
        io::messages.add("evaluate_comparison: could not evaluate "+expr+"!",
		     "In_Parameter", io::message::error);
        return false;        
    } else {
        return true;
    }
  }
 
  template < class T >
  bool evaluate_logical_expression(const T &var, const std::string &expr) {
    bool result = false;
    std::vector<std::string> substrings_or = io::split(expr, "||");
    for (unsigned int i=0; i<substrings_or.size(); i++) {
        bool iresult=true;
        std::vector<std::string> substrings_and = io::split(io::trimmed(substrings_or[i]), "&&");
        for (unsigned int j=0; j<substrings_and.size(); j++) {
            iresult = iresult && evaluate_comparison(var, io::trimmed(substrings_and[j]));
        }
        result = result || iresult;
    }
    return result;
  }
  
  
  /**
   * @class Block
   * reads parameters from a GROMOS block in form of a vector of strings
   */
  class Block {
  private:
    /*
     * default constructor
     */
    Block();
  
  public:


    /*
     * Constructors
     */
    Block(std::string name, std::string exampleblock="") : _blockname(name), _block_error(0), _numlines(0), _exampleblock(exampleblock)  {} 
    Block(std::string name, std::vector<std::string> &buffer, std::string exampleblock="") : _blockname(name), _block_error(0), _exampleblock(exampleblock) { read_buffer(buffer); } 

    /**
     * Destructor.
     */
    ~Block() {}
    
    
    /**
     * check if there is something in the block
     * read block into istringstream
     * @param buffer vector of strings, the first of which is the block name, the last the "END"
     */
   int read_buffer(std::vector<std::string> &buffer, bool required=false);
   
   /**
    * get current content of _block_parameters
    * so it can be added to error messages
    */
    std::string get_readparameters();
   
   /**
    * read values and check their types and ranges
    * @param var the variable
    * @param expr e.g. >=0 && < 1
    * @param allowed comma-separated string of allowed values
    */
    template <class T > 
    int get_next_parameter(std::string varname, T &var, std::string expr="", std::string allowed="", bool allow_missing=false);

   /**
    * print the read parameters if there were errors
    * else check if there are leftover parameters and warn
    */
    void get_final_messages(bool print_read=true);
    
    /**
     * set the block example string
     */
    void set_exampleblock(std::stringstream exampleblock) { _exampleblock=exampleblock.str(); };
    
    std::string name() { return _blockname; }
    bool error() { return _block_error; }
    int numlines() { return _numlines; }
    void reset_error() { _block_error = 0; }
    
    /**
     * get the block example
     */
    std::string exampleblock() { return _exampleblock; }
    
   private:
     std::istringstream _lineStream;
     std::vector<std::string> _par_names, _par_values;
     std::string _blockname;
     int _block_error; 
     int _numlines;   
     std::string _exampleblock;
  };
  
  /**
   * if unsuccessful, var is just not set!! this was also the previous behaviour
   */
  template < class T > 
  int io::Block::get_next_parameter(std::string  varname, T &var, std::string expr, std::string allowed, bool allow_missing) {
    _par_names.push_back(varname);
    std::string tmp_string;
    _lineStream >> tmp_string;
  
    if (_lineStream.eof()) {
      if (!allow_missing) {
      io::messages.add(_blockname + " block reached END before "+varname
          +" could be read!",
		    "BlockInput", io::message::error);
      }
      _par_values.push_back(tmp_string);
      _block_error=1;
      return 1;  
    }

    // check if we have a valid data type
    if (io::is_valid_type(var, tmp_string.c_str())) {
      std::istringstream iss(tmp_string);
      iss >> var;
    } else {
      io::messages.add(_blockname + " block: wrong value type for "+varname+": "
                       +tmp_string,  //+" (required: "+typeid(var).name()+")"  NOTE: I could write also which it should be with typeid(var).name(), but this is a bit cryptic as the output is compiler/system? dependent - I get i and d for int and double
		     "In_Parameter", io::message::error);
      _par_values.push_back(tmp_string);
      _block_error=1;
      return 1;
    }

    _par_values.push_back(io::to_string(var));
    io::trim(expr);
    io::trim(allowed);
    bool ok=false;
    if (allowed == "" && expr == "") ok=true;

    std::ostringstream errormsg;
    errormsg << _blockname << " block: " <<varname<< " parameter is " << var << ", should be: ";
    if (expr != "") {
       ok = evaluate_logical_expression(var, expr);
       errormsg << " " << expr;
    }
  
    // even if the variable does not fit with expr, 
    // it is ok if it is in allowed
    if (allowed != "") {
        if (expr != "") errormsg << "\n\t  or\n";
        errormsg << " one of: ";
        std::vector<std::string> substrings=io::split(allowed, ",");
        for (unsigned int i=0; i<substrings.size(); i++) {
          T value;
          std::stringstream stream(substrings[i]);
          stream >> value;
          if (var == value) {
            ok=true;
          }
      }
      errormsg << allowed;
    }
  
    if (!ok) {
      _block_error=1;
      io::messages.add(errormsg.str(), "In_Parameter", io::message::error);
      return 1;
    }
    return 0;
  }

}

#endif

