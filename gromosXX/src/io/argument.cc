/**
 * @file argument.cc
 * implementation of the Argument class.
 */

#include <stdheader.h>

#include "argument.h"
#include "blockinput.h"
#include <util/error.h>
#include <util/usage.h>

namespace io{

  typedef std::multimap<std::string,std::string>::value_type argType;

  // checks if an argument is known
  static unsigned int isKnown(const std::string str, std::set<std::string> known_args){
    return unsigned(known_args.count(str));
  }

  Argument::Argument()
    : std::multimap<std::string,std::string>(),
      d_usage(""),
      d_prog(""),
      d_known(),
      empty("")
  {
  }
  
  int Argument::parse(int argc, char **argv,
		      util::Known &known, bool empty_ok) 
  {

    if(argc) d_prog = argv[0];

    if (argc == 1 && !empty_ok) return E_USAGE;

    // for(int i=0; i<known.size(); ++i)
    // d_known.insert(std::string(known[i]));
    d_known = known;

    std::string s("");
  
    for(int i=1;i<argc;i++){
    
      if(std::string(argv[i])=="@f"){
	// input file
	++i;
	std::ifstream inp(argv[i]);
	if (parse_line(inp)){
	  inp.close();
	  return E_USAGE;
	}
	
	inp.close();
      }
      else
	s += std::string(argv[i])+' ';
    }

    std::istringstream is(s.c_str());
    return parse_line(is);
  }

  Argument::~Argument() {
  }

  /**
   * read in the arguments.
   */
  int Argument::parse_line(std::istream &istr)
  {
    // get away the comments
    std::string buff;
    std::string s("");
  
    while(istr.good()) {
      io::getline(istr, buff);

      s += buff;
      s += '\n';
    }

    std::istringstream is(s);

    std::string str, last;
  
    if(!(is>>last))
      return 0;
    if(last[0]!='@'){
      std::cerr << "argument does not begin with @ : " << last << std::endl;
      return E_USAGE;
    }
    
    last=last.substr(1);

    if(find(last)!=end())
      erase(lower_bound(last), upper_bound(last));

    if(!isKnown(last, d_known)) {
      std::cerr << "unknown argument : " << last << std::endl;
      return E_USAGE;
    }

    while(is>>str){
      if(str[0] == '@'){
	if(find(last) == end())
	  insert(argType(last,""));
	last=str.substr(1);
	
	if(find(last)!=end())
	  for(Argument::iterator l=lower_bound(last);
	      l!=upper_bound(last);++l)
	    erase(l);

	if(!isKnown(last, d_known)){
	  std::cerr << "unknown argument: " << last << std::endl;
	  return E_USAGE;
	}

	continue;
      }
      else
	insert(argType(last,str));
    }

    // insert the last one without value
    if(find(last)==end())
      insert(argType(last,""));
  
    return 0;
  }
  


  const std::string &Argument::operator[](const std::string &str)const
  {
    const_iterator f = find(str);
    if(f == end()){
      return empty;
    }
  
    else return find(str)->second;
  }

  int Argument::check(const std::string &str, int num_args)const
  {
    if(find(str) == end())
      return -1;
    int num=0;
    for(const_iterator l=lower_bound(str), u=upper_bound(str);
	l!=u;++l)
      if (l->second!="")++num;
    if(num<num_args)
      return -1;
    return 0;
  }

  int Argument::count(const std::string &str)const
  {
    if(find(str)==end())
      return -1;
    int num=0;
    for(const_iterator l=lower_bound(str), u=upper_bound(str);
	l!=u;++l)
      if (l->second!="")++num;
    return num;
  }

}
