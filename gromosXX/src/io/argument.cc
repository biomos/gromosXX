/**
 * @file argument.cc
 * implementation of the Argument class.
 */

#include <algorithm>
#include <sstream>
#include <fstream>
#include <set>
#include <map>

#include "argument.h"

namespace io{

  typedef std::multimap<std::string,std::string>::value_type argType;

  // checks if an argument is known
  static int isKnown(const std::string str, std::set<std::string> known_args){
    /*
      std::set<std::string>::const_iterator fr, to=known_args.end();
      for(fr=known_args.begin();fr!=to; fr++)
      if(*fr==str)break;
      if(fr==to)return 0;
      return 1;
    */
    return known_args.count(str);
  }

  Argument::Argument(int argc, char **argv, int nknown, 
		       char **known, const std::string &usage)
    : std::multimap<std::string,std::string>(),
      d_usage(""),
      d_prog(""),
      d_known()
  {

    if(argc) d_prog = argv[0];
    d_usage="\n\n"+usage;

    for(int i=0;i<nknown;++i)
      d_known.insert(std::string(known[i]));

    std::string s("");
  
    for(int i=1;i<argc;i++){
    
      if(std::string(argv[i])=="@f"){
	// input file
	++i;
	std::ifstream inp(argv[i]);
	inp >> *this;
	inp.close();
      }
      else
	s += std::string(argv[i])+' ';
    }
  
    std::istringstream is(s.c_str());
    is >> *this;
  }

  Argument::~Argument() {
  }

  /**
   * read in the arguments.
   */
  std::istream &io::operator>>(std::istream &istr, Argument &args)
  {
    // get away the comments
    char buff[1000];
    std::string s("");
  
    while(istr.good()&&istr.getline(buff,1000)){
      s+=std::string(buff);
      if(s.find("#")<=s.length())s=s.erase(s.find("#"));
      s+='\n';
    }
    std::istringstream is(s.c_str());
  
    std::string str, last;
  
    if(!(is>>last))
      return istr;
    if(last[0]!='@')
      throw std::string(args.d_usage);
  
    last=last.substr(1);
    if(args.find(last)!=args.end())
      args.erase(args.lower_bound(last), args.upper_bound(last));

    if(!isKnown(last, args.d_known))
      throw std::string(args.d_usage);
  
    while(is>>str){
      if(str[0] == '@'){
	if(args.find(last) == args.end())
	  args.insert(argType(last,""));
	last=str.substr(1);
	if(args.find(last)!=args.end())
	  for(Argument::iterator l=args.lower_bound(last);
	      l!=args.upper_bound(last);++l)
	    args.erase(l);

	if(!isKnown(last, args.d_known))	
	  throw std::string(args.d_usage);
      
	continue;
      }
      else
	args.insert(argType(last,str));
    }

    // insert the last one without value
    if(args.find(last)==args.end())
      args.insert(argType(last,""));
  
    return istr;
  }
  


  const std::string &Argument::operator[](const std::string &str)const
  {
    const_iterator f = find(str);
    if(f == end()){
      throw std::string(d_usage);
    }
  
    else return find(str)->second;
  }

  int Argument::check(const std::string &str, int num_args)const
  {
    if(find(str) == end())
      throw std::string(d_usage);
    int num=0;
    for(const_iterator l=lower_bound(str), u=upper_bound(str);
	l!=u;++l)
      if (l->second!="")++num;
    if(num<num_args)
      throw std::string(d_usage);
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
