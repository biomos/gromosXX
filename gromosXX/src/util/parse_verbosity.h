/**
 * @file parse_verbosity.h
 */

namespace util
{
  int parse_verbosity(io::Argument &args, std::string flag = "verb", 
		      std::ostream &os = std::cout);
}
