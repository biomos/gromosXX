

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <string>
#include <map>
#include <vector>

using namespace std;

vector<string> pairlist, nonbonded;
vector<string> makefile;
vector<string> header;
vector<string> code;

map<string, vector<pair<string, vector<string> > > > parameter;

void create_current(map<string, vector<pair<string, vector<string> > > >::const_iterator it,
		    map<string, vector<pair<string, vector<string> > > >::const_iterator const to,
		    map<string, string> &current,
		    vector<pair<pair<string, string>, pair<string, string> > > const &rule);

void write_solution(map<string, string> &current);
void write_header();
void write_code();
void add_to_code(string name, map<string, string> &current);
void write_makefile();

int main(int argc, char *argv[])
{
  if (argc != 2){
    cout << "usage:\n\t"
	 << argv[0]
	 << " template_file\n\n";
    return 1;
  }

  ifstream templ(argv[1]);

  if (!templ.is_open()){
    cerr << "could not open template file \"split.dat\"\n" << endl;
    return 1;
  }
  
  string s;
  
  templ >> s;
  if (s != "PAIRLIST"){
    cerr << "wrong file format! PAIRLIST block expected!\n" << endl;
    return 1;
  }
  
  templ >> s;
  while (s != "END"){
    pairlist.push_back(s);
    templ >> s;
  }

  templ >> s;
  if (s != "NONBONDED"){
    cerr << "wrong file format! NONBONDED block expected!\n" << endl;
    return 1;
  }
  
  templ >> s;
  while(s != "END"){
    nonbonded.push_back(s);
    templ >> s;
  }
  
  string name, value;
  vector<string> condition;
  vector<pair<string, vector<string> > > state;
  vector<pair<pair<string, string>, pair<string, string> > > rule;
  string r1, r2;

  while (!templ.eof()){
    templ >> s;
    if (templ.eof()) break;
    
    if (s == "PARAMETER"){
      templ >> name;
      
      templ >> s;
      while(s != "END"){
	condition.clear();
	while(s != "@"){
	  condition.push_back(s);
	  templ >> s;
	}
	
	templ >> s;
	state.push_back(pair<string, vector<string> >(s, condition));
	
	templ >> s;
      }
      parameter[name] = state;
      state.clear();
      
    }
    else if (s == "RULE"){
      
      while (true){
	templ >> r1;
	if (r1 == "END") break;
	templ >> r2;
	
	string::size_type sep1 = r1.find('@');
	string::size_type sep2 = r2.find('@');
	
	if (sep1 == string::npos || sep2 == string::npos){
	  cerr << "error in rule definition " << r1 << "\t" << r2 << endl;
	  return 1;
	}
	string r1a = r1.substr(0, sep1);
	string r1b = r1.substr(sep1+1, string::npos);
	string r2a = r2.substr(0, sep2);
	string r2b = r2.substr(sep2+1, string::npos);

	cout << "Rule: " << r1a << " @ " << r1b << "  &&  "
	     << r2a << " @ " << r2b << endl;
	
	rule.push_back(pair<pair<string, string>, pair<string, string> >
		       (pair<string, string>(r1a, r1b),
			pair<string, string>(r2a, r2b)));
      }
    }
    
  }

  // init the makefile vector
  makefile.push_back("split_files =");
  
  map<string, string> current;
  
  map<string, vector<pair<string, vector<string> > > >::const_iterator
    it = parameter.begin(),
    to = parameter.end();

  create_current(it, to, current, rule);

  write_header();
  write_code();
  write_makefile();
  
  return 0;
  
}

void create_current(map<string, vector<pair<string, vector<string> > > >::const_iterator it,
		    map<string, vector<pair<string, vector<string> > > >::const_iterator const to,
		    map<string, string> &current,
		    vector<pair<pair<string, string>, pair<string, string> > > const &rule)
{
  
  if (it == to){
    // we have a solution
    // check rules
    vector<pair<pair<string, string>, pair<string, string> > >::const_iterator
      r_it = rule.begin(),
      r_to = rule.end();
    
    for( ; r_it != r_to; ++r_it){
      if ((current[r_it->first.first] == r_it->first.second) &&
	  (current[r_it->second.first] == r_it->second.second))
	return;
    }

    // ok, rules fulfilled
    map<string, string>::const_iterator sol_it = current.begin(),
      sol_to = current.end();
    
    for( ; sol_it != sol_to; ++sol_it)
      cout << sol_it->first << " = " << sol_it->second << "\n";

    cout << "\n"; 

    write_solution(current);

    return;
  }
  
  vector<pair<string, vector<string> > >::const_iterator 
    v_it = it->second.begin(),
    v_to = it->second.end();
  
  map<string, vector<pair<string, vector<string> > > >::const_iterator  next = it;
  ++next;
  
  for( ; v_it != v_to; ++v_it){
    current[it->first] = v_it->first;

    create_current(next, to, current, rule);
  }
}

string & operator++(string & s)
{
  return s = s + "  ";
}

string & operator--(string & s)
{
  return s = s.substr(0, s.length() - 2);
}

void write_solution(map<string, string> &current)
{
  static int count = 0;
  ++count;

  string indent = "";
  
  stringstream ss, name;
  ss << "split/split_" << count << ".cc";
  name << "split_" << count;
  
  makefile.push_back(" \\\n\tnonbonded/" + ss.str());

  add_to_code("create_" + name.str(), current);

  string decl = "void create_" + name.str() + "(topology::Topology const & topo,\n"
    "\tconfiguration::Configuration const & conf,\n"
    "\tsimulation::Simulation const & sim,\n"
    "\tinteraction::Forcefield & ff,\n"
    "\tio::In_Topology & it)";
  
  header.push_back(decl + ";\n\n");

  ofstream splitf(ss.str().c_str());
  
  // header
  splitf << "// " << name.str() << "\n\n"
	 << "#include \"nonbonded_split.h\"\n\n"
	 << "using namespace interaction;\n"
	 << "using namespace math;\n"
	 << "using namespace std;\n\n"
	 << decl
	 << "\n{\n";
  
  // body
  ++indent;

  // some output of what we are doing
  splitf << indent << "// some output of what we are doing\n";
  
  map<string, string>::const_iterator sol_it = current.begin(),
    sol_to = current.end();
    
  for( ; sol_it != sol_to; ++sol_it){

    splitf << indent << "cout << setw(40) << left << \"" 
	   << sol_it->first.substr(2, string::npos)
	   << "\" << setw(30) << left << \""
	   << sol_it->second
	   << "\" << \"\\n\";\n";
  }
  splitf << "\n";
  
  // nonbonded interaction type
  splitf << indent << "typedef ";

  for(size_t i=0; i<nonbonded.size(); ++i){
    if (current.count(nonbonded[i]))
      splitf << current[nonbonded[i]];
    else{
      if (nonbonded[i] == ">"){
	--indent; 
	splitf << "\n" << indent; 
      }
      
      splitf << nonbonded[i];

      if (nonbonded[i] == "<"){
	++indent;
	splitf << "\n" << indent;
      }
      if (nonbonded[i] == ",") splitf << "\n" << indent;
    }
  }

  splitf << " nonbonded_interaction_type;\n\n";

  // pairlist type
  splitf << indent << "typedef ";

  for(size_t i=0; i<pairlist.size(); ++i){
    if (current.count(pairlist[i]))
      splitf << current[pairlist[i]];
    else{
      if (pairlist[i] == ">"){
	--indent; 
	splitf << "\n" << indent; 
      }
      
      splitf << pairlist[i];

      if (pairlist[i] == "<"){
	++indent;
	splitf << "\n" << indent;
      }
      if (pairlist[i] == ",") splitf << "\n" << indent;
    }
  }

  splitf << " pairlist_algorithm_type;\n";
  
  // create the pairlist
  splitf << indent << "pairlist_algorithm_type * pa = new pairlist_algorithm_type;\n\n"
	 << indent << "nonbonded_interaction_type * ni = new nonbonded_interaction_type(pa);\n\n"
	 << indent << "it.read_lj_parameter(ni->lj_parameter());\n"
	 << indent << "ni->initialize(topo, conf, sim);\n\n"
	 << indent << "ff.push_back(ni);\n\n";

  splitf << "}\n\n";

  splitf.close();
  
}

void write_header()
{
  ofstream hf("split_header.h");
  if (!hf.is_open()){
    cerr << "could not open split_header.h for writing" << endl;
    exit(1);
  }
  
  for(size_t i=0; i < header.size(); ++i){
    hf << header[i];
  }
  
  hf.close();

}

void add_to_code(string name, map<string, string> &current)
{

  string indent="  ";
  
  stringstream ss;
  
  ss << " if(";
  indent+="        ";
  
  map<string, string>::const_iterator
    it = current.begin(),
    to = current.end();

  for( ; it != to; ++it){
    if(it != current.begin())
      ss << " && \n" << indent;
    
    vector<pair<string, vector<string> > > & p = parameter[it->first];
    for(size_t i=0; i < p.size(); ++i){
      if (p[i].first == it->second){
	
	ss << "( ";
	for(size_t j=0; j<p[i].second.size(); j++){
	  ss << p[i].second[j] << " ";
	  
	}
	ss << ")";
      }
    }
  }
  ss << "){\n\n";

  --indent;
  --indent;
  --indent;
  
  ss << indent << name << "(topo, conf, sim, ff, it);\n\n";
  --indent;
  ss << indent << "}\n\n"
     << indent << "else";
  
  code.push_back(ss.str());
  
}

void write_code()
{
  ofstream cd("create_nonbonded.cc");
  if (!cd.is_open()){
    cerr << "could not open create_nonbonded.cc for writing" << endl;
    exit(1);
  }

  cd << "/**\n"
     << " * @file create_nonbonded.cc\n"
     << " * create the nonbonded interaction.\n"
     << " */\n\n"
     << "#include \"create_header.h\"\n"
     << "#include \"split_header.h\"\n\n"
     << "#undef MODULE\n"
     << "#undef SUBMODULE\n"
     << "#define MODULE interaction\n"
     << "#define SUBMODULE nonbonded\n\n"
     << "using namespace interaction;\n"
     << "using namespace math;\n\n";
  
  cd << "int interaction::create_g96_nonbonded(interaction::Forcefield & ff,\n"
     << "\t\ttopology::Topology const & topo,\n"
     << "\t\tsimulation::Simulation const & sim,\n"
     << "\t\tconfiguration::Configuration const & conf,\n"
     << "\t\tio::In_Topology & it,\n"
     << "\t\tbool quiet)\n"
     << "{\n\n ";
  
  
  string indent="  ";
  
  for(size_t i=0; i < code.size(); ++i){
    cd << code[i];
  }

  cd << indent << " {\n";
  ++indent;
  
  cd << indent << "io::messages.add(\"case not supported!\",\n"
     << indent << "                 \"create_nonbonded\",\n"
     << indent << "                 io::message::error);\n"
     << indent << "return 1;\n";
  --indent;
  cd << indent << "}\n";
  
  cd << indent << "return 0;\n";
  --indent;
  cd << indent << "}\n\n";
  
  cd.close();

}

void write_makefile()
{
  for(size_t i=0; i<makefile.size(); ++i)
    cout << makefile[i];
  cout << "\n\n";
}
