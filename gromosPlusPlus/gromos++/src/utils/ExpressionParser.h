/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// ExpressionParser.h

#ifndef INCLUDED_UTILS_EXPRESSIONPARSER
#define INCLUDED_UTILS_EXPRESSIONPARSER

#include <stack>
#include <map>
#include <stdexcept>
#include <sstream>

#include "../bound/Boundary.h"
#include "Value.h"

namespace utils
{
  enum operation_enum{
    // functions
    // - basic
    op_abs = 10,
    op_sqrt = 11,
    // - trigonometric:
    op_sin = 20,
    op_cos = 21,
    op_tan = 22,
    op_asin = 23,
    op_acos = 24,
    op_atan = 25,
    // - transcendent:
    op_exp = 30,
    op_ln = 31,
    // - vector:
    op_abs2 = 51,
    op_dot = 52,
    op_cross = 53,
    op_ni = 54,
    // unary expressions
    op_unary = 100,
    op_umin = 101,
    op_uplus = 102,
    op_not = 103,
    // binary expressions
    op_binary = 200,
    op_mul = 201,
    op_div = 202,
    op_add = 203,
    op_sub = 204,
    // logical expressions
    op_logical = 300,
    op_eq = 301,
    op_neq = 302,
    op_less = 303,
    op_greater = 304,
    op_lesseq = 305,
    op_greatereq = 306,
    op_and = 400,
    op_or = 401,
    // ternary operator
    op_ternary = 500,
    op_condask = 501,
    op_condition = 502,
    // sequence operator
    op_comma   = 999,
    // no op
    op_undef = 1000
  };

  template<typename T>
  class ExpressionParser;

  template<typename T>
  class ValueTraits;
  
  template<>
  class ValueTraits<int>;
  
  template<>
  class ValueTraits<Value>;

  /**
   * @section ExpressionParser Expression Parser
   * parse expressions
   * @todo document more, since it is referenced in xray_map.cc
   **/
  template<typename T>
  class ExpressionParser
  {
  public:

    enum expr_enum{
      expr_value,
      expr_operator,
      expr_function,
      expr_variable
    };
    
    struct expr_struct
    {
      expr_struct(T value) 
	: type(expr_value), value(value), op(op_undef), name("") {}
      expr_struct(operation_enum op) 
	: type(expr_operator), value(T()), op(op), name("") 
      {
	if (op < op_unary) type = expr_function;
      }
      expr_struct(std::string s) 
	: type(expr_variable), value(T()), op(op_undef), name(s) {}
      

      expr_enum type;
      T value;
      operation_enum op;
      std::string name;
      
      std::string toString()
      {
	std::ostringstream os;
	
	switch(type){
	  case expr_value: os <<    "value: " << value; break;
	  case expr_operator: os << "op:    " << op; break;
	  case expr_function: os << "f:     " << op; break;
	  case expr_variable: os << "var:   " << name; break;
	  default: os << "-- unknown --";
	}
	return os.str();
      }
    };

    /**
     * Constructor
     */
    ExpressionParser();

    /**
     * Constructor
     * (for types that need a System and a Boundary)
     */
    ExpressionParser(gcore::System * sys,
		     bound::Boundary * pbc);
    
    /**
     * parse an expression and calculate the result
     * @param s expression
     * @param var map of parameters and known variables
     */
    T parse_expression(std::string s,
		       std::map<std::string, T> & var);

    /**
     * parse an expression and store it in the
     * UPN stack expr
     * @param s expression string
     * @param var map of known variables and parameters
     * (the parameters are needed as some might be used by special types; see value traits)
     * @param expr returns resulting expression
     */
    void parse_expression(std::string s,
			  std::map<std::string, T> & var,
			  std::vector<expr_struct> & expr);

    /**
     * parse a named expression s (name = expr) and store it in the
     * map of UPN expressions expr
     * @param s expression string
     * @param var map of known variables and parameters
     * (the parameters are needed as some might be used by special types; see value traits)
     * @param expr map of expressions
     */
    void parse_expression(std::string s,
			  std::map<std::string, T> & var,
			  std::map<std::string, std::vector<expr_struct> > & expr);


    /**
     * parse a vector of named expressions (name = expr) and store it in the
     * map of UPN expressions expr
     * @param s expression string
     * @param var map of known variables and parameters
     * (the parameters are needed as some might be used by special types; see value traits)
     * @param expr map of expressions
     */
    void parse_expression(std::vector<std::string> s,
			  std::map<std::string, T> & var,
			  std::map<std::string, std::vector<expr_struct> > & expr);

    /**
     * get the result from the expression expr
     * using parameters and known variables from
     * the map var.
     * @param expr expression
     * @param var  variables and parameters
     */
    T calculate(std::vector<expr_struct> & expr,
		std::map<std::string, T> & var);

    /**
     * get the result of a named expression from
     * a given map of expressions
     * @param name expression name (variable)
     * @param var parameters and known variables
     * @param expr map of expressions (var = expr)
     */
    void calculate(std::string name,
		   std::map<std::string, std::vector<expr_struct> > & expr,
		   std::map<std::string, T> & var);

    /**
     * get results from a map of expressions using
     * a map of parameters and known variables
     * @param expr expressions
     * @param var parameters and known variables
     */
    void calculate(std::map<std::string, std::vector<expr_struct> > & expr,
		   std::map<std::string, T> & var);

    /**
     * clear the variables associated with expressions from
     * the parameter and known variables map, causing
     * recalculation.
     * @param expr variables associated with expressions
     * @param var parameters and known variables
     */
    void clear(std::map<std::string, std::vector<expr_struct> > & expr,
	       std::map<std::string, T> & var);

    /**
     * print the UPN (reverse polnish notation) stack
     * of an expression
     */
    void print_expression(std::vector<expr_struct> & expr,
			  std::ostream & os = std::cout);

    /**
     * general functions, called from value traits that support them
     */
    void do_general_function(operation_enum op);
    /**
     * trigonomeric functions, called from value traits that support them
     */
    void do_trigonometric_function(operation_enum op);
    /**
     * transcendent functions, called from value traits that support them
     */
    void do_transcendent_function(operation_enum op);
    /**
     * vector functions, called from value traits that support them
     */
    void do_vector_function(operation_enum op);
    /**
     * comparisons and logical operators, called from value traits that support them
     */
    void do_logical_operation(operation_enum op);

  private:
    /**
     * init operations
     * (map operation tokens to the op_enum)
     */
    void init_op();
    
    /**
     * parse an expression
     * @param arg1 current argument
     * @param op1 current operator
     * @param rest rest of the expression
     * @param var parameters and known variables
     */
    T _parse_expression(T arg1, operation_enum op1,
			std::string rest,
			std::map<std::string, T> & var);

    /**
     * parse a single token
     * @param s expression
     * @param it current position in expression
     * @param var parameters and known variables
     * @param expr resulting expression UPN stack
     */
    operation_enum _parse_token(operation_enum op,
		      std::string s,
		      std::string::size_type & it,
		      std::map<std::string, T> & var,
		      std::vector<expr_struct> & expr);

    /**
     * parse a function
     * @param s expression
     * @param it current position in expression
     * @param var parameters and known variables
     * @param expr resulting expression UPN stack
     */
    bool _parse_function(std::string s,
			 std::string::size_type & it,
			 std::map<std::string, T> & var,
			 std::vector<expr_struct> & expr);

    /**
     * parse a value
     * @param s expression
     * @param it current position in expression
     * @param var parameters and known variables
     * @param expr resulting expression UPN stack
     */
    void _parse_value(std::string s,
		     std::string::size_type & it,
		     std::map<std::string, T> & var,
		     std::vector<expr_struct> & expr);
    
    /**
     * commit a value to the resulting expression
     * @param s expression
     * @param it current position in expression
     * @param var parameters and known variables
     * @param expr resulting expression UPN stack
     */
    void _commit_value(std::string s,
		      std::string::size_type & it,
		      std::map<std::string, T> & var,
		      std::vector<expr_struct> & expr);
    
    /**
     * commit an operator to the resulting expression
     * @param expr resulting expression UPN stack
     */
    void _commit_operator(operation_enum op,
			  std::vector<expr_struct> & expr);

    /**
     * parse an unary operator
     * @param s expression
     * @param it current position in expression
     */
    operation_enum _parse_unary_operator(std::string s,
					 std::string::size_type & it);

    /**
     * parse an operator
     * @param s expression
     * @param it current position in expression
     */
    operation_enum _parse_operator(std::string s,
				   std::string::size_type & it);

    /**
     * execute an operation (during calculation)
     * using the current stack (internal)
     */
    T do_operation(operation_enum op);
    
    /**
     * exectute a function (during calculation)
     * using the current stack (internal)
     */
    void do_function(operation_enum op);

    /**
     * known operators and their string representation
     */
    std::map<std::string, operation_enum> d_op;
    /**
     * characters recognized as operators
     */
    std::string d_op_string;
    /**
     * internal stack (used during calculation)
     */
    std::stack<T> d_stack;
    /**
     * Value trait. Specifies which operations are
     * allowed for a given value type T.
     */
    ValueTraits<T> d_value_traits;
  };

  
  //////////////////////////////////////////////////
  // traits class for values
  //////////////////////////////////////////////////
  template<typename T>
  class ValueTraits
  {
  public:
    ValueTraits(gcore::System * sys,
		bound::Boundary *pbc) : d_sys(sys), d_pbc(pbc) {}
    
    ValueTraits() : d_sys(NULL), d_pbc(NULL) {}

    T parseValue(std::string s, std::map<std::string, T> & var)
    {
      std::istringstream is(s);
      T t;
      if (!(is >> t))
	throw gromos::Exception("expression parser", "could not read value (" + s + ")");
      return t;
    }
    
    static void do_function(operation_enum op,
			    ExpressionParser<T> & ep)
    {
      ep.do_general_function(op);
      ep.do_trigonometric_function(op);
      ep.do_transcendent_function(op);
      // ep.do_vector_function(op);
    }

    static void do_operation(operation_enum op,
			     ExpressionParser<T> & ep)
    {
      ep.do_logical_operation(op);
    }
    

  private:
    gcore::System * d_sys;
    bound::Boundary * d_pbc;
  };

  // integer specialisation
  template<>
  class ValueTraits<int>
  {
  public:
    ValueTraits() {}

    int parseValue(std::string s, std::map<std::string, int> & var)
    {
      std::istringstream is(s);
      int t;
      if (!(is >> t))
	throw gromos::Exception("expression parser", "could not read value (" + s + ")");
      return t;
    }
    
    static void do_function(operation_enum op,
			    ExpressionParser<int> & ep)
    {
      ep.do_general_function(op);
    }

    static void do_operation(operation_enum op,
			     ExpressionParser<int> & ep)
    {
      ep.do_logical_operation(op);
    }

  };

  // and finally the one for Values!
  template<>
  class ValueTraits<Value>
  {
  public:
    ValueTraits(gcore::System *sys,
		bound::Boundary *pbc) : d_sys(sys), d_pbc(pbc) {}
    
    // ValueTraits() : d_sys(NULL), d_pbc(NULL) {}

    Value parseValue(std::string s, std::map<std::string, Value> & var)
    {
      Value v;
      v.parse(s, var, *d_sys, d_pbc);
      return v;
    }
    
    static void do_function(operation_enum op,
			    ExpressionParser<Value> & ep)
    {
      ep.do_general_function(op);
      ep.do_trigonometric_function(op);
      ep.do_transcendent_function(op);
      ep.do_vector_function(op);
    }

    static void do_operation(operation_enum op,
			     ExpressionParser<Value> & ep)
    {
      ep.do_logical_operation(op);
    }

    bound::Boundary * pbc()
    {
      if (d_pbc == NULL) throw gromos::Exception("expression parser", "pbc is NULL");
      return d_pbc;
    }

    gcore::System * sys()
    {
      if (d_sys == NULL) throw gromos::Exception("expression parser", "sys is NULL");
      return d_sys;
    }

  private:
    gcore::System * d_sys;
    bound::Boundary * d_pbc;
  };

}

#include "ExpressionParser.cc"

#endif
