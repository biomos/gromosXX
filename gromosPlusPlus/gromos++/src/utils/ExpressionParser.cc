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

// ExpressionParser.cc
// This is not a good solution!
#ifndef INCLUDED_UTILS_EXPRESSIONPARSER
#include "ExpressionParser.h"
#endif

#include "Value.h"
#include "parse.h"
#include "../gcore/System.h"
#include "../bound/Boundary.h"
#include "../gromos/Exception.h"

#include <cmath>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace utils
{

  //////////////////////////////////////////////////

  template<typename T>
  ExpressionParser<T>::ExpressionParser()
  {
    init_op();
  }
  
  template<typename T>
  ExpressionParser<T>::ExpressionParser(gcore::System * sys,
					bound::Boundary * pbc)
    : d_value_traits(sys, pbc)
  {
    init_op();
  }

  template<typename T>
  void ExpressionParser<T>::init_op()
  {
    d_op_string = "*/+-,=!><&|?:";

    d_op["sin"] = op_sin;
    d_op["cos"] = op_cos;
    d_op["tan"] = op_tan;
    d_op["asin"] = op_asin;
    d_op["acos"] = op_acos;
    d_op["atan"] = op_atan;
    d_op["exp"] = op_exp;
    d_op["ln"] = op_ln;
    d_op["abs"] = op_abs;
    d_op["abs2"] = op_abs2;
    d_op["dot"] = op_dot;
    d_op["cross"] = op_cross;
    d_op["ni"] = op_ni;
    d_op["sqrt"] = op_sqrt;

    d_op["*"] = op_mul;
    d_op["/"] = op_div;
    d_op["+"] = op_add;
    d_op["-"] = op_sub;
    
    d_op["!"] = op_not;
    d_op["=="] = op_eq;
    d_op["!="] = op_neq;
    d_op["<"] = op_less;
    d_op[">"] = op_greater;
    d_op["<="] = op_lesseq;
    d_op[">="] = op_greatereq;
    d_op["&&"] = op_and;
    d_op["||"] = op_or;

    d_op["?"] = op_condask;
    d_op[":"] = op_condition;
    
    d_op[","] = op_comma;
    
  }
  
  template<typename T>
  T ExpressionParser<T>::parse_expression(std::string s,
					  std::map<std::string, T> & var)
  {
    std::string::size_type it = 0;
    std::vector<expr_struct> expr;

    operation_enum op = op_undef;
    do{
      op = _parse_token(op, s, it, var, expr);
    }while (op != op_undef);
    
    try{
      T res = calculate(expr, var);
      return res;
    }
    catch(std::string n){
      throw gromos::Exception("expression parser", "Variable '" + n + "' unknown!");
    }
  }

  template<typename T>
  void ExpressionParser<T>::parse_expression(std::string s,
					     std::map<std::string, T> & var,
					     std::vector<expr_struct> & expr)
  {
    std::string::size_type it = 0;

    operation_enum op = op_undef;
    do{
      op = _parse_token(op, s, it, var, expr);
    }while (op != op_undef);
    // _parse_token(op_undef, s, it, var, expr);
  }

  template<typename T>
  void ExpressionParser<T>::parse_expression(std::string s,
					     std::map<std::string, T> & var,
					     std::map<std::string, std::vector<expr_struct> > & expr)
  {
    std::string::size_type it = 0;
    it = s.find('=');
    
    if (it == std::string::npos)
      throw gromos::Exception("expression parser", "expected '=' not found in expression");
    
    std::string namestr = s.substr(0, it-1);
    std::istringstream is(namestr);
    std::string name;
    is >> name;
    ++it;
    
    // std::cerr << "name = " << name << std::endl;

    std::vector<expr_struct> single_expr;

    operation_enum op = op_undef;
    do{
      op = _parse_token(op, s, it, var, expr);
    }while (op != op_undef);

    // _parse_token(op_undef, s, it, var, single_expr);
    
    expr[name] = single_expr;
  }

  template<typename T>
  void ExpressionParser<T>::parse_expression(std::vector<std::string> s,
					     std::map<std::string, T> & var,
					     std::map<std::string, std::vector<expr_struct> > & expr)
  {
    for(int i=0; i<s.size(); ++i){
      parse_expression(s[i], var, expr);
    }
  }

  template<typename T>
  void ExpressionParser<T>::print_expression(std::vector<expr_struct> & expr,
					     std::ostream & os)
  {
    for(unsigned int i=0; i < expr.size(); ++i){
      switch(expr[i].type){
	case expr_value:
	  os << "value     \t" << expr[i].value << "\n";
	  break;
	case expr_operator:
	  os << "operator  \t" << expr[i].op << "\n";
	  break;
	case expr_function:
	  os << "function  \t" << expr[i].op << "\n";
	  break;
	case expr_variable:
	  os << "variable  \t" << expr[i].name << "\n";
	  break;
	default:
	  os << "error\n";
      }
    }
  }

  template<typename T>
  operation_enum ExpressionParser<T>::_parse_token(operation_enum op,
					 std::string s,
					 std::string::size_type & it,
					 std::map<std::string, T> & var,
					 std::vector<expr_struct> & expr)
  {
    // function, unary operator, value
    _parse_value(s, it, var, expr);
    
    operation_enum op2 = op_undef;
    if (it != std::string::npos)
      op2 = _parse_operator(s, it);

    // std::cerr << "considering order op1=" << op << " op2=" << op2 << std::endl;
    
    if (op <= op2){
      return op2;
    }
    else{

      do{
	if (op2 == op_undef) break;
	
	operation_enum op3 = _parse_token(op2, s, it, var, expr);
	_commit_operator(op2, expr);
	
	op2 = op3;
	// std::cerr << "considering order op=" << op << " op2=" << op2 << std::endl;

      } while(op2 < op);
      
      // std::cerr << "and backtracking" << std::endl;
      
      return op2;
    }
  }

  template<typename T>
  bool ExpressionParser<T>::_parse_function
  (
   std::string s, 
   std::string::size_type & it,
   std::map<std::string, T> & var,
   std::vector<expr_struct> & expr
   )
  {
    // std::cerr << "parsing function " 
    // << s.substr(it, std::string::npos) << std::endl;
    
    std::string::size_type bra = s.find_first_not_of(" ", it);
    if (bra == std::string::npos) return false;
    
    operation_enum fop = op_undef;
        bool found = false;
    
    // just a bracket
    if (s[bra] == '('){
      bra += 1;
      found = true;
    }
    else{
      // or a 'real' function
      std::map<std::string,operation_enum>::const_iterator
	it = d_op.begin(),
	to = d_op.end();
      
      for( ; it!=to; ++it){
	std::string::size_type len = it->first.length() + 1;
	if (s.substr(bra, len) == it->first + "("){
	  bra += len;
	  fop = it->second;
	  found = true;
	  break;
	}
      }
    }
    if (!found){
      return false;
    }
    
    std::string::size_type ket = find_matching_bracket(s, '(', bra);
    if (ket == std::string::npos){
      throw gromos::Exception("expression parser", "could not find matching bracket");
    }

    std::string::size_type fit = 0;

    operation_enum op = op_undef;
    do{
      op = _parse_token(op, s.substr(bra, ket-bra-1), fit, var, expr);
    }while (op != op_undef);

    if (fop != op_undef){
      _commit_operator(fop, expr);
    }
    
    if (s.length() > ket)
      it = ket;
    else it = std::string::npos;
    
    return true;
  }
  
  template<typename T>
  operation_enum
  ExpressionParser<T>::_parse_unary_operator
  (
   std::string s, 
   std::string::size_type & it
   )
  {
    // std::cerr << "parsing unary: " << s.substr(it, std::string::npos)
    // << std::endl;

    std::string::size_type bra = s.find_first_not_of(" ", it);
    if (bra == std::string::npos) return op_undef;
  
    // std::cerr << "\tunary: " << s[bra] << std::endl;
  
    if (s[bra] == '-'){ it = bra + 1; return op_umin; }
    if (s[bra] == '+'){ it = bra + 1; return op_uplus; }   

    return op_undef;
  }

  template<typename T>
  operation_enum
  ExpressionParser<T>::_parse_operator
  (
   std::string s, 
   std::string::size_type & it
   )
  {
    // std::cerr << "parse_operator: " << s.substr(it, std::string::npos)
    // << std::endl;

    std::string::size_type bra = s.find_first_of(d_op_string, it);

    if (bra == std::string::npos) return op_undef;
    
    // check if its a two character operator
    int op_len = 1;
    std::string::size_type bra2 = s.find_first_of(d_op_string, bra+1);
    if (bra2 == bra+1) op_len = 2;

    std::string ops = s.substr(bra, op_len);
    // std::cerr << "ops: " << ops << std::endl;
    
    if (d_op.count(ops) == 0)
      throw gromos::Exception("expression parser", "operator undefined! (" + ops + ")");

    it = bra + op_len;
    return d_op[ops];
  }

  template<typename T>
  void ExpressionParser<T>::_parse_value
  (
   std::string s, 
   std::string::size_type & it,
   std::map<std::string, T> & var,
   std::vector<expr_struct> & expr
   )
  {
    // std::cerr << "parse value: " << s.substr(it, std::string::npos) << std::endl;
    
    if (!_parse_function(s, it, var, expr)){

      // std::cerr << "no, then try unary operator" << std::endl;
      
      operation_enum u_op = op_undef;
      u_op = _parse_unary_operator(s, it);
      
      if (u_op != op_undef){
	_parse_value(s, it, var, expr);
	_commit_operator(u_op, expr);
	return;
      }

      _commit_value(s, it, var, expr);
    }
  }
  
  template<typename T>
  void ExpressionParser<T>::_commit_value
  (
   std::string s, 
   std::string::size_type & it,
   std::map<std::string, T> & var,
   std::vector<expr_struct> & expr
   )
  {
    std::string::size_type vit = find_par(s, d_op_string, it, "(", ")");

    // std::cerr << "_commit_value: " << s.substr(it, vit-it) 
    // << " from original " << s << std::endl;
    
    try{
      expr_struct e(d_value_traits.parseValue(s.substr(it, vit - it), var));
      expr.push_back(e);
    }
    catch(gromos::Exception e){
      std::istringstream is(s.substr(it, vit-it));
      std::string name;
      is >> name;

      if (name == "")
	throw gromos::Exception("expression parser", "argument missing!");
      
      expr_struct ex(name);
      expr.push_back(ex);
    }
    
    it = vit;
    if (it >=  s.length()) it = std::string::npos;
  }

  template<typename T>
  void ExpressionParser<T>::_commit_operator
  (
   operation_enum op,
   std::vector<expr_struct> & expr
   )
  {
    // std::cerr << "committing operator " << op << std::endl;
    expr.push_back(expr_struct(op));
  }

  template<typename T>
  T ExpressionParser<T>::calculate
  (
   std::vector<expr_struct> & expr,
   std::map<std::string, T> & var
   )
  {
    for(unsigned int i=0; i<expr.size(); ++i){
      switch(expr[i].type){
	case expr_value: d_stack.push(expr[i].value); break;
	case expr_variable: {
	  if (var.count(expr[i].name) == 0)
	    throw gromos::Exception("expression parser", "bad token: " + expr[i].name);
	  d_stack.push(var[expr[i].name]);
	  break;
	}
	case expr_function: do_function(expr[i].op); break;
	case expr_operator: do_operation(expr[i].op); break;
	default: throw gromos::Exception("expression parser", "wrong expression type");
      }
    }
    
    if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
    T res = d_stack.top();
    d_stack.pop();
    return res;
  }

  template<typename T>
  void ExpressionParser<T>::calculate
  (
   std::string name,
   std::map<std::string, std::vector<expr_struct> > & expr,
   std::map<std::string, T> & var
   )
  {
    if (var.count(name) == 0){
      
      while(true){

	try{
	  T res = calculate(expr[name], var);
	  var[name] = res;
	  break;
	}
	catch(std::string n){
	  if (n == name)
	    throw gromos::Exception("expression parser", "Implicit expression for '" + n + "'");
	  
	  if (expr.count(n) == 0)
	    throw gromos::Exception("expression parser", "No expression to calculate '" + n + "'");
	  
	  calculate(n, expr, var);
	}
      }
    }
  }


  template<typename T>
  void ExpressionParser<T>::calculate
  (
   std::map<std::string, std::vector<expr_struct> > & expr,
   std::map<std::string, T> & var
   )
  {
    typename std::map<std::string, std::vector<expr_struct> >::const_iterator
      it = expr.begin(),
      to = expr.end();
    
    for( ; it != to; ++it){
      
	calculate(it->first, expr, var);
    }
  }

  template<typename T>
  void ExpressionParser<T>::clear
  (
   std::map<std::string, std::vector<expr_struct> > & expr,
   std::map<std::string, T> & var
   )
  {
    typename std::map<std::string, std::vector<expr_struct> >::const_iterator
      it = expr.begin(),
      to = expr.end();
    
    for( ; it != to; ++it){
      
      if (var.count(it->first)){
	// std::cerr << "\tclearing " << it->first << " = " << var[it->first] << std::endl;
	var.erase(it->first);
      }
      
    }
  }


  template<typename T>
  T ExpressionParser<T>::do_operation(operation_enum op)
  {
    T res;
    if (op < op_unary) throw gromos::Exception("expression parser", "operator is function");
    if (op < op_binary){
      if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
      T arg = d_stack.top();
      d_stack.pop();
      switch(op){
	case op_uplus: res = arg; break;
	case op_umin: res = -arg; break;
	default: throw gromos::Exception("expression parser", "unknown unary operator");
      }
    }
    else if (op < op_logical){
      if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
      T arg2 = d_stack.top();
      d_stack.pop();
      if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
      T arg1 = d_stack.top();
      d_stack.pop();
      switch(op){
	case op_add: res = arg1 + arg2; break;
	case op_sub: res = arg1 - arg2; break;
	case op_mul: res = arg1 * arg2; break;
	case op_div: res = arg1 / arg2; break;
	default: throw gromos::Exception("expression parser", "unknown binary operator");
      }
    }
    else if (op == op_comma){
      res = d_stack.top();
      d_stack.pop();
    }
    else if (op < op_undef){
      // delegate the logical operators, not everything might
      // be defined with those...
      ValueTraits<T>::do_operation(op, *this);
      return d_stack.top();
    }
    else{
      throw gromos::Exception("expression parser", "unknown / undef operator");
    }
    
    d_stack.push(res);
    return res;
  }

  template<typename T>
  void ExpressionParser<T>::do_function(operation_enum op)
  {
    ValueTraits<T>::do_function(op, *this);
  }

  template<typename T>
  void ExpressionParser<T>::do_logical_operation(operation_enum op)
  {
    if (op == op_undef) return;
    
    if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");

    if (op == op_not){
      T arg = d_stack.top();
      T res = ! arg;
      op = op_undef;
      d_stack.pop();
      d_stack.push(res);
      return;
    }

    if (op < op_ternary){
      T arg2 = d_stack.top();
      d_stack.pop();
      if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
      T arg1 = d_stack.top();
      d_stack.pop();
      
      T res;
      switch(op){
	case op_eq:        res = (arg1 == arg2); break;
	case op_neq:       res = (arg1 != arg2); break;
	case op_less:      res = (arg1 < arg2); break;
	case op_greater:   res = (arg1 > arg2); break;
	case op_lesseq:    res = (arg1 <= arg2); break;
	case op_greatereq: res = (arg1 >= arg2); break;
	case op_and:       res = (arg1 && arg2); break;
	case op_or:        res = (arg1 || arg2); break;
	default: return;
      }

      op = op_undef;
      d_stack.push(res);
      return;
    }

    if (op == op_condask){
      // do nothing! leave the arguments for the op_condition
      return;
    }
    
    if (op == op_condition){
      T arg3 = d_stack.top();
      d_stack.pop();
      if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
      T arg2 = d_stack.top();
      d_stack.pop();
      if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
      T arg1 = d_stack.top();
      d_stack.pop();
      
      T res = (arg1 == true) ? arg2 : arg3;
      op = op_undef;
      d_stack.push(res);
      return;
    }
  }
  
  template<typename T>
  void ExpressionParser<T>::do_general_function(operation_enum op)
  {
    if (op == op_undef) return;
    
    T res;
    if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
    T arg = d_stack.top();
    
    switch(op){
      case op_abs: res = T(abs(arg)); break;
      case op_sqrt: res = T(sqrt(arg)); break;
      default: return;
    }
    
    op = op_undef;
    d_stack.pop();
    d_stack.push(res);
  }
  
  template<typename T>
  void ExpressionParser<T>::do_trigonometric_function(operation_enum op)
  {
    if (op == op_undef) return;
    
    T res;
    if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
    T arg = d_stack.top();

    switch(op){
      case op_sin: res = sin(arg); break;
      case op_cos: res = cos(arg); break;
      case op_tan: res = tan(arg); break;
      case op_asin: res = asin(arg); break;
      case op_acos: res = acos(arg); break;
      case op_atan: res = atan(arg); break;
	
      default: return;
    }
    
    op = op_undef;
    d_stack.pop();
    d_stack.push(res);
  }

  template<typename T>
  void ExpressionParser<T>::do_transcendent_function(operation_enum op)
  {
    if (op == op_undef) return;
    
    T res;
    if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
    T arg = d_stack.top();
    
    switch(op){
      case op_exp: res = exp(arg); break;
      case op_ln:  res = log(arg); break;
      default: return;
    }
    
    op = op_undef;
    d_stack.pop();
    d_stack.push(res);
  }

  template<typename T>
  void ExpressionParser<T>::do_vector_function(operation_enum op)
  {
    if (op == op_undef) return;
    
    T res;
    
    switch(op){
      case op_abs2:
	{
	  if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
	  T arg = d_stack.top();
	  d_stack.pop();
	  res = abs2(arg);
	  break;
	}
      case op_dot:
	{
	  if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
	  T arg1 = d_stack.top();
	  d_stack.pop();
	  if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
	  T arg2 = d_stack.top();
	  d_stack.pop();
	  res = dot(arg1, arg2);
	  break;
	}
      case op_cross:
	{
	  if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
	  T arg1 = d_stack.top();
	  d_stack.pop();
	  if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
	  T arg2 = d_stack.top();
	  d_stack.pop();
	  res = cross(arg1, arg2);
	  break;
	}
      case op_ni:
	{
	  if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
	  T arg1 = d_stack.top();
	  d_stack.pop();
	  if (d_stack.empty()) throw gromos::Exception("expression parser", "too few arguments on stack");
	  T arg2 = d_stack.top();
	  d_stack.pop();
	  res = Value(d_value_traits.pbc()->nearestImage
		      (arg1.vec(), arg2.vec(), d_value_traits.sys()->box()));
	  break;
	}
      default: return;
    }

    op = op_undef;
    d_stack.push(res);
  }
}
