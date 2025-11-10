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

// gmath_Expression

#ifndef INCLUDED_GMATH_EXPRESSION
#define INCLUDED_GMATH_EXPRESSION

#include <string>
#include <vector>

namespace gmath
{
  /**
   * Class Expression
   * A class to do pocket calculator expressions
   * 
   * This class allows you to calculate expressions that have been
   * specified using a string. It knows how to use basic operators
   * +, -, * and / and the order in which to apply them. Or at least that
   * * and / come before + and -. Brackets and the functions cos, sin, log 
   * and exp are also implemented
   *
   * @class Expression
   * @author B.C. Oostenbrink
   * @ingroup gmath
   */
  class Expression
  {
    std::vector<double> d_val;
    std::vector<std::string> d_op;
    int d_new;
    double d_result;

  public:
    /**
     * Expression constructor
     * @param s The string contains an expression consisting of
     *                 numbers, the tokens *, /, + and -, brackets ( and ), 
     *                 operators sin, cos, log and exp
     *                 and variables, indicated by a1, a2, .. an
     *                 before calculation, these need to be set by a call to
     *                 setValue giving a vector<double> of length n.
     */
    Expression(std::string s);
    /**
     * Expression constructor
     * @param s The string contains an expression consisting of
     *                 numbers, the tokens *, /, + and -, brackets ( and ),
     *                 operators sin, cos, log and exp
     *                 and variables, indicated by a1, a2, .. an
     * @param v This vector should be of length n and contain
     *                 the variables that are needed in the Expression
     */
    Expression(std::string s, std::vector<double>& v);
    /**
     * Method to re-define the expression. Be carefull that the number of 
     * variables that are required might change.
     * @param s The string contains an expression consisting of
     *                 numbers, the tokens *, /, + and -, brackets ( and ),
     *                 operators sin, cos and, log exp
     *                 and variables, indicated by a1, a2, .. an
     *                 before calculation, these need to be set by a call to
     *                 setValue giving a vector<double> of length n.
     */
    void setExpression(std::string s);
    /**
     * Method to set the variable that are needed to evaluate the expression
     * @param v A vector of doubles that should contain at 
     *                 least the number of elements of the highest variable
     *                 that is specified in the expression.
     */
    void setValues(std::vector<double>& v);
    /**
     * Method that evaluates the expression with the current set of variables
     */
    double value();
    /**
     * Method that writes the expression to an ostream. Good for debugging.
     */
    void writeExpression(std::ostream& os);
    /**
     * Method that writes the expression to an ostream, replacing variable
     * names with their current value.
     */
    void writeExpressionValue(std::ostream& os);
    
  private:
    void Tokenize(const std::string& str,
		  std::vector<std::string>& tokens,
		  const std::string& delimiters);
    double calc(std::vector<std::string>& op, std::vector<double>& val, 
		int f, int t);
    bool allowed_token(std::string s);
    bool is_function(std::string s, double on, double& res);
    bool is_function(std::string s);
    bool is_operator(std::string s);
    bool is_variable(std::string s);
    bool is_bracket(std::string s);
    bool is_number(std::string s);
    void find_bracket(std::vector<std::string>&op, int &first, int &last);
    void special_functions(std::vector<std::string>& op,
			   std::vector<double>& val,
			   int f, int& t);
  };
}
#endif
