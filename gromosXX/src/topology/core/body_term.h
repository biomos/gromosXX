/**
 * @file body_term.h
 * define the structures to hold topology information.
 */

#ifndef INCLUDED_BODY_TERM_H
#define INCLUDED_BODY_TERM_H

namespace topology
{
  /**
   * @struct one_body_term_struct
   * one body terms (topological information).
   */
  struct one_body_term_struct
  {
    /**
     * Constructor.
     * @param i atom i.
     * @param t interaction type.
     */
    one_body_term_struct(size_t i, size_t t) : i(i), type(t) {};
    
    /**
     * atom i.
     */
    size_t i;
    /**
     * interaction type.
     */
    size_t type;
    /**
     * equal operator
     */
    bool operator==(one_body_term_struct const & b){
      return (b.i==i && b.type==type);
    }
    
  };

  /**
   * @struct perturbed_one_body_term_struct
   * perturbed one body terms (topological information).
   */
  struct perturbed_one_body_term_struct
  {
    /**
     * Constructor.
     * @param i atom i.
     * @param t_A interaction type.
     * @param t_B interaction type.
     */
    perturbed_one_body_term_struct(size_t i, size_t t_A, size_t t_B) 
      : i(i), A_type(t_A), B_type(t_B) {};
    
    /**
     * atom i.
     */
    size_t i;
    /**
     * interaction type for state A.
     */
    size_t A_type;

    /**
     * interaction type for state B.
     */
    size_t B_type;
    /**
     * equal operator
     */
    bool operator==(perturbed_one_body_term_struct const & b){
      return (b.i==i && b.A_type==A_type && b.B_type==B_type);
    }

  };
  
  /**
   * @struct two_body_term_struct
   * two body term interaction topological information.
   */
  struct two_body_term_struct : public one_body_term_struct
  {
    /**
     * Constructor.
     * @param i atom i.
     * @param j atom j.
     * @param t interaction type.
     */
    two_body_term_struct(size_t i, size_t j, size_t t) 
      : one_body_term_struct(i, t), j(j) {};
    
    /**
     * atom j.
     */
    size_t j;
    /**
     * equal operator
     */
    bool operator==(two_body_term_struct const & b){
      return (one_body_term_struct::operator==(b) && b.j==j);
    }

  };

  /**
   * @struct perturbed_two_body_term_struct
   * perturbed two body term interaction topological information.
   */
  struct perturbed_two_body_term_struct : public perturbed_one_body_term_struct
  {
    /**
     * Constructor.
     * @param i atom i.
     * @param j atom j.
     * @param t_A interaction type for state A.
     * @param t_B interaction type for state B.
     */
    perturbed_two_body_term_struct(size_t i, size_t j, size_t t_A, size_t t_B) 
      : perturbed_one_body_term_struct(i, t_A, t_B), j(j) {};
    
    /**
     * atom j.
     */
    size_t j;
    /**
     * equal operator
     */
    bool operator==(perturbed_two_body_term_struct const & b){
      return (perturbed_one_body_term_struct::operator==(b) && b.j==j);
    }
  };

  /**
   * @struct three_body_term_struct
   * three body term interaction topological information.
   */
  struct three_body_term_struct : public two_body_term_struct
  {
    /**
     * Constructor.
     * @param i atom i.
     * @param j atom j.
     * @param k atom k.
     * @param t interaction type.
     */
    three_body_term_struct(size_t i, size_t j, size_t k, size_t t) 
      : two_body_term_struct(i, j, t), k(k) {};
    
    /**
     * atom k.
     */
    size_t k;
    /**
     * equal operator
     */
    bool operator==(three_body_term_struct const & b){
      return (two_body_term_struct::operator==(b) && b.k==k);
    }

  };

  /**
   * @struct perturbed_three_body_term_struct
   * perturbed three body term interaction topological information.
   */
  struct perturbed_three_body_term_struct : public perturbed_two_body_term_struct
  {
    /**
     * Constructor.
     * @param i atom i.
     * @param j atom j.
     * @param k atom k.
     * @param t_A interaction type for state A.
     * @param t_B interaction type for state B.
     */
    perturbed_three_body_term_struct(size_t i, size_t j, size_t k, 
				     size_t t_A, size_t t_B) 
      : perturbed_two_body_term_struct(i, j, t_A, t_B), k(k) {};
    
    /**
     * atom k.
     */
    size_t k;
    /**
     * equal operator
     */
    bool operator==(perturbed_three_body_term_struct const & b){
      return (perturbed_two_body_term_struct::operator==(b) && b.k==k);
    }
  };

  /**
   * @struct four_body_term_struct
   * four body term interaction topological information.
   */
  struct four_body_term_struct : public three_body_term_struct
  {
    /**
     * Constructor.
     * @param i atom i.
     * @param j atom j.
     * @param k atom k.
     * @param l atom l.
     * @param t interaction type.
     */
    four_body_term_struct(size_t i, size_t j, size_t k, size_t l, size_t t) 
      : three_body_term_struct(i, j, k, t), l(l) {};
    
    /**
     * atom l.
     */
    size_t l;
    /**
     * equal operator
     */
    bool operator==(four_body_term_struct const & b){
      return (three_body_term_struct::operator==(b) && b.l==l);
    }

  };

  /**
   * @struct perturbed_four_body_term_struct
   * perturbed four body term interaction topological information.
   */
  struct perturbed_four_body_term_struct : public perturbed_three_body_term_struct
  {
    /**
     * Constructor.
     * @param i atom i.
     * @param j atom j.
     * @param k atom k.
     * @param l atom l.
     * @param t_A interaction type for state A.
     * @param t_B interaction type for state B.
     */
    perturbed_four_body_term_struct(size_t i, size_t j, size_t k, size_t l, 
				    size_t t_A, size_t t_B) 
      : perturbed_three_body_term_struct(i, j, k, t_A, t_B), l(l) {};
    
    /**
     * atom l.
     */
    size_t l;
    /**
     * equal operator
     */
    bool operator==(perturbed_four_body_term_struct const & b){
      return (perturbed_three_body_term_struct::operator==(b) && b.l==l);
    }
  };

  struct position_restraint_struct
  {
    /**
     * Constructor.
     */
    position_restraint_struct(size_t seq, math::Vec pos, 
			       double bfactor = 1.0)
      : seq(seq), pos(pos), bfactor(bfactor)
    {}
    
    /**
     * sequence number.
     */
    size_t seq;
    /**
     * position.
     */
    math::Vec pos;
    /**
     * atomic b-factor.
     */
    double bfactor;
  };
    
}

#endif
