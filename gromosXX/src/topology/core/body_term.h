/**
 * @file body_term.h
 * define the structures to hold topology information
 * of the n-body terms (like bonds, angles, dihedrals).
 */

#ifndef INCLUDED_BODY_TERM_H
#define INCLUDED_BODY_TERM_H

namespace topology
{
  /**
   * @enum functional_form
   * functional form
   */
  enum functional_form{
    /**
     * full harmonic (attractive and repulsive)
     */
    harmonic = 1,
    /**
     * half harmonic, attractive
     */
    attractive = 2,
    /**
     * half harmonic, repulsive
     */
    repulsive = 3
  };
  
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
    one_body_term_struct(unsigned int i, unsigned int t) : i(i), type(t) {};
    
    /**
     * atom i.
     */
    unsigned int i;
    /**
     * interaction type.
     */
    unsigned int type;
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
    perturbed_one_body_term_struct(unsigned int i, unsigned int t_A, unsigned int t_B) 
      : i(i), A_type(t_A), B_type(t_B) {};
    
    /**
     * atom i.
     */
    unsigned int i;
    /**
     * interaction type for state A.
     */
    unsigned int A_type;

    /**
     * interaction type for state B.
     */
    unsigned int B_type;
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
    two_body_term_struct(unsigned int i, unsigned int j, unsigned int t) 
      : one_body_term_struct(i, t), j(j) {};
    
    /**
     * atom j.
     */
    unsigned int j;
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
    perturbed_two_body_term_struct(unsigned int i, unsigned int j, unsigned int t_A, unsigned int t_B) 
      : perturbed_one_body_term_struct(i, t_A, t_B), j(j) {};
    
    /**
     * atom j.
     */
    unsigned int j;
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
    three_body_term_struct(unsigned int i, unsigned int j, unsigned int k, unsigned int t) 
      : two_body_term_struct(i, j, t), k(k) {};
    
    /**
     * atom k.
     */
    unsigned int k;
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
    perturbed_three_body_term_struct(unsigned int i, unsigned int j, unsigned int k, 
				     unsigned int t_A, unsigned int t_B) 
      : perturbed_two_body_term_struct(i, j, t_A, t_B), k(k) {};
    
    /**
     * atom k.
     */
    unsigned int k;
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
    four_body_term_struct(unsigned int i, unsigned int j, unsigned int k, unsigned int l, unsigned int t) 
      : three_body_term_struct(i, j, k, t), l(l) {};
    
    /**
     * atom l.
     */
    unsigned int l;
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
    perturbed_four_body_term_struct(unsigned int i, unsigned int j, unsigned int k, unsigned int l, 
				    unsigned int t_A, unsigned int t_B) 
      : perturbed_three_body_term_struct(i, j, k, t_A, t_B), l(l) {};
    
    /**
     * atom l.
     */
    unsigned int l;
    /**
     * equal operator
     */
    bool operator==(perturbed_four_body_term_struct const & b){
      return (perturbed_three_body_term_struct::operator==(b) && b.l==l);
    }
  };

  /**
   * Position restraints information.
   */
  struct position_restraint_struct
  {
    /**
     * Constructor.
     */
    position_restraint_struct(unsigned int seq, math::Vec pos, 
			       double bfactor = 1.0)
      : seq(seq), pos(pos), bfactor(bfactor)
    {}
    
    /**
     * sequence number.
     */
    unsigned int seq;
    /**
     * position.
     */
    math::Vec pos;
    /**
     * atomic b-factor.
     */
    double bfactor;
  };
    
  /**
   * J-Value restraints.
   */
  struct jvalue_restraint_struct
  {
    /**
     * Constructor.
     */
    jvalue_restraint_struct(int i, int j, int k, int l,
			    double K, double J0,
			    double a, double b, double c,
			    double delta,
			    functional_form H)
      : i(i), j(j), k(k), l(l),
	K(K), J0(J0),
	a(a), b(b), c(c), delta(delta),
	H(H)
    {
    }
    
    /**
     * atom sequence numbers.
     */
    unsigned int i, j, k, l;
    /**
     * force constant
     */
    double K;
    /**
     * J0
     */
    double J0;
    /**
     * Karplus parameter.
     */
    double a, b, c;
    /**
     * phase shift
     */
    double delta;
    /**
     * functional form.
     * half harmonic attractive,
     * half harmonic repulsive,
     * harmonic
     */
    functional_form H;
  };

}

#endif
