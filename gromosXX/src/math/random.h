/**
 * @file random.h
 * generate random numbers using GSL or the GROMOS96 algorithm
 */

#ifndef INCLUDED_RANDOM_H
#define INCLUDED_RANDOM_H

#include <gsl/gsl_rng.h>

namespace math {
  /** 
   * @class RandomGenerator
   * interface to a generic random number generator
   */
  class RandomGenerator {
    public :
    /**
     * default constructor
     */
    RandomGenerator() : gauss_m(0.0), gauss_s(1.0) {}
    /** 
     * default destructor
     */  
    virtual ~RandomGenerator() {};

    /**
     * accessor to the data needed to restart the random number generator
     * in very trivial generators, this is the seed.
     */
    virtual std::string seed() = 0;
    /**
     * sets the state of the random number genertor to the one of 
     * definied by the seed string.
     *
     * if it fails, a std::runtime_error is thrown.
     */
    virtual void seed(std::string s) = 0;

    /**
     * samples a random number between [0;1]
     */
    virtual double get() = 0;
    /**
     * samples a gaussian distributed random number
     */
    virtual double get_gauss() = 0;
 
    /**
     * gets the mean of the gaussian distribution
     */ 
    double & mean() { return gauss_m; }
    /**
     * sets the mean of the gaussian distribution
     */
    void mean(double m) { gauss_m = m; }
    /**
     * gets the standard deviation of the gaussian distribution
     */
    double & stddev() { return gauss_s; }
    /**
     * sets the standard deviation of the gaussian distribution
     */
    void stddev(double s) { gauss_s = s;}

    protected:
    /**
     * mean of the gaussian distribution
     */
    double gauss_m;
    /**
     * standard deviation of the gaussian distribution
     */
    double gauss_s;
  };

  /**
   * @class RandomGeneratorG96
   * a GROMOS96 style random number generator
   */
  class RandomGeneratorG96 : public RandomGenerator {
    public:
    /**
     * create an instance with zero seed.
     */
    RandomGeneratorG96() { 
      seed(0);
    }
    /**
     * crate an instance with given seed.
     */
    RandomGeneratorG96(std::string s) { 
      seed(s);
    }
    /**
     * destroy instance
     */
    ~RandomGeneratorG96() {}

    std::string seed();
    /**
     * sets the seed. The seed must contain an integer value or 
     * std::runtime_error is thrown.
     */
    void seed(std::string s);

    /**
     * use the GROMOS96 algorithm to get a random number between [0;1]
     * This uses a pretty default UNIX random number generator
     *
     * IRAND = (IRAND * B + 1) % A
     *
     * with pretty non default values for B (31415821) and A (100000000)
     */
    double get();
    /**
      * gets a gaussian distributed random number create by the Box-Mueller
      * method
      */
    double get_gauss();
 
    protected:
    /**
     * internal seed
     */
    int ig;
  };

  /**
   * @class RandomGeneratorGSL
   * GSL random number generator
   *
   * This random number genertor can be tuned using environment variables and
   * has the big advantage that the random numbers are far more sophisticated
   * than in the G96 algorithm.
   *
   * See: http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html
   */
  class RandomGeneratorGSL : public RandomGenerator {
    public:

    /**
     * creates an instance
     */
    RandomGeneratorGSL();
    /**
     * creates an instance unsing the given seed
     */
    RandomGeneratorGSL(std::string s);

    /**
     * destroys the instance
     */
    ~RandomGeneratorGSL();

    /**
     * gets the internal state of the random number generator (binary data)
     * as a base64 encoded string. There is no possibility to express the
     * internal state as a simple integer like in simpler algorithms.
     */
    std::string seed();
    /**
     * sets the internal state of the random number generator. If a base64
     * encoded string is given, the state is set accordingly. If an integer
     * number is given, the latter is used as seed.
     *
     * If everything goes wrong, a std::runtime_error is thrown.
     */
    void seed(std::string s);
 
    /**
     * samples a uniformly distributed random number between [0;1]
     */
    double get();
    /**
     * samples a gaussian distributed random number using the Box-Muller
     * method (GSL implementation)
     */
    double get_gauss();

    protected:
    /*
     * creates the GSL random number generator
     */
    void create_instance();
    /*
     * the GSL random number generator
     */
    gsl_rng *rng;
  };
}

#endif

