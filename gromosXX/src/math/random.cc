/**
 * @file random.cc
 * implementation of random number generator
 */

#include <stdheader.h>
#include <gromosXX/math/random.h>
#include <gromosXX/util/coding.h>

#include <gsl/gsl_randist.h>

std::string math::RandomGeneratorG96::seed() {
  std::ostringstream buf; buf << ig;
  return buf.str();
}

void math::RandomGeneratorG96::seed(std::string s) {
  std::istringstream str(s);
  if (!(str >> ig))
    throw std::runtime_error("invalid seed");
}

double math::RandomGeneratorG96::get() {

  // this code was more or less copied from GROMOS96

  const int m    = 100000000;
  const int m1   = 10000;
  const int mult = 31425821;

  int irand = std::abs(ig) % m;

//  MULTIPLY IRAND BY MULT, BUT TAKE INTO ACCOUNT THAT OVERFLOW
//  MUST BE DISCARDED, AND DO NOT GENERATE AN ERROR.

  int irandh = int(irand / m1);
  int irandl = irand % m1;
  int multh  = int(mult / m1);
  int multl  = mult % m1;

  irand = ((irandh*multl+irandl*multh) % m1) * m1 + irandl+multl;
  irand = (irand + 1) % m;

//  CONVERT IRAND TO A REAL RANDOM NUMBER BETWEEN 0 AND 1.

  double r = int(irand / 10) * 10.0;
  r = r / m;
  if ((r <= 0.0) || (r > 1.0))
      r = 0.0;
  ig = irand;
  return r;
}

double math::RandomGeneratorG96::get_gauss() {
  // Box-Muller from GROMOS96
  double w1, w2, r;
  do {
    w1 = 2.0 * get() - 1.0;
    w2 = 2.0 * get() - 1.0;
 
    r = w1 * w1 + w2 * w2;
  } while( r > 1.0 || r == 0.0);

  return mean() + stddev() * w1 * sqrt(-2.0 * log(r) / r);
}

math::RandomGeneratorGSL::RandomGeneratorGSL() {
  create_instance();
  std::ostringstream buf; buf << gsl_rng_default_seed;
  seed(buf.str());
}

math::RandomGeneratorGSL::RandomGeneratorGSL(std::string s) {
  create_instance();
  seed(s);
}

math::RandomGeneratorGSL::~RandomGeneratorGSL() {
  // release the memory of the GSL rng
  gsl_rng_free(rng);
}

void math::RandomGeneratorGSL::create_instance() {
  // get a random number generator type
  const gsl_rng_type *rng_type;
  //enable control via environment variables
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;
  // get the rundom number generator
  rng = gsl_rng_alloc(rng_type);
}

void math::RandomGeneratorGSL::seed(std::string s) {
  std::istringstream str(s);
  unsigned long int init;

  // try to decode it to a string of binary data
  std::string decoded = util::base64_decode(s);
  if (decoded != "") {
    // copy the binary data to the gsl state
    memcpy(gsl_rng_state(rng), decoded.c_str(), decoded.length());
  } else if ((str >> init))  // try to read an integer
    gsl_rng_set(rng, init);
  else
    throw std::runtime_error("invalid seed");
}

std::string math::RandomGeneratorGSL::seed() {
  // get the GSL rng's state pointer and size
  void * state = gsl_rng_state (rng);
  size_t n = gsl_rng_size (rng);

  // encode with base64 to string
  return util::base64_encode(reinterpret_cast<const unsigned char*>(state), n);
}

double math::RandomGeneratorGSL::get() {
  // sample uniform
  return gsl_rng_uniform(rng);
}

double math::RandomGeneratorGSL::get_gauss() {
  // sample gaussian
  return mean() + gsl_ran_gaussian(rng, stddev());
}

