/**
 * @file random.cc
 * implementation of random number generator
 */

#include "../stdheader.h"
#include <cstring>
#include "../math/random.h"
#include "../util/coding.h"
#include "../simulation/multibath.h"
#include "../simulation/parameter.h"

#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE math
#define SUBMODULE math

#include <gsl/gsl_randist.h>

math::RandomGenerator * math::RandomGenerator::create(
const simulation::Parameter &param, std::string s) {
  switch(param.rng.rng) {
    case simulation::random_g96 :
      return new RandomGeneratorG96(s);
    case simulation::random_gsl :
      return new RandomGeneratorGSL(s, param.rng.gsl_rng);
    default :
      throw std::runtime_error("invalid random number algorithm");
  }
}

bool math::RandomGenerator::check(const simulation::Parameter &param) {
  try {
    math::RandomGenerator::create(param, "0");
    return true;
  } catch(std::runtime_error &e) {
    io::messages.add(e.what(), "RandomGenerator", io::message::error);  
  }
  return false;
}

std::string math::RandomGeneratorG96::seed() {
  std::ostringstream buf; buf << ig;
  ig = std::abs(ig) % m;
  stored = false;
  return buf.str();
}

void math::RandomGeneratorG96::seed(std::string s) {
  std::istringstream str(s);
  if (!(str >> ig))
    throw std::runtime_error("invalid seed");
}

double math::RandomGeneratorG96::get() {

  // this code was more or less copied from GROMOS96
  int irand = std::abs(ig) % m;

//  MULTIPLY IRAND BY MULT, BUT TAKE INTO ACCOUNT THAT OVERFLOW
//  MUST BE DISCARDED, AND DO NOT GENERATE AN ERROR.

  int irandh = irand / m1;
  int irandl = irand % m1;
  int multh  = mult / m1;
  int multl  = mult % m1;

  irand = ((irandh*multl+irandl*multh) % m1) * m1 + irandl*multl;
  irand = (irand + 1) % m;

//  CONVERT IRAND TO A REAL RANDOM NUMBER BETWEEN 0 AND 1.

  double r = int(irand / 10) * 10.0 / double(m);
  if ((r <= 0.0f) || (r > 1.0f))
      r = 0.0;
  ig = irand;
  return r;
}

double math::RandomGeneratorG96::get_gauss() {
  // Box-Muller from GROMOS96
  
  if (stored) { // then just return the stored value
    stored = false;
    return mean() + stddev() * stored_gaussian;
  }
  
  double w1 = 0.0, w2 = 0.0, r = 0.0;
  do {
    w1 = 2.0 * get() - 1.0;
    w2 = 2.0 * get() - 1.0;
 
    r = w1 * w1 + w2 * w2;
  } while( r > 1.0 || r == 0.0);

  // store the second gaussian
  stored = true;
  // only store gassian(0.0, 1.0) part. stddev and mean may change!
  stored_gaussian = w2 * sqrt(-2.0 * log(r) / r);
  
  return mean() + stddev() * w1 * sqrt(-2.0 * log(r) / r);
}

std::string math::RandomGeneratorG96::description() {
  return "GROMOS random number generator";
}


math::RandomGeneratorGSL::RandomGeneratorGSL(const int algorithm) {
  create_instance(algorithm);
  std::ostringstream buf; buf << gsl_rng_default_seed;
  seed(buf.str());
}

math::RandomGeneratorGSL::RandomGeneratorGSL(std::string s, const int algorithm) {
  create_instance(algorithm);
  seed(s);
}

math::RandomGeneratorGSL::~RandomGeneratorGSL() {
  // release the memory of the GSL rng
  gsl_rng_free(rng);
}

void math::RandomGeneratorGSL::create_instance(const int algorithm) {
  gsl_rng_env_setup();
  const gsl_rng_type **t = nullptr, **t0 = nullptr;
          
  t0 = gsl_rng_types_setup();
         
  bool hasDefault = false; 
  int i = 0;
  for(i = 0, t = t0; *t != 0; t++, i++) {
    std::string name((*t)->name);
    if (name == "mt19937" && algorithm == -1) {
      hasDefault = true;
      rng = gsl_rng_alloc(*t);
      return;
    }
    
    if (i == algorithm) {
      rng = gsl_rng_alloc(*t);
      return;     
    }
  }
  
  if (algorithm == -1 && !hasDefault) {
    throw std::runtime_error("GSL default algorithm (mt19937) not available");
  }
  
  throw std::runtime_error("GSL rng algorithm not found");
}

void math::RandomGeneratorGSL::seed(std::string s) {
  std::istringstream str(s);
  unsigned int init = 0;
  // try to decode it to a string of binary data
  std::string decoded = util::base64_decode(s);
  if (decoded != "") {
    DEBUG(10, "\trandom: setting state to: " << decoded);
    // copy the binary data to the gsl state
    memcpy(gsl_rng_state(rng), decoded.c_str(), decoded.length());
  } else if ((str >> init)) {  // try to read an integer
    DEBUG(10, "\trandom: setting seed to: " << init);
    gsl_rng_set(rng, init);
  } else {
    throw std::runtime_error("invalid seed");
  }
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

std::string math::RandomGeneratorGSL::description() {
  std::ostringstream desc;
  desc << "GSL random number generator (" << gsl_rng_name(rng) << ")";
  return desc.str();
}

