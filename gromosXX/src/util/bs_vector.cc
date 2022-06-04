#include "../stdheader.h"

//#include "../math/periodicity.h"
#include "../io/message.h"

#include "../util/template_split.h"
#include "../util/debug.h"
#include "../util/bs_vector.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE bs_leus

/**************************************************
 * BS_VECTOR
 */
double util::BS_Vector::abs2(){
  double length2 = 0;
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++){
    length2 += *it * *it;
  }
  return length2;
}

double util::BS_Vector::normalize(){
  double length = sqrt(this->abs2());
  if (length > 0){
    double inverseLength = 1 / length;
    this->scale(inverseLength);
  }
  return length;
}

void util::BS_Vector::nullify(){
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++){
    *it = 0;
  }
}

void util::BS_Vector::scale(const double scalar){
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++){
    DEBUG(10, "Value * scalar =  " << *it << " * " << scalar);
    *it *= scalar;
    DEBUG(10, "\t = " << *it)
  }
}

void util::BS_Vector::minus(const BS_Vector& subtrahend, BS_Vector& result){
  assert(this->size() == subtrahend.size());
  if (this->size() != subtrahend.size()){
    io::messages.add("Two BS_Vectors with different sizes (minus)!", "BS_Vector",
            io::message::critical);
    DEBUG(5, "Two BS_Vectors with different sizes (minus)!")
    return;
  }
  result.clear();
  // minuend
  BS_Vector::iterator m_i = this->begin(),
          m_end = this->end();
  BS_Vector::const_iterator s_i = subtrahend.begin();
  
  double diff = 0.0;
  for (; m_i != m_end; m_i++, s_i++){
    diff = *m_i - *s_i;
    DEBUG(12, "BS_Vector.minus: Value: " << diff);
    result.push_back(diff);
  }
}

util::BS_Vector util::BS_Vector::operator *(const double scalar){
  BS_Vector result;
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++){
    result.push_back(*it * scalar);
  }
  return result;
}

util::BS_Vector util::BS_Vector::operator +(const BS_Vector& summand) {
  assert(this->size() == summand.size());
  if (this->size() != summand.size()){
    io::messages.add("Two BS_Vectors with different sizes (plus)!", "BS_Vector",
            io::message::critical);
    DEBUG(5, "Two BS_Vectors with different sizes (plus)!")
  }
  BS_Vector result;
  BS_Vector::const_iterator s_i = summand.begin();
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++, s_i++){
    result.push_back(*it + *s_i);
  }
  return result;
}
util::BS_Vector util::BS_Vector::operator -(const BS_Vector &subtrahend){
  assert(this->size() == subtrahend.size());
  BS_Vector result;
  for (unsigned int i = 0; i < this->size(); i++){
    result.push_back((*this)[i] - subtrahend[i]);
  }
  return result;
}

void util::BS_Vector::operator +=(const BS_Vector& summand) {
  assert(this->size() == summand.size());
  if (this->size() != summand.size()){
    io::messages.add("Two BS_Vectors with different sizes (+=)!", "BS_Vector",
            io::message::critical);
    DEBUG(5, "Two BS_Vectors with different sizes (+=)! (" << this->size() << " != " << summand.size() << ")");
  }
  BS_Vector::const_iterator s_i = summand.begin();
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++, s_i++){
    *it += *s_i;
  }
}

double util::BS_Vector::dot(const BS_Vector &other){
  assert(this->size() == other.size());
  if (this->size() != other.size()){
    io::messages.add("Two BS_Vectors with different sizes (dot)!", "BS_Vector",
            io::message::critical);
    DEBUG(5, "Two BS_Vectors with different sizes (dot)!")
  }
  double dotProduct = 0;
  BS_Vector::const_iterator other_i = other.begin();
  BS_Vector::iterator it = this->begin(),
                      to = this->end();
  for (; it != to; it++, other_i++){
    dotProduct += *it * *other_i;
  }
  return dotProduct;
}

void util::BS_Vector::create(std::vector<double>& values)
{
  this->clear();
  std::vector<double>::iterator val_i = values.begin(),
          val_end = values.end();
  for (; val_i != val_end; val_i++){
    this->push_back(*val_i);
  }
}

std::string util::BS_Vector::str(){
  std::ostringstream os;
  os << "[ ";
  BS_Vector::const_iterator vec_i = this->begin(),
          vec_end = this->end();
  
  if (vec_i != vec_end){
    os << *vec_i;
    vec_i++;
  }
  else {
    os << "'empty'";
  }
  for (; vec_i != vec_end; vec_i++){
    os << ", " << *vec_i;
  }
  os << " ]";
  return os.str();
}
