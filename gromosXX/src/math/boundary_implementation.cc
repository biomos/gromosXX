/**
 * @file boundary_implementation.cc
 * implementation of the periodic boundary condition functions.
 */


#undef MODULE
#undef SUBMODULE
#define MODULE math

#ifdef WIN32
// Converts a floating point value to an integer, very fast.
inline int rint(float param)
{
	// Uses the FloatToInt functionality
	int a;
	int *int_pointer = &a;

	__asm  fld  param
	__asm  mov  edx,int_pointer
	__asm  FRNDINT
	__asm  fistp dword ptr [edx];

	return a;
}
#endif

/*
 * Constructor : vacuum
 */
inline math::Boundary_Implementation<math::vacuum>
::Boundary_Implementation(math::Box const & b)
  : m_box(b)
{
}

/**
 * Constructor : rectangular
 */
inline math::Boundary_Implementation<math::rectangular>
::Boundary_Implementation(math::Box const & b)
 : m_box(b)
{
  for(int i=0; i<3; ++i)
    m_half_box(i) = 0.5 * abs(m_box(i));
}

/**
 * Constructor : triclinic
 */
inline math::Boundary_Implementation<math::triclinic>
::Boundary_Implementation(math::Box const & b)
  : m_box(b)
{
  double volume = dot(cross(m_box(K), m_box(L)), m_box(M));
  
  assert(volume != 0);
  m_cross_K_L_M(0) = cross(m_box(L), m_box(M)) / -volume;
  m_cross_K_L_M(1) = cross(m_box(K), m_box(M)) / volume;
  m_cross_K_L_M(2) = cross(m_box(K), m_box(L)) / -volume;
}

////////////////////////////////////////////////////////////////////////////////
// box / shift vector accessors
////////////////////////////////////////////////////////////////////////////////

/**
 * const box accessor : vacuum
 */
inline math::Box const math::Boundary_Implementation<math::vacuum>::box()const
{
  return m_box;
}

/**
 * const box accessor : rectangular
 */
inline math::Box const & math::Boundary_Implementation<math::rectangular>::box()const
{
  return m_box;
}

/**
 * const box accessor : triclinic
 */
inline math::Box const & math::Boundary_Implementation<math::triclinic>::box()const
{
  return m_box;
}

/**
 * box element accessor (d1,d2) : vacuum
 */
inline double math::Boundary_Implementation<math::vacuum>
::box(unsigned int d1, unsigned int d2)const
{
  return m_box(d1)(d2);
}

/**
 * box element accessor (d1,d2) : rectangular
 */
inline double math::Boundary_Implementation<math::rectangular>
::box(unsigned int d1, unsigned int d2)const
{
  return m_box(d1)(d2);
}

/**
 * box element accessor (d1,d2) : triclinic
 */
inline double math::Boundary_Implementation<math::triclinic>
::box(unsigned int d1, unsigned int d2)const
{
  return m_box(d1)(d2);
}

/**
 * shift struct accessor : vacuum
 */
inline math::Boundary_Implementation<math::vacuum>::shift_struct &
math::Boundary_Implementation<math::vacuum>
::shift(unsigned int i)
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * const shift struct accessor : vacuum
 */
inline math::Boundary_Implementation<math::vacuum>::shift_struct const &
math::Boundary_Implementation<math::vacuum>
::shift(unsigned int i)const
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * shift struct accessor : rectangular
 */
inline math::Boundary_Implementation<math::rectangular>::shift_struct &
math::Boundary_Implementation<math::rectangular>
::shift(unsigned int i)
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * const shift struct accessor : rectangular
 */
inline math::Boundary_Implementation<math::rectangular>::shift_struct const &
math::Boundary_Implementation<math::rectangular>
::shift(unsigned int i)const
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * shift struct accessor : triclininc
 */
inline math::Boundary_Implementation<math::triclinic>::shift_struct &
math::Boundary_Implementation<math::triclinic>
::shift(unsigned int i)
{
  assert(27 > i);
  return m_shift[i];
}

/**
 * const shift struct accessor : triclinic
 */
inline math::Boundary_Implementation<math::triclinic>::shift_struct const &
math::Boundary_Implementation<math::triclinic>
::shift(unsigned int i)const
{
  assert(27 > i);
  return m_shift[i];
}

////////////////////////////////////////////////////////////////////////////////
// nearest image functions
////////////////////////////////////////////////////////////////////////////////

/**
 * nearest image : vacuum
 */
inline int math::Boundary_Implementation<math::vacuum>
::nearest_image(Vec const &v1,
		Vec const &v2,
		Vec &nim)const
{
  nim = v1 - v2;
  return 0;
}

/**
 * nearest image : rectangular
 */
inline int math::Boundary_Implementation<math::rectangular>
::nearest_image(Vec const &v1, Vec const &v2,
		Vec &nim)const
{
  // nim = v1 - v2;

  int na = 0, nb = 0, nc = 0;
  
  nim[0] = v1[0] - v2[0];
  while (nim[0] > m_half_box[0]) {
    nim[0] -= 2. * m_half_box[0];
    na = -9;
  }
  while (nim[0] < -m_half_box[0]) {
    nim[0] += 2. * m_half_box[0];
    na = 9;
  }
  nim[1] = v1[1] - v2[1];
  while (nim[1] > m_half_box[1]) {
    nim[1] -= 2. * m_half_box[1];
    nb = -3;
  }
  while (nim[1] < -m_half_box[1]) {
    nim[1] += 2. * m_half_box[1];
    nb = 3;
  }
  nim[2] = v1[2] - v2[2];
  while (nim[2] > m_half_box[2]) {
    nim[2] -= 2. * m_half_box[2];
    nc = -1;
  }
  while (nim[2] < -m_half_box[2]) {
    nim[2] += 2. * m_half_box[2];
    nc = 1;
  }
  
  return na + nb + nc;

  /*
  for(int d=0; d<3; ++d){
    nim(d) = v1(d) - v2(d);

    if (fabs(nim(d)) >= m_half_box(d)){
      nim(d) -= abs(m_box(d)) * rint(nim(d)/abs(m_box(d)));

    }
  }
   */
}      

/**
 * nearest image : triclinic
 */
inline int math::Boundary_Implementation<math::triclinic>
::nearest_image(Vec const &v1,
		Vec const &v2,
		Vec &nim)const
{
  // nim has to be out of the loop here!
  nim = v1 - v2;
  //for(int d=0; d<3; ++d){
  while (nim(2) >= 0.5*m_box(2,2)) {
    nim -= m_box(2);
  }
  while (nim(2) <= -0.5*m_box(2,2)) {
    nim += m_box(2);
  }
  
  double nim_y = nim(1) - m_box(2,1)/m_box(2,2)*nim(2);
  while (nim_y >= 0.5*m_box(1,1)) {
    nim -= m_box(1);
    nim_y = nim(1) - m_box(2,1)/m_box(2,2)*nim(2);
  }
  while (nim_y <= -0.5*m_box(1,1)) {
    nim += m_box(1);
    nim_y = nim(1) - m_box(2,1)/m_box(2,2)*nim(2);
  }

  double nim_x = nim(0)- 
                 m_box(1,0)*(nim(1)-m_box(2,1)*nim(2)/m_box(2,2))/m_box(1,1) -
                 m_box(2,0)*nim(2)/m_box(2,2);
  while (nim_x >= 0.5*m_box(0,0)) {
    nim -= m_box(0);
    nim_x = nim(0)- 
            m_box(1,0)*(nim(1)-m_box(2,1)*nim(2)/m_box(2,2))/m_box(1,1) -
            m_box(2,0)*nim(2)/m_box(2,2);
  }
  while (nim_x <= -0.5*m_box(0,0)) {
    nim += m_box(0);
    nim_x = nim(0)- 
            m_box(1,0)*(nim(1)-m_box(2,1)*nim(2)/m_box(2,2))/m_box(1,1) -
            m_box(2,0)*nim(2)/m_box(2,2);
  }
   
  return 0;
   
  /*
  //for(int d=0; d<3; ++d){
    // i think the if statement might be wrong for really 
    // triclinic cases! - > agree
    // now we are in the rotated frame of the box: 
    // a along x, b in the x-y plane, c arbitrary
    // - > triangular matrix - > trivially solvable set of equations
    // only c has components along z direction
    if(nim(2)*nim(2) >= 0.25*m_box(2,2)*m_box(2,2))
      nim -= m_box(2) * rint(nim(2)/fabs(m_box(2,2)));  
    //b is along x and y
   double nim_y=nim(1)-m_box(2,1)/m_box(2,2)*nim(2);
    if (nim_y*nim_y >= 0.25*m_box(1,1)*m_box(1,1))
      nim -= m_box(1) * rint(nim_y/fabs(m_box(1,1)));
    //a is along x
   double nim_x=nim(0)- 
            m_box(1,0)*(nim(1)-m_box(2,1)*nim(2)/m_box(2,2))/m_box(1,1)
            - m_box(2,0)*nim(2)/m_box(2,2);
    if (nim_x*nim_x >= 0.25*m_box(0,0)* m_box(0,0))
      nim -= m_box(0) * rint(nim_x/fabs(m_box(0,0)));
 // }
   return 0;
  */
}

#ifdef __AVX2__
/**
 * nearest image : rectangular
 */
inline void math::Boundary_Implementation<math::rectangular>
::nearest_image(const double *v1, const double *v2, double *nim
              , const unsigned *ids1, const unsigned *ids2
              , int *ret, const unsigned N) const
/** Takes 4*N vectors and calculates their nims
 */
{
  // nim = v1 - v2;
  static const int ns[3] = {9, 3, 1};
  static const unsigned stride[4] = {0, 8*3, 8*6, 8*9};
  __m128i vindex = _mm_loadu_si128((__m128i*)stride);
  for (unsigned i = 0; i < N; ++i) {
    __m256d m256_v1, m256_v2, m256_gt, m256_lt;
    __m256i m256i_n = _mm256_set1_epi64x(0);
    __m256d m256_nims[3];
    for (unsigned j = 0; j < 3; ++j) {
      // initialize box vectors
      __m256d m256_half_box = _mm256_set1_pd(m_half_box[i]);
      __m256d m256_2half_box = m256_half_box * 2;
      // initialize vectors
      m256_v1 = _mm256_i32gather_pd(v1, vindex, 1);
      m256_v2 = _mm256_i32gather_pd(v2, vindex, 1);
      m256_nims[j] = _mm256_sub_pd(m256_v1, m256_v2);

      //greater-than, ordered, signaling
      m256_gt = _mm256_cmp_pd(m256_nims[j], m256_half_box, _CMP_GT_OS);
      //less-than, ordered, signaling
      m256_lt = _mm256_cmp_pd(m256_nims[j], -m256_half_box, _CMP_LT_OS);
      m256i_n = _mm256_add_epi64(
                  m256i_n
                , _mm256_and_si256(
                    _mm256_castpd_si256(m256_gt)
                  , _mm256_set1_epi64x(-ns[j])
                  )
                );
      m256i_n = _mm256_add_epi64(
                  m256i_n
                , _mm256_and_si256(
                    _mm256_castpd_si256(m256_lt)
                  , _mm256_set1_epi64x(ns[j])
                  )
                );
      __m256d m256_mask = _mm256_or_pd(m256_gt, m256_lt);
      while (_mm256_movemask_pd(m256_mask) != 0) {
        __m256d m256_shift_down = _mm256_and_pd(m256_gt, -m256_2half_box);
        __m256d m256_shift_up = _mm256_and_pd(m256_lt, m256_2half_box);
        __m256d m256_shift = _mm256_or_pd(m256_shift_up, m256_shift_down);
        m256_nims[j] = _mm256_add_pd(m256_nims[j], m256_shift);
        m256_gt = _mm256_cmp_pd(m256_nims[j], m256_half_box, _CMP_GT_OS);
        m256_lt = _mm256_cmp_pd(m256_nims[j], -m256_half_box, _CMP_LT_OS);
        m256_mask = _mm256_or_pd(m256_gt, m256_lt);
      }
      m256i_n = _mm256_add_epi64(
                  m256i_n
                , _mm256_and_si256(
                    _mm256_castpd_si256(m256_lt)
                  , _mm256_set1_epi64x(ns[j])
                  )
                );
    }
    for (unsigned j = 0; j < 3; ++j) {
      // TODO: Permute
      // current alignment nim1_x, ... , nim4_x, nim1_y, ... , nim4_y, nim1_z, ... , nim4_z
      // should be changed to nim1_x, nim1_y, nim1_z, ..., nim4_x, nim4_y, nim4_z,
      _mm256_storeu_pd(&nim[4*(3*i+j)], m256_nims[j]);
    }
    _mm256_storeu_si256((__m256i*)ret, m256i_n);
  }

/*
  io::messages.add("AVX2 instructions not supported in this build",
                  "Boundary_Implementation", io::message::critical);
  return;
*/
}

/**
 * nearest image : rectangular
 */
inline void math::Boundary_Implementation<math::rectangular>
::nearest_image(const __m256d m256_v1[3], const __m256d m256_v2[3], __m256d nim[3]
              , __m128i *ret) const
/** Takes 4*N vectors and calculates their nims
 * other possible versions 
 * direct loadu, would need rearrangement of data using _mm256_permute4x64_pd and _mm256_permute2f128_pd
 * or use _mm256_i64gather_pd for striding (cool!)
 *::nearest_image(VArray const &v1, VArray const &v2, VArray &nim)const
 * single vector with 4 others
 *::nearest_image(Vec const &v1, VArray const &v2, VArray &nim)const
 */
{
  // nim = v1 - v2;
  static const int ns[3] = {9, 3, 1};
  __m256d m256_gt, m256_lt;
  __m256i m256i_n = _mm256_set1_epi64x(0);
  __m256d m256_half_box[3];
  __m256d m256_2half_box[3];
  for (unsigned i = 0; i < 3; ++i) {
    m256_half_box[i] = _mm256_set1_pd(m_half_box[i]);
    m256_2half_box[i] = m256_half_box[i] * 2;
  }

  for (unsigned i = 0; i < 3; ++i) {
    nim[i] = _mm256_sub_pd(m256_v1[i], m256_v2[i]);

    //greater-than, ordered, signaling
    m256_gt = _mm256_cmp_pd(nim[i], m256_half_box[i], _CMP_GT_OS);
    //less-than, ordered, signaling
    m256_lt = _mm256_cmp_pd(nim[i], -m256_half_box[i], _CMP_LT_OS);
    m256i_n = _mm256_sub_epi32(
                m256i_n
              , _mm256_and_si256(
                  _mm256_castpd_si256(m256_gt)
                , _mm256_set1_epi32(ns[i])
                )
              );
    m256i_n = _mm256_add_epi32(
                m256i_n
              , _mm256_and_si256(
                  _mm256_castpd_si256(m256_lt)
                , _mm256_set1_epi32(ns[i])
                )
              );
    __m256d m256_mask = _mm256_or_pd(m256_gt, m256_lt);
    while (_mm256_movemask_pd(m256_mask) != 0) {
      const __m256d m256_shift_down = _mm256_and_pd(m256_gt, -m256_2half_box[i]);
      const __m256d m256_shift_up = _mm256_and_pd(m256_lt, m256_2half_box[i]);
      const __m256d m256_shift = _mm256_or_pd(m256_shift_up, m256_shift_down);
      nim[i] = _mm256_add_pd(nim[i], m256_shift);
      m256_gt = _mm256_cmp_pd(nim[i], m256_half_box[i], _CMP_GT_OS);
      m256_lt = _mm256_cmp_pd(nim[i], -m256_half_box[i], _CMP_LT_OS);
      m256_mask = _mm256_or_pd(m256_gt, m256_lt);
    }
  }
  static const __m256i shuffle = {(long long)2<<32|0, (long long)6<<32|4, 0, 0};
  *ret = _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(m256i_n, shuffle));
}
#endif // __AVX2__

        
////////////////////////////////////////////////////////////////////////////////
// grid stuff
////////////////////////////////////////////////////////////////////////////////


/**
 * calculate box components of vector v : vacuum
 * (lattice vector multipliers)
 */
inline void math::Boundary_Implementation<math::vacuum>
::box_components(Vec const &v, Vec & n)const
{
  n = 0;
}

/**
 * calculate box components of vector v : rectangular
 * (lattice vector multipliers)
 */
inline void math::Boundary_Implementation<math::rectangular>
::box_components(Vec const &v, Vec & n)const
{
  for(int d=0; d<3; ++d){
    n(d) = v(d) / m_box(d)(d);
  }
}

/**
 * calculate box components of vector v : triclinic
 * (lattice vector multipliers)
 */
inline void math::Boundary_Implementation<math::triclinic>
::box_components(Vec const &v, Vec & n)const
{
  for(int d=0; d<3; ++d){
    n(d) = -dot(m_cross_K_L_M(d), v);
  }

}

/**
 * recalc shift vectors
 * and also update cell index shifts : rectangular
 */
inline void math::Boundary_Implementation<math::rectangular>
::recalc_shift_vectors(unsigned int num_cells[3])
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].cell[0] = k * num_cells[0];
	m_shift[index].cell[1] = l * num_cells[1];
	m_shift[index].cell[2] = m * num_cells[2];

	m_shift[index].pos = 
	  k * box()(0) +
	  l * box()(1) +
	  m * box()(2);
	
      }
    }
  }  
}

/**
 * recalc shift vectors : rectangular
 */
inline void math::Boundary_Implementation<math::rectangular>
::recalc_shift_vectors()
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){
	m_shift[index].pos(0) = k * box(0,0);
	m_shift[index].pos(1) = l * box(1,1);
	m_shift[index].pos(2) = m * box(2,2);
      }
    }
  }  
}

/**
 * recalc shift vectors
 * and also update cell index shifts : triclinic
 */
inline void math::Boundary_Implementation<math::triclinic>
::recalc_shift_vectors(unsigned int num_cells[3])
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].cell[0] = k * num_cells[0];
	m_shift[index].cell[1] = l * num_cells[1];
	m_shift[index].cell[2] = m * num_cells[2];

	m_shift[index].pos = 
	  k * box()(0) +
	  l * box()(1) +
	  m * box()(2);
	
      }
    }
  }  
}

/**
 * recalc shift vectors : triclinic
 */
inline void math::Boundary_Implementation<math::triclinic>
::recalc_shift_vectors()
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].pos = 
	  k * box()(0) +
	  l * box()(1) +
	  m * box()(2);
	
      }
    }
  }  
}

/**
 * recalc shift vectors: vacuum
 * (everything is 0)
 */
inline void math::Boundary_Implementation<math::vacuum>
::recalc_shift_vectors(unsigned int num_cells[3])
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].cell[0] = 0;
	m_shift[index].cell[1] = 0;
	m_shift[index].cell[2] = 0;

	m_shift[index].pos = 
	  0 * box()(0) +
	  0 * box()(1) +
	  0 * box()(2);
	
      }
    }
  }  
}

/**
 * recalc shift vector : vacuum
 * (everything 0)
 */
inline void math::Boundary_Implementation<math::vacuum>
::recalc_shift_vectors()
{
  int index=0;
  for(int k=-1; k<2; ++k){
    for(int l=-1; l<2; ++l){
      for(int m=-1; m<2; ++m, ++index){

	m_shift[index].pos = 0;
	
      }
    }
  }  
}

