dnl configure.in for GromosXX

AC_PREREQ(2.58)

AC_INIT([md++], [1.5.0], [biomos@igc.phys.chem.ethz.ch], [md++], [https//www.gromos.net])
dnl this is a safety check to see wheter --srcdir exists
AC_CONFIG_SRCDIR(src/math/math.cc)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
dnl get canonical target and host information
AC_CANONICAL_TARGET
AC_CANONICAL_HOST
AC_MSG_RESULT([Configure MD++ for $target])

dnl set C++ as main language
AC_LANG(C++)

dnl include version information
sinclude([VERSION])

MD_VERSION=$MD_MAJOR_VERSION.$MD_MINOR_VERSION.$MD_MICRO_VERSION
AC_DEFINE_UNQUOTED(MD_VERSION, "$MD_VERSION", [MD version info])

CONF_DATE=`date`
AC_DEFINE_UNQUOTED(MD_DATE, "$CONF_DATE", [MD date (configure)])

echo ""
echo "GROMOSXX molecular dynamics program (MD++)"
echo "version ${MD_MAJOR_VERSION}.${MD_MINOR_VERSION}.${MD_MICRO_VERSION}"
echo "${CONF_DATE}, ETH Zuerich"
echo ""
echo ""

AM_INIT_AUTOMAKE([no-define subdir-objects])
AM_CONFIG_HEADER(config.h)

dnl option for production builds
AC_ARG_ENABLE(debug,
  [  --enable-debug          compile with debug options],
  ,enable_debug=no)    

AC_ARG_ENABLE(openmp,
  [  --enable-openmp         compile with openMP support],
  , enable_openmp=no)

AC_ARG_ENABLE(mpi,
  [  --enable-mpi            compile with MPI support],
  , enable_mpi=no)

AC_ARG_ENABLE(repex,
  [  --disable-repex         compile without replica exchange support],
  , enable_repex=yes)

AC_ARG_ENABLE(profile,
  [  --enable-profile        compile with profiling on],
  , enable_profile=no)

AC_ARG_ENABLE(profileopt,
  [  --enable-profileopt     use profile information],
  , enable_profileopt=no)

AC_ARG_ENABLE(gprof,
  [  --enable-gprof          compile with profiling on],
  , enable_gprof=no)

AC_ARG_ENABLE(allstatic, 
  [ --enable-allstatic       compile static-only],
  , enable_allstatic=no)

AC_ARG_ENABLE(forcegroups, 
  [ --enable-forcegroups    compile with force group support],
  ,enable_forcegroups=no)

AC_ARG_ENABLE(heaviside, 
  [ --enable-heaviside    compile with Heaviside truncation],
  ,enable_heaviside=no)

AC_ARG_ENABLE(shared,
[  --enable-shared         compile with shared libraries],
    , enable_shared=no)

dnl use empty default flags but allow user to override
: ${CFLAGS=""}
: ${CXXFLAGS=""}
: ${LDFLAGS=""}

if eval "test x$enable_allstatic = xyes"; then
  LDFLAGS="$LDFLAGS -all-static"
  echo "configuring for static build"
fi

if eval "test x$enable_forcegroups = xyes"; then
  AC_DEFINE_UNQUOTED(XXFORCEGROUPS, 1, [enable force groups])
fi

if eval "test x$enable_heaviside = xyes"; then
  AC_DEFINE_UNQUOTED(XXHEAVISIDE, 1, [enable Heaviside truncation])
fi

dnl check for required programs
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET
AC_PROG_INSTALL

if eval "test x$enable_openmp = xyes"; then
  if eval "test x$enable_mpi = xyes"; then
    AC_MSG_ERROR([Either enable OpenMP or MPI, but not both at the same time.])
  fi
fi


dnl detect the compiler
compiler="none"
dnl need to test for intel first, as it looks nowadays
dnl suspiciously like GNU g++
if eval "test x$compiler = xnone"; then AC_PROG_CXX_INTELCC ; fi
if eval "test x$compiler = xnone"; then
  if eval "test x$GXX = xyes"; then
    compiler=gnucc
    echo "using GNU c++ compiler"
  fi
fi
if eval "test x$compiler = xnone"; then AC_PROG_CXX_SUNCC   ; fi
if eval "test x$compiler = xnone"; then AC_PROG_CXX_PGI     ; fi

case "$compiler" in
  intel*)             # Intel C++.  tested with version 8.0beta (precompiled headers)
    echo "Configuring for Intel C++"
    CXX_PICKY_FLAGS="-w1"
    dnl consider adding -xP, -xN or -xW for Intel processors
    dnl maybe some inline control: -Qoption,c,-ip_ninl_max_stats=2000
    dnl or maybe -ipo if your lucky...
    dnl CXX_OPTIMIZE_FLAGS="-O3 -ipo -restrict -align -DNDEBUG -vec_report"
    CXX_OPTIMIZE_FLAGS="-O3 -ip -restrict -align -DNDEBUG -vec_report"
    CXX_NDEBUG_FLAGS="-DNDEBUG -g"
    CXX_DEBUG_FLAGS="-g"
    CXXLD="xild"

    dnl OpenMP
    if eval "test x$enable_openmp = xyes"; then
      echo "   OpenMP enabled"
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -openmp -openmp_report -DOMP"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -openmp -openmp_report -DOMP"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -openmp -openmp_report -DOMP"
    fi

    dnl MPI
    if eval "test x$enable_mpi = xyes"; then
      echo "   MPI enabled"
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
    fi
       
    dnl profiling
    if eval "test x$enable_profile = xyes"; then
      echo "   profiling enabled"
      mkdir -p profile
      PROFDIR=`pwd`
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -prof-genx -prof-dir=${PROFDIR}/profile"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -prof-genx -prof-dir=${PROFDIR}/profile"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -prof-genx -prof-dir=${PROFDIR}/profile"
    fi

    dnl profiled optimisation
    if eval "test x$enable_profileopt = xyes"; then
      echo "   using profile data for optimisation"
      PROFDIR=`pwd`
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -prof-use -prof-dir=${PROFDIR}/profile"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -prof-use -prof-dir=${PROFDIR}/profile"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -prof-use -prof-dir=${PROFDIR}/profile"
    fi

    if eval "test x$enable_debug = xyes"; then
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS $CXX_PICKY_FLAGS"
      MY_CXXFLAGS="$CXX_DEBUG_FLAGS"
      echo "   configuring for debug build"
    else
      MY_CXXFLAGS="$CXX_OPTIMIZE_FLAGS"
      echo "   configuring for release build"
    fi

    AC_DEFINE_UNQUOTED(COMPILER_ICC, 1, [Intel C++ compiler])
  ;;
  gnu*)               # EGCS  http://egcs.cygnus.com/
    echo "Configuring for GNU C++"
    CXX_PICKY_FLAGS="-Wall"
    CXX_OPTIMIZE_FLAGS="-O3 -DNDEBUG"
    CXX_NDEBUG_FLAGS="-DNDEBUG -g"
    CXX_DEBUG_FLAGS="-g"
    CXXLD="$CXX"

    #      gprof support...
    if eval "test x$enable_gprof = xyes"; then
      echo "   gprof profiling enabled"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -pg"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -pg"
    fi

    #      profile support...       
    if eval "test x$enable_profile = xyes"; then
      echo "   profiling enabled"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -fprofile-generate"
      LDFLAGS="$LDFLAGS -fprofile-generate"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -fprofile-generate"
    fi

    if eval "test x$enable_gprof = xno -a x$enable-profile = xno"; then
      # this opt flag is incompatible with profiling
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -fomit-frame-pointer"
    fi

    #      profile optimization
    if eval "test x$enable_profileopt = xyes"; then
      echo "   profiling based optimizations enabled"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -fprofile-use"
      LDFLAGS="$LDFLAGS -fprofile-use"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -fprofile-use"
    fi

    dnl OpenMP
    if eval "test x$enable_openmp = xyes"; then
      echo "   OpenMP enabled"
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -fopenmp -DOMP"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -fopenmp -DOMP"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -fopenmp -DOMP"
    fi

    if eval "test x$enable_mpi = xyes"; then
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
    fi

    #      debug
    if eval "test x$enable_debug = xyes"; then
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS $CXX_PICKY_FLAGS"
      MY_CXXFLAGS="$CXX_DEBUG_FLAGS"
       echo "configuring for debug build"
    else
      MY_CXXFLAGS="$CXX_OPTIMIZE_FLAGS"
      echo "configuring for release build"
    fi

    dnl Why do we need this?
    F77=g77
    F77_OPTIMIZE_FLAGS="-O2 -funroll-loops"
    FORTLIBS=-lg2c
    AC_DEFINE_UNQUOTED(COMPILER_GCC, 1, [GNU C++ compiler])
  ;;
  suncc*)             # SUN studio (tested with SUN studio 10)
    echo "Configuring for SUN Studio"
    dnl CXX_PICKY_FLAGS="+w2 -xcheck -xsb"
    CXX_PICKY_FLAGS="+w2"
    CXX_OPTIMIZE_FLAGS="-fast -xdepend -xvector -DNDEBUG"
    CXX_NDEBUG_FLAGS="-DNDEBUG -g"
    CXX_DEBUG_FLAGS="-g"
    CXXLD="$CXX"

    if eval "test x$enable_openmp = xyes"; then
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -xopenmp -DOMP"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -xopenmp -DOMP"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -xopenmp -DOMP"
    fi

    if eval "test x$enable_mpi = xyes"; then
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
    fi

    if eval "test x$enable_debug = xyes"; then
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS $CXX_PICKY_FLAGS"
      MY_CXXFLAGS="$CXX_DEBUG_FLAGS"
      echo "configuring for debug build"
    else
      MY_CXXFLAGS="$CXX_OPTIMIZE_FLAGS"
      echo "configuring for release build"
    fi

    AC_DEFINE_UNQUOTED(COMPILER_SUNCC, 1, [SUN C++ compiler])
  ;;
  pgicc*)             # PGI compilers
    echo "Configuring for PGI"
    CXX_PICKY_FLAGS="-Minform=warn"
    CXX_OPTIMIZE_FLAGS="-fast -fastsse -DNDEBUG"
    CXX_NDEBUG_FLAGS="-DNDEBUG -g"
    CXX_DEBUG_FLAGS="-g"
    CXXLD="$CXX"

    if eval "test x$enable_openmp = xyes"; then
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -DOMP"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -DOMP"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -DOMP"
    else
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -Mnoopenmp"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -Mnoopenmp"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -Mnoopenmp"
    fi

    if eval "test x$enable_mpi = xyes"; then
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
      CXX_NDEBUG_FLAGS="$CXX_NDEBUG_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -DXXMPI -DMPICH_IGNORE_CXX_SEEK"
    fi

    if eval "test x$enable_debug = xyes"; then
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS $CXX_PICKY_FLAGS"
      MY_CXXFLAGS="$CXX_DEBUG_FLAGS"
      echo "configuring for debug build"
    else
      MY_CXXFLAGS="$CXX_OPTIMIZE_FLAGS"
      echo "configuring for release build"
    fi

    AC_DEFINE_UNQUOTED(COMPILER_PGICC, 1, [PGI compiler])
  ;;
  *)
    CXXLD="$CXX"
    AC_MSG_WARN([Compiler not recognized])
  ;;
esac

dnl set c++11 standard
MY_CXXFLAGS="-std=c++11 ${MY_CXXFLAGS}"

dnl subst compiler flags
AC_SUBST(MY_CXXFLAGS)
AC_SUBST(CXX_NDEBUG_FLAGS)

AC_PREFIX_DEFAULT(`pwd`)

AM_PATH_GSL
AM_PATH_CCP4_CLIPPER
AM_PATH_FFTW3
AM_PATH_HOOMD
AM_PATH_CUDA
AM_PATH_XTB
AM_WITH_SCHNETPACK

dnl check header files
AC_HEADER_STDC
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

AC_CHECK_FUNCS(hypot)
AC_CHECK_FUNCS(expm1)
AC_CHECK_FUNCS(acosh)
AC_CHECK_FUNCS(asinh)
AC_CHECK_FUNCS(atanh)
AC_CHECK_FUNCS(system)
AC_CHECK_FUNCS(mkstemp)
AC_CHECK_FUNCS(mkdtemp)
AC_CHECK_FUNCS(getenv)
AC_CHECK_FUNCS(symlink)
AC_CHECK_FUNCS(unlink)
AC_CHECK_FUNCS(chdir)
AC_CHECK_FUNCS(getcwd)

AC_MSG_CHECKING([for isinf with <cmath>])
AC_TRY_LINK([#include <cmath>], [float f = 0.0; std::isinf(f)],[
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_ISINF],[1],[has isinf])],
  AC_MSG_RESULT(no))

AC_MSG_CHECKING([for finite with <cmath>])
AC_TRY_LINK([#include <cmath>], [float f = 0.0; std::finite(f)],[
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_FINITE],[1],[have finite])
  ac_cv_func_finite=yes],
  AC_MSG_RESULT(no))

if test "x${ac_cv_func_finite}" != "xyes" ; then
  AC_MSG_CHECKING([for isfinite with <cmath>])
  AC_TRY_LINK([#include <cmath>], [float f = 0.0; std::isfinite(f)],[
    AC_MSG_RESULT(yes)
    AC_DEFINE([HAVE_ISFINITE],[1],[have isfinite])],
    AC_MSG_RESULT(no))
fi

AC_MSG_CHECKING([for isnan with <cmath>])
AC_TRY_LINK([#include <cmath>], [float f = 0.0; std::isnan(f)],[
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_ISNAN],[1],[have isnan])],
  AC_MSG_RESULT(no))

dnl OpenBSD has a broken implementation of log1p.
case "$host" in
  *-*-*openbsd*)
    AC_MSG_RESULT([avoiding OpenBSD system log1p - using gsl version])
  ;;
  *)
    AC_CHECK_FUNCS(log1p)
  ;;
esac

dnl Checks for libraries.
AC_CHECK_LIB(m, sqrt)
AC_CHECK_LIB(socket, listen)
AC_CHECK_LIB(nsl, uname)
AC_CHECK_HEADER([zlib.h],[],AC_MSG_ERROR([zlib.h header file missing.]))
AC_CHECK_LIB(z, zlibVersion)

if eval "test x$enable_openmp = xyes"; then
  AC_CHECK_LIB([fftw3_omp], [${PREFIX}init_threads],, AC_MSG_ERROR([FFTW3 threads library missing.]), [-fopenmp -lfftw3 -lm])
fi

dnl check if full STL support
AC_CV_CXX_VERSION_OK
if test "$ac_cv_cxx_version_ok" = no; then
  echo "Your C++ compiler does not support the STL"
  echo "to the right extent. Check your C++ compiler..."
  exit
fi

if eval "test x$enable_repex = xno"; then
  echo "   replica-exchange disabled"
  REPEX=no
else
  REPEX=yes
  AC_SUBST(REPEX)
  AC_DEFINE_UNQUOTED(REPEX, "$REPEX", [replica-exchange simulation])
fi


AC_CACHE_SAVE

AC_SUBST(MD_MAJOR_VERSION)
AC_SUBST(MD_MINOR_VERSION)
AC_SUBST(MD_MICRO_VERSION)
AC_SUBST(MD_VERSION)

AC_SUBST(compiler)

AC_DEFINE_UNQUOTED(TOP_SOURCE_DIR, "$srcdir", [top src directory])

AC_CONFIG_FILES([
	src/Makefile
	src/cukernel/Makefile
	src/math/Makefile
	src/topology/Makefile
	src/configuration/Makefile
	src/algorithm/Makefile
	src/simulation/Makefile
	src/interaction/Makefile
	src/io/Makefile
	src/replicaExchange/Makefile
	src/util/Makefile
	src/check/Makefile
	program/Makefile
	contrib/Makefile
	Makefile
])
	
AC_OUTPUT
