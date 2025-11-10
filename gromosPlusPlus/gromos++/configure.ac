dnl configure.in for GROMOS++
dnl Written by Roland Buergi, modified for gsl by M. Kastenholz
dnl and again by Markus Christen
dnl 12/20/00

AC_PREREQ(2.53)

AC_INIT([gromos++], [1.6.0], [biomos@igc.phys.chem.ethz.ch], [gromos++], [https//www.gromos.net])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_LANG(C++)

dnl include version information
sinclude([VERSION])

GROMOS_VERSION=$GROMOS_MAJOR_VERSION.$GROMOS_MINOR_VERSION.$GROMOS_MICRO_VERSION
AC_DEFINE_UNQUOTED(GROMOS_VERSION, "$GROMOS_VERSION", [Gromos version info])
AC_DEFINE_UNQUOTED(GROMOS_MAJOR_VERSION, $GROMOS_MAJOR_VERSION, [Gromos major version info])
AC_DEFINE_UNQUOTED(GROMOS_MINOR_VERSION, $GROMOS_MINOR_VERSION, [Gromos minor version info])
AC_DEFINE_UNQUOTED(GROMOS_MICRO_VERSION, $GROMOS_MICRO_VERSION, [Gromos micro version info])

CONF_DATE=`date`
AC_DEFINE_UNQUOTED(GROMOS_DATE, "$CONF_DATE", [MD date (configure)])

dnl large file support
LFS=`getconf LFS_CFLAGS`
AC_SUBST(LFS)

echo ""
echo "GROMOSXX analysis package (GROMOS++)"
echo "version ${GROMOS_VERSION}"
echo ""
echo ""

dnl option for production builds
AC_ARG_ENABLE(debug,
  [  --enable-debug          compile with debug options],
  ,enable_debug=no)

AC_ARG_ENABLE(openmp,
  [  --enable-openmp         compile with openMP support],
  , enable_openmp=no)

AC_ARG_ENABLE(profile,
  [  --enable-profile        compile with profiling on],
  , enable_profile=no)

AC_ARG_ENABLE(profileopt,
  [  --enable-profileopt     use profile information],
  , enable_profileopt=no)

AC_ARG_ENABLE(gprof,
  [  --enable-gprof          compile with profiling on],
  , enable_gprof=no)

AC_ARG_ENABLE(shared,
[  --enable-shared          compile with shared libraries],
    , enable_shared=no)


AM_INIT_AUTOMAKE([no-define subdir-objects])
AM_CONFIG_HEADER(config.h)

dnl Initialise libtool
AM_PROG_LIBTOOL

dnl Check if all programmes are there
AC_PROG_CXX
AC_PROG_MAKE_SET

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

dnl set some defaults to CXXFLAGS and LDFLAGS
extra_lib=""

dnl take out -O and -g from CXXFLAGS, CFLAGS, as these are set
dnl automagically

CXXFLAGS=`echo $CXXFLAGS | sed 's/\-O.//g' | sed 's/\-g//g'`
CFLAGS=`echo $CFLAGS | sed 's/\-O.//g' | sed 's/\-g//g'`
LDFLAGS=`echo $LDFLAGS | sed 's/\-O.//g' | sed 's/\-g//g'`

case "$compiler" in
  intel*)             # Intel C++.  tested with version 8.0beta (precompiled headers)
    echo "Configuring for Intel C++"
    CXX_PICKY_FLAGS="-w1"
    dnl consider adding -xP, -xN or -xW for Intel processors
    dnl maybe some inline control: -Qoption,c,-ip_ninl_max_stats=2000
    dnl or maybe -ipo if your lucky...
    dnl CXX_OPTIMIZE_FLAGS="-O3 -ipo -restrict -align -DNDEBUG -vec_report"
    CXX_OPTIMIZE_FLAGS="-O3 -ip -restrict -align -DNDEBUG -vec_report"
    CXX_DEBUG_FLAGS="-g"
    CXXLD="xild"

    dnl OpenMP
    if eval "test x$enable_openmp = xyes"; then
      echo "   OpenMP enabled"
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -openmp -openmp_report -DOMP"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -openmp -openmp_report -DOMP"
    fi

    dnl profiling
    if eval "test x$enable_profile = xyes"; then
      echo "   profiling enabled"
      mkdir -p profile
      PROFDIR=`pwd`
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -prof-genx -prof-dir=${PROFDIR}/profile"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -prof-genx -prof-dir=${PROFDIR}/profile"
    fi

    dnl profiled optimisation
    if eval "test x$enable_profileopt = xyes"; then
      echo "   using profile data for optimisation"
      PROFDIR=`pwd`
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -prof-use -prof-dir=${PROFDIR}/profile"
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
    CXX_OPTIMIZE_FLAGS="-g -O3 -funroll-loops -fno-gcse -DNDEBUG"
    CXX_DEBUG_FLAGS="-g"
    CXXLD="$CXX"

    #      gprof support...
    if eval "test x$enable_gprof = xyes"; then
      echo "   gprof profiling enabled"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -pg"
    fi

    #      profile support...
    if eval "test x$enable_profile = xyes"; then
      echo "   profiling enabled"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -fprofile-generate"
      LDFLAGS="$LDFLAGS -fprofile-generate"
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
    fi

    dnl OpenMP
    if eval "test x$enable_openmp = xyes"; then
      echo "   OpenMP enabled"
      CXX_DEBUG_FLAGS="$CXX_DEBUG_FLAGS -fopenmp -DOMP"
      CXX_OPTIMIZE_FLAGS="$CXX_OPTIMIZE_FLAGS -fopenmp -DOMP"
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
  *)
    CXXLD="$CXX"
    AC_MSG_WARN([Compiler not recognized])
  ;;
esac

dnl subst compiler flags
AC_SUBST(MY_CXXFLAGS)

CXXFLAGS="${CXXFLAGS} -std=c++11 "

dnl Checks for header files.
AC_HEADER_STDC
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

AC_CHECK_FUNCS(hypot)
AC_CHECK_FUNCS(expm1)
AC_CHECK_FUNCS(acosh)
AC_CHECK_FUNCS(asinh)
AC_CHECK_FUNCS(atanh)

AC_CHECK_FUNCS(isinf, ,[
    AC_MSG_CHECKING([for isinf with <math.h>])
    AC_TRY_LINK([#include <math.h>], [float f = 0.0; isinf(f)],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_ISINF)],
        AC_MSG_RESULT(no))])

AC_CHECK_FUNCS(finite, ,[
    AC_MSG_CHECKING([for finite with <math.h>])
    AC_TRY_LINK([#include <math.h>], [float f = 0.0; finite(f)],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_FINITE)
        ac_cv_func_finite=yes],
        AC_MSG_RESULT(no))])

if test "x${ac_cv_func_finite}" != "xyes" ; then
  AC_CHECK_FUNCS(isfinite, ,[
      AC_MSG_CHECKING([for isfinite with <math.h>])
      AC_TRY_LINK([#include <math.h>], [float f = 0.0; isfinite(f)],[
          AC_MSG_RESULT(yes)
          AC_DEFINE(HAVE_ISFINITE)],
          AC_MSG_RESULT(no))])
fi

AC_CHECK_FUNCS(isnan, ,[
    AC_MSG_CHECKING([for isnan with <math.h>])
    AC_TRY_LINK([#include <math.h>], [float f = 0.0; isnan(f)],[
        AC_MSG_RESULT(yes)
        AC_DEFINE(HAVE_ISNAN)],
        AC_MSG_RESULT(no))])

dnl OpenBSD has a broken implementation of log1p.
case "$host" in
    *-*-*openbsd*)
       AC_MSG_RESULT([avoiding OpenBSD system log1p - using gsl version])
       ;;
    *)
       AC_CHECK_FUNCS(log1p)
       ;;
esac

dnl
dnl FIXME: on Solaris we will need a prototype for hypot, since it is
dnl missing from their math.h according to Jerry.
dnl
dnl extern double hypot(double, double);
dnl

AC_PREFIX_DEFAULT(`pwd`)

dnl Checks for libraries.
dnl Replace `main' with a function in -lm:
AC_CHECK_LIB(m, sqrt)
AC_CHECK_HEADER([zlib.h],[],AC_MSG_ERROR([zlib.h header file missing.]))
AC_CHECK_LIB(z, zlibVersion)

AC_CXX_LIB_GSL([required])
AM_PATH_CCP4_CLIPPER([])
AM_PATH_FFTW3([])
AM_PATH_MDPP([])
AM_PATH_GROMACS([])
AM_PATH_VMD([])

dnl check if full STL support
AC_CV_CXX_VERSION_OK
if test "$ac_cv_cxx_version_ok" = no; then
  echo "Your C++ compiler does not support the STL"
  echo "to the right extent. Check your C++ comipler..."
  exit
fi

AC_CHECK_HEADER([readline/readline.h],
    [AC_DEFINE_UNQUOTED([HAVE_READLINE],[],[Have readline header])],
    [AC_MSG_WARN([Cannot find readline header files. Falling back to stdin])])
AC_CHECK_LIB(readline, readline)
AC_CHECK_HEADER([wordexp.h],
    [AC_DEFINE_UNQUOTED([HAVE_WORDEXP_H],[],[Have wordexp header])],
    [AC_MSG_WARN([Cannot find wordexp.])])
AC_CHECK_FUNCS([wordexp])
AC_CHECK_FUNCS([wordfree])

AC_CACHE_SAVE

AC_SUBST(GROMOS_MAJOR_VERSION)
AC_SUBST(GROMOS_MINOR_VERSION)
AC_SUBST(GROMOS_MICRO_VERSION)
AC_SUBST(GROMOS_VERSION)
AC_SUBST(LSF)

AC_CONFIG_FILES([
	src/Makefile
	src/gmath/Makefile
	src/gromos/Makefile
	src/gcore/Makefile
	src/gio/Makefile
	src/bound/Makefile
	src/fit/Makefile
	src/args/Makefile
	src/utils/Makefile
	src/pb/Makefile
	programs/Makefile
	contrib/Makefile
	data/Makefile
	Makefile
])

AC_OUTPUT

echo ""
echo "configure finished successfully"
echo "type make and make install next"
echo ""

