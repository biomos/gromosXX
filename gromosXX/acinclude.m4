dnl Local function, as GROMOS++ depends
dnl on the STL. -> Does the C++ compiler
dnl support the STL to the degree necessary?
dnl 
AC_DEFUN(AC_CV_CXX_VERSION_OK,
  [AC_CACHE_CHECK(whether the compiler supports the STL,
   ac_cv_cxx_version_ok,
     [AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([#include <vector>],[
	std::vector<int> v; v.push_back(0);return 0;],
         ac_cv_cxx_version_ok=yes, ac_cv_cxx_version_ok=no)
      AC_LANG_RESTORE
  ])
])

AC_DEFUN(AC_PROG_CXX_MWERKS,
[AC_CACHE_CHECK(whether we are using Metrowerks Codewarrior C++, MWERKS_CXX,
[cat > conftest.c <<EOF
#ifdef __MWERKS__
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  MWERKS_CXX=yes
  compiler=mwerks
else
  MWERKS_CXX=no
fi])])

AC_DEFUN(AC_PROG_CXX_MSVCPP,
[AC_CACHE_CHECK(whether we are using Visual C++, MSVCPP_CXX,
[cat > conftest.c <<EOF
#if defined(_MSC_VER) && !defined(__MWERKS__) && !defined(__ICL)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  MSVCPP_CXX=yes
  compiler=msvcpp
else
  MSVCPP_CXX=no
fi])])

AC_DEFUN(AC_PROG_CXX_SGICC,
[AC_CACHE_CHECK(whether we are using SGI MIPSpro C++, SGI_CXX,
[cat > conftest.c <<EOF
# if defined(__sgi) && !defined(__GNUC__)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  SGI_CXX=yes
  compiler=sgicc
else
  SGI_CXX=no
fi])])

AC_DEFUN(AC_PROG_CXX_SUNCC,
[AC_CACHE_CHECK(whether we are using Sun C++, SUN_CXX,
[cat > conftest.c <<EOF
# if defined(__SUNPRO_CC) 
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  SUN_CXX=yes
  compiler=suncc
else
  SUN_CXX=no
fi])])

AC_DEFUN(AC_PROG_CXX_INTELCC,
[AC_CACHE_CHECK(whether we are using Intel C++, INTEL_CXX,
[cat > conftest.c <<EOF
# if defined(__ICC)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  INTEL_CXX=yes
  compiler=intelcc
else
  INTEL_CXX=no
fi])])

AC_DEFUN(AC_PROG_CXX_KAICC,
[AC_CACHE_CHECK(whether we are using KAI C++, KAI_CXX,
[cat > conftest.c <<EOF
# if defined(__KCC)
  yes;
#endif
EOF
if AC_TRY_COMMAND(${CXX} -E conftest.c) | egrep yes >/dev/null 2>&1; then
  KAI_CXX=yes
  compiler=kaicc
else
  KAI_CXX=no
fi])])

AC_DEFUN(AC_BZ_SET_COMPILER,
  [cxxwith=`echo $1 | sed -e 's/ /@/'`
   case "$cxxwith" in
     *:*@*)                 # Full initialization syntax
       CXX=`echo "$cxxwith" | sed  -n -e 's/.*:\(.*\)@.*/\1/p'`
       CXXFLAGS=`echo "$cxxwith" | sed  -n -e 's/.*:.*@\(.*\)/\1/p'`
     ;;
     *:*)                   # Simple initialization syntax
       CXX=`echo "$cxxwith" | sed  -n -e 's/.*:\(.*\)/\1/p'`
       CXXFLAGS=$3
     ;;
     *)                     # Default values
       CXX=$2
       CXXFLAGS=$3
     ;;
   esac])

dnl Determine a Fortran 77 compiler to use.  If `F77' is not already set
dnl in the environment, check for `g77', `f77' and `f2c', in that order.
dnl Set the output variable `F77' to the name of the compiler found.
dnl 
dnl If using `g77' (the GNU Fortran 77 compiler), then `AC_PROG_F77'
dnl will set the shell variable `G77' to `yes', and empty otherwise.  If
dnl the output variable `FFLAGS' was not already set in the environment,
dnl then set it to `-g -02' for `g77' (or `-O2' where `g77' does not
dnl accept `-g').  Otherwise, set `FFLAGS' to `-g' for all other Fortran
dnl 77 compilers.
dnl 
dnl AC_PROG_F77()
AC_DEFUN(AC_MTL_PROG_F77,
[AC_BEFORE([$0], [AC_PROG_CPP])dnl
if test -z "$F77"; then
  AC_CHECK_PROGS(F77, g77 f77 f2c)
    test -z "$F77" && AC_MSG_WARN([no acceptable Fortran 77 compiler found in \$PATH])
fi

AC_PROG_F77_WORKS
AC_PROG_F77_GNU

if test $ac_cv_prog_g77 = yes; then
  G77=yes
dnl Check whether -g works, even if FFLAGS is set, in case the package
dnl plays around with FFLAGS (such as to build both debugging and
dnl normal versions of a library), tasteless as that idea is.
  ac_test_FFLAGS="${FFLAGS+set}"
  ac_save_FFLAGS="$FFLAGS"
  FFLAGS=
  AC_PROG_F77_G
  if test "$ac_test_FFLAGS" = set; then
    FFLAGS="$ac_save_FFLAGS"
  elif test $ac_cv_prog_f77_g = yes; then
    FFLAGS="-g -O2"
  else
    FFLAGS="-O2"
  fi
else
  G77=
  test "${FFLAGS+set}" = set || FFLAGS="-g"
fi
])

dnl @synopsis AC_CXX_LIB_BLITZ([optional-string "required"])
dnl
dnl Check whether Blitz++ is installed.
dnl Blitz++ is available at http://oonumerics.org/blitz.
dnl
dnl   Set the path for Blitz++  with the option
dnl      --with-blitz[=DIR]
dnl   Blitz headers should be under DIR/includes
dnl   Blitz library should be under DIR/lib
dnl   Then try to compile and run a simple program with a Blitz Array
dnl   Optional argument `required' triggers an error if Blitz++ not installed
dnl 
dnl @version $Id$
dnl @author Patrick Guio <patrick.guio@matnat.uio.no>
dnl
AC_DEFUN([AC_MSG_ERROR_BLITZ],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the Blitz++ template library
available at http://oonumerics.org/blitz
When installed give the directory of installation with the option
  --with-blitz@<:@=DIR@:>@
])])


AC_DEFUN([AC_CXX_LIB_BLITZ],[

AC_ARG_WITH(blitz,
AS_HELP_STRING([--with-blitz@<:@=DIR@:>@],[Set the path for Blitz++]),
[],[withval='yes'])

if test "$1" = required -a "$withval" = no ; then
	AC_MSG_ERROR_BLITZ
fi

if test "$withval" != no ; then

	saveCPPFLAGS=$CPPFLAGS
	saveLDFLAGS=$LDFLAGS
	saveLIBS=$LIBS

	if test "$withval" != 'yes'; then
		CPPFLAGS="-I$withval/include"
		LDFLAGS="-L$withval/lib"
	fi
	LIBS="-lblitz"

	AC_CACHE_CHECK([whether Blitz++ is installed],ac_cxx_lib_blitz,
	[AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[
#include <blitz/array.h>
]],[[
blitz::Array<int,1> x(10);
x = blitz::tensor::i;
	]])],[ac_cxx_lib_blitz=yes],[ac_cxx_lib_blitz=no])
	AC_LANG_RESTORE
	])

	CPPFLAGS=$saveCPPFLAGS
	LDFLAGS=$saveLDFLAGS
	LIBS=$saveLIBS

	if test "$ac_cxx_lib_blitz" = yes ; then
		if test "$withval" != yes ; then
			CPPFLAGS="-I$withval/include"
			BLITZ_LDFLAGS="-L$withval/lib"
dnl			AC_SUBST(BLITZ_CPPFLAGS)
			AC_SUBST(BLITZ_LDFLAGS)
		fi
		BLITZ_LIB="-lblitz"
		AC_SUBST(BLITZ_LIB)
	else
		if test "$1" = required ; then
			AC_MSG_ERROR_BLITZ
		fi
	fi

fi

])

