dnl Local function, as GROMOS++ depends
dnl on the STL. -> Does the C++ compiler
dnl support the STL to the degree necessary?
dnl 
AC_DEFUN([AC_CV_CXX_VERSION_OK],
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

dnl @synopsis AC_CXX_LIB_GSL([optional-string "required"])
dnl
dnl Check whether Gnu Scientific Library (GSL) is installed.
dnl GSL is available from
dnl www.gnu.org
dnl
dnl   Set the path for GSL  with the option
dnl      --with-gsl[=DIR]
dnl   GSL headers should be under DIR/include
dnl   GSL library should be under DIR/lib
dnl   Then try to compile and run a simple program with a gsl random number
dnl   Optional argument `required' triggers an error if GSL not installed
dnl 
dnl @version $Id$
dnl @author Patrick Guio <patrick.guio@matnat.uio.no>
dnl
AC_DEFUN([AC_MSG_ERROR_GSL],[
AC_MSG_ERROR([
$PACKAGE_STRING requires the Gnu Scientific Library (GSL)
When installed give the directory of installation with the option
  --with-gsl@<:@=DIR@:>@
])])


AC_DEFUN([AC_CXX_LIB_GSL],[

AC_ARG_WITH(gsl,
AS_HELP_STRING([--with-gsl@<:@=DIR@:>@],[Set the path for GSL]),
[],[withval='yes'])

if test "$1" = required -a "$withval" = no ; then
	AC_MSG_ERROR_GSL
fi

if test "$withval" != no ; then

	saveCPPFLAGS=$CPPFLAGS
	saveLDFLAGS=$LDFLAGS
	saveLIBS=$LIBS

	if test "$withval" != 'yes'; then
		CPPFLAGS="-I$withval/include"
		LDFLAGS="-L$withval/lib -Wl,-rpath,$withval/lib"
	fi
	LIBS="-lgsl -lgslcblas"

	AC_CACHE_CHECK([whether Gnu Scientific Library is installed],ac_cv_cxx_lib_gsl,
	[AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[
#include <gsl/gsl_matrix.h>
]],[[
gsl_matrix * mat = gsl_matrix_alloc(3,3);
gsl_matrix_set_zero(mat);
gsl_matrix_free(mat);
	]])],[ac_cxx_lib_gsl=yes],[ac_cxx_lib_gsl=no])
	AC_LANG_RESTORE
	])

	CPPFLAGS=$saveCPPFLAGS
	LDFLAGS=$saveLDFLAGS
	LIBS=$saveLIBS

	if test "$ac_cxx_lib_gsl" = yes ; then
		if test "$withval" != yes ; then
			CPPFLAGS="$CPPFLAGS -I$withval/include"
			GSL_LDFLAGS="-L$withval/lib -Wl,-rpath,$withval/lib"
			AC_SUBST(GSL_LDFLAGS)
		fi
		GSL_LIB="-lgsl -lgslcblas"
		AC_SUBST(GSL_LIB)

   		AC_DEFINE_UNQUOTED([HAVE_GSL],[],[Gnu Scientific Library])

	else
		if test "$1" = required ; then
			AC_MSG_ERROR_GSL
		fi
	fi

fi

])
AC_DEFUN([AC_PROG_CXX_SUNCC],
[AC_CACHE_CHECK(whether we are using Sun C++, SUN_cv_CXX,
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

AC_DEFUN([AC_PROG_CXX_INTELCC],
[AC_CACHE_CHECK(whether we are using Intel C++, INTEL_cv_CXX,
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

dnl check for lib CCP4/Clipper
AC_DEFUN([AM_PATH_CCP4_CLIPPER],[
  dnl allow for ccp4 lib directory specification
  AC_ARG_WITH(ccp4,
    [  --with-ccp4=DIR         CCP4 library directory to use],
    [
      [CXXFLAGS="$CXXFLAGS -I${withval}/include -L${withval}/lib"]
      [LDFLAGS="$LDFLAGS -L${withval}/lib"]
    ],
    [
      AC_MSG_WARN([Assuming default paths for CCP4])
    ])
  AC_ARG_WITH(clipper,
    [  --with-clipper=DIR      clipper library directory to use],
    [
      [CXXFLAGS="$CXXFLAGS -I${withval}/include -L${withval}/lib"]
      [LDFLAGS="$LDFLAGS -L${withval}/lib"]
      [CLIPPER_LIB="-lclipper-ccp4 -lccp4c -lclipper-contrib -lclipper-core -lrfftw -lfftw -lpthread"]
      AC_DEFINE_UNQUOTED([HAVE_CLIPPER],[],[Have clipper x-ray library])
      [have_clipper=yes]
    ],
    [
      AC_MSG_WARN([clipper path was not specified. Disabling clipper support])
      [CLIPPER_LIB=""]
      [have_clipper=no]
    ]
  )
  dnl check for lib with these settings. To be implemented
  AC_SUBST(CLIPPER_LIB)
])

dnl check for lib FFTW3
AC_DEFUN([AM_PATH_FFTW3],[
  PREFIX='fftw_'
  if eval "test x$have_clipper = xyes"; then
    PREFIX='fftw3_'
  fi
  AC_DEFINE_UNQUOTED([FFTW_PREFIX], [${PREFIX}], [prefix for fftw3 library])
  AC_ARG_WITH(fftw,
    [  --with-fftw=DIR         fftw library directory to use],
    [
      [CXXFLAGS="$CXXFLAGS -I${withval}/include -L${withval}/lib"]
      [LDFLAGS="$LDFLAGS -L${withval}/lib"]
    ],
    [
      AC_MSG_WARN([fftw path was not specified. Trying default paths...])
    ]
  )
  dnl check for lib with these settings and add flags automatically
  AC_CHECK_LIB([fftw3], [${PREFIX}version],,, [-lm])
])


dnl check for md++
AC_DEFUN([AM_PATH_MDPP],[
  AC_ARG_WITH(mdpp,
    [  --with-mdpp=DIR         GROMOS MD++ directory],
    [
      [CXXFLAGS="$CXXFLAGS -I${withval}/include -L${withval}/lib"]
      [LDFLAGS="$LDFLAGS -L${withval}/lib"]
      [MDPP_LIB="-lmdpp"]
      AC_DEFINE_UNQUOTED([HAVE_MDPP],[],[Have md++ lib])
      AC_DEFINE_UNQUOTED([GROMOSXX],[],[MD++])
    ],
    [
      AC_MSG_WARN([GROMOS MD++ path was not specified. Disabling MD++ support.])
      [MDPP_LIB=""]
    ]
  )
  AC_SUBST(MDPP_LIB)
])

dnl check for gromacs
AC_DEFUN([AM_PATH_GROMACS],[
  AC_ARG_WITH(gromacs,
    [  --with-gromacs=DIR      Gromacs directory],
    [
      [CXXFLAGS="$CXXFLAGS -I${withval}/include -I${withval}/include/gromacs -L${withval}/lib"]
      [LDFLAGS="$LDFLAGS -L${withval}/lib"]
      [GMX_LIB="-lgmx"]
      AC_DEFINE_UNQUOTED([HAVE_GMX],[],[Have Gromacs lib])
    ],
    [
      AC_MSG_WARN([Gromacs path was not specified. Disabling Gromacs support.])
      [GMX_LIB=""]
    ]
  )
  AC_SUBST(GMX_LIB)
])

dnl check for VMD
AC_DEFUN([AM_PATH_VMD],[
  AC_ARG_WITH(vmd,
    [  --with-vmd=DIR      VMD plugin sources directory],
    [
      [CXXFLAGS="$CXXFLAGS -I${withval}/include "]
      [CPPFLAGS="$CPPFLAGS -I${withval}/include "]
    ],[
    AC_MSG_WARN([Searching default include paths for VMD plugin headers...])
    ]
  )
  AC_CHECK_HEADER([molfile_plugin.h],
    [AC_DEFINE_UNQUOTED([HAVE_VMD],[],[Have VMD plugin header])],
    [AC_MSG_WARN([Cannont find VMD header files. Disabling VMD plugin])])
])
