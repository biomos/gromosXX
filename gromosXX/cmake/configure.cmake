# check compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(COMPILER_PGICC 1)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(COMPILER_GCC 1)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(COMPILER_ICC 1)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "SunPro")
    set(COMPILER_SUNCC 1)
endif()

# check if clipper was found
if(CLIPPER_FOUND)
    set(HAVE_CLIPPER 1)
endif()

# set fftw3 prefix
set(FFTW_PREFIX fftw_)

# set version
set(MD_VERSION "\"${MD_MAJOR_VERSION}.${MD_MINOR_VERSION}.${MD_MICRO_VERSION}\"")

# set date
string(TIMESTAMP MD_DATE "\"%Y-%m-%d %H:%M:%S UTC\"" UTC)

# set top directory
set(TOP_SOURCE_DIR \"${PROJECT_SOURCE_DIR}\")

# check for headers
include(CheckIncludeFileCXX)
check_include_file_cxx(dlfcn.h HAVE_DLFCN_H)
check_include_file_cxx(inttypes.h HAVE_INTTYPES_H)
check_include_file_cxx(memory.h HAVE_MEMORY_H)
check_include_file_cxx(stdint.h HAVE_STDINT_H)
check_include_file_cxx(stdlib.h HAVE_STDLIB_H)
check_include_file_cxx(strings.h HAVE_STRINGS_H)
check_include_file_cxx(string.h HAVE_STRING_H)
check_include_file_cxx(sys/stat.h HAVE_SYS_STAT_H)
check_include_file_cxx(sys/types.h HAVE_SYS_TYPES_H)
check_include_file_cxx(unistd.h HAVE_UNISTD_H)

include(CheckFunctionExists)
#set(CMAKE_REQUIRED_INCLUDES math.h)
#set(CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(system HAVE_SYSTEM)
check_function_exists(tmpnam HAVE_TMPNAM)
check_function_exists(symlink HAVE_SYMLINK)
check_function_exists(unlink HAVE_UNLINK)
check_function_exists(chdir HAVE_CHDIR)
check_function_exists(getcwd HAVE_GETCWD)
check_function_exists(remove HAVE_REMOVE)
check_function_exists(mkstemp HAVE_MKSTEMP)
check_function_exists(mkdtemp HAVE_MKDTEMP)
check_function_exists(getenv HAVE_GETENV)

# check if xtb is properly installed
if(XTB)
    # global includes / libraries vs target specific ones
    
    set(CMAKE_REQUIRED_INCLUDES ${XTB_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${XTB_LIBRARIES})

    # check header file
    check_include_file_cxx(xtb.h HAVE_XTB_H)
    if(NOT HAVE_XTB_H)
         message(FATAL_ERROR "xtb.h not found - provided path to xtb invalid")
    endif()

    # check simple function call
    check_function_exists(xtb_getAPIVersion HAVE_XTBFUNCTION)
    if(NOT HAVE_XTBFUNCTION)
        message(FATAL_ERROR "xtb_getAPIVersion not found - provided path to xtb invalid")
    endif()
endif()

#INCLUDE (CheckSymbolExists)
#CHECK_SYMBOL_EXISTS (acosh cmath.h ACOSH_FOUND)

# check for math functions
set(MATH_FUNCTIONS
        expm1
        acosh
        asinh
        atanh
        )
include(CheckCXXSourceCompiles)
foreach(func ${MATH_FUNCTIONS})
    string(TOUPPER ${func} func_upper)
    check_cxx_source_compiles(
            " #include <cmath>
	    int main (int argc, char* argv[]) { ${func} (0.); }
	  "
            HAVE_${func_upper}
    )
endforeach()

check_cxx_source_compiles(
        " #include <cmath>
	int main (int argc, char* argv[]) { hypot(0.,0.);}
	"
        HAVE_HYPOT
)
check_cxx_source_compiles(
        " #include <cmath>
	int main (int argc, char* argv[]) { std::isinf(0.);}
	"
        HAVE_ISINF
)
check_cxx_source_compiles(
        " #include <cmath>
	int main (int argc, char* argv[]) { std::finite(0.);}
	"
        HAVE_FINITE
)
check_cxx_source_compiles(
        " #include <cmath>
	int main (int argc, char* argv[]) { std::isfinite(0.);}
	"
        HAVE_ISFINITE
)
check_cxx_source_compiles(
        " #include <cmath>
	int main (int argc, char* argv[]) { std::isnan(0.);}
	"
        HAVE_ISNAN
)

# check if size_t is defined in sys/types.h
if(HAVE_SYS_TYPES_H)
    check_cxx_source_compiles(
            " #include <sys/types.h>
		int main (int argc, char* argv[]) {size_t x = 0;}
		"
            SIZE_T_DEFINED
    )
    if(NOT SIZE_T_DEFINED)
        set(size_t unsigned int)
    endif()
endif()

# set variables for configuration depending on options
if(FORCEGROUPS)
    set(XXFORCEGROUPS 1)
endif()
if(HEAVISIDE)
    set(XXHEAVISIDE 1)
endif()

# generate config.h
configure_file(${PROJECT_SOURCE_DIR}/cmake/config.h.in ${PROJECT_BINARY_DIR}/config.h)
