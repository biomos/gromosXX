/**
 * @file stdheader.h
 * file for later use with precompiled headers.
 */

#ifndef INCLUDED_STDHEADER_H
#define INCLUDED_STDHEADER_H

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <stdexcept>
#include <cassert>

#include <algorithm>
#include <typeinfo>
#include <math/gmath.h>
#include <util/debug.h>

#include <random/normal.h>

#include <io/message.h>

#ifdef COMPILER_GCC
#include <cxxabi.h>
using namespace ranlib;
#endif

#endif

