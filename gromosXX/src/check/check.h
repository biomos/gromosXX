/**
 * @file check.h
 * functions and macros to make testing easier.
 */

#define to_str(name) # name
#define STR2(name) to_str(name)

#define CHECK_MODULE(name, result) \
{ \
  std::cout << "\n\033[22;34mentering " << STR2(name) \
            << "\033[0m" << std::endl; \
  result += name (); \
}

#define CHECK_RESULT(res, total) \
if (res){ \
  std::cout << std::setw(10) << std::right \
  << "\033[1;31mfailed\033[22;0m" << std::endl; \
  total += res; \
} else { \
  std::cout << std::setw(10) << std::right \
  << "\033[1;32mok\033[22;0m" << std::endl; \
}

#define CHECK_EQUAL(x1, x2, res) \
if (x1 != x2){ \
  std::cout << "\n\t" << STR(x1) << " =\t" << x1 \
            << "\n\t" << STR(x2) << " =\t" << x2 \
            << "\n\t" << std::setw(60) << " "; \
  ++res; \
}
 /*   std::cout.precision(6); \
    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield); \
    std::cout << "\n\t" << std::setw(15) << STR(x1) << " =\t" << x1 \
              << "\n\t" << std::setw(15) << STR(x2) << " =\t" << x2 \
              << "\n\t" << std::setw(15) << "rel error" << " =\t" << diff \
              << "\n\t" << std::setw(60) << " "; \
   */
#define CHECK_APPROX_EQUAL(x1, x2, eps, res) \
{ \
  double diff; \
  if (x1 != 0) diff = fabs((x1 - x2) / x1); \
  else diff = x2; \
  if (diff > eps){ \
    ++res; \
  } \
}

#define CHECK_APPROX_EQUAL_RMSFERR(x1, x2, eps, f, res) \
{ \
  double diff; \
  if (x1 != 0) diff = fabs(max(fabs(x1 - x2) - f, 0.0) / x1); \
  else diff = max(x2 - f, 0.0); \
  if (diff > eps){ \
    ++res; \
  } \
}

#define CHECKING(name, res) \
{ \
std::cout << "\033[34mchecking\033[0m " << std::setw(70) \
          << std::left << name << std::flush; \
res = 0; \
}

#define RESULT(last, total) \
if(last == 0) {\
  std::cout << std::setw(10) << std::right \
  << "\033[1;32mok\033[22;0m" << std::endl; \
} else { \
  std::cout << std::setw(10) << std::right \
  << "\033[1;31mfailed\033[22;0m" << std::endl; \
  total += last; \
}

#define GETFILE(ifs, name) \
{ \
string s = "../../" TOP_SOURCE_DIR "src/check/data/" name; \
ifs.open(s.c_str()); \
}

#define GETFILEPATH(p, name, subdir) \
{ \
  p = "../../" TOP_SOURCE_DIR "/" subdir name; \
}
