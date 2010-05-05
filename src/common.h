// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __HERMES2D_COMMON_H
#define __HERMES2D_COMMON_H

#include "config.h"

// common headers
#include <stdexcept>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> //allows to use offsetof
#include <string.h>
#include <cstdarg>
#include <assert.h>
#include <pthread.h>
#include <math.h>

#include <float.h>

// STL stuff
#include <vector>
#include <string>
#include <map>
#include <set>
#include <queue>
#include <sstream>
#include <fstream>

// platform compatibility stuff
#include "compat.h"

// others
#include <Judy.h>
#include "auto_local_array.h"
#include "common_time_period.h"

// Enabling second derivatives in weak forms. Turned off by default. Second
// derivatives are employed, among others, by stabilization methods for
// transport equations. For usage see the example linear-convection-diffusion.
#define H2D_SECOND_DERIVATIVES_ENABLED

enum // node types
{
  TYPE_VERTEX = 0,
  TYPE_EDGE = 1
};

enum // element modes
{
  MODE_TRIANGLE = 0,
  MODE_QUAD = 1
};


const int ANY = -1234;

enum Edges_flag
{
	ANY_BOUNDARY_EDGE = ANY,
	ANY_EDGE = -12345,
	ANY_INNER_EDGE = -123456,
};


// how many bits the order number takes
const int order_bits = 5;
const int order_mask = (1 << order_bits) - 1;


// macros for combining quad horizontal and vertical orders
#define make_quad_order(h_order, v_order) (((v_order) << order_bits) + (h_order))
#define get_h_order(order) ((order) & order_mask)
#define get_v_order(order) ((order) >> order_bits)
extern HERMES2D_API const std::string get_quad_order_str(const int quad_order); ///< Returns string representation of the quad order: used for debugging purposses.

#ifdef COMPLEX
  #include <complex>
  typedef std::complex<double> cplx;
  typedef cplx scalar;
  typedef cplx complex2[2];
#else
  typedef double scalar;
#endif


typedef int int2[2];
typedef int int3[3];
typedef int int4[4];
typedef int int5[5];

typedef double double2[2];
typedef double double3[3];
typedef double double4[4];
typedef double double2x2[2][2];
typedef double double3x2[3][2];

typedef scalar scalar2[2];
typedef scalar scalar3[3];


inline int sqr(int x) { return x*x; }
inline double sqr(double x) { return x*x; }
#ifdef COMPLEX
inline double sqr(cplx x) { return std::norm(x); }
#endif

inline double magn(double x) { return fabs(x); }
#ifdef COMPLEX
inline double magn(cplx x) { return std::abs(x); }
#endif

inline double conj(double a) {  return a; }
#ifdef COMPLEX
inline cplx conj(cplx a) { return std::conj(a); }
#endif

#define is_int(x) ((int) (x) == (x))

/* basic logging functions */
struct HERMES2D_API __h2d_log_info { ///< Info of a log record. Used for output log function.
  const char code;
  const char* log_file;
  const char* src_function;
  const char* src_file;
  const int src_line;
  __h2d_log_info(const char code, const char* log_file, const char* src_function, const char* src_file, const int src_line)
    : code(code), log_file(log_file), src_function(src_function), src_file(src_file), src_line(src_line) {};
};
#define __LOG_INFO(__event) __h2d_log_info(__event, HERMES2D_LOG_FILE, __CURRENT_FUNCTION, __FILE__, __LINE__)
extern HERMES2D_API void __h2d_exit_if(bool cond, int code = -1); ///< Exits with the given code if condition if met.
extern HERMES2D_API bool __h2d_log_message_if(bool cond, const __h2d_log_info& info, const char* msg, ...); ///< Logs a message if condition is met.

/* log file */
#ifdef HERMES2D_REPORT_NO_FILE
#  define HERMES2D_LOG_FILE NULL
#else
# ifdef HERMES2D_REPORT_FILE
#  define HERMES2D_LOG_FILE HERMES2D_REPORT_FILE
# else
#  ifndef HERMES2D_TEST
#    define HERMES2D_LOG_FILE "hermes2d.log" // default filename for a library
#  else
#    define HERMES2D_LOG_FILE "test.log" // default filename for a library test
#  endif
# endif
#endif

/* function name */
#ifdef _WIN32 //Win32
# ifdef __MINGW32__ //MinGW
#   define __CURRENT_FUNCTION __func__
# else //MSVC and other compilers
#   define __CURRENT_FUNCTION __FUNCTION__
# endif
#else //Linux and Mac
# define __CURRENT_FUNCTION __PRETTY_FUNCTION__
#endif

/* event codes */
#define H2D_EC_ERROR 'E' /* errors */
#define H2D_EC_ASSERT 'X' /* asserts */
#define H2D_EC_WARNING 'W' /* warnings */
#define H2D_EC_INFO 'I' /* info about results */
#define H2D_EC_VERBOSE 'V' /* more details for info */
#define H2D_EC_TRACE 'R' /* execution tracing */
#define H2D_EC_TIME 'T' /* time measurements */
#define H2D_EC_DEBUG 'D' /* general debugging messages */

/* error and assert macros */
#define error(...) __h2d_exit_if(__h2d_log_message_if(true, __LOG_INFO(H2D_EC_ERROR), __VA_ARGS__))
#define error_if(__cond, ...) __h2d_exit_if(__h2d_log_message_if(__cond, __LOG_INFO(H2D_EC_ERROR), __VA_ARGS__))
#ifndef NDEBUG
# define assert_msg(__cond, ...) assert(!__h2d_log_message_if(!(__cond), __LOG_INFO(H2D_EC_ASSERT), __VA_ARGS__))
#else
# define assert_msg(__cond, ...)
#endif

/* reporting macros */
#ifdef HERMES2D_REPORT_ALL
# define HERMES2D_REPORT_WARNING
# define HERMES2D_REPORT_INFO
# define HERMES2D_REPORT_VERBOSE
# define HERMES2D_REPORT_TRACE
# define HERMES2D_REPORT_TIME
#endif
#ifdef HERMES2D_REPORT_RUNTIME_CONTROL
# define H2D_RCTR(__var) __var /* reports will be controled also by runtime report control variables */
extern HERMES2D_API bool __h2d_report_warn;
extern HERMES2D_API bool __h2d_report_info;
extern HERMES2D_API bool __h2d_report_verbose;
extern HERMES2D_API bool __h2d_report_trace;
extern HERMES2D_API bool __h2d_report_time;
extern HERMES2D_API bool __h2d_report_debug;
#else
# define H2D_RCTR(__var) true /* reports will be controled strictly by preprocessor directives */
#endif

#if defined(HERMES2D_REPORT_WARNING) || defined(HERMES2D_REPORT_RUNTIME_CONTROL)
# define warn(...) __h2d_log_message_if(true && H2D_RCTR(__h2d_report_warn), __LOG_INFO(H2D_EC_WARNING), __VA_ARGS__)
# define warn_if(__cond, ...) __h2d_log_message_if((__cond) && H2D_RCTR(__h2d_report_warn), __LOG_INFO(H2D_EC_WARNING), __VA_ARGS__)
#else
# define warn(...)
# define warn_if(__cond, ...)
#endif
#if defined(HERMES2D_REPORT_INFO) || defined(HERMES2D_REPORT_RUNTIME_CONTROL)
# define info(...) __h2d_log_message_if(true  && H2D_RCTR(__h2d_report_info), __LOG_INFO(H2D_EC_INFO), __VA_ARGS__)
# define info_if(__cond, ...) __h2d_log_message_if((__cond) && H2D_RCTR(__h2d_report_warn), __LOG_INFO(H2D_EC_INFO), __VA_ARGS__)
#else
# define info(...)
# define info_if(__cond, ...)
#endif
#if defined(HERMES2D_REPORT_VERBOSE) || defined(HERMES2D_REPORT_RUNTIME_CONTROL)
# define verbose(...) __h2d_log_message_if(true && H2D_RCTR(__h2d_report_verbose), __LOG_INFO(H2D_EC_VERBOSE), __VA_ARGS__)
#else
# define verbose(...)
#endif
#if defined(HERMES2D_REPORT_TRACE) || defined(HERMES2D_REPORT_RUNTIME_CONTROL)
# define trace(...) __h2d_log_message_if(true && H2D_RCTR(__h2d_report_trace), __LOG_INFO(H2D_EC_TRACE), __VA_ARGS__)
#else
# define trace(...)
#endif
#if defined(HERMES2D_REPORT_TIME) || defined(HERMES2D_REPORT_RUNTIME_CONTROL)
# define report_time(...) __h2d_log_message_if(true && H2D_RCTR(__h2d_report_time), __LOG_INFO(H2D_EC_TIME), __VA_ARGS__)
#else
# define report_time(...)
#endif
#if defined(_DEBUG) || !defined(NDEBUG) || defined(HERMES2D_REPORT_RUNTIME_CONTROL)
# define debug_log(...) __h2d_log_message_if(true && H2D_RCTR(__h2d_report_debug), __LOG_INFO(H2D_EC_DEBUG), __VA_ARGS__)
#else
# define debug_log(...)
#endif

/* file operations */
void __hermes2d_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream, const __h2d_log_info& err_info);
void __hermes2d_fread(void* ptr, size_t size, size_t nitems, FILE* stream, const __h2d_log_info& err_info);

#define hermes2d_fwrite(ptr, size, nitems, stream) \
      __hermes2d_fwrite((ptr), (size), (nitems), (stream), __LOG_INFO(H2D_EC_ERROR))

#define hermes2d_fread(ptr, size, nitems, stream) \
      __hermes2d_fread((ptr), (size), (nitems), (stream), __LOG_INFO(H2D_EC_ERROR))

/* python support */
extern HERMES2D_API void throw_exception(char *text);

#endif

