/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#include "defines.h"

#include <iostream>
#include <sstream>

#ifdef DWITHGPERFTOOLS
// using longer path to avoid collisions
#include "install_gperftools/include/profiler.h"
#endif


#ifdef __linux__

#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <libintl.h>
#include <stdio.h>
#include <sys/mman.h>
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle


// This function produces a stack backtrace with demangled function & method names.
// from https://gist.github.com/fmela/591333,
// Copyright (c) 2009-2017, Farooq Mela
void print_callstack(int skip = 3)
{
    void *callstack[128];
    const int nMaxFrames = sizeof(callstack) / sizeof(callstack[0]);
    char buf[1024];
    int nFrames = backtrace(callstack, nMaxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);

    std::ostringstream trace_buf;
    for (int i = skip; i < nFrames; i++) {
        //printf("%s\n", symbols[i]);

        Dl_info info;
        if (dladdr(callstack[i], &info) && info.dli_sname) {
            char *demangled = NULL;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i],
                     status == 0 ? demangled :
                     info.dli_sname == 0 ? symbols[i] : info.dli_sname,
                     (char *)callstack[i] - (char *)info.dli_saddr);
            free(demangled);
        } else {
            snprintf(buf, sizeof(buf), "%-3d %*p %s\n",
                     i, int(2 + sizeof(void*) * 2), callstack[i], symbols[i]);
        }
        trace_buf << buf;
    }
    free(symbols);
    if (nFrames == nMaxFrames)
        trace_buf << "[truncated]\n";

    std::cerr << trace_buf.str() << "\n";
}

extern const char *__progname;

/* This function, when passed a string containing an asserted
   expression, a filename, and a line number, prints a message
   on the standard error stream of the form:
  a.c:10: foobar: Assertion `a == b' failed.
   It then aborts program execution via a call to `abort'.  */

#ifdef FATAL_PREPARE_INCLUDE
# include FATAL_PREPARE_INCLUDE
#endif

// from https://github.com/lattera/glibc/blob/master/assert/assert.c
void
mcell__assert_fail_base (const char *fmt, const char *assertion, const char *file,
        unsigned int line, const char *function)
{
  char *str;

  std::cerr << "-------- ASSERT MESSAGE --------\n\n";

  int total;
  if (__asprintf (&str, fmt,
      __progname, __progname[0] ? ": " : "",
      file, line,
      function ? function : "", function ? ": " : "",
      assertion, &total) >= 0)
    {
      /* Print the message.  */
      fputs(str, stderr);
      fflush (stderr);
      free (str);
    }
  else {
    /* At least print a minimal message.  */
    static const char errstr[] = "Unexpected error during assert.\n";
    fputs(str, stderr);
    fflush (stderr);
  }

  std::cerr << "\n-------- ASSERT STACK TRACE --------\n";
  print_callstack();

  abort ();
}

void
mcell__assert_fail (const char *assertion, const char *file, unsigned int line,
         const char *function)
{
  mcell__assert_fail_base ("%s%s%s:%u: %s%sAssertion `%s' failed.\n%n",
      assertion, file, line, function);
}
#endif

using namespace std;

namespace MCell {

std::ostream & operator<<(std::ostream &out, const vec3_t &a) {
  out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return out;
}


string vec3_t::to_string() const {
  stringstream ss;
  ss << *this;
  return ss.str();
}


void vec3_t::dump(const std::string extra_comment, const std::string ind) const {
  cout << ind << extra_comment << *this << "\n";
}


std::ostream & operator<<(std::ostream &out, const vec2_t &a) {
  out << "(" << a.u << ", " << a.v << ")";
  return out;
}


string vec2_t::to_string() const {
  stringstream ss;
  ss << *this;
  return ss.str();
}


void vec2_t::dump(const std::string extra_comment, const std::string ind) const {
  cout << ind << extra_comment << *this << "\n";
}


void SimulationStats::dump() {
  cout << "Total number of ray-subvolume intersection tests (number of ray_trace calls): " << ray_voxel_tests << "\n";
  cout << "Total number of ray-polygon intersection tests: " << ray_polygon_tests << "\n";
  cout << "Total number of ray-polygon intersections: " << ray_polygon_colls << "\n";
  cout << "Total number of molecule moves between walls: " << mol_moves_between_walls << "\n";
}


void SimulationConfig::dump() {
  cout << "time_unit: \t\t" << time_unit << " [float_t] \t\t\n";
  cout << "length_unit: \t\t" << length_unit << " [float_t] \t\t\n";
  cout << "rx_radius_3d: \t\t" << rx_radius_3d << " [float_t] \t\t\n";
  cout << "partition_edge_length: \t\t" << partition_edge_length << " [float_t] \t\t\n";
  cout << "subpartitions_per_partition_dimension: \t\t" << subpartitions_per_partition_dimension << " [uint] \t\t\n";
  cout << "subpartition_edge_length: \t\t" << subpartition_edge_length << " [float_t] \t\t\n";
}

} // namespace mcell
