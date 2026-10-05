#pragma once

#include <cmath>

// mu0 and h_bar both need pi. M_PI is a POSIX extension rather than standard
// C++ (MSVC hides it behind _USE_MATH_DEFINES), so define it here rather than
// leaving each caller to guess. g++ already defines _GNU_SOURCE for C++, so on
// glibc this is a no-op and the fallback only fires where M_PI is genuinely
// absent.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double mu0 = 4.0 * M_PI * 1.0e-7;
const double kB = 1.3806503e-23;
const double gamma0 = 2.210173e5;
const double e_charge = -1.60217733e-19;
const double h_bar = 6.6262e-34 / (2.0 * M_PI);
