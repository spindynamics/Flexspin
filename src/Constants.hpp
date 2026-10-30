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

// Vacuum magnetic permeability (mu_0), in N/A^2 (or H/m)
const double mu0 = 4.0 * M_PI * 1.0e-7;

// Boltzmann constant (kB), in J/K
const double kB = 1.3806503e-23;

// Gyromagnetic ratio (often used in micromagnetics/OOMMF), in m/(A·s)
const double gamma0 = 2.210173e5;

// Elementary charge of an electron (e), in Coulombs (C)
const double e_charge = -1.60217733e-19;

// Reduced Planck constant (hbar = h / 2π), in J·s
const double h_bar = 6.6262e-34 / (2.0 * M_PI);

