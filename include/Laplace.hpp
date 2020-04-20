#ifndef LAPLACE_HPP_ 
#define LAPLACE_HPP_

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "MatrixOperations.hpp"

void siteindex(
  unsigned int* site_ip, unsigned int* site_im, 
  unsigned int* site_jp, unsigned int* site_jm, const unsigned int& L);

void siteindex(const SiteIndices& SITES, const unsigned int& L);

#endif
