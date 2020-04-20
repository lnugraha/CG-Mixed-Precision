#include "Laplace.hpp"


void siteindex(
  unsigned int* site_ip, unsigned int* site_im, 
  unsigned int* site_jp, unsigned int* site_jm, 
  const unsigned int& L)    // lattice size  L x L
{
  int ip, im, jp, jm, site;
  for (int i=0; i<L; i++){
    im = (i-1+L)%L;
    ip = (i+1)%L;
    for (int j=0; j<L; j++){
      site = i*L + j;
      jp = (j+1)%L;
      jm = (j-1+L)%L;
      site_ip[site] = ip*L + j;
      site_im[site] = im*L + j;
      site_jp[site] = i*L + jp;
      site_jm[site] = i*L + jm;
    } // END-FOR j
  } // END-FOR i
}

void siteindex(const SiteIndices& SITES, const unsigned int& L){
  int ip, im, jp, jm, site;
  for (int i=0; i<L; i++){
    im = (i-1+L)%L;
    ip = (i+1)%L;
    for (int j=0; j<L; j++){
      site = i*L + j;
      jp = (j+1)%L;
      jm = (j-1+L)%L;
      SITES.site_ip[site] = ip*L + j;
      SITES.site_im[site] = im*L + j;
      SITES.site_jp[site] = i*L + jp;
      SITES.site_jm[site] = i*L + jm;
    } // END-FOR j
  } // END-FOR i


}
