#include "Laplace.hpp"

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

void siteindex(const SiteIndices& SITES)
{
  int ip, im, jp, jm, site;
  for (int i=0; i<SITES.L; i++){
    im = (i-1+SITES.L)%SITES.L;
    ip = (i+1)%SITES.L;
    for (int j=0; j<SITES.L; j++){
      site = i*SITES.L + j;
      jp = (j+1)%SITES.L;
      jm = (j-1+SITES.L)%SITES.L;
      SITES.site_ip[site] = ip*SITES.L + j;
      SITES.site_im[site] = im*SITES.L + j;
      SITES.site_jp[site] = i*SITES.L + jp;
      SITES.site_jm[site] = i*SITES.L + jm;
    } // END-FOR j
  } // END-FOR i
}
