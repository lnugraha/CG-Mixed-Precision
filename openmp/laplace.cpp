#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


// set up the nearest neighbor site indices for the laplacian operator on 2d lattice

void siteindex(int* site_ip, int* site_im, int* site_jp, int* site_jm, const int& L)    // lattice size  L x L
{
  int i,j;
  int ip,im,jp,jm;
  int site;

  for (i=0; i<L; i++) {
    im = (i-1+L)%L;
    ip = (i+1)%L;
    for (j=0; j<L; j++) {
      site = i*L + j;
      jp = (j+1)%L;
      jm = (j-1+L)%L;
      site_ip[site] = ip*L + j;
      site_im[site] = im*L + j;
      site_jp[site] = i*L + jp;
      site_jm[site] = i*L + jm;
    }
  }
  return;
}
