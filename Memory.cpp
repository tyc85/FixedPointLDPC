#include "Memory.h"
#include "MemoryMacro.h"
#include <stdio.h>
#include <stdlib.h>
#include <new.h>
#include <iostream>

using std::cerr;
using std::cout;
using std::cin;



void nrerror(const char error_text[])
{
	fprintf(stdout,"Numerical Recipes run-time error...\n");
	fprintf(stdout,"%s\n",error_text);
	fprintf(stdout,"...now exiting to system...\n");
	exit(1);
}


int *ivector(int nl,int nh)
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

unsigned char *cvector(int nl,int nh)
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned) (nh-nl+1)*sizeof(unsigned char));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(int nl,int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}


int **imatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,**m;
  
  m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
  if (!m) nrerror("allocation failure 1 in imatrix()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
    if (!m[i]) nrerror("allocation failure 2 in imatrix()");
    m[i] -= ncl;
  }
	//printf("%d * %d matrix allocated at %x", nrh-nrl+1, nch-ncl+1, m);
  return m;
}
double **dmatrix_2(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;
	double *temp;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;
	temp = (double *) malloc((unsigned) (nrh-nrl+1)*(nch-ncl+1)*sizeof(double));
	if (!temp) nrerror("allocation failure 1 in dmatrix()");
	for(i=nrl;i<=nrh;i++) {
		m[i]=&temp[(nrh-nrl+1)*i];
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}
double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}
void free_dvector(double *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
	for(i=nrh;i>=nrl;i--) 
		free( (char*)(m[i]+ncl));
	free((char*) (m+nrl));
}

void free_ivector(int *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
	for(i=nrh;i>=nrl;i--) 
		free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}
