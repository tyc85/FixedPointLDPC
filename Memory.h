#ifndef MEMORY_H
#define MEMORY_H

void free_ivector(int*, int, int);
void free_imatrix(int**, int, int, int, int);
void free_dmatrix(double **,int, int, int, int);
void free_dvector(double *, int, int);

void nrerror(const char*);
int *ivector(int, int);
int **imatrix(int, int, int, int);
unsigned char *cvector(int, int);
double *dvector(int,int);
double **dmatrix(int,int,int,int);
double **dmatrix_2(int,int,int,int);

#endif