#ifndef SIMPLEX_H
#define SIMPLEX_H

double simplex(int m, int n, double **a, double *b, double *c, double *x, double y);

typedef struct simplex_t{
	int m, n; 			
	int* var;
	double** a;		 
	double* b;		
	double* c;
	double* x;
	double y;		
} simplex_t;
#endif
