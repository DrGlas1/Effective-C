#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#define EPSILON 1e-8

typedef struct simplex_t{
	int m; 			
	int n; 			
	int* var;
	double** a;		 
	double* b;		
	double* c;
	double* x;
	double y;		
} simplex_t;

void print(simplex_t* s);

int init(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
	int i;
	int k = 0;
	s->m = m;
    	s->n = n;
    	s->var = var;
    	s->a = a;
    	s->b = b;
    	s->c = c;
    	s->x = x;
    	s->y = y;
	if (s->var == NULL) {	
		s->var = malloc(sizeof(int) * (m + n + 1));
		for(i = 0; i < m + n; i++) {
			s->var[i] = i;
		}
	}
	for (i = 1; i < m; i++) {
		if (b[i] < b[k]) {
			k = i;
		}
	} 
	return k;

}

bool initial(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
	int i, j, k;
	k = init(s, m, n, a, b, c, x, y, var);
	if (b[k] > 0) {
		return true;
	}
	// TODO - handle this case
	return false;
}



void prepare(simplex_t* s, int k) {
	int m = s->m;
	int n = s->n;
	int i;
	for(i = m + n; i > n; i--) {
		s->var[i] = s->var[i - 1];
	}
	s->var[n] = m + n;
	n++;
	for(i = 0; i < m; i++) {
		s->a[i][n - 1] = -1;
	}

}

int select_nonbasic(simplex_t* s) {
	for(int i = 0; i < s->n; i++) {
		if (s->c[i] > EPSILON) {
			return i;
		}
}
	return -1;
}


void pivot(simplex_t* s, int row, int col) {
	double** a = s->a;
	double* b = s->b;
	double* c = s->c;
	int m = s->m;
	int n = s->n;
	int t = s->var[col];
	int i;
	int j;
	s->var[col] = s->var[n + row];
	s->var[n + row] = t;
	s->y += c[col] * b[row] / a[row][col];
	for(i = 0; i < n; i++) {
		if (i != col) {
			c[i] -= c[col] * a[row][i] / a[row][col];
		}
	}
	c[col] = -c[col] / a[row][col];
	for(i = 0; i < m; i++) {
		if (i != row) {
			b[i] -= a[i][col] * b[row] / a[row][col];
		}
	}
	for(i = 0; i < m; i++) {
		if (i != row) {
			for(j = 0; j < n; j++) {
				if (j != col) {
					a[i][j] -= a[i][col] * a[row][j] / a[row][col];
				}
			}
		}
	}
	for(i = 0; i < m; i++) {
		if (i != row) {
			a[i][col] = -a[i][col] / a[row][col];
		}
	}
	for(i = 0; i < n; i++) {
		if (i != col) {
			a[row][i] = -a[row][i] / a[row][col];
		}
	}
	b[row] = b[row] / a[row][col];
	a[row][col] = 1 / a[row][col];
}


double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h) {
	simplex_t* s = malloc(sizeof(simplex_t));
	int i, row, col;
	if (!initial(s, m, n, a, b, c, x, y, var)) {
		free(s->var);
		return NAN;
	}
	while((col = select_nonbasic(s)) >= 0) {
		row = -1;
		for(i = 0; i < m; i++) {
			if (a[i][col] > EPSILON && (row < 0 || b[i] / a[i][col] < b[row]/a[row][col])) {
				row = i;
			}
		}
		if (row < 0) {
			free(s->var);
			return INFINITY;
		}
		pivot(s, row, col);
		print(s);
	}
	if (h == 0) {
		for(i = 0; i < n; i++) {
			if (s->var[i] < n) {
				x[s->var[i]] = 0;
			}
		}
		for(i = 0; i < m; i++) {
			if (s->var[n + i] < n) {
				x[s->var[n + i]] = s->b[i];
			}
		}
		free(s->var);
	} else {
		for(i = 0; i < n; i++) {
			x[i] = 0;
		}
		for(i = n; i < n+m; i++) {
			x[i] = s->b[i-n];
		}
	}
	double res = s->y;
	free(s);
	return res;
}

double simplex(int m, int n, double** a, double* b, double* c, double* x, double y) {
	return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

void print(simplex_t* s) {
	double** a = s->a;
	double* b = s->b;
	double* c = s->c;
	int* var = s->var;
	int m = s->m;
	int n = s->n;
	int i, j;
	printf("-----------------------------------------------\n");
	printf("max z = ");
	for(i = 0; i < m; i++) {
		printf("%lfx_%d ", c[i],var[i]);
		printf("+ ");
	}
	printf("%lf\n", s->y);

	for(i = 0; i < m; i++) {
		printf("x_%d = -(", var[i+n]);
		for(j = 0; j < n; j++) {
			printf("%lfx_%d", a[i][j], var[j]);
			if (j != n - 1) {
				printf(" + ");
			}
		}
		printf(")\n");
	}
	printf("\n");
}

int main() {
	int m, n;

	scanf("%d %d", &m, &n);

	printf("%d, %d\n", m, n);

	double** a;
	double b[m];
	double c[n];
	double x[n + 1];
	
	int i, j;
	
	a = (double**)malloc(sizeof(double) * m);
	for(i = 0; i < m; i++) {
		a[i] = (double*)malloc(sizeof(double) * (n+1));
	}
		
	for(i = 0; i < n; i++) {
		scanf("%lf", &c[i]);
	}

	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			scanf("%lf", &a[i][j]);
		}
	}

	for(i = 0; i < m; i++) {
		scanf("%lf", &b[i]);
	}

	printf("max z = ");
	for(i = 0; i < n; i++) {
		printf("%lfx_%d ", c[i], i);
		if (i != n - 1) {
			printf("+ ");
		}
	}
	printf("\n");

	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			printf("%lfx_%d ", a[i][j], i);
			if (j != n - 1) {
				printf("+ ");
			} else {
				printf("\u2264 %lf", b[i]);
			}
		}
		printf("\n");
	}
	double res = simplex(m, n, a, b, c, x, 0); 
	printf("%lf\n", res);
	
	for(i = 0; i < m; i++) {
		free(a[i]);
	}

	free(a);

	return 0;
}
