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

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);

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
	for(i = 0; i < n; i++) {
		printf("%lfx_%d ", c[i],var[i]);
		printf("+ ");
	}
	printf("%lf\n", s->y);

	for(i = 0; i < m; i++) {
		printf("x_%d = %lf -(", var[i+n], b[i]);
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
	free(s->x);
	free(s->c);
	s->x = malloc(sizeof(double)*(m+n));
	s->c = malloc(sizeof(double)*n);
	s->n = n;
	pivot(s, k, n - 1);
}

int select_nonbasic(simplex_t* s) {
	for(int i = 0; i < s->n; i++) {
		if (s->c[i] > EPSILON) {
			return i;
		}
	}
	return -1;
}

bool initial(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
	int i, j, k;
	double w;
	k = init(s, m, n, a, b, c, x, y, var);
	if (b[k] > 0) {
		return true;
	}
	printf("Prepare s");
	prepare(s,k);
	print(s);
	n = s->n;
	s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);
	for(i = 0; i < m + n; i++) {
		if(s->var[i] == m+n-1) {
			if(fabs(s->x[i]) > EPSILON) {
				free(s->x);
				free(s->c);
				return false;
			} else {
				break;
			
			}
		}
	}		
	if (i >= n) {
		j = 0;
		for(j = 0; k < n; k++) {
			if(fabs(s->a[i-n][k]) > fabs(s->a[i-n][j])) {
				j = k;
			}
		}
		pivot(s, i-n, j);
		i = j;	
	} else if (i < n - 1) {
		k = s->var[i];
		s->var[i] = s->var[n-1];
		s->var[n-1] = k;
	}
	free(s->c);
	s->c = c;
	s->y = y;
	for(k = n-1; k < n + m - 1; k++) {
		s->var[k] = s->var[k+1];
	}
	n = s->n;
	s->n--;
	double t[n];
        for(k = 0; k < n; k++) {
		for(j = 0; j < n; j++) {
			if(k == s->var[j]) {
				t[j] = t[j] + s->c[k];
				goto next_k;
			}
		}
		for(j = 0; j < m; j++) {
			if(s->var[n+j] == k) {
				break;
			}	
		}
		printf("Questionable part");
		s->y += s->c[k] * s->b[j];
		for(i = 0; i < n; i++) {
			t[i] -= s->c[k] * s->a[i][j];
		}
		next_k:;
	}
	for(i = 0; i < n; i++) {
		s->c[i] = t[i];
	}
	free(s->x);
	return true;
}


double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h) {
	simplex_t* s = malloc(sizeof(simplex_t));
	int i, row, col;
	if (!initial(s, m, n, a, b, c, x, y, var)) {
		free(s->var);
		return NAN;
	}
	print(s);
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

int main() {
	int m, n;

	scanf("%d %d", &m, &n);

	printf("%d, %d\n", m, n);

	double** a = (double**)malloc(sizeof(double) * m);
	double* b = (double*)malloc(sizeof(double) * m);
	double* c = (double*)malloc(sizeof(double) * n);
	double* x =(double*)malloc(sizeof(double) * (n+1)); 
	
	int i, j;
	
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
