#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#define EPSILON 1e-6

typedef struct simplex_t{
	int m, n; 			
	int* var;
	double** a;		 
	double* b;		
	double* c;
	double* x;
	double y;		
} simplex_t;

typedef struct node_t {
	int m, n, k, h;
	double xh, ak, bk;
	double* min;
	double* max;
	double** a;
	double* b;
	double* x;
	double* c;
	double z;
} node_t;

typedef struct queue {
	struct queue* succ;
	node_t* node;
} queue;

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);


void free_node(node_t* node) {
	free(node->min);
	free(node->max);
	free(node->b);
	free(node->x);
	free(node->c);
	for(int i = 0; i < node->m + 1; i++) {
		free(node->a[i]);
	}
	free(node->a);
	free(node);
}

queue* push(queue* head, node_t* node) {
	queue* temp = head;
	queue* q = (queue*)malloc(sizeof(queue));
	q->node = node;
	q->succ = NULL;
	if (temp == NULL) {
		temp = q;
		return temp;
	}
	while(temp->succ != NULL) {
		temp = temp->succ;
	}
	temp->succ = q;
	return head;
}

node_t* pop(queue** head) {
	queue* node = *head;
	queue* prev = NULL;
	while (node->succ != NULL) {
		prev = node;
		node = node->succ;
	}
	if (prev == NULL) {
		*head = NULL;  
	} else {
		prev->succ = NULL;  
	}
	node_t* popped_node = node->node;

	free(node);
	return popped_node;
}

void delete(queue** head, node_t* node) {
	queue* temp = *head, *prev;
	if (temp != NULL && temp->node == node) {
		*head = temp->succ;
		free_node(temp->node);
		free(temp);
		return;
	}
	while(temp != NULL && temp->node != node) {
		prev = temp;
		temp = temp->succ;
	}

	if (temp == NULL) return;

	prev->succ = temp->succ;
	free_node(temp->node);
	free(temp);
}

void prune_nodes(queue** head, double z) {
	queue* temp = *head, *prev;
	if (temp != NULL && temp->node->z < z) {
		*head = temp->succ;
		free_node(temp->node);
		free(temp);
		return;
	}
	while(temp != NULL) {
		if (temp->node->z < z) {
			prev->succ = temp->succ;
			free_node(temp->node);
			free(temp);
			temp = prev->succ;
		} else {
			prev = temp;
			temp = temp->succ;
		}
		
	}
}

void print_length(queue* h) {
	int l = 0;
	return;
	while(h != NULL) {
		l++;
		h = h->succ;
	}
	printf("%d\n",l);
}

void print_queue(queue* h) {
	while(h != NULL && h->node != NULL) {
		printf("%lf ", h->node->z);
		h = h->succ;
	}
	if (h != NULL && h->node == NULL) {
		printf("\nMADNESS");
	}
	printf("\n");
}

void print(simplex_t* s) {
	double** a = s->a;
	double* b = s->b;
	double* c = s->c;
	double* x = s->x;
	int* var = s->var;
	int m = s->m;
	int n = s->n;
	int i, j;
	printf("-----------------------------------------------\n");
	printf("maximize: ");
	for(i = 0; i < n; i++) {
		printf("%.2fx_%d ", c[i],var[i]);
		printf("+ ");
	}
	
	printf("%.2f\n", s->y);
	printf("\nsubject to: \n");
	for(i = 0; i < m; i++) {
		printf("x_%d = %.2f -(", var[i+n], b[i]);
		for(j = 0; j < n; j++) {
			printf("%.2fx_%d", a[i][j], var[j]);
			if (j != n - 1) {
				printf(" + ");
			}
		}
		printf(")\n");
	}
	printf("\n");
}

int init(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
	int i, k;
	s->m = m;
    	s->n = n;
    	s->var = var;
    	s->a = a;
    	s->b = b;
    	s->c = c;
    	s->x = x;
    	s->y = y;
	printf("heloooeloo %d %d\n",m,n);
	if (s->var == NULL) {	
		s->var = malloc((m+n+1) * sizeof(int));
		for(i = 0; i < m + n; i++) {
			s->var[i] = i;
		}
	}
	for (k = 0,i = 1; i < m; i++) {
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
	int i,j,t;
	t = s->var[col];
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
			a[row][i] = a[row][i] / a[row][col];
		}
	}
	b[row] = b[row] / a[row][col];
	a[row][col] = 1 / a[row][col];
	printf("pivot row=%d col=%d\n", row, col);
	print(s);
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
	s->x = malloc((m+n) * sizeof(double));
	s->c = malloc(n * sizeof(double));
	s->c[n-1] = -1;
	s->n = n;
	printf("S after prepare\n");
	print(s);
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
	printf("Hello\n");
	print(s);
	if (b[k] >= 0) {
		return true;
	}
	printf("b[%d]=%lf\n", k, b[k]);
	printf("initial\n");
	print(s);
	prepare(s,k);
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
		for(k = 0; k < n; k++) {
			if(fabs(s->a[i-n][k]) > fabs(s->a[i-n][j])) {
				j = k;
			}
		}
		pivot(s, i-n, j);
		i = j;	
	}
       	if (i < n - 1) {
		k = s->var[i];
		s->var[i] = s->var[n-1];
		s->var[n-1] = k;
		for(k = 0; k < m; k++) {
			w = s->a[k][n-1];
			s->a[k][n-1] = s->a[k][i];
			s->a[k][i] = w;
		}
	}
	free(s->c);
	s->c = c;
	s->y = y;
	for(k = n-1; k < n + m - 1; k++) {
		s->var[k] = s->var[k+1];
	}
	n = --s->n;
	//Look at this
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
		s->y += s->c[k] * s->b[j];
		for(i = 0; i < n; i++) {
			t[i] -= s->c[k] * s->a[j][i];
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
	printf("xsimplex\n");
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


node_t* initial_node(int m, int n, double** a, double* b, double* c) {
	node_t* p = (node_t*)malloc(sizeof(node_t));
	p->a = (double**)malloc((m+1) * sizeof(double));
	int i, j;
	for(i = 0; i < m + 1; i++) {
		p->a[i] = (double*)malloc((n+1) * sizeof(double));
	}
	p->b = (double*)malloc((m+1) * sizeof(double));
	p->c = (double*)malloc((n+1) * sizeof(double));
	p->x = (double*)malloc((n+1) * sizeof(double));
	p->min = (double*)malloc(n * sizeof(double));
	p->max = (double*)malloc(n * sizeof(double));
	p->m = m;
	p->n = n;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n + 1; j++) {
			p->a[i][j] = a[i][j];
		}
		p->b[i] = b[i];
	}
	for(i = 0; i < n; i++) {
		p->c[i] = c[i];
	}
	for(int i = 0; i < n; i++) { 
		p->min[i] = -INFINITY;
		p->max[i] = INFINITY;
	}
	return p;

}

node_t* extend(node_t* p, int m, int n, double** a, double* b, double* c, int k, double ak, double bk) {
	node_t* q = (node_t*)malloc(sizeof(node_t));
	int i,j;
	q->k = k;
	q->ak = ak;
	q->bk = bk;
	if (ak > 0 && p->max[k] < INFINITY) {
		q->m = p->m;
	} else if (ak < 0 && p->min[k] > 0) {
		q->m = p->m;
	} else {
		q->m = p->m + 1;
	}
	q->n = p->n;
	q->h = -1;
	q->a = (double**)malloc((q->m+1) * sizeof(double));
	for(i = 0; i < q->m + 1; i++) {
		q->a[i] = (double*)malloc((q->n+1) * sizeof(double));
	}
	q->b = (double*)malloc((q->m+1) * sizeof(double));
	q->c = (double*)malloc((q->n+1) * sizeof(double));
	q->x = (double*)malloc((q->n+1) * sizeof(double));
	q->min = (double*)malloc(n * sizeof(double));
	q->max = (double*)malloc(n * sizeof(double));
	for(i = 0; i < p->n; i++) {
		q->min[i] = p->min[i];
		q->max[i] = p->max[i];
	}
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			q->a[i][j] = a[i][j];
		}
		q->b[i] = b[i];
	}
	for(i = 0; i < n; i++) {
		q->c[i] = c[i];
	}
	if (ak > 0) {
		if (q->max[k] == INFINITY || bk < q->max[k]) {
			q->max[k] = bk;
		} else if (q->min[k] == -INFINITY || -bk > q->min[k]) {
			q->min[k] = -bk;
		}
		
	}
	for(i = m, j = 0; j < n; j++) {
		if (q->min[j] > -INFINITY) {
			q->a[i][j] = -1;
			q->b[i] = -q->min[j];
			i++;
		}
		if (q->max[j] < INFINITY) {
			q->a[i][j] = 1;
			q->b[i] = q->max[j];
			i++;
		}
	}
	return q;
}

bool is_integer(double* xp) {
	double x = *xp;
	double r = lround(x);
	if (fabs(r - x) < EPSILON) {
		*xp = r;
		return true;
	}
	return false;
}

bool integer(node_t* p) {
	int i;
	for(i = 0; i < p->n; i++) {
		if (!is_integer(&p->x[i])) {
			return 0;
		}
	}
	return 1;
}

void bound(node_t* p, queue* h, double* zp, double* x) {
	if(p == NULL) return; //TODO - remove
	if (p->z > *zp) {
		*zp = p->z;
		memcpy(x, p->x, sizeof(double) * p->n);
		
		print_queue(h);
		printf("Pruning nodes\n");
		prune_nodes(&h, p->z);
		printf("%lf\n", p->z);
		print_queue(h);
	}
}

bool is_finite(double x) {
	return x != NAN && fabs(x) != INFINITY;
}

bool branch(node_t* q, double z) {
	double min, max;
	if (q->z < z) return false;
	for(int h = 0; h < q->n; h++) {
		if(!is_integer(&q->x[h])) {
			if (q->min[h] == -INFINITY) {
				min = 0;
			} else {
				min = q->min[h];
			}
			max = q->max[h];
			if (floor(q->x[h]) < min || ceil(q->x[h]) > max) continue;
			q->h = h;
			q->xh = q->x[h];
			return true;
		}
	}
	return false;
}

void succ(node_t* p, queue* h, int m, int n, double** a, double* b, double* c, int k, double ak, double bk, double* zp, double* x) {
	node_t* q = extend(p, m, n, a, b, c, k, ak, bk);
	if (q == NULL) return;
	q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
	if (is_finite(q->z)) {
		if (integer(q)) {
			bound(q, h, zp, x);
		} else if (branch(q, *zp)) {
			push(h, q);
			return;
		}
	}
	free_node(q);
}

double intopt(int m, int n, double** a, double* b, double* c, double* x) {
	node_t* p = initial_node(m, n, a, b, c), *node;
	queue* h = (queue*)malloc(sizeof(queue));
	h->succ = NULL;
	h->node = p;
	double z = -INFINITY;
	p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);
	if (integer(p) || !is_finite(p->z)) {
		z = p->z;
		if (integer(p)) {
			for(int i = 0; i < p->n; i++) {
				x[i] = p->x[i];
			}
		}
		queue* temp = h, *prev;
		while(temp != NULL) {
			free_node(h->node);
			prev = temp;
			temp = temp->succ;
			free(prev);
		}
		free_node(p);
		return z;
	}
	branch(p, z);
	while(h != NULL) {
		printf("Before pop\n");
		print_queue(h);
		node = pop(&h);
		printf("After pop\n");
		print_queue(h);
		succ(node, h, m, n, a, b, c, node->h, 1, floor(node->xh), &z, x);
		succ(node, h, m, n, a, b, c, node->h, -1, -ceil(node->xh), &z, x);
		free_node(node);
	}
	if (z == -INFINITY) return NAN;
	return z;
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
		printf("%.2fx_%d ", c[i], i);
		if (i != n - 1) {
			printf("+ ");
		}
	}
	printf("\n");

	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			printf("%.2fx_%d ", a[i][j], i);
			if (j != n - 1) {
				printf("+ ");
			} else {
				printf("\u2264 %lf", b[i]);
			}
		}
		printf("\n");
	}
	printf("\n");
	double res = 0;
	res = intopt(m, n, a, b, c, x); 
	printf("%lf\n", res);
	
	for(i = 0; i < m; i++) {
		free(a[i]);
	}

	free(a);
	free(b);
	free(c);
	free(x);
	return 0;
}
