#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "simplex.h"
#define EPSILON 1e-6

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
	for(i = 0; i < m + 1; i++) {
		for(j = 0; j < n + 1; j++) {
			p->a[i][j] = a[i][j];
		}
		p->b[i] = b[i];
	}
	for(i = 0; i < n + 1; i++) {
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
	for(i = 0; i < m + 1; i++) {
		q->a[i] = (double*)malloc((q->n+1) * sizeof(double));
	}
	q->b = (double*)malloc((q->m+1) * sizeof(double));
	q->c = (double*)malloc((q->n+1) * sizeof(double));
	q->x = (double*)malloc((q->n+1) * sizeof(double));
	q->min = (double*)malloc(q->n * sizeof(double));
	q->max = (double*)malloc(q->n * sizeof(double));
	for(i = 0; i < n; i++) {
		q->min[i] = p->min[i];
		q->max[i] = p->max[i];
		q->c[i] = p->c[i];
	}
	for(i = 0; i < m + 1; i++) {
		for(j = 0; j < n + 1; j++) {
			q->a[i][j] = p->a[i][j];
		}
		q->b[i] = p->b[i];
	}
	if (ak > 0) {
		if (q->max[k] == INFINITY || bk < q->max[k]) {
			q->max[k] = bk;
		} else if (q->min[k] == -INFINITY || -bk > q->min[k]) {
			q->min[k] = -bk;
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
	}
	return q;
}

bool is_integer(double* xp) {
	double x = *xp;
	double r = round(x);
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
	if (p->z > *zp) {
		*zp = p->z;
		for(int i = 0; i < p->n; i++) {
			x[i] = p->x[i];
		}
		print_queue(h);
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
			if (floor(q->x[h]) < min || floor(q->x[h]) > max) continue;
			q->h = h;
			q->xh = q->x[h];
			return true;
		}
	}
	return false;
}

void succ(node_t* p, queue* h, int n, int m, double** a, double* b, double* c, int k, double ak, double bk, double* zp, double* x) {
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
	node_t* p = initial_node(m, n, a, b, c);
	queue* h = (queue*)malloc(sizeof(queue));
	h->succ = NULL;
	h->node = p;
	double z = -INFINITY;
	p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);
	if (integer(p) || !is_finite(p->z)) {
		z = p->z;
		if (integer(p)) {
			for(int i = 0; i < n; i++) {
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
		return z;
	}
	branch(p, z);
	while(h != NULL) {
		printf("Before pop\n");
		print_queue(h);
		node_t* node = pop(&h);
		printf("After pop\n");
		print_queue(h);
		succ(p, h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
		succ(p, h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);
		free_node(node);
	}
	if (z == -INFINITY) return NAN;
	return z;
}
