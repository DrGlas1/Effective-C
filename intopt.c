#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simplex.h"

double EPSILON = 1e-6;

typedef struct node_t node_t;
struct node_t {
  int m, n, k, h;
  double xh, ak, bk;
  double *min;
  double *max;
  double **a;
  double *b;
  double *x;
  double *c;
  double z;
};

typedef struct set_t set_t;
struct set_t {
  set_t *succ;
  node_t *node;
};

void free_node(node_t *node);

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


void print_set(set_t* h) {
  int l = 0;
  while(h != NULL) {
    printf("%d ", l);
    if (h->node == NULL) {
      printf("NULL\n");
    }
    else {
      printf("%lf\n",h->node->z);
    }
    l++;
    h = h->succ;
  }
}

bool node_in_set(set_t* h, node_t* node) {
	while(h != NULL) {
		if (h->node == node) return true;
		h = h->succ;
	}
	return false;
}

void free_node(node_t *node) {
  free(node->min);
  free(node->max);
  free(node->b);
  free(node->x);
  free(node->c);
  for (int i = 0; i < node->m + 1; i++) {
    free(node->a[i]);
  }
  free(node->a);
  free(node);
}

bool isempty(set_t *h) {
  while (h != NULL) {
    if (h->node != NULL) {
      return false;
    }
    h = h->succ;
  }
  return true;
}

set_t *new_set(node_t *node) {
  set_t *set = calloc(1, sizeof(set_t));

  set->succ = NULL;
  set->node = node;

  return set;
}

bool add(set_t **set, node_t *node) {
  set_t* temp = *set;
  if (temp == NULL) {
  	*set = new_set(node);
	return true;
  }
  while (temp->succ != NULL) {
    temp = temp->succ;
  }
  temp->succ = new_set(node);
  return true;
}

void prune_nodes(set_t** head, double z) {
  set_t* temp = *head;
  set_t *prev = NULL;
  while (temp != NULL) {
      if (temp != NULL && temp->node->z + EPSILON < z) {
        if (prev == NULL) {
              *head = temp->succ;
        } else {
            prev->succ = temp->succ;
      }
      free_node(temp->node);
      free(temp);
	    temp = (prev == NULL) ? *head : prev->succ;
      } else {
          prev = temp;
          temp = temp->succ;
      }
  }
}

node_t* pop(set_t** head) {
    if (*head == NULL) {
        return NULL;
    }

    set_t* temp = *head;
    set_t* prev = NULL;

    while (temp != NULL && temp->node == NULL) {
        prev = temp;
        temp = temp->succ;
    }

    if (temp == NULL) {
        return NULL;
    }

    node_t* node = temp->node;

    if (prev == NULL) {
        *head = temp->succ;
    } else {
        prev->succ = temp->succ;
    }
    free(temp);
    return node;
}

node_t *initial_node(int m, int n, double **a, double *b, double *c) {
  node_t *p = malloc(sizeof(node_t));
  int i, j;
  p->a = calloc((m + 1), sizeof(double *));
  for (i = 0; i < m + 1; i += 1) {
    p->a[i] = calloc(n + 1, sizeof(double));
  }
  p->b = malloc((m + 1)* sizeof(double));
  p->c = malloc((n + 1)* sizeof(double));
  p->x = malloc((n + 1) *sizeof(double));
  p->min = malloc(n * sizeof(double));
  p->max = malloc(n * sizeof(double));
  p->m = m;
  p->n = n;

  for (i = 0; i < m; i += 1) {
    for (j = 0; j < n; j += 1) {
      p->a[i][j] = a[i][j];
    }
  }
  memcpy(p->b, b, sizeof(double) * m);
  p->b[m] = 0;
  memcpy(p->c, c, sizeof(double) * n);
  p->c[n] = 0;

  for (i = 0; i < n; i += 1) {
    p->min[i] = -INFINITY;
    p->max[i] = INFINITY;
  }
  return p;
}

node_t *extend(node_t *p, int m, int n, double **a, double *b, double *c, int k,
               double ak, double bk) {
  node_t *q = calloc(1, sizeof(node_t));
  int i, j;
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
  q->a = (double**)calloc((q->m + 1), sizeof(double));
  for (i = 0; i < q->m + 1; i += 1) {
    q->a[i] = calloc((q->n + 1), sizeof(double));
  }
  q->b = malloc((q->m + 1) * sizeof(double));
  q->c = malloc((q->n + 1) * sizeof(double));
  q->x = malloc((q->n + 1) * sizeof(double));
  q->min = malloc(n * sizeof(double));
  q->max = malloc(n * sizeof(double));
  for (i = 0; i < p->n; i += 1) {
    q->min[i] = p->min[i];
    q->max[i] = p->max[i];
  }
  for (i = 0; i < m; i += 1) {
    for (j = 0; j < n; j += 1) {
      q->a[i][j] = a[i][j];
    }
    q->b[i] = b[i];
  }
  for (i = 0; i < n; i += 1) {
    q->c[i] = c[i];
  }
  if (ak > 0) {
    if (q->max[k] == INFINITY || bk < q->max[k]) {
      q->max[k] = bk;
    }
  } else if (q->min[k] == -INFINITY || -bk > q->min[k]) {
    q->min[k] = -bk;
  }
  for (i = m, j = 0; j < n; j += 1) {
    if (q->min[j] > -INFINITY) {
      q->a[i][j] = -1;
      q->b[i] = -q->min[j];
      i += 1;
    }
    if (q->max[j] < INFINITY) {
      q->a[i][j] = 1;
      q->b[i] = q->max[j];
      i += 1;
    }
  }
  return q;
}

int is_integer(double *xp) {

  double x = *xp;
  double r = lround(x); 
  if (fabs(r - x) < EPSILON) {
    *xp = r;
    return 1;
  } else {
    return 0;
  }
}

int integer(node_t *p) {
  int i;
  for (i = 0; i < p->n; i += 1) {
    if (!is_integer(&p->x[i])) {
      return 0;
    }
  }
  return 1;
}

void bound(node_t* p, double* zp, double* x) {
	if(p == NULL) return;
	if (p->z > *zp) {
		*zp = p->z;
		memcpy(x, p->x, sizeof(double) * p->n);
	}
}

int branch(node_t *q, double z) {
  double min, max;
  if (q->z < z) {
    return 0;
  }

  int h, i;

  for (h = 0; h < q->n; h += 1) {
    if (!is_integer(&q->x[h])) {
      if (q->min[h] == -INFINITY) {
        min = 0;
      } else {
        min = q->min[h];
      }
      max = q->max[h];

      if (floor(q->x[h]) < min || ceil(q->x[h]) > max) {
        continue;
      }

      q->h = h;
      q->xh = q->x[h];

      return 1;
    }
  }
  return 0;
}

void succ(node_t *p, set_t **h, int m, int n, double **a, double *b, double *c,
          int k, double ak, double bk, double *zp, double *x) {
  node_t *q = extend(p, m, n, a, b, c, k, ak, bk);
  if (q == NULL) {
    return;
  }
  q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
  if (isfinite(q->z)) {
    if (integer(q)) {
      bound(q, zp, x);
    } else if (branch(q, *zp)) {
      add(h, q);
      return;
    }
  }
  free_node(q);
}

double intopt(int m, int n, double **a, double *b, double *c, double *x) {
  node_t *p = initial_node(m, n, a, b, c);

  set_t *h = new_set(p);

  double z = -INFINITY; 
  p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);
  if (integer(p) || !isfinite(p->z)) {
    z = p->z;
    if (integer(p)) {
      memcpy(x, p->x, sizeof(double) * p->n);
    }

    while (h != NULL) {
      if (h->node != NULL) {
        free_node(h->node);
      }
      set_t *temp = h->succ;
      free(h);
      h = temp;
    }

    return z;
  }
  branch(p, z);
  double prev_max;
  while (h != NULL) {
    p = pop(&h);
    prev_max = z;
    succ(p, &h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
    succ(p, &h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);
    if (prev_max < z) prune_nodes(&h, z);
    free_node(p);
  }

  while (h != NULL) {
    if (h->node != NULL) {
      free_node(h->node);
    }
    set_t *temp = h->succ;
    free(h);
    h = temp;
  }

  if (z == -INFINITY) {
    return NAN;
  } else {
    return z;
  }
}

