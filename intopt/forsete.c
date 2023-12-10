#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EPSILON 1e-6
#define ALPHA 0.5
#define INIT_P_COST 1.0

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
  double xh, ak, bk, z;
  double *min;
  double *max;
  double **a;
  double *b;
  double *x;
  double *c;
} node_t;

typedef struct p_queue p_queue;
struct p_queue {
  p_queue *succ;
  node_t *node;
};

double temp, recip, x_recip;

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

p_queue *new_set(node_t *node) {
  p_queue *set = malloc(sizeof(p_queue));
  set->succ = NULL;
  set->node = node;
  return set;
}

void add(p_queue **set, node_t *node) {
  p_queue* n = new_set(node);
  if (*set == NULL || (*set)->node->z < node->z) {
    n->succ = *set;
    *set = n;
    return;
  }
  p_queue* temp = *set;
  while (temp->succ != NULL && temp->succ->node->z > node->z) {
    temp = temp->succ;
  }
  n->succ = temp->succ;
  temp->succ = n;
}

void prune_nodes(p_queue** head, double z) {
  p_queue* temp = *head;
  p_queue *prev = NULL;
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

node_t* pop(p_queue** head) {
  if (*head == NULL) {
      return NULL;
  }
  p_queue* temp = *head;
  p_queue* prev = NULL;
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



bool initial(simplex_t *s, int m, int n, double **a, double *b, double *c,
            double *x, double y, int *var);

int init(simplex_t *s, int m, int n, double **a, double *b, double *c,
         double *x, double y, int *var) {
  int i, k;
  s->m = m;
  s->n = n;
  s->var = var;
  s->a = a;
  s->b = b;
  s->x = x;
  s->c = c;
  s->y = y;

  if (s->var == NULL) {
    s->var = malloc((m + n + 1) * sizeof(int));
    for (i = 0; i < m + n; i += 1) {
      s->var[i] = i;
    }
  }
  for (k = 0, i = 1; i < m; i += 1) {
    if (b[i] < b[k]) {
      k = i;
    }
  }
  return k;
}

int select_nonbasic(simplex_t *s) {
  int i, max = -1;
  temp = EPSILON;
  for (i = 0; i < s->n; i++) {
    if(s->c[i] > temp) {
      max = i, temp = s->c[i];
    }
  }
  return max;
}

void pivot(simplex_t *s, int row, int col) {
  double **a = s->a;
  double *b = s->b;
  double *c = s->c;

  recip = 1 / a[row][col];
  x_recip = c[col] * recip;


  int m = s->m;
  int n = s->n;
  int i, j, t;
  t = s->var[col];
  s->var[col] = s->var[n + row];
  s->var[n + row] = t;
  s->y = s->y + b[row] * x_recip;

  temp = c[col];
  for (i = 0; i < n; i += 1) {
    c[i]-= a[row][i] * x_recip;
  }
  c[col] = temp;
  c[col] = -x_recip;
  
  x_recip = recip * b[row];

  temp = b[row];
  for (i = 0; i < m; i += 1) {
    b[i] -= a[i][col] * x_recip; 
    if (i != row) {
      for (j = 0; j < n; j += 1) {
        if (j != col) {
          a[i][j] -= a[i][col] * a[row][j] * recip;
        }
      }
    }
  }
  b[row] = temp;
  for (i = 0; i < m; i += 1) {
    a[i][col] *= -recip;
  }
  for (i = 0; i < n; i += 1) {
    a[row][i] *= recip;
  }

  b[row] = x_recip;
  a[row][col] = recip;
}

void prepare(simplex_t *s, int k) {
  int m = s->m;
  int n = s->n;
  int i;

  for (i = m + n; i > n; i -= 1) {
    s->var[i] = s->var[i - 1];
  }
  s->var[n] = m + n;
  n += 1;

  for (i = 0; i < m; i += 1) {
    s->a[i][n - 1] = -1;
  }
  s->x = malloc((m + n) * sizeof(double));
  s->c = calloc(n, sizeof(double));
  s->c[n - 1] = -1;
  s->n = n;
  pivot(s, k, n - 1);
}

double xsimplex(int m, int n, double **a, double *b, double *c, double *x,
                double y, int *var, int h) {
  simplex_t s;
  int i, row, col;
  if (!initial(&s, m, n, a, b, c, x, y, var)) {
    free(s.var);
    return NAN;
  }

  while ((col = select_nonbasic(&s)) >= 0) {
    row = -1;
    for (i = 0; i < m; i += 1) {
      if (a[i][col] > EPSILON &&
          (row < 0 || b[i] * a[row][col] < b[row] * a[i][col])) {
        row = i;
      }
    }
    if (row < 0) {
      free(s.var);
      return INFINITY;
    }
    pivot(&s, row, col);
  }

  if (h == 0) {
    for (i = 0; i < n; i += 1) {
      if (s.var[i] < n) {
        x[s.var[i]] = 0;
      }
    }
    for (i = 0; i < m; i += 1) {
      if (s.var[n + i] < n) {
        x[s.var[n + i]] = s.b[i];
      }
    }
    free(s.var);
  } else {
    for (i = 0; i < n; i += 1) {
      x[i] = 0;
    }
    for (i = n; i < n + m; i += 1) {
      x[i] = s.b[i - n];
    }
  }
  return s.y;
}

bool initial(simplex_t *s, int m, int n, double **a, double *b, double *c,
            double *x, double y, int *var) {
  int i, j, k;
  double w;
  k = init(s, m, n, a, b, c, x, y, var);
  if (s->b[k] >= 0) {
    return 1;
  }
  prepare(s, k);
  n = s->n;
  s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

  for (i = 0; i < m + n; i += 1) {
    if (s->var[i] == m + n - 1) {
      if (fabs(s->x[i]) > EPSILON) {
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
    for (k = 0; k < n; k += 1) {
      if (fabs(s->a[i - n][k]) > fabs(s->a[i - n][j])) {
        j = k;
      }
    }
    pivot(s, i - n, j);
    i = j;
  }

  if (i < n - 1) {
    k = s->var[i];
    s->var[i] = s->var[n - 1];
    s->var[n - 1] = k;
    for (k = 0; k < m; k += 1) {
      w = s->a[k][n - 1];
      s->a[k][n - 1] = s->a[k][i];
      s->a[k][i] = w;
    }
  }

  free(s->c);
  s->c = c;
  s->y = y;
  for (k = n - 1; k < n + m - 1; k += 1) {
    s->var[k] = s->var[k + 1];
  }

  n = --s->n;
  double *t = calloc(n, sizeof(double));

  for (k = 0; k < n; k += 1) {
    for (j = 0; j < n; j += 1) {
      if (k == s->var[j]) {
        t[j] = t[j] + s->c[k];
        goto next_k;
      }
    }
    for (j = 0; j < m; j += 1) {
      if (s->var[n + j] == k) {
        break;
      }
    }
    s->y = s->y + s->c[k] * s->b[j];
    for (i = 0; i < n; i += 1) {
      t[i] = t[i] - s->c[k] * s->a[j][i];
    }
  next_k:;
  }
  for (i = 0; i < n; i += 1) {
    s->c[i] = t[i];
  }
  free(t);
  free(s->x);
  return true;
}

double simplex(int m, int n, double **a, double *b, double *c, double *x,
               double y) {
  return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

node_t *initial_node(int m, int n, double **a, double *b, double *c) {
  node_t *p = malloc(sizeof(node_t));
  int i, j;
  p->a =  (double**)malloc((m + 1) * sizeof(double *));
  for (i = 0; i < m + 1; i += 1) {
    p->a[i] =  (double*)malloc((n + 1) * sizeof(double));
  }
  p->b = (double*)malloc((m + 1)* sizeof(double));
  p->c = (double*) malloc((n + 1)* sizeof(double));
  p->x = (double*)malloc((n + 1) *sizeof(double));
  p->min = (double*)malloc(n * sizeof(double));
  p->max = (double*)malloc(n * sizeof(double));
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
  node_t *q = malloc(sizeof(node_t));
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
  q->a = (double**)malloc((q->m + 1) * sizeof(double));
  for (i = 0; i < q->m + 1; i += 1) {
    q->a[i] = (double*)calloc((q->n + 1), sizeof(double));
  }
  q->b =(double*) malloc((q->m + 1) * sizeof(double));
  q->c = (double*)malloc((q->n + 1) * sizeof(double));
  q->x = (double*)malloc((q->n + 1) * sizeof(double));
  q->min = (double*)malloc(n * sizeof(double));
  q->max = (double*)malloc(n * sizeof(double));
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
  if (q->z < z) {
    return 0;
  }
  double min, max;
  int h, i;
  for(h = 0; h < q->n; h += 1) {
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

void update_p_cost(double* p_up, double* p_down, int k, double delta) {
  p_up[k] = (1-ALPHA) * p_up[k] + ALPHA * delta;
  p_down[k] = (1-ALPHA) * p_down[k] - ALPHA * delta;
}

void succ(node_t *p, p_queue **h, int m, int n, double **a, double *b, double *c,
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
  p_queue *h = new_set(p);

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
      p_queue *temp = h->succ;
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
    p_queue *temp = h->succ;
    free(h);
    h = temp;
  }

  if (z == -INFINITY) {
    return NAN;
  } else {
    return z;
  }
}

int main() {
  int i, j;
  int m;
  int n;
  scanf("%d %d", &m, &n);
  double *c;
  double **a;
  double *b;

  c = (double*)malloc(n * sizeof(double));
  a = (double**)malloc(m * sizeof(double *));
  for (i = 0; i < m; i += 1) {
    a[i] = (double*)malloc((n + 1) * sizeof(double));
  }
  b = (double*)malloc(m * sizeof(double));

  for (i = 0; i < n; i += 1) {
    scanf("%lf", &c[i]);
  }

  for (i = 0; i < m; i += 1) {
    for (j = 0; j < n; j += 1) {
      scanf("%lf", &a[i][j]);
    }
  }

  for (i = 0; i < m; i += 1) {
    scanf("%lf", &b[i]);
  }

  double *x = malloc(m * sizeof(double));
  double res = intopt(m, n, a, b, c, x); // Call either intopt or simplex
  printf("ANSWER IS: %lf \n", res);
  free(x);
  for (i = 0; i < m; i += 1) {
    free(a[i]);
  }
  free(a);
  free(b);
  free(c);
  return 0;
}
