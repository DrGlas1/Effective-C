#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "simplex.h"
#define EPSILON 1e-6



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
  int i;
  for (i = 0; i < s->n; i += 1) {
    if (s->c[i] > EPSILON) {
      return i;
    }
  }
  return -1;
}

void pivot(simplex_t *s, int row, int col) {
  double **a = s->a;
  double *b = s->b;
  double *c = s->c;
  int m = s->m;
  int n = s->n;
  int i, j, t;
  t = s->var[col];
  s->var[col] = s->var[n + row];
  s->var[n + row] = t;
  s->y = s->y + c[col] * b[row] / a[row][col];

  for (i = 0; i < n; i += 1) {
    if (i != col) {
      c[i] = c[i] - c[col] * a[row][i] / a[row][col];
    }
  }
  c[col] = -c[col] / a[row][col];

  for (i = 0; i < m; i += 1) {
    if (i != row) {
      b[i] = b[i] - a[i][col] * b[row] / a[row][col];
    }
  }
  for (i = 0; i < m; i += 1) {
    if (i != row) {
      for (j = 0; j < n; j += 1) {
        if (j != col) {
          a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
        }
      }
    }
  }
  for (i = 0; i < m; i += 1) {
    if (i != row) {
      a[i][col] = -a[i][col] / a[row][col];
    }
  }
  for (i = 0; i < n; i += 1) {
    if (i != col) {
      a[row][i] = a[row][i] / a[row][col];
    }
  }
  b[row] = b[row] / a[row][col];
  a[row][col] = 1 / a[row][col];
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
          (row < 0 || b[i] / a[i][col] < b[row] / a[row][col])) {
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
  double *t;
  t = calloc(n, sizeof(double));

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
