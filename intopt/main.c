#include <stdio.h>
#include <stdlib.h>
#include "intopt.h"
#include "simplex.h"

int main() {
  int i, j;
  int m;
  int n;
  scanf("%d %d", &m, &n);
  double *c;
  double **a;
  double *b;

  c = malloc(n * sizeof(double));
  a = malloc(m * sizeof(double *));
  for (i = 0; i < m; i += 1) {
    a[i] = malloc((n + 1) * sizeof(double));
  }
  b = malloc(m * sizeof(double));

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
