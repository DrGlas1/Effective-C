#include <stdio.h>
#include <stdlib.h>

int main() {
{
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
		a[i] = (double*)malloc(sizeof(double) * n);
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
	
	for(i = 0; i < m; i++) {
		free(a[i]);
	}
	free(a);

	return 0;
}
}
