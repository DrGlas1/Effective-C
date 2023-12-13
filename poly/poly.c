#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct poly_t {
	unsigned char size;
	int* exp;
	int val[];
} poly_t;


void swap(int* a, int* b) {  
	int t = *a;  
	*a = *b;  
	*b = t;  
}  

int partition(int arr[], int comp[], int low, int high) {  
    int pivot = arr[high];  
    int i = (low - 1);  
  
	for (int j = low; j <= high - 1; j++) {  
		if (arr[j] > pivot) {  
			i++;  
			swap(&arr[i], &arr[j]);
			swap(&comp[i], &comp[j]);
		}  
	}  
	swap(&arr[i + 1], &arr[high]);  
	return (i + 1);  
}  

void quicksort(int arr[], int comp[], int low, int high) {  
	if (low < high) {  
		int pi = partition(arr, comp, low, high);  
		quicksort(arr, comp, low, pi - 1);  
		quicksort(arr, comp, pi + 1, high);  
	}  
}

char is_digit(char c) {
	return '0' <= c && c <= '9';
}

poly_t* new_poly_from_string(const char* var) {
	int len = strlen(var), exp = 0, val;
	char negative = 0;
	unsigned char i;
	unsigned char size = 1;
	for(i = 0; i < len; i++) {
		if (var[i] == 'x') size++;
	}
	poly_t* poly = (poly_t*)calloc(1,sizeof(poly_t) + size * sizeof(int));
	poly->exp = (int*)calloc(1,sizeof(int));
	poly->size = size;
	poly->val[0] = 1;
	size = 0;
	for(i = 0; i < len; i++) {
		if (var[i] == ' ' || var[i] == '+') continue;
		if (var[i] == '-') negative = 1;
		if (is_digit(var[i])) {
			val = var[i++] - '0';
			while(i < len && is_digit(var[i])) {
				val *= 10;
				val += var[i] - '0';
			}
			if(negative) {
				poly->val[size] = -val;
				negative = 0;
			} else	poly->val[size] = val;
		}
		if (var[i] == 'x') {
			i++;
			if (var[i] == '^') i++;
			while(i < len && is_digit(var[i])) {
				exp *= 10;
				exp += var[i++] - '0';
			}
			if (exp == 0) exp++;
			poly->exp[size] = exp;
			exp = 0;
			size++;
		}
	}
	return poly;
}

void free_poly(poly_t* poly) {
	free(poly->exp);
	free(poly);
}

poly_t* mul(poly_t* poly1, poly_t* poly2) {
	unsigned char size = poly1->size * poly2->size;
	poly_t* res = (poly_t*)calloc(1,sizeof(poly_t) + size * sizeof(int));
	res->exp = (int*)calloc(1,sizeof(int));
	unsigned char k = 0;
	int prev_exp = -1, exp;
	for(unsigned char i = 0; i < poly1->size; i++) {
		for(unsigned char j = 0; j < poly2->size; j++) {
			exp = poly1->exp[i] + poly2->exp[j]; 
			if (prev_exp == exp) {
				k--;
				res->val[k] += poly1->val[i] * poly2->val[j];
			} else {
				res->exp[k] = (prev_exp = exp);
				res->val[k] = poly1->val[i] * poly2->val[j];	
			}		
			k++;
		}
	}
	quicksort(res->exp, res->val, 0, k-1);
	res->size = k;
	return res;
}


void print_poly(poly_t* poly) {
    for (int i = 0; i < poly->size; i++) {
        if (poly->val[i] == 0)continue;
        if (poly->val[i] < 0 || i > 0) printf("%s", (poly->val[i] < 0) ? "- " : "+ ");
	if (abs(poly->val[i]) != 1 || poly->exp[i] == 0) printf("%d", abs(poly->val[i]));
        if (poly->exp[i] > 0) {
            printf("x");
            if (poly->exp[i] > 1) {
                printf("^%d", poly->exp[i]);
            }
        }
        if (i != poly->size - 1) printf(" ");
    }
    printf("\n");
}

