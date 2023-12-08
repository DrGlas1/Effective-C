#include <stdio.h>
#define SIZE 10

int stack[SIZE];

char is_digit(char ch) {
	return '0' <= ch && ch <= '9';
}

int main() {
	char ch;
	char err = '-';
	int val;
	int line = 1;
	signed char sp = -1;
	while ((ch = getchar()) != EOF) {
next:
		if(is_digit(ch)) {
			int v = ch - '0';
			while((is_digit(ch = getchar()))) {
				v = v * 10 + ch - '0';
			}
			if (sp == 10);
			else stack[++sp] = v;
			goto next;
		} else if(ch == ' ') continue;
		else if(ch == '\n') {
			if (sp == -1) {
				printf("line %d: error at \\n\n", line);
			}
			else if (err != '-') {
				printf("line %d: error at %c\n", line, err);	
				err ='-'; 
			} else {
				printf("line %d: %d\n", line, stack[sp]);
			}
			line++;
			continue;
		}
		else {
			switch(ch) {
				case '+':
					if (sp < 1) {
						err = '+';
						continue;
					}
					val = stack[sp] + stack[sp-1];
					stack[--sp] = val;
					break;
				case '-':
					if (sp < 1) {
						err = '-';
						continue;
					}
					val = stack[sp-1] - stack[sp];
					stack[--sp] = val;
					break;
				case '*':
					if (sp < 1) {
						err = '*';
						continue;
					}
					val = stack[sp] * stack[sp-1];
					stack[--sp] = val;
					break;
				case '/':
					if (sp < 1 || stack[sp] == 0) {
						err = '/';
						continue;
					}
					val = stack[sp-1] / stack[sp];
					stack[--sp] = val;
					break;
				default:
					fprintf(stderr, "Error: Invalid char\n");
			}
		}
	}
}
