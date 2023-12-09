#include <ctype.h>
#include <stdio.h>
#define SIZE 10
int stack[SIZE];
int main() {
	int ch, val, sp = -1, line = 1;
	char err = '-';
	while ((ch = getchar()) != EOF) {
next:
		if(isdigit(ch)) {
			int v = ch - '0';
			while((isdigit(ch = getchar()))) {
				v = v * 10 + ch - '0';
			}
			if (sp == SIZE - 1) {
				err = (char)v + '0';
				goto eol;
			}
			else stack[++sp] = v;
			goto next;
		} else if(ch == ' ') continue;
		else if(ch == '\n') {
eol:
			if (err != '-') {
				printf("line %d: error at %c\n", line++, err);	
				while(ch != '\n') {
					ch = getchar();
				};
				err ='-', sp = -1; 
			} else if (sp != 0) {
				printf("line %d: error at \\n\n", line++);
			} else {
				val = stack[sp--];
				printf("line %d: %d\n", line++, val);
			}
		}
		else {
			switch(ch) {
				case '+':
					if (sp < 1) {
						err = '+';
						goto eol;
					}
					val = stack[sp] + stack[sp-1];
					stack[--sp] = val;
					break;
				case '-':
					if (sp < 1) {
						err = '-';
						goto eol;
					}
					val = stack[sp-1] - stack[sp];
					stack[--sp] = val;
					break;
				case '*':
					if (sp < 1) {
						err = '*';
						goto eol;
					}
					val = stack[sp] * stack[sp-1];
					stack[--sp] = val;
					break;
				case '/':
					if (sp < 1 || stack[sp] == 0) {
						err = '/';
						goto eol;
					}
					val = stack[sp-1] / stack[sp];
					stack[--sp] = val;
					break;
				default:
					err = ch;
					goto eol;
			}
		}
	}
}
