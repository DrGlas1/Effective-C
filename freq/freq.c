#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUF_SIZE 1024

typedef struct list list;
struct list{
	list* next;
	char* word;
	int count;
};

int isprime(int num) {
	int i, temp = 0;
	for(i = 2; i <= num / 2; i++) {
		if (num % i == 0) {
			temp++;
			break;
		}
	}
	return temp == 0 && num != 1;
}

int add(list** head, const char* word) {
	if (*head == NULL) {
		*head = (list*)malloc(sizeof(list));
	        (*head)->next = NULL;
		(*head)->word = strdup(word);
		(*head)->count = 1;
		return 1;
	}

	list* current = *head;
	while (current != NULL) {
		if (strcmp(current->word, word) == 0) {
			current->count++;
			return 0;
		}
		if (current->next != NULL) {
			current = current->next;
		} else {
			break;
		}
	}

	current->next = (list*)malloc(sizeof(list));
	current = current->next;
	current->next = NULL;
	current->word = strdup(word);
	current->count = 1;
	return 1;
}

int delete(list** head, const char* word) {
	if(*head == NULL)
		return 0;
	list* prev = NULL;
	list* current = *head;
	while(current != NULL) {
		if (strcmp(current->word, word)== 0) {
			if (prev != NULL) {
				prev->next = current->next;
			} else {
				*head = current->next;
			}
			free(current);
			return 1;
		}
		prev = current;
		current = current->next;
	}
	return 0;
}

void print_most_common(list* head) {
	int count = 0;
	char* word;
	while(head != NULL) {
		if (head->count > count) {
			count = head->count;
			word = head->word;
		}
		head = head->next;
	}
	printf("result: %s %d\n", word, count);
}

void free_list(list* head) {
	if (head == NULL) return;
	list* curr;
	list* next = head->next;
	free(head);
	while(next != NULL) {
		curr = next;
		next = next->next;
		free(curr);
	}
}

int main() {
	int line = 1;
	char buffer[BUF_SIZE];
	list* words = NULL;
	while(fgets(buffer, BUF_SIZE, stdin) != NULL) {
		buffer[strcspn(buffer, "\n")] = '\0';
		if(isprime(line)) {
			printf("trying to delete %s: ", buffer);
			if (delete(&words, buffer)) {
				printf("deleted\n");
			} else {
				printf("not found\n");
			}
		}else {
			if (add(&words, buffer)) {
				printf("added %s\n", buffer);
			} else {
				printf("counted %s\n", buffer);
			}
		}
		line++;
	}
	print_most_common(words);
	free_list(words);
	return 0;
}
