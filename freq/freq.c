#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUF_SIZE 256

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
		(*head)->word = (char*)malloc(strlen(word) + 1);
		strcpy((*head)->word, word);
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
	current->word = (char*)malloc(strlen(word) + 1);
	strcpy(current->word, word);
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
			free(current->word);
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
	char word[BUF_SIZE];
	while(head != NULL) {
		if (head->count > count) {
			count = head->count;
			strcpy(word, head->word);
		}
		head = head->next;
	}
	printf("result: %s %d\n", word, count);
}

void free_list(list* head) {
	if (head == NULL) return;
	list* current;
	list* next = head->next;
	free(head->word);
	free(head);
	while(next != NULL) {
		current = next;
		next = next->next;
		free(current->word);
		free(current);
	}
}

int main() {
	char buffer[BUF_SIZE];
	int line = 1;
	list* words = NULL;
	while(fscanf(stdin, "%s", buffer) != EOF) {
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
