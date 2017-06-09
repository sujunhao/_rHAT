#ifndef HAFFMAN_H_
#define HAFFMAN_H_

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define MX_HCODE_L 255



typedef struct {
	int l, r;
	double freq;
	int code_n;	
}hnode;

int cmpfunc1 (const void * a, const void * b);
int cmpfunc2 (const void * a, const void * b);

int *init_pQ(int n, int o);
int pop_pQ(int *q);
void push_pQ(int *q, int n, int o);
void create_htree(int r, char *s, int n, int o);
int creat_hnode(double *freq, int n, int o);
//enter symbol index return haffman code
char *int_get_hcode(int s, int o);
//enter haffman code return symbol index
int hcode_get_int(int r, char *s, int n, int o);
const char *byte_to_binary(uint64_t x);
//let a most 64 bit tar write in array with the last pos is l;
void set_HAbit(size_t *l, uint64_t *p, int tar, int o);
//let a most 64 bit tar write in array with the last pos is l;
//for have code 0000 should alway read
void get_HAbit(size_t *l, uint64_t *p, uint64_t *tar, int o);
int create_HAFF();


#endif