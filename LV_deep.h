#ifndef LV_DEEP_H
#define LV_DEEP_H


#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SNP_vector.h"
#include "Haffman.h"
#include "limits.h"
#define CountTrailingZeroes(x, ans) do{ans = __builtin_ctzll(x);} while(0)



//text string tree node unit
typedef struct {
	int id;
	int deep;
	size_t begin, len;
} deep_node;


//whole deep info
//all include all node in deep string
typedef struct {
	size_t s_l;
	char *deep_s;

	int n_l;
	int the_deep;
	deep_node *all_node;
} DI;


// int ***LV;
// LV = (int ***)malloc(sizeof(int **) * (di->n_l));
// LV save the LV info. for each block, error and diagonal
// db table
// init all as mark_no_vis

// init Far_reach mark if a diagoal have reach the end
// Far_reach[b][d] b: in which block, d: in which diagoal
// 0 if not reach end, 1 if have reach
// int **Far_reach;

//to save the boundry of deep
//to save succeed the lv info. between different block(deep)
// int **mark_bit;

typedef struct {
	int ***LV;
	int **Far_reach;
	int **mark_bit;

	int n_l;
	int max_error;
	int lp;
	int deep;
} LV_ENTITY;


// typedef enum {
// 	mat = 'M',
//     del = 'D',
//     ins = 'I',
// 	mis = 'X',
// 	// mat = '0',
//     // del = '1',
//     // ins = '2',
// 	// mis = '3',
// } Cigar;

//init the_deep = -1;
void init_node_space(DI* di, size_t s_l, int n_l);
void free_node_space(DI* di);

//unpdate node and di the_deep info.
void mark_node_id(DI *di, int node_index, int deep, size_t begin, size_t len);

void parse_target_string(DI *di, size_t *ll, size_t *rr, ref_s *rs, var_s *vs);

void log_DI(DI *di); 

void parse_pattern_string(char *p, size_t l, char *new_p);

void print_bit4(char o);
void print_bit4_(char o);

void init_lv_space(LV_ENTITY *lv, int n_l, int max_error, int lp, int deep);
void free_lv_space(LV_ENTITY *lv);

int init_lv(LV_ENTITY *lv, int n_l, int max_error, int lp, int deep, int mark_no_vis);

void run_deep_lv(LV_ENTITY *lv, DI *di, char *p, size_t lp, int max_error);

void get_Cigar(LV_ENTITY *lv, DI *di, char *pattern, int lp, int max_error, \
		int global_reach, int the_b, int the_e, int the_d);
#endif