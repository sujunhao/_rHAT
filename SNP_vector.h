#ifndef SNP_VECTOR_H_
#define SNP_VECTOR_H_

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <stdbool.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include "func.h"



#define ENCODE_REF_BIT 2
#define ENCODE_REF_BLOCK 64

//ref malloc max size default
#define init_size 1e8
#define unit_snp_len 4
#define unit_ref_clen 2


KSEQ_INIT(gzFile, gzread)


extern int map_bit[], map_bit4[];
extern char map_c[], map_c4[];


//ref seq unit
typedef struct {
	//len, index, name, comment
	size_t l, i;
	char *n, *c;
} ref_type;

//whole ref seq
typedef struct {
	size_t ref_s_l;
	uint64_t *ref_s;

	size_t uref_l;
	ref_type *uref;
} ref_s;

//variant unit
typedef struct {
	char *rsid;
	char *chr;
	size_t pos;
	size_t n_allele;
	char **allele;
	bool if_snp;
} v_type;

//variant type have two bit vector for snp and mnp
typedef struct {
	//N(ref)
	uint64_t s_m;
	//snp type bit vector and sampling
	uint64_t *bit_v1;
	uint64_t *bv1_s;

	//mnp type bit vector and sampling
	uint64_t *bit_v2;
	uint64_t *bv2_s;

	//variant
	size_t v_l;

	//snp vector
	size_t snp_l;
	uint64_t *snp;

	size_t mnp_l;
	size_t last_p;
	uint64_t *mnp;
} var_s;

void init_ref(ref_s *rs, char *ref_fn);
void free_ref(ref_s *rs);
void read_ref(ref_s *rs, char *ref_fn);

//read ATCGX map to 0-3
int get_bit(char x);
//store string s in rs_b in 64 int format, start in b_l in rs_b; get ref char
void ref_set2int(char *s, size_t s_l, uint64_t *rs_b, size_t *b_l);
char get_ref_c(ref_s *rs, size_t p);
void print_ref_in_char(ref_s *rs, char *fn);

void init_var(ref_s *rs, var_s *vs, char *vcf_fn);
void free_var(var_s *vs);
void read_var(ref_s *rs, var_s *vs, char *vcf_fn);
void read_var_stat(char *vcf_fn);

//save_bit1(bit_v1, pos, is_snp, bv1_s, snp_i);
//save_bit2(bit_v2, pos, is_mnp, bv2_s, mnp_i);
//save_sv(snp_l, SV, variant_unit)
//save_mv(mnp_l, MV, variant_unit)
void save_bit1(var_s *vs, size_t pos, size_t *last_p);
void save_bit2(var_s *vs, size_t pos, size_t *last_p);

void save_sv(var_s *vs, char **al, size_t al_n);
void save_mv(var_s *vs, char **al, size_t al_n);

//mark bit vector bv pos p to 1  
void mark_mbit(uint64_t *bv, size_t p);
//get bit vector pos p value (0 or 1)
int get_mbit(uint64_t *bv, size_t p);

//for var snp vector set pos [l] to store snp info; get snp info in int
void mark_snp(size_t l, uint64_t *snp, char **al, size_t snp_l);
int get_snp(size_t l, uint64_t *snp);

// uint64_t int2bit(uint64_t n, uint64_t len);
void mark_mnp(size_t *l, uint64_t *mnp, char **al, size_t al_n);

//print to check if read vcf right
void log_vcf_v(var_s *vs, ref_s *rs, char *fn);

//query var result store in rns
void get_var(size_t *l, size_t *r, ref_s *rs, var_s *vs, char *rns, size_t *rns_l);

void set_Abit(size_t *l, uint64_t *p, uint64_t tar, size_t tl);
void get_Abit(size_t *l, uint64_t *p, uint64_t *tar, size_t tl);


unsigned long string_hash(unsigned char *str);

int get_pos_from_chrom(size_t *l, size_t *r, ref_s *rs, char *chrom);


void parse_pattern_string_4bit(char *p, size_t l, char *new_p);


#endif