#include "SNP_vector.h"
#include "Haffman.h"


const uint64_t last_mark = 0x3;
//map A:00, C:01, G:10, T:11
int map_bit[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3};
//map A:0x0001, C:0x0010, G:0x0100, T:0x1000 N:0x1111
int map_bit4[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 
	0, 4, 0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 0, 8, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 4, 
	0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 0, 8};
//map 0:A, 1:C, 2:G, 3:T
char map_c[]={'A', 'C', 'G', 'T'};
//map 0x0001:A, 0x0010:C, 0x0100:G, 0x1000:T, 0x1111:N
char map_c4[]={0, 'A', 'C', 0, 'G', 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 'N'};

size_t ref_buffer=0, snp_buffer=0, mnp_buffer=0, bitv_buffer=0;
size_t same_pos=0;
size_t alt_n_mx=0, alt_unit_mx=0;
const size_t max_n = 1000;
size_t alt_n_cnt[1000], alt_u_cnt[1000];
size_t alt_n_mx_allow = 19, alt_u_mx_allow = 50;

int alt_n_l = 6, alt_u_l = 6;


inline int get_bit(char x)
{
	return (x =='N') ? rand()%4 : map_bit[x];
}


void init_ref(ref_s *rs, char *ref_fn)
{
	size_t rl=1, nl=1;
	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen(ref_fn, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		printf("%2lu %13s len: %lu i:%lu\n", nl, seq->name.s, seq->seq.l, rl);
		// if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		// printf("seq: %s\n", seq->seq.s);
		// if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
		rl += seq->seq.l;
		++nl;
	}
	printf("total Len:%lu\n", rl);
	rs->ref_s_l=0;
	rs->ref_s = (uint64_t *)malloc(sizeof(uint64_t) * (rl/32+1));
	rs->uref_l=0;
	rs->uref = (ref_type *)malloc(sizeof(ref_type) * nl);
	
	ref_buffer += (sizeof(uint64_t) * (rl/32+1) + sizeof(ref_type) * nl);
	kseq_destroy(seq); // STEP 5: destroy seq
	gzclose(fp);
}

void free_ref(ref_s *rs)
{
	free(rs->ref_s);
	for (size_t i = 0; i<(rs->uref_l); ++i)
	{
		free((rs->uref)[i].n);
		free((rs->uref)[i].c);
	}
	free((rs->uref));
}

void read_ref(ref_s *rs, char *ref_fn)
{
	gzFile fp;
	kseq_t *seq;
	int l;
	uint64_t unit_number = 0;
	fp = gzopen(ref_fn, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		// printf("name: %s\n", seq->name.s);
		// if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		// printf("seq: %s\n", seq->seq.s);
		// if (seq->qual.l) printf("qual: %s\n", seq->qual.s);

		//set ference value
		ref_type urt;
		char *seq_n = (char *)malloc(sizeof(char) * (seq->name.l + 1));
		strcpy(seq_n, seq->name.s);
		char *seq_c = (char *)malloc(sizeof(char) * (seq->comment.l + 1));
		strcpy(seq_c, seq->comment.s);
		urt.n = seq_n;
		urt.c = seq_c;
		urt.l = seq->seq.l;
		urt.i = rs->ref_s_l;

		ref_set2int(seq->seq.s, seq->seq.l, rs->ref_s, &(rs->ref_s_l));

		rs->uref[rs->uref_l++] = urt;
	}
	// printf("return value: %d\n", l);
	kseq_destroy(seq); // STEP 5: destroy seq
	gzclose(fp);
}


void ref_set2int(char *s, size_t s_l, uint64_t *rs_b, size_t *b_l)
{
	size_t _l = (*b_l), _i = 0;

	_i = _l / 32;

	size_t tmp_l = (_l % (32)) * 2;
	size_t tmp = (rs_b[_i]) >> (64-tmp_l);

	// // printf("%lu %lu\n", tmp_l, ref_b[b_i]);
	for (size_t i=0; i<s_l; ++i)
	{
		tmp = (tmp<<2) | get_bit(s[i]);
		tmp_l += 2;

		if (tmp_l == 64)
		{
			rs_b[_i++] = tmp;
			tmp_l = 0;
		}
	}
	if (tmp_l != 0)
	{
		rs_b[_i++] = (tmp<<(64-tmp_l));
	}	

	// // printf("%lu\n", ref_b[b_i-1]);
	(*b_l) += s_l;
} 

void init_var(ref_s *rs, var_s *vs, char *vcf_fn)
{
	memset(alt_n_cnt, 0, sizeof(alt_n_cnt));
	memset(alt_u_cnt, 0, sizeof(alt_u_cnt));
	size_t vl=0, snpl=0;
	htsFile *fp    = hts_open(vcf_fn,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec    = bcf_init();

    while ( bcf_read(fp, hdr, rec)>=0 )
    {
    	if (bcf_is_snp(rec)) ++snpl;
    	++vl;
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",vcf_fn,ret);
        exit(ret);
    }

	vs->s_m = 64 * 5;
 	vs->bit_v1 = (uint64_t *)malloc(sizeof(uint64_t) * (rs->ref_s_l / 64 + 1));
 	vs->bit_v2 = (uint64_t *)malloc(sizeof(uint64_t) * (rs->ref_s_l / 64 + 1));
 	memset(vs->bit_v1, 0, sizeof(uint64_t) * (rs->ref_s_l / 64 + 1));
 	memset(vs->bit_v2, 0, sizeof(uint64_t) * (rs->ref_s_l / 64 + 1));
	vs->bv1_s = (uint64_t *)malloc(sizeof(uint64_t) * (rs->ref_s_l / vs->s_m + 1));
	vs->bv2_s = (uint64_t *)malloc(sizeof(uint64_t) * (rs->ref_s_l / vs->s_m + 1));
	
	vs->v_l = 0;
	vs->snp_l = 0;
	vs->snp = (uint64_t *)malloc(sizeof(uint64_t) * (snpl/16 + 1));
 	memset(vs->snp, 0, sizeof(uint64_t) * (snpl/16 + 1));
	
	vs->mnp_l = 0;
	vs->last_p = 0;
	// vs->mnp = (uint64_t *)malloc(sizeof(uint64_t) * ((vl-snpl)/2 + 1));
 // 	memset(vs->mnp, 0, sizeof(uint64_t) * ((vl-snpl)/2 + 1));
 	vs->mnp = (uint64_t *)malloc(sizeof(uint64_t) * ((vl-snpl) + 1));
 	memset(vs->mnp, 0, sizeof(uint64_t) * ((vl-snpl) + 1));

 	bitv_buffer += 2 * ((sizeof(uint64_t) * (rs->ref_s_l / 64 + 1)) + sizeof(uint64_t) * (rs->ref_s_l / vs->s_m + 1));
 	snp_buffer += (sizeof(uint64_t) * (snpl/16 + 1));
 	mnp_buffer += (sizeof(uint64_t) * ((vl-snpl)/2 + 1));
}

void free_var(var_s *vs)
{
	free(vs->bit_v1);
	free(vs->bit_v2);
	free(vs->bv1_s);
	free(vs->bv2_s);
	free(vs->snp);
	free(vs->mnp);
}

void read_var(ref_s *rs, var_s *vs, char *vcf_fn)
{
	htsFile *fp    = hts_open(vcf_fn,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec    = bcf_init();

    int ref_i = -1;
    size_t ref_i_begin=0;
    size_t same_c = 0, tmp_l=0;;
    char chrom_n[150]="", chrom_l[150]="";
    bool is_snp=false;
    size_t v_pos, last_p1=0, last_p2=0, last_p=0;

    while ( bcf_read(fp, hdr, rec)>=0 )
    {
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);

        //POS:rec->pos IS_SNP:bcf_is_snp(rec) CHR:bcf_hdr_id2name(hdr, rec->rid)
        //RSID:rec->d.id REF ALT:rec->n_allele, rec->allele[i]
		is_snp = bcf_is_snp(rec);
        strcpy(chrom_n, bcf_hdr_id2name(hdr, rec->rid));
        
        

        // check CHROM
        if (strcmp(chrom_l, chrom_n) != 0)
        {
        	// printf("chrom switch: %s to %s\n", chrom_n, vt.chr);
        	strcpy(chrom_l, chrom_n);
        	++ref_i;
        	if (ref_i >= rs->uref_l)
        	{
        		printf("ref number out of range\n");
        		break;
        	}
        	ref_i_begin = rs->uref[ref_i].i;
        }
       	else if (rec->pos < last_p)
       	{
       		//if pos appear again ignore
       		fprintf(stderr, "\rsame POS confilct: %s %lu %lu\n", chrom_n, last_p, (unsigned long)rec->pos);
       		++same_c;
       		continue;
       	}
        last_p = rec->pos + strlen(rec->d.allele[0]);
        v_pos = ref_i_begin + rec->pos;
        if (v_pos >= rs->ref_s_l) 
        {
        	printf("vpos out range,vpos:%lu res_s_l:%lu ref_i:%lu pos:%lu\n", v_pos, rs->ref_s_l, 1UL*ref_i_begin, 1UL*rec->pos);
        	continue;
        }

        
    	if (rec->d.allele[0][0] != 'N' && (get_ref_c(rs, v_pos) != rec->d.allele[0][0]))
    	{
    		printf("ref error in chrom:%s pos:%lu ref:%c REF:%s\n", chrom_n, (unsigned long)last_p, get_ref_c(rs, v_pos), rec->d.allele[0]);
    		continue;
    	}

        //save var
		// vs->v[vs->v_l] = vt;
		// printf("%lu %s %d %d %s\n", v_pos, vt.rsid, is_snp, get_mbit(ss->bit_v1, v_pos), rec->d.allele[0]);
		
        if (is_snp)
        {
        	//save bit vector 1
        	save_bit1(vs, v_pos, &last_p1);
        	//save snp_vector
        	save_sv(vs, rec->d.allele,rec->n_allele);
    		++(vs->snp_l);
        }
        else
        {
        	//save bit vector 2
        	save_bit2(vs, v_pos, &last_p2);
        	//save mnp_vector
        	save_mv(vs, rec->d.allele,rec->n_allele);
    		++(vs->mnp_l);
        }

        ++(vs->v_l);
    }


    same_pos = same_c;
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",vcf_fn,ret);
        exit(ret);
    }
}

void read_var_stat(char *vcf_fn)
{
	htsFile *fp    = hts_open(vcf_fn,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec    = bcf_init();

    size_t same_c = 0, tmp_l=0;
    char chrom_n[150]="", chrom_l[150]="";
    bool is_snp=false;
    size_t v_pos, last_p1=0, last_p2=0, last_p=0;

    while ( bcf_read(fp, hdr, rec)>=0 )
    {
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);

        //POS:rec->pos IS_SNP:bcf_is_snp(rec) CHR:bcf_hdr_id2name(hdr, rec->rid)
        //RSID:rec->d.id REF ALT:rec->n_allele, rec->allele[i]
		is_snp = bcf_is_snp(rec);
        strcpy(chrom_n, bcf_hdr_id2name(hdr, rec->rid));
        
        if (!is_snp)
        {
        	alt_n_mx = MX(alt_n_mx, rec->n_allele);
        	if (rec->n_allele > max_n) 
        	{
        		fprintf(stderr, "number of vcf ref len out of range: %d max enable: %lu\n", rec->n_allele, max_n);
        		rec->n_allele = max_n;
        	}
        	++alt_n_cnt[rec->n_allele];
	        for (int i=0; i<rec->n_allele; ++i)
	        {
	        	tmp_l = strlen(rec->d.allele[i]);
	        	alt_unit_mx = MX(alt_unit_mx, tmp_l);
	        	if (tmp_l > max_n)
	        	{
	        		fprintf(stderr, "number of vcf alt len out of range: %lu max enable: %lu\n", tmp_l, max_n);
	        		tmp_l = max_n;
	        	}
	        	++alt_u_cnt[tmp_l];
	        }
        }

    }


    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",vcf_fn,ret);
        exit(ret);
    }
}

void save_bit1(var_s *vs, size_t v_pos, size_t *last_p)
{
	mark_mbit(vs->bit_v1, v_pos);
    while ( ((*last_p) * vs->s_m) <= v_pos)
    {
    	vs->bv1_s[(*last_p)++] = vs->snp_l;
    }
}

void save_bit2(var_s *vs, size_t v_pos, size_t *last_p)
{
	mark_mbit(vs->bit_v2, v_pos);
    while ( ((*last_p) * vs->s_m) <= v_pos)
    {
    	vs->bv2_s[(*last_p)++] = vs->last_p;
    }
}

void save_sv(var_s *vs, char **al, size_t al_n)
{
    mark_snp(vs->snp_l, vs->snp, al, al_n);
}

void save_mv(var_s *vs, char **al, size_t al_n)
{
	mark_mnp(&(vs->last_p), vs->mnp, al, al_n);
}

void mark_mbit(uint64_t *bv, size_t p)
{
	uint64_t z = 1UL;
	size_t n=p/64, i=p%64;
	while (i<63) 
	{
		++i;
		z = z << 1;
	}
	bv[n] = (bv[n] | z);
}

int get_mbit(uint64_t *bv, size_t p)
{
	size_t n=p/64, i=p%64;
	uint64_t z = bv[n];
	while (i<63) 
	{
		++i;
		z = z >> 1;
	}
	return (z & 1);
}


int get_snp(size_t l, uint64_t *snp)
{
	size_t n=l/(64/unit_snp_len), i=l%(64/unit_snp_len);
	uint64_t t = snp[n];
	while (i!=0)
	{
		t = t << 4;
		--i;
	}
	return (int)(t>>60);
}

//let a most 64 bit tar write in array with the last pos is l;
inline void set_Abit(size_t *l, uint64_t *p, uint64_t tar, size_t tl)
{
	// printf("%lu %lu\n", tar, tl);
	int m = (*l)/64, n=(*l)%64;
	if (n + tl <= 64)
	{
		p[m] = p[m] | (tar << (64 - n - tl));
		// printf("%s\n", byte_to_binary(p[m]));
	}
	else 
	{
		p[m] = p[m] | (tar >> (tl - (64 - n)));
		p[m+1] = (tar << (128 - tl - n));
		// printf("%s\n", byte_to_binary(p[m]));
		// printf("%s\n", byte_to_binary(p[m+1]));
	}
	(*l) += tl;
}

//let a most 64 bit tar write in array with the last pos is l;
inline void get_Abit(size_t *l, uint64_t *p, uint64_t *tar, size_t tl)
{
	*tar = 0;
	int m = (*l)/64, n=(*l)%64;
	if (n + tl <= 64)	
		(*tar) = (p[m] << n >> n >> (64 - n - tl));
	else 
	{
		(*tar) = (p[m] << n >> n << (tl - (64 - n)));
		(*tar) = (*tar) | (p[m+1] >> (128 - tl - n));
	}
	(*l) += tl;
}

void mark_mnp(size_t *l, uint64_t *mnp, char **al, size_t al_n)
{
	size_t tmp=0;
	if (al_n > alt_n_mx_allow)
	{
		printf("have alt_n_mx_allow out of range\n");
		al_n = alt_n_mx_allow;
	}
	// set_Abit(l, mnp, al_n, alt_n_l);
	set_HAbit(l, mnp, al_n, 0);
	for (int i=0; i<al_n; ++i)
	{
		tmp = strlen(al[i]);
		if (tmp > alt_u_mx_allow)
		{
			tmp = alt_u_mx_allow;
		}
		// 	// set_Abit(l, mnp, 0, alt_u_l);
		// 	set_HAbit(l, mnp, 0, 1);
		// else
		{
			set_HAbit(l, mnp, tmp, 1);
			for (int j=0; j<tmp; ++j)
			{
				set_Abit(l, mnp, get_bit(al[i][j]), 2);
			}
		}
	}
	
}



void mark_snp(size_t l, uint64_t *snp, char **al, size_t snp_l)
{
	uint64_t t = 0;
	for (int i=1; i<snp_l; ++i)
	{
		t = t | map_bit4[al[i][0]];
	}

	size_t n=l/(64/unit_snp_len), i=(64/unit_snp_len) - 1 - l%(64/unit_snp_len);
	while (i!=0)
	{
		t = t << 4;
		--i;
	}
	snp[n] = snp[n]|t;
	return;
}


//64 bit to 32 len string
void int2str(uint64_t x, char *s, size_t l)
{
	int t;
	uint64_t tx=x;
	if (l==32)
	for (int i=l-1; i>=0; --i)
	{
		t = (tx>>(i*2)) & last_mark;
		s[31-i] = map_c[t];
	}
	else
	{
		tx = x >> (64 - l*2);
		for (int i=l-1; i>=0; --i)
		{
			t = (tx>>(i*2)) & last_mark;
			s[l-1-i] = map_c[t];
		}
	}
	s[l] = '\0';
	// printf("%s\n", s);
}

void print_ref_in_char(ref_s *rs, char *fn)
{
	size_t l = rs->ref_s_l / 32;
	size_t last = rs->ref_s_l % 32;

	FILE  *out;
	out = fopen(fn, "w+");

	char s[33];

	for (size_t i = 0; i<l; ++i)
	{
		int2str(rs->ref_s[i], s, 32);
		fprintf(out, "%s", s);
		// if (i%2) fprintf(out, "\n");
	}

	if (last > 0)
	{
		int2str(rs->ref_s[l], s, last);
		fprintf(out, "%s", s);
	}

	fprintf(out, "\n");

	for (size_t i = 0; i<rs->uref_l; ++i)
	{
		fprintf(out, "%lu ", rs->uref[i].l);
		fprintf(out, "%s\n", rs->uref[i].n);
		// if (i%2) fprintf(out, "\n");
	}
}

inline char get_ref_c(ref_s *rs, size_t p)
{
	size_t l = p / 32;
	size_t la = p % 32, t=31;

	uint64_t tmp = rs->ref_s[l];
	while (t > la)
	{
		tmp = (tmp >> 2);
		--t;
	}
	tmp = (tmp & last_mark);
	return map_c[tmp];
}


void log_vcf_v(var_s *vs, ref_s *rs, char *fn)
{
	FILE  *out;
	out = fopen(fn, "w+");

	size_t a = 0, b = rs->ref_s_l, tmp, tar_n, tar_u, ref_i=0, the_i = 0;
	size_t cnt=0, sum=0;
	a = a / vs->s_m * vs->s_m;
	b = MN(b, rs->ref_s_l);
	size_t s_cnt = vs->bv1_s[a / vs->s_m], m_cnt = vs->bv2_s[a / vs->s_m];
	fprintf(out, "%lu\n", b);

	printf("ref_i                 name          i        cnt        sum\n");
	printf("%5lu %20s %10lu ", ref_i, rs->uref[ref_i].n, 0UL);
	for (size_t i=a; i<b; ++i)
	{
		while ((rs->uref[ref_i]).l + (rs->uref[ref_i]).i <= i)
		{
			++ref_i;
			the_i = (rs->uref[ref_i]).i;
			sum += cnt;
			printf("%10lu %10lu\n%5lu %20s %10lu ", cnt, sum, ref_i, rs->uref[ref_i].n, i);
			fprintf(out, "new chr %s\n", rs->uref[ref_i].n);
			cnt = 0;
		} 
		if (get_mbit(vs->bit_v1, i))
		{
			++cnt;
			fprintf(out, "%lu\t%lu\t", ref_i+1, i-the_i+1);
			fprintf(out, "%c\t", get_ref_c(rs, i));
			bool is_ref=false;
			int o = get_snp(s_cnt, vs->snp);
			if (o & 1) 
			{
				fprintf(out, is_ref?",A":"A");
				is_ref = true;
			}
			if (o & 2) 
			{
				fprintf(out, is_ref?",C":"C");
				is_ref = true;
			}
			if (o & 4)
			{
				fprintf(out, is_ref?",G":"G");
				is_ref = true;
			}
			if (o & 8)
			{
				fprintf(out, is_ref?",T":"T");
				is_ref = true;
			}
			fprintf(out, "\n");

			++s_cnt;
		}

		if (get_mbit(vs->bit_v2, i))
		{
			++cnt;
			fprintf(out, "%lu\t%lu\t", ref_i+1, i-the_i+1);
			// get_Abit(&m_cnt, vs->mnp, &tar_n, alt_n_l);
			get_HAbit(&m_cnt, vs->mnp, &tar_n, 0);
			// fprintf(out, "%lu ", tar_n);
			for (int j=0; j<tar_n; ++j)
			{
				// get_Abit(&m_cnt, vs->mnp, &tar_u, alt_u_l);
				get_HAbit(&m_cnt, vs->mnp, &tar_u, 1);
				// fprintf(out, "%lu ", tar_u);

				for (int k=0; k<tar_u; ++k)
				{
					get_Abit(&m_cnt, vs->mnp, &tmp, 2);
					// fprintf(out, "%lu ", tmp);
					fprintf(out, "%c", map_c[tmp]);
				}
				fprintf(out, j==0||j==tar_n-1?"\t":",");
			}
			fprintf(out, "\n");
		}
		if (i==b-1)
		{
			sum += cnt;
			printf("%10lu %10lu\n", cnt, sum);
		}
	}
	// fclose(out);
}

void get_var(size_t *ll, size_t *rr, ref_s *rs, var_s *vs, char *rns, size_t *rns_l)
{
	size_t l = *ll, r = *rr;
	if (r<=l) return;
	l = ( l < rs->ref_s_l ) ? l : rs->ref_s_l;
	l = l/(vs->s_m)*(vs->s_m);
	r = ( r < rs->ref_s_l ) ? r : rs->ref_s_l;
	// *ll = l;
	*rr = r;
	size_t tt=(*rns_l), tar_n, tar_u, tmp;
	size_t s_cnt = vs->bv1_s[l / vs->s_m], m_cnt = vs->bv2_s[l / vs->s_m];
	for (size_t k = l; k < *ll; ++k)
	{
		if (get_mbit(vs->bit_v1, k))
		{
			++s_cnt;
		}
		else if (get_mbit(vs->bit_v2, k))
		{
			get_HAbit(&m_cnt, vs->mnp, &tar_n, 0);
			for (int i=0; i<tar_n; ++i)
			{
				get_HAbit(&m_cnt, vs->mnp, &tar_u, 1);

				
				for (int j=0; j<tar_u; ++j)
				{
					get_Abit(&m_cnt, vs->mnp, &tmp, 2);
				}

			}
		}

	}

	for (size_t k = *ll; k < r; ++k)
	{
		if (get_mbit(vs->bit_v1, k))
		{
			rns[tt++] = get_ref_c(rs, k);
			rns[tt++]='[';
			
			int o = get_snp(s_cnt, vs->snp);
			if (o & 1) rns[tt++]='A';
			if (o & 2) rns[tt++]='C';
			if (o & 4) rns[tt++]='G';
			if (o & 8) rns[tt++]='T';
			// rns[tt++] = o;
			++s_cnt;
			rns[tt++]=']';
		}
		else if (get_mbit(vs->bit_v2, k))
		{
			// get_Abit(&m_cnt, vs->mnp, &tar_n, alt_n_l);
			get_HAbit(&m_cnt, vs->mnp, &tar_n, 0);
			for (int i=0; i<tar_n; ++i)
			{
				// get_Abit(&m_cnt, vs->mnp, &tar_u, alt_u_l);
				get_HAbit(&m_cnt, vs->mnp, &tar_u, 1);

				if (i==0)
				{
					rns[tt++]='(';
					//should sure that len(ref) > 0
					k += MX((tar_u-1), 0);
				}
				else if (i==1) rns[tt++]='[';
				for (int j=0; j<tar_u; ++j)
				{
					get_Abit(&m_cnt, vs->mnp, &tmp, 2);
					// fprintf(out, "%lu ", tmp);
					rns[tt++]=map_c[tmp];
				}
				if (i==0) rns[tt++]=')';
				else if (i==tar_n-1) rns[tt++]=']';
				else rns[tt++]=',';

			}
		}
		else rns[tt++] = get_ref_c(rs, k);

	}
	*rns_l = tt;
}

unsigned long string_hash(unsigned char *str)
{
    unsigned long hash = 5381;
    int c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}


int get_pos_from_chrom(size_t *l, size_t *r, ref_s *rs, char* chrom)
{
	for (size_t i=0; i<(rs->uref_l); ++i)
	{
		if (strcmp(chrom, (rs->uref)[i].n) == 0)
		{
			// if ((*l < 0 || *l >=(rs->uref)[i].l) || \
			// 	(*r < 0 || *r >=(rs->uref)[i].l))
			// {
			// 	printf("get chrom error, position out of range\n");
			// 	*l = MN(*l, (rs->uref)[i].l);
			// 	*r = MN(*r, (rs->uref)[i].l);
			// }
			*l = MX(0, *l);
			*r = MX(0, *r);
			*l = MN(*l, (rs->uref)[i].l) + (rs->uref)[i].i;
			*r = MN(*r, (rs->uref)[i].l) + (rs->uref)[i].i;
			// printf("query chrom:%s, match:%s, index:%lu, len:%lu\n", chrom,\
			// 				(rs->uref)[i].n, (rs->uref)[i].i, (rs->uref)[i].l);
			// printf("begin:%lu, end:%lu len:%lu\n", *l, *r, *r-*l);
			return 1;
		}
	}
	return 0;
}

// parse pattern string from char  to bit4. ie A to 0001
void parse_pattern_string_4bit(char *p, size_t l, char *new_p)
{
	for (size_t i=0; i<l; ++i)
	{
		// if (p[i] == 'N')	new_p[i] = 0;
		// else 	new_p[i] = map_bit4[p[i]];
		new_p[i] = map_bit4[p[i]];
	}
	for (size_t i=l; i<(l+9); ++i)
	{
		new_p[i] = 0;
	}
	return;
}