#include "Haffman.h"

//1-hcode_len is the code
int hcode_len[2] = {0, 0};
int mark_max[2] = {0, 0};
char hcode[2][MX_HCODE_L][500];
char hsymbol[2][MX_HCODE_L][20];
int theroot[2];

double freq[2][MX_HCODE_L];
char tmpcode[500];


hnode *thehnode[2];
//the haffman tree store in hnode,
//the theroot get get by creat_hnode()
//map symbol to haffman code by int_get_hcode()
//map haffman code to symbol by hcode_get_int()
//symbol store is hsmybol





int cmpfunc1 (const void * a, const void * b)
{
   return ( thehnode[0][*(int*)a].freq > thehnode[0][*(int*)b].freq );
}
int cmpfunc2 (const void * a, const void * b)
{
   return ( thehnode[1][*(int*)a].freq > thehnode[1][*(int*)b].freq );
}

int *init_pQ(int n, int o)
{
    static int *q;
    q = (int *)malloc(sizeof(int) * (2*n-1));
    for (int i=0; i<2*n-1; ++i) 
    	if (i<n) q[i] = i;
    	else q[i] = -1;

    if (o==0)	qsort(q, n, sizeof(int), cmpfunc1);
    else qsort(q, n, sizeof(int), cmpfunc2);

    // printf("ppp\n");
    return q;
}


int pop_pQ(int *q)
{
	// printf("pop\n");
	int tmp = q[0], i=1;
	while(q[i] != -1)
	{
		q[i-1] = q[i];
		++i;
	}
	q[i-1] = q[i];
	// printf("pop\n");
	return tmp;
}

void push_pQ(int *q, int n, int o)
{
	// printf("push\n");

	int i=0, tmp, l = n;
	while (q[i] != -1 && thehnode[o][q[i]].freq<thehnode[o][n].freq) ++i;
	// printf("push\n");
	while (q[i] != -1)
	{
		tmp = q[i];
		q[i] = l;
		l = tmp;
		++i;
	}
	q[i] = l;
	// printf("push\n");

}

void create_htree(int r, char *s, int n, int o)
{
	if (thehnode[o][r].l == -1)
	{
		s[n] = '\0';
		strcpy(hcode[o][thehnode[o][r].code_n], s);
		// printf("%d %s\n", thehnode[r].code_n, s);
	}
	else
	{
		s[n] = '1';
		create_htree(thehnode[o][r].l, s, n+1, o);
		s[n] = '0';
		create_htree(thehnode[o][r].r, s, n+1, o);
	}
	return;
}

int creat_hnode(double *freq, int n, int o)
{
	thehnode[o] = (hnode *)malloc(sizeof(hnode) * (2*n-1));
	for (int i=0; i<n; ++i)
	{
		// thehnode[o][i].code_n = i+1;
		thehnode[o][i].code_n = atoi(hsymbol[o][i]);
		thehnode[o][i].freq = freq[i];
		thehnode[o][i].l = -1;
		thehnode[o][i].r = -1;
	}
	int thei = n;
	int *pQ = init_pQ(n, o);
	// for (int j=0; pQ[j]!=-1; ++j) printf("%1.8lf ", thehnode[pQ[j]].freq);
	// 	printf("\n");
	int l, r, m;
	for (int i=0; i<n-1; ++i)
	{
		// printf("%d %d\n", i, n-1);
		// for (int j=0; pQ[j]!=-1; ++j) printf("%.8lf ", thehnode[pQ[j]].freq);
		// printf("\n");
		l = pop_pQ(pQ);
		r = pop_pQ(pQ);
		thehnode[o][thei].l = l;
		thehnode[o][thei].r = r;
		thehnode[o][thei].freq = thehnode[o][l].freq + thehnode[o][r].freq;
		push_pQ(pQ, thei, o);
		// for (int j=0; pQ[j]!=-1; ++j) printf("%.8lf ", thehnode[pQ[j]].freq);
		// printf("\n");
		++thei;
	}

	char tmpc[500];
	create_htree(pQ[0], tmpc, 0, o);
	theroot[o] = pQ[0];
	free(pQ);
	return theroot[o];
}


//enter symbol index return haffman code
char *int_get_hcode(int s, int o)
{
	return hcode[o][s];
}


//enter haffman code return symbol index
int hcode_get_int(int r, char *s, int n, int o)
{
	if (s[n]=='\0') return thehnode[o][r].code_n;
	if (thehnode[o][r].l != -1)
	{
		if (s[n]=='1') return hcode_get_int(thehnode[o][r].l, s, n+1, o);
		else if (s[n]=='0') return hcode_get_int(thehnode[o][r].r, s, n+1, o);	
	}
	else return -1;
}


const char *byte_to_binary(uint64_t x)
{
    static char b[65];
    b[0] = '\0';

    uint64_t z;
    for (z = (1UL << 63); z > 0; z >>= 1)
    {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }

    return b;
}


//let a most 64 bit tar write in array with the last pos is l;
void set_HAbit(size_t *l, uint64_t *p, int tar, int o)
{
	// if (tar >= mark_max[o])
	// {
	// 	printf("%d\n", tar);
	// 	tar = mark_max[o];
	// }
	if (tar > mark_max[o]) return;
	// printf("%lu %lu\n", tar, tl);
	char *th = int_get_hcode(tar, o);
	int tl = strlen(th);

	int m = (*l)/64, n = (*l)%64;
	uint64_t xmark = 1UL << (64 - n - 1);
	if (n + tl <= 64)
	{
		for (int i=0; i<tl; ++i)
		{
			if (th[i]=='1') p[m] = p[m] | xmark;
			xmark = xmark >> 1;
		}
		// printf("%s\n", byte_to_binary(p[m]));
	}
	else 
	{
		for (int i=0; i<tl; ++i)
		{
			if (n + i == 64)
			{
				xmark = (1UL)<<63;
				++m;
			}
			if (th[i]=='1') p[m] = p[m] | xmark;
			xmark = xmark >> 1;
		}
		// printf("%s\n", byte_to_binary(p[m]));
		// printf("%s\n", byte_to_binary(p[m+1]));
	}
	(*l) += tl;
}

//let a most 64 bit tar write in array with the last pos is l;
//for have code 0000 should alway read
void get_HAbit(size_t *l, uint64_t *p, uint64_t *tar, int o)
{
	*tar = 0;
	int m = (*l)/64, n=(*l)%64;

	int have_s = 0;
	int rn = theroot[o], tl=0;
	uint64_t xmark = 1UL << (64 - n - 1);
	while (!have_s)
	{
		if (thehnode[o][rn].l == -1)
		{
			if (thehnode[o][rn].freq < -1e16)	return;
			have_s = 1;
			*tar = thehnode[o][rn].code_n;
		}
		else
		{
			if (n + tl == 64)
			{
				xmark = 1UL << (64 - 1);
				++m;
			}
			if ((p[m] & xmark) != 0) rn = thehnode[o][rn].l;
			else rn = thehnode[o][rn].r;
			xmark = xmark >> 1;
			++tl;
		}
	}
	(*l) += tl;
}



int create_HAFF()
{
	FILE *fref;

	fref = fopen("freq.stat", "r");
	if (!fref)
	{
		fprintf(stderr, "can't find Haffman stat file: freq.stat\n");
		return -1;
	}

	int sn = 0;
	double sum_f[2];
	int root=0;

	for (int o=0; o<2; ++o)
	{
		sn = 0;
		fscanf(fref, "%d\n", &sn);
		hcode_len[o] = sn;
		sum_f[o] = 0;
		for (int i=0; i<sn; ++i) 
		{
			fscanf(fref, "%s %lf\n", hsymbol[o][i], &freq[o][i]);
			// printf("%s %lf\n", hsymbol[o][i], freq[o][i]);
			sum_f[o] += freq[o][i];
			if (i==sn-1)
			{
				mark_max[o] = atoi(hsymbol[o][i]);
				// printf("%d\n", mark_max[o]);
			}
		}
		for (int i=0; i<sn; ++i)
		{
			freq[o][i] = freq[o][i] / sum_f[o];
			// printf("%2.10lf\n", freq[o][i]);
		}

		root =  creat_hnode(freq[o], sn, o);

		// for (int i=0; i<sn; ++i)
		// 	printf("%d %s\n", atoi(hsymbol[o][i]), int_get_hcode(atoi(hsymbol[o][i]), o));
		// char qtmp[]="11111111111";
		// printf("%s %d\n", qtmp, hcode_get_int(root, qtmp, 0, o));
		
		// double cnt_h_len = 0;
		// for (int i=0; i<sn; ++i)
		// {
		// 	cnt_h_len += (freq[o][i]*sum_f[o]*strlen(int_get_hcode(i, o)));
		// }
		// printf("6 code len: %.0lf haf len:%.0lf save:%2.2lf%%\n", sum_f[o] * 6 , cnt_h_len, (sum_f[o] * 6 - cnt_h_len)/(sum_f[o] * 6) * 100);
		

		// uint64_t p[5]={0, 0, 0, 0, 0};
		// size_t l=0, t;
		// set_HAbit(&l, p, 10, o);
		// set_HAbit(&l, p, 1, o);
		// set_HAbit(&l, p, 10, o);
		// set_HAbit(&l, p, 10, o);
		// set_HAbit(&l, p, 100, o);
		

		// printf("%s\n", byte_to_binary(p[0]));
		// printf("%s\n", byte_to_binary(p[1]));
		// printf("%s\n", byte_to_binary(p[2]));

		// l = 0;
		// get_HAbit(&l, p, &t, o);
		// printf("%lu %lu\n", l, t);
		// get_HAbit(&l, p, &t, o);
		// printf("%lu %lu\n", l, t);
		// get_HAbit(&l, p, &t, o);
		// printf("%lu %lu\n", l, t);
		// get_HAbit(&l, p, &t, o);
		// printf("%lu %lu\n", l, t);
		// printf("\n");
	}

	printf("created Haffman code\n");
	fclose(fref);
	return 1;
}