#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "alignment.h"
#include "RHT.h"
#include "ksw.h"
#define PX(X) std::cout << X << std::endl
using namespace std;

#define PRINTALN

uint32_t t_wait = 1024;
const uint8_t seq_nt4_tablet[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int8_t  mat[25];
int8_t match=2, mism=-5, gape=2, gapo=1;
const char correspondTable[] = "MIDNSHP=X";
uint32_t *cigar;
int n_cigar = 0;
uint8_t  readqry[2048];
uint8_t  refqry[2048];
uint8_t *readqry_;
uint8_t *refqry_;


void cout_string_in_char(char* read, size_t l, size_t len, ofstream &outt)
{
    for (size_t i = 0; i<len; ++i) outt << read[i+l];
}
inline void transIntoDec(uint8_t *transtr,char *str, int length)
{
    for (int i=0;i<length;++i) {
        transtr[i] = seq_nt4_tablet[str[i]];
    }
}


typedef struct window_node {
    uint32_t index_of_W, index_of_R;
    size_t len;
}WINDOW_NODE;

class DAG
{
// private:
public:
	WINDOW_NODE *wd;
	size_t **V;
	uint32_t *dp;
	size_t node_i;
	size_t *path;
	size_t p_index;

public:
	DAG(size_t max_num)
	{
		wd = new WINDOW_NODE[max_num];
		node_i = 0;
	}
	~DAG()
	{
		if (V != NULL)
		{
			for (size_t i=0; i<node_i; ++i) delete [] V[i];
	        delete [] V;
		}
		if (dp != NULL)	delete [] dp;
		if (path != NULL) delete [] path;
        delete [] wd;
	}

    void clear()
    {
        node_i = 0;
    }
	size_t add_node(uint32_t w, uint32_t r, size_t l)
	{
		wd[node_i].index_of_W = w;
		wd[node_i].index_of_R = r;
		wd[node_i].len = l;
		return ++node_i;
	}

    bool check_node(size_t i, size_t ri, size_t node_i_)
    {
        if (node_i_ > 1 && i > wd[node_i_-1].index_of_W && (i - wd[node_i_-1].index_of_W) == wd[node_i_-1].len - PointerListLen + 1\
            && ri > wd[node_i_-1].index_of_R && (ri - wd[node_i_-1].index_of_R) == wd[node_i_-1].len - PointerListLen + 1)
        {
            ++wd[node_i_-1].len;
            return true;
        }
        return false;
    }

    bool check_node(size_t i, size_t ri, size_t node_i_, int z)
    {
        if (node_i_ > 1 && i > wd[node_i_-1].index_of_W && (i - wd[node_i_-1].index_of_W) == wd[node_i_-1].len - PointerListLen + 1\
            && ri > wd[node_i_-1].index_of_R && (ri - wd[node_i_-1].index_of_R) == wd[node_i_-1].len - PointerListLen + 1)
        {
            return true;
        }
        return false;
    }


    void print_log(char* dna_f, char* read, ofstream &outt)
    {
        for (size_t o = p_index, i; o>0; --o)
        {
            if (o==p_index) continue;
            i = path[o-1];
            outt << wd[i].index_of_W - PointerListLen + 1 << " " << wd[i].index_of_R - PointerListLen + 1<< " " << wd[i].len << "\n";
            cout_string_in_char(dna_f, wd[i].index_of_W + 1 - PointerListLen, wd[i].len, outt);
            outt << endl;
            cout_string_in_char(read, wd[i].index_of_R + 1 - PointerListLen, wd[i].len, outt);
            outt << endl; 
            outt << endl; 
        }
        // for (size_t i = 0; i<node_i; ++i)
     //    {
     //        outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
     //        // outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << endl; 
     //    }
     //    outt << endl;
    }

    void print_log(string& dna_f, string& read, ofstream &outt)
    {
        for (size_t o = p_index, i; o>0; --o)
        {
            if (o==p_index) continue;
            i = path[o-1];
            outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len) << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << endl; 
        }
        outt<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    	// for (size_t i = 0; i<node_i; ++i)
     //    {
     //        outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
     //        // outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << endl; 
     //    }
     //    outt << endl;
    }

    void find_path()
    {
    	dp = new uint32_t[node_i];
    	// memset(dp, 0, sizeof(dp));
        for (size_t i=0; i<node_i; ++i) dp[i] = 0;
    	for (size_t i=0; i<node_i-1; ++i)
    	{
    		for (size_t j=i+1; j<node_i; ++j)
    		{
                  // printf("%lu %lu %lu | ", (unsigned long)wd[j].len, (unsigned long)dp[i], (unsigned long)dp[j]);
                if (V[i][j]) dp[j] = (dp[j] > wd[j].len+dp[i] ? dp[j] : wd[j].len+dp[i]);
    		      // printf("%lu %lu %lu | ", (unsigned long)i, (unsigned long)j, (unsigned long)dp[j]);
            }
            // printf("\n");
        }
    	size_t thenode = node_i - 1;

    	//path from 1 to p_index) is the index of path string
    	size_t p_tmp;
    	p_index = 0;
    	path = new size_t[node_i];
    	while (thenode != 0)
    	{
    		p_tmp = thenode;
    		for (size_t i = 0; i< p_tmp; ++i)
    		{
    			if (i==thenode || V[i][thenode]==0) continue;
    			if (dp[thenode] == wd[thenode].len + dp[i])
    			{
    				thenode = i;
    				// PX(i);
    				path[p_index++] = i;
    				break;
    			}
    		}
    		if(thenode == p_tmp) break;
    	}
    }
    
    double do_alignment(string& dna_f, size_t window_up, size_t window_down, string&  read, string&ss, ofstream& outt)
    {
        // for (size_t o = p_index, i; o>0; --o)
        // {
        //     if (o==p_index) continue;
        //     i = path[o-1];
        //     outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        // }
        // outt<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

        size_t last_w=window_up + PointerListLen - 1, last_r=PointerListLen-1;
        double score=0;
        size_t index_ss=0;
        char snum[100];
        memset(snum, 0, sizeof(snum));
        size_t i, w, r, l;
        // printf("%lu %lu\n", (unsigned long)(p_index), (unsigned long)(window_down));
        if (p_index >= 3)
        {
            i = path[p_index-2];
            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            if (w!=last_w && r!=last_r)
            {
                // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));
                score+=get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, w-last_w, r-last_r, w-last_w, r-last_r, ss);
            }
            last_w = w;
            last_r = r;
        }
        // printf("2333\n");

        for (size_t o = p_index; o>0; --o)
        {
            if (o == p_index) continue;
            i = path[o-1];

            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            l = wd[i].len;


            if (w!=last_w || r!=last_r)
            {
                if (w!=last_w && r!=last_r)
                {
                    // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w), (unsigned long)(w), (unsigned long)(last_r), (unsigned long)(r));
                    // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));
                    score+=get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, 0, 0, w-last_w, r-last_r, ss);
                }
                else if(w==last_w)
                {
                    score+=(r-last_r)*scg;
                    sprintf(snum, "%lu", (unsigned long)(r-last_r));
                    ss.append(snum);
                    sprintf(snum, "I");
                    ss.append(snum);
                }
                else if(r==last_r)
                {
                    score+=(w-last_w)*scg;
                    sprintf(snum, "%lu", (unsigned long)(w-last_w));
                    ss.append(snum);
                    sprintf(snum, "D");
                    ss.append(snum);
                }
            }
            last_w = w+l;
            last_r = r+l;
            score+=(l)*scy;
            sprintf(snum, "%lu", (unsigned long)l);
            ss.append(snum);
            sprintf(snum, "M");
            ss.append(snum);
            // outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        }

        // printf("2333\n");
        if (p_index >= 3)
        {
            w = window_down;
            r = read.size();
            if (w!=last_w && r!=last_r)
            {
                score += get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, 0, 0, 0, 0, ss);
            }
            last_w = w;
            last_r = r;
        }
        // outt << ss << endl;
        printf("%lf\n", score);
        return score;
    }

    double do_alignment(char* dna_f, size_t window_up, size_t window_down, char*  read, size_t read_len, Cigar *ss, FILE* outt)
    {
        // for (size_t o = p_index, i; o>0; --o)
        // {
        //     if (o==p_index) continue;
        //     i = path[o-1];
        //     outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        // }
        // outt<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

        size_t last_w=window_up + PointerListLen - 1, last_r=PointerListLen-1;
        double score=0;
        size_t index_ss=0;
        char snum[100];
        memset(snum, 0, sizeof(snum));
        size_t i, w, r, l;
        // printf("%lu %lu\n", (unsigned long)(p_index), (unsigned long)(window_down));
        // if (p_index >= 3)
        {
            i = path[p_index-2];
            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            // PX(i);
            //     printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));

            if (w!=last_w && r!=last_r)
            {
                // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));
                // score+=get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, w-last_w, r-last_r, w-last_w, r-last_r, ss);
                score+=get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, w-last_w, r-last_r, w-last_w, r-last_r, ss);
            }
            last_w = w;
            last_r = r;
        }

        #ifdef PRINTALN
        cout << "A: ";
        for (size_t i = 0; i<(ss->length); ++i)
        {
            cout << cigar_int_to_len(ss->seq[i]) << cigar_int_to_op(ss->seq[i]);
        }
        cout << endl;
        #endif

        // printf("2333\n");

        // if (p_index >= 3)
        for (size_t o = p_index; o>0; --o)
        {
            if (o == p_index) continue;
            i = path[o-1];

            // PX(i);
            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            l = wd[i].len;

            // printf("%d %d %d %d\n, ", (int)last_w, (int)w, (int)last_r, (int)r);

            if (w!=last_w || r!=last_r)
            {
                if (w!=last_w && r!=last_r)
                {
                    printf("%lu %lu %lu %lu\n", (unsigned long)(last_w), (unsigned long)(w), (unsigned long)(last_r), (unsigned long)(r));
                    printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));
                    score+=get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, 0, 0, w-last_w, r-last_r, ss);
                }
                else if(w==last_w)
                {
                    score+=(r-last_r)*scg;
                    // sprintf(snum, "%lu", (unsigned long)(r-last_r));
                    // strcat(ss, snum);
                    // strcat(ss, "I");
                    ss->seq[(ss->length)++] = to_cigar_int((uint32_t)(r-last_r), 'I');
                }
                else if(r==last_r)
                {
                    score+=(w-last_w)*scg;
                    // sprintf(snum, "%lu", (unsigned long)(w-last_w));
                    // strcat(ss, snum);
                    // strcat(ss, "D");
                    ss->seq[(ss->length)++] = to_cigar_int((uint32_t)(w-last_w), 'D');

                }

                #ifdef PRINTALN
                for (size_t i = 0; i<(ss->length); ++i)
                {
                    cout << cigar_int_to_len(ss->seq[i]) << cigar_int_to_op(ss->seq[i]);
                }
                cout << endl;
                #endif



            }
            last_w = w+l;
            last_r = r+l;
            // printf("%d %d %d %d\n, ", (int)last_w, (int)w, (int)last_r, (int)r);

            score+=(l)*scy;
            // sprintf(snum, "%lu", (unsigned long)l);
            // strcat(ss, snum);
            // strcat(ss, "M");
            ss->seq[(ss->length)++] = to_cigar_int((uint32_t)(l), 'M');

            // outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
            #ifdef PRINTALN
            cout << "B: ";
            for (size_t i = 0; i<(ss->length); ++i)
            {
                cout << cigar_int_to_len(ss->seq[i]) << cigar_int_to_op(ss->seq[i]);
            }
            cout << endl;
            #endif
        }

     


        // // printf("2333\n");
        // if (p_index >= 3)
        {
            w = window_down + PointerListLen - 1;
            r = read_len + PointerListLen - 1;
            #ifdef PRINTALN
            printf("%d %d %d %d\n, ", (int)last_w, (int)w, (int)last_r, (int)r);
            #endif
            if (w!=last_w && r!=last_r)
            {
                score += get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, 0, 0, 0, r-last_r, ss);
            }
            last_w = w;
            last_r = r;
        }
        #ifdef PRINTALN
        cout << "C: ";
        for (size_t i = 0; i<(ss->length); ++i)
        {
            cout << cigar_int_to_len(ss->seq[i]) << cigar_int_to_op(ss->seq[i]);
        }
        cout << endl;
        #endif
        refine_cigar(ss);
        #ifdef PRINTALN
        for (size_t i = 0; i<(ss->length); ++i)
        {
            cout << cigar_int_to_len(ss->seq[i]) << cigar_int_to_op(ss->seq[i]);
        }
        cout << endl;
        #endif

        // outt << ss << endl;
        // printf("%lf\n", score);
        return score;

    }
    double do_alignment(char* dna_f, size_t window_up, size_t window_down, char*  read, size_t read_len, char* ss, FILE* outt)
    {
        // for (size_t o = p_index, i; o>0; --o)
        // {
        //     if (o==p_index) continue;
        //     i = path[o-1];
        //     outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        // }
        // outt<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

        size_t last_w=window_up + PointerListLen - 1, last_r=PointerListLen-1;
        double score=0;
        size_t index_ss=0;
        char snum[100];
        memset(snum, 0, sizeof(snum));
        size_t i, w, r, l;
        printf("%lu %lu\n", (unsigned long)(p_index), (unsigned long)(window_down));
        if (p_index >= 3)
        {
            i = path[p_index-2];
            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            if (w!=last_w && r!=last_r)
            {
                // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));
                score+=get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, w-last_w, r-last_r, w-last_w, r-last_r, ss);
            }
            last_w = w;
            last_r = r;
        }
        // printf("2333\n");

        for (size_t o = p_index; o>0; --o)
        {
            if (o == p_index) continue;
            i = path[o-1];

            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            l = wd[i].len;


            if (w!=last_w || r!=last_r)
            {
                if (w!=last_w && r!=last_r)
                {
                    // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w), (unsigned long)(w), (unsigned long)(last_r), (unsigned long)(r));
                    // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));
                    score+=get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, 0, 0, w-last_w, r-last_r, ss);
                }
                else if(w==last_w)
                {
                    score+=(r-last_r)*scg;
                    sprintf(snum, "%lu", (unsigned long)(r-last_r));
                    strcat(ss, snum);
                    strcat(ss, "I");
                }
                else if(r==last_r)
                {
                    score+=(w-last_w)*scg;
                    sprintf(snum, "%lu", (unsigned long)(w-last_w));
                    strcat(ss, snum);
                    strcat(ss, "D");
                }
            }
            last_w = w+l;
            last_r = r+l;
            score+=(l)*scy;
            sprintf(snum, "%lu", (unsigned long)l);
            strcat(ss, snum);
            strcat(ss, "M");
            // outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        }

        // printf("2333\n");
        if (p_index >= 3)
        {
            w = window_down;
            r = read_len;
            if (w!=last_w && r!=last_r)
            {
                score += get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, 0, 0, 0, 0, ss);
            }
            last_w = w;
            last_r = r;
        }
        // outt << ss << endl;
        printf("%lf\n", score);
        return score;
    } 
    double do_alignment(char* ref, size_t window_up, size_t window_down, char*  read, size_t read_l, string& thecigar, FILE* outt)
    {
        
        int k=0;
        for (int i=0; i<5; ++i)
        {
            for (int j=0; j<5; ++j)
            {
                if (i<4 && j<4) mat[k++] = i == j? match : mism;
                else mat[k++]=0;
            }
        }


        


        int read_len;
        int ref_len;
        int w_;

        uint32_t countM = 0;
        char     trans_cigar[500];
        int startPosCigar = 0;

        int qlen = 0;
        int tlen = 0;


        double score=0;

        size_t last_w=window_up + PointerListLen - 1, last_r=PointerListLen-1;
        size_t i, w, r, l;

        if (p_index < 3) return 0;
        if (p_index >= 3)
        {
            i = path[p_index-2];
            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            l = wd[i].len;
            ref_len = w-last_w;
            read_len = r-last_r;

            if (ref_len > PointerListLen / 2)
            {
                size_t tmp = ref_len - PointerListLen/2;
                ref_len = PointerListLen / 2;
                last_w += tmp;
            }
            if (read_len > PointerListLen / 2)
            {
                size_t tmp = read_len - PointerListLen/2;
                read_len = PointerListLen / 2;
                last_r += tmp;
            }

            // cout << string(read, last_r - PointerListLen + 1, read_len) << " " << string(ref, last_w - PointerListLen + 1, ref_len) << endl;
            if (w!=last_w && r!=last_r)
            {
                transIntoDec(readqry,read + last_r - PointerListLen + 1,read_len);
                transIntoDec(refqry,ref + last_w - PointerListLen + 1,ref_len);
                readqry_ = readqry;
                refqry_ = refqry;
                // cout << r << " " << last_r << " " << ref_len << " " << read_len << endl;
                w_ = max(ref_len, read_len);
                n_cigar = 0;
                // score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
                score += ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len , &qlen, &tlen, &cigar, &n_cigar) - read_len;
                if (n_cigar) {
                    for (int z=0;z<n_cigar-1;++z) {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    }

                    if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
                        countM = cigar[n_cigar-1] >> 4;
                        //startPosCigar += sprintf(trans_cigar,"%u%c",countM,'M');
                        //sams[countSam].cigar.append(trans_cigar);
                    } else {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
                        thecigar.append(trans_cigar);
                    }

                } 
                free(cigar);
            }
            score+=(l)*match;
            // cout << "---" << countM << endl;
            countM+=l;
            // cout << "---" << countM << endl;
            last_w = w-PointerListLen+l+1;
            last_r = r-PointerListLen+l+1;
        }

        for (size_t o = p_index-2; o>0; --o)
        {
            if (o == p_index) continue;
            i = path[o-1];

            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            l = wd[i].len;

            

            ref_len = w-last_w-PointerListLen+1;
            read_len = r-last_r-PointerListLen+1;


            // cout << last_r << " " << r << " " << last_w << " " << w << " " << l << endl;
            // cout << "----\n" << string(read, last_r, read_len) << " " << string(ref, last_w, ref_len) << endl;

            if (w>last_w+PointerListLen-1|| r>last_r+PointerListLen-1)
            {
                if (ref_len > 2048 || read_len > 2048) 
                {
                    return -INF;
                    cout << "hehe\t" << ref_len << " " << read_len << endl;
                }

                transIntoDec(readqry,read + last_r,read_len);
                transIntoDec(refqry,ref + last_w, ref_len);
                readqry_ = readqry;
                refqry_ = refqry;
                w_ = max(ref_len, read_len);
                n_cigar = 0;
                score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
                if (n_cigar - 1) 
                {
                    if (correspondTable[cigar[0]&0xf] == 'M') {
                        countM += (cigar[0] >> 4);
                        startPosCigar += sprintf(trans_cigar,"%uM",countM);
                        thecigar.append(trans_cigar);
                    }
                    else
                    {
                        startPosCigar += sprintf(trans_cigar,"%uM",countM);
                        thecigar.append(trans_cigar);
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    } 
                    countM = 0;
                    for (int z=1;z<n_cigar-1;++z) {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    }

                    if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
                        countM = cigar[n_cigar-1] >> 4;
                        //startPosCigar += sprintf(trans_cigar,"%u%c",countM,'M');
                        //sams[countSam].cigar.append(trans_cigar);
                    } else {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
                        thecigar.append(trans_cigar);
                    }

                } 
                else if (n_cigar>0)
                {
                    if (correspondTable[cigar[0]&0xf] == 'M') 
                        countM += (cigar[0] >> 4);
                    else {
                        startPosCigar += sprintf(trans_cigar,"%uM",countM);
                        thecigar.append(trans_cigar);
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                        countM = 0;

                    }
                }
                free(cigar);
            }

            score+=(l)*match;
            // cout << "---" << countM << endl;
            countM+=l;
            // cout << "---" << countM << endl;
            last_w = w-PointerListLen+l+1;
            last_r = r-PointerListLen+l+1;
            // outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        }

        // printf("****\n");
        if (p_index >= 3)
        {
            w = window_down;
            r = read_l;

            ref_len = w-last_w;
            read_len = r-last_r;

            if (ref_len > PointerListLen/2) ref_len = PointerListLen / 2;
            if (read_len > PointerListLen/2) read_len = PointerListLen / 2;

            // cout << last_r << " " << r << " " << last_w << " " << w << " " << l << endl;
            // cout << "----\n" << string(read, last_r, read_len) << " " << string(ref, last_w, ref_len) << endl;

            if (w>last_w && r>last_r)
            {
                transIntoDec(readqry,read + last_r,read_len);
                transIntoDec(refqry,ref + last_w, ref_len);
                readqry_ = readqry;
                refqry_ = refqry;
                w_ = max(ref_len, read_len);
                n_cigar = 0;
                qlen = 0;
                tlen = 0;
                score += ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len , &qlen, &tlen, &cigar, &n_cigar) - read_len;
                // score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
                if (n_cigar) {
                    if (correspondTable[cigar[0]&0xf] == 'M') {
                        countM += (cigar[0] >> 4);
                        startPosCigar += sprintf(trans_cigar,"%uM",countM);
                        thecigar.append(trans_cigar);
                    }
                    else
                    {
                        startPosCigar += sprintf(trans_cigar,"%uM",countM);
                        thecigar.append(trans_cigar);
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    } 
                    countM = 0;
                    for (int z=1;z<n_cigar;++z) {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    }

                } 
                free(cigar);
            }
        }

        // printf("%lf\n", score);
        return score;
    }

    
    void create_matrix()
    {
        V = new size_t*[node_i];
        for (size_t i=0; i<node_i; ++i) V[i] = new size_t[node_i];
        uint32_t iw, ir, jw, jr;
    	size_t il;

        // for (size_t i=0; i<node_i; ++i)
        // {
        //     iw = wd[i].index_of_W;
        //     ir = wd[i].index_of_R;
        //     il = wd[i].len;
        //     printf("%lu %lu %lu ", (unsigned long)iw, (unsigned long)ir, (unsigned long)il);
        //     printf("\n");
        // }


        for (size_t i=0; i<node_i-1; ++i)
        {
            iw = wd[i].index_of_W;
            ir = wd[i].index_of_R;
            il = wd[i].len;
            for (size_t j=i+1; j<node_i; ++j)
            {
                V[i][j]=0;
                jw = wd[j].index_of_W;
                jr = wd[j].index_of_R;
                // printf("%lu %lu %lu ", (unsigned long)jw, (unsigned long)(il + iw ), (unsigned long)iw);
                // printf("%lu %lu %lu ", (unsigned long)jr, (unsigned long)(il + ir ), (unsigned long)ir);

                if (jw >= il + iw && jr >= il  + ir && jr <= ir + il + t_wait)
                {
                    V[i][j] = 1;
                }
                // printf("(%lu %lu) %lu ", (unsigned long)i, (unsigned long)j, (unsigned long)V[i][j]);
            }
            // printf("\n");
        }
    }
};