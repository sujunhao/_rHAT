#ifndef DAG_H
#define DAG_H

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

// #define PRINTALN

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
    int32_t **V;
	char **Link;
	int32_t *dp;
	size_t node_i;
	size_t *path;
	size_t p_index;
    size_t _max_num;

public:
	DAG(size_t max_num)
	{
        _max_num = max_num;
        wd = (WINDOW_NODE *)malloc(max_num * sizeof(WINDOW_NODE));

        V = (int32_t **)malloc(max_num * sizeof(int32_t *));
        for (size_t i=0; i<max_num; ++i) V[i] = (int32_t *)malloc(max_num * sizeof(int32_t));

        Link = (char **)malloc(max_num * sizeof(char *));
        for (size_t i=0; i<max_num; ++i) Link[i] = (char *)malloc(max_num * sizeof(char));

        dp = (int32_t *)malloc(max_num * sizeof(int32_t));

        path = (size_t *)malloc(max_num * sizeof(size_t));

        
		node_i = 0;
	}
	~DAG()
	{
		if (V != NULL)
		{
			for (size_t i=0; i<_max_num; ++i) free(V[i]);
	        free(V);
		}
        if (Link != NULL)
        {
            for (size_t i=0; i<_max_num; ++i) free(Link[i]);
            free(Link);
        }
		if (dp != NULL)	free(dp);
		if (path != NULL) free(path);
  //       free(wd);
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
  
    

    double do_alignment(char* dna_f, size_t window_up, size_t window_down, char*  read, size_t read_len, Cigar *ss, size_t *first_pos, FILE* outt)
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
        // size_t index_ss=0;
        // char snum[100];
        // memset(snum, 0, sizeof(snum));
        size_t i, w, r, l;
        // printf("%lu %lu\n", (unsigned long)(p_index), (unsigned long)(window_down));
        // if (p_index >= 3)

        i = path[p_index-2];
        w = wd[i].index_of_W;
        r = wd[i].index_of_R;
        // PX(i);
        //     printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));


        if (w!=last_w && r!=last_r)
        {
            // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));
            // score+=get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, w-last_w, r-last_r, w-last_w, r-last_r, ss);
            score+=get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, w-last_w, 0, w-last_w, r-last_r, ss);
        }
        last_w = w;
        last_r = r;

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
                    // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w), (unsigned long)(w), (unsigned long)(last_r), (unsigned long)(r));
                    // printf("%lu %lu %lu %lu\n", (unsigned long)(last_w - PointerListLen + 1), (unsigned long)(w-last_w), (unsigned long)(last_r - PointerListLen + 1), (unsigned long)(r-last_r));
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

     

        size_t pos_len = 0;
        for (size_t i = 0; i<(ss->length); ++i)
        {
            if (cigar_int_to_op(ss->seq[i]) == 'I') continue;
            pos_len += cigar_int_to_len(ss->seq[i]);
        }
        (*first_pos) = last_w - PointerListLen + 1 - pos_len;

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
    
    

    
    // void create_matrix()
    // {
    //     V = new double*[node_i];
    //     for (size_t i=0; i<node_i; ++i) V[i] = new double[node_i];
    //     uint32_t iw, ir, jw, jr;
    // 	size_t il, jl;
    //     double tw, tr, dis;
    //     // for (size_t i=0; i<node_i; ++i)
    //     // {
    //     //     iw = wd[i].index_of_W;
    //     //     ir = wd[i].index_of_R;
    //     //     il = wd[i].len;
    //     //     printf("%lu %lu %lu ", (unsigned long)iw, (unsigned long)ir, (unsigned long)il);
    //     //     printf("\n");
    //     // }


    //     for (size_t i=0; i<node_i-1; ++i)
    //     {
    //         iw = wd[i].index_of_W;
    //         ir = wd[i].index_of_R;
    //         il = wd[i].len;
    //         for (size_t j=i+1; j<node_i; ++j)
    //         {
    //             V[i][j]=0;
    //             jw = wd[j].index_of_W;
    //             jr = wd[j].index_of_R;
    //             jl = wd[j].len;

    //             // printf("%lu %lu %lu ", (unsigned long)jw, (unsigned long)(il + iw ), (unsigned long)iw);
    //             // printf("%lu %lu %lu ", (unsigned long)jr, (unsigned long)(il + ir ), (unsigned long)ir);

    //             // if (jw >= il + iw && jr >= il  + ir && jr <= ir + il + t_wait)
    //             // {
    //             //     V[i][j] = 1;
    //             // }
    //             if (i==0 || j==(node_i-1))
    //             {
    //                 V[i][j] = 1 + jl;
    //             }
    //             else if (jw >= il + iw && jr >= il  + ir)
    //             {
    //                 tw = abs(1.0 * jw - iw);
    //                 tr = abs(1.0 * jr - ir);
    //                 dis = abs(tw - tr);
    //                 V[i][j] = 1 + jl + scg * (1 + dis);
    //             }

    //             // printf("(%lu %lu) %lu ", (unsigned long)i, (unsigned long)j, (unsigned long)V[i][j]);
    //         }
    //         // printf("\n");
    //     }
    // }

    void find_path()
    {
        uint32_t iw, ir, jw, jr;
        size_t il, jl;
        int32_t tw, tr, dis, weight;

        // memset(dp, 0, sizeof(dp));
        for (size_t i=0; i<node_i; ++i) dp[i] = 0;

        // printf("node len: %lu\n", (unsigned long)node_i);
        for (size_t i=0; i<node_i-1; ++i)
        {
            iw = wd[i].index_of_W;
            ir = wd[i].index_of_R;
            il = wd[i].len;
            for (size_t j=i+1; j<node_i; ++j)
            {
                V[i][j]=0;

                //for V[i][j] may == 0, so should have link to mark the point is linked
                Link[i][j] = 0;
                jw = wd[j].index_of_W;
                jr = wd[j].index_of_R;
                jl = wd[j].len;
                if (i==0 || j==(node_i-1))
                {
                    // V[i][j] = 1 + jl;
                    Link[i][j] = 1;
                    weight = (int32_t)(scy * jl);
                    V[i][j] = weight;
                    dp[j] = ((dp[j] > weight+dp[i]) ? dp[j] : weight+dp[i]);
                }
                else if (jw >= il + iw && jr >= il  + ir)
                {
                    Link[i][j] = 1;
                    tw = abs(1.0 * jw - iw);
                    tr = abs(1.0 * jr - ir);
                    dis = abs(tw - tr);
                    weight = (int32_t)(1 + scy * jl + scg * (1 + dis)) + 1;
                    V[i][j] = weight;
                    dp[j] = (dp[j] > weight+dp[i] ? dp[j] : weight+dp[i]);
                }
                // printf("(%lu,%lu) %d ", (unsigned long)i, (unsigned long)j, weight);
            }
            // printf("---\n");
        }
        // printf("%d\n", dp[node_i-1]);

        // for (size_t i=0; i<node_i-1; ++i)
        // {
        //     for (size_t j=i+1; j<node_i; ++j)
        //     {
        //           // printf("%lu %lu %lu | ", (unsigned long)wd[j].len, (unsigned long)dp[i], (unsigned long)dp[j]);
        //         if (V[i][j]) dp[j] = (dp[j] > wd[j].len+dp[i] ? dp[j] : wd[j].len+dp[i]);
        //           // printf("%lu %lu %lu | ", (unsigned long)i, (unsigned long)j, (unsigned long)dp[j]);
        //     }
        //     // printf("\n");
        // }
        size_t thenode = node_i - 1;
        //path from 1 to p_index) is the index of path string
        size_t p_tmp;
        p_index = 0;
        // path = new size_t[node_i];
        while (thenode != 0)
        {
            p_tmp = thenode;
            for (size_t i = 0; i< p_tmp; ++i)
            {
                if (i==thenode || Link[i][thenode]==0) continue;
                if (dp[thenode] == V[i][thenode] + dp[i])
                {
                    thenode = i;
                    // PX(i);
                    path[p_index++] = i;
                    break;
                }
            }
            if(thenode == p_tmp) break;
        }


        // for (size_t i=0; i<node_i; ++i) free(V[i]);
        // free(V);
        // free(dp);
    }

};




#endif