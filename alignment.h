#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <cstring>
#include <cmath>
#include <queue>
#include <unistd.h>
#include <cstring>
#include <stdint.h>
#include "cigar.h"
#include "LV_deep.h"
using namespace std;

double scy = 1, scn = -2, scg = -3;
double INF = 0xFFFFFFFF;

double zero = 1e-16;
char sstate[]="MDI";

size_t max_dp_ref_len = 15000, max_dp_read_len = 15000;

const uint32_t mx_len=15000;

// double dp[mx_len][mx_len];
// int v[mx_len];
double **dp;
int *v;


// int max_error = 50, max_pattern=5000;
LV_ENTITY *lv;


//return score 
//input the first string, string len, the second, the start point the end point the alignment out string and its len

double get_alignment(char* s1, size_t last_w, size_t lw, char* s2, size_t last_r, size_t lr, size_t si, size_t sj, size_t ei, size_t ej, Cigar *ss, int log=0)
{
    // printf(">%d %d\n", (int)lw, (int)lr);
    if (log)
        printf(">%.*s\n>%.*s\n", (int)lw, s1+last_w, (int)lr, s2+last_r);
    if (lw >= max_dp_ref_len || lr >= max_dp_read_len) return 0;
    // if (lw >= mx_len || lr>=mx_len) return 0;
    size_t offset_1 = last_w, offset_2 = last_r;
    size_t n_s1 = lw, n_s2 = lr, n_s;
    size_t max_n = max(n_s1, n_s2) * 2;

    size_t xi=0, xj=0, index_v=0;
    double dp_max=-INF;


    double tmp = 0.0;
    
    for (size_t i = 0; i <= n_s1; ++i)
    {
        for (size_t j = 0; j <= n_s2; ++j)
        {
            //set the init dp score
            if (i<=si && j<=sj) dp[i][j] = 0;
            else dp[i][j] = -INF;

            //set_dp
            tmp = dp[i][j];
            if (i) tmp = std::max(tmp, dp[i-1][j] + scg);
            if (j) tmp = std::max(tmp, dp[i][j-1] + scg);
            if (i > 0 && j > 0) tmp = std::max(tmp, dp[i-1][j-1] + ((s1[offset_1 + i-1]==s2[offset_2 + j-1]) ? scy : scn));
            dp[i][j] = tmp;
            // if (log)
            //     printf("%3.0lf ", dp[i][j]);

            // find the last dp position to trace
            if (i>=ei && j>= ej && dp[i][j]>dp_max)
            {
                dp_max = dp[i][j];
                xi = i;
                xj = j;
            }
        }
        // if (log)
        //     printf("\n");
    }


    while (abs(dp[xi][xj])>1e-16 || xi>si || xj>sj)
    {
        if (xi && xj && dp[xi][xj] == (dp[xi-1][xj-1] + ((s1[offset_1+xi-1]==s2[offset_2+xj-1]) ? scy : scn))) 
        {
            --xi;
            --xj;
            v[index_v++] = 0;
            continue;
        }
        if (xj && dp[xi][xj] == dp[xi][xj-1] + scg)
        {
            --xj;
            v[index_v++] = 2;
            continue;
        }
        if (xi && dp[xi][xj] == dp[xi-1][xj] + scg)
        {
            --xi;
            v[index_v++] = 1;
            continue;
        }
    }

    // printf("%lu 2333sadasd\n", (unsigned long)index_v);
    //if index_v == 0 will error
    if (index_v==0) return 0;
    size_t i=index_v-1, m=offset_1+xi, n=offset_2+xj;
    int state=-1;
    uint32_t cnt=0;

    for (; i>=0; --i) 
    {
        // printf(i==0?"%d\n":"%d ", v[i]);

        if (state==v[i]) ++cnt;
        else
        {
            if (state>=0)
            {
                ss->seq[(ss->length)++] = to_cigar_int(cnt, sstate[state]);
            }
            state=v[i];
            cnt=1;
        }
        if (i==0) break;
    }
    if (state>=0)
        ss->seq[(ss->length)++] = to_cigar_int(cnt, sstate[state]);    

    // printf("%f\n", dp_max);
    return dp_max;
} 


void run_deep_lv_cigar(LV_ENTITY *lv, DI *di, char *p, size_t lp, int max_error, Cigar *ss);


//
double get_alignment_lv(ref_s *rs, var_s *vs, size_t last_w, size_t lw, char* s2_4bit, size_t last_r, size_t lr, size_t si, size_t sj, size_t ei, size_t ej, Cigar *ss)
{
    //shoud init LV
    ////num of node, m_e, m_p, num of deep
    //init_lv_space(&lv, 100, max_error, max_pattern, 100);

    // printf(">%d %d\n", (int)lw, (int)lr);
    
    int error = 20;
    double score = 0;

    DI di;
    size_t begin = last_w, end = last_w+lw;
    run_deep_lv_cigar(lv, &di, s2_4bit+last_r, lr, error, ss);
    free_node_space(&di);


    return score;

} 


void run_deep_lv_cigar(LV_ENTITY *lv, DI *di, char *p, size_t lp, int max_error, Cigar *ss)
{
    ;
}