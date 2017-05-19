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
using namespace std;

double scy = 1, scn = -2, scg = -2;
double INF = 0xFFFFFFFF;
double zero = 1e-16;
const uint32_t mx_len=15000;
double dp[mx_len][mx_len];
int v[mx_len];


#define MAPSTR "MIDNSHP=X"
#ifndef BAM_CIGAR_SHIFT
#define BAM_CIGAR_SHIFT 4u
#endif


typedef struct {
    uint32_t* seq;
    int32_t length;
} Cigar;

/*!     @function               Produce CIGAR 32-bit unsigned integer from CIGAR operation and CIGAR length
        @param  length          length of CIGAR
        @param  op_letter       CIGAR operation character ('M', 'I', etc)
        @return                 32-bit unsigned integer, representing encoded CIGAR operation and length
*/
uint32_t to_cigar_int (uint32_t length, char op_letter)
{
        switch (op_letter) {
                case 'M': /* alignment match (can be a sequence match or mismatch */
                default:
                        return length << BAM_CIGAR_SHIFT;
                case 'S': /* soft clipping (clipped sequences present in SEQ) */
                        return (length << BAM_CIGAR_SHIFT) | (4u);
                case 'D': /* deletion from the reference */
                        return (length << BAM_CIGAR_SHIFT) | (2u);
                case 'I': /* insertion to the reference */
                        return (length << BAM_CIGAR_SHIFT) | (1u);
                case 'H': /* hard clipping (clipped sequences NOT present in SEQ) */
                        return (length << BAM_CIGAR_SHIFT) | (5u);
                case 'N': /* skipped region from the reference */
                        return (length << BAM_CIGAR_SHIFT) | (3u);
                case 'P': /* padding (silent deletion from padded reference) */
                        return (length << BAM_CIGAR_SHIFT) | (6u);
                case '=': /* sequence match */
                        return (length << BAM_CIGAR_SHIFT) | (7u);
                case 'X': /* sequence mismatch */
                        return (length << BAM_CIGAR_SHIFT) | (8u);
        }
        return (uint32_t)-1; // This never happens
}

/*! @function       Extract CIGAR operation character from CIGAR 32-bit unsigned integer
    @param  cigar_int   32-bit unsigned integer, representing encoded CIGAR operation and length
    @return         CIGAR operation character ('M', 'I', etc)
*/
//char cigar_int_to_op (uint32_t cigar_int);
static inline char cigar_int_to_op(uint32_t cigar_int) 
{
    return (cigar_int & 0xfU) > 8 ? 'M': MAPSTR[cigar_int & 0xfU];
}


/*! @function       Extract length of a CIGAR operation from CIGAR 32-bit unsigned integer
    @param  cigar_int   32-bit unsigned integer, representing encoded CIGAR operation and length
    @return         length of CIGAR operation
*/
//uint32_t cigar_int_to_len (uint32_t cigar_int);
static inline uint32_t cigar_int_to_len (uint32_t cigar_int)
{
    return cigar_int >> BAM_CIGAR_SHIFT;
}


void refine_cigar(Cigar *cg)
{
    uint32_t cnt = 0, tmpn;
    char c = cigar_int_to_op((cg->seq)[0]), tmpc;
    int32_t len = cg->length, nl = 0;
    for (int i=0; i<len; ++i)
    {
        tmpn = cigar_int_to_len((cg->seq)[i]);
        tmpc = cigar_int_to_op((cg->seq)[i]);
        if (tmpc == c)  cnt += tmpn;
        else 
        {
            (cg->seq)[nl++] = to_cigar_int(cnt, c);
            c = tmpc;
            cnt = tmpn;
        }
        // printf("%c", d_path[i]);
    }
    (cg->seq)[nl++] = to_cigar_int(cnt, c);
    (cg->length) = nl;
}
//return score 
//input the first string, string len, the second, the start point the end point the alignment out string and its len
double get_alignment(string &s1, size_t last_w, size_t lw, string &s2, size_t last_r, size_t lr, size_t si, size_t sj, size_t ei, size_t ej, string &ss)
{
    // std::cout << s1.substr(last_w, lw) << '\n' << s2.substr(last_r, lr) << '\n';

    size_t offset_1 = last_w, offset_2 = last_r;
    size_t n_s1 = lw, n_s2 = lr, n_s;
    double **dp;
    int *v;
    size_t max_n = max(n_s1, n_s2) * 2;
    dp = new double*[max_n];
    for (size_t i=0; i<max_n; ++i) dp[i]=new double[max_n];
    v = new int[max_n];

    //set the init dp score
    for (size_t i = 0; i <= n_s1; ++i)
        for (size_t j = 0; j <= n_s2; ++j)
            if (i<=si && j<=sj) dp[i][j] = 0;
            else dp[i][j] = -INF;
    //set_dp
    double tmp = 0.0;
    for (size_t i = 0; i <= n_s1; ++i)
    {
        for (size_t j = 0; j <= n_s2; ++j)
        {
            tmp = dp[i][j];
            if (i) tmp = std::max(tmp, dp[i-1][j] + scg);
            if (j) tmp = std::max(tmp, dp[i][j-1] + scg);
            if (i > 0 && j > 0) tmp = std::max(tmp, dp[i-1][j-1] + ((s1[offset_1 + i-1]==s2[offset_2 + j-1]) ? scy : scn));
            dp[i][j] = tmp;
            // printf("%3.0lf ", dp[i][j]);
        }
        // printf("\n");
    }
    // find the last dp position to trace
    size_t xi=0, xj=0, index_v=0;
    double dp_max=-INF;
    for (size_t i=ei; i<=n_s1; ++i)
        for (size_t j=ej; j<=n_s2; ++j)
            if (dp[i][j]>dp_max)
            {
                dp_max = dp[i][j];
                xi = i;
                xj = j;
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
    char snum[100];
    char sstate[]="MDI";
    unsigned long cnt=0;

    for (; i>=0; --i) 
    {
        // printf(i==0?"%d\n":"%d ", v[i]);

        if (state==v[i]) ++cnt;
        else
        {
            if (state>=0)
            {
                sprintf(snum, "%lu", cnt);
                ss.append(snum);
                ss.append(sstate+state, 1);
            }
            state=v[i];
            cnt=1;
        }
        if (i==0) break;
        // printf("%d ", v[i]);
        // printf("%lu 2333sadasd\n", (unsigned long)i);
    }
    if (state>=0)
    {
        sprintf(snum, "%lu", cnt);
        ss.append(snum);
        ss.append(sstate+state, 1);
    }
    
    // string s_1(index_v, 0), s_2(index_v, 0);
    // s_1.clear();
    // s_2.clear();
    // i=index_v-1; m=offset_1+xi; n=offset_2+xj;
    // for (; i>=0; --i) 
    // {
    //     printf(i==0?"%d\n":"%d ", v[i]);
    //     switch(v[i])
    //     {
    //         case 0:
    //             s_1.append(1, s1[m++]);
    //             s_2.append(1, s2[n++]);
    //             break;
    //         case 1:
    //             s_1.append(1, s1[m++]);
    //             s_2.append("_");
    //             break;
    //         case 2:
    //             s_1.append("_");
    //             s_2.append(1, s2[n++]);
    //             break;
    //     }
    //     if (i==0) break;
    //     // printf("%d ", v[i]);
    // }
    // std::cout << s_1 << '\n' << s_2 << '\n';

    for (size_t i=0; i<max_n; ++i) delete [] dp[i];
    delete [] dp;
    delete [] v;
    printf("%f\n", dp_max);
    return dp_max;
}

double get_alignment(char* s1, size_t last_w, size_t lw, char* s2, size_t last_r, size_t lr, size_t si, size_t sj, size_t ei, size_t ej, char* ss)
{
    // std::cout << s1.substr(last_w, lw) << '\n' << s2.substr(last_r, lr) << '\n';
    if (lw >= mx_len || lr>=mx_len) return 0;
    size_t offset_1 = last_w, offset_2 = last_r;
    size_t n_s1 = lw, n_s2 = lr, n_s;
    size_t max_n = max(n_s1, n_s2) * 2;

    //set the init dp score
    for (size_t i = 0; i <= n_s1; ++i)
        for (size_t j = 0; j <= n_s2; ++j)
            if (i<=si && j<=sj) dp[i][j] = 0;
            else dp[i][j] = -INF;
    //set_dp
    double tmp = 0.0;
    for (size_t i = 0; i <= n_s1; ++i)
    {
        for (size_t j = 0; j <= n_s2; ++j)
        {
            tmp = dp[i][j];
            if (i) tmp = std::max(tmp, dp[i-1][j] + scg);
            if (j) tmp = std::max(tmp, dp[i][j-1] + scg);
            if (i > 0 && j > 0) tmp = std::max(tmp, dp[i-1][j-1] + ((s1[offset_1 + i-1]==s2[offset_2 + j-1]) ? scy : scn));
            dp[i][j] = tmp;
            // printf("%3.0lf ", dp[i][j]);
        }
        // printf("\n");
    }
    // find the last dp position to trace
    size_t xi=0, xj=0, index_v=0;
    double dp_max=-INF;
    for (size_t i=ei; i<=n_s1; ++i)
        for (size_t j=ej; j<=n_s2; ++j)
            if (dp[i][j]>dp_max)
            {
                dp_max = dp[i][j];
                xi = i;
                xj = j;
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
    char snum[100];
    char sstate[]="MDI";
    unsigned long cnt=0;

    for (; i>=0; --i) 
    {
        // printf(i==0?"%d\n":"%d ", v[i]);

        if (state==v[i]) ++cnt;
        else
        {
            if (state>=0)
            {
                sprintf(snum, "%lu", cnt);
                strcat(ss, snum);
                strncat(ss, sstate+state, 1);
            }
            state=v[i];
            cnt=1;
        }
        if (i==0) break;
        // printf("%d ", v[i]);
        // printf("%lu 2333sadasd\n", (unsigned long)i);
    }
    if (state>=0)
    {
        sprintf(snum, "%lu", cnt);
        strcat(ss, snum);
        strncat(ss, sstate+state, 1);
    }
    
    // string s_1(index_v, 0), s_2(index_v, 0);
    // s_1.clear();
    // s_2.clear();
    // i=index_v-1; m=offset_1+xi; n=offset_2+xj;
    // for (; i>=0; --i) 
    // {
    //     printf(i==0?"%d\n":"%d ", v[i]);
    //     switch(v[i])
    //     {
    //         case 0:
    //             s_1.append(1, s1[m++]);
    //             s_2.append(1, s2[n++]);
    //             break;
    //         case 1:
    //             s_1.append(1, s1[m++]);
    //             s_2.append("_");
    //             break;
    //         case 2:
    //             s_1.append("_");
    //             s_2.append(1, s2[n++]);
    //             break;
    //     }
    //     if (i==0) break;
    //     // printf("%d ", v[i]);
    // }
    // std::cout << s_1 << '\n' << s_2 << '\n';

    // printf("%f\n", dp_max);
    return dp_max;
} 


double get_alignment(char* s1, size_t last_w, size_t lw, char* s2, size_t last_r, size_t lr, size_t si, size_t sj, size_t ei, size_t ej, Cigar *ss)
{
    // std::cout << s1.substr(last_w, lw) << '\n' << s2.substr(last_r, lr) << '\n';
    if (lw >= mx_len || lr>=mx_len) return 0;
    size_t offset_1 = last_w, offset_2 = last_r;
    size_t n_s1 = lw, n_s2 = lr, n_s;
    size_t max_n = max(n_s1, n_s2) * 2;

    //set the init dp score
    for (size_t i = 0; i <= n_s1; ++i)
        for (size_t j = 0; j <= n_s2; ++j)
            if (i<=si && j<=sj) dp[i][j] = 0;
            else dp[i][j] = -INF;
    //set_dp
    double tmp = 0.0;
    for (size_t i = 0; i <= n_s1; ++i)
    {
        for (size_t j = 0; j <= n_s2; ++j)
        {
            tmp = dp[i][j];
            if (i) tmp = std::max(tmp, dp[i-1][j] + scg);
            if (j) tmp = std::max(tmp, dp[i][j-1] + scg);
            if (i > 0 && j > 0) tmp = std::max(tmp, dp[i-1][j-1] + ((s1[offset_1 + i-1]==s2[offset_2 + j-1]) ? scy : scn));
            dp[i][j] = tmp;
            // printf("%3.0lf ", dp[i][j]);
        }
        // printf("\n");
    }
    // find the last dp position to trace
    size_t xi=0, xj=0, index_v=0;
    double dp_max=-INF;
    for (size_t i=ei; i<=n_s1; ++i)
        for (size_t j=ej; j<=n_s2; ++j)
            if (dp[i][j]>dp_max)
            {
                dp_max = dp[i][j];
                xi = i;
                xj = j;
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
    char snum[100];
    char sstate[]="MDI";
    uint32_t cnt=0;

    for (; i>=0; --i) 
    {
        // printf(i==0?"%d\n":"%d ", v[i]);

        if (state==v[i]) ++cnt;
        else
        {
            if (state>=0)
            {
                // sprintf(snum, "%lu", cnt);
                // strcat(ss, snum);
                // strncat(ss, sstate+state, 1);
                ss->seq[(ss->length)++] = to_cigar_int(cnt, sstate[state]);
                // printf("%d %c ", cnt, sstate[state]);
            }
            state=v[i];
            cnt=1;
        }
        if (i==0) break;
        // printf("%d ", v[i]);
        // printf("%lu 2333sadasd\n", (unsigned long)i);
    }
    if (state>=0)
    {
        // sprintf(snum, "%lu", cnt);
        // strcat(ss, snum);
        // strncat(ss, sstate+state, 1);
        ss->seq[(ss->length)++] = to_cigar_int(cnt, sstate[state]);


    }
    
    // string s_1(index_v, 0), s_2(index_v, 0);
    // s_1.clear();
    // s_2.clear();
    // i=index_v-1; m=offset_1+xi; n=offset_2+xj;
    // for (; i>=0; --i) 
    // {
    //     printf(i==0?"%d\n":"%d ", v[i]);
    //     switch(v[i])
    //     {
    //         case 0:
    //             s_1.append(1, s1[m++]);
    //             s_2.append(1, s2[n++]);
    //             break;
    //         case 1:
    //             s_1.append(1, s1[m++]);
    //             s_2.append("_");
    //             break;
    //         case 2:
    //             s_1.append("_");
    //             s_2.append(1, s2[n++]);
    //             break;
    //     }
    //     if (i==0) break;
    //     // printf("%d ", v[i]);
    // }
    // std::cout << s_1 << '\n' << s_2 << '\n';

    // printf("%f\n", dp_max);
    return dp_max;
} 
