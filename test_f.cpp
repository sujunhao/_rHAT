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

double scy = 1, scn = -3, scg = -2;
double INF = 0xFFFFFFFF;
double zero = 1e-16;


//return score 
//input the first string, string len, the second, the start point the end point the alignment out string and its len
double get_alignment(char* s1, size_t last_w, size_t lw, char* s2, size_t last_r, size_t lr, size_t si, size_t sj, size_t ei, size_t ej, char* ss)
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

    for (size_t i=0; i<max_n; ++i) delete [] dp[i];
    delete [] dp;
    delete [] v;
    printf("%f\n", dp_max);
    return dp_max;
}

int main()
{
    char s1[]="ATGCACGT", s2[]="ATGCCC";
    char s3[]="asdasdasdasdasdasdasdasdasd", s4[]="asdasdsadasasdsaddasdasdsadasdasdasd";
    string ss1, ss2, sss;
    char ss[1000];
    memset(ss, 0, sizeof(ss));
    size_t n=0;
    size_t k;

    while (scanf("%s%s", s1, s2))
    {
        get_alignment(s1, 0, strlen(s1), s2, 0, strlen(s2), 0, 0, strlen(s1), strlen(s2), ss);
        puts(ss);
        // get_alignment(s1, strlen(s1), s2, strlen(s2), 0, 0, 0, 0, ss, n);
        // get_alignment(s1, strlen(s1), s2, strlen(s2), strlen(s1), strlen(s2), strlen(s1), strlen(s2), ss, n);
    }

    // while (std::cin >> ss1 >> ss2)
    // {
    //     // get_alignment(ss1, ss2, 0, 0, ss1.size(), ss2.size(), sss);
    //     // get_alignment(ss1, ss2, 0, 0, 0, 0, sss);
    //     // get_alignment(ss1, ss2, ss1.size(), ss2.size(), ss1.size(), ss2.size(), sss);
    //     sss.clear();
    //     get_alignment(ss1, 0, ss1.size(), ss2, 0, ss2.size(), 0, 0, ss1.size(), ss2.size(), sss);

    // }
    
    return 0;

}