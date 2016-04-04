#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <queue>
#include <unistd.h>
#include <cstring>
#include <stdint.h>
using namespace std;

double scy = 1, scn = -1, scg = -2;

class ALIGNMENT
{
    private:
        double cy, cn, cg, INF;
        size_t L[2];
        char *s1, *s2;
        int *v;
        double **dp;
        size_t n_s1, n_s2, index_v;
        enum{LOCAL, GLOBAL, SEMI} alignment_state;
        bool semi_state[2];
    public:
        ALIGNMENT()
        {
            cy = scy; cn = scn; cg = scg;
            INF = 0xFFFFFFFF;
        }
        ~ALIGNMENT()
        {

        }

        size_t get_global_(char ss1[], char ss2[], char s3[], char s4[], size_t n_ss1, size_t n_ss2)
        {
            s1 = ss1;
            s2 = ss2;
            n_s1 = n_ss1;
            n_s2 = n_ss2;

            alignment_state = GLOBAL;

            size_t n_max = std::max(n_s1, n_s2) + 1;

            dp = new double*[n_max];
            for (size_t i=0; i<n_max; ++i) dp[i] = new double[n_max];
            set_dp_();

            for (size_t i=0; i<=n_s1; ++i)
            {
                for (size_t j=0; j<=n_s2; ++j)
                {
                    printf(j==n_s2?"%3.0lf\n":"%3.0lf ", dp[i][j]);
                }
            }

            v = new int[n_max<<1];
            index_v = 0;
            get_dp_(n_s1, n_s2);

            for (size_t i=0; i<index_v; ++i)
                printf(i==index_v-1?"%d----\n":"%d ", v[i]);

            for (size_t i=0, m=L[0], n=L[1] ; i<index_v; ++i)
            {
                switch(v[i])
                {
                    case 0:
                        s3[i] = s1[m++];
                        s4[i] = s2[n++];
                        break;
                    case 1:
                        s3[i] = s1[m++];
                        s4[i] = '_';
                        break;
                    case 2:
                        s3[i] = '_';
                        s4[i] = s2[n++];
                        break;

                }
            }
            s3[index_v] = '\0';
            s4[index_v] = '\0';
            puts(s3);
            puts(s4);

            return index_v;
        }
        size_t get_semi_(char ss1[], char ss2[], char s3[], char s4[], size_t n_ss1, size_t n_ss2)
        {

            s1 = ss1;
            s2 = ss2;
            n_s1 = n_ss1;
            n_s2 = n_ss2;

            alignment_state = SEMI;

            size_t n_max = std::max(n_s1, n_s2) + 1;

            dp = new double*[n_max];
            for (size_t i=0; i<n_max; ++i) dp[i] = new double[n_max];
            set_dp_();

            for (size_t i=0; i<=n_s1; ++i)
            {
                for (size_t j=0; j<=n_s2; ++j)
                {
                    printf(j==n_s2?"%3.0lf\n":"%3.0lf ", dp[i][j]);
                }
            }

            v = new int[n_max<<1];
            index_v = 0;
            size_t a = n_s1, b = n_s2;
            double k = dp[a][b];
            for (size_t i = 0; i < n_s1; ++i)
            {
                if (dp[i][n_s2] > k) 
                {
                    k = dp[i][n_s2];
                    a = i;
                    b = n_s2;
                }
            }
            for (size_t i = 0; i < n_s2; ++i)
            {
                if (dp[n_s1][i] > k) 
                {
                    k = dp[n_s1][i];
                    a = n_s1;
                    b = i;
                }
            }
            semi_state[0]=true;
            semi_state[1]=true;
            if (a<n_s1) semi_state[1]=false;
            if (b<n_s2) semi_state[0]=false;
            // cout << endl << k << endl;
            get_dp_(a, b);
            // get_dp_(n_s1-1, n_s2-1);

            for (size_t i=0; i<index_v; ++i)
                printf(i==index_v-1?"%d----\n":"%d ", v[i]);

            for (size_t i=0, m=L[0], n=L[1] ; i<index_v; ++i)
            {
                switch(v[i])
                {
                    case 0:
                        s3[i] = s1[m++];
                        s4[i] = s2[n++];
                        break;
                    case 1:
                        s3[i] = s1[m++];
                        s4[i] = '_';
                        break;
                    case 2:
                        s3[i] = '_';
                        s4[i] = s2[n++];
                        break;

                }
            }
            s3[index_v] = '\0';
            s4[index_v] = '\0';
            puts(s3);
            puts(s4);

            return index_v;
            
        }
        void set_dp_()
        {
            double tmp = 0.0;
            for (size_t i = 0; i <= n_s1; ++i)
            {
                for (size_t j = 0; j <= n_s2; ++j)
                {
                    if (alignment_state == GLOBAL)
                    {
                        dp[i][j] = -INF;
                    }
                    else if (alignment_state == SEMI)
                    {
                        if (i == 0) dp[i][j] = 0;
                        else if (j == 0) dp[i][j] = 0;
                        else dp[i][j] = -INF;
                    }
                    else
                    {
                        dp[i][j] = 0;
                    }
                    if (i==0 && j==0) dp[i][j] = 0;
                    tmp = dp[i][j];
                    if (i) tmp = std::max(tmp, dp[i-1][j] + cg);
                    if (j) tmp = std::max(tmp, dp[i][j-1] + cg);
                    if (i > 0 && j > 0) tmp = std::max(tmp, dp[i-1][j-1] + ((s1[i-1]==s2[j-1]) ? cy : cn));
                    dp[i][j] = tmp;
                    // printf("%3.0lf ", dp[i][j]);
                }
                // printf("\n");
            }
        }

        void get_dp_(size_t i, size_t j)
        {
            // printf("->%lu %lu %d\n", (unsigned long)i, (unsigned long)j, v[index_v]);
            if (alignment_state == GLOBAL)
            {
                if (i==0 && j==0)
                {
                    L[0] = i;
                    L[1] = j;
                    return;
                }
            }
            else if (alignment_state == SEMI)
            {
                if (i==0 && semi_state[1]==true)
                {
                    L[0] = i;
                    L[1] = j;
                    return;
                }
                if (j==0 && semi_state[0]==true)
                {
                    L[0] = i;
                    L[1] = j;
                    return;
                }
            }
            else
            {
                if (dp[i][j] == 0) 
                {
                    L[0] = i;
                    L[1] = j;
                    return;
                }
            }

            if (i && j && dp[i][j] == (dp[i-1][j-1] + ((s1[i-1]==s2[j-1]) ? cy : cn))) 
            {
                get_dp_(i-1, j-1);
                v[index_v++] = 0;
                return;
            }
            if (j && dp[i][j] == dp[i][j-1] + cg)
            {
                get_dp_(i, j-1);
                v[index_v++] = 2;
                return;
            }
            if (i && dp[i][j] == dp[i-1][j] + cg)
            {
                get_dp_(i-1, j);
                v[index_v++] = 1;
                return;
            }
        }
};

int main()
{
    char s1[]="ATGCACGT", s2[]="ATGCCC";
    char s3[]="asdasdasdasdasdasdasdasdasd", s4[]="asdasdsadasasdsaddasdasdsadasdasdasd";
    size_t n;

    while (scanf("%s%s", s1, s2))
    {
        ALIGNMENT A;
        n = A.get_global_(s1, s2, s3, s4, strlen(s1), strlen(s2));
        n = A.get_semi_(s1, s2, s3, s4, strlen(s1), strlen(s2));
        
    }
    
    return 0;

}