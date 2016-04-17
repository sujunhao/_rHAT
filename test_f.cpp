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
#include "ksw.h"
using namespace std;


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