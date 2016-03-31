#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "RHT.h"
using namespace std;

// #define PRINT_WINDOW_INDEX
uint32_t z[]={0, 2, 2, 5, 6, 34, 55, 56, 67, 100};
uint32_t search(uint32_t p_index)
{
    //binary search
    uint32_t l=0, r=9, m;
    while (l < r)
    {
        m = (l+r)/2;
        if (p_index <= (z[m]))   r = m;
        else l = m + 1;
    }
    // printf("%lu %lu %lu\n", (unsigned long)l, (unsigned long)m, (unsigned long)r);
    if (z[r] == p_index) return r+1;
    return 0;
}
int main(int argc, char** argv) 
{

    unsigned long k;
    for (size_t i=0; i<1000; ++i)
        if (search(i)) printf("%lu\n", i);
    // while (~scanf("%lu", &k))
    // {
    //     uint32_t f = k;
    //     printf("%lu\n", (unsigned long)search(f));
    // }


    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);

    return 0;
    
}
