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

int main(int argc, char** argv) 
{

    
    ifstream inf;
    inf.open("E.coli.fa");

    string dna_name, dna_s;
    const size_t tmp = 1024;
    char dna_w[tmp];
    string dna_w()

    uint64_t dna_ref_b = 0;
    size_t  index_s = 0;
    size_t  index_w = 0;

    getline(inf, dna_name);

    while (inf >> dna_s)
    {
        index_s = 0;
        while (index_s < dna_s.size())
        {
            ++dna_ref_b;
            dna_w[index_w] = dna_s[index_s]

            if (dna_ref_b > WindowListLen / 2 && ((dna_ref_b - 1 - WindowListLen / 2) % WindowListLen + 1 >= PointerListLen))
            {
                if ((dna_ref_b - 1) % WindowListLen + 1>= PointerListLen)
                {
                    if ((dna_ref_b / WindowListLen) * 2 < 2 * ((dna_ref_b - WindowListLen / 2 ) / WindowListLen) + 1)
                    {
                        db.link_string(dna_w.substr(dna_w.size() - PointerListLen, PointerListLen), ((dna_ref_b - 1) / WindowListLen) * 2);
                        db.link_string(dna_w.substr(dna_w.size() - PointerListLen, PointerListLen), 2 * ((dna_ref_b - 1- WindowListLen / 2 ) / WindowListLen) + 1);
                    }
                    else
                    {
                        db.link_string(dna_w.substr(dna_w.size() - PointerListLen, PointerListLen), 2 * ((dna_ref_b - 1 - WindowListLen / 2 ) / WindowListLen) + 1);
                        db.link_string(dna_w.substr(dna_w.size() - PointerListLen, PointerListLen), ((dna_ref_b - 1) / WindowListLen) * 2);
                    }
                }
                else 
                    db.link_string(dna_w.substr(dna_w.size() - PointerListLen, PointerListLen), 2 * ((dna_ref_b - 1 - WindowListLen / 2 ) / WindowListLen) + 1);
            }
            else if ((dna_ref_b - 1) % WindowListLen + 1>= PointerListLen)
            {
                    db.link_string(dna_w.substr(dna_w.size() - PointerListLen, PointerListLen), ((dna_ref_b - 1) / WindowListLen) * 2);
            }
            ++index_s;
            ++index_w;
            index_w %= tmp;
        }
    }



    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);

    return 0;
    
}

