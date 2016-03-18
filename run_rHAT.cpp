#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "RHT.h"
using namespace std;


typedef struct window {
    uint32_t index_of_W;
    uint32_t count;
    size_t index_of_L;
}WINDOW;


int main(int argc, char** argv) 
{
	int c;
    opterr = 0;
    while ((c = getopt (argc, argv, "p:w:")) != -1) 
    {
        switch (c)
        {
            case 'p':
                if (optarg)
                    PointerListLen = atol(optarg);
                break;
            case 'w':
                if (optarg)
                    WindowListLen = atol(optarg);
                break;
            default:
                abort ();
        }
    }
    std::cout << PointerListLen << " " << WindowListLen << "\n";

    ifstream inf, inRead, inRHT;
    inf.open("GCF_000005845.2_ASM584v2_genomic.fna");
    inRHT.open("out_RHT");
    inRead.open("dna_read");

    ofstream out;
    out.open("out");


    // ofstream outRHT;
    // outRHT.open("out_R");

    string dna_read;
    size_t index_r, num_w;
    vector<uint32_t> L_2;

    inRead >> dna_read;

    dna_bitset db(PointerListLen);
    db.read_hash_in(inRHT);
    // db.write_hash_out(outRHT);


    //dna_read is the read length = L/2
    index_r = 0;
    while (index_r < dna_read.size()) 
    {
        ++index_r;
        if (index_r >= PointerListLen)
        {
            L_2.push_back(to_bit(dna_read.substr(index_r - PointerListLen, index_r)));
            num_w += db[L_2.back()].size();
        }
    }

    // for (size_t i=0; i<L_2.size(); ++i)
    // {
    //    out<<L_2[i]<< " " << to_string(L_2[i]) << "\n";
    // }




    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
	return 0;
}