#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "RHT.h"
using namespace std;

#define PRINTLOG

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


    ifstream inf;
    inf.open("E.coli.fa");

    FILE *pout;
    pout = fopen("out_RHT", "wb");

    #ifdef PRINTLOG
    ofstream out;
    out.open("out_R");
    #endif


    string dna_name, dna_s, dna_w;

    getline(inf, dna_name);
    while (inf >> dna_s) dna_w.append(dna_s);

    uint32_t tmp = 0, tar;

    //init mask
    MASK = MASK >> (32 - PointerListLen * 2);

    // std::cout << dna_w.size() << endl;
    
    //init rht table with dna_size
    RHT db(dna_w.size());


    for (size_t i=0; i<dna_w.size(); ++i)
    {   
        tmp = (tmp << 2) | get_c[dna_w[i]];
        tar = tmp & MASK;
        if (i>=PointerListLen - 1)
        {
            //if i > window NO.1
            if (i > WindowListLen / 2 && ((i - WindowListLen / 2) % WindowListLen + 1 >= PointerListLen))
            {
                if (i % WindowListLen + 1>= PointerListLen)
                {
                    if ((i /  WindowListLen) * 2 < 2 * ((i - WindowListLen / 2 ) / WindowListLen) + 1)
                    {
                        db.link_string(tar, (i / WindowListLen) * 2);
                        db.link_string(tar, 2 * ((i- WindowListLen / 2 ) / WindowListLen) + 1);
                    }
                    else
                    {
                        db.link_string(tar, 2 * ((i - WindowListLen / 2 ) / WindowListLen) + 1);
                        db.link_string(tar, (i / WindowListLen) * 2);
                    }
                }
                else 
                    db.link_string(tar,2 * ((i - WindowListLen / 2 ) / WindowListLen) + 1);
            }
            else if (i % WindowListLen + 1>= PointerListLen)
            {
                    db.link_string(tar, (i / WindowListLen) * 2);
            }

        }
    } 

    db.create_p_w();
    db.write_hash(pout);

    #ifdef PRINTLOG
    db.write_hash_test(out);
    out.close();
    #endif


    fclose(pout);
    inf.close();

    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);

    return 0;
    
}
