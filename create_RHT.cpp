#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "RHT.h"
using namespace std;

// #define PRINTLOG

int main(int argc, char** argv) 
{

    int c;
    opterr = 0;
    char d_n[100]="E.coli.fa";
    //set argument,p for PointerListLen,w for WindowListLen,f for dna_read file(if the -pwf is provided, argument is required)
    while ((c = getopt (argc, argv, "p:w:d:h")) != -1) 
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
            case 'd':
                if (optarg)
                    strcpy(d_n, optarg);
                break;
            case 'h':
                printf("\tuse dna file to create a RHT table to file out_RHT\n");
                printf("\tcreate_RHT [option]\n");
                printf("\t-d \tdna file           [%s]\n", d_n);
                printf("\t-w \tset WindowListLen  [%lu]\n", (unsigned long)WindowListLen);
                printf("\t-p \tset PointerListLen [%lu]\n", (unsigned long)PointerListLen);
                printf("\t-h \thelp\n");
                return 0;
                break;
            default:
                abort ();
        }
    }

    std::cout << "PointerListLen: " << PointerListLen << " " 
              << "| WindowListLen: " << WindowListLen << " "
              << "| dna file: " << string(d_n) << endl;


    ifstream inf;
    inf.open(d_n);

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
