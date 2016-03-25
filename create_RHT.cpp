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

    //receive PointerListLength and windowListLength
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
    pout = fopen("out_RHT", "w");
    
    // ofstream pout;
    // pout.open("out_RHT");

    string dna_name, dna_s, dna_w;

    dna_bitset db(PointerListLen);
    // std::cout << to_string(to_bit("AACTAATTGGT")) << std::endl;

    uint64_t dna_ref_b = 0;
    size_t  index_s = 0;

    getline(inf, dna_name);
    // std::cout << dna_name  << endl;

    #ifdef PRINT_WINDOW_INDEX
        ofstream outf;
        outf.open("out_window");
    #endif

    size_t k, dna_s_size;

// SCANF_WINDOW
    while (inf >> dna_s)
    {

        dna_w.append(dna_s);
        index_s = 0;
        k = 0;
        dna_s_size = dna_s.size();
        while (index_s < dna_s_size) 
        {
            ++dna_ref_b;
            k = dna_s_size - index_s - 1;
            ++index_s;
           if (dna_ref_b > WindowListLen / 2 && ((dna_ref_b - 1 - WindowListLen / 2) % WindowListLen + 1 >= PointerListLen))
            {
                if ((dna_ref_b - 1) % WindowListLen + 1>= PointerListLen)
                {
                    if ((dna_ref_b / WindowListLen) * 2 < 2 * ((dna_ref_b - WindowListLen / 2 ) / WindowListLen) + 1)
                    {
                        db.link_string(dna_w.substr(dna_w.size() - k - PointerListLen, PointerListLen), ((dna_ref_b - 1) / WindowListLen) * 2);
                        db.link_string(dna_w.substr(dna_w.size() - k  - PointerListLen, PointerListLen), 2 * ((dna_ref_b - 1- WindowListLen / 2 ) / WindowListLen) + 1);
                    }
                    else
                    {
                        db.link_string(dna_w.substr(dna_w.size() - k  - PointerListLen, PointerListLen), 2 * ((dna_ref_b - 1 - WindowListLen / 2 ) / WindowListLen) + 1);
                        db.link_string(dna_w.substr(dna_w.size() - k  - PointerListLen, PointerListLen), ((dna_ref_b - 1) / WindowListLen) * 2);
                    }
                }
                else 
                    db.link_string(dna_w.substr(dna_w.size() - k  - PointerListLen, PointerListLen), 2 * ((dna_ref_b - 1 - WindowListLen / 2 ) / WindowListLen) + 1);
            }
            else if ((dna_ref_b - 1) % WindowListLen + 1>= PointerListLen)
            {
                    db.link_string(dna_w.substr(dna_w.size() - k  - PointerListLen, PointerListLen), ((dna_ref_b - 1) / WindowListLen) * 2);
            }

            #ifdef PRINT_WINDOW_INDEX
                if (dna_ref_b >= WindowListLen)
                    if (dna_ref_b % (WindowListLen / 2) == 0)
                        if (dna_ref_b % (WindowListLen) == 0)
                            outf << dna_w.substr(dna_w.size() - WindowListLen, WindowListLen) << " " <<  (dna_ref_b / WindowListLen) * 2 - 2<< "| \n";
                        else
                            outf << dna_w.substr(dna_w.size() - WindowListLen, WindowListLen) << " " <<  2 * ((dna_ref_b - WindowListLen / 2 ) / WindowListLen) - 1<< "| \n";
            #endif
        }

        // #ifndef PRINT_WINDOW_INDEX
        // if (dna_w.size() > 2 * WindowListLen) dna_w = dna_w.substr(dna_w.size()-PointerListLen*2, PointerListLen*2);
        // #endif

        #ifdef PRINT_WINDOW_INDEX
        if (dna_w.size() > 5 * WindowListLen) dna_w = dna_w.substr(dna_w.size()-WindowListLen, WindowListLen);
        #endif
    }

    #ifdef PRINT_WINDOW_INDEX
        if (dna_ref_b % WindowListLen != 0)
        if (dna_ref_b % WindowListLen >= WindowListLen / 2) 
        {
            outf << dna_w.substr(dna_w.size() - dna_ref_b % WindowListLen, dna_ref_b % WindowListLen) << " " <<  (dna_ref_b / WindowListLen) * 2 << "| \n";
            outf << dna_w.substr(dna_w.size() - dna_ref_b % WindowListLen + WindowListLen / 2, dna_ref_b % WindowListLen - WindowListLen / 2) << " " <<  2 * ((dna_ref_b - WindowListLen / 2 ) / WindowListLen) + 1<< "| \n";
        }
        else
        {
            outf << dna_w.substr(dna_w.size() - dna_ref_b % WindowListLen - WindowListLen / 2, dna_ref_b % WindowListLen + WindowListLen / 2) << " " <<  2 * ((dna_ref_b - WindowListLen / 2 ) / WindowListLen) + 1<< "| \n";
            outf << dna_w.substr(dna_w.size() - dna_ref_b % WindowListLen, dna_ref_b % WindowListLen) << " " <<  (dna_ref_b / WindowListLen) * 2 << "| \n";
        }

        outf.close();
    #endif


    db.write_hash_out(pout);
    
    fclose(pout);
    inf.close();


    RHT rht(dna_ref_b);
    for (size_t i = 0; i < dna_ref_b; ++i)
    {
        if (dna_ref_b > WindowListLen / 2 && ((dna_ref_b - 1 - WindowListLen / 2) % WindowListLen + 1 >= PointerListLen))
        {
            if ((dna_ref_b - 1) % WindowListLen + 1>= PointerListLen)
            {
                if ((dna_ref_b / WindowListLen) * 2 < 2 * ((dna_ref_b - WindowListLen / 2 ) / WindowListLen) + 1)
                {
                    rht.link(dna_w.substr(dna_ref_b - PointerListLen, PointerListLen), ((dna_ref_b - 1) / WindowListLen) * 2);
                    rht.link(dna_w.substr(dna_ref_b  - PointerListLen, PointerListLen), 2 * ((dna_ref_b - 1- WindowListLen / 2 ) / WindowListLen) + 1);
                }
                else
                {
                    rht.link(dna_w.substr(dna_ref_b  - PointerListLen, PointerListLen), 2 * ((dna_ref_b - 1 - WindowListLen / 2 ) / WindowListLen) + 1);
                    rht.link(dna_w.substr(dna_ref_b  - PointerListLen, PointerListLen), ((dna_ref_b - 1) / WindowListLen) * 2);
                }
            }
            else 
                rht.link(dna_w.substr(dna_ref_b  - PointerListLen, PointerListLen), 2 * ((dna_ref_b - 1 - WindowListLen / 2 ) / WindowListLen) + 1);
        }
        else if ((dna_ref_b - 1) % WindowListLen + 1>= PointerListLen)
        {
                rht.link(dna_w.substr(dna_ref_b  - PointerListLen, PointerListLen), ((dna_ref_b - 1) / WindowListLen) * 2);
        }
    }
    FILE *oo;
    oo = fopen("oo", "w+");
    rht.print_hash(oo);
    fclose(oo);

    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);

    return 0;
    
}

