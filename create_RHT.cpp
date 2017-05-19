#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <queue>
#include <unistd.h>
#include <algorithm>
#include <sstream>
#include <zlib.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "kseq.h"
#include "RHT.h"
using namespace std;

KSEQ_INIT(gzFile, gzread)

// #define PRINTLOG


void to_upper(char *s)
{
  size_t l = strlen(s);
  for (size_t i=0; i<l; ++i)
    s[i] = toupper(s[i]);
}

void print_usage()
{
    printf("\tuse dna file to create a RHT table to file out_RHT\n");
    printf("\tcreate_RHT [option]\n");
    printf("\t-f \tdna file           [FILE]\n");
    printf("\t-w \tset WindowListLen  [%lu]\n", (unsigned long)WindowListLen);
    printf("\t-p \tset PointerListLen [%lu]\n", (unsigned long)PointerListLen);
    printf("\t-h \thelp\n");
}

int main(int argc, char** argv) 
{

    int c;
    opterr = 0;
    char d_n[100]="";
    //set argument,p for PointerListLen,w for WindowListLen,f for dna_read file(if the -pwf is provided, argument is required)
    while ((c = getopt (argc, argv, "p:w:f:h")) != -1) 
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
            case 'f':
                if (optarg)
                    strcpy(d_n, optarg);
                break;
            case 'h':
                print_usage();
                return 0;
                break;
            default:
                abort ();
        }
    }

    if (d_n[0] == '\0')
    {
        print_usage();
        return 1;
    }

    std::cout << "PointerListLen: " << PointerListLen << " " 
              << "| WindowListLen: " << WindowListLen << " "
              << "| dna file: " << string(d_n) << endl;

    //read ref
    char *res_s = (char *)malloc(10 * sizeof *res_s);;
    size_t ref_len = 0;
    int l = 0;
    gzFile fp1 = gzopen(d_n, "r");
    if (fp1)
    {
        kseq_t *seq1 = kseq_init(fp1);
        while ((l = kseq_read(seq1)) >= 0) {
            res_s = (char *)realloc(res_s, (ref_len + (seq1->seq.l+1)) * sizeof(char*));
            // // ref_num = (int8_t *)malloc(sizeof(int8_t) * seq1->seq.l);
            to_upper(seq1->seq.s);
            strncpy(res_s + ref_len, seq1->seq.s, seq1->seq.l);
            // // ref = seq1->seq.s;
            ref_len += (seq1->seq.l);
            res_s[ref_len] = '\0';
            // ref_num = char_create_num(ref, seq1->seq.l, nt_table);
            // for (m = 0; m < ref_len; ++m) ref_num[m] = nt_table[(int)ref[m]];
        }
        // printf("%s\n", res_s);
        kseq_destroy(seq1);
    }
    else
    {
        fprintf(stderr, "read ref file %s error\n", d_n);
    }
    gzclose(fp1);

    // return 0;

    // ifstream inf;
    // inf.open(d_n);

    FILE *pout;
    strcat(d_n, ".out_RHT");
    printf("write hash to %s\n", d_n);
    pout = fopen(d_n, "wb+");

    #ifdef PRINTLOG
    ofstream out;
    out.open("out_R");
    #endif


    string dna_name, dna_s, dna_w;

    // getline(inf, dna_name);
    // while (inf >> dna_s) dna_w.append(dna_s);

    uint32_t tmp = 0, tar;

    //init mask
    MASK = MASK >> (32 - PointerListLen * 2);

    // std::cout << dna_w.size() << endl;
    
    //init rht table with dna_size
    // RHT db(dna_w.size());
    RHT db(ref_len);


    for (size_t i=0; i<ref_len; ++i)
    // for (size_t i=0; i<dna_w.size(); ++i)
    {   
        tmp = (tmp << 2) | get_c[res_s[i]];
        // tmp = (tmp << 2) | get_c[dna_w[i]];
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
    // inf.close();

    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);

    return 0;
    
}
