#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <sstream>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include <zlib.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "kseq.h"
#include "DAG.h"
using namespace std;

KSEQ_INIT(gzFile, gzread)


#define PRINTLOG

void to_upper(char *s)
{
  size_t l = strlen(s);
  for (size_t i=0; i<l; ++i)
    s[i] = toupper(s[i]);
}


typedef struct window_cnt {
    uint32_t index_of_W;
    size_t cnt;
    bool operator<(const window_cnt& w) const
    {
        return cnt > w.cnt;
    }
}WINDOW_CNT;

typedef struct wq{
    size_t _i, i, n;
    uint32_t d, tar;
}WQ;

void print_usage()
{
    printf("\tuse out_RHT to run rHAT\n");
    printf("\trun_rHAT [option]\n");
    printf("\t-f \tdna file           [FILE]\n");
    printf("\t-r \tread file          [FILE]\n");
    printf("\t-o \toutput file        [out]\n");
    printf("\t-w \tset WindowListLen  [%lu]\n", (unsigned long)WindowListLen);
    printf("\t-p \tset PointerListLen [%lu]\n", (unsigned long)PointerListLen);
    printf("\t-h \thelp\n");
}

int main(int argc, char** argv) 
{ 
    int c;
    opterr = 0;
    char d_n[100]="";
    char r_n[100]="";
    char o_n[100]="out";
    while ((c = getopt (argc, argv, "p:w:f:r:o:h")) != -1) 
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
            case 'r':
                if (optarg)
                    strcpy(r_n, optarg);
                break;
            case 'o':
                if (optarg)
                    strcpy(o_n, optarg);  
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
              << "| dna file: " << string(d_n) << " " 
              << "| read file: " << string(r_n) << endl;


    //read ref
    char *ref_s = (char *)malloc(10 * sizeof(char));;
    size_t ref_len = 0;
    int l = 0;
    gzFile fp1 = gzopen(d_n, "r");
    if (fp1)
    {
        kseq_t *seq1 = kseq_init(fp1);
        while ((l = kseq_read(seq1)) >= 0) {
            ref_s = (char *)realloc(ref_s, (ref_len + (seq1->seq.l+1)) * sizeof(char*));
            // // ref_num = (int8_t *)malloc(sizeof(int8_t) * seq1->seq.l);
            to_upper(seq1->seq.s);
            strncpy(ref_s + ref_len, seq1->seq.s, seq1->seq.l);
            // // ref = seq1->seq.s;
            ref_len += (seq1->seq.l);
            ref_s[ref_len] = '\0';
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

    printf("load ref succeed, total len: %lu\n", (unsigned long)ref_len);
    // for (size_t i = 0; i<200; ++i)
    // {
    //     printf("%c", ref_s[i]);
    // }
    // printf("\n");
    // for (size_t i = ref_len-200; i<ref_len; ++i)
    // {
    //     printf("%c", ref_s[i]);
    // }
    // printf("\n");


    // ifstream inRef, inRead;
    FILE *out;
    FILE *inRHT;
    // inRef.open(d_n);
    // inRead.open(r_n);
    out = fopen(o_n, "wb");

    strcat(d_n, ".out_RHT");
    printf("read hash from %s\n", d_n);
    inRHT = fopen(d_n, "rb");

    //-----------------------------------get dna_ref string
    // string dna_name, dna_s, dna_f;
    // getline(inRef, dna_name);
    // while (inRef >> dna_s) dna_f.append(dna_s);

    //-----------------------------------get hash table
    // RHT rht(dna_f.size());
    RHT rht(ref_len);
    rht.read_hash(inRHT);
    cout << "succeed load RHT\n";
    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    
    ofstream outt;
    outt.open("check.out");
    // rht.write_hash_test(outt);

    //-----------------------------------get read
    // string read_m, read, read_x;
    //init mask
    MASK = MASK >> (32 - PointerListLen * 2);
    // //the hit window number that will go to alignment of each read
    const uint32_t full = 3;
    uint32_t hit_w[full];
    uint32_t tmp = 0, tar;

    // const char *c_dna_f = dna_f.c_str();
    const char *c_dna_f = ref_s;

    char sss[10000];
    char thesss[10000];


    //#2
    size_t WC[10000];
    memset(WC, 0, sizeof(WC));
    size_t tmp_wc, max_wc=0;
    //#
    RHT rrht(30000);
    bool vis[30000];


    char *read_s;
    size_t read_len = 0;
    gzFile fp2 = gzopen(r_n, "r");
    if (fp2)
    {
        kseq_t *seq2 = kseq_init(fp2);
        while ((l = kseq_read(seq2)) >= 0) 
        {
            to_upper(seq2->seq.s);
            // // ref = seq1->seq.s;
            read_len = (seq2->seq.l);
            // const char *c_read = read.c_str();
            char *c_read = seq2->seq.s;
            bool reverse = false;
            // size_t read_len = read.size();
            
            
            double score = -INF, score_enough=100;
            double score_too_low=-100;
            size_t windown_hit_too_low=10, hit_too_many=500;
            string thess(read_len+WindowListLen, 0);
            thess.clear();

            Cigar cg;
            cg.length = 0;
            cg.seq =  (uint32_t *)malloc(100 * sizeof(uint32_t));


            while (true)
            {
                //-----------------------------------step 1 find hit windows 
                //Hwin is the first heap to set window priority_queue

                //reverse read
                if (reverse)
                {
                    for (size_t i = 0; i<read_len/2; ++i)
                    {   
                        char tmp = c_read[i];
                        c_read[i] = c_read[read_len - 1 - i];
                        c_read[read_len - 1 - i] = tmp;
                    }

                    for (size_t i = 0; i<read_len; ++i)
                    {
                        switch(c_read[i])
                        {
                            case 'A':
                                c_read[i] = 'T';
                                break;
                            case 'T':
                                c_read[i] = 'A';
                                break;
                            case 'G':
                                c_read[i] = 'C';
                                break;
                            case 'C':
                                c_read[i] = 'G';
                                break;    
                        }
                    }
                }
                #ifdef PRINTLOG
                    // printf("%s %d %s\n", seq2->name.s, reverse, c_read);
                    outt << seq2->name.s << " " << reverse << " " << c_read << endl;
                #endif

        


                //#1 split the center read
                size_t len_up_stream = 0, len_down_stream = 0;
                size_t theLen = WindowListLen / 2;
                //from len_up_stream to theLen is the center string of read
                if (read_len * 2 > WindowListLen)
                {
                    len_up_stream = read_len / 2 - theLen / 2;
                    len_down_stream = read_len - len_up_stream - theLen;
                }
                else
                {   
                    len_up_stream = 0;
                    len_down_stream = 0;
                    theLen = read_len;
                }
                // outt << read.substr(len_up_stream, theLen) << endl;

                //print log
                // #ifdef PRINTLOG
                // printf("the center read is:\n");
                // for (size_t i=0; i<theLen; ++i)
                // {
                //     printf("%c", c_read[i+len_up_stream]);
                // }
                // printf("\n");
                // #endif
                //#


                //#2 get the most hitted window
                tmp = 0;
                max_wc = 0;
                for (size_t i=len_up_stream, k; i < len_up_stream + theLen; ++i)
                {
                    tmp = (tmp << 2) | get_c[c_read[i]];
                    tar = tmp & MASK;
                    if (i - len_up_stream >= PointerListLen - 1)
                    {
                        k = rht.P[tar];
                        while (k && rht.P[tar+1]>k)
                        {
                            tmp_wc = rht.W[k - 1];
                            max_wc = max(max_wc, tmp_wc);
                            ++WC[tmp_wc];
                            ++k;
                        }
                    }
                }

                //Hwin_cnt is the second heap to count window hit number
                priority_queue<WINDOW_CNT> Hwin_cnt;
                WINDOW_CNT wc;
                for (size_t i=0, l=0; i<=max_wc; ++i)
                {
                    if (WC[i] > l)
                    {
                        wc.index_of_W = i;
                        wc.cnt = WC[i];
                        if (wc.cnt < hit_too_many) 
                        {
                            Hwin_cnt.push(wc);
                            if (Hwin_cnt.size() > full) Hwin_cnt.pop();
                            l = Hwin_cnt.top().cnt;
                        }
                    }
                    WC[i]=0;
                }
                
                uint32_t z = full, mx;
                while(!Hwin_cnt.empty())
                {                
                    #ifdef PRINTLOG
                        //print the hit window and their hit times
                        // outt << Hwin_cnt.top().cnt << " " << Hwin_cnt.top().index_of_W << endl;
                        outt << Hwin_cnt.top().cnt << " " << Hwin_cnt.top().index_of_W << endl;
                    #endif
                    hit_w[--z] = Hwin_cnt.top().index_of_W;
                    if (z==0) 
                    {
                        mx = Hwin_cnt.top().cnt;
                    }
                    Hwin_cnt.pop();
                }
                
                if (mx < windown_hit_too_low)
                {
                    if (reverse) break;
                    reverse = true;
                    continue;
                }
                //#


                //-----------------------------------step 2 create read hash table
                tmp = 0;
                rrht.clear();
                for (size_t i=0; i < read_len; ++i)
                {
                    tmp = (tmp << 2) | get_c[c_read[i]];
                    tar = tmp & MASK;
                    if (i >= PointerListLen - 1)
                    {
                        rrht.link_string(tar, i);
                    }
                }
                rrht.sort_pw();

                //-----------------------------------step 3 for each window hit do alignment
                uint32_t last_z=hit_w[0]+5;
                for (size_t z=0; z<full; ++z)
                {
                    //-----------------------------------struct window and read main body
                    outt  << last_z << " " << hit_w[z] << endl;
                    
                    if (abs(1.0*last_z-hit_w[z]) <= 1) break;
                    // if (abs(hit_w[z]-last_z) <= 1) break;
                    last_z = hit_w[z];
                    size_t window_up, window_down;
                    window_up = ( hit_w[z] * (WindowListLen / 2) >= len_up_stream ) ? ( hit_w[z] * (WindowListLen / 2) - len_up_stream ) : (0);
                    window_down = ( hit_w[z] * (WindowListLen / 2) + WindowListLen + len_down_stream <= ref_len ) ? ( hit_w[z] * (WindowListLen / 2) + WindowListLen + len_down_stream ) : (ref_len);
                    #ifdef PRINTLOG
                    outt  << last_z << endl;
                    outt  << window_up << " " << window_down << endl;
                    // outt << dna_f.substr(window_up, window_down - window_up) << endl;

                    // printf("%lu %lu\n", (unsigned long)window_up, (unsigned long)window_down);
                    for (size_t i = window_up; i < window_down; ++i)
                    {
                        outt << ref_s[i];
                    }
                    outt<< "\n";
                    #endif
                    // struct graph node include the first and last node
                    DAG dag(window_down - window_up + 3);
                    dag.clear();
                    //the first node
                    dag.add_node(window_up + PointerListLen - 1, PointerListLen - 1, 0);
                    
                    outt << "> ";
                    cout_string_in_char(c_read, (dag.wd)[(dag.node_i-1)].index_of_R - PointerListLen, (dag.wd)[(dag.node_i-1)].len, outt);
                    outt << endl; 

                    tmp = 0;
                    memset(vis, 0, sizeof(vis));
                    queue<WQ> Q;
                    WQ qt;
                    for (size_t i=window_up; i < window_down; ++i)
                    {
                        tmp = (tmp << 2) | get_c[c_dna_f[i]];
                        tar = tmp & MASK;
                        // outt << "**" << endl;

                        if (i - window_up >= PointerListLen - 1)
                        {
                            //if find window & read match l-mer
                            uint32_t d = 0;
                            d = rrht.search(tar);
                            size_t pw, n;
                            uint32_t pr;
                            // printf("%lu\n", (unsigned long)d);
                            if (d)
                            {
                                // outt << to_string(tar) << " " << d << endl;
                                while ((rrht.PW[d-1] >> 32) == tar && !vis[(uint32_t)rrht.PW[d-1]])
                                {
                                    // outt << "asd" << endl;
                                    pw = i;
                                    pr = rrht.PW[d-1];
                                    n = dag.add_node(pw, pr, PointerListLen);

                                    qt.i = i;
                                    qt.d = d;
                                    qt.tar = tar;
                                    qt.n = n;
                                    // outt << "+ " << pr << " "<< i - window_up << " ";

                                    //if overlap
                                    vis[pr]=true;
                                    ++pw;
                                    ++pr;
                                    while(!vis[pr] && c_dna_f[pw]==c_read[pr] && dag.check_node(pw, pr, n))
                                    {
                                        // vis[pr]=true;
                                        ++pw;
                                        ++pr;
                                        // outt << (unsigned long)pw - window_up << "/" << pr << " ";
                                    }
                                    // for (size_t t = 1; t<PointerListLen-1; ++t) vis[pr + t] = true;
                                    qt._i = pw - 1;
                                    Q.push(qt);
                                    ++d;
                                    // outt << wd[node_i].index_of_W << " " << wd[node_i].index_of_R << " " << wd[node_i].len << " ";
                                    outt << "+ " << pr << " "<< qt._i - window_up << " ";
                                    cout_string_in_char(c_read, (dag.wd)[(dag.node_i-1)].index_of_R - PointerListLen + 1, (dag.wd)[(dag.node_i-1)].len, outt);
                                    outt << endl; 
                                }
                            }

                            


                            if (i==Q.front()._i)
                            {
                                // outt << i - window_up << "> \n" ;
                                memset(vis, 0, sizeof(vis));
                                // qt = Q.front();
                                // Q.pop();
                                // d = qt.d;
                                // while ((rrht.PW[d-1] >> 32) == qt.tar)
                                // {

                                //     pw = qt.i;
                                //     pr = rrht.PW[d-1];
                                //     //if overlap
                                //     vis[pr]=false;
                                //     ++pw;
                                //     ++pr;
                                //     while(vis[pr] && c_dna_f[pw]==c_read[pr] && dag.check_node(pw, pr, qt.n, 1))
                                //     {
                                //         vis[pr]=false;
                                //         ++pw;
                                //         ++pr;
                                //     }
                                //     // for (size_t t = 1; t<PointerListLen-2; ++t) vis[pr + t] = false;
                                //     ++d;
                                //     // outt << wd[node_i].index_of_W << " " << wd[node_i].index_of_R << " " << wd[node_i].len << " " << read.substr(wd[node_i].index_of_R, wd[node_i].len) << endl; 
                                //     // outt << "> ";
                                //     // cout_string_in_char(c_read, (dag.wd)[(dag.node_i-1)].index_of_R - PointerListLen + 1, (dag.wd)[(dag.node_i-1)].len, outt);
                                //     // outt << endl; 
                                // }

                            }
                        }
                    }
                    dag.add_node(window_down + PointerListLen - 1, read_len + PointerListLen - 1, 0);

                    // dag.print_log(outt);

                    //linking nodes edge
                    //node matrix 
                    dag.create_matrix();
                    //do use node matrix to run dp
                    dag.find_path();

                    // -----------------------------------use global & semiglobal alignment to struct alignment and get sroce
                    #ifdef PRINTLOG
                        // outt << read_m << endl;
                        dag.print_log(const_cast<char*>(c_dna_f), const_cast<char*>(c_read), outt);
                    #endif

                    double t;
                    // string ss;
                    // // dag.do_alignment(dna_f, read, outt);
                    // t = dag.do_alignment(dna_f, window_up, window_down, read, ss, outt);

                    // strcpy(sss, "");
                    // t = dag.do_alignment(const_cast<char*>(c_dna_f), window_up, window_down, const_cast<char*>(c_read), read_len, sss, out);
                    
                    cg.length = 0;
                    t = dag.do_alignment(const_cast<char*>(c_dna_f), window_up, window_down, const_cast<char*>(c_read), read_len, &cg, out);
                    // t = dag.do_alignment(const_cast<char*>(c_dna_f), window_up, window_down, const_cast<char*>(c_read), read_len, ss, out);

                    for (size_t i = 0; i<(cg.length); ++i)
                    {
                        outt << cigar_int_to_len(cg.seq[i]) << cigar_int_to_op(cg.seq[i]);
                    }
                    outt << endl;

                    outt<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";


                    if ((t) > score)
                    {
                        score = t;

                        memset(thesss, 0, sizeof(thesss));
                        strcpy(thesss, sss);
                        
                        // thess = ss;
                    }
                    else break;

                    if (score < score_too_low) break;
                    // dag.do_alignment(c_dna_f, window_up, window_down - PointerListLen + 1, c_read, out);

                    if (t <= score_too_low) break;
                    if (score >= score_enough) break;   
                }

                if (score >= score_enough) break;
                break;
                if (reverse) break;
                reverse = true;
            }



            //-----------------------------------out put the read best alignment
            // fprintf(out, "%s\n", read_m.c_str());
            // fprintf(out, "%s\n", seq2->name.s);
            // fprintf(out, "%d\n", reverse);
            // fprintf(out, "%.0lf\n", score);
            // fprintf(out, "%s\n", thesss);
            
            // outt << score << endl;
            // outt << thesss << endl;
            // outt << endl;
           
            // printf("%.0lf\n", score);
            // printf("%s\n", thess.c_str());


        }
        kseq_destroy(seq2);
    }
    else
    {
        fprintf(stderr, "read file %s error\n", r_n);
    }
    gzclose(fp2);
    
    

    

    // inRef.close();
    // inRead.close();
    outt.close();
    fclose(out);
    fclose(inRHT);


        


    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}