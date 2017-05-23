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


// #define PRINTLOG
#define READINFO

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
    size_t _i, t;
    bool operator<(const wq& w) const
    {
        return _i > w._i;
    }
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
    uint32_t hit_times[full];
    uint32_t tmp = 0, tar;

    // const char *c_dna_f = dna_f.c_str();
    const char *c_dna_f = ref_s;

    const size_t max_read = 500, max_window = 4096, min_window = 512, max_ref=1UL<<25;
    //#2
    size_t WC[max_ref / min_window];
    memset(WC, 0, sizeof(WC));
    size_t tmp_wc, max_wc=0;
    //#
    RHT rrht(max_read);
    bool vis[max_read + max_window];
    // size_t vis_m[30000];

    DAG dag(WindowListLen + 3);
    Cigar cg, thecg;
    cg.length = 0;
    cg.seq =  (uint32_t *)malloc(max_read * sizeof(uint32_t));
    thecg.length = 0;
    thecg.seq =  (uint32_t *)malloc(max_read * sizeof(uint32_t));
    
    char *read_s;
    size_t read_len = 0;
    gzFile fp2 = gzopen(r_n, "r");
    if (fp2)
    {
        kseq_t *seq2 = kseq_init(fp2);
        while ((l = kseq_read(seq2)) >= 0) 
        {
            #ifdef READINFO
            int max_pattern=max_read;
            char *seq_n = (char *)malloc(sizeof(char) * (100 + 1));
            char *new_p = (char *)malloc(sizeof(char) * (max_pattern + 1));
            strcpy(seq_n, seq2->name.s);
            char *seq_c = (char *)malloc(sizeof(char) * (max_pattern * 2 + 100));
            strcpy(seq_c, seq2->comment.s);
            char *p = strtok(seq_c, " =");
            char *p1, *p2;
            char *config, *chrom, *hap_ref, *edit_s;
            size_t begin, end;

            size_t extend_time = 0;

            while (p != NULL)
            {
                p1 = p;
                p2 = strtok(NULL, " =");
                p = strtok(NULL, " =");
                // printf("%s %s\n", p1, p2);
                if (strcmp(p1, "chrom") == 0)
                {
                    chrom = p2;
                }
                else if (strcmp(p1, "orig_begin") == 0)
                {
                    begin = atoi(p2);
                }
                else if (strcmp(p1, "orig_end") == 0)
                {
                    end = atoi(p2);
                    end += 20;
                }
                else if (strcmp(p1, "haplotype_infix") == 0)
                {
                    hap_ref = p2;
                }
                else if (strcmp(p1, "edit_string") == 0)
                {
                    edit_s = p2;
                }
            }
            // outt << seq_n << " " << chrom << " " << begin << " " << end << " " << edit_s << " ";
            outt << seq_n << " " << begin << " " << edit_s << " ";

            #endif
            to_upper(seq2->seq.s);
            // // ref = seq1->seq.s;
            read_len = (seq2->seq.l);
            // const char *c_read = read.c_str();
            char *c_read = seq2->seq.s;
            bool reverse = false;
            // size_t read_len = read.size();
            
            
            double score = -INF, score_enough=read_len*scy*0.95;
            double score_too_low=(-0.8)*read_len*scy;
            size_t windown_hit_too_low=20, hit_too_many=1000;
            string thess(read_len+WindowListLen, 0);
            thess.clear();

            


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
                    if (reverse) outt << "reversed ";
                    outt << "\nread seq:\n" << c_read << endl;
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
                
                uint32_t z = full, mx = 0;
                z = min((int)z, (int)Hwin_cnt.size());
                #ifdef PRINTLOG
                    outt << "match block:\ntimes    bolck    offset\n";
                #endif
                while(!Hwin_cnt.empty())
                {                
                    #ifdef PRINTLOG
                        //print the hit window and their hit times
                        // outt << Hwin_cnt.top().cnt << " " << Hwin_cnt.top().index_of_W << endl;
                        outt << Hwin_cnt.top().cnt << " " << Hwin_cnt.top().index_of_W << " " << Hwin_cnt.top().index_of_W * (WindowListLen / 2) << endl;
                    #endif
                    hit_w[--z] = Hwin_cnt.top().index_of_W;
                    hit_times[z] = Hwin_cnt.top().cnt;
                    if (z==0) 
                    {
                        mx = Hwin_cnt.top().cnt;
                    }
                    Hwin_cnt.pop();
                }
                
                if (mx < windown_hit_too_low)
                {
                    #ifdef PRINTLOG
                        //print the hit window and their hit times
                        // outt << Hwin_cnt.top().cnt << " " << Hwin_cnt.top().index_of_W << endl;
                        outt << mx << " window hit too low\n";
                    #endif
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
                        // outt << to_string(tar) << " " << i << endl;
                    }
                }
                rrht.sort_pw();
                #ifdef PRINTLOG
                    outt << "created read hash table\n";
                #endif
                // rrht.print_pw(outt);

                //-----------------------------------step 3 for each window hit do alignment
                uint32_t last_z=hit_w[0]+5;
                for (size_t z=0; z<full; ++z)
                {
                    //-----------------------------------struct window and read main body
                    // outt  << last_z << " " << hit_w[z] << endl;
                    
                    if (abs(1.0*last_z-hit_w[z]) <= 1) break;
                    
                    last_z = hit_w[z];
                    ++extend_time;
                    size_t window_up, window_down;
                    size_t ref_block_len = WindowListLen / 4;
                    //just target ref block to width at least WindowListLen + ref_block_len
                    window_up = ( hit_w[z] * (WindowListLen / 2) >= ((ref_block_len / 2) + len_up_stream) ) ? ( hit_w[z] * (WindowListLen / 2) - ((ref_block_len / 2) + len_up_stream) ) : (0);
                    window_down = ( hit_w[z] * (WindowListLen / 2) + WindowListLen + ((ref_block_len / 2) + len_down_stream) <= ref_len ) ? ( hit_w[z] * (WindowListLen / 2) + WindowListLen + ((ref_block_len / 2) + len_down_stream) ) : (ref_len);
                    #ifdef PRINTLOG
                        outt << "check block: " << hit_w[z] << " ";
                        // outt  << last_z << endl;
                        outt  << "offset from " << window_up << " to " << window_down << endl;
                        // outt << dna_f.substr(window_up, window_down - window_up) << endl;

                        // printf("%lu %lu\n", (unsigned long)window_up, (unsigned long)window_down);
                        outt << "target ref string:\n";
                        for (size_t i = window_up; i < window_down; ++i)
                        {
                            outt << ref_s[i];
                        }
                        outt<< "\n";
                    #endif
                    // struct graph node include the first and last node
                    dag.clear();
                    //the first node
                    dag.add_node(window_up + PointerListLen - 1, PointerListLen - 1, 0);
                    
                    #ifdef PRINTLOG
                        outt << "add  node begin\n";
                        // cout_string_in_char(c_read, (dag.wd)[(dag.node_i-1)].index_of_R - PointerListLen, (dag.wd)[(dag.node_i-1)].len, outt);
                        // outt << endl; 
                    #endif

                    tmp = 0;
                    // memset(vis, 0, sizeof(vis));
                    // memset(vis_m, 0, sizeof(vis_m));
                    priority_queue<WQ> Q;
                    WQ qt;
                    size_t s_ref_len = window_down - window_up;
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
                            // continue;
                            size_t pw, n, ttmp;
                            uint32_t pr;
                            // printf("%lu\n", (unsigned long)d);
                            // outt << "|";

                            if (d)
                            {   
                                while ((rrht.PW[d-1] >> 32) == tar)
                                {
                                    // outt << "-";
                                    // outt << ">" << (uint32_t)rrht.PW[d-1] << " " << to_string((rrht.PW[d-1] >> 32)) << endl;
                                    // outt << (i - window_up + read_len - (uint32_t)rrht.PW[d-1]) << endl;

                                    if (!vis[(i - window_up + read_len - (uint32_t)rrht.PW[d-1])])
                                    {
                                    // outt << ">" << (uint32_t)rrht.PW[d-1] << " " << to_string((rrht.PW[d-1] >> 32)) << endl;
                                        // outt << "asd" << endl;
                                        pw = i;
                                        pr = rrht.PW[d-1];
                                        n = dag.add_node(pw, pr, PointerListLen);

                                        ttmp = (i - window_up + read_len - pr);

                                        //if overlap
                                        vis[ttmp]=true;
                                        // outt << "B: " << pr << " "<<(i - window_up + read_len - pr) << " ";
                                        ++pw;
                                        ++pr;
                                        while(c_dna_f[pw]==c_read[pr] && dag.check_node(pw, pr, n))
                                        {
                                            // vis[(pw - window_up + read_len - pr)]=true;
                                            ++pw;
                                            ++pr;
                                            // outt << (unsigned long)pw - window_up << "/" << pr << " ";
                                        }
                                        // for (size_t t = 1; t<PointerListLen-1; ++t) vis[pr + t] = true;
                                        qt.t = ttmp;
                                        qt._i = pw - 1;
                                        // vis_m[pw-window_up] = 1 + ttmp;
                                        Q.push(qt);
                                        // outt << wd[node_i].index_of_W << " " << wd[node_i].index_of_R << " " << wd[node_i].len << " ";
                                        #ifdef PRINTLOG
                                        // outt << "+ " << pr << " "<< qt._i - window_up << " ";
                                        // cout_string_in_char(c_read, (dag.wd)[(dag.node_i-1)].index_of_R - PointerListLen + 1, (dag.wd)[(dag.node_i-1)].len, outt);
                                        // outt << endl; 
                                        #endif
                                    }
                                    ++d;
                                }
                            }

                            
                            while (!Q.empty() && Q.top()._i == i)
                            {
                                // outt << "F " << c_dna_f[i] << " " << i - window_up << "> \n" ;

                                qt = Q.top();
                                Q.pop();
                                vis[qt.t] = false;
                            }

                            // if (vis_m[i-window_up+1])
                            // {
                            //     vis[vis_m[i-window_up+1] - 1] = 0;
                            //     // outt << "F " << vis_m[i-window_up+1] - 1 << " " << i - window_up << "> \n" ;
                            //     // memset(vis, 0, sizeof(vis));
                            //     qt = Q.front();
                            //     Q.pop();
                            //     d = qt.d;
                            //     while ((rrht.PW[d-1] >> 32) == qt.tar)
                            //     {

                            //         pw = qt.i;
                            //         pr = rrht.PW[d-1];
                            //         //if overlap
                            //         vis[pr]=false;
                            //         ++pw;
                            //         ++pr;
                            //         while(vis[pr] && c_dna_f[pw]==c_read[pr] && dag.check_node(pw, pr, qt.n, 1))
                            //         {
                            //             vis[pr]=false;
                            //             ++pw;
                            //             ++pr;
                            //         }
                            //         // for (size_t t = 1; t<PointerListLen-2; ++t) vis[pr + t] = false;
                            //         ++d;
                            //         // outt << wd[node_i].index_of_W << " " << wd[node_i].index_of_R << " " << wd[node_i].len << " " << read.substr(wd[node_i].index_of_R, wd[node_i].len) << endl; 
                            //         // outt << "> ";
                            //         // cout_string_in_char(c_read, (dag.wd)[(dag.node_i-1)].index_of_R - PointerListLen + 1, (dag.wd)[(dag.node_i-1)].len, outt);
                            //         // outt << endl; 
                            //     // }

                            // }
                        }
                    }
                    dag.add_node(window_down + PointerListLen - 1, read_len + PointerListLen - 1, 0);
                    #ifdef PRINTLOG
                         outt << "add  node end\n";
                    #endif

                    // dag.print_log(outt);

                    //linking nodes edge
                    //node matrix 
                    // dag.create_matrix();
                    //do use node matrix to run dp
                     #ifdef PRINTLOG
                         outt << "print the node path:\noffset   node_index\n";
                    #endif
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

                    #ifdef PRINTLOG
                        outt << "cigar: ";
                        for (size_t i = 0; i<(cg.length); ++i)
                        {
                            outt << cigar_int_to_len(cg.seq[i]) << cigar_int_to_op(cg.seq[i]);
                        }
                        outt << endl;
                    #endif


                    #ifdef PRINTLOG
                        outt << "score: ";
                        outt << t << endl;
                        outt<< "~~~~~~~~~~~~~~~~~\n";
                    #endif


                    if ((t) > score)
                    {
                        score = t;

                        // strcpy(thesss, sss);
                        // thecg.length = cg.length;
                        // // memset(thecg.seq, 0, sizeof(thecg.seq));
                        // memcpy(thecg.seq, cg.seq, sizeof(cg.seq));
                        copy_cigar(&thecg, &cg);
                        
                        // thess = ss;
                    }
                    else if (z > 0 && hit_times[z] > hit_times[z-1])
                    {
                        #ifdef PRINTLOG
                            outt << hit_times[z] << " " << hit_times[z-1] << "\n";
                            outt << "break for score no increse\n";
                        #endif
                        break;
                    }

                    if (score < score_too_low)
                    {
                        #ifdef PRINTLOG
                            outt << "break for best score too low\n";
                        #endif
                        break;
                    }
                    // dag.do_alignment(c_dna_f, window_up, window_down - PointerListLen + 1, c_read, out);

                    if (t <= score_too_low) 
                    {
                        #ifdef PRINTLOG
                            outt << "break for score too low\n";
                        #endif
                        break;
                    }
                    if (score >= score_enough) 
                    {
                        #ifdef PRINTLOG
                            outt << "break for score_enough\n";
                        #endif
                        break;
                    }   
                }

                if (score >= score_enough)
                {
                    #ifdef PRINTLOG
                        outt << "no reverse for score_enough\n";
                    #endif
                    break;
                }
                
                // break;
                if (reverse) break;
                reverse = true;
            }


            outt << extend_time << " ";
            for (size_t i = 0; i<(thecg.length); ++i)
            {
                outt << cigar_int_to_len(thecg.seq[i]) << cigar_int_to_op(thecg.seq[i]);
            }
            outt << endl;

            #ifdef PRINTLOG
                outt << "\nfinal cigar: ";
                for (size_t i = 0; i<(thecg.length); ++i)
                {
                    outt << cigar_int_to_len(thecg.seq[i]) << cigar_int_to_op(thecg.seq[i]);
                }
                outt << endl;
                outt << "final score: ";
                outt << score << endl;
                outt<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

            #endif


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