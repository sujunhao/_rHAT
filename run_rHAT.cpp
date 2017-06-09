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
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "kseq.h"
#include "DAG.h"
// #include "LV_deep.h"
#include "SNP_vector.h"
#include "Haffman.h"
using namespace std;

// KSEQ_INIT(gzFile, gzread)


// #define PRINTLOG
// #define READINFO


void char2char(char *s)
{
  size_t l = strlen(s);
  for (size_t i=0; i<l; ++i)
    s[i] = map_c[(s[i] =='N') ? rand()%4 : map_bit[s[i]]];
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
    printf("\t-v \tVCF file           [FILE]\n");
    printf("\t-r \tread file          [FILE]\n");
    printf("\t-o \toutput file        [out]\n");
    printf("\t-w \tset WindowListLen  [%lu]\n", (unsigned long)WindowListLen);
    printf("\t-p \tset PointerListLen [6-16][%lu]\n", (unsigned long)PointerListLen);
    printf("\t-h \thelp\n");
}


extern double **dp;
extern int *v;
extern size_t max_dp_ref_len, max_dp_read_len;
// int max_error = 50, max_pattern=5000;
extern LV_ENTITY *lv;




int main(int argc, char** argv) 
{ 
    char d_n[500]="";
    char r_n[500]="";
    char v_n[500]="";
    char r_f[500]="";
    char r_h[500]="";
    char o_n[500]="out.sam";
    
    int c;
    opterr = 0;
    while ((c = getopt (argc, argv, "v:p:w:f:r:o:h")) != -1) 
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
            case 'v':
                if (optarg)
                    strcpy(v_n, optarg);
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

    if (PointerListLen<6 || PointerListLen > 16)
    {
        std::cout << "PointerListLen out of range[6-16]: " << PointerListLen << endl;
        return 1;
    }
    std::cout << "PointerListLen: " << PointerListLen << " " 
              << "| WindowListLen: " << WindowListLen << " "
              << "| dna file: " << string(d_n) << " " 
              << "| read file: " << string(r_n) << endl;


    //read ref

    FILE *out;
    out = fopen(o_n, "wb");

    int sam = 1;
    if (sam)    fprintf(out, "@HD\tVN:0.1.1\n");


    ref_s rs_entity;
    ref_s *rs = &rs_entity;

    //init ref space
    size_t rl=1, nl=1;
    int l;
    gzFile fp;
    kseq_t *seq;
    strcpy(r_f, d_n);
    strcat(r_f, ".hash_ref");
    fp = gzopen(r_f, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        if (sam) fprintf(out, "@SQ\tSN:%s\tLN:%d\n", seq->name.s, (int32_t)seq->seq.l);
        printf("%2lu %13s len: %lu i:%lu\n", nl, seq->name.s, seq->seq.l, rl);
        // if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
        // printf("seq: %s\n", seq->seq.s);
        // if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        rl += seq->seq.l;
        ++nl;
    }
    printf("total Len:%lu\n", rl);
    rs->ref_s_l=0;
    rs->ref_s = (uint64_t *)malloc(sizeof(uint64_t) * (rl/32+1));
    rs->uref_l=0;
    rs->uref = (ref_type *)malloc(sizeof(ref_type) * nl);
    char *ref_s = (char *)malloc((rl+1) * sizeof(char));
    
    // ref_buffer += (sizeof(uint64_t) * (rl/32+1) + sizeof(ref_type) * nl);
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp);

    // //read ref
    // // read_ref(&rs, d_n);
    uint64_t unit_number = 0;
    fp = gzopen(r_f, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        //set ference value
        ref_type urt;
        char *seq_n = (char *)malloc(sizeof(char) * (seq->name.l + 1));
        strcpy(seq_n, seq->name.s);
        char *seq_c = (char *)malloc(sizeof(char) * (seq->comment.l + 1));
        strcpy(seq_c, seq->comment.s);
        urt.n = seq_n;
        urt.c = seq_c;
        urt.l = seq->seq.l;
        urt.i = rs->ref_s_l;

        // char2char(seq->seq.s);
        strncpy(ref_s+rs->ref_s_l, seq->seq.s, seq->seq.l);

        ref_set2int(seq->seq.s, seq->seq.l, rs->ref_s, &(rs->ref_s_l));

        rs->uref[rs->uref_l++] = urt;
    }
    ref_s[rs->ref_s_l] = '\0';
    // printf("return value: %d\n", l);
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp);

    char log_ref_file[] = "log_ref";
    print_ref_in_char(rs, log_ref_file);

    
    size_t ref_len = rs->ref_s_l;

    if (sam)
    {
        fprintf(out, "@PG\tID:run_rHAT\tVN:0.1\tCL:");
        for (size_t i=0; i<argc; ++i)
        {
            fprintf(out, i==(argc-1) ? "%s\n":"%s ", argv[i]);
        }
    }

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
    // return 1;


    //read VCF
    var_s vs_entity;
    var_s *vs = &vs_entity;

    int have_vcf = 0;

    if (v_n[0] != '\0') have_vcf = 1;
    if (have_vcf)
    {
        create_HAFF();

        printf("read VCF file from:%s\n", v_n);
        init_var(rs, vs, v_n);
        read_var(rs, vs, v_n);
        
        char log_vcf[50]="vcf.log";
        log_vcf_v(vs, rs, log_vcf);
    }

    // ifstream inRef, inRead;
    
    FILE *inRHT;
    // inRef.open(d_n);
    // inRead.open(r_n);

    strcpy(r_h, d_n);
    strcat(r_h, ".out_RHT");
    printf("read hash from %s\n", d_n);
    inRHT = fopen(r_h, "rb");

    //-----------------------------------get dna_ref string
    // string dna_name, dna_s, dna_f;
    // getline(inRef, dna_name);
    // while (inRef >> dna_s) dna_f.append(dna_s);

    //-----------------------------------get hash table
    // RHT rht(dna_f.size());
    RHT rht(ref_len);
    if (rht.read_hash(inRHT)) return 1;
    cout << "succeed load RHT\n";
    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    
    ofstream outt;
    outt.open("check.out");
    // rht.write_hash_test(outt);

    //-----------------------------------get read
    // string read_m, read, read_x;
    //init mask mask = 0xffffffff
    MASK = MASK >> (32 - PointerListLen * 2);

    // //the hit window number that will go to alignment of each read
    uint32_t full = 3;
    uint32_t *hit_w = (uint32_t *)malloc(full*sizeof(uint32_t));
    uint32_t *hit_times = (uint32_t *)malloc(full*sizeof(uint32_t));
    uint32_t tmp = 0, tar;

    // const char *c_dna_f = dna_f.c_str();
    const char *c_dna_f = ref_s;


    //init max_read length
    size_t max_read = 100;
    fp = gzopen(r_n, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        // if (seq->seq.l > max_read) printf("%s\n", seq->name.s);
        max_read = MX(max_read, seq->seq.l);

    }
    // max_read += 10;
    printf("max read length: %lu\n", (unsigned long)max_read);
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp);

    max_dp_ref_len = WindowListLen + max_read + 5;
    max_dp_read_len = max_read + 5;

    size_t max_window_index = 1 + (ref_len) / WindowListLen * 2 ;
    //#2

    //init dp space
    dp = (double **)malloc(sizeof(double *) * max_dp_ref_len);
    for(size_t i=0; i<max_dp_ref_len; ++i)
    {
        dp[i] = (double *)malloc(sizeof(double) * (max_dp_read_len));
    }
    v = (int *)malloc(sizeof(int) * 2 * (MX(max_dp_read_len, max_dp_ref_len)));

    int max_error = 50, max_node = 500, max_deep = 500;
    if (have_vcf)
    {
        //len node error read deep
        lv = (LV_ENTITY *)malloc(sizeof(LV_ENTITY));
        init_lv_space(lv, max_node, max_error, max_dp_read_len, max_deep);
    }

    int windown_hit_too_low=20, hit_too_many=WindowListLen / 2;

    //ref window count array
    int *WC = (int *)calloc(max_window_index, sizeof(int));

    size_t tmp_wc, max_wc=0;

    //read hash table
    RHT rrht(max_read);
    bool vis[max_read*2 + WindowListLen];

    DAG dag(max_read + WindowListLen + 3);

    Cigar cg, thecg;
    cg.length = 0;
    cg.seq =  (uint32_t *)malloc(max_read * sizeof(uint32_t));
    thecg.length = 0;
    thecg.seq =  (uint32_t *)malloc(max_read * sizeof(uint32_t));
    
    char *read_s;
    size_t read_len = 0;
    gzFile fp2 = gzopen(r_n, "r");
            
    char *seq_c = (char *)malloc(sizeof(char) * (max_read * 2 + 100));

    char *c_read_4bit;
    if (have_vcf)
    {
        c_read_4bit = (char *)malloc(sizeof(char) * (max_read * 2 + 100));
    }
    // PX("asd");
    // cout << max_read + WindowListLen + 3 << endl;
    // return 1;
    int best_if_reverse = 0, find_a = 0;
    size_t best_pos = 0;
    if (!fp2) fprintf(stderr, "read file %s error\n", r_n);
    else
    {
        kseq_t *seq2 = kseq_init(fp2);
        while ((l = kseq_read(seq2)) >= 0) 
        {

            outt << seq2->name.s << " ";


            //print read info from simulated read data
            #ifdef READINFO
            strcpy(seq_c, seq2->comment.s);
            char *p = strtok(seq_c, " =");
            char *p1, *p2;
            char *config, *chrom, *hap_ref, *edit_s;
            size_t begin, end;
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
            outt << begin << " " << edit_s << " ";
            #endif


            size_t extend_time = 0;
            char2char(seq2->seq.s);
            // // ref = seq1->seq.s;
            read_len = (seq2->seq.l);
            // const char *c_read = read.c_str();

            char *c_read = seq2->seq.s;
            char *re_read;
            bool reverse = false;
            // size_t read_len = read.size();
            
            
            double score = -INF, score2, score_enough=read_len*scy*0.95;
            double score_too_low=(-0.8)*read_len*scy;
            
            best_if_reverse = 0;
            find_a = 0;

            while (true)
            {
                //-----------------------------------step 1 find hit windows 
                //Hwin is the first heap to set window priority_queue

                //reverse read
                if (reverse)
                {
                    re_read =  (char *)malloc((read_len+1) * sizeof(char));
                    strcpy(re_read, seq2->seq.s);

                    for (size_t i = 0; i<read_len/2; ++i)
                    {   
                        char tmp = re_read[i];
                        re_read[i] = re_read[read_len - 1 - i];
                        re_read[read_len - 1 - i] = tmp;
                    }

                    for (size_t i = 0; i<read_len; ++i)
                    {
                        switch(re_read[i])
                        {
                            case 'A':
                                re_read[i] = 'T';
                                break;
                            case 'T':
                                re_read[i] = 'A';
                                break;
                            case 'G':
                                re_read[i] = 'C';
                                break;
                            case 'C':
                                re_read[i] = 'G';
                                break;    
                        }
                    }
                    c_read = re_read;
                }

                if (have_vcf)
                {
                    parse_pattern_string_4bit(c_read, read_len, c_read_4bit);
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
                        outt << mx << " window hit too low reverse read or end read alignment\n";
                    #endif
                    if (reverse) break;
                    reverse = true;
                    continue;
                }
                //#


                //-----------------------------------step 2 create read hash table
                tmp = 0;
                rrht.clear_PW();
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
                // rrht.create_p_w();

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

                    // memset(vis, 0, sizeof(vis));
                    // memset(vis_m, 0, sizeof(vis_m));

                    size_t pw, n, ttmp;
                    uint32_t pr;
                    uint32_t d = 0;
                    
                    priority_queue<WQ> Q;
                    WQ qt;
                    tmp = 0;
                    for (size_t i=window_up; i < window_down; ++i)
                    {
                        tmp = (tmp << 2) | get_c[c_dna_f[i]];
                        tar = tmp & MASK;
                        // outt << "**" << endl;

                        if (i - window_up >= PointerListLen - 1)
                        {
                            //if find window & read match l-mer
                            // d = 0;
                            d = rrht.search(tar);
                            // continue;
                            // printf("%lu\n", (unsigned long)d);
                            // outt << "|";

                            while (d && (rrht.PW[d-1] >> 32) == tar)
                            {   
                                // outt << "-";
                                // outt << ">" << (uint32_t)rrht.PW[d-1] << " " << to_string((rrht.PW[d-1] >> 32)) << endl;
                                // outt << (i - window_up + read_len - (uint32_t)rrht.PW[d-1]) << endl;

                                pw = i;
                                pr = rrht.PW[d-1];

                                if (!vis[(i - window_up + read_len - pr)])
                                {
                                // outt << ">" << (uint32_t)rrht.PW[d-1] << " " << to_string((rrht.PW[d-1] >> 32)) << endl;
                                    // outt << "asd" << endl;
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


                            // d = rrht.P[tar];
                            // while (d && rrht.P[tar+1]>d)
                            // {
                            //     pr = rrht.W[d-1];
                            //     pw = i;
                            //     if (!vis[(i - window_up + read_len - pr)])
                            //     {
                            //     // outt << ">" << (uint32_t)rrht.PW[d-1] << " " << to_string((rrht.PW[d-1] >> 32)) << endl;
                            //         // outt << "asd" << endl;
                            //         n = dag.add_node(pw, pr, PointerListLen);

                            //         ttmp = (i - window_up + read_len - pr);

                            //         //if overlap
                            //         vis[ttmp]=true;
                            //         // outt << "B: " << pr << " "<<(i - window_up + read_len - pr) << " ";
                            //         ++pw;
                            //         ++pr;
                            //         while(c_dna_f[pw]==c_read[pr] && dag.check_node(pw, pr, n))
                            //         {
                            //             // vis[(pw - window_up + read_len - pr)]=true;
                            //             ++pw;
                            //             ++pr;
                            //             // outt << (unsigned long)pw - window_up << "/" << pr << " ";
                            //         }
                            //         // for (size_t t = 1; t<PointerListLen-1; ++t) vis[pr + t] = true;
                            //         qt.t = ttmp;
                            //         qt._i = pw - 1;
                            //         // vis_m[pw-window_up] = 1 + ttmp;
                            //         Q.push(qt);
                            //         // outt << wd[node_i].index_of_W << " " << wd[node_i].index_of_R << " " << wd[node_i].len << " ";
                            //         #ifdef PRINTLOG
                            //         // outt << "+ " << pr << " "<< qt._i - window_up << " ";
                            //         // cout_string_in_char(c_read, (dag.wd)[(dag.node_i-1)].index_of_R - PointerListLen + 1, (dag.wd)[(dag.node_i-1)].len, outt);
                            //         // outt << endl; 
                            //         #endif
                            //     }
                            //     ++d;
                            // }
                            //free match pair block
                            while (!Q.empty() && Q.top()._i == i)
                            {
                                // outt << "F " << c_dna_f[i] << " " << i - window_up << "> \n" ;

                                qt = Q.top();
                                Q.pop();
                                vis[qt.t] = false;
                            }

                        }
                    }
                    dag.add_node(window_down + PointerListLen - 1, read_len + PointerListLen - 1, 0);
                    #ifdef PRINTLOG
                         outt << "add  node end\n";
                    #endif

                    double t = -INF;
                    // dag.print_log(outt);

                    // linking nodes edge
                    // node matrix 
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

                    
                    cg.length = 0;
                    size_t pp = 0;
                    t = dag.do_alignment(const_cast<char*>(c_dna_f), window_up, window_down, const_cast<char*>(c_read), read_len, &cg, &pp, out);

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


                    if (t > score)
                    {
                        score2 = score;

                        find_a = 1;
                        best_pos = pp;
                        if (reverse == 1)   best_if_reverse = 1;
                        score = t;

                        copy_cigar(&thecg, &cg);
                    }

                    if (score < score_too_low)
                    {
                        #ifdef PRINTLOG
                            outt << "break for best score too low\n";
                        #endif
                        break;
                    }
                    else if (z > 0 && hit_times[z] > hit_times[z-1])
                    {
                        #ifdef PRINTLOG
                            outt << hit_times[z] << " " << hit_times[z-1] << "\n";
                            outt << "break for score no increse\n";
                        #endif
                        break;
                    }
                    
                    

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
                        outt << "no reverse need for score_enough\n";
                    #endif
                    break;
                }
                
                // break;
                if (reverse) break;
                reverse = true;
            }


            //print read alignment result
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

            if (sam)
            {
                //read name
                fprintf(out, "%s\t", seq2->name.s);
                //flag
                int sam_flag = 0;
                if (!find_a)
                {
                    fprintf(out, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
                    continue;
                }
                if (best_if_reverse)    sam_flag += 16;
                fprintf(out, "%d\t", sam_flag);
                //ref name
                for (size_t i = (rs->uref_l)-1; i>=0; --i)
                {
                    if ((rs->uref)[i].i <= best_pos)
                    {
                        fprintf(out, "%s\t", (rs->uref)[i].n);
                        break;
                    }
                }
                // if (best_pos == 0) best_pos = (rs->uref)[0].n;
                //pos 1 base
                fprintf(out, "%lu\t", (unsigned long)(best_pos+1));
                //maq
                /*
                MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest
                integer. A value 255 indicates that the mapping quality is not available.
                */
                
                // int tot_n=0, m_n=0, t_n;
                // for (size_t i = 0; i<(thecg.length); ++i)
                // {
                //     t_n = cigar_int_to_len(thecg.seq[i]);
                //     tot_n += t_n;
                //     if ( cigar_int_to_op(thecg.seq[i]) == 'M')
                //         m_n += t_n;
                // }
                // if (tot_n == m_n)
                //     fprintf(out, "255\t");
                // else
                //     fprintf(out, "%d\t", (int)((-10)*log10((tot_n - m_n)*1.0/tot_n)));
                
                // uint32_t mapq = -4.343 * log(1 - (double)abs(score - score2)/(double)score);
                // mapq = (uint32_t) (mapq + 4.99);
                // mapq = mapq < 254 ? mapq : 254;
                // fprintf(out, "%d\t", mapq);

                fprintf(out, "%d\t", 255);
                //cigar
                for (size_t i = 0; i<(thecg.length); ++i)
                {
                    fprintf(out, "%d%c", cigar_int_to_len(thecg.seq[i]), \
                        cigar_int_to_op(thecg.seq[i]));
                }
                fprintf(out, "\t");
                //ref next
                //position next
                //observed Template LENgth
                fprintf(out, "*\t0\t0\t");
                //segment SEQuence
                if (best_if_reverse)    fprintf(out, "%s\t", re_read);
                else fprintf(out, "%s\t", seq2->seq.s);
                //QUAL 
                if (seq2->qual.s) fprintf(out, "%s", seq2->qual.s);
                else fprintf(out, "*");
            
                if (score) fprintf(out, "\tAS:i:%lu", (unsigned long)score);

                int have_v = 0;
                if (have_v)
                {
                    fprintf(out, "\tVP:A:pos_vtype");
                }
            //     if (a->score2 > 0) fprintf(stdout, "ZS:i:%d\n", a->score2);
            //     else fprintf(stdout, "\n");

                fprintf(out, "\n");
            }
            // if (re_read) free(re_read);
        }
        kseq_destroy(seq2);
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