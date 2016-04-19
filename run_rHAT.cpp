#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <sstream>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "DAG.h"
using namespace std;

// #define PRINTLOG

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


int main(int argc, char** argv) 
{ 
    int c;
    opterr = 0;
    char d_n[100]="E.coli.fa";
    char r_n[100]="dna_read";
    char o_n[100]="out";
    while ((c = getopt (argc, argv, "p:w:d:r:o:h")) != -1) 
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
            case 'r':
                if (optarg)
                    strcpy(r_n, optarg);
            case 'o':
                if (optarg)
                    strcpy(o_n, optarg);   
            case 'h':
                printf("\tuse out_RHT to run rHAT\n");
                printf("\trun_rHAT [option]\n");
                printf("\t-d \tdna file           [%s]\n", d_n);
                printf("\t-r \tread file          [%s]\n", r_n);
                printf("\t-o \toutput file        [%s]\n", o_n);
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
              << "| dna file: " << string(d_n) << " " 
              << "| read file: " << string(r_n) << endl;

    ifstream inRef, inRead;
    ofstream outt;
    FILE *out;
    FILE *inRHT;
    inRef.open(d_n);
    inRead.open(r_n);
    outt.open("outt");
    out = fopen(o_n, "wb");
    inRHT = fopen("out_RHT", "rb");

    //-----------------------------------get dna_ref string
    string dna_name, dna_s, dna_f;
    getline(inRef, dna_name);
    while (inRef >> dna_s) dna_f.append(dna_s);

    //-----------------------------------get hash table
    RHT rht(dna_f.size());
    rht.read_hash(inRHT);
    cout << "succeed load RHT\n";
    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    // rht.write_hash_test(outt);

    //-----------------------------------get read
    string read_m, read, read_x;
    //init mask
    MASK = MASK >> (32 - PointerListLen * 2);
    //the hit window number that will go to alignment of each read
    const uint32_t full = 3;
    uint32_t hit_w[full];
    uint32_t tmp = 0, tar;

    const char *c_dna_f = dna_f.c_str();

    char sss[10000];
    char thesss[10000];


    //#2
    size_t WC[10000];
    memset(WC, 0, sizeof(WC));
    size_t tmp_wc, max_wc=0;
    //#
    RHT rrht(30000);
    bool vis[30000];


    while (inRead >> read_m)
    {
        //-----------------------------------get each read seq, mark, and info
        inRead >> read;
        inRead >> read_m;
        inRead >> read_x;
        // if (read.size() > 10000)
        // cout << read.size() << endl;
        // outt << read << endl;

        const char *c_read = read.c_str();
        bool reverse = false;
        size_t read_len = read.size();
        
        
        double score = -INF, score_enough=100;
        double score_too_low=-100;
        size_t windown_hit_too_low=150, hit_too_many=500;
        string thess(read_len+WindowListLen, 0);
        thess.clear();

        while (true)
        {
            //-----------------------------------step 1 find hit windows 
            //Hwin is the first heap to set window priority_queue

            //reverse read
            if (reverse)
            {
                std::reverse(read.begin(),read.end());
                for (size_t z = 0; z < read.size(); ++z)
                {
                    switch(read[z])
                    {
                        case 'A':
                            read[z] = 'T';
                            break;
                        case 'T':
                            read[z] = 'A';
                            break;
                        case 'G':
                            read[z] = 'C';
                            break;
                        case 'C':
                            read[z] = 'G';
                            break;    
                    }
                }
            }
            #ifdef PRINTLOG
                outt << read << endl;
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
            //#


            //#2 get the most hitted window
            tmp = 0;
            max_wc = 0;
            for (size_t i=len_up_stream, k; i < len_up_stream + theLen; ++i)
            {
                tmp = (tmp << 2) | get_c[read[i]];
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
                tmp = (tmp << 2) | get_c[read[i]];
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
                if (abs(last_z-hit_w[z]) <= 1) break;
                last_z = hit_w[z];
                size_t window_up, window_down;
                window_up = ( hit_w[z] * (WindowListLen / 2) >= len_up_stream ) ? ( hit_w[z] * (WindowListLen / 2) - len_up_stream ) : (0);
                window_down = ( hit_w[z] * (WindowListLen / 2) + WindowListLen + len_down_stream <= dna_f.size() ) ? ( hit_w[z] * (WindowListLen / 2) + WindowListLen + len_down_stream ) : (dna_f.size());
                #ifdef PRINTLOG
                outt << "\n" << window_up << endl;
                outt << dna_f.substr(window_up, window_down - window_up) << endl;
                #endif
                // struct graph node include the first and last node
                // dag.clear();
                DAG dag(window_down - window_up + 3);
                dag.clear();
                //the first node
                dag.add_node(window_up + PointerListLen - 1, PointerListLen - 1, 0);

                tmp = 0;
                memset(vis, 0, sizeof(vis));
                queue<WQ> Q;
                WQ qt;
                for (size_t i=window_up; i < window_down; ++i)
                {
                    tmp = (tmp << 2) | get_c[dna_f[i]];
                    tar = tmp & MASK;
                    if (i - window_up >= PointerListLen - 1)
                    {
                        //if find window & read match l-mer
                        uint32_t d;
                        d = rrht.search(tar);
                        size_t pw, n;
                        uint32_t pr;
                        if (d)
                        {
                            // outt << to_string(tar) << " " << d << endl;
                            while ((rrht.PW[d-1] >> 32) == tar && !vis[(uint32_t)rrht.PW[d-1]])
                            {
                                pw = i;
                                pr = rrht.PW[d-1];
                                n = dag.add_node(pw, pr, PointerListLen);

                                qt.i = i;
                                qt.d = d;
                                qt.tar = tar;
                                qt.n = n;

                                //if overlap
                                vis[pr]=true;
                                ++pw;
                                ++pr;
                                while(!vis[pr] && dna_f[pw]==read[pr] && dag.check_node(pw, pr, n))
                                {
                                    vis[pr]=true;
                                    ++pw;
                                    ++pr;
                                    // printf("%lu\n", (unsigned long)pw);
                                }
                                qt._i = pw;
                                Q.push(qt);
                                ++d;
                                // outt << wd[node_i].index_of_W << " " << wd[node_i].index_of_R << " " << wd[node_i].len << " " << read.substr(wd[node_i].index_of_R, wd[node_i].len) << endl; 
                            }
                        }

                        


                        if (i==Q.front()._i)
                        {
                            qt = Q.front();
                            Q.pop();
                            d = qt.d;
                            while ((rrht.PW[d-1] >> 32) == qt.tar)
                            {

                                pw = qt.i;
                                pr = rrht.PW[d-1];
                                //if overlap
                                vis[pr]=false;
                                ++pw;
                                ++pr;
                                while(vis[pr] && dna_f[pw]==read[pr] && dag.check_node(pw, pr, qt.n, 1))
                                {
                                    vis[pr]=false;
                                    ++pw;
                                    ++pr;
                                }
                                ++d;
                                // outt << wd[node_i].index_of_W << " " << wd[node_i].index_of_R << " " << wd[node_i].len << " " << read.substr(wd[node_i].index_of_R, wd[node_i].len) << endl; 
                            }

                        }
                    }
                }
                dag.add_node(window_down, read_len, 0);

                //dag.print_log(outt);

                //linking nodes edge
                //node matrix 
                dag.create_matrix();
                //do use node matrix to run dp
                dag.find_path();

                //-----------------------------------use global & semiglobal alignment to struct alignment and get sroce
                #ifdef PRINTLOG
                    outt << read_m << endl;
                    dag.print_log(dna_f, read, outt);
                #endif

                double t;
                string ss;
                // // dag.do_alignment(dna_f, read, outt);
                // t = dag.do_alignment(dna_f, window_up, window_down, read, ss, outt);

                // strcpy(sss, "");
                // t = dag.do_alignment(const_cast<char*>(c_dna_f), window_up, window_down, const_cast<char*>(c_read), read_len, sss, out);
                t = dag.do_alignment(const_cast<char*>(c_dna_f), window_up, window_down, const_cast<char*>(c_read), read_len, ss, out);

                if ((t) > score)
                {
                    score = t;

                    // memset(thesss, 0, sizeof(thesss));
                    // strcpy(thesss, sss);
                    
                    thess = ss;
                }
                else break;

                if (score < score_too_low) break;
                // dag.do_alignment(c_dna_f, window_up, window_down - PointerListLen + 1, c_read, out);

                if (t <= score_too_low) break;
                if (score >= score_enough) break;   
            }

            if (score >= score_enough) break;
            if (reverse) break;
            reverse = true;
        }



        //-----------------------------------out put the read best alignment
        fprintf(out, "%s\n", read_m.c_str());
        fprintf(out, "%d\n", reverse);
        fprintf(out, "%.0lf\n", score);
        // fprintf(out, "%s\n", thesss);
        fprintf(out, "%s\n", thess.c_str());


        // printf("%.0lf\n", score);
        // printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);

    }

    
    

    

    inRef.close();
    inRead.close();
    outt.close();
    fclose(out);
    fclose(inRHT);


        


    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}