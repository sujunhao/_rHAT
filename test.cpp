#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <sstream>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "RHT.h"
#include "alignment.h"
using namespace std;

// #define PRINTLOG

typedef struct window {
    uint32_t index_of_W;
    size_t a, b;

    bool operator<(const window& w) const
    {
        return index_of_W >= w.index_of_W;
    }
}WINDOW;

typedef struct window_cnt {
    uint32_t index_of_W;
    size_t cnt;
    bool operator<(const window_cnt& w) const
    {
        return cnt > w.cnt;
    }
}WINDOW_CNT;

typedef struct window_node {
    uint32_t index_of_W, index_of_R;
    size_t len;
}WINDOW_NODE;


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

    ifstream inRef, inRead;
    ofstream outt;
    FILE *out;
    FILE *inRHT;
    inRef.open("E.coli.fa");
    inRead.open("dna_read");
    outt.open("outz");
    out = fopen("out", "wb");
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

    while (inRead >> read_m)
    {
        //-----------------------------------get each read seq, mark, and info
        inRead >> read;
        inRead >> read_m;
        inRead >> read_x;
        // outt << read << endl;

        RHT rrht(read.size());
        bool reverse = 0;
        while (true)
        {
            //-----------------------------------step 1 find hit windows  (tot 8s)
            size_t read_len = read.size();
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

            //Hwin is the first heap to set window priority_queue
            priority_queue<WINDOW> Hwin;
            WINDOW w, w_tmp;
            tmp = 0;
            for (size_t i=len_up_stream; i < len_up_stream + theLen; ++i)
            {
                tmp = (tmp << 2) | get_c[read[i]];
                tar = tmp & MASK;
                if (i - len_up_stream >= PointerListLen - 1)
                {
                    if (rht.P[tar] && rht.P[tar+1] && rht.P[tar+1] > rht.P[tar])
                    {
                        w.index_of_W = rht.W[rht.P[tar] - 1];
                        // outt << read_m << " " << w.index_of_W << endl;
                        w.a = rht.P[tar] + 1;
                        w.b = rht.P[tar+1];
                        Hwin.push(w);
                    }
                }
            }
            //Hwin_cnt is the second heap to count window hit number
            priority_queue<WINDOW_CNT> Hwin_cnt;
            WINDOW_CNT wc;
            uint32_t last = w.index_of_W, thecnt = 0;
            while (!Hwin.empty())
            {
                w = Hwin.top();
                //do Hwin pop process
                Hwin.pop();
                if (w.b > w.a)
                {
                    w_tmp.index_of_W = rht.W[w.a - 1];
                    w_tmp.a = w.a + 1;
                    w_tmp.b = w.b;
                    Hwin.push(w_tmp);
                }

                //if a new window num appear, push into hwin_cnt
                if (w.index_of_W != last)
                {
                    wc.index_of_W = last;
                    wc.cnt = thecnt;
                    Hwin_cnt.push(wc);
                    if (Hwin_cnt.size() > full) Hwin_cnt.pop();
                    last = w.index_of_W;
                    thecnt = 1;
                }
                else 
                {
                    ++thecnt;
                }
            }
            uint32_t z = full;
            while(!Hwin_cnt.empty())
            {
                hit_w[--z] = Hwin_cnt.top().index_of_W;
                
                #ifdef PRINTLOG
                //print the hit window and their hit times
                    outt << Hwin_cnt.top().cnt << " " << Hwin_cnt.top().index_of_W << endl;
                #endif
                Hwin_cnt.pop();
            }

           
            //-----------------------------------step 2 create read hash table
            tmp = 0;
            rrht.clear();
            for (size_t i=0; i < read.size(); ++i)
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
            for (size_t z=0; z<full && z<1; ++z)
            {
                //-----------------------------------struct window and read main body
                size_t window_up, window_down;
                window_up = ( hit_w[z] * (WindowListLen / 2) >= len_up_stream ) ? ( hit_w[z] * (WindowListLen / 2) - len_up_stream ) : (0);
                window_down = ( hit_w[z] * (WindowListLen / 2) + WindowListLen + len_down_stream <= dna_f.size() ) ? ( hit_w[z] * (WindowListLen / 2) + WindowListLen + len_down_stream ) : (dna_f.size());
                // outt << dna_f.substr(window_up, window_down - window_up) << endl;

                // struct graph node include the first and last node
                uint32_t d;
                size_t node_i=0;
                WINDOW_NODE *wd, wd_tmp;
                wd = new WINDOW_NODE[window_down - window_up + 3];
                //the first node
                wd_tmp.index_of_W = window_up;
                wd_tmp.index_of_R = PointerListLen - 1;
                wd_tmp.len = 0;
                wd[node_i++] = wd_tmp;
                tmp = 0;
                for (size_t i=window_up; i < window_down; ++i)
                {
                    tmp = (tmp << 2) | get_c[dna_f[i]];
                    tar = tmp & MASK;
                    if (i - window_up >= PointerListLen - 1)
                    {
                        //if find window & read match l-mer
                        d = rrht.search(tar);
                        if (d)
                        {
                            // outt << to_string(tar) << " " << d << endl;
                            while ((rrht.PW[d-1] >> 32) == tar)
                            {
                                //i is tha match window string last index
                                wd_tmp.index_of_W = i;
                                wd_tmp.index_of_R = rrht.PW[d-1];
                                wd_tmp.len = PointerListLen;
                                ++d;
                                //if overlap
                                if (node_i > 1 && i > wd[node_i-1].index_of_W && (i - wd[node_i-1].index_of_W) == wd[node_i-1].len - PointerListLen + 1)
                                {
                                    ++wd[node_i-1].len;
                                    continue;
                                }
                                wd[node_i] = wd_tmp;
                                // outt << wd[node_i].index_of_W << " " << wd[node_i].index_of_R << " " << wd[node_i].len << " " << read.substr(wd[node_i].index_of_R, wd[node_i].len) << endl; 
                                ++node_i;
                            }
                        }
                    }
                }
                wd_tmp.index_of_W = window_down;
                wd_tmp.index_of_R = read.size();
                wd_tmp.len = 0;
                wd[node_i++] = wd_tmp;
                // for (size_t i = 0; i<node_i; ++i)
                // {
                //     outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
                // }
                // outt << endl;


                //linking nodes edge
                // uint32_t t_wait = 1024;
                // uint32_t iw, ir, il, jw, jr;
                // for (size_t i=0; i<node_i-1; ++i)
                // {
                //     iw = wd[i].index_of_W;
                //     ir = wd[i].index_of_R;
                //     il = wd[i].len;
                //     for (size_t j=i+1; j<node_i; ++j)
                //     {
                //         jw = wd[j].index_of_W;
                //         jr = wd[j].index_of_R;
                //         if (jw >= il + iw && jr >= il + ir && jr <= ir + il + t_wait)
                //         {

                //         }
                //     }
                // }

                //-----------------------------------use global & semiglobal alignment to struct alignment and get sroce
                
                delete [] wd;
            }



            if (reverse == 1) break;
            reverse = 1;
        }

        //-----------------------------------out put the read best alignment

    }

    
    

    

    inRef.close();
    inRead.close();
    outt.close();
    fclose(out);
    fclose(inRHT);


        


    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}