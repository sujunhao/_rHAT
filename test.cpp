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

#define PRINTLOG

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
    outt.open("outt");
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

    while (inRead >> read_m)
    {
        //-----------------------------------get each read seq, mark, and info
        inRead >> read;
        inRead >> read_m;
        inRead >> read_x;
        // outt << read << endl;

        bool reverse = 0;
        while (true)
        {
            //-----------------------------------step 1 find hit windows
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

            uint32_t tmp = 0, tar;
            priority_queue<WINDOW> Hwin;
            WINDOW w, w_tmp;
            for (size_t i=len_up_stream; i < len_up_stream + theLen; ++i)
            {
                tmp = (tmp << 2) | get_c[read[i]];
                tar = tmp & MASK;
                if (i - len_up_stream >= PointerListLen - 1)
                {
                    if (rht.P[tar] && rht.P[tar+1] && rht.P[tar] > rht.P[tar+1])
                    {
                        w.index_of_W = rht.W[rht.P[tar]];
                        w.a = rht.P[tar] + 1;
                        w.b = rht.P[tar+1];
                        Hwin.push(w);
                    }
                }
            }
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
                    w_tmp.index_of_W = w.index_of_W;
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

            //-----------------------------------for each window hit, process to struct main alignment body

                //-----------------------------------use global & semiglobal alignment to struct alignment and get sroce


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