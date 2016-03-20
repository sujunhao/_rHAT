#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "RHT.h"
using namespace std;

#define PRINTLOG

typedef struct window {
    uint32_t index_of_W;
    size_t index_of_L;
    size_t count_L;

    bool operator<(const window& w) const
    {
        return index_of_W >= w.index_of_W;
    }

}WINDOW;

typedef struct window_cnt {
    uint32_t index_of_W;
    size_t count;

    bool operator<(const window_cnt& w) const
    {
        return count > w.count;
    }

}WINDOW_CNT;


typedef struct window_string {
    uint32_t index_of_W;
    size_t count;
    string s;
    uint64_t index_up;
    size_t len;
}WINDOW_STRING;

typedef struct to_sort 
{ 
    uint32_t w;
    size_t index;  
    bool operator<(const to_sort& b) const
    {
        return w < b.w;
    }
}TO_SORT;


bool cp_window(WINDOW &a, const window& w) 
{
    a.index_of_L = w.index_of_L;
    a.index_of_W = w.index_of_W;
    a.count_L = w.count_L;
}

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

    ifstream inf, inRead, inRHT;
    inf.open("E.coli.fa");
    inRHT.open("out_RHT");
    inRead.open("dna_read");

    ofstream out;
    out.open("out");


    // ofstream outRHT;
    // outRHT.open("out_R");

    string dna_read;
    size_t index_r, num_w;
    uint32_t L_read[WindowListLen];
    uint32_t last_W_index;

    inRead >> dna_read;
    //dna_read is the read length = L/2

    dna_bitset db(PointerListLen);
    db.read_hash_in(inRHT);
    cout << "succeed load RHT\n";
    // db.write_hash_out(outRHT);

    priority_queue<WINDOW> Qwin;
    priority_queue<WINDOW_CNT> QWin_C;
    WINDOW w, w1;
    WINDOW_CNT wc, wc1;
    num_w = 0;


    //splite the read center sequence
    uint32_t read_len = dna_read.size();
    uint32_t len_up_stream = 0, len_down_stream = 0;
    uint32_t theLen = WindowListLen / 2;
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

    string center_read = dna_read.substr(len_up_stream, theLen);

    #ifdef PRINTLOG
    out << center_read << endl;
    #endif

    index_r = 0;
    while (index_r < center_read.size()) 
    {
        ++index_r;
        if (index_r >= PointerListLen)
        {
            L_read[index_r - PointerListLen] = to_bit(center_read.substr(index_r - PointerListLen, index_r));
            if (db.get_count(L_read[index_r - PointerListLen]) > 0)
            {
                w.index_of_L = index_r - PointerListLen;
                w.count_L = 1;
                w.index_of_W = db.get_window(L_read[w.index_of_L], 1);
                num_w += db.get_count(L_read[w.index_of_L]) - 1;
                Qwin.push(w);
            }

        }
    }


    uint32_t wc_num = 0;
    const uint32_t wc_full = 5;
    uint32_t hit_w[wc_full];
    last_W_index = Qwin.top().index_of_W;
    wc.index_of_W = Qwin.top().index_of_W;
    wc.count = 0;

    while (num_w || Qwin.empty())
    {
        cp_window(w, Qwin.top());
        #ifdef PRINTLOG
        out << w.index_of_L << " " << to_string(L_read[w.index_of_L]) << " " << w.count_L << " " << w.index_of_W << endl;
        #endif

        if (last_W_index == w.index_of_W)
        {
            ++wc.count;
        }
        else
        {
            ++wc_num;
            QWin_C.push(wc);
            if (wc_num > wc_full)
            {
                --wc_num;
                QWin_C.pop();
            }
            last_W_index = w.index_of_W;
            wc.index_of_W = w.index_of_W;;
            wc.count = 1;
        }
        Qwin.pop();
        if (w.count_L < db.get_count(L_read[w.index_of_L]))
        {
            w1.index_of_L = w.index_of_L;
            w1.index_of_W = db.get_window(L_read[w.index_of_L], w.count_L + 1);
            w1.count_L = w.count_L + 1;
            Qwin.push(w1);
            --num_w;
        }

        if (Qwin.empty())
        {
            ++wc_num;
            QWin_C.push(wc);
            if (wc_num > wc_full)
            {
                --wc_num;
                QWin_C.pop();
            }
        }
        
    }

    while (!Qwin.empty())
    {
        cp_window(w, Qwin.top());
        #ifdef PRINTLOG
        out << w.index_of_L << " " << to_string(L_read[w.index_of_L]) << " " << w.count_L << " " << w.index_of_W << endl;
        #endif

        if (last_W_index == w.index_of_W)
        {
            ++wc.count;
        }
        else
        {
            ++wc_num;
            QWin_C.push(wc);
            if (wc_num > wc_full)
            {
                --wc_num;
                QWin_C.pop();
            }
            last_W_index = w.index_of_W;
            wc.index_of_W = w.index_of_W;;
            wc.count = 1;
        }

        Qwin.pop();
        if (Qwin.empty())
        {
            ++wc_num;
            QWin_C.push(wc);
            if (wc_num > wc_full)
            {
                --wc_num;
                QWin_C.pop();
            }
        }
    }

    //the top 5 hit window index store in hit_w
    uint32_t wcfull = min(wc_full, wc_num), k=wcfull;

    // cout << wcfull << " " << k << endl;
    WINDOW_STRING ws[wcfull];


    //ws store the most high hit window string, hit time count, etc
    while (!QWin_C.empty())
    {
        #ifdef PRINTLOG
        out << QWin_C.top().count << " " << QWin_C.top().index_of_W << endl;
        #endif
        // hit_w[--k] = QWin_C.top().index_of_W;
        ws[--k].count = QWin_C.top().count;
        ws[k].index_of_W = QWin_C.top().index_of_W;
        QWin_C.pop();
    }



    //ws already have index_w and count info., this function is to find ws array string, index_up and len
    TO_SORT K[wcfull];
    for (size_t i=0; i < wcfull; ++i)
    {
            K[i].w = ws[i].index_of_W * (WindowListLen / 2) > len_up_stream ? ws[i].index_of_W * (WindowListLen / 2) - len_up_stream : 0;
            K[i].index = i;
    }
    sort(K, K+wcfull);
    // for (size_t i=0; i < wcfull; ++i)
    // {
    //     std::cout << K[i].w << " " << K[i].index << endl;
    // }
    
    string dna_name, dna_s, dna_w;

    uint64_t dna_ref_b = 0;
    size_t  index_s = 0, ix_window = 0, thegetwindowlen = len_up_stream + WindowListLen + len_down_stream;

    getline(inf, dna_name);
    std::cout << dna_name  << endl;

// SCANF_WINDOW
    while (inf >> dna_s)
    {
        index_s = 0;
        while (index_s < dna_s.size()) 
        {
            ++dna_ref_b;
            dna_w.append(dna_s, index_s++, 1);
            while (ix_window < wcfull && K[ix_window].w + thegetwindowlen == dna_ref_b)
            {
                ws[K[ix_window].index].s = dna_w.substr(dna_w.size() - thegetwindowlen, thegetwindowlen);
                ws[K[ix_window].index].len = thegetwindowlen;
                ws[K[ix_window].index].index_up = K[ix_window].w;
                ++ix_window;
            }
        }

        if (ix_window >= wcfull) break;
        if (dna_w.size() > 2 * thegetwindowlen) dna_w = dna_w.substr(dna_w.size()-thegetwindowlen, thegetwindowlen);
    }


    #ifdef PRINTLOG
    out << "index_W s count  index_up len\n";
    for (uint32_t i=0; i<wcfull; i++)
    {
        out << ws[i].index_of_W \
        << " " << ws[i].s \
        << " " << ws[i].count \
        << " " << ws[i].index_up \
        << " " << ws[i].len \
        << endl;
    }
    #endif






    /*for (size_t i=0; i <= index_r - PointerListLen; ++i)
    {
       out << L_read[i] << " " << to_string(L_read[i]) << "\n";
    }*/




    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
	return 0;
}