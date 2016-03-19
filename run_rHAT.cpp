#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "RHT.h"
using namespace std;


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
    inf.open("GCF_000005845.2_ASM584v2_genomic.fna");
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

    index_r = 0;
    while (index_r < dna_read.size()) 
    {
        ++index_r;
        if (index_r >= PointerListLen)
        {
            L_read[index_r - PointerListLen] = to_bit(dna_read.substr(index_r - PointerListLen, index_r));
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
        out << w.index_of_L << " " << to_string(L_read[w.index_of_L]) << " " << w.count_L << " " << w.index_of_W << endl;

        if (last_W_index == w.index_of_W)
        {
            ++wc.count;
        }
        else
        {
            ++wc_num;
            QWin_C.push(wc);
            if (wc_num> wc_full)
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
    }

    while (!Qwin.empty())
    {
        cp_window(w, Qwin.top());
        out << w.index_of_L << " " << to_string(L_read[w.index_of_L]) << " " << w.count_L << " " << w.index_of_W << endl;

        if (last_W_index == w.index_of_W)
        {
            ++wc.count;
        }
        else
        {
            ++wc_num;
            QWin_C.push(wc);
            if (wc_num> wc_full)
            {
                --wc_num;
                QWin_C.pop();
            }
            last_W_index = w.index_of_W;
            wc.index_of_W = w.index_of_W;;
            wc.count = 1;
        }

        Qwin.pop();
    }

    //the top 5 hit window index store in hit_w
    uint32_t k=wc_full;
    while (!QWin_C.empty())
    {
        out << QWin_C.top().count << " " << QWin_C.top().index_of_W << endl;
        hit_w[--k] = QWin_C.top().index_of_W;
        QWin_C.pop();

    }

    for (uint32_t i=0; i<wc_full; i++)
        out << hit_w[i] << endl;


    /*for (size_t i=0; i <= index_r - PointerListLen; ++i)
    {
       out << L_read[i] << " " << to_string(L_read[i]) << "\n";
    }*/




    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
	return 0;
}