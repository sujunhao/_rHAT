#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "RHT.h"
#include "global_alignment.h"
using namespace std;

// #define PRINTLOG

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

typedef struct read_hash
{ 
    vector<size_t> index;
}READ_HASH;


bool cp_window(WINDOW &a, const window& w) 
{
    a.index_of_L = w.index_of_L;
    a.index_of_W = w.index_of_W;
    a.count_L = w.count_L;
}


typedef struct match_point
{
    size_t len;
    uint64_t index_of_W;
    uint32_t index_of_R;
    vector<size_t> U;
}MATCH_POINT;

vector<uint32_t> dp;
vector<bool> vis;

uint32_t find_path_mpv(vector<MATCH_POINT> mpv, size_t v)
{
    if (vis[v]) return dp[v];
    uint32_t value = 0;
    for (size_t i = 0; i < mpv[v].U.size(); ++i)
    {
        if (vis[mpv[v].U[i]]) value = max(value, dp[mpv[v].U[i]] + (uint32_t)mpv[v].len);
        else
        {
            value = max(value, find_path_mpv(mpv, mpv[v].U[i]) + (uint32_t)mpv[v].len);
        }
    }
    vis[v] = true;
    // std::cout << v << " " << value << endl;
    return dp[v] = value;
}

void get_path_mpv(vector<MATCH_POINT> mpv, size_t v, vector<size_t>& out_path)
{
    for (size_t i = 0; i < mpv[v].U.size(); ++i)
    {
        if (dp[mpv[v].U[i]] + (uint32_t)mpv[v].len == dp[v])
        {
            get_path_mpv(mpv, mpv[v].U[i], out_path);
            out_path.push_back(v);
            return;
        }
    }
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


    dna_bitset db(PointerListLen);
    db.read_hash_in(inRHT);
    cout << "succeed load RHT\n";
    // db.write_hash_out(outRHT);

    priority_queue<WINDOW> Qwin;
    priority_queue<WINDOW_CNT> QWin_C;
    WINDOW w, w1;
    WINDOW_CNT wc, wc1;
    num_w = 0;

    inRead >> dna_read;
    out << dna_read << endl << endl;
    //create_read_hash_table
    const size_t readListlen = 11;
    size_t read_hash_len = 1 << (2 * readListlen);
    vector<READ_HASH> read_h(read_hash_len);

    index_r = 0;
    while (index_r < dna_read.size())
    {
        ++index_r;
        if (index_r >= readListlen)
        {
            read_h[to_bit(dna_read.substr(index_r - readListlen, readListlen))].index.push_back(index_r - readListlen);
            // out << dna_read.substr(index_r - readListlen, readListlen) << " " << to_string(to_bit(dna_read.substr(index_r - readListlen, readListlen))) << " " << index_r - readListlen << endl;
        }
    }
    //inside the read_h store vector if some k_mer string appear in read

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

    for (uint32_t o = 0; o < wcfull; ++o)
    {
        if (o>0) break;
        string ss = ws[o].s;
        out << ss << endl;
        size_t count = ws[o].count, window_len = ws[o].len;
        uint64_t index_up = ws[o].index_up;
        uint32_t index_w = ws[o].index_of_W;
        vector<MATCH_POINT> mpv;
        MATCH_POINT mp_tmp;

        for (size_t i = 0; i < ss.size(); ++i)
        {
            if (i >= readListlen)
            {
                if (read_h[to_bit(ss.substr(i - readListlen, readListlen))].index.size() > 0)
                {
                    size_t k = 0;
                    vector<size_t> vv;
                    vv = read_h[to_bit(ss.substr(i - readListlen, readListlen))].index;
                    while (k < vv.size())
                    {
                        size_t match_len = readListlen;
                        for (size_t j = 0; i + j < ss.size(); ++j)
                        {
                            if (dna_read[vv[k] + readListlen + j] == ss[i + j])
                                ++match_len;
                            else 
                                break;
                        }
                        mp_tmp.len = match_len;
                        mp_tmp.index_of_W = i - readListlen;
                        mp_tmp.index_of_R = vv[k];
                        // out << mp_tmp.len \
                        // << " " << mp_tmp.index_of_W \
                        // << " " <<  ss.substr(i - readListlen, mp_tmp.len)\
                        // << " " << mp_tmp.index_of_R \
                        // << endl;
                        mpv.push_back(mp_tmp);
                        ++k;
                    }
                }
            }
        }

        //init V_start, V_end
        MATCH_POINT Vs, Ve;
        Vs.len = Ve.len = 0;
        Vs.index_of_R = Vs.index_of_W = 0;
        Ve.index_of_W = window_len;
        Ve.index_of_R = dna_read.size();
        mpv.push_back(Vs);
        mpv.push_back(Ve);

        out << mpv.size() << " " << Ve.index_of_W << " " << Ve.index_of_R << endl << endl;
        const uint32_t t_wait = 1024;
        for (size_t i = 0; i < mpv.size(); ++i)
        {
            for (size_t j = 0; j < mpv.size(); ++j)
            {
                if (i != j)
                if(mpv[i].index_of_R + mpv[i].len <= mpv[j].index_of_R && mpv[i].index_of_R + mpv[i].len + t_wait>= mpv[j].index_of_R  && mpv[i].index_of_W + mpv[i].len <= mpv[j].index_of_W)
                {
                    mpv[j].U.push_back(i);
                    // out << " i i_r i_w" << " " << i \
                    // << " | " << mpv[i].index_of_R  \
                    // << " | " << mpv[i].index_of_W  \
                    // << " | " << j \
                    // << " | " << mpv[j].index_of_R  \
                    // << " | " << mpv[j].index_of_W  \
                    // << " | " << mpv[i].len  \
                    << endl;
                }
            }
        }

        dp.clear();
        vis.clear();
        for (size_t i = 0; i < mpv.size(); ++i) 
        {
            dp.push_back(0);
            vis.push_back(false);
        }
        find_path_mpv(mpv, mpv.size() - 1);
        // for (size_t i = 0; i < mpv.size(); ++i) 
        // {
        //     out << i << " " << dp[i] << endl;
        // }

        vector<size_t> out_path;
        get_path_mpv(mpv, mpv.size() - 1, out_path);

        //S store the output alignment
        global_alignment ga;
        string S1, S2, s1, s2, s3, s4;
        size_t last_r = 0, last_w = 0;
        for (size_t i = 0; i < out_path.size(); ++i)
        {
            size_t z = out_path[i];
            s1.clear();
            s2.clear();
            if (mpv[z].index_of_R > last_r) s1 = dna_read.substr(last_r, mpv[z].index_of_R - last_r);
            if (mpv[z].index_of_W > last_w) s2 = ss.substr(last_w, mpv[z].index_of_W - last_w);
            if (s1.size() || s2.size())
            {
                if (s1.size() != 0 || last_w != 0)
                {
                    ga.set(s1, s2, s3, s4);
                    S1.append(s3);
                    S2.append(s4);
                }
            }
            S1.append(dna_read.substr(mpv[z].index_of_R, mpv[z].len));
            S2.append(ss.substr(mpv[z].index_of_W, mpv[z].len));
            last_r = mpv[z].index_of_R + mpv[z].len;
            last_w = mpv[z].index_of_W + mpv[z].len;
            out << "s1----------\n" << s1 << "\ns2---------\n" << s2 << "\ns3--------\n" << s3 << "\ns4--------\n" << s4 << "\nr----------\n" << dna_read.substr(mpv[z].index_of_R, mpv[z].len) << "\nw---------\n" << ss.substr(mpv[z].index_of_W, mpv[z].len) << "\nS1---------\n" << S1 << "\nS2--------\n" << S2 << "\n---------------\n";
            // out << z
            // << mpv[z].len \
            // << endl \
            // << mpv[z].index_of_R \
            // << " " << dna_read.substr(mpv[z].index_of_R, mpv[z].len) \
            // << endl \
            // << mpv[z].index_of_W \
            // << " " << ss.substr(mpv[z].index_of_W, mpv[z].len) \
            // out << endl << S1 << endl << S2 << endl;
        }
        s1.clear();
        s2.clear();
        if (last_r < dna_read.size()) s1 = dna_read.substr(last_r, dna_read.size() - last_r);
        if (last_w < ss.size()) s2 = ss.substr(last_w, ss.size() - last_w);
        if (s1.size()!= 0)
        {
            ga.set(s1, s2, s3, s4);
            S1.append(s3);
            S2.append(s4);
        }

        for (size_t i = 0; i < S1.size(); ++i)
        {
            if (S1[i] != '_')
            {
                S1 = S1.substr(i, S1.size() - i);
                S2 = S2.substr(i, S2.size() - i);
                break;
            }
        }

        for (size_t i = S1.size() - 1; i >= 0; --i)
        {
            if (S1[i] != '_')
            {
                S1 = S1.substr(0, i+1);
                S2 = S2.substr(0, i+1);
                break;
            }
        }
        out << endl << S1 << endl;
        out << endl << S2 << endl;

    }





    /*for (size_t i=0; i <= index_r - PointerListLen; ++i)
    {
       out << L_read[i] << " " << to_string(L_read[i]) << "\n";
    }*/




    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
	return 0;
}