#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <queue>
#include <unistd.h>
#include <stdint.h>
#include "alignment.h"
#define PX(X) std::cout << X << std::endl
using namespace std;

uint32_t t_wait = 1024;


typedef struct window_node {
    uint32_t index_of_W, index_of_R;
    size_t len;
}WINDOW_NODE;

class DAG
{
private:
	WINDOW_NODE *wd;
	size_t **V;
	uint32_t *dp;
	size_t node_i;
	size_t *path;
	size_t p_index;
public:
	DAG(size_t max_num)
	{
		wd = new WINDOW_NODE[max_num];
		node_i = 0;
	}
	~DAG()
	{
		if (V != NULL)
		{
			for (size_t i=0; i<node_i; ++i) delete [] V[i];
	        delete [] V;
		}
		if (dp != NULL)	delete [] dp;
		if (path != NULL) delete [] path;
        delete [] wd;
	}

	void add_node(uint32_t w, uint32_t r, size_t l)
	{
		wd[node_i].index_of_W = w;
		wd[node_i].index_of_R = r;
		wd[node_i].len = l;
		++node_i;
	}
	bool check_node(size_t i)
    {
    	if (node_i > 1 && i > wd[node_i-1].index_of_W && (i - wd[node_i-1].index_of_W) == wd[node_i-1].len - PointerListLen + 1)
    	{
    		++wd[node_i-1].len;
    		return true;
    	}
        return false;
    }

    void print_log(string& dna_f, string& read, ofstream &outt)
    {
    	for (size_t i = 0; i<node_i; ++i)
        {
            outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
            // outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << endl; 
        }
        outt << endl;
    }

    void find_path()
    {
    	dp = new uint32_t[node_i];
    	memset(dp, 0, sizeof(dp));
    	for (size_t i=0; i<node_i-1; ++i)
    	{
    		for (size_t j=i+1; j<node_i; ++j)
    		{
    			if (V[i][j]) dp[j] = (dp[j] > wd[j].len+dp[i] ? dp[j] : wd[j].len+dp[i]);
    		}
    	}
    	size_t thenode = node_i - 1;

    	//path from 1 to p_index) is the index of path string
    	size_t p_tmp;
    	p_index = 0;
    	path = new size_t[node_i];
    	while (thenode != 0)
    	{
    		p_tmp = thenode;
    		for (size_t i = 0; i< node_i; ++i)
    		{
    			if (i==thenode) continue;
    			if (dp[thenode] == wd[thenode].len + dp[i])
    			{
    				thenode = i;
    				// PX(i);
    				path[p_index++] = i;
    				break;
    			}
    		}
    		if(thenode == p_tmp) break;
    	}
    }

    uint32_t do_alignment(string& dna_f, string&  read, ofstream& outt)
    {
    	for (size_t o = p_index, i; o>0; --o)
    	{
            if (o==p_index) continue;
    		i = path[o-1];
            outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
    	}
    	outt<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

        // size_t last_w=PointerListLen - 1, last_r=PointerListLen-1;
        // for (size_t o = p_index, i; o>0; --o)
        // {
        //     if (o==p_index) continue;
        //     i = path[o-1];

        //     w = wd[i].index_of_W;
        //     r = wd[i].index_of_R;
        //     l = wd[i].len;

        //     if (w!=last_w || r!=last_r)
        //     {
        //         if (w!=last_w && r!=last_r)
        //         {

        //         }
        //         else if(w==last_w)
        //         {

        //         }
        //         else if(r==last_r)
        //         {

        //         }
        //     }
        //     outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        // }
    }
    
    uint32_t do_alignment(string& dna_f, size_t window_up, size_t window_down, string&  read, string&ss, ofstream& outt)
    {
        // for (size_t o = p_index, i; o>0; --o)
        // {
        //     if (o==p_index) continue;
        //     i = path[o-1];
        //     outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        // }
        // outt<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

        size_t last_w=window_up + PointerListLen - 1, last_r=PointerListLen-1;
        size_t index_ss=0;
        char snum[100];
        memset(snum, 0, sizeof(snum));
        for (size_t o = p_index, i, r, w, l; o>0; --o)
        {
            if (o==p_index) continue;
            i = path[o-1];

            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            l = wd[i].len;

            if (w!=last_w || r!=last_r)
            {
                if (w!=last_w && r!=last_r)
                {
                    get_alignment(dna_f, last_w - PointerListLen + 1, w-last_w, read, last_r - PointerListLen + 1, r-last_r, 0, 0, w-last_w, r-last_r, ss);
                }
                else if(w==last_w)
                {
                    sprintf(snum, "%lu", (unsigned long)(r-last_r));
                    ss.append(snum);
                    sprintf(snum, "I");
                    ss.append(snum);
                }
                else if(r==last_r)
                {
                    sprintf(snum, "%lu", (unsigned long)(w-last_w));
                    ss.append(snum);
                    sprintf(snum, "D");
                    ss.append(snum);
                }
            }
            last_w = w+l;
            last_r = r+l;
            sprintf(snum, "%lu", (unsigned long)l);
            ss.append(snum);
            sprintf(snum, "M");
            ss.append(snum);
            // outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        }

        outt << ss << endl;
    }
    void create_matrix()
    {
    	V = new size_t*[node_i];
        for (size_t i=0; i<node_i; ++i) V[i] = new size_t[node_i];
        uint32_t iw, ir, jw, jr;
    	size_t il;

        for (size_t i=0; i<node_i-1; ++i)
        {
            iw = wd[i].index_of_W;
            ir = wd[i].index_of_R;
            il = wd[i].len;
            for (size_t j=i+1; j<node_i; ++j)
            {
                V[i][j]=0;
                jw = wd[j].index_of_W;
                jr = wd[j].index_of_R;
                if (jw >= il + iw && jr >= il + ir && jr <= ir + il + t_wait)
                {
                    V[i][j] = 1;
                }
            }
        }
    }
};