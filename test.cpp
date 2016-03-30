#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <bitset>
#include <queue>
#include <map>
#include <unistd.h>
#include <stdint.h>
// #include "global_alignment.h"
using namespace std;

size_t PointerListLen = 11;
size_t WindowListLen = 2048;

uint32_t get_c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    3};

uint32_t MASK = 0xffffffff;
#define BASE_MASK 0x3

enum
{
    BASE_A = 0x0,
    BASE_C = 0x1,
    BASE_G = 0x2,
    BASE_T = 0x3,
};

string to_string(uint32_t bb)
{
    uint32_t bitnum = bb, shift, base, mask, bit_len = PointerListLen;
    string s(bit_len, 0);
    for (size_t i = 0; i < bit_len; i++)
    {
        shift = 2*bit_len-2-2*(i%(bit_len));
        mask = BASE_MASK << shift;
        base = (bitnum & mask) >> shift;

        switch (base)
        {
            case BASE_A:
                s[i] = 'A';
                break;
            case BASE_C:
                s[i] = 'C';
                break;
            case BASE_G:
                s[i] = 'G';
                break;
            case BASE_T:
                s[i] = 'T';
                break;
            default:
                std::cout << "invalid bitnum\n";
                return "";
        }
    }
    return s;
}

class RHT
{
    private:
        size_t pw_len, p_len, w_len;
        uint32_t *P, *W;
        uint64_t *PW;
        uint32_t index, w_index;
        uint64_t tmp;
    public:
        RHT(size_t k)
        {
            p_len = (size_t)1 << (PointerListLen << 1) + 1;
            w_len = k << 1;
            pw_len = w_len;

            P = new uint32_t[p_len];
            W = new uint32_t[w_len];
            PW = new uint64_t[pw_len];

            index = 0;
            w_index = 0;
        }
        ~RHT()
        {
            delete [] P;
            delete [] W;
            delete [] PW;
        }

        void link_string(uint32_t p_index, uint32_t w_index)
        {
            tmp = p_index;
            tmp = tmp << 32;
            tmp = tmp | w_index;
            // std::cout << std::bitset<64>(tmp) << " " << std::bitset<32>(p_index) << " " << std::bitset<32>(w_index) << std::endl;
            PW[index++] = tmp;
        }

        void create_p_w()
        {
            std::sort(PW, PW+index);
            // memset(P, 0, sizeof(P));
            uint32_t p_tmp, last;
            w_index=0;
            p_tmp = PW[0] >> 32;
            last = p_tmp;
            P[p_tmp] = w_index + 1;
            W[w_index] = PW[0];
            ++w_index;
            // std::cout << std::bitset<64>(PW[0]) << " " << std::bitset<32>(p_tmp) << " " << std::bitset<32>(W[w_index])  << std::endl;
            for (size_t i = 1; i<index; ++i)
            {
                if (PW[i]==PW[i-1]) continue;
                p_tmp = PW[i] >> 32;
                if (p_tmp == last) 
                {
                    W[w_index++] = PW[i];
                    continue;
                }
                
                P[last+1] = w_index+1;
                P[p_tmp] = w_index+1;
                last = p_tmp;
                W[w_index++] = PW[i];
                
            }
            P[last+1] = w_index+1;
        }

        void write_hash(FILE *out)
        {
            // fprintf(out, "%lu|", (unsigned long)w_index);
            fwrite(P, sizeof(uint32_t), p_len, out);
            fwrite(&w_index, sizeof(uint32_t), 1, out);
            std::cout << w_index << std::endl;
            fwrite(W, sizeof(uint32_t), w_index, out);
        }

        void write_hash_test(std::ofstream &out)
        {
            for (size_t i = 0; i<p_len-1; ++i)
            {
                if (P[i] > 0 && P[i+1] > 0 && P[i+1] > P[i])
                {
                    // out << std::bitset<32>(i);
                    out << to_string((uint32_t)i);

                    for (size_t j = P[i]-1; j < P[i+1]-1; ++j)
                    {
                        out << " " << W[j];
                    }
                    out << "\n";
                }
            }
        }
        void read_hash(FILE *in)
        {
            // unsigned long d;
            // fscanf(in, "%lu|", &d);
            // w_index = d;
            fread(P, sizeof(uint32_t), p_len, in);
            fread(&w_index, sizeof(uint32_t), 1, in);
            std::cout << w_index << std::endl;
            fread(W, sizeof(uint32_t), w_index, in);
        }
};

int main()
{

	RHT rht(4938920);
	FILE *pout;
    pout = fopen("out_R", "rb");
    ofstream out;
    out.open("outtt");

	rht.read_hash(pout);
	rht.write_hash_test(out);


	out.close();
	fclose(pout);

	return 0;
}