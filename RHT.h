#include <iostream>
#include <string>
#include <fstream>
#include <bitset>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <stdint.h>
using std::string;

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


//transformat form string to bitset
uint32_t to_bit(string s)
{
    uint32_t bitnum = 0, shift, bit_len = PointerListLen;
    for (size_t i = 0; i < bit_len; i++) 
    {
        shift = 2*bit_len-2-2*(i%(bit_len));
        switch (s[i])
        {
            case 'A':
                bitnum |= BASE_A << shift;
                break;
            case 'C':
                bitnum |= BASE_C << shift;
                break;
            case 'G':
                bitnum |= BASE_G << shift;
                break;
            case 'T':
                bitnum |= BASE_T << shift;
                break;
            default:
                std::cout << s << " invalid DNA base\n";
                return 0;

        }
    }
    return bitnum;
}

//transformat form bitset to string
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


//keep the RHT
class dna_bitset
{
    public:
        dna_bitset (size_t len) 
        {
            bit_len = len;
            m_len = 1 << (2 * bit_len);
            m_data = new std::vector<std::pair<uint32_t, uint32_t> >[m_len];
        }

        ~dna_bitset ()
        {
            delete [] m_data;
        }

        size_t get_count(uint32_t n)
        {
            return m_data[n].size();
        }

        uint32_t get_window(uint32_t n, size_t v)
        {
            return m_data[n][v-1].first;
        }

        void link_string(string s, uint32_t dna_ref_index)
        {
            uint32_t n = to_bit(s), v = m_data[n].size();
            if (v > 0)
            {
                if (m_data[n][v - 1].first == dna_ref_index)
                    ++m_data[n][v - 1].second;
                else if (v > 1 && m_data[n][v - 2].first == dna_ref_index)
                    ++m_data[n][v - 2].second;
                else if (m_data[n][v-1].first > dna_ref_index)
                    m_data[n].insert(m_data[n].end() - 1, std::make_pair(dna_ref_index, 1));
                else
                    m_data[n].push_back(std::make_pair(dna_ref_index, 1));
            }
            else
                m_data[n].push_back(std::make_pair(dna_ref_index, 1));            
            //std::cout << s << " " << n << " " << v << " " << dna_ref_index << std::endl;
        }

        void link_string(uint32_t n, uint32_t dna_ref_index)
        {
            uint32_t v = m_data[n].size();
            if (v > 0)
            {
                if (m_data[n][v - 1].first == dna_ref_index)
                    ++m_data[n][v - 1].second;
                else if (v > 1 && m_data[n][v - 2].first == dna_ref_index)
                    ++m_data[n][v - 2].second;
                else if (m_data[n][v-1].first > dna_ref_index)
                    m_data[n].insert(m_data[n].end() - 1, std::make_pair(dna_ref_index, 1));
                else
                    m_data[n].push_back(std::make_pair(dna_ref_index, 1));
            }
            else
                m_data[n].push_back(std::make_pair(dna_ref_index, 1));            
            // std::cout << s << " " << n << " " << v << " " << dna_ref_index << std::endl;
        }



        //write hash table to fstream
        void write_hash_out(std::ofstream &oo)
        {
            string ts;
            for (uint32_t i=0; i < m_len; ++i) 
            {
                ts = to_string(i);
                if (m_data[i].size() > 0)
                {
                    //int string
                    oo << ts;
                    // fprintf(oo, "%s", ts.c_str());

                    //in bitnum style
                    // oo << to_bit(ts);

                    for (size_t j=0; j < m_data[i].size(); ++j)
                    {
                        //write the hash table to file, the pointer string, the windowList index and occure time
                        // if (m_data[i][j].second  > 1)
                        // fprintf(oo, " %d %d", m_data[i][j].first, m_data[i][j].second);
                        // oo << " " << m_data[i][j].first <<  " " << m_data[i][j].second;

                        //the pointer string, the windowList index
                        // fprintf(oo, " %d", m_data[i][j].first);
                        oo << " " << m_data[i][j].first;
                    }
                    // fprintf(oo, "\n");
                    oo << std::endl;  
               }
            }
        }

        void write_hash_out(FILE *oo)
        {
            string ts;
            for (uint32_t i=0; i < m_len; ++i) 
            {
                ts = to_string(i);
                if (m_data[i].size() > 0)
                {
                    //int string
                    // oo << ts;
                    fprintf(oo, "%s", ts.c_str());

                    //in bitnum style
                    // oo << to_bit(ts);

                    for (size_t j=0; j < m_data[i].size(); ++j)
                    {
                        //write the hash table to file, the pointer string, the windowList index and occure time
                        // if (m_data[i][j].second  > 1)
                        // fprintf(oo, " %d %d", m_data[i][j].first, m_data[i][j].second);
                        // oo << " " << m_data[i][j].first <<  " " << m_data[i][j].second;

                        //the pointer string, the windowList index
                        fprintf(oo, " %d", m_data[i][j].first);
                        // oo << " " << m_data[i][j].first;
                    }
                    fprintf(oo, "\n");
                    // oo << std::endl;  
               }
            }
        }


        void read_hash_in(std::ifstream &ii)
        {
            string ss, line;
            uint32_t n, t;
            while (getline(ii, line))
            {
                std::istringstream ll(line);
                ll >> ss;
                n = to_bit(ss);
                while (ll >> t)
                {
                    m_data[n].push_back(std::make_pair(t, 1));
                }
            }
        }
    private:
        std::vector<std::pair<uint32_t, uint32_t> >* m_data;
        //the table list size  
        uint32_t m_len;
        //the table item size 11
        size_t bit_len;
};


class RHT
{
    private:
        size_t pw_len, p_len, w_len;
        uint64_t tmp;
    public:
        uint32_t index, w_index;
        uint64_t *PW;
        uint32_t *P, *W;

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

        void clear()
        {
            memset(PW, sizeof(PW), 0);
            // memset(P, sizeof(P), 0);
            index = 0;
        }

        void link_string(uint32_t p_index, uint32_t w_index)
        {
            // P[p_index] = index;
            tmp = p_index;
            tmp = tmp << 32;
            tmp = tmp | w_index;
            // std::cout << std::bitset<64>(tmp) << " " << std::bitset<32>(p_index) << " " << std::bitset<32>(w_index) << std::endl;
            PW[index++] = tmp;
        }

        uint32_t search(uint32_t p_index)
        {
            //binary search
            uint32_t l=0, r=index-1, m;
            while (l < r)
            {
                m = (l+r)/2;
                if (p_index <= (PW[m]>>32)) r = m;
                else l = m + 1;
            }
            if ((PW[r]>>32) == p_index) return r+1;
            return 0;
        }

        void sort_pw()
        {
            std::sort(PW, PW+index);
            // uint32_t p_tmp;
            // for (size_t i = 0; i < index; ++i)
            // {
            //     p_tmp = PW[i]>>32;
            //     if (!P[p_tmp])
            //     P[p_tmp] = i+1;
            // }
        }
        void create_p_w()
        {
            std::sort(PW, PW+index);
            memset(P, 0, sizeof(P));
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
            // std::cout << w_index << std::endl;
            fwrite(W, sizeof(uint32_t), w_index, out);
        }

        void read_hash(FILE *in)
        {
            fread(P, sizeof(uint32_t), p_len, in);
            fread(&w_index, sizeof(uint32_t), 1, in);
            // std::cout << w_index << std::endl;
            fread(W, sizeof(uint32_t), w_index, in);
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
};