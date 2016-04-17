#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <cstring>
#include <cmath>
#include <queue>
#include <unistd.h>
#include <cstring>
#include <stdint.h>
#include "ksw.h"

using namespace std;
size_t PointerListLen = 11;

const uint8_t seq_nt4_tablet[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
void prseq(char *str, int len, bool enter)
{
    for (int i=0;i<len;++i)  {
        cout<<str[i];
    }
    if (enter)
        cout<<endl;

}
inline void transIntoDec(uint8_t *transtr,char *str, int length)
{
    for (int i=0;i<length;++i) {
        transtr[i] = seq_nt4_tablet[str[i]];
        // printf("%c ", str[i]);
    }
}


int8_t  mat[25];
int8_t match=1, mism=-1, gape=2, gapo=-1;
const char correspondTable[] = "MIDNSHP=X";



uint32_t path_d[]={32, 27, 15, 10};
uint32_t path_r[]={31, 24, 11, 10};
uint32_t path_len[]={0, 11, 11, 0};
size_t p_index=4;

    double do_alignment(char* ref, size_t window_up, size_t window_down, char*  read, size_t read_l, string& thecigar, FILE* outt)
    {
        
        int k=0;
        for (int i=0; i<5; ++i)
        {
            for (int j=0; j<5; ++j)
            {
                if (i<4 && j<4) mat[k++] = i == j? match : mism;
                else mat[k++]=0;
            }
        }


        uint32_t *cigar;
        int n_cigar = 0;
        uint8_t  readqry[1024];
        uint8_t  refqry[1024];
        uint8_t *readqry_;
        uint8_t *refqry_;


        int read_len;
        int ref_len;
        int w_;

        uint32_t countM = 0;
        char     trans_cigar[50];
        int startPosCigar = 0;


        double score=0;

        size_t last_w=window_up + PointerListLen - 1, last_r=PointerListLen-1;
        size_t i, w, r, l;

        if (p_index >= 3)
        {
            i = path[p_index-2];
            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            l = wd[i].len;
            ref_len = w-last_w;
            read_len = r-last_r;
            cout << string(read, last_r - PointerListLen + 1, read_len) << " " << string(ref, last_w - PointerListLen + 1, ref_len) << endl;
            if (w!=last_w && r!=last_r)
            {
                transIntoDec(readqry,read + last_r - PointerListLen + 1,read_len);
                transIntoDec(refqry,ref + last_w - PointerListLen + 1,ref_len);
                readqry_ = readqry;
                refqry_ = refqry;
                // cout << r << " " << last_r << " " << ref_len << " " << read_len << endl;
                w_ = max(ref_len, read_len);
                n_cigar = 0;
                score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
                if (n_cigar) {
                    for (int z=0;z<n_cigar-1;++z) {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    }

                    if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
                        countM = cigar[n_cigar-1] >> 4;
                        //startPosCigar += sprintf(trans_cigar,"%u%c",countM,'M');
                        //sams[countSam].cigar.append(trans_cigar);
                    } else {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
                        thecigar.append(trans_cigar);
                    }

                } 
                free(cigar);
            }
            score+=(l)*match;
            cout << "---" << countM << endl;
            countM+=l;
            cout << "---" << countM << endl;
            last_w = w;
            last_r = r;
        }

        for (size_t o = p_index-2; o>1; --o)
        {
            if (o == p_index) continue;
            i = path[o-1];

            w = wd[i].index_of_W;
            r = wd[i].index_of_R;
            l = wd[i].len;

            

            ref_len = w-last_w-l;
            read_len = r-last_r-l;

            cout << last_r << " " << r << " " << last_w << " " << w << " " << l << endl;
            cout << string(read, last_r+1, read_len) << " " << string(ref, last_w+1, ref_len) << endl;

            if (w!=last_w || r!=last_r)
            {
                transIntoDec(readqry,read + last_r + 1,read_len);
                transIntoDec(refqry,ref + last_w + 1,ref_len);
                readqry_ = readqry;
                refqry_ = refqry;
                w_ = max(ref_len, read_len);
                n_cigar = 0;
                score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
                if (n_cigar) {
                    if (correspondTable[cigar[0]&0xf] == 'M') {
                        countM += (cigar[0] >> 4);
                        startPosCigar += sprintf(trans_cigar,"%uM",countM);
                        thecigar.append(trans_cigar);
                    }
                    else
                    {
                        startPosCigar += sprintf(trans_cigar,"%uM",countM);
                        thecigar.append(trans_cigar);
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    } 
                    countM = 0;
                    for (int z=1;z<n_cigar-1;++z) {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    }

                    if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
                        countM = cigar[n_cigar-1] >> 4;
                        //startPosCigar += sprintf(trans_cigar,"%u%c",countM,'M');
                        //sams[countSam].cigar.append(trans_cigar);
                    } else {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
                        thecigar.append(trans_cigar);
                    }

                } 
                free(cigar);
            }

            score+=(l)*match;
            cout << "---" << countM << endl;
            countM+=l;
            cout << "---" << countM << endl;
            last_w = w;
            last_r = r;
            // outt << wd[i].index_of_W << " " << wd[i].index_of_R << " " << wd[i].len << " " << read.substr(wd[i].index_of_R + 1 - PointerListLen, wd[i].len) << " " << dna_f.substr(wd[i].index_of_W + 1 - PointerListLen, wd[i].len)<< endl; 
        }

        printf("****\n");
        if (p_index >= 3)
        {
            w = window_down;
            r = read_l;
            ref_len = w-last_w-1;
            read_len = r-last_r-1;

            cout << last_r << " " << r << " " << last_w << " " << w << " " << l << endl;
            cout << string(read, last_r+1, read_len) << " " << string(ref, last_w+1, ref_len) << endl;

            if (w!=last_w || r!=last_r)
            {
                transIntoDec(readqry,read + last_r + 1,read_len);
                transIntoDec(refqry,ref + last_w + 1,ref_len);
                readqry_ = readqry;
                refqry_ = refqry;
                w_ = max(ref_len, read_len);
                n_cigar = 0;
                score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
                if (n_cigar) {
                    if (correspondTable[cigar[0]&0xf] == 'M') {
                        countM += (cigar[0] >> 4);
                        startPosCigar += sprintf(trans_cigar,"%uM",countM);
                        thecigar.append(trans_cigar);
                    }
                    else
                    {
                        startPosCigar += sprintf(trans_cigar,"%uM",countM);
                        thecigar.append(trans_cigar);
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    } 
                    countM = 0;
                    for (int z=1;z<n_cigar;++z) {
                        startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                        //sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
                        //++startPosCigar;
                        thecigar.append(trans_cigar);
                    }
                } 
                free(cigar);
            }
        }

        printf("%lf\n", score);
        return score;
    } 

    

int main()
{
    //initiate ksw parameters

    char read[]="ACCCCCCCCCCCAATTTTTTTTTTTAAAGG";
    char ref[]="AAAAACCCCCCCCCCCATTTTTTTTTTTGAGG";
    cout << string(ref) << endl << string(read) << endl;
    string thecigar;
    do_alignment(const_cast<char*>(ref), 0, strlen(ref), const_cast<char*>(read), strlen(read), thecigar, stdout);
    cout << thecigar << endl;
    
    return 0;

}