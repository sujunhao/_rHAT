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
void transIntoDec(uint8_t *transtr,char *str, int length)
{
    for (int i=0;i<length;++i) {
        transtr[i] = seq_nt4_tablet[str[i]];
    }
}


int8_t match=1, mism=-1, gape=2, gapo=-1;
int8_t  mat[25];
string cigar;

const   char    correspondTable[] = "MIDNSHP=X";

string thecigar;
int main()
{
    //initiate ksw parameters
    int k=0;
    for (int i=0; i<5; ++i)
    {
        for (int j=0; j<5; ++j)
        {
            if (i<4 && j<4) mat[k++] = i == j? match : mism;
            else mat[k++]=0;
        }
    }

    int n_cigar = 0;
    uint32_t *cigar;
    uint8_t  readqry[1024];
    uint8_t  refqry[1024];
    uint8_t *readqry_;
    uint8_t *refqry_;

    char read[]="ATTGC";
    char ref[]="ACTTGAC";

    int read_len = strlen(read);
    int ref_len = strlen(ref);
    int w;

    transIntoDec(readqry,read,read_len);
    transIntoDec(refqry,ref,ref_len);

    readqry_ = readqry;
    refqry_ = refqry;
    //prseq(readStartP,read_len,true);
    //prseq(refStartP,ref_len,true);

    w = read_len > ref_len ? read_len : ref_len;

    n_cigar = 0;

    double score = 0;
    score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w,&n_cigar,&cigar);

    //fprintf(stderr,"%d %d %d %d %d %d\n",order[i-1],order[i],read_len,ref_len,score, n_cigar);
    //first cigar[0] is M
    uint32_t countM = 0;
    char     trans_cigar[50];
    int startPosCigar = 0;
    if (n_cigar) {
        if (correspondTable[cigar[0]&0xf] == 'M') {
            countM += (cigar[0] >> 4);
            startPosCigar += sprintf(trans_cigar,"%uM",countM);
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
    //cigarbuflen -= usedCigarsize;
    //cigarStartP += usedCigarsize;
    
            

    cout << thecigar <<endl;
    return 0;

}