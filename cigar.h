#ifndef CIGAR_H_
#define CIGAR_H_


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
using namespace std;


#define MAPSTR "MIDNSHP=X"
#ifndef BAM_CIGAR_SHIFT
#define BAM_CIGAR_SHIFT 4u
#endif



typedef struct {
    uint32_t* seq;
    int32_t length;
} Cigar;



/*!     @function               Produce CIGAR 32-bit unsigned integer from CIGAR operation and CIGAR length
        @param  length          length of CIGAR
        @param  op_letter       CIGAR operation character ('M', 'I', etc)
        @return                 32-bit unsigned integer, representing encoded CIGAR operation and length
*/
uint32_t to_cigar_int (uint32_t length, char op_letter)
{
        switch (op_letter) {
                case 'M': /* alignment match (can be a sequence match or mismatch */
                default:
                        return length << BAM_CIGAR_SHIFT;
                case 'S': /* soft clipping (clipped sequences present in SEQ) */
                        return (length << BAM_CIGAR_SHIFT) | (4u);
                case 'D': /* deletion from the reference */
                        return (length << BAM_CIGAR_SHIFT) | (2u);
                case 'I': /* insertion to the reference */
                        return (length << BAM_CIGAR_SHIFT) | (1u);
                case 'H': /* hard clipping (clipped sequences NOT present in SEQ) */
                        return (length << BAM_CIGAR_SHIFT) | (5u);
                case 'N': /* skipped region from the reference */
                        return (length << BAM_CIGAR_SHIFT) | (3u);
                case 'P': /* padding (silent deletion from padded reference) */
                        return (length << BAM_CIGAR_SHIFT) | (6u);
                case '=': /* sequence match */
                        return (length << BAM_CIGAR_SHIFT) | (7u);
                case 'X': /* sequence mismatch */
                        return (length << BAM_CIGAR_SHIFT) | (8u);
        }
        return (uint32_t)-1; // This never happens
}

/*! @function       Extract CIGAR operation character from CIGAR 32-bit unsigned integer
    @param  cigar_int   32-bit unsigned integer, representing encoded CIGAR operation and length
    @return         CIGAR operation character ('M', 'I', etc)
*/
//char cigar_int_to_op (uint32_t cigar_int);
static inline char cigar_int_to_op(uint32_t cigar_int) 
{
    return (cigar_int & 0xfU) > 8 ? 'M': MAPSTR[cigar_int & 0xfU];
}

//copy cigar 2 to 1
static inline char copy_cigar(Cigar *cg1, Cigar *cg2)
{
    size_t len = cg2->length;
    cg1->length = len;
    for (size_t i=0; i<len; ++i) cg1->seq[i] = cg2->seq[i];
} 

/*! @function       Extract length of a CIGAR operation from CIGAR 32-bit unsigned integer
    @param  cigar_int   32-bit unsigned integer, representing encoded CIGAR operation and length
    @return         length of CIGAR operation
*/
//uint32_t cigar_int_to_len (uint32_t cigar_int);
static inline uint32_t cigar_int_to_len (uint32_t cigar_int)
{
    return cigar_int >> BAM_CIGAR_SHIFT;
}


void refine_cigar(Cigar *cg)
{
    uint32_t cnt = 0, tmpn;
    char c = cigar_int_to_op((cg->seq)[0]), tmpc;
    int32_t len = cg->length, nl = 0;
    for (int i=0; i<len; ++i)
    {
        tmpn = cigar_int_to_len((cg->seq)[i]);
        tmpc = cigar_int_to_op((cg->seq)[i]);
        if (tmpc == c)  cnt += tmpn;
        else 
        {
            (cg->seq)[nl++] = to_cigar_int(cnt, c);
            c = tmpc;
            cnt = tmpn;
        }
        // printf("%c", d_path[i]);
    }
    (cg->seq)[nl++] = to_cigar_int(cnt, c);
    (cg->length) = nl;
}


#endif