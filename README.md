# rHAT

[ref][1]


## usage:
```
make

//help
./create_RHT -h
./run_rHAT -h


//use E.coli.fa file to create RHT
./create_RHT -d E.coli.fa


//-d for dna file, -r for read file, -o for ouput file
./run_rHAT -d E.coli.fa -r dna_read -o out

```


-------

The rHAT program has two different parts.The first part is create_RHT. The second part is run_rHAT.

The first part create_RHT is aimed to use a DNA array to create an RHT table.
The second part run_rHAT is for using an RHT table to run rHAT process.

To run rHAT should set the pointerlistlen length and window length. They are store in PointerListLen(size_t [11]) and WindowListLen(size_t [2048])

##The main process of rHAT program 
```
1.0 create_RHT process:
1.1 first of all, the DNA string is store in a std::string named dna_w;
1.2 for each char in dna_w, convert into a two bit numnber(A->00, C->01, G->10, T->11) and then append to a variant named tmp;
1.3 take the last 2*PointerListLen bit nums in tmp to indicate pointer indices;
1.4 take the dna_w index to count the window index;
1.5 use pointer index and window index to bind together and store into a array named PW(uint64_t);
1.6 sort(std::sort) the PW array and then separate it into a pointer lists and a window lists as the RHT table;
1.7 write the RHT table to a file called out_RHT;


2.0 run_rHAT process:
2.1 get the DNA and read arrays form file and store in dna_f(std::string) and read(std::string);
2.2 get the RHT table from out_RHT(see 1.7);
2.3 for each read, split the center string of read. Does this process by counting the center read len(size_t theLen, if read < WindowlistLen/2 then theLen=read len), the len from read's start postion to read-center's start postion(size_t len_up_stream) and the len from read-center's end postion to read's end position(size_t len_down_stream);
2.4 in read center array do 1.1 to 1.2 process, use pointer index and RHT table loaded before to find the window index;
2.5 create an array to count the each window index's appearing time of 2.4 process;
2.6 use a minimum heap to get the k highest hitted windew index in the array;
2.8 use the whole read array to create a read l-mer RHT table;
2.8.1 use read to do 1.2 and 1.3 process;
2.8.2 take the read index as the read RHT window index, and do 1.5 process;
2.8.3 sort the read PW array;
2.9 for each k highest hitted window index in 2.6 do create-DAG and do_alignment process;

3.0 DAG process:
3.1 use window index to get DNA window array of DNA string and extend upstream len_upstream(see 2.3) char and extend downstream(see 2.3) char to create DNA window array;
3.2 use DNA window arrays(3.1) to run 1.2 and 1.3 process and use the pointer index and read's RHT table to find the read hitted index;
3.3 set the DNA window index(3.2) and read hitted inde(3.2) and match length (default PointerListLen) to create in a directed acyclic graph;
3.4 if two point p1 and p2 in DAG where p1.index_w = p2.index_w-1, p1.index_r=p2.index_r-1 and p1.len=p2.len+1, there is no need to store p2 and as a result the p2 can be deleted;
3.5 add the start(a, b, 0)(a is the window start postion in 3.2->1.2 pocess and b is the read start postion in 2.8.1->1.2 process) and end(windowlen, readlen, 0) point to DAG;
3.6 for each point in DAG, if two of them satisfy the condition "VRi + VLi <= VRj <= VRi + VLi + Twait and VSi + VLi <= VSj"(more detail see [ref][1] 2.4 3), add an edge(in a 2D array) to link that two pointers;
3.7 use the point and edge information to find a path from start point to end point(see 3.5) that having the maximum score(see [ref][1] 2.4 4, by a dynamic programming process).
3.8 for the path(3.7) in DAG, the point in the path is the window&read array that matched, then the path between two points is the gap array of window&read, so cut the gap array of window&read to do alignment process(4.0);

4.0 do alignment process:
4.1 get the first and last gaps in window&read to do semi-global alignment using ksw;
4.2 get the remaining gaps to do global alignment process;

5.0 print the result:
5.1 for a read, there have k-highest hitted window indexes, for each hitted index, there have a max alignment score from 3.0 and 4.0. print the max score in that alignment and print out corebonding alignment;
```

--------------------------

##Main data structure used in create_RHT.cpp

class RHT to create RHT table(see 1.5)
```
uint32_t *P;
uint32_t *W;
uint64_t *PW;
```
PW[] is a 64 bit number array. for each unit in PW[], the first 32 bit is the pointer index and the last 32 bit is the window index in DNA string;(see 1.6)
P[] is use to store first 32 bit info in PW[]
W[] is th last 32 bit array in PW[];



##Main data structure used in run_RHT.cpp

in the process of getting the k-highest window hitted index, the counting process is work by a mini heap:

```
priority_queue<WINDOW_CNT>
the WINDOW_CNT is 
typedef struct window_cnt {
    uint32_t index_of_W;
    size_t cnt;
    bool operator<(const window_cnt& w) const
    {
        return cnt > w.cnt;
    }
}WINDOW_CNT;
```
the WINDOW_CNT stores the window index and count number, and will be compared by count number.

the DAG pointer is some instant of 
```
typedef struct window_node {
    uint32_t index_of_W, index_of_R;
    size_t len;
}WINDOW_NODE;
```
which have window index, read index and length information


[1]: (http://bioinformatics.oxfordjournals.org/content/early/2016/01/09/bioinformatics.btv662.full#abstract-1)
