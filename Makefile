all: rHAT

rHAT: 
	g++ run_rHAT.cpp ksw.c -lz -o run_rHAT
	g++ create_RHT.cpp -lz -o create_RHT