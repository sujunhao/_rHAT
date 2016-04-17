all: rHAT

rHAT: 
	g++ run_rHAT.cpp ksw.c -o run_rHAT
	g++ create_RHT.cpp -o create_RHT

