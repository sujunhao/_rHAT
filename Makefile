all: rHAT

rHAT: 
	g++ ksw.c -c -O3 -o ksw.o
	g++ run_rHAT.cpp -o run_rHAT
	g++ create_RHT.cpp -o create_RHT

