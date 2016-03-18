all: rHAT

rHAT: 
	g++ run_rHAT.cpp -o run_rHAT
	g++ create_RHT.cpp -o create_RHT

