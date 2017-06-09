# CC=	g++
# SV_OBJS=SNP_vector.o
# LV_OBJS=LV_deep.o
# HF_OBJS=Haffman.o
# RUN_R=run_rHAT.o

# all: rHAT

# rHAT: $(SV_OBJS) $(HF_OBJS) $(LV_OBJS) $(RUN_R)
# 	$(CC) create_RHT.cpp -lz -O3 -o create_RHT 
# 	$(CC) -pthread $(RUN_R) $(HF_OBJS) $(SV_OBJS) $(LV_OBJS) ksw.c -O3 -lz -o run_rHAT ../htslib/libhts.a

# $(RUN_R):run_rHAT.cpp

# $(LV_OBJS):LV_deep.h

# $(SV_OBJS): SNP_vector.h 

# $(HF_OBJS):Haffman.h

# clean:
# 	$(RM) *.o test


all: rHAT

rHAT:
	g++ create_RHT.cpp -lz -O3 -o create_RHT 
	g++ -c -o Haffman.o Haffman.c
	g++ -c -o SNP_vector.o SNP_vector.c
	g++ -c -o LV_deep.o LV_deep.c
	g++ -pthread run_rHAT.cpp LV_deep.o SNP_vector.o Haffman.o ksw.c -O3 -lz -o run_rHAT ../htslib/libhts.a
	# g++ -pthread run_rHAT.cpp SNP_vector.o Haffman.o ksw.c -O3 -lz -o run_rHAT ../htslib/libhts.a

