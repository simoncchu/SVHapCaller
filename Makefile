BAMTOOLS=../bamtools-master/include
BAMTOOLS_LD=../bamtools-master/lib

CC = g++
CFLAGS =  -O3 -Wall -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)

all: SVHapCaller

public_func.o :	public_func.cpp public_func.h
	$(CC) $(CFLAGS) -c public_func.cpp

bam_parse.o : bam_parse.cpp bam_parse.h
	$(CC) $(CFLAGS) -c bam_parse.cpp

reads_parse.o : reads_parse.cpp reads_parse.h
	$(CC) $(CFLAGS) -c reads_parse.cpp

vcf_parse.o : vcf_parse.cpp vcf_parse.h
	$(CC) $(CFLAGS) -c vcf_parse.cpp

SvSnpLinkage.o : SvSnpLinkage.cpp SvSnpLinkage.h
	$(CC) $(CFLAGS) -c SvSnpLinkage.cpp

main.o: main.cpp 
	$(CC) $(CFLAGS) -c main.cpp

SVHapCaller: public_func.o bam_parse.o reads_parse.o vcf_parse.o SvSnpLinkage.o main.o 
	$(CC) $(CFLAGS) -o SVHapCaller_1 public_func.o bam_parse.o reads_parse.o vcf_parse.o SvSnpLinkage.o main.o \
	-lbamtools -lz -lm
