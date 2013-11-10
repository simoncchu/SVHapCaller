#ifndef _CORE_PAIR_END_READS_H_
#define _CORE_PAIR_END_READS_H_
#include<vector>
#include<utility>
#include<string>
#include"bam_alignment_record.h"

/*
Description: Parse reads information from each alignment records in the sam file. 
*/
class ReadsParser
{

public:
	ReadsParser(){};
	ReadsParser(int len_read);
	ReadsParser(int len_read, int insert_size, int std_variation);
	ReadsParser(std::vector<BamAlignmentRecord*>* records,int len_read, int insert_size, int std_variation, int cnt_base_diff);

public:
	void setReadLenth(int readlenth);
	void setInsertSize(int insertsize);
	void setStdVariation(int stdvariation);

public:
	void sortAlignRecords();// sort alignment records according to pos field in each alignment 
	std::string getCigarByPos(int pos);//get the cigar of a read by the pos of the read, and the reads are sorted. 
	void cntPEReadsSpanDel(long brkpnt1, long brkpnt2, int& discon_cnt, int& con_cnt);  //get number of concordant and disconcordant
	                                                                                    //paired-end reads span a deletion. 
	void cntOverlapReadsOverDel(long brkpnt1, long brkpnt2, int& overlap_cnt);//get number of overlapped reads in a deletion
	void cntSplitReadsOverDel(long brkpnt1,long brkpnt2, int& split_cnt, int& unsplit_cnt);//Description: Count the number of split reads and unsplit read	                                                                                       
	void cntSplitMappedReads(long brkpnt1,long brkpnt2, int& split_cnt, int& unsplit_cnt);// Count the # of reads split mapped over the deletion.
	void cntMapUnMapReads(long brkpnt1,long brkpnt2, int& mapped_cnt, int& unmapped_cnt);//Count the number of mapped reads and unmapped reads
	                                                                                     //between [brkpnt1, brkpnt2]. 

public:
	bool isPEReadEncompassDel(int read_pos, int mate_pos, int lbrkpnt, int rbrkpnt);//whether a pair-end reads encompassing a deletion.
	bool isReadSplitDel(int read_pos,std::vector< std::pair<std::string,int> >& cigar, int lbrkpnt, int rbrkpnt,int slack=3);//whether a read split-mapped at a deletion breakpoint.
	bool inReadOverDelFullMap(int read_pos, int lbrkpnt, int rbrkpnt);//whether a read fully mapped betweeen two breakpoints.

public:
	int getReadLength();
	int getInsertSize();
	int getStdVariation();

private:
	int readlength; //length of the read in the sam file 
	std::vector<BamAlignmentRecord*>* vrecords; // alignment records of the sam file. 
	int insert_size; // insert size of the paired-end reads
	int std_variation;// standard variation of the paired-end reads 
	std::vector<int> vmappedreads;// save the mapped reads of a deletion

	int cnt_base_diff;//when reads split mapped,how many bases are different from the reference between[brkpnt1,brkpnt2]
};

#endif
