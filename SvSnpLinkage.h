#ifndef _H_RESRAPETIS_
#define _H_RESRAPETIS_

#include<string>
#include<vector>
#include<utility>
#include"bam_alignment_record.h" 
#include"reads_parse.h"

class SvSnpLinkage
{
public:
	SvSnpLinkage();
	//SvSnpLinkage();
	SvSnpLinkage(int readlength, int insertsize, int deviation, std::string bamfile, std::string svfile, std::string snpfile);

public:
	void findLinkage();//call the overlapped reads.
	void getHitReads(std::vector<BamAlignmentRecord*>& bam_aln_records);//separately get the hit reads for SV and SNP respectively.

private:
	void getSVHitDiscdntReads(BamAlignmentRecord*& barecord);
	void getSVHitSplitReads(BamAlignmentRecord*& barecord);
	void getSVHitMBtwnReads(BamAlignmentRecord*& barecord);
	void getSNPHitReads(BamAlignmentRecord*& barecord);

	void parseSV(int chrom);//load all the sv position of specified chromosome into memory 
	void parseSNP(int chrom);//load all the snp position of specified chromosome into memory 

public:
	std::vector<int>vsv;//contain all the hit reads for sv
	std::vector<std::pair<int,int> > vsv_discordant; //save pair-end reads positions contain information of deletion, only for discordant pair.
	std::vector<int> vsv_split;//save split-mapped reads positions.
	std::vector<int> vsv_btwn;//save fully mapped reads positions that mapped between breakpoints.
	std::vector<int> vsnp;//save reads position contain snp information.

private:
	void combineAllHitReadsSV();//combine all the hit reads for sv

private:
	ReadsParser rp;
	
private:
	std::string bamfile;
	std::string vcf_sv_file;
	std::string vcf_snp_file;
	std::vector<std::pair<int,int> > vsv_pos;//to save the parsed positions of sv.
	std::vector<int> vsnp_pos;//to save the parsed positions of snp.
};

#endif


