#ifndef BAM_PARSE_H_
#define BAM_PARSE_H_

#include <api/BamReader.h>
#include "api/BamAux.h"
using namespace BamTools;

#include"bam_alignment_record.h"
#include"bam_header_record.h"

#include<vector>
#include<string>
using namespace std;

class BamParse
{

public:
	BamParse(string filename);
	BamParse(vector<string>*& listname);
	~BamParse();

public:
	void parseHeader();
	bool parseAlignment(int chrom_begin=-1, int start=-1, int chrom_end=-1, int end=-1);

	void sortAlignmentByChrom();//sort all the alignment, in the priority of Chrom first, then read Name.

private:
	void setAlignmentRecord(BamAlignment& al,BamAlignmentRecord*& bar);
	void parseCigar(vector<CigarOp>& cgdata, BamAlignmentRecord*& bar);//parse cigar 
	bool loadIndex(BamReader& br);
public:	
	std::vector<BamAlignmentRecord*> bam_aln_records;
	std::vector<BamHeaderRecord*> bam_head_records;
	int min_read_length;
	int max_read_length;
	int insert_size;
	double std_variation;

private:
	string filename;//bam file name
	string filename_index;//bam index file name 
	vector<string>* listfname; 

};

#endif

