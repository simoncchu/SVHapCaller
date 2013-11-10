#include"bam_parse.h"

#include<iostream>
#include <algorithm>

using namespace std;

BamParse::BamParse(string filename)
{
	this->filename=filename;
	this->filename_index=filename+".bai";
	min_read_length=500;
	max_read_length=0;
	insert_size=0;
	std_variation=0.0;
}

BamParse::BamParse(vector<string>*& listname)
{
	this->listfname=listname;
	min_read_length=500;
	max_read_length=0;
	insert_size=0;
	std_variation=0.0;
}

BamParse::~BamParse()
{
	int size=0;

	size=bam_head_records.size();
	for(int i=0;i<size;i++)
		delete bam_head_records[i];

	size=bam_aln_records.size();
	for(int i=0;i<size;i++)
		delete bam_aln_records[i];

}

void BamParse::parseCigar(vector<CigarOp>& cgdata, BamAlignmentRecord*& bar)//parse cigar
{
	pair<string,int> ptemp;
	for(int i=0;i<cgdata.size();i++)
	{
		ptemp.first=cgdata[i].Type;
		
		ptemp.second=cgdata[i].Length;

		bar->cigar.push_back(ptemp);
	}
}

void BamParse::parseHeader()
{
	
}

bool BamParse::loadIndex(BamReader& reader)
{
	bool bsignal = reader.OpenIndex(this->filename_index);
	if(bsignal==false) 
	{
		cerr << "Index file not found, now create it!!!"<<endl;
		bool bci=reader.CreateIndex();
		if(bci==false)
		{
			cerr << "Index file cannot be created!!!"<<endl;
			return false;
		}
		else
		{
			bool blct=reader.LocateIndex();
			if(blct==false)
			{
				cerr << "Index file cannot be located!!!"<<endl;
				return false;
			}
		}

	}
	return true;
}

/*
Description:
	Load all the bam into memory at one time if no parameters set, otherwise load the needed part of the bam.
	Save the parsed info into vector. 
*/
bool BamParse::parseAlignment(int chrom1, int chrom1_begin, int chrom2, int chrom2_end)
{
	BamReader reader;
    if ( !reader.Open(filename) ) {
        cerr << "Bamtools ERROR: could not open input BAM file: " << filename << endl;
        return false;
    }
		
	//check whether need to set a region.
	if(chrom1>-1 && chrom1_begin>-1 && chrom2>-1 && chrom2_end>-1)
	{
		this->loadIndex(reader);
		BamRegion br(chrom1,chrom1_begin,chrom2,chrom2_end);
		bool is_set=reader.SetRegion(br);
		if(is_set==false)
		{
			return false;//cannot set the region.
		}
	}

	//process input data
    BamAlignment al;   
	while ( reader.GetNextAlignment(al) )
	{
		if(al.Position<0) continue;

		BamAlignmentRecord* bar=new BamAlignmentRecord();
		setAlignmentRecord(al,bar);
		bam_aln_records.push_back(bar);
	}

	reader.Close();
	return true;
}



void BamParse::setAlignmentRecord(BamAlignment& al,BamAlignmentRecord*& bar)
{
	// mandatory fields 
	//parse QNAME 
	bar->qName=al.Name;
	//parse FLAG
	bar->flag=al.AlignmentFlag;
	//parse REF SEQ ID 
	bar->rID=al.RefID;// start from 0, for human_g1k_v37, 0 represent chr1, and 23 represent chrx, and so on.
	//parse 1-based mappping position, POS
	bar->pos=al.Position;
	//parse MAPQ
	bar->mapQ=al.MapQuality;
	//parse bin
	bar->bin=al.Bin;
	//parse CIGAR
	parseCigar(al.CigarData, bar);
	//parse RNEXT ID
	bar->rNextID=al.MateRefID;
	//parse PNEXT, 1-based
	bar->pNext=al.MatePosition;
	//parse TLEN
	bar->tLen=al.Length;
	//parse SEQ
	bar->seq="*";
	//bar->seq=al.QueryBases;
	//parse QUAL
	bar->qual=al.Qualities;
		
	/*option fields*/
	//string stroptfileds="";
	//for(int i=10;i<vemp.size();i++)
	//	stroptfileds+=vemp[i];
}

bool cmp(BamAlignmentRecord* bar1, BamAlignmentRecord* bar2)
{
	if(bar1->rID==bar2->rID)
	{
		if(bar1->pos==bar2->pos)
			return bar1->qName < bar2->qName;
		else 
			return bar1->pos < bar2->pos;
	}
	else 
		return bar1->rID < bar2->rID;
}

/*
Description:
Sort all the alignment, in the priority of Chrom first, then read Name.
*/
void BamParse::sortAlignmentByChrom()
{
	sort(bam_aln_records.begin(),bam_aln_records.end(),cmp);
}

