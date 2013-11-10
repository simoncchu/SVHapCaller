#include"SvSnpLinkage.h"
#include<fstream>
#include<iostream>
#include<algorithm>

using namespace std;

//const int INSERT_SIZE=500;
//const int STDVAR=50;
//const int READLEN=100;

const int FLAG_UNMAP=4;
const int FLAG_PAIRMAP=3;


SvSnpLinkage::SvSnpLinkage(int readlength, int insertsize, int deviation, std::string bamfile, std::string svfile, std::string snpfile)
{
	rp.setReadLenth(readlength);
	rp.setInsertSize(insertsize);
	rp.setStdVariation(deviation);

	this->bamfile=bamfile;
	this->vcf_sv_file=svfile;
	this->vcf_snp_file=snpfile;
}

/*
Description:
	Call out the useful reads, by analyzing the two vectors.
*/
void SvSnpLinkage::findLinkage()
{
	combineAllHitReadsSV();
	int size_sv=vsv.size();
	int size_snp=vsnp.size();
	int i=0,j=0;
	int cnt=0;
	while(i<size_sv && j<size_snp)
	{
		if(vsnp[j]==vsv[i])
			cnt++;
		else if(vsnp[j]<vsv[i])
			j++;
		else 
			i++;
	}
	cout<<cnt<<endl;
}

/*
Description: 
	Get all the hit reads.
Note:
	For the bam files, need to check the "flag" part, whether correct. 
*/
void SvSnpLinkage::getHitReads(std::vector<BamAlignmentRecord*>& bam_aln_records)//get the hit reads for SV respectively.
{
	int sv_size=vsv_pos.size();
	int snp_size=vsnp_pos.size();
	int ipsv=0, ipsnp=0;
	int isv_start = vsv_pos[0].first;
	int isv_end = vsv_pos[0].second;
	int isnp_pos = vsnp_pos[0];

	//check "Flag" of bam/sam file, and decide use which function. 
	int size=bam_aln_records.size();
	for(int i=0;i<size;i++)
	{
		int iread_pos=bam_aln_records[i]->pos;
		int imate_pos=bam_aln_records[i]->pNext;
		if(imate_pos<=iread_pos) continue;

		//find a proper deletion position.
		while(isv_start<iread_pos)
		{
			ipsv++;
			isv_start = vsv_pos[ipsv].first;
		}
		isv_end=vsv_pos[ipsv].second;

		//find a proper snp position
		while(isnp_pos<iread_pos)
		{
			ipsnp++;
			isnp_pos = vsnp_pos[ipsnp];
		}

		/*---For Deletion---------------------------------------------------------*/
		int flag=bam_aln_records[i]->flag;
		int iunmap=flag & FLAG_UNMAP;
		int ipairmap=flag & FLAG_PAIRMAP;
		if(iunmap==FLAG_UNMAP)
		{//read is unmapped
			//check whether split at breakpnt
			if(rp.isReadSplitDel(iread_pos,bam_aln_records[i]->cigar,isv_start,isv_end)==true)
			{
				getSVHitSplitReads(bam_aln_records[i]);
			}
		}
		else if(iunmap==FLAG_PAIRMAP)
		{//both reads of a pair if fully mapped
			//check whether encompassing and discordant.
			if(rp.isPEReadEncompassDel(iread_pos,imate_pos,isv_start,isv_end)==true)
			{
				getSVHitDiscdntReads(bam_aln_records[i]);
			}
		}
		else
		{//read is mapped, but mate is unmapped
			//check whether fully mapped between brkpnts.
			if(rp.inReadOverDelFullMap(iread_pos,isv_start,isv_end)==true)
			{
				getSVHitMBtwnReads(bam_aln_records[i]);
			}
		}
		
		/*-----For SNP----------------------------------------------------------------*/
		if(isnp_pos>=iread_pos && isnp_pos<=(iread_pos + this->rp.getReadLength()))
			getSNPHitReads(bam_aln_records[i]);

	}//end of for.

}
/*
Description:
	Get the reads positions for split reads.
*/
void SvSnpLinkage::getSVHitDiscdntReads(BamAlignmentRecord*& barecord)
{
	vsv_discordant.push_back(std::make_pair(barecord->pos,barecord->pNext));
}

void SvSnpLinkage::getSVHitSplitReads(BamAlignmentRecord*& barecord)
{
	vsv_split.push_back(barecord->pos);
}

void SvSnpLinkage::getSVHitMBtwnReads(BamAlignmentRecord*& barecord)
{
	vsv_btwn.push_back(barecord->pos);
}

void SvSnpLinkage::getSNPHitReads(BamAlignmentRecord*& barecord)//get the hit reads for  SNP
{
	vsnp.push_back(barecord->pos);
}

/*
Description:
	Parse the needed SV according to given chrom. And load into memory. 
*/
void SvSnpLinkage::parseSV(int spcfy_chrom)
{
	ifstream fsv;
	fsv.open(this->vcf_sv_file.c_str(),ifstream::in);
	int chrom,sv_start, sv_end;

	
	while(fsv>>chrom>>sv_start>>sv_end)
	{
		if(chrom>spcfy_chrom)
			break;
		else if(chrom==spcfy_chrom)
		{
			vsv_pos.push_back(std::make_pair(sv_start,sv_end));
		}
		//discard the genotype information
		std::string strdiscard;
		getline(fsv,strdiscard);
	}

	fsv.close();
}

/*
Description:
	Parse the needed SNP according to given chrom. And load into memory.
*/
void SvSnpLinkage::parseSNP(int spcfy_chrom)
{
	ifstream fsnp;
	fsnp.open(this->vcf_snp_file.c_str(),ifstream::in);
	int chrom,snp_pos;
	std::string name;
	
	while(fsnp>>chrom>>name>>snp_pos)
	{
		if(chrom>spcfy_chrom)
			break;
		else if(chrom==spcfy_chrom)
		{
			vsnp_pos.push_back(snp_pos);
		}
		//discard the genotype information
		std::string strdiscard;
		getline(fsnp,strdiscard);
	}

	fsnp.close();
}


void SvSnpLinkage::combineAllHitReadsSV()//combine all the hit reads for sv
{
	int size=vsv_discordant.size();
	for(int i=0;i<size;i++)
	{
		vsv.push_back(vsv_discordant[i].first);
	}

	size=vsv_split.size();
	for(int i=0;i<size;i++)
	{
		vsv.push_back(vsv_split[i]);
	}

	size=vsv_btwn.size();
	for(int i=0;i<size;i++)
	{
		vsv.push_back(vsv_btwn[i]);
	}
	sort(vsv.begin(),vsv.end());
}