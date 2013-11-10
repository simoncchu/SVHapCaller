#ifndef _CORE_VCF_PARSE_H_
#define _CORE_VCF_PARSE_H_

#include<vector>
#include<utility>
#include<string>

class VCFParse
{
public:
	VCFParse(char* buffer, int len, int indvd_num);
	void parseBrkpntGenotype(); //parse breakpoints positions and genotypes in the vcf file. 

public:
	std::vector<std::pair<int,int> > brkpnts; // save the breakpoints in the vcf file.  
	std::vector<std::vector<std::pair<int,int> > > genotypes; // save the genotypes of each individual for each breakpoints. 
	                                                          //and genotypes[i] denotes the genotypes of the ith deletion
private:
	char* buffer; // pointer to the buffer store the vcf file 
	int len; // file length of the vcf file 
	int indvdu_num;// number of individuals 

};

#endif

