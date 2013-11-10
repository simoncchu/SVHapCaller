
#include"vcf_parse.h"
#include <stdlib.h>
#include<iostream>

VCFParse::VCFParse(char* buffer, int len, int indvdu_num)
{
	this->buffer = buffer;
	this->len = len;
	this->indvdu_num = indvdu_num;
}

/*
Description: parse the breakpoints and genotype in the vcf file.
*/
void VCFParse::parseBrkpntGenotype()
{// Because the vcf file in this project doesn't follow the "VCF (Variant Call Format) version 4.0", so doesn't use VCFRecord class here. 
	
	std::string strbp1="";
	std::string strbp2="";
	int igntp1;
	int igntp2;
	std::string temp="";

	int ibp1=0;
	int ibp2=0;
	int cnt_tab=0; 
	//int cnt_deletion=0; // tell the current number of deletions 

	std::vector<std::pair<int,int> > vtemp; 
	vtemp.clear();

	for(int i=0;i<len;i++)
	{
		if(buffer[i]=='\n')
		{
			ibp1=atoi(strbp1.c_str());
			ibp2=atoi(strbp2.c_str());

			brkpnts.push_back(std::make_pair(ibp1,ibp2)); 
			cnt_tab=0;

			/* //depends whether there is a 'space' at end of each line 
			char c1 = temp.at(0);
			char c2 = temp.at(2);
			igntp1 = c1-'0';
			igntp2 = c2-'0'; 
			vtemp.push_back(std::make_pair(igntp1,igntp2));*/

			genotypes.push_back(vtemp);//each deletion a line, so move to another deletion 
			vtemp.clear();
		}
		else if(buffer[i]=='\t' || buffer[i]=='	' || buffer[i]==' ')
		{
			cnt_tab++;
			if(cnt_tab==3)
			{
				strbp1 = temp;
			}
			else if(cnt_tab==4)
			{
				strbp2 = temp;
			}
			else if(cnt_tab>4)
			{
				char c1 = temp.at(0);
				char c2 = temp.at(2);
				igntp1 = c1-'0';
				igntp2 = c2-'0'; 

				vtemp.push_back(std::make_pair(igntp1,igntp2));
			}
			temp="";
		}
		else
		{
			temp+=buffer[i]; 
		}
	}
	
}
