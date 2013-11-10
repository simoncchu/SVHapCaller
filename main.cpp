#include<iostream>
#include<string>
#include"bam_parse.h"

using namespace std;

int main(int argc, char* argv[])
{
	
	int chromid=0;
	string bam_name = argv[1];
	BamParse bp(bam_name);
	//for(int i=1;i<23;i++)
	//{
		bool b=bp.parseAlignment(0,0,0,249250621);
		if(b==false)
			cout<<"set region failed!!"<<endl;
		//chromid=i;
	//}

	cout<<bp.bam_aln_records.size()<<endl;
	//each chromosome by chromosome
	//for bam
	//also for vcf files 
	return 0;
}

