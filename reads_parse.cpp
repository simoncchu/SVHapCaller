#include<iostream>
#include<cstdio>
#include<algorithm>
#include"reads_parse.h"
#include"bam_alignment_record.h"
#include"public_func.h"

ReadsParser::ReadsParser(int len_read)
{
	this->readlength=len_read;
}

ReadsParser::ReadsParser(int len_read, int insert_size, int std_variation)
{
	this->readlength=len_read;
	this->insert_size=insert_size;
	this->std_variation=std_variation;
}

ReadsParser::ReadsParser(std::vector<BamAlignmentRecord*>* records,int len_read, int insert_size, int std_variation, int cnt_base_diff)
{
	this->vrecords=records;
	this->readlength= len_read;
	this->insert_size=insert_size;
	this->std_variation=std_variation;
	this->cnt_base_diff=cnt_base_diff;
}

void ReadsParser::setReadLenth(int readlength)
{
	this->readlength=readlength;
}

void ReadsParser::setInsertSize(int insertsize)
{
	this->insert_size=insertsize;
}

void ReadsParser::setStdVariation(int stdvariation)
{
	this->std_variation=stdvariation;
}

/*
Description: Parse sam/bam file and count the paired-end reads span the deletion, 
			 with given breakpoints, insert size, standard variation.
Return value: discon_cnt, save the number of discordant paired-end reads span the deletion. 
			  con_cnt, save the number of concordant paired-end reads span the deletion.
*/
void ReadsParser::cntPEReadsSpanDel(long brkpnt1,long brkpnt2, int& discon_cnt, int &con_cnt)
{
	discon_cnt=0;
	con_cnt=0;

	int len_deletion = brkpnt2 - brkpnt1;
	int leftmost_position = brkpnt1 - (insert_size + 3*std_variation);
	if(leftmost_position<0) leftmost_position=0;
	int rightmost_position=brkpnt2 + (insert_size + 3*std_variation);

	int size=vrecords->size(); 


	for(int i=0;i<size;i++) 
	{
		int pos=vrecords->at(i)->pos;

		if(pos>brkpnt2) break; // prerequirment is the sam file is sorted.
		
		int mate_pos = vrecords->at(i)->pNext;
		
		if(mate_pos>rightmost_position)
		{
			continue;
		}
		else if( pos < leftmost_position)
		{
			continue;
		}
		else if(pos > brkpnt1)
		{
			continue;
		}
		else if(mate_pos < brkpnt2)
		{
			continue;
		}
		else if(mate_pos < pos)
		{	
			continue;
		}
		else
		{
			
			std::string stemp = vrecords->at(i)->cigar.at(0).first; 
			int itemp = vrecords->at(i)->cigar.at(0).second;

			if(stemp!="M" || itemp!=readlength) // to check whether it is perfect match with the reference
			{
				//std::cout<<"6 "<<readlength<<" ; "<<itemp<<stemp<<" "<<pos<<" "<<mate_pos<<std::endl; /////////////////////////////////////////////////////////
				continue;
			}
			else
			{
				//check whether the mate-read is match perfectly 
				std::string cigar_mate = ReadsParser::getCigarByPos(mate_pos);
				std::string tmplt = PubFuncs::cvtInt2Str(readlength) + "M"; 
				
				//std::cout<<"mate cigar:"<<cigar_mate<<" temp cigar"<<tmplt<<std::endl;////////////////////////////////////////////////////////////////////////////////
				
				if(cigar_mate != tmplt) //mate read doesn't match well 
				{
					//std::cout<<"9 mate cigar"<<cigar_mate<<std::endl; //////////////////////////////////////////////////////
					continue;
				}
				int discord_distant = insert_size + 3*std_variation;
				if((mate_pos-pos) > discord_distant)
				{
					//std::cout<<"7 "<<pos<<" "<<mate_pos<<std::endl; //////////////////////////////////////////////////////
					discon_cnt++;
				}	
				else
				{
					//std::cout<<"8 "<<pos<<" "<<mate_pos<<std::endl; //////////////////////////////////////////////////////
					con_cnt++;
				}

			}

		}
	}
	

	return;
}

//void ReadsParser::cntPEReadsSpanDel(long brkpnt1,long brkpnt2, int& discon_cnt, int &con_cnt)
//{
//	
//	discon_cnt=0;
//	con_cnt=0;
//
//	//int len_deletion = brkpnt2 - brkpnt1;
//	//int leftmost_position = brkpnt1 - (insert_size + 3*std_variation + len_deletion); 
//	//if(leftmost_position<0) leftmost_position=0; 
//
//	int size=vrecords->size(); 
//	for(int i=0;i<size;i++) 
//	{
//		int pos=vrecords->at(i)->pos;
//		
//		if(pos>brkpnt1) break; // prerequirment is the sam file is sorted. 
//		
//		int mate_pos = vrecords->at(i)->pNext; 
//		
//		if( ((pos+readlength) < brkpnt1) && (pos< mate_pos) && (mate_pos>=brkpnt2)) 
//		{
//			std::string stemp = vrecords->at(i)->cigar.at(0).first; 
//			int itemp = vrecords->at(i)->cigar.at(0).second;
//
//			if((stemp=="M") && (itemp==readlength)) // to check whether it is perfect match with the reference 
//			{
//				//check whether the mate-read is match perfectly    
//				std::string cigar_mate = ReadsParser::getCigarByPos(mate_pos);
//				std::string tmplt = PubFuncs::cvtInt2Str(readlength) + "M"; 
//				if(cigar_mate != tmplt) //mate read doesn't match well
//					continue;
//
//				int discord_distant = insert_size + 3*std_variation;
//				if((mate_pos-pos) > discord_distant)
//				{
//					std::cout<<brkpnt1<<" "<<brkpnt2<<" "<<pos<<" "<<mate_pos<<std::endl; 
//					discon_cnt++;
//				}	
//				else
//				{
//					con_cnt++;
//				}
//
//			}
//		}
//	}
//	
//	return;
//}

/*
Description: Count the number of split-mapped reads. 
			And in this function, we only consider part of the read split at the breakpoint, but doesn't consider the other part.
*/
void ReadsParser::cntSplitReadsOverDel(long brkpnt1,long brkpnt2, int& split_cnt, int& unsplit_cnt) 
{
	split_cnt=0;
	unsplit_cnt=0;

	int vsize = vrecords->size();
	std::string temp_cigar;
	int pos;
	
	for(int i=0;i<vsize;i++)
	{
		pos=vrecords->at(i)->pos;
		if((pos > (brkpnt1-readlength) && pos < brkpnt1) || (((pos+readlength)>brkpnt2) && (pos<brkpnt2) ))
		{			
			int cg_size = vrecords->at(i)->cigar.size();
			if(cg_size==0) 
			{
				continue;
			}

			int itemp1=vrecords->at(i)->cigar.at(0).second;
			std::string strtemp1=vrecords->at(i)->cigar.at(0).first; 

			if(strtemp1=="*") 
			{
				continue;
			}

			if(cg_size==1) //totally mapped or totally unmapped.  
			{
				if(strtemp1!="M")
				{
					split_cnt++;
				}
				else 
				{
					unsplit_cnt++;
				}
			}
			else if(pos>=(brkpnt1-readlength) && pos<brkpnt1) // reads split at the left breakpoints 
			{
				int itemp2=vrecords->at(i)->cigar.at(cg_size-1).second; //last cigar length
				std::string strtemp2=vrecords->at(i)->cigar.at(cg_size-1).first; //last cigar 

				//the part surpass brkpnt1 mapped on the ref, and length > cnt_base_diff, then unsplit_cnt++
				//1. get part [brkpnt1,end of reads] length and cigar 
				
				int len;//save how many bases haven't been check between [brkpnt1,pos+readlength]
				int imapped=0;// save how many bases mapped between [brkpnt1,pos+readlength] 
				int index1=cg_size-1;

				if((pos+readlength)<=brkpnt2)  
				{
					len=pos+readlength-brkpnt1-itemp2;
				
					if(strtemp2=="M")
					{
						imapped+=itemp2;
					}
					else
					{
						imapped=0;
					}	
				}
				else // two sides of the reads overlap with sequence 
				{// move forward from right to left until reach brkpnt2. 
					int temp_start=pos+readlength; 
					int itemp4;
					std::string strtemp4;
					while(temp_start>brkpnt2)
					{
						itemp4=vrecords->at(i)->cigar.at(index1).second;
						strtemp4=vrecords->at(i)->cigar.at(index1).first; 

						temp_start-=itemp4;
						index1--;
					}
					if(strtemp4=="M")
					{
						imapped+=(brkpnt2-temp_start);
					}
					len=temp_start-brkpnt1;
				}
							

				while( len>=0 )
				{
					index1--;
					if(index1<0)
					{
						break;
					}
					int itemp3=vrecords->at(i)->cigar.at(index1).second;
					std::string strtemp3=vrecords->at(i)->cigar.at(index1).first; 
					
					len-=itemp3;
					if(strtemp3=="M") 
					{
						if(len>=0)
							imapped+=itemp3;
						else
						{
							imapped+=(itemp3+len);
						}
					}
					
				}
				if( imapped>= cnt_base_diff) //check whether unsplit is caused by the more than cnt_base_differ bases
				{
					unsplit_cnt++;
				}
				//the part surpass brkpnt1 not mapped 
				else//check whether split is caused by the more than cnt_base_differ bases
				{
					split_cnt++; 
				}
			}
			else // reads split at the right breakpnts 
			{
				bool b1=true;
				int len=brkpnt2-pos;
				int index1=0;
				
				int imapped;
				int temp_pos=pos; 
				if(strtemp1=="M")
				{
					if(itemp1<=len)
						imapped=itemp1;
					else
						imapped=len;
				}
				else
				{
					imapped=0;
				}
				temp_pos+=itemp1; 

				while( imapped<len )
				{
					index1++;
					if(index1>=cg_size)
					{
						break;
					}
					int itemp3=vrecords->at(i)->cigar.at(index1).second;
					std::string strtemp3=vrecords->at(i)->cigar.at(index1).first; 
					
					temp_pos+=itemp3;
					if(strtemp3=="M") 
					{
						if(temp_pos<=brkpnt2)  
							imapped+=itemp3; 
						else
						{
							imapped+=(itemp3-(temp_pos-brkpnt2)); 
							break;
						}
					}
					if(temp_pos>brkpnt2) break;
				}
				if(imapped>= cnt_base_diff) //check whether unsplit is caused by the more than cnt_base_differ bases
				{
					unsplit_cnt++;
				}
				//the part surpass brkpnt1 not mapped 
				else//check whether split is caused by the more than cnt_base_differ bases
				{
					split_cnt++; 
				}

			}// end of "else if"  

		}//end of pos
	}
}

/*
Description: Count the number of split-mapped reads. 
			And in this function, we strictly ask them map or split
*/
void ReadsParser::cntSplitMappedReads(long brkpnt1,long brkpnt2, int& split_cnt, int& unsplit_cnt) 
{
	split_cnt=0;
	unsplit_cnt=0;

	int vsize = vrecords->size();
	std::string temp_cigar;
	int pos;
	
	int interval=3; //slack value at the breakpoint

	for(int i=0;i<vsize;i++)
	{
		pos=vrecords->at(i)->pos;
		if((pos > (brkpnt1-readlength) && pos < brkpnt1) || (((pos+readlength)>brkpnt2) && (pos<brkpnt2) ))
		{			
			int cg_size = vrecords->at(i)->cigar.size();
			if(cg_size==0) 
			{
				continue;
			}

			int itemp1=vrecords->at(i)->cigar.at(0).second;
			std::string strtemp1=vrecords->at(i)->cigar.at(0).first; 

			if(strtemp1=="*") 
			{
				continue;
			}

			if(cg_size==1) //totally mapped or totally unmapped. 
			{
				if(strtemp1!="M")
				{
					//split_cnt++; 
				}
				else 
				{
					unsplit_cnt++;
				}
			}
			else if(pos>=(brkpnt1-readlength) && pos<brkpnt1) // reads split at the left breakpoints  
			{	
				int ileft=brkpnt1-pos;
				if((itemp1>ileft-interval && itemp1<ileft+interval) && strtemp1=="M")
					split_cnt++;

			}
			else // reads split at the right breakpnts 
			{
				int itemp2=vrecords->at(i)->cigar.at(cg_size-1).second; //last cigar length
				std::string strtemp2=vrecords->at(i)->cigar.at(cg_size-1).first; //last cigar 

				int iright=pos-brkpnt2;
				if((itemp2>iright-interval && itemp2<iright+interval) && strtemp2=="M")
					split_cnt++;
			}// end of "else if"  

		}//end of if pos
	}
}

/*
Description: Get number of overlapped reads over a deletion. 
// how about calc the reads coverage? 
*/
//void ReadsParser::cntOverlapReadsOverDel(long brkpnt1, long brkpnt2, int& overlap_cnt)
//{	
//	/*
//	1. first, reads should mapped; this step is not finished.....
//	2. then, reads overlaped with each other. 
//	*/
//
//	overlap_cnt=0;
//	int vsize = vrecords->size();
//	int flag=0;
//	int lastpos=brkpnt1;
//	int pos;
//	for(int i=0;i<vsize;i++)
//	{
//		pos=vrecords->at(i)->pos;
//		if(pos >= (brkpnt1-readlength) && pos <= brkpnt2)
//		{
//			if(flag==0)
//			{
//				flag=1;
//			}
//			else if(pos-lastpos<readlength)
//			{
//				overlap_cnt++;
//			}
//			else
//			{
//				flag=0;
//			}
//			
//		}
//	}
//	
//}



/*
Description:
Count the number of mapped reads and unmapped reads between [brkpnt1, brkpnt2]. 
*/

void ReadsParser::cntMapUnMapReads(long brkpnt1,long brkpnt2, int& mapped_cnt, int& unmapped_cnt)
{
	mapped_cnt=0;
	int vsize = vrecords->size();
	std::string temp_cigar;
	int pos;
	std::string tmplt = PubFuncs::cvtInt2Str(readlength) + "M"; 
	
	
	//if(brkpnt1==81094047 || brkpnt1==81072325) std::cout<<"size:"<<vsize<<std::endl; /////////////////////////////////////////////////////////////////////
	for(int i=0;i<vsize;i++)
	{
		pos=vrecords->at(i)->pos;
		
		if(pos >= brkpnt1 && pos <= (brkpnt2-readlength))
		{
			//if((brkpnt1==81094047 || brkpnt1==81072325) && i==1181) std::cout<<"pos:"<<pos<<std::endl; /////////////////////////////////////////////////////////////////////
			temp_cigar = getCigarByPos(pos);
			
			//if(brkpnt1==81094047) std::cout<<"temp cigar:"<<temp_cigar<<std::endl; /////////////////////////////////////////////////////////////////////
			
			if(temp_cigar==tmplt) // here is strict constrain 
			{
				mapped_cnt++;
			}
			else
				unmapped_cnt++;
		}
		else if(pos>brkpnt2)
		{
			//if(brkpnt1==81094047) std::cout<<"break"<<std::endl; 
			break;
		}
	}
}

/*
Description: Get the # of overlapped reads in a deletion.

And this function must excuted after function cntMapUnMapReads, because we need the positions of mapped reads
*/
void ReadsParser::cntOverlapReadsOverDel(long brkpnt1, long brkpnt2, int& overlap_cnt)
{
	overlap_cnt=0;
	int vsize = vmappedreads.size();
	std::string temp_cigar;
	int pos1,pos2;
	std::string tmplt = PubFuncs::cvtInt2Str(readlength) + "M"; 
	
	for(int i=0;i<vsize;i++)
	{
		pos1=vmappedreads[i];
		for(int j=i+1;j<vsize;j++)
		{
			pos2=vmappedreads[j];
			if((pos1+readlength)>pos2)
				overlap_cnt++;
		}
		
	}
	
} 

/*
Description: Get the cigar of a read by the position of the read, 
             and the reads are sorted according to pos in ascendent order. 
             using binary search to find the read. 
*/
std::string ReadsParser::getCigarByPos(int pos) 
{
	int vsize = vrecords->size();
	int left=0;
	int right=vsize-1;
	int mid=(left+right)/2;


	while(true)
	{
		if(left>=right) break; 
		int temppos = vrecords->at(mid)->pos;

		//if(pos==81180399) std::cout<<"left:"<<left<<" mid:"<<mid<<" right:"<<right<<" "<<vrecords->at(mid)->pos<<std::endl;//////////////////////////////////

		if(vrecords->at(mid)->pos == pos || (right-left==1))
		{
			if(right-left==1)
			{
				if(vrecords->at(left)->pos == pos)
					mid=left;
				else if(vrecords->at(right)->pos == pos)
					mid=right;
				else
					return "no cigar";
			}
			//if(pos==81180399) std::cout<<"equal:"<<vrecords->at(mid)->pos<<std::endl;//////////////////////////////////
			std::string cigar="";
			int cg_size = vrecords->at(mid)->cigar.size();
			for(int i=0;i<cg_size;i++)
			{
				cigar += PubFuncs::cvtInt2Str(vrecords->at(mid)->cigar.at(i).second);
				cigar += vrecords->at(mid)->cigar.at(i).first;
			}
			return cigar;
		}
		else if(vrecords->at(mid)->pos > pos)
		{
			//if(pos==81180399) std::cout<<"smaller"<<vrecords->at(mid)->pos<<std::endl;//////////////////////////////////
			right = mid;
		}
		else
		{
			//if(pos==81180399) std::cout<<"larger"<<vrecords->at(mid)->pos<<", left:"<< left<<" right:"<<right <<std::endl;//////////////////////////////////
			left = mid;
		}
		mid = (left+right)/2;

	}
	return "no cigar";
}


/*
Description: Compare two BamAlignmentRecords according to position 
			 It should be global function.
*/
bool cmpAR(BamAlignmentRecord* ba1, BamAlignmentRecord* ba2)
{
	return (ba1->pos < ba2->pos);
}

/*
Description: Sort alignment records according to pos field in each alignment.
*/
void ReadsParser::sortAlignRecords()
{
	std::sort(vrecords->begin(),vrecords->end(),cmpAR);
}



/*
Description:
	Check whether pair-end reads encompassing deletion and is discordant.
	And make sure both pair-end are fully mapped, this can be down by samtools.
*/
bool ReadsParser::isPEReadEncompassDel(int read_pos, int mate_pos, int lbrkpnt, int rbrkpnt)//whether a pair-end reads encompassing a deletion.
{

	//encompassing and discordant, then return true
	int criterion = this->insert_size + 3*this->std_variation;
	bool bdiscordant = (rbrkpnt - lbrkpnt)>criterion;
	bool bleft = (read_pos + this->readlength) < lbrkpnt;
	bool bright = (mate_pos - this->readlength) > rbrkpnt;
	if(bdiscordant && bleft && bright)
	{
		return true;
	}
	else
		return false;
	
}

/*
Desription:
	Check whether a read split at breakpoint.
	Make sure the reads is not fully mapped.
*/
bool ReadsParser::isReadSplitDel(int read_pos,std::vector< std::pair<std::string,int> >& cigar,int lbrkpnt,int rbrkpnt, int slack)
{
	int cg_size=cigar.size();
	bool blsplit = (read_pos > (lbrkpnt-readlength)) && (read_pos < lbrkpnt);
	bool brsplit = ((read_pos + this->readlength) > rbrkpnt) && (read_pos<rbrkpnt) ;
	if(blsplit==true && brsplit==true)
	{			
		int cg_size = cigar.size();
		if(cg_size<=1)//fully map or wrong 
		{
			return false;
		}

		int itemp1=cigar.at(0).second;
		std::string strtemp1=cigar.at(0).first; 
		if(strtemp1=="*") 
		{
			return false;
		}

		if(blsplit==true) // reads split at the left breakpoints.
		{	
			int ileft=lbrkpnt-read_pos;
			if((itemp1>(ileft-slack) && itemp1<(ileft+slack)) && strtemp1=="M")
				return true;

		}
		else if(brsplit==true) // reads split at the right breakpnts 
		{
			int itemp2=cigar.at(cg_size-1).second; //last cigar length
			std::string strtemp2=cigar.at(cg_size-1).first; //last cigar 

			int iright=read_pos - rbrkpnt + this->readlength;
			if((itemp2>(iright-slack) && itemp2<(iright+slack)) && strtemp2=="M")
				return true;
		}

	}
	
	return false;
}

/*
Description:
	Check whether a read is fully mapped between two brkpnt.
	Make sure the reads is fully mapped.
*/
bool ReadsParser::inReadOverDelFullMap(int read_pos, int lbrkpnt, int rbrkpnt)//whether a read fully mapped betweeen two breakpoints.
{
	if(read_pos>lbrkpnt && read_pos<rbrkpnt)
		return true;
	else
		return false;
}


int ReadsParser::getReadLength()
{
	return this->readlength;
}

int ReadsParser::getInsertSize()
{
	return this->insert_size;
}

int ReadsParser::getStdVariation()
{
	return this->std_variation;
}