/*
 * Mass-spectrometry Filtering Project
 * This program is property of University of California at San Diego
 *
 * Supervisor: Pavel Pevzner
 * coded by Dumitru Brinza
 *
 * ppevzner@cs.ucsd.edu
 * dima@cs.ucsd.edu
 *
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "PepnovoTags.h"

using namespace std;

int TTag::NumberOfFlankingMassesToMatch = 1;

//---------------------------------------------------------------------
string itos(int i)	// convert int to string
{
		stringstream s;
		s << i;
		return s.str();
}
//----------------------------------------------------------------------
/*
int main(int argc, char *argv[])
{
	std::vector<std::vector<TTag> > tags;
	string fname = argv[1];
	for(int i=0;i<=3;i++) LoadTags(fname+itos(i)+".txt", tags, 1);
	std::vector<unsigned int> temp;
	for(int i=0;i<10;i++)temp.push_back(i);
	IntersectTags(tags,0,temp,1,1);
	for(int i=0;i<temp.size();i++) cout << temp[i] << endl;
}
*/
//---------------------------------------------------------------------
void IntersectTags(std::vector<std::vector<TTag> > &tags, unsigned int spectrumId, std::vector<unsigned int> & candidates, int numberOfTagsToMatch, int numberOfFlankingMassesToMatch)
{
 TTag::NumberOfFlankingMassesToMatch = numberOfFlankingMassesToMatch;

 std::vector<unsigned int> subsetOfCandidates;
 int candidateId;
 for(int c=0;c<candidates.size();c++)
 {
   candidateId = candidates[c];
   int overlap = 0;
   for(int x=0;x<tags[spectrumId].size()&&overlap<numberOfTagsToMatch;x++)
   for(int y=0;y<tags[candidateId].size()&&overlap<numberOfTagsToMatch;y++)
   if(tags[spectrumId][x] == tags[candidateId][y]) overlap++;
   if(overlap>=numberOfTagsToMatch)subsetOfCandidates.push_back(candidateId);
 }

 candidates.swap(subsetOfCandidates);
}

//---------------------------------------------------------------------

void UpdateTagList(std::vector<TTag> &tag_list, TTag & tag )
{
 int start = -1;
 for(int i=0;i<tag_list.size();i++)
 {
  if(tag_list[i]==tag)
    {
     if(tag_list[i].confedence<tag.confedence) { start=i; tag_list[start].set(tag); }
     else return;
     break;
    }
 }

 if(start == -1)
	 {
	 start = tag_list.size()-1;
	 if(tag_list.size()<20) {tag_list.push_back(tag); start++;}
     else tag_list[start].set(tag);
	 }

 for(int i=start-1;i>=0;i--)
   if(tag_list[i].confedence<tag_list[i+1].confedence)
    {
     TTag temp; temp.set(tag_list[i]);
     tag_list[i].set(tag_list[i+1]);
     tag_list[i+1].set(temp);
    }
    else break;
  }
//-------------------------------------------------------------------
void LoadTags(string FileName, std::vector<std::vector<TTag> > &tags, int numberOfFlankingMassesToMatch){

	 TTag::NumberOfFlankingMassesToMatch = numberOfFlankingMassesToMatch;

	ifstream infile;
    infile.open(FileName.c_str(), ios::binary);
	string line;

   if (!infile) {
        cerr << "Unable to read " <<  FileName << endl;
        exit(1);
    }

   int num_peptides = -1, nn, i, prev_index = 0, cur_index, extra_spectra;
   string::size_type loc1,loc2;

   while(getline(infile, line))
    {
	 if(line.length()>2 && line.substr(0,2)==">>")
	 {
      std::vector<TTag> temp_list_of_tags;
     cur_index = atoi(line.substr(5).c_str())+1;
	  extra_spectra = (cur_index - prev_index + 10000) % 10000;
	  prev_index = cur_index;

          for(i=1;i<extra_spectra;i++)
    	  {
    	  num_peptides++;
    	  tags.push_back(temp_list_of_tags);
    	  }
	 num_peptides++;
	 if(tags.size()<=num_peptides) tags.push_back(temp_list_of_tags);

	 getline(infile, line);
	 while(getline(infile, line)&&line.length()>2)
   	    {
	    TTag temp_tag;
            loc1 = line.find( "\t", 0 );
	        loc2 = line.find( "\t", loc1+1 );
	        temp_tag.confedence = atof(line.substr(loc1+1,loc2-loc1-1).c_str());
	        loc1 = line.find( "\t", loc2+1 );
	        loc2 = line.find( "\t", loc1+1 );
	        temp_tag.n_gap=(unsigned short)(10*atof(line.substr(loc1+1,loc2-loc1-1).c_str()));
	        loc1 = line.find( "\t", loc2+1 );
	        temp_tag.c_gap=(unsigned short)(10*atof(line.substr(loc2+1,loc1-loc2-1).c_str()));
	        loc1 = line.find( "\t", loc1+1 );
	        loc2 = line.find( "\t", loc1+1 );
	        temp_tag.charge = (char)atoi(line.substr(loc1+1,loc2-loc1-1).c_str());
  	        memcpy(temp_tag.sequence,&line[loc2+1],3);
   	        if(temp_tag.sequence[0]=='Q')temp_tag.sequence[0]='K';
	        if(temp_tag.sequence[1]=='Q')temp_tag.sequence[1]='K';
	        if(temp_tag.sequence[2]=='Q')temp_tag.sequence[2]='K';

	        if(tags[num_peptides].size()<20 || tags[num_peptides][tags[num_peptides].size()-1].confedence<temp_tag.confedence) UpdateTagList(tags[num_peptides],temp_tag);
        }
   }
 }
   infile.close();
   cout << endl;

}
//---------------------------------------------------------------------
