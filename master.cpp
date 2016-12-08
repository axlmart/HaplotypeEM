#include <iostream> 
#include <fstream> 
#include <map> 
#include <vector> 
#include <string>
#include "Matrix615.h"
#include <cmath>
#include <set>
#include <map>
using namespace std;

//this has a problem in the enumeration of our haplotypes when doing 3 bps (because we enumeration at 4)
// idea! add indicator variable that would increase starting index to 4 if locus_size==3
// something like int blah=(locus_size==3 ? 1 : 0)*4
// then int i = 0 + blah would be 0 if locus_size==2 , 4 if locus_size==3
// we need this distinction to properly update hap_prob, print the correct observed counts table, and maybe
// something else needs it too

void init_EM( map<int, std::set <int> > pos_hapmap , map<int, int> obs_map , Matrix615<int>&
	my_matrix , map <int , vector<int> >& pos_hap_all_markers , int locus_num){
	
	//pos_hapmap is the map of possible haplotypes
	//obs_map counts the number of observed haplotypes
	//my_matrix is the data file
	//pos_hap_all_markers gets updated with possible haplotypes at each 2 or 3 bp "locus"
	//locus_num is the set of 2 or 3 bps we are looking at. Locus num ranges from o to floor(#markers)
	
	int locus_size = 2; //start out assuming we are looking at 2 markers
	if ((my_matrix.colNums()-locus_num*2)==3){ locus_size=3;} //if we are looking at the last locus
								//and it is of size 3, account for it
	int index_adj=(locus_size==3 ? 1 : 0)*4;
	
	//should allocate and delete this memory
	map<int, double> hap_prob;
	for (int i=0+index_adj ; i<pow(2,locus_size)+index_adj; ++i){
		hap_prob[i]= 1 / double(pow(2,locus_size));
	}
	
	for (int n=0; n<my_matrix.rowNums(); ++n){//for N rows in locus matrix 
		if (locus_size==2){// if num markers in the locus =2
			for(int m=0; m<locus_size; ++m){ //for each marker in the locus of size 2
				if(m==0){
					if(my_matrix.data[n][locus_num*2 + m]==0){//if we have two 0s at the first marker
						pos_hapmap[n].erase(2); pos_hapmap[n].erase(3); }//remove 2,3 from pos_hapmap[n]
					else if(my_matrix.data[n][locus_num*2 + m]==11){//if we have two ones at the first marker
						pos_hapmap[n].erase(0); pos_hapmap[n].erase(1); } //remove 0,1
				} else{
					if(my_matrix.data[n][locus_num*2 + m]==0){//if we have two  0s at the second marker
						pos_hapmap[n].erase(1); pos_hapmap[n].erase(3); } //remove 1 and 3
					else if(my_matrix.data[n][locus_num*2 + m]==11){ //if we have two 1s at the second marker
						pos_hapmap[n].erase(0); pos_hapmap[n].erase(2); } //remove 0 and 2
				}
			}
		}
		
		else if (locus_size==3){// if num markers in the locus = 3
			for (int m=0; m<locus_size; ++m){ //for each marker in the locus of size 3
				if (m==0){
					if (my_matrix.data[n][locus_num*2 + m]==0){
						pos_hapmap[n].erase(8); pos_hapmap[n].erase(9);
						pos_hapmap[n].erase(10); pos_hapmap[n].erase(11);}
					else if(my_matrix.data[n][locus_num*2 + m]==11){
						pos_hapmap[n].erase(4); pos_hapmap[n].erase(5);
						pos_hapmap[n].erase(6); pos_hapmap[n].erase(7);}
				} else if (m==1){
					if (my_matrix.data[n][locus_num*2 + m]==0){
						pos_hapmap[n].erase(6); pos_hapmap[n].erase(7);
						pos_hapmap[n].erase(10); pos_hapmap[n].erase(11);}
					else if(my_matrix.data[n][locus_num*2 + m]==11){
						pos_hapmap[n].erase(4); pos_hapmap[n].erase(5);
						pos_hapmap[n].erase(8); pos_hapmap[n].erase(9);}
				} else{
					if (my_matrix.data[n][locus_num*2 + m]==0){
						pos_hapmap[n].erase(5); pos_hapmap[n].erase(7);
						pos_hapmap[n].erase(9); pos_hapmap[n].erase(11);}
					else if(my_matrix.data[n][locus_num*2 + m]==11){
						pos_hapmap[n].erase(4); pos_hapmap[n].erase(6);
						pos_hapmap[n].erase(8); pos_hapmap[n].erase(10);}
				}
			}
		}
		else{}
		
		std::set<int>::iterator it;
		
		if (pos_hapmap[n].size()==2){//if size pos_hapmap[n]==2 iterate through set using *iter and  ++obs_map[*it]
			for(it=pos_hapmap[n].begin(); it!=pos_hapmap[n].end(); ++it){
				++obs_map[*it];
			}
		}
		else if (pos_hapmap[n].size()==1){ //if size pos_hapmap[n]==1 obs_map[*it]=obs_map[*it]+2
			it=pos_hapmap[n].begin();
			obs_map[*it]=obs_map[*it]+2;
		}
		else{}
		
	}
			
	//add something to print the observed haps
	
	std::set<int>::iterator it2;
	
	cout << "Possible Haplotypes" << endl << endl;
	
	for (int i=0; i<my_matrix.rowNums(); i++){
		cout << "Person " << i+1 << ": ";
		for (it2=pos_hapmap[i].begin(); it2!=pos_hapmap[i].end(); ++it2){
			cout<< *it2 << " ";
		}
		cout<<endl;
	}
	
	cout << endl << endl << "Observed Counts" << endl << endl;
	
	for (int i=0+index_adj; i<pow(2,locus_size)+index_adj; ++i){
		cout << "Haplotype " << i << ": ";
		cout<< obs_map[i];
		cout << endl;
	}	
		
	//Run EM udating observed map and hap_prob (still needs to be added)
	//push back hap probs >0 into pos_hap_all_markers[locus_num]
	
	for (int i=0+index_adj; i<pow(2,locus_size)+index_adj; ++i){
		if (hap_prob[i] > 0.000001){ //change to fit whatever we decide
			pos_hap_all_markers[locus_num].push_back(i);
		}
	}
	
	
	
}
		
	
			
	
	
	
	
	


int main(int argc, char** argv){
	Matrix615<int> my_matrix;
	my_matrix.readFromFile(argv[1]);


// Divide the matrix in submatrices of 2 columns including all rows.
// Call outside EM function on each subproblem. Hi Axel!

// Note for now (at least) the entries are intergers (because changing the Matrix.h 
// to allow for strings seemed complicated). Thus each value represents a different genotype :
// integer value --> genotype 

// 0 --> 00
// 1 --> 01
// 11 --> 11
//
//
// Note we will need to modify this loop to account for the possibility of odd # of markers.

	int N = my_matrix.rowNums();
	int M = my_matrix.colNums();

/*create a map of possible haplotypes where each key is a person. 
0 corresponds to <0,0>
1 corresponds to <0,1>
2 = <1,0>
3 = <1,1>
4 = <0,0,0>
5 = <0,0,1>
6 = <0,1,0>
7 = <0,1,1>
8 = <1,0,0>
9 = <1,0,1>
10 = <1,1,0>
11 = <1,1,1>
*/
 
	int hap_nums_2bp[]={0,1,2,3}; //arrays to be turned into sets
	int hap_nums_3bp[]={4,5,6,7,8,9,10,11};
	std::set<int> hap_nums_2bp_set (hap_nums_2bp, hap_nums_2bp+4);
	std::set<int> hap_nums_3bp_set (hap_nums_3bp, hap_nums_3bp+8);

	map<int, std::set <int> > master_hapmap_2bp;

	for(int n=0; n<N; ++n){
		master_hapmap_2bp[n] = hap_nums_2bp_set;
	}

	map<int , std::set <int> > master_hapmap_3bp;

	for( int n=0; n<N; ++n){
		master_hapmap_3bp[n] = hap_nums_3bp_set;
	}
	
	map<int, int> observed_map_2bp; //these keep count of the number of haplotypes we observed
	map<int, int> observed_map_3bp;
	
	map<int, vector<int> > pos_hap_all_markers;
	
	
	//add for loop to go through every 2/3 bp hap
	
	for (int i=0; i< floor((double)M/2); i++){
		cout << "Marker Pair " << i << endl;
		if ((M-i*2)==3){
			init_EM(master_hapmap_3bp, observed_map_3bp, my_matrix,
				pos_hap_all_markers, i);
		}
		else{
			init_EM(master_hapmap_2bp, observed_map_2bp, my_matrix,
				pos_hap_all_markers, i);
		}
		cout << endl;
	}
	
	
	//print pos_haps (should be all at this point since we haven't added the EM part)
	
	for (int i=0; i< floor((double)M/2); ++i){
		cout << "Marker Pair " << i << endl;
		for (int j=0; j<pos_hap_all_markers[i].size(); j++){
			cout << pos_hap_all_markers[i][j]<< " ";
		}
		cout<<endl;
	}

	return 0;	

}

	
	
	


 





